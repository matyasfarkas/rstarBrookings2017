##########################################################################################
## CASE 2C: MATCH EA CONVENIENCE-YIELD DYNAMICS IN THE US MODEL
##########################################################################################
# Goal:
#   Build a US conditional forecast where the US convenience yield follows the EA
#   convenience-yield dynamics after 1998Q4. We keep the US initial state and the
#   US non-CY shocks fixed. The only shocks solved period by period are the four
#   CY shocks:
#       b_liqtil_sh, b_liqp_sh, b_safetil_sh, b_safep_sh
#
# Convenience yield is always defined as the sum of the four CY states:
#       b_liqtil_t + b_liqp_t + b_safetil_t + b_safep_t
##########################################################################################

using DSGE, HDF5, Plots, StatsPlots
using DataFrames, Dates
using LinearAlgebra

default(
    titlefontsize = 18,
    guidefontsize = 13,
    tickfontsize = 13,
    legendfontsize = 12,
    foreground_color_text = :black,
    foreground_color_axis = :black,
    foreground_color_guide = :black,
    foreground_color_border = :black,
)

##########################################################################################
## SETTINGS
##########################################################################################

use_FG_in_EA = true
smoother = :durbin_koopman
anchor_date = quartertodate("1998-Q4")
first_panel_legend_pos = (0.25, 0.25)

script_dir = @__DIR__
basepath = dirname(script_dir)
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")
figroot = joinpath(saveroot, "Final Paper", "Figures")

cy_state_names = [:b_liqtil_t, :b_liqp_t, :b_safetil_t, :b_safep_t]
cy_shock_names = [:b_liqtil_sh, :b_liqp_sh, :b_safetil_sh, :b_safep_sh]

ea_spec = "ss24"
ea_vintage = use_FG_in_EA ? "250115" : "250116"
ea_mode_file = joinpath(
    saveroot,
    "output_data",
    "m1010",
    ea_spec,
    "estimate",
    "raw",
    "paramsmode_vint=$(ea_vintage).h5",
)

us_spec = "ss20"
us_vintage = "250825"
us_mode_file = joinpath(
    saveroot,
    "output_data",
    "m1010",
    us_spec,
    "estimate",
    "raw",
    "paramsmode_vint=$(us_vintage).h5",
)

##########################################################################################
## SMALL HELPERS
##########################################################################################

function set_model_settings!(m, dataroot, saveroot, vintage, smoother)
    m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
    m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
    m <= DSGE.Setting(:data_vintage, vintage)
    m <= DSGE.Setting(:reoptimize, false)
    m <= DSGE.Setting(:calculate_hessian, false)
    m <= DSGE.Setting(:date_mainsample_start, quartertodate("1970-Q3"))
    m <= DSGE.Setting(:date_presample_start, quartertodate("1970-Q2"))
    m <= DSGE.Setting(:date_forecast_start, quartertodate("2024-Q3"))
    m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))
    m <= DSGE.Setting(:forecast_smoother, smoother)
end

function load_and_smooth_model(spec, vintage, mode_file, dataroot, saveroot, smoother)
    m = Model1010(spec)
    set_model_settings!(m, dataroot, saveroot, vintage, smoother)
    df = load_data(m; check_empty_columns = false)
    specify_mode!(m, mode_file)
    system = DSGE.compute_system(m)
    states_sm, shocks_sm, pseudo_sm = DSGE.smooth(m, df, system; draw_states = false)
    T = size(states_sm, 2)
    df_aligned = df[end-T+1:end, :]
    dates = df_aligned.date
    return m, df_aligned, system, states_sm, shocks_sm, pseudo_sm, dates
end

function model_series(m, states_mat, obs_mat, pseudo_mat, sym::Symbol)
    if haskey(m.endogenous_states, sym)
        return vec(states_mat[m.endogenous_states[sym], :])
    elseif haskey(m.endogenous_states_augmented, sym)
        return vec(states_mat[m.endogenous_states_augmented[sym], :])
    elseif haskey(m.observables, sym)
        obs_mat === nothing && error("Observable $(sym) requested, but no observables matrix was provided.")
        return vec(obs_mat[m.observables[sym], :])
    elseif haskey(m.pseudo_observables, sym)
        return vec(pseudo_mat[m.pseudo_observables[sym], :])
    else
        error("Series $(sym) not found in model.")
    end
end

function baseline_series(m, df_aligned::DataFrame, states_sm, pseudo_sm, sym::Symbol)
    if haskey(m.observables, sym)
        sym in names(df_aligned) || error("Observable $(sym) not found in data.")
        return vec(df_aligned[!, sym])
    else
        return model_series(m, states_sm, nothing, pseudo_sm, sym)
    end
end

function convenience_yield(m, states_mat, pseudo_mat, cy_state_names)
    # Total convenience yield is the sum of the four smoothed CY state variables.
    missing_states = [
        sym for sym in cy_state_names
        if !(haskey(m.endogenous_states, sym) || haskey(m.endogenous_states_augmented, sym))
    ]
    isempty(missing_states) || error("Convenience-yield state(s) not found: $(missing_states)")

    cy = model_series(m, states_mat, nothing, pseudo_mat, cy_state_names[1])
    for sym in cy_state_names[2:end]
        cy = cy .+ model_series(m, states_mat, nothing, pseudo_mat, sym)
    end
    return cy
end

function resolve_long_rate(m)
    candidates = [
        :obs_longrate, :obs_long_rate, :longrate, :LongRate,
        :LongTermRate, :LongTermNominalRate, :tenyearrate, :obs_tenyearrate
    ]
    for sym in candidates
        if haskey(m.endogenous_states, sym) ||
           haskey(m.endogenous_states_augmented, sym) ||
           haskey(m.observables, sym) ||
           haskey(m.pseudo_observables, sym)
            return sym
        end
    end
    error("Could not locate the long-term interest rate series.")
end

function replace_shock_block(base_shocks, m, shock_names, shock_vals)
    size(base_shocks, 2) == size(shock_vals, 1) || error("Shock replacement horizon mismatch.")
    size(shock_vals, 2) == length(shock_names) || error("Shock replacement width mismatch.")

    shocks_new = copy(base_shocks)
    for (i, shock_name) in enumerate(shock_names)
        shock_ind = m.exogenous_shocks[shock_name]
        shocks_new[shock_ind, :] .= shock_vals[:, i]
    end
    return shocks_new
end

function solve_cy_target_shocks(system_US, m_US, anchor_state, base_shocks_tail, cy_target_tail, shock_names; shock_prior_vals)
    # Starting from the US state and US non-CY shocks, adjust only the CY shocks
    # so that the next-period sum of CY states exactly equals the target path.
    horizon_tail = size(base_shocks_tail, 2)
    length(cy_target_tail) == horizon_tail || error("CY target horizon mismatch.")

    TTT = system_US[:TTT]
    RRR = system_US[:RRR]
    CCC = system_US[:CCC]

    cy_inds = [
        haskey(m_US.endogenous_states, sym) ? m_US.endogenous_states[sym] : m_US.endogenous_states_augmented[sym]
        for sym in cy_state_names
    ]
    shock_inds = [m_US.exogenous_shocks[shock_name] for shock_name in shock_names]
    cy_loading = vec(sum(RRR[cy_inds, shock_inds], dims = 1))
    loading_norm2 = dot(cy_loading, cy_loading)
    loading_norm2 > eps(eltype(cy_loading)) || error("Selected shocks cannot move convenience yield.")

    shocks_conditional = copy(base_shocks_tail)
    s_prev = collect(anchor_state)

    for t in 1:horizon_tail
        shocks_conditional[shock_inds, t] .= shock_prior_vals[t, :]

        # Forecast once using the prior shocks, then close the gap to the CY target.
        s_before_targeting = CCC + TTT * s_prev + RRR * shocks_conditional[:, t]
        cy_gap = cy_target_tail[t] - sum(s_before_targeting[cy_inds])

        shocks_conditional[shock_inds, t] .+= (cy_gap / loading_norm2) .* cy_loading
        s_prev = CCC + TTT * s_prev + RRR * shocks_conditional[:, t]
    end

    return shocks_conditional
end

function splice_forecast_tail(baseline, actual_tail, counterfactual_tail; anchor_idx)
    length(actual_tail) == length(counterfactual_tail) || error("Tail length mismatch.")
    length(baseline) == anchor_idx + length(actual_tail) || error("Baseline/tail length mismatch.")

    out = collect(baseline)
    out[anchor_idx+1:end] .= baseline[anchor_idx+1:end] .+ (counterfactual_tail .- actual_tail)
    return out
end

function four_quarter_sum(x)
    y = fill(NaN, length(x))
    for t in 4:length(x)
        y[t] = x[t] + x[t-1] + x[t-2] + x[t-3]
    end
    return y
end

function make_panel(dates_plot, baseline, red, blue, title, xticks; show_legend = false)
    p = plot(
        dates_plot,
        baseline,
        color = :black,
        lw = 2,
        linestyle = :solid,
        title = title,
        label = show_legend ? "US baseline" : "",
        xticks = xticks,
    )
    plot!(
        p,
        dates_plot,
        red,
        color = :red,
        lw = 2,
        linestyle = :dashdot,
        label = show_legend ? "US excluding CY shocks" : "",
    )
    plot!(
        p,
        dates_plot,
        blue,
        color = :blue,
        lw = 2,
        linestyle = :dash,
        label = show_legend ? "Counterfactual: EA CY dynamics" : "",
    )
    show_legend && plot!(p, legend = first_panel_legend_pos)
    return p
end

##########################################################################################
## LOAD AND SMOOTH EA AND US MODELS
##########################################################################################

m_EA, df_EA, system_EA, states_EA, shocks_EA, pseudo_EA, dates_EA =
    load_and_smooth_model(ea_spec, ea_vintage, ea_mode_file, dataroot, saveroot, smoother)

m_US, df_US, system_US, states_US, shocks_US, pseudo_US, dates_US =
    load_and_smooth_model(us_spec, us_vintage, us_mode_file, dataroot, saveroot, smoother)

common_dates = intersect(dates_EA, dates_US)
isempty(common_dates) && error("No overlapping dates between EA and US samples.")

idx_EA = findall(in(common_dates), dates_EA)
idx_US = findall(in(common_dates), dates_US)
dates_common_EA = dates_EA[idx_EA]
dates_common_US = dates_US[idx_US]
dates_common_EA == dates_common_US || error("EA/US date alignment failed.")
dates_common = dates_common_EA

anchor_idx = findfirst(==(anchor_date), dates_common)
anchor_idx === nothing && error("Anchor date $(anchor_date) not found in common sample.")
anchor_idx < 4 && error("Need at least 3 quarters before the anchor to construct y/y inflation.")

mask_plot = dates_common .>= anchor_date

df_EA_common = df_EA[idx_EA, :]
states_EA_common = states_EA[:, idx_EA]
pseudo_EA_common = pseudo_EA[:, idx_EA]

df_US_common = df_US[idx_US, :]
states_US_common = states_US[:, idx_US]
pseudo_US_common = pseudo_US[:, idx_US]

##########################################################################################
## BUILD BASELINE SERIES AND EA CY TARGET
##########################################################################################

longrate_sym = resolve_long_rate(m_US)

cy_US_raw = convenience_yield(m_US, states_US_common, pseudo_US_common, cy_state_names)
cy_EA_raw = convenience_yield(m_EA, states_EA_common, pseudo_EA_common, cy_state_names)

long_US_raw = baseline_series(m_US, df_US_common, states_US_common, pseudo_US_common, longrate_sym)
policy_US_raw = baseline_series(m_US, df_US_common, states_US_common, pseudo_US_common, :obs_nominalrate)
infl_US_raw = baseline_series(m_US, df_US_common, states_US_common, pseudo_US_common, :obs_corepce)
output_US_raw = baseline_series(m_US, df_US_common, states_US_common, pseudo_US_common, :y_t)
rstar_US_raw = baseline_series(m_US, df_US_common, states_US_common, pseudo_US_common, :Forward5YearRealNaturalRate)

cy_target_raw = copy(cy_US_raw)
# Use the EA change in CY after 1998Q4, but start from the US 1998Q4 level.
cy_target_raw[anchor_idx+1:end] .=
    cy_US_raw[anchor_idx] .+ (cy_EA_raw[anchor_idx+1:end] .- cy_EA_raw[anchor_idx])

##########################################################################################
## CONDITIONAL FORECASTS
##########################################################################################

tail_dates_idx_US = idx_US[anchor_idx+1:end]
tail_dates_idx_EA = idx_EA[anchor_idx+1:end]
tail_horizon = length(tail_dates_idx_US)

anchor_state_US = states_US_common[:, anchor_idx]
us_shocks_tail = copy(shocks_US[:, tail_dates_idx_US])

zero_cy_shocks = zeros(tail_horizon, length(cy_shock_names))
# Red path: remove US CY shocks, leaving all other US shocks unchanged.
red_shocks_tail = replace_shock_block(us_shocks_tail, m_US, cy_shock_names, zero_cy_shocks)

ea_cy_shock_prior = Float64.(hcat(
    [vec(shocks_EA[m_EA.exogenous_shocks[shock_name], tail_dates_idx_EA]) for shock_name in cy_shock_names]...
))
# Blue path: begin from EA CY shocks, then solve the exact CY-shock sequence
# needed for the US model to reproduce the EA CY target dynamics.
blue_shocks_tail = solve_cy_target_shocks(
    system_US,
    m_US,
    anchor_state_US,
    us_shocks_tail,
    cy_target_raw[anchor_idx+1:end],
    cy_shock_names;
    shock_prior_vals = ea_cy_shock_prior,
)

states_actual, obs_actual, pseudo_actual = forecast(system_US, collect(anchor_state_US), us_shocks_tail)
states_red, obs_red, pseudo_red = forecast(system_US, collect(anchor_state_US), red_shocks_tail)
states_blue, obs_blue, pseudo_blue = forecast(system_US, collect(anchor_state_US), blue_shocks_tail)

cy_actual_tail = convenience_yield(m_US, states_actual, pseudo_actual, cy_state_names)
cy_red_tail = convenience_yield(m_US, states_red, pseudo_red, cy_state_names)
cy_blue_tail = convenience_yield(m_US, states_blue, pseudo_blue, cy_state_names)

cy_match_error = maximum(abs.(cy_blue_tail .- cy_target_raw[anchor_idx+1:end]))
cy_match_error <= 1e-8 || @warn "Case2c CY target match error" cy_match_error

cy_red_raw = splice_forecast_tail(cy_US_raw, cy_actual_tail, cy_red_tail; anchor_idx = anchor_idx)
cy_blue_raw = cy_target_raw

long_red_raw = splice_forecast_tail(
    long_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, longrate_sym),
    model_series(m_US, states_red, obs_red, pseudo_red, longrate_sym);
    anchor_idx = anchor_idx,
)
long_blue_raw = splice_forecast_tail(
    long_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, longrate_sym),
    model_series(m_US, states_blue, obs_blue, pseudo_blue, longrate_sym);
    anchor_idx = anchor_idx,
)

policy_red_raw = splice_forecast_tail(
    policy_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :obs_nominalrate),
    model_series(m_US, states_red, obs_red, pseudo_red, :obs_nominalrate);
    anchor_idx = anchor_idx,
)
policy_blue_raw = splice_forecast_tail(
    policy_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :obs_nominalrate),
    model_series(m_US, states_blue, obs_blue, pseudo_blue, :obs_nominalrate);
    anchor_idx = anchor_idx,
)

infl_red_raw = splice_forecast_tail(
    infl_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :obs_corepce),
    model_series(m_US, states_red, obs_red, pseudo_red, :obs_corepce);
    anchor_idx = anchor_idx,
)
infl_blue_raw = splice_forecast_tail(
    infl_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :obs_corepce),
    model_series(m_US, states_blue, obs_blue, pseudo_blue, :obs_corepce);
    anchor_idx = anchor_idx,
)

output_red_raw = splice_forecast_tail(
    output_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :y_t),
    model_series(m_US, states_red, obs_red, pseudo_red, :y_t);
    anchor_idx = anchor_idx,
)
output_blue_raw = splice_forecast_tail(
    output_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :y_t),
    model_series(m_US, states_blue, obs_blue, pseudo_blue, :y_t);
    anchor_idx = anchor_idx,
)

rstar_red_raw = splice_forecast_tail(
    rstar_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :Forward5YearRealNaturalRate),
    model_series(m_US, states_red, obs_red, pseudo_red, :Forward5YearRealNaturalRate);
    anchor_idx = anchor_idx,
)
rstar_blue_raw = splice_forecast_tail(
    rstar_US_raw,
    model_series(m_US, states_actual, obs_actual, pseudo_actual, :Forward5YearRealNaturalRate),
    model_series(m_US, states_blue, obs_blue, pseudo_blue, :Forward5YearRealNaturalRate);
    anchor_idx = anchor_idx,
)

##########################################################################################
## CONVERT TO PLOTTED UNITS
##########################################################################################

cy_US = 4 .* cy_US_raw
cy_red = 4 .* cy_red_raw
cy_blue = 4 .* cy_blue_raw
cy_target = 4 .* cy_target_raw

long_US = 4 .* long_US_raw
long_red = 4 .* long_red_raw
long_blue = 4 .* long_blue_raw

policy_US = 4 .* policy_US_raw
policy_red = 4 .* policy_red_raw
policy_blue = 4 .* policy_blue_raw

infl_US = four_quarter_sum(infl_US_raw)
infl_red = four_quarter_sum(infl_red_raw)
infl_blue = four_quarter_sum(infl_blue_raw)

output_US = output_US_raw
output_red = output_red_raw
output_blue = output_blue_raw

rstar_US = 4 .* rstar_US_raw
rstar_red = 4 .* rstar_red_raw
rstar_blue = 4 .* rstar_blue_raw

##########################################################################################
## PLOT
##########################################################################################

first_tick_year = Dates.year(anchor_date)
last_tick_year = Dates.year(dates_common[end])
tick_step = 5
year_tick_years = collect(first_tick_year:tick_step:last_tick_year)
year_tick_dates = [quartertodate("$(y)-Q1") for y in year_tick_years]
year_tick_labels = string.(year_tick_years)
xticks_common = (year_tick_dates, year_tick_labels)

plot_dates = dates_common[mask_plot]

p_cy = make_panel(
    plot_dates,
    cy_US[mask_plot],
    cy_red[mask_plot],
    cy_blue[mask_plot],
    "Convenience yield (APR)",
    xticks_common;
    show_legend = true,
)
p_long = make_panel(plot_dates, long_US[mask_plot], long_red[mask_plot], long_blue[mask_plot], "Long-term interest rate (APR)", xticks_common)
p_pol = make_panel(plot_dates, policy_US[mask_plot], policy_red[mask_plot], policy_blue[mask_plot], "Policy rate (APR)", xticks_common)
p_inf = make_panel(plot_dates, infl_US[mask_plot], infl_red[mask_plot], infl_blue[mask_plot], "Core PCE inflation (%, yoy)", xticks_common)
p_out = make_panel(plot_dates, output_US[mask_plot], output_red[mask_plot], output_blue[mask_plot], "Output deviation from trend (%)", xticks_common)
p_rstar = make_panel(plot_dates, rstar_US[mask_plot], rstar_red[mask_plot], rstar_blue[mask_plot], "r* (Forward 5-year real natural rate, APR)", xticks_common)

p_3x2 = plot(
    p_cy, p_long,
    p_pol, p_inf,
    p_out, p_rstar,
    layout = (3, 2),
    size = (1100, 1200),
)

display(p_3x2)

if use_FG_in_EA
    savefig(joinpath(figroot, "Case2c Exorbitant privilege 3x2 baseline_plus_delta.pdf"))
    savefig(joinpath(figroot, "Case2c Exorbitant privilege 3x2 baseline_plus_delta.png"))
else
    savefig(joinpath(figroot, "Case2c Exorbitant privilege without FG shocks in EA 3x2 baseline_plus_delta.pdf"))
    savefig(joinpath(figroot, "Case2c Exorbitant privilege without FG shocks in EA 3x2 baseline_plus_delta.png"))
end
