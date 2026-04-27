##########################################################################################
## CASE 2D: EA NO EXORBITANT PRIVILEGE
##########################################################################################
#
# This is the EA version of the Case2B "No Exorbitant Privilege" chart.
# Black is the smoothed EA baseline. Red is the EA conditional forecast obtained
# by setting the four convenience-yield shocks to zero after 1998Q4:
#     b_liqtil_sh, b_liqp_sh, b_safetil_sh, b_safep_sh
# The output panel also removes the contribution of the EA mu shock, stored in
# this model as the Greek-symbol shock name represented below by "\u03bc_sh".
#
# The plotted convenience yield matches the pseudo-observable definition:
#     b_liq_t + b_safe_t + 100*log(lnb_liq) + 100*log(lnb_safe)
# which is equivalent to the four convenience-yield component states plus the
# model-specific drift term.
##########################################################################################

using DSGE, HDF5, Plots, StatsPlots
using DataFrames, Dates
using Statistics

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

cy_component_state_names = [:b_liqtil_t, :b_liqp_t, :b_safetil_t, :b_safep_t]
cy_shock_names = [:b_liqtil_sh, :b_liqp_sh, :b_safetil_sh, :b_safep_sh]
mu_shock_candidates = [:mu_sh, Symbol("\u03bc_sh")]

##########################################################################################
## HELPERS
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

function load_and_smooth_ea_model(spec, vintage, mode_file, dataroot, saveroot, smoother)
    m = Model1010(spec)
    set_model_settings!(m, dataroot, saveroot, vintage, smoother)
    df = load_data(m; check_empty_columns = false)
    specify_mode!(m, mode_file)
    system = DSGE.compute_system(m)
    states_sm, shocks_sm, pseudo_sm = DSGE.smooth(m, df, system; draw_states = false)
    T = size(states_sm, 2)
    df_aligned = df[end-T+1:end, :]
    return m, df_aligned, system, states_sm, shocks_sm, pseudo_sm, df_aligned.date
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

function first_available_series(m, df_aligned::DataFrame, candidates::Vector{Symbol})
    for sym in candidates
        if haskey(m.observables, sym) && sym in names(df_aligned)
            return sym
        elseif haskey(m.endogenous_states, sym) ||
               haskey(m.endogenous_states_augmented, sym) ||
               haskey(m.pseudo_observables, sym)
            return sym
        end
    end
    error("Could not locate any of these series: $(candidates)")
end

function first_available_shock(m, candidates::Vector{Symbol})
    for sym in candidates
        haskey(m.exogenous_shocks, sym) && return sym
    end
    error("Could not locate any of these shocks: $(candidates)")
end

has_model_series(m, sym::Symbol) =
    haskey(m.endogenous_states, sym) || haskey(m.endogenous_states_augmented, sym)

function convenience_yield_drift(m)
    return 100 * log(m[:lnb_liq]) + 100 * log(m[:lnb_safe])
end

function convenience_yield(m, states_mat, pseudo_mat)
    if has_model_series(m, :b_liq_t) && has_model_series(m, :b_safe_t)
        cy = model_series(m, states_mat, nothing, pseudo_mat, :b_liq_t)
        cy .+= model_series(m, states_mat, nothing, pseudo_mat, :b_safe_t)
    else
        missing_states = [
            sym for sym in cy_component_state_names if !has_model_series(m, sym)
        ]
        isempty(missing_states) || error("Convenience-yield state(s) not found: $(missing_states)")

        cy = model_series(m, states_mat, nothing, pseudo_mat, cy_component_state_names[1])
        for sym in cy_component_state_names[2:end]
            cy .+= model_series(m, states_mat, nothing, pseudo_mat, sym)
        end
    end

    return cy .+ convenience_yield_drift(m)
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

function apply_counterfactual_tail(baseline, actual_tail, counterfactual_tail; anchor_idx)
    length(actual_tail) == length(counterfactual_tail) || error("Tail series length mismatch.")
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

function make_panel(dates_plot, baseline, zero_cy, title, xticks; show_legend = false)
    p = plot(
        dates_plot,
        baseline,
        color = :black,
        lw = 2,
        linestyle = :solid,
        title = title,
        label = show_legend ? "EA - Actual" : "",
        xticks = xticks,
    )
    plot!(
        p,
        dates_plot,
        zero_cy,
        color = :red,
        lw = 2,
        linestyle = :dashdot,
        label = show_legend ? "EA - No CY innovations" : "",
    )
    show_legend && plot!(p, legend = first_panel_legend_pos)
    return p
end

##########################################################################################
## LOAD EA MODEL AND BUILD THE ZERO-CY-SHOCK FORECAST
##########################################################################################

m_EA, df_EA, system_EA, states_EA, shocks_EA, pseudo_EA, dates_EA =
    load_and_smooth_ea_model(ea_spec, ea_vintage, ea_mode_file, dataroot, saveroot, smoother)

anchor_idx = findfirst(==(anchor_date), dates_EA)
anchor_idx === nothing && error("Anchor date $(anchor_date) not found in the EA sample.")
anchor_idx < 4 && error("Need at least 3 quarters before the anchor to construct y/y inflation.")

mask_plot = dates_EA .>= anchor_date
tail_idx = anchor_idx+1:length(dates_EA)
tail_horizon = length(tail_idx)

longrate_sym = first_available_series(
    m_EA,
    df_EA,
    [:obs_longrate, :obs_long_rate, :longrate, :LongRate, :LongTermRate, :LongTermNominalRate, :tenyearrate, :obs_tenyearrate],
)
inflation_sym = first_available_series(m_EA, df_EA, [:obs_corepce, :obs_gdpdeflator, :pi_t])
mu_shock_name = first_available_shock(m_EA, mu_shock_candidates)
mu_shock_names = [mu_shock_name]

cy_base_raw = convenience_yield(m_EA, states_EA, pseudo_EA)
long_base_raw = baseline_series(m_EA, df_EA, states_EA, pseudo_EA, longrate_sym)
policy_base_raw = baseline_series(m_EA, df_EA, states_EA, pseudo_EA, :obs_nominalrate)
infl_base_raw = baseline_series(m_EA, df_EA, states_EA, pseudo_EA, inflation_sym)
output_base_raw = baseline_series(m_EA, df_EA, states_EA, pseudo_EA, :y_t)
rstar_base_raw = baseline_series(m_EA, df_EA, states_EA, pseudo_EA, :Forward5YearRealNaturalRate)

anchor_state = states_EA[:, anchor_idx]
actual_shocks_tail = copy(shocks_EA[:, tail_idx])
zero_cy_shocks_tail = replace_shock_block(
    actual_shocks_tail,
    m_EA,
    cy_shock_names,
    zeros(tail_horizon, length(cy_shock_names)),
)
zero_mu_shocks_tail = replace_shock_block(
    actual_shocks_tail,
    m_EA,
    mu_shock_names,
    zeros(tail_horizon, length(mu_shock_names)),
)
zero_cy_zero_mu_shocks_tail = replace_shock_block(
    zero_cy_shocks_tail,
    m_EA,
    mu_shock_names,
    zeros(tail_horizon, length(mu_shock_names)),
)

states_actual, obs_actual, pseudo_actual =
    forecast(system_EA, collect(anchor_state), actual_shocks_tail)
states_zero_cy, obs_zero_cy, pseudo_zero_cy =
    forecast(system_EA, collect(anchor_state), zero_cy_shocks_tail)
states_zero_mu, obs_zero_mu, pseudo_zero_mu =
    forecast(system_EA, collect(anchor_state), zero_mu_shocks_tail)
states_zero_cy_zero_mu, obs_zero_cy_zero_mu, pseudo_zero_cy_zero_mu =
    forecast(system_EA, collect(anchor_state), zero_cy_zero_mu_shocks_tail)

cy_actual_tail = convenience_yield(m_EA, states_actual, pseudo_actual)
cy_zero_cy_tail = convenience_yield(m_EA, states_zero_cy, pseudo_zero_cy)
cy_zero_cy_raw = apply_counterfactual_tail(cy_base_raw, cy_actual_tail, cy_zero_cy_tail; anchor_idx = anchor_idx)

long_zero_cy_raw = apply_counterfactual_tail(
    long_base_raw,
    model_series(m_EA, states_actual, obs_actual, pseudo_actual, longrate_sym),
    model_series(m_EA, states_zero_cy, obs_zero_cy, pseudo_zero_cy, longrate_sym);
    anchor_idx = anchor_idx,
)
policy_zero_cy_raw = apply_counterfactual_tail(
    policy_base_raw,
    model_series(m_EA, states_actual, obs_actual, pseudo_actual, :obs_nominalrate),
    model_series(m_EA, states_zero_cy, obs_zero_cy, pseudo_zero_cy, :obs_nominalrate);
    anchor_idx = anchor_idx,
)
infl_zero_cy_raw = apply_counterfactual_tail(
    infl_base_raw,
    model_series(m_EA, states_actual, obs_actual, pseudo_actual, inflation_sym),
    model_series(m_EA, states_zero_cy, obs_zero_cy, pseudo_zero_cy, inflation_sym);
    anchor_idx = anchor_idx,
)
output_base_ex_mu_raw = apply_counterfactual_tail(
    output_base_raw,
    model_series(m_EA, states_actual, obs_actual, pseudo_actual, :y_t),
    model_series(m_EA, states_zero_mu, obs_zero_mu, pseudo_zero_mu, :y_t);
    anchor_idx = anchor_idx,
)
output_zero_cy_ex_mu_raw = apply_counterfactual_tail(
    output_base_raw,
    model_series(m_EA, states_actual, obs_actual, pseudo_actual, :y_t),
    model_series(m_EA, states_zero_cy_zero_mu, obs_zero_cy_zero_mu, pseudo_zero_cy_zero_mu, :y_t);
    anchor_idx = anchor_idx,
)
rstar_zero_cy_raw = apply_counterfactual_tail(
    rstar_base_raw,
    model_series(m_EA, states_actual, obs_actual, pseudo_actual, :Forward5YearRealNaturalRate),
    model_series(m_EA, states_zero_cy, obs_zero_cy, pseudo_zero_cy, :Forward5YearRealNaturalRate);
    anchor_idx = anchor_idx,
)

##########################################################################################
## CONVERT TO CASE2B PLOTTED UNITS
##########################################################################################

cy_base = 4 .* cy_base_raw
cy_zero_cy = 4 .* cy_zero_cy_raw

long_base = 4 .* long_base_raw
long_zero_cy = 4 .* long_zero_cy_raw

policy_base = 4 .* policy_base_raw
policy_zero_cy = 4 .* policy_zero_cy_raw

infl_base = four_quarter_sum(infl_base_raw)
infl_zero_cy = four_quarter_sum(infl_zero_cy_raw)

# Shift the baseline output path to have zero mean over the plotted sample.
output_level_shift = -mean(output_base_ex_mu_raw[mask_plot])
output_base = output_base_ex_mu_raw .+ output_level_shift
output_zero_cy = output_zero_cy_ex_mu_raw .+ output_level_shift

rstar_base = 4 .* rstar_base_raw
rstar_zero_cy = 4 .* rstar_zero_cy_raw

##########################################################################################
## PLOT EA CASE2B CHART
##########################################################################################

first_tick_year = Dates.year(anchor_date)
last_tick_year = Dates.year(dates_EA[end])
tick_step = 5
year_tick_years = collect(first_tick_year:tick_step:last_tick_year)
year_tick_dates = [quartertodate("$(y)-Q1") for y in year_tick_years]
year_tick_labels = string.(year_tick_years)
xticks_common = (year_tick_dates, year_tick_labels)

plot_dates = dates_EA[mask_plot]

p_cy = make_panel(
    plot_dates,
    cy_base[mask_plot],
    cy_zero_cy[mask_plot],
    "EA convenience yield (APR)",
    xticks_common;
    show_legend = true,
)
p_long = make_panel(plot_dates, long_base[mask_plot], long_zero_cy[mask_plot], "EA long-term interest rate (APR)", xticks_common)
p_pol = make_panel(plot_dates, policy_base[mask_plot], policy_zero_cy[mask_plot], "EA policy rate (APR)", xticks_common)
p_inf = make_panel(plot_dates, infl_base[mask_plot], infl_zero_cy[mask_plot], "EA inflation (%, yoy)", xticks_common)
p_out = make_panel(plot_dates, output_base[mask_plot], output_zero_cy[mask_plot], "EA output deviation from trend (%)", xticks_common)
p_rstar = make_panel(plot_dates, rstar_base[mask_plot], rstar_zero_cy[mask_plot], "EA r* (APR)", xticks_common)

p_3x2 = plot(
    p_cy, p_long,
    p_pol, p_inf,
    p_out, p_rstar,
    layout = (3, 2),
    size = (1100, 1200),
)

display(p_3x2)

if use_FG_in_EA
    savefig(p_3x2, joinpath(figroot, "EA Case2B no exorbitant privilege 3x2.pdf"))
    savefig(p_3x2, joinpath(figroot, "EA Case2B no exorbitant privilege 3x2.png"))
else
    savefig(p_3x2, joinpath(figroot, "EA Case2B no exorbitant privilege without FG 3x2.pdf"))
    savefig(p_3x2, joinpath(figroot, "EA Case2B no exorbitant privilege without FG 3x2.png"))
end
