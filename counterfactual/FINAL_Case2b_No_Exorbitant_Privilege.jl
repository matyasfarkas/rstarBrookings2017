##########################################################################################
## SETUP
##########################################################################################
using DSGE, ClusterManagers, HDF5, Plots, StatsPlots
using DataFrames, CSV, Dates

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

first_panel_legend_pos = (0.22, 0.25) # relative position within the first panel (convenience yield)



##############
# Load the EA model with FG 
##############
use_FG_in_EA  = true  # set to false to load the EA model without FG shocks


# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss24")

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")


m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
if use_FG_in_EA 
    m <= DSGE.Setting(:data_vintage, "250115")
else
    m <= DSGE.Setting(:data_vintage, "250116")
end
m <= DSGE.Setting(:reoptimize, false)
m <= DSGE.Setting(:calculate_hessian, false)
m <= DSGE.Setting(:date_mainsample_start,  quartertodate("1970-Q3"))
m <= DSGE.Setting(:date_presample_start,  quartertodate("1970-Q2"))

# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q3"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))


df = load_data(m; check_empty_columns = false)
if use_FG_in_EA
mode_file = joinpath(saveroot,"output_data/m1010/ss24/estimate/raw/" ,  "paramsmode_vint=250115.h5")
else
mode_file = joinpath(saveroot, "output_data/m1010/ss24/estimate/raw/" ,  "paramsmode_vint=250116.h5")
end
specify_mode!(m, mode_file)
system = DSGE.compute_system(m)

states = Dict{Symbol, Matrix{Float64}}()
shocks = Dict{Symbol, Matrix{Float64}}()
pseudo = Dict{Symbol, Matrix{Float64}}()

shock_labels = [key for (key, _) in sort(collect(m.exogenous_shocks), by = x -> x[2])]


combined = DSGE.OrderedDict{Symbol, Int64}()
# First insert all entries from endogenous_states
for (k, v) in m.endogenous_states
combined[k] = v
end

# Then insert entries from endogenous_states_augmented
for (k, v) in m.endogenous_states_augmented
combined[k] = v
end
state_labels = [key for (key, _) in sort(collect(combined), by = x -> x[2])]
pseudo_labels = [key for (key, _) in sort(collect(m.pseudo_observables), by = x -> x[2])]

system = DSGE.compute_system(m)

states = Dict{Symbol, Matrix{Float64}}()
shocks = Dict{Symbol, Matrix{Float64}}()
pseudo = Dict{Symbol, Matrix{Float64}}()

shock_labels = [key for (key, _) in sort(collect(m.exogenous_shocks), by = x -> x[2])]


combined = DSGE.OrderedDict{Symbol, Int64}()
# First insert all entries from endogenous_states
for (k, v) in m.endogenous_states
combined[k] = v
end

# Then insert entries from endogenous_states_augmented
for (k, v) in m.endogenous_states_augmented
combined[k] = v
end
state_labels = [key for (key, _) in sort(collect(combined), by = x -> x[2])]
pseudo_labels = [key for (key, _) in sort(collect(m.pseudo_observables), by = x -> x[2])]

system = DSGE.compute_system(m)

states_df = Dict{Symbol, DataFrame}()
shocks_df = Dict{Symbol, DataFrame}()
pseudo_df = Dict{Symbol, DataFrame}()

smoother = :durbin_koopman #:hamilton, :koopman, :carter_kohn, 
m <= DSGE.Setting(:forecast_smoother, smoother)

states[smoother], shocks[smoother], pseudo[smoother] = DSGE.smooth(m, df, system; draw_states = false)

dates = df.date[end-size(states[smoother],2)+1:end]
mat = states[smoother]'  # transpose to 259×91
states_df[smoother] = DataFrame(hcat(dates, mat), [:date; state_labels])

mat = shocks[smoother]'  # transpose to 259×29
shocks_df[smoother] = DataFrame(hcat(dates, mat), [:date; shock_labels])

mat = pseudo[smoother]'  # transpose to 259×22
pseudo_df[smoother] = DataFrame(hcat(dates, mat), [:date; pseudo_labels])

# Collect from the dataframe the respective shocks

#shocks_df[smoother][:, [:date; :b_liqtil_sh; :b_liqp_sh; :b_safetil_sh; :b_safep_sh]]
privilege_shock_names = [:b_liqtil_sh; :b_liqp_sh; :b_safetil_sh; :b_safep_sh]
privilege_shock_vals  = convert(Matrix,shocks_df[smoother][:, [:b_liqtil_sh; :b_liqp_sh; :b_safetil_sh; :b_safep_sh]])
 
horizon= size(privilege_shock_vals,1)
nshocks = size(system[:RRR], 2)
nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)
shocks = zeros(nshocks, horizon)

m1 = Model1010("ss20")
mode_file = joinpath(saveroot,"output_data/m1010/ss20/estimate/raw" ,  "paramsmode_vint=250825.h5")
specify_mode!(m1, mode_file)
system_US = DSGE.compute_system(m1)

for (i, shock_name) in enumerate(privilege_shock_names)
    shock_ind_US = m1.exogenous_shocks[shock_name]
    shocks[shock_ind_US, :] .= privilege_shock_vals[:, i]
end

# Compute IRF for a unit shock at time t for the specified shock
        states_privilege, obs_privilege, pseudo_privilege = forecast(system, s_0, shocks)
               
plotstart = 114 # 1999-Q1
horizon = size(privilege_shock_vals,1) - plotstart + 1
using Plots
p1 = plot(dates[plotstart:end],obs_privilege[m.observables[:obs_nominalrate],plotstart:end],color=:blue, lw=2,title="Policy rate")
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(dates[plotstart:end],obs_privilege[m.observables[:obs_corepce],plotstart:end],color=:blue, lw=2,title="Inflation")
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(dates[plotstart:end],states_privilege[m.endogenous_states[:y_t],plotstart:end],color=:blue, lw=2 , title="Output")#
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(dates[plotstart:end],pseudo_privilege[m.pseudo_observables[:Forward5YearRealNaturalRate],plotstart:end],color=:blue, lw=2,title="r* (Forward 5-year real natural rate)")#
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
# p6 = plot(1:horizon,pseudo[m.pseudo_observables[:RealNaturalRate],:, m.exogenous_shocks[:rm_sh]],title="Real natural rate")#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")


plot(p1, p2, p3, p4,layout=(2,2), legend=false)
plot!(size=(960,540))
if use_FG_in_EA 
    savefig( joinpath(saveroot, "Final Paper","Figures", "Exorbitant privilege.pdf") )   # saves the plot from p as a .pdf vector graphic
else
    savefig( joinpath(saveroot, "Final Paper","Figures","Exorbitant privilege without FG shocks in EA.pdf"))   # saves the plot from p as a .pdf vector graphic
end



# --- Write plotted series to CSV ---
using CSV, DataFrames
df_out = DataFrame(
    Date = dates[plotstart:end],
    PolicyRate = obs_privilege[m.observables[:obs_nominalrate], plotstart:end],
    Inflation = obs_privilege[m.observables[:obs_corepce], plotstart:end],
    Output = states_privilege[m.endogenous_states[:y_t], plotstart:end],
    Forward5YearRealNaturalRate = pseudo_privilege[m.pseudo_observables[:Forward5YearRealNaturalRate], plotstart:end]
)
if use_FG_in_EA
    CSV.write("Exorbitant_privilege.csv", df_out)
else
    CSV.write("Exorbitant_privilege_without_FG_shocks_in_EA.csv", df_out)
end

system10 = compute_system(m)
horzion = 40
states_irf10, obs_irf10, pseudo_irf10 = impulse_responses(system10, horizon)

# output
p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]],title="Permanent liquidity shock", label=["Basline model"])
plot!(legend=:bottomright)

p2 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_safep_sh]]] ,title="Permanent safety shock", label=["Basline model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:zp_sh ]] ] ,title="Permanent technology shock", label=["Basline model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:rm_shl6 ]] ] ,title="FG6 shock", label=["Basline model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)


# 
p1 = plot(1:horizon,-[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]./ minimum(obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]])],title="Output", label=["Basline model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,-[obs_irf10[m.observables[:obs_corepce],:, m.exogenous_shocks[:b_liqp_sh]]./ minimum(obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]])] ,title="Inflation", label=["Basline model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon, -[obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]] ./ minimum(obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]])], title="Policy rate", label=["Basline model"])
plot!(legend=:bottomright)
p4=  plot(1:horizon,-[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh ]]./ minimum(obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]]) ] ,title="r* (Forward 5-year real natural rate)", label=["Basline model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)

    savefig( joinpath(saveroot, "Final Paper","Figures", "IRF_to_permanet_liquidity_shock_scaled.pdf") )       # saves the plot from p as a .pdf vector graphic

##########################################################################################
## APPENDED: 3x2 PANEL WITH CONVENIENCE YIELD + LONG RATE
##########################################################################################

# -------------------------
# Helper functions
# -------------------------
function _get_series(m, states_mat, obs_mat, pseudo_mat, sym::Symbol)
    if haskey(m.endogenous_states, sym)
        return vec(states_mat[m.endogenous_states[sym], :])
    elseif haskey(m.endogenous_states_augmented, sym)
        return vec(states_mat[m.endogenous_states_augmented[sym], :])
    elseif haskey(m.observables, sym)
        return vec(obs_mat[m.observables[sym], :])
    elseif haskey(m.pseudo_observables, sym)
        return vec(pseudo_mat[m.pseudo_observables[sym], :])
    else
        error("Series $(sym) not found in states / observables / pseudo-observables.")
    end
end

function _find_first_symbol(m, candidates::Vector{Symbol})
    for sym in candidates
        if haskey(m.endogenous_states, sym) ||
           haskey(m.endogenous_states_augmented, sym) ||
           haskey(m.observables, sym) ||
           haskey(m.pseudo_observables, sym)
            return sym
        end
    end
    return nothing
end

const convenience_yield_state_names = [:b_liqtil_t, :b_liqp_t, :b_safetil_t, :b_safep_t]

_has_state_series(m, sym::Symbol) =
    haskey(m.endogenous_states, sym) || haskey(m.endogenous_states_augmented, sym)

function _get_convenience_yield_drift(m)
    return 100 * log(m[:lnb_liq]) + 100 * log(m[:lnb_safe])
end

function _sum_state_series(m, states_mat, obs_mat, pseudo_mat, syms::Vector{Symbol})
    sum_series = _get_series(m, states_mat, obs_mat, pseudo_mat, syms[1])
    for sym in syms[2:end]
        sum_series .+= _get_series(m, states_mat, obs_mat, pseudo_mat, sym)
    end
    return sum_series
end

function _get_convenience_yield_series(m, states_mat, obs_mat, pseudo_mat)
    if _has_state_series(m, :b_liq_t) && _has_state_series(m, :b_safe_t)
        cy_series = _get_series(m, states_mat, obs_mat, pseudo_mat, :b_liq_t)
        cy_series .+= _get_series(m, states_mat, obs_mat, pseudo_mat, :b_safe_t)
    else
        missing_states = [
            sym for sym in convenience_yield_state_names if !_has_state_series(m, sym)
        ]
        isempty(missing_states) || error("Convenience-yield state(s) not found in model: $(missing_states)")

        cy_series = _sum_state_series(m, states_mat, obs_mat, pseudo_mat, convenience_yield_state_names)
    end
    return cy_series .+ _get_convenience_yield_drift(m)
end

function _get_long_rate_series(m, states_mat, obs_mat, pseudo_mat)
    sym = _resolve_long_rate(m)
    return _get_series(m, states_mat, obs_mat, pseudo_mat, sym)
end

function _resolve_long_rate(m)
    candidates = [
        :obs_longrate, :obs_long_rate, :longrate, :LongRate,
        :LongTermRate, :LongTermNominalRate, :tenyearrate, :obs_tenyearrate
    ]
    sym = _find_first_symbol(m, candidates)
    sym === nothing && error("""
Could not locate the long-term interest rate series automatically.
Please update candidate names in _resolve_long_rate(...) to match your model.
""")
    return sym
end

function _get_baseline_series(m, df_aligned::DataFrame, states_sm, pseudo_sm, sym::Symbol)
    if haskey(m.observables, sym)
        sym in names(df_aligned) || error("Observable $(sym) not found in df_aligned.")
        return vec(df_aligned[!, sym])
    elseif haskey(m.endogenous_states, sym)
        return vec(states_sm[m.endogenous_states[sym], :])
    elseif haskey(m.endogenous_states_augmented, sym)
        return vec(states_sm[m.endogenous_states_augmented[sym], :])
    elseif haskey(m.pseudo_observables, sym)
        return vec(pseudo_sm[m.pseudo_observables[sym], :])
    else
        error("Series $(sym) not found in model.")
    end
end

function _replace_shock_block!(
    shocks_mat::AbstractMatrix,
    m_US,
    shock_names::Vector{Symbol},
    shock_vals::AbstractMatrix
)
    size(shocks_mat, 2) == size(shock_vals, 1) || error("Shock replacement horizon mismatch.")
    size(shock_vals, 2) == length(shock_names) || error("Shock replacement width mismatch.")

    for (i, shock_name) in enumerate(shock_names)
        shock_ind_US = m_US.exogenous_shocks[shock_name]
        shocks_mat[shock_ind_US, :] .= shock_vals[:, i]
    end

    return shocks_mat
end

function _forecast_us_from_anchor(system_US, anchor_state::AbstractVector, shocks_tail::AbstractMatrix)
    return forecast(system_US, collect(anchor_state), shocks_tail)
end

function _apply_counterfactual_tail(
    baseline::AbstractVector,
    actual_tail::AbstractVector,
    cf_tail::AbstractVector;
    anchor_idx::Int = 1
)
    length(actual_tail) == length(cf_tail) || error("Tail series length mismatch.")
    length(baseline) == anchor_idx + length(actual_tail) || error("Baseline/tail length mismatch.")

    out = collect(baseline)
    out[anchor_idx+1:end] .= baseline[anchor_idx+1:end] .+ (cf_tail .- actual_tail)
    return out
end

# -------------------------
# Load and smooth the US baseline model for the black convenience-yield line
# -------------------------
m1 <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m1 <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m1 <= DSGE.Setting(:data_vintage, "250825")
m1 <= DSGE.Setting(:reoptimize, false)
m1 <= DSGE.Setting(:calculate_hessian, false)
m1 <= DSGE.Setting(:date_mainsample_start, quartertodate("1970-Q3"))
m1 <= DSGE.Setting(:date_presample_start, quartertodate("1970-Q2"))
m1 <= DSGE.Setting(:date_forecast_start, quartertodate("2024-Q3"))
m1 <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))
m1 <= DSGE.Setting(:forecast_smoother, smoother)

df_US = load_data(m1; check_empty_columns = false)
states_US_sm, shocks_US_sm, pseudo_US_sm = DSGE.smooth(m1, df_US, system_US; draw_states = false)
df_US_aligned = df_US[end-size(states_US_sm, 2)+1:end, :]
dates_US = df_US.date[end-size(states_US_sm, 2)+1:end]

# -------------------------
# Build the additional plotted series
# -------------------------
# Align sample start with your existing plotstart
plot_start_date = dates[plotstart]
mask_CF = dates .>= plot_start_date
mask_US = dates_US .>= plot_start_date
common_dates = intersect(dates[mask_CF], dates_US[mask_US])
isempty(common_dates) && error("No overlapping dates between EA and US samples in the earlier 3x2 block.")

idx_CF = findall(in(common_dates), dates)
idx_US = findall(in(common_dates), dates_US)

dates_plot = dates[idx_CF]
dates_plot == dates_US[idx_US] || error("Earlier 3x2 block date alignment failed: date order mismatch.")

df_US_common = df_US_aligned[idx_US, :]
states_US_common = states_US_sm[:, idx_US]
pseudo_US_common = pseudo_US_sm[:, idx_US]
privilege_shock_vals_EA = convert(Matrix, shocks_df[smoother][idx_CF, privilege_shock_names])

longrate_sym = _resolve_long_rate(m1)

cy_US = _get_convenience_yield_series(m1, states_US_common, nothing, pseudo_US_common)
longrate_US = _get_baseline_series(m1, df_US_common, states_US_common, pseudo_US_common, longrate_sym)
policy_US = _get_baseline_series(m1, df_US_common, states_US_common, pseudo_US_common, :obs_nominalrate)
inflation_US = _get_baseline_series(m1, df_US_common, states_US_common, pseudo_US_common, :obs_corepce)
output_US = _get_baseline_series(m1, df_US_common, states_US_common, pseudo_US_common, :y_t)
rstar_US = _get_baseline_series(m1, df_US_common, states_US_common, pseudo_US_common, :Forward5YearRealNaturalRate)

anchor_state_US = states_US_common[:, 1]
us_shocks_tail = copy(shocks_US_sm[:, idx_US[2:end]])
zero_privilege_tail = _replace_shock_block!(
    copy(us_shocks_tail),
    m1,
    privilege_shock_names,
    zeros(length(idx_US) - 1, length(privilege_shock_names))
)
ea_privilege_tail = _replace_shock_block!(
    copy(us_shocks_tail),
    m1,
    privilege_shock_names,
    privilege_shock_vals_EA[2:end, :]
)

states_actual_US, obs_actual_US, pseudo_actual_US = _forecast_us_from_anchor(system_US, anchor_state_US, us_shocks_tail)
states_remove_US, obs_remove_US, pseudo_remove_US = _forecast_us_from_anchor(system_US, anchor_state_US, zero_privilege_tail)
states_cf_US, obs_cf_US, pseudo_cf_US = _forecast_us_from_anchor(system_US, anchor_state_US, ea_privilege_tail)

cy_actual_tail = _get_convenience_yield_series(m1, states_actual_US, obs_actual_US, pseudo_actual_US)
cy_remove_tail = _get_convenience_yield_series(m1, states_remove_US, obs_remove_US, pseudo_remove_US)
cy_cf_tail = _get_convenience_yield_series(m1, states_cf_US, obs_cf_US, pseudo_cf_US)
cy_remove_US = _apply_counterfactual_tail(cy_US, cy_actual_tail, cy_remove_tail)
cy_CF = _apply_counterfactual_tail(cy_US, cy_actual_tail, cy_cf_tail)

longrate_actual_tail = _get_series(m1, states_actual_US, obs_actual_US, pseudo_actual_US, longrate_sym)
longrate_remove_tail = _get_series(m1, states_remove_US, obs_remove_US, pseudo_remove_US, longrate_sym)
longrate_cf_tail = _get_series(m1, states_cf_US, obs_cf_US, pseudo_cf_US, longrate_sym)
longrate_remove_US = _apply_counterfactual_tail(longrate_US, longrate_actual_tail, longrate_remove_tail)
longrate_CF = _apply_counterfactual_tail(longrate_US, longrate_actual_tail, longrate_cf_tail)

policy_actual_tail = _get_series(m1, states_actual_US, obs_actual_US, pseudo_actual_US, :obs_nominalrate)
policy_remove_tail = _get_series(m1, states_remove_US, obs_remove_US, pseudo_remove_US, :obs_nominalrate)
policy_cf_tail = _get_series(m1, states_cf_US, obs_cf_US, pseudo_cf_US, :obs_nominalrate)
policy_remove_US = _apply_counterfactual_tail(policy_US, policy_actual_tail, policy_remove_tail)
policy_CF = _apply_counterfactual_tail(policy_US, policy_actual_tail, policy_cf_tail)

inflation_actual_tail = _get_series(m1, states_actual_US, obs_actual_US, pseudo_actual_US, :obs_corepce)
inflation_remove_tail = _get_series(m1, states_remove_US, obs_remove_US, pseudo_remove_US, :obs_corepce)
inflation_cf_tail = _get_series(m1, states_cf_US, obs_cf_US, pseudo_cf_US, :obs_corepce)
inflation_remove_US = _apply_counterfactual_tail(inflation_US, inflation_actual_tail, inflation_remove_tail)
inflation_CF = _apply_counterfactual_tail(inflation_US, inflation_actual_tail, inflation_cf_tail)

output_actual_tail = _get_series(m1, states_actual_US, obs_actual_US, pseudo_actual_US, :y_t)
output_remove_tail = _get_series(m1, states_remove_US, obs_remove_US, pseudo_remove_US, :y_t)
output_cf_tail = _get_series(m1, states_cf_US, obs_cf_US, pseudo_cf_US, :y_t)
output_remove_US = _apply_counterfactual_tail(output_US, output_actual_tail, output_remove_tail)
output_CF = _apply_counterfactual_tail(output_US, output_actual_tail, output_cf_tail)

rstar_actual_tail = _get_series(m1, states_actual_US, obs_actual_US, pseudo_actual_US, :Forward5YearRealNaturalRate)
rstar_remove_tail = _get_series(m1, states_remove_US, obs_remove_US, pseudo_remove_US, :Forward5YearRealNaturalRate)
rstar_cf_tail = _get_series(m1, states_cf_US, obs_cf_US, pseudo_cf_US, :Forward5YearRealNaturalRate)
rstar_remove_US = _apply_counterfactual_tail(rstar_US, rstar_actual_tail, rstar_remove_tail)
rstar_CF = _apply_counterfactual_tail(rstar_US, rstar_actual_tail, rstar_cf_tail)

# -------------------------
# X-axis ticks: show years only
# -------------------------
first_tick_year = Dates.year(plot_start_date)
last_tick_year  = Dates.year(dates[end])

# If yearly labels are too crowded, change 1 to 2 or 4
tick_step = 5

year_tick_years  = collect(first_tick_year:tick_step:last_tick_year)
year_tick_dates  = [quartertodate("$(y)-Q1") for y in year_tick_years]
year_tick_labels = string.(year_tick_years)

# -------------------------
# Build 3x2 panel
# -------------------------
p_cy = plot(
    dates_plot, cy_US,
    color = :black, lw = 2, linestyle = :solid,
    title = "Convenience yield",
    label = "US - Actual",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_cy, dates_plot, cy_remove_US, color = :red, lw = 2, linestyle = :dashdot, label = "US excluding matched privilege shocks")
plot!(p_cy, dates_plot, zeros(length(dates_plot)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

plot!(
    p_cy,
    dates_plot, cy_CF,
    color = :blue, lw = 2,
    label = "Counterfactual - EA CY innovations",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_cy, legend = first_panel_legend_pos)

p_long = plot(
    dates_plot, longrate_US,
    color = :black, lw = 2, linestyle = :solid,
    title = "Long-term interest rate",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_long, dates_plot, longrate_remove_US, color = :red, lw = 2, linestyle = :dashdot, label = "")
plot!(p_long, dates_plot, longrate_CF, color = :blue, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))
plot!(p_long, dates_plot, zeros(length(dates_plot)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_pol = plot(
    dates_plot, policy_US,
    color = :black, lw = 2, linestyle = :solid,
    title = "Policy rate",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_pol, dates_plot, policy_remove_US, color = :red, lw = 2, linestyle = :dashdot, label = "")
plot!(p_pol, dates_plot, policy_CF, color = :blue, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))
plot!(p_pol, dates_plot, zeros(length(dates_plot)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_inf = plot(
    dates_plot, inflation_US,
    color = :black, lw = 2, linestyle = :solid,
    title = "Inflation",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_inf, dates_plot, inflation_remove_US, color = :red, lw = 2, linestyle = :dashdot, label = "")
plot!(p_inf, dates_plot, inflation_CF, color = :blue, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))
plot!(p_inf, dates_plot, zeros(length(dates_plot)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_out = plot(
    dates_plot, output_US,
    color = :black, lw = 2, linestyle = :solid,
    title = "Output",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_out, dates_plot, output_remove_US, color = :red, lw = 2, linestyle = :dashdot, label = "")
plot!(p_out, dates_plot, output_CF, color = :blue, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))
plot!(p_out, dates_plot, zeros(length(dates_plot)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_rstar = plot(
    dates_plot, rstar_US,
    color = :black, lw = 2, linestyle = :solid,
    title = "r* (Forward 5-year real natural rate)",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_rstar, dates_plot, rstar_remove_US, color = :red, lw = 2, linestyle = :dashdot, label = "")
plot!(p_rstar, dates_plot, rstar_CF, color = :blue, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))
plot!(p_rstar, dates_plot, zeros(length(dates_plot)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_3x2 = plot(
    p_cy, p_long,
    p_pol, p_inf,
    p_out, p_rstar,
    layout = (3, 2),
    size = (1100, 1200)
)

display(p_3x2)

if use_FG_in_EA
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege 3x2.pdf"))
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege 3x2.png"))
else
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege without FG shocks in EA 3x2.pdf"))
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege without FG shocks in EA 3x2.png"))
end

# -------------------------
# Export the plotted series
# -------------------------
df_plot_3x2 = DataFrame(
    Date = dates_plot,

    ConvenienceYield_Baseline = cy_US,
    ConvenienceYield_USExcludingMatchedPrivilege = cy_remove_US,
    ConvenienceYield_Counterfactual = cy_CF,

    LongTermRate_Baseline = longrate_US,
    LongTermRate_USExcludingMatchedPrivilege = longrate_remove_US,
    LongTermRate_Counterfactual = longrate_CF,

    PolicyRate_Baseline = policy_US,
    PolicyRate_USExcludingMatchedPrivilege = policy_remove_US,
    PolicyRate_Counterfactual = policy_CF,

    Inflation_Baseline = inflation_US,
    Inflation_USExcludingMatchedPrivilege = inflation_remove_US,
    Inflation_Counterfactual = inflation_CF,

    Output_Baseline = output_US,
    Output_USExcludingMatchedPrivilege = output_remove_US,
    Output_Counterfactual = output_CF,

    Forward5YearRealNaturalRate_Baseline = rstar_US,
    Forward5YearRealNaturalRate_USExcludingMatchedPrivilege = rstar_remove_US,
    Forward5YearRealNaturalRate_Counterfactual = rstar_CF
)

if use_FG_in_EA
    CSV.write(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant_privilege_3x2.csv"), df_plot_3x2)
else
    CSV.write(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant_privilege_without_FG_shocks_in_EA_3x2.csv"), df_plot_3x2)
end
##########################################################################################
## APPENDED: REWORKED 3x2 PANEL WITH COUNTERFACTUAL US PATHS
## - Common-date alignment between EA and US samples
## - Black solid = US baseline
## - Red dash-dot = US path with privilege shocks set to zero from 1999Q1 onward
## - Blue dashed = US path with EA privilege shocks from 1999Q1 onward
## - Counterfactual paths share the observed 1998Q4 starting point
## - Inflation = obs_corepce shown as y/y (4-quarter sum)
## - Convenience yield, LT rate, policy rate, r* shown as APR (*4)
##########################################################################################

# -------------------------
# Helper functions (unique names to avoid collisions)
# -------------------------
function _ep_find_first_symbol(m, candidates::Vector{Symbol})
    for sym in candidates
        if haskey(m.endogenous_states, sym) ||
           haskey(m.endogenous_states_augmented, sym) ||
           haskey(m.observables, sym) ||
           haskey(m.pseudo_observables, sym)
            return sym
        end
    end
    return nothing
end

function _ep_get_baseline_series(m, df_aligned::DataFrame, states_sm, pseudo_sm, sym::Symbol)
    if haskey(m.observables, sym)
        sym in names(df_aligned) || error("Observable $(sym) not found in df_aligned.")
        return vec(df_aligned[!, sym])
    elseif haskey(m.endogenous_states, sym)
        return vec(states_sm[m.endogenous_states[sym], :])
    elseif haskey(m.endogenous_states_augmented, sym)
        return vec(states_sm[m.endogenous_states_augmented[sym], :])
    elseif haskey(m.pseudo_observables, sym)
        return vec(pseudo_sm[m.pseudo_observables[sym], :])
    else
        error("Series $(sym) not found in model.")
    end
end

function _ep_get_delta_series(m, states_delta, obs_delta, pseudo_delta, sym::Symbol)
    if haskey(m.observables, sym)
        return vec(obs_delta[m.observables[sym], :])
    elseif haskey(m.endogenous_states, sym)
        return vec(states_delta[m.endogenous_states[sym], :])
    elseif haskey(m.endogenous_states_augmented, sym)
        return vec(states_delta[m.endogenous_states_augmented[sym], :])
    elseif haskey(m.pseudo_observables, sym)
        return vec(pseudo_delta[m.pseudo_observables[sym], :])
    else
        error("Series $(sym) not found in model.")
    end
end

function _ep_resolve_long_rate(m)
    candidates = [
        :obs_longrate, :obs_long_rate, :longrate, :LongRate,
        :LongTermRate, :LongTermNominalRate, :tenyearrate, :obs_tenyearrate
    ]
    sym = _ep_find_first_symbol(m, candidates)
    sym === nothing && error("""
Could not locate the long-term interest rate series automatically.
Please update candidate names in _ep_resolve_long_rate(...) to match your model.
""")
    return sym
end

function _ep_four_quarter_sum(x::AbstractVector)
    y = fill(NaN, length(x))
    for t in 4:length(x)
        y[t] = x[t] + x[t-1] + x[t-2] + x[t-3]
    end
    return y
end

# -------------------------
# Smooth the US baseline model
# -------------------------
m1 <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m1 <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m1 <= DSGE.Setting(:data_vintage, "250825")
m1 <= DSGE.Setting(:reoptimize, false)
m1 <= DSGE.Setting(:calculate_hessian, false)
m1 <= DSGE.Setting(:date_mainsample_start, quartertodate("1970-Q3"))
m1 <= DSGE.Setting(:date_presample_start, quartertodate("1970-Q2"))
m1 <= DSGE.Setting(:date_forecast_start, quartertodate("2024-Q3"))
m1 <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))
m1 <= DSGE.Setting(:forecast_smoother, smoother)

df_US_ep = load_data(m1; check_empty_columns = false)
states_US_sm_ep, shocks_US_sm_ep, pseudo_US_sm_ep = DSGE.smooth(m1, df_US_ep, system_US; draw_states = false)

T_US_ep = size(states_US_sm_ep, 2)
df_US_aligned_ep = df_US_ep[end-T_US_ep+1:end, :]
dates_US_ep = df_US_aligned_ep.date

# -------------------------
# Align EA and US on COMMON DATES
# -------------------------
common_dates_ep = intersect(dates, dates_US_ep)
isempty(common_dates_ep) && error("No overlapping dates between EA and US samples.")

idx_EA_ep = findall(in(common_dates_ep), dates)
idx_US_ep = findall(in(common_dates_ep), dates_US_ep)

dates_common_EA_ep = dates[idx_EA_ep]
dates_common_US_ep = dates_US_ep[idx_US_ep]
dates_common_EA_ep == dates_common_US_ep || error("Common-date extraction failed: date order mismatch.")

dates_plot_ep = dates_common_EA_ep

# -------------------------
# Collect EA privilege shocks on the common sample
# -------------------------
privilege_shock_vals_EA_ep = convert(Matrix, shocks_df[smoother][idx_EA_ep, privilege_shock_names])
# -------------------------
# Resolve plotted variables
# -------------------------
longrate_sym_ep = _ep_resolve_long_rate(m1)

# -------------------------
# Start date and anchor normalization date: 1998Q4
# -------------------------
plot_start_date_ep = quartertodate("1998-Q4")
anchor_idx_ep = findfirst(==(plot_start_date_ep), dates_plot_ep)
anchor_idx_ep === nothing && error("Anchor date 1998Q4 not found in dates_plot_ep.")
anchor_idx_ep < 4 && error("Need at least 3 quarters before 1998Q4 to construct y/y inflation.")

mask_plot_ep = dates_plot_ep .>= plot_start_date_ep

# -------------------------
# Align US baseline objects to common dates
# -------------------------
df_US_common_ep = df_US_aligned_ep[idx_US_ep, :]
states_US_common_ep = states_US_sm_ep[:, idx_US_ep]
pseudo_US_common_ep = pseudo_US_sm_ep[:, idx_US_ep]

# -------------------------
# Build RAW baseline and counterfactual series first
# -------------------------

# Convenience yield
cy_base_raw_ep = _get_convenience_yield_series(m1, states_US_common_ep, nothing, pseudo_US_common_ep)

# Long rate
long_base_raw_ep = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, longrate_sym_ep)

# Policy rate
policy_base_raw_ep = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :obs_nominalrate)

# Inflation (quarterly obs_corepce)
infl_base_raw_ep = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :obs_corepce)

# Output gap
output_base_raw_ep = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :y_t)

# r*
rstar_base_raw_ep = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :Forward5YearRealNaturalRate)

anchor_state_US_ep = states_US_common_ep[:, anchor_idx_ep]
us_shocks_tail_ep = copy(shocks_US_sm_ep[:, idx_US_ep[anchor_idx_ep+1:end]])
zero_privilege_tail_ep = _replace_shock_block!(
    copy(us_shocks_tail_ep),
    m1,
    privilege_shock_names,
    zeros(length(idx_US_ep) - anchor_idx_ep, length(privilege_shock_names))
)
ea_privilege_tail_ep = _replace_shock_block!(
    copy(us_shocks_tail_ep),
    m1,
    privilege_shock_names,
    privilege_shock_vals_EA_ep[anchor_idx_ep+1:end, :]
)

states_actual_ep, obs_actual_ep, pseudo_actual_ep = _forecast_us_from_anchor(system_US, anchor_state_US_ep, us_shocks_tail_ep)
states_remove_US_ep, obs_remove_US_ep, pseudo_remove_US_ep = _forecast_us_from_anchor(system_US, anchor_state_US_ep, zero_privilege_tail_ep)
states_cf_ep, obs_cf_ep, pseudo_cf_ep             = _forecast_us_from_anchor(system_US, anchor_state_US_ep, ea_privilege_tail_ep)

cy_actual_raw_tail_ep = _get_convenience_yield_series(m1, states_actual_ep, obs_actual_ep, pseudo_actual_ep)
cy_remove_raw_tail_ep = _get_convenience_yield_series(m1, states_remove_US_ep, obs_remove_US_ep, pseudo_remove_US_ep)
cy_cf_raw_tail_ep = _get_convenience_yield_series(m1, states_cf_ep, obs_cf_ep, pseudo_cf_ep)

cy_remove_raw_ep = _apply_counterfactual_tail(
    cy_base_raw_ep,
    cy_actual_raw_tail_ep,
    cy_remove_raw_tail_ep;
    anchor_idx = anchor_idx_ep
)
cy_cf_raw_ep = _apply_counterfactual_tail(
    cy_base_raw_ep,
    cy_actual_raw_tail_ep,
    cy_cf_raw_tail_ep;
    anchor_idx = anchor_idx_ep
)

long_actual_raw_tail_ep = _ep_get_delta_series(m1, states_actual_ep, obs_actual_ep, pseudo_actual_ep, longrate_sym_ep)
long_remove_raw_tail_ep = _ep_get_delta_series(m1, states_remove_US_ep, obs_remove_US_ep, pseudo_remove_US_ep, longrate_sym_ep)
long_cf_raw_tail_ep = _ep_get_delta_series(m1, states_cf_ep, obs_cf_ep, pseudo_cf_ep, longrate_sym_ep)
long_remove_raw_ep = _apply_counterfactual_tail(long_base_raw_ep, long_actual_raw_tail_ep, long_remove_raw_tail_ep; anchor_idx = anchor_idx_ep)
long_cf_raw_ep = _apply_counterfactual_tail(long_base_raw_ep, long_actual_raw_tail_ep, long_cf_raw_tail_ep; anchor_idx = anchor_idx_ep)

policy_actual_raw_tail_ep = _ep_get_delta_series(m1, states_actual_ep, obs_actual_ep, pseudo_actual_ep, :obs_nominalrate)
policy_remove_raw_tail_ep = _ep_get_delta_series(m1, states_remove_US_ep, obs_remove_US_ep, pseudo_remove_US_ep, :obs_nominalrate)
policy_cf_raw_tail_ep = _ep_get_delta_series(m1, states_cf_ep, obs_cf_ep, pseudo_cf_ep, :obs_nominalrate)
policy_remove_raw_ep = _apply_counterfactual_tail(policy_base_raw_ep, policy_actual_raw_tail_ep, policy_remove_raw_tail_ep; anchor_idx = anchor_idx_ep)
policy_cf_raw_ep = _apply_counterfactual_tail(policy_base_raw_ep, policy_actual_raw_tail_ep, policy_cf_raw_tail_ep; anchor_idx = anchor_idx_ep)

infl_actual_raw_tail_ep = _ep_get_delta_series(m1, states_actual_ep, obs_actual_ep, pseudo_actual_ep, :obs_corepce)
infl_remove_raw_tail_ep = _ep_get_delta_series(m1, states_remove_US_ep, obs_remove_US_ep, pseudo_remove_US_ep, :obs_corepce)
infl_cf_raw_tail_ep = _ep_get_delta_series(m1, states_cf_ep, obs_cf_ep, pseudo_cf_ep, :obs_corepce)
infl_remove_raw_ep = _apply_counterfactual_tail(infl_base_raw_ep, infl_actual_raw_tail_ep, infl_remove_raw_tail_ep; anchor_idx = anchor_idx_ep)
infl_cf_raw_ep = _apply_counterfactual_tail(infl_base_raw_ep, infl_actual_raw_tail_ep, infl_cf_raw_tail_ep; anchor_idx = anchor_idx_ep)

output_actual_raw_tail_ep = _ep_get_delta_series(m1, states_actual_ep, obs_actual_ep, pseudo_actual_ep, :y_t)
output_remove_raw_tail_ep = _ep_get_delta_series(m1, states_remove_US_ep, obs_remove_US_ep, pseudo_remove_US_ep, :y_t)
output_cf_raw_tail_ep = _ep_get_delta_series(m1, states_cf_ep, obs_cf_ep, pseudo_cf_ep, :y_t)
output_remove_raw_ep = _apply_counterfactual_tail(output_base_raw_ep, output_actual_raw_tail_ep, output_remove_raw_tail_ep; anchor_idx = anchor_idx_ep)
output_cf_raw_ep = _apply_counterfactual_tail(output_base_raw_ep, output_actual_raw_tail_ep, output_cf_raw_tail_ep; anchor_idx = anchor_idx_ep)

rstar_actual_raw_tail_ep = _ep_get_delta_series(m1, states_actual_ep, obs_actual_ep, pseudo_actual_ep, :Forward5YearRealNaturalRate)
rstar_remove_raw_tail_ep = _ep_get_delta_series(m1, states_remove_US_ep, obs_remove_US_ep, pseudo_remove_US_ep, :Forward5YearRealNaturalRate)
rstar_cf_raw_tail_ep = _ep_get_delta_series(m1, states_cf_ep, obs_cf_ep, pseudo_cf_ep, :Forward5YearRealNaturalRate)
rstar_remove_raw_ep = _apply_counterfactual_tail(rstar_base_raw_ep, rstar_actual_raw_tail_ep, rstar_remove_raw_tail_ep; anchor_idx = anchor_idx_ep)
rstar_cf_raw_ep = _apply_counterfactual_tail(rstar_base_raw_ep, rstar_actual_raw_tail_ep, rstar_cf_raw_tail_ep; anchor_idx = anchor_idx_ep)

# -------------------------
# Convert to plotted units
# -------------------------
# APR variables: *4
cy_base_ep = 4 .* cy_base_raw_ep
cy_remove_cf_ep = 4 .* cy_remove_raw_ep
cy_cf_ep = 4 .* cy_cf_raw_ep

long_base_ep = 4 .* long_base_raw_ep
long_remove_cf_ep = 4 .* long_remove_raw_ep
long_cf_ep = 4 .* long_cf_raw_ep

policy_base_ep = 4 .* policy_base_raw_ep
policy_remove_cf_ep = 4 .* policy_remove_raw_ep
policy_cf_ep = 4 .* policy_cf_raw_ep

rstar_base_ep = 4 .* rstar_base_raw_ep
rstar_remove_cf_ep = 4 .* rstar_remove_raw_ep
rstar_cf_ep = 4 .* rstar_cf_raw_ep

# Inflation: year-over-year = 4-quarter sum of quarterly inflation
infl_base_ep = _ep_four_quarter_sum(infl_base_raw_ep)
infl_remove_cf_ep = _ep_four_quarter_sum(infl_remove_raw_ep)
infl_cf_ep = _ep_four_quarter_sum(infl_cf_raw_ep)

# Output gap: leave as-is
output_base_ep = output_base_raw_ep
output_remove_cf_ep = output_remove_raw_ep
output_cf_ep = output_cf_raw_ep

# -------------------------
# X-axis ticks: show years only
# -------------------------
first_tick_year_ep = Dates.year(plot_start_date_ep)
last_tick_year_ep  = Dates.year(dates_plot_ep[end])
tick_step_ep = 5

year_tick_years_ep  = collect(first_tick_year_ep:tick_step_ep:last_tick_year_ep)
year_tick_dates_ep  = [quartertodate("$(y)-Q1") for y in year_tick_years_ep]
year_tick_labels_ep = string.(year_tick_years_ep)

# -------------------------
# Build 3x2 panel
# -------------------------
p_cy_ep = plot(
    dates_plot_ep[mask_plot_ep], cy_base_ep[mask_plot_ep],
    color = :black, lw = 2, linestyle = :solid,
    title = "Convenience yield (APR)",
    label = "US - Actual",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_cy_ep,
    dates_plot_ep[mask_plot_ep], cy_remove_cf_ep[mask_plot_ep],
    color = :red, lw = 2, linestyle = :dashdot,
    label = "US - No CY innovations"
)
plot!(
    p_cy_ep,
    dates_plot_ep[mask_plot_ep], cy_cf_ep[mask_plot_ep],
    color = :blue, lw = 2, linestyle = :dash,
    label = "Counterfactual - EA CY innovations"
)
plot!(p_cy_ep, legend = first_panel_legend_pos)

p_long_ep = plot(
    dates_plot_ep[mask_plot_ep], long_base_ep[mask_plot_ep],
    color = :black, lw = 2, linestyle = :solid,
    title = "Long-term interest rate (APR)",
    label = "",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_long_ep,
    dates_plot_ep[mask_plot_ep], long_remove_cf_ep[mask_plot_ep],
    color = :red, lw = 2, linestyle = :dashdot,
    label = ""
)
plot!(
    p_long_ep,
    dates_plot_ep[mask_plot_ep], long_cf_ep[mask_plot_ep],
    color = :blue, lw = 2, linestyle = :dash,
    label = ""
)

p_pol_ep = plot(
    dates_plot_ep[mask_plot_ep], policy_base_ep[mask_plot_ep],
    color = :black, lw = 2, linestyle = :solid,
    title = "Policy rate (APR)",
    label = "",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_pol_ep,
    dates_plot_ep[mask_plot_ep], policy_remove_cf_ep[mask_plot_ep],
    color = :red, lw = 2, linestyle = :dashdot,
    label = ""
)
plot!(
    p_pol_ep,
    dates_plot_ep[mask_plot_ep], policy_cf_ep[mask_plot_ep],
    color = :blue, lw = 2, linestyle = :dash,
    label = ""
)

p_inf_ep = plot(
    dates_plot_ep[mask_plot_ep], infl_base_ep[mask_plot_ep],
    color = :black, lw = 2, linestyle = :solid,
    title = "Core PCE inflation (%, yoy)",
    label = "",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_inf_ep,
    dates_plot_ep[mask_plot_ep], infl_remove_cf_ep[mask_plot_ep],
    color = :red, lw = 2, linestyle = :dashdot,
    label = ""
)
plot!(
    p_inf_ep,
    dates_plot_ep[mask_plot_ep], infl_cf_ep[mask_plot_ep],
    color = :blue, lw = 2, linestyle = :dash,
    label = ""
)

p_out_ep = plot(
    dates_plot_ep[mask_plot_ep], output_base_ep[mask_plot_ep].+2,
    color = :black, lw = 2, linestyle = :solid,
    title = "Output deviation from trend (%)",
    label = "",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_out_ep,
    dates_plot_ep[mask_plot_ep], output_remove_cf_ep[mask_plot_ep].+2,
    color = :red, lw = 2, linestyle = :dashdot,
    label = ""
)
plot!(
    p_out_ep,
    dates_plot_ep[mask_plot_ep], output_cf_ep[mask_plot_ep].+2,
    color = :blue, lw = 2, linestyle = :dash,
    label = ""
)

p_rstar_ep = plot(
    dates_plot_ep[mask_plot_ep], rstar_base_ep[mask_plot_ep],
    color = :black, lw = 2, linestyle = :solid,
    title = "r* (APR)",
    label = "",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_rstar_ep,
    dates_plot_ep[mask_plot_ep], rstar_remove_cf_ep[mask_plot_ep],
    color = :red, lw = 2, linestyle = :dashdot,
    label = ""
)
plot!(
    p_rstar_ep,
    dates_plot_ep[mask_plot_ep], rstar_cf_ep[mask_plot_ep],
    color = :blue, lw = 2, linestyle = :dash,
    label = ""
)

p_3x2_ep = plot(
    p_cy_ep, p_long_ep,
    p_pol_ep, p_inf_ep,
    p_out_ep, p_rstar_ep,
    layout = (3, 2),
    size = (1100, 1200)
)

display(p_3x2_ep)

if use_FG_in_EA
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege 3x2 baseline_plus_delta.pdf"))
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege 3x2 baseline_plus_delta.png"))
else
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege without FG shocks in EA 3x2 baseline_plus_delta.pdf"))
    savefig(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant privilege without FG shocks in EA 3x2 baseline_plus_delta.png"))
end

##########################################################################################
## APPENDED: US ONLY 3x2 PANEL, BLACK AND RED LINES
##########################################################################################

function _us_black_red_panel(dates_plot, baseline, no_cy, title, xticks; show_legend = false)
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
        no_cy,
        color = :red,
        lw = 2,
        linestyle = :dashdot,
        label = show_legend ? "US - No CY shocks" : "",
    )
    show_legend && plot!(p, legend = first_panel_legend_pos)
    return p
end

p_cy_us_only = _us_black_red_panel(
    dates_plot_ep[mask_plot_ep],
    cy_base_ep[mask_plot_ep],
    cy_remove_cf_ep[mask_plot_ep],
    "Convenience yield (APR)",
    (year_tick_dates_ep, year_tick_labels_ep);
    show_legend = true,
)
p_long_us_only = _us_black_red_panel(
    dates_plot_ep[mask_plot_ep],
    long_base_ep[mask_plot_ep],
    long_remove_cf_ep[mask_plot_ep],
    "Long-term interest rate (APR)",
    (year_tick_dates_ep, year_tick_labels_ep),
)
p_pol_us_only = _us_black_red_panel(
    dates_plot_ep[mask_plot_ep],
    policy_base_ep[mask_plot_ep],
    policy_remove_cf_ep[mask_plot_ep],
    "Policy rate (APR)",
    (year_tick_dates_ep, year_tick_labels_ep),
)
p_inf_us_only = _us_black_red_panel(
    dates_plot_ep[mask_plot_ep],
    infl_base_ep[mask_plot_ep],
    infl_remove_cf_ep[mask_plot_ep],
    "Core PCE inflation (%, yoy)",
    (year_tick_dates_ep, year_tick_labels_ep),
)
p_out_us_only = _us_black_red_panel(
    dates_plot_ep[mask_plot_ep],
    output_base_ep[mask_plot_ep] .+ 2,
    output_remove_cf_ep[mask_plot_ep] .+ 2,
    "Output deviation from trend (%)",
    (year_tick_dates_ep, year_tick_labels_ep),
)
p_rstar_us_only = _us_black_red_panel(
    dates_plot_ep[mask_plot_ep],
    rstar_base_ep[mask_plot_ep],
    rstar_remove_cf_ep[mask_plot_ep],
    "r* (APR)",
    (year_tick_dates_ep, year_tick_labels_ep),
)

p_3x2_us_only = plot(
    p_cy_us_only, p_long_us_only,
    p_pol_us_only, p_inf_us_only,
    p_out_us_only, p_rstar_us_only,
    layout = (3, 2),
    size = (1100, 1200)
)

display(p_3x2_us_only)

if use_FG_in_EA
    savefig(p_3x2_us_only, joinpath(saveroot, "Final Paper", "Figures", "US no exorbitant privilege 3x2 baseline_plus_delta.pdf"))
    savefig(p_3x2_us_only, joinpath(saveroot, "Final Paper", "Figures", "US no exorbitant privilege 3x2 baseline_plus_delta.png"))
else
    savefig(p_3x2_us_only, joinpath(saveroot, "Final Paper", "Figures", "US no exorbitant privilege without FG shocks in EA 3x2 baseline_plus_delta.pdf"))
    savefig(p_3x2_us_only, joinpath(saveroot, "Final Paper", "Figures", "US no exorbitant privilege without FG shocks in EA 3x2 baseline_plus_delta.png"))
end

# -------------------------
# Export plotted series
# -------------------------
df_plot_3x2_ep = DataFrame(
    Date = dates_plot_ep[mask_plot_ep],

    ConvenienceYield_Baseline_APR = cy_base_ep[mask_plot_ep],
    ConvenienceYield_USExcludingMatchedPrivilege_APR = cy_remove_cf_ep[mask_plot_ep],
    ConvenienceYield_Counterfactual_APR = cy_cf_ep[mask_plot_ep],

    LongRate_Baseline_APR = long_base_ep[mask_plot_ep],
    LongRate_USExcludingMatchedPrivilege_APR = long_remove_cf_ep[mask_plot_ep],
    LongRate_Counterfactual_APR = long_cf_ep[mask_plot_ep],

    PolicyRate_Baseline_APR = policy_base_ep[mask_plot_ep],
    PolicyRate_USExcludingMatchedPrivilege_APR = policy_remove_cf_ep[mask_plot_ep],
    PolicyRate_Counterfactual_APR = policy_cf_ep[mask_plot_ep],

    CorePCEInflation_Baseline_YY = infl_base_ep[mask_plot_ep],
    CorePCEInflation_USExcludingMatchedPrivilege_YY = infl_remove_cf_ep[mask_plot_ep],
    CorePCEInflation_Counterfactual_YY = infl_cf_ep[mask_plot_ep],

    OutputGap_Baseline = output_base_ep[mask_plot_ep],
    OutputGap_USExcludingMatchedPrivilege = output_remove_cf_ep[mask_plot_ep],
    OutputGap_Counterfactual = output_cf_ep[mask_plot_ep],

    Forward5YearRealNaturalRate_Baseline_APR = rstar_base_ep[mask_plot_ep],
    Forward5YearRealNaturalRate_USExcludingMatchedPrivilege_APR = rstar_remove_cf_ep[mask_plot_ep],
    Forward5YearRealNaturalRate_Counterfactual_APR = rstar_cf_ep[mask_plot_ep]
)

if use_FG_in_EA
    CSV.write(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant_privilege_3x2_baseline_plus_delta.csv"), df_plot_3x2_ep)
else
    CSV.write(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant_privilege_without_FG_shocks_in_EA_3x2_baseline_plus_delta.csv"), df_plot_3x2_ep)
end

##########################################################################################
## APPENDED: EXPORT FULL-SAMPLE CONVENIENCE YIELD SERIES FOR US AND EA
##########################################################################################

# EA full-sample convenience yield
df_EA_full_ep = df[end-size(states[smoother], 2)+1:end, :]
states_EA_full_ep = states[smoother]
pseudo_EA_full_ep = pseudo[smoother]
dates_EA_full_ep = df_EA_full_ep.date

cy_EA_full_ep = _get_convenience_yield_series(m, states_EA_full_ep, nothing, pseudo_EA_full_ep)
ea_b_liqtil_full_ep  = _get_series(m, states_EA_full_ep, nothing, nothing, :b_liqtil_t)
ea_b_liqp_full_ep    = _get_series(m, states_EA_full_ep, nothing, nothing, :b_liqp_t)
ea_b_safetil_full_ep = _get_series(m, states_EA_full_ep, nothing, nothing, :b_safetil_t)
ea_b_safep_full_ep   = _get_series(m, states_EA_full_ep, nothing, nothing, :b_safep_t)
ea_liquidity_state_total_full_ep = ea_b_liqtil_full_ep .+ ea_b_liqp_full_ep
ea_safety_state_total_full_ep = ea_b_safetil_full_ep .+ ea_b_safep_full_ep
ea_cy_state_total_full_ep = ea_liquidity_state_total_full_ep .+ ea_safety_state_total_full_ep
ea_liquidity_trend_full_ep = fill(100 * log(m[:lnb_liq]), length(dates_EA_full_ep))
ea_safety_trend_full_ep = fill(100 * log(m[:lnb_safe]), length(dates_EA_full_ep))
ea_cy_trend_full_ep = ea_liquidity_trend_full_ep .+ ea_safety_trend_full_ep
ea_liquidity_cy_full_ep = ea_liquidity_state_total_full_ep .+ ea_liquidity_trend_full_ep
ea_safety_cy_full_ep = ea_safety_state_total_full_ep .+ ea_safety_trend_full_ep

# US full-sample convenience yield
cy_US_full_ep = _get_convenience_yield_series(m1, states_US_sm_ep, nothing, pseudo_US_sm_ep)
us_b_liqtil_full_ep  = _get_series(m1, states_US_sm_ep, nothing, nothing, :b_liqtil_t)
us_b_liqp_full_ep    = _get_series(m1, states_US_sm_ep, nothing, nothing, :b_liqp_t)
us_b_safetil_full_ep = _get_series(m1, states_US_sm_ep, nothing, nothing, :b_safetil_t)
us_b_safep_full_ep   = _get_series(m1, states_US_sm_ep, nothing, nothing, :b_safep_t)
us_liquidity_state_total_full_ep = us_b_liqtil_full_ep .+ us_b_liqp_full_ep
us_safety_state_total_full_ep = us_b_safetil_full_ep .+ us_b_safep_full_ep
us_cy_state_total_full_ep = us_liquidity_state_total_full_ep .+ us_safety_state_total_full_ep
us_liquidity_trend_full_ep = fill(100 * log(m1[:lnb_liq]), length(dates_US_ep))
us_safety_trend_full_ep = fill(100 * log(m1[:lnb_safe]), length(dates_US_ep))
us_cy_trend_full_ep = us_liquidity_trend_full_ep .+ us_safety_trend_full_ep
us_liquidity_cy_full_ep = us_liquidity_state_total_full_ep .+ us_liquidity_trend_full_ep
us_safety_cy_full_ep = us_safety_state_total_full_ep .+ us_safety_trend_full_ep

# Optional APR versions
cy_EA_full_apr_ep = 4 .* cy_EA_full_ep
cy_US_full_apr_ep = 4 .* cy_US_full_ep
ea_cy_state_total_full_apr_ep = 4 .* ea_cy_state_total_full_ep
us_cy_state_total_full_apr_ep = 4 .* us_cy_state_total_full_ep
ea_cy_trend_full_apr_ep = 4 .* ea_cy_trend_full_ep
us_cy_trend_full_apr_ep = 4 .* us_cy_trend_full_ep
ea_liquidity_cy_full_apr_ep = 4 .* ea_liquidity_cy_full_ep
ea_safety_cy_full_apr_ep = 4 .* ea_safety_cy_full_ep
us_liquidity_cy_full_apr_ep = 4 .* us_liquidity_cy_full_ep
us_safety_cy_full_apr_ep = 4 .* us_safety_cy_full_ep
ea_liquidity_state_total_full_apr_ep = 4 .* ea_liquidity_state_total_full_ep
us_liquidity_state_total_full_apr_ep = 4 .* us_liquidity_state_total_full_ep
ea_liquidity_trend_full_apr_ep = 4 .* ea_liquidity_trend_full_ep
us_liquidity_trend_full_apr_ep = 4 .* us_liquidity_trend_full_ep
ea_b_liqtil_full_apr_ep = 4 .* ea_b_liqtil_full_ep
us_b_liqtil_full_apr_ep = 4 .* us_b_liqtil_full_ep
ea_b_liqp_full_apr_ep = 4 .* ea_b_liqp_full_ep
us_b_liqp_full_apr_ep = 4 .* us_b_liqp_full_ep
ea_safety_state_total_full_apr_ep = 4 .* ea_safety_state_total_full_ep
us_safety_state_total_full_apr_ep = 4 .* us_safety_state_total_full_ep
ea_safety_trend_full_apr_ep = 4 .* ea_safety_trend_full_ep
us_safety_trend_full_apr_ep = 4 .* us_safety_trend_full_ep
ea_b_safetil_full_apr_ep = 4 .* ea_b_safetil_full_ep
us_b_safetil_full_apr_ep = 4 .* us_b_safetil_full_ep
ea_b_safep_full_apr_ep = 4 .* ea_b_safep_full_ep
us_b_safep_full_apr_ep = 4 .* us_b_safep_full_ep

# Export as one CSV with full samples preserved, APR-only. Within each model
# block, the grouped columns add up cleanly and the overall total is last.
df_cy_EA_ep = DataFrame(Date = dates_EA_full_ep,
                        LiquidityTransitory_EA_APR = ea_b_liqtil_full_apr_ep,
                        LiquidityPermanent_EA_APR = ea_b_liqp_full_apr_ep,
                        LiquidityStateTotal_EA_APR = ea_liquidity_state_total_full_apr_ep,
                        LiquidityTrend_EA_APR = ea_liquidity_trend_full_apr_ep,
                        LiquidityConvenienceYield_EA_APR = ea_liquidity_cy_full_apr_ep,
                        SafetyTransitory_EA_APR = ea_b_safetil_full_apr_ep,
                        SafetyPermanent_EA_APR = ea_b_safep_full_apr_ep,
                        SafetyStateTotal_EA_APR = ea_safety_state_total_full_apr_ep,
                        SafetyTrend_EA_APR = ea_safety_trend_full_apr_ep,
                        SafetyConvenienceYield_EA_APR = ea_safety_cy_full_apr_ep,
                        ConvenienceYieldStateTotal_EA_APR = ea_cy_state_total_full_apr_ep,
                        ConvenienceYieldTrend_EA_APR = ea_cy_trend_full_apr_ep,
                        ConvenienceYield_EA_APR = cy_EA_full_apr_ep)

df_cy_US_ep = DataFrame(Date = dates_US_ep,
                        LiquidityTransitory_US_APR = us_b_liqtil_full_apr_ep,
                        LiquidityPermanent_US_APR = us_b_liqp_full_apr_ep,
                        LiquidityStateTotal_US_APR = us_liquidity_state_total_full_apr_ep,
                        LiquidityTrend_US_APR = us_liquidity_trend_full_apr_ep,
                        LiquidityConvenienceYield_US_APR = us_liquidity_cy_full_apr_ep,
                        SafetyTransitory_US_APR = us_b_safetil_full_apr_ep,
                        SafetyPermanent_US_APR = us_b_safep_full_apr_ep,
                        SafetyStateTotal_US_APR = us_safety_state_total_full_apr_ep,
                        SafetyTrend_US_APR = us_safety_trend_full_apr_ep,
                        SafetyConvenienceYield_US_APR = us_safety_cy_full_apr_ep,
                        ConvenienceYieldStateTotal_US_APR = us_cy_state_total_full_apr_ep,
                        ConvenienceYieldTrend_US_APR = us_cy_trend_full_apr_ep,
                        ConvenienceYield_US_APR = cy_US_full_apr_ep)

df_cy_full_ep = join(df_cy_EA_ep, df_cy_US_ep, on = :Date,  kind = :outer)

CSV.write(joinpath(saveroot, "Final Paper", "Figures", "ConvenienceYield_fullsample_US_EA.csv"), df_cy_full_ep)


