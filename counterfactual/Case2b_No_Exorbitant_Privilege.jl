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

first_panel_legend_pos = (0.25, 0.17)



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

function _get_convenience_yield_series(m, states_mat, obs_mat, pseudo_mat)
    # First try a total convenience-yield variable
    total_candidates = [
        :ConvenienceYield, :convenience_yield, :convenienceyield,
        :cy, :cy_t, :privilege, :privilege_t
    ]
    total_sym = _find_first_symbol(m, total_candidates)
    if total_sym !== nothing
        return _get_series(m, states_mat, obs_mat, pseudo_mat, total_sym)
    end

    # Otherwise try summing liquid + safety components
    liq_candidates  = [:lnb_liq, :lnb_liq_t, :b_liq, :b_liq_t, :liq_convenience_yield]
    safe_candidates = [:lnb_safe, :lnb_safe_t, :b_safe, :b_safe_t, :safe_convenience_yield]

    liq_sym  = _find_first_symbol(m, liq_candidates)
    safe_sym = _find_first_symbol(m, safe_candidates)

    if liq_sym !== nothing && safe_sym !== nothing
        return _get_series(m, states_mat, obs_mat, pseudo_mat, liq_sym) .+
               _get_series(m, states_mat, obs_mat, pseudo_mat, safe_sym)
    end

    error("""
Could not locate a convenience yield series automatically.
Please update candidate names in _get_convenience_yield_series(...) to match your model.
""")
end

function _get_long_rate_series(m, states_mat, obs_mat, pseudo_mat)
    candidates = [
        :obs_longrate, :obs_long_rate, :longrate, :LongRate,
        :LongTermRate, :LongTermNominalRate, :tenyearrate, :obs_tenyearrate
    ]
    sym = _find_first_symbol(m, candidates)
    sym === nothing && error("""
Could not locate the long-term interest rate series automatically.
Please update candidate names in _get_long_rate_series(...) to match your model.
""")
    return _get_series(m, states_mat, obs_mat, pseudo_mat, sym)
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
dates_US = df_US.date[end-size(states_US_sm, 2)+1:end]

# -------------------------
# Build the additional plotted series
# -------------------------
cy_US = _get_convenience_yield_series(m1, states_US_sm, nothing, pseudo_US_sm)
cy_CF = _get_convenience_yield_series(m1, states_privilege, obs_privilege, pseudo_privilege)

longrate_CF = _get_long_rate_series(m1, states_privilege, obs_privilege, pseudo_privilege)

# Existing counterfactual series from your current script
policy_CF    = vec(obs_privilege[m1.observables[:obs_nominalrate], :])
inflation_CF = vec(obs_privilege[m1.observables[:obs_corepce], :])
output_CF    = vec(states_privilege[m1.endogenous_states[:y_t], :])
rstar_CF     = vec(pseudo_privilege[m1.pseudo_observables[:Forward5YearRealNaturalRate], :])

# Align sample start with your existing plotstart
plot_start_date = dates[plotstart]
mask_CF = dates .>= plot_start_date
mask_US = dates_US .>= plot_start_date

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
    dates_US[mask_US], cy_US[mask_US],
    color = :black, lw = 2,
    title = "Convenience yield",
    label = "US",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_cy, dates[mask_US], zeros(sum(mask_US)), lc = :black, lw = 2, label = "",xticks = (year_tick_dates, year_tick_labels))

plot!(
    p_cy,
    dates[mask_CF], cy_CF[mask_CF],
    color = :blue, lw = 2,
    label = "Counterfactual",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_cy, legend = first_panel_legend_pos)

p_long = plot(
    dates[mask_CF], longrate_CF[mask_CF],
    color = :blue, lw = 2,
    title = "Long-term interest rate",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_long, dates[mask_CF], zeros(sum(mask_CF)), lc = :black, lw = 2, label = "",xticks = (year_tick_dates, year_tick_labels))

p_pol = plot(
    dates[mask_CF], policy_CF[mask_CF],
    color = :blue, lw = 2,
    title = "Policy rate",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_pol, dates[mask_CF], zeros(sum(mask_CF)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_inf = plot(
    dates[mask_CF], inflation_CF[mask_CF],
    color = :blue, lw = 2,
    title = "Inflation",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_inf, dates[mask_CF], zeros(sum(mask_CF)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_out = plot(
    dates[mask_CF], output_CF[mask_CF],
    color = :blue, lw = 2,
    title = "Output",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_out, dates[mask_CF], zeros(sum(mask_CF)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_rstar = plot(
    dates[mask_CF], rstar_CF[mask_CF],
    color = :blue, lw = 2,
    title = "r* (Forward 5-year real natural rate)",
    label = "",
    xticks = (year_tick_dates, year_tick_labels)
)
plot!(p_rstar, dates[mask_CF], zeros(sum(mask_CF)), lc = :black, lw = 2, label = "", xticks = (year_tick_dates, year_tick_labels))

p_3x2 = plot(
    p_cy, p_long,
    p_pol, p_inf,
    p_out, p_rstar,
    layout = (3, 2),
    legend = false,
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

# # -------------------------
# # Export the plotted series
# # -------------------------
# df_plot_cf = DataFrame(
#     Date = dates[mask_CF],
#     ConvenienceYieldCounterfactual = cy_CF[mask_CF],
#     LongTermRateCounterfactual = longrate_CF[mask_CF],
#     PolicyRateCounterfactual = policy_CF[mask_CF],
#     InflationCounterfactual = inflation_CF[mask_CF],
#     OutputCounterfactual = output_CF[mask_CF],
#     Forward5YearRealNaturalRateCounterfactual = rstar_CF[mask_CF]
# )

# df_plot_us = DataFrame(
#     Date = dates_US[mask_US],
#     ConvenienceYieldUS = cy_US[mask_US]
# )

# df_plot_3x2 = leftjoin(df_plot_cf, df_plot_us, on = :Date)

# if use_FG_in_EA
#     CSV.write(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant_privilege_3x2.csv"), df_plot_3x2)
# else
#     CSV.write(joinpath(saveroot, "Final Paper", "Figures", "Exorbitant_privilege_without_FG_shocks_in_EA_3x2.csv"), df_plot_3x2)
# end
##########################################################################################
## APPENDED: REWORKED 3x2 PANEL WITH LEVEL-MATCHED COUNTERFACTUALS
## - Common-date alignment between EA and US samples
## - Black solid = US baseline
## - Blue dashed = US baseline + normalized (EA-US) convenience-yield-shock delta
## - Start date = 1998Q4
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

function _ep_resolve_convenience_yield(m)
    total_candidates = [
        :ConvenienceYield, :convenience_yield, :convenienceyield,
        :cy, :cy_t, :privilege, :privilege_t
    ]
    total_sym = _ep_find_first_symbol(m, total_candidates)
    if total_sym !== nothing
        return (:single, total_sym, nothing)
    end

    liq_candidates  = [:lnb_liq, :lnb_liq_t, :b_liq, :b_liq_t, :liq_convenience_yield]
    safe_candidates = [:lnb_safe, :lnb_safe_t, :b_safe, :b_safe_t, :safe_convenience_yield]

    liq_sym  = _ep_find_first_symbol(m, liq_candidates)
    safe_sym = _ep_find_first_symbol(m, safe_candidates)

    if liq_sym !== nothing && safe_sym !== nothing
        return (:sum, liq_sym, safe_sym)
    end

    error("""
Could not locate a convenience yield series automatically.
Please update candidate names in _ep_resolve_convenience_yield(...) to match your model.
""")
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

_ep_normalize_delta(delta::AbstractVector, anchor_idx::Int) = delta .- delta[anchor_idx]

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
# Build EA-US convenience-yield shock DIFFERENCE on common sample
# -------------------------
privilege_shock_vals_EA_ep = convert(Matrix, shocks_df[smoother][idx_EA_ep, privilege_shock_names])

privilege_shock_vals_US_ep = hcat(
    [vec(shocks_US_sm_ep[m1.exogenous_shocks[sh], idx_US_ep]) for sh in privilege_shock_names]...
)

delta_privilege_shock_vals_ep = privilege_shock_vals_EA_ep .- privilege_shock_vals_US_ep

# -------------------------
# Feed the delta shocks through the US model
# -------------------------
horizon_cf_ep = size(delta_privilege_shock_vals_ep, 1)
nshocks_US_ep = size(system_US[:RRR], 2)
nstates_US_ep = size(system_US[:TTT], 1)

s_0_US_ep = zeros(nstates_US_ep)
shocks_delta_ep = zeros(nshocks_US_ep, horizon_cf_ep)

for (i, shock_name) in enumerate(privilege_shock_names)
    shock_ind_US = m1.exogenous_shocks[shock_name]
    shocks_delta_ep[shock_ind_US, :] .= delta_privilege_shock_vals_ep[:, i]
end

system_US_delta_ep = DSGE.zero_system_constants(system_US)
states_delta_ep, obs_delta_ep, pseudo_delta_ep = forecast(system_US_delta_ep, s_0_US_ep, shocks_delta_ep)

# -------------------------
# Resolve plotted variables
# -------------------------
cy_kind_ep, cy_sym1_ep, cy_sym2_ep = _ep_resolve_convenience_yield(m1)
longrate_sym_ep = _ep_resolve_long_rate(m1)

# -------------------------
# Align US baseline objects to common dates
# -------------------------
df_US_common_ep = df_US_aligned_ep[idx_US_ep, :]
states_US_common_ep = states_US_sm_ep[:, idx_US_ep]
pseudo_US_common_ep = pseudo_US_sm_ep[:, idx_US_ep]

# -------------------------
# Build RAW baseline and delta series first
# -------------------------

# Convenience yield
if cy_kind_ep == :single
    cy_base_raw_ep  = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, cy_sym1_ep)
    cy_delta_raw_ep = _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, cy_sym1_ep)
else
    cy_base_raw_ep  = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, cy_sym1_ep) .+
                      _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, cy_sym2_ep)
    cy_delta_raw_ep = _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, cy_sym1_ep) .+
                      _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, cy_sym2_ep)
end

# Long rate
long_base_raw_ep  = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, longrate_sym_ep)
long_delta_raw_ep = _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, longrate_sym_ep)

# Policy rate
policy_base_raw_ep  = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :obs_nominalrate)
policy_delta_raw_ep = _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, :obs_nominalrate)

# Inflation (now obs_corepce)
infl_base_raw_ep  = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :obs_corepce)
infl_delta_raw_ep = _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, :obs_corepce)

# Output gap
output_base_raw_ep  = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :y_t)
output_delta_raw_ep = _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, :y_t)

# r*
rstar_base_raw_ep  = _ep_get_baseline_series(m1, df_US_common_ep, states_US_common_ep, pseudo_US_common_ep, :Forward5YearRealNaturalRate)
rstar_delta_raw_ep = _ep_get_delta_series(m1, states_delta_ep, obs_delta_ep, pseudo_delta_ep, :Forward5YearRealNaturalRate)

# -------------------------
# Start date and anchor normalization date: 1998Q4
# -------------------------
plot_start_date_ep = quartertodate("1998-Q4")
anchor_idx_ep = findfirst(==(plot_start_date_ep), dates_plot_ep)
anchor_idx_ep === nothing && error("Anchor date 1998Q4 not found in dates_plot_ep.")
anchor_idx_ep < 4 && error("Need at least 3 quarters before 1998Q4 to construct y/y inflation.")

mask_plot_ep = dates_plot_ep .>= plot_start_date_ep

# -------------------------
# Convert to plotted units
# -------------------------
# APR variables: *4
cy_base_ep     = 4 .* cy_base_raw_ep
cy_delta_ep    = 4 .* cy_delta_raw_ep

long_base_ep   = 4 .* long_base_raw_ep
long_delta_ep  = 4 .* long_delta_raw_ep

policy_base_ep = 4 .* policy_base_raw_ep
policy_delta_ep = 4 .* policy_delta_raw_ep

rstar_base_ep  = 4 .* rstar_base_raw_ep
rstar_delta_ep = 4 .* rstar_delta_raw_ep

# Inflation: year-over-year = 4-quarter sum of quarterly inflation
infl_base_ep   = _ep_four_quarter_sum(infl_base_raw_ep)
infl_delta_ep  = _ep_four_quarter_sum(infl_delta_raw_ep)

# Output gap: leave as-is
output_base_ep  = output_base_raw_ep
output_delta_ep = output_delta_raw_ep

# -------------------------
# Normalize deltas to start from the same level at 1998Q4
# -------------------------
cy_cf_ep     = cy_base_ep     .+ _ep_normalize_delta(cy_delta_ep, anchor_idx_ep)
long_cf_ep   = long_base_ep   .+ _ep_normalize_delta(long_delta_ep, anchor_idx_ep)
policy_cf_ep = policy_base_ep .+ _ep_normalize_delta(policy_delta_ep, anchor_idx_ep)
infl_cf_ep   = infl_base_ep   .+ _ep_normalize_delta(infl_delta_ep, anchor_idx_ep)
output_cf_ep = output_base_ep .+ _ep_normalize_delta(output_delta_ep, anchor_idx_ep)
rstar_cf_ep  = rstar_base_ep  .+ _ep_normalize_delta(rstar_delta_ep, anchor_idx_ep)

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
    label = "US baseline",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_cy_ep,
    dates_plot_ep[mask_plot_ep], cy_cf_ep[mask_plot_ep],
    color = :blue, lw = 2, linestyle = :dash,
    label = "Counterfactual"
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
    dates_plot_ep[mask_plot_ep], infl_cf_ep[mask_plot_ep],
    color = :blue, lw = 2, linestyle = :dash,
    label = ""
)

p_out_ep = plot(
    dates_plot_ep[mask_plot_ep], output_base_ep[mask_plot_ep],
    color = :black, lw = 2, linestyle = :solid,
    title = "Output deviation from trend (%)",
    label = "",
    xticks = (year_tick_dates_ep, year_tick_labels_ep)
)
plot!(
    p_out_ep,
    dates_plot_ep[mask_plot_ep], output_cf_ep[mask_plot_ep],
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

# -------------------------
# Export plotted series
# -------------------------
df_plot_3x2_ep = DataFrame(
    Date = dates_plot_ep[mask_plot_ep],

    ConvenienceYield_Baseline_APR = cy_base_ep[mask_plot_ep],
    ConvenienceYield_Counterfactual_APR = cy_cf_ep[mask_plot_ep],

    LongRate_Baseline_APR = long_base_ep[mask_plot_ep],
    LongRate_Counterfactual_APR = long_cf_ep[mask_plot_ep],

    PolicyRate_Baseline_APR = policy_base_ep[mask_plot_ep],
    PolicyRate_Counterfactual_APR = policy_cf_ep[mask_plot_ep],

    CorePCEInflation_Baseline_YY = infl_base_ep[mask_plot_ep],
    CorePCEInflation_Counterfactual_YY = infl_cf_ep[mask_plot_ep],

    OutputGap_Baseline = output_base_ep[mask_plot_ep],
    OutputGap_Counterfactual = output_cf_ep[mask_plot_ep],

    Forward5YearRealNaturalRate_Baseline_APR = rstar_base_ep[mask_plot_ep],
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

# Resolve convenience-yield symbols in each model
cy_kind_EA_ep, cy_sym1_EA_ep, cy_sym2_EA_ep = _ep_resolve_convenience_yield(m)
cy_kind_US_ep, cy_sym1_US_ep, cy_sym2_US_ep = _ep_resolve_convenience_yield(m1)

# EA full-sample convenience yield
df_EA_full_ep = df[end-size(states[smoother], 2)+1:end, :]
states_EA_full_ep = states[smoother]
pseudo_EA_full_ep = pseudo[smoother]
dates_EA_full_ep = df_EA_full_ep.date

if cy_kind_EA_ep == :single
    cy_EA_full_ep = _ep_get_baseline_series(m, df_EA_full_ep, states_EA_full_ep, pseudo_EA_full_ep, cy_sym1_EA_ep)
else
    cy_EA_full_ep = _ep_get_baseline_series(m, df_EA_full_ep, states_EA_full_ep, pseudo_EA_full_ep, cy_sym1_EA_ep) .+
                    _ep_get_baseline_series(m, df_EA_full_ep, states_EA_full_ep, pseudo_EA_full_ep, cy_sym2_EA_ep)
end

# US full-sample convenience yield
if cy_kind_US_ep == :single
    cy_US_full_ep = _ep_get_baseline_series(m1, df_US_aligned_ep, states_US_sm_ep, pseudo_US_sm_ep, cy_sym1_US_ep)
else
    cy_US_full_ep = _ep_get_baseline_series(m1, df_US_aligned_ep, states_US_sm_ep, pseudo_US_sm_ep, cy_sym1_US_ep) .+
                    _ep_get_baseline_series(m1, df_US_aligned_ep, states_US_sm_ep, pseudo_US_sm_ep, cy_sym2_US_ep)
end

# Optional APR versions
cy_EA_full_apr_ep = 4 .* cy_EA_full_ep
cy_US_full_apr_ep = 4 .* cy_US_full_ep

# Export as one CSV with full samples preserved
df_cy_EA_ep = DataFrame(Date = dates_EA_full_ep,
                        ConvenienceYield_EA = cy_EA_full_ep,
                        ConvenienceYield_EA_APR = cy_EA_full_apr_ep)

df_cy_US_ep = DataFrame(Date = dates_US_ep,
                        ConvenienceYield_US = cy_US_full_ep,
                        ConvenienceYield_US_APR = cy_US_full_apr_ep)

df_cy_full_ep = join(df_cy_EA_ep, df_cy_US_ep, on = :Date,  kind = :outer)

CSV.write(joinpath(saveroot, "Final Paper", "Figures", "ConvenienceYield_fullsample_US_EA.csv"), df_cy_full_ep)