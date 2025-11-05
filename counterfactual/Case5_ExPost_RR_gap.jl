using DSGE, Dates, DataFrames,OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics, ModelConstructors, LinearAlgebra
# using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics
using Measures  # for mm margins

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl and the test forecast driver file
# 25 August 2025
#############################

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")
# params_mode = load_draws(m, :mode)
# # Switch off the pi_target shock
# params_mode[m.param_index[:σ_π_target_sh]] = 0.0
# DSGE.update!(m, params_mode)

# Settings for data, paths, etc.
m <= DSGE.Setting(:data_vintage, "161223")
# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))
#  mode_file = rawpath(m, "estimate", "paramsmode.h5")
#         #mode_file = replace(mode_file, "ss20", "ss18")
#         DSGE.update!(m, h5read(mode_file, "params"))


system = DSGE.compute_system(m)



nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)


""""
obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, state_ind::Int, shock_inds::Vector{Int},
                                                    system::System{Float64})

Given a desired path `x` for state `state_ind` over `horizon` periods, and a vector of shock indices
`shock_inds` (one per period), back out the necessary values of shocks over those periods such that
s^i_{1:h} = x_{1:h}, using the model's forecast function.

Returns a matrix of required shocks of size (nshocks, horizon).
"""
function obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, m::AbstractDSGEModel,var_name::Symbol, shock_inds::Vector{Int},
                                                         system::System{Float64})


    horizon = length(x)
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    s_0 = zeros(nstates)
    shocks = zeros(nshocks, horizon)
    var_names, var_class, peg_ind =
    if var_name in keys(m.endogenous_states)
        m.endogenous_states, :states, m.endogenous_states[var_name]

    else
        m.observables, :obs,  m.observables[var_name]
    end


    for t in 1:horizon
        
        # Compute IRF for a unit shock at time t for the specified shock
        test_shocks = zeros(nshocks, horizon)
        test_shocks[shock_inds[t], t] = 1.0
        states, obs, _ = forecast(system, s_0, test_shocks)
                if var_class == :states
                    irf = states[peg_ind, t] # Impact of a unit shock at t on state at t
                else 
                    irf = obs[peg_ind,t] # Impact of a unit shock at t on obs at t
                end
        # Compute effect of previous shocks
        prev_effect = 0.0
        if t > 1
            prev_shocks = shocks[:, 1:t-1]
            prev_states, prev_obs, _ = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
            if var_class == :states
                prev_effect = prev_states[peg_ind, t]
            else 
                prev_effect = prev_obs[peg_ind,t] # Impact of a unit shock at t on obs at t
            end
        end

        # Required shock at time t to achieve desired value
        shocks[shock_inds[t], t] = (x[t] - prev_effect) / irf
    end

    return shocks
end

# Example test (assuming you have a model m and system):
# x = [1.0, 2.0]
# shock_inds = [m.exogenous_shocks[:shock1] m.exogenous_shocks[:shock1];
#               m.exogenous_shocks[:shock2] m.exogenous_shocks[:shock2]]
# shocks = obtain_shocks_from_desired_state_path_iterative(x, m, :obs_nominalrate, shock_inds, system)


# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]

## Load in HLW real time estiamtes of R*
csv_path = joinpath(basepath, "Main results","US", "Ex_post_real_rate_gaps.csv")
US_dataset = DataFrame(CSV.File(csv_path))

valid_idx = findall(row -> !ismissing(row[:date]) && !ismissing(row[:Target]), eachrow(US_dataset))
dates= US_dataset.date[valid_idx]

##### Alternative if realrate gap change of HLW was implemented using policy rate shocks alone
desired_path = skipmissing(US_dataset.Target[valid_idx]) |> collect  # Desired path for the state variable
desired_path = -desired_path
var_name =:obs_nominalrate
horizon = size(desired_path, 1)
# desired_path = [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
shock_name =   :rm_sh # MP shock implements the real rate gap change
# shock_syms = collect(Iterators.filter(k -> startswith(String(k), "rm_sh"), keys(m.exogenous_shocks)))

# Ensure `shock_inds` is a 1D array
shock_ind = 10


var_names, var_class =
    if var_name in keys(m.endogenous_states)
        m.endogenous_states, :states
    else
        m.observables, :obs
    end
    exo          = m.exogenous_shocks
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 1)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    states = zeros(nstates, horizon, nshocks)
    obs    = zeros(nobs,    horizon, nshocks)
    pseudo = zeros(npseudo, horizon, nshocks)

    # Set constant system matrices to 0
    system = DSGE.zero_system_constants(system)

    s_0 = zeros(nstates)

    # Isolate single shock
    shocks = zeros(nshocks, horizon)
    for t = 1:horizon
        var_value = desired_path[t]
        if var_class == :states
            var_value_att = var_value - obs[m.endogenous_states[var_name],t, m.exogenous_shocks[shock_name]]
            shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_state_value(var_value_att,
                                                                        var_names[var_name],
                                                                        exo[shock_name],
                                                                        system[:RRR])
        else # == :obs
            var_value_att = var_value - obs[m.observables[var_name],t, m.exogenous_shocks[shock_name]]
            shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(var_value_att,
                                                                        var_names[var_name],
                                                                        exo[shock_name],
                                                                        system[:ZZ],
                                                                        system[:RRR])
        end
    
    # Iterate state space forward
    states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)
    end

    # states, obs, pseudo = forecast(system, s_0, shocks)

# --- Step 1: Compute IRFs for each shock ---
plotvars = [:Forward5YearRealNaturalRate,:obs_nominalrate,  :pi_t,  :y_t,:ExAnteRealRate,:RealNaturalRate] # Output, Inflation, Policy Rate, R*
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

using Plots
p1 = plot(plotdates,desired_path,title="Targeted Path")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:,exo[shock_name]].*4,title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
# p3 = plot(plotdates, pseudo[m.pseudo_observables[:π_t], :, shock_idx].*4,title="Inflation")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:,exo[shock_name]].*4,title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:,exo[shock_name]],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:ExAnteRealRate],:,exo[shock_name]],title="Ex-ante real rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:,exo[shock_name]],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "Main results/What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# --- Plot side by side the real rate gaps (updated for 3D arrays, no baseline) ---
using Plots

titlelist = ["Ex-ante real rate"]
plotvars = [:ExAnteRealRate]
zeroline = zeros(horizon)
plots_ar1 = Vector{Any}(undef, 1)
shock_idx = exo[:rm_sh]  # 10 for your model
for (i, v) in enumerate(plotvars)
    # Get counterfactual series for this variable from pseudo 3D array
    if v == :pi_t
        cf_series = pseudo[m.pseudo_observables[:π_t], :, shock_idx]
    elseif v in keys(m.observables)
        cf_series = obs[m.observables[v], :, shock_idx]
    elseif v in keys(m.endogenous_states)
        cf_series = states[m.endogenous_states[v], :, shock_idx]
    elseif v in keys(m.pseudo_observables)
        cf_series = pseudo[m.pseudo_observables[v], :, shock_idx]
    else
        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
        plots_ar1[i] = plot(title=string(v), legend=false)
        continue
    end
    cf_series[1] = 0.0 # Align first value to zero
    cf_short = cf_series[end-length(plotdates)+1:end]
    years = unique(year.(plotdates))
    year_tick_dates = [findfirst(d -> year(d) == y, plotdates) !== nothing ? plotdates[findfirst(d -> year(d) == y, plotdates)] : Date(string(y)*"-01-01") for y in years]
    year_tick_labels = [string(y) for y in years]
    xtick_tuple = (year_tick_dates, year_tick_labels)
    # Plot: counterfactual (red), zero line (black)
    p1 = plot(plotdates, cf_short, label="Counterfactual HLW r", color=:blue, lw=2, xticks=xtick_tuple, ylim=(-5.5,4), title="Counterfactual HLW Ex-ante real rate")
    plot!(p1, plotdates, zeroline, lc=:black, lw=2, label="")
    plot!(p1, legend=false)
    plots_ar1[i] = p1
end

plt = plot(plots_ar1[1], layout=(1,1), legend=false)
plot!(plt, size=(960,540))

savefig(plt, "Main results/compare_baseline_vs_HLW_change_only_ex_ante_realrate.pdf")

# --- Plot 2x1 inflation and output grid: counterfactual only (updated for 3D arrays, no baseline) ---
titlelist = ["Inflation", "Output"]
plotvars = [:obs_gdpdeflator, :y_t]
zeroline = zeros(horizon)
plots_arr = Vector{Any}(undef, length(plotvars))
for (i, v) in enumerate(plotvars)
    if v == :pi_t
        cf_series = pseudo[m.pseudo_observables[:π_t], :, shock_idx]
    elseif v in keys(m.observables)
        cf_series = obs[m.observables[v], :, shock_idx].*4*10
    elseif v in keys(m.endogenous_states)
        cf_series = states[m.endogenous_states[v], :, shock_idx]
    elseif v in keys(m.pseudo_observables)
        cf_series = pseudo[m.pseudo_observables[v], :, shock_idx]
    else
        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
        plots_arr[i] = plot(title=string(v), legend=false)
        continue
    end
    cf_series[1] = 0.0
    cf_short = cf_series[end-length(plotdates)+1:end]
    years = unique(year.(plotdates))
    year_tick_dates = [findfirst(d -> year(d) == y, plotdates) !== nothing ? plotdates[findfirst(d -> year(d) == y, plotdates)] : Date(string(y)*"-01-01") for y in years]
    year_tick_labels = [string(y) for y in years]
    xtick_tuple = (year_tick_dates, year_tick_labels)
    p = plot(plotdates, cf_short, label="Counterfactual HLW Real Rate Gap", color=:blue, lw=2, title=titlelist[i], xticks=xtick_tuple)
    plot!(p, plotdates, zeroline, lc=:black, lw=2, label="")
    plot!(p, legend=false)
    plots_arr[i] = p
end
while length(plots_arr) < 2
    push!(plots_arr, plot(title="", legend=false))
end
plt2 = plot(plots_arr[1], plots_arr[2], layout=(1,2), legend=false)
plot!(plt2, size=(960,540))
savefig(plt2, "Main results/compare_baseline_vs_HLW_change_only_inflation_and_output.pdf")

