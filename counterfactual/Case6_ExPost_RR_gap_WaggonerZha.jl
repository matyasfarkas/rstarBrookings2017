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
m <= DSGE.Setting(:data_vintage, "250825")
# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))
 mode_file = rawpath(m, "estimate", "paramsmode.h5")
        #mode_file = replace(mode_file, "ss20", "ss18")
        DSGE.update!(m, h5read(mode_file, "params"))


system = DSGE.compute_system(m)



nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)

##############
# Function Setup
##############


function obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, m::AbstractDSGEModel, var_name::Symbol, shock_inds::Matrix{Int}, system::System{Float64})
    # shock_inds: matrix of size (nshocks, horizon), each column is the set of shock indices for that period
    horizon = length(x)
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    s_0 = zeros(nstates)
    shocks = zeros(nshocks, horizon)
    var_names, var_class, peg_ind =
        if var_name in keys(m.endogenous_states)
            m.endogenous_states, :states, m.endogenous_states[var_name]
        elseif var_name in keys(m.observables)
            m.observables, :obs,  m.observables[var_name]
        elseif var_name in keys(m.pseudo_observables)
            m.pseudo_observables, :pseudo,  m.pseudo_observables[var_name]
        else
            error("Variable $var_name not found in endogenous states, observables, or pseudo-observables.")
            return
        end

    for t in 1:horizon
        # Build IRF vector for all shocks at time t
        n_shocks_t = size(shock_inds, 1)
        IRFvec = zeros(n_shocks_t)
        for j in 1:n_shocks_t
            test_shocks = zeros(nshocks, horizon)
            test_shocks[shock_inds[j, t], t] = 1.0
            states, obs, pseudo = forecast(system, s_0, test_shocks)
            if var_class == :states
                IRFvec[j] = states[peg_ind, t]
            elseif var_class == :obs
                IRFvec[j] = obs[peg_ind, t]
            elseif var_class == :pseudo
                IRFvec[j] = pseudo[peg_ind, t]
            end
        end
        # Compute effect of previous shocks
        prev_effect = 0.0
        if t > 1
            prev_shocks = shocks[:, 1:t-1]
            prev_states, prev_obs, prev_pseudo = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
            if var_class == :states
                prev_effect = prev_states[peg_ind, t]
            elseif var_class == :obs
                prev_effect = prev_obs[peg_ind, t]
            elseif var_class == :pseudo
                prev_effect = prev_pseudo[peg_ind, t]
            end
        end
        # Solve for required shocks: IRFvec * shocks_t = x[t] - prev_effect
        # If underdetermined, use least squares
        shocks_t = pinv(IRFvec') * (x[t] - prev_effect)
        for j in 1:n_shocks_t
            shocks[shock_inds[j, t], t] = shocks_t[j]
        end
    end
    return shocks
end

# Example test (assuming you have a model m and system):
# x = [1.0, 2.0]
# shock_inds = [m.exogenous_shocks[:shock1] m.exogenous_shocks[:shock1];
#               m.exogenous_shocks[:shock2] m.exogenous_shocks[:shock2]]
# shocks = obtain_shocks_from_desired_state_path_iterative(x, m, :obs_nominalrate, shock_inds, system)


# DSGE.Settings for data, paths, etc.
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
# desired_path = -desired_path
var_name =:RealRateGap

# hock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh ] # Convenience yield shocks causing the difference
# shock_syms = [  :zp_sh ] # Convenience yield shocks causing the difference
shock_syms =keys(m.exogenous_shocks)
#  collect(Iterators.filter(k -> !startswith(String(k), "rm_sh"), keys(m.exogenous_shocks)))

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:Forward5YearRealNaturalRate,:obs_nominalrate,  :pi_t,  :y_t,:ExAnteRealRate,:RealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Target")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealRateGap],:],title="Real rate gap")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "Main results/What_if_real_rate_gap_change_HLW_using_all_shocks.pdf")   # saves the plot from p as a .pdf vector graphic

# --- Plot side by side the real rate gaps (updated for 3D arrays, no baseline) ---
using Plots

titlelist = ["Real real rate"]
plotvars = [:RealRateGap]
zeroline = zeros(horizon)
plots_ar1 = Vector{Any}(undef, 1)
for (i, v) in enumerate(plotvars)
    # Get counterfactual series for this variable from pseudo 3D array
    if v == :pi_t
        cf_series = pseudo[m.pseudo_observables[:π_t], :]
    elseif v in keys(m.observables)
        cf_series = obs[m.observables[v], :]
    elseif v in keys(m.endogenous_states)
        cf_series = states[m.endogenous_states[v], :]
    elseif v in keys(m.pseudo_observables)
        cf_series = pseudo[m.pseudo_observables[v], :]
    else
        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
        plots_ar1[i] = plot(title=string(v), legend=false)
        continue
    end
    cf_series[1] = 0.0 # Align first value to zero
    cf_short = cf_series[end-length(plotdates)+1:end].*4
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

savefig(plt, "Main results/compare_baseline_vs_HLW_change_allshocks_ex_ante_realrate.pdf")

# --- Plot 2x1 inflation and output grid: counterfactual only (updated for 3D arrays, no baseline) ---
titlelist = ["Inflation", "Output"]
plotvars = [:obs_gdpdeflator, :y_t]
zeroline = zeros(horizon)
plots_arr = Vector{Any}(undef, length(plotvars))
for (i, v) in enumerate(plotvars)
    if v == :pi_t
        cf_series = pseudo[m.pseudo_observables[:π_t], :]
    elseif v in keys(m.observables)
        cf_series = obs[m.observables[v], :]
    elseif v in keys(m.endogenous_states)
        cf_series = states[m.endogenous_states[v], :]
    elseif v in keys(m.pseudo_observables)
        cf_series = pseudo[m.pseudo_observables[v], :]
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
savefig(plt2, "Main results/compare_baseline_vs_HLW_change_only_inflation_and_output_allshocks.pdf")

