using DSGE, Dates, DataFrames,OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics, ModelConstructors
# using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl and the test forecast driver file
# 25 August 2025
#############################


##############
# Function Setup
##############



function obtain_shocks_from_desired_state_path_iterative(x::Matrix{Float64}, m::AbstractDSGEModel, var_names::Vector{Symbol}, shock_inds::Matrix{Int}, system::System{Float64})
    # x: (ntargets, horizon), var_names: vector of symbols
    horizon = size(x, 2)
    ntargets = size(x, 1)
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    s_0 = zeros(nstates)
    shocks = zeros(nshocks, horizon)
    # For each target, get (var_class, peg_ind)
    var_classes = Vector{Symbol}(undef, ntargets)
    peg_inds = Vector{Int}(undef, ntargets)
    for i in 1:ntargets
        v = var_names[i]
        if v in keys(m.endogenous_states)
            var_classes[i] = :states
            peg_inds[i] = m.endogenous_states[v]
        elseif v in keys(m.observables)
            var_classes[i] = :obs
            peg_inds[i] = m.observables[v]
        elseif v in keys(m.pseudo_observables)
            var_classes[i] = :pseudo
            peg_inds[i] = m.pseudo_observables[v]
        else
            error("Variable $v not found in endogenous states, observables, or pseudo-observables.")
            return
        end
    end

    for t in 1:horizon
        # Build IRF matrix for all shocks and all targets at time t
        n_shocks_t = size(shock_inds, 1)
        IRFmat = zeros(ntargets, n_shocks_t)
        for j in 1:n_shocks_t
            test_shocks = zeros(nshocks, horizon)
            test_shocks[shock_inds[j, t], t] = 1.0
            states, obs, pseudo = forecast(system, s_0, test_shocks)
            for i in 1:ntargets
                if var_classes[i] == :states
                    IRFmat[i, j] = states[peg_inds[i], t]
                elseif var_classes[i] == :obs
                    IRFmat[i, j] = obs[peg_inds[i], t]
                elseif var_classes[i] == :pseudo
                    IRFmat[i, j] = pseudo[peg_inds[i], t]
                end
            end
        end
        # Compute effect of previous shocks for all targets
        prev_effect = zeros(ntargets)
        if t > 1
            prev_shocks = shocks[:, 1:t-1]
            prev_states, prev_obs, prev_pseudo = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
            for i in 1:ntargets
                if var_classes[i] == :states
                    prev_effect[i] = prev_states[peg_inds[i], t]
                elseif var_classes[i] == :obs
                    prev_effect[i] = prev_obs[peg_inds[i], t]
                elseif var_classes[i] == :pseudo
                    prev_effect[i] = prev_pseudo[peg_inds[i], t]
                end
            end
        end
        # Solve for required shocks: IRFmat * shocks_t = x[:, t] - prev_effect
        # If underdetermined, use least squares
        shocks_t = pinv(IRFmat) * (x[:, t] - prev_effect)
        for j in 1:n_shocks_t
            shocks[shock_inds[j, t], t] = shocks_t[j]
        end
    end
    return shocks
end

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")

## Load in HLW real time estiamtes of R*
csv_path = joinpath(basepath, "Main results", "DSGE_vs_HLW.csv")
hlw_rstar = DataFrame(CSV.File(csv_path))

valid_idx = findall(row -> !ismissing(row[:date]) && !ismissing(row[:HLW]) && !ismissing(row[:mean]), eachrow(hlw_rstar))

rstar_diff = hlw_rstar.HLW[valid_idx] .- hlw_rstar.mean[valid_idx]
dates= hlw_rstar.date[valid_idx]

desired_path = rstar_diff#rstar_diff[end-16:end] # Desired path for the state variable
var_name =:Forward5YearRealNaturalRate

shock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh, :zp_sh ] # Convenience yield shocks causing the difference

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

