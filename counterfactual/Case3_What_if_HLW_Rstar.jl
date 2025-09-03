using DSGE, Dates, DataFrames,OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics, ModelConstructors, LinearAlgebra
# using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl and the test forecast driver file
# 25 August 2025
#############################



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
            m.observables, :obs,  m.observables[var_name]
        elseif var_name in keys(m.pseudo_observables)
            m.pseudo_observables, :pseudo,  m.pseudo_observables[var_name]
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
#               m.exogenous_shocks[:shock2] m.exogenous_shocks[:shock2]]
# shocks = obtain_shocks_from_desired_state_path_iterative(x, m, :obs_nominalrate, shock_inds, system)


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

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="HLW r* minus DSGE r*")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "Main results/rstar_what_if_HLW.pdf")   # saves the plot from p as a .pdf vector graphic





# Alternative if r* did not increase post COVID19
desired_path =hlw_rstar.mean[end-20:259] .- hlw_rstar.mean[end-20]  #rstar_diff[end-16:end] # Desired path for the state variable
# desired_path = -desired_path
var_name =:Forward5YearRealNaturalRate

# shock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh ] # Convenience yield shocks causing the difference
shock_syms = [  :zp_sh ] # Convenience yield shocks causing the difference

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Change in r* since end of COVID19")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

savefig( "Main results/rstar_had_not_increased.pdf")   # saves the plot from p as a .pdf vector graphic


# df = load_data(m; check_empty_columns = false)

# ## Smooth for R* and structural shocks
# states = Dict{Symbol, Matrix{Float64}}()
# shocks = Dict{Symbol, Matrix{Float64}}()
# pseudo = Dict{Symbol, Matrix{Float64}}()

# shock_labels = [key for (key, _) in sort(collect(m.exogenous_shocks), by = x -> x[2])]
# combined = DSGE.OrderedDict{Symbol, Int64}()
# # First insert all entries from endogenous_states
# for (k, v) in m.endogenous_states
# combined[k] = v
# end
# # Then insert entries from endogenous_states_augmented
# for (k, v) in m.endogenous_states_augmented
# combined[k] = v
# end
# state_labels = [key for (key, _) in sort(collect(combined), by = x -> x[2])]
# pseudo_labels = [key for (key, _) in sort(collect(m.pseudo_observables), by = x -> x[2])]
# states_df = Dict{Symbol, DataFrame}()
# shocks_df = Dict{Symbol, DataFrame}()
# pseudo_df = Dict{Symbol, DataFrame}()
# smoother = :durbin_koopman #:hamilton, :koopman, :carter_kohn, 
# m <= DSGE.Setting(:forecast_smoother, smoother)
# states[smoother], shocks[smoother], pseudo[smoother] =
#         DSGE.smooth(m, df, system; draw_states = false)
# dates = df.date[end-size(states[smoother],2)+1:end]
# mat = states[smoother]'  # transpose to 259×91
# states_df[smoother] = DataFrame(hcat(dates, mat), [:date; state_labels])
# mat = shocks[smoother]'  # transpose to 259×29
# shocks_df[smoother] = DataFrame(hcat(dates, mat), [:date; shock_labels])
# mat = pseudo[smoother]'  #
# pseudo_df[smoother] = DataFrame(hcat(dates, mat), [:date; pseudo_labels])

# ## Compute the delta to get HLW R* = :ExpectedAvg5YearRealNaturalRate
# rstar_smoothed_obs  = pseudo_df[smoother][:, [:date, :Forward5YearRealNaturalRate]]
