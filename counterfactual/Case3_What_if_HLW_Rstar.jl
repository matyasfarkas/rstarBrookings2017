using DSGE, Dates, DataFrames,OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics, ModelConstructors
# using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl and the test forecast driver file
# 25 August 2025
#############################

#############################
# To use:
# Just run in the Julia REPL
# include("run_default.jl")
# Note that the estimation
# step will take 2-3 hours.
#############################

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '/', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")

## Load in HLW real time estiamtes of R*
hlw_rstar = CSV.read(joinpath(basepath, "Main results/DSGE_vs_HLW.csv"), DataFrame)

valid_idx = findall(row -> !ismissing(row[:date]) && !ismissing(row[:HLW]) && !ismissing(row[:mean]), eachrow(hlw_rstar))

diff = hlw_rstar.mean[valid_idx] .- hlw_rstar.HLW[valid_idx]
#print(diff)

##############
# Model Setup
##############
# Instantiate the FRBNY DSGE model object
m = Model1010("ss20")
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "250826")
do_not_run_estimation = true


shock_name = :rm_sh # Select MP to implement the specific path in state variable 
var_name = :obs_nominalrate # Select the targeted state variable
var_value = -1.0  # Select the depth of the path
peg_horizon =6;


system = DSGE.compute_system(m)

function obtain_shock_from_desired_pseudo_obs_value(obs_value::Float64, pseudo_ind::Int,
                                             shock_ind::Int, ZZ::Matrix{Float64},
                                             RRR::Matrix{Float64})
    return obs_value/dot(ZZ[pseudo_ind, :], RRR[:, shock_ind])
end


# Setup - copied from impulse_responses.jl
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
    for t = 1:peg_horizon
            var_value_att = var_value - pseudo[m.pseudo_observables[var_name],t, m.exogenous_shocks[shock_name]]
            shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(var_value_att,
                                                                        var_names[var_name],
                                                                        exo[shock_name],
                                                                        system[:ZZ],
                                                                        system[:RRR])
        end
    
    # Iterate state space forward
    states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)
    end


using Plots
p1 = plot(1:horizon,states[m.endogenous_states[:rm_t],:, m.exogenous_shocks[:rm_sh]],title="Monetary policy shock")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]],title="Policy rate")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]],title="Inflation")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]],title="Output")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh]],title="r* (Forward 5-year real natural rate)")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:RealNaturalRate],:, m.exogenous_shocks[:rm_sh]],title="Real natural rate")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,p5,p6,layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "irf/FTPL_Equilibrium_IRF_Policy_rate_with_MP_shock.pdf")   # saves the plot from p as a .pdf vector graphic


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
