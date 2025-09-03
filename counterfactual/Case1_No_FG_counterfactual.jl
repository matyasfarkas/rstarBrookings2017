using DSGE, Dates, DataFrames,OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics, ModelConstructors
# using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl and the test forecast driver file
# 25 August 2025
#############################



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

# SET UP MODEL SYSTEM FOR SMOOTHING
m1 = Model1010("ss20")
mode_file = joinpath("dsge/output_data/m1010/ss20/estimate/raw" ,  "paramsmode_vint=161223.h5")
specify_mode!(m1, mode_file)
system = DSGE.compute_system(m1)

## Initialize model object
m = Model1010("ss20")

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()),"dsge", "input_data")
saveroot = dirname(@__FILE__())
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "250825")
m <= DSGE.Setting(:use_population_forecast, false)

# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))

df = load_data(m; check_empty_columns = false)



var_name =:Forward5YearRealNaturalRate
shock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh ] # Convenience yield shocks causing the difference

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Targeted rstar change")
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




# Alternative if r* did not increase
desired_path =hlw_rstar.mean[end-20:259] .- hlw_rstar.mean[end-20]  #rstar_diff[end-16:end] # Desired path for the state variable
var_name =:Forward5YearRealNaturalRate

shock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh ] # Convenience yield shocks causing the difference

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Targeted r* difference")
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
