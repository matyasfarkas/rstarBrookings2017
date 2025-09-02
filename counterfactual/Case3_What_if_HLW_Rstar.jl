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
hlw_rstar = CSV.read(joinpath(dataroot, "raw/Laubach_Williams_real_time_estimates.csv"), DataFrame)


##############
# Model Setup
##############
# Instantiate the FRBNY DSGE model object
m = Model1010("ss20")
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "250826")
do_not_run_estimation = true

system = DSGE.compute_system(m)
df = load_data(m; check_empty_columns = false)

## Smooth for R* and structural shocks
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
states_df = Dict{Symbol, DataFrame}()
shocks_df = Dict{Symbol, DataFrame}()
pseudo_df = Dict{Symbol, DataFrame}()
smoother = :durbin_koopman #:hamilton, :koopman, :carter_kohn, 
m <= DSGE.Setting(:forecast_smoother, smoother)
states[smoother], shocks[smoother], pseudo[smoother] =
        DSGE.smooth(m, df, system; draw_states = false)
dates = df.date[end-size(states[smoother],2)+1:end]
mat = states[smoother]'  # transpose to 259×91
states_df[smoother] = DataFrame(hcat(dates, mat), [:date; state_labels])
mat = shocks[smoother]'  # transpose to 259×29
shocks_df[smoother] = DataFrame(hcat(dates, mat), [:date; shock_labels])
mat = pseudo[smoother]'  #
pseudo_df[smoother] = DataFrame(hcat(dates, mat), [:date; pseudo_labels])

## Compute the delta to get HLW R* = :ExpectedAvg5YearRealNaturalRate
rstar_smoothed_obs  = pseudo_df[smoother][:, [:date, :Forward5YearRealNaturalRate]].+system.pseudo_measurement.DD_pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate]]

using Dates

# 1. Extract and clean HLW R* data
hlw_dates = hlw_rstar[6:end, 1]
hlw_rstar_vals = hlw_rstar[6:end, 3]

valid_idx = findall(!ismissing, hlw_dates)
hlw_dates_clean = hlw_dates[valid_idx]
hlw_rstar_vals_clean = hlw_rstar_vals[valid_idx]

# Convert HLW dates to Date type
hlw_dates_parsed = Date.(hlw_dates_clean, dateformat"m/d/y")

# Convert HLW R* values to Float64
hlw_rstar_vals_num = parse.(Float64, hlw_rstar_vals_clean)

# Convert HLW dates to quarter string, e.g., "1960Q1"
hlw_quarters = string.(year.(hlw_dates_parsed)) .* "Q" .* string.(ceil.(Int, month.(hlw_dates_parsed) ./ 3))

hlw_df = DataFrame(quarter = hlw_quarters, hlw_rstar = hlw_rstar_vals_num)

# 2. Extract and clean model R* data
model_dates = rstar_smoothed_obs[:, :date]
model_rstar = rstar_smoothed_obs[:, :Forward5YearRealNaturalRate]

# Convert model dates to quarter string, e.g., "1960Q1"
model_quarters = string.(year.(model_dates)) .* "Q" .* string.(ceil.(Int, month.(model_dates) ./ 3))

model_df = DataFrame(quarter = model_quarters, model_rstar = model_rstar)

# 3. Inner join on quarter
joined = innerjoin(hlw_df, model_df, on = :quarter)

# 4. Compute difference
joined.diff = joined.model_rstar .- joined.hlw_rstar

