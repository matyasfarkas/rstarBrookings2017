using DSGE, ModelConstructors, Distributed
using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl amd the test forecast driver file
# 09 April 2025
#############################

#############################
# To use:
# Just run in the Julia REPL
# include("run_default.jl")
# Note that the estimation
# step will take 2-3 hours.
#############################

##############
# Model Setup
##############
# Instantiate the FRBNY DSGE model object
m = Model1010("ss20")

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "250331")
m <= Setting(:date_forecast_start, quartertodate("2025-Q1"))

# The following settings ensure that this script runs in
# a short amount of time. To properly estimate and
# forecast, we recommend either using the default settings
# (i.e. comment out the settings below) or
# changing the settings yourself.
m <= Setting(:n_mh_simulations, 500) # Do 500 MH steps during estimation
m <= Setting(:n_mh_blocks, 10) # Do 10 blocks

# m <= Setting(:use_population_forecast, false) # Population forecast not available as data to turn off
m <= Setting(:forecast_block_size, 5) # adjust block size to run on small number of estimations
do_not_run_estimation = false
#############
# Estimation
#############
# Reoptimize parameter vector, compute Hessian at mode, and full posterior
# parameter sampling.
#
# Note some columns will have missing data because not all our data
# is publicly available. By default, `load_data` will error
# if any of the columns in the loaded DataFrame is empty,
# but we turn this feature off by setting the keyword `check_empty_columns = false`.
# Warnings will be still thrown after calling load_data indicating which columns
# are empty. However, estimate will still run when data is missing.
# @time begin
# # Run estimation
# if do_not_run_estimation
#     # Start from full sample mode
#    mode_file = rawpath(m, "estimate", "paramsmode.h5")
#    #mode_file = replace(mode_file, "ss20", "ss18")
#    DSGE.update!(m, h5read(mode_file, "params"))
#    hessian_file = rawpath(m, "estimate", "hessian.h5")
#    DSGE.specify_hessian!(m, hessian_file)

# else
#     df = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
#         # Define the COVID19 date range
#         start_date = Date(2020, 3, 31)
#         end_date = Date(2020, 9, 30) 
#         # Replace rows within the date range with missing values
#         allowmissing!(df)
#         df[(df.date .>= start_date) .& (df.date .<= end_date), 2:end] .= missing
#     data = df_to_matrix(m, df)
#     DSGE.estimate(m, data)

# end
# end

        df = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
        # Define the COVID19 date range
        COVID_start_date = Date(2020, 3, 31)
        COVID_end_date = Date(2020, 9, 30) 
        # Replace rows within the date range with missing values
        allowmissing!(df)
        df[(df.date .>= COVID_start_date) .& (df.date .<= COVID_end_date), 2:end] .= missing

for i_year = 2002:2025

    m <= Setting(:date_forecast_start, quartertodate(string(i_year)*"-Q1"))
    data = df_to_matrix(m, df[df.date.<Date(i_year, 3, 31),:])
    DSGE.estimate(m, data)

# Load posterior draws
params_draws = load_draws(m, :full)
params_mode = load_draws(m, :mode)
shock_labels = [key for (key, _) in sort(collect(m.exogenous_shocks), by = x -> x[2])]

obs_labels = [key for (key, _) in sort(collect(m.observables), by = x -> x[2])]

combined = OrderedDict{Symbol, Int64}()
# First insert all entries from endogenous_states
for (k, v) in m.endogenous_states
combined[k] = v
end

# Then insert entries from endogenous_states_augmented
for (k, v) in m.endogenous_states_augmented
combined[k] = v
end
state_labels = [key for (key, _) in sort(collect(combined), by = x -> x[2])]
kkk=1
for i in 1:4:4000
    params = params_draws[kkk,:]

    DSGE.update!(m, params)

    system = compute_system(m, tvis = false)
    # Unpack system
    TTT    = system[:TTT]
    RRR    = system[:RRR]
    CCC    = system[:CCC]
    QQ     = system[:QQ]
    ZZ     = system[:ZZ]
    DD     = system[:DD]
    EE     = system[:EE]

 CSV.write("statespace/" * string(i_year)*"_"* string(kkk) * "_TTT.csv",DataFrame(TTT,state_labels))
 CSV.write("statespace/" * string(i_year)*"_"* string(kkk) *"_RRR.csv",DataFrame(RRR,shock_labels))
 CSV.write("statespace/" * string(i_year)*"_"* string(kkk) *"_CCC.csv",DataFrame(CCC',state_labels))
 CSV.write("statespace/" * string(i_year)*"_"* string(kkk) *"_QQ.csv",DataFrame(QQ,shock_labels))
 CSV.write("statespace/" * string(i_year)*"_"* string(kkk) *"_ZZ.csv",DataFrame(ZZ',obs_labels))
 CSV.write("statespace/" * string(i_year)*"_"* string(kkk) *"_DD.csv",DataFrame(DD',obs_labels))
 CSV.write("statespace/" * string(i_year)*"_"* string(kkk) *"_EE.csv",DataFrame(EE,obs_labels))
kkk = kkk+1
 
#  @inline Φ(s_t1::Vector{S}, ϵ_t::Vector{S}) = TTT*s_t1 + RRR*ϵ_t + CCC
#  @inline Ψ(s_t::Vector{S}) = ZZ*s_t + DD

#  # Define shock and measurement error distributions
#  nshocks = size(QQ, 1)
#  nobs    = size(EE, 1)
#  F_ϵ = Distributions.MvNormal(zeros(nshocks), QQ)
#  F_u = Distributions.MvNormal(zeros(nobs),    EE)

#  return Φ, Ψ, F_ϵ, F_u
end
end

output_vars = [ :histobs, :forecastobs ]

forecast_one(m,:mode, :none, output_vars; df=df, check_empty_columns = false, verbose = :high)
output_files = get_forecast_output_files(m, :mode, :none, output_vars)

cond_type = :none
write_meansbands_tables_all(m, :mode, cond_type, output_vars,vars = [:obs_gdp])


# Load posterior draws
params_mode = load_draws(m, :mode)

DSGE.update!(m, params_mode)
shock_labels = [key for (key, _) in sort(collect(m.exogenous_shocks), by = x -> x[2])]

obs_labels = [key for (key, _) in sort(collect(m.observables), by = x -> x[2])]

combined = OrderedDict{Symbol, Int64}()
# First insert all entries from endogenous_states
for (k, v) in m.endogenous_states
combined[k] = v
end

# Then insert entries from endogenous_states_augmented
for (k, v) in m.endogenous_states_augmented
combined[k] = v
end
state_labels = [key for (key, _) in sort(collect(combined), by = x -> x[2])]
system = compute_system(m, tvis = false)
# Unpack system
TTT    = system[:TTT]
RRR    = system[:RRR]
CCC    = system[:CCC]
QQ     = system[:QQ]
ZZ     = system[:ZZ]
DD     = system[:DD]
EE     = system[:EE]


CSV.write("statespace/mode_TTT.csv",DataFrame(TTT,state_labels))
CSV.write("statespace/mode_RRR.csv",DataFrame(RRR,shock_labels))
CSV.write("statespace/mode_CCC.csv",DataFrame(CCC',state_labels))
CSV.write("statespace/mode_QQ.csv",DataFrame(QQ,shock_labels))
CSV.write("statespace/mode_ZZ.csv",DataFrame(ZZ',obs_labels))
CSV.write("statespace/mode_DD.csv",DataFrame(DD',obs_labels))
CSV.write("statespace/mode_EE.csv",DataFrame(EE,obs_labels))



## DSGE VAR
dsgevar_λ         = 0.5    # What λ do you want to use?
lags = 1

dsgevar= DSGEVAR(m) 
DSGE.update!(dsgevar;observables = collect(keys(m.observables))[1:14], lags = lags, λ = dsgevar_λ)
data = df_to_matrix(m, df)
estimate(dsgevar, data[1:14,130:end]; run_csminwel = true)
ss