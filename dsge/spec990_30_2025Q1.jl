using DSGE, ClusterManagers, HDF5

using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics

cd("C:/Matyi/rstarBrookings2017-master/rstarBrookings2017-master")
##########################################################################################
## SETUP
##########################################################################################

# Initialize model object
m = Model990()
# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= DSGE.Setting(:data_vintage, "250331")
m <= DSGE.Setting(:date_forecast_start, quartertodate("2025-Q1"))
m <= DSGE.Setting(:n_mh_simulations, 500) # Do 500 MH steps during estimation
m <= DSGE.Setting(:n_mh_blocks, 10) # Do 10 blocks

df = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
        # Define the COVID19 date range
        start_date = Date(2020, 3, 31)
        end_date = Date(2020, 9, 30) 
        # Replace rows within the date range with missing values
        allowmissing!(df)
        df[(df.date .>= start_date) .& (df.date .<= end_date), 2:end] .= missing
    data = df_to_matrix(m, df)
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
# Save mode
params = params_mode

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

 CSV.write("statespace_vanila_DNGS/990init_TTT.csv",DataFrame(TTT,state_labels))
 CSV.write("statespace_vanila_DNGS/990init_RRR.csv",DataFrame(RRR,shock_labels))
 CSV.write("statespace_vanila_DNGS/990init_CCC.csv",DataFrame(CCC',state_labels))
 CSV.write("statespace_vanila_DNGS/990init_QQ.csv",DataFrame(QQ,shock_labels))
 CSV.write("statespace_vanila_DNGS/990init_ZZ.csv",DataFrame(ZZ',obs_labels))
 CSV.write("statespace_vanila_DNGS/990init_DD.csv",DataFrame(DD',obs_labels))
 CSV.write("statespace_vanila_DNGS/990init_EE.csv",DataFrame(EE,obs_labels))


kkk=1
for i in 50050:50:100000
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

 CSV.write("statespace_vanila_DNGS/" * string(kkk) * "_TTT.csv",DataFrame(TTT,state_labels))
 CSV.write("statespace_vanila_DNGS/" * string(kkk) *"_RRR.csv",DataFrame(RRR,shock_labels))
 CSV.write("statespace_vanila_DNGS/" * string(kkk) *"_CCC.csv",DataFrame(CCC',state_labels))
 CSV.write("statespace_vanila_DNGS/" * string(kkk) *"_QQ.csv",DataFrame(QQ,shock_labels))
 CSV.write("statespace_vanila_DNGS/" * string(kkk) *"_ZZ.csv",DataFrame(ZZ',obs_labels))
 CSV.write("statespace_vanila_DNGS/" * string(kkk) *"_DD.csv",DataFrame(DD',obs_labels))
 CSV.write("statespace_vanila_DNGS/" * string(kkk) *"_EE.csv",DataFrame(EE,obs_labels))
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

# ##########################################################################################
# ## RUN
# ##########################################################################################

# # Run estimation
# if run_estimation

#     if reoptimize(m)
#         # Start from ss18 mode
#         mode_file = rawpath(m, "estimate", "paramsmode.h5")
#         #mode_file = replace(mode_file, "ss20", "ss18")
#         DSGE.update!(m, h5read(mode_file, "params"))
#     else
#         # Use calculated ss18 mode
#         mode_file = joinpath(dataroot, "user", "paramsmode_vint=240324.h5")
#         specify_mode!(m, mode_file)
#     end

#     # Use calculated ss18 hessian
#     if !calculate_hessian(m)
#         hessian_file = joinpath(dataroot, "user", "hessian_vint=240324.h5")
#         specify_hessian(m, hessian_file)
#     end
#     df = DSGE.load_data(m)
#     data = df_to_matrix(m, df)
#     estimate(m, data; verbose=:low)

#     # Print tables of estimated parameter moments
#     groupings = DSGE.parameter_groupings(m)
#     moment_tables(m, groupings = groupings)
# end

# # Forecast step: produces smoothed histories and shock decompositions
# if run_modal_forecast || run_full_forecast

#     # what do we want to produce?
#     output_vars = [:histpseudo, :forecastpseudo]#, :shockdecpseudo]

#     # conditional type
#     cond_type = :none

#     # Forecast label: all forecast output filenames will contain this string
#     forecast_string = ""

#     # Modal forecast
#     if run_modal_forecast
#         # run modal forecasts and save all draws
#         forecast_one(m, :mode, cond_type, output_vars; verbose = :high)

#         # compute means and bands
#         compute_meansbands(m, :mode, cond_type, output_vars)

#                 # print history means and bands tables to csv
#                 table_vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
#                 :RealNaturalRate, :Forward5YearRealNaturalRate,
#                 :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
#                 :Forward30YearRealNaturalRate]
  
#                 write_meansbands_tables_all(m, :mode, cond_type, [:histpseudo], forecast_string = forecast_string,
#                               vars = table_vars)

#     end

#     # Full-distribution forecast
#     if run_full_forecast
#         #my_procs = DSGE.addprocsfcn(nworkers)
#         ClusterManagers.@everywhere using DSGE

#         DSGE.forecast_one(m, :full, cond_type, output_vars; verbose = :high, forecast_string = forecast_string)
#         rstar_bands = [0.68, 0.95]
#         DSGE.compute_meansbands(m, :full, cond_type, output_vars; verbose = :high, density_bands = rstar_bands,
#                            forecast_string = forecast_string)
#         #rmprocs(my_procs)

#         DSGE.meansbands_to_matrix(m, :full, cond_type, output_vars; forecast_string = forecast_string)

#         # print history means and bands tables to csv
#         table_vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
#                       :RealNaturalRate, :Forward5YearRealNaturalRate,
#                       :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
#                       :Forward30YearRealNaturalRate]
#         DSGE.write_meansbands_tables_all(m, :full, cond_type, [:histpseudo], forecast_string = forecast_string,
#                                     vars = table_vars)

#         # print shockdec means and bands tables to csv
#         if any(x->contains(string(x), "shockdec"), output_vars)
#             shockdec_vars = [:RealNaturalRate, :Forward30YearRealNaturalRate]

#             DSGE.write_meansbands_tables_all(m, :full, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo],
#                                         vars = shockdec_vars,
#                                         forecast_string = forecast_string)

#         end
#     end
# end

# nothing
