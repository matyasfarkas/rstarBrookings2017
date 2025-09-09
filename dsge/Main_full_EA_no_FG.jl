using DSGE, ClusterManagers, HDF5, Plots, StatsPlots

##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = false 
run_modal_forecast = true 
run_full_forecast  = false

# What do you want to do?
run_estimation     = true 
run_modal_forecast = true 
run_full_forecast  = false

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss23")
# This is the Euro Area estimation specification:
# function ss23!(m::Model1010)
#     # ss20 with altered prior for EA dataset. 
#     # SS labour growth is mean zero:
#     m <= parameter(:Lmean, 0.91, (-1000., 1000.), (-1e3, 1e3), ModelConstructors.Untransformed(), Normal(0.91, 1.), fixed=false,
#                    description="Lmean: Mean level of hours.",
#                    tex_label="\\bar{L}")
#     # SS growth rate of technology is 1.2:  
#     m <= parameter(:γ, 0.1823, (-5.0, 5.0), (-5., 5.), ModelConstructors.Untransformed(), Normal(0.1823, 0.1), fixed=false, scaling = x -> x/100,
#                    description="γ: The log of the steady-state growth rate of technology.",
#                    tex_label="100\\gamma")
# end
# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "250114")
m <= DSGE.Setting(:use_population_forecast, false)

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC

m <= DSGE.Setting(:reoptimize, false)
m <= DSGE.Setting(:calculate_hessian, false)
m <= DSGE.Setting(:date_mainsample_start,  quartertodate("1972-Q2"))
m <= DSGE.Setting(:date_presample_start,  quartertodate("1970-Q2"))

# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q3"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))


df = load_data(m; check_empty_columns = false)
output_vars = Vector{Symbol}(undef,0)
do_histforecast= true
do_shockdecs   = true
    if do_histforecast
        # Write data to create historical and forecast output
        output_vars = vcat(output_vars, [:histpseudo, :histobs, :histstdshocks,
                                         :hist4qpseudo, :hist4qobs, :histutpseudo,
                                         :forecastpseudo, :forecastobs, :forecastutpseudo,
                                         :forecast4qpseudo, :forecast4qobs, :forecaststdshocks])
    end

    if do_shockdecs
        # Shock decompositions of forecasts
        output_vars = vcat(output_vars, [:dettrendobs, :dettrendpseudo, :trendobs,
                                         :trendpseudo, :shockdecpseudo, :shockdecobs])
    end
    
# usual_model_forecast(m, :mode, :none, output_vars,     forecast_string = "",                         density_bands = [.5, .6, .68, .7, .8, .9],                         check_empty_columns = false)
# sections = [:estimation, :forecast]
# output_vars = [:forecastobs, :forecastpseudo,:shockdecobs, :shockdecpseudo]
# plot_standard_model_packet(m, :mode, :none, output_vars,
#                                forecast_string = "",
#                                sections = sections)
#    write_standard_model_packet(m, :mode, :none, output_vars,
                                # sections = sections, forecast_string = "")
# moment_tables(m)


# using CSV

# function save_shock_decomposition_to_csv(m, var, class, input_type, cond_type; forecast_string = "", groups = shock_groupings(m), file_path = "shock_decomposition.csv")
#     # Read in MeansBands
#     output_vars = [Symbol(prod, class) for prod in [:shockdec, :trend, :dettrend, :hist, :forecast]]
#     mbs = map(output_var -> read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string), output_vars)

#     # Prepare the shock decomposition table
#     df = DSGE.prepare_means_table_shockdec(mbs[1], mbs[2], mbs[3], var, mb_hist = mbs[4], mb_forecast = mbs[5], detexify_shocks = false, groups = groups)

#     # Save to CSV
#     CSV.write(file_path, df)
# end
# save_shock_decomposition_to_csv(m, :obs_gdpdeflator, :obs, :mode, :none; file_path = "woFG_inflation_shock_decomposition.csv")
# save_shock_decomposition_to_csv(m, :obs_nominalrate, :obs, :mode, :none; file_path = "woFG_FFR_shock_decomposition.csv")
# save_shock_decomposition_to_csv(m, :obs_gdp, :obs, :mode, :none; file_path = "woFG_gdp_shock_decomposition.csv")
# save_shock_decomposition_to_csv(m, :Forward5YearRealNaturalRate, :pseudo, :mode, :none; file_path = "woFG_rstar_shock_decomposition.csv")

forecast_one(m, :mode, :none, output_vars; verbose = :high , check_empty_columns = false)

forecast_string =""
cond_type = :none
                table_vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                :RealNaturalRate, :Forward5YearRealNaturalRate,
                :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                :Forward30YearRealNaturalRate]
write_meansbands_tables_all(m, :mode, cond_type, [:histpseudo], forecast_string = forecast_string,vars = table_vars)

shockdec_vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                :RealNaturalRate, :Forward5YearRealNaturalRate,
                :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                :Forward30YearRealNaturalRate]

DSGE.write_meansbands_tables_all(m, :mode, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo],
                                        vars = shockdec_vars,
                                        forecast_string = forecast_string)

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
#     df = DSGE.load_data(m,try_disk = false, check_empty_columns = false, summary_statistics = :none)
#     data = df_to_matrix(m, df)
#     estimate(m, data; verbose=:low)

#     # Print tables of estimated parameter moments
#     groupings = DSGE.parameter_groupings(m)
#     moment_tables(m, groupings = groupings)
# end

# # Forecast step: produces smoothed histories and shock decompositions
# if run_modal_forecast || run_full_forecast

#     # what do we want to produce?
#     output_vars = [:histpseudo, :forecastpseudo,:histobs]#, :shockdecpseudo]

#     # conditional type
#     cond_type = :none

#     # Forecast label: all forecast output filenames will contain this string
#     forecast_string = ""

#     # Modal forecast
#     if run_modal_forecast
#         # run modal forecasts and save all draws
#         forecast_one(m, :mode, cond_type, output_vars; verbose = :high,check_empty_columns = false)

#         # compute means and bands
#         compute_meansbands(m, :mode, cond_type, output_vars)


#                 # print history means and bands tables to csv
#                 table_vars = [:obs_nominalrate,:obs_gdp,:obs_longrate]
  
#                 write_meansbands_tables_all(m, :mode, cond_type, [:histobs], forecast_string = forecast_string,
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
