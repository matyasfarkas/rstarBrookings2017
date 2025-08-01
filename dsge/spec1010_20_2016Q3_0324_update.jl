using DSGE, ClusterManagers, HDF5

##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = false 
run_modal_forecast = false 
run_full_forecast  = true

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "240324")
m <= DSGE.Setting(:use_population_forecast, false)

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC
m <= DSGE.Setting(:reoptimize, true)
m <= DSGE.Setting(:calculate_hessian, true)

# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2016-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2016-Q4"))

m <= DSGE.Setting(:forecast_block_size,  50)
nworkers = 20
addprocsfcn = addprocs_sge # choose to work with your scheduler; see ClusterManagers.jl

##########################################################################################
## RUN
##########################################################################################

# Run estimation
if run_estimation

    if reoptimize(m)
        # Start from ss18 mode
        mode_file = rawpath(m, "estimate", "paramsmode.h5")
        #mode_file = replace(mode_file, "ss20", "ss18")
        DSGE.update!(m, h5read(mode_file, "params"))
    else
        # Use calculated ss18 mode
        mode_file = joinpath(dataroot, "user", "paramsmode_vint=240324.h5")
        specify_mode!(m, mode_file)
    end

    # Use calculated ss18 hessian
    if !calculate_hessian(m)
        hessian_file = joinpath(dataroot, "user", "hessian_vint=240324.h5")
        specify_hessian(m, hessian_file)
    end
    df = DSGE.load_data(m)
    data = df_to_matrix(m, df)
    estimate(m, data; verbose=:low)

    # Print tables of estimated parameter moments
    groupings = DSGE.parameter_groupings(m)
    moment_tables(m, groupings = groupings)
end

# Forecast step: produces smoothed histories and shock decompositions
if run_modal_forecast || run_full_forecast

    # what do we want to produce?
    output_vars = [:histpseudo, :forecastpseudo]#, :shockdecpseudo]

    # conditional type
    cond_type = :none

    # Forecast label: all forecast output filenames will contain this string
    forecast_string = ""

    # Modal forecast
    if run_modal_forecast
        # run modal forecasts and save all draws
        forecast_one(m, :mode, cond_type, output_vars; verbose = :high)

        # compute means and bands
        compute_meansbands(m, :mode, cond_type, output_vars)

                # print history means and bands tables to csv
                table_vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                :RealNaturalRate, :Forward5YearRealNaturalRate,
                :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                :Forward30YearRealNaturalRate]

                write_meansbands_tables_all(m, :mode, cond_type, [:histpseudo], forecast_string = forecast_string,
                              vars = table_vars)

    end

    # Full-distribution forecast
    if run_full_forecast
        #my_procs = DSGE.addprocsfcn(nworkers)
        ClusterManagers.@everywhere using DSGE

        DSGE.forecast_one(m, :full, cond_type, output_vars; verbose = :high, forecast_string = forecast_string)
        rstar_bands = [0.68, 0.95]
        DSGE.compute_meansbands(m, :full, cond_type, output_vars; verbose = :high, density_bands = rstar_bands,
                           forecast_string = forecast_string)
        #rmprocs(my_procs)

        DSGE.meansbands_to_m  atrix(m, :full, cond_type, output_vars; forecast_string = forecast_string)

        # print history means and bands tables to csv
        table_vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                      :RealNaturalRate, :Forward5YearRealNaturalRate,
                      :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                      :Forward30YearRealNaturalRate]
        DSGE.write_meansbands_tables_all(m, :full, cond_type, [:histpseudo], forecast_string = forecast_string,
                                    vars = table_vars)

        # print shockdec means and bands tables to csv
        if any(x->contains(string(x), "shockdec"), output_vars)
            shockdec_vars = [:RealNaturalRate, :Forward30YearRealNaturalRate]

            DSGE.write_meansbands_tables_all(m, :full, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo],
                                        vars = shockdec_vars,
                                        forecast_string = forecast_string)

        end
    end
end

nothing
