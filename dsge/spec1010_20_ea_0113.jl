using DSGE, ClusterManagers, HDF5, Distributions, DataFrames 
using Distributed


# TO RUN NARRATIVE SYSTEM PRIORS ADD Dataframes to the packages and in DSGE package's src/estimate/posterior.jl function in the section of # Return total log-likelihood, excluding the presample before the return command the following lines:

# logSP_minus = get_setting(m, :logSP_minus)
# logSPweight = get_setting(m, :logSPweight_minus)
# threshold = get_setting(m,:threshold)
# simnum = get_setting(m,:SP_simnum)
# SP_shock = :g_sh
# SP_obs = :obs_consumption
# SP_first = "2020-Q2"
# SP_last = "2021-Q2"
# SP_vint = "2020-Q2"
# simnum = 100
# # Declare containers
# HVD_c_obs_percent = zeros(1)
# logSP = zeros(1)
# logSPweight  =zeros(1)
# HVD_c_obs_percent_sim = zeros(simnum)
# # df = DSGE.load_data(m)
# start_year = 1959
# start_quarter = 3
# end_year = 2023
# end_quarter = 4

# quarters = []
# for year in start_year:end_year
#     for quarter in 1:4
#         if (year == start_year && quarter < start_quarter) || (year == end_year && quarter > end_quarter)
#             continue
#         end
#         push!(quarters, DSGE.quartertodate("$year-Q$quarter"))
#     end
# end

# datenums = quarters
# df = DataFrame(date=datenums)
# variable_names = string.(keys(m.observables))
# for (i, var) in enumerate(variable_names)
#     df[!,var] = data[i, :]
# end
# _, histshocks, _, _ = DSGE.smooth(m, df, system)
# m <= DSGE.Setting(:shockdec_startdate, nothing)
# _, shockobs, _ = DSGE.shock_decompositions(m, system, histshocks)
# try 
#     HVD_c_obs_percent = sum(abs(shockobs[m.observables[SP_obs], (findfirst(row -> row[1] == quartertodate(SP_first), eachrow(df))-1) : (findfirst(row -> row[1] == quartertodate(SP_last), eachrow(df))-1) ,m.exogenous_shocks[SP_shock]])) / sum(sum(abs.(shockobs[m.observables[SP_obs],(findfirst(row -> row[1] == quartertodate(SP_first), eachrow(df))-1) : (findfirst(row -> row[1] == quartertodate(SP_last), eachrow(df))-1) , [i for i in 1:size(shockobs,3) if i != m.exogenous_shocks[SP_shock]]] )))
# catch
#     HVD_c_obs_percent = abs(shockobs[m.observables[SP_obs], findfirst(row -> row[1] == quartertodate(SP_vint), eachrow(df))-1 ,m.exogenous_shocks[SP_shock]])/sum(abs.(shockobs[m.observables[SP_obs], findfirst(row -> row[1] == quartertodate(SP_vint), eachrow(df))-1 , [i for i in 1:size(shockobs,3) if i != m.exogenous_shocks[SP_shock]]] ))
# end
# logSP = log(pdf(DSGE.BetaAlt(0.10, 0.05), HVD_c_obs_percent))
# println("Log system prior is:", logSP)
# # Check if the change in logSP is small

# if abs(logSP - logSP_minus) > threshold
#     println("Change in logSP is significant. Running simulation loop...")
    
#     # Run the simulation loop
#     HVD_c_obs_percent_sim = zeros(simnum, 1)
#         for i = 1:simnum
#         sim_data = DSGE.simulate_observables(m; burnin=500, n_periods=size(df, 1))'
#         sim_df = copy(df)

#         for j in 1:ncol(df)-1
#             non_nan_idx = findall(!isnan, df[:, j+1])
#             sim_df[non_nan_idx, j+1] .= sim_data[non_nan_idx, j]
#         end

#         _, histshocks_sim, _, _ = DSGE.smooth(m, sim_df, system; cond_type = :full)
#         shockstates, shockobssim, shockpseudo = DSGE.shock_decompositions(m, system, histshocks_sim)
#         try 
#             HVD_c_obs_percent_sim[i] = sum(abs(shockobs[m.observables[SP_obs], (findfirst(row -> row[1] == quartertodate(SP_first), eachrow(df))-1) : (findfirst(row -> row[1] == quartertodate(SP_last), eachrow(df))-1) ,m.exogenous_shocks[SP_shock]])) / sum(sum(abs.(shockobs[m.observables[SP_obs],(findfirst(row -> row[1] == quartertodate(SP_first), eachrow(df))-1) : (findfirst(row -> row[1] == quartertodate(SP_last), eachrow(df))-1) , [i for i in 1:size(shockobs,3) if i != m.exogenous_shocks[SP_shock]]] )))
#         catch
#             HVD_c_obs_percent_sim[i] = abs(shockobssim[m.observables[SP_obs], findfirst(row -> row[1] == quartertodate(SP_vint), eachrow(df))-1, m.exogenous_shocks[SP_shock]]) / sum(abs.(shockobssim[m.observables[SP_obs], findfirst(row -> row[1] == quartertodate(SP_vint), eachrow(df))-1, [i for i in 1:size(shockobssim, 3) if i != m.exogenous_shocks[SP_shock]]]))
#         end
# end

#     logSPweight = mean(log.(pdf(DSGE.BetaAlt(0.10, 0.05), HVD_c_obs_percent_sim)))
#     println("Log system prior weight is:", logSPweight)
# else
#     println("Change in logSP is smaller than ", threshold, ". Skipping simulation loop...")
#     logSPweight =  try get_setting(m, :logSPweight_minus) catch; return end
#     # Reuse the previous weight
# end
# # Update global variables
# m <= DSGE.Setting(:logSP_minus,  logSP)
# m <= DSGE.Setting(:logSPweight_minus, logSPweight)

# return ψ_l * sum(filter_likelihood(m, data, system;
#                                    include_presample = false, tol = tol)) +
#                                        ψ_p * penalty + logSPweight + logSP


##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = true 
run_modal_forecast = true 
run_full_forecast  = false

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "250113")
m <= DSGE.Setting(:use_population_forecast, false)

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC

m <= DSGE.Setting(:reoptimize, true)
m <= DSGE.Setting(:calculate_hessian, true)
m <= DSGE.Setting(:date_mainsample_start,  quartertodate("1970-Q2"))
m <= DSGE.Setting(:date_presample_start,  quartertodate("1970-Q2"))

# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q3"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))

m <= DSGE.Setting(:forecast_block_size,  50)
nworkers = 20
addprocsfcn = addprocs_sge # choose to work with your scheduler; see ClusterManagers.jl


m <= DSGE.Setting(:logSP_minus,  0.0)
m <= DSGE.Setting(:logSPweight_minus, 0.0)
m <= DSGE.Setting(:threshold,  10^-4)
m <= DSGE.Setting(:SP_simnum, 1000)


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
        mode_file = joinpath(dataroot, "user", "paramsmode_vint=250113.h5")
        specify_mode!(m, mode_file)
    end

    # Use calculated ss18 hessian
    if !calculate_hessian(m)
        hessian_file = joinpath(dataroot, "user", "hessian_vint=250113.h5")
        DSGE.specify_hessian!(m, hessian_file)
    end
    df = DSGE.load_data(m)
    data = df_to_matrix(m, df)
    
    max_attempts = 500
    attempt = 1
    estimation_ran = false
    while attempt <= max_attempts && !estimation_ran
        try
            DSGE.estimate(m, data; verbose=:low)
            estimation_ran = true
        catch e
            println("Estimation failed on attempt $attempt: $e. Retrying...")
            attempt += 1
        end
    end
    if estimation_ran==false
        error("Estimation failed after $max_attempts attempts.")
    end

    # Print tables of estimated parameter moments
    groupings = DSGE.parameter_groupings(m)
    moment_tables(m, groupings = groupings)
end



# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC

m <= DSGE.Setting(:reoptimize, false)
m <= DSGE.Setting(:calculate_hessian, false)
m <= DSGE.Setting(:date_mainsample_start,  quartertodate("1970-Q2"))
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
        usual_model_forecast(m, :mode, :none, output_vars,     forecast_string = "",                         density_bands = [.5, .6, .68, .7, .8, .9],                         check_empty_columns = false)
    sections = [:estimation, :forecast]
output_vars = [:forecastobs, :forecastpseudo,:shockdecobs, :shockdecpseudo]



    plot_standard_model_packet(m, :mode, :none, output_vars,
                               forecast_string = "",
                               sections = sections)

                                   write_standard_model_packet(m, :mode, :none, output_vars,
                                sections = sections, forecast_string = "")
    moment_tables(m)
using CSV

function save_shock_decomposition_to_csv(m, var, class, input_type, cond_type; forecast_string = "", groups = shock_groupings(m), file_path = "shock_decomposition.csv")
    # Read in MeansBands
    output_vars = [Symbol(prod, class) for prod in [:shockdec, :trend, :dettrend, :hist, :forecast]]
    mbs = map(output_var -> read_mb(m, input_type, cond_type, output_var, forecast_string = forecast_string), output_vars)

    # Prepare the shock decomposition table
    df = DSGE.prepare_means_table_shockdec(mbs[1], mbs[2], mbs[3], var, mb_hist = mbs[4], mb_forecast = mbs[5], detexify_shocks = false, groups = groups)

    # Save to CSV
    CSV.write(file_path, df)
end
save_shock_decomposition_to_csv(m, :obs_gdpdeflator, :obs, :mode, :none; file_path = "inflation_shock_decomposition.csv")
save_shock_decomposition_to_csv(m, :obs_nominalrate, :obs, :mode, :none; file_path = "FFR_shock_decomposition.csv")
save_shock_decomposition_to_csv(m, :obs_gdp, :obs, :mode, :none; file_path = "gdp_shock_decomposition.csv")

save_shock_decomposition_to_csv(m, :Forward5YearRealNaturalRate, :pseudo, :mode, :none; file_path = "rstar_shock_decomposition.csv")
forecast_string =""
cond_type = :none
                table_vars = [:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                :RealNaturalRate, :Forward5YearRealNaturalRate,
                :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                :Forward30YearRealNaturalRate]
write_meansbands_tables_all(m, :mode, cond_type, [:histpseudo], forecast_string = forecast_string,vars = table_vars)


# Forecast step: produces smoothed histories and shock decompositions
if run_modal_forecast || run_full_forecast

    # what do we want to produce?
    output_vars = [:histpseudo, :forecastpseudo]#, :shockdecpseudo]

    # conditional type
    cond_type = :none

    # Forecast label: all forecast output filenames will contain this string
    forecast_string = ""
    m <= DSGE.Setting(:date_mainsample_start,  quartertodate("1985-Q1"))
    m <= DSGE.Setting(:date_presample_start,  quartertodate("1987-Q1"))
    
    # Modal forecast
    if run_modal_forecast
        # run modal forecasts and save all draws
        forecast_one(m, :mode, cond_type, output_vars; verbose = :low)

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

        DSGE.meansbands_to_matrix(m, :full, cond_type, output_vars; forecast_string = forecast_string)

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
