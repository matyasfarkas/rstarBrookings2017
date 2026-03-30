using DSGE, ClusterManagers, HDF5, Plots, StatsPlots

##########################################################################################
## SETUP
##########################################################################################

# What do you want to do?
run_estimation     = false 
run_modal_forecast = true 

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")
# params_mode = load_draws(m, :mode)
# # Switch off the pi_target shock
# params_mode[m.param_index[:σ_π_target_sh]] = 0.0
# DSGE.update!(m, params_mode)

# Settings for data, paths, etc.
dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "250825")
m <= DSGE.Setting(:use_population_forecast, false)

# Settings for estimation
# set to false => will load pre-computed mode and hessian before MCMC
m <= DSGE.Setting(:reoptimize, false)
m <= DSGE.Setting(:calculate_hessian, false)

# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2025-Q3"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2025-Q3"))

m <= DSGE.Setting(:forecast_block_size,  50)
m <= DSGE.Setting(:optimization_iterations, 100,"Number of iterations the optimizer should run for")
m <= DSGE.Setting(:n_mh_simulations, 10000,"Number of draws saved (after thinning) per block in Metropolis-Hastings")
m <= DSGE.Setting(:mh_adaptive_accpt, false,"Whether to use adaptive acceptance rate in Metropolis-Hastings")
m <= DSGE.Setting(:n_mh_blocks, 2,"Number of blocks for Metropolis-Hastings")
m <= DSGE.Setting(:mh_c, 0.75,"Step size used for adaptive acceptance rate in Metropolis-Hastings")
m <= DSGE.Setting(:n_mh_burn, 1,"Number of blocks to use as burn-in in Metropolis-Hastings")
m <= DSGE.Setting(:mh_thin, 5,"Metropolis-Hastings thinning step")
m <= DSGE.Setting(:mh_cc, 0.09,"Jump size for Metropolis-Hastings (after initialization)")
m <= DSGE.Setting(:mh_cc0, 0.01,"Jump size for initialization of Metropolis-Hastings")
m <= DSGE.Setting(:mh_α, 1.0,"Mixture proportion for adaptive acceptance rate in Metropolis-Hastings")
m <= DSGE.Setting(:shockdec_startdate,  quartertodate("1990-Q1"))

nworkers = 20
addprocsfcn = addprocs_sge # choose to work with your scheduler; see ClusterManagers.jl

df = load_data(m; check_empty_columns = false)
data = df_to_matrix(m, df)
    if !calculate_hessian(m)
        hessian_file = joinpath(saveroot, "output_data", "m1010", "ss20", "estimate", "raw", "hessian_vint=250825.h5")
        DSGE.specify_hessian!(m, hessian_file)
    end

if reoptimize(m)
    estimate(m, data; verbose=:low)
    groupings = DSGE.parameter_groupings(m)
    moment_tables(m, groupings = groupings)
else
    params_mode = load_draws(m, :mode)
    DSGE.update!(m, params_mode)
    DSGE.steadystate!(m)
end

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
output_vars = [:forecastobs, :forecastpseudo,:shockdecobs, :shockdecpseudo, :histstdshocks]
# plot_standard_model_packet(m, :mode, :none, output_vars,
#                                forecast_string = "",
#                                sections = sections)
# write_standard_model_packet(m, :mode, :none, output_vars,
#                                 sections = sections, forecast_string = "")                                
usual_model_forecast(m, :mode, :none, output_vars,     forecast_string = "",                         density_bands = [.5, .6, .68, .7, .8, .9],                         check_empty_columns = false)
# sections = [:estimation, :forecast]
# output_vars = [:forecastobs, :forecastpseudo,:shockdecobs, :shockdecpseudo]
# plot_standard_model_packet(m, :mode, :none, output_vars,
#                                forecast_string = "",
#                                sections = sections)
# write_standard_model_packet(m, :mode, :none, output_vars,
#                                 sections = sections, forecast_string = "")                                
# moment_tables(m)


cond_type = :none
forecast_string =""

forecast_one(m, :mode, cond_type, output_vars; verbose = :high)

# compute means and bands
compute_meansbands(m, :mode, cond_type, output_vars)

                # print history means and bands tables to csv
shockdec_vars = [:π_t, :y_t, :rm_t, :rm_tl1,:rm_tl2,:rm_tl3,:rm_tl4,:rm_tl5,:rm_tl6,:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                :RealRateGap,:Forward5YearRateGap,:ExpectedAvg5YearRateGap,:RealNaturalRate, :Forward5YearRealNaturalRate,
                :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                :Forward30YearRealNaturalRate]

                
write_meansbands_tables_all(m, :mode, cond_type, [:histpseudo,:shockdecpseudo], forecast_string = forecast_string,
                              vars = shockdec_vars)



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

shockdec_vars = [:π_t, :y_t, :rm_t, :rm_tl1,:rm_tl2,:rm_tl3,:rm_tl4,:rm_tl5,:rm_tl6,:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                :RealRateGap,:Forward5YearRateGap,:ExpectedAvg5YearRateGap,:RealNaturalRate, :Forward5YearRealNaturalRate,
                :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                :Forward30YearRealNaturalRate]

                DSGE.write_meansbands_tables_all(m, :mode, cond_type, [:shockdecpseudo, :trendpseudo, :dettrendpseudo],
                                        vars = shockdec_vars,
                                        forecast_string = forecast_string)

# save_shock_decomposition_to_csv(m, :obs_gdpdeflator, :obs, :mode, :none; file_path = "wFG_inflation_shock_decomposition.csv")
# save_shock_decomposition_to_csv(m, :obs_nominalrate, :obs, :mode, :none; file_path = "wFG_FFR_shock_decomposition.csv")
# save_shock_decomposition_to_csv(m, :obs_gdp, :obs, :mode, :none; file_path = "wFG_gdp_shock_decomposition.csv")
# save_shock_decomposition_to_csv(m, :Forward5YearRealNaturalRate, :pseudo, :mode, :none; file_path = "wFG_rstar_shock_decomposition.csv")
forecast_string =""
cond_type = :none
shockdec_vars = [:π_t, :y_t, :rm_t, :rm_tl1,:rm_tl2,:rm_tl3,:rm_tl4,:rm_tl5,:rm_tl6,:ExAnteRealRate, :Forward5YearRealRate, :Forward10YearRealRate,
                :RealRateGap,:Forward5YearRateGap,:ExpectedAvg5YearRateGap,:RealNaturalRate, :Forward5YearRealNaturalRate,
                :Forward10YearRealNaturalRate, :Forward20YearRealNaturalRate,
                :Forward30YearRealNaturalRate]
write_meansbands_tables_all(m, :mode, cond_type, [:histpseudo], forecast_string = forecast_string,vars = shockdec_vars)



shockdec_vars = [:obs_gdpdeflator, :obs_nominalrate, :obs_gdp,:obs_corepce]

DSGE.write_meansbands_tables_all(m, :mode, cond_type, [:shockdecobs, :trendobs, :dettrendobs],
                                        vars = shockdec_vars,
                                        forecast_string = forecast_string)


################################################################################
# Write historical standardized monetary policy shocks to CSV
################################################################################
using JLD2, CSV, DataFrames, Dates, Statistics

input_type = :mode
cond_type  = :none

histstd_file = DSGE.get_forecast_filename(m, input_type, cond_type, :histstdshocks)

# Load saved histstdshocks object
histstd_obj = JLD2.jldopen(histstd_file, "r") do f
    Dict(
        "arr"           => f["arr"],
        "shock_indices" => f["shock_indices"],
        "date_indices"  => f["date_indices"]
    )
end

A             = histstd_obj["arr"]
shock_indices = histstd_obj["shock_indices"]
date_inds     = histstd_obj["date_indices"]

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Get shock index whether keys are Symbols or Strings
function shock_index(shock_indices, sh::Symbol)
    sh_str = String(sh)

    if haskey(shock_indices, sh)
        return shock_indices[sh]
    elseif haskey(shock_indices, sh_str)
        return shock_indices[sh_str]
    else
        return nothing
    end
end

# Convert date metadata into sorted dates and time indices
# Works for Dict{Date,Int} or Dict{Int,Date}
function sorted_dates_and_tinds(date_inds)
    ks = collect(keys(date_inds))
    vs = collect(values(date_inds))

    if !isempty(ks) && first(ks) isa Date
        dates = sort(ks)
        tinds = [date_inds[d] for d in dates]
        return dates, tinds
    elseif !isempty(vs) && first(vs) isa Date
        pairs_vec = sort(collect(pairs(date_inds)), by = x -> x.first)
        tinds = [p.first for p in pairs_vec]
        dates = [p.second for p in pairs_vec]
        return dates, tinds
    else
        error("date_indices does not appear to contain Dates.")
    end
end

# ------------------------------------------------------------------------------
# Collect monetary-policy shocks present in this specification
# ------------------------------------------------------------------------------

wanted_shocks = Symbol[:rm_sh]

for i in 1:DSGE.n_mon_anticipated_shocks(m)
    push!(wanted_shocks, Symbol("rm_shl$i"))
end

# Optional AIT anticipated shocks
try
    if haskey(m.settings, :add_ait_rm) && get_setting(m, :add_ait_rm)
        for i in DSGE.mon_anticipated_ait_shocks(m)
            push!(wanted_shocks, Symbol("rm_ait_shl$i"))
        end
    end
catch
    # skip if not relevant for this spec
end

# Keep only shocks that actually exist in the saved file
wanted_shocks = [sh for sh in wanted_shocks if !isnothing(shock_index(shock_indices, sh))]

if isempty(wanted_shocks)
    error("No monetary-policy shocks found in shock_indices.")
end

# ------------------------------------------------------------------------------
# Build dataframe
# ------------------------------------------------------------------------------

dates, tinds = sorted_dates_and_tinds(date_inds)
df = DataFrame(date = dates)

for sh in wanted_shocks
    sidx = shock_index(shock_indices, sh)

    series = if ndims(A) == 2
        # mode run: nshocks × T
        vec(A[sidx, tinds])
    elseif ndims(A) == 3
        # full run: ndraws × nshocks × T
        vec(dropdims(median(A[:, sidx, tinds], dims=1), dims=1))
    else
        error("Unexpected dimension of histstdshocks array: $(size(A))")
    end

df[!, sh] = series
end

# ------------------------------------------------------------------------------
# Write CSV
# ------------------------------------------------------------------------------

csv_path = joinpath(saveroot, "Final Paper", "US_histstd_monetary_shocks.csv")
CSV.write(csv_path, df)

println("Saved historical standardized monetary shocks to: ", csv_path)
println("Shocks written: ", join(string.(wanted_shocks), ", "))

# ##########################################################################################
# ## RUN
# ##########################################################################################
# using DSGE, ClusterManagers, HDF5, Plots, StatsPlots
# # What do you want to do?
# run_estimation     = true 
# run_modal_forecast = false 
# run_full_forecast  = true

# m = Model1010("ss20")
# # Settings for data, paths, etc.
# dataroot = joinpath(dirname(@__FILE__()), "input_data")
# saveroot = dirname(@__FILE__())
# m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
# m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
# m <= DSGE.Setting(:data_vintage, "250825")
# m <= DSGE.Setting(:use_population_forecast, false)
# m <= DSGE.Setting(:reoptimize, false)
# m <= DSGE.Setting(:calculate_hessian, false)

# # Settings for forecast dates
# m <= DSGE.Setting(:date_forecast_start,  quartertodate("2025-Q3"))
# m <= DSGE.Setting(:date_conditional_end, quartertodate("2025-Q3"))

# m <= DSGE.Setting(:forecast_block_size,  1000)
# m <= DSGE.Setting(:n_mh_blocks, 10,"Number of blocks for Metropolis-Hastings")

# # Run estimation
# if run_estimation

#     if reoptimize(m)
#         # Start from ss20 mode
#         mode_file = rawpath(m, "estimate", "paramsmode.h5")
#         #mode_file = replace(mode_file, "ss20", "ss18")
#         DSGE.update!(m, h5read(mode_file, "params"))
#     else
#         # Use calculated ss20 mode
        mode_file = joinpath(dataroot, "user", "paramsmode_vint=250825.h5")
        specify_mode!(m, mode_file)
    end

    # Use calculated hessian
    if !calculate_hessian(m)
        hessian_file = joinpath(saveroot, "output_data", "m1010", "ss20", "estimate", "raw", "hessian_vint=250825.h5")
        DSGE.specify_hessian!(m, hessian_file)
    end
    df = DSGE.load_data(m,try_disk = true, check_empty_columns = false, summary_statistics = :none)
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
                table_vars = [:Forward5YearRealNaturalRate]
  
                write_meansbands_tables_all(m, :mode, cond_type, [:histpseudo], forecast_string = forecast_string,
                              vars = table_vars)

    end

    # Full-distribution forecast
    if run_full_forecast
        #my_procs = DSGE.addprocsfcn(nworkers)
        ClusterManagers.@everywhere using DSGE

        DSGE.forecast_one(m, :full, cond_type, output_vars; verbose = :high, forecast_string = forecast_string,check_empty_columns = false)
        rstar_bands = [0.68, 0.95]
        DSGE.compute_meansbands(m, :full, cond_type, output_vars; verbose = :high, density_bands = rstar_bands,
                           forecast_string = forecast_string)
        #rmprocs(my_procs)

        DSGE.meansbands_to_matrix(m, :full, cond_type, output_vars; forecast_string = forecast_string)

        # print history means and bands tables to csv
        table_vars = [:Forward5YearRealNaturalRate]
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