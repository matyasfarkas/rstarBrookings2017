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
@time begin
# Run estimation
if do_not_run_estimation
    # Start from full sample mode
   mode_file = rawpath(m, "estimate", "paramsmode.h5")
   #mode_file = replace(mode_file, "ss20", "ss18")
   DSGE.update!(m, h5read(mode_file, "params"))
   hessian_file = rawpath(m, "estimate", "hessian.h5")
   DSGE.specify_hessian!(m, hessian_file)

else
    df = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
        # Define the COVID19 date range
        start_date = Date(2020, 3, 31)
        end_date = Date(2020, 9, 30) 
        # Replace rows within the date range with missing values
        allowmissing!(df)
        df[(df.date .>= start_date) .& (df.date .<= end_date), 2:end] .= missing
    data = df_to_matrix(m, df)
    DSGE.estimate(m, data)

end
end
# produce LaTeX tables of parameter moments
# moment_tables(m)

############
# Forecasts
############
# (1) forecast_one produces forecasts of the variables specified in output_vars
# and stores the forecasts in model units (log deviations from steady state).
# (2) compute_meansbands stores the forecast in a MeansBands type, defined in DSGE.jl,
# which transforms the forecast from model units to "human readable" units, i.e.
# the GDP forecast is transformed from annualized quarterly log per-capita growth rates
# to annualized quarterly aggregate percent change.

# The variables we want to forecast. In this case, all of the model observables
output_vars = [ :histobs, :forecastobs ]

# Modal forecast (point forecast)
forecast_one(m, :mode, :none, output_vars; check_empty_columns = false)
m <= Setting(:use_population_forecast, false) # disable population forecasts
compute_meansbands(m, :mode, :none, output_vars; check_empty_columns = false)

#  Modal conditional forecast
# m <= Setting(:cond_vintage, "250331")
# forecast_one(m, :mode, :semi, output_vars; check_empty_columns = false)
# compute_meansbands(m, :mode, :semi, output_vars; check_empty_columns = false)
# table_vars = [:obs_nominalrate]
# write_meansbands_tables_all(m, :mode, :semi, [:forecastobs], forecast_string = "", vars = table_vars)


# Full-distribution forecast (point forecast (mean) and uncertainty bands)
# Optionally add 10 processes to run the forecast in parallel (uncomment the 3 lines below).
# Alternatively, you can load the ClusterManagers package and add processes
# using one of the schedulers such as SGE or Slurm.
# addprocs(4)
# @everywhere using DSGE
# m <= Setting(:use_parallel_workers, true)

# Comment out the line below if you did not run the forecast in parallel.
# rmprocs(procs())
# plot_history_and_forecast(m, :obs_nominalrate, :obs, :mode, :semi,
#                               use_bdd = :bdd_and_unbdd,
#                               start_date = DSGE.quartertodate("2007-Q1"),
#                               end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
#                               verbose = :none)

## Conditional forecast  - works!!!!
# forecast_one(m, :full, :none, output_vars; check_empty_columns = false)
# compute_meansbands(m, :full, :none, output_vars; check_empty_columns = false)
# m <= Setting(:use_population_forecast, true) # disable population forecasts 
# plot_history_and_forecast(m, :obs_gdp, :obs, :mode, :full,
#                               use_bdd = :bdd_and_unbdd,
#                               start_date = DSGE.quartertodate("2007-Q1"),
#                               end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
#                               verbose = :none)

                          
                                                       
## Run the rolling window estimation with pseudo real time out of sample forecasts

finvars = [:date, :obs_nominalrate,:obs_BBBspread, :obs_AAAspread,:obs_longrate ]

#  We define the following financial variables - Based on the m1010 observables in DSGE.jl

    ############################################################################
    ## 1. Nominal short-term interest rate (3 months)
    ############################################################################
        # FROM: Nominal effective federal funds rate (aggregate daily data at a
        #       quarterly frequency at an annual rate)
        # TO:   Nominal effective fed funds rate, at a quarterly rate


    ############################################################################
    ## 2. Spread: BAA-20yr TBill
    ############################################################################
        # FROM: Baa corporate bond yield (percent annualized), and 20-year
        #       treasury note yield (percent annualized)
        # TO:   Baa yield - 20T yield spread at a quarterly rate
        # Note: Moody's corporate bond yields on the H15 are based on corporate
        #       bonds with remaining maturities of at least 20 years.
        # Another note: There were no 20-year treasuries between
        #       1987-1993. Therefore, we use the LTGOVTBD series until 2000, then splice in the GS20 series.

    ############################################################################
    ## 3. Spread: AAA-20yr TBill
    ############################################################################
        # FROM: AAA corporate bond yield (percent annualized), and 20-year
        #       treasury note yield (percent annualized)
        # TO:   AAA yield - 20T yield spread at a quarterly rate
        # Note: Moody's corporate bond yields on the H15 are based on corporate
        #       bonds with remaining maturities of at least 20 years.
        # Another note: There were no 20-year treasuries between
        #       1987-1993. Therefore, we use the LTGOVTBD series until 2000, then splice in the GS20 series.
#  We condition them on the following shocks:  

finshocks = [:b_liqtil_sh, :b_liqp_sh,:b_safetil_sh, :b_safep_sh,:rm_sh] #:γ_sh, :σ_ω_sh,:μ_e_sh, :b_liqtil_sh, :b_liqp_sh,:b_safetil_sh, :b_safep_sh,:rm_sh



    
k =1
for i_year in 2001:2024
    for i_quarter in 1:4

        #  Construct model used for forecasting based on full sample estimated DSGE
                mnc = m
                mnc <= Setting(:data_vintage, "250331")
                m <= Setting(:date_forecast_start, quartertodate("2025-Q1"))
                df_final = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)

        # TIming           
                date_last_obs = quartertodate(string(i_year)*"-Q"*string(i_quarter))
                vint_last_obs = string(date_last_obs)
                vint_last_obs = vint_last_obs[3:4]*vint_last_obs[6:7]*vint_last_obs[9:10]
           
        #  NOWCASTS
            if i_quarter ==4
                date_nc_cond = quartertodate(string(i_year+1)*"-Q"*string(1))
                vint_nc_cond= string(date_nc_cond)
                mnc <= Setting(:date_forecast_start,date_nc_cond)
            else
                date_nc_cond = quartertodate(string(i_year)*"-Q"*string(i_quarter+1))
                vint_nc_cond= string(date_nc_cond)
                mnc <= Setting(:date_forecast_start,date_nc_cond)
            end
            vint_nc_cond = vint_nc_cond[3:4]*vint_nc_cond[6:7]*vint_nc_cond[9:10]

        #  1 YEAR AHEAD
            if i_quarter ==4
                    date_t1_cond = quartertodate(string(i_year+2)*"-Q"*string(1))
                    vint_t1_cond= string(date_t1_cond)
                    mnc <= Setting(:date_forecast_start,date_t1_cond)
            else
                    date_t1_cond = quartertodate(string(i_year+1)*"-Q"*string(i_quarter+1))
                    vint_t1_cond= string(date_t1_cond)
                    mnc <= Setting(:date_forecast_start,date_t1_cond)
            end
            vint_t1_cond = vint_t1_cond[3:4]*vint_t1_cond[6:7]*vint_t1_cond[9:10]
            date_t1_cond_end = quartertodate(string(i_year+1)*"-Q"*string(i_quarter))
        # 2 YEARS AHEAD 
        if i_quarter ==4
            date_t2_cond = quartertodate(string(i_year+3)*"-Q"*string(1))
            vint_t2_cond= string(date_t2_cond)
            mnc <= Setting(:date_forecast_start,date_t2_cond)
        else
            date_t2_cond = quartertodate(string(i_year+2)*"-Q"*string(i_quarter+1))
            vint_t2_cond= string(date_t2_cond)
            mnc <= Setting(:date_forecast_start,date_t2_cond)
        end
        vint_t2_cond = vint_t2_cond[3:4]*vint_t2_cond[6:7]*vint_t2_cond[9:10]
        date_t2_cond_end = quartertodate(string(i_year+2)*"-Q"*string(i_quarter))

       #  Compute unconditional forecasts
                mnc <= Setting(:date_forecast_start,date_last_obs)
                df = load_data(mnc, try_disk = false, check_empty_columns = false, summary_statistics = :none)
                
                forecast_one(mnc, :full, :none, output_vars; df=df, check_empty_columns = false, verbose = :high)
                # plot_history_and_forecast(mnc, finvars[2:end] , :obs, :mode, :none,
                # use_bdd = :bdd_and_unbdd,
                # start_date = DSGE.quartertodate(string(i_year-10)*"-Q"*string(1)),
                # end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
                # verbose = :none)
                
                output_files = get_forecast_output_files(mnc, :full, :none, output_vars)
                out =  Dict{Symbol, Array{Float64}}()
                for var in keys(output_files)
                    out[var] = FileIO.load(output_files[var], "arr")
                end
                # fc = read_forecast_output(mnc,:full,:none,:forecastobs,:obs_gdp)

                temp_indices = [mnc.observables[k] for k in finvars[2:end]]
                if date_t2_cond > df_final.date[end]
                    maxhoriz =size(df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t2_cond),: ],1)
                    t2_delta = cat(out[:forecastobs][:,temp_indices,1:maxhoriz] -  PermutedDimsArray(reshape(repeat(convert(Matrix,df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= df_final.date[end]), finvars[2:end]])', 1, 1,size( out[:forecastobs][:,temp_indices,1:maxhoriz],1)), size( out[:forecastobs][:,temp_indices,1:maxhoriz],2), size( out[:forecastobs][:,temp_indices,1:maxhoriz],3), size( out[:forecastobs][:,temp_indices,1:maxhoriz],1)),(3, 1, 2)), zeros(size( out[:forecastobs][:,temp_indices,1:maxhoriz],1),size(temp_indices,1),9-maxhoriz),dims=3)                    
                else
                    t2_delta = out[:forecastobs][:,temp_indices,1:9] -  PermutedDimsArray(reshape(repeat(convert(Matrix,df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t2_cond), finvars[2:end]])', 1, 1,size( out[:forecastobs][:,temp_indices,1:9],1)), size( out[:forecastobs][:,temp_indices,1:9],2), size( out[:forecastobs][:,temp_indices,1:9],3), size( out[:forecastobs][:,temp_indices,1:9],1)),(3, 1, 2)) # This results in ndraws x ntargets x horizon matrix!
                end
                
                if date_t1_cond > df_final.date[end]
                        maxhoriz =size(df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t1_cond),: ],1)
                        t1_delta = cat(out[:forecastobs][:,temp_indices,1:maxhoriz] -  PermutedDimsArray(reshape(repeat(convert(Matrix,df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= df_final.date[end]), finvars[2:end]])', 1, 1,size( out[:forecastobs][:,temp_indices,1:maxhoriz],1)), size( out[:forecastobs][:,temp_indices,1:maxhoriz],2), size( out[:forecastobs][:,temp_indices,1:maxhoriz],3), size( out[:forecastobs][:,temp_indices,1:maxhoriz],1)),(3, 1, 2)), zeros(size( out[:forecastobs][:,temp_indices,1:maxhoriz],1),size(temp_indices,1),5-maxhoriz),dims=3)                    
                else
                    t1_delta = out[:forecastobs][:,temp_indices,1:5] -  PermutedDimsArray(reshape(repeat(convert(Matrix,df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t1_cond), finvars[2:end]])', 1, 1,size( out[:forecastobs][:,temp_indices,1:5],1)), size( out[:forecastobs][:,temp_indices,1:5],2), size( out[:forecastobs][:,temp_indices,1:5],3), size( out[:forecastobs][:,temp_indices,1:5],1)),(3, 1, 2)) # This results in ndraws x ntargets x horizon matrix!
                end
                
                if date_nc_cond > df_final.date[end]
                    nc_delta =  zeros(size( out[:forecastobs][:,temp_indices,1:maxhoriz],1),size(temp_indices,1),1)                    
                else

                    nc_delta = out[:forecastobs][:,temp_indices,1] - repeat(convert(Matrix,  df_final[(df_final.date .==date_nc_cond), finvars[2:end]]),size(out[:forecastobs][:,temp_indices,1],1),1)
                end
                    # df_f= DataFrame(out[:forecastobs]',names(df[:,2:end]))
            
                # df_target_nc=insertcols!(DataFrame(df_f[1,finvars[2:end]]), 1, :date => date_nc_cond)  .- filter!(row -> row.date == date_nc_cond, df_final[:,finvars[1:end]])
                # df_target_t1=insertcols!(DataFrame(df_f[1:5,finvars[2:end]]), 1, :date => df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t1_cond), 1])  .-  df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t1_cond), finvars]
                # df_target_t2=insertcols!(DataFrame(df_f[1:9,finvars[2:end]]), 1, :date => df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t2_cond), 1])  .-  df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t2_cond), finvars]

                file_path_nc = joinpath(get_setting(mnc, :dataroot), "scenarios", "nowcast_"* vint_last_obs *".jld2")
                file_path_t1 = joinpath(get_setting(mnc, :dataroot), "scenarios", "forecastt1_"* vint_last_obs *".jld2")
                file_path_t2 = joinpath(get_setting(mnc, :dataroot), "scenarios", "forecastt2_"* vint_last_obs *".jld2")
                
                nowcast = Scenario(:nowcast, "Nowcast Scenario",finvars[2:end], finshocks, vint_last_obs ; shock_scaling = 1.0, draw_states = false)
                forecastt1 = Scenario(:forecastt1, "One-Year Ahead Conditional Scenario",finvars[2:end], finshocks, vint_last_obs ; shock_scaling = 1.0, draw_states = false)
                forecastt2 = Scenario(:forecastt2, "Two-Year Ahead Conditional Scenario",finvars[2:end], finshocks, vint_last_obs ; shock_scaling = 1.0, draw_states = false)

                
                target_indices = OrderedDict(s => i for (i, s) in enumerate(finvars[2:end]))

                save(file_path_nc, "arr", nc_delta,"target_indices", target_indices)
                save(file_path_t1, "arr", t1_delta,"target_indices", target_indices)
                save(file_path_t2, "arr", t2_delta,"target_indices", target_indices)

                conditional_forecast_nc = forecast_scenario(mnc, nowcast, verbose =:high)
                conditional_forecast_t1 = forecast_scenario(mnc, forecastt1, verbose =:high)
                conditional_forecast_t2 = forecast_scenario(mnc, forecastt2, verbose =:high)
                # Reconstructing the unconditonal forecast with conditional forecast + delta = unconditional forecast in the key dimensions
                # tempp = conditional_forecast_t1[:forecastobs][:,temp_indices,1:5]  + PermutedDimsArray(reshape(repeat(convert(Matrix,df_final[(df_final.date .>=date_nc_cond) .& (df_final.date .<= date_t1_cond), finvars[2:end]])', 1, 1,size( out[:forecastobs][:,temp_indices,1:5],1)), size( out[:forecastobs][:,temp_indices,1:5],2), size( out[:forecastobs][:,temp_indices,1:5],3), size( out[:forecastobs][:,temp_indices,1:5],1)),(3, 1, 2))
                # sum(tempp - out[:forecastobs][:,temp_indices,1:5])
                # PASS: THIS SHOULD BE ZERO!!!
                
                final_forecast_numbers_nc = conditional_forecast_nc[:forecastobs][:,:,:] + out[:forecastobs][:,:,:]
                final_forecast_numbers_t1 = conditional_forecast_t1[:forecastobs][:,:,:] + out[:forecastobs][:,:,:]
                final_forecast_numbers_t2 = conditional_forecast_t2[:forecastobs][:,:,:] + out[:forecastobs][:,:,:]
                
                table_vars = [:obs_gdp;  :obs_hours;       :obs_wages        ;       :obs_gdpdeflator  ;       :obs_corepce      ;       :obs_nominalrate  ;       :obs_consumption  ;       :obs_investment   ;       :obs_BBBspread    ;       :obs_longinflation;        :obs_longrate     ;        :obs_tfp          ;        :obs_gdi          ;        :obs_AAAspread    ]
                
                    # Quantiles of interest
                    QQ = [0.16, 0.25, 0.5, 0.75, 0.84]
                    quantile_labels = ["q16", "q25", "q50", "q75", "q84"]

                    # Loop through variables
                    for i_var in 1:length(table_vars)
                        # Select the relevant forecast data for the variable across horizons and draws
                        var_fc_nc = final_forecast_numbers_nc[:, i_var, :]
                        var_fc_t1 = final_forecast_numbers_t1[:, i_var, :]
                        var_fc_t2 = final_forecast_numbers_t2[:, i_var, :]

                        # Compute quantiles for each horizon (i.e., over the third dimension)
                        quantiles_nc = reduce(hcat, map(x -> quantile(x, QQ), eachslice(var_fc_nc, dims=2)))'
                        quantiles_t1 = reduce(hcat, map(x -> quantile(x, QQ), eachslice(var_fc_t1, dims=2)))'
                        quantiles_t2 = reduce(hcat, map(x -> quantile(x, QQ), eachslice(var_fc_t2, dims=2)))'

                        horizons = 1:size(var_fc_nc, 2)

                        # Create DataFrames
                        df_nc = DataFrame(horizon=horizons)
                        df_t1 = DataFrame(horizon=horizons)
                        df_t2 = DataFrame(horizon=horizons)

                        for (i, label) in enumerate(quantile_labels)
                            df_nc[!, label] = quantiles_nc[:, i]
                            df_t1[!, label] = quantiles_t1[:, i]
                            df_t2[!, label] = quantiles_t2[:, i]
                        end

                        # Construct file paths
                        varname = string(table_vars[i_var])
                        outpath_nc = tablespath(mnc, "forecast") * "\\Finshocks\\forecast_" * varname * "_cond=semi_para=full_vint=" * vint_last_obs *"_nc.csv"
                        outpath_t1 = tablespath(mnc, "forecast") * "\\Finshocks\\forecast_" * varname * "_cond=semi_para=full_vint=" * vint_last_obs *"_t1.csv"
                        outpath_t2 = tablespath(mnc, "forecast") * "\\Finshocks\\forecast_" * varname * "_cond=semi_para=full_vint=" * vint_last_obs *"_t2.csv"

                        # Save to CSV
                        CSV.write(outpath_nc, df_nc)
                        CSV.write(outpath_t1, df_t1)
                        CSV.write(outpath_t2, df_t2)
                    end

                # Possibly cross check with these commands...
                # scenario_means_bands(mnc, nowcast, verbose = :high,density_bands = [0.05, 0.1, 0.3, 0.5])
                # mb = read_scenario_mb(mnc, nowcast, :forecastobs)
                # mb.bands[:obs_gdp]

                
            end
           
end
  

