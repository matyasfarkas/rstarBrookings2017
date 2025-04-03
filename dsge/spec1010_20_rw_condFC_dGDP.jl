using DSGE, ModelConstructors, Distributed
using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl
# 31 March 2025
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
do_not_run_estimation = true
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
# addprocs(8)
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

finvars = [:date, :obs_nominalrate, :obs_BBBspread, :obs_AAAspread]

for i_year in 2020:2024
    for i_quarter in 1:4
        #  Construct model used for forecasting based on full sample estimated DSGE
                mnc = m
                mnc <= Setting(:data_vintage, "250331")
                mnc <=  Setting(:cond_semi_names, finvars,"Observables used in semiconditional forecasts")
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
                
                date_last_obs = quartertodate(string(i_year)*"-Q"*string(i_quarter))
                vint_last_obs = string(date_last_obs)
                vint_last_obs = vint_last_obs[3:4]*vint_last_obs[6:7]*vint_last_obs[9:10]

        #  Create population forecasts
                filename = inpath(m, "raw", "population_data_levels_" * "2024q4_full" * ".csv") # Load full population data 
                popdf = CSV.read(filename, DataFrame, copycols = true)
                DSGE.format_dates!(:date, popdf)
                sort!(popdf, :date)
                rename!(popdf,:CNP16OV => :POPULATION)
                popdf_forecast_to_write=filter!(row -> row.date < date_nc_cond, popdf)
                file_path = joinpath(joinpath(get_setting(mnc, :dataroot), "raw"), "population_forecast_" * vint_nc_cond * ".csv")
                CSV.write(file_path, popdf_forecast_to_write; missingstring="")
        
       
       #  Create nowcasts
                mnc <= Setting(:date_forecast_start,date_nc_cond)
                mnc <= Setting(:date_conditional_end,date_nc_cond)

                dfnc = load_data(mnc, try_disk = false, check_empty_columns = false, summary_statistics = :none)
                conditional_data_path = joinpath(get_setting(mnc, :dataroot),"data")
                conditional_path =  joinpath(get_setting(mnc, :dataroot), "cond")
                dfnc[end,Not(finvars)] .= missing
                df_to_write = DataFrame(dfnc[end,:])
                file_path = joinpath(conditional_path, "cond_cdid=02_cdvt="* vint_nc_cond *".csv")
                CSV.write(file_path, df_to_write; missingstring="")
        # #  Create population forecasts
        #         filename = inpath(m, "raw", "population_data_levels_" * "2024q4_full" * ".csv") # Load full population data 
        #         popdf = CSV.read(filename, DataFrame, copycols = true)
        #         DSGE.format_dates!(:date, popdf)
        #         sort!(popdf, :date)
        #         rename!(popdf,:CNP16OV => :POPULATION)
        #         popdf_forecast_to_write=filter!(row -> row.date < date_nc_cond, popdf)
        #         file_path = joinpath(joinpath(get_setting(mnc, :dataroot), "raw"), "population_forecast_" * vint_nc_cond * ".csv")
        #         CSV.write(file_path, popdf_forecast_to_write; missingstring="")

        # Create dataset to append the conditioning - Depreciated as we can pass the dataframe directly.
                # mnc <= Setting(:cond_vintage, vint_nc_cond)
                # mnc <= Setting(:date_mainsample_end,date_last_obs)
                # mnc <= Setting(:date_forecast_start,date_last_obs)
                # dfrt = load_data(mnc, try_disk = false, check_empty_columns = false, summary_statistics = :none)
                # This creates the file boguosly nowcasts in data_cdid=02_cdvt=vint_nc_cond_dsid=04_vint=250331.csv has the wrong data at the end...
                # forecast_one(mnc, :mode, :semi, output_vars; check_empty_columns = false, verbose = :high) 
                # Need to manually update the conditioning from the conditioning file... 
                # file_path = joinpath(conditional_data_path, "data_cdid=02_cdvt="* vint_nc_cond *"_dsid=04_vint=250331.csv")
                # CSV.write(file_path, dfnc; missingstring="")
        
        # Compute and store forecasts
                mnc <= Setting(:forecast_block_size, 100) # adjust block size to run on smaller number of posterior draws
                mnc <= Setting(:date_forecast_start,date_last_obs)
                mnc <= Setting(:date_conditional_end,date_last_obs)
                mnc <= Setting(:date_mainsample_end,date_last_obs)
                forecast_one(mnc, :full, :semi, output_vars; df= dfnc, check_empty_columns = false, verbose = :high)
                compute_meansbands(mnc, :full, :semi, output_vars; df= dfnc, check_empty_columns = false)
                plot_history_and_forecast(mnc, :obs_gdp , :obs, :full, :semi,
                                            use_bdd = :bdd_and_unbdd,
                                            start_date = DSGE.quartertodate(string(i_year-10)*"-Q"*string(1)),
                                            end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
                                            verbose = :none)
                mv(figurespath(mnc,"forecast" ) *"\\forecast_obs_gdp_cond=semi_para=full_vint=250331.pdf", figurespath(mnc,"forecast" ) *"\\obs_gdp_forecast_density_semi_nowcast_cond="* vint_last_obs *".pdf",force=true)

                table_vars = [:obs_gdp;  :obs_hours;       :obs_wages        ;       :obs_gdpdeflator  ;       :obs_corepce      ;       :obs_nominalrate  ;       :obs_consumption  ;       :obs_investment   ;       :obs_BBBspread    ;       :obs_longinflation;        :obs_longrate     ;        :obs_tfp          ;        :obs_gdi          ;        :obs_AAAspread    ;        :obs_nominalrate1 ;        :obs_nominalrate2 ;        :obs_nominalrate3 ;        :obs_nominalrate4 ;        :obs_nominalrate5 ;        :obs_nominalrate6]
                write_meansbands_tables_all(mnc, :full, :semi,output_vars, forecast_string ="" , vars = table_vars)
                for i_var in size(table_vars,1)
                    mv(tablespath(mnc,"forecast" ) *"\\forecast_" * string(table_vars[i_var]) * "_cond=semi_para=full_vint=250331.csv", tablespath(mnc,"forecast") * "\\" * string(table_vars[i_var]) * "_forecast_density_semi_nowcast_cond="* vint_last_obs  *".csv",force=true)
                end            

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


        # Generate 1Y conditioning Dataframe
                mnc <= Setting(:date_forecast_start,date_t1_cond)
                mnc <= Setting(:date_conditional_end,date_t1_cond)
                mnc <= Setting(:cond_cdid,"03")

                dft1 = load_data(mnc, try_disk = false, check_empty_columns = false, summary_statistics = :none)
                dft1[end-4:end,Not(finvars)] .= missing
                df_to_write = DataFrame(dft1[end-4:end,:])
                file_path = joinpath(conditional_path, "cond_cdid=03_cdvt="* vint_nc_cond *".csv")
                CSV.write(file_path, df_to_write; missingstring="")                
        # #   Create population forecasts
        #         filename = inpath(m, "raw", "population_data_levels_" * "2024q4_full" * ".csv") # Load full population data 
        #         popdf = CSV.read(filename, DataFrame, copycols = true)
        #         DSGE.format_dates!(:date, popdf)
        #         sort!(popdf, :date)
        #         rename!(popdf,:CNP16OV => :POPULATION)
        #         popdf_forecast_to_write=filter!(row -> row.date < date_t1_cond, popdf)
        #         file_path = joinpath(joinpath(get_setting(mnc, :dataroot), "raw"), "population_forecast_" * vint_t1_cond * ".csv")
        #         CSV.write(file_path, popdf_forecast_to_write; missingstring="")

        # Compute and store forecasts
                mnc <= Setting(:forecast_block_size, 100) # adjust block size to run on smaller number of posterior draws
                mnc <= Setting(:date_forecast_start,dft1[end-5,:date])
                mnc <= Setting(:date_mainsample_end,dft1[end-5,:date])
                mnc <= Setting(:date_conditional_end,date_t1_cond_end)
                forecast_one(mnc, :full, :semi, output_vars; df= dft1, check_empty_columns = false, verbose = :high)
                compute_meansbands(mnc, :full, :semi, output_vars; df= dft1, check_empty_columns = false)
                plot_history_and_forecast(mnc, :obs_gdp , :obs, :full, :semi,
                                            use_bdd = :bdd_and_unbdd,
                                            start_date = DSGE.quartertodate(string(i_year-10)*"-Q"*string(1)),
                                            end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
                                            verbose = :none)
                mv(figurespath(mnc,"forecast" ) *"\\forecast_obs_gdp_cond=semi_para=full_vint=250331.pdf", figurespath(mnc,"forecast" ) *"\\obs_gdp_forecast_density_semi_1year_h5_cond="* vint_last_obs *".pdf",force=true)

                table_vars = [:obs_gdp;  :obs_hours;       :obs_wages        ;       :obs_gdpdeflator  ;       :obs_corepce      ;       :obs_nominalrate  ;       :obs_consumption  ;       :obs_investment   ;       :obs_BBBspread    ;       :obs_longinflation;        :obs_longrate     ;        :obs_tfp          ;        :obs_gdi          ;        :obs_AAAspread    ;        :obs_nominalrate1 ;        :obs_nominalrate2 ;        :obs_nominalrate3 ;        :obs_nominalrate4 ;        :obs_nominalrate5 ;        :obs_nominalrate6]
                write_meansbands_tables_all(mnc, :full, :semi,output_vars, forecast_string ="" , vars = table_vars)
                for i_var in size(table_vars,1)
                    mv(tablespath(mnc,"forecast" ) *"\\forecast_" * string(table_vars[i_var]) * "_cond=semi_para=full_vint=250331.csv", tablespath(mnc,"forecast") * "\\" * string(table_vars[i_var]) * "_forecast_density_semi_1year_cond="* vint_last_obs * ".csv",force=true) 
                end    
                    
                
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

        # Generate 2Y conditioning Dataframe
                mnc <= Setting(:date_forecast_start,date_t2_cond)
                mnc <= Setting(:date_conditional_end,date_t2_cond)
                mnc <= Setting(:cond_cdid,"04")

                dft2 = load_data(mnc, try_disk = false, check_empty_columns = false, summary_statistics = :none)
                dft2[end-8:end,Not(finvars)] .= missing
                df_to_write = DataFrame(dft1[end-8:end,:])
                file_path = joinpath(conditional_path, "cond_cdid=04_cdvt="* vint_nc_cond *".csv")
                CSV.write(file_path, df_to_write; missingstring="")                
        # #   Create population forecasts
        #         filename = inpath(m, "raw", "population_data_levels_" * "2024q4_full" * ".csv") # Load full population data 
        #         popdf = CSV.read(filename, DataFrame, copycols = true)
        #         DSGE.format_dates!(:date, popdf)
        #         sort!(popdf, :date)
        #         rename!(popdf,:CNP16OV => :POPULATION)
        #         popdf_forecast_to_write=filter!(row -> row.date < date_t1_cond, popdf)
        #         file_path = joinpath(joinpath(get_setting(mnc, :dataroot), "raw"), "population_forecast_" * vint_t1_cond * ".csv")
        #         CSV.write(file_path, popdf_forecast_to_write; missingstring="")

        # Compute and store forecasts
                mnc <= Setting(:forecast_block_size, 100) # adjust block size to run on smaller number of posterior draws
                mnc <= Setting(:date_forecast_start,dft2[end-9,:date])
                mnc <= Setting(:date_mainsample_end,dft2[end-9,:date])
                mnc <= Setting(:date_conditional_end,dft2[end,:date])
                forecast_one(mnc, :full, :semi, output_vars; df= dft2, check_empty_columns = false, verbose = :high)
                compute_meansbands(mnc, :full, :semi, output_vars; df= dft1, check_empty_columns = false)
                plot_history_and_forecast(mnc, :obs_gdp , :obs, :full, :semi,
                                            use_bdd = :bdd_and_unbdd,
                                            start_date = DSGE.quartertodate(string(i_year-10)*"-Q"*string(1)),
                                            end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
                                            verbose = :none)
                mv(figurespath(mnc,"forecast" ) *"\\forecast_obs_gdp_cond=semi_para=full_vint=250331.pdf", figurespath(mnc,"forecast" ) *"\\obs_gdp_forecast_density_semi_2years_h9_cond="* vint_last_obs *".pdf",force=true)

                table_vars = [:obs_gdp;  :obs_hours;       :obs_wages        ;       :obs_gdpdeflator  ;       :obs_corepce      ;       :obs_nominalrate  ;       :obs_consumption  ;       :obs_investment   ;       :obs_BBBspread    ;       :obs_longinflation;        :obs_longrate     ;        :obs_tfp          ;        :obs_gdi          ;        :obs_AAAspread    ;        :obs_nominalrate1 ;        :obs_nominalrate2 ;        :obs_nominalrate3 ;        :obs_nominalrate4 ;        :obs_nominalrate5 ;        :obs_nominalrate6]
                write_meansbands_tables_all(mnc, :full, :semi,output_vars, forecast_string ="" , vars = table_vars)
                for i_var in 1:size(table_vars,1)
                    mv(tablespath(mnc,"forecast" ) *"\\forecast_" * string(table_vars[i_var]) * "_cond=semi_para=full_vint=250331.csv", tablespath(mnc,"forecast") * "\\" * string(table_vars[i_var]) * "_forecast_density_semi_2years_cond="* vint_last_obs *".csv",force=true)
                end            
            end
           
end
  

