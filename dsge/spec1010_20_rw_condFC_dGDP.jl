using DSGE, ModelConstructors, Distributed
using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl
# 28 March 2025
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
if run_estimation
   
    # Start from full sample mode
   mode_file = rawpath(m, "estimate", "paramsmode.h5")
   #mode_file = replace(mode_file, "ss20", "ss18")
   DSGE.update!(m, h5read(mode_file, "params"))
else
    df = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
    data = df_to_matrix(m, df)
    DSGE.estimate(m, data)

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
m <= Setting(:cond_vintage, "250331")
forecast_one(m, :mode, :full, output_vars; check_empty_columns = false)
compute_meansbands(m, :mode, :full, output_vars; check_empty_columns = false)

table_vars = [:obs_gdp]
write_meansbands_tables_all(m, :mode, :full, [:forecastobs], forecast_string = "", vars = table_vars)


# Full-distribution forecast (point forecast (mean) and uncertainty bands)
# Optionally add 10 processes to run the forecast in parallel (uncomment the 3 lines below).
# Alternatively, you can load the ClusterManagers package and add processes
# using one of the schedulers such as SGE or Slurm.
# addprocs(10)
# @everywhere using DSGE
# m <= Setting(:use_parallel_workers, true)

# Comment out the line below if you did not run the forecast in parallel.
# rmprocs(procs())
plot_history_and_forecast(m, :obs_gdp, :obs, :mode, :full,
                              use_bdd = :bdd_and_unbdd,
                              start_date = DSGE.quartertodate("2007-Q1"),
                              end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
                              verbose = :none)

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

for i_year in 2000:2001
    for i_quarter in 1:4
        #  Construct conditioning datasets
        mnc = m

        
        #  Create Nowcast cond series
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
       mnc <= Setting(:date_forecast_start,date_nc_cond)
       mnc <= Setting(:date_conditional_end,date_nc_cond)

       dfnc = load_data(mnc, try_disk = false, check_empty_columns = false, summary_statistics = :none)

        conditional_path =  joinpath(get_setting(mnc, :dataroot), "cond")
        dfnc[end,Not(finvars)] .= missing
        df_to_write = DataFrame(dfnc[end,:])
        file_path = joinpath(conditional_path, "cond_cdid=02_cdvt=" * vint_nc_cond * ".csv")
        CSV.write(file_path, df_to_write; missingstring="")


        date_last_obs = quartertodate(string(i_year)*"-Q"*string(i_quarter))
        mnc <= Setting(:date_forecast_start,date_last_obs)

        dfnc = load_data(mnc, try_disk = false, check_empty_columns = false, summary_statistics = :none)


        #  Create population forecasts
        filename = inpath(m, "raw", "population_data_levels_" * "2024q4_full" * ".csv") # Load full population data 
        popdf = CSV.read(filename, DataFrame, copycols = true)
        DSGE.format_dates!(:date, popdf)
        sort!(popdf, :date)
        rename!(popdf,:CNP16OV => :POPULATION)
        popdf_forecast_to_write=filter!(row -> row.date < date_nc_cond, popdf)
        file_path = joinpath(joinpath(get_setting(mnc, :dataroot), "raw"), "population_forecast_" * vint_nc_cond * ".csv")
        CSV.write(file_path, popdf_forecast_to_write; missingstring="")


        # Update the forecast start to be the looper in the model
        mnc <= Setting(:forecast_block_size, 100) # adjust block size to run on small number of estimations
        forecast_one(mnc, :full, :full, output_vars; check_empty_columns = false, verbose = :low)
        compute_meansbands(mnc, :full, :full, output_vars; check_empty_columns = false)
        plot_history_and_forecast(mnc, :obs_gdp, :obs, :mode, :full,
                                      use_bdd = :bdd_and_unbdd,
                                      start_date = DSGE.quartertodate(string(i_year-10)*"-Q"*string(1)),
                                      end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
                                      verbose = :none)
        table_vars = [:obs_gdp]
        write_meansbands_tables_all(mnc, :mode, :full, [:forecastobs], forecast_string ="" , vars = table_vars)
        mv(tablespath(mnc,"forecast")*"\\forecast_obs_gdp_cond=full_para=mode_vint=250331.csv", tablespath(mnc,"forecast")*"\\forecast_obs_gdp_nowcast_cond="* string(date_last_obs) *".csv",force=true)
        # #
        # vint_1yearahead_cond= string(quartertodate(string(i_year+1)*"-Q"*string(i_quarter)))
        # vint_1yearahead_cond = vint_1yearahead_cond[3:4]*vint_1yearahead_cond[6:7]*vint_1yearahead_cond[9:10]

        # vint_2yearsahead_cond= string(quartertodate(string(i_year+2)*"-Q"*string(i_quarter)))
        # vint_2yearsahead_cond = vint_2yearsahead_cond[3:4]*vint_2yearsahead_cond[6:7]*vint_2yearsahead_cond[9:10]

    
        # m <= DSGE.Setting(:data_vintage, vint)

        # forecast_one(m, :full, :full, output_vars; check_empty_columns = false)
        # compute_meansbands(m, :full, :none, output_vars; check_empty_columns = false)
        # plot_history_and_forecast(m, :obs_gdp, :obs, :full, :none,
        #                             use_bdd = :bdd_and_unbdd,
        #                             start_date = DSGE.quartertodate(string(i_year-5)*"-Q"*string(i_quarter)),
        #                             end_date = DSGE.iterate_quarters(date_forecast_start(m), 10),
        #                             verbose = :none)

    end
           
end
  

