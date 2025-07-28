using DSGE, ModelConstructors, Distributed, Random
using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics
using Plots.PlotMeasures
#############################
# Written by: Matyas Farkas, IMF
# Based on the run_dsgevar.jl script of DSGE.jl
# 01 June 2025
#############################

# What do you want to do?
estimate_dsgevar  = true   # Estimate a DSGEVAR using SMC
use_estim_output  = true   # Use posterior from SMC estimation for following code. Otherwise, use default parameters.
get_VAR_system    = true   # Compute VAR coefficients and innovations-covariance matrix
do_modal_irf      = true   # Compute IRFs using modal parameters
compare_modal_irf = true   # Plot the modal DSGEVAR λ = ∞ rotation IRF and the actual DSGE IRF for observables
do_full_band_irf  = true   # Compute IRFs using parameters drawn from a distribution
create_meansbands = false  # Save full_band_irfs to MeansBands
do_parallel       = false  # Use parallel workers
n_workers         = 10
dsgevar_λ         = 0.5    # What λ do you want to use?

Random.seed!(20201793)

##############
# Model Setup
##############
# Instantiate the FRBNY DSGE model object
m = Model1010("ss20")

# estimate as of 2015-Q3 using the default data vintage from 2015 Nov 27
m <= Setting(:data_vintage, "250331")
m <= Setting(:date_forecast_start, quartertodate("2025-Q1"))
m <= Setting(:date_mainsample_start, DSGE.quartertodate("1960-Q1")) # Set bounds on the periods of
m <= Setting(:date_presample_start, DSGE.quartertodate("1959-Q3"))  # data being loaded so we don't

df = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
   mode_file = rawpath(m, "estimate", "paramsmode.h5")
      
# Switch off measurement errors 
m[:e_AAA_t]     =   0.
m[:e_BBB_t]     =   0.
m[:e_corepce_t] =   0.
m[:e_gdi_t]     =   0.
m[:e_gdi_t1]    =   0.
m[:e_gdp_t]     =   0.
m[:e_gdp_t1]    =   0.
m[:e_gdpdef_t]  =   0.

modal_paras = map(x -> x.value, m.parameters)
forecast_string = "test_dsgevar" # Change this to an empty string if you don't want an identifier for saved output
m <= Setting(:sampling_method, :SMC)
m <= Setting(:impulse_response_horizons, 20)

m <= Setting(:n_particles, 20)          # a small number just to make the estimation and IRFs faster.
m <= Setting(:use_fixed_schedule, true) # Fix schedule to fewer iterations
m <= Setting(:n_Φ, 10)                  # to speed up the estimation further

m <= Setting(:date_mainsample_start, DSGE.quartertodate("1960-Q1")) # Set bounds on the periods of
m <= Setting(:date_presample_start, DSGE.quartertodate("1959-Q3"))  # data being loaded so we don't
m <= Setting(:date_forecast_start, DSGE.quartertodate("2025-Q1"))   # get NaNs.
data = df_to_matrix(m, load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none))    # the rows of the data MUST correspond
                                        # to the exact order of the observables specified by the DSGEVAR.

# Construct a DSGEVAR from AnSchorfheide
dsgevar = DSGEVAR(m)
DSGE.update!(dsgevar, observables = collect(keys(m.observables)), shocks = collect(keys(m.exogenous_shocks)), lags = 4, λ = dsgevar_λ)

# Estimate the DSGEVAR
if estimate_dsgevar
    if do_parallel
        my_procs = addprocs(n_workers)
        @everywhere using DSGE, OrderedCollections
    end

    estimate(dsgevar, data; run_csminwel = false, verbose = :none)

    if do_parallel
        rmprocs(my_procs)
    end
end

