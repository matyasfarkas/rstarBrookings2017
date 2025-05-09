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
# do_not_run_estimation = false
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
# # Load posterior draws

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

drawnum = 100;
simnum = 100
n_periods=size(df, 1)
cor_sim = Matrix{Float64}(undef, simnum,drawnum)

robs_sim = Array{Float64}(undef, n_periods,simnum,drawnum);
rstar_sim = Array{Float64}(undef, n_periods,simnum,drawnum);
ϵ_sim = Array{Float64}(undef,29,262,simnum,drawnum);
drawi=1

for drawloop in 1:40:4000
    params = params_draws[drawloop,:]

    DSGE.update!(m, params)
    system = compute_system(m, tvis = false)
    # Unpack system
    TTT    = system[:TTT]
    RRR    = system[:RRR]
    CCC    = system[:CCC]
    QQ     = system[:QQ]
    

    for simi = 1:simnum
                # sim_data = DSGE.simulate_observables(m; burnin=500, n_periods=size(df, 1))'
                ϵ_draw = rand(DegenerateMvNormal(zeros(29), sqrt.(QQ)), 500 + size(df, 1))

                s, ϵ = DSGE.simulate_states(TTT,RRR,CCC,QQ; burnin =500,n_periods = 262,ϵ=ϵ_draw)
                pseudo = DSGE.pseudo_measurement(m, TTT, RRR,CCC) 
                pseudo_y = Matrix{Float64}(undef, size(pseudo.ZZ_pseudo,1), n_periods)

                for i in 1:n_periods
                    pseudo_y[:, i] = pseudo.DD_pseudo + pseudo.ZZ_pseudo*s[:, i + 500] 
                end

                y = Matrix{Float64}(undef, length(system[:DD]), n_periods)

                for i in 1:n_periods
                    y[:, i] = system[:DD] + system[:ZZ]*s[:, i + 500] 
                end

                rstar = pseudo_y[m.pseudo_observables[:Forward5YearRealNaturalRate], :];
                robs = y[m.observables[:obs_nominalrate],:];
                cor_sim[simi,drawi]=cor(robs,rstar);
                robs_sim[:,simi,drawi]=robs;
                rstar_sim[:,simi,drawi]=rstar;
                ϵ_sim[:,:,simi,drawi] = ϵ[:,501:end];

    end
    println("Posterior draw " * string(drawi) * " done...")
    drawi= drawi + 1;

end

using Plots
plot(1:262,[rstar_sim[:,100,100] ,  robs_sim[:,100,100]],title="NY-FED simulations", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate"])
savefig( "obsnominalrate_vs_5yforwardrealrate.pdf")   # saves the plot from p as a .pdf vector graphic

CSV.write("update/simulated_correlation.csv",DataFrame(reshape(cor_sim,10000)',:auto))
CSV.write("update/simulated_correlation_description.csv",describe(DataFrame((cor_sim),:auto)))
# Now the shock decomposition
        drawi=1
        robs_HVD = Array{Float64}(undef,262,29,simnum,drawnum);
        rstar_HVD = Array{Float64}(undef,262,29,simnum,drawnum);

        for simi = 1:simnum
            drawi=1
            for drawloop in 1:40:4000
            params = params_draws[drawloop,:]
            DSGE.update!(m, params)
            system = compute_system(m, tvis = false)
            stateHVD, obsHVD , pseudoHVD = shock_decompositions(system, 1, ϵ_sim[:,:,simi,drawi], 1, 262)
            robs_HVD[:,:,simi,drawi] =obsHVD[ m.observables[:obs_nominalrate],:,:]
            rstar_HVD[:,:,simi,drawi] =pseudoHVD[ m.pseudo_observables[:Forward5YearRealNaturalRate],:,:]
            # println("Posterior draw " * string(drawi) * " done...")
            drawi = drawi +1
            end
            println("HVD simulation draw " * string(simi) * " done...")

        end

robs_HVD_towrite = reshape(robs_HVD, 262, 29, 100 * 100)
rstar_HVD_towrite = reshape(rstar_HVD, 262, 29, 100 * 100)

hvdcorr = Array{Float64}(undef,10000, 29);

for simi = 1:10000
    for shocki = 1:29
    hvdcorr[simi,shocki] = cor(robs_HVD_towrite[:,shocki,simi],rstar_HVD_towrite[:,shocki,simi])
    end
end

mean_robs_HVD_df_towrite = DataFrame(mean(robs_HVD_towrite,dims=3)[:,:,1],collect(keys(m.exogenous_shocks)))
mean_rstar_HVD_df_towrite = DataFrame(mean(rstar_HVD_towrite,dims=3)[:,:,1],collect(keys(m.exogenous_shocks)))

CSV.write("update/simulated_HVD_robs.csv",mean_robs_HVD_df_towrite)
CSV.write("update/simulated_HVD_rstar.csv",mean_rstar_HVD_df_towrite)

corr_summary = describe(DataFrame(collect(keys(m.exogenous_shocks))[1] => vec(hvdcorr[:,1])))

for shocki =1:29
    try
    CSV.write("update/simulated_HVD_corr_"* string(collect(keys(m.exogenous_shocks))[shocki]) *".csv",describe(DataFrame(collect(keys(m.exogenous_shocks))[shocki] => vec(hvdcorr[:,shocki]))))
    newrow = describe(DataFrame(collect(keys(m.exogenous_shocks))[shocki] => vec(hvdcorr[:,shocki])))
    corr_summary = vcat(corr_summary, newrow)

    catch

        println("Structural shock " * string(collect(keys(m.exogenous_shocks))[shocki]) * " skipped as it had been switched off...")
    end
end


CSV.write("update/simulated_HVD_corr_allshocks.csv",(corr_summary))
