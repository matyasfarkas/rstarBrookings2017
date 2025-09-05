using DSGE, ClusterManagers, HDF5, Plots, StatsPlots
using DataFrames, CSV

##########################################################################################
## SETUP
##########################################################################################
##############
# Load the EA model with FG 
##############
use_FG_in_EA  = true  # set to false to load the EA model without FG shocks


# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")


m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
if use_FG_in_EA 
    m <= DSGE.Setting(:data_vintage, "250113")
else
    m <= DSGE.Setting(:data_vintage, "250114")
end
m <= DSGE.Setting(:reoptimize, false)
m <= DSGE.Setting(:calculate_hessian, false)
m <= DSGE.Setting(:date_mainsample_start,  quartertodate("1970-Q3"))
m <= DSGE.Setting(:date_presample_start,  quartertodate("1970-Q2"))

# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q3"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))


df = load_data(m; check_empty_columns = false)
if use_FG_in_EA
mode_file = joinpath("dsge/output_data/m1010/ss20/estimate/raw/" ,  "paramsmode_vint=250113.h5")
else
mode_file = joinpath("dsge/output_data/m1010/ss20/estimate/raw/" ,  "paramsmode_vint=250114.h5")
end
specify_mode!(m, mode_file)
system = DSGE.compute_system(m)

states = Dict{Symbol, Matrix{Float64}}()
shocks = Dict{Symbol, Matrix{Float64}}()
pseudo = Dict{Symbol, Matrix{Float64}}()

shock_labels = [key for (key, _) in sort(collect(m.exogenous_shocks), by = x -> x[2])]


combined = DSGE.OrderedDict{Symbol, Int64}()
# First insert all entries from endogenous_states
for (k, v) in m.endogenous_states
combined[k] = v
end

# Then insert entries from endogenous_states_augmented
for (k, v) in m.endogenous_states_augmented
combined[k] = v
end
state_labels = [key for (key, _) in sort(collect(combined), by = x -> x[2])]
pseudo_labels = [key for (key, _) in sort(collect(m.pseudo_observables), by = x -> x[2])]

system = DSGE.compute_system(m)

states_df = Dict{Symbol, DataFrame}()
shocks_df = Dict{Symbol, DataFrame}()
pseudo_df = Dict{Symbol, DataFrame}()

smoother = :durbin_koopman #:hamilton, :koopman, :carter_kohn, 
m <= DSGE.Setting(:forecast_smoother, smoother)

states[smoother], shocks[smoother], pseudo[smoother] = DSGE.smooth(m, df, system; draw_states = false)

dates = df.date[end-size(states[smoother],2)+1:end]
mat = states[smoother]'  # transpose to 259×91
states_df[smoother] = DataFrame(hcat(dates, mat), [:date; state_labels])

mat = shocks[smoother]'  # transpose to 259×29
shocks_df[smoother] = DataFrame(hcat(dates, mat), [:date; shock_labels])

mat = pseudo[smoother]'  # transpose to 259×22
pseudo_df[smoother] = DataFrame(hcat(dates, mat), [:date; pseudo_labels])

# Collect from the dataframe the respective shocks

#shocks_df[smoother][:, [:date; :b_liqtil_sh; :b_liqp_sh; :b_safetil_sh; :b_safep_sh]]
privilege_shock_names = [:b_liqtil_sh; :b_liqp_sh; :b_safetil_sh; :b_safep_sh]
privilege_shock_vals  = convert(Matrix,shocks_df[smoother][:, [:b_liqtil_sh; :b_liqp_sh; :b_safetil_sh; :b_safep_sh]])
 
horizon= size(privilege_shock_vals,1)
nshocks = size(system[:RRR], 2)
nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)
shocks = zeros(nshocks, horizon)

m1 = Model1010("ss20")
mode_file = joinpath("dsge/output_data/m1010/ss20/estimate/raw" ,  "paramsmode_vint=161223.h5")
specify_mode!(m1, mode_file)
system_US = DSGE.compute_system(m1)

for (i, shock_name) in enumerate(privilege_shock_names)
    shock_ind_US = m1.exogenous_shocks[shock_name]
    shocks[shock_ind_US, :] .= privilege_shock_vals[:, i]
end

# Compute IRF for a unit shock at time t for the specified shock
        states_privilege, obs_privilege, pseudo_privilege = forecast(system, s_0, shocks)
               
plotstart = 114 # 1999-Q1
horizon = size(privilege_shock_vals,1) - plotstart + 1
using Plots
p1 = plot(dates[plotstart:end],obs_privilege[m.observables[:obs_nominalrate],plotstart:end],title="Policy rate")
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(dates[plotstart:end],obs_privilege[m.observables[:obs_gdpdeflator],plotstart:end],title="Inflation")
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(dates[plotstart:end],states_privilege[m.endogenous_states[:y_t],plotstart:end],title="Output")#
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(dates[plotstart:end],pseudo_privilege[m.pseudo_observables[:Forward5YearRealNaturalRate],plotstart:end],title="r* (Forward 5-year real natural rate)")#
plot!(dates[plotstart:end],zeros(horizon,1),lc=:black,lw=2,label="")
# p6 = plot(1:horizon,pseudo[m.pseudo_observables[:RealNaturalRate],:, m.exogenous_shocks[:rm_sh]],title="Real natural rate")#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,layout=(2,2), legend=false)
plot!(size=(960,540))
if use_FG_in_EA 
    savefig( "Main results/Exorbitant privilege.pdf")   # saves the plot from p as a .pdf vector graphic

else
    savefig( "Main results/Exorbitant privilege without FG shocks in EA.pdf")   # saves the plot from p as a .pdf vector graphic
end


system10 = compute_system(m)
horzion = 40
states_irf10, obs_irf10, pseudo_irf10 = impulse_responses(system10, horizon)

# output
p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]],title="Permanent liquidity shock", label=["Basline model"])
plot!(legend=:bottomright)

p2 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_safep_sh]]] ,title="Permanent safety shock", label=["Basline model"])
plot!(legend=:bottomright)

p3 = plot(1:horizon,[ states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:zp_sh ]] ] ,title="Permanent technology shock", label=["Basline model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:rm_shl6 ]] ] ,title="FG6 shock", label=["Basline model" ])

plot(p1, p2, p3, p4, layout=(2,2), legend=false)


# 
p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]],title="Output", label=["Basline model"])
plot!(legend=:bottomright)

p2 = plot(1:horizon,[obs_irf10[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:b_liqp_sh]]] ,title="Inflation", label=["Basline model"])
plot!(legend=:bottomright)

p3 = plot(1:horizon,[ obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]] ] ,title="Policy rate", label=["Basline model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh ]] ] ,title="R* (Forward 5 year real natural rate)", label=["Basline model" ])

plot(p1, p2, p3, p4, layout=(2,2), legend=false)

