
using DSGE, Dates, DataFrames,OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics, ModelConstructors, LinearAlgebra
# using Nullables, DataFrames, OrderedCollections, Dates,HDF5, CSV, JLD2, FileIO, Statistics
using Measures  # for mm margins

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl and the test forecast driver file
# 25 August 2025
#############################

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")
# params_mode = load_draws(m, :mode)
# # Switch off the pi_target shock
# params_mode[m.param_index[:σ_π_target_sh]] = 0.0
# DSGE.update!(m, params_mode)

# Settings for data, paths, etc.
m <= DSGE.Setting(:data_vintage, "161223")
# Settings for forecast dates
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))
 mode_file = rawpath(m, "estimate", "paramsmode.h5")
        #mode_file = replace(mode_file, "ss20", "ss18")
        DSGE.update!(m, h5read(mode_file, "params"))


system = DSGE.compute_system(m)



nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)

##############
# Function Setup
##############


function obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, m::AbstractDSGEModel, var_name::Symbol, shock_inds::Matrix{Int}, system::System{Float64})
    # shock_inds: matrix of size (nshocks, horizon), each column is the set of shock indices for that period
    horizon = length(x)
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    s_0 = zeros(nstates)
    shocks = zeros(nshocks, horizon)
    var_names, var_class, peg_ind =
        if var_name in keys(m.endogenous_states)
            m.endogenous_states, :states, m.endogenous_states[var_name]
        elseif var_name in keys(m.observables)
            m.observables, :obs,  m.observables[var_name]
        elseif var_name in keys(m.pseudo_observables)
            m.pseudo_observables, :pseudo,  m.pseudo_observables[var_name]
        else
            error("Variable $var_name not found in endogenous states, observables, or pseudo-observables.")
            return
        end

    for t in 1:horizon
        # Build IRF vector for all shocks at time t
        n_shocks_t = size(shock_inds, 1)
        IRFvec = zeros(n_shocks_t)
        for j in 1:n_shocks_t
            test_shocks = zeros(nshocks, horizon)
            test_shocks[shock_inds[j, t], t] = 1.0
            states, obs, pseudo = forecast(system, s_0, test_shocks)
            if var_class == :states
                IRFvec[j] = states[peg_ind, t]
            elseif var_class == :obs
                IRFvec[j] = obs[peg_ind, t]
            elseif var_class == :pseudo
                IRFvec[j] = pseudo[peg_ind, t]
            end
        end
        # Compute effect of previous shocks
        prev_effect = 0.0
        if t > 1
            prev_shocks = shocks[:, 1:t-1]
            prev_states, prev_obs, prev_pseudo = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
            if var_class == :states
                prev_effect = prev_states[peg_ind, t]
            elseif var_class == :obs
                prev_effect = prev_obs[peg_ind, t]
            elseif var_class == :pseudo
                prev_effect = prev_pseudo[peg_ind, t]
            end
        end
        # Solve for required shocks: IRFvec * shocks_t = x[t] - prev_effect
        # If underdetermined, use least squares
        shocks_t = pinv(IRFvec') * (x[t] - prev_effect)
        for j in 1:n_shocks_t
            shocks[shock_inds[j, t], t] = shocks_t[j]
        end
    end
    return shocks
end

# Example test (assuming you have a model m and system):
# x = [1.0, 2.0]
# shock_inds = [m.exogenous_shocks[:shock1] m.exogenous_shocks[:shock1];
#               m.exogenous_shocks[:shock2] m.exogenous_shocks[:shock2]]
# shocks = obtain_shocks_from_desired_state_path_iterative(x, m, :obs_nominalrate, shock_inds, system)


# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]

## Load in HLW real time estiamtes of R*
csv_path = joinpath(basepath, "Main results", "DSGE_vs_HLW.csv")
hlw_rstar = DataFrame(CSV.File(csv_path))

valid_idx = findall(row -> !ismissing(row[:date]) && !ismissing(row[:HLW]) && !ismissing(row[:mean]), eachrow(hlw_rstar))

rstar_diff = hlw_rstar.HLW[valid_idx] .- hlw_rstar.mean[valid_idx]
dates= hlw_rstar.date[valid_idx]

desired_path = rstar_diff#rstar_diff[end-16:end] # Desired path for the state variable
var_name =:Forward5YearRealNaturalRate

shock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh, :zp_sh ] # Convenience yield shocks causing the difference

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="HLW r* minus DSGE r*")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "Main results/rstar_what_if_HLW.pdf")   # saves the plot from p as a .pdf vector graphic





# Alternative if r* had been the HLW post COVID
desired_path = hlw_rstar.HLW[end-20:259] .- hlw_rstar.mean[end-20:259]

var_name =:Forward5YearRealNaturalRate

# hock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh ] # Convenience yield shocks causing the difference
# shock_syms = [  :zp_sh ] # Convenience yield shocks causing the difference
shock_syms = collect(Iterators.filter(k -> !startswith(String(k), "rm_sh"), keys(m.exogenous_shocks)))

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:Forward5YearRealNaturalRate,:obs_nominalrate,  :pi_t,  :y_t,:ExAnteRealRate,:RealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Post-Covid difference of HLW vs DSGE ")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.pseudo_observables[:π_t],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.pseudo_observables[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))
savefig("Main results/what_if_rstar_had_been_HLW_post_COVID.pdf")   # saves the plot from p as a .pdf vector graphic

# --- Write plotted series to CSV ---
using CSV, DataFrames
df_out_norstarincrease_WZ = DataFrame(
        Date = plotdates,
        HLW_minus_DSGE = desired_path,
        PolicyRate = obs[m.observables[:obs_nominalrate], :],
        Inflation = obs[m.observables[:obs_gdpdeflator], :],
        Output = states[m.endogenous_states[:y_t], :],
        ExAnteRealRate = pseudo[m.pseudo_observables[:ExAnteRealRate], :],
        RealNaturalRate = pseudo[m.pseudo_observables[:RealNaturalRate], :]
)
CSV.write("Main results/what_if_rstar_had_been_HLW_post_COVID.csv", df_out_norstarincrease_WZ)


# === Baseline vs Counterfactual Plotting Section ===
using CSV, DataFrames

# Directory and vintage for DSGE smoothed series
dsge_table_dir = joinpath(basepath, "dsge", "output_data", "m1010", "ss20", "forecast", "tables")
vintage = "250825"  # Update if needed to match your DSGE output
cond = "none"
para = "mode"

# Helper to get the correct filename for a variable
function get_hist_filename(var::Symbol)
        return "hist_" * String(var) * "_cond=" * cond * "_para=" * para * "_vint=" * vintage * ".csv"
end


# Load all baseline series into a Dict{Symbol, DataFrame} using CSV.File
baseline_dfs = Dict{Symbol, DataFrame}()
for v in plotvars
        fname = get_hist_filename(v)
        fpath = joinpath(dsge_table_dir, fname)
        if isfile(fpath)
                baseline_dfs[v] = DataFrame(CSV.File(fpath))
        else
                @warn "File not found for variable $(v): $(fpath)"
        end
end

# --- Plot all variables in a 3x2 grid: baseline (blue) vs counterfactual (red) ---
using Plots

titlelist = ["r* (Forward 5-year real natural rate)","Policy rate", "Inflation", "Output", "Ex-ante real rate", "Real natural rate (flex-price real rate)"]
zeroline = zeros(length(plotdates))
plots_arr = Vector{Any}(undef, length(plotvars))
for (i, v) in enumerate(plotvars)

        if !haskey(baseline_dfs, v)
                plots_arr[i] = plot(title=string(v), legend=false) # empty plot if missing
                continue
        end
        df_base = baseline_dfs[v]
        varcol = names(df_base)[names(df_base) .!= :date][end]
        # Filter to plotdates and ensure order matches plotdates
        df_base_short = DataFrames.filter(row -> row.date in plotdates, df_base)
        # If not already sorted, sort by date
        sort!(df_base_short, :date)
        # Get counterfactual series for this variable
        if v == :pi_t
                cf_series = pseudo[m.pseudo_observables[:π_t], :]
                else
                if v in keys(m.observables)
                        cf_series = obs[m.observables[v], :]
                elseif v in keys(m.endogenous_states)
                        cf_series = states[m.endogenous_states[v], :]
                elseif v in keys(m.pseudo_observables)
                        cf_series = pseudo[m.pseudo_observables[v], :]
                else
                        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
                        plots_arr[i] = plot(title=string(v), legend=false)
                        continue
                end
        end

        cf_short = cf_series[end-length(plotdates)+1:end]
        # Format x-ticks to show only the year
        # Format x-ticks to show only the year (yyyy)
        years = unique(year.(plotdates))
        # Find the first date in plotdates for each year
        year_tick_dates = [findfirst(d -> year(d) == y, plotdates) !== nothing ? plotdates[findfirst(d -> year(d) == y, plotdates)] : Date(string(y)*"-01-01") for y in years]
        year_tick_labels = [string(y) for y in years]
        xtick_tuple = (year_tick_dates, year_tick_labels)
        # Plot: baseline (blue), counterfactual (red), zero line (black)
        p = plot(df_base_short.date, df_base_short[!, varcol], label="DSGE Baseline", color=:blue, lw=2, title=titlelist[i], xticks=xtick_tuple)
        plot!(p, plotdates, cf_short+df_base_short[!, varcol], label="Counterfactual HLW r*", color=:red, lw=2)
        plot!(p, plotdates, zeroline, lc=:black, lw=2, label="")
        plot!(p, legend=false)
        plots_arr[i] = p
end

# Fill up to 6 plots if plotvars < 6
while length(plots_arr) < 6
        push!(plots_arr, plot(title="", legend=false))
end

plt = plot(plots_arr[1], plots_arr[2], plots_arr[3], plots_arr[4], plots_arr[5], plots_arr[6], layout=(3,2), legend=false)
plot!(plt, size=(960,540))
savefig(plt, "Main results/compare_baseline_vs_HLW_grid.pdf")


# Alternative if r* did not increase post COVID19
desired_path =hlw_rstar.mean[end-20:259] .- hlw_rstar.mean[end-20]  #rstar_diff[end-16:end] # Desired path for the state variable
# desired_path = -desired_path
var_name =:Forward5YearRealNaturalRate

# hock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh ] # Convenience yield shocks causing the difference
# shock_syms = [  :zp_sh ] # Convenience yield shocks causing the difference
shock_syms = collect(Iterators.filter(k -> !startswith(String(k), "rm_sh"), keys(m.exogenous_shocks)))

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:Forward5YearRealNaturalRate,:obs_nominalrate,  :pi_t,  :y_t,:ExAnteRealRate,:RealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Change in r* since end of COVID19")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

savefig( "Main results/rstar_had_not_increased_woFG.pdf")   # saves the plot from p as a .pdf vector graphic

# === Baseline vs Counterfactual Plotting Section ===
using CSV, DataFrames

# Directory and vintage for DSGE smoothed series
dsge_table_dir = joinpath(basepath, "dsge", "output_data", "m1010", "ss20", "forecast", "tables")
vintage = "250825"  # Update if needed to match your DSGE output
cond = "none"
para = "mode"

# Helper to get the correct filename for a variable
function get_hist_filename(var::Symbol)
        return "hist_" * String(var) * "_cond=" * cond * "_para=" * para * "_vint=" * vintage * ".csv"
end


# Load all baseline series into a Dict{Symbol, DataFrame} using CSV.File
baseline_dfs = Dict{Symbol, DataFrame}()
for v in plotvars
        fname = get_hist_filename(v)
        fpath = joinpath(dsge_table_dir, fname)
        if isfile(fpath)
                baseline_dfs[v] = DataFrame(CSV.File(fpath))
        else
                @warn "File not found for variable $(v): $(fpath)"
        end
end

# --- Plot all variables in a 3x2 grid: baseline (blue) vs counterfactual (red) ---
using Plots

titlelist = ["r* (Forward 5-year real natural rate)","Policy rate", "Inflation", "Output", "Ex-ante real rate", "Real natural rate (flex-price real rate)"]
zeroline = zeros(length(plotdates))
plots_arr = Vector{Any}(undef, length(plotvars))
for (i, v) in enumerate(plotvars)

        if !haskey(baseline_dfs, v)
                plots_arr[i] = plot(title=string(v), legend=false) # empty plot if missing
                continue
        end
        df_base = baseline_dfs[v]
        varcol = names(df_base)[names(df_base) .!= :date][end]
        # Filter to plotdates and ensure order matches plotdates
        df_base_short = DataFrames.filter(row -> row.date in plotdates, df_base)
        # If not already sorted, sort by date
        sort!(df_base_short, :date)
        # Get counterfactual series for this variable
        if v == :pi_t
                cf_series = pseudo[m.pseudo_observables[:π_t], :]
                else
                if v in keys(m.observables)
                        cf_series = obs[m.observables[v], :]
                elseif v in keys(m.endogenous_states)
                        cf_series = states[m.endogenous_states[v], :]
                elseif v in keys(m.pseudo_observables)
                        cf_series = pseudo[m.pseudo_observables[v], :]
                else
                        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
                        plots_arr[i] = plot(title=string(v), legend=false)
                        continue
                end
        end

        cf_short = cf_series[end-length(plotdates)+1:end]
        # Format x-ticks to show only the year
        # Format x-ticks to show only the year (yyyy)
        years = unique(year.(plotdates))
        # Find the first date in plotdates for each year
        year_tick_dates = [findfirst(d -> year(d) == y, plotdates) !== nothing ? plotdates[findfirst(d -> year(d) == y, plotdates)] : Date(string(y)*"-01-01") for y in years]
        year_tick_labels = [string(y) for y in years]
        xtick_tuple = (year_tick_dates, year_tick_labels)
        # Plot: baseline (blue), counterfactual (red), zero line (black)
        p = plot(df_base_short.date, df_base_short[!, varcol], label="DSGE Baseline", color=:blue, lw=2, title=titlelist[i], xticks=xtick_tuple)
        plot!(p, plotdates, -cf_short+df_base_short[!, varcol], label="Counterfactual HLW r*", color=:red, lw=2)
        plot!(p, plotdates, zeroline, lc=:black, lw=2, label="")
        plot!(p, legend=false)
        plots_arr[i] = p
end

# Fill up to 6 plots if plotvars < 6
while length(plots_arr) < 6
        push!(plots_arr, plot(title="", legend=false))
end

plt = plot(plots_arr[1], plots_arr[2], plots_arr[3], plots_arr[4], plots_arr[5], plots_arr[6], layout=(3,2), legend=false)
plot!(plt, size=(960,540))
savefig(plt, "Main results/compare_baseline_vs_no_increase_grid.pdf")

# Alternative if r* did not increase post COVID19
desired_path_wFG =hlw_rstar[end-20:259,10] .- hlw_rstar[end-20,10]  #rstar_diff[end-16:end] # Desired path for the state variable
shocks_path_wFG = obtain_shocks_from_desired_state_path_iterative(desired_path_wFG,m, var_name, shock_inds, system)
statesFG, obsFG, pseudoFG = forecast(system, s_0, shocks_path_wFG)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path_wFG,title="Change in r* since end of COVID19")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obsFG[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obsFG[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,statesFG[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudoFG[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudoFG[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "Main results/rstar_had_not_increased_wFG.pdf")   # saves the plot from p as a .pdf vector graphic


using Plots
# Create each subplot WITHOUT a legend
p1 = plot(plotdates, desired_path, title="Change in r* since end of COVID19", label="", legend=false)
plot!(plotdates, desired_path_wFG, lc=:red, label="")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p2 = plot(plotdates, obs[m.observables[:obs_nominalrate],:], title="Policy rate", label="", legend=false)
plot!(plotdates, obsFG[m.observables[:obs_nominalrate],:], lc=:red, label="")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p3 = plot(plotdates, obs[m.observables[:obs_gdpdeflator],:], title="Inflation", label="", legend=false)
plot!(plotdates, obsFG[m.observables[:obs_gdpdeflator],:], lc=:red, label="")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p4 = plot(plotdates, states[m.endogenous_states[:y_t],:], title="Output", label="", legend=false)
plot!(plotdates, statesFG[m.endogenous_states[:y_t],:], lc=:red, label="")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p5 = plot(plotdates, pseudo[m.pseudo_observables[:ExAnteRealRate],:], title="Ex-ante real rate", label="", legend=false)
plot!(plotdates, pseudoFG[m.pseudo_observables[:ExAnteRealRate],:], lc=:red, label="")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p6 = plot(plotdates, pseudo[m.pseudo_observables[:RealNaturalRate],:], title="Real natural rate", label="", legend=false)
plot!(plotdates, pseudoFG[m.pseudo_observables[:RealNaturalRate],:], lc=:red, label="")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

# Combine subplots and add a single legend above the first row, centered
plt = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), legend=:top, legendfontsize=10, size=(Int(960*1.5), Int(540*1.5)))

# # Add dummy series to the first subplot to show both legend entries
# plot!(plt[1], plotdates, desired_path, label="Without FG")
# plot!(plt[1], plotdates, desired_path_wFG, label="With FG")

savefig(plt, "Main results/rstar_had_not_increased_comparison_TFP.pdf")




# --- build your six plots (no legends anywhere in the grid)
zeroline = zeros(length(plotdates))

p1 = plot(plotdates, desired_path,  title="Change in r* since end of COVID-19",
          label="", legend=false)
plot!(p1, plotdates, desired_path_wFG, lc=:red, label="")
plot!(p1, plotdates, zeroline, lc=:black, lw=2, label="")

p2 = plot(plotdates, obs[m.observables[:obs_nominalrate],:], title="Policy rate",
          label="", legend=false)
plot!(p2, plotdates, obsFG[m.observables[:obs_nominalrate],:], lc=:red, label="")
plot!(p2, plotdates, zeroline, lc=:black, lw=2, label="")

p3 = plot(plotdates, obs[m.observables[:obs_gdpdeflator],:], title="Inflation",
          label="", legend=false)
plot!(p3, plotdates, obsFG[m.observables[:obs_gdpdeflator],:], lc=:red, label="")
plot!(p3, plotdates, zeroline, lc=:black, lw=2, label="")

p4 = plot(plotdates, states[m.endogenous_states[:y_t],:], title="Output",
          label="", legend=false)
plot!(p4, plotdates, statesFG[m.endogenous_states[:y_t],:], lc=:red, label="")
plot!(p4, plotdates, zeroline, lc=:black, lw=2, label="")

p5 = plot(plotdates, pseudo[m.pseudo_observables[:ExAnteRealRate],:], title="Ex-ante real rate",
          label="", legend=false)
plot!(p5, plotdates, pseudoFG[m.pseudo_observables[:ExAnteRealRate],:], lc=:red, label="")
plot!(p5, plotdates, zeroline, lc=:black, lw=2, label="")

p6 = plot(plotdates, pseudo[m.pseudo_observables[:RealNaturalRate],:], title="Real natural rate",
          label="", legend=false)
plot!(p6, plotdates, pseudoFG[m.pseudo_observables[:RealNaturalRate],:], lc=:red, label="")
plot!(p6, plotdates, zeroline, lc=:black, lw=2, label="")


# --- legend-only subplot
plegend = plot(legend = :top, framestyle = :none, grid = false,
               xticks = false, yticks = false, xlabel = "", ylabel = "",
               margin = 0mm, legendtitle = "")

# add dummy series for legend entries (NaN avoids plotting, only legend appears)
plot!(plegend, [0.0], [NaN], label = "Without FG",  foreground_color_legend = nothing)
plot!(plegend, [0.0], [NaN], label = "With FG",    lc = :red, foreground_color_legend = nothing)

# --- combine with your 3×2 grid of real plots
plt = plot(plegend, p1, p2, p3, p4, p5, p6;
           layout = @layout([a{0.10h}; grid(3,2)]),
           size   = (Int(960*1.5), Int(540*1.5)),
           top_margin = 4mm, bottom_margin = 4mm)

display(plt)
savefig(plt, "Main results/rstar_had_not_increased_comparison_TFP_legend.pdf")


# Alternative if r* did not increase post COVID19
desired_path =hlw_rstar.mean[end-20:259] .- hlw_rstar.mean[end-20]  #rstar_diff[end-16:end] # Desired path for the state variable
# desired_path = -desired_path
var_name =:Forward5YearRealNaturalRate


# Exclude all shocks whose names start with :rm_sh (including :rm_sh, :rm_shl1, ...)
shock_syms = collect(Iterators.filter(k -> !startswith(String(k), "rm_sh"), keys(m.exogenous_shocks)))

shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Change in r* since end of COVID19")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

savefig( "Main results/rstar_had_not_increased_WaggonerZha.pdf")   # saves the plot from p as a .pdf vector graphic



# --- Write plotted series to CSV ---
using CSV, DataFrames
df_out_norstarincrease_WZ = DataFrame(
        Date = plotdates,
        HLW_minus_DSGE = desired_path,
        PolicyRate = obs[m.observables[:obs_nominalrate], :],
        Inflation = obs[m.observables[:obs_gdpdeflator], :],
        Output = states[m.endogenous_states[:y_t], :],
        ExAnteRealRate = pseudo[m.pseudo_observables[:ExAnteRealRate], :],
        RealNaturalRate = pseudo[m.pseudo_observables[:RealNaturalRate], :]
)
CSV.write("Main results/what_if_rstar_had_not_increased_WaggonerZha.csv", df_out_norstarincrease_WZ)




# Alternative if r* did not increase post COVID19

desired_path =(hlw_rstar[end-20:259,10] .- hlw_rstar[end-20,10]) .- (hlw_rstar.mean[end-20:259] .- hlw_rstar.mean[end-20])  #rstar_diff[end-16:end] # Desired path for the state variable
# desired_path = -desired_path
var_name =:Forward5YearRealNaturalRate

# shock_syms = [  :b_liqtil_sh,   :b_liqp_sh,  :b_safetil_sh,  :b_safep_sh ] # Convenience yield shocks causing the difference
# shock_syms = keys(m.exogenous_shocks) # All shocks causing the difference

shock_syms = collect(Iterators.filter(k -> !startswith(String(k), "rm_sh"), keys(m.exogenous_shocks)))
shock_inds = repeat(reshape([m.exogenous_shocks[shock_name] for shock_name in shock_syms], :, 1), 1, length(desired_path))

shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo = forecast(system, s_0, shocks_path)
# --- Step 1: Compute IRFs for each shock ---
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate , :Forward5YearRealNaturalRate] # Output, Inflation, Policy Rate, R*
horizon = size(shocks_path, 2)
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

horizon = size(shocks_path, 2)
using Plots
p1 = plot(plotdates,desired_path,title="Change in r* since end of COVID19")
#p1 = plot(plotdates,states[m.endogenous_states[:b_liq_t],:],title="Combined liquidity shocks")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(plotdates,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(plotdates,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(plotdates,states[m.endogenous_states[:y_t],:],title="Output")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(plotdates,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate")#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(plotdates,pseudo[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate",)#
plot!(plotdates,zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

savefig( "Main results/rstar_increase_impact_of_FG_WaggonerZha.pdf")   # saves the plot from p as a .pdf vector graphic


