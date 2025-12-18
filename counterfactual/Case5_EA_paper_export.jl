using DSGE, Dates, DataFrames, OrderedCollections, HDF5, CSV, JLD2, FileIO, Statistics, ModelConstructors, LinearAlgebra
using Measures  # for mm margins
using Plots

#############################
# Written by: Matyas Farkas, IMF
# Based on the run_default.jl script of DSGE.jl and the test forecast driver file
# 25 August 2025
#############################

# =========================
# Export helpers (CSV next to PDF, same basename)
# =========================

function data_path_like_pdf(pdf_path::AbstractString; ext::String = ".csv")::String
    dir = dirname(pdf_path)
    base = splitext(basename(pdf_path))[1]
    return joinpath(dir, String(base) * ext)
end

# Make names valid Symbols for older DataFrames versions
function _colsym(name::AbstractString)::Symbol
    s = String(name)
    s = replace(s, r"[^A-Za-z0-9_]+" => "_")
    s = strip(s, '_')
    isempty(s) && (s = "series")
    # avoid starting with a digit
    if !isempty(s) && isdigit(s[1])
        s = "v_" * s
    end
    return Symbol(s)
end

function export_chart_data(pdf_path::AbstractString,
                           x,
                           series::Dict{String,<:AbstractVector})
    out = data_path_like_pdf(pdf_path; ext=".csv")
    mkpath(dirname(out))

    df = DataFrame()
    df[!, :x] = collect(x)

    for (k, v) in series
        df[!, _colsym(k)] = collect(v)
    end

    CSV.write(out, df)
    return out
end

# =========================
# Model setup (UNCHANGED)
# =========================

m = Model1010("ss20")

m <= DSGE.Setting(:data_vintage, "250113")
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))

mode_file = rawpath(m, "estimate", "paramsmode.h5")
DSGE.update!(m, h5read(mode_file, "params"))

system = DSGE.compute_system(m)

nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]

## Load in HLW real time estimates of R*
csv_path = joinpath(basepath, "Main results", "EA", "Ex_post_real_rate_gaps_EA.csv")
US_dataset = DataFrame(CSV.File(csv_path))

valid_idx = findall(row -> !ismissing(row[:date]) && !ismissing(row[:Target]), eachrow(US_dataset))
dates = US_dataset.date[valid_idx]

##### Alternative if realrate gap change of HLW was implemented using policy rate shocks alone
desired_path = skipmissing(US_dataset.Target[valid_idx]) |> collect   # Desired path for the state variable
desired_path = -desired_path
var_name = :obs_nominalrate
horizon = size(desired_path, 1)

shock_name = :rm_sh # MP shock implements the real rate gap change
shock_ind = 10

var_names, var_class =
    if var_name in keys(m.endogenous_states)
        m.endogenous_states, :states
    else
        m.observables, :obs
    end

exo     = m.exogenous_shocks
nshocks = size(system[:RRR], 2)
nstates = size(system[:TTT], 1)
nobs    = size(system[:ZZ], 1)
npseudo = size(system[:ZZ_pseudo], 1)

states = zeros(nstates, horizon, nshocks)
obs    = zeros(nobs,    horizon, nshocks)
pseudo = zeros(npseudo, horizon, nshocks)

# Set constant system matrices to 0
system = DSGE.zero_system_constants(system)

s_0 = zeros(nstates)

# Isolate single shock (UNCHANGED)
shocks = zeros(nshocks, horizon)
for t = 1:horizon
    var_value = desired_path[t]
    if var_class == :states
        var_value_att = var_value - obs[m.endogenous_states[var_name], t, m.exogenous_shocks[shock_name]]
        shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_state_value(var_value_att,
                                                                                var_names[var_name],
                                                                                exo[shock_name],
                                                                                system[:RRR])
    else # == :obs
        var_value_att = var_value - obs[m.observables[var_name], t, m.exogenous_shocks[shock_name]]
        shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(var_value_att,
                                                                              var_names[var_name],
                                                                              exo[shock_name],
                                                                              system[:ZZ],
                                                                              system[:RRR])
    end

    # Iterate state space forward (UNCHANGED)
    states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)
end

# =========================
# Figure 1 + export (UNCHANGED plotting, add export only)
# =========================

plotvars = [:Forward5YearRealNaturalRate, :obs_nominalrate, :pi_t, :y_t, :ExAnteRealRate, :RealNaturalRate]
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

p1 = plot(plotdates, desired_path, title="Targeted Path")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")
p2 = plot(plotdates, obs[m.observables[:obs_nominalrate], :, exo[shock_name]], title="Policy rate")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")
p3 = plot(plotdates, obs[m.observables[:obs_gdpdeflator].*4, :, exo[shock_name]], title="Inflation")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")
p4 = plot(plotdates, states[m.endogenous_states[:y_t], :, exo[shock_name]], title="Output")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")
p5 = plot(plotdates, pseudo[m.pseudo_observables[:ExAnteRealRate], :, exo[shock_name]], title="Ex-ante real rate")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")
p6 = plot(plotdates, pseudo[m.pseudo_observables[:RealNaturalRate], :, exo[shock_name]], title="Real natural rate")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

plot(p1, p2, p3, p4, p5, p6, layout=(3,2), legend=false)
plot!(size=(960,540))

pdf_path = "Main results/EA_What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf"
savefig(pdf_path)

export_chart_data(
    pdf_path,
    plotdates,
    Dict(
        "Targeted Path"     => vec(desired_path),
        "Zero line"         => vec(zeros(horizon,1)),
        "Policy rate"       => vec(obs[m.observables[:obs_nominalrate], :, exo[shock_name]]),
        "Inflation"         => vec(obs[m.observables[:obs_gdpdeflator].*4, :, exo[shock_name]]),
        "Output"            => vec(states[m.endogenous_states[:y_t], :, exo[shock_name]]),
        "Ex-ante real rate" => vec(pseudo[m.pseudo_observables[:ExAnteRealRate], :, exo[shock_name]]),
        "Real natural rate" => vec(pseudo[m.pseudo_observables[:RealNaturalRate], :, exo[shock_name]]),
    )
)

# =========================
# Figure 2: Ex-ante real rate only + export
# (Fix: store cf_short inside loop; do NOT reference it outside)
# =========================

titlelist = ["Ex-ante real rate"]
plotvars = [:ExAnteRealRate]
zeroline = zeros(horizon)
plots_ar1 = Vector{Any}(undef, 1)
shock_idx = exo[:rm_sh]  # 10 for your model

series_store_exante = Dict{String,Vector{Float64}}()  # ADDED: store plotted series

for (i, v) in enumerate(plotvars)
    if v == :pi_t
        cf_series = pseudo[m.pseudo_observables[:π_t], :, shock_idx]
    elseif v in keys(m.observables)
        cf_series = obs[m.observables[v], :, shock_idx]
    elseif v in keys(m.endogenous_states)
        cf_series = states[m.endogenous_states[v], :, shock_idx]
    elseif v in keys(m.pseudo_observables)
        cf_series = pseudo[m.pseudo_observables[v], :, shock_idx]
    else
        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
        plots_ar1[i] = plot(title=string(v), legend=false)
        continue
    end

    cf_series[1] = 0.0 # Align first value to zero
    cf_short = cf_series[end-length(plotdates)+1:end].*4

    # Store for export (ADDED ONLY)
    series_store_exante["Counterfactual HLW Ex-ante real rate"] = collect(cf_short)

    years = unique(year.(plotdates))
    year_tick_dates = [findfirst(d -> year(d) == y, plotdates) !== nothing ? plotdates[findfirst(d -> year(d) == y, plotdates)] : Date(string(y)*"-01-01") for y in years]
    year_tick_labels = [string(y) for y in years]
    xtick_tuple = (year_tick_dates, year_tick_labels)

    p1 = plot(plotdates, cf_short, label="Counterfactual HLW r", color=:blue, lw=2,
              xticks=xtick_tuple, ylim=(-5.5,4), title="Counterfactual HLW Ex-ante real rate")
    plot!(p1, plotdates, zeroline, lc=:black, lw=2, label="")
    plot!(p1, legend=false)
    plots_ar1[i] = p1
end

plt = plot(plots_ar1[1], layout=(1,1), legend=false)
plot!(plt, size=(960,540))

pdf_path = "Main results/compare_EA_baseline_vs_HLW_change_only_ex_ante_realrate.pdf"
savefig(plt, pdf_path)

# Export (uses stored dict; never references cf_short directly)
export_chart_data(
    pdf_path,
    plotdates,
    merge(series_store_exante, Dict("Zero line" => vec(zeroline)))
)

# =========================
# Figure 3: Inflation and Output + export
# (Already storing cf_short inside loop, keep that)
# =========================

titlelist = ["Inflation", "Output"]
plotvars = [:obs_gdpdeflator, :y_t]
zeroline = zeros(horizon)
plots_arr = Vector{Any}(undef, length(plotvars))
series_store = Dict{String,Vector{Float64}}()

for (i, v) in enumerate(plotvars)
    if v == :pi_t
        cf_series = pseudo[m.pseudo_observables[:π_t], :, shock_idx]
    elseif v in keys(m.observables)
        cf_series = obs[m.observables[v], :, shock_idx].*4
    elseif v in keys(m.endogenous_states)
        cf_series = states[m.endogenous_states[v], :, shock_idx]./2
    elseif v in keys(m.pseudo_observables)
        cf_series = pseudo[m.pseudo_observables[v], :, shock_idx]
    else
        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
        plots_arr[i] = plot(title=string(v), legend=false)
        continue
    end

    cf_series[1] = 0.0
    cf_short = cf_series[end-length(plotdates)+1:end]

    # Store for export (ADDED ONLY)
    series_store[titlelist[i]] = collect(cf_short)

    years = unique(year.(plotdates))
    year_tick_dates = [findfirst(d -> year(d) == y, plotdates) !== nothing ? plotdates[findfirst(d -> year(d) == y, plotdates)] : Date(string(y)*"-01-01") for y in years]
    year_tick_labels = [string(y) for y in years]
    xtick_tuple = (year_tick_dates, year_tick_labels)

    p = plot(plotdates, cf_short, label="Counterfactual HLW Real Rate Gap", color=:blue, lw=2,
             title=titlelist[i], xticks=xtick_tuple)
    plot!(p, plotdates, zeroline, lc=:black, lw=2, label="")
    plot!(p, legend=false)
    plots_arr[i] = p
end

while length(plots_arr) < 2
    push!(plots_arr, plot(title="", legend=false))
end

plt2 = plot(plots_arr[1], plots_arr[2], layout=(1,2), legend=false)
plot!(plt2, size=(960,540))

pdf_path = "Main results/compare_EA_baseline_vs_HLW_change_only_inflation_and_output.pdf"
savefig(plt2, pdf_path)

export_chart_data(
    pdf_path,
    plotdates,
    merge(series_store, Dict("Zero line" => vec(zeroline)))
)
