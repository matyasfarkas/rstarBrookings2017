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
# Original code (UNCHANGED computation)
# =========================

m = Model1010("ss20")

m <= DSGE.Setting(:data_vintage, "161223")
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))

system = DSGE.compute_system(m)

nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)

""""
obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, state_ind::Int, shock_inds::Vector{Int},
                                                    system::System{Float64})

Given a desired path `x` for state `state_ind` over `horizon` periods, and a vector of shock indices
`shock_inds` (one per period), back out the necessary values of shocks over those periods such that
s^i_{1:h} = x_{1:h}, using the model's forecast function.

Returns a matrix of required shocks of size (nshocks, horizon).
"""
function obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, m::AbstractDSGEModel, var_name::Symbol, shock_inds::Vector{Int},
                                                         system::System{Float64})
    horizon = length(x)
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    s_0 = zeros(nstates)
    shocks = zeros(nshocks, horizon)

    var_names, var_class, peg_ind =
        if var_name in keys(m.endogenous_states)
            m.endogenous_states, :states, m.endogenous_states[var_name]
        else
            m.observables, :obs,  m.observables[var_name]
        end

    for t in 1:horizon
        test_shocks = zeros(nshocks, horizon)
        test_shocks[shock_inds[t], t] = 1.0
        states, obs, _ = forecast(system, s_0, test_shocks)

        irf = (var_class == :states) ? states[peg_ind, t] : obs[peg_ind, t]

        prev_effect = 0.0
        if t > 1
            prev_shocks = shocks[:, 1:t-1]
            prev_states, prev_obs, _ = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
            prev_effect = (var_class == :states) ? prev_states[peg_ind, t] : prev_obs[peg_ind, t]
        end

        shocks[shock_inds[t], t] = (x[t] - prev_effect) / irf
    end

    return shocks
end

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]

# Load HLW real time estimates of R*
csv_path = joinpath(basepath, "Main results", "US", "Ex_post_real_rate_gaps.csv")
US_dataset = DataFrame(CSV.File(csv_path))

valid_idx = findall(row -> !ismissing(row[:date]) && !ismissing(row[:Target]), eachrow(US_dataset))
dates = US_dataset.date[valid_idx]

desired_path = skipmissing(US_dataset.Target[valid_idx]) |> collect
desired_path = -desired_path
var_name = :obs_nominalrate
horizon = size(desired_path, 1)

shock_name = :rm_sh
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

system = DSGE.zero_system_constants(system)
s_0 = zeros(nstates)

shocks = zeros(nshocks, horizon)
for t = 1:horizon
    var_value = desired_path[t]
    if var_class == :states
        var_value_att = var_value - obs[m.endogenous_states[var_name], t, m.exogenous_shocks[shock_name]]
        shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_state_value(var_value_att,
                                                                                var_names[var_name],
                                                                                exo[shock_name],
                                                                                system[:RRR])
    else
        var_value_att = var_value - obs[m.observables[var_name], t, m.exogenous_shocks[shock_name]]
        shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(var_value_att,
                                                                              var_names[var_name],
                                                                              exo[shock_name],
                                                                              system[:ZZ],
                                                                              system[:RRR])
    end

    states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)
end

plotvars = [:Forward5YearRealNaturalRate, :obs_nominalrate, :pi_t, :y_t, :ExAnteRealRate, :RealNaturalRate]
plotdates = Date.(dates[end-horizon+1:end], dateformat"mm/dd/yyyy")

# =========================
# Figure 1 + export
# =========================

p1 = plot(plotdates, desired_path, title="Targeted Path")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p2 = plot(plotdates, obs[m.observables[:obs_nominalrate], :, exo[shock_name]].*4, title="Policy rate")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p3 = plot(plotdates, obs[m.observables[:obs_gdpdeflator], :, exo[shock_name]].*4, title="Inflation")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p4 = plot(plotdates, states[m.endogenous_states[:y_t], :, exo[shock_name]], title="Output")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p5 = plot(plotdates, pseudo[m.pseudo_observables[:ExAnteRealRate], :, exo[shock_name]], title="Ex-ante real rate")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

p6 = plot(plotdates, pseudo[m.pseudo_observables[:RealNaturalRate], :, exo[shock_name]], title="Real natural rate")
plot!(plotdates, zeros(horizon,1), lc=:black, lw=2, label="")

plot(p1, p2, p3, p4, p5, p6, layout=(3,2), legend=false)
plot!(size=(960,540))

pdf_path = "Main results/What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf"
savefig(pdf_path)

export_chart_data(
    pdf_path,
    plotdates,
    Dict(
        "Targeted Path"     => vec(desired_path),
        "Zero line"         => vec(zeros(horizon,1)),
        "Policy rate"       => vec(obs[m.observables[:obs_nominalrate], :, exo[shock_name]].*4),
        "Inflation"         => vec(obs[m.observables[:obs_gdpdeflator], :, exo[shock_name]].*4),
        "Output"            => vec(states[m.endogenous_states[:y_t], :, exo[shock_name]]),
        "Ex-ante real rate" => vec(pseudo[m.pseudo_observables[:ExAnteRealRate], :, exo[shock_name]]),
        "Real natural rate" => vec(pseudo[m.pseudo_observables[:RealNaturalRate], :, exo[shock_name]]),
    )
)

# =========================
# Figure 2 + export (fix scoping via dict store)
# =========================

titlelist = ["Ex-ante real rate"]
plotvars = [:ExAnteRealRate]
zeroline = zeros(horizon)
plots_ar1 = Vector{Any}(undef, 1)
shock_idx = exo[:rm_sh]

series_store_exante = Dict{String,Vector{Float64}}()

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

    cf_series[1] = 0.0
    cf_short = cf_series[end-length(plotdates)+1:end].*4

    series_store_exante["Counterfactual HLW Ex-ante real rate"] = collect(cf_short)

    years = unique(year.(plotdates))
    year_tick_dates = [findfirst(d -> year(d) == y, plotdates) !== nothing ? plotdates[findfirst(d -> year(d) == y, plotdates)] : Date(string(y)*"-01-01") for y in years]
    year_tick_labels = [string(y) for y in years]
    xtick_tuple = (year_tick_dates, year_tick_labels)

    p = plot(plotdates, cf_short, label="Counterfactual HLW r", color=:blue, lw=2,
             xticks=xtick_tuple, ylim=(-5.5,4), title="Counterfactual HLW Ex-ante real rate")
    plot!(p, plotdates, zeroline, lc=:black, lw=2, label="")
    plot!(p, legend=false)
    plots_ar1[i] = p
end

plt = plot(plots_ar1[1], layout=(1,1), legend=false)
plot!(plt, size=(960,540))

pdf_path = "Main results/compare_baseline_vs_HLW_change_only_ex_ante_realrate.pdf"
savefig(plt, pdf_path)

export_chart_data(
    pdf_path,
    plotdates,
    merge(series_store_exante, Dict("Zero line" => vec(zeroline)))
)

# =========================
# Figure 3 + export (store both series during loop)
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
        cf_series = states[m.endogenous_states[v], :, shock_idx]
    elseif v in keys(m.pseudo_observables)
        cf_series = pseudo[m.pseudo_observables[v], :, shock_idx]
    else
        @warn "Variable $(v) not found in model observables/states/pseudo-observables."
        plots_arr[i] = plot(title=string(v), legend=false)
        continue
    end

    cf_series[1] = 0.0
    cf_short = cf_series[end-length(plotdates)+1:end]

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

pdf_path = "Main results/compare_baseline_vs_HLW_change_only_inflation_and_output.pdf"
savefig(plt2, pdf_path)

export_chart_data(
    pdf_path,
    plotdates,
    merge(series_store, Dict("Zero line" => vec(zeroline)))
)
