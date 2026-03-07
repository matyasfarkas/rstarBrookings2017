using DSGE, Dates, DataFrames, OrderedCollections, CSV, HDF5, JLD2, FileIO
using Statistics, ModelConstructors, LinearAlgebra
using Plots
using Measures

#------------------------------------------------------------------------------
# Helpers
#------------------------------------------------------------------------------

"""
    var_loc(m, v::Symbol) -> (cls::Symbol, idx::Int)

Return where `v` lives and its index:
- `:states`  => m.endogenous_states[v]
- `:obs`     => m.observables[v]
- `:pseudo`  => m.pseudo_observables[v]
"""
function var_loc(m, v::Symbol)
    if v in keys(m.endogenous_states)
        return (:states, m.endogenous_states[v])
    elseif v in keys(m.observables)
        return (:obs, m.observables[v])
    elseif v in keys(m.pseudo_observables)
        return (:pseudo, m.pseudo_observables[v])
    else
        error("""
        Variable $(v) not found.

        Check exact symbol spelling, e.g. :pi_t versus :π_t.

        Matching state keys:
        $(filter(x -> occursin(lowercase(String(v)), lowercase(String(x))), collect(keys(m.endogenous_states))))

        Matching observable keys:
        $(filter(x -> occursin(lowercase(String(v)), lowercase(String(x))), collect(keys(m.observables))))

        Matching pseudo-observable keys:
        $(filter(x -> occursin(lowercase(String(v)), lowercase(String(x))), collect(keys(m.pseudo_observables))))
        """)
    end
end

"""
    get_series(m, states, obs, pseudo, v::Symbol; shock_idx=nothing) -> Vector

Extract a 1D series for variable `v` from the appropriate container.
Handles both 2D arrays [var, t] and 3D arrays [var, t, shock].
"""
function get_series(m, states, obs, pseudo, v::Symbol; shock_idx::Union{Nothing,Int}=nothing)
    cls, idx = var_loc(m, v)
    A = cls == :states ? states : cls == :obs ? obs : pseudo

    if ndims(A) == 2
        return vec(A[idx, :])
    elseif ndims(A) == 3
        shock_idx === nothing && error("Array for $(v) is 3D; pass shock_idx=...")
        return vec(A[idx, :, shock_idx])
    else
        error("Unsupported array dimension $(ndims(A)) for $(v)")
    end
end

"""
    year_xticks(dts::Vector{Date}) -> (tick_dates, tick_labels)
"""
function year_xticks(dts::Vector{Date})
    yrs = unique(year.(dts))
    tick_dates = Date[]
    for y in yrs
        j = findfirst(d -> year(d) == y, dts)
        push!(tick_dates, dts[j])
    end
    return (tick_dates, string.(yrs))
end

"""
    plot_with_zero(dts, y; title="", ylim=nothing, xticks=nothing)

Make a plot with a zero line.
"""
function plot_with_zero(dts::Vector{Date}, y::AbstractVector;
                        title::AbstractString = "",
                        ylim = nothing,
                        xticks = nothing)
    p = plot(dts, y; title=title, lw=2, legend=false)
    plot!(p, dts, zeros(length(dts)); lc=:black, lw=2, label="")
    ylim !== nothing && plot!(p; ylim=ylim)
    xticks !== nothing && plot!(p; xticks=xticks)
    return p
end

"""
    save_pdf_and_csv(fig, pdf_path, dts; series)

Save figure as PDF and export plotted series to CSV with same stem name.
"""
function save_pdf_and_csv(fig, pdf_path::AbstractString, dts::Vector{Date};
                          series::AbstractDict{Symbol,<:AbstractVector})
    savefig(fig, pdf_path)

    csv_path = replace(pdf_path, r"\.pdf$" => ".csv")
    if csv_path == pdf_path
        csv_path *= ".csv"
    end

    df = DataFrame(Date = dts)
    for (k, v) in series
        @assert length(v) == length(dts) "Column $(k) length $(length(v)) ≠ dates length $(length(dts))"
        df[!, k] = collect(v)
    end
    CSV.write(csv_path, df)
    return nothing
end

#------------------------------------------------------------------------------
# Shock backout
#------------------------------------------------------------------------------


""""
obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, state_ind::Int, shock_inds::Vector{Int},
                                                    system::System{Float64})

Given a desired path `x` for state `state_ind` over `horizon` periods, and a vector of shock indices
`shock_inds` (one per period), back out the necessary values of shocks over those periods such that
s^i_{1:h} = x_{1:h}, using the model's forecast function.

Returns a matrix of required shocks of size (nshocks, horizon).
"""
function obtain_shocks_from_desired_state_path_iterative(x::Vector{Float64}, m::AbstractDSGEModel,var_name::Symbol, shock_inds::Vector{Int},
                                                         system::System{Float64})


    horizon = length(x)
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    s_0 = zeros(nstates)
    shocks = zeros(nshocks, horizon)
    var_names, var_class, peg_ind =
    if var_name in keys(m.endogenous_states)
        m.endogenous_states, :states, m.endogenous_states[var_name]

    elseif var_name in keys(m.observables)
        m.observables, :obs,  m.observables[var_name]
    elseif var_name in keys(m.pseudo_observables)
        m.pseudo_observables, :pseudo, m.pseudo_observables[var_name]
    end


    for t in 1:horizon
        
        # Compute IRF for a unit shock at time t for the specified shock
        test_shocks = zeros(nshocks, horizon)
        test_shocks[shock_inds[t], t] = 1.0
        states, obs, pseudo = forecast(system, s_0, test_shocks)
                if var_class == :states
                    irf = states[peg_ind, t] # Impact of a unit shock at t on state at t
                elseif var_class == :obs
                    irf = obs[peg_ind,t] # Impact of a unit shock at t on obs at t
                elseif var_class == :pseudo
                    irf = pseudo[peg_ind,t] # Impact of a unit shock at t on pseudo-obs at t
                end
        # Compute effect of previous shocks
        prev_effect = 0.0
        if t > 1
            prev_shocks = shocks[:, 1:t-1]
            prev_states, prev_obs, prev_pseudo = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
            if var_class == :states
                prev_effect = prev_states[peg_ind, t]
            elseif var_class == :obs
                prev_effect = prev_obs[peg_ind,t] # Impact of a unit shock at t on obs at t
            elseif var_class == :pseudo
                prev_effect = prev_pseudo[peg_ind,t] # Impact of a unit shock at t on pseudo-obs at t   
            end
        end

        # Required shock at time t to achieve desired value
        shocks[shock_inds[t], t] = (x[t] - prev_effect) / irf
    end

    return shocks
end
#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------



# Model init
m = Model1010("ss20")
m <= DSGE.Setting(:data_vintage, "250825")

mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")
datafolder = joinpath(basepath, "dsge", "output_data")
mode_file = joinpath(datafolder, "m1010","ss20","estimate","raw", "paramsmode_vint=250825.h5")
specify_mode!(m, mode_file)




m <= DSGE.Setting(:date_forecast_start, quartertodate("2024-Q4"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))

system = DSGE.compute_system(m)
system = DSGE.zero_system_constants(system)

nstates = size(system[:TTT], 1)
s0 = zeros(nstates)

# Data path
basepath = abspath(joinpath(saveroot, ".."))
csv_path = joinpath(basepath, "Main results", "US", "Ex_post_real_rate_gaps.csv")
US = DataFrame(CSV.File(csv_path))

valid = findall(r -> !ismissing(r[:date]) && !ismissing(r[:actual_MP_stance]), eachrow(US))
dates_raw = US.date[valid]
desired_path = -collect(skipmissing(US.actual_MP_stance[valid]))/4
# desired_path= [-0.5347660219210864, -0.42580683141933096, -0.2848610819447991, -0.1902461026228856, -0.13293768204673204, -0.09838621322188786, -0.0771225475192518, -0.06375336981838424, -0.05522381954465688, -0.04973500673305435, -0.046172815925899724, -0.04381636832118916, -0.04218579829272206, -0.04095821113356902, -0.03991776319380584, -0.03892372752154872, -0.037888760326541024, -0.03676344848609328, -0.03552502983049402]

var_name   = :ExAnteRealRate # :obs_nominalrate, :ExAnteRealRate, :y_t
shock_name = :rm_sh
horizon    = length(desired_path)

# Parse dates
dates_sub = dates_raw[end-horizon+1:end]
plotdates = eltype(dates_sub) <: Date ? collect(dates_sub) :
            Date.(dates_sub, dateformat"mm/dd/yyyy")

# Back out shocks
exo = m.exogenous_shocks
shock_inds = fill(exo[shock_name], horizon)

shocks = obtain_shocks_from_desired_state_path_iterative(
    Float64.(desired_path),
    m,
    var_name,
    shock_inds,
    system
)


states, obs, pseudo = forecast(system, s0, shocks)
# shockks = zeros(29, 19)
# shockks[exo[shock_name], 1] = -1
# states, obs, pseudo = forecast(system, s0, shockks)

#------------------------------------------------------------------------------
# Figure 1: 3x2 grid
#------------------------------------------------------------------------------

targeted   = Float64.(desired_path)
policy     = get_series(m, states, obs, pseudo, :obs_nominalrate)
infl       = 4 .* get_series(m, states, obs, pseudo, :obs_corepce)
output     = get_series(m, states, obs, pseudo, :y_t)
exante_ann = get_series(m, states, obs, pseudo, :ExAnteRealRate)
obs_gdp       = 4 .* get_series(m, states, obs, pseudo, :obs_gdp)

p1 = plot_with_zero(plotdates, targeted,   title="Targeted Path")
p2 = plot_with_zero(plotdates, policy,     title="Policy rate")
p3 = plot_with_zero(plotdates, infl,   title="Inflation")
p4 = plot_with_zero(plotdates, output,     title="Output")
p5 = plot_with_zero(plotdates, exante_ann, title="Ex-ante real rate")
p6 = plot_with_zero(plotdates, obs_gdp,    title="GDP growth (%, qoq annualized)")

fig1 = plot(p1, p2, p3, p4, p5, p6; layout=(3,2), legend=false, size=(960,540))

pdf_path1 = joinpath(saveroot, "Final Paper","Counterfactual", "What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf")
save_pdf_and_csv(fig1, pdf_path1, plotdates; series=OrderedDict(
    :TargetedPath                => targeted,
    :PolicyRate                  => policy,
    :Inflation                   => infl,
    :Output                      => output,
    :ExAnteRealRate_Annualized   => exante_ann,
    :obs_gdp                     => obs_gdp
))

#------------------------------------------------------------------------------
# Figure 2: Ex-ante real rate only
#------------------------------------------------------------------------------

xt = year_xticks(plotdates)
y2 = copy(exante_ann)
y2[1] = 0.0

fig2 = plot_with_zero(
    plotdates,
    y2;
    title  = "Counterfactual HLW Ex-ante real rate",
    ylim   = (-5.5, 4),
    xticks = xt
)
plot!(fig2; size=(960,540))

pdf_path2 = joinpath(saveroot, "Final Paper","Counterfactual", "compare_baseline_vs_HLW_change_only_ex_ante_realrate.pdf")
save_pdf_and_csv(fig2, pdf_path2, plotdates; series=OrderedDict(
    :ExAnteRealRate_Annualized => y2
))

#------------------------------------------------------------------------------
# Figure 3: Inflation + Output
#------------------------------------------------------------------------------

y_inf = copy(infl); y_inf[1] = 0.0
y_out = copy(output);   y_out[1] = 0.0

p_inf = plot_with_zero(plotdates, y_inf; title="Inflation (%, yoy)", xticks=xt)
p_out = plot_with_zero(plotdates, y_out; title="Output (% dev from SS)", xticks=xt)

fig3 = plot(p_inf, p_out; layout=(1,2), legend=false, size=(960,540))

pdf_path3 = joinpath(saveroot, "Final Paper","Counterfactual", "Paper_Compare_baseline_vs_HLW_change_only_inflation_and_output.pdf")
save_pdf_and_csv(fig3, pdf_path3, plotdates; series=OrderedDict(
    :InflationYoY => y_inf,
    :Output       => y_out
))