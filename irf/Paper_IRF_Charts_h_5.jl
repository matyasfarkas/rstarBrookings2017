using DSGE
using Plots
using XLSX
using Dates

# =========================
# Excel export
# =========================
path = dirname(@__FILE__())

function sheet_from_pdf(pdf_path::AbstractString; maxlen::Int=31)::String
    base = basename(pdf_path)
    name = splitext(base)[1]
    name = String(name)
    name = replace(name, r"[:\\\/\?\*\[\]]" => "_")
    name = strip(name)
    isempty(name) && (name = "Sheet")
    if lastindex(name) > maxlen
        name = first(name, maxlen)
    end
    return String(name)
end

function export_series_xlsx(filepath::AbstractString,
                            sheetname_in::AbstractString,
                            t,
                            series::Dict{String,<:AbstractVector})
    mkpath(dirname(filepath))
    sheetname = String(sheetname_in)

    colnames = ["t"; collect(keys(series))]
    cols = Any[collect(t)]
    for nm in colnames[2:end]
        push!(cols, series[nm])
    end

    mode = isfile(filepath) ? "rw" : "w"

    XLSX.openxlsx(filepath, mode=mode) do xf
        if !(sheetname in XLSX.sheetnames(xf))
            XLSX.addsheet!(xf, sheetname)
        end
        sh = xf[sheetname]
        XLSX.writetable!(sh, cols, colnames; anchor_cell=XLSX.CellRef("A1"))
        sh[1, length(colnames)+2] = "exported_at"
        sh[1, length(colnames)+3] = string(now())
    end

    return filepath
end

# =========================
# Paths / model setup
# =========================
saveroot = dirname(@__FILE__)
mypath = @__DIR__

idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")
datafolder = joinpath(basepath, "dsge", "output_data")

horizon     = 20
peg_horizon = 5                     # <-- CHANGED: 5-quarter peg
tgrid       = collect(1:horizon)
xlsx_out    = joinpath(path, "irf", "Paper IRFs.xlsx")

m = Model1010("ss20")
m <= DSGE.Setting(:data_vintage, "250825")

params_mode = load_draws(m, :mode)
DSGE.update!(m, params_mode)
DSGE.steadystate!(m)

mode_file = joinpath(datafolder, "m1010", "ss20", "estimate", "raw", "paramsmode_vint=250825.h5")
specify_mode!(m, mode_file)

system = compute_system(m)
system = DSGE.zero_system_constants(system)

nstates = size(system[:TTT], 1)
s_0 = zeros(nstates)

# =========================
# Plot helpers
# =========================
my_xticks = 0:4:20

panel_ylim(i::Int) =
    i == 1 ? (-3.0, 3.0)  :
    i == 2 ? (-1.5, 1.5)  :
    i == 3 ? (-0.1, 0.05) :
    i == 4 ? (-1.0, 1.0)  :
    i == 5 ? (-0.1, 0.1)  :
    i == 6 ? (-1.5, 1.5)  :
    error("No calibrated ylim for panel $i")

function scale_for_plot(var_sym::Symbol, x::AbstractVector)
    if var_sym in (:obs_nominalrate, :obs_corepce,
                   :Forward5YearRealNaturalRate, :RealNaturalRate,
                   :ExAnteRealRate)
        return 4 .* x
    else
        return x
    end
end

function get_series(m, states, obs, pseudo, var_sym::Symbol, PlotT::Int)
    if haskey(m.observables, var_sym)
        return vec(obs[m.observables[var_sym], 1:PlotT])
    elseif haskey(m.endogenous_states, var_sym)
        return vec(states[m.endogenous_states[var_sym], 1:PlotT])
    elseif haskey(m.pseudo_observables, var_sym)
        return vec(pseudo[m.pseudo_observables[var_sym], 1:PlotT])
    else
        error("Variable $(var_sym) not found in observables, states, or pseudo observables.")
    end
end

function collect_irfs(m, system, s_0, plotvars::Vector{Symbol},
                      shock_syms::Vector{Symbol}, PlotT::Int;
                      shock_scales::Dict{Symbol,Float64}=Dict{Symbol,Float64}())
    nvars   = length(plotvars)
    nshocks = length(shock_syms)
    irfmat  = zeros(PlotT, nvars, nshocks)
    shockmat = zeros(size(system[:RRR], 2), PlotT, nshocks)

    for (j, shock_sym) in enumerate(shock_syms)
        shocks = zeros(size(system[:RRR], 2), PlotT)
        shocks[m.exogenous_shocks[shock_sym], 1] = get(shock_scales, shock_sym, 1.0)
        states, obs, pseudo, _ = forecast(system, s_0, shocks)

        for (i, var_sym) in enumerate(plotvars)
            irfmat[:, i, j] .= get_series(m, states, obs, pseudo, var_sym, PlotT)
        end
        shockmat[:, :, j] .= shocks
    end

    return irfmat, shockmat
end

function build_fg_solution(irfmatR::Array{Float64,3},
                           irfmatFG::Array{Float64,3},
                           peg_len::Int;
                           peg_level::Float64 = -1.0/4)

    PlotT = size(irfmatFG, 1)
    nvars = size(irfmatFG, 2)

    # Target path for the policy rate net of the contemporaneous rm_sh contribution
    FG_vec = fill(peg_level, peg_len) .- irfmatR[1:peg_len, 1, 1]

    # Policy-rate response matrix to rm_shl1...rm_shl5
    R_mp_mat = zeros(peg_len, peg_len)
    for s in 1:peg_len
        R_mp_mat[:, s] .= irfmatFG[1:peg_len, 1, s]
    end

    fg_weights = R_mp_mat \ FG_vec

    # Total IRF = contemporaneous rm_sh + weighted FG shocks
    total_irf = copy(irfmatR[:, :, 1])
    for s in 1:peg_len
        total_irf .+= fg_weights[s] .* irfmatFG[:, :, s]
    end

    return total_irf, fg_weights
end

function plot_star_panel!(plt, impact_weight::Float64, fg_weights::AbstractVector, PlotT::Int)
    # quarter 1: contemporaneous rm_sh
    scatter!(plt, [1], [4 * impact_weight],
        marker=:star5,
        markersize=7,
        markercolor=:black,
        markerstrokecolor=:black,
        label="",
        xticks=my_xticks
    )

    # quarters 2..6: rm_shl1..rm_shl5
    scatter!(plt, collect(2:length(fg_weights)+1), 4 .* fg_weights,
        marker=:star5,
        markersize=7,
        markercolor=:red,
        markerstrokecolor=:red,
        label="",
        xticks=my_xticks
    )

    plot!(plt, 1:PlotT, zeros(PlotT), lc=:black, lw=1, label="")
end

function build_star_export_dict(impact_weight::Float64,
                                fg_weights::AbstractVector,
                                total_irf::Matrix{Float64},
                                plotvars::Vector{Symbol},
                                PlotT::Int)
    d = Dict{String,Vector{Float64}}()

    # shocks shown in panel 1
    d["Shock q1 rm_sh"] = vcat([impact_weight], zeros(PlotT-1))
    for j in 1:length(fg_weights)
        series = zeros(PlotT)
        series[j+1] = fg_weights[j]   # q2..q6
        d["Shock q$(j+1) rm_shl$(j)"] = series
    end

    # plotted series
    for (i, var_sym) in enumerate(plotvars)
        nm = string(var_sym)
        d[nm] = vec(scale_for_plot(var_sym, total_irf[:, i]))
    end

    d["Zero line"] = zeros(PlotT)
    return d
end

function draw_six_panel_star_figure(plotvars::Vector{Symbol},
                                    total_irf::Matrix{Float64},
                                    impact_weight::Float64,
                                    fg_weights::AbstractVector;
                                    titles::Vector{String},
                                    PlotT::Int)

    p = plot(layout=(3,2), size=(1200,800), legend=false)
    TT = 1:PlotT

    for i in 1:(length(plotvars)+1)
        if i == 1
            plot_star_panel!(p[i], impact_weight, fg_weights, PlotT)
        else
            var_sym = plotvars[i-1]
            y = vec(scale_for_plot(var_sym, total_irf[:, i-1]))
            plot!(p[i], TT, y, lw=2, label="", xticks=my_xticks)
            plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
        end

        title!(p[i], titles[i])
        ylabel!(p[i], "%")
        xlabel!(p[i], "Quarter")
        ylims!(p[i], panel_ylim(i))
    end

    return p
end

function run_star_fg_exercise(m, system, s_0;
                              plotvars::Vector{Symbol},
                              rm_scale::Float64,
                              pdf_file::String,
                              PlotT::Int=horizon,
                              peg_len::Int=peg_horizon)

    # contemporaneous rm_sh
    irfmatR, shockmatR = collect_irfs(
        m, system, s_0,
        plotvars,
        [:rm_sh],
        PlotT;
        shock_scales=Dict(:rm_sh => rm_scale)
    )

    # FG shocks: rm_shl1..rm_shl5  --> stars at quarters 2..6
    fg_syms = [:rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5]
    irfmatFG, _ = collect_irfs(
        m, system, s_0,
        plotvars,
        fg_syms,
        PlotT
    )

    total_irf, fg_weights = build_fg_solution(irfmatR, irfmatFG, peg_len; peg_level=-1.0/4)

    impact_weight = shockmatR[m.exogenous_shocks[:rm_sh], 1, 1]

    titles = ["Anticipated Policy Innovations (APR)",
              "Policy rate (APR)",
              "Inflation (%, qoq annualized)",
              "Output (% dev from SS)",
              "r* (Forward 5-year real natural rate, APR)",
              "Ex-ante real rate (APR)"]

    p = draw_six_panel_star_figure(plotvars, total_irf, impact_weight, fg_weights;
                                   titles=titles, PlotT=PlotT)

    pdf_path = joinpath(saveroot, "Final Paper", "figures", pdf_file)
    savefig(p, pdf_path)

    export_series_xlsx(
        xlsx_out,
        sheet_from_pdf(pdf_path),
        collect(1:PlotT),
        build_star_export_dict(impact_weight, fg_weights, total_irf, plotvars, PlotT)
    )

    return p, total_irf, fg_weights
end

# =========================
# 1) Contemporaneous MP innovations only
#    5-quarter peg with rm_sh in quarters 1:5
# =========================
begin
    shock_name = :rm_sh
    var_name   = :obs_nominalrate
    var_value  = -1.0/4

    exo     = m.exogenous_shocks
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    nobs    = size(system[:ZZ], 1)
    npseudo = size(system[:ZZ_pseudo], 1)

    states = zeros(nstates, horizon, nshocks)
    obs    = zeros(nobs,    horizon, nshocks)
    pseudo = zeros(npseudo, horizon, nshocks)

    shocks = zeros(nshocks, horizon)

    for t in 1:peg_horizon
        var_value_att = var_value - obs[m.observables[var_name], t, exo[shock_name]]
        shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(
            var_value_att,
            m.observables[var_name],
            exo[shock_name],
            system[:ZZ],
            system[:RRR]
        )

        states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ =
            forecast(system, s_0, shocks)
    end

    y = shocks[m.exogenous_shocks[:rm_sh], :]
    x = 1:horizon
    idx = 1:peg_horizon

    p1 = plot(
        x[idx],
        4 .* y[idx],
        seriestype=:scatter,
        marker=:star5,
        markersize=6,
        markercolor=:black,
        title="Contemporaneous Policy Innovations (APR)",
        label="",
        xticks=my_xticks
    )
    ylabel!(p1, "%")
    xlabel!(p1, "Quarter")
    plot!(p1, x, zeros(horizon), lc=:black, lw=2, label="")
    ylims!(p1, panel_ylim(1))

    p2 = plot(1:horizon, 4 .* obs[m.observables[:obs_nominalrate], :, m.exogenous_shocks[:rm_sh]],
        title="Policy rate (APR)", ylims=panel_ylim(2), xticks=my_xticks)
    ylabel!(p2, "%"); xlabel!(p2, "Quarter"); plot!(p2, zeros(horizon), lc=:black, lw=2, label="")

    p3 = plot(1:horizon, 4 .* obs[m.observables[:obs_corepce], :, m.exogenous_shocks[:rm_sh]],
        title="Inflation (%, qoq annualized)", ylims=panel_ylim(3), xticks=my_xticks)
    ylabel!(p3, "%"); xlabel!(p3, "Quarter"); plot!(p3, zeros(horizon), lc=:black, lw=2, label="")

    p4 = plot(1:horizon, states[m.states[:y_t], :, m.exogenous_shocks[:rm_sh]],
        title="Output (% dev from SS)", ylims=panel_ylim(4), xticks=my_xticks)
    ylabel!(p4, "%"); xlabel!(p4, "Quarter"); plot!(p4, zeros(horizon), lc=:black, lw=2, label="")

    p5 = plot(1:horizon, 4 .* pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate], :, m.exogenous_shocks[:rm_sh]],
        title="r* (Forward 5-year real natural rate, APR)", ylims=panel_ylim(5), xticks=my_xticks)
    ylabel!(p5, "%"); xlabel!(p5, "Quarter"); plot!(p5, zeros(horizon), lc=:black, lw=2, label="")

    p6 = plot(1:horizon, 4 .* pseudo[m.pseudo_observables[:ExAnteRealRate], :, m.exogenous_shocks[:rm_sh]],
        title="Ex-ante real rate (APR)", ylims=panel_ylim(6), xticks=my_xticks)
    ylabel!(p6, "%"); xlabel!(p6, "Quarter"); plot!(p6, zeros(horizon), lc=:black, lw=2, label="")

    p = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(1200,800), legend=false)

    pdf_path = joinpath(saveroot, "Final Paper", "figures", "Interest_rate_peg_IRF_Policy_rate_with_MP_shock.pdf")
    savefig(p, pdf_path)

    export_series_xlsx(
        xlsx_out,
        sheet_from_pdf(pdf_path),
        tgrid,
        Dict(
            "Monetary policy shock (rm_sh)"        => vec(shocks[m.exogenous_shocks[:rm_sh], :]),
            "Policy rate (obs_nominalrate)"       => vec(4 .* obs[m.observables[:obs_nominalrate], :, m.exogenous_shocks[:rm_sh]]),
            "Inflation (obs_corepce)"             => vec(4 .* obs[m.observables[:obs_corepce], :, m.exogenous_shocks[:rm_sh]]),
            "Output (obs_gdp)"                    => vec(states[m.states[:y_t], :, m.exogenous_shocks[:rm_sh]]),
            "r* (Forward5YearRealNaturalRate)"    => vec(4 .* pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate], :, m.exogenous_shocks[:rm_sh]]),
            "Ex-ante real rate (ExAnteRealRate)"  => vec(4 .* pseudo[m.pseudo_observables[:ExAnteRealRate], :, m.exogenous_shocks[:rm_sh]]),
            "Zero line"                           => zeros(horizon),
        )
    )
end

# =========================
# Common plotvars for the star-panel FG exercises
# =========================
plotvars = [:obs_nominalrate,
            :obs_corepce,
            :obs_gdp,
            :Forward5YearRealNaturalRate,
            :ExAnteRealRate]

# =========================
# 2) FG exercise:
#    rm_sh at t=1 is nonzero, plus 5-quarter peg from rm_shl1..rm_shl5
#    stars in first panel at quarters 1..6
# =========================
begin
    rm_scale = -1 / 0.45
    run_star_fg_exercise(
        m, system, s_0;
        plotvars=plotvars,
        rm_scale=rm_scale,
        pdf_file="interest_rate_peg_with_contemporaneous_innovation_1_and_rest_FG_no_increase.pdf",
        PlotT=horizon,
        peg_len=peg_horizon
    )
end

# =========================
# 3) FG exercise:
#    rm_sh at t=1 is zero, peg implemented by rm_shl1..rm_shl5 only
#    stars in first panel at quarters 1..6
# =========================
begin
    rm_scale = 0.0
    run_star_fg_exercise(
        m, system, s_0;
        plotvars=plotvars,
        rm_scale=rm_scale,
        pdf_file="interest_rate_peg_with_contemporaneous_innovation_0_and_rest_FG.pdf",
        PlotT=horizon,
        peg_len=peg_horizon
    )
end

# =========================
# 4) FG exercise:
#    choose rm_sh so that impact policy-rate response is +0.25 (quarterly model units),
#    then peg next 5 quarters with rm_shl1..rm_shl5
#    stars in first panel at quarters 1..6
# =========================
begin
    # First compute the policy-rate response to a -1 shock, then scale to hit +0.25 on impact
    tmp_plotvars = [:obs_nominalrate]
    irf_unit, _ = collect_irfs(
        m, system, s_0,
        tmp_plotvars,
        [:rm_sh],
        horizon;
        shock_scales=Dict(:rm_sh => -1.0)
    )

    unit_policy_impact = irf_unit[1, 1, 1]
    rm_scale = 0.25 / unit_policy_impact

    run_star_fg_exercise(
        m, system, s_0;
        plotvars=plotvars,
        rm_scale=rm_scale,
        pdf_file="interest_rate_peg_with_contemporaneous_innovation_target_plus_025_and_rest_FG.pdf",
        PlotT=horizon,
        peg_len=peg_horizon
    )
end

# =========================
# 5) Permanent liquidity shock peg
#    kept as a 5-quarter peg as well
# =========================
begin
    shock_name   = :b_liqp_sh
    var_name     = :obs_nominalrate
    var_value    = -1.0/4
    desired_path = vcat(fill(var_value, peg_horizon), zeros(horizon - peg_horizon))

    exo     = m.exogenous_shocks
    nshocks = size(system[:RRR], 2)
    nstates = size(system[:TTT], 1)
    nobs    = size(system[:ZZ], 1)
    npseudo = size(system[:ZZ_pseudo], 1)

    states = zeros(nstates, horizon, nshocks)
    obs    = zeros(nobs,    horizon, nshocks)
    pseudo = zeros(npseudo, horizon, nshocks)

    shocks = zeros(nshocks, horizon)

    for t in 1:horizon
        var_value_att = desired_path[t] - obs[m.observables[var_name], t, exo[shock_name]]
        shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(
            var_value_att,
            m.observables[var_name],
            exo[shock_name],
            system[:ZZ],
            system[:RRR]
        )

        states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ =
            forecast(system, s_0, shocks)
    end

    y = shocks[m.exogenous_shocks[:b_liqp_sh], :]
    x = 1:horizon
    idx = findall(abs.(y) .> 1e-10)

    p1 = plot(
        x[idx],
        4 .* y[idx],
        seriestype=:scatter,
        marker=:star5,
        markersize=6,
        markercolor=:blue,
        title="Permanent Liquidity Innovations (APR)",
        label="",
        xticks=my_xticks,
        ylims=panel_ylim(1)
    )
    ylabel!(p1, "%")
    xlabel!(p1, "Quarter")
    plot!(p1, x, zeros(horizon), lc=:black, lw=2, label="")

    p2 = plot(1:horizon, 4 .* obs[m.observables[:obs_nominalrate], :, m.exogenous_shocks[:b_liqp_sh]],
        title="Policy rate (APR)", xticks=my_xticks, ylims=panel_ylim(2))
    ylabel!(p2, "%"); xlabel!(p2, "Quarter"); plot!(p2, zeros(horizon), lc=:black, lw=2, label="")

    p3 = plot(1:horizon, 4 .* obs[m.observables[:obs_corepce], :, m.exogenous_shocks[:b_liqp_sh]],
        title="Inflation (%, qoq annualized)", xticks=my_xticks, ylims=panel_ylim(3))
    ylabel!(p3, "%"); xlabel!(p3, "Quarter"); plot!(p3, zeros(horizon), lc=:black, lw=2, label="")

    p4 = plot(1:horizon, states[m.states[:y_t], :, m.exogenous_shocks[:b_liqp_sh]],
        title="Output (% dev from SS)", xticks=my_xticks, ylims=panel_ylim(4))
    ylabel!(p4, "%"); xlabel!(p4, "Quarter"); plot!(p4, zeros(horizon), lc=:black, lw=2, label="")

    p5 = plot(1:horizon, 4 .* pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate], :, m.exogenous_shocks[:b_liqp_sh]],
        title="r* (Forward 5-year real natural rate, APR)", xticks=my_xticks, ylims=panel_ylim(5))
    ylabel!(p5, "%"); xlabel!(p5, "Quarter"); plot!(p5, zeros(horizon), lc=:black, lw=2, label="")

    p6 = plot(1:horizon, 4 .* pseudo[m.pseudo_observables[:ExAnteRealRate], :, m.exogenous_shocks[:b_liqp_sh]],
        title="Ex-ante real rate (APR)", xticks=my_xticks, ylims=panel_ylim(6))
    ylabel!(p6, "%"); xlabel!(p6, "Quarter"); plot!(p6, zeros(horizon), lc=:black, lw=2, label="")

    p = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(1200,800), legend=false)

    pdf_path = joinpath(saveroot, "Final Paper", "figures", "IRF_rate_peg_with_permanent_liquidity_shock.pdf")
    savefig(p, pdf_path)

    export_series_xlsx(
        xlsx_out,
        sheet_from_pdf(pdf_path),
        tgrid,
        Dict(
            "Permanent liquidity shock (b_liqp_sh)" => vec(shocks[m.exogenous_shocks[:b_liqp_sh], :]),
            "Policy rate (obs_nominalrate)"        => vec(4 .* obs[m.observables[:obs_nominalrate], :, m.exogenous_shocks[:b_liqp_sh]]),
            "Inflation (obs_corepce)"              => vec(4 .* obs[m.observables[:obs_corepce], :, m.exogenous_shocks[:b_liqp_sh]]),
            "Output (obs_gdp)"                     => vec(states[m.states[:y_t], :, m.exogenous_shocks[:b_liqp_sh]]),
            "r* (Forward5YearRealNaturalRate)"     => vec(4 .* pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate], :, m.exogenous_shocks[:b_liqp_sh]]),
            "Ex-ante real rate (ExAnteRealRate)"   => vec(4 .* pseudo[m.pseudo_observables[:ExAnteRealRate], :, m.exogenous_shocks[:b_liqp_sh]]),
            "Zero line"                            => zeros(horizon),
        )
    )
end



#####################
# Permanent liquidity shock #
#####################

m = Model1010("ss20");
system = compute_system(m)
m <= DSGE.Setting(:data_vintage, "250825")
 params_mode = load_draws(m, :mode)
DSGE.update!(m, params_mode)
DSGE.steadystate!(m)
system = compute_system(m)

shock_name = :b_liqp_sh # Select MP to implement the specific path in state variable 
var_name = :obs_nominalrate # Select the targeted state variable
var_value = -1.0/4  # Select the depth of the path
peg_horizon = 5;
desired_path = vcat(fill(var_value, peg_horizon), zeros(horizon - peg_horizon))

# Setup - copied from impulse_responses.jl
    var_names, var_class =
    if var_name in keys(m.endogenous_states)
        m.endogenous_states, :states
    else
        m.observables, :obs
    end
    exo          = m.exogenous_shocks
    nshocks      = size(system[:RRR], 2)
    nstates      = size(system[:TTT], 1)
    nobs         = size(system[:ZZ], 1)
    npseudo      = size(system[:ZZ_pseudo], 1)

    states = zeros(nstates, horizon, nshocks)
    obs    = zeros(nobs,    horizon, nshocks)
    pseudo = zeros(npseudo, horizon, nshocks)

    # Set constant system matrices to 0
    system = DSGE.zero_system_constants(system)

    s_0 = zeros(nstates)

    # Isolate single shock
    shocks = zeros(nshocks, horizon)
    for t = 1:horizon
        if var_class == :states
            var_value_att = desired_path[t] - obs[m.endogenous_states[var_name],t, m.exogenous_shocks[shock_name]]
            shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_state_value(var_value_att,
                                                                        var_names[var_name],
                                                                        exo[shock_name],
                                                                        system[:RRR])
        else # == :obs
            var_value_att = desired_path[t] - obs[m.observables[var_name],t, m.exogenous_shocks[shock_name]]
            shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(var_value_att,
                                                                        var_names[var_name],
                                                                        exo[shock_name],
                                                                        system[:ZZ],
                                                                        system[:RRR])
        end
    
    # Iterate state space forward
    states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)
    end


using Plots
#y = states[m.endogenous_states[:b_liqp_t],:, m.exogenous_shocks[:b_liqp_sh]]
y = shocks[m.exogenous_shocks[:b_liqp_sh],:]
x = 1:horizon

# Keep only non-zero entries
idx = [1,2,3,4,5,6]   # example horizons
#idx = abs.(y) .> 10^-6

p1 = plot(
    x[idx],
    y[idx]*4,
    seriestype = :scatter,
    marker = :star5,
    markersize = 6,
    markercolor = :blue,
    title = "Permanent Liquidity Innovations (APR)",
    label = "", xticks = my_xticks, ylims = panel_ylim(1)./3
)
    ylabel!(p1, "%")
    xlabel!(p1, "Quarter")
plot!(p1, x, zeros(horizon), lc=:black, lw=2, label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Policy rate (APR)", xticks = my_xticks, ylims = panel_ylim(2))#
    ylabel!(p2, "%")
    xlabel!(p2, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_corepce],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Inflation (%, qoq annualized)", xticks = my_xticks, ylims = panel_ylim(3).*6)#
    ylabel!(p3, "%")
    xlabel!(p3, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,states[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Output (% dev from SS)", xticks = my_xticks, ylims = panel_ylim(4).*3)#
    ylabel!(p4, "%")
    xlabel!(p4, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="r* (Forward 5-year real natural rate, APR)", ylims = panel_ylim(5).*15, xticks = my_xticks)#
    ylabel!(p5, "%")
    xlabel!(p5, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Ex-ante real rate (APR)", xticks = my_xticks, ylims = panel_ylim(6))#
    ylabel!(p6, "%")
    xlabel!(p6, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

    p = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(1200,800), legend=false)
pdf_path = joinpath(saveroot,"Final Paper","figures","IRF_rate_peg_with_permanent_liquidity_shock.pdf")
savefig(pdf_path)   # saves the plot from p as a .pdf vector graphic

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    tgrid,
    Dict(
        "Monetary policy shock (rm_t)"            => vec(states[m.endogenous_states[:rm_t],:, m.exogenous_shocks[:rm_sh]]),
        "Policy rate (obs_nominalrate)"          => vec(obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]]),
        "Inflation (obs_corepce)"            => vec(obs[m.observables[:obs_corepce],:, m.exogenous_shocks[:rm_sh]]),
        "Output (obs_gdp)"                       => vec(obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]]),
        "r* (Forward5YearRealNaturalRate)"       => vec(pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh]]),
        "Ex-ante real rate (ExAnteRealRate)"     => vec(pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:rm_sh]]),
        "Zero line"                              => zeros(horizon),
    )
)
# ============================