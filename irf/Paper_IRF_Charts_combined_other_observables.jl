using DSGE
using Plots

# Standalone companion to Paper_IRF_Charts_combined_updated.jl.
# It repeats the 6-quarter policy-rate peg exercise for other observables.

repo_root = normpath(joinpath(@__DIR__, ".."))
dataroot = joinpath(repo_root, "dsge", "input_data")
datafolder = joinpath(repo_root, "dsge", "output_data")
saveroot = joinpath(repo_root, "dsge")
figures_dir = joinpath(saveroot, "Final Paper", "figures")
mkpath(figures_dir)

horizon = 20
peg_horizon = 6
target_policy_rate = -1.0 / 4.0
plot_quarters = 1:horizon
xticks_20q = [1, 4, 8, 12, 16, 20]

plotvars = [
    :obs_nominalrate,
    :obs_nominalrate6,
    :obs_AAAspread,
    :obs_BBBspread,
]

panel_titles = [
    "Policy rate (APR)",
    "OIS 6-quarter forward rate (APR)",
    "AAA spread (APR)",
    "BAA spread (APR)",
]

excel_dark_blue = "#1F497D"
excel_dark_red = "#C0504D"
excel_orange = "#F79646"

function build_ss20_mode_model()
    m = Model1010("ss20")
    m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
    m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
    m <= DSGE.Setting(:data_vintage, "250825")
    params_mode = load_draws(m, :mode)
    DSGE.update!(m, params_mode)
    DSGE.steadystate!(m)
    mode_file = joinpath(datafolder, "m1010", "ss20", "estimate", "raw", "paramsmode_vint=250825.h5")
    specify_mode!(m, mode_file)
    return m
end

function assert_model_symbol(m, var::Symbol)
    if !(haskey(m.observables, var) || haskey(m.endogenous_states, var) || haskey(m.pseudo_observables, var))
        error("Model variable $(var) was not found in observables, states, or pseudo-observables.")
    end
end

function assert_shock_symbol(m, shock::Symbol)
    haskey(m.exogenous_shocks, shock) || error("Model shock $(shock) was not found in exogenous_shocks.")
end

function model_series(m, states, obs, pseudo, var::Symbol)
    if haskey(m.observables, var)
        return vec(obs[m.observables[var], 1:horizon])
    elseif haskey(m.endogenous_states, var)
        return vec(states[m.endogenous_states[var], 1:horizon])
    elseif haskey(m.pseudo_observables, var)
        return vec(pseudo[m.pseudo_observables[var], 1:horizon])
    else
        error("Model variable $(var) was not found in observables, states, or pseudo-observables.")
    end
end

function forecast_with_shocks(system, shocks)
    s_0 = zeros(size(system[:TTT], 1))
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    return (states = states, obs = obs, pseudo = pseudo, shocks = shocks)
end

function path_matrix(m, solution)
    paths = zeros(horizon, length(plotvars))
    for (i, var) in enumerate(plotvars)
        paths[:, i] .= model_series(m, solution.states, solution.obs, solution.pseudo, var)
    end
    return paths
end

apr(paths::AbstractMatrix) = paths .* 4.0

function padded_ylim(values; min_pad::Float64 = 0.05)
    vals = collect(skipmissing(vec(values)))
    vmin = minimum(vals)
    vmax = maximum(vals)
    span = vmax - vmin
    pad = max(min_pad, 0.12 * span)
    return (vmin - pad, vmax + pad)
end

function solve_single_shock_policy_peg(m, system, shock_name::Symbol; solve_periods::Int = peg_horizon)
    assert_shock_symbol(m, shock_name)

    nshocks = size(system[:RRR], 2)
    shocks = zeros(nshocks, horizon)
    exo_idx = m.exogenous_shocks[shock_name]
    policy_idx = m.observables[:obs_nominalrate]
    solution = forecast_with_shocks(system, shocks)

    for t in 1:solve_periods
        remaining_gap = target_policy_rate - solution.obs[policy_idx, t]
        shocks[exo_idx, t] = DSGE.obtain_shock_from_desired_obs_value(
            remaining_gap,
            policy_idx,
            exo_idx,
            system[:ZZ],
            system[:RRR],
        )
        solution = forecast_with_shocks(system, shocks)
    end

    return solution
end

function unit_shock_paths(m, system, shock_names::Vector{Symbol})
    for shock in shock_names
        assert_shock_symbol(m, shock)
    end

    nshocks_total = size(system[:RRR], 2)
    paths = zeros(horizon, length(plotvars), length(shock_names))

    for (j, shock) in enumerate(shock_names)
        shocks = zeros(nshocks_total, horizon)
        shocks[m.exogenous_shocks[shock], 1] = 1.0
        solution = forecast_with_shocks(system, shocks)
        paths[:, :, j] .= path_matrix(m, solution)
    end

    return paths
end

function solve_news_policy_peg(m, system; base_paths = zeros(horizon, length(plotvars)))
    news_shocks = [:rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6]
    irfmat = unit_shock_paths(m, system, news_shocks)

    policy_response = zeros(peg_horizon, peg_horizon)
    for j in 1:peg_horizon
        policy_response[:, j] .= irfmat[1:peg_horizon, 1, j]
    end

    target_path = fill(target_policy_rate, peg_horizon) .- base_paths[1:peg_horizon, 1]
    weights = policy_response \ target_path

    paths = copy(base_paths)
    for j in 1:peg_horizon
        paths .+= weights[j] .* irfmat[:, :, j]
    end

    return (paths = paths, weights = weights)
end

function compute_policy_instrument_paths(m, system)
    surprise_solution = solve_single_shock_policy_peg(m, system, :rm_sh)
    surprise_paths = path_matrix(m, surprise_solution)

    anticipated = solve_news_policy_peg(m, system)

    mp_initial_solution = solve_single_shock_policy_peg(m, system, :rm_sh; solve_periods = 1)
    mp_initial_paths = path_matrix(m, mp_initial_solution)
    mixed = solve_news_policy_peg(m, system; base_paths = mp_initial_paths)

    return [
        (label = "Surprise shocks", color = excel_dark_blue, paths = apr(surprise_paths)),
        (label = "Anticipated shocks", color = excel_dark_red, paths = apr(anticipated.paths)),
        (label = "Mix of surprise and anticipated", color = excel_orange, paths = apr(mixed.paths)),
    ]
end

function compute_liquidity_instrument_paths(m, system)
    permanent_solution = solve_single_shock_policy_peg(m, system, :b_liqp_sh)
    transitory_solution = solve_single_shock_policy_peg(m, system, :b_liqtil_sh)

    return [
        (label = "Permanent liquidity shock", color = excel_dark_red, paths = apr(path_matrix(m, permanent_solution))),
        (label = "Transitory liquidity shock", color = excel_orange, paths = apr(path_matrix(m, transitory_solution))),
    ]
end

function make_2x2_other_observables_figure(series_specs; title_prefix::String = "")
    p = plot(layout = (2, 2), size = (1150, 760))

    for panel in 1:length(plotvars)
        panel_values = Float64[]

        for spec in series_specs
            vals = spec.paths[:, panel]
            append!(panel_values, vals)
            plot!(p[panel], plot_quarters, vals;
                color = spec.color,
                lw = 2.5,
                label = panel == 1 ? spec.label : "",
                legend = panel == 1 ? :topright : false,
                xticks = xticks_20q,
            )
        end

        plot!(p[panel], plot_quarters, zeros(horizon);
            color = :black,
            lw = 1,
            label = "",
            legend = panel == 1 ? :topright : false,
            xticks = xticks_20q,
        )

        title_text = isempty(title_prefix) ? panel_titles[panel] : string(title_prefix, panel_titles[panel])
        title!(p[panel], title_text)
        ylabel!(p[panel], "%")
        xlabel!(p[panel], "Quarter")
        xlims!(p[panel], (1, horizon))
        ylims!(p[panel], padded_ylim(vcat(panel_values, [0.0])))
    end

    plot!(p[1]; legend = :topright)
    return p
end

function save_pdf_and_png(p, basename_no_ext::AbstractString)
    pdf_path = joinpath(figures_dir, string(basename_no_ext, ".pdf"))
    png_path = joinpath(figures_dir, string(basename_no_ext, ".png"))
    savefig(p, pdf_path)
    savefig(p, png_path)
    println("Saved ", pdf_path)
    println("Saved ", png_path)
end

m = build_ss20_mode_model()
system = DSGE.zero_system_constants(compute_system(m))

for var in plotvars
    assert_model_symbol(m, var)
end

policy_series = compute_policy_instrument_paths(m, system)
policy_fig = make_2x2_other_observables_figure(policy_series)
save_pdf_and_png(policy_fig, "interest_rate_peg_combined_other_observables")

liquidity_series = compute_liquidity_instrument_paths(m, system)
liquidity_fig = make_2x2_other_observables_figure(liquidity_series)
save_pdf_and_png(liquidity_fig, "interest_rate_peg_liquidity_other_observables")
