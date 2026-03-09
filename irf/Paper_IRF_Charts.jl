using DSGE;
using Plots # no need for `using Plots` as that is reexported here

# =========================
# Excel export (ADDED ONLY)
# =========================
using XLSX
using Dates
saveroot = dirname(@__FILE__())
path = dirname(@__FILE__)

# Sheet name = PDF basename; return a *String* (not SubString)
function sheet_from_pdf(pdf_path::AbstractString; maxlen::Int=31)::String
    base = basename(pdf_path)
    name = splitext(base)[1]              # might be SubString internally
    name = String(name)                  # FORCE String

    # Excel forbids: : \ / ? * [ ]
    name = replace(name, r"[:\\\/\?\*\[\]]" => "_")
    name = strip(name)
    isempty(name) && (name = "Sheet")

    if lastindex(name) > maxlen
        name = first(name, maxlen)
    end
    return String(name)                  # ensure String again
end

"""
Write/overwrite a sheet in an Excel workbook.

- Uses XLSX.writetable! (stable across XLSX.jl versions).
- If workbook exists: open in rw and write sheet
- If workbook does not exist: create it via openxlsx(mode="w")
"""
function export_series_xlsx(filepath::AbstractString,
                            sheetname_in::AbstractString,
                            t,
                            series::Dict{String,<:AbstractVector})

    mkpath(dirname(filepath))
    sheetname = String(sheetname_in)  # FORCE String

    # Assemble columns (exactly as provided; no computation changes)
    colnames = ["t"; collect(keys(series))]
    cols = Any[collect(t)]
    for nm in colnames[2:end]
        push!(cols, series[nm])
    end

    mode = isfile(filepath) ? "rw" : "w"

    XLSX.openxlsx(filepath, mode=mode) do xf
        # If sheet exists, delete+recreate is version-fragile; instead, just write over A1.
        # Ensure sheet exists:
        if !(sheetname in XLSX.sheetnames(xf))
            XLSX.addsheet!(xf, sheetname)   # now safe because sheetname is String
        end

        sh = xf[sheetname]

        # Clear is not needed; writetable! overwrites starting at A1.
        XLSX.writetable!(sh, cols, colnames; anchor_cell=XLSX.CellRef("A1"))

        # Stamp (optional, but helpful)
        sh[1, length(colnames)+2] = "exported_at"
        sh[1, length(colnames)+3] = string(now())
    end

    return filepath
end

# =========================
# Update model parameters with the mode - note that specify_mode! is needed to update the system matrices for the IRF computations below
# =========================

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "input_data")
saveroot = joinpath(basepath, "dsge")

idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
datafolder = joinpath(basepath, "dsge", "output_data")

horizon  = 20
peg_horizon = 6 # Length of the peg in periods
m = Model1010("ss20");
m <= DSGE.Setting(:data_vintage, "250825")
params_mode = load_draws(m, :mode)
DSGE.update!(m, params_mode)
DSGE.steadystate!(m)
mode_file = joinpath(datafolder, "m1010","ss20","estimate","raw", "paramsmode_vint=250825.h5")
specify_mode!(m, mode_file)

var_name = :obs_nominalrate # Select the targeted state variable
system = compute_system(m)

# Output workbook (ADDED ONLY)
xlsx_out = joinpath(path, "irf", "Paper IRFs.xlsx")
tgrid = collect(1:horizon)

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

    else
        m.observables, :obs,  m.observables[var_name]
    end


    for t in 1:horizon-1
        
        # Compute IRF for a unit shock at time t for the specified shock
        test_shocks = zeros(nshocks, horizon)
        test_shocks[shock_inds[t], t] = 1.0
        states, obs, _ = forecast(system, s_0, test_shocks)
                if var_class == :states
                    irf = states[peg_ind, t] # Impact of a unit shock at t on state at t
                else 
                    irf = obs[peg_ind,t] # Impact of a unit shock at t on obs at t
                end
        # Compute effect of previous shocks
        prev_effect = 0.0
        if t > 1
            prev_shocks = shocks[:, 1:t-1]
            prev_states, prev_obs, _ = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
            if var_class == :states
                prev_effect = prev_states[peg_ind, t]
            else 
                prev_effect = prev_obs[peg_ind,t] # Impact of a unit shock at t on obs at t
            end
        end

        # Required shock at time t to achieve desired value
        shocks[shock_inds[t], t] = (x[t] - prev_effect) / irf
    end

    return shocks
end


#####################   
# Contemporaneous MP innovations #
#####################


shock_name = :rm_sh # Select MP to implement the specific path in state variable 
var_name = :obs_nominalrate # Select the targeted state variable
var_value = -1.0/4  # Select the depth of the path
peg_horizon = 6;

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
    for t = 1:peg_horizon
        if var_class == :states
            var_value_att = var_value - obs[m.endogenous_states[var_name],t, m.exogenous_shocks[shock_name]]
            shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_state_value(var_value_att,
                                                                        var_names[var_name],
                                                                        exo[shock_name],
                                                                        system[:RRR])
        else # == :obs
            var_value_att = var_value - obs[m.observables[var_name],t, m.exogenous_shocks[shock_name]]
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
y = shocks[m.exogenous_shocks[:rm_sh],:]
x = 1:horizon

# Keep only non-zero entries
idx = [1,2,3,4,5,6]   # example horizons
#idx = abs.(y) .> 10^-6
my_xticks = 0:4:20

p1 = plot(
    x[idx],
    y[idx]*4,
    seriestype = :scatter,
    marker = :star5,
    markersize = 6,
    markercolor = :black,
    title = "Contemporanous Policy Innovations (APR)",
    label = "",
    xticks = my_xticks

)
    ylabel!(p1, "%")
    xlabel!(p1, "Quarter")
plot!(p1, x, zeros(horizon), lc=:black, lw=2, label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]]*4,title="Policy rate (APR)", xticks = my_xticks)
    ylabel!(p2, "%")
    xlabel!(p2, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_corepce],:, m.exogenous_shocks[:rm_sh]]*4,title="Inflation (%, qoq annualized)",    xticks = my_xticks)
    ylabel!(p3, "%")
    xlabel!(p3, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,states[m.endogenous_states[:y_t],:, m.exogenous_shocks[:rm_sh]],title="Output (% dev from SS)",    xticks = my_xticks)
    ylabel!(p4, "%")
    xlabel!(p4, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh]]*4,title="r* (Forward 5-year real natural rate, APR)", ylims = (-0.1, 0.1),    xticks = my_xticks)
    ylabel!(p5, "%")
    xlabel!(p5, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:rm_sh]]*4,title="Ex-ante real rate (APR)",    xticks = my_xticks)
    ylabel!(p6, "%")
    xlabel!(p6, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,p5,p6,layout=(3,2), legend=false)
plot!(layout=(3,2), size=(1200,800))
pdf_path = joinpath(saveroot,"Final Paper","figures","Interest_rate_peg_IRF_Policy_rate_with_MP_shock.pdf")
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
# ================================================


# ########################################################################
# Proper forward guidance implementation by Matyas Farkas, IMF 19/08/2025 based on the Matlab code from Jesper Linde and Zoltan Jakab
# ########################################################################

# using LinearAlgebra
# using Plots

# # --- Step 1: Compute IRFs for each shock ---
PlotT = horizon # Total IRF horizon to plot
# plotvars = [:y_t, :obs_corepce, :obs_nominalrate] # Output, Inflation, Policy Rate
# shock_syms = [:rm_sh, :rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6] # MP + 1-6 FG shocks

# nvars = length(plotvars)
# nshocks = length(shock_syms)

# # Store IRFs: irfmat[t, var, shock]
# irfmat = zeros(PlotT, nvars, nshocks)
# for (j, shock_sym) in enumerate(shock_syms)
#     shocks = zeros(size(system[:RRR], 2), PlotT)
#     shocks[m.exogenous_shocks[shock_sym], 1] = 1.0
#     states, obs, pseudo, _ = forecast(system, s_0, shocks)
#     for (i, var_sym) in enumerate(plotvars)
#         if var_sym in keys(m.observables)
#             irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
#         elseif var_sym in keys(m.endogenous_states)
#             irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
#         elseif var_sym in keys(m.pseudo_observables)
#             irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
#         end
#     end
# end

# # --- Step 2: Loop over FG horizons and compute weights and IRFs ---
# FGhorz = 1:(peg_horizon-1) # Try 1 to 5 horizon pegs
# FGplotmat = zeros(PlotT, nvars, length(FGhorz))
# shk_weights_store = zeros(PlotT, length(FGhorz))

# for (hidx, horz) in enumerate(FGhorz)
#     FGdur = horz + 1 # Duration of FG in periods (Matlab uses +1)
#     FG_vec = fill(-1.0, FGdur)/4 # Desired policy rate path

#     # Build IRF matrix for policy rate
#     R_mp_mat = zeros(FGdur, FGdur)
#     for s = 1:FGdur
#         R_mp_mat[:, s] .= irfmat[1:FGdur, 3, s] # 3rd var is policy rate
#     end

#     # Solve for shock weights
#     shk_weights = R_mp_mat \ FG_vec
#     shk_weights_store[1:FGdur, hidx] .= shk_weights

#     # Construct total IRFs for each variable
#     for s = 1:FGdur
#         FGplotmat[:, :, hidx] .+= shk_weights[s] .* irfmat[:, :, s]
#     end
# end


# # Propagate for full PlotT
# # --- Step 3: Plot results ---
# titles = ["Output", "Inflation", "Nominal Policy Rate"]
# TT = 1:PlotT

# p = plot(layout=(1,3), size=(1200,400))
# for i = 1:nvars
#     plot!(p[i], TT, FGplotmat[:, i, peg_horizon-1], lw=2, label="FG IRF (horizon=$peg_horizon)")
#     plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
#     title!(p[i], titles[i])
#     ylabel!(p[i], "Percent")
#     xlabel!(p[i], "Quarter")
# end
# plot!(p)
# pdf_path = joinpath(saveroot,"Final Paper","figures","Vanila_FG_6horizon_policy_rate_output_inflation.pdf")
# savefig(pdf_path)

# # ===== ADDED ONLY: export underlying plotted data =====
# export_series_xlsx(
#     xlsx_out,
#     sheet_from_pdf(pdf_path),
#     collect(1:PlotT),
#     Dict(
#         "Output"       => vec(FGplotmat[:, 1, peg_horizon-1]),
#         "Inflation"    => vec(FGplotmat[:, 2, peg_horizon-1]*4),
#         "Policy Rate"  => vec(FGplotmat[:, 3, peg_horizon-1]*4),
#         "Zero line"    => zeros(PlotT),
#     )
# )
# # ================================================

# # Optional: plot weights for each horizon
# pw = plot(1:PlotT, shk_weights_store[:, peg_horizon-1], lw=2, label="Shock weights")
# xlabel!("Quarter ahead shocks")
# ylabel!("Weight")
# title!("Shock Weights for FG horizon $peg_horizon")
# pdf_path = "FG_6horizon_shock_weights.pdf"
# savefig(pdf_path)

# # ===== ADDED ONLY: export underlying plotted data =====
# export_series_xlsx(
#     xlsx_out,
#     sheet_from_pdf(pdf_path),
#     collect(1:PlotT),
#     Dict(
#         "Shock weights" => vec(shk_weights_store[:, peg_horizon-1]),
#         "Zero line"     => zeros(PlotT),
#     )
# )
# # ================================================

#########################################################################
# Adding the other variables to plot
#########################################################################
horizon = 20
plotvars = [ :obs_nominalrate,:obs_corepce,  :obs_gdp, :Forward5YearRealNaturalRate, :ExAnteRealRate] 

titles = ["Combined monetary policy shocks","Policy rate", "Inflation", "Output", "r* (Forward 5-year real natural rate)", "Ex-ante real rate"]
nvars = length(plotvars)


shock_syms = [:rm_sh, :rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6] # MP + 1-6 FG shocks

nvars = length(plotvars)
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)
for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1.0
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    for (i, var_sym) in enumerate(plotvars)
        if var_sym in keys(m.observables)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif var_sym in keys(m.endogenous_states)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif var_sym in keys(m.pseudo_observables)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end


# --- Step 2: Loop over FG horizons and compute weights and IRFs ---
FGhorz = 1:(peg_horizon-1) # Try 1 to 5 horizon pegs
FGplotmat = zeros(PlotT, nvars, length(FGhorz))
shk_weights_store = zeros(PlotT, length(FGhorz))

for (hidx, horz) in enumerate(FGhorz)
    FGdur = horz + 1 # Duration of FG in periods (Matlab uses +1)
    FG_vec = fill(-0.25, FGdur) # Desired policy rate path

    # Build IRF matrix for policy rate
    R_mp_mat = zeros(FGdur, FGdur)
    for s = 1:FGdur
        R_mp_mat[:, s] .= irfmat[1:FGdur, 1, s] # 3rd var is policy rate
    end

    # Solve for shock weights
    shk_weights = R_mp_mat \ FG_vec
    shk_weights_store[1:FGdur, hidx] .= shk_weights

    # Construct total IRFs for each variable
    for s = 1:FGdur
        FGplotmat[:, :, hidx] .+= shk_weights[s] .* irfmat[:, :, s]
    end
end

# --- Step 3: Plot results ---
TT = 1:PlotT

p = plot(layout=(3,2), size=(1200,800))
for i = 1:nvars+1
    if i == 1
        plot!(p[i], TT, shk_weights_store[:, peg_horizon-1]*4, lw=2, label="", xticks=my_xticks)
    elseif i ==34
        if plotvars[i-1] in keys(m.observables)
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", xticks=my_xticks)
        end
    else
        if plotvars[i-1] in keys(m.observables)
            plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1]*4, lw=2, label="", xticks=my_xticks)
        elseif plotvars[i-1] in keys(m.pseudo_observables)
            if plotvars[i-1] == :Forward5YearRealNaturalRate || plotvars[i-1] == :RealNaturalRate
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", ylims=(-0.1, 0.1), xticks=my_xticks)
            elseif plotvars[i-1] == :ExAnteRealRate
                 plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1]*4, lw=2, label="", xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", xticks=my_xticks)
            end
        end
    end
    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "Percent")
    xlabel!(p[i], "Quarter")
end
plot!(p)
# pdf_path = joinpath(saveroot,"Final Paper","figures","Vanilla_interest_rate_peg_with_MP_and_FG_stars.pdf")
# savefig(pdf_path)

# ===== ADDED ONLY: export underlying plotted data =====
# export_series_xlsx(
#     xlsx_out,
#     sheet_from_pdf(pdf_path),
#     collect(1:PlotT),
#     Dict(
#         "Combined monetary policy shocks" => vec(shk_weights_store[:, peg_horizon-1]),
#         "Policy rate"                     => vec(FGplotmat[:, 1, peg_horizon-1]),
#         "Inflation"                       => vec(FGplotmat[:, 2, peg_horizon-1]),
#         "Output"                          => vec(FGplotmat[:, 3, peg_horizon-1]),
#         "r* (Forward 5-year real natural rate)" => vec(FGplotmat[:, 4, peg_horizon-1]),
#         "Ex-ante real rate"               => vec(FGplotmat[:, 5, peg_horizon-1]),
#         "Zero line"                       => zeros(PlotT),
#     )
# )
# ================================================


#########################################################################
# FG exercise: "rm_sh = 1 (1.738std) at t=0" + peg for periods 1..h with FG shocks only (rm_shl1..rm_shl6)
#########################################################################

# --- User knobs ---
PlotT      = horizon                  # IRF length

plotvars = [:obs_nominalrate, :obs_corepce, :obs_gdp,
            :Forward5YearRealNaturalRate, :ExAnteRealRate]

titles = ["Anticipated Policy Innovations (APR)",
          "Policy rate (APR)",
          "Inflation (%, qoq annualized)",
          "Output (% dev from SS)",
          "r* (Forward 5-year real natural rate, APR)",
          "Ex-ante real rate (APR)"]

nvars = length(plotvars)

# Include contemporaneous MP shock at t=0 plus 1..6 FG (news) shocks
shock_syms = [:rm_sh]
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)

for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = -1 #     shocks[m.exogenous_shocks[shock_sym], 1] = -1
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1/obs[m.observables[:obs_nominalrate], 1] # Normalize to 1 unit shock on policy rate
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    for (i, var_sym) in enumerate(plotvars)
        if haskey(m.observables, var_sym)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif haskey(m.endogenous_states, var_sym)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif haskey(m.pseudo_observables, var_sym)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end

irfmatR = irfmat # IRF of policy rate to contemporaneous MP shock (used for peg target and FG weights)

plotvars = [ :obs_nominalrate,:obs_corepce,  :obs_gdp, :Forward5YearRealNaturalRate, :ExAnteRealRate] 

titles = ["Anticipated Policy Innovations (APR)","Policy rate (APR)", "Inflation (%, qoq annualized)", "Output (% dev from SS)", "r* (Forward 5-year real natural rate, APR)", "Ex-ante real rate (APR)"]
nvars = length(plotvars)


shock_syms = [:rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6] #  1-6 FG shocks

nvars = length(plotvars)
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)
for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1.0
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    for (i, var_sym) in enumerate(plotvars)
        if var_sym in keys(m.observables)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif var_sym in keys(m.endogenous_states)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif var_sym in keys(m.pseudo_observables)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end


# --- Step 2: Loop over FG horizons and compute weights and IRFs ---
FGhorz = 1:(peg_horizon-1) # Try 1 to 5 horizon pegs
FGplotmat = zeros(PlotT, nvars, length(FGhorz))
shk_weights_store = zeros(PlotT, length(FGhorz))

for (hidx, horz) in enumerate(FGhorz)
    FGdur = horz + 1 # Duration of FG in periods (Matlab uses +1)
    FG_vec = (fill(-1.0, FGdur)-irfmatR[1:FGdur, 1, 1]) # Desired policy rate path

    # Build IRF matrix for policy rate
    R_mp_mat = zeros(FGdur, FGdur)
    for s = 1:FGdur
        R_mp_mat[:, s] .= irfmat[1:FGdur, 1, s] # 3rd var is policy rate
    end

    # Solve for shock weights
    shk_weights = R_mp_mat \ FG_vec
    shk_weights_store[1:FGdur, hidx] .= shk_weights

    # Construct total IRFs for each variable
    for s = 1:FGdur
        FGplotmat[:, :, hidx] .+= shk_weights[s] .* irfmat[:, :, s]
    end
end
FGplotmat[:, 1, peg_horizon-1] .+= irfmatR[:, 1, 1] # Add contemporaneous MP shock back to policy rate IRF  

# --- Step 3: Plot results ---
TT = 1:PlotT

p = plot(layout=(3,2), size=(1200,800))

for i = 1:nvars+1
    if i == 1
        # shk_weights_store[t, h] = weight on the t-th shock (MP/news shock at that timing)
        # for the forward-guidance implementation of length FGdur = h+1 (i.e., horizon h uses 2..(h+1) shocks).
        # Visualization: "triangular build-up" — at t=1 plot 2 stars, at t=2 plot 3, ..., at t=5 plot 6.
        Tmax = 5 # min(5, PlotT)                                  # show periods 1..5
        Hmax = 6             # up to 6 horizons/columns (=> up to 6 stars)
        alphas = collect(range(1.0, 0.15, length=Hmax))        # longer horizon => more transparent

        for t in Tmax #1:Tmax
            nh = Hmax #min(t + 1, Hmax)                              # t=1 -> 2 stars, ..., t=5 -> 6 stars
            for hidx in 1:nh
                y = shk_weights_store[hidx, t]
                # (Optional) suppress numerical zeros:
                # if abs(y) <= 1e-12; continue; end
                if hidx == 1
                    plot!(p[i], [hidx], [irfmatR[1, 1, 1]];
                        seriestype = :scatter,
                        marker = :star5,
                        markersize = 7,
                        markercolor = RGBA(0.0,0.0,0.0,1.0), # fully opaque red for contemporaneous shock
                        markerstrokecolor = RGBA(0.0,0.0,0.0,1.0),
                        label = "",
                        xticks=my_xticks
                    )
                end
                plot!(p[i], [hidx+1], [y];
                    seriestype = :scatter,
                    marker = :star5,
                    markersize = 7,
                    markercolor = RGBA(1.0, 0.0, 0.0,1.0), # fully transparent red for marker color (invisible star)
                    markerstrokecolor = RGBA(1.0, 0.0, 0.0, 1.0),
                    label = "",
                    xticks=my_xticks
                )
            end
        end

    else
        if plotvars[i-1] in keys(m.observables)
            if i == 2 # Policy rate: add contemporaneous MP shock to FG response
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",ylims=(-1.5, 1.5), xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        elseif plotvars[i-1] in keys(m.pseudo_observables)
            if plotvars[i-1] == :Forward5YearRealNaturalRate || plotvars[i-1] == :RealNaturalRate
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", ylims=(-0.1, 0.1),xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        end
    end

    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "%")
    xlabel!(p[i], "Quarter")
end

plot!(p)

pdf_path = joinpath(saveroot,"Final Paper","figures","interest_rate_peg_with_contemporaneous_innovation_1_and_rest_FG.pdf")
savefig(pdf_path)

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    collect(1:PlotT),
    Dict(
        "Combined monetary policy shocks" => vec(shk_weights_store[:, peg_horizon-1]),
        "Policy rate"                     => vec(FGplotmat[:, 1, peg_horizon-1]),
        "Inflation"                       => vec(FGplotmat[:, 2, peg_horizon-1]),
        "Output"                          => vec(FGplotmat[:, 3, peg_horizon-1]),
        "r* (Forward 5-year real natural rate)" => vec(FGplotmat[:, 4, peg_horizon-1]),
        "Ex-ante real rate"               => vec(FGplotmat[:, 5, peg_horizon-1]),
        "Zero line"                       => zeros(PlotT),
    )
)
# ================================================


#########################################################################
# FG exercise: "rm_sh = 0 at t=0" + peg for periods 1..h with FG shocks only (rm_shl1..rm_shl6)
#########################################################################

# --- User knobs ---
PlotT      = horizon                  # IRF length

plotvars = [:obs_nominalrate, :obs_corepce, :obs_gdp,
            :Forward5YearRealNaturalRate, :ExAnteRealRate]

titles = ["Anticipated Policy Innovations (APR)",
          "Policy rate (APR)",
          "Inflation (%, yoy)",
          "Output (% dev from SS)",
          "r* (Forward 5-year real natural rate, APR)",
          "Ex-ante real rate (APR)"]

nvars = length(plotvars)

# Include contemporaneous MP shock at t=0 plus 1..6 FG (news) shocks
shock_syms = [:rm_sh]
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)

for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = 0.0
    states, obs, pseudo, _ = forecast(system, s_0, shocks)

    for (i, var_sym) in enumerate(plotvars)
        if haskey(m.observables, var_sym)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif haskey(m.endogenous_states, var_sym)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif haskey(m.pseudo_observables, var_sym)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end

irfmatR = irfmat # IRF of policy rate to contemporaneous MP shock (used for peg target and FG weights)

plotvars = [ :obs_nominalrate,:obs_corepce,  :obs_gdp, :Forward5YearRealNaturalRate, :ExAnteRealRate] 

titles = ["Anticipated Policy Innovations (APR)","Policy rate (APR)", "Inflation (%, qoq annualized)", "Output (% dev from SS)", "r* (Forward 5-year real natural rate, APR)", "Ex-ante real rate (APR)"]
nvars = length(plotvars)


shock_syms = [:rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6] #  1-6 FG shocks

nvars = length(plotvars)
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)
for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1.0
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1/obs[m.observables[:obs_nominalrate], 1] # Normalize to 1 unit shock on policy rate
    states, obs, pseudo, _ = forecast(system, s_0, shocks)

    for (i, var_sym) in enumerate(plotvars)
        if var_sym in keys(m.observables)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif var_sym in keys(m.endogenous_states)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif var_sym in keys(m.pseudo_observables)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end


# --- Step 2: Loop over FG horizons and compute weights and IRFs ---
FGhorz = 1:(peg_horizon-1) # Try 1 to 5 horizon pegs
FGplotmat = zeros(PlotT, nvars, length(FGhorz))
shk_weights_store = zeros(PlotT, length(FGhorz))

for (hidx, horz) in enumerate(FGhorz)
    FGdur = horz + 1 # Duration of FG in periods (Matlab uses +1)
    FG_vec = (fill(-1.0, FGdur)-irfmatR[1:FGdur, 1, 1]) # Desired policy rate path

    # Build IRF matrix for policy rate
    R_mp_mat = zeros(FGdur, FGdur)
    for s = 1:FGdur
        R_mp_mat[:, s] .= irfmat[1:FGdur, 1, s] # 3rd var is policy rate
    end

    # Solve for shock weights
    shk_weights = R_mp_mat \ FG_vec
    shk_weights_store[1:FGdur, hidx] .= shk_weights

    # Construct total IRFs for each variable
    for s = 1:FGdur
        FGplotmat[:, :, hidx] .+= shk_weights[s] .* irfmat[:, :, s]
    end
end
FGplotmat[:, 1, peg_horizon-1] .+= irfmatR[:, 1, 1] # Add contemporaneous MP shock back to policy rate IRF  

# --- Step 3: Plot results ---
TT = 1:PlotT

p = plot(layout=(3,2), size=(1200,800))

for i = 1:nvars+1
    if i == 1
        # shk_weights_store[t, h] = weight on the t-th shock (MP/news shock at that timing)
        # for the forward-guidance implementation of length FGdur = h+1 (i.e., horizon h uses 2..(h+1) shocks).
        # Visualization: "triangular build-up" — at t=1 plot 2 stars, at t=2 plot 3, ..., at t=5 plot 6.
        Tmax = 5 # min(5, PlotT)                                  # show periods 1..5
        Hmax = 6             # up to 6 horizons/columns (=> up to 6 stars)
        alphas = collect(range(1.0, 0.15, length=Hmax))        # longer horizon => more transparent

        for t in Tmax #1:Tmax
            nh = Hmax #min(t + 1, Hmax)                              # t=1 -> 2 stars, ..., t=5 -> 6 stars
            for hidx in 1:nh
                y = shk_weights_store[hidx, t]
                # (Optional) suppress numerical zeros:
                # if abs(y) <= 1e-12; continue; end
                if hidx == 1
                    plot!(p[i], [hidx], [irfmatR[1, 1, 1]];
                        seriestype = :scatter,
                        marker = :star5,
                        markersize = 7,
                        markercolor = RGBA(0.0,0.0,0.0,1.0), # fully opaque red for contemporaneous shock
                        markerstrokecolor = RGBA(0.0,0.0,0.0,1.0),
                        label = "",
                        xticks=my_xticks
                    )
                end
                plot!(p[i], [hidx+1], [y];
                    seriestype = :scatter,
                    marker = :star5,
                    markersize = 7,
                    markercolor = RGBA(1.0, 0.0, 0.0,1.0), # fully transparent red for marker color (invisible star)
                    markerstrokecolor = RGBA(1.0, 0.0, 0.0, 1.0),
                    label = "",
                    xticks=my_xticks
                )
            end
        end

    else
        if plotvars[i-1] in keys(m.observables)
            if i == 2 # Policy rate: add contemporaneous MP shock to FG response
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",ylims=(-1.5, 1.5), xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        elseif plotvars[i-1] in keys(m.pseudo_observables)
            if plotvars[i-1] == :Forward5YearRealNaturalRate || plotvars[i-1] == :RealNaturalRate
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", ylims=(-0.1, 0.1),xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        end
    end

    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "%")
    xlabel!(p[i], "Quarter")
end

plot!(p)

pdf_path = joinpath(saveroot,"Final Paper","figures","interest_rate_peg_with_contemporaneous_innovation_0_and_rest_FG.pdf")
savefig(pdf_path)

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    collect(1:PlotT),
    Dict(
        "Combined monetary policy shocks" => vec(shk_weights_store[:, peg_horizon-1]),
        "Policy rate"                     => vec(FGplotmat[:, 1, peg_horizon-1]),
        "Inflation"                       => vec(FGplotmat[:, 2, peg_horizon-1]),
        "Output"                          => vec(FGplotmat[:, 3, peg_horizon-1]),
        "r* (Forward 5-year real natural rate)" => vec(FGplotmat[:, 4, peg_horizon-1]),
        "Ex-ante real rate"               => vec(FGplotmat[:, 5, peg_horizon-1]),
        "Zero line"                       => zeros(PlotT),
    )
)
# ================================================



#########################################################################
# FG exercise: "rm_sh = 1.738 at t=0" that implements period 0 policy rate to be 0.25 and peg for periods 1..h with FG shocks only (rm_shl1..rm_shl6)
#########################################################################

# --- User knobs ---
PlotT      = horizon                  # IRF length

plotvars = [:obs_nominalrate, :obs_corepce, :obs_gdp,
            :Forward5YearRealNaturalRate, :ExAnteRealRate]

titles = ["Anticipated Policy Innovations (APR)",
          "Policy rate (APR)",
          "Inflation (%, qoq annualized)",
          "Output (% dev from SS)",
          "r* (Forward 5-year real natural rate, APR)",
          "Ex-ante real rate (APR)"]

nvars = length(plotvars)

# Include contemporaneous MP shock at t=0 plus 1..6 FG (news) shocks
shock_syms = [:rm_sh]
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)

for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = -1 #     shocks[m.exogenous_shocks[shock_sym], 1] = -1
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1/obs[m.observables[:obs_nominalrate], 1] # Normalize to 1 unit shock on policy rate
    states, obs, pseudo, _ = forecast(system, s_0, shocks)

    for (i, var_sym) in enumerate(plotvars)
        if haskey(m.observables, var_sym)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif haskey(m.endogenous_states, var_sym)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif haskey(m.pseudo_observables, var_sym)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end

irfmatR = irfmat # IRF of policy rate to contemporaneous MP shock (used for peg target and FG weights)

plotvars = [ :obs_nominalrate,:obs_corepce,  :obs_gdp, :Forward5YearRealNaturalRate, :ExAnteRealRate] 

titles = ["Anticipated Policy Innovations (APR)","Policy rate (APR)", "Inflation (%, qoq annualized)", "Output (% dev from SS)", "r* (Forward 5-year real natural rate, APR)", "Ex-ante real rate (APR)"]
nvars = length(plotvars)


shock_syms = [:rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6] #  1-6 FG shocks

nvars = length(plotvars)
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)
for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1.0
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    for (i, var_sym) in enumerate(plotvars)
        if var_sym in keys(m.observables)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif var_sym in keys(m.endogenous_states)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif var_sym in keys(m.pseudo_observables)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end


# --- Step 2: Loop over FG horizons and compute weights and IRFs ---
FGhorz = 1:(peg_horizon-1) # Try 1 to 5 horizon pegs
FGplotmat = zeros(PlotT, nvars, length(FGhorz))
shk_weights_store = zeros(PlotT, length(FGhorz))

for (hidx, horz) in enumerate(FGhorz)
    FGdur = horz + 1 # Duration of FG in periods (Matlab uses +1)
    FG_vec = (fill(-1.0, FGdur)-irfmatR[1:FGdur, 1, 1]) # Desired policy rate path

    # Build IRF matrix for policy rate
    R_mp_mat = zeros(FGdur, FGdur)
    for s = 1:FGdur
        R_mp_mat[:, s] .= irfmat[1:FGdur, 1, s] # 3rd var is policy rate
    end

    # Solve for shock weights
    shk_weights = R_mp_mat \ FG_vec
    shk_weights_store[1:FGdur, hidx] .= shk_weights

    # Construct total IRFs for each variable
    for s = 1:FGdur
        FGplotmat[:, :, hidx] .+= shk_weights[s] .* irfmat[:, :, s]
    end
end
FGplotmat[:, 1, peg_horizon-1] .+= irfmatR[:, 1, 1] # Add contemporaneous MP shock back to policy rate IRF  

# --- Step 3: Plot results ---
TT = 1:PlotT

p = plot(layout=(3,2), size=(1200,800))

for i = 1:nvars+1
    if i == 1
        # shk_weights_store[t, h] = weight on the t-th shock (MP/news shock at that timing)
        # for the forward-guidance implementation of length FGdur = h+1 (i.e., horizon h uses 2..(h+1) shocks).
        # Visualization: "triangular build-up" — at t=1 plot 2 stars, at t=2 plot 3, ..., at t=5 plot 6.
        Tmax = 5 # min(5, PlotT)                                  # show periods 1..5
        Hmax = 6             # up to 6 horizons/columns (=> up to 6 stars)
        alphas = collect(range(1.0, 0.15, length=Hmax))        # longer horizon => more transparent

        for t in Tmax #1:Tmax
            nh = Hmax #min(t + 1, Hmax)                              # t=1 -> 2 stars, ..., t=5 -> 6 stars
            for hidx in 1:nh
                y = shk_weights_store[hidx, t]
                # (Optional) suppress numerical zeros:
                # if abs(y) <= 1e-12; continue; end
                if hidx == 1
                    plot!(p[i], [hidx], [irfmatR[1, 1, 1]];
                        seriestype = :scatter,
                        marker = :star5,
                        markersize = 7,
                        markercolor = RGBA(0.0,0.0,0.0,1.0), # fully opaque red for contemporaneous shock
                        markerstrokecolor = RGBA(0.0,0.0,0.0,1.0),
                        label = "",
                        xticks=my_xticks
                    )
                end
                plot!(p[i], [hidx+1], [y];
                    seriestype = :scatter,
                    marker = :star5,
                    markersize = 7,
                    markercolor = RGBA(1.0, 0.0, 0.0,1.0), # fully transparent red for marker color (invisible star)
                    markerstrokecolor = RGBA(1.0, 0.0, 0.0, 1.0),
                    label = "",
                    xticks=my_xticks
                )
            end
        end

    else
        if plotvars[i-1] in keys(m.observables)
            if i == 2 # Policy rate: add contemporaneous MP shock to FG response
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",ylims=(-1.5, 1.5), xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        elseif plotvars[i-1] in keys(m.pseudo_observables)
            if plotvars[i-1] == :Forward5YearRealNaturalRate || plotvars[i-1] == :RealNaturalRate
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", ylims=(-0.1, 0.1),xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        end
    end

    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "%")
    xlabel!(p[i], "Quarter")
end

plot!(p)

pdf_path = joinpath(saveroot,"Final Paper","figures","interest_rate_peg_with_contemporaneous_innovation_1_and_rest_FG.pdf")
savefig(pdf_path)

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    collect(1:PlotT),
    Dict(
        "Combined monetary policy shocks" => vec(shk_weights_store[:, peg_horizon-1]),
        "Policy rate"                     => vec(FGplotmat[:, 1, peg_horizon-1]),
        "Inflation"                       => vec(FGplotmat[:, 2, peg_horizon-1]),
        "Output"                          => vec(FGplotmat[:, 3, peg_horizon-1]),
        "r* (Forward 5-year real natural rate)" => vec(FGplotmat[:, 4, peg_horizon-1]),
        "Ex-ante real rate"               => vec(FGplotmat[:, 5, peg_horizon-1]),
        "Zero line"                       => zeros(PlotT),
    )
)
# ================================================


#########################################################################
# FG exercise: "rm_sh = 0 at t=0" + peg for periods 1..h with FG shocks only (rm_shl1..rm_shl6)
#########################################################################

# --- User knobs ---
PlotT      = horizon                  # IRF length

plotvars = [:obs_nominalrate, :obs_corepce, :obs_gdp,
            :Forward5YearRealNaturalRate, :ExAnteRealRate]

titles = ["Anticipated Policy Innovations (APR)",
          "Policy rate (APR)",
          "Inflation (%, yoy)",
          "Output (% dev from SS)",
          "r* (Forward 5-year real natural rate, APR)",
          "Ex-ante real rate (APR)"]

nvars = length(plotvars)

# Include contemporaneous MP shock at t=0 plus 1..6 FG (news) shocks
shock_syms = [:rm_sh]
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)

for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = 0.0
    states, obs, pseudo, _ = forecast(system, s_0, shocks)

    for (i, var_sym) in enumerate(plotvars)
        if haskey(m.observables, var_sym)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif haskey(m.endogenous_states, var_sym)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif haskey(m.pseudo_observables, var_sym)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end

irfmatR = irfmat # IRF of policy rate to contemporaneous MP shock (used for peg target and FG weights)

plotvars = [ :obs_nominalrate,:obs_corepce,  :obs_gdp, :Forward5YearRealNaturalRate, :ExAnteRealRate] 

titles = ["Anticipated Policy Innovations (APR)","Policy rate (APR)", "Inflation (%, qoq annualized)", "Output (% dev from SS)", "r* (Forward 5-year real natural rate, APR)", "Ex-ante real rate (APR)"]
nvars = length(plotvars)


shock_syms = [:rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6] #  1-6 FG shocks

nvars = length(plotvars)
nshocks = length(shock_syms)

# Store IRFs: irfmat[t, var, shock]
irfmat = zeros(PlotT, nvars, nshocks)
for (j, shock_sym) in enumerate(shock_syms)
    shocks = zeros(size(system[:RRR], 2), PlotT)
    shocks[m.exogenous_shocks[shock_sym], 1] = 1.0
    states, obs, pseudo, _ = forecast(system, s_0, shocks)
    for (i, var_sym) in enumerate(plotvars)
        if var_sym in keys(m.observables)
            irfmat[:, i, j] .= obs[m.observables[var_sym], 1:PlotT]
        elseif var_sym in keys(m.endogenous_states)
            irfmat[:, i, j] .= states[m.endogenous_states[var_sym], 1:PlotT]
        elseif var_sym in keys(m.pseudo_observables)
            irfmat[:, i, j] .= pseudo[m.pseudo_observables[var_sym], 1:PlotT]
        end
    end
end


# --- Step 2: Loop over FG horizons and compute weights and IRFs ---
FGhorz = 1:(peg_horizon-1) # Try 1 to 5 horizon pegs
FGplotmat = zeros(PlotT, nvars, length(FGhorz))
shk_weights_store = zeros(PlotT, length(FGhorz))

for (hidx, horz) in enumerate(FGhorz)
    FGdur = horz + 1 # Duration of FG in periods (Matlab uses +1)
    FG_vec = (fill(-1.0, FGdur)-irfmatR[1:FGdur, 1, 1]) # Desired policy rate path

    # Build IRF matrix for policy rate
    R_mp_mat = zeros(FGdur, FGdur)
    for s = 1:FGdur
        R_mp_mat[:, s] .= irfmat[1:FGdur, 1, s] # 3rd var is policy rate
    end

    # Solve for shock weights
    shk_weights = R_mp_mat \ FG_vec
    shk_weights_store[1:FGdur, hidx] .= shk_weights

    # Construct total IRFs for each variable
    for s = 1:FGdur
        FGplotmat[:, :, hidx] .+= shk_weights[s] .* irfmat[:, :, s]
    end
end
FGplotmat[:, 1, peg_horizon-1] .+= irfmatR[:, 1, 1] # Add contemporaneous MP shock back to policy rate IRF  

# --- Step 3: Plot results ---
TT = 1:PlotT

p = plot(layout=(3,2), size=(1200,800))

for i = 1:nvars+1
    if i == 1
        # shk_weights_store[t, h] = weight on the t-th shock (MP/news shock at that timing)
        # for the forward-guidance implementation of length FGdur = h+1 (i.e., horizon h uses 2..(h+1) shocks).
        # Visualization: "triangular build-up" — at t=1 plot 2 stars, at t=2 plot 3, ..., at t=5 plot 6.
        Tmax = 5 # min(5, PlotT)                                  # show periods 1..5
        Hmax = 6             # up to 6 horizons/columns (=> up to 6 stars)
        alphas = collect(range(1.0, 0.15, length=Hmax))        # longer horizon => more transparent

        for t in Tmax #1:Tmax
            nh = Hmax #min(t + 1, Hmax)                              # t=1 -> 2 stars, ..., t=5 -> 6 stars
            for hidx in 1:nh
                y = shk_weights_store[hidx, t]
                # (Optional) suppress numerical zeros:
                # if abs(y) <= 1e-12; continue; end
                if hidx == 1
                    plot!(p[i], [hidx], [irfmatR[1, 1, 1]];
                        seriestype = :scatter,
                        marker = :star5,
                        markersize = 7,
                        markercolor = RGBA(0.0,0.0,0.0,1.0), # fully opaque red for contemporaneous shock
                        markerstrokecolor = RGBA(0.0,0.0,0.0,1.0),
                        label = "",
                        xticks=my_xticks
                    )
                end
                plot!(p[i], [hidx+1], [y];
                    seriestype = :scatter,
                    marker = :star5,
                    markersize = 7,
                    markercolor = RGBA(1.0, 0.0, 0.0,1.0), # fully transparent red for marker color (invisible star)
                    markerstrokecolor = RGBA(1.0, 0.0, 0.0, 1.0),
                    label = "",
                    xticks=my_xticks
                )
            end
        end

    else
        if plotvars[i-1] in keys(m.observables)
            if i == 2 # Policy rate: add contemporaneous MP shock to FG response
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",ylims=(-1.5, 1.5), xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        elseif plotvars[i-1] in keys(m.pseudo_observables)
            if plotvars[i-1] == :Forward5YearRealNaturalRate || plotvars[i-1] == :RealNaturalRate
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", ylims=(-0.1, 0.1),xticks=my_xticks)
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="",xticks=my_xticks)
            end
        end
    end

    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "%")
    xlabel!(p[i], "Quarter")
end

plot!(p)

pdf_path = joinpath(saveroot,"Final Paper","figures","interest_rate_peg_with_contemporaneous_innovation_0_and_rest_FG.pdf")
savefig(pdf_path)

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    collect(1:PlotT),
    Dict(
        "Combined monetary policy shocks" => vec(shk_weights_store[:, peg_horizon-1]),
        "Policy rate"                     => vec(FGplotmat[:, 1, peg_horizon-1]),
        "Inflation"                       => vec(FGplotmat[:, 2, peg_horizon-1]),
        "Output"                          => vec(FGplotmat[:, 3, peg_horizon-1]),
        "r* (Forward 5-year real natural rate)" => vec(FGplotmat[:, 4, peg_horizon-1]),
        "Ex-ante real rate"               => vec(FGplotmat[:, 5, peg_horizon-1]),
        "Zero line"                       => zeros(PlotT),
    )
)
# ================================================





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
peg_horizon = 6;
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
idx = [1,2,3,4,5,6,7,8]   # example horizons
#idx = abs.(y) .> 10^-6

p1 = plot(
    x[idx],
    y[idx]*4,
    seriestype = :scatter,
    marker = :star5,
    markersize = 6,
    markercolor = :blue,
    title = "Permanent Liquidity Innovations (APR)",
    label = "", xticks = my_xticks
)
    ylabel!(p1, "%")
    xlabel!(p1, "Quarter")
plot!(p1, x, zeros(horizon), lc=:black, lw=2, label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Policy rate (APR)", xticks = my_xticks)#
    ylabel!(p2, "%")
    xlabel!(p2, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_corepce],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Inflation (%, qoq annualized)", xticks = my_xticks)#
    ylabel!(p3, "%")
    xlabel!(p3, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,states[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Output (% dev from SS)", xticks = my_xticks)#
    ylabel!(p4, "%")
    xlabel!(p4, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="r* (Forward 5-year real natural rate, APR)", ylims = (-1.1, 1.1), xticks = my_xticks)#
    ylabel!(p5, "%")
    xlabel!(p5, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:b_liqp_sh]]*4,title="Ex-ante real rate (APR)", xticks = my_xticks)#
    ylabel!(p6, "%")
    xlabel!(p6, "Quarter")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,p5,p6,layout=(3,2), legend=false)
plot!(size=(960,540))
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