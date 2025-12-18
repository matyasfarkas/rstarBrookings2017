using DSGE;
using Plots # no need for `using Plots` as that is reexported here

# =========================
# Excel export (ADDED ONLY)
# =========================
using XLSX
using Dates

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
# Your original script
# =========================

path = dirname(@__FILE__)
horizon  = 40
peg_horizon = 6 # Length of the peg in periods
m = Model1010("ss20");

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


    for t in 1:horizon
        
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
# Standard MP shock #
#####################

m = Model1010("ss20");
system = compute_system(m)
shock_name = :rm_sh # Select MP to implement the specific path in state variable 
var_name = :obs_nominalrate # Select the targeted state variable
var_value = -1.0  # Select the depth of the path
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
p1 = plot(1:horizon,states[m.endogenous_states[:rm_t],:, m.exogenous_shocks[:rm_sh]],title="Monetary policy shock")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]],title="Policy rate")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]],title="Inflation")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]],title="Output")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh]],title="r* (Forward 5-year real natural rate)", ylims = (-0.1, 0.1))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:rm_sh]],title="Ex-ante real rate")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,p5,p6,layout=(3,2), legend=false)
plot!(size=(960,540))
pdf_path = "irf/FTPL_Equilibrium_IRF_Policy_rate_with_MP_shock.pdf"
savefig(pdf_path)   # saves the plot from p as a .pdf vector graphic

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    tgrid,
    Dict(
        "Monetary policy shock (rm_t)"            => vec(states[m.endogenous_states[:rm_t],:, m.exogenous_shocks[:rm_sh]]),
        "Policy rate (obs_nominalrate)"          => vec(obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]]),
        "Inflation (obs_gdpdeflator)"            => vec(obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]]),
        "Output (obs_gdp)"                       => vec(obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]]),
        "r* (Forward5YearRealNaturalRate)"       => vec(pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh]]),
        "Ex-ante real rate (ExAnteRealRate)"     => vec(pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:rm_sh]]),
        "Zero line"                              => zeros(horizon),
    )
)
# ================================================

desired_path = vec(var_value *ones(peg_horizon)) # Desired path for the state variable
shock_inds = [m.exogenous_shocks[:rm_sh],m.exogenous_shocks[:rm_sh],m.exogenous_shocks[:rm_sh],m.exogenous_shocks[:rm_sh],m.exogenous_shocks[:rm_sh],m.exogenous_shocks[:rm_sh]] # This replicates the only MP path case
shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
states, obs, pseudo, _ = forecast(system, s_0, hcat(shocks_path, zeros(size(collect(m.exogenous_shocks),1), horizon-peg_horizon)))
using Plots
p1 = plot(1:horizon,states[m.endogenous_states[:rm_t],:],title="Combined monetary policy shocks")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,obs[m.observables[:obs_gdp],:],title="Output")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)", ylims = (-0.1, 0.1))#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate", ylims = (-0.1, 0.1))#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

########################################################################
# An ever changing shock for the end of the peg that implements the FG path #
########################################################################
desired_path = vec(var_value *ones(peg_horizon)) # Desired path for the state variable
var_name =:obs_nominalrate
mp_ind = m.exogenous_shocks[:rm_sh]
fg_inds = [m.exogenous_shocks[:rm_shl6],
           m.exogenous_shocks[:rm_shl5],
           m.exogenous_shocks[:rm_shl4],
           m.exogenous_shocks[:rm_shl3],
           m.exogenous_shocks[:rm_shl2],
           m.exogenous_shocks[:rm_shl1]] # This uses FG1-FG6 shocks
shock_inds = vcat(fg_inds,mp_ind) # Using mp and 1-6 horizon FG shocks
shocks_path = obtain_shocks_from_desired_state_path_iterative(desired_path,m, var_name, shock_inds, system)
print("Shocks path: ",shocks_path)
states, obs, pseudo, _ = forecast(system, s_0, hcat(shocks_path, zeros(size(collect(m.exogenous_shocks),1), horizon-peg_horizon)))

using Plots
p1 = plot(1:horizon,states[m.endogenous_states[:rm_t],:],title="Combined monetary policy shocks")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:],title="Policy rate")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:],title="Inflation")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,obs[m.observables[:obs_gdp],:],title="Output")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)", ylims = (-0.1, 0.1))#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate")#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

pdf_path = "irf/FTPL_Equilibrium_IRF_Policy_rate_with_promise_for_period.pdf"
savefig(pdf_path)   # saves the plot from p as a .pdf vector graphic

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    tgrid,
    Dict(
        "Combined monetary policy shocks (rm_t)" => vec(states[m.endogenous_states[:rm_t],:]),
        "Policy rate (obs_nominalrate)"         => vec(obs[m.observables[:obs_nominalrate],:]),
        "Inflation (obs_gdpdeflator)"           => vec(obs[m.observables[:obs_gdpdeflator],:]),
        "Output (obs_gdp)"                      => vec(obs[m.observables[:obs_gdp],:]),
        "r* (Forward5YearRealNaturalRate)"      => vec(pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:]),
        "Ex-ante real rate (ExAnteRealRate)"    => vec(pseudo[m.pseudo_observables[:ExAnteRealRate],:]),
        "Zero line"                             => zeros(horizon),
    )
)
# ================================================


# ########################################################################
# Proper forward guidance implementation by Matyas Farkas, IMF 19/08/2025 based on the Matlab code from Jesper Linde and Zoltan Jakab
# ########################################################################

using LinearAlgebra
using Plots

# --- Step 1: Compute IRFs for each shock ---
PlotT = horizon # Total IRF horizon to plot
plotvars = [:obs_gdp, :obs_gdpdeflator, :obs_nominalrate] # Output, Inflation, Policy Rate
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
    FG_vec = fill(-1.0, FGdur) # Desired policy rate path

    # Build IRF matrix for policy rate
    R_mp_mat = zeros(FGdur, FGdur)
    for s = 1:FGdur
        R_mp_mat[:, s] .= irfmat[1:FGdur, 3, s] # 3rd var is policy rate
    end

    # Solve for shock weights
    shk_weights = R_mp_mat \ FG_vec
    shk_weights_store[1:FGdur, hidx] .= shk_weights

    # Construct total IRFs for each variable
    for s = 1:FGdur
        FGplotmat[:, :, hidx] .+= shk_weights[s] .* irfmat[:, :, s]
    end
end


# Propagate for full PlotT
# --- Step 3: Plot results ---
titles = ["Output", "Inflation", "Nominal Policy Rate"]
TT = 1:PlotT

p = plot(layout=(1,3), size=(1200,400))
for i = 1:nvars
    plot!(p[i], TT, FGplotmat[:, i, peg_horizon-1], lw=2, label="FG IRF (horizon=$peg_horizon)")
    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "Percent")
    xlabel!(p[i], "Quarter")
end
plot!(p)
pdf_path = "irf/FG_6horizon_policy_rate_output_inflation.pdf"
savefig(pdf_path)

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    collect(1:PlotT),
    Dict(
        "Output"       => vec(FGplotmat[:, 1, peg_horizon-1]),
        "Inflation"    => vec(FGplotmat[:, 2, peg_horizon-1]),
        "Policy Rate"  => vec(FGplotmat[:, 3, peg_horizon-1]),
        "Zero line"    => zeros(PlotT),
    )
)
# ================================================

# Optional: plot weights for each horizon
pw = plot(1:PlotT, shk_weights_store[:, peg_horizon-1], lw=2, label="Shock weights")
xlabel!("Quarter ahead shocks")
ylabel!("Weight")
title!("Shock Weights for FG horizon $peg_horizon")
pdf_path = "irf/FG_6horizon_shock_weights.pdf"
savefig(pdf_path)

# ===== ADDED ONLY: export underlying plotted data =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    collect(1:PlotT),
    Dict(
        "Shock weights" => vec(shk_weights_store[:, peg_horizon-1]),
        "Zero line"     => zeros(PlotT),
    )
)
# ================================================

#########################################################################
# Adding the other variables to plot
#########################################################################

plotvars = [ :obs_nominalrate,:obs_gdpdeflator,  :obs_gdp, :Forward5YearRealNaturalRate, :ExAnteRealRate] 

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
    FG_vec = fill(-1.0, FGdur) # Desired policy rate path

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
        plot!(p[i], TT, shk_weights_store[:, peg_horizon-1], lw=2, label="")
    else
        if plotvars[i-1] in keys(m.observables)
            plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="")
        elseif plotvars[i-1] in keys(m.pseudo_observables)
            if plotvars[i-1] == :Forward5YearRealNaturalRate || plotvars[i-1] == :RealNaturalRate
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", ylims=(-0.1, 0.1))
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="")
            end
        end
    end
    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "Percent")
    xlabel!(p[i], "Quarter")
end
plot!(p)
pdf_path = "irf/FG_6horizon_policy_rate_output_inflation_and_rstar.pdf"
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
# Adding the other variables to plot
#########################################################################

plotvars = [ :obs_nominalrate,:obs_gdpdeflator,  :obs_gdp, :Forward5YearRealNaturalRate, :ExAnteRealRate] 

titles = ["Combined monetary policy shocks","Policy rate", "Inflation", "Output", "r* (Forward 5-year real natural rate)", "Ex-ante real rate"]
nvars = length(plotvars)


shock_syms = [:rm_shl1, :rm_shl2, :rm_shl3, :rm_shl4, :rm_shl5, :rm_shl6] # MP + 1-6 FG shocks

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
    FG_vec = fill(-1.0, FGdur) # Desired policy rate path

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
        plot!(p[i], TT, shk_weights_store[:, peg_horizon-1], lw=2, label="")
    else
        if plotvars[i-1] in keys(m.observables)
            plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="")
        elseif plotvars[i-1] in keys(m.pseudo_observables)
            if plotvars[i-1] == :Forward5YearRealNaturalRate || plotvars[i-1] == :RealNaturalRate
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="", ylims=(-0.1, 0.1))
            else
                plot!(p[i], TT, FGplotmat[:, i-1, peg_horizon-1], lw=2, label="")
            end
        end
    end
    plot!(p[i], TT, zeros(PlotT), lc=:black, lw=1, label="")
    title!(p[i], titles[i])
    ylabel!(p[i], "Percent")
    xlabel!(p[i], "Quarter")
end
plot!(p)
pdf_path = "irf/FG_6horizon_policy_rate_output_inflation_and_rstar_FISHERIAN.pdf"
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


# (Your remaining commented-out blocks and OLD CODE remain unchanged below.)
#  OLD CODE

#####################
# Standard MP shock #
#####################

m = Model1010("ss20");
system = compute_system(m)
shock_name = :rm_sh # Select MP to implement the specific path in state variable 
var_name = :obs_nominalrate # Select the targeted state variable
var_value = -1.0  # Select the depth of the path
peg_horizon =6;

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
p1 = plot(1:horizon,states[m.endogenous_states[:rm_t],:, m.exogenous_shocks[:rm_sh]],title="Monetary policy shock")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]],title="Policy rate")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]],title="Inflation")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]],title="Output")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh]],title="r* (Forward 5-year real natural rate)", ylims = (-0.1, 0.1))#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:RealNaturalRate],:, m.exogenous_shocks[:rm_sh]],title="Real natural rate", ylims = (-0.1, 0.1))#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,p5,p6,layout=(3,2), legend=false)
plot!(size=(960,540))

pdf_path = "irf/FTPL_Equilibrium_IRF_Policy_rate_with_MP_shock.pdf"
savefig(pdf_path)   # saves the plot from p as a .pdf vector graphic

# ===== ADDED ONLY: export underlying plotted data (same PDF name sheet; overwrites sheet) =====
export_series_xlsx(
    xlsx_out,
    sheet_from_pdf(pdf_path),
    tgrid,
    Dict(
        "Monetary policy shock (rm_t)"            => vec(states[m.endogenous_states[:rm_t],:, m.exogenous_shocks[:rm_sh]]),
        "Policy rate (obs_nominalrate)"          => vec(obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]]),
        "Inflation (obs_gdpdeflator)"            => vec(obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]]),
        "Output (obs_gdp)"                       => vec(obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]]),
        "r* (Forward5YearRealNaturalRate)"       => vec(pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh]]),
        "Real natural rate (RealNaturalRate)"    => vec(pseudo[m.pseudo_observables[:RealNaturalRate],:, m.exogenous_shocks[:rm_sh]]),
        "Zero line"                              => zeros(horizon),
    )
)
# ================================================
