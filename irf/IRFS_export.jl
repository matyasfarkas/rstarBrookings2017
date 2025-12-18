using DSGE;
using Plots # no need for `using Plots` as that is reexported here

# --- ADDED (export only; does not affect computations) ---
using XLSX
using Dates

"""
Append/replace a sheet in an Excel workbook, writing time + named series.
This does NOT compute anything; it only writes the vectors you pass in.
"""
function export_series_xlsx(filepath::AbstractString,
                            sheetname::AbstractString,
                            t,
                            series::Dict{String,<:AbstractVector})

    mkpath(dirname(filepath))
    mode = isfile(filepath) ? "rw" : "w"

    XLSX.openxlsx(filepath, mode=mode) do xf
        # Get or create sheet (no deletes required)
        sh =
            if sheetname in XLSX.sheetnames(xf)
                xf[sheetname]
            else
                XLSX.addsheet!(xf, sheetname)
            end

        # Write header
        sh[1,1] = "t"
        colnames = collect(keys(series))
        for (j, name) in enumerate(colnames)
            sh[1, j+1] = name
        end

        # Write data
        T = length(t)
        for i in 1:T
            sh[i+1, 1] = t[i]
            for (j, name) in enumerate(colnames)
                sh[i+1, j+1] = series[name][i]
            end
        end

        # Optional: stamp
        sh[1, length(colnames)+3] = "exported_at"
        sh[1, length(colnames)+4] = string(now())
    end

    return filepath
end
# --- END ADDED ---

path = dirname(@__FILE__)
horizon  = 20
peg_horizon = 6 # Length of the peg in periods
m = Model1010("ss20");

var_name = :obs_nominalrate # Select the targeted state variable
system = compute_system(m)

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

p1 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]],title="Policy rate",ylims = (-1.25, 0.75))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]],title="Inflation",ylims = (-0.4, 0.2))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]],title="Output",ylims = (-3,3))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:rm_sh]],title="Ex-ante real rate",ylims = (-1.2,0.8))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,layout=(2,2), legend=false)
plot!(size=(960,540))

# --- ADDED: export underlying series for this figure ---
xlsx_out = joinpath(path, "irf", "IRF_underlying_data.xlsx")
tgrid = collect(1:horizon)
export_series_xlsx(
    xlsx_out,
    "MPshock_rm_sh",
    tgrid,
    Dict(
        "Policy rate"       => vec(obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]]),
        "Inflation"         => vec(obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]]),
        "Output"            => vec(obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]]),
        "Ex-ante real rate" => vec(pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:rm_sh]]),
        "Zero line"         => zeros(horizon),
    )
)
# --- END ADDED ---

savefig( "irf/Presentation_FTPL_Equilibrium_IRF_Policy_rate_with_MP_shock.pdf")   # saves the plot from p as a .pdf vector graphic


########################################################################
# A constant shock for the end of the peg that implements the desired path #
########################################################################
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

# --- ADDED: export underlying series for this figure ---
export_series_xlsx(
    xlsx_out,
    "ConstShockPath_3x2",
    tgrid,
    Dict(
        "Combined monetary policy shocks (rm_t)" => vec(states[m.endogenous_states[:rm_t],:]),
        "Policy rate"                            => vec(obs[m.observables[:obs_nominalrate],:]),
        "Inflation"                              => vec(obs[m.observables[:obs_gdpdeflator],:]),
        "Output"                                 => vec(obs[m.observables[:obs_gdp],:]),
        "r* (Forward 5-year real natural rate)"   => vec(pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:]),
        "Ex-ante real rate"                       => vec(pseudo[m.pseudo_observables[:ExAnteRealRate],:]),
        "Zero line"                               => zeros(horizon),
    )
)
# --- END ADDED ---


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
p1 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:],title="Policy rate",ylims = (-1.25, 0.75))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:],title="Inflation",ylims = (-0.4, 0.2))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdp],:],title="Output",ylims = (-3,3))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:],title="Ex-ante real rate",ylims = (-1.2,0.8))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
plot!(size=(960,540))

# --- ADDED: export underlying series for this figure ---
export_series_xlsx(
    xlsx_out,
    "FGpath_2x2",
    tgrid,
    Dict(
        "Policy rate"       => vec(obs[m.observables[:obs_nominalrate],:]),
        "Inflation"         => vec(obs[m.observables[:obs_gdpdeflator],:]),
        "Output"            => vec(obs[m.observables[:obs_gdp],:]),
        "Ex-ante real rate" => vec(pseudo[m.pseudo_observables[:ExAnteRealRate],:]),
        "Zero line"         => zeros(horizon),
    )
)
# --- END ADDED ---

savefig( "irf/Presentation_FTPL_Equilibrium_IRF_Policy_rate_with_promise_for_period.pdf")   # saves the plot from p as a .pdf vector graphic


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
p1 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:],title="Policy rate",ylims = (-1.25, 0.75))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:],title="Inflation", ylims = (-0.5, 0.25))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdp],:],title="Output", ylims = (-3,1))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)", ylims = (-0.1, 0.1))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
plot!(size=(960,540))

# --- ADDED: export underlying series for this figure ---
export_series_xlsx(
    xlsx_out,
    "FGpath_with_rstar_2x2",
    tgrid,
    Dict(
        "Policy rate"                          => vec(obs[m.observables[:obs_nominalrate],:]),
        "Inflation"                            => vec(obs[m.observables[:obs_gdpdeflator],:]),
        "Output"                               => vec(obs[m.observables[:obs_gdp],:]),
        "r* (Forward 5-year real natural rate)"=> vec(pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:]),
        "Zero line"                            => zeros(horizon),
    )
)
# --- END ADDED ---

savefig( "irf/Presentation_vs_pliq_FTPL_Equilibrium_IRF_Policy_rate_with_promise_for_period.pdf")   # saves the plot from p as a .pdf vector graphic
##############################


#####################
# Permanent liquidity shock #
#####################


mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "output_data")
saveroot = joinpath(basepath, "dsge")

var_name = :obs_nominalrate # Select the targeted state variable
m = Model1010("ss20");
mode_file = joinpath(dataroot, "m1010","ss20","estimate","raw", "paramsmode_vint=161223.h5")
specify_mode!(m, mode_file)
system = compute_system(m)

states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)
shock_name = :b_liqp_sh # Select MP to implement the specific path in state variable 
var_name = :obs_nominalrate # Select the targeted state variable
var_value = -1.0  # Select the depth of the path
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
p1 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]],title="Policy rate",ylims = (-1.25, 0.75))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:b_liqp_sh]],title="Inflation",ylims = (-0.5, 0.25))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,states[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]],title="Output",ylims = (-3,1))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]],title="r* (Forward 5-year real natural rate)")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,layout=(2,2), legend=false)
plot!(size=(960,540))

# --- ADDED: export underlying series for this figure ---
export_series_xlsx(
    xlsx_out,
    "PermLiquidity_b_liqp_sh_2x2",
    tgrid,
    Dict(
        "Policy rate"                          => vec(obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]]),
        "Inflation"                            => vec(obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:b_liqp_sh]]),
        "Output (y_t state)"                    => vec(states[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]),
        "r* (Forward 5-year real natural rate)" => vec(pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]]),
        "Zero line"                             => zeros(horizon),
    )
)
# --- END ADDED ---

savefig( "irf/Presentation_IRF_rate_peg_with_permanent_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic







using DSGE;
using Plots # no need for `using Plots` as that is reexported here

path = dirname(@__FILE__)
horizon  = 20
peg_horizon = 6 # Length of the peg in periods
m = Model1010("ss20");

var_name = :obs_nominalrate # Select the targeted state variable
system = compute_system(m)

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
var_value = 1.0  # Select the depth of the path
peg_horizon = 1;

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

p1 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]],title="Policy rate",ylims = (-1.25, 0.75))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]],title="Inflation",ylims = (-0.4, 0.2))
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]],title="Output",ylims = (-3,3))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:rm_sh]],title="Ex-ante real rate",ylims = (-1.2,0.8))#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,layout=(2,2), legend=false)
plot!(size=(960,540))
