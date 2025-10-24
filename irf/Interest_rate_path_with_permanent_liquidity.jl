using DSGE;
using Plots # no need for `using Plots` as that is reexported here

mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "output_data")
saveroot = joinpath(basepath, "dsge")

horizon  = 40
peg_horizon = 6 # Length of the peg in periods

var_name = :obs_nominalrate # Select the targeted state variable
m = Model1010("ss20");
mode_file = joinpath(dataroot, "m1010","ss20","estimate","raw", "paramsmode_vint=161223.h5")
specify_mode!(m, mode_file)
system = compute_system(m)

states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)
p1 = plot(1:horizon,[states_irf[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]],title="Output", label=["US model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,[obs_irf[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:b_liqp_sh]]] ,title="Inflation", label=["US model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ obs_irf[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]] ] ,title="Policy rate", label=["US model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh ]] ] ,title="r*", label=["US model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)


#     nshocks = size(system[:RRR], 2)
#     nstates = size(system[:TTT], 1)
#     s_0 = zeros(nstates)
#     shocks = zeros(nshocks, horizon)


# shock_name = :b_safep_sh # Select MP to implement the specific path in state variable 
# var_name = :obs_nominalrate # Select the targeted state variable
# var_value = -1.0  # Select the depth of the path
# peg_horizon = 6;

# for t = 1:peg_horizon
#             var_value_att = var_value - obs[m.observables[var_name],t, m.exogenous_shocks[shock_name]]
#             shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(var_value_att,
#                                                                         var_names[var_name],
#                                                                         exo[shock_name],
#                                                                         system[:ZZ],
#                                                                         system[:RRR])

#         states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)
# end


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
                    irf = states[peg_ind, 1] # Impact of a unit shock at t on state at t
                else 
                    irf = obs[peg_ind,1] # Impact of a unit shock at t on obs at t
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
# Permanent liquidity shock #
#####################

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
p1 = plot(1:horizon,states[m.endogenous_states[:b_liqp_t],:, m.exogenous_shocks[:b_liqp_sh]],title="Permanent liquidity shock")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]],title="Policy rate")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p3 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:b_liqp_sh]],title="Inflation")
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p4 = plot(1:horizon,states[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]],title="Output")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p5 = plot(1:horizon,pseudo[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]],title="r* (Forward 5-year real natural rate)")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,pseudo[m.pseudo_observables[:ExAnteRealRate],:, m.exogenous_shocks[:b_liqp_sh]],title="Ex-ante real rate")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")

plot(p1, p2, p3, p4,p5,p6,layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "irf/IRF_rate_peg_with_permanent_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic
