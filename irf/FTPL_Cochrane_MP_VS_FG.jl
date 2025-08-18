using DSGE;
using Plots # no need for `using Plots` as that is reexported here

path = dirname(@__FILE__)
horizon  = 80

"""
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
p5 = plot(1:horizon,obs[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,obs[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))
savefig( "irf/FTPL_Equilibrium_IRF_Policy_rate_with_MP_shock.pdf")   # saves the plot from p as a .pdf vector graphic


########################################################################
# An ever changing shock for the end of the peg that implements the FG #
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
p5 = plot(1:horizon,obs[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,obs[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

savefig( "irf/FTPL_Equilibrium_IRF_Policy_rate_with_promise_for_period.pdf")   # saves the plot from p as a .pdf vector graphic

#  OLD CODE

#####################
# Standard MP shock #
#####################

# m = Model1010("ss20");
# system = compute_system(m)
# shock_name = :rm_sh # Select MP to implement the specific path in state variable 
# var_name = :obs_nominalrate # Select the targeted state variable
# var_value = -1.0  # Select the depth of the path
# peg_horizon =6;

# # Setup - copied from impulse_responses.jl
#     var_names, var_class =
#     if var_name in keys(m.endogenous_states)
#         m.endogenous_states, :states
#     else
#         m.observables, :obs
#     end
#     exo          = m.exogenous_shocks
#     nshocks      = size(system[:RRR], 2)
#     nstates      = size(system[:TTT], 1)
#     nobs         = size(system[:ZZ], 1)
#     npseudo      = size(system[:ZZ_pseudo], 1)

#     states = zeros(nstates, horizon, nshocks)
#     obs    = zeros(nobs,    horizon, nshocks)
#     pseudo = zeros(npseudo, horizon, nshocks)

#     # Set constant system matrices to 0
#     system = DSGE.zero_system_constants(system)

#     s_0 = zeros(nstates)

#     # Isolate single shock
#     shocks = zeros(nshocks, horizon)
#     for t = 1:peg_horizon
#         if var_class == :states
#             var_value_att = var_value - obs[m.endogenous_states[var_name],t, m.exogenous_shocks[shock_name]]
#             shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_state_value(var_value_att,
#                                                                         var_names[var_name],
#                                                                         exo[shock_name],
#                                                                         system[:RRR])
#         else # == :obs
#             var_value_att = var_value - obs[m.observables[var_name],t, m.exogenous_shocks[shock_name]]
#             shocks[exo[shock_name], t] = DSGE.obtain_shock_from_desired_obs_value(var_value_att,
#                                                                         var_names[var_name],
#                                                                         exo[shock_name],
#                                                                         system[:ZZ],
#                                                                         system[:RRR])
#         end
    
#     # Iterate state space forward
#     states[:, :, exo[shock_name]], obs[:, :, exo[shock_name]], pseudo[:, :, exo[shock_name]], _ = forecast(system, s_0, shocks)
#     end


# using Plots
# p1 = plot(1:horizon,states[m.endogenous_states[:rm_t],:, m.exogenous_shocks[:rm_sh]],title="Monetary policy shock")
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
# p2 = plot(1:horizon,obs[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh]],title="Policy rate")
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
# p3 = plot(1:horizon,obs[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]],title="Inflation")
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
# p4 = plot(1:horizon,obs[m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_sh]],title="Output")#
# plot!(zeros(horizon,1),lc=:black,lw=2,label="")
# plot(p1, p2, p3, p4, layout=(2,2), legend=false)
# plot!(size=(960,540))

# savefig( "irf/FTPL_Equilibrium_IRF_Policy_rate_with_MP_shock.pdf")   # saves the plot from p as a .pdf vector graphic
