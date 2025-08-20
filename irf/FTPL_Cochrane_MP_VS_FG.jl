using DSGE;
using Plots # no need for `using Plots` as that is reexported here

path = dirname(@__FILE__)
horizon  = 20
peg_horizon = 6 # Length of the peg in periods

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
desired_path = vec(-1.0 * ones(peg_horizon)) # Desired path for the state variable

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
p5 = plot(1:horizon,obs[m.pseudo_observables[:Forward5YearRealNaturalRate],:],title="r* (Forward 5-year real natural rate)")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
p6 = plot(1:horizon,obs[m.pseudo_observables[:RealNaturalRate],:],title="Real natural rate")#
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot(p1, p2, p3, p4,p5,p6, layout=(3,2), legend=false)
plot!(size=(960,540))

savefig( "irf/FTPL_Equilibrium_IRF_Policy_rate_with_promise_for_period.pdf")   # saves the plot from p as a .pdf vector graphic


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
savefig("irf/FG_5horizon_policy_rate_output_inflation.pdf")

# Optional: plot weights for each horizon
pw = plot(1:PlotT, shk_weights_store[:, peg_horizon-1], lw=2, label="Shock weights")
xlabel!("Quarter ahead shocks")
ylabel!("Weight")
title!("Shock Weights for FG horizon $peg_horizon")
savefig("irf/FG_5horizon_shock_weights.pdf")

#########################################################################




# ########################################################################
# # Expectation errors on a length of a MP peg - Agents anticipate peg to be implemented  and are suprized consequtively that it gets extended                          #
# ########################################################################



# """
#     obtain_shocks_from_desired_multi_state_path_iterative(x::Matrix{Float64}, m::AbstractDSGEModel, var_names::Vector{Symbol},
#                                                           shock_inds::Vector{Vector{Int}}, system::System{Float64})

# Given a desired path `x` for multiple variables (size k × horizon), a vector of variable names `var_names` (length k),
# and a vector of shock index vectors `shock_inds` (length horizon, each element is a vector of k shock indices),
# back out the necessary values of shocks for each period such that each target is met, using the model's forecast function.

# Returns a matrix of required shocks of size (nshocks, horizon).
# """
# function obtain_shocks_from_desired_multi_state_path_iterative(x::Matrix{Union{Float64,Missing}}, m::AbstractDSGEModel, var_names::Vector{Symbol},
#                                                                shock_inds::Vector{Vector{Int}}, system::System{Float64})
#     horizon = size(x, 2)
#     k = size(x, 1)
#     nshocks = size(system[:RRR], 2)
#     nstates = size(system[:TTT], 1)
#     s_0 = zeros(nstates)
#     shocks = zeros(nshocks, horizon)

#     # Get indices and classes for each target variable
#     var_info = Vector{Tuple{Int, Symbol}}(undef, k)
#     for i in 1:k
#         if var_names[i] in keys(m.endogenous_states)
#             var_info[i] = (m.endogenous_states[var_names[i]], :states)
#         else
#             var_info[i] = (m.observables[var_names[i]], :obs)
#         end
#     end

#     for t in 1:horizon
#         # Find non-missing targets for this period
#         target_inds = findall(!ismissing, x[:, t])
#         n_targets = length(target_inds)
#         if n_targets == 0
#             continue
#         end

#         # Build IRF matrix for this period (n_targets × n_targets shocks)
#         IRFmat = zeros(n_targets, n_targets)
#         for i in 1:n_targets
#             for j in 1:n_targets
#                 test_shocks = zeros(nshocks, horizon)
#                 test_shocks[shock_inds[t][target_inds[j]], t] = 1.0
#                 states, obs, _ = forecast(system, s_0, test_shocks)
#                 idx, cls = var_info[target_inds[i]]
#                 IRFmat[i, j] = cls == :states ? states[idx, t] : obs[idx, t]
#             end
#         end

#         # Compute effect of previous shocks
#         prev_effect = zeros(n_targets)
#         if t > 1
#             prev_shocks = shocks[:, 1:t-1]
#             prev_states, prev_obs, _ = forecast(system, s_0, hcat(prev_shocks, zeros(nshocks, horizon-t+1)))
#             for i in 1:n_targets
#                 idx, cls = var_info[target_inds[i]]
#                 prev_effect[i] = cls == :states ? prev_states[idx, t] : prev_obs[idx, t]
#             end
#         end

#         # Solve for required shocks at time t
#         shocks_t = IRFmat \ (collect(skipmissing(x[:, t])) - prev_effect)
#         for j in 1:n_targets
#             shocks[shock_inds[t][target_inds[j]], t] = shocks_t[j]
#         end
#     end

#     return shocks
# end


# # Setup
# horizon_topeg = 7 # Enough to show extension
# nstates = size(system[:TTT], 1)
# nshocks = size(system[:RRR], 2)
# s_0 = zeros(nstates)

# # OrderedDict indices for expected short-term rates
# rate_syms = [:rm_t, :rm_tl1, :rm_tl2, :rm_tl3, :rm_tl4, :rm_tl5]
# rate_inds = [m.endogenous_states[sym] for sym in rate_syms]

# # For period 0: target -1 for current and next 6 expected rates
# targets_0 = fill(-1.0, horizon_topeg-1)
# var_names_0 = rate_syms
# shock_inds_0 = [m.exogenous_shocks[:rm_sh], m.exogenous_shocks[:rm_shl1], m.exogenous_shocks[:rm_shl2],
#                 m.exogenous_shocks[:rm_shl3], m.exogenous_shocks[:rm_shl4], m.exogenous_shocks[:rm_shl5]]

# # For period 1: extend horizon to 7 ahead (add :rm_tl6)
# rate_syms_ext = [:rm_t, :rm_tl1, :rm_tl2, :rm_tl3, :rm_tl4, :rm_tl5, :rm_tl6]
# targets_1 = fill(-1.0, horizon_topeg)
# var_names_1 = rate_syms_ext
# shock_inds_1 = [m.exogenous_shocks[:rm_sh], m.exogenous_shocks[:rm_shl1], m.exogenous_shocks[:rm_shl2],
#                 m.exogenous_shocks[:rm_shl3], m.exogenous_shocks[:rm_shl4], m.exogenous_shocks[:rm_shl5],
#                 m.exogenous_shocks[:rm_shl6]]

# # Build target matrix and shock index vector for each period
# x_targets = Matrix{Union{Float64,Missing}}(missing, 7, horizon)
# shock_inds = Vector{Vector{Int}}(undef, horizon)
# var_names = rate_syms_ext

# # Period 0: target first 6 rates, last entry is missing
# x_targets[1:6, 1] .= -1.0
# x_targets[7, 1] = missing
# shock_inds[1] = shock_inds_0

# # Period 1: target all 7 rates
# x_targets[:, 2] .= -1.0
# shock_inds[2] = shock_inds_1

# # Periods 3+: set all entries to missing (or set some to -1.0 if you want to target)
# for t in 3:horizon
#     x_targets[:, t] .= missing
#     shock_inds[t] = shock_inds_1
# end
# # Compute shocks
# shocks = obtain_shocks_from_desired_multi_state_path_iterative(x_targets, m, var_names, shock_inds, system)

# # Simulate model
# states, obs, pseudo, _ = forecast(system, s_0, shocks)

# # Plot results
# using Plots

# p1 = plot(1:horizon, states[m.endogenous_states[:rm_t], :], title="Combined monetary policy shocks")
# plot!(zeros(horizon), lc=:black, lw=2, label="")

# p2 = plot(1:horizon, obs[m.observables[:obs_nominalrate], :], title="Policy rate")
# plot!(zeros(horizon), lc=:black, lw=2, label="")

# p3 = plot(1:horizon, obs[m.observables[:obs_gdpdeflator], :], title="Inflation")
# plot!(zeros(horizon), lc=:black, lw=2, label="")

# p4 = plot(1:horizon, obs[m.observables[:obs_gdp], :], title="Output")
# plot!(zeros(horizon), lc=:black, lw=2, label="")

# p5 = plot(1:horizon, obs[m.pseudo_observables[:Forward5YearRealNaturalRate], :], title="r* (Forward 5-year real natural rate)")
# plot!(zeros(horizon), lc=:black, lw=2, label="")

# p6 = plot(1:horizon, obs[m.pseudo_observables[:RealNaturalRate], :], title="Real natural rate")
# plot!(zeros(horizon), lc=:black, lw=2, label="")

# plot(p1, p2, p3, p4, p5, p6, layout=(3,2), legend=false)
# plot!(size=(960,540))

# savefig("irf/FTPL_Equilibrium_IRF_Extended_Expectations.pdf")




# rate_syms_ext = [:rm_t, :rm_tl1, :rm_tl2, :rm_tl3, :rm_tl4, :rm_tl5, :rm_tl6,:R_t1 ]

# plots = []
# for sym in rate_syms_ext
#     idx = m.endogenous_states[sym]
#     p = plot(1:horizon, states[idx, :], title="Expectation: $(string(sym))", lw=2)
#     plot!(p, zeros(horizon), lc=:black, lw=1, label="")
#     push!(plots, p)
# end

# plot(plots..., layout=(4,2), legend=false, size=(960,540))
# savefig("irf/FTPL_Equilibrium_IRF_Policy_Expectations.pdf")


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
