using DSGE;
path = dirname(@__FILE__)


m = Model1010("ss20");

system = compute_system(m)

shock_names = [:b_safetil_sh, :b_safep_sh, :b_liqtil_sh, :b_liqp_sh];
shock_values = ones(4,1);
horizon = 80;
variable_names = [:y_t,:y_f_t,:Rstarn,:rstar,:i_t];
titlelist = keys(m.observable_mappings);

var_name = :y_t;
shock_name = :b_safetil_sh;
var_value = 10.
shock_values = 1.

# states_irf, obs_irf, pseudo_irf = impulse_responses(m, system, horizon, shock_name , var_name, var_value)


shock_names = [:b_safetil_sh]
shock_values = [1.0]
states_irf, obs_irf, pseudo_irf  = impulse_responses(m, system, horizon, shock_names, shock_values)




using Plots # no need for `using Plots` as that is reexported here
plot(obs_irf[1,:,4],label="IRF of output gap to permanent safety shock",lw=4) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("Impulse Response Function of Output \n in Del Negro et al. (2018)")
xlabel!("% deviation from SS")
ylabel!("Horizon")

plot(obs_irf[6,:,4],label="IRF of policy rate to permanent safety shock",lw=4) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("Impulse Response Function of policy rate \n in Del Negro et al. (2018)")
xlabel!("% deviation from SS")
ylabel!("Horizon")



var_name1 = :y_t;
shock_name1 = :b_liqtil_sh;
var_value = 10.
states_irf, obs_irf1, pseudo_irf = impulse_responses(m, system, horizon, shock_name1 , var_name1, var_value)

plot(obs_irf1[1,:,2],label="IRF of output gap to permanent liquidity shock",lw=4,lc=:red) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("Impulse Response Function of Output \n in Del Negro et al. (2018)")
xlabel!("% deviation from SS")
ylabel!("Horizon")


plot(obs_irf1[6,:,2],label="IRF of policy rate  to permanent liquidity shock",lw=4,lc=:red) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("Impulse Response Function of policy rate \n in Del Negro et al. (2018)")
xlabel!("% deviation from SS")
ylabel!("Horizon")


