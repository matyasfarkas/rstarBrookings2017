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
plot(obs_irf[1,:,4],label="IRF of GDP per capita to permanent safety shock",lw=4) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("IRF of gdp per capita to permanent safety shock \n in Del Negro et al. (2017)")
xlabel!("% deviation from SS")
ylabel!("Horizon")

plot(obs_irf[6,:,4],label="IRF of policy rate to permanent safety shock",lw=4) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("IRF of policy rate to permanent safety shock\n in Del Negro et al. (2017)")
xlabel!("% deviation from SS")
ylabel!("Horizon")



var_name1 = :y_t;
shock_name1 = :b_liqtil_sh;
var_value = 10.
states_irf, obs_irf1, pseudo_irf = impulse_responses(m, system, horizon, shock_name1 , var_name1, var_value)

plot(obs_irf1[1,:,2],label="IRF of gdp per capita to permanent liquidity shock",lw=4,lc=:red) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("IRF of output gap to permanent liquidity shock\n in Del Negro et al. (2018)")
xlabel!("% deviation from SS")
ylabel!("Horizon")


plot(obs_irf1[6,:,2],label="IRF of policy rate  to permanent liquidity shock",lw=4,lc=:red) # This is the permanent preference shock's impact on output - states and observables are in declaration order.
plot!(zeros(horizon,1),lc=:black,lw=2,label="")
plot!(size=(600,600))
title!("IRF of policy rate  to permanent liquidity shock \n in Del Negro et al. (2018)")
xlabel!("% deviation from SS")
ylabel!("Horizon")


system = compute_system(m)

states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon*4)

plot(1:horizon*4,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]]  ],title="IRF to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate"])
plot!(legend=:bottomright)

savefig( "IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic


plot(1:horizon*4,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_shl6]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_shl6]]  ],title="IRF to 6 period-ahead forward guidance shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate"])
plot!(legend=:bottomright)
savefig( "IRF_to_FG6_shock.pdf")   # saves the plot from p as a .pdf vector graphic