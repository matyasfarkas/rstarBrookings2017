using DSGE;
using Plots # no need for `using Plots` as that is reexported here

path = dirname(@__FILE__)
m = Model1010("ss20");

system = compute_system(m)
horizon =80
states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1010/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1010/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_shl6]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_shl6]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_shl6]] ],title="IRFs to 6 period-ahead forward guidance shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate"  "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1010/IRF_to_FG6_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1010/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1010/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

using DSGE;
path = dirname(@__FILE__)
m = Model1011("ss20");

system = compute_system(m)
horizon =80
states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1011/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1011/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1011/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1011/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1011/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1011/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

using DSGE;
path = dirname(@__FILE__)
m = Model1012("ss20");

system = compute_system(m)
horizon =80
states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

