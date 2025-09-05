using DSGE;
using Plots # no need for `using Plots` as that is reexported here

# DSGE.Settings for data, paths, etc.
mypath = @__DIR__
idx = findlast(c -> c == '\\', mypath)
basepath = mypath[1:idx]
dataroot = joinpath(basepath, "dsge", "output_data")
saveroot = joinpath(basepath, "dsge")

horizon  = 40



# Standard model
m10 = Model1010("ss20");
mode_file = joinpath(dataroot, "m1010","ss20","estimate","raw", "paramsmode_vint=161223.h5")
specify_mode!(m10, mode_file)
system10 = compute_system(m10)
states_irf10, obs_irf10, pseudo_irf10 = impulse_responses(system10, horizon)


# EA model
m11 = Model1010("ss20");
use_FG_in_EA = true
if use_FG_in_EA
mode_file = joinpath(dataroot, "m1010","ss20","estimate","raw", "paramsmode_vint=250113.h5")
else
mode_file =joinpath(dataroot, "m1010","ss20","estimate","raw", "paramsmode_vint=250114.h5")
end
specify_mode!(m11, mode_file)
system11 = DSGE.compute_system(m11)
states_irf11, obs_irf11, pseudo_irf11 = impulse_responses(system11, horizon)



plot(1:horizon,[pseudo_irf10[m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:b_liqp_sh]], pseudo_irf11[m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:b_liqp_sh]] ],title="IRFs of r* to a permanent liquidity shock", label=["US model" "EA model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[pseudo_irf10[m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:b_safep_sh]],pseudo_irf11[m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:b_safep_sh]]  ],title="IRFs of r*  to a permanent safety shock", label=["US model" "EA model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:zp_sh ]],pseudo_irf11[m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:zp_sh]]  ],title="IRFs of r* to permanent technology shock", label=["US model" "EA model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:z_sh ]],pseudo_irf11[m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:z_sh]]  ],title="IRFs of r* to transitory technology shock", label=["US model" "EA model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_transitory_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:μ_sh ]],pseudo_irf11[m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:μ_sh]]  ],title="IRFs of r* to investment specific technology shock", label=["US model" "EA model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_investment_specific_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:σ_ω_sh ]],pseudo_irf11[m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:σ_ω_sh]]  ],title="IRFs of r* to risk shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_risk_shock.pdf")   # saves the plot from p as a .pdf vector graphic





p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_liqp_sh]]],title="Output", label=["US model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,[obs_irf10[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:b_liqp_sh]]] ,title="Inflation", label=["US model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh ]] ] ,title="Policy rate", label=["US model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh ]] ] ,title="r*", label=["US model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
savefig( "irf/all/IRF_output_inflation_FFR_to_Liquidity.pdf")   # saves the plot from p as a .pdf vector graphic


p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:b_safep_sh]]],title="Output", label=["US model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,[obs_irf10[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:b_safep_sh]]] ,title="Inflation", label=["US model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh ]] ] ,title="Policy rate", label=["US model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh ]] ] ,title="r*", label=["US model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
savefig( "irf/all/IRF_output_inflation_FFR_to_safety.pdf")   # saves the plot from p as a .pdf vector graphic


p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:zp_sh]]],title="Output", label=["US model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,[obs_irf10[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:zp_sh]]] ,title="Inflation", label=["US model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] ] ,title="Policy rate", label=["US model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]] ] ,title="r*", label=["US model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
savefig( "irf/all/IRF_output_inflation_FFR_to_permanentTFP.pdf")   # saves the plot from p as a .pdf vector graphic



p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:rm_shl6]]],title="Output", label=["US model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,[obs_irf10[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_shl6]]] ,title="Inflation", label=["US model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_shl6 ]] ] ,title="Policy rate", label=["US model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_shl6 ]] ] ,title="r*", label=["US model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
savefig( "irf/all/IRF_output_inflation_FFR_to_FG6.pdf")   # saves the plot from p as a .pdf vector graphic

p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:rm_sh]]],title="Output", label=["US model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,[obs_irf10[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_sh]]] ,title="Inflation", label=["US model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_sh ]] ] ,title="Policy rate", label=["US model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_sh ]] ] ,title="r*", label=["US model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
savefig( "irf/all/IRF_output_inflation_FFR_to_MP.pdf")   # saves the plot from p as a .pdf vector graphic

p1 = plot(1:horizon,[states_irf10[m.endogenous_states[:y_t],:, m.exogenous_shocks[:g_sh]]],title="Output", label=["US model"])
plot!(legend=:bottomright)
p2 = plot(1:horizon,[obs_irf10[m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:g_sh]]] ,title="Inflation", label=["US model"])
plot!(legend=:bottomright)
p3 = plot(1:horizon,[ obs_irf10[m.observables[:obs_nominalrate],:, m.exogenous_shocks[:g_sh ]] ] ,title="Policy rate", label=["US model" ])
plot!(legend=:bottomright)
p4=  plot(1:horizon,[ pseudo_irf10[m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:g_sh ]] ] ,title="r*", label=["US model" ])
plot(p1, p2, p3, p4, layout=(2,2), legend=false)
savefig( "irf/all/IRF_output_inflation_FFR_to_Gshock.pdf")   # saves the plot from p as a .pdf vector graphic

# # policy rate
# plot(1:horizon,[obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:b_liqp_sh]], obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:b_liqp_sh]] ],title="IRFs of the policy rate to permanent liquidity shock", label=["Basline model" "No FG model"])
# plot!(legend=:bottomright)
# savefig( "irf/all/IRF_robs_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:b_safep_sh]], obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:b_safep_sh]] ],title="IRFs of  the policy rate to permanent safety shock", label=["Basline model" "No FG model"])
# plot!(legend=:bottomright)
# savefig( "irf/all/IRF_robs_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:zp_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:zp_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:zp_sh ]] ],title="IRFs of  the policy rateto permanent technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
# plot!(legend=:bottomright)
# savefig( "irf/all/IRF_robs_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:z_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:z_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:z_sh ]] ],title="IRFs of  the policy rateto transitory technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
# plot!(legend=:bottomright)
# savefig( "irf/all/IRF_robs_to_transitory_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:μ_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:μ_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:μ_sh ]] ],title="IRFs of  the policy rateto investment specific technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
# plot!(legend=:bottomright)
# savefig( "irf/all/IRF_robs_to_investment_specific_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:σ_ω_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:σ_ω_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:σ_ω_sh ]] ],title="IRFs of  the policy rateto risk shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
# plot!(legend=:bottomright)
# savefig( "irf/all/IRF_robs_to_risk_shock.pdf")   # saves the plot from p as a .pdf vector graphic



# obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]]
# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic


# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_shl6]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_shl6]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_shl6]] ],title="IRFs to 6 period-ahead forward guidance shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate"  "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_FG6_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# using DSGE;
# path = dirname(@__FILE__)
# m = Model1011("ss20");

# system = compute_system(m)
# horizon =80
# states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# using DSGE, Plots;
# path = dirname(@__FILE__)
# m = Model1012("ss20");

# system = compute_system(m)
# horizon =80
# states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

# # plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# # plot!(legend=:bottomright)
# # savefig( "irf/m1012/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# # plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# # plot!(legend=:bottomright)
# # savefig( "irf/m1012/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1012/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1012/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1012/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1012/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic


# horizon =80
# states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" "pi:GDP Deflator"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_inflation_target_shock.pdf")   # saves the plot from p as a .pdf vector graphic
# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_shl6]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_shl6]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_shl6]],obs_irf[ m.observables[:obs_gdpdeflator],:, m.exogenous_shocks[:rm_shl6]]  ],title="IRFs to FG6 shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" "pi:GDP Deflator"])
# plot!(legend=:bottomright)


