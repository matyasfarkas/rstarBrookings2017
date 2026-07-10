##########################################################################################
## Posterior-mode convenience yield extraction in the policy-rate domain
##########################################################################################

using DSGE, HDF5
using DataFrames, CSV, Dates

repo_root = dirname(dirname(@__DIR__))
dataroot = joinpath(repo_root, "dsge", "input_data")
saveroot = joinpath(repo_root, "dsge")
figroot = joinpath(saveroot, "Final Paper", "Figures", "CY_rate_domain")
mkpath(figroot)

smoother = :durbin_koopman
cy_state_names = [:b_liqtil_t, :b_liqp_t, :b_safetil_t, :b_safep_t]

function configure_model!(m, vintage)
    m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
    m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
    m <= DSGE.Setting(:data_vintage, vintage)
    m <= DSGE.Setting(:reoptimize, false)
    m <= DSGE.Setting(:calculate_hessian, false)
    m <= DSGE.Setting(:date_mainsample_start, quartertodate("1970-Q3"))
    m <= DSGE.Setting(:date_presample_start, quartertodate("1970-Q2"))
    m <= DSGE.Setting(:date_forecast_start, quartertodate("2024-Q3"))
    m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q3"))
    m <= DSGE.Setting(:forecast_smoother, smoother)
    return m
end

function smooth_at_mode(spec, vintage, mode_file)
    m = configure_model!(Model1010(spec), vintage)
    df = load_data(m; check_empty_columns = false)
    specify_mode!(m, mode_file)
    system = DSGE.compute_system(m)
    states, shocks, pseudo = DSGE.smooth(m, df, system; draw_states = false)
    dates = df.date[end-size(states, 2)+1:end]
    return m, dates, states, pseudo
end

function state_series(m, states_mat, sym::Symbol)
    if haskey(m.endogenous_states, sym)
        return vec(states_mat[m.endogenous_states[sym], :])
    elseif haskey(m.endogenous_states_augmented, sym)
        return vec(states_mat[m.endogenous_states_augmented[sym], :])
    else
        error("State $(sym) not found in model.")
    end
end

function cy_drift(m)
    return 100 * log(m[:lnb_liq]) + 100 * log(m[:lnb_safe])
end

function policy_rate_coefficient_in_euler(m)
    sigma_c = m[Symbol(Char(0x03c3), :_c)]
    habit_term = m[:h] * exp(-m[:z_star])
    return (1.0 - habit_term) / (sigma_c * (1.0 + habit_term))
end

cy_rate_domain_scale(m) = -1.0 / policy_rate_coefficient_in_euler(m)

function cy_components(m, states_mat)
    b_liqtil = state_series(m, states_mat, :b_liqtil_t)
    b_liqp = state_series(m, states_mat, :b_liqp_t)
    b_safetil = state_series(m, states_mat, :b_safetil_t)
    b_safep = state_series(m, states_mat, :b_safep_t)

    liquidity_state = b_liqtil .+ b_liqp
    safety_state = b_safetil .+ b_safep
    state_total = liquidity_state .+ safety_state

    liquidity_trend = fill(100 * log(m[:lnb_liq]), size(states_mat, 2))
    safety_trend = fill(100 * log(m[:lnb_safe]), size(states_mat, 2))
    trend_total = liquidity_trend .+ safety_trend

    return (
        b_liqtil = b_liqtil,
        b_liqp = b_liqp,
        liquidity_state = liquidity_state,
        liquidity_trend = liquidity_trend,
        liquidity_total = liquidity_state .+ liquidity_trend,
        b_safetil = b_safetil,
        b_safep = b_safep,
        safety_state = safety_state,
        safety_trend = safety_trend,
        safety_total = safety_state .+ safety_trend,
        state_total = state_total,
        trend_total = trend_total,
        total = state_total .+ trend_total,
    )
end

function rate_domain_dataframe(label, dates, m, states_mat)
    c = cy_components(m, states_mat)
    scale = cy_rate_domain_scale(m)
    to_rate_domain_apr(x) = 4 .* scale .* x

    return DataFrame(
        :Date => dates,
        Symbol("CYRateDomainScale_", label) => fill(scale, length(dates)),
        Symbol("LiquidityTransitory_", label, "_RateDomain_APR") => to_rate_domain_apr(c.b_liqtil),
        Symbol("LiquidityPermanent_", label, "_RateDomain_APR") => to_rate_domain_apr(c.b_liqp),
        Symbol("LiquidityStateTotal_", label, "_RateDomain_APR") => to_rate_domain_apr(c.liquidity_state),
        Symbol("LiquidityTrend_", label, "_RateDomain_APR") => to_rate_domain_apr(c.liquidity_trend),
        Symbol("LiquidityConvenienceYield_", label, "_RateDomain_APR") => to_rate_domain_apr(c.liquidity_total),
        Symbol("SafetyTransitory_", label, "_RateDomain_APR") => to_rate_domain_apr(c.b_safetil),
        Symbol("SafetyPermanent_", label, "_RateDomain_APR") => to_rate_domain_apr(c.b_safep),
        Symbol("SafetyStateTotal_", label, "_RateDomain_APR") => to_rate_domain_apr(c.safety_state),
        Symbol("SafetyTrend_", label, "_RateDomain_APR") => to_rate_domain_apr(c.safety_trend),
        Symbol("SafetyConvenienceYield_", label, "_RateDomain_APR") => to_rate_domain_apr(c.safety_total),
        Symbol("ConvenienceYieldStateTotal_", label, "_RateDomain_APR") => to_rate_domain_apr(c.state_total),
        Symbol("ConvenienceYieldTrend_", label, "_RateDomain_APR") => to_rate_domain_apr(c.trend_total),
        Symbol("ConvenienceYield_", label, "_RateDomain_APR") => to_rate_domain_apr(c.total),
    )
end

ea_mode_file = joinpath(saveroot, "output_data", "m1010", "ss24", "estimate", "raw", "paramsmode_vint=250115.h5")
us_mode_file = joinpath(saveroot, "output_data", "m1010", "ss20", "estimate", "raw", "paramsmode_vint=250825.h5")

m_EA, dates_EA, states_EA, pseudo_EA = smooth_at_mode("ss24", "250115", ea_mode_file)
m_US, dates_US, states_US, pseudo_US = smooth_at_mode("ss20", "250825", us_mode_file)

df_EA = rate_domain_dataframe("EA", dates_EA, m_EA, states_EA)
df_US = rate_domain_dataframe("US", dates_US, m_US, states_US)
df_out = join(df_EA, df_US, on = :Date, kind = :outer)

out_path = joinpath(figroot, "ConvenienceYield_fullsample_US_EA_rate_domain.csv")
CSV.write(out_path, df_out)

println("EA CY rate-domain scale (-1/Xi_r): ", cy_rate_domain_scale(m_EA))
println("US CY rate-domain scale (-1/Xi_r): ", cy_rate_domain_scale(m_US))
println("Saved ", out_path)
