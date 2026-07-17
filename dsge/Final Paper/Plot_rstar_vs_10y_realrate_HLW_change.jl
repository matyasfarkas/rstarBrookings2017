using CSV
using DataFrames
using Dates
using Plots
using Measures

# Compare 5-year-horizon and 10-year-horizon r-star estimates.
# The 5-year estimate carries the posterior confidence bands used in the paper.

repo_root = abspath(joinpath(@__DIR__, "..", ".."))
dsge_root = joinpath(repo_root, "dsge")
figures_dir = joinpath(@__DIR__, "Figures")
mkpath(figures_dir)

us_rstar_path = joinpath(dsge_root, "output_data", "m1010", "ss20", "forecast", "tables",
                         "hist_Forward5YearRealNaturalRate_cond=none_para=full_vint=250825.csv")
us_ten_year_path = joinpath(dsge_root, "output_data", "m1010", "ss20", "forecast", "tables",
                            "hist_Forward10YearRealNaturalRate_cond=none_para=mode_vint=250825.csv")
ea_rstar_path = joinpath(dsge_root, "output_data", "m1010", "ss24", "forecast", "tables",
                         "hist_Forward5YearRealNaturalRate_cond=none_para=full_vint=250115.csv")
ea_ten_year_path = joinpath(dsge_root, "output_data", "m1010", "ss24", "forecast", "tables",
                            "hist_Forward10YearRealNaturalRate_cond=none_para=mode_vint=250115.csv")

decimal_year(d::Date) = year(d) + (dayofyear(d) - 1) / 365.25

function parse_date_value(x)
    x isa Date && return x
    s = String(x)
    try
        return Date(s, dateformat"yyyy-mm-dd")
    catch
        return Date(s, dateformat"mm/dd/yyyy")
    end
end

function read_table(path::AbstractString)
    isfile(path) || error("Missing input file: $(path)")
    return DataFrame(CSV.File(path; normalizenames = false))
end

function load_rstar_bands(path::AbstractString)
    df = read_table(path)
    rename!(df, Dict(
        "date"     => :Date,
        "68.0% LB" => :RStar5Y_LB68,
        "95.0% LB" => :RStar5Y_LB95,
        "95.0% UB" => :RStar5Y_UB95,
        "68.0% UB" => :RStar5Y_UB68,
        "mean"     => :RStar5Y_Mean
    ))
    df.Date = parse_date_value.(df.Date)
    df.Year = decimal_year.(df.Date)
    return df
end

function load_mode_series(path::AbstractString, outname::Symbol)
    df = read_table(path)
    rename!(df, Dict("date" => :Date, "mean" => outname))
    df.Date = parse_date_value.(df.Date)
    return df[:, [:Date, outname]]
end

function attach_ten_year_rstar(rstar_df::DataFrame, ten_year_df::DataFrame)
    ten_year_lookup = Dict{Date, Float64}()
    for i in 1:nrow(ten_year_df)
        ten_year_lookup[ten_year_df.Date[i]] = ten_year_df.RStar10Y_Mean[i]
    end

    keep = [haskey(ten_year_lookup, d) for d in rstar_df.Date]
    out = rstar_df[keep, :]
    out.RStar10Y_Mean = [ten_year_lookup[d] for d in out.Date]
    return out
end

function add_country_column(df::DataFrame, country::AbstractString)
    out = copy(df)
    out.Country = fill(country, nrow(out))
    return out
end

function rstar_panel(df::DataFrame, title_text::AbstractString, xlims_tuple, xtick_years, ylims_tuple, legend_pos)
    bg = RGB(0.93, 0.93, 0.93)
    outer_band = RGB(0.80, 0.80, 0.80)
    inner_band = RGB(0.70, 0.70, 0.70)

    five_label = "r* - E_t[r*_{t+20}] - 5-year horizon r*"
    ten_label = "E_t[r*_{t+40}] - 10-year horizon r*"

    p = plot(
        framestyle = :box,
        background_color = bg,
        background_color_inside = bg,
        foreground_color_border = :black,
        foreground_color_axis = :black,
        foreground_color_text = :black,
        xlims = xlims_tuple,
        ylims = ylims_tuple,
        xticks = (xtick_years, string.(xtick_years)),
        xlabel = "Year",
        ylabel = "%",
        title = title_text,
        legend = legend_pos
    )

    plot!(p, df.Year, zeros(nrow(df)); color = :black, lw = 1, label = "")
    plot!(p, df.Year, df.RStar5Y_UB95; fillrange = df.RStar5Y_LB95,
          fillcolor = outer_band, fillalpha = 1.0, linealpha = 0, label = "")
    plot!(p, df.Year, df.RStar5Y_UB68; fillrange = df.RStar5Y_LB68,
          fillcolor = inner_band, fillalpha = 1.0, linealpha = 0, label = "")
    plot!(p, df.Year, df.RStar5Y_Mean; color = :black, lw = 2.5,
          linestyle = :dash, label = five_label)
    plot!(p, df.Year, df.RStar10Y_Mean; color = :red, lw = 2.5,
          label = ten_label)
    return p
end

gr()

default(
    fontfamily = "Times New Roman",
    grid = false,
    guidefontsize = 18,
    tickfontsize = 13,
    titlefontsize = 20,
    legendfontsize = 11,
    left_margin = 8mm,
    right_margin = 5mm,
    bottom_margin = 7mm,
    top_margin = 7mm
)

us_levels = attach_ten_year_rstar(load_rstar_bands(us_rstar_path),
                                  load_mode_series(us_ten_year_path, :RStar10Y_Mean))
ea_levels = attach_ten_year_rstar(load_rstar_bands(ea_rstar_path),
                                  load_mode_series(ea_ten_year_path, :RStar10Y_Mean))

us_plot = rstar_panel(us_levels, "United States: 5-year and 10-year horizon r-star",
                      (2008.0, 2025.5), collect(2008:2:2025), (-1.5, 3.5), (0.42, 0.78))
ea_plot = rstar_panel(ea_levels, "Euro area: 5-year and 10-year horizon r-star",
                      (1999.0, 2025.5), collect(1999:2:2025), (-2.5, 2.5), (0.42, 0.78))

combined_plot = plot(
    us_plot,
    ea_plot;
    layout = (2, 1),
    size = (1500, 1100)
)

combined_pdf = joinpath(figures_dir, "rstar_5y_vs_10y_horizon_US_EA.pdf")
combined_png = joinpath(figures_dir, "rstar_5y_vs_10y_horizon_US_EA.png")
legacy_current_pdf = joinpath(figures_dir, "rstar_current_5y_10y_horizon_US_EA.pdf")
legacy_current_png = joinpath(figures_dir, "rstar_current_5y_10y_horizon_US_EA.png")
legacy_hlw_pdf = joinpath(figures_dir, "rstar_5y_vs_10y_realrate_HLW_change_US_EA.pdf")
legacy_hlw_png = joinpath(figures_dir, "rstar_5y_vs_10y_realrate_HLW_change_US_EA.png")
savefig(combined_plot, combined_pdf)
savefig(combined_plot, combined_png)
savefig(combined_plot, legacy_current_pdf)
savefig(combined_plot, legacy_current_png)
savefig(combined_plot, legacy_hlw_pdf)
savefig(combined_plot, legacy_hlw_png)

level_export = vcat(
    add_country_column(us_levels[:, [:Date, :RStar5Y_LB68, :RStar5Y_LB95, :RStar5Y_UB95,
                                     :RStar5Y_UB68, :RStar5Y_Mean, :RStar10Y_Mean]], "US"),
    add_country_column(ea_levels[:, [:Date, :RStar5Y_LB68, :RStar5Y_LB95, :RStar5Y_UB95,
                                     :RStar5Y_UB68, :RStar5Y_Mean, :RStar10Y_Mean]], "EA")
)

CSV.write(joinpath(figures_dir, "rstar_5y_vs_10y_horizon_US_EA.csv"), level_export)
CSV.write(joinpath(figures_dir, "rstar_current_5y_10y_horizon_US_EA.csv"), level_export)

println("Saved: $(combined_pdf)")
println("Saved: $(combined_png)")
println("Saved level export in: $(figures_dir)")