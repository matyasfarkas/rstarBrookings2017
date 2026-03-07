using CSV
using DataFrames
using Dates
using Plots
using Measures

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------
saveroot = dirname(@__FILE__())

path = joinpath(saveroot, "output_data", "m1010", "ss24", "forecast", "tables", "hist_Forward5YearRealNaturalRate_cond=none_para=full_vint=250115.csv")

df = DataFrame(CSV.File(
    path;
    delim = ',',
    normalizenames = false,
    types = Dict("date" => Date),
    dateformat = dateformat"yyyy-mm-dd"
))

println(names(df))

rename!(df, Dict(
    "date"     => :date,
    "68.0% LB" => :lb68,
    "95.0% LB" => :lb95,
    "95.0% UB" => :ub95,
    "68.0% UB" => :ub68,
    "mean"     => :mean
))

# Keep only 1960–2016 if you want to match the original figure
# df = filter(:date => d -> d <=/.ate(2016, 12, 31), df)

# ------------------------------------------------------------------------------
# Convert Date to decimal year
# ------------------------------------------------------------------------------

decimal_year(d::Date) = year(d) + (dayofyear(d) - 1) / 365.25
df.year = decimal_year.(df.date)

# ------------------------------------------------------------------------------
# Plot styling
# ------------------------------------------------------------------------------

gr()

default(
    fontfamily = "Times New Roman",
    legend = false,
    grid = false,
    foreground_color = :black,
    guidefontsize = 20,
    tickfontsize = 15,
    titlefontsize = 22
)

bg = RGB(0.93, 0.93, 0.93)
outer_band = RGB(0.80, 0.80, 0.80)
inner_band = RGB(0.70, 0.70, 0.70)

p = plot(
    size = (1920, 1080),
    framestyle = :box,
    background_color = bg,
    background_color_inside = bg,
    foreground_color_border = :black,
    foreground_color_axis = :black,
    foreground_color_text = :black,
    left_margin = 10mm,
    right_margin = 10mm,
    bottom_margin = 10mm,
    top_margin = 26mm,
    xlims = (1999.0, 2025.0),
    ylims = (-2.5, 2.5),
    xticks = (collect(1999:1:2025), string.(1999:1:2025)),
    yticks = -2.5:0.5:2.5,
    xlabel = "Year",
    ylabel = "%",
    # title = "r*_t"
)
plot!(
    p,
    df.year,
    zeros(size(df.year));
    color = :black,
    lw = 1,
)

# 95% band
plot!(
    p,
    df.year,
    df.ub95;
    fillrange = df.lb95,
    fillcolor = outer_band,
    fillalpha = 1.0,
    linealpha = 0
)

# 68% band
plot!(
    p,
    df.year,
    df.ub68;
    fillrange = df.lb68,
    fillcolor = inner_band,
    fillalpha = 1.0,
    linealpha = 0
)

# Mean line
plot!(
    p,
    df.year,
    df.mean;
    color = :black,
    lw = 2.5,
    linestyle = :dash
)

display(p)
savepath_png = joinpath(saveroot, "Final Paper", "Figures", "EArstar_wCB_replication.png")
savefig(p, savepath_png)
savepath_pdf = joinpath(saveroot, "Final Paper", "Figures", "EArstar_wCB_replication.pdf")
savefig(p, savepath_pdf)