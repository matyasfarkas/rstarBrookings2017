using Dates, CSV, DataFrames, HTTP, Statistics


"""
Download multiple FRED series as a single CSV using fredgraph.csv.
No API key required.
"""
function fredgraph_csv(ids::Vector{String})
    url = "https://fred.stlouisfed.org/graph/fredgraph.csv?id=" * join(ids, ",")
    resp = HTTP.get(url; readtimeout=60)
    return String(resp.body)
end

"""
Parse FRED CSV where missing values are ".".
Returns DataFrame with Date column and Float64 columns (missing allowed).
"""
function parse_fred_csv(csv_str::String)
    io = IOBuffer(csv_str)
    df = DataFrame(CSV.File(IOBuffer(csv_str)))
    rename!(df, Symbol("observation_date") => :date)
    df.date = Date.(df.date)
for c in names(df)
    c == :date && continue
    df[!, c] = map(df[!, c]) do x
        if ismissing(x) || x == "."
            missing
        elseif x isa Number
            Float64(x)
        else
            parse(Float64, String(x))
        end
    end
end

    return df
end

"""
Given yields (annual percent) at maturities in years, build log discount factors
using continuous-compounding approximation:
D(T)=exp(-y_c*T), y_c=log(1+y/100)
Interpolate linearly in ln D.
Return lnD(t) for requested t-grid.
"""
function ln_discount_curve(yields_pct::Vector{Float64}, mats::Vector{Float64}, tgrid::Vector{Float64})
    # include t=0
    mats2 = vcat(0.0, mats)
    lnD2  = vcat(0.0, [-log(1 + y/100)*T for (y,T) in zip(yields_pct, mats)])
    # linear interpolation in lnD
    out = similar(tgrid)
    for (i,t) in pairs(tgrid)
        if t ≤ mats2[1]
            out[i] = lnD2[1]
        elseif t ≥ mats2[end]
            out[i] = lnD2[end]
        else
            j = findlast(m -> m ≤ t, mats2)
            t0,t1 = mats2[j], mats2[j+1]
            y0,y1 = lnD2[j], lnD2[j+1]
            w = (t - t0)/(t1 - t0)
            out[i] = (1-w)*y0 + w*y1
        end
    end
    return out
end

"""
Compute average forward rate over each quarter interval using ln discounts.
For interval [t0,t1], continuous forward f = (lnD(t0)-lnD(t1))/(t1-t0).
Convert to annual simple percent: (exp(f)-1)*100.
Then to percent-per-quarter: /4.
"""
function quarterly_forward_rates(lnD::Vector{Float64}, tgrid::Vector{Float64})
    @assert length(lnD) == length(tgrid)
    rates = Float64[]
    for h in 1:6
        t0 = (h-1)*0.25
        t1 = h*0.25
        i0 = findfirst(==(t0), tgrid)
        i1 = findfirst(==(t1), tgrid)
        f  = (lnD[i0] - lnD[i1])/(t1 - t0)              # cont. forward
        ann_pct = (exp(f) - 1)*100                      # annual % (simple)
        push!(rates, ann_pct/4)                         # % per quarter
    end
    return rates
end

# --- MAIN ---

ids = ["DGS1MO","DGS3MO","DGS6MO","DGS1","DGS2"]  # public H.15 series
raw = fredgraph_csv(ids)
df  = parse_fred_csv(raw)

# maturities in years corresponding to the ids above
mats = [1/12, 0.25, 0.5, 1.0, 2.0]

# grid of quarter endpoints out to 6 quarters
tgrid = collect(0.0:0.25:1.5)

# daily forward rates (6 horizons)
for h in 1:6
    df[!, Symbol("obs_nominalrate$(h)")] = Vector{Union{Missing,Float64}}(missing, nrow(df))
end

for i in 1:nrow(df)
    y = [df[i, Symbol(id)] for id in ids]
    if any(ismissing, y)
        continue
    end
    lnD = ln_discount_curve(collect(Float64, y), mats, tgrid)
    qr  = quarterly_forward_rates(lnD, tgrid)
    for h in 1:6
        df[i, Symbol("obs_nominalrate$(h)")] = qr[h]
    end
end

# Convert to quarterly frequency: take last available observation in each quarter
df_q = sort(df, :date)
df_q[!, :qend] = Date.(year.(df_q.date), month.(df_q.date) .+ (3 .- ((month.(df_q.date).-1) .% 3 .+ 1)), 1) .+ Month(1) .- Day(1)
g = groupby(df_q, :qend)
out = combine(g) do sdf
    sdf[end, [:qend, Symbol("obs_nominalrate1"), Symbol("obs_nominalrate2"), Symbol("obs_nominalrate3"),
                   Symbol("obs_nominalrate4"), Symbol("obs_nominalrate5"), Symbol("obs_nominalrate6")]]
end
rename!(out, :qend => :date)
sort!(out, :date)

CSV.write("yield_curve_dataset_treasury_from_fred.csv", out)
