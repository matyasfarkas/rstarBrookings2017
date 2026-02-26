using Dates, CSV, DataFrames, HTTP

# -------------------------
# FRED download helpers
# -------------------------
function fredgraph_csv(ids::Vector{String})
    url = "https://fred.stlouisfed.org/graph/fredgraph.csv?id=" * join(ids, ",")
    resp = HTTP.get(url; readtimeout=60)
    return String(resp.body)
end

function parse_fred_csv(csv_str::String)
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
                y = tryparse(Float64, String(x))
                y === nothing ? missing : y
            end
        end
    end
    return df
end

# -------------------------
# Curve + forward rates
# -------------------------
function ln_discount_curve(yields_pct::Vector{Float64}, mats::Vector{Float64}, tgrid::Vector{Float64})
    mats2 = vcat(0.0, mats)
    lnD2  = vcat(0.0, [-log(1 + y/100) * T for (y, T) in zip(yields_pct, mats)])

    out = similar(tgrid)
    for (i, t) in pairs(tgrid)
        if t <= mats2[1]
            out[i] = lnD2[1]
        elseif t >= mats2[end]
            out[i] = lnD2[end]
        else
            j = findlast(m -> m <= t, mats2)
            t0, t1 = mats2[j], mats2[j+1]
            y0, y1 = lnD2[j], lnD2[j+1]
            w = (t - t0) / (t1 - t0)
            out[i] = (1 - w) * y0 + w * y1
        end
    end
    return out
end

function quarterly_forward_rates_percent_per_quarter(lnD::Vector{Float64})
    @assert length(lnD) == 7  # tgrid = 0:0.25:1.5
    rates = Vector{Float64}(undef, 6)
    for h in 1:6
        i0 = h
        i1 = h + 1
        f  = (lnD[i0] - lnD[i1]) / 0.25
        ann_pct = (exp(f) - 1) * 100
        rates[h] = ann_pct / 4
    end
    return rates
end

# -------------------------
# Ragged-edge mask logic
# -------------------------
"""
Apply your desired ragged-edge structure to ant1..ant6.

Rules:
- Keep only dates from 2008-12-31 onward.
- 2008Q4–2013Q4: keep ant1..ant6
- 2014Q1: ant6 missing
- 2014Q2: ant5..ant6 missing
- 2014Q3: ant4..ant6 missing
- 2014Q4: ant3..ant6 missing
- 2015Q1: ant2..ant6 missing
- 2015Q2–2019Q4: all ant1..ant6 missing
- 2020Q1–2025Q2: keep ant1..ant6 (full path)
- Outside these windows: drop rows (so dataset ends 2025Q2)
"""
function apply_ragged_edges!(out::DataFrame)
    # restrict to [2008Q4, 2025Q2]
    startd = Date(2008, 12, 31)
    endd   = Date(2025, 6, 30)
    filter!(row -> (row.date >= startd) && (row.date <= endd), out)

    # helper: set a range of ants missing
    function set_missing!(df, d::Date, ants::Vector{Int})
        idx = findall(df.date .== d)
        isempty(idx) && return
        for h in ants
            df[idx, Symbol("ant$(h)")] .= missing
        end
    end

    # 2014 raggeding
    set_missing!(out, Date(2014,3,31),  [6])
    set_missing!(out, Date(2014,6,30),  [5,6])
    set_missing!(out, Date(2014,9,30),  [4,5,6])
    set_missing!(out, Date(2014,12,31), [3,4,5,6])

    # 2015Q1 only ant1
    set_missing!(out, Date(2015,3,31),  [2,3,4,5,6])

    # 2015Q2–2019Q4 all missing
    for d in out.date
        if d >= Date(2015,6,30) && d <= Date(2019,12,31)
            for h in 1:6
                out[out.date .== d, Symbol("ant$(h)")] .= missing
            end
        end
    end

    # 2020Q1–2025Q2 full path: do nothing (keep whatever computed)
    # But ensure those rows exist; if FRED ends earlier, they will be missing anyway.

    return out
end

# -------------------------
# MAIN
# -------------------------
ids  = ["DGS1MO","DGS3MO","DGS6MO","DGS1","DGS2"]
mats = [1/12,     0.25,    0.5,     1.0,   2.0]

raw = fredgraph_csv(ids)
df  = parse_fred_csv(raw)

tgrid = collect(0.0:0.25:1.5)

# compute daily ant1..ant6
for h in 1:6
    df[!, Symbol("ant$(h)")] = Vector{Union{Missing,Float64}}(missing, nrow(df))
end

for i in 1:nrow(df)
    y = [df[i, Symbol(id)] for id in ids]
    if any(ismissing, y)
        continue
    end
    lnD = ln_discount_curve(Float64.(y), mats, tgrid)
    qr  = quarterly_forward_rates_percent_per_quarter(lnD)
    for h in 1:6
        df[i, Symbol("ant$(h)")] = qr[h]
    end
end

# -------------------------------------------------------------------
# Quarterly aggregation: pick LAST day in quarter with complete ant1..ant6
# -------------------------------------------------------------------
df_q = sort(df, :date)
df_q.qend = lastdayofquarter.(df_q.date)

antcols = Symbol.("ant" .* string.(1:6))

# A row is "valid" if ant1..ant6 are all non-missing
valid = completecases(df_q[:, antcols])

# For each quarter, find the last index with valid=true
g = groupby(df_q, :qend)

last_valid_idx = Vector{Union{Missing,Int}}(missing, length(g))
for (k, sdf) in enumerate(g)
    # sdf is a SubDataFrame; use its row indices in the parent
    parent_rows = parentindices(sdf)[1]   # vector of row indices in df_q
    ok = valid[parent_rows]
    if any(ok)
        last_valid_idx[k] = parent_rows[findlast(ok)]
    else
        last_valid_idx[k] = missing
    end
end

# Build output quarter list (one row per quarter)
qends = [keys(g)[k].qend for k in 1:length(g)]
out = DataFrame(date = qends)

# Fill ant1..ant6 from last valid row when available, else missing
for c in antcols
    out[!, c] = Vector{Union{Missing,Float64}}(missing, nrow(out))
end

for k in 1:nrow(out)
    idx = last_valid_idx[k]
    if !ismissing(idx)
        for c in antcols
            out[k, c] = df_q[idx, c]
        end
    end
end

sort!(out, :date)


# add ant7..ant13
for h in 7:13
    out[!, Symbol("ant$(h)")] = Vector{Union{Missing,Float64}}(missing, nrow(out))
end

# apply your ragged-edge structure and truncate window
apply_ragged_edges!(out)

# reorder columns
select!(out, [:date; Symbol.("ant" .* string.(1:13))...])

# format date and write with NaN tokens
out.date = Dates.format.(out.date, dateformat"dd-mm-yy")

# OUTPUT PATH (edit if needed)
#
CSV.write("C:/Mac/Home/Documents/GitHub/rstarBrookings2017/dsge/input_data/raw/ois_250825.csv", out; missingstring="NaN")
#CSV.write("yield_curve_dataset_treasury_from_fred.csv", out; missingstring="NaN")

println("Wrote yield curve file with columns date, ant1..ant13")
println("Rows: ", nrow(out), "  Cols: ", ncol(out))
