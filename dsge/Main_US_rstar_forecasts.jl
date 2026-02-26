using DSGE, ClusterManagers, HDF5, Plots, StatsPlots, CSV, DataFrames, Dates,FileIO

# Initialize model object
# Note that the default for m1010 uses 6 anticipated shocks
m = Model1010("ss20")

dataroot = joinpath(dirname(@__FILE__()), "input_data")
saveroot = dirname(@__FILE__())
m <= DSGE.Setting(:dataroot, dataroot, "Input data directory path")
m <= DSGE.Setting(:saveroot, saveroot, "Output data directory path")
m <= DSGE.Setting(:data_vintage, "250825")
m <= DSGE.Setting(:use_population_forecast, false)
m <= DSGE.Setting(:date_forecast_start,  quartertodate("2025-Q2"))
m <= DSGE.Setting(:date_conditional_end, quartertodate("2025-Q2"))

params_mode = load_draws(m, :mode)
DSGE.update!(m, params_mode)
DSGE.steadystate!(m)

## THIS SECTION UPDATES THE DATA VINTAGE IF NEEDED
# m <= DSGE.Setting(:data_vintage, "250825")
# m <= DSGE.Setting(:use_population_forecast, false)
# params_mode = load_draws(m, :mode)
# DSGE.update!(m, params_mode)

# # Settings for estimation
# # set to false => will load pre-computed mode and hessian before MCMC
# m <= DSGE.Setting(:reoptimize, true)
# m <= DSGE.Setting(:calculate_hessian, true)

# # Settings for forecast dates
# m <= DSGE.Setting(:date_forecast_start,  quartertodate("2024-Q4"))
# m <= DSGE.Setting(:date_conditional_end, quartertodate("2024-Q4"))
#                 m <= DSGE.Setting(:date_forecast_start, quartertodate("2025-Q1"))
#                 df_final = load_data(m, try_disk = false, check_empty_columns = false, summary_statistics = :none)
#                 df_yc = DataFrame(CSV.File("C:/Mac/Home/Documents/GitHub/rstarBrookings2017/dsge/input_data/raw/ois_250825.csv", missingstring="NaN"))
#                 df_final[198:end,end-5:end] = df_yc[1:65,2:7]

df_final = load_data(m; check_empty_columns = false)
maxhoriz = 20
full_start_date = df_final[1, :date]
full_end_date   = DSGE.iterate_quarters(df_final[end, :date], maxhoriz)
all_dates       = DSGE.get_quarter_ends(full_start_date, full_end_date)
output_vars = [:histobs, :forecastobs, :histpseudo,:forecastpseudo]

# Initialize the Master DataFrame
master_df = DataFrame(date = all_dates)

# Get variable index for r-star
outvar = :Forward5YearRealNaturalRate
var_idx = m.pseudo_observables[outvar]


mnc = m

for i_year in 2001:2025
    for i_quarter in 1:4
        
        # Timing & Slicing
        date_last_obs = quartertodate(string(i_year)*"-Q"*string(i_quarter))
        if date_last_obs > df_final[end, :date]
            continue
        end

        # Define Vintage Label (e.g., "01Q1")
        vint_label = string(i_year)[3:4] * "Q" * string(i_quarter)
        col_name   = Symbol("vint_" * vint_label)

        # Update Model Settings
        mnc <= DSGE.Setting(:date_mainsample_end, date_last_obs)
        mnc <= DSGE.Setting(:date_forecast_start, DSGE.iterate_quarters(date_last_obs, 1))
        mnc <= DSGE.Setting(:date_conditional_end, date_last_obs)

        # Slice Data
        idx = findfirst(d -> d == date_last_obs, df_final[:, :date])
        df_slice = df_final[1:idx, :]

        # --- 2. RUN FORECAST ---
        forecast_one(mnc, :mode, :none, output_vars; 
                     df = df_slice, check_empty_columns = false, verbose = :low)
        
        # --- 3. EXTRACT AND ALIGN DATA ---
        out_files = get_forecast_output_files(mnc, :mode, :none, output_vars)
        hist_p    = FileIO.load(out_files[:histpseudo], "arr")
        fore_p    = FileIO.load(out_files[:forecastpseudo], "arr")

        h_vals = hist_p[var_idx, :]
        f_vals = fore_p[var_idx, 1:min(size(fore_p, 2), maxhoriz)]

        # Determine exact dates for the values provided by the model
        # History dates: Count back from date_last_obs based on result length
        h_dates = [DSGE.iterate_quarters(date_last_obs, -i) for i in (length(h_vals)-1):-1:0]
        
        # Forecast dates: Start from the quarter after last_obs
        f_dates = [DSGE.iterate_quarters(date_last_obs, i) for i in 1:length(f_vals)]

        # Concatenate dates and values
        combined_dates = vcat(h_dates, f_dates)
        combined_vals  = vcat(Vector{Union{Float64, Missing}}(h_vals), 
                              Vector{Union{Float64, Missing}}(f_vals))

        # --- 4. FIXING DIMENSION MISMATCH & JOINING ---
        # Initialize a temp dataframe with the full timeline as Missing
        temp_vint_df = DataFrame(date = all_dates)
        temp_vint_df[!, :value] = Vector{Union{Float64, Missing}}(missing, length(all_dates))

        # Create a small dataframe of current results and update the temp one
        current_results = DataFrame(date = combined_dates, result_val = combined_vals)
        
        # Left join current results onto the full timeline to ensure alignment
        temp_vint_df = join(temp_vint_df, current_results, on = :date, makeunique=true, kind = :left)
        
        # Move the values over and clean up
        temp_vint_df[!, :value] = temp_vint_df.result_val
        select!(temp_vint_df, [:date, :value])

        # Join into the Master DataFrame
        master_df = join(master_df, temp_vint_df, on = :date, makeunique=true, kind = :left)
        rename!(master_df, :value => col_name)

        println("Processed Vintage: $vint_label")
    end
end

# --- 5. FINAL EXPORT ---
output_path = joinpath(saveroot, "rstar_pseudo_realtime_results.csv")
CSV.write(output_path, master_df)
println("Process complete. File saved to $output_path")

# 1. Define zoom window and specific Q4 tick dates
zoom_start = Date(2008, 1, 1)
zoom_end   = Date(2026, 12, 31)

# Generate a list of dates for every year's Q4 (Dec 31st)
tick_dates = [Date(y, 12, 31) for y in 2008:2026]
# Create the corresponding labels: ["2008Q4", "2009Q4", ...]
tick_labels = [string(y) * "Q4" for y in 2008:2026]

# 2. Initialize the plot
p = plot(title = "",
         xlabel = "Year", 
         ylabel = "Natural Rate of Interest (r*), (%)",
         legend = :topright,
         size = (1000, 600),
         grid = :both,
         minorgrid = true)

# 3. Identify vintage columns
vint_cols = DataFrames.filter(col -> startswith(string(col), "vint_"), names(master_df))

# 4. Plot background vintages (Red, transparent)
for col in vint_cols[1:end-1]
    if !all(ismissing.(master_df[!, col]))
        plot!(p, master_df.date, master_df[!, col]*4, # Annualize the quarterly values for better visual comparison
              color = :red, linealpha = 0.5, linewidth = 0.5, label = "") 
    end
end

# 5. Highlight the LATEST vintage (Black, solid)
latest_vint = vint_cols[end]
plot!(p, master_df.date, master_df[!, latest_vint]*4, 
      color = :black, linewidth = 2.5, label = "")

# 6. ADD ZERO LINE
# We use hline! for a horizontal line. color=:black and linestyle=:dash for clarity.
hline!(p, [0.0], color = :black,  linewidth = 1, label = "")

# 7. FIX X-AXIS: Limits and Manual Labels
# We pass a tuple of (Numeric_Positions, String_Labels) to xticks!
# We must use Dates.value() for positions to match the coordinate system
xlims!(p, (Dates.value(zoom_start), Dates.value(zoom_end)))
xticks!(p, (Dates.value.(tick_dates), tick_labels), xrotation = 45)

# 8. Add Recession Shading
vspan!(p, [Dates.value(Date(2007, 12, 1)), Dates.value(Date(2009, 6, 30))], 
       color=:gray, alpha=0.2, label="")
vspan!(p, [Dates.value(Date(2020, 1, 1)), Dates.value(Date(2020, 6, 30))], 
       color=:gray, alpha=0.2, label="")
# vspan!(p, [Dates.value(Date(2025,1, 1)), Dates.value(Date(2026, 12, 31))], 
#        color=:blue, alpha=0.6, label="")


# 9. Save and Display
savefig(joinpath(saveroot, "US_rstar_final_chart.png"))
display(p)


# 1. Identify the Latest Vintage for the Anchor Line
latest_vint = vint_cols[end]

# Extract history from the latest vintage (The Black Body)
anchor_mask = .!ismissing.(master_df[!, latest_vint])
anchor_dates = master_df.date[anchor_mask]
anchor_vals  = master_df[anchor_mask, latest_vint] # Annualize the quarterly values for better visual comparison

# 2. Initialize the Plot
p_starfish = plot(title = "", #title = "Natural Rate (r*) Real-Time Forecasts",
         xlabel = "Year", ylabel = "Natural Rate of Interest (r*), (%)",
         legend = :none, # Legend dropped as requested
         size = (1000, 600), grid = :both)

# 3. Plot the RED Forecasts (Tethered to the Black Line)
for col in vint_cols[1:end-1]
    # Get the date this vintage 'started'
    v_str = string(col)
    v_year = 2000 + parse(Int, v_str[6:7])
    v_qrt  = parse(Int, v_str[9:9])
    v_start_date = quartertodate("$(v_year)-Q$(v_qrt)")
    
    # 3a. Find the "Tether Point": The value of the BLACK line at this start date
    tether_idx = findfirst(d -> d == v_start_date, anchor_dates)
    
    # If the anchor line doesn't exist yet for this vintage, skip it
    if tether_idx === nothing
        continue
    end
    
    tether_val = anchor_vals[tether_idx]

    # 3b. Mask the FORECAST data (Dates strictly after the start date)
    fcast_mask = (master_df.date .> v_start_date) .& (.!ismissing.(master_df[!, col]))
    
    if any(fcast_mask)
        # We prepend the tether point to the forecast series to ensure connection
        f_dates = vcat(v_start_date, master_df.date[fcast_mask])
        f_vals  = vcat(tether_val, master_df[fcast_mask, col])
        
        plot!(p_starfish, f_dates, f_vals*4, # Annualize the quarterly values for better visual comparison
              color = :red, linealpha = 0.6, linewidth = 1.0)
    end
end

# 4. Plot the Latest Vintage's OWN forecast (The Red Dash at the end)
# This connects the very end of the black line to the current future projections
last_hist_date = anchor_dates[end-20]
curr_fcast_mask = (master_df.date.> last_hist_date) 
if any(curr_fcast_mask)
    plot!(p_starfish, vcat(last_hist_date, master_df.date[curr_fcast_mask]), 
          vcat(anchor_vals[end-20], master_df[curr_fcast_mask, latest_vint])*4, # Annualize the quarterly values for better visual comparison
          color = :red, linealpha = 0.6, linewidth = 1.0)
end

# 5. Plot the SINGLE Black Anchor Line (History)
plot!(p_starfish, anchor_dates[1:end-20], anchor_vals[1:end-20]*4, # Annualize the quarterly values for better visual comparison
      color = :black, linewidth = 3.0)

# 6. Formatting (Zero line, limits, ticks)
hline!(p_starfish, [0.0], color = :black, linewidth = 1)

# 7. Shading
vspan!(p_starfish, [Dates.value(Date(2007,12,1)), Dates.value(Date(2009,6,30))], color=:gray, alpha=0.15)
vspan!(p_starfish, [Dates.value(Date(2020,1,1)), Dates.value(Date(2020,6,30))], color=:gray, alpha=0.15)
xlims!(p_starfish, (Dates.value(zoom_start), Dates.value(zoom_end)))
xticks!(p_starfish, (Dates.value.(tick_dates), tick_labels), xrotation = 45)

# 8. Save
savefig(joinpath(saveroot, "US_rstar_starfish_chart_wFG.png"))
display(p_starfish)