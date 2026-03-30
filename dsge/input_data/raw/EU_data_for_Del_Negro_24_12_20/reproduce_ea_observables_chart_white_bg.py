from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

csv_path = Path("data_dsid=04_vint=250115.csv")
output_path = Path("ea_observables_white_bg.png")

df = pd.read_csv(csv_path)
df["date"] = pd.to_datetime(df["date"])

series_order = [
    "obs_gdp", "obs_hours", "obs_wages", "obs_gdpdeflator",
    "obs_corepce", "obs_nominalrate", "obs_consumption", "obs_investment",
    "obs_BBBspread", "obs_longinflation", "obs_longrate", "obs_tfp",
    "obs_gdi", "obs_AAAspread", "obs_nominalrate1", "obs_nominalrate2",
    "obs_nominalrate3", "obs_nominalrate4", "obs_nominalrate5", "obs_nominalrate6",
]

pretty_names = {
    "obs_gdp": "gdp",
    "obs_hours": "hours",
    "obs_wages": "wages",
    "obs_gdpdeflator": "gdpdeflator",
    "obs_corepce": "corepce",
    "obs_nominalrate": "nominalrate",
    "obs_consumption": "consumption",
    "obs_investment": "investment",
    "obs_BBBspread": "BBBspread",
    "obs_longinflation": "longinflation",
    "obs_longrate": "longrate",
    "obs_tfp": "tfp",
    "obs_gdi": "gdi",
    "obs_AAAspread": "AAAspread",
    "obs_nominalrate1": "nominalrate1",
    "obs_nominalrate2": "nominalrate2",
    "obs_nominalrate3": "nominalrate3",
    "obs_nominalrate4": "nominalrate4",
    "obs_nominalrate5": "nominalrate5",
    "obs_nominalrate6": "nominalrate6",
}

notes = (
    "Notes: obs_gdp=real output growth; obs_hours=hours worked; obs_wages=real wage growth; "
    "obs_gdpdeflator=GDP deflator inflation; obs_corepce=core PCE inflation; obs_nominalrate=policy rate; "
    "obs_consumption=consumption growth; obs_investment=investment growth; obs_BBBspread=BBB (Baa) credit spread; "
    "obs_longinflation=long-run inflation expectations; obs_longrate=long-run bond yield; obs_tfp=TFP growth; "
    "obs_gdi=gross domestic income growth; obs_AAAspread=AAA credit spread; nominalrate1-6=expected policy rates "
    "1 to 6 quarters ahead. All series are quarterly and expressed in percent."
)

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.8,
    "grid.color": "#D0D0D0",
    "axes.edgecolor": "black",
})

fig, axes = plt.subplots(5, 4, figsize=(16, 20), facecolor="white")
axes = axes.flatten()

for i, col in enumerate(series_order):
    ax = axes[i]
    data = df.loc[df[col].notna(), ["date", col]]

    ax.set_facecolor("white")
    ax.plot(data["date"], data[col], linewidth=1.0, color="#1f77b4")
    ax.set_title(pretty_names[col], pad=4)

    span_years = (data["date"].max() - data["date"].min()).days / 365.25
    if span_years >= 35:
        ax.xaxis.set_major_locator(mdates.YearLocator(10))
    elif span_years >= 20:
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
    else:
        ax.xaxis.set_major_locator(mdates.YearLocator(2))

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.tick_params(axis="x", rotation=0)

fig.suptitle("Euro Area Observables (all variables)", fontsize=18, y=0.985)
fig.text(0.01, 0.035, notes, ha="left", va="bottom", fontsize=8, wrap=True)
fig.text(0.5, 0.02, "Appendix B.1 for measurement equations.", ha="center", va="bottom", fontsize=8)

fig.tight_layout(rect=[0.02, 0.07, 0.99, 0.96])
fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
print(f"Saved chart to {output_path.resolve()}")
