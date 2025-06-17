import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logging
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake
        snakemake = mock_snakemake("evaluate_general_run", run="ref-constr50")

    configure_logging(snakemake)

    capex_df = pd.read_csv(snakemake.input.capex_csv, index_col=0)
    opex_df = pd.read_csv(snakemake.input.opex_csv, index_col=0)
    planning_horizons = snakemake.params.planning_horizons

    # Merge CAPEX and OPEX completely
    cap = capex_df.copy()
    opex = opex_df.copy()
    year_cols = [str(y) for y in planning_horizons]

    # Outer join to preserve all rows
    combined_df = cap[["component", "bus", "carrier", "country"]].copy()
    combined_df = combined_df.join(opex[["component", "bus", "carrier", "country"]], how="outer", lsuffix="_cap",
                                   rsuffix="_opex")

    # Prioritise non-null metadata
    for col in ["component", "bus", "carrier", "country"]:
        combined_df[col] = combined_df[f"{col}_cap"].combine_first(combined_df[f"{col}_opex"])
        combined_df.drop(columns=[f"{col}_cap", f"{col}_opex"], inplace=True)

    # Merge year-wise values
    cap_vals = cap[year_cols]
    opex_vals = opex[year_cols]
    combined_vals = cap_vals.add(opex_vals, fill_value=0)

    # Final assembly
    combined_df = combined_df.join(combined_vals)

    # --- Define carrier grouping rules ---
    carrier_group_map = {
        "fossil plants": [
            "gas", "coal", "oil", "lignite", "CCGT", "Open-Cycle Gas", "Combined-Cycle Gas",
            "urban central gas CHP", "urban central oil CHP", "urban central coal CHP",
            "urban central lignite CHP", "CHP"
        ],
        "sustainable biomass": [
            "biogas", "solid biomass", "biomass to liquid"
        ],
        "unsustainable biomass": [
            "unsustainable", "unsustainable biogas", "unsustainable solid biomass", "unsustainable bioliquids"
        ],
        "onshore wind": ["onshore wind"],
        "offshore wind": ["offshore wind"],
        "solar": ["solar"],
        "solar rooftop": ["solar rooftop"],
        "hydro power": ["run of river", "reservoir", "pumped hydro storage", "dam", "hydro"],
        "heat systems": ["heat", "boiler", "resistive heater", "air heat pump", "water tanks", "vent"],
        "storage": ["battery", "H2 store", "thermal storage", "home battery", "water tanks", "hydrogen storage"],
        "hydrogen": [
            "H2", "hydrogen", "electrolysis", "fuel cell", "SMR", "sabatier"#, "methanol", "methanolisation"
        ],
        "nuclear": ["nuclear", "uranium"],
        "CO2 removal": ["DAC", "co2 sequestered", "co2 stored", "process emissions"],
        "electricity grid": ["electricity distribution grid"]#, "AC", "DC"],
    }


    def assign_group(carrier):
        for group, keywords in carrier_group_map.items():
            if any(k.lower() in carrier.lower() for k in keywords):
                return group
        return "other"

    combined_df["carrier_grouped"] = combined_df["carrier"].apply(assign_group)

    # --- Log carriers grouped as 'other' ---
    other_carriers = combined_df[combined_df["carrier_grouped"] == "other"]["carrier"].unique()
    logger.info(f"Carriers grouped under 'other': {sorted(other_carriers.tolist())}")

    # --- Group and prepare data ---
    def prepare_plot_data(df, countries=None):
        if countries:
            df = df[df["country"].isin(countries)]
        grouped = df.groupby("carrier_grouped")[[str(y) for y in planning_horizons]].sum()
        grouped = grouped / 1e9  # Convert to Bn EUR
        return grouped

    all_data = prepare_plot_data(combined_df)
    de_data = prepare_plot_data(combined_df, countries=["DE"])

    # --- Prepare sorted carrier list (same order left and right, bottom-up) ---
    sorted_carriers = all_data.sum(axis=1).sort_values(ascending=False).index.tolist()

    # --- Plotting with custom legend, x-ticks, y-axes ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=False)

    # Store plotted colors for shared legend
    color_map = {}

    for ax, data, title in zip(axes, [all_data, de_data], ["All countries", "Germany only"]):
        # Reindex to ensure same stack order
        data = data.reindex(sorted_carriers).fillna(0)
        values = data.values
        polys = ax.stackplot(planning_horizons, values, labels=data.index)

        # Save color mapping from left plot
        if title == "All countries":
            color_map = dict(zip(data.index, [p.get_facecolor()[0] for p in polys]))

        # For right plot, set colors to match left
        else:
            for p, key in zip(polys, data.index):
                p.set_color(color_map[key])

        ax.set_title(title)
        ax.set_xlabel("Year")
        ax.set_xticks(planning_horizons)
        ax.set_xlim(min(planning_horizons), max(planning_horizons))
        ax.set_ylabel("CAPEX + OPEX (Bn. Euro)" if title == "All countries" else None)

    # --- Shared legend (outside plot) with matching order ---
    legend_handles = [mpatches.Patch(label=k, color=color_map[k]) for k in reversed(sorted_carriers)]
    fig.legend(handles=legend_handles,
               loc="center left",
               bbox_to_anchor=(0.85, 0.5),
               title="Technology Group")

    # Add run label as top-right title
    fig.text(0.98, 0.98, f"Run: {snakemake.wildcards.run}",
             ha="right", va="top", fontsize=12)

    plt.tight_layout(rect=[0, 0, 0.8, 1])
    plt.show()
    plt.savefig(snakemake.output.total_and_DE_capex_opex_graph)
    plt.close()