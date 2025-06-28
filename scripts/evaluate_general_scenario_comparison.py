import os
import re
import pypsa
import logging
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scripts._evaluation_helpers import load_networks_from_path_list, load_csvs_from_path_list, compare_value, plot_line_comparison
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)




if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_general_scenario_comparison",
        )

    configure_logging(snakemake)

    logger.info("loading data.")
    #nn = load_networks_from_path_list(snakemake.input.networks)
    cc_capex = load_csvs_from_path_list(snakemake.input.statistics_capex_buscarrier_csv)
    cc_opex = load_csvs_from_path_list(snakemake.input.statistics_opex_buscarrier_csv)
    cc_optimalcapacity = load_csvs_from_path_list(snakemake.input.statistics_optimalcapacity_buscarrier_csv)

    planning_horizons = snakemake.params.planning_horizons

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX
    # ─────────────────────────────────────────────────────────────────────────────

    cc_capex_opex = {}
    cc_capex_opex_de = {}

    for run in cc_capex:
        df_cap = cc_capex[run]
        df_ope = cc_opex[run]
        year_cols = [str(y) for y in planning_horizons]

        # Gesamt
        df_all = df_cap[year_cols].add(df_ope[year_cols], fill_value=0)
        cc_capex_opex[run] = df_all

        # Nur DE
        df_cap_de = df_cap[df_cap["country"] == "DE"]
        df_ope_de = df_ope[df_ope["country"] == "DE"]

        df_de = df_cap_de[year_cols].add(df_ope_de[year_cols], fill_value=0)
        cc_capex_opex_de[run] = df_de


    def get_capex_opex_series(df_run):
        s = df_run.sum(axis=0)
        s.index = s.index.astype(int)
        return s.sort_index() * 1e-9

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    plot_line_comparison(
        cc=cc_capex_opex,
        title="CAPEX + OPEX Total (Bn. Euro)",
        expr=get_capex_opex_series,
        #output=snakemake.output.total_capexopex_graph,
        ax=ax1,
    )

    plot_line_comparison(
        cc=cc_capex_opex_de,
        title="CAPEX + OPEX DE (Bn. Euro)",
        expr=get_capex_opex_series,
        #output=snakemake.output.DE_capexopex_graph,
        ax=ax2,
    )

    fig.tight_layout()
    fig.savefig(snakemake.output.total_and_DE_capexopex_graph)
    logger.info(f"Created plot Total and DE capexopex graph and saved to {snakemake.output.total_and_DE_capexopex_graph}")

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX Stackplot Raster
    # ─────────────────────────────────────────────────────────────────────────────

    # --- carrier_group_map wie in deinem zweiten Skript ---
    carrier_group_map = {
        "fossil plants": [
            "gas", "coal", "oil", "lignite", "CCGT", "Open-Cycle Gas", "Combined-Cycle Gas",
            "urban central gas CHP", "urban central oil CHP", "urban central coal CHP",
            "urban central lignite CHP", "CHP", "gas primary", "oil primary", "coal for industry",
            "gas for industry", "gas for industry CC", "oil refining"
        ],
        "sustainable biomass": [
            "biogas", "solid biomass", "biomass to liquid",
            "biogas to gas", "biogas to gas CC", "biomass to liquid CC", "solid biomass for industry",
            "solid biomass for industry CC"
        ],
        "unsustainable biomass": [
            "unsustainable", "unsustainable biogas", "unsustainable solid biomass", "unsustainable bioliquids"
        ],
        "onshore wind": ["onshore wind"],
        "offshore wind": ["offshore wind", "offshore wind (AC)", "offshore wind (DC)"],
        "solar": ["solar"],
        "solar rooftop": ["solar rooftop"],
        "hydro power": ["run of river", "reservoir", "dam", "hydro", "reservoir & dam", "Pumped Hydro Storage"],

        "boiler": [
            "boiler",
            "residential rural biomass boiler", "residential rural gas boiler", "residential rural oil boiler",
            "residential urban decentral biomass boiler", "residential urban decentral gas boiler",
            "residential urban decentral oil boiler",
            "services rural biomass boiler", "services rural gas boiler", "services rural oil boiler",
            "services urban decentral biomass boiler", "services urban decentral gas boiler",
            "services urban decentral oil boiler",
            "urban central gas boiler"
        ],
        "heat pump": [
            "air heat pump", "ground heat pump",
            "residential rural air heat pump", "residential rural ground heat pump",
            "residential urban decentral air heat pump",
            "services rural air heat pump", "services rural ground heat pump",
            "services urban decentral air heat pump",
            "urban central air heat pump"
        ],
        # "other heating": [
        #     "heat", "vent", "resistive heater", "water tanks",
        #     "residential rural heat vent", "residential rural resistive heater", "residential rural water tanks",
        #     "residential rural water tanks charger", "residential rural water tanks discharger",
        #     "residential urban decentral heat vent", "residential urban decentral resistive heater",
        #     "residential urban decentral water tanks", "residential urban decentral water tanks charger",
        #     "residential urban decentral water tanks discharger",
        #     "services rural heat vent", "services rural resistive heater", "services rural water tanks",
        #     "services rural water tanks charger", "services rural water tanks discharger",
        #     "services urban decentral heat vent", "services urban decentral resistive heater",
        #     "services urban decentral water tanks", "services urban decentral water tanks charger",
        #     "services urban decentral water tanks discharger",
        #     "urban central heat vent", "urban central resistive heater", "urban central water tanks",
        #     "urban central water tanks charger", "urban central water tanks discharger"
        # ],

        "hydrogen": [
            "H2", "hydrogen", "electrolysis", "fuel cell", "sabatier",
            "H2 Electrolysis", "H2 Fuel Cell", "H2 OCGT", "SMR", "SMR CC",
            "H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted", "H2 Store"
        ],
        "nuclear": ["nuclear", "uranium"],
        # "CO2 removal": [
        #     "DAC", "co2 sequestered", "co2 stored", "process emissions",
        #     "co2", "process emissions CC"
        # ],
        "electricity grid": [
            "electricity distribution grid", "AC", "DC"
        ],

        "storage": [
            "battery", "thermal storage", "home battery", "hydrogen storage",
            "Battery Storage", "battery charger", "battery discharger",
            "home battery charger", "home battery discharger"
        ],
        # "mobility": [
        #     "BEV charger", "land transport oil", "kerosene for aviation", "shipping methanol", "shipping oil",
        #     "agriculture machinery oil"
        # ],
        # "methanol": [
        #     "methanol", "methanolisation", "industry methanol"
        # ],
        # "HVC": ["HVC to air"],
        # "waste": ["waste CHP", "waste CHP CC"],
        # "renewable fuels": ["renewable gas", "renewable oil", "Fischer-Tropsch"],
        # "naphtha": ["naphtha for industry"]
    }


    def assign_group(carrier):
        c = carrier.lower().strip()
        for group, keywords in carrier_group_map.items():
            for k in keywords:
                if c == k.lower().strip():
                    return group
        return "other"


    def extract_constraint_from_run(run):
        match = re.search(r"constr(\d+)", run)
        if match:
            return int(match.group(1))
        else:
            return 100  # default fallback for runs without constraint


    # --- Prepare input data (CAPEX + OPEX in 2050) ---
    def get_combined_2050(cc_capex, cc_opex):
        combined = {}

        for run in cc_capex:
            df_cap = cc_capex[run][["carrier", "country", "2050"]].rename(
                columns={"2050": "capex"}
            )
            df_ope = cc_opex[run][["carrier", "country", "2050"]].rename(
                columns={"2050": "opex"}
            )

            # Outer merge, damit auch Carrier ohne capex (nur opex) und umgekehrt drinbleiben
            df = pd.merge(
                df_cap,
                df_ope,
                on=["carrier", "country"],
                how="outer"
            )

            # Fehlende Werte auf 0 setzen
            df["capex"] = df["capex"].fillna(0)
            df["opex"] = df["opex"].fillna(0)

            # Gesamtkosten
            df["total"] = df["capex"] + df["opex"]

            # Gruppierung und Metadaten
            df["carrier_grouped"] = df["carrier"].apply(assign_group)
            df["constraint"] = extract_constraint_from_run(run)
            df["run"] = run

            combined[run] = df

        return combined


    combined_by_run = get_combined_2050(cc_capex, cc_opex)

    # --- Prepare DataFrames for plotting ---
    df_all = pd.concat(combined_by_run.values(), ignore_index=True)
    df_all["constraint"] = df_all["constraint"].astype(int)


    # --- Helper for filtering ---
    def filter_and_group(df, run_prefix, countries):
        subset = df[df["run"].str.startswith(run_prefix)]
        if countries is not None:
            subset = subset[subset["country"].isin(countries)]
        grouped = subset.groupby(["constraint", "carrier_grouped"])["total"].sum().unstack().fillna(0) / 1e9
        return grouped


    reference_all = filter_and_group(df_all, "ref", None)
    reference_de = filter_and_group(df_all, "ref", ["DE"])
    clever_all = filter_and_group(df_all, "cle", None)
    clever_de = filter_and_group(df_all, "cle", ["DE"])

    # --- Sort carrier groups for consistent stack order ---
    all_carriers = pd.concat([reference_all, reference_de, clever_all, clever_de])
    carrier_order = all_carriers.sum().sort_values(ascending=False).index.tolist()

    # --- Plotting: 2×2 Grid ---
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey='row')

    # Geänderte Rasterlogik: Oben = All, Unten = DE; Links = ref, Rechts = clever
    plot_data = [
        (reference_all, "Reference – All Countries", axs[0, 0]),
        (clever_all, "Clever – All Countries", axs[0, 1]),
        (reference_de, "Reference – Germany", axs[1, 0]),
        (clever_de, "Clever – Germany", axs[1, 1]),
    ]

    color_map = {}

    for data, title, ax in plot_data:
        data = data.reindex(columns=carrier_order, fill_value=0)
        x = data.index
        y = data.values.T
        polys = ax.stackplot(x, y, labels=data.columns)

        if not color_map:
            color_map = {name: poly.get_facecolor()[0] for name, poly in zip(data.columns, polys)}
        else:
            for poly, name in zip(polys, data.columns):
                poly.set_color(color_map[name])

        ax.set_title(title)
        ax.set_xlabel("Constraint (%)")
        ax.set_ylabel("CAPEX + OPEX in 2050 (Bn. EUR)")
        ax.set_xticks([1, 25, 50, 100])
        ax.set_xlim(105, 0)  # ← ACHTUNG: Achse umgedreht

    # Gemeinsame Legende (reversed für passende Stapelreihenfolge)
    handles = [mpatches.Patch(label=k, color=color_map[k]) for k in reversed(carrier_order)]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.85, 0.5), title="Technology Group")

    fig.tight_layout(rect=[0, 0, 0.85, 1])
    #plt.show()
    fig.savefig(snakemake.output.total_and_DE_capexopex_stackplot_2050)
    logger.info(f"Saved stackplot of 2050 technologies to {snakemake.output.total_and_DE_capexopex_stackplot_2050}")


    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX Relative Lineplot
    # ─────────────────────────────────────────────────────────────────────────────

    # 1) Compute relative series
    def compute_relative(df, ref100):
        return df.div(ref100, axis=1) * 100


    ref_all_100 = reference_all.loc[100]
    clever_all_100 = clever_all.loc[100]
    ref_de_100 = reference_de.loc[100]
    clever_de_100 = clever_de.loc[100]

    rel_ref_all = compute_relative(reference_all, ref_all_100)
    rel_clever_all = compute_relative(clever_all, clever_all_100)
    rel_ref_de = compute_relative(reference_de, ref_de_100)
    rel_clever_de = compute_relative(clever_de, clever_de_100)

    # 2) Define paired groups for shared colors & linestyles
    pair_groups = [
        ("onshore wind", "offshore wind"),
        ("electricity grid", "storage"),
        ("solar", "solar rooftop"),
        ("heat pump", "boiler"),
        ("sustainable biomass", "unsustainable biomass"),
        ("nuclear", "fossil plants"),
    ]

    paired = {c for pair in pair_groups for c in pair}
    others = [c for c in carrier_order if c not in paired]

    cmap = plt.get_cmap("tab20")
    n_colors = len(pair_groups) + len(others)
    colors = [cmap(i % cmap.N) for i in range(n_colors)]

    color_map = {}
    for idx, (c1, c2) in enumerate(pair_groups):
        color_map[c1] = colors[idx]
        color_map[c2] = colors[idx]
    for j, c in enumerate(others, start=len(pair_groups)):
        color_map[c] = colors[j]

    linestyle_map = {}
    for c1, c2 in pair_groups:
        linestyle_map[c1] = "-"
        linestyle_map[c2] = "--"
    for c in others:
        linestyle_map[c] = "-"

    # 3) Plotting: 2×2 grid
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    panels = [
        (rel_ref_all, "Reference – All Countries", axs[0, 0]),
        (rel_clever_all, "Clever – All Countries", axs[0, 1]),
        (rel_ref_de, "Reference – Germany", axs[1, 0]),
        (rel_clever_de, "Clever – Germany", axs[1, 1]),
    ]

    for df_rel, title, ax in panels:
        df = df_rel.reindex(columns=carrier_order, fill_value=0)
        for carrier in carrier_order:
            ax.plot(
                df.index,
                df[carrier],
                color=color_map[carrier],
                linestyle=linestyle_map[carrier],
                label=carrier
            )
        ax.set_title(title)
        ax.set_xlabel("Constraint (%)")
        ax.set_xticks([1, 25, 50, 100])
        ax.set_xlim(105, 0)
        ax.set_ylim(0, 250)
        ax.grid(True)

    # Common y-label
    axs[0, 0].set_ylabel("Relative to 100% baseline (%)")

    # Build legend handles with pairs grouped together
    legend_order = []
    for c1, c2 in pair_groups:
        legend_order += [c1, c2]
    legend_order += others

    handles = [
        mpl.lines.Line2D([0], [0],
                         color=color_map[carrier],
                         linestyle=linestyle_map[carrier],
                         lw=2)
        for carrier in legend_order
    ]

    fig.legend(
        handles,
        legend_order,
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        title="Technology Group"
    )

    fig.tight_layout(rect=[0, 0, 0.85, 1])
    #plt.show()
    fig.savefig(snakemake.output.total_and_DE_capexopex_relative_change_2050)

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX Absolute Change Lineplot
    # ─────────────────────────────────────────────────────────────────────────────

    import matplotlib.pyplot as plt
    import matplotlib as mpl


    # 1) Compute absolute change series (Bn EUR) vs. 100% constraint baseline
    def compute_absolute_change(df, base100):
        """
        Subtract the 100% constraint values (base100) from each series.
        Result is absolute change in Bn EUR.
        """
        return df.sub(base100, axis=1)


    # Baselines at 100% constraint
    base_ref_all = reference_all.loc[100]
    base_clever_all = clever_all.loc[100]
    base_ref_de = reference_de.loc[100]
    base_clever_de = clever_de.loc[100]

    # Absolute change DataFrames
    abs_ref_all = compute_absolute_change(reference_all, base_ref_all)
    abs_clever_all = compute_absolute_change(clever_all, base_clever_all)
    abs_ref_de = compute_absolute_change(reference_de, base_ref_de)
    abs_clever_de = compute_absolute_change(clever_de, base_clever_de)

    # 2) Paired groups for shared colors & linestyles (reuse from relative plot)
    pair_groups = [
        ("onshore wind", "offshore wind"),
        ("electricity grid", "storage"),
        ("solar", "solar rooftop"),
        ("heat pump", "boiler"),
        ("sustainable biomass", "unsustainable biomass"),
        ("nuclear", "fossil plants"),
    ]

    paired = {c for pair in pair_groups for c in pair}
    others = [c for c in carrier_order if c not in paired]

    cmap = plt.get_cmap("tab20")
    n_colors = len(pair_groups) + len(others)
    colors = [cmap(i % cmap.N) for i in range(n_colors)]

    color_map = {}
    for idx, (c1, c2) in enumerate(pair_groups):
        color_map[c1] = colors[idx]
        color_map[c2] = colors[idx]
    for j, c in enumerate(others, start=len(pair_groups)):
        color_map[c] = colors[j]

    linestyle_map = {}
    for c1, c2 in pair_groups:
        linestyle_map[c1] = "-"
        linestyle_map[c2] = "--"
    for c in others:
        linestyle_map[c] = "-"

    # 3) Plotting: 2×2 grid of absolute-change lineplots
    fig, axs = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    panels = [
        (abs_ref_all, "Reference – All Countries", axs[0, 0]),
        (abs_clever_all, "Clever – All Countries", axs[0, 1]),
        (abs_ref_de, "Reference – Germany", axs[1, 0]),
        (abs_clever_de, "Clever – Germany", axs[1, 1]),
    ]

    for df_abs, title, ax in panels:
        df = df_abs.reindex(columns=carrier_order, fill_value=0)
        for carrier in carrier_order:
            ax.plot(
                df.index,
                df[carrier],
                color=color_map[carrier],
                linestyle=linestyle_map[carrier],
                label=carrier
            )
        ax.set_title(title)
        ax.set_xlabel("Constraint (%)")
        ax.set_xticks([1, 25, 50, 100])
        ax.set_xlim(105, 0)
        ax.grid(True)

    # Common y-label for absolute change
    axs[0, 0].set_ylabel("Absolute change from 100% baseline (Bn EUR)")

    # Build and place legend with paired groups together
    legend_order = []
    for c1, c2 in pair_groups:
        legend_order += [c1, c2]
    legend_order += others

    handles = [
        mpl.lines.Line2D(
            [0], [0],
            color=color_map[carrier],
            linestyle=linestyle_map[carrier],
            lw=2
        )
        for carrier in legend_order
    ]

    fig.legend(
        handles,
        legend_order,
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        title="Technology Group"
    )

    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()
    fig.savefig(snakemake.output.total_and_DE_capexopex_absolute_change_2050)


    # ─────────────────────────────────────────────────────────────────────────────
    #  Generation
    # ─────────────────────────────────────────────────────────────────────────────

    # Solar
    def get_solar_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Solar")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Solar capacity (GW)",
        expr=get_solar_capacity_series,
        output=snakemake.output.gen_solar_graph,
    )


    #  Onshore Wind
    def get_onwind_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Onshore Wind")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Onshore Wind capacity (GW)",
        expr=get_onwind_capacity_series,
        output=snakemake.output.gen_onwind_graph,
    )


    #  Offshore Wind (AC)
    def get_offwind_ac_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Offshore Wind (AC)")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Offshore Wind AC capacity (GW)",
        expr=get_offwind_ac_capacity_series,
        output=snakemake.output.gen_offwind_ac_graph,
    )


    #  Offshore Wind (DC)
    def get_offwind_dc_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Offshore Wind (DC)")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Offshore Wind DC capacity (GW)",
        expr=get_offwind_dc_capacity_series,
        output=snakemake.output.gen_offwind_dc_graph,
    )

    logger.info("All general scenario comparisons done.")