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


    # ── Globale Einstellungen ──
    plt.rcParams['font.size'] = 7
    plt.rcParams['savefig.dpi'] = 300

    # ── Figure mit 5.46" Breite und altem 12:5-Verhältnis ──
    width = 5.46
    height = width * (5 / 8)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, height))

    # ── Linker Plot: alle Länder ──
    plot_line_comparison(
        cc=cc_capex_opex,
        title="All countries",
        expr=get_capex_opex_series,
        ax=ax1,
    )
    ax1.set_ylabel("CAPEX + OPEX in Bn. Euro")
    ax1.set_xticks([2020, 2030, 2040, 2050])
    ax1.set_xlabel("Year")
    ax1.grid(True, linestyle='--', linewidth=0.5)  # gestricheltes Grid

    # ── Rechter Plot: nur Deutschland ──
    plot_line_comparison(
        cc=cc_capex_opex_de,
        title="Germany",
        expr=get_capex_opex_series,
        ax=ax2,
    )
    ax2.set_ylabel("CAPEX + OPEX in Bn. Euro")
    ax2.set_xticks([2020, 2030, 2040, 2050])
    ax2.set_xlabel("Year")
    ax2.grid(True, linestyle='--', linewidth=0.5)

    # ── Linienbreiten halbieren ──
    for ax in (ax1, ax2):
        for line in ax.get_lines():
            line.set_linewidth(line.get_linewidth() / 2)

    # ── Gemeinsame Legende ──
    # erst Handles/Labels einsammeln, dann einzelne Legenden entfernen
    handles, labels = ax1.get_legend_handles_labels()
    ax1.get_legend().remove()
    ax2.get_legend().remove()

    fig.legend(
        handles, labels,
        loc='lower center',
        ncol=3,  # <— hier 3 Spalten erzwingen
        bbox_to_anchor=(0.5, 0.0),
        columnspacing=5.0,  # optional: Abstand zwischen Spalten
    )

    # ── Layout anpassen und speichern ──
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.3)  # Platz für Legende
    #plt.show()
    fig.savefig(
        snakemake.output.total_and_DE_capexopex_graph,
        dpi=300
    )
    logger.info(
        f"Created plot Total and DE capexopex graph and saved to "
        f"{snakemake.output.total_and_DE_capexopex_graph}"
    )

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX Stackplot Raster
    # ─────────────────────────────────────────────────────────────────────────────

    # ── Globale Einstellungen ──
    plt.rcParams['font.size'] = 7
    plt.rcParams['savefig.dpi'] = 300

    # ── Figure mit 5.46" Breite und 1:1-Verhältnis ──
    width = 5.46
    height = width * (7 / 5)
    fig, axs = plt.subplots(2, 2, figsize=(width, height), sharex=True, sharey=True)

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

    # --- Plotting: 2×2 Grid Stackplots ---
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

        # Stackplot zeichnen
        polys = ax.stackplot(x, y, labels=data.columns)

        # Farben konsistent setzen
        if not color_map:
            color_map = {name: poly.get_facecolor()[0] for name, poly in zip(data.columns, polys)}
        else:
            for poly, name in zip(polys, data.columns):
                poly.set_color(color_map[name])

        # Titel & Achsen
        ax.set_title(title)
        ax.set_xticks([1, 25, 50, 100])
        ax.set_xlim(105, 0)
        ax.grid(True, linestyle='--', linewidth=0.5)

    # Obere Reihe: keine X-Achsen-Beschriftung
    for ax in axs[0, :]:
        ax.tick_params(labelbottom=False)
        ax.set_xlabel("")

    # Gemeinsame X-Beschriftung unten
    for ax in axs[1, :]:
        ax.set_xlabel("Constraint (%)")

    # Gemeinsame Y-Beschriftung links unten
    axs[1, 0].set_ylabel("CAPEX + OPEX in 2050 (Bn. EUR)")

    # Legend unten zentriert mit 4 Spalten wie beim Lineplot
    handles = [mpatches.Patch(label=k, color=color_map[k]) for k in reversed(carrier_order)]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        title="Technology Group (from top to bottom)",
        ncol=1,
        columnspacing=2.0
    )

    # Layout anpassen und speichern
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    fig.subplots_adjust(bottom=0.4)
    fig.savefig(snakemake.output.total_and_DE_capexopex_stackplot_2050)

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX Relative Change Lineplot
    # ─────────────────────────────────────────────────────────────────────────────

    # ── Globale Einstellungen ──
    plt.rcParams['font.size'] = 7
    plt.rcParams['savefig.dpi'] = 300

    # ── Figure mit 5.46" Breite und 1:1-Verhältnis ──
    width = 5.46
    height = width * (5 / 5)
    fig, axs = plt.subplots(2, 2, figsize=(width, height), sharex=True, sharey=True)


    # 1) Compute relative series vs. 100% baseline
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

    # 2) Farben & Linestyles vorbereiten (wie oben beim Absolute-Plot)
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

    # 3) Plotting: 2×2 Grid, relative Werte
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
                lw=1
            )
        ax.set_title(title)
        # X-Ticks & Achsenbegrenzung
        ax.set_xticks([1, 25, 50, 100])
        ax.set_xlim(105, 0)
        ax.set_ylim(0,300)
        # gestricheltes Raster
        ax.grid(True, linestyle='--', linewidth=0.5)

    # Obere Reihe: keine X-Achsenbeschriftung
    for ax in axs[0, :]:
        ax.tick_params(labelbottom=False)
        ax.set_xlabel("")

    # Gemeinsame Y-Beschriftung links unten
    axs[1, 0].set_ylabel("Relative to 100% baseline (%)")

    # 4) Legende wie beim Absolute-Plot
    legend_order = []
    for c1, c2 in pair_groups:
        legend_order += [c1, c2]
    legend_order += others

    handles = [
        mpl.lines.Line2D([0], [0],
                         color=color_map[carrier],
                         linestyle=linestyle_map[carrier],
                         lw=1)
        for carrier in legend_order
    ]

    fig.legend(
        handles,
        legend_order,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        title="Technology Group",
        ncol=4,
        columnspacing=2.0
    )

    # Layout & Speichern
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    fig.subplots_adjust(bottom=0.25)
    fig.savefig(snakemake.output.total_and_DE_capexopex_relative_change_2050)

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX Absolute Change Lineplot
    # ─────────────────────────────────────────────────────────────────────────────


    # ── Globale Einstellungen ──
    plt.rcParams['font.size'] = 7
    plt.rcParams['savefig.dpi'] = 300

    # ── Figure mit 5.46" Breite und altem 12:5-Verhältnis ──
    width = 5.46
    height = width * (5 / 5)
    fig, axs = plt.subplots(2, 2, figsize=(width, height), sharex=True, sharey=True)


    # 1) Compute absolute change series (Bn EUR) vs. 100% constraint baseline
    def compute_absolute_change(df, base100):
        return df.sub(base100, axis=1)


    base_ref_all = reference_all.loc[100]
    base_clever_all = clever_all.loc[100]
    base_ref_de = reference_de.loc[100]
    base_clever_de = clever_de.loc[100]

    abs_ref_all = compute_absolute_change(reference_all, base_ref_all)
    abs_clever_all = compute_absolute_change(clever_all, base_clever_all)
    abs_ref_de = compute_absolute_change(reference_de, base_ref_de)
    abs_clever_de = compute_absolute_change(clever_de, base_clever_de)

    # 2) Farben & Linestyles vorbereiten
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
        # gestricheltes Raster
        ax.grid(True, linestyle='--', linewidth=0.5)

    for ax in axs[0, :]:
        ax.tick_params(labelbottom=False)
        ax.set_xlabel("")

    # halbe Linienbreite in allen Subplots
    for ax in axs.flat:
        for line in ax.get_lines():
            line.set_linewidth(line.get_linewidth() / 2)

    axs[1, 0].set_ylabel(
        "Absolute change from\n100% baseline in Bn. Euro",
        #labelpad=8  # optional: Abstand zum Plot-Rand
    )

    # Build and place legend
    legend_order = []
    for c1, c2 in pair_groups:
        legend_order += [c1, c2]
    legend_order += others

    handles = [
        mpl.lines.Line2D(
            [0], [0],
            color=color_map[carrier],
            linestyle=linestyle_map[carrier],
            lw=1  # passt zur halbierten Linienbreite
        )
        for carrier in legend_order
    ]

    fig.legend(
        handles,
        legend_order,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        title="Technology Group",
        ncol=4,
        columnspacing=2.0
    )

    fig.tight_layout(rect=[0, 0, 0.95, 1])
    fig.subplots_adjust(bottom=0.25)
    #plt.show()
    fig.savefig(snakemake.output.total_and_DE_capexopex_absolute_change_2050, dpi=300)


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