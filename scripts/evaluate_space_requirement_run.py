import os
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import logging
from scripts._evaluation_helpers import load_networks_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def space_req_DLU_DE(df_DLU_country, run_label=""):
    total_area_de = 357.022  # 1000 km²
    planning_horizons = df_DLU_country.columns

    # Daten vorbereiten: summe über alle Länder und nur DE
    df_all = df_DLU_country.groupby("carrier").sum()
    df_de = df_DLU_country.loc["DE"]

    # Nur Carrier mit > 0 Flächenverbrauch
    df_all = df_all.loc[(df_all != 0).any(axis=1)]
    df_de = df_de.loc[(df_de != 0).any(axis=1)]

    # Gleiche Carrier-Reihenfolge
    sorted_carriers = sorted(set(df_all.index) | set(df_de.index))
    df_all = df_all.reindex(sorted_carriers).fillna(0)
    df_de = df_de.reindex(sorted_carriers).fillna(0)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    color_map = {}

    # --- Linker Plot: alle Länder ---
    vals_all = df_all.values
    polys_all = ax1.stackplot(planning_horizons, vals_all, labels=df_all.index)
    color_map = dict(zip(df_all.index, [p.get_facecolor()[0] for p in polys_all]))

    ax1.set_title("DLU (all countries)")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("1000 km²")
    ax1.set_xticks(planning_horizons)
    ax1.set_xlim(min(planning_horizons), max(planning_horizons))

    # --- Rechter Plot: nur DE ---
    vals_de = df_de.values
    polys_de = ax2.stackplot(planning_horizons, vals_de, labels=df_de.index)
    for p, k in zip(polys_de, df_de.index):
        p.set_color(color_map[k])

    if is_constrained:
        max_limits_DE = max_limits_df.loc["DE"]
        max_limit_series = max_limits_DE.copy()
        max_limit_series = max_limit_series * 1e-9
        max_limit_series.loc[2020] = df_de.sum().loc[2020]
        max_limit_series.index = max_limit_series.index.astype(int)
        max_limit_series = max_limit_series.sort_index()
        ax2.plot(planning_horizons, max_limit_series, linestyle="--", color="black", label="Max limit")

        # --- Max limit total (alle Länder) ---
        max_limit_total = max_limits_df.sum() * 1e-9
        max_limit_total.loc[2020] = df_all.sum().loc[2020]
        max_limit_total.index = max_limit_total.index.astype(int)
        max_limit_total = max_limit_total.sort_index()
        ax1.plot(planning_horizons, max_limit_total, linestyle="--", color="black", label="Max limit")

    ax2.set_title("DLU DE")
    ax2.set_xlabel("Year")
    ax2.set_xticks(planning_horizons)
    ax2.set_xlim(min(planning_horizons), max(planning_horizons))

    # Sekundärachse: Anteil DE-Fläche
    def abs_to_rel(x): return x / total_area_de
    def rel_to_abs(x): return x * total_area_de
    secax = ax2.secondary_yaxis("right", functions=(abs_to_rel, rel_to_abs))
    secax.set_ylabel("% of DE area")
    secax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))

    # Legende
    legend_handles = [mpatches.Patch(label=k, color=color_map[k]) for k in reversed(sorted_carriers)]
    if is_constrained:
        legend_handles.append(mpatches.Patch(label="Max limit", color="black"))
    fig.legend(handles=legend_handles,
               loc="center left",
               bbox_to_anchor=(0.75, 0.5),
               title="Carrier")

    # Run-Label (optional)
    if run_label:
        fig.text(0.98, 0.98, f"Run: {run_label}", ha="right", va="top", fontsize=12)

    plt.tight_layout(rect=[0, 0, 0.75, 1])
    plt.show()
    fig.savefig(snakemake.output.DLU_DE, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_space_requirement_run",
            run = "ref-constr25",
        )

    configure_logging(snakemake)
    DLU_config=snakemake.params.DLU_config
    is_constrained = DLU_config["constraint"]["enable"]

    df_DLU = pd.read_csv(snakemake.input.space_requirements_DLU_csv)
    df_DLU_country = (
        df_DLU.drop(columns=["bus"])
        .groupby(["country", "carrier"], as_index=True)
        .sum(numeric_only=True)
    )
    df_DLU_country = df_DLU_country / 1e9  # m² → 1000 km²
    df_DLU_country.columns = df_DLU_country.columns.astype(int)

    if "EU" in df_DLU_country.index.get_level_values(0):
        df_DLU_country = df_DLU_country.drop(index="EU", level="country")


    if is_constrained:
        max_limits = DLU_config["constraint"].get("max_limit")
        max_limits_df = pd.DataFrame.from_dict(max_limits, orient='index').apply(pd.to_numeric)
        logger.info("constraints loaded.")
    else:
        logger.info("No space requirement constraints defined.")

    space_req_DLU_DE(df_DLU_country, run_label=snakemake.wildcards.run)