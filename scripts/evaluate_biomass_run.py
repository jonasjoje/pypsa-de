import os
import pandas as pd
import matplotlib.pyplot as plt
import logging
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def create_total_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path):

    # Jahre als Strings für DataFrame-Spalten
    year_columns = [str(y) for y in planning_horizons]

    # Relevante Carrier
    biomass_carriers = [
        "biogas", "unsustainable biogas",
        "solid biomass", "unsustainable solid biomass",
        "unsustainable bioliquids"
    ]

    # Filter auf relevante Carrier
    supply_biomass = supply_df[supply_df["carrier"].isin(biomass_carriers)]
    e_sum_max_biomass = e_sum_max_df[e_sum_max_df["carrier"].isin(biomass_carriers)]
    e_sum_min_biomass = e_sum_min_df[e_sum_min_df["carrier"].isin(biomass_carriers)]

    # Gesamterzeugung pro Jahr berechnen
    supply_total = supply_biomass[year_columns].sum()
    e_max = e_sum_max_biomass[year_columns].sum()
    e_min = e_sum_min_biomass[year_columns].sum()

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))

    # Potenzialbandbreite als Fläche
    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")

    # Tatsächliche Erzeugung als Linie
    ax.plot(year_columns, supply_total, label="Actual biomass generation", color="green", linewidth=2)

    # Achsenbeschriftung und Layout
    ax.set_title("Biomass Generation and Potential")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()
    plt.tight_layout()
    #plt.show()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")


def create_total_biocrops_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path):

    # Jahre als Strings
    year_columns = [str(y) for y in planning_horizons]

    # Nur unsustainable Carrier
    biocrop_carriers = [
        "unsustainable biogas",
        "unsustainable solid biomass",
        "unsustainable bioliquids"
    ]

    # Filter
    supply_biocrops = supply_df[supply_df["carrier"].isin(biocrop_carriers)]
    e_sum_max_biocrops = e_sum_max_df[e_sum_max_df["carrier"].isin(biocrop_carriers)]
    e_sum_min_biocrops = e_sum_min_df[e_sum_min_df["carrier"].isin(biocrop_carriers)]

    # Summe je Jahr
    supply_total = supply_biocrops[year_columns].sum()
    e_max = e_sum_max_biocrops[year_columns].sum()
    e_min = e_sum_min_biocrops[year_columns].sum()

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(year_columns, supply_total, label="Actual biocrop generation", color="darkred", linewidth=2)

    ax.set_title("Generation from Biocrops (Unsustainable Biomass)")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()
    plt.tight_layout()
    plt.show()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")



if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake
        snakemake = mock_snakemake(
            "evaluate_biomass_run",
            run="ref-constr50",
        )

    configure_logging(snakemake)

    # Load inputs
    e_sum_min_df = pd.read_csv(snakemake.input.e_sum_min_csv, index_col=0)
    e_sum_max_df = pd.read_csv(snakemake.input.e_sum_max_csv, index_col=0)
    supply_stats_df = pd.read_csv(snakemake.input.statistics_supply_csv, index_col=0)
    planning_horizons = snakemake.params.planning_horizons

    # Call evaluation functions
    create_total_biomass_graph(
        e_sum_min_df,
        e_sum_max_df,
        supply_stats_df,
        snakemake.output.total_biomass
    )

    create_total_biocrops_graph(
        e_sum_min_df,
        e_sum_max_df,
        supply_stats_df,
        snakemake.output.total_biocrops
    )

    # generation stackplot

    # unsustainable solid biomass
