import os
import pandas as pd
import matplotlib.pyplot as plt
import logging
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)


def calculate_e_min_unsus_solid_DE(e_sum_max_df, biomass_config, planning_horizons):
    year_columns = [str(y) for y in planning_horizons]
    max_DE = e_sum_max_df[
        (e_sum_max_df["country"] == "DE") &
        (e_sum_max_df["carrier"] == "unsustainable solid biomass")
    ][year_columns].sum()
    share = pd.Series({str(y): biomass_config["share_unsustainable_min"].get(y, 0) for y in planning_horizons})
    e_min_usb = max_DE.astype(float) * share.astype(float)
    return e_min_usb


def create_total_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    biomass_carriers = ["biogas", "unsustainable biogas", "solid biomass", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply_biomass = supply_df[supply_df["carrier"].isin(biomass_carriers)]
    e_sum_max_biomass = e_sum_max_df[e_sum_max_df["carrier"].isin(biomass_carriers)]
    e_sum_min_biomass = e_sum_min_df[e_sum_min_df["carrier"].isin(biomass_carriers)]
    supply_total = pd.to_numeric(supply_biomass[year_columns].sum(), errors="coerce")
    e_max = pd.to_numeric(e_sum_max_biomass[year_columns].sum(), errors="coerce")
    e_min = pd.to_numeric(e_sum_min_biomass[year_columns].sum(), errors="coerce")
    #e_min += pd.to_numeric(globalconstraint_unsus_solid_biomass_min[year_columns])
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(year_columns, supply_total, label="Actual biomass generation", color="green", linewidth=2)
    ax.set_title("Biomass Generation and Potential")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()
    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()


def create_DE_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    biomass_carriers = ["biogas", "unsustainable biogas", "solid biomass", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply_biomass = supply_df[supply_df["carrier"].isin(biomass_carriers)]
    e_sum_max_biomass = e_sum_max_df[e_sum_max_df["carrier"].isin(biomass_carriers)]
    e_sum_min_biomass = e_sum_min_df[e_sum_min_df["carrier"].isin(biomass_carriers)]
    supply_total = pd.to_numeric(supply_biomass[year_columns].sum(), errors="coerce")
    e_max = pd.to_numeric(e_sum_max_biomass[year_columns].sum(), errors="coerce")
    e_min = pd.to_numeric(e_sum_min_biomass[year_columns].sum(), errors="coerce")
    #e_min += e_min_usb_DE
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(year_columns, supply_total, label="Actual biomass generation", color="green", linewidth=2)
    ax.set_title("Biomass Generation and Potential (DE)")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()
    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()


def create_total_biocrops_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    biocrop_carriers = ["unsustainable biogas", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply_biocrops = supply_df[supply_df["carrier"].isin(biocrop_carriers)]
    e_sum_max_biocrops = e_sum_max_df[e_sum_max_df["carrier"].isin(biocrop_carriers)]
    e_sum_min_biocrops = e_sum_min_df[e_sum_min_df["carrier"].isin(biocrop_carriers)]
    supply_total = pd.to_numeric(supply_biocrops[year_columns].sum(), errors="coerce")
    e_max = pd.to_numeric(e_sum_max_biocrops[year_columns].sum(), errors="coerce")
    e_min = pd.to_numeric(e_sum_min_biocrops[year_columns].sum(), errors="coerce")
    #e_min += pd.to_numeric(globalconstraint_unsus_solid_biomass_min[year_columns])
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(year_columns, supply_total, label="Actual biocrop generation", color="darkred", linewidth=2)
    ax.set_title("Generation from Biocrops (Unsustainable Biomass)")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()
    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()


def create_DE_biocrops_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    biocrop_carriers = ["unsustainable biogas", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply_biocrops = supply_df[supply_df["carrier"].isin(biocrop_carriers)]
    e_sum_max_biocrops = e_sum_max_df[e_sum_max_df["carrier"].isin(biocrop_carriers)]
    e_sum_min_biocrops = e_sum_min_df[e_sum_min_df["carrier"].isin(biocrop_carriers)]
    supply_total = pd.to_numeric(supply_biocrops[year_columns].sum(), errors="coerce")
    e_max = pd.to_numeric(e_sum_max_biocrops[year_columns].sum(), errors="coerce")
    e_min = pd.to_numeric(e_sum_min_biocrops[year_columns].sum(), errors="coerce")
    #e_min += e_min_usb_DE
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(year_columns, supply_total, label="Actual biocrop generation", color="darkred", linewidth=2)
    ax.set_title("Biocrop Generation (Unsustainable Biomass, DE)")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()
    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()


def create_total_biomass_stack_graph(supply_df, planning_horizons, output_path=None, ax=None):
    biomass_carriers = ["unsustainable solid biomass", "unsustainable biogas", "unsustainable bioliquids", "solid biomass", "biogas"]
    year_columns = [str(y) for y in planning_horizons]
    supply_biomass = supply_df[supply_df["carrier"].isin(biomass_carriers)].copy()
    supply_biomass["carrier"] = pd.Categorical(supply_biomass["carrier"], categories=biomass_carriers, ordered=True)
    supply_stacked = supply_biomass.groupby("carrier")[year_columns].sum().loc[biomass_carriers].T
    supply_stacked = supply_stacked.apply(pd.to_numeric, errors="coerce")
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    supply_stacked.plot.area(ax=ax)
    ax.set_title("Biomass Generation by Carrier")
    ax.set_ylabel("Generation [MWh]")
    ax.set_xlabel("Year")
    ax.legend(title="Carrier", loc="upper left", bbox_to_anchor=(1.0, 1.0))
    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()


def create_DE_biomass_stack_graph(supply_df, planning_horizons, output_path=None, ax=None):
    biomass_carriers = ["unsustainable solid biomass", "unsustainable biogas", "unsustainable bioliquids", "solid biomass", "biogas"]
    year_columns = [str(y) for y in planning_horizons]
    supply_biomass = supply_df[supply_df["carrier"].isin(biomass_carriers)].copy()
    supply_biomass["carrier"] = pd.Categorical(supply_biomass["carrier"], categories=biomass_carriers, ordered=True)
    supply_stacked = supply_biomass.groupby("carrier")[year_columns].sum().loc[biomass_carriers].T
    supply_stacked = supply_stacked.apply(pd.to_numeric, errors="coerce")
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    supply_stacked.plot.area(ax=ax)
    ax.set_title("Biomass Generation by Carrier (DE)")
    ax.set_ylabel("Generation [MWh]")
    ax.set_xlabel("Year")
    ax.legend(title="Carrier", loc="upper left", bbox_to_anchor=(1.0, 1.0))
    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()

def create_total_unsustainable_solid_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None,
                                                   ax=None):
    year_columns = [str(y) for y in planning_horizons]
    biocrop_carriers = ["unsustainable solid biomass"]

    supply_biocrops = supply_df[supply_df["carrier"].isin(biocrop_carriers)]
    e_sum_max_biocrops = e_sum_max_df[e_sum_max_df["carrier"].isin(biocrop_carriers)]
    e_sum_min_biocrops = e_sum_min_df[e_sum_min_df["carrier"].isin(biocrop_carriers)]

    supply_total = pd.to_numeric(supply_biocrops[year_columns].sum(), errors="coerce")
    e_max = pd.to_numeric(e_sum_max_biocrops[year_columns].sum(), errors="coerce")
    e_min = pd.to_numeric(e_sum_min_biocrops[year_columns].sum(), errors="coerce")
    #e_min += pd.to_numeric(globalconstraint_unsus_solid_biomass_min[year_columns])

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(year_columns, supply_total, label="Actual biocrop generation", color="darkred", linewidth=2)
    ax.set_title("Generation from Woody Biocrops (Unsustainable Solid Biomass)")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()

    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()


def create_DE_unsustainable_solid_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    biocrop_carriers = ["unsustainable solid biomass"]
    supply_biocrops = supply_df[supply_df["carrier"].isin(biocrop_carriers)]
    e_sum_max_biocrops = e_sum_max_df[e_sum_max_df["carrier"].isin(biocrop_carriers)]
    supply_total = pd.to_numeric(supply_biocrops[year_columns].sum(), errors="coerce")
    e_max = pd.to_numeric(e_sum_max_biocrops[year_columns].sum(), errors="coerce")
    e_min = e_min_usb_DE
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(year_columns, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(year_columns, supply_total, label="Actual biocrop generation", color="darkred", linewidth=2)
    ax.set_title("Unsustainable Solid Biomass (DE)")
    ax.set_ylabel("Generation [MWh]")
    ax.legend()
    if output_path:
        fig = ax.get_figure()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    # plt.show()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake
        snakemake = mock_snakemake("evaluate_biomass_run", run="ref-constr50")

    configure_logging(snakemake)

    e_sum_min_df = pd.read_csv(snakemake.input.e_sum_min_csv, index_col=0)
    e_sum_max_df = pd.read_csv(snakemake.input.e_sum_max_csv, index_col=0)
    supply_stats_df = pd.read_csv(snakemake.input.statistics_supply_csv, index_col=0)
    global_constraints_constant_df = pd.read_csv(snakemake.input.globalconstraints_constant_csv, index_col=0)
    planning_horizons = snakemake.params.planning_horizons
    biomass_config = snakemake.params.biomass

    globalconstraint_unsus_solid_biomass_min = global_constraints_constant_df.loc["unsustainable biomass min"]
    #e_min_usb_DE = calculate_e_min_unsus_solid_DE(e_sum_max_df, biomass_config, planning_horizons)
    e_min_usb_DE = 0 # if biomass_spatial is deactivated in config. else activate previous line instead.

    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(20, 10))
    axes = axes.flatten()

    # Total plots (obere Zeile)
    create_total_biomass_graph(e_sum_min_df, e_sum_max_df, supply_stats_df, snakemake.output.total_biomass, ax=axes[0])
    create_total_biomass_stack_graph(supply_stats_df, planning_horizons, snakemake.output.total_biomass_stack,
                                     ax=axes[1])
    create_total_biocrops_graph(e_sum_min_df, e_sum_max_df, supply_stats_df, snakemake.output.total_biocrops,
                                ax=axes[2])
    create_total_unsustainable_solid_biomass_graph(
        e_sum_min_df,
        e_sum_max_df,
        supply_stats_df,
        snakemake.output.total_unsustainable_solid_biomass,
        ax=axes[3]
    )

    # DE plots (untere Zeile)
    create_DE_biomass_graph(
        e_sum_min_df.query('country == "DE"'),
        e_sum_max_df.query('country == "DE"'),
        supply_stats_df.query('country == "DE"'),
        snakemake.output.DE_biomass,
        ax=axes[4]
    )
    create_DE_biomass_stack_graph(
        supply_stats_df.query('country == "DE"'),
        planning_horizons,
        snakemake.output.DE_biomass_stack,
        ax=axes[5]
    )
    create_DE_biocrops_graph(
        e_sum_min_df.query('country == "DE"'),
        e_sum_max_df.query('country == "DE"'),
        supply_stats_df.query('country == "DE"'),
        snakemake.output.DE_biocrops,
        ax=axes[6]
    )
    create_DE_unsustainable_solid_biomass_graph(
        e_sum_min_df.query('country == "DE"'),
        e_sum_max_df.query('country == "DE"'),
        supply_stats_df.query('country == "DE"'),
        snakemake.output.DE_unsustainable_solid_biomass,
        ax=axes[7]
    )

    plt.tight_layout()
    plt.show()
    fig.savefig(snakemake.output.raster_biomass, dpi=300, bbox_inches="tight")
