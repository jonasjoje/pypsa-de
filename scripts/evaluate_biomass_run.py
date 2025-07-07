#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import logging
from scripts._helpers import configure_logging

# Set global font size
plt.rcParams['font.size'] = 7


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
    carriers = ["biogas", "unsustainable biogas", "solid biomass", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply = supply_df[supply_df["carrier"].isin(carriers)]

    # Convert to TWh
    e_max = pd.to_numeric(
        e_sum_max_df[e_sum_max_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    e_min = pd.to_numeric(
        e_sum_min_df[e_sum_min_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    supply_total = pd.to_numeric(supply[year_columns].sum(), errors="coerce") / 1e6

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    xs = [int(y) for y in year_columns]
    ax.fill_between(xs, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(xs, supply_total, label="Actual biomass generation", linewidth=2)
    ax.set_title("Biomass\nAll Countries", fontsize=7)

    # Disable scientific notation on y-axis
    ax.ticklabel_format(style='plain', axis='y')
    # Legend below plot
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


def create_total_biomass_stack_graph(supply_df, planning_horizons, output_path=None, ax=None):
    carriers = ["unsustainable solid biomass", "unsustainable biogas", "unsustainable bioliquids", "solid biomass", "biogas"]
    year_columns = [str(y) for y in planning_horizons]
    supply = supply_df[supply_df["carrier"].isin(carriers)].copy()
    supply["carrier"] = pd.Categorical(supply["carrier"], categories=carriers, ordered=True)

    stacked = (
        supply.groupby("carrier")[year_columns]
        .sum()
        .loc[carriers]
        .T
        .apply(pd.to_numeric, errors="coerce")
        / 1e6  # TWh
    )

    # Convert index to int for proper ticks
    stacked.index = stacked.index.astype(int)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    stacked.plot.area(ax=ax)
    ax.set_title("Biomass by Carrier\nAll Countries", fontsize=7)

    # Disable scientific notation
    ax.ticklabel_format(style='plain', axis='y')

    # Legend below plot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(labels))

    # Only show x-ticks for specific years
    ax.set_xticks([2020, 2030, 2040, 2050])
    ax.set_xticklabels(['2020', '2030', '2040', '2050'])

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


def create_total_biocrops_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    carriers = ["unsustainable biogas", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply = supply_df[supply_df["carrier"].isin(carriers)]

    e_max = pd.to_numeric(
        e_sum_max_df[e_sum_max_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    e_min = pd.to_numeric(
        e_sum_min_df[e_sum_min_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    supply_total = pd.to_numeric(supply[year_columns].sum(), errors="coerce") / 1e6

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    xs = [int(y) for y in year_columns]
    ax.fill_between(xs, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(xs, supply_total, label="Actual biocrop generation", linewidth=2)
    ax.set_title("Biocrops\nAll Countries", fontsize=7)

    ax.ticklabel_format(style='plain', axis='y')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


def create_total_unsustainable_solid_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    carriers = ["unsustainable solid biomass"]
    supply = supply_df[supply_df["carrier"].isin(carriers)]

    e_max = pd.to_numeric(
        e_sum_max_df[e_sum_max_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    e_min = pd.to_numeric(
        e_sum_min_df[e_sum_min_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    supply_total = pd.to_numeric(supply[year_columns].sum(), errors="coerce") / 1e6

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    xs = [int(y) for y in year_columns]
    ax.fill_between(xs, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(xs, supply_total, label="Actual biocrop generation", linewidth=2)
    ax.set_title("Woody Biocrops\nAll Countries", fontsize=7)

    ax.ticklabel_format(style='plain', axis='y')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


def create_DE_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    carriers = ["biogas", "unsustainable biogas", "solid biomass", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply = supply_df[supply_df["carrier"].isin(carriers)]

    e_max = pd.to_numeric(
        e_sum_max_df[e_sum_max_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    e_min = pd.to_numeric(
        e_sum_min_df[e_sum_min_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    supply_total = pd.to_numeric(supply[year_columns].sum(), errors="coerce") / 1e6

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    xs = [int(y) for y in year_columns]
    ax.fill_between(xs, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(xs, supply_total, label="Actual generation", linewidth=2)
    ax.set_title("Biomass\nGermany", fontsize=7)

    ax.ticklabel_format(style='plain', axis='y')
    ax.legend(loc='upper left', bbox_to_anchor=(-0.2, -0.2), ncol=1)

    ax.set_xticks([2020, 2030, 2040, 2050])
    ax.set_xticklabels(['2020', '2030', '2040', '2050'])

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


def create_DE_biomass_stack_graph(supply_df, planning_horizons, output_path=None, ax=None):
    carriers = ["unsustainable solid biomass", "unsustainable biogas", "unsustainable bioliquids", "solid biomass", "biogas"]
    year_columns = [str(y) for y in planning_horizons]
    supply = supply_df[supply_df["carrier"].isin(carriers)].copy()
    supply["carrier"] = pd.Categorical(supply["carrier"], categories=carriers, ordered=True)

    stacked = (
        supply.groupby("carrier")[year_columns]
        .sum()
        .loc[carriers]
        .T
        .apply(pd.to_numeric, errors="coerce")
        / 1e6  # TWh
    )
    stacked.index = stacked.index.astype(int)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    stacked.plot.area(ax=ax)
    ax.set_title("Biomass by Carrier\nGermany", fontsize=7)
    ax.ticklabel_format(style='plain', axis='y')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(-0.1, -0.2), ncol=1) #len(labels))

    ax.set_xticks([2020, 2030, 2040, 2050])
    ax.set_xticklabels(['', '', '', ''])

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


def create_DE_biocrops_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    carriers = ["unsustainable biogas", "unsustainable solid biomass", "unsustainable bioliquids"]
    supply = supply_df[supply_df["carrier"].isin(carriers)]

    e_max = pd.to_numeric(
        e_sum_max_df[e_sum_max_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    e_min = pd.to_numeric(
        e_sum_min_df[e_sum_min_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    supply_total = pd.to_numeric(supply[year_columns].sum(), errors="coerce") / 1e6

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    xs = [int(y) for y in year_columns]
    ax.fill_between(xs, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(xs, supply_total, label="Actual generation", linewidth=2)
    ax.set_title("Biocrops\nGermany", fontsize=7)

    ax.ticklabel_format(style='plain', axis='y')
    #ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1)

    ax.set_xticks([2020, 2030, 2040, 2050])
    ax.set_xticklabels(['', '', '', ''])

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


def create_DE_unsustainable_solid_biomass_graph(e_sum_min_df, e_sum_max_df, supply_df, output_path=None, ax=None):
    year_columns = [str(y) for y in planning_horizons]
    carriers = ["unsustainable solid biomass"]
    supply = supply_df[supply_df["carrier"].isin(carriers)]

    e_max = pd.to_numeric(
        e_sum_max_df[e_sum_max_df["carrier"].isin(carriers)][year_columns].sum(),
        errors="coerce"
    ) / 1e6
    e_min = calculate_e_min_unsus_solid_DE(e_sum_max_df, biomass_config, planning_horizons) / 1e6
    supply_total = pd.to_numeric(supply[year_columns].sum(), errors="coerce") / 1e6

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46/2))

    xs = [int(y) for y in year_columns]
    ax.fill_between(xs, e_min, e_max, color="gray", alpha=0.3, label="Potential range")
    ax.plot(xs, supply_total, label="Actual biocrop generation", linewidth=2)
    ax.set_title("Woody Biocrops\nGermany", fontsize=7)

    ax.ticklabel_format(style='plain', axis='y')
    #ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1)

    ax.set_xticks([2020, 2030, 2040, 2050])
    ax.set_xticklabels(['', '', '', ''])

    if output_path:
        ax.get_figure().savefig(output_path, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake
        snakemake = mock_snakemake("evaluate_biomass_run", run="reference")

    configure_logging(snakemake)

    # Data loading
    e_sum_min_df = pd.read_csv(snakemake.input.e_sum_min_csv, index_col=0)
    e_sum_max_df = pd.read_csv(snakemake.input.e_sum_max_csv, index_col=0)
    supply_stats_df = pd.read_csv(snakemake.input.statistics_supply_csv, index_col=0)
    planning_horizons = snakemake.params.planning_horizons
    biomass_config = snakemake.params.biomass

    # Compute DE-specific minimum
    e_min_usb_DE = calculate_e_min_unsus_solid_DE(e_sum_max_df, biomass_config, planning_horizons) if biomass_config.get('biomass_spatial', True) else 0

    # Create 2x4 grid, share y-axis within rows
    fig, axes = plt.subplots(
        nrows=2,
        ncols=4,
        figsize=(5.46, 4.5),
        sharey='row',
        sharex='col',
    )
    axes = axes.flatten()

    # Upper row: total plots (no legends)
    create_total_biomass_graph(e_sum_min_df, e_sum_max_df, supply_stats_df, snakemake.output.total_biomass, ax=axes[0])
    create_total_biomass_stack_graph(supply_stats_df, planning_horizons, snakemake.output.total_biomass_stack, ax=axes[1])
    create_total_biocrops_graph(e_sum_min_df, e_sum_max_df, supply_stats_df, snakemake.output.total_biocrops, ax=axes[2])
    create_total_unsustainable_solid_biomass_graph(e_sum_min_df, e_sum_max_df, supply_stats_df, snakemake.output.total_unsustainable_solid_biomass, ax=axes[3])

    # Lower row: DE plots (with legends)
    create_DE_biomass_graph(e_sum_min_df.query('country == "DE"'), e_sum_max_df.query('country == "DE"'), supply_stats_df.query('country == "DE"'), snakemake.output.DE_biomass, ax=axes[4])
    create_DE_biomass_stack_graph(supply_stats_df.query('country == "DE"'), planning_horizons, snakemake.output.DE_biomass_stack, ax=axes[5])
    create_DE_biocrops_graph(e_sum_min_df.query('country == "DE"'), e_sum_max_df.query('country == "DE"'), supply_stats_df.query('country == "DE"'), snakemake.output.DE_biocrops, ax=axes[6])
    create_DE_unsustainable_solid_biomass_graph(e_sum_min_df.query('country == "DE"'), e_sum_max_df.query('country == "DE"'), supply_stats_df.query('country == "DE"'), snakemake.output.DE_unsustainable_solid_biomass, ax=axes[7])

    # Remove legends for upper row
    for ax in axes[:4]:
        if ax.get_legend():
            ax.get_legend().remove()

    # Only y-tick labels on left column
    for idx, ax in enumerate(axes):
        if idx % 4 != 0:
            ax.tick_params(labelleft=False)

    # Y-label on left column, X-label on bottom row
    for ax in axes:
        ax.set_xlabel('')
        ax.set_ylabel('')
    for idx in (0, 4):
        axes[idx].set_ylabel('Generation in TWh')
    for ax in axes[4:]:
        ax.set_xlabel('')

    fig.text(
        0.95, 0.05,
        f"Scenario: {snakemake.wildcards.run}",
        ha='right', va='bottom',
        fontsize=7,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.5)
    )

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4, wspace=0.2, bottom=0.15)
    fig.savefig(snakemake.output.raster_biomass, dpi=300, bbox_inches="tight")
