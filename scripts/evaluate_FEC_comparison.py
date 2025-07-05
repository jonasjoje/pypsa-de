#!/usr/bin/env python3
import os
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 7

from scripts._evaluation_helpers import (
    load_csvs_from_path_list,
)
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    # Snakemake stub
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake
        snakemake = mock_snakemake("evaluate_FEC_comparison")

    configure_logging(snakemake)
    logger.info("Loading data")

    # load withdrawal statistics
    cc_withdrawal = load_csvs_from_path_list(
        snakemake.input.statistics_withdrawal_csvs
    )

    # filter scenarios to only 'reference' and 'clever'
    desired = {"reference", "clever"}
    cc_withdrawal = {
        name: df for name, df in cc_withdrawal.items()
        if name in desired
    }

    planning_horizons = snakemake.params.planning_horizons
    years = list(map(str, planning_horizons))
    base_dir = os.path.dirname(snakemake.output.total_and_DE_FEC_graph)

    # helper functions
    def get_total_FEC(df):
        s = df.loc[:, df.columns.intersection(years)].sum(axis=0)
        s.index = [int(y) for y in s.index]
        return s.sort_index() * 1e-6

    def get_country_FEC(df, country='DE'):
        mask = df['country'] == country
        s = df.loc[mask, years].sum(axis=0)
        s.index = [int(y) for y in s.index]
        return s.sort_index() * 1e-6

    def get_sector_FEC(df, sector, sector_map):
        mask = df['carrier'].map(sector_map) == sector
        s = df.loc[mask, years].sum(axis=0)
        s.index = [int(y) for y in s.index]
        return s.sort_index() * 1e-6

    def get_country_sector_FEC(df, sector, country, sector_map):
        mask = (df['country'] == country) & (df['carrier'].map(sector_map) == sector)
        s = df.loc[mask, years].sum(axis=0)
        s.index = [int(y) for y in s.index]
        return s.sort_index() * 1e-6

    # sector map from original script
    sector_map = {
        'electricity': 'Residential',
        'residential rural heat': 'Residential',
        'urban central heat': 'Residential',
        'residential urban decentral heat': 'Residential',
        'land transport EV': 'Domestic transport',
        'land transport fuel cell': 'Domestic transport',
        'land transport oil': 'Domestic transport',
        'kerosene for aviation': 'International transport',
        'shipping methanol': 'International transport',
        'shipping oil': 'International transport',
        'industry electricity': 'Industry',
        'low-temperature heat for industry': 'Industry',
        'H2 for industry': 'Industry',
        'coal for industry': 'Industry',
        'gas for industry': 'Industry',
        'industry methanol': 'Industry',
        'naphtha for industry': 'Industry',
        'solid biomass for industry': 'Industry',
        'services rural heat': 'Tertiary',
        'services urban decentral heat': 'Tertiary',
        'agriculture electricity': 'Agriculture',
        'agriculture heat': 'Agriculture',
        'agriculture machinery oil': 'Agriculture',
    }
    unique_sectors = sorted(set(sector_map.values()))

    # Combined Total and DE FEC
    logger.info("Plotting combined Total and DE FEC")
    fig, axes = plt.subplots(1, 2, figsize=(5.46, 2.5),
                         sharey=True)   # Ratio 12:5

    for label, df in cc_withdrawal.items():
        total = get_total_FEC(df)
        de = get_country_FEC(df, country='DE')
        axes[0].plot(total.index, total.values, label=label)
        axes[1].plot(de.index, de.values, label=label)

    axes[0].set_title("All Countries")
    axes[1].set_title("Germany")

    for ax in axes:
        ax.set_xlabel("Year")
        ax.set_ylabel("FEC in TWh")
        ax.set_ylim(bottom=0)
        ax.set_xticks([2020, 2030, 2040, 2050])
        ax.grid(True, linestyle='--')

    # Gemeinsame Legende unten
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='lower center',
               ncol=len(labels),
               frameon=False,
               bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.25)  # Platz für die Legende
    plt.savefig(snakemake.output.total_and_DE_FEC_graph, dpi=300)
    plt.close(fig)

    # Combined Sector FEC plots (All Countries vs Germany)
    logger.info("Plotting combined sector FEC comparisons")
    for sector in unique_sectors:
        fig, axes = plt.subplots(1, 2, figsize=(5.46, 2.5),
                         sharey=True)

        for label, df in cc_withdrawal.items():
            tot_sec = get_sector_FEC(df, sector, sector_map)
            de_sec = get_country_sector_FEC(df, sector, 'DE', sector_map)
            axes[0].plot(tot_sec.index, tot_sec.values, label=label)
            axes[1].plot(de_sec.index, de_sec.values, label=label)

        axes[0].set_title(f"All Countries - {sector}")
        axes[1].set_title("Germany")

        for ax in axes:
            ax.set_xlabel("Year")
            ax.set_ylabel(f"{sector} FEC in TWh")
            ax.set_ylim(bottom=0)
            ax.set_xticks([2020, 2030, 2040, 2050])
            ax.grid(True, linestyle='--')

        # Gemeinsame Legende unten
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels,
                   loc='lower center',
                   ncol=len(labels),
                   frameon=False,
                   bbox_to_anchor=(0.5, 0.0))

        fig.tight_layout()
        fig.subplots_adjust(bottom=0.25)
        out_path = os.path.join(base_dir, f"FEC_comparison_{sector}.png")
        plt.savefig(out_path, dpi=300)
        plt.close(fig)
        logger.info(f"Saved combined plot for sector: {sector}")

    logger.info("All combined FEC comparisons done.")
