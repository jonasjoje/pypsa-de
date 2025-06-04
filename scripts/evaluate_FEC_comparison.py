import os
import re
import pypsa
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts._evaluation_helpers import load_networks_from_path_list, compare_value, plot_line_comparison, filter_statistics_by_country, load_csvs_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_FEC_comparison",
        )

    configure_logging(snakemake)

    logger.info("loading data")
    #nn = load_networks_from_path_list(snakemake.input.networks)
    cc_withdrawal = load_csvs_from_path_list(snakemake.input.statistics_withdrawal_csvs)

    planning_horizons = snakemake.params.planning_horizons


    # ─────────────────────────────────────────────────────────────────────────────
    #  Total FEC from statistics
    # ─────────────────────────────────────────────────────────────────────────────

    def get_total_FEC_in_twh(df_run):
        s = (
            df_run
            .loc[:, df_run.columns.intersection(map(str, planning_horizons))]
            .sum(axis=0)
        )
        return s * 1e-6
    plot_line_comparison(
                    cc = cc_withdrawal,
                    title="Total FEC in TWh",
                    expr=get_total_FEC_in_twh,
                    output=snakemake.output.total_FEC_graph)

    # ─────────────────────────────────────────────────────────────────────────────
    #  Sector FEC from statistics
    # ─────────────────────────────────────────────────────────────────────────────
    sector_map = {
        # Residential
        'electricity': 'residential',
        'rural heat': 'residential',
        'urban central heat': 'residential',
        'urban decentral heat': 'residential',
        # Transport domestic
        'land transport EV': 'transport domestic',
        'land transport fuel cell': 'transport domestic',
        'land transport oil': 'transport domestic',
        # Transport international
        'kerosene for aviation': 'transport international',
        'shipping methanol': 'transport international',
        'shipping oil': 'transport international',
        # Industry
        'industry electricity': 'industry',
        'low-temperature heat for industry': 'industry',
        'H2 for industry': 'industry',
        'coal for industry': 'industry',
        'gas for industry': 'industry',
        'industry methanol': 'industry',
        'naphtha for industry': 'industry',
        'solid biomass for industry': 'industry',
        'process emissions': 'industry',
        # Agriculture
        'agriculture electricity': 'agriculture',
        'agriculture heat': 'agriculture',
        'agriculture machinery oil': 'agriculture',
    }
    unique_sectors = set(sector_map.values())

    base_dir = os.path.dirname(snakemake.output.total_FEC_graph)


    def get_sector_FEC_series(df_run, sector):
        mask = df_run["carrier"].map(sector_map) == sector
        yrs = list(map(str, planning_horizons))
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(y) for y in s.index]
        return s.sort_index() * 1e-6


    for sector in unique_sectors:
        expr = lambda df_run, sector=sector: get_sector_FEC_series(df_run, sector)
        plot_line_comparison(
            cc=cc_withdrawal,
            title=f"{sector} FEC in TWh",
            expr=expr,
            output=os.path.join(base_dir, f"total_FEC_{sector}_graph.png"),
        )
        #print(f"{sector} done")

    # ─────────────────────────────────────────────────────────────────────────────
    #  DE FEC from statistics
    # ─────────────────────────────────────────────────────────────────────────────
    country = 'DE'

    def get_country_FEC_series(df_run, country):
        mask = df_run["country"] == country
        yrs = list(map(str, planning_horizons))
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(y) for y in s.index]
        return s.sort_index() * 1e-6

    expr = lambda df_run, country=country: get_country_FEC_series(df_run, country)
    plot_line_comparison(
        cc=cc_withdrawal,
        title=f"{country} FEC in TWh",
        expr=expr,
        output=snakemake.output.DE_FEC_graph,
    )


    # ─────────────────────────────────────────────────────────────────────────────
    #  DE Sector FEC from statistics
    # ─────────────────────────────────────────────────────────────────────────────
    def get_country_sector_FEC_series(df_run, sector, country):
        mask_country = df_run["country"] == country
        mask_sector = df_run.loc[mask_country, "carrier"].map(sector_map) == sector
        yrs = list(map(str, planning_horizons))
        sub = df_run.loc[mask_country & mask_sector, yrs]
        s = sub.sum(axis=0)
        s.index = [int(y) for y in s.index]
        return s.sort_index() * 1e-6

    country = 'DE'

    for sector in unique_sectors:
        expr = lambda df_run, sector=sector, country=country: get_country_sector_FEC_series(df_run, sector, country)

        plot_line_comparison(
            cc=cc_withdrawal,
            title=f"{country} {sector} FEC in TWh",
            expr=expr,
            output=os.path.join(base_dir, f"{country}_FEC_{sector}_graph.png"),
        )
        print(f"{sector} done")


    logger.info("All FEC comparisons done.")