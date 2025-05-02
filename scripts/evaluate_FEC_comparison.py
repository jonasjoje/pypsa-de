import os
import re
import pypsa
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts._evaluation_helpers import load_networks_from_path_list, compare_value, plot_line_comparison, filter_statistics_by_country
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_FEC_comparison",
        )

    configure_logging(snakemake)

    logger.info("loading networks")
    nn = load_networks_from_path_list(snakemake.input.networks)


    # ─────────────────────────────────────────────────────────────────────────────
    #  Total FEC from statistics
    # ─────────────────────────────────────────────────────────────────────────────

    expr = lambda n: n.statistics.withdrawal(comps="Load").sum() * 1e-6
    plot_line_comparison(
                    nn = nn,
                    title="Total FEC in TWh",
                    expr=expr,
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

    def compute_sector_FEC(n, sector):
        df = n.statistics().loc["Load"].Withdrawal
        df_grouped = df.groupby(df.index.map(sector_map)).sum()
        value = df_grouped.loc[sector]*1e-6  # to TWh/a
        return value

    for sector in unique_sectors:
        print(sector)
        expr = lambda n: compute_sector_FEC(n, sector)
        plot_line_comparison(
            nn=nn,
            title=f"{sector} FEC in TWh",
            expr=expr,
            output=os.path.join(base_dir, f"total_FEC_{sector}_graph.png"),)
        print(f"{sector} done")

    # ─────────────────────────────────────────────────────────────────────────────
    #  DE FEC from statistics
    # ─────────────────────────────────────────────────────────────────────────────
    country = 'DE'

    expr = lambda n: filter_statistics_by_country(n, country).loc["Load"].Withdrawal.sum() * 1e-6
    plot_line_comparison(
                    nn = nn,
                    title=f"{country} FEC in TWh",
                    expr=expr,
                    output=snakemake.output.DE_FEC_graph)


    # ─────────────────────────────────────────────────────────────────────────────
    #  DE Sector FEC from statistics
    # ─────────────────────────────────────────────────────────────────────────────
    def compute_sector_FEC_country(n, sector, country):
        df = filter_statistics_by_country(n, country).loc["Load"].Withdrawal
        df_grouped = df.groupby(df.index.map(sector_map)).sum()
        value = df_grouped.loc[sector]*1e-6  # to TWh/a
        return value

    country = 'DE'

    for sector in unique_sectors:
        print(sector)
        expr = lambda n: compute_sector_FEC_country(n, sector, country)
        plot_line_comparison(
            nn=nn,
            title=f"{country} {sector} FEC in TWh",
            expr=expr,
            output=os.path.join(base_dir, f"{country}_FEC_{sector}_graph.png"),)
        print(f"{sector} done")


    logger.info("All FEC comparisons done.")