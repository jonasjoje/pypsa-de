import os
import re
import pypsa
import logging
import matplotlib.pyplot as plt
from scripts._evaluation_helpers import load_networks_from_path_list, compare_value, plot_line_comparison
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



    logger.info("All FEC comparisons done.")