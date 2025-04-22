import os
import re
import pypsa
import logging
import matplotlib.pyplot as plt
from scripts._evaluation_helpers import load_networks_from_path_list, compare_value
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "general_scenario_comparison",
        )

    configure_logging(snakemake)

    logger.info("loading networks")
    nn = load_networks_from_path_list(snakemake.input.networks)

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEEX
    # ─────────────────────────────────────────────────────────────────────────────

    title = "CAPEX + OPEX (Bn. Euro)"
    logger.info(f"Create plot {title}")
    expr = lambda n: n.statistics.capex().sum() + n.statistics.opex().sum()
    df_objective = compare_value(expr, nn) * 1e-9

    ax = df_objective.plot.line()
    plt.ylabel(title)
    fig = ax.get_figure()
    fig.savefig(snakemake.output.objective_graph)
    logger.info(f"Created plot {title} and saved to {snakemake.output.objective_graph}")


    logger.info("All general scenario comparisons done.")