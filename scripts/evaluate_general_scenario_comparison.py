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
            "evaluate_general_scenario_comparison",
        )

    configure_logging(snakemake)

    logger.info("loading networks")
    nn = load_networks_from_path_list(snakemake.input.networks)

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX
    # ─────────────────────────────────────────────────────────────────────────────

    expr = lambda n: (n.statistics.capex().sum() + n.statistics.opex().sum()) * 1e-9
    plot_line_comparison(
                    nn = nn,
                    title="CAPEX + OPEX (Bn. Euro)",
                    expr=expr,
                    output=snakemake.output.total_capexopex_graph)

    # ─────────────────────────────────────────────────────────────────────────────
    #  Generation
    # ─────────────────────────────────────────────────────────────────────────────

    # Solar
    expr = lambda n: n.statistics.optimal_capacity().loc[("Generator", "Solar")] * 1e-3
    plot_line_comparison(
                    nn = nn,
                    title="Solar capacity (GW)",
                    expr=expr,
                    output=snakemake.output.gen_solar_graph)

    # Onwind
    expr = lambda n: n.statistics.optimal_capacity().loc[("Generator", "Onshore Wind")] * 1e-3
    plot_line_comparison(
                    nn = nn,
                    title="Onshore Wind capacity (GW)",
                    expr=expr,
                    output=snakemake.output.gen_onwind_graph)

    # Offshore Wind AC
    expr = lambda n: n.statistics.optimal_capacity().loc[("Generator", "Offshore Wind (AC)")] * 1e-3
    plot_line_comparison(
                    nn = nn,
                    title="Offshore Wind AC capacity (GW)",
                    expr=expr,
                    output=snakemake.output.gen_offwind_ac_graph
                )

    # Offshore Wind DC
    expr = lambda n: n.statistics.optimal_capacity().loc[("Generator", "Offshore Wind (DC)")] * 1e-3
    plot_line_comparison(
                    nn = nn,
                    title="Offshore Wind DC capacity (GW)",
                    expr=expr,
                    output=snakemake.output.gen_offwind_dc_graph
                )




    logger.info("All general scenario comparisons done.")