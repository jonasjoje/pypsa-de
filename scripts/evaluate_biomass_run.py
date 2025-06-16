import os
import pandas as pd
import matplotlib.pyplot as plt
import logging
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def evaluate_biomass(e_sum_min_df, e_sum_max_df, supply_stats_df, planning_horizons):

    # Initial checks
    logger.debug(f"e_sum_min_df shape: {e_sum_min_df.shape}")
    logger.debug(f"e_sum_max_df shape: {e_sum_max_df.shape}")
    logger.debug(f"supply_stats_df shape: {supply_stats_df.shape}")
    logger.debug(f"Planning horizons: {planning_horizons}")

    # TODO: Add plotting logic here
    logger.info("Plotting not yet implemented.")
    # Placeholder: empty plot
    fig, ax = plt.subplots()
    ax.set_title("Total Biomass Placeholder")
    plt.tight_layout()
    plt.show()
    fig.savefig(snakemake.output.total_biomass, dpi=300, bbox_inches='tight')

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

    # Call evaluation function
    evaluate_biomass(e_sum_min_df, e_sum_max_df, supply_stats_df, planning_horizons)
