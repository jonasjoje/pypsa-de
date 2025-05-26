import os
import pypsa
import logging
from scripts._evaluation_helpers import load_networks_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def space_requirement_csvs(network_list):
    dpath = {year: next(p for p in network_list if str(year) in p)
             for year in planning_horizons}
    for path in network_list:
        print(f"Evaluating {path}")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_run_csvs",
            run = "reference",
        )

        planning_horizons = snakemake.params.planning_horizons

        configure_logging(snakemake)

        space_requirement_csvs(snakemake.input.network_list)