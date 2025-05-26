import os
import pypsa
import pandas as pd
import logging
from scripts._evaluation_helpers import load_networks_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_space_requirement_run",
            run = "reference",
        )

        configure_logging(snakemake)

        print(snakemake.input.space_requirements_DLU_csv)
        print(snakemake.output.map)
        print("done")