import os
import re
import pypsa

from scripts._evaluation_helpers import load_networks_from_path_list

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "general_scenario_comparison",
        )

    nn = load_networks_from_path_list(snakemake.input.networks)

    # Zusammenfassung in test.txt schreiben
    with open(snakemake.output.test, "w") as f:
        for run, years in nn.items():
            for year in sorted(years):
                f.write(f"{run}_{year}\n")


    print("fertig")