import os
import pypsa

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "general_scenario_comparison",
        )


    with open(snakemake.output.test, "w") as f:
        for p in snakemake.input.networks:
            f.write(p + "\n")


    print("fertig")