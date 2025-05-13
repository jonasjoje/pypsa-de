
import os
import pandas as pd
import logging
from pathlib import Path
import pypsa


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "get_FEC_reference",
            run="reference",
        )


    print("get FEC from reference scenario.")

    net_file = Path(snakemake.input.network)
    out_path = Path(snakemake.output.FEC_value)

    n = pypsa.Network(net_file)

    value = n.statistics().loc["Load"].Withdrawal.sum()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(f"{value}\n")

    print("FEC_reference.txt exported.")