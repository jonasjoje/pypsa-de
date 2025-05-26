import os
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import logging
from scripts._evaluation_helpers import load_networks_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def space_req_DLU_DE(df_DLU_country):
    # todo: stackplot carrier + max-limit
    DLU_DE = df_DLU_country.loc["DE"]
    df_plot = DLU_DE.to_frame(name="DLU_DE")
    df_plot.index = df_plot.index.astype(int)

    if is_constrained:
        max_limits_DE = max_limits_df.loc["DE"]
        df_plot["max_limit_DE"] = max_limits_DE

    df_plot = df_plot * 1e-9
    df_plot.loc[2020, "max_limit_DE"] = df_plot.loc[2020, "DLU_DE"]

    total_area_de = 357.022  # 1000 km²

    fig, ax = plt.subplots(figsize=(8, 5))
    df_plot.plot(ax=ax, kind="line", ylim=(0, None))

    def abs_to_rel(x):
        return x / total_area_de

    def rel_to_abs(x):
        return x * total_area_de

    secax = ax.secondary_yaxis(
        "right",
        functions=(abs_to_rel, rel_to_abs)
    )
    secax.set_ylabel("% of country area")
    secax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))

    ax.set_xticks(df_plot.index.astype(int))
    ax.set_title("DLU DE")
    ax.set_ylabel("1000 km²")
    ax.legend()
    plt.tight_layout()
    #plt.show()
    fig.savefig(snakemake.output.DLU_DE,
                dpi=300,
                bbox_inches='tight')

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_space_requirement_run",
            run = "ref-constr50",
        )

        configure_logging(snakemake)
        DLU_config=snakemake.params.DLU_config
        is_constrained = DLU_config["constraint"]["enable"]

        df_DLU = pd.read_csv(snakemake.input.space_requirements_DLU_csv, index_col=0)
        df_DLU_country = df_DLU.drop(columns=['bus', 'carrier']).groupby('country', as_index=True).sum().drop(index="EU")


        if is_constrained:
            max_limits = DLU_config["constraint"].get("max_limit")
            max_limits_df = pd.DataFrame.from_dict(max_limits, orient='index').apply(pd.to_numeric)
            logger.info("constraints loaded.")
        else:
            logger.info("No space requirement constraints defined.")

        space_req_DLU_DE(df_DLU_country)