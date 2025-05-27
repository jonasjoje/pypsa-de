import os
import pypsa
import pandas as pd
import logging
from scripts._evaluation_helpers import load_networks_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def space_requirements_DLU_csvs(network_list):
    # map year to its file path
    dpath = {y: next(p for p in network_list if str(y) in p)
             for y in planning_horizons}

    space_data = {}
    attr_records = []

    # load each network once, collect data
    for year, path in dpath.items():
        logger.info(f"Loading {year}: {path}")
        net = pypsa.Network(path)
        space_data[year] = net.generators['space_req_DLU_opt']
        attrs = net.generators[['bus', 'carrier']].copy()
        attrs['year'] = year
        attr_records.append(attrs)

    # build DataFrame of space requirements
    df_space = pd.DataFrame(space_data).fillna(0)
    df_space.index.name = 'generator_index'

    # combine and dedupe attributes
    df_attr = pd.concat(attr_records)
    df_attr = df_attr[~df_attr.index.duplicated(keep='first')]
    df_attr['country'] = df_attr['bus'].str[:2]
    static = df_attr[['bus', 'carrier', 'country']]

    # final DataFrame
    df = static.join(df_space, how='outer')
    return df

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_run_csvs",
            run = "ref-constr20",
        )

    configure_logging(snakemake)

    planning_horizons = snakemake.params.planning_horizons

    DLU_df = space_requirements_DLU_csvs(snakemake.input.network_list)
    DLU_df.to_csv(snakemake.output.space_requirements_DLU_csv)
    logger.info("Exported DLU csv to: snakemake.output.space_requirements_DLU_csv")