
import os
import pandas as pd
import logging
from pathlib import Path
import pypsa
from scripts._helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)

def get_reference_FEC(n):
    df = n.loads[['bus', 'p_set', 'carrier']].copy()
    df['FEC_const'] = df.p_set.clip(lower=0) * n.snapshot_weightings.generators.sum()
    df['FEC_t'] = n.loads_t.p_set.clip(lower=0).mul(n.snapshot_weightings.generators, axis=0).sum()
    df['FEC'] = df['FEC_const'] + df['FEC_t'].fillna(0)

    sector_map = {
        # Residential
        'electricity': 'Residential',
        'residential rural heat': 'Residential',
        'urban central heat': 'Residential',
        'residential urban decentral heat': 'Residential',
        # Domestic transport
        'land transport EV': 'Domestic transport',
        'land transport fuel cell': 'Domestic transport',
        'land transport oil': 'Domestic transport',
        # International transport
        'kerosene for aviation': 'International transport',
        'shipping methanol': 'International transport',
        'shipping oil': 'International transport',
        # Industry
        'industry electricity': 'Industry',
        'low-temperature heat for industry': 'Industry',
        'H2 for industry': 'Industry',
        'coal for industry': 'Industry',
        'gas for industry': 'Industry',
        'industry methanol': 'Industry',
        'naphtha for industry': 'Industry',
        'solid biomass for industry': 'Industry',
        # Tertiary
        'services rural heat': 'Tertiary',
        'services urban decentral heat': 'Tertiary',
        # Agriculture
        'agriculture electricity': 'Agriculture',
        'agriculture heat': 'Agriculture',
        'agriculture machinery oil': 'Agriculture',
    }

    df['sector'] = df['carrier'].map(sector_map)
    df['country'] = df['bus'].str[:2]
    df_FEC = (
        df
        .groupby(['country', 'sector'])['FEC']
        .sum()
        .unstack(fill_value=0)
    )

    df_FEC.drop(index='EU', inplace=True)

    return df_FEC

def get_FEC_ratios(countries, years):
    ratio_list = []
    for country in countries:
        file_code = 'UK' if country == 'GB' else country
        file_path = os.path.join(os.getcwd(),f"data/CLEVER/ChartData_{file_code}.xlsx")
        df = pd.read_excel(file_path, sheet_name='Chart 10', header=2)
        df.rename(columns={'Unnamed: 0': 'Year'}, inplace=True)
        df.set_index('Year', inplace=True)
        df = df.loc[years]

        df['Industry'] = df['Industry'] + df['Non-energy']
        df.drop(columns=['Maritime bunkers', 'Energy', 'Non-energy'], errors='ignore', inplace=True)
        df.rename(columns={'Aviation bunkers': 'International transport'}, inplace=True)

        df = df.div(df.iloc[0], axis=1)

        df.index = pd.MultiIndex.from_product(
            [[country], df.index],
            names=['Country', 'Year']
        )
        ratio_list.append(df)

    return pd.concat(ratio_list)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "get_FEC_reference",
            run="reference",
        )
        os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    years = snakemake.config["scenario"]["planning_horizons"]
    countries = snakemake.config["countries"]
    #countries = ['DE']

    ref_FEC = get_reference_FEC(pypsa.Network(snakemake.input.network))

    multi_df = get_FEC_ratios(countries, years)

    multi_df_FEC = multi_df.mul(
        ref_FEC,
        axis=0,
        level='Country'
    )

    year_dfs = {
        year: multi_df_FEC.xs(year, level='Year')
        for year in years
    }

    for year, out_path in zip(years, snakemake.output.FEC_files):
        df = year_dfs[year]
        df.to_csv(out_path, index=True)

    logger.info("Exported FEC references.")
