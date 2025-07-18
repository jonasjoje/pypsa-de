# retrieve from data set from danish energy agency
# https://ens.dk/en/analyses-and-statistics/technology-data-generation-electricity-and-district-heating
# brings data in form of costs.csv
# creates separate csv for every available year and for mean, optimist, pessimist

# todo: fix error with retrieve-parameter

import urllib.request
import os
import pandas as pd
import logging

from scripts._helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

def download_file(url, path, retrieve_flag):
    if retrieve_flag:
        if os.path.exists(path):
            logger.warning(f"File exists at {path}, will be overwritten")
        urllib.request.urlretrieve(url, path)
        logger.info(f"File downloaded and saved to {path}")
    else:
        if os.path.exists(path):
            logger.info(f"Retrieve disabled, using existing file at {path}")
        else:
            raise RuntimeError(f"File {path} does not exist and retrieve is disabled")


def adjust_space_requirement_values(sr):
    sr['val'] = pd.to_numeric(sr['val'], errors='coerce')

    # adjust to m2/MW
    sr.loc[sr['val'].notna(), 'val'] *= 1e3
    sr.loc[:, 'unit'] = sr['unit'].str.replace(r'^1000', '', regex=True)

    return sr


def adjust_onwind_values(sr, data):
    generating_capacity = data[data['par'] == "Generating capacity for one unit [MW_e]"]
    generating_capacity = generating_capacity.rename(columns={'val': 'capacity_value'})[
        ['Technology', 'capacity_value', 'year', 'est']]
    sr = sr.merge(generating_capacity, on=['Technology', 'year', 'est'], how='left')
    sr.loc[sr['Technology'] == 'Onshore wind turbine, utility - renewable power - wind - large', 'val'] = (
            (50 * 50) / sr['capacity_value'])
    sr = sr.drop(columns=['capacity_value'])
    return sr

def adjust_flh_calculation(sr, flh):
    flh['val'] = pd.to_numeric(flh['val'], errors='coerce')

    # compute to MW nominal power (DEA assumes MW_e as continuously provided power)
    merged = sr.merge(
        flh[['ws', 'est', 'year', 'val']],
        on=['ws', 'est', 'year'],
        how='left',
        suffixes=('', '_flh')
    )
    # Berechne den angepassten Wert nur für Zeilen, bei denen ein flh-Wert vorhanden ist:
    # val_neu = val_sr * 8760 / val_flh
    mask = merged['val_flh'].notna()
    merged.loc[mask, 'val'] = merged.loc[mask, 'val'] * 8760 / merged.loc[mask, 'val_flh']

    sr = merged.drop(columns=['val_flh'])

    return sr


def cleanup_dataframe(sr, mapping, year, power_specific_generators):
    # Drop unnecessary columns
    sr = sr.drop(columns=['ref', 'note', 'priceyear', 'cat', 'ws'])
    # Rename columns
    sr = sr.rename(columns={
        'Technology': 'technology',
        'par': 'parameter',
        'val': 'value',
        'unit': 'unit'
    })
    # Remove units from the parameter column
    sr['parameter'] = sr['parameter'].str.replace(r'\s\[.*\]$', '', regex=True)
    # Filter for the specific year and only 'ctrl' values
    sr = sr[(sr['year'] == int(year)) & (sr['est'] == 'ctrl')]
    # Check if all power-specific generators exist in mapping
    unmapped_generators = set(power_specific_generators) - set(mapping.values())
    if unmapped_generators:
        raise RuntimeError(f"These power-specific generators are not in the DEA-mapping: {unmapped_generators}")
    # Apply mapping, filter
    sr['technology'] = sr['technology'].map(mapping)
    sr = sr.loc[sr['technology'].isin(power_specific_generators)]
    # Drop 'year' and 'est' columns
    sr = sr.drop(columns=['year', 'est'])
    # Add extra columns
    sr['source'] = "Danish Energy Agency, technology_data_for_el_and_dh.xlsx"
    sr['further description'] = None
    sr['currency_year'] = None
    # Reorder columns
    sr = sr[['technology', 'parameter', 'value', 'unit', 'source', 'further description', 'currency_year']]
    return sr

def add_energy_specific_technologies(sr, energy_specific_generators):
    # Convert dictionary to DataFrame
    energy_df = pd.DataFrame([
        {
            'technology': tech,
            'parameter': 'Space requirement',
            'value': value,
            'unit': 'm2/MWh',
            'source': 'Manually defined in config',
            'further description': None,
            'currency_year': None
        }
        for tech, value in energy_specific_generators.items()
    ])

    # Append to the existing DataFrame
    return pd.concat([sr, energy_df], ignore_index=True)


def apply_overwrite_values(df, overwrite_config, planning_year):
    for tech, year_overrides in overwrite_config.items():
        if planning_year in year_overrides:
            override_val = year_overrides[planning_year]
            mask = df['technology'] == tech

            if mask.any():
                # Existierende Zeilen überschreiben
                df.loc[mask, 'value'] = override_val
                df.loc[mask, 'source'] = (
                    f"Config override: {override_val} m2/MW for {tech}"
                )
            else:
                # Neue Zeile anfügen
                new_row = {
                    'technology': tech,
                    'parameter': 'Space requirement',
                    'value': override_val,
                    'unit': 'm2/MW',
                    'source': f"Config override: {override_val} m2/MW for {tech}",
                    'further description': None,
                    'currency_year': None
                }
                df = pd.concat([df, pd.DataFrame([new_row])],
                               ignore_index=True)
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake, get_scenarios

        snakemake = mock_snakemake(
            "retrieve_space_requirement_data",
            planning_horizons="2030",
            space_req_type="dist",
            run="reference",
        )
        rootpath = ".."
    else:
        from scripts._helpers import get_scenarios

        rootpath = "."

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(snakemake.params.dea_sheet_path), exist_ok=True)

    # Download file depending on 'retrieve'
    download_file(snakemake.params.url, snakemake.params.dea_sheet_path, snakemake.params.retrieve)

    # Read the specific sheet "alldata_flat"
    data = pd.read_excel(snakemake.params.dea_sheet_path, sheet_name='alldata_flat')

    # Filter rows where the "par" column contains "Space requirement" and "cat" equals "Energy/technical data"
    sr = data[(data['par'].str.contains("Space requirement", na=False)) & (data['cat'] == "Energy/technical data")]

    # Adjust the space requirement values and unit strings
    sr = adjust_space_requirement_values(sr)
    # Overwrite wind values based on generating capacity
    sr = adjust_onwind_values(sr, data)

    #sr = adjust_flh_calculation(sr, flh)

    # Define mapping dictionary (todo: map more)
    technology_mapping = {
        'PV - renewable power - solar - utility-scale, ground mounted': 'solar',
        'Onshore wind turbine, utility - renewable power - wind - large': 'onwind'
    }
    # Cleanup, filter for year and mapping
    year = snakemake.wildcards.planning_horizons
    power_specific_generators = snakemake.params.power_specific_generators
    sr = cleanup_dataframe(sr, technology_mapping, year, power_specific_generators)

    # Overwrite power specific generators based on config
    overwrite_values = snakemake.params.overwrite_values
    if overwrite_values:
        planning_year = int(snakemake.wildcards.planning_horizons)
        sr = apply_overwrite_values(sr, overwrite_values, planning_year)

    # Add energy specific technologies
    if snakemake.params.energy_specific_generators:
        sr = add_energy_specific_technologies(sr, snakemake.params.energy_specific_generators)

    # Save final dataframe as CSV
    sr.to_csv(snakemake.output.csv_file, index=False)
    logger.info(f"Saved dataframe to {snakemake.output.csv_file}")