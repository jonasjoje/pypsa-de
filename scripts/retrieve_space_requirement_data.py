# todo: documentation
# retrieve from data set from danish energy agency
# https://ens.dk/en/analyses-and-statistics/technology-data-generation-electricity-and-district-heating
# brings data in form of costs.csv
# creates separate csv for every available year and for mean, optimist, pessimist

# prio 1
# todo: als snakemake skript: logging
# prio 2
# todo: vergleichen mit anderen retrieve skripts
# todo: mehr technologien mappen
# prio 3
# todo: verallgemeinerung auf andere Datensätze?
# todo: aktuellste version des Datensatzes erkennen
# todo: progressbar bei download

import urllib.request
import os
import pandas as pd
import logging

from scripts._helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake, get_scenarios

        snakemake = mock_snakemake(
            "retrieve_space_requirement_data",
            year="2020",
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

    # Download the file only if it does not already exist
    if not os.path.exists(snakemake.params.dea_sheet_path):
        urllib.request.urlretrieve(snakemake.params.url, snakemake.params.dea_sheet_path)
        print(f"File downloaded and saved to {snakemake.params.dea_sheet_path}")
    else:
        print(f"File already exists at {snakemake.params.dea_sheet_path}")

    # Extract rows with 'Space requirement' from the 'alldata_flat' sheet
    if os.path.exists(snakemake.params.dea_sheet_path):
        # Read the specific sheet "alldata_flat"
        data = pd.read_excel(snakemake.params.dea_sheet_path, sheet_name='alldata_flat')

        # Filter rows where the "par" column contains the string "Space requirement"
        # and "cat" column contains the value "Energy/technical data"
        sr = data[(data['par'].str.contains("Space requirement", na=False)) & (data['cat'] == "Energy/technical data")]

        # Adjust space requirements for onwind based on generating capacity
        generating_capacity = data[(data['par'] == "Generating capacity for one unit [MW_e]")]
        generating_capacity = generating_capacity.rename(columns={'val': 'capacity_value'})[['Technology', 'capacity_value', 'year', 'est']]

        # Merge generating capacity with space requirements for onwind
        sr = sr.merge(generating_capacity, on=['Technology', 'year', 'est'], how='left')
        sr.loc[sr['Technology'] == 'Onshore wind turbine, utility - renewable power - wind - large', 'val'] \
            = (50 * 50) / sr['capacity_value'] * 1e-3

        # Drop the merged capacity column
        sr = sr.drop(columns=['capacity_value'])

        # Drop unnecessary columns
        sr = sr.drop(columns=['ref', 'note', 'priceyear', 'cat', 'ws'])  # todo: priceyear wichtig?

        # Rename and reorder columns
        sr = sr.rename(columns={
            'Technology': 'technology',
            'par': 'parameter',
            'val': 'value',
            'unit': 'unit'
        })

        # Remove units from the parameter column
        sr['parameter'] = sr['parameter'].str.replace(r'\s\[.*\]$', '', regex=True)

        # Add columns
        sr['source'] = "Danish Energy Agency, technology_data_for_el_and_dh.xlsx"
        sr['further description'] = None
        sr['currency_year'] = None

        # Reorder columns
        sr = sr[
            ['technology', 'parameter', 'value', 'unit', 'est', 'year', 'source', 'further description', 'currency_year']]

        # Map technologies to standard names todo: map more
        technology_mapping = {
            'PV - renewable power - solar - utility-scale, ground mounted': 'solar',
            'Onshore wind turbine, utility - renewable power - wind - large': 'onwind'
        }

        # Apply the mapping and drop rows with no mapping
        sr['technology'] = sr['technology'].map(technology_mapping)
        sr = sr.dropna(subset=['technology'])

        # Filter for the specific year and only 'mean' values (ctrl)
        year = int(snakemake.wildcards.planning_horizons)
        year_df = sr[(sr['year'] == year) & (sr['est'] == 'ctrl')]

        if not year_df.empty:
            # Drop the 'year' and 'est' columns
            year_df = year_df.drop(columns=['year', 'est'])
            year_df.to_csv(snakemake.output.csv_file, index=False)
            print(f"Saved dataframe to {snakemake.output.csv_file}")

        print("Finished")