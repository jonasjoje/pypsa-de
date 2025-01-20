# todo: documentation
# retrieve from data set from danish energy agency
# https://ens.dk/en/analyses-and-statistics/technology-data-generation-electricity-and-district-heating
# brings data in form of costs.csv
# creates separate csv for every available year and for mean, optimist, pessimist

# prio 1
# todo: als snakemake skript: wildcards, outputs, logging
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
            #cost_horizon="mean"
        )
        rootpath = ".."
    else:
        from scripts._helpers import get_scenarios
        rootpath = "."

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    # URL of the latest Excel file
    download_url = "https://ens.dk/media/5795/download"

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(snakemake.output.dea_sheet), exist_ok=True)

    # Download the file only if it does not already exist
    if not os.path.exists(snakemake.output.dea_sheet):
        urllib.request.urlretrieve(download_url, snakemake.output.dea_sheet)
        print(f"File downloaded and saved to {snakemake.output.dea_sheet}")
    else:
        print(f"File already exists at {snakemake.output.dea_sheet}")

    # Extract rows with 'Space requirement' from the 'alldata_flat' sheet
    if os.path.exists(snakemake.output.dea_sheet):
        # Read the specific sheet "alldata_flat"
        data = pd.read_excel(snakemake.output.dea_sheet, sheet_name='alldata_flat')

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
            'PV - renewable power - solar - utility-scale, ground mounted': 'solar-utility',
            'Onshore wind turbine, utility - renewable power - wind - large': 'onwind'
        }

        # Apply the mapping and drop rows with no mapping
        sr['technology'] = sr['technology'].map(technology_mapping)
        sr = sr.dropna(subset=['technology'])

        # Split the dataframe by 'est' and 'year'
        est_values = ['ctrl', 'upper', 'lower']
        est_labels = {'ctrl': 'mean', 'upper': 'optimist', 'lower': 'pessimist'}

        resources_dir = os.path.join(rootpath, "resources", "space_requirements")

        for est in est_values:
            est_dir = os.path.join(resources_dir, est_labels[est])
            os.makedirs(est_dir, exist_ok=True)

            est_df = sr[sr['est'] == est]
            for year in sr['year'].unique():
                year_df = est_df[est_df['year'] == year]
                if not year_df.empty:
                    # Drop the 'year' and 'est' columns
                    year_df = year_df.drop(columns=['year', 'est'])
                    file_name = f"space_requirement_{year}_{est_labels[est]}.csv"
                    file_path = os.path.join(est_dir, file_name)
                    year_df.to_csv(file_path, index=False)
                    print(f"Saved dataframe to {file_path}")

        print("Finished")