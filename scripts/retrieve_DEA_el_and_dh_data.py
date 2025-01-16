# todo: documentation
# retrieve from data set from danish energy agency
# https://ens.dk/en/analyses-and-statistics/technology-data-generation-electricity-and-district-heating

# prio 1
# todo: als snakemake skript
# prio 2
# todo: vergleichen mit anderen retrieve skripts
# prio 3
# todo: verallgemeinerung auf andere Datensätze?
# todo: aktuellste version des Datensatzes erkennen

import urllib.request
import os
import pandas as pd

# URL of the latest Excel file
download_url = "https://ens.dk/media/5795/download"

# Ensure the script's directory is used as the base path
script_dir = os.path.dirname(os.path.abspath(__file__))
path_out = os.path.join(script_dir, "..", "data", "DEA_electricity_district_heat_data_sheet.xlsx")

# Ensure the output directory exists
#os.makedirs(os.path.dirname(path_out), exist_ok=True)

# Download the file only if it does not already exist
if not os.path.exists(path_out):
    urllib.request.urlretrieve(download_url, path_out)
    print(f"File downloaded and saved to {path_out}")
else:
    print(f"File already exists at {path_out}")


# Extract rows with 'Space requirement' from the 'alldata_flat' sheet
# and bring in costs.csv form
if os.path.exists(path_out):
    # Read the specific sheet "alldata_flat"
    data = pd.read_excel(path_out, sheet_name='alldata_flat')

    # Filter rows where the "par" column contains the string "Space requirement"
    # and "cat" column contains the value "Energy/technical data"
    sr = data[(data['par'].str.contains("Space requirement", na=False)) & (data['cat'] == "Energy/technical data")]

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

    # Add empty columns
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

    resources_dir = os.path.join(script_dir, "..", "resources", "space_requirements")

    for est in est_values:
        est_dir = os.path.join(resources_dir, est_labels[est])
        os.makedirs(est_dir, exist_ok=True)

        est_df = sr[sr['est'] == est]
        years = est_df['year'].unique()
        for year in years:
            year_df = est_df[est_df['year'] == year]
            # Drop the 'year' and 'est' columns
            year_df = year_df.drop(columns=['year', 'est'])
            file_name = f"space_requirement_{year}.csv"
            file_path = os.path.join(est_dir, file_name)
            year_df.to_csv(file_path, index=False)
            print(f"Saved dataframe to {file_path}")

    print("Finished")