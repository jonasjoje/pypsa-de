# todo: documentation
# retrieve from data set from danish energy agency
# https://ens.dk/en/analyses-and-statistics/technology-data-generation-electricity-and-district-heating

# todo: als snakemake skript
# todo: überprüfen ob (aktuellste) datei vorhanden ist
# todo: vergleichen mit anderen retrieve skripts

import urllib.request
import os

# URL of the latest Excel file
download_url = "https://ens.dk/media/5795/download"

# Ensure the script's directory is used as the base path
script_dir = os.path.dirname(os.path.abspath(__file__))
path_out = os.path.join(script_dir, "..", "data", "DEA_electricity_district_heat_data_sheet.xlsx")

# Ensure the output directory exists
#os.makedirs(os.path.dirname(path_out), exist_ok=True)

# Download and save the file
urllib.request.urlretrieve(download_url, path_out)

print('Download complete.')