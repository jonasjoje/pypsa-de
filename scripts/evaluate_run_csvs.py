import os
import pypsa
import logging
from scripts._evaluation_helpers import load_networks_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

def space_requirement_csvs(network_list):
    dpath = {year: next(p for p in network_list if str(year) in p)
             for year in planning_horizons}
    df_list = []
    for year, path in dpath.items():
        print(f"Lade Netzwerk für {year}: {path}")
        n = pypsa.Network(path)
        # Wähle Generator-Attribute aus
        df_year = n.generators[['space_req_opt', 'bus', 'carrier']].copy()
        # Extrahiere Land aus dem Bus-Namen (Beispiel: letzte Komponente nach '_')
        df_year['country'] = df_year['bus'].apply(lambda b: b.split('_')[-1])
        # Spalten für dieses Jahr mit MultiIndex (Attribut, Jahr)
        df_year.columns = pd.MultiIndex.from_product(
            [df_year.columns, [year]],
            names=['attribute', 'year']
        )
        df_list.append(df_year)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_run_csvs",
            run = "reference",
        )

        planning_horizons = snakemake.params.planning_horizons

        configure_logging(snakemake)

        space_requirement_csvs(snakemake.input.network_list)