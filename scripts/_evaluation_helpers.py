import os
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger(__name__)

def load_networks_from_path_list(path_list):
    nn = {}
    for filepath in path_list:
        parts = filepath.split(os.sep)
        run = parts[parts.index("networks") - 1]
        year = int(filepath.split("_")[-1].split(".")[0])
        net = pypsa.Network(filepath)
        nn.setdefault(run, {})[year] = net
    return nn

def load_csvs_from_path_list(path_list):
    cc = {}
    for filepath in path_list:
        parts = filepath.split(os.sep)
        idx_csvs = parts.index("csvs")
        run = parts[idx_csvs - 1]
        df = pd.read_csv(filepath, index_col=0)
        cc[run] = df
    return cc

def compare_value_networks(get_value, nn):
    """
    Vergleicht einen bestimmten Wert aus den Netzwerken.

    Parameter:
      get_value: Funktion, die ein pypsa.Network-Objekt entgegennimmt und den gewünschten Wert zurückgibt.
      nn: Dictionary mit der Struktur { run_name: { planning_horizon: pypsa.Network } }

    Rückgabe:
      Ein DataFrame mit planning_horizons als Zeilen und run_names als Spalten.
    """
    results = {}
    for run_name, net_dict in nn.items():
        for planning_horizon, net in net_dict.items():
            value = get_value(net)
            if run_name not in results:
                results[run_name] = {}
            results[run_name][planning_horizon] = value
    df = pd.DataFrame(results)
    df = df.sort_index()  # Sortiere nach planning_horizons (Zeilen)
    return df

def compare_value(get_series, cc):
    """
    Vergleicht eine Zeitreihe, die aus jedem Run-DataFrame extrahiert wird.

    Parameter:
      get_series: Funktion, die ein DataFrame (für einen Run) entgegennimmt
                  und eine pandas.Series zurückgibt:
                    - Index = Jahr (z. B. [2020, 2025, 2030, …])
                    - Werte = die Kennzahl für dieses Jahr (float).

      cc: Dictionary { run_name: DataFrame, … }

    Rückgabe:
      DataFrame mit
        - Index = Jahr (z. B. 2020, 2025, 2030, …)
        - Spalten = run_names (z. B. "reference", "scenarioB", …)
      Jede Zelle = Wert, den get_series(df_run) für das jeweilige Jahr zurückliefert.
    """
    results = {}
    for run_name, df_run in cc.items():
        # get_series(df_run) muss eine Series liefern, deren Index Jahre sind:
        #   z. B. pd.Series(index=[2020,2025,2030], values=[..., ..., ...])
        series = get_series(df_run)
        results[run_name] = series

    # Aus dem Dict aus Series wird ein DataFrame: Index=Jahre, Spalten=run_names
    df_compare = pd.DataFrame(results)
    df_compare.sort_index(inplace=True)
    return df_compare


def plot_line_comparison(cc, title, expr, output):
    logger.info(f"Create  plot {title}")
    df = compare_value(expr, cc)
    ax = df.plot.line()
    plt.ylabel(title)
    fig = ax.get_figure()
    fig.savefig(output)
    logger.info(f"Created plot {title} and saved to {output}")

def filter_statistics_by_country(n, country):
    df = n.statistics(groupby=["bus","carrier"])
    df.index.set_names(["component", "bus", "carrier"], inplace=True)
    df=df.reset_index().query(f"bus.str.contains(@country)")
    df=df.groupby(["component","carrier"]).sum().drop(columns="bus")
    return df