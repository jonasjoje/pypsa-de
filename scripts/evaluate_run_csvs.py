import os
import pypsa
import pandas as pd
import logging
from scripts._evaluation_helpers import load_networks_from_path_list
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_run_csvs",
            run = "ref-constr50",
        )

    configure_logging(snakemake)

    planning_horizons = snakemake.params.planning_horizons
    network_list = snakemake.input.network_list

    # --------------------------------------------------------------
    # B) Konfigurierbare DataFrame-Definitionen
    # --------------------------------------------------------------
    expr_withdrawal = lambda net: (
        net.statistics
        .withdrawal(groupby=["bus", "carrier"], comps=["Load"])
        .to_frame("value")
        .reset_index(level=0, drop=True)
        .reset_index()
        .assign(idx=lambda x: x["bus"] + "_" + x["carrier"])
        .set_index("idx")
    )
    dataframes_config = [
        {
            "name": "space_requirements_DLU",
            "index_func": lambda net: net.generators.index,
            "static": {
                "bus": lambda net: net.generators["bus"],
                "carrier": lambda net: net.generators["carrier"],
                # Beispiel für berechnetes Attribut
                "country": lambda net: net.generators["bus"].str[:2],
            },
            "variable": {
                "space_req_DLU_opt": lambda net: net.generators["space_req_DLU_opt"],
            },
            "fillna": 0,
            "output_path": snakemake.output.space_requirements_DLU_csv,
        },
        {
            "name": "statistics_withdrawal_load_buscarrier",
            # hier „expr_withdrawal(net)“ aufrufen, um das DataFrame zu kriegen:
            "index_func": lambda net: expr_withdrawal(net).index,
            "static": {
                "bus": lambda net: expr_withdrawal(net)["bus"],
                "carrier": lambda net: expr_withdrawal(net)["carrier"],
                "country": lambda net: expr_withdrawal(net)["bus"].str[:2],
            },
            "variable": {
                "s_nom": lambda net: expr_withdrawal(net)["value"],
            },
            "fillna": 0,
            "output_path": snakemake.output.statistics_withdrawal_load_buscarrier_csv,
        },
        # Hier kannst Du beliebig weitere DataFrame-Definitionen anhängen:
        # {
        #     "name": "another_table",
        #     "index_func": lambda net: net.lines.index,
        #     "static": {
        #         "bus0": lambda net: net.lines["bus0"],
        #         "bus1": lambda net: net.lines["bus1"],
        #     },
        #     "variable": {
        #         "s_nom": lambda net: net.lines["s_nom"],
        #     },
        #     "fillna": float("nan"),
        #     "output_path": snakemake.output.lines_overview_csv,
        # },
    ]

    # --------------------------------------------------------------
    # C) Mapping YEAR → Pfad finden (damit das korrekte Netz geladen wird)
    # --------------------------------------------------------------
    dpath = {
        str(year): next(p for p in network_list if str(year) in os.path.basename(p))
        for year in planning_horizons
    }

    # --------------------------------------------------------------
    # D) Speicherstrukturen für Zwischenergebnisse aufbauen
    # --------------------------------------------------------------
    # Für jedes DataFrame in dataframes_config speichern wir statische und variable Daten je Jahr.
    # static_data["space_requirements_DLU"]["2020"] ist später ein DataFrame mit statischen Spalten für 2020.
    static_data = {conf["name"]: {} for conf in dataframes_config}
    variable_data = {conf["name"]: {} for conf in dataframes_config}

    # --------------------------------------------------------------
    # E) Hauptschleife: JEDE NETZWERKDATEI einmalig laden, Werte extrahieren
    # --------------------------------------------------------------
    for year in planning_horizons:
        year_str = str(year)
        path = dpath[year_str]
        logger.info(f"Load network for year {year_str} from {path}")
        net = pypsa.Network(path)

        for conf in dataframes_config:
            name = conf["name"]
            idx = conf["index_func"](net)

            # 1) Statische Attribute aus dem Netzwerk auslesen:
            #    Wir bauen ein DataFrame mit Index = idx und allen Konfig-spalten
            df_stat = pd.DataFrame(
                {col: func(net) for col, func in conf["static"].items()},
                index=idx,
            )
            # Speichern in static_data[name][year_str]
            static_data[name][year_str] = df_stat

            # 2) Variable Attribute auslesen:
            #    Jedes Attribut wird als pd.Series (Index = idx) zurückgeliefert
            #    Wir speichern für jedes DataFrame und für jedes Jahr ein Dictionary
            #    variable_data[name][year_str] = Series (Index = idx)
            var_dict = {}
            for var_name, func in conf["variable"].items():
                series = func(net)
                # Sicherstellen, dass die Series genau auf idx ausgerichtet ist
                series = series.reindex(idx)
                var_dict[var_name] = series
            variable_data[name][year_str] = var_dict

        # Netz freigeben, damit nur ein Netzwerk im Speicher liegt
        del net

    # --------------------------------------------------------------
    # F) EIN DataFrame pro Konfiguration zusammenbauen (Wide-Format)
    # --------------------------------------------------------------
    final_dataframes = {}
    for conf in dataframes_config:
        name = conf["name"]
        fillval = conf["fillna"]

        # 1) Statische Daten über alle Jahre zusammenführen:
        #    Wir concat'en alle df_stat per Jahr untereinander und entfernen Duplikate basierend auf dem Index.
        static_years = [static_data[name][str(year)] for year in planning_horizons]
        df_static_all = pd.concat(static_years, axis=0)
        # Falls ein Index (z.B. Generator-ID) in mehreren Jahren existiert,
        # behalten wir den ersten Eintrag. (Alle statischen Werte sind über Jahre identisch.)
        df_static_unique = df_static_all[~df_static_all.index.duplicated(keep='first')]

        # 2) Variable Daten über alle Jahre im Wide-Format zusammenführen:
        var_year_dfs = []
        for year in planning_horizons:
            year_str = str(year)
            var_dict = variable_data[name][year_str]
            df_var = pd.DataFrame(var_dict)
            df_var = df_var.rename(columns=lambda c: year_str)
            var_year_dfs.append(df_var)
        df_var_wide = pd.concat(var_year_dfs, axis=1)

        # 3) Statisch + Variabel zusammenführen (outer join):
        #    df_static_unique hat Index = gesamter Satz aller statischen Einträge.
        #    df_var_wide hat Index = gesamter Satz aller Indizes (aus allen Jahren).
        #    Mit outer join kombinieren wir beide Index-Mengen.
        df_final = df_static_unique.join(df_var_wide, how="outer")

        # 4) Fehlende Werte auffüllen
        df_final = df_final.fillna(fillval)

        # 5) Index-Name setzen
        df_final.index.name = "index"

        final_dataframes[name] = df_final

    # Am Ende von Teil 2 haben wir in final_dataframes[namen] jeweils ein komplettes pd.DataFrame.
    # Z. B. final_dataframes["space_requirements_DLU"].columns =
    # ["bus_2020","carrier_2020","country_2020","bus_2030",…,"space_req_DLU_opt_2020",…]

    # --------------------------------------------------------------
    # 3. Teil: Exportieren aller DataFrames als CSV
    # --------------------------------------------------------------
    for conf in dataframes_config:
        name = conf["name"]
        outpath = conf["output_path"]
        df_to_save = final_dataframes[name]

        # Pfad anlegen, falls Verzeichnisse fehlen
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        df_to_save.to_csv(outpath)
        logger.info(f"Export DataFrame '{name}' as CSV → {outpath}")