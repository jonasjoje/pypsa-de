#!/usr/bin/env python3
import yaml
import pandas as pd
import matplotlib.pyplot as plt

# 1) Hard-gecodete Liste mit Deinen Prozent-Szenarien
PERCENTS = [0.5, 0.25, 0.01]  # z.B. 50%, 25%, 1%

# 2) Referenzwerte für 2020
REFERENCE_2020 = {
    "AT": 24488350460.855003,
    "BE": 7082851499.127298,
    "CH": 40923400.566239,
    "CZ": 18277166048.033268,
    "DE": 80970593113.51149,
    "DK": 8124059585.565379,
    "ES": 28928916463.091045,
    "FR": 57422254142.49341,
    "GB": 28680213620.23215,
    "IT": 41332969667.49164,
    "LU": 604067339.0603615,
    "NL": 10431340832.62824,
    "NO": 4950592377.9087925,
    "PL": 47391021188.88941,
    "SE": 48804456986.293304,
}

YEARS = [2030, 2040, 2050]

def interpolate(v0, v50, year):
    """Lineare Interpolation: 2020->v0, 2050->v50."""
    if year == 2020:
        return v0
    if year == 2050:
        return v50
    t = (year - 2020) / 30.0
    return v0 + (v50 - v0) * t

def build_scenario_block(percent):
    """Erzeugt den max_limit-Block für ein einzelnes Prozent-Szenario."""
    block = {}
    for country, v0 in REFERENCE_2020.items():
        v50 = v0 * percent
        block[country] = {
            year: interpolate(v0, v50, year)
            for year in YEARS
        }
    return block

def main():
    scenarios = {}
    for p in PERCENTS:
        label = f"{int(p*100)}%"
        scenarios[label] = {
            "land_use_module": {
                "types": {
                    "DLU": {
                        "constraint": {
                            "enable": True,
                            "max_limit": build_scenario_block(p)
                        }
                    }
                }
            }
        }

    # 3) Schreibe alles in eine YAML-Datei
    with open("../resources/space_req_max_limits.yaml", "w") as f:
        yaml.dump(scenarios, f, sort_keys=False)

    print("max_limits.yaml erzeugt mit Szenarien:", ", ".join(scenarios.keys()))


    ### for testing
    # land = "FR"
    #
    # # scenarios ist Dein Dict aus dem Script
    # data = {
    #     label: {
    #         year: scenarios[label]["land_use_module"]["types"]["DLU"]["constraint"]["max_limit"][land][year]
    #         for year in YEARS
    #     }
    #     for label in scenarios
    # }
    #
    # # DataFrame: Zeilen = Szenario-Labels (z.B. "50%"), Spalten = Jahre
    # df = pd.DataFrame.from_dict(data, orient="index")
    # print(df)
    #
    # # Für den Plot transponieren, damit Jahre auf der x-Achse stehen
    # df.T.plot(marker="o")
    # plt.xlabel("Year")
    # plt.ylabel("Max limit")
    # plt.title(f"{land} max_limit across scenarios")
    # plt.legend(title="Scenario")
    # plt.show()
    # print("testplot")

if __name__ == "__main__":
    main()
