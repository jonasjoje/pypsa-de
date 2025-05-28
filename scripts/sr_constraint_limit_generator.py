#!/usr/bin/env python3
import yaml
import pandas as pd
import matplotlib.pyplot as plt

# 1) Hard-gecodete Liste mit Deinen Prozent-Szenarien
PERCENTS = [0.5, 0.25, 0.20, 0.15, 0.10, 0.05, 0.01]  # z.B. 50%, 25%, 1%

# 2) Referenzwerte für 2020 – genau die von Dir genannten
REFERENCE_2020 = {
    "AT": 2.448579e10,
    "BE": 7.082864e9,
    "CH": 4.094392e7,
    "CZ": 1.827723e10,
    "DE": 9.913593e10,
    "DK": 8.120990e9,
    "ES": 2.891575e10,
    "FR": 5.739774e10,
    "GB": 2.885746e10,
    "IT": 2.588361e10,
    "LU": 6.040852e8,
    "NL": 7.595576e9,
    "NO": 4.950605e9,
    "PL": 4.737675e10,
    "SE": 4.880445e10,
}

YEARS = [2020, 2030, 2040, 2050]

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
