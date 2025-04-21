# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule evaluation_all:
    input:
        expand("tmp/{rule}.done", rule=config["evaluation"]["enable"])
    output:
        temp("tmp/all_evaluations.done")


rule test:
    input:
        data = "results/20250411_reference/reference/networks/base_s_adm__none_2050.nc"
    output:
        temp("tmp/test.done")

rule general_scenario_comparison:
    input:
        networks = expand(
            RESULTS + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            clusters=config["scenario"]["clusters"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=config["scenario"]["planning_horizons"],
            run=config["run"]["name"]
        )
    output:
        test = "results/" + run["prefix"] + "/EVALUATION/general_comparison/test.txt",
        done = temp("tmp/general_scenario_comparison.done")
    script:
        "../scripts/evaluate_general_scenario_comparison.py"