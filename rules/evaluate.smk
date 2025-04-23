# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import os
os.makedirs(".tmp", exist_ok=True)

EVALUATION = "results/" + run["prefix"] + "/EVALUATION/"

rule evaluate_all:
    input:
        expand(".tmp/{rule}.done", rule=config["evaluation"]["enable"])
    output:
        touch(".tmp/all_evaluations.done")


rule test:
    input:
        data = "results/20250411_reference/reference/networks/base_s_adm__none_2050.nc"
    output:
        touch(".tmp/test.done")


GENERAL_COMPARISON = EVALUATION + "general_comparison/"
rule evaluate_general_scenario_comparison:
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
        total_capexopex_graph = GENERAL_COMPARISON + "total_capexopex_graph.png",
        gen_solar_graph = GENERAL_COMPARISON +"gen_solar_graph.png",
        gen_onwind_graph = GENERAL_COMPARISON +"gen_onwind_graph.png",
        gen_offwind_ac_graph= GENERAL_COMPARISON +"gen_offwind_ac_graph.png",
        gen_offwind_dc_graph= GENERAL_COMPARISON +"gen_offwind_dc_graph.png",
        done = touch(".tmp/evaluate_general_scenario_comparison.done")
    resources:
        mem = 30000
    log:
        EVALUATION + "logs/evaluate_general_comparison.log"
    script:
        "../scripts/evaluate_general_scenario_comparison.py"

rule evaluate_FEC_comparison:
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
        total_FEC_graph = GENERAL_COMPARISON + "total_FEC_graph.png",
        done = touch(".tmp/evaluate_FEC_comparison.done")
    resources:
        mem = 30000
    log:
        EVALUATION + "logs/evaluate_FEC_comparison.log"
    script:
        "../scripts/evaluate_FEC_comparison.py"