# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import os
os.makedirs(".tmp", exist_ok=True)

EVALUATION = "results/" + run["prefix"] + "/EVALUATION/"

rule evaluate_all:
    input:
        expand(".tmp/{rule}.done", rule=config["evaluation"]["enable"]),
        expand(".tmp/evaluate_space_requirement_run_{run}.done", run=config["run"]["name"])
    output:
        touch(".tmp/all_evaluations.done")


rule test:
    input:
        data = "results/20250411_reference/reference/networks/base_s_adm__none_2050.nc"
    output:
        touch(".tmp/test.done")


GENERAL_COMPARISON = EVALUATION + "general_comparison/"
rule evaluate_general_scenario_comparison:
    params:
        planning_horizons=config_provider("scenario","planning_horizons"),
    input:
        # networks = expand(
        #     RESULTS + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
        #     clusters=config["scenario"]["clusters"],
        #     opts=config["scenario"]["opts"],
        #     sector_opts=config["scenario"]["sector_opts"],
        #     planning_horizons=config["scenario"]["planning_horizons"],
        #     run=config["run"]["name"]
        # ),
        statistics_capex_buscarrier_csv= expand(RESULTS + "csvs/statistics_capex_buscarrier.csv", run=config["run"]["name"]),
        statistics_opex_buscarrier_csv= expand(RESULTS + "csvs/statistics_opex_buscarrier.csv", run=config["run"]["name"]),
        statistics_optimalcapacity_buscarrier_csv= expand(RESULTS + "csvs/statistics_optimalcapacity_buscarrier.csv", run=config["run"]["name"]),
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
    params:
        planning_horizons = config_provider("scenario","planning_horizons"),
    input:
        networks = expand(
            RESULTS + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            clusters=config["scenario"]["clusters"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=config["scenario"]["planning_horizons"],
            run=config["run"]["name"]
        ),
        statistics_withdrawal_csvs = expand(
            RESULTS + "csvs/statistics_withdrawal_load_buscarrier.csv",
            run=config["run"]["name"]
        )
    output:
        total_FEC_graph = GENERAL_COMPARISON + "total_FEC_graph.png",
        DE_FEC_graph = GENERAL_COMPARISON + "DE_FEC_graph.png",
        done = touch(".tmp/evaluate_FEC_comparison.done")
    resources:
        mem = 30000
    log:
        EVALUATION + "logs/evaluate_FEC_comparison.log"
    script:
        "../scripts/evaluate_FEC_comparison.py"

rule evaluate_run_csvs:
    params:
        planning_horizons = config_provider("scenario","planning_horizons")
    input:
        network_list = expand(
            RESULTS
            + "networks/base_s_{clusters}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"],
            allow_missing=True,
        )
    output:
        space_requirements_DLU_csv = RESULTS + "csvs/space_requirements_DLU.csv",
        statistics_withdrawal_load_buscarrier_csv = RESULTS + "csvs/statistics_withdrawal_load_buscarrier.csv",
        statistics_capex_buscarrier_csv = RESULTS + "csvs/statistics_capex_buscarrier.csv",
        statistics_opex_buscarrier_csv = RESULTS + "csvs/statistics_opex_buscarrier.csv",
        statistics_optimalcapacity_buscarrier_csv= RESULTS + "csvs/statistics_optimalcapacity_buscarrier.csv",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        RESULTS + "logs/evaluate_run_csvs.log",
    script:
        "../scripts/evaluate_run_csvs.py"

rule evaluate_space_requirement_run:
    params:
        DLU_config=config_provider("land_use_module","types","DLU")
    input:
        space_requirements_DLU_csv = RESULTS + "csvs/space_requirements_DLU.csv"
    output:
        DLU_DE = RESULTS + "graphs/DLU_DE.png",
        # map = expand(RESULTS + "maps/space_requirement_map_{planning_horizons}.png",
        #     planning_horizons = config["scenario"]["planning_horizons"],
        #     allow_missing=True)
        done = touch(".tmp/evaluate_space_requirement_run_{run}.done")
    resources:
        mem_mb=10000,
    log:
        RESULTS + "logs/evaluate_space_requirement_run.log"
    script:
        "../scripts/evaluate_space_requirement_run.py"
