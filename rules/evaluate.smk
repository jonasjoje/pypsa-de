# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import os
os.makedirs(".tmp", exist_ok=True)

EVALUATION = "results/" + run["prefix"] + "/EVALUATION/"

rule evaluate_all:
    input:
        expand(".tmp/{rule}.done", rule=config["evaluation"]["enable"]),
        ".tmp/evaluate_space_requirement_comparison.done",
        expand(".tmp/evaluate_general_run_{run}.done", run=config["run"]["name"]),
        expand(".tmp/evaluate_space_requirement_run_{run}.done", run=config["run"]["name"]),
        expand(".tmp/evaluate_biomass_run_{run}.done", run=config["run"]["name"]),
    output:
        touch(".tmp/all_evaluations.done")


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
        space_requirements_dist_csv = RESULTS + "csvs/space_requirements_dist.csv",
        statistics_withdrawal_load_buscarrier_csv = RESULTS + "csvs/statistics_withdrawal_load_buscarrier.csv",
        statistics_supply_generator_buscarrier_csv= RESULTS + "csvs/statistics_supply_generator_buscarrier.csv",
        statistics_capex_buscarrier_csv = RESULTS + "csvs/statistics_capex_buscarrier.csv",
        statistics_opex_buscarrier_csv = RESULTS + "csvs/statistics_opex_buscarrier.csv",
        statistics_optimalcapacity_buscarrier_csv= RESULTS + "csvs/statistics_optimalcapacity_buscarrier.csv",
        generators_p_nom_opt= RESULTS + "csvs/generators_p_nom_opt.csv",
        generators_p_nom_min= RESULTS + "csvs/generators_p_nom_min.csv",
        generators_p_nom_max= RESULTS + "csvs/generators_p_nom_max.csv",
        generators_e_sum_min = RESULTS + "csvs/generators_e_sum_min.csv",
        generators_e_sum_max = RESULTS + "csvs/generators_e_sum_max.csv",
        globalconstraints_constant = RESULTS + "csvs/globalconstraints_constant.csv",
        globalconstraints_mu = RESULTS + "csvs/globalconstraints_mu.csv",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        RESULTS + "logs/evaluate_run_csvs.log",
    script:
        "../scripts/evaluate_run_csvs.py"


GENERAL_COMPARISON = EVALUATION + "general_comparison/"
rule evaluate_general_scenario_comparison:
    params:
        planning_horizons=config_provider("scenario","planning_horizons"),
    input:
        statistics_capex_buscarrier_csv= expand(RESULTS + "csvs/statistics_capex_buscarrier.csv", run=config["run"]["name"]),
        statistics_opex_buscarrier_csv= expand(RESULTS + "csvs/statistics_opex_buscarrier.csv", run=config["run"]["name"]),
        statistics_optimalcapacity_buscarrier_csv= expand(RESULTS + "csvs/statistics_optimalcapacity_buscarrier.csv", run=config["run"]["name"]),
    output:
        #total_capexopex_graph = GENERAL_COMPARISON + "total_capexopex_graph.png",
        #DE_capexopex_graph= GENERAL_COMPARISON + "DE_capexopex_graph.png",
        total_and_DE_capexopex_graph= GENERAL_COMPARISON + "total_and_DE_capexopex_graph.png",
        total_and_DE_capexopex_stackplot_2050 = GENERAL_COMPARISON + "total_and_DE_capexopex_stackplot_2050.png",
        total_and_DE_capexopex_relative_change_2050 = GENERAL_COMPARISON + "total_and_DE_capexopex_relative_change_2050.png",
        total_and_DE_capexopex_absolute_change_2050 = GENERAL_COMPARISON + "total_and_DE_capexopex_absolute_change_2050.png",
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
        statistics_withdrawal_csvs = expand(
            RESULTS + "csvs/statistics_withdrawal_load_buscarrier.csv",
            run=config["run"]["name"]
        )
    output:
        total_and_DE_FEC_graph = GENERAL_COMPARISON + "total_and_DE_FEC_graph.png",
        done = touch(".tmp/evaluate_FEC_comparison.done")
    resources:
        mem = 30000
    log:
        EVALUATION + "logs/evaluate_FEC_comparison.log"
    script:
        "../scripts/evaluate_FEC_comparison.py"

rule evaluate_space_requirement_comparison:
    params:
        planning_horizons=config_provider("scenario","planning_horizons"),
    input:
        space_requirements_DLU_csv = expand(RESULTS + "csvs/space_requirements_DLU.csv",run=config["run"]["name"]),
        space_requirements_dist_csv = expand(RESULTS + "csvs/space_requirements_dist.csv",run=config["run"]["name"]),
        constraint_mu_DLU_csv = expand(RESULTS + "csvs/globalconstraints_mu.csv",run=config["run"]["name"]),
    output:
        DLU_vs_constraint_stack = GENERAL_COMPARISON + "DLU_vs_constraint_stack.png",
        DLU_mu_vs_constraint = GENERAL_COMPARISON + "DLU_mu_vs_constraint.png",
        dist_vs_constraint_stack= GENERAL_COMPARISON + "dist_vs_constraint_stack.png",
        done = touch(".tmp/evaluate_space_requirement_comparison.done")
    resources:
        mem = 30000
    log:
        EVALUATION + "logs/evaluate_space_requirement_comparison.log"
    script:
        "../scripts/evaluate_space_requirement_comparison.py"



rule evaluate_general_run:
    params:
        planning_horizons = config_provider("scenario","planning_horizons"),
    input:
        capex_csv = RESULTS + "csvs/statistics_capex_buscarrier.csv",
        opex_csv= RESULTS + "csvs/statistics_opex_buscarrier.csv",
        generators_p_nom_max_csv = RESULTS + "csvs/generators_p_nom_max.csv",
        generators_p_nom_min_csv = RESULTS + "csvs/generators_p_nom_min.csv",
        generators_p_nom_ops_csv = RESULTS + "csvs/generators_p_nom_opt.csv",
        globalconstraints_constant_csv= RESULTS + "csvs/globalconstraints_constant.csv",
    output:
        total_and_DE_capex_opex_graph = RESULTS + "graphs/total_and_DE_capex_opex.png",
        total_and_DE_solar_capacity_graph = RESULTS + "graphs/total_and_DE_solar_capacity.png",
        total_and_DE_onwind_capacity_graph = RESULTS + "graphs/total_and_DE_onwind_capacity.png",
        done= touch(".tmp/evaluate_general_run_{run}.done"),
    log:
        RESULTS + "logs/evaluate_general_run_{run}.log"
    script:
        "../scripts/evaluate_general_run.py"


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
        RESULTS + "logs/evaluate_space_requirement_run_{run}.log"
    script:
        "../scripts/evaluate_space_requirement_run.py"

rule evaluate_biomass_run:
    params:
        planning_horizons = config_provider("scenario","planning_horizons"),
        biomass = config_provider("biomass")
    input:
        e_sum_min_csv = RESULTS + "csvs/generators_e_sum_min.csv",
        e_sum_max_csv = RESULTS + "csvs/generators_e_sum_max.csv",
        statistics_supply_csv = RESULTS + "csvs/statistics_supply_generator_buscarrier.csv",
        globalconstraints_constant_csv = RESULTS + "csvs/globalconstraints_constant.csv",
    output:
        total_biomass=RESULTS + "graphs/total_biomass.png",
        DE_biomass=RESULTS + "graphs/DE_biomass.png",
        total_biocrops=RESULTS + "graphs/total_biocrops.png",
        DE_biocrops=RESULTS + "graphs/DE_biocrops.png",
        total_biomass_stack=RESULTS + "graphs/total_biomass_stack.png",
        DE_biomass_stack=RESULTS + "graphs/DE_biomass_stack.png",
        total_unsustainable_solid_biomass=RESULTS + "graphs/total_unsustainable_solid_biomass.png",
        DE_unsustainable_solid_biomass=RESULTS + "graphs/DE_unsustainable_solid_biomass.png",
        raster_biomass = RESULTS + "graphs/raster_biomass.png",

        done = touch(".tmp/evaluate_biomass_run_{run}.done"),
    log:
        RESULTS + "logs/evaluate_biomass_run_{run}.log"
    script:
        "../scripts/evaluate_biomass_run.py"

