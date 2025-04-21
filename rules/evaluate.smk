# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule evaluation_all:
    input:
        expand("tmp/{rule}.done", rule=config["evaluation"]["enable"])
    output:
        "results/all_evaluations.done"


rule test:
    input:
        data = "results/20250411_reference/reference/networks/base_s_adm__none_2050.nc"
    output:
        temp("tmp/test.done")

# rule evaluate_space_requirement:
#     params:
#         planning_horizons=config_provider("scenario", "planning_horizons")
#     input:
#         expand(
#             RESULTS + "postnetworks/base_s_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
#             allow_missing=True,
#             **config["scenario"],
#         )
#     output:
#         report=directory(RESULTS + "space_requirement_summary")
#     script:
#         "../scripts/evaluate_space_requirement.py"