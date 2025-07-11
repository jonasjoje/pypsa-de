# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script modifies the industrial production values to match the FORECAST
model.

This includes
- Production|Non-Metallic Minerals|Cement
- Production|Steel
- Production|Chemicals|Ammonia
- Production|Chemicals|Methanol
- Production|Non-Ferrous Metals
- Production|Pulp and Paper
"""

import logging

logger = logging.getLogger(__name__)
import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        snakemake = mock_snakemake(
            "modify_industry_demand",
            simpl="",
            clusters=27,
            opts="",
            ll="v1.0",
            sector_opts="None",
            run="8Gt_Bal_v3",
            planning_horizons=2020,
        )

    configure_logging(snakemake)
    # leitmodell for industry demand
    leitmodell = "FORECAST v1.0"

    year = snakemake.input.industrial_production_per_country_tomorrow.split("_")[
        -1
    ].split(".")[0]
    if snakemake.params.db_name == "ariadne2_intern" and year == "2020":
        logger.warning(f"Assuming {leitmodell} uses 2021 as base year instead of 2020.")
        year = "2021"
    existing_industry = pd.read_csv(
        snakemake.input.industrial_production_per_country_tomorrow, index_col=0
    )

    # read in ariadne database
    ariadne = (
        pd.read_csv(
            snakemake.input.ariadne,
            index_col=["model", "scenario", "region", "variable", "unit"],
        )
        .loc[
            leitmodell,
            snakemake.params.reference_scenario,
            "Deutschland",
            :,
            "Mt/yr",
        ]
        .multiply(1000)
    )

    logger.info(
        "German industry demand before modification",
    )
    logger.info(
        existing_industry.loc[
            "DE",
            [
                "Cement",
                "Electric arc",
                "Integrated steelworks",
                "DRI + Electric arc",
                "Ammonia",
                "Methanol",
                "Pulp production",
                "Paper production",
                "Ceramics & other NMM",
            ],
        ],
    )

    # write Cement, Ammonia and Methanol directly to dataframe
    existing_industry.loc["DE", "Cement"] = ariadne.loc[
        "Production|Non-Metallic Minerals|Cement", year
    ]
    existing_industry.loc["DE", "Ammonia"] = ariadne.loc[
        "Production|Chemicals|Ammonia", year
    ]
    existing_industry.loc["DE", "Methanol"] = ariadne.loc[
        "Production|Chemicals|Methanol", year
    ]

    # get ratio of pulp and paper production
    pulp_ratio = existing_industry.loc["DE", "Pulp production"] / (
        existing_industry.loc["DE", "Pulp production"]
        + existing_industry.loc["DE", "Paper production"]
    )

    existing_industry.loc["DE", "Pulp production"] = (
        ariadne.loc["Production|Pulp and Paper", year] * pulp_ratio
    )
    existing_industry.loc["DE", "Paper production"] = ariadne.loc[
        "Production|Pulp and Paper", year
    ] * (1 - pulp_ratio)

    # non-metallic minerals
    existing_industry.loc["DE", "Ceramics & other NMM"] = (
        ariadne.loc["Production|Non-Metallic Minerals", year]
        - ariadne.loc["Production|Non-Metallic Minerals|Cement", year]
    )

    # get steel ratios from existing_industry
    steel = existing_industry.loc[
        "DE", ["Electric arc", "Integrated steelworks", "DRI + Electric arc"]
    ]
    ratio = steel / steel.sum()

    # multiply with steel production including primary and secondary steel since distinguishing is taken care of later
    existing_industry.loc[
        "DE", ["Electric arc", "Integrated steelworks", "DRI + Electric arc"]
    ] = (ratio * ariadne.loc["Production|Steel", year])

    logger.info("German demand after modification")
    logger.info(
        existing_industry.loc[
            "DE",
            [
                "Cement",
                "Electric arc",
                "Integrated steelworks",
                "DRI + Electric arc",
                "Ammonia",
                "Methanol",
                "Pulp production",
                "Paper production",
                "Ceramics & other NMM",
            ],
        ],
    )

    existing_industry.to_csv(
        snakemake.output.industrial_production_per_country_tomorrow
    )
