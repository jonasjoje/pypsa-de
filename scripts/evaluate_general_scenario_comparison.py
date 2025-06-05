import os
import re
import pypsa
import logging
import matplotlib.pyplot as plt
from scripts._evaluation_helpers import load_networks_from_path_list, load_csvs_from_path_list, compare_value, plot_line_comparison
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)




if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "evaluate_general_scenario_comparison",
        )

    configure_logging(snakemake)

    logger.info("loading data.")
    #nn = load_networks_from_path_list(snakemake.input.networks)
    cc_capex = load_csvs_from_path_list(snakemake.input.statistics_capex_buscarrier_csv)
    cc_opex = load_csvs_from_path_list(snakemake.input.statistics_opex_buscarrier_csv)
    cc_optimalcapacity = load_csvs_from_path_list(snakemake.input.statistics_optimalcapacity_buscarrier_csv)

    planning_horizons = snakemake.params.planning_horizons

    # ─────────────────────────────────────────────────────────────────────────────
    #  CAPEX + OPEX
    # ─────────────────────────────────────────────────────────────────────────────

    cc_capex_opex = {}
    for run, df_cap in cc_capex.items():
        df_ope = cc_opex[run]
        year_cols = [str(y) for y in planning_horizons]
        df_sum = df_cap[year_cols].copy()
        for yr in year_cols:
            df_sum[yr] = df_cap[yr] + df_ope[yr]
        cc_capex_opex[run] = df_sum


    def get_capex_opex_series(df_run):
        s = df_run.sum(axis=0)
        s.index = s.index.astype(int)
        return s.sort_index() * 1e-9

    plot_line_comparison(
        cc=cc_capex_opex,
        title="CAPEX + OPEX (Bn. Euro)",
        expr=get_capex_opex_series,
        output=snakemake.output.total_capexopex_graph,
    )

    # ─────────────────────────────────────────────────────────────────────────────
    #  Generation
    # ─────────────────────────────────────────────────────────────────────────────

    # Solar
    def get_solar_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Solar")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Solar capacity (GW)",
        expr=get_solar_capacity_series,
        output=snakemake.output.gen_solar_graph,
    )


    #  Onshore Wind
    def get_onwind_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Onshore Wind")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Onshore Wind capacity (GW)",
        expr=get_onwind_capacity_series,
        output=snakemake.output.gen_onwind_graph,
    )


    #  Offshore Wind (AC)
    def get_offwind_ac_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Offshore Wind (AC)")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Offshore Wind AC capacity (GW)",
        expr=get_offwind_ac_capacity_series,
        output=snakemake.output.gen_offwind_ac_graph,
    )


    #  Offshore Wind (DC)
    def get_offwind_dc_capacity_series(df_run):
        mask = (df_run["component"] == "Generator") & (df_run["carrier"] == "Offshore Wind (DC)")
        yrs = [str(y) for y in planning_horizons]
        s = df_run.loc[mask, yrs].sum(axis=0)
        s.index = [int(year) for year in s.index]
        return s.sort_index() * 1e-3
    plot_line_comparison(
        cc=cc_optimalcapacity,
        title="Offshore Wind DC capacity (GW)",
        expr=get_offwind_dc_capacity_series,
        output=snakemake.output.gen_offwind_dc_graph,
    )

    logger.info("All general scenario comparisons done.")