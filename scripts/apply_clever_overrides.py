#!/usr/bin/env python3
# scripts/apply_clever_overrides.py

"""
Load the standard energy, heat-share, CO₂ and transport tables,
and—if in a “Clever” run—apply all of the same Clever overrides
that used to live inline in build_energy_totals.py.
"""

import pandas as pd

# 1) Read in the “vanilla” outputs
energy_df          = pd.read_csv(snakemake.input.energy,          index_col=[0,1])
district_heat_df   = pd.read_csv(snakemake.input.district_heat,   index_col=0)
co2_df             = pd.read_csv(snakemake.input.co2,             index_col=0)
transport_df       = pd.read_csv(snakemake.input.transport,       index_col=[0,1])

is_clever  = snakemake.config["clever"]

if is_clever:
    # 2) Load all Clever inputs
    clever_residential = pd.read_csv(snakemake.input.clever_residential, index_col=0)
    clever_tertiary    = pd.read_csv(snakemake.input.clever_tertiary,    index_col=0)
    clever_transport   = pd.read_csv(snakemake.input.clever_transport,   index_col=0)
    clever_agriculture = pd.read_csv(snakemake.input.clever_agriculture, index_col=0)
    clever_macro       = pd.read_csv(snakemake.input.clever_macro,       index_col=0)
    clever_afolub      = pd.read_csv(snakemake.input.clever_afolub,      index_col=0)

    year = snakemake.params.energy["energy_totals_year"]
    countries = snakemake.params.countries

    # 3) Override energy_df exactly as before
    for country in countries:
        energy_df.loc[(country, year), 'total road']                  = clever_transport.loc[country, 'Total_Road']
        energy_df.loc[(country, year), 'electricity road']            = clever_transport.loc[country, 'Electricity_Road']
        energy_df.loc[(country, year), 'total passenger cars']        = clever_transport.loc[country, 'Total final energy consumption in passenger road mobility']
        energy_df.loc[(country, year), 'electricity passenger cars']  = clever_transport.loc[country, 'Final electricity consumption for passenger road mobility']
        energy_df.loc[(country, year), 'total rail']                  = clever_transport.loc[country, 'Total_rail']
        energy_df.loc[(country, year), 'electricity rail']            = clever_transport.loc[country, 'Electricity_rail']
        energy_df.loc[(country, year), 'total rail passenger']        = clever_transport.loc[country, 'Total final energy consumption in rail passenger transport']
        energy_df.loc[(country, year), 'electricity rail passenger']  = clever_transport.loc[country, 'Final electricity consumption in rail passenger transport']
        energy_df.loc[(country, year), 'total rail freight']          = clever_transport.loc[country, 'Total final energy consumption in rail freight transport']
        energy_df.loc[(country, year), 'electricity rail freight']    = clever_transport.loc[country, 'Final electricity consumption in rail freight transport']
        energy_df.loc[(country, year), 'total aviation passenger']    = clever_transport.loc[country, 'Total final energy consumption for air travel']
        energy_df.loc[(country, year), 'total international aviation'] = clever_transport.loc[country, 'Total final energy consumption for air travel']
        energy_df.loc[(country, year), 'total domestic navigation']    = clever_transport.loc[country, 'Final energy consumption from liquid fuels in national water freight transport']
        energy_df.loc[(country, year), 'total international navigation']= clever_transport.loc[country, 'Final energy consumption from liquid fuels in international water freight transport']

        energy_df.loc[(country, year), 'total residential space']      = clever_residential.loc[country, 'Total final energy consumption for space heating in the residential sector']
        energy_df.loc[(country, year), 'electricity residential space']= clever_residential.loc[country, 'Final electricity consumption for space heating in the residential sector']
        energy_df.loc[(country, year), 'total residential water']      = clever_residential.loc[country, 'Total final energy consumption for domestic hot water']
        energy_df.loc[(country, year), 'electricity residential water']= clever_residential.loc[country, 'Final electricity consumption for domestic hot water']
        energy_df.loc[(country, year), 'total residential cooking']    = clever_residential.loc[country, 'Total final energy consumption for domestic cooking']
        energy_df.loc[(country, year), 'electricity residential cooking']= clever_residential.loc[country, 'Final electricity consumption for domestic cooking']
        energy_df.loc[(country, year), 'total residential']            = clever_residential.loc[country, 'Total final energy consumption in the residential sector']
        energy_df.loc[(country, year), 'electricity residential']      = clever_residential.loc[country, 'Final electricity consumption in the residential sector']
        energy_df.loc[(country, year), 'derived heat residential']     = clever_residential.loc[country, 'Final energy consumption from heating networks in the residential sector']
        energy_df.loc[(country, year), 'thermal uses residential']     = clever_residential.loc[country, 'Thermal_uses_residential']

        energy_df.loc[(country, year), 'total services space']         = clever_tertiary.loc[country, 'Total final energy consumption for space heating in the tertiary sector (with climatic corrections) ']
        energy_df.loc[(country, year), 'electricity services space']   = clever_tertiary.loc[country, 'Final electricity consumption for space heating in the tertiary sectorr']
        energy_df.loc[(country, year), 'total services water']         = clever_tertiary.loc[country, 'Total final energy consumption for hot water in the tertiary sector']
        energy_df.loc[(country, year), 'electricity services water']   = clever_tertiary.loc[country, 'Final electricity consumption for hot water in the tertiary sector']
        energy_df.loc[(country, year), 'total services cooking']       = clever_tertiary.loc[country, 'Total Final energy consumption for cooking in the tertiary sector']
        energy_df.loc[(country, year), 'electricity services cooking'] = clever_tertiary.loc[country, 'Final electricity consumption for cooking in the tertiary sector']
        energy_df.loc[(country, year), 'total services']               = clever_tertiary.loc[country, 'Total final energy consumption in the tertiary sector']
        energy_df.loc[(country, year), 'electricity services']         = clever_tertiary.loc[country, 'Final electricity consumption in the tertiary sector']
        energy_df.loc[(country, year), 'derived heat services']        = clever_tertiary.loc[country, 'Final energy consumption from heating networks in the tertiary sector']
        energy_df.loc[(country, year), 'thermal uses services']        = clever_tertiary.loc[country, 'Thermal_uses_tertiary']

        energy_df.loc[(country, year), 'total agriculture']            = clever_agriculture.loc[country, 'Total Final energy consumption in agriculture']
        energy_df.loc[(country, year), 'total agriculture electricity']= clever_agriculture.loc[country, 'Final electricity consumption in agriculture']
        energy_df.loc[(country, year), 'total agriculture machinery']  = clever_agriculture.loc[country, 'Final oil consumption in agriculture']
        energy_df.loc[(country, year), 'total agriculture heat']       = clever_agriculture.loc[country, 'Total_agriculture_heat']

        # override population if needed (affects later transport step)
        # … etc …

    # 4) Override CO₂ with Clever_AFOLUB
    for country in co2_df.index:
        co2_df.at[country, 'agriculture'] = clever_afolub.loc[country, 'Total CO2 emissions from agriculture']
        co2_df.at[country, 'LULUCF']       = clever_afolub.loc[country, 'Total CO2 emissions from the LULUCF sector']

    # 5) Recompute transport “number cars” if needed
    if is_clever:
        for country in transport_df.index.get_level_values(0).unique():
            year0 = snakemake.params.transport_year if hasattr(snakemake.params, "transport_year") else year
            ppl_per_vehicle = clever_transport.loc[country, 'Average number of people per vehicle']
            stock = transport_df.loc[(country, year0), 'number cars'].sum()
            transport_df.loc[(country, year0), 'number cars'] = stock / ppl_per_vehicle

# 6) Write out the new files
energy_df.to_csv(snakemake.output.energy_C)
district_heat_df.to_csv(snakemake.output.district_heat_C)
co2_df.to_csv(snakemake.output.co2_C)
transport_df.to_csv(snakemake.output.transport_C)
