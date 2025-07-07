#!/usr/bin/env python3
import os
import re
import logging
import pandas as pd
import matplotlib.pyplot as plt
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake
        snakemake = mock_snakemake("evaluate_sensitivity_analysis_comparison")

    configure_logging(snakemake)

    # Szenarien, die wir vergleichen wollen
    wanted = [
        "reference",
        "ref-constr01",
        "sens1.5-reference",
        "sens1.5-ref-constr01",
        "sens2.0-reference",
        "sens2.0-ref-constr01",
    ]

    # CSVs einlesen
    df_list = []
    for path in snakemake.input.space_requirements_DLU_csv:
        run = os.path.basename(os.path.dirname(os.path.dirname(path)))
        if run not in wanted:
            continue
        df = pd.read_csv(path)
        df["run"] = run
        df_list.append(df)
    if not df_list:
        logger.error("Keine Daten für die gewünschten Szenarien gefunden!")
        raise SystemExit(1)
    df_all = pd.concat(df_list, ignore_index=True)

    # Jahre-Spalten ins Long-Format und Umrechnung in km²
    year_cols = [c for c in df_all.columns if re.fullmatch(r"\d{4}", c)]
    df_m = (
        df_all
        .melt(
            id_vars=["run", "country", "carrier"],
            value_vars=year_cols,
            var_name="year",
            value_name="area_m2"
        )
        .assign(
            year=lambda d: d["year"].astype(int),
            area_km2=lambda d: d["area_m2"] * 1e-6
        )
    )

    # Nur die drei interessierenden Technologien
    techs = ["unsustainable solid biomass", "unsustainable bioliquids", "unsustainable biogas", "solar", "onwind"]
    df_m = df_m[df_m["carrier"].isin(techs)].copy()

    # Alle Jahre einmal sortiert
    years = sorted(df_m["year"].unique())

    # Feste Farben (gleich für dashed+solid) und Linienstile
    color_map = {
        "reference":           "#1b5e20",  # dunkelgrün
        "ref-constr01":        "#1b5e20",
        "sens1.5-reference":   "#f9a825",  # goldgelb
        "sens1.5-ref-constr01":"#f9a825",
        "sens2.0-reference":   "#b71c1c",  # dunkelrot
        "sens2.0-ref-constr01":"#b71c1c",
    }
    linestyle_map = {
        "reference":            "--",
        "sens1.5-reference":    "--",
        "sens2.0-reference":    "--",
        "ref-constr01":         "-",
        "sens1.5-ref-constr01": "-",
        "sens2.0-ref-constr01": "-",
    }

    # Plot-Defaults wie im Beispiel
    plt.rcParams.update({
        "font.size": 7,
        "axes.titlesize": 7,
        "lines.linewidth": 0.5,
    })
    width = 5.46
    height = width * (10/30) * len(techs)
    fig, axes = plt.subplots(len(techs), 2, figsize=(width, height), sharex=True)

    regions = [
        ("All Countries",      lambda df: df),
        ("Germany",            lambda df: df[df["country"] == "DE"])
    ]

    # Pro Technologie (Reihe) und Region (Spalte) plotten
    for i, tech in enumerate(techs):
        for j, (region_label, region_filter) in enumerate(regions):
            ax = axes[i, j]
            sub = region_filter(df_m[df_m["carrier"] == tech])
            for run, grp in sub.groupby("run"):
                yearly = (
                    grp.groupby("year")["area_km2"]
                       .sum()
                       .reindex(years, fill_value=0)
                       .reset_index()
                )
                ax.plot(
                    yearly["year"],
                    yearly["area_km2"],
                    label=run,
                    color=color_map[run],
                    linestyle=linestyle_map[run],
                )
            # Titel und Achsenbeschriftung
            title = tech[0].upper() + tech[1:] + "\n" + region_label
            ax.set_title(title)
            ax.set_ylabel("DLU in km²")
            ax.grid(True, linestyle="--", linewidth=0.3)
            ax.set_xlim(years[0], years[-1])
            if i == len(techs)-1:
                ax.set_xticks(years)
                ax.set_xlabel("Year")

    # Gemeinsame Legende unten
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0],
               color=color_map[r],
               linestyle=linestyle_map[r],
               label=r)
        for r in wanted
    ]
    fig.legend(
        handles=handles,
        title="Scenarios",
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.15),
        frameon=False
    )

    fig.subplots_adjust(
        bottom=-0.05,
        left=0.08,
        right=0.98,
        top=0.95,
        hspace=0.4,
        wspace=0.4
    )

    # Speichern
    fig.savefig(
        snakemake.output.sensitivity_analysis_plot,
        dpi=300,
        bbox_inches="tight"
    )
    logger.info(f"Saved sensitivity-analysis comparison plot to {snakemake.output.sensitivity_analysis_plot}")

    # ─────────────────────────────────────────────────────────────────────────────
    # Zweite Abbildung: Gesamte DLU über alle Technologien
    # ─────────────────────────────────────────────────────────────────────────────
    fig2, axes2 = plt.subplots(1, 2, figsize=(width, width * (10 / 14)), sharex=True)

    for j, (region_label, region_filter) in enumerate(regions):
        ax = axes2[j]
        sub_all = region_filter(df_m)  # alle Technologieträger zusammen
        for run, grp in sub_all.groupby("run"):
            total = (
                grp
                .groupby("year")["area_km2"]
                .sum()
                .reindex(years, fill_value=0)
                .reset_index()
            )
            ax.plot(
                total["year"],
                total["area_km2"]*1e-3,
                label=run,
                color=color_map[run],
                linestyle=linestyle_map[run],
            )
        ax.set_title(f"Total DLU\n{region_label}")
        ax.set_ylabel("DLU in 1000 km²")
        ax.grid(True, linestyle="--", linewidth=0.3)
        ax.set_xlim(years[0], years[-1])
        if j == 1:
            ax.set_xticks(years)
            ax.set_xlabel("Year")

    # Gemeinsame Legende unter der zweiten Abbildung
    from matplotlib.lines import Line2D

    handles2 = [
        Line2D([0], [0], color=color_map[r], linestyle=linestyle_map[r], label=r)
        for r in wanted
    ]
    fig2.legend(
        handles=handles2,
        title="Scenario",
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.2),
        frameon=False
    )

    fig2.subplots_adjust(
        bottom=0.02,
        left=0.08,
        right=0.98,
        top=0.90,
        wspace=0.3
    )

    fig2.savefig(
        snakemake.output.sensitivity_analysis_plot.replace(".png", "_total.png"),
        dpi=300,
        bbox_inches="tight"
    )
    logger.info(
        f"Saved total-DLU comparison plot to {snakemake.output.sensitivity_analysis_plot.replace('.png', '_total.png')}")

