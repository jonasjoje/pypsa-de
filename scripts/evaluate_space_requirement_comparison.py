import os
import re
import logging
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scripts._helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake
        snakemake = mock_snakemake("evaluate_space_requirement_comparison")

    configure_logging(snakemake)

    # ─────────────────────────────────────────────────────────────────────────────
    # Load data
    # ─────────────────────────────────────────────────────────────────────────────
    df_list = []
    for path in snakemake.input.space_requirements_DLU_csv:
        run_name = os.path.basename(os.path.dirname(os.path.dirname(path)))
        df = pd.read_csv(path)
        df["run"] = run_name
        df_list.append(df)

    df_all = pd.concat(df_list, ignore_index=True)

    # Filter only rows with area > 0 in 2050
    df_all = df_all[df_all["2050"] > 0].copy()

    # Rename "2050" to a consistent column name for plotting
    df_all = df_all.rename(columns={"2050": "area_km2"})
    # area_km2 war in m² – erst in km², dann in 1000 km² umrechnen:
    df_all['area_km2'] = df_all['area_km2'] * 1e-6  # m² → km²


    # Extract constraint value from run name
    def extract_constraint(run):
        match = re.search(r"constr(\d+)", run)
        if match:
            return int(match.group(1))
        return 100


    df_all["constraint"] = df_all["run"].apply(extract_constraint)

    # Filter carriers with nonzero total area
    carriers_with_area = (
        df_all.groupby("carrier")["area_km2"].sum().loc[lambda s: s > 0].index.tolist()
    )
    df_all = df_all[df_all["carrier"].isin(carriers_with_area)]


    # ─────────────────────────────────────────────────────────────────────────────
    # Subsetting for plot panels
    # ─────────────────────────────────────────────────────────────────────────────
    def subset(df, prefix, countries):
        sub = df[df["run"].str.startswith(prefix)]
        if countries is not None:
            sub = sub[sub["country"].isin(countries)]
        grouped = sub.groupby(["constraint", "carrier"])["area_km2"].sum().unstack().fillna(0) / 1e3
        return grouped


    ref_all = subset(df_all, "ref", None)
    cle_all = subset(df_all, "cle", None)
    ref_de = subset(df_all, "ref", ["DE"])
    cle_de = subset(df_all, "cle", ["DE"])

    # Determine carrier plotting order by total area
    all_data = pd.concat([ref_all, cle_all, ref_de, cle_de])
    carrier_order = all_data.sum().sort_values(ascending=False).index.tolist()

    # ─────────────────────────────────────────────────────────────────────────────
    # Plotting
    # ─────────────────────────────────────────────────────────────────────────────
    # Globale Plot-Einstellungen
    plt.rcParams.update({
        'font.size': 7,
        'axes.titlesize': 7,
        'lines.linewidth': 0.5,  # Linien nur halb so dick
    })

    new_width = 5.46
    new_height = new_width * 10 / 14  # erhaltenes Seitenverhältnis: ~3.9
    fig, axs = plt.subplots(2, 2, figsize=(new_width, new_height), sharex='col')

    plot_data = [
        (ref_all, "Reference – All Countries", axs[0, 0]),
        (cle_all, "Clever – All Countries", axs[0, 1]),
        (ref_de, "Reference – Germany", axs[1, 0]),
        (cle_de, "Clever – Germany", axs[1, 1]),
    ]

    color_map = {}

    for data, title, ax in plot_data:
        data = data.reindex(columns=carrier_order, fill_value=0)
        x = data.index
        y = data.values.T
        polys = ax.stackplot(x, y, labels=data.columns)

        if not color_map:
            color_map = {name: poly.get_facecolor()[0] for name, poly in zip(data.columns, polys)}
        else:
            for poly, name in zip(polys, data.columns):
                poly.set_color(color_map[name])

        ax.set_title(title)
        ax.set_xlabel("Constraint (%)")
        ax.set_ylabel("Land Area in 1000 km²")
        ax.ticklabel_format(style='plain', axis='y')  # keine wissenschaftliche Notation

        ax.set_xlim(105, 0)
        ax.set_xticks(data.index)

    # Legende unter den Plots
    handles = [mpatches.Patch(label=k, color=color_map[k]) for k in reversed(carrier_order)]
    fig.legend(
        handles=handles,
        title="Carrier",
        loc='lower center',
        ncol=5,
        bbox_to_anchor=(0.5, -0.18),
        frameon=False
    )

    # Mehr Platz für Legende und X-Achsentitel
    fig.subplots_adjust(
        bottom=0.1,  # schiebt die gesamte Grafik nach oben
        left=0.08,
        right=0.98,
        top=0.95,
        hspace=0.3,  # vertikaler Abstand zwischen den Subplots
        wspace=0.25  # horizontaler Abstand
    )

    # Speichern mit tight bbox, damit auch die Legende außen mitgedruckt wird
    fig.savefig(
        snakemake.output.DLU_vs_constraint_stack,
        dpi=300,
        bbox_inches='tight'
    )
    logger.info(f"Saved plot to {snakemake.output.DLU_vs_constraint_stack}")



    # Globale Einstellungen
    plt.rcParams['font.size'] = 7
    plt.rcParams['axes.titlesize'] = 7
    plt.rcParams['lines.linewidth'] = 0.5  # Linien nur halb so dick

    # ─────────────────────────────────────────────────────────────────────────────
    # Plot mu over constraint
    # ─────────────────────────────────────────────────────────────────────────────
    mu_list = []
    for path in snakemake.input.constraint_mu_DLU_csv:
        run_name = os.path.basename(os.path.dirname(os.path.dirname(path)))
        df = pd.read_csv(path)
        df["run"] = run_name
        mu_list.append(df)

    df_mu = pd.concat(mu_list, ignore_index=True)

    # Filter only TotalSpaceRequirement_DLU constraints
    df_mu = df_mu[df_mu["index"].str.contains("TotalSpaceRequirement_DLU")].copy()

    # Extract constraint, country, scenario
    df_mu["constraint"] = df_mu["run"].apply(extract_constraint)
    df_mu["scenario"] = df_mu["run"].str.extract(r"^(ref|cle)")[0]
    df_mu["country"] = df_mu["index"].str.extract(r"DLU_([A-Z]{2})")[0].fillna("ALL")

    # Keep only year 2050 and drop NaN mu
    df_mu = df_mu[["scenario", "constraint", "country", "2050"]].rename(columns={"2050": "mu"})
    df_mu = df_mu.dropna(subset=["mu"])
    df_mu["mu"] = -df_mu["mu"]

    df_mu = df_mu[df_mu["country"] != "CH"]  # Ausreißer rausnehmen

    # Prepare consistent color mapping
    countries = sorted(df_mu["country"].unique())
    color_cycle = plt.cm.get_cmap("tab20", len(countries))
    color_map = {country: color_cycle(i) for i, country in enumerate(countries)}

    # Plotting
    # Breite 5.46 inch, Höhe beibehalten
    old_width, old_height = 14, 6
    new_width = 5.46
    new_height = new_width * old_height / old_width

    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(new_width, new_height),
        sharey=True
    )

    for ax, scen in zip((ax1, ax2), ("ref", "cle")):
        sub = df_mu[df_mu["scenario"] == scen]
        for country, group in sub.groupby("country"):
            ax.plot(group["constraint"], group["mu"], label=country, color=color_map[country])
        ax.set_title(f"{'Reference' if scen == 'ref' else 'Clever'} Scenarios")
        ax.set_xlabel("Constraint (%)")
        ax.set_xlim(55, 0)
        ax.set_xticks([50, 25, 1])
        ax.set_yscale('log')
        #ax.set_ylim(0, 2)
        ax.grid(True, linestyle='--', linewidth=0.3)

    ax1.set_ylabel("Marginal Value of Land (2050)\nin million EUR/km²")

    # Legende unterhalb, 5 Spalten, außerhalb des Graphen
    handles = [plt.Line2D([0], [0], color=color_map[c], label=c) for c in countries]
    fig.legend(
        handles=handles,
        title="Country",
        loc='lower center',
        ncol=7,
        bbox_to_anchor=(0.5, -0.25),
        frameon=False
    )

    # Mehr Platz unten lassen
    fig.tight_layout(rect=[0, 0.0, 1, 1])

    # Speichern in 300 dpi
    fig.savefig(
        snakemake.output.DLU_mu_vs_constraint,
        dpi=300,
        bbox_inches='tight'
    )

    logger.info(f"Saved mu plot to {snakemake.output.DLU_mu_vs_constraint}")

    # ─────────────────────────────────────────────────────────────────────────────
    # Disturbed
    # ─────────────────────────────────────────────────────────────────────────────
    df_list = []
    for path in snakemake.input.space_requirements_dist_csv:
        run_name = os.path.basename(os.path.dirname(os.path.dirname(path)))
        df = pd.read_csv(path)
        df["run"] = run_name
        df_list.append(df)

    df_all = pd.concat(df_list, ignore_index=True)

    # Filter only rows with area > 0 in 2050
    df_all = df_all[df_all["2050"] > 0].copy()

    # Rename "2050" to a consistent column name for plotting
    df_all = df_all.rename(columns={"2050": "area_km2"})
    # area_km2 war in m² – in km² umrechnen
    df_all["area_km2"] = df_all["area_km2"] * 1e-6


    # Extract constraint value from run name
    def extract_constraint(run):
        match = re.search(r"constr(\d+)", run)
        if match:
            return int(match.group(1))
        return 100


    df_all["constraint"] = df_all["run"].apply(extract_constraint)

    # Filter carriers with nonzero total area
    carriers_with_area = (
        df_all.groupby("carrier")["area_km2"].sum().loc[lambda s: s > 0].index.tolist()
    )
    df_all = df_all[df_all["carrier"].isin(carriers_with_area)]


    # Subsetting for plot panels
    def subset(df, prefix, countries):
        sub = df[df["run"].str.startswith(prefix)]
        if countries is not None:
            sub = sub[sub["country"].isin(countries)]
        grouped = sub.groupby(["constraint", "carrier"])["area_km2"].sum().unstack().fillna(0) / 1e3
        return grouped


    ref_all = subset(df_all, "ref", None)
    cle_all = subset(df_all, "cle", None)
    ref_de = subset(df_all, "ref", ["DE"])
    cle_de = subset(df_all, "cle", ["DE"])

    # Determine carrier plotting order by total area
    all_data = pd.concat([ref_all, cle_all, ref_de, cle_de])
    carrier_order = all_data.sum().sort_values(ascending=False).index.tolist()


    # Plotting
    # Plotting im gleichen Stil wie der mu-Plot
    new_width = 5.46
    new_height = new_width * 10 / 14
    fig, axs = plt.subplots(2, 2, figsize=(new_width, new_height), sharex='col')

    plot_data = [
        (ref_all, "Reference – All Countries", axs[0, 0]),
        (cle_all, "Clever – All Countries", axs[0, 1]),
        (ref_de, "Reference – Germany", axs[1, 0]),
        (cle_de, "Clever – Germany", axs[1, 1]),
    ]

    color_map = {}

    for data, title, ax in plot_data:
        data = data.reindex(columns=carrier_order, fill_value=0)
        x = data.index
        y = data.values.T
        polys = ax.stackplot(x, y, labels=data.columns)

        if not color_map:
            color_map = {
                name: poly.get_facecolor()[0]
                for name, poly in zip(data.columns, polys)
            }
        else:
            for poly, name in zip(polys, data.columns):
                poly.set_color(color_map[name])

        ax.set_title(title)
        ax.set_xlabel("Constraint (%)")
        ax.set_ylabel("Disturbed Area in 1000 km²")
        ax.set_xlim(105, 0)
        ax.set_xticks(data.index)
        ax.ticklabel_format(style='plain', axis='y')  # keine wissenschaftliche Notation

    # Legende unter den Plots
    handles = [
        mpatches.Patch(label=k, color=color_map[k])
        for k in reversed(carrier_order)
    ]
    fig.legend(
        handles=handles,
        title="Carrier",
        loc='lower center',
        ncol=5,
        bbox_to_anchor=(0.5, -0.18),
        frameon=False
    )

    # Mehr Platz für Legende und X-Achsentitel
    fig.subplots_adjust(
        bottom=0.1,  # schiebt die Grafik nach oben
        left=0.08,
        right=0.98,
        top=0.95,
        hspace=0.3,
        wspace=0.25
    )

    # Speichern mit tight bbox
    fig.savefig(
        snakemake.output.dist_vs_constraint_stack,
        dpi=300,
        bbox_inches='tight'
    )
    logger.info(f"Saved plot to {snakemake.output.dist_vs_constraint_stack}")

