import os
import re
import pypsa
import matplotlib.pyplot as plt
from scripts._evaluation_helpers import load_networks_from_path_list, compare_value

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "general_scenario_comparison",
        )

    nn = load_networks_from_path_list(snakemake.input.networks)

    # Datenvergleich der Zielfunktion
    expr = lambda n: n.objective  # todo: richtige objective
    # DataFrame mit Werten in Mrd. Euro
    df_objective = compare_value(expr, nn) * 1e-9

    # Plot erstellen
    ax = df_objective.plot.line()
    plt.ylabel("Zielfunktion (Mrd. Euro)")
    fig = ax.get_figure()

    # Ausgabedatei speichern
    fig.savefig(snakemake.output.objective_graph)


    print("fertig")