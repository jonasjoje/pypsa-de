import os
import pypsa
import pandas as pd

# Snakemake input is provided as a list
networks = snakemake.input.networks
planning_horizons = snakemake.params.planning_horizons

# Print the loaded networks
print("Loaded networks:")
for network in networks:
    print(f"- {network}")

# Dictionary to store space requirements per year
space_requirements = {}

# Process each network file
for network_file, year in zip(networks, planning_horizons):
    n = pypsa.Network(network_file)
    total_space = n.generators["space_req_opt"].sum() * 1e-3  # Convert from 1000 m² to km²
    space_requirements[year] = total_space

# Ensure the output directory exists
output_dir = snakemake.output.report
os.makedirs(output_dir, exist_ok=True)

# Write results to a file
output_file = os.path.join(output_dir, "space_requirement_report.txt")
with open(output_file, "w") as f:
    f.write("Loaded networks:\n")
    for network in networks:
        f.write(f"- {network}\n")
    f.write("\nSpace Requirements:\n")
    for year, space in space_requirements.items():
        f.write(f"{year}: {space:.2f} km²\n")