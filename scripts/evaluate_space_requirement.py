import os

# Snakemake input is provided as a list
networks = snakemake.input.networks

# Print the loaded networks
print("Loaded networks:")
for network in networks:
    print(f"- {network}")

# Ensure the output directory exists
output_dir = snakemake.output.report
os.makedirs(output_dir, exist_ok=True)

# Write results to a file
output_file = os.path.join(output_dir, "space_requirement_report.txt")
with open(output_file, "w") as f:
    f.write("Loaded networks test:\n")
    for network in networks:
        f.write(f"- {network}\n")