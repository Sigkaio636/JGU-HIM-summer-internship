import pathlib
import yaml
import csv
import re

# Get directory names
script_dir = pathlib.Path(__file__).parent.resolve()
output_dir = script_dir / "outputs"
output_files = list(output_dir.iterdir())
print("output_files", output_files)

# Read yaml
with open(f"{script_dir}/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Extract Vc_values
Vc_param = config.get("sweep", {}).get("Vc", [-80, 2, -70])
Vc_values = []
while Vc_param[0] <= Vc_param[2]:
    Vc_values.append(round(Vc_param[0], 10))
    Vc_param[0] += Vc_param[1]
# Vc_values =equiv np.arange(start=-80, stop=-70, step=2)

# For each output_file read the binding energy 
results_save = []
for i, output_file in enumerate(output_files):
    # Obtain base name (input_001.dat → 001)
    output_path = output_dir / output_file
    base_name = pathlib.Path(output_path).stem.split("_")[0]

    with output_path.open("r") as f:
        content = f.read()
        print(content)

        # Use regular expression to extract binding energy
        match = re.search(r"E=\s*([-+]?\d+\.\d+)\s+MeV", content)
        if match:
            energy = float(match.group(1))
        else:
            energy = None

        results_save.append({"archivo": base_name, "energia": energy})

# Order results and add Vc to data tuple: ("archivo", "energia", "Vc")
results_sorted = sorted(results_save, key=lambda x: x["archivo"])
for result, vc in zip(results_sorted, Vc_values):
    result["Vc"] = vc
print("results_sorted", results_sorted)

# Save CSV
csv_path = pathlib.Path(f"{script_dir}/resultados.csv")
with open(csv_path, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["archivo", "Vc", "energia"])
    writer.writeheader()
    writer.writerows(results_sorted)
