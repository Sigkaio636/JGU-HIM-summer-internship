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
Vc_param = config.get("sweep", {}).get("Vc", [-80, 5, -70])
Vc_values = []
while Vc_param[0] <= Vc_param[2]:
    Vc_values.append(round(Vc_param[0], 10))
    Vc_param[0] += Vc_param[1]
# Vc_values =equiv np.arange(start=-80, stop=-70, step=5)


# Extract Vso_values
Vso_param = config.get("sweep", {}).get("Vs", [-10, 5, 10])  # default por si no está
Vso_values = []
while Vso_param[0] <= Vso_param[2]:
    Vso_values.append(round(Vso_param[0], 10))
    Vso_param[0] += Vso_param[1]
# Vc_values =equiv np.arange(start=-10, stop=10, step=5)

# For each output_file read the binding energy 
results_save = []
for i, output_file in enumerate(output_files):
    # Obtain base name (input_001.dat → 001)
    output_path = output_dir / output_file
    base_name = pathlib.Path(output_path).stem.split("_")[0]

    with output_path.open("r") as f:
        content = f.read()

        # Use regular expression to extract binding energy
        match = re.search(r"E=\s*([-+]?\d+\.\d+)\s+MeV", content)
        if match:
            energy = float(match.group(1))
        else:
            energy = None  # o puedes saltarlo con `continue`

        results_save.append({"archivo": base_name, "energia": energy})

# Order results and add Vc,Vso to data tuple: ("archivo", "energia", "Vc", "Vso")
results_sorted = sorted(results_save, key=lambda x: x["archivo"])
i = 0
for vc in Vc_values:
    for vso in Vso_values:
        results_sorted[i]["Vc"] = vc
        results_sorted[i]["Vso"] = vso
        i += 1

# Save CSV
csv_path = pathlib.Path(f"{script_dir}/resultados.csv")
with open(csv_path, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["archivo", "Vc", "Vso", "energia"])
    writer.writeheader()
    writer.writerows(results_sorted)
