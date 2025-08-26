import pathlib
import yaml
import csv
import re

script_dir = pathlib.Path(__file__).parent.resolve()
output_dir = script_dir / "outputs"
output_files = list(output_dir.iterdir())
print("output_files", output_files)

# Lee el yaml
with open(f"{script_dir}/config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Extrae Vc_values
Vc_param = config.get("sweep", {}).get("Vc", [-80, 2, -70])  # default por si no está
Vc_values = []
while Vc_param[0] <= Vc_param[2]:
    Vc_values.append(round(Vc_param[0], 10))
    Vc_param[0] += Vc_param[1]

results_save = []
for i, output_file in enumerate(output_files):
    # Obtener nombre base (input_001.dat → 001)
    output_path = output_dir / output_file
    base_name = pathlib.Path(output_path).stem.split("_")[0]

    with output_path.open("r") as f:
        content = f.read()
        print(content)

        # Usa una expresión regular para extraer un valor numérico
        match = re.search(r"E=\s*([-+]?\d+\.\d+)\s+MeV", content)
        if match:
            energy = float(match.group(1))
        else:
            energy = None  # o puedes saltarlo con `continue`

        results_save.append({"archivo": base_name, "energia": energy})

results_sorted = sorted(results_save, key=lambda x: x["archivo"])
for result, vc in zip(results_sorted, Vc_values):
    result["Vc"] = vc
print("results_sorted", results_sorted)

# Guardar CSV
csv_path = pathlib.Path(f"{script_dir}/resultados.csv")
with open(csv_path, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["archivo", "Vc", "energia"])
    writer.writeheader()
    writer.writerows(results_sorted)
