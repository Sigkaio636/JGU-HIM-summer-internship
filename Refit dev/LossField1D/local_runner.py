import time
import pathlib
import yaml
from fabric import Connection


# 1) Load configuration
cfg = yaml.safe_load(open("config.yaml"))
host = cfg["hpc"]["host"]
user = cfg["hpc"]["user"]
remote_base = cfg["hpc"]["remote_base"]
Vc_list = cfg["sweep"]["Vc"]

# 2) Generate a unique run_id
timestamp = time.strftime("%Y%m%dT%H%M%S")
run_id = f"run-{timestamp}"
remote_dir = f"{remote_base}/{run_id}"

# 3) SSH connection (assume that your key is in  ~/.ssh/id_rsa)
key_file = str(pathlib.Path.home() / ".ssh" / "id_ed25519")
conn = Connection(host=host, user=user, connect_kwargs={"key_filename": key_file})

# 4) Create remote directory
conn.run(f"mkdir -p {remote_dir}")

# 5) Upload scripts and configuration
conn.put("create_inpfile.py", remote=remote_dir + "/create_inpfile.py")
conn.put("config.yaml", remote=remote_dir + "/config.yaml")
conn.put("create_outcsv.py", remote=remote_dir + "/create_outcsv.py")

# 6) Build the remote call
cmd = f"cd {remote_dir} && " f"python3 create_inpfile.py "

# 7) Execute and log
print("▶ Launching input generation on the cluster...")
result = conn.run(cmd, pty=True)
print("✔ stdout:")
print(result.stdout)
print("✖ stderr:")
print(result.stderr)

# 8) Obtain the input list remotly
result = conn.run(f"ls {remote_dir}/inputs/*_input.dat", hide=True)
input_files = result.stdout.strip().splitlines()

# 9) Create the output directory 
conn.run(f"mkdir -p {remote_dir}/outputs")

# 10) Execute program for each input
for input_path in input_files:
    # Obtain base name (input_001.dat → 001)
    base_name = pathlib.Path(input_path).stem.split("_")[0]
    # print(input_path, base_name)

    # Execute the binary with I/O redirection
    print(f"▶ Lanzando computo {base_name} por Boscos...")
    result = conn.run(
        f"{remote_base}/Boscos.out < {input_path} > {remote_dir}/outputs/{base_name}_output.txt"
    )
    # print("✔ stdout:")
    # print(result.stdout)
    # print("✖ stderr:")
    # print(result.stderr)

# 11) Execute output reading and csv generation
print("▶ Launching results reading and CSV generation...")
conn.run(f"cd {remote_dir} && " f"python3 create_outcsv.py")
print("✔ stdout:")
print(result.stdout)
print("✖ stderr:")
print(result.stderr)

# 12) Obtain csv with final results
local_path = f"C:/Users/ekici/Desktop/Mainz HIM/Refit dev/LossField1D/res_{run_id}.csv"
conn.get(remote=f"{remote_dir}/resultados.csv", local=local_path)

# x) Close connection
conn.close()
print(":3")
