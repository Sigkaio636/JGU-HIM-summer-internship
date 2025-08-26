import time
import pathlib
import yaml
from fabric import Connection
from invoke.exceptions import UnexpectedExit

# 1) Load configuration
cfg = yaml.safe_load(open("config.yaml"))
host = cfg["hpc"]["host"]
user = cfg["hpc"]["user"]
remote_base = cfg["hpc"]["remote_base"]

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
conn.put("remote_runner.py", remote=remote_dir + "/remote_runner.py")
conn.put("config.yaml", remote=remote_dir + "/config.yaml")

# 6) Create output directory
conn.run(f"mkdir -p {remote_dir}/outputs")

# 7) Initialize remote executer
print("▶ Lunching remote Runner  ...")
result = conn.run(f"cd {remote_dir} && " f"python3 remote_runner.py")
# print("✔ stdout:")
# print(result.stdout)
# print("✖ stderr:")
# print(result.stderr)

# 8) Obtain the input list remotly
local_path = f"C:/Users/ekici/Desktop/Mainz HIM/Refit dev/VcVso_Opt/res_{run_id}.csv"

try:
    conn.run(f"test -f {remote_dir}/run_log.csv", hide=True)
    file_exists = True
except UnexpectedExit:
    file_exists = False

if file_exists:
    conn.get(remote=f"{remote_dir}/run_log.csv", local=local_path)
else:
    print("run_log.csv not generated")

# x) Close connection
conn.close()
print(":3")
