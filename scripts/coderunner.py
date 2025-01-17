import subprocess
import re
import csv

def get_total_time(command):
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[Error] Command failed: {command}")
            print(result.stderr)
            return None

        search = re.search(r"Total time:\s+(\S+)", result.stdout)
        if search:
            return float(search.group(1))
        return None
    except Exception as e:
        print(f"[Exception] An error occurred while executing {command}: {str(e)}")
        return None

configurations = [
    {"threads": 1, "cores": 1, "command": "./sequential inputs/large.in"},
    {"threads": 2, "cores": 1, "command": "./openmp inputs/large.in 2"},
    {"threads": 4, "cores": 1, "command": "./openmp inputs/large.in 4"},
    {"threads": 8, "cores": 1, "command": "./openmp inputs/large.in 8"},
    {"threads": 1, "cores": 2, "command": "mpirun -np 2 ./mpi inputs/large.in"},
    {"threads": 1, "cores": 4, "command": "mpirun -np 4 ./mpi inputs/large.in"},
    {"threads": 1, "cores": 8, "command": "mpirun -np 8 --oversubscribe ./mpi inputs/large.in"},
    {"threads": 2, "cores": 2, "command": "mpirun -np 2 ./hybrid inputs/large.in 2"},
    {"threads": 4, "cores": 2, "command": "mpirun -np 2 ./hybrid inputs/large.in 4"},
    {"threads": 8, "cores": 2, "command": "mpirun -np 2 ./hybrid inputs/large.in 8"},
    {"threads": 2, "cores": 4, "command": "mpirun -np 4 ./hybrid inputs/large.in 2"},
    {"threads": 4, "cores": 4, "command": "mpirun -np 4 ./hybrid inputs/large.in 4"},
    {"threads": 8, "cores": 4, "command": "mpirun -np 4 ./hybrid inputs/large.in 8"},
    {"threads": 2, "cores": 8, "command": "mpirun -np 8 --oversubscribe ./hybrid inputs/large.in 2"},
    {"threads": 4, "cores": 8, "command": "mpirun -np 8 --oversubscribe ./hybrid inputs/large.in 4"},
    {"threads": 8, "cores": 8, "command": "mpirun -np 8 --oversubscribe ./hybrid inputs/large.in 8"}
]

# Prepare to save the results
with open('performance_data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Threads", "Cores", "TotalTime"])

    for cfg in configurations:
        time = get_total_time(cfg["command"])
        if time is not None:
            writer.writerow([cfg["threads"], cfg["cores"], time])

print("Data collection complete and saved to performance_data.csv.")
