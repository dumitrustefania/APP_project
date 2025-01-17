import subprocess
import matplotlib.pyplot as plt
import os
import re

# Define input files and corresponding sizes
input_path = "inputs/small.in"

# Run the program with each input file and measure the runtime
def run_program():
    runtimes = {}  # Initialize a dictionary for all executables
    
    #  for i in 1,2 4, 8, i*2
    for th_num in (1, 2, 4, 8, 16):
        runtimes[th_num] = {}
        for proc_num in  (1, 2, 4, 8, 16):
            print(f"Running with {th_num} threads and {proc_num} processes")
            if proc_num >= 8:
                result = subprocess.run(
                    ["mpirun", "-np", str(proc_num), "--oversubscribe", "./hybrid", input_path, str(th_num)],
                    capture_output=True,
                    text=True
                )
            else:
                result = subprocess.run(
                        ["mpirun", "-np", str(proc_num), "./hybrid", input_path, str(th_num)],
                        capture_output=True,
                        text=True
                    )
          
            runtime_match = re.search(r"Total time: (\d+\.\d+)", result.stdout)
            if runtime_match:
                runtime = float(runtime_match.group(1))
                runtimes[th_num][proc_num] = runtime
                print(f"Runtime for {th_num} threads and {proc_num} processes: {runtime} seconds")
            else:
                print(f"Error: Could not extract runtime for {th_num} threads and {proc_num} processes")
                runtimes[th_num][proc_num] = None
    return runtimes

def plot_runtimes(runtimes):
    # Prepare data for plotting
    threads, processes, runtime_values = [], [], []

    for thread_num, proc_map in runtimes.items():
        for proc_num, runtime in proc_map.items():
            threads.append(thread_num)
            processes.append(proc_num)
            runtime_values.append(runtime)
    
    # Normalize runtime values for color mapping
    norm = plt.Normalize(min(runtime_values), max(runtime_values))
    colors = plt.cm.coolwarm(norm(runtime_values))  # Coolwarm colormap for visualization

    # Create 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(threads, processes, runtime_values, c=colors, marker='o', s=50)

    # Add colorbar
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='coolwarm'), ax=ax)
    cbar.set_label('Runtime (seconds)')

    # Set axis labels and ticks
    ax.set_xlabel('Threads')
    ax.set_ylabel('Processes')
    ax.set_zlabel('Runtime (seconds)')

    # Set ticks for threads and processes (Ox and Oy)
    tick_values = [1, 2, 4, 8, 16]
    ax.set_xticks(tick_values)
    ax.set_yticks(tick_values)

    # Title and show plot
    plt.title('Runtime vs Threads and Processes')
    plt.show()

if __name__ == "__main__":
    try:
        runtimes = run_program()
        # Ensure that all runtime data is valid (not None)
        if all(all(r is not None for r in exec_runtimes) for exec_runtimes in runtimes.values()):
            plot_runtimes(runtimes)
        else:
            print("Some runtimes could not be measured. Skipping plot.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred during compilation or execution: {e}")
