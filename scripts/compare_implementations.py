import subprocess
import matplotlib.pyplot as plt
import os
import re

# Define input files and corresponding sizes
input_files = ["small.in", "medium.in", "large.in"]
input_sizes = ["Small (700 x 700)", "Medium (1100 x 1100)", "Large (1500 x 1500)"]
input_dir = "inputs"
executables = ["sequential", "mpi", "openmp", "hybrid"]

# Run the program with each input file and measure the runtime
def run_program():
    runtimes = {exec_name: [] for exec_name in executables}  # Initialize a dictionary for all executables
    
    for input_file in input_files:
        input_path = os.path.join(input_dir, input_file)
        print(f"Running with input file: {input_path}")
        
        for exec_name in executables:
            print(f"Running {exec_name} with input file: {input_path}")
            # Execute the program and capture the output
            if exec_name == "mpi" or exec_name == "hybrid":
                result = subprocess.run(
                    ["mpirun", "-np", "4", f"./{exec_name}", input_path],
                    capture_output=True,
                    text=True
                )
            else:
                result = subprocess.run(
                    [f"./{exec_name}", input_path],
                    capture_output=True,
                    text=True
                )

            # Extract the runtime from the output using regex
            runtime_match = re.search(r"Total time: (\d+\.\d+)", result.stdout)
            if runtime_match:
                runtime = float(runtime_match.group(1))
                runtimes[exec_name].append(runtime)
                print(f"Runtime for {exec_name} with {input_file}: {runtime} seconds")
            else:
                print(f"Error: Could not extract runtime for {exec_name} with {input_file}")
                runtimes[exec_name].append(None)

    return runtimes

# Plot the runtimes using a bar chart
def plot_runtimes(runtimes, filename="comparison.png"):
    # Number of input files (Small, Medium, Large)
    n = len(input_sizes)
    
    # Set up the bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Width of the bars
    bar_width = 0.2
    
    # Define positions for each set of bars
    indices = range(n)
    offset = -bar_width

    # Plot the bars for each executable (sequential, mpi, openp)
    for i, exec_name in enumerate(executables):
        ax.bar(
            [index + offset + bar_width * i for index in indices],
            runtimes[exec_name],
            bar_width,
            label=exec_name.capitalize(),  # Capitalize for better display
        )
    
    # Set labels, title, and other plot elements
    ax.set_xlabel("Input Size")
    ax.set_ylabel("Runtime (seconds)")
    ax.set_title("Runtime Comparison of Different Executables for Various Input Sizes")
    ax.set_xticks([index + bar_width for index in indices])  # Set x-ticks in the center of the grouped bars
    ax.set_xticklabels(input_sizes)
    ax.legend(title="Executables")
    ax.grid(axis='y')
    plt.tight_layout()
    
    # Save the plot with a specified filename
    output_path = os.path.join("./photos", filename)
    plt.savefig(output_path)
    print(f"Plot saved as {output_path}")
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
