import subprocess
import matplotlib.pyplot as plt
import os
import re

# Define input files and corresponding sizes
input_files = ["small.in", "medium.in", "large.in"]
input_sizes = ["Small (700 x 700)", "Medium (1100 x 1100)", "Large (1500 x 1500)"]
input_dir = "inputs"
target = "sequential"

# Run the C++ program with each input file and measure the runtime
def run_program():
    runtimes = []
    for input_file in input_files:
        input_path = os.path.join(input_dir, input_file)
        print(f"Running with input file: {input_path}")
        
        # Execute the program and capture the output
        result = subprocess.run(
            [f"./{target}", input_path],
            capture_output=True,
            text=True
        )

        # Extract the runtime from the output using regex
        runtime_match = re.search(r"Total time: (\d+\.\d+)", result.stdout)
        if runtime_match:
            runtime = float(runtime_match.group(1))
            runtimes.append(runtime)
            print(f"Runtime for {input_file}: {runtime} seconds")
        else:
            print(f"Error: Could not extract runtime for {input_file}")
            runtimes.append(None)

    return runtimes

# Plot the runtimes using a bar chart
def plot_runtimes(runtimes, filename="sequential_run_times.png"):
    plt.figure(figsize=(10, 6))
    plt.bar(input_sizes, runtimes, color='skyblue')
    plt.xlabel("Input Size")
    plt.ylabel("Runtime (seconds)")
    plt.title("Runtime of C++ Program for Different Input Sizes")
    plt.grid(axis='y')
    plt.tight_layout()
    
    # Save the plot with a specified filename
    output_path = os.path.join("./photos", filename)
    plt.savefig(output_path)
    print(f"Plot saved as {output_path}")
    plt.show()

if __name__ == "__main__":
    try:
        runtimes = run_program()
        if all(r is not None for r in runtimes):
            plot_runtimes(runtimes)
        else:
            print("Some runtimes could not be measured. Skipping plot.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred during compilation or execution: {e}")
