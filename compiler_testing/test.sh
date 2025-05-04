#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u
# Pipestatus: exit code of a pipeline is the status of the last command to exit with non-zero status, or zero if no command exited with a non-zero status
set -o pipefail

# --- Configuration ---

# Directory containing the C++ source files (relative to this script's location)
CPP_SRC_DIR="cpp_impl"
# Source files
SOURCE_FILES=("$CPP_SRC_DIR/main.cpp" "$CPP_SRC_DIR/comp.cpp" "$CPP_SRC_DIR/perf.cpp")
# Output executable name (will be placed in CPP_SRC_DIR)
OUTPUT_EXEC="$CPP_SRC_DIR/SparseGEMM.out"
# Output file for results (will be placed in the same directory as the script)
RESULTS_FILE="compiler_testing/compiler_results.txt"
# Arguments to pass to the executable
RUN_ARGS="-M 32 -K 1024 -N 4096 -s 4"
# Required preprocessor define
DEFINE_FLAGS="-DPMU"

# --- Compiler Configurations ---
# List of commands and titles. Adjust g++-14 if your version differs.
# Format: "Title" "Compiler Command Base"
# The script will add source files, -o target, and defines automatically.
declare -a configs=(
    "Clang O0"              "clang++ -O0 -fopenmp"
    "Clang O2 Native"       "clang++ -O2 -march=native -fopenmp"
    "Clang O3 Native"       "clang++ -O3 -march=native -fopenmp"
    "Clang Ofast Native"    "clang++ -Ofast -march=native -fopenmp" # Verify correctness
    "Clang O3 Native LTO"   "clang++ -O3 -march=native -fopenmp -flto"
    # --- GCC Section ---
    # Check your GCC version (e.g., g++-13, g++-14) and update command if needed
    "GCC O0"                "g++-14 -O0 -fopenmp"
    "GCC O2 Native"         "g++-14 -O2 -march=native -fopenmp"
    "GCC O3 Native"         "g++-14 -O3 -march=native -fopenmp"
    "GCC Ofast Native"      "g++-14 -Ofast -march=native -fopenmp" # Verify correctness
    "GCC O3 Native LTO"     "g++-14 -O3 -march=native -fopenmp -flto"
)

# --- Script Logic ---

echo "Starting compiler tests..."
# Clear previous results file
> "$RESULTS_FILE"
echo "Results will be saved to: $RESULTS_FILE"

# Temporary file for capturing output
TEMP_OUT=$(mktemp)
# Ensure temp file is cleaned up on exit
trap 'rm -f "$TEMP_OUT"' EXIT

# Check if sudo is needed (PMU often requires elevated privileges)
# Running the whole script with sudo might be necessary if the executable needs it.
if [[ $EUID -ne 0 ]]; then
   echo "WARNING: PMU access (-DPMU) might require root privileges." >&2
   echo "WARNING: Consider running this script with 'sudo bash $0'" >&2
fi

# Loop through configurations
for (( i=0; i<${#configs[@]}; i+=2 )); do
    title="${configs[i]}"
    compile_base_cmd="${configs[i+1]}"

    echo "---------------------------------------------"
    echo "Testing: $title"

    # Construct the full compile command
    compile_cmd="$compile_base_cmd ${SOURCE_FILES[*]} -o $OUTPUT_EXEC $DEFINE_FLAGS"

    echo "Compile command: $compile_cmd"

    # Compile the code
    if ! $compile_cmd; then
        echo "ERROR: Compilation failed for '$title'. Skipping." >&2
        echo "$title: COMPILE FAILED" >> "$RESULTS_FILE"
        continue # Skip to the next configuration
    fi
    echo "Compilation successful."

    # Construct the run command
    # NOTE: sudo is NOT added here. Run the whole script with sudo if needed.
    run_cmd="$OUTPUT_EXEC $RUN_ARGS"
    echo "Run command: $run_cmd"

    # Run the executable and capture output
    if ! $run_cmd > "$TEMP_OUT" 2>&1; then
        echo "ERROR: Execution failed for '$title'. Skipping." >&2
        cat "$TEMP_OUT" # Print output to see the error
        echo "$title: EXECUTION FAILED" >> "$RESULTS_FILE"
        rm -f "$OUTPUT_EXEC" # Clean up executable
        continue # Skip to the next configuration
    fi

    # Extract the cycle count (assumes first number on line containing "cycles")
    # Using awk for robustness: prints the first field ($1) of the line containing "cycles"
    cycles=$(grep 'cycles' "$TEMP_OUT" | awk '{print $1}')

    if [[ -z "$cycles" ]]; then
        echo "ERROR: Could not extract cycle count for '$title'." >&2
        cat "$TEMP_OUT" # Print output to help debug
        echo "$title: FAILED TO EXTRACT CYCLES" >> "$RESULTS_FILE"
    else
        echo "Cycles: $cycles"
        # Append result to file
        echo "$title: $cycles" >> "$RESULTS_FILE"
    fi

    # Clean up executable
    rm -f "$OUTPUT_EXEC"

done

echo "---------------------------------------------"
echo "Testing complete. Results saved in $RESULTS_FILE"

exit 0