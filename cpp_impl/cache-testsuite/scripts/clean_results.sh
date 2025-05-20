#!/bin/bash

# --- Default values ---
# These can be overridden by CLI arguments
TO_BE_SUBTRACTED_DEFAULT="./measurements/raw/gemm.txt"
BASE_RESULTS_DEFAULT="./measurements/raw/base.txt"
OUTPUT_FILE_DEFAULT="./measurements/clean/clean.txt"

# --- Script paths (assuming they are in these locations relative to the script or in PATH) ---
# If these scripts are not in the same directory or in PATH, provide full paths.
SUB_SCRIPT="./scripts/sub_cache_res.py"
CACHE_ANLZ="./scripts/cache_miss_analyzer.py"

# --- Temporary file ---
# Using mktemp for a more robust temporary file creation
TEMP_FILE=$(mktemp)
# Ensure TEMP_FILE is removed on script exit (normal or error)
trap 'rm -f "$TEMP_FILE"' EXIT

# --- Function to display usage ---
usage() {
    echo "Usage: $0 [--sub <subtracted_file>] [--base <base_file>] [--out <output_file>]"
    echo "Options:"
    echo "  --sub FILE    Path to the file from which base results will be subtracted (default: $TO_BE_SUBTRACTED_DEFAULT)"
    echo "  --base FILE   Path to the base results file (default: $BASE_RESULTS_DEFAULT)"
    echo "  --out FILE    Path to the output file (default: $OUTPUT_FILE_DEFAULT)"
    echo "  -h, --help    Display this help message"
    exit 1
}

# --- Parse Command Line Arguments ---
# Initialize variables with default values
sub_file="$TO_BE_SUBTRACTED_DEFAULT"
base_file="$BASE_RESULTS_DEFAULT"
output_file="$OUTPUT_FILE_DEFAULT"

# Using a loop to handle long options like --sub
# This is a common way to handle long options in bash without `getopt` (the GNU version)
# which is not always available or behaves differently.
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sub)
            sub_file="$2"
            shift # past argument
            shift # past value
            ;;
        --base)
            base_file="$2"
            shift # past argument
            shift # past value
            ;;
        --out)
            output_file="$2"
            shift # past argument
            shift # past value
            ;;
        -h|--help)
            usage
            ;;
        *)
            # Unknown option
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# --- Validate that the Python scripts exist and are executable ---
if [ ! -f "$SUB_SCRIPT" ] ; then
    echo "Error: Subtraction script '$SUB_SCRIPT' not found or not executable."
    exit 1
fi

if [ ! -f "$CACHE_ANLZ" ] ; then
    echo "Error: Cache analysis script '$CACHE_ANLZ' not found or not executable."
    exit 1
fi

# --- Outputting the configuration being used ---
echo "--- Configuration ---"
echo "Subtraction Script: $SUB_SCRIPT"
echo "Cache Analyzer Script: $CACHE_ANLZ"
echo "File to be Subtracted (Input): $sub_file"
echo "Base Results File (Input): $base_file"
echo "Temporary File: $TEMP_FILE"
echo "Final Output File: $output_file"
echo "---------------------"
echo ""

# --- Check if input files exist ---
if [ ! -f "$sub_file" ]; then
    echo "Error: Subtracted file '$sub_file' not found."
    exit 1
fi

if [ ! -f "$base_file" ]; then
    echo "Error: Base file '$base_file' not found."
    exit 1
fi

# --- Ensure output directory exists ---
output_dir=$(dirname "$output_file")
if [ ! -d "$output_dir" ]; then
    echo "Output directory '$output_dir' does not exist. Creating it."
    mkdir -p "$output_dir"
    if [ $? -ne 0 ]; then
        echo "Error: Could not create output directory '$output_dir'."
        exit 1
    fi
fi


# --- Main script execution ---
echo "Step 1: Running subtraction script..."
# The arguments --fields "hits" "miss" --excl "compul" "pref" are kept as in your example.
# If these also need to be configurable, the script would need further modification.
python3 "$SUB_SCRIPT" "$sub_file" "$base_file" --fields "hits" "miss" --excl "compul" "pref" > "$TEMP_FILE"

if [ $? -ne 0 ]; then
    echo "Error: Subtraction script failed."
    # rm -f "$TEMP_FILE" # Already handled by trap
    exit 1
fi
echo "Subtraction script completed. Intermediate results in $TEMP_FILE"

echo "Step 2: Running cache analyzer script..."
python3 "$CACHE_ANLZ" "$TEMP_FILE" > "$output_file"

if [ $? -ne 0 ]; then
    echo "Error: Cache analyzer script failed."
    # rm -f "$TEMP_FILE" # Already handled by trap
    exit 1
fi
echo "Cache analyzer script completed. Final output in $output_file"

# Temporary file is automatically removed by the trap
# rm -f "$TEMP_FILE"
echo "Script finished successfully."

exit 0
