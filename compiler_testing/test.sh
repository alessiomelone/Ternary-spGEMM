#!/usr/bin/env bash
# RUN WITH SUDO

set -e
set -u
set -o pipefail

RUN_ARGS="-M 32 -K 1024 -N 4096 -s 4"


CPP_SRC_DIR="cpp_impl"
SOURCE_FILES=("$CPP_SRC_DIR/main.cpp" "$CPP_SRC_DIR/comp.cpp" "$CPP_SRC_DIR/perf.cpp")
OUTPUT_EXEC="$CPP_SRC_DIR/SparseGEMM.out"
RESULTS_FILE="compiler_testing/compiler_results.txt"
DEFINE_FLAGS="-DPMU"

COMPILERS=("g++")
FLAG_SETS=(
    "-O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG"
    "-O3 -march=native -mtune=native -fstrict-aliasing -DNDEBUG"
    "-Ofast -march=native -mtune=native -ffast-math -funroll-loops -fstrict-aliasing -DNDEBUG"

    "-O3 -march=native -mtune=native -flto=auto -fstrict-aliasing -DNDEBUG"
    "-O3 -march=native -mtune=native -flto=auto -funroll-loops -fstrict-aliasing -DNDEBUG"
    "-Ofast -march=native -mtune=native -flto=auto -ffast-math -funroll-loops -fstrict-aliasing -DNDEBUG"
)

echo "Starting compiler tests..."
> "$RESULTS_FILE"
echo "Results will be saved to: $RESULTS_FILE"

best_cycles=""
best_title=""

TEMP_OUT=$(mktemp)
trap 'rm -f "$TEMP_OUT"' EXIT

if [[ $EUID -ne 0 ]]; then
   echo "WARNING: PMU access (-DPMU) might require root privileges." >&2
   echo "WARNING: Consider running this script with 'sudo bash $0'" >&2
fi

for compiler in "${COMPILERS[@]}"; do
    for flags in "${FLAG_SETS[@]}"; do
        title="$compiler ${flags// /_}"
        compile_base_cmd="$compiler $flags"

        echo "---------------------------------------------"
        echo "Testing: $title"

        compile_cmd="$compile_base_cmd ${SOURCE_FILES[*]} -o $OUTPUT_EXEC $DEFINE_FLAGS"

        echo "Compile command: $compile_cmd"

        if ! $compile_cmd > "$TEMP_OUT" 2>&1; then
            echo "$title: COMPILE FAILED" >> "$RESULTS_FILE"
            echo "Compiler output for $title:" >> "$RESULTS_FILE"
            cat "$TEMP_OUT" >> "$RESULTS_FILE"
            echo "-----------------------" >> "$RESULTS_FILE"
            continue
        fi

        echo "Compilation successful."

        run_cmd="$OUTPUT_EXEC $RUN_ARGS"
        echo "Run command: $run_cmd"

        if ! $run_cmd > "$TEMP_OUT" 2>&1; then
            echo "ERROR: Execution failed for '$title'. Skipping." >&2
            cat "$TEMP_OUT"
            echo -e "$title:\nEXECUTION FAILED" >> "$RESULTS_FILE"
            rm -f "$OUTPUT_EXEC"
            echo "-----------------------" >> "$RESULTS_FILE"
            continue
        fi

        vals=( $(grep 'cycles' "$TEMP_OUT" | awk '{print $1}') )
        if [[ ${#vals[@]} -lt 2 ]]; then
            cycles="${vals[0]}"
        else
            cycles=$(awk -v a="${vals[0]}" -v b="${vals[1]}" 'BEGIN{if (a<b) print a; else print b}')
        fi

        if [[ -z "$cycles" ]]; then
            echo "ERROR: Could not extract cycle count for '$title'." >&2
            cat "$TEMP_OUT"
            echo "$title: FAILED TO EXTRACT CYCLES" >> "$RESULTS_FILE"
        else
            echo "Cycles: $cycles"
            echo -e "$title:\n$cycles" >> "$RESULTS_FILE"

            if [[ -z "$best_cycles" ]]; then
                best_cycles="$cycles"
                best_title="$title"
            else
                is_less=$(awk -v a="$cycles" -v b="$best_cycles" 'BEGIN{print (a<b)}')
                if [[ "$is_less" -eq 1 ]]; then
                    best_cycles="$cycles"
                    best_title="$title"
                fi
            fi
        fi

        rm -f "$OUTPUT_EXEC"
        echo "-----------------------" >> "$RESULTS_FILE"
    done
done

echo
echo -e "Fastest configuration:\n$best_title with $best_cycles cycles"
echo -e "Fastest configuration:\n$best_title with $best_cycles cycles" >> "$RESULTS_FILE"

echo "---------------------------------------------"
echo "Testing complete. Results saved in $RESULTS_FILE"

exit 0