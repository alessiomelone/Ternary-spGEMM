#!/usr/bin/env bash

set -e
set -u
set -o pipefail

CPP_SRC_DIR="cpp_impl"
SOURCE_FILES=("$CPP_SRC_DIR/main.cpp" "$CPP_SRC_DIR/comp.cpp" "$CPP_SRC_DIR/perf.cpp")
OUTPUT_EXEC="$CPP_SRC_DIR/SparseGEMM.out"
RESULTS_FILE="compiler_testing/compiler_results.txt"
RUN_ARGS="-M 32 -K 1024 -N 4096 -s 4"
DEFINE_FLAGS="-DPMU"

declare -a configs=(
    "AppleClang O3 M2"               "clang++ -O3 -mcpu=apple-m2 -mtune=apple-m2 -fstrict-aliasing -DNDEBUG"
    "AppleClang Ofast M2"            "clang++ -Ofast -mcpu=apple-m2 -mtune=apple-m2 -ffast-math -fstrict-aliasing -DNDEBUG"
    "AppleClang O3 M2 LTO"           "clang++ -O3 -mcpu=apple-m2 -mtune=apple-m2 -flto=thin -fstrict-aliasing -DNDEBUG"
    "LLVM Clang-15 O2 M2 OpenMP"     "clang++-15 -O2 -mcpu=apple-m2 -mtune=apple-m2 -fopenmp -fstrict-aliasing -DNDEBUG"
    "LLVM Clang-15 O3 M2 OpenMP"     "clang++-15 -O3 -mcpu=apple-m2 -mtune=apple-m2 -fopenmp -flto=thin -ffast-math"
    "LLVM Clang-15 Ofast M2 OpenMP"  "clang++-15 -Ofast -mcpu=apple-m2 -mtune=apple-m2 -fopenmp -ffast-math"
    "GCC-14 O3 M2 OpenMP"            "g++-14 -O3 -march=armv8.4-a -fopenmp -flto=auto -fstrict-aliasing -DNDEBUG"
    "GCC-14 Ofast M2 OpenMP"         "g++-14 -Ofast -march=armv8.4-a -fopenmp -funroll-loops -ffast-math"
    "LLVM Clang-15 PGO Generate"     "clang++-15 -O2 -mcpu=apple-m2 -mtune=apple-m2 -fopenmp -fprofile-generate"
    "LLVM Clang-15 PGO Use"          "clang++-15 -O2 -mcpu=apple-m2 -mtune=apple-m2 -fopenmp -fprofile-use"
    "Clang AddressSanitizer"         "clang++ -O1 -g -fsanitize=address -fno-omit-frame-pointer"
    "Clang UndefinedSanitizer"       "clang++ -O1 -g -fsanitize=undefined -fno-omit-frame-pointer"
)

echo "Starting compiler tests..."
> "$RESULTS_FILE"
echo "Results will be saved to: $RESULTS_FILE"

TEMP_OUT=$(mktemp)
trap 'rm -f "$TEMP_OUT"' EXIT

if [[ $EUID -ne 0 ]]; then
   echo "WARNING: PMU access (-DPMU) might require root privileges." >&2
   echo "WARNING: Consider running this script with 'sudo bash $0'" >&2
fi

for (( i=0; i<${#configs[@]}; i+=2 )); do
    title="${configs[i]}"
    compile_base_cmd="${configs[i+1]}"

    echo "---------------------------------------------"
    echo "Testing: $title"

    compile_cmd="$compile_base_cmd ${SOURCE_FILES[*]} -o $OUTPUT_EXEC $DEFINE_FLAGS"

    echo "Compile command: $compile_cmd"

    if ! $compile_cmd; then
        echo "ERROR: Compilation failed for '$title'. Skipping." >&2
        echo "$title: COMPILE FAILED" >> "$RESULTS_FILE"
        continue
    fi
    echo "Compilation successful."

    run_cmd="$OUTPUT_EXEC $RUN_ARGS"
    echo "Run command: $run_cmd"

    if ! $run_cmd > "$TEMP_OUT" 2>&1; then
        echo "ERROR: Execution failed for '$title'. Skipping." >&2
        cat "$TEMP_OUT"
        echo "$title: EXECUTION FAILED" >> "$RESULTS_FILE"
        rm -f "$OUTPUT_EXEC"
        continue
    fi

    cycles=$(grep 'cycles' "$TEMP_OUT" | awk '{print $1}')

    if [[ -z "$cycles" ]]; then
        echo "ERROR: Could not extract cycle count for '$title'." >&2
        cat "$TEMP_OUT"
        echo "$title: FAILED TO EXTRACT CYCLES" >> "$RESULTS_FILE"
    else
        echo "Cycles: $cycles"
        echo "$title: $cycles" >> "$RESULTS_FILE"
    fi

    rm -f "$OUTPUT_EXEC"

done

echo "---------------------------------------------"
echo "Testing complete. Results saved in $RESULTS_FILE"

exit 0