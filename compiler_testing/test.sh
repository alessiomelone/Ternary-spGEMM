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

declare -a configs=(
    "AppleClang O3 M2"               "clang++ -O3 -mcpu=native -fstrict-aliasing -DNDEBUG"
    "AppleClang Ofast M2"            "clang++ -Ofast -mcpu=native -ffast-math -fstrict-aliasing -DNDEBUG"
    "AppleClang O3 M2 LTO"           "clang++ -O3 -mcpu=native -flto=thin -fstrict-aliasing -DNDEBUG"
    # "LLVM Clang-15 O2 M2 OpenMP"     "clang++-15 -O2 -mcpu=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -fstrict-aliasing -DNDEBUG"
    # "LLVM Clang-15 O3 M2 OpenMP"     "clang++-15 -O3 -mcpu=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -flto=thin -ffast-math"
    # "LLVM Clang-15 Ofast M2 OpenMP"  "clang++-15 -Ofast -mcpu=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -ffast-math"
    # "GCC-14 O3 M2 OpenMP"            "g++-14 -O3 -march=armv8.4-a -fopenmp -flto=auto -fstrict-aliasing -DNDEBUG"
    # "GCC-14 Ofast M2 OpenMP"         "g++-14 -Ofast -march=armv8.4-a -fopenmp -funroll-loops -ffast-math"
    # "LLVM Clang-15 PGO Generate"     "clang++-15 -O2 -mcpu=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -fprofile-generate"
    # "LLVM Clang-15 PGO Use"          "clang++-15 -O2 -mcpu=native -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -fprofile-use"
    "Clang AddressSanitizer"         "clang++ -O1 -g -fsanitize=address -fno-omit-frame-pointer"
    "Clang UndefinedSanitizer"       "clang++ -O1 -g -fsanitize=undefined -fno-omit-frame-pointer"
    "AppleClang O3 apple-m2"               "clang++ -O3 -mcpu=native -ffp-contract=fast -fstrict-aliasing -DNDEBUG"
    "AppleClang Ofast apple-m2 Unroll"     "clang++ -Ofast -mcpu=native -funroll-loops -ffast-math -fstrict-aliasing -DNDEBUG"
    # "AppleClang O3 ThinLTO apple-m2"       "clang++ -O3 -mcpu=native -flto=thin -fvectorize -fstrict-aliasing -DNDEBUG"
    # "AppleClang O3 OpenMP apple-m2"        "clang++ -O3 -mcpu=native -fopenmp -ffp-contract=fast -fstrict-aliasing -DNDEBUG"
    # "AppleClang Ofast LTO OpenMP"          "clang++ -Ofast -mcpu=apple-m2 -fopenmp -flto=thin -ffast-math -funroll-loops"
    # "AppleClang PGO Generate"              "clang++ -O2 -mcpu=apple-m2 -fprofile-instr-generate -fcoverage-mapping -fopenmp"
    # "AppleClang PGO Use ThinLTO"           "clang++ -O3 -mcpu=apple-m2 -flto=thin -fprofile-instr-use=default.profdata -fopenmp -ffp-contract=fast"
    "AppleClang Bolt-prep"                 "clang++ -O3 -mcpu=apple-m2 -flto=auto -fdata-sections -ffunction-sections -Wl,-dead_strip -DNDEBUG"
    "AppleClang Bolt-optimised"            "llvm-bolt a.out -o a.bolt --deterministic --lite=0 -reorder-blocks=cache+ --reorder-functions=hfsort"
    # "LLVM17 O3 apple-m2 Vectorize"         "clang++-17 -O3 -mcpu=apple-m2 -moutline-atomics -ffp-contract=fast -fstrict-aliasing -DNDEBUG"
    # "LLVM17 Ofast OpenMP apple-m2"         "clang++-17 -Ofast -mcpu=apple-m2 -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -ffast-math"
    # "GCC14 O3 armv8.5-a OpenMP"            "g++-14 -O3 -march=armv8.5-a+crc+crypto+dotprod -fopenmp -flto=auto -fstrict-aliasing -DNDEBUG"
    # "GCC14 Ofast Unroll OpenMP"            "g++-14 -Ofast -march=armv8.5-a+crc+crypto+dotprod -fopenmp -funroll-loops -ffast-math"
    "Clang ASan apple-m2 FastDebug"        "clang++ -O1 -gline-tables-only -fsanitize=address -mcpu=native -fomit-frame-pointer"
    "Clang UBSan apple-m2"                 "clang++ -O1 -g -fsanitize=undefined -mcpu=native -fno-omit-frame-pointer"
)

echo "Starting compiler tests..."
> "$RESULTS_FILE"
echo "Results will be saved to: $RESULTS_FILE"

# Initialize best result tracking
best_cycles=""
best_title=""

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

    # Extract cycle counts and choose the minimum of the two runs
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
        # Track fastest cycle count
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

# Summarize fastest result
echo
echo -e "Fastest configuration:\n$best_title with $best_cycles cycles"
echo -e "Fastest configuration:\n$best_title with $best_cycles cycles" >> "$RESULTS_FILE"

echo "---------------------------------------------"
echo "Testing complete. Results saved in $RESULTS_FILE"

exit 0