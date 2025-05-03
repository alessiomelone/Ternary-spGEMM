#!/bin/bash

declare -a test_cases=(
    "1 512 2048"
    "1 1024 4096"
    "1 2048 8192"
    "1 4096 16384"
    "256 512 2048"
    "256 1024 4096"
    "256 2048 8192"
    "256 4096 16384"
    "1000 4096 16384"
    "4000 4096 16384"
)

# Sparsity value
s=2

CSV_FILE="results/sparse_basic_$(date +%Y-%m-%d).csv"

mkdir -p results

echo "M,K,N,nonZero,clock_cycles,clock_freq_mhz,clock_seconds,timeofday_seconds,vct_cycles,vct_seconds,vct_freq_mhz,pmu_cycles,pmu_seconds,pmu_freq_mhz,pmu_instructions,pmu_branches,pmu_branch_misses,pmu_ipc" > "$CSV_FILE"

for test_case in "${test_cases[@]}"; do
    # Parse the test case
    read m k n <<< "$test_case"
    
    echo "Running benchmark with M=$m, K=$k, N=$n, nonZero=$s"
    
    # Compile the code
    OUTPUT=$(gcc -DPMU -arch arm64 -DVALIDATE -O3 ./c_impl/sparse_format.c ./c_impl/TestImpl.c -o ./c_impl/TestImpl.out)

    # Run the executable and capture output
    OUTPUT=$(sudo ./c_impl/TestImpl.out -M "$m" -K "$k" -N "$n" -s "$s" )
    
    # Check if test case failed
    if echo "$OUTPUT" | grep -qi "test case failed"; then
        echo "TESTCASE FAILED for M=$m, K=$k, N=$n, nonZero=$s"
        echo "$OUTPUT"
        continue  # Skip to next iteration instead of exiting
    fi

    # Parse output values
    clock_cycles=$(echo "$OUTPUT" | grep "clock_cycles=" | cut -d'=' -f2)
    clock_freq_mhz=$(echo "$OUTPUT" | grep "clock_freq_mhz=" | cut -d'=' -f2)
    clock_seconds=$(echo "$OUTPUT" | grep "clock_seconds=" | cut -d'=' -f2)
    timeofday_seconds=$(echo "$OUTPUT" | grep "timeofday_seconds=" | cut -d'=' -f2)
    vct_cycles=$(echo "$OUTPUT" | grep "vct_cycles=" | cut -d'=' -f2)
    vct_seconds=$(echo "$OUTPUT" | grep "vct_seconds=" | cut -d'=' -f2)
    vct_freq_mhz=$(echo "$OUTPUT" | grep "vct_freq_mhz=" | cut -d'=' -f2)
    pmu_cycles=$(echo "$OUTPUT" | grep "pmu_cycles=" | cut -d'=' -f2)
    pmu_seconds=$(echo "$OUTPUT" | grep "pmu_seconds=" | cut -d'=' -f2)
    pmu_freq_mhz=$(echo "$OUTPUT" | grep "pmu_freq_mhz=" | cut -d'=' -f2)
    pmu_instructions=$(echo "$OUTPUT" | grep "pmu_instructions=" | cut -d'=' -f2)
    pmu_branches=$(echo "$OUTPUT" | grep "pmu_branches=" | cut -d'=' -f2)
    pmu_branch_misses=$(echo "$OUTPUT" | grep "pmu_branch_misses=" | cut -d'=' -f2)
    pmu_ipc=$(echo "$OUTPUT" | grep "pmu_ipc=" | cut -d'=' -f2)
    
    # Check if any of the required values are missing
    if [ -z "$clock_cycles" ] || [ -z "$clock_freq_mhz" ] || [ -z "$clock_seconds" ] || [ -z "$timeofday_seconds" ]; then
        echo "WARNING: Missing required output values for M=$m, K=$k, N=$n, nonZero=$s"
        echo "$OUTPUT"
        continue  # Skip to next iteration instead of exiting
    fi

    echo "$m,$k,$n,$s,$clock_cycles,$clock_freq_mhz,$clock_seconds,$timeofday_seconds,$vct_cycles,$vct_seconds,$vct_freq_mhz,$pmu_cycles,$pmu_seconds,$pmu_freq_mhz,$pmu_instructions,$pmu_branches,$pmu_branch_misses,$pmu_ipc" >> "$CSV_FILE"
    
    echo "Completed benchmark for M=$m, K=$k, N=$n, nonZero=$s"
done

echo "Benchmarking complete. Results saved to $CSV_FILE."
