#!/bin/bash

M=(1 16 64 256 1000 4000 16000 64000)
K=(512 1024 2048 4096 2048 4096 8192 16384)
N=(2048 4096 8192 16384 512 1024 2048 4096)
nonZero=(2 4 8 16)

CSV_FILE="results/sparse_fixed_$(date +%Y-%m-%d).csv"

mkdir -p results

echo "M,K,N,nonZero,clock_cycles,clock_freq_mhz,clock_seconds,timeofday_seconds,vct_cycles,vct_seconds,vct_freq_mhz,pmu_cycles,pmu_seconds,pmu_freq_mhz,pmu_instructions,pmu_branches,pmu_branch_misses,pmu_ipc" > "$CSV_FILE"

# Loop over all 8*8 permutations
for m_idx in {0..7}; do
    kn_idx=0
    m=${M[$m_idx]}
    k=${K[$kn_idx]}
    n=${N[$kn_idx]} 
    
    for s in "${nonZero[@]}"; do
        echo "Running benchmark with M=$m, K=$k, N=$n, nonZero=$s"
        
        # Run the executable and capture output
        OUTPUT=$(sudo ./c_impl/TestImpl.out -M "$m" -K "$k" -N "$n" -s "$s")
        
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
done

echo "Benchmarking complete. Results saved to $CSV_FILE."
