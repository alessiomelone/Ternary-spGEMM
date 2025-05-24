import subprocess
import re
import json
import argparse
import os

def ensure_make():
    print("Running make to ensure latest build...")
    try:
        subprocess.run(["make"], check=True)
        print("Make completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running make: {e}")
        raise

def run_and_parse_benchmark(save_results=False, raw_output=False):
    # Ensure make has been run
    ensure_make()
    
    test_cases = [
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 205dd08 (BlockedCSC is still wrong, I need some sleep)
        (   1,  1,  8),
        (   1,  2,  8),
        (   1,  4,  8),
        (   1,  8,  8),
        (   1,  16,  8),
        (   1,  32,  8),
        (   1,  64,  8),
        (   1,  128,  8),
        (   1,  256,  8),
        (   1,  512,  8),
        (   1, 1024,  8),
        (   1, 2048,  8),
        (   1, 4096, 8),
        (   1, 8192,  8),
        (   1, 16384,  8),
        (   1, 32768,  8),
        (   1, 65536,  8),
        (   1, 131072,  8),
        (   1, 262144,  8),
        (   1, 524288,  8),
        (   1, 1048576,  8),
        (   1, 2097152,  8),
        (   1, 4194304,  8),
        (   1, 8388608,  8),
        (   1, 16777216,  8),
        (   1, 33554432,  8),
        (   1, 67108864,  8),
        (   1, 134217728,  8),
<<<<<<< HEAD
=======
        (   1,  512,  2048),
        (   1, 1024,  4096),
        (   1, 2048,  8192),
        (   1, 4096, 16384),
        ( 256,  512,  2048),
>>>>>>> a4043d8 (added calibrate as flag (Default CALIBRATE), removed Baraq's unrolled 5 (templated))
=======
>>>>>>> 205dd08 (BlockedCSC is still wrong, I need some sleep)
        # ( 256, 1024,  4096),
        # ( 256, 2048,  8192)
    ]
    non_zero_values = [2, 4, 8, 16]
    executable_path = "./SparseGEMM.out" 
    base_command = ["sudo", executable_path]

    all_results = []

    print(f"Running benchmark with executable: {executable_path}")
    print(f"Testing non-zero values: {non_zero_values}\n")

    for non_zero_s in non_zero_values:
        print(f"\n=== Testing with non-zero value: {non_zero_s} ===\n")
        
        for m_val, k_val, n_val in test_cases:
            print(f"--- Running test case: M={m_val}, K={k_val}, N={n_val}, s={non_zero_s} ---")
            command = base_command + [
                "-M", str(m_val),
                "-K", str(k_val),
                "-N", str(n_val),
                "-s", str(non_zero_s)
            ]

            try:
                process = subprocess.run(command, capture_output=True, text=True, check=False) # NO TIMEOUT

                if process.returncode != 0:
                    print(f"ERROR: Benchmark run failed for M={m_val}, K={k_val}, N={n_val}")
                    print(f"Return code: {process.returncode}")
                    print(f"Stderr:\n{process.stderr}")
                    all_results.append({
                        "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": non_zero_s},
                        "error": process.stderr,
                        "results": {}
                    })
                    continue

                stdout_output = process.stdout

                if raw_output:
                    print("\nRaw output:")
                    print(stdout_output)
                    print("-" * 80)
                    continue

                # Regex to find test case correctness (passed/failed)
                correctness_matches = re.findall(r"Test case (.*?) (passed|failed)!", stdout_output)
                correctness_status = {fn.strip(): status for fn, status in correctness_matches}
 
                # Regex to find function name and its cycle count
                matches = re.findall(r"Running: (.*?)\s*\n\s*([\d\.eE+-]+) cycles", stdout_output)
                
                current_test_results = {}
                if matches:
                    for func_name, cycles_str in matches:
                        stripped_func_name = func_name.strip()
                        try:
                            cycles = float(cycles_str)
                            current_test_results[stripped_func_name] = cycles
                            print(f"  {stripped_func_name}: {cycles:.2e} cycles")
                            if correctness_status.get(stripped_func_name) == "failed":
                                print(f"  WARNING: {stripped_func_name} failed correctness check!")
                        except ValueError:
                            print(f"  Could not parse cycle count for {stripped_func_name}: {cycles_str}")
                            current_test_results[stripped_func_name] = "Error parsing cycles"
                else:
                    print("  No performance results found in output.")


                all_results.append({
                    "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": non_zero_s},
                    "results": current_test_results
                })

            except subprocess.TimeoutExpired:
                print(f"ERROR: Benchmark run timed out for M={m_val}, K={k_val}, N={n_val}")
                all_results.append({
                    "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": non_zero_s},
                    "error": "TimeoutExpired",
                    "results": {}
                })
            except Exception as e:
                print(f"An unexpected error occurred for M={m_val}, K={k_val}, N={n_val}: {e}")
                all_results.append({
                    "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": non_zero_s},
                    "error": str(e),
                    "results": {}
                })
            print("-" * (len(f"--- Running test case: M={m_val}, K={k_val}, N={n_val}, s={non_zero_s} ---")) + "\n")

    # Optionally, save all results to a JSON file
    if save_results:
        output_file = "benchmark_results.json"
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=4)
        print(f"\nAll benchmark results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run SparseGEMM benchmarks')
    parser.add_argument('-s', '--save', action='store_true', help='Save benchmark results to JSON file')
    parser.add_argument('-r', '--raw', action='store_true', help='Show raw output of test cases without parsing')
    args = parser.parse_args()
    run_and_parse_benchmark(save_results=args.save, raw_output=args.raw) 