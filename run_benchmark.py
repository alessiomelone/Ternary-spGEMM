import subprocess
import re
import json

def run_and_parse_benchmark():
    test_cases = [
        (   1,  512,  2048),
        (   1, 1024,  4096),
        (   1, 2048,  8192),
        (   1, 4096, 16384),
        ( 256,  512,  2048),
        ( 256, 1024,  4096),
        ( 256, 2048,  8192),
        ( 256, 4096, 16384),
    ]
    non_zero_s = 2
    # Assuming SparseGEMM.out is in the cpp_impl subdirectory relative to this script
    executable_path = "./cpp_impl/SparseGEMM.out" 
    base_command = ["sudo", executable_path]

    all_results = []

    print(f"Running benchmark with executable: {executable_path}")
    print(f"Using -s (nonZero) parameter: {non_zero_s}\n")

    for m_val, k_val, n_val in test_cases:
        print(f"--- Running test case: M={m_val}, K={k_val}, N={n_val}, s={non_zero_s} ---")
        command = base_command + [
            "-M", str(m_val),
            "-K", str(k_val),
            "-N", str(n_val),
            "-s", str(non_zero_s)
        ]

        try:
            # print(f"Executing command: {' '.join(command)}") # Uncomment for debugging command
            process = subprocess.run(command, capture_output=True, text=True, check=False, timeout=300) # 5 min timeout

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
            # print(f"Raw stdout:\n{stdout_output[:500]}...") # Uncomment for debugging output

            # Regex to find function name and its cycle count
            # Handles various number formats including scientific notation
            # and skips the "Starting program. X functions registered." and "Test case ... passed/failed!" lines
            matches = re.findall(r"Running: (.*?)\s*\n\s*([\d\.eE+-]+) cycles", stdout_output)
            
            # Regex to find test case correctness (passed/failed)
            correctness_matches = re.findall(r"Test case (.*?) (passed|failed)!", stdout_output)
            correctness_status = {fn.strip(): status for fn, status in correctness_matches}

            current_test_results = {}
            if matches:
                for func_name, cycles_str in matches:
                    stripped_func_name = func_name.strip()
                    try:
                        cycles = float(cycles_str) # Convert cycles to float
                        current_test_results[stripped_func_name] = cycles
                        print(f"  {stripped_func_name}: {cycles:.2e} cycles")
                        if correctness_status.get(stripped_func_name) == "failed":
                            print(f"  WARNING: {stripped_func_name} failed correctness check!")
                    except ValueError:
                        print(f"  Could not parse cycle count for {stripped_func_name}: {cycles_str}")
                        current_test_results[stripped_func_name] = "Error parsing cycles"
            else:
                print("  No performance results found in output.")
                # print(f"Full stdout for M={m_val}, K={k_val}, N={n_val}:\n{stdout_output}")


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
    output_file = "benchmark_results.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=4)
    print(f"\nAll benchmark results saved to {output_file}")

if __name__ == "__main__":
    run_and_parse_benchmark() 