import subprocess
import re
import json
import argparse

def strip_ansi_codes(text):
    """Remove ANSI escape sequences (color codes) from text"""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)

def run_and_parse_benchmark(save_results=False):
    # Test case configurations from the provided arrays
    M_values = [1, 16, 64, 256, 1000, 4000, 16000, 64000]
    K_values = [512, 1024, 2048, 4096, 2048, 4096, 8192, 16384]
    N_values = [2048, 4096, 8192, 16384, 512, 1024, 2048, 4096]
    
    # Generate test cases: each M with each corresponding (K,N) pair
    test_cases = []
    for m_val in M_values:
        for k_val, n_val in zip(K_values, N_values):
            test_cases.append((m_val, k_val, n_val))
    
    non_zero_s = 2
    executable_path = "./sparseGEMM.out" 
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
            
            # Strip ANSI color codes from output before parsing
            clean_output = strip_ansi_codes(stdout_output)

            # Regex to find test case correctness (passed/failed)
            correctness_matches = re.findall(r"Test case (.*?) (passed|failed)!", clean_output)
            correctness_status = {fn.strip(): status for fn, status in correctness_matches}
 
            # Regex to find function name, cycle count, and speedup
            matches = re.findall(r"Running: (.*?)\s*\n\s*([\d\.eE+-]+) cycles\s*\n\s*Speedup is: ([\d\.eE+-]+)", clean_output)
            
            current_test_results = {}
            if matches:
                for func_name, cycles_str, speedup_str in matches:
                    stripped_func_name = func_name.strip()
                    try:
                        cycles = float(cycles_str)
                        speedup = float(speedup_str)
                        current_test_results[stripped_func_name] = {"cycles": cycles, "speedup": speedup}
                        print(f"  {stripped_func_name}: {speedup:.4f}x speedup")
                        if correctness_status.get(stripped_func_name) == "failed":
                            print(f"  WARNING: {stripped_func_name} failed correctness check!")
                    except ValueError:
                        print(f"  Could not parse results for {stripped_func_name}: cycles={cycles_str}, speedup={speedup_str}")
                        current_test_results[stripped_func_name] = "Error parsing results"
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
    args = parser.parse_args()
    run_and_parse_benchmark(save_results=args.save) 