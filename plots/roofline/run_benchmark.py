import subprocess
import re
import json
import argparse

def run_and_parse_benchmark(save_results=False):
    test_cases = [
        (   128,  512,  2048),
        (   128, 1024,  4096),
        (   128, 2048,  8192),
        (   128, 4096, 16384),
        ( 256,  512,  2048),
        ( 256, 1024,  4096),
        ( 256, 2048,  8192),
        ( 256, 4096, 16384),
    ]
    non_zero_s = 8
    executable_path = "../../SparseGEMM.out" 
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

            # Regex to find test case correctness (passed/failed)
            correctness_matches = re.findall(r"Test case (.*?) (passed|failed)!", stdout_output)
            correctness_status = {fn.strip(): status for fn, status in correctness_matches}
 
            # Regex to find function name and its cycle count
            # regex = r"Running: (.*?)\s*\n\s*([\d\.eE+-]+) cycles\s*\n\s*Performance: ([\d\.eE+-]+) flops/cycle"

            # Regex to capture the values after "Running:", "Performance:", and "Operational Intensity:"
            # It uses re.DOTALL to make '.' match newlines, allowing us to skip intermediate lines easily.
            regex = re.compile(
        r"Running:\s*(.*?)\s*\n"           # Capture Group 1: Text after "Running:"
        r".*?"                             # Non-greedy match for any lines in between
        r"Performance:\s*([\d\.eE+-]+)"    # Capture Group 2: Number after "Performance:"
        r".*?"                             # Non-greedy match for any lines in between
        r"Operational Intensity:\s*([\d\.eE+-]+)", # Capture Group 3: Number after "Operational Intensity:"
        re.DOTALL  # Make '.' match newline characters as well
    )
            matches = re.findall(regex, stdout_output)


            current_test_results = {}
            if matches:
                for func_name, fpc, oi in matches:
                    stripped_func_name = func_name.strip()
                    stripped_func_name = stripped_func_name.split('31m')[1]
                    stripped_func_name = stripped_func_name.split('\u001b')[0]
                    try:
                        fpc = float(fpc)
                        oi = float(oi)
                        current_test_results[stripped_func_name] = { 'operational_intensity' : oi, 'performance' : fpc }
                        print(f"  {stripped_func_name}: {fpc} flops/cycle , {oi} flops/Byte")
                        if correctness_status.get(stripped_func_name) == "failed":
                            print(f"  WARNING: {stripped_func_name} failed correctness check!")
                    except ValueError:
                        print(f"  Could not parse cycle count for {stripped_func_name}: {oi}")
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
    args = parser.parse_args()
    run_and_parse_benchmark(save_results=args.save) 