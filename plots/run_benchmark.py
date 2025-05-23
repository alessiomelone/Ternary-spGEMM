import subprocess
import re
import json
import argparse

def run_and_parse_benchmark(vary='M', sparsity=8, save_results=False, outname=''):
    reference_values = [1, 16, 64, 256, 1000, 4000]#, 16000, 64000]

    test_cases = []
    if vary == 'M':
        for m in reference_values:
            test_cases.append((m, 1024, 1024))
    elif vary == 'K':
        for k in reference_values:
            test_cases.append((1024, k, 1024))
    else:  # vary == 'N'
        for n in reference_values:
            test_cases.append((1024, 1024, n))

    executable_path = "./SparseGEMM.out" 
    base_command = ["sudo", executable_path]

    all_results = []

    print(f"Running benchmark with executable: {executable_path}")
    print(f"Using sparsity {sparsity}\n")

    for m_val, k_val, n_val in test_cases:
        print(f"--- Running test case: M={m_val}, K={k_val}, N={n_val}, s={sparsity} ---")
        command = base_command + [
            "-M", str(m_val),
            "-K", str(k_val),
            "-N", str(n_val),
            "-s", str(sparsity)
        ]

        try:
            process = subprocess.run(command, capture_output=True, text=True, check=False) # NO TIMEOUT

            if process.returncode != 0:
                print(f"ERROR: Benchmark run failed for M={m_val}, K={k_val}, N={n_val}")
                print(f"Return code: {process.returncode}")
                print(f"Stderr:\n{process.stderr}")
                all_results.append({
                    "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": sparsity},
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
        r"Total Input Size:\s*([\d\.eE+-]+)"    # Capture Group 2: Number after "Performance:"
        r".*?"                             # Non-greedy match for any lines in between
        r"Operational Intensity:\s*([\d\.eE+-]+)", # Capture Group 3: Number after "Operational Intensity:"
        re.DOTALL  # Make '.' match newline characters as well
    )
            matches = re.findall(regex, stdout_output)


            current_test_results = {}
            if matches:
                for func_name, fpc, total_input_size, oi in matches:
                    stripped_func_name = func_name.strip()
                    stripped_func_name = stripped_func_name.split('31m')[1]
                    stripped_func_name = stripped_func_name.split('\u001b')[0]
                    try:
                        fpc = float(fpc)
                        oi = float(oi)
                        total_input_sz = int(float(total_input_size))
                        current_test_results[stripped_func_name] = { 'total_input_size' : total_input_sz, 'operational_intensity' : oi, 'performance' : fpc }
                        print(f"  {stripped_func_name}: {total_input_sz} Bytes input, {fpc} flops/cycle, {oi} flops/Byte")
                        if correctness_status.get(stripped_func_name) == "failed":
                            print(f"  WARNING: {stripped_func_name} failed correctness check!")
                    except ValueError:
                        print(f"  Could not parse cycle count for {stripped_func_name}: {oi} : {total_input_size}")
                        current_test_results[stripped_func_name] = "Error parsing cycles"
            else:
                print("  No performance results found in output.")


            all_results.append({
                "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": sparsity},
                "results": current_test_results
            })

        except subprocess.TimeoutExpired:
            print(f"ERROR: Benchmark run timed out for M={m_val}, K={k_val}, N={n_val}")
            all_results.append({
                "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": sparsity},
                "error": "TimeoutExpired",
                "results": {}
            })
        except Exception as e:
            print(f"An unexpected error occurred for M={m_val}, K={k_val}, N={n_val}: {e}")
            all_results.append({
                "test_case": {"M": m_val, "K": k_val, "N": n_val, "s": sparsity},
                "error": str(e),
                "results": {}
            })
        print("-" * (len(f"--- Running test case: M={m_val}, K={k_val}, N={n_val}, s={sparsity} ---")) + "\n")

    # Optionally, save all results to a JSON file
    if save_results:
        output_file = outname
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=4)
        print(f"\nAll benchmark results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run SparseGEMM benchmarks')
    parser.add_argument('-s', '--save', action='store_true', help='Save benchmark results to JSON file')
    parser.add_argument('--output', type=str, default='benchmark_results.json', help='Name of JSON file')
    parser.add_argument('--vary', choices=['M', 'K', 'N'], default='M',
                        help='Matrix dimension to vary; others fixed at 1024')
    parser.add_argument('--sparsity', type=int, default=8,
                        help='Non‑zero elements per row passed to SparseGEMM via -s')
    args = parser.parse_args()
    run_and_parse_benchmark(args.vary, args.sparsity, args.save, args.output)