import subprocess
import re
import json
import argparse

def run_and_parse_benchmark(save_results=False, outname='', varyonly=None):
    # If only one dimension should vary, define dimension lists and defaults
    list_M = [1, 16, 64, 256, 1000, 4000, 16000, 64000]
    list_K = [512, 1024, 2048, 4096, 8192, 16384]
    list_N = [512, 1024, 2048, 4096, 8192, 16384]
    default_val = 1024

    if varyonly == 'M':
        test_cases = [(m, default_val, default_val) for m in list_M]
    elif varyonly == 'K':
        test_cases = [(default_val, k, default_val) for k in list_K]
    elif varyonly == 'N':
        test_cases = [(default_val, default_val, n) for n in list_N]
    else:
        # Original full list of test cases
        test_cases = [
            (1, 512, 2048),
            (16, 1024, 4096),
            (64, 2048, 8192),
            (256, 4096, 16384),
            (1000, 2048, 512),
            (4000, 4096, 1024),
            (16000, 8192, 2048),
            (64000, 16384, 4096)
        ]

    list_s = [2, 4, 8, 16]
    executable_path = "./SparseGEMM.out" 
    base_command = ["sudo", executable_path]

    all_results = []

    print(f"Running benchmark with executable: {executable_path}")
    print(f"Using -s (nonZero) parameter: {list_s}\n")

    for m_val, k_val, n_val in test_cases:
        case_result = {
            "test_case": {"M": m_val, "K": k_val, "N": n_val},
            "results": {}
        }
        header = f"--- Running test case: M={m_val}, K={k_val}, N={n_val} ---"
        print(header)
        for s in list_s:
            command = base_command + [
                "-M", str(m_val),
                "-K", str(k_val),
                "-N", str(n_val),
                "-s", str(s)
            ]
            try:
                process = subprocess.run(command, capture_output=True, text=True, check=False)
                if process.returncode != 0:
                    print(f"  ERROR: Benchmark run failed for M={m_val}, K={k_val}, N={n_val}, s={s}")
                    print(f"  Return code: {process.returncode}")
                    print(f"  Stderr:\n{process.stderr}")
                    continue

                stdout_output = process.stdout

                correctness_matches = re.findall(r"Test case (.*?) (passed|failed)!", stdout_output)
                correctness_status = {fn.strip(): status for fn, status in correctness_matches}

                regex = re.compile(
                    r"Running:\s*(.*?)\s*\n"
                    r".*?"
                    r"Performance:\s*([\d\\.eE+-]+)"
                    r".*?"
                    r"Total Input Size:\s*([\d\\.eE+-]+)"
                    r".*?"
                    r"Operational Intensity:\s*([\d\\.eE+-]+)",
                    re.DOTALL
                )
                matches = re.findall(regex, stdout_output)

                if matches:
                    for func_name, fpc, total_input_size, oi in matches:
                        stripped_func_name = re.sub(r'\x1b\[[0-9;]*m', '', func_name).strip()
                        try:
                            fpc_val = float(fpc)
                            oi_val = float(oi)
                            total_input_sz = int(float(total_input_size))
                            # Determine size to save based on varyonly parameter
                            if varyonly == 'M':
                                size_to_save = m_val
                                # size_label = 'M'
                            elif varyonly == 'K':
                                size_to_save = k_val
                                # size_label = 'K'
                            elif varyonly == 'N':
                                size_to_save = n_val
                                # size_label = 'N'
                            else:
                                size_to_save = total_input_sz
                            size_label = 'total_input_size'
                            print(f"  {stripped_func_name} (Sparsity {s}): {size_label}={size_to_save}, {fpc_val} flops/cycle, {oi_val} flops/Byte")
                            case_result["results"][f"{stripped_func_name} (Sparsity {s})"] = {
                                size_label: size_to_save,
                                "operational_intensity": oi_val,
                                "performance": fpc_val
                            }
                            if correctness_status.get(stripped_func_name) == "failed":
                                print(f"    WARNING: {stripped_func_name} failed correctness check!")
                        except ValueError:
                            print(f"  Could not parse results for {stripped_func_name} at sparsity {s}")
                else:
                    print(f"  No performance results found for sparsity {s}.")
            except subprocess.TimeoutExpired:
                print(f"  ERROR: Benchmark run timed out for M={m_val}, K={k_val}, N={n_val}, s={s}")
            except Exception as e:
                print(f"  An unexpected error occurred for M={m_val}, K={k_val}, N={n_val}, s={s}: {e}")
        print("-" * len(header) + "\n")
        all_results.append(case_result)

    if save_results:
        output_file = outname
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=4)
        print(f"\nAll benchmark results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run SparseGEMM benchmarks')
    parser.add_argument('-s', '--save', action='store_true', help='Save benchmark results to JSON file')
    parser.add_argument('--output', type=str, default='benchmark_results.json', help='Name of JSON file')
    parser.add_argument('--varyonly', type=str, choices=['M','K','N'], help='Only vary one dimension: M, K, or N')
    args = parser.parse_args()
    run_and_parse_benchmark(args.save, args.output, args.varyonly)