import json
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Compute average speedup from a JSON file for a specific K value")
    parser.add_argument("json_file", help="Path to the JSON file containing performance data")
    parser.add_argument("K", type=int, help="The K value to filter results by")
    args = parser.parse_args()
    k_to_select = args.K

    try:
        with open(args.json_file, "r") as f:
            data = json.load(f)
            # Filter entries to only those matching the specified K
            filtered_data = []
            for entry in data:
                tc = entry.get("test_case", {})
                if tc.get("K") == k_to_select:
                    filtered_data.append(entry)
            data = filtered_data
    except (IOError, json.JSONDecodeError) as e:
        print(f"Error reading JSON: {e}", file=sys.stderr)
        sys.exit(1)

    speedups = []
    for entry in data:
        results = entry.get("results", {})
        base_perf = None
        unrolled_perf = None
        for func_name, metrics in results.items():
            perf = metrics.get("performance")
            if perf is None:
                continue
            if "unrolled" in func_name.lower():
                unrolled_perf = perf
            else:
                base_perf = perf
        if base_perf is None or unrolled_perf is None:
            continue
        speedups.append(unrolled_perf / base_perf)

    if not speedups:
        print(f"No valid performance entries found for K = {k_to_select}.", file=sys.stderr)
        sys.exit(1)

    avg = sum(speedups) / len(speedups)
    print(f"Average speedup over all K: {avg:.2f}")

if __name__ == "__main__":
    main()