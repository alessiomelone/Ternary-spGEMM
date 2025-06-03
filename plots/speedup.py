import json
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Compute average speedup from a JSON file")
    parser.add_argument("json_file", help="Path to the JSON file containing performance data")
    args = parser.parse_args()

    try:
        with open(args.json_file, "r") as f:
            data = json.load(f)
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
        print("No valid performance entries found.", file=sys.stderr)
        sys.exit(1)

    avg = sum(speedups) / len(speedups)
    print(f"Average speedup over all K: {avg:.2f}")

if __name__ == "__main__":
    main()