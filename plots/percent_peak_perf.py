

import json
import argparse

def compute_peak_percentages(data, peak):
    max_perf = {}
    for test in data:
        for func, stats in test['results'].items():
            perf = stats.get('performance', 0)
            if func not in max_perf or perf > max_perf[func]:
                max_perf[func] = perf

    percentages = {func: (perf / peak) * 100 for func, perf in max_perf.items()}
    return percentages

def main():
    parser = argparse.ArgumentParser(description='Compute peak performance percentages.')
    parser.add_argument('json_file', help='Path to JSON file with performance data.')
    parser.add_argument('--peak', type=float, default=4.0, help='Peak flops per cycle (default: 4).')
    args = parser.parse_args()

    with open(args.json_file, 'r') as f:
        data = json.load(f)

    percentages = compute_peak_percentages(data, args.peak)

    for func, pct in sorted(percentages.items()):
        print(f'{func}: {pct:.2f}%')

if __name__ == '__main__':
    main()