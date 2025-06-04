#!/usr/bin/env python3

import json
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Calibri', 'Helvetica', 'Gill Sans MT', 'Arial']
plt.rcParams['axes.titlesize'] = 24
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
Xaxis_fontsize = 16
functionnamesize = 16
plt.rcParams['axes.linewidth'] = 0.8
import matplotlib.ticker
import pandas as pd
import os

def round_size(size):
    millions = size / 1_000_000.0
    ret = int(millions * 10) / 10
    dec = ret - int(ret)
    dec *= 10
    if dec < 3:
        return int(ret)
    elif dec > 7:
        return int(ret) + 1
    else:
        return int(ret) + 0.5

def format_size_to_human_readable_M(value_in_units):
    if not isinstance(value_in_units, (int, float)):
        return str(value_in_units)
    millions = value_in_units / 1_000_000.0
    return f"{millions:.1f}"

def create_ratio_plot(json_filepaths, title, outname, xlabel, inline_labels=False):
    """
    Reads multiple JSON files, each containing a baseline and an unrolled version of a function.
    Computes, for each input size, the ratio: (unrolled_performance / baseline_performance).
    Plots those ratios vs input size on a single log‐scale X axis, with styling matching the original.
    """
    all_plot_data = []
    for json_filepath in json_filepaths:
        try:
            with open(json_filepath, 'r') as f:
                data = json.load(f)
        except FileNotFoundError:
            print(f"Error: The file {json_filepath} was not found. Skipping.")
            continue
        except json.JSONDecodeError:
            print(f"Error: Could not decode JSON from {json_filepath}. Skipping.")
            continue
        except Exception as e:
            print(f"Unexpected error reading {json_filepath}: {e}. Skipping.")
            continue

        # Determine baseline vs unrolled key names
        # Assume JSON 'results' dict contains exactly two functions:
        # one whose name contains "Unrolled" (case‐insensitive), and one baseline.
        # If there are exactly two keys, assign accordingly; otherwise warn.
        sample_item = None
        for item in data:
            if item.get('results'):
                sample_item = item
                break
        if not sample_item:
            print(f"Warning: No valid 'results' entries in {json_filepath}. Skipping.")
            continue

        func_names = list(sample_item['results'].keys())
        if len(func_names) != 2:
            print(f"Warning: Expected exactly 2 functions in {json_filepath} but found {len(func_names)}. Attempting best guess.")
        # Identify unrolled by "unrolled" substring
        baseline_name = None
        unrolled_name = None
        for name in func_names:
            if 'unroll' in name.lower():
                unrolled_name = name
            else:
                baseline_name = name
        # If pattern matching failed, default to first two
        if not baseline_name or not unrolled_name:
            if len(func_names) >= 2:
                baseline_name, unrolled_name = func_names[0], func_names[1]
            else:
                print(f"Warning: Couldn't identify baseline/unrolled in {json_filepath}. Skipping.")
                continue

        # For labeling the line, use the JSON filename (without extension)
        filename_label = os.path.splitext(os.path.basename(json_filepath))[0]

        # Collect (input_size, ratio) pairs
        for item_index, item in enumerate(data):
            results = item.get('results', {})
            if baseline_name not in results or unrolled_name not in results:
                print(f"Warning: Missing '{baseline_name}' or '{unrolled_name}' in item {item_index} of {json_filepath}. Skipping that item.")
                continue

            base_entry = results[baseline_name]
            unrolled_entry = results[unrolled_name]

            # Extract and validate input size
            input_size = base_entry.get('total_input_size')
            if input_size is None:
                print(f"Warning: Missing 'total_input_size' for '{baseline_name}' in item {item_index} of {json_filepath}. Skipping.")
                continue
            try:
                input_size = float(input_size)
            except ValueError:
                print(f"Warning: Invalid 'total_input_size' ({input_size}) for '{baseline_name}' in item {item_index} of {json_filepath}. Skipping.")
                continue

            # Extract and validate baseline performance
            base_perf = base_entry.get('performance')
            if base_perf is None:
                print(f"Warning: Missing 'performance' for '{baseline_name}' at input_size {input_size} in {json_filepath}. Skipping.")
                continue
            try:
                base_perf = float(base_perf)
            except ValueError:
                print(f"Warning: Invalid 'performance' ({base_perf}) for '{baseline_name}' at input_size {input_size} in {json_filepath}. Skipping.")
                continue
            if base_perf == 0:
                print(f"Warning: Zero baseline performance for '{baseline_name}' at input_size {input_size} in {json_filepath}. Skipping.")
                continue

            # Extract and validate unrolled performance
            unrolled_perf = unrolled_entry.get('performance')
            if unrolled_perf is None:
                print(f"Warning: Missing 'performance' for '{unrolled_name}' at input_size {input_size} in {json_filepath}. Skipping.")
                continue
            try:
                unrolled_perf = float(unrolled_perf)
            except ValueError:
                print(f"Warning: Invalid 'performance' ({unrolled_perf}) for '{unrolled_name}' at input_size {input_size} in {json_filepath}. Skipping.")
                continue

            ratio = unrolled_perf / base_perf
            all_plot_data.append({
                'Input Size': input_size,
                'Performance Ratio': ratio,
                'Label': filename_label
            })

    if not all_plot_data:
        print("No valid data across all JSONs. Exiting without plotting.")
        return

    df = pd.DataFrame(all_plot_data)

    plt.figure(figsize=(18, 8))
    ax = plt.gca()
    ax.set_facecolor('#dddddd')
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color='white')
    ax.xaxis.grid(False)
    for spine in ['left', 'top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)

    raw_labels = df['Label'].unique()
    avg_ratios = {lbl: df[df['Label'] == lbl]['Performance Ratio'].mean() for lbl in raw_labels}
    labels_sorted = sorted(raw_labels, key=lambda l: avg_ratios[l], reverse=True)
    palette = ['red', 'orange', 'purple', 'green', 'gray', 'blue']
    for idx, lbl in enumerate(labels_sorted):
        subdf = df[df['Label'] == lbl].sort_values(by='Input Size')
        if subdf.empty:
            continue
        color = palette[idx % len(palette)]
        lw = 3 if idx == 0 else 1.8
        line, = ax.plot(subdf['Input Size'],
                        subdf['Performance Ratio'],
                        marker='o',
                        linestyle='-',
                        linewidth=lw,
                        color=color,
                        label=lbl)
        if inline_labels:
            mid_index = len(subdf) // 2
            mid_x = subdf['Input Size'].iloc[mid_index]
            mid_y = subdf['Performance Ratio'].iloc[mid_index]
            ax.annotate(
                lbl,
                xy=(mid_x, mid_y),
                xytext=(5, 5),
                textcoords='offset points',
                va='bottom',
                ha='left',
                fontsize=functionnamesize,
                fontweight='bold',
                color=color,
                zorder=10
            )
    if not inline_labels:
        ax.legend(fontsize=functionnamesize, loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_xlabel(xlabel, fontsize=Xaxis_fontsize)
    ax.yaxis.set_label_coords(0, 1.02)
    escaped_title = title.replace(" ", r"\ ")
    ax.set_title(r"$\bf{" + escaped_title + r"}$" + "\n" + "(ratio over baseline)", loc='left')
    ax.grid(axis='y', color='white')
    ax.set_xscale('log')

    unique_input_sizes = sorted(df['Input Size'].unique())
    if unique_input_sizes:
        tick_positions = unique_input_sizes
        tick_labels = [str(int(s)) for s in tick_positions]
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(tick_positions))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(tick_labels))
        ax.tick_params(axis='x', which='major', labelrotation=0)
        plt.setp(ax.get_xticklabels(which='major'), ha='right')
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        print(f"Info: Set X-axis ticks to all {len(unique_input_sizes)} unique input sizes.")
    else:
        print("Info: No unique input sizes found to set as X-axis ticks.")

    plt.tight_layout(rect=[0, 0, 0.8, 1])
    ax.set_ylim(bottom=0)

    # Save output
    dt = datetime.now()
    ts = int(datetime.timestamp(dt) * 100)
    output_dir = os.path.join("plots", "plot_images", "performance_ratio")
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, f"{outname}.png")
    try:
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {plot_filename}")
    except Exception as e:
        print(f"Error saving plot: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate a plot of performance increase (ratio over baseline) vs input size for multiple JSON files."
    )
    parser.add_argument(
        "json_files",
        nargs='+',
        help="One or more JSON files, each containing a baseline function and its unrolled version."
    )
    parser.add_argument("--title", default="Unrolling Speedup per Foramt on Apple M2, 3.49 GHz")
    parser.add_argument("--outname", default="performance_ratio_plot")
    parser.add_argument("--xlabel", default="K")
    parser.add_argument(
        "--inline-labels",
        action="store_true",
        help="If set, annotate each line with its label at midpoint."
    )
    args = parser.parse_args()

    # If outname is default and multiple JSONs provided, append joined filenames
    outname = args.outname
    if outname == "performance_ratio_plot":
        basename_list = [
            os.path.splitext(os.path.basename(fp))[0] for fp in args.json_files
        ]
        joined = "_".join(basename_list)
        outname = outname + "_" + joined

    create_ratio_plot(
        args.json_files,
        args.title,
        outname,
        args.xlabel,
        args.inline_labels
    )