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
gflopsfontsize = 16
functionnamesize = 16
plt.rcParams['axes.linewidth'] = 0.8
import matplotlib.ticker
import pandas as pd
import os

def round_size(size):
    millions = size / 1_000_000.0
    ret = int(millions * 10) / 10
    dec = ret - int(ret)
    # print(dec)
    dec *= 10
    if dec < 3:
        return int(ret)
    elif dec > 7:
        return int(ret) + 1
    else:
        return int(ret) + 0.5
    return ret

def format_size_to_human_readable_M(value_in_units):
    if not isinstance(value_in_units, (int, float)):
        return str(value_in_units)
    
    millions = value_in_units / 1_000_000.0
    return f"{millions:.1f}"

def create_performance_plot(json_filepath, title, outname, xlabel, inline_labels=False):
    try:
        with open(json_filepath, 'r') as f:
            data = json.load(f)
        print(f"Successfully loaded data from {json_filepath}")
    except FileNotFoundError:
        print(f"Error: The file {json_filepath} was not found.")
        return
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from the file {json_filepath}. Please check the file format.")
        return
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        return

    plot_data = []
    for item_index, item in enumerate(data):
        test_case_info = item.get('test_case', {})
        results_dict = item.get('results', {})

        if not results_dict:
            print(f"Warning: No 'results' found in item {item_index}. Skipping.")
            continue

        for func_name, results_data in results_dict.items():
            input_size = results_data.get('total_input_size')
            performance = results_data.get('performance')

            if input_size is None:
                print(f"Warning: Missing 'total_input_size' for '{func_name}' in item {item_index}. Skipping.")
                continue
            try:
                input_size = float(input_size)
            except ValueError:
                print(f"Warning: Invalid 'total_input_size' for '{func_name}' (value: {input_size}) in item {item_index}. Skipping.")
                continue

            if performance is None:
                print(f"Warning: Missing 'performance' for '{func_name}' (input_size: {input_size}) in item {item_index}. Skipping.")
                continue
            try:
                performance = float(performance)
            except ValueError:
                print(f"Warning: Invalid 'performance' for '{func_name}' (value: {performance}) in item {item_index}. Skipping.")
                continue

            plot_data.append({
                'Input Size': input_size,
                'Function Name': func_name,
                'Performance (Flops/Cycle)': performance
            })

    if not plot_data:
        print("No data to plot. Check JSON structure, 'total_input_size', and 'performance'.")
        return

    df = pd.DataFrame(plot_data)

    plt.figure(figsize=(18, 8))
    ax = plt.gca()
    ax.set_facecolor('#dddddd')
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color='white')
    ax.xaxis.grid(False)
    for spine in ['left', 'top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)

    raw_function_names = df['Function Name'].unique()
    avg_perfs = {func: df[df['Function Name'] == func]['Performance (Flops/Cycle)'].mean() for func in raw_function_names}
    function_names = sorted(raw_function_names, key=lambda f: avg_perfs[f], reverse=True)
    palette = ['red', 'orange', 'purple', 'green', 'gray', 'blue']
    for idx, func_name in enumerate(function_names):
        func_df = df[df['Function Name'] == func_name].sort_values(by='Input Size')
        if func_df.empty:
            continue
        color = palette[idx % len(palette)]
        lw = 3 if idx == 0 else 1.8
        line, = ax.plot(func_df['Input Size'],
                        func_df['Performance (Flops/Cycle)'],
                        marker='o',
                        linestyle='-',
                        linewidth=lw,
                        color=color,
                        label=func_name)
        if inline_labels:
            # Annotate in the middle of the line
            mid_index = len(func_df) // 2
            mid_x = func_df['Input Size'].iloc[mid_index]
            mid_y = func_df['Performance (Flops/Cycle)'].iloc[mid_index]
            ax.annotate(
                func_name,
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
    ax.set_title(r"$\bf{" + escaped_title + r"}$" + "\n" + "(flops/cycle)", loc='left')
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
        
        print(f"Info: Set X-axis ticks to ALL {len(unique_input_sizes)} unique input sizes, formatted as 'XX.X M'.")
    else:
        print("Info: No unique input sizes found to set as X-axis ticks.")

    plt.tight_layout(rect=[0, 0, 0.8, 1])
    
    ax.set_ylim(bottom=0)
    dt = datetime.now()
    ts = datetime.timestamp(dt) * 100
    ts = int(ts)
    output_dir = os.path.join("plots", "plot_images", "performance")
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, f"{outname}.png")
    try:
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {plot_filename}")
    except Exception as e:
        print(f"Error saving plot: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate performance plot with X-axis labels like '44.5 M'.")
    parser.add_argument("json_file")
    parser.add_argument("--title", default='Performance of DFT 2ⁿ on Apple M2, 3.49 GHz')
    parser.add_argument("--outname", default='performance_plot')
    parser.add_argument("--xlabel", default='Total Input Size (MB)')
    parser.add_argument(
        "--inline-labels",
        action="store_true",
    )
    args = parser.parse_args()
    outname = args.outname
    if outname == 'performance_plot':
        to_app = args.json_file.split('/')[-1]
        to_app = to_app.split('.')[0]
        outname += '_' + to_app 
    create_performance_plot(
        args.json_file,
        args.title,
        outname,
        args.xlabel,
        args.inline_labels
    )