import json
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker # For FixedLocator, FixedFormatter, NullLocator
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
    # print(ret)
    return ret

def format_size_to_human_readable_M(value_in_units):
    """
    Formats a numerical value into a string representation in millions (M),
    rounded to one decimal place.
    Example: 44509294 -> "44.5 M"
             500000   -> "0.5 M"
             44000000 -> "44.0 M"
    """
    if not isinstance(value_in_units, (int, float)):
        # Fallback for non-numeric types, though input should be numeric
        return str(value_in_units)
    
    millions = value_in_units #/ 1_000_000.0
    return f"{millions:.1f}"
    # return f"{int(millions)} M"

def create_performance_plot(json_filepath, title, outname):
    """
    Reads performance data from a JSON file, processes it,
    and generates a performance plot. X-axis ticks are set to all
    unique 'total_input_size' values and formatted in millions (e.g., "44.5 M").
    """
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
                input_size = float(input_size) # Keep as float for calculations
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
                'Input Size': round_size(input_size),
                'Function Name': func_name,
                'Performance (Flops/Cycle)': performance
            })

    if not plot_data:
        print("No data to plot. Check JSON structure, 'total_input_size', and 'performance'.")
        return

    df = pd.DataFrame(plot_data)

    plt.figure(figsize=(14, 8))
    ax = plt.gca()

    function_names = df['Function Name'].unique()
    for func_name in function_names:
        func_df = df[df['Function Name'] == func_name].sort_values(by='Input Size')
        if not func_df.empty:
            ax.plot(func_df['Input Size'], func_df['Performance (Flops/Cycle)'], marker='o', linestyle='-', label=func_name)

    ax.set_xlabel('Total Input Size (MB)')
    ax.set_ylabel('Performance (Flops/Cycle)')
    ax.set_title(title)
    ax.legend(title='Function Name', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True)
    ax.set_xscale('log')

    # --- X-axis tick annotation: Format labels as "XX.X M" ---
    unique_input_sizes = sorted(df['Input Size'].unique())
    
    if unique_input_sizes:
        tick_positions = unique_input_sizes
        # Use the new helper function to format labels
        tick_labels = [format_size_to_human_readable_M(s) for s in tick_positions]

        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(tick_positions))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(tick_labels))
        
        ax.tick_params(axis='x', which='major', labelrotation=0)
        plt.setp(ax.get_xticklabels(which='major'), ha='right')

        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        
        print(f"Info: Set X-axis ticks to ALL {len(unique_input_sizes)} unique input sizes, formatted as 'XX.X M'.")
        if len(unique_input_sizes) > 15: # Adjusted threshold for this format
            print("Warning: A large number of unique input sizes are being used as X-axis ticks.")
            print("         This may lead to a cluttered X-axis even with formatted labels.")
    else:
        print("Info: No unique input sizes found to set as X-axis ticks.")
    # --- End of X-axis tick annotation enhancement ---

    plt.tight_layout(rect=[0, 0, 0.85, 1])

    dt = datetime.now()
    ts = datetime.timestamp(dt) * 100
    ts = int(ts)
    output_dir = os.path.join("plots", "plot_images")
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, f"{outname}_{ts}.png")
    try:
        plt.savefig(plot_filename, dpi=300)
        print(f"Plot saved as {plot_filename}")
    except Exception as e:
        print(f"Error saving plot: {e}")

    # csv_filename = "performance_data.csv"
    # try:
    #     df.to_csv(csv_filename, index=False)
    #     print(f"Processed data saved as {csv_filename}")
    # except Exception as e:
    #     print(f"Error saving CSV data: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate performance plot with X-axis labels like '44.5 M'.")
    parser.add_argument("json_file", help="Path to input JSON file.")
    parser.add_argument("--title", default='Function Performance vs. Total Input Size', help="Path to input JSON file.")
    parser.add_argument("--outname", default='performance_plot')
    args = parser.parse_args()
    outname = args.outname
    if outname == 'performance_plot':
        to_app = args.json_file.split('/')[-1]
        to_app = to_app.split('.')[0]
        outname += '_' + to_app
    create_performance_plot(args.json_file, args.title, outname)