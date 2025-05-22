import matplotlib.pyplot as plt
import numpy as np
import sys
import json
import csv
import argparse
import os # Added for directory and path operations
from collections import defaultdict

def plot_roofline(measured_op_intensities, measured_performances,
                  beta_memory_bw, pi_max_cpu_perf,
                  point_labels=None, plot_title="Roofline Plot",
                  oi_label="Operational Intensity (Flops/Byte)",
                  perf_label="Performance (Flops/Cycle)"):
    plot_points = True
    if not measured_op_intensities or not measured_performances or \
       len(measured_op_intensities) == 0 or len(measured_performances) == 0:
        print("Warning: No measured data points loaded. Plotting only ceilings.")
        _default_ridge_oi = pi_max_cpu_perf / beta_memory_bw if beta_memory_bw > 0 else 1.0
        # Create dummy OIs for oi_range calculation if no points
        measured_op_intensities = np.array([_default_ridge_oi / 10 if _default_ridge_oi > 0 else 0.1, 
                                            _default_ridge_oi * 10 if _default_ridge_oi > 0 else 10.0])
        measured_performances = np.array([pi_max_cpu_perf / 2, pi_max_cpu_perf / 2]) # Not plotted, for range calc
        plot_points = False


    if plot_points and (len(measured_op_intensities) != len(measured_performances)):
        # This check should ideally be handled by the CSV parser, but as a safeguard:
        raise ValueError("Operational intensity and performance lists must have the same length.")
    if plot_points and point_labels and len(point_labels) != len(measured_op_intensities):
        # This implies some labels were None, and the list was then filtered if not all None.
        # The load_data_from_csv ensures point_labels is None if it's not usable.
        print("Warning: Point labels length mismatch, disabling individual text labels.")
        point_labels = None # Fallback

    # Ensure these are numpy arrays for calculations
    measured_op_intensities = np.array(measured_op_intensities, dtype=float)
    measured_performances = np.array(measured_performances, dtype=float)


    ridge_point_oi = pi_max_cpu_perf / beta_memory_bw if beta_memory_bw > 0 else float('inf')
    
    min_oi_data = np.min(measured_op_intensities) if len(measured_op_intensities) > 0 else (ridge_point_oi / 10 if ridge_point_oi != float('inf') else 0.1)
    max_oi_data = np.max(measured_op_intensities) if len(measured_op_intensities) > 0 else (ridge_point_oi * 10 if ridge_point_oi != float('inf') else 10.0)
    
    _min_ois_to_consider = [min_oi_data / 5]
    if ridge_point_oi != float('inf') : _min_ois_to_consider.append(ridge_point_oi / 10)
    min_oi = min(_min_ois_to_consider + [0.01]) 

    _max_ois_to_consider = [max_oi_data * 5]
    if ridge_point_oi != float('inf') : _max_ois_to_consider.append(ridge_point_oi * 10)
    max_oi = max(_max_ois_to_consider + [min_oi * 1000 if min_oi > 0 else 100.0]) 

    if min_oi <= 0: min_oi = 0.001 
    if max_oi <= min_oi : max_oi = min_oi * 1000 

    oi_range = np.logspace(np.log10(min_oi), np.log10(max_oi), num=200)
    roofline_perf = np.minimum(pi_max_cpu_perf, beta_memory_bw * oi_range)
    unbounded_beta_perf = beta_memory_bw * oi_range

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor('#f0f0f0') # A light gray color

    ax.plot(oi_range, unbounded_beta_perf, linestyle=':', color='dimgray', lw=1.5, label=f"Memory Bandwidth (β={beta_memory_bw:.2f})")
    ax.axhline(y=pi_max_cpu_perf, color='gray', linestyle='--', lw=1.5, label=f"Peak Performance (π={pi_max_cpu_perf:.2f})")
    ax.plot(oi_range, roofline_perf, 'k-', lw=2.5, label='Attainable Performance Ceiling')

    if plot_points:
        # If point_labels is a list (even with Nones), it means there was a label column.
        # We will try to annotate points if labels exist for them.
        # For the legend, we create one entry if labels are used for annotation, or the generic one.
        has_any_actual_label = point_labels is not None and any(l is not None for l in point_labels)
        
        if has_any_actual_label:
            # Plot points individually to handle potentially None labels for legend
            for i in range(len(measured_op_intensities)):
                current_label = point_labels[i] if point_labels[i] else None # Legend label only if string
                ax.plot(measured_op_intensities[i], measured_performances[i], 'o', markersize=9, label=current_label)
                if point_labels[i] is not None: # Text annotation only if string
                    ax.text(measured_op_intensities[i] * 1.1, measured_performances[i], point_labels[i], fontsize=9, verticalalignment='center')
        else: # No valid labels at all, or point_labels was None from start
            ax.plot(measured_op_intensities, measured_performances, 'o', markersize=9, label='Measured Kernels/Points')


    if ridge_point_oi != float('inf') and ridge_point_oi <= max_oi and ridge_point_oi >= min_oi:
        ax.plot(ridge_point_oi, pi_max_cpu_perf, 'rX', markersize=12, label=f'Ridge Point ({ridge_point_oi:.2f} F/B)')
        ax.text(ridge_point_oi * 0.95, pi_max_cpu_perf * 0.90, f'({ridge_point_oi:.2f}, {pi_max_cpu_perf:.2f})',
                 color='red', fontsize=9, horizontalalignment='right', verticalalignment='top')

    # Annotation for Peak Performance (π)
    idx_pi_label = int(len(oi_range) * 0.85)
    if idx_pi_label >= len(oi_range): idx_pi_label = len(oi_range) - 1
    if idx_pi_label < 0: idx_pi_label = 0
    annot_pi_x_pos = oi_range[idx_pi_label]
    if not (np.isfinite(annot_pi_x_pos) and annot_pi_x_pos >= min_oi and annot_pi_x_pos <= max_oi) :
        annot_pi_x_pos = (min_oi + max_oi) * 0.8
    try: pi_units = perf_label.split('(')[-1].split(')')[0]
    except Exception: pi_units = "Perf Units"
    ax.text(annot_pi_x_pos, pi_max_cpu_perf * 1.08,
             f"Peak Performance (π) = {pi_max_cpu_perf:.2f} {pi_units}",
             color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=10)

    if beta_memory_bw > 0: # Only attempt to label if beta is meaningful
        # Choose an x-position on the slope, closer to the start of oi_range
        # (Using annot_beta_x_idx as it was in the user-provided script's context for this label)
        annot_beta_x_idx = max(5, int(len(oi_range) * 0.1)) 
        if annot_beta_x_idx >= len(oi_range): annot_beta_x_idx = len(oi_range) -1 
        
        annot_beta_x_pos = oi_range[annot_beta_x_idx]
        annot_beta_y_pos = beta_memory_bw * annot_beta_x_pos
        
        if ridge_point_oi != float('inf') and annot_beta_y_pos > pi_max_cpu_perf * 0.9 : 
            annot_beta_x_pos = min(ridge_point_oi * 0.5, oi_range[max(0, int(len(oi_range) * 0.2))]) 
            annot_beta_y_pos = beta_memory_bw * annot_beta_x_pos

        # Angle calculation for text rotation (as in user's provided script)
        # p1_display_for_angle = plt.gca().transData.transform_point((annot_beta_x_pos, annot_beta_y_pos))
        
        next_idx = min(annot_beta_x_idx + 1, len(oi_range) -1)
        if len(oi_range) > next_idx and oi_range[next_idx] == annot_beta_x_pos and next_idx + 1 < len(oi_range): 
          next_idx +=1
        
        annot_beta_x_pos_p2 = oi_range[next_idx] if len(oi_range) > next_idx else annot_beta_x_pos # handle very short oi_range
        annot_beta_y_pos_p2 = beta_memory_bw * annot_beta_x_pos_p2

        # if annot_beta_x_pos_p2 == annot_beta_x_pos : 
        #     angle_deg = 45 
        # else:
        #     p2_display_for_angle = plt.gca().transData.transform_point((annot_beta_x_pos_p2, annot_beta_y_pos_p2))
        #     # Using p1_display_for_angle and p2_display_for_angle for arctan2
        #     angle_rad = np.arctan2(p2_display_for_angle[1] - p1_display_for_angle[1], p2_display_for_angle[0] - p1_display_for_angle[0])
        #     angle_deg = np.degrees(angle_rad)
        #     if angle_deg > 75 : angle_deg = 75 
        #     if angle_deg < 0 : angle_deg = 0 
        
        dy = annot_beta_y_pos_p2 - annot_beta_y_pos
        dx = annot_beta_x_pos_p2 - annot_beta_x_pos
        angle = np.rad2deg(np.arctan2(dy, dx))

        try: 
            oi_byte_unit = oi_label.split('(')[-1].split('/')[1].split(')')[0]
            perf_cycle_unit = perf_label.split('(')[-1].split('/')[1].split(')')[0]
            beta_units = f"{oi_byte_unit}/{perf_cycle_unit}"
        except IndexError:
            beta_units = "BW Units"

        plt.text(annot_beta_x_pos, annot_beta_y_pos, # Original positioning from user script
                 f"Memory Bandwidth (β) = {beta_memory_bw:.2f} {beta_units}",
                 color='dimgray',
                 horizontalalignment='left', 
                 verticalalignment='bottom', 
                 fontsize=10, 
                 rotation=angle,
                 rotation_mode='anchor',
                 transform_rotates_text=True) # <<< KEY ADDITION HERE

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(oi_label, fontsize=14)
    ax.set_ylabel(perf_label, fontsize=14)
    ax.set_title(plot_title, fontsize=16, pad=20)

    all_x_points = list(measured_op_intensities if plot_points else []) + [min_oi, max_oi]
    if ridge_point_oi != float('inf'): all_x_points.append(ridge_point_oi)
    all_y_points = list(measured_performances if plot_points else []) + [pi_max_cpu_perf]
    roof_y_points_for_limits = np.minimum(pi_max_cpu_perf, beta_memory_bw * np.array([min_oi, max_oi, ridge_point_oi if ridge_point_oi != float('inf') and ridge_point_oi > min_oi and ridge_point_oi < max_oi else min_oi]))
    all_y_points.extend(list(roof_y_points_for_limits))
    all_y_points = [y for y in all_y_points if y > 0 and np.isfinite(y)]
    all_x_points = [x for x in all_x_points if x > 0 and np.isfinite(x)]

    if not all_x_points or not all_y_points: 
        final_min_oi, final_max_oi = 0.001, 100.0
        final_min_perf, final_max_perf = 0.001, 100.0
    else:
        final_min_perf = np.min(all_y_points) / 2.0
        final_max_perf = np.max(all_y_points) * 2.0
        final_min_oi = np.min(all_x_points) / 2.0
        final_max_oi = np.max(all_x_points) * 2.0
        if final_min_perf <= 0: final_min_perf = 0.0001 
        if final_max_perf <= final_min_perf: final_max_perf = final_min_perf * 1000
        if final_min_oi <= 0: final_min_oi = 0.0001
        if final_max_oi <= final_min_oi: final_max_oi = final_min_oi * 1000
    
    ax.set_ylim(bottom=final_min_perf, top=final_max_perf)
    ax.set_xlim(left=final_min_oi, right=final_max_oi)

    ax.legend(loc='best', fontsize=10, framealpha=0.7)
    ax.grid(True, which="both", ls="-", alpha=0.3)
    plt.tight_layout()
    plt.show()


def load_data_from_csv(filepath):
    """
    Loads roofline parameters and data points from a CSV file.
    """
    beta = None
    pi_cpu = None
    ois = []
    perfs = []
    labels = [] # Will store strings or None
    
    parsing_data_points = False
    header_map = {} # To map column names to indices

    # Define flexible header names
    oi_keys = {"oi", "operational_intensity", "operational intensity (flops/byte)", "intensity (flop/b)"}
    perf_keys = {"performance", "flops/cycle", "perf", "perf (flop/c)", "performance (flops/cycle)"}
    label_keys = {"label", "point_label", "name", "kernel", "kernel_name"}

    try:
        with open(filepath, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            for row_num, row in enumerate(reader):
                if not row or not row[0].strip() or row[0].strip().startswith('#'):
                    continue
                
                clean_row = [cell.strip() for cell in row]
                param_key_csv = clean_row[0].lower().replace('_', ' ') # Normalize key

                if not parsing_data_points:
                    if param_key_csv == "beta memory bw" or param_key_csv == "beta":
                        if len(clean_row) > 1: beta = float(clean_row[1])
                        else: raise ValueError(f"CSV format error: '{row[0]}' missing value on line {row_num+1}")
                    elif param_key_csv == "pi max cpu perf" or param_key_csv == "pi":
                        if len(clean_row) > 1: pi_cpu = float(clean_row[1])
                        else: raise ValueError(f"CSV format error: '{row[0]}' missing value on line {row_num+1}")
                    else:
                        # Check if this row is a potential data header
                        temp_header = [h.lower().strip() for h in clean_row]
                        is_oi_header_present = any(key in temp_header for key in oi_keys)
                        is_perf_header_present = any(key in temp_header for key in perf_keys)

                        if is_oi_header_present and is_perf_header_present:
                            for idx, h_name in enumerate(temp_header):
                                header_map[h_name] = idx
                            
                            # Find actual indices
                            for key_set, target_var_name in [(oi_keys, "oi_idx"), (perf_keys, "perf_idx"), (label_keys, "label_idx")]:
                                found = False
                                for key in key_set:
                                    if key in header_map:
                                        globals()[target_var_name] = header_map[key]
                                        found = True
                                        break
                                if not found and target_var_name != "label_idx": # OI and Perf are mandatory
                                    raise ValueError(f"CSV Data header on line {row_num+1} missing mandatory column for {target_var_name.replace('_idx','')} (e.g., oi, performance)")
                                elif not found and target_var_name == "label_idx":
                                    globals()[target_var_name] = -1 # Label is optional
                            
                            parsing_data_points = True
                else: # Parsing data points
                    if not header_map: # Should not happen if parsing_data_points is True
                        raise ValueError("CSV format error: Data points found before a valid header.")
                    
                    current_oi, current_perf, current_label = None, None, None
                    try:
                        current_oi = float(clean_row[globals()["oi_idx"]])
                        current_perf = float(clean_row[globals()["perf_idx"]])
                        if globals()["label_idx"] != -1 and len(clean_row) > globals()["label_idx"] and clean_row[globals()["label_idx"]]:
                            current_label = clean_row[globals()["label_idx"]]
                        
                        ois.append(current_oi)
                        perfs.append(current_perf)
                        labels.append(current_label) # Appends None if no label
                    except (IndexError, ValueError) as e:
                        print(f"Warning: Skipping data row {row_num+1} due to error ({e}): {row}")
                        continue
        
        if beta is None or pi_cpu is None:
            raise ValueError("CSV file must define 'beta_memory_bw' and 'pi_max_cpu_perf'.")
        
        if not parsing_data_points and (len(ois) == 0) : # No data header found or no data points
            print("Warning: No data points or data header found in CSV. Only ceilings will be plotted.")
            # plot_roofline handles empty ois/perfs

        # If all labels are None (e.g. no label column or all entries empty), pass None for point_labels
        if all(l is None for l in labels):
            labels = None

    except FileNotFoundError:
        print(f"Error: CSV file not found at {filepath}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error parsing CSV file '{filepath}': {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading CSV '{filepath}': {e}", file=sys.stderr)
        sys.exit(1)
        
    return beta, pi_cpu, ois, perfs, labels

def generate_roofline_csvs_from_json(input_json_path, output_directory, beta_bw, pi_perf, title):
    """
    Reads data from the input JSON file, processes it, and writes one
    CSV file per test case to the specified output directory.
    Each CSV is formatted for the roofline plotting script.

    Args:
        input_json_path (str): Path to the input JSON file.
        output_directory (str): Path to the directory where output CSV files will be saved.
        beta_bw (float): Beta (memory bandwidth) value.
        pi_perf (float): Pi (max CPU performance) value.
    """
    try:
        with open(input_json_path, 'r', encoding='utf-8') as f_json:
            data = json.load(f_json)
    except FileNotFoundError:
        print(f"Error: Input JSON file not found at '{input_json_path}'", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{input_json_path}'. Please check its format.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading JSON: {e}", file=sys.stderr)
        sys.exit(1)

    if not isinstance(data, list):
        print("Error: Expected JSON data to be a list of test cases.", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    try:
        os.makedirs(output_directory, exist_ok=True)
    except OSError as e:
        print(f"Error: Could not create output directory '{output_directory}': {e}", file=sys.stderr)
        sys.exit(1)

    # Process each item (test_case block) in the JSON data
    for item_index, item in enumerate(data):
        if not isinstance(item, dict) or "test_case" not in item or "results" not in item:
            print(f"Warning: Skipping invalid item (index {item_index}) in JSON: {item}", file=sys.stderr)
            continue
        
        test_case = item.get("test_case", {})
        results = item.get("results", {})

        try:
            m = test_case.get("M")
            k = test_case.get("K")
            n = test_case.get("N")
            s_val = test_case.get("s") # Now using 's' for filename

            if None in [m, k, n, s_val]: # Ensure all needed components for filename and OI are present
                print(f"Warning: Skipping item (index {item_index}) due to missing M, K, N, or s in test_case: {test_case}", file=sys.stderr)
                continue
            
            # Cast to int for filename and OI calculation
            m, k, n, s_val = int(m), int(k), int(n), int(s_val)

            # Construct filename for this test case
            output_csv_filename = f"roofline_M{m}_K{k}_N{n}_s{s_val}.csv"
            title = f"Roofline plot for test case M:{m}, K:{k}, N:{n} - s:{s_val}"
            output_csv_filename = f'temp.csv'
            output_csv_filepath = os.path.join(output_directory, output_csv_filename)

            with open(output_csv_filepath, 'w', newline='', encoding='utf-8') as f_csv:
                writer = csv.writer(f_csv)

                # Write machine parameters
                writer.writerow(["parameter", "value"])
                writer.writerow(["beta_memory_bw", beta_bw])
                writer.writerow(["pi_max_cpu_perf", pi_perf])
                writer.writerow([]) # Empty line for separation
                writer.writerow([f"# Kernel Performance Data for M={m}, K={k}, N={n}, s={s_val}"])

                # Write data header
                writer.writerow(["oi", "performance", "label"])

                # Calculate "oi" for this test case (same for all results in this file)
                memory_CHANGE_ME = (m * k + k * n + n) * 4 # Bytes

                if not results:
                    print(f"Info: No 'results' entries for test case M={m},K={k},N={n},s={s_val}. CSV file '{output_csv_filepath}' will contain parameters only.", file=sys.stdout)

                for result_key, res_dict in results.items():
                    # if performance_value is None:
                    #     print(f"Warning: Skipping result '{result_key}' for M={m},K={k},N={n},s={s_val} due to null performance value.", file=sys.stderr)
                    #     continue
                    
                    try:
                        cycles_val_float = float(res_dict['cycles'])
                        perf_val_float = float(res_dict['fpc'])
                    except (ValueError, TypeError):
                        print(f"Warning: Skipping result '{result_key}' for M={m},K={k},N={n},s={s_val} due to non-numeric performance value: '{res_dict['fpc']}'.", file=sys.stderr)
                        continue

                    oi_value = (cycles_val_float * perf_val_float) / (memory_CHANGE_ME / 1e9) # Cycles / GB

                    label = f"{result_key}_{m}x{k}x{n}"
                    writer.writerow([oi_value, perf_val_float, label])
            
            # print(f"Successfully created '{output_csv_filepath}'")

            beta, pi_cpu, ois, perfs, point_labels = load_data_from_csv(output_csv_filepath)
            os.remove(output_csv_filepath)

            plot_roofline(ois, perfs, 
                        beta, pi_cpu, 
                        point_labels=point_labels, 
                        plot_title=title,
                        oi_label='Operational Intensity (Flops/GB)',
                        perf_label='Performance (Flops/Cycle)')

        except (TypeError, ValueError) as e:
            print(f"Warning: Skipping item (index {item_index}) due to data type error for M,K,N,s or results: {e} (test_case: {test_case})", file=sys.stderr)
            continue
        except IOError:
            print(f"Error: Could not write to output CSV file '{output_csv_filepath}'", file=sys.stderr)
            # Continue to next item if one file fails
        except Exception as e:
            print(f"An unexpected error occurred while processing item (index {item_index}): {e}", file=sys.stderr)
            # Continue to next item


def aggregate_and_write_prefix_csvs(input_json_path, output_directory, beta_bw, pi_perf):
    """
    Reads data from the input JSON, aggregates results by algorithm prefix across all test cases,
    and writes one CSV file per prefix to the specified output directory.
    Each CSV is formatted for the roofline plotting script.

    Args:
        input_json_path (str): Path to the input JSON file.
        output_directory (str): Path to the directory where output CSV files will be saved.
        beta_bw (float): Beta (memory bandwidth) value.
        pi_perf (float): Pi (max CPU performance) value.
    """
    try:
        with open(input_json_path, 'r', encoding='utf-8') as f_json:
            all_test_data = json.load(f_json)
    except FileNotFoundError:
        print(f"Error: Input JSON file not found at '{input_json_path}'", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{input_json_path}'. Please check its format.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading JSON: {e}", file=sys.stderr)
        sys.exit(1)

    if not isinstance(all_test_data, list):
        print("Error: Expected JSON data to be a list of test cases.", file=sys.stderr)
        sys.exit(1)

    try:
        os.makedirs(output_directory, exist_ok=True)
    except OSError as e:
        print(f"Error: Could not create output directory '{output_directory}': {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 1: Aggregate data by prefix ---
    aggregated_data_by_prefix = defaultdict(list)

    for item_index, item in enumerate(all_test_data):
        if not isinstance(item, dict) or "test_case" not in item or "results" not in item:
            print(f"Warning: Skipping invalid item (index {item_index}) in JSON: {item}", file=sys.stderr)
            continue
        
        test_case = item.get("test_case", {})
        results_dict = item.get("results", {})

        try:
            m_val = test_case.get("M")
            k_val = test_case.get("K")
            n_val = test_case.get("N")
            s_val = test_case.get("s")

            if None in [m_val, k_val, n_val, s_val]:
                print(f"Warning: Skipping item (index {item_index}) due to missing M, K, N, or s in test_case: {test_case}", file=sys.stderr)
                continue
            
            m, k, n, s = int(m_val), int(k_val), int(n_val), int(s_val)
            
            memory_CHANGE_ME = (m * k + k * n + n) * 4
            product_m_n_k = m * n * k

            for original_key, result_data in results_dict.items():
                if not isinstance(result_data, dict) or "fpc" not in result_data:
                    print(f"Warning: Skipping result '{original_key}' for M{m}K{k}N{n}s{s} due to missing 'fpc' or invalid format.", file=sys.stderr)
                    continue
                
                fpc = result_data.get("fpc")
                if fpc is None:
                    print(f"Warning: Skipping result '{original_key}' for M{m}K{k}N{n}s{s} due to null 'fpc' value.", file=sys.stderr)
                    continue
                
                try:
                    fpc_float = float(fpc)
                except (ValueError, TypeError):
                    print(f"Warning: Skipping result '{original_key}' for M{m}K{k}N{n}s{s} due to non-numeric 'fpc': '{fpc}'.", file=sys.stderr)
                    continue

                prefix = original_key.split('_')[0]
                # New label format: <original_key>(prod:<M*N*K>_s:<s>)
                label = f"{original_key}(prod:{product_m_n_k}_s:{s})"
                
                aggregated_data_by_prefix[prefix].append({
                    "oi": memory_CHANGE_ME,
                    "performance": fpc_float,
                    "label": label
                })

        except (TypeError, ValueError) as e:
            print(f"Warning: Skipping item (index {item_index}) due to data type error for M,K,N,s: {e} (test_case: {test_case})", file=sys.stderr)
            continue
        except Exception as e_item:
             print(f"An unexpected error occurred processing item (index {item_index}): {e_item}", file=sys.stderr)

    # --- Step 2: Write one CSV per prefix ---
    if not aggregated_data_by_prefix:
        print("No valid data was aggregated from the JSON file. No CSV files will be created.", file=sys.stdout)
        return

    for prefix, data_rows_list in aggregated_data_by_prefix.items():
        # output_csv_filename = f"roofline_{prefix}.csv"
        # output_csv_filepath = os.path.join(output_directory, output_csv_filename)
        output_csv_filepath = './temp.csv'

        try:
            with open(output_csv_filepath, 'w', newline='', encoding='utf-8') as f_csv:
                writer = csv.writer(f_csv)
                # Write machine parameters
                writer.writerow(["parameter", "value"])
                writer.writerow(["beta_memory_bw", beta_bw])
                writer.writerow(["pi_max_cpu_perf", pi_perf])
                writer.writerow([]) 
                writer.writerow([f"# Kernel Performance Data for Prefix: {prefix} (aggregated across all test cases)"])
                
                # Write data header
                writer.writerow(["oi", "performance", "label"])
                
                for row_data in data_rows_list:
                    writer.writerow([row_data["oi"], row_data["performance"], row_data["label"]])

            beta, pi_cpu, ois, perfs, point_labels = load_data_from_csv(output_csv_filepath)
            os.remove(output_csv_filepath)

            plot_roofline(ois, perfs, 
                        beta, pi_cpu, 
                        point_labels=point_labels, 
                        plot_title=f'Roofline plot for: {prefix}',
                        oi_label='Operational Intensity (Flops/GB)',
                        perf_label='Performance (Flops/Cycle)')
            
            # print(f"Successfully created '{output_csv_filepath}' for prefix '{prefix}'")
        except IOError:
            print(f"Error: Could not write to output CSV file '{output_csv_filepath}'", file=sys.stderr)
        except Exception as e_csv:
            print(f"An unexpected error occurred while writing CSV '{output_csv_filepath}': {e_csv}", file=sys.stderr)


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(
#         description="Convert JSON results to multiple CSV files (one per algorithm prefix, aggregated across test cases) for roofline plotting.",
#         formatter_class=argparse.RawTextHelpFormatter
#     )
#     parser.add_argument("input_json", help="Path to the input JSON file.")
#     parser.add_argument("output_dir", help="Path to the directory where output CSV files will be saved.")
#     parser.add_argument("--beta", type=float, required=True, help="Beta (memory bandwidth) value (e.g., Bytes/Cycle).")
#     parser.add_argument("--pi", type=float, required=True, help="Pi (max CPU performance) value (e.g., Flops/Cycle).")

#     args = parser.parse_args()

#     aggregate_and_write_prefix_csvs(args.input_json, args.output_dir, args.beta, args.pi)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a Roofline plot from data in a CSV file.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""\
CSV File Format Example:
--------------------------
# Machine Parameters (keys are case-insensitive, space or underscore)
beta_memory_bw, 1.0
pi_max_cpu_perf, 2.0

# Optional comments can be used

# Kernel Performance Data
# Columns header: Must include OI and Performance. Label is optional.
# Recognized OI headers: oi, operational_intensity, intensity (flop/b)
# Recognized Perf headers: performance, flops/cycle, perf, perf (flop/c)
# Recognized Label headers: label, point_label, name, kernel
oi,performance,label
0.1,0.08,Kernel_A
0.5,0.4,Kernel_B
1.0,0.7
2.0,1.2,Kernel_D
"""
    )
    parser.add_argument("json_filepath", help="Path to the JSON file containing roofline data.")
    # parser.add_argument("--title", help="Optional title for the plot.", default="Roofline Plot")
    # parser.add_argument("--oi_label", help="Optional label for the OI axis.", default="Operational Intensity (Flops/GB)")
    # parser.add_argument("--perf_label", help="Optional label for the Performance axis.", default="Performance (Flops/Cycle)")
    parser.add_argument("--beta", type=float, default=100.0, help="Beta (memory bandwidth) value (e.g., Bytes/Cycle).")
    parser.add_argument("--pi", type=float, default=4.0, help="Pi (max CPU performance) value (e.g., Flops/Cycle).")
    parser.add_argument("--per_test", action='store_true', help="Group results by test.")
    args = parser.parse_args()

    if args.per_test:
        generate_roofline_csvs_from_json(args.json_filepath, './', args.beta, args.pi, '')
    else:
        aggregate_and_write_prefix_csvs(args.json_filepath, './', args.beta, args.pi)