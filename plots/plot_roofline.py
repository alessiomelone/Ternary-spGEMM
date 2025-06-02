import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import sys
import json
import csv
import argparse
import os
from collections import defaultdict

# --- plot_roofline_for_prefix function (MODIFIED for annotation position and background) ---
def plot_roofline_for_prefix(
    algorithms_data, # Dict: { "FullName1": [{"oi":val, "perf":val, "annotation_label":str}, ...], ... }
    beta_memory_bw,
    pi_max_cpu_perf,
    plot_title="Roofline Plot",
    oi_label="Operational Intensity (Flops/Byte)",
    perf_label="Performance (Flops/Cycle)"
):
    """
    Generates a roofline plot for a specific prefix, showing multiple algorithm variants.
    Point annotations are 'total_input_size', positioned on top with transparent background.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor('#f0f0f0')

    all_ois_for_plot_range = []
    if algorithms_data:
        for algo_name, points in algorithms_data.items():
            if points:
                all_ois_for_plot_range.extend([p['oi'] for p in points])
    
    if not all_ois_for_plot_range:
        _default_ridge_oi = pi_max_cpu_perf / beta_memory_bw if beta_memory_bw > 0 else 1.0
        min_oi = _default_ridge_oi / 10 if _default_ridge_oi > 0 else 0.01
        max_oi = _default_ridge_oi * 10 if _default_ridge_oi > 0 else 10.0
        if min_oi <= 0: min_oi = 0.001
        if max_oi <= min_oi : max_oi = min_oi * 1000
    else:
        min_oi_data = np.min(all_ois_for_plot_range)
        max_oi_data = np.max(all_ois_for_plot_range)
        ridge_point_oi_temp = pi_max_cpu_perf / beta_memory_bw if beta_memory_bw > 0 else float('inf')
        _min_ois_to_consider = [min_oi_data / 5]
        if ridge_point_oi_temp != float('inf') : _min_ois_to_consider.append(ridge_point_oi_temp / 10)
        min_oi = min(_min_ois_to_consider + [0.01])
        _max_ois_to_consider = [max_oi_data * 5]
        if ridge_point_oi_temp != float('inf') : _max_ois_to_consider.append(ridge_point_oi_temp * 10)
        max_oi = max(_max_ois_to_consider + [min_oi * 1000 if min_oi > 0 else 100.0])
        if min_oi <= 0: min_oi = 0.001
        if max_oi <= min_oi : max_oi = min_oi * 1000

    oi_range = np.logspace(np.log10(min_oi), np.log10(max_oi), num=200)
    roofline_perf_ceiling = np.minimum(pi_max_cpu_perf, beta_memory_bw * oi_range)
    unbounded_beta_perf_line = beta_memory_bw * oi_range

    ax.plot(oi_range, unbounded_beta_perf_line, linestyle=':', color='dimgray', lw=1.5, label=f"Unbounded Mem BW (β={beta_memory_bw:.2f})")
    ax.axhline(y=pi_max_cpu_perf, color='gray', linestyle='--', lw=1.5, label=f"Peak Perf (π={pi_max_cpu_perf:.2f})")
    ax.plot(oi_range, roofline_perf_ceiling, 'k-', lw=2.5, label='Attainable Performance Ceiling')

    if algorithms_data:
        for full_algorithm_name, points_list in algorithms_data.items():
            if not points_list:
                continue
            current_ois = np.array([p['oi'] for p in points_list], dtype=float)
            current_performances = np.array([p['performance'] for p in points_list], dtype=float)
            current_annotations = [p['annotation_label'] for p in points_list] 
            current_annotations_int = [int(input_size) for input_size in current_annotations]
            current_annotations_int = sorted(current_annotations_int)
            min_annot_val = min(current_annotations_int)
            max_annot_val = max(current_annotations_int)

            if len(current_ois) > 0:
                sorted_indices = np.argsort(current_ois)
                sorted_ois = current_ois[sorted_indices]
                sorted_perfs = current_performances[sorted_indices]
                sorted_annotations = [current_annotations[i] for i in sorted_indices]
                
                ax.plot(sorted_ois, sorted_perfs, 
                        marker='o', linestyle='-', markersize=7, label=full_algorithm_name) 

                for i in range(len(sorted_ois)):
                    if sorted_annotations[i] is not None:
                        if int(sorted_annotations[i]) == min_annot_val or int(sorted_annotations[i]) == max_annot_val:
                            try: # Format annotation for readability
                                annot_val = float(sorted_annotations[i])
                                if annot_val >= 1e6: annot_text = f"{annot_val/1e6:.2f}MB"
                                elif annot_val >= 1e3: annot_text = f"{annot_val/1e3:.2f}KB"
                                elif isinstance(annot_val, float) and annot_val.is_integer(): annot_text = str(int(annot_val))
                                else: annot_text = f"{int(annot_val)}" if isinstance(annot_val, float) else str(sorted_annotations[i])
                            except ValueError: annot_text = str(sorted_annotations[i]) 
                        else:
                            annot_text = ''

                        # # Position text strictly on top
                        # y_text_position = sorted_perfs[i] * 1.05 # Slight nudge up on log scale
                        # if y_text_position > np.max(ax.get_ylim())*0.95 : # Prevent going too high off screen
                        #     y_text_position = sorted_perfs[i] * 0.95
                        #     v_align = 'top'
                        # else:
                        #     v_align = 'bottom'

                        y_text_position = sorted_perfs[i] * 1.05 # Slight nudge up on log scale
                        v_align = 'bottom'
                        # Basic check to prevent label going too far off-screen if point is at the very top
                        current_ylim = ax.get_ylim() # Get current limits for this check
                        if y_text_position > current_ylim[1]*0.85 : # If getting too close to top limit
                            y_text_position = sorted_perfs[i] * 0.95 # Place below
                            v_align = 'top'
                        elif y_text_position < current_ylim[0]*1.15 and v_align == 'bottom': # If too close to bottom limit
                            y_text_position = sorted_perfs[i] * 1.05 # Ensure it stays above
                            v_align = 'bottom'

                        if int(sorted_annotations[i]) == min_annot_val:
                            v_align = 'top'
                            horalign='left'
                        elif int(sorted_annotations[i]) == max_annot_val:
                            v_align = 'bottom'
                            horalign='right'
                        else:
                            horalign='center'
                        ax.text(sorted_ois[i], # X: Same as point
                                y_text_position, 
                                annot_text, 
                                fontsize=6, 
                                color='black', # Explicitly black text
                                horizontalalignment=horalign, # Center text horizontally over the point
                                verticalalignment=v_align, # Bottom of text starts at y_text_position
                                # No bbox argument for transparent background
                                )
    
    ridge_point_oi = pi_max_cpu_perf / beta_memory_bw if beta_memory_bw > 0 else float('inf')
    if ridge_point_oi != float('inf') and ridge_point_oi <= max_oi and ridge_point_oi >= min_oi: # Show ridge point
        ax.plot(ridge_point_oi, pi_max_cpu_perf, 'rX', markersize=10, label=f'Ridge Pt ({ridge_point_oi:.2f} F/B)')
        ax.text(ridge_point_oi * 0.95, pi_max_cpu_perf * 0.90, f'({ridge_point_oi:.2f}, {pi_max_cpu_perf:.2f})',
                 color='red', fontsize=9, horizontalalignment='right', verticalalignment='top')

    idx_pi_label = int(len(oi_range) * 0.85)
    annot_pi_x_pos = oi_range[min(idx_pi_label, len(oi_range)-1)] 
    if not (np.isfinite(annot_pi_x_pos) and annot_pi_x_pos >= min_oi and annot_pi_x_pos <= max_oi) :
        annot_pi_x_pos = (min_oi + max_oi) * 0.8
    try: pi_units = perf_label.split('(')[-1].split(')')[0]
    except Exception: pi_units = "Perf Units"
    ax.text(annot_pi_x_pos, pi_max_cpu_perf * 1.08,
             f"Peak Perf (π) = {pi_max_cpu_perf:.2f} {pi_units}",
             color='black', horizontalalignment='center', verticalalignment='bottom', fontsize=10)

    if beta_memory_bw > 0: 
        if ridge_point_oi != float('inf') and ridge_point_oi > min_oi * 1.2:
            anchor_oi_beta = np.sqrt(min_oi * ridge_point_oi)
            anchor_oi_beta = max(min_oi * 1.2, anchor_oi_beta)
            anchor_oi_beta = min(ridge_point_oi * 0.6, anchor_oi_beta)
        else:
            anchor_oi_beta = oi_range[max(0,int(len(oi_range)*0.1))]
        anchor_oi_beta = max(min_oi, anchor_oi_beta)
        anchor_perf_beta = beta_memory_bw * anchor_oi_beta
        if anchor_perf_beta > pi_max_cpu_perf * 0.8 and pi_max_cpu_perf > 0 : 
             anchor_perf_beta = pi_max_cpu_perf * 0.5
             anchor_oi_beta = anchor_perf_beta / beta_memory_bw if beta_memory_bw > 0 else min_oi
             anchor_oi_beta = max(min_oi, anchor_oi_beta)

        idx1_angle = max(0, int(len(oi_range) * 0.1))
        idx2_angle = min(len(oi_range) - 1, int(len(oi_range) * 0.3))
        if idx1_angle >= idx2_angle:
            if len(oi_range) > 1: idx1_angle, idx2_angle = 0, 1
            else: idx1_angle, idx2_angle = 0, 0 
        oi_p1_angle, oi_p2_angle = oi_range[idx1_angle], oi_range[idx2_angle]
        p1_data_angle = (oi_p1_angle, beta_memory_bw * oi_p1_angle)
        p2_data_angle = (oi_p2_angle, beta_memory_bw * oi_p2_angle)
        p1_display_angle = ax.transData.transform_point(p1_data_angle)
        p2_display_angle = ax.transData.transform_point(p2_data_angle)
        if np.all(np.isclose(p1_display_angle, p2_display_angle)):
            angle_rad = np.deg2rad(45)
        else:
            angle_rad = np.arctan2(p2_display_angle[1] - p1_display_angle[1], p2_display_angle[0] - p1_display_angle[0])
        angle_deg = np.degrees(angle_rad)

        offset_points = 7 
        perp_angle_rad = angle_rad + np.pi / 2 
        dx_display = offset_points * np.cos(perp_angle_rad)
        dy_display = offset_points * np.sin(perp_angle_rad)
        anchor_beta_display = ax.transData.transform_point((anchor_oi_beta, anchor_perf_beta))
        text_pos_beta_display = (anchor_beta_display[0] + dx_display, anchor_beta_display[1] + dy_display)
        text_pos_beta_data = ax.transData.inverted().transform_point(text_pos_beta_display)
        try:
            oi_den_unit = oi_label.split('(')[-1].split('/')[1].replace(')','').strip()
            perf_den_unit = perf_label.split('(')[-1].split('/')[1].replace(')','').strip()
            beta_units_str = f"{oi_den_unit}/{perf_den_unit}"
        except Exception: beta_units_str = "BW Units"
        ax.text(text_pos_beta_data[0], text_pos_beta_data[1], 
                 f"Mem BW (β) = {beta_memory_bw:.2f} {beta_units_str}",
                 color='dimgray', horizontalalignment='center', verticalalignment='center', 
                 fontsize=10, rotation=angle_deg, rotation_mode='anchor',
                 transform_rotates_text=True)

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(oi_label, fontsize=14); ax.set_ylabel(perf_label, fontsize=14)
    ax.set_title(plot_title, fontsize=16, pad=20)
    
    all_plot_y_points = [pi_max_cpu_perf]
    if algorithms_data:
        for algo_name, points in algorithms_data.items():
            if points: all_plot_y_points.extend([p['performance'] for p in points])
    all_plot_y_points.extend(list(roofline_perf_ceiling)) 
    valid_y_points = [y for y in all_plot_y_points if y > 0 and np.isfinite(y)]
    if not valid_y_points: valid_y_points = [0.001, pi_max_cpu_perf if pi_max_cpu_perf > 0 else 1.0]
    final_min_perf = np.min(valid_y_points) / 2.0
    final_max_perf = np.max(valid_y_points) * 2.0
    if final_min_perf <= 0: final_min_perf = 0.0001 
    if final_max_perf <= final_min_perf: final_max_perf = final_min_perf * 1000
    
    ax.set_ylim(bottom=final_min_perf, top=final_max_perf)
    ax.set_xlim(left=min_oi, right=max_oi) 

    ax.legend(loc='best', fontsize=9, framealpha=0.7) 
    ax.grid(True, which="both", ls="-", alpha=0.3)
    plt.tight_layout()
    # plt.show()
    dt = datetime.now()
    ts = datetime.timestamp(dt) * 100
    ts = int(ts)
    output_dir = os.path.join("plots", "plot_images")
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, f"roofline_plot_{plot_title.split(' ')[-1]}_{ts}.png")
    try:
        plt.savefig(plot_filename, dpi=300)
        print(f"Plot saved as {plot_filename}")
    except Exception as e:
        print(f"Error saving plot: {e}")


# --- load_data_from_csv (as provided) ---
def load_data_from_csv(filepath):
    beta = None; pi_cpu = None; ois = []; perfs = []; labels = []
    parsing_data_points = False; header_map = {}
    oi_keys = {"oi", "operational_intensity", "operational intensity (flops/byte)", "intensity (flop/b)"}
    perf_keys = {"performance", "flops/cycle", "perf", "perf (flop/c)", "performance (flops/cycle)"}
    label_keys = {"label", "point_label", "name", "kernel", "kernel_name"}
    _temp_indices = {"oi_idx": -1, "perf_idx": -1, "label_idx": -1}
    try:
        with open(filepath, 'r', newline='', encoding='utf-8') as f:
            reader = csv.reader(f)
            for rnum, row in enumerate(reader):
                if not row or not row[0].strip() or row[0].strip().startswith('#'): continue
                clean_row = [c.strip() for c in row]; pk_csv = clean_row[0].lower().replace('_',' ')
                if not parsing_data_points:
                    if pk_csv=="beta memory bw" or pk_csv=="beta": beta=float(clean_row[1])
                    elif pk_csv=="pi max cpu perf" or pk_csv=="pi": pi_cpu=float(clean_row[1])
                    else:
                        th = [h.lower().strip() for h in clean_row]
                        if any(k in th for k in oi_keys) and any(k in th for k in perf_keys):
                            for i,hn in enumerate(th): header_map[hn]=i
                            for ks,tvn in [(oi_keys,"oi_idx"),(perf_keys,"perf_idx"),(label_keys,"label_idx")]:
                                found=False
                                for k_csv in ks:
                                    if k_csv in header_map: _temp_indices[tvn]=header_map[k_csv];found=True;break
                                if not found and tvn!="label_idx": raise ValueError(f"Header missing {tvn}")
                            parsing_data_points=True
                else:
                    if not header_map: raise ValueError("Data before header")
                    oi_idx,perf_idx,lbl_idx = _temp_indices["oi_idx"],_temp_indices["perf_idx"],_temp_indices["label_idx"]
                    try:
                        ois.append(float(clean_row[oi_idx])); perfs.append(float(clean_row[perf_idx]))
                        lbl = None
                        if lbl_idx!=-1 and len(clean_row)>lbl_idx and clean_row[lbl_idx]: lbl=clean_row[lbl_idx]
                        labels.append(lbl)
                    except (IndexError,ValueError) as e: print(f"Warn: Skip row {rnum+1} ({e}): {row}")
        if beta is None or pi_cpu is None: raise ValueError("beta/pi params missing")
        if not parsing_data_points and not ois: print("Warn: No data points/header")
        if all(l is None for l in labels): labels=None
    except FileNotFoundError: print(f"Err: No CSV {filepath}", file=sys.stderr); sys.exit(1)
    except ValueError as e: print(f"Err: Parse CSV {filepath}: {e}", file=sys.stderr); sys.exit(1)
    except Exception as e: print(f"Err: Read CSV {filepath}: {e}", file=sys.stderr); sys.exit(1)
    return beta, pi_cpu, ois, perfs, labels


# --- generate_plots_by_prefix function (MODIFIED to extract total_input_size for annotation) ---
def generate_plots_by_prefix(input_json_path, beta_bw, pi_perf, cli_oi_label, cli_perf_label):
    """
    Reads JSON, aggregates data by algorithm prefix, and generates one plot per prefix.
    Point annotations will now be the "total_input_size".
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

    data_grouped_by_prefix_then_full_key = defaultdict(lambda: defaultdict(list))

    for item_index, item in enumerate(all_test_data):
        if not isinstance(item, dict) or "test_case" not in item or "results" not in item:
            continue
        
        test_case = item.get("test_case", {})
        results_dict = item.get("results", {})

        try:
            # M, K, N, s are available from test_case if needed for other things,
            # but not directly used for the annotation label if "total_input_size" is primary.
            # m_val, k_val, n_val, s_val = test_case.get("M"), test_case.get("K"), test_case.get("N"), test_case.get("s")
            # if None in [m_val, k_val, n_val, s_val]:
            #     print(f"Warning: Skipping item (index {item_index}) for M,K,N,s: one or more are missing. Test case: {test_case}", file=sys.stderr)
            #     continue
            # m, k, n, s = int(m_val), int(k_val), int(n_val), int(s_val)

            for original_key, result_data in results_dict.items():
                if not isinstance(result_data, dict):
                    print(f"Warning: Result data for '{original_key}' in test_case index {item_index} is not a dict. Skipping.", file=sys.stderr)
                    continue

                oi_value = result_data.get("operational_intensity")
                perf_value = result_data.get("performance")
                total_input_size = result_data.get("total_input_size") # Get the new field

                if oi_value is None or perf_value is None or total_input_size is None:
                    print(f"Warning: Result '{original_key}' in test_case index {item_index} missing 'operational_intensity', 'performance', or 'total_input_size'. Data: {result_data}", file=sys.stderr)
                    continue
                
                try:
                    oi_float = float(oi_value)
                    perf_float = float(perf_value)
                    # Annotation label is now total_input_size (as string)
                    annotation_label = str(total_input_size) 
                except (ValueError, TypeError) as e:
                    print(f"Warning: Non-numeric OI/Perf/total_input_size for '{original_key}' in test_case index {item_index}. OI:'{oi_value}', Perf:'{perf_value}', Size:'{total_input_size}'. Error: {e}", file=sys.stderr)
                    continue

                prefix = original_key.split('_')[0]
                data_grouped_by_prefix_then_full_key[prefix][original_key].append({
                    "oi": oi_float,
                    "performance": perf_float,
                    "annotation_label": annotation_label # This is now total_input_size
                })
        except (TypeError, ValueError) as e:
            print(f"Warning: Error processing M,K,N,s for item index {item_index}: {e}", file=sys.stderr)
        except Exception as e_item:
             print(f"An unexpected error processing item (index {item_index}): {e_item}", file=sys.stderr)


    if not data_grouped_by_prefix_then_full_key:
        print("No valid data was aggregated from the JSON file. No plots will be generated.", file=sys.stdout)
        return

    for prefix, all_algos_data_for_this_prefix in data_grouped_by_prefix_then_full_key.items():
        if not all_algos_data_for_this_prefix:
            continue
            
        print(f"\nGenerating plot for prefix: {prefix}")
        plot_roofline_for_prefix(
            all_algos_data_for_this_prefix, 
            beta_bw, 
            pi_perf, 
            plot_title=f'Roofline Plot for Prefix: {prefix}',
            oi_label=cli_oi_label,
            perf_label=cli_perf_label)


# --- generate_plot_all (new function) ---
def generate_plot_all(input_json_path, beta_bw, pi_perf, cli_oi_label, cli_perf_label):
    """
    Reads JSON and generates a single roofline plot with all algorithm variants.
    """
    try:
        with open(input_json_path, 'r', encoding='utf-8') as f_json:
            all_test_data = json.load(f_json)
    except FileNotFoundError:
        print(f"Error: Input JSON file not found at '{input_json_path}'", file=sys.stderr); sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{input_json_path}'. Please check its format.", file=sys.stderr); sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading JSON: {e}", file=sys.stderr); sys.exit(1)

    if not isinstance(all_test_data, list):
        print("Error: Expected JSON data to be a list of test cases.", file=sys.stderr); sys.exit(1)

    from collections import defaultdict
    data_grouped = defaultdict(list)
    for item in all_test_data:
        if not isinstance(item, dict) or "results" not in item:
            continue
        results_dict = item.get("results", {})
        for original_key, result_data in results_dict.items():
            if not isinstance(result_data, dict):
                continue
            oi_value = result_data.get("operational_intensity")
            perf_value = result_data.get("performance")
            total_input_size = result_data.get("total_input_size")
            if oi_value is None or perf_value is None or total_input_size is None:
                continue
            try:
                oi_float = float(oi_value)
                perf_float = float(perf_value)
                annotation_label = str(total_input_size)
            except (ValueError, TypeError):
                continue
            data_grouped[original_key].append({
                "oi": oi_float,
                "performance": perf_float,
                "annotation_label": annotation_label
            })

    if not data_grouped:
        print("No valid data to plot.", file=sys.stdout)
        return

    print("Generating combined plot for all algorithms")
    plot_roofline_for_prefix(
        data_grouped,
        beta_bw,
        pi_perf,
        plot_title="Roofline Plot for All Algorithms",
        oi_label=cli_oi_label,
        perf_label=cli_perf_label
    )

# --- generate_roofline_csvs_from_json (for --per_test mode, kept as is but consider annotation change if desired) ---
def generate_roofline_csvs_from_json(input_json_path, output_directory, beta_bw, pi_perf, title_placeholder, cli_oi_label, cli_perf_label):
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

    for item_index, item in enumerate(all_test_data): 
        if not isinstance(item, dict) or "test_case" not in item or "results" not in item:
            print(f"Warning: Skipping invalid item (index {item_index}) in JSON: {item}", file=sys.stderr)
            continue
        
        test_case = item.get("test_case", {})
        results = item.get("results", {})

        try:
            m, k, n, s_val = (test_case.get(c) for c in ("M", "K", "N", "s"))
            if None in [m, k, n, s_val]:
                print(f"Warning: Skipping item (index {item_index}) due to missing M,K,N,s in test_case: {test_case}", file=sys.stderr)
                continue
            m, k, n, s_val = int(m), int(k), int(n), int(s_val)
            
            plot_specific_title = f"Roofline M{m}K{k}N{n}s{s_val}" 
            temp_csv_filename = f'temp_roofline_M{m}_K{k}_N{n}_s{s_val}_{os.getpid()}.csv' 
            temp_csv_filepath = os.path.join(output_directory, temp_csv_filename)
            
            with open(temp_csv_filepath, 'w', newline='', encoding='utf-8') as f_csv:
                writer = csv.writer(f_csv)
                writer.writerow(["parameter", "value"])
                writer.writerow(["beta_memory_bw", beta_bw])
                writer.writerow(["pi_max_cpu_perf", pi_perf])
                writer.writerow([]) 
                writer.writerow([f"# Data for M={m},K={k},N={n},s{s_val}"])
                writer.writerow(["oi", "performance", "label"]) 
                
                if not results:
                    print(f"Info: No 'results' for M{m}K{k}N{n}s{s_val}. Plot will show ceilings only.", file=sys.stdout)
                
                for result_key, res_dict in results.items():
                    if not isinstance(res_dict, dict): continue 
                    oi_val = res_dict.get("operational_intensity")
                    perf_val = res_dict.get("performance") 
                    # For --per_test mode, label is typically the algorithm name for the legend
                    # The annotation would be this label. If "total_input_size" should also be used
                    # here, the _original_plot_roofline would need to handle it.
                    # For now, keeping it as result_key for the 'label' column in CSV.
                    annotation_for_this_point = str(res_dict.get("total_input_size", result_key)) # Use total_input_size if available

                    if oi_val is None or perf_val is None:
                        continue
                    try:
                        # The CSV 'label' column will be used by _original_plot_roofline's `point_labels`
                        # This determines legend entries and default text annotation.
                        # If you want total_input_size as annotation, _original_plot_roofline needs
                        # to be modified to expect/use it. For now, it expects algo keys as labels.
                        writer.writerow([float(oi_val), float(perf_val), result_key]) 
                    except (ValueError, TypeError):
                        continue 
            
            _, _, ois_loaded, perfs_loaded, labels_loaded = load_data_from_csv(temp_csv_filepath)
            
            try: os.remove(temp_csv_filepath)
            except OSError as e_rm: print(f"Warning: Could not remove temp file '{temp_csv_filepath}': {e_rm}", file=sys.stderr)

            _original_plot_roofline(ois_loaded, perfs_loaded, 
                        beta_bw, pi_perf, 
                        point_labels=labels_loaded, 
                        plot_title=plot_specific_title,
                        oi_label=cli_oi_label,
                        perf_label=cli_perf_label)
        except (TypeError, ValueError) as e:
            print(f"Warning: Error in test_case M,K,N,s for item index {item_index}: {e}", file=sys.stderr)
        except Exception as e_item:
            print(f"An unexpected error processing item (index {item_index}): {e_item}", file=sys.stderr)

# --- _original_plot_roofline (Placeholder - PASTE YOUR PREVIOUS VERSION HERE) ---
def _original_plot_roofline(measured_op_intensities, measured_performances,
                  beta_memory_bw, pi_max_cpu_perf,
                  point_labels=None, plot_title="Roofline Plot",
                  oi_label="Operational Intensity (Flops/Byte)",
                  perf_label="Performance (Flops/Cycle)"):
    # This is the plot_roofline function version that was used with the --per_test mode.
    # It plots individual points, and 'point_labels' are algorithm names for legend/annotation.
    # If you want 'total_input_size' as annotation here too, this function needs modification.
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor('#f0f0f0')
    ridge_point_oi = pi_max_cpu_perf / beta_memory_bw if beta_memory_bw > 0 else float('inf')
    # Simplified OI range calculation for this placeholder
    min_oi = 0.01; max_oi = 100.0
    if measured_op_intensities is not None and len(measured_op_intensities) > 0:
        min_oi_data = np.min(measured_op_intensities)
        max_oi_data = np.max(measured_op_intensities)
        min_oi = min(min_oi_data / 5, ridge_point_oi / 10 if ridge_point_oi != float('inf') else min_oi_data / 5, 0.01)
        max_oi = max(max_oi_data * 5, ridge_point_oi * 10 if ridge_point_oi != float('inf') else max_oi_data * 5, min_oi * 1000)
    if min_oi <= 0: min_oi = 0.001 
    if max_oi <= min_oi : max_oi = min_oi * 1000 

    oi_range = np.logspace(np.log10(min_oi), np.log10(max_oi), num=200)
    roofline_perf = np.minimum(pi_max_cpu_perf, beta_memory_bw * oi_range)
    unbounded_beta_perf = beta_memory_bw * oi_range
    ax.plot(oi_range, unbounded_beta_perf, linestyle=':', color='dimgray', lw=1.5, label=f"Unbounded Mem BW (β={beta_memory_bw:.2f})")
    ax.axhline(y=pi_max_cpu_perf, color='gray', linestyle='--', lw=1.5, label=f"Peak Perf (π={pi_max_cpu_perf:.2f})")
    ax.plot(oi_range, roofline_perf, 'k-', lw=2.5, label='Attainable Performance Ceiling')
    
    if measured_op_intensities is not None and len(measured_op_intensities) > 0:
        plot_points_present = True
        # In --per_test mode, point_labels are the algorithm names.
        # Each is a distinct series for the legend. No connecting lines between different algos.
        if point_labels and len(point_labels) == len(measured_op_intensities):
            unique_algo_keys = sorted(list(set(l for l in point_labels if l))) # Get unique algo names for consistent coloring
            cmap = plt.get_cmap('tab10', max(1, len(unique_algo_keys)))
            color_map = {key: cmap(i) for i, key in enumerate(unique_algo_keys)}

            for i in range(len(measured_op_intensities)):
                algo_key_label = point_labels[i]
                color = color_map.get(algo_key_label, 'k') # Default color if somehow not in map
                ax.plot(measured_op_intensities[i], measured_performances[i], 
                        marker='o', linestyle='None', markersize=9, color=color, label=algo_key_label if algo_key_label else "Point")
                # Annotation for --per_test is typically the algorithm key itself if not further specified
                # If you want total_input_size here, it needs to be passed in point_labels or another arg.
                if algo_key_label:
                    ax.text(measured_op_intensities[i] * 1.1, measured_performances[i], algo_key_label, fontsize=7, va='center')
        else:
            ax.plot(measured_op_intensities, measured_performances, 'o', markersize=9, label='Measured Points')

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(oi_label); ax.set_ylabel(perf_label); ax.set_title(plot_title)
    # Remove duplicate legend entries if any caused by plotting individual points with same label
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles)) # Create a dict to ensure unique labels
    ax.legend(by_label.values(), by_label.keys(), loc='best', fontsize=8)
    ax.grid(True, which="both", ls="-", alpha=0.3); plt.tight_layout(); plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Roofline plots from JSON data.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""\
For prefix aggregation (default): Generates one plot per algorithm prefix.
  Legend shows full algorithm names. Point annotations show total_input_size.

For per-test-case plots: Use --per_test. Generates one plot for each
  test case (M,K,N,s combination) showing all its algorithms.
  Point annotations are currently the algorithm names (can be adapted if total_input_size is written to temp CSV for this mode).

JSON 'results' for all modes should contain 'operational_intensity', 
'performance' (fpc), and 'total_input_size' for each algorithm variant.
"""
    )
    parser.add_argument("json_filepath", help="Path to the JSON file containing roofline data.")
    parser.add_argument("--output_dir", default='./', help="Output directory (used for temporary files in --per_test mode). Default: ./")
    parser.add_argument("--beta", type=float, default=24.0, help="Beta (memory bandwidth) value (e.g., Bytes/Cycle). Default: 24.0")
    parser.add_argument("--pi", type=float, default=4.0, help="Pi (max CPU performance) value (e.g., Flops/Cycle). Default: 4.0")
    parser.add_argument("--per_test", action='store_true', help="Generate one plot per test case. If not set, aggregates by algorithm prefix.")
    parser.add_argument("--oi_label", help="Label for the OI axis.", default="Operational Intensity (Flops/Byte)")
    parser.add_argument("--perf_label", help="Label for the Performance axis.", default="Performance (Flops/Cycle)")
    parser.add_argument("-a", "--all", action="store_true", help="Generate a single plot containing all algorithms in the JSON.")
    
    args = parser.parse_args()
    
    if args.all:
        print("Operating in all mode: Generating one plot with all algorithms.")
        generate_plot_all(args.json_filepath, args.beta, args.pi, args.oi_label, args.perf_label)
    elif args.per_test:
        print("Operating in --per_test mode: Generating one plot per M,K,N,s combination.")
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir, exist_ok=True)
        generate_roofline_csvs_from_json(
            args.json_filepath,
            args.output_dir,
            args.beta,
            args.pi,
            "",
            args.oi_label,
            args.perf_label
        )
    else:
        print("Operating in default mode: Aggregating by algorithm prefix, one plot per prefix.")
        generate_plots_by_prefix(
            args.json_filepath,
            args.beta,
            args.pi,
            args.oi_label,
            args.perf_label
        )