import sys
import re
import argparse

def subtract_files_line_by_line(results_filepath, base_filepath, visible_fields=None, excl=None):
    """
    Reads two cache simulation files line by line, subtracts numerical data
    from the base file from the results file, and prints the Simultaneous results.
    Percentage fields are taken directly from the results file.
    Only specified 'visible_fields' are printed if the option is used.
    Assumes both files have the same number of lines and structure for subtraction.
    """
    try:
        with open(results_filepath, 'r') as results_file, \
             open(base_filepath, 'r') as base_file:

            for line_res in results_file:
                line_base = base_file.readline()

                if not line_base:
                    # Base file is shorter, print rest of results file as is
                    print(line_res, end="")
                    continue

                # Attempt to parse the line as a field with a label and value
                # Regex breakdown:
                # ^(?P<indent>\s*)                  # Capture leading whitespace (indentation)
                # (?P<label_full>[A-Za-z0-9\s_.-]+?:\s*) # Capture label text (non-greedy) up to colon, and trailing spaces
                # (?P<value_spaces>\s*)            # Capture spaces between label and value (if any, usually caught by label_full's \s*)
                # (?P<value>[\d.]+)                # Capture numerical value (digits and decimal point)
                # (?P<percent>%?)                  # Capture optional percentage sign
                # (?P<suffix>.*)$                   # Capture any remaining characters on the line (suffix)
                match_res = re.match(r"^(?P<indent>\s*)(?P<label_full>[A-Za-z0-9\s_.,()-]+?:\s*)(?P<value_spaces>\s*)(?P<value>[\d.]+)(?P<percent>%?)(?P<suffix>.*)$", line_res.rstrip('\n'))

                if match_res:
                    res_parts = match_res.groupdict()
                    original_value_str_res = res_parts['value']
                    is_percentage = bool(res_parts['percent'])
                    
                    # Clean the label by removing the colon and stripping whitespace
                    # e.g., "  Cache Hits:   " -> "Cache Hits"
                    label_text = res_parts['label_full'].strip().removesuffix(':').strip()
                    should_process_and_print = False

                    if not is_percentage and visible_fields is not None: # visible_fields is not None means --fields was used
                        for f in visible_fields:
                            if f.lower() in label_text.lower():
                                should_process_and_print = True
                        if excl is not None:
                            for f in excl:
                                if f.lower() in label_text.lower():
                                    should_process_and_print = False

                    # print(label_text, should_process_and_print)
                    if should_process_and_print:
                        if is_percentage:
                            print(line_res, end="")
                        else:
                            # Now try to match the base line with a similar structure to extract its value
                            match_base = re.match(r"^(?P<indent>\s*)(?P<label_full>[A-Za-z0-9\s_.,()-]+?:\s*)(?P<value_spaces>\s*)(?P<value>[\d.]+)(?P<percent>%?)(?P<suffix>.*)$", line_base.rstrip('\n'))
                            
                            if match_base:
                                base_parts = match_base.groupdict()
                                original_value_str_base = base_parts['value']

                                try:
                                    if '.' in original_value_str_res or '.' in original_value_str_base:
                                        val_res = float(original_value_str_res)
                                        val_base = float(original_value_str_base)
                                        diff = val_res - val_base
                                        # Try to preserve original number of decimal places from results file
                                        if '.' in original_value_str_res:
                                            decimals = len(original_value_str_res.split('.')[-1])
                                            diff_str = f"{diff:.{decimals}f}"
                                        elif '.' in original_value_str_base: # use base if res is int but base is float
                                            decimals = len(original_value_str_base.split('.')[-1])
                                            diff_str = f"{diff:.{decimals}f}"
                                        else: # both looked like ints but were parsed as float
                                            diff_str = str(int(diff))
                                    else:
                                        val_res = int(original_value_str_res)
                                        val_base = int(original_value_str_base)
                                        diff = val_res - val_base
                                        diff_str = str(diff)

                                    prefix_for_print = f"{res_parts['indent']}{res_parts['label_full']}{res_parts['value_spaces']}"
                                    suffix_for_print = res_parts['suffix'] if res_parts['suffix'] else ''
                                    print(f"{prefix_for_print}{diff_str}{res_parts['percent']}{suffix_for_print}")

                                except ValueError:
                                    # Problem converting numbers, print result line as is
                                    print(line_res, end="")
                            else:
                                # Base line didn't match, print result line as is
                                print(line_res, end="")
                    # else: field was filtered out, do nothing
                else:
                    # Line doesn't match the expected field format (e.g., a header, title, or empty line)
                    # Print it as is, regardless of filtering options for fields
                    print(line_res, end="")
            
            # Print any remaining lines from results_file if base_file was shorter
            for line_res_remaining in results_file:
                print(line_res_remaining, end="")


    except FileNotFoundError:
        print(f"Error: One or both files not found. Searched for '{results_filepath}' and '{base_filepath}'.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

# --- Main execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Subtracts numerical data from a base cache simulation file from a results file. "
                    "Allows specifying which leaf fields (e.g., 'Cache Hits') should be visible. "
                    "Percentage fields and lines not matching the field structure are always shown."
    )
    parser.add_argument("results_filepath", help="Path to the results cache simulation file.")
    parser.add_argument("base_filepath", help="Path to the base cache simulation file.")
    parser.add_argument(
        "--fields",
        nargs="*",  # 0 or more arguments
        default=None, # Default to None if --fields is not provided
        help="Optional list of exact leaf field names to display (e.g., \"Cache Hits\" \"L1 Misses\"). "
             "Field names are matched after removing the trailing colon and extra spaces. "
             "If this option is provided but no names are listed, no numeric fields will be shown. "
             "If this option is omitted, all numeric fields are shown."
    )
    parser.add_argument(
        "--excl",
        nargs="*",  # 0 or more arguments
        default=None, # Default to None if --fields is not provided
        help="Optional list of exact leaf field names to display (e.g., \"Cache Hits\" \"L1 Misses\"). "
             "Field names are matched after removing the trailing colon and extra spaces. "
             "If this option is provided but no names are listed, no numeric fields will be shown. "
             "If this option is omitted, all numeric fields are shown."
    )

    args = parser.parse_args()

    # The regex for label_full was updated to be more inclusive: [A-Za-z0-9\s_.,()-]+?:\s*
    # This helps capture labels like "L1-dcache hits:", "Total cycles (non-spec.)" etc.
    # The cleaning `label_text = res_parts['label_full'].strip().removesuffix(':').strip()`
    # is designed to normalize this for comparison.

    subtract_files_line_by_line(args.results_filepath, args.base_filepath, args.fields, args.excl)

    # To redirect the output of this script to a file, you would do it from your shell:
    # python your_script_name.py results.txt base.txt --fields "Field Name 1" "Another Field" > diff_output.txt