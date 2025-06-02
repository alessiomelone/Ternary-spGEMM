import re
import sys

def print_calculated_miss_ratio(header_indent, cache_name, cache_type, hits, misses, child_hits):
    """
    Calculates and prints the cache miss ratio for a given cache section.
    This function is called after all lines for that cache section have been printed.

    Args:
        header_indent (str): The indentation string of the cache section's header line.
        cache_name (str): The name of the cache (e.g., "L1I0", "LL").
        cache_type (str): The type of cache ("L1" or "LL").
        hits (int): Number of cache hits.
        misses (int): Number of cache misses.
        child_hits (int): Number of child hits (relevant for LL cache).
    """
    denominator = 0
    if cache_type == "LL":
        denominator = hits + child_hits + misses
    elif cache_type == "L1":
        denominator = hits + misses
    else: # Should not happen with expected format
        miss_ratio_str = "N/A (unknown cache type)"
        print(f"{header_indent}  Calculated Miss Ratio: {miss_ratio_str}")
        return

    if denominator == 0:
        miss_ratio_str = "N/A (division by zero: H={hits}, M={misses}, CH={child_hits})"
    else:
        miss_ratio = (misses / denominator) * 100
        miss_ratio_str = f"{miss_ratio:.2f}% (M={misses}, H={hits}"
        if cache_type == "LL":
            miss_ratio_str += f", CH={child_hits}"
        miss_ratio_str += ")"


    # The "Calculated Miss Ratio:" line should be indented similarly to
    # the "Hits:", "Misses:" lines within that section.
    # These stat lines are typically indented 2 spaces more than their header.
    summary_indent = header_indent + "  "
    print(f"{summary_indent}Calculated Miss Ratio: {miss_ratio_str}")

def process_log_file(filepath):
    """
    Processes the cache simulation log file to calculate and print miss ratios.

    Args:
        filepath (str): The path to the log file.
    """
    # Regex to identify the start of a cache stats section header.
    # Captures:
    #   group(1): Indentation of the header line.
    #   group(2): Cache name (e.g., "L1I0", "LL").
    cache_header_regex = re.compile(r"^(\s*)([A-Za-z0-9]+)\s*\(.*stats:$")
    
    # Regex to capture Hits, Misses, Child hits lines.
    # Captures:
    #   group(1): Label ("Hits", "Misses", or "Child hits").
    #   group(2): Numerical value.
    stats_line_regex = re.compile(r"^\s*(Hits|Misses|Child hits):\s*([\d]+)")

    # Variables to store stats for the current cache section being processed
    current_cache_header_indent = ""
    current_cache_name = None
    current_cache_type = None # "L1" or "LL"
    current_hits = 0
    current_misses = 0
    current_child_hits = 0 # Specific to LL cache
    last_lines = []
    skip_lines = ['Dummy', 'Attempting to fill', 'application exited']
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if '0 thread' in line:
                    last_lines.append(line)
                    continue
                skip = False
                for s in skip_lines:
                    if s in line:
                        skip = True
                        break
                if skip:
                    continue

                match_header = cache_header_regex.match(line)
                match_stats = stats_line_regex.match(line)

                if match_header:
                    # This line is a cache section header.
                    # If we were processing a previous cache section, print its miss ratio now.
                    if current_cache_name:
                        print_calculated_miss_ratio(current_cache_header_indent, 
                                                    current_cache_name, current_cache_type,
                                                    current_hits, current_misses, current_child_hits)

                    # Start tracking the new cache section
                    current_cache_header_indent = match_header.group(1)
                    current_cache_name = match_header.group(2)
                    
                    if current_cache_name.startswith("L1"):
                        current_cache_type = "L1"
                    elif current_cache_name == "LL": # Assuming "LL" is the exact name for Last Level Cache
                        current_cache_type = "LL"
                    else:
                        # This case should ideally not be reached if the log format is consistent
                        # and only L1x and LL caches are present.
                        current_cache_type = "Unknown" 
                    
                    # Reset stats for the new section
                    current_hits = 0
                    current_misses = 0
                    current_child_hits = 0
                
                elif current_cache_name and match_stats:
                    # This line is a stats line (Hits, Misses, Child hits)
                    # within the current cache section.
                    label = match_stats.group(1)
                    try:
                        value = int(match_stats.group(2))
                        if label == "Hits":
                            current_hits = value
                        elif label == "Misses":
                            current_misses = value
                        elif label == "Child hits":
                            # Only assign to child_hits if it's an LL cache,
                            # otherwise it's an unexpected field for L1.
                            if current_cache_type == "LL":
                                current_child_hits = value
                            else:
                                sys.stderr.write(f"Warning: 'Child hits' found for non-LL cache '{current_cache_name}' in line: {line.strip()}\n")
                    except ValueError:
                        sys.stderr.write(f"Warning: Could not parse numeric value in line: {line.strip()}\n")
                
                # If the line is not a cache header and not a stats line for the current cache,
                # it's either a line outside any cache section (like "Core #0...") or
                # a descriptive line within a cache section that isn't a stat we parse.
                # These lines are already printed by `print(line, end="")`.
                # The crucial part is that a new cache_header_regex match signals the end
                # of the previous section, or EOF signals the end of the last section.
                print(line, end="") # Print the current line from the log

            # After the loop, if there's an unfinished cache section (the last one in the file),
            # print its calculated miss ratio.
            if current_cache_name:
                print_calculated_miss_ratio(current_cache_header_indent, 
                                            current_cache_name, current_cache_type,
                                            current_hits, current_misses, current_child_hits)
            for l in last_lines:
                print(l, end='')

    except FileNotFoundError:
        sys.stderr.write(f"Error: Log file not found at '{filepath}'\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An unexpected error occurred: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <log_filepath>", file=sys.stderr)
        sys.exit(1) # Exit with an error code for incorrect usage
    
    log_filepath = sys.argv[1]
    process_log_file(log_filepath)
