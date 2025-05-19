#!/bin/bash

# Default values
DEFAULT_PROFILER_ABS_PATH=/home/kali-user-1/Downloads/DynamoRIO-AArch64-Linux-11.3.0-1/bin64/drrun
DEFAULT_OUT_FILE_PATH=../bin/a.out
DEFAULT_CUSTOM_CFLAGS=""
DEFAULT_PROFILER_OUTPUT_FILE="../measurements/prof_out.txt" # Default to current directory for profiler output

# Initialize variables with default values
PROFILER_ABS_PATH="$DEFAULT_PROFILER_ABS_PATH"
OUT_FILE_PATH="$DEFAULT_OUT_FILE_PATH"
CUSTOM_CFLAGS="$DEFAULT_CUSTOM_CFLAGS"
PROFILER_OUTPUT_FILE="$DEFAULT_PROFILER_OUTPUT_FILE"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --ccflags)
        CUSTOM_CFLAGS="$2"
        shift # past argument
        shift # past value
        ;;
        --out)
        OUT_FILE_PATH="$2"
        shift # past argument
        shift # past value
        ;;
        --prof)
        PROFILER_ABS_PATH="$2"
        shift # past argument
        shift # past value
        ;;
        --prof_out_file)
        PROFILER_OUTPUT_FILE="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        echo "Unknown option: $1"
        echo "Usage: $0 [--ccflags \"<flags>\"] [--out <path>] [--prof <path>] [--prof_out_file <dir_path>]"
        exit 1
        ;;
    esac
done

# Output the configuration being used
echo "Compiler Flags: ${CUSTOM_CFLAGS:-<none>}"
echo "Output File Path: $OUT_FILE_PATH"
echo "Profiler Absolute Path: $PROFILER_ABS_PATH"
echo "Profiler Output Directory: $PROFILER_OUTPUT_FILE"
echo ""

# Compilation command
COMPILE_COMMAND="g++ -O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG $CUSTOM_CFLAGS -I../.. ../../comp.cpp ../src/main_test_cache.cpp -o $OUT_FILE_PATH"
echo "Executing compilation command:"
echo "$COMPILE_COMMAND"
eval "$COMPILE_COMMAND"

# Check if compilation was successful
if [ $? -ne 0 ]; then
  echo "Compilation failed."
  exit 1
fi
echo "Compilation successful."
echo ""

# Profiler command
# Note: drmemtrace typically creates a subdirectory within the specified output directory.
# The actual trace files will be inside a folder like $PROFILER_OUTPUT_FILE/drmemtrace.<app_name>.<pid>.<tid>.dir/
PROFILER_COMMAND="$PROFILER_ABS_PATH -t drmemtrace -- \"$OUT_FILE_PATH\""
echo "Executing profiler command:"
echo "$PROFILER_COMMAND"
eval "$PROFILER_COMMAND" > "$PROFILER_OUTPUT_FILE" 2>&1

if [ $? -ne 0 ]; then
  echo "Profiler execution failed."
  exit 1
fi
echo "Profiler execution finished. Output should be in '$PROFILER_OUTPUT_FILE'."
exit 0