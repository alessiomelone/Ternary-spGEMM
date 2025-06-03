
# Determine output JSON filename from first argument, default to timestamp if not provided
if [[ -n "$1" ]]; then
    FILE_NAME="$1"
    shift
    # Remove .json extension if provided
    if [[ "$FILE_NAME" == *.json ]]; then
        BASENAME="${FILE_NAME%.json}"
    else
        BASENAME="$FILE_NAME"
    fi
    OUTPUT_JSON="plots/output_jsons/${BASENAME}.json"
else
    DATESTRING=$(date +"%A, %-d.%-m %-I:%M%P")
    OUTPUT_JSON="plots/output_jsons/$DATESTRING.json"
fi

# Parse command-line flags
VARYONLY=""
SPARSITYONLY=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --varyonly)
            VARYONLY="$2"
            shift
            shift
            ;;
        --sparsityonly)
            SPARSITYONLY="$2"
            shift
            shift
            ;;
        *)
            shift
            ;;
    esac
done

# Validate VARYONLY
if [[ "$VARYONLY" != "M" && "$VARYONLY" != "K" && "$VARYONLY" != "N" ]]; then
    VARYONLY=""
fi

if [ -n "$VARYONLY" ]; then
    VAR_ARG="--varyonly $VARYONLY"
else
    VAR_ARG=""
fi

# Build sparsity argument
if [ -n "$SPARSITYONLY" ]; then
    SPARSE_ARG="--sparsityonly $SPARSITYONLY"
else
    SPARSE_ARG=""
fi

make clean
make INSTRUMENT=1
sudo python3 plots/run_benchmark.py -s --output "$OUTPUT_JSON" $VAR_ARG $SPARSE_ARG

if [ "$VARYONLY" = "M" ]; then
    X_LABEL="M"
    TITLE="Performance vs. M"
elif [ "$VARYONLY" = "K" ]; then
    X_LABEL="K"
    TITLE="Performance vs. K"
elif [ "$VARYONLY" = "N" ]; then
    X_LABEL="N"
    TITLE="Performance vs. N"
else
    X_LABEL="Total Input Size (MB)"
    TITLE="Performance vs. Total Input Size"
fi

TITLE="$TITLE on Apple M2, 3.49 GHz"

python plots/plot_perf.py "$OUTPUT_JSON" --xlabel "$X_LABEL" --title "$TITLE"
python plots/plot_roofline.py "$OUTPUT_JSON" -a