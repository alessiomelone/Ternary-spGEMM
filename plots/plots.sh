DATESTRING=$(date +"%A, %-d.%-m %-I:%M%P")
OUTPUT_JSON="plots/output_jsons/$DATESTRING.json"

if [ "$1" = "--varyonly" ]; then
    shift
    VARYONLY=$1
else
    VARYONLY=$1
fi

if [[ "$VARYONLY" != "M" && "$VARYONLY" != "K" && "$VARYONLY" != "N" ]]; then
    VARYONLY=""
fi

if [ -n "$VARYONLY" ]; then
    VAR_ARG="--varyonly $VARYONLY"
else
    VAR_ARG=""
fi

make clean
make INSTRUMENT=1
sudo python plots/run_benchmark.py -s --output "$OUTPUT_JSON" $VAR_ARG

if [ "$VARYONLY" = "M" ]; then
    X_LABEL="M"
    TITLE="Function Performance vs. M"
elif [ "$VARYONLY" = "K" ]; then
    X_LABEL="K"
    TITLE="Function Performance vs. K"
elif [ "$VARYONLY" = "N" ]; then
    X_LABEL="N"
    TITLE="Function Performance vs. N"
else
    X_LABEL="Total Input Size (MB)"
    TITLE="Function Performance vs. Total Input Size"
fi
python plots/plot_perf.py "$OUTPUT_JSON" --xlabel "$X_LABEL" --title "$TITLE"
python plots/plot_roofline.py "$OUTPUT_JSON" -a