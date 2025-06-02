# Default parameters
VARY='M'
SPARSITY=8

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --vary) VARY="$2"; shift ;;
    --sparsity) SPARSITY="$2"; shift ;;
    *) echo "Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

DATESTRING=$(date +"%A, %-d.%-m %-I:%M%P")
OUTPUT_JSON="plots/roofline/output_jsons/$DATESTRING.json"

make clean
make INSTRUMENT=1
sudo python plots/run_benchmark.py --vary "$VARY" --sparsity "$SPARSITY" -s --output "$OUTPUT_JSON"

# OUTPUT_JSON="Saturday, 24.5 5:56P.json"
python plots/roofline/plot_perf.py "$OUTPUT_JSON"
python plots/roofline/plot_roofline.py "$OUTPUT_JSON" -a