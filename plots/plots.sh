DATESTRING=$(date +"%A, %-d.%-m %-I:%M%P")
OUTPUT_JSON="plots/roofline/output_jsons/$DATESTRING.json"

make clean
make INSTRUMENT=1
sudo python plots/roofline/run_benchmark.py -s --output "$OUTPUT_JSON"

# OUTPUT_JSON="Saturday, 24.5 5:56P.json"
python plots/roofline/plot_perf.py "$OUTPUT_JSON"
python plots/roofline/plot_roofline.py "$OUTPUT_JSON" -a