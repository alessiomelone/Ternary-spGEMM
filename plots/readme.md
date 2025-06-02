# Prerun:
Edit `run_benchmark.py` to modify the test cases (choices of M, K, N, s) you wanna measure.

# Run:
```
$ sudo python3 run_benchmark.py -s --output benchmarks_whatevername.json
$ python3 plot_roofline.py benchmarks_whatevername.json
$ python3 plot_perf.py benchmarks_whatevername.json --title="Performance vs Input Size"
```
# TODO: 
* Maybe fix axis points
* Save roofline plots