# Run:
```
$ sudo python3 run_benchmark.py -s --output benchmarks_whatevername.json
$ python3 plot_roofline.py benchmarks_whatevername.json
$ python3 plot_perf.py benchmarks_whatevername.json --title="Performance vs Input Size"
```
# TODO: 
* Maybe fix axis points
* Save roofline plots
<!-- Right now it's the hardcoded value for what should be the total memory requirements for each data structure. -->