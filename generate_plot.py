import pandas as pd
import os
import glob
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
#  Configuration
# ---------------------------------------------------------------------------
results_dir = "./results/"
list_of_files = glob.glob(os.path.join(results_dir, '*.csv')) # Get all csv files in the directory
if not list_of_files:
    raise FileNotFoundError(f"No CSV files found in {results_dir}")
csv_path = max(list_of_files, key=os.path.getmtime) # Find the latest file based on modification time
print(f"Using data from: {csv_path}")
K_FIXED, N_FIXED    = 512, 2048   # fixed dimensions
NONZERO_FIXED       = 4           # sparsity level to keep
# ---------------------------------------------------------------------------

# 1.  Load the data
df = pd.read_csv(csv_path)

# 2.  Filter the slice we care about
df_slice = df[
    (df["K"] == K_FIXED) &
    (df["N"] == N_FIXED) &
    (df["nonZero"] == NONZERO_FIXED)
].copy()

# 3.  Compute FLOPs and FLOPs/cycle
df_slice["flops"] = (
    df_slice["M"] *
    df_slice["N"] *
    (1 + df_slice["K"] / df_slice["nonZero"])
)
df_slice["flops_per_cycle"] = df_slice["flops"] / df_slice["pmu_cycles"]

# 4.  Find an optimisation‑variant column (if any)
variant_column = None
for col in ["compiler", "variant", "flags", "opt", "optimization", "config"]:
    if col in df_slice.columns:
        variant_column = col
        break
if variant_column is None:            # no such column? invent one
    variant_column = "variant"
    df_slice[variant_column] = "run"

# 5.  Average repeats → one point per (M, variant)
df_plot = (
    df_slice
    .groupby(["M", variant_column], as_index=False)["flops_per_cycle"]
    .mean()
    .sort_values("M")
)

# 6.  Plot
plt.figure(figsize=(8, 4))
markers = ["o", "s", "^", "d", "v", ">", "<", "p", "h", "x", "+"]

for idx, (variant, block) in enumerate(df_plot.groupby(variant_column)):
    block = block.sort_values("M")
    plt.plot(
        block["M"],
        block["flops_per_cycle"],
        marker=markers[idx % len(markers)],
        linewidth=2,
        label=variant,
    )
plt.xscale('log', base=2)
plt.yscale('log', base=2)
plt.title("Performance (flops/cycle)")
plt.xlabel("M")
plt.ylabel("flops/cycle")
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.show()
