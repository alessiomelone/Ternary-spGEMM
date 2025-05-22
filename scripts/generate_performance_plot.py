import argparse, csv, itertools, os, re, subprocess, sys, time
from pathlib import Path
from typing import List, Tuple, Dict
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

# ---------- constants --------------------------------------------------------
M_LIST = [1, 16, 64, 256, 1000, 4000, 16000, 64000]
K_LIST = [512, 1024, 2048, 4096, 2048, 4096, 8192, 16384]
N_LIST = [2048, 4096, 8192, 16384, 512, 1024, 2048, 4096]
SPARSITY_TO_SFLAG = {   # non-zeros = 1 / 2**sflag
    0.50: 1,   # 50 %  nnz  → -s 1
    0.25: 2,   # 75 % sparsity
    0.125: 3,  # 87 % sparsity
    0.0625: 4  # 93 % sparsity
}
DATE = (lambda b: next(f"{b}{'' if i==1 else f'({i})'}" for i in range(1,100) if not os.path.exists(f"results/{b if i==1 else f'{b}({i})'}")))(datetime.now().strftime("%A %-d.%-m, %-I:%M%p").lower())
DEFAULT_FUNCTION = "CSC_base"
BIN_PATH = Path("cpp_impl/SparseGEMM.out")
CSV_PATH = Path(f"results/{DATE}/perf_results.csv")
PLOT_PATH = Path(f"results/{DATE}/performance.png")

PERF_RE = re.compile(r"Running:\s*(?P<name>\S+).*?Performance:\s*(?P<flops>[0-9.]+)",
                     re.S)

# ---------- helpers ----------------------------------------------------------
def build_binary(force: bool = False):
    if BIN_PATH.exists() and not force:
        return
    print("[build] compiling SparseGEMM …")
    cmd = (
        "g++ -O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG "
        "cpp_impl/main.cpp cpp_impl/perf.cpp "
        "-o cpp_impl/SparseGEMM.out -DPMU"
    )
    subprocess.run(cmd.split(), check=True)
    print("       done.")

def run_case(M: int, K: int, N: int, s_flag: int, kernel: str) -> float:
    """Return flops/cycle for given run or np.nan if failed/not found."""
    cmd = ["sudo", str(BIN_PATH), "-M", str(M), "-K", str(K), "-N", str(N), "-s", str(s_flag)]
    try:
        out = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(f"[warn] run failed: {e}")
        return np.nan

    match = PERF_RE.search(out)
    while match:
        if match.group("name") == kernel:
            return float(match.group("flops"))
        # keep searching (multiple kernels in same output)
        out = out[match.end():]
        match = PERF_RE.search(out)
    return np.nan

def n_elements(M: int, K: int, N: int) -> int:
    """Total elements in X and W matrices (M*K + K*N)."""
    return M * K + K * N

def parse_range(rng: str, max_len: int):
    """Parse “a-b” into (a-1 .. b-1) inclusive zero-based range."""
    lo, hi = map(int, rng.split('-', 1))
    if not (1 <= lo <= hi <= max_len):
        raise ValueError(f"range must be between 1 and {max_len}")
    return range(lo-1, hi)

def main():
    parser = argparse.ArgumentParser(description="SparseGEMM benchmark sweep")
    parser.add_argument("m_range", nargs="?", help="1-based index range for M/K/N groups, e.g. '1-5'")
    parser.add_argument("s_range", nargs="?", help="sparsity-exponent range, e.g. '2-3' (1/2**2..1/2**3 nnz)")
    parser.add_argument("--kernel", default=DEFAULT_FUNCTION, help="name as printed after 'Running:'")
    parser.add_argument("--rebuild", action="store_true", help="force binary rebuild")
    parser.add_argument("--skip_errors", action="store_true", help="ignore shapes that crash")
    args = parser.parse_args()

    # build if needed
    build_binary(args.rebuild)

    # determine which M/K/N groups to sweep
    m_idxs = parse_range(args.m_range, len(M_LIST)) if args.m_range else range(len(M_LIST))

    # determine which sparsity flags to sweep
    all_sflags = sorted(SPARSITY_TO_SFLAG.values())
    if args.s_range:
        s_idxs = parse_range(args.s_range, len(all_sflags))
        sflags = [all_sflags[i] for i in s_idxs]
    else:
        sflags = all_sflags

    # precompute full-span n_elements for fixed x-axis limits
    full_n = [n_elements(M_LIST[i], K_LIST[i], N_LIST[i]) for i in range(len(M_LIST))]
    min_n, max_n = min(full_n), max(full_n)

    results: List[Tuple[int, float, float]] = []

    for i in m_idxs:
        M = M_LIST[i]; K = K_LIST[i]; N = N_LIST[i]
        for sflag in sflags:
            sparsity = next(p for p,q in SPARSITY_TO_SFLAG.items() if q == sflag)
            print(f"[run] i={i} M={M:<6} K={K:<6} N={N:<6} sparsity={sparsity:>6.3f}", end="  ")
            flops = run_case(M, K, N, sflag, args.kernel)
            if np.isnan(flops):
                print("— skipped")
                if not args.skip_errors:
                    print("   (use --skip_errors to continue)")
                    sys.exit(1)
            else:
                print(f"-> {flops:.3f} flops/cycle")
                results.append((n_elements(M, K, N), flops, sparsity))

    # dump CSV
    if results:
        CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
        with CSV_PATH.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["n_elements", "flops_per_cycle", "sparsity"])
            writer.writerows(results)
        print(f"[save] wrote {CSV_PATH}")

    # make plot
    if results:
        results.sort(key=lambda x: x[0])
        xs = np.array([r[0] for r in results])
        ys = np.array([r[1] for r in results])
        ss = np.array([r[2] for r in results])

        plt.figure(figsize=(8,6), dpi=150)
        for sparsity in sorted(SPARSITY_TO_SFLAG):
            mask = ss == sparsity
            plt.plot(xs[mask], ys[mask], marker="o",
                     label=f"{int((1-sparsity)*100)}% zeros")

        plt.xscale("log", base=2)
        plt.xlim(min_n, max_n)
        plt.xlabel("Total elements in X and W (M*K + K*N)")
        plt.ylabel("flops / cycle")
        plt.title(f"SparseGEMM – {args.kernel}")
        plt.grid(True, which="both", ls=":")
        plt.legend()
        plt.tight_layout()
        PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(PLOT_PATH)
        print(f"[save] wrote {PLOT_PATH}")


if __name__ == "__main__":
    main()
