#!/usr/bin/env python3
"""
Safe multi-process simulation launcher

Usage examples:

# External executable that accepts --seed and --out arguments
python3 scripts/run_simulations.py \
  --executable ./bin/AvalancheRunner \
  --n 20 --parallel 6 \
  --exe-arg "--steps=10000" --exe-arg "--size=100" \
  --seed-start 1000 \
  --out-plot combined.png --out-csv combined.csv

# If executable writes CSV to stdout, use --capture-stdout to capture it
python3 scripts/run_simulations.py --executable ./sim_bin --capture-stdout --n 10

Requirements:
- Python 3.7+
- pip install pandas matplotlib

What it does:
- Launches up to --parallel simulations concurrently (each as its own OS process)
- Each simulation writes its data to a unique temporary CSV file
  (either by providing an output path via --out-arg to the executable or by capturing stdout)
- After all simulations finish, reads the CSV files into pandas, aligns them on a common time grid if needed,
  computes mean and standard deviation across runs, and writes a single combined CSV and a single plot
  (mean line + shaded std band).

Safety features:
- Uses a temporary run directory for outputs and atomically moves final outputs into place
- Handles KeyboardInterrupt by terminating child processes
- Retries failed runs (configurable retries)
- Validates that at least one CSV was produced before aggregating

Adapt the --out-arg name (default: --out) and --seed-arg name (default: --seed) to match the simulation binary's CLI.
"""

from __future__ import annotations
import argparse
import concurrent.futures
import csv
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import List, Optional, Tuple

# Lazy imports for plotting/aggregation
try:
    import pandas as pd
    import matplotlib.pyplot as plt
except Exception:
    pd = None  # type: ignore
    plt = None  # type: ignore


def ensure_deps():
    if pd is None or plt is None:
        print("Missing dependency: pandas and matplotlib are required. Install with: pip install pandas matplotlib", file=sys.stderr)
        sys.exit(2)


def atomic_move(src: Path, dest: Path) -> None:
    """Move file atomically when possible."""
    dest_parent = dest.parent
    dest_parent.mkdir(parents=True, exist_ok=True)
    temp_dest = dest_parent / (dest.name + ".tmp")
    shutil.move(str(src), str(temp_dest))
    os.replace(str(temp_dest), str(dest))


class ChildProcessManager:
    """Manage launched subprocesses so they can be terminated on KeyboardInterrupt."""

    def __init__(self):
        self._procs: List[subprocess.Popen] = []

    def add(self, p: subprocess.Popen):
        self._procs.append(p)

    def terminate_all(self):
        for p in list(self._procs):
            if p.poll() is None:
                try:
                    p.terminate()
                except Exception:
                    pass

    def kill_all(self):
        for p in list(self._procs):
            if p.poll() is None:
                try:
                    p.kill()
                except Exception:
                    pass


CHILD_MANAGER = ChildProcessManager()


def run_single_simulation(
    exe: str,
    idx: int,
    out_path: Optional[Path],
    seed: Optional[int],
    exe_args: List[str],
    out_arg_name: str,
    seed_arg_name: str,
    capture_stdout: bool,
    timeout: Optional[float],
    retries: int = 0,
) -> Tuple[int, int, Optional[Path], int, Optional[str]]:
    """
    Returns tuple (idx, attempts, out_path_or_none, returncode, stderr_text)
    """
    attempts = 0
    last_err: Optional[str] = None
    while attempts <= retries:
        attempts += 1
        cmd = [exe] + exe_args[:]
        if seed is not None:
            cmd += [f"{seed_arg_name}", str(seed)]
        if out_path is not None and not capture_stdout:
            cmd += [f"{out_arg_name}", str(out_path)]

        try:
            if capture_stdout:
                # Capture stdout and write to out_path ourselves
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                CHILD_MANAGER.add(p)
                stdout, stderr = p.communicate(timeout=timeout)
                rc = p.returncode
                stderr_text = stderr.decode(errors="replace") if stderr else None
                if rc == 0 and out_path is not None:
                    out_path.write_bytes(stdout)
                    return (idx, attempts, out_path, rc, stderr_text)
                elif rc == 0 and out_path is None:
                    # caller didn't ask for output file
                    return (idx, attempts, None, rc, stderr_text)
                else:
                    last_err = stderr_text
            else:
                # Let the child write the file directly
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                CHILD_MANAGER.add(p)
                _, stderr = p.communicate(timeout=timeout)
                rc = p.returncode
                stderr_text = stderr.decode(errors="replace") if stderr else None
                if rc == 0:
                    return (idx, attempts, out_path, rc, stderr_text)
                else:
                    last_err = stderr_text
        except subprocess.TimeoutExpired:
            try:
                p.kill()
            except Exception:
                pass
            last_err = f"Timeout after {timeout} seconds"
        except KeyboardInterrupt:
            # propagate so parent can handle
            raise
        except Exception as e:
            last_err = str(e)

        # small backoff before retry
        time.sleep(0.2 * attempts)

    return (idx, attempts, out_path, -1, last_err)


def aggregate_csvs(csv_paths: List[Path], time_col: Optional[str] = None) -> pd.DataFrame:
    """Read CSVs and aggregate. If a time_col is given, align by that.

    Returns DataFrame with columns: time (if any), mean, std, and optionally count.
    """
    dfs = []
    for p in csv_paths:
        try:
            df = pd.read_csv(p)
        except Exception:
            # try whitespace delimited
            df = pd.read_csv(p, delim_whitespace=True)
        dfs.append(df)

    if len(dfs) == 0:
        raise RuntimeError("No CSV files to aggregate")

    # If user provided a time column name or the first column looks like time
    if time_col is None:
        # heuristics: if first column named 'time' or 't' or 'step'
        for c in dfs[0].columns:
            if c.lower() in ("time", "t", "step", "steps"):
                time_col = c
                break
        if time_col is None:
            # fallback: use index positions (0..n-1) assuming same length
            lengths = {len(d) for d in dfs}
            if len(lengths) != 1:
                # align by a generated integer time index up to max length using NaN fill
                max_len = max(len(d) for d in dfs)
                aligned = pd.DataFrame()
                for i, d in enumerate(dfs):
                    col_name = f"run_{i}"
                    # if dataframe has a single data column, use it; else attempt to use the last column
                    if d.shape[1] == 1:
                        aligned[col_name] = d.iloc[:, 0]
                    else:
                        aligned[col_name] = d.iloc[:, -1]
                stats = pd.DataFrame()
                stats["mean"] = aligned.mean(axis=1)
                stats["std"] = aligned.std(axis=1)
                stats.index.name = "index"
                return stats.reset_index()
            else:
                # same length, use integer index
                aligned = pd.DataFrame({f"run_{i}": d.iloc[:, -1] if d.shape[1] > 1 else d.iloc[:, 0] for i, d in enumerate(dfs)})
                stats = pd.DataFrame()
                stats["mean"] = aligned.mean(axis=1)
                stats["std"] = aligned.std(axis=1)
                stats.index.name = "index"
                return stats.reset_index()

    # Align by time_col
    # Reindex each df to the union of time points using interpolation
    all_times = pd.Index(sorted(set().union(*[set(d[time_col]) for d in dfs])))
    reindexed = []
    for d in dfs:
        single = d.set_index(time_col)
        # pick the last numeric column as the value column if multiple
        val_col = single.select_dtypes(include=["number"]).columns[-1]
        series = single[val_col].reindex(all_times)
        # interpolate missing values
        series = series.interpolate().ffill().bfill()
        reindexed.append(series)
    stacked = pd.concat(reindexed, axis=1)
    stacked.columns = [f"run_{i}" for i in range(len(reindexed))]
    stats = pd.DataFrame({"time": all_times, "mean": stacked.mean(axis=1), "std": stacked.std(axis=1)})
    return stats


def plot_stats(df_stats: pd.DataFrame, out_plot: Path, time_col_name: str = "time") -> None:
    ensure_deps()
    plt.figure(figsize=(8, 5))
    x = df_stats[time_col_name]
    y = df_stats["mean"]
    yerr = df_stats["std"]
    plt.plot(x, y, label="mean")
    plt.fill_between(x, y - yerr, y + yerr, alpha=0.3, label="std")
    plt.xlabel(time_col_name)
    plt.ylabel("value")
    plt.legend()
    plt.tight_layout()
    tmp = out_plot.with_suffix(out_plot.suffix + ".tmp")
    plt.savefig(tmp)
    # atomic move
    os.replace(tmp, out_plot)
    plt.close()


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run many independent simulations in parallel and aggregate outputs")
    parser.add_argument("--executable", required=True, help="Path to simulation executable")
    parser.add_argument("--n", type=int, required=True, help="Number of independent simulations to run")
    parser.add_argument("--parallel", type=int, default=4, help="Max concurrent simulations")
    parser.add_argument("--seed-start", type=int, default=1, help="Starting seed (incremented by 1 for each run)")
    parser.add_argument("--out-arg", default="--out", help="Name of the output file argument for the executable")
    parser.add_argument("--seed-arg", default="--seed", help="Name of the seed argument for the executable")
    parser.add_argument("--exe-arg", action="append", default=[], help="Additional fixed args to pass to the executable (can be repeated)")
    parser.add_argument("--capture-stdout", action="store_true", help="If set, capture the executable's stdout and save as CSV")
    parser.add_argument("--timeout", type=float, default=300.0, help="Per-run timeout in seconds")
    parser.add_argument("--retries", type=int, default=1, help="Retries per failed run")
    parser.add_argument("--out-plot", default="combined.png", help="Final combined plot path")
    parser.add_argument("--out-csv", default="combined.csv", help="Final combined CSV path")
    parser.add_argument("--time-col", default=None, help="Name of the time column to align on (auto-detected if omitted)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary per-run outputs for debugging")
    args = parser.parse_args(argv)

    exe = os.path.abspath(args.executable)
    if not os.path.exists(exe):
        print(f"Executable not found: {exe}", file=sys.stderr)
        return 2

    # Prepare a temporary run directory
    run_dir = Path(tempfile.mkdtemp(prefix="sim_run_"))
    print(f"Run directory: {run_dir}")

    # prepare graceful shutdown on SIGINT/SIGTERM
    def _handle_sig(signum, frame):
        print("Received signal, terminating child processes...", file=sys.stderr)
        CHILD_MANAGER.terminate_all()

    signal.signal(signal.SIGINT, _handle_sig)
    signal.signal(signal.SIGTERM, _handle_sig)

    jobs = []
    for i in range(args.n):
        seed = args.seed_start + i
        out_path = Path(run_dir) / f"sim_{i:04d}.csv"
        jobs.append((exe, i, out_path, seed, args.exe_arg, args.out_arg, args.seed_arg, args.capture_stdout, args.timeout, args.retries))

    # Use a ThreadPoolExecutor to start subprocesses concurrently (subprocesses are separate OS processes)
    results = []
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.parallel) as ex:
            futs = [ex.submit(lambda job: run_single_simulation(*job), job) for job in jobs]
            for f in concurrent.futures.as_completed(futs):
                try:
                    res = f.result()
                except KeyboardInterrupt:
                    print("Interrupted by user", file=sys.stderr)
                    CHILD_MANAGER.terminate_all()
                    raise
                except Exception as e:
                    print("Run exception:", e, file=sys.stderr)
                    res = None
                if res is not None:
                    idx, attempts, outp, rc, stderr_text = res
                    print(f"Run {idx}: rc={rc} attempts={attempts}")
                    if rc != 0:
                        print(f"  STDERR: {stderr_text}", file=sys.stderr)
                    results.append(res)
    except KeyboardInterrupt:
        print("User aborted, cleaning up...", file=sys.stderr)
        CHILD_MANAGER.terminate_all()
        if not args.keep_temp:
            shutil.rmtree(run_dir, ignore_errors=True)
        return 130

    # collect produced CSV paths
    csv_paths = [r[2] for r in results if r[2] is not None and Path(r[2]).exists()]
    if len(csv_paths) == 0:
        print("No CSV outputs produced; aborting aggregation", file=sys.stderr)
        if not args.keep_temp:
            shutil.rmtree(run_dir, ignore_errors=True)
        return 3

    # Aggregate
    ensure_deps()
    csv_paths = [Path(p) for p in csv_paths]
    try:
        stats = aggregate_csvs(csv_paths, time_col=args.time_col)
    except Exception as e:
        print(f"Aggregation failed: {e}", file=sys.stderr)
        if not args.keep_temp:
            shutil.rmtree(run_dir, ignore_errors=True)
        return 4

    # Write combined CSV atomically
    out_csv_path = Path(args.out_csv).resolve()
    tmp_csv = out_csv_path.with_suffix(out_csv_path.suffix + ".tmp")
    stats.to_csv(tmp_csv, index=False)
    atomic_move(tmp_csv, out_csv_path)
    print(f"Wrote combined CSV: {out_csv_path}")

    # Plot
    out_plot_path = Path(args.out_plot).resolve()
    plot_stats(stats, out_plot_path, time_col_name=(args.time_col or stats.columns[0]))
    print(f"Wrote combined plot: {out_plot_path}")

    if not args.keep_temp:
        shutil.rmtree(run_dir, ignore_errors=True)
    else:
        print(f"Temporary run directory left at: {run_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
