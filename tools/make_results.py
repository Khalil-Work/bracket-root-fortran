#!/usr/bin/env python3
"""Regenerate figures, animation and the run summary from the Fortran CSVs.

Reads build/csv/<equation>_<method>.csv and build/csv/summary.csv (written by
examples/compare.f90 via `make run`) and produces:

  docs/results/convergence_comparison.png   bracket half-width vs iteration
  docs/results/residual_histories.png       |f(x)| vs iteration
  docs/results/bracket_animation.gif        bisection vs safeguarded secant
  docs/results/run_summary.json             compiler + metrics digest

The numerical solvers are pure Fortran; this script only visualises and
aggregates their recorded output.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import os
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_DIR = os.path.join(ROOT, "build", "csv")
OUT_DIR = os.path.join(ROOT, "docs", "results")

METHODS = ["bisection", "illinois", "safeguarded_secant"]
METHOD_LABEL = {
    "bisection": "bisection",
    "illinois": "Illinois",
    "safeguarded_secant": "safeguarded secant",
}
EQUATIONS = [
    ("cubic", r"$f(x) = x^3 - 2x - 5$"),
    ("cosx", r"$f(x) = \cos x - x$"),
    ("sinx", r"$f(x) = \sin x - \frac{1}{2}$"),
    ("omega", r"$f(x) = x\,e^{x} - 1$"),
]
CUBIC_ROOT = 2.0945514815423265

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.grid": False,
    "font.size": 9,
})


def read_history(equation: str, method: str) -> list[dict]:
    path = os.path.join(CSV_DIR, f"{equation}_{method}.csv")
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append({
                "iteration": int(row["iteration"]),
                "a": float(row["a"]),
                "b": float(row["b"]),
                "x": float(row["x"]),
                "fx": float(row["fx"]),
                "halfwidth": float(row["halfwidth"]),
            })
    return rows


def read_summary() -> list[dict]:
    with open(os.path.join(CSV_DIR, "summary.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


def make_still_figures() -> None:
    # --- bracket half-width vs iteration, log scale, no grid --------------
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0))
    for ax, (eq, label) in zip(axes.flat, EQUATIONS):
        for method in METHODS:
            rows = read_history(eq, method)
            k = [r["iteration"] for r in rows]
            w = [max(r["halfwidth"], 1e-18) for r in rows]
            n = rows[-1]["iteration"]
            ax.semilogy(k, w, marker="o", markersize=3, linewidth=1.2,
                        label=f"{METHOD_LABEL[method]}: {n} its")
        ax.set_title(label)
        ax.set_xlabel("iteration $k$")
        ax.set_ylabel("bracket half-width $(b_k-a_k)/2$")
        ax.grid(False)
        ax.legend(fontsize=8)
    fig.suptitle("Convergence comparison: bracket half-width per iteration", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "convergence_comparison.png"), dpi=150)
    plt.close(fig)

    # --- |f(x_k)| vs iteration, log scale, no grid ------------------------
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0))
    for ax, (eq, label) in zip(axes.flat, EQUATIONS):
        for method in METHODS:
            rows = read_history(eq, method)
            k = [r["iteration"] for r in rows]
            r_abs = [max(abs(r["fx"]), 1e-18) for r in rows]
            n = rows[-1]["iteration"]
            ax.semilogy(k, r_abs, marker="o", markersize=3, linewidth=1.2,
                        label=f"{METHOD_LABEL[method]}: {n} its")
        ax.set_title(label)
        ax.set_xlabel("iteration $k$")
        ax.set_ylabel(r"$|f(x_k)|$")
        ax.grid(False)
        ax.legend(fontsize=8)
    fig.suptitle("Residual histories: $|f(x_k)|$ per iteration", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "residual_histories.png"), dpi=150)
    plt.close(fig)


def cubic(x):
    return x ** 3 - 2 * x - 5


def draw_panel(ax, row, name, final):
    a, b, x, fx = row["a"], row["b"], row["x"], row["fx"]
    w = max(b - a, 1e-15)
    pad = 0.45 * w
    lo, hi = a - pad, b + pad
    xs = np.linspace(lo, hi, 800)
    ys = cubic(xs)
    ax.plot(xs, ys, color="#1f77b4", lw=1.4, label="$f(x)$")
    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.axvspan(a, b, color="tab:orange", alpha=0.15)
    ax.axvline(a, color="tab:orange", ls="--", lw=1.0)
    ax.axvline(b, color="tab:orange", ls="--", lw=1.0)
    ax.axvline(CUBIC_ROOT, color="tab:green", ls=":", lw=1.0)
    ax.plot([x], [fx], "o", color="tab:red", ms=5, zorder=5)
    ymin, ymax = float(ys.min()), float(ys.max())
    span = (ymax - ymin) or 1.0
    ax.set_ylim(ymin - 0.15 * span, ymax + 0.15 * span)
    ax.set_xlim(lo, hi)
    status = " (converged)" if final else ""
    ax.set_title(f"{name} -- iteration {row['iteration']}{status}", fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    ax.grid(False)
    ax.text(0.02, 0.04, f"(b-a)/2 = {row['halfwidth']:.2e}",
            transform=ax.transAxes, fontsize=8, va="bottom")


def render_frame(t, bis, saf):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), dpi=100)
    draw_panel(axes[0], bis[min(t, len(bis) - 1)], "bisection", t >= len(bis) - 1)
    draw_panel(axes[1], saf[min(t, len(saf) - 1)], "safeguarded secant", t >= len(saf) - 1)
    fig.suptitle("Bracket Root: bracket [a, b] on  x^3 - 2x - 5 = 0"
                 "   (root x* = 2.0945514815423265)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    img = Image.fromarray(buf[:, :, :3].copy())
    plt.close(fig)
    return img


def make_animation() -> int:
    bis = read_history("cubic", "bisection")
    saf = read_history("cubic", "safeguarded_secant")
    n = max(len(bis), len(saf))
    repeat = max(1, math.ceil(30 / n))          # ensure at least 30 frames
    frame_idx = [t for t in range(n) for _ in range(repeat)]
    if len(frame_idx) > 80:                     # honour the 80-frame cap
        step = math.ceil(len(frame_idx) / 80)
        frame_idx = frame_idx[::step][:80]

    frames = [render_frame(t, bis, saf) for t in frame_idx]
    gif_path = os.path.join(OUT_DIR, "bracket_animation.gif")
    frames[0].save(gif_path, save_all=True, append_images=frames[1:],
                   duration=400, loop=0)
    for f in frames:
        f.close()
    return len(frame_idx)


def write_run_summary(frames: int) -> None:
    fc = os.environ.get("FC", "gfortran")
    gfortran = shutil.which(fc)
    if gfortran:
        out = subprocess.run([gfortran, "--version"], capture_output=True, text=True)
        compiler = out.stdout.splitlines()[0].strip() if out.stdout else "unknown"
    else:
        compiler = "gfortran not found on PATH"
    payload = {
        "project": "bracket-root-fortran",
        "generated_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "compiler": compiler,
        "fortran_flags": "-std=f2018 -Wall -Wextra -O2 -fall-intrinsics",
        "driver_tolerances": {"xtol": 1.0e-14, "ftol": 1.0e-14, "maxiter": 200},
        "animation": {"file": "docs/results/bracket_animation.gif",
                      "frames": frames, "frame_ms": 400},
        "results": [
            {
                "method": r["method"],
                "equation": r["equation"],
                "iterations": int(r["iterations"]),
                "final_residual": float(r["final_residual"]),
                "flag": int(r["flag"]),
                "flag_name": r["flag_name"],
                "root": float(r["root"]),
                "reference_root": float(r["reference_root"]),
                "abs_error": float(r["abs_error"]),
            }
            for r in read_summary()
        ],
    }
    with open(os.path.join(OUT_DIR, "run_summary.json"), "w") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    make_still_figures()
    frames = make_animation()
    write_run_summary(frames)
    gif = os.path.join(OUT_DIR, "bracket_animation.gif")
    with Image.open(gif) as im:
        print(f"animation: {im.n_frames} frames, {im.size[0]}x{im.size[1]} px")
    print(f"wrote figures and run_summary.json to {OUT_DIR}")


if __name__ == "__main__":
    main()
