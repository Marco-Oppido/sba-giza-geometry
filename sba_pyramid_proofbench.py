#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GizaPyramids-SBA proofbench (streamlined)

README (reproducibility)
- Command:
  python3 sba_theory_proofbench.py --outdir outputs_giza_only --samples 50000
- Expected SHA-256 (submission tracking):
  giza_proofbank_report.txt  = 3d101f2d32475de7afe3ceb95c897a8cc31b7327ab3457c3903fbbdeb97507d8
  giza_proofbank_results.json = ffde4eed07201de5174bd60102cb65417b6e08e37500e0e06f05ed4cb0f451e8
  giza_proofbank_table.csv   = 5c7c83fe9ae8d9340a115ea9ddd0ac802517910ca97d7dd3a505c7293dbc5c00
  parameters.json            = 6f96619241c21bec0a49e547f0805de65301b52dafea5d1fd49bd2f70be678e7

This script runs only the pyramid-theory validation pipeline and writes:
  - giza_proofbank_table.csv
  - giza_proofbank_results.json
  - giza_proofbank_report.txt
  - parameters.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import math
import os
import sys
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np

OUTPUT_DIR = "outputs"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("giza_proofbench")


@dataclass
class Config:
    OUTPUT_DIR: str = OUTPUT_DIR
    RUN_ID: str = ""

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def ensure_output_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_csv_dicts(path: str, rows: List[Dict[str, object]]) -> None:
    if not rows:
        with open(path, "w", encoding="utf-8", newline=""):
            return
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)



def generate_required_figures(result: Dict[str, object], figures_dir: str = "figures") -> None:
    """Generate the two LaTeX-required figures using matplotlib."""
    os.makedirs(figures_dir, exist_ok=True)

    ratio_table = result["ratio_table"]
    monte_carlo = result["monte_carlo"]
    sqrt5 = float(result["sqrt5"])
    phi5 = float(result["phi5"])

    labels = ["User values", "Published values"]
    base_vals = [
        float(ratio_table[0]["r_base_khufu_menkaure"]),
        float(ratio_table[1]["r_base_khufu_menkaure"]),
    ]
    vol_vals = [
        float(ratio_table[0]["r_volume_khufu_menkaure"]),
        float(ratio_table[1]["r_volume_khufu_menkaure"]),
    ]

    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(6.6, 3.6), dpi=220)
    ax.bar(x - width / 2, base_vals, width, label="Base ratio Khufu/Menkaure", color="#4C78A8")
    ax.bar(x + width / 2, vol_vals, width, label="Volume ratio Khufu/Menkaure", color="#F58518")
    ax.axhline(sqrt5, color="#4C78A8", linestyle="--", linewidth=1.4, label=r"Target $\sqrt{5}$")
    ax.axhline(phi5, color="#F58518", linestyle="--", linewidth=1.4, label=r"Target $\varphi^5$")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Ratio")
    ax.set_title("Giza endpoint ratios vs SBA targets")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=7, loc="upper left", ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(figures_dir, "giza_endpoint_ratios.png"), bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 3.6), dpi=220)
    y_positions = [1.0, 0.0]
    means = [float(monte_carlo["base_ratio_mean"]), float(monte_carlo["volume_ratio_mean"])]
    q025 = [float(monte_carlo["base_ratio_q025"]), float(monte_carlo["volume_ratio_q025"])]
    q975 = [float(monte_carlo["base_ratio_q975"]), float(monte_carlo["volume_ratio_q975"])]

    for i, (mean, low, high) in enumerate(zip(means, q025, q975)):
        ax.plot([low, high], [y_positions[i], y_positions[i]], color="black", linewidth=2)
        ax.plot(mean, y_positions[i], "o", color="#2CA02C", markersize=6)

    rng = np.random.default_rng(20260328)
    base_sigma = (q975[0] - q025[0]) / (2 * 1.96)
    vol_sigma = (q975[1] - q025[1]) / (2 * 1.96)
    base_samples = rng.normal(means[0], base_sigma, 120)
    vol_samples = rng.normal(means[1], vol_sigma, 120)
    ax.scatter(base_samples, np.full_like(base_samples, 1.0) + rng.normal(0, 0.02, 120), s=8, alpha=0.25, color="#1f77b4")
    ax.scatter(vol_samples, np.full_like(vol_samples, 0.0) + rng.normal(0, 0.02, 120), s=8, alpha=0.25, color="#d62728")

    ax.axvline(sqrt5, color="#4C78A8", linestyle="--", linewidth=1.3, label=r"$\sqrt{5}$")
    ax.axvline(phi5, color="#F58518", linestyle="--", linewidth=1.3, label=r"$\varphi^5$")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(["Base ratio", "Volume ratio"])
    ax.set_xlabel("Ratio")
    ax.set_title("Monte Carlo intervals and residual scatter proxy")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(figures_dir, "giza_monte_carlo_intervals.png"), bbox_inches="tight")
    plt.close(fig)

def verify_giza_proofbank(outdir: str, mc_samples: int = 50000) -> Dict[str, object]:
    """Run Giza ratio checks and Monte Carlo uncertainty diagnostics."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    sqrt5 = math.sqrt(5.0)
    phi5 = phi**5

    datasets: List[Dict[str, float]] = [
        {
            "dataset": "user_working_values",
            "b_menkaure": 104.6,
            "b_khafre": 215.3,
            "b_khufu": 230.3,
            "v_menkaure": 237000.0,
            "v_khafre": 2211000.0,
            "v_khufu": 2592000.0,
        },
        {
            "dataset": "published_common_values",
            "b_menkaure": (104.6 + 102.2) / 2.0,
            "b_khafre": 215.25,
            "b_khufu": 230.33,
            "v_menkaure": 235183.0,
            "v_khafre": 2211096.0,
            "v_khufu": 2600000.0,
        },
    ]

    ratio_table: List[Dict[str, object]] = []
    for row in datasets:
        r_base = row["b_khufu"] / row["b_menkaure"]
        r_vol = row["v_khufu"] / row["v_menkaure"]
        ratio_table.append(
            {
                "dataset": row["dataset"],
                "r_base_khufu_menkaure": r_base,
                "target_sqrt5": sqrt5,
                "base_relative_error_pct": 100.0 * abs(r_base - sqrt5) / sqrt5,
                "r_volume_khufu_menkaure": r_vol,
                "target_phi5": phi5,
                "volume_relative_error_pct": 100.0 * abs(r_vol - phi5) / phi5,
                "r_base_khafre_menkaure": row["b_khafre"] / row["b_menkaure"],
                "r_base_khufu_khafre": row["b_khufu"] / row["b_khafre"],
                "r_volume_khafre_menkaure": row["v_khafre"] / row["v_menkaure"],
                "r_volume_khufu_khafre": row["v_khufu"] / row["v_khafre"],
            }
        )

    rng = np.random.default_rng(20260328)
    b_m = rng.uniform(102.2, 104.6, size=mc_samples)
    b_kh = rng.normal(230.33, 0.20, size=mc_samples)
    v_m = rng.normal(235183.0, 3500.0, size=mc_samples)
    v_kh = rng.normal(2600000.0, 30000.0, size=mc_samples)

    base_ratio = b_kh / b_m
    volume_ratio = v_kh / v_m

    base_q025 = float(np.quantile(base_ratio, 0.025))
    base_q975 = float(np.quantile(base_ratio, 0.975))
    vol_q025 = float(np.quantile(volume_ratio, 0.025))
    vol_q975 = float(np.quantile(volume_ratio, 0.975))

    monte_carlo = {
        "samples": int(mc_samples),
        "base_ratio_mean": float(np.mean(base_ratio)),
        "base_ratio_q025": base_q025,
        "base_ratio_q975": base_q975,
        "volume_ratio_mean": float(np.mean(volume_ratio)),
        "volume_ratio_q025": vol_q025,
        "volume_ratio_q975": vol_q975,
        "sqrt5_inside_base_interval": int(base_q025 <= sqrt5 <= base_q975),
        "phi5_inside_volume_interval": int(vol_q025 <= phi5 <= vol_q975),
    }

    result: Dict[str, object] = {
        "phi": phi,
        "sqrt5": sqrt5,
        "phi5": phi5,
        "ratio_table": ratio_table,
        "monte_carlo": monte_carlo,
    }

    write_csv_dicts(os.path.join(outdir, "giza_proofbank_table.csv"), ratio_table)
    with open(os.path.join(outdir, "giza_proofbank_results.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    lines = [
        "Pyramids of Giza proofbank",
        "Target tests: Khufu/Menkaure base ~ sqrt(5), Khufu/Menkaure volume ~ phi^5",
        f"phi={phi:.15f}, sqrt5={sqrt5:.15f}, phi^5={phi5:.15f}",
        "",
        "Deterministic ratio checks:",
    ]
    for row in ratio_table:
        lines.append(
            f"[{row['dataset']}] R_base={float(row['r_base_khufu_menkaure']):.6f} "
            f"(err={float(row['base_relative_error_pct']):.3f}%), "
            f"R_vol={float(row['r_volume_khufu_menkaure']):.6f} "
            f"(err={float(row['volume_relative_error_pct']):.3f}%)"
        )
        lines.append(
            f"   Khafre bridge: base(Khafre/Menkaure)={float(row['r_base_khafre_menkaure']):.6f}, "
            f"base(Khufu/Khafre)={float(row['r_base_khufu_khafre']):.6f}, "
            f"vol(Khafre/Menkaure)={float(row['r_volume_khafre_menkaure']):.6f}, "
            f"vol(Khufu/Khafre)={float(row['r_volume_khufu_khafre']):.6f}"
        )

    lines.extend(
        [
            "",
            "Uncertainty envelope (Monte Carlo):",
            f"samples={monte_carlo['samples']}",
            f"Base ratio 95% interval=[{base_q025:.6f}, {base_q975:.6f}], mean={float(monte_carlo['base_ratio_mean']):.6f}",
            f"Volume ratio 95% interval=[{vol_q025:.6f}, {vol_q975:.6f}], mean={float(monte_carlo['volume_ratio_mean']):.6f}",
            f"sqrt5 inside base interval: {int(monte_carlo['sqrt5_inside_base_interval'])}",
            f"phi^5 inside volume interval: {int(monte_carlo['phi5_inside_volume_interval'])}",
            "",
            "Interpretation: endpoint compatibility can be quantified; causal SBA execution mapping remains hypothesis-layer.",
            "Artifacts: giza_proofbank_table.csv, giza_proofbank_results.json, giza_proofbank_report.txt",
        ]
    )

    with open(os.path.join(outdir, "giza_proofbank_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    return result


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="GizaPyramids-SBA proofbench (streamlined): runs only Giza pyramid-theory validation."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=OUTPUT_DIR,
        help="Output directory for Giza proofbank artifacts.",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=50000,
        help="Monte Carlo sample count for uncertainty intervals (default: 50000).",
    )
    args = parser.parse_args(argv)

    config = Config()
    config.OUTPUT_DIR = args.outdir
    config.RUN_ID = hashlib.sha256(f"giza|{args.samples}|{config.OUTPUT_DIR}".encode("utf-8")).hexdigest()[:12]

    ensure_output_dir(config.OUTPUT_DIR)
    logger.info("Proofbench output directory: %s", config.OUTPUT_DIR)

    with open(os.path.join(config.OUTPUT_DIR, "parameters.json"), "w", encoding="utf-8") as f:
        json.dump(config.to_dict(), f, indent=2)
    logger.info("Wrote parameters.json")

    logger.info("Running streamlined Giza-only proofbank validations...")
    result = verify_giza_proofbank(outdir=config.OUTPUT_DIR, mc_samples=args.samples)
    generate_required_figures(result, figures_dir="figures")
    logger.info("Giza proofbank artifacts written to giza_proofbank_* files")
    logger.info("Required LaTeX figures written to figures/giza_endpoint_ratios.png and figures/giza_monte_carlo_intervals.png")
    logger.info("Proofbench run completed. Inspect %s for artifacts.", config.OUTPUT_DIR)


if __name__ == "__main__":
    main()