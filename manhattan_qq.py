
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# GWAS Manhattan + QQ (publication-ready)
# Usage:
#   python code/manhattan_qq.py --in data/gwas_sumstats_template.csv --sep , --build GRCh38 --out out
#
# Required columns: SNP, CHR, BP, P

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True, help="Input CSV/TSV with SNP,CHR,BP,P")
    ap.add_argument("--sep", default=",", help="Separator: ',' for CSV, '\\t' for TSV (default: ',')")
    ap.add_argument("--build", default="GRCh38", help="Genome build label for title")
    ap.add_argument("--out", default=None, help="Output basename (default: derived from timestamp)")
    ap.add_argument("--include-chrs", default="1-22", help="Chromosomes to include (e.g., '1-22' or '1-22,23')")
    ap.add_argument("--sugg", type=float, default=1e-5, help="Suggestive p-value threshold")
    ap.add_argument("--gwas", type=float, default=5e-8, help="Genome-wide p-value threshold")
    ap.add_argument("--annotate-top", type=int, default=5, help="Annotate top-N hits by p-value (0=off)")
    ap.add_argument("--highlight", default=None, help="Optional file with SNP IDs to highlight (one per line)")
    args = ap.parse_args()

    INFILE = Path(args.infile)
    SEP = args.sep
    BUILD = args.build
    OUT_BASENAME = args.out or f"gwas_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    SUGG = args.sugg
    GWAS = args.gwas
    ANNOTATE_TOP_N = args.annotate_top

    # Parse chromosomes to include
    include = []
    for part in args.include_chrs.split(","):
        part = part.strip()
        if "-" in part:
            a,b = part.split("-")
            include.extend(range(int(a), int(b)+1))
        else:
            include.append(int(part))
    INCLUDE_CHRS = sorted(set(include))

    # Load data
    df = pd.read_csv(INFILE, sep=SEP)
    required = {"SNP", "CHR", "BP", "P"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {missing} in {INFILE}")

    df = df[list(required)].dropna().copy()
    df["CHR"] = (df["CHR"].astype(str).str.replace("X", "23", regex=False)).astype(float).astype(int)
    df = df[df["CHR"].isin(INCLUDE_CHRS)].copy()

    eps = np.nextafter(0, 1)
    df.loc[df["P"] <= 0, "P"] = eps

    # Prep positions
    df = df.sort_values(["CHR", "BP"]).reset_index(drop=True)
    chr_max_bp = df.groupby("CHR")["BP"].max().cumsum()
    offsets = {c: int(chr_max_bp.shift(fill_value=0).loc[c]) for c in chr_max_bp.index}
    df["pos"] = df.apply(lambda r: r["BP"] + offsets[r["CHR"]], axis=1)
    centers = df.groupby("CHR")["pos"].median()
    df["neglog10p"] = -np.log10(df["P"])

    # Read highlight list if provided
    highlight_snps = set()
    if args.highlight:
        hp = Path(args.highlight)
        if hp.exists():
            with open(hp) as f:
                for line in f:
                    s = line.strip()
                    if s:
                        highlight_snps.add(s)

    # Manhattan
    fig = plt.figure(figsize=(14, 6))
    ax = plt.gca()

    for i, chr_id in enumerate(centers.index, start=1):
        if i % 2 == 0:
            start = df.loc[df["CHR"] == chr_id, "pos"].min()
            end   = df.loc[df["CHR"] == chr_id, "pos"].max()
            ax.axvspan(start, end, alpha=0.05)

    for c, sub in df.groupby("CHR"):
        mask_sig  = sub["P"] < GWAS
        mask_sug  = (sub["P"] < SUGG) & ~mask_sig
        mask_rest = ~mask_sig & ~mask_sug

        ax.scatter(sub.loc[mask_rest, "pos"], sub.loc[mask_rest, "neglog10p"],
                   s=6, alpha=0.7, linewidths=0)
        ax.scatter(sub.loc[mask_sug, "pos"], sub.loc[mask_sug, "neglog10p"],
                   s=10, alpha=0.8, linewidths=0)
        ax.scatter(sub.loc[mask_sig, "pos"], sub.loc[mask_sig, "neglog10p"],
                   s=18, alpha=0.95, linewidths=0)

    # Optional highlight
    if highlight_snps:
        sub_hl = df[df["SNP"].astype(str).isin(highlight_snps)]
        if not sub_hl.empty:
            ax.scatter(sub_hl["pos"], sub_hl["neglog10p"], s=28, alpha=0.98, linewidths=0)

    ax.axhline(-np.log10(SUGG), linestyle="--", linewidth=1.5)
    ax.axhline(-np.log10(GWAS), linestyle="-",  linewidth=1.5)

    if ANNOTATE_TOP_N and ANNOTATE_TOP_N > 0:
        for _, r in df.nsmallest(ANNOTATE_TOP_N, "P").iterrows():
            ax.text(r["pos"], r["neglog10p"] + 0.2, str(r["SNP"]), fontsize=8,
                    ha="center", va="bottom")

    ymax = max(9, float(np.nanpercentile(df["neglog10p"], 99.9) + 0.8))
    ax.set_ylim(0, ymax)
    ax.set_xlim(df["pos"].min(), df["pos"].max())
    ax.set_xticks(centers.values)
    ax.set_xticklabels(centers.index, fontsize=10)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(f"Manhattan Plot (Build: {BUILD})", fontsize=15, weight="bold", pad=10)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    out_dir = INFILE.parent if INFILE.parent.name else Path(".")
    png_path = out_dir / f"{OUT_BASENAME}_manhattan.png"
    pdf_path = out_dir / f"{OUT_BASENAME}_manhattan.pdf"
    plt.tight_layout()
    plt.savefig(png_path, dpi=220)
    plt.savefig(pdf_path)
    plt.close()

    # QQ plot
    pvals = np.sort(df["P"].values)
    n = len(pvals)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(pvals)

    fig = plt.figure(figsize=(6, 6))
    ax = plt.gca()
    ax.scatter(exp, obs, s=6, alpha=0.8, linewidths=0)
    maxv = float(max(exp.max(), obs.max()))
    ax.plot([0, maxv], [0, maxv], linestyle="--")
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title("QQ Plot", fontsize=13, weight="bold")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(axis="both", linestyle=":", alpha=0.3)
    png_path2 = out_dir / f"{OUT_BASENAME}_qq.png"
    pdf_path2 = out_dir / f"{OUT_BASENAME}_qq.pdf"
    plt.tight_layout()
    plt.savefig(png_path2, dpi=220)
    plt.savefig(pdf_path2)
    plt.close()

    print("Saved:", png_path, pdf_path, png_path2, pdf_path2)

if __name__ == "__main__":
    main()
