[README.md](https://github.com/user-attachments/files/22732414/README.md)
# GWAS Manhattan + QQ (Python)

This repo contains a **publication-ready Manhattan plot** and **QQ plot** script for GWAS summary statistics.

## What you get
- `code/manhattan_qq.py`: Python script to generate Manhattan + QQ plots
- `data/gwas_sumstats_template.csv`: Example summary stats (SNP, CHR, BP, P)
- Output figures: PNG + vector PDF

## Requirements
- Python 3.8+
- Packages: `pandas`, `numpy`, `matplotlib`

Install:
```bash
pip install pandas numpy matplotlib
```

## Usage
```bash
python code/manhattan_qq.py --in data/gwas_sumstats_template.csv --sep , --build GRCh38 --out out
```

Options:
- `--in`: input CSV/TSV path with columns `SNP,CHR,BP,P`
- `--sep`: separator (`,` for CSV, `\t` for TSV)
- `--build`: genome build label for the figure title (e.g., GRCh37/GRCh38)
- `--out`: output basename (default: timestamp)
- `--include-chrs`: chromosomes to include (default `1-22`, can be `1-22,23` to include X)
- `--sugg`: suggestive p-value threshold (default `1e-5`)
- `--gwas`: genome-wide p-value threshold (default `5e-8`)
- `--annotate-top`: annotate top-N SNPs by p-value (default `5`)
- `--highlight`: file with SNPs to highlight (one per line)

## Example
```bash
python code/manhattan_qq.py --in data/gwas_sumstats_template.csv --sep , --build GRCh38 --out demo
```

This will produce:
- `demo_manhattan.png` and `.pdf`
- `demo_qq.png` and `.pdf`

## How to publish
1. Create a new repo on GitHub (e.g., `gwas-manhattan-qq`).
2. Upload this folder (or push via git):
```bash
git init
git add .
git commit -m "Add GWAS Manhattan + QQ plotting kit"
git branch -M main
git remote add origin https://github.com/YourUser/gwas-manhattan-qq.git
git push -u origin main
```
3. Add a short LinkedIn post describing the application and link to the repo.

## License
MIT
