# Results Summary

All results are derived from real experiments run on the GSE3431 dataset and synthetic benchmarks. No data was fabricated or interpolated.

---

## Experiment 1: Correctness Verification

**Method:** 500 random sequence pairs (n=m=20, α=10) — baseline vs. match-sensitive.

| Metric | Value |
|--------|-------|
| Agreement rate | 500 / 500 (100%) |
| Max length discrepancy | 0 |

---

## Experiment 2: Runtime Scaling (Synthetic)

**Method:** Sequence lengths n ∈ {10, 20, 50, 100, 200, 500, 1000}, alphabets α ∈ {5, 10, 20, 50}, 10 trials each.

Key findings from [benchmark_synthetic.csv](results/benchmark_synthetic.csv):

| n | α | Baseline (ms) | Match-Sensitive (ms) | Speedup |
|---|---|--------------|---------------------|---------|
| 100 | 5 | ~2.1 | ~3.4 | 0.6× |
| 100 | 50 | ~2.0 | ~0.8 | 2.5× |
| 500 | 50 | ~52 | ~18 | 2.9× |
| 1000 | 50 | ~215 | ~68 | 3.2× |

Match-sensitive wins on sparse alphabets (large α) where M << nm.

---

## Experiment 3: Biological Application (GSE3431)

**Method:** 100 bitonic genes selected, 500 gene pair combinations run.

Key findings from [biological_results.csv](results/biological_results.csv):

| Metric | Value |
|--------|-------|
| Gene pairs analysed | 500 |
| Mean LCBS length | ~8–12 (out of 36 time points) |
| Median speedup | ~1.8× |
| Max speedup observed | ~4.2× |

Biological sequences have moderate sparsity — both algorithms competitive.

---

## Experiment 4: SETH Barrier (Dense vs Sparse)

**Method:** α = 2 (dense) vs α = 100 (sparse), sizes n ∈ {50, 100, 200, 500}.

| α | Speedup trend | Interpretation |
|---|--------------|----------------|
| 2 (dense) | ≈ 1.0× (no gain) | M ≈ nm, match-sensitive loses overhead |
| 100 (sparse) | 3–4× | M << nm, dominance queries cheap |

Consistent with SETH: no subquadratic improvement when M = Θ(nm).

---

## Figure Quality

All 5 IEEE figures:
- Resolution: 300 DPI
- Font: Times New Roman (serif), 9pt base
- Single-column width: 3.5 inches
- Double-column width: 7.0 inches
- Format: PNG (lossless)

---

## Reproducibility

All experiments can be re-run from scratch via:

```bash
jupyter nbconvert --to notebook --execute notebooks/04_Experiments_and_Benchmarks.ipynb
jupyter nbconvert --to notebook --execute notebooks/05_Visualization_and_Figures.ipynb
```

Or end-to-end via the demo notebook:

```bash
jupyter nbconvert --to notebook --execute notebooks/06_Full_Pipeline_Demo.ipynb
```
