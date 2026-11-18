"""Run the heat-equity analysis and write a short report with figures.

    python scripts/run_analysis.py \
        --heat data/heat_greenness_tracts_GA_059.csv \
        --acs data/acs_GA_059.csv \
        --out results/

Inputs:
  --heat  CSV exported by gee/heat_greenness_by_tract.js
  --acs   CSV written by scripts/get_acs.py
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import heatequity as he  # noqa: E402

INDICATORS = {"pct_poverty": 1, "pct_minority": 1, "pct_age65": 1, "median_income": -1}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--heat", required=True, type=Path)
    p.add_argument("--acs", required=True, type=Path)
    p.add_argument("--out", default=Path("results"), type=Path)
    p.add_argument("--k", type=int, default=6, help="neighbours for Moran's I")
    a = p.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    df = he.load_and_join(a.heat, a.acs)
    df["vulnerability"] = he.vulnerability_index(df, INDICATORS)
    df["tree_pct"] = 100 * df["tree_frac"]
    df = df.dropna(subset=["vulnerability", "LST_C", "NDVI", "tree_pct"])

    q = he.quintile_table(df, "vulnerability", ["LST_C", "NDVI", "tree_pct", "pct_poverty"])
    gap = q.loc["Q5", "LST_C"] - q.loc["Q1", "LST_C"]
    gap_words = f"{abs(gap):.1f} °C {'hotter' if gap >= 0 else 'cooler'}"

    model = he.ols(df, "LST_C", ["tree_pct", "vulnerability"])
    w = he.knn_weights(he.lonlat_to_km(df["lon"], df["lat"]), k=min(a.k, len(df) - 1))
    i_lst = he.morans_i(df["LST_C"], w)
    i_res = he.morans_i(model.residuals, w)

    # ---- figures
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.bar(q.index, q["LST_C"], color=plt.cm.YlOrRd(np.linspace(0.3, 0.9, 5)))
    ax.set_ylim(q["LST_C"].min() - 1.5, q["LST_C"].max() + 0.8)
    ax.set_xlabel("Social vulnerability quintile (Q1 = least, Q5 = most)")
    ax.set_ylabel("Mean summer land surface temperature (°C)")
    ax.set_title(f"Most vulnerable tracts are {gap_words}\nthan the least vulnerable")
    for x, v in zip(q.index, q["LST_C"]):
        ax.text(x, v + 0.1, f"{v:.1f}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(a.out / "lst_by_vulnerability.png", dpi=150)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    sc = ax.scatter(df["tree_pct"], df["LST_C"], c=df["vulnerability"], cmap="magma_r", s=28, edgecolor="k", linewidth=0.3)
    xs = np.linspace(df["tree_pct"].min(), df["tree_pct"].max(), 50)
    ys = model.params["const"] + model.params["tree_pct"] * xs + model.params["vulnerability"] * df["vulnerability"].mean()
    ax.plot(xs, ys, color="k", lw=1)
    ax.set_xlabel("Tree cover (% of tract)")
    ax.set_ylabel("Summer land surface temperature (°C)")
    fig.colorbar(sc, label="Vulnerability index")
    ax.set_title(f"Each +10 pts tree cover: {10 * model.params['tree_pct']:+.2f} °C")
    fig.tight_layout()
    fig.savefig(a.out / "lst_vs_tree_cover.png", dpi=150)

    # ---- report
    s = model.summary().round(3)
    report = f"""# Results

Tracts analysed: **{len(df)}** (population {df['population'].sum():,.0f})

## Heat and greenness by social vulnerability quintile

{q.round(2).to_markdown()}

The most vulnerable fifth of tracts is **{gap_words}** than the least vulnerable fifth
(population-weighted mean summer land surface temperature).

## Regression: LST ~ tree cover + vulnerability (HC1 robust SE)

{s.to_markdown()}

R² = {model.r_squared:.3f}, n = {model.n}

## Spatial autocorrelation (Moran's I, {len(df)} tracts, k = {min(a.k, len(df) - 1)} nearest neighbours)

| Variable | I | E[I] | pseudo p |
| --- | --- | --- | --- |
| Land surface temperature | {i_lst[0]:.3f} | {i_lst[1]:.3f} | {i_lst[2]:.3f} |
| Regression residuals | {i_res[0]:.3f} | {i_res[1]:.3f} | {i_res[2]:.3f} |

Significant clustering in the residuals means a spatial regression model
(spatial lag or error) would be more appropriate than OLS.

![](lst_by_vulnerability.png)
![](lst_vs_tree_cover.png)
"""
    (a.out / "results.md").write_text(report)
    df.to_csv(a.out / "tracts_joined.csv")
    print(report)


if __name__ == "__main__":
    main()
