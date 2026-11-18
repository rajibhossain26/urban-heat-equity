import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import heatequity as he

ROOT = Path(__file__).resolve().parents[1]


def synthetic_county(n_side=12, seed=1):
    """Tracts on a grid where heat follows a smooth spatial gradient, trees
    cool, and vulnerability is higher in the hot east."""
    rng = np.random.default_rng(seed)
    gx, gy = np.meshgrid(np.arange(n_side), np.arange(n_side))
    lon = -83.5 + gx.ravel() * 0.02
    lat = 33.9 + gy.ravel() * 0.02
    east = gx.ravel() / (n_side - 1)
    n = n_side**2
    tree = np.clip(0.6 - 0.4 * east + rng.normal(0, 0.05, n), 0, 1)
    poverty = np.clip(10 + 25 * east + rng.normal(0, 3, n), 0, 100)
    lst = 38 - 8 * tree + 2 * east + rng.normal(0, 0.3, n)
    geoid = [f"13059{i:06d}" for i in range(n)]
    heat = pd.DataFrame({"GEOID": geoid, "lon": lon, "lat": lat, "LST_C": lst,
                         "NDVI": 0.2 + 0.6 * tree, "tree_frac": tree, "built_frac": 1 - tree})
    acs = pd.DataFrame({"GEOID": geoid, "population": rng.integers(800, 6000, n),
                        "median_income": 80000 - 900 * poverty, "pct_poverty": poverty,
                        "pct_minority": np.clip(20 + 40 * east, 0, 100), "pct_age65": 15 + rng.normal(0, 2, n)})
    return heat, acs


def test_load_and_join_pads_geoids_and_filters(tmp_path):
    heat, acs = synthetic_county()
    heat.loc[0, "GEOID"] = heat.loc[0, "GEOID"].lstrip("1")  # lost leading digits
    acs.loc[0, "GEOID"] = acs.loc[0, "GEOID"].lstrip("1")
    acs.loc[1, "population"] = 10
    heat.loc[2, "LST_C"] = np.nan
    heat.to_csv(tmp_path / "h.csv", index=False)
    acs.to_csv(tmp_path / "a.csv", index=False)
    df = he.load_and_join(tmp_path / "h.csv", tmp_path / "a.csv")
    assert len(df) == len(heat) - 2
    assert all(len(g) == 11 for g in df.index)


def test_vulnerability_index_direction():
    df = pd.DataFrame({"pov": [0.0, 10.0, 20.0], "inc": [90.0, 60.0, 30.0]})
    v = he.vulnerability_index(df, {"pov": 1, "inc": -1})
    assert v.is_monotonic_increasing
    assert v.mean() == pytest.approx(0.0)
    with pytest.raises(ValueError):
        he.vulnerability_index(df, {"pov": 2})


def test_quintiles_weighted_means():
    df = pd.DataFrame({"v": np.arange(10.0), "y": np.arange(10.0), "population": [1, 3] * 5})
    q = he.quintile_table(df, "v", ["y"])
    assert list(q.index) == ["Q1", "Q2", "Q3", "Q4", "Q5"]
    assert q.loc["Q1", "y"] == pytest.approx(0.75)  # (0*1 + 1*3) / 4
    assert q["tracts"].sum() == 10


def test_ols_recovers_coefficients():
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"a": rng.normal(size=500), "b": rng.normal(size=500)})
    df["y"] = 2 + 3 * df.a - 1.5 * df.b + rng.normal(0, 0.1, 500)
    r = he.ols(df, "y", ["a", "b"])
    assert r.params.round(1).tolist() == [2.0, 3.0, -1.5]
    assert r.r_squared > 0.99
    assert (r.std_err > 0).all()
    classic = he.ols(df, "y", ["a", "b"], robust=False)
    assert classic.params.equals(r.params)


def test_knn_weights_row_standardized():
    w = he.knn_weights(np.array([[0, 0], [1, 0], [0, 1], [5, 5]]), k=2)
    assert np.allclose(w.sum(axis=1), 1)
    assert np.all(np.diag(w) == 0)
    assert w[0, 1] == 0.5 and w[0, 2] == 0.5


def test_morans_i_detects_clustering_and_randomness():
    n_side = 15
    gx, gy = np.meshgrid(np.arange(n_side), np.arange(n_side))
    coords = np.column_stack([gx.ravel(), gy.ravel()]).astype(float)
    w = he.knn_weights(coords, k=4)
    smooth = gx.ravel() + gy.ravel()
    i, e, p = he.morans_i(smooth, w, permutations=199)
    assert i > 0.8 and p < 0.01 and e == pytest.approx(-1 / 224)
    noise = np.random.default_rng(3).normal(size=n_side**2)
    i2, _, p2 = he.morans_i(noise, w, permutations=199)
    assert abs(i2) < 0.15 and p2 > 0.01


def test_lonlat_to_km_scale():
    xy = he.lonlat_to_km([0, 0], [0, 1])
    assert np.hypot(*(xy[1] - xy[0])) == pytest.approx(111.19, abs=0.05)


def test_run_analysis_end_to_end(tmp_path):
    heat, acs = synthetic_county()
    heat.to_csv(tmp_path / "heat.csv", index=False)
    acs.to_csv(tmp_path / "acs.csv", index=False)
    out = tmp_path / "results"
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "run_analysis.py"),
         "--heat", str(tmp_path / "heat.csv"), "--acs", str(tmp_path / "acs.csv"), "--out", str(out)],
        check=True, capture_output=True, text=True,
    )
    report = (out / "results.md").read_text()
    assert "Tracts analysed: **144**" in report
    assert (out / "lst_by_vulnerability.png").stat().st_size > 10_000
    joined = pd.read_csv(out / "tracts_joined.csv")
    # the synthetic county is built so vulnerable tracts are hotter
    q = he.quintile_table(joined, "vulnerability", ["LST_C"])
    assert q.loc["Q5", "LST_C"] > q.loc["Q1", "LST_C"]


def test_acs_tidy_handles_missing_codes():
    sys.path.insert(0, str(ROOT / "scripts"))
    import get_acs

    rec = {v: "100" for v in get_acs.VARIABLES}
    rec.update({"B19013_001E": "-666666666", "B17001_002E": "25", "B03002_003E": "60",
                "state": "13", "county": "059", "tract": "000100"})
    for v in get_acs.AGE_65_PLUS:
        rec[v] = "1"
    row = get_acs.tidy([rec])[0]
    assert row["GEOID"] == "13059000100"
    assert row["median_income"] is None
    assert row["pct_poverty"] == 25.0
    assert row["pct_minority"] == 40.0
    assert row["pct_age65"] == 12.0


def test_summary_file_parsing_matches_api_layout():
    sys.path.insert(0, str(ROOT / "scripts"))
    import get_acs

    assert get_acs.summary_column("B01001_020E") == "B01001_E020"
    lines = [
        "GEO_ID|B17001_E001|B17001_M001|B17001_E002|B17001_M002\n",
        "0500000US13059|120000|50|30000|900\n",
        "1400000US13057090100|1000|10|100|5\n",
        "1400000US13059000100|800|20|200|9\n",
        "1400000US13059000200|500|15|-666666666|-222222222\n",
        "1400000US13061950100|400|10|40|4\n",
    ]
    recs = get_acs.parse_summary_lines(lines, "13", "059", ["B17001_001E", "B17001_002E"])
    assert sorted(recs) == ["13059000100", "13059000200"]
    assert recs["13059000100"] == {"state": "13", "county": "059", "tract": "000100",
                                   "B17001_001E": "800", "B17001_002E": "200"}
    rows = get_acs.tidy(recs.values())
    assert rows[0]["pct_poverty"] == 25.0
    assert rows[1]["pct_poverty"] is None
