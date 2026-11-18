"""Who lives in the hottest, least green neighborhoods?

Tools to join tract-level heat and greenness (from Google Earth Engine) with
census vulnerability indicators and test for disparities, including spatial
autocorrelation. Only numpy and pandas are required.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


def _geoid(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(11)


def load_and_join(heat_csv, acs_csv, min_population: float = 100) -> pd.DataFrame:
    """Read both CSVs, join on the 11-digit tract GEOID, and drop tracts with
    very few residents or no satellite data."""
    heat = pd.read_csv(heat_csv, dtype={"GEOID": str})
    acs = pd.read_csv(acs_csv, dtype={"GEOID": str})
    heat["GEOID"], acs["GEOID"] = _geoid(heat["GEOID"]), _geoid(acs["GEOID"])
    df = heat.merge(acs, on="GEOID", how="inner", validate="one_to_one")
    df = df[(df["population"] >= min_population) & df["LST_C"].notna()].copy()
    return df.set_index("GEOID")


def vulnerability_index(df: pd.DataFrame, indicators: Mapping[str, int]) -> pd.Series:
    """Mean z-score of the indicators. ``indicators`` maps a column to +1 if a
    higher value means more vulnerable (e.g. pct_poverty) or -1 if lower means
    more vulnerable (e.g. median_income). Missing values are ignored per tract."""
    z = []
    for col, sign in indicators.items():
        if sign not in (1, -1):
            raise ValueError("indicator direction must be +1 or -1")
        x = df[col].astype(float)
        sd = x.std(ddof=0)
        z.append(sign * (x - x.mean()) / (sd if sd > 0 else 1.0))
    return pd.concat(z, axis=1).mean(axis=1, skipna=True).rename("vulnerability")


def quintile_table(df: pd.DataFrame, by: str, columns: Sequence[str], weight: str | None = "population") -> pd.DataFrame:
    """Population-weighted mean of ``columns`` in each quintile of ``by``
    (Q1 = lowest, Q5 = highest)."""
    q = pd.qcut(df[by].rank(method="first"), 5, labels=[f"Q{i}" for i in range(1, 6)])
    w = df[weight] if weight else pd.Series(1.0, index=df.index)
    rows = {}
    for label, idx in df.groupby(q, observed=True).groups.items():
        sub, ww = df.loc[idx], w.loc[idx]
        rows[label] = {c: np.average(sub[c], weights=ww) for c in columns}
        rows[label]["tracts"] = len(idx)
        rows[label]["population"] = ww.sum() if weight else np.nan
    return pd.DataFrame(rows).T


@dataclass
class OLSResult:
    params: pd.Series
    std_err: pd.Series
    t_values: pd.Series
    r_squared: float
    n: int
    residuals: np.ndarray

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame({"coef": self.params, "std_err": self.std_err, "t": self.t_values})


def ols(df: pd.DataFrame, y: str, x: Sequence[str], robust: bool = True) -> OLSResult:
    """Ordinary least squares with an intercept. ``robust`` uses HC1
    heteroskedasticity-consistent standard errors."""
    data = df[[y, *x]].dropna()
    Y = data[y].to_numpy(float)
    X = np.column_stack([np.ones(len(data)), data[list(x)].to_numpy(float)])
    n, k = X.shape
    if n <= k:
        raise ValueError("not enough observations")
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    xtx_inv = np.linalg.pinv(X.T @ X)
    if robust:
        meat = X.T @ (X * resid[:, None] ** 2)
        cov = xtx_inv @ meat @ xtx_inv * n / (n - k)
    else:
        cov = xtx_inv * (resid @ resid) / (n - k)
    se = np.sqrt(np.diag(cov))
    names = ["const", *x]
    r2 = 1 - (resid @ resid) / np.sum((Y - Y.mean()) ** 2)
    return OLSResult(
        params=pd.Series(beta, index=names),
        std_err=pd.Series(se, index=names),
        t_values=pd.Series(beta / se, index=names),
        r_squared=float(r2),
        n=n,
        residuals=resid,
    )


def knn_weights(coords: np.ndarray, k: int = 6) -> np.ndarray:
    """Row-standardized k-nearest-neighbour spatial weights (dense, for the
    few hundred tracts of a county or metro)."""
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if not 0 < k < n:
        raise ValueError("k must be between 1 and n - 1")
    d = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    nn = np.argsort(d, axis=1)[:, :k]
    w = np.zeros((n, n))
    w[np.repeat(np.arange(n), k), nn.ravel()] = 1.0 / k
    return w


def lonlat_to_km(lon, lat) -> np.ndarray:
    """Equirectangular projection to km around the data's mean latitude:
    accurate enough for neighbour search within a metro area."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    x = np.radians(lon) * 6371.0 * np.cos(np.radians(lat.mean()))
    y = np.radians(lat) * 6371.0
    return np.column_stack([x, y])


def morans_i(values, w: np.ndarray, permutations: int = 999, seed: int | None = 0):
    """Global Moran's I with a permutation pseudo p-value.

    Returns (I, expected I under no autocorrelation, pseudo p-value). Positive I
    means similar values cluster in space (hot tracts next to hot tracts).
    """
    z = np.asarray(values, dtype=float)
    z = z - z.mean()
    n, s0 = len(z), w.sum()
    denom = z @ z
    if denom == 0:
        raise ValueError("values are constant")
    stat = lambda v: (n / s0) * (v @ w @ v) / denom
    observed = stat(z)
    rng = np.random.default_rng(seed)
    sims = np.array([stat(rng.permutation(z)) for _ in range(permutations)])
    expected = -1.0 / (n - 1)
    # one-sided pseudo p-value in the direction of the observed value (as in PySAL)
    extreme = np.sum(sims >= observed) if observed >= expected else np.sum(sims <= observed)
    p = (extreme + 1) / (permutations + 1)
    return float(observed), expected, float(p)
