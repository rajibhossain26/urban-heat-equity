"""Download tract-level social vulnerability indicators from the Census API.

    python scripts/get_acs.py --state 13 --county 059 --year 2023 --out data/acs_GA_059.csv

Uses the ACS 5-year detailed tables. No API key is needed for one county;
add --key YOUR_KEY for heavier use (https://api.census.gov/data/key_signup.html).
"""

import argparse
import csv
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

AGE_65_PLUS = [f"B01001_{i:03d}E" for i in list(range(20, 26)) + list(range(44, 50))]
VARIABLES = {
    "B01003_001E": "population",
    "B19013_001E": "median_income",
    "B17001_001E": "_poverty_universe",
    "B17001_002E": "_below_poverty",
    "B03002_001E": "_race_total",
    "B03002_003E": "_nh_white",
    "B01001_001E": "_age_total",
    **{v: "_age65" for v in AGE_65_PLUS},
}


def fetch(state: str, county: str, year: int, key: str | None = None):
    params = {"get": ",".join(VARIABLES), "for": "tract:*", "in": f"state:{state} county:{county}"}
    if key:
        params["key"] = key
    url = f"https://api.census.gov/data/{year}/acs/acs5?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as resp:
        rows = json.load(resp)
    header, data = rows[0], rows[1:]
    return [dict(zip(header, r)) for r in data]


def to_number(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v < 0 else v  # the API uses large negative codes for missing values


def tidy(records):
    out = []
    for r in records:
        num = {k: to_number(r[k]) for k in VARIABLES}
        age65 = sum(num[v] or 0 for v in AGE_65_PLUS)

        def pct(a, b):
            return None if not b or a is None else round(100 * a / b, 2)

        out.append(
            {
                "GEOID": r["state"] + r["county"] + r["tract"],
                "population": num["B01003_001E"],
                "median_income": num["B19013_001E"],
                "pct_poverty": pct(num["B17001_002E"], num["B17001_001E"]),
                "pct_minority": pct(
                    None if num["B03002_003E"] is None or num["B03002_001E"] is None
                    else num["B03002_001E"] - num["B03002_003E"],
                    num["B03002_001E"],
                ),
                "pct_age65": pct(age65, num["B01001_001E"]),
            }
        )
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--state", required=True, help="2-digit state FIPS, e.g. 13")
    p.add_argument("--county", required=True, help="3-digit county FIPS, e.g. 059")
    p.add_argument("--year", type=int, default=2023, help="last year of the ACS 5-year release")
    p.add_argument("--key", default=None)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)

    rows = tidy(fetch(a.state, a.county, a.year, a.key))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} tracts to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
