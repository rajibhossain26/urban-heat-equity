"""Download tract-level social vulnerability indicators from the ACS 5-year data.

    python scripts/get_acs.py --state 13 --county 059 --year 2023 --out data/acs_GA_059.csv

By default this reads the Census Bureau's table-based summary files
(https://www2.census.gov/programs-surveys/acs/summary_file/), which need no
API key. They are national files (about 400 MB in total for the five tables
used here), so the download streams them and keeps only the county's tracts;
expect a few minutes.

The Census Data API now requires a key for every request. If you have one
(https://api.census.gov/data/key_signup.html), pass --key YOUR_KEY to use the
API instead, which is much faster.
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

SUMMARY_FILE_URL = ("https://www2.census.gov/programs-surveys/acs/summary_file/{year}/"
                    "table-based-SF/data/5YRData/acsdt5y{year}-{table}.dat")


def fetch_api(state: str, county: str, year: int, key: str):
    params = {"get": ",".join(VARIABLES), "for": "tract:*", "in": f"state:{state} county:{county}", "key": key}
    url = f"https://api.census.gov/data/{year}/acs/acs5?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as resp:
        rows = json.load(resp)
    header, data = rows[0], rows[1:]
    return [dict(zip(header, r)) for r in data]


def summary_column(api_var: str) -> str:
    """API name -> summary-file column name, e.g. B01001_020E -> B01001_E020."""
    table, rest = api_var.split("_")
    return f"{table}_E{rest[:3]}"


def parse_summary_lines(lines, state: str, county: str, variables, records=None):
    """Read one summary-file table (an iterable of text lines, header first)
    and add the county's tract values to ``records`` (dict keyed by GEOID).
    Rows are sorted by GEO_ID, so reading stops after the county's block."""
    records = {} if records is None else records
    prefix = f"1400000US{state}{county}"
    it = iter(lines)
    header = next(it).rstrip("\n").split("|")
    cols = {v: header.index(summary_column(v)) for v in variables if summary_column(v) in header}
    seen = False
    for line in it:
        if not line.startswith(prefix):
            if seen:
                break
            continue
        seen = True
        fields = line.rstrip("\n").split("|")
        fips = fields[0][len("1400000US"):]
        rec = records.setdefault(fips, {"state": fips[:2], "county": fips[2:5], "tract": fips[5:]})
        for v, i in cols.items():
            rec[v] = fields[i]
    return records


def fetch_summary_files(state: str, county: str, year: int):
    records: dict = {}
    for table in sorted({v.split("_")[0] for v in VARIABLES}):
        wanted = [v for v in VARIABLES if v.startswith(table + "_")]
        url = SUMMARY_FILE_URL.format(year=year, table=table.lower())
        print(f"Reading {url}", file=sys.stderr)
        req = urllib.request.Request(url, headers={"User-Agent": "urban-heat-equity"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            lines = (raw.decode("utf-8") for raw in resp)
            parse_summary_lines(lines, state, county, wanted, records)
    return [records[k] for k in sorted(records)]


def fetch(state: str, county: str, year: int, key: str | None = None):
    return fetch_api(state, county, year, key) if key else fetch_summary_files(state, county, year)


def to_number(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v < 0 else v  # Census uses large negative codes for missing values


def tidy(records):
    out = []
    for r in records:
        num = {k: to_number(r.get(k)) for k in VARIABLES}
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
    p.add_argument("--key", default=None, help="Census API key (optional; faster than the summary files)")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)

    rows = tidy(fetch(a.state, a.county, a.year, a.key))
    if not rows:
        raise SystemExit(f"No tracts found for state {a.state}, county {a.county}.")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} tracts to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
