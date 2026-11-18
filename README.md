# Urban Heat & Greenness Equity

**Who lives in the hottest, least green neighborhoods?** This is an open, reproducible pipeline. It measures summer land surface temperature (LST), greenness and tree cover for every census tract with Google Earth Engine. It then tests whether socially vulnerable communities carry more of the heat burden, using regression and spatial statistics.

[![tests](https://github.com/rajibhossain26/urban-heat-equity/actions/workflows/tests.yml/badge.svg)](https://github.com/rajibhossain26/urban-heat-equity/actions/workflows/tests.yml)
![Google Earth Engine](https://img.shields.io/badge/Google%20Earth%20Engine-4285F4?logo=google-earth&logoColor=white)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)

Extreme heat is the deadliest weather hazard in the US. Exposure is not shared equally: neighborhoods with less tree canopy and more pavement can be several degrees hotter. This pipeline turns that question into numbers for any US county.

## Pipeline

```
Landsat 8/9 (summer, cloud-masked) ─┐
ESA WorldCover (trees, built, water) ├─► GEE zonal stats by tract ─► heat_greenness.csv ─┐
US Census TIGER 2020 tracts ────────┘                                                  ├─► run_analysis.py ─► report + figures
ACS 5-year summary files (poverty, income, race, age 65+) ─► acs.csv ──────────┘
```

| Step | File | What it does |
| --- | --- | --- |
| 1 | [`gee/heat_greenness_by_tract.js`](gee/heat_greenness_by_tract.js) | Builds a median summer composite of Landsat Collection 2 surface temperature (°C) and NDVI, masking clouds, shadows and open water. It summarizes these by tract, along with tree and built-up fractions, and exports a CSV. |
| 2 | [`scripts/get_acs.py`](scripts/get_acs.py) | Downloads tract-level ACS indicators from the Census summary files, no API key needed (standard library only; pass `--key` to use the faster Census API). |
| 3 | [`scripts/run_analysis.py`](scripts/run_analysis.py) | Joins the two, builds a vulnerability index (mean z-score), and compares heat across vulnerability quintiles. It fits `LST ~ tree cover + vulnerability` with robust standard errors and tests spatial clustering with Moran's I. |

## Run it

```bash
git clone https://github.com/rajibhossain26/urban-heat-equity.git
cd urban-heat-equity
pip install -e ".[dev]"

# 1. Paste gee/heat_greenness_by_tract.js into the Earth Engine Code Editor, set the state and county FIPS,
#    run it, export the CSV from the Tasks tab, and save it to data/
# 2. Census data for the same county
python scripts/get_acs.py --state 13 --county 059 --year 2023 --out data/acs_GA_059.csv
# 3. Analysis
python scripts/run_analysis.py --heat data/heat_greenness_tracts_GA_059.csv --acs data/acs_GA_059.csv --out results/
```

`results/results.md` gets a quintile table, regression output and Moran's I. The figures are `lst_by_vulnerability.png` and `lst_vs_tree_cover.png`.

## Results: Athens–Clarke County, GA

<!-- TODO (Rajib): run the pipeline and paste the two figures and 2-3 key numbers here, e.g.
     "The most vulnerable fifth of tracts is X °C hotter than the least vulnerable fifth;
      each +10 percentage points of tree cover is associated with -X °C." -->

*Coming soon.*

## Methods notes

- **LST vs. air temperature.** Landsat measures the temperature of the ground surface in mid-morning, when the satellite passes over. It is a well-established proxy for heat exposure, but it is not the air temperature people feel. Results should be read as relative differences between tracts.
- **Median composite** across several summers reduces the effect of single hot days and leftover cloud edges.
- **Water is masked** before averaging, so tracts with lakes or rivers are not artificially cool.
- **Vulnerability index** is the mean z-score of % below poverty, % minority (all except non-Hispanic white), % age 65+, and median household income (reversed). You can edit `INDICATORS` in `run_analysis.py`.
- **Spatial autocorrelation.** If the OLS residuals are clustered (significant Moran's I), the standard errors are too optimistic and a spatial lag or error model is the next step.
- **Tests:** `pytest` runs the statistics on synthetic data with known answers, plus an end-to-end run of the analysis script.

## Possible extensions

- [ ] Heat-related ED visits or CDC PLACES health outcomes as the dependent variable
- [ ] Spatial error and lag models (PySAL `spreg`)
- [ ] Multiple counties or the whole Atlanta metro, with a county fixed effect
- [ ] Street-level greenery from street-view imagery

## Author

**Rajib Hossain**, PhD researcher in GIScience & environmental health, University of Georgia · [LinkedIn](https://www.linkedin.com/in/hossainrajib)

## Data sources

- Landsat 8–9 Collection 2 Level-2, USGS (via Google Earth Engine `LANDSAT/LC08/C02/T1_L2`, `LANDSAT/LC09/C02/T1_L2`)
- ESA WorldCover 10 m 2021 v200 (`ESA/WorldCover/v200`)
- US Census Bureau TIGER/Line 2020 tracts (`TIGER/2020/TRACT`) and American Community Survey 5-year estimates

MIT License
