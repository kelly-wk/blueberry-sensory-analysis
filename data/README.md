# Data directory

- `blueberry_2013_sensory.csv` is the 50-row analysis-ready subset derived from Gilbert et al. (2015) S2 under CC BY 4.0.
- `source_manifest.json` records source URL, DOI, hash, license, and transformation rules.
- `raw/` is gitignored. `make reproduce` downloads the 94,712-byte S2 workbook there, verifies its SHA-256, and rebuilds the CSV.

See `../DATA_CARD.md` for full attribution and limitations.

