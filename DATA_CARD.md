# Data card: 2013 blueberry sensory subset

## Dataset summary

| Item | Value |
|---|---|
| Unit of observation | One genotype sampled in one location-by-harvest environment |
| Rows | 50 |
| Year | 2013 |
| Locations | C, H, W (codes retained from the publication) |
| Environment blocks | 9 location-by-harvest combinations |
| Genotypes | Emerald, Endura, Farthing, Meadowlark, Primadonna, Scintilla |
| Sensory fields | overall liking, texture, sweetness, sourness, flavor |
| Personal data | None; no panelist-level records are present |

The five outcomes are published sample-level sensory summaries. Overall liking and texture are hedonic measures; sweetness, sourness, and flavor are perceived intensities. Consult the source article for the general Labeled Magnitude Scale design and panel procedures.

## Source and attribution

Gilbert JL, Guthart MJ, Gezan SA, Pisaroglo de Carvalho M, Schwieterman ML, Colquhoun TA, et al. (2015). “Identifying Breeding Priorities for Blueberry Flavor Using Biochemical, Sensory, and Genotype by Environment Analyses.” *PLOS ONE* 10(9): e0138494.

- Article DOI: <https://doi.org/10.1371/journal.pone.0138494>
- S2 supplement DOI: <https://doi.org/10.1371/journal.pone.0138494.s002>
- Article data-availability statement: all relevant data are in the paper and supporting information.
- License: Creative Commons Attribution 4.0 International (CC BY 4.0).
- Verified S2 XLSX SHA-256: `fade6546de831bf15c5b4ef0fa234b13d25764a599ae87a04f906bbf051acab6`.

The 94,712-byte source workbook is downloaded to `data/raw/` and deliberately ignored by Git. The small transformed CSV is committed under the same CC BY terms so the verified analysis can run offline. Attribution must be preserved in redistributions.

## Transformation

`scripts/prepare_data.R` performs the complete transformation:

1. download S2 and verify its pinned SHA-256;
2. read worksheet `S2 Table`, with three descriptive rows above the header;
3. retain rows whose `gID` is non-missing and `Year == 2013`;
4. select `gID`, year, location, harvest, genotype, Overall Liking, Texture, Sweetness, Sourness, and Flavor;
5. require complete values for all retained fields;
6. round the five published sensory values to their reported one-decimal precision;
7. sort deterministically by location, harvest, and genotype.

The script asserts 50 unique sample IDs, nine environments, and six named genotypes. `data/source_manifest.json` records the executable provenance.

## Repeated-measures structure and imbalance

These are not 50 fully exchangeable biological replicates. Genotypes recur across location-by-harvest environments. The main PCA bootstrap therefore resamples whole environment blocks. Cluster stability removes whole environments, and external predictive validation leaves out one whole environment at a time.

The design is mildly incomplete:

- Emerald, Endura, Farthing, and Meadowlark each have 9 observations;
- Primadonna has 8 observations because H3 is absent;
- Scintilla has 6 observations because all W environments are absent.

This imbalance limits separation of genotype and environment and is explicitly retained rather than imputed.

## Appropriate uses

- teaching or benchmarking PCA and clustering diagnostics;
- studying scaling sensitivity and resampling stability;
- exploratory comparison of sample-level sensory profiles;
- demonstrating group-aware external validation.

## Inappropriate uses

- inferring individual consumer preferences or demographic effects;
- making health, food-safety, or commercial cultivar recommendations;
- claiming genotype causes a sensory outcome from this observational subset;
- treating an exploratory cluster as a validated biological class;
- reconstructing panelist-level data that the supplement does not provide.

## Known limitations

- only one year is analyzed;
- only six repeatedly sampled genotypes are in this subset;
- only nine environment blocks are available for resampling;
- sensory values are already aggregated, so consumer-level uncertainty cannot be rebuilt;
- location codes are preserved but farm management and weather are not modeled here;
- no missing values are imputed, and the incomplete design can influence validation.

