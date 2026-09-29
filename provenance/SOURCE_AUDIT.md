# Public Source and Analysis Audit

## Canonical public source

The project uses the article's S2 supporting-information workbook:

- article: <https://doi.org/10.1371/journal.pone.0138494>
- supplement: <https://doi.org/10.1371/journal.pone.0138494.s002>
- pinned SHA-256: `fade6546de831bf15c5b4ef0fa234b13d25764a599ae87a04f906bbf051acab6`
- source size: 94,712 bytes
- retrieval date used for this reconstruction: 2026-09-29
- license: CC BY 4.0, with attribution to Gilbert et al.

PLOS states that all relevant data are available in the article and supporting information. The analysis uses only the 2013 five-trait sensory subset required for the stated research questions.

## Analytical risks addressed

1. **Scaling sensitivity.** All five fields use related sensory scales, but their empirical variances differ. The project compares standardized and unscaled PCA and distance calculations explicitly.
2. **Cluster-count leakage.** The number of clusters is selected without using the six known genotype labels; genotype is retained only for external checks.
3. **Repeated environments.** Genotypes recur across locations and harvests, so uncertainty and validation respect environment blocks rather than treating every row as independent.
4. **Single-run instability.** Gap statistics, silhouettes, method agreement, environment subsampling, and scale sensitivity are all reported.
5. **In-sample label agreement.** The workflow adds leave-one-environment-out validation and environment-constrained permutations.

## Publication boundary

Safe to publish with attribution:

- executable downloader and transformation code;
- the 50-row derived public-data subset;
- public-source and artifact hashes;
- aggregate results, original figures, tests, and report prose.

Excluded:

- panelist-level records, which are not present in the analyzed supplement;
- unrelated biochemical and volatile columns not required for this project;
- signed download URLs, account metadata, personal information, and credentials.
