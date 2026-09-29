# ブルーベリー官能評価の多変量解析 / Multivariate Analysis of Blueberry Sensory Data

## Abstract

This study analyzes the open S2 supporting data from Gilbert et al. (2015): 50 location-by-harvest-by-genotype samples from 2013, 9 environment blocks, 6 genotypes, and five sensory ratings. The first two standardized principal components explain 89.2% of total variance, compared with 91.0% without scaling. Environment-block resampling, scale sensitivity, gap statistics, silhouettes, cluster subsampling, and leave-one-environment-out validation distinguish reproducible continuous structure from weak evidence for discrete groups.

The central result is conservative: the gap statistic's first-SE rule selects `k=1`, so the six known genotypes must not be equated with six natural sensory clusters. When at least two clusters are forced, mean silhouette selects an exploratory `k=4` solution, but this is only a descriptive partition. Separately, genotype labels that never enter feature construction achieve above-permutation accuracy in leave-one-environment-out nearest-centroid validation. The five ratings therefore contain genotype-related signal that generalizes across environments, but that signal does not naturally form six compact groups.

## 1. Research questions

1. What are the main continuous axes of the five sensory ratings in 2013?
2. Does scaling materially change the PCA or clustering conclusions?
3. Are PCA loadings stable under environment-block resampling?
4. How many discrete groups does the data support, and how stable are they across methods, environment subsets, and scaling choices?
5. Without using genotype labels to select `k`, can sensory features identify genotype after holding out complete environments?

This is an exploratory measurement study. It estimates no treatment effect and makes no causal claim.

## 2. Data source and license

The public data are from Jessica L. Gilbert et al. (2015), *Identifying Breeding Priorities for Blueberry Flavor Using Biochemical, Sensory, and Genotype by Environment Analyses*, PLOS ONE 10(9): e0138494, DOI: [10.1371/journal.pone.0138494](https://doi.org/10.1371/journal.pone.0138494). The article and supporting information are published under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). `scripts/download_source.R` pins the S2 XLSX URL and SHA-256; `scripts/prepare_data.R` selects 2013 and the five sensory fields from `S2 Table`.

The workflow treats scaling as an explicit sensitivity question, selects the cluster count without genotype labels, and quantifies loading and clustering stability before genotype is used for external checks.

## 3. Design and methods

### 3.1 Repeated structure

The 50 observations come from 9 environment blocks formed by location C/H/W and harvest. Genotypes recur across environments, so an ordinary row bootstrap would overstate independent information. PCA uncertainty uses 1000 environment-block bootstrap replicates: each replicate samples nine complete blocks with replacement, refits the mean, standard deviation, and PCA, and aligns the first two axes by permutation and sign.

### 3.2 PCA and scaling

The primary analysis z-standardizes each of the five traits. A sensitivity analysis retains the original scales: the fields use related rating scales, but their empirical variances differ. The report includes eigenvector loadings, variable-component correlation loadings, explained variance, score correlations, and angles between the first two subspaces.

### 3.3 Cluster count, stability, and external labels

The gap statistic uses 200 reference bootstraps for `k=1…8`; mean silhouettes are computed for k-means and Ward.D2 over `k=2…8`. The primary rule is the gap first-SE criterion. If it returns `k=1`, a silhouette-selected solution with `k≥2` is reported for visualization only. Stability uses 500 environment subsamples, each retaining 8/9 complete blocks without replacement (`ceiling(0.8 × 9)`). ARI also measures agreement between k-means and Ward and between standardized and unscaled solutions.

Genotype labels never enter PCA, cluster-count selection, or clustering. External checks are: (1) cluster-genotype ARI with 4999 genotype-label permutations constrained within environment blocks; and (2) leave-one-environment-out nearest-centroid classification, with scaling and genotype centroids estimated only from the remaining environments and 999 similarly constrained permutations.

## 4. Results

### 4.1 PCA

Standardized PC1 and PC2 explain 56.5% and 32.7%, or 89.2% cumulatively; the first two unscaled axes explain 46.0% and 45.0%.

| Sensory variable | PC1 correlation loading | PC2 correlation loading |
|---|---:|---:|
| overall liking | 0.921 | -0.270 |
| texture | 0.744 | 0.489 |
| sweetness | 0.828 | -0.502 |
| sourness | -0.004 | 0.973 |
| flavor | 0.859 | 0.354 |

PC1 aligns overall liking, texture, sweetness, and flavor, while sourness is near zero; PC2 is dominated by positive sourness variation. In the environment-block bootstrap, median cosines between the resampled and full-sample PC1/PC2 axes are 0.997 and 0.993. The median maximum angle between the first two subspaces is 4.3° (95% interval 1.3°–11.4°). This quantifies axis uncertainty more directly than a single biplot.

The first two standardized and unscaled subspace angles are 11.6° and 13.5°; corresponding score correlations are 0.822 and 0.822. The raw scale gives sourness substantially more influence on the first axis, so scaling cannot be treated as inconsequential.

### 4.2 Evidence for discrete clusters

For standardized data, the first-SE rule selects `k=1`, so the primary analysis does not support a discrete partition. The corresponding rule selects `k=1` without scaling. The primary evidence therefore favors a continuous sensory space rather than clearly separated natural clusters.

After imposing `k≥2`, mean silhouette selects `k=4` for standardized data. K-means and Ward agree with ARI 0.735. Across environment subsamples, median ARI relative to the full-sample solution is 0.952 for k-means (95% interval 0.721–1.000) and 0.864 for Ward (0.622–1.000). Standardized/unscaled agreement is 0.793 for k-means and 0.769 for Ward.

The exploratory `k=4` clusters agree weakly with the six genotype labels: k-means ARI=0.069 (within-environment permutation p=0.0042) and Ward ARI=0.057 (p=0.0106). The clusters therefore cannot be named as genotype groups, and the known count of six genotypes cannot be used to justify six clusters.

### 4.3 Genotype signal across environments

Leave-one-environment-out nearest-centroid accuracy is 46.0% and balanced accuracy is 47.7%. Constrained permutations have mean accuracy 16.8%, a 95% interval of 6.0%–28.1%, and Monte Carlo p=0.0010. The sensory vector therefore contains genotype information that generalizes across environments, but this may be an overlapping continuous displacement and does not require six separated clusters.

## 5. Conclusion

The 2013 blueberry sensory data have a clear, reproducible, low-dimensional continuous structure: one combined liking/sweetness/flavor axis and one sourness-dominated axis. Genotypes show sensory differences that predict across environments, but unsupervised evidence does not support a natural partition into six genotype clusters. The defensible conclusion is continuous genotype-related sensory signal; discrete clusters are a weaker exploratory summary that is sensitive to method, environment, and scaling.

## 6. Limitations

- Only nine environment blocks are available, so the tails of block-bootstrap intervals are themselves imprecise.
- Genotypes are not fully balanced in 2013: Primadonna lacks H3 and Scintilla lacks three W environments.
- The table contains sample-level sensory summaries, not individual ratings from all 217 consumers, so participant-level random effects and demographic heterogeneity cannot be reconstructed.
- Nearest-centroid validation establishes signal in the five ratings, not genetic causality; location, harvest time, and their interactions may still affect ratings.
- Gap statistics, silhouettes, and ARI answer different questions; exploratory `k=4` is not an estimate of a true cluster count.
- Results cover 2013 and six genotypes repeated in this subset, not all 19 genotypes or other years.

## 7. Reproducibility and artifacts

Rebuild completely from the public source:

```bash
make reproduce
```

Re-run offline from the committed, attributed 50-row derived snapshot:

```bash
make all
```

The random seed is `20260929`. `results/summary.json` provides machine-readable conclusions; every table, PNG, session record, and SHA-256 manifest is generated by the script. See `DATA_CARD.md` and `provenance/SOURCE_AUDIT.md` for source and publication boundaries.
