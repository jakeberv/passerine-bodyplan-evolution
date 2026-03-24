# Supplementary Code Repository

## Rates of passerine body plan evolution in time and space (*in press*)

**Authors:** Jacob S. Berv¹² [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-5962-0621), Charlotte M. Probst¹ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-4394-6931), Santiago Claramunt³⁴ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-8926-5974), J. Ryan Shipley⁵ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-9864-2498), Matt Friedman²⁶⁷ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-0114-7384), Stephen A. Smith⁸ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-2035-9531), David F. Fouhey⁹¹⁰ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-5028-5161), Brian C. Weeks¹ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-2967-2970)

**Affiliations:**  
¹ School for Environment and Sustainability, University of Michigan (Ann Arbor, MI, USA)  
² Michigan Institute for Data Science and AI in Society, University of Michigan (Ann Arbor, MI, USA)  
³ Department of Biological Sciences, University of New Orleans (New Orleans, LA, USA)  
⁴ Department of Ecology and Evolutionary Biology, University of Toronto (Toronto, ON, Canada)  
⁵ Department of Fish Ecology and Evolution, Swiss Federal Institute of Aquatic Science and Technology (Eawag) (Kastanienbaum, Switzerland)  
⁶ Museum of Paleontology, University of Michigan (Ann Arbor, MI, USA)  
⁷ Department of Earth and Environmental Sciences, University of Michigan (Ann Arbor, MI, USA)  
⁸ Department of Ecology and Evolutionary Biology, University of Michigan (Ann Arbor, MI, USA)  
⁹ Department of Computer Science, Courant Institute of Mathematical Sciences, New York University (New York, NY, USA)  
¹⁰ Department of Electrical and Computer Engineering, Tandon School of Engineering, New York University (New York, NY, USA)

**Corresponding authors:** Jacob S. Berv; Brian C. Weeks

> [!NOTE]
> **Note on reuse and analysis**
>
> This repository primarily serves as an archival record of the analyses used in this study.  
> Readers interested in running these methods on their own datasets should use the **bifrost** R package, which provides a supported and generalizable implementation of this analysis pipeline:
> [CRAN](https://cran.r-project.org/package=bifrost) · [GitHub](https://github.com/jakeberv/bifrost)

---

### At a glance

- `TemporalAnalyses.R` + `TemporalAnalyses-functions.R` → rate–shift search, simulations, and summary plots (Figures 1A–2).
- `SpatialAnalyses.R` + `SpatialAnalyses-functions.R` → species matching, range/grid processing, climate joins, spatial models, and figures (Figures 1B–D, 3–4).
- Heavy steps are cached as `.RDS` files (archived on Zenodo); scripts read those by default.

This repository provides the code to replicate the primary analysis pipeline in our study, which integrates phylogenetic comparative methods with spatial statistics. The workflow is divided into two primary R scripts and two function libraries. Each primary script loads dependencies, checks for the `data/` directory, and executes the analysis pipeline; some optional blocks reuse helper definitions across the two function files.

---

## Table of Contents

- [Quickstart](#quickstart)
- [Analytical Workflow](#analytical-workflow)
- [Data and Phylogeny](#data-and-phylogeny)
- [Citation](#citation)
- [Reproducibility Details](#reproducibility-details)
- [License](#license)
- [Contact](#contact)

---

## Quickstart

### Known caveats

- Cached inputs are used for most heavy steps, but one live GBIF rematch call remains uncommented in `SpatialAnalyses.R`.
- The scripts contain macOS/Unix-specific calls (`quartz()`, `pbmclapply`/`mclapply`), so Linux/macOS is the supported runtime target for full execution.
- `TemporalAnalyses.R` sources a pinned remote helper (`fastdivrate` `dr.R`) at runtime.

### Cached-Only Execution (Current Default)

This repository is configured for cached-mode execution. Heavy recomputation blocks are intentionally disabled, and scripts load precomputed objects from `data/`.

- `TemporalAnalyses.R` expects cached search/simulation/false-negative outputs under `data/temporal/`.
- `SpatialAnalyses.R` expects cached matching/range/grid outputs under `data/spatial/`.
- `saveRDS(...)` lines tied to heavy cached outputs are commented out.
- If recomputing from raw inputs, re-enable the relevant compute/save blocks and regenerate caches into the same `data/` paths.

### System requirements

-   **R ≥ 4.2**
-   Geospatial system libraries (needed for `sf`, `terra`, `lwgeom`) and **V8** (for `h3`/`h3jsr`):

macOS (Homebrew):

``` bash
brew install gdal geos proj pkg-config v8
```

Ubuntu/Debian:

``` bash
sudo apt-get update
sudo apt-get install -y gdal-bin libgdal-dev libgeos-dev libproj-dev libudunits2-dev libv8-dev
```

Windows:

- Install **Rtools** (matching your R version). `sf`, `terra`, `lwgeom`, and `V8` can install from CRAN binaries, but this repository's scripts use macOS/Unix-specific calls (`quartz()`, `pbmclapply`/`mclapply`) in active/optional blocks, so Linux/macOS is the supported runtime target for full execution.

> If `cairo_pdf()` fails on Linux, install Cairo headers: `sudo apt-get install libcairo2-dev`.

---

### 1) Clone

``` bash
git clone https://github.com/jakeberv/passerine-bodyplan-evolution
cd passerine-bodyplan-evolution
```

### 2) Install R packages


``` r

# CRAN packages used by active script paths
# (h3 is installed from GitHub below; parallel is base R)
cran_pkgs <- c(
  "ape","phytools","mvMORPH",
  "future","future.apply","viridis",
  "palaeoverse","boot","RColorBrewer","classInt","pbmcapply","phylolm","scales",
  "tidyverse","patchwork","readxl",
  "univariateML","evd","fitdistrplus","HDInterval",
  "rgbif","lwgeom","sf","data.table","pbapply","progressr",
  "h3jsr","geomorph","igraph","ggplot2","rnaturalearth","rnaturalearthdata",
  "tidyterra","spdep","terra","spatialreg",
  "MASS","dplyr","stringr","stargazer",
  "TeachingDemos","phangorn","ggrain","ggsignif","circlize","gridExtra",
  "corpcor","fastmatch","e1071","dispRity","caper","gplots","RRphylo",
  "matrixStats","geiger"
)

need <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)

# GitHub-only package
if (!requireNamespace("h3", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("crazycapivara/h3-r")
}

# Bioconductor package used by active figure code paths
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
}

# parallel is part of base R; no install required
message("Done. If there were no errors above, all requested packages are installed.")

```

> **External utility pinned from GitHub:**\
> `TemporalAnalyses.R` sources `fastdivrate`’s `dr.R`:
>
> ``` r
> source('https://raw.githubusercontent.com/jonchang/fastdivrate/ee8fcd0b7e623253d4d1225dadd10ffec966b8dc/R/dr.R')
> ```

---

### 3) Prepare input data

Run scripts from the repository root. Cached inputs are used for most heavy steps, and scripts expect a populated `data/` directory.

Download the published supplementary cached inputs archive from Zenodo, unpack it, and place the resulting `data/` folder at the repository root. The archive includes its own `data/README.md`, which serves as the canonical data manifest for the cached inputs. The minimum required structure is:

```         
data/
  temporal/
    01_core_inputs/external_skelevision/
      dat.mvgls.RDS
      sampled_cv.RDS
      ythlida_supertree.rescale.tre
    02_shift_search/new_bifrost/
      min10.ic20.{bic,gic}.RDS
      min10.ic40.{bic,gic}.RDS
      min20.ic20.{bic,gic}.RDS
      min20.ic40.{bic,gic}.RDS
      min30.ic20.{bic,gic}.RDS
      min30.ic40.{bic,gic}.RDS
    03_simulations/
      simresults_*_{BIC,GIC}{2,10}.RDS
    04_false_negative_outputs/
      output_FN.*.{bic,gic}.RDS
      output_FN.*.{bic,gic}.corr.RDS
    05_climate_reference/
      aba6853_tables_s8_s34.xlsx
  spatial/
    00_raw_range_source/
      ranges_4-16-22_multipolygons
    01_taxonomy_matching/
      gbif_matches.RDS
      ranges_gbif_matches.RDS
    02_ranges_and_filters/
      passeriformes.filtered*.RDS
      resident_species*_filter_res3.RDS
    03_grid_climate_weights/
      spatial_coords.RDS
      spatial_coords.12.RDS
      trMat.1000.idw.RDS
    04_temporal_bridge_inputs/
      sampled_cv_shift_metrics_8_08_25.RDS
```

---

**External datasets:**
- **WorldClim 2.1** bioclim rasters (5 arc-min) are staged under `data/spatial/05_bioclim_rasters/wc2.1_5m_bio/` for reruns of climate extraction blocks.
- **GBIF** backbone is mostly cached via `.RDS` files, but one live rematch call remains uncommented in `SpatialAnalyses.R` (`rgbif::name_backbone_checklist(...)`).
- Country boundaries via `{rnaturalearth}` are retrieved at runtime.

Published supplementary cached inputs archive on Zenodo:

- DOI: <https://doi.org/10.5281/zenodo.19198393>

### Troubleshooting

-   **Install errors for `sf`/`terra`/`lwgeom`:** ensure GDAL/GEOS/PROJ and, on Linux, `libudunits2-dev` are installed (see commands above).
-   **`h3` / `h3jsr` errors:** install system **V8** (`brew install v8` or `sudo apt-get install libv8-dev`) then reinstall the packages.
-   **Long runtimes / memory:** reduce `{future}` workers or rely on the cached `.RDS` files rather than recomputing from raw sources.
-   **`ggrain` missing:** install `ggrain` or skip the raincloud plot block that calls `create_raincloud_with_wilcox()`; that plotting section will fail without it.
-   **Strict cached/offline spatial run:** comment out the live GBIF rematch line in `SpatialAnalyses.R` and use cached matching objects.
-   **Temporal clade annotation helper:** the `extract_clade_info()` block in `TemporalAnalyses.R` uses `name_backbone_checklist()`; source `SpatialAnalyses-functions.R` in the same session or replace with `rgbif::name_backbone_checklist()`.

---

## Analytical Workflow

### 1. Temporal Analysis: Detecting Evolutionary Rate Shifts

In this stage, we identify shifts in patterns of multivariate body plan evolution across the \~50-million-year history of passerine birds. We developed a custom heuristic search algorithm to find the best-fitting multi-regime multivariate Brownian Motion (mvBM) model for high-dimensional skeletal data, with body mass included as a covariate to account for phylogenetic variation in size. The algorithm iteratively evaluates regime shifts on different clades across the phylogeny, guided by an information criterion (GIC or BIC), to identify an optimal configuration. We validated this approach through extensive simulations, which confirmed a low false-positive rate and high accuracy in recovering true shifts under various scenarios. This analysis produces the plot in **Figure 1A** and the data on shift timing and direction used in **Figure 2**.

#### Heuristic Search Algorithm

The core of the temporal analysis is the `searchOptimalConfiguration` function. This parallelized wrapper executes a greedy, stepwise search by: 1. **Generating Candidate Models:** The `generatePaintedTrees` function identifies all clades eligible for a rate shift and creates a list of candidate `simmap` trees, each representing a single potential shift. 2. **Evaluating Candidates:** Each candidate model is fit in parallel. The `fitMvglsAndExtractGIC.formula` (or `...BIC.formula`) function is used to fit a multi-regime mvBM model and calculate the corresponding information criterion. 3. **Building the Aggregated Model:** The algorithm iteratively adds the best-supported shift to the model using the `addShiftToModel` helper function, accepting it only if it improves the overall model fit beyond a predefined threshold. 4. **Validating Shifts:** After the search is complete, the `removeShiftFromTree` function is used to temporarily remove each accepted shift to calculate **IC Weights**, which quantify the statistical support for each individual shift.

### 2. Spatial Analysis: Linking Rates to Geography

In this stage, we project the evolutionary rates onto a global map and test for correlations with environmental variables.

#### Data Preparation and Metric Aggregation

We begin by harmonizing the species in the phylogeny with their geographic range polygons from BirdLife International using a modified `name_backbone_checklist` function. We then process these ranges using `merge_species_shapes_planar_parallel` and `collapse_species_data_parallel` to create a single, valid geometry for each species. We convert these polygons into a global grid using the H3 geospatial indexing system via the `h3jsr` R package. We then calculate a key species-specific metric, a **time-weighted lineage rate**, using the `compute_shift_rates.branches` function. This metric summarizes the rate history from the tree's root to each tip, giving more weight to recent evolutionary events. Finally, we use `aggregate_h3_metrics` to calculate **range-weighted means** of these rates, local species richness, and other variables for the assemblage of birds in each grid cell.

#### Investigating Phenotypic Variation

Within the spatial analysis, we analyze the structure of phenotypic variation within each grid cell to investigate various hypotheses for observed rate gradients. For the assemblage of species within each H3 grid cell, we calculated the following metrics: - **Phenotypic Disparity vs. Divergence:** We use `summarize_disparity_vs_divergence` to calculate a standardized residual by regressing pairwise Mahalanobis distances against pairwise phylogenetic divergence times. This estimates whether species are more or less different from each other than expected given their shared history (**Figure 4B**). - **Covariance Structure:** We use `local_vcv` to fit a local multivariate Brownian motion model for each cell. From the resulting trait covariance matrix, we estimate **Modularity** (the degree of trait clustering) and **Effective Dimensionality** (the number of independent trait axes) (**Figures 4A & 4C**). - **Phylogenetic Signal:** We use `aggregate_phylo_signal` to compute **Kmult** for each cell's assemblage to estimate the strength of phylogenetic signal in the multivariate data (**Figure 4D**).

#### Spatial Statistical Modeling

We use **Spatial Autoregressive Models** (via `spatialreg::lagsarlm`) to investigate the relationship between evolutionary rates and predictors like latitude, temperature seasonality (BioClim4), and local species richness. We explicitly account for spatial autocorrelation by constructing a spatial weights matrix using `spdep::dnearneigh` and `spdep::nb2listw`, based on an inverse-distance scheme. This analysis generates the data for the scatterplots in **Figures 1B, 1C, & 1D** and the global maps in **Figure 3**.

---

### Data and Phylogeny

The analyses in this study are based on three primary data sources:

1.  **Phenotypic Data:** A dataset of 12 linear skeletal measurements from 14,419 voucher specimens, representing 2,057 species of passerine birds. Missing values were imputed using a multivariate evolutionary model (Weeks et al. 2025).

2.  **Phylogeny:** A time-calibrated supertree of Passeriformes based on the phylogenetic framework of Claramunt et al. (2025), with 101 additional taxa grafted into the topology using Matrix Representation with Parsimony.

3.  **Geographic Ranges:** Global species distribution polygons from BirdLife International, used to define spatial assemblages.

-   Weeks, B. C., Zhou, Z., Probst, C. M., Berv, J. S., O’Brien, B., Benz, B. W., Skeen, H. R., Ziebell, M., Bodt, L. & Fouhey, D. F. (2025). *Skeletal trait measurements for thousands of bird species*. **Scientific Data**, 12, Article 884. <https://doi.org/10.1038/s41597-025-05234-y>
-   Claramunt, S., Sheard, C., Brown, J. W., Cortés-Ramírez, G., Cracraft, J., Su, M. M., Weeks, B. C. & Tobias, J. A. (2025). *A new time tree of birds reveals the interplay between dispersal, geographic range size, and diversification*. **Current Biology**. Advance online publication. <https://doi.org/10.1016/j.cub.2025.07.004>

---

## Citation

If you use the code in this repository, please cite:

> Berv, J. S., Probst, C. M., Claramunt, S., Shipley, J. R., Friedman, M., Smith, S. A., Fouhey, D. F., & Weeks, B. C. (2026). *Supplementary code for Rates of passerine body plan evolution in time and space* [Computer software]. GitHub. <https://github.com/jakeberv/passerine-bodyplan-evolution>

If you use the cached supplementary data archive, please also cite:

> Berv, J., Probst, C., Claramunt, S., Shipley, J. R., Friedman, M., Smith, S., Fouhey, D., & Weeks, B. (2026). *Supplementary data archive for Rates of passerine body plan evolution in time and space* (v1.0.0) [Data set]. Zenodo. <https://doi.org/10.5281/zenodo.19198393>

The manuscript citation will be added here once the in-press paper has final bibliographic details.

## Reproducibility Details

Complete session/version provenance and full references are preserved in:

- [docs/reproducibility.md](docs/reproducibility.md)

---

## License

-   **Code:** GNU General Public License, version 2 or later. See `LICENSE`.
-   **Data:** Governed by providers’ licenses (e.g., WorldClim, GBIF). Verify terms before redistribution.

---

## Contact

**Maintainer:** Jacob S. Berv — [jberv@umich.edu](mailto:jberv@umich.edu) · [jacob.berv@gmail.com](mailto:jacob.berv@gmail.com)
