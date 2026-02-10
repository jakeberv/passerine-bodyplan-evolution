# Supplemental Code Repository

## Rates of passerine body plan evolution in space and time

**Authors:** Jacob S. Berv¹² [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-5962-0621), Charlotte M. Probst¹ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-4394-6931), Santiago Claramunt³⁴ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-8926-5974), J. Ryan Shipley⁵ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-9864-2498), Matt Friedman²⁶⁷ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-0114-7384), Stephen A. Smith⁸ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-2035-9531), David F. Fouhey⁹¹⁰ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-5028-5161), Brian C. Weeks¹ [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-2967-2970)

**Affiliations:**
- ¹ School for Environment and Sustainability, University of Michigan (Ann Arbor, MI, USA)
- ² Michigan Institute for Data Science and AI in Society, University of Michigan (Ann Arbor, MI, USA)
- ³ Department of Biological Sciences, University of New Orleans (New Orleans, LA, USA)
- ⁴ Department of Ecology and Evolutionary Biology, University of Toronto (Toronto, ON, Canada)
- ⁵ Department of Fish Ecology and Evolution, Swiss Federal Institute of Aquatic Science and Technology (Eawag) (Kastanienbaum, Switzerland)
- ⁶ Museum of Paleontology, University of Michigan (Ann Arbor, MI, USA)
- ⁷ Department of Earth and Environmental Sciences, University of Michigan (Ann Arbor, MI, USA)
- ⁸ Department of Ecology and Evolutionary Biology, University of Michigan (Ann Arbor, MI, USA)
- ⁹ Department of Computer Science, Courant Institute of Mathematical Sciences, New York University (New York, NY, USA)
- ¹⁰ Department of Electrical and Computer Engineering, Tandon School of Engineering, New York University (New York, NY, USA)

**Corresponding authors:** Jacob S. Berv; Brian C. Weeks


------------------------------------------------------------------------

### At a glance

- `TemporalAnalyses.R` + `TemporalAnalyses-functions.R` → rate–shift search, simulations, and summary plots (Figures 1A–2).
- `SpatialAnalyses.R` + `SpatialAnalyses-functions.R` → species matching, range/grid processing, climate joins, spatial models, and figures (Figures 1B–D, 3–4).
- Heavy steps are cached as `.RDS` files (archived on Zenodo); scripts read those by default.

This repository provides the code to replicate the primary analysis pipeline in our study, which integrates phylogenetic comparative methods with spatial statistics. The workflow is divided into two primary R scripts, and two scripts containing new functions relevant to each. Each primary script loads dependencies, sets up a working directory, and then executes the analytical pipeline.

*To facilitate reproducibility, most of the “heavy” execution lines in the scripts are commented out, allowing the script to load preprocessed `.RDS` objects. `.RDS` files are a standard format for storing R objects, enabling the preservation of data states across sessions.*

> [!NOTE]
> **Note on reuse and analysis**
>
> This repository primarily serves as an archival record of the analyses used in this study.  
> Users interested in running these methods on their own datasets should use the **bifrost** R package, which provides a supported and generalizable implementation of the analysis pipeline:  
> [CRAN](https://cran.r-project.org/package=bifrost) · [GitHub](https://github.com/jakeberv/bifrost)

------------------------------------------------------------------------

## Quickstart

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

Windows: - Install **Rtools** (matching your R version). `sf`, `terra`, `lwgeom`, and `V8` will install from CRAN binaries.

> If `cairo_pdf()` fails on Linux, install Cairo headers: `sudo apt-get install libcairo2-dev`.

------------------------------------------------------------------------

### 1) Clone

``` bash
git clone https://github.com/jakeberv/passerine-bodyplan-evolution
cd passerine-bodyplan-evolution
```

### 2) Install R packages


``` r

# CRAN packages (everything you load except h3 + parallel)
cran_pkgs <- c(
  "ape","phytools","mvMORPH",
  "future","future.apply","viridis",
  "palaeoverse","boot","RColorBrewer","classInt","pbmcapply","phylolm","scales",
  "tidyverse","patchwork","readxl",
  "univariateML","evd","fitdistrplus","HDInterval",
  "rgbif","lwgeom","sf","data.table","pbapply","progressr",
  "h3jsr","geomorph","igraph","ggplot2","rnaturalearth","rnaturalearthdata",
  "tidyterra","spdep","terra","spatialreg",
  "MASS","dplyr","stringr","stargazer"
)

need <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)

# GitHub-only package
if (!requireNamespace("h3", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("crazycapivara/h3-r")
}

# parallel is part of base R; no install required
message("Done. If there were no errors above, all requested packages are installed.")

```

> **External utility pinned from GitHub:**\
> `TemporalAnalyses.R` sources `fastdivrate`’s `dr.R`:
>
> ``` r
> source('https://raw.githubusercontent.com/jonchang/fastdivrate/refs/heads/master/R/dr.R')
> ```

------------------------------------------------------------------------

### 3) Prepare input data

Create a `data/` folder at the repo root and place/copy the cached inputs there. These are the objects the scripts expect (relative to repo root):

```         
data/
  # multivariate trait objects used in temporal analyses
  dat.mvgls.RDS
  sampled_cv.RDS
  sampled_cv_shift_metrics_8_08_25.RDS

  # results from the shift search (consumed by both workflows)
  new_bifrost/
    min10.ic20.{bic,gic}.RDS
    min10.ic40.{bic,gic}.RDS
    min20.ic20.{bic,gic}.RDS
    min20.ic40.{bic,gic}.RDS
    min30.ic20.{bic,gic}.RDS
    min30.ic40.{bic,gic}.RDS

  # GBIF name matching caches
  gbif_matches.RDS
  ranges_gbif_matches.RDS

  # passerine range/grid objects & helpers
  passeriformes.filtered.1_res3.RDS
  passeriformes.filtered.12_res3.RDS
  passeriformes.filtered.merged.deduplicated.1.RDS
  passeriformes.filtered.merged.deduplicated.12.RDS
  resident_species1_filter_res3.RDS
  resident_species12_filter_res3.RDS
  spatial_coords.RDS
  spatial_coords.12.RDS
  trMat.1000.idw.RDS
```

------------------------------------------------------------------------

**External datasets:** - **WorldClim 2.1** bioclim rasters (5 arc‑min), especially **BioClim 4 (temperature seasonality)**. If recomputing climate summaries, download from WorldClim and set a directory path in the spatial functions; by default the script uses cached `spatial_coords*.RDS`. - GBIF backbone is accessed via `{rgbif}`; we rely on the cached `.RDS` matches by default. - Country boundaries via `{rnaturalearth}` are retrieved at runtime.

DOIs/links to archives (e.g., Zenodo) will be added when paper is published.

### Troubleshooting

-   **Install errors for `sf`/`terra`/`lwgeom`:** ensure GDAL/GEOS/PROJ and, on Linux, `libudunits2-dev` are installed (see commands above).
-   **`h3` / `h3jsr` errors:** install system **V8** (`brew install v8` or `sudo apt-get install libv8-dev`) then reinstall the packages.
-   **Long runtimes / memory:** reduce `{future}` workers or rely on the cached `.RDS` files rather than recomputing from raw sources.
-   **`ggrain` missing:** comment out `library(ggrain)` in `TemporalAnalyses-functions.R` if installation fails; the remaining code will run.

------------------------------------------------------------------------

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

------------------------------------------------------------------------

### Data and Phylogeny

The analyses in this study are based on three primary data sources:

1.  **Phenotypic Data:** A dataset of 12 linear skeletal measurements from 14,419 voucher specimens, representing 2,057 species of passerine birds. Missing values were imputed using a multivariate evolutionary model (Weeks et al. 2025).

2.  **Phylogeny:** A time-calibrated supertree of Passeriformes based on the phylogenetic framework of Claramunt et al. (2025), with 101 additional taxa grafted into the topology using Matrix Representation with Parsimony.

3.  **Geographic Ranges:** Global species distribution polygons from BirdLife International, used to define spatial assemblages.

-   Weeks, B. C., Zhou, Z., Probst, C. M., Berv, J. S., O’Brien, B., Benz, B. W., Skeen, H. R., Ziebell, M., Bodt, L. & Fouhey, D. F. (2025). *Skeletal trait measurements for thousands of bird species*. **Scientific Data**, 12, Article 884. <https://doi.org/10.1038/s41597-025-05234-y>
-   Claramunt, S., Sheard, C., Brown, J. W., Cortés-Ramírez, G., Cracraft, J., Su, M. M., Weeks, B. C. & Tobias, J. A. (2025). *A new time tree of birds reveals the interplay between dispersal, geographic range size, and diversification*. **Current Biology**. Advance online publication. <https://doi.org/10.1016/j.cub.2025.07.004>

------------------------------------------------------------------------

## Citation

If you use the code or data from this repository, please cite our paper:

> Berv, J.S., Probst, C.M., Claramunt, S., Shipley, J.R., Friedman, M., Smith, S.A., Fouhey, D.F., & Weeks, B.C. (Year). Rates of passerine body plan evolution in space and time. *Journal Placeholder*, Volume(Issue), pages. [DOI placeholder]

------------------------------------------------------------------------

# System and Package Versions (Temporal Analyses)

Analyses were conducted using the R Statistical language (version 4.4.2; R Core Team, 2024) on macOS Sequoia 15.6, using the packages boot (version 1.3.31; Angelo Canty, Ripley, 2024), maps (version 3.4.2.1; Becker OScbRA et al., 2024), future (version 1.34.0; Bengtsson H, 2021), future.apply (version 1.11.3; Bengtsson H, 2021), classInt (version 0.4.10; Bivand R, 2023), intervals (version 0.15.5; Bourgon R, 2024), rgbif (version 3.8.1.2; Chamberlain S et al., 2025), mvMORPH (version 1.2.1; Clavel J et al., 2015), fitdistrplus (version 1.2.2; Delignette-Muller ML, Dutang C, 2015), viridisLite (version 0.4.2; Garnier et al., 2023), viridis (version 0.6.5; Garnier et al., 2024), lubridate (version 1.9.3; Grolemund G, Wickham H, 2011), phylolm (version 2.6.5; Ho LST, Ane C, 2014), palaeoverse (version 1.4.0; Jones LA et al., 2023), subplex (version 1.9; King AA, Rowan T, 2024), pbmcapply (version 1.5.1; Kuang K et al., 2022), report (version 0.6.1; Makowski D et al., 2023), HDInterval (version 0.2.4; Meredith M, Kruschke J, 2022), univariateML (version 1.5.0; Moss J, 2019), tibble (version 3.2.1; Müller K, Wickham H, 2023), RColorBrewer (version 1.1.3; Neuwirth E, 2022), ape (version 5.8; Paradis E, Schliep K, 2019), patchwork (version 1.3.0.9000; Pedersen T, 2024), phytools (version 2.3.0; Revell L, 2024), corpcor (version 1.6.10; Schafer J et al., 2021), phangorn (version 2.12.1; Schliep K, 2011), TeachingDemos (version 2.13; Snow G, 2024), evd (version 2.3.7.1; Stephenson AG, 2002), survival (version 3.7.0; Therneau T, 2024), MASS (version 7.3.61; Venables WN, Ripley BD, 2002), ggplot2 (version 3.5.1; Wickham H, 2016), forcats (version 1.0.0; Wickham H, 2023), stringr (version 1.5.1; Wickham H, 2023), tidyverse (version 2.0.0; Wickham H et al., 2019), readxl (version 1.4.3; Wickham H, Bryan J, 2023), dplyr (version 1.1.4; Wickham H et al., 2023), purrr (version 1.0.2; Wickham H, Henry L, 2023), readr (version 2.1.5; Wickham H et al., 2024), scales (version 1.3.0; Wickham H et al., 2023) and tidyr (version 1.3.1; Wickham H et al., 2024).

## References

-   Angelo Canty, B. D. Ripley (2024). *boot: Bootstrap R (S-Plus) Functions*. R package version 1.3-31. A. C. Davison, D. V. Hinkley (1997). *Bootstrap Methods and Their Applications*. Cambridge University Press, Cambridge. ISBN 0-521-57391-2, <doi:10.1017/CBO9780511802843>.
-   Becker OScbRA, Minka ARWRvbRBEbTP, Deckmyn. A (2024). *maps: Draw Geographical Maps*. R package version 3.4.2.1, <https://CRAN.R-project.org/package=maps>.
-   Bengtsson H (2021). “A Unifying Framework for Parallel and Distributed Processing in R using Futures.” *The R Journal*, *13*(2), 208-227. <doi:10.32614/RJ-2021-048> <https://doi.org/10.32614/RJ-2021-048>, <https://doi.org/10.32614/RJ-2021-048>.
-   Bengtsson H (2021). “A Unifying Framework for Parallel and Distributed Processing in R using Futures.” *The R Journal*, *13*(2), 208-227. <doi:10.32614/RJ-2021-048> <https://doi.org/10.32614/RJ-2021-048>, <https://doi.org/10.32614/RJ-2021-048>.
-   Bivand R (2023). *classInt: Choose Univariate Class Intervals*. R package version 0.4-10, <https://CRAN.R-project.org/package=classInt>.
-   Bourgon R (2024). *intervals: Tools for Working with Points and Intervals*. R package version 0.15.5, <https://CRAN.R-project.org/package=intervals>.
-   Chamberlain S, Barve V, Mcglinn D, Oldoni D, Desmet P, Geffert L, Ram K (2025). *rgbif: Interface to the Global Biodiversity Information Facility API*. R package version 3.8.1.2, <https://CRAN.R-project.org/package=rgbif>. Chamberlain S, Boettiger C (2017). “R Python, and Ruby clients for GBIF species occurrence data.” *PeerJ PrePrints*. <https://doi.org/10.7287/peerj.preprints.3304v1>.
-   Clavel J, Escarguel G, Merceron G (2015). “mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data.” *Methods in Ecology and Evolution*, *6*, 1311-1319.
-   Delignette-Muller ML, Dutang C (2015). “fitdistrplus: An R Package for Fitting Distributions.” *Journal of Statistical Software*, *64*(4), 1-34. <doi:10.18637/jss.v064.i04> <https://doi.org/10.18637/jss.v064.i04>.
-   Garnier, Simon, Ross, Noam, Rudis, Robert, Camargo, Pedro A, Sciaini, Marco, Scherer, Cédric (2023). *viridis(Lite) - Colorblind-Friendly Color Maps for R*. <doi:10.5281/zenodo.4678327> <https://doi.org/10.5281/zenodo.4678327>, viridisLite package version 0.4.2, <https://sjmgarnier.github.io/viridis/>.
-   Garnier, Simon, Ross, Noam, Rudis, Robert, Camargo, Pedro A, Sciaini, Marco, Scherer, Cédric (2024). *viridis(Lite) - Colorblind-Friendly Color Maps for R*. <doi:10.5281/zenodo.4679423> <https://doi.org/10.5281/zenodo.4679423>, viridis package version 0.6.5, <https://sjmgarnier.github.io/viridis/>.
-   Grolemund G, Wickham H (2011). “Dates and Times Made Easy with lubridate.” *Journal of Statistical Software*, *40*(3), 1-25. <https://www.jstatsoft.org/v40/i03/>.
-   Ho LST, Ane C (2014). “A linear-time algorithm for Gaussian and non-Gaussian trait evolution models.” *Systematic Biology*, *63*, 397-408.
-   Jones LA, Gearty W, Allen BJ, Eichenseer K, Dean CD, Galván S, Kouvari M, Godoy PL, Nicholl CSC, Buffan L, Dillon EM, Flannery-Sutherland JT, Chiarenza AA (2023). “palaeoverse: A community-driven R package to support palaeobiological analysis.” *Methods in Ecology and Evolution*, *14(9)*, 2205-2215. <doi:10.1111/2041-210X.14099> <https://doi.org/10.1111/2041-210X.14099>.
-   King AA, Rowan T (2024). *subplex: Unconstrained Optimization using the Subplex Algorithm*. R package version 1.9, <https://CRAN.R-project.org/package=subplex>.
-   Kuang K, Kong Q, Napolitano F (2022). \_pbmcapply: Tracking the Progress of Mc\*pply with Progress Bar\_. R package version 1.5.1, <https://CRAN.R-project.org/package=pbmcapply>.
-   Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” *CRAN*. <https://easystats.github.io/report/>.
-   Meredith M, Kruschke J (2022). *HDInterval: Highest (Posterior) Density Intervals*. R package version 0.2.4, <https://CRAN.R-project.org/package=HDInterval>.
-   Moss J (2019). “univariateML: An R package for maximum likelihood estimation of univariate densities.” *Journal of Open Source Software*, *4*(44), 1863. <doi:10.21105/joss.01863> <https://doi.org/10.21105/joss.01863>, <https://doi.org/10.21105/joss.01863>.
-   Müller K, Wickham H (2023). *tibble: Simple Data Frames*. R package version 3.2.1, <https://CRAN.R-project.org/package=tibble>.
-   Neuwirth E (2022). *RColorBrewer: ColorBrewer Palettes*. R package version 1.1-3, <https://CRAN.R-project.org/package=RColorBrewer>.
-   Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.” *Bioinformatics*, *35*, 526-528. <doi:10.1093/bioinformatics/bty633> <https://doi.org/10.1093/bioinformatics/bty633>.
-   Pedersen T (2024). *patchwork: The Composer of Plots*. R package version 1.3.0.9000, commit 2695a9f0200b7fd73f295d5c8a3e13e3943078c5, <https://github.com/thomasp85/patchwork>.
-   R Core Team (2024). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
-   Revell L (2024). “phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things).” *PeerJ*, *12*, e16505. <doi:10.7717/peerj.16505> <https://doi.org/10.7717/peerj.16505>.
-   Schafer J, Opgen-Rhein R, Zuber V, Ahdesmaki M, Silva APD, Strimmer. K (2021). *corpcor: Efficient Estimation of Covariance and (Partial) Correlation*. R package version 1.6.10, <https://CRAN.R-project.org/package=corpcor>.
-   Schliep K (2011). “phangorn: phylogenetic analysis in R.” *Bioinformatics*, *27*(4), 592-593. <doi:10.1093/bioinformatics/btq706> <https://doi.org/10.1093/bioinformatics/btq706>. Schliep K, Potts A, Morrison D, Grimm G (2017). “Intertwining phylogenetic trees and networks.” *Methods in Ecology and Evolution*, *8*(10), 1212-1220.
-   Snow G (2024). *TeachingDemos: Demonstrations for Teaching and Learning*. R package version 2.13, <https://CRAN.R-project.org/package=TeachingDemos>.
-   Stephenson AG (2002). “evd: Extreme Value Distributions.” *R News*, *2*(2), 31-32. <https://CRAN.R-project.org/doc/Rnews/>.
-   Therneau T (2024). *A Package for Survival Analysis in R*. R package version 3.7-0, <https://CRAN.R-project.org/package=survival>. Terry M. Therneau, Patricia M. Grambsch (2000). *Modeling Survival Data: Extending the Cox Model*. Springer, New York. ISBN 0-387-98784-3.
-   Venables WN, Ripley BD (2002). *Modern Applied Statistics with S*, Fourth edition. Springer, New York. ISBN 0-387-95457-0, <https://www.stats.ox.ac.uk/pub/MASS4/>.
-   Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. ISBN 978-3-319-24277-4, <https://ggplot2.tidyverse.org>.
-   Wickham H (2023). *forcats: Tools for Working with Categorical Variables (Factors)*. R package version 1.0.0, <https://CRAN.R-project.org/package=forcats>.
-   Wickham H (2023). *stringr: Simple, Consistent Wrappers for Common String Operations*. R package version 1.5.1, <https://CRAN.R-project.org/package=stringr>.
-   Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” *Journal of Open Source Software*, *4*(43), 1686. <doi:10.21105/joss.01686> <https://doi.org/10.21105/joss.01686>.
-   Wickham H, Bryan J (2023). *readxl: Read Excel Files*. R package version 1.4.3, <https://CRAN.R-project.org/package=readxl>.
-   Wickham H, François R, Henry L, Müller K, Vaughan D (2023). *dplyr: A Grammar of Data Manipulation*. R package version 1.1.4, <https://CRAN.R-project.org/package=dplyr>.
-   Wickham H, Henry L (2023). *purrr: Functional Programming Tools*. R package version 1.0.2, <https://CRAN.R-project.org/package=purrr>.
-   Wickham H, Hester J, Bryan J (2024). *readr: Read Rectangular Text Data*. R package version 2.1.5, <https://CRAN.R-project.org/package=readr>.
-   Wickham H, Pedersen T, Seidel D (2023). *scales: Scale Functions for Visualization*. R package version 1.3.0, <https://CRAN.R-project.org/package=scales>.
-   Wickham H, Vaughan D, Girlich M (2024). *tidyr: Tidy Messy Data*. R package version 1.3.1, <https://CRAN.R-project.org/package=tidyr>.

# System and Package Versions (Spatial Analyses)

Analyses were conducted using the R Statistical language (version 4.4.2; R Core Team, 2024) on macOS Sequoia 15.6, using the packages geomorph (version 4.0.10; Baken E et al., 2021), data.table (version 1.17.0; Barrett T et al., 2025), Matrix (version 1.7.1; Bates D et al., 2024), maps (version 3.4.2.1; Becker OScbRA et al., 2024), future (version 1.34.0; Bengtsson H, 2021), future.apply (version 1.11.3; Bengtsson H, 2021), progressr (version 0.15.1; Bengtsson H, 2024), spatialreg (version 1.3.6; Bivand R et al., 2021), spData (version 2.3.3; Bivand R et al., 2024), spdep (version 1.3.8; Bivand R, Wong D, 2018), rgbif (version 3.8.1.2; Chamberlain S et al., 2025), mvMORPH (version 1.2.1; Clavel J et al., 2015), igraph (version 2.1.1; Csardi G, Nepusz T, 2006), viridisLite (version 0.4.2; Garnier et al., 2023), viridis (version 0.6.5; Garnier et al., 2024), tidyterra (version 0.6.1; Hernangómez D, 2023), terra (version 1.8.54; Hijmans R, 2025), stargazer (version 5.2.3; Hlavac M, 2022), subplex (version 1.9; King AA, Rowan T, 2024), pbmcapply (version 1.5.1; Kuang K et al., 2022), h3 (version 3.7.2; Kuethe S, 2022), RRPP (version 2.1.2; Collyer Adams, 2024), report (version 0.6.1; Makowski D et al., 2023), rnaturalearth (version 1.0.1; Massicotte P, South A, 2023), rgl (version 1.3.12; Murdoch D, Adler D, 2024), h3jsr (version 1.3.1; O'Brien L, 2023), ape (version 5.8; Paradis E, Schliep K, 2019), lwgeom (version 0.2.14; Pebesma E, 2024), sf (version 1.0.21; Pebesma E, Bivand R, 2023), patchwork (version 1.3.0.9000; Pedersen T, 2024), phytools (version 2.3.0; Revell L, 2024), corpcor (version 1.6.10; Schafer J et al., 2021), pbapply (version 1.7.2; Solymos P, Zawadzki Z, 2023), rnaturalearthdata (version 1.0.0; South A et al., 2024), MASS (version 7.3.61; Venables WN, Ripley BD, 2002), ggplot2 (version 3.5.1; Wickham H, 2016), stringr (version 1.5.1; Wickham H, 2023), dplyr (version 1.1.4; Wickham H et al., 2023) and scales (version 1.3.0; Wickham H et al., 2023).

## References

-   Baken E, Collyer M, Kaliontzopoulou A, Adams D (2021). “geomorph v4.0 and gmShiny: enhanced analytics and a new graphical interface for a comprehensive morphometric experience.” *Methods in Ecology and Evolution*, *12*, 2355-2363. Adams D, Collyer M, Kaliontzopoulou A, Baken E (2025). “Geomorph: Software for geometric morphometric analyses. R package version 4.0.10.” <https://cran.r-project.org/package=geomorph>. Collyer M, Adams D (2024). “RRPP: Linear Model Evaluation with Randomized Residuals in a Permutation Procedure, R package version 2.0.4.” <https://cran.r-project.org/package=RRPP>. Collyer M, Adams D (2018). “RRPP: An R package for fitting linear models to high‐dimensional data using residual randomization.”
-   Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T, Schwendinger B, Krylov I (2025). *data.table: Extension of `data.frame`*. R package version 1.17.0, <https://CRAN.R-project.org/package=data.table>.
-   Bates D, Maechler M, Jagan M (2024). *Matrix: Sparse and Dense Matrix Classes and Methods*. R package version 1.7-1, <https://CRAN.R-project.org/package=Matrix>.
-   Becker OScbRA, Minka ARWRvbRBEbTP, Deckmyn. A (2024). *maps: Draw Geographical Maps*. R package version 3.4.2.1, <https://CRAN.R-project.org/package=maps>.
-   Bengtsson H (2021). “A Unifying Framework for Parallel and Distributed Processing in R using Futures.” *The R Journal*, *13*(2), 208-227. <doi:10.32614/RJ-2021-048> <https://doi.org/10.32614/RJ-2021-048>, <https://doi.org/10.32614/RJ-2021-048>.
-   Bengtsson H (2021). “A Unifying Framework for Parallel and Distributed Processing in R using Futures.” *The R Journal*, *13*(2), 208-227. <doi:10.32614/RJ-2021-048> <https://doi.org/10.32614/RJ-2021-048>, <https://doi.org/10.32614/RJ-2021-048>.
-   Bengtsson H (2024). *progressr: An Inclusive, Unifying API for Progress Updates*. R package version 0.15.1, <https://CRAN.R-project.org/package=progressr>.
-   Bivand R, Millo G, Piras G (2021). “A Review of Software for Spatial Econometrics in R.” *Mathematics*, *9*(11). <doi:10.3390/math9111276> <https://doi.org/10.3390/math9111276>, <https://www.mdpi.com/2227-7390/9/11/1276>. Bivand R, Piras G (2015). “Comparing Implementations of Estimation Methods for Spatial Econometrics.” *Journal of Statistical Software*, *63*(18), 1-36. <doi:10.18637/jss.v063.i18> <https://doi.org/10.18637/jss.v063.i18>. Bivand R, Hauke J, Kossowski T (2013). “Computing the Jacobian in Gaussian spatial autoregressive models: An illustrated comparison of available methods.” *Geographical Analysis*, *45*(2), 150-179. <doi:10.1111/gean.12008> <https://doi.org/10.1111/gean.12008>. Bivand R, Pebesma E, Gómez-Rubio V (2013). *Applied spatial data analysis with R, Second edition*. Springer, NY. <https://asdar-book.org/>. Pebesma E, Bivand R (2023). *Spatial Data Science With Applications in R*. Chapman & Hall. <https://r-spatial.org/book/>.
-   Bivand R, Nowosad J, Lovelace R (2024). *spData: Datasets for Spatial Analysis*. R package version 2.3.3, <https://CRAN.R-project.org/package=spData>.
-   Bivand R, Wong D (2018). “Comparing implementations of global and local indicators of spatial association.” *TEST*, *27*(3), 716-748. <doi:10.1007/s11749-018-0599-x> <https://doi.org/10.1007/s11749-018-0599-x>. Roger Bivand (2022). “R Packages for Analyzing Spatial Data: A Comparative Case Study with Areal Data.” *Geographical Analysis*, *54*(3), 488-518. <doi:10.1111/gean.12319> <https://doi.org/10.1111/gean.12319>. Bivand R, Pebesma E, Gómez-Rubio V (2013). *Applied spatial data analysis with R, Second edition*. Springer, NY. <https://asdar-book.org/>. Pebesma E, Bivand R (2023). *Spatial Data Science With Applications in R*. Chapman & Hall. <https://r-spatial.org/book/>.
-   Chamberlain S, Barve V, Mcglinn D, Oldoni D, Desmet P, Geffert L, Ram K (2025). *rgbif: Interface to the Global Biodiversity Information Facility API*. R package version 3.8.1.2, <https://CRAN.R-project.org/package=rgbif>. Chamberlain S, Boettiger C (2017). “R Python, and Ruby clients for GBIF species occurrence data.” *PeerJ PrePrints*. <https://doi.org/10.7287/peerj.preprints.3304v1>.
-   Clavel J, Escarguel G, Merceron G (2015). “mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data.” *Methods in Ecology and Evolution*, *6*, 1311-1319.
-   Csardi G, Nepusz T (2006). “The igraph software package for complex network research.” *InterJournal*, *Complex Systems*, 1695. <https://igraph.org>. Csárdi G, Nepusz T, Traag V, Horvát S, Zanini F, Noom D, Müller K (2025). *igraph: Network Analysis and Visualization in R*. <doi:10.5281/zenodo.7682609> <https://doi.org/10.5281/zenodo.7682609>, R package version 2.1.1, <https://CRAN.R-project.org/package=igraph>.
-   Garnier, Simon, Ross, Noam, Rudis, Robert, Camargo, Pedro A, Sciaini, Marco, Scherer, Cédric (2023). *viridis(Lite) - Colorblind-Friendly Color Maps for R*. <doi:10.5281/zenodo.4678327> <https://doi.org/10.5281/zenodo.4678327>, viridisLite package version 0.4.2, <https://sjmgarnier.github.io/viridis/>.
-   Garnier, Simon, Ross, Noam, Rudis, Robert, Camargo, Pedro A, Sciaini, Marco, Scherer, Cédric (2024). *viridis(Lite) - Colorblind-Friendly Color Maps for R*. <doi:10.5281/zenodo.4679423> <https://doi.org/10.5281/zenodo.4679423>, viridis package version 0.6.5, <https://sjmgarnier.github.io/viridis/>.
-   Hernangómez D (2023). “Using the tidyverse with terra objects: the tidyterra package.” *Journal of Open Source Software*, *8*(91), 5751. ISSN 2475-9066, <doi:10.21105/joss.05751> <https://doi.org/10.21105/joss.05751>, <https://doi.org/10.21105/joss.05751>.
-   Hijmans R (2025). *terra: Spatial Data Analysis*. R package version 1.8-54, <https://CRAN.R-project.org/package=terra>.
-   Hlavac M (2022). *stargazer: Well-Formatted Regression and Summary Statistics Tables*. Social Policy Institute, Bratislava, Slovakia. R package version 5.2.3, <https://CRAN.R-project.org/package=stargazer>.
-   King AA, Rowan T (2024). *subplex: Unconstrained Optimization using the Subplex Algorithm*. R package version 1.9, <https://CRAN.R-project.org/package=subplex>.
-   Kuang K, Kong Q, Napolitano F (2022). \_pbmcapply: Tracking the Progress of Mc\*pply with Progress Bar\_. R package version 1.5.1, <https://CRAN.R-project.org/package=pbmcapply>.
-   Kuethe S (2022). *h3: R Bindings for H3*. R package version 3.7.2, commit 6b658e832f6907581d9d5c5296d611c4e4cf372a, <https://github.com/crazycapivara/h3-r>.
-   M. L. Collyer D. C. Adams (2024). *RRPP: Linear Model Evaluation with Randomized Residuals in a Permutation Procedure. R package version 2.1.2.*. <https://CRAN.R-project.org/package=RRPP>. M. L. Collyer D. C. Adams (2018). “RRPP: An R package for fitting linear models to high‐dimensional data using residual randomization.” *Methods in Ecology and Evolution*, *9*, 1772-1779. <https://doi.org/10.1111/2041-210X.13029>.
-   Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” *CRAN*. <https://easystats.github.io/report/>.
-   Massicotte P, South A (2023). *rnaturalearth: World Map Data from Natural Earth*. R package version 1.0.1, <https://CRAN.R-project.org/package=rnaturalearth>.
-   Murdoch D, Adler D (2024). *rgl: 3D Visualization Using OpenGL*. R package version 1.3.12, <https://CRAN.R-project.org/package=rgl>.
-   O'Brien L (2023). *h3jsr: Access Uber's H3 Library*. R package version 1.3.1, <https://CRAN.R-project.org/package=h3jsr>.
-   Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.” *Bioinformatics*, *35*, 526-528. <doi:10.1093/bioinformatics/bty633> <https://doi.org/10.1093/bioinformatics/bty633>.
-   Pebesma E (2024). *lwgeom: Bindings to Selected 'liblwgeom' Functions for Simple Features*. R package version 0.2-14, <https://CRAN.R-project.org/package=lwgeom>.
-   Pebesma E, Bivand R (2023). *Spatial Data Science: With applications in R*. Chapman and Hall/CRC. <doi:10.1201/9780429459016> <https://doi.org/10.1201/9780429459016>, <https://r-spatial.org/book/>. Pebesma E (2018). “Simple Features for R: Standardized Support for Spatial Vector Data.” *The R Journal*, *10*(1), 439-446. <doi:10.32614/RJ-2018-009> <https://doi.org/10.32614/RJ-2018-009>, <https://doi.org/10.32614/RJ-2018-009>.
-   Pedersen T (2024). *patchwork: The Composer of Plots*. R package version 1.3.0.9000, commit 2695a9f0200b7fd73f295d5c8a3e13e3943078c5, <https://github.com/thomasp85/patchwork>.
-   R Core Team (2024). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
-   Revell L (2024). “phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things).” *PeerJ*, *12*, e16505. <doi:10.7717/peerj.16505> <https://doi.org/10.7717/peerj.16505>.
-   Schafer J, Opgen-Rhein R, Zuber V, Ahdesmaki M, Silva APD, Strimmer. K (2021). *corpcor: Efficient Estimation of Covariance and (Partial) Correlation*. R package version 1.6.10, <https://CRAN.R-project.org/package=corpcor>.
-   Solymos P, Zawadzki Z (2023). \_pbapply: Adding Progress Bar to '\*apply' Functions\_. R package version 1.7-2, <https://CRAN.R-project.org/package=pbapply>.
-   South A, Michael S, Massicotte P (2024). *rnaturalearthdata: World Vector Map Data from Natural Earth Used in 'rnaturalearth'*. R package version 1.0.0, <https://CRAN.R-project.org/package=rnaturalearthdata>.
-   Venables WN, Ripley BD (2002). *Modern Applied Statistics with S*, Fourth edition. Springer, New York. ISBN 0-387-95457-0, <https://www.stats.ox.ac.uk/pub/MASS4/>.
-   Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. ISBN 978-3-319-24277-4, <https://ggplot2.tidyverse.org>.
-   Wickham H (2023). *stringr: Simple, Consistent Wrappers for Common String Operations*. R package version 1.5.1, <https://CRAN.R-project.org/package=stringr>.
-   Wickham H, François R, Henry L, Müller K, Vaughan D (2023). *dplyr: A Grammar of Data Manipulation*. R package version 1.1.4, <https://CRAN.R-project.org/package=dplyr>.
-   Wickham H, Pedersen T, Seidel D (2023). *scales: Scale Functions for Visualization*. R package version 1.3.0, <https://CRAN.R-project.org/package=scales>. \>

## License

-   **Code:** (add a `LICENSE` file).
-   **Data:** Governed by providers’ licenses (e.g., WorldClim, GBIF). Verify terms before redistribution.

------------------------------------------------------------------------

## Contact

**Maintainer:** Jacob S. Berv — [jberv\@umich.edu](mailto:jberv@umich.edu) \ [jacob.berv\@gmail.com](mailto:jacob.berv@gmail.com)
