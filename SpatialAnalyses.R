# Master R code for executing analyses for
# Rates of Avian Body Plan Evolution in Space and Time
# Spatial Analyses (primarily)

#load necessary packages
library(rgbif)
library(lwgeom)
library(sf)
library(data.table)
library(pbapply)
library(future.apply)
library(pbmcapply)
library(progressr)
library(phytools)
library(h3jsr)
library(geomorph)
library(mvMORPH)
library(future)
library(igraph)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(parallel)
library(viridis)
library(tidyterra)
library(spdep)
library(terra)  # For handling raster files
library(h3)
library(spatialreg)
library(MASS)
library(dplyr)
library(stringr)
library(stargazer)
library(patchwork)

setwd('/Users/cotinga/jacob.berv@gmail.com/Code/passerine-bodyplan-evolution')

#read in function definitions
source("SpatialAnalyses-functions.R")
source("TemporalAnalyses-functions.R")
#length(lsf.str(envir = .GlobalEnv))

#system report and package references
{
  #require(report)
  #session <- sessionInfo()
  #r <- report(session)
  #report saved to README.md (08/08/2025)
}

#these next sections run the search pipeline
#this code uses bracket notation {} to delineate discrete sections
#each section is labeled with a brief comment, 
#section contents are described within each section
#code is formatted to read .RDS files representing 
#intermediate/processed data objects-- uncommenting required to re-run


# Step 1: Load the ranges dataset and prepare it
load("/Users/cotinga/Downloads/ranges_4-16-22_multipolygons")

# Inspect the dataset
class(ranges)
str(ranges)

#Load and clean data
{
# Create a new column with underscores instead of spaces
ranges$sci_name_underscored <- gsub(" ", "_", ranges$sci_name)

# Step 2: Perform GBIF name matching for sampled_cv$phylo
# -------------------------------------------------------
# Assume sampled_cv is your existing dataset
# Inspect the structure of sampled_cv

sampled_cv_shift_metrics<-readRDS('sampled_cv_shift_metrics_8_08_25.RDS')
sampled_cv<-sampled_cv_shift_metrics #adding due to data loading simplification
#sum(sampled_cv_shift_metrics$phylo == sampled_cv$phylo)

# Perform GBIF name matching for `sampled_cv$phylo`
# gbif_matches <- name_backbone_checklist(sampled_cv$phylo)
# saveRDS(gbif_matches, file='gbif_matches.RDS')
gbif_matches <- readRDS('gbif_matches.RDS')

# Add the GBIF canonical names to the sampled_cv dataset
sampled_cv$phylo_matched <- gbif_matches$canonicalName

# Step 3: Inspect FUZZY, HIGHERRANK, and NONE matches
# ----------------------------------------------------
# Inspect FUZZY matches
fuzzy_matches <- gbif_matches[gbif_matches$matchType == "FUZZY", ]
cat("Number of FUZZY matches:", nrow(fuzzy_matches), "\n")
print(fuzzy_matches)

# Inspect HIGHERRANK matches
higherrank_matches <- gbif_matches[gbif_matches$matchType == "HIGHERRANK", ]
cat("Number of HIGHERRANK matches:", nrow(higherrank_matches), "\n")
print(higherrank_matches)

# Inspect NONE matches
none_matches <- gbif_matches[gbif_matches$matchType == "NONE", ]
cat("Number of NONE matches:", nrow(none_matches), "\n")
print(none_matches)

# Step 4: Handle Problematic Matches (FUZZY, HIGHERRANK, and NONE)
# ----------------------------------------------------------------
# Define replacements for HIGHERRANK and NONE matches
replacement_mapping <- data.frame(
  original_name = c("Tangara_cyanoptera", "Camaroptera_undosa", "Bernieria_madagascariensis"),
  suggested_name = c("Stilpnia_cyanoptera", "Calamonastes_undosus", "Phyllastrephus_madagascariensis"),
  stringsAsFactors = FALSE
)

# Update the `sampled_cv$phylo` list with replacements
name_mapping <- data.frame(
  original_name = sampled_cv$phylo,
  suggested_name = sampled_cv$phylo,  # Default to original names
  stringsAsFactors = FALSE
)

for (i in 1:nrow(replacement_mapping)) {
  name_mapping$suggested_name[name_mapping$original_name == replacement_mapping$original_name[i]] <- replacement_mapping$suggested_name[i]
}

# Check how many names were replaced
replaced_count <- sum(name_mapping$original_name != name_mapping$suggested_name)
cat("Number of names replaced:", replaced_count, "\n")
print(head(name_mapping[name_mapping$original_name != name_mapping$suggested_name, ]))

# Step 5: Perform GBIF name matching for the updated suggested names
# -------------------------------------------------------------------
gbif_matches <- rgbif::name_backbone_checklist(name_mapping$suggested_name)
#saveRDS(gbif_matches, file="updated_gbif_matches.RDS")
sum(is.na(gbif_matches$verbatim_name))
length(unique(gbif_matches$speciesKey))
}

#GBIF processing/filering
{
# Step 6: Perform GBIF name matching for ranges dataset
# ------------------------------------------------------
# Convert ranges$sci_name to character for GBIF matching
# ranges_gbif_matches <- name_backbone_checklist(as.character(ranges$sci_name))
# saveRDS(ranges_gbif_matches, file='ranges_gbif_matches.RDS')
ranges_gbif_matches <- readRDS(file='ranges_gbif_matches.RDS')

table(ranges_gbif_matches$matchType)
length(unique(ranges_gbif_matches$canonicalName))
length(unique(gbif_matches$canonicalName))
sum(gbif_matches$canonicalName %in% ranges_gbif_matches$canonicalName)
sum(gbif_matches$speciesKey %in% ranges_gbif_matches$speciesKey)

unmatched<-gbif_matches[!gbif_matches$speciesKey %in% ranges_gbif_matches$speciesKey,]
#unmatched$genus[1] %in% ranges_gbif_matches$genus
all(ranges$sci_name== ranges_gbif_matches$verbatim_name)

# Add speciesKey from ranges_gbif_matches to ranges
ranges_merged<-cbind(ranges, ranges_gbif_matches)

# Filter ranges to include only species matched in GBIF from sampled_cv
filtered_rows <- ranges_merged[ranges_merged$speciesKey %in% gbif_matches$speciesKey, ]

# Inspect the number of unique species in the filtered dataset
cat("Number of unique species in the final dataset:", length(unique(filtered_rows$speciesKey)), "\n")

# Filter for resident species only (seasonal = 1)
resident_species.1 <- filtered_rows[filtered_rows$seasonal == 1, ]
resident_species.12 <- filtered_rows[filtered_rows$seasonal %in% c(1, 2), ]
}

#capturing/processing the ranges data for the dataset
{
#year round ranges
{
  resident_species.1.merged <- merge_species_shapes_planar_parallel(resident_species.1, workers = 8)
  #valid<-unlist(pbmclapply(resident_species.1.merged$Shape, st_is_valid, mc.cores=8))
  resident_species.1.merged.deduplicated <- collapse_species_data_parallel(resident_species.1.merged, workers=8, verbose=T)
  #valid<-unlist(pbmclapply(resident_species.1.merged.deduplicated$Shape, st_is_valid, mc.cores=8))
  resident_species.1 <- resident_species.1.merged.deduplicated #overwrite
}
#year round ranges + breeding ranges
{
  resident_species.12.merged <- merge_species_shapes_planar_parallel(resident_species.12, workers = 8)
  #valid<-unlist(pbmclapply(resident_species.1.merged$Shape, st_is_valid, mc.cores=8))
  resident_species.12.merged.deduplicated <- collapse_species_data_parallel(resident_species.12.merged, workers=8, verbose=T)
  #valid<-unlist(pbmclapply(resident_species.1.merged.deduplicated$Shape, st_is_valid, mc.cores=8))
  resident_species.12 <- resident_species.12.merged.deduplicated
}

#code below for capturing/processing the ranges for Passeriformes
{
  passeriformes <- ranges_merged[ranges_merged$order == "Passeriformes",]
  passeriformes <- passeriformes[!is.na(passeriformes$speciesKey),] #filter out any rows with NA speciesKey at this stage
  passeriformes.filtered.1 <- passeriformes[passeriformes$seasonal == 1, ] #filter only for 'year round ranges'
  passeriformes.filtered.12 <- passeriformes[passeriformes$seasonal %in% c(1, 2), ] #filter only for 'year round + breeding ranges'
  # passeriformes.filtered.merged.1 <- merge_species_shapes_planar_parallel(passeriformes.filtered.1, workers = 8)
  # passeriformes.filtered.merged.deduplicated.1 <- collapse_species_data_parallel(passeriformes.filtered.merged.1, workers=8, verbose=T)
  # passeriformes.filtered.merged.12 <- merge_species_shapes_planar_parallel(passeriformes.filtered.12, workers = 8)
  # passeriformes.filtered.merged.deduplicated.12 <- collapse_species_data_parallel(passeriformes.filtered.merged.12, workers=8, verbose=T)
  #saveRDS(object = passeriformes.filtered.merged.deduplicated.1, file='passeriformes.filtered.merged.deduplicated.1.RDS')
  #saveRDS(object = passeriformes.filtered.merged.deduplicated.12, file='passeriformes.filtered.merged.deduplicated.12.RDS')
  passeriformes.filtered.merged.deduplicated.1 <- readRDS(file='passeriformes.filtered.merged.deduplicated.1.RDS')
  passeriformes.filtered.merged.deduplicated.12 <- readRDS(file='passeriformes.filtered.merged.deduplicated.12.RDS')
}

#cleaning up memory usage
plan(sequential)
gc()
plan(sequential)
gc()

}

#assign h3 cells to the dataset
{
  # # Configure the number of parallel workers
  # plan(multisession, workers = 8)  # Adjust the number of cores
  # # Apply the fallback function to all geometries in parallel
  # resident_species.1$h3_cells <- future_lapply(
  #   resident_species.1$Shape,
  #   function(geometry) assign_h3_approx(geometry, 3),
  #   future.seed = TRUE  # Ensures reproducible parallel-safe random numbers
  # )
  # saveRDS(resident_species.1, file='resident_species1_filter_res3.RDS')
  resident_species.1 <- readRDS(file='resident_species1_filter_res3.RDS')
  # Reset the parallel plan to sequential (optional cleanup step)
  #plan(sequential)
  
  # # # Configure the number of parallel workers
  # plan(multisession, workers = 8)  # Adjust the number of cores
  # # Apply the fallback function to all geometries in parallel
  # resident_species.12$h3_cells <- future_lapply(
  #   resident_species.12$Shape,
  #   function(geometry) assign_h3_approx(geometry, 3),
  #   future.seed = TRUE  # Ensures reproducible parallel-safe random numbers
  # )
  # saveRDS(resident_species.12, file='resident_species12_filter_res3.RDS')
  resident_species.12 <- readRDS(file='resident_species12_filter_res3.RDS')
  # Reset the parallel plan to sequential (optional cleanup step)
  # plan(sequential)
}

#assign h3 cells for all passerines
{
  # # # # Configure the number of parallel workers
  # plan(multisession, workers = 8)  # Adjust the number of cores
  # # Apply the fallback function to all geometries in parallel
  # passeriformes.filtered.merged.deduplicated.1$h3_cells <- future_lapply(
  #   passeriformes.filtered.merged.deduplicated.1$Shape,
  #   function(geometry) assign_h3_approx(geometry, 3),
  #   future.seed = TRUE  # Ensures reproducible parallel-safe random numbers
  # )
  # saveRDS(passeriformes.filtered.merged.deduplicated.1, file='passeriformes.filtered.1_res3.RDS')
  passeriformes.filtered.merged.deduplicated.1 <- readRDS(file='passeriformes.filtered.1_res3.RDS')
  # Reset the parallel plan to sequential (optional cleanup step)
  # plan(sequential)
  
  # # # Configure the number of parallel workers
  # plan(multisession, workers = 8)  # Adjust the number of cores
  # # Apply the fallback function to all geometries in parallel
  # passeriformes.filtered.merged.deduplicated.12$h3_cells <- future_lapply(
  #   passeriformes.filtered.merged.deduplicated.12$Shape,
  #   function(geometry) assign_h3_approx(geometry, 3),
  #   future.seed = TRUE  # Ensures reproducible parallel-safe random numbers
  # )
  # saveRDS(passeriformes.filtered.merged.deduplicated.12, file='passeriformes.filtered.12_res3.RDS')
  passeriformes.filtered.merged.deduplicated.12 <- readRDS(file='passeriformes.filtered.12_res3.RDS')
  # Reset the parallel plan to sequential (optional cleanup step)
  # plan(sequential)
  
  
}

#data checks
{
  any(is.na(resident_species.1$h3_cells))
  unique(unlist(resident_species.1$h3_cells))
  
  any(is.na(resident_species.12$h3_cells))
  unique(unlist(resident_species.12$h3_cells))
  
  sampled_cv_shift_metrics$log_tip_rate #looks good
  sampled_cv_shift_metrics$phylo #looks good 
  all(sort(name_mapping$original_name) == sort(sampled_cv_shift_metrics$phylo)) #the names match
}

# Check speciesKey mapping and generate lookup tables for branch metrics
{
  
  # Step 1: Map speciesKey to sampled_cv_shift_metrics using name_mapping and gbif_matches
  sampled_cv_shift_metrics$speciesKey <- gbif_matches$speciesKey[
    match(sampled_cv_shift_metrics$phylo, name_mapping$original_name)
  ]
  
  cat("Check speciesKeys in sampled_cv_shift_metrics:\n")
  print(head(sampled_cv_shift_metrics$speciesKey))
  cat("Number of NA speciesKeys in sampled_cv_shift_metrics:", sum(is.na(sampled_cv_shift_metrics$speciesKey)), "\n")
  
  # Step 2: Verify overlap between resident_species$speciesKey and sampled_cv_shift_metrics$speciesKey
  cat("Number of matching speciesKeys between resident_species and sampled_cv_shift_metrics:", 
      sum(resident_species.1$speciesKey %in% sampled_cv_shift_metrics$speciesKey), "\n")
  
  cat("Number of matching speciesKeys between resident_species and sampled_cv_shift_metrics:", 
      sum(resident_species.12$speciesKey %in% sampled_cv_shift_metrics$speciesKey), "\n")
  
  # Step 3: Create a lookup table for speciesKey to Tip_Phenotypic_Rate
  lookup_table.terminal <- setNames((sampled_cv_shift_metrics$log_tip_rate), sampled_cv_shift_metrics$speciesKey)
  #alternatively use weighted rates from shift events
  lookup_table.lineage.shifts <- setNames((sampled_cv_shift_metrics$log_lineage_rate_shifts), sampled_cv_shift_metrics$speciesKey)
  #alternatively use weighted rates from node events
  lookup_table.lineage.branches <- setNames((sampled_cv_shift_metrics$log_lineage_rate_branches), sampled_cv_shift_metrics$speciesKey)
  #lookup table for tipDR
  lookup_table.tipDR <- setNames(log(sampled_cv_shift_metrics$tipDR), sampled_cv_shift_metrics$speciesKey)
  
  #lookup table for residuals
  lookup_table.residuals <- setNames(
    lapply(seq_len(nrow(sampled_cv_shift_metrics)), function(i) {
      cols <- grep("^residuals\\.", names(sampled_cv_shift_metrics), value = TRUE)
      setNames(
        as.numeric(sampled_cv_shift_metrics[i, cols]),
        cols
      )
    }),
    sampled_cv_shift_metrics$speciesKey
  )
  
  # Debugging Step: Check the lookup table
  cat("Lookup table for log_tip_rate:\n")
  print(head(lookup_table.terminal))
}

# Assign values to spatial data frames based on lookup tables
{
resident_species.1$log_tip_rate <- lookup_table.terminal[as.character(resident_species.1$speciesKey)]
#resident_species.1$log_lineage_rate <- lookup_table.lineage.shifts[as.character(resident_species.1$speciesKey)] #change to lineage_rate
resident_species.1$log_lineage_rate <- lookup_table.lineage.branches[as.character(resident_species.1$speciesKey)] #change to lineage_rate
#for year round + breeding
resident_species.12$log_tip_rate <- lookup_table.terminal[as.character(resident_species.12$speciesKey)]
#resident_species.12$log_lineage_rate <- lookup_table.lineage.shifts[as.character(resident_species.12$speciesKey)] #change to lineage_rate
resident_species.12$log_lineage_rate <- lookup_table.lineage.branches[as.character(resident_species.12$speciesKey)] #change to lineage_rate

# Step 4a: Assign tipDR to resident_species based on speciesKey
resident_species.1$log_tipDR <- lookup_table.tipDR[as.character(resident_species.1$speciesKey)]
resident_species.12$log_tipDR <- lookup_table.tipDR[as.character(resident_species.12$speciesKey)]

#step5 : assign trait residuals to resident_species based on speciesKey
resident_species.1$residuals <- lookup_table.residuals[as.character(resident_species.1$speciesKey)]
resident_species.12$residuals <- lookup_table.residuals[as.character(resident_species.12$speciesKey)]
}

# Debugging Step: Check if log_tip_rate was correctly assigned
{
  cat("Check log_tip_rate in resident_species.1:\n")
  print(head(resident_species.1$log_tip_rate))
  
  cat("Check log_tip_rate in resident_species.123:\n")
  print(head(resident_species.12$log_tip_rate))
  
  cat("Number of NA log_tip_rate values in resident_species.1:", sum(is.na(resident_species.1$log_tip_rate)), "\n")
  cat("Number of NA log_tip_rate values in resident_species.12:", sum(is.na(resident_species.12$log_tip_rate)), "\n")
  
  # Step 5: Verify unmatched speciesKeys (if any)
  unmatched_keys.1 <- resident_species.1$speciesKey[is.na(resident_species.1$log_tip_rate)]
  unmatched_keys.12 <- resident_species.12$speciesKey[is.na(resident_species.12$log_tip_rate)]
  
  cat("Unmatched speciesKeys in resident_species.1:\n")
  print(unmatched_keys.1)
  
  cat("Unmatched speciesKeys in resident_species.1:\n")
  print(unmatched_keys.12)
}

# Verify that log_tip_rate values match for overlapping speciesKeys
{
  verification.1 <- data.frame(
    speciesKey = resident_species.1$speciesKey,
    resident_rate = resident_species.1$log_tip_rate,
    sampled_rate = lookup_table.terminal[as.character(resident_species.1$speciesKey)], 
    row.names = NULL
  )
  # Check if all assigned rates match the lookup rates
  mismatched.1 <- verification.1[verification.1$resident_rate.1 != verification.1$sampled_rate, ]
  cat("Number of mismatched log_tip_rate values:", nrow(mismatched.1), "\n")
  print(mismatched.1)
  length(unique(unlist(resident_species.1$h3_cells)))
  
  # Verify that Tip_Phenotypic_Rate values match for overlapping speciesKeys
  verification.12 <- data.frame(
    speciesKey = resident_species.12$speciesKey,
    resident_rate = resident_species.12$log_tip_rate,
    sampled_rate = lookup_table.terminal[as.character(resident_species.12$speciesKey)],
    row.names = NULL
  )
  # Check if all assigned rates match the lookup rates
  mismatched.12 <- verification.12[verification.12$resident_rate.12 != verification.12$sampled_rate, ]
  cat("Number of mismatched Tip_Phenotypic_Rate values:", nrow(mismatched.12), "\n")
  print(mismatched.12)
  length(unique(unlist(resident_species.12$h3_cells)))
}

#options(future.globals.maxSize= 4000*1024^2)
#aggregate metrics across h3 cells
{
resident_species.1.h3.metrics <- aggregate_h3_metrics(resident_species.1, cores=3, use_sqrt = F)
resident_species.12.h3.metrics <- aggregate_h3_metrics(resident_species.12, cores=3, use_sqrt= F)
}

#calculate dispersion/variance metrics per grid cell
{
# Set memory limit for parallel futures if needed
options(future.globals.maxSize = 4 * 1024^3)
resident_species.1.h3.mahalanobis <- aggregate_h3_disparity.mahalanobis(
  resident_species.1,
  cores = 6,
  method_order = c("solve")
)

# Set memory limit for parallel futures if needed
options(future.globals.maxSize = 4 * 1024^3)
resident_species.1.h3.centroid_dispersion <- aggregate_h3_dispersion_to_centroid(
  resident_species = resident_species.1,
  cores = 6,  # Adjust this to match the number of cores you want to use
  method_order = c("solve")  # Or use c("shrink", "solve") if corpcor is installed
)

# load fitted object data
min10.ic20.gic <- readRDS(file='new_bifrost/min10.ic20.gic.RDS')
min10.ic20.gic$model_fit_history$fits<-NULL #wipe out the history to save memory

resident_species.1.h3.divtimes <- aggregate_h3_divergence_times(
  resident_species = resident_species.1,
  phy = min10.ic20.gic$tree_no_uncertainty_untransformed,
  key_to_tip_map = sampled_cv_shift_metrics[, c("speciesKey", "phylo")],
  cores = 6
)

disparity.divergence.weighted <- summarize_disparity_vs_divergence(
  mahalanobis_data = resident_species.1.h3.mahalanobis,
  divergence_data = resident_species.1.h3.divtimes,
  plan_strategy = "multisession",
  cores = 8, bootstrap = T, n_boot=100, use_weighted = T, standardize = F, keep_model_fit = FALSE
)

disparity.divergence <- summarize_disparity_vs_divergence(
  mahalanobis_data = resident_species.1.h3.mahalanobis,
  divergence_data = resident_species.1.h3.divtimes,
  plan_strategy = "multisession",
  cores = 8, bootstrap = T, n_boot=100, use_weighted = F, standardize = F
)

#cleanup
plan(sequential)
gc()

kz_per_cell <- aggregate_phylo_signal(
  resident_species = resident_species.1,
  phy = min10.ic20.gic$tree_no_uncertainty_untransformed,
  key_to_tip_map = sampled_cv_shift_metrics[, c("speciesKey", "phylo")],
  trait_column = "residuals",
  min_species = 2,
  iter = 1000,
  lambda_method = 'front', cores = 8
)

results_vcv <- local_vcv(
  resident_species = resident_species.1,                  # same input as before
  phy = min10.ic20.gic$tree_no_uncertainty_untransformed, # same phylogeny
  key_to_tip_map = sampled_cv_shift_metrics[, c("speciesKey", "phylo")], # same mapping
  n_min_species = 2,                                      # optional threshold for species per cell
  plan_strategy = "multisession",                         # parallel backend
  cores = 8                                               # same number of cores
)
hist(results_vcv$metrics_df$modularity)



}

#check weights
{
first_h3_cell <- resident_species.1.h3.metrics$weighting_df$h3_cell[1]
print(first_h3_cell)
subset_df <- resident_species.1.h3.metrics$weighting_df %>% filter(h3_cell == first_h3_cell)
}

#generate metrics for all passeriformes
{
#using 3 cores to use less memory
passeriformes.filtered.1.metrics <- aggregate_h3_metrics_counts(passeriformes.filtered.merged.deduplicated.1, cores = 2, use_sqrt = F)
gc(); plan(sequential);
passeriformes.filtered.12.metrics <- aggregate_h3_metrics_counts(passeriformes.filtered.merged.deduplicated.12, cores = 2, use_sqrt = F)
gc(); plan(sequential);
passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted_log <- log(passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted)
gc(); plan(sequential);
passeriformes.filtered.12.metrics$metrics_df$species_key_count_range_weighted_log <- log(passeriformes.filtered.12.metrics$metrics_df$species_key_count_range_weighted)
gc(); plan(sequential);

#now add sampling fraction to dataset (overwriting)
resident_species.1.h3.metrics.backup <- resident_species.1.h3.metrics
{
  resident_species.1.h3.metrics$metrics_df <- merge(resident_species.1.h3.metrics$metrics_df, 
                                                    passeriformes.filtered.1.metrics$metrics_df[, c("h3_cell", "species_key_count", "species_key_count_range_weighted")], 
                                                    by = "h3_cell", all.x = TRUE, suffixes = c("", ".passerines"))
  resident_species.1.h3.metrics$metrics_df$sampling_frac <- resident_species.1.h3.metrics$metrics_df$species_key_count/resident_species.1.h3.metrics$metrics_df$species_key_count.passerines
  resident_species.1.h3.metrics$metrics_df$sampling_frac_asin <- asin(sqrt(resident_species.1.h3.metrics$metrics_df$sampling_frac))
  resident_species.1.h3.metrics$metrics_df$sampling_frac.range_weighted <- resident_species.1.h3.metrics$metrics_df$species_key_count_range_weighted/resident_species.1.h3.metrics$metrics_df$species_key_count_range_weighted.passerines
  resident_species.1.h3.metrics$metrics_df$sampling_frac.range_weighted_asin <- asin(sqrt(resident_species.1.h3.metrics$metrics_df$sampling_frac.range_weighted))
}
resident_species.12.h3.metrics.backup <- resident_species.12.h3.metrics
{
  resident_species.12.h3.metrics$metrics_df <- merge(resident_species.12.h3.metrics$metrics_df, 
                                                     passeriformes.filtered.12.metrics$metrics_df[, c("h3_cell", "species_key_count", "species_key_count_range_weighted")], 
                                                     by = "h3_cell", all.x = TRUE, suffixes = c("", ".passerines"))
  resident_species.12.h3.metrics$metrics_df$sampling_frac <- resident_species.12.h3.metrics$metrics_df$species_key_count/resident_species.12.h3.metrics$metrics_df$species_key_count.passerines
  resident_species.12.h3.metrics$metrics_df$sampling_frac_asin <- asin(sqrt(resident_species.12.h3.metrics$metrics_df$sampling_frac))
  resident_species.12.h3.metrics$metrics_df$sampling_frac.range_weighted <- resident_species.12.h3.metrics$metrics_df$species_key_count_range_weighted/resident_species.12.h3.metrics$metrics_df$species_key_count_range_weighted.passerines
  resident_species.12.h3.metrics$metrics_df$sampling_frac.range_weighted_asin <- asin(sqrt(resident_species.12.h3.metrics$metrics_df$sampling_frac.range_weighted))
}
}

# Precompute polygons
{
polygons_sf.1 <- precompute_h3_polygons(resident_species.1.h3.metrics, num_cores = 8)
polygons_sf.12 <- precompute_h3_polygons(resident_species.12.h3.metrics, num_cores = 8)
polygons_passerines_sf.1 <- precompute_h3_polygons(passeriformes.filtered.1.metrics, num_cores = 8)
polygons_passerines_sf.12 <- precompute_h3_polygons(passeriformes.filtered.12.metrics, num_cores = 8) ##check this -- might be a bug

#cleanup
plan(sequential)
gc()
plan(sequential)
gc()
}

#test plots
{
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.metrics,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "mean_tip_rate"    # Column to visualize
  )
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.metrics,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "tip_rate_range_weighted_mean"    # Column to visualize
  )
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.metrics,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "mean_lineage_rate"    # Column to visualize
  )
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.metrics,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "mean_lineage_rate_range_weighted"    # Column to visualize
  )
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.metrics,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "mean_tipDR"    # Column to visualize
  )
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.metrics,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "tipDR_range_weighted_mean"    # Column to visualize
  )
  
  # # Plot the precomputed polygons
  # plot_h3_metrics_global(
  #   results = resident_species.1.h3.vcv.disparity,                # Results containing metrics_df
  #   polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
  #   metric_column = "vcv_trace_unweighted"    # Column to visualize
  # )
  # 
  # # Plot the precomputed polygons
  # plot_h3_metrics_global(
  #   results = resident_species.1.h3.vcv.disparity,                # Results containing metrics_df
  #   polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
  #   metric_column = "vcv_trace_range_weighted",    # Column to visualize
  #   log_transform=F)
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = results_vcv,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "det_Pinv",    # Column to visualize
    log_transform=T)
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.mahalanobis,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "mean_mahalanobis_disparity", log_transform = F    # Column to visualize
  ) /
    plot_h3_metrics_global(
      results = resident_species.1.h3.mahalanobis,                # Results containing metrics_df
      polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
      metric_column = "median_mahalanobis_disparity", log_transform = F    # Column to visualize
    ) /
    plot_h3_metrics_global(
      results = resident_species.1.h3.centroid_dispersion,                # Results containing metrics_df
      polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
      metric_column = "mean_mahalanobis_distance_to_centroid", log_transform = F    # Column to visualize
    )
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.divtimes,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "max_divergence_time", log_transform = F    # Column to visualize
  )
  
  # Plot the precomputed polygons
  plot_h3_metrics_global(
    results = resident_species.1.h3.divtimes,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "mean_divergence_time", log_transform = F    # Column to visualize
  )
  
  disparity.divergence$metrics_df <- disparity.divergence$summary_df
  plot_h3_metrics_global(
    results = disparity.divergence,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "standardized_residual", log_transform = F    # Column to visualize
  )
  
  # plot_h3_metrics_global(
  #   results = results_mean,                # Results containing metrics_df
  #   polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
  #   metric_column = "empirical_disp", log_transform = T    # Column to visualize
  # )
  
  disparity.divergence.weighted$metrics_df <- disparity.divergence.weighted$summary_df
  plot_h3_metrics_global(
    results = disparity.divergence.weighted,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "prop_extreme_residuals", log_transform = F    # Column to visualize
  )
  
  kz_per_cell$metrics_df<-kz_per_cell
  #k_per_cell$metrics_df[k_per_cell$metrics_df$Z < 0, ]
  plot_h3_metrics_global(
    results = kz_per_cell,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "Z", log_transform = F    # Column to visualize
  )
  plot_h3_metrics_global(
    results = kz_per_cell,                # Results containing metrics_df
    polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
    metric_column = "K", log_transform = F    # Column to visualize
  )
  
  
}

#test plots for disparity and phylogenetic signal
{
  
  # Example call
  plot_disparity_metrics_by_latitude(
    disparity_results = disparity.divergence,
    kz_results = kz_per_cell,
    add_kz_metrics = TRUE,
    dispersion_results = resident_species.1.h3.mahalanobis$metrics_df,
    add_dispersion_metrics = TRUE,
    metrics = c("standardized_residual", "K", "median_mahalanobis_disparity"),
    zero_line_metrics = c("Z"),
    one_line_metrics  = c("K"),
    point_size = 0.75,
    add_loess = FALSE,
    bin_size = 15,
    step_size = 5,
    p_threshold = 0.05, color_by_region = T, point_alpha = 0.25, legend_alpha = 0.8
  )
  
  # Example call
  plot_disparity_metrics_by_latitude(
    disparity_results = disparity.divergence,
    kz_results = kz_per_cell,
    add_kz_metrics = TRUE,
    dispersion_results = resident_species.1.h3.mahalanobis$metrics_df,
    add_dispersion_metrics = TRUE,
    metrics = c("K"),
    zero_line_metrics = c("Z"),
    one_line_metrics  = c("K"),
    point_size = 0.75,
    add_loess = FALSE,
    bin_size = 15,
    step_size = 5,
    p_threshold = 0.05, color_by_region = T, point_alpha = 0.25, legend_alpha = 0.8,
    n_species_threshold = 10, low_n_species_color = 'yellow'
  ) #+ geom_vline(xintercept = 68, linetype = "solid", color = "black")
    #geom_vline(xintercept = - 40, linetype = "solid", color = "black")
  
}

#plotting test for vcv metrics
{
  
  
  plot_vcv_metrics_by_latitude_v2(
    results_vcv,
    metrics         = c("trace", "effective_dimensionality", "modularity"),
    source          = "estimated",
    add_loess       = F,
    point_size      = 0.75,
    point_alpha     = 0.25,
    zero_line_metrics = c("effective_dim_per_species"),
    one_line_metrics  = c("lambda1_prop"),
    color_by_region = TRUE,
    legend_alpha = 0.8
  )
  
  
}

#build disparity and vcv patchworks
{
  dp_plot  <- plot_disparity_metrics_by_latitude(
    disparity_results     = disparity.divergence,
    kz_results            = kz_per_cell,
    add_kz_metrics        = TRUE,
    dispersion_results    = resident_species.1.h3.mahalanobis$metrics_df,
    add_dispersion_metrics= TRUE,
    metrics               = c("standardized_residual", "median_mahalanobis_disparity", "K"),
    zero_line_metrics     = c("Z"),
    one_line_metrics      = c("K"),
    point_size            = 0.75,
    add_loess             = FALSE,
    bin_size              = 10,
    step_size             = 2.5,
    p_threshold           = 0.05,
    color_by_region       = TRUE,
    point_alpha           = 0.25,
    legend_alpha          = 0.8, p_adjust_method = 'BY',
    n_species_threshold = 10, low_n_species_color = 'yellow'
  )
  
  vcv_plot <- plot_vcv_metrics_by_latitude_v2(
    results_vcv           = results_vcv,
    metrics               = c("modularity", "effective_dimensionality", "trace"),
    source                = "estimated",
    add_loess             = FALSE,
    point_size            = 0.75,
    point_alpha           = 0.25,
    zero_line_metrics     = c("effective_dim_per_species"),
    one_line_metrics      = c("lambda1_prop"),
    color_by_region       = TRUE,
    legend_alpha          = 0.8
  ) & 
    theme(legend.position = "none") 
  
  # 2) re-wrap vcv sub‐plots, keeping legend only on the first one
  vcv_plot <- wrap_plots(
    vcv_plot[[1]] + xlim(-85,85),
    vcv_plot[[2]] + xlim(-85,85),
    vcv_plot[[3]] + xlim(-85,85),
    ncol = 3
  )
  
  # 2) re-wrap vcv sub‐plots, keeping legend only on the first one
  dp_plot <- wrap_plots(
    dp_plot[[1]] + theme(legend.position = "none") + xlim(-85,85),
    dp_plot[[2]] + theme(legend.position = "none") + xlim(-85,85),
    dp_plot[[3]] + theme(legend.position = "right") + xlim(-85,85),
    ncol = 3
  )
  
}

#plotting for figure 4/supp fig 16
{
quartz(width=11*1.1, height=5*1.1, type='pdf', file = 'extinction.pdf')
dp_plot / vcv_plot
dev.off()

quartz(height=3*2.1, width=(9*1.5)/1.6, type='pdf', file = 'Figure4.pdf')
(
  ((vcv_plot[[1]] + xlim(-75, 75) + theme_classic() + theme(legend.position = "none")) |
     (dp_plot[[1]] + xlim(-75, 75) + theme_classic() + theme(legend.position = "right"))) /
    ((vcv_plot[[2]] + xlim(-75, 75) + theme_classic()+  theme(legend.position = "none")) |
       (dp_plot[[3]] + xlim(-75, 75) + theme_classic() + theme(legend.position = "none")))
) + 
  plot_annotation(tag_levels = "A")
dev.off()

quartz(height=3*1.15, width=(9*1.5)/1.65, type='pdf', file = 'SuppFig16.pdf')
(
  
  (dp_plot[[2]] + xlim(-75, 75) + theme_classic() + theme(legend.position = "none")) |
    (vcv_plot[[3]] + xlim(-75, 75) + theme_classic()+  theme(legend.position = "right"))
) + 
  plot_annotation(tag_levels = "A")
dev.off()
}

#more test plots
{
#lambert cylindrical equal area
# Plot the precomputed polygons
plot_h3_metrics_global(
  results = resident_species.1.h3.metrics,                # Results containing metrics_df
  polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
  metric_column = "mean_lineage_rate_range_weighted" ,   # Column to visualize
  crs = "+proj=cea +lat_ts=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
)

#robinson
plot_h3_metrics_global(
  results = resident_species.1.h3.metrics,                # Results containing metrics_df
  polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
  metric_column = "mean_lineage_rate_range_weighted" ,   # Column to visualize
  crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  ocean_color = alpha('lightblue', 0.25), border_thickness = 0.25
)

#robinson
plot_h3_metrics_global(
  results = resident_species.1.h3.metrics,                # Results containing metrics_df
  polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
  metric_column = "mean_tipDR" ,   # Column to visualize
  crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  ocean_color = alpha('lightblue', 0.25), border_thickness = 0.25
)
}

#passerine diversity plots (testing)
{
  
#robinson (all passerine diversity)
plot_h3_metrics_global(
  results = passeriformes.filtered.1.metrics,                # Results containing metrics_df
  polygons_sf = polygons_passerines_sf.1,     # Precomputed H3 polygons
  metric_column = "species_key_count" ,   # Column to visualize
  crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  ocean_color = alpha('lightblue', 0.25), border_thickness = 0.25
)
  
#robinson (all passerine diversity)
tmp1<-plot_h3_metrics_global(
  results = passeriformes.filtered.1.metrics,                # Results containing metrics_df
  polygons_sf = polygons_passerines_sf.1,     # Precomputed H3 polygons
  metric_column = "species_key_count_range_weighted" ,   # Column to visualize
  crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  ocean_color = alpha('lightblue', 0.25), border_thickness = 0.25
)
#robinson (all passerine diversity)
tmp2<-plot_h3_metrics_global(
  results = passeriformes.filtered.1.metrics,                # Results containing metrics_df
  polygons_sf = polygons_passerines_sf.1,     # Precomputed H3 polygons
  metric_column = "species_key_count_range_weighted_log" ,   # Column to visualize
  crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  ocean_color = alpha('lightblue', 0.25), border_thickness = 0.25
)
plot(tmp1 / tmp2)

median(passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted)
max(passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted)
min(passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted)

#robinson (all passerine diversity)
sampling_frac.1<-plot_h3_metrics_global(
  results = resident_species.1.h3.metrics,                # Results containing metrics_df
  polygons_sf = polygons_sf.1,     # Precomputed H3 polygons
  metric_column = "sampling_frac_asin" ,   # Column to visualize
  crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  ocean_color = alpha('lightblue', 0.25), border_thickness = 0.25
)

#robinson (all passerine diversity)
sampling_frac.12 <- plot_h3_metrics_global(
  results = resident_species.12.h3.metrics,                # Results containing metrics_df
  polygons_sf = polygons_sf.12,     # Precomputed H3 polygons
  metric_column = "sampling_frac_asin" ,   # Column to visualize
  crs = "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs",
  ocean_color = alpha('lightblue', 0.25), border_thickness = 0.25
)

sampling_frac.1 / sampling_frac.12
 
  
}

#setting up baseline plots
#year round
{
  
  #mean_tip_rate.1 
  {
    
    custom_palette.mean_tip_rate.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_tip_rate,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_tip_rate.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_tip_rate",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_tip_rate.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_tip_rate.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_tip_rate,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_tip_rate.1$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_tip_rate.1$breaks,  # Breaks from custom_palette
      xlab = "Mean tip rate",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    ) + ylim(0, 1050) + xlim(-9.5, -7)
    
    mean_tip_rate.1.combined <- mean_tip_rate.1 / ((plot_spacer() | mean_tip_rate.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Mean Tip Rates")
    
  }
  
  #tip_rate_range_weighted_mean.1
  {
    custom_palette.tip_rate_range_weighted_mean.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$tip_rate_range_weighted_mean,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    tip_rate_range_weighted_mean.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "tip_rate_range_weighted_mean",         # The column you want to visualize
      custom_palette = alpha(custom_palette.tip_rate_range_weighted_mean.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    tip_rate_range_weighted_mean.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$tip_rate_range_weighted_mean,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.tip_rate_range_weighted_mean.1$colors,  # Colors from custom_palette
      breaks = custom_palette.tip_rate_range_weighted_mean.1$breaks,  # Breaks from custom_palette
      xlab = "Mean tip rate, range weighted",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    ) + ylim(0, 1150) + xlim(-10, -6)
    
    tip_rate_range_weighted_mean.1.combined <- tip_rate_range_weighted_mean.1 / ((plot_spacer() | tip_rate_range_weighted_mean.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Tip Rate Range Weighted Mean")
    
    
  }
  
  #mean_lineage_rate.1
  {
    custom_palette.mean_lineage_rate.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_lineage_rate,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_lineage_rate.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_lineage_rate",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_lineage_rate.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_lineage_rate.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_lineage_rate,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_lineage_rate.1$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_lineage_rate.1$breaks,  # Breaks from custom_palette
      xlab = "Mean lineage rate",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    ) + ylim(0, 1000)
    
    mean_lineage_rate.1.combined <- mean_lineage_rate.1 / ((plot_spacer() | mean_lineage_rate.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Mean Lineage Rates")
    
  }
  
  #mean_lineage_rate_range_weighted.1
  {
    custom_palette.mean_lineage_rate_range_weighted.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_lineage_rate_range_weighted,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_lineage_rate_range_weighted.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_lineage_rate_range_weighted",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_lineage_rate_range_weighted.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_lineage_rate_range_weighted.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_lineage_rate_range_weighted,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_lineage_rate_range_weighted.1$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_lineage_rate_range_weighted.1$breaks,  # Breaks from custom_palette
      xlab = "Mean lineage rate, range weighted",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    ) + ylim(0, 1000)
    
    mean_lineage_rate_range_weighted.1.combined <- mean_lineage_rate_range_weighted.1 / ((plot_spacer() | mean_lineage_rate_range_weighted.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Mean Lineage Rate Range Weighted")
    
    
  }
  
  #mean_tipDR.1
  {
    custom_palette.mean_tipDR.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_tipDR,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_tipDR.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_tipDR",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_tipDR.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_tipDR.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$mean_tipDR,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_tipDR.1$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_tipDR.1$breaks,  # Breaks from custom_palette
      xlab = "Mean tipDR",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    ) + ylim(0,1150) + xlim(-3, 1.5)
    
    mean_tipDR.1.combined <- mean_tipDR.1 / ((plot_spacer() | mean_tipDR.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of tipDR")
    mean_tipDR.1.combined
    
    
  }
  
  #tipDR_range_weighted_mean.1
  {
    custom_palette.tipDR_range_weighted_mean.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$tipDR_range_weighted_mean,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    tipDR_range_weighted_mean.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "tipDR_range_weighted_mean",         # The column you want to visualize
      custom_palette = alpha(custom_palette.tipDR_range_weighted_mean.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    tipDR_range_weighted_mean.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$tipDR_range_weighted_mean,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.tipDR_range_weighted_mean.1$colors,  # Colors from custom_palette
      breaks = custom_palette.tipDR_range_weighted_mean.1$breaks,  # Breaks from custom_palette
      xlab = "Mean tipDR, range weighted",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    ) + ylim(0, 1300) + xlim(-3, 1)
    
    tipDR_range_weighted_mean.1.combined <- tipDR_range_weighted_mean.1 / ((plot_spacer() | tipDR_range_weighted_mean.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of tipDR, Range Weighted")
    tipDR_range_weighted_mean.1.combined
    
    
  }
  
  #species_key_count.1
  {
    
    custom_palette.species_key_count.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$species_key_count,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    # Use the modified function with the custom palette
    sp.count.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "species_key_count",         # The column you want to visualize
      custom_palette = alpha(custom_palette.species_key_count.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    sp.count.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$species_key_count,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.species_key_count.1$colors,  # Colors from custom_palette
      breaks = custom_palette.species_key_count.1$breaks,  # Breaks from custom_palette
      ylab = "Cell Counts",          # X-axis label
      xlab = "Species Counts",                 # Y-axis label
      title = "Histogram of Species Counts",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.75, cex.main = 1.75, cex.lab = 1.5         # Do not use a log scale
    ) + ylim(0, 1500) + xlim(0, 250)
    
    sp.count.1.combined <- sp.count.1 / ((plot_spacer() | sp.count.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Species Count")
    
    
  }
  
  #sampling_frac.1
  {
    custom_palette.sampling_frac.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$sampling_frac,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "linear",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    # Use the modified function with the custom palette
    sampling_frac.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "sampling_frac",         # The column you want to visualize
      custom_palette = alpha(custom_palette.sampling_frac.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    sampling_frac.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$sampling_frac,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.sampling_frac.1$colors,  # Colors from custom_palette
      breaks = custom_palette.sampling_frac.1$breaks,  # Breaks from custom_palette
      ylab = "Cell Counts",          # X-axis label
      xlab = "Species Sampling Fraction",                 # Y-axis label
      title = "Histogram of Species Counts",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.75, cex.main = 1.75, cex.lab = 1.5         # Do not use a log scale
    ) + ylim(0, 2000) + xlim(0, 1)
    
    sampling_frac.1.combined <- sampling_frac.1 / ((plot_spacer() | sampling_frac.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Sampling Fraction")
    
    
  }
  
  #sampling_frac_asin.1
  {
    
    custom_palette.sampling_frac_asin.1 <- assignRateColors(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$sampling_frac_asin,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    sampling_frac_asin.1<-plot_h3_metrics_global.custom.raster(
      results = resident_species.1.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.1,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "sampling_frac_asin",         # The column you want to visualize
      custom_palette = alpha(custom_palette.sampling_frac_asin.1$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 15, 
      lon_graticule_spacing = 15, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    sampling_frac_asin.1.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.1.h3.metrics$metrics_df$sampling_frac_asin,
        resident_species.1.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.sampling_frac_asin.1$colors,  # Colors from custom_palette
      breaks = custom_palette.sampling_frac_asin.1$breaks,  # Breaks from custom_palette
      ylab = "Cell Counts",          # X-axis label
      xlab = "Species Sampling Fraction\n(arcsine square root)",                 # Y-axis label
      title = "Histogram of Species Counts",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.75, cex.main = 1.75, cex.lab = 1.5         # Do not use a log scale
    ) + ylim(0, 1700) + xlim(0.1, 1.6)
    
    sampling_frac_asin.1.combined <- sampling_frac_asin.1 / ((plot_spacer() | sampling_frac_asin.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Sampling Fraction (Arcsine Transformed)")
    
  }
  
  
  plot(resident_species.1.h3.metrics$metrics_df$tip_rate_range_weighted_mean~
       resident_species.1.h3.metrics$metrics_df$mean_lineage_rate_range_weighted, 
       ylab = "Tip-rate grid cell mean estimates",
       xlab = "Lineage-rate grid cell mean estimates", main = "Tip vs Lineage-rate in spatial analyses")
  
  quartz(file = "Lineage_vs_Tip_rates.pdf", height = 8, width = 8, type = 'pdf')
  with(resident_species.1.h3.metrics$metrics_df, {
    # initial fit to get studentized residuals (lineage ~ tip)
    fit_all <- lm(mean_lineage_rate_range_weighted ~ tip_rate_range_weighted_mean)
    stud    <- rstudent(fit_all)
    keep    <- abs(stud) <= 3
    
    # plot filtered data: x = tip-rate, y = lineage-rate
    plot(
      tip_rate_range_weighted_mean[keep],
      mean_lineage_rate_range_weighted[keep],
      xlab = "Tip-rate grid cell mean estimates",
      ylab = "Lineage-rate grid cell mean estimates",
      main = "Lineage vs Tip-rate in spatial analyses\n(| rstudent | > 3 removed)",
      pch  = 21,
      bg   = adjustcolor("grey", alpha.f = 0.50),
      col  = "black",
      lwd  = 0.01,
      cex  = 2.00
    )
    
    # refit on filtered data and add regression line (lineage ~ tip)
    fit <- lm(mean_lineage_rate_range_weighted[keep] ~ tip_rate_range_weighted_mean[keep])
    abline(fit, col = "red", lwd = 2, lty = 2)
    
    # correlation and R^2 on filtered data
    r   <- cor(tip_rate_range_weighted_mean[keep],
               mean_lineage_rate_range_weighted[keep],
               use = "complete.obs")
    r2  <- summary(fit)$r.squared
    a   <- coef(fit)[1]
    b   <- coef(fit)[2]
    
    eqn_text <- bquote(
      y == .(round(a, 3)) + .(round(b, 3)) * x ~
        "(" * R^2 == .(round(r2, 3)) * "," ~ r == .(round(r, 3)) * ")"
    )
    
    legend("topleft", legend = as.expression(eqn_text), bty = "n")
  })
  dev.off()
  
}

#patchwork arangement of baseline plots
{
  
  pdf(file = 'mean_tip_rates.combined.pdf', height=10, width=7.5*2)
  (mean_tip_rate.1.combined | tip_rate_range_weighted_mean.1.combined) +
    plot_annotation(tag_levels = 'A') &
    theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  dev.off()
  
  pdf(file = 'mean_lineage_rates.combined.pdf', height=10, width=7.5*2)
  (mean_lineage_rate.1.combined | mean_lineage_rate_range_weighted.1.combined) +
    plot_annotation(tag_levels = 'A') &
    theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  dev.off()
  
  
  pdf(file = 'mean_tipDRs.combined.pdf', height=10, width=7.5*2)
  (mean_tipDR.1.combined | tipDR_range_weighted_mean.1.combined) +
    plot_annotation(tag_levels = 'A') &
    theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  dev.off()
  
  
  pdf(file = 'sampling_fraction.combined.pdf', height=10, width=7.5*3)
  (sp.count.1.combined | sampling_frac.1.combined | sampling_frac_asin.1.combined) +
    plot_annotation(tag_levels = 'A') &
    theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  dev.off()
  
  
}

#year round + breeding 
{
  
  #mean_tip_rate.12
  {
    
    custom_palette.mean_tip_rate.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_tip_rate,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_tip_rate.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_tip_rate",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_tip_rate.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_tip_rate.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_tip_rate,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_tip_rate.12$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_tip_rate.12$breaks,  # Breaks from custom_palette
      xlab = "Mean tip rate",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    mean_tip_rate.12.combined <- mean_tip_rate.12 / ((plot_spacer() | mean_tip_rate.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Mean Tip Rates")
    
  }
  
  #tip_rate_range_weighted_mean.12
  {
    custom_palette.tip_rate_range_weighted_mean.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$tip_rate_range_weighted_mean,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    tip_rate_range_weighted_mean.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "tip_rate_range_weighted_mean",         # The column you want to visualize
      custom_palette = alpha(custom_palette.tip_rate_range_weighted_mean.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    tip_rate_range_weighted_mean.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$tip_rate_range_weighted_mean,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.tip_rate_range_weighted_mean.12$colors,  # Colors from custom_palette
      breaks = custom_palette.tip_rate_range_weighted_mean.12$breaks,  # Breaks from custom_palette
      xlab = "Mean tip rate, range weighted",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    tip_rate_range_weighted_mean.12.combined <- tip_rate_range_weighted_mean.12 / ((plot_spacer() | tip_rate_range_weighted_mean.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Tip Rate Range Weighted Mean")
    
    
  }
  
  #mean_lineage_rate.12
  {
    custom_palette.mean_lineage_rate.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_lineage_rate,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_lineage_rate.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_lineage_rate",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_lineage_rate.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_lineage_rate.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_lineage_rate,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_lineage_rate.12$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_lineage_rate.12$breaks,  # Breaks from custom_palette
      xlab = "Mean lineage rate",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    mean_lineage_rate.12.combined <- mean_lineage_rate.12 / ((plot_spacer() | mean_lineage_rate.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Mean Lineage Rates")
    
  }
  
  #mean_lineage_rate_range_weighted.12
  {
    
    custom_palette.mean_lineage_rate_range_weighted.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_lineage_rate_range_weighted,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_lineage_rate_range_weighted.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_lineage_rate_range_weighted",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_lineage_rate_range_weighted.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_lineage_rate_range_weighted.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_lineage_rate_range_weighted,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_lineage_rate_range_weighted.12$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_lineage_rate_range_weighted.12$breaks,  # Breaks from custom_palette
      xlab = "Mean lineage rate, range weighted",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    mean_lineage_rate_range_weighted.12.combined <- mean_lineage_rate_range_weighted.12 / ((plot_spacer() | mean_lineage_rate_range_weighted.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Mean Lineage Rate Range Weighted")
    
    
  }
  
  #mean_tipDR.12
  {
    custom_palette.mean_tipDR.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_tipDR,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    mean_tipDR.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "mean_tipDR",         # The column you want to visualize
      custom_palette = alpha(custom_palette.mean_tipDR.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    mean_tipDR.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$mean_tipDR,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.mean_tipDR.12$colors,  # Colors from custom_palette
      breaks = custom_palette.mean_tipDR.12$breaks,  # Breaks from custom_palette
      xlab = "Mean tipDR",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    mean_tipDR.12.combined <- mean_tipDR.12 / ((plot_spacer() | mean_tipDR.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of tipDR")
    mean_tipDR.12.combined
    mean_tipDR.1.combined
    
    
  }
  
  #tipDR_range_weighted_mean.1
  {
    custom_palette.tipDR_range_weighted_mean.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$tipDR_range_weighted_mean,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    tipDR_range_weighted_mean.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "tipDR_range_weighted_mean",         # The column you want to visualize
      custom_palette = alpha(custom_palette.tipDR_range_weighted_mean.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    tipDR_range_weighted_mean.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$tipDR_range_weighted_mean,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.tipDR_range_weighted_mean.12$colors,  # Colors from custom_palette
      breaks = custom_palette.tipDR_range_weighted_mean.12$breaks,  # Breaks from custom_palette
      xlab = "Mean tipDR, range weighted",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Rate Variation",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    tipDR_range_weighted_mean.12.combined <- tipDR_range_weighted_mean.12 / ((plot_spacer() | tipDR_range_weighted_mean.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of tipDR, Range Weighted")
    tipDR_range_weighted_mean.12.combined
    
    
  }
  
  #species_key_count.12
  {
    
    custom_palette.species_key_count.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$species_key_count,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    # Use the modified function with the custom palette
    sp.count.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "species_key_count",         # The column you want to visualize
      custom_palette = alpha(custom_palette.species_key_count.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    sp.count.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$species_key_count,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.species_key_count.12$colors,  # Colors from custom_palette
      breaks = custom_palette.species_key_count.12$breaks,  # Breaks from custom_palette
      xlab = "Counts",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Species Counts",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    sp.count.12.combined <- sp.count.12 / ((plot_spacer() | sp.count.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Species Count")
    
    
  }
  
  #sampling_frac.12
  {
    
    custom_palette.sampling_frac.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$sampling_frac,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "linear",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    # Use the modified function with the custom palette
    sampling_frac.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "sampling_frac",         # The column you want to visualize
      custom_palette = alpha(custom_palette.sampling_frac.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    sampling_frac.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$sampling_frac,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.sampling_frac.12$colors,  # Colors from custom_palette
      breaks = custom_palette.sampling_frac.12$breaks,  # Breaks from custom_palette
      xlab = "Sampling Fraction",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Sampling Fraction",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    sampling_frac.12.combined <- sampling_frac.12 / ((plot_spacer() | sampling_frac.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Sampling Fraction")
    
    
  }
  
  #sampling_frac_asin.12
  {
    
    custom_palette.sampling_frac_asin.12 <- assignRateColors(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$sampling_frac_asin,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 20,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    sampling_frac_asin.12<-plot_h3_metrics_global.custom.raster(
      results = resident_species.12.h3.metrics,  # Your results object
      polygons_sf = polygons_sf.12,             # The polygons_sf object (precomputed H3 polygons)
      metric_column = "sampling_frac_asin",         # The column you want to visualize
      custom_palette = alpha(custom_palette.sampling_frac_asin.12$colors, 1),  # Use the colors from your custom_palette
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
      border_thickness = 0.25, lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
    )
    
    sampling_frac_asin.12.hist<-ratesHistogramGG(
      rates = setNames(
        resident_species.12.h3.metrics$metrics_df$sampling_frac_asin,
        resident_species.12.h3.metrics$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.sampling_frac_asin.12$colors,  # Colors from custom_palette
      breaks = custom_palette.sampling_frac_asin.12$breaks,  # Breaks from custom_palette
      xlab = "Sampling Fraction (asin)",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Histogram of Sampling Fraction (asin)",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    sampling_frac_asin.12.combined <- sampling_frac_asin.12 / ((plot_spacer() | sampling_frac_asin.12.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Sampling Fraction (Arcsine Transformed)")
    
  }
  
}

#patchwork arangement
{
  mean_tip_rate.12.combined
  tip_rate_range_weighted_mean.12.combined
  mean_lineage_rate.12.combined
  mean_lineage_rate_range_weighted.12.combined
  sp.count.12.combined
  sampling_frac.12.combined
  sampling_frac_asin.12.combined
  
}

#all passerines baseline plots
{
  
  #year round ranges  
  {
    #species_key_count_range_weighted_log.passerines.1
    { 
      custom_palette.species_key_count_range_weighted_log.passerines.1 <- assignRateColors(
        rates = setNames(
          passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted_log,
          passeriformes.filtered.1.metrics$metrics_df$h3_cell
        ),
        breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
        logcolor = FALSE,           # Do not apply log transformation
        palette = "RdYlBu",         # Specify the color palette
        nbreaks = 20,             # Let assignRateColors determine the number of breaks
        reverse = TRUE,              # Reverse the palette
        largeN = 20000,
        samp_prop = 1
      )
      
      # Use the modified function with the custom palette
      species_key_count_range_weighted_log.passerines.1<-plot_h3_metrics_global.custom.raster(
        results = passeriformes.filtered.1.metrics,  # Your results object
        polygons_sf = polygons_passerines_sf.1,             # The polygons_sf object (precomputed H3 polygons)
        metric_column = "species_key_count_range_weighted_log",         # The column you want to visualize
        custom_palette = alpha(custom_palette.species_key_count_range_weighted_log.passerines.1$colors, 1),  # Use the colors from your custom_palette
        crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
        ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
        border_thickness = 0.25, lat_graticule_spacing = 15, 
        lon_graticule_spacing = 15, graticule_color = 'white',
        graticule_linetype = 1, 
        graticule_thickness = 0.1,
        raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
      )
      
      
      species_key_count_range_weighted_log.passerines.1.hist <-ratesHistogramGG(
        rates = setNames(
          passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted_log,
          passeriformes.filtered.1.metrics$metrics_df$h3_cell
        ),                               # Named vector of rates
        colors = custom_palette.species_key_count_range_weighted_log.passerines.1$colors,  # Colors from custom_palette
        breaks = custom_palette.species_key_count_range_weighted_log.passerines.1$breaks,  # Breaks from custom_palette
        xlab = "Log range-weighted\n species counts",          # X-axis label
        ylab = "Cell Counts",                 # Y-axis label
        title = "Histogram of Richness",  # Title
        useDensity = FALSE,              # Plot counts (not density)
        logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
      ) + ylim(0,1050) + xlim(-9, 10)
      
      richness.combined <- species_key_count_range_weighted_log.passerines.1 / ((plot_spacer() | species_key_count_range_weighted_log.passerines.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Local Species Richnes")
      
      
    }
    
    #species_key_count_range_weighted.passerines.1 
    {
      custom_palette.species_key_count_range_weighted.passerines.1 <- assignRateColors(
        rates = setNames(
          passeriformes.filtered.1.metrics$metrics_df$species_key_count_range_weighted,
          passeriformes.filtered.1.metrics$metrics_df$h3_cell
        ),
        breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
        logcolor = FALSE,           # Do not apply log transformation
        palette = "RdYlBu",         # Specify the color palette
        nbreaks = 20,             # Let assignRateColors determine the number of breaks
        reverse = TRUE,              # Reverse the palette
        largeN = 20000,
        samp_prop = 1
      )
      
      # Use the modified function with the custom palette
      species_key_count_range_weighted.passerines.1 <-plot_h3_metrics_global.custom.raster(
        results = passeriformes.filtered.1.metrics,  # Your results object
        polygons_sf = polygons_passerines_sf.1,             # The polygons_sf object (precomputed H3 polygons)
        metric_column = "species_key_count_range_weighted",         # The column you want to visualize
        custom_palette = alpha(custom_palette.species_key_count_range_weighted.passerines.1$colors, 1),  # Use the colors from your custom_palette
        crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
        ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
        border_thickness = 0.25, lat_graticule_spacing = 15, 
        lon_graticule_spacing = 15, graticule_color = 'white',
        graticule_linetype = 1, 
        graticule_thickness = 0.1,
        raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
      )
      
    }
    
    #species_key_count.passerines.1
    {
      custom_palette.species_key_count.passerines.1 <- assignRateColors(
        rates = setNames(
          passeriformes.filtered.1.metrics$metrics_df$species_key_count,
          passeriformes.filtered.1.metrics$metrics_df$h3_cell
        ),
        breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
        logcolor = FALSE,           # Do not apply log transformation
        palette = "RdYlBu",         # Specify the color palette
        nbreaks = 20,             # Let assignRateColors determine the number of breaks
        reverse = TRUE              # Reverse the palette
      )
      
      # Use the modified function with the custom palette
      species_key_count.passerines.1<-plot_h3_metrics_global.custom.raster(
        results = passeriformes.filtered.1.metrics,  # Your results object
        polygons_sf = polygons_passerines_sf.1,             # The polygons_sf object (precomputed H3 polygons)
        metric_column = "species_key_count",         # The column you want to visualize
        custom_palette = alpha(custom_palette.species_key_count.passerines.1$colors, 1),  # Use the colors from your custom_palette
        crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
        ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
        border_thickness = 0.25, lat_graticule_spacing = 15, 
        lon_graticule_spacing = 15, graticule_color = 'white',
        graticule_linetype = 1, 
        graticule_thickness = 0.1,
        raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
      )
      
      
      species_key_count.passerines.1.hist <-ratesHistogramGG(
        rates = setNames(
          passeriformes.filtered.1.metrics$metrics_df$species_key_count,
          passeriformes.filtered.1.metrics$metrics_df$h3_cell
        ),                               # Named vector of rates
        colors = custom_palette.species_key_count.passerines.1$colors,  # Colors from custom_palette
        breaks = custom_palette.species_key_count.passerines.1$breaks,  # Breaks from custom_palette
        xlab = "Species Counts",          # X-axis label
        ylab = "Cell Counts",                 # Y-axis label
        title = "Histogram of Species Counts",  # Title
        useDensity = FALSE,              # Plot counts (not density)
        logscale = FALSE, lwd=0.1, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
      ) + ylim(0, 1300) + xlim(0, 320)
      
      counts.combined <- species_key_count.passerines.1 / ((plot_spacer() | species_key_count.passerines.1.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of Species Counts")
      
      
    }
    
  }
  
  #year round ranges + breeding ranges
  {
    #species_key_count_range_weighted_log.passerines.12
    { 
      custom_palette.species_key_count_range_weighted_log.passerines.12 <- assignRateColors(
        rates = setNames(
          passeriformes.filtered.12.metrics$metrics_df$species_key_count_range_weighted_log,
          passeriformes.filtered.12.metrics$metrics_df$h3_cell
        ),
        breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
        logcolor = FALSE,           # Do not apply log transformation
        palette = "RdYlBu",         # Specify the color palette
        nbreaks = 20,             # Let assignRateColors determine the number of breaks
        reverse = TRUE,              # Reverse the palette
        largeN = 20000,
        samp_prop = 1
      )
      
      # Use the modified function with the custom palette
      species_key_count_range_weighted_log.passerines.12<-plot_h3_metrics_global.custom.raster(
        results = passeriformes.filtered.12.metrics,  # Your results object
        polygons_sf = polygons_passerines_sf.12,             # The polygons_sf object (precomputed H3 polygons)
        metric_column = "species_key_count_range_weighted_log",         # The column you want to visualize
        custom_palette = alpha(custom_palette.species_key_count_range_weighted_log.passerines.12$colors, 1),  # Use the colors from your custom_palette
        crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
        ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
        border_thickness = 0.25, lat_graticule_spacing = 20, 
        lon_graticule_spacing = 20, graticule_color = 'white',
        graticule_linetype = 1, 
        graticule_thickness = 0.1,
        raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
      )
    }
    
    #species_key_count_range_weighted.passerines.12
    {
      custom_palette.species_key_count_range_weighted.passerines.12 <- assignRateColors(
        rates = setNames(
          passeriformes.filtered.12.metrics$metrics_df$species_key_count_range_weighted,
          passeriformes.filtered.12.metrics$metrics_df$h3_cell
        ),
        breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
        logcolor = FALSE,           # Do not apply log transformation
        palette = "RdYlBu",         # Specify the color palette
        nbreaks = 20,             # Let assignRateColors determine the number of breaks
        reverse = TRUE,              # Reverse the palette
        largeN = 20000,
        samp_prop = 1
      )
      
      # Use the modified function with the custom palette
      species_key_count_range_weighted.passerines.12 <-plot_h3_metrics_global.custom.raster(
        results = passeriformes.filtered.12.metrics,  # Your results object
        polygons_sf = polygons_passerines_sf.12,             # The polygons_sf object (precomputed H3 polygons)
        metric_column = "species_key_count_range_weighted",         # The column you want to visualize
        custom_palette = alpha(custom_palette.species_key_count_range_weighted.passerines.12$colors, 1),  # Use the colors from your custom_palette
        crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
        ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
        border_thickness = 0.25, lat_graticule_spacing = 20, 
        lon_graticule_spacing = 20, graticule_color = 'white',
        graticule_linetype = 1, 
        graticule_thickness = 0.1,
        raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
      )
      
    }
    
    #species_key_count.passerines.12
    {
      custom_palette.species_key_count.passerines.12 <- assignRateColors(
        rates = setNames(
          passeriformes.filtered.12.metrics$metrics_df$species_key_count,
          passeriformes.filtered.12.metrics$metrics_df$h3_cell
        ),
        breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
        logcolor = FALSE,           # Do not apply log transformation
        palette = "RdYlBu",         # Specify the color palette
        nbreaks = 20,             # Let assignRateColors determine the number of breaks
        reverse = TRUE,              # Reverse the palette
        largeN = 20000,
        samp_prop = 1
      )
      
      # Use the modified function with the custom palette
      species_key_count.passerines.12<-plot_h3_metrics_global.custom.raster(
        results = passeriformes.filtered.12.metrics,  # Your results object
        polygons_sf = polygons_passerines_sf.12,             # The polygons_sf object (precomputed H3 polygons)
        metric_column = "species_key_count",         # The column you want to visualize
        custom_palette = alpha(custom_palette.species_key_count.passerines.12$colors, 1),  # Use the colors from your custom_palette
        crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",                       # Coordinate Reference System
        ocean_color = alpha('#87B7DB', 0.5),               # Color for the ocean
        border_thickness = 0.25, lat_graticule_spacing = 20, 
        lon_graticule_spacing = 20, graticule_color = 'white',
        graticule_linetype = 1, 
        graticule_thickness = 0.1,
        raster_basemap = NULL, raster_alpha = 1 # Thickness of the ocean border
      )
      
    }
  }
  
  species_key_count.passerines.1
  species_key_count_range_weighted.passerines.1 /
    species_key_count_range_weighted_log.passerines.1
  
  species_key_count.passerines.12
  species_key_count_range_weighted.passerines.12 /
    species_key_count_range_weighted_log.passerines.12
  
  
  pdf(file = 'passerine_diversity.combined.pdf', height=10, width=7.5*2)
  (counts.combined | richness.combined) +
    plot_annotation(tag_levels = 'A') &
    theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  dev.off()
  
}

#starting to set up spatial lag models
{

spatial_coords <- h3_to_geo(resident_species.1.h3.metrics$metrics_df$h3_cell)
spatial_coords <- data.frame(spatial_coords, 
                             skelevision_counts = resident_species.1.h3.metrics$metrics_df$species_key_count, 
                             h3 = resident_species.1.h3.metrics$metrics_df$h3_cell,  
                             mean_lineage_rate_range_weighted = resident_species.1.h3.metrics$metrics_df$mean_lineage_rate_range_weighted,
                             mean_lineage_rate = resident_species.1.h3.metrics$metrics_df$mean_lineage_rate,
                             mean_tip_rate = resident_species.1.h3.metrics$metrics_df$mean_tip_rate,
                             tip_rate_range_weighted_mean = resident_species.1.h3.metrics$metrics_df$tip_rate_range_weighted_mean,
                             mean_tipDR = resident_species.1.h3.metrics$metrics_df$mean_tipDR,
                             tipDR_range_weighted_mean = resident_species.1.h3.metrics$metrics_df$tipDR_range_weighted_mean,
                             sampling_frac = resident_species.1.h3.metrics$metrics_df$sampling_frac,
                             passerines_counts = resident_species.1.h3.metrics$metrics_df$species_key_count.passerines,
                             passerines_counts_range_weighted = resident_species.1.h3.metrics$metrics_df$species_key_count_range_weighted.passerines
)

spatial_coords.12 <- h3_to_geo(resident_species.12.h3.metrics$metrics_df$h3_cell)
spatial_coords.12 <- data.frame(spatial_coords.12, 
                                skelevision_counts = resident_species.12.h3.metrics$metrics_df$species_key_count, 
                                h3 = resident_species.12.h3.metrics$metrics_df$h3_cell,  
                                mean_lineage_rate_range_weighted = resident_species.12.h3.metrics$metrics_df$mean_lineage_rate_range_weighted,
                                mean_lineage_rate = resident_species.12.h3.metrics$metrics_df$mean_lineage_rate,
                                mean_tip_rate = resident_species.12.h3.metrics$metrics_df$mean_tip_rate,
                                tip_rate_range_weighted_mean = resident_species.12.h3.metrics$metrics_df$tip_rate_range_weighted_mean,
                                mean_tipDR = resident_species.12.h3.metrics$metrics_df$mean_tipDR,
                                tipDR_range_weighted_mean = resident_species.12.h3.metrics$metrics_df$tipDR_range_weighted_mean,
                                sampling_frac = resident_species.12.h3.metrics$metrics_df$sampling_frac,
                                passerines_counts = resident_species.12.h3.metrics$metrics_df$species_key_count.passerines,
                                passerines_counts_range_weighted = resident_species.12.h3.metrics$metrics_df$species_key_count_range_weighted.passerines
)

  
  
}

#setting up bioclim
{
  
  # # # # Run the function
  # spatial_coords <- load_all_bioclim(
  #   bioclim_dir = "/Users/cotinga/Downloads/wc2.1_5m_bio",
  #   spatial_coords = spatial_coords,
  #   polygons_sf = polygons_sf.1,
  #   summary_function_name = 'mean'
  # )
  # saveRDS(spatial_coords, file='spatial_coords.RDS')
  spatial_coords<- readRDS('spatial_coords.RDS')
  
  # # # # Run the function
  # spatial_coords.12 <- load_all_bioclim(
  #   bioclim_dir = "/Users/cotinga/Downloads/wc2.1_5m_bio",
  #   spatial_coords = spatial_coords.12,
  #   polygons_sf = polygons_sf.12,
  #   summary_function_name = 'mean'
  # )
  # saveRDS(spatial_coords.12, file='spatial_coords.12.RDS')
  spatial_coords.12<- readRDS('spatial_coords.12.RDS')
  
  length(spatial_coords$lat)
  length(spatial_coords.12$lat)
  
}

#merged datasets with bioclim PCA (testing)
{
  # Step 1: Remove rows with missing values
  climate_vars.cleaned <- spatial_coords[complete.cases(spatial_coords), ]
  #generating PC scores for bioclim -- only including those variables linked to patterns of variability
  pca.clim <- prcomp(climate_vars.cleaned[, c("wc2.1_5m_bio_2_mean",
                                              "wc2.1_5m_bio_4_mean",
                                              "wc2.1_5m_bio_7_mean",
                                              "wc2.1_5m_bio_15_mean")], scale. = TRUE)
  
  #plot(pca.clim$x[,1] ~ climate_vars.cleaned$wc2.1_5m_bio_4_mean)
  climate_vars.cleaned<-cbind(climate_vars.cleaned, pca.clim$x)
  climate_vars.cleaned <- climate_vars.cleaned[!rownames(climate_vars.cleaned) %in% c("8113", "9003"), ] #remove 2 cells because they have no neighbors at the 1000 radius
  rownames(climate_vars.cleaned)<-NULL
  
  # Step 1: Remove rows with missing values
  climate_vars.cleaned.12 <- spatial_coords.12[complete.cases(spatial_coords.12), ]
  #generating PC scores for bioclim
  pca.clim.12 <- prcomp(climate_vars.cleaned.12[, c("wc2.1_5m_bio_2_mean",
                                                    "wc2.1_5m_bio_4_mean",
                                                    "wc2.1_5m_bio_7_mean",
                                                    "wc2.1_5m_bio_15_mean")], scale. = TRUE)
  #plot(pca.clim$x[,1] ~ climate_vars.cleaned$wc2.1_5m_bio_4_mean)
  climate_vars.cleaned.12<-cbind(climate_vars.cleaned.12, pca.clim.12$x)
  
}

#Compute spatial weights
{
  #1000km
  {
    #Compute spatial weights
    weights_idw.1000.clim <- nb2listw(
      dnearneigh(as.matrix(climate_vars.cleaned[, c("lng", "lat")]), d1 = 0, d2 = 1000, longlat = TRUE),
      glist = lapply(
        nbdists(
          dnearneigh(as.matrix(climate_vars.cleaned[, c("lng", "lat")]), d1 = 0, d2 = 1000, longlat = TRUE),
          as.matrix(climate_vars.cleaned[, c("lng", "lat")]),
          longlat = TRUE
        ),
        function(x) 1 / x^2
      ),
      style = "W",
      zero.policy = TRUE
    )
    
    #testing more efficient computaton -- it is identical 
    {
      
      # Compute spatial weights (linear inverse distance weighting)
      neighbors <- dnearneigh(as.matrix(climate_vars.cleaned[, c("lng", "lat")]), d1 = 0, d2 = 1000, longlat = TRUE)
      distances <- nbdists(neighbors, as.matrix(climate_vars.cleaned[, c("lng", "lat")]), longlat = TRUE)
      
      weights_idw.1000.clim.tmp <- nb2listw(
        neighbors,
        glist = lapply(distances, function(x) 1 / x ^2),
        style = "W",
        zero.policy = TRUE
      )
      
      {
    
        
        # Extract longitude and latitude
        coords <- as.matrix(climate_vars.cleaned[, c("lng", "lat")])
        
        # Compute neighbors within 1000 km
        nb <- dnearneigh(coords, d1 = 0, d2 = 1000, longlat = TRUE)
        
        # Compute inverse distance squared weights (avoiding division by zero)
        dist_list <- nbdists(nb, coords, longlat = TRUE)
        dist_list <- lapply(dist_list, function(x) 1 / x^2)  # Avoid division by zero
        
        # Convert to spatial weights matrix
        weights_idw.1000.clim <- nb2listw(nb, glist = dist_list, style = "W", zero.policy = TRUE)
        
      }
      {
   
        
        # Extract longitude and latitude
        coords <- as.matrix(climate_vars.cleaned[, c("lng", "lat")])
        
        # Compute neighbors within 1000 km
        nb <- dnearneigh(coords, d1 = 0, d2 = 1000, longlat = TRUE)
        
        # Compute inverse distance squared weights (avoiding division by zero)
        dist_list <- nbdists(nb, coords, longlat = TRUE)
        dist_list <- lapply(dist_list, function(x) 1 / x)  # Avoid division by zero
        
        # Convert to spatial weights matrix
        weights_idwl.1000.clim <- nb2listw(nb, glist = dist_list, style = "W", zero.policy = TRUE)
        
      }
      
      
    }
    
    
    #Compute spatial weights (linear)
    weights_idwl.1000.clim <- nb2listw(
      dnearneigh(as.matrix(climate_vars.cleaned[, c("lng", "lat")]), d1 = 0, d2 = 1000, longlat = TRUE),
      glist = lapply(
        nbdists(
          dnearneigh(as.matrix(climate_vars.cleaned[, c("lng", "lat")]), d1 = 0, d2 = 1000, longlat = TRUE),
          as.matrix(climate_vars.cleaned[, c("lng", "lat")]),
          longlat = TRUE
        ),
        function(x) 1 / x
      ),
      style = "W",
      zero.policy = TRUE
    )
    
    
    
  }
  
}

#standardizing vars
{
  #now trying with standardized coefficients
  climate_vars.cleaned$wc2.1_5m_bio_4_mean_z <- scale(climate_vars.cleaned$wc2.1_5m_bio_4_mean)
  climate_vars.cleaned$lat_z <- scale(climate_vars.cleaned$lat, center=F) #no centering because 0 mean equator
  climate_vars.cleaned$tipDR_range_weighted_mean_z <- scale(climate_vars.cleaned$tipDR_range_weighted_mean)
  climate_vars.cleaned$passerines_counts_range_weighted_z <- scale(log(climate_vars.cleaned$passerines_counts_range_weighted))
  climate_vars.cleaned$passerines_counts_z <- scale(log(climate_vars.cleaned$passerines_counts))
  
}

#fitting SAR-Lag models (aggregate) and leave one out analyses
{
  slm_model.step9c.2c <- lagsarlm(
    mean_lineage_rate_range_weighted ~ 
      (poly(wc2.1_5m_bio_4_mean, 2) + poly(log(passerines_counts_range_weighted), 2, raw=T)) * poly(lat_z, 2, raw=T) + 
      asin(sqrt(sampling_frac)) + tipDR_range_weighted_mean_z,
    data = climate_vars.cleaned, 
    listw = weights_idw.1000.clim, 
    zero.policy = TRUE,
    method = "eigen"
  )
  
  slm_model.step9c.2c.noDR <- lagsarlm(
    mean_lineage_rate_range_weighted ~ 
      (poly(wc2.1_5m_bio_4_mean, 2) + poly(log(passerines_counts_range_weighted), 2, raw=T)) * poly(lat_z, 2, raw=T) + 
      asin(sqrt(sampling_frac)),
    data = climate_vars.cleaned, 
    listw = weights_idw.1000.clim, 
    zero.policy = TRUE,
    method = "eigen"
  )
  
  #leave one out model weights for sar-lag aggregated
  step9c.2.loo <- loo_model_weights(slm_model.step9c.2c, climate_vars.cleaned, weights_idw.1000.clim, recompute_full_model = T)
  #step9c.2.loo <- loo_model_weights_from_obj(slm_model.step9c.2c, weights_idw.1000.clim)
  #alternative using the model object based approach
  
  #generate stargazer latex table of loo
  {

    
    # Create a cleaned-up copy
    step9c.2.loo.metrics_clean <- step9c.2.loo$metrics %>%
      rename(
        "Term" = Term,
        "AIC" = AIC,
        "BIC" = BIC,
        "ΔAIC" = DeltaAIC,
        "ΔBIC" = DeltaBIC,
        "AIC Weight" = AIC_Weight,
        "BIC Weight" = BIC_Weight
      ) %>%
      mutate(
        Term = str_replace_all(Term, 
                               "poly\\(log\\(passerines_counts_range_weighted\\), 2, raw = T\\)1", 
                               "Species richness"),
        Term = str_replace_all(Term, 
                               "poly\\(log\\(passerines_counts_range_weighted\\), 2, raw = T\\)2", 
                               "Species richness²"),
        Term = str_replace_all(Term, 
                               "poly\\(wc2.1_5m_bio_4_mean, 2\\)1", 
                               "Temperature seasonality"),
        Term = str_replace_all(Term, 
                               "poly\\(wc2.1_5m_bio_4_mean, 2\\)2", 
                               "Temperature seasonality²"),
        Term = str_replace_all(Term, 
                               "poly\\(lat_z, 2, raw = T\\)1", 
                               "Latitude"),
        Term = str_replace_all(Term, 
                               "poly\\(lat_z, 2, raw = T\\)2", 
                               "Latitude²"),
        Term = str_replace_all(Term, 
                               "Full Model", 
                               "Full model"),
        Term = str_replace_all(Term, 
                               ":", 
                               " × ")
      )
    
    # Print cleaned LaTeX table
    stargazer(
      step9c.2.loo.metrics_clean,
      summary = FALSE,
      rownames = FALSE,
      title = "LOO Model Weights",
      label = "tab:model_metrics_clean",
      digits = 2,
      font.size = "small"
    )
    
  }
  
  #plots of predicted values under aggregate models
  {
    
    predicted_rates_df <- list(
      metrics_df = data.frame(
        h3_cell = climate_vars.cleaned$h3,  # H3 cell IDs
        #predicted_rates.2 = fitted(slm_model.step2),
        #predicted_rates.7b = fitted(slm_model.step7b),
        #predicted_rates.9 = fitted(slm_model.step9),
        predicted_rates.9c.2c = fitted(slm_model.step9c.2c),
        predicted_rates.9c.2c.residuals = residuals(slm_model.step9c.2c)#,
        #predicted_rates.9c.2c.durbin = fitted(slm_model.step9c.2c.durbin),
        #predicted_rates.9c.2c.durbin.residuals = residuals(slm_model.step9c.2c.durbin)
      )
    )
    
    
    #9c.2c
    {
      #predicted
      {
        custom_palette.predicted_rates.9c.2c <- assignRateColors(
          rates = setNames(
            predicted_rates_df$metrics_df$predicted_rates.9c.2c,
            predicted_rates_df$metrics_df$h3_cell
          ),
          breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
          logcolor = FALSE,           # Do not apply log transformation
          palette = "RdYlBu",         # Specify the color palette
          nbreaks = 30,             # Let assignRateColors determine the number of breaks
          reverse = TRUE,              # Reverse the palette
          largeN = 20000,
          samp_prop = 1
        )
        
        # Step 3: Visualize the residuals using the custom plotting function
        predicted_rates_plot.9c.2c <- plot_h3_metrics_global.custom.raster(
          results = predicted_rates_df,                          # The constructed results object
          polygons_sf = polygons_sf.1,               # The polygons with H3 indices
          metric_column = "predicted_rates.9c.2c",               # The column for residuals
          custom_palette = custom_palette.predicted_rates.9c.2c$colors,    # Use a viridis palette for residuals
          crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",  # Coordinate Reference System
          ocean_color = alpha('#87B7DB', 0.5),       # Color for the ocean
          border_thickness = 0.25, 
          lat_graticule_spacing = 15, 
          lon_graticule_spacing = 15, 
          graticule_color = 'white',
          graticule_linetype = 1, 
          graticule_thickness = 0.1,
          raster_basemap = NULL, 
          raster_alpha = 1                            # Thickness of the ocean border
        )
        
        predicted_rates.hist.9c.2c <-ratesHistogramGG(
          rates = setNames(
            predicted_rates_df$metrics_df$predicted_rates.9c.2c,
            predicted_rates_df$metrics_df$h3_cell
          ),                               # Named vector of rates
          colors = custom_palette.predicted_rates.9c.2c$colors,  # Colors from custom_palette
          breaks = custom_palette.predicted_rates.9c.2c$breaks,  # Breaks from custom_palette
          xlab = "Predicted Rates",          # X-axis label
          ylab = "Counts",                 # Y-axis label
          title = "Predicted Mean Rates",  # Title
          useDensity = FALSE,              # Plot counts (not density)
          logscale = FALSE, lwd=0.01, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
        )
      }
      
      valid_methods <- c("linear", "quantile", "jenks", "sd", "equal", 
                         "pretty", "kmeans", "hclust", "bclust", "fisher", 
                         "dpih", "headtails", "maximum", "box")
      
      
      #residuals
      {
        custom_palette.residuals.9c.2c <- assignRateColors(
          rates = setNames(
            predicted_rates_df$metrics_df$predicted_rates.9c.2c.residuals,
            predicted_rates_df$metrics_df$h3_cell
          ),
          breaksmethod = "linear",    # Use "fisher" for better performance with large datasets
          logcolor = FALSE,           # Do not apply log transformation
          palette = "RdYlBu",         # Specify the color palette
          nbreaks = 30,             # Let assignRateColors determine the number of breaks
          reverse = TRUE,              # Reverse the palette
          largeN = 20000,
          samp_prop = 1
        )
        
        # Step 3: Visualize the residuals using the custom plotting function
        residual_rates_plot.9c.2c <- plot_h3_metrics_global.custom.raster(
          results = predicted_rates_df,                          # The constructed results object
          polygons_sf = polygons_sf.1,               # The polygons with H3 indices
          metric_column = "predicted_rates.9c.2c.residuals",               # The column for residuals
          custom_palette = custom_palette.residuals.9c.2c$colors,    # Use a viridis palette for residuals
          crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",  # Coordinate Reference System
          ocean_color = alpha('#87B7DB', 0.5),       # Color for the ocean
          border_thickness = 0.25, 
          lat_graticule_spacing = 15, 
          lon_graticule_spacing = 15, 
          graticule_color = 'white',
          graticule_linetype = 1, 
          graticule_thickness = 0.1,
          raster_basemap = NULL, 
          raster_alpha = 1                            # Thickness of the ocean border
        )
        
        residuals.hist.9c.2c.residuals <-ratesHistogramGG(
          rates = setNames(
            predicted_rates_df$metrics_df$predicted_rates.9c.2c.residuals,
            predicted_rates_df$metrics_df$h3_cell
          ),                               # Named vector of rates
          colors = custom_palette.residuals.9c.2c$colors,  # Colors from custom_palette
          breaks = custom_palette.residuals.9c.2c$breaks,  # Breaks from custom_palette
          xlab = "Residuals",          # X-axis label
          ylab = "Counts (pseudo-log)",                 # Y-axis label
          title = "Residual Variability",  # Title
          useDensity = FALSE,              # Plot counts (not density)
          logscale = FALSE, lwd=0.01, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
        )+
          scale_y_continuous(
            trans  = pseudo_log_trans(base = 10, sigma = 20)
          ) +
          coord_cartesian(
            ylim   = c(-5, NA),   # force y to start at 0, let ggplot pick the top
            expand = FALSE
          ) 
        
        #+ xlim(-1.5,1.5)
        #generate_color_bar(break_colors = custom_palette.residuals.9c.2c$break_colors)
        
        
      }
      
      pdf(file="predicted_vs_residuals.pdf", height=9*2, width=4.5*2)
      predicted_rates_plot.9c.2c /
        residual_rates_plot.9c.2c
      dev.off()
      
    }

    
  }
  
  
  mean_lineage_rate_range_weighted.1 /
    predicted_rates_plot.9c.2c /
    residual_rates_plot.9c.2c
  
  predicted_rates_plot.9c.2c.durbin / 
    predicted_rates_plot.9c.2c
  
  plot(fitted(slm_model.step9c.2c) ~ slm_model.step9c.2c$y)
  summary(lm(fitted(slm_model.step9c.2c) ~ slm_model.step9c.2c$y))
  
  plot(fitted(slm_model.step9c.2c.durbin) ~ slm_model.step9c.2c$y)
  summary(lm(fitted(slm_model.step9c.2c.durbin) ~ slm_model.step9c.2c$y))
  
  
  
}

#fitting independent effects SAR-lag
{
  #richness (idw)
  {
    #linear fit - without and with controls
    {
      # slm_model.richness.1000idw.lin <- lagsarlm(
      #   mean_lineage_rate_range_weighted ~ 
      #     passerines_counts_range_weighted_z,
      #   data = climate_vars.cleaned, 
      #   listw = weights_idw.1000.clim, 
      #   zero.policy = F,
      #   method = "eigen"
      # )
      # 
      # summary(slm_model.richness.1000idw.lin, Nagelkerke=T)
      # AIC(slm_model.richness.1000idw.lin) #-12885.73
      # BIC(slm_model.richness.1000idw.lin) #-12856.71
      # moran.test(residuals(slm_model.richness.1000idw.lin), weights_idw.1000.clim)
      # 0.06889648
      
      
      slm_model.richness.1000idw.lin.con <- lagsarlm(
        mean_lineage_rate_range_weighted ~ 
          passerines_counts_range_weighted_z * I(lat^2) + asin(sqrt(sampling_frac)),
        data = climate_vars.cleaned, 
        listw = weights_idw.1000.clim, 
        zero.policy = F,
        method = "eigen"
      )
      
      summary(slm_model.richness.1000idw.lin.con, Nagelkerke=T)
      AIC(slm_model.richness.1000idw.lin.con) #-12883.01
      BIC(slm_model.richness.1000idw.lin.con) #-12832.22
      moran.test(residuals(slm_model.richness.1000idw.lin.con), weights_idw.1000.clim)
      0.06942038
      
      slm_model.richness.1000idw.lin.con.loo <- loo_model_weights(model = slm_model.richness.1000idw.lin.con, data = climate_vars.cleaned, listw = weights_idw.1000.clim, recompute_full_model = T )
      
      
      # slm_model.richness.1000idw.lin.con.poly <- lagsarlm(
      #   mean_lineage_rate_range_weighted ~ 
      #     passerines_counts_range_weighted_z * poly(lat, 2, raw=T) + asin(sqrt(sampling_frac)),
      #   data = climate_vars.cleaned, 
      #   listw = weights_idw.1000.clim, 
      #   zero.policy = F,
      #   method = "eigen"
      # )
      
      
      #non-spatial model
      lm_model.richness.1000idw.lin.con <- lm(
        mean_lineage_rate_range_weighted ~ 
          passerines_counts_range_weighted_z * I(lat^2) + asin(sqrt(sampling_frac)),
        data = climate_vars.cleaned)
      moran.test(residuals(lm_model.richness.1000idw.lin.con), weights_idw.1000.clim)
      
      #decrease in spatial autocorrelation
      moran.richness<-((0.5016752 - 0.06942038) / 0.5016752) * 100
      
      
    }
    
  }
  
  #bio4 (idw)
  {

    slm_model.bio4_z.1000idw.con <- lagsarlm(
      mean_lineage_rate_range_weighted ~ 
        poly(wc2.1_5m_bio_4_mean_z, 2, raw=T) * poly(lat_z, 2, raw=T) + asin(sqrt(sampling_frac)),
      data = climate_vars.cleaned, 
      listw = weights_idw.1000.clim, 
      zero.policy = F,
      method = "eigen"
    )
    summary(slm_model.bio4_z.1000idw.con)
    AIC(slm_model.bio4_z.1000idw.con) # -12884.73
    BIC(slm_model.bio4_z.1000idw.con) #-12797.66
    moran.test(residuals(slm_model.bio4_z.1000idw.con), weights_idw.1000.clim)
    6.902988e-02
    
    slm_model.bio4_z.1000idw.con.loo <- loo_model_weights(model = slm_model.bio4_z.1000idw.con, data = climate_vars.cleaned, listw = weights_idw.1000.clim, recompute_full_model = T )
    
    
    #non-spatial model
    lm_model.bio4.2.1000idw <- lm(
      mean_lineage_rate_range_weighted ~ 
        poly(wc2.1_5m_bio_4_mean, 2) * poly(lat, 2) + asin(sqrt(sampling_frac)),
      data = climate_vars.cleaned)
    moran.test(residuals(lm_model.bio4.2.1000idw), weights_idw.1000.clim)
    
    
    #decrease in spatial autocorrelation
    moran.bio4<-((0.4316325 - 0.06902988) / 0.4316325) * 100
    #84.00726
    
    
  }
  
  #lat2 (idw)
  {
    
    slm_model.lat2.1000idw.con <- lagsarlm(
      mean_lineage_rate_range_weighted ~ 
        I(lat^2) + asin(sqrt(sampling_frac)),  # Preserve quadratic interpretation
      data = climate_vars.cleaned, 
      listw = weights_idw.1000.clim, 
      zero.policy = TRUE,
      method = "eigen"
    )
    summary(slm_model.lat2.1000idw.con, Nagelkerke=T) #0.74171 
    AIC(slm_model.lat2.1000idw.con)  #-12882.1
    BIC(slm_model.lat2.1000idw.con) #-12845.82
    moran.test(residuals(slm_model.lat2.1000idw.con), weights_idw.1000.clim)
    0.06936483
    
    slm_model.lat2.1000idw.con.loo <- loo_model_weights(model = slm_model.lat2.1000idw.con, data = climate_vars.cleaned, listw = weights_idw.1000.clim, recompute_full_model = T )
    
    
    #non-spatial model
    lm_model.lat2.1000idw <- lm(
      mean_lineage_rate_range_weighted ~ 
        I(lat^2) + asin(sqrt(sampling_frac)),  # Preserve quadratic interpretation
      data = climate_vars.cleaned)
    moran.test(residuals(lm_model.lat2.1000idw), weights_idw.1000.clim)
    moran.lat<-((0.5014061 - 0.06936483) / 0.5021143) * 100
    
    
  }
  
}

#computing spatial impacts
{
  # Compute traces of the weights matrix
  W <- as(weights_idw.1000.clim, "CsparseMatrix")
  #trMat <- trW(W, type = "mult")  # Compute trace-based approximation
  #saveRDS(trMat, file='trMat.1000.idw.RDS')
  trMat<-readRDS(trMat, file='trMat.1000.idw.RDS')
  
  #independent models
  #latitude + samling fraction
  impacts_approx.slm_model.lat2.1000idw.con <- impacts(slm_model.lat2.1000idw.con, tr = trMat, R=1000)
  impacts_approx.summary.slm_model.lat2.1000idw.con <- summary(impacts_approx.slm_model.lat2.1000idw.con, zstats=T, short=F)
  
  #bio4 + lat + interaction + sampling fraction
  impacts_approx.slm_model.bio4_z.1000idw.con <- impacts(slm_model.bio4_z.1000idw.con, tr = trMat, R=1000)
  impacts_approx.summary.slm_model.bio4_z.1000idw.con <- summary(impacts_approx.slm_model.bio4_z.1000idw.con, zstats=T, short=F)
  
  #richness + lat2 + sampling fraction + interaction between lat and counts
  impacts_approx.slm_model.richness.1000idw.lin.con <- impacts(slm_model.richness.1000idw.lin.con, tr = trMat, R=1000)
  impacts_approx.summary.slm_model.richness.1000idw.lin.con <- summary(impacts_approx.slm_model.richness.1000idw.lin.con, zstats=T, short=F)
  
  #impacts scores for aggregate models
  impacts_approx.step9c.2c <- impacts(slm_model.step9c.2c, tr = trMat, R=1000)
  impacts_approx.summary.step9c.2c <- summary(impacts_approx.step9c.2c, zstats=T, short=F)
  summary(slm_model.step9c.2c, Nagelkerke=T)
  
  impacts_approx.step9c.2c.noDR <- impacts(slm_model.step9c.2c.noDR, tr = trMat, R=1000)
  impacts_approx.summary.step9c.2c.noDR <- summary(impacts_approx.step9c.2c.noDR, zstats=T, short=F)
  summary(slm_model.step9c.2c, Nagelkerke=T)
  
 
  
  
}

#generate and add palettes to data
{latitude_data <- list(
  metrics_df = data.frame(
    h3_cell = climate_vars.cleaned$h3,  # H3 cell IDs
    lat = climate_vars.cleaned$lat,
    bio4 = climate_vars.cleaned$wc2.1_5m_bio_4_mean
  )
)
  
  custom_palette.lat <- assignRateColors(
    rates = setNames(
      latitude_data$metrics_df$lat,
      latitude_data$metrics_df$h3_cell
    ),
    breaksmethod = "linear",    # Use "fisher" for better performance with large datasets
    logcolor = FALSE,           # Do not apply log transformation
    palette = "RdYlBu",         # Specify the color palette
    nbreaks = 30,             # Let assignRateColors determine the number of breaks
    reverse = TRUE,              # Reverse the palette
    largeN = 20000,
    samp_prop = 1
  )
  
  custom_palette.abs.lat <- assignRateColors(
    rates = setNames(
      abs(latitude_data$metrics_df$lat),
      latitude_data$metrics_df$h3_cell
    ),
    breaksmethod = "linear",    # Use "fisher" for better performance with large datasets
    logcolor = FALSE,           # Do not apply log transformation
    palette = "RdYlBu",         # Specify the color palette
    nbreaks = 30,             # Let assignRateColors determine the number of breaks
    reverse = TRUE,              # Reverse the palette
    largeN = 20000,
    samp_prop = 1
  )
  
  custom_palette.lat2 <- assignRateColors(
    rates = setNames(
      (latitude_data$metrics_df$lat^2),
      latitude_data$metrics_df$h3_cell
    ),
    breaksmethod = "linear",    # Use "fisher" for better performance with large datasets
    logcolor = FALSE,           # Do not apply log transformation
    palette = "RdYlBu",         # Specify the color palette
    nbreaks = 30,             # Let assignRateColors determine the number of breaks
    reverse = TRUE,              # Reverse the palette
    largeN = 20000,
    samp_prop = 1
  )
  
  
  custom_palette.bio4 <- assignRateColors(
    rates = setNames(
      latitude_data$metrics_df$bio4,
      latitude_data$metrics_df$h3_cell
    ),
    breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
    logcolor = FALSE,           # Do not apply log transformation
    palette = "RdYlBu",         # Specify the color palette
    nbreaks = 30,             # Let assignRateColors determine the number of breaks
    reverse = TRUE,              # Reverse the palette
    largeN = 20000,
    samp_prop = 1
  )
  
  
  # Add a new column for colors based on the 'h3' mapping
  climate_vars.cleaned$custom_palette.bio4 <- custom_palette.bio4$colors[climate_vars.cleaned$h3]
  climate_vars.cleaned$custom_palette.lat <- custom_palette.lat$colors[climate_vars.cleaned$h3]
  climate_vars.cleaned$custom_palette.abs.lat <- custom_palette.abs.lat$colors[climate_vars.cleaned$h3]
  climate_vars.cleaned$custom_palette.lat2 <- custom_palette.lat2$colors[climate_vars.cleaned$h3]
}

#visualizing bioclim datasets
{
  
  bioclim_mean <- list(
    metrics_df = data.frame(
      h3_cell = climate_vars.cleaned$h3,  # H3 cell IDs
      #bio1 = climate_vars.cleaned$wc2.1_5m_bio_1_mean,
      bio2 = climate_vars.cleaned$wc2.1_5m_bio_2_mean,
      #bio3 = climate_vars.cleaned$wc2.1_5m_bio_3_mean,
      bio4 = climate_vars.cleaned$wc2.1_5m_bio_4_mean,
      #bio5 = climate_vars.cleaned$wc2.1_5m_bio_5_mean,
      #bio6 = climate_vars.cleaned$wc2.1_5m_bio_6_mean,
      bio7 = climate_vars.cleaned$wc2.1_5m_bio_7_mean,
      #bio8 = climate_vars.cleaned$wc2.1_5m_bio_8_mean,
      #bio9 = climate_vars.cleaned$wc2.1_5m_bio_9_mean,
      #bio10 = climate_vars.cleaned$wc2.1_5m_bio_10_mean,
      #bio11 = climate_vars.cleaned$wc2.1_5m_bio_11_mean,
      #bio12 = climate_vars.cleaned$wc2.1_5m_bio_12_mean,
      #bio13 = climate_vars.cleaned$wc2.1_5m_bio_13_mean,
      #bio14 = climate_vars.cleaned$wc2.1_5m_bio_14_mean,
      bio15 = climate_vars.cleaned$wc2.1_5m_bio_15_mean,
      #bio16 = climate_vars.cleaned$wc2.1_5m_bio_16_mean,
      #bio17 = climate_vars.cleaned$wc2.1_5m_bio_17_mean,
      #bio18 = climate_vars.cleaned$wc2.1_5m_bio_18_mean,
      #bio19 = climate_vars.cleaned$wc2.1_5m_bio_19_mean,
      PC1 = climate_vars.cleaned$PC1,
      PC2 = climate_vars.cleaned$PC2,
      PC3 = climate_vars.cleaned$PC3
    )
  )
  
  
  #Visualizing bioclim 4
  {
    custom_palette.bio4 <- assignRateColors(
      rates = setNames(
        bioclim_mean$metrics_df$bio4,
        bioclim_mean$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 30,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    # Step 3: Visualize the residuals using the custom plotting function
    bio4_plot <- plot_h3_metrics_global.custom.raster(
      results = bioclim_mean,                          # The constructed results object
      polygons_sf = polygons_sf.1,               # The polygons with H3 indices
      metric_column = "bio4",               # The column for residuals
      custom_palette = custom_palette.bio4$colors,    # Use a viridis palette for residuals
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",  # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),       # Color for the ocean
      border_thickness = 0.25, 
      lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, 
      graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, 
      raster_alpha = 1                            # Thickness of the ocean border
    )
    
    bio4.hist<-ratesHistogramGG(
      rates = setNames(
        bioclim_mean$metrics_df$bio4,
        bioclim_mean$metrics_df$h3_cell
      ),                               # Named vector of rates
      colors = custom_palette.bio4$colors,  # Colors from custom_palette
      breaks = custom_palette.bio4$breaks,  # Breaks from custom_palette
      xlab = "BioClim4",          # X-axis label
      ylab = "Counts",                 # Y-axis label
      title = "Temperature Variability",  # Title
      useDensity = FALSE,              # Plot counts (not density)
      logscale = FALSE, lwd=0.01, cex.axis = 1.25, cex.main = 1.5         # Do not use a log scale
    )
    
    bio4_plot.combined <- bio4_plot / ((plot_spacer() | bio4.hist | plot_spacer()) + plot_layout(widths = c(1, 10, 1))) + plot_layout(heights = c(4, 2)) + plot_annotation(title = "Global distribution of BioClim4")
    
    
    # Step 4: Display the plot
    print(bio4_plot)
  }
  bio4_plot
  
  #Visualizing bioclim 15
  {
    
    custom_palette.bio15 <- assignRateColors(
      rates = setNames(
        bioclim_mean$metrics_df$bio15,
        bioclim_mean$metrics_df$h3_cell
      ),
      breaksmethod = "fisher",    # Use "fisher" for better performance with large datasets
      logcolor = FALSE,           # Do not apply log transformation
      palette = "RdYlBu",         # Specify the color palette
      nbreaks = 30,             # Let assignRateColors determine the number of breaks
      reverse = TRUE,              # Reverse the palette
      largeN = 20000,
      samp_prop = 1
    )
    
    # Step 3: Visualize the residuals using the custom plotting function
    bio15_plot <- plot_h3_metrics_global.custom.raster(
      results = bioclim_mean,                          # The constructed results object
      polygons_sf = polygons_sf.1,               # The polygons with H3 indices
      metric_column = "bio15",               # The column for residuals
      custom_palette = custom_palette.bio15$colors,    # Use a viridis palette for residuals
      crs = "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs",  # Coordinate Reference System
      ocean_color = alpha('#87B7DB', 0.5),       # Color for the ocean
      border_thickness = 0.25, 
      lat_graticule_spacing = 20, 
      lon_graticule_spacing = 20, 
      graticule_color = 'white',
      graticule_linetype = 1, 
      graticule_thickness = 0.1,
      raster_basemap = NULL, 
      raster_alpha = 1                            # Thickness of the ocean border
    )
    
    # Step 4: Display the plot
    print(bio15_plot)
  }
  bio15_plot
  
  
}

#generating scatterplots to show spatial effects
{

  #lineage rate v richness
  {
    lineage_rate_v_richness <- generate_scatter_with_colorbar(
      data = climate_vars.cleaned[(order(abs(climate_vars.cleaned$lat))), ],
      x_var = "passerines_counts_range_weighted_z",
      y_var = "mean_lineage_rate_range_weighted", 
      title =  "Effect of species richness on lineage rates",
      color_var = 'custom_palette.bio4', x_label = 'Z-Score local species richness \n (log range-weighted counts)',
      color_palette = custom_palette.bio4$break_colors, xmin_perc = 0.72,
      point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', zero_intercept=F,
    )
    lineage_rate_v_richness <- lineage_rate_v_richness + geom_smooth(method = "lm", formula = y ~ x, color = "red", se = F, level=0.95)
    #lineage_rate_v_richness <- lineage_rate_v_richness + geom_smooth(method = "lm", formula = y ~ poly(x, 1), color = "red", se = T, level=0.95)
    
    #predicted
    lineage_rate_v_richness.predicted <- generate_scatter_with_colorbar(
      data = transform(climate_vars.cleaned, predicted = predict(slm_model.richness.1000idw.lin.con, type="TS"))[(order(abs(climate_vars.cleaned$lat))), ],
      x_var = "passerines_counts_range_weighted_z",
      y_var = "predicted.fit", 
      title =  "SAR-Lag model fit",
      color_var = 'custom_palette.bio4', x_label = 'Z-Score local species richness \n (log range-weighted counts)',
      color_palette = custom_palette.bio4$break_colors, xmin_perc = 0.72,
      point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', zero_intercept=F, y_label = "Predicted range weighted lineage rate", colorbar=F
    ) + geom_smooth(method = "lm", formula = y ~ x, color = "red") +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = ggplot_build(lineage_rate_v_richness)$layout$panel_params[[1]]$y.range) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = "atop(bold(hat(y)) == bold(rho * W * y) + bold(beta[0]) + bold(beta[1] * R) + beta[2] * L^2, 
           + beta[3] * S + beta[4] * R * L^2)",
        parse = TRUE,
        family = "CMU Serif",
        size = 5,
        hjust = 1.05,
        vjust = 1.2
      )
    
    lineage_rate_v_richness | lineage_rate_v_richness.predicted
    
    
  }
  
  # Generate stargazer LaTeX table of LOO metrics for richness model
  {
    
    # Create a cleaned-up copy
    richness_loo_metrics_clean <- slm_model.richness.1000idw.lin.con.loo$metrics %>%
      rename(
        "Term"        = Term,
        "AIC"         = AIC,
        "BIC"         = BIC,
        "ΔAIC"        = DeltaAIC,
        "ΔBIC"        = DeltaBIC,
        "AIC Weight"  = AIC_Weight,
        "BIC Weight"  = BIC_Weight
      ) %>%
      mutate(
        Term = str_replace_all(Term, 
                               "passerines_counts_range_weighted_z", 
                               "Species richness (z)"),
        Term = str_replace_all(Term, 
                               "I\\(lat\\^2\\)", 
                               "Latitude²"),
        Term = str_replace_all(Term, 
                               "Full Model", 
                               "Full model"),
        Term = str_replace_all(Term, 
                               ":", 
                               " × ")
      )
    
    # Print cleaned LaTeX table
    stargazer(
      richness_loo_metrics_clean,
      summary   = FALSE,
      rownames  = FALSE,
      title     = "LOO Model Weights for Richness SAR-Lag Model",
      label     = "tab:richness_model_metrics",
      digits    = 2,
      font.size = "small"
    )
  }
  
  #lineage rate v seasonality
  {
    lineage_rate_v_seasonality <- generate_scatter_with_colorbar(
      data = climate_vars.cleaned[rev(order((climate_vars.cleaned$lat))), ],
      x_var = "wc2.1_5m_bio_4_mean",
      y_var = "mean_lineage_rate_range_weighted", 
      title =  "Effect of Temperature Seasonality on lineage rates",
      color_var = 'custom_palette.lat',
      color_palette = custom_palette.lat$break_colors,
      point_size = 2, point_alpha = 0.5, bar_title = 'Latitude', 
      xmin_perc = 0.85, ymin_perc = 0.8, x_label = 'BioClim4', y_nudge = 0.2,
    ) + ylim(-8.5,-6.2) +
      geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "red", se = F) + ylim(-8.5,-6.2)
    #adding prediction interval
    # lineage_rate_v_seasonality <- lineage_rate_v_seasonality + 
    #   stat_predict(method = "lm", formula = y ~ I(x^2),
    #                geom = "ribbon", fill = alpha('grey',0.1), colour = "grey", linewidth = 0.5, linetype = "solid", level=0.95)  # Confidence interval as solid lines only
    # 
    # lineage_rate_v_seasonality<- ggMarginal(
    #   lineage_rate_v_seasonality,
    #   type = "boxplot", outliers=F,
    #   margins = "both",
    #   size = 20,  # Adjust density plot size
    #   fill = NA,  # No fill, making it fully transparent
    #   color = "black",  # Thin black outline for density
    #   alpha = 1,  # Ensure full opacity for the line
    #   linewidth = 0.1  # Thin density lines
    # )
    
    lineage_rate_v_seasonality.predicted <- generate_scatter_with_colorbar(
      data = transform(climate_vars.cleaned, predicted = predict(slm_model.bio4_z.1000idw.con, type="TS"))[rev(order((climate_vars.cleaned$lat))), ],
      x_var = "wc2.1_5m_bio_4_mean",
      y_var = "predicted.fit", 
      title =  "SAR-Lag model fit",
      color_var = 'custom_palette.lat',
      color_palette = custom_palette.lat$break_colors,
      point_size = 2, point_alpha = 0.5, bar_title = 'Latitude', 
      xmin_perc = 0.85, ymin_perc = 0.8, x_label = 'BioClim4', y_nudge = 0.2, y_label = "Predicted range weighted lineage rate", colorbar = F
    ) + ylim(-8.5,-6.2) + 
      geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "red", se = F) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = "atop(bold(hat(y)) == bold(rho * W * y) + bold(beta[0]) + bold(beta[1] * T) + bold(beta[2] * T^2) + bold(beta[3] * L) + bold(beta[4] * L^2) + beta[5] * S, 
               bold(beta[6] * T * L) + bold(beta[7] * T^2 * L) + beta[8] * T * L^2 + bold(beta[9] * T^2 * L^2))",
        parse = TRUE,
        family = "CMU Serif",
        size = 5,
        hjust = 1.05,
        vjust = 1.2
      )
    
    lineage_rate_v_seasonality | lineage_rate_v_seasonality.predicted
    
    
  }
  
  # Generate stargazer LaTeX table of LOO metrics for temperature seasonality model
  {
    
    # Create a cleaned-up copy
    bio4_loo_metrics_clean <- slm_model.bio4_z.1000idw.con.loo$metrics %>%
      rename(
        "Term"        = Term,
        "AIC"         = AIC,
        "BIC"         = BIC,
        "ΔAIC"        = DeltaAIC,
        "ΔBIC"        = DeltaBIC,
        "AIC Weight"  = AIC_Weight,
        "BIC Weight"  = BIC_Weight
      ) %>%
      mutate(
        Term = str_replace_all(Term, 
                               "poly\\(wc2.1_5m_bio_4_mean_z, 2, raw = T\\)1", 
                               "Temperature seasonality"),
        Term = str_replace_all(Term, 
                               "poly\\(wc2.1_5m_bio_4_mean_z, 2, raw = T\\)2", 
                               "Temperature seasonality²"),
        Term = str_replace_all(Term, 
                               "poly\\(lat_z, 2, raw = T\\)1", 
                               "Latitude"),
        Term = str_replace_all(Term, 
                               "poly\\(lat_z, 2, raw = T\\)2", 
                               "Latitude²"),
        Term = str_replace_all(Term, 
                               "Full Model", 
                               "Full model"),
        Term = str_replace_all(Term, 
                               ":", 
                               " × ")
      )
    
    # Print cleaned LaTeX table
    stargazer(
      bio4_loo_metrics_clean,
      summary   = FALSE,
      rownames  = FALSE,
      title     = "LOO Model Weights for Temperature Seasonality SAR-Lag Model",
      label     = "tab:bio4_model_metrics",
      digits    = 2,
      font.size = "small"
    )
  }
  
  #lineage rate v latitude
  {
    set.seed(1)
    lineage_rate_v_latitude <- generate_scatter_with_colorbar(
      data = climate_vars.cleaned[sample(order((climate_vars.cleaned$lat))), ],
      x_var = "lat",
      y_var = "mean_lineage_rate_range_weighted", 
      title =  "Effect of latitude on lineage rates",
      color_var = 'custom_palette.bio4',
      color_palette = custom_palette.bio4$break_colors,
      point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', 
      x_label = 'Latitude', ymin_perc = 0.8, y_nudge = 0.3, xmin_perc = 0.715,
      x_limits = c(-70,95), zero_intercept = T
    )
    #additional tweaks
    lineage_rate_v_latitude <- lineage_rate_v_latitude + ylim(-8.5,-6.2) + geom_smooth(method = "lm", formula = y ~ I(x^2), color = "red", se = F, level=0.95)
    
    set.seed(1)
    lineage_rate_v_latitude.predicted <- generate_scatter_with_colorbar(
      data = transform(climate_vars.cleaned, predicted = predict(slm_model.lat2.1000idw.con, type="TS"))[sample(order((climate_vars.cleaned$lat))), ],
      x_var = "lat",
      y_var = "predicted.fit", 
      title =  "SAR-Lag model fit",
      color_var = 'custom_palette.bio4',
      color_palette = custom_palette.bio4$break_colors,
      point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', 
      x_label = 'Latitude', ymin_perc = 0.8, y_nudge = 0.3, xmin_perc = 0.715,
      x_limits = c(-70,95), zero_intercept = T, y_label = "Predicted range weighted lineage rate", colorbar = F
    ) + ylim(-8.5,-6.2) + geom_smooth(method = "lm", formula = y ~ I(x^2), color = "red", se = F) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = "bold(hat(y)) == bold(rho * W * y) + bold(beta[0]) + bold(beta[1] * L^2) + beta[2] * S",
        parse = TRUE,
        family = "CMU Serif",
        size = 5,
        hjust = 1.05,
        vjust = 1.2
      )
    
    lineage_rate_v_latitude | lineage_rate_v_latitude.predicted
    
  }

  # Generate stargazer LaTeX table of LOO metrics for latitude model
  {
    
    # Create a cleaned-up copy
    lat_loo_metrics_clean <- slm_model.lat2.1000idw.con.loo$metrics %>%
      rename(
        "Term"        = Term,
        "AIC"         = AIC,
        "BIC"         = BIC,
        "ΔAIC"        = DeltaAIC,
        "ΔBIC"        = DeltaBIC,
        "AIC Weight"  = AIC_Weight,
        "BIC Weight"  = BIC_Weight
      ) %>%
      mutate(
        Term = str_replace_all(Term, 
                               "I\\(lat\\^2\\)", 
                               "Latitude²"),
        Term = str_replace_all(Term, 
                               "Full Model", 
                               "Full model")
      )
    
    # Print cleaned LaTeX table
    stargazer(
      lat_loo_metrics_clean,
      summary   = FALSE,
      rownames  = FALSE,
      title     = "LOO Model Weights for Latitude SAR-Lag Model",
      label     = "tab:latitude_model_metrics",
      digits    = 2,
      font.size = "small"
    )
  }
  
  #lineage rate v latitude- aggregate
  {
    
    set.seed(1)
    lineage_rate_v_latitude.aggregate.predicted <- generate_scatter_with_colorbar(
      data = transform(climate_vars.cleaned, predicted = predict(slm_model.step9c.2c, type="TS"))[sample(order((climate_vars.cleaned$lat))), ],
      x_var = "lat",
      y_var = "predicted.fit", 
      title =  "SAR-Lag model fit",
      color_var = 'custom_palette.bio4',
      color_palette = custom_palette.bio4$break_colors,
      point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', 
      x_label = 'Latitude', ymin_perc = 0.8, y_nudge = 0.3, xmin_perc = 0.715,
      x_limits = c(-70,95), zero_intercept = T, y_label = "Predicted range weighted lineage rate", colorbar = F
    ) + ylim(-8.5,-6.2) + geom_smooth(method = "lm", formula = y ~ I(x^2), color = "red", se = F) 
    # annotate(
    #   "text",
    #   x = Inf, y = Inf,
    #   label = "bold(hat(y)) == bold(rho * W * y) + bold(beta[0]) + bold(beta[1] * L^2) + beta[2] * S",
    #   parse = TRUE,
    #   family = "CMU Serif",
    #   size = 5,
    #   hjust = 1.05,
    #   vjust = 1.2
    # )
    
    lineage_rate_v_latitude | lineage_rate_v_latitude.aggregate.predicted
    
    lineage_rate_v_latitude.predicted | lineage_rate_v_latitude.aggregate.predicted
  }
  
  
  quartz(file='scatter.pdf', width = (10/(1.75))*2, height=25/1.75, type = 'pdf')
  (lineage_rate_v_latitude | lineage_rate_v_latitude.predicted) /
    (lineage_rate_v_seasonality | lineage_rate_v_seasonality.predicted) /
    (lineage_rate_v_richness | lineage_rate_v_richness.predicted)
  dev.off()
  
  lineage_rate_v_latitude
  
  quartz(file='latitude_pars.pdf', width = (10/(1.75))*1.9, height=(25/1.75)/2.75, type = 'pdf')
  (lineage_rate_v_latitude + scale_x_continuous(limits = c(-70, 95), breaks = c(-50, 0, 50), labels = c("-50", "0", "50")) | lineage_rate_v_latitude.predicted) + scale_x_continuous(limits = c(-70, 95), breaks = c(-50, 0, 50), labels = c("-50", "0", "50")) +
    plot_annotation(tag_levels = 'A')
  dev.off()
  
  quartz(file='seasonality_pars.pdf', width = (10/(1.75))*2.1, height=(25/1.75)/2.65, type = 'pdf')
  (lineage_rate_v_seasonality | lineage_rate_v_seasonality.predicted)+
    plot_annotation(tag_levels = 'A')
  dev.off()
  
  quartz(file='richness_pars.pdf', width = (10/(1.75))*2, height=(25/1.75)/2.5, type = 'pdf')
  (lineage_rate_v_richness | lineage_rate_v_richness.predicted)+
    plot_annotation(tag_levels = 'A')
  dev.off()
  
  #richness  v latitude
  {
    set.seed(1)
    richness_v_latitude <- generate_scatter_with_colorbar(
      data = climate_vars.cleaned[sample(order((climate_vars.cleaned$lat))), ],
      x_var = "lat",
      y_var = "passerines_counts_range_weighted_z", 
      title =  "Effect of latitude on local species richness",
      color_var = 'custom_palette.bio4',
      color_palette = custom_palette.bio4$break_colors,
      point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', 
      x_label = 'Latitude', ymin_perc = -.5, y_nudge = 0.3, xmin_perc = 0.625,
      x_limits = c(-70,95), zero_intercept = T, y_label = "Log range weighted counts"
    )
    richness_v_latitude <- richness_v_latitude + geom_smooth(method = "lm", formula = y ~ I(x^2), color = "red", se = F, level=0.95)
    
  }
  
  #alt plots with color bar scaled differently for horizontal plotting
  {
    
    #lineage rate v richness
    {
      lineage_rate_v_richness <- generate_scatter_with_colorbar(
        data = climate_vars.cleaned[(order(abs(climate_vars.cleaned$lat))), ],
        x_var = "passerines_counts_range_weighted_z",
        y_var = "mean_lineage_rate_range_weighted", 
        title =  "Effect of species richness on lineage rates",
        color_var = 'custom_palette.bio4', x_label = 'Z-Score local species richness \n (log range-weighted counts)',
        color_palette = custom_palette.bio4$break_colors, xmin_perc = 0.665, y_nudge = 0.1,
        point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', zero_intercept=F,
      )
      lineage_rate_v_richness <- lineage_rate_v_richness + geom_smooth(method = "lm", formula = y ~ x, color = "red", se = F, level=0.95)
      #lineage_rate_v_richness <- lineage_rate_v_richness + geom_smooth(method = "lm", formula = y ~ poly(x, 1), color = "red", se = T, level=0.95)
      
      
    }
    
    #lineage rate v seasonality
    {
      lineage_rate_v_seasonality <- generate_scatter_with_colorbar(
        data = climate_vars.cleaned[rev(order((climate_vars.cleaned$lat))), ],
        x_var = "wc2.1_5m_bio_4_mean",
        y_var = "mean_lineage_rate_range_weighted", 
        title =  "Effect of Temperature Seasonality on lineage rates",
        color_var = 'custom_palette.lat',
        color_palette = custom_palette.lat$break_colors,
        point_size = 2, point_alpha = 0.5, bar_title = 'Latitude', 
        xmin_perc = 0.815, ymin_perc = 0.81, x_label = 'BioClim4', y_nudge = 0.3
      ) + ylim(-8.5,-6.2) +
        geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "red", se = F) + ylim(-8.5,-6.2)
      
      
    }
    
    #lineage rate v latitude
    {
      set.seed(1)
      lineage_rate_v_latitude <- generate_scatter_with_colorbar(
        data = climate_vars.cleaned[sample(order((climate_vars.cleaned$lat))), ],
        x_var = "lat",
        y_var = "mean_lineage_rate_range_weighted", 
        title =  "Effect of latitude on lineage rates",
        color_var = 'custom_palette.bio4',
        color_palette = custom_palette.bio4$break_colors,
        point_size = 2, point_alpha = 0.5, bar_title = 'Bio4', 
        x_label = 'Latitude', ymin_perc = 0.81, y_nudge = 0.3, xmin_perc = 0.625,
        x_limits = c(-70,95), zero_intercept = T
      )
      #additional tweaks
      lineage_rate_v_latitude <- lineage_rate_v_latitude + ylim(-8.5,-6.2) + geom_smooth(method = "lm", formula = y ~ I(x^2), color = "red", se = F, level=0.95)
      
      
    }
    
  }
  
  
  quartz(file='scatter_long.pdf', height = 4*1, width=14*1.25, type = 'pdf')
  (lineage_rate_v_latitude | lineage_rate_v_seasonality | lineage_rate_v_richness | richness_v_latitude) +
    plot_layout(guides = "collect") & 
    theme(plot.margin = margin(5.5, 15, 5.5, 15))
  dev.off()
  
  
}

#figure 3
{
  # -------------------------------------------------------------------
  # Themes
  # -------------------------------------------------------------------
  map_theme <- theme(
    plot.margin = margin(0, 0, 0, 0)
  )
  
  hist_theme <- theme(
    plot.margin = margin(t = 40, r = 40, b = 20, l = 40)
  )
  
  # -------------------------------------------------------------------
  # Histogram: 98% tall, centered with thin top/bottom spacers
  hist_centered <- function(plot) {
    plot_spacer() / (plot + hist_theme) / plot_spacer() +
      plot_layout(heights = c(0.01, 0.98, 0.01))
  }
  
  # -------------------------------------------------------------------
  # Row 1
  # -------------------------------------------------------------------
  
  row1 <- (
    (predicted_rates_plot.9c.2c + ggtitle(NULL) + map_theme) |
      hist_centered(predicted_rates.hist.9c.2c + ylim(0,620))  
  ) + plot_layout(widths = c(3, 1))
  

  # -------------------------------------------------------------------
  # Row 2
  # -------------------------------------------------------------------
  row2 <- (
    (bio4_plot + ggtitle(NULL) + map_theme) |
      hist_centered(bio4.hist + ylim(0, 950) + xlim(NA,2500))
  ) + plot_layout(widths = c(3, 1))
  
  # -------------------------------------------------------------------
  # Row 3
  # -------------------------------------------------------------------
  row3 <- (
    (residual_rates_plot.9c.2c + ggtitle(NULL) + map_theme) |
      hist_centered(residuals.hist.9c.2c.residuals + xlim(-1.5,1.5)) 
  ) + plot_layout(widths = c(3, 1))
  
  # -------------------------------------------------------------------
  # Combine All Rows Vertically
  # -------------------------------------------------------------------
  combined_plot <- (row1 / row2 / row3) +
    plot_layout(
      heights = c(1, 1, 1),  # Equal height per row
      guides = "collect"     # Shared legend
    )
  
  # -------------------------------------------------------------------
  # Render
  # -------------------------------------------------------------------
  combined_plot
  
  quartz(file='climate_plot_right.pdf', height=11, width=11, type='pdf')
  combined_plot
  dev.off()
  
}

#quantile diagnostics
{
  # =====================================================================================
  # STREAMLINED PRODUCTION SCRIPT (Quantile tails + diagnostics + FINAL 25%-only figure)
  #
  # Goal:
  #   - Keep ALL diagnostic calculations/prints (10% + 25% tails)
  #   - Only build the PLOTS/OBJECTS needed for the FINAL output figure:
  #       Count maps (25%):    map_fast25, map_slow25, leg_fast25, leg_slow25
  #       Fraction maps (25%): map_hot25_frac, map_cold25_frac, leg_hot25_frac, leg_cold25_frac
  #       Diagnostics (25%):   p1_25, p2_25, p3_25, p4_25
  #   - Still computes metrics.fast10/slow10 so diagnostics can run for 10% too.
  #
  # Requires existing objects/functions in your workspace:
  #   Data:
  #     resident_species.1
  #     polygons_sf.1
  #     resident_species.1.h3.metrics$metrics_df with columns:
  #       h3_cell, mean_lineage_rate_range_weighted, species_key_count
  #   Functions:
  #     aggregate_h3_metrics_counts()
  #     assignRateColors()
  #     plot_h3_metrics_global.custom.raster()
  #     generate_color_bar()
  #     generate_scatter_with_colorbar()
  # =====================================================================================
  
    cat("\n=== Running streamlined quantile maps + diagnostics + final figure (25% only output) ===\n")
    
    # ---------------------------
    # Libraries
    # ---------------------------
    library(ggplot2)
    library(patchwork)
    library(scales)
    
    # ---------------------------
    # Sanity checks (core inputs)
    # ---------------------------
    need_objs <- c("resident_species.1", "polygons_sf.1", "resident_species.1.h3.metrics")
    missing <- need_objs[!vapply(need_objs, exists, logical(1))]
    if (length(missing) > 0) stop("Missing required object(s): ", paste(missing, collapse = ", "))
    
    need_fns <- c("aggregate_h3_metrics_counts", "assignRateColors", "plot_h3_metrics_global.custom.raster",
                  "generate_color_bar", "generate_scatter_with_colorbar")
    missing_fns <- need_fns[!vapply(need_fns, exists, logical(1))]
    if (length(missing_fns) > 0) stop("Missing required function(s): ", paste(missing_fns, collapse = ", "))
    
    needed_cols <- c("h3_cell", "mean_lineage_rate_range_weighted", "species_key_count")
    if (!all(needed_cols %in% names(resident_species.1.h3.metrics$metrics_df))) {
      stop("resident_species.1.h3.metrics$metrics_df missing one of: ", paste(needed_cols, collapse = ", "))
    }
    
    # ---------------------------
    # Parameters (single source of truth)
    # ---------------------------
    crs_moll   <- "+proj=moll +lat_0=0 +datum=WGS84 +units=m +no_defs"
    ocean_col  <- alpha("#87B7DB", 0.5)
    
    pal_args <- list(
      breaksmethod = "fisher", logcolor = FALSE, palette = "RdYlBu",
      nbreaks = 20, reverse = TRUE, largeN = 20000, samp_prop = 1
    )
    
    map_leg_widths <- c(14, 0.75)                      # maps wider, legend narrower
    map_side_pad   <- 0.03                             # horizontal padding around each map-panel
    leg_vpad_heights <- c(0.15, 0.70, 0.15)            # shrink colorbars vertically ~30%
    
    diag_left_sp  <- 0.18
    diag_right_sp <- 0.45
    diag_mid      <- 1.00
    diag_row_heights <- c(1, 1, 1.171875, 1.171875)
    diag_pad <- theme(plot.margin = margin(10, 10, 10, 10))
    
    # ---------------------------
    # Small local helpers (no external side effects)
    # ---------------------------
    make_pal <- function(rates_named) {
      do.call(assignRateColors, c(list(rates = rates_named), pal_args))
    }
    
    make_leg <- function(break_colors, bar_title) {
      generate_color_bar(break_colors, bar_title = bar_title, n_ticks = 5) +
        theme(plot.margin = margin(0, 0, 0, 0))
    }
    
    vshrink_legend <- function(leg) {
      (plot_spacer() / leg / plot_spacer()) + plot_layout(heights = leg_vpad_heights)
    }
    
    wrap_map_panel <- function(map_plot, legend_plot_v) {
      core <- (map_plot | legend_plot_v) + plot_layout(widths = map_leg_widths)
      wrap_elements(
        full = ((plot_spacer() | core | plot_spacer()) +
                  plot_layout(widths = c(map_side_pad, 1, map_side_pad))),
        ignore_tag = FALSE, clip = FALSE
      )
    }
    
    wrap_diag_panel <- function(p) {
      wrap_elements(
        full = ((plot_spacer() | p | plot_spacer()) +
                  plot_layout(widths = c(diag_left_sp, diag_mid, diag_right_sp))),
        ignore_tag = FALSE, clip = FALSE
      )
    }
    
    # ===================================================================================
    # A) Quantile tails: compute species-level mean rates and define 10% / 25% tail keys
    # ===================================================================================
    cat("\n--- A) Computing species tail sets (10% + 25%) ---\n")
    
    resident_species.1$speciesKey <- as.character(resident_species.1$speciesKey)
    
    rate_by_key <- tapply(
      resident_species.1$log_lineage_rate,
      resident_species.1$speciesKey,
      mean,
      na.rm = TRUE
    )
    rate_by_key <- rate_by_key[!is.na(rate_by_key)]
    
    q25 <- quantile(rate_by_key, probs = c(0.25, 0.75), na.rm = TRUE)
    q10 <- quantile(rate_by_key, probs = c(0.10, 0.90), na.rm = TRUE)
    
    keys_slow25 <- names(rate_by_key)[rate_by_key <= q25[[1]]]
    keys_fast25 <- names(rate_by_key)[rate_by_key >= q25[[2]]]
    keys_slow10 <- names(rate_by_key)[rate_by_key <= q10[[1]]]
    keys_fast10 <- names(rate_by_key)[rate_by_key >= q10[[2]]]
    
    cat(sprintf("25%% tails: %d slow | %d fast species\n", length(keys_slow25), length(keys_fast25)))
    cat(sprintf("10%% tails: %d slow | %d fast species\n", length(keys_slow10), length(keys_fast10)))
    
    # ===================================================================================
    # B) H3 aggregation (counts): build metrics objects for 25% and 10% tails
    # ===================================================================================
    cat("\n--- B) Aggregating tail richness per H3 cell (counts) ---\n")
    
    resident.slow25 <- resident_species.1[resident_species.1$speciesKey %in% keys_slow25, ]
    resident.fast25 <- resident_species.1[resident_species.1$speciesKey %in% keys_fast25, ]
    resident.slow10 <- resident_species.1[resident_species.1$speciesKey %in% keys_slow10, ]
    resident.fast10 <- resident_species.1[resident_species.1$speciesKey %in% keys_fast10, ]
    
    metrics.slow25 <- aggregate_h3_metrics_counts(resident.slow25, cores = 4, use_sqrt = FALSE, verbose = FALSE)
    metrics.fast25 <- aggregate_h3_metrics_counts(resident.fast25, cores = 4, use_sqrt = FALSE, verbose = FALSE)
    metrics.slow10 <- aggregate_h3_metrics_counts(resident.slow10, cores = 4, use_sqrt = FALSE, verbose = FALSE)
    metrics.fast10 <- aggregate_h3_metrics_counts(resident.fast10, cores = 4, use_sqrt = FALSE, verbose = FALSE)
    
    # ===================================================================================
    # C) COUNT maps (25% only): create map_fast25/map_slow25 + legends
    # ===================================================================================
    cat("\n--- C) Building 25% COUNT maps + legends (objects for final figure) ---\n")
    
    rates_fast25 <- setNames(metrics.fast25$metrics_df$species_key_count, metrics.fast25$metrics_df$h3_cell)
    rates_slow25 <- setNames(metrics.slow25$metrics_df$species_key_count, metrics.slow25$metrics_df$h3_cell)
    
    pal_fast25 <- make_pal(rates_fast25)
    pal_slow25 <- make_pal(rates_slow25)
    
    map_fast25 <- plot_h3_metrics_global.custom.raster(
      results = metrics.fast25, polygons_sf = polygons_sf.1,
      metric_column = "species_key_count", custom_palette = alpha(pal_fast25$colors, 1),
      crs = crs_moll, ocean_color = ocean_col,
      border_thickness = 0.25, graticule_color = "white", graticule_thickness = 0.1
    ) + ggtitle("Top 25% rapidly evolving lineages")
    
    map_slow25 <- plot_h3_metrics_global.custom.raster(
      results = metrics.slow25, polygons_sf = polygons_sf.1,
      metric_column = "species_key_count", custom_palette = alpha(pal_slow25$colors, 1),
      crs = crs_moll, ocean_color = ocean_col,
      border_thickness = 0.25, graticule_color = "white", graticule_thickness = 0.1
    ) + ggtitle("Bottom 25% slowly evolving lineages")
    
    leg_fast25 <- make_leg(pal_fast25$break_colors, "Species count")
    leg_slow25 <- make_leg(pal_slow25$break_colors, "Species count")
    
    # ===================================================================================
    # D) FRACTION maps (25% only): create map_hot25_frac/map_cold25_frac + legends
    # ===================================================================================
    cat("\n--- D) Building 25% FRACTION maps + legends (objects for final figure) ---\n")
    
    total_df <- resident_species.1.h3.metrics$metrics_df[, c("h3_cell","species_key_count")]
    names(total_df)[2] <- "total_richness"
    
    hot25_df <- metrics.fast25$metrics_df[, c("h3_cell","species_key_count")]
    names(hot25_df)[2] <- "hot25_count"
    cold25_df <- metrics.slow25$metrics_df[, c("h3_cell","species_key_count")]
    names(cold25_df)[2] <- "cold25_count"
    
    frac25_df <- merge(total_df, hot25_df, by = "h3_cell", all.x = TRUE)
    frac25_df <- merge(frac25_df, cold25_df, by = "h3_cell", all.x = TRUE)
    frac25_df$hot25_count[is.na(frac25_df$hot25_count)] <- 0
    frac25_df$cold25_count[is.na(frac25_df$cold25_count)] <- 0
    
    frac25_df$hot25_fraction  <- with(frac25_df, ifelse(total_richness > 0, hot25_count  / total_richness, NA_real_))
    frac25_df$cold25_fraction <- with(frac25_df, ifelse(total_richness > 0, cold25_count / total_richness, NA_real_))
    frac25_df$hot25_fraction[is.na(frac25_df$hot25_fraction)] <- 0
    frac25_df$cold25_fraction[is.na(frac25_df$cold25_fraction)] <- 0
    
    metrics.hot25_fraction <- metrics.fast25
    metrics.cold25_fraction <- metrics.slow25
    metrics.hot25_fraction$metrics_df  <- frac25_df[, c("h3_cell","hot25_fraction")]
    metrics.cold25_fraction$metrics_df <- frac25_df[, c("h3_cell","cold25_fraction")]
    
    rates_hot25_frac  <- setNames(metrics.hot25_fraction$metrics_df$hot25_fraction,  metrics.hot25_fraction$metrics_df$h3_cell)
    rates_cold25_frac <- setNames(metrics.cold25_fraction$metrics_df$cold25_fraction, metrics.cold25_fraction$metrics_df$h3_cell)
    
    pal_hot25_frac  <- make_pal(rates_hot25_frac)
    pal_cold25_frac <- make_pal(rates_cold25_frac)
    
    map_hot25_frac <- plot_h3_metrics_global.custom.raster(
      results = metrics.hot25_fraction, polygons_sf = polygons_sf.1,
      metric_column = "hot25_fraction", custom_palette = alpha(pal_hot25_frac$colors, 1),
      crs = crs_moll, ocean_color = ocean_col,
      border_thickness = 0.25, graticule_color = "white", graticule_thickness = 0.1
    ) + ggtitle("Top 25% lineages (fraction of assemblage)")
    
    map_cold25_frac <- plot_h3_metrics_global.custom.raster(
      results = metrics.cold25_fraction, polygons_sf = polygons_sf.1,
      metric_column = "cold25_fraction", custom_palette = alpha(pal_cold25_frac$colors, 1),
      crs = crs_moll, ocean_color = ocean_col,
      border_thickness = 0.25, graticule_color = "white", graticule_thickness = 0.1
    ) + ggtitle("Bottom 25% lineages (fraction of assemblage)")
    
    leg_hot25_frac  <- make_leg(pal_hot25_frac$break_colors,  "Fraction")
    leg_cold25_frac <- make_leg(pal_cold25_frac$break_colors, "Fraction")
    
    # ===================================================================================
    # E) Diagnostics (KEEP ALL calculations): run for BOTH 10% and 25% tails
    #     - Produces p1..p4 (10%) and p1_25..p4_25 (25%) objects
    # ===================================================================================
    cat("\n--- E) Diagnostics for Reviewer 2.4 (10% + 25%) ---\n")
    
    grid_df <- resident_species.1.h3.metrics$metrics_df[, c(
      "h3_cell", "mean_lineage_rate_range_weighted", "species_key_count"
    )]
    names(grid_df) <- c("h3_cell", "mean_rate_rw", "total_richness")
    
    hot10_df <- metrics.fast10$metrics_df[, c("h3_cell", "species_key_count")]
    names(hot10_df)[2] <- "hot10_count"
    cold10_df <- metrics.slow10$metrics_df[, c("h3_cell", "species_key_count")]
    names(cold10_df)[2] <- "cold10_count"
    
    hot25_df <- metrics.fast25$metrics_df[, c("h3_cell", "species_key_count")]
    names(hot25_df)[2] <- "hot25_count"
    cold25_df <- metrics.slow25$metrics_df[, c("h3_cell", "species_key_count")]
    names(cold25_df)[2] <- "cold25_count"
    
    diag_df <- merge(grid_df, hot10_df, by = "h3_cell", all.x = TRUE)
    diag_df <- merge(diag_df, cold10_df, by = "h3_cell", all.x = TRUE)
    diag_df <- merge(diag_df, hot25_df, by = "h3_cell", all.x = TRUE)
    diag_df <- merge(diag_df, cold25_df, by = "h3_cell", all.x = TRUE)
    
    diag_df$hot10_count[is.na(diag_df$hot10_count)] <- 0
    diag_df$cold10_count[is.na(diag_df$cold10_count)] <- 0
    diag_df$hot25_count[is.na(diag_df$hot25_count)] <- 0
    diag_df$cold25_count[is.na(diag_df$cold25_count)] <- 0
    
    diag_df <- diag_df[
      is.finite(diag_df$mean_rate_rw) &
        is.finite(diag_df$total_richness) &
        diag_df$total_richness > 0,
    ]
    
    thr90 <- quantile(diag_df$mean_rate_rw, 0.90, na.rm = TRUE)
    diag_df$high_mean10pct <- diag_df$mean_rate_rw >= thr90
    
    diag_df$custom_palette.black <- "black"
    dummy_break_colors <- list(breaks = c(0, 1), colors = c("black"))
    
    box_theme <- theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.line.y.right = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
      )
    
    for (tail in c(10, 25)) {
      
      exp_prop <- tail / 100
      hot_count_col  <- paste0("hot",  tail, "_count")
      cold_count_col <- paste0("cold", tail, "_count")
      hot_prop_col   <- paste0("hot",  tail, "_prop")
      cold_prop_col  <- paste0("cold", tail, "_prop")
      
      diag_df[[hot_prop_col]]  <- diag_df[[hot_count_col]]  / diag_df$total_richness
      diag_df[[cold_prop_col]] <- diag_df[[cold_count_col]] / diag_df$total_richness
      diag_df[[paste0("hot",  tail, "_enrich")]]  <- diag_df[[hot_prop_col]]  - exp_prop
      diag_df[[paste0("cold", tail, "_enrich")]] <- diag_df[[cold_prop_col]] - exp_prop
      
      cat("\n============================================================\n")
      cat(sprintf("DIAGNOSTICS FOR %d%% TAILS (expected fraction = %.2f)\n", tail, exp_prop))
      cat("============================================================\n")
      
      cat("\n=== A) Do mean rates track HOT/COLD counts? (counts can be richness-driven) ===\n")
      cat(sprintf("Spearman cor(mean_rate_rw, %s)  = %.4f\n", hot_count_col,
                  cor(diag_df$mean_rate_rw, diag_df[[hot_count_col]], method = "spearman")))
      cat(sprintf("Spearman cor(mean_rate_rw, %s) = %.4f\n", cold_count_col,
                  cor(diag_df$mean_rate_rw, diag_df[[cold_count_col]], method = "spearman")))
      
      cat("\n=== B) Do mean rates track HOT/COLD proportions? (key swamping diagnostic) ===\n")
      cat(sprintf("Spearman cor(mean_rate_rw, %s)   = %.4f\n", hot_prop_col,
                  cor(diag_df$mean_rate_rw, diag_df[[hot_prop_col]], method = "spearman")))
      cat(sprintf("Spearman cor(mean_rate_rw, %s)  = %.4f\n", cold_prop_col,
                  cor(diag_df$mean_rate_rw, diag_df[[cold_prop_col]], method = "spearman")))
      
      cat("\n=== C) Are the highest mean-rate cells enriched for hot lineages / depleted of cold? ===\n")
      summ <- data.frame(
        group = c("bottom 90% mean-rate cells", "top 10% mean-rate cells"),
        n = c(sum(!diag_df$high_mean10pct), sum(diag_df$high_mean10pct)),
        mean_hot_prop = c(mean(diag_df[[hot_prop_col]][!diag_df$high_mean10pct]),
                          mean(diag_df[[hot_prop_col]][diag_df$high_mean10pct])),
        median_hot_prop = c(median(diag_df[[hot_prop_col]][!diag_df$high_mean10pct]),
                            median(diag_df[[hot_prop_col]][diag_df$high_mean10pct])),
        mean_cold_prop = c(mean(diag_df[[cold_prop_col]][!diag_df$high_mean10pct]),
                           mean(diag_df[[cold_prop_col]][diag_df$high_mean10pct])),
        median_cold_prop = c(median(diag_df[[cold_prop_col]][!diag_df$high_mean10pct]),
                             median(diag_df[[cold_prop_col]][diag_df$high_mean10pct]))
      )
      print(summ)
      
      cat("\n=== D) Sensitivity: correlations after filtering out very species-poor cells ===\n")
      for (minN in c(1, 10, 20)) {
        sub <- diag_df[diag_df$total_richness >= minN, ]
        cat(sprintf("\n-- min richness >= %d (n=%d) --\n", minN, nrow(sub)))
        cat(sprintf("Spearman cor(mean_rate_rw, %s)  = %.4f\n", hot_prop_col,
                    cor(sub$mean_rate_rw, sub[[hot_prop_col]], method = "spearman")))
        cat(sprintf("Spearman cor(mean_rate_rw, %s) = %.4f\n", cold_prop_col,
                    cor(sub$mean_rate_rw, sub[[cold_prop_col]], method = "spearman")))
      }
      
      # Plots (objects) for this tail set
      set.seed(1)
      p1_tmp <- generate_scatter_with_colorbar(
        data = diag_df[sample(seq_len(nrow(diag_df))), ],
        x_var = hot_prop_col,
        y_var = "mean_rate_rw",
        color_var = "custom_palette.black",
        color_palette = dummy_break_colors,
        title = sprintf("Mean rates vs fast%d proportion (%.2f = global expectation)", tail, exp_prop),
        x_label = sprintf("Fraction of local assemblage in top %d%% fast lineages", tail),
        y_label = "Mean lineage rate (range-weighted; log scale)",
        point_size = 2, point_alpha = 0.5,
        bar_title = NULL,
        x_limits = c(0, 1),
        xmin_perc = 0.80, ymin_perc = 0.80, y_nudge = 0.3,
        zero_intercept = FALSE,
        colorbar = FALSE
      ) +
        geom_vline(xintercept = exp_prop, linetype = 2) +
        geom_smooth(
          data = diag_df,
          mapping = aes(x = .data[[hot_prop_col]], y = mean_rate_rw),
          method = "lm", formula = y ~ x,
          color = "red", se = FALSE,
          inherit.aes = FALSE
        )
      
      set.seed(1)
      p2_tmp <- generate_scatter_with_colorbar(
        data = diag_df[sample(seq_len(nrow(diag_df))), ],
        x_var = cold_prop_col,
        y_var = "mean_rate_rw",
        color_var = "custom_palette.black",
        color_palette = dummy_break_colors,
        title = sprintf("Mean rates vs slow%d proportion (%.2f = global expectation)", tail, exp_prop),
        x_label = sprintf("Fraction of local assemblage in bottom %d%% slow lineages", tail),
        y_label = "Mean lineage rate (range-weighted; log scale)",
        point_size = 2, point_alpha = 0.5,
        bar_title = NULL,
        x_limits = c(0, 1),
        xmin_perc = 0.80, ymin_perc = 0.80, y_nudge = 0.3,
        zero_intercept = FALSE,
        colorbar = FALSE
      ) +
        geom_vline(xintercept = exp_prop, linetype = 2) +
        geom_smooth(
          data = diag_df,
          mapping = aes(x = .data[[cold_prop_col]], y = mean_rate_rw),
          method = "lm", formula = y ~ x,
          color = "red", se = FALSE,
          inherit.aes = FALSE
        )
      
      p3_tmp <- ggplot(diag_df, aes(x = factor(high_mean10pct), y = .data[[hot_prop_col]])) +
        geom_boxplot(outlier.size = 0.3, color = "black", fill = "white") +
        geom_hline(yintercept = exp_prop, linetype = 2) +
        labs(
          x = "Cell in top 10% of mean lineage rate?",
          y = sprintf("fast%d proportion", tail),
          title = sprintf("Are high mean-rate cells enriched for fast%d lineages?", tail)
        ) +
        box_theme
      
      p4_tmp <- ggplot(diag_df, aes(x = factor(high_mean10pct), y = .data[[cold_prop_col]])) +
        geom_boxplot(outlier.size = 0.3, color = "black", fill = "white") +
        geom_hline(yintercept = exp_prop, linetype = 2) +
        labs(
          x = "Cell in top 10% of mean lineage rate?",
          y = sprintf("slow%d proportion", tail),
          title = sprintf("Are high mean-rate cells depleted of slow%d lineages?", tail)
        ) +
        box_theme
      
      if (tail == 10) {
        p1 <- p1_tmp; p2 <- p2_tmp; p3 <- p3_tmp; p4 <- p4_tmp
      } else {
        p1_25 <- p1_tmp; p2_25 <- p2_tmp; p3_25 <- p3_tmp; p4_25 <- p4_tmp
      }
    }
    
    cat("\nDone: diagnostics (10% + 25%) and required plot objects created.\n")
    
    # ===================================================================================
    # F) FINAL 25%-only summary figure (your tuned layout)
    # ===================================================================================
    cat("\n--- F) Writing FINAL 25%-only summary figure PDF ---\n")
    
    # Shrink legends vertically (but keep the original legend objects intact)
    leg_fast25_v <- vshrink_legend(leg_fast25)
    leg_slow25_v <- vshrink_legend(leg_slow25)
    leg_hot25_v  <- vshrink_legend(leg_hot25_frac)
    leg_cold25_v <- vshrink_legend(leg_cold25_frac)
    
    panel_ct_hot25  <- wrap_map_panel(map_fast25,     leg_fast25_v)
    panel_ct_cold25 <- wrap_map_panel(map_slow25,     leg_slow25_v)
    panel_fr_hot25  <- wrap_map_panel(map_hot25_frac, leg_hot25_v)
    panel_fr_cold25 <- wrap_map_panel(map_cold25_frac,leg_cold25_v)
    
    # Diagnostics label wrapping (E/F: title + x + y; G/H: title only)
    p1_25_adj <- p1_25 +
      labs(
        title = "Mean rates vs fast25 proportion\n(0.25 = global expectation)",
        x = "Fraction of local assemblage\nin top 25% fast lineages",
        y = "Mean lineage rate\n(range-weighted; log scale)"
      ) + diag_pad
    
    p2_25_adj <- p2_25 +
      labs(
        title = "Mean rates vs slow25 proportion\n(0.25 = global expectation)",
        x = "Fraction of local assemblage\nin bottom 25% slow lineages",
        y = "Mean lineage rate\n(range-weighted; log scale)"
      ) + diag_pad
    
    # Convert p3_25 / p4_25 to vioplot-like violins WITHOUT rebuilding from scratch
    p3_25_vio <- p3_25
    p3_25_vio$layers <- Filter(
      function(l) !(inherits(l$geom, "GeomBoxplot") || inherits(l$geom, "GeomHline") || inherits(l$geom, "GeomViolin")),
      p3_25_vio$layers
    )
    p3_25_vio <- p3_25_vio +
      geom_violin(trim = TRUE, scale = "width", adjust = 3, color = "black", fill = "grey92") +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "black", color = "black") +
      stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, fill = "white", color = "white") +
      geom_hline(yintercept = 0.25, linetype = 2) +
      coord_cartesian(ylim = c(0, 1))
    
    p4_25_vio <- p4_25
    p4_25_vio$layers <- Filter(
      function(l) !(inherits(l$geom, "GeomBoxplot") || inherits(l$geom, "GeomHline") || inherits(l$geom, "GeomViolin")),
      p4_25_vio$layers
    )
    p4_25_vio <- p4_25_vio +
      geom_violin(trim = TRUE, scale = "width", adjust = 3, color = "black", fill = "grey92") +
      geom_boxplot(width = 0.12, outlier.shape = NA, fill = "black", color = "black") +
      stat_summary(fun = median, geom = "point", shape = 21, size = 2.2, fill = "white", color = "white") +
      geom_hline(yintercept = 0.25, linetype = 2) +
      coord_cartesian(ylim = c(0, 1))
    
    p3_25_adj <- p3_25_vio + labs(title = "Are high mean-rate cells enriched\nfor fast25 lineages?") + diag_pad
    p4_25_adj <- p4_25_vio + labs(title = "Are high mean-rate cells depleted\nof slow25 lineages?") + diag_pad
    
    panel_p1_25 <- wrap_diag_panel(p1_25_adj)
    panel_p2_25 <- wrap_diag_panel(p2_25_adj)
    panel_p3_25 <- wrap_diag_panel(p3_25_adj)
    panel_p4_25 <- wrap_diag_panel(p4_25_adj)
    
    final_fig <-
      (panel_ct_hot25 | panel_ct_cold25) /
      (panel_fr_hot25 | panel_fr_cold25) /
      (panel_p1_25 | panel_p2_25) /
      (panel_p3_25 | panel_p4_25)
    
    final_fig <- final_fig +
      plot_layout(heights = diag_row_heights) +
      plot_annotation(tag_levels = "A") &
      theme(
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0.01, 0.99)
      )
    
    pdf("SuppFig17_summary_25only_counts_fractions_diagnostics.pdf", height = 16, width = 16)
    print(final_fig)
    dev.off()
    
    cat("\nSaved: SuppFig17_summary_25only_counts_fractions_diagnostics.pdf\n")
    cat("=== DONE ===\n")
  }
