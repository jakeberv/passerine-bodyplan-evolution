# Master R code for executing analyses for
# Rates of Avian Body Plan Evolution in Space and Time
# Temporal Analyses (primarily)

#load necessary packages
library(ape)
library(phytools)
library(mvMORPH)
library(future)
library(future.apply)
library(viridis)
library(parallel)
library(palaeoverse)
library(boot)
library(RColorBrewer)
library(classInt)
library(pbmcapply)
library(phylolm)
library(scales)
library(tidyverse)
library(patchwork)
library(readxl)
library(univariateML)
library(evd)
library(fitdistrplus)
library(HDInterval)

setwd('/Users/cotinga/jacob.berv@gmail.com/Code/passerine-bodyplan-evolution')

#read in function definitions
source("TemporalAnalyses-functions.R")

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

#fitting main models
#BIC/BIC of 20,40, with minimum clade sizes of 30,20,10
{
#IC 40
{
  ### read in various datasets
  supertree <- read.tree(file='/Users/cotinga/jacob.berv@gmail.com/Code/skelevision_analysis/ythlida_supertree.rescale.tre')
  dat.mvgls <- readRDS(file='/Users/cotinga/jacob.berv@gmail.com/Code/skelevision_analysis/dat.mvgls.RDS')
  sampled_cv <- readRDS(file= '/Users/cotinga/jacob.berv@gmail.com/Code/skelevision_analysis/sampled_cv.RDS')
  superdat.shuffle<-as.matrix(dat.mvgls[supertree$tip.label,]) #reorder the dataset
  #need to reorder the sampled_cv dataset
  rownames(sampled_cv) <- sampled_cv$phylo
  sampled_cv <- sampled_cv[supertree$tip.label,]
  superdat.shuffle.factor <- cbind(superdat.shuffle, tlevel = sampled_cv$tlevel, tniche = sampled_cv$tniche)
  
  #options(future.globals.maxSize = 20 * 1024^3)
  # min30.ic40.gic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 30,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 40,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 40,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="GIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #     )
  # saveRDS(min30.ic40.gic, file='min30.ic40.gic.RDS')
  min30.ic40.gic <- readRDS(file='new_bifrost/min30.ic40.gic.RDS')
  min30.ic40.gic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min30.ic40.gic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 20 * 1024^3)
  # min30.ic40.bic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 30,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 40,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 40,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="BIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min30.ic40.bic, file='min30.ic40.bic.RDS')
  min30.ic40.bic <- readRDS(file='new_bifrost/min30.ic40.bic.RDS')
  min30.ic40.bic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min30.ic40.bic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 30 * 1024^3)
  # min20.ic40.gic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 20,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 40,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 40,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="GIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min20.ic40.gic, file='min20.ic40.gic.RDS')
  min20.ic40.gic<-readRDS(file='new_bifrost/min20.ic40.gic.RDS')
  min20.ic40.gic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min20.ic40.gic$model_fit_history$ic_acceptance_matrix)
  #searchConvergentShifts(input_model_results = min10.ic40.bic, trait_data = superdat.shuffle, formula = 'trait_data[,c(-13)]~trait_data[,c(13)]', convergence_threshold = 1, IC = 'BIC', plot = T, model = "BMM", method='LL', error=T, REML=T, distance_metric = 'AverageVariance')
  
  #options(future.globals.maxSize = 30 * 1024^3)
  # min20.ic40.bic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 20,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 40,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 40,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="BIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min20.ic40.bic, file='min20.ic40.bic.RDS')
  min20.ic40.bic<-readRDS(file='new_bifrost/min20.ic40.bic.RDS')
  min20.ic40.bic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min20.ic40.bic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 35 * 1024^3)
  # min10.ic40.gic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 10,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 40,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 40,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="GIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min10.ic40.gic, file='min10.ic40.gic.RDS')
  min10.ic40.gic <- readRDS(file='new_bifrost/min10.ic40.gic.RDS')
  min10.ic40.gic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min10.ic40.gic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 35 * 1024^3)
  # min10.ic40.bic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 10,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 40,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 40,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="BIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min10.ic40.bic, file='min10.ic40.bic.RDS')
  min10.ic40.bic <- readRDS(file='new_bifrost/min10.ic40.bic.RDS')
  min10.ic40.bic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min10.ic40.bic$model_fit_history$ic_acceptance_matrix)
  
  
}

#IC 20
{
  ### read in various datasets
  supertree <- read.tree(file='/Users/cotinga/jacob.berv@gmail.com/Code/skelevision_analysis/ythlida_supertree.rescale.tre')
  dat.mvgls <- readRDS(file='/Users/cotinga/jacob.berv@gmail.com/Code/skelevision_analysis/dat.mvgls.RDS')
  sampled_cv <- readRDS(file= '/Users/cotinga/jacob.berv@gmail.com/Code/skelevision_analysis/sampled_cv.RDS')
  superdat.shuffle<-as.matrix(dat.mvgls[supertree$tip.label,]) #reorder the dataset
  #need to reorder the sampled_cv dataset
  rownames(sampled_cv) <- sampled_cv$phylo
  sampled_cv <- sampled_cv[supertree$tip.label,]
  superdat.shuffle.factor <- cbind(superdat.shuffle, tlevel = sampled_cv$tlevel, tniche = sampled_cv$tniche)
  
  #options(future.globals.maxSize = 20 * 1024^3)
  # min30.ic20.gic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 30,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 20,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 20,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="GIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min30.ic20.gic, file='min30.ic20.gic.RDS')
  min30.ic20.gic <- readRDS(file='new_bifrost/min30.ic20.gic.RDS')
  min30.ic20.gic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min30.ic20.gic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 20 * 1024^3)
  # min30.ic20.bic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 30,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 20,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 20,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="BIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min30.ic20.bic, file='min30.ic20.bic.RDS')
  min30.ic20.bic <- readRDS(file='new_bifrost/min30.ic20.bic.RDS')
  min30.ic20.bic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min30.ic20.bic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 30 * 1024^3)
  # min20.ic20.gic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 20,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 20,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 20,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="GIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min20.ic20.gic, file='min20.ic20.gic.RDS')
  min20.ic20.gic<-readRDS(file='new_bifrost/min20.ic20.gic.RDS')
  min20.ic20.gic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min20.ic20.gic$model_fit_history$ic_acceptance_matrix)
  #searchConvergentShifts(input_model_results = min10.ic20.bic, trait_data = superdat.shuffle, formula = 'trait_data[,c(-13)]~trait_data[,c(13)]', convergence_threshold = 1, IC = 'BIC', plot = T, model = "BMM", method='LL', error=T, REML=T, distance_metric = 'AverageVariance')
  
  #options(future.globals.maxSize = 30 * 1024^3)
  # min20.ic20.bic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 20,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 20,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 20,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="BIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min20.ic20.bic, file='min20.ic20.bic.RDS')
  min20.ic20.bic<-readRDS(file='new_bifrost/min20.ic20.bic.RDS')
  min20.ic20.bic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min20.ic20.bic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 35 * 1024^3)
  # min10.ic20.gic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 10,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 20,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 20,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="GIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min10.ic20.gic, file='min10.ic20.gic.RDS')
  min10.ic20.gic <- readRDS(file='new_bifrost/min10.ic20.gic.RDS')
  min10.ic20.gic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min10.ic20.gic$model_fit_history$ic_acceptance_matrix)
  
  #options(future.globals.maxSize = 35 * 1024^3)
  # min10.ic20.bic <-
  #   searchOptimalConfiguration.mclapply(
  #     baseline_tree = paintSubTree(
  #       reorder(as.phylo(supertree), order='postorder'),
  #       node = length(supertree$tip.label) + 1,
  #       state = 0),
  #     trait_data = superdat.shuffle,
  #     min_descendant_tips = 10,
  #     num_cores = 20,
  #     ic_uncertainty_threshold = 20,
  #     uncertainty = F,
  #     shift_acceptance_threshold = 20,
  #     postorder_traversal = F,
  #     uncertaintyweights = F,
  #     uncertaintyweights_par = T,
  #     IC="BIC",
  #     plot=F,
  #     formula='trait_data[,c(-13)]~trait_data[,c(13)]', 
  #     method='LL', 
  #     error=T
  #   )
  # saveRDS(min10.ic20.bic, file='min10.ic20.bic.RDS')
  min10.ic20.bic <- readRDS(file='new_bifrost/min10.ic20.bic.RDS')
  min10.ic20.bic$model_fit_history$fits<-NULL #wipe out the history to save memory
  plot_ic_acceptance_matrix(matrix_data = min10.ic20.bic$model_fit_history$ic_acceptance_matrix)
  
  
}

}
  
#store a baseline model for reference later
{
#fitting BM model to the full dataset
GIC(mvgls(
  formula = superdat.shuffle[,c(-13)]~superdat.shuffle[,c(13)],
  tree =  paintSubTree(
    reorder(as.phylo(supertree), order = 'postorder'),
    node = length(supertree$tip.label) + 1,
    state = 0
  ), model='BM', method='LL', error=T
)) #GIC: -141939.1, BIC -141495.8

global.model.mvBM<-mvgls(
  formula = superdat.shuffle[,c(-13)]~superdat.shuffle[,c(13)],
  tree =  paintSubTree(
    reorder(as.phylo(supertree), order = 'postorder'),
    node = length(supertree$tip.label) + 1,
    state = 0
  ), model='BM', method='LL', error=T
)
}

#simulations for model performance
{
#estimate mean vars and means covar for simulations
{
  tmp<-mvgls(
    formula = superdat.shuffle[,c(-13)]~superdat.shuffle[,c(13)],
    tree =  paintSubTree(
      reorder(as.phylo(supertree), order = 'postorder'),
      node = length(supertree$tip.label) + 1,
      state = 0
    ), model='BM', method='LL'
  )
  meanvar<-mean(diag(tmp$sigma$Pinv)) #0.0002500118
  sdvar <- sd(diag(tmp$sigma$Pinv)) #0.0001372484
  diag(tmp$sigma$Pinv)<-NA
  meancov<-mean(tmp$sigma$Pinv[upper.tri(tmp$sigma$Pinv)]) #2.248166e-05
  sdcov <- sd(tmp$sigma$Pinv[upper.tri(tmp$sigma$Pinv)]) #0.000121608
}

#simulations for false positives
{
  
  set.seed(seed=1)
  simdata<-simulate_traits_BM1(n_species = 100, n_traits = 10)
  
  
  #100 tips
  {
    
    #GIC threshold of 2
    {
      # simresults_100tip_10D_GIC2 <- run_FP_shift_inference(
      #   base_tree = supertree,
      #   n_datasets     = 100,
      #   n_traits       = 10,
      #   tree_tip_count = 100,
      #   num_cores      = 8,
      #   search_options = list(
      #     formula                     = "trait_data~1",
      #     min_descendant_tips         = 10,
      #     ic_uncertainty_threshold    = 10,
      #     shift_acceptance_threshold  = 2,
      #     uncertainty                 = FALSE,
      #     uncertaintyweights          = FALSE,
      #     uncertaintyweights_par      = FALSE,
      #     postorder_traversal         = FALSE,
      #     plot                        = FALSE,
      #     IC                          = "GIC",
      #     store_model_fit_history     = FALSE,
      #     num_cores                   = 1, 
      #     method = 'LL'
      #   ),
      #   simulation_options = list(
      #     variance_mean    = 0.00025,
      #     variance_sd      = 0.0001372484,
      #     covariance_mean  = 2.248166e-05,
      #     covariance_sd    = 0.000121608
      #   ),
      #   seed = 5
      # )
      # saveRDS(simresults_100tip_10D_GIC2, file="simresults_100tip_10D_GIC2.RDS")
      simresults_100tip_10D_GIC2 <- readRDS(file = "simresults_100tip_10D_GIC2.RDS")
      
      #what is the FP rate?
      FP_list.1.gic<-list()
      for(i in 1:length(simresults_100tip_10D_GIC2)){
        FP_list.1.gic[i]<-length(simresults_100tip_10D_GIC2[[i]]$shift_nodes_no_uncertainty)/simresults_100tip_10D_GIC2[[i]]$num_candidates
      }
      mean(unlist(FP_list.1.gic))
      
      FP.1.gic <- plot_fp_histogram(FP_list.1.gic,
                                    title = "False Positives. 100 tips, 10D, GIC: 2, min: 10",
                                    xlab = "FP rate",
                                    ylab = "Tree resamples")
      
      FP.1.gic
      
    }
    
    #GIC threshold of 10
    {
      # simresults_100tip_10D_GIC10 <- run_FP_shift_inference(
      #   base_tree = supertree,
      #   n_datasets = 100,
      #   n_traits = 10,
      #   tree_tip_count = 100,
      #   num_cores = 8,
      #   search_options = list(
      #     formula = "trait_data~1",
      #     min_descendant_tips = 10,
      #     ic_uncertainty_threshold = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty = FALSE,
      #     uncertaintyweights = FALSE,
      #     uncertaintyweights_par = FALSE,
      #     postorder_traversal = FALSE,
      #     plot = FALSE,
      #     IC = "GIC",
      #     store_model_fit_history = FALSE,
      #     num_cores = 1, 
      #     method = 'LL'
      #   ),
      #   simulation_options = list(
      #     variance_mean = 0.00025,
      #     variance_sd   = 0.0001372,
      #     covariance_mean = 0.0000225,
      #     covariance_sd   = 0.0001216
      #   ),
      #   seed = 5
      # )
      # saveRDS(simresults_100tip_10D_GIC10, file="simresults_100tip_10D_GIC10.RDS")
      simresults_100tip_10D_GIC10<- readRDS(file="simresults_100tip_10D_GIC10.RDS")
      
      #what is the FP rate?
      FP_list.2.gic<-list()
      for(i in 1:length(simresults_100tip_10D_GIC10)){
        FP_list.2.gic[i]<-length(simresults_100tip_10D_GIC10[[i]]$shift_nodes_no_uncertainty)/simresults_100tip_10D_GIC10[[i]]$num_candidates
      }
      mean(unlist(FP_list.2.gic))
      
      FP.2.gic <- plot_fp_histogram(FP_list.2.gic,
                                    title = "False Positives. 100 tips, 10D, GIC: 10, min: 10",
                                    xlab = "FP rate",
                                    ylab = "Tree resamples")
      
      FP.2.gic
      
      
    }
    
  }
  
  #simulation pipeline, 200 tips
  {
    # simresults_200tip_10D_GIC2 <- run_FP_shift_inference(
    #   base_tree = supertree,
    #   n_datasets        = 100,
    #   n_traits          = 10,
    #   tree_tip_count    = 200,
    #   num_cores         = 8,
    #   search_options    = list(
    #     formula = "trait_data~1",
    #     min_descendant_tips = 10,
    #     ic_uncertainty_threshold = 10,
    #     shift_acceptance_threshold = 2,
    #     uncertainty = FALSE,
    #     uncertaintyweights = FALSE,
    #     uncertaintyweights_par = FALSE,
    #     postorder_traversal = FALSE,
    #     plot = FALSE,
    #     IC = "GIC",
    #     store_model_fit_history = FALSE,
    #     num_cores = 1, 
    #     method = 'LL'
    #   ),
    #   simulation_options = list(
    #     variance_mean   = 0.00025,
    #     variance_sd     = 0.00014,
    #     covariance_mean = 2.25e-05,
    #     covariance_sd   = 0.00012
    #   ),
    #   seed = 5
    # )
    # saveRDS(simresults_200tip_10D_GIC2, file="simresults_200tip_10D_GIC2.RDS")
    simresults_200tip_10D_GIC2 <- readRDS(file="simresults_200tip_10D_GIC2.RDS")
    
    #what is the FP rate?
    FP_list.3.gic<-list()
    for(i in 1:length(simresults_200tip_10D_GIC2)){
      FP_list.3.gic[i]<-length(simresults_200tip_10D_GIC2[[i]]$shift_nodes_no_uncertainty)/simresults_200tip_10D_GIC2[[i]]$num_candidates
    }
    mean(unlist(FP_list.3.gic))
    
    FP.3.gic <- plot_fp_histogram(FP_list.3.gic,
                                  title = "False Positives. 200 tips, 10D, GIC: 2, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.3.gic
    
    
    # simresults_200tip_10D_GIC10 <- run_FP_shift_inference(
    #   base_tree = supertree,
    #   n_datasets        = 100,
    #   n_traits          = 10,
    #   tree_tip_count    = 200,
    #   num_cores         = 8,
    #   search_options    = list(
    #     formula = "trait_data~1",
    #     min_descendant_tips = 10,
    #     ic_uncertainty_threshold = 10,
    #     shift_acceptance_threshold = 10,
    #     uncertainty = FALSE,
    #     uncertaintyweights = FALSE,
    #     uncertaintyweights_par = FALSE,
    #     postorder_traversal = FALSE,
    #     plot = FALSE,
    #     IC = "GIC",
    #     store_model_fit_history = FALSE,
    #     num_cores = 1, 
    #     method = 'LL'
    #   ),
    #   simulation_options = list(
    #     variance_mean   = 0.00025,
    #     variance_sd     = 0.00014,
    #     covariance_mean = 2.25e-05,
    #     covariance_sd   = 0.00012
    #   ),
    #   seed = 5
    # )
    # saveRDS(simresults_200tip_10D_GIC10, file='simresults_200tip_10D_GIC10.RDS')
    simresults_200tip_10D_GIC10<-readRDS(file='simresults_200tip_10D_GIC10.RDS')
    
    #what is the FP rate?
    FP_list.4.gic<-list()
    for(i in 1:length(simresults_200tip_10D_GIC10)){
      FP_list.4.gic[i]<-length(simresults_200tip_10D_GIC10[[i]]$shift_nodes_no_uncertainty)/simresults_200tip_10D_GIC10[[i]]$num_candidates
    }
    mean(unlist(FP_list.4.gic))
    
    FP.4.gic <- plot_fp_histogram(FP_list.4.gic,
                                  title = "False Positives. 200 tips, 10D, GIC: 10, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.4.gic
    
  }
  
  #simulation pipeline, 300 tips
  {
    
    # simresults_300tip_10D_GIC2 <- run_FP_shift_inference(
    #   base_tree = supertree,
    #   n_datasets        = 100,
    #   n_traits          = 10,
    #   tree_tip_count    = 300,
    #   num_cores         = 8,
    #   search_options    = list(
    #     formula = "trait_data~1",
    #     min_descendant_tips = 10,
    #     ic_uncertainty_threshold = 10,
    #     shift_acceptance_threshold = 2,
    #     uncertainty = FALSE,
    #     uncertaintyweights = FALSE,
    #     uncertaintyweights_par = FALSE,
    #     postorder_traversal = FALSE,
    #     plot = FALSE,
    #     IC = "GIC",
    #     store_model_fit_history = FALSE,
    #     num_cores = 1, 
    #     method = 'LL'
    #   ),
    #   simulation_options = list(
    #     variance_mean   = 0.00025,
    #     variance_sd     = 0.00014,
    #     covariance_mean = 2.25e-05,
    #     covariance_sd   = 0.00012
    #   ),
    #   seed = 5
    # )
    # saveRDS(simresults_300tip_10D_GIC2, file="simresults_300tip_10D_GIC2.RDS")
    simresults_300tip_10D_GIC2<-readRDS('simresults_300tip_10D_GIC2.RDS')
    
    #what is the FP rate?
    FP_list.5.gic<-list()
    for(i in 1:length(simresults_300tip_10D_GIC2)){
      FP_list.5.gic[i]<-length(simresults_300tip_10D_GIC2[[i]]$shift_nodes_no_uncertainty)/simresults_300tip_10D_GIC2[[i]]$num_candidates
    }
    mean(unlist(FP_list.5.gic))
    
    FP.5.gic <- plot_fp_histogram(FP_list.5.gic,
                                  title = "False Positives. 300 tips, 10D, GIC: 2, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.5.gic
    
    
    # simresults_300tip_10D_GIC10 <- run_FP_shift_inference(
    #   base_tree = supertree,
    #   n_datasets        = 100,
    #   n_traits          = 10,
    #   tree_tip_count    = 300,
    #   num_cores         = 8,
    #   search_options    = list(
    #     formula = "trait_data~1",
    #     min_descendant_tips = 10,
    #     ic_uncertainty_threshold = 10,
    #     shift_acceptance_threshold = 10,
    #     uncertainty = FALSE,
    #     uncertaintyweights = FALSE,
    #     uncertaintyweights_par = FALSE,
    #     postorder_traversal = FALSE,
    #     plot = FALSE,
    #     IC = "GIC",
    #     store_model_fit_history = FALSE,
    #     num_cores = 1, 
    #     method = 'LL'
    #   ),
    #   simulation_options = list(
    #     variance_mean   = 0.00025,
    #     variance_sd     = 0.00014,
    #     covariance_mean = 2.25e-05,
    #     covariance_sd   = 0.00012
    #   ),
    #   seed = 5
    # )
    # saveRDS(simresults_300tip_10D_GIC10, file="simresults_300tip_10D_GIC10.RDS")
    simresults_300tip_10D_GIC10 <- readRDS(file="simresults_300tip_10D_GIC10.RDS")
    
    #what is the FP rate?
    FP_list.6.gic<-list()
    for(i in 1:length(simresults_300tip_10D_GIC10)){
      FP_list.6.gic[i]<-length(simresults_300tip_10D_GIC10[[i]]$shift_nodes_no_uncertainty)/simresults_300tip_10D_GIC10[[i]]$num_candidates
    }
    mean(unlist(FP_list.6.gic))
    
    FP.6.gic <- plot_fp_histogram(FP_list.6.gic,
                                  title = "False Positives. 300 tips, 10D, GIC: 10, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.6.gic
    
    
  }
  
  quartz(width = 10, height = 10, type = 'pdf', file = "FP_histograms.GIC.pdf")
  
  # Combine with alphabetical panel tags
  (
    (FP.1.gic | FP.2.gic) /
      (FP.3.gic | FP.4.gic) /
      (FP.5.gic | FP.6.gic)
  ) + plot_annotation(tag_levels = "A")
  
  dev.off()
  
  #BIC
  #100 tips
  {
    
    #BIC threshold of 2
    {
      simresults_100tip_10D_BIC2 <- run_FP_shift_inference(
        base_tree = supertree,
        n_datasets     = 100,
        n_traits       = 10,
        tree_tip_count = 100,
        num_cores      = 8,
        search_options = list(
          formula                     = "trait_data~1",
          min_descendant_tips         = 10,
          ic_uncertainty_threshold    = 10,
          shift_acceptance_threshold  = 2,
          uncertainty                 = FALSE,
          uncertaintyweights          = FALSE,
          uncertaintyweights_par      = FALSE,
          postorder_traversal         = FALSE,
          plot                        = FALSE,
          IC                          = "BIC",
          store_model_fit_history     = FALSE,
          num_cores                   = 1,
          method = 'LL'
        ),
        simulation_options = list(
          variance_mean    = 0.00025,
          variance_sd      = 0.0001372484,
          covariance_mean  = 2.248166e-05,
          covariance_sd    = 0.000121608
        ),
        seed = 5
      )
      saveRDS(simresults_100tip_10D_BIC2, file="simresults_100tip_10D_BIC2.RDS")
      simresults_100tip_10D_BIC2 <- readRDS("simresults_100tip_10D_BIC2.RDS")
      
      #what is the FP rate?
      FP_list.1.bic<-list()
      for(i in 1:length(simresults_100tip_10D_GIC2)){
        FP_list.1.bic[i]<-length(simresults_100tip_10D_BIC2[[i]]$shift_nodes_no_uncertainty)/simresults_100tip_10D_BIC2[[i]]$num_candidates
      }
      mean(unlist(FP_list.1.bic))
      
      FP.1.bic <- plot_fp_histogram(FP_list.1.bic,
                                    title = "False Positives. 100 tips, 10D, BIC: 2, min: 10",
                                    xlab = "FP rate",
                                    ylab = "Tree resamples")
      
      FP.1.bic
      
    }
    
    #BIC threshold of 10
    {
      simresults_100tip_10D_BIC10 <- run_FP_shift_inference(
        base_tree = supertree,
        n_datasets = 100,
        n_traits = 10,
        tree_tip_count = 100,
        num_cores = 8,
        search_options = list(
          formula = "trait_data~1",
          min_descendant_tips = 10,
          ic_uncertainty_threshold = 10,
          shift_acceptance_threshold = 10,
          uncertainty = FALSE,
          uncertaintyweights = FALSE,
          uncertaintyweights_par = FALSE,
          postorder_traversal = FALSE,
          plot = FALSE,
          IC = "BIC",
          store_model_fit_history = FALSE,
          num_cores = 1, 
          method = 'LL'
        ),
        simulation_options = list(
          variance_mean = 0.00025,
          variance_sd   = 0.0001372,
          covariance_mean = 0.0000225,
          covariance_sd   = 0.0001216
        ),
        seed = 5
      )
      saveRDS(simresults_100tip_10D_BIC10, file="simresults_100tip_10D_BIC10.RDS")
      simresults_100tip_10D_BIC10<- readRDS(file='simresults_100tip_10D_BIC10.RDS')
      
      #what is the FP rate?
      FP_list.2.bic<-list()
      for(i in 1:length(simresults_100tip_10D_BIC10)){
        FP_list.2.bic[i]<-length(simresults_100tip_10D_BIC10[[i]]$shift_nodes_no_uncertainty)/simresults_100tip_10D_BIC10[[i]]$num_candidates
      }
      mean(unlist(FP_list.2.bic))
      
      FP.2.bic <- plot_fp_histogram(FP_list.2.bic,
                                    title = "False Positives. 100 tips, 10D, BIC: 10, min: 10",
                                    xlab = "FP rate",
                                    ylab = "Tree resamples")
      
      FP.2.bic
      
      
    }
    
  }
  
  #simulation pipeline, 200 tips
  {
    simresults_200tip_10D_BIC2 <- run_FP_shift_inference(
      base_tree = supertree,
      n_datasets        = 100,
      n_traits          = 10,
      tree_tip_count    = 200,
      num_cores         = 8,
      search_options    = list(
        formula = "trait_data~1",
        min_descendant_tips = 10,
        ic_uncertainty_threshold = 10,
        shift_acceptance_threshold = 2,
        uncertainty = FALSE,
        uncertaintyweights = FALSE,
        uncertaintyweights_par = FALSE,
        postorder_traversal = FALSE,
        plot = FALSE,
        IC = "BIC",
        store_model_fit_history = FALSE,
        num_cores = 1,
        method = 'LL'
      ),
      simulation_options = list(
        variance_mean   = 0.00025,
        variance_sd     = 0.00014,
        covariance_mean = 2.25e-05,
        covariance_sd   = 0.00012
      ),
      seed = 5
    )
    saveRDS(simresults_200tip_10D_BIC2, file="simresults_200tip_10D_BIC2.RDS")
    simresults_200tip_10D_BIC2 <- readRDS("simresults_200tip_10D_BIC2.RDS")
    
    #what is the FP rate?
    FP_list.3.bic<-list()
    for(i in 1:length(simresults_200tip_10D_BIC2)){
      FP_list.3.bic[i]<-length(simresults_200tip_10D_BIC2[[i]]$shift_nodes_no_uncertainty)/simresults_200tip_10D_BIC2[[i]]$num_candidates
    }
    mean(unlist(FP_list.3.bic))
    
    FP.3.bic <- plot_fp_histogram(FP_list.3.bic,
                                  title = "False Positives. 200 tips, 10D, BIC: 2, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.3.bic
    
    
    simresults_200tip_10D_BIC10 <- run_FP_shift_inference(
      base_tree = supertree,
      n_datasets        = 100,
      n_traits          = 10,
      tree_tip_count    = 200,
      num_cores         = 8,
      search_options    = list(
        formula = "trait_data~1",
        min_descendant_tips = 10,
        ic_uncertainty_threshold = 10,
        shift_acceptance_threshold = 10,
        uncertainty = FALSE,
        uncertaintyweights = FALSE,
        uncertaintyweights_par = FALSE,
        postorder_traversal = FALSE,
        plot = FALSE,
        IC = "BIC",
        store_model_fit_history = FALSE,
        num_cores = 1,
        method = 'LL'
      ),
      simulation_options = list(
        variance_mean   = 0.00025,
        variance_sd     = 0.00014,
        covariance_mean = 2.25e-05,
        covariance_sd   = 0.00012
      ),
      seed = 5
    )
    saveRDS(simresults_200tip_10D_BIC10, file='simresults_200tip_10D_BIC10.RDS')
    simresults_200tip_10D_BIC10<-readRDS(file='simresults_200tip_10D_BIC10.RDS')
    
    #what is the FP rate?
    FP_list.4.bic<-list()
    for(i in 1:length(simresults_200tip_10D_BIC10)){
      FP_list.4.bic[i]<-length(simresults_200tip_10D_BIC10[[i]]$shift_nodes_no_uncertainty)/simresults_200tip_10D_BIC10[[i]]$num_candidates
    }
    mean(unlist(FP_list.4.bic))
    
    FP.4.bic <- plot_fp_histogram(FP_list.4.bic,
                                  title = "False Positives. 200 tips, 10D, BIC: 10, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.4.bic
    
  }
  
  #simulation pipeline, 300 tips
  {
    
    simresults_300tip_10D_BIC2 <- run_FP_shift_inference(
      base_tree = supertree,
      n_datasets        = 100,
      n_traits          = 10,
      tree_tip_count    = 300,
      num_cores         = 8,
      search_options    = list(
        formula = "trait_data~1",
        min_descendant_tips = 10,
        ic_uncertainty_threshold = 10,
        shift_acceptance_threshold = 2,
        uncertainty = FALSE,
        uncertaintyweights = FALSE,
        uncertaintyweights_par = FALSE,
        postorder_traversal = FALSE,
        plot = FALSE,
        IC = "BIC",
        store_model_fit_history = FALSE,
        num_cores = 1,
        method = 'LL'
      ),
      simulation_options = list(
        variance_mean   = 0.00025,
        variance_sd     = 0.00014,
        covariance_mean = 2.25e-05,
        covariance_sd   = 0.00012
      ),
      seed = 5
    )
    saveRDS(simresults_300tip_10D_BIC2, file="simresults_300tip_10D_BIC2.RDS")
    simresults_300tip_10D_BIC2<-readRDS(file='simresults_300tip_10D_BIC2.RDS')
    
    #what is the FP rate?
    FP_list.5.bic<-list()
    for(i in 1:length(simresults_300tip_10D_BIC2)){
      FP_list.5.bic[i]<-length(simresults_300tip_10D_BIC2[[i]]$shift_nodes_no_uncertainty)/simresults_300tip_10D_BIC2[[i]]$num_candidates
    }
    mean(unlist(FP_list.5.bic))
    
    FP.5.bic <- plot_fp_histogram(FP_list.5.bic,
                                  title = "False Positives. 300 tips, 10D, BIC: 2, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.5.bic
    
    
    simresults_300tip_10D_BIC10 <- run_FP_shift_inference(
      base_tree = supertree,
      n_datasets        = 100,
      n_traits          = 10,
      tree_tip_count    = 300,
      num_cores         = 8,
      search_options    = list(
        formula = "trait_data~1",
        min_descendant_tips = 10,
        ic_uncertainty_threshold = 10,
        shift_acceptance_threshold = 10,
        uncertainty = FALSE,
        uncertaintyweights = FALSE,
        uncertaintyweights_par = FALSE,
        postorder_traversal = FALSE,
        plot = FALSE,
        IC = "BIC",
        store_model_fit_history = FALSE,
        num_cores = 1,
        method = 'LL'
      ),
      simulation_options = list(
        variance_mean   = 0.00025,
        variance_sd     = 0.00014,
        covariance_mean = 2.25e-05,
        covariance_sd   = 0.00012
      ),
      seed = 5
    )
    saveRDS(simresults_300tip_10D_BIC10, file="simresults_300tip_10D_BIC10.RDS")
    simresults_300tip_10D_BIC10 <- readRDS(file='simresults_300tip_10D_BIC10.RDS')
    
    #what is the FP rate?
    FP_list.6.bic<-list()
    for(i in 1:length(simresults_300tip_10D_BIC10)){
      FP_list.6.bic[i]<-length(simresults_300tip_10D_BIC10[[i]]$shift_nodes_no_uncertainty)/simresults_300tip_10D_BIC10[[i]]$num_candidates
    }
    mean(unlist(FP_list.6.bic))
    
    FP.6.bic <- plot_fp_histogram(FP_list.6.bic,
                                  title = "False Positives. 300 tips, 10D, BIC: 10, min: 10",
                                  xlab = "FP rate",
                                  ylab = "Tree resamples")
    
    FP.6.bic
    
    
    
  }
  
  # Start quartz PDF device
  quartz(width = 10, height = 10, type = 'pdf', file = "FP_histograms.BIC.pdf")
  
  # Generate patchwork layout with labels
  (
    (FP.1.bic | FP.2.bic) /
      (FP.3.bic | FP.4.bic) /
      (FP.5.bic | FP.6.bic)
  ) +
    plot_annotation(tag_levels = "A")
  
  # Close device
  dev.off()
  
  
  # ##########
  # #tring EB -- not complete
  # {
  # 
  # set.seed(seed=1)
  # simdata<-simulate_traits_EB(n_species = 250, n_traits = 2, beta = -5)
  # 
  # sim_search1<-searchOptimalConfiguration(
  #   baseline_tree = paintSubTree(
  #     reorder(as.phylo(simdata$tree), order = 'postorder'),
  #     node = length(simdata$tree$tip.label) + 1,
  #     state = 0
  #   ),
  #   trait_data = simdata$data,
  #   formula = 'trait_data~1',
  #   min_descendant_tips = 5,
  #   num_cores = 3,
  #   ic_uncertainty_threshold = 10,
  #   shift_acceptance_threshold = 5,
  #   uncertainty = F,
  #   uncertaintyweights = F,
  #   uncertaintyweights_par = F,
  #   postorder_traversal = F,
  #   plot = T,
  #   IC = 'GIC')
  # 
  # plotSimmap(
  #   sim_search1$tree_no_uncertainty_untransformed,
  #   ftype='off',
  #   direction='upwards',
  #   colors = generateViridisColorScale.fine(log(sim_search1$model_no_uncertainty$param))$NamedColors
  # )
  # generateLegendLabel(min30.ic40.gic, inset_y = 0.1)
  # 
  # 
  # ###
  # set.seed(seed=1)
  # simdata<-simulate_traits_EB(n_species = 200, n_traits = 10, beta = -0.5)
  # 
  # sim_search1<-searchOptimalConfiguration(
  #   baseline_tree = paintSubTree(
  #     reorder(as.phylo(simdata$tree), order = 'postorder'),
  #     node = length(simdata$tree$tip.label) + 1,
  #     state = 0
  #   ),
  #   trait_data = simdata$data,
  #   formula = 'trait_data~1',
  #   min_descendant_tips = 10,
  #   num_cores = 3,
  #   ic_uncertainty_threshold = 10,
  #   shift_acceptance_threshold = 1,
  #   uncertainty = F,
  #   uncertaintyweights = F,
  #   uncertaintyweights_par = F,
  #   postorder_traversal = F,
  #   plot = T,
  #   IC = 'GIC')
  # 
  # plotSimmap(
  #   sim_search1$tree_no_uncertainty_untransformed,
  #   ftype='off',
  #   direction='upwards',
  #   colors = generateViridisColorScale(sim_search1$model_no_uncertainty$param)$NamedColors
  # )
  # generateLegendLabel(min30.ic40.gic, inset_y = 0.1)
  # }
  #   
}

#false negatives
#proportional scaling
{
  #5 shifts
  {
    #output_FN.1.gic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: GIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.1.gic <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 250,
      #   num_cores       = 8,
      #   seed            = 5,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 5,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE, 
      #     scaleCovariance = F
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "GIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.1.gic, file='output_FN.1.gic.RDS')
      output_FN.1.gic <- readRDS('output_FN.1.gic.RDS')
      
      metrics.output_FN.1.gic <- evaluate_shift_recovery(
        simdata    = output_FN.1.gic$simdata,
        simresults = output_FN.1.gic$results,
        fuzzyDist  = 2, weighted = T
      )
      
      
    }
    #output_FN.1.bic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: BIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.2.bic <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 250,
      #   num_cores       = 8,
      #   seed            = 5,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 5,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE  # Optional visualization of simmap
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "BIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.2.bic, file='output_FN.2.bic.RDS')
      output_FN.2.bic<-readRDS(file='output_FN.2.bic.RDS')
      
      metrics.output_FN.2.bic <- evaluate_shift_recovery(
        simdata    = output_FN.2.bic$simdata,
        simresults = output_FN.2.bic$results,
        fuzzyDist  = 2, weighted = T
      )
      
    }
    
    
  }
  
  #10 shifts
  {
    #output_FN.1.gic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: GIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.3.gic <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 350,
      #   num_cores       = 8,
      #   seed            = 1,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 10,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE  # Optional visualization of simmap
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "GIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.3.gic, file='output_FN.3.gic.RDS')
      output_FN.3.gic <- readRDS(file='output_FN.3.gic.RDS')
      
      metrics.output_FN.3.gic <- evaluate_shift_recovery(
        simdata    = output_FN.3.gic$simdata,
        simresults = output_FN.3.gic$results,
        fuzzyDist  = 2, weighted = T
      )
      
    }
    
    #output_FN.1.bic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: BIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.4.bic <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 350,
      #   num_cores       = 8,
      #   seed            = 1,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 10,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE  # Optional visualization of simmap
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "BIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.4.bic, file='output_FN.4.bic.RDS')
      output_FN.4.bic <- readRDS(file='output_FN.4.bic.RDS')
      
      metrics.output_FN.4.bic <- evaluate_shift_recovery(
        simdata    = output_FN.4.bic$simdata,
        simresults = output_FN.4.bic$results,
        fuzzyDist  = 2, weighted = T
      )
      
    }
    
    
  }
  
  make_performance_table(caption = 'Shift Detection Metrics - Proportional Scaling',
                         list(
                           "5 Shifts (250 tips) — GIC"  = metrics.output_FN.1.gic,
                           "5 Shifts (250 tips) — BIC"  = metrics.output_FN.2.bic,
                           "10 Shifts (350 tips) — GIC" = metrics.output_FN.3.gic,
                           "10 Shifts (350 tips) — BIC" = metrics.output_FN.4.bic
                         )
  )
  
}

#false negatives
#correlation scaling
{
  #5 shifts
  {
    #output_FN.1.gic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: GIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.1.gic.corr <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 250,
      #   num_cores       = 8,
      #   seed            = 5,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 5,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE, 
      #     scaleCovariance = T
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "GIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.1.gic.corr, file='output_FN.1.gic.corr.RDS')
      output_FN.1.gic.corr<-readRDS(file='output_FN.1.gic.corr.RDS')
      
      metrics.output_FN.1.gic.corr <- evaluate_shift_recovery(
        simdata    = output_FN.1.gic.corr$simdata,
        simresults = output_FN.1.gic.corr$results,
        fuzzyDist  = 2, weighted = T
      )
      
      
    }
    
    #output_FN.1.bic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: BIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.2.bic.corr <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 250,
      #   num_cores       = 8,
      #   seed            = 5,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 5,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE, 
      #     scaleCovariance = T
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "BIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.2.bic.corr, file='output_FN.2.bic.corr.RDS')
      output_FN.2.bic.corr<-readRDS(file='output_FN.2.bic.corr.RDS')
      
      metrics.output_FN.2.bic.corr <- evaluate_shift_recovery(
        simdata    = output_FN.2.bic.corr$simdata,
        simresults = output_FN.2.bic.corr$results,
        fuzzyDist  = 2, weighted = T
      )
      
    }
    
    
  }
  
  #10 shifts
  {
    #output_FN.1.gic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: GIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.3.gic.corr <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 350,
      #   num_cores       = 8,
      #   seed            = 1,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 10,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE, 
      #     scaleCovariance = T
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "GIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.3.gic.corr, file='output_FN.3.gic.corr.RDS')
      output_FN.3.gic.corr <- readRDS(file='output_FN.3.gic.corr.RDS')
      
      metrics.output_FN.3.gic.corr <- evaluate_shift_recovery(
        simdata    = output_FN.3.gic.corr$simdata,
        simresults = output_FN.3.gic.corr$results,
        fuzzyDist  = 2, weighted = T
      )
      
    }
    
    #output_FN.1.bic. simulation: 100 replciates, 250 tips. min10-20. emp params. 5 shifts. 
    #inference: BIC, min 10, threshold of 10
    {
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # 1. Call the FN simulation wrapper with updated search settings
      # –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      # output_FN.4.bic.corr <- run_FN_shift_inference(
      #   base_tree = supertree,
      #   n_datasets      = 100,
      #   tree_tip_count  = 350,
      #   num_cores       = 5,
      #   seed            = 1,
      #   simulation_options = list(
      #     minCladeSize    = 10,
      #     maxCladeSize    = 20,
      #     dimensions      = 10,
      #     variance        = 0.00025,
      #     variance_sd     = 0.00014,
      #     covariance      = 2.25e-05,
      #     covariance_sd   = 0.00012,
      #     scaleFactor     = "sample",
      #     scaleFactorRange = c(0.1, 2.0),
      #     excludeRange     = c(0.5, 1.5),
      #     numShifts       = 10,
      #     buffer          = 3,
      #     maxAttempts     = 100,
      #     plot            = FALSE, 
      #     scaleCovariance = T
      #   ),
      #   search_options = list(
      #     formula                    = "trait_data ~ 1",
      #     min_descendant_tips        = 10,
      #     ic_uncertainty_threshold   = 10,
      #     shift_acceptance_threshold = 10,
      #     uncertainty                = FALSE,             # Matches your previous script
      #     uncertaintyweights         = FALSE,
      #     uncertaintyweights_par     = TRUE,              # <- Key difference
      #     postorder_traversal        = FALSE,
      #     plot                       = TRUE,
      #     IC                         = "BIC",
      #     store_model_fit_history    = FALSE,
      #     num_cores                  = 1,
      #     method                     = "LL"
      #   )
      # )
      # saveRDS(output_FN.4.bic.corr, file='output_FN.4.bic.corr.RDS')
      output_FN.4.bic.corr <- readRDS(file='output_FN.4.bic.corr.RDS')
      
      metrics.output_FN.4.bic.corr <- evaluate_shift_recovery(
        simdata    = output_FN.4.bic.corr$simdata,
        simresults = output_FN.4.bic.corr$results,
        fuzzyDist  = 2, weighted = T
      )
      
    }
    
    
  }
  
  make_performance_table(caption = 'Shift Detection Metrics - Correlation Scaling',
                         list(
                           "5 Shifts (250 tips) — GIC"  = metrics.output_FN.1.gic.corr,
                           "5 Shifts (250 tips) — BIC"  = metrics.output_FN.2.bic.corr,
                           "10 Shifts (350 tips) — GIC" = metrics.output_FN.3.gic.corr,
                           "10 Shifts (350 tips) — BIC" = metrics.output_FN.4.bic.corr
                         ), increment_counter = 1
  )
  
  
  
  
  
}
}

#create a plot of all the runs (exploratory, not in manuscript)
{
  #pdf(file = "runs.pdf", height = 6, width = 6, encoding = "UTF-8")
  quartz(type = "pdf", file = "runs.pdf", height = 16, width = 20)
  # Create a layout matrix for a 6x2 grid, filled column-wise
  layout_matrix <- matrix(c(1,2,3,4,5,6,7,8, 9, 10, 11, 12), nrow=6, ncol=2, byrow=FALSE)
  # Apply the layout
  layout(layout_matrix)
  #GIC column
  {
    
    plotSimmap(
      min30.ic40.gic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction='upwards',
      colors = generateViridisColorScale(min30.ic40.gic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min30.ic40.gic, inset_y = 0.1)
    
    plotSimmap(
      min30.ic20.gic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction='upwards',
      colors = generateViridisColorScale(min30.ic20.gic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min30.ic20.gic, inset_y = 0.1)
    
    plotSimmap(
      min20.ic40.gic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction='upwards',
      colors = generateViridisColorScale(min20.ic40.gic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min20.ic40.gic, inset_y = 0.1)
    
    
    plotSimmap(
      min20.ic20.gic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction='upwards',
      colors = generateViridisColorScale(log(min20.ic20.gic$model_no_uncertainty$param))$NamedColors
    )
    generateLegendLabel(min20.ic20.gic, inset_y = 0.1)
    
    
    plotSimmap(
      min10.ic40.gic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction = 'upwards',
      colors = generateViridisColorScale(min10.ic40.gic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min10.ic40.gic, inset_y = 0.1)
    
    #quartz(type = "pdf", file = "M10IC20.pdf", height = 6, width = 30)
    plotSimmap(
      min10.ic20.gic$tree_no_uncertainty_untransformed,
      direction = 'upwards', ftype='off', #, fsize=0.09, offset=2, lwd=1,
      colors = generateViridisColorScale(min10.ic20.gic$model_no_uncertainty$param)$NamedColors, lwd=0.5
    )
    generateLegendLabel(min10.ic20.gic, inset_y = 0.1)
    #dev.off()
    #nodelabels(node = as.numeric(names(getStates(min10.ic20.gic$tree_no_uncertainty_untransformed, type='nodes'))), text = getStates(min10.ic20.gic$tree_no_uncertainty_untransformed, type='nodes'), frame='none', cex=0.5)
    
  }
  #BIC column
  {
    plotSimmap(
      min30.ic40.bic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction='upwards',
      colors = generateViridisColorScale(min30.ic40.bic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min30.ic40.bic, inset_y = 0.1)
    
    plotSimmap(
      min30.ic20.bic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction='upwards',
      colors = generateViridisColorScale(min30.ic20.bic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min30.ic20.bic, inset_y = 0.1)
    
    
    plotSimmap(
      min20.ic40.bic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction = 'upwards',
      colors = generateViridisColorScale(min20.ic40.bic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min20.ic40.bic, inset_y = 0.1)
    
    
    plotSimmap(
      min20.ic20.bic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction = 'upwards',
      colors = generateViridisColorScale(min20.ic20.bic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min20.ic20.bic, inset_y = 0.1)
    
    
    plotSimmap(
      min10.ic40.bic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction = 'upwards',
      colors = generateViridisColorScale(min10.ic40.bic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min10.ic40.bic, inset_y = 0.1)
    
    
    plotSimmap(
      min10.ic20.bic$tree_no_uncertainty_untransformed,
      ftype='off',
      direction = 'upwards',
      colors = generateViridisColorScale(min10.ic20.bic$model_no_uncertainty$param)$NamedColors
    )
    generateLegendLabel(min10.ic20.bic, inset_y = 0.1)
    
  }
  dev.off()
}

#plot of acceptance curves individually
{
  quartz(type = "pdf", file = "runsIC.pdf", height = 16, width = 8)
  # Create a layout matrix for a 6x2 grid, filled column-wise
  layout_matrix <- matrix(c(1,2,3,4,5,6,7,8, 9, 10, 11, 12), nrow=6, ncol=2, byrow=FALSE)
  # Apply the layout
  layout(layout_matrix)
  #GIC column
  {
    plot_ic_acceptance_matrix(
      matrix_data = min30.ic40.gic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min30.ic40.gic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min30.ic20.gic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min30.ic20.gic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min20.ic40.gic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min20.ic40.gic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min20.ic20.gic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min20.ic20.gic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min10.ic40.gic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min10.ic40.gic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min10.ic20.gic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min10.ic20.gic Acceptance Plot"
    )
  }
  #BIC column
  {
    plot_ic_acceptance_matrix(
      matrix_data = min30.ic40.bic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min30.ic40.bic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min30.ic20.bic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min30.ic20.bic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min20.ic40.bic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min20.ic40.bic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min20.ic20.bic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min20.ic20.bic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min10.ic40.bic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min10.ic40.bic Acceptance Plot"
    )
    plot_ic_acceptance_matrix(
      matrix_data = min10.ic20.bic$model_fit_history$ic_acceptance_matrix,
      plot_title = "min10.ic20.bic Acceptance Plot"
    )
  }
  dev.off()
}

#plotting the curves on top of each other
{
  quartz(type = "pdf", file = "IC_overlay.pdf", height = 5, width = 10)
  par(mfrow=c(1,2))
  {
    # Plot the curves together GIC
    plot_multiple_ic_acceptance_curves(
      list(
        "min30.ic40.gic" = min30.ic40.gic,
        "min30.ic20.gic" = min30.ic20.gic,
        "min20.ic40.gic" = min20.ic40.gic,
        "min20.ic20.gic" = min20.ic20.gic,
        "min10.ic40.gic" = min10.ic40.gic,
        "min10.ic20.gic" = min10.ic20.gic
      ),
      plot_title = "GIC Decline Comparison Across Models"
    )
    
    # Plot the curves together BIC
    plot_multiple_ic_acceptance_curves(
      list(
        "min30.ic40.bic" = min30.ic40.bic,
        "min30.ic20.bic" = min30.ic20.bic,
        "min20.ic40.bic" = min20.ic40.bic,
        "min20.ic20.bic" = min20.ic20.bic,
        "min10.ic40.bic" = min10.ic40.bic,
        "min10.ic20.bic" = min10.ic20.bic
      ),
      plot_title = "BIC Decline Comparison Across Models"
    )
  }
  dev.off()
}

#store all of the runs in one object for downstream processing
{
runs<-list(min30.ic40.gic=min30.ic40.gic, 
           min30.ic20.gic=min30.ic20.gic, 
           min20.ic40.gic=min20.ic40.gic, 
           min20.ic20.gic=min20.ic20.gic, 
           min10.ic40.gic=min10.ic40.gic, 
           min10.ic20.gic=min10.ic20.gic, 
           min30.ic40.bic=min30.ic40.bic, 
           min30.ic20.bic=min30.ic20.bic, 
           min20.ic40.bic=min20.ic40.bic, 
           min20.ic20.bic=min20.ic20.bic, 
           min10.ic40.bic=min10.ic40.bic, 
           min10.ic20.bic=min10.ic20.bic)

#store the runs just GIC
runs.gic<-list(min30.ic40.gic=min30.ic40.gic,
               min30.ic20.gic=min30.ic20.gic,
               min20.ic40.gic=min20.ic40.gic,
               min20.ic20.gic=min20.ic20.gic,
               min10.ic40.gic=min10.ic40.gic,
               min10.ic20.gic=min10.ic20.gic)

#store the runs just BIC
runs.bic<-list(
  min30.ic40.bic=min30.ic40.bic,
  min30.ic20.bic=min30.ic20.bic,
  min20.ic40.bic=min20.ic40.bic,
  min20.ic20.bic=min20.ic20.bic,
  min10.ic40.bic=min10.ic40.bic,
  min10.ic20.bic=min10.ic20.bic)
}

#plotting search behavior for supp
{
quartz(file = 'search_behavior.pdf', width = 4*1.5, height = 6*1.5, type = 'pdf')
#par(mfrow = c(6, 2), mar = c(2, 2, 2, 2))  # Fill row-wise: [GIC left, BIC right] pairs

par(mfrow = c(2, 1))  # Fill row-wise: [GIC left, BIC right] pairs

#plot_ic_acceptance_matrix(runs$min30.ic40.gic$model_fit_history$ic_acceptance_matrix, plot_title = 'min30.ic40.gic Search Behavior')
plot_ic_acceptance_matrix(runs$min30.ic40.bic$model_fit_history$ic_acceptance_matrix, plot_title = 'min30.ic40.bic Search Behavior')
#plot_ic_acceptance_matrix(runs$min30.ic20.gic$model_fit_history$ic_acceptance_matrix, plot_title = 'min30.ic20.gic Search Behavior')
#plot_ic_acceptance_matrix(runs$min30.ic20.bic$model_fit_history$ic_acceptance_matrix, plot_title = 'min30.ic20.bic Search Behavior')
#plot_ic_acceptance_matrix(runs$min20.ic40.gic$model_fit_history$ic_acceptance_matrix, plot_title = 'min20.ic40.gic Search Behavior')
#plot_ic_acceptance_matrix(runs$min20.ic40.bic$model_fit_history$ic_acceptance_matrix, plot_title = 'min20.ic40.bic Search Behavior')
#plot_ic_acceptance_matrix(runs$min20.ic20.gic$model_fit_history$ic_acceptance_matrix, plot_title = 'min20.ic20.gic Search Behavior')
#plot_ic_acceptance_matrix(runs$min20.ic20.bic$model_fit_history$ic_acceptance_matrix, plot_title = 'min20.ic20.bic Search Behavior')
#plot_ic_acceptance_matrix(runs$min10.ic40.gic$model_fit_history$ic_acceptance_matrix, plot_title = 'min10.ic40.gic Search Behavior')
#plot_ic_acceptance_matrix(runs$min10.ic40.bic$model_fit_history$ic_acceptance_matrix, plot_title = 'min10.ic40.bic Search Behavior')
plot_ic_acceptance_matrix(runs$min10.ic20.gic$model_fit_history$ic_acceptance_matrix, plot_title = 'min10.ic20.gic Search Behavior')
#plot_ic_acceptance_matrix(runs$min10.ic20.bic$model_fit_history$ic_acceptance_matrix, plot_title = 'min10.ic20.bic Search Behavior')

dev.off()
}

#create plots of increases/decreases in rate by frequency
{
run_deltas<-lapply(runs, rateSummary)
relative_frequencies <- do.call(rbind, lapply(run_deltas, calculateRelativeFrequencies))
relative_counts <- do.call(rbind, lapply(run_deltas, calculateRelativeCounts))
boxplot(relative_frequencies)

analyze_paired_counts(as.data.frame(relative_counts))

quartz(file='run_rates.pdf', width=12,height=4, type='pdf')
create_raincloud_with_wilcox(relative_frequencies, horizontal_spacing = 0.75, title = 'Shift frequencies (GIC+BIC)', 
                             group_colors = c("decrease" = "#e41a1c", "increase" = "#377eb8"),
                             y_limits=c(0,1)) + geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") |
  create_raincloud_with_wilcox(relative_frequencies[!grepl("bic", rownames(relative_frequencies)), ], horizontal_spacing = 0.75, title = 'Shift frequencies (GIC)', 
                               group_colors = c("decrease" = "#e41a1c", "increase" = "#377eb8"),
                               y_limits=c(0,1)) + geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") |
  create_raincloud_with_wilcox(relative_frequencies[!grepl("gic", rownames(relative_frequencies)), ], horizontal_spacing = 0.75, title = 'Shift frequencies (BIC)', 
                               group_colors = c("decrease" = "#e41a1c", "increase" = "#377eb8"),
                               y_limits=c(0,1)) + geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") 
dev.off()
}

#process runs
{
processruns<-processTransitionRates(runs)
processruns.gic<-processTransitionRates(runs.gic)
processruns.bic<-processTransitionRates(runs.bic)

rowMeans(rbind(sapply(processruns, function(df) min(abs(df$percentage_change), na.rm = TRUE)),
               sapply(processruns, function(df) median(abs(df$percentage_change), na.rm = TRUE)),
               sapply(processruns, function(df) max(abs(df$percentage_change), na.rm = TRUE))))
}

#estimate modal increase/decrease
{
  
  # Function to estimate mode using kernel density
  density_mode <- function(x) {
    d <- density(x, na.rm = TRUE)
    d$x[which.max(d$y)]
  }
  
  # Merge all data frames in processruns
  all_changes <- do.call(rbind, processruns)
  
  # Estimate mode from density for increases and decreases
  mode_increase <- density_mode(all_changes$percentage_change[all_changes$rate_change == "increase"])
  mode_decrease <- density_mode(all_changes$percentage_change[all_changes$rate_change == "decrease"])
  
  # Print results
  cat("Density Mode % Change (Increase):", round(mode_increase, 2), "%\n")
  cat("Density Mode % Change (Decrease):", round(mode_decrease, 2), "%\n")
  
  
  bootstrap_mode_test <- function(x, y, n_boot = 1000, bw = NULL, plot = FALSE) {
    # KDE mode estimator
    estimate_mode <- function(z) {
      d <- density(z, na.rm = TRUE)
      d$x[which.max(d$y)]
    }
    
    # Smoothed resampler
    smoothed_sample <- function(data, n, bw) {
      sample(data, n, replace = TRUE) + rnorm(n, 0, sd = bw)
    }
    
    # Mode difference
    m_x <- estimate_mode(x)
    m_y <- estimate_mode(y)
    d_obs <- abs(m_x - m_y)
    
    # Pooled sample
    pooled <- c(x, y)
    n_x <- length(x)
    n_y <- length(y)
    bw <- bw %||% bw.nrd0(pooled)
    
    # Naive bootstrap
    d_null_naive <- replicate(n_boot, {
      resample <- sample(pooled, n_x + n_y, replace = TRUE)
      abs(estimate_mode(resample[1:n_x]) - estimate_mode(resample[(n_x + 1):(n_x + n_y)]))
    })
    
    # Smoothed bootstrap
    d_null_smoothed <- replicate(n_boot, {
      x_star <- smoothed_sample(pooled, n_x, bw)
      y_star <- smoothed_sample(pooled, n_y, bw)
      abs(estimate_mode(x_star) - estimate_mode(y_star))
    })
    
    # P-values
    p_naive <- mean(d_null_naive >= d_obs)
    p_smoothed <- mean(d_null_smoothed >= d_obs)
    
    # Output object
    result <- list(
      observed_mode_difference = d_obs,
      naive = list(
        p_value = p_naive,
        null_distribution = d_null_naive
      ),
      smoothed = list(
        p_value = p_smoothed,
        null_distribution = d_null_smoothed,
        bandwidth = bw
      )
    )
    
    # Summary print
    cat("Bootstrap Mode Difference Test (BMDT)\n")
    cat("Observed mode difference:", round(d_obs, 3), "\n")
    cat("Naive bootstrap p-value  :", signif(p_naive, 3), "\n")
    cat("Smoothed bootstrap p-value:", signif(p_smoothed, 3), "\n")
    cat("Bandwidth used for smoothing:", signif(bw, 3), "\n")
    
    # Plotting
    if (plot) {
      op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
      hist(d_null_smoothed, breaks = 30, col = "lightblue",
           main = "Smoothed Null Distribution", xlab = "Mode Difference")
      abline(v = d_obs, col = "red", lwd = 2)
      legend("topright", legend = paste("Observed =", round(d_obs, 2),
                                        "\nP =", signif(p_smoothed, 3)),
             bty = "n")
      
      hist(d_null_naive, breaks = 30, col = "lightgray",
           main = "Naive Null Distribution", xlab = "Mode Difference")
      abline(v = d_obs, col = "red", lwd = 2)
      legend("topright", legend = paste("Observed =", round(d_obs, 2),
                                        "\nP =", signif(p_naive, 3)),
             bty = "n")
      par(op)
    }
    
    invisible(result)
  }
  
  bootstrap_mode_results<-bootstrap_mode_test(log(abs(all_changes$rate_delta[all_changes$rate_change == "decrease"])), 
                      log(abs(all_changes$rate_delta[all_changes$rate_change == "increase"])), 
                      n_boot=10000, plot = T)
  
}

#plot rate differentials
{
quartz(file='rateDeltas.pdf', type='pdf', width=12, height=4)
par(mfrow=c(1,3))
plotRateChangeDensities(processruns, use_log = T, text_y=0.35, ylim=c(0, 0.4), xlim=c(-10.5,-5), bty='n', test_type = 'ks', y_nudge = 8, x_nudge = 0.5, legend_x_nudge = -0.021, title = "Rate Δ (GIC + BIC)")
plotRateChangeDensities(processruns.gic, use_log = T, text_y=0.35, ylim=c(0, 0.4), xlim=c(-10.5,-5), bty='n', test_type = 'ks', y_nudge = 8, x_nudge = 0.5, legend_x_nudge = -0.021, title = "Rate Δ (GIC)")
plotRateChangeDensities(processruns.bic, use_log = T, text_y=0.35, ylim=c(0, 0.4), xlim=c(-10.5,-5), bty='n', test_type = 'ks', y_nudge = 8, x_nudge = 0.5, legend_x_nudge = -0.021, title = "Rate Δ (BIC)")
dev.off()
}

#plotting for Figure 2/ westerhold data
#westerhold data
{
  
  westerhold<-read_xlsx(path="/Users/cotinga/jacob.berv@gmail.com/Code/passerine-bodyplan-evolution/other_data/aba6853_tables_s8_s34.xlsx", skip = 1, sheet = 'Table S33')
  westerhold<-westerhold[complete.cases(westerhold$age_tuned), ]
  westerhold.long<-read_xlsx(path="/Users/cotinga/jacob.berv@gmail.com/Code/passerine-bodyplan-evolution/other_data/aba6853_tables_s8_s34.xlsx", skip = 1, sheet = 'Table S34', trim_ws = TRUE)
  westerhold.long<-westerhold.long[complete.cases(westerhold.long$age_tuned), ]
  
  #new plot Figure 2
  {
    pdf(file='climate_plot-tmp.pdf', height=(6*1.5)*0.9, width=(6*1.4)*0.9)
    
    #d18 plot
    {
      
      #set up plotting area (but dont plot the data)
      # Call the function with all customizations
      label_plot_with_annotations(
        data = westerhold, 
        data_long = westerhold.long, 
        x_coords = NULL, 
        x_col = "age_tuned", 
        y_col = "benthic d18O VPDB CorrAdjusted", 
        y_col_long = "ISOBENd18oLOESSsmoothLongTerm", 
        plot_cex = 0.1, 
        plot_type = 'b', 
        plot_xlim = c(45.24011, 1.74),#1.74), #min(westerhold$age_tuned)
        plot_ylim = c(6, -1.5), 
        point_col = scales::alpha('#377eb8',  0.25), 
        y_buffer = 0.4, 
        label_cex = 0.8,
        line_color = "red",       
        text_color = "black",
        line_width = 1.2,         
        line_type = 1,            
        text_size = 1, 
        title = expression(delta^18 * O * " Temperature Proxy"), 
        xlab = "Ma", 
        ylab = expression(delta^18 * O),
        suppress_x_axis = T, plot_data = F
      )
      
      # 1. Get all epochs except Holocene
      epochs <- GTS2020[GTS2020$rank == "epoch" & GTS2020$interval_name != "Holocene", ]
      
      # 3. Custom abbreviations for Pliocene and Pleistocene
      epochs$abbr <- ifelse(
        epochs$interval_name %in% c("Pliocene", "Pleistocene"),
        c("Plio.", "Plei."),
        epochs$interval_name
      )
      
      axis_geo(
        side = 1,
        intervals = epochs,
        height = 0.05,
        fill = c("white", "grey90", "white", "grey90", "white"),
        lab = TRUE,
        lab_size = 0.7,
        tick_labels = TRUE,
        bord_col = "black", lwd=list(1),
        abbr = TRUE,  # this now uses the customized `abbr` column
        skip = NULL, exact=F, round=1, neg = F, bkgd = NA
      )
      
      usr <- par("usr")
      x_min <- 0; x_max <- 46
      minor_ticks <- seq(10*floor(x_min/10), 10*ceiling(x_max/10), by = 1)
      minor_ticks <- minor_ticks[minor_ticks >= x_min & minor_ticks <= x_max]  # clamp to geo range
      nudge <- -0.05 * diff(usr[3:4])
      segments(
        x0 = minor_ticks, x1 = minor_ticks,
        y0 = usr[3] + nudge, y1 = usr[3] - 0.0075 * diff(usr[3:4]) + nudge,
        col = "black", lwd = 0.8, xpd = NA
      )
      
      # op <- par(xpd = NA)                  # allow drawing outside plot
      # usr <- par("usr")                    # get current plot dimensions
      # 
      # # Draw outside the normal plot region
      # rect(
      #   xleft = -0,
      #   xright = usr[2],
      #   ybottom = usr[3] - 1 * diff(usr[3:4]),
      #   ytop = usr[4] + 2 * diff(usr[3:4]),
      #   col = "white", border = NA
      # )
      # 
      # # Restore clipping and xpd
      # par(op)                               # reset xpd and other settings
      # do.call("clip", as.list(usr))         # reset clipping to plot area
      # 
      #box()
      #axis_geo(side = 1, intervals = "epoch", tick_labels = T, 
      #        lwd=1, height=0.05, fill = c("white", "grey90", "white", "grey90", "white"))
      
      #par(xpd=NA)
      add_event_interval(
        start = 40.5, end = 40.0,
        label = "Middle Eocene Climate Optimum",
        y = 0,
        interval_height = 0.05,
        label_offset = -0.3,
        wrap_words = 2,
        lwd = 1.5,
        col = "black",
        label_cex = 0.75
      )
      #abline(v=mean(c(40.5,40.0)), ylim=c(2,4))
      
      add_event_interval(
        start = 34.1, end = 33.6,
        label = "Eocene Oligocene Transition",
        y = 0.75,
        interval_height = 0.05,
        label_offset = -0.3,
        wrap_words = 2,
        lwd = 1.5,
        col = "black",
        label_cex = 0.75
      )
      
      add_event_interval(
        start = 17.0, end = 14.7,
        label = "Miocene Climate Optimum",
        y = 1,
        interval_height = 0.05,
        label_offset = -0.3,
        wrap_words = 2,
        lwd = 1.5,
        col = "black",
        label_cex = 0.75
      )
      
      # add_event_interval(
      #   start = 14.2, end = 13.8,
      #   label = "Middle Miocene Climatic Transition",
      #   y = 1,
      #   interval_height = 0.05,
      #   label_offset = -0.3,
      #   wrap_words = 2,
      #   lwd = 1.5,
      #   col = "black",
      #   label_cex = 0.75
      # )
      # 
      
      
      # add_event_interval(
      #   start = 15, end = 13,
      #   label = "Middle Miocene Climate Transition",
      #   y = 2,
      #   interval_height = 0.05,
      #   label_offset = -0.3,
      #   wrap_words = 2,
      #   lwd = 1.5,
      #   col = "black",
      #   label_cex = 0.75
      # )
      
      # add_event_interval(
      #   start = 7.2, end = 5.5,
      #   label = "Late Miocene Cooling",
      #   y = 2.25,
      #   interval_height = 0.05,
      #   label_offset = -0.3,
      #   wrap_words = 2,
      #   lwd = 1.5,
      #   col = "black",
      #   label_cex = 0.75
      # )
      
      # par(xpd=T)
      
      par(new=T)
      # Call the function with all customizations
      label_plot_with_annotations(
        data = westerhold, 
        data_long = westerhold.long, 
        x_coords = NULL, 
        x_col = "age_tuned", 
        y_col = "benthic d18O VPDB CorrAdjusted", 
        y_col_long = "ISOBENd18oLOESSsmoothLongTerm", 
        plot_cex = 0.1, 
        plot_type = 'b', 
        plot_xlim = c(45.24011, 1.74),#1.74), #min(westerhold$age_tuned)
        plot_ylim = c(6, -1.5), 
        point_col = scales::alpha('#377eb8',  0.25), 
        y_buffer = 0.4, 
        label_cex = 0.8,
        line_color = "red",       
        text_color = "black",
        line_width = 1.2,         
        line_type = 1,            
        text_size = 1, 
        title = expression(delta^18 * O * " Temperature Proxy"), 
        xlab = "Ma", 
        ylab = expression(delta^18 * O),
        suppress_x_axis = T
      )
      
    }
    usr <- par("usr"); y0 <- usr[3]-0.05*diff(usr[3:4])
    segments(c(usr[1],usr[2]), y0, c(usr[1],usr[2]), c(0,2), xpd=NA)
    #add minor ticks on left y
    usr <- par("usr")
    minor_ticks <- seq(0, 6, by = diff(pretty(usr[3:4], 5))[1] / 10)
    segments(usr[1], minor_ticks, usr[1] - 0.015 * diff(usr[1:2]), minor_ticks, col = "black", lwd = 0.8, xpd = NA)
    
    #d18 variability
    {
      # tmp<-calculate_variability(westerhold, "age_tuned", "benthic d18O VPDB CorrAdjusted", bin_size = 5, slide_step = 0.1, metric = "acf", lag = 5)
      # plot(tmp$variability ~ tmp$bin_midpoint, xlim=c(45.24011, min(westerhold$age_tuned)), type='b', cex=0.5, col='rosybrown')
      # tmp<-calculate_variability(westerhold, "age_tuned", "benthic d18O VPDB CorrAdjusted", bin_size = 5, slide_step = 0.1, metric = "sd")
      # points(tmp$variability ~ tmp$bin_midpoint, xlim=c(45.24011, min(westerhold$age_tuned)), type='b', cex=0.5, col='skyblue')
      par(new=T)
      plot_with_dual_y_axis(westerhold, 
                            x_col = "age_tuned", 
                            y_col1 = "benthic d18O VPDB CorrAdjusted", 
                            y_col2 = "benthic d18O VPDB CorrAdjusted",
                            xlim=c(45.24011, 1.74), #min(westerhold$age_tuned)), 
                            plot_type = "spline", 
                            line_width1 = 5, 
                            line_width2 = 2,
                            lag = 2, 
                            plot_metric2 = F, 
                            metric1 = 'sd',
                            metric2 = 'sd',
                            col1 = make.transparent(color='darkgrey', alpha=0.75), 
                            slide_step = 0.5,
                            bin_size = 5, 
                            ylab1 = 'SD',
                            ylab2= 'SD',
                            metric1_on_right = T,
                            main = "",
                            xlab = "",
                            suppress_x_axis = T, ylim1 = c(-1.65,0.65), suppress_y_axis = T)
      axis(side=4, at = c(0.2,0.4,0.6), col = 'darkgrey', col.ticks = 'darkgrey', col.axis='darkgrey')
      # Place label in user coords (x = right margin, y = center of axis)
      usr <- par("usr")
      text(x = usr[2] -5,   # push a bit outside plot (adjust as needed)
           y = mean(c(0.2, 0.6)),  # center along the axis range
           labels = "SD", srt = 90, cex = 1, col = "darkgrey", xpd = TRUE)
      
      par(new=T)
      plot_with_dual_y_axis(westerhold, 
                            x_col = "age_tuned", 
                            y_col1 = "benthic d18O VPDB CorrAdjusted", 
                            y_col2 = "benthic d18O VPDB CorrAdjusted",
                            xlim=c(45.24011, 1.74), #min(westerhold$age_tuned)), 
                            plot_type = "spline", 
                            line_width1 = 1, 
                            line_width2 = 1.5, 
                            lag = 2, 
                            plot_metric2 = F, 
                            metric1 = 'sd',
                            metric2 = 'sd',
                            col1 = 'black', 
                            slide_step = 0.5, 
                            bin_size = 5,
                            ylab1 = 'SD',
                            ylab2= 'SD',
                            metric1_on_right = T,
                            main = "",
                            xlab = "",
                            suppress_x_axis = T,
                            suppress_y_axis = T,
                            line_type1 = 2, ylim1 = c(-1.65,0.65)
      )
      
    }
    
    #spline summary
    {
      # splines <- processTreeDataFromObjects.filter(
      #   runs,
      #   bin_size = 5,
      #   slide_step = 0.5,
      #   legend = F,
      #   xlim = c(45.24011, min(westerhold$age_tuned)),
      #   filter_by = NULL,
      #   normalization = 'log',
      #   ylim = c(0, 3.5), central_tendency = 'mean'
      # )
      #abline(v=15, lty=2)
      #abline(v=32.5, lty=2)
      #abline(v=40, lty=2)
      #abline(v=24.5, lty=2)
      # splines <- processTreeDataFromObjects.filter(
      #   runs,
      #   bin_size = 5,
      #   slide_step = 0.5,
      #   legend = F,
      #   xlim = c(45.24011, min(westerhold$age_tuned)),
      #   filter_by = 'increase',
      #   normalization = 'log',
      #   ylim = c(0, 1.5)
      # )
      par(new=T)
      splines.all <- processTreeDataFromObjects.filter(
        runs,
        bin_size = 5,
        slide_step = 0.5,
        legend = T,
        xlim = c(45.24011, 1.74), #min(westerhold$age_tuned)),
        filter_by = 'both',
        normalization = 'log',
        ylim = c(0, 4.5), 
        legend_position = 'topright', 
        hdi_prob = 0.90, 
        central_tendency = 'mean',
        bootstrap_hdi = T, y_axis_side = 'right', show_x_axis = F, show_y_axis = F, title = '')
      axis(side=4, at = c(0,1,2,3), col = 'black', col.ticks = 'black', col.axis='black')
      
      #add_ratio_curve(splines.all, 
      #                xlim = c(45.24011-3, min(westerhold$age_tuned+4)), 
      #                lwd = 4, lty = 1)
      #abline(h=1)
      
    }
    
    dev.off()
    
    
  }
  
}

#calculating branch/shift metrics
{

shift_rate_metrics.shifts <- compute_shift_rates.shifts(min10.ic20.gic$tree_no_uncertainty_untransformed, 
                                                        n_cores = 8, weighting = "ES", 
                                                        verbose = T, decay_base = 2, 
                                                        param = min10.ic20.gic$model_no_uncertainty$param, 
                                                        use_absolute_magnitudes = T, normalize_weights = T, T = 5)

shift_rate_metrics.branches <- compute_shift_rates.branches(min10.ic20.gic$tree_no_uncertainty_untransformed, 
                                                            n_cores = 8, weighting = "ES", 
                                                            verbose = T, decay_base = 2, 
                                                            param = min10.ic20.gic$model_no_uncertainty$param, 
                                                            use_absolute_magnitudes = T, normalize_weights = T, T = 5)

Total_Shift_Counts <- shift_rate_metrics.shifts$Unweighted_Shift_Count
names(Total_Shift_Counts) <- shift_rate_metrics.shifts$Tip

#store the values for lineage rate (shift events focux)
log_lineage_rate.shifts <- log(shift_rate_metrics.shifts$Weighted_Phenotypic_Rate)
names(log_lineage_rate.shifts) <- shift_rate_metrics.shifts$Tip
log_tip_rate.shifts <- log(shift_rate_metrics.shifts$Tip_Phenotype_Rate)
names(log_tip_rate.shifts) <- shift_rate_metrics.shifts$Tip

#store the values for lineage rate (branch events focus)
log_lineage_rate.branches <- log(shift_rate_metrics.branches$Weighted_Phenotypic_Rate)
names(log_lineage_rate.branches) <- shift_rate_metrics.branches$Tip
log_tip_rate.branches <- log(shift_rate_metrics.branches$Tip_Phenotype_Rate)
names(log_tip_rate.shifts) <- shift_rate_metrics.branches$Tip


#check the waiting time calculations
waiting.times.global<-summarize_waiting_times(ladderize(untangle(min10.ic20.gic$tree_no_uncertainty_untransformed)))
waiting.times<-summarize_tree_shift_metrics(min10.ic20.gic, verbose = T)
waiting.times.temporal<-analyze_shift_waiting_times(min10.ic20.gic, verbose = T)
waiting.times.lineage<-analyze_shift_waiting_times_by_branch(min10.ic20.gic, verbose = T)

#add tipDR (as a control) using jon chang's code
source('https://raw.githubusercontent.com/jonchang/fastdivrate/refs/heads/master/R/dr.R')
tipDR <- DR_statistic(min10.ic20.gic$tree_no_uncertainty_untransformed)

sampled_cv_shift_metrics<-cbind(sampled_cv, Shift_Counts = Total_Shift_Counts, log_lineage_rate_branches = log_lineage_rate.branches, log_lineage_rate_shifts = log_lineage_rate.shifts, log_tip_rate = log_tip_rate.shifts, tipDR = tipDR) 

#model residuals (for spatial analysis of disparity)
#shift.residuals<-residuals(min10.ic20.gic$model_no_uncertainty, type = 'response')[sampled_cv$phylo,]
shift.residuals<-residuals(global.model.mvBM, type = 'response')[sampled_cv$phylo,]
sampled_cv_shift_metrics<-cbind(sampled_cv_shift_metrics, residuals = shift.residuals) #this is ok bc names match

#try reordering to match plot order
sampled_cv_shift_metrics<-sampled_cv_shift_metrics[
  jntools::get_tips_in_ape_plot_order(ladderize(untangle(
    min10.ic20.gic$tree_no_uncertainty_untransformed
  ))),]

#save the object to use for spatial analysis script
saveRDS(sampled_cv_shift_metrics, file='sampled_cv_shift_metrics_8_08_25.RDS')
sampled_cv_shift_metrics<-readRDS('sampled_cv_shift_metrics_8_08_25.RDS')
  
#cor(log(shift_rate_metrics.branches$Tip_Phenotype_Rate), 
#log(shift_rate_metrics.branches$Weighted_Phenotypic_Rate))

#cor(log(shift_rate_metrics.shifts$Tip_Phenotype_Rate), 
#    log(shift_rate_metrics.shifts$Weighted_Phenotypic_Rate))


}

#distribution fitting/testing
{
  
  #define the color palette for figure 1 histogram/branch rates
  rate_colors <- assignRateColors(min10.ic20.gic$model_no_uncertainty$param, 
                                  breaksmethod = 'fisher', logcolor = T, 
                                  color.interval = NULL, palette = "RdYlBu", 
                                  nbreaks = 10, reverse = T, jitter_on_failure = T)
  
  
  # Load and log-transform evolutionary rate data
  rates <- log(min10.ic20.gic$model_no_uncertainty$param)
  
  # Fit Gumbel and GEV distributions
  fit_gumbel <- mlgumbel(rates)
  fit_gumbel.boot <- bootstrapml(fit_gumbel)
  gumbel_kl <- compute_kl_bootstrap(fit_gumbel, data = rates, n_boot = 10000)
  
  #dists<-model_select(rates - min(rates) + 1, type = "continuous", return="all")
  dists<-model_select(rates, type = "continuous", return="all")
  
  #identify waiting times
  global.waiting.times<-model_select(waiting.times.temporal$waiting_times$TimeToNext[-1], type = "continuous", return="all", models = c(
    "exp",       # Baseline: memoryless process, simplest null model
    "gamma",     # Flexible shape; handles overdispersion and skew
    "weibull",   # Varying hazard rate; increasing or decreasing failure risk
    "lnorm",     # Skewed positive values; models multiplicative processes well
    "invgauss"   # First-passage time processes (e.g. Brownian motion-based),
  ))
  
  global.lineage.sum.times<-model_select((do.call(rbind, waiting.times.lineage$lineage_waiting_times)$TimeToNext), type = "continuous", return="all", models = c(
    "exp",       # Baseline: memoryless process, simplest null model
    "gamma",     # Flexible shape; handles overdispersion and skew
    "weibull",   # Varying hazard rate; increasing or decreasing failure risk
    "lnorm",     # Skewed positive values; models multiplicative processes well
    "invgauss"   # First-passage time processes (e.g. Brownian motion-based)
  ))
  
  separate.lineage.waiting.times <- fit_lineage_distributions(waiting.times.lineage$lineage_waiting_times, num_cores = 6, 
                                                              models = c(
                                                                "exp",       # Baseline: memoryless process, simplest null model
                                                                "gamma",     # Flexible shape; handles overdispersion and skew
                                                                "weibull",   # Varying hazard rate; increasing or decreasing failure risk
                                                                "lnorm",     # Skewed positive values; models multiplicative processes well
                                                                "invgauss"   # First-passage time processes (e.g. Brownian motion-based)
                                                              ))

  lambda_df <- extract_exp_lambda(separate.lineage.waiting.times)
  summary(lambda_df$Lambda, na.rm=T)
  mean(lambda_df$Lambda, na.rm=T)
  hdi(lambda_df$Lambda, na.rm=T)
  
  
  boot_gumbel_raw <- bootstrapml(
    fit_gumbel, 
    reps = 1000, 
    map = identity,
    reducer = identity
  )
  
  ratesHistogram(
    min10.ic20.gic$model_no_uncertainty$param,
    rate_colors$colors,
    rate_colors$breaks,
    plotBrks = TRUE,
    useDensity = TRUE,
    logscale = TRUE,
    title = "Break method: Fisher",
    xBuffer = 0.02,
    yBuffer = 0.03,
    brksCol = 'darkgrey',
    lwd = 0.5,
    cex.axis = 0.5,
    cex.lab = 0.75,
    cex.main = 0.75,
    x.axis.label.line = 1.25,
    y.axis.label.line = 1.5,
    x.tick.label.line = 0.25,
    y.tick.label.line = 0.55, plotOverlay = 
      {overlay_gumbel_fit(
        fit_gumbel,          # The main MLE fit
        boot_gumbel_raw,     # The entire bootstrap param matrix
        data = rates,        # Optional raw data for density
        addDataDensity = F,
        addLegend =  T,
        col = "red",
        lty = 2,
        lwd = 2, bandFill = alpha('red', 0.1)
      )}
  )
}

#plotting figure 1 stuff
{

  #par(mfrow = c(1,1))
  cairo_pdf(file='fanPlot.pdf', height=11*1.2, width=8.5*1.2)
  #top 1/2
  #grid.newpage()
  {
    library(grid)
    library(ggplot2)
    library(TeachingDemos)
 
    # 1) Plot your arc tree
    plotFanTree.wTraits(
      tree = (min10.ic20.gic$tree_no_uncertainty_untransformed),
      X = cbind(scale(sampled_cv_shift_metrics$log_lineage_rate_branches)),
      colorvec = rate_colors$colors,
      type = 'arc', ftype = 'off', part = 0.5, 
      arc_height = 1.25, 
      color_palette = hcl.colors, trait_scale = 0.15, plot=FALSE
    )
    
    # #plot((sampled_cv_shift_metrics$log_lineage_rate_branches) ~ reconstructed_lat)
    # plotFanTree.wTraits(
    #   tree = (min10.ic20.gic$tree_no_uncertainty_untransformed),
    #   X = cbind(abs(reconstructed_lat)),
    #   colorvec = rate_colors$colors,
    #   type = 'arc', ftype = 'off', part = 0.5, 
    #   arc_height = 1.25, 
    #   color_palette = hcl.colors, trait_scale = 0.15, plot=FALSE
    # )
    # Assuming `tree` is already plotted and `min10.ic20.gic$ic_weights` is loaded
    lowconf<-plotLowConfidenceShiftNodes(tree, min10.ic20.gic$ic_weights, dot_cex = 2, threshold = 0.9, dot_color = alpha('black', alpha = 0.5), pch_select=16)
    
    subplot.gg(
      plot_object = generate_color_bar(
        break_colors = list(
          breaks = seq(
            min((sampled_cv_shift_metrics$log_lineage_rate_branches)),
            max((sampled_cv_shift_metrics$log_lineage_rate_branches)),
            length.out = 101
          ),
          colors = hcl.colors(100, "Viridis")
        ),
        bar_title = "Log Lineage Rate", n_ticks = 4
      ),
      bbox = c(95, 75-10, 110, 115-10),  # these are user coordinates (xmin, ymin, xmax, ymax)
      units = "npc"  # bbox coordinates are converted internally
    )
    
    
    # Call function to add node labels
    keys<-plotNodeLabelsByRateChange(
      tree = min10.ic20.gic$tree_no_uncertainty_untransformed,
      rate_data = runs$min10.ic20.gic,
      filter_by = "increase",  # Can be "increase", "decrease", or "all"
      color_increase = "white", 
      color_decrease = "blue", 
      base_cex = 1.5, 
      scale_factor = 0.08, 
      alpha = 0.75, 
      transform_method = 'square_root', 
      letter_cex = 0.75, plot_letters = T
    )
    (keys$PercentageChange)
    
    merge_keys_weights<-merge(
      keys,
      min10.ic20.gic$ic_weights,
      by.x = "Node",
      by.y = "node",
      all.x = TRUE
    )
    
    with(density(keys$PercentageChange), x[which.max(y)])
    tmp<-extractMaxAgeOfRegimesWithRateChanges(runs$min10.ic20.gic)
    with(density(tmp[tmp$rate_change == "decrease", ]$percentage_change), x[which.max(y)])
    median(tmp[tmp$rate_change == "decrease", ]$percentage_change)
    
    
    # Add small black dots for low-confidence shift nodes
    # low_conf_nodes <- plotLowConfidenceShiftNodes(
    #   tree = min10.ic20.gic$tree_no_uncertainty_untransformed,
    #   ic_weights = min10.ic20.gic$ic_weights,
    #   threshold = 0.9,
    #   dot_cex = 0.5,
    #   dot_color = "red"
    # )
    
    # nodelabels(node=c(3807, 3153, 3924, 3900, 2311, 2073, 2603))
    # low_conf_nodes <- plotLowConfidenceShiftNodes(
    #   tree = min10.ic20.gic$tree_no_uncertainty_untransformed,
    #   ic_weights = min10.ic20.gic$ic_weights,
    #   threshold = 0.9,
    #   dot_cex = 0.5,
    #   dot_color = "red"
    # )
    
    # 2) Inset histogram
    TeachingDemos::subplot(
      ratesHistogram(
        min10.ic20.gic$model_no_uncertainty$param,
        rate_colors$colors,
        rate_colors$breaks,
        plotBrks = TRUE,
        useDensity = TRUE,
        logscale = TRUE,
        title = "Break method: Fisher",
        xBuffer = 0.02,
        yBuffer = 0.03,
        brksCol = 'darkgrey',
        lwd = 0.25,
        cex.axis = 0.5,
        cex.lab = 0.75,
        cex.main = 0.75,
        x.axis.label.line = 1.25,
        y.axis.label.line = 1.5,
        x.tick.label.line = 0.25,
        y.tick.label.line = 0.55, plotOverlay = 
          {overlay_gumbel_fit(
            fit_gumbel,          # The main MLE fit
            boot_gumbel_raw,     # The entire bootstrap param matrix
            data = rates,        # Optional raw data for density
            addDataDensity = F,
            addLegend =  T,
            col = "red",
            lty = 2,
            lwd = 2, bandFill = alpha('red', 0.25), legend_cex = 0.65
          )}
      ),
      x = "center",
      y = "center",
      size = c(2.4, 1.25),
      hadj = 0.35,
      vadj = 1.25
    )
    
  }
  
  dev.off()  # CLOSE PDF DEVICE (ensures only one page)
  
  
  #get clade information (for taxonomy identification)
  {
    keynodes <- extractMaxAgeOfRegimesWithRateChanges(min10.ic20.gic, filter_by = 'increase')$max_age_node[-1]
    perc_changes <- extractMaxAgeOfRegimesWithRateChanges(min10.ic20.gic, filter_by = 'increase')$percentage_change[-1]
    #modal accelaration is 120-166%
    
    clade_info <- extract_clade_info(min10.ic20.gic, keynodes, cores = 4)
    
    extractMaxAgeOfRegimesWithRateChanges(min10.ic20.gic, filter_by = 'increase')
    
    # Given a list of taxa for a focal node and a list of taxa for its sister node from a phylogenetic tree, produce a precise and concise taxonomic assessment that:
    #   1.	Identifies the narrowest named monophyletic group that fully encompasses all taxa within the focal node.
    # 2.	Determines whether the focal node corresponds exactly to a recognized taxonomic rank, or if it represents an incomplete subset due to missing or unsampled members.
    # 3.	Compares the composition of the focal node to its sister node to evaluate whether the focal node forms a distinct, named taxonomic rank or an unranked clade.
    # 4.	Explicitly states whether the focal node is taxonomically equivalent to a known rank and provides justification for any mismatches.
    # 5.	Identifies any taxa missing from the focal node that would be required for it to fully match a recognized taxonomic rank.
    # 6.	Assigns a concise, publication-ready descriptor to the focal node based on its closest taxonomic match (e.g., a formal rank, an informal grouping, or a relevant clade name).
    # 7.	Notes any informal, biogeographic, or alternative classifications that provide additional evolutionary context.
    # 
    # The response should be precise, publication-ready, and supported by the taxonomic structure of the provided lists. We want the best rank identifier considering the phylogenetic context and what is included in the sample. Eg, we may have some missing lineages known to be a member of the rank, but as long as the sister node includes different ranks, we can assign the focal rank to the most appropriate rank that is compatible with our sample. Assume missing lineages in a focal sample (and their absence in the sister sample) are due to incomplete sampling. Use terms like sensu stricto or sensu lato, only if helpful in clarifing your analysis. Lastly, for you final determination, provide the reference for the suggested taxonomic rank, for example, in this format:  "Ohlson et al (2013)" 
    
    #clade 1 - Core Sylvioidea -\n Babbler group -- possible ref? 10.1016/j.ympev.2005.05.015
    {
      clade_info$unique_families[1]$Node_2845
      unique(clade_info$clade_members[1]$Node_2845$genus)
      unique(clade_info$clade_members[1]$Node_2845$verbatim_name)
      clade_info$sister_unique_families[1]$Node_2845
      unique(clade_info$sister_clade_members[1]$Node_2845$genus)
      #unique(clade_info$sister_clade_members[1]$Node_2845$verbatim_name)
      #getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2845, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2845), cex=0.5)
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2845)
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2844)
      
      #“Core Sylvioidea (excluding Hirundinidae and Phylloscopidae).”
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2844), node = 178, state = 1, anc.state = 0), fsize=0.4, type='arc')
      #text(x = 0, y = -1, labels = "Core Sylvioidea -\n Babbler group", cex = 2)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2843), node = 190, state = 1, anc.state = 0), fsize=0.4)
      
      
      data.frame(unique(clade_info$clade_members[1]$Node_2845$verbatim_name), clade_info$clade_members[1]$Node_2845$family)
      
    }
    
    #clade 2 - Tody-tyrants clade \n (Todirostrinae: Tyrannidae) -- source unclear
    {
      clade_info$unique_families[2]$Node_3790
      unique(clade_info$clade_members[2]$Node_3790$genus)
      unique(clade_info$clade_members[2]$Node_3790$verbatim_name)
      
      clade_info$sister_unique_families[2]$Node_3790
      unique(clade_info$sister_clade_members[2]$Node_3790$genus)
      unique(clade_info$sister_clade_members[2]$Node_3790$verbatim_name)
      
      #getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3790, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3790), cex=0.5)
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3790)
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3789)
      
      # tody-tyrants clade (Todirostrinae: Tyrannidae)
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3789), node = 36 , state = 1, anc.state = 0), fsize=0.5, type='arc')
      text(x = 0, y = -1, labels = "Tody-tyrants clade \n (Todirostrinae: Tyrannidae)", cex = 2)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3788), node = 60 , state = 1, anc.state = 0), fsize=0.5)
      
      
    }
    
    #clade 3 - Core Passeroidea \n Passerida? Ericsson 2003? Reference: Ericson et al. (2003). Evolution, Biogeography, and Patterns of Diversification in Passerine Birds.
    {
      clade_info$unique_families[3]$Node_2065
      unique(clade_info$clade_members[3]$Node_2065$genus)
      unique(clade_info$clade_members[3]$Node_2065$family)
      #unique(clade_info$clade_members[3]$Node_2065$verbatim_name)
      clade_info$sister_unique_families[3]$Node_2065
      unique(clade_info$sister_clade_members[3]$Node_2065$genus)
      #unique(clade_info$clade_members[3]$Node_2065$verbatim_name)
      #getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2065, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2064), cex=0.1)
      #“Core Passeroidea” (a major oscine radiation encompassing multiple avian families, distinct from the basal Petroicidae lineage).
      #Passerida
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2064)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2064), node = 1292 , state = 1, anc.state = 0), fsize=0.05, type='arc', lwd=0.5)
      text(x = 0, y = -1, labels = "Core Passeroidea \n Passerida?", cex = 2)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2063), node = 1498 , state = 1, anc.state = 0), fsize=0.05, type='arc', lwd=0.5)
      
      
      data.frame(clade_info$clade_members[3]$Node_2065$verbatim_name, clade_info$clade_members[3]$Node_2065$family)
      
    }
    
    #clade 4, core lonchura
    {
      clade_info$unique_families[4]$Node_2708
      unique(clade_info$clade_members[4]$Node_2708$genus)
      unique(clade_info$clade_members[4]$Node_2708$verbatim_name)
      clade_info$sister_unique_families[4]$Node_2708
      unique(clade_info$sister_clade_members[4]$Node_2708$genus)
      unique(clade_info$sister_clade_members[4]$Node_2708$verbatim_name)
      #a subset of Lonchura and Euodice within Estrildidae, core Lonchura
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2708, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2708))
      #Core Lonchura / excluding Lonchura nana, within Estrildidae -- 
      #lonchura nana is NOT lonchura any more -- so this is equivalent to Lonchura
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2707)
      
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2706), node = 20 , state = 1, anc.state = 0), fsize=0.5)
      text(x = 0, y = -1, labels = "Core Passeroidea \n Passerida?", cex = 2)
      
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2706))
      
      
    }
    
    #clade 5 - core vidua -- indigobirds
    {
      clade_info$unique_families[5]$Node_2740 
      unique(clade_info$clade_members[5]$Node_2740$genus)
      unique(clade_info$clade_members[5]$Node_2740$verbatim_name)
      clade_info$sister_unique_families[5]$Node_2740 
      unique(clade_info$sister_clade_members[5]$Node_2740$genus)
      unique(clade_info$sister_clade_members[5]$Node_2740$verbatim_name)
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2740, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2740))
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2736)
      
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2735))
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2735), node = 25 , state = 1, anc.state = 0), fsize=0.5)
      
      #“Vidua species group (excluding Vidua regia), an unranked subset of Vidua within Viduidae.”
    }
    
    #clade 6 - Acanthizidae
    {
      clade_info$unique_families[6]$Node_3594 
      unique(clade_info$clade_members[6]$Node_3594$genus)
      unique(clade_info$clade_members[6]$Node_3594$verbatim_name)
      clade_info$sister_unique_families[6]$Node_3594 
      unique(clade_info$sister_clade_members[6]$Node_3594$genus)
      unique(clade_info$sister_clade_members[6]$Node_3594$verbatim_name)
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3594, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3594))
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3593)
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3560))
      #Acanthizidae -- some lineages missing -- A tribe within Acanthizidae (either existing or newly designated).
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3560), node = 81 , state = 1, anc.state = 0), fsize=0.5)
      
      
    }
    
    #clade 7 - crow-magpie-nutcracker clade within Corvidae (“core corvids” (Corvinae sensu lato) Jønsson et al. (2016).)
    {
      clade_info$unique_families[7]$Node_3369 
      unique(clade_info$clade_members[7]$Node_3369$genus)
      unique(clade_info$clade_members[7]$Node_3369$verbatim_name)
      clade_info$sister_unique_families[7]$Node_3369 
      unique(clade_info$sister_clade_members[7]$Node_3369$genus)
      unique(clade_info$sister_clade_members[7]$Node_3369$verbatim_name)
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3369, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3369))
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3368)
      #the crow-magpie-nutcracker clade within Corvidae, historically recognized Corvinae sensu lato, 
      #traditionally distinguished from the jay-dominated sister clade (Cyanocoracinae). Core Corvids
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3369-1), node = 54, state = 1, anc.state = 0), fsize=0.5, type='arc')
      text(x = 0, y = -1, labels = "crow-magpie-nutcracker clade \n within Corvidae", cex = 1.5)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3367), node = 57 , state = 1, anc.state = 0), fsize=0.5)
      
      
    }
    
    #clade 8 - Sylviidae (sensu stricto) -  revised Sylviidae concept as proposed by Ohlson et al. (2013).
    {
      clade_info$unique_families[8]$Node_2895 
      unique(clade_info$clade_members[8]$Node_2895$genus)
      unique(clade_info$clade_members[8]$Node_2895$verbatim_name)
      clade_info$sister_unique_families[8]$Node_2895 
      unique(clade_info$sister_clade_members[8]$Node_2895$genus)
      unique(clade_info$sister_clade_members[8]$Node_2895$verbatim_name)
      #getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2895, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2895))
      
      #represents our sample of Sylviidae
      #Sylviidae (sensu stricto), an incompletely sampled but monophyletic subset of the family Sylviidae.
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2847)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2846), node = 140 , state = 1, anc.state = 0), fsize=0.5)
      
      
    }
    
    #clade 9 - new world jays, Ericson et al. (2005) and Ohlson et al. (2013)?
    {
      clade_info$unique_families[9]$Node_3398 
      unique(clade_info$clade_members[9]$Node_3398$genus)
      unique(clade_info$clade_members[9]$Node_3398$verbatim_name)
      clade_info$sister_unique_families[9]$Node_3398 
      unique(clade_info$sister_clade_members[9]$Node_3398$genus)
      unique(clade_info$sister_clade_members[9]$Node_3398$verbatim_name)
      #getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3398, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3368))
      #get_mrca_of_set(min10.ic20.gic$tree_no_uncertainty_untransformed, descendants = c(3369, 3398))
      #it aligns most closely with Cyanocoracinae sensu lato and can be described informally as the New World Jays.
      # jay clade within Corvidae
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3367)
      
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3366))
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3366), node = 59, state = 1, anc.state = 0), fsize=0.5)
      #text(x = 0, y = -1, labels = "Cyanocoracinae \n new world jays?", cex = 1.5)
      
    }
    
    #clade 10 - Emberizoidea (superfamily), sensu Ohlson et al. 2013, Ohlson, J. I., Fjeldså, J., & Ericson, P. G. P. (2013). “Phylogeny and classification of the New World suboscines (Aves, Passeriformes).” Zoologica Scripta, 42(1), 20-37.
    {
      clade_info$unique_families[10]$Node_2072 
      unique(clade_info$clade_members[10]$Node_2072$genus)
      #unique(clade_info$clade_members[10]$Node_2072$verbatim_name)
      clade_info$sister_unique_families[10]$Node_2072 
      unique(clade_info$sister_clade_members[10]$Node_2072$genus)
      #unique(clade_info$sister_clade_members[10]$Node_2072$verbatim_name)
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2072, mode=c('label', 'number'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2072))
      #Emberizoidea (superfamily),  this clade can confidently be referred to as an incomplete yet representative subset of Emberizoidea.
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2071)
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2070))
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2070), node = 595, state = 1, anc.state = 0), fsize=0.2, type='arc', lwd=0.75)
      text(x = 0, y = -1, labels = "Emberizoidea \n superfamily?", cex = 1.5)
      
    }
    
    #clade 11 -- cotingidae sensu lato, bulk of frugivorous, lek-forming cotingas, often considered part of the “true” cotingas in a broad sense. Ohlson et al. (2013)
    {
      clade_info$unique_families[11]$Node_3870 
      unique(clade_info$clade_members[11]$Node_3870$genus)
      unique(clade_info$clade_members[11]$Node_3870$verbatim_name)
      clade_info$sister_unique_families[11]$Node_3870 
      unique(clade_info$sister_clade_members[11]$Node_3870$genus)
      unique(clade_info$sister_clade_members[11]$Node_3870$verbatim_name)
      #Core Cotingidae Clade, encompassing a diverse group of genera while excluding 
      #Rupicola and Snowornis, which represent a separate lineage within the family.
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3870, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3870))
      ##this is one node away from cotinga MRCA
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3868)
      
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3868))
      3868
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3868), node = 31, state = 1, anc.state = 0), fsize=0.5)
      
      
    }
    
    #clade 12 - Core Thamnophiloidea (Thamnophilidae + Conopophagidae) sensu Ohlson et al. (2013) -- antbird radiation, “Core Thamnophiloidea (Thamnophilidae + Conopophagidae) sensu Ohlson et al. (2013)”
    {
      clade_info$unique_families[12]$Node_4062 
      unique(clade_info$clade_members[12]$Node_4062$genus)
      unique(clade_info$clade_members[12]$Node_4062$verbatim_name)
      clade_info$sister_unique_families[12]$Node_4062 
      unique(clade_info$sister_clade_members[12]$Node_4062$genus)
      unique(clade_info$sister_clade_members[12]$Node_4062$verbatim_name)
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 4062, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 4062), cex=0.5)
      #Thamnophiloidea #Thamnophilidae" "Conopophagidae" -- #Thamnophilidae radiation within the infraorder Furnariida (??Thamnophiloidea)
      #Thamnophilida (Infraorder, Tyranni)
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 4062)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3923), node = 323, state = 1, anc.state = 0), fsize=0.35)
      text(x = 0, y = -1, labels = "Thamnophiloidea \n superfamily?", cex = 1.5)
      
    }
    
    #clade 13 - core wren clade? Henicorhina + allies clade?
    {
      clade_info$unique_families[13]$Node_3100 
      #unique(clade_info$clade_members[13]$Node_3100$genus)
      unique(clade_info$clade_members[13]$Node_3100$verbatim_name)
      
      clade_info$sister_unique_families[13]$Node_3100 
      #unique(clade_info$sister_clade_members[13]$Node_3100$genus)
      unique(clade_info$sister_clade_members[13]$Node_3100$verbatim_name)
      
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3100, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3100))
      #“Core wren clade (Troglodytidae: Henicorhina-Cyphorhinus group)”, Thryothorus group
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3098), node = 22, state = 1, anc.state = 0), fsize=0.7)
      text(x = 0, y = -1, labels = "Core wren clade", cex = 1.5)
      
      
    }
    
    #clade 14 --core ploceus weaver clade
    {
      clade_info$unique_families[14]$Node_2759
      unique(clade_info$clade_members[14]$Node_2759$genus)
      unique(clade_info$clade_members[14]$Node_2759$verbatim_name)
      clade_info$sister_unique_families[14]$Node_2759
      unique(clade_info$sister_clade_members[14]$Node_2759$genus)
      unique(clade_info$sister_clade_members[14]$Node_2759$verbatim_name)
      #Core Ploceus within Plocdidae
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2759, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2759))
      #Core Ploceus within Plocdidae -- close to the root (2 nodes)-- african ploceus weaver clade
      #African Ploceus weaver clade????
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2758)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2757), node = 20, state = 1, anc.state = 0), fsize=0.7)
      
      text(x = 0, y = -1, labels = "Core African Ploceus \n weaver clade", cex = 1.5)
      
      
    }
    
    #clade 15, Tangara-Thraupis clade within Thraupidae, mostly tangara
    {
      clade_info$unique_families[15]$Node_2447
      unique(clade_info$clade_members[15]$Node_2447$genus)
      unique(clade_info$clade_members[15]$Node_2447$verbatim_name)
      
      clade_info$sister_unique_families[15]$Node_2447
      unique(clade_info$sister_clade_members[15]$Node_2447$genus)
      unique(clade_info$sister_clade_members[15]$Node_2447$verbatim_name)
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 2447, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2447))
      #Tangara-Thraupis clade within Thraupidae, mostly tangara
      
      getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2446)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 2445), node = 38, state = 1, anc.state = 0), fsize=0.7)
      
      
      
    }
    
    #clade 16 Furnariinae (sensu lato), Derryberry et al. (2011), “core ovenbirds”  Furnariinae Gray, 1840
    {
      clade_info$unique_families[16]$Node_3934
      unique(clade_info$clade_members[16]$Node_3934$genus)
      unique(clade_info$clade_members[16]$Node_3934$verbatim_name)
      
      clade_info$sister_unique_families[16]$Node_3934
      unique(clade_info$sister_clade_members[16]$Node_3934$genus)
      unique(clade_info$sister_clade_members[16]$Node_3934$verbatim_name)
      
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3934, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3934), cex=0.5)
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = getParent(as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), 3934)), cex=0.5)
      #Core Furnariini (Furnariidae)-- nested group within furnariidae
      
      #plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3933), node = 68, state = 1, anc.state = 0), fsize=0.25, type='arc')
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3932), node = 88, state = 1, anc.state = 0), fsize=0.75)
      #text(x = 0, y = -1, labels = "Core ovenbirds \n Furnariinae (sensu lato)", cex = 1.5)
      
      Synallaxini
      Furnariini
      
      
    }
    
    #clade 17, Furnariida radiation, Ericson et al. (2006) – Diversification of the South American suboscine birds: The Furnariida-Thamnophilida split and their evolutionary history in Tyranni.
    {
      clade_info$unique_families[17]$Node_3924
      unique(clade_info$clade_members[17]$Node_3924$genus)
      unique(clade_info$clade_members[17]$Node_3924$verbatim_name)
      
      clade_info$sister_unique_families[17]$Node_3924
      unique(clade_info$sister_clade_members[17]$Node_3924$genus)
      unique(clade_info$sister_clade_members[17]$Node_3924$verbatim_name)
      
      #getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3924, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3924-1), cex=0.5)
      #“Core Furnariida Clade
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3924-1), node = 185, state = 1, anc.state = 0), fsize=0.4)
      text(x = 0, y = -1, labels = "Core Furnariida Clade \n(infraorder) ", cex = 1.5)
      
      
    }
    
    #clade 18 - Core Pipridae clade (Pipridae) -- Pipra-Manacus clade
    {
      clade_info$unique_families[18]$Node_3900
      unique(clade_info$clade_members[18]$Node_3900$genus)
      unique(clade_info$clade_members[18]$Node_3900$verbatim_name)
      clade_info$sister_unique_families[18]$Node_3900
      unique(clade_info$sister_clade_members[18]$Node_3900$genus)
      unique(clade_info$sister_clade_members[18]$Node_3900$verbatim_name)
      getSisters(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3900, mode=c('label'))
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3900), cex=0.5)
      #Core Pipridae clade (Pipridae) -- Pipra-Manacus clade
      
      getParent(min10.ic20.gic$tree_no_uncertainty_untransformed, node = 3899)
      
      plot(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3898), cex=0.5)
      
      plot(paintSubTree(extract.clade(phy=as.phylo(min10.ic20.gic$tree_no_uncertainty_untransformed), node = 3898), node = 21, state = 1, anc.state = 0), fsize=1)
      
      
    }
    
    --
      replicated early bursts link cliamate and body plan evolution---
      --
      
      
      
  }
  

}

#Post-hoc analysis of phenotypic integration
{

# Step 1: Fit posthoc models
runs_with_posthoc <- fitPosthocModels(
  run_list = runs,
  user_formula = 'trait_data[,c(-13)] ~ trait_data[,c(13)]',
  input_data = superdat.shuffle,
  num_cores = 8
)

# Step 2: Assemble vars_cors for each run
all_vars_cors <- generateVarsCorsList(runs_with_posthoc, 
                                      remove_high_corr = T, 
                                      corr_threshold = 0.95)

posthoc_integration_pgls<-phylolm(
  scale(log(rate)) ~ scale(log(vars)) + scale(fisher_z_transform(corrs)),
  data = all_vars_cors$min20.ic20.gic,
  phy = as.phylo(
    collapsePhylogenyByStates(min20.ic20.gic$tree_no_uncertainty_untransformed)
  ), model='BM'
)

posthoc_integration_gls<-lm(
  log(rate) ~ scale(log(vars)) + scale(fisher_z_transform(corrs)),
  data = all_vars_cors$min20.ic20.gic
)

posthoc_integration_plots<-plotVarsCors(all_vars_cors, 
                                        resid_sd_threshold_corrs = 2, 
                                        resid_sd_threshold_vars = 2, 
                                        point_alpha = 0.5, 
                                        ci_level = 0.99, point_size = 3.5)


quartz(width=11*0.9, height=5*0.9, file='posthoc_integration_plot.pdf', type='pdf')
print(posthoc_integration_plots$combined)
dev.off()

#looking at covariances in PC space -- generating heatmaps
{
  ## 1. Vectorize post-hoc covariance matrices and run PCA
  res_vec <- vectorize_posthoc_cov(
    runs_with_posthoc$min10.ic20.gic$posthoc,
    use_correlation = TRUE,
    min_n = 10
  )
  
  Sigma_mat <- res_vec$Sigma_mat
  reg_ids   <- res_vec$reg_ids
  reg_age   <- res_vec$reg_age   # numeric vector of ages, named by regime ID
  
  Sigma_pca <- prcomp(Sigma_mat, center = TRUE, scale. = TRUE)
  scores    <- Sigma_pca$x
  
  #output which taxa are on the top and bottom of PC1
  {
    ## pc1 scores named by regime ID
    pc1_scores <- setNames(scores[, 1], reg_ids)
    
    ## top 5 and bottom 5 regimes by PC1
    top_ids    <- names(sort(pc1_scores, decreasing = TRUE))[1:5]
    bottom_ids <- names(sort(pc1_scores, decreasing = FALSE))[1:5]
    
    ## helper to get tip labels for a given regime ID
    get_tips <- function(reg_id) {
      runs_with_posthoc$min10.ic20.gic$posthoc[[reg_id]]$corrSt$phy$tip.label
    }
    
    ## print results
    cat("Top 5 regimes on PC1:\n")
    for (id in top_ids) {
      cat("\nRegime", id, "(PC1 =", pc1_scores[id], ")\n")
      print(get_tips(id))
    }
    
    cat("\n\nBottom 5 regimes on PC1:\n")
    for (id in bottom_ids) {
      cat("\nRegime", id, "(PC1 =", pc1_scores[id], ")\n")
      print(get_tips(id))
    }
  }
  
  ## 2. PC1–PC2 plot colored by regime age
  age_palette <- colorRampPalette(c("blue", "red"))
  age_rank    <- rank(reg_age, ties.method = "average")
  cols        <- age_palette(length(reg_age))[age_rank]
  
  pdf(file='pc1_pc2_diagnostic.pdf', height=5*1.2, width=7*1.2)
  {
 
    par(mfrow = c(2, 3))
    
    ## PC1 vs PC2 (no r here)
    #r_pc1_pc2 <- cor(scores[, 2], scores[, 1])
    plot(
      scores[, 2], scores[, 1],
      ylab = "PC1",
      xlab = "PC2",
      pch  = 21,
      col  = "black",
      bg   = alpha("grey", 0.5),
      cex  = 1.5,
      main = "PCA of regime-specific\ncorrelation structure"
    )
    #abline(lm(scores[, 1] ~ scores[, 2]), col = "red", lwd = 2, lty = 2)
    #legend("topleft", legend = paste0("r = ", round(r_pc1_pc2, 3)), bty = "n", 
    #       inset  = c(-0.01, 0.01))
    #text(scores[, 1], scores[, 2], labels = reg_ids, pos = 3, cex = 0.6)
    
    ## PC1 vs age
    r_pc1_age <- cor(reg_age, scores[, 1])
    plot(
      reg_age, scores[, 1],
      xlab = "Regime age (Ma)",
      ylab = "PC1 score (covariance structure)",
      pch  = 21,
      col  = "black",
      bg   = alpha("grey", 0.5),
      cex  = 1.5,
      main = "PC1 vs age",
      xlim = c(0, 40)
    )
    abline(lm(scores[, 1] ~ reg_age), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste0("r = ", round(r_pc1_age, 3)), bty = "n", 
           inset  = c(-0.01, 0.01))
    
    ## PC2 vs age
    r_pc2_age <- cor(reg_age, scores[, 2])
    plot(
      reg_age, scores[, 2],
      xlab = "Regime age (Ma)",
      ylab = "PC2 score (covariance structure)",
      pch  = 21,
      col  = "black",
      bg   = alpha("grey", 0.5),
      cex  = 1.5,
      main = "PC2 vs age",
      xlim = c(0, 40)
    )
    abline(lm(scores[, 2] ~ reg_age), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste0("r = ", round(r_pc2_age, 3)), bty = "n",
           inset  = c(-0.01, 0.01))
    
    ## PC3 vs age
    r_pc3_age <- cor(reg_age, scores[, 3])
    plot(
      reg_age, scores[, 3],
      xlab = "Regime age (Ma)",
      ylab = "PC3 score (covariance structure)",
      pch  = 21,
      col  = "black",
      bg   = alpha("grey", 0.5),
      cex  = 1.5,
      main = "PC3 vs age",
      xlim = c(0, 40)
    )
    abline(lm(scores[, 3] ~ reg_age), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste0("r = ", round(r_pc3_age, 3)), bty = "n",
           inset  = c(-0.01, 0.01))
    
    ## PC4 vs age
    r_pc4_age <- cor(reg_age, scores[, 4])
    plot(
      reg_age, scores[, 4],
      xlab = "Regime age (Ma)",
      ylab = "PC4 score (covariance structure)",
      pch  = 21,
      col  = "black",
      bg   = alpha("grey", 0.5),
      cex  = 1.5,
      main = "PC4 vs age",
      xlim = c(0, 40)
    )
    abline(lm(scores[, 4] ~ reg_age), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste0("r = ", round(r_pc4_age, 3)), bty = "n",
           inset  = c(-0.01, 0.01))
    
    par(mfrow = c(1, 1))
    
    
  }
  dev.off()

  ## 4. Helper to rebuild PC loadings into a trait x trait matrix
  example_sigma <- runs_with_posthoc$min10.ic20.gic$posthoc[[ reg_ids[1] ]]$sigma$Pinv
  trait_names   <- rownames(example_sigma)

  # Overwrite with cleaner anatomical labels
  trait_names <- c(
    tarsus          = "Tibiotarsus",
    metatarsus      = "Tarsometatarsus",
    femur           = "Femur",
    humerus         = "Humerus",
    ulna            = "Ulna",
    radius          = "Radius",
    carpometacarpus = "Carpometacarpus",
    second_digit    = "2nd digit phalanx",
    cv.keel.1       = "Keel",
    cv.furcula.1    = "Furcula",
    sclerotic_ring  = "Sclerotic ring",
    cv.skull.1      = "Skull–bill length"
  )[trait_names]
  
  # Inspect PCs 1–5
  plot_pc_loadings_heatmap(Sigma_pca, trait_names, pcs = 1)
  
  pdf(file = "Supp4_b.pdf", height= 7*0.8, width = 18*0.8)
  plot_pc_loadings_heatmap_CH(Sigma_pca, trait_names, pcs = 1:4)
  dev.off()
  
  pdf(file = "Supp4_b_2.pdf", height= 7*0.55, width = 18*0.8)
  plot_pc_loadings_heatmap_CH_local(Sigma_pca, trait_names, pcs = 1:4, show_legend = F, show_dend = F, dend_size_mm = 10)
  dev.off()
  
  sort(
  scores[, 1]
  )
  
  
  runs_with_posthoc$min10.ic20.gic$posthoc$`4`$corrSt$phy$tip.label
  runs_with_posthoc$min10.ic20.gic$posthoc$`275`$corrSt$phy$tip.label
  runs_with_posthoc$min10.ic20.gic$posthoc$`203`$corrSt$phy$tip.label
  
  runs_with_posthoc$min10.ic20.gic$posthoc$`408`$corrSt$phy$tip.label
  runs_with_posthoc$min10.ic20.gic$posthoc$`241`$corrSt$phy$tip.label
  runs_with_posthoc$min10.ic20.gic$posthoc$`167`$corrSt$phy$tip.label
  
  
}

## 5. Diagnostic: does PC1 align with within-wing vs cross-module correlations?
{
  # Define mapping from internal names to pretty names (same as above)
  name_map <- c(
    tarsus          = "Tibiotarsus",
    metatarsus      = "Tarsometatarsus",
    femur           = "Femur",
    humerus         = "Humerus",
    ulna            = "Ulna",
    radius          = "Radius",
    carpometacarpus = "Carpometacarpus",
    second_digit    = "2nd digit phalanx",
    cv.keel.1       = "Keel",
    cv.furcula.1    = "Furcula",
    sclerotic_ring  = "Sclerotic ring",
    cv.skull.1      = "Skull–bill length"
  )
  
  # Define module membership using the pretty names
  wing_traits     <- c("Carpometacarpus", "2nd digit phalanx",
                       "Radius", "Ulna", "Humerus", "Keel", "Furcula")
  hindcran_traits <- c("Tibiotarsus", "Tarsometatarsus", "Femur",
                       "Skull–bill length", "Sclerotic ring")
  
  # Storage vectors aligned with reg_ids / scores[, 1]
  within_wing_mean  <- numeric(length(reg_ids))
  cross_module_mean <- numeric(length(reg_ids))
  
  for (i in seq_along(reg_ids)) {
    reg <- reg_ids[i]
    # Extract covariance, convert to correlation
    Sigma_reg <- runs_with_posthoc$min10.ic20.gic$posthoc[[reg]]$sigma$Pinv
    C_reg     <- cov2cor(Sigma_reg)
    
    # Map internal names to pretty names
    internal_names <- rownames(C_reg)
    pretty_names   <- name_map[internal_names]
    
    # Re-label the matrix with pretty names
    dimnames(C_reg) <- list(pretty_names, pretty_names)
    
    # Indices for modules
    wing_idx     <- pretty_names %in% wing_traits
    hindcran_idx <- pretty_names %in% hindcran_traits
    
    # Within-wing mean correlation (upper triangle only, no diagonal)
    if (sum(wing_idx) > 1) {
      W <- C_reg[wing_idx, wing_idx, drop = FALSE]
      within_wing_mean[i] <- mean(W[upper.tri(W, diag = FALSE)], na.rm = TRUE)
    } else {
      within_wing_mean[i] <- NA
    }
    
    # Cross-module mean correlation: hindlimb+cranial vs wing/flight
    if (sum(wing_idx) > 0 && sum(hindcran_idx) > 0) {
      HxW <- C_reg[hindcran_idx, wing_idx, drop = FALSE]
      cross_module_mean[i] <- mean(HxW, na.rm = TRUE)
    } else {
      cross_module_mean[i] <- NA
    }
  }
  
  # PC1 scores
  pc1_scores <- scores[, 1]
  
  # Correlations between PC1 and the two summary measures
  r_within   <- cor(pc1_scores, within_wing_mean,  use = "complete.obs")
  r_cross    <- cor(pc1_scores, cross_module_mean, use = "complete.obs")
  cat("\nDiagnostic correlations for PC1:\n")
  cat("cor(PC1, within-wing mean)      = ", r_within, "\n")
  cat("cor(PC1, hindcran <-> wing mean) = ", r_cross,  "\n")
  
  # Scatterplots
  pdf(file = "pc1_diagnostic.pdf", width=(6*1.4)*1.1, height=(3*1.5)*1.1)
  {
  par(mfrow = c(1, 2))
  
  plot(pc1_scores, within_wing_mean,
       xlab = "PC1 score",
       ylab = "Mean within-wing correlation",
       main = "PC1 vs within-wing",
       pch  = 21,
       col  = "black",
       bg   = adjustcolor("grey", alpha.f = 0.5),
       cex  = 1.5)
  abline(lm(within_wing_mean ~ pc1_scores), col = "red", lwd = 2, lty = 2)
  legend("topleft", legend = paste0("r = ", round(r_within, 3)),
         bty = "n", inset = c(-0.01, 0.01))
  
  plot(pc1_scores, cross_module_mean,
       xlab = "PC1 score",
       ylab = "Mean hindlimb+cranial ↔ wing correlation",
       main = "PC1 vs cross-module",
       pch  = 21,
       col  = "black",
       bg   = adjustcolor("grey", alpha.f = 0.5),
       cex  = 1.5)
  abline(lm(cross_module_mean ~ pc1_scores), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste0("r = ", round(r_cross, 3)),
         bty = "n", inset = c(-0.01, 0.01))
  
  par(mfrow = c(1, 1))
  }
  dev.off()
  
  
}

## Diagnostic for PC2: head–hindlimb vs head–wing correlations
{
  ## Diagnostic for PC2: head–hindlimb vs head–wing correlations
  
  # Re-use name_map from the PC1 block (assumed already defined)
  # Define cranial, hindlimb, and wing sets using the same pretty names
  cranial_traits <- c("Skull–bill length", "Sclerotic ring")
  hind_traits    <- c("Tibiotarsus", "Tarsometatarsus", "Femur")
  wing_traits    <- c("Carpometacarpus", "2nd digit phalanx",
                      "Radius", "Ulna", "Humerus", "Keel", "Furcula")
  
  # Storage vectors aligned with reg_ids / scores[, 2]
  cran_hind_mean <- numeric(length(reg_ids))  # head–hindlimb
  cran_wing_mean <- numeric(length(reg_ids))  # head–wing
  
  for (i in seq_along(reg_ids)) {
    reg <- reg_ids[i]
    
    # Extract covariance, convert to correlation
    Sigma_reg <- runs_with_posthoc$min10.ic20.gic$posthoc[[reg]]$sigma$Pinv
    C_reg     <- cov2cor(Sigma_reg)
    
    # Map internal names to pretty names
    internal_names <- rownames(C_reg)
    pretty_names   <- name_map[internal_names]
    
    # Re-label the matrix with pretty names
    dimnames(C_reg) <- list(pretty_names, pretty_names)
    
    # Indices for modules
    cran_idx <- pretty_names %in% cranial_traits
    hind_idx <- pretty_names %in% hind_traits
    wing_idx <- pretty_names %in% wing_traits
    
    # Mean head–hindlimb correlation
    if (sum(cran_idx) > 0 && sum(hind_idx) > 0) {
      CH <- C_reg[cran_idx, hind_idx, drop = FALSE]
      cran_hind_mean[i] <- mean(CH, na.rm = TRUE)
    } else {
      cran_hind_mean[i] <- NA
    }
    
    # Mean head–wing correlation
    if (sum(cran_idx) > 0 && sum(wing_idx) > 0) {
      CW <- C_reg[cran_idx, wing_idx, drop = FALSE]
      cran_wing_mean[i] <- mean(CW, na.rm = TRUE)
    } else {
      cran_wing_mean[i] <- NA
    }
  }
  
  # PC2 scores
  pc2_scores <- scores[, 2]
  
  # Correlations between PC2 and the two summary measures
  r_head_hind <- cor(pc2_scores, cran_hind_mean, use = "complete.obs")
  r_head_wing <- cor(pc2_scores, cran_wing_mean, use = "complete.obs")
  
  cat("\nDiagnostic correlations for PC2:\n")
  cat("cor(PC2, head–hindlimb mean) = ", r_head_hind, "\n")
  cat("cor(PC2, head–wing mean)      = ", r_head_wing, "\n")
  
  # Scatterplots (matching PC1 style)
  pdf(file = "pc2_diagnostic.pdf", width = (6*1.4)*1.1, height = (3*1.5)*1.1)
  {
    par(mfrow = c(1, 2))
    
    plot(pc2_scores, cran_hind_mean,
         xlab = "PC2 score",
         ylab = "Mean head–hindlimb correlation",
         main = "PC2 vs head–hindlimb",
         pch  = 21,
         col  = "black",
         bg   = adjustcolor("grey", alpha.f = 0.5),
         cex  = 1.5)
    abline(lm(cran_hind_mean ~ pc2_scores), col = "red", lwd = 2, lty = 2)
    legend("topleft", legend = paste0("r = ", round(r_head_hind, 3)),
           bty = "n", inset = c(-0.01, 0.01))
    
    plot(pc2_scores, cran_wing_mean,
         xlab = "PC2 score",
         ylab = "Mean head–wing correlation",
         main = "PC2 vs head–wing",
         pch  = 21,
         col  = "black",
         bg   = adjustcolor("grey", alpha.f = 0.5),
         cex  = 1.5)
    abline(lm(cran_wing_mean ~ pc2_scores), col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste0("r = ", round(r_head_wing, 3)),
           bty = "n", inset = c(-0.01, 0.01))
    
    par(mfrow = c(1, 1))
  }
  dev.off()
}

## Diagnostic: regime rate vs PCs 1–4
{
  # 1. Build a named vector of rates by regime/state
  rate_by_state <- setNames(all_vars_cors$min10.ic20.gic$rate,
                            all_vars_cors$min10.ic20.gic$State)
  
  # 2. Align rates with the PCA regimes (reg_ids) and PC scores
  rate_vec <- rate_by_state[reg_ids]        # in same order as scores rows
  keep     <- !is.na(rate_vec)
  
  rate_use <- rate_vec[keep]
  log_rate <- log(rate_use)
  
  pc_scores_use <- scores[keep, , drop = FALSE]  # subset PCs to same regimes
  
  # 3. Compute correlations for PCs 1–4
  cor_pcs <- sapply(1:4, function(k) {
    cor(pc_scores_use[, k], log_rate, use = "complete.obs")
  })
  names(cor_pcs) <- paste0("PC", 1:4)
  
  cat("\nDiagnostic correlations for log(rate) vs PCs 1–4:\n")
  print(cor_pcs)
  
  # 4. Scatterplots: log(rate) vs PC1–PC4
  par(mfrow = c(2, 2))
  for (k in 1:4) {
    x <- pc_scores_use[, k]
    plot(
      x, log_rate,
      xlab = paste0("PC", k, " score (covariance structure)"),
      ylab = "Log regime rate",
      pch  = 16,
      main = paste0("Rate vs PC", k)
    )
    abline(lm(log_rate ~ x), col = "red", lwd = 2)
    # optional: add the correlation in the corner
    legend("topleft",
           legend = paste0("r = ", round(cor_pcs[k], 3)),
           bty = "n")
  }
  par(mfrow = c(1, 1))
}


}


