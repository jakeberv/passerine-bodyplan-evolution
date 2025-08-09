#Set of primary helper functions and primary search functions 
{
#' Generate Painted Subtrees for Eligible Nodes in a Phylogenetic Tree
#'
#' This function identifies internal nodes in a rooted phylogenetic tree that have 
#' at least a specified minimum number of descendant tips and generates a list of 
#' trees with subtrees painted (annotated) at those nodes.
#' Each painted tree represents a potential sub-model for further analysis.
generatePaintedTrees <- function(tree, min_tips, state = "shift") {
  if (!is.rooted(tree)) {
    tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
  }
  
  getEligibleNodes <- function(tree, min_tips) {
    eligible_nodes <- c()
    for (node in 1:(Nnode(tree))) {
      internal_node <- node + Ntip(tree)
      descendants <- getDescendants(tree, internal_node)
      tip_descendants <- tree$tip.label[descendants[descendants <= Ntip(tree)]]
      if (length(tip_descendants) >= min_tips) {
        eligible_nodes <- c(eligible_nodes, internal_node)
      }
    }
    return(eligible_nodes)
  }
  
  eligible_nodes <- getEligibleNodes(tree, min_tips)
  cat(paste(length(eligible_nodes), "eligible nodes are detected", '\n'))
  
  painted_trees <- list()
  
  for (node in eligible_nodes) {
    tree_copy <- tree
    tree_copy <- paintSubTree(tree_copy, node, state)
    painted_trees[[paste("Node", node)]] <- tree_copy
  }
  
  cat(paste(length(painted_trees), "sub-models generated", '\n'))
  return(painted_trees)
}


#' Fit mvgls Model to a Painted Tree and Extract GIC Score
#'
#' This function fits a multivariate generalized least squares (mvgls) model using 
#' a SIMMAP-formatted (painted) phylogenetic tree and a matrix of trait data, and 
#' returns the fitted model along with its Generalized Information Criterion (GIC) score.
fitMvglsAndExtractGIC <- function(painted_tree, trait_data) {
  # Ensure trait_data is a matrix
  if (!is.matrix(trait_data)) {
    stop("trait_data must be a matrix.")
  }
  
  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }
  
  # Fit the mvgls model directly using the matrix
  model <- mvgls(trait_data ~ 1, tree = painted_tree, model = "BMM", method='LL')
  gic_value <- GIC(model)
  
  # Return a list containing the model and the GIC
  return(list(model = model, GIC = gic_value))
}

#' Fit mvgls Model to a Painted Tree and Extract BIC Score
#'
#' This function fits a multivariate generalized least squares (mvgls) model using 
#' a SIMMAP-formatted (painted) phylogenetic tree and a matrix of trait data, and 
#' returns the fitted model along with its Bayesian Information Criterion (BIC) score.
fitMvglsAndExtractBIC <- function(painted_tree, trait_data) {
  # Ensure trait_data is a matrix
  if (!is.matrix(trait_data)) {
    stop("trait_data must be a matrix.")
  }
  
  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }
  
  # Fit the mvgls model directly using the matrix
  model <- mvgls(trait_data ~ 1, tree = painted_tree, model = "BMM", method='LL')
  bic_value <- BIC(model)
  
  # Return a list containing the model and the GIC
  return(list(model = model, BIC = bic_value))
}


#' Fit mvgls Model Using a Formula and Extract GIC Score
#'
#' This function fits a multivariate generalized least squares (mvgls) model to a 
#' SIMMAP-formatted (painted) phylogenetic tree using a user-specified formula and 
#' trait data, and returns the fitted model along with its Generalized Information Criterion (GIC) score.
#' The function automatically switches between "BM" and "BMM" models depending on the number of painted regimes.
fitMvglsAndExtractGIC.formula <- function(formula, painted_tree, trait_data, ...) {
  # Ensure trait_data is a matrix
  #if (!is.matrix(trait_data)) {
  #  stop("trait_data must be a matrix.")
  #}
  
  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }
  
  # Validate that formula is provided and is a character
  if (missing(formula) || !is.character(formula)) {
    stop("A character formula must be provided.")
  }
  
  # Convert the string formula to an actual formula object
  formula_obj <- as.formula(formula)
  
  # Fit the mvgls model using the user-defined formula
  
  # On the baseline tree, we have to switch to model = BM bc there is only one painted regime
  # Then we switch back to BMM
  
  if(length(unique(getStates(tree=painted_tree))) == 1){
    model <- mvgls(formula_obj, tree = painted_tree, model = "BM", ...)
  } else {
    model <- mvgls(formula_obj, tree = painted_tree, model = "BMM", ...)
  }
  gic_value <- GIC(model)
  
  # Return a list containing the model and the GIC
  return(list(model = model, GIC = gic_value))
}

#' Fit mvgls Model Using a Formula and Extract BIC Score
#'
#' This function fits a multivariate generalized least squares (mvgls) model to a 
#' SIMMAP-formatted (painted) phylogenetic tree using a user-specified formula and 
#' trait data, and returns the fitted model along with its Bayesian Information Criterion (BIC) score.
#' The model type is automatically selected: "BM" is used if there is only one regime in the tree, otherwise "BMM" is applied.
fitMvglsAndExtractBIC.formula <- function(formula, painted_tree, trait_data, ...) {
  # # Ensure trait_data is a matrix
  # if (!is.matrix(trait_data)) {
  #   stop("trait_data must be a matrix.")
  # }
  
  # Make sure the row names match the tip labels of the painted_tree
  if (!identical(rownames(trait_data), painted_tree$tip.label)) {
    stop("Row names of trait_data must exactly match the tip labels of the tree.")
  }
  
  # Validate that formula is provided and is a character
  if (missing(formula) || !is.character(formula)) {
    stop("A character formula must be provided.")
  }
  
  # Convert the string formula to an actual formula object
  formula_obj <- as.formula(formula)
  
  # Fit the mvgls model using the user-defined formula
  
  # On the baseline tree, we have to switch to model = BM bc there is only one painted regime
  # Then we switch back to BMM
  
  if(length(unique(getStates(tree=painted_tree))) == 1){
    model <- mvgls(formula_obj, tree = painted_tree, model = "BM", ...)
  } else {
    model <- mvgls(formula_obj, tree = painted_tree, model = "BMM", ...)
  }
  bic_value <- BIC(model)
  
  # Return a list containing the model and the GIC
  return(list(model = model, BIC = bic_value))
}

#' Calculate Delta GIC Scores Relative to a Baseline Model
#'
#' This function computes the difference in Generalized Information Criterion (GIC) scores 
#' between a baseline model (assumed to be the first in the list) and a series of alternative 
#' models fitted to painted phylogenetic trees. It returns a named vector of delta GIC values, 
#' representing the improvement or decline in model fit relative to the baseline.
calculateAllDeltaGIC <- function(model_results, painted_tree_list) {
  # Check if the painted trees have names
  if (!all(sapply(painted_tree_list, function(x) !is.null(names(x))))) {
    stop("All trees in 'painted_tree_list' must have names.")
  }
  
  # Extract the names from the painted_tree_list
  tree_list_names <- names(painted_tree_list)
  
  # Retrieve the baseline GIC from the first model result
  # Ensure we are accessing the numeric GIC value correctly
  baseline_gic <- model_results[[1]]$GIC$GIC
  if (!is.numeric(baseline_gic)) {
    stop("The baseline GIC value must be numeric.")
  }
  
  # Use lapply to iterate over the list of model_results and calculate the delta GIC
  delta_gic <- lapply(seq_along(model_results), function(i) {
    # Ensure we are accessing the numeric GIC value correctly
    current_gic <- model_results[[i]]$GIC$GIC
    if (!is.numeric(current_gic)) {
      stop(paste("The GIC value for model", i, "must be numeric."))
    }
    # Calculate the difference in GIC between the baseline and the shift model
    return(baseline_gic - current_gic)
  })
  
  # Assign names to the delta GIC values and convert it to a named vector
  delta_gic <- setNames(unlist(delta_gic), tree_list_names)
  
  return(delta_gic)
}

#' Paint a Subtree in a Phylogenetic Tree with Optional Selective Overwriting
#'
#' This function modifies a phylogenetic tree (of class "phylo") by painting a specified subtree 
#' starting at a given node with a new state, optionally preserving existing state mappings unless overwritten.
#' It returns a SIMMAP-style tree with updated edge mappings and supports both full and selective painting, 
#' as well as optional stem painting from the parent edge.
paintSubTree_mod <- function(tree, node, state, anc.state="1", stem=FALSE, overwrite=TRUE) {
  if (!inherits(tree, "phylo")) stop("tree should be an object of class \"phylo\".")
  if (stem == 0 && node <= length(tree$tip)) stop("stem must be TRUE for node <= N")
  if (is.null(tree$edge.length)) tree <- compute.brlen(tree)
  
  if (is.null(tree$maps)) {
    maps <- as.list(tree$edge.length)
    for (i in 1:length(maps)) names(maps[[i]]) <- anc.state
  } else {
    maps <- tree$maps
  }
  
  if (overwrite) {
    # Original behavior: Overwrite entire subtree
    desc <- getDescendants(tree, node)
    z <- which(tree$edge[,2] %in% desc)
    for (i in z) {
      maps[[i]] <- sum(maps[[i]])
      names(maps[[i]]) <- state
    }
  } else {
    # Modified behavior: Selective overwriting
    target_state <- if (node > length(tree$tip)) names(maps[[which(tree$edge[,2] == node)]]) else anc.state
    desc <- getDescendants(tree, node)
    z <- which(tree$edge[,2] %in% desc)
    for (i in z) {
      if (names(maps[[i]]) == target_state) {
        maps[[i]] <- sum(maps[[i]])
        names(maps[[i]]) <- state
      }
    }
  }
  
  if (stem && node > length(tree$tip)) {
    stem_edge <- which(tree$edge[,2] == node)
    maps[[stem_edge]] <- sum(maps[[stem_edge]]) * c(1 - stem, stem)
    names(maps[[stem_edge]]) <- c(anc.state, state)
  }
  
  s <- vector()
  for (i in 1:nrow(tree$edge)) s <- c(s, names(maps[[i]]))
  s <- unique(s)
  mapped.edge <- matrix(0, length(tree$edge.length), length(s), dimnames=list(edge=apply(tree$edge, 1, function(x) paste(x, collapse=",")), state=s))
  for (i in 1:length(maps)) {
    for (j in 1:length(maps[[i]])) {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }
  
  tree$mapped.edge <- mapped.edge
  tree$maps <- maps
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  
  return(tree)
}

#' Remove a Painted Shift from a Subtree in a SIMMAP-Formatted Phylogenetic Tree
#'
#' This function removes a previously painted shift (regime change) from a specified node and its descendant branches 
#' in a SIMMAP-style phylogenetic tree by selectively overwriting branches painted with the shift state.
#' The original ancestral state (inherited from the parent node) is restored, with optional stem painting of the edge leading to the node.
paintSubTree_removeShift <- function(tree, shift_node, stem=FALSE) {
  if (!inherits(tree, "phylo")) stop("tree should be an object of class 'phylo'.")
  if (is.null(tree$edge.length)) tree <- compute.brlen(tree)
  
  if (is.null(tree$maps)) {
    maps <- as.list(tree$edge.length)
    for (i in 1:length(maps)) names(maps[[i]]) <- "1"  # Assuming '1' is the default ancestral state
  } else {
    maps <- tree$maps
  }
  
  # Get parent node and its state
  parent_node <- phytools::getParent(tree, shift_node)
  parent_state <- if (!is.na(parent_node)) getStates(tree, type = "nodes")[as.character(parent_node)] else "1"  # Default to '1' if parent is NA
  
  # Handle stem painting if applicable
  if (stem && shift_node > length(tree$tip)) {
    stem_edge <- which(tree$edge[,2] == shift_node)
    maps[[stem_edge]] <- sum(maps[[stem_edge]]) * c(1 - stem, stem)
    names(maps[[stem_edge]]) <- c(parent_state, parent_state)
  }
  
  # Get descendants and selectively overwrite branches
  desc <- getDescendants(tree, shift_node)
  shift_node_state <- getStates(tree, type = "nodes")[as.character(shift_node)]
  
  z <- which(tree$edge[,2] %in% desc)
  for (i in z) {
    if (names(maps[[i]]) == shift_node_state) {
      maps[[i]] <- sum(maps[[i]])
      names(maps[[i]]) <- parent_state
    }
  }
  
  # Update the tree with the new maps
  s <- vector()
  for (i in 1:nrow(tree$edge)) s <- c(s, names(maps[[i]]))
  s <- unique(s)
  mapped.edge <- matrix(0, length(tree$edge.length), length(s), dimnames = list(edge = apply(tree$edge, 1, function(x) paste(x, collapse = ",")), state = s))
  for (i in 1:length(maps)) {
    for (j in 1:length(maps[[i]])) {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }
  
  tree$mapped.edge <- mapped.edge
  tree$maps <- maps
  class(tree) <- c("simmap", setdiff(class(tree), "simmap"))
  
  return(tree)
}

#' Add a New Shift to a Phylogenetic Tree Model
#'
#' This function adds a new evolutionary regime (shift) to a SIMMAP-style phylogenetic tree 
#' by painting the subtree starting at a specified node with a new unique state identifier.
#' The shift ID is incremented and returned for tracking multiple shifts in downstream modeling.
addShiftToModel <- function(tree, shift_node, current_shift_id) {
  # Update the shift ID
  next_shift_id <- current_shift_id + 1
  
  # Paint the subtree with the new regime/shift id
  painted_tree <- paintSubTree_mod(tree, node = shift_node, state = as.character(next_shift_id), overwrite=F, stem = F)
  
  # Return a list with the updated tree and the new shift ID
  return(list(tree = painted_tree, shift_id = next_shift_id))
}

#' Remove a Shift from a Phylogenetic Tree by Reverting to Parent State
#'
#' This function removes a regime shift from a SIMMAP-style phylogenetic tree by repainting the subtree 
#' starting at a specified node with the state of its parent node. It uses `paintSubTree_removeShift()` 
#' to selectively revert the painted shift without affecting unrelated branches, with optional stem edge handling.
removeShiftFromTree <- function(tree, shift_node, stem=F) {
  #print(paste("Removing shift from node:", shift_node))
  
  # Retrieve node states; names of this vector are node indices
  node_states <- getStates(tree, type="nodes")
  #print("Current node states:")
  #print(node_states)
  
  # Using phytools' getParent to find the parent node index
  parent_node <- phytools::getParent(tree, shift_node)
  #print(paste("Parent node of", shift_node, "is:", parent_node))
  
  # Indicate the state of the shift node
  shift_state <- node_states[as.character(shift_node)]
  #print(paste("State of shift node", shift_node, "is:", shift_state))
  
  # Check if parent node index is valid
  if (!is.na(parent_node) && parent_node %in% names(node_states)) {
    # Get the state of the parent node
    parent_state <- node_states[as.character(parent_node)]
    #print(paste("State of parent node", parent_node, "is:", parent_state))
    
    # Check if parent state is NA
    if (!is.na(parent_state)) {
      # Paint the subtree at the shift node with the parent's state, without overwriting descendants
      #print(paste("Painting subtree at node", shift_node, "with state", parent_state, "to remove shift"))
      tree <- paintSubTree_removeShift(tree, shift_node, stem=stem)  # Using the specialized function for shift removal
    } else {
      #print(paste("State of parent node", parent_node, "is NA. Cannot remove shift."))
    }
  } else {
    #print("Invalid parent node. Cannot remove shift.")
  }
  
  return(tree)
}

#' Identify Nodes Representing Shifts in a SIMMAP Tree
#'
#' This function identifies internal nodes in a SIMMAP-style phylogenetic tree 
#' that represent regime shifts, based on the most recent common ancestors (MRCAs) 
#' of tips sharing the same painted state.
whichShifts <- function(tree) {
  tip_states <- getStates(tree, type = "tips")
  unique_states <- unique(tip_states)
  
  shift_nodes <- c()
  for (state in unique_states) {
    tips_with_state <- names(tip_states[tip_states == state])
    if (length(tips_with_state) > 1) {
      # Use ape::getMRCA to find the most recent common ancestor
      mrca_node <- ape::getMRCA(tree, tips_with_state)
      shift_nodes <- c(shift_nodes, mrca_node)
    }
  }
  
  return(unique(shift_nodes))
}

#' Extract Regime-Specific Variance-Covariance Matrices from a BMM mvgls Model
#'
#' This function extracts and scales the regime-specific variance-covariance (VCV) matrices 
#' from a mvgls model fitted under the BMM (Brownian Motion with multiple regimes) framework.
#' It returns a named list of VCV matrices, each corresponding to an evolutionary regime.
extractRegimeVCVs <- function(model_output) {
  # Ensure the required components are in the model_output
  if (!"param" %in% names(model_output) || !"sigma" %in% names(model_output) || !"Pinv" %in% names(model_output$sigma)) {
    return(NULL)
    stop("model_output does not contain the required components.")
  }
  
  # Extract the precision matrix (Pinv) for the first regime (base VCV)
  base_Pinv <- model_output$sigma$Pinv
  
  # Get the parameter for the first regime (base rate)
  base_param <- model_output$param[1]
  
  # List to store VCVs for each regime
  vcv_list <- list()
  
  # Iterate through the parameters and calculate VCV for each regime
  param_names <- names(model_output$param)
  for (i in seq_along(param_names)) {
    regime_name <- param_names[i]
    regime_param <- model_output$param[regime_name]
    
    # Scale the precision matrix to get the regime's covariance matrix
    # For the first regime, no scaling is needed
    if (i == 1) {
      regime_vcv <- base_Pinv
    } else {
      regime_vcv <- base_Pinv * (regime_param / base_param)
    }
    
    # Add the VCV matrix to the list, using the regime's parameter name
    vcv_list[[regime_name]] <- regime_vcv
  }
  
  return(vcv_list)
}

# # Supporting function getDescendants (if not already defined in the environment)
#' Get All Descendants of a Node in a Phylogenetic Tree
#'
#' Recursively retrieves all descendant nodes (internal and/or tips) from a specified node 
#' in a phylogenetic tree of class `"phylo"`. Optionally includes the starting node in the output.
getDescendants <- function(tree, node, include.node = FALSE) {
  # Function to recursively find all descendants of a node
  descendants <- numeric(0)
  for (i in which(tree$edge[,1] == node)) {
    descendants <- c(descendants, tree$edge[i,2], getDescendants(tree, tree$edge[i,2]))
  }
  if (include.node) descendants <- c(node, descendants)
  return(unique(descendants))
}

#' Search for the Optimal Shift Configuration in a Phylogenetic Model
#'
#' This function performs a stepwise search to identify the optimal configuration of evolutionary regime shifts 
#' on a phylogenetic tree, using information criteria (GIC or BIC) to guide model selection. It evaluates candidate 
#' shifts in parallel, fits models at each stage, optionally accounts for model uncertainty, and computes model weights.
searchOptimalConfiguration <-
  function(baseline_tree,
           trait_data,
           formula = 'trait_data~1',
           min_descendant_tips,
           num_cores = 2,
           ic_uncertainty_threshold = 1.0,
           shift_acceptance_threshold = 1.0,
           uncertainty = F,
           uncertaintyweights = F,
           uncertaintyweights_par = F,
           postorder_traversal = F,
           plot = T,
           IC = 'GIC', 
           store_model_fit_history = TRUE, ...) {
    
    # Capture user input
    user_input <- as.list(match.call())
    
    #generate initial set of painted candidate trees with shifts at each sub-node
    cat('Generating candidate shift models...\n')
    candidate_trees <- generatePaintedTrees(baseline_tree, min_descendant_tips)
    candidate_trees_shifts <- candidate_trees[-1]
    
    #fit the initial baseline model to the baseline tree with a global regime (state=0)
    cat('Fitting baseline model...\n')
    
    #select which information criterion to use
    if (IC != "GIC" && IC != "BIC") {
      stop("IC must be GIC or BIC")
    }
    if(IC=="GIC"){
      baseline_model <- fitMvglsAndExtractGIC.formula(formula, candidate_trees[[1]], trait_data, ...)
      baseline_ic <- baseline_model$GIC$GIC
    }
    if(IC=="BIC"){
      baseline_model <- fitMvglsAndExtractBIC.formula(formula, candidate_trees[[1]], trait_data, ...)
      baseline_ic <- baseline_model$BIC$BIC
    }
    cat(paste('Baseline ', IC, ':', round(baseline_ic, digits = 2), '\n'))
    
    #evaluate all of the candidate trees under GIC or BIC
    # Capture additional arguments into a list
    args_list <- list(...)
    
    if (.Platform$OS.type == "unix") {
      plan(multicore, workers = num_cores)
    } else {
      plan(multisession, workers = num_cores)
    }
    cat('Fitting sub-models in parallel...\n')
    candidate_results <-
      future.apply::future_lapply(candidate_trees_shifts, function(tree) {
        if (IC == "GIC") {
          do.call(fitMvglsAndExtractGIC.formula, c(list(formula, tree, trait_data), args_list))
        } else if (IC == "BIC") {
          do.call(fitMvglsAndExtractBIC.formula, c(list(formula, tree, trait_data), args_list))
        }
      }, future.seed = TRUE, future.scheduling = TRUE)
    plan(sequential)
    
    #generate the delta IC lists
    if (IC == "GIC"){
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$GIC$GIC)
    } else if (IC == "BIC") {
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$BIC$BIC)
    }
    
    cat('Sorting and evaluating shifts...\n')
    sorted_candidates <- candidate_trees_shifts[order(delta_ic_list, decreasing=T)]
    current_best_tree <- baseline_tree
    current_best_ic <- baseline_ic
    shift_id <- 0
    
    #plot current shift tree
    if(plot==T){
      plotSimmap(current_best_tree, ftype='off')
    }
    
    #in case the user wants postorder traversal for shift searching (default off)
    if(postorder_traversal){
      print('Candidate shifts are sorted in Postorder')
      sorted_candidates <- (candidate_trees_shifts)
      current_best_tree <- baseline_tree
      current_best_ic <- baseline_ic
      shift_id <- 0
    }
    
    shift_vec <- list() #initialize shift_vec
    model_with_shift_no_uncertainty<-NULL #initialize output
    best_tree_no_uncertainty<-NULL #initialize output
    # Initialize the list to collect warning messages
    warnings_list <- list()
    model_fit_history<-list() #new object to capture the full history of the search
    
    #Run the primary shift configuration search
    
    # for (i in seq_along(sorted_candidates)) {
    #   shift_node_name <- names(sorted_candidates)[i]
    #   shift_node_number <- as.integer(sub("Node ", "", shift_node_name))
    #   percent_complete <- round((i / length(sorted_candidates)) * 100, 2)
    #   print(paste('Evaluating shift at node:', shift_node_number, '-', percent_complete, '% complete'))
    #   
    #   add_shift_result <- addShiftToModel(current_best_tree, shift_node_number, shift_id)
    #   shifted_tree <- add_shift_result$tree
    #   shift_id <- add_shift_result$shift_id
    #   
    #   if(plot==T){
    #     nodelabels(text = shift_id, node=shift_node_number)
    #   }
    #   
    #   # tryCatch({
    #   #   if (IC == "GIC") {
    #   #     model_with_shift <- fitMvglsAndExtractGIC.formula(formula, shifted_tree, trait_data, ...)
    #   #     new_ic <- model_with_shift$GIC$GIC
    #   #   } else if (IC == "BIC") {
    #   #     model_with_shift <- fitMvglsAndExtractBIC.formula(formula, shifted_tree, trait_data, ...)
    #   #     new_ic <- model_with_shift$BIC$BIC
    #   #   }
    #   #   
    #   #   delta_ic <- current_best_ic - new_ic
    #   #   
    #   #   if (delta_ic >= shift_acceptance_threshold) {
    #   #     current_best_tree <- shifted_tree
    #   #     current_best_ic <- new_ic
    #   #     print(paste('Shift at node', shift_node_number, 'accepted. Updated', IC, ':', current_best_ic, 'Delta', IC, ':', delta_ic))
    #   #     
    #   #     shift_vec[[length(shift_vec) + 1]] <- shift_node_number
    #   #     
    #   #     best_tree_no_uncertainty <- current_best_tree
    #   #     model_with_shift_no_uncertainty <- model_with_shift
    #   #   } else {
    #   #     print(paste('Shift at node', shift_node_number, 'rejected. Delta', IC, ':', delta_ic, 'is less than threshold:', shift_acceptance_threshold))
    #   #   }
    #   # }, error = function(e) {
    #   #   warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
    #   #   warning(warning_message)
    #   #   warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #   # }, warning = function(w) {
    #   #   warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
    #   #   warning(warning_message)
    #   #   warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #   # })
    #   tryCatch({
    #     if (IC == "GIC") {
    #       model_with_shift <- withCallingHandlers(
    #         fitMvglsAndExtractGIC.formula(formula, shifted_tree, trait_data, ...),
    #         warning = function(w) {
    #           # Log the warning but continue execution
    #           warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
    #           warning(warning_message)
    #           warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #           invokeRestart("muffleWarning")  # Prevent the warning from propagating further
    #         }
    #       )
    #       new_ic <- model_with_shift$GIC$GIC
    #     } else if (IC == "BIC") {
    #       model_with_shift <- withCallingHandlers(
    #         fitMvglsAndExtractBIC.formula(formula, shifted_tree, trait_data, ...),
    #         warning = function(w) {
    #           # Log the warning but continue execution
    #           warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
    #           warning(warning_message)
    #           warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #           invokeRestart("muffleWarning")  # Prevent the warning from propagating further
    #         }
    #       )
    #       new_ic <- model_with_shift$BIC$BIC
    #     }
    #     
    #     # Calculate delta IC
    #     delta_ic <- current_best_ic - new_ic
    #     
    #     # Decision logic
    #     if (delta_ic >= shift_acceptance_threshold) {
    #       current_best_tree <- shifted_tree
    #       current_best_ic <- new_ic
    #       print(paste('Shift at node', shift_node_number, 'accepted. Updated', IC, ':', current_best_ic, 'Delta', IC, ':', delta_ic))
    #       
    #       shift_vec[[length(shift_vec) + 1]] <- shift_node_number
    #       
    #       best_tree_no_uncertainty <- current_best_tree
    #       model_with_shift_no_uncertainty <- model_with_shift
    #     } else {
    #       print(paste('Shift at node', shift_node_number, 'rejected. Delta', IC, ':', delta_ic, 'is less than threshold:', shift_acceptance_threshold))
    #     }
    #   }, error = function(e) {
    #     # Handle errors
    #     warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
    #     warning(warning_message)
    #     warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #   })
    #   
    #   
    #   if(plot==T){
    #     colorvec <- setNames(object = c('black', rainbow(length(unique(getStates(shifted_tree, type = 'both')))-1)),
    #                          nm = sort(as.numeric(unique(getStates(shifted_tree, type = 'both')))))
    #     plotSimmap(current_best_tree, colors=colorvec, fsize=0.0001)
    #   }
    # }
    
    for (i in seq_along(sorted_candidates)) {
      shift_node_name <- names(sorted_candidates)[i]
      shift_node_number <- as.integer(sub("Node ", "", shift_node_name))
      percent_complete <- round((i / length(sorted_candidates)) * 100, 2)
      cat(paste('Evaluating shift at node:', shift_node_number, '-', percent_complete, '% complete', '\n'))
      
      add_shift_result <- addShiftToModel(current_best_tree, shift_node_number, shift_id)
      shifted_tree <- add_shift_result$tree
      shift_id <- add_shift_result$shift_id
      
      if(plot==T){
        nodelabels(text = shift_id, node=shift_node_number)
      }
      
      tryCatch({
        if (IC == "GIC") {
          model_with_shift <- withCallingHandlers(
            fitMvglsAndExtractGIC.formula(formula, shifted_tree, trait_data, ...),
            warning = function(w) {
              warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
              warning(warning_message)
              warnings_list[[length(warnings_list) + 1]] <<- warning_message
              invokeRestart("muffleWarning") 
            }
          )
          new_ic <- model_with_shift$GIC$GIC
        } else if (IC == "BIC") {
          model_with_shift <- withCallingHandlers(
            fitMvglsAndExtractBIC.formula(formula, shifted_tree, trait_data, ...),
            warning = function(w) {
              warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
              warning(warning_message)
              warnings_list[[length(warnings_list) + 1]] <<- warning_message
              invokeRestart("muffleWarning") 
            }
          )
          new_ic <- model_with_shift$BIC$BIC
        }
        
        # Calculate delta IC
        delta_ic <- current_best_ic - new_ic
        
        # Store model fit and acceptance status (including delta_ic)
        if (store_model_fit_history) {
          model_fit_history[[length(model_fit_history) + 1]] <- list(
            model = model_with_shift,
            accepted = delta_ic >= shift_acceptance_threshold,
            delta_ic = delta_ic 
          )
        }
        
        # Decision logic (unchanged)
        if (delta_ic >= shift_acceptance_threshold) {
          current_best_tree <- shifted_tree
          current_best_ic <- new_ic
          cat(paste('Shift at node', shift_node_number, 'accepted. Updated', IC, ':', round(current_best_ic, digits = 2), 'Delta', IC, ':', round(delta_ic, digits = 2), '\n')) 
          
          shift_vec[[length(shift_vec) + 1]] <- shift_node_number
          
          best_tree_no_uncertainty <- current_best_tree
          model_with_shift_no_uncertainty <- model_with_shift
        } else {
          cat(paste('Shift at node', shift_node_number, 'rejected. Delta', IC, ':', round(delta_ic, digits = 2), 'is less than threshold:', shift_acceptance_threshold, '\n')) 
        }
      }, error = function(e) {
        # Handle errors (unchanged)
        warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
        warning(warning_message)
        warnings_list[[length(warnings_list) + 1]] <<- warning_message
        
        # Also store the error in the model fit history
        if (store_model_fit_history) {
          model_fit_history[[length(model_fit_history) + 1]] <- list(
            model = NULL, 
            accepted = FALSE,
            delta_ic = NA,
            error = e$message
          )
        }
        
      })
      
      if(plot==T){
        colorvec <- setNames(object = c('black', rainbow(length(unique(getStates(shifted_tree, type = 'both')))-1)),
                             nm = sort(as.numeric(unique(getStates(shifted_tree, type = 'both')))))
        plotSimmap(current_best_tree, colors=colorvec, ftype='off')
      }
    }
    
    #print(paste(shift_vec))
    shifts_no_uncertainty<-unlist(shift_vec)
    cat(paste("Shifts detected at nodes:", paste(shift_vec, collapse = ", "), '\n'))
    
    # If activated, this section removes shifts after re-evaluating
    model_without_shift <- NULL # Initialize as NULL
    if (uncertainty) {
      cat('Post-search re-evaluation to reduce overfitting...')
      shift_nodes <- unlist(shift_vec)
      print(paste("Re-evaluating nodes", shift_nodes))
      shift_vec_uncertainty <- shift_nodes # Tracker (to remove shifts from)
      
      root_node <- Ntip(baseline_tree) + 1
      for (shift_node_number in shift_nodes) {
        if (shift_node_number != root_node) {
          cat(paste('Re-evaluating shift at node:', shift_node_number))
          tree_without_shift <- removeShiftFromTree(current_best_tree, shift_node_number)
          
          tryCatch({
            # Evaluate the model without the current shift using the selected IC
            if (IC == "GIC") {
              temp_model <- fitMvglsAndExtractGIC.formula(formula, tree_without_shift, trait_data, ...)
              ic_without_shift <- temp_model$GIC$GIC
            } else if (IC == "BIC") {
              temp_model <- fitMvglsAndExtractBIC.formula(formula, tree_without_shift, trait_data, ...)
              ic_without_shift <- temp_model$BIC$BIC
            }
            
            if (abs(current_best_ic - ic_without_shift) <= ic_uncertainty_threshold) {
              current_best_tree <- tree_without_shift
              current_best_ic <- ic_without_shift
              cat(paste('Shift at node', shift_node_number, 'removed. Updated', IC, ':', round(current_best_ic, digits=2)))
              shift_vec_uncertainty <- shift_vec_uncertainty[!shift_vec_uncertainty == shift_node_number]
              model_without_shift <- temp_model # Update only if the condition is met
            }
          }, error = function(e) {
            warning(paste("Error in re-evaluating shift at node", shift_node_number, ":", e$message))
          })
        } else {
          cat(paste('Skipping re-evaluation for the root node:', shift_node_number))
        }
      }
    }
    
    # New Section for Calculating Information Criterion Weights Post Optimization
    ic_weights_df <- NA  # Initialize the results vector
    if (xor(uncertaintyweights, uncertaintyweights_par)) {
      if (uncertaintyweights) {
        if (length(unlist(shift_vec)) > 0) {
          cat('Calculating IC weights for initially identified shifts...\n')
          ic_weights_df <- data.frame(node = integer(), ic_with_shift = numeric(), ic_without_shift = numeric(), delta_ic = numeric(), ic_weight = numeric())
          
          # Retrieve the IC of the optimized model before uncertainty analysis
          original_ic <- if (IC == "GIC") model_with_shift_no_uncertainty$GIC$GIC else model_with_shift_no_uncertainty$BIC$BIC
          
          cat(paste("Considering", length(shift_vec), "shifts in the candidate set",'\n'))
          cat(paste("There are", length(unique(getStates(best_tree_no_uncertainty, type = 'both'))) - 1, "shifts in the mapped tree",'\n'))
          
          for (shift_node_number in unlist(shift_vec)) {
            cat(paste('Re-estimating model without shift at node:', shift_node_number, '\n'))
            
            # Remove the shift temporarily from best_tree_no_uncertainty and re-estimate the IC
            tree_without_current_shift <- removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
            model_without_current_shift_function <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
            model_without_current_shift <- model_without_current_shift_function(formula, tree_without_current_shift, trait_data, ...)
            ic_without_current_shift <- if (IC == "GIC") model_without_current_shift$GIC$GIC else model_without_current_shift$BIC$BIC
            
            # Calculate the difference in IC
            delta_ic <- original_ic - ic_without_current_shift
            
            # Use aicw function to calculate the IC weight
            ic_weights <- aicw(c(original_ic, ic_without_current_shift))$aicweights
            ic_weight <- ic_weights[1]  # First element is the weight for the model with the shift
            
            cat(paste("IC weight for the shift is", round(ic_weight, digits=2)))
            ic_weights_df <- rbind(
              ic_weights_df,
              data.frame(
                node = shift_node_number,
                ic_with_shift = original_ic,
                ic_without_shift = ic_without_current_shift,
                delta_ic = delta_ic,
                ic_weight = ic_weight
              )
            )
          }
        } else {
          cat("No shifts were detected in the initial search, skipping IC weights calculation.\n")
          ic_weights_df <- NA
        }
      }
      
      if (uncertaintyweights_par) {
        if (length(unlist(shift_vec)) > 0) {
          cat('Calculating IC weights for initially identified shifts in parallel...\n')
          ic_weights_df <- data.frame(node = integer(), ic_with_shift = numeric(), ic_without_shift = numeric(), delta_ic = numeric(), ic_weight = numeric())
          
          # Retrieve the IC of the optimized model before uncertainty analysis
          original_ic <- if (IC == "GIC") model_with_shift_no_uncertainty$GIC$GIC else model_with_shift_no_uncertainty$BIC$BIC
          
          # Prepare a list of trees with shifts removed
          shift_removed_trees <- lapply(unlist(shift_vec), function(shift_node_number) {
            return(removeShiftFromTree(best_tree_no_uncertainty, shift_node_number))
          })
          
          # Enable parallel processing
          if (.Platform$OS.type == "unix") {
            plan(multicore, workers = num_cores)
          } else {
            plan(multisession, workers = num_cores)
          }
          # Capture additional arguments into a list
          args_list <- list(...)
          # Calculate IC weights in parallel
          ic_results <- future.apply::future_lapply(shift_removed_trees, function(tree) {
            model_function <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
            model_without_shift <- do.call(model_function, c(list(formula, tree, trait_data), args_list))
            ic_without_shift <- if (IC == "GIC") model_without_shift$GIC$GIC else model_without_shift$BIC$BIC
            delta_ic <- original_ic - ic_without_shift
            ic_weights <- aicw(c(original_ic, ic_without_shift))$aicweights
            return(c(ic_weight_withshift = ic_weights[1], ic_weight_withoutshift = ic_weights[2], delta_ic = delta_ic))
          }, future.seed = TRUE, future.scheduling = T)
          
          # Reset to sequential plan
          future::plan(future::sequential)
          
          # Add results to the dataframe
          for (i in seq_along(shift_removed_trees)) {
            shift_node_number <- unlist(shift_vec)[i]
            ic_res <- ic_results[[i]]
            ic_weights_df <- rbind(ic_weights_df, data.frame(
              node = shift_node_number,
              ic_with_shift = original_ic,
              ic_without_shift = original_ic - ic_res['delta_ic'],
              delta_ic = ic_res['delta_ic'],
              ic_weight_withshift = ic_res['ic_weight_withshift'],
              ic_weight_withoutshift = ic_res['ic_weight_withoutshift'],
              evidence_ratio = ic_res['ic_weight_withshift'] / ic_res['ic_weight_withoutshift']
            ))
          }
        } else {
          cat("No shifts were detected in the initial search, skipping IC weights calculation.\n")
          ic_weights_df <- NA
        }
      }
    } else {
      print("Only one of uncertaintyweights or uncertaintyweights_par can be set to TRUE")
    }
    
    # Print statements for the optimal configuration and delta GIC/BIC
    if (IC == "GIC") {
      cat(paste('Optimal configuration found with GIC:', round(current_best_ic, digits=2), '\n'))
      cat(paste('Global Delta GIC:', round(baseline_ic, digits=2) - round(current_best_ic, digits=2), '\n'))
    } else if (IC == "BIC") {
      cat(paste('Optimal configuration found with BIC:', round(current_best_ic, digits=2), '\n'))
      cat(paste('Global Delta BIC:', round(baseline_ic, digits=2) - round(current_best_ic, digits=2), '\n'))
    }
    
    # #Assembling the results
    # {
    #   { #taking directly from the model fit to ensure the correct tree is transferred to the results
    #     if(uncertainty) {
    #       opt_uncertainty_transformed <- model_without_shift$model$corrSt$phy
    #       opt_uncertainty_untransformed <- model_without_shift$model$corrSt$phy
    #       opt_uncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length #transfer the edge lengths
    #     }
    #     #taking directly from the model fit to ensure the correct tree is transferred to the results
    #     opt_nouncertainty_transformed <- model_with_shift_no_uncertainty$model$corrSt$phy
    #     opt_nouncertainty_untransformed <- model_with_shift_no_uncertainty$model$corrSt$phy
    #     opt_nouncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length #transfer the edge lengths
    #   }
    # 
    # 
    # # Create the main list that will always be returned
    # result_list <- list(
    #   user_input = user_input,
    #   tree_no_uncertainty_transformed = opt_nouncertainty_transformed,
    #   tree_no_uncertainty_untransformed = opt_nouncertainty_untransformed,
    #   model_no_uncertainty = model_with_shift_no_uncertainty$model,
    #   shift_nodes_no_uncertainty = shifts_no_uncertainty,
    #   optimal_ic = current_best_ic, 
    #   baseline_ic = baseline_ic,
    #   IC_used = IC,
    #   num_candidates = length(sorted_candidates),
    #   model_fit_history = model_fit_history
    # )
    # 
    # #generate the VCVs per regime from the overall model fit
    # {
    #   model_output<-result_list$model_no_uncertainty
    #   result_list$VCVs <- extractRegimeVCVs(model_output)
    # }
    # 
    # # Add the ic_weights to the list conditionally
    # if (uncertaintyweights | uncertaintyweights_par) {
    #   result_list$ic_weights <- ic_weights_df
    # }
    # 
    # # Add the uncertainty outputs to the list conditionally
    # if (uncertainty) {
    #   result_list$tree_uncertainty_transformed <- opt_uncertainty_transformed
    #   result_list$tree_uncertainty_untransformed <- opt_uncertainty_untransformed
    #   result_list$model_uncertainty <- model_without_shift$model
    #   result_list$shift_nodes_uncertainty <- shift_vec_uncertainty
    # }
    # 
    # # Add the warnings to the output list conditionally
    # if(length(warnings_list) > 0) {
    #   result_list$warnings <- warnings_list
    # }
    # 
    # }
    
    # Assembling the results
    {
      # Taking directly from the model fit to ensure the correct tree is transferred to the results
      if(uncertainty) {
        opt_uncertainty_transformed <- model_without_shift$model$corrSt$phy
        opt_uncertainty_untransformed <- model_without_shift$model$corrSt$phy
        opt_uncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length 
      }
      
      opt_nouncertainty_transformed <- model_with_shift_no_uncertainty$model$corrSt$phy
      opt_nouncertainty_untransformed <- model_with_shift_no_uncertainty$model$corrSt$phy
      opt_nouncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length 
      
      # Create the main list that will always be returned
      result_list <- list(
        user_input = user_input,
        tree_no_uncertainty_transformed = opt_nouncertainty_transformed,
        tree_no_uncertainty_untransformed = opt_nouncertainty_untransformed,
        model_no_uncertainty = model_with_shift_no_uncertainty$model,
        shift_nodes_no_uncertainty = shifts_no_uncertainty,
        optimal_ic = current_best_ic, 
        baseline_ic = baseline_ic,
        IC_used = IC,
        num_candidates = length(sorted_candidates),
        model_fit_history = model_fit_history
      )
      
      # Create the IC and acceptance matrix from the model fit history
      if (store_model_fit_history) {
        ic_acceptance_matrix <- do.call(rbind, lapply(model_fit_history, function(x) {
          if (is.null(x$model)) {
            c(NA, x$accepted)
          } else {
            if (IC == "GIC") c(x$model$GIC$GIC, x$accepted) else c(x$model$BIC$BIC, x$accepted)
          }
        }))
        
        result_list$model_fit_history <- list(
          fits = model_fit_history,
          ic_acceptance_matrix = ic_acceptance_matrix
        )
      }
      
      # Generate the VCVs per regime from the overall model fit
      model_output <- result_list$model_no_uncertainty
      result_list$VCVs <- extractRegimeVCVs(model_output)
      
      # Add the ic_weights to the list conditionally
      if (uncertaintyweights | uncertaintyweights_par) {
        result_list$ic_weights <- ic_weights_df
      }
      
      # Add the uncertainty outputs to the list conditionally
      if (uncertainty) {
        result_list$tree_uncertainty_transformed <- opt_uncertainty_transformed
        result_list$tree_uncertainty_untransformed <- opt_uncertainty_untransformed
        result_list$model_uncertainty <- model_without_shift$model
        result_list$shift_nodes_uncertainty <- shift_vec_uncertainty
      }
      
      # Add the warnings to the output list conditionally
      if(length(warnings_list) > 0) {
        result_list$warnings <- warnings_list
      }
      
    }
    
    return(result_list)
    
  }
#... is optional arguments to be passed to mvgls


#' Search for the Optimal Shift Configuration Using `mclapply` for Parallelization
#'
#' Identical to `searchOptimalConfiguration()`, this function searches for the optimal configuration
#' of evolutionary shifts in a phylogenetic model, but uses `parallel::mclapply()` instead of `future_lapply()`
#' for parallel model fitting. This version is recommended only for Unix-based systems (Linux/macOS),
#' as `mclapply()` is not supported on Windows.
searchOptimalConfiguration.mclapply <-
  function(baseline_tree,
           trait_data,
           formula = 'trait_data~1',
           min_descendant_tips,
           num_cores = 2,
           ic_uncertainty_threshold = 1.0,
           shift_acceptance_threshold = 1.0,
           uncertainty = F,
           uncertaintyweights = F,
           uncertaintyweights_par = F,
           postorder_traversal = F,
           plot = T,
           IC = 'GIC', store_model_fit_history = TRUE, ...) {
    
    # Capture user input
    user_input <- as.list(match.call())
    
    #generate initial set of painted candidate trees with shifts at each sub-node
    cat('Generating candidate shift models...\n')
    candidate_trees <- generatePaintedTrees(baseline_tree, min_descendant_tips)
    candidate_trees_shifts <- candidate_trees[-1]
    
    #fit the initial baseline model to the baseline tree with a global regime (state=0)
    cat('Fitting baseline model...\n')
    
    #select which information criterion to use
    if (IC != "GIC" && IC != "BIC") {
      stop("IC must be GIC or BIC")
    }
    if(IC=="GIC"){
      baseline_model <- fitMvglsAndExtractGIC.formula(formula, candidate_trees[[1]], trait_data, ...)
      baseline_ic <- baseline_model$GIC$GIC
    }
    if(IC=="BIC"){
      baseline_model <- fitMvglsAndExtractBIC.formula(formula, candidate_trees[[1]], trait_data, ...)
      baseline_ic <- baseline_model$BIC$BIC
    }
    cat(paste('Baseline ', IC, ':', round(baseline_ic, digits = 2), '\n'))
    
    #evaluate all of the candidate trees under GIC or BIC
    # Capture additional arguments into a list
    args_list <- list(...)
    
    cat('Fitting sub-models in parallel...\n')
    
    # Temporarily assign variables to the global environment
    assign("formula", formula, envir = .GlobalEnv)
    assign("trait_data", trait_data, envir = .GlobalEnv)
    assign("args_list", args_list, envir = .GlobalEnv)
    assign("IC", IC, envir = .GlobalEnv)
    
    #plan(multisession, workers = num_cores)
    candidate_results <-
      mclapply(candidate_trees_shifts, function(tree) {
        if (IC == "GIC") {
          do.call(fitMvglsAndExtractGIC.formula, c(list(formula, tree, trait_data), args_list))
        } else if (IC == "BIC") {
          do.call(fitMvglsAndExtractBIC.formula, c(list(formula, tree, trait_data), args_list))
        }
      }, mc.cores=num_cores)
    #plan(sequential)
    
    # Clean up the global environment by removing temporary variables
    rm(formula, trait_data, args_list, IC, envir = .GlobalEnv)
    
    #generate the delta IC lists
    if (IC == "GIC"){
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$GIC$GIC)
    } else if (IC == "BIC") {
      delta_ic_list <- sapply(candidate_results, function(res) baseline_ic - res$BIC$BIC)
    }
    
    cat('Sorting and evaluating shifts...\n')
    sorted_candidates <- candidate_trees_shifts[order(delta_ic_list, decreasing=T)]
    current_best_tree <- baseline_tree
    current_best_ic <- baseline_ic
    shift_id <- 0
    
    #plot current shift tree
    if(plot==T){
      plotSimmap(current_best_tree, ftype='off')
    }
    
    #in case the user wants postorder traversal for shift searching (default off)
    if(postorder_traversal){
      print('Candidate shifts are sorted in Postorder')
      sorted_candidates <- (candidate_trees_shifts)
      current_best_tree <- baseline_tree
      current_best_ic <- baseline_ic
      shift_id <- 0
    }
    
    shift_vec <- list() #initialize shift_vec
    model_with_shift_no_uncertainty<-NULL #initialize output
    best_tree_no_uncertainty<-NULL #initialize output
    # Initialize the list to collect warning messages
    warnings_list <- list()
    model_fit_history<-list() #new object to capture the full history of the search
    
    #Run the primary shift configuration search
    
    # for (i in seq_along(sorted_candidates)) {
    #   shift_node_name <- names(sorted_candidates)[i]
    #   shift_node_number <- as.integer(sub("Node ", "", shift_node_name))
    #   percent_complete <- round((i / length(sorted_candidates)) * 100, 2)
    #   print(paste('Evaluating shift at node:', shift_node_number, '-', percent_complete, '% complete'))
    #   
    #   add_shift_result <- addShiftToModel(current_best_tree, shift_node_number, shift_id)
    #   shifted_tree <- add_shift_result$tree
    #   shift_id <- add_shift_result$shift_id
    #   
    #   if(plot==T){
    #     nodelabels(text = shift_id, node=shift_node_number)
    #   }
    #   
    #   # tryCatch({
    #   #   if (IC == "GIC") {
    #   #     model_with_shift <- fitMvglsAndExtractGIC.formula(formula, shifted_tree, trait_data, ...)
    #   #     new_ic <- model_with_shift$GIC$GIC
    #   #   } else if (IC == "BIC") {
    #   #     model_with_shift <- fitMvglsAndExtractBIC.formula(formula, shifted_tree, trait_data, ...)
    #   #     new_ic <- model_with_shift$BIC$BIC
    #   #   }
    #   #   
    #   #   delta_ic <- current_best_ic - new_ic
    #   #   
    #   #   if (delta_ic >= shift_acceptance_threshold) {
    #   #     current_best_tree <- shifted_tree
    #   #     current_best_ic <- new_ic
    #   #     print(paste('Shift at node', shift_node_number, 'accepted. Updated', IC, ':', current_best_ic, 'Delta', IC, ':', delta_ic))
    #   #     
    #   #     shift_vec[[length(shift_vec) + 1]] <- shift_node_number
    #   #     
    #   #     best_tree_no_uncertainty <- current_best_tree
    #   #     model_with_shift_no_uncertainty <- model_with_shift
    #   #   } else {
    #   #     print(paste('Shift at node', shift_node_number, 'rejected. Delta', IC, ':', delta_ic, 'is less than threshold:', shift_acceptance_threshold))
    #   #   }
    #   # }, error = function(e) {
    #   #   warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
    #   #   warning(warning_message)
    #   #   warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #   # }, warning = function(w) {
    #   #   warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
    #   #   warning(warning_message)
    #   #   warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #   # })
    #   tryCatch({
    #     if (IC == "GIC") {
    #       model_with_shift <- withCallingHandlers(
    #         fitMvglsAndExtractGIC.formula(formula, shifted_tree, trait_data, ...),
    #         warning = function(w) {
    #           # Log the warning but continue execution
    #           warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
    #           warning(warning_message)
    #           warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #           invokeRestart("muffleWarning")  # Prevent the warning from propagating further
    #         }
    #       )
    #       new_ic <- model_with_shift$GIC$GIC
    #     } else if (IC == "BIC") {
    #       model_with_shift <- withCallingHandlers(
    #         fitMvglsAndExtractBIC.formula(formula, shifted_tree, trait_data, ...),
    #         warning = function(w) {
    #           # Log the warning but continue execution
    #           warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
    #           warning(warning_message)
    #           warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #           invokeRestart("muffleWarning")  # Prevent the warning from propagating further
    #         }
    #       )
    #       new_ic <- model_with_shift$BIC$BIC
    #     }
    #     
    #     # Calculate delta IC
    #     delta_ic <- current_best_ic - new_ic
    #     
    #     # Decision logic
    #     if (delta_ic >= shift_acceptance_threshold) {
    #       current_best_tree <- shifted_tree
    #       current_best_ic <- new_ic
    #       print(paste('Shift at node', shift_node_number, 'accepted. Updated', IC, ':', current_best_ic, 'Delta', IC, ':', delta_ic))
    #       
    #       shift_vec[[length(shift_vec) + 1]] <- shift_node_number
    #       
    #       best_tree_no_uncertainty <- current_best_tree
    #       model_with_shift_no_uncertainty <- model_with_shift
    #     } else {
    #       print(paste('Shift at node', shift_node_number, 'rejected. Delta', IC, ':', delta_ic, 'is less than threshold:', shift_acceptance_threshold))
    #     }
    #   }, error = function(e) {
    #     # Handle errors
    #     warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
    #     warning(warning_message)
    #     warnings_list[[length(warnings_list) + 1]] <<- warning_message
    #   })
    #   
    #   
    #   if(plot==T){
    #     colorvec <- setNames(object = c('black', rainbow(length(unique(getStates(shifted_tree, type = 'both')))-1)),
    #                          nm = sort(as.numeric(unique(getStates(shifted_tree, type = 'both')))))
    #     plotSimmap(current_best_tree, colors=colorvec, fsize=0.0001)
    #   }
    # }
    
    for (i in seq_along(sorted_candidates)) {
      shift_node_name <- names(sorted_candidates)[i]
      shift_node_number <- as.integer(sub("Node ", "", shift_node_name))
      percent_complete <- round((i / length(sorted_candidates)) * 100, 2)
      cat(paste('Evaluating shift at node:', shift_node_number, '-', percent_complete, '% complete', '\n'))
      
      add_shift_result <- addShiftToModel(current_best_tree, shift_node_number, shift_id)
      shifted_tree <- add_shift_result$tree
      shift_id <- add_shift_result$shift_id
      
      if(plot==T){
        nodelabels(text = shift_id, node=shift_node_number)
      }
      
      tryCatch({
        if (IC == "GIC") {
          model_with_shift <- withCallingHandlers(
            fitMvglsAndExtractGIC.formula(formula, shifted_tree, trait_data, ...),
            warning = function(w) {
              warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
              warning(warning_message)
              warnings_list[[length(warnings_list) + 1]] <<- warning_message
              invokeRestart("muffleWarning") 
            }
          )
          new_ic <- model_with_shift$GIC$GIC
        } else if (IC == "BIC") {
          model_with_shift <- withCallingHandlers(
            fitMvglsAndExtractBIC.formula(formula, shifted_tree, trait_data, ...),
            warning = function(w) {
              warning_message <- paste("Warning in evaluating shift at node", shift_node_number, ":", w$message)
              warning(warning_message)
              warnings_list[[length(warnings_list) + 1]] <<- warning_message
              invokeRestart("muffleWarning") 
            }
          )
          new_ic <- model_with_shift$BIC$BIC
        }
        
        # Calculate delta IC
        delta_ic <- current_best_ic - new_ic
        
        # Store model fit and acceptance status (including delta_ic)
        if (store_model_fit_history) {
          model_fit_history[[length(model_fit_history) + 1]] <- list(
            model = model_with_shift,
            accepted = delta_ic >= shift_acceptance_threshold,
            delta_ic = delta_ic 
          )
        }
        
        # Decision logic (unchanged)
        if (delta_ic >= shift_acceptance_threshold) {
          current_best_tree <- shifted_tree
          current_best_ic <- new_ic
          cat(paste('Shift at node', shift_node_number, 'accepted. Updated', IC, ':', round(current_best_ic, digits = 2), 'Delta', IC, ':', round(delta_ic, digits = 2), '\n')) 
          
          shift_vec[[length(shift_vec) + 1]] <- shift_node_number
          
          best_tree_no_uncertainty <- current_best_tree
          model_with_shift_no_uncertainty <- model_with_shift
        } else {
          cat(paste('Shift at node', shift_node_number, 'rejected. Delta', IC, ':', round(delta_ic, digits = 2), 'is less than threshold:', shift_acceptance_threshold, '\n')) 
        }
      }, error = function(e) {
        # Handle errors (unchanged)
        warning_message <- paste("Error in evaluating shift at node", shift_node_number, ":", e$message)
        warning(warning_message)
        warnings_list[[length(warnings_list) + 1]] <<- warning_message
        
        # Also store the error in the model fit history
        if (store_model_fit_history) {
          model_fit_history[[length(model_fit_history) + 1]] <- list(
            model = NULL,
            accepted = FALSE,
            delta_ic = NA,
            error = e$message
          )
        }
      })
      
      if(plot==T){
        colorvec <- setNames(object = c('black', rainbow(length(unique(getStates(shifted_tree, type = 'both')))-1)),
                             nm = sort(as.numeric(unique(getStates(shifted_tree, type = 'both')))))
        plotSimmap(current_best_tree, colors=colorvec, ftype='off')
      }
    }
    
    #print(paste(shift_vec))
    shifts_no_uncertainty<-unlist(shift_vec)
    cat(paste("Shifts detected at nodes:", paste(shift_vec, collapse = ", "), '\n'))
    
    # If activated, this section removes shifts after re-evaluating
    model_without_shift <- NULL # Initialize as NULL
    if (uncertainty) {
      cat('Post-search re-evaluation to reduce overfitting...')
      shift_nodes <- unlist(shift_vec)
      print(paste("Re-evaluating nodes", shift_nodes))
      shift_vec_uncertainty <- shift_nodes # Tracker (to remove shifts from)
      
      root_node <- Ntip(baseline_tree) + 1
      for (shift_node_number in shift_nodes) {
        if (shift_node_number != root_node) {
          cat(paste('Re-evaluating shift at node:', shift_node_number))
          tree_without_shift <- removeShiftFromTree(current_best_tree, shift_node_number)
          
          tryCatch({
            # Evaluate the model without the current shift using the selected IC
            if (IC == "GIC") {
              temp_model <- fitMvglsAndExtractGIC.formula(formula, tree_without_shift, trait_data, ...)
              ic_without_shift <- temp_model$GIC$GIC
            } else if (IC == "BIC") {
              temp_model <- fitMvglsAndExtractBIC.formula(formula, tree_without_shift, trait_data, ...)
              ic_without_shift <- temp_model$BIC$BIC
            }
            
            if (abs(current_best_ic - ic_without_shift) <= ic_uncertainty_threshold) {
              current_best_tree <- tree_without_shift
              current_best_ic <- ic_without_shift
              cat(paste('Shift at node', shift_node_number, 'removed. Updated', IC, ':', round(current_best_ic, digits=2)))
              shift_vec_uncertainty <- shift_vec_uncertainty[!shift_vec_uncertainty == shift_node_number]
              model_without_shift <- temp_model # Update only if the condition is met
            }
          }, error = function(e) {
            warning(paste("Error in re-evaluating shift at node", shift_node_number, ":", e$message))
          })
        } else {
          cat(paste('Skipping re-evaluation for the root node:', shift_node_number))
        }
      }
    }
    
    # New Section for Calculating Information Criterion Weights Post Optimization
    ic_weights_df <- NA  # Initialize the results vector
    if (xor(uncertaintyweights, uncertaintyweights_par)) {
      if (uncertaintyweights) {
        if (length(unlist(shift_vec)) > 0) {
          cat('Calculating IC weights for initially identified shifts...\n')
          ic_weights_df <- data.frame(node = integer(), ic_with_shift = numeric(), ic_without_shift = numeric(), delta_ic = numeric(), ic_weight = numeric())
          
          # Retrieve the IC of the optimized model before uncertainty analysis
          original_ic <- if (IC == "GIC") model_with_shift_no_uncertainty$GIC$GIC else model_with_shift_no_uncertainty$BIC$BIC
          
          cat(paste("Considering", length(shift_vec), "shifts in the candidate set",'\n'))
          cat(paste("There are", length(unique(getStates(best_tree_no_uncertainty, type = 'both'))) - 1, "shifts in the mapped tree",'\n'))
          
          for (shift_node_number in unlist(shift_vec)) {
            cat(paste('Re-estimating model without shift at node:', shift_node_number, '\n'))
            
            # Remove the shift temporarily from best_tree_no_uncertainty and re-estimate the IC
            tree_without_current_shift <- removeShiftFromTree(best_tree_no_uncertainty, shift_node_number)
            model_without_current_shift_function <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
            model_without_current_shift <- model_without_current_shift_function(formula, tree_without_current_shift, trait_data, ...)
            ic_without_current_shift <- if (IC == "GIC") model_without_current_shift$GIC$GIC else model_without_current_shift$BIC$BIC
            
            # Calculate the difference in IC
            delta_ic <- original_ic - ic_without_current_shift
            
            # Use aicw function to calculate the IC weight
            ic_weights <- aicw(c(original_ic, ic_without_current_shift))$aicweights
            ic_weight <- ic_weights[1]  # First element is the weight for the model with the shift
            
            cat(paste("IC weight for the shift is", round(ic_weight, digits=2)))
            ic_weights_df <- rbind(
              ic_weights_df,
              data.frame(
                node = shift_node_number,
                ic_with_shift = original_ic,
                ic_without_shift = ic_without_current_shift,
                delta_ic = delta_ic,
                ic_weight = ic_weight
              )
            )
          }
        } else {
          cat("No shifts were detected in the initial search, skipping IC weights calculation.\n")
          ic_weights_df <- NA
        }
      }
      
      if (uncertaintyweights_par) {
        if (length(unlist(shift_vec)) > 0) {
          cat('Calculating IC weights for initially identified shifts in parallel...\n')
          ic_weights_df <- data.frame(node = integer(), ic_with_shift = numeric(), ic_without_shift = numeric(), delta_ic = numeric(), ic_weight = numeric())
          
          # Retrieve the IC of the optimized model before uncertainty analysis
          original_ic <- if (IC == "GIC") model_with_shift_no_uncertainty$GIC$GIC else model_with_shift_no_uncertainty$BIC$BIC
          
          # Prepare a list of trees with shifts removed
          shift_removed_trees <- lapply(unlist(shift_vec), function(shift_node_number) {
            return(removeShiftFromTree(best_tree_no_uncertainty, shift_node_number))
          })
          
          # Enable parallel processing
          # Capture additional arguments into a list
          args_list <- list(...)
          
          # Temporary assignment to the global environment for this mclapply section
          assign("formula", formula, envir = .GlobalEnv)
          assign("IC", IC, envir = .GlobalEnv)
          assign("trait_data", trait_data, envir = .GlobalEnv)
          assign("args_list", args_list, envir = .GlobalEnv)
          assign("original_ic", original_ic, envir = .GlobalEnv)
          
          # Calculate IC weights in parallel
          ic_results <- mclapply(shift_removed_trees, function(tree) {
            model_function <- if (IC == "GIC") fitMvglsAndExtractGIC.formula else fitMvglsAndExtractBIC.formula
            model_without_shift <- do.call(model_function, c(list(formula, tree, trait_data), args_list))
            ic_without_shift <- if (IC == "GIC") model_without_shift$GIC$GIC else model_without_shift$BIC$BIC
            delta_ic <- original_ic - ic_without_shift
            ic_weights <- aicw(c(original_ic, ic_without_shift))$aicweights
            return(c(ic_weight_withshift = ic_weights[1], ic_weight_withoutshift = ic_weights[2], delta_ic = delta_ic))
          }, mc.cores=num_cores)
          
          # Clean up the global environment by removing temporary variables
          rm(formula, IC, trait_data, args_list, original_ic, envir = .GlobalEnv)
          
          # Add results to the dataframe
          for (i in seq_along(shift_removed_trees)) {
            shift_node_number <- unlist(shift_vec)[i]
            ic_res <- ic_results[[i]]
            ic_weights_df <- rbind(ic_weights_df, data.frame(
              node = shift_node_number,
              ic_with_shift = original_ic,
              ic_without_shift = original_ic - ic_res['delta_ic'],
              delta_ic = ic_res['delta_ic'],
              ic_weight_withshift = ic_res['ic_weight_withshift'],
              ic_weight_withoutshift = ic_res['ic_weight_withoutshift'],
              evidence_ratio = ic_res['ic_weight_withshift'] / ic_res['ic_weight_withoutshift']
            ))
          }
        } else {
          cat("No shifts were detected in the initial search, skipping IC weights calculation.\n")
          ic_weights_df <- NA
        }
      }
    } else {
      print("Only one of uncertaintyweights or uncertaintyweights_par can be set to TRUE")
    }
    
    # Print statements for the optimal configuration and delta GIC/BIC
    if (IC == "GIC") {
      cat(paste('Optimal configuration found with GIC:', round(current_best_ic, digits=2), '\n'))
      cat(paste('Global Delta GIC:', round(baseline_ic, digits=2) - round(current_best_ic, digits=2), '\n'))
    } else if (IC == "BIC") {
      cat(paste('Optimal configuration found with BIC:', round(current_best_ic, digits=2), '\n'))
      cat(paste('Global Delta BIC:', round(baseline_ic, digits=2) - round(current_best_ic, digits=2), '\n'))
    }
    
    # #Assembling the results
    # {
    #   { #taking directly from the model fit to ensure the correct tree is transferred to the results
    #     if(uncertainty) {
    #       opt_uncertainty_transformed <- model_without_shift$model$corrSt$phy
    #       opt_uncertainty_untransformed <- model_without_shift$model$corrSt$phy
    #       opt_uncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length #transfer the edge lengths
    #     }
    #     #taking directly from the model fit to ensure the correct tree is transferred to the results
    #     opt_nouncertainty_transformed <- model_with_shift_no_uncertainty$model$corrSt$phy
    #     opt_nouncertainty_untransformed <- model_with_shift_no_uncertainty$model$corrSt$phy
    #     opt_nouncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length #transfer the edge lengths
    #   }
    # 
    # 
    # # Create the main list that will always be returned
    # result_list <- list(
    #   user_input = user_input,
    #   tree_no_uncertainty_transformed = opt_nouncertainty_transformed,
    #   tree_no_uncertainty_untransformed = opt_nouncertainty_untransformed,
    #   model_no_uncertainty = model_with_shift_no_uncertainty$model,
    #   shift_nodes_no_uncertainty = shifts_no_uncertainty,
    #   optimal_ic = current_best_ic, 
    #   baseline_ic = baseline_ic,
    #   IC_used = IC,
    #   num_candidates = length(sorted_candidates),
    #   model_fit_history = model_fit_history
    # )
    # 
    # #generate the VCVs per regime from the overall model fit
    # {
    #   model_output<-result_list$model_no_uncertainty
    #   result_list$VCVs <- extractRegimeVCVs(model_output)
    # }
    # 
    # # Add the ic_weights to the list conditionally
    # if (uncertaintyweights | uncertaintyweights_par) {
    #   result_list$ic_weights <- ic_weights_df
    # }
    # 
    # # Add the uncertainty outputs to the list conditionally
    # if (uncertainty) {
    #   result_list$tree_uncertainty_transformed <- opt_uncertainty_transformed
    #   result_list$tree_uncertainty_untransformed <- opt_uncertainty_untransformed
    #   result_list$model_uncertainty <- model_without_shift$model
    #   result_list$shift_nodes_uncertainty <- shift_vec_uncertainty
    # }
    # 
    # # Add the warnings to the output list conditionally
    # if(length(warnings_list) > 0) {
    #   result_list$warnings <- warnings_list
    # }
    # 
    # }
    
    # Assembling the results
    {
      # Taking directly from the model fit to ensure the correct tree is transferred to the results
      if(uncertainty) {
        opt_uncertainty_transformed <- model_without_shift$model$corrSt$phy
        opt_uncertainty_untransformed <- model_without_shift$model$corrSt$phy
        opt_uncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length 
      }
      
      opt_nouncertainty_transformed <- model_with_shift_no_uncertainty$model$corrSt$phy
      opt_nouncertainty_untransformed <- model_with_shift_no_uncertainty$model$corrSt$phy
      opt_nouncertainty_untransformed$edge.length <- best_tree_no_uncertainty$edge.length 
      
      # Create the main list that will always be returned
      result_list <- list(
        user_input = user_input,
        tree_no_uncertainty_transformed = opt_nouncertainty_transformed,
        tree_no_uncertainty_untransformed = opt_nouncertainty_untransformed,
        model_no_uncertainty = model_with_shift_no_uncertainty$model,
        shift_nodes_no_uncertainty = shifts_no_uncertainty,
        optimal_ic = current_best_ic, 
        baseline_ic = baseline_ic,
        IC_used = IC,
        num_candidates = length(sorted_candidates),
        model_fit_history = model_fit_history
      )
      
      # Create the IC and acceptance matrix from the model fit history
      if (store_model_fit_history) {
        ic_acceptance_matrix <- do.call(rbind, lapply(model_fit_history, function(x) {
          if (is.null(x$model)) {
            c(NA, x$accepted)
          } else {
            if (IC == "GIC") {
              c(x$model$GIC$GIC, x$accepted)
            } else {
              c(x$model$BIC$BIC, x$accepted)
            }
          }
        }))
        
        result_list$model_fit_history <- list(
          fits = model_fit_history, 
          ic_acceptance_matrix = ic_acceptance_matrix
        )
      }
      
      # Add the matrix to the model_fit_history
      result_list$model_fit_history <- list(
        fits = result_list$model_fit_history, 
        ic_acceptance_matrix = ic_acceptance_matrix
      )
      
      # Generate the VCVs per regime from the overall model fit
      model_output <- result_list$model_no_uncertainty
      result_list$VCVs <- extractRegimeVCVs(model_output)
      
      # Add the ic_weights to the list conditionally
      if (uncertaintyweights | uncertaintyweights_par) {
        result_list$ic_weights <- ic_weights_df
      }
      
      # Add the uncertainty outputs to the list conditionally
      if (uncertainty) {
        result_list$tree_uncertainty_transformed <- opt_uncertainty_transformed
        result_list$tree_uncertainty_untransformed <- opt_uncertainty_untransformed
        result_list$model_uncertainty <- model_without_shift$model
        result_list$shift_nodes_uncertainty <- shift_vec_uncertainty
      }
      
      # Add the warnings to the output list conditionally
      if(length(warnings_list) > 0) {
        result_list$warnings <- warnings_list
      }
      
    }
    
    return(result_list)
    
  }
#... is optional arguments to be passed to mvgls
}

## helper functions for downstream analysis below ##

#' Extract Maximum Ages of Regimes in a SIMMAP Tree
#'
#' This function identifies the oldest (most ancestral) node associated with each 
#' unique regime (state) in a SIMMAP-style phylogenetic tree and returns a summary 
#' of the maximum ages per regime.
extractMaxAgeOfRegimes <- function(simmap_tree) {
  if (!inherits(simmap_tree, "phylo")) {
    stop("The input must be a 'phylo' object.")
  }
  
  # Convert simmap tree to phylo tree for compatibility with tree.age
  phylo_tree <- as.phylo(simmap_tree)
  
  # Get the ages of all nodes
  ages_info <- dispRity::tree.age(phylo_tree, order = 'past')
  # Filter out the tip ages
  ages_info <- ages_info[!ages_info$ages == 0,]
  
  # Get the states for all nodes
  states <- phytools::getStates(simmap_tree, type = 'both')
  
  # Initialize an empty list to store results
  results <- list()
  
  # Iterate over unique states to find the maximum age node
  for (state in unique(states)) {
    nodes_in_state <- names(states[states == state])
    ages_in_state <- ages_info$ages[ages_info$elements %in% nodes_in_state]
    
    if (length(ages_in_state) > 0) {
      max_age <- max(ages_in_state, na.rm = TRUE)
      max_age_index <- which.max(ages_in_state)
      max_age_node <- nodes_in_state[max_age_index]
    } else {
      max_age <- NA
      max_age_node <- NA
    }
    
    results[[state]] <- c(state = state, max_age_node = max_age_node, max_age = max_age)
  }
  
  # Construct the data frame from the list of results
  max_ages_df <- do.call(rbind, results)
  max_ages_df <- as.data.frame(max_ages_df, stringsAsFactors = FALSE)
  max_ages_df$max_age_node <- as.integer(max_ages_df$max_age_node)
  max_ages_df$max_age <- as.numeric(max_ages_df$max_age)
  
  return(max_ages_df)
}

#functions to slice the tree into time bins and count the number of lineages in each bin
{
  #' Calculate Average Lineage Counts in Sliding Time Bins
  #'
  #' Uses an LTT (lineages through time) plot to count the number of lineages 
  #' within sliding time bins across a phylogenetic tree.
  calculateAverageLineagesFromTree <- function(tree, bin_size, slide_step) {
    if (!is.numeric(bin_size) || bin_size <= 0 || !is.numeric(slide_step) || slide_step <= 0) {
      stop("bin_size and slide_step must be positive numbers.")
    }
    
    ltt_data <- ape::ltt.plot.coords(tree)
    ltt_data[, 1] <- abs(ltt_data[, 1])
    max_time <- max(ltt_data[, 1])
    
    bin_starts <- seq(0, max_time - bin_size, by = slide_step)
    bin_ends <- bin_starts + bin_size
    bin_midpoints <- bin_starts + bin_size / 2
    
    average_lineages <- data.frame(
      bin_start = bin_starts,
      bin_end = bin_ends,
      bin_midpoint = bin_midpoints,
      average = numeric(length(bin_starts))
    )
    
    for (i in seq_along(bin_starts)) {
      bin_data <- ltt_data[ltt_data[, 1] >= bin_starts[i] & ltt_data[, 1] < bin_ends[i], ]
      if (length(bin_data[, 2]) > 0) {
        average_lineages$average[i] <- max(max(bin_data[, 2]), 2)
      } else {
        average_lineages$average[i] <- 2
      }
    }
    
    return(average_lineages)
  }
  
  
  #' Count Values Within Sliding Time Bins
  #'
  #' Counts how many values in a numeric vector (e.g., node ages or event times) fall 
  #' within each of a set of sliding time bins.
  countAgesInBins <- function(ages, bin_size, slide_step) {
    if (!is.numeric(bin_size) || bin_size <= 0 || !is.numeric(slide_step) || slide_step <= 0) {
      stop("bin_size and slide_step must be positive numbers.")
    }
    
    if (!is.numeric(ages)) {
      stop("Ages must be a numeric vector.")
    }
    
    max_age <- max(ages)
    bin_starts <- seq(0, max_age - bin_size, by = slide_step)
    bin_ends <- bin_starts + bin_size
    bin_midpoints <- bin_starts + bin_size / 2
    
    age_counts <- data.frame(
      bin_start = bin_starts,
      bin_end = bin_ends,
      bin_midpoint = bin_midpoints,
      count = numeric(length(bin_starts))
    )
    
    for (i in seq_along(bin_starts)) {
      age_counts$count[i] <- sum(ages >= bin_starts[i] & ages < bin_ends[i])
    }
    
    return(age_counts)
  }
  
  #' Calculate Total Branch Length per Time Bin Using Tree Slices
  #'
  #' Slices a phylogenetic tree at sliding time intervals and computes the total 
  #' branch length within each bin using overlapping tree fragments.
  calculateTotalBranchLengthWithTreeSlice <- function(tree, bin_size, slide_step, cores = parallel::detectCores() - 1) {
    # Convert the tree to a phylo object
    tree <- as.phylo(tree)
    
    # Calculate maximum height of the tree
    tree_height <- max(ape::node.depth.edgelength(tree))
    
    # Generate bin start and end points
    bin_starts <- seq(0, tree_height - bin_size, by = slide_step)
    bin_ends <- bin_starts + bin_size
    bin_midpoints <- (bin_starts + bin_ends) / 2
    
    # Helper function to add branch positions
    add_branch_positions <- function(tree) {
      node_depths <- ape::node.depth.edgelength(tree)
      branches <- data.frame(
        X1 = tree$edge[, 1],
        X2 = tree$edge[, 2],
        start = node_depths[tree$edge[, 1]],
        end = node_depths[tree$edge[, 1]] + tree$edge.length,
        length = tree$edge.length
      )
      return(branches)
    }
    
    # Helper function to calculate overlap contributions
    calculate_contributions <- function(branches, slice_start, slice_end) {
      sapply(seq_len(nrow(branches)), function(i) {
        branch <- branches[i, ]
        overlap <- max(0, min(slice_end, branch$end) - max(slice_start, branch$start))
        return(overlap)
      })
    }
    
    # Parallel computation for each bin
    calculate_bin_total <- function(slice_start, slice_end) {
      # Slice the tree rootwards at the end of the bin
      sliced_tree <- tryCatch(
        phytools::treeSlice(tree, slice = slice_end, trivial = TRUE, orientation = "rootwards"),
        error = function(e) {
          cat("Error during slicing at", slice_end, ":", e$message, "\n")
          return(NULL)
        }
      )
      
      if (is.null(sliced_tree)) {
        return(NA)
      }
      
      # Add branch positions
      branches <- add_branch_positions(sliced_tree)
      
      # Calculate contributions for the current bin
      contributions <- calculate_contributions(branches, slice_start, slice_end)
      
      # Sum up the contributions
      return(sum(contributions))
    }
    
    # Perform parallel slicing and length computation
    cat("Starting parallel computation for bins...\n")
    total_branch_lengths <- pbmcapply::pbmcmapply(
      calculate_bin_total,
      slice_start = bin_starts,
      slice_end = bin_ends,
      SIMPLIFY = TRUE,
      mc.cores = cores
    )
    
    # Construct the results data frame
    result <- data.frame(
      bin_start = bin_starts,
      bin_end = bin_ends,
      bin_midpoint = bin_midpoints,
      total_branch_length = total_branch_lengths
    )
    
    return(result)
  }
  
  #' Calculate Total Branch Length in Sliding Time Bins
  #'
  #' Calculates the amount of phylogenetic branch length in each sliding time bin 
  #' by slicing the tree rootward and computing cumulative differences.
  calculateTotalBranchLengthInBins <- function(tree, bin_size, slide_step, cores = 2) {
    # Convert tree to phylo object
    tree <- as.phylo(tree)
    
    # Calculate the maximum height of the tree (time from root to tips)
    tree_height <- max(ape::node.depth.edgelength(tree))
    
    # Define bin starts and ends
    bin_starts <- seq(0, tree_height - bin_size, by = slide_step)  # Start from 0
    bin_ends <- bin_starts + bin_size                             # End at bin_size
    bin_midpoints <- bin_starts + bin_size / 2                    # Midpoint of bins
    
    # Function to calculate cumulative branch length from the root to a given slice point
    calculate_rootward_length <- function(slice_point) {
      sliced_tree <- tryCatch(
        phytools::treeSlice(tree, slice = slice_point, trivial = TRUE, orientation = "rootwards"),
        error = function(e) {
          cat("Error during slicing at", slice_point, ":", e$message, "\n")
          return(0)
        }
      )
      # Sum branch lengths of the sliced tree(s)
      if (inherits(sliced_tree, "phylo")) {
        return(sum(sliced_tree$edge.length))
      } else if (inherits(sliced_tree, "multiPhylo")) {
        return(sum(sapply(sliced_tree, function(subtree) sum(subtree$edge.length))))
      } else {
        return(0)
      }
    }
    
    # Generate slice points for all bins
    slice_points <- unique(c(bin_starts, bin_ends))
    slice_points <- sort(slice_points[slice_points <= tree_height])  # Ensure valid slice points
    
    # Compute rootward branch lengths at all slice points in parallel
    rootward_branch_lengths <- pbmclapply(slice_points, calculate_rootward_length, mc.cores = cores)
    rootward_branch_lengths <- unlist(rootward_branch_lengths)
    
    # Map slice points to their corresponding branch lengths
    rootward_length_map <- setNames(rootward_branch_lengths, slice_points)
    
    # Compute total branch length for each bin
    total_branch_length_per_bin <- sapply(seq_along(bin_starts), function(i) {
      if (bin_starts[i] == 0) {
        # For the first bin: use only the cumulative length at bin_end
        return(rootward_length_map[as.character(bin_ends[i])])
      } else {
        # For other bins: compute the difference between bin_end and bin_start
        start_length <- rootward_length_map[as.character(bin_starts[i])]
        end_length <- rootward_length_map[as.character(bin_ends[i])]
        return(end_length - start_length)
      }
    })
    
    # Reverse the total branch lengths to reflect the descending order from tips to root
    total_branch_length_per_bin <- rev(total_branch_length_per_bin)
    
    # Create the result data frame
    result <- data.frame(
      bin_start = bin_starts,
      bin_end = bin_ends,
      bin_midpoint = bin_midpoints,
      total_branch_length = total_branch_length_per_bin,  # Descending total branch lengths
      tree_age_at_bin_start = bin_starts,  # Starts at 0 and increases
      row.names = seq_len(length(bin_starts))  # Sequential row names
    )
    
    return(result)
  }
  
}


#' Integrative mvgls Fitting on Regime-Defined Subtrees
#'
#' This function takes a SIMMAP tree (typically output from `searchOptimalConfiguration`),
#' identifies unique regimes based on mapped tip states (`getStates(type = "tips")`),
#' and fits separate `mvgls` models for each regime.
#'
#' Each subtree (corresponding to a regime) is extracted and passed to `mvgls()` along with
#' a matching subset of the trait data. The fitting can be run in parallel using `future_lapply`.
{
  
  #' Fit mvgls Models on Subtrees Defined by Tip Regimes
  #'
  #' Fits separate `mvgls` models to subtrees extracted from a SIMMAP tree,
  #' where each subtree contains tips with a unique regime (state).
  fitMvglsOnSubtrees <- function(user_formula, simmap_tree, input_data, model = 'BM', ...) {
    # Initialize the list to collect mvgls models with placeholders for each state
    unique_states <- unique(getStates(simmap_tree, type='tips'))
    mvgls_models <- setNames(vector("list", length(unique_states)), unique_states)
    
    # Fit models for each state
    for (state in unique_states) {
      print(paste("Processing state:", state))
      tips <- names(getStates(simmap_tree, type='tips'))[getStates(simmap_tree, type='tips') == state]
      
      if (length(tips) > 2) {
        subtree <- keep.tip(simmap_tree, tips)
        # Ensure the subset of data matches the tip labels in the subtree
        trait_data_subset <- input_data[match(subtree$tip.label, rownames(input_data)), , drop = FALSE]
        
        if (is.matrix(trait_data_subset) && nrow(trait_data_subset) > 2) {
          # Print number of rows to help debugging
          print(paste("Number of rows in filtered data:", nrow(trait_data_subset)))
          
          # Create a local environment to evaluate the formula
          eval_env <- new.env()
          assign("trait_data", trait_data_subset, envir = eval_env)
          
          # Parse and evaluate the formula in the local environment
          formula_obj <- as.formula(user_formula)
          environment(formula_obj) <- eval_env
          
          tryCatch({
            # Fit the model using the parsed formula
            model_fit <- mvgls(formula_obj, tree = subtree, model = model, ...)
            mvgls_models[[state]] <- model_fit
          }, error = function(e) {
            warning(paste("Failed to fit model for state", state, ":", e$message))
            mvgls_models[[state]] <- paste("Error:", e$message)
          })
        } else {
          warning(paste("Insufficient data for fitting model in state", state))
          mvgls_models[[state]] <- paste("Insufficient data: only", nrow(trait_data_subset), "rows available")
        }
      } else {
        warning(paste("No tips found for state", state))
        mvgls_models[[state]] <- "No tips found"
      }
    }
    
    return(mvgls_models)
  }
  
  #' Parallelized mvgls Fitting by Tip Regimes
  #'
  #' Similar to `fitMvglsOnSubtrees()` but utilizes parallel processing
  #' via `future.apply` to accelerate fitting across regimes.
  fitMvglsOnSubtrees.parallel <- function(user_formula, simmap_tree, input_data, model = 'BM', num_cores = parallel::detectCores(), ...) {
    # Set up parallel processing
    plan(multisession, workers = num_cores)
    
    # Get unique states and prepare a list with NULLs for each state
    unique_states <- unique(getStates(simmap_tree, type='tips'))
    mvgls_models <- setNames(vector("list", length(unique_states)), unique_states)
    
    # Prepare data for each state
    for (state in unique_states) {
      tips <- names(getStates(simmap_tree, type='tips'))[getStates(simmap_tree, type='tips') == state]
      
      if (length(tips) >= 2) {  # Only proceed if more >= 2 tips
        subtree <- keep.tip(simmap_tree, tips)
        trait_data_subset <- input_data[match(subtree$tip.label, rownames(input_data)), , drop = FALSE]
        
        # Check and convert data to appropriate type
        if (is.list(trait_data_subset) || !is.matrix(trait_data_subset)) {
          trait_data_subset <- as.matrix(trait_data_subset)
        }
        
        # Create a local environment to evaluate the formula
        eval_env <- new.env()
        assign("trait_data", trait_data_subset, envir = eval_env)
        formula_obj <- as.formula(user_formula)
        environment(formula_obj) <- eval_env
        
        # Try fitting the model, catch errors
        mvgls_models[[state]] <- tryCatch({
          mvgls(formula_obj, tree = subtree, model = model, ...)
        }, error = function(e) {
          warning(paste("Failed to fit model for state", state, ":", e$message))
          NULL  # Return NULL in case of an error
        })
      } else {
        warning(paste("Not enough tips for state", state, ": found", length(tips), "tips"))
        mvgls_models[[state]] <- NULL  # Explicitly set to NULL
      }
    }
    
    # Reset future's plan to sequential
    plan(sequential)
    
    return(mvgls_models)
  }
  
  #' Extract Penalized Precision Matrices from mvgls Model List
  extractPinvMatrices <- function(model_list) {
    if (!is.list(model_list)) {
      stop("The input must be a list of model objects.")
    }
    
    pinv_matrices <- list()  # Initialize an empty list to store Pinv matrices
    
    # Iterate over each model in the list
    for (state in names(model_list)) {
      model <- model_list[[state]]
      if (!is.null(model) && "sigma" %in% names(model)) {
        if ("Pinv" %in% names(model$sigma)) {
          pinv_matrices[[state]] <- model$sigma$Pinv
        } else {
          pinv_matrices[[state]] <- NULL  # Or you can choose not to add it at all if Pinv does not exist
        }
      } else {
        pinv_matrices[[state]] <- NULL  # Model is NULL or doesn't have a sigma component
      }
    }
    
    if (length(pinv_matrices) == 0) {
      warning("No valid Pinv matrices found in any models.")
    }
    
    return(pinv_matrices)
  }
  
  #' Compute Mean Trait Correlations from VCV Matrices
  calculateMeanCorrelation <- function(mat_list, absolute = FALSE) {
    means <- sapply(mat_list, function(mat) {
      cor_mat <- cov2cor(mat)
      if (absolute) cor_mat <- abs(cor_mat)
      # Use upper.tri() with diag = FALSE to select only off-diagonal elements in the upper triangle
      mean(cor_mat[upper.tri(cor_mat, diag = FALSE)], na.rm = TRUE)
    })
    return(means)
  }
  
  #' Compute Mean Trait Variances from VCV Matrices
  calculateMeanVariance <- function(vcv_list) {
    mean_variances <- sapply(vcv_list, function(vcv) {
      # Extract the diagonal elements representing the variances
      variances <- diag(vcv)
      # Calculate the mean of these variances
      mean(variances, na.rm = TRUE)
    })
    return(mean_variances)
  }
  
  #' Compute Median Trait Correlations from VCV Matrices
  calculateMedianCorrelation <- function(mat_list, absolute = FALSE) {
    medians <- sapply(mat_list, function(mat) {
      cor_mat <- cov2cor(mat)
      if (absolute) cor_mat <- abs(cor_mat)
      # Use upper.tri() with diag = FALSE to select only off-diagonal elements in the upper triangle
      median(cor_mat[upper.tri(cor_mat, diag = FALSE)], na.rm = TRUE)
    })
    return(medians)
  }
  
  #' Compute Median Trait Variances from VCV Matrices
  calculateMedianVariance <- function(vcv_list) {
    median_variances <- sapply(vcv_list, function(vcv) {
      # Extract the diagonal elements representing the variances
      variances <- diag(vcv)
      # Calculate the mean of these variances
      median(variances, na.rm = TRUE)
    })
    return(median_variances)
  }
  
  #' Match Correlation Estimates with Model Parameters
  matchCorrelationsWithParams <- function(correlation_means, params) {
    # Extract the names (labels) from both the correlation means and the parameters
    correlation_labels <- names(correlation_means)
    param_labels <- names(params)
    
    # Find common labels to ensure we only match existing ones
    common_labels <- intersect(correlation_labels, param_labels)
    
    # Create a data frame to store the matched results
    matched_data <- data.frame(
      State = common_labels,
      Correlation_Mean = correlation_means[common_labels],
      Param = params[common_labels],
      row.names = common_labels # Setting row names to match state labels
    )
    
    return(matched_data)
  }
  
  #' Match Variance Estimates with Model Parameters
  matchVarsWithParams <- function(var_means, params) {
    # Extract the names (labels) from both the correlation means and the parameters
    var_labels <- names(var_means)
    param_labels <- names(params)
    
    # Find common labels to ensure we only match existing ones
    common_labels <- intersect(var_labels, param_labels)
    
    # Create a data frame to store the matched results
    matched_data <- data.frame(
      State = common_labels,
      Var_Mean = var_means[common_labels],
      Param = params[common_labels],
      row.names = common_labels # Setting row names to match state labels
    )
    
    return(matched_data)
  }
  
  #' Parallel Post-hoc mvgls Fitting by Node-defined Clades
  #' Fits `mvgls` models to foreground and background clades defined by a list of nodes.
  fitMvglsByNodeRegions.parallel <- function(user_formula, tree, input_data, nodes, model = 'BM', num_cores = parallel::detectCores(), ...) {
    # Ensure required package is available
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required. Please install it with install.packages('future.apply').")
    }
    
    # Set up parallel processing
    future::plan(future::multisession, workers = num_cores)
    
    # Helper function to safely fit a pair of mvgls models for a node
    fit_pair_for_node <- function(node) {
      # Try to extract foreground clade; fail gracefully if invalid
      foreground_tree <- tryCatch({
        extract.clade(tree, node)
      }, error = function(e) {
        warning(paste("Could not extract clade for node", node, ":", e$message))
        return(NULL)
      })
      
      if (is.null(foreground_tree)) return(NULL)
      
      # Drop foreground tips to define background tree
      background_tree <- drop.tip(tree, foreground_tree$tip.label)
      
      # Match trait data
      foreground_data <- input_data[match(foreground_tree$tip.label, rownames(input_data)), , drop = FALSE]
      background_data <- input_data[match(background_tree$tip.label, rownames(input_data)), , drop = FALSE]
      
      # Safe model fitting (reused logic)
      fit_model_safe <- function(subtree, trait_data_subset) {
        if (is.list(trait_data_subset) || !is.matrix(trait_data_subset)) {
          trait_data_subset <- as.matrix(trait_data_subset)
        }
        
        if (nrow(trait_data_subset) <= 2) {
          return(NULL)
        }
        
        eval_env <- new.env()
        assign("trait_data", trait_data_subset, envir = eval_env)
        formula_obj <- as.formula(user_formula)
        environment(formula_obj) <- eval_env
        
        tryCatch({
          mvgls(formula_obj, tree = subtree, model = model, ...)
        }, error = function(e) {
          warning(paste("Failed to fit model at node", node, ":", e$message))
          return(NULL)
        })
      }
      
      # Fit both models
      fg_model <- fit_model_safe(foreground_tree, foreground_data)
      bg_model <- fit_model_safe(background_tree, background_data)
      
      # Return structured result
      return(list(
        node = node,
        foreground = fg_model,
        background = bg_model
      ))
    }
    
    # Run parallel model fitting
    results <- future.apply::future_lapply(nodes, fit_pair_for_node)
    
    # Label the list by node number
    names(results) <- paste0("node_", nodes)
    
    # Reset future plan
    future::plan(future::sequential)
    
    return(results)
  }
  
  #' Compute and Match Variances and Correlations to Regime Rates
  #'
  #' Uses `Pinv` matrices to derive average variances and correlations, and aligns them with evolutionary rates.
  generateVarsCors <- function(posthoc_models, param_table) {
    # Extract penalized covariance matrices
    mat <- extractPinvMatrices(posthoc_models)
    
    # Compute mean correlation (absolute) and mean variance for each regime
    corrs <- calculateMeanCorrelation(mat, absolute = TRUE)
    vars <- calculateMeanVariance(mat)
    
    # Match metrics to model parameters (e.g., evolutionary rates)
    matchcor <- matchCorrelationsWithParams(correlation_means = corrs, params = param_table)
    matchvars <- matchVarsWithParams(var_means = vars, params = param_table)
    
    # Match param values to their state names from param_table
    match_states <- sapply(matchvars$Param, function(p) {
      matched_names <- names(param_table)[which(param_table == p)]
      if (length(matched_names) == 1) {
        return(matched_names)
      } else {
        return(NA)  # ambiguous or unmatched
      }
    })
    
    # Combine into final dataframe
    vars_cors <- data.frame(
      rate = matchvars$Param,
      vars = matchvars$Var_Mean,
      corrs = matchcor$Correlation_Mean,
      State = match_states
    )
    
    rownames(vars_cors)<- vars_cors$State
    
    return(vars_cors)
  }
  
  #' Batch Post-hoc mvgls Fitting on Multiple Run Configurations
  #'
  #' Automates post-hoc fitting of `mvgls` models for each run in a run list.
  fitPosthocModels <- function(run_list, user_formula, input_data, num_cores = 3, model = 'BM', error = TRUE,
                               tree_element = "tree_no_uncertainty_untransformed") {
    for (run_name in names(run_list)) {
      cat("Fitting posthoc model for:", run_name, "\n")
      run <- run_list[[run_name]]
      
      # Fit the model
      posthoc_model <- fitMvglsOnSubtrees.parallel(
        user_formula = user_formula,
        num_cores = num_cores,
        simmap_tree = run[[tree_element]],
        input_data = input_data,
        model = model,
        error = error
      )
      
      # Store the result
      run_list[[run_name]]$posthoc <- posthoc_model
    }
    
    return(run_list)
  }
  
  #' Generate Variance-Correlation Summaries Across Model Runs
  #'
  #' Processes multiple runs and summarizes regime-specific trait variance, correlation, and rates.
  generateVarsCorsList <- function(
    run_list,
    posthoc_name = "posthoc",
    param_path = c("model_no_uncertainty", "param"),
    remove_high_corr = FALSE,
    corr_threshold = 0.99
  ) {
    results <- list()
    
    for (run_name in names(run_list)) {
      cat("Generating vars_cors for:", run_name, "\n")
      run <- run_list[[run_name]]
      
      posthoc_model <- run[[posthoc_name]]
      param_table <- run[[param_path[1]]][[param_path[2]]]
      
      vars_cors <- generateVarsCors(posthoc_model, param_table)
      
      if (remove_high_corr) {
        original_n <- nrow(vars_cors)
        vars_cors <- vars_cors[abs(vars_cors$corrs) <= corr_threshold, ]
        removed_n <- original_n - nrow(vars_cors)
        
        cat(sprintf("  Removed %d rows with |corrs| > %.2f\n", removed_n, corr_threshold))
      }
      
      results[[run_name]] <- vars_cors
    }
    
    return(results)
  }
  
  #' Plot Rate vs. Variance and Trait Integration (Correlation)
  #'
  #' Generates ggplot-based panels showing the relationship between evolutionary rate and
  #' post-hoc regime variance or trait correlation. Includes bootstrap CIs and residual filtering.
  plotVarsCors <- function(
    vars_cors_list,
    point_alpha = 0.8,
    point_size = 3,
    resid_sd_threshold_vars = 2,
    resid_sd_threshold_corrs = 2,
    n_boot = 1000,
    ci_level = 0.99
  ) {
    library(ggplot2)
    library(dplyr)
    library(patchwork)
    
    fisher_z_transform <- function(r) {
      0.5 * log((1 + r) / (1 - r))
    }
    
    # Combine and preprocess
    combined <- dplyr::bind_rows(
      lapply(names(vars_cors_list), function(name) {
        df <- vars_cors_list[[name]]
        df$model <- name
        df$fisher_z_corr <- fisher_z_transform(df$corrs)
        df$log_rate <- log(df$rate)
        df$log_vars <- log(df$vars)
        return(df)
      }),
      .id = "run_id"
    )
    
    # Residual filtering: Variance model
    lm_vars <- lm(log_rate ~ log_vars, data = combined)
    combined$vars_resid <- rstudent(lm_vars)
    combined_vars_clean <- combined %>% filter(abs(vars_resid) <= resid_sd_threshold_vars)
    removed_vars <- combined %>% filter(abs(vars_resid) > resid_sd_threshold_vars)
    if (nrow(removed_vars) > 0) {
      cat("\nRemoved from variance plot (studentized residuals > ", resid_sd_threshold_vars, "):\n")
      print(removed_vars[, c("model", "log_vars", "log_rate", "vars_resid")])
    }
    
    # Residual filtering: Correlation model
    lm_corrs <- lm(log_rate ~ fisher_z_corr, data = combined)
    combined$corrs_resid <- rstudent(lm_corrs)
    combined_corrs_clean <- combined %>% filter(abs(corrs_resid) <= resid_sd_threshold_corrs)
    removed_corrs <- combined %>% filter(abs(corrs_resid) > resid_sd_threshold_corrs)
    if (nrow(removed_corrs) > 0) {
      cat("\nRemoved from correlation plot (studentized residuals > ", resid_sd_threshold_corrs, "):\n")
      print(removed_corrs[, c("model", "fisher_z_corr", "log_rate", "corrs_resid")])
    }
    
    # Axis limits
    xlim_vars <- range(combined_vars_clean$log_vars, na.rm = TRUE)
    ylim_vars <- range(combined_vars_clean$log_rate, na.rm = TRUE)
    xlim_corrs <- range(combined_corrs_clean$fisher_z_corr, na.rm = TRUE)
    ylim_corrs <- range(combined_corrs_clean$log_rate, na.rm = TRUE)
    
    # Bootstrap CI: Variance
    x_seq_vars <- seq(xlim_vars[1], xlim_vars[2], length.out = 100)
    preds_vars <- replicate(n_boot, {
      boot_sample <- combined_vars_clean[sample(nrow(combined_vars_clean), replace = TRUE), ]
      predict(lm(log_rate ~ log_vars, data = boot_sample), newdata = data.frame(log_vars = x_seq_vars))
    })
    ci_vars <- apply(preds_vars, 1, quantile, probs = c((1 - ci_level)/2, 1 - (1 - ci_level)/2))
    mean_vars <- rowMeans(preds_vars)
    ci_df_vars <- data.frame(x = x_seq_vars, y = mean_vars, ymin = ci_vars[1, ], ymax = ci_vars[2, ])
    
    # Bootstrap CI: Correlation
    x_seq_corrs <- seq(xlim_corrs[1], xlim_corrs[2], length.out = 100)
    preds_corrs <- replicate(n_boot, {
      boot_sample <- combined_corrs_clean[sample(nrow(combined_corrs_clean), replace = TRUE), ]
      predict(lm(log_rate ~ fisher_z_corr, data = boot_sample), newdata = data.frame(fisher_z_corr = x_seq_corrs))
    })
    ci_corrs <- apply(preds_corrs, 1, quantile, probs = c((1 - ci_level)/2, 1 - (1 - ci_level)/2))
    mean_corrs <- rowMeans(preds_corrs)
    ci_df_corrs <- data.frame(x = x_seq_corrs, y = mean_corrs, ymin = ci_corrs[1, ], ymax = ci_corrs[2, ])
    
    # Regression captions
    lm_clean_vars <- lm(log_rate ~ log_vars, data = combined_vars_clean)
    sv <- summary(lm_clean_vars)
    caption_vars <- paste0(
      "y = ", round(coef(lm_clean_vars)[2], 3),
      "x + ", round(coef(lm_clean_vars)[1], 3),
      " | R² = ", round(sv$r.squared, 3),
      " | p = ", ifelse(sv$coefficients[2, 4] < 0.001, "< 0.001", format(round(sv$coefficients[2, 4], 3), nsmall = 3))
    )
    
    lm_clean_corrs <- lm(log_rate ~ fisher_z_corr, data = combined_corrs_clean)
    sc <- summary(lm_clean_corrs)
    caption_corrs <- paste0(
      "y = ", round(coef(lm_clean_corrs)[2], 3),
      "x + ", round(coef(lm_clean_corrs)[1], 3),
      " | R² = ", round(sc$r.squared, 3),
      " | p = ", ifelse(sc$coefficients[2, 4] < 0.001, "< 0.001", format(round(sc$coefficients[2, 4], 3), nsmall = 3))
    )
    
    # Theme
    custom_theme <- theme_minimal(base_size = 13) +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "right",
        plot.caption = element_text(hjust = 0.5)
      )
    
    # Plot 1: Variance
    p1 <- ggplot(combined_vars_clean, aes(x = log_vars, y = log_rate, color = model)) +
      geom_point(alpha = point_alpha, size = point_size) +
      geom_line(data = ci_df_vars, aes(x = x, y = y), inherit.aes = FALSE, color = "black", linewidth = 1) +
      geom_ribbon(data = ci_df_vars, aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE, alpha = 0.2, fill = "gray50") +
      coord_cartesian(xlim = xlim_vars, ylim = ylim_vars) +
      labs(
        x = "Log(Mean Regime Variance, \nPost-hoc Fit)",
        y = "Log(Mean Regime Variance, \nOriginal Model)",
        title = "Rate vs. Variance",
        caption = caption_vars
      ) +
      custom_theme
    
    # Plot 2: Correlation
    p2 <- ggplot(combined_corrs_clean, aes(x = fisher_z_corr, y = log_rate, color = model)) +
      geom_point(alpha = point_alpha, size = point_size) +
      geom_line(data = ci_df_corrs, aes(x = x, y = y), inherit.aes = FALSE, color = "black", linewidth = 1) +
      geom_ribbon(data = ci_df_corrs, aes(x = x, ymin = ymin, ymax = ymax), inherit.aes = FALSE, alpha = 0.2, fill = "gray50") +
      coord_cartesian(xlim = xlim_corrs, ylim = ylim_corrs) +
      labs(
        x = "Fisher Z(Trait Correlation, \nPost-hoc Fit)",
        y = "Log(Mean Regime Variance, \nOriginal Model)",
        title = "Rate vs. Trait Integration",
        caption = caption_corrs
      ) +
      custom_theme
    
    # Combine and return
    combined_plot <- p1 + p2 + plot_layout(guides = "collect")
    
    print(combined_plot)
    
    return(list(
      combined = combined_plot,
      rate_vs_variance = p1,
      rate_vs_correlation = p2
    ))
  }
  
  
}

#' Generate Scaled Viridis Color Palette for Rate Parameters
#'
#' This function creates a color mapping for rate parameters (e.g., evolutionary rates),
#' using the `viridis` color palette. Parameters are first sorted and normalized from 0 to 1,
#' then mapped to evenly spaced colors for clear visual differentiation.
generateViridisColorScale <- function(params) {
  # Sort parameters and keep their names
  sorted_indices <- order(params)
  sorted_params <- params[sorted_indices]
  
  # Normalize the sorted parameter values to a range from 0 to 1
  normalized_sorted_params <- (sorted_params - min(sorted_params)) / (max(sorted_params) - min(sorted_params))
  
  # Use the normalized values to get colors from the viridis palette
  colors <- viridis(length(normalized_sorted_params))
  
  # Associate each color with its state (name), using the sorted order
  named_sorted_colors <- setNames(colors, names(sorted_params))
  
  # Create a second list for the actual parameter values and their associated colors, also in sorted order
  param_color_mapping <- setNames(sorted_params, names(sorted_params))
  
  # Return both the original named sorted colors and the param-color mapping
  return(list("NamedColors" = named_sorted_colors, "ParamColorMapping" = param_color_mapping))
}

#' Generate High-Resolution Viridis Color Scale for Rate Parameters
#'
#' This function maps uniquely named rate parameters to colors in a high-resolution
#' viridis color gradient. Useful for fine-grained visualizations where subtle
#' differences in parameter values should be reflected in color.
generateViridisColorScale.fine <- function(params, num_colors = 100) {
  # Ensure params are unique and have names
  if(length(unique(params)) != length(params) || is.null(names(params))) {
    stop("Each parameter value must be unique and named.")
  }
  
  # Normalize the parameter values to the range of 1 to num_colors
  normalized_params <- (params - min(params)) / (max(params) - min(params))
  color_indices <- 1 + (num_colors - 1) * normalized_params
  
  # Round the indices and ensure they stay within bounds
  color_indices <- round(color_indices)
  color_indices[color_indices < 1] <- 1
  color_indices[color_indices > num_colors] <- num_colors
  
  # Generate a set of colors from the viridis palette
  all_colors <- viridis(num_colors)
  
  # Select colors based on the rounded indices
  selected_colors <- all_colors[color_indices]
  
  # Associate each color with its state (name)
  named_colors <- setNames(selected_colors, names(params))
  
  # Create a second list for the actual parameter values and their associated colors
  param_color_mapping <- setNames(params, names(params))
  
  # Return both the original named colors and the param-color mapping
  return(list("NamedColors" = named_colors, "ParamColorMapping" = param_color_mapping))
}


## function from Liam's blog for plotting
## http://blog.phytools.org/2024/03/even-more-about-plotting-discrete.html?m=1
plotFanTree.wTraits<-function(tree, X, colorvec, type=c("arc","fan"), color_palette=hcl.colors, trait_scale=0.07,
                              ...){
  X<-if(is.vector(X)) as.matrix(X) else X
  h<-max(nodeHeights(tree))
  d<-min(ncol(X)*trait_scale*h,h)
  type<-type[1]
  if(!(type%in%c("arc","fan"))) type<-"fan"
  ftype<-if(hasArg(ftype)) list(...)$ftype else "i"
  fsize<-if(hasArg(fsize)) list(...)$fsize else 0.5
  part<-if(hasArg(part)) list(...)$part else 
    min(0.99,(Ntip(tree)-2)/Ntip(tree))
  arc_height<-if(hasArg(arc_height)) list(...)$arc_height else 
    0.7
  spacer<-if(hasArg(spacer)) list(...)$spacer else 0.025
  spacer<-spacer*(2*pi*part/(Ntip(tree)-1))/2
  xlim<-if(hasArg(xlim)) list(...)$xlim else NULL
  ylim<-if(hasArg(ylim)) list(...)$ylim else NULL
  if(hasArg(colors)) colors<-list(...)$colors
  else {
    colors<-list()
    for(i in 1:ncol(X)){
      if(is.numeric(X[,i])){
        colors[[i]]<-setNames(color_palette(n=100), 1:100)
      } else {
        if(!is.factor(X[,i])) X[,i]<-as.factor(X[,i])
        colors[[i]]<-setNames(
          palette.colors(n=length(levels(X[,i]))),
          levels(X[,i]))
      }
    }
  }
  tt<-tree
  tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
    tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+d
  plotTree(tt,type=type,ftype=ftype,fsize=fsize,
           part=part,color="transparent",
           arc_height=arc_height*h/max(nodeHeights(tt)),
           xlim=xlim,ylim=ylim)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  outer_rad<-max(pp$xx)
  plotSimmap(tree,type=type,part=part,
             lwd=1,add=TRUE,xlim=pp$x.lim,ylim=pp$y.lim,
             arc_height=arc_height,ftype="off", colors=colorvec)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  inner_rad<-max(pp$xx)
  par(lend=3)
  for(i in 1:ncol(X)){
    if(is.numeric(X[,i])){
      x_seq<-seq(min(X[,i]),max(X[,i]),length.out=100)
      x_ind<-sapply(X[,i],function(x,y) which.min((x-y)^2),
                    y=x_seq)
      colors[[i]]<-colorRampPalette(colors[[i]])(n=100)
      cols<-colors[[i]][x_ind]
    } else {
      cols<-colors[[i]][X[tree$tip.label,i]]  
    }
    for(j in 1:Ntip(tree)){
      start<-if(pp$xx[j]>0) 
        (i-1)*(d/ncol(X))+(2/7)*(d/ncol(X)) else 
          -((i-1)*(d/ncol(X))+(2/7)*(d/ncol(X)))
      end<-if(pp$xx[j]>0) i*d/ncol(X) else -i*d/ncol(X)
      th<-atan(pp$yy[j]/pp$xx[j])
      theta<-(2*pi*part/(Ntip(tree)-1))/2-spacer
      sign<-if(pp$xx[j]>0) 1 else -1
      H1<-(sign*inner_rad+start)/cos(theta)
      H2<-(sign*inner_rad+end)/cos(theta)
      th_up<-th+theta
      th_down<-th-theta
      x<-c(H1*cos(th_down),H2*cos(th_down),
           H2*cos(th_up),H1*cos(th_up))
      y<-c(H1*sin(th_down),H2*sin(th_down),
           H2*sin(th_up),H1*sin(th_up))
      polygon(x,y,col=cols[j],border=FALSE)
    }
  }
  invisible(colors)
}

#Color assignment and visualization tools for rate data, inspired by BAMM-style rate maps.
{
  
  #' Assign discrete colors to rate values using simple binning
  #'
  #' Assigns colors to named rate values by dividing them into bins based on one of three methods:
  #' "linear", "quantile", or "jenks". Supports log transformation, clamping to a value range, and
  #' simple color palette selection.
  assignRateColors.simple <- function(rates, breaksmethod = "linear", logcolor = FALSE, 
                                      color.interval = NULL, palette = "temperature", 
                                      nbreaks = NULL) {
    if (is.null(nbreaks)) {
      nbreaks <- max(2, min(10, round(sqrt(length(rates)))))  # Example heuristic
      cat(paste("Input data have been automatically binned into ", nbreaks, "categories."))
    }
    if (!is.numeric(rates) || is.null(names(rates))) {
      stop("`rates` must be a named numeric vector.")
    }
    if (!breaksmethod %in% c("linear", "quantile", "jenks")) {
      stop("`breaksmethod` must be one of 'linear', 'quantile', or 'jenks'.")
    }
    if (logcolor && any(rates <= 0)) {
      stop("All rates must be positive for log transformation.")
    }
    if (!is.null(color.interval) && (length(color.interval) != 2 || !is.numeric(color.interval))) {
      stop("`color.interval` must be a numeric vector of length 2.")
    }
    
    if (logcolor) {
      rates <- log(rates)
    }
    if (!is.null(color.interval)) {
      rates <- pmin(pmax(rates, color.interval[1]), color.interval[2])
    }
    
    breaks <- NULL
    if (breaksmethod == "linear") {
      breaks <- seq(min(rates), max(rates), length.out = nbreaks + 1)
    } else if (breaksmethod == "quantile") {
      breaks <- quantile(rates, probs = seq(0, 1, length.out = nbreaks + 1))
    } else if (breaksmethod == "jenks") {
      if (!requireNamespace("classInt", quietly = TRUE)) {
        stop("The `classInt` package is required for Jenks natural breaks.")
      }
      if (length(unique(rates)) < nbreaks) {
        nbreaks <- length(unique(rates))
        warning(sprintf("Adjusted nbreaks to %d due to limited unique values.", nbreaks))
      }
      breaks <- classInt::classIntervals(rates, n = nbreaks, style = "jenks")$brks
    }
    
    getPalette <- function(pal, nbreaks) {
      if (length(pal) >= 3) {
        return(colorRampPalette(pal, space = 'Lab')(nbreaks))
      } else if (tolower(pal) == "viridis") {
        if (!requireNamespace("viridis", quietly = TRUE)) {
          stop("The `viridis` package is required for the viridis palette.")
        }
        return(viridis::viridis(nbreaks))
      } else if (tolower(pal) == "temperature") {
        if (!requireNamespace("gplots", quietly = TRUE)) {
          stop("The `gplots` package is required for the temperature palette.")
        }
        return(gplots::rich.colors(nbreaks))
      } else if (tolower(pal) == "terrain") {
        return(terrain.colors(nbreaks))
      } else {
        stop("Unsupported palette. Use a valid custom palette or 'viridis', 'temperature', 'terrain'.")
      }
    }
    
    color_generator <- getPalette(palette, nbreaks)
    color_bins <- cut(rates, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    assigned_colors <- color_generator[color_bins]
    
    names(assigned_colors) <- names(rates)
    return(list(colors = assigned_colors, breaks = breaks))
  }
  #' Flexible rate-to-color mapping with advanced binning and palette options
  #'
  #' Provides robust, highly configurable binning and color assignment for rate data.
  #' Supports over a dozen break methods (e.g., "dpih", "box", "jenks"), customizable palettes,
  #' jitter handling for unstable breaks, and reversed scales.
  assignRateColors <- function(rates, breaksmethod = "linear", logcolor = FALSE, 
                               color.interval = NULL, palette = "temperature", 
                               nbreaks = NULL, jitter_on_failure = TRUE, 
                               reverse = FALSE, ...) {
    if (is.null(nbreaks)) {
      nbreaks <- max(2, min(10, round(sqrt(length(rates)))))  # Example heuristic
      cat(paste("Input data have been automatically binned into ", nbreaks, "categories.\n"))
    }
    if (!is.numeric(rates) || is.null(names(rates))) {
      stop("`rates` must be a named numeric vector.")
    }
    # Update breaksmethod validation to include all supported styles
    valid_methods <- c("linear", "quantile", "jenks", "sd", "equal", 
                       "pretty", "kmeans", "hclust", "bclust", "fisher", 
                       "dpih", "headtails", "maximum", "box")
    if (!breaksmethod %in% valid_methods) {
      stop(paste("`breaksmethod` must be one of:", paste(valid_methods, collapse = ", ")))
    }
    if (logcolor && any(rates <= 0)) {
      stop("All rates must be positive for log transformation.")
    }
    if (!is.null(color.interval) && (length(color.interval) != 2 || !is.numeric(color.interval))) {
      stop("`color.interval` must be a numeric vector of length 2.")
    }
    
    # Log transformation
    if (logcolor) {
      rates <- log(rates)
    }
    # Clamp rates to the specified color interval
    if (!is.null(color.interval)) {
      rates <- pmin(pmax(rates, color.interval[1]), color.interval[2])
    }
    
    # Predefined palettes
    predefined_palettes <- list(
      "RdYlBu" = RColorBrewer::brewer.pal(11, "RdYlBu"),
      "BrBG" = RColorBrewer::brewer.pal(11, "BrBG"),
      "PiYG" = RColorBrewer::brewer.pal(11, "PiYG"),
      "PRGn" = RColorBrewer::brewer.pal(11, "PRGn"),
      "PuOr" = RColorBrewer::brewer.pal(11, "PuOr"),
      "RdBu" = RColorBrewer::brewer.pal(11, "RdBu"),
      "BuOr" = c("blue", "white", "orange"),
      "BuOrRd" = c("blue", "orange", "red"),
      "DkRdBu" = c("darkred", "white", "darkblue"),
      "BuDkOr" = c("blue", "darkorange"),
      "GnPu" = c("green", "purple"),
      "RdYlGn" = RColorBrewer::brewer.pal(11, "RdYlGn"),
      "Spectral" = RColorBrewer::brewer.pal(11, "Spectral"),
      "temperature" = gplots::rich.colors(64),
      "terrain" = terrain.colors(nbreaks),
      "grayscale" = gray.colors(nbreaks),
      "revgray" = rev(gray.colors(nbreaks))
    )
    
    # Compute breaks based on the selected breaksmethod
    breaks <- NULL
    if (breaksmethod == "linear") {
      breaks <- seq(min(rates), max(rates), length.out = nbreaks + 1)
    } else if (breaksmethod == "quantile") {
      breaks <- quantile(rates, probs = seq(0, 1, length.out = nbreaks + 1))
    } else if (breaksmethod == "dpih") {
      tryCatch({
        # Attempt to compute breaks using dpih
        breaks <- classInt::classIntervals(rates, n = nbreaks, style = "dpih", ...)$brks
        
        # Ensure unique breaks
        breaks <- unique(breaks)
        
        # Check for sufficient number of unique breaks
        if (length(breaks) < 2) {
          stop("Insufficient unique breaks for dpih method. Cannot proceed.")
        }
        
        # Ensure breaks span the full range of rates
        if (min(breaks) > min(rates) || max(breaks) < max(rates)) {
          warning("Extending breaks to include full range of data for dpih method.")
          breaks <- c(min(rates), breaks, max(rates))
          breaks <- unique(breaks)
        }
      }, error = function(e) {
        # If dpih fails, stop execution with an informative error
        stop("Error with dpih method: Unable to compute valid breaks.")
      })
    } else if (breaksmethod == "box") {
      tryCatch({
        breaks <- classInt::classIntervals(rates, n = nbreaks, style = "box", ...)$brks
        breaks[is.infinite(breaks)] <- min(rates) - 0.05 * abs(min(rates))
      }, error = function(e) {
        warning("Box method failed. Falling back to linear breaks.")
        breaks <- seq(min(rates), max(rates), length.out = nbreaks + 1)
      })
    } else if (breaksmethod == "jenks") {
      # Handle jenks method with progressive jitter
      success <- FALSE
      jitter_amount <- 0  # Default jitter amount is 0
      for (jitter_attempt in 1:100) {
        if (jitter_attempt == 1) {
          # First attempt without jitter
          try_breaks <- tryCatch(
            classInt::classIntervals(rates, n = nbreaks, style = "jenks", ...)$brks,
            error = function(e) NULL
          )
        } else {
          # Add progressive jitter
          jitter_amount <- jitter_attempt * 0.01
          jittered_rates <- rates + runif(length(rates), -jitter_amount, jitter_amount)
          try_breaks <- tryCatch(
            classInt::classIntervals(jittered_rates, n = nbreaks, style = "jenks", ...)$brks,
            error = function(e) NULL
          )
        }
        
        # Check if breaks are unique
        if (!is.null(try_breaks) && length(unique(try_breaks)) == length(try_breaks)) {
          breaks <- try_breaks
          
          # Ensure breaks span the full range of rates
          if (min(breaks) > min(rates)) {
            breaks <- c(min(rates), breaks)
          }
          if (max(breaks) < max(rates)) {
            breaks <- c(breaks, max(rates))
          }
          breaks <- unique(breaks)  # Ensure uniqueness after modification
          
          cat(sprintf("Jitter applied for 'jenks': %e after %d attempts\n", jitter_amount, jitter_attempt))
          success <- TRUE
          break
        }
      }
      
      # If all attempts fail, stop with an error
      if (!success) {
        stop("Failed to find unique breaks with 'jenks' after 100 attempts.")
      }
    } else {
      if (!requireNamespace("classInt", quietly = TRUE)) {
        stop("The `classInt` package is required for advanced breaks methods.")
      }
      if (length(unique(rates)) < nbreaks) {
        nbreaks <- length(unique(rates))
        warning(sprintf("Adjusted nbreaks to %d due to limited unique values.", nbreaks))
      }
      breaks <- classInt::classIntervals(rates, n = nbreaks, style = breaksmethod, ...)$brks
    }
    
    # Generate colors
    getPalette <- function(pal, NCOLORS) {
      if (pal %in% names(predefined_palettes)) {
        color_list <- predefined_palettes[[pal]]
        if (reverse) {
          color_list <- rev(color_list)  # Reverse the palette if reverse = TRUE
        }
        return(colorRampPalette(color_list)(NCOLORS))
      } else if (tolower(pal) == "viridis") {
        if (!requireNamespace("viridis", quietly = TRUE)) {
          stop("The `viridis` package is required for the viridis palette.")
        }
        colors <- viridis::viridis(NCOLORS)
        if (reverse) {
          colors <- rev(colors)  # Reverse if requested
        }
        return(colors)
      } else {
        stop("Unsupported palette. Use a predefined palette or 'viridis'.")
      }
    }
    
    color_generator <- getPalette(palette, length(breaks) - 1)
    color_bins <- cut(rates, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    assigned_colors <- color_generator[color_bins]
    
    # Store numeric breaks and their corresponding color
    break_colors <- list(
      breaks = breaks,  # The breakpoints used for binning
      colors = color_generator  # The colors assigned to each break
    )
    
    names(assigned_colors) <- names(rates)
    return(list(colors = assigned_colors, breaks = breaks, break_colors = break_colors))
  }
  
  #' Draw a vertical color bar legend for a rate palette
  #'
  #' Uses ggplot2 to render a legend-like bar showing the relationship between rate bins and
  #' assigned colors. Typically used to accompany rate-mapped visualizations.
  generate_color_bar <- function(break_colors, bar_title = "Temperature Seasonality", n_ticks = 5) {
    # Validate input
    if (!is.list(break_colors) || !all(c("breaks", "colors") %in% names(break_colors))) {
      stop("Input must be a list with 'breaks' and 'colors'.")
    }
    if (length(break_colors$breaks) != length(break_colors$colors) + 1) {
      stop("Number of colors must be one less than the number of breaks.")
    }
    
    # Create bin data frame
    color_bar_df <- data.frame(
      ymin = head(break_colors$breaks, -1),
      ymax = tail(break_colors$breaks, -1),
      fill = break_colors$colors
    )
    
    # Build color bar plot
    color_bar <- ggplot(color_bar_df) +
      geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = ymin, ymax = ymax, fill = fill)) +
      scale_fill_identity() +
      scale_y_continuous(breaks = pretty(range(break_colors$breaks), n = n_ticks)) +
      labs(title = bar_title, fill = NULL) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_text(angle = 90),
        plot.title = element_text(hjust = 1, size = 10, face = "bold", margin = margin(b = 5)),
        legend.position = "none"
      ) +
      coord_cartesian(clip = "off")
    
    return(color_bar)
  }
  
  #' Embed a ggplot figure as an inset in a base R plot
  #'
  #' Inserts a ggplot2 plot object into an existing base plot using specified coordinates,
  #' corner keywords (e.g., "topright"), or a bounding box. Useful for embedding legends or
  #' summary charts.
  subplot.gg <- function(plot_object, x = "topright", y = NULL, size = NULL, bbox = NULL, 
                         units = "npc", just = c(0.5, 0.5)) {
    # plot_object: a ggplot2 object that you wish to inset
    # x, y: if bbox is NULL, x can be a numeric coordinate (with y) in the base plot's user space,
    #       or a character string like "topright", "topleft", "bottomright", or "bottomleft".
    # size: a vector c(width, height) specifying the size of the inset (in units defined by 'units').
    # bbox: alternatively, a vector c(xmin, ymin, xmax, ymax) in user coordinates that defines the inset box.
    # units: units for the size when using 'size' (default "npc"; could be "inches" if you prefer).
    # just: justification of the inset relative to the specified center (default is centered).
    
    # If a bounding box is provided, calculate center and size from it:
    if (!is.null(bbox)) {
      if (length(bbox) != 4)
        stop("bbox must be a vector of length 4: c(xmin, ymin, xmax, ymax)")
      xleft   <- grconvertX(bbox[1], from = "user", to = "npc")
      ybottom <- grconvertY(bbox[2], from = "user", to = "npc")
      xright  <- grconvertX(bbox[3], from = "user", to = "npc")
      ytop    <- grconvertY(bbox[4], from = "user", to = "npc")
      width  <- xright - xleft
      height <- ytop - ybottom
      x_center <- xleft + width/2
      y_center <- ybottom + height/2
    } else {
      # Determine center from x, y.
      if (is.character(x)) {
        usr <- par("usr")  # user coordinates of the current plot
        pos <- tolower(x)
        if (pos == "topright") {
          x_center <- usr[2]
          y_center <- usr[4]
          just <- c(1, 1)
        } else if (pos == "topleft") {
          x_center <- usr[1]
          y_center <- usr[4]
          just <- c(0, 1)
        } else if (pos == "bottomright") {
          x_center <- usr[2]
          y_center <- usr[3]
          just <- c(1, 0)
        } else if (pos == "bottomleft") {
          x_center <- usr[1]
          y_center <- usr[3]
          just <- c(0, 0)
        } else {
          x_center <- (usr[1] + usr[2]) / 2
          y_center <- (usr[3] + usr[4]) / 2
          just <- c(0.5, 0.5)
        }
      } else {
        if (is.null(y))
          stop("If x is numeric, y must also be provided")
        x_center <- x
        y_center <- y
      }
      # Convert center from user coordinates to npc:
      x_center <- grconvertX(x_center, from = "user", to = "npc")
      y_center <- grconvertY(y_center, from = "user", to = "npc")
      # If size is not provided, use a default (e.g., 0.15 by 0.4 in npc units)
      if (is.null(size)) {
        width <- 0.15
        height <- 0.4
      } else {
        width  <- size[1]
        height <- size[2]
      }
    }
    
    # Create a viewport at the calculated center with the given size.
    vp <- viewport(x = x_center, y = y_center, width = unit(width, units), height = unit(height, units), 
                   just = just)
    # Print the ggplot object inside this viewport.
    print(plot_object, vp = vp)
  }
  
  #' Base R histogram with rate-based color mapping
  #'
  #' Generates a histogram of rates using base graphics, coloring each bin according to the
  #' associated rate values. Includes options for density display, axis label styling, and break overlays.
  ratesHistogram <- function(rates, colors, breaks, plotBrks = TRUE, 
                             xlab = NULL, ylab = NULL, title = NULL, 
                             useDensity = FALSE, logscale = FALSE, 
                             lwd = 0.2, lty = 1, brksCol = "black", 
                             xBuffer = 0.05, yBuffer = 0.05, 
                             cex.axis = 0.75, cex.lab = 1, cex.main = 1.1, 
                             x.axis.label.line = 2, y.axis.label.line = 2, 
                             x.tick.label.line = 0.5, y.tick.label.line = 0.5, 
                             plotOverlay = NULL, ...) {
    if (!is.numeric(rates) || is.null(names(rates))) {
      stop("`rates` must be a named numeric vector.")
    }
    if (length(colors) != length(rates) || !all(names(colors) == names(rates))) {
      stop("`colors` must be a named vector of the same length as `rates`, with matching names.")
    }
    if (is.null(breaks) || length(breaks) < 2) {
      stop("`breaks` must be a numeric vector of at least two breakpoints.")
    }
    
    if (logscale) {
      if (any(rates <= 0)) {
        stop("Log transformation is not possible for non-positive rates.")
      }
      rates <- log(rates)
    }
    
    hist_data <- hist(rates, breaks = breaks, plot = FALSE)
    midpoints <- hist_data$mids
    heights <- if (useDensity) hist_data$density else hist_data$counts
    
    if (is.null(ylab)) {
      ylab <- if (useDensity) "Density" else "Counts"
    }
    if (is.null(xlab)) {
      xlab <- if (logscale) "Rates - log scale" else "Rates"
    }
    
    bin_indices <- cut(rates, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    per_bin_color <- sapply(seq_along(midpoints), function(b) {
      inds <- which(bin_indices == b)
      if (length(inds) > 0) {
        return(colors[inds[1]])
      } else {
        return("gray")
      }
    })
    
    x_range <- range(breaks)
    y_max <- max(heights, 0)
    x_min_with_buffer <- x_range[1] - (x_range[2] - x_range[1]) * xBuffer
    y_min_with_buffer <- -y_max * yBuffer
    
    plot.new()
    plot.window(xlim = c(x_min_with_buffer, x_range[2]), 
                ylim = c(y_min_with_buffer, y_max), 
                xaxs = "i", yaxs = "i")
    box(bty = 'L')
    
    for (i in seq_along(midpoints)) {
      rect(xleft = breaks[i], xright = breaks[i + 1],
           ybottom = 0, ytop = heights[i],
           col = per_bin_color[i], border = NA)
    }
    
    # Draw x-axis ticks and labels
    par(mgp = c(3, x.tick.label.line, 0))  # Set tick label distance for x-axis
    x_ticks <- signif(seq(x_min_with_buffer, x_range[2], length.out = 5), 2)
    x_ticks <- x_ticks[-1]
    axis(1, at = x_ticks, cex.axis = cex.axis, lwd = 1)
    
    # Draw y-axis ticks and labels
    par(mgp = c(3, y.tick.label.line, 0))  # Set tick label distance for y-axis
    y_ticks <- seq(0, y_max, length.out = 5)
    if (useDensity) {
      y_ticks <- signif(y_ticks, 2)
    }
    axis(2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = cex.axis, lwd = 1)
    
    # Add axis labels
    mtext(xlab, side = 1, line = x.axis.label.line, cex = cex.lab, ...)
    mtext(ylab, side = 2, line = y.axis.label.line, cex = cex.lab, ...)
    
    # Add title
    if (!is.null(title)) {
      title(main = title, line = 0.5, cex.main = cex.main)
    }
    
    # Add vertical lines for breaks
    if (plotBrks) {
      abline(v = breaks, lwd = lwd, lty = lty, col = brksCol)
    }
    
    # Call the overlay function if provided
    if (!is.null(plotOverlay) && is.function(plotOverlay)) {
      plotOverlay()
    }
    
  }
  
  #' Simulate multi-modal rate distributions
  #'
  #' Generates synthetic rate values from a mixture of distributions
  #' (normal, exponential, uniform). Useful for testing color binning and visualization.
  simulateRates <- function(num_modes = 4, total_samples = 2000,
                            distribution_types = c("normal", "exponential", "uniform"),
                            range_normal = c(0.5, 0.1), range_exponential = c(1, 2),
                            range_uniform = c(0, 1), seed = NULL) {
    if (!is.null(seed)) set.seed(seed)    # Set seed if provided
    
    # Distribute total samples equally among modes
    samples_per_mode <- round(total_samples / num_modes)
    
    rates <- c()                          # Initialize vector to store rates
    
    # Generate data for each mode
    for (i in 1:num_modes) {
      dist_type <- sample(distribution_types, 1)  # Randomly select a distribution type
      
      if (dist_type == "normal") {
        mean <- runif(1, range_normal[1] - range_normal[2], range_normal[1] + range_normal[2])
        sd <- runif(1, 0.05, range_normal[2])
        rates <- c(rates, rnorm(samples_per_mode, mean = mean, sd = sd))
        
      } else if (dist_type == "exponential") {
        rate <- runif(1, range_exponential[1], range_exponential[2])
        rates <- c(rates, rexp(samples_per_mode, rate = rate))
        
      } else if (dist_type == "uniform") {
        min_val <- runif(1, range_uniform[1], range_uniform[2] - 0.1)
        max_val <- runif(1, min_val + 0.1, range_uniform[2])
        rates <- c(rates, runif(samples_per_mode, min = min_val, max = max_val))
      }
    }
    
    # Add names to rates
    names(rates) <- paste0("Rate_", seq_along(rates))
    
    # Shuffle the rates to avoid clustering by mode
    rates <- sample(rates)
    
    return(rates)
  }
  
  #' ggplot2-based histogram of rates with bin-specific colors
  #'
  #' Similar to ratesHistogram(), but implemented using ggplot2 for more control and styling.
  #' Colors bins based on rate values and optionally overlays break lines.
  ratesHistogramGG <- function(rates, colors, breaks,
                               plotBrks = TRUE, 
                               xlab = NULL, ylab = NULL, title = NULL, 
                               useDensity = FALSE, logscale = FALSE, 
                               lwd = 0.2, lty = 1, brksCol = "black", 
                               xBuffer = 0.05, yBuffer = 0.05, 
                               cex.axis = 0.75, cex.lab = 1, cex.main = 1.1,
                               x.axis.label.line = 2, y.axis.label.line = 2,
                               x.tick.label.line = 0.5, y.tick.label.line = 0.5,
                               ...) {
    # 1. Basic checks
    if (!is.numeric(rates) || is.null(names(rates))) {
      stop("`rates` must be a named numeric vector.")
    }
    if (length(colors) != length(rates) || !all(names(colors) == names(rates))) {
      stop("`colors` must be a named vector of the same length as `rates`, with matching names.")
    }
    if (is.null(breaks) || length(breaks) < 2) {
      stop("`breaks` must be a numeric vector of at least two breakpoints.")
    }
    
    # 2. Log-scale transformation
    if (logscale) {
      if (any(rates <= 0)) {
        stop("Log transformation is not possible for non-positive rates.")
      }
      rates <- log(rates)
    }
    
    # 3. Compute histogram info
    hist_data <- hist(rates, breaks = breaks, plot = FALSE)
    bin_counts  <- hist_data$counts
    bin_density <- hist_data$density
    heights     <- if (useDensity) bin_density else bin_counts
    midpoints   <- hist_data$mids
    
    # 4. Axis labels
    if (is.null(ylab)) {
      ylab <- if (useDensity) "Density" else "Counts"
    }
    if (is.null(xlab)) {
      xlab <- if (logscale) "Rates - log scale" else "Rates"
    }
    
    # 5. Per-bin color assignment: first observation in each bin
    bin_indices <- cut(rates, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    per_bin_color <- sapply(seq_along(midpoints), function(i) {
      inds <- which(bin_indices == i)
      if (length(inds) > 0) {
        first_name <- names(rates)[inds[1]]
        colors[first_name]
      } else {
        "gray"
      }
    })
    
    # 6. Create data frame for plotting each bin as a rectangle
    bin_df <- data.frame(
      xleft    = head(breaks, -1),
      xright   = tail(breaks, -1),
      ybottom  = 0,
      ytop     = heights,
      fillcolor = per_bin_color
    )
    
    # 7. Determine axis limits (with buffer)
    x_range <- range(breaks)
    y_max   <- max(heights, 0)
    x_min_with_buffer <- x_range[1] - (x_range[2] - x_range[1]) * xBuffer
    x_max_with_buffer <- x_range[2]
    y_min_with_buffer <- 0 - y_max * yBuffer
    
    # 8. Manual tick positions
    x_ticks <- signif(seq(x_min_with_buffer, x_max_with_buffer, length.out = 5), 2)
    if (length(x_ticks) > 1) {
      x_ticks <- x_ticks[-1]
    }
    y_ticks <- seq(0, y_max, length.out = 5)
    if (useDensity) {
      y_ticks <- signif(y_ticks, 2)
    }
    
    # 9. Build the ggplot
    p <- ggplot(bin_df) +
      geom_rect(
        aes(xmin = xleft, xmax = xright, ymin = ybottom, ymax = ytop, fill = fillcolor),
        color = NA
      ) +
      scale_fill_identity() +
      coord_cartesian(
        xlim = c(x_min_with_buffer, x_max_with_buffer),
        ylim = c(y_min_with_buffer, y_max),
        expand = FALSE
      ) +
      scale_x_continuous(
        breaks = x_ticks,
        labels = x_ticks
      ) +
      scale_y_continuous(
        breaks = y_ticks,
        labels = y_ticks
      ) +
      labs(
        x = xlab,
        y = ylab,
        title = title
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid          = element_blank(),
        panel.border        = element_blank(),
        axis.line.x.bottom  = element_line(color = "black", linewidth = 0.5),
        axis.line.y.left    = element_line(color = "black", linewidth = 0.5),
        axis.line.x.top     = element_blank(),
        axis.line.y.right   = element_blank(),
        axis.ticks          = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length   = unit(0.2, "lines"),
        axis.text.x         = element_text(
          size = rel(cex.axis),
          margin = margin(t = 5.5 * x.tick.label.line, unit = "pt")
        ),
        axis.text.y         = element_text(
          size = rel(cex.axis),
          margin = margin(r = 5.5 * y.tick.label.line, unit = "pt"),
          angle = 0
        ),
        axis.title.x        = element_text(
          size = rel(cex.lab),
          margin = margin(t = 5.5 * x.axis.label.line, unit = "pt")
        ),
        axis.title.y        = element_text(
          size = rel(cex.lab),
          margin = margin(r = 5.5 * y.axis.label.line, unit = "pt"),
          angle = 90
        ),
        plot.title          = element_text(
          hjust = 0.5,
          size = rel(cex.main),
          margin = margin(b = 5.5 * 0.5, unit = "pt")
        ),
        # Default R plot margins: c(5.1, 4.1, 4.1, 2.1)
        plot.margin         = margin(5.1 * 5.5, 4.1 * 5.5, 4.1 * 5.5, 2.1 * 5.5, unit = "pt")
      )
    
    if (plotBrks) {
      p <- p + geom_vline(
        xintercept = breaks,
        color = brksCol,
        linetype = lty,
        linewidth = lwd
      )
    }
    
    return(p)
  }
  
}


#' Scatterplot with embedded color legend for rate visualization
#'
#' Creates a customizable ggplot2 scatterplot where each point is colored using a precomputed color mapping
#' (e.g., based on evolutionary rates). An optional vertical reference line at x = 0 can be included,
#' and a compact color bar legend is embedded into the top-right of the plot.
generate_scatter_with_colorbar <- function(
    data,
    x_var,
    y_var,
    color_var,  # Name of the column containing precomputed colors
    color_palette,
    title = "Effect of Temperature Variability on Evolution Rates",
    x_label = "Local species richness (log weighted counts)",
    y_label = "Range weighted lineage rate",
    point_size = 2,   # Default point size
    point_alpha = 0.75, # Default transparency
    bar_title = NULL,
    xmin_perc = 0.8,
    ymin_perc = 0.85,
    y_nudge = 0,
    x_limits = NULL,
    zero_intercept = FALSE,  # NEW: Toggle for vertical line at x = 0
    colorbar = TRUE
) {
  # 1. Ensure the color column exists in the dataset
  if (!color_var %in% colnames(data)) {
    stop(paste("Column", color_var, "not found in dataset."))
  }
  
  # 2. Generate the scatterplot
  scatter_plot <- ggplot(
    data, 
    aes(x = .data[[x_var]], 
        y = .data[[y_var]], 
        color = .data[[color_var]])  
  )
  
  # **NEW: Conditionally add the vertical line (BEHIND points)**
  if (zero_intercept) {
    scatter_plot <- scatter_plot + geom_vline(xintercept = 0, linetype = 2, color = 'grey')
  }
  
  # Add scatter points
  scatter_plot <- scatter_plot + 
    geom_point(alpha = point_alpha, size = point_size) +  # Uses user-defined size & transparency
    scale_color_identity() +  
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = 14) +
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
  
  # Apply X-axis limits if provided by user
  if (!is.null(x_limits)) {
    scatter_plot <- scatter_plot + xlim(x_limits)
  }
  
  # 4. Generate the color bar
  color_bar <- generate_color_bar(color_palette, bar_title = bar_title)
  
  # 5. Convert the color bar plot into a grob
  color_bar_grob <- ggplotGrob(color_bar)
  
  # 6. Overlay the color bar in the top-right corner
  final_xmin <- if (!is.null(x_limits)) x_limits[2] * xmin_perc else max(data[[x_var]]) * xmin_perc
  final_xmax <- if (!is.null(x_limits)) x_limits[2] else max(data[[x_var]])
  
  if(colorbar==T){
    final_plot <- scatter_plot +
      annotation_custom(
        color_bar_grob, 
        xmin = final_xmin, 
        xmax = final_xmax, 
        ymin = min(data[[y_var]]) * ymin_perc,
        ymax = max(data[[y_var]]) + y_nudge
      )
  } else {
    final_plot <- scatter_plot
  }
  
  
  return(final_plot)
}


#' Generate multi-line legend summary for IC model comparison
#'
#' Constructs a text-based summary comparing baseline and optimal information criterion (IC) scores
#' from a model configuration search object. Designed to be optionally drawn as a multi-line
#' legend on a plot, summarizing baseline score, optimal score, and the delta (ΔIC).
generateLegendLabel <- function(object, x = "bottomleft", plot = TRUE, inset_y = 0.02) {
  # Automatically get the name of the object
  object_name <- deparse(substitute(object))
  
  # Validate the structure and types of the input object
  if (!is.list(object) || 
      !"IC_used" %in% names(object) || 
      !"optimal_ic" %in% names(object) || 
      !"baseline_ic" %in% names(object) || 
      !is.numeric(object$optimal_ic[[1]]) || 
      !is.numeric(object$baseline_ic[[1]])) {
    stop("Input object does not have the correct structure or types")
  }
  
  # Retrieve the type of IC used (BIC or GIC)
  ic_type <- object$IC_used
  
  # Calculate the absolute difference between baseline and optimal IC
  delta_ic <- abs(object$baseline_ic[[1]] - object$optimal_ic[[1]])
  
  # Create the legend text as individual lines
  legend_text <- c(
    sprintf("%s", object_name),
    sprintf("Baseline %s is %.2f", ic_type, object$baseline_ic[[1]]),
    sprintf("Optimal %s is %.2f", ic_type, object$optimal_ic[[1]]),
    sprintf("Δ%s is %.2f", ic_type, delta_ic)  # Delta symbol
  )
  
  # Plot the legend only if plot is TRUE
  if (plot) {
    # Using `legend()` to create multi-line legend
    legend(x, legend = legend_text, bty = "n", cex = 1.5, inset = c(0, inset_y))
  }
  
  return(legend_text)
}


#' Count and summarize state transitions on a node-based SIMMAP tree
#'
#' Analyzes a node-mapped phylogenetic tree in SIMMAP format to compute:
#' - Number of state transitions at internal nodes
#' - Time spent in each state across branches
#' - Transition matrix summarizing directional changes
countNodeTransitionsSimmap <- function(tree) {
  if (!inherits(tree, "phylo") || is.null(tree$maps) || is.null(tree$edge.length)) {
    stop("Tree must be a 'phylo' object and must contain 'maps' and 'edge.length' attributes")
  }
  
  # Extract all unique states from the maps
  all_states <- unique(unlist(lapply(tree$maps, names)))
  
  # Initialize the transition matrix and time in states
  num_states <- length(all_states)
  transition_matrix <- matrix(0, nrow = num_states, ncol = num_states,
                              dimnames = list(all_states, all_states))
  time_in_states <- setNames(vector("numeric", length(all_states)), all_states)
  
  # Navigate through each branch to accumulate time spent in each state
  for (i in seq_along(tree$maps)) {
    map <- tree$maps[[i]]
    branch_length <- tree$edge.length[i]
    if (length(map) > 1) {
      # Calculate segment lengths assuming equal division among segments
      segment_lengths <- rep(branch_length / length(map), length(map))
    } else {
      # If only one segment, its length is the whole branch
      segment_lengths <- branch_length
    }
    
    # Accumulate time in each state
    for (j in seq_along(map)) {
      state <- names(map)[j]
      time_in_states[state] <- time_in_states[state] + segment_lengths[j]
    }
  }
  
  # Calculate transitions and state times at nodes
  node_indices <- (length(tree$tip.label) + 1):length(tree$edge[,1])
  for (node in node_indices) {
    parent_edges <- which(tree$edge[,2] == node)
    child_edges <- which(tree$edge[,1] == node)
    
    # States at the end of parent branches
    parent_states <- sapply(parent_edges, function(edge) {
      map <- tree$maps[[edge]]
      if (length(map) > 0) tail(names(map), n = 1) else NA
    })
    
    # States at the start of child branches
    child_states <- sapply(child_edges, function(edge) {
      map <- tree$maps[[edge]]
      if (length(map) > 0) head(names(map), n = 1) else NA
    })
    
    # If all parent states are the same and different from child states, count a transition
    if (length(unique(parent_states)) == 1 && !is.na(parent_states[1])) {
      for (child_state in unique(child_states)) {
        if (parent_states[1] != child_state && !is.na(child_state)) {
          transition_matrix[parent_states[1], child_state] <- 
            transition_matrix[parent_states[1], child_state] + 1
        }
      }
    }
  }
  
  # Print the transition matrix
  print(transition_matrix)
  trans_df <- as.data.frame(as.table(transition_matrix))
  names(trans_df) <- c("state 1", "state 2", "count")
  trans_df <- trans_df[trans_df$count > 0, ]
  trans_df$count <- as.integer(trans_df$count)
  cat("List of all transitions and their counts:\n")
  apply(trans_df, 1, function(x) cat(sprintf('%s -> %s, %s\n', x[1], x[2], x[3])))
  
  # Calculate total time and proportions
  total_time <- sum(time_in_states)
  proportions <- time_in_states / total_time
  
  # Print time spent in each state and proportions
  cat("\nTime spent in each state (absolute and proportion):\n")
  print(data.frame(Time = time_in_states, Proportion = proportions))
  
  #convert the transition details
  trans_df[] <- lapply(trans_df, function(x) if(is.factor(x)) as.character(x) else x)
  
  # Return all results as a list
  return(list(transition_matrix = transition_matrix, transition_details = trans_df, time_in_states = data.frame(time = time_in_states, proportion = proportions)))
}

#' Annotate state transitions with evolutionary rate information
#'
#' Enhances a transition table (typically from a SIMMAP transition analysis) by appending
#' evolutionary rate estimates from a list of variance-covariance (VCV) matrices. Rates are summarized
#' as average variances from the diagonal of each state's VCV matrix.
addRatesToTransitions <- function(transition_details, VCVs) {
  # Verify that the necessary components are present
  if (is.null(transition_details) || !is.data.frame(transition_details)) {
    stop("transition_details must be a data frame.")
  }
  if (is.null(VCVs) || !is.list(VCVs)) {
    stop("VCVs must be a list of matrices.")
  }
  
  # Check if all states in transition_details are available in VCVs
  all_states <- unique(c(transition_details$`state 1`, transition_details$`state 2`))
  if (!all(all_states %in% names(VCVs))) {
    stop("Not all states in transition_details are present in VCVs.")
  }
  
  # Function to calculate the mean variance of the diagonal elements of a matrix
  getAverageVariance <- function(state) {
    if (exists(state, where = VCVs) && is.matrix(VCVs[[state]])) {
      return(mean(diag(VCVs[[state]])))
    } else {
      return(NA)  # Return NA if the state does not exist or is not properly formatted
    }
  }
  
  # Map average variances for each state based on VCVs
  transition_details$rate_1 <- sapply(transition_details$`state 1`, getAverageVariance)
  transition_details$rate_2 <- sapply(transition_details$`state 2`, getAverageVariance)
  
  # Calculate the rate delta
  transition_details$rate_delta <- transition_details$rate_2 - transition_details$rate_1
  
  # Determine the direction of rate change
  transition_details$rate_change <- ifelse(transition_details$rate_delta > 0, "increase",
                                           ifelse(transition_details$rate_delta < 0, "decrease", "no change"))
  
  # Calculate percentage rate change
  transition_details$percentage_change <- with(transition_details, (rate_delta / rate_1) * 100)
  # Handling potential division by zero
  transition_details$percentage_change[is.infinite(transition_details$percentage_change)] <- NA
  
  # Calculate log ratio if both rates are greater than 0
  transition_details$log_ratio <- with(transition_details, ifelse(rate_1 > 0 & rate_2 > 0,
                                                                  log(rate_2 / rate_1), NA))
  
  # Return the updated dataframe
  return(transition_details)
}

#' Summarize directional evolutionary rate changes across transitions
#'
#' Computes a high-level summary of how evolutionary rates change across discrete state transitions,
#' using a SIMMAP tree and corresponding VCV matrices from a search object.
rateSummary <- function(search_object) {
  # Check for necessary components
  if (is.null(search_object$tree_no_uncertainty_untransformed) || is.null(search_object$VCVs)) {
    stop("Each object must contain a 'tree_no_uncertainty_untransformed' and 'VCVs' attributes.")
  }
  
  # Process the transitions in the tree
  transition_details <- countNodeTransitionsSimmap(search_object$tree_no_uncertainty_untransformed)$transition_details
  
  # Add variance rates to the transition details
  updated_transition_details <- addRatesToTransitions(transition_details, search_object$VCVs)
  
  # Summarize rate changes and return
  rate_change_summary <- table(updated_transition_details$rate_change)
  
  return(rate_change_summary)
}

#' Convert absolute counts to relative frequencies
#' Calculates the proportion of each element in a numeric vector relative to the total sum.
calculateRelativeFrequencies <- function(item) {
  total <- sum(item)
  rel_freq <- item / total
  return(rel_freq)
}

#' Return original counts without transformation
#'
#' Simply returns the input vector of counts without modification. Likely used as a placeholder
#' or for interface consistency with other functions (e.g., `calculateRelativeFrequencies()`).
calculateRelativeCounts <- function(item) {
  total <- sum(item)
  rel_count <- item
  return(rel_count)
}

#' Summarize directional rate changes across transitions
#'
#' Computes descriptive statistics (mean, variance, standard deviation) separately for
#' evolutionary rate increases and decreases across regime transitions. Focuses on
#' `rate_delta` values and the categorical `rate_change` direction.
summarizeRateChanges <- function(data) {
  # Check if the necessary column exists
  if (!("rate_delta" %in% names(data)) || !("rate_change" %in% names(data))) {
    stop("Data must contain 'rate_delta' and 'rate_change' columns.")
  }
  
  # Subset data for increases and decreases
  increase_data <- subset(data, rate_change == "increase")
  decrease_data <- subset(data, rate_change == "decrease")
  
  # Calculate mean, variance, and standard deviation for increases
  increase_mean <- mean(increase_data$rate_delta, na.rm = TRUE)
  increase_variance <- var(increase_data$rate_delta, na.rm = TRUE)
  increase_sd <- sqrt(increase_variance)  # Standard deviation is the square root of variance
  
  # Calculate mean, variance, and standard deviation for decreases
  decrease_mean <- mean(decrease_data$rate_delta, na.rm = TRUE)
  decrease_variance <- var(decrease_data$rate_delta, na.rm = TRUE)
  decrease_sd <- sqrt(decrease_variance)
  
  # Create a summary dataframe
  summary_df <- data.frame(
    Change = c("Increase", "Decrease"),
    Mean = c(increase_mean, decrease_mean),
    Variance = c(increase_variance, decrease_variance),
    SD = c(increase_sd, decrease_sd)
  )
  
  # Return the summary dataframe
  return(summary_df)
}

#' Summarize rate changes and test differences between increase vs. decrease
#'
#' Provides a detailed summary of rate changes including descriptive statistics and
#' results of t-tests comparing the magnitude of increases and decreases.
summarizeRateChangesWithTest <- function(data) {
  if (!("rate_delta" %in% names(data)) || !("rate_change" %in% names(data))) {
    stop("Data must contain 'rate_delta' and 'rate_change' columns.")
  }
  
  # Subset data for increases and decreases
  increase_data <- subset(data, rate_change == "increase")
  decrease_data <- subset(data, rate_change == "decrease")
  
  # Calculate mean, variance, and standard deviation for increases
  increase_mean <- mean(increase_data$rate_delta, na.rm = TRUE)
  increase_variance <- var(increase_data$rate_delta, na.rm = TRUE)
  increase_sd <- sqrt(increase_variance)
  
  # Calculate mean, variance, and standard deviation for decreases
  decrease_mean <- mean(decrease_data$rate_delta, na.rm = TRUE)
  decrease_variance <- var(decrease_data$rate_delta, na.rm = TRUE)
  decrease_sd <- sqrt(decrease_variance)
  
  # Perform t-tests on raw and absolute rate deltas
  t_test_raw <- t.test(increase_data$rate_delta, decrease_data$rate_delta, var.equal = FALSE)
  t_test_abs <- t.test(abs(increase_data$rate_delta), abs(decrease_data$rate_delta), var.equal = FALSE)
  
  # Create a summary dataframe for statistics
  summary_stats <- data.frame(
    Change = c("Increase", "Decrease"),
    Mean = c(increase_mean, decrease_mean),
    Variance = c(increase_variance, decrease_variance),
    SD = c(increase_sd, decrease_sd)
  )
  
  # Create a list for t-test results
  t_test_results <- list(
    raw = list(
      statistic = t_test_raw$statistic,
      p_value = t_test_raw$p.value
    ),
    absolute = list(
      statistic = t_test_abs$statistic,
      p_value = t_test_abs$p.value
    )
  )
  
  # Return all results as a list containing both the summary stats and t-test results
  return(list(
    summary_stats = summary_stats,
    t_test_results = t_test_results
  ))
}


#' Identify maximum regime ages and annotate with rate shift direction
#'
#' Extracts the oldest node associated with each unique regime (state) in a SIMMAP tree
#' and annotates each with the direction and magnitude of evolutionary rate change
#' (increase, decrease, or root). Optionally filters results by direction.
extractMaxAgeOfRegimesWithRateChanges <- function(sub_object, filter_by = NULL) {
    # Ensure the sub_object contains necessary components
    if (is.null(sub_object$tree_no_uncertainty_untransformed) || is.null(sub_object$VCVs)) {
      stop("The sub_object must contain a tree and VCVs.")
    }
    
    # Extract the SIMMAP tree
    simmap_tree <- sub_object$tree_no_uncertainty_untransformed
    
    # Step 1: Generate transition details using countNodeTransitionsSimmap and addRatesToTransitions
    transition_details <- addRatesToTransitions(
      transition_details = countNodeTransitionsSimmap(simmap_tree)$transition_details,
      VCVs = sub_object$VCVs
    )
    
    # Convert SIMMAP tree to phylo for compatibility with dispRity::tree.age
    phylo_tree <- as.phylo(simmap_tree)
    
    # Step 2: Get the ages of all nodes using dispRity::tree.age
    ages_info <- dispRity::tree.age(phylo_tree, order = 'past')
    ages_info <- ages_info[!ages_info$ages == 0,]  # Exclude tip ages
    
    # Step 3: Get the states for all nodes using phytools::getStates (using SIMMAP tree)
    states <- phytools::getStates(simmap_tree, type = 'both')
    
    # Step 4: Initialize an empty list to store results
    results <- list()
    
    # Iterate over unique states to find the maximum age node
    for (state in unique(states)) {
      nodes_in_state <- names(states[states == state])
      ages_in_state <- ages_info$ages[ages_info$elements %in% nodes_in_state]
      
      if (length(ages_in_state) > 0) {
        max_age <- max(ages_in_state, na.rm = TRUE)
        max_age_index <- which.max(ages_in_state)
        max_age_node <- nodes_in_state[max_age_index]
      } else {
        max_age <- NA
        max_age_node <- NA
      }
      
      # Append max age and node info to the results list
      results[[state]] <- c(state = state, max_age_node = max_age_node, max_age = max_age)
    }
    
    # Construct the data frame from the list of results
    max_ages_df <- do.call(rbind, results)
    max_ages_df <- as.data.frame(max_ages_df, stringsAsFactors = FALSE)
    max_ages_df$max_age_node <- as.integer(max_ages_df$max_age_node)
    max_ages_df$max_age <- as.numeric(max_ages_df$max_age)
    
    # # Step 5: Add rate change information (increase/decrease) to the max_ages_df
    # # Match the "state 2" from the transition details with the states in max_ages_df
    # rate_changes <- transition_details[, c("state 2", "rate_change")]
    # colnames(rate_changes) <- c("state", "rate_change")
    # 
    # # Ensure state is a character to match with max_ages_df state
    # max_ages_df$state <- as.character(max_ages_df$state)
    # rate_changes$state <- as.character(rate_changes$state)
    # 
    # # Merge the rate change information into the max_ages_df
    # max_ages_df <- merge(max_ages_df, rate_changes, by = "state", all.x = TRUE)
    # 
    # # Assign the label "root" to the root state (state 0)
    # max_ages_df$rate_change[is.na(max_ages_df$rate_change) & max_ages_df$state == "0"] <- "root"
    
    # Step 5: Add rate change information (increase/decrease) to the max_ages_df
    rate_changes <- transition_details[, c("state 2", "rate_change", "rate_delta", "percentage_change")]
    colnames(rate_changes) <- c("state", "rate_change", "rate_delta", "percentage_change")
    
    # Ensure state is a character to match max_ages_df
    max_ages_df$state <- as.character(max_ages_df$state)
    rate_changes$state <- as.character(rate_changes$state)
    
    # Merge the rate change information into max_ages_df
    max_ages_df <- merge(max_ages_df, rate_changes, by = "state", all.x = TRUE)
    
    # Assign "root" to the root state's rate change and handle missing rate data
    max_ages_df$rate_change[is.na(max_ages_df$rate_change) & max_ages_df$state == "0"] <- "root"
    max_ages_df$rate_delta[is.na(max_ages_df$rate_delta) & max_ages_df$state == "0"] <- NA
    max_ages_df$percentage_change[is.na(max_ages_df$percentage_change) & max_ages_df$state == "0"] <- NA
    
    # Ensure numeric format for rate delta and percentage change
    max_ages_df$rate_delta <- as.numeric(max_ages_df$rate_delta)
    max_ages_df$percentage_change <- as.numeric(max_ages_df$percentage_change)
    
    # Step 6: Filter based on user input (increase, decrease), keep the root state always
    if (!is.null(filter_by)) {
      max_ages_df <- max_ages_df[max_ages_df$rate_change == filter_by | max_ages_df$rate_change == "root", ]
    }
    
    return(max_ages_df)
}


#Simulation functions


#' Simulate multivariate traits under Brownian motion with a random covariance matrix
#'
#' Generates a random phylogenetic tree (or uses a provided one), simulates multivariate
#' trait data under a BM1 model using a randomly constructed positive-definite covariance
#' matrix, and returns the full simulation.
simulate_traits_BM1 <- function(n_species, n_traits,
                                variance_mean = 0.1, variance_sd = 0.05, 
                                covariance_mean = 0.0, covariance_sd = 0.25, 
                                tree = NULL, max_attempts = 100) {
  if (is.null(tree)) {
    tree <- pbtree(n = n_species)
  }
  
  attempt <- 1
  repeat {
    # Generate random variances (diagonal)
    variances <- abs(rnorm(n_traits, mean = variance_mean, sd = variance_sd))
    
    # Generate lower triangular matrix
    L <- matrix(rnorm(n_traits * n_traits, mean = covariance_mean, sd = covariance_sd), 
                ncol = n_traits)
    L[upper.tri(L)] <- 0
    
    # Construct symmetric covariance matrix
    sigma <- t(L) %*% L
    diag(sigma) <- variances + diag(sigma)
    
    # Check for positive definiteness via eigenvalues
    if (all(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values > 0)) {
      break
    }
    
    attempt <- attempt + 1
    if (attempt > max_attempts) {
      stop("Failed to generate a valid positive definite covariance matrix after multiple attempts.")
    }
  }
  
  # Simulation parameters
  theta <- rep(0, n_traits)
  params <- list(ntraits = n_traits, sigma = sigma, theta = theta)
  
  # Simulate data
  simulated_data <- mvSIM(tree, nsim = 1, model = "BM1", param = params)
  
  return(list(tree = tree, data = simulated_data, covariance_matrix = sigma))
}

#' Simulate multivariate traits under an Early Burst (EB) evolutionary model
#'
#' Generates a random tree and simulates multivariate traits evolving under an EB model,
#' using a randomly constructed positive-definite covariance matrix and a user-defined
#' or default \code{beta} parameter controlling the acceleration or deceleration of rates over time.
simulate_traits_EB <- function(n_species, n_traits,
                               variance_mean = 0.1, variance_sd = 0.05,
                               covariance_mean = 0, covariance_sd = 0.05,
                               beta = 0.1) {  # Default beta value set to 0.1
  # Generating a random tree with specified number of species
  tree <- pbtree(n = n_species)
  
  # Generating random variances (diagonal elements)
  variances <- abs(rnorm(n_traits, mean = variance_mean, sd = variance_sd))  # Ensure non-negative values
  
  # Create a lower triangular matrix L from normal distribution
  L <- matrix(rnorm(n_traits * n_traits, mean = covariance_mean, sd = covariance_sd), ncol = n_traits)
  L[upper.tri(L)] <- 0  # Zero out upper triangle to make it lower triangular
  
  # Constructing a positive definite covariance matrix
  sigma <- t(L) %*% L  # Multiply L by its transpose
  diag(sigma) <- variances + diag(sigma)  # Add variances to the diagonal elements
  
  # Parameters setup for the simulation
  theta <- rep(0, n_traits)  # Zero mean for the ancestral state
  params <- list(ntraits = n_traits, sigma = sigma, theta = theta, beta = beta)
  
  # Simulate the dataset using the Early Burst model
  simulated_data <- mvSIM(tree, nsim = 1, model = "EB", param = params)
  
  # Return a list containing the tree, the simulated data, the covariance matrix, and the beta value
  return(list(tree = tree, data = simulated_data, covariance_matrix = sigma, beta = beta))
}

#' Plot normalized shift frequencies over time across multiple evolutionary models
#'
#' Processes a list of model output objects containing SIMMAP trees, extracts the timing
#' of regime shifts, normalizes them by lineage richness over time, and visualizes these 
#' dynamics using smoothed spline plots. Produces one curve per model and overlays average 
#' and median trends.
processTreeDataFromObjects <- function(object_list, bin_size, slide_step, legend=T, xlim=NULL) {
  if (!is.list(object_list)) {
    stop("The input 'object_list' must be a list of objects containing trees.")
  }
  
  if (!is.numeric(bin_size) || bin_size <= 0) {
    stop("bin_size must be a positive number.")
  }
  
  if (!is.numeric(slide_step) || slide_step <= 0) {
    stop("slide_step must be a positive number.")
  }
  
  results <- list()
  spline_data <- list() # List to store individual spline data
  average_data <- list(x = numeric(), y = numeric()) # List to store average spline data
  median_data <- list(x = numeric(), y = numeric()) # List to store median spline data
  
  
  colors <- viridis(length(object_list)) # Assign a unique color to each tree
  plot_labels <- names(object_list) # Use the names from the list as labels
  
  x_all <- numeric() # Collect all x-values for averaging and plotting limits
  y_all <- list() # Collect all y-values for averaging
  
  for (i in seq_along(object_list)) {
    object <- object_list[[i]]
    if (!is.null(object$tree_no_uncertainty_untransformed)) {
      tree <- object$tree_no_uncertainty_untransformed
      
      # Extract maximum ages of regimes
      ages <- extractMaxAgeOfRegimes(simmap_tree = tree)$max_age
      counts <- countAgesInBins(ages = ages, bin_size = bin_size, slide_step = slide_step)
      
      # Calculate average lineages
      lineages <- calculateAverageLineagesFromTree(tree, bin_size = bin_size, slide_step = slide_step)
      
      # Store results in a list
      results[[i]] <- list(counts = counts, lineages = lineages)
      
      # Prepare for spline plot
      x <- counts$bin_midpoint
      #if ()
      y <- counts$count / log(lineages$average)
      
      #plot(counts$count ~ lineages$average)
      
      # Collect data for averaging
      x_all <- unique(c(x_all, x))
      y_all[[i]] <- approx(x, y, x_all, rule = 2)$y # Interpolate/extrapolate to common x-values
      
      # Calculate spline for individual data
      spline_fit <- smooth.spline(x = x, y = y)
      spline_data[[i]] <- list(x = spline_fit$x, y = spline_fit$y, type = 'l', col = colors[i], lwd = 2)
    } else {
      warning(paste("Tree not found in object at index", i))
      results[[i]] <- NULL
    }
  }
  
  # Calculate average spline from interpolated y-values
  if (length(y_all) > 0) {
    average_data$y <- rowMeans(do.call(cbind, y_all), na.rm = TRUE)
    average_data$x <- x_all
    average_spline <- smooth.spline(average_data$x, average_data$y)
  }
  
  # Calculate median spline from interpolated y-values
  if (length(y_all) > 0) {
    median_data$y <- matrixStats::rowMedians(do.call(cbind, y_all), na.rm = TRUE)
    median_data$x <- x_all
    median_spline <- smooth.spline(median_data$x, median_data$y)
  }
  
  # Plot all splines together
  if (length(spline_data) > 0) {
    plot(NA, NA, type = 'n', xlim = xlim, ylim = range(sapply(spline_data, function(s) range(s$y))),
         xlab = "Time (bin midpoint)", ylab = "Normalized Event Frequency", main = "Comparative Splines for Shifts")
    for (spline in spline_data) {
      lines(spline$x, spline$y, type = spline$type, col = make.transparent(spline$col, 0.25), lwd = spline$lwd)
    }
    # Plot the average and median splines
    lines(median_spline, type = 'l', col = "grey", lwd = 5, lty=1, lend=2)
    lines(average_spline, type = 'l', col = "black", lwd = 5, lty=1, lend=2)
    if(legend){
      legend("topright", c(plot_labels, "Average", "Median"), col = c(colors, "black", "grey"), lwd = 2, title = "Legend", cex=0.5)
    }
  } else {
    print("No valid data to plot.")
  }
  #old xlim was rev(range(x_all))
  return(list(results = results, plots = spline_data))
}


#' Visualize regime shift frequencies over time with statistical filtering and uncertainty
#'
#' Processes a list of model objects containing SIMMAP trees and extracts regime shift
#' timing. Allows optional filtering by rate change direction (increase or decrease),
#' normalization by lineage richness or branch length, and visualizes central trends with 
#' uncertainty (HDIs).
processTreeDataFromObjects.filter <- function(object_list, bin_size, slide_step, legend = TRUE,
                                              xlim = NULL, ylim = NULL, filter_by = NULL,
                                              normalization = "log", hdi_prob = 0.95,
                                              central_tendency = "mean", legend_position = "topright",
                                              bootstrap_hdi = FALSE, n_boot = 100,
                                              y_axis_side = "left",           # "left" or "right"
                                              show_x_axis = TRUE, 
                                              show_y_axis = TRUE,
                                              title = NULL) {           # logical
  
  if (!is.list(object_list)) stop("The input 'object_list' must be a list.")
  if (!is.numeric(bin_size) || bin_size <= 0) stop("bin_size must be a positive number.")
  if (!is.numeric(slide_step) || slide_step <= 0) stop("slide_step must be a positive number.")
  if (!normalization %in% c("log", "raw", "absolute", "branch_length", "log_branch_length"))
    stop("Invalid 'normalization'.")
  if (!is.numeric(hdi_prob) || hdi_prob <= 0 || hdi_prob >= 1)
    stop("'hdi_prob' must be between 0 and 1.")
  if (!central_tendency %in% c("mean", "median", "mode"))
    stop("'central_tendency' must be one of 'mean', 'median', or 'mode'.")
  if (!y_axis_side %in% c("left", "right"))
    stop("'y_axis_side' must be either 'left' or 'right'.")
  if (!is.logical(show_x_axis))
    stop("'show_x_axis' must be TRUE or FALSE.")
  
  get_hdi_from_density <- function(values, prob) {
    d <- density(values, na.rm = TRUE, n = 512)
    dx <- d$x
    dy <- d$y
    ord <- order(dy, decreasing = TRUE)
    cum_prob <- cumsum(dy[ord]) / sum(dy)
    idx <- which(cum_prob <= prob)
    range(dx[ord][idx])
  }
  
  get_mode_from_density <- function(values) {
    d <- density(values, na.rm = TRUE, n = 512)
    d$x[which.max(d$y)]
  }
  
  process_filtered_type <- function(filter_type, color) {
    x_all <- numeric()
    y_all <- list()
    branch_lengths <- NULL
    
    for (i in seq_along(object_list)) {
      object <- object_list[[i]]
      if (!is.null(object$tree_no_uncertainty_untransformed)) {
        tree <- object$tree_no_uncertainty_untransformed
        
        if (is.null(branch_lengths) && normalization %in% c("branch_length", "log_branch_length")) {
          branch_lengths <- calculateTotalBranchLengthInBins(tree, bin_size, slide_step)
        }
        
        ages <- if (is.null(filter_type)) {
          extractMaxAgeOfRegimes(simmap_tree = tree)$max_age
        } else {
          extractMaxAgeOfRegimesWithRateChanges(sub_object = object, filter_by = filter_type)$max_age
        }
        
        counts <- countAgesInBins(ages, bin_size, slide_step)
        
        if (!normalization %in% c("branch_length", "log_branch_length")) {
          lineages <- calculateAverageLineagesFromTree(tree, bin_size, slide_step)
        }
        
        y <- switch(normalization,
                    "log" = counts$count / log(lineages$average),
                    "absolute" = counts$count / lineages$average,
                    "branch_length" = counts$count / branch_lengths$total_branch_length,
                    "log_branch_length" = counts$count / log(branch_lengths$total_branch_length),
                    "raw" = counts$count)
        
        x <- counts$bin_midpoint
        x_all <- unique(c(x_all, x))
        y_all[[i]] <- approx(x, y, x_all, rule = 2)$y
      }
    }
    
    x_all <- sort(x_all)
    y_matrix <- do.call(cbind, lapply(y_all, function(y) approx(x_all, y, x_all, rule = 2)$y))
    
    center_y <- switch(central_tendency,
                       "mean" = rowMeans(y_matrix, na.rm = TRUE),
                       "median" = matrixStats::rowMedians(y_matrix, na.rm = TRUE),
                       "mode" = apply(y_matrix, 1, get_mode_from_density))
    
    if (bootstrap_hdi) {
      hdi_bounds <- t(apply(y_matrix, 1, function(values) {
        boots <- replicate(n_boot, {
          samp <- sample(values, replace = TRUE)
          get_hdi_from_density(samp, prob = hdi_prob)
        })
        apply(boots, 1, median, na.rm = TRUE)
      }))
    } else {
      hdi_bounds <- t(apply(y_matrix, 1, get_hdi_from_density, prob = hdi_prob))
    }
    
    lower_spline <- smooth.spline(x_all, hdi_bounds[, 1])
    upper_spline <- smooth.spline(x_all, hdi_bounds[, 2])
    center_spline <- smooth.spline(x_all, center_y)
    
    list(center = center_spline, lower = lower_spline, upper = upper_spline,
         color = color, x_all = x_all, y_matrix  = y_matrix)
  }
  
  color_decrease <- "steelblue"
  color_increase <- "firebrick"
  
  if (!is.null(filter_by) && filter_by == "both") {
    decrease <- process_filtered_type("decrease", color_decrease)
    increase <- process_filtered_type("increase", color_increase)
    plot_title <- paste0(central_tendency, " & ", round(hdi_prob * 100), "% HDI for Rate Increases & Decreases")
  } else {
    result <- process_filtered_type(filter_by, "black")
    filter_text <- if (is.null(filter_by)) "Model Shifts" else filter_by
    plot_title <- paste0(central_tendency, " & ", round(hdi_prob * 100), "% HDI for ", filter_text)
  }
  
  if (!is.null(title)) {plot_title = title}
  
  ylab_text <- switch(normalization,
                      "log" = "Event Frequency - log(lineages)",
                      "absolute" = "Event Frequency - lineages",
                      "branch_length" = "Event Frequency - total branch length",
                      "log_branch_length" = "Event Frequency - log(total branch length)",
                      "raw" = "Shift Counts")
  
  all_y <- if (!is.null(filter_by) && filter_by == "both") {
    c(decrease$lower$y, decrease$upper$y, increase$lower$y, increase$upper$y)
  } else {
    c(result$lower$y, result$upper$y)
  }
  if (is.null(ylim)) ylim <- range(all_y, na.rm = TRUE)
  
  # --- NEW BOX/AXIS HANDLING -------------------------------------------------
  bty_plot  <- "n"                            # no box: we’ll add only what we need
  xaxt_plot <- ifelse(show_x_axis, "s", "n")  # default x-axis behaviour
  ylab_plot <- ifelse(y_axis_side == "left", ylab_text, "")
  
  # Dynamically suppress x-axis label if show_x_axis is FALSE
  xlab_plot <- ifelse(show_x_axis, "Time (bin midpoint)", "")
  
  plot(NA, NA, type = 'n', xlim = xlim, ylim = ylim,
       xlab = xlab_plot, ylab = ylab_plot, main = plot_title,
       bty = bty_plot, xaxt = xaxt_plot, yaxt = 'n')
  
  ## y-axis ticks (only if requested)
  if (show_y_axis) {
    axis(ifelse(y_axis_side == "left", 2, 4))
  }
  ## Always show the y-axis label
  mtext(ylab_text, side = ifelse(y_axis_side == "left", 2, 4), line = 2.5)
  
  ## ensure bottom baseline is present even if x-axis suppressed
  if (!show_x_axis) axis(1, labels = FALSE, tick = FALSE)
  # ---------------------------------------------------------------------------
  
  if (!is.null(filter_by) && filter_by == "both") {
    polygon(c(decrease$upper$x, rev(decrease$lower$x)),
            c(decrease$upper$y, rev(decrease$lower$y)),
            col = adjustcolor(decrease$color, alpha.f = 0.2), border = NA)
    lines(decrease$center, col = decrease$color, lwd = 3)
    
    polygon(c(increase$upper$x, rev(increase$lower$x)),
            c(increase$upper$y, rev(increase$lower$y)),
            col = adjustcolor(increase$color, alpha.f = 0.2), border = NA)
    lines(increase$center, col = increase$color, lwd = 3)
    
    if (legend) {
      legend(legend_position, legend = c("Rate Decrease", "Rate Increase"),
             col = c(color_decrease, color_increase), lwd = 2, bty = "n")
    }
    
    return(list(
      decrease = decrease,
      increase = increase,
      hdi_prob = hdi_prob,
      central_tendency = central_tendency,
      bootstrap_hdi = bootstrap_hdi,
      n_boot = n_boot
    ))
    
  } else {
    polygon(c(result$upper$x, rev(result$lower$x)),
            c(result$upper$y, rev(result$lower$y)),
            col = adjustcolor(result$color, alpha.f = 0.5), border = NA)
    lines(result$center, col = result$color, lwd = 2)
    
    if (legend) {
      legend(legend_position, legend = c(central_tendency, paste0(round(hdi_prob * 100), "% HDI")),
             col = c("black", "grey70"), lwd = c(2, NA), pch = c(NA, 15),
             pt.cex = 2, bty = "n")
    }
    
    return(list(
      result = result,
      hdi_prob = hdi_prob,
      central_tendency = central_tendency,
      bootstrap_hdi = bootstrap_hdi,
      n_boot = n_boot
    ))
  }
}

#' Compute temporal variability metrics in sliding time windows
#'
#' Calculates a variety of variability metrics (e.g., standard deviation, coefficient of variation, 
#' entropy, etc.) across overlapping sliding windows along a time axis in a dataset.
calculate_variability <- function(data, time_column, value_column, bin_size, slide_step, metric = "sd", num_bins = 10, lag = 1) {
  
  # Ensure bin_size and slide_step are positive numbers
  if (!is.numeric(bin_size) || bin_size <= 0) {
    stop("bin_size must be a positive number.")
  }
  
  if (!is.numeric(slide_step) || slide_step <= 0) {
    stop("slide_step must be a positive number.")
  }
  
  # Valid metrics now include 'sd', 'cv', 'entropy', 'variance', 'first_diff_sd', 'total_variation', and 'acf'
  if (!metric %in% c("sd", "cv", "entropy", "variance", "first_diff_sd", "total_variation", "acf")) {
    stop("Invalid metric. Choose from 'sd', 'cv', 'entropy', 'variance', 'first_diff_sd', 'total_variation', or 'acf'.")
  }
  
  # Get the range of time
  time_range <- range(data[[time_column]], na.rm = TRUE)
  
  # Initialize lists to store the results
  bin_midpoints <- numeric()
  variability <- numeric()
  
  # Start the sliding window process
  for (start_time in seq(time_range[1], time_range[2] - bin_size, by = slide_step)) {
    
    # Define the window range
    end_time <- start_time + bin_size
    
    # Subset the data within this window
    window_data <- subset(data, data[[time_column]] >= start_time & data[[time_column]] < end_time)
    
    # Remove any rows with NA values in the variable of interest
    window_data <- window_data[!is.na(window_data[[value_column]]), ]
    
    # Only calculate if there are enough points
    if (nrow(window_data) > 1) {
      bin_midpoint <- mean(c(start_time, end_time))
      
      if (metric == "sd") {
        # Standard deviation
        variability_value <- sd(window_data[[value_column]], na.rm = TRUE)
        
      } else if (metric == "cv") {
        # Coefficient of Variation (CV = sd / mean)
        window_mean <- mean(window_data[[value_column]], na.rm = TRUE)
        window_sd <- sd(window_data[[value_column]], na.rm = TRUE)
        
        # Avoid division by zero in case of a mean close to 0
        if (window_mean != 0) {
          variability_value <- window_sd / window_mean
        } else {
          variability_value <- NA  # Undefined CV when mean is 0
        }
        
      } else if (metric == "entropy") {
        # Shannon entropy
        # Bin the isotope values into discrete categories
        window_values <- window_data[[value_column]]
        
        # Create a histogram of values in the specified number of bins
        hist_counts <- hist(window_values, breaks = num_bins, plot = FALSE)$counts
        
        # Calculate the probabilities (relative frequencies)
        prob <- hist_counts / sum(hist_counts)
        
        # Remove zero probabilities to avoid log(0)
        prob <- prob[prob > 0]
        
        # Calculate Shannon entropy
        entropy_value <- -sum(prob * log(prob))
        variability_value <- entropy_value
        
      } else if (metric == "variance") {
        # Variance
        variability_value <- var(window_data[[value_column]], na.rm = TRUE)
        
      } else if (metric == "first_diff_sd") {
        # Standard deviation of the first differences
        differences <- diff(window_data[[value_column]], lag = 1)
        
        # Only compute if we have more than one difference
        if (length(differences) > 1) {
          variability_value <- sd(differences, na.rm = TRUE)
        } else {
          variability_value <- NA
        }
        
      } else if (metric == "total_variation") {
        # Total Variation (sum of absolute first differences)
        differences <- diff(window_data[[value_column]], lag = 1)
        
        # Only compute if we have more than one difference
        if (length(differences) > 1) {
          variability_value <- sum(abs(differences), na.rm = TRUE)
        } else {
          variability_value <- NA
        }
        
      } else if (metric == "acf") {
        # Autocorrelation at a specified lag
        acf_value <- acf(window_data[[value_column]], lag.max = lag, plot = FALSE)$acf[lag + 1]
        variability_value <- acf_value
      }
      
      # Store the results
      bin_midpoints <- c(bin_midpoints, bin_midpoint)
      variability <- c(variability, variability_value)
    }
  }
  
  # Return the results as a data frame
  return(data.frame(bin_midpoint = bin_midpoints, variability = variability))
}

#' Apply Fisher's Z-transformation to correlation coefficients
#'
#' Transforms a numeric vector of Pearson correlation coefficients into Z-scores using 
#' Fisher's Z-transformation. This transformation stabilizes the variance of correlations
#' and is commonly used in meta-analysis and hypothesis testing.
fisher_z_transform <- function(correlations) {
  # Check if all values are within the valid range for correlation coefficients
  if (any(correlations < -1 | correlations > 1)) {
    stop("All correlation coefficients must be between -1 and 1.")
  }
  
  # Apply the Fisher Z-transformation
  z_transformed <- 0.5 * log((1 + correlations) / (1 - correlations))
  
  return(z_transformed)
}

#Helper functions for post-hoc analysis of shift patterns
{
#' Collapse phylogeny by monophyletic state groups
#'
#' Simplifies a SIMMAP phylogenetic tree by collapsing monophyletic groups of tips 
#' that share the same state into a single representative tip, and removing non-monophyletic 
#' state groups entirely.
collapsePhylogenyByStates <- function(simmap_tree) {
  # Ensure the tree is of the correct class
  if (!inherits(simmap_tree, "simmap")) {
    stop("Input tree must be of class 'simmap'.")
  }
  
  print("Starting to process the tree based on state information.")
  
  # Extract states from the tree tips
  tip_states <- getStates(simmap_tree, type = 'tips')
  unique_states <- unique(tip_states)
  
  # Clone the original tree to modify
  modified_tree <- simmap_tree
  
  # Iterate over each unique state
  for (state in unique_states) {
    tips_with_state <- names(tip_states[tip_states == state])
    print(paste("Processing state:", state, "- Number of tips:", length(tips_with_state)))
    
    # Check if the group of tips is monophyletic
    if (is.monophyletic(phy = modified_tree, tips = tips_with_state)) {
      print(paste("State", state, "is monophyletic. Collapsing to one tip."))
      # If monophyletic, drop all tips but one
      tips_to_drop <- tips_with_state[-1]  # Keep the first tip
      modified_tree <- drop.tip(modified_tree, tips_to_drop)
      # Rename the remaining tip to reflect the state name
      modified_tree$tip.label[modified_tree$tip.label == tips_with_state[1]] <- state
    } else {
      # If not monophyletic, remove all tips with that state
      print(paste("State", state, "is not monophyletic. Removing all associated tips."))
      modified_tree <- drop.tip(modified_tree, tips_with_state)
    }
  }
  
  print("Finished processing the tree.")
  
  return(modified_tree)
}

#' Two-pass filter for collapsing phylogeny by state
#'
#' Applies a two-pass filtering procedure to simplify a SIMMAP phylogenetic tree 
#' by first collapsing monophyletic state groups and then removing non-monophyletic ones.z
collapsePhylogenyByStates.2 <- function(simmap_tree) {
  # Ensure the tree is of the correct class
  if (!inherits(simmap_tree, "simmap")) {
    stop("Input tree must be of class 'simmap'.")
  }
  
  print("Starting to process the tree based on state information.")
  
  # Extract states from the tree tips
  tip_states <- getStates(simmap_tree, type = 'tips')
  unique_states <- unique(tip_states)
  
  # Clone the original tree to modify
  modified_tree <- simmap_tree
  
  # First pass: collapse monophyletic states
  for (state in unique_states) {
    tips_with_state <- names(tip_states[tip_states == state])
    if (is.monophyletic(phy = modified_tree, tips = tips_with_state)) {
      print(paste("State", state, "is monophyletic. Collapsing to one tip."))
      # If monophyletic, drop all tips but one
      tips_to_drop <- tips_with_state[-1]  # Keep the first tip
      modified_tree <- drop.tip(modified_tree, tips_to_drop)
      # Rename the remaining tip to reflect the state name
      remaining_tip <- tips_with_state[1]
      modified_tree$tip.label[modified_tree$tip.label == remaining_tip] <- state
    }
  }
  
  # Update states after modifications
  tip_states <- getStates(modified_tree, type = 'tips')
  unique_states <- unique(tip_states)
  
  # Second pass: remove non-monophyletic states
  for (state in unique_states) {
    tips_with_state <- names(tip_states[tip_states == state])
    if (!is.monophyletic(phy = modified_tree, tips = tips_with_state)) {
      print(paste("State", state, "is not monophyletic after first pass. Removing all associated tips."))
      modified_tree <- drop.tip(modified_tree, tips_with_state)
    }
  }
  
  print("Finished processing the tree.")
  
  return(modified_tree)
}

#' Sample One Tip per State from a SIMMAP Tree
#'
#' Randomly samples one tip per unique discrete state from a SIMMAP-formatted phylogenetic tree. 
#' Optionally repeats the sampling process multiple times.
sampleOneTipPerState <- function(simmap_tree, n = 1) {
  # Check input
  if (!inherits(simmap_tree, "simmap")) {
    stop("Input must be a 'simmap' object.")
  }
  
  get_sample <- function() {
    tip_states <- getStates(simmap_tree, type = 'tips')
    unique_states <- unique(tip_states)
    
    tips_to_keep <- character()
    new_tip_labels <- character()
    
    for (state in unique_states) {
      state_tips <- names(tip_states[tip_states == state])
      chosen_tip <- sample(state_tips, 1)
      tips_to_keep <- c(tips_to_keep, chosen_tip)
      new_tip_labels <- c(new_tip_labels, state)
    }
    
    # Drop all but selected tips
    sampled_tree <- drop.tip(simmap_tree, setdiff(simmap_tree$tip.label, tips_to_keep))
    # Rename tips to state names
    sampled_tree$tip.label <- new_tip_labels
    return(sampled_tree)
  }
  
  if (n == 1) {
    return(get_sample())
  } else {
    return(replicate(n, get_sample(), simplify = FALSE))
  }
}

}


#' Plot Densities of Rate Changes with Statistical Comparisons
#'
#' Visualizes the distribution of rate change magnitudes (increases vs. decreases) using smoothed density plots,
#' and performs statistical tests (KS, Wilcoxon) with optional bootstrapping to assess distributional differences.
plotRateChangeDensities <- function(data_input, use_log = FALSE, test_type = "all", text_y = 0.5,
                                    y_nudge = 0, x_nudge = 0, legend_x_nudge = 0,
                                    title = "Rate Change Densities", cores = 4, ...) {
  if (is.data.frame(data_input)) {
    data_list <- list(data_input)
  } else if (is.list(data_input) && all(sapply(data_input, is.data.frame))) {
    data_list <- data_input
  } else {
    stop("Input must be a data frame or list of data frames.")
  }
  
  all_increases <- numeric(0)
  all_decreases <- numeric(0)
  
  for (df in data_list) {
    if (!all(c("rate_delta", "rate_change") %in% names(df))) {
      stop("Each data frame must contain 'rate_delta' and 'rate_change' columns.")
    }
    
    increases <- df$rate_delta[df$rate_change == "increase"]
    decreases <- df$rate_delta[df$rate_change == "decrease"]
    
    if (use_log) {
      increases <- log(abs(increases))
      decreases <- log(abs(decreases))
    } else {
      increases <- abs(increases)
      decreases <- abs(decreases)
    }
    
    all_increases <- c(all_increases, increases)
    all_decreases <- c(all_decreases, decreases)
  }
  
  combined_data <- c(all_increases, all_decreases)
  group_labels <- c(rep(1, length(all_increases)), rep(2, length(all_decreases)))
  
  ks_test_summary <- ""
  wilcox_test_summary <- ""
  modal_summary <- ""
  
  if (test_type == "ks" || test_type == "all") {
    ks_result <- ks.test(all_increases, all_decreases, simulate.p.value = TRUE, B = 10000)
    cat("KS Test Results:\n")
    print(ks_result)
    
    ks_boot_pvalue <- function(data, indices) {
      data1 <- data[indices[group_labels[indices] == 1]]
      data2 <- data[indices[group_labels[indices] == 2]]
      ks.test(data1, data2, simulate.p.value = TRUE, B = 10000)$p.value
    }
    
    boot_results_ks <- boot(data = combined_data, statistic = ks_boot_pvalue, R = 100,
                            strata = group_labels, parallel = 'multicore', ncpus = cores)
    p_values_ks <- boot_results_ks$t
    mean_p_value_ks <- mean(p_values_ks)
    
    cat("Bootstrap KS P-Value Mean:", mean_p_value_ks, "\n")
    
    ks_test_summary <- sprintf("KS Test: D = %.3g,\np = %.3g\nBootstrap Mean p = %.3g",
                               ks_result$statistic, ks_result$p.value, mean_p_value_ks)
  }
  
  if (test_type == "wilcox" || test_type == "all") {
    wilcox_result <- wilcox.test(all_increases, all_decreases, paired = FALSE, exact = FALSE, conf.int = TRUE)
    cat("Regular Wilcoxon Test Results:\n")
    print(wilcox_result)
    
    wilcox_boot_pvalue <- function(data, indices) {
      data1 <- data[indices[group_labels[indices] == 1]]
      data2 <- data[indices[group_labels[indices] == 2]]
      wilcox.test(data1, data2, paired = FALSE, exact = FALSE)$p.value
    }
    
    boot_results_wilcox <- boot(data = combined_data, statistic = wilcox_boot_pvalue, R = 1000, strata = group_labels)
    p_values_wilcox <- boot_results_wilcox$t
    mean_p_value_wilcox <- mean(p_values_wilcox)
    
    cat("Bootstrap Wilcoxon P-Value Mean:", mean_p_value_wilcox, "\n")
    
    wilcox_test_summary <- sprintf("Wilcoxon Test: W = %.3g,\np = %.3g\nBootstrap Mean p = %.3g",
                                   wilcox_result$statistic, wilcox_result$p.value, mean_p_value_wilcox)
  }
  
  density_increases <- density(all_increases)
  density_decreases <- density(all_decreases)
  
  count_increases <- length(all_increases)
  count_decreases <- length(all_decreases)
  frequency_increases <- count_increases / (count_increases + count_decreases)
  frequency_decreases <- count_decreases / (count_increases + count_decreases)
  
  area_increases <- sum(density_increases$y * diff(density_increases$x[1:2]))
  area_decreases <- sum(density_decreases$y * diff(density_decreases$x[1:2]))
  density_increases$y <- density_increases$y / area_increases * frequency_increases
  density_decreases$y <- density_decreases$y / area_decreases * frequency_decreases
  
  
  if (length(all_increases) > 0 && length(all_decreases) > 0) {
    plot(density_increases, main = title, xlab = "Magnitude of Rate Change",
         col = "#377eb8", ylab = 'Relative Frequencies', bty = 'o', ...)
    polygon(density_increases, col = rgb(0, 0, 1, 0.25))
    lines(density_decreases, col = "#e41a1c")
    polygon(density_decreases, col = rgb(1, 0, 0, 0.25))
    box(bty = "L")
    
    abline(v = density_increases$x[which.max(density_increases$y)], col = "#377eb8", lty = 2)
    abline(v = density_decreases$x[which.max(density_decreases$y)], col = "#e41a1c", lty = 2)
    
    text_y <- if (is.null(text_y)) par("usr")[3] + 0.95 * (par("usr")[4] - par("usr")[3]) else text_y
    text(density_increases$x[which.max(density_increases$y)] + 0.5, text_y,
         labels = paste("N: ", count_increases), pos = 3, col = "#377eb8")
    text(density_decreases$x[which.max(density_decreases$y)] - 0.5, text_y,
         labels = paste("N: ", count_decreases), pos = 3, col = "#e41a1c")
    
    modal_x_increases <- density_increases$x[which.max(density_increases$y)]
    modal_x_decreases <- density_decreases$x[which.max(density_decreases$y)]
    ratio_log_space <- exp(modal_x_increases) / exp(modal_x_decreases)
    
    cat("Modal value for increases:", modal_x_increases, "\n")
    cat("Modal value for decreases:", modal_x_decreases, "\n")
    cat("Ratio in linear space:", ratio_log_space, "\n")
    
    modal_summary <- sprintf("Modal (inc): %.*g\nModal (dec): %.*g\nRatio (linear): %.*g",
                             3, modal_x_increases, 3, modal_x_decreases, 3, ratio_log_space)
    
    legend_text <- c()
    if (test_type == "ks" || test_type == "all") {
      legend_text <- c(legend_text, ks_test_summary)
    }
    if (test_type == "wilcox" || test_type == "all") {
      legend_text <- c(legend_text, wilcox_test_summary)
    }
    
    legend("topright", legend = c("Increases", "Decreases"), col = c("#377eb8", "#e41a1c"), lty = 1,
           fill = c(rgb(0, 0, 1, 0.25), rgb(1, 0, 0, 0.25)), title = "Rate Change Type", bty = 'n', inset = c(legend_x_nudge, 0))
    
    mtext(paste(c(legend_text, modal_summary), collapse = "\n"),
          side = 3, line = -3 - y_nudge, adj = 1,
          cex = 0.75, at = par("usr")[2] - 0.5 + x_nudge)
  } else {
    print("No data to plot.")
  }
}

#' Process Transition Rates Across Multiple Phylogenetic Trees
#'
#' Computes transition details and associated evolutionary rate changes for a list of input objects,
#' each containing a SIMMAP tree and a corresponding variance-covariance matrix (VCV).
processTransitionRates <- function(object_list) {
  # This list will store the results
  results_list <- list()
  
  # Loop through each object in the list
  for (obj_name in names(object_list)) {
    obj <- object_list[[obj_name]]
    
    # Perform the transition counting and rate adding
    transition_details <- countNodeTransitionsSimmap(obj$tree_no_uncertainty_untransformed)$transition_details
    results <- addRatesToTransitions(
      transition_details = transition_details,
      VCVs = obj$VCVs
    )
    
    # Store the result with the name of the object
    results_list[[obj_name]] <- results
    print(paste("Processed", obj_name))
  }
  
  return(results_list)
}

# Function to compare two lists and print a summary of the counts
compare_lists <- function(list1, list2) {
  
  # Find the common elements
  common_elements <- intersect(list1, list2)
  
  # Find unique elements in each list
  unique_list1 <- setdiff(list1, list2)
  unique_list2 <- setdiff(list2, list1)
  
  # Create a summary table
  summary_table <- data.frame(
    Element = c(common_elements, unique_list1, unique_list2),
    Source = c(rep("Common", length(common_elements)), 
               rep("Unique to List 1", length(unique_list1)), 
               rep("Unique to List 2", length(unique_list2)))
  )
  
  # Print the summary counts
  cat("Summary of Counts:\n")
  cat("Total elements in list1:", length(list1), "\n")
  cat("Total elements in list2:", length(list2), "\n")
  cat("Common elements:", length(common_elements), "\n")
  cat("Unique elements in list1:", length(unique_list1), "\n")
  cat("Unique elements in list2:", length(unique_list2), "\n")
  
  # Return the summary table
  return(summary_table)
}

#' Plot Dual-Metric Time Series with Independent Y-Axes
#'
#' Visualizes two metrics (e.g., standard deviation, autocorrelation) over time using a dual y-axis layout.
#' Supports smoothed or point-based plotting, customizable axis positioning, and dynamic suppression of plot elements.
plot_with_dual_y_axis <- function(data, x_col, y_col1, y_col2 = NULL, 
                                  metric1 = "sd", metric2 = "acf", 
                                  bin_size = 5, slide_step = 0.1, lag = 2, 
                                  xlim = NULL, ylim1 = NULL, ylim2 = NULL, 
                                  col1 = 'skyblue', col2 = 'rosybrown', 
                                  ylab1 = "SD", ylab2 = "ACF", 
                                  xlab = "Time", main = "Temperature Variability",
                                  pch1 = 16, pch2 = 17, 
                                  plot_type = "points",  # "points" or "spline"
                                  line_col1 = NULL, line_col2 = NULL, 
                                  line_width1 = 1, line_width2 = 1, 
                                  line_type1 = 1, line_type2 = 1,  # Add line type options for each metric
                                  plot_metric2 = TRUE,  # Flag to control second metric plot
                                  metric1_on_right = FALSE,  # Flag to control the side of metric1's y-axis
                                  suppress_x_axis = FALSE, 
                                  suppress_y_axis = FALSE  # New parameter to suppress y-axis,
) {
  
  # Set default line colors to match point colors if not specified
  if (is.null(line_col1)) line_col1 <- col1
  if (is.null(line_col2)) line_col2 <- col2
  
  # Manually set margins regardless of the value of plot_metric2
  par(mar = c(5, 4, 4, 5) + 0.1)  # Bottom, Left, Top, Right margins with space for right axis
  
  # First calculation for the first metric
  tmp1 <- calculate_variability(data, x_col, y_col1, bin_size = bin_size, slide_step = slide_step, metric = metric1, lag = lag)
  
  # Define x-axis limits if not provided
  if (is.null(xlim)) {
    xlim <- c(min(tmp1$bin_midpoint), max(tmp1$bin_midpoint))
  }
  
  # Check if the x-axis is reversed (the first element in xlim is larger than the second)
  is_reversed <- xlim[1] > xlim[2]
  
  # Plot first metric (ACF or other) based on user choice for "points" or "spline"
  if (plot_type == "points") {
    plot(tmp1$variability ~ tmp1$bin_midpoint, xlim = xlim, ylim = ylim1, 
         type = 'b', cex = 0.7, col = col1, 
         ylab = "",  # Suppress default y-axis label
         xlab = ifelse(suppress_x_axis, "", xlab), main = main, 
         pch = pch1, axes = FALSE, lty = line_type1, bty='n')  # Suppress default axes and set line type
  } else {
    plot(tmp1$bin_midpoint, tmp1$variability, type = 'n', xlim = xlim, ylim = ylim1, 
         ylab = "", xlab = ifelse(suppress_x_axis, "", xlab), main = main, axes = FALSE, bty='n')
    lines(smooth.spline(tmp1$bin_midpoint, tmp1$variability), col = line_col1, lwd = line_width1, lty = line_type1)
  }
  
  # Conditionally add x-axis
  if (!suppress_x_axis) {
    axis(1)
  }
  
  # Check if the y-axis should be suppressed
  if (!suppress_y_axis) {
    # Check if the second metric is being plotted
    if (plot_metric2) {
      axis(2, col = col1, col.axis = col1)  # Left-side y-axis for first metric
      mtext(ylab1, side = 2, line = 3, col = col1)  # Left-side y-axis label
    } else {
      # If plot_metric2 is FALSE, allow metric1_on_right option
      if (metric1_on_right) {
        axis(4, col = 'darkgrey', col.axis = 'darkgrey')
        #axis(4, col = col1, col.axis = col1)  # Right-side y-axis for first metric
        #mtext(ylab1, side = 4, line = 3, col = col1)  # Right-side y-axis label
        mtext(ylab1, side = 4, line = 3, col = 'darkgrey')  # Right-side y-axis label
      } else {
        axis(2, col = col1, col.axis = col1)  # Left-side y-axis for first metric
        mtext(ylab1, side = 2, line = 3, col = col1)  # Left-side y-axis label
      }
    }
  }
  
  # Only plot the second metric if the user wants it
  if (plot_metric2) {
    # Second calculation for the second metric
    tmp2 <- calculate_variability(data, x_col, y_col2, bin_size = bin_size, slide_step = slide_step, metric = metric2)
    
    # Overlay second plot with no axes
    par(new = TRUE)
    
    # Plot second metric (SD or other) based on user choice for "points" or "spline"
    if (plot_type == "points") {
      plot(tmp2$variability ~ tmp2$bin_midpoint, xlim = xlim, ylim = ylim2, 
           type = 'b', cex = 0.7, col = col2, 
           axes = FALSE, xlab = "", ylab = "", 
           pch = pch2, lty = line_type2, bty='n')  # Set line type for the second metric
    } else {
      plot(tmp2$bin_midpoint, tmp2$variability, type = 'n', xlim = xlim, ylim = ylim2, 
           axes = FALSE, xlab = "", ylab = "", bty='n')
      lines(smooth.spline(tmp2$bin_midpoint, tmp2$variability), col = line_col2, lwd = line_width2, lty = line_type2)
    }
    
    # Add the second y-axis on the right (use col2 for the axis color)
    if (!suppress_y_axis) {
      axis(4, col = col2, col.axis = col2)  # Right-side y-axis for second metric
      #axis(4, col = 'darkgrey', col.axis = 'darkgrey')  # Right-side y-axis for second metric
      mtext(ylab2, side = 4, line = 3, col = col2)  # Right-side y-axis label
    }
  }
  
  # Add the box around the plot
  #box()
}

#' Plot Isotope Data with Custom Annotations
#'
#' Generates a scatter or line plot of benthic isotope data with optional point annotations
#' using arrows and custom text labels positioned at user-defined angles and offsets.
label_plot_with_annotations <- function(data, data_long, x_coords, 
                                        x_col = "age_tuned", 
                                        y_col = "benthic d18O VPDB CorrAdjusted", 
                                        y_col_long = "ISOBENd18oLOESSsmoothLongTerm", 
                                        plot_cex = 0.1, 
                                        plot_type = 'b', 
                                        plot_xlim = c(45.24011, min(data[[x_col]])), 
                                        plot_ylim = c(5.5, -0.1), 
                                        point_col = 'lightblue',  # Color for points in data
                                        point_pch = 20,           # pch for points in data
                                        long_point_col = 'black',    # Color for points in data_long
                                        long_point_pch = 20,       # pch for points in data_long
                                        y_buffer = 0.3, 
                                        label_cex = 0.8, 
                                        line_color = "red", 
                                        text_color = "blue",
                                        line_width = 1.2, 
                                        line_type = 1, 
                                        text_size = 0.7, 
                                        title = "Plot Title", 
                                        xlab = "X Axis Label", 
                                        ylab = "Y Axis Label",
                                        suppress_x_axis = FALSE,
                                        plot_data = T) {  # New argument to suppress x-axis
  
  # Helper function to wrap text with a maximum of two words per line, treating hyphenated words as two
  wrap_text_by_words <- function(label, max_words_per_line = 2) {
    words <- strsplit(label, " ")[[1]]  # Split the label by spaces only
    lines <- c()  # To store the lines of text
    current_line <- ""  # To accumulate words into the current line
    current_word_count <- 0  # To keep track of the word count per line
    
    for (word in words) {
      word_count <- ifelse(grepl("-", word), 2, 1)  # Count hyphenated words as two words
      if (current_word_count + word_count > max_words_per_line) {
        lines <- c(lines, current_line)
        current_line <- word
        current_word_count <- word_count
      } else {
        if (current_line == "") {
          current_line <- word
        } else {
          current_line <- paste(current_line, word)
        }
        current_word_count <- current_word_count + word_count
      }
    }
    if (nchar(current_line) > 0) {
      lines <- c(lines, current_line)
    }
    return(paste(lines, collapse = "\n"))
  }
  
  # Function to label specific points with user-defined angle, line length, line color, and text color
  label_points_with_user_inputs <- function(data, x_coords, y_buffer = 0.3, cex = 0.8, 
                                            line_color = "red", text_color = "blue", 
                                            line_width = 1.2, line_type = 1, text_size = 0.7) {
    for (x_name in names(x_coords)) {
      x_val <- x_coords[[x_name]]$x
      angle <- x_coords[[x_name]]$angle
      line_length <- x_coords[[x_name]]$line_length
      
      closest_idx <- which.min(abs(data[[x_col]] - x_val))
      y_val <- data[[y_col_long]][closest_idx]
      
      x_end <- data[[x_col]][closest_idx] + line_length * cos(angle * pi / 180)
      y_end <- y_val + line_length * sin(angle * pi / 180)
      
      # Draw a line from the point to the calculated end point with user-specified color, width, and type
      segments(data[[x_col]][closest_idx], y_val, x_end, y_end, col = line_color, lwd = line_width, lty = line_type)
      
      label_y_end <- ifelse(y_end > y_val, y_end + y_buffer, y_end - y_buffer)
      
      wrapped_label <- wrap_text_by_words(x_name, max_words_per_line = 2)
      
      # Add the text label with user-specified color and size
      text(x_end, label_y_end, labels = wrapped_label, adj = c(0.5, NA), cex = text_size, col = text_color)
    }
  }
  
  # Adjust plot margins if necessary
  par(mar = c(5, 4, 4, 5) + 0.1)
  
  if (suppress_x_axis) {
    plot(NA, NA, 
         xlim = plot_xlim, 
         ylim = plot_ylim, 
         main = title,
         xlab = "",         # Suppress x-axis label
         ylab = ylab,
         xaxt = "n",        
         bty = 'n')
    axis(2)
    axis(1, labels = FALSE, tick = FALSE)
  } else {
    plot(NA, NA,
         xlim = plot_xlim, 
         ylim = plot_ylim, 
         main = title,     
         xlab = xlab,      
         ylab = ylab,      
         bty = 'n')
  }
  
  # Only plot data if requested
  if (plot_data) {
    points(data[[y_col]] ~ data[[x_col]], 
           type = plot_type, 
           cex = plot_cex, 
           col = point_col,  
           pch = point_pch)
    points(data_long[[y_col_long]] ~ data_long[[x_col]], 
           type = plot_type, 
           cex = plot_cex, 
           col = long_point_col,  
           pch = long_point_pch)
  }
  
  # Label the points with user-defined line length, angle, and wrapped text labels
  label_points_with_user_inputs(data_long, x_coords, y_buffer = y_buffer, 
                                cex = label_cex, line_color = line_color, 
                                text_color = text_color, line_width = line_width, 
                                line_type = line_type, text_size = text_size)
}

#' Add Annotated Event Interval to Plot
#'
#' Draws a horizontal event interval with optional shaded background, vertical ticks, and a label.
#' Useful for annotating plots with geologic stages, climate events, or other time-bounded periods.
add_event_interval <- function(start,
                               end,
                               label,
                               y = 0,
                               interval_height = 0.05,
                               label_offset = 0.015,
                               wrap_words = NULL,
                               lwd = 2,
                               col = "black",
                               fill_col = "grey",           # Fill color for shaded bar
                               alpha = 0.5,                # Transparency for shaded bar
                               label_cex = 0.8,
                               label_srt = 0,
                               ...) {
  
  # Wrap label if requested
  if (!is.null(wrap_words)) {
    words <- unlist(strsplit(label, " "))
    wrapped_lines <- split(words, ceiling(seq_along(words) / wrap_words))
    label <- sapply(wrapped_lines, paste, collapse = " ")
    label <- paste(label, collapse = "\n")  # Final wrapped label with newlines
  }
  
  # --- Compute bar boundaries ---
  y_bottom <- y - interval_height / 2
  y_top <- y + interval_height / 2
  
  # --- 1. SHADED RECTANGLE (no border) ---
  usr <- par("usr")        # [x1, x2, y1, y2]
  y_axis_bottom <- usr[3]  # bottom of plot (x-axis)
  
  rect(xleft = start, 
       xright = end, 
       ybottom = y_axis_bottom, 
       ytop = y,                     # <-- Use bar y (not bar handle top)
       col = adjustcolor(fill_col, alpha.f = alpha), 
       border = NA)   # <-- No border
  
  # --- 2. SHORT BAR + LABEL ---
  # Short vertical bars
  segments(x0 = start, y0 = y_bottom, x1 = start, y1 = y_top, col = col, lwd = lwd, ...)
  segments(x0 = end,   y0 = y_bottom, x1 = end,   y1 = y_top, col = col, lwd = lwd, ...)
  
  # Horizontal connector bar
  segments(x0 = start, y0 = y, x1 = end, y1 = y, col = col, lwd = lwd, ...)
  
  # Label above the short bar
  text(x = mean(c(start, end)),
       y = y_top + label_offset,
       labels = label,
       cex = label_cex,
       srt = label_srt)
}


#' Plot IC Acceptance Matrix with Optional Rate of Improvement Overlay
#'
#' Visualizes information criterion (IC) scores across sub-model evaluations, highlighting accepted
#' and rejected shifts. Optionally overlays a secondary y-axis showing the rate of improvement 
#' (first difference of IC scores) with shaded markers and lines.
plot_ic_acceptance_matrix <- function(matrix_data, plot_title = "IC Acceptance Matrix Scatter Plot", 
                                      plot_rate_of_improvement = T) {
  # Adjust margins for balanced spacing
  par(mar = c(5, 5.5, 4, 6), mgp = c(3, 0.6, 0))  # Adjust mgp to move tick labels closer
  
  # Extract y-values and category values
  y_values <- matrix_data[, 1]
  categories <- matrix_data[, 2]
  
  # Calculate rate of improvement (differences between consecutive IC scores)
  rate_of_improvement <- diff(y_values)
  
  # Define x-axis and y-axis limits using pretty ticks with padding
  x_values <- seq_along(y_values)
  x_ticks <- pretty(x_values)
  y_ticks <- pretty(y_values)
  x_limits <- range(x_ticks)
  y_limits <- range(y_ticks)
  
  # Define limits for the rate of improvement (for the secondary y-axis)
  rate_ticks <- pretty(range(rate_of_improvement))
  rate_limits <- c(-400, 150) #range(rate_ticks)
  
  # Identify the baseline IC as the first IC score
  baseline_ic <- y_values[1]
  
  # Plot the rate of improvement optionally
  if (plot_rate_of_improvement) {
    plot(
      x_values[-1], rate_of_improvement,  # Rate of improvement (x values shifted for diff())
      col = NA,  # Suppress default plotting
      type = "n", lty = "solid", lwd = 1,  # Set up the plot environment
      xlab = "", ylab = "",  # Suppress axis labels for overlay
      xlim = x_limits, ylim = rate_limits,  # Secondary y-axis scaling
      xaxt = "n", yaxt = "n", bty = "n"  # Suppress axes for overlay
    )
    
    # Add a black horizontal line at y = 0 with restricted x-range
    lines(
      x = c(min(x_values), max(x_ticks)+20),  # Extend from data limit to axis limit
      y = c(0, 0),  # Horizontal line at y = 0
      col = rgb(0, 0, 0, alpha = 0.5), lty = 1, lwd = 0.7)
    
    # Plot the rate of improvement curve
    lines(
      x_values[-1], rate_of_improvement, 
      col = "grey",  # Semi-transparent black line
      lty = "solid", lwd = 0.8  # Thin solid line
    )
    
    # Add small black dots on the rate curve for accepted shifts
    accepted_x <- which(categories[-1] == 1)  # Accepted shifts correspond to categories == 1
    points(
      x_values[accepted_x + 1], rate_of_improvement[accepted_x],  # Offset by 1 for diff()
      col = rgb(0, 0, 0, alpha = 0.5), pch = 16, cex = 0.3  # Small black dots
    )
    
    # Add the secondary y-axis for rate of improvement with transparency
    axis(
      4, at = rate_ticks, labels = rate_ticks, 
      las = 1, cex.axis = 0.75, tck = -0.02, 
      col = rgb(0, 0, 0, alpha = 0.5),         # Color of ticks matches line transparency
      col.axis = rgb(0, 0, 0, alpha = 0.5)    # Color of tick labels matches line transparency
    )
  }
  
  # Plot the IC scores on top
  par(new = TRUE)  # Enable overlaying
  plot(
    x_values, y_values,
    col = NA,  # Suppress default point plotting; add manually
    xlab = "Sub-model evaluated", 
    ylab = "",  # Leave blank; we'll use mtext for the y-axis label
    type = "n",  # Suppress plotting, just set up the environment
    xlim = x_limits,  # Adjust x-axis to tick-aligned limits
    ylim = y_limits,  # Adjust y-axis to IC score limits
    xaxt = "n",  # Suppress default x-axis
    yaxt = "n",  # Suppress default y-axis
    cex.lab = 1,  # Ensure consistent font size for axis labels
    bty = "n"  # Remove the box around the plot
  )
  
  # Custom title placement
  title(
    main = plot_title,
    line = 2,   # Adjust this value to raise/lower the title
    cex.main = 1.0 # Optional: Adjust title size
  )
  
  # Manually add x-axis
  axis(
    1, at = x_ticks, labels = x_ticks, 
    cex.axis = 1, tck = -0.02  # Ticks and labels for x-axis
  )
  
  # Manually add y-axis for IC scores
  axis(
    2, at = y_ticks, labels = y_ticks, 
    las = 1, cex.axis = 0.75, tck = -0.02  # Horizontal labels and ticks for y-axis
  )
  
  # Add the IC Score y-axis label using mtext for custom positioning
  mtext("IC Score", side = 2, line = 3.5, cex = 0.6)  # Customize 'line' as needed
  
  # Plot the rejected shifts (red dots, smaller size) first (behind accepted shifts)
  points(
    x_values[categories == 0], y_values[categories == 0],
    col = rgb(1, 0, 0, alpha = 1), pch = 3, cex = 0.4, lwd=0.3  # Smaller red dots for rejected shifts
  )
  
  # Combine all blue dots (excluding the baseline IC) for a continuous line
  blue_x <- x_values[categories == 1]
  blue_y <- y_values[categories == 1]
  
  # Plot the blue line connecting all accepted shifts (BEHIND the dots)
  lines(
    blue_x, blue_y,
    col = "blue", type = "l", lwd = 1.1  # Thinner blue line
  )
  
  # Plot the blue dots for accepted shifts with a fine black outline (ON TOP of the line)
  points(
    x_values[-1][categories[-1] == 1], y_values[-1][categories[-1] == 1],
    col = "black", bg = "blue", pch = 21, cex = 0.8, lwd = 0.3  # Blue fill with hairline black outline
  )
  
  # Plot the baseline IC value as a larger black dot (after lines/dots for layering)
  points(
    x = x_values[1], y = baseline_ic, col = "black", pch = 19, cex = 1.0  # Black dot for baseline IC
  )
  
  # Add a label for the baseline IC
  text(
    x = x_values[1] + 2, y = baseline_ic, labels = paste0(round(baseline_ic, 2)),
    pos = 4, col = "black", cex = 0.6
  )
  
  # Add the label for the minimum accepted IC score
  min_accepted_index <- which(categories == 1 & y_values == min(y_values[categories == 1]))
  min_accepted_value <- y_values[min_accepted_index]
  
  text(
    x = min_accepted_index, 
    y = min_accepted_value - diff(range(y_values)) * 0.02,  # Slight vertical offset
    labels = paste0(round(min_accepted_value, 2)), 
    pos = 1, col = "black", cex = 0.6
  )
  
  # Add a clean legend for IC scores, rate of improvement, and baseline IC
  legend(
    "topright",
    inset = c(0.04, -0.10),
    legend = c("Rejected shift", "Accepted shift", "IC Change", "Baseline IC"), 
    col = c("red", "blue", rgb(0, 0, 0, alpha = 0.5), "black"), 
    lty = c(NA, 1, ifelse(plot_rate_of_improvement, 1, NA), NA),  # Line for "Accepted shift" and "IC Change"
    pch = c(3, 21, ifelse(plot_rate_of_improvement, NA, 19), 19), # Cross (3), Dot with line (21), "Baseline IC" (19)
    pt.bg = c(NA, "blue", NA, NA),  # Background for "Accepted shift"
    pt.lwd = c(NA, 0.5, NA, 0),     # Fine outline for "Accepted shift"
    cex = 0.65, bty = "n", xpd = T  # Slightly smaller legend text and allow margin overlap
  )
}


#' Plot Multiple IC Acceptance Curves for Model Comparison
#'
#' Generates overlaid IC (Information Criterion) decline curves for multiple models,
#' highlighting accepted shifts, baseline ICs, and the globally optimal model fit.
plot_multiple_ic_acceptance_curves <- function(model_list, plot_title = "IC Decline Comparison") {
  # Ensure the input is a list
  if (!is.list(model_list)) {
    stop("Input must be a list of model objects.")
  }
  
  # Extract the IC acceptance matrices and baselines
  ic_matrices <- lapply(model_list, function(model) model$model_fit_history$ic_acceptance_matrix)
  baselines <- sapply(ic_matrices, function(mat) mat[1, 1])  # First IC score for each model
  
  # Determine global x and y axis limits
  x_raw_limits <- range(sapply(ic_matrices, function(mat) seq_along(mat[, 1])))
  y_raw_limits <- range(sapply(ic_matrices, function(mat) mat[, 1]))
  
  # Add padding to the y-axis lower limit (10%)
  y_padding <- diff(y_raw_limits) * 0.1
  y_padded_limits <- c(y_raw_limits[1] - y_padding, y_raw_limits[2])
  y_final_limits <- pretty(y_padded_limits)
  
  # Add padding to the x-axis upper limit (5%)
  x_padding <- diff(x_raw_limits) * 0.05
  x_padded_limits <- c(x_raw_limits[1], x_raw_limits[2] + x_padding)
  x_final_limits <- pretty(x_padded_limits)
  
  # Find the globally minimum IC value and its model
  global_min_value <- min(sapply(ic_matrices, function(mat) min(mat[, 1])))
  global_min_model <- which(sapply(ic_matrices, function(mat) any(mat[, 1] == global_min_value)))
  global_min_matrix <- ic_matrices[[global_min_model]]
  global_min_index <- which(global_min_matrix[, 1] == global_min_value)
  
  # Generate unique colors using viridis
  colors <- viridis(length(model_list))
  
  # Adjust margins for shared axes and title
  par(mar = c(5, 5.5, 4, 2), mgp = c(3, 0.6, 0))  # Adjust mgp to move tick labels closer
  
  # Create the base plot
  plot(
    1, type = "n",  # Empty plot to set up the environment
    xlim = range(x_final_limits), ylim = range(y_final_limits), 
    xlab = "Sub-model evaluated", 
    ylab = "IC Score", 
    main = plot_title, 
    xaxt = "n", yaxt = "n", bty = "n"
  )
  
  # Add x and y axes
  axis(1, at = x_final_limits, cex.axis = 1, tck = -0.02)
  axis(2, at = y_final_limits, las = 1, cex.axis = 0.75, tck = -0.02)
  
  # Add IC decline curves for each model
  for (i in seq_along(ic_matrices)) {
    mat <- ic_matrices[[i]]
    y_values <- mat[, 1]
    categories <- mat[, 2]
    x_values <- seq_along(y_values)
    baseline_ic <- baselines[i]
    
    # Plot the line connecting accepted shifts
    lines(
      x_values[categories == 1], y_values[categories == 1],
      col = colors[i], type = "l", lwd = 1.1  # Line for accepted shifts
    )
    
    # Plot the blue dots for accepted shifts (ON TOP of rejected shifts)
    points(
      x_values[categories == 1], y_values[categories == 1],
      col = "black", bg = colors[i], pch = 21, cex = 0.8, lwd = 0.3  # Fine black outline
    )
    
    # Plot the baseline IC value as a larger black dot
    points(
      x = x_values[1], y = baseline_ic, col = "black", pch = 19, cex = 1.0
    )
    
    # Add a label for the baseline IC
    text(
      x = x_values[1] + 2, y = baseline_ic, labels = paste0(round(baseline_ic, 2)),
      pos = 4, col = "black", cex = 0.6
    )
  }
  
  # Add the label for the globally minimum accepted IC score
  global_min_x <- global_min_index
  global_min_y <- global_min_value
  text(
    x = global_min_x, 
    y = global_min_y - diff(range(y_final_limits)) * 0.02,  # Slight vertical offset
    labels = paste0(round(global_min_y, 2)), 
    pos = 1, col = "black", cex = 0.6
  )
  
  # Add a legend to identify each curve
  legend(
    "topright",
    inset = c(0.04, -0.10),
    legend = names(model_list),
    col = colors,
    lty = 1,  # Solid line for each curve
    pch = 21,
    pt.bg = colors,
    pt.lwd = 0.5,
    cex = 0.65, bty = "n", xpd = T
  )
}

#Functions to study the shift counts and the waiting times
{
  
  #' Compute State Shift Rates and Phenotypic Rates from a SIMMAP Tree
  #'
  #' Calculates unweighted and weighted state shift rates, phenotypic rates, and associated metrics 
  #' for each tip in a SIMMAP-formatted phylogenetic tree. Designed to trace lineage state transitions 
  #' and quantify the evolutionary dynamics using optional weighting schemes and state-dependent parameters.
  compute_shift_rates <- function(mapped_tree, weighting = "ES", verbose = TRUE, n_cores = 1, decay_base = 2, param = NULL, use_absolute_magnitudes = F, normalize_weights = F) {
      # Check weighting scheme validity
      if (!is.null(weighting) && !weighting %in% c("ES")) {
        stop("Invalid weighting scheme. Use 'ES' for Equal Splits or NULL for no weighting.")
      }
      
      # Validate param
      if (!is.null(param) && (!is.numeric(param) || is.null(names(param)))) {
        stop("'param' must be a named numeric vector where names correspond to state values.")
      }
      
      # Retrieve states at nodes and tips
      states <- phytools::getStates(mapped_tree, type = "both")
      root_node <- max(mapped_tree$edge[, 1])
      Ntip <- length(mapped_tree$tip.label)
      
      # Extract states for tips
      tip_states <- states[mapped_tree$tip.label]
      
      # Calculate total tree height (maximum distance from root to any tip)
      total_tree_height <- max(ape::node.depth.edgelength(mapped_tree))
      
      # Calculate total times using root-to-tip distances
      root_to_tip_distances <- diag(vcv(mapped_tree))
      
      if (verbose && n_cores == 1) cat("Starting computation for", Ntip, "tips...\n")
      
      # Helper function to log verbose messages
      log_verbose <- function(message, ...) {
        if (verbose && n_cores == 1) {
          cat(sprintf(message, ...), "\n")
        }
      }
      
      # Helper function to detect shifts
      detect_shift <- function(parent_node, child_node, edge_length, cumulative_time, node_age) {
        parent_state <- states[as.character(parent_node)]
        child_state <- states[as.character(child_node)]
        shift_detected <- FALSE
        shift_node_age <- NA
        shift_magnitude <- NA
        
        # Log the parent and child states
        log_verbose("  Parent state: %s, Child state: %s", parent_state, child_state)
        
        if (!is.na(parent_state) && !is.na(child_state) && parent_state != child_state) {
          shift_detected <- TRUE
          shift_node_age <- max(nodeHeights(mapped_tree)) - node_age
          
          # Calculate shift magnitude if param is provided
          if (!is.null(param) && all(c(parent_state, child_state) %in% names(param))) {
            shift_magnitude <- param[parent_state] - param[child_state] #reversed the order bc the child/parent states are mixed up somehow
            if (use_absolute_magnitudes) {
              shift_magnitude <- abs(shift_magnitude)
            }
          }
          
          # Log details about the shift
          log_verbose("    Shift detected: Parent state (%s) -> Child state (%s)", parent_state, child_state)
          log_verbose("    Edge length: %f, Cumulative time: %f", edge_length, cumulative_time)
          log_verbose("    Shift node age (from present): %f", shift_node_age)
          if (!is.na(shift_magnitude)) {
            log_verbose("    Shift magnitude (delta): %f", shift_magnitude)
          }
        }
        
        return(list(shift_detected = shift_detected, shift_node_age = shift_node_age, shift_magnitude = shift_magnitude))
      }
      
      # Helper function to calculate rates
      calculate_rates <- function(unweighted_shift_count, weighted_shift_count, total_time, speciation_events) {
        # Rates per time
        unweighted_shift_rate_per_time <- ifelse(total_time > 0, unweighted_shift_count / total_time, NA)
        weighted_shift_rate_per_time <- if (!is.na(weighted_shift_count) && total_time > 0) weighted_shift_count / total_time else NA
        
        # Rates per speciation event
        unweighted_shift_rate_per_speciation <- ifelse(speciation_events > 0, unweighted_shift_count / speciation_events, NA)
        weighted_shift_rate_per_speciation <- if (!is.na(weighted_shift_count) && speciation_events > 0) weighted_shift_count / speciation_events else NA
        
        # Combined rates (time × speciation events)
        combined_unweighted_rate <- ifelse(total_time > 0 && speciation_events > 0,
                                           unweighted_shift_count / (total_time * speciation_events), NA)
        combined_weighted_rate <- ifelse(!is.na(weighted_shift_count) && total_time > 0 && speciation_events > 0,
                                         weighted_shift_count / (total_time * speciation_events), NA)
        
        return(list(
          unweighted_shift_rate_per_time = unweighted_shift_rate_per_time,
          weighted_shift_rate_per_time = weighted_shift_rate_per_time,
          unweighted_shift_rate_per_speciation = unweighted_shift_rate_per_speciation,
          weighted_shift_rate_per_speciation = weighted_shift_rate_per_speciation,
          combined_unweighted_rate = combined_unweighted_rate,
          combined_weighted_rate = combined_weighted_rate
        ))
      }
      
      # Helper function to calculate Total Weighted Magnitude
      calculate_total_weighted_magnitude <- function(shift_magnitudes, shift_ages, decay_base, normalize_weights) {
        weights <- 1 / (decay_base ^ shift_ages)
        # Normalize weights if required
        if (normalize_weights) {
          weights <- weights / sum(weights, na.rm = TRUE)
        }
        weighted_magnitudes <- shift_magnitudes * weights
        total_weighted_magnitude <- sum(weighted_magnitudes, na.rm = TRUE)
        return(total_weighted_magnitude)
      }
      
      # Helper function to calculate Weighted Phenotype Rates
      calculate_weighted_phenotypic_rates <- function(states, shift_ages, param, decay_base, normalize_weights) {
        # Retrieve parameter values for all states
        param_values <- param[states]
        
        # Calculate weights using decay_base and shift ages
        weights <- 1 / (decay_base ^ shift_ages)
        
        # Normalize weights if required
        if (normalize_weights) {
          weights <- weights / sum(weights, na.rm = TRUE)
        }
        
        # Compute weighted parameter values
        weighted_values <- param_values * weights
        
        # Sum the weighted values to get the total
        total_weighted_rate <- sum(weighted_values, na.rm = TRUE)
        
        return(total_weighted_rate)
      }
      
      # Function to process each tip
      process_tip <- function(k) {
        log_verbose("\nProcessing tip: %s (Tip ID: %d)", mapped_tree$tip.label[k], k)
        
        # Get lineage from tip to root, including root
        lineage_nodes <- c(k, phangorn::Ancestors(mapped_tree, k, type = "all"))
        log_verbose("Lineage nodes (including root): %s", paste(lineage_nodes, collapse = ", "))
        
        shift_positions <- c()
        shift_ages <- c()
        shift_magnitudes <- c()
        cumulative_time <- 0
        total_time <- root_to_tip_distances[k]  # Use precomputed root-to-tip distance for this tip
        speciation_events <- 0  # Count of speciation events along the lineage
        
        for (i in seq_along(lineage_nodes)) {
          node <- lineage_nodes[i]
          
          # Skip the tip node itself when processing (special case handling)
          if (i == 1) next
          
          # Count speciation events (exclude tip and root)
          if (node != k && node != root_node) {
            speciation_events <- speciation_events + 1
          }
          
          if (node != root_node) {
            # Non-root nodes: compare with parent node
            parent_node <- phangorn::Ancestors(mapped_tree, node, type = "parent")
            if (length(parent_node) == 0) next
            parent_node <- parent_node[1]
            
            # Edge index between parent and node
            edge_index <- which((mapped_tree$edge[, 1] == parent_node) & (mapped_tree$edge[, 2] == node))
            if (length(edge_index) == 0) next
            edge_length <- mapped_tree$edge.length[edge_index]
            cumulative_time <- cumulative_time + edge_length
            
            shift_result <- detect_shift(parent_node, node, edge_length, cumulative_time, nodeheight(mapped_tree, node=node))
            if (shift_result$shift_detected) {
              shift_positions <- c(shift_positions, i)
              shift_ages <- c(shift_ages, shift_result$shift_node_age)
              shift_magnitudes <- c(shift_magnitudes, shift_result$shift_magnitude)
            }
          } else {
            # Root node: compare with child node
            if (i < length(lineage_nodes)) {
              child_node <- lineage_nodes[i + 1]
              
              # Edge index between root and child node
              edge_index <- which((mapped_tree$edge[, 1] == node) & (mapped_tree$edge[, 2] == child_node))
              if (length(edge_index) == 0) next
              edge_length <- mapped_tree$edge.length[edge_index]
              cumulative_time <- cumulative_time + edge_length
              
              shift_result <- detect_shift(node, child_node, edge_length, cumulative_time, nodeheight(mapped_tree, node=child_node))
              if (shift_result$shift_detected) {
                shift_positions <- c(shift_positions, i)
                shift_ages <- c(shift_ages, shift_result$shift_node_age)
                shift_magnitudes <- c(shift_magnitudes, shift_result$shift_magnitude)
              }
            }
          }
        }
        
        # Shift count
        unweighted_shift_count <- length(shift_positions)
        
        # Weighted shift count (if applicable)
        if (!is.null(weighting) && weighting == "ES") {
          weights <- 1 / (decay_base^shift_ages)
          # Normalize weights if required
          if (normalize_weights) {
            weights <- weights / sum(weights, na.rm = TRUE)
          }
          log_verbose("  Shift ages: %s", paste(shift_ages, collapse = ", "))
          log_verbose("  Weights for shifts: %s", paste(weights, collapse = ", "))
          weighted_shift_count <- sum(weights)
          inverse_weighted_shift_count <- 1 / weighted_shift_count  # Add inverse weighted count
        } else {
          weighted_shift_count <- NA
          inverse_weighted_shift_count <- NA
        }
        
        # Total Weighted Magnitude
        total_weighted_magnitude <- NA
        if (length(shift_magnitudes) > 0) {
          total_weighted_magnitude <- calculate_total_weighted_magnitude(shift_magnitudes, shift_ages, decay_base, normalize_weights)
        }
        log_verbose("  Total Weighted Magnitude: %f", total_weighted_magnitude)
        
        # **Weighted Phenotypic Rates**
        weighted_phenotypic_rate <- calculate_weighted_phenotypic_rates(
          states = states[lineage_nodes[shift_positions]],  # Updated: only shifted nodes
          shift_ages = shift_ages,
          param = param,
          decay_base = decay_base,
          normalize_weights
        )
        
        log_verbose("  Shift states: %s", paste(states[lineage_nodes[shift_positions]], collapse = ", "))
        log_verbose("  Weighted Phenotypic Rate: %f", weighted_phenotypic_rate)
        
        # Calculate rates
        rates <- calculate_rates(unweighted_shift_count, weighted_shift_count, total_time, speciation_events)
        
        return(c(unweighted_shift_count, weighted_shift_count, inverse_weighted_shift_count, total_time,
                 rates$unweighted_shift_rate_per_time, rates$weighted_shift_rate_per_time,
                 rates$unweighted_shift_rate_per_speciation, rates$weighted_shift_rate_per_speciation,
                 rates$combined_unweighted_rate, rates$combined_weighted_rate, total_weighted_magnitude,
                 weighted_phenotypic_rate))  # Return the new metric here
      }
      
      # Determine whether to use parallel processing
      if (n_cores > 1) {
        results <- pbmclapply(1:Ntip, process_tip, mc.cores = n_cores)
      } else {
        results <- lapply(1:Ntip, process_tip)
      }
      
      # Compile results
      unweighted_shift_counts <- sapply(results, function(x) x[1])
      weighted_shift_counts <- sapply(results, function(x) x[2])
      inverse_weighted_shift_counts <- sapply(results, function(x) x[3])
      total_times <- sapply(results, function(x) x[4])
      unweighted_shift_rate_per_time <- sapply(results, function(x) x[5])
      weighted_shift_rate_per_time <- sapply(results, function(x) x[6])
      unweighted_shift_rate_per_speciation <- sapply(results, function(x) x[7])
      weighted_shift_rate_per_speciation <- sapply(results, function(x) x[8])
      combined_unweighted_rate <- sapply(results, function(x) x[9])
      combined_weighted_rate <- sapply(results, function(x) x[10])
      total_weighted_magnitudes <- sapply(results, function(x) x[11])
      # Extract Weighted Phenotypic Rates
      weighted_phenotypic_rates <- sapply(results, function(x) x[12])  # New metric (12th element)
      
      # Compute Tip_Phenotype_Rate
      tip_phenotype_rate <- if (!is.null(param)) param[tip_states] else rep(NA, Ntip)
      
      # Add to results data frame
      results_df <- data.frame(
        Tip = mapped_tree$tip.label,
        Unweighted_Shift_Count = unweighted_shift_counts,
        Total_Time = total_times,
        Unweighted_Shift_Rate_Per_Time = unweighted_shift_rate_per_time,
        Unweighted_Shift_Rate_Per_Speciation = unweighted_shift_rate_per_speciation,
        #Combined_Unweighted_Rate = combined_unweighted_rate,
        #Total_Weighted_Magnitude = total_weighted_magnitudes,
        Weighted_Phenotypic_Rate = weighted_phenotypic_rates,
        Tip_Phenotype_Rate = tip_phenotype_rate  # Add the new metric here
      )
      
      # Add weighted metrics if applicable
      if (!is.null(weighting) && weighting == "ES") {
        results_df$Weighted_Shift_Count <- weighted_shift_counts
        results_df$Inverse_Weighted_Shift_Count <- inverse_weighted_shift_counts
        results_df$Weighted_Shift_Rate_Per_Time <- weighted_shift_rate_per_time
        results_df$Weighted_Shift_Rate_Per_Speciation <- weighted_shift_rate_per_speciation
        #results_df$Combined_Weighted_Rate <- combined_weighted_rate
      }
      
      # Print metric units explanations only if verbose
      if (verbose) {
        cat("\nUnits of each metric:\n")
        cat("  Tip: The name of the species/tip in the tree.\n")
        cat("  Unweighted_Shift_Count: Number of state shifts along the lineage.\n")
        if (!is.null(weighting) && weighting == "ES") {
          cat("  Weighted_Shift_Count: Weighted number of state shifts along the lineage (using decay parameter weighting).\n")
          cat("    Equation: Weighted_Shift_Count = sum(1 / (decay_base^shift_ages))\n")
          cat("  Inverse_Weighted_Shift_Count: Inverse of the weighted shift count (mirrors DR statistic).\n")
          cat("    Equation: Inverse_Weighted_Shift_Count = 1 / Weighted_Shift_Count\n")
        }
        cat("  Total_Time: Total evolutionary time from root to tip.\n")
        cat("    Equation: Total_Time = sum(edge lengths along lineage)\n")
        cat("  Unweighted_Shift_Rate_Per_Time: Number of shifts per unit time.\n")
        cat("    Equation: Unweighted_Shift_Rate_Per_Time = Unweighted_Shift_Count / Total_Time\n")
        if (!is.null(weighting) && weighting == "ES") {
          cat("  Weighted_Shift_Rate_Per_Time: Weighted number of shifts per unit time (using decay parameter weighting).\n")
          cat("    Equation: Weighted_Shift_Rate_Per_Time = Weighted_Shift_Count / Total_Time\n")
        }
        cat("  Unweighted_Shift_Rate_Per_Speciation: Number of shifts per speciation event.\n")
        cat("    Equation: Unweighted_Shift_Rate_Per_Speciation = Unweighted_Shift_Count / Speciation_Events\n")
        if (!is.null(weighting) && weighting == "ES") {
          cat("  Weighted_Shift_Rate_Per_Speciation: Weighted number of shifts per speciation event (using decay parameter weighting).\n")
          cat("    Equation: Weighted_Shift_Rate_Per_Speciation = Weighted_Shift_Count / Speciation_Events\n")
        }
        #cat("  Combined_Unweighted_Rate: Shifts per unit time per speciation event.\n")
        #cat("    Equation: Combined_Unweighted_Rate = Unweighted_Shift_Count / (Total_Time * Speciation_Events)\n")
        #if (!is.null(weighting) && weighting == "ES") {
        #  cat("  Combined_Weighted_Rate: Weighted shifts per unit time per speciation event (using decay parameter weighting).\n")
        #  cat("    Equation: Combined_Weighted_Rate = Weighted_Shift_Count / (Total_Time * Speciation_Events)\n")
        #}
        #cat("  Total_Weighted_Magnitude: Sum of weighted magnitudes of shifts.\n")
        #cat("    Equation: Total_Weighted_Magnitude = sum(abs(magnitudes) * (1 / (decay_base^shift_ages)))\n")
        cat("  Weighted_Phenotypic_Rate: Weighted rate of phenotypic evolution based on param values and shift ages.\n")
        cat("    Equation: Weighted_Phenotypic_Rate = sum(param[state] * (1 / (decay_base^shift_ages)))\n")
        cat("  Tip_Phenotype_Rate: The param value associated with the tip's final state.\n")
        cat("    Equation: Tip_Phenotype_Rate = param[state_of_tip]\n")
      }
      
      if (verbose && n_cores == 1) {
        cat("\nComputation completed. Summary of results:\n")
        print(head(results_df))
      }
      
      return(results_df)
    }

  #' Compute Shift Rates Based on Discrete State Transitions Along Lineages (Tip-wise)
  #'
  #' This function computes unweighted and weighted shift rates along tip-wise lineages 
  #' of a SIMMAP tree. It calculates rates of state shifts and optionally incorporates 
  #' user-defined phenotype parameters to produce phenotypic evolution metrics. 
  compute_shift_rates.shifts <- function(mapped_tree, weighting = "ES", verbose = TRUE, n_cores = 1, decay_base = 2, param = NULL, use_absolute_magnitudes = F, normalize_weights = F, T = NULL) {
    # Check weighting scheme validity
    if (!is.null(weighting) && !weighting %in% c("ES")) {
      stop("Invalid weighting scheme. Use 'ES' for Equal Splits or NULL for no weighting.")
    }
    
    # Validate param
    if (!is.null(param) && (!is.numeric(param) || is.null(names(param)))) {
      stop("'param' must be a named numeric vector where names correspond to state values.")
    }
    
    # Retrieve states at nodes and tips
    states <- phytools::getStates(mapped_tree, type = "both")
    root_node <- ape::Ntip(mapped_tree)+1
    Ntip <- length(mapped_tree$tip.label)
    
    # Extract states for tips
    tip_states <- states[mapped_tree$tip.label]
    
    # Calculate total times using root-to-tip distances
    root_to_tip_distances <- diag(vcv(mapped_tree))
    
    if (verbose && n_cores == 1) cat("Starting computation for", Ntip, "tips...\n")
    
    # Helper function to log verbose messages
    log_verbose <- function(message, ...) {
      if (verbose && n_cores == 1) {
        cat(sprintf(message, ...), "\n")
      }
    }
    
    # Helper function to detect shifts
    detect_shift <- function(parent_node, child_node, edge_length, cumulative_time, node_age) {
      parent_state <- states[as.character(parent_node)]
      child_state <- states[as.character(child_node)]
      shift_detected <- FALSE
      shift_node_age <- NA
      shift_magnitude <- NA
      
      # Log the parent and child states
      log_verbose("  Parent state: %s, Child state: %s", parent_state, child_state)
      
      if (!is.na(parent_state) && !is.na(child_state) && parent_state != child_state) {
        shift_detected <- TRUE
        shift_node_age <- max(nodeHeights(mapped_tree)) - node_age
        
        # Calculate shift magnitude if param is provided
        if (!is.null(param) && all(c(parent_state, child_state) %in% names(param))) {
          shift_magnitude <- param[parent_state] - param[child_state] #reversed the order bc the child/parent states are mixed up somehow
          if (use_absolute_magnitudes) {
            shift_magnitude <- abs(shift_magnitude)
          }
        }
        
        # Log details about the shift
        log_verbose("    Shift detected: Parent state (%s) -> Child state (%s)", parent_state, child_state)
        log_verbose("    Edge length: %f, Cumulative time: %f", edge_length, cumulative_time)
        log_verbose("    Shift node age (from present): %f", shift_node_age)
        if (!is.na(shift_magnitude)) {
          log_verbose("    Shift magnitude (delta): %f", shift_magnitude)
        }
      }
      
      return(list(shift_detected = shift_detected, shift_node_age = shift_node_age, shift_magnitude = shift_magnitude))
    }
    
    # Helper function to calculate rates
    calculate_rates <- function(unweighted_shift_count, weighted_shift_count, total_time, speciation_events) {
      # Rates per time
      unweighted_shift_rate_per_time <- ifelse(total_time > 0, unweighted_shift_count / total_time, NA)
      weighted_shift_rate_per_time <- if (!is.na(weighted_shift_count) && total_time > 0) weighted_shift_count / total_time else NA
      
      # Rates per speciation event
      unweighted_shift_rate_per_speciation <- ifelse(speciation_events > 0, unweighted_shift_count / speciation_events, NA)
      weighted_shift_rate_per_speciation <- if (!is.na(weighted_shift_count) && speciation_events > 0) weighted_shift_count / speciation_events else NA
      
      # Combined rates (time × speciation events)
      combined_unweighted_rate <- ifelse(total_time > 0 && speciation_events > 0,
                                         unweighted_shift_count / (total_time * speciation_events), NA)
      combined_weighted_rate <- ifelse(!is.na(weighted_shift_count) && total_time > 0 && speciation_events > 0,
                                       weighted_shift_count / (total_time * speciation_events), NA)
      
      return(list(
        unweighted_shift_rate_per_time = unweighted_shift_rate_per_time,
        weighted_shift_rate_per_time = weighted_shift_rate_per_time,
        unweighted_shift_rate_per_speciation = unweighted_shift_rate_per_speciation,
        weighted_shift_rate_per_speciation = weighted_shift_rate_per_speciation,
        combined_unweighted_rate = combined_unweighted_rate,
        combined_weighted_rate = combined_weighted_rate
      ))
    }
    
    calculate_weighted_phenotypic_rates <- function(states, shift_ages, param, decay_base, normalize_weights) {
      # Retrieve parameter values for all states
      param_values <- param[states]
      cat(">> Parameter values for states:\n")
      cat("   States: ", paste(states, collapse = ", "), "\n")
      cat("   Parameter values: ", paste(param_values, collapse = ", "), "\n")
      
      # Calculate weights using decay_base and shift ages
      
      if (!is.null(T)) {
        weights <- 1 / (decay_base ^ (shift_ages / T))  # Half-life decay
      } else {
        weights <- 1 / (decay_base ^ shift_ages)  # Default decay
      }
      
      cat(">> Initial weights (without normalization):\n")
      cat("   Shift ages: ", paste(shift_ages, collapse = ", "), "\n")
      cat("   Weights: ", paste(weights, collapse = ", "), "\n")
      
      # Normalize weights if required
      if (normalize_weights) {
        weights <- weights / sum(weights, na.rm = TRUE)
        cat(">> Normalized weights:\n")
        cat("   Weights: ", paste(weights, collapse = ", "), "\n")
      }
      
      # Compute weighted parameter values
      weighted_values <- param_values * weights
      cat(">> Weighted parameter values:\n")
      cat("   Weighted values: ", paste(weighted_values, collapse = ", "), "\n")
      
      # Sum the weighted values to get the total
      total_weighted_rate <- sum(weighted_values, na.rm = TRUE)
      cat(">> Total weighted phenotypic rate: ", total_weighted_rate, "\n")
      
      return(total_weighted_rate)
    }
    
    # Function to process each tip
    process_tip <- function(k) {
      log_verbose("\nProcessing tip: %s (Tip ID: %d)", mapped_tree$tip.label[k], k)
      
      #Get root age
      root_age <- max(nodeHeights(mapped_tree))
      
      # Get lineage from tip to root, including root
      lineage_nodes <- c(k, phangorn::Ancestors(mapped_tree, k, type = "all"))
      log_verbose("Lineage nodes (including root): %s", paste(lineage_nodes, collapse = ", "))
      
      shift_positions <- c()
      shift_ages <- c()
      shift_states <- c()
      shift_magnitudes <- c()
      cumulative_time <- 0
      total_time <- root_to_tip_distances[k]  # Use precomputed root-to-tip distance for this tip
      speciation_events <- 0  # Count of speciation events along the lineage
      
      for (i in seq_along(lineage_nodes)) {
        node <- lineage_nodes[i]
        
        # Skip the tip node itself when processing (special case handling)
        if (i == 1) next
        
        # Count speciation events (exclude tip and root)
        if (node != k && node != root_node) {
          speciation_events <- speciation_events + 1
        }
        
        if (node != root_node) {
          # Non-root nodes: compare with parent node
          parent_node <- phangorn::Ancestors(mapped_tree, node, type = "parent")
          if (length(parent_node) == 0) next
          parent_node <- parent_node[1]
          
          # Edge index between parent and node
          edge_index <- which((mapped_tree$edge[, 1] == parent_node) & (mapped_tree$edge[, 2] == node))
          if (length(edge_index) == 0) next
          edge_length <- mapped_tree$edge.length[edge_index]
          cumulative_time <- cumulative_time + edge_length
          
          shift_result <- detect_shift(parent_node, node, edge_length, cumulative_time, nodeheight(mapped_tree, node=node))
          if (shift_result$shift_detected) {
            shift_positions <- c(shift_positions, i)
            shift_ages <- c(shift_ages, shift_result$shift_node_age)
            shift_magnitudes <- c(shift_magnitudes, shift_result$shift_magnitude)
            shift_states <- c(shift_states, states[as.character(node)])
          }
        } else {
          # Root node: compare with child node
          if (i < length(lineage_nodes)) {
            child_node <- lineage_nodes[i + 1]
            
            # Edge index between root and child node
            edge_index <- which((mapped_tree$edge[, 1] == node) & (mapped_tree$edge[, 2] == child_node))
            if (length(edge_index) == 0) next
            edge_length <- mapped_tree$edge.length[edge_index]
            cumulative_time <- cumulative_time + edge_length
            
            shift_result <- detect_shift(node, child_node, edge_length, cumulative_time, nodeheight(mapped_tree, node=child_node))
            if (shift_result$shift_detected) {
              shift_positions <- c(shift_positions, i)
              shift_ages <- c(shift_ages, shift_result$shift_node_age)
              shift_magnitudes <- c(shift_magnitudes, shift_result$shift_magnitude)
              shift_states <- c(shift_states, states[as.character(child_node)])
            }
          }
        }
      }
      
      # Shift count
      unweighted_shift_count <- length(shift_positions)
      
      # Weighted shift count (if applicable)
      if (!is.null(weighting) && weighting == "ES") {
        weights <- 1 / (decay_base^shift_ages)
        # Normalize weights if required
        if (normalize_weights) {
          weights <- weights / sum(weights, na.rm = TRUE)
        }
        log_verbose("  Shift ages: %s", paste(shift_ages, collapse = ", "))
        log_verbose("  Weights for shifts: %s", paste(weights, collapse = ", "))
        weighted_shift_count <- sum(weights)
        inverse_weighted_shift_count <- 1 / weighted_shift_count  # Add inverse weighted count
      } else {
        weighted_shift_count <- NA
        inverse_weighted_shift_count <- NA
      }
      
      
      all_shift_ages <- c(shift_ages, root_age)
      all_shift_states <- c(shift_states, states[as.character(root_node)])
      
      # **Weighted Phenotypic Rates**
      weighted_phenotypic_rate <- calculate_weighted_phenotypic_rates(
        states = all_shift_states,
        shift_ages = all_shift_ages,
        param = param,
        decay_base = decay_base,
        normalize_weights
      )
      
      log_verbose("  Combined shift ages (including root): %s", paste(all_shift_ages, collapse = ", "))
      log_verbose("  Combined shift states (including root): %s", paste(all_shift_states, collapse = ", "))
      log_verbose("  Weighted Phenotypic Rate: %e", weighted_phenotypic_rate)
      
      # Calculate rates
      rates <- calculate_rates(unweighted_shift_count, weighted_shift_count, total_time, speciation_events)
      
      return(c(unweighted_shift_count, weighted_shift_count, inverse_weighted_shift_count, total_time,
               rates$unweighted_shift_rate_per_time, rates$weighted_shift_rate_per_time,
               rates$unweighted_shift_rate_per_speciation, rates$weighted_shift_rate_per_speciation,
               rates$combined_unweighted_rate, rates$combined_weighted_rate,
               weighted_phenotypic_rate))
    }
    
    # Determine whether to use parallel processing
    if (n_cores > 1) {
      results <- pbmclapply(1:Ntip, process_tip, mc.cores = n_cores)
    } else {
      results <- lapply(1:Ntip, process_tip)
    }
    
    # Compile results
    unweighted_shift_counts <- sapply(results, function(x) x[1])
    weighted_shift_counts <- sapply(results, function(x) x[2])
    inverse_weighted_shift_counts <- sapply(results, function(x) x[3])
    total_times <- sapply(results, function(x) x[4])
    unweighted_shift_rate_per_time <- sapply(results, function(x) x[5])
    weighted_shift_rate_per_time <- sapply(results, function(x) x[6])
    unweighted_shift_rate_per_speciation <- sapply(results, function(x) x[7])
    weighted_shift_rate_per_speciation <- sapply(results, function(x) x[8])
    combined_unweighted_rate <- sapply(results, function(x) x[9])
    combined_weighted_rate <- sapply(results, function(x) x[10])
    weighted_phenotypic_rates <- sapply(results, function(x) x[11])
    
    # Compute Tip_Phenotype_Rate
    tip_phenotype_rate <- if (!is.null(param)) param[tip_states] else rep(NA, Ntip)
    
    # Add to results data frame
    results_df <- data.frame(
      Tip = mapped_tree$tip.label,
      Unweighted_Shift_Count = unweighted_shift_counts,
      Total_Time = total_times,
      Unweighted_Shift_Rate_Per_Time = unweighted_shift_rate_per_time,
      Unweighted_Shift_Rate_Per_Speciation = unweighted_shift_rate_per_speciation,
      Weighted_Phenotypic_Rate = weighted_phenotypic_rates,
      Tip_Phenotype_Rate = tip_phenotype_rate  # Add the new metric here
    )
    
    # Add weighted metrics if applicable
    if (!is.null(weighting) && weighting == "ES") {
      #results_df$Weighted_Shift_Count <- weighted_shift_counts
      #results_df$Inverse_Weighted_Shift_Count <- inverse_weighted_shift_counts
      results_df$Weighted_Shift_Rate_Per_Time <- weighted_shift_rate_per_time
      results_df$Weighted_Shift_Rate_Per_Speciation <- weighted_shift_rate_per_speciation
    }
    
    # Print metric units explanations only if verbose
    if (verbose) {
      cat("\nUnits of each metric:\n")
      cat("  Tip: The name of the species/tip in the tree.\n")
      cat("  Unweighted_Shift_Count: Number of state shifts along the lineage.\n")
      if (!is.null(weighting) && weighting == "ES") {
        #cat("  Weighted_Shift_Count: Weighted number of state shifts along the lineage (using decay parameter weighting).\n")
        #cat("    Equation: Weighted_Shift_Count = sum(1 / (decay_base^shift_ages))\n")
        #cat("  Inverse_Weighted_Shift_Count: Inverse of the weighted shift count (mirrors DR statistic).\n")
        #cat("    Equation: Inverse_Weighted_Shift_Count = 1 / Weighted_Shift_Count\n")
      }
      cat("  Total_Time: Total evolutionary time from root to tip.\n")
      cat("    Equation: Total_Time = sum(edge lengths along lineage)\n")
      cat("  Unweighted_Shift_Rate_Per_Time: Number of shifts per unit time.\n")
      cat("    Equation: Unweighted_Shift_Rate_Per_Time = Unweighted_Shift_Count / Total_Time\n")
      if (!is.null(weighting) && weighting == "ES") {
        cat("  Weighted_Shift_Rate_Per_Time: Weighted number of shifts per unit time (using decay parameter weighting).\n")
        cat("    Equation: Weighted_Shift_Rate_Per_Time = Weighted_Shift_Count / Total_Time\n")
      }
      cat("  Unweighted_Shift_Rate_Per_Speciation: Number of shifts per speciation event.\n")
      cat("    Equation: Unweighted_Shift_Rate_Per_Speciation = Unweighted_Shift_Count / Speciation_Events\n")
      if (!is.null(weighting) && weighting == "ES") {
        cat("  Weighted_Shift_Rate_Per_Speciation: Weighted number of shifts per speciation event (using decay parameter weighting).\n")
        cat("    Equation: Weighted_Shift_Rate_Per_Speciation = Weighted_Shift_Count / Speciation_Events\n")
      }
      cat("  Weighted_Phenotypic_Rate: Weighted rate of phenotypic evolution based on param values and shift ages.\n")
      cat("    Equation: Weighted_Phenotypic_Rate = sum(param[state] * (1 / (decay_base^shift_ages)))\n")
      cat("  Tip_Phenotype_Rate: The param value associated with the tip's final state.\n")
      cat("    Equation: Tip_Phenotype_Rate = param[state_of_tip]\n")
    }
    
    if (verbose && n_cores == 1) {
      cat("\nComputation completed. Summary of results:\n")
      print(head(results_df))
    }
    
    return(results_df)
  }
  
  #' Compute Branch-Based Shift and Phenotypic Rates for Each Tip
  #'
  #' Calculates unweighted and weighted rates of state shifts and phenotypic evolution along 
  #' entire branches (not just at shift points) leading to each tip of a SIMMAP-formatted tree.
  #' It uses full node-level state histories to weight parameters across the entire ancestral lineage.
  compute_shift_rates.branches <- function(mapped_tree, weighting = "ES", verbose = TRUE, n_cores = 1, decay_base = 2, param = NULL, use_absolute_magnitudes = F, normalize_weights = F, T = NULL) {
    # Check weighting scheme validity
    if (!is.null(weighting) && !weighting %in% c("ES")) {
      stop("Invalid weighting scheme. Use 'ES' for Equal Splits or NULL for no weighting.")
    }
    
    # Validate param
    if (!is.null(param) && (!is.numeric(param) || is.null(names(param)))) {
      stop("'param' must be a named numeric vector where names correspond to state values.")
    }
    
    # Retrieve states at nodes and tips
    states <- phytools::getStates(mapped_tree, type = "both")
    root_node <- ape::Ntip(mapped_tree)+1
    Ntip <- length(mapped_tree$tip.label)
    
    # Extract states for tips
    tip_states <- states[mapped_tree$tip.label]
    
    # Calculate total times using root-to-tip distances
    root_to_tip_distances <- diag(vcv(mapped_tree))
    
    if (verbose && n_cores == 1) cat("Starting computation for", Ntip, "tips...\n")
    
    # Helper function to log verbose messages
    log_verbose <- function(message, ...) {
      if (verbose && n_cores == 1) {
        cat(sprintf(message, ...), "\n")
      }
    }
    
    # Helper function to detect shifts
    detect_shift <- function(parent_node, child_node, edge_length, cumulative_time, node_age) {
      parent_state <- states[as.character(parent_node)]
      child_state <- states[as.character(child_node)]
      shift_detected <- FALSE
      shift_node_age <- NA
      shift_magnitude <- NA
      
      # Log the parent and child states
      log_verbose("  Parent state: %s, Child state: %s", parent_state, child_state)
      
      if (!is.na(parent_state) && !is.na(child_state) && parent_state != child_state) {
        shift_detected <- TRUE
        shift_node_age <- max(nodeHeights(mapped_tree)) - node_age
        
        # Calculate shift magnitude if param is provided
        if (!is.null(param) && all(c(parent_state, child_state) %in% names(param))) {
          shift_magnitude <- param[parent_state] - param[child_state] #reversed the order bc the child/parent states are mixed up somehow
          if (use_absolute_magnitudes) {
            shift_magnitude <- abs(shift_magnitude)
          }
        }
        
        # Log details about the shift
        log_verbose("    Shift detected: Parent state (%s) -> Child state (%s)", parent_state, child_state)
        log_verbose("    Edge length: %f, Cumulative time: %f", edge_length, cumulative_time)
        log_verbose("    Shift node age (from present): %f", shift_node_age)
        if (!is.na(shift_magnitude)) {
          log_verbose("    Shift magnitude (delta): %f", shift_magnitude)
        }
      }
      
      return(list(shift_detected = shift_detected, shift_node_age = shift_node_age, shift_magnitude = shift_magnitude))
    }
    
    # Helper function to calculate rates
    calculate_rates <- function(unweighted_shift_count, weighted_shift_count, total_time, speciation_events) {
      # Rates per time
      unweighted_shift_rate_per_time <- ifelse(total_time > 0, unweighted_shift_count / total_time, NA)
      weighted_shift_rate_per_time <- if (!is.na(weighted_shift_count) && total_time > 0) weighted_shift_count / total_time else NA
      
      # Rates per speciation event
      unweighted_shift_rate_per_speciation <- ifelse(speciation_events > 0, unweighted_shift_count / speciation_events, NA)
      weighted_shift_rate_per_speciation <- if (!is.na(weighted_shift_count) && speciation_events > 0) weighted_shift_count / speciation_events else NA
      
      # Combined rates (time × speciation events)
      combined_unweighted_rate <- ifelse(total_time > 0 && speciation_events > 0,
                                         unweighted_shift_count / (total_time * speciation_events), NA)
      combined_weighted_rate <- ifelse(!is.na(weighted_shift_count) && total_time > 0 && speciation_events > 0,
                                       weighted_shift_count / (total_time * speciation_events), NA)
      
      return(list(
        unweighted_shift_rate_per_time = unweighted_shift_rate_per_time,
        weighted_shift_rate_per_time = weighted_shift_rate_per_time,
        unweighted_shift_rate_per_speciation = unweighted_shift_rate_per_speciation,
        weighted_shift_rate_per_speciation = weighted_shift_rate_per_speciation,
        combined_unweighted_rate = combined_unweighted_rate,
        combined_weighted_rate = combined_weighted_rate
      ))
    }
    
    calculate_weighted_phenotypic_rates <- function(states, shift_ages, param, decay_base, normalize_weights) {
      # Retrieve parameter values for all states
      param_values <- param[states]
      cat(">> Parameter values for states:\n")
      cat("   States: ", paste(states, collapse = ", "), "\n")
      cat("   Parameter values: ", paste(param_values, collapse = ", "), "\n")
      
      # Calculate weights using decay_base and shift ages
      if (!is.null(T)) {
        weights <- 1 / (decay_base ^ (shift_ages / T))  # Half-life decay
      } else {
        weights <- 1 / (decay_base ^ shift_ages)  # Default decay
      }
      
      cat(">> Initial weights (without normalization):\n")
      cat("   Shift ages: ", paste(shift_ages, collapse = ", "), "\n")
      cat("   Weights: ", paste(weights, collapse = ", "), "\n")
      
      # Normalize weights if required
      if (normalize_weights) {
        weights <- weights / sum(weights, na.rm = TRUE)
        cat(">> Normalized weights:\n")
        cat("   Weights: ", paste(weights, collapse = ", "), "\n")
      }
      
      # Compute weighted parameter values
      weighted_values <- param_values * weights
      cat(">> Weighted parameter values:\n")
      cat("   Weighted values: ", paste(weighted_values, collapse = ", "), "\n")
      
      # Sum the weighted values to get the total
      total_weighted_rate <- sum(weighted_values, na.rm = TRUE)
      cat(">> Total weighted phenotypic rate: ", total_weighted_rate, "\n")
      
      return(total_weighted_rate)
    }
    
    # Function to process each tip
    process_tip <- function(k) {
      log_verbose("\nProcessing tip: %s (Tip ID: %d)", mapped_tree$tip.label[k], k)
      
      #Get root age
      root_age <- max(nodeHeights(mapped_tree))
      
      # Get lineage from tip to root, including root
      lineage_nodes <- c(k, phangorn::Ancestors(mapped_tree, k, type = "all"))
      log_verbose("Lineage nodes (including root): %s", paste(lineage_nodes, collapse = ", "))
      
      shift_positions <- c()
      shift_ages <- c()
      shift_states <- c()
      shift_magnitudes <- c()
      cumulative_time <- 0
      total_time <- root_to_tip_distances[k]  # Use precomputed root-to-tip distance for this tip
      speciation_events <- 0  # Count of speciation events along the lineage
      
      for (i in seq_along(lineage_nodes)) {
        node <- lineage_nodes[i]
        
        # Skip the tip node itself when processing (special case handling)
        if (i == 1) next
        
        # Count speciation events (exclude tip and root)
        if (node != k && node != root_node) {
          speciation_events <- speciation_events + 1
        }
        
        if (node != root_node) {
          # Non-root nodes: compare with parent node
          parent_node <- phangorn::Ancestors(mapped_tree, node, type = "parent")
          if (length(parent_node) == 0) next
          parent_node <- parent_node[1]
          
          # Edge index between parent and node
          edge_index <- which((mapped_tree$edge[, 1] == parent_node) & (mapped_tree$edge[, 2] == node))
          if (length(edge_index) == 0) next
          edge_length <- mapped_tree$edge.length[edge_index]
          cumulative_time <- cumulative_time + edge_length
          
          shift_result <- detect_shift(parent_node, node, edge_length, cumulative_time, nodeheight(mapped_tree, node=node))
          if (shift_result$shift_detected) {
            shift_positions <- c(shift_positions, i)
            shift_ages <- c(shift_ages, shift_result$shift_node_age)
            shift_magnitudes <- c(shift_magnitudes, shift_result$shift_magnitude)
            shift_states <- c(shift_states, states[as.character(node)])
          }
        } else {
          # Root node: compare with child node
          if (i < length(lineage_nodes)) {
            child_node <- lineage_nodes[i + 1]
            
            # Edge index between root and child node
            edge_index <- which((mapped_tree$edge[, 1] == node) & (mapped_tree$edge[, 2] == child_node))
            if (length(edge_index) == 0) next
            edge_length <- mapped_tree$edge.length[edge_index]
            cumulative_time <- cumulative_time + edge_length
            
            shift_result <- detect_shift(node, child_node, edge_length, cumulative_time, nodeheight(mapped_tree, node=child_node))
            if (shift_result$shift_detected) {
              shift_positions <- c(shift_positions, i)
              shift_ages <- c(shift_ages, shift_result$shift_node_age)
              shift_magnitudes <- c(shift_magnitudes, shift_result$shift_magnitude)
              shift_states <- c(shift_states, states[as.character(child_node)])
            }
          }
        }
      }
      
      # Shift count
      unweighted_shift_count <- length(shift_positions)
      
      # Weighted shift count (if applicable)
      if (!is.null(weighting) && weighting == "ES") {
        weights <- 1 / (decay_base^shift_ages)
        # Normalize weights if required
        if (normalize_weights) {
          weights <- weights / sum(weights, na.rm = TRUE)
        }
        log_verbose("  Shift ages: %s", paste(shift_ages, collapse = ", "))
        log_verbose("  Weights for shifts: %s", paste(weights, collapse = ", "))
        weighted_shift_count <- sum(weights)
        inverse_weighted_shift_count <- 1 / weighted_shift_count  # Add inverse weighted count
      } else {
        weighted_shift_count <- NA
        inverse_weighted_shift_count <- NA
      }
      
      
      #all_shift_ages <- c(shift_ages, root_age)
      #all_shift_states <- c(shift_states, states[as.character(root_node)])
      path_nodes <- phangorn::Ancestors(mapped_tree, k, type = "all")
      
      all_node_states <- sapply(path_nodes, function(n) {
        states[as.character(n)]
      })
      
      all_node_ages <- sapply(path_nodes, function(n) {
        root_age - nodeheight(mapped_tree, node = n)
      })
      
      # # **Weighted Phenotypic Rates**
      # weighted_phenotypic_rate <- calculate_weighted_phenotypic_rates(
      #   states = all_shift_states,
      #   shift_ages = all_shift_ages,
      #   param = param,
      #   decay_base = decay_base,
      #   normalize_weights
      # )
      
      weighted_phenotypic_rate <- calculate_weighted_phenotypic_rates(
        states = all_node_states,
        shift_ages = all_node_ages,
        param = param,
        decay_base = decay_base,
        normalize_weights = normalize_weights
      )
      
      
      log_verbose("  Combined shift ages (including root): %s", paste(all_shift_ages, collapse = ", "))
      log_verbose("  Combined shift states (including root): %s", paste(all_shift_states, collapse = ", "))
      log_verbose("  Weighted Phenotypic Rate: %e", weighted_phenotypic_rate)
      
      #stop()
      
      # Calculate rates
      rates <- calculate_rates(unweighted_shift_count, weighted_shift_count, total_time, speciation_events)
      
      return(c(unweighted_shift_count, weighted_shift_count, inverse_weighted_shift_count, total_time,
               rates$unweighted_shift_rate_per_time, rates$weighted_shift_rate_per_time,
               rates$unweighted_shift_rate_per_speciation, rates$weighted_shift_rate_per_speciation,
               rates$combined_unweighted_rate, rates$combined_weighted_rate,
               weighted_phenotypic_rate))
    }
    
    # Determine whether to use parallel processing
    if (n_cores > 1) {
      results <- pbmclapply(1:Ntip, process_tip, mc.cores = n_cores)
    } else {
      results <- lapply(1:Ntip, process_tip)
    }
    
    # Compile results
    unweighted_shift_counts <- sapply(results, function(x) x[1])
    weighted_shift_counts <- sapply(results, function(x) x[2])
    inverse_weighted_shift_counts <- sapply(results, function(x) x[3])
    total_times <- sapply(results, function(x) x[4])
    unweighted_shift_rate_per_time <- sapply(results, function(x) x[5])
    weighted_shift_rate_per_time <- sapply(results, function(x) x[6])
    unweighted_shift_rate_per_speciation <- sapply(results, function(x) x[7])
    weighted_shift_rate_per_speciation <- sapply(results, function(x) x[8])
    combined_unweighted_rate <- sapply(results, function(x) x[9])
    combined_weighted_rate <- sapply(results, function(x) x[10])
    weighted_phenotypic_rates <- sapply(results, function(x) x[11])
    
    # Compute Tip_Phenotype_Rate
    tip_phenotype_rate <- if (!is.null(param)) param[tip_states] else rep(NA, Ntip)
    
    # Add to results data frame
    results_df <- data.frame(
      Tip = mapped_tree$tip.label,
      Unweighted_Shift_Count = unweighted_shift_counts,
      Total_Time = total_times,
      Unweighted_Shift_Rate_Per_Time = unweighted_shift_rate_per_time,
      Unweighted_Shift_Rate_Per_Speciation = unweighted_shift_rate_per_speciation,
      Weighted_Phenotypic_Rate = weighted_phenotypic_rates,
      Tip_Phenotype_Rate = tip_phenotype_rate  # Add the new metric here
    )
    
    # Add weighted metrics if applicable
    if (!is.null(weighting) && weighting == "ES") {
      #results_df$Weighted_Shift_Count <- weighted_shift_counts
      #results_df$Inverse_Weighted_Shift_Count <- inverse_weighted_shift_counts
      results_df$Weighted_Shift_Rate_Per_Time <- weighted_shift_rate_per_time
      results_df$Weighted_Shift_Rate_Per_Speciation <- weighted_shift_rate_per_speciation
    }
    
    # Print metric units explanations only if verbose
    if (verbose) {
      cat("\nUnits of each metric:\n")
      cat("  Tip: The name of the species/tip in the tree.\n")
      cat("  Unweighted_Shift_Count: Number of state shifts along the lineage.\n")
      if (!is.null(weighting) && weighting == "ES") {
        #cat("  Weighted_Shift_Count: Weighted number of state shifts along the lineage (using decay parameter weighting).\n")
        #cat("    Equation: Weighted_Shift_Count = sum(1 / (decay_base^shift_ages))\n")
        #cat("  Inverse_Weighted_Shift_Count: Inverse of the weighted shift count (mirrors DR statistic).\n")
        #cat("    Equation: Inverse_Weighted_Shift_Count = 1 / Weighted_Shift_Count\n")
      }
      cat("  Total_Time: Total evolutionary time from root to tip.\n")
      cat("    Equation: Total_Time = sum(edge lengths along lineage)\n")
      cat("  Unweighted_Shift_Rate_Per_Time: Number of shifts per unit time.\n")
      cat("    Equation: Unweighted_Shift_Rate_Per_Time = Unweighted_Shift_Count / Total_Time\n")
      if (!is.null(weighting) && weighting == "ES") {
        cat("  Weighted_Shift_Rate_Per_Time: Weighted number of shifts per unit time (using decay parameter weighting).\n")
        cat("    Equation: Weighted_Shift_Rate_Per_Time = Weighted_Shift_Count / Total_Time\n")
      }
      cat("  Unweighted_Shift_Rate_Per_Speciation: Number of shifts per speciation event.\n")
      cat("    Equation: Unweighted_Shift_Rate_Per_Speciation = Unweighted_Shift_Count / Speciation_Events\n")
      if (!is.null(weighting) && weighting == "ES") {
        cat("  Weighted_Shift_Rate_Per_Speciation: Weighted number of shifts per speciation event (using decay parameter weighting).\n")
        cat("    Equation: Weighted_Shift_Rate_Per_Speciation = Weighted_Shift_Count / Speciation_Events\n")
      }
      cat("  Weighted_Phenotypic_Rate: Weighted rate of phenotypic evolution based on param values and shift ages.\n")
      cat("    Equation: Weighted_Phenotypic_Rate = sum(param[state] * (1 / (decay_base^shift_ages)))\n")
      cat("  Tip_Phenotype_Rate: The param value associated with the tip's final state.\n")
      cat("    Equation: Tip_Phenotype_Rate = param[state_of_tip]\n")
    }
    
    if (verbose && n_cores == 1) {
      cat("\nComputation completed. Summary of results:\n")
      print(head(results_df))
    }
    
    return(results_df)
  }
  
  #' Summarize Global Shift Metrics for Rate-Transition Models
  #'
  #' This function summarizes rate shift dynamics from a tree model object (or a list of such objects),
  #' computing global shift rates, waiting times, and half-lives based on regime transitions.
  #' It supports both individual model objects and named lists of model objects.
  summarize_tree_shift_metrics <- function(model_object, verbose = FALSE) {
    # Helper function to process a single model object
    process_single_model <- function(single_model_object, object_name = NA, total_objects = NA) {
      if (!is.na(object_name)) {
        cat(sprintf("Processing object '%s'...\n", object_name))
      }
      
      # Extract tree and validate
      tree <- single_model_object$tree_no_uncertainty_untransformed
      if (is.null(tree)) stop("Tree not found in model object.")
      
      if (verbose) cat("Tree structure retrieved. Checking tree object...\n")
      
      # Extract transitions and suppress all console outputs, warnings, and messages
      transitions <- suppressWarnings(suppressMessages({
        capture.output(
          result <- extractMaxAgeOfRegimesWithRateChanges(single_model_object),
          file = NULL
        )
        result # Return the actual output of the function
      }))
      
      if (verbose) {
        cat("Extracting transitions and identifying increases/decreases...\n")
        print(head(transitions))
      }
      
      # Validate transitions
      if (is.null(transitions) || !is.data.frame(transitions)) {
        stop("Invalid transitions data extracted from the model object.")
      }
      
      # Total branch length
      total_branch_length <- sum(tree$edge.length)
      
      # Calculate shift counts
      total_shifts <- nrow(transitions)
      increase_shifts <- sum(transitions$rate_change == "increase")
      decrease_shifts <- sum(transitions$rate_change == "decrease")
      
      # Calculate shift rates
      shift_rate_all <- total_shifts / total_branch_length
      shift_rate_increase <- increase_shifts / total_branch_length
      shift_rate_decrease <- decrease_shifts / total_branch_length
      
      # Mean waiting times
      mean_waiting_time_all <- total_branch_length / total_shifts
      mean_waiting_time_increase <- if (increase_shifts > 0) total_branch_length / increase_shifts else NA
      mean_waiting_time_decrease <- if (decrease_shifts > 0) total_branch_length / decrease_shifts else NA
      
      # Half-lives
      half_life_all <- log(2) / shift_rate_all
      half_life_increase <- if (increase_shifts > 0) log(2) / shift_rate_increase else NA
      half_life_decrease <- if (decrease_shifts > 0) log(2) / shift_rate_decrease else NA
      
      # Format values for output
      results_table <- data.frame(
        Metric = c(
          "Total Branch Length",
          "Total Shifts",
          "Shift Rate (All)",
          "Shift Rate (Increase)",
          "Shift Rate (Decrease)",
          "Mean Waiting Time (All)",
          "Mean Waiting Time (Increase)",
          "Mean Waiting Time (Decrease)",
          "Half-Life (All)",
          "Half-Life (Increase)",
          "Half-Life (Decrease)"
        ),
        Value = c(
          total_branch_length,
          total_shifts,
          shift_rate_all,
          shift_rate_increase,
          shift_rate_decrease,
          mean_waiting_time_all,
          mean_waiting_time_increase,
          mean_waiting_time_decrease,
          half_life_all,
          half_life_increase,
          half_life_decrease
        ),
        Units = c(
          "Ma",
          "Count",
          "Shifts per Ma",
          "Shifts per Ma",
          "Shifts per Ma",
          "Ma",
          "Ma",
          "Ma",
          "Ma",
          "Ma",
          "Ma"
        )
      )
      
      # Print human-readable output for single-object processing
      if (verbose) {
        cat("\n--- Results Summary for", object_name, "---\n\n")
        print(transform(results_table, Value = formatC(Value, format = "f", digits = 5)), row.names = FALSE)
      }
      
      # Return results
      return(list(
        total_branch_length_Ma = total_branch_length,
        total_shifts = total_shifts,
        shift_rates_per_Ma = list(all = shift_rate_all, increase = shift_rate_increase, decrease = shift_rate_decrease),
        mean_waiting_times_Ma = list(all = mean_waiting_time_all, increase = mean_waiting_time_increase, decrease = mean_waiting_time_decrease),
        half_lives_Ma = list(all = half_life_all, increase = half_life_increase, decrease = half_life_decrease),
        results_table = results_table
      ))
    }
    
    # Main function logic
    if (!is.null(model_object$tree_no_uncertainty_untransformed)) {
      # Input is a single object, process directly
      return(tryCatch(
        process_single_model(model_object, object_name = "Single Model"),
        error = function(e) stop(paste("Error processing the model object:", e$message))
      ))
    } else if (is.list(model_object) && length(model_object) > 0) {
      # Input is a list, process each object
      object_names <- names(model_object)
      if (is.null(object_names)) stop("The input list must have named objects for proper organization.")
      
      results <- setNames(lapply(object_names, function(name) {
        tryCatch(
          process_single_model(model_object[[name]], object_name = name),
          error = function(e) {
            warning(sprintf("Error processing object '%s': %s", name, e$message))
            return(NULL)
          }
        )
      }), object_names)
      
      # Filter out NULL results
      valid_results <- results[!sapply(results, is.null)]
      
      # Generate summary table if results are valid
      if (length(valid_results) > 0) {
        all_metrics <- lapply(valid_results, function(res) res$results_table)
        metric_names <- all_metrics[[1]]$Metric
        metric_units <- all_metrics[[1]]$Units
        values_matrix <- do.call(cbind, lapply(all_metrics, function(df) df$Value))
        
        mean_values <- rowMeans(values_matrix, na.rm = TRUE)
        median_values <- apply(values_matrix, 1, median, na.rm = TRUE)
        sd_values <- apply(values_matrix, 1, sd, na.rm = TRUE)
        
        summary_table <- data.frame(
          Metric = metric_names,
          Units = metric_units,
          Mean = mean_values,
          Median = median_values,
          SD = sd_values
        )
        
        # Print the summary
        cat("\n--- Summary Table (Aggregated Metrics) ---\n\n")
        print(summary_table, row.names = FALSE)
        
        return(list(
          individual_results = valid_results,
          summary_table = summary_table
        ))
      } else {
        warning("No valid results to summarize.")
        return(NULL)
      }
    } else {
      stop("Invalid input: Expected a single model object or a named list of model objects.")
    }
  }
  
  #' Summarize Waiting Times Between Evolutionary State Shifts
  #'
  #' Traverses all root-to-tip paths in a simmap-formatted tree to calculate the distribution of
  #' waiting times (branch lengths) between evolutionary state transitions. Returns summary
  #' statistics including the mean, median, and variance of waiting times.
  summarize_waiting_times <- function(mapped_tree, verbose = TRUE) {
    library(phangorn)
    library(phytools)
    
    # Retrieve states for nodes
    states <- phytools::getStates(mapped_tree, type = "both")
    root_node <- max(mapped_tree$edge[, 1])
    Ntip <- length(mapped_tree$tip.label)
    waiting_times <- numeric()
    
    # Process each tip to extract all waiting times
    for (tip_idx in seq_len(Ntip)) {
      # Get the lineage from tip to root
      lineage_nodes <- c(tip_idx, Ancestors(mapped_tree, tip_idx, type = "all"))
      
      # Traverse the lineage and record waiting times
      for (i in seq_along(lineage_nodes)[-1]) {  # Skip the tip itself
        parent_node <- lineage_nodes[i]
        child_node <- lineage_nodes[i - 1]
        
        # Compare states between parent and child
        parent_state <- states[as.character(parent_node)]
        child_state <- states[as.character(child_node)]
        if (is.na(parent_state) || is.na(child_state) || parent_state == child_state) {
          next
        }
        
        # Record the branch length as a waiting time
        edge_index <- which((mapped_tree$edge[, 1] == parent_node) & (mapped_tree$edge[, 2] == child_node))
        if (length(edge_index) == 0) next
        branch_length <- mapped_tree$edge.length[edge_index]
        
        if (!is.na(branch_length) && branch_length > 0) {
          waiting_times <- c(waiting_times, branch_length)
        }
      }
    }
    
    # If no waiting times are detected, return NA
    if (length(waiting_times) == 0) {
      if (verbose) cat("No shifts detected in the tree. Returning NA for all metrics.\n")
      return(list(
        mean_waiting_time = NA,
        median_waiting_time = NA,
        variance_waiting_time = NA,
        all_waiting_times = numeric()
      ))
    }
    
    # Summarize waiting times
    summary_results <- list(
      mean_waiting_time = mean(waiting_times, na.rm = TRUE),
      median_waiting_time = median(waiting_times, na.rm = TRUE),
      variance_waiting_time = var(waiting_times, na.rm = TRUE),
      all_waiting_times = waiting_times
    )
    
    # Verbose output
    if (verbose) {
      cat("Summary of Waiting Times Across All Root-to-Tip Paths:\n")
      cat("Mean Waiting Time:", summary_results$mean_waiting_time, "\n")
      cat("Median Waiting Time:", summary_results$median_waiting_time, "\n")
      cat("Variance of Waiting Times:", summary_results$variance_waiting_time, "\n")
    }
    
    return(summary_results)
  }
  
  #' Analyze Waiting Times Between Successive Evolutionary Shifts
  #'
  #' Computes and categorizes the waiting times between consecutive rate shifts (e.g., increases or decreases)
  #' in a phylogenetic model. Designed for `bayou`-like model objects that contain a transition table
  #' and an associated tree. The function handles both expected and alternative column formats.
  analyze_shift_waiting_times <- function(model_object, verbose = FALSE) {
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package 'ape' is required but not installed.")
    }
    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop("Package 'phytools' is required but not installed.")
    }
    
    # Step 1: Extract Transitions
    if (verbose) cat("Extracting transitions from model object...\n")
    transitions <- suppressWarnings(suppressMessages({
      capture.output(
        result <- extractMaxAgeOfRegimesWithRateChanges(model_object),
        file = NULL
      )
      result
    }))
    
    if (!is.data.frame(transitions)) stop("Transitions extraction failed. Not a valid data frame.")
    cat(sprintf("Transitions table has %d rows and %d columns.\n", nrow(transitions), ncol(transitions)))
    
    # Step 2: Detect and Map Column Names
    expected_cols <- c("node", "type")
    alternative_cols <- c("max_age_node", "rate_change")
    
    if (all(expected_cols %in% colnames(transitions))) {
      if (verbose) cat("Found expected column names: node, type.\n")
    } else if (all(alternative_cols %in% colnames(transitions))) {
      colnames(transitions)[colnames(transitions) == "max_age_node"] <- "node"
      colnames(transitions)[colnames(transitions) == "rate_change"] <- "type"
      if (verbose) cat("Renamed columns: max_age_node -> node, rate_change -> type.\n")
    } else {
      cat("Available columns in transitions:\n")
      print(colnames(transitions))
      stop("Transitions data frame must contain one of the following column sets:\n",
           "1. node, type\n",
           "2. max_age_node, rate_change\n")
    }
    
    # Step 3: Extract Tree and Compute Node Heights
    if (verbose) cat("Extracting tree and computing node heights...\n")
    tree <- model_object$tree_no_uncertainty_untransformed
    if (is.null(tree)) stop("Tree not found in model object.")
    
    node_heights_matrix <- phytools::nodeHeights(tree)
    
    # Compute node heights for all nodes
    all_nodes <- 1:(ape::Ntip(tree) + tree$Nnode)
    node_height_vector <- numeric(length = length(all_nodes))
    for (i in 1:nrow(node_heights_matrix)) {
      child <- tree$edge[i, 2]
      node_height_vector[child] <- node_heights_matrix[i, 2]
    }
    
    # Create a data frame for node heights
    node_heights_df <- data.frame(
      node = all_nodes,
      height = node_height_vector,
      stringsAsFactors = FALSE
    )
    
    # Step 4: Merge Shifts with Node Heights
    if (verbose) cat("Merging shifts with node heights...\n")
    shifts_with_time <- merge(transitions, node_heights_df, by = "node")
    
    if (any(is.na(shifts_with_time$height))) {
      warning("Some shifts have missing node heights. These shifts will be excluded.")
      shifts_with_time <- shifts_with_time[!is.na(shifts_with_time$height), ]
    }
    cat(sprintf("Total shifts after merging with node heights: %d\n", nrow(shifts_with_time)))
    
    # Step 5: Sort Shifts Chronologically
    if (verbose) cat("Sorting shifts chronologically...\n")
    shifts_sorted <- shifts_with_time[order(shifts_with_time$height), ]
    
    # Add a synthetic shift for the root node
    root_node <- max(tree$edge[, 1])
    if (!any(shifts_sorted$node == root_node)) {
      root_height <- 0  # Root always at height 0
      root_state <- shifts_sorted$type[1]  # First observed state
      root_shift <- data.frame(
        node = root_node,
        type = "root",
        height = root_height,
        stringsAsFactors = FALSE
      )
      
      # Ensure the synthetic root row matches the column structure of shifts_sorted
      missing_columns <- setdiff(names(shifts_sorted), names(root_shift))
      for (col in missing_columns) {
        root_shift[[col]] <- NA
      }
      
      shifts_sorted <- rbind(root_shift, shifts_sorted)
    }
    
    # Step 6: Calculate Waiting Times Between Consecutive Shifts
    if (verbose) cat("Calculating waiting times between consecutive shifts...\n")
    if (nrow(shifts_sorted) < 2) {
      stop("Not enough shifts to compute waiting times.")
    }
    
    waiting_times <- data.frame(
      PreviousShift = shifts_sorted$type[-nrow(shifts_sorted)],
      NextShift = shifts_sorted$type[-1],
      TimeToNext = diff(shifts_sorted$height),
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("Calculated %d waiting times.\n", nrow(waiting_times)))
    
    # Step 7: Categorize Shift Sequences with All Combinations
    if (verbose) cat("Categorizing shift sequences...\n")
    waiting_times$Category <- ifelse(
      waiting_times$PreviousShift == "increase" & waiting_times$NextShift == "increase", "Increase → Increase",
      ifelse(
        waiting_times$PreviousShift == "increase" & waiting_times$NextShift == "decrease", "Increase → Decrease",
        ifelse(
          waiting_times$PreviousShift == "decrease" & waiting_times$NextShift == "increase", "Decrease → Increase",
          ifelse(
            waiting_times$PreviousShift == "decrease" & waiting_times$NextShift == "decrease", "Decrease → Decrease",
            "Other"
          )
        )
      )
    )
    
    # Step 8: Summarize Waiting Times Including All Categories
    if (verbose) cat("Summarizing waiting times...\n")
    categories <- unique(waiting_times$Category)
    summary_stats <- data.frame(
      Category = character(),
      Mean = numeric(),
      Variance = numeric(),
      Count = integer(),
      stringsAsFactors = FALSE
    )
    
    for (cat in categories) {
      times <- waiting_times$TimeToNext[waiting_times$Category == cat]
      summary_stats <- rbind(
        summary_stats,
        data.frame(
          Category = cat,
          Mean = mean(times, na.rm = TRUE),
          Variance = var(times, na.rm = TRUE),
          Count = length(times),
          stringsAsFactors = FALSE
        )
      )
    }
    
    # Display summary
    cat("\n--- Waiting Time Analysis ---\n\n")
    print(summary_stats, row.names = FALSE)
    
    # Return results
    return(list(
      waiting_times = waiting_times,
      stats = summary_stats
    ))
  }
  
  #' Analyze Waiting Times Between Shifts Along Individual Lineages
  #'
  #' Computes waiting times between consecutive rate shifts (e.g., increases or decreases)
  #' along individual root-to-tip paths (lineages) in a phylogenetic tree. This branch-level
  #' analysis offers a finer resolution compared to global waiting time summaries.
  analyze_shift_waiting_times_by_branch <- function(model_object, verbose = FALSE) {
    # ---------------------------------------------------------
    # Function to Analyze Waiting Times Between Shifts
    # within Each Lineage (Branch-by-Branch) in a Phylogenetic Tree
    #
    # Parameters:
    # - model_object: The phylogenetic model object containing the tree and shift data.
    # - verbose: Logical flag to print detailed processing information.
    #
    # Returns:
    # A list containing:
    # - lineage_waiting_times: A list of data frames with waiting times per lineage.
    # - all_waiting_times: A combined data frame of all waiting times across lineages.
    # - stats: Summary statistics for each shift category.
    # ---------------------------------------------------------
    
    # -------------------------------
    # Step 0: Load Required Packages
    # -------------------------------
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package 'ape' is required but not installed.")
    }
    if (!requireNamespace("phytools", quietly = TRUE)) {
      stop("Package 'phytools' is required but not installed.")
    }
    
    # -------------------------------
    # Helper Function: Identify Root Node
    # -------------------------------
    get_root_node <- function(tree) {
      # Identifies the root node as the node that is a parent but not a child
      all_parent_nodes <- unique(tree$edge[,1])
      all_child_nodes <- unique(tree$edge[,2])
      root_node <- setdiff(all_parent_nodes, all_child_nodes)
      if(length(root_node) !=1){
        stop("Cannot uniquely identify the root node.")
      }
      return(root_node)
    }
    
    # -------------------------------
    # Step 1: Extract Transitions
    # -------------------------------
    if (verbose) cat("Extracting transitions from model object...\n")
    transitions <- suppressWarnings(suppressMessages({
      capture.output(
        result <- extractMaxAgeOfRegimesWithRateChanges(model_object),
        file = NULL
      )
      result
    }))
    
    if (!is.data.frame(transitions)) stop("Transitions extraction failed. Not a valid data frame.")
    cat(sprintf("Transitions table has %d rows and %d columns.\n", nrow(transitions), ncol(transitions)))
    
    # -------------------------------
    # Step 2: Detect and Map Column Names
    # -------------------------------
    # Define possible column name sets
    expected_cols <- c("node", "type")
    alternative_cols <- c("max_age_node", "rate_change")
    
    if (all(expected_cols %in% colnames(transitions))) {
      if (verbose) cat("Found expected column names: node, type.\n")
    } else if (all(alternative_cols %in% colnames(transitions))) {
      # Rename columns to standard names
      colnames(transitions)[colnames(transitions) == "max_age_node"] <- "node"
      colnames(transitions)[colnames(transitions) == "rate_change"] <- "type"
      if (verbose) cat("Renamed columns: max_age_node -> node, rate_change -> type.\n")
    } else {
      # Print available columns for debugging
      cat("Available columns in transitions:\n")
      print(colnames(transitions))
      stop("Transitions data frame must contain one of the following column sets:\n",
           "1. node, type\n",
           "2. max_age_node, rate_change\n")
    }
    
    # -------------------------------
    # Step 3: Extract Tree and Compute Node Heights
    # -------------------------------
    if (verbose) cat("Extracting tree and computing node heights...\n")
    tree <- model_object$tree_no_uncertainty_untransformed
    if (is.null(tree)) stop("Tree not found in model object.")
    
    # Compute node heights using phytools::nodeHeights
    node_heights_matrix <- phytools::nodeHeights(tree)
    
    # Compute node heights for all nodes
    all_nodes <- 1:(ape::Ntip(tree) + tree$Nnode)
    node_height_vector <- numeric(length = length(all_nodes))
    for (i in 1:nrow(node_heights_matrix)) {
      child <- tree$edge[i, 2]
      node_height_vector[child] <- node_heights_matrix[i, 2]
    }
    
    # Create a data frame mapping nodes to their heights
    node_heights_df <- data.frame(
      node = all_nodes,
      height = node_height_vector,
      stringsAsFactors = FALSE
    )
    
    # -------------------------------
    # Step 4: Merge Shifts with Node Heights
    # -------------------------------
    if (verbose) cat("Merging shifts with node heights...\n")
    shifts_with_time <- merge(transitions, node_heights_df, by = "node")
    
    # Exclude shifts originating from the root (type "root")
    shifts_with_time <- shifts_with_time[shifts_with_time$type != "root", ]
    
    # Exclude shifts with missing heights
    if (any(is.na(shifts_with_time$height))) {
      warning("Some shifts have missing node heights. These shifts will be excluded.")
      shifts_with_time <- shifts_with_time[!is.na(shifts_with_time$height), ]
    }
    cat(sprintf("Total shifts after merging with node heights: %d\n", nrow(shifts_with_time)))
    
    # -------------------------------
    # Step 5: Extract Root-to-Tip Paths
    # -------------------------------
    if (verbose) cat("Extracting root-to-tip paths...\n")
    root_node <- get_root_node(tree)
    paths <- lapply(1:ape::Ntip(tree), function(tip) ape::nodepath(tree, from = root_node, to = tip))
    
    # -------------------------------
    # Step 6: Calculate Waiting Times Within Each Lineage
    # -------------------------------
    if (verbose) cat("Calculating waiting times for each lineage...\n")
    lineage_waiting_times <- list()
    
    for (i in seq_along(paths)) {
      path <- paths[[i]]
      # Select shifts that are within the current path
      path_shifts <- shifts_with_time[shifts_with_time$node %in% path, ]
      
      # Proceed only if there are at least two shifts in the path
      if (nrow(path_shifts) >= 2) {
        # Sort shifts by height to ensure chronological order
        path_shifts <- path_shifts[order(path_shifts$height), ]
        
        # Calculate waiting times between consecutive shifts
        waiting_times <- data.frame(
          Lineage = i,
          PreviousShift = path_shifts$type[-nrow(path_shifts)],
          NextShift = path_shifts$type[-1],
          TimeToNext = diff(path_shifts$height),
          stringsAsFactors = FALSE
        )
        
        # Append to the list
        lineage_waiting_times[[i]] <- waiting_times
      } else if (verbose) {
        cat(sprintf("Lineage %d has less than 2 shifts and will be excluded.\n", i))
      }
    }
    
    # -------------------------------
    # Step 7: Combine Waiting Times Across Lineages
    # -------------------------------
    if (length(lineage_waiting_times) > 0) {
      all_waiting_times <- do.call(rbind, lineage_waiting_times)
    } else {
      stop("No waiting times were calculated. Check if your data contains sufficient shifts.")
    }
    cat(sprintf("Calculated waiting times for %d lineages.\n", length(lineage_waiting_times)))
    
    # -------------------------------
    # Step 8: Categorize Shift Sequences
    # -------------------------------
    if (verbose) cat("Categorizing shift sequences...\n")
    all_waiting_times$Category <- ifelse(
      all_waiting_times$PreviousShift == "increase" & all_waiting_times$NextShift == "increase", "Increase → Increase",
      ifelse(
        all_waiting_times$PreviousShift == "increase" & all_waiting_times$NextShift == "decrease", "Increase → Decrease",
        ifelse(
          all_waiting_times$PreviousShift == "decrease" & all_waiting_times$NextShift == "increase", "Decrease → Increase",
          ifelse(
            all_waiting_times$PreviousShift == "decrease" & all_waiting_times$NextShift == "decrease", "Decrease → Decrease",
            "Other"  # Captures any unexpected shift sequences
          )
        )
      )
    )
    
    # -------------------------------
    # Step 9: Summarize Waiting Times
    # -------------------------------
    if (verbose) cat("Summarizing waiting times...\n")
    categories <- unique(all_waiting_times$Category)
    summary_stats <- data.frame(
      Category = character(),
      Mean = numeric(),
      Variance = numeric(),
      Count = integer(),
      stringsAsFactors = FALSE
    )
    
    for (cat in categories) {
      times <- all_waiting_times$TimeToNext[all_waiting_times$Category == cat]
      summary_stats <- rbind(
        summary_stats,
        data.frame(
          Category = cat,
          Mean = mean(times, na.rm = TRUE),
          Variance = var(times, na.rm = TRUE),
          Count = length(times),
          stringsAsFactors = FALSE
        )
      )
    }
    
    # -------------------------------
    # Step 10: Display Summary
    # -------------------------------
    cat("\n--- Waiting Time Analysis (Branch-by-Branch) ---\n\n")
    print(summary_stats, row.names = FALSE)
    
    # -------------------------------
    # Step 11: Return Results
    # -------------------------------
    return(list(
      lineage_waiting_times = lineage_waiting_times,  # List of waiting times per lineage
      all_waiting_times = all_waiting_times,          # Combined waiting times across all lineages
      stats = summary_stats                           # Summary statistics
    ))
  }
  
  
}

#' Create Paired Raincloud Plot with Wilcoxon Test Annotation
#'
#' Generates a paired raincloud plot comparing two related distributions (e.g., "decrease" vs "increase"),
#' with statistical annotation from a Wilcoxon signed-rank test. The plot combines half-violin, boxplot,
#' and dot plot elements to illustrate both distribution and individual paired observations.
create_raincloud_with_wilcox <- function(data,
                                         alternative = "two.sided",
                                         title = "Shift Frequencies",
                                         palette = "Dark2",
                                         horizontal_spacing = 0.1,
                                         y_limits = NULL,
                                         group_colors = NULL) {
  # Load necessary libraries
  library(ggplot2)
  library(ggrain)    # Ensure this package is installed
  library(ggsignif)
  library(tidyr)
  
  # ==============================
  # 1. Data Validation and Preparation
  # ==============================
  
  data <- as.data.frame(data)
  
  required_cols <- c("decrease", "increase")
  if (!all(required_cols %in% names(data))) {
    stop(paste(
      "Data must contain the following columns:",
      paste(required_cols, collapse = ", ")
    ))
  }
  
  if (any(is.na(data$decrease)) || any(is.na(data$increase))) {
    stop("Data contains missing values. Please handle them before plotting.")
  }
  
  data$id <- seq_len(nrow(data))
  
  data_long <- pivot_longer(
    data,
    cols = c("decrease", "increase"),
    names_to = "Group",
    values_to = "Frequency"
  )
  
  data_long$Group <- factor(data_long$Group, levels = c("decrease", "increase"))
  
  # ==============================
  # 2. Statistical Testing
  # ==============================
  
  wilcox_test <- wilcox.test(
    x = data$decrease,
    y = data$increase,
    paired = TRUE,
    exact = FALSE,
    alternative = alternative
  )
  
  test_statistic <- signif(wilcox_test$statistic, 3)
  p_value <- signif(wilcox_test$p.value, 3)
  annotation_text <- paste0("W = ", test_statistic, ", p = ", p_value)
  
  # ==============================
  # 3. Dynamic y_position Calculation
  # ==============================
  
  max_freq <- max(data_long$Frequency, na.rm = TRUE)
  y_pos <- max_freq + 0.15 * max_freq
  
  # ==============================
  # 4. Plot Generation
  # ==============================
  
  plot <- ggplot(data_long,
                 aes(
                   x = Group,
                   y = Frequency,
                   fill = Group,
                   color = Group
                 )) +
    geom_rain(
      rain.side = 'f1x1',
      id.long.var = "id",
      line.args = list(
        color = "black",
        alpha = 0.5,
        linewidth = 0.25
      ),
      boxplot.args = list(color = "black", outlier.shape = NA),
      violin.args = list(alpha = 0.4, color = NA),
      point.args = list(size = 3, alpha = 0.8)
    ) +
    geom_signif(
      comparisons = list(c("decrease", "increase")),
      annotations = annotation_text,
      y_position = y_pos * 0.9,
      tip_length = 0.01,
      textsize = 4,
      color = 'black'
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5),
      plot.background = element_blank()
    ) +
    guides(fill = 'none', color = 'none') +
    labs(title = title, x = "Group", y = "Shift Frequency") +
    scale_y_continuous(limits = y_limits) +
    scale_x_discrete(expand = expansion(mult = c(horizontal_spacing, horizontal_spacing)))
  
  # ==============================
  # 5. Apply Color Scales
  # ==============================
  
  if (!is.null(group_colors)) {
    plot <- plot +
      scale_fill_manual(values = group_colors) +
      scale_color_manual(values = group_colors)
  } else {
    plot <- plot +
      scale_fill_brewer(palette = palette) +
      scale_color_brewer(palette = palette)
  }
  
  # ==============================
  # 6. Return
  # ==============================
  
  print(wilcox_test)
  return(plot)
}

#' Analyze Paired Count Data with Statistical Tests
#'
#' Performs a suite of statistical analyses on paired count data representing two conditions 
#' (e.g., 'decrease' and 'increase'). Includes Wilcoxon signed-rank test, paired t-tests on 
#' raw counts and proportions, and binomial proportion tests both per replicate and in total.
analyze_paired_counts <- function(data) {
  # Ensure the data has the required columns
  if (!all(c("decrease", "increase") %in% names(data))) {
    stop("The data frame must have 'decrease' and 'increase' columns.")
  }
  
  # Wilcoxon Signed-Rank Test
  wilcox_result <- wilcox.test(data$decrease, data$increase, paired = TRUE, exact = FALSE)
  cat("Wilcoxon Signed-Rank Test:\n")
  print(wilcox_result)
  
  # Paired t-Test on Counts
  t_test_result <- t.test(data$decrease, data$increase, paired = TRUE)
  cat("\nPaired t-Test on Counts:\n")
  print(t_test_result)
  
  # Compute proportions
  prop_decrease <- data$decrease / (data$decrease + data$increase)
  prop_increase <- data$increase / (data$decrease + data$increase)
  
  # Paired t-Test on Proportions
  prop_t_test <- t.test(prop_decrease, prop_increase, paired = TRUE)
  cat("\nPaired t-Test on Proportions:\n")
  print(prop_t_test)
  
  # Binomial Tests for Each Replicate
  cat("\nBinomial Tests for Each Replicate:\n")
  for (i in 1:nrow(data)) {
    cat(sprintf("Replicate: %s\n", rownames(data)[i]))
    
    # Create a matrix for prop.test
    test_matrix <- matrix(c(data$decrease[i], data$increase[i]), nrow = 1, byrow = TRUE)
    
    # Perform a binomial proportion test using matrix
    binom_test <- prop.test(test_matrix, correct = FALSE)
    print(binom_test)
  }
  
  # Overall Binomial Test
  cat("\nSummary of Total Counts:\n")
  total_counts_matrix <- matrix(c(sum(data$decrease), sum(data$increase)), nrow = 1)
  overall_test <- prop.test(total_counts_matrix, correct = FALSE)
  print(overall_test)
}


#Functions to estimate mode and hdi
estimate_mode <- function(x) {
  # Remove NA values
  x <- na.omit(x)
  
  # Kernel density estimation
  density_est <- density(x)
  
  # Find the mode as the x-value corresponding to the maximum density
  mode_est <- density_est$x[which.max(density_est$y)]
  
  return(mode_est)
}
estimate_mode_and_hdi<- function(x, prob = 0.95) {
  x <- sort(na.omit(x))  # Sort the data and remove NA values
  n <- length(x)
  interval_size <- floor(prob * n)  # Number of data points in the interval
  
  # Find the shortest interval
  min_width <- Inf
  hdi_start <- 0
  hdi_end <- 0
  
  for (i in 1:(n - interval_size)) {
    width <- x[i + interval_size] - x[i]
    if (width < min_width) {
      min_width <- width
      hdi_start <- x[i]
      hdi_end <- x[i + interval_size]
    }
  }
  
  # Find the mode (most frequent value or midpoint of the data distribution)
  density_est <- density(x)
  mode_est <- density_est$x[which.max(density_est$y)]
  
  return(list(mode = mode_est, hdi = c(hdi_start, hdi_end)))
}


#' Compute KL Divergence with Bootstrap Confidence Interval
#'
#' Calculates the Kullback-Leibler (KL) divergence between an empirical distribution (via kernel density estimation)
#' and a fitted parametric distribution (Gumbel or Lognormal). Uses parametric bootstrap to estimate confidence intervals.
compute_kl_bootstrap <- function(fit, data, n_boot = 1000, conf_level = 0.95) {
  library(evd)
  library(univariateML)
  
  # Determine the model type
  model_type <- attr(fit, "model")
  
  if (model_type == "Gumbel") {
    mu_hat <- fit["mu"]
    sigma_hat <- fit["sigma"]
  } else if (model_type == "Lognormal") {
    meanlog_hat <- fit["meanlog"]
    sdlog_hat <- fit["sdlog"]
  } else {
    stop("Unsupported distribution type. Only Gumbel and Lognormal are supported.")
  }
  
  # Compute empirical density using kernel density estimation (KDE)
  empirical_density <- density(data)
  
  # Function to compute KL divergence given parameters
  compute_kl <- function(param1, param2) {
    if (model_type == "Gumbel") {
      theoretical_density <- dgumbel(empirical_density$x, loc = param1, scale = param2)
    } else if (model_type == "Lognormal") {
      theoretical_density <- dlnorm(empirical_density$x, meanlog = param1, sdlog = param2)
    }
    
    # Normalize densities
    dx <- diff(empirical_density$x)[1]  # Approximate bin width
    empirical_density$y <- empirical_density$y / sum(empirical_density$y * dx)
    theoretical_density <- theoretical_density / sum(theoretical_density * dx)
    
    # Compute KL divergence
    kl_div <- sum(empirical_density$y * log(empirical_density$y / theoretical_density), na.rm = TRUE) * dx
    return(kl_div)
  }
  
  # Compute KL divergence for the MLE estimates
  if (model_type == "Gumbel") {
    kl_mle <- compute_kl(mu_hat, sigma_hat)
  } else if (model_type == "Lognormal") {
    kl_mle <- compute_kl(meanlog_hat, sdlog_hat)
  }
  
  # Compute bootstrap samples for KL divergence
  set.seed(123)  # Ensure reproducibility
  boot_samples <- bootstrapml(fit, reps = n_boot, map = identity, reducer = identity)
  
  # Extract bootstrapped parameter values
  param1_boot <- boot_samples[1, ]
  param2_boot <- boot_samples[2, ]
  
  # Apply KL computation to each bootstrap sample
  kl_boot_values <- mapply(compute_kl, param1_boot, param2_boot)
  
  # Compute confidence interval
  alpha <- (1 - conf_level) / 2
  kl_CI <- quantile(kl_boot_values, probs = c(alpha, 1 - alpha))
  
  # Print results
  cat("Kullback-Leibler Divergence (KL):", round(kl_mle, 6), "\n")
  cat(conf_level * 100, "% Confidence Interval for KL:", 
      round(kl_CI[1], 6), "-", round(kl_CI[2], 6), "\n")
  
  # Return results as a list
  return(list(
    KL_MLE = kl_mle,
    KL_CI = kl_CI
  ))
}


#' Overlay Gumbel Fit and Bootstrap Confidence Band on Density Plot
#'
#' Plots the maximum likelihood Gumbel distribution fit over an optional data density,
#' and overlays a pointwise 95% confidence band based on bootstrapped parameter estimates.
#' Optionally includes a legend with fitted parameter values and KL divergence results.
overlay_gumbel_fit <- function(fit_gumbel,
                               # univariateML object with final MLE
                               boot_gumbel_raw,
                               # matrix of shape (2, reps) from bootstrapml(..., map=identity, reducer=identity)
                               x_vals = NULL,
                               col = "red",
                               lty = 2,
                               lwd = 2,
                               addLegend = TRUE,
                               # Optional data for raw density
                               data = NULL,
                               addDataDensity = FALSE,
                               dataCol = "gray",
                               dataLty = 1,
                               dataLwd = 2,
                               # New arguments to allow customization of the CI polygon
                               bandFill = rgb(1, 0, 0, 0.2),
                               bandBorder = NA,
                               legend_cex = 1) {
  # We use 'evd' for dgumbel()
  library(evd)
  
  ## Helper: check if a plot region is active
  is_plot_active <- function() {
    if (length(dev.list()) == 0) return(FALSE)
    usr_vals <- par("usr")
    if (is.null(usr_vals)) return(FALSE)
    # If usr is c(0,1,0,1), it likely means no high-level plot
    if (all(usr_vals == c(0, 1, 0, 1))) return(FALSE)
    TRUE
  }
  
  # 1) Extract the best-fit Gumbel parameters (mu, sigma)
  mu_hat    <- fit_gumbel["mu"]
  sigma_hat <- fit_gumbel["sigma"]
  
  # 2) Optional raw data density
  density_data <- NULL
  if (addDataDensity && !is.null(data)) {
    density_data <- density(data)
  }
  
  # 3) Determine x-range for the Gumbel fit
  #    Default around best-fit mu +/- 4*sigma
  default_range <- c(mu_hat - 4 * sigma_hat, mu_hat + 4 * sigma_hat)
  data_xrange   <- if (!is.null(density_data)) range(density_data$x) else c(Inf, -Inf)
  
  # If x_vals not provided, unify default_range & data_xrange or use existing plot
  if (is.null(x_vals)) {
    if (is_plot_active()) {
      usr_vals <- par("usr")  # c(x1, x2, y1, y2)
      x_vals <- seq(usr_vals[1], usr_vals[2], length.out = 200)
    } else {
      x_min <- min(default_range[1], data_xrange[1])
      x_max <- max(default_range[2], data_xrange[2])
      if (x_min == Inf || x_max == -Inf) {
        x_min <- default_range[1]
        x_max <- default_range[2]
      }
      x_vals <- seq(x_min, x_max, length.out = 200)
    }
  }
  
  # 4) Compute the best-fit Gumbel density
  y_vals_gumbel <- dgumbel(x_vals, loc = mu_hat, scale = sigma_hat)
  
  # 5) Build the pointwise confidence band from boot_gumbel_raw
  #    boot_gumbel_raw is shape (2, reps) => each col = (mu_i, sigma_i)
  #    We'll compute density for each iteration, store in a matrix.
  
  # Number of bootstrap reps
  reps <- ncol(boot_gumbel_raw)
  
  # Matrix to hold densities: rows = length(x_vals), cols = reps
  dens_matrix <- matrix(NA, nrow = length(x_vals), ncol = reps)
  
  # For each bootstrap sample, compute the entire Gumbel density
  for (b in seq_len(reps)) {
    mu_b    <- boot_gumbel_raw[1, b]
    sigma_b <- boot_gumbel_raw[2, b]
    dens_matrix[, b] <- dgumbel(x_vals, loc = mu_b, scale = sigma_b)
  }
  
  # For each x, compute the 2.5% and 97.5% quantiles across the reps
  # => a pointwise confidence band
  dens_lower <- apply(dens_matrix, 1, quantile, probs = 0.025)
  dens_upper <- apply(dens_matrix, 1, quantile, probs = 0.975)
  
  # 6) If no plot is active, create one sized to fit everything
  plot_active <- is_plot_active()
  if (!plot_active) {
    # We'll figure out y_max from best-fit + band + optional data
    y_max_candidates <- c(y_vals_gumbel, dens_lower, dens_upper)
    if (!is.null(density_data)) {
      y_max_candidates <- c(y_max_candidates, density_data$y)
    }
    y_max <- max(y_max_candidates, na.rm = TRUE)
    
    plot(
      range(x_vals),
      c(0, y_max),
      type = "n",
      xlab = "Log Evolutionary Rate",
      ylab = "Density",
      main = "Gumbel Fit with 95% CI (Pointwise)"
    )
  }
  
  # 7) Draw optional data density
  if (!is.null(density_data)) {
    lines(
      density_data$x,
      density_data$y,
      col = dataCol,
      lty = dataLty,
      lwd = dataLwd
    )
  }
  
  # 8) Draw the pointwise confidence band
  polygon(
    c(x_vals, rev(x_vals)),
    c(dens_lower, rev(dens_upper)),
    col = bandFill,
    border = bandBorder
  )
  
  # 9) Overlay the best-fit Gumbel line
  lines(
    x_vals,
    y_vals_gumbel,
    col = col,
    lty = lty,
    lwd = lwd
  )
  
  # 10) Legend if requested
  if (addLegend) {
    # Compute KL Divergence
    gumbel_kl <- compute_kl_bootstrap(fit_gumbel, data = data, n_boot = 10000)
    
    if (!is.null(density_data)) {
      # Legend labels: Data, MLE, 95% MLE CI, (half-space gap), KLD (with D subscript)
      legend_labels <- c(
        "Data",
        sprintf("MLE: \u03BC = %.2f, \u03C3 = %.2f", mu_hat, sigma_hat),
        "95% MLE CI",
        " ",  # 1/2 a space inserted here
        bquote(KL[D] ~ ":" ~ .(sprintf("%.2f", gumbel_kl$KL_MLE)) * " [" *
                 .(sprintf("%.2f", gumbel_kl$KL_CI[1])) * " - " *
                 .(sprintf("%.2f", gumbel_kl$KL_CI[2])) * "]")
      )
      
      # Corresponding styles
      # 1) Data: line only
      # 2) MLE: line only
      # 3) 95% MLE CI: small colored square
      # 4) Gap: no symbol/line (for spacing)
      # 5) KLD: line with no symbol
      legend_cols  <- c(dataCol, col, bandFill, NA, "black")
      legend_ltys  <- c(dataLty, lty, NA, NA, NA)
      legend_lwds  <- c(dataLwd, lwd, NA, NA, NA)
      legend_pchs  <- c(NA, NA, 15, NA, NA)  # Removed symbol for KLD
      legend_fills <- c(NA, NA, NA, NA, NA)  # No patch fill in legend
    } else {
      # Legend labels: MLE, 95% MLE CI, (half-space gap), KLD (with D subscript)
      legend_labels <- c(
        sprintf("MLE: \u03BC = %.2f, \u03C3 = %.2f", mu_hat, sigma_hat),
        "95% MLE CI",
        " ",  # 1/2 a space inserted here
        bquote(bold(D[KL]) ~ ":" ~ .(sprintf("%.2f", gumbel_kl$KL_MLE)) * " [" *
                 .(sprintf("%.2f", gumbel_kl$KL_CI[1])) * " - " *
                 .(sprintf("%.2f", gumbel_kl$KL_CI[2])) * "]")
      )
      
      # Corresponding styles
      legend_cols  <- c(col, bandFill, NA, "black")
      legend_ltys  <- c(lty, NA, NA, NA)
      legend_lwds  <- c(lwd, NA, NA, NA)
      legend_pchs  <- c(NA, 15, NA, NA)  # Removed symbol for KLD
      legend_fills <- c(NA, NA, NA, NA)
    }
    
    legend("topright",
           legend  = legend_labels,
           col     = legend_cols,
           lty     = legend_ltys,
           lwd     = legend_lwds,
           pch     = legend_pchs,
           fill    = legend_fills,
           pt.cex  = 1.2,    # Adjust symbol size if desired
           bty     = "n",    # No box around legend
           border  = NA,
           cex = legend_cex)
  }
  
  
}


#' Visualize Rate Change Nodes on a Phylogenetic Tree with Scaled Node Labels
#'
#' Plots node labels on a phylogenetic tree using circles scaled by the magnitude of rate changes 
#' (e.g., evolutionary rate increases or decreases). Optional letter annotations are added to identify nodes.
#' Useful for visualizing where and how strongly rates shift in a tree based on results from 
#' `extractMaxAgeOfRegimesWithRateChanges()`.
plotNodeLabelsByRateChange <- function(
    tree, 
    rate_data, 
    filter_by = "increase",  # Options: "increase", "decrease", "all"
    color_increase = "red", 
    color_decrease = "blue", 
    base_cex = 2,  # Base size of nodes
    scale_factor = 0.5,  # Controls scaling effect
    alpha = 0.7,  # Transparency of node labels
    transform_method = "log1p",  # Options: "log1p" or "square_root"
    plot_letters = TRUE,  # Toggle letter overlay (default = TRUE)
    letter_cex = 0.8,  # Size of the letter labels
    letter_color = "black"  # Color of the letter labels
) {
  # Step 1: Extract Rate Change Data
  tmp <- extractMaxAgeOfRegimesWithRateChanges(rate_data)
  
  # Step 2: Filter Data Based on User Choice
  if (filter_by != "all") {
    tmp <- tmp[tmp$rate_change == filter_by, ]
  }
  
  # Step 3: Apply Transformation for Scaling
  if (transform_method == "log1p") {
    tmp$transformed_size <- base_cex + (log1p(abs(tmp$percentage_change)) * scale_factor)
  } else if (transform_method == "square_root") {
    tmp$transformed_size <- base_cex + ((abs(tmp$percentage_change))^(1/2) * scale_factor)
  } else {
    stop("Invalid transform_method. Choose 'log1p' or 'square_root'.")
  }
  
  # Step 4: Define Colors Based on Rate Change Direction
  tmp$color <- ifelse(tmp$rate_change == "increase", color_increase, color_decrease)
  
  # Step 5: Plot Node Labels with Scaled Circles
  nodelabels(
    node = tmp$max_age_node,
    pch = 21,  
    cex = tmp$transformed_size,  
    bg = adjustcolor(tmp$color, alpha.f = alpha)  
  )
  
  # Step 6: Overlay Letters on Top of the Circles (if enabled)
  if (plot_letters) {
    # Assign unique letters to each node
    n <- nrow(tmp)
    if (n > length(LETTERS)) {
      stop("More nodes than available unique letters. Consider a different assignment method.")
    }
    tmp$letter <- LETTERS[1:n]
    
    # Overlay letters by calling nodelabels again
    nodelabels(
      node = tmp$max_age_node,
      text = tmp$letter,
      frame = "none",
      cex = letter_cex,
      col = letter_color
    )
  }
  
  # Step 7: Create and print a report mapping node numbers to details
  if (plot_letters) {
    report <- data.frame(
      Node = tmp$max_age_node,
      Letter = tmp$letter,
      PercentageChange = tmp$percentage_change,
      ScaledSize = tmp$transformed_size
    )
    print(report)
    return(report)
  } else {
    return(NULL)
  }
}


#' Highlight Low-Confidence Shift Nodes on a Phylogenetic Tree
#'
#' Plots nodes on a phylogenetic tree that have low confidence in rate shifts,
#' based on a threshold applied to `ic_weight_withshift`. These nodes are visually 
#' marked using the `nodelabels()` function for further inspection.
plotLowConfidenceShiftNodes <- function(
    tree,
    ic_weights,
    threshold = 0.9,
    dot_cex = 0.5,
    dot_color = "grey",
    pch_select = 19
) {
  low_conf_nodes <- ic_weights[ic_weights$ic_weight_withshift < threshold, "node"]
  
  if (length(low_conf_nodes) == 0) {
    message("No nodes found with ic_weight_withshift below threshold.")
    return(invisible(NULL))
  }
  
  nodelabels(
    node = low_conf_nodes,
    pch = pch_select,
    cex = dot_cex,
    col = dot_color
  )
  
  return(low_conf_nodes)
}


#' Extract Clade and Sister Clade Information from Key Nodes in a Phylogenetic Tree
#'
#' For a given set of internal nodes (key nodes) in a phylogenetic tree, this function extracts:
#' - Tip members of the clade
#' - Unique families within the clade (if `family` column is present)
#' - Tip members and families of the sister clade (if it exists)
extract_clade_info <- function(tree_data, keynodes, cores = 2) {
  results <- pbmclapply(keynodes, function(node) {
    # Extract clade members for the key node
    clade_members <- caper::clade.members(phy = tree_data$tree_no_uncertainty_untransformed, x = node, tip.labels = TRUE)
    clade_members_df <- as.data.frame(name_backbone_checklist(clade_members))
    
    # Extract unique families for the key node
    unique_families <- if ("family" %in% colnames(clade_members_df)) {
      unique(na.omit(clade_members_df$family))
    } else {
      character(0)
    }
    
    # Get the sister node
    sister_node <- phytools::getSisters(tree_data$tree_no_uncertainty_untransformed, node, mode = "number")
    
    # Initialize empty variables in case no sister exists
    sister_clade_members_df <- data.frame()
    sister_unique_families <- character(0)
    
    # Process the sister node if it exists
    if (!is.null(sister_node) && length(sister_node) > 0) {
      sister_clade_members <- caper::clade.members(phy = tree_data$tree_no_uncertainty_untransformed, x = sister_node, tip.labels = TRUE)
      sister_clade_members_df <- as.data.frame(name_backbone_checklist(sister_clade_members))
      
      sister_unique_families <- if ("family" %in% colnames(sister_clade_members_df)) {
        unique(na.omit(sister_clade_members_df$family))
      } else {
        character(0)
      }
    }
    
    # Return results as a list for each node
    list(
      clade_members = clade_members_df,
      unique_families = unique_families,
      sister_clade_members = sister_clade_members_df,
      sister_unique_families = sister_unique_families
    )
  }, mc.cores = cores)
  
  # Assign names to the results
  names(results) <- paste0("Node_", keynodes)
  
  # Split results into separate lists
  clade_members_list <- lapply(results, `[[`, "clade_members")
  unique_families_list <- lapply(results, `[[`, "unique_families")
  sister_clade_members_list <- lapply(results, `[[`, "sister_clade_members")
  sister_unique_families_list <- lapply(results, `[[`, "sister_unique_families")
  
  return(list(
    clade_members = clade_members_list,
    unique_families = unique_families_list,
    sister_clade_members = sister_clade_members_list,
    sister_unique_families = sister_unique_families_list
  ))
}


#Functions for Leave-one-out spatial analyses
#' Leave-One-Term-Out Model Comparison for SAR Models
#'
#' Performs model comparison by iteratively removing individual predictor terms from the full
#' spatial autoregressive (SAR) model and evaluating changes in AIC, BIC, and model weights.
#' Supports optional recomputation of the full model for consistency.
{
loo_model_weights <- function(model, data, listw, 
                              control_terms = c("(Intercept)", 
                                                "asin(sqrt(sampling_frac))", 
                                                "tipDR_range_weighted_mean_z"),
                              response_var = "mean_lineage_rate_range_weighted",
                              n_terms = NULL,
                              recompute_full_model = FALSE) {
  # Build the design matrix from the original model and data
  mm <- model.matrix(formula(model), data = data)
  
  # Get expanded predictor names
  all_terms_expanded <- colnames(mm)
  
  # Define candidate columns: all columns except control terms
  candidate_cols <- setdiff(all_terms_expanded, control_terms)
  
  # Optionally limit candidate columns to first n_terms if specified
  if (!is.null(n_terms)) {
    candidate_cols <- candidate_cols[1:min(n_terms, length(candidate_cols))]
  }
  
  # Extract the response variable from the data
  resp <- data[[response_var]]
  
  # Optionally refit the full model to ensure consistency
  if (recompute_full_model) {
    # Construct full model data the same way as reduced models
    full_data <- as.data.frame(mm)  # Use design matrix
    full_data <- full_data[, setdiff(colnames(full_data), "(Intercept)"), drop = FALSE]  # Drop intercept
    full_data[[response_var]] <- resp  # Add response variable
    
    full_model <- lagsarlm(as.formula(paste(response_var, "~ .")),
                           data = full_data,
                           listw = listw,
                           zero.policy = TRUE,
                           method = "eigen")
    full_AIC <- AIC(full_model)
    full_BIC <- BIC(full_model)
  } else {
    full_AIC <- AIC(model)
    full_BIC <- BIC(model)
    full_model <- model  # Keep the original model if not recomputed
  }
  
  # Helper function: remove one column, refit the model, and return a list with AIC, BIC, and model object
  fit_after_removal <- function(term) {
    new_cols <- setdiff(colnames(mm), term)
    new_data <- as.data.frame(mm[, new_cols, drop = FALSE])
    # Remove the intercept column (it is added automatically by the formula)
    new_data <- new_data[, setdiff(colnames(new_data), "(Intercept)"), drop = FALSE]
    new_data[[response_var]] <- resp
    updated_model <- lagsarlm(as.formula(paste(response_var, "~ .")),
                              data = new_data,
                              listw = listw,
                              zero.policy = TRUE,
                              method = "eigen")
    list(AIC = AIC(updated_model), BIC = BIC(updated_model), model = updated_model)
  }
  
  # Use pblapply (parallel with progress bar) to iterate over candidate columns
  res_list <- pblapply(candidate_cols, fit_after_removal)
  
  # Extract AIC and BIC values for candidate models
  aic_vals <- sapply(res_list, function(x) x$AIC)
  bic_vals <- sapply(res_list, function(x) x$BIC)
  
  results <- data.frame(
    Term = candidate_cols,
    AIC = aic_vals,
    BIC = bic_vals,
    stringsAsFactors = FALSE
  )
  
  # Compute ΔAIC and ΔBIC relative to the full model's AIC and BIC
  results$DeltaAIC <- results$AIC - full_AIC
  results$DeltaBIC <- results$BIC - full_BIC
  
  # Compute weights using the aicw function (for both AIC and BIC)
  results$AIC_Weight <- sapply(results$AIC, function(candidate_AIC) 
    aicw(c(full_AIC, candidate_AIC))$aicweights[1]
  )
  results$BIC_Weight <- sapply(results$BIC, function(candidate_BIC) 
    aicw(c(full_BIC, candidate_BIC))$aicweights[1]
  )
  
  # Create a named list of the updated model objects
  model_list <- lapply(res_list, function(x) x$model)
  names(model_list) <- candidate_cols
  
  # Prepend the recomputed full model (if applicable)
  if (recompute_full_model) {
    model_list <- c(list("Full Model" = full_model), model_list)
  }
  
  # Add full model metrics to the results for reference
  full_metrics <- data.frame(
    Term = "Full Model",
    AIC = full_AIC,
    BIC = full_BIC,
    DeltaAIC = 0,
    DeltaBIC = 0,
    AIC_Weight = 1,  # Full model always gets weight 1 in its own comparison
    BIC_Weight = 1,
    stringsAsFactors = FALSE
  )
  
  results <- rbind(full_metrics, results)
  results <- results[order(results$AIC), ]
  
  list(metrics = results, models = model_list)
}



loo_model_weights_from_obj <- function(mod, listw, 
                                       control_terms = c("(Intercept)", 
                                                         "asin(sqrt(sampling_frac))", 
                                                         "tipDR_range_weighted_mean_z"),
                                       response_var = "y",  # assumes mod$y holds the response
                                       n_terms = NULL) {
  # Extract the design matrix and response from the model object
  X <- mod$X
  y <- mod$y
  
  # Get all column names from the design matrix
  all_terms <- colnames(X)
  
  # Define candidate columns by removing the control terms
  candidate_cols <- setdiff(all_terms, control_terms)
  
  # Optionally limit candidate columns to the first n_terms if specified
  if (!is.null(n_terms)) {
    candidate_cols <- candidate_cols[1:min(n_terms, length(candidate_cols))]
  }
  
  # Full model metrics (from the passed model object)
  full_AIC <- AIC(mod)
  full_BIC <- BIC(mod)
  full_model <- mod
  
  # Helper function: remove one candidate column, then refit the model
  fit_after_removal <- function(term) {
    new_cols <- setdiff(all_terms, term)
    new_data <- as.data.frame(X[, new_cols, drop = FALSE])
    # Remove the intercept column if present (it will be added automatically)
    new_data <- new_data[, setdiff(colnames(new_data), "(Intercept)"), drop = FALSE]
    new_data[[response_var]] <- y  # Add the response vector
    updated_model <- lagsarlm(as.formula(paste(response_var, "~ .")),
                              data = new_data,
                              listw = listw,
                              zero.policy = mod$zero.policy,
                              method = mod$method)
    list(AIC = AIC(updated_model), BIC = BIC(updated_model), model = updated_model)
  }
  
  # Use pblapply to iterate over candidate columns with a progress bar
  res_list <- pbapply::pblapply(candidate_cols, fit_after_removal)
  
  # Extract AIC and BIC values from each reduced model
  aic_vals <- sapply(res_list, function(x) x$AIC)
  bic_vals <- sapply(res_list, function(x) x$BIC)
  
  results <- data.frame(
    Term = candidate_cols,
    AIC = aic_vals,
    BIC = bic_vals,
    stringsAsFactors = FALSE
  )
  
  # Compute ΔAIC and ΔBIC relative to the full model's AIC/BIC
  results$DeltaAIC <- results$AIC - full_AIC
  results$DeltaBIC <- results$BIC - full_BIC
  
  # Compute weights using the aicw function (ensure aicw() is defined/loaded)
  results$AIC_Weight <- sapply(results$AIC, function(candidate_AIC) 
    aicw(c(full_AIC, candidate_AIC))$aicweights[1]
  )
  results$BIC_Weight <- sapply(results$BIC, function(candidate_BIC) 
    aicw(c(full_BIC, candidate_BIC))$aicweights[1]
  )
  
  # Create a named list of the updated model objects from the reduced fits
  model_list <- lapply(res_list, function(x) x$model)
  names(model_list) <- candidate_cols
  
  # Prepend the full model to the list of models
  model_list <- c(list("Full Model" = full_model), model_list)
  
  # Add full model metrics for reference
  full_metrics <- data.frame(
    Term = "Full Model",
    AIC = full_AIC,
    BIC = full_BIC,
    DeltaAIC = 0,
    DeltaBIC = 0,
    AIC_Weight = 1,
    BIC_Weight = 1,
    stringsAsFactors = FALSE
  )
  
  results <- rbind(full_metrics, results)
  results <- results[order(results$AIC), ]
  
  list(metrics = results, models = model_list)
}
}


#More simulation functions
{

#' Subsample a Phylogenetic Tree to a Specified Number of Tips
#'
#' Randomly samples a specified number of tips from a base phylogenetic tree and returns
#' a new tree containing only those tips. The tree is reordered in postorder to ensure consistent structure.
treesampler <- function(base_tree, sample){keep.tip(reorder(as.phylo(base_tree), order = 'postorder'), tip = sample(base_tree$tip.label, sample))}


#' Simulate and Paint Multiple Trait Shifts on a Phylogenetic Tree
#'
#' Simulates multiple evolutionary rate shifts in a phylogenetic tree by identifying eligible clades,
#' painting regimes, constructing ancestral and derived covariance matrices, and simulating multivariate trait data
#' under a multi-rate Brownian Motion model (BMM).
simulateAndPaintMultipleShifts <- function(
  tree,
  minCladeSize        = 10,
  maxCladeSize        = 50,
  dimensions          = 10,
  variance            = 1,       # mean for diagonal variance
  covariance          = 0.5,     # mean for off‐diagonal covariance
  variance_sd         = 0.1,     # SD around the variance mean
  covariance_sd       = 0.25,    # SD around the covariance mean
  scaleFactor         = "sample",
  scaleFactorRange    = c(0.1, 2.0),
  excludeRange        = c(0.9, 1.1),
  numShifts           = 2,
  maxAttempts         = 100,
  plot                = FALSE,
  buffer              = 5,
  scaleCovariance     = FALSE
) {
    if (!inherits(tree, "phylo")) {
      stop("The tree must be an object of class 'phylo'.")
    }
    
    require(RRphylo)  # Ensure the RRphylo package is loaded
    
    outerAttempt <- 0
    while (outerAttempt < maxAttempts) {
      # ──────────────────────────────────────────────────────────────────────────
      # 1. Paint the “root” ancestral regime
      # ──────────────────────────────────────────────────────────────────────────
      numTips <- length(tree$tip.label)
      root <- numTips + 1
      tree <- paintSubTree(tree, node = root, state = "ancestral")
      
      # ──────────────────────────────────────────────────────────────────────────
      # 2. Identify all valid candidate nodes (clade size between min/max)
      # ──────────────────────────────────────────────────────────────────────────
      validCandidates <- setdiff(tree$edge[, 1], c(tree$tip.label, root))
      validCandidates <- Filter(function(node) {
        descendants <- getDescendants(tree, node)
        numDescendants <- length(descendants[descendants <= numTips])
        return(numDescendants >= minCladeSize && numDescendants <= maxCladeSize)
      }, validCandidates)
      
      if (length(validCandidates) < numShifts) {
        cat(sprintf("Attempt %d: Not enough valid nodes meet the clade size requirements. Restarting...\n",
                    outerAttempt + 1))
        outerAttempt <- outerAttempt + 1
        next
      }
      
      # ──────────────────────────────────────────────────────────────────────────
      # 3. Randomly sample `numShifts` distinct shift nodes, enforcing no overlap
      # ──────────────────────────────────────────────────────────────────────────
      shiftNodes <- integer(0)
      for (shiftIndex in seq_len(numShifts)) {
        if (length(validCandidates) == 0) {
          cat("No valid candidates left for further shifts. Restarting...\n")
          outerAttempt <- outerAttempt + 1
          next
        }
        randomNode <- sample(validCandidates, 1)
        shiftNodes <- c(shiftNodes, randomNode)
        
        descendants <- getDescendants(tree, randomNode)
        ancestors   <- getMommy(tree, randomNode)
        validCandidates <- setdiff(validCandidates, c(descendants, ancestors, randomNode))
      }
      
      # ──────────────────────────────────────────────────────────────────────────
      # 4. Enforce minimum “buffer” (phylogenetic distance) between any two shift nodes
      # ──────────────────────────────────────────────────────────────────────────
      allDistancesSufficient <- TRUE
      if (length(shiftNodes) > 1) {
        for (i in seq_len(length(shiftNodes) - 1)) {
          for (j in (i + 1):length(shiftNodes)) {
            distanceInfo <- distNodes(tree, node = c(shiftNodes[i], shiftNodes[j]), clus = 0)
            if (distanceInfo$node < buffer) {
              cat(sprintf(
                "Buffer check failed: Distance between node %d and node %d is %d, less than buffer %d. Restarting...\n",
                shiftNodes[i], shiftNodes[j], distanceInfo$node, buffer
              ))
              allDistancesSufficient <- FALSE
              break
            }
          }
          if (!allDistancesSufficient) {
            break
          }
        }
      }
      
      if (!allDistancesSufficient) {
        outerAttempt <- outerAttempt + 1
        next
      }
      
      # ──────────────────────────────────────────────────────────────────────────
      # 5. Build a positive‐definite “ancestral” covariance matrix
      # ──────────────────────────────────────────────────────────────────────────
      ntraits <- dimensions
      ancestralSigma <- NULL
      vcvAttempt <- 0
      
      while (vcvAttempt < maxAttempts) {
        # 5a. Sample diagonal variances
        sampledVariances <- abs(rnorm(ntraits, mean = variance, sd = variance_sd))
        # 5b. Build a random “base” for off‐diagonal elements
        L <- matrix(rnorm(ntraits * ntraits, mean = covariance, sd = covariance_sd),
                    ncol = ntraits)
        L[upper.tri(L)] <- 0
        sigmaBase <- t(L) %*% L
        diag(sigmaBase) <- sigmaBase[cbind(seq_len(ntraits), seq_len(ntraits))] + sampledVariances
        
        ev <- eigen(sigmaBase, symmetric = TRUE, only.values = TRUE)$values
        if (all(ev > 0)) {
          ancestralSigma <- sigmaBase
          break
        }
        vcvAttempt <- vcvAttempt + 1
      }
      
      if (is.null(ancestralSigma)) {
        cat("Failed to generate positive definite ancestral covariance matrix after", maxAttempts, "attempts. Restarting...\n")
        outerAttempt <- outerAttempt + 1
        next
      }
      
      # ──────────────────────────────────────────────────────────────────────────
      # 6. Initialize lists of VCVs (ancestral + derived)
      # ──────────────────────────────────────────────────────────────────────────
      VCVs <- list(ancestral = ancestralSigma)
      sigmaList <- list(ancestral = ancestralSigma)
      sampledScaleFactors <- numeric(numShifts)
      
      # ──────────────────────────────────────────────────────────────────────────
      # 7. For each shift, generate a “derived” covariance matrix
      # ──────────────────────────────────────────────────────────────────────────
      failedAll <- FALSE
      for (shiftCount in seq_len(numShifts)) {
        derivedStateName <- paste("derived", shiftCount, sep = "_")
        
        # 7a. Draw a scale‐factor `s` (unless fixed)
        if (scaleFactor == "sample") {
          part1 <- runif(1, min = scaleFactorRange[1], max = excludeRange[1])
          part2 <- runif(1, min = excludeRange[2], max = scaleFactorRange[2])
          currentScaleFactor <- sample(c(part1, part2), 1)
        } else {
          currentScaleFactor <- scaleFactor
        }
        sampledScaleFactors[shiftCount] <- currentScaleFactor
        
        # 7b. Attempt to construct a valid “derived” Sigma, retrying up to maxAttempts
        vcvDerivedAttempt <- 0
        repeat {
          if (scaleCovariance) {
            # Step i: extract ancestral correlation matrix
            corMat <- cov2cor(ancestralSigma)
            # Step ii: multiply only off‐diagonals by s
            scaledCor <- corMat
            scaledCor[lower.tri(scaledCor)] <- corMat[lower.tri(corMat)] * currentScaleFactor
            scaledCor[upper.tri(scaledCor)] <- t(scaledCor)[upper.tri(corMat)]
            # Step iii: adjust variances by 1/s
            sdVec <- sqrt(diag(ancestralSigma))
            newVars <- (sdVec^2) / currentScaleFactor
            # Step iv: rebuild the derived covariance
            derivedSigma <- diag(sqrt(newVars)) %*% scaledCor %*% diag(sqrt(newVars))
          } else {
            # “Standard” BMM: simply multiply entire ancestral Sigma by s
            derivedSigma <- currentScaleFactor * ancestralSigma
          }
          
          # Check positive definiteness
          ev <- eigen(derivedSigma, symmetric = TRUE, only.values = TRUE)$values
          if (all(ev > 0)) {
            break
          }
          
          vcvDerivedAttempt <- vcvDerivedAttempt + 1
          if (vcvDerivedAttempt >= maxAttempts) {
            cat(sprintf("Failed to generate valid derived matrix at shift %d after %d attempts. Restarting...\n",
                        shiftCount, maxAttempts))
            outerAttempt <- outerAttempt + 1
            failedAll <- TRUE
            break
          }
        }
        
        if (failedAll) {
          break
        }
        
        # 7c. Record the valid derived Sigma
        VCVs[[derivedStateName]]    <- derivedSigma
        sigmaList[[derivedStateName]] <- derivedSigma
        
        # 7d. Paint the subtree under the new regime
        tree <- paintSubTree(tree, node = shiftNodes[shiftCount],
                             state = derivedStateName, anc.state = "ancestral")
      }
      
      if (failedAll) {
        next
      }
      
      # ──────────────────────────────────────────────────────────────────────────
      # 8. Verify we have exactly (numShifts + 1) distinct regimes
      # ──────────────────────────────────────────────────────────────────────────
      if (length(unique(getStates(tree))) != numShifts + 1) {
        cat("State validation failed. Restarting the entire simulation procedure...\n")
        outerAttempt <- 0
        next
      }
      
      # ──────────────────────────────────────────────────────────────────────────
      # 9. Optionally plot the final painted tree
      # ──────────────────────────────────────────────────────────────────────────
      if (plot) {
        plot(tree, ftype = "off")
      }
      
      # ──────────────────────────────────────────────────────────────────────────
      # 10. Simulate multivariate trait data under a multi‐rate BM (“BMM”) model
      # ──────────────────────────────────────────────────────────────────────────
      simulatedData <- mvSIM(
        tree,
        nsim  = 1,
        model = "BMM",
        param = list(sigma = sigmaList, theta = rep(0, ntraits))
      )
      
      # ──────────────────────────────────────────────────────────────────────────
      # 11. Return everything
      # ──────────────────────────────────────────────────────────────────────────
      return(list(
        paintedTree         = tree,
        shiftNodes          = shiftNodes,
        simulatedData       = simulatedData,
        VCVs                = VCVs,
        sampledScaleFactors = if (scaleFactor == "sample") sampledScaleFactors else NULL
      ))
    }
    
    stop("Failed to place all required shifts after maximum attempts.")
}

  
#' Plot Histogram of False Positive (FP) Rates
#'
#' Creates a histogram of false positive (FP) rates from a list of numeric vectors, 
#' with a dashed vertical line and label indicating the mean FP rate.
plot_fp_histogram <- function(FP_list, 
                              title = "Histogram of FP", 
                              xlab = "FP rate", 
                              ylab = "Count",
                              mean_label_y_frac = 0.8,
                              x_offset_frac = 0.1) {
  
  data_vec <- unlist(FP_list)
  mean_val <- mean(data_vec)
  
  formatted_mean <- if (mean_val < 0.01) {
    formatC(mean_val, format = "e", digits = 2)
  } else {
    formatC(round(mean_val, 2), format = "f", digits = 2)
  }
  
  hist_data <- ggplot(mapping = aes(x = data_vec)) +
    geom_histogram(bins = 10)
  gb <- ggplot_build(hist_data)
  max_count <- max(gb$data[[1]]$count)
  label_y_pos <- mean_label_y_frac * max_count
  
  x_range <- range(data_vec)
  x_offset <- x_offset_frac * diff(x_range)
  label_x_pos <- mean_val + x_offset
  
  ggplot(mapping = aes(x = data_vec)) +
    geom_histogram(bins = 10, fill = "steelblue", color = "black") +
    geom_vline(xintercept = mean_val, 
               linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", 
             x = label_x_pos, 
             y = label_y_pos, 
             label = paste0("Mean = ", formatted_mean),
             vjust = -0.5, hjust = 0, size = 4, fontface = "italic") +
    labs(title = title, x = xlab, y = ylab) +
    coord_cartesian(clip = "off", expand = TRUE) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 12),  # reduced font size here
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12)
    )
}

#' Run False Positive Shift Inference on Simulated Datasets
#'
#' Simulates multiple trait datasets and runs a shift inference routine to estimate false positive rates
#' using parallel processing. Each dataset is generated with a subsampled tree and traits simulated under
#' a Brownian motion model with no shifts.
run_FP_shift_inference <- function(
    base_tree = NULL,
    n_datasets        = 100,
    n_traits          = 10,
    tree_tip_count    = 100,
    search_options    = list(),
    simulation_options = list(),  # New: args to pass to simulate_traits_BM1
    num_cores         = 2,
    seed              = 5
) {
  # Suppress fork warnings in RStudio
  Sys.setenv(R_PARALLELLY_SUPPORTSMULTICORE_UNSTABLE = "quiet")
  
  # Load required packages
  require(future)
  require(future.apply)
  require(progressr)
  require(ape)
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Normalize n_traits to length == n_datasets
  if (length(n_traits) != n_datasets) {
    if (length(n_traits) == 1) {
      n_traits <- rep(n_traits, n_datasets)
    } else {
      n_traits <- n_traits[seq_len(n_datasets)]
      warning("n_traits was truncated to match n_datasets.")
    }
  }
  
  # Simulate data (trees + traits)
  message("Simulating data...")
  simdata <- lapply(seq_len(n_datasets), function(i) {
    sim_args <- modifyList(simulation_options, list(
      n_species = NULL,
      n_traits  = n_traits[i],
      tree      = treesampler(base_tree, sample = tree_tip_count)
    ))
    do.call(simulate_traits_BM1, sim_args)
  })
  
  # Plan for multisession (RStudio-safe)
  future::plan(future::multisession, workers = num_cores)
  
  # Register a simple text progress bar
  progressr::handlers(progressr::handler_txtprogressbar)
  
  # Run inference with progress tracking
  message("Starting parallel inference...")
  results <- with_progress({
    p <- progressr::progressor(along = simdata)
    
    future_lapply(seq_along(simdata), function(i) {
      baseline_tree <- paintSubTree(
        reorder(as.phylo(simdata[[i]]$tree), order = "postorder"),
        node  = length(simdata[[i]]$tree$tip.label) + 1,
        state = 0
      )
      
      search_args <- modifyList(search_options, list(
        baseline_tree = baseline_tree,
        trait_data    = simdata[[i]]$data
      ))
      
      result <- {
        null_con <- if (.Platform$OS.type == "windows") {
          file("NUL", open = "w")
        } else {
          file("/dev/null", open = "w")
        }
        
        sink(null_con)
        sink(null_con, type = "message")
        out <- NULL
        try({
          out <- do.call(searchOptimalConfiguration, search_args)
        }, silent = TRUE)
        sink(type = "message")
        sink()
        close(null_con)
        out
      }
      
      p()
      result
    }, future.seed = TRUE)
  })
  
  return(results)
}

#' Run False Negative Shift Inference on Simulated Datasets with Known Shifts
#'
#' Simulates datasets with multiple rate shifts on phylogenetic trees using 
#' `simulateAndPaintMultipleShifts()` and performs inference to detect shifts 
#' using `searchOptimalConfiguration()`. Enables parallel computation.
run_FN_shift_inference <- function(
    base_tree = NULL,
    n_datasets         = 100,
    tree_tip_count     = 200,
    simulation_options = list(),
    search_options     = list(),
    num_cores          = 4,
    seed               = 5
) {
  #───────────────────────────────────────────────────────────────────────────────
  # 1. Suppress fork warnings in RStudio & load packages
  #───────────────────────────────────────────────────────────────────────────────
  Sys.setenv(R_PARALLELLY_SUPPORTSMULTICORE_UNSTABLE = "quiet")
  require(future)
  require(future.apply)
  require(progressr)
  require(ape)
  require(RRphylo)
  
  #───────────────────────────────────────────────────────────────────────────────
  # 2. Set seed for reproducibility
  #───────────────────────────────────────────────────────────────────────────────
  set.seed(seed)
  
  #───────────────────────────────────────────────────────────────────────────────
  # 3. Prepare default simulation parameters for simulateAndPaintMultipleShifts()
  #───────────────────────────────────────────────────────────────────────────────
  default_sim_opts <- list(
    minCladeSize     = 10,
    maxCladeSize     = 50,
    dimensions       = 10,
    variance         = 0.00025,
    covariance       = 2.25e-05,
    variance_sd      = 0.00014,
    covariance_sd    = 0.00012,
    scaleFactor      = "sample",
    scaleFactorRange = c(0.1, 2.0),
    excludeRange     = c(0.5, 1.5),
    numShifts        = 2,
    buffer           = 3,
    maxAttempts      = 100,
    plot             = FALSE
  )
  sim_opts <- modifyList(default_sim_opts, simulation_options)
  
  #───────────────────────────────────────────────────────────────────────────────
  # 4. SIMULATION: parallel generation of n_datasets of paintedTree + simulatedData
  #───────────────────────────────────────────────────────────────────────────────
  message("Simulating data with multiple rate shifts…")
  future::plan(future::multisession, workers = num_cores)
  
  simdata <- with_progress({
    p <- progressr::progressor(along = 1:n_datasets)
    
    future_lapply(seq_len(n_datasets), function(i) {
      this_tree <- treesampler(base_tree, sample = tree_tip_count)
      args      <- modifyList(sim_opts, list(tree = this_tree))
      result    <- do.call(simulateAndPaintMultipleShifts, args)
      p()
      result
    }, future.seed = TRUE)
  })
  
  # Reset plan to avoid conflicts between phases
  future::plan(future::sequential)
  
  #───────────────────────────────────────────────────────────────────────────────
  # 5. Prepare default inference parameters for searchOptimalConfiguration()
  #───────────────────────────────────────────────────────────────────────────────
  default_search_opts <- list(
    formula                    = "trait_data ~ 1",
    min_descendant_tips        = 10,
    num_cores                  = 1,
    ic_uncertainty_threshold   = 10,
    shift_acceptance_threshold = 10,
    uncertainty                = FALSE,
    uncertaintyweights         = FALSE,
    uncertaintyweights_par     = FALSE,
    postorder_traversal        = FALSE,
    plot                       = FALSE,
    IC                         = "GIC",
    store_model_fit_history    = FALSE
  )
  search_opts_template <- modifyList(default_search_opts, search_options)
  
  #───────────────────────────────────────────────────────────────────────────────
  # 6. INFERENCE: parallel inference on each simulated dataset
  #───────────────────────────────────────────────────────────────────────────────
  message("Running inference on simulated datasets…")
  future::plan(future::multisession, workers = num_cores)
  
  results <- with_progress({
    p <- progressr::progressor(along = simdata)
    
    future_lapply(seq_along(simdata), function(i) {
      painted_tree   <- simdata[[i]]$paintedTree
      simulated_data <- simdata[[i]]$simulatedData
      
      baseline_tree <- paintSubTree(
        reorder(as.phylo(painted_tree), order = "postorder"),
        node  = length(painted_tree$tip.label) + 1,
        state = 0
      )
      
      search_args <- modifyList(
        search_opts_template,
        list(
          baseline_tree = baseline_tree,
          trait_data    = simulated_data
        )
      )
      
      result <- {
        null_con <- if (.Platform$OS.type == "windows") file("NUL", open = "w") else file("/dev/null", open = "w")
        sink(null_con); sink(null_con, type = "message")
        out <- NULL
        try({ out <- do.call(searchOptimalConfiguration, search_args) }, silent = TRUE)
        sink(type = "message"); sink()
        close(null_con)
        out
      }
      
      p()
      result
    }, future.seed = TRUE)
  })
  
  future::plan(future::sequential)
  
  #───────────────────────────────────────────────────────────────────────────────
  # 7. RETURN: both the simulated datasets and the inference results
  #───────────────────────────────────────────────────────────────────────────────
  return(list(
    simdata = simdata,
    results = results
  ))
}


#' Evaluate Accuracy of Inferred Shifts Against Simulated Ground Truth
#'
#' Compares inferred shift nodes to known true shift nodes across multiple datasets
#' using strict (exact match) and fuzzy (within phylogenetic distance) matching.
#' Computes common metrics including precision, recall, F1 score, specificity, false positive rate,
#' and balanced accuracy. Optionally supports weighting based on information criterion weights.
evaluate_shift_recovery <- function(
    simdata, simresults,
    fuzzyDist = 2,
    verbose   = TRUE,
    weighted  = FALSE
) {
  require(RRphylo)
  require(ape)          # for nodepath()
  
  # ─────────────────────────────── helper utilities ───────────────────────────
  safe_div <- function(num, den) ifelse(den == 0, NA, num / den)
  harm     <- function(p, r)
    ifelse(is.na(p + r) || p + r == 0, NA, 2 * p * r / (p + r))
  
  # cumulative contingency tables
  strict <- c(TP = 0, FP = 0, FN = 0, TN = 0)
  fuzzy  <- c(TP = 0, FP = 0, FN = 0, TN = 0)
  
  # weighted sums (weights apply **only to predicted nodes**)
  w_strict <- c(TP = 0, FP = 0)
  w_fuzzy  <- c(TP = 0, FP = 0)
  
  # ───────────────────────────── iterate over datasets ─────────────────────────
  for (k in seq_along(simdata)) {
    true_nodes     <- simdata[[k]]$shiftNodes
    inferred_nodes <- simresults[[k]]$shift_nodes_no_uncertainty
    cand_k         <- simresults[[k]]$num_candidates
    tree_k         <- simdata[[k]]$paintedTree
    
    ## weights attached to **predicted** nodes (if available)
    weights <- NULL
    if (weighted && !is.null(simresults[[k]]$ic_weights)) {
      w_df   <- simresults[[k]]$ic_weights
      weights <- setNames(w_df$ic_weight_withshift, w_df$node)
    }
    
    # ────────── 1 · STRICT counts (exact‑node matches) ──────────
    TP_nodes <- intersect(true_nodes, inferred_nodes)
    FP_nodes <- setdiff(inferred_nodes, true_nodes)
    FN_nodes <- setdiff(true_nodes, inferred_nodes)
    TNk      <- cand_k - length(TP_nodes) - length(FP_nodes) - length(FN_nodes)
    
    strict <- strict + c(TP = length(TP_nodes),
                         FP = length(FP_nodes),
                         FN = length(FN_nodes),
                         TN = TNk)
    
    if (weighted && !is.null(weights)) {
      w_strict["TP"] <- w_strict["TP"] +
        sum(weights[as.character(TP_nodes)], na.rm = TRUE)
      w_strict["FP"] <- w_strict["FP"] +
        sum(weights[as.character(FP_nodes)], na.rm = TRUE)
    }
    
    # ────────── 2 · FUZZY counts (≤ fuzzyDist, one‑to‑one) ──────────
    if (length(true_nodes) == 0 || length(inferred_nodes) == 0) {
      
      ## trivial edge cases
      TPf <- 0
      FPf <- length(inferred_nodes)
      FNf <- length(true_nodes)
      
    } else {
      ## distance matrix (Inf if beyond fuzzyDist)
      dist_mat <- matrix(Inf,
                         nrow = length(inferred_nodes),
                         ncol = length(true_nodes))
      
      for (i in seq_along(inferred_nodes))
        for (j in seq_along(true_nodes)) {
          d <- tryCatch(length(nodepath(tree_k,
                                        inferred_nodes[i],
                                        true_nodes[j])) - 1,
                        error = function(e) Inf)
          if (d <= fuzzyDist) dist_mat[i, j] <- d
        }
      
      ## greedy one‑to‑one assignment (fine for small matrices)
      matched_inf  <- rep(FALSE, length(inferred_nodes))
      matched_true <- rep(FALSE, length(true_nodes))
      while (TRUE) {
        min_d <- min(dist_mat, na.rm = TRUE)
        if (!is.finite(min_d)) break                 # no more pairs within radius
        idx <- which(dist_mat == min_d, arr.ind = TRUE)[1, ]  # first minimum
        matched_inf[idx[1]]  <- TRUE
        matched_true[idx[2]] <- TRUE
        dist_mat[idx[1], ]   <- Inf                 # drop that inferred row
        dist_mat[ , idx[2]]  <- Inf                 # drop that true column
      }
      
      TPf <- sum(matched_inf)                       # one TP per true shift
      FPf <- length(inferred_nodes) - TPf
      FNf <- length(true_nodes) - sum(matched_true)
    }
    
    TNf <- cand_k - TPf - FPf - FNf
    fuzzy <- fuzzy + c(TP = TPf, FP = FPf, FN = FNf, TN = TNf)
    
    if (weighted && !is.null(weights)) {
      w_fuzzy["TP"] <- w_fuzzy["TP"] +
        sum(weights[as.character(inferred_nodes[matched_inf])], na.rm = TRUE)
      w_fuzzy["FP"] <- w_fuzzy["FP"] +
        sum(weights[as.character(inferred_nodes[!matched_inf])], na.rm = TRUE)
    }
  }
  
  # ────────────────────────── 3 · aggregate metrics ───────────────────────────
  ## strict
  prec <- safe_div(strict["TP"], strict["TP"] + strict["FP"])
  rec  <- safe_div(strict["TP"], strict["TP"] + strict["FN"])
  spec <- safe_div(strict["TN"], strict["TN"] + strict["FP"])
  f1   <- harm(prec, rec)
  fpr  <- safe_div(strict["FP"], strict["FP"] + strict["TN"])
  bal  <- mean(c(rec, spec), na.rm = TRUE)
  
  ## fuzzy
  f_prec <- safe_div(fuzzy["TP"], fuzzy["TP"] + fuzzy["FP"])
  f_rec  <- safe_div(fuzzy["TP"], fuzzy["TP"] + fuzzy["FN"])
  f_spec <- safe_div(fuzzy["TN"], fuzzy["TN"] + fuzzy["FP"])
  f_f1   <- harm(f_prec, f_rec)
  f_fpr  <- safe_div(fuzzy["FP"], fuzzy["FP"] + fuzzy["TN"])
  f_bal  <- mean(c(f_rec, f_spec), na.rm = TRUE)
  
  ## weighted (weights only on predictions)
  if (weighted) {
    wp_strict <- safe_div(w_strict["TP"], w_strict["TP"] + w_strict["FP"])
    wr_strict <- safe_div(w_strict["TP"], strict["TP"] + strict["FN"])
    wf1_strict <- harm(wp_strict, wr_strict)
    
    wp_fuzzy <- safe_div(w_fuzzy["TP"], w_fuzzy["TP"] + w_fuzzy["FP"])
    wr_fuzzy <- safe_div(w_fuzzy["TP"], fuzzy["TP"] + fuzzy["FN"])
    wf1_fuzzy <- harm(wp_fuzzy, wr_fuzzy)
  }
  
  # ────────────────────────────────── 4 · output ──────────────────────────────
  if (verbose) {
    cat("\nStrict Performance Metrics\n-------------------------\n")
    cat(sprintf("Precision        : %.3f\n", prec))
    cat(sprintf("Recall           : %.3f\n", rec))
    cat(sprintf("F1 Score         : %.3f\n", f1))
    cat(sprintf("Specificity      : %.3f\n", spec))
    cat(sprintf("False Pos. Rate  : %.3f\n", fpr))
    cat(sprintf("Balanced Accuracy: %.3f\n\n", bal))
    
    cat(sprintf("Fuzzy Matching (distance ≤ %d)\n----------------------------------------\n",
                fuzzyDist))
    cat(sprintf("Fuzzy Precision        : %.3f\n", f_prec))
    cat(sprintf("Fuzzy Recall           : %.3f\n", f_rec))
    cat(sprintf("Fuzzy F1 Score         : %.3f\n", f_f1))
    cat(sprintf("Fuzzy Specificity      : %.3f\n", f_spec))
    cat(sprintf("Fuzzy False Pos. Rate  : %.3f\n", f_fpr))
    cat(sprintf("Fuzzy Balanced Accuracy: %.3f\n\n", f_bal))
    
    if (weighted) {
      cat("Weighted Metrics (weights on predicted nodes)\n"
          ,"-------------------------------------------\n", sep = "")
      cat(sprintf("Weighted Precision (strict): %.3f\n", wp_strict))
      cat(sprintf("Weighted Recall    (strict): %.3f\n", wr_strict))
      cat(sprintf("Weighted F1 Score  (strict): %.3f\n", wf1_strict))
      cat(sprintf("Weighted Precision (fuzzy) : %.3f\n", wp_fuzzy))
      cat(sprintf("Weighted Recall    (fuzzy) : %.3f\n", wr_fuzzy))
      cat(sprintf("Weighted F1 Score  (fuzzy) : %.3f\n\n", wf1_fuzzy))
    }
  }
  
  invisible(list(
    strict = list(precision = prec, recall = rec, f1 = f1,
                  specificity = spec, fpr = fpr, balanced_accuracy = bal),
    fuzzy  = list(precision = f_prec, recall = f_rec, f1 = f_f1,
                  specificity = f_spec, fpr = f_fpr, balanced_accuracy = f_bal),
    weighted = if (weighted) list(
      strict = list(precision = wp_strict, recall = wr_strict, f1 = wf1_strict),
      fuzzy  = list(precision = wp_fuzzy,  recall = wr_fuzzy,  f1 = wf1_fuzzy)
    ),
    counts = list(strict = strict, fuzzy = fuzzy)
  ))
}


#' Generate LaTeX Table of Performance Metrics from Evaluation Results
#'
#' Formats and prints a LaTeX table of performance metrics (precision, recall, F1 score, etc.)
#' from one or more shift recovery evaluations using the \code{stargazer} package.
#' Supports both strict and fuzzy matching metrics, with and without IC-based weighting.
make_performance_table <- function(metrics,
                                   caption = "Shift Detection Performance Metrics",
                                   label   = "tab:performance",
                                   increment_counter = NULL) {
  
  if (!is.list(metrics[[1]])) metrics <- list(Main = metrics)
  if (is.null(names(metrics)) || any(names(metrics) == ""))
    names(metrics) <- LETTERS[seq_along(metrics)]
  
  r2 <- function(x) round(x, 2)           # helper → 2‑decimal rounding
  
  build_block <- function(m){
    data.frame(
      Metric = c("Precision","Recall","F1 Score",
                 "Specificity","False Positive Rate","Balanced Accuracy"),
      Strict            = r2(c(m$strict$precision,
                               m$strict$recall,
                               m$strict$f1,
                               m$strict$specificity,
                               m$strict$fpr,
                               m$strict$balanced_accuracy)),
      `Strict—Weighted` = r2(c(m$weighted$strict$precision,
                               m$weighted$strict$recall,
                               m$weighted$strict$f1,
                               NA, NA, NA)),
      Fuzzy             = r2(c(m$fuzzy$precision,
                               m$fuzzy$recall,
                               m$fuzzy$f1,
                               m$fuzzy$specificity,
                               m$fuzzy$fpr,
                               m$fuzzy$balanced_accuracy)),
      `Fuzzy—Weighted`  = r2(c(m$weighted$fuzzy$precision,
                               m$weighted$fuzzy$recall,
                               m$weighted$fuzzy$f1,
                               NA, NA, NA)),
      check.names = FALSE
    )
  }
  bigdf <- do.call(rbind, lapply(metrics, build_block))
  
  # capture stargazer output
  tex <- capture.output(
    stargazer::stargazer(
      bigdf,
      type            = "latex",
      summary         = FALSE,
      rownames        = FALSE,
      title           = caption,
      label           = label,
      digits          = 2,              # two‑decimal display
      float           = TRUE,
      font.size       = "small",
      table.placement = "!htbp"
    )
  )
  
  # insert rule + label in front of every block
  starts <- grep("^Precision\\s*&", tex)
  for (i in rev(seq_along(starts))) {
    tag <- names(metrics)[i]
    tex <- append(tex,
                  c("\\hline",
                    sprintf("\\multicolumn{5}{l}{\\textbf{%s}} \\\\", tag),
                    "\\hline"),
                  after = starts[i] - 1)
  }
  
  # optionally insert \setcounter
  if (!is.null(increment_counter)) {
    set_line <- sprintf("\\setcounter{table}{%d}", increment_counter)
    tex <- append(tex, set_line, after = which(grepl("^\\\\begin\\{table\\}", tex)))
  }
  
  cat(tex, sep = "\n")
  invisible(tex)
}


#' Count Unique States in Simulated Trees to Assess False Negatives
#'
#' This function processes a list of simulation results and counts the number
#' of unique painted states (regimes) in each result's untransformed tree.
#' If a result is missing a tree (i.e., \code{NULL}), it assumes no shifts
#' were detected and assigns a count of 1.
#'
#' This is typically used to assess the presence of false negatives in shift detection.
processFalseNegSimulations <- function(simResults) {
    # Initialize a vector to store the results
    uniqueStatesCounts <- vector("list", length(simResults))
    
    # Loop through the entire list of simulation results
    for (i in seq_along(simResults)) {
      if (!is.null(simResults[[i]]$tree_no_uncertainty_untransformed)) {
        # Calculate the number of unique states in the tree
        numUniqueStates <- length(unique(getStates(simResults[[i]]$tree_no_uncertainty_untransformed)))
        uniqueStatesCounts[[i]] <- numUniqueStates
      } else {
        # If the tree object is NULL, assume no shifts were detected, and assign 1
        uniqueStatesCounts[[i]] <- 1
        cat(sprintf("Element at index %d is NULL. No shifts detected, assigning 1 state.\n", i))
      }
    }
    
    return(uniqueStatesCounts)
  }



}

#' Fit and Compare Statistical Distributions for Lineage Waiting Times
#'
#' Fits a set of candidate probability distributions to lineage-specific waiting times
#' and computes AIC and BIC weights to assess relative model support.
#' Designed for parallel execution to improve speed.
fit_lineage_distributions <- function(waiting_times_list, 
                                      models = c(
                                        "exp",        # Simple, memoryless baseline
                                        "gamma",      # Flexible two-parameter generalization of exponential
                                        "weibull",    # Standard for increasing/decreasing hazard
                                        "lnorm",      # Captures skew from multiplicative processes
                                        "invgauss",   # Diffusion/first-passage-based timing
                                        "invweibull", # Heavy-tailed variant of Weibull
                                        "llogis",     # Non-monotonic hazard; rise then fall
                                        "pareto",     # Rare extreme waiting times; power-law tail
                                        "betapr"      # Flexible, heavy-tailed and positive-only
                                      ),
                                      num_cores = 2) {
  message("Fitting models for ", length(waiting_times_list), " lineages...")
  
  lineage_fits <- pbmclapply(seq_along(waiting_times_list), function(i) {
    wt <- waiting_times_list[[i]]$TimeToNext
    if (length(wt) < 2 || anyNA(wt)) return(list(models = NULL, AICw = NULL, BICw = NULL))
    
    res <- tryCatch({
      fits <- model_select(
        x = wt,
        models = models,
        type = "continuous",
        return = "all"
      )
      
      # Calculate AIC weights
      aic_vec <- setNames(fits$AIC, fits$ml)
      bic_vec <- setNames(fits$BIC, fits$ml)
      
      aic_weights <- geiger::aicw(aic_vec)
      bic_weights <- geiger::aicw(bic_vec)
      
      list(
        models = fits,
        AICw = aic_weights,
        BICw = bic_weights
      )
    }, error = function(e) {
      list(models = NULL, AICw = NULL, BICw = NULL)
    })
    
    return(res)
  }, mc.cores = num_cores)
  
  names(lineage_fits) <- paste0("Lineage_", seq_along(waiting_times_list))
  message("Done. Model fits completed.")
  return(lineage_fits)
}

#' Plot Model Weights per Lineage
#'
#' Creates a stacked bar plot of model weights (AIC or BIC) for each lineage
#' from fitted distribution results, and summarizes the best-performing models.
plot_model_weights_per_lineage <- function(fit_results,
                                           weight_type = c("aic", "bic"),
                                           candidate_models = NULL) {
  weight_type <- match.arg(weight_type)
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  extract_weights <- function(entry, lineage_name, candidate_models) {
    fitted_models <- entry$models$ml
    if (!all(candidate_models %in% fitted_models)) return(NULL)
    
    weight_matrix <- switch(
      weight_type,
      aic = entry$AICw,
      bic = entry$BICw
    )
    
    if (is.null(weight_matrix)) return(NULL)
    
    weight_matrix <- weight_matrix[rownames(weight_matrix) %in% candidate_models, , drop = FALSE]
    
    data.frame(
      Lineage = lineage_name,
      Model = rownames(weight_matrix),
      Weight = weight_matrix[, "w"],
      stringsAsFactors = FALSE
    )
  }
  
  if (is.null(candidate_models)) {
    stop("Please specify the candidate_models argument (vector of model names to require).")
  }
  
  lineage_names <- names(fit_results)
  total_lineages <- length(lineage_names)
  
  weights_list <- mapply(
    extract_weights,
    fit_results,
    lineage_names,
    MoreArgs = list(candidate_models = candidate_models),
    SIMPLIFY = FALSE
  )
  
  included <- sum(!sapply(weights_list, is.null))
  excluded <- total_lineages - included
  
  message(sprintf("Included %d lineages with complete fits. Excluded %d lineages missing at least one model.",
                  included, excluded))
  
  if (included == 0) {
    warning("No lineages met the inclusion criteria. Plot will not be generated.")
    return(invisible(NULL))
  }
  
  weights_df <- do.call(rbind, weights_list)
  
  ggplot(weights_df, aes(x = Lineage, y = Weight, fill = Model)) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("Model Weights per Lineage (", toupper(weight_type), ")", sep = ""),
      x = "Lineage",
      y = "Model Weight",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
  # Get the best model (highest weight) for each lineage
  best_models <- weights %>%
    group_by(Lineage) %>%
    slice_max(order_by = Weight, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Count how many times each model is best
  summary <- best_models %>%
    count(Model) %>%
    mutate(Percentage = 100 * n / sum(n))
  
  print(summary)
  return(weights_df)
}

#' Extract λ (Rate Parameter) from Exponential Fits
#'
#' Retrieves the estimated rate parameter (\eqn{\lambda}) from exponential
#' distribution fits for each lineage.
extract_exp_lambda <- function(results.fits) {
  lineage_names <- names(results.fits)
  
  lambda_vals <- sapply(lineage_names, function(ln) {
    model_list <- results.fits[[ln]]$models$univariateML
    if (!is.null(model_list$mlexp)) {
      return(model_list$mlexp[[1]])  # the rate parameter λ
    } else {
      return(NA)  # if exp wasn't fit or failed
    }
  })
  
  # Return as data.frame
  data.frame(Lineage = lineage_names, Lambda = lambda_vals)
}
