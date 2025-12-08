#' Phylogenetic rate–environment time-series correlation
#'
#' Computes bin-wise phylogenetically corrected mean nodal rate derivatives and correlates them with a binned environmental time series.
#'
#' @param tree A `phylo` (simmap) object with mapped regimes.
#' @param param_by_state Named numeric vector of per-regime rate scalars (e.g. BMM rate multipliers).
#' @param env_time Numeric vector of environmental ages (Ma before present, 0 = present).
#' @param env_value Numeric vector of environmental values (e.g. δ18O) matching `env_time`.
#' @param bin_width Bin width in Ma for discretizing time.
#' @param t_min,t_max Optional numeric limits (Ma before present) for the analysis time window.
#' @param verbose Logical; if TRUE, print diagnostics.
#' @param parallel Logical; if TRUE, build the tip × bin matrix in parallel using futures.
#' @param n_workers Integer; number of workers to use when `parallel = TRUE`.
#' @param show_progress Logical; if TRUE, show a progress bar while building the matrix.
#' @param center_param Logical; if TRUE (default), center and scale param_means to get param_std
#'   (z-score). If FALSE, only scale (no centering), so param_std = param_means / sd(param_means)
#'   and 0 corresponds to zero derivative.
#'
#' @return A list with bin midpoints (`t_mid`), parameter means (`param_means`), environmental means (`env_means`), standardized vectors (`param_std`, `env_std`), correlation (`correlation`), model summary (`lm_summary`), the tip × bin matrix (`X`), bin breaks (`breaks`), and the time window used.
#' @export
phylo_env_timeseries_corr_scalar <- function(tree,
                                             param_by_state,
                                             env_time,
                                             env_value,
                                             bin_width     = 2,
                                             t_min         = NULL,
                                             t_max         = NULL,
                                             verbose       = TRUE,
                                             parallel      = FALSE,
                                             n_workers     = NULL,
                                             show_progress = TRUE,
                                             center_param  = TRUE) {
  ## ---------- basic tree sizes ----------
  n_tips  <- Ntip(tree)
  n_nodes <- Nnode(tree)
  
  ## ---------- combined states (tips + nodes) ----------
  tip_states  <- phytools::getStates(tree, type = "tips")
  node_states <- phytools::getStates(tree, type = "nodes")
  
  names(tip_states)  <- as.character(1:n_tips)
  names(node_states) <- as.character((n_tips + 1):(n_tips + n_nodes))
  
  states_combined <- c(tip_states, node_states)
  
  ## ---------- node heights and ages before present ----------
  all_ids <- as.integer(names(states_combined))
  node_height <- sapply(all_ids, function(nd) phytools::nodeheight(tree, nd))
  names(node_height) <- names(states_combined)
  
  tree_height <- max(node_height, na.rm = TRUE)
  node_age    <- tree_height - node_height   # age before present (Ma)
  
  if (verbose) {
    cat("Tree height (Ma):", tree_height, "\n")
  }
  
  ## ---------- determine time window and bins ----------
  env_time_max <- max(env_time, na.rm = TRUE)
  overlap_max  <- min(tree_height, env_time_max)
  
  if (is.null(t_min)) {
    t_min_use <- 0
  } else {
    t_min_use <- max(0, t_min)
  }
  
  if (is.null(t_max)) {
    t_max_use <- overlap_max
  } else {
    t_max_use <- min(overlap_max, t_max)
  }
  
  if (t_min_use >= t_max_use) {
    stop("Invalid time window: t_min >= t_max after overlap adjustment.")
  }
  
  breaks <- seq(t_min_use, t_max_use, by = bin_width)
  if (tail(breaks, 1) < t_max_use) {
    breaks <- c(breaks, t_max_use)
  }
  n_int <- length(breaks) - 1
  if (n_int < 1) stop("'breaks' must have at least two elements.")
  
  if (verbose) {
    cat("Analysis time window (Ma): [", t_min_use, ", ", t_max_use, "]\n", sep = "")
    cat("Bin width (Ma):", bin_width, "\n")
    cat("Number of time bins:", n_int, "\n")
  }
  
  ## ---------- edge-level param changes (sanity check) ----------
  edges      <- tree$edge
  parent_ids <- edges[, 1]
  child_ids  <- edges[, 2]
  
  s_parent <- states_combined[as.character(parent_ids)]
  s_child  <- states_combined[as.character(child_ids)]
  
  has_param <- s_parent %in% names(param_by_state) & s_child %in% names(param_by_state)
  s_parent  <- s_parent[has_param]
  s_child   <- s_child[has_param]
  child_ids <- child_ids[has_param]
  
  if (length(s_parent) == 0) {
    warning("No edges where both parent and child states have parameters; returning NA.")
    t_mid <- (breaks[-length(breaks)] + breaks[-1]) / 2
    return(list(
      t_mid       = t_mid,
      param_means = rep(NA_real_, n_int),
      env_means   = rep(NA_real_, n_int),
      param_std   = NULL,
      env_std     = NULL,
      correlation = NA_real_,
      lm_summary  = NULL,
      X           = NULL,
      breaks      = breaks,
      time_window = c(t_min_use, t_max_use)
    ))
  }
  
  param_parent <- param_by_state[s_parent]
  param_child  <- param_by_state[s_child]
  
  change_idx_all  <- which(param_parent != param_child)
  num_changes_all <- length(change_idx_all)
  
  change_child_ids_all <- child_ids[change_idx_all]
  change_ages_all      <- node_age[as.character(change_child_ids_all)]
  
  in_window <- !is.na(change_ages_all) &
    change_ages_all >= t_min_use &
    change_ages_all <= t_max_use
  change_ages <- change_ages_all[in_window]
  num_changes <- length(change_ages)
  
  if (verbose) {
    cat("Total edges with parent!=child parameter (all times):", num_changes_all, "\n")
    cat("Edges with parameter changes within window:", num_changes, "\n")
    if (num_changes > 0) {
      cat("Change ages range within window (Ma before present):",
          range(change_ages, na.rm = TRUE), "\n")
      change_hist <- hist(change_ages, breaks = breaks, plot = FALSE)$counts
      cat("Number of change events per time bin (within window):\n")
      print(change_hist)
    }
  }
  
  if (num_changes == 0) {
    warning("No edges with parameter changes within the specified time window; returning NA.")
    t_mid <- (breaks[-length(breaks)] + breaks[-1]) / 2
    return(list(
      t_mid       = t_mid,
      param_means = rep(NA_real_, n_int),
      env_means   = rep(NA_real_, n_int),
      param_std   = NULL,
      env_std     = NULL,
      correlation = NA_real_,
      lm_summary  = NULL,
      X           = NULL,
      breaks      = breaks,
      time_window = c(t_min_use, t_max_use)
    ))
  }
  
  ## ---------- paths from root to each tip ----------
  paths <- ape::nodepath(tree)
  
  ## ---------- helper to compute one X row (tip x bins) ----------
  compute_X_row <- function(i, paths, states_combined, breaks, param_by_state, node_age) {
    path_nodes <- paths[[i]]
    n_int      <- length(breaks) - 1
    row_vec    <- numeric(n_int)
    
    if (length(path_nodes) < 2) return(row_vec)
    
    # root -> ... -> tip; node_age decreases along the path
    for (j in seq_len(n_int)) {
      t0 <- breaks[j]
      t1 <- breaks[j + 1]
      
      for (k in 2:(length(path_nodes) - 1)) {
        node_k <- path_nodes[k]  # child (younger)
        age_k  <- node_age[as.character(node_k)]
        
        if (!is.na(age_k) && age_k < t1 && age_k >= t0) {
          node_parent <- path_nodes[k - 1]  # parent (older)
          
          s_k    <- states_combined[as.character(node_k)]
          s_prev <- states_combined[as.character(node_parent)]
          
          if (is.na(s_k) || is.na(s_prev)) next
          if (!s_k    %in% names(param_by_state)) next
          if (!s_prev %in% names(param_by_state)) next
          
          val_k    <- param_by_state[s_k]
          val_prev <- param_by_state[s_prev]
          
          # Only consider actual parameter changes:
          if (val_k == val_prev) next
          
          age_prev <- node_age[as.character(node_parent)]
          # parent must be older (larger age) than child
          if (is.na(age_prev) || age_prev <= age_k) next
          
          # change per unit forward time (parent -> child):
          # ages are "before present", so dt = age_parent - age_child > 0
          deriv <- (val_k - val_prev) / (age_prev - age_k)
          
          # Accumulate (sum) derivatives within this bin
          row_vec[j] <- row_vec[j] + deriv
        }
      }
    }
    row_vec
  }
  
  ## ---------- optional: set up parallel + progress infrastructure ----------
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE) ||
        !requireNamespace("progressr", quietly = TRUE)) {
      stop("Packages 'future', 'future.apply', and 'progressr' are required for parallel = TRUE.")
    }
    # Set future plan based on n_workers
    if (!is.null(n_workers)) {
      future::plan(future::multisession, workers = n_workers)
    } else {
      # default plan if none set; this is conservative
      future::plan(future::multisession)
    }
    
    if (show_progress) {
      progressr::handlers(global = TRUE)
      progressr::handlers("txtprogressbar")
    }
  }
  
  ## ---------- build X (tip x time-bin derivatives), with optional parallel + progress ----------
  if (!parallel) {
    # sequential
    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = n_tips, style = 3)
      X_list <- vector("list", n_tips)
      for (i in seq_len(n_tips)) {
        X_list[[i]] <- compute_X_row(
          i,
          paths           = paths,
          states_combined = states_combined,
          breaks          = breaks,
          param_by_state  = param_by_state,
          node_age        = node_age
        )
        setTxtProgressBar(pb, i)
      }
      close(pb)
      X <- do.call(rbind, X_list)
    } else {
      X <- t(sapply(
        X   = seq_len(n_tips),
        FUN = compute_X_row,
        paths           = paths,
        states_combined = states_combined,
        breaks          = breaks,
        param_by_state  = param_by_state,
        node_age        = node_age
      ))
    }
  } else {
    # parallel with future.apply
    if (show_progress) {
      X <- progressr::with_progress({
        p <- progressr::progressor(along = seq_len(n_tips))
        t(future.apply::future_sapply(
          X   = seq_len(n_tips),
          FUN = function(i, paths, states_combined, breaks, param_by_state, node_age) {
            p()
            compute_X_row(
              i,
              paths           = paths,
              states_combined = states_combined,
              breaks          = breaks,
              param_by_state  = param_by_state,
              node_age        = node_age
            )
          },
          paths           = paths,
          states_combined = states_combined,
          breaks          = breaks,
          param_by_state  = param_by_state,
          node_age        = node_age,
          future.seed     = TRUE
        ))
      })
    } else {
      X <- t(future.apply::future_sapply(
        X   = seq_len(n_tips),
        FUN = compute_X_row,
        paths           = paths,
        states_combined = states_combined,
        breaks          = breaks,
        param_by_state  = param_by_state,
        node_age        = node_age,
        future.seed     = TRUE
      ))
    }
  }
  
  rownames(X) <- tree$tip.label
  
  if (verbose) {
    cat("Non-zero derivatives per bin (colSums(X != 0)):\n")
    print(colSums(X != 0))
  }
  
  ## ---------- phylogenetically corrected mean per bin ----------
  param_means <- rep(NA_real_, n_int)
  names(param_means) <- paste0("bin_", seq_len(n_int))
  
  cols_with_data <- which(colSums(X != 0) > 0)
  if (length(cols_with_data) == 0) {
    warning("No time intervals with non-zero derivatives across tip paths; returning NA correlation.")
    t_mid <- (breaks[-length(breaks)] + breaks[-1]) / 2
    return(list(
      t_mid       = t_mid,
      param_means = param_means,
      env_means   = rep(NA_real_, n_int),
      param_std   = NULL,
      env_std     = NULL,
      correlation = NA_real_,
      lm_summary  = NULL,
      X           = X,
      breaks      = breaks,
      time_window = c(t_min_use, t_max_use)
    ))
  }
  
  for (j in cols_with_data) {
    xj   <- X[, j]
    taxa <- which(xj != 0)
    if (length(taxa) < 2) next
    
    ones <- matrix(1, length(taxa), 1)
    C_little <- ape::vcv.phylo(ape::keep.tip(tree, tree$tip.label[taxa]))
    
    C_inv <- tryCatch(solve(C_little), error = function(e) NULL)
    if (is.null(C_inv)) next
    
    num <- t(ones) %*% C_inv %*% xj[taxa]
    den <- sum(C_inv)
    param_means[j] <- as.numeric(num / den)
  }
  
  if (verbose) {
    cat("phylo-corrected param_means:\n")
    print(param_means)
  }
  
  ## ---------- bin environmental series ----------
  env_means <- rep(NA_real_, n_int)
  for (j in seq_len(n_int)) {
    t0 <- breaks[j]
    t1 <- breaks[j + 1]
    in_bin <- env_time >= t0 & env_time < t1
    if (any(in_bin)) {
      env_means[j] <- mean(env_value[in_bin], na.rm = TRUE)
    }
  }
  
  if (verbose) {
    cat("env_means per bin:\n")
    print(env_means)
  }
  
  ## ---------- correlate ----------
  t_mid <- (breaks[-length(breaks)] + breaks[-1]) / 2
  
  valid <- which(!is.na(param_means) & !is.na(env_means))
  if (length(valid) < 3) {
    warning("Fewer than 3 intervals with non-NA parameter and environmental means; returning NA correlation.")
    return(list(
      t_mid       = t_mid,
      param_means = param_means,
      env_means   = env_means,
      param_std   = NULL,
      env_std     = NULL,
      correlation = NA_real_,
      lm_summary  = NULL,
      X           = X,
      breaks      = breaks,
      time_window = c(t_min_use, t_max_use)
    ))
  }
  
  # y: either centered+scaled (z-score) or scaled-only (no centering),
  # depending on center_param flag
  param_std <- scale(param_means[valid], center = center_param, scale = TRUE)[, 1]
  # x: still z-scored (centered+scaled)
  env_std   <- scale(env_means[valid],   center = TRUE,        scale = TRUE)[, 1]
  
  r      <- cor(param_std, env_std, use = "complete.obs")
  lm_fit <- lm(param_std ~ env_std)
  
  list(
    t_mid       = t_mid,
    param_means = param_means,
    env_means   = env_means,
    param_std   = param_std,
    env_std     = env_std,
    correlation = r,
    lm_summary  = summary(lm_fit),
    lm          = lm_fit,
    X           = X,
    breaks      = breaks,
    time_window = c(t_min_use, t_max_use)
  )
}


#' Plot phylogenetic rate–environment time-series results
#'
#' Produces a 2×2-style diagnostic plot from `phylo_env_timeseries_corr_scalar()` output,
#' including a nodal-derivative heatmap, binned environmental series, and a standardized
#' rate–environment correlation panel with optional annotations.
#'
#' @param res List returned by `phylo_env_timeseries_corr_scalar()`.
#' @param env_time Numeric vector of environmental ages (Ma before present).
#' @param env_value Numeric vector of environmental values (e.g. δ18O) matching `env_time`.
#' @param max_tips_to_plot Maximum number of tips to display in the heatmap (for readability).
#' @param main_prefix Optional character prefix for panel titles.
#' @param reverse_time Logical; if TRUE, plot time with older ages on the left and 0 Ma on the right.
#' @param flip_env_y Logical; if TRUE, flip the δ18O y-axis visually (higher δ18O plotted lower),
#'   matching common paleo-style plots. Underlying values are unchanged.
#' @param meta_text Optional character string to override the default annotation of bin width and
#'   time interval in the environmental panel.
#' @param annotate_axes Logical; if TRUE, add user-specified text labels near the x and y axes
#'   of the correlation panel to clarify the meaning of positive/negative values.
#' @param x_annot_left Optional character string placed near the left side of the correlation
#'   x-axis (e.g. "warmer bins (lower δ18O)").
#' @param x_annot_right Optional character string placed near the right side of the correlation
#'   x-axis (e.g. "colder bins (higher δ18O)").
#' @param y_annot_bottom Optional character string placed near the bottom of the correlation
#'   y-axis (e.g. "net rate decreases").
#' @param y_annot_top Optional character string placed near the top of the correlation
#'   y-axis (e.g. "net rate increases").
#' @param shade_quadrants Logical; if TRUE, lightly shade the warm–rate-down (x < 0, y < 0)
#'   and cold–rate-up (x > 0, y > 0) quadrants in the correlation panel.
#' @param quad_col_warm_down Color used to shade the warm–rate-down quadrant.
#' @param quad_col_cold_up Color used to shade the cold–rate-up quadrant.
#' @param quad_alpha Numeric alpha (0–1) controlling the transparency of quadrant shading.
#'
#' @return Invisibly returns NULL; called for its plotting side effects.
#' @export
plot_phylo_env_timeseries <- function(res,
                                      env_time,
                                      env_value,
                                      max_tips_to_plot = 200,
                                      main_prefix = NULL,
                                      reverse_time = FALSE,
                                      flip_env_y = FALSE,
                                      meta_text = NULL,
                                      annotate_axes = FALSE,
                                      x_annot_left = NULL,
                                      x_annot_right = NULL,
                                      y_annot_bottom = NULL,
                                      y_annot_top = NULL,
                                      shade_quadrants = FALSE,
                                      quad_col_warm_down = "lightblue",
                                      quad_col_cold_up  = "lightcoral",
                                      quad_alpha        = 0.3) {
  if (is.null(res$t_mid)) {
    stop("Result object must have t_mid (bin midpoints). Did you pass the output of phylo_env_timeseries_corr_scalar()?")
  }
  
  t_mid       <- res$t_mid
  env_means   <- res$env_means
  param_std   <- res$param_std
  env_std     <- res$env_std
  X           <- res$X
  r           <- res$correlation
  lm_sum      <- res$lm_summary
  
  # Title prefix
  if (is.null(main_prefix)) {
    main_prefix <- ""
  } else if (nzchar(main_prefix)) {
    main_prefix <- paste0(main_prefix, " - ")
  }
  
  # Time limits (for whole plotting window)
  xlim_time <- range(t_mid, na.rm = TRUE)
  if (reverse_time) {
    xlim_time <- rev(xlim_time)  # visually reverse time axis
  }
  
  # Env limits
  env_range <- range(c(env_means, env_value), na.rm = TRUE)
  ylim_env  <- env_range
  if (flip_env_y) {
    ylim_env <- rev(env_range)  # visually flip δ18O axis
  }
  
  # Y-axis label for env panel
  env_ylab <- expression(paste("Binned ", delta^{18}, "O"))
  if (flip_env_y) {
    env_ylab <- expression(paste("Binned ", delta^{18}, "O (colder down)"))
  }
  
  # Meta annotation text (bin width + time window) - always shown
  if (!is.null(meta_text)) {
    meta_str <- meta_text
  } else {
    # Infer bin width from breaks if available
    bw <- NA_real_
    if (!is.null(res$breaks)) {
      diffs <- diff(res$breaks)
      if (length(diffs) > 0) {
        # use median bin size to ignore a slightly shorter last bin
        bw <- stats::median(diffs, na.rm = TRUE)
      }
    }
    # Time window from result if present, otherwise from t_mid
    if (!is.null(res$time_window)) {
      tw <- res$time_window
    } else {
      tw <- range(t_mid, na.rm = TRUE)
    }
    meta_str <- sprintf("Bin width ≈ %.1f Ma\nInterval = %.1f–%.1f Ma",
                        bw, tw[1], tw[2])
  }
  
  # 2x2 layout where the right column spans both rows
  layout(matrix(c(1, 3,
                  2, 3), nrow = 2, byrow = TRUE))
  
  # layout(
  #   matrix(c(1, 0,
  #            1, 3,
  #            2, 3,
  #            2, 0),
  #          nrow = 4, byrow = TRUE),
  #   widths  = c(1, 1),   # left column wider than right
  #   heights = c(1, 1, 1, 1)
  # )
  
  ## ------------------------------------------------------------------------
  ## (1) Heatmap of absolute nodal derivatives (tip x time-bin)
  ## ------------------------------------------------------------------------
  if (!is.null(X)) {
    X_abs  <- abs(X)
    n_tips <- nrow(X_abs)
    
    # Subsample tips for plotting if needed
    if (n_tips > max_tips_to_plot) {
      tip_idx <- sort(unique(round(seq(1, n_tips, length.out = max_tips_to_plot))))
    } else {
      tip_idx <- 1:n_tips
    }
    
    z <- t(X_abs[tip_idx, , drop = FALSE])  # bins x tips
    
    image(x = t_mid,
          y = seq_along(tip_idx),
          z = z,
          xlab = "Time (Ma before present)",
          ylab = "Tips (subset)",
          xlim = xlim_time,
          main = paste0(main_prefix, "Abs. nodal derivatives"))
  } else {
    plot.new()
    text(0.5, 0.5, "No X matrix in result", cex = 1.2)
  }
  
  ## ------------------------------------------------------------------------
  ## (2) Binned environmental series + raw curve
  ## ------------------------------------------------------------------------
  plot(t_mid, env_means,
       type = "b",
       xlab = "Time (Ma before present)",
       ylab = env_ylab,
       xlim = xlim_time,
       ylim = ylim_env,
       main = paste0(main_prefix, "Binned environmental series"))
  
  # Add raw env in background
  ord <- order(env_time)
  lines(env_time[ord], env_value[ord],
        lty = 3)
  
  legend("topright",
         legend = c("Binned δ18O", "Raw δ18O"),
         lty    = c(1, 3),
         pch    = c(1, NA),
         bty    = "n")
  
  # Add meta annotation (place inside plot, top-left, nudged in)
  usr <- par("usr")  # c(xmin, xmax, ymin, ymax) for current panel
  x_meta <- usr[1] + 0.02 * (usr[2] - usr[1])  # a bit right of left edge
  y_meta <- usr[4] - 0.02 * (usr[4] - usr[3])  # a bit below top edge
  text(x   = x_meta,
       y   = y_meta,
       labels = meta_str,
       adj   = c(0, 1),  # left-aligned, top-aligned
       cex   = 0.8)
  
  ## ------------------------------------------------------------------------
  ## (3) Variance-normalized means vs variance-normalized env series
  ## ------------------------------------------------------------------------
  if (!is.null(param_std) && !is.null(env_std)) {
    # Basic scatter first
    plot(env_std, param_std,
         xlab = "Environmental series (z-score)",
         ylab = "Phylo-corrected mean derivative (scaled)",
         main = paste0(
           main_prefix,
           "Rate–environment relationship (r = ",
           ifelse(is.na(r), "NA", round(r, 2)),
           ")"
         ),
         pch  = 21, col='black', bg='grey', cex = 2)
    
    usr_cor <- par("usr")  # c(xmin, xmax, ymin, ymax)
    
    # Optional quadrant shading (drawn behind lines/points conceptually,
    # but here with alpha so points/lines remain visible)
    if (shade_quadrants) {
      # Warm + rate-down quadrant (x < 0, y < 0)
      if (usr_cor[1] < 0 && usr_cor[3] < 0) {
        rect(xleft   = usr_cor[1],
             ybottom = usr_cor[3],
             xright  = 0,
             ytop    = 0,
             border  = NA,
             col     = grDevices::adjustcolor(quad_col_warm_down, alpha.f = quad_alpha))
      }
      # Cold + rate-up quadrant (x > 0, y > 0)
      if (usr_cor[2] > 0 && usr_cor[4] > 0) {
        rect(xleft   = 0,
             ybottom = 0,
             xright  = usr_cor[2],
             ytop    = usr_cor[4],
             border  = NA,
             col     = grDevices::adjustcolor(quad_col_cold_up, alpha.f = quad_alpha))
      }
    }
    
    # Re-draw points on top (in case shading obscured them slightly)
    points(env_std, param_std, pch  = 21, col='black', bg='grey', cex = 2)
    
    # Zero lines to show the quadrants explicitly
    abline(h = 0, v = 0, lty = 3, col = "grey60")
    
    # Regression line
    abline(lm(param_std ~ env_std), lwd = 2)
    # 1:1 line for reference
    abline(0, 1, lty = 2, col = "grey50")
    
    # Add regression stats if available, bottom-right, nudged in
    if (!is.null(lm_sum)) {
      slope <- lm_sum$coefficients["env_std", "Estimate"]
      p_val <- lm_sum$coefficients["env_std", "Pr(>|t|)"]
      R2    <- lm_sum$r.squared
      
      stats_str <- sprintf("slope = %.2f\nR² = %.2f\np = %.3g",
                           slope, R2, p_val)
      
      usr_cor <- par("usr")  # update in case abline expanded limits
      x_pos <- usr_cor[2] - 0.02 * (usr_cor[2] - usr_cor[1])  # a bit left of right edge
      y_pos <- usr_cor[3] + 0.02 * (usr_cor[4] - usr_cor[3])  # a bit above bottom
      
      text(x   = x_pos,
           y   = y_pos,
           labels = stats_str,
           adj   = c(1, 0),  # right-aligned, bottom-aligned
           cex   = 1.5)
    }
    
    # Extra axis-meaning annotations, if requested
    if (annotate_axes) {
      ## X-axis annotations: keep using mtext (already horizontal and well-placed)
      if (!is.null(x_annot_left)) {
        mtext(x_annot_left, side = 1, line = 2, adj = 0, cex = 0.8)
      }
      if (!is.null(x_annot_right)) {
        mtext(x_annot_right, side = 1, line = 2, adj = 1, cex = 0.8)
      }
      
      ## Y-axis annotations: draw horizontal labels just left of the axis
      usr_cor <- par("usr")  # c(xmin, xmax, ymin, ymax)
      dx <- usr_cor[2] - usr_cor[1]
      dy <- usr_cor[4] - usr_cor[3]
      
      # allow drawing outside the plot region
      op <- par(xpd = NA)
      on.exit(par(op), add = TRUE)
      
      # x position a bit left of the y-axis (closer to axis: 0.03 * dx)
      x_lab <- usr_cor[1] - 0.03 * dx
      
      # bottom label near lower part of axis
      if (!is.null(y_annot_bottom)) {
        y_bot <- usr_cor[3] + 0.15 * dy
        text(x   = x_lab,
             y   = y_bot,
             labels = y_annot_bottom,
             adj   = c(1, 0.5),  # right-aligned, centered vertically
             cex   = 0.9)        # slightly larger
      }
      
      # top label near upper part of axis
      if (!is.null(y_annot_top)) {
        y_top <- usr_cor[4] - 0.15 * dy
        text(x   = x_lab,
             y   = y_top,
             labels = y_annot_top,
             adj   = c(1, 0.5),
             cex   = 0.9)
      }
    }
    
    legend("topleft",
           inset  = c(0.01, 0.01),         # tiny nudge inwards from the corner
           legend = c("Bins", "Regression", "1:1 line"),
           pch    = c(16, NA, NA),
           lty    = c(NA, 1, 2),
           col    = c("black", "black", "grey50"),
           bty    = "n",
           cex    = 1.2)                   # increase legend text/points size
  } else {
    plot.new()
    text(0.5, 0.5, "Too few bins for standardized correlation", cex = 1.2)
  }
  
  invisible(NULL)
}


## ============================================================================
##  Setup
## ============================================================================

library(ape)
library(phytools)
library(future)
library(future.apply)
library(progressr)

## NOTE: phylo_env_timeseries_corr_scalar() and plot_phylo_env_timeseries()
## are assumed to be defined in your environment.


## ============================================================================
##  Data extraction
## ============================================================================

# mvgls model and tree
model <- runs$min10.ic20.gic$model_no_uncertainty
tree  <- runs$min10.ic20.gic$tree_no_uncertainty_untransformed
param_by_state <- model$param

# Westerhold environmental series
env_time  <- westerhold.long$age_tuned                     # 0 = present, increasing into past
env_value <- westerhold.long$ISOBENd18oLOESSsmoothLongTerm

# Quick checks (optional)
cat("env_time range (Ma before present):\n")
print(range(env_time, na.rm = TRUE))
cat("env_value summary:\n")
print(summary(env_value))


## ============================================================================
##  Analysis parameters
## ============================================================================

# bin_width      <- 4      # 3-Ma bins
# t_min          <- 0      # start at present
# t_max          <- NULL   # full overlap (no restriction)
# use_parallel   <- TRUE
# n_workers      <- 6
# show_progress  <- TRUE

bin_width      <- 4      # 3-Ma bins
t_min          <- 0      # start at present
t_max          <- NULL   # full overlap (no restriction)
use_parallel   <- TRUE
n_workers      <- 6
show_progress  <- TRUE
center_param   <- FALSE

## ============================================================================
##  Run analysis: full overlap window
## ============================================================================

cat("\n================ FULL OVERLAP WINDOW ================\n")

res_full <- phylo_env_timeseries_corr_scalar(
  tree           = tree,
  param_by_state = param_by_state,
  env_time       = env_time,
  env_value      = env_value,
  bin_width      = bin_width,
  t_min          = t_min,      # 0 Ma
  t_max          = t_max,      # NULL => up to overlap_max
  verbose        = TRUE,
  parallel       = use_parallel,
  n_workers      = n_workers,
  show_progress  = show_progress, 
  center_param   = center_param
)

cat("\n===== RESULT SUMMARY (full window) =====\n")
cat("Time window (Ma):", paste(res_full$time_window, collapse = " – "), "\n")
cat("Correlation (param_std vs env_std):", res_full$correlation, "\n\n")
print(res_full$lm_summary)


## ============================================================================
##  Plotting: full window
## ============================================================================

cat("\nPlotting full-window results...\n")


quartz(file='timeseries.pdf', height=8, width=10, type='pdf')
plot_phylo_env_timeseries(
  res       = res_full,
  env_time  = env_time,
  env_value = env_value,
  max_tips_to_plot = 2000,
  main_prefix      = "Full window",
  reverse_time     = TRUE,
  flip_env_y       = TRUE,
  annotate_axes    = TRUE,
  x_annot_left     = expression(
    atop(delta^18*O < bar(delta^18*O),
         "Warmer")
  ),
  x_annot_right    = expression(
    atop(delta^18*O > bar(delta^18*O),
         "Colder")
  ),
  y_annot_bottom   = "Rate\ndecrease ↓",
  y_annot_top      = "Rate\nincrease ↑",
  shade_quadrants  = TRUE,
  quad_col_warm_down = "lightblue",
  quad_col_cold_up  = "lightcoral",
  quad_alpha        = 0.3
)
dev.off()


## ============================================================================
##  (Optional) Inspect key outputs
## ============================================================================

cat("\nBin midpoints (t_mid):\n")
print(res_full$t_mid)

cat("\nparam_means:\n")
print(res_full$param_means)

cat("\nenv_means:\n")
print(res_full$env_means)

# res_full$X contains the full tip × bin derivative matrix.

