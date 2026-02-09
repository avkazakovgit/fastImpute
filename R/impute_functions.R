#' Fast Kalman Imputation
#'
#' Imputes missing values using a fast C++ Kalman Smoother with robust parameter estimation.
#'
#' @param DT A data.table containing the panel data.
#' @param vars Character vector of variable names to impute.
#' @param id_var Character string of the ID variable name.
#' @param time_var Character string of the time variable name.
#' @param params_sample Numeric (0-1). Fraction of IDs to use for parameter estimation.
#' @param max_n_ids Integer. Max IDs to use for parameter estimation.
#' @param out_suff Character string or NULL. Suffix for new columns. If NULL, overwrites in place.
#' @param max_hole Integer or Inf. Max length of internal NA gap to fill.
#' @param max_endpoint Integer or Inf. Max distance from data to fill at endpoints.
#' @param degree Integer (1 or 2). 1 = Level (Random Walk), 2 = Trend (Local Linear).
#' @param inplace Logical. If TRUE, modifies DT by reference.
#' @param verbose Logical. Print progress.
#'
#' @import data.table
#' @import KFAS
#' @export
kalman_fast <- function(
    DT, vars, id_var, time_var,
    params_sample = 0.1,
    max_n_ids = 1000,
    out_suff = "_kal",     # NULL = Overwrite, String = New Cols
    max_hole = Inf,        # Max length of internal NA gap to fill
    max_endpoint = Inf,    # Max distance from data to fill at starts/ends
    degree = 2L,
    inplace = TRUE,
    verbose = TRUE
) {
  # Note: require() is replaced by @import in packages
  # We assume inputs are checked by the package environment

  if (!data.table::is.data.table(DT)) {
    DT <- data.table::as.data.table(DT)
    inplace <- TRUE
  }
  if (!inplace) DT <- data.table::copy(DT)

  # Ensure sorted for time-series validity
  data.table::setorderv(DT, c(id_var, time_var))

  # --- Helper: Mask Calculation for Holes/Endpoints ---
  calc_fill_mask <- function(x, mh, me) {
    is_na <- is.na(x)
    if (!any(is_na)) return(rep(FALSE, length(x)))

    n <- length(x)
    idx_obs <- which(!is_na)
    if (length(idx_obs) == 0) return(rep(FALSE, n)) # All NA case

    first_obs <- idx_obs[1]
    last_obs  <- idx_obs[length(idx_obs)]
    mask <- rep(FALSE, n)

    # 1. Internal Holes
    if (first_obs < last_obs - 1) {
      r <- rle(is_na)
      # Flag runs that are NA (TRUE) and short enough
      ok_runs <- r$values & (r$lengths <= mh)
      internal_mask <- rep(ok_runs, r$lengths)

      # Zero out endpoints from this logic (they obey 'me', not 'mh')
      if (first_obs > 1) internal_mask[1:(first_obs-1)] <- FALSE
      if (last_obs < n)  internal_mask[(last_obs+1):n]  <- FALSE

      mask <- mask | internal_mask
    }

    # 2. Endpoints
    if (is.infinite(me)) {
      if (first_obs > 1) mask[1:(first_obs-1)] <- TRUE
      if (last_obs < n)  mask[(last_obs+1):n]  <- TRUE
    } else if (me > 0) {
      if (first_obs > 1) {
        start <- max(1, first_obs - me)
        mask[start:(first_obs-1)] <- TRUE
      }
      if (last_obs < n) {
        end <- min(n, last_obs + me)
        mask[(last_obs+1):end] <- TRUE
      }
    }
    return(mask)
  }

  # --- 1. Helper: Transformation Logic ---
  apply_transform <- function(x) {
    min_val <- min(x, na.rm = TRUE)
    if (min_val > 0) return(list(vals = log(x), type = "log", offset = 0))
    return(list(vals = asinh(x), type = "asinh", offset = 0))
  }

  # --- 2. Robust Parameter Estimator ---
  estimate_params_robust <- function(dt_source, var_name) {
    u_ids <- unique(dt_source[[id_var]])
    n_total <- length(u_ids)
    n_target <- max(50, min(ceiling(n_total * params_sample), max_n_ids, n_total))
    sel_ids <- if (n_target >= n_total) u_ids else sample(u_ids, n_target)

    if(verbose) message(sprintf("  Estimating params for '%s' using %d IDs...", var_name, length(sel_ids)))

    # Extract & Transform
    sub_dt <- dt_source[dt_source[[id_var]] %in% sel_ids, .SD, .SDcols = c(id_var, var_name)]
    raw_vals <- as.numeric(sub_dt[[var_name]])

    tf <- apply_transform(raw_vals)
    sub_dt[, val_tf := tf$vals]

    # Scale per ID (Standardize)
    sub_dt[, `:=`(mu = mean(val_tf, na.rm=TRUE), sig = sd(val_tf, na.rm=TRUE)), by = c(id_var)]
    sub_dt[is.na(sig) | sig < 1e-8, sig := 1]
    sub_dt[, y_scaled := (val_tf - mu) / sig]

    # Format for KFAS
    y_list <- split(sub_dt$y_scaled, sub_dt[[id_var]])
    y_long <- unlist(lapply(y_list, function(x) c(x, rep(NA_real_, 5))))

    if (sum(!is.na(y_long)) < 10) return(NULL)

    tryCatch({
      Q_init <- replicate(degree, matrix(NA_real_), simplify = FALSE)
      mod <- KFAS::SSModel(y ~ KFAS::SSMtrend(degree = degree, Q = Q_init), data = data.frame(y = y_long), H = matrix(NA_real_))
      fit <- KFAS::fitSSM(mod, inits = rep(log(0.1), degree + 1), method = "BFGS", control=list(maxit=50, trace=0))

      h <- fit$model$H[1,1,1]
      q <- numeric(degree)
      q[1] <- fit$model$Q[1,1,1]
      if(degree >= 2) q[2] <- fit$model$Q[2,2,1]

      return(list(H = h, Q = q, init_var = h*10, type = tf$type, valid = TRUE))
    }, error = function(e) {
      return(list(H = 0.1, Q = rep(0.01, degree), init_var = 1.0, type = tf$type, valid = FALSE))
    })
  }

  # --- 3. Main Loop ---
  for (v in vars) {
    # 3.1 Estimate Global Params
    pars <- estimate_params_robust(DT, v)
    if (is.null(pars)) next

    H_val <- pars$H
    Q_val <- pars$Q
    init_P <- pars$init_var
    tf_type <- pars$type

    if (verbose) message(sprintf("  Processing %s | Type: %s | H: %.4f | Q: %.4f", v, tf_type, H_val, Q_val[1]))

    # 3.2 Temp Columns for Calculation
    temp_tf_col <- paste0("temp_tf_", v)

    raw_v <- as.numeric(DT[[v]])
    if (tf_type == "log") {
      DT[, (temp_tf_col) := log(raw_v)]
    } else {
      DT[, (temp_tf_col) := asinh(raw_v)]
    }

    # 3.3 Apply Smoother
    smooth_col_name <- paste0("temp_smooth_", v)

    DT[, (smooth_col_name) := {
      val <- .SD[[1]]
      mu <- mean(val, na.rm=TRUE)
      sig <- sd(val, na.rm=TRUE)
      if (is.na(sig) || sig < 1e-8) sig <- 1
      y_scaled <- (val - mu) / sig

      if (all(is.na(y_scaled))) {
        as.numeric(NA)
      } else {
        first_obs <- na.omit(y_scaled)[1]
        # CALL C++ FUNCTION HERE
        res_scaled <- kalman_smoother_cpp(y_scaled, degree, H_val, Q_val, first_obs, init_P)
        res_tf <- (res_scaled * sig) + mu

        # Reverse Transform
        res_final <- if (tf_type == "log") exp(res_tf) else sinh(res_tf)

        # --- NEW: APPLY HOLE/ENDPOINT FILTER ---
        if (!is.infinite(max_hole) || !is.infinite(max_endpoint)) {
          # Use the RAW transformed data 'val' (which has NAs) to detect holes
          mask <- calc_fill_mask(val, max_hole, max_endpoint)
          # Set imputed values to NA where mask is FALSE
          res_final[!mask] <- NA
        }
        res_final
      }
    }, by = id_var, .SDcols = temp_tf_col]

    # 3.4 SAFE MERGE (The "Imputation" Step)
    if (is.null(out_suff)) {
      # Overwrite in-place
      DT[, (v) := data.table::fcoalesce(as.numeric(get(v)), get(smooth_col_name))]
    } else {
      # Create new column with suffix
      out_col <- paste0(v, out_suff)
      DT[, (out_col) := data.table::fcoalesce(as.numeric(get(v)), get(smooth_col_name))]
    }

    # 3.5 Cleanup Temp Columns
    DT[, c(temp_tf_col, smooth_col_name) := NULL]
  }

  DT
}


#' Fast Factor-Based Imputation (Robust)
#'
#' Imputes missing values using an EM-SVD approach.
#' Automatically applies log/asinh transformation and Double Centering (Median-based)
#' to prevent outliers from corrupting the imputation.
#'
#' @export
fbi_fast <- function(
    DT,
    vars,
    id_var,
    time_var,
    rank = 4L,
    thresh = 1e-5,
    max_iter = 20L,
    out_suff = "_fbi",      # NULL = Overwrite
    max_hole = Inf,
    max_endpoint = Inf,
    verbose = TRUE
) {
  if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)
  data.table::setorderv(DT, c(id_var, time_var))

  # --- Helper: Mask Calculation ---
  calc_fill_mask <- function(x, mh, me) {
    is_na <- is.na(x)
    if (!any(is_na)) return(rep(FALSE, length(x)))
    n <- length(x)
    idx_obs <- which(!is_na)
    if (length(idx_obs) == 0) return(rep(FALSE, n))
    first_obs <- idx_obs[1]; last_obs <- idx_obs[length(idx_obs)]
    mask <- rep(FALSE, n)
    if (first_obs < last_obs - 1) {
      r <- rle(is_na)
      mask <- mask | rep(r$values & (r$lengths <= mh), r$lengths)
      if (first_obs > 1) mask[1:(first_obs-1)] <- FALSE
      if (last_obs < n)  mask[(last_obs+1):n]  <- FALSE
    }
    if (is.infinite(me)) {
      if (first_obs > 1) mask[1:(first_obs-1)] <- TRUE
      if (last_obs < n)  mask[(last_obs+1):n]  <- TRUE
    } else if (me > 0) {
      if (first_obs > 1) mask[max(1, first_obs - me):(first_obs-1)] <- TRUE
      if (last_obs < n)  mask[(last_obs+1):min(n, last_obs + me)] <- TRUE
    }
    return(mask)
  }

  for (v in vars) {
    if (verbose) message(sprintf("Processing '%s' (Rank: %d)...", v, rank))

    # --- STEP 0: CHECK & TRANSFORM ---
    raw_vals_global <- as.numeric(DT[[v]])
    min_val  <- min(raw_vals_global, na.rm = TRUE)
    use_log <- (min_val > 0)

    if (use_log) {
      DT[, temp_v_trans := log(as.numeric(get(v)))]
    } else {
      DT[, temp_v_trans := asinh(as.numeric(get(v)))]
    }

    # --- STEP A: Pivot ---
    form <- as.formula(paste(time_var, "~", id_var))
    wide_dt <- data.table::dcast(DT, form, value.var = "temp_v_trans")

    time_index <- wide_dt[[time_var]]
    X <- as.matrix(wide_dt[, -1])

    T_dim <- nrow(X)
    N_dim <- ncol(X)
    curr_rank <- min(rank, T_dim - 1, N_dim - 1)

    if (curr_rank < 1) {
      DT[, temp_v_trans := NULL]
      if(verbose) message("  -> Skipped: Dimensions too small.")
      next
    }

    # --- STEP B: PRE-PROCESSING (Double Centering) ---
    mu_unit <- colMeans(X, na.rm = TRUE)
    sd_unit <- apply(X, 2, sd, na.rm = TRUE)
    sd_unit[is.na(sd_unit) | sd_unit < 1e-8] <- 1
    mu_unit[is.na(mu_unit)] <- 0

    X_std <- t((t(X) - mu_unit) / sd_unit)

    # Use MEDIAN for time centering to kill outliers
    mu_time <- apply(X_std, 1, median, na.rm = TRUE)
    mu_time[is.na(mu_time)] <- 0
    X_cent <- X_std - mu_time

    # --- STEP C: EM-SVD ALGORITHM ---
    na_idx <- which(is.na(X_cent))

    if (length(na_idx) > 0) {
      X_filled <- X_cent
      X_filled[na_idx] <- 0

      for (iter in 1:max_iter) {
        svd_fit <- try(svd(X_filled, nu = curr_rank, nv = curr_rank), silent=TRUE)
        if (inherits(svd_fit, "try-error")) break

        d_mat <- if(curr_rank==1) matrix(svd_fit$d[1],1,1) else diag(svd_fit$d[1:curr_rank])
        X_rec <- svd_fit$u[, 1:curr_rank, drop=FALSE] %*% d_mat %*% t(svd_fit$v[, 1:curr_rank, drop=FALSE])

        curr_vals <- X_filled[na_idx]
        new_vals  <- X_rec[na_idx]
        X_filled[na_idx] <- new_vals

        rss <- sum((curr_vals - new_vals)^2) / (sum(curr_vals^2) + 1e-10)
        if (rss < thresh) break
      }

      X_final_std <- X_filled + mu_time
      X_final <- t(t(X_final_std) * sd_unit + mu_unit)

    } else {
      X_final <- X
    }

    # --- STEP E: REVERSE TRANSFORMATION ---
    if (use_log) {
      X_final <- exp(X_final)
    } else {
      X_final <- sinh(X_final)
    }

    # --- STEP F: Reconstruct ---
    res_wide <- as.data.table(X_final)
    res_wide[, (time_var) := time_index]

    res_long <- data.table::melt(
      res_wide, id.vars = time_var, variable.name = id_var, value.name = "imp_val",
      variable.factor = FALSE
    )

    # ID Restoration
    target_class <- class(DT[[id_var]])[1]
    if ("integer64" %in% class(DT[[id_var]])) {
      res_long[, (id_var) := bit64::as.integer64(as.character(get(id_var)))]
    } else if (is.numeric(DT[[id_var]])) {
      res_long[, (id_var) := as.numeric(as.character(get(id_var)))]
    } else {
      res_long[, (id_var) := as.character(get(id_var))]
    }

    data.table::setkeyv(res_long, c(id_var, time_var))
    data.table::setkeyv(DT, c(id_var, time_var))

    DT[res_long, temp_imp := i.imp_val]

    # --- FIXED: Use get(v) inside the group-by operation ---
    if (!is.infinite(max_hole) || !is.infinite(max_endpoint)) {
      DT[, temp_mask := calc_fill_mask(get(v), max_hole, max_endpoint), by = c(id_var)]
      DT[temp_mask == FALSE, temp_imp := NA]
      DT[, temp_mask := NULL]
    }

    if (is.null(out_suff)) {
      DT[, (v) := data.table::fcoalesce(as.numeric(get(v)), temp_imp)]
    } else {
      out_col <- paste0(v, out_suff)
      DT[, (out_col) := data.table::fcoalesce(as.numeric(get(v)), temp_imp)]
    }

    DT[, c("temp_imp", "temp_v_trans") := NULL]
    rm(wide_dt, X, X_cent, X_filled, res_wide, res_long, time_index)
    gc()
  }

  return(DT)
}

#' Hybrid Imputation (FBI + Kalman)
#'
#' Runs a cascade of imputation: first FBI (Market Factors), then Kalman (Local Trend)
#' on any remaining NAs.
#'
#' @param DT A data.table.
#' @param vars Character vector of variables.
#' @param id_var Character string.
#' @param time_var Character string.
#' @param fbi_rank Integer. Rank for FBI.
#' @param kalman_degree Integer. Degree for Kalman.
#' @param max_hole Integer or Inf.
#' @param max_endpoint Integer or Inf.
#' @param out_suff Character string or NULL.
#' @param verbose Logical.
#'
#' @import data.table
#' @export
impute_hybrid <- function(
    DT, vars, id_var, time_var,
    fbi_rank = 4L,
    kalman_degree = 1L,
    max_hole = 3L,
    max_endpoint = 1L,
    out_suff = "_hybrid", # NULL = Overwrite
    verbose = TRUE
) {
  if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)

  # Create a temporary working copy to avoid messy intermediate overwrites
  working_dt <- if (!is.null(out_suff)) data.table::copy(DT) else DT

  # --- STEP 1: FBI (The Heavy Lifter) ---
  if (verbose) message(">>> STEP 1: Running FBI (Market Structure)...")

  # We run FBI with out_suff = NULL so it fills the gaps in working_dt
  working_dt <- fbi_fast(
    DT = working_dt, vars = vars, id_var = id_var, time_var = time_var,
    rank = fbi_rank,
    max_hole = max_hole,
    max_endpoint = max_endpoint,
    out_suff = NULL,
    verbose = verbose
  )

  # --- STEP 2: Kalman (The Cleanup Crew) ---
  if (verbose) message(">>> STEP 2: Running Kalman on Residual NAs...")

  # Only run Kalman on columns that still have NAs
  vars_with_na <- vars[sapply(working_dt[, ..vars], function(x) any(is.na(x)))]

  if (length(vars_with_na) > 0) {
    working_dt <- kalman_fast(
      DT = working_dt, vars = vars_with_na, id_var = id_var, time_var = time_var,
      degree = kalman_degree,
      max_hole = max_hole,
      max_endpoint = max_endpoint,
      out_suff = NULL,
      verbose = verbose
    )
  }

  # --- STEP 3: Final Column Management ---
  if (!is.null(out_suff)) {
    for (v in vars) {
      new_v <- paste0(v, out_suff)
      # Transfer the filled column to the ORIGINAL DT
      DT[, (new_v) := data.table::fcoalesce(as.numeric(get(v)), working_dt[[v]])]
    }
    return(DT)
  }

  return(working_dt)
}

#' Validation Protocol
#'
#' Masks a percentage of data and tests the accuracy of imputation functions.
#'
#' @param DT Data.table.
#' @param vars Vector of variable names.
#' @param id_var ID variable.
#' @param time_var Time variable.
#' @param impute_fn The function to test (e.g., kalman_fast).
#' @param mask_rate Fraction of data to mask (0-1).
#' @param seed Random seed.
#' @param plot_worst_n Integer. Number of diagnostic plots to return.
#' @param ... Arguments passed to impute_fn.
#'
#' @import data.table
#' @import ggplot2
#' @export
validate_imputation <- function(DT,
                                vars,
                                id_var,
                                time_var,
                                impute_fn,
                                mask_rate = 0.2,
                                seed = 123,
                                plot_worst_n = 4,
                                ...) {

  set.seed(seed)
  message("--- Starting Validation Protocol ---")

  if (!data.table::is.data.table(DT)) DT <- data.table::as.data.table(DT)

  # 1. Create Testing Set
  DT_test <- data.table::copy(DT)
  mask_map <- list()

  for (v in vars) {
    truth_col <- paste0("truth_", v)
    DT_test[, (truth_col) := get(v)]

    # Identify valid indices to drop
    valid_idxs <- DT_test[which(!is.na(get(v))), .I]
    n_drop <- floor(length(valid_idxs) * mask_rate)

    if (n_drop > 0) {
      drop_idx <- sample(valid_idxs, n_drop)
      mask_map[[v]] <- drop_idx
      DT_test[drop_idx, (v) := NA]
    }
  }

  message(sprintf("Masking complete (%.0f%%). Executing Imputation Strategy...", mask_rate * 100))

  # 2. Run Imputation
  start_time <- Sys.time()
  DT_imputed <- impute_fn(DT = DT_test, vars = vars, id_var = id_var, time_var = time_var, ...)

  if (!data.table::is.data.table(DT_imputed)) {
    stop("The 'impute_fn' must return a data.table object.")
  }

  duration <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
  message(sprintf("Imputation finished in %s seconds.", duration))

  # 3. Calculate Stats
  results_list <- list()

  for (v in vars) {
    idx <- mask_map[[v]]
    if (length(idx) == 0) next

    vec_truth <- DT_test[idx, get(paste0("truth_", v))]
    vec_imputed <- DT_imputed[idx, get(v)]

    # Check for complete failure (All NAs)
    if (all(is.na(vec_imputed))) {
      warning(paste0("Variable '", v, "' returned ALL NAs."))
      results_list[[v]] <- data.table(
        Variable = v, RMSE = NA, MAE = NA, NRMSE_Pct = NA,
        Status = "FAILED (All NAs)", Error_Type = "Check Args"
      )
      next
    }

    diffs <- vec_truth - vec_imputed

    # Robust Stats Calculation
    rmse <- sqrt(mean(diffs^2, na.rm = TRUE))
    mae <- mean(abs(diffs), na.rm = TRUE)
    range_val <- max(vec_truth, na.rm = TRUE) - min(vec_truth, na.rm = TRUE)

    nrmse_pct <- if (!is.na(range_val) && range_val > 0) (rmse / range_val) * 100 else NA

    # Status Logic
    status <- "Good"
    if (!is.na(nrmse_pct)) {
      if (nrmse_pct > 20) status <- "CRITICAL"
      else if (nrmse_pct > 10) status <- "Check"
    }

    # Shape Logic (Robust Fix Here)
    shape_flag <- "Stable"

    # FIX: Check is.finite() to catch Inf, and is.na() on the final ratio
    if (is.finite(mae) && mae > 0 && is.finite(rmse)) {
      ratio <- rmse / mae
      if (!is.na(ratio) && ratio > 2) {
        shape_flag <- "Spiky Errors"
      }
    } else if (is.infinite(rmse) || is.infinite(mae)) {
      # Specifically flag divergence
      shape_flag <- "DIVERGED (Inf)"
      status <- "CRITICAL"
    }

    results_list[[v]] <- data.table(
      Variable = v,
      RMSE = round(rmse, 4),
      MAE = round(mae, 4),
      NRMSE_Pct = round(nrmse_pct, 2),
      Status = status,
      Error_Type = shape_flag
    )
  }

  stats_dt <- rbindlist(results_list)
  if (nrow(stats_dt) > 0) data.table::setorder(stats_dt, -NRMSE_Pct)

  # 4. Generate Diagnostic Plots
  plots_out <- list()
  if (nrow(stats_dt) > 0) {
    vars_to_plot <- head(stats_dt$Variable, plot_worst_n)

    for (v in vars_to_plot) {
      masked_idxs <- mask_map[[v]]
      affected_ids <- unique(DT_test[masked_idxs, get(id_var)])

      if (length(affected_ids) > 0) {
        demo_id <- affected_ids[1]

        # Slices
        truth_slice <- DT_test[get(id_var) == demo_id]
        imp_slice   <- DT_imputed[get(id_var) == demo_id]
        mask_points <- truth_slice[, .I] %in% masked_idxs

        # Plot
        p <- ggplot() +
          geom_line(data = truth_slice, aes(x = get(time_var), y = get(paste0("truth_", v))),
                    color = "grey70", size = 0.8, alpha = 0.6) +
          geom_line(data = imp_slice, aes(x = get(time_var), y = get(v)),
                    color = "steelblue", size = 0.8) +
          geom_point(data = truth_slice[mask_points], aes(x = get(time_var), y = get(paste0("truth_", v))),
                     color = "red", size = 3, shape = 4, stroke = 2) +
          labs(
            title = paste0("Diagnostic: ", v, " (ID: ", demo_id, ")"),
            subtitle = paste0("NRMSE: ", stats_dt[Variable == v, NRMSE_Pct], "% | Status: ", stats_dt[Variable == v, Status]),
            y = "Value", x = "Time"
          ) +
          theme_minimal()

        plots_out[[v]] <- p
      }
    }
  }

  return(list(stats = stats_dt, plots = plots_out))
}
