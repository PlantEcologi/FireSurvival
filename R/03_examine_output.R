extract_beta_draws <- function(fit_obj) {
  engine <- fit_obj$engine
  fit <- fit_obj$fit

  if (engine == "cmdstanr") {
    draws_array <- fit$draws(variables = "beta", format = "draws_array")
    draws_matrix <- fit$draws(variables = "beta", format = "matrix")
    return(list(array = draws_array, matrix = draws_matrix))
  }

  if (engine == "rstan") {
    arr <- rstan::extract(fit, pars = "beta", permuted = FALSE)
    # rstan returns iter x chains x params for a single parameter vector
    dimnames(arr)[[3]] <- paste0("beta[", seq_len(dim(arr)[3]), "]")
    mat <- as.matrix(fit, pars = "beta")
    return(list(array = arr, matrix = mat))
  }

  stop("Unknown engine: ", engine)
}

summarise_beta <- function(draws_matrix, coef_names) {
  q <- apply(draws_matrix, 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
  data.frame(
    coefficient = coef_names,
    mean = colMeans(draws_matrix),
    sd = apply(draws_matrix, 2, stats::sd),
    q2.5 = q[1, ],
    q50 = q[2, ],
    q97.5 = q[3, ],
    row.names = NULL,
    check.names = FALSE
  )
}

plot_trace <- function(draws_array, coef_names, out_file) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    message("Package 'coda' not installed; skipping trace plot.")
    return(invisible(NULL))
  }

  n_iter <- dim(draws_array)[1]
  n_chain <- dim(draws_array)[2]
  n_coef <- dim(draws_array)[3]

  ml <- coda::mcmc.list(lapply(seq_len(n_chain), function(ch) {
    m <- matrix(draws_array[, ch, ], nrow = n_iter, ncol = n_coef)
    colnames(m) <- coef_names
    coda::mcmc(m)
  }))

  grDevices::pdf(out_file, width = 12, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  coda::traceplot(ml)
}

plot_posterior_hist <- function(draws_matrix, coef_names, out_file) {
  grDevices::pdf(out_file, width = 12, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  n_coef <- length(coef_names)
  n_col <- ceiling(sqrt(n_coef))
  n_row <- ceiling(n_coef / n_col)
  par(mfrow = c(n_row, n_col), mar = c(3, 3, 2, 1))
  for (j in seq_along(coef_names)) {
    hist(draws_matrix[, j], main = coef_names[j], xlab = "", col = "grey80", border = "white")
    abline(v = mean(draws_matrix[, j]), col = "firebrick", lwd = 2)
  }
}

plot_forest <- function(summary_df, out_file) {
  ord <- seq_len(nrow(summary_df))
  grDevices::pdf(out_file, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(
    summary_df$q50[ord], ord,
    xlim = range(c(summary_df$q2.5, summary_df$q97.5)),
    yaxt = "n", ylab = "", xlab = "Coefficient (median and 95% CrI)", pch = 16
  )
  axis(2, at = ord, labels = summary_df$coefficient[ord], las = 2)
  segments(summary_df$q2.5[ord], ord, summary_df$q97.5[ord], ord, lwd = 2)
  abline(v = 0, lty = 2)
}

examine_fire_model <- function(
  fit_file = "results/model/fire_stan_fit.rds",
  out_dir = "results/summary"
) {
  fit_obj <- readRDS(fit_file)
  draws <- extract_beta_draws(fit_obj)

  coef_names <- fit_obj$metadata$x_colnames
  if (is.null(coef_names)) {
    coef_names <- colnames(draws$matrix)
  }

  # Ensure matrix columns align with model coefficient order.
  if (ncol(draws$matrix) == length(coef_names)) {
    colnames(draws$matrix) <- coef_names
  }

  summary_df <- summarise_beta(draws$matrix, coef_names)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(summary_df, file.path(out_dir, "beta_summary.csv"), row.names = FALSE)

  original_order <- c(
    "intercept", "TAAve", "TSAve", "TSHotweek",
    "PATot", "PSTot", "aao", "pptconc", "spring", "summer", "fall"
  )
  if (length(coef_names) != length(original_order)) {
    warning(
      "Coefficient length mismatch: expected ", length(original_order),
      ", got ", length(coef_names), ". Marking matches as NA."
    )
    stan_order_aligned <- rep(NA_character_, length(original_order))
    n_fill <- min(length(coef_names), length(original_order))
    stan_order_aligned[seq_len(n_fill)] <- coef_names[seq_len(n_fill)]
    compare_df <- data.frame(
      original_order = original_order,
      stan_order = stan_order_aligned,
      matched = NA,
      row.names = NULL,
      check.names = FALSE
    )
  } else {
    compare_df <- data.frame(
      original_order = original_order,
      stan_order = coef_names,
      matched = original_order == coef_names,
      row.names = NULL,
      check.names = FALSE
    )
  }
  utils::write.csv(compare_df, file.path(out_dir, "coefficient_order_comparison.csv"), row.names = FALSE)

  # Acceptance-rate summaries
  accept_df <- NULL
  if (fit_obj$engine == "cmdstanr") {
    sampler_diag_mat <- fit_obj$fit$sampler_diagnostics(format = "matrix")
    accept_df <- data.frame(mean_accept_stat = mean(sampler_diag_mat[, "accept_stat__"]))
  } else if (fit_obj$engine == "rstan") {
    sp <- rstan::get_sampler_params(fit_obj$fit, inc_warmup = FALSE)
    accepts <- unlist(lapply(sp, function(x) x[, "accept_stat__"]))
    accept_df <- data.frame(mean_accept_stat = mean(accepts))
  }

  if (!is.null(accept_df)) {
    utils::write.csv(accept_df, file.path(out_dir, "acceptance_summary.csv"), row.names = FALSE)
  }

  if (!is.null(fit_obj$diagnostics)) {
    utils::write.csv(fit_obj$diagnostics, file.path(out_dir, "diagnostics_rhat_ess.csv"), row.names = FALSE)
  }

  plot_trace(draws$array, coef_names, file.path(out_dir, "traceplots_beta.pdf"))
  plot_forest(summary_df, file.path(out_dir, "forestplot_beta.pdf"))
  plot_posterior_hist(draws$matrix, coef_names, file.path(out_dir, "posterior_hist_beta.pdf"))

  message("Saved summaries and plots to: ", out_dir)
  invisible(list(summary = summary_df, comparison = compare_df, acceptance = accept_df))
}

if (sys.nframe() == 0) {
  examine_fire_model()
}
