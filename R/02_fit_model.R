fit_fire_model <- function(
  prepared_data_file = "data/processed/fire_stan_data.rds",
  stan_file = "stan/firemodel.stan",
  output_file = "results/model/fire_stan_fit.rds",
  chains = 3,
  iter_warmup = 2000,
  iter_sampling = 10000,
  thin = 10,
  seed = 1234,
  parallel_chains = chains
) {
  prepared <- readRDS(prepared_data_file)
  stan_data <- prepared$stan_data

  if (!file.exists(stan_file)) {
    stop("Stan model not found: ", stan_file)
  }

  engine <- NULL
  fit <- NULL

  cmdstan_ready <- FALSE
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    cmdstan_ready <- tryCatch({
      !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
    }, error = function(e) {
      FALSE
    })
  }

  if (cmdstan_ready) {
    engine <- "cmdstanr"
    message("Using cmdstanr backend")
    mod <- cmdstanr::cmdstan_model(stan_file)
    fit <- mod$sample(
      data = stan_data,
      seed = seed,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      thin = thin,
      refresh = 200
    )

    beta_summary <- fit$summary(variables = "beta")
    diagnostics <- beta_summary[, c("variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail")]
  } else if (requireNamespace("rstan", quietly = TRUE)) {
    engine <- "rstan"
    message("Using rstan backend")
    fit <- rstan::sampling(
      object = rstan::stan_model(file = stan_file),
      data = stan_data,
      chains = chains,
      warmup = iter_warmup,
      iter = (iter_warmup + iter_sampling)/thin,
      seed = seed,
      refresh = 200
    )

    s <- rstan::summary(fit, pars = "beta", probs = c(0.025, 0.5, 0.975))$summary
    diagnostics <- data.frame(
      variable = rownames(s),
      mean = s[, "mean"],
      sd = s[, "sd"],
      rhat = s[, "Rhat"],
      ess_bulk = s[, "n_eff"],
      ess_tail = NA_real_,
      row.names = NULL,
      check.names = FALSE
    )
  } else {
    stop("Neither cmdstanr (with CmdStan installed) nor rstan is available.")
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  out <- list(
    engine = engine,
    fit = fit,
    diagnostics = diagnostics,
    metadata = prepared$metadata,
    prepared_data_file = prepared_data_file,
    stan_file = stan_file
  )
  saveRDS(out, output_file)

  message("Saved fitted model to: ", output_file)
  print(diagnostics)
  invisible(out)
}

if (sys.nframe() == 0) {
  fit_fire_model()
}
