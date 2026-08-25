release_urls <- c(
  conc = "https://github.com/AdamWilsonLab/FireSurvival2010/releases/download/Data/conc.csv",
  data = "https://github.com/AdamWilsonLab/FireSurvival2010/releases/download/Data/data.csv",
  static = "https://github.com/AdamWilsonLab/FireSurvival2010/releases/download/Data/static.csv"
)

zscore <- function(x) {
  as.numeric(scale(x))
}

download_release_data <- function(data_dir = "data/raw", overwrite = FALSE) {
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(release_urls)) {
    out_file <- file.path(data_dir, paste0(nm, ".csv"))
    if (!file.exists(out_file) || isTRUE(overwrite)) {
      message("Downloading ", nm, " -> ", out_file)
      utils::download.file(release_urls[[nm]], out_file, mode = "wb", quiet = FALSE)
    }
  }

  invisible(
    list(
      conc = file.path(data_dir, "conc.csv"),
      data = file.path(data_dir, "data.csv"),
      static = file.path(data_dir, "static.csv")
    )
  )
}

prepare_fire_data <- function(
  data_dir = "data/raw",
  output_file = "data/processed/fire_stan_data.rds",
  left_censor_length = 112L,
  scale_pptconc = TRUE
) {
  paths <- download_release_data(data_dir = data_dir)

  data <- utils::read.csv(paths$data)
  static <- utils::read.csv(paths$static)
  conc <- utils::read.csv(paths$conc)

  # Legacy script dropped first column in conc.csv; do that only when appropriate.
  if (ncol(conc) > 1 && !"pid" %in% names(conc) && "pid" %in% names(conc[, -1, drop = FALSE])) {
    conc <- conc[, -1, drop = FALSE]
  }

  static <- static[order(static$gid), ]
  data <- data[order(data$gid, data$sid), ]

  static2 <- merge(static, conc, by = "pid", all.x = TRUE, sort = FALSE)
  static2 <- static2[order(static2$gid), ]

  data <- merge(
    data,
    unique(static2[c("gid", "pptconc")]),
    by = "gid",
    all.x = TRUE,
    sort = FALSE
  )
  data <- data[order(data$gid, data$sid), ]

  if (isTRUE(scale_pptconc)) {
    data$pptconc <- zscore(data$pptconc)
  }

  data$aao[is.na(data$aao)] <- mean(data$aao, na.rm = TRUE)

  # Legacy season coding from source data: 1 = fall, 3 = spring, 4 = summer.
  data$spring <- as.integer(data$season == 3)
  data$summer <- as.integer(data$season == 4)
  data$fall <- as.integer(data$season == 1)

  x_cols <- c(
    "intercept", "TAAve", "TSAve", "TSHotweek",
    "PATot", "PSTot", "aao", "pptconc", "spring", "summer", "fall"
  )

  data$intercept <- 1
  X <- as.matrix(data[, x_cols])

  gid_unique <- sort(unique(data$gid))
  n_grid <- length(gid_unique)
  N <- nrow(data)
  seasons_per_grid <- as.integer(table(factor(data$gid, levels = gid_unique)))
  n_seas <- as.integer(max(seasons_per_grid))
  if (length(unique(seasons_per_grid)) != 1L) {
    warning("Grid cells have unequal numbers of seasons; using max seasons per grid in metadata.")
  }

  gid_index <- split(seq_len(N), data$gid)

  fire <- as.integer(data$fire)

  # Left-censor handling (legacy non-spatial model): first 112 seasons were pre-1980.
  left_censor_index <- rep(FALSE, N)
  for (g in names(gid_index)) {
    idx <- gid_index[[g]]
    take_n <- min(length(idx), left_censor_length)
    left_censor_index[idx[seq_len(take_n)]] <- TRUE
  }

  # Exclude the season immediately after any observed fire (legacy likelihood behavior).
  post_fire_exclude <- rep(FALSE, N)
  for (g in names(gid_index)) {
    idx <- gid_index[[g]]
    if (length(idx) <= 1L) next
    fire_local <- which(fire[idx] == 1L)
    fire_local <- fire_local[fire_local < length(idx)]
    if (length(fire_local) > 0L) {
      post_fire_exclude[idx[fire_local + 1L]] <- TRUE
    }
  }

  include_row <- !(left_censor_index | post_fire_exclude)

  # Keep a first-observed-fire index for reference/diagnostics.
  first_observed_fire <- vapply(gid_unique, function(g) {
    idx <- which(data$gid == g)
    fire_idx <- idx[fire[idx] == 1L]
    if (length(fire_idx) == 0L) {
      return(NA_integer_)
    }
    min(fire_idx)
  }, integer(1))

  stan_data <- list(
    N = N,
    K = ncol(X),
    X = X,
    y = fire,
    include_row = as.integer(include_row)
  )

  prepared <- list(
    stan_data = stan_data,
    metadata = list(
      x_colnames = colnames(X),
      gid_unique = gid_unique,
      n_grid = n_grid,
      n_seas = n_seas,
      seasons_per_grid = seasons_per_grid,
      left_censor_length = left_censor_length,
      first_observed_fire_row = first_observed_fire,
      include_row_count = sum(include_row),
      exclude_row_count = sum(!include_row)
    ),
    data = data
  )

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(prepared, output_file)
  message("Saved prepared data to: ", output_file)

  invisible(prepared)
}

if (sys.nframe() == 0) {
  prepare_fire_data()
}
