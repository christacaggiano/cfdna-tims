suppressPackageStartupMessages({
  library(data.table)
  library(bigstatsr)
  library(caret)
  library(pROC)
  library(dplyr)
})

# ------------------------- Config -------------------------

ucla_meth_path <- "ucla.txt"
uq_meth_path   <- "uq.txt"
ucla_cov_path  <- "ucla_covariates.txt"
uq_cov_path    <- "uq_covariates.txt"
meta_path      <- "metadata_final.txt"

manual_remove <- c("45", "96", "109")

covars_ucla <- c("sex female", "age", "starting", "concentration", "race white")
covars_uq   <- c("sex female", "age", "starting", "concentration")

# ------------------------ Helpers -------------------------

read_meth_as_samples <- function(path) {
  df <- fread(path, data.table = FALSE)
  if (!is.numeric(df[[1]])) df <- df[, -1, drop = FALSE]  # drop CpG index column if present
  as.data.frame(t(df))                                    # rows = samples, cols = CpGs
}

to_binary01 <- function(x) {
  x <- as.vector(unlist(x))
  ifelse(x %in% c("case","ALS","TRUE","1",1,TRUE), 1, 0)
}

minmax_scale <- function(Xtr, Xte) {
  mins <- apply(Xtr, 2, min)
  maxs <- apply(Xtr, 2, max)
  rng  <- pmax(maxs - mins, .Machine$double.eps)
  Xtr2 <- sweep(sweep(Xtr, 2, mins, "-"), 2, rng, "/")
  Xte2 <- sweep(sweep(Xte, 2, mins, "-"), 2, rng, "/")
  list(train = Xtr2, test = Xte2)
}

safe_FBM <- function(M, path="~/temp_matrix") {
  system("rm -f ~/temp_matrix.bk")
  as_FBM(as.matrix(M), backingfile = path)
}

report_auc <- function(y, score, label) {
  a <- as.numeric(pROC::auc(y, score))
  cat(sprintf("%s AUC: %.4f\n", label, a))
  a
}

reduce_covariates <- function(C_train) {
  if (is.null(C_train) || ncol(C_train) == 0) return(list(C = NULL, idx = integer(0)))
  keep <- which(apply(C_train, 2, function(v) length(unique(v)) > 1))
  if (length(keep) == 0) return(list(C = NULL, idx = integer(0)))
  list(C = C_train[, keep, drop = FALSE], idx = keep)
}

build_cohort <- function(cohort, meth_path, cov_df, meta_df, covars_keep) {
  md <- meta_df %>%
    filter(group == cohort) %>%
    filter(!sample_num %in% manual_remove) %>%
    filter(sample_type != "PLS")

  X <- read_meth_as_samples(meth_path)

  if (all(grepl("^V[0-9]+$", rownames(X)))) {
    rn <- as.character(cov_df$sample_num)
    n  <- min(nrow(X), length(rn))
    rownames(X)[seq_len(n)] <- rn[seq_len(n)]
  }

  keep_ids <- intersect(as.character(md$sample_num), as.character(rownames(X)))
  X <- X[keep_ids, , drop = FALSE]

  cov_df <- cov_df[match(keep_ids, as.character(cov_df$sample_num)), , drop = FALSE]
  C <- cov_df[, intersect(c("sample_num", covars_keep), names(cov_df)), drop = FALSE]
  C <- C[, setdiff(names(C), "sample_num"), drop = FALSE]

  y <- to_binary01(md$sample_type[match(keep_ids, md$sample_num)])

  stopifnot(nrow(X) == length(y), nrow(C) == length(y))

  list(X = X, C = C, y = y, ids = keep_ids)
}

xcohort_auc <- function(train, test, label, alphas = c(1e-4, 0.01, 0.1, 0.5)) {
  scaled <- minmax_scale(train$X, test$X)

  red <- reduce_covariates(train$C)
  Ctr <- red$C
  Cte <- if (length(red$idx)) test$C[, red$idx, drop = FALSE] else NULL

  fit <- big_spLogReg(
    X = safe_FBM(scaled$train),
    y01.train = train$y,
    alphas = alphas,
    covar.train = if (!is.null(Ctr)) data.matrix(Ctr) else NULL,
    warn = FALSE
  )

  preds <- predict(fit, safe_FBM(scaled$test),
                   covar.row = if (!is.null(Cte)) data.matrix(Cte) else NULL)

  report_auc(test$y, preds, label)
}

within_cv_auc <- function(dat, label, k = 10, seed = 1, alphas = c(1e-4, 0.01, 0.1, 0.5)) {
  set.seed(seed)
  folds <- createFolds(dat$y, k = k, list = TRUE)
  pred <- rep(NA_real_, length(dat$y))

  for (i in seq_along(folds)) {
    idx_te <- folds[[i]]
    idx_tr <- setdiff(seq_along(dat$y), idx_te)

    scaled <- minmax_scale(dat$X[idx_tr, , drop = FALSE],
                           dat$X[idx_te, , drop = FALSE])

    red <- reduce_covariates(dat$C[idx_tr, , drop = FALSE])
    Ctr <- red$C
    Cte <- if (length(red$idx)) dat$C[idx_te, red$idx, drop = FALSE] else NULL

    fit <- big_spLogReg(
      X = safe_FBM(scaled$train),
      y01.train = dat$y[idx_tr],
      alphas = alphas,
      covar.train = if (!is.null(Ctr)) data.matrix(Ctr) else NULL,
      warn = FALSE
    )

    pred[idx_te] <- predict(
      fit,
      safe_FBM(scaled$test),
      covar.row = if (!is.null(Cte)) data.matrix(Cte) else NULL
    )
  }

  report_auc(dat$y, pred, label)
}

combined_cv_auc <- function(ucla, uq, label = "Combined (CV)", k = 10, seed = 1,
                            alphas = c(1e-4, 0.01, 0.1, 0.5)) {
  X <- rbind(ucla$X, uq$X)
  C <- rbind(ucla$C, uq$C)
  y <- c(ucla$y, uq$y)

  set.seed(seed)
  folds <- createFolds(y, k = k, list = TRUE)
  pred <- rep(NA_real_, length(y))

  for (i in seq_along(folds)) {
    idx_te <- folds[[i]]
    idx_tr <- setdiff(seq_along(y), idx_te)

    scaled <- minmax_scale(X[idx_tr, , drop = FALSE],
                           X[idx_te, , drop = FALSE])

    red <- reduce_covariates(C[idx_tr, , drop = FALSE])
    Ctr <- red$C
    Cte <- if (length(red$idx)) C[idx_te, red$idx, drop = FALSE] else NULL

    fit <- big_spLogReg(
      X = safe_FBM(scaled$train),
      y01.train = y[idx_tr],
      alphas = alphas,
      covar.train = if (!is.null(Ctr)) data.matrix(Ctr) else NULL,
      warn = FALSE
    )

    pred[idx_te] <- predict(
      fit,
      safe_FBM(scaled$test),
      covar.row = if (!is.null(Cte)) data.matrix(Cte) else NULL
    )
  }

  report_auc(y, pred, label)
}

# ------------------------ Load data ------------------------

ucla_meth <- read_meth_as_samples(ucla_meth_path)
uq_meth   <- read_meth_as_samples(uq_meth_path)
ucla_cov  <- fread(ucla_cov_path, data.table = FALSE)
uq_cov    <- fread(uq_cov_path,   data.table = FALSE)
meta      <- fread(meta_path,     data.table = FALSE)

ucla_cov$sample_num <- as.character(ucla_cov$sample_num)
uq_cov$sample_num   <- as.character(uq_cov$sample_num)
meta$sample_num     <- as.character(meta$sample_num)

# --------------------- Build cohorts ----------------------

ucla <- build_cohort("UCLA", ucla_meth_path, ucla_cov, meta, covars_ucla)
uq   <- build_cohort("UQ",   uq_meth_path,   uq_cov,   meta, covars_uq)

cat(sprintf("UCLA: %d samples, %d CpGs, %d covariates\n",
            nrow(ucla$X), ncol(ucla$X), ncol(ucla$C)))
cat(sprintf("UQ:   %d samples, %d CpGs, %d covariates\n",
            nrow(uq$X),   ncol(uq$X),   ncol(uq$C)))

stopifnot(all(ucla$y %in% c(0,1)), all(uq$y %in% c(0,1)))

# ---------------------- Experiments -----------------------

auc_ucla_to_uq <- xcohort_auc(ucla, uq,   "UCLA→UQ")
auc_uq_to_ucla <- xcohort_auc(uq,   ucla, "UQ→UCLA")

auc_ucla_cv <- within_cv_auc(ucla, "UCLA (CV)")
auc_uq_cv   <- within_cv_auc(uq,   "UQ (CV)")

auc_combined <- combined_cv_auc(ucla, uq, "Combined (CV)")

cat("\n===== SUMMARY =====\n")
cat(sprintf("UCLA→UQ:   %.4f\n", auc_ucla_to_uq))
cat(sprintf("UQ→UCLA:   %.4f\n", auc_uq_to_ucla))
cat(sprintf("UCLA (CV): %.4f\n", auc_ucla_cv))
cat(sprintf("UQ (CV):   %.4f\n", auc_uq_cv))
cat(sprintf("Combined:   %.4f\n", auc_combined))
