#' EstimFDRcontrol
#'
#' Estimate loss of FDR control for scDD using pooled cells and the 
#' Kolmogorov-Smirnov test.
#' Quality controlled, library sized normalized (and not log transformed)
#' values should be used.
#' 
#' When measurement comming from different batches are treated 
#' equally independent (pooling cells), then methods for detecting 
#' differentially distributed genes lose control over the false discovery rate
#' (FDR).
#' scDD is moderately robust to this effect, yet quite powerful.
#' This function uses your dataset to predict how much scDD underestimates
#' the FDR.
#' The predicted true FDR is returned (corrected for scDD's FDR control loss).
#' 
#' scDD is one of the most robust method with regard to batch effects
#' (when pooling cells).
#' If it's loss of FDR control is to high then you need to resort to drastic
#' methods such as proposed by A. Lun, Biostatistics, 2017.
#' Here, cells are summed up within batches and then edgeR is used with these
#' sums.
#' That means you effectively have to degrade your single cell experiment
#' to a bulk experiment.
#' 
#' If scDD's FDR control is within an acceptable range, we encourage using
#' the above mentioned method with edgeR (Lun, 2017) and scDD together.
#' scDD is quite good at detecting different shapes in distributions while
#' edgeR is good at detecting different means.
#' 
#' @note A linear mixed model is fitted iteratively to each gene, so this
#' function can take a while.
#'
#' @param dt        \code{num} matrix or data.frame of (not log transformed)
#'                  expression values with rows as genes and columns as cells
#' @param batches   vector identifying batches (\code{chr}, or \code{int})
#' @param groups    vector identifying batches (\code{chr}, or \code{int})
#' @param FDR       \code{num} scalar of predicted FDR
#'
#' @return \code{num} for true (corrected) FDR
#'
#' @export
#' 
#' @examples 
#' ds <- SimulateData()
#' trueFDR <- EstimFDRcontrol(ds$table, ds$pData$batch, ds$pData$group, 0.1)
#' trueFDR
#'
EstimFDRcontrol <- function(dt, batches, groups, FDR) {
  stopifnot(inherits(dt, "data.frame") || inherits(dt, "matrix"))
  stopifnot(length(batches) == length(groups) && length(groups) == ncol(dt))
  batches <- as.factor(batches)
  groups <- as.factor(groups)
  ubatches <- unique(batches)
  ugroups <- unique(groups)
  if (length(ugroups) != 2) stop("Only compare 2 groups.")
  idxs <- list(which(groups == ugroups[1]), which(groups == ugroups[2]))
  if (length(intersect(idxs[[1]], idxs[[2]])) != 0L) {
    stop("There are batches which appear in both groups.
         Congratulations, your experiment is not completely confounded.
         You don't need this function.")
  }
  if (length(unique(batches[idxs[[1]]])) < 2 || 
      length(unique(batches[idxs[[2]]])) < 2) {
    stop("For at least one group there is less than 1 batch.
         That means I cannot estimate the batch effect.")
  }
  
  # cells/batch
  n_hat <- ncol(dt) / length(ubatches)
  
  # within batch var
  cv <- function(x) sd(x) / mean(x)
  meanCVs <- numeric(length(ubatches))
  for (i in seq_along(ubatches)) {
    x <- apply(dt[, batches == ubatches[i]], 1, cv)
    x[!is.finite(x)] <- NA
    meanCVs[i] <- mean(x, na.rm = TRUE)
  }
  cv_hat_within <- mean(meanCVs)
  
  # between batch var
  CVs <- numeric(length(ugroups))
  for (i in seq_along(idxs)) {
    l <- lapply(unique(batches[idxs[[i]]]), function(batch) {
      rowMeans(dt[, batches == batch])
    })
    CVs[i] <- mean(apply(do.call(cbind, l), 1, cv), na.rm = TRUE)
  }
  cv_hat_between <- mean(CVs)
  
  # LMM
  corrs <- limma::duplicateCorrelation(
    log(dt + 1), 
    design = model.matrix(~ 1 + groups), 
    block = batches
  )
  batch_corrs <- corrs$consensus.correlation
  
  # test scDD
  df <- data.frame(
    n_hat = n_hat,
    cv_hat_within = cv_hat_within,
    cv_hat_between = cv_hat_between,
    batch_corrs = batch_corrs,
    predFDR = FDR
  )
  dm <- xgboost::xgb.DMatrix(data = as.matrix(df))
  scDDpred <- predict(model.scDD, dm)
  scDDpred + FDR
}
