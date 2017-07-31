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
#' 
#' It returns a bool of whether scDD's predicted FDR can be trusted or not.
#' This is decided by a classifier based on a flexible discriminant analysis.
#' 
#' If scDD's FDR control is within an acceptable range, we encourage using
#' the above mentioned method with edgeR and summed up batches (Lun, 2017) 
#' and scDD together.
#' scDD is quite good at detecting different shapes in distributions while
#' edgeR is good at detecting different means.
#'
#' @param dt        \code{num} matrix or data.frame of (not log transformed)
#'                  expression values with rows as genes and columns as cells
#' @param batches   vector identifying batches (\code{chr}, or \code{int})
#' @param groups    vector identifying batches (\code{chr}, or \code{int})
#' @param FDR       \code{num} scalar of predicted FDR
#'
#' @return list of \code{bool} whether scDD's FDR can be trusted
#'         and \code{num} of posterior probability from the flexible 
#'         discriminant analysis
#'
#' @export
#' 
#' @examples 
#' ds <- SimulateData()
#' res <- EstimFDRcontrol(ds$table, ds$pData$batch, ds$pData$group, 0.1)
#' res
#' ds <- SimulateData(nCells = 50)
#' res <- EstimFDRcontrol(ds$table, ds$pData$batch, ds$pData$group, 0.1)
#' res
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
  n_hat <- as.vector(ncol(dt) / length(ubatches))
  
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
  
  # mean count
  mean_count <- mean(as.vector(as.matrix(dt)))
  
  # test scDD
  X <- data.frame(
    n_hat = n_hat,
    cv_hat_within = cv_hat_within,
    cv_hat_between = cv_hat_between,
    mean_count = mean_count,
    predFDR = FDR
  )
  post <- as.vector(predict(model, X, type = "posterior")[, 2])
  lab <- post > 0.0008
  list(use_scDD = lab, posterior = post)
}
