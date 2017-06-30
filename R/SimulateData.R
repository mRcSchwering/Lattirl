#' PoiBeta
#'
#' Create Poisson-Beta distributions.
#'
#' The Poisson-Beta distribution has mean and variance of
#' \itemize{
#'   \item Mu = k * a / (a + b)
#'   \item Sig^2 = k^2 * a * b / ((a + b + 1) * (a + b)^2)
#' }
#' If k, a, or b have length greater 1, a matrix with n columns is created.
#'
#' @param k \code{numeric} vector for k.
#'        if length > 1 then several rows are computed
#' @param a \code{numeric} vector for a
#' @param b \code{numeric} vector for b
#' @param n \code{integer} number of values (number of columns in value)
#' @param d \code{numeric} vector for probability of zero counts
#'
#' @return \code{numeric} matrix of length k, a, or b rows and n columns
#'
#'
PoiBeta <- function(k, a, b, n = 100, d = 0) {
  if (a<=0 || b<=0 || k<=0 || n<=0) {
    stop("All parameters must be positive!")
  }
  
  #Check if we have vectors and that they are the same length
  kLen <- length(k)
  aLen <- length(a)
  bLen <- length(b)
  if ((kLen==1 && aLen==1 && bLen==1) ||
      (kLen>1 && aLen==1 && bLen==1) ||
      (kLen==1 && aLen>1 && bLen==1) ||
      (kLen==1 && aLen==1 && bLen>1) ||
      (kLen==bLen && aLen==bLen) ||
      (kLen==bLen && aLen==1) ||
      (kLen==aLen && bLen==1) ||
      (kLen==1 && aLen==bLen)) {
    m <- max(c(kLen, aLen, bLen))
    x <- mat.or.vec(m, n)
    if (kLen==1) k <- rep(k, m)
    if (aLen==1) a <- rep(a, m)
    if (bLen==1) b <- rep(b, m)
    if (length(d) == 1) d <- rep(d, m)
    
    # go through features
    for (i in 1:m) {
      
      #First generate Beta random variables
      y <- rbeta(n, a[i], b[i])
      
      #Then use for the Poisson intensities
      x[i,] <- rpois(n, k[i]*y)
      
      # add dropouts
      Nz <- floor(d[i] * n) - sum(x[i,] == 0)
      if ( Nz > 0 ) {
        x[i, sample(which(x[i,] > 0), Nz)] <- 0
      }
    }
    return(x)
    
  } else {
    stop("Array lengths of input parameters inconsistent!")
  }
}




#' SimulateData
#'
#' Generate a simulated scRNA-seq dataset of 2 biological groups with several
#' batches per group, a defined variability within and between batches, and
#' a defined proportion and amount of differentially distributed genes.
#' 
#' Count tables, phenotypic- and feature-information tables are generated
#' to resemble a scRNA-seq experiment.
#' Count distributions are modelled using a 3 parameter Beta-Poisson
#' distributions: \emph{Pois(k * Beta(a, b))}.
#'
#' For low average expression (genewise) additional 0's are introduced to better
#' resemble real dropout rates.
#' This is done according to a Hill equation, which defines the propability of
#' a zero (dropout) from the average expression by
#' \emph{P_zero = x^n / (Km + x^n)} where x is the average expression.
#'
#' The defaults resemble the \emph{Tung} 2016 dataset.
#' In general, these parameters fit well to datasets which were produced using
#' UMI's and the fluidigm platform.
#' For raw readcounts the absolute count values and variabilities are higher.
#' Count tables produced using Drop-Seq usually have more cells and a much
#' larger dropout rate.
#'
#' @param k,a,b       \code{numeric} scalar for average BetaPoisson parameters
#'                    in log2
#' @param _CV         \code{numeric} scalar for coefficient of variation of
#'                    a parameter with a batch (in log2)
#' @param _CV2        \code{numeric} scalar for coefficient of variation of
#'                    average parameters between batches (in log2)
#' @param Km,n        \code{numeric} parameter scalar defining Hill function
#'                    by which 0's (dropouts) will be introduced
#' @param pZero_SD    \code{numeric} scalar for Std of introduced dropouts
#' @param nBatches    \code{integer} scalar number of batches per group
#' @param nGenes      \code{integer} scalar number of genes
#' @param nCells      \code{integer} scalar average number of cells per batch
#' @param nCells_SD   \code{numeric} scalar coefficient of deviation for
#'                    numbers of cells in batches
#' @param ddP         \code{numeric} scalar for proportion of differentially
#'                    distributed genes
#' @param ddM         \code{numeric} scalar controlling amount of differential
#'                    distribution (3 would be quite strong)
#'
#' @return list of generated count table, phenotypic- (cell-) and feature-
#'         (gene-) information tables, and a list which holds the Parameters
#'         for each batch and each gene.
#'
#' @export
#' 
#' @examples 
#' ds <- ds <- SimulateData()
#' str(ds)
#'
SimulateData <- function(
  k = 4.2, k_CV = 0.25, k_CV2 = 0.08,
  a = 1.9, a_CV = 0.35, a_CV2 = 0.07,
  b = 3.5, b_CV = 0.40, b_CV2 = 0.09,
  Km = 4.5, n = 3.4,
  pZero_SD = 0.10,
  nBatches = 3,
  nGenes = 10000,
  nCells = 100,
  nCells_CV = 0.3,
  ddP = 0, ddM = 1
) {
  # create fData
  fData <- data.table::data.table(
    gene = paste0("Gene", 1:nGenes),
    label = "random"
  )
  # define differentially distributed genes
  ddN <- round(nGenes * ddP)
  if (ddN > 10) {
    steps <- round(seq(0, ddN, length.out = 5))
    fData[1:steps[2], label := "mean"]
    fData[(steps[2] + 1):steps[3], label := "shape"]
    fData[(steps[3] + 1):steps[4], label := "both_a"]
    fData[(steps[4] + 1):steps[5], label := "both_b"]
  }
  
  # generate PoiBeta parameter distributions
  ks <- rnorm(nGenes, k, k_CV * k)
  as <- rnorm(nGenes, a, a_CV * a)
  bs <- rnorm(nGenes, b, b_CV * b)
  
  # create factors for differential distros
  if (ddM <= 0) ddM <- 1
  kM <- rep(1, nGenes)
  aM <- rep(1, nGenes)
  bM <- rep(1, nGenes)
  kM[fData[, label == "mean"]] <- ddM
  aM[fData[, label == "shape"]] <- ddM * 1 / 10^(ddM - 1)
  bM[fData[, label == "shape"]] <- ddM * 1 / 10^(ddM - 1)
  aM[fData[, label == "both_a"]] <- ddM * 1.5
  bM[fData[, label == "both_b"]] <- ddM * 2
  
  # generate parameters for each batch
  nSmpls <- nBatches * 2
  pars <- vector("list", nSmpls)
  for (i in 1:nSmpls) {
    
    # add batch effects
    ks_btch <- ks + rnorm(nGenes, 0, abs(k_CV2 * ks))
    as_btch <- as + rnorm(nGenes, 0, abs(a_CV2 * as))
    bs_btch <- bs + rnorm(nGenes, 0, abs(b_CV2 * bs))
    
    # transform params back
    ks_btch <- 2^ks_btch
    as_btch <- 2^as_btch
    bs_btch <- 2^bs_btch
    
    # add differential distributions
    if (i > nBatches) {
      ks_btch <- ks_btch * kM
      as_btch <- as_btch * aM
      bs_btch <- bs_btch * bM
    }
    
    # calculate dropout rates
    mus <- ks_btch * as_btch / (as_btch + bs_btch)
    if (Km < 0 || n < 0) {
      pZeros <- 0
    } else {
      pZeros <- 1 - mus^n / (Km + mus^n)
      pZeros <- pZeros + rnorm(nGenes, 0, pZero_SD)
      pZeros[pZeros > 1] <- 1
      pZeros[pZeros < 0] <- 0
    }
    pars[[i]] <- list(ks = ks_btch, as = as_btch, bs = bs_btch, pZeros = pZeros)
  }
  
  # generate batch sample sizes
  # nbinom/binom/pois dep on mu and sigma
  nCells_var <- (nCells_CV * nCells)^2
  if (nCells_var > nCells) {
    s <- nCells^2 / (nCells_var - nCells)
    smplSizes <- rnbinom(nSmpls, mu = nCells, size = s)
  } else if (nCells_var < nCells) {
    p <- 1 - nCells_var / nCells
    s <- nCells / p
    smplSizes <- rbinom(nSmpls, size = s, prob = p)
  } else {
    smplSizes <- rpois(nSmpls, lambda = nCells)
  }
  
  # generate readcounts
  batches <- lapply(1:nSmpls, function(x){
    PoiBeta(pars[[x]]$ks, pars[[x]]$as, pars[[x]]$bs,
            n = smplSizes[x], d = pars[[x]]$pZeros)
  })
  ds <- do.call(cbind, batches)
  
  # create pData
  btchLabs <- paste0("Batch", 1:nSmpls)
  pData <- data.table::data.table(
    sample = paste0("Cell", 1:ncol(ds)),
    batch = unlist(
      lapply(1:nSmpls,function(x) rep(btchLabs[x], smplSizes[x]))
    ),
    group = c( rep("Group1", sum(smplSizes[1:nBatches])),
               rep("Group2", sum(smplSizes[-nBatches:-1])) )
  )
  
  names(pars) <- btchLabs
  return(list(pData = pData, fData = fData, table = ds, params = pars))
}