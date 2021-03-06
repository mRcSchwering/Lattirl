% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulateData.R
\name{SimulateData}
\alias{SimulateData}
\title{SimulateData}
\usage{
SimulateData(k = 4.2, k_CV = 0.25, k_CV2 = 0.08, a = 1.9, a_CV = 0.35,
  a_CV2 = 0.07, b = 3.5, b_CV = 0.4, b_CV2 = 0.09, Km = 4.5,
  n = 3.4, pZero_SD = 0.1, nBatches = 3, nGenes = 10000, nCells = 100,
  nCells_CV = 0.3, ddP = 0, ddM = 1)
}
\arguments{
\item{k, a, b}{\code{numeric} scalar for average BetaPoisson parameters
in log2}

\item{Km, n}{\code{numeric} parameter scalar defining Hill function
by which 0's (dropouts) will be introduced}

\item{pZero_SD}{\code{numeric} scalar for Std of introduced dropouts}

\item{nBatches}{\code{integer} scalar number of batches per group}

\item{nGenes}{\code{integer} scalar number of genes}

\item{nCells}{\code{integer} scalar average number of cells per batch}

\item{ddP}{\code{numeric} scalar for proportion of differentially
distributed genes}

\item{ddM}{\code{numeric} scalar controlling amount of differential
distribution (3 would be quite strong)}

\item{_CV}{\code{numeric} scalar for coefficient of variation of
a parameter with a batch (in log2)}

\item{_CV2}{\code{numeric} scalar for coefficient of variation of
average parameters between batches (in log2)}

\item{nCells_SD}{\code{numeric} scalar coefficient of deviation for
numbers of cells in batches}
}
\value{
list of generated count table, phenotypic- (cell-) and feature-
        (gene-) information tables, and a list which holds the Parameters
        for each batch and each gene.
}
\description{
Generate a simulated scRNA-seq dataset of 2 biological groups with several
batches per group, a defined variability within and between batches, and
a defined proportion and amount of differentially distributed genes.
}
\details{
Count tables, phenotypic- and feature-information tables are generated
to resemble a scRNA-seq experiment.
Count distributions are modelled using a 3 parameter Beta-Poisson
distributions: \emph{Pois(k * Beta(a, b))}.

For low average expression (genewise) additional 0's are introduced to better
resemble real dropout rates.
This is done according to a Hill equation, which defines the propability of
a zero (dropout) from the average expression by
\emph{P_zero = x^n / (Km + x^n)} where x is the average expression.

The defaults resemble the \emph{Tung} 2016 dataset.
In general, these parameters fit well to datasets which were produced using
UMI's and the fluidigm platform.
For raw readcounts the absolute count values and variabilities are higher.
Count tables produced using Drop-Seq usually have more cells and a much
larger dropout rate.
}
\examples{
ds <- ds <- SimulateData()
str(ds)

}
