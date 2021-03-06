% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLML.R
\name{MLML}
\alias{MLML}
\title{MLE (maximum likelihood estimates) of 5-mC and 5-hmC levels.}
\usage{
MLML(U.matrix = NULL, T.matrix = NULL, G.matrix = NULL,
  H.matrix = NULL, L.matrix = NULL, M.matrix = NULL,
  iterative = FALSE, tol = 1e-05)
}
\arguments{
\item{U.matrix}{Converted cytosines (T counts or U channel) from standard BS-conversion (True 5-C).}

\item{T.matrix}{Unconverted cytosines (C counts or M channel) from standard BS-conversion (reflecting 5-mC+5-hmC).}

\item{G.matrix}{Converted cytosines (T counts or U channel) from TAB-conversion (reflecting 5-C + 5-mC).}

\item{H.matrix}{Unconverted cytosines (C counts or M channel) from TAB-conversion(reflecting True 5-hmC).}

\item{L.matrix}{Converted cytosines (T counts or U channel) from oxBS-conversion (reflecting 5-C + 5-hmC).}

\item{M.matrix}{Unconverted cytosines (C counts or M channel) from oxBS-conversion (reflecting True 5-mC).}

\item{iterative}{logical. If iterative=TRUE EM-algorithm is used. For the combination of
two methods, iterative=FALSE returns the exact constrained MLE using the
pool-adjacent-violators algorithm (PAVA). When all three methods are combined,
iterative=FALSE returns the constrained MLE using Lagrange multiplier.}

\item{tol}{convergence tolerance; considered only if iterative=TRUE}
}
\value{
The returned value is a list with the following components.

\item{mC}{maximum likelihood estimate for the proportion of methylation.}

\item{hmC}{maximum likelihood estimate for the proportion of hydroxymethylation.}

\item{C}{maximum likelihood estimate for the proportion of unmethylation.}

\item{methods}{the conversion methods used to produce the MLE}
}
\description{
MLE (maximum likelihood estimates) of 5-mC and 5-hmC levels.
}
\details{
The function returns MLE estimates (binomial model assumed).
 The function assumes that the order of the rows and columns in the input matrices
 are consistent. In addition, all the input matrices must have the same dimension.
 Usually, rows represent CpG loci and columns are the samples.
}
\examples{
# load the example datasets from BS, oxBS and TAB methods
data(C_BS_sim)
data(C_OxBS_sim)
data(T_BS_sim)
data(T_OxBS_sim)
data(C_TAB_sim)
data(T_TAB_sim)

# obtain MLE via EM-algorithm for BS+oxBS:
results_em <- MLML(T.matrix = C_BS_sim , U.matrix = T_BS_sim,
L.matrix = T_OxBS_sim, M.matrix = C_OxBS_sim,iterative=TRUE)

# obtain constrained exact MLE for BS+oxBS:
results_exact <- MLML(T.matrix = C_BS_sim , U.matrix = T_BS_sim,
L.matrix = T_OxBS_sim, M.matrix = C_OxBS_sim)

# obtain MLE via EM-algorithm for BS+TAB:
results_em <- MLML(T.matrix = C_BS_sim , U.matrix = T_BS_sim,
G.matrix = T_TAB_sim, H.matrix = C_TAB_sim,iterative=TRUE)

# obtain constrained exact MLE for BS+TAB:
results_exact <- MLML(T.matrix = C_BS_sim , U.matrix = T_BS_sim,
G.matrix = T_TAB_sim, H.matrix = C_TAB_sim)

# obtain MLE via EM-algorithm for oxBS+TAB:
results_em <- MLML(L.matrix = T_OxBS_sim, M.matrix = C_OxBS_sim,
G.matrix = T_TAB_sim, H.matrix = C_TAB_sim,iterative=TRUE)

# obtain constrained exact MLE for oxBS+TAB:
results_exact <- MLML(L.matrix = T_OxBS_sim, M.matrix = C_OxBS_sim,
G.matrix = T_TAB_sim, H.matrix = C_TAB_sim)

# obtain MLE via EM-algorithm for BS+oxBS+TAB:
results_em <- MLML(T.matrix = C_BS_sim , U.matrix = T_BS_sim,
L.matrix = T_OxBS_sim, M.matrix = C_OxBS_sim,
G.matrix = T_TAB_sim, H.matrix = C_TAB_sim,iterative=TRUE)

#' # obtain MLE via Lagrange multiplier for BS+oxBS+TAB:
results_exact <- MLML(T.matrix = C_BS_sim , U.matrix = T_BS_sim,
L.matrix = T_OxBS_sim, M.matrix = C_OxBS_sim,
G.matrix = T_TAB_sim, H.matrix = C_TAB_sim)


# Example of datasets with zero counts and missing values:

C_BS_sim[1,1] <- 0
C_OxBS_sim[1,1] <- 0
C_TAB_sim[1,1] <- 0
T_BS_sim[1,1] <- 0
T_OxBS_sim[1,1] <- 0
T_TAB_sim[1,1] <- 0

C_BS_sim[2,2] <- NA
C_OxBS_sim[2,2] <- NA
C_TAB_sim[2,2] <- NA
T_BS_sim[2,2] <- NA
T_OxBS_sim[2,2] <- NA
T_TAB_sim[2,2] <- NA



}
\references{
Kiihl SF, Martinez-Garrido MJ, Domingo-Relloso A, Bermudez J, Tellez-Plaza M. MLML2R: an
  R package for maximum likelihood estimation of DNA methylation and hydroxymethylation
  proportions. Statistical Applications in Genetics and Molecular Biology. 2019;18(1).
  doi:10.1515/sagmb-2018-0031.

  Qu J, Zhou M, Song Q, Hong EE, Smith AD. MLML: consistent simultaneous estimates of
  DNA methylation and hydroxymethylation. Bioinformatics. 2013;29(20):2645-2646.
  doi:10.1093/bioinformatics/btt459.

  Ayer M, Brunk HD, Ewing GM, Reid WT, Silverman E. An Empirical Distribution Function
  for Sampling with Incomplete Information. Ann. Math. Statist. 1955, 26(4), 641-647.
  doi:10.1214/aoms/1177728423.

  Zongli Xu, Jack A. Taylor, Yuet-Kin Leung, Shuk-Mei Ho, Liang Niu;
  oxBS-MLE: an efficient method to estimate 5-methylcytosine and 5-hydroxymethylcytosine
  in paired bisulfite and oxidative bisulfite treated DNA,
  Bioinformatics, 2016;32(23):3667-3669.
}
\author{
Samara Kiihl samara@ime.unicamp.br;
Maria Jose Garrido;
Arce Domingo-Relloso;
Jose Bermudez;
Maria Tellez-Plaza.
}
