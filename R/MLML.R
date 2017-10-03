#' MLE (maximum likelihood estimates) of 5-mC and 5-hmC levels.
#'
#' @param G Unmethylated channel (intensities/counts) from TAB-conversion (5-C + 5-mC).
#' @param H Methylated channel (intensities/counts) from TAB-conversion  (True 5-hmC).
#' @param L Unmethylated channel (intensities/counts) from oxBS-conversion (5-C + 5-hmC).
#' @param M Methylated channel (intensities/counts) from oxBS-conversion (True 5-mC).
#' @param T Methylated channel (intensities/counts) from standard BS-conversion (5-mC+5-hmC).
#' @param U Unmethylated channel (intensities/counts) from standard BS-conversion (True 5-C).
#' @param tol convergence tolerance.
#' @return The returned value is a list with the following components.
#' @return \item{mC}{maximum likelihood estimate for the proportion of methylation.}
#' @return \item{hmC}{maximum likelihood estimate for the proportion of hydroxymethylation.}
#' @return \item{C}{maximum likelihood estimate for the proportion of unmethylation.}
#' @return \item{methods}{the conversion methods used to produce the MLE}
#' @references
#'   Qu et al. MLML: consistent simultaneous estimates of DNA methylation and hydroxymethylation.
#'   Bioinformatics 2013, 29(20) pages 2645--2646.

MLML <- function(
  G       = NULL,
  H       = NULL,
  L       = NULL,
  M       = NULL,
  T       = NULL,
  U       = NULL,
  tol = 0.0001)
{

  g <- if (is.null(G)) {
    NULL
  } else {
    round(G)
  }
  h <- if (is.null(H)) {
    NULL
  } else {
    round(H)
  }
  l <- if (is.null(L)) {
    NULL
  } else {
    round(L)
  }
  m <- if (is.null(M)) {
    NULL
  } else {
    round(M)
  }
  u <- if (is.null(U)) {
    NULL
  } else {
    round(U)
  }
  t <- if (is.null(T)) {
    NULL
  } else {
    round(T)
  }

  pme <- 0.5
  phe <- 0.5
  diff <- 1

  pm <- matrix()
  ph <- matrix()

  if (!is.null(G) &
      !is.null(H) &
      !is.null(L) & !is.null(M) &
      !is.null(T) & !is.null(U)) {
    ######### 3 methods

    while (diff > tol) {
      pme_ant <- pme
      phe_ant <- phe
      k <- t * (pme / (pme + phe))
      j <- g * (pme / (1 - phe))
      pme <- (m + j + k) / (m + l + h + g + u + t)
      phe <- ((h - k + t) / (h + g + t + u - j - k)) * (1 - pme)
      diff_pme <- abs(max(pme-pme_ant))
      diff_phe <- abs(max(phe-phe_ant))
      diff <- abs(max(diff_pme,diff_phe))
    }

    methods <- c("TAB-conversion +  oxBS-conversion + standard BS-conversion")


  } else if (is.null(G) || is.null(H)) {
    ##### oxBS-seq + BS-seq

    while (diff > tol) {
      pme_ant <- pme
      phe_ant <- phe
      k <- t * (pme / (pme + phe))
      pme <- (m + k) / (m + l + u + t)
      phe <- ((-k + t) / (t + u - k)) * (1 - pme)
      diff_pme <- abs(max(pme-pme_ant))
      diff_phe <- abs(max(phe-phe_ant))
      diff <- abs(max(diff_pme,diff_phe))
    }
    
    methods <- c("oxBS-conversion + standard BS-conversion")


  } else if (is.null(M) || is.null(L)) {
    ##### TAB-seq + BS-seq

    while (diff > tol) {
      pme_ant <- pme
      phe_ant <- phe
      k <- t * (pme / (pme + phe))
      j <- g * (pme / (1 - phe))

      pme <- (j + k) / (g + h + t + u)
      phe <- ((h - k + t) / (g + h + t + u - k - j)) * (1 - pme)
      diff_pme <- abs(max(pme-pme_ant))
      diff_phe <- abs(max(phe-phe_ant))
      diff <- abs(max(diff_pme,diff_phe))
    }
    
    methods <- c("TAB-conversion + standard BS-conversion")

  } else {
    ##### TAB-seq + Ox-seq

    while (diff > tol) {
      pme_ant <- pme
      phe_ant <- phe
      j <- g * (pme / (1 - phe))

      pme <- (j + m) / (g + h + m + l)
      phe <- ((h) / (h + g - j)) * (1 - pme)
      diff_pme <- abs(max(pme-pme_ant))
      diff_phe <- abs(max(phe-phe_ant))
      diff <- abs(max(diff_pme,diff_phe))
    }

    methods <- c("TAB-conversion +  oxBS-conversion")

  }

  pm <- pme
  ph <- phe
  pu <- (1 - pm - ph)

    proportions <- list()
    proportions$mC <- pm
    proportions$hmC <- ph
    proportions$C <- pu
    proportions$methods <- methods

    row.names(proportions[[1]]) <- row.names(T)
    colnames(proportions[[1]]) <- colnames(T)
    row.names(proportions[[2]]) <- row.names(T)
    colnames(proportions[[2]]) <- colnames(T)
    row.names(proportions[[3]]) <- row.names(T)
    colnames(proportions[[3]]) <- colnames(T)

    return(proportions)
  }
