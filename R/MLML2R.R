#' MLE estimates of 5-mC and 5-hmC levels.
#'
#' @param G no 5-hmC from TAB-conversion.
#' @param H 5-hmC from TAB-conversion.
#' @param L no 5-mC from oxBS-conversion.
#' @param M 5-mC from oxBS-conversion.
#' @param T 5-mC+5-hmC from BS-conversion.
#' @param U no methylation from BS-conversion.
#' @param tol convergent tolerance.
#' @return The returned value is a list with the following components.
#' @return \item{mC}{estimate for the proportion of methylation.}
#' @return \item{hmC}{estimate for the proportion of hydroxymethylation.}
#' @return \item{C}{estimate for the proportion of unmethylation.}

MLML2R <- function(
  G       = NULL,
  H       = NULL,
  L       = NULL,
  M       = NULL,
  T       = NULL,
  U       = NULL,
  tol = 0.001)
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
    ######### 3 métodos

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

  } else if (is.null(G) || is.null(H)) {
    ##### ox-seq + bs-seq

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

  } else if (is.null(M) || is.null(L)) {
    ##### tab-seq + bs-seq

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

  } else {
    ##### tab-seq + ox-seq

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
  }

  #pm <- round(pme, 4)
  #ph <- round(phe, 4)
  pm <- pme
  ph <- phe
  pu <- (1 - pm - ph)

    proportions <- list()
    proportions$mC <- pm
    proportions$hmC <- ph
    proportions$C <- pu

    row.names(proportions[[1]]) <- row.names(T)
    colnames(proportions[[1]]) <- colnames(T)
    row.names(proportions[[2]]) <- row.names(T)
    colnames(proportions[[2]]) <- colnames(T)
    row.names(proportions[[3]]) <- row.names(T)
    colnames(proportions[[3]]) <- colnames(T)

    return(proportions)
  }
