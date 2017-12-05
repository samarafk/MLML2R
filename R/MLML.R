#' MLE (maximum likelihood estimates) of 5-mC and 5-hmC levels.
#'
#' @param Gc Unmethylated channel (intensities/counts) from TAB-conversion (5-C + 5-mC).
#' @param Hc Methylated channel (intensities/counts) from TAB-conversion  (True 5-hmC).
#' @param Lc Unmethylated channel (intensities/counts) from oxBS-conversion (5-C + 5-hmC).
#' @param Mc Methylated channel (intensities/counts) from oxBS-conversion (True 5-mC).
#' @param Tc Methylated channel (intensities/counts) from standard BS-conversion (5-mC+5-hmC).
#' @param Uc Unmethylated channel (intensities/counts) from standard BS-conversion (True 5-C).
#' @param iterative logical. If iterative=TRUE EM-algorithm is used. For the combination of
#'  two methods, iterative=FALSE returns the exact constrained MLE using the
#'  pool-adjacent-violators algorithm (PAVA). When all three methods are combined,
#'  iterative=FALSE returns the constrained MLE using Lagrange multiplier.
#' @param tol convergence tolerance; considered only if iterative=TRUE
#' @details The function returns MLE estimates (binomial model assumed).
#'  The function assumes that the order of the rows and columns in the input matrices
#'  are consistent. In addition, all the input matrices must have the same dimension.
#'  Usually, rows represent CpG loci and columns are the samples.
#' @return The returned value is a list with the following components.
#' @return \item{mC}{maximum likelihood estimate for the proportion of methylation.}
#' @return \item{hmC}{maximum likelihood estimate for the proportion of hydroxymethylation.}
#' @return \item{C}{maximum likelihood estimate for the proportion of unmethylation.}
#' @return \item{methods}{the conversion methods used to produce the MLE}
#' @examples
#' # load the example datasets from BS, oxBS and TAB methods
#' data(MethylatedBS_sim)
#' data(MethylatedOxBS_sim)
#' data(UnMethylatedBS_sim)
#' data(UnMethylatedOxBS_sim)
#' data(MethylatedTAB_sim)
#' data(UnMethylatedTAB_sim)
#'
#' # obtain MLE via EM-algorithm for BS+oxBS:
#' results_em <- MLML(Tc = MethylatedBS_sim , Uc = UnMethylatedBS_sim,
#' Lc = UnMethylatedOxBS_sim, Mc = MethylatedOxBS_sim,iterative=TRUE,tol=0.0001)
#'
#' # obtain constrained exact MLE for BS+oxBS:
#' results_exact <- MLML(Tc = MethylatedBS_sim , Uc = UnMethylatedBS_sim,
#' Lc = UnMethylatedOxBS_sim, Mc = MethylatedOxBS_sim)
#'
#' # obtain MLE via EM-algorithm for BS+TAB:
#' results_em <- MLML(Tc = MethylatedBS_sim , Uc = UnMethylatedBS_sim,
#' Gc = UnMethylatedTAB_sim, Hc = MethylatedTAB_sim,tol=0.0001)
#'
#' # obtain constrained exact MLE for BS+TAB:
#' results_exact <- MLML(Tc = MethylatedBS_sim , Uc = UnMethylatedBS_sim,
#' Gc = UnMethylatedTAB_sim, Hc = MethylatedTAB_sim)
#'
#' # obtain MLE via EM-algorithm for oxBS+TAB:
#' results_em <- MLML(Lc = UnMethylatedOxBS_sim, Mc = MethylatedOxBS_sim,
#' Gc = UnMethylatedTAB_sim, Hc = MethylatedTAB_sim,iterative=TRUE,tol=0.0001)
#'
#' # obtain constrained exact MLE for oxBS+TAB:
#' results_exact <- MLML(Lc = UnMethylatedOxBS_sim, Mc = MethylatedOxBS_sim,
#' Gc = UnMethylatedTAB_sim, Hc = MethylatedTAB_sim)
#'
#' # obtain MLE via EM-algorithm for BS+oxBS+TAB:
#' results_em <- MLML(Tc = MethylatedBS_sim , Uc = UnMethylatedBS_sim,
#' Lc = UnMethylatedOxBS_sim, Mc = MethylatedOxBS_sim,
#' Gc = UnMethylatedTAB_sim, Hc = MethylatedTAB_sim,iterative=TRUE,tol=0.0001)
#'
#' #' # obtain MLE via Lagrange multiplier for BS+oxBS+TAB:
#' results_exact <- MLML(Tc = MethylatedBS_sim , Uc = UnMethylatedBS_sim,
#' Lc = UnMethylatedOxBS_sim, Mc = MethylatedOxBS_sim,
#' Gc = UnMethylatedTAB_sim, Hc = MethylatedTAB_sim)
#'
#' @author
#' Samara Kiihl samara@ime.unicamp.br;
#' Maria Jose Garrido;
#' Arce Domingo-Relloso;
#' Jose Bermudez;
#' Maria Tellez-Plaza.
#'
#' @references
#'   Qu J, Zhou M, Song Q, Hong EE, Smith AD. MLML: consistent simultaneous estimates of
#'   DNA methylation and hydroxymethylation. Bioinformatics. 2013;29(20):2645-2646.
#'   doi:10.1093/bioinformatics/btt459.
#'
#'   Ayer M, Brunk HD, Ewing GM, Reid WT, Silverman E. An Empirical Distribution Function
#'   for Sampling with Incomplete Information. Ann. Math. Statist. 1955, 26(4), 641–647.
#'   doi:10.1214/aoms/1177728423.
#'
#'   Zongli Xu, Jack A. Taylor, Yuet-Kin Leung, Shuk-Mei Ho, Liang Niu;
#'   oxBS-MLE: an efficient method to estimate 5-methylcytosine and 5-hydroxymethylcytosine
#'   in paired bisulfite and oxidative bisulfite treated DNA,
#'   Bioinformatics, 2016;32(23):3667–3669.

MLML <- function(Gc       = NULL,
                 Hc       = NULL,
                 Lc       = NULL,
                 Mc       = NULL,
                 Tc       = NULL,
                 Uc       = NULL,
                 iterative = FALSE,
                 tol=0.00001)
{
  g <- if (is.null(Gc)) {
    NULL
  } else {
    round(Gc)
  }
  h <- if (is.null(Hc)) {
    NULL
  } else {
    round(Hc)
  }
  l <- if (is.null(Lc)) {
    NULL
  } else {
    round(Lc)
  }
  m <- if (is.null(Mc)) {
    NULL
  } else {
    round(Mc)
  }
  u <- if (is.null(Uc)) {
    NULL
  } else {
    round(Uc)
  }
  t <- if (is.null(Tc)) {
    NULL
  } else {
    round(Tc)
  }

  four <-function(x1,x2,x3,x4){
    return(identical(x1,x2) & identical(x1,x3) & identical(x1,x4))
  }

  six <-function(x1,x2,x3,x4,x5,x6){
    return(identical(x1,x2) & identical(x1,x3) & identical(x1,x4) &
           identical(x1,x5) & identical(x1,x6))
  }

  pme <- 0.3
  phe <- 0.5
  diff <- 1

  pm <- matrix()
  ph <- matrix()

  if (!is.null(Gc) &
      !is.null(Hc) &
      !is.null(Lc) & !is.null(Mc) &
      !is.null(Tc) & !is.null(Uc)) {
    ######### 3 methods

      if (!six(rownames(Lc),rownames(Mc),rownames(Tc),
                        rownames(Uc),rownames(Gc),rownames(Hc)))
          {stop("Row names are inconsistent")}
    else {
      if (!iterative)
      {
        pm0 <- pm <- m/(m+l)
        ph0 <- ph <- h/(h+g)
        pc0 <- pc <- u/(u+t)
        Vm <- pm*(1-pm)/(m+l)
        Vh <- ph*(1-ph)/(h+g)
        Vc <- pc*(1-pc)/(u+t)
        lam <- (1-pm0-ph0-pc0)/(Vm+Vh+Vc)
        pm <- pm0 + lam*Vm
        ph <- ph0 + lam*Vh
        pc <- pc0 + lam*Vc
        Vm <- pm*(1-pm)/(m+l)
        Vh <- ph*(1-ph)/(h+g)
        Vc <- pc*(1-pc)/(u+t)
        lam <- (1-pm0-ph0-pc0)/(Vm+Vh+Vc)
        pme <- pm0 + lam*Vm
        phe <- ph0 + lam*Vh
      } else {

    while (diff > tol) {
      pme_ant <- pme
      phe_ant <- phe
      k <- t * (pme / (pme + phe))
      j <- g * (pme / (1 - phe))
      pme <- (m + j + k) / (m + l + h + g + u + t)
      phe <- ((h - k + t) / (h + g + t + u - j - k)) * (1 - pme)
      diff_pme <- abs(max(pme - pme_ant))
      diff_phe <- abs(max(phe - phe_ant))
      diff <- abs(max(diff_pme, diff_phe))
    }
      }
}

    methods <-
      c("TAB-conversion +  oxBS-conversion + standard BS-conversion")


  } else if (is.null(Gc) || is.null(Hc)) {
    ##### oxBS-seq + BS-seq

if (!four(rownames(Lc),rownames(Mc),rownames(Tc),
                        rownames(Uc))){stop("Row names are inconsistent")}
    else {

    if (!iterative)
    {
      pme = ifelse(t / (t + u) >= m / (m + l) , m / (m + l) ,
          (m + t) / (m + l + u + t))
      phe = ifelse(t / (t + u) >= m / (m + l), t / (t + u) - m / (m + l) , 0)
    } else {
      while (diff > tol) {
        pme_ant <- pme
        phe_ant <- phe
        k <- t * (pme / (pme + phe))
        pme <- (m + k) / (m + l + u + t)
        phe <- ((-k + t) / (t + u - k)) * (1 - pme)
        diff_pme <- abs(max(pme - pme_ant))
        diff_phe <- abs(max(phe - phe_ant))
        diff <- abs(max(diff_pme, diff_phe))
      }

    }

}

    methods <- c("oxBS-conversion + standard BS-conversion")


  } else if (is.null(Mc) || is.null(Lc)) {
    ##### TAB-seq + BS-seq

       if (!four(rownames(Gc),rownames(Hc),rownames(Tc),
                        rownames(Uc))){stop("Row names are inconsistent")}
    else {

    if (!iterative)
    {
      pme = ifelse(t / (t + u) >= h / (g + h) , t / (t + u) - h / (g + h) , 0)
      phe = ifelse(t / (t + u) >= h / (g + h) , h / (g + h),
          (h + t) / (g + h + t + u))
    } else {
      while (diff > tol) {
        pme_ant <- pme
        phe_ant <- phe
        k <- t * (pme / (pme + phe))
        j <- g * (pme / (1 - phe))

        pme <- (j + k) / (g + h + t + u)
        phe <- ((h - k + t) / (g + h + t + u - k - j)) * (1 - pme)
        diff_pme <- abs(max(pme - pme_ant))
        diff_phe <- abs(max(phe - phe_ant))
        diff <- abs(max(diff_pme, diff_phe))
      }
    }
}
    methods <- c("TAB-conversion + standard BS-conversion")


  } else {
    ##### TAB-seq + Ox-seq

      if (!four(rownames(Gc),rownames(Hc),rownames(Mc),
                        rownames(Lc))){stop("Row names are inconsistent")}
    else {

    if (!iterative)
    {
      pme = ifelse(m/(m+l)<=g/(h+g),m/(m+l),(g+m)/(g+h+m+l))
      phe = ifelse(m/(m+l)<=g/(h+g),h/(h+g),(h+l)/(g+h+m+l))
    } else {
      while (diff > tol) {
        pme_ant <- pme
        phe_ant <- phe
        j <- g * (pme / (1 - phe))

        pme <- (j + m) / (g + h + m + l)
        phe <- ((h) / (h + g - j)) * (1 - pme)
        diff_pme <- abs(max(pme - pme_ant))
        diff_phe <- abs(max(phe - phe_ant))
        diff <- abs(max(diff_pme, diff_phe))
      }
    }
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

  return(proportions)
}
