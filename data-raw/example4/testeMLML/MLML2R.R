MLML2 <- function(G.matrix       = NULL,
                 H.matrix       = NULL,
                 L.matrix       = NULL,
                 M.matrix       = NULL,
                 T.matrix       = NULL,
                 U.matrix       = NULL,
                 iterative = FALSE,
                 tol=0.00001)
{
  g <- if (is.null(G.matrix)) {
    NULL
  } else {
    round(G.matrix)
  }
  h <- if (is.null(H.matrix)) {
    NULL
  } else {
    round(H.matrix)
  }
  l <- if (is.null(L.matrix)) {
    NULL
  } else {
    round(L.matrix)
  }
  m <- if (is.null(M.matrix)) {
    NULL
  } else {
    round(M.matrix)
  }
  u <- if (is.null(U.matrix)) {
    NULL
  } else {
    round(U.matrix)
  }
  t <- if (is.null(T.matrix)) {
    NULL
  } else {
    round(T.matrix)
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
  
  if (!is.null(G.matrix) &
      !is.null(H.matrix) &
      !is.null(L.matrix) & !is.null(M.matrix) &
      !is.null(T.matrix) & !is.null(U.matrix)) {
    ######### 3 methods
    
    if (!six(rownames(L.matrix),rownames(M.matrix),rownames(T.matrix),
             rownames(U.matrix),rownames(G.matrix),rownames(H.matrix)))
    {stop("Row names are inconsistent")}
    else {
      if (!iterative)
      {
        pm0 <- pm <- m/(m+l)
        ph0 <- ph <- h/(h+g)
        pc0 <- pc <- u/(u+t)
        ss <- pm0+ph0+pc0+1e-08
        pm0 <- pm <- ifelse( (m == 0 | h == 0 | u == 0 | g==0 | l==0 | t==0)  & (m/(m+l) + h/(h+g) + u/(u+t)) >= 1, pm0/ss , pm0)
        ph0 <- ph <- ifelse( (m == 0 | h == 0 | u == 0 | g==0 | l==0 | t==0)  & (m/(m+l) + h/(h+g) + u/(u+t)) >= 1, ph0/ss , ph0)
        pc0 <- pc <- ifelse( (m == 0 | h == 0 | u == 0 | g==0 | l==0 | t==0)  & (m/(m+l) + h/(h+g) + u/(u+t)) >= 1, pc0/ss , pc0)
        
        Vm <- pm*(1-pm)/(m+l)
        Vh <- ph*(1-ph)/(h+g)
        Vc <- pc*(1-pc)/(u+t)
        lam <- (1-pm0-ph0-pc0)/(Vm+Vh+Vc+1e-08)
        pm <- pm0 + lam*Vm
        ph <- ph0 + lam*Vh
        pc <- pc0 + lam*Vc
        Vm <- pm*(1-pm)/(m+l)
        Vh <- ph*(1-ph)/(h+g)
        Vc <- pc*(1-pc)/(u+t)
        lam <- (1-pm0-ph0-pc0)/(Vm+Vh+Vc+1e-08)
        pme <- pm0 + lam*Vm
        phe <- ph0 + lam*Vh
      } else {
        
        while (diff > tol) {
          pme_ant <- pme
          phe_ant <- phe
          k <- t * (pme / (pme + phe + 1e-08))
          j <- g * (pme / (1 - phe + 1e-08))
          pme <- (m + j + k) / (m + l + h + g + u + t)
          phe <- ((h - k + t) / (h + g + t + u - j - k)) * (1 - pme)
          diff_pme <- abs(max(pme - pme_ant,na.rm=TRUE))
          diff_phe <- abs(max(phe - phe_ant,na.rm=TRUE))
          diff <- abs(max(diff_pme, diff_phe,na.rm=TRUE))
        }
      }
    }
    
    methods <-
      c("TAB-conversion +  oxBS-conversion + standard BS-conversion")
    
    
  } else if (is.null(G.matrix) || is.null(H.matrix)) {
    ##### oxBS-seq + BS-seq
    
    if (!four(rownames(L.matrix),rownames(M.matrix),rownames(T.matrix),
              rownames(U.matrix))){stop("Row names are inconsistent")}
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
          k <- t * (pme / (pme + phe + 1e-08))
          pme <- (m + k) / (m + l + u + t)
          phe <- ((-k + t) / (t + u - k)) * (1 - pme)
          diff_pme <- abs(max(pme - pme_ant,na.rm=TRUE))
          diff_phe <- abs(max(phe - phe_ant,na.rm=TRUE))
          diff <- abs(max(diff_pme, diff_phe,na.rm=TRUE))
        }
        
      }
      
    }
    
    methods <- c("oxBS-conversion + standard BS-conversion")
    
    
  } else if (is.null(M.matrix) || is.null(L.matrix)) {
    ##### TAB-seq + BS-seq
    
    if (!four(rownames(G.matrix),rownames(H.matrix),rownames(T.matrix),
              rownames(U.matrix))){stop("Row names are inconsistent")}
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
          k <- t * (pme / (pme + phe + 1e-08))
          j <- g * (pme / (1 - phe + 1e-08))
          
          pme <- (j + k) / (g + h + t + u)
          phe <- ((h - k + t) / (g + h + t + u - k - j)) * (1 - pme)
          diff_pme <- abs(max(pme - pme_ant,na.rm=TRUE))
          diff_phe <- abs(max(phe - phe_ant,na.rm=TRUE))
          diff <- abs(max(diff_pme, diff_phe,na.rm=TRUE))
        }
      }
    }
    methods <- c("TAB-conversion + standard BS-conversion")
    
    
  } else {
    ##### TAB-seq + Ox-seq
    
    if (!four(rownames(G.matrix),rownames(H.matrix),rownames(M.matrix),
              rownames(L.matrix))){stop("Row names are inconsistent")}
    else {
      
      if (!iterative)
      {
        pme = ifelse(m/(m+l)<=g/(h+g),m/(m+l),(g+m)/(g+h+m+l))
        phe = ifelse(m/(m+l)<=g/(h+g),h/(h+g),(h+l)/(g+h+m+l))
      } else {
        while (diff > tol) {
          pme_ant <- pme
          phe_ant <- phe
          j <- g * (pme / (1 - phe + 1e-08))
          
          pme <- (j + m) / (g + h + m + l)
          phe <- ((h) / (h + g - j)) * (1 - pme)
          diff_pme <- abs(max(pme - pme_ant,na.rm=TRUE))
          diff_phe <- abs(max(phe - phe_ant,na.rm=TRUE))
          diff <- abs(max(diff_pme, diff_phe,na.rm=TRUE))
        }
      }
    }
    methods <- c("TAB-conversion +  oxBS-conversion")
    
    
    
  }
  
  pm <- pme
  ph <- phe
  
  pm[is.nan(pm)] <- NA
  ph[is.nan(ph)] <- NA
  
  pu <- (1 - pm - ph)
  
  proportions <- list()
  proportions$mC <- pm
  proportions$hmC <- ph
  proportions$C <- pu
  proportions$methods <- methods
  
  return(proportions)
}