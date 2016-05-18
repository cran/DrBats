# Aim: Apply dimension reduction methods to a matrix
# Persons : Gabrielle Weinrott [cre, aut]
# Date : 27/10/2015

##' Perform a PCA using Deville's method
##'
##' @param X a data matrix
##' @param t a matrix of observation times corresponding to X
##' @param t.range the range of observation times in vector form (ex. t.range = c(0, 1000))
##' @param breaks integer number of histogram windows
##'
##' @return X.histo the matrix projected onto the histogram basis
##' @return U.histo a matrix of eigenvectors in the histogram basis
##' @return Cp a matrix of principal components
##' @return lambda a vector of eigenvalues
##' @return perc.lambda a vector of the percentage of total inertia explained by each principal component
##'
##' @examples
##' res <- drbats.simul(N = 5, P = 100, t.range = c(5, 100), breaks = 8)
##' res.pca <- pca.Deville(res$X, res$t.simul, t.range = c(5, 100), breaks = 8)
##' res.pca
##'
##' @author Gabrielle Weinrott
##'
##' @references JC Deville, "Methodes statisiques et numeriques de l'analyse harmonique", Annales de l'INSEE, 1974.
##' @export
##'
pca.Deville <- function(X, t, t.range, breaks){
  X <- as.matrix(X)
  t <- as.matrix(t)

  if(nrow(X) != nrow(t) | ncol(X) != ncol(t)){
    stop("X and t dimensions must match")
  }

  if(!is.vector(t.range) | length(t.range) != 2){
    stop("t.range must be a vector of length 2")
  }
  breaks <- as.integer(breaks)
  if(breaks <= 0){
    stop("breaks must be a positive integer")
  }

  histo <- histoProj(X, t, t.range, breaks)
  X.histo <- histo$X.proj

  weight.C <- scale(X.histo, center = TRUE, scale = FALSE)

  step <- histo$windows[2] - histo$windows[1]

  cov.weight.C <- (breaks/step)*stats::cov(weight.C)

  eg <- eigen(cov.weight.C)
  U <- eg$vectors
  lambda <- eg$values
  perc.lambda <- eg$values/sum(eg$values)
  U.histo <- U*sqrt(step/breaks)

  Cp <- weight.C%*%U.histo

  res <- list(X.histo = X.histo,
              U.histo = U.histo,
              Cp = Cp,
              lambda = lambda,
              perc.lambda = perc.lambda)

  return(res)
}


##' Perform Interval-PCA using Deville's method
##'
##' @param X a data matrix
##' @param t a matrix of observation times corresponding to X
##' @param t.range the range of observation times in vector form (ex. t.range = c(0, 1000))
##' @param breaks integer number of histogram windows
##'
##' @return scores
##' @return scores.min
##' @return scores.max
##'
##' @author Gabrielle Weinrott
##' @references JC Deville, "Methodes statisiques et numeriques de l'analyse harmonique", Annales de l'INSEE, 1974.
##'
##' @export
##'
interval.pca.Deville <- function(X, t, t.range, breaks){

  X <- as.matrix(X)
  t <- as.matrix(t)

  if(nrow(X) != nrow(t) | ncol(X) != ncol(t)){
    stop("X and t dimensions must match")
  }

  if(!is.vector(t.range) | length(t.range) != 2){
    stop("t.range must be a vector of length 2")
  }
  breaks <- as.integer(breaks)
  if(breaks <= 0){
    stop("breaks must be a positive integer")
  }

  proj <- histoProj(X, t, t.range, breaks)

  Xmin <- proj$X.min
  Xmax <- proj$X.max
  X.histo <- proj$X.proj

  Mx <- apply(X.histo, 2, mean) # keep scaling parameters for min and max

  scale.X <- matrix(nrow = nrow(X.histo), ncol = ncol(X.histo))
  for(i in 1:ncol(X.histo)){
    scale.X[ , i] <- X.histo[ , i] - Mx[i]
  }

  step <- proj$windows[2] - proj$windows[1]
  cov.weight.X <- (breaks/step)*stats::cov(scale.X)

  eg <- eigen(cov.weight.X)
  U <- eg$vectors
  U.histo <- U*sqrt(step/breaks)
  scores <- scale.X%*%U.histo

  scale.Xmin <- matrix(nrow = nrow(Xmin), ncol = ncol(Xmin))
  scale.Xmax <- matrix(nrow = nrow(Xmin), ncol = ncol(Xmin))
  for(i in 1:ncol(Xmin)){
    scale.Xmin[ , i] <- Xmin[ , i] - Mx[i]
    scale.Xmax[ , i] <- Xmax[ , i] - Mx[i]
  }

  Cp.min <- scale.Xmin%*%U.histo
  Cp.max <- scale.Xmax%*%U.histo

  res <- list(scores = scores, scores.MIN = Cp.min, scores.MAX = Cp.max)
  return(res)
}


##' Perform a weighted PCA using Deville's method
##' on a data matrix X that we project
##' onto a histogram basis and weighted
##'
##' @param X a data matrix
##' @param t a matrix of observation times corresponding to X
##' @param t.range the range of observation times in vector form (ex. t.range = c(a, b))
##' @param breaks integer number of histogram windows
##' @param Qp a matrix of weights, if Qp = NULL the function specifies a diagonal weight matrix
##'
##' @return X.histo the matrix projected onto the histogram basis
##' @return U.histo a matrix of eigenvectors in the histogram basis
##' @return Cp a matrix of principal components
##' @return lambda a vector of eigenvalues
##' @return perc.lambda a vector of the percentage of total inertia explained by each principal component
##'
##' @author Gabrielle Weinrott

##' @examples
##' res <- drbats.simul(N = 5, P = 100, t.range = c(5, 100), breaks = 8)
##' res.weighted <- weighted.Deville(res$X, res$t.simul, t.range = c(5, 100), breaks = 8, Qp = NULL)
##' res.weighted
##'
##' @export
##'
weighted.Deville <- function(X, t, t.range, breaks, Qp = NULL){

  X <- as.matrix(X)
  t <- as.matrix(t)

  if(nrow(X) != nrow(t) | ncol(X) != ncol(t)){
    stop("X and t dimensions must match")
  }

  if(!is.vector(t.range) | length(t.range) != 2){
    stop("t.range must be a vector of length 2")
  }
  breaks <- as.integer(breaks)
  if(breaks <= 0){
    stop("breaks must be a positive integer")
  }
  if(!is.null(Qp)){
    if(ncol(Qp)!= nrow(Qp)){
    stop("Qp is not square")
    }
  }

  histo <- histoProj(X, t, t.range, breaks)
  X.histo <- histo$X.proj
  X.count <- histo$X.count

  if(is.null(Qp)){
    X.count.col <- apply(X.count, MARGIN = 2, FUN = sum)
    X.count.sum <- sum(X.count)
    Qp <- diag(breaks*X.count.col/X.count.sum)
  }

  if(ncol(Qp) != ncol(X.histo) | nrow(Qp) != ncol(X.histo)){
    stop("invalid matrix Qp!!")
  }

  weight.X.histo <- X.histo%*%Qp

  weight.C <- scale(weight.X.histo, center = TRUE, scale = FALSE)

  step <- histo$windows[2] - histo$windows[1]

  cov.weight.C <- (breaks/step)*stats::cov(weight.C)

  eg <- eigen(cov.weight.C)
  U <- eg$vectors
  lambda <- eg$values
  perc.lambda <- eg$values/sum(eg$values)
  U.histo <- U*sqrt(step/breaks)

  Cp <- weight.C%*%U.histo

  res <- list(X.histo = X.histo,
              U.histo = U.histo,
              Cp = Cp,
              lambda = lambda,
              perc.lambda = perc.lambda)

  return(res)
}


##' Perform Coinertia Analysis on the PCA
##' of the Weighted PCA and Deville's PCA
##'
##' @param X.histo the data matrix projected onto the histogram basis
##' @param Qp a matrix of weights, if Qp = NULL the function specifies a diagonal weight matrix
##' @param X a data matrix, if X.histo is NULL and needs to be built
##' @param t a matrix of observation times, if X.histo is NULL and needs to be built
##' @param t.range the range of observation times in vector form,
##' if X.histo is NULL and needs to be built (default: t.range = c(0, 1000))
##' @param breaks integer number of histogram windows
##'
##' @return co_weight the co-inertia object
##'
##' @author Gabrielle Weinrott
##'
##' @examples
##' res <- drbats.simul(N = 5, P = 100, t.range = c(5, 100), breaks = 8)
##' res.coinertia <- coinertia.drbats(X = res$X, t = res$t.simul, t.range = c(5, 100), breaks = 8)
##' res.coinertia
##'
##' @export
##'
coinertia.drbats <- function(X.histo = NULL, Qp = NULL, X = NULL, t = NULL,
                             t.range = c(0, 1000), breaks){

  if(is.null(X.histo) & is.null(X) & is.null(t) & is.null(t.range)){
    stop("No data specified")
  }

  if(!is.null(X.histo) & is.null(Qp)){
    stop("If you specify X.histo you must also specify Qp")
  }

  if(!is.null(X.histo) & !is.null(Qp)){
    if(ncol(Qp) != ncol(X.histo) | nrow(Qp) != ncol(X.histo)){
      stop("invalid Qp matrix !")
    }
  }

  if(is.null(X.histo) & is.null(Qp)){
    if(is.null(t.range)){
      stop("Must specify t.range")
    }
    if(!is.vector(t.range) | length(t.range) != 2){
      stop("t.range must be a vector of length 2")
    }
    breaks <- as.integer(breaks)
    if(breaks <= 0){
      stop("breaks must be a positive integer")
    }

    histo <- histoProj(X, t, t.range, breaks)
    X.histo <- histo$X.proj
    X.count <- histo$X.count
    if(is.null(Qp)){
      X.count.col <- apply(X.count, MARGIN = 2, FUN = sum)
      X.count.sum <- sum(X.count)
      Qp <- diag(breaks*X.count.col/X.count.sum)
    }
    w.histo <- data.frame(X.histo%*%Qp)
    X.histo <- data.frame(X.histo)
    noweight <- ade4::dudi.pca(X.histo, scale = FALSE, scannf = F, nf = 3)
    weight <- ade4::dudi.pca(w.histo, scale = FALSE, scannf = F, nf = 3)

    co_weight <- ade4::coinertia(noweight, weight, scan = F, nf = 2)
  }
  else{
    w.histo <- data.frame(X.histo%*%Qp)
    X.histo <- data.frame(X.histo)
    noweight <- ade4::dudi.pca(X.histo, scale = FALSE, scannf = F, nf = 3)
    weight <- ade4::dudi.pca(w.histo, scale = FALSE, scannf = F, nf = 3)

    co_weight <- ade4::coinertia(noweight, weight, scan = F, nf = 2)
  }

  return(co_weight = co_weight)
}

