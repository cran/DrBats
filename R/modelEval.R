# Aim: Evaluate a Bayesian Latent Factor Model
# Persons : Gabrielle Weinrott [cre, aut]
# Date: January 10th, 2016

##' Convert a STAN objet to MCMC list
##'
##'@param stanfit a STAN object
##'
##'@return codafit an mcmc.list
##'
##'@author Gabrielle Weinrott
##'
##'@examples
##'data(stanfit) # output of modelFit or main.modelFit
##'coda.fit <- coda.obj(stanfit)
##'head(coda.fit)
##'
##'@export
##'
coda.obj <- function(stanfit){
  if(is.null(stanfit)){
    stop("STAN object not specified")
  }
  codafit <- coda::mcmc.list(lapply(1:ncol(stanfit), function(x) coda::mcmc(as.array(stanfit)[,x,])))
  return(codafit)
}

theta <- function(v1, v2){
  theta <- acos( sum(v1*v2) / ( sqrt(sum(v1 * v1)) * sqrt(sum(v2 * v2)) ) )
  return(theta)
}


reflect <- function(N, P, D, coda.fit, rotation, real.W, real.B){

  Q <- (P*D)-(D*(D-1)/2)

  nchains <- length(coda.fit)
  niter <- nrow(coda.fit[[1]])

  W <- list()
  tR <- list()
  B <- list()
  for(k in 1:nchains){

    tR[[k]] <- matrix(coda.fit[[k]][ , (N*D+Q+3):(N*D+Q+3+P*D-1)], nrow = niter, ncol = P*D)
    B[[k]] <- matrix(coda.fit[[k]][ , 1:(N*D)], nrow = niter, ncol = N*D)

    Wk <- matrix(nrow = niter, ncol = P*D)
    Bk <- matrix(nrow = niter, ncol = N*D)

    for(i in 1:niter){
      # work on the kth chain, ith iteration
      tRki <- matrix(tR[[k]][i, ], nrow = P, ncol = D)
      Bki <- matrix(B[[k]][i, ], nrow = D, ncol = N)

      # in terms of W
      Wki <- tRki%*%t(rotation)

      for(j in 1:P){
        for(l in 1:D){
          if(sign(Wki[j, l]) != sign(real.W[j, l])){
            Wki[j, l] = -(Wki[j, l])
          }
        }
      }
      Wk[i, ] <- as.vector(Wki)

      for(l in 1:D){
        for(n in 1:N){
          if(sign(Bki[l, n]) != sign(real.B[l, n])){
            Bki[l, n] = -(Bki[l, n])
          }
        }
      }
      Bk[i, ] <- as.vector(Bki)
    }
    W[[k]] <- Wk
    B[[k]] <- Bk
  }
  res <- list(W = W, B = B)
  return(res)
}

##' Post-process an MCMC list with reflection issues
##'
##' @param N number of individuals
##' @param P number of variables
##' @param D number of latent factors
##' @param coda.fit an MCMC list
##' @param rotation a DxD rotation matrix
##' @param real.W a reference latent factor matrix
##' @param real.B a reference factor loadings matrix
##'
##' @return mc.simu a clean MCMC list corrected for reflection issues
##'
##' @author Gabrielle Weinrott
##'
##' @examples
##'data(toydata) # simulated data
##'data(stanfit) # output of modelFit or main.modelFit
##'coda.fit <- coda.obj(stanfit)
##'
##'data.simul <- toydata$Y.simul$Y
##'N = nrow(data.simul)
##'D = toydata$wlu$D
##'P = ncol(data.simul)
##' ## PCA in the histogram basis
##' obs <- toydata$X
##' times <- toydata$t
##' pca.data <- pca.Deville(obs, times, t.range = c(min(times), max(times)), breaks = 15)
##' ## Post-processing landmark information
##' rotation <- toydata$wlu$Q # rotation matrix
##' real.W <- toydata$wlu$W # PCA-determined latent factors
##' real.B <- t(pca.data$Cp[, 1:(toydata$wlu$D)]) # PCA-determined scores
##' codafit.clean <- clean.mcmc(N, P, D, coda.fit, rotation, real.W, real.B)
##' head(codafit.clean)
##'
##' @export
##'
clean.mcmc <- function(N, P, D, coda.fit, rotation, real.W, real.B){

  N <- as.integer(N)
  if(N <= 0){
    stop("N must be a positive integer")
  }
  P <- as.integer(P)
  if(P <= 0){
    stop("P must be a positive integer")
  }
  D <- as.integer(D)
  if(D <= 0){
    stop("D must be a positive integer")
  }
  if(is.null(coda.fit)){
    stop("No MCMC.list specified")
  }
  if(is.null(rotation)){
    stop("No rotation specified")
  }
  if(is.null(real.W)){
    stop("No reference latent factors specified")
  }
  if(is.null(real.B)){
    stop("No reference scores specified")
  }

  Q <- (P*D)-(D*(D-1)/2)

  nchains <- length(coda.fit)
  niter <- nrow(coda.fit[[1]])

  ref <- reflect(N, P, D, coda.fit, rotation, real.W, real.B)

  simu <- list()

  for(k in 1:nchains){
    B <- ref$B[[k]]
    B.names <- c()
    for(l in 1:N){
      for(n in 1:D){
        B.names.i <- paste("B_", n, ",", l, sep = "")
        B.names <- c(B.names, B.names.i)
      }
    }

    W <- ref$W[[k]]
    W.names <- c()
    for(l in 1:D){
      for(j in 1:P){
        W.names.i <- paste("W_", j, ",", l, sep = "")
        W.names <- c(W.names, W.names.i)
      }
    }

    sigma <- matrix(coda.fit[[k]][ , N*D+Q+1], nrow = niter, ncol = 1)
    tau <- matrix(coda.fit[[k]][ , N*D+Q+2], nrow = niter, ncol = 1)

    together <- cbind(B, W, sigma, tau)
    colnames(together) <- c(B.names, W.names, "sigma2", "tau2")

    together <- coda::mcmc(together)
    simu[[k]] <- together
  }

  mc.simu <- coda::mcmc.list(simu)

  return(mc.simu)
}

