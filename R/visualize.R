# Aim : calculate posterior densities and produce corresponding graphics
# Author : Gabrielle Weinrott
# Date : January 15th, 2016

##' Calculate the log likelihood of the model
##'
##' @param Y is an NxP matrix
##' @param W is a PxD matrix
##' @param B is a DxN matrix
##' @param sigma2 is a real number
##' @param tau2 is a real number
##' @param a the first parameter of the inverse gamma distribution
##' @param b the second parameter of the inverse gamma distribution
##'
##' @return loglik the likelihood
##'
##' @author Gabrielle Weinrott
##'
##' @keywords internal
##'
calc.loglik <- function(Y, W, B, sigma2, tau2, a = 0.001, b = 0.001){
  if(is.null(Y)){
    stop("No data matrix Y specified")
  }
  if(is.null(W)){
    stop("No latent factor matrix W specified")
  }
  if(is.null(B)){
    stop("No scores matrix B specified")
  }
  if(is.null(sigma2)){
    stop("No sigma2 specified")
  }
  if(is.null(tau2)){
    stop("No tau2 specified")
  }
  a <- as.numeric(a)
  b <- as.numeric(b)

  sum.Y = 0
  sum.B = 0
  for(i in 1:nrow(Y)){
    mult.Y <- t(Y[i, ] - W%*%B[ , i])%*%(Y[i, ] - W%*%B[ , i])
    sum.Y = sum.Y + mult.Y
    mult.B <- t(B[ , i])%*%(B[ , i])
    sum.B = sum.B + mult.B
  }
  YB.post <- -1/2*(sum.Y/sigma2 + sum.B)

  sum.w = 0
  for(j in 1:nrow(W)){
    for(k in 1:ncol(W)){
      mult.w <- W[j, k]^2
      sum.w = sum.w + mult.w
    }
  }
  W.post <- -1/2*(sum.w/tau2)
  IG <- -b*(tau2 + sigma2)/(tau2*sigma2)
  expo <- YB.post + W.post + IG
  coeff <- 2*log(b^a/gamma(a)) + log(sqrt(2*pi)) - (nrow(Y) + a + 1)*log(sigma2) + (nrow(W)*ncol(W) + a + 1)*log(tau2)
  loglik <- expo*coeff

  return(loglik)
}

##' Calculate the unnormalized posterior density of the model
##'
##' @param mcmc.output an mcmc list as produced by clean.mcmc
##' @param Y the data matrix
##' @param D the number of latent factors
##' @param chain the chain to plot (default = 1)
##'
##' @return post a vector containing the posterior density at each iteration##' @examples
##' @examples
##' require(DrBats)
##' data("toydata")
##' data("stanfit")
##' codafit <- coda.obj(stanfit)
##' Y <- toydata$Y.simul$Y
##' N = nrow(Y)
##' D = toydata$wlu$D
##' P = ncol(Y)
##' ## PCA in the histogram basis
##' obs <- toydata$X
##' times <- toydata$t
##' pca.data <- pca.Deville(obs, times, t.range = c(min(times), max(times)), breaks = 15)
##' ## Post-processing landmark information
##' rotation <- toydata$wlu$Q # rotation matrix
##' real.W <- toydata$wlu$W # PCA-determined latent factors
##' real.B <- t(pca.data$Cp[, 1:(toydata$wlu$D)]) # PCA-determined scores
##' mcmc.output <- clean.mcmc(N, P, D, codafit, rotation, real.W, real.B)
##' post <- postdens(mcmc.output, Y, D, chain = 1)
##'
##' ## plot the density
##' hist(post)
##'
##' @export
##' @author Gabrielle Weinrott
##'
postdens <- function(mcmc.output, Y, D, chain = 1){
  if(is.null(mcmc.output)){
    stop("No mcmc.output specified")
  }
  if(is.null(Y)){
    stop("No data matrix Y specified")
  }
  if(is.null(D)){
    stop("No number of latent factors D specified")
  }
  chain <- as.integer(chain)
  if(chain <= 0){
    stop("Number of chains must be a positive integer")
  }

  N <- nrow(Y)
  P <- ncol(Y)
  Q <- (P*D)-(D*(D-1)/2)
  post <- c()
  for(i in 1:nrow(mcmc.output[[chain]])){
    B <- matrix(mcmc.output[[chain]][i, 1:(N*D)], nrow = D, ncol = N)
    W <- matrix(mcmc.output[[chain]][i, (N*D+1):(N*D+P*D)], nrow = P, ncol = D)

    sigma2 <- mcmc.output[[chain]][i, N*D+P*D+1]
    tau2 <- mcmc.output[[chain]][i, N*D+P*D+2]
    post.i <- calc.loglik(Y, W, B, sigma2, tau2)
    post <- c(post, post.i)
  }
  return(post)
}

##' Format scores output for visualization
##'
##' @param mcmc.output an mcmc list as produced by clean.mcmc
##' @param Y the matrix of data
##' @param D the number of latent factors
##' @param chain the chain to use (default = 1)
##' @param axes the axes to use (default = c(1, 2))
##' @param quant a vector of quantiles to retain (default = NULL)
##'
##' @return mean.df are the MCMC estimates for the parmeters
##' @return points.df contains all of the estimates of the chain
##' @return contour.df contains the exterior points of the convex hull of the cloud of estimates
##'
##' @examples
##' require(DrBats)
##' data("toydata")
##' data("stanfit")
##' codafit <- coda.obj(stanfit)
##' Y <- toydata$Y.simul$Y
##' N = nrow(Y)
##' D = toydata$wlu$D
##' P = ncol(Y)
##' ## PCA in the histogram basis
##' obs <- toydata$X
##' times <- toydata$t
##' pca.data <- pca.Deville(obs, times, t.range = c(min(times), max(times)), breaks = 15)
##' ## Post-processing landmark information
##' rotation <- toydata$wlu$Q # rotation matrix
##' real.W <- toydata$wlu$W # PCA-determined latent factors
##' real.B <- t(pca.data$Cp[, 1:(toydata$wlu$D)]) # PCA-determined scores
##' mcmc.output <- clean.mcmc(N, P, D, codafit, rotation, real.W, real.B)
##' beta.res <- visbeta(mcmc.output, Y, D, chain = 1, axes = c(1, 2), quant = c(0.05, 0.95))
##'
##' ggplot2::ggplot() +
##' ggplot2::geom_path(data = beta.res$contour.df, ggplot2::aes(x = x, y = y, colour = ind)) +
##' ggplot2::geom_point(data = beta.res$mean.df, ggplot2::aes(x = x, y = y, colour = ind))
##'
##' @export
##'
##' @author Gabrielle Weinrott
##'
visbeta <- function(mcmc.output, Y, D, chain = 1, axes = c(1, 2), quant = NULL){
  if(is.null(mcmc.output)){
    stop("No mcmc.output specified")
  }
  if(is.null(Y)){
    stop("No data matrix Y specified")
  }
  if(is.null(D)){
    stop("No number of latent factors D specified")
  }
  chain <- as.integer(chain)
  if(chain <= 0){
    stop("Number of chains must be a positive integer")
  }
  if(!is.null(quant)){
    for(i in 1:length(quant)){
      if(quant[i] <= 0 | quant[i] >= 1)
        stop("quant must be a vector of probabilites")
    }
  }

  N = nrow(Y)
  P = ncol(Y)

  # order densities
  post <- postdens(mcmc.output, Y, D, chain)
  decreasing.post <- order(post)

  # calculate quantiles
  if(!is.null(quant)){
    q <- stats::quantile(post[decreasing.post], probs = quant, na.rm = T)
    indi <- which(post > q[1] & post < q[2])
  }
  if(is.null(quant)){
    indi = 1:nrow(mcmc.output[[chain]])
  }

  B.chain <- mcmc.output[[chain]][indi, 1:(N*D)]
  # for data frame formatting and gggplot2
  x.cols <- seq(axes[1], N*D, by = D)
  y.cols <- seq(axes[2], N*D, by = D)

  # mean estimates
  B.mean <- apply(B.chain, 2, mean)
  B.x <- c(); B.y = c(); B.ind = c()
  for(i in 1:length(x.cols)){
    B.x <- c(B.x, B.mean[x.cols[i]])
    B.y <- c(B.y, B.mean[y.cols[i]])
    B.ind <- c(B.ind, i)
  }
  mean.df <- data.frame(x = B.x, y = B.y, ind = B.ind)
  mean.df$ind <- as.factor(mean.df$ind)


  # cloud of points
  B.x <- c(); B.y = c(); B.ind = c()
  for(i in 1:length(x.cols)){
    B.x <- c(B.x, B.chain[ , x.cols[i]])
    B.y <- c(B.y, B.chain[ , y.cols[i]])
    B.ind <- c(B.ind, rep(i, nrow(B.chain)))
  }
  B.df <- data.frame(x = B.x, y = B.y, ind = B.ind)
  B.df$ind <- as.factor(B.df$ind)

  # contour
  contour.x = c() ; contour.y = c() ; contour.ind = c()
  for(i in 1:length(x.cols)){
    hpts.i <- grDevices::chull(B.chain[ , x.cols[i]:y.cols[i]])
    hpts.i <- c(hpts.i, hpts.i[1])
    contour.x <- c(contour.x, B.chain[hpts.i, x.cols[i]:y.cols[i]][ , 1])
    contour.y <- c(contour.y, B.chain[hpts.i, x.cols[i]:y.cols[i]][ , 2])
    contour.ind <- c(contour.ind, rep(i, length(hpts.i)))
  }

  contour.df <- data.frame(x = contour.x, y = contour.y, ind = contour.ind)
  contour.df$ind <- as.factor(contour.df$ind)

  res.beta = list(mean.df = mean.df, points.df = B.df, contour.df = contour.df)
  return(res.beta)

}

##' Plot the estimates for the latent factors
##'
##' @param mcmc.output an mcmc list as produced by clean.mcmc
##' @param Y the matrix of data
##' @param D the number of latent factors
##' @param chain the chain to plot (default = 1)
##' @param factors a vector indicating the factors to plot (default = c(1, 2))
##'
##' @author Gabrielle Weinrott
##'
##' @return res.W a data frame containing the estimates for the factors, and their lower
##' and upper bounds
##' @return Inertia the percentage of total inertia captured by each of the factors
##'
##'
##' @examples
##' require(DrBats)
##' data("toydata")
##' data("stanfit")
##' codafit <- coda.obj(stanfit)
##' Y <- toydata$Y.simul$Y
##' N = nrow(Y)
##' D = toydata$wlu$D
##' P = ncol(Y)
##' ## PCA in the histogram basis
##' obs <- toydata$X
##' times <- toydata$t
##' pca.data <- pca.Deville(obs, times, t.range = c(min(times), max(times)), breaks = 15)
##' ## Post-processing landmark information
##' rotation <- toydata$wlu$Q # rotation matrix
##' real.W <- toydata$wlu$W # PCA-determined latent factors
##' real.B <- t(pca.data$Cp[, 1:(toydata$wlu$D)]) # PCA-determined scores
##' mcmc.output <- clean.mcmc(N, P, D, codafit, rotation, real.W, real.B)
##' W.res <- visW(mcmc.output, Y, D, chain = 1, factors = c(1, 2))
##'
##' @export
##'
##'
visW <- function(mcmc.output, Y, D, chain = 1, factors = c(1, 2)){
  if(is.null(mcmc.output)){
    stop("No mcmc.output specified")
  }
  if(is.null(Y)){
    stop("No data matrix Y specified")
  }
  if(is.null(D)){
    stop("No number of latent factors D specified")
  }
  chain <- as.integer(chain)
  if(chain <= 0){
    stop("Number of chains must be a positive integer")
  }
  for(i in 1:length(factors)){
    if(factors[i] > D)
      stop("Invalid factors")
  }



  N = nrow(Y)
  P = ncol(Y)
  Q <- (P*D)-(D*(D-1)/2)
  if(length(mcmc.output) == 1){
    W.med <- matrix(summary(mcmc.output)[[2]][(N*D+1):(N*D+P*D), 3], nrow = P, ncol = D)
    W.lower <- matrix(summary(mcmc.output)[[2]][(N*D+1):(N*D+P*D), 1], nrow = P, ncol = D)
    W.upper <- matrix(summary(mcmc.output)[[2]][(N*D+1):(N*D+P*D), 5], nrow = P, ncol = D)
  }
  else{
    W.med <- matrix(summary(mcmc.output[[chain]])[[2]][(N*D+1):(N*D+P*D), 3], nrow = P, ncol = D)
    W.lower <- matrix(summary(mcmc.output[[chain]])[[2]][(N*D+1):(N*D+P*D), 1], nrow = P, ncol = D)
    W.upper <- matrix(summary(mcmc.output[[chain]])[[2]][(N*D+1):(N*D+P*D), 5], nrow = P, ncol = D)
  }

  lower = c() ; est = c() ; upper = c() ; fact = c() ; inertia = c()
  for(i in 1:length(factors)){
    lower = c(lower, W.lower[ , factors[i]])
    est = c(est, W.med[ , factors[i]])
    inertia = c(inertia, stats::var(W.med[ , factors[i]]))
    upper = c(upper, W.upper[ , factors[i]])
    fact = c(fact, rep(i, nrow(W.med)))
  }

  res.W <- data.frame(Lower.est = lower,
                      Estimation = est,
                      Upper.est = upper,
                      Factor = as.factor(fact))

  Inertia = inertia/sum(inertia)

  res <- list(res.W = res.W, Inertia = Inertia)

  return(res)
}
