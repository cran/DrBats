# Aim: simulate data to test a Bayesian Latent Factor Model
# Persons : Gabrielle Weinrott [cre, aut]

# Function to generate a set of signals
#
# N : number of individuals
# P : number of observations
# b.range : range of location of first gaussian
# c.range : range of location for second gaussian
# b.sd : standard dev of first gaussian
# c.sd : standard dev of second gaussian
# a.range : range of values for slope
# y.range : range of loaction for y-intercept
# amp : amplitude of cosinus function
# per : period of cosinus function
# data.type : type of data to generate (options: sparse, sparse.tend, sparse.tend.cos)
#
gen_X <- function(N = 10,
                  P = 150,
                  b.range = c(0.2, 0.35),
                  c.range = c(0.65, 0.8),
                  b.sd = 2,
                  c.sd = 2,
                  a.range = c(-0.4, 0.4),
                  y.range = c(0, 10),
                  amp = 10,
                  per = 12,
                  data.type = "sparse"){

  theta =  c(0, 10, 1) # parameters for the O-U process

  t = 0:P
  X <- matrix(nrow = N, ncol = length(t))

  amp.sample_b = c(-12:-8, 8:12)
  amp.sample_c = c(-8:-4, 4:8)

  if(data.type == "sparse"){
    for(i in 1:N){
      b_i <- stats::runif(1, min = b.range[1]*max(t), max = b.range[2]*max(t)) # location of first peak
      c_i <- stats::runif(1, min = c.range[1]*max(t), max = c.range[2]*max(t)) # location of second peak
      gauss_b <- sample(amp.sample_b, 1)*stats::dnorm(t, mean = b_i, sd = b.sd)
      gauss_c <- sample(amp.sample_c, 1)*stats::dnorm(t, mean = c_i, sd = c.sd)
      X[i, ] <- gauss_b*t+ gauss_c*t+ sde::rsOU(n = t, theta)
    }
  }

  else if(data.type == "sparse.tend"){
    for(i in 1:N){
      y_i <- stats::runif(1, min = y.range[1], max = y.range[2]) # intercept
      a_i <- stats::runif(1, min = a.range[1], max = a.range[2]) # slope
      b_i <- stats::runif(1, min = b.range[1]*max(t), max = b.range[2]*max(t)) # location of first peak
      c_i <- stats::runif(1, min = c.range[1]*max(t), max = c.range[2]*max(t)) # location of first peak
      gauss_b <- sample(amp.sample_b, 1)*stats::dnorm(t, mean = b_i, sd = b.sd)
      gauss_c <- sample(amp.sample_c, 1)*stats::dnorm(t, mean = c_i, sd = c.sd)
      X[i, ] <- y_i + a_i*t + gauss_b*t+ gauss_c*t + sde::rsOU(n = t, theta)
    }
  }

  else if(data.type == "sparse.tend.cos"){
    for(i in 1:N){
      y_i <- stats::runif(1, min = y.range[1], max = y.range[2]) # intercept
      a_i <- stats::runif(1, min = a.range[1], max = a.range[2]) # slope
      b_i <- stats::runif(1, min = b.range[1]*max(t), max = b.range[2]*max(t)) # location of first peak
      c_i <- stats::runif(1, min = c.range[1]*max(t), max = c.range[2]*max(t)) # location of first peak
      gauss_b <- sample(amp.sample_b, 1)*stats::dnorm(t, mean = b_i, sd = b.sd)
      gauss_c <- sample(amp.sample_c, 1)*stats::dnorm(t, mean = c_i, sd = c.sd)
      X[i, ] <- y_i + a_i*t + gauss_b*t+ gauss_c*t + amp*cos(t/per) + sde::rsOU(n = t, theta)
    }
  }

  return(X = X[ , -1])
}


# Simulate a matrix of observation times
#
# N : number of individuals
# P : number of observations
# t.range : range of times for the P observations
# b.range : range of location of first gaussian
# c.range : range of location for second gaussian
#
simul.matrix.t <- function(N = 10, P = 150, t.range = c(0, 1000),
                           b.range = c(0.2, 0.4), c.range = c(0.6, 0.8)){

  t <- data.frame()

  for(i in 1:N){
    t.1 <- sort(round(stats::runif(3*P/10, t.range[1], b.range[1]*t.range[2]), 2), decreasing = FALSE)
    t.b <- sort(round(stats::runif(P/20, b.range[1]*t.range[2], b.range[2]*t.range[2]), 2), decreasing = FALSE)
    t.12 <- sort(round(stats::runif(3*P/10, b.range[2]*t.range[2], c.range[1]*t.range[2]), 2), decreasing = FALSE)
    t.c <- sort(round(stats::runif(P/20, c.range[1]*t.range[2], c.range[2]*t.range[2]), 2), decreasing = FALSE)

    # correct for the number of times simulated depending on P
    if(P/20 != floor(P/20)){
      t.2 <- sort(round(stats::runif(3*P/10 + 1, c.range[2]*t.range[2], t.range[2]), 2), decreasing = FALSE)
    }
    else{
      t.2 <- sort(round(stats::runif(3*P/10, c.range[2]*t.range[2], t.range[2]), 2), decreasing = FALSE)
    }
    ti <- c(t.1, t.b, t.12, t.c, t.2)
    t <- rbind(t, ti)
  }

  return(t)
}

##' PCA data projected onto a histogram basis
##'
##' @param X the data matrix
##' @param t the matrix of observation times
##' @param t.range a vector specifying the observation time range (default : c(0, 1000))
##' @param breaks the number of breaks in the histogram basis (default : breaks = 15)
##' @return Xt.proj a matrix of projected observations
##' @return  U a matrix of eigenvectors
##' @return lambda a vector of eigenvalues
##' @return lambda.perc the percentage of inertia captured by each axis
##'
##' @examples
##' res <- drbats.simul(N = 5, P = 100, t.range = c(5, 100), breaks = 8)
##' pca.proj.Xt(res$X, res$t.simul, t.range = c(0, 100), breaks = 8)
##'
##' @export
##'
##' @author Gabrielle Weinrott
##'
pca.proj.Xt <- function(X,
                        t,
                        t.range = c(0, 1000),
                        breaks = 15){

  Xt.proj <- histoProj(X, t, t.range, breaks)
  X.histo <- Xt.proj$X.proj

  YC <- scale(X.histo, center = TRUE, scale = FALSE)
  cov.YC <- stats::cov(YC)
  eg <- eigen(cov.YC, symmetric = T)
  lambda <- eg$values
  lambda.perc <- lambda/sum(lambda)
  U <- eg$vectors

  res <- list(Xt.proj = Xt.proj,
              U = U,
              lambda = lambda,
              lambda.perc = lambda.perc)
  return(res)
}

##' Build and decompose a low-rank matrix W
##'
##' @description Build and decompose a low-rank matrix from
##' a matrix of eigenvectors and eigenvalues
##' from principal component analysis
##'
##'@param U a matrix of eigenvectors
##'@param lambda a vector of corresponding eigenvalues
##'
##'@return W a low-rank matrix
##'@return D the number of latent factors
##'@return Q the orthogonal matrix of the W = QR matrix decomposition
##'@return R the upper triangular matrix of the W = QR matrix decomposition
##'
##'@examples
##' res <- drbats.simul(N = 5, P = 100, t.range = c(5, 100), breaks = 8)
##' res.pca <- pca.Deville(res$X, res$t.simul, t.range = c(5, 100), breaks = 8)
##' Wres.pca <- W.QR(res.pca$U, res.pca$lambda)
##' Wres.pca
##'
##'@author Gabrielle Weinrott
##'@export
##'
W.QR <- function(U, lambda){
  if(!is.data.frame(U) && !is.matrix(U)){
    stop("U must be a matrix or data frame")
  }
  if(ncol(U) != nrow(U)){
    stop("U must have the same number of rows and columns")
  }
  if(!is.vector(lambda)){
    stop("lambda must be a vector")
  }
  if(length(lambda) != nrow(U)){
    stop("lambda must be of same length as each column of U")
  }

  U <- as.matrix(U)
  W <- matrix(nrow = nrow(U), ncol = ncol(U))
  biglambda = c()

  count = c()
  # get rid of tiny eigenvalues
  for(i in 1:length(lambda)){
    if(lambda[i]/sum(lambda) < 1/length(lambda)){
      count[i] = i
      lambda[i] = 0
    }
  }

  D <- min(which(lambda == 0)) - 1
  if(D < 2) D = 2

  biglambda <- as.matrix(lambda*diag(1, ncol(U)))

  W <- U%*%biglambda[ , 1:D]

  qrW <- Matrix::qr(t(W))
  Q <- qr.Q(qrW) # rotation matrix
  R <- qr.R(qrW) # upper triangular

  res <- list(W = W, D = D, Q = Q, R = R)
  return(res)
}

# Simulate a data matrix Y
#
# N : number of individuals
# W : a matrix of latent factors
# D : the number of latent factors (ncol(W))
# R : a rotation matrix
# sigma2 : variance of the error
#
simul.matrix.Y <- function(N = 10,
                           W,
                           D,
                           R,
                           sigma2 = 0.2){
  W.lower <- t(R)

  Y.simul <- matrix(nrow = N, ncol = nrow(W))
  beta.simul <- matrix(nrow = D, ncol = N)
  epsilon.simul <- matrix(nrow = N, ncol = nrow(W))

  for(i in 1:N){
    epsilon <- MASS::mvrnorm(nrow(W), 0, sigma2)
    beta <- MASS::mvrnorm(D, 0, 1)
    Y.simul.i <- t(W.lower%*%beta + epsilon)
    Y.simul[i, ] <- Y.simul.i

    beta.simul[ , i] <- beta
    epsilon.simul[i , ] <- epsilon
  }

  res = list(Y = Y.simul, beta = beta.simul, epsilon = epsilon.simul)
  return(res)
}

##' Main simulation function
##'
##' @param N integer number of functions to simulate (default = 10)
##' @param P a number of observation times (default = 150)
##' @param t.range a range of times in which to place the P observations (default = c(1, 1000))
##' @param b.range a vector giving the range of values for the mean of the first mode (default b.range = c(0.2, 0.4))
##' @param c.range a vector giving the range of values for the mean of the second mode (default c.range = c(0.6, 0.8))
##' @param b.sd the standard deviation for the first mode (default b.sd = 2)
##' @param c.sd the standard deviation for the second mode (default c.sd = 2)
##' @param a.range a vector giving the range of values for the slope (default a.range = c(-0.4, 0.4))
##' @param y.range a vector giving the range of values for the intercept (default y.range = c(0, 10))
##' @param amp the amplitude of the cosine function (default = 10)
##' @param per the periodicity of the cosine function (default = 12)
##' @param data.type string indicating type of functions (options :sparse, sparse.tend, sparse.tend.cos)
##' @param breaks number of breaks in the histogram basis
##' @param sigma2 the precision of the error terms (default = 0.2)
##' @param seed integer specification of a seed (default = NULL)
##'
##' @return Y.simul a list containing a matrix Y, a matrix beta, and a matrix epsilon
##' @return t.simul a matrix of simulated observation times
##' @return X the underlying signal to build the data, see DataSimulationandProjection vignette
##' @return proj.pca the outputs of the function pca.proj.Xt
##' @return wlu the outputs of the function W.QR
##'
##' @examples
##' res <- drbats.simul(N = 5, P = 100, t.range = c(5, 100), breaks = 8)
##' X <- res$X
##' t <- res$t.simul
##' # To plot the observations, ie the rows
##' matplot(t(t), t(X), type = 'l', xlab = "Time", ylab = "X")
##'
##' @author Gabrielle Weinrott
##' @export
##'
drbats.simul <- function(N = 10,
                         P = 150,
                         t.range = c(0, 1000),
                         b.range = c(0.2, 0.4),
                         c.range = c(0.6, 0.8),
                         b.sd = 2,
                         c.sd = 2,
                         a.range = c(-0.4, 0.4),
                         y.range = c(0, 10),
                         amp = 10,
                         per = 12,
                         data.type = "sparse",
                         breaks = 15,
                         sigma2 = 0.2,
                         seed = NULL){

  N <- as.integer(N)
  if(N <= 1){
    stop("N must be greater than 1")
  }
  P <- as.integer(P)
  if(P <= 0){
    stop("P must be a positive integer")
  }
  if(!is.vector(t.range) | length(t.range) != 2){
    stop("t.range must be a vector of length 2")
  }
  if(!is.vector(b.range) | length(b.range) != 2){
    stop("b.range must be a vector of length 2")
  }
  if(!is.vector(c.range) | length(c.range) != 2){
    stop("c.range must be a vector of length 2")
  }
  b.sd <- as.numeric(b.sd)
  c.sd <- as.numeric(c.sd)
  if(b.sd <= 0 | c.sd <= 0){
    stop("b.sd and c.sd must be positive values")
  }

  if(!is.vector(a.range) | length(a.range) != 2){
    stop("a.range must be a vector of length 2")
  }
  if(!is.vector(y.range) | length(y.range) != 2){
    stop("y.range must be a vector of length 2")
  }
  amp <- as.integer(amp)
  per <- as.integer(per)
  if(amp <= 0 | per <= 0){
    stop("amp and per must be positive integers")
  }
  if(data.type != "sparse" & data.type != "sparse.tend" &
     data.type != "sparse.tend.cos"){
    stop("invalid data.type")
  }
  breaks <- as.integer(breaks)
  if(breaks <= 0){
    stop("breaks must be a positive integer")
  }
  if(!is.numeric(sigma2) | sigma2 <= 0){
    stop("sigma2 must be a positive number")
  }

  set.seed(seed)

  X <- gen_X(N, P, b.range, c.range, b.sd, c.sd, a.range, y.range, amp, per, data.type)
  t <- as.matrix(simul.matrix.t(N, P, t.range, b.range, c.range))

  proj.pca <- pca.proj.Xt(X, t, t.range, breaks)
  wlu <- W.QR(proj.pca$U, proj.pca$lambda)

  Y.simul <- simul.matrix.Y(N, wlu$W, wlu$D, wlu$R, sigma2)

  res <- list(Y.simul = Y.simul, t.simul = t, X = X, proj.pca = proj.pca, wlu = wlu)

  return(res)
}

# Aim: project longitudinal data onto a histogram basis

# Use the rectangle method to calculate the area under a
# curve Xi(ti)
#
# Xi : a vector of observations
# ti : a vector of observation times corresponding to Xi
# windows : the windows of time on which to use the rectangle method
#
intRec <- function(Xi,
                   ti,
                   windows){

  step <- windows[2] - windows[1]

  area = NULL # a vector of length (length(windows) - 1)
  max = NULL
  min = NULL
  count = NULL

  for(i in 1:(length(windows) - 1)){
    mask = (ti >= windows[i]) & (ti <= windows[i+1])
    t.int <- ti[mask]
    X.int <- Xi[mask]
    count[i] = length(X.int)
    if(count[i] >= 1){
      area.i = sum(X.int)/count[i]
      max.i = max(X.int)
      min.i = min(X.int)
    }
    if(count[i] < 1){
      area.i = NA
      max.i = NA
      min.i = NA
    }
    area[i] = area.i
    max[i] = max.i
    min[i] = min.i
  }

  # if the first element is NA
  if(is.na(area[1])){
    if(!is.na(area[2])){
      area[1] <- area[2]
      max[1] = max[2]
      min[1] = min[2]
      warning("the first two intervals have been merged")
    }
    else if(!is.na(area[3])) {
      area[1] <- area[3]
      area[2] <- area[3]
      max[1] = max[3]
      max[2] = max[3]
      min[1] = min[3]
      min[2] = min[3]
      warning("the first three intervals have been merged")
    }
    else if(!is.na(area[4])) {
      area[1] <- area[4]
      area[2] <- area[4]
      area[3] <- area[4]
      max[1] = max[4]
      max[2] = max[4]
      max[3] = max[4]
      min[1] = min[4]
      min[2] = min[4]
      min[3] = min[4]
      warning("the first four intervals have been merged")
    }
    else if(!is.na(area[5])) {
      area[1] <- area[5]
      area[2] <- area[5]
      area[3] <- area[5]
      area[4] <- area[5]
      max[1] <- max[5]
      max[2] <- max[5]
      max[3] <- max[5]
      max[4] <- max[5]
      min[1] <- min[5]
      min[2] <- min[5]
      min[3] <- min[5]
      min[4] <- min[5]
      warning("the first four intervals have been merged")
    }
  }

  if(is.na(area[1])) stop("reduce number of breaks")

  # if there are NAs elsewhere in the vector
  for(i in 2:(length(area))){
    if(is.na(area[i]) && !is.na(area[i+1]) && !is.na(area[i-1])){
      area[i] = (area[i-1] + area[i+1])/2
      max[i] = (max[i-1] + max[i+1])/2
      min[i] = (min[i-1] + min[i+1])/2
    }
    if(is.na(area[i]) && is.na(area[i+1]) && !is.na(area[i-1]) && !is.na(area[i+2])){
      area[i] <- (area[i-1] + area[i+2])/2
      area[i+1] <- (area[i-1] + area[i+2])/2
      max[i] <- (max[i-1] + max[i+2])/2
      max[i+1] <- (max[i-1] + max[i+2])/2
      min[i] <- (min[i-1] + min[i+2])/2
      min[i+1] <- (min[i-1] + min[i+2])/2
    }
    else if(is.na(area[i]) && is.na(area[i+1]) && is.na(area[i+2])){
      warning("increase window size")
    }
  }

  # if the last element is NA
  if(is.na(area[length(area)])){
    if(!is.na(area[length(area)-1])){
      area[length(area)] <- area[length(area) - 1]
      max[length(max)] <- max[length(max) - 1]
      min[length(min)] <- min[length(min) - 1]
      warning("the last two intervals have been merged")
    }
    else if(is.na(area[length(area)-1]) && !is.na(length(area)[i-2])){
      area[length(area)] <- area[length(area) - 2]
      area[length(area) - 1] <- area[length(area) - 3]
      max[length(max)] <- max[length(max) - 2]
      max[length(max) - 1] <- max[length(max) - 3]
      min[length(min)] <- min[length(min) - 2]
      min[length(min) - 1] <- min[length(min) - 3]
      warning("the last three intervals have been merged")
    }
    else warning("increase window size")
  }

  if(sum(is.na(area)) > 0){
    stop("increase window size")
  }

  res <- list(area = area, count = count, max = max, min = min)

  return(res)
}


##' Project a set of curves onto a histogram basis
##'
##' @param X a matrix
##' @param t a matrix of observation times
##' @param t.range a range of times in which to place the P projections (default = c(0, 1000))
##' @param breaks the number of intervals in the histogram basis
##'
##' @return X.proj the matrix X after projection
##' @return X.count a matrix containing the number of observations used to build the projection onto the histogram basis
##' @return windows a vector containing the first time of each window of the histogram intervals
##' @return X.max the matrix of minimum values in each window
##' @return X.min the matrix of maximum values in each window
##'
##' @examples
##' res <- drbats.simul(N = 5, P = 100, t.range = c(5, 100), breaks = 8)
##' res.proj <- histoProj(res$X, res$t.simul, t.range = c(5, 100), breaks = 8)
##' res.proj
##'
##' @author Gabrielle Weinrott
##' @export
##'
histoProj <- function(X, t, t.range, breaks){

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

  windows <- seq(min(t.range), max(t.range), length.out = breaks)

  X.proj <- matrix(nrow = nrow(X), ncol = length(windows) - 1)
  X.count <- matrix(nrow = nrow(X), ncol = length(windows) - 1)
  X.max <- matrix(nrow = nrow(X), ncol = length(windows) - 1)
  X.min <- matrix(nrow = nrow(X), ncol = length(windows) - 1)

  for(i in 1:nrow(X)){
    integ.X <- intRec(X[i, ], t[i, ], windows)
    X.proj[i, ] <- integ.X$area
    X.count[i, ] <- integ.X$count
    X.max[i, ] <- integ.X$max
    X.min[i, ] <- integ.X$min
  }

  res <- list(X.proj = X.proj, X.count = X.count, windows = windows,
              X.max = X.max, X.min = X.min)

  return(res)
}
