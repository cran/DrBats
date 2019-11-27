# Aim: Fit a Bayesian Latent Factor Model
# Persons : Gabrielle Weinrott [cre, aut]

##' Fit a Bayesian Latent Factor to a data set
##' using STAN
##'
##' @param model a string indicating the type of model ("PLT", or sparse", default = "PLT")
##' @param var.prior the family of priors to use for the variance parameters ("IG" for inverse gamma, or "cauchy")
##' @param prog a string indicating the MCMC program to use (default = "stan")
##' @param parallel true or false, whether or not to parelleize (done using the package "parallel")
##' @param Xhisto matrix of simulated data (projected onto the histogram basis)
##' @param nchains number of chains (default = 2)
##' @param nthin the number of thinned interations (default = 1)
##' @param niter number of iterations (default = 1e4)
##' @param R rotation matrix of the same dimension as the number of desired latent factors
##'
##' @return stanfit, a STAN object
##'
##' @references The Stan Development Team Stan Modeling Language User's Guide and Reference Manual. http://mc-stan.org/
##' @author Gabrielle Weinrott
##'
##'
##' @export
##' @import rstan
modelFit <- function(model = "PLT",
                     var.prior = "IG",
                     prog = "stan",
                     parallel = TRUE,
                     Xhisto = NULL,
                     nchains = 4,
                     nthin = 10,
                     niter = 10000,
                     R = NULL){

  if(is.null(Xhisto)){
    stop("No data specified!")
  }
  if(model != "PLT" & model != "sparse"){
    stop("Invalid model type")
  }
  if(prog != "stan"){
    warning("Invalid program type, defaulting to stan")
    prog = "stan"
  }
  if(!is.null(parallel) & parallel != TRUE & parallel != FALSE){
    parallel = FALSE
    warning("Invalid parallel input (must be TRUE or FALSE), defaulting to FALSE")
  }
  nchains <- as.integer(nchains)
  if(nchains <= 0){
    stop("Number of chains must be a positive integer")
  }
  nthin <- as.integer(nthin)
  if(nthin <= 0){
    stop("Number of thinning iterations must be a positive integer")
  }
  niter <- as.integer(niter)
  if(niter <= 0){
    stop("Number of iterations must be a positive integer")
  }
  if(is.null(R)){
    warning("No rotation matrix specified, using the identity matrix of dimension 3")
    R <- diag(1, 3)
  }
  if(var.prior != "IG" & var.prior != "cauchy"){
    stop("Invalid variance prior family, must select either IG or cauchy")
  }

  Xhisto <- scale(Xhisto, center = TRUE, scale = FALSE)

  rstan::rstan_options(auto_write = TRUE)

  if(parallel == TRUE){
    options(mc.cores = parallel::detectCores())
  }

  N <- dim(Xhisto)[1]
  P <- dim(Xhisto)[2]
  D <- nrow(R)

  if(D >= P)
    stop("D must be smaller than ncol(Xhisto)")
  Q <- P*D-(D*(D-1)/2)


  if(model == "PLT"){
    stan_data <- list(P=P, N=N, D=D, Q=Q, Xhisto = Xhisto, R=R)

    if(var.prior == "IG"){
      scode <- "data {
      int<lower=1> N; // observations
      int<lower=1> P;  // variables
      int<lower=1> D; // latent variables
      int<lower=1> Q;  // number of off-diagonal elements

      vector[P]  Xhisto[N]; // data matrix
      matrix[D, D] R; // rotation matrix
      }

      parameters {
      vector[D] B[N]; // factor loadings
      vector[Q] offdiag;
      real<lower=0> sigma2;
      real<lower=0> tau2;
      }

      transformed parameters {
      matrix[P, D]  tL;
      matrix[P, D] W;
      {
      int index;
      for (j in 1:D) {
      index <- index + 1;
      tL[j,j] <- offdiag[index];
      for (i in (j+1):P) {
      index <- index + 1;
      tL[i,j] <- offdiag[index];
      }
      }
      for(i in 1:(D-1)){
      for(j in (i+1):D){
      tL[i,j] <- 0;
      }
      }
      }
      W <- tL*R;
      }

      model {
      offdiag ~ normal(0, tau2); // priors of the loadings
      tau2 ~ inv_gamma(0.001, 0.001);
      sigma2 ~ inv_gamma(0.001, 0.001);

      for (n in 1:N){
      B[n] ~ normal(0, 1); // factor constraints
      Xhisto[n] ~ normal(W*B[n], sigma2); //the likelihood
      }
      }
      "
    }

    if(var.prior == "cauchy"){
      scode <- "data {
    int<lower=1> N; // observations
      int<lower=1> P;  // variables
      int<lower=1> D; // latent variables
      int<lower=1> Q;  // number of off-diagonal elements

      vector[P]  Xhisto[N]; // data matrix
      matrix[D, D] R; // rotation matrix
      }

      parameters {
      vector[D] B[N]; // factors
      vector[Q] offdiag;
      real<lower=0> sigma;
      real<lower=0> tau;
      }

      transformed parameters {
      matrix[P, D]  tL;
      matrix[P, D] W;
      {
      int index;
      for (j in 1:D) {
      index <- index + 1;
      tL[j,j] <- offdiag[index];
      for (i in (j+1):P) {
      index <- index + 1;
      tL[i,j] <- offdiag[index];
      }
      }
      for(i in 1:(D-1)){
      for(j in (i+1):D){
      tL[i,j] <- 0;
      }
      }
      }
      W <- tL*R ;
      }

      model {
      offdiag ~ normal(0, tau^2); // priors of the loadings
      tau ~ cauchy(0, 5);
      sigma ~ cauchy(0, 5);

      for (n in 1:N){
      B[n] ~ normal(0, 1); // factor constraints
      Xhisto[n] ~ normal(W*B[n], sigma^2); //the likelihood
      }
      }
      "
    }
    }


  if(model == "sparse"){
    stan_data <- list(P=P, N=N, D=D, Xhisto = Xhisto)

    if(var.prior == "IG"){
      scode <- "data {
      int<lower=1> N; // observations
      int<lower=1> P;  // variables
      int<lower=1> D; // latent variables

      vector[P]  Xhisto[N]; // data matrix
      }

      parameters {
      vector[D] B[N]; // factor loadings
      matrix[P, D] W; // latent factors
      real<lower=0> sigma2;
      vector[D] tau2;
      }

      model {
      sigma2 ~ inv_gamma(0.001, 0.001);

      for(i in 1:D){
      tau2[i] ~ inv_gamma(0.001, 0.001);
      W[ ,i] ~ double_exponential(0, tau2[i]);
      }

      for (n in 1:N){
      B[n] ~ normal(0, 1); // factor constraints
      Xhisto[n] ~ normal(W*B[n], sigma2); //the likelihood
      }
      }
      "
    }

    if(var.prior == "cauchy"){
      scode <- "data {
      int<lower=1> N; // observations
      int<lower=1> P;  // variables
      int<lower=1> D; // latent variables

      vector[P]  Xhisto[N]; // data matrix
      }

      parameters {
      vector[D] B[N]; // factor loadings
      matrix[P, D] W; // latent factors
      real<lower=0> sigma;
      vector[D] tau;
      }

      model {
      sigma ~ cauchy(0, 5);

      for(i in 1:D){
      tau[i] ~ cauchy(0, 5);
      W[ ,i] ~ double_exponential(0, tau[i]);
      }

      for (n in 1:N){
      B[n] ~ normal(0, 1); // factor constraints
      Xhisto[n] ~ normal(W*B[n], sigma^2); //the likelihood
      }
      }
      "
    }

  }

  stanfit <- rstan::stan(model_code = scode, data = stan_data, chains = nchains,
                         thin = nthin, iter = niter)

  return(stanfit)
}
