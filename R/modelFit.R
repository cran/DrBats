# Aim: Fit a Bayesian Latent Factor Model
# Persons : Gabrielle Weinrott [cre, aut]

##' Fit a Bayesian Latent Factor to a data set
##' using STAN
##'
##' @param model a string indicating the type of model ("igPLT", "cauchyPLT" or "igsparse", default = "PLT")
##' @param prog a string indicating the MCMC program to use (default = "stan")
##' @param parallel true or false, whether or not to parelleize (done using the package "parallel")
##' @param Xhisto matrix of simulated data (projected onto the histogram basis)
##' @param nchains number of chains (default = 2)
##' @param nthin the number of thinned interations (default = 1)
##' @param niter number of iterations (default = 1e4)
##' @param D number of latent factors
##'
##' @return stanfit, a STAN object
##'
##' @references The Stan Development Team Stan Modeling Language User's Guide and Reference Manual. http://mc-stan.org/
##' @author Gabrielle Weinrott
##'
##' @export
##'
modelFit <- function(model = "igPLT",
                     prog = "stan",
                     parallel = TRUE,
                     Xhisto = NULL,
                     nchains = 4,
                     nthin = 10,
                     niter = 10000,
                     D){

  if(is.null(Xhisto)){
    stop("No data specified!")
  }
  if(model != "igPLT" & model != "cauchyPLT" & model != "igsparse"){
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
  D <- as.integer(D)
  if(is.null(D)){
    stop("No number of latent factors specified")
  }
  if(D <= 0){
    stop("D must be a positive integer")
  }

  Xhisto <- scale(Xhisto, center = TRUE, scale = FALSE)

  rstan::rstan_options(auto_write = TRUE)

  if(parallel == TRUE){
    options(mc.cores = parallel::detectCores())
  }

  N <- dim(Xhisto)[1]
  P <- dim(Xhisto)[2]

  if(D >= P)
    stop("D must be smaller than ncol(Xhisto)")
  Q <- P*D-(D*(D-1)/2)


  if(model == "igPLT"){
    stan_data <- list(P=P, N=N, D=D, Q=Q, Xhisto = Xhisto)

    scode <- "data {
    int<lower=1> N; // observations
    int<lower=1> P;  // variables
    int<lower=1> D; // latent variables
    int<lower=1> Q;  // number of off-diagonal elements

    vector[P]  Xhisto[N]; // data matrix
    }

    parameters {
    vector[D] B[N]; // factor loadings
    vector[Q] offdiag;
    real<lower=0> sigma2;
    real<lower=0> tau2;
    }

    transformed parameters {
    matrix[P, D]  tR;
    {
    int index;
    for (j in 1:D) {
    index <- index + 1;
    tR[j,j] <- offdiag[index];
    for (i in (j+1):P) {
    index <- index + 1;
    tR[i,j] <- offdiag[index];
    }
    }
    for(i in 1:(D-1)){
    for(j in (i+1):D){
    tR[i,j] <- 0;
    }
    }
    }
    }

    model {
    offdiag ~ normal(0, tau2); // priors of the loadings
    tau2 ~ inv_gamma(0.001, 0.001);
    sigma2 ~ inv_gamma(0.001, 0.001);

    for (n in 1:N){
    B[n] ~ normal(0, 1); // factor constraints
    Xhisto[n] ~ normal(tR*B[n], sigma2); //the likelihood
    }
    }
    "
  }
  if(model == "cauchyPLT"){
    stan_data <- list(P=P, N=N, D=D, Q=Q, Xhisto = Xhisto)

    scode <- "data {
    int<lower=1> N; // observations
    int<lower=1> P;  // variables
    int<lower=1> D; // latent variables
    int<lower=1> Q;  // number of off-diagonal elements

    vector[P]  Xhisto[N]; // data matrix
    }

    parameters {
    vector[D] B[N]; // factors
    vector[Q] offdiag;
    real<lower=0> sigma;
    real<lower=0> tau;
    }

    transformed parameters {
    matrix[P, D]  tR;
    {
    int index;
    for (j in 1:D) {
    index <- index + 1;
    tR[j,j] <- offdiag[index];
    for (i in (j+1):P) {
    index <- index + 1;
    tR[i,j] <- offdiag[index];
    }
    }
    for(i in 1:(D-1)){
    for(j in (i+1):D){
    tR[i,j] <- 0;
    }
    }
    }
    }

    model {
    offdiag ~ normal(0, tau^2); // priors of the loadings
    tau ~ cauchy(0, 5);
    sigma ~ cauchy(0, 5);

    for (n in 1:N){
    B[n] ~ normal(0, 1); // factor constraints
    Xhisto[n] ~ normal(tR*B[n], sigma^2); //the likelihood
    }
    }
    "
  }

  if(model == "igsparse"){
    stan_data <- list(P=P, N=N, D=D, Xhisto = Xhisto)

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

  stanfit <- rstan::stan(model_code = scode, data = stan_data, chains = nchains,
                         thin = nthin, iter = niter, refresh = -1)

  return(stanfit)
}


##' Main function to fit the Bayesian Latent Factor model
##'
##' @param params a list of parameters for the modelFit function, see modelFit help for details
##' @param out.file the path to the file that will be saved
##' @param verbose binary variable (default = 1)
##' @param obj.version output version (eg. "1.1.plt")
##'
##' @return nothing
##'
##' @author Gabrielle Weinrott
##'
##' @keywords internal
main.modelFit <- function(params, out.file, verbose = 1, obj.version){

  if(is.null(out.file)){
    stop("out.file not specified")
  }

  if(is.null(obj.version)){
    obj.version = 0.1
  }

  if(verbose > 0){
    write("fitting the STAN model, this might take a minute...", stdout())
  }

  stanfit <- modelFit(params$model,
                  params$prog,
                  params$parallel,
                  params$Xhisto,
                  params$nchains,
                  params$nthin,
                  params$niter,
                  params$D)

  if(verbose > 0){
    write(paste(Sys.time(), "model fit for the", params$model, "model, with ",
                params$prog, ";", params$nchains, " chains, parallelizing = ", params$parallel,
                params$niter, " iterations, ",
                params$nthin, "thinned iterations"), stdout())
  }

  if(verbose > 0){
    write("writing an Rda file...", stdout())
  }

  save(stanfit, file = paste(out.file, "stanfit", obj.version, ".RData", sep = ""))

  if(verbose > 0){
    write("done!", stdout())
  }
}
