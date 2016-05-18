## ---- echo = F-----------------------------------------------------------
rm(list=ls())

## ---- messages = F-------------------------------------------------------
require(DrBats)
mycol<-c("#ee204d", "#1f75fe", "#1cac78", "#ff7538", "#b4674d", "#926eae",
                 "#fce883", "#000000", "#78dbe2", "#6e5160", "#ff43a4")

data("toydata")
data("stanfit")

## ------------------------------------------------------------------------
codafit <- coda.obj(stanfit)

data.simul <- toydata$Y.simul$Y

N = nrow(data.simul)
D = toydata$wlu$D
P = ncol(data.simul)

## PCA in the histogram basis
obs <- toydata$X
times <- toydata$t
pca.data <- pca.Deville(obs, times, t.range = c(min(times), max(times)), breaks = 15)

## Post-processing landmark information
rotation <- toydata$wlu$Q # rotation matrix
real.W <- toydata$wlu$W # PCA-determined latent factors
real.B <- t(pca.data$Cp[, 1:(toydata$wlu$D)]) # PCA-determined scores

codafit.clean <- clean.mcmc(N, P, D, codafit, rotation, real.W, real.B)

## ------------------------------------------------------------------------
post <- postdens(codafit.clean, Y = toydata$Y.simul$Y, D = toydata$wlu$D, chain = 1)
hist(post, main = "Histogram of the posterior density", xlab = "Density")

## ------------------------------------------------------------------------
beta.res <- visbeta(codafit.clean, toydata$Y.simul$Y, toydata$wlu$D, chain = 1, axes = c(1, 2), quant = c(0.05, 0.95))

ggplot2::ggplot(beta.res$mean.df, ggplot2::aes(x = x, y = y, colour = ind)) +
  ggplot2::geom_point(ggplot2::aes(x = x, y = y, colour = ind)) +
  ggplot2::ggtitle("HMC estimate of the scores")

## ------------------------------------------------------------------------
ggplot2::ggplot() +
  ggplot2::geom_point(data = beta.res$points.df, ggplot2::aes(x = x, y = y, colour = ind)) +
  ggplot2::geom_point(data = beta.res$mean.df, ggplot2::aes(x = x, y = y, colour = ind)) +
  ggplot2::ggtitle("Cloud of HMC estimates of the scores")

## ------------------------------------------------------------------------
ggplot2::ggplot() +
  ggplot2::geom_path(data = beta.res$contour.df, ggplot2::aes(x = x, y = y, colour = ind)) +
  ggplot2::geom_point(data = beta.res$mean.df, ggplot2::aes(x = x, y = y, colour = ind)) +
  ggplot2::ggtitle("Convex hull of HMC estimates of the scores")

