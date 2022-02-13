# DrBats

Feed longitudinal data into a Bayesian Latent Factor Model to obtain a low-rank representation. Parameters are estimated using a Hamiltonian Monte Carlo algorithm with STAN. See G. Weinrott, B. Fontez, N. Hilgert and S. Holmes, "Bayesian Latent Factor Model for Functional Data Analysis", Actes des JdS 2016.
  
# Installation

To install the **DrBats** package, the easiest is to install it directly from Gitlab. Open an R session and run the following commands:

```R
library(remotes) 
XXXX
```

# Usage

Once the package is installed on your computer, it can be loaded into a R session:

```R
library(DrBats)
help(package="DrBats")
```

# Citation

As a lot of time and effort were spent in creating the **DrBats** method, please cite it when using it for data analysis:

G. Weinrott, B. Fontez, N. Hilgert and S. Holmes, Bayesian Latent Factor Model for Functional Data Analysis", Actes des JdS 2016.

You should also cite the **DrBats** package:

```R
citation("DrBats")
```

See also citation() for citing R itself.
