# VBel

Variational Bayes for fast and accurate empirical likelihood (AEL) inference

# About this package

description

This package allows you to run AEL on a data set in R, and can use C++ for faster computation 
(1s for R, 0.2s for Rcpp and 0.1s for Rcpp with pre-z calculation).

* * *

# Pre-installation instructions (Mac Users Only)
To install this package in Mac requires a Fortran compiler (through its RcppEigen dependency).
Chances are, your current Fortran compiler is not up-to-date. To update your Fortran compiler, simply follow the steps here: <br />
&nbsp;

1. In your Mac App Store, search "Xcode" and install. <br />
2. Open Terminal application. Type in

```{eval=FALSE}
xcode-select --install
```
&nbsp; &nbsp;&nbsp;
and follow the instructions.<br />
&nbsp; &nbsp;&nbsp;
3. Click on the link [here](https://github.com/fxcoudert/gfortran-for-macOS/releases). Download the gfortan dmg file according to your MacOS version. <br />
&nbsp; &nbsp;&nbsp;
4. Open the dmg file, run the gfortran installer, follow all the instructions.

An alternative recommended method is to use the packet manager [Homebrew](https://docs.brew.sh/Installation):

&nbsp;
1. Check if you have homebrew with
```{eval=FALSE}
$ brew doctor
```
&nbsp; &nbsp;&nbsp;
If you don't have it installed, use the following code from the Homebrew webiste. Check the website that it hasn't changed since.
It will ask for your user password (you won't see characters as you type). Follow the instructions.
```{eval=FALSE}
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
&nbsp; &nbsp;&nbsp;
2. Install GFortran using gcc (contains GFortran).
```{eval=FALSE}
brew install gcc
```

* * *

# Installation
```{r}
# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("jlimrasc/VBel")
```

* * *

# Toy example
```{r}
library(VaDA)
set.seed(1)

# Set up variables
x    <- runif(30, min = -5, max = 5) # Random Uniform Distribution
elip <- rnorm(30, mean = 0, sd = 1)  # Random Normal Distribution
y    <- 0.75 - x + elip
lam0 <- matrix(c(0,0), nrow = 2)
th   <- matrix(c(0.8277, -1.0050), nrow = 2)
a <- 0.001
z    <- cbind(x, y)
h    <- function(z, th) {
    xi <- z[1]
    yi <- z[2]
    h_zith <- c(yi - th[1] - th[2] * xi)
    h_zith[2] <- xi * h_zith[1]
    matrix(h_zith, nrow = 2)
}

# Excecute function
compute_AEL_R(th, h, lam0, a, z)
```

* * *