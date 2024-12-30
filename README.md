---
editor_options: 
  markdown: 
    wrap: 72
---

# VBel

Variational Bayes for fast and accurate empirical likelihood inference

# About this package

This package allows you to run GVA on a data set in R and C++ for faster
computation (for 10,000 iterations of GVA: 40.23s for partial R and cpp
and 28.5s for purely cpp computation).

This package also allows you to run AEL on a data set in R and C++ for
faster computation (for 500 iterations of AEL: 0.2s for purely cpp and
0.1s for R and cpp with pre-z calculation).

------------------------------------------------------------------------

# Pre-installation instructions (Mac Users Only)

To install this package in Mac requires a Fortran compiler (through its
RcppEigen dependency). Chances are, your current Fortran compiler is not
up-to-date. To update your Fortran compiler, simply follow the steps
here: <br />  

1.  In your Mac App Store, search "Xcode" and install. <br />
2.  Open Terminal application. Type in

``` {eval="FALSE"}
xcode-select --install
```

     and follow the instructions.<br />      3. Click on the link
[here](https://github.com/fxcoudert/gfortran-for-macOS/releases).
Download the gfortan dmg file according to your MacOS version. <br />  
   4. Open the dmg file, run the gfortran installer, follow all the
instructions.

An alternative recommended method is to use the packet manager
[Homebrew](https://docs.brew.sh/Installation):

  1. Check if you have homebrew with

``` {eval="FALSE"}
$ brew doctor
```

     If you don't have it installed, use the following code from the
Homebrew webiste. Check the website that it hasn't changed since. It
will ask for your user password (you won't see characters as you type).
Follow the instructions.

``` {eval="FALSE"}
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

     2. Install GFortran using gcc (contains GFortran).

``` {eval="FALSE"}
brew install gcc
```

------------------------------------------------------------------------

# Installation

```{r}
# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("jlimrasc/VBel")
```

------------------------------------------------------------------------

# Toy example

```{r}
library(VBel)

# -----------------------------
# Generate toy variables
# -----------------------------
# Generating 30 data points from a simple linear-regression model
set.seed(1)
x    <- runif(30, min = -5, max = 5)
vari <- rnorm(30, mean = 0, sd = 1)
y    <- 0.75 - x + vari
lam0 <- matrix(c(0,0), nrow = 2)
th   <- matrix(c(0.8277, -1.0050), nrow = 2)
z    <- cbind(x, y)

# Specify moment condition functions for linear regression and its corresponding derivative
h    <- function(z, th) {
 xi     <- z[1]
 yi     <- z[2]
 h_zith <- c(yi - th[1] - th[2] * xi, xi*(yi - th[1] - th[2] * xi))
 matrix(h_zith, nrow = 2)
}

delthh <- function(z, th) {
 xi <- z[1]
 matrix(c(-1, -xi, -xi, -xi^2), 2, 2)
}

# Specify derivative of log prior
delth_logpi <- function(theta) { -0.0001 * mu0 }

# Specify AEL constant and Newton-Rhapson iteration
a         <- 0.00001
AEL_iters <- 500

# Specify initial values for GVA mean vector and Cholesky
reslm <- lm(y ~ x)
mu0   <- matrix(unname(reslm$coefficients),2,1)
C0    <- unname(t(chol(vcov(reslm))))

# Specify details for ADADELTA (Stochastic Gradient-Descent)
SGD_iters <- 10000
epsil     <- 10^-5
rho       <- 0.9

# -----------------------------
# Excecute functions
# -----------------------------
ansAELRcpp <- compute_AEL(th, h, lam0, a, z, AEL_iters)

resultGVA <-compute_GVA(mu, C0, h, delthh, delth_logpi, z, lam0, rho, epsil, a, 
SGD_iters, AEL_iters)

diagnostic_plot(resultGVA) # Plot the results to check for convergence
```

------------------------------------------------------------------------
