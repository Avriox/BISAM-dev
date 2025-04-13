<div align="center">
  <img src="img/logo.png" alt="BISAM Logo" width="200"/>

  <h1>BISAM</h1>
  <p><strong>Bayesian Indicator Saturation Model in R</strong></p>

  <a href="https://github.com/Avriox/BISAM">
    <img src="https://img.shields.io/badge/Status-Under%20Development-orange.svg" alt="Development Status"/>
  </a>
  <a href="https://github.com/Avriox/BISAM/actions/workflows/R-CMD-check.yaml">
    <img src="https://github.com/Avriox/BISAM/actions/workflows/R-CMD-check.yaml/badge.svg?branch=candidate" alt="R CMD Check Status"/>
  </a>
  <a href="https://cran.r-project.org/package=BISAM">
    <img src="https://www.r-pkg.org/badges/version/BISAM" alt="CRAN Status"/>
  </a>
  <a href="LICENSE"> <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT"/> </a>
  <img src="https://img.shields.io/badge/Depends-Rcpp-blue.svg" alt="Depends Rcpp"/>

</div>

---

> **⚠️ Disclaimer:** BISAM is currently under **active development**. The API may change without notice, and the package is **not yet available on CRAN**. The most up-to-date development version is available in the `candidate` branch. Please use with caution in production environments.

---

## Overview

**BISAM** provides an R implementation of Bayesian Indicator Saturation Models. This statistical technique is designed to robustly identify and model multiple structural breaks, outliers, or other forms of parameter instability within regression models or time series data.

The core computations are accelerated using **Rcpp** for improved performance, making it suitable for analysing larger datasets or more complex models.

---

## 🚀 Installation

You can install BISAM using one of the following methods:

**1. From GitHub (Development Version):**

The latest development version (currently recommended) is available on the `candidate` branch. You'll need the `remotes` package.

```r
# Install remotes if you haven't already
# install.packages("remotes")

# Install BISAM from the candidate branch
remotes::install_github("Avriox/BISAM@candidate")

# Note: Once the 'main' branch is stable, you might use:
# remotes::install_github("Avriox/BISAM")