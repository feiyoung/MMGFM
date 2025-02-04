# MMGFM
High-dimensional multi-study multi-modality covariate-augmented generalized factor model

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/MMGFM)](https://cran.r-project.org/package=MMGFM)
[![](https://cranlogs.r-pkg.org/badges/MMGFM?color=orange)](https://cran.r-project.org/package=MMGFM)
[![](https://cranlogs.r-pkg.org/badges/grand-total/MMGFM?color=orange)](https://cran.r-project.org/package=MMGFM)
<!-- badges: end -->


Latent factor models that integrate data from multiple sources/studies or modalities have garnered considerable attention across various disciplines. However, existing methods predominantly focus either on multi-study integration or multi-modality integration, rendering them insufficient for analyzing the diverse modalities measured across multiple studies. To address this limitation and cater to practical needs, we introduce a high-dimensional generalized factor model that seamlessly integrates  multi-modality data from multiple studies, while also accommodating additional covariates.  



Check out  [Package Website](https://feiyoung.github.io/MMGFM/index.html) for a more complete description of the methods and analyses. 

# Installation
"MMGFM" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}
## Method 1：
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/MMGFM")

## Method 2: install from CRAN
install.packages("MMGFM")
```



## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Simulated example1](https://feiyoung.github.io/MMGFM/articles/simu_MMGFM.html)
* [Simulated example2](https://feiyoung.github.io/MMGFM/articles/simu2_MMGFM.html)

## Simulated codes
For the codes in simulation study, check the `simu_code` directory of the repo.


## News

MMGFM version 1.1 released! (2024-09-17) 


