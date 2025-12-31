
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fkcentroids

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fkcentroids)](https://CRAN.R-project.org/package=fkcentroids)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/fkcentroids)](https://CRAN.R-project.org/package=fkcentroids)
<!-- badges: end -->

Functional K-Centroids Clustering Using Phase and Amplitude Components

## What is fkcentroids?

`fkcentroids` aims to conduct functional *k*-means and *k*-medians
clustering by jointly considering the phase and amplitude variations of
functions. The routines first synchronizes each function, thereby
separating phase and amplitude components. For the phase components, a
Bayes Hilbert space structure is imposed to enable effective analysis,
while for the amplitude components, the usual *L*<sup>2</sup> geometry
is employed. Routines for function synchronization, as well as
clustering procedures based on the extracted phase and amplitude
components, are provided.

## How can I get fkcentroids?

<!-- Version 0.4.0 of this package is available on [CRAN](https://cran.r-project.org/package=fkcentroids): -->
<!-- ``` r -->
<!-- install.packages("fkcentroids") -->
<!-- library(fkcentroids) -->
<!-- ``` -->

You can install the development version (version 0.0.0.9000) of
`fkcentroids` via:

``` r
devtools::install_github("seungwoo-stat/fkcentroids")
library(fkcentroids)
```

## How do I use it?

Using the Berkeley growth curve functions from the **fda** package, the
functions can be clustered using `fkmeans()` (or `fkmedians()`) after
decomposing them into their phase and amplitude components.

``` r
library(fkcentroids)
library(fda)
#> Loading required package: splines
#> Loading required package: fds
#> Loading required package: rainbow
#> Loading required package: MASS
#> Loading required package: pcaPP
#> Loading required package: RCurl
#> Loading required package: deSolve
#> 
#> Attaching package: 'fda'
#> The following object is masked from 'package:graphics':
#> 
#>     matplot
#> The following object is masked from 'package:datasets':
#> 
#>     gait
data(growth)
Ytilde <- cbind(growth$hgtm, growth$hgtf)
x <- growth$age
fkmeans(Ytilde, x, t = seq(0, 1, length.out = 100),
        sync_map = "auc", sync_args = 1, k = 2, nstart = 10)
#> Functional k-means clustering with 2 clusters of sizes 42, 51 
#> 
#> Clustering vector: 
#>  boy01  boy02  boy03  boy04  boy05  boy06  boy07  boy08  boy09  boy10  boy11 
#>      1      1      1      1      1      1      1      1      1      1      1 
#>  boy12  boy13  boy14  boy15  boy16  boy17  boy18  boy19  boy20  boy21  boy22 
#>      1      1      1      1      1      1      1      1      1      1      1 
#>  boy23  boy24  boy25  boy26  boy27  boy28  boy29  boy30  boy31  boy32  boy33 
#>      1      1      1      1      1      1      1      1      1      1      1 
#>  boy34  boy35  boy36  boy37  boy38  boy39 girl01 girl02 girl03 girl04 girl05 
#>      1      1      1      1      1      1      2      2      2      2      2 
#> girl06 girl07 girl08 girl09 girl10 girl11 girl12 girl13 girl14 girl15 girl16 
#>      2      2      2      2      2      2      2      2      2      2      2 
#> girl17 girl18 girl19 girl20 girl21 girl22 girl23 girl24 girl25 girl26 girl27 
#>      2      2      2      2      2      2      2      2      1      2      2 
#> girl28 girl29 girl30 girl31 girl32 girl33 girl34 girl35 girl36 girl37 girl38 
#>      2      2      2      2      2      2      2      2      2      2      2 
#> girl39 girl40 girl41 girl42 girl43 girl44 girl45 girl46 girl47 girl48 girl49 
#>      2      2      2      2      2      2      2      2      2      2      1 
#> girl50 girl51 girl52 girl53 girl54 
#>      2      1      2      2      2 
#> 
#> Within cluster sum of squares by cluster: 
#> [1] 27.38324 28.68932
#>  (between_SS / total_SS = 35.7 %)
#> 
#> Available components: 
#>  [1] "cluster"       "centers.Xclrv" "centers.Y"     "totss"        
#>  [5] "withinss"      "tot.withinss"  "betweenss"     "size"         
#>  [9] "iter"          "ifault"        "alpha0"
```

Refer to the package documentation for details.

``` r
?`fkcentroids-package`
```

## Where can I learn more?

Visit [this repo](https://github.com/seungwoo-stat/fkcentroids-paper)
for code to reproduce the figures and analysis from the paper Kang and
Oh (2026).

## References

Seungwoo Kang and Hee-Seok Oh. Multiview Representation and Clustering
of Functional Data. *Unpublished Manuscript*, 2026.
