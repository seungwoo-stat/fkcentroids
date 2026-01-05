
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
is employed. By weighting phase and amplitude variation differently,
clustering results can be obtained from multiple perspectives. Routines
for function synchronization, as well as clustering procedures based on
the extracted phase and amplitude components, are provided.

## How can I get fkcentroids?

<!-- Version 0.4.0 of this package is available on [CRAN](https://cran.r-project.org/package=fkcentroids): -->
<!-- ``` r -->
<!-- install.packages("fkcentroids") -->
<!-- library(fkcentroids) -->
<!-- ``` -->

You can install the development version (version 0.0.0.9001) of
`fkcentroids` via:

``` r
devtools::install_github("seungwoo-stat/fkcentroids")
library(fkcentroids)
```

## How do I use it?

Using the Berkeley growth curves from the **fda** package, the functions
can be clustered using `fkmeans()` (or `fkmedians()`) after decomposing
them into their phase and amplitude components. Weights between phase
and amplitude variations can be controlled by the argument
`alpha_scale`, which is set to 1 by default.

``` r
library(fkcentroids)
data(growth, package = "fda")
Ytilde <- cbind(growth$hgtm, growth$hgtf)
x <- growth$age
(res_pa <- fkmeans(Ytilde, x, t = seq(0, 1, length.out = 100),
                   sync_map = "auc", sync_args = 1, alpha_scale = 1, 
                   k = 2, nstart = 10))
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
(res_raw <- fkmeans(Ytilde, x, t = seq(0, 1, length.out = 100),
                    sync_map = "none", k = 2, nstart = 10))
#> Functional k-means clustering with 2 clusters of sizes 51, 42 
#> 
#> Clustering vector: 
#>  boy01  boy02  boy03  boy04  boy05  boy06  boy07  boy08  boy09  boy10  boy11 
#>      2      1      1      2      1      1      1      2      1      2      2 
#>  boy12  boy13  boy14  boy15  boy16  boy17  boy18  boy19  boy20  boy21  boy22 
#>      2      2      2      2      2      2      2      1      2      1      1 
#>  boy23  boy24  boy25  boy26  boy27  boy28  boy29  boy30  boy31  boy32  boy33 
#>      1      1      1      2      2      1      2      2      2      2      1 
#>  boy34  boy35  boy36  boy37  boy38  boy39 girl01 girl02 girl03 girl04 girl05 
#>      1      2      2      2      2      1      1      1      2      2      1 
#> girl06 girl07 girl08 girl09 girl10 girl11 girl12 girl13 girl14 girl15 girl16 
#>      1      1      2      1      2      1      1      1      1      2      1 
#> girl17 girl18 girl19 girl20 girl21 girl22 girl23 girl24 girl25 girl26 girl27 
#>      1      2      1      2      2      1      1      1      2      1      1 
#> girl28 girl29 girl30 girl31 girl32 girl33 girl34 girl35 girl36 girl37 girl38 
#>      1      1      1      1      1      2      1      2      1      1      2 
#> girl39 girl40 girl41 girl42 girl43 girl44 girl45 girl46 girl47 girl48 girl49 
#>      2      2      1      1      2      1      1      1      1      1      2 
#> girl50 girl51 girl52 girl53 girl54 
#>      1      2      1      2      2 
#> 
#> Within cluster sum of squares by cluster: 
#> [1] 9265.661 9449.935
#>  (between_SS / total_SS = 47.7 %)
#> 
#> Available components: 
#>  [1] "cluster"        "centers.Ytilde" "totss"          "withinss"      
#>  [5] "tot.withinss"   "betweenss"      "size"           "iter"          
#>  [9] "ifault"         "alpha0"
```

It can be checked that the correct classification rate using phase and
amplitude components of each curve is 90/93, whereas functional
*k*-means clustering conducted directly on the observed curves yields
correct classification rate of 58/93. Refer to the package documentation
for details.

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
