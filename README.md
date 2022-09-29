
<!-- README.md is generated from README.Rmd. Please edit that file -->

# revpathsensitivity

<!-- badges: start -->
<!-- badges: end -->

This code implements reverse sensitivity testing on compound Poisson
processes. This is joint work with Silvana M. Pesenti and Sebastian
Jaimungal.

## Installation

You can install the development version of revpathsensitivity from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("emmakroell/revpathsensitivity")
```

## Example

``` r
library(revpathsensitivity)
```

In this example, we set the jump size distribution to
![\Gamma(\alpha, \beta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CGamma%28%5Calpha%2C%20%5Cbeta%29 "\Gamma(\alpha, \beta)")
where
![\alpha = 2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha%20%3D%202 "\alpha = 2")
and
![\beta = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta%20%3D%201 "\beta = 1").
We first create a class containing all the information needed on this
random variable.

``` r
# jump size is Gamma(2,1) distributed
# Set up distribution class
gamma_2_1 <- new_RPS_dist(dist_fun = function(x,parms) pgamma(x,shape=parms$alpha,
                                                              rate=parms$beta),
                          dens_fun = function(x,parms) dgamma(x,shape=parms$alpha,
                                                              rate=parms$beta),
                          sim_fun = function(x,parms) rgamma(x,shape=parms$alpha,
                                                             rate=parms$beta),
                          char_fun = function(x,parms) (1 - 1i*x/parms$beta)^(-parms$alpha),
                          min = 0,
                          max = qgamma(1-1e-10,shape=2,rate=1),
                          mean_fun = function(parms) parms$alpha/parms$beta,
                          parms = list(alpha = 2, beta=1))
```

We decide to stress the 90% VaR of the model up 10%. This code imposes
the stress and simulates under the model.

``` r
ex_gamma_VaR <- stressed_sim(kappa = 5, jump_dist = gamma_2_1, 
                             stress_type = "VaR",
                             stress_parms = list(c=0.9,VaR_stress=1.1),
                             Npaths=1e4, end_time=1, dt=2e-3)
#> No stress time specified. Applying stress at terminal time. 
#> eta: 0.4833721
```

Here we plot the output of the stressed model. We compare the paths
under the stressed model to those under the baseline model, as well as
to how the intensity of the jumps have changed under the stressed
probability measure.

``` r
plot_paths(ex_gamma_VaR, Npaths=15, quantiles=list(lower=0.1,upper=0.9))
#> $X
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

    #> 
    #> $kappa

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

### CVaR constraint

Now we want to impose an additional constraint, corresponding to an 8%
increase in CVaR (also known as Expected Shortfall).

``` r
ex_gamma_CVaR <- stressed_sim(kappa = 5, jump_dist = gamma_2_1,
                              stress_type = "CVaR", 
                              stress_parms = list(c=0.9,
                                                 VaR_stress=1.1,
                                                 CVaR_stress=1.08),
                              Npaths=1e4, end_time=1, dt=2e-3)
#> No stress time specified. Applying stress at terminal time. 
#> eta: 0.3326332 -0.04110227
```

Here we plot the paths of the stressed model as well as the jump
intensity over time.

``` r
plot_paths(ex_gamma_CVaR, Npaths=15, quantiles=list(lower=0.1,upper=0.9))
#> $X
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

    #> 
    #> $kappa

<img src="man/figures/README-unnamed-chunk-6-2.png" width="100%" />

### Bivariate example

We can also examine the case where
![X_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_t "X_t")
is bivariate and a stress is applied to only the first component,
![X^1_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X%5E1_t "X^1_t").
In this example, we set the severity distributions to both be
![\Gamma(2,1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CGamma%282%2C1%29 "\Gamma(2,1)")
with a dependence relationship given by the
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")-copula
with correlation 0.8 and 3 degrees of freedom. The implementation of a
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")-copula
with fixed marginals is done using the `copula` package.

``` r
gamma_t <- new_RPS_dist_biv(copula=copula::tCopula(0.8,df=3,dim=2),
                            margins=c("gamma", "gamma"),
                            sim_fun = function(x,parms) rgamma(x,shape=parms$alpha1,
                                                           rate=parms$beta1),
                            dens_fun = function(x,parms) dgamma(x,shape=parms$alpha1,
                                                            rate=parms$beta1),
                            char_fun = function(x,parms) (1 - 1i*x/parms$beta1)^
                              (-parms$alpha1),
                            mean_fun = function(parms) parms$alpha1/parms$beta1,
                            parms = list(alpha1 = 2, beta1=1, alpha2=2, beta2=1))
```

Here we examine the effect of a 10% increase in
![\text{VaR}\_{0.9}(X^1_T)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctext%7BVaR%7D_%7B0.9%7D%28X%5E1_T%29 "\text{VaR}_{0.9}(X^1_T)")
at the 90% level.

``` r
gamma_t_ex <- stressed_sim_biv(kappa = 5, jump_dist = gamma_t,
                               stress_type = "VaR",
                               stress_parms = list(c=0.9,VaR_stress=1.1),
                               Npaths=1e4, end_time=1, dt=2e-3)
#> No stress time specified. Applying stress at terminal time. 
#> eta: 0.4833721
```

``` r
plot_paths_biv(gamma_t_ex, Npaths=15, quantiles=list(lower=0.1,upper=0.9))
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

We can also plot the contours of the copula density under
![\mathbb{P}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BP%7D "\mathbb{P}")
and
![\mathbb{Q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BQ%7D "\mathbb{Q}").
Here we plot it at the terminal time,
![t=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t%3D1 "t=1"):

``` r
plot_copula_dens_contour(gamma_t_ex,time=1)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />
