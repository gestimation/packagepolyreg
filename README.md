
<!-- README.md is generated from README.Rmd. Please edit that file -->

# packagepolyreg

<!-- badges: start -->
<!-- badges: end -->

The goal of packagepolyreg is to test polyreg function.

## Installation

You can install the development version of packagepolyreg from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("gestimation/packagepolyreg")
```

## Example

This is a basic example of direct polynomial regression analysis.
Dataset bmt is a well-known competing risks dataset that is included in
package mets.

``` r
library(packagepolyreg)
data(bmt)
result <- polyreg(nuisance.model = Event(time, cause)~age+tcell, exposure = 'platelet',
cens.model = Event(time,cause==0)~+1, data = bmt, 
effect.measure1='RR', effect.measure2='RR', time.point=24, outcome.type='COMPETINGRISK')
```
