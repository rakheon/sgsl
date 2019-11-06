# d2wlasso: data-driven weighted lasso
The R package `sgsl` implements the sparse group-subgroup lasso.

The reference for the original pliable lasso can be found at:
* [Identification of important regressor groups, subgroups
and individuals via regularization methods: application to
gut microbiome data](https://doi.org/10.1093/bioinformatics/btt608) by Tanya P. Garcia et al (2014).

## Installation

```
devtools::install_github("rakheon/sgsl", force = TRUE)
```

## Example

```
# data generation for linear models
x = matrix(rnorm(100*5, 0, 1),100,5)
z = matrix(rbinom(100, 1, 0.5),100,1)
y = matrix(z[,1] + 2*x[,1] - 2*x[,2] + rnorm(100, 0, 1), 100)

# variable selection with d2wlasso for linear models



```
