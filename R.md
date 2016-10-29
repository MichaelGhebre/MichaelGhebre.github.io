---
layout: page
title: R
permalink: /R/
---

We wrote R function to implement the model-based clustering from scratch. The EM-algorithm is initialized using several heuristic-clustering techniques, such as kmeans, kmedoids and hierarchical, here is the code. The singularity issue was dealt by adding small value to the diagonal of the variance-covariance matrix if the variance covariance is not positive definite, here is the code.

```r
sigma <- function(x, y){
z <- 2*x + 3*y + x*y
return(z)
}
```
