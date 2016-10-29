---
layout: page
title: R
permalink: /R/
---

We wrote R function to implement the model-based clustering from scratch. The EM-algorithm is initialized using several heuristic-clustering techniques, such as kmeans, kmedoids and hierarchical, here is the code. The singularity issue was dealt by adding small value to the diagonal of the variance-covariance matrix if the variance covariance is not positive definite, here is the code.

```r
#'  gmmEM 
#'
#' This is model-based clustering technique using EM algorithm which initialized using k-means or k-medoids algorithms
#' @param x a numeric  matrix or dataframe  for clusterng; c is  number of clusters;  
#' @return  Estiamted means and variance-covariance matrix in each cluster, and BIC, loglikelihood, number of iteration, clusters (class) and number of parameters
#' @export
#' @examples
#' gmmEM(iris[,1:4], c=2,initialize = c('kmeans'))
#' 


gmmEM <- function(x, c = 1:10, initialize = c("kmeans", "kmedoids", "fuzzykmeans","hierarchical"),           iter = NULL, tol = NULL, ...) {
    
    require(clusterSim)  # for initialization the kmeans clustering algorithm
    require(FactMixtAnalysis)  #  for calculating confusion matrix
    require(mvtnorm)  # for multivariate density function
    require(e1071)
    

    x <- as.matrix(x)
    
    if (is.null(x)) 
        return(NULL)
    
    # X is non-missing
    if (any(is.na(x))) {
        warning("NA's in the dataframe")
        return(NULL)
    }
    
    
    if (nrow(x) == 1) {
        warning(" The input should be more than a single observation")
        return(NULL)
    }
    
    
    if (any(!is.numeric(x))) {
        warning("The input is not a numeric dataframe or matrix")
        return(NULL)
    }
    
    
    
    n <- nrow(x)
    d <- ncol(x)
    
    
    #  hierarchical clustering 
    
    hhclust <- function(x, method ="complete", cluster = 1:20){
      
      tree <- hclust(dist(x), method = "complete")
      group <- cutree(tree,c)
      
      return(list(Call=tree, cluster=group))
    }
    
    
    
    if (c == 1 && d == 1) {
        j <- 1  # number of clusters
        
        # initialise the algorithm using k-means or k-medoids
        
        if (initialize == "kmedoids") {
            k <- pam(x, c, do.swap = F)
       } else if (initialize == "fuzzykmeans") {
          k <- cmeans(x, c, iter.max = 100,  dist = "manhattan", method="cmeans")
          
        } else if (initialize == "hierarchical") {
          
          k <- hhclust(x, cluster=c)
          
        } else {
            k <- kmeans(x, c, nstart = 1)
        }
        
        mu <- mean(x[k$cluster == 1])
        sigma <- sd(x[k$cluster == 1])
        
        p <- sum(k$cluster == 1)/length(x)
        
        loglik <- sum(log(p * (dnorm(x, mu, sigma))))
        
        p <- (2 * j - 2 + 3 * j * d + j * d^2)/2  # parameters or degree of freedom ; j = number of clusters
        bic <- -2 * loglik + p * log(n)
        
        return(list(mean = mu, sigma = sigma, loglikhood = loglik, n = n, 
            BIC = bic, df = p, clusters = j))
        
    } else if (c >= 2 & d == 1) {
        
        # initialize the EM-algorithm
        if (initialize == "kmedoids") {
            k <- pam(x, c, do.swap = F)
        } else if (initialize == "fuzzykmeans") {
          k <- cmeans(x, c, iter.max = 100,  dist = "manhattan", method="cmeans")
          
        } else if (initialize == "hierarchical") {
          k <- hhclust(x, cluster=c)
          
        } else {
            k <- kmeans(x, c, nstart = 25)
        }
        
        j <- c
        mu <- lapply(1:c, function(i) mean(x[k$cluster == i, ]))
        mu <- lapply(mu, function(x) {
            replace(x, is.na(x), .Machine$double.eps)
        })
        
        sigma <- lapply(1:c, function(i) sd(x[k$cluster == i, ]))
        sigma <- lapply(sigma, function(x) {
            replace(x, is.na(x), .Machine$double.eps)
        })
        
        
        # sigma <- lapply(sigma,function(x)
        # replace(sigma,x==0,.Machine$double.eps)) # replacing zeros sigma
        
        p <- lapply(1:c, function(i) sum(k$cluster == i)/nrow(x))
        p <- lapply(p, function(x) {
            replace(x, is.na(x), .Machine$double.eps)
        })
        
        
        Ntau <- lapply(1:c, function(i) p[[i]] * (dnorm(x, mu[[i]], sigma[[i]])))
        Ntau <- sapply(Ntau, cbind)
        SUMtau <- apply(Ntau, 1, sum)
        loglik <- sum(log(SUMtau))
        
        tol = 1e-06
        iter <- 0
        tau <- 0
        
        
        # for loop
        for (i in 1:1000) {
            
            iter <- iter + 1
            
            Ntau <- lapply(1:c, function(i) p[[i]] * (dnorm(x, mu[[i]], 
                sigma[[i]])))
            Ntau <- sapply(Ntau, cbind)
            SUMtau <- apply(Ntau, 1, sum)
            
            tau <- Ntau/SUMtau
            p <- lapply(1:c, function(i) sum(tau[, i])/nrow(x))
            p <- lapply(p, function(x) {
                replace(x, is.na(x), .Machine$double.eps)
            })
            
            mu <- lapply(1:c, function(i) sum(tau[, i] * x)/sum(tau[, i]))
            mu <- lapply(mu, function(x) {
                replace(x, is.na(x), .Machine$double.eps)
            })
            
            sigma <- lapply(1:c, function(i) sqrt(sum(tau[, i] * (x - mu[[i]])^2)/sum(tau[, 
                i])))
            sigma <- lapply(sigma, function(x) {
                replace(x, is.na(x), .Machine$double.eps)
            })
            
            loglik_0 <- loglik
            
            Ntau <- lapply(1:c, function(i) p[[i]] * (dnorm(x, mu[[i]], 
                sigma[[i]])))
            Ntau <- sapply(Ntau, cbind)
            SUMtau <- apply(Ntau, 1, sum)
            loglik <- sum(log(SUMtau))
            
            del_loglik <- loglik - loglik_0
            
            if ((abs(del_loglik) < tol) || (is.nan(abs(del_loglik)))) {
                
                break
                
            }
        }
        
        
        prob <- as.data.frame(tau)
        class <- apply(prob, 1, which.max)
        
        par <- (2 * j - 2 + 3 * j * d + j * d^2)/2  # parameters or degree of freedom; j is number of clusters
        
        bic <- -2 * loglik + par * log(n)
        
        return(list(lambda = p, mu = mu, sigma = sigma, loglikhood = loglik, 
            n = n, BIC = bic, df = par, number_of_iteration = iter, clusters = j, 
            class = class))
        
    } else if (c >= 2 & d >= 2) {
      
        # initialise the algorithm using k-means or kmedoids
        if (initialize == "kmedoids") {
          k <- pam(x, c, do.swap = F)
          
        } else if (initialize == "fuzzykmeans") {
          k <- cmeans(x, c, iter.max = 100,  dist = "manhattan", method="cmeans")
          
        } else if (initialize == "hierarchical") {
          k <- hhclust(x, cluster=c)
          
        } else {
          k <- kmeans(x, x[initial.Centers(x, c), ])
        }
        
         j <- length(unique(k$cluster))
        
    
  
        mu <- lapply(1:c, function(i) apply(x[k$cluster == i, ], 2, mean))
        sigma <- lapply(1:c, function(i) cov(x[k$cluster == i, ]))
        sigma <- lapply(sigma, function(x) sigmaFixer(x))  # fixing singularity 
        
        p <- lapply(1:c, function(i) sum(k$cluster == i)/nrow(x))
        
        tol <- 1e-06
        iter <- 0
        tau <- 0
        
        Ntau <- lapply(1:c, function(i) p[[i]] * (dmvnorm(x, mu[[i]], sigma[[i]])))
        Ntau <- sapply(Ntau, cbind)
        SUMtau <- apply(Ntau, 1, sum)
        
        loglik <- sum(log(SUMtau))
        
        # for loop
        for (i in 1:1000) {
            
            iter <- iter + 1
            
            Ntau <- lapply(1:c, function(i) p[[i]] * (dmvnorm(x, mu[[i]], 
                sigma[[i]])))
            Ntau <- sapply(Ntau, cbind)
            SUMtau <- apply(Ntau, 1, sum)
            
            tau <- Ntau/SUMtau
            p <- lapply(1:c, function(i) sum(tau[, i])/nrow(x))
            
            mu <- lapply(1:c, function(i) apply(x, 2, function(x) sum(tau[, 
                i] * x)/sum(tau[, i])))
            
            x_m <- lapply(1:c, function(i) sweep(x, 2, mu[[i]], "-"))
            
            sigma <- lapply(1:c, function(i) (t(as.matrix(x_m[[i]])) %*% 
                (as.matrix(x_m[[i]]) * tau[, i]))/sum(tau[, i]))
            sigma <- lapply(sigma, function(x) sigmaFixer(x))  # fixing singularity 
            
            loglik_0 <- loglik
            Ntau <- lapply(1:c, function(i) p[[i]] * (dmvnorm(x, mu[[i]], 
                sigma[[i]])))
            Ntau <- sapply(Ntau, cbind)
            SUMtau <- apply(Ntau, 1, sum)
            loglik <- sum(log(SUMtau))
            
            del_loglik <- loglik - loglik_0
            
            
            if ((abs(del_loglik) < tol) || (is.nan(abs(del_loglik)))) {
                
                break
                
            }
            
        }
        
        prob <- as.data.frame(tau)
        class <- apply(prob, 1, which.max)
        par <- (2 * j - 2 + 3 * j * d + j * d^2)/2  # parameters or degree of freedom
        
        bic <- -2 * loglik + par * log(n)
        
        return(list(lambda = p, mu = mu, sigma = sigma, loglikhood = loglik, 
            number_of_iteration = iter, n = n, BIC = bic, df = par, class = class))
        
    }
} 

```

#Singularity Issues in Model-based Clustering 

If the variance-covariance is not positive definite, a small value is added to the diagonal of the matrix till it is positive definite using the following R function. 

```r
sigmaFixer <- function(sigma, ...) {
    sigma <- as.matrix(sigma)
    r <- dim(sigma)[1]
    fixrate <- 0.01
    covfixmat <- array(1, c(r, r)) + fixrate * diag(1, r)
    min_limit <- .Machine$double.eps * 10
    
    if (!all(is.finite(sigma))) {
        warning("covariance matrix is not finite")
    }
    
    # Enforces squareness and symmetricity
    nsigma <- sigma - ((sigma - t(sigma))/2)
    iter <- 0
    
    # Checking the covariance matrix is not positive definite
    while (postDef(nsigma) == 0 & (iter < 10000)) {
        
        iter <- iter + 1
        d <- diag(nsigma)
        if (any(d <= min_limit)) {
            m <- max(abs(d)) * fixrate
            neg <- min(d)
            
            if (neg < 0) {
                addit <- (m - neg) * diag(1, r)
            } else {
                if (m < min_limit) {
                  m <- min_limit
                }
                addit <- m * diag(1, r)
            }
            nsigma <- nsigma + addit
        } else {
            # Increase the diagonal values by 1%
            nsigma <- nsigma * covfixmat
        }
    }
    # return(list(nsigma, iter))
    return(nsigma)
}


# testing for postive definite
postDef <- function(Sigma) {
    Sigma <- as.matrix(Sigma)
    q <- try(chol(Sigma)[2], TRUE)
    if (q == 0) {
        t <- 1
    } else {
        t <- 0
    }
    return(t)
} 
```
