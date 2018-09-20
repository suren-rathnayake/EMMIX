EMMIX <- function (dat, g, distr = "mvn", ncov = 3, clust = NULL, 
            init = NULL, itmax = 1000, epsilon = 1e-06, nkmeans = 0, 
            nrandom = 10, nhclust = FALSE, debug = FALSE, initloop = 20) {

  dat <- as.matrix(dat)
  obj <- list(loglik = -Inf)

  if (!is.null(init) | !missing(init)) {
  
    if (((distr == "mvn") || (distr == "msn")) && (!exists("init$dof"))) 
      init$dof <- rep(100, g)

    if (((distr == "msn") || (distr == "mst")) && (!exists("init$delta"))) 
      init$delta <- rep(100, g) 

    obj <- emmixfit2(dat, g, init, distr, ncov, itmax, epsilon)
  } 

  if (!is.null(clust)) {
    obj_cl <- emmixfit1(dat, g, clust, distr, ncov, itmax, epsilon, initloop) 

    if (obj_cl$error <= 1) {
       if (obj_cl$loglik > obj$loglik)
        obj <- obj_cl
    }
  }

  if ((nkmeans >= 1) || (nrandom >= 1) || (nhclust >= 1)) {

    init <- init.mix(dat, g, distr, ncov, nkmeans, nrandom, nhclust, initloop)

    if (!is.null(init)) {

      obj_cl <- emmixfit2(dat, g, init, distr, ncov, itmax, epsilon)

      if (obj_cl$error <= 1) {

        if (obj_cl$loglik > obj$loglik)
          obj <- obj_cl
      }
    } else {

      warning("not find initial values")
      obj <- list()
      obj$error = 20
    }

  }
  
  error <- obj$error
  ret <- NULL
  msg <- switch(tolower(error), 
    `1` = paste("stopped at (did not converge within) ",
        itmax, " iterations"), 
    `2` = paste("density fails at initial steps! "),
    `3` = paste("allocation fails at initial steps"), 
    `12` = paste("density fails at estps! "),
    `13` = paste("allocation fails at esteps"), 
    `20` = paste("not find initials")
   )
   
   if (error > 1) {
        cat("\n-----------------------\n")
        warning("error code=", error, "\n", msg, "\n")
        cat("\n-----------------------\n\n")
   }
   
   if (error <= 1) {
        ICL <- getICL(dat, nrow(dat), ncol(dat), g, distr, ncov,
            obj$pro, obj$mu, obj$sigma, obj$dof, obj$delta, obj$clust)
        ret <- obj
        ret$icl <- ICL$ICL
        if (debug) {
            msg <- switch(tolower(distr), mvn = paste(g, "- Component Multivariate Normal Mixture Model"),
                mvt = paste(g, "- Component Multivariate t-Mixture Model"),
                msn = paste(g, "- Component Multivariate Skew Normal Mixture Model"),
                mst = paste(g, "- Component Multivariate Skew-t Mixture Model"))
            cat("\n-----------------------\n\n")
            cat(msg, "\n")
            cat("\n-----------------------\n\n")
            switch(tolower(distr), mvn = print(obj[1:8]), mvt = print(obj[1:9]),
                msn = print(obj[c(1:8, 10)]), mst = print(obj[1:10]))
            print(ICL)
            cat("\n-----------------------\n")
        }
    }
    ret$g <- g
    ret$distr <- distr
    ret$ncov <- ncov
    
    if (distr == "mvn") {

        ret$dof <- ret$delta <- NULL
    }

    if (distr == "mvt") {
        
        ret$delta <- NULL
    }

    if (distr == "msn") {
        
        ret$dof <- NULL
    }
 
    ret$call <- match.call()
    class(ret) <- c(class(ret), "emmix", distr)
    ret
}
