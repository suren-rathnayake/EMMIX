EMMIX <- function (dat, g, distr = "mvn", ncov = 3, clust = NULL, init = NULL,
    itmax = 1000, epsilon = 1e-06, nkmeans = 0, nrandom = 10,
    nhclust = FALSE, debug = FALSE, initloop = 20)
{
    dat <- as.matrix(dat)

    if (!is.null(init) | !missing(init)) {

        obj <- emmixfit2(dat, g, init, distr, ncov, itmax, epsilon)
    } else {

        if (is.null(clust) | missing(clust)) {

            init <- try(init.mix(dat, g, distr, ncov, nkmeans, nrandom, 
                        nhclust, initloop))

            if (!is.null(init)) {

                obj <- emmixfit2(dat, g, init, distr, ncov, itmax, epsilon)
            } else {

                warning("not find initial values")
                obj <- list()
                obj$error = 20
            }
        } else {

            obj <- emmixfit1(dat, g, clust, distr, ncov, itmax,
                epsilon, initloop)
        }
    }

    error <- obj$error
    ret <- NULL
    msg <- switch(tolower(error), `1` = paste("stopped at (did not converge within) ",
        itmax, " iterations"), `2` = paste("density fails at initial steps! "),
        `3` = paste("allocation fails at initial steps"), `12` = paste("density fails at estps! "),
        `13` = paste("allocation fails at esteps"), `20` = paste("not find initials"))
    if (error > 1) {
        cat("\n-----------------------\n")
        warning("error code=", error, "\n", msg, "\n")
        cat("\n-----------------------\n\n")
    }
    if (error <= 1) {
        ICL <- getICL(dat, nrow(dat), ncol(dat), g, distr, ncov,
            obj$pro, obj$mu, obj$sigma, obj$dof, obj$delta, obj$clust)
        ret <- obj
        ret$ICL <- ICL$ICL
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
    class(ret) <- c(class(ret), "emmix")
    ret
}
