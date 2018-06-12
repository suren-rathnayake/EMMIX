EMMIX2 <- function (dat, g, init, nvcov = 0, neq = 0, itmax = 1000, epsilon = 1e-06,
    debug = TRUE)
{
    dat <- as.matrix(dat)
    obj <- emmixfit3(dat, g, init, nvcov, neq, itmax, epsilon)
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
        ret <- obj
        if (debug) {
            msg <- paste("Restricted ", g, "- Component Multivariate Normal Mixture Model")
            cat("\n-----------------------\n\n")
            cat(msg, "\n")
            cat("\n-----------------------\n\n")
            print(obj[1:8])
        }
    }
    ret
}
