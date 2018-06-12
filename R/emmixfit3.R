emmixfit3 <- function (dat, g, init, nvcov = 0, neq = 0, itmax, epsilon)
{
    if ((neq != 0) && length(nvcov) != neq)
        stop("the length of nvcov (specified equal cov components ids) should be same as neq!")
    dat <- as.matrix(dat)
    n <- nrow(dat)
    p <- ncol(dat)
    if (is.null(init) | missing(init))
        stop("init should be provided")
    pro <- init$pro
    mu <- init$mu
    sigma <- init$sigma
    obj <- .C("emmixfit3", PACKAGE = "EMMIX", as.double(dat),
        as.integer(n), as.integer(p), as.integer(g), as.integer(nvcov),
        as.integer(neq), pro = as.double(pro), mu = as.double(mu),
        sigma = as.double(sigma), tau = double(n * g), sumtau = double(g),
        loglik = double(1), lk = double(itmax), aic = double(1),
        bic = double(1), clust = integer(n), error = integer(1),
        as.integer(itmax), as.double(epsilon))[c(7:18)]
    lk <- obj$lk
    lk <- lk[lk != 0]
    list(distr = "mvn", error = obj$error, loglik = obj$loglik,
        bic = obj$bic, aic = obj$aic, pro = obj$pro, mu = array(obj$mu,
            c(p, g)), sigma = array(obj$sigma, c(p, p, g)), clust = obj$clust,
        tau = array(obj$tau, c(n, g)), lk = lk)
}
