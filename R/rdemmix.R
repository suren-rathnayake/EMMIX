rdemmix <- function (nvect, distr, mu, sigma, dof = NULL, delta = NULL)
{

  g <- dim(sigma)[3]
  if (g != dim(mu)[2])
    stop("Number of components in mu and sigma are not the same.")

  p <- dim(sigma)[1]
  if (p != dim(mu)[1])
    stop("Number of dimensions in mu and sigma are not the same.")


  if ((distr == "mvt") || (distr == "mst")) {
    if (is.null(dof)) {
      stop("The dof needs to specified for mvt and mst mixture models.")
    }
  }

  if ((distr == "mst") || (distr == "msn")) {
    if (is.null(delta)) {
      stop("The delta needs to specified for mvt and mst mixture models.")
    }
  }

  if (length(c(nvect)) != g)
      stop("nvect should be a vector")
  ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3,
      mst = 4, 5)
  if (ndist > 4)
      stop("the model specified is not available yet")

    if (is.null(dof))
       dof <- rep(4, g)
    if (is.null(delta))
       delta <- array(0, c(p, g))

  if (length(c(mu)) != (p * g))
    stop(paste("mu should be a ", p, "by", g, "matrix!"))
  if (length(c(sigma)) != (p * p * g))
    stop(paste("sigma should be a ", p, "by", p, "by", g,
            " array!"))
  #if (length(c(dof)) != g)
#    stop(paste("dof should be a ", g, " vector!"))
  #if (length(c(delta)) != (p * g))
#    stop(paste("delta should be a ", p, "by", g, " array!"))
  mu <- array(mu, c(p, g))
  sigma <- array(sigma, c(p, p, g))
  delta <- array(delta, c(p, g))
  dat <- array(0, c(10, p))
    mvrand <- function(n, p, ndist, mean, cov, nu, del) {
        switch(ndist, `1` = rdmvn(n, mean = mean, cov = cov),
            `2` = rdmvt(n, mean = mean, cov = cov, nu = nu),
            `3` = rdmsn(n, mean = mean, cov = cov, delta = del),
            `4` = rdmst(n, mean = mean, cov = cov, nu = nu, delta = del))
    }
    if (g >= 1)
        for (h in 1:g) {
            if (nvect[h] > 0)
                dat <- rbind(dat, mvrand(nvect[h], p, ndist,
                  mu[, h], sigma[, , h], dof[h], delta[, h, drop = F]))
        }
    dat[-(1:10), ]
}


rdemmix2 <- function (n, distr, pro, mu, sigma, dof = NULL, delta = NULL)
{
  g <- dim(sigma)[3]
  if ((g != dim(mu)[2]) || (g != length(pro)))
    stop("Number of components in mu and sigma are not the same.")

  p <- dim(sigma)[1]
  if (p != dim(mu)[1])
    stop("Number of dimensions in mu and sigma are not the same.")

  n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
  nn <- n0
    if (length(nn) < g) {
        nn <- rep(0, g)
        for (i in as.numeric(names(n0))) nn[i] <- n0[paste(i)]
    }
    names(nn) <- NULL
    rdemmix(nn, distr, mu, sigma, dof, delta)
}



rdemmix3 <- function (n, distr, pro, mu, sigma, dof = NULL, delta = NULL)
{
    g <- dim(sigma)[3]
    if (g != dim(mu)[2])
      stop("Number of components in mu and sigma are not the same.")

    p <- dim(sigma)[1]
    if (p != dim(mu)[1])
      stop("Number of dimensions in mu and sigma are not the same.")

    if (length(pro) != g)
        stop(paste("pro should be a ", g, " vector!"))
    n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
    nn <- n0
    if (length(nn) < g) {
        nn <- rep(0, g)
        for (i in as.numeric(names(n0))) nn[i] <- n0[paste(i)]
    }
    names(nn) <- NULL
    dat <- rdemmix(nn, distr, mu, sigma, dof, delta)
    list(data = dat, cluster = rep(1:g, nn))
}
