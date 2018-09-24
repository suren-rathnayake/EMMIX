rdmvn <- function (n, mean, cov) {

  cov <- as.matrix(cov)
  p <- nrow(cov)

  if (nrow(cov) != ncol(cov)) {
    stop("cov must be a square matrix")
  }
  if (length(mean) != nrow(cov)) {
    stop("mean and cov have non-conforming size")
  }
  rmvnorm(n, mean = mean, sigma = cov, method = "chol")
}

rdmvt <- function (n, mean, cov, nu) {

  cov <- as.matrix(cov)
  p <- nrow(cov)
  u <- rgamma(n, nu/2, nu/2)
  t(t(rdmvn(n, mean = rep(0, p), cov = cov)/sqrt(u)) + mean)
}

rdmsn <- function (n, mean, cov, del) {

  x <-  rdmvn(n, mean, cov)
  z <-  abs(rnorm(n))
  as.matrix(z %*% t(del) + x)
}

rdmst <- function (n, mean, cov, nu, del) {

  p <- nrow(cov)
  u <- rgamma(n, nu/2, nu/2)
  x <- t(t(rdmvn(n, mean = rep(0, p), cov = cov) /sqrt(u)) + mean)
  z <- abs(rnorm(n) / sqrt(u))
  as.matrix(z %*% t(del) + x)
}

ddmvn <- function (dat, mean, cov) {

  p <- nrow(cov)
  exp(ddmix(dat, 1, "mvn", mean, cov, 0, rep(0, p)))
}

ddmvt <- function (dat, mean, cov, nu) {

  p <- nrow(cov)
  exp(ddmix(dat, 1, "mvt", mean, cov, nu, rep(0, p)))
}

ddmsn <- function (dat, mean, cov, del) {

  exp(ddmix(dat, 1, "msn", mean, cov, 0, del))
}

ddmst <-function (dat, mean, cov, nu, del) {

  exp(ddmix(dat, 1, "mst", mean, cov, nu, del))
}
