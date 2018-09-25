dmix <- function (dat, model, ...) {

  n <- nrow(dat)
  p <- ncol(dat)

  mu <- model$mu 
  sigma <- model$sigma 
  pro <- model$pro 

  g <- model$g
  if (!exists("ncov", model)) {
    stop("The model parameter must contain `ncov` which describe \
      structure of mixture model. See help(EMMIX).")
  }
  ncov <- model$ncov

  if (!exists("distr", model)) {
    stop("The model parameter must contain `distr` which describe \
      the type of component distributions. See help(EMMIX).")
  }  
  distr <- model$distr

  if (exists("dof", model)) {

    dof <- model$dof
  } else {

    dof <- rep(1000, g)
  }

  if (exists("delta", model)) {

    delta <- model$delta
  } else {

    delta <- array(0, c(p, g))
  }

  ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3,
                   mst = 4, 5)
  if (ndist > 4)
    stop("the model specified is not available yet")
    
  dat <- as.matrix(dat)
    
  if (n == 1 & (ncol(dat) == 1))
    dat <- t(dat)
    
  if (nrow(dat) != n | ncol(dat) != p)
    stop("dat does not match n and p.")
    
  if (length(c(mu)) != (p * g))
    stop(paste("mu should be a ", p, "by", g, "matrix!"))
    
  if (length(c(sigma)) != (p * p * g))
    stop(paste("sigma should be a ", p, "by", p, "by", g, " array!"))
    
  if (length(c(dof)) != g)
    stop(paste("dof should be a ", g, " vector!"))

  if (length(c(delta)) != (p * g))
    stop(paste("delta should be a ", p, "by", g, " array!"))
    
  obj <- .C("ddmix", PACKAGE = "EMMIX", as.double(dat), as.integer(n),
        as.integer(p), as.integer(g), as.integer(ndist), as.double(mu),
        as.double(sigma), as.double(dof), as.double(delta), den = double(n *
            g), error = integer(1))[10:11]
    
  if (obj$error)
    stop("error")
    
  den <- matrix(obj$den, ncol = g)
  den 
}
