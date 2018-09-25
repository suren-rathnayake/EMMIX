rmix <- function (n, model, ...) {

  g <- model$g
  p <- dim(model$sigma)[1]

  mu <- model$mu 
  sigma <- model$sigma 
  pro <- model$pro 

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

  if (g != dim(mu)[2])
    stop("Number of components in mu is not model$g.")

  if (g != dim(sigma)[3])
    stop("Number of components in sigma is not model$g.")

  if (p != dim(mu)[1])
    stop("Number of dimensions in mu and sigma are not the same.")

  if (length(pro) != g)
        stop(paste("pro should be a ", g, " vector!"))

  n0 <- table(sample(1 : g, n, replace = TRUE, prob = pro))
  nn <- n0
  if (length(nn) < g) {
    nn <- rep(0, g)
    for (i in as.numeric(names(n0))) nn[i] <- n0[paste(i)]
  }
  
  names(nn) <- NULL
  dat <- rdemmix(nn, distr, mu, sigma, dof, delta)

  if (labels == TRUE) {
    return(list(data = dat, labels = rep(1 : g, nn)))
  } else {
    return(dat)
  }
}