bootstrap <- function (x, model, B = 99, replace = TRUE,
    itmax = 1000, epsilon = 1e-05, ...) {

  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  g <- model$g
  if (!exists("ncov", model)) {
    error("The model parameter must contain `ncov` which describe \
      structure of mixture model. See help(EMMIX).")
  }
  ncov <- model$ncov

  if (!exists("distr", model)) {
    error("The model parameter must contain `distr` which describe \
      the type of component distributions. See help(EMMIX).")
  }  
  distr <- model$distr

  if (missing(model))
     stop("Please provide `model` parameters.")
  
  counter <- 0
  if (distr == "mvn") {

    ret <- matrix(NA, nrow = B, ncol = g * (1 + p + p * p))
  }  else if (distr == "mvt") {

    ret <- matrix(NA, nrow = B, ncol = g * (2 + p + p * p))
  } else {

    ret <- matrix(NA, nrow = B, ncol = g * (2 + 2*p + p * p))
  }

  for (i in 1 : (2 * B)) {

    if (replace) {

       dat <- x[sample(1:n, n, replace = TRUE), ]
    } else {

      dat <- rdemmix3(n, distr, model$pro, model$mu,
                      model$sigma, model$dof, model$delta)
    }
    
    obj <- emmixfit2(dat, g, model, distr, ncov, itmax, epsilon)
    
    if (obj$error > 1)
      next
    
    counter <- counter + 1

    if (distr == "mvn") {    

        ret[counter, ] <- c(obj$pro, obj$mu, obj$sigma)
    }

    if (distr == "mvt") {    

        ret[counter, ] <- c(obj$pro, obj$mu, obj$sigma, obj$dof)
    }

    if (distr == "msn") {

        obj$dof <- rep(NA, g)    
        ret[counter, ] <- c(obj$pro, obj$mu, obj$sigma, obj$dof, obj$delta)
    }

    if (distr == "mst")  {
      
        ret[counter, ] <- c(obj$pro, obj$mu, obj$sigma, obj$dof, obj$delta)
    }

    
    if (counter >= B)
      break
  }

  std <- sqrt(apply(ret[1 : counter, ], MARGIN = 2, FUN = "var"))
  
  se_pi <- std[1 : g]

  se_mu <- matrix(std[g + 1 : (p * g)], nrow = p, ncol = g)

  se_sigma <- array(NA, c(p, p, g))
  for (i in 1 : g) 
    se_sigma[,, i] <- matrix(std[(g + p * g) + (p*p*(i-1)) + 1 : (p * p)], 
                              ncol = p)
  se_err <- list(se_pi = se_pi, se_mu = se_mu, se_sigma = se_sigma)
  
  if ((distr == "mvt") || (distr == "mst")) 
    se_err$se_dof <- std[((g + p * g) + (p * p * g)) + 1 : g] 

  if ((distr == "msn") || (distr == "mst"))
    se_err$se_delta <- matrix(std[((g + p * g) + (p * p * g) + g) + 
                                1 : (g * p)], nrow = p, ncol = g) 

  se_err 
}
