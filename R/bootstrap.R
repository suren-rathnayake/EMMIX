bootstrap <- function (x, model, B = 99, replace = TRUE,
    itmax = 1000, epsilon = 1e-05, ...) {

  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  g <- model$g
  ncov <- model$ncov
  distr <- model$distr

  if (missing(model))
     stop("please provide the function model parameters.")
  
  counter <- 0
  nnn <- g * (1 + p + p * p + 1 + p)
  ret <- array(0, c(B, nnn))

  dimnames(ret) <- list(1 : B, c(paste("pi", 1 : g, sep = ""),
                    paste("mu", rep(1 : p, g), rep(paste(1 : g, sep = ""), 
                    rep(p, g)), sep = ""), 
                    paste("sigma", rep(paste(rep(1 : p, rep(p, p)), 
                    rep(1 : p, p), sep = ""), g), 
                    rep(paste(",", 1 : g, sep = ""), rep(p * p, g)), sep = ""), 
                    paste("dof",   1 : g, sep = ""), 
                    paste("delta", rep(1 : p, g), rep(paste(1 : g, sep = ""), 
                    rep(p, g)), sep = "")))

  for (i in 1 : (2 * B)) {

    if (replace)
       dat <- x[sample(1:n, n, replace = TRUE), ]
    else 
      dat <- rdemmix3(n, distr, model$pro, model$mu,
                      model$sigma, model$dof, model$delta)
    
    obj <- emmixfit2(dat, g, model, distr, ncov, itmax, epsilon)
    
    if (obj$error > 1)
      next
    
    counter <- counter + 1
    ret[counter, ] <- c(obj$pro, obj$mu, obj$sigma, obj$dof, obj$delta)
    
    if (counter >= B)
      break
  }
  std <- sqrt(apply(ret[1:counter, ], MARGIN = 2, FUN = "var"))
  names(std) <- dimnames(ret)[[2]]
  #std
  
  std_pi <- std[1 : g]
  std_mu <- matrix(std[g + 1 : (p * g)], nrow = p, ncol = g)
  std_sigma <- array(NA, c(p, p, g))
  for (i in 1 : g) 
    std_sigma[,, i] <- matrix(std[(g + p * g) + (p*p*(i-1)) + 1 : (p * p)], 
                              ncol = p)
  std_err <- list(std_pi = std_pi, std_mu = std_mu, std_sigma = std_sigma)
  
  if ((distr == "mvt") || (distr == "mst")) 
    std_err$std_dof <- std[((g + p * g) + (p * p * g)) + 1 : g] 

  if ((distr == "msn") || (distr == "mst"))
    std_err$std_delta <- matrix(std[((g + p * g) + (p * p * g) + g) + 
                                1 : (g * p)], nrow = p, ncol = g) 

  std_err 
}
