fwsim <- function(g = 10, k = 20000, r = 7, alpha = 1, mu = 0.003, trace = TRUE, ...) { 
  if (!is.numeric(g) || length(g) != 1) stop("The number of generations, g, must be an integer")
  if (!is.numeric(k) || length(k) != 1) stop("The initial population size, k, must be an integer")
  if (!is.numeric(r) || length(r) != 1) stop("The number of loci, r, must be an integer")
  
  if (!is.numeric(alpha) || (length(alpha) != 1 && length(alpha) != g)) 
    stop("The growth rate, alpha, must be a numeric vector of size 1 or g")

  if (!is.logical(trace) || length(trace) != 1) stop("trace must be a logical")
  
  if (!is.numeric(mu)) stop("The mutation rate, mu, must be a numeric")
  
  if (is.matrix(mu)) {
    if (!(nrow(mu) == 2 && ncol(mu) == r)) {
      stop("The mutation is specified as a matrix, but with the wrong dimensions. Expected a 2xr matrix. ")
    }
  } else if (length(mu) == r) {
    mu <- matrix(rep(mu, 2)/2, nrow = 2, byrow = TRUE)
  } else if (length(mu) == 1) {
    mu <- matrix(rep(mu, 2*r)/2, nrow = 2, byrow = TRUE)
  } else {
    stop("The mutation specification is not understood. Please read the manual at ?fwsim.")
  }

  if (any(mu < 0 | mu > 1)) {
    stop("Please only use mutation rates between 0 and 1.")
  }
  
  ans <- .Call("fwsim", 
    as.integer(g), 
    as.integer(k), 
    as.integer(r), 
    as.numeric(alpha), 
    as.numeric(mu), 
    as.integer(trace),
    as.integer(length(alpha)),
    PACKAGE = "fwsim")
  
  ans$haplotypes <- data.frame(ans$haplotypes)
  
  newr <- ncol(ans$haplotypes) - 1
  
  colnames(ans$haplotypes)[1:newr] <- paste("Locus", 1:newr, sep = "")
  colnames(ans$haplotypes)[newr + 1] <- "N"
  
  if (sum(ans$haplotypes$N) != ans$sizes[g+1]) stop("sum(ans$haplotypes$N) != ans$sizes[g+1]")
  
  return(ans)
}

