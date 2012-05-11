fwsim <- function(g = 10, k = 20000, r = 7, alpha = 1, mu = 0.003, save.gs = NULL, trace = TRUE, ...) { 
  if (!is.numeric(g) || length(g) != 1) stop("The number of generations, g, must be an integer")
  if (!is.numeric(k) || length(k) != 1) stop("The initial population size, k, must be an integer")
  if (!is.numeric(r) || length(r) != 1) stop("The number of loci, r, must be an integer")
  
  if (!is.numeric(alpha) || (length(alpha) != 1 && length(alpha) != g)) 
    stop("The growth rate, alpha, must be a numeric vector of size 1 or g")

  if (!is.logical(trace) || length(trace) != 1) stop("trace must be a logical")
  
  if (!is.numeric(mu)) stop("The mutation rate, mu, must be a numeric")
  
  if (is.matrix(mu)) {
    if (!(nrow(mu) == 2 && ncol(mu) == r)) {
      stop("The mutation is specified as a matrix, but with the wrong dimensions. Expected a 2xr matrix.")
    }
  } else if (length(mu) == r) {
    mu <- matrix(rep(mu, 2)/2, nrow = 2, byrow = TRUE)
  } else if (length(mu) == 1) {
    mu <- matrix(rep(mu, 2*r)/2, nrow = 2, byrow = TRUE)
  } else {
    stop("The mutation specification cannot be inferred. Please read the manual at ?fwsim.")
  }

  if (any(mu < 0 | mu > 1)) {
    stop("Please only use mutation rates between 0 and 1.")
  }
  
  if (!is.null(save.gs) && length(save.gs) > 0) {
    if (!is.numeric(save.gs)) stop("save.gs must be an integer.")
    save.gs <- sort(round(save.gs))
    if (any(save.gs <= 0 | save.gs >= g)) stop("If not null, save.gs must be numbers between 0 and g, both excluded.")
    
    new.gs <- rep(0, g-1)
    new.gs[save.gs] <- 1
    save.gs <- c(-1, new.gs, -1)
  } else {
    save.gs <- rep(-1, g+1) # For easier handling in C
  }
  
  # Let's let the GC clean up before the run
  gc()
  
  ans <- .Call("fwsim", 
    as.integer(g), 
    as.integer(k), 
    as.integer(r), 
    as.numeric(alpha), 
    as.numeric(mu), 
    as.integer(save.gs), 
    as.integer(trace),
    as.integer(length(alpha)),
    PACKAGE = "fwsim")
  
  #print(ans$haplotypes)
  if (!is.null(ans$haplotypes)) {
    ans$haplotypes <- data.frame(ans$haplotypes)  
    colnames(ans$haplotypes)[1:r] <- paste("Locus", 1:r, sep = "")
    colnames(ans$haplotypes)[r + 1] <- "N"
  } else {
    warning("Population died out")
  }

  if (!is.null(ans$intermediate.haplotypes)) {
    ans$intermediate.haplotypes <- lapply(ans$intermediate.haplotypes, function(l) {
      if (is.null(l)) return(NULL)
      
      newl <- data.frame(l)
      colnames(newl)[1:r] <- paste("Locus", 1:r, sep = "")
      colnames(newl)[r + 1] <- "N"
      return(newl)
    })
  }
    
  if ((is.null(ans$haplotypes) && ans$sizes[g+1] != 0) || 
      (!is.null(ans$haplotypes) && sum(ans$haplotypes$N) != ans$sizes[g+1])) {
    
    #if (save.gs[2] != -1) {
    #  is <- tail(which(save.gs == 1))
    #  print(ans$intermediate.haplotypes[is])
    #  print(lapply(ans$intermediate.haplotypes[is], function(hap) sum(hap$N)))
    #}

    #print(ans$haplotypes)
        
    #print(tail(ans$sizes))
    stop(sum(ans$haplotypes$N), " = sum(ans$haplotypes$N) != ans$sizes[g+1] = ", ans$sizes[g+1])
  }
  
  # Let's let the GC clean up again
  gc()
  
  return(ans)
}

