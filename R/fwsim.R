fwsim <- function(g = 10, k = 20000, r = 7, alpha = 1, mu = 0.003, mutmodel = 1, trace = TRUE, trace.loc.mut = 4, ...) { 
  ans <- .Call("fwsim", 
    as.integer(g), 
    as.integer(k), 
    as.integer(r), 
    as.numeric(alpha), 
    as.numeric(mu), 
    as.integer(trace),
    as.integer(trace.loc.mut),
    PACKAGE = "fwsim")
  
  ans$haplotypes <- data.frame(ans$haplotypes)
  
  newr <- ncol(ans$haplotypes) - 1
  
  colnames(ans$haplotypes)[1:newr] <- paste("Locus", 1:newr, sep = "")
  colnames(ans$haplotypes)[newr + 1] <- "N"
  
  if (sum(ans$haplotypes$N) != ans$sizes[g+1]) stop("sum(ans$haplotypes$N) != tail(ans$sizes, 1)")
  
  return(ans)
}

