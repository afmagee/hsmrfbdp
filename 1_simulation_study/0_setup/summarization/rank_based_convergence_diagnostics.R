# Compute the stan crew's new rank-based convergence diagnostics
# Input: list of data.frames, matrices, or mcmc class objects of MCMC output
# Output: depends on return.both
#         return.both=TRUE,  returns both the folded and unfolded rank-based splitPSRF measures as named matrix
#         return.both=FALSE, returns the maximum of each diagnostic per parameter as named matrix with one column

diagnoseConvergence <- function(chains, return.both=TRUE) {
  # Naming of variables follows Gelman et al. (2013)
  
  m <- 2*length(chains)
  mcmcdim <- do.call(rbind,lapply(chains,dim))
  
  # We're assuming these are all runs of equal length of the same model, all MCMC outputs should be the same size
  if ( !all(apply(mcmcdim,2,function(x){length(unique(x))}) == 1) ) {
    stop("MCMC logs are not of the same size")
  }
  
  if ( mcmcdim[1,1] %% 2 == 1 ) {
    # Chains have odd numbers of samples, n will be variable and we must assign the odd sample to one half of each chain
    # n_vec <- rep(c(floor(mcmcdim[1,1]/2),ceiling(mcmcdim[1,1]/2)),2)
    warning("Chains have odd numbers of samples, discarding last sample in calculation of diagnostics")
    chains <- lapply(chains,function(x){x[-dim(x)[1],]})
    n <- dim(chains[[1]])[1]/2
  } else {
    # n_vec <- rep(mcmcdim[,1]/2,2)
    n <- mcmcdim[1,1]/2
  }
  ends <- cumsum(rep(n,m))
  starts <- 1 + cumsum(c(0,rep(n,m-1)))

  # recover()
  ## Concatenated MCMC for ranks
  all_chains <- do.call(rbind,chains)
  
  ## Compute rank-based R hat
  unfolded_ranks <- apply(all_chains,2,rank,ties.method="average")
  
  # Split chains in half
  chains <- lapply(1:m,function(i){unfolded_ranks[starts[i]:ends[i],]})
  
  # Compute some stuff!
  # Chains now in columns, variables in rows
  psi_bar <- do.call(cbind,lapply(chains,function(x){colMeans(x)}))
  s_sq <- do.call(cbind,lapply(chains,function(x){apply(x,2,var)}))
  
  W <- rowSums(s_sq)/m
  B <- (n/(m-1)) * rowSums((psi_bar - rowMeans(psi_bar))^2)
  
  v_hat <- (n-1)/n*W + 1/n*B
  
  zsRhat <- sqrt(v_hat/W)
  
  ## Compute folded rank-based R hat
  folded_ranks <- apply(apply(all_chains,2,function(x){abs(x - median(x))}),2,rank,ties.method="average")
  
  # Split chains in half
  chains <- lapply(1:m,function(i){folded_ranks[starts[i]:ends[i],]})
  
  # Compute some stuff!
  # Chains now in columns, variables in rows
  psi_bar <- do.call(cbind,lapply(chains,function(x){colMeans(x)}))
  s_sq <- do.call(cbind,lapply(chains,function(x){apply(x,2,var)}))
  
  W <- rowSums(s_sq)/m
  B <- (n/(m-1)) * rowSums((psi_bar - rowMeans(psi_bar))^2)
  
  v_hat <- (n-1)/n*W + 1/n*B
  
  fzsRhat <- sqrt(v_hat/W)
  # recover()

  both <- cbind(zsRhat,fzsRhat)
  
  if ( return.both ) {
    return(both)
  } else {
    return(cbind(apply(both,1,max)))
  }
}

rankESS <- function(chains) {
  # recover()
  ## Concatenated MCMC for ranks
  all_chains <- do.call(rbind,chains)
  
  unfolded_ranks <- apply(all_chains,2,rank,ties.method="average")
  
  ess <- effectiveSize(unfolded_ranks)
  
  names(ess) <- colnames(chains[[1]])
  
  return(ess)  
}
