source("1_simulation_study/0_setup/simulation/zeta_helper_functions.R")

# prior.n.shifts: the desired prior expectation of the number of shifts
# ncell: the number of cells/time windows in the BDP (1 + number of possible shift points)
# U: the definition of a shift is theta < -U or theta > U
# jump=log(x) corresponds to an x-fold increase or decrease of the variable, assuming the HSMRF is on the log-scale
# We treat the change between each grid cell as a Bernoulli RV, so the collection of changes becomes binomial
# From this we can calculate the expected number of cells where a jump occurs
setZetaHSMRFExpectedNumberOfJumps <- function(prior.n.shifts,ncell=100,jump=log(2)) {
  # recover()
  # Probability of a jump for a value of zeta
  # We average the conditional p(jump | gamma) over p(gamma)
  quants <- seq(0.0001,0.9999,length.out=2000)
  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants)/2.0
  probs <- 1/length(quants) # we're using quantiles, each gamma is equally likely
  # Function to optimize
  fn <- function(zeta) {
    # Grid of gammas
    gammas <- qcauchy(quants,0,zeta)
    # Number of expected jumps for each value of sigma
    num_expected_jumps <- sapply(gammas,function(x) {
      p_jump_one_cell_this_gamma <- pRightTailHorseshoeGrid(jump,x,grid.size=2000)/0.5
      return(p_jump_one_cell_this_gamma * (ncell-1))
    })
    # Average the per-sigma E(n_jumps) over p(sigma) to get overall expectation given zeta
    this_expected_num_jumps <- sum(probs * num_expected_jumps)
    return( (log(this_expected_num_jumps) - log(prior.n.shifts))^2 ) # Distance to target
  }
  # Find best value of zeta
  opts <- optimize(fn,c(0,1))
  zeta <- opts$minimum
  # Compute the prior on number of jumps for this zeta (to show user how well we approximated the target)
  gammas <- qcauchy(quants,0,zeta)
  num_expected_jumps <- sapply(gammas,function(x) {
    p_jump_one_cell_this_gamma <- pRightTailHorseshoeGrid(jump,x,grid.size=2000)/0.5
    return(p_jump_one_cell_this_gamma * (ncell-1))
  })
  computed_num_expected_jumps <- sum(probs * num_expected_jumps)
  return(list(zeta=zeta,E.n=computed_num_expected_jumps))
}

# prior.n.shifts: the desired prior expectation of the number of shifts
# ncell: the number of cells/time windows in the BDP (1 + number of possible shift points)
# U: the definition of a shift is theta < -U or theta > U
# jump=log(x) corresponds to an x-fold increase or decrease of the variable, assuming the HSMRF is on the log-scale
# We treat the change between each grid cell as a Bernoulli RV, so the collection of changes becomes binomial
# From this we can calculate the expected number of cells where a jump occurs
setZetaHSMRFExpectedNumberOfJumps <- function(prior.n.shifts,ncell=100,jump=log(2)) {
  # recover()
  # Probability of a jump for a value of zeta
  # We average the conditional p(jump | gamma) over p(gamma)
  quants <- seq(0.0001,0.9999,length.out=2000)
  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants)/2.0
  probs <- 1/length(quants) # we're using quantiles, each gamma is equally likely
  # Function to optimize
  fn <- function(zeta) {
    # Grid of gammas
    gammas <- qcauchy(quants,0,zeta)
    # Number of expected jumps for each value of sigma
    num_expected_jumps <- sapply(gammas,function(x) {
      p_jump_one_cell_this_gamma <- pRightTailHorseshoeGrid(jump,x,grid.size=2000)/0.5
      return(p_jump_one_cell_this_gamma * (ncell-1))
    })
    # Average the per-sigma E(n_jumps) over p(sigma) to get overall expectation given zeta
    this_expected_num_jumps <- sum(probs * num_expected_jumps)
    return( (log(this_expected_num_jumps) - log(prior.n.shifts))^2 ) # Distance to target
  }
  # Find best value of zeta
  opts <- optimize(fn,c(0,1))
  zeta <- opts$minimum
  # Compute the prior on number of jumps for this zeta (to show user how well we approximated the target)
  gammas <- qcauchy(quants,0,zeta)
  num_expected_jumps <- sapply(gammas,function(x) {
    p_jump_one_cell_this_gamma <- pRightTailHorseshoeGrid(jump,x,grid.size=2000)/0.5
    return(p_jump_one_cell_this_gamma * (ncell-1))
  })
  computed_num_expected_jumps <- sum(probs * num_expected_jumps)
  return(list(zeta=zeta,E.n=computed_num_expected_jumps))
}

# prior.n.shifts: the desired prior expectation of the number of shifts
# ncell: the number of cells/time windows in the BDP (1 + number of possible shift points)
# U: the definition of a shift is theta < -U or theta > U
# jump=log(x) corresponds to an x-fold increase or decrease of the variable, assuming the HSMRF is on the log-scale
# We treat the change between each grid cell as a Bernoulli RV, so the collection of changes becomes binomial
# From this we can calculate the expected number of cells where a jump occurs
setZetaGMRFExpectedNumberOfJumps <- function(prior.n.shifts,ncell=100,jump=log(2)) {
  # recover()
  # Probability of a jump for a value of zeta
  # We average the conditional p(jump | sigma) over p(sigma)
  quants <- seq(0.0001,0.9999,length.out=2000)
  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants)/2.0
  probs <- 1/length(quants) # we're using quantiles, each gamma is equally likely
  # Function to optimize
  fn <- function(zeta) {
    # Grid of sigmas
    sigmas <- qcauchy(quants,0,zeta)
    # Number of expected jumps for each value of sigma
    num_expected_jumps <- sapply(sigmas,function(x) {
      p_jump_one_cell_this_sigma <- pnorm(jump,0,x,lower.tail=FALSE)/0.5
      return(p_jump_one_cell_this_sigma * (ncell-1))
    })
    # Average the per-sigma E(n_jumps) over p(sigma) to get overall expectation given zeta
    this_expected_num_jumps <- sum(probs * num_expected_jumps)
    return( (log(this_expected_num_jumps) - log(prior.n.shifts))^2 ) # Distance to target
  }
  # Find best value of zeta
  opts <- optimize(fn,c(0,1))
  zeta <- opts$minimum
  # Compute the prior on number of jumps for this zeta (to show user how well we approximated the target)
  sigmas <- qcauchy(quants,0,zeta)
  num_expected_jumps <- sapply(sigmas,function(x) {
    p_jump_one_cell_this_sigma <- pnorm(jump,0,x,lower.tail=FALSE)/0.5
    return(p_jump_one_cell_this_sigma * (ncell-1))
  })
  computed_num_expected_jumps <- sum(probs * num_expected_jumps)
  return(list(zeta=zeta,E.n=computed_num_expected_jumps))
}