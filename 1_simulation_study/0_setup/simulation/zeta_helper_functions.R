##############################################################
# Helper functions and functions for wrangling the horseshoe #
##############################################################

# The density of Faulkner and Minin's approximation to the Horseshoe
dApproxHorseshoe <- function(x, gamma=1) {
  return((2 * pi^3 * gamma^2)^(-1/2) * ( (sqrt(pi) - 2)/(2 * sqrt(2) - 4) * log(1 + (4*gamma^2)/(x^2)) + (sqrt(2) - sqrt(pi))/(sqrt(2) - 2) * log(1 + (2*gamma^2)/(x^2)) ))
}

# Inverse cotangent (for computing integral of approximate density)
acot <- function(x) {
  if (x < 0) {
    return(-0.5 * pi - atan(x))
  } else {
    return(0.5 * pi - atan(x))
  }
}

# Integral of approximate HS density, for computing tail probabilities
FApproxHorseshoe <- function(x, gamma=1) {
  # Thanks to Wolfram-Alpha for integrating this
  (x * (2 * (sqrt(2) - sqrt(pi)) * log((2*gamma^2) / (x^2) + 1) + (sqrt(pi) - 2) * log((4*gamma^2) / (x^2) + 1)) - 4 * (sqrt(2*pi) - 2) * gamma * atan(x / (gamma*sqrt(2))) + 4 * (sqrt(pi) - 2) * gamma * acot(2 * gamma / x)) /
    (2 * sqrt(2) * (sqrt(2) - 2) * pi^(3/2) * sqrt(gamma^2))
}

# Right tail probability of the appriximate HS
# Computes F(Inf) - F(x), where Inf is machine's maximum double
pRightTailApproxHorseshoe <- function(x, gamma=1) {
  FApproxHorseshoe(.Machine$double.xmax,gamma=gamma) - FApproxHorseshoe(x,gamma=gamma)
}

# Right tail probability of the appriximate HS
# Integrates the density function
pRightTailHorseshoeIntegral <- function(x, gamma=1) {
  fn <- function(y) {
    pnorm(x,0,y,lower.tail=FALSE)*(dcauchy(y,0,gamma)/0.5)
  }
  integrate(fn,x,Inf,subdivisions=10000,abs.tol=1e-5)$value
}

# Right tail probability of the HS
# Integrates the density function via grid, so is slow
# This produces values that most closely match those simulated in R
# This does not guarantee correctness, but neither does using the approximation or calling integrate()
# This way we at least avoid headaches when trying to match our analytical results to simulations
pRightTailHorseshoeGrid <- function(x, gamma=1, grid.size=5000) {
  quants <- seq(1e-10,1-1e-10,length.out=grid.size)
  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants)/2.0
  probs <- 1/length(quants) # we're using quantiles, each gamma is equally likely
  sigmas <- qcauchy(quants,0,gamma)
  sum(pnorm(x,0,sigmas,lower.tail=FALSE) * probs)
}

# Draw a random HS(gamma) variable using hierarchical formulation
rHS <- function(n,gamma=1) {
  sigma <- abs(rcauchy(n,0,gamma))
  theta <- rnorm(n,0,sigma)
  return(theta)
}

# Simulate a trajectory from a 1-D HSMRF, given the number of cells and the global shrinkage parameter
simulateHSTrajectory <- function(zeta,n.cells) {
  # Draw gamma
  gamma <- abs(rcauchy(1,0,zeta))
  return(cumsum(rHS(n.cells-1,gamma)))
}

# Simulate and plot a realization of a HS trajectory
# log.scale = TRUE if the process is on the log scale of the real variables
plotSimulatedHSTrajectory <- function(zeta,n.cells,log.scale=TRUE) {
  gamma <- abs(rcauchy(1,0,zeta))
  traj <- c(0,cumsum(rHSzeta(n.cells-1,gamma)))
  if ( log.scale ) {
    traj <- exp(traj)
  }
  plot(x=c(1:n.cells),y=traj,xlab="t",ylab="X(t)",type="l")
}

# Simulate a trajectory from a 1-D GMRF, given the number of cells and the global shrinkage parameter
simulateBMTrajectory <- function(zeta,n.cells) {
  # Draw sigma
  sigma <- abs(rcauchy(1,0,zeta))
  return(cumsum(rnorm(n.cells-1,0,sigma)))
}

# Simulate and plot a realization of a HS trajectory
# log.scale = TRUE if the process is on the log scale of the real variables
plotSimulatedBMTrajectory <- function(zeta,n.cells,log.scale=TRUE) {
  sigma <- abs(rcauchy(1,0,zeta))
  traj <- c(0,cumsum(rHSzeta(n.cells-1,sigma)))
  if ( log.scale ) {
    traj <- exp(traj)
  }
  plot(x=c(1:n.cells),y=traj,xlab="t",ylab="X(t)",type="l")
}

# Calculate total variation norm of a trajectory
TVN <- function(traj) {
  sum(abs(traj[2:length(traj)] - traj[1:(length(traj)-1)]))
}

