## Functions for finding simulation parameters for variable-rate models

# We simulate with a speciation rate that is constant in the beginning and end and linear in the middle
# This gives us a function where we can control the magnitude of shift and the time over which it occurs
# Controlling the time window is also possible
# This model spans from a 2-epoch model with an instantaneous jump (radius=0), to a model with a decline over the whole tree (radius=tree age)
# It also includes the constant-rate case as a trivial case, but choosing parameters for those is much simpler 

# Parameterization is in terms of
# 1) Shift magnitude (rate_new/rate_old), so 2.0 is a two-fold increase, 0.5 a two-fold decrease
# 2) Shift center
# 3) Shift radius
# So the rate is rate[1] from t_inf to t_shift - radius, and changes linearly to rate[2], which proceeds until t_0
# Shorthanding radius as r, and using []/() notation, we have the intervals
#    [0,t_shift - r), [t_shift - r, t_shift + r) , [t_shift + r, infinity)
# This notation doesn't technically matter, since the actual rates are continuous, but it helps us define them in a way that accomodates radius=0 cases when simulating

# A function to find a value of lambda that produces trees of a given size for a constant-linear-constant setup
findParams <- function(t.shift,shift.radius,fold.change,extinction,target.taxa,sampling.fraction,tree.age=100) {
  # recover()
  mu <- function(x){return(extinction)}
  
  fn <- function(par) {
    rate_old <- par
    rate_new <- par * fold.change
    shift_delta <- rate_new - rate_old # absolute value of the change
    l <- function(x){
      if ( x < t.shift - shift.radius) {
        return(rate_old)
      } else if ( x >= t.shift + shift.radius ) {
        return(rate_new)
      } else {
        # Calculate rate by figuring out the proportion of the way along the shift the time is, multiplied by the magnitude of the shift
        prop <- (x - (t.shift - shift.radius))/(2 * shift.radius)
        return(rate_old + prop * shift_delta)
      }
    }
    nt <- tess.nTaxa.expected(begin=0,t=tree.age,end=tree.age,lambda=l,mu=mu,samplingProbability=rho)
    return((nt-target.taxa)^2)
  }
  # Avoid optimizer getting lost out at infinity
  grid_fn <- sapply(seq(-10,1,1),function(i){fn(10^i)})
  ub <- 10^(seq(-10,1,1)[max(which(is.finite(grid_fn)))])
  opts <- optimize(fn,interval=c(0,ub))
  rate_old <- opts$minimum
  rate_new <- fold.change * rate_old
  res <- c(rate_old,rate_new)
  names(res) <- c("rate.old","rate.new")
  return( res )
}

# A function for obtaining the speciation rate at any point, for plotting and simulating
speciationRate <- function(t,t.shift,shift.radius,rate.old,rate.new) {
  if ( t < t.shift - shift.radius) {
    return(rate.old)
  } else if ( t >= t.shift + shift.radius ) {
    return(rate.new)
  } else {
    # Calculate rate by figuring out the proportion of the way along the shift the time is, multiplied by the magnitude of the shift
    prop <- (t - (t.shift - shift.radius))/(2 * shift.radius)
    return(rate.old + prop * (rate.new - rate.old))
  }
}
