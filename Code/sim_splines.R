# ------------------------------- generate splines
### B splines
make.Bsplines <- function(ordered_obs, start = 0, knots, knots_position, degree){
  ##### input ###
  ## ordered_obs : ordered observed event times.
  ## start : minimim of observed times.
  ## end : maximum of observed times.
  ## knots : the number of knots.
  ## knots_position : position of knots length of (knots + 1)
  ## degree : degree of the splines. (set degree = 3 for cubic splines)
  ##### output ###
  ## basis : B of length (knots + 4).
  
  n = length(ordered_obs) 
  degree = 3 # cubic basis function
  basis = matrix(0, nrow = length(ordered_obs), ncol =(knots + degree + 1))
  knotlist = rep(0, (knots + 8))
  knotlist[1:3] = seq(start, ordered_obs[1], ordered_obs[1] / 3)[1:3]
  knotlist[4] = ordered_obs[1]
  knotlist[(5:(4+ knots))] =  knots_position
  knotlist[(4+knots) + 1] = ordered_obs[n]
  knotlist[(4+knots + 2)] = ordered_obs[n] + 1 
  knotlist[(4+knots + 3)] = ordered_obs[n] + 2
  knotlist[(4+knots + 4)] = ordered_obs[n] + 3
  
  for(t in 1:length(ordered_obs)){
    for(k in 1:ncol(basis)){
      tmp <- 0
      for(ll in k:(k+4)){
        denom <- prod(c(knotlist[k:(k+4)] - knotlist[ll])[knotlist[k:(k+4)] - knotlist[ll] !=0] )
        tmp <- tmp + max( knotlist[ll] - ordered_obs[t], 0)^3 / denom
      }
      basis[t, k] <- (knotlist[(k+4)] - knotlist[k]) * tmp
    }
  }
  
  basis <- round(basis, 6)
  return(basis)
}

### I splines
make.Isplines <- function(ordered_obs, start, knots, knots_position, degree){
  ##### input ###
  ## ordered_obs : ordered observed event times.
  ## start : minimim of observed times.
  ## end : maximum of observed times.
  ## knots : the number of knots.
  ## knots_position : position of knots length of (knots + 1)
  ## degree : degree of the splines. (set degree = 3 for cubic splines)
  ##### output ###
  ## basis : I of length (knots + 4).
  
  n = length(ordered_obs) 
  degree <- 3
  basis <- matrix(0, nrow = length(ordered_obs), ncol =(knots + degree + 1))
  knotlist <- rep(0, (knots + 8))
  knotlist[1:3] <- seq(start, ordered_obs[1], ordered_obs[1] / 3)[1:3]
  knotlist[4] = ordered_obs[1]
  knotlist[(5:(4+ knots))] <- knots_position
  knotlist[(4+knots) + 1] <- ordered_obs[n]
  knotlist[(4+knots + 2)] <- ordered_obs[n] + 1 
  knotlist[(4+knots + 3)] <- ordered_obs[n] + 2
  knotlist[(4+knots + 4)] <- ordered_obs[n] + 3
  
  for(t in 1:length(ordered_obs)){
    for(k in 1:ncol(basis)){
      tmp1 <- 0 ; tmp2 <- 0
      for(ll in k:(k+4)){
        denom <- prod(c(knotlist[k:(k+4)] - knotlist[ll])[knotlist[k:(k+4)] - knotlist[ll] !=0] )
        tmp1 <- tmp1 + max(knotlist[ll], 0)^4 / denom
        tmp2 <- tmp2 +  max(knotlist[ll] - ordered_obs[t], 0)^4 / denom
      }
      basis[t, k] <- ((knotlist[(k+4)] - knotlist[k])/4) *  ( tmp1 - tmp2)
    }
  }
  basis <- round(basis, 6)
  return(basis)
}