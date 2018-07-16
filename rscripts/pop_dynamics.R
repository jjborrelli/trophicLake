# Dynamic Models
## Libraries
library(deSolve)

# Rosenzweig-MacArthur

romac <- function(time, parms, state){
  with(as.list(c(state, parms)),{
    dx <- b * x * (1 - (x / K)) - ((c * x)/(a + x)) * y 
    dy <- e * ((c * x)/(a + x)) * y - d * y
    return(list(c(dx, dy)))
  })
}

par.RM <- list(
  b = 1.2,
  K = 2,
  a = .25, 
  c = .81, 
  e = .8,
  d = .25
)

# LKE model 
## Loladze et al. 2000

lke <- function(time, parms, state){
  with(as.list(c(state, parms)), {
    dx <- b * x * min(1 - (x/K), 1 - (x/((P - theta * y)/q))) - ((c * x)/(a + x)) * y
    dy <- e * min(1, ((P - theta * y)/x)/theta) * ((c * x)/(a + x)) * y - d * y
    
    return(list(c(dx, dy)))
  })
}

par.LKE <- list(
  b = 1.2,
  K = 1,
  a = .25,
  theta = 0.03,
  q = 0.0038,
  P = 2,
  c = .81, 
  e = .8,
  d = .25
)

## LKE with stochasticity

lke.stoch <- function(time, parms, state){
  with(as.list(c(state, parms)), {
    P.t <- P[ceiling(time)]
    dx <- b[ceiling(time)] * x * min(1 - (x/K[ceiling(time)]), 1 - (x/((P.t - theta * y)/q))) - 
      ((c * x)/(a + x)) * y
    dy <- e * min(1, ((P.t - theta * y)/x)/theta) * ((c * x)/(a + x)) * y - 
      d[ceiling(time)] * y
    
    return(list(c(dx, dy)))
  })
}

par.LKEstoch <- list(
  b = rnorm(1001, 1.2, .12),
  K = seq(1, .1, length.out = 1001),
  a = .25,
  theta = 0.03,
  q = 0.0038,
  P = abs(rnorm(1001, seq(2, .25, length.out = 1001), .0025)),
  c = .81, 
  e = .8,
  d = rnorm(1001, .25, .025)
)

# LKE 2004 
## 3 spp with predator competition

lke04 <- function(time, parms, state){
  with(as.list(c(state, parms)), {
    dx <- r * x * (1 - x/min(K, ((P.t - s1 * y1 - s2 * y2)/q))) - 
      ((c1 * x)/(a1 + x)) * y1 - ((c2 * x)/(a2 + x)) * y2
    
    dy1 <- e1 * min(1, (((P.t - s1 * y1 - s2 * y2)/x)/s1)) * 
      ((c1 * x)/(a1 + x)) * y1 - d1 * y1
    
    dy2 <- e2 * min(1, (((P.t - s1 * y1 - s2 * y2)/x)/s2)) * 
      ((c2 * x)/(a2 + x)) * y2 - d1 * y2
    
    return(list(c(dx, dy1, dy2)))
  })
}


par.LKE04 <- list(
  r = 0.93,
  K = .5, # 0 - 1.75
  c1 = 0.7,
  c2 = 0.8, 
  a1 = 0.3, 
  a2 = 0.2,
  e1 = 0.72,
  e2 = 0.76,
  d1 = 0.23,
  d2 = 0.2,
  s1 = 0.032,
  s2 = 0.05,
  q = 0.004,
  P.t = 0.03
)

# Peace13

peace13 <- function(time, parms, state){
  with(as.list(c(state, parms)), {
    dx <- b * x * (1 - x/min(K, (P - theta * y)/q)) - 
      min(((c * x)/(a + x)), (fhat * theta)/((P - theta * y)/x)) * y
    
    dy <- min(e * ((c * x)/(a + x)), 
              ((P - theta * y)/x)/theta * ((c * x)/(a + x)),
              e * fhat * (theta/((P - theta * y)/x))) * y - d * y
    
    return(list(c(dx, dy)))
  })
}


par.P13 <- list(
  b = 1.2,
  K = 1.5,
  a = 0.25,
  theta = 0.03,
  q = 0.0038,
  P = 0.05, # .03 - .2
  c = 0.81, 
  e = 0.8,
  d = 0.25
)

par.P13$fhat <- (par.P13$c * 1e8)/(par.P13$a + 1e8)


# Chen17

chen17 <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dx <- b * x * (1 - x/min(K, ((P.t - theta.y * y - theta.z * z)/q))) - 
      ((c.y * x)/(a.y + x)) * y
    
    dy <- e.y * min(1, ((P.t - theta.y * y - theta.z * z)/q)/ theta.y) * 
      ((c.y * x)/(a.y + x)) * y - ((c.z * y)/(a.z + y)) * z - d.y * y
    
    dz <- e.z * min(1, theta.y/theta.z) * ((c.z * y)/(a.z + y)) * z - d.z * z
    
    return(list(c(dx, dy, dz)))
  })
}


par.C17 <- list(
  b = 1.2,
  K = 2, # 0-10 
  P.t = 0.09,
  theta.y = 0.03,
  theta.z = 0.013,
  q = 0.0038,
  c.y = 0.81,
  a.y = 0.25,
  e.y = 0.8,
  d.y = 0.25,
  c.z = 0.03,
  a.z = 0.75,
  e.z = 0.75,
  d.z = 0.003
)
