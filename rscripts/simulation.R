# Simulate

# get pop dynamics model functions
source("rscripts/pop_dynamics.R")

# get migration functions
source("rscripts/migration.R")


###############################
# Set up the lake

# Horizontal dimension
xdim <- 5
# Vertical dimension
ydim <- 5
# Number of simulation timesteps (days)
timesteps <- 120


# set up array to store biomass of prey
lake.x <- array(0, dim = c(ydim, xdim, timesteps))
# set up array to store biomass of predator
lake.y <- array(0, dim = c(ydim, xdim, timesteps))
# set up array to store biomass of top pred
lake.z <- array(0, dim = c(ydim, xdim, timesteps))
# set up lake phosphorus concentration
lake.P <- matrix(runif(xdim*ydim, 0.03, 0.2), nrow = ydim, ncol = xdim)

# initial biomass for both species
lake.x[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
lake.y[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
#lake.z[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
lake.z[,,1] <- matrix((lake.P - par.P14fixed$theta * lake.y[,,1])/lake.x[,,1])

###############################
# Simulate dynamics

t0 <- Sys.time()
for(t in 1:(timesteps - 1)){
  for(i in 1:nrow(lake.x)){
    for(j in 1:ncol(lake.x)){
      
      state = c(x = lake.x[i,j,t], y = lake.y[i,j,t], Q = lake.z[i,j,t])
      
      out <- ode(state, 1:2, P14_fixed, par.P14fixed)
      
      lake.x[i,j,t+1] <- out[2, 2]
      lake.y[i,j,t+1] <- out[2, 3]
      lake.z[i,j,t+1] <- out[2, 4]
      
    }
  }
  
  # pick migration
  ## Random
  #lake.x <- migration.random(lake.x, t)
  #lake.y <- migration.random(lake.y, t)
  
  ## Density
  lake.y[,,t+1] <- migration.density(lake.y[,,t], lake.x[,,t], max.move = 1, prop.migrant = .1)
  #lake.z[,,t+1] <- migration.density(lake.z[,,t], lake.y[,,t], max.move = 2, prop.migrant = .1)
  
  if(((t + 1) %% 5) == 0){cat("Day ", t+1, "of ", timesteps, "simulated\n")}
}
t1 <- Sys.time()
t1-t0


###############################
# Visualize results

par(mfrow = c(3,1), mar = c(1,2,1,.2))
matplot(t(lake.z[,1,]), typ = "l")
matplot(t(lake.y[,1,]), typ = "l")
matplot(t(lake.x[,1,]), typ = "l")

pheatmap::pheatmap(lake.x[,,180], cluster_rows = F, cluster_cols = F, cellwidth = 20, cellheight = 20)

lake.y[1,1,]

all <- cbind(apply(lake.y, 3, sum), apply(lake.z, 3, sum))
matplot(all, typ = "l")
