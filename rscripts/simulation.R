# Simulate

# get pop dynamics model functions
source("rscripts/pop_dynamics.R")

# get migration functions
source("rscripts/migration.R")

# get plotting functions
source("rscripts/plotting.R")


###############################
# Set up the lake

# Horizontal dimension
xdim <- 10
# Vertical dimension
ydim <- 10
# Number of simulation timesteps (days)
timesteps <- 300


# set up array to store biomass of prey
lake.x <- array(0, dim = c(ydim, xdim, timesteps))
# set up array to store biomass of predator
lake.y <- array(0, dim = c(ydim, xdim, timesteps))
# set up array to store biomass of top pred
lake.z <- array(0, dim = c(ydim, xdim, timesteps))
# set up lake phosphorus concentration
# gradient
lake.P <- matrix((runif(xdim*ydim, .03, 1)), nrow = ydim, ncol = xdim)
lake.P <- lake.P[sample(1:nrow(lake.P)),]

# initial biomass for both species
lake.x[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
lake.y[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
#lake.z[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
lake.z[,,1] <- matrix((lake.P - par.P14fixed$theta * lake.y[,,1])/lake.x[,,1])

###############################
# Simulate dynamics

errors.t <- vector(length = timesteps)
t0 <- Sys.time()
for(t in 1:(timesteps - 1)){
  errors.i <- matrix(nrow = nrow(lake.x), ncol = 3)
  for(i in 1:nrow(lake.x)){
    for(j in 1:ncol(lake.x)){
      lake.x[,,t][lake.x[,,t] <= 10^-10] <- 0
      lake.y[,,t][lake.y[,,t] <= 10^-10] <- 0
      lake.x[,,t][lake.x[,,t] == 0] <- 10^-4
      
      
      par.P14fixed$P <- lake.P[i,j]
      state = c(x = lake.x[i,j,t], y = lake.y[i,j,t], Q = lake.z[i,j,t])
      
      out <- ode(state, 1:2, P14_fixed, par.P14fixed)
      
      errors <- c(is.nan(out[2,2]), is.nan(out[2,3]), is.nan(out[2,4])) 
      
      lake.x[i,j,t+1] <- ifelse(is.nan(out[2, 2]), 10^-6, out[2,2])
      lake.y[i,j,t+1] <- ifelse(is.nan(out[2, 3]), 0, out[2,3])
      lake.z[i,j,t+1] <- ifelse(is.nan(out[2, 4]), 0, out[2,4])
      
    }
    errors.i[i, ] <- errors
  }
  errors.t[t] <- any(is.na(errors.i))
  
  # pick migration
  ## Random
  #lake.x <- migration.random(lake.x, t)
  #lake.y[,,t+1] <- migration.random(lake.y[,,t+1])
  
  ## Density
  lake.y[,,t+1] <- migration.density(lake.y[,,t+1], lake.x[,,t+1], max.move = 1, prop.migrant = .1)
  #lake.x[,,t+1] <- migration.random(lake.x[,,t+1], max.move = 1)
  #lake.P <- migration.random(lake.P, max.move = 1)
  
  if(((t + 1) %% 10) == 0){cat("Day ", t+1, "of ", timesteps, "simulated\n")}
}
t1 <- Sys.time()
t1-t0


###############################
# Visualize results

threePlot(3, lake.z, lake.y, lake.x)

#par(mfrow = c(2,1), mar = c(1,2,1,.2))
pheatmap(apply(lake.x[,,250:timesteps], c(1,2), mean), cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15)
pheatmap(apply(lake.y[,,250:timesteps], c(1,2), mean), cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15)
pheatmap(apply(lake.z[,,250:timesteps], c(1,2), mean), cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 15)

pheatmap(lake.P, cluster_rows = F, cluster_cols = F, cellwidth = 5, cellheight = 5)



lake.y[1,1,]

all2 <- cbind(apply(lake.x, 3, sum), apply(lake.y, 3, sum))
matplot(all2, typ = "l")
matplot(all)

spatdynplot(lake.x, t = timesteps, interval = .2, name = "xlake.gif")
spatdynplot(lake.y, t = timesteps, interval = .2, name = "ylake.gif")
dev.off()

matplot(t(apply(lake.x, 3, colMeans)), typ = "l")
apply(lake.y, 3, colMeans)

