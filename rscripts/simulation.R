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
timesteps <- 500


# set up array to store biomass of prey
lake.x <- array(0, dim = c(ydim, xdim, timesteps))
# set up array to store biomass of predator
lake.y <- array(0, dim = c(ydim, xdim, timesteps))
# set up array to store biomass of top pred
lake.z <- array(0, dim = c(ydim, xdim, timesteps))
# set up lake phosphorus concentration
# gradient
lake.P <- array(0, dim = c(ydim, xdim, timesteps))
lake.P[,,1] <- matrix(sort(runif(xdim*ydim, .08, .08)), nrow = ydim, ncol = xdim)
lake.P[,,1] <- lake.P[sample(1:nrow(lake.P)),,1]
#cols <- rev(sort(sample(1:ncol(lake.P), ncol(lake.P)/2)))
#lake.P <- (lake.P[,c(cols, (1:nrow(lake.P))[-cols])])
# initial biomass for both species
lake.x[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
lake.y[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
#lake.z[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
lake.z[,,1] <- matrix((lake.P[,,1] - par.P14fixed$theta * lake.y[,,1])/lake.x[,,1])

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
      
      
      par.P14fixed$P <- lake.P[i,j,t]
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
  
  lake.x[,,t+1][lake.x[,,t+1] < 0] <- 0
  lake.y[,,t+1][lake.y[,,t+1] < 0] <- 0

  ## Zoop Migration
  mig <- migration(lake.y[,,t+1], lake.x[,,t+1], lake.z[,,t+1], lake.P[,,t], 
                   par.P14fixed$theta,
                             max.move = 1, prop.migrant = .1, method = "quota")
  lake.y[,,t+1] <- mig[[1]]
  lake.P[,,t+1] <- mig[[2]]
  
  if(((t + 1) %% 10) == 0){cat("Day ", t+1, "of ", timesteps, "simulated\n")}
}
t1 <- Sys.time()
t1-t0


#run1 <- list(lake.x, lake.y, lake.z, lake.P) # random
#run2 <- list(lake.x, lake.y, lake.z, lake.P) # density
#run3 <- list(lake.x, lake.y, lake.z, lake.P) # quality 
#run4 <- list(lake.x, lake.y, lake.z, lake.P) # quota

run1.1 <- list(lake.x, lake.y, lake.z, lake.P)
#run2.1 <- list(lake.x, lake.y, lake.z, lake.P)
#run3.1 <- list(lake.x, lake.y, lake.z, lake.P)

#run1.2 <- list(lake.x, lake.y, lake.z, lake.P)
run2.2 <- list(lake.x, lake.y, lake.z, lake.P)
#run3.2 <- list(lake.x, lake.y, lake.z, lake.P)


###############################
# Visualize results
test <- run1.2
threePlot(3, test[[3]], test[[2]], test[[1]])
threePlot(5, lake.z[,,], lake.y[,,], lake.x[,,])

#par(mfrow = c(2,1), mar = c(1,2,1,.2))

cell <- 5
iter <- 90
pheatmap(apply(lake.x[,,iter], c(1,2), mean), cluster_rows = F, cluster_cols = F, cellwidth = cell, cellheight = cell)
pheatmap(apply(lake.y[,,iter], c(1,2), mean), cluster_rows = F, cluster_cols = F, cellwidth = cell, cellheight = cell)
pheatmap(apply(lake.z[,,iter], c(1,2), mean), cluster_rows = F, cluster_cols = F, cellwidth = cell, cellheight = cell)

pheatmap(apply(lake.P[,,iter], c(1,2), mean), cluster_rows = F, cluster_cols = F, cellwidth = cell, cellheight = cell)

#pheatmap(lake.P, cluster_rows = F, cluster_cols = F, cellwidth = 5, cellheight = 5)

#saveRDS(list(lake.x, lake.y, lake.z, lake.P), "sides_gradientP_run_100.rds")
#run <- readRDS("gradientP_run_100.rds")

lake.y[1,1,]


all2 <- cbind(apply(run1[[1]], 3, mean), apply(run1[[2]], 3, mean))
matplot(all2, typ = "l")
all2 <- cbind(apply(run2[[1]], 3, mean), apply(run2[[2]], 3, mean))
matplot(all2, typ = "l")
all2 <- cbind(apply(run3[[1]], 3, mean), apply(run3[[2]], 3, mean))
matplot(all2, typ = "l")

all2 <- cbind(apply(run1.1[[1]], 3, mean), apply(run1.1[[2]], 3, mean))
matplot(all2, typ = "l")
all2 <- cbind(apply(run2.1[[1]], 3, mean), apply(run2.1[[2]], 3, mean))
matplot(all2, typ = "l")
all2 <- cbind(apply(run3.1[[1]], 3, mean), apply(run3.1[[2]], 3, mean))
matplot(all2, typ = "l")
#matplot(all2)

spatdynplot(lake.x, t = timesteps, interval = .2, name = "xlake.html")
spatdynplot(lake.y, t = timesteps, interval = .2, name = "ylake.html")
dev.off()

matplot(t(apply(lake.x, 3, colMeans)), typ = "l")
apply(lake.y, 3, colMeans)

matplot(cbind(t(lake.P[1,,]),t(lake.P[2,,]),t(lake.P[3,,]),t(lake.P[4,,])), typ = "l")



library(reshape2)
library(ggplot2)
library(gganimate)

lx <- data.frame(melt(lake.x), lake = "producer")
ly <- data.frame(melt(lake.y), lake = "grazer")

ggplot(dplyr::filter(rbind(lx, ly), Var3 %in% 1:100), aes(x = Var1, y = Var2, fill = value)) + geom_tile() + facet_wrap(~lake) + 
  transition_states(Var3,transition_length = 1, state_length = 2)
