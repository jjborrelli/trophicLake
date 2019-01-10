# Lake set up

make_lakes <- function(xdim, ydim, timesteps, pmin = 0.8, pmax = 0.8, theta, gradient = FALSE){
  
  # set up array to store biomass of prey
  lake.x <- array(0, dim = c(xdim, ydim, timesteps))
  # set up array to store biomass of predator
  lake.y <- array(0, dim = c(xdim, ydim, timesteps))
  # set up array to store biomass of top pred
  lake.z <- array(0, dim = c(xdim, ydim, timesteps))
  # set up lake phosphorus concentration
  # gradient
  lake.P <- array(0, dim = c(xdim, ydim, timesteps))
  
  if(gradient){
    lake.P[,,1] <- matrix(sort(runif(xdim*ydim, pmin, pmax)), nrow = ydim, ncol = xdim)
    lake.P[,,1] <- lake.P[sample(1:nrow(lake.P)),,1]
  }else{
    lake.P[,,1] <- matrix((runif(xdim*ydim, pmin, pmax)), nrow = ydim, ncol = xdim)
  }
  
  # initial biomass for both species
  lake.x[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
  lake.y[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
  #lake.z[,,1] <- matrix(runif(ydim * xdim), ydim, xdim)
  lake.z[,,1] <- matrix((lake.P[,,1] - theta * lake.y[,,1])/lake.x[,,1])
  
  return(list(lake.x, lake.y, lake.z, lake.P))
}


## Simulation

run_dynamics <- function(lakes, pars, mig.method = c("random", "density", "quality", "quota"), p.mig = .1){
  
  if(length(mig.method) > 1){
    mig.method <- mig.method[1]
    warning("Only select one migration method. Using first value.")
  }
  
  lake.x <- lakes[[1]] 
  lake.y <- lakes[[2]]
  lake.z <- lakes[[3]]
  lake.P <- lakes[[4]]
  
  timesteps <- dim(lake.x)[3]
  
  t0 <- Sys.time()
  for(t in 1:(timesteps - 1)){
    
    for(i in 1:nrow(lake.x)){
      for(j in 1:ncol(lake.x)){
        lake.x[,,t][lake.x[,,t] <= 10^-10] <- 0
        lake.y[,,t][lake.y[,,t] <= 10^-10] <- 0
        lake.x[,,t][lake.x[,,t] == 0] <- 10^-4
        
        
        pars$P <- lake.P[i,j,t]
        state = c(x = lake.x[i,j,t], y = lake.y[i,j,t], Q = lake.z[i,j,t])
        #state = c(x = lake.x[i,j,t], y = lake.y[i,j,t])
        
        #out <- ode(state, 1:2, P14_fixed_alt, pars)
        out <- ode(state, 1:2, P14_fixed, pars)
        
        
        lake.x[i,j,t+1] <- ifelse(is.nan(out[2, 2]), 10^-6, out[2,2])
        lake.y[i,j,t+1] <- ifelse(is.nan(out[2, 3]), 0, out[2,3])
        lake.z[i,j,t+1] <- ifelse(is.nan(out[2, 4]), 0, out[2,4])
        
      }
      
    }
    
    
    lake.x[,,t+1][lake.x[,,t+1] < 0] <- 0
    lake.y[,,t+1][lake.y[,,t+1] < 0] <- 0
    
    ## Zoop Migration
    mig <- migration_diffuse(lake.y[,,t+1], lake.x[,,t+1], lake.z[,,t+1], lake.P[,,t], 
                     pars$theta,
                     max.move = 1, prop.migrant = p.mig, method = mig.method)
    lake.y[,,t+1] <- mig[[1]]
    lake.P[,,t+1] <- mig[[2]]
    
    if(((t + 1) %% 10) == 0){cat("Day ", t+1, "of ", timesteps, "simulated\n")}
  }
  t1 <- Sys.time()
  cat("Elapsed Time: ", t1-t0, "\n")
  
  return(list(lake.x, lake.y, lake.z, lake.P))
}


#filenames <- c("x100_m-8_p-03.rds","x100_m-8_p-08.rds","x100_m-8_p-13.rds","x100_m-8_p-18.rds")
filenames <- c("diff_x100_m-8_p-03.rds","diff_x100_m-8_p-08.rds","diff_x100_m-8_p-13.rds","diff_x100_m-8_p-18.rds")
#filenames <- c("x100_m-8_p-03_alt.rds","x100_m-8_p-18_alt.rds")
totP <- c(.03, .08, .13, .18)
#totP <- c(.03, .18)
timing <-c()
for(x in 1:2){
  ## Example 
  tA <- Sys.time()
  mlak <- make_lakes(100, 100, 300, totP[x], totP[x], par.P14fixed$theta)
  mig.types <- c("random", "density", "quality", "quota")
  dyn <- list()
  for(i in 1:4){
    dyn[[i]] <- run_dynamics(mlak, par.P14fixed, mig.method = mig.types[i], p.mig = .8)
  }
  
  saveRDS(dyn, filenames[x])
  rm(dyn)
  tB <- Sys.time()
  
  timing[x] <- tB - tA
}


dyn <- readRDS("x100_m-8_p-03.rds")

par(mfcol = c(4,4), mar = c(1,2,1,.2))
for(i in 1:4){
  matplot(t(apply(dyn[[i]][[4]], 3, as.vector))[,1:10], typ = "l")
  matplot(t(apply(dyn[[i]][[3]], 3, as.vector))[,1:10], typ = "l")
  matplot(t(apply(dyn[[i]][[2]], 3, as.vector))[,1:10], typ = "l")
  matplot(t(apply(dyn[[i]][[1]], 3, as.vector))[,1:10], typ = "l")
}

library(ggplot2)
library(animation)
library(gganimate)
library(patchwork)

t1 <- 93

dfphyto <- (reshape2::melt(dyn[[2]][[1]][,,t1]))
dfphyto$spp <- "Phytoplankton"
pp <- ggplot(dfphyto, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradient(limits = c(0, 1.5), low = "lightblue", high = "darkgreen") + theme_void() + facet_wrap(~spp)


dfzoop <- (reshape2::melt(dyn[[2]][[2]][,,t1]))
dfzoop$spp <- "Zooplankton"
zp <- ggplot(dfzoop, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_gradient(limits = c(0, 6), low = "lightblue", high = "darkblue") + theme_void() + facet_wrap(~spp)

pp + zp + plot_layout(ncol = 1, heights = c(1, 1))

boxplot(dyn[[2]][[2]])
max(dfzoop$value)

dyn <- readRDS("x100_m-8_p-03.rds")

dfphyto <- (reshape2::melt(dyn[[1]][[1]]))
anim <- ggplot(dfphyto, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
  scale_fill_continuous(low = "lightblue", high = "darkgreen") + 
  labs(title = "Phytoplankton in timestep: {closest_state}", fill = "Biomass") + 
  theme_void() + 
  transition_states(Var3, transition_length = 1, state_length = 1) + 
  enter_fade() + exit_fade()

animate(anim, fps = 20, duration = 30)
