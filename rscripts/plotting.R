# Plotting
library(animation)

## plot the nth column of each "lake" through time

threePlot <- function(n, lake1, lake2, lake3){
  par(mfrow = c(3,1), mar = c(1,2,1,.2))
  matplot(t(lake1[,n,]), typ = "l")
  matplot(t(lake2[,n,]), typ = "l")
  matplot(t(lake3[,n,]), typ = "l")
}


## plot heatmap of the dynamics

spatdynplot <- function(lake, t, interval, name){
  
  saveGIF({
    for(i in 1:t){
      (pheatmap(lake[,,i], cluster_rows = F, cluster_cols = F, main = i))
    }
  }, interval = interval, movie.name = name)
  
}
