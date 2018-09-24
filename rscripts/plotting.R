# Plotting
library(animation)
library(pheatmap)

## plot the nth column of each "lake" through time

threePlot <- function(n, lake1, lake2, lake3){
  par(mfrow = c(3,1), mar = c(1,2,1,.2))
  matplot(t(lake1[,n,]), typ = "l")
  matplot(t(lake2[,n,]), typ = "l")
  matplot(t(lake3[,n,]), typ = "l")
}

lakeDynPlot <- function(dyn){
  par(mfrow = c(4,1), mar = c(1,2,1,.2))
  matplot(t(apply(dyn[[4]], 3, as.vector))[-c(1:5),], typ = "l")
  matplot(t(apply(dyn[[3]], 3, as.vector))[-c(1:5),], typ = "l")
  matplot(t(apply(dyn[[2]], 3, as.vector))[-c(1:5),], typ = "l")
  matplot(t(apply(dyn[[1]], 3, as.vector))[-c(1:5),], typ = "l")
} 


## plot heatmap of the dynamics

spatdynplot <- function(lake, t, interval, name){
  
  saveHTML({
    for(i in 1:t){
      (pheatmap(lake[,,i], cluster_rows = F, cluster_cols = F, main = i))
    }
  }, interval = interval, movie.name = name)
  
}
