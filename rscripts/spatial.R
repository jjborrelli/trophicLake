library(gstat)
library(sp)
library(reshape2)
library(dplyr)

phyto <- apply(dyn[[1]][[1]], 3, melt)

for(i in 1:length(phyto)){
  coordinates(phyto[[i]]) <- ~Var1+Var2
}


t1 <- Sys.time()
vgms <- lapply(phyto, function(x) variogram(log(value) ~ 1, data = x))
t2 <- Sys.time()
t2 - t1
vgms[[2]]
plot(vgms[[1]])

x = 300
plot(vgms[[x]]$dist, vgms[[x]]$gamma, typ = "o")

tA <- Sys.time()
fv <- list()
fvr <- c()
for(i in 1:length(vgms)){
  fv[[i]] <- fit.variogram(vgms[[i]], vgm(c("Lin")))
  fvr[i] <- fv[[i]]$range[2]
}
tB <- Sys.time()
tB - tA

plot(vgms[[100]], fv[[100]])
plot(fvr)


library(lctools)
phyto <- apply(dyn[[1]][[1]], 3, melt)

tm1 <- Sys.time()
mIlist <- list()
times <- seq(25, 300, 25)
for(i in 1:length(times)){
  mIlist[[i]] <- moransI(phyto[[times[i]]][,1:2],6,phyto[[times[i]]][,3])[-1]
  print(i)
}
tm2 <- Sys.time()
tm2 - tm1


par(mfrow = c(4,4))
pheatmap::pheatmap(dyn[[1]][[1]][,,200], cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(dyn[[2]][[1]][,,200], cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(dyn[[3]][[1]][,,200], cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(dyn[[4]][[1]][,,200], cluster_rows = F, cluster_cols = F)


x <- cbind(1:98, 2:99, 3:100)
allmd <- list()
for(i in 1:nrow(x)){
  md <- c()
  for(j in 1:nrow(x)){
    md[j] <- mean(dist(dyn[[1]][[1]][x[i,],x[j,],22]))
  }
  allmd[[i]] <- md
}
mean(unlist(allmd))
