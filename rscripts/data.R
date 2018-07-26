# Data
library(data.table)


# compile 3d array into 2d dataframe
collect_bio <- function(lake.x, lake.y){
  
  xlake <- list()
  ylake <- list()
  for(i in 1:xdim){
    for(j in 1:ydim){
      xlake[[i]] <- data.frame(lake = "x", col = i, row = j, time = 1:length(lake.x[j,i,]),
                               biomass = lake.x[j,i,])
      ylake[[i]] <- data.frame(lake = "y", col = i, row = j, time = 1:length(lake.x[j,i,]),
                               biomass = lake.y[j,i,])
    }
  }
  

  lakebio <- rbind(rbindlist(xlake), rbindlist(ylake))
  return(lakebio)
}



lakebio2 <- collect_bio(lake.x, lake.y)