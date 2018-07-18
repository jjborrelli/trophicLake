# Migration

## Completely randomized migration

migration.random <- function(lake, max.move = 1){
  
  new.lake <- lake
  for(i in 1:nrow(lake)){
    for(j in 1:ncol(lake)){
      migrants <- lake[i,j] * runif(1,0,.1)
      new.lake[i,j] <- lake[i,j] - migrants
      
      xdir <- seq(j-max.move, j+max.move, 1)
      xdir <- xdir[xdir %in% 1:ncol(lake)]
      
      ydir <- seq(i-max.move, i+max.move, 1)
      ydir <- ydir[ydir %in% 1:nrow(lake)]
      
      
      newcell <- c(sample(xdir, 1), sample(ydir, 1))
      
      new.lake[newcell[2], newcell[1]] <- new.lake[newcell[2], newcell[1]] + migrants
      
    }
  }
  return(new.lake)
}


## Migration based on probability derived by prey density

migration.density <- function(ylake, xlake, max.move = 1, prop.migrant){
  new.lake <- ylake
  for(i in 1:nrow(ylake)){
    for(j in 1:ncol(ylake)){
      #migrants <- ylake[i,j] * prop.migrant
      #new.lake[i,j] <- ylake[i,j] - migrants
      
      xdir <- seq(j-max.move, j+max.move, 1)
      xdir <- xdir[xdir %in% 1:ncol(xlake)]
      
      ydir <- seq(i-max.move, i+max.move, 1)
      ydir <- ydir[ydir %in% 1:nrow(xlake)]
      
      high <- reshape2::melt((xlake[ydir, xdir]))
      high$Var1 <- xdir[high$Var1]
      high$Var2 <- ydir[high$Var2]
      high$prob <- order(high$value)/sum(order(high$value))
      
      #newcell <- high[sample(1:nrow(high), 1, prob = high$prob),]
      
      for(x in 1:nrow(high)){
        new.lake[high$Var2[x], high$Var1[x]] <- new.lake[high$Var2[x], high$Var1[x]] + 
          ylake[i,j] * high$prob[x]
      }
      new.lake[i,j] <- new.lake[i,j] - ylake[i,j]
      
      
    }
  }
  return(new.lake)
}


migration.diffusion <- function(ylake, xlake, max.move = 1){
  new.lake <- ylake
  for(i in 1:nrow(ylake)){
    for(j in 1:ncol(ylake)){
      #migrants <- ylake[i,j] * prop.migrant
      #new.lake[i,j] <- ylake[i,j] - migrants
      
      xdir <- seq(i-max.move, i+max.move, 1)
      xdir <- xdir[xdir %in% 1:ncol(xlake)]
      
      ydir <- seq(j-max.move, j+max.move, 1)
      ydir <- ydir[ydir %in% 1:nrow(xlake)]
      
      high <- reshape2::melt(xlake[ydir, xdir])
      high$Var1 <- xdir[high$Var1]
      high$Var2 <- ydir[high$Var2]
      high$prob <- order(high$value)/sum(order(high$value))
      
      #newcell <- high[sample(1:nrow(high), 1, prob = high$prob),]
      
      for(x in 1:nrow(high)){
        new.lake[high$Var1[x], high$Var2[x]] <- new.lake[high$Var1[x], high$Var2[x]] + 
          ylake[i,j] * high$prob[x]
      }
      new.lake[i,j] <- new.lake[i,j] - ylake[i,j]
      
      
    }
  }
  return(new.lake)
}


migration.density1 <- function(ylake, xlake, max.move = 1, prop.migrant){
  new.lake <- ylake
  for(i in 1:nrow(ylake)){
    for(j in 1:ncol(ylake)){
      migrants <- ylake[i,j] * prop.migrant
      new.lake[i,j] <- ylake[i,j] - migrants
      
      xdir <- seq(i-max.move, i+max.move, 1)
      xdir <- xdir[xdir %in% 1:ncol(xlake)]
      
      ydir <- seq(j-max.move, j+max.move, 1)
      ydir <- ydir[ydir %in% 1:nrow(xlake)]
      
      high <- reshape2::melt(xlake[xdir, ydir])
      high$Var1 <- xdir[high$Var1]
      high$Var2 <- ydir[high$Var2]
      high$prob <- order(high$value)/sum(order(high$value))
      
      newcell <- high[sample(1:nrow(high), 1, prob = high$prob),]
      
      new.lake[newcell$Var1, newcell$Var2] <- new.lake[newcell$Var1, newcell$Var2] + migrants
      
      
    }
  }
  return(new.lake)
}
