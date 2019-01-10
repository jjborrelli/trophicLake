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


migration <- function(ylake, xlake, qlake, plake, theta, prop.migrant = 0.1, method = "random", max.move = 1){
  
  new.lake.P <- plake
  # create new lake to reflect migration
  new.lake <- ylake
  # how much migration from each cell
  lake.migrant <- ylake * prop.migrant
  
  for(i in 1:nrow(new.lake)){
    for(j in 1:ncol(new.lake)){
      
      # determine all possible vertical moves
      rdir <- seq(i-max.move, i+max.move, 1)
      rdir <- rdir[rdir %in% 1:nrow(xlake)]
      #determine all possible horizontal moves
      cdir <- seq(j-max.move, j+max.move, 1)
      cdir <- cdir[cdir %in% 1:ncol(xlake)]
      
      if(method == "random"){
        xsampl <- sample(rdir, 1)
        ysampl <- sample(cdir, 1)
      }
      
      if(method == "density"){
        high <- reshape2::melt(xlake[rdir, cdir])
        high$Var1 <- rdir[high$Var1]
        high$Var2 <- cdir[high$Var2]
        
        if(sum(high$value > 0)){
          high$prob <- (high$value)/sum((high$value))
        }else{
          high$prob <- .25
        }
        
        newcell <- high[sample(1:nrow(high), 1, prob = high$prob),]
        xsampl <- newcell$Var1
        ysampl <- newcell$Var2
      }
      
      if(method == "quality"){
        high <- reshape2::melt(xlake[rdir, cdir] * qlake[rdir, cdir])
        high$Var1 <- rdir[high$Var1]
        high$Var2 <- cdir[high$Var2]
        if(sum(high$value > 0)){
          high$prob <- (high$value)/sum((high$value))
        }else{
          high$prob <- .25
        }
        newcell <- high[sample(1:nrow(high), 1, prob = high$prob),]
        xsampl <- newcell$Var1
        ysampl <- newcell$Var2
      }
      
      if(method == "quota"){
        high <- reshape2::melt(qlake[rdir, cdir])
        high$Var1 <- rdir[high$Var1]
        high$Var2 <- cdir[high$Var2]
        if(sum(high$value > 0)){
          high$prob <- (high$value)/sum((high$value))
        }else{
          high$prob <- .25
        }
        newcell <- high[sample(1:nrow(high), 1, prob = high$prob),]
        xsampl <- newcell$Var1
        ysampl <- newcell$Var2
      }
      
      pmoved <- lake.migrant * theta
      
      new.lake[xsampl, ysampl] <- new.lake[xsampl, ysampl] + lake.migrant[i,j]
      new.lake[i,j] <- new.lake[i,j] - lake.migrant[i,j]
      
      new.lake.P[xsampl, ysampl] <- new.lake.P[xsampl, ysampl] + pmoved[i,j]
      new.lake.P[i,j] <- new.lake.P[i,j] - pmoved[i,j]
    }
  }
  
  return(list(new.lake, new.lake.P))
}

migration_diffuse <- function(ylake, xlake, qlake, plake, theta, prop.migrant = 0.1, method = "random", max.move = 1){
  
  # how much migration from each cell
  lake.migrant <- ylake * prop.migrant

  # create new lake to reflect migration
  new.lake <- ylake - lake.migrant
  
  new.lake.P <- plake - (lake.migrant * theta)
  
  for(i in 1:nrow(new.lake)){
    for(j in 1:ncol(new.lake)){
      
      # determine all possible vertical moves
      rdir <- seq(i-max.move, i+max.move, 1)
      rdir <- rdir[rdir %in% 1:nrow(xlake)]
      #determine all possible horizontal moves
      cdir <- seq(j-max.move, j+max.move, 1)
      cdir <- cdir[cdir %in% 1:ncol(xlake)]
      
      if(method == "random"){
        cells <- expand.grid(rdir, cdir)
        zmove <- diff(c(0, sort(runif(nrow(cells)-1, 0, lake.migrant[i,j])), lake.migrant[i,j]))
        new.lake[rdir, cdir] <- new.lake[rdir, cdir] + zmove
        new.lake.P[rdir, cdir] <- new.lake.P[rdir, cdir] + (zmove * theta)
      }
      
      if(method == "density"){
        zmove <- (xlake[rdir,cdir]/sum(xlake[rdir,cdir]) * lake.migrant[i,j])
        new.lake[rdir, cdir] <- new.lake[rdir, cdir] + zmove
        new.lake.P[rdir, cdir] <- new.lake.P[rdir, cdir] + (zmove * theta)
      }
      
      if(method == "quality"){
        qualake <- xlake * qlake
        zmove <- (qualake[rdir,cdir]/sum(qualake[rdir,cdir]) * lake.migrant[i,j])
        new.lake[rdir, cdir] <- new.lake[rdir, cdir] + zmove
        new.lake.P[rdir, cdir] <- new.lake.P[rdir, cdir] + (zmove * theta)
      }
      
      if(method == "quota"){
        zmove <- (qlake[rdir,cdir]/sum(qlake[rdir,cdir]) * lake.migrant[i,j])
        new.lake[rdir, cdir] <- new.lake[rdir, cdir] + zmove
        new.lake.P[rdir, cdir] <- new.lake.P[rdir, cdir] + (zmove * theta)
      }
      
    }
  }
  
  return(list(new.lake, new.lake.P))
}
