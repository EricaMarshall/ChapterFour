Developing_targeted <- function(baseRaster, weightRaster=NULL, devArea=100, repeats=1,
                                 bufZone=1500, maxvalue=FALSE, timestep = FALSE, 
                                 Threshold = 0,  mainDir = "Z:/Erica/CHAPTER_3",
                                 subDir = "outputDirectory", SDM = NULL,
                                 Scenario){
  require(dplyr)
  require(dismo)
  require(spdep)
  require(raster)
  require(progress)
  if(is.null(weightRaster)){
    newRaster <- baseRaster
  } else{
    newRaster <- baseRaster * weightRaster
  }
  fullRaster <- fullRaster2 <- newRaster
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    
  }
  newRaster <- calc(newRaster, function(x){x[x < Threshold] <- NA; return(x)})
  # newRaster <- na.omit(newRaster)
  message("Pre-processing is done!")
  dd2 <- c()
  ids <- c() # save IDs of the dev sites
  for(t in 1:repeats){
    npixle <- 0
    while(npixle < devArea * 1.1){ # find buffers  with more than the minimum number of cells needed
      oldw <- getOption("warn")
      options(warn = -1)
      firstrand <- dismo::randomPoints(newRaster, 1)
      buf <- raster::buffer(SpatialPoints(firstrand), bufZone) # buffer size!
      samllRaster <- newRaster %>% 
        crop(buf) %>% 
        mask(buf)
      point2 <- rasterToPoints(samllRaster, spatial=TRUE, na.rm = TRUE)
      names(point2)[1] <- "suitability"
      if(!is.null(weightRaster)){
        point2$base <- raster::extract(baseRaster, point2)
      }
      npixle <- nrow(point2)
      options(warn = oldw)
    }
    par(mfrow=c(1,2))
    plot(newRaster, zlim=c(0,1), legend.width=1.2)
    plot(buf, add=T)
    plot(samllRaster, zlim=c(0,1), legend.width=1.2)
    neighbour <- spdep::knearneigh(point2, k=8)
    mp <- as.data.frame(firstrand)
    indices <- c()
    indices[1] <- which(neighbour$x[,1] == mp[1,1] & neighbour$x[,2] == mp[1,2])
    pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                     total=devArea-1, clear=FALSE, width=75) # add progress bar
    # browser()
    for(j in 2:devArea){
      # browser()
      nn <- neighbour$nn[indices,]
      nn <- nn[!nn %in% indices] # remove the previous indices from the neighbours
      if(maxvalue==TRUE){
        npoints <- point2[nn, ]
        npmax <- which(npoints@data[,"suitability"] == max(npoints$suitability))
        mp3 <- as.data.frame(npoints[npmax,]@coords)
        randomnn <- which(neighbour$x[,1] == mp3[1,1] & neighbour$x[,2] == mp3[1,2])
        indices[j] <- randomnn
        thePoint <- neighbour$x[randomnn,]
        points(thePoint[1], thePoint[2], pch=16, cex=0.4)
      } else{
        if(length(nn) > 1){ # when lenght(nn) < 1 it consider it 1:nn
          randomnn <- sample(nn, size=1)#, prob=point2$suitability[nn])
          indices[j] <- randomnn
          thePoint <- neighbour$x[randomnn,]
          points(thePoint[1], thePoint[2], pch=16, cex=0.4)
        } else if(length(nn) == 1){
          randomnn <- nn
          indices[j] <- randomnn
          thePoint <- neighbour$x[randomnn,]
          points(thePoint[1], thePoint[2], pch=16, cex=0.4)
        } else{
          message("A new patch is being search for neighbourhood!")
          searchSpace <- 1:nrow(neighbour$nn)
          searchSpace <- searchSpace[!searchSpace %in% indices]
          randNew <- sample(searchSpace, 1)
          indices[j] <- randNew
          thePoint <- neighbour$x[randNew,]
          points(thePoint[1], thePoint[2], pch=16, cex=0.4)
          
        }
      }
      pb$tick()
    }
    message("Some development happened! :/")
    if(length(devArea) > 10 || !is.null(devArea)){
      dd <- c()
      dd <- raster::cellFromXY(newRaster, neighbour$x[indices, ])
      dd2 <- append(dd2, dd)
      # newdd <- dd2
      ids <- append(ids, rep(t+1, length(dd)))
      if(repeats > 1 && t < repeats){
        raster::values(newRaster)[dd] <- NA
      }
      if(timestep == TRUE){
        raster::values(fullRaster)[dd2] <- NA
        Dev_t <- fullRaster
        writeRaster(Dev_t, paste0(Scenario,"Dev_", t, ".asc"), format = "ascii", overwrite = TRUE)
      }
    }
    print(paste("The number of development cells removed:", length(dd2)))
    message(paste("Repeat", t, "is done!", "Hooray! :))"))
  }
  if(anyNA(dd2)){
    print(paste("Nummber of NAs:", length(which(is.na(dd2)))))
    warning("The cells have NAs")
    ids <- ids[- which(is.na(dd2))]
    dd2 <- dd2[- which(is.na(dd2))]
  }
  # print(dd2)
  # View(dd2)
  # browser()
  raster::values(fullRaster)[dd2] <- NA
  raster::values(fullRaster2) <- NA
  raster::values(fullRaster2)[dd2] <- ids
  names(fullRaster) <- "Development"
  names(fullRaster2) <- "IDs"
  return(stack(fullRaster, fullRaster2)) 
}

#debugonce(Developing_Threshold)
