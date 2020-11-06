# Authors: Erica Marshall and Roozbeh Valavi
# contact: marshalle@student.unimelb.edu.au / valavi.r@gmail.com
# Date : December 2019
# Version 0.3
# incrementCell: the number of cells that increase around the 
# previous buffer each time - the higher the faster
# but too much high might make a patch disconnected
offsetting_timestep <- function(baseRaster, weightRaster, devRaster, bufferZone=10000, 
                                mainDir = "Z:/Erica/CHAPTER_3",
                                subDir = "outputDirectory", offPatches=10, 
                                type="value", threshold=0.5,
                                offsetMultiplier=1, incrementCell=3){
  require(dplyr)
  require(sf)
  require(dismo)
  require(fasterize)
  newRaster <- baseRaster * weightRaster
  fullRaster <- fullRaster2 <- theOriginalRaster <- newRaster
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
  }
  # Convert each unique development site into its own polygon so that developments may be offset individually
  devRaster[values(devRaster)== 0] <- NA
  message(paste("The number of development cells is", length(which(!is.na(values(devRaster))))))
  # Convert each unique development site into its own polygon so that developments may be offset individually
  dev <- raster::rasterToPolygons(devRaster, dissolve = TRUE) %>%
    sf::st_as_sf()
  # sf::st_cast("MULTIPOLYGON")
  dev$area <- st_area(dev)
  # remove small polygons
  dev <- dev[which(dev$area > prod(raster::res(devRaster)) * 20), ] %>% 
    sf::as_Spatial()
  # make a buffer to target offsetting starting point
  # devBuff <- buffer(dev, bufferZone) %>%
  #   crop(newRaster)
  devValue <- sum(values(newRaster)[which(!is.na(values(devRaster)))]) # total dev values
  if(type == "area"){
    offsetvalueTotal <- developedArea * offsetMultiplier
  } else if(type == "value"){
    offsetvalueTotal <- devValue * offsetMultiplier
  }
  values(newRaster)[which(!is.na(values(devRaster)))] <- NA
  values(baseRaster)[which(!is.na(values(devRaster)))] <- NA
  if(type == "area"){
    restorationPotential <- length(newRaster[values(newRaster) < threshold])
  }else if (type == "value"){
    restorationPotential <- sum(baseRaster[], na.rm = TRUE)-sum(newRaster[], na.rm = TRUE)
  }
  if(offsetvalueTotal > restorationPotential){
    message("The offset requirements are greater than the restoration potential of the landscape and this offset will fail")
  }
  par(mfrow=c(1,2))
  plot(newRaster)
  plot(dev, add=TRUE)
  # plot(devBuff, add=T, border="red")
  plot(newRaster)
  # plot(devBuff, add=T, border="red")
  # devBuff <- devBuff - dev
  # buffRaster <- fasterize(sf::st_as_sf(devBuff), devRaster)
  # plot(buffRaster)
  oo <- c() # save all the offset indices
  gg <- c() # save all the offset gain (values)
  fullRaster <- newRaster
  if(sum(length(dev)) < 2){
    x <- offPatches
  } else{
    x <- seq_len(nrow(dev))
    message("The number of habitat patches: ", nrow(dev))
  }
  randpoint <- rasterToPoints(newRaster, spatial = TRUE)
  message("Pre-processing is done!")
  for(t in x){
    # the desired offset values for the patch
    if(sum(x) < 2){
      offsetvalue <- offsetvalueTotal / offPatches
    } else{
      offsetvalue <- sum(unlist(raster::extract(theOriginalRaster, dev[t,])), na.rm = TRUE) * offsetMultiplier
    }
    if(length(which(!is.na(values(newRaster)))) < 10){
      message(cat("There is no cell left in the landscape", "\n",
                  "Of habitat was not offset", offsetvalueTotal - offsetvalue))
      dtf <- data.frame(OffsetRequired = offsetvalueTotal, 
                        ValueAchieved = offsetvalue)
      break
    }
    if(sum(which(!is.na(values(newRaster)))) < offsetvalue){
      message(cat("There is no cell left in the landscape", "\n",
                  "Of habitat was not offset", offsetvalueTotal - offsetvalue))
      break
    }
    of_indices <- c()
    highPoints <- c()
    offsetGain <- c() # to count the offset gain we need to add to the raster later
    e <- 0
    s <- 0
    n <- 0
    rfp <- NA
    rfp2 <- NA
    # browser()
    while(is.na(rfp) || rfp2 < threshold){
      # nn <- dismo::randomPoints(buffRaster, 1)
      nn <- randpoint[sample(1:nrow(randpoint), 1),]
      rfp <- values(newRaster)[cellFromXY(newRaster, nn)]
      rfp2 <- values(baseRaster)[cellFromXY(baseRaster, nn)]
    }
    points(nn)
    # calculate minimum buffer
    if(type == "area"){
      resRaster <- prod(res(newRaster))
      rr <- sqrt((resRaster * offsetvalue) / 3.141593)
    } else if(type == "value"){
      rr <- offsetvalue 
    }
    while(s < offsetvalue){
      b <- buffer(SpatialPoints(nn), rr + e)
      # plot(b, add=TRUE)
      e <- e + (res(newRaster)[1] * incrementCell)
      v <- unlist(cellFromPolygon(newRaster, b))
      v <- v[which(!is.na(values(newRaster)[v]))]
      if(length(v) > 0){
        for(i in v){
          if(is.null(weightRaster)){
            n <- n + 1   # Count it as an offset point
            of_indices[n] <- i # store it as an a already considered indice
            offsetGain[n] <- threshold
            if(type == "value"){
              s <- s + (threshold - values(newRaster)[i]) # condition value
            } else if(type == "area"){
              s <- n
            } 
          } else{
            if(values(baseRaster)[i] >= threshold){ # value of suitability or combined raster
              if(values(newRaster)[i] <= values(baseRaster)[i]){
                n <- n + 1
                of_indices[n] <- i
                offsetGain[n] <- values(baseRaster)[i]
                if(type == "value"){
                  s <- s + (values(baseRaster)[i] - values(newRaster)[i])
                } else if(type == "area"){
                  s <- n
                } 
              } else{
                highPoints <- append(highPoints, i)
              }
            } else{
              highPoints <- append(highPoints, i)
            }
          }
          if(s >= offsetvalue){
            break
          } 
        }
      }
      searchedPoints <- unique(append(of_indices, highPoints))
      points(xyFromCell(newRaster, of_indices), col="blue", cex=0.1)
      values(newRaster)[searchedPoints] <- NA
      values(baseRaster)[searchedPoints] <- NA
      # values(buffRaster)[searchedPoints] <- NA
      if(length(which(!is.na(values(newRaster)))) < 10){
        message(paste("There is no more available habitat for restoration"))
        break
      }
      # if(all(which(!is.na(values(newRaster)))) >= all(which(!is.na(values(baseRaster))))){
      #   message(print("There is no more available habitat for restoration"))
      # break
      # }
      if(type == "area"){
        progress <- length(which(!is.na(newRaster)))
      }else if (type == "value"){
        progress <- sum(which(!is.na(newRaster[])))
      }
      if(progress == restorationPotential){
        break
        message("There is no more available habitat for restoration")
      }
    }
    message("Biodiversity offsetting recovered the habitat quality! :)")
    if(length(of_indices) > 1 || !is.null(of_indices)){
      gg <- append(gg, offsetGain)
      oo <- append(oo, of_indices)
      # save the raster of each patch
      values(theOriginalRaster)[unlist(cellFromPolygon(theOriginalRaster, dev[t,]))] <- NA
      theOriginalRaster[of_indices] <- offsetGain
      writeRaster(theOriginalRaster, paste0("offset_patch_", t, ".tif"), overwrite = TRUE)
      message(paste("Repeat", t, "is done!", "Hooray! :))"))
    } else{
      values(theOriginalRaster)[unlist(cellFromPolygon(theOriginalRaster, dev[t,]))] <- NA
      writeRaster(theOriginalRaster, paste0("offset_patch_", t, ".tif"), overwrite = TRUE)
      message(paste("Repeat", t, "is done!", "Hooray! :))"))
      message(paste("The offset of repeat", t, "has not been taken"))
      message(paste("The offset value is", offsetvalue))
    }
  }
  print(paste("The number of cells used for offsetting:", length(oo)))
  print(paste("Sum of the development value previously taken:", devValue))
  print(paste("Sum of the offset gain added:", sum(gg) - sum(values(fullRaster)[oo])))
  fullRaster[oo] <- gg
  return(fullRaster)
}

