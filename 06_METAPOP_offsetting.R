offsetting_Metapop <- function(baseRaster,newRaster, weightRaster, devRaster, 
                                 mainDir = "Z:/Erica/CHAPTER_3",
                                 subDir = "outputDirectory", offPatches=10, 
                                 type="value",threshold=0.5,SDM,bufferZone =5000,
                                 offsetMultiplier=1, incrementCell=3){
  require(dplyr)
  require(sf)
  require(dismo)
  require(fasterize)
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    
  }
  theOriginalRaster <- theOriginalRaster2 <-  newRaster
  differenceRaster <- baseRaster - newRaster
  fullRaster <- fullRaster2 <- SDM*weightRaster
  # Convert each unique development site into its 
  # count the number of dev pixels
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
  devBuff <- buffer(dev, bufferZone) %>%
    crop(newRaster)
  devValue <- sum(values(newRaster)[which(!is.na(values(devRaster)))]) # total dev values
  values(differenceRaster)[which(!is.na(values(devRaster)))] <- NA
  values(fullRaster)[which(!is.na(values(devRaster)))] <- NA
  values(theOriginalRaster2)[which(!is.na(values(devRaster)))] <- NA
  par(mfrow=c(1,2))
  plot(differenceRaster)
  plot(dev, add=TRUE)
  plot(devBuff, add=T, border="red")
  plot(differenceRaster)
  plot(devBuff, add=T, border="red")
  devBuff <- devBuff - dev
  buffRaster <- fasterize(sf::st_as_sf(devBuff), devRaster)
  # plot(buffRaster)
  if(type == "area"){
    offsetvalueTotal <- developedArea * offsetMultiplier
  } else if(type == "value"){
    offsetvalueTotal <- devValue * offsetMultiplier
  }
  oo <- c() # save all the offset indices
  gg <- c() # save all the offset gain (values)
  hh <- c()
  Metapop <- c()
  Repeat <- c()
  if(sum(length(dev)) < 2){
    x <- offPatches
  } else{
    x <- seq_len(nrow(dev))
    message("The number of habitat patches: ", nrow(dev))
  }
  message("Pre-processing is done!")
  for(t in x){
    # the desired offset values for the patch
    if(sum(x) < 2){
      offsetvalue <- offsetvalueTotal / offPatches
    } else{
      offsetvalue <- sum(unlist(raster::extract(theOriginalRaster, dev[t,])), na.rm = TRUE) * offsetMultiplier
    }
    if(length(which(!is.na(values(newRaster)))) < 1){
      message(cat("There is no cell left in the landscape", "\n",
                  "Of habitat was not offset", offsetvalueTotal - offsetvalue))
      break
    }
    of_indices <- c()
    highPoints <- c()
    offsetGain <- c() # to count the offset gain we need to add to the raster later
    HSGain <- c()
    e <- 0
    s <- 0
    n <- 0
    rfp <- NA
    rfp2 <- NA
    # browser()
    while(is.na(rfp) || rfp2 < threshold){
      nn <- dismo::randomPoints(buffRaster, 1)
      rfp <- values(newRaster)[cellFromXY(newRaster, nn)]
      rfp2 <- values(SDM)[cellFromXY(SDM, nn)]
    }
    points(nn)
    # calculate minimum buffer
    if(type == "area"){
      resRaster <- prod(res(newRaster))
      rr <- sqrt((resRaster * offsetvalue) / 3.141593)
    } else if(type == "value"){
      rr <- offsetvalue / 2
    }
    while(s < offsetvalue){
      if(offsetvalue == 0){
        break
      } else {
        b <- buffer(SpatialPoints(nn), rr + e)
        # plot(b, add=TRUE)
        e <- e + (res(differenceRaster)[1] * incrementCell)
        v <- unlist(cellFromPolygon(differenceRaster, b))
        v <- v[which(!is.na(values(differenceRaster)[v]))]
        # if(length(v) > 0){
        for (j in v) {
          n <- n + 1
          of_indices[n] <- j
          if(type == "value"){
            s <- s + values(differenceRaster)[j]
          } else if(type == "area"){
            s <- n
          }
          if(s >= offsetvalue){
            break
          } 
        }
        searchedPoints <- unique(append(of_indices, highPoints))
        points(xyFromCell(differenceRaster, of_indices), col="blue", cex=0.1)
        values(buffRaster)[searchedPoints] <- NA
        values(differenceRaster)[searchedPoints] <- NA
        if(length(which(!is.na(values(newRaster)))) < 10){
          break
        }
        if(length(which(!is.na(values(differenceRaster)))) < 1){
          a <- offsetvalue-s
          Metapop <- append(Metapop, a)
          Repeat <- append(Repeat, t)
          break
        }
      }
    }
    message("Biodiversity offsetting recovered the habitat quality! :)")
    # offsetGain[of_indices] <- values(baseRaster)[of_indices]
    # HSGain[of_indices] <- values(SDM)[of_indices]
    if(length(of_indices) > 10 || !is.null(of_indices)){
      # hh <- append(hh, HSGain)
      # gg <- append(gg, offsetGain)
      oo <- append(oo, of_indices)
      values(theOriginalRaster)[unlist(cellFromPolygon(theOriginalRaster, dev[t,]))] <- NA
      theOriginalRaster[of_indices] <- offsetGain[of_indices]
      values(fullRaster2)[unlist(cellFromPolygon(fullRaster2, dev[t,]))] <- NA
      fullRaster2[of_indices] <- HSGain[of_indices]
      writeRaster(theOriginalRaster, paste0("Metapop_patch_", t, ".tif"), overwrite = TRUE)
      writeRaster(fullRaster2, paste0("HS_Patch_", t, ".tif"), overwrite = TRUE)
      message(paste("Repeat", t, "is done!", "Hooray! :))"))
    } else{
      message(paste("The offset of repeat", t, "has not been taken"))
      message(paste("The offset value is", offsetvalue))
      values(theOriginalRaster)[unlist(cellFromPolygon(theOriginalRaster, dev[t,]))] <- NA
      values(fullRaster2)[unlist(cellFromPolygon(fullRaster2, dev[t,]))] <- NA
      writeRaster(theOriginalRaster, paste0("Metapop_patch_", t, ".tif"), overwrite = TRUE)
      writeRaster(fullRaster2, paste0("HS_Patch_", t, ".tif"), overwrite = TRUE)
      message(paste("Repeat", t, "is done!", "Hooray! :))"))
    }
    # save the raster of each patch
    
  }
  print(paste("The number of cells used for offsetting:", length(oo)))
  print(paste("Sum of the development value previously taken:", devValue))
  print(paste("Sum of the offset gain added:", sum(values(differenceRaster)[oo])))
  dtf <- data.frame(Run = Repeat,
                    Unoffset = Metapop)
  saveRDS(dtf, file ="Metaopop_notoffset.rds")
  fullRaster[oo] <- values(SDM)[oo]
  theOriginalRaster2[oo] <- values(baseRaster)[oo]
  my_rasters <- list(fullRaster, theOriginalRaster2)
  return(my_rasters)
}
