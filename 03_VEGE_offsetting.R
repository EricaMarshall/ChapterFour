# Authors: Erica Marshall and Roozbeh Valavi
# contact: marshalle@student.unimelb.edu.au / valavi.r@gmail.com
# Date : December 2019
# Version 0.3
# incrementCell: the number of cells that increase around the 
# previous buffer each time - the higher the faster
# but too much high might make a patch disconnected
offsetting_Veg <- function(baseRaster, weightRaster, devRaster, vegRaster, 
                           mainDir = "Z:/Erica/CHAPTER_3",
                           subDir = "outputDirectory",  
                           type="value", threshold=0.5,
                           offsetMultiplier=1, incrementCell=3, csvdir){
  require(dplyr)
  require(ggplot2)
  require(sf)
  require(dismo)
  require(fasterize)
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
  # count the number of dev pixels
  theOriginalRaster <- newRaster # to retrive the dev values for each patch
  # Convert each unique development site into its own polygon so that developments may be offset individually
  devRaster[values(devRaster)== 0] <- NA
  message(paste("The number of development cells is", length(which(!is.na(values(devRaster))))))
  # Convert each unique development site into its own polygon so that developments may be offset individually
  
  dev <- raster::rasterToPolygons(devRaster, dissolve = TRUE) %>%
    sf::st_as_sf()
  dev$area <- st_area(dev)
  # remove small polygons
  dev <- dev[which(dev$area > prod(raster::res(devRaster)) * 20), ] %>% 
    sf::as_Spatial()
  devValue <- sum(values(newRaster)[which(!is.na(values(devRaster)))]) # total dev values
  if(type == "area"){
    offsetvalueTotal <- sum(dev$area) * offsetMultiplier
  } else if(type == "value"){
    offsetvalueTotal <- devValue * offsetMultiplier
  }
  if(type == "area"){
    restorationPotential <- length(newRaster[values(newRaster) < threshold])
  }else if (type == "value"){
    restorationPotential <- sum(values(baseRaster), na.rm = TRUE)-sum(values(newRaster), na.rm = TRUE)
  }
  if(offsetvalueTotal > restorationPotential){
    message("The offset requirements are greater than the restoration potential of the landscape and this offset will fail")
  }
  values(newRaster)[which(!is.na(values(devRaster)))] <- NA
  values(baseRaster)[which(!is.na(values(devRaster)))] <- NA
  plot(newRaster)
  plot(dev, add=TRUE)
  oo <- c() # save all the offset indices
  gg <- c() # save all the offset gain (values)
  message("The number of habitat patches: ", nrow(dev))
  message("Pre-processing is done!")
  # for ---------------------------------------------------------------------
  for(z in seq_len(nrow(dev))){
    devRaster <- fasterize(sf::st_as_sf(dev[z,]), devRaster)
    lostVeg <- vegRaster * devRaster
    lostCondition <- fullRaster * devRaster
    df <- as.data.frame(setNames(stack(lostVeg, lostCondition), c("class", "condition")), na.rm = TRUE)
    head(df)
    VegCount <- df %>% 
      group_by(class) %>% 
      summarise(sum(condition)) %>% 
      left_join(df %>% group_by(class) %>% count()) %>% 
      setNames(c("class", "condition", "count")) %>% 
      as.data.frame()
    ggbar <- ggplot(VegCount, aes(x = as.factor(class), y = count, fill = count)) +
      geom_bar(stat = "identity") +
      theme_bw()
    print(ggbar)
    values(newRaster)[which(!is.na(values(devRaster)))] <- NA
    values(baseRaster)[which(!is.na(values(devRaster)))] <- NA
    values(vegRaster)[which(!is.na(values(devRaster)))] <- NA
    if(length(which(!is.na(values(newRaster)))) <= 1){
      message(cat("There is no cell left in the landscape", "\n",
                  "Of habitat was not offset", offsetvalueTotal - offsetvalue))
      break
    }
    # message("Pre-processing is done!")
    oo <- c()
    gg <- c()
    # class loop --------------------------------------------------------------
    for(t in VegCount$class){
      
      if(type == "area"){
        offsetvalue <- VegCount$count[which(VegCount$class == t)]
      } else if(type == "value"){
        offsetvalue <- VegCount$condition[which(VegCount$class == t)]
      }
      offsetvalue <- offsetvalue * offsetMultiplier
      classVeg <- newRaster # the selected veg class raster
      values(classVeg)[which(values(vegRaster) != t)] <- NA
      classBase <- baseRaster
      values(classBase)[which(values(vegRaster) != t)] <- NA
      if(all(which(!is.na(values(classVeg)))) < threshold){
        message(paste("None of the cells in vege class", t, "are above the threshold for offsets and these impacts are unoffsetable"))
        break
      }
      if(all(which(!is.na(values(classBase)))) < threshold){
        message(paste("None of the cells in vege class", t, "are above the threshold for offsets and these impacts are unoffsetable"))
        break()
      }
      if(length(which(!is.na(values(classVeg)) <= 1))){
        message(cat("There are no cells left in the landscape from class ", t, "\n"))
        break()
      }
      par(mfrow = c(1,2))
      plot(vegRaster, col = viridis::viridis(nrow(VegCount)), axes = FALSE, box = FALSE)
      plot(classVeg, axes = FALSE, box = FALSE, zlim = c(0,1), main = paste("Class", t))
      
      of_indices <- c()
      unoffsetable <- c()
      restorationVal <- c() # to count the offset gain we need to add to the raster later
      e <- 0
      s <- 0
      n <- 0
      rfp <- NA
      rfp2 <- 0
      randpoint <- rasterToPoints(classVeg, spatial = TRUE)
      # 
      # if(t == 5){
      #   browser()
      # }
      randpoint <- randpoint[which(randpoint$layer >= threshold), ]
      if(nrow(randpoint) < 2){
        next
      }
      # if(all(which(!is.na(randpoint$layer))) < threshold){
      #   message(paste("None of the cells in vege class", t, "are above the threshold for offsets and these impacts are unoffsetable"))
      #   next
      # }
      
      while(is.na(rfp) || rfp2 < threshold){
        nn <- randpoint[sample(1:nrow(randpoint), 1),]
        rfp <- values(classVeg)[cellFromXY(classVeg, nn)]
        rfp2 <- values(baseRaster)[cellFromXY(baseRaster, nn)]
      }
      points(nn)
      # calculate minimum buffer
      if(type == "area"){
        resRaster <- prod(res(classVeg))
        rr <- sqrt((resRaster * offsetvalue) / 3.141593)
      } else if(type == "value"){
        rr <- offsetvalue / 2
      }
      while(s < offsetvalue){
        b <- raster::buffer(nn, rr + e)
        # plot(b, add=TRUE)
        e <- e + (res(classVeg)[1] * incrementCell)
        v <- unlist(cellFromPolygon(classVeg, b))
        v <- v[which(!is.na(values(classVeg)[v]))]
        # if(all(which(!is.na(values(classVeg)[v]))) < threshold){
        #   message(paste("None of the cells in vege class", t, "are above the threshold for offsets and these impacts are unoffsetable"))
        #   break
        # }
        if(length(v) > 0){
          for(i in v){
            if(is.null(weightRaster)){
              n <- n + 1   # Count it as an offset point
              of_indices[n] <- i # store it as an a already considered index
              restorationVal[n] <- threshold
              if(type == "value"){
                s <- s + (threshold - values(classVeg)[i]) # condition value
              } else if(type == "area"){
                s <- n
              } 
            } else{
              if(values(baseRaster)[i] >= threshold){ # value of suitability or combined raster
                if(values(classVeg)[i] <= values(baseRaster)[i]){
                  n <- n + 1
                  of_indices[n] <- i
                  restorationVal[n] <- values(baseRaster)[i]
                  if(type == "value"){
                    s <- s + (values(baseRaster)[i] - values(classVeg)[i])
                  } else if(type == "area"){
                    s <- n
                  } 
                } else{
                  unoffsetable <- append(unoffsetable, i)
                }
              } else{
                unoffsetable <- append(unoffsetable, i)
              } 
            }
            if(s >= offsetvalue){
              break
            } 
          }
        }
        searchedPoints <- unique(append(of_indices, unoffsetable))
        points(xyFromCell(classVeg, of_indices), col="blue", cex=0.1)
        values(newRaster)[searchedPoints] <- NA
        values(baseRaster)[searchedPoints] <- NA
        values(classVeg)[searchedPoints] <- NA
        if(length(which(!is.na(values(classVeg)))) <= 1){ # we may not lose all the cells because of classes
          message(paste("The offset of class", t, "has not been taken"))
          message(paste("The total remaining offset value required is ", offsetvalue - s))
          break
        }
      }
      message("Biodiversity offsetting recovered the habitat quality! :)")
      if(length(of_indices) > 1|| !is.null(of_indices)){
        gg <- append(gg, restorationVal)
        oo <- append(oo, of_indices)
      } else{
        message(paste("The offset of class", t, "has not been taken"))
        message(paste("The total remaining offset value required is ", offsetvalue - s))
      }
      message(paste("Class", t, "has been offset!", "Hooray! :))"))
    }
    message("Biodiversity offsetting recovered the habitat quality! :)")
    if(length(of_indices) > 1 || !is.null(of_indices)){
      gg <- append(gg, restorationVal)
      oo <- append(oo, of_indices)
      values(theOriginalRaster)[unlist(cellFromPolygon(theOriginalRaster, dev[z,]))] <- NA
      theOriginalRaster[of_indices] <- restorationVal
      writeRaster(theOriginalRaster, paste0("offset_patch_", z, ".tif"), overwrite = TRUE)
      message(paste("Repeat", z, "is done!", "Hooray! :))"))
    } else{
      message(paste("The offset of repeat", z, "has not been taken"))
      message(paste("The offset value is", offsetvalue))
      values(theOriginalRaster)[unlist(cellFromPolygon(theOriginalRaster, dev[z,]))] <- NA
      writeRaster(theOriginalRaster, paste0("offset_patch_", z, ".tif"), overwrite = TRUE)
      message(paste("Repeat", z, "is done!", "Hooray! :))"))
    }
  }
  cat("The number of cells used for offsetting:", length(oo), "\n")
  cat("Sum of the development value previously taken:", devValue, "\n")
  cat("Total restoration value added:", sum(gg) - sum(values(theOriginalRaster)[oo]), "\n")
  fullRaster[oo] <- gg
  # restorationRaster <- finalRaster - theOriginalRaster
  restorationRaster <- theOriginalRaster - fullRaster
  veg1 <- vegRaster
  full1 <- restorationRaster
  values(veg1) <- NA
  values(full1) <- NA
  veg1[oo] <- values(vegRaster)[oo]
  full1[oo] <- values(restorationRaster)[oo]
  
  df <- as.data.frame(setNames(stack(veg1, full1), c("class", "restored")), na.rm = TRUE)
  reportTable <- df %>% 
    group_by(class) %>% 
    summarise(sum(restored)) %>% 
    left_join(df %>% group_by(class) %>% count()) %>% 
    setNames(c("class", "restored_value", "restored_area")) %>% 
    left_join(setNames(VegCount, c("class", "required_value", "required_area"))) %>% 
    as.data.frame()
  reportTable$required_value <- reportTable$required_value * offsetMultiplier
  reportTable$required_area <- reportTable$required_area * offsetMultiplier
  reportTable$unoffset_value <- reportTable$restored_value - reportTable$required_value
  reportTable$unoffset_area <- reportTable$restored_area - reportTable$required_area
  cat("The offseting type is", type, "\n")
  if(type == "area"){
    reportTable <- reportTable[, c(1, 3, 5, 7)]
  } else{
    reportTable <- reportTable[, c(1, 2, 4, 6)]
  }
  print(reportTable)
  write.csv(reportTable, csvdir)
  # for ends ----------------------------------------------------------------
  
  print(paste("The number of cells used for offsetting:", length(oo)))
  print(paste("Sum of the development value previously taken:", devValue))
  print(paste("Sum of the offset gain added:", sum(gg) - sum(values(fullRaster)[oo])))
  # fullRaster[oo] <- gg
  return(fullRaster)
}

  