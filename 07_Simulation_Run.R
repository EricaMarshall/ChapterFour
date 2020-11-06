source("01_DEVELOPMENT_TIME_function.R")
source("01_DEVELOPMENT_TIME_AVOID_function.R")
source("02_OFFSET_TIME_function.R")
source("03_VEGE_offsetting.R")
source("04_ABUNDANCE_offsetting.R")
source("05_RICHNESS_offsetting.R")
source("06_METAPOP_offsetting.R")
library(raster)

#### Load data ##### 
Species_Map <- raster()
Condition_layer <- raster()
Vegetation_layer <- raster()
Abundance_layer <- raster()
Richness_layer <- raster()
Metapop_layer <- raster()

#### Run development scenario ##### 
for (i in 1:50) {
  raster_dev <- Developing_targeted(baseRaster = Species_Map, 
                                    weightRaster = Condition_layer,
                                    devArea = 5000, 
                                    repeats = 20,
                                    bufZone = 10000, 
                                    Threshold = 0.331) 
  # Threshold values can be species specific so developments target areas of high importance to species
  # or related to landuse types or can be 0 if you want random development
  writeRaster(raster_dev, paste0("Outputs/Developments/Development_repeat_", i, ".tif"))
}

#### Load the development files in for offsetting #### 
load_devs <- function(devRaster, n){
  require(raster)
  raster_list <- list()
  for (i in c(1:n)) {
    # browser()
    raster_list[[i]] <- raster(print(paste0(devRaster,i,".asc", sep = "")))
  }
  Dev <- stack(raster_list)
  return(Dev)
}
Development_rasters <- load_devs("Outputs/Developments/Development_repeat_", n= 50)

##### Offset these development impacts ##### 
#Vegetation area and vegetation condition
for (i in 1:50){
  raster_off <- offsetting_Veg(baseRaster = Species_Map,
                               vegRaster = Vegetation_layer,
                               weightRaster = ConLayer,
                               devRaster = Development_rasters[[i]],
                               type = "value",
                               offPatches = 20,
                               threshold = 0.331,
                               offsetMultiplier = 1,
                               incrementCell = 100)
  writeRaster(raster_off[[1]], paste0("Vegetation_Condition", i ,".tif"), sep="_", format="GTiff", overwrite=TRUE)
  
}
# habitat suitability 
for (i  in 1:50) {
  raster_off <- offsetting_timestep(baseRaster = Species_Map, # Area X SDM and Condition XSDM metrics used both baseraster and weightraster
                                    weightRaster = Condition_layer, # if you don't have a species map and are just offsetting vegetation condition 
                                    # weightraster = NULL and condition layer is the baseraster 
                                    devRaster = Development_rasters[[i]], 
                                    bufferZone = 10000,
                                    offPatches = 12, # Number of development impacts in your landscape and the number of offsets required
                                    type = "value", # This determines whether offsets will be calculated based on area or condition 
                                    # value = condition, area = area
                                    threshold = 0.331, # this can be any value, here a species specific threshold has been used to target offsets to areas of high value to the species. 
                                    #But can also be zero so offset happen randomly 
                                    offsetMultiplier = 1, incrementCell=100)
  writeRaster(raster_off, "HabitatSuitability", i, ".tif")
  
}

# Abundance 
for (i in 1:50) {
  raster_off <- offsetting_Abundance(baseRaster = Abundance_layer,
                                     newRaster = Abundance_layer*ConLayer,
                                     weightRaster = ConLayer,
                                     devRaster = Development_rasters[[i]],
                                     SDM = Species_Map, 
                                     type = "value",
                                     offPatches = 20,
                                     threshold = 0.331,
                                     offsetMultiplier = 1,
                                     incrementCell = 100)
  writeRaster(raster_off[[1]], paste0("Abundance_HS", i ,".tif"), sep="_", format="GTiff", overwrite=TRUE)
  writeRaster(raster_off[[2]], paste0("Abundance", i ,".tif"), sep="_", format="GTiff", overwrite=TRUE)
}


# Metapopulation connectivity
for (i in 1:50) {raster_off <- offsetting_Metapop(baseRaster = Metapop_layer,
                                                  newRaster = Metapop_layer*Conlayer,
                                                  weightRaster = Conlayer,
                                                  devRaster = Development_rasters[[i]],
                                                  SDM = Species_Map,
                                                  bufferZone = 10000,
                                                  offPatches = 20,
                                                  type = "value",
                                                  offsetMultiplier = 1,
                                                  threshold = 0.331,
                                                  incrementCell = 100)
writeRaster(raster_off[[1]], paste0("Metapop_Connect_HS", i, ".tif"), overwrite = TRUE)
writeRaster(raster_off[[2]], paste0("Metapop_Connect", i, ".tif"), overwrite = TRUE)
}

# Richness 
for (i in 1:50) {raster_off <- offsetting_Richness(baseRaster = Richness_layer,
                                                   weightRaster = ConLayer,
                                                   devRaster = Development_rasters[[i]],
                                                   bufferZone = 10000, 
                                                   SDM = Species_Map,
                                                   bufferZone = 10000,
                                                   offPatches = 20,
                                                   type = "value",
                                                   offsetMultiplier = 1,
                                                   threshold = 0.331,
                                                   incrementCell = 100)
writeRaster(raster_off[[1]], paste0("Richness_HS", i ,".tif"), sep="_", format="GTiff", overwrite=TRUE)
writeRaster(raster_off[[2]], paste0("Richness", i ,".tif"), sep="_", format="GTiff", overwrite=TRUE)
}

