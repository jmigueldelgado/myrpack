getRaster <- function(name,ms)
  {
    require(rgdal)
    if (nchar(Sys.getenv("GISRC")) > 0) {
      oechoCmd <- get.echoCmdOption()       
      set.echoCmdOption(TRUE)
      print(paste("\nImporting raster from mapset",ms))
      SpGrid <- readRAST6(name,mapset = "SOTER_DB")
      set.echoCmdOption(oechoCmd)
        }
    return(SpGrid)
  }

writeRaster <- function(SpGrid,ms)
  {
    require(rgdal)
    execGRASS("g.region",region = c("bregion@PERMANENT"))
    execGRASS("r.mask",flag = c("o"), input= c("basin@PERMANENT"))
    if (nchar(Sys.getenv("GISRC")) > 0) {
      oechoCmd <- get.echoCmdOption()       
      set.echoCmdOption(TRUE)
      writeRAST6(SpGrid, colnames(SpGrid@data),overwrite=TRUE)
      set.echoCmdOption(oechoCmd)
    }
  }

rasterExists <- function(rast_name,ms)
    {
        require(rgdal)
        if (nchar(Sys.getenv("GISRC")) > 0) {
            oechoCmd <- get.echoCmdOption()       
            set.echoCmdOption(TRUE)
            rlist <- execGRASS("g.list",type = "rast",mapset = ms,intern = TRUE)
            set.echoCmdOption(oechoCmd)
        }
        yesorno <- grep(rast_name,rlist)
        return(length(yesorno))
    }

importSoter <- function(soterName)
    {
        msOrig <- "SOTER_EUdb"
        msDest <- "SOTER_DB"
        execGRASS("g.mapset", mapset=msDest)
        execGRASS("g.region",region = c("bregion@PERMANENT"))
        execGRASS("r.mask", flag = c("r"))
        execGRASS("r.proj", flag = c("overwrite"), input = soterName, location = c("EuropeanGrid"), mapset = msOrig) 
        execGRASS("r.mask",flag = c("o"), input= c("basin@PERMANENT"))
        execGRASS("r.resample",flags = c("overwrite"),input=soterName,output=paste0(soterName,"_clip"))
    }

soterExists <- function(soterName)
  {
    yes <- rasterExists(soterName,"SOTER_DB")
    if(yes)
      {
        print(paste(soterName,"exists in the database"))
      }
    else stop(paste("consider importing soil rasters to the database:",soterName))
    yes <- rasterExists(paste0(soterName,"_clip"),"SOTER_DB")
    if(yes)
      {
        print(paste0(soterName,"_clip"," exists in the database"))
      } else
        {
          ms <- "SOTER_DB"
          source("GRASS_defs.R")
          initGRASS( gisBase = gis.Base,gisDbase = gis.DBase, home = Home, location = loc, mapset = ms, override = TRUE)
          if (nchar(Sys.getenv("GISRC")) > 0) {
            oechoCmd <- get.echoCmdOption()       
            set.echoCmdOption(TRUE)
            execGRASS("g.region",region = c("bregion@PERMANENT"))
            execGRASS("r.mask",flag = c("o"), input= c("basin@PERMANENT"))
            execGRASS("r.resample",flags = c("overwrite"),input=soterName,output=paste0(soterName,"_clip"))
            set.echoCmdOption(oechoCmd) 
          }
        }
  }


get_subbasins_centroids <- function(subbasins)
    {
        COORcentroids <- coordinates(subbasins)
        rownames(COORcentroids) <- seq(1,nrow(COORcentroids))
        colnames(COORcentroids) <- c("X","Y")
        SP_centroids <- SpatialPointsDataFrame(COORcentroids, data.frame(CENTROID=subbasins@data$cat) ,proj4string=CRS(local_proj), bbox = NULL)
        return(SP_centroids)
    }



###############################
##### Interpolate over    #####
#####   Subbasins centroids ###
#####  with nearest neighbor ##
###############################

interpl_NearestNeigh <- function(subbasins,ST)
    {
        #returns spacetime dataframe with timeseries in each subbasin
        SP_centroids <- get_subbasins_centroids(subbasins)
        
        SP_st <- ST@sp
        temp <- apply(coordinates(SP_centroids),1,function(SPx) {spDistsN1(SP_st,SPx,longlat=FALSE)})
        ind <- apply(temp,2,which.min)
        ST_centroids <- ST[ind,,1]
        return(ST_centroids)
    }

