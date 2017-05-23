check_st_ids <- function(X)
{
    #' X must be a spacetime object
    f <- file("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.echam/ECHAM/51/pdd/P.2005-01-01.csv")
    l <- readLines(con=f,n=1)
    close(con=f)
    st <- unlist(strsplit(l,split="[,]"))
    st <- st[-1:-3]
    echamxds_ids <- as.numeric(st)

    X_ids <- as.numeric(rownames(X@sp@coords))


    ## checksum
    if(sum(X_ids %in% echamxds_ids)==length(echamxds_ids))
        cat("all ids found in station list\n")
    else
    {
        stop("error: loaded ids do not match other echamxds downscaling")
    }

    return(X[as.character(echamxds_ids),])
    
}

write_xds_input <- function(Xname,Xst)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Castanhao/StationData/",Xname,".csv")
    df <- as.data.frame(Xst)

    Y <- year(df$time)
    M <- month(df$time)
    D <- day(df$time)

    df <- cbind(Y,M,D,df)

    dfcast <- dcast(df,Y+M+D ~ sp.ID,value.var="values")
    write.table(dfcast,address,quote=FALSE,sep=",",na="-9999",col.names=TRUE,row.names=FALSE)    
}

load_xds_input <- function(Xname)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/",Xname,".csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% Xname,]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(x=Xmeta$lon[Xmeta$id %in% colnames(X)],y=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    
    X <- X[,!is.na(XY[,1])]
    XY <- XY[!is.na(XY[,1]),]
    
    X <- X[,!is.na(XY[,2])]
    XY <- XY[!is.na(XY[,2]),]
    
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "values"
    
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="UTC-3")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

loadPxds_input <- function()
{
    Xname <- "P"
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/",Xname,".csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% Xname,]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(x=Xmeta$lon[Xmeta$id %in% colnames(X)],y=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    
    X <- X[,!is.na(XY[,1])]
    XY <- XY[!is.na(XY[,1]),]
    
    X <- X[,!is.na(XY[,2])]
    XY <- XY[!is.na(XY[,2]),]
    
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "values"
    
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="UTC-3")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


loadPecham <- function(member,yeari)
{
    member <- as.numeric(member)
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Data/ECHAM/hindcast/PCP-RAD/PCP-RAD-",member,"/pcp-daily-echam46-hind8110-en51-jan",yeari,"_198101-08.nc")
    nc <- open.nc(address)
    ntime <- dim.inq.nc(nc,0)$length
    nlat <- dim.inq.nc(nc,2)$length
    nlon <- dim.inq.nc(nc,1)$length
    
    time <- as.POSIXct(var.get.nc(nc,0,start=1,count=ntime)*60*60,origin="1981-01-01 00:00:00", tz="BRT") ###hours since day 1 of january 1981 00:00:00
    lon <- var.get.nc(nc,1,start=1,count=nlon)
    lat <- var.get.nc(nc,2,start=1,count=nlat)
                                        #members <- var.get.nc(nc,2,start=1,count=nmembers)
    xloc <- var.get.nc(nc,3,start=c(1,1,1),count=c(nlon,nlat,1)) ####  note that time, lat and lon are in the wrong places...
    
    dimnames(xloc)[[1]] <- as.character(lon)
    dimnames(xloc)[[2]] <- as.character(lat)
    xloc <- melt(xloc)
    colnames(xloc) <- c("lon","lat","values")
    
    coo <- coordinates(cbind(xloc[,1],xloc[,2]))
    SP <- SpatialPoints(coords=coo, proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
                                        #    SPix <- SpatialPixels(SP)
    
    x <- var.get.nc(nc,3,start=c(1,1,1),count=c(nlon,nlat,ntime))
    
    
    dimnames(x)[[1]] <- as.character(lon)
    dimnames(x)[[2]] <- as.character(lat)
                                        #       dimnames(x)[[3]] <- as.character(members)
    dimnames(x)[[3]] <- as.character(time)
    xmelt <- melt(x)
    colnames(xmelt) <- c("lon","lat","time","value")
    xmelt$value[xmelt$value<0.1] <- 0
    
    xx <- dcast(xmelt,lon+lat  ~ time)
    cnames <- colnames(xx[3:(2+ntime)])
    xxx <- xx[,4:(2+ntime)]-xx[,3:(1+ntime)]
    xxx <- cbind(xxx,xxx[,ncol(xxx)])
    colnames(xxx) <- cnames
    xxx <- cbind(xx[,1:2],xxx)
    xxmelt <- melt(xxx,id.vars=c("lon","lat"))
    colnames(xxmelt) <- colnames(xmelt)
    xxmelt[,1] <- as.numeric(xxmelt[,1])
    xxmelt[,2] <- as.numeric(xxmelt[,2])
    
    xxxx <- xxmelt[with(xxmelt,order(time,-lat,lon)),]

    SP@coords[SP@coords[,1]>180,1] <- SP@coords[SP@coords[,1]>180,1]-360
    
    stObj <- STFDF(sp=SP, time=time, data=data.frame(values=xxxx[,4]))
    return(stObj)
}

#' @export
loadPechamxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.echam/ECHAM/",memberi,"/pdd/P.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "P",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "values"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

#' @export
loadPrsm <- function(memberi,yeari)
{
                                        #
    #"~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/echamrsm/pcp-daily-rsm97-hind8110-1981-2015/"

    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/echamrsm/pcp-daily-rsm97-hind8110-1981-2015/pcp.daily.rsm97.hind8110.jan",yeari,".",yeari,"01-",yeari,"08.nc")
#    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    nc <- open.nc(con=address,write=FALSE)

#    print.nc(nc)
    lat <- var.get.nc(nc,variable="lat")
    lon <- var.get.nc(nc,variable="lon")
    time <- var.get.nc(nc,variable="time")
    ntime <- length(time)
#    ensemble <- var.get.nc(nc,variable="ensemble")
    t0 <- att.get.nc(nc,variable="time",attribute="units")
    t0 <- paste(strsplit(t0,split=" ")[[1]][3],strsplit(t0,split=" ")[[1]][4])    
    t0 <- as.POSIXct(t0,tz="BRT")
    var <- "pcp"
    tt <- t0+days(time)    


    #var.inq.nc(nc,"lat")
    xloc <- var.get.nc(nc,variable=var,c(1,1,1,1),c(length(lon),length(lat),1,1)) #### ensemble member and time switched places,,,


    dimnames(xloc)[[1]] <- as.character(lon)
    dimnames(xloc)[[2]] <- as.character(lat)
    xloc <- melt(xloc)
    colnames(xloc) <- c("lon","lat","values")
    
    coo <- coordinates(cbind(xloc[,1],xloc[,2]))
    SP <- SpatialPoints(coords=coo, proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
                                        #    SPix <- SpatialPixels(SP)
    
    x <- var.get.nc(nc,variable=var,c(1,1,1,memberi),c(length(lon),length(lat),length(time),1)) #### ensemble member and time switched places,,,
    
    
    dimnames(x)[[1]] <- as.character(lon)
    dimnames(x)[[2]] <- as.character(lat)
                                        #       dimnames(x)[[3]] <- as.character(members)
    dimnames(x)[[3]] <- as.character(time)
    xmelt <- melt(x,id.vars=c("lon","lat"))
    colnames(xmelt) <- c("lon","lat","time","value")
        
    xx <- xmelt[with(xmelt,order(time,-lat,lon)),]

    SP@coords[SP@coords[,1]>180,1] <- SP@coords[SP@coords[,1]>180,1]-360
    
    stObj <- STFDF(sp=SP, time=tt, data=data.frame(values=xx[,4]))

    close.nc(nc)

    return(stObj)

}

#' @export
loadPmmfs <- function(member,yeari)
{
    member <- as.numeric(member)
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/mmsf/",yeari,"/","pf.01.sfc.nc")
    nc <- open.nc(address)
    
    ntime <- dim.inq.nc(nc,3)$length
    nmembers <- dim.inq.nc(nc,2)$length
    nlat <- dim.inq.nc(nc,1)$length
    nlon <- dim.inq.nc(nc,0)$length
    
    time <- as.POSIXct(var.get.nc(nc,3,start=1,count=ntime)*60*60,origin="1900-01-01 00:00:00", tz="BRT") ###hours since day 1 of january 1900 00:00:00
    lon <- var.get.nc(nc,0,start=1,count=nlon)
    lat <- var.get.nc(nc,1,start=1,count=nlat)
                                        #members <- var.get.nc(nc,2,start=1,count=nmembers)
    
    xloc <- (var.get.nc(nc,4,start=c(1,1,1,1),count=c(nlon,nlat,1,1)))
    
    dimnames(xloc)[[1]] <- as.character(lon)
    dimnames(xloc)[[2]] <- as.character(lat)
    xloc <- melt(xloc)
    colnames(xloc) <- c("lon","lat","values")
    
    coo <- coordinates(cbind(xloc[,1],xloc[,2]))
    SP <- SpatialPoints(coords=coo, proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
                                        #    SPix <- SpatialPixels(SP)
    
    x <- (var.get.nc(nc,6,start=c(1,1,member,1),count=c(nlon,nlat,1,ntime))*att.get.nc(nc,6,0) + att.get.nc(nc,6,1))*1000  ### unpacked_value = add_offset + ( packed_value * scale_factor )  from http://www.ecmwf.int/en/why-do-my-netcdf-data-only-contain-integers
    close.nc(nc)
    
    
    dimnames(x)[[1]] <- as.character(lon)
    dimnames(x)[[2]] <- as.character(lat)
                                        #       dimnames(x)[[3]] <- as.character(members)
    dimnames(x)[[3]] <- as.character(time)
    xmelt <- melt(x)
    colnames(xmelt) <- c("lon","lat","time","value")
    
    xx <- dcast(xmelt,lon+lat  ~ time)
    cnames <- colnames(xx[3:(2+ntime)])
    xxx <- xx[,4:(2+ntime)]-xx[,3:(1+ntime)]
    xxx <- cbind(xxx,xxx[,ncol(xxx)])
    colnames(xxx) <- cnames
    xxx <- cbind(xx[,1:2],xxx)
    xxmelt <- melt(xxx,id.vars=c("lon","lat"))
    colnames(xxmelt) <- colnames(xmelt)
    xxmelt[,1] <- as.numeric(xxmelt[,1])
    xxmelt[,2] <- as.numeric(xxmelt[,2])
    
    xxxx <- xxmelt[with(xxmelt,order(time,-lat,lon)),]

    SP@coords[SP@coords[,1]>180,1] <- SP@coords[SP@coords[,1]>180,1]-360
    
    stObj <- STFDF(sp=SP, time=time, data=data.frame(values=xxxx[,4]))
    return(stObj)
}

#' @export
loadTechamxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.echam/ECHAM/",memberi,"/pdd/T.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

#' @export
loadRechamxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.echam/ECHAM/",memberi,"/pdd/R.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


#' @export
loadHechamxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.echam/ECHAM/",memberi,"/pdd/H.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


#' @export
loadTechameqm <- function(memberi,yeari)
{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/echam/",memberi,"/T.",yeari,".csv")

    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

#' @export
loadHechameqm <- function(memberi,yeari)
{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/echam/",memberi,"/H.",yeari,".csv")

    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

#' @export
loadRechameqm <- function(memberi,yeari)

{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/echam/",memberi,"/R.",yeari,".csv")

    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


#' @export
loadTmmsfeqm <- function(memberi,yeari)
{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/ecmwf_mmsf/",memberi,"/T.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

#' @export
loadHmmsfeqm <- function(memberi,yeari)
{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/ecmwf_mmsf/",memberi,"/H.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


#' @export
loadRmmsfeqm <- function(memberi,yeari)
{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/ecmwf_mmsf/",memberi,"/R.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


#' @export
loadTmmsfxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.ecmwf/mmsf.LR.bc/",yeari,"/pdd/",memberi,"/T.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

#' @export
loadHmmsfxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.ecmwf/mmsf.LR.bc/",yeari,"/pdd/",memberi,"/H.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

#' @export
loadRmmsfxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.ecmwf/mmsf.LR.bc/",yeari,"/pdd/",memberi,"/R.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "T",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "value"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


#' @export
loadPmmsfxds <- function(memberi,yeari)
{
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/xds.ecmwf/mmsf.LR.bc/",yeari,"/pdd/",memberi,"/P.",yeari,"-01-01.csv")
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "P",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "values"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}


#' @export
loadHmmsf <- function(member,yeari)
{
                                        #    member <- as.numeric(member)
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/mmsf/",yeari,"/","pf.01.pl.nc")
    nc <- open.nc(address)

    
    vari <- 7 #### variable i is relative humidity
    vartemp <- 6
    nlevel <- dim.inq.nc(nc,2)$length
    leveli <- var.get.nc(nc,2,start=1,count=nlevel)
    leveli <- which(leveli==1000)
    ntime <- dim.inq.nc(nc,4)$length
    nmembers <- dim.inq.nc(nc,3)$length
    memberi <- var.get.nc(nc,3,start=1,count=nmembers)
    memberi <- which(memberi==(as.numeric(member)-1))
    nlat <- dim.inq.nc(nc,1)$length
    nlon <- dim.inq.nc(nc,0)$length
    



    time <- as.POSIXct(var.get.nc(nc,4,start=1,count=ntime)*60*60,origin="1900-01-01 00:00:00", tz="BRT") ###hours since day 1 of january 1900 00:00:00
    lon <- var.get.nc(nc,0,start=1,count=nlon)
    lat <- var.get.nc(nc,1,start=1,count=nlat)
                                        #members <- var.get.nc(nc,2,start=1,count=nmembers)
    
    xloc <- (var.get.nc(nc,vari,start=c(1,1,leveli,1,1),count=c(nlon,nlat,1,1,1)))
    
    dimnames(xloc)[[1]] <- as.character(lon)
    dimnames(xloc)[[2]] <- as.character(lat)
    xloc <- melt(xloc)
    colnames(xloc) <- c("lon","lat","values")
    
    coo <- coordinates(cbind(xloc[,1],xloc[,2]))
    SP <- SpatialPoints(coords=coo, proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
                                        #    SPix <- SpatialPixels(SP)

    #### this is the specific humidity
    x <- (var.get.nc(nc,vari,start=c(1,1,leveli,memberi,1),count=c(nlon,nlat,1,1,ntime))*att.get.nc(nc,vari,0) + att.get.nc(nc,vari,1))  ### unpacked_value = add_offset + ( packed_value * scale_factor )  from http://www.ecmwf.int/en/why-do-my-netcdf-data-only-contain-integers

    #### we need temperature to calculate the saturation vapor pressure
    temp <- (var.get.nc(nc,vartemp,start=c(1,1,leveli,memberi,1),count=c(nlon,nlat,1,1,ntime))*att.get.nc(nc,vartemp,0) + att.get.nc(nc,vartemp,1))  ### unpacked_value = add_offset + ( packed_value * scale_fact    

    #### saturation vapor pressure
    e_sat <- exp(77.3450 + 0.0057*temp - 7235/temp)/(temp^8.2) ### in Pa

    #### vapor pressure
    e <- 100000*x/(x+0.622)

    hum <- 100*e/e_sat
    dimnames(hum)[[1]] <- as.character(lon)
    dimnames(hum)[[2]] <- as.character(lat)
                                        #       dimnames(x)[[3]] <- as.character(members)
    dimnames(hum)[[3]] <- as.character(time)
    xmelt <- melt(hum)
    colnames(xmelt) <- c("lon","lat","time","value")
            
    xxxx <- xmelt[with(xmelt,order(time,-lat,lon)),]
    SP@coords[SP@coords[,1]>180,1] <- SP@coords[SP@coords[,1]>180,1]-360

    stObj <- STFDF(sp=SP, time=time, data=data.frame(values=xxxx[,4]))
    close.nc(nc)
    return(stObj)
}

#' @export
loadTmmsf <- function(member,yeari)
{
                                        #    member <- as.numeric(member)
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/mmsf/",yeari,"/","pf.01.pl.nc")
    nc <- open.nc(address)

    
    vari <- 6 #### variable i is temperature
    nlevel <- dim.inq.nc(nc,2)$length
    leveli <- var.get.nc(nc,2,start=1,count=nlevel)
    leveli <- which(leveli==1000)
    ntime <- dim.inq.nc(nc,4)$length
    nmembers <- dim.inq.nc(nc,3)$length
    memberi <- var.get.nc(nc,3,start=1,count=nmembers)
    memberi <- which(memberi==(as.numeric(member)-1))
    nlat <- dim.inq.nc(nc,1)$length
    nlon <- dim.inq.nc(nc,0)$length
    



    time <- as.POSIXct(var.get.nc(nc,4,start=1,count=ntime)*60*60,origin="1900-01-01 00:00:00", tz="BRT") ###hours since day 1 of january 1900 00:00:00
    lon <- var.get.nc(nc,0,start=1,count=nlon)
    lat <- var.get.nc(nc,1,start=1,count=nlat)
                                        #members <- var.get.nc(nc,2,start=1,count=nmembers)
    
    xloc <- (var.get.nc(nc,vari,start=c(1,1,leveli,1,1),count=c(nlon,nlat,1,1,1)))
    
    dimnames(xloc)[[1]] <- as.character(lon)
    dimnames(xloc)[[2]] <- as.character(lat)
    xloc <- melt(xloc)
    colnames(xloc) <- c("lon","lat","values")
    
    coo <- coordinates(cbind(xloc[,1],xloc[,2]))
    SP <- SpatialPoints(coords=coo, proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
                                        #    SPix <- SpatialPixels(SP)

    #### this is air temperature
    x <- (var.get.nc(nc,vari,start=c(1,1,leveli,memberi,1),count=c(nlon,nlat,1,1,ntime))*att.get.nc(nc,vari,0) + att.get.nc(nc,vari,1))  ### unpacked_value = add_offset + ( packed_value * scale_factor )  from http://www.ecmwf.int/en/why-do-my-netcdf-data-only-contain-integers

    dimnames(x)[[1]] <- as.character(lon)
    dimnames(x)[[2]] <- as.character(lat)
                                        #       dimnames(x)[[3]] <- as.character(members)
    dimnames(x)[[3]] <- as.character(time)
    xmelt <- melt(x)
    colnames(xmelt) <- c("lon","lat","time","value")
            
    xxxx <- xmelt[with(xmelt,order(time,-lat,lon)),]
    SP@coords[SP@coords[,1]>180,1] <- SP@coords[SP@coords[,1]>180,1]-360

    stObj <- STFDF(sp=SP, time=time, data=data.frame(values=xxxx[,4]))
    close.nc(nc)
    return(stObj)
}


#' @export
loadWmmsf <- function(member,yeari)
{
                                        #    member <- as.numeric(member)
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/mmsf/",yeari,"/","pf.01.pl.nc")
    nc <- open.nc(address)

    var_u <- 8
    var_v <- 9

    nlevel <- dim.inq.nc(nc,2)$length
    leveli <- var.get.nc(nc,2,start=1,count=nlevel)
    leveli <- which(leveli==1000)
    ntime <- dim.inq.nc(nc,4)$length
    nmembers <- dim.inq.nc(nc,3)$length
    memberi <- var.get.nc(nc,3,start=1,count=nmembers)
    memberi <- which(memberi==(as.numeric(member)-1))
    nlat <- dim.inq.nc(nc,1)$length
    nlon <- dim.inq.nc(nc,0)$length
    



    time <- as.POSIXct(var.get.nc(nc,4,start=1,count=ntime)*60*60,origin="1900-01-01 00:00:00", tz="BRT") ###hours since day 1 of january 1900 00:00:00
    lon <- var.get.nc(nc,0,start=1,count=nlon)
    lat <- var.get.nc(nc,1,start=1,count=nlat)
                                        #members <- var.get.nc(nc,2,start=1,count=nmembers)
    
    xloc <- (var.get.nc(nc,var_u,start=c(1,1,leveli,1,1),count=c(nlon,nlat,1,1,1)))
    
    dimnames(xloc)[[1]] <- as.character(lon)
    dimnames(xloc)[[2]] <- as.character(lat)
    xloc <- melt(xloc)
    colnames(xloc) <- c("lon","lat","values")
    
    coo <- coordinates(cbind(xloc[,1],xloc[,2]))
    SP <- SpatialPoints(coords=coo, proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
                                        #    SPix <- SpatialPixels(SP)

    #### this is U velocity
    xu <- (var.get.nc(nc,var_u,start=c(1,1,leveli,memberi,1),count=c(nlon,nlat,1,1,ntime))*att.get.nc(nc,var_u,0) + att.get.nc(nc,var_u,1))  ### unpacked_value = add_offset + ( packed_value * scale_factor )  from http://www.ecmwf.int/en/why-do-my-netcdf-data-only-contain-integers

    #### this is V velocity
    xv <- (var.get.nc(nc,var_v,start=c(1,1,leveli,memberi,1),count=c(nlon,nlat,1,1,ntime))*att.get.nc(nc,var_v,0) + att.get.nc(nc,var_v,1))  ### unpacked_value = add_offset + ( packed_value * scale_factor )  from http://www.ecmwf.int/en/why-do-my-netcdf-data-only-contain-integers

    x <- (xv^2+xu^2)^0.5
    dimnames(x)[[1]] <- as.character(lon)
    dimnames(x)[[2]] <- as.character(lat)
                                        #       dimnames(x)[[3]] <- as.character(members)
    dimnames(x)[[3]] <- as.character(time)
    xmelt <- melt(x)
    colnames(xmelt) <- c("lon","lat","time","value")
            
    xxxx <- xmelt[with(xmelt,order(time,-lat,lon)),]
    SP@coords[SP@coords[,1]>180,1] <- SP@coords[SP@coords[,1]>180,1]-360

    stObj <- STFDF(sp=SP, time=time, data=data.frame(values=xxxx[,4]))
    close.nc(nc)
    return(stObj)
}


#' @export
loadRmmsf <- function(member,yeari)
{
    member <- as.numeric(member)
    address <- paste0("/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/mmsf/",yeari,"/","pf.rad.01.sfc.nc")
    nc <- open.nc(address)

        ntime <- dim.inq.nc(nc,3)$length
        nmembers <- dim.inq.nc(nc,2)$length
        nlat <- dim.inq.nc(nc,1)$length
        nlon <- dim.inq.nc(nc,0)$length

        time <- as.POSIXct(var.get.nc(nc,3,start=1,count=ntime)*60*60,origin="1900-01-01 00:00:00", tz="BRT") ###hours since day 1 of january 1900 00:00:00
        lon <- var.get.nc(nc,0,start=1,count=nlon)
        lat <- var.get.nc(nc,1,start=1,count=nlat)
        #members <- var.get.nc(nc,2,start=1,count=nmembers)
        
    xloc <- (var.get.nc(nc,4,start=c(1,1,1,1),count=c(nlon,nlat,1,1)))
        
        dimnames(xloc)[[1]] <- as.character(lon)
        dimnames(xloc)[[2]] <- as.character(lat)
        xloc <- melt(xloc)
        colnames(xloc) <- c("lon","lat","values")
        
    coo <- coordinates(cbind(xloc[,1],xloc[,2]))
    SP <- SpatialPoints(coords=coo, proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    #    SPix <- SpatialPixels(SP)
        
    x <- (var.get.nc(nc,4,start=c(1,1,member,1),count=c(nlon,nlat,1,ntime))*att.get.nc(nc,4,0) + att.get.nc(nc,4,1))/(3600*24)  ### unpacked_value = add_offset + ( packed_value * scale_factor )  from http://www.ecmwf.int/en/why-do-my-netcdf-data-only-contain-integers
        
        
        dimnames(x)[[1]] <- as.character(lon)
        dimnames(x)[[2]] <- as.character(lat)
#       dimnames(x)[[3]] <- as.character(members)
        dimnames(x)[[3]] <- as.character(time)
        xmelt <- melt(x)
        colnames(xmelt) <- c("lon","lat","time","value")

        xx <- dcast(xmelt,lon+lat  ~ time)
        cnames <- colnames(xx[3:(2+ntime)])
        xxx <- xx[,4:(2+ntime)]-xx[,3:(1+ntime)]
        xxx <- cbind(xxx,xxx[,ncol(xxx)])
        colnames(xxx) <- cnames
        xxx <- cbind(xx[,1:2],xxx)
        xxmelt <- melt(xxx,id.vars=c("lon","lat"))
        colnames(xxmelt) <- colnames(xmelt)
        xxmelt[,1] <- as.numeric(xxmelt[,1])
        xxmelt[,2] <- as.numeric(xxmelt[,2])
        
        xxxx <- xxmelt[with(xxmelt,order(time,-lat,lon)),]

        
        stObj <- STFDF(sp=SP, time=time, data=data.frame(values=xxxx[,4]))
            close.nc(nc)
    return(stObj)
}

#' @export
loadPmmsfeqm <- function(memberi,yeari)
{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/ecmwf_mmsf/",memberi,"/P.",yeari,"-01-01.csv")
 ### loop over months   
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "P",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.numeric(t(X)))
    colnames(Vec) <- "values"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}



#' @export
loadPechameqm <- function(memberi,yeari)
{
    address <- paste0("~/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/eqm/pred/echam/"
,memberi,"/P.",yeari,".csv")
 ### loop over months   
    X <- read.table(address,sep=",",header=TRUE,check.names=FALSE,fill=TRUE,na.strings="-9999")
    addressmeta <- "/home/delgado/SESAM/sesam_data/DFG_Erkenntnis_Transfer/Climate_Prediction/xds/Ceara/StationData/"
    
    Xmeta <- read.table(paste0(addressmeta,"meta.txt"),sep=",",header=TRUE, check.names=FALSE,fill=TRUE,na.strings="-9999")
    Xmeta <- Xmeta[Xmeta$var %in% "P",]
    
    T <- X[,1:3]
    X <- X[,-1:-3]
    
    X <- X[,order(colnames(X))]
    Vec <- as.data.frame(as.vector(t(X)))
    colnames(Vec) <- "values"
                                        #        mtrx <- cbind(mtrx,Vec)
    
    
    XY <- data.frame(lon=Xmeta$lon[Xmeta$id %in% colnames(X)],lat=Xmeta$lat[Xmeta$id %in% colnames(X)],row.names=Xmeta$id[Xmeta$id %in% colnames(X)])
    XY <- XY[order(rownames(XY)),]
    ll <- SpatialPoints(XY,proj4string=CRS("+proj=longlat +ellps=WGS84 +no_defs"))
    
    t <- as.POSIXct(paste(T$Y,T$M,T$D,sep="-",tz="BRT")) 
    stObj <- STFDF(sp=ll,time=t,data=Vec)
    return(stObj)
}

