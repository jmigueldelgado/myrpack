

#' @export
startGRASS <- function(loc,ms,version)
    {

        nodename <- Sys.info()["nodename"]
        switch(nodename,
               MEKONG={
                   gis.Base6 <- "/usr/local/grass-6.4.4svn"
                   gis.Base7 <- "/usr/lib/grass70/"
                   gisDB <- "/home/delgado/grassdata"
               },
               mondego={
                   gis.Base6 <-  "/home/delgado/local/grass-6.4.6svn/"
                   gis.Base7 <-  "/usr/lib/grass70/"
                   gisDB="/home/delgado/grassdata"
               },
               vouga={
                   gis.Base6 <-  "/home/delgado/local/grass-6.4.6svn/"
                   gis.Base7 <-  "/usr/lib/grass70/"
                   gisDB="/home/delgado/grassdata"
               },
               {
                   gis.Base6 <- "/usr/lib/grass64"
                   gis.Base7 <-  "/usr/lib/grass70/"
                   gisDB <- "/home/delgado/grassdata"
               }
               )
        
        if(version>=7) {
            require(rgrass7)
                                        #            require(spgrass6)
            gis.Base <- gis.Base7
        } else
        {
            require(spgrass6);
 #           require(rgrass7)
            gis.Base <- gis.Base6
        }

        initGRASS(gisBase=gis.Base,
                  home="/home/delgado",                                            # The directory in which to create the .gisrc file; can usually be set to tempdir()
                  location=loc,                                          # GRASS location
                  mapset=ms,                                             # corresp. mapset
                  gisDbase=gisDB,  # path to grass data directory containing the location specified above and all corresp. data
                  override=TRUE)
    }


load_INMET_T <- function(id_name,path)
{
    ID <- id_name$id
    name <- id_name$name
    x <- read.table(paste0(path,'/',name,'_MG_',ID,'.txt'),skip=16,header=T,sep=';',dec=".",colClasses=c('numeric','character','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
    d <- strptime(paste(x[,2],x[,3],sep=" "),"%d/%m/%Y %H%M",tz="UTC")
    i0 <- hour(d)==0
    i12 <- hour(d)==12
  #  Tmax <- xts(x$TempMaxima[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  Tmin <- xts(x$TempMinima[i12],order.by=ymd(paste(year(d[i12]),month(d[i12]),day(d[i12]),sep="-"))) ####################
    X <- xts(x$Temp.Comp.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  H <- xts(x$Umidade.Relativa.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  V <- xts(x$Velociade.do.Vento.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  In <- xts(x$Insolacao[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  xx <- merge(Tmin,Tmax,Tmean,H,V,In,all=TRUE)
#    xx <- as.data.frame(xx)
#    xx$id <- rownames(xx)
#    xx <- melt(xx,id.vars="id")
    #xx$value <- as.numeric(xx$value,)
#    xx <- xx[!is.na(xx$value),]
    #dayofmonth <- as.numeric(substr(xx$variable,6,7))
    #xx$variable <- dayofmonth
    #xx$dates <- xx$dates+24*60*60*(xx$variable-1)
    #xx <- xx[,c(1,3)]
#    xx <- xx[with(xx,order(id)),]
                                        #    X_xts <- xts(x=xx$value,order.by=xx$id)
    colnames(X) <- ID
    return(X)
}

load_INMET_Tmin <- function(id_name,path)
{
    ID <- id_name$id
    name <- id_name$name
    x <- read.table(paste0(path,'/',name,'_MG_',ID,'.txt'),skip=16,header=T,sep=';',dec=".",colClasses=c('numeric','character','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
    d <- strptime(paste(x[,2],x[,3],sep=" "),"%d/%m/%Y %H%M",tz="UTC")
    i0 <- hour(d)==0
    i12 <- hour(d)==12
  #  Tmax <- xts(x$TempMaxima[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
    X <- xts(x$TempMinima[i12],order.by=ymd(paste(year(d[i12]),month(d[i12]),day(d[i12]),sep="-"))) ####################
#    X <- xts(x$Temp.Comp.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  H <- xts(x$Umidade.Relativa.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  V <- xts(x$Velociade.do.Vento.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  In <- xts(x$Insolacao[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  xx <- merge(Tmin,Tmax,Tmean,H,V,In,all=TRUE)
#    xx <- as.data.frame(xx)
#    xx$id <- rownames(xx)
#    xx <- melt(xx,id.vars="id")
    #xx$value <- as.numeric(xx$value,)
#    xx <- xx[!is.na(xx$value),]
    #dayofmonth <- as.numeric(substr(xx$variable,6,7))
    #xx$variable <- dayofmonth
    #xx$dates <- xx$dates+24*60*60*(xx$variable-1)
    #xx <- xx[,c(1,3)]
#    xx <- xx[with(xx,order(id)),]
                                        #    X_xts <- xts(x=xx$value,order.by=xx$id)
    colnames(X) <- ID
    return(X)
}

load_INMET_H <- function(id_name,path)
{
    ID <- id_name$id
    name <- id_name$name
    x <- read.table(paste0(path,'/',name,'_MG_',ID,'.txt'),skip=16,header=T,sep=';',dec=".",colClasses=c('numeric','character','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
    d <- strptime(paste(x[,2],x[,3],sep=" "),"%d/%m/%Y %H%M",tz="UTC")
    i0 <- hour(d)==0
    i12 <- hour(d)==12
  #  Tmax <- xts(x$TempMaxima[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
#    X <- xts(x$TempMinima[i12],order.by=ymd(paste(year(d[i12]),month(d[i12]),day(d[i12]),sep="-"))) ####################
#    X <- xts(x$Temp.Comp.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
    X <- xts(x$Umidade.Relativa.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  V <- xts(x$Velociade.do.Vento.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  In <- xts(x$Insolacao[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  xx <- merge(Tmin,Tmax,Tmean,H,V,In,all=TRUE)
#    xx <- as.data.frame(xx)
#    xx$id <- rownames(xx)
#    xx <- melt(xx,id.vars="id")
    #xx$value <- as.numeric(xx$value,)
#    xx <- xx[!is.na(xx$value),]
    #dayofmonth <- as.numeric(substr(xx$variable,6,7))
    #xx$variable <- dayofmonth
    #xx$dates <- xx$dates+24*60*60*(xx$variable-1)
    #xx <- xx[,c(1,3)]
#    xx <- xx[with(xx,order(id)),]
                                        #    X_xts <- xts(x=xx$value,order.by=xx$id)
    colnames(X) <- ID
    return(X)
}

load_INMET_In <- function(id_name,path)
{
    ID <- id_name$id
    name <- id_name$name
    x <- read.table(paste0(path,'/',name,'_MG_',ID,'.txt'),skip=16,header=T,sep=';',dec=".",colClasses=c('numeric','character','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
    d <- strptime(paste(x[,2],x[,3],sep=" "),"%d/%m/%Y %H%M",tz="UTC")
    i0 <- hour(d)==0
    i12 <- hour(d)==12
  #  Tmax <- xts(x$TempMaxima[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
    #X <- xts(x$TempMinima[i12],order.by=ymd(paste(year(d[i12]),month(d[i12]),day(d[i12]),sep="-"))) ####################
#    X <- xts(x$Temp.Comp.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  H <- xts(x$Umidade.Relativa.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  V <- xts(x$Velociade.do.Vento.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
    X <- xts(x$Insolacao[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  xx <- merge(Tmin,Tmax,Tmean,H,V,In,all=TRUE)
#    xx <- as.data.frame(xx)
#    xx$id <- rownames(xx)
#    xx <- melt(xx,id.vars="id")
    #xx$value <- as.numeric(xx$value,)
#    xx <- xx[!is.na(xx$value),]
    #dayofmonth <- as.numeric(substr(xx$variable,6,7))
    #xx$variable <- dayofmonth
    #xx$dates <- xx$dates+24*60*60*(xx$variable-1)
    #xx <- xx[,c(1,3)]
#    xx <- xx[with(xx,order(id)),]
                                        #    X_xts <- xts(x=xx$value,order.by=xx$id)
    colnames(X) <- ID
    return(X)
}


load_INMET_V <- function(id_name,path)
{
    ID <- id_name$id
    name <- id_name$name
    x <- read.table(paste0(path,'/',name,'_MG_',ID,'.txt'),skip=16,header=T,sep=';',dec=".",colClasses=c('numeric','character','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
    d <- strptime(paste(x[,2],x[,3],sep=" "),"%d/%m/%Y %H%M",tz="UTC")
    i0 <- hour(d)==0
    i12 <- hour(d)==12
  #  Tmax <- xts(x$TempMaxima[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
   # X <- xts(x$TempMinima[i12],order.by=ymd(paste(year(d[i12]),month(d[i12]),day(d[i12]),sep="-"))) ####################
#    X <- xts(x$Temp.Comp.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  H <- xts(x$Umidade.Relativa.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
    X <- xts(x$Velociade.do.Vento.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  In <- xts(x$Insolacao[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  xx <- merge(Tmin,Tmax,Tmean,H,V,In,all=TRUE)
#    xx <- as.data.frame(xx)
#    xx$id <- rownames(xx)
#    xx <- melt(xx,id.vars="id")
    #xx$value <- as.numeric(xx$value,)
#    xx <- xx[!is.na(xx$value),]
    #dayofmonth <- as.numeric(substr(xx$variable,6,7))
    #xx$variable <- dayofmonth
    #xx$dates <- xx$dates+24*60*60*(xx$variable-1)
    #xx <- xx[,c(1,3)]
#    xx <- xx[with(xx,order(id)),]
                                        #    X_xts <- xts(x=xx$value,order.by=xx$id)
    colnames(X) <- ID
    return(X)
}

load_INMET_Tmax <- function(id_name,path)
{
    ID <- id_name$id
    name <- id_name$name
    x <- read.table(paste0(path,'/',name,'_MG_',ID,'.txt'),skip=16,header=T,sep=';',dec=".",colClasses=c('numeric','character','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
    d <- strptime(paste(x[,2],x[,3],sep=" "),"%d/%m/%Y %H%M",tz="UTC")
    i0 <- hour(d)==0
    i12 <- hour(d)==12
    X <- xts(x$TempMaxima[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
#    X <- xts(x$TempMinima[i12],order.by=ymd(paste(year(d[i12]),month(d[i12]),day(d[i12]),sep="-"))) ####################
#    X <- xts(x$Temp.Comp.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  H <- xts(x$Umidade.Relativa.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  V <- xts(x$Velociade.do.Vento.Media[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  In <- xts(x$Insolacao[i0],order.by=ymd(paste(year(d[i0]),month(d[i0]),day(d[i0]),sep="-"))) ####################
  #  xx <- merge(Tmin,Tmax,Tmean,H,V,In,all=TRUE)
#    xx <- as.data.frame(xx)
#    xx$id <- rownames(xx)
#    xx <- melt(xx,id.vars="id")
    #xx$value <- as.numeric(xx$value,)
#    xx <- xx[!is.na(xx$value),]
    #dayofmonth <- as.numeric(substr(xx$variable,6,7))
    #xx$variable <- dayofmonth
    #xx$dates <- xx$dates+24*60*60*(xx$variable-1)
    #xx <- xx[,c(1,3)]
#    xx <- xx[with(xx,order(id)),]
                                        #    X_xts <- xts(x=xx$value,order.by=xx$id)
    colnames(X) <- ID
    return(X)
}
