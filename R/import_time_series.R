get_tower_meteo <- function(filename)
    {
        require("xts")
        X <- list()
        header <- readLines("/home/delgado/GEOECOLOGY/user/delgado/EVAP_PT/tower/header", n = 1, encoding = "ISO-8859-15")
        X$df <- read.table(filename, fileEncoding = "ISO-8859-15", header = FALSE, fill = TRUE, blank.lines.skip = FALSE, skip = 1, sep=",", dec = ".")
        g3 <- grep('^[0-9]{3}$',X$df[,4])
        g2 <- grep('^[0-9]{2}$',X$df[,4])
        X$df[g3,4] <- paste0(0,X$df[g3,4])
        X$df[g2,4] <- paste0(0,X$df[g2,4])
        thetime <- as.POSIXct(strptime(apply(X$df[,2:4],1,toString),format="%Y,%j,%H%M"))
        colnames(X$df) <-  strsplit(header,",")[[1]]
        X$ts <- as.xts(x = X$df[,-2:-4], order.by = thetime)
        return(X$ts)
    }


getTFriedrichs <- function(filename)
    {
        #filename <- "/home/delgado/GEOECOLOGY/user/delgado/EVAP_PT/HS/COMBILOG5.LOG"
        require("xts")
        X <- list()
        header <- readLines(filename, n = 8, encoding = "ISO-8859-15")
        X$df <- read.table(filename, fileEncoding = "ISO-8859-15", header = FALSE, fill = TRUE, blank.lines.skip = FALSE, skip = 8, sep="\t", dec = ",")
        thetime <- as.POSIXct(X$df[,2], format = "%d.%m.%Y %H:%M:%S")
        header[6] <- gsub(" ","_",header[6])
        colnames(X$df) <-  strsplit(header[6],"\t")[[1]]
        X$ts <- as.xts(x = X$df[,-1:-2], order.by = thetime)
        #X$header <- rbind(header[6],header[8])
        return(X$ts)
    }

getSNIRH <- function(filename)
    {
      require("xts")
      X <- list()
      X$df <- read.table(filename, fileEncoding = "ISO-8859-15", header = TRUE, fill = TRUE, blank.lines.skip = FALSE, skip = 3, sep="\t")
      X$df <- X$df[-(match(c(TRUE),X$df[,1]==""):nrow(X$df)),]
      thetime <- as.POSIXct(X$df[,1], format = "%d/%m/%Y %H:%M")
      X$ts <- xts(x = X$df[,seq(2,ncol(X$df)-1,2)], order.by = thetime)
    }

getWASA <- function(filename)
    {
        require("xts")
        X <- list()
        X$df <- read.table(filename, header = TRUE, fill = TRUE, blank.lines.skip = FALSE, skip = 1)
        if(grepl("X",colnames(X$df)[3])) #daily data
           {
               thetime <- as.POSIXct(strptime(apply(X$df[,1:2],1,toString),format="%Y, %j"))
               X$ts <- xts(x = X$df[,3], order.by = thetime)
           }
           else { #hourly data
               thetime <- as.POSIXct(strptime(apply(X$df[,1:3],1,toString),format="%Y, %j, %H"))
               X$ts <- xts(x = X$df[,4], order.by = thetime)
           }
return(X$ts)
    }
