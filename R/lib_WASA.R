
######################333 define locals!!!!!


locals <- function(cenario,nome,wasa_home,lump_home)
    {
        l <- list()
        l$cenario <- cenario
        l$cenario0 <- "default"
        l$wasa_home <- wasa_home
        l$lump_home <- lump_home
        l$lump_out <- paste(l$lump_home,cenario,nome,"LUMP_out",sep="/")
        l$lump_out0 <- paste(l$lump_home,l$cenario0,nome,"LUMP_out",sep="/")
        l$soil_data <- paste(l$lump_home,"out_t",sep="/")
        l$db_out <- paste(l$lump_home,cenario,nome,"db_out",sep="/")
        l$dbname <- "mysqlitedb"

        l$wasa_input_dir <- paste(l$wasa_home,cenario,nome,"Input/",sep="/")

        l$wasa_output_dir <- paste(l$wasa_home,cenario,nome,"Output/",sep="/")
        return(l)
    }

##################3
##################
#################

edit_dodat <- function(wasa_input_dir,wasa_output_dir)
    {
        dodat <- paste(wasa_input_dir,"do.dat",sep="/")
        f <- file(dodat,open="r")
        text <- readLines(f)
        close(f)
        
        text[2] <- wasa_input_dir
        text[3] <- wasa_output_dir
        text[4] <- "2003     //tstart (start year of simulation)"
        text[5] <- "2006     //tstop (end year of simulation)"
        text[6] <- "8     //mstart (start month of simulation)"
        text[7] <- "7     //mstop (end month of simulation)"
        text[21] <- ".t.     //dohour: do hourly version"
        text[30] <- "1     //dt: timestep in hours"
        text[36] <- ".t.     //save state of storage  to files  after siulation period (optional)"
        text[37] <- ".t.     //save state of storage  to files  after siulation period (optional)"
        
### write new do.dat
        f <- file(dodat,open="w+")
        write(text,file=f,sep="\n")
        close(f)
    }


get_tc_area <- function(l,sbi=1)
    {
        wasa_input_dir <- l$wasa_input_dir
        LUMP_tables <- l$lump_out
        lu <- getLUinfo(input_dir)
        tc <- get_tc_properties(input_dir)
        tclist <- unique(tc["tc_id"])
        lulist <- unique(tc["lu_id"])

        horizons <- read.table(paste0(LUMP_tables,"/horizons.dat"),header=TRUE) ### read soil properties and obtain porosity of soilto calculate available water holding capacity of LU
        theta_s_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        volume_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        area_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        
#        sbi <- 1
        
        for(tci in tclist[,1])
            {
                soil_in_tc <- tc$soil_id[tc$tc_id %in% tci]
                TCarea <- tc$TC_area[tc$tc_id %in% tci]
                soil_area <- tc$frac[tc$tc_id %in% tci]*TCarea
                lui <- unique(tc$lu_id[tc$tc_id %in% tci])
                theta_s_in_soil <- numeric()
                depth_in_soil <- numeric()
                for(soili in soil_in_tc)
                    {
                        hor_in_soil <- horizons[horizons$soil_id %in% soil_in_tc,]
                        theta_s_in_soil <- c(theta_s_in_soil,sum(hor_in_soil$theta_s*hor_in_soil$thickness)/sum(hor_in_soil$thickness))
                        depth_in_soil <- c(depth_in_soil,sum(hor_in_soil$thickness))
                    }
                theta_s_in_tc <- sum(tc$frac[tc$tc_id %in% tci]*theta_s_in_soil)
                depth_in_tc <- sum(tc$frac[tc$tc_id %in% tci]*depth_in_soil)
                
                position <- tc$TCposition[tc$tc_id %in% tci]

                switch(
                    position[1],
                    theta_s_lu$cabeceira[theta_s_lu$pid==lui] <- theta_s_in_tc,
                    theta_s_lu$encosta[theta_s_lu$pid==lui] <- theta_s_in_tc,
                    theta_s_lu$vale[theta_s_lu$pid==lui] <- theta_s_in_tc
                    )
                switch(
                    position[1],
                    volume_lu$cabeceira[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000,
                    volume_lu$encosta[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000,
                    volume_lu$vale[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000 #### in m3
                    )
                switch(
                    position[1],
                    area_lu$cabeceira[area_lu$pid==lui] <- TCarea[1]*1000000,
                    area_lu$encosta[area_lu$pid==lui] <- TCarea[1]*1000000,
                    area_lu$vale[area_lu$pid==lui] <- TCarea[1]*1000000 ### in m2
                    )
                
            }
        
        volume_lu$total <- apply(volume_lu[c("vale","encosta","cabeceira")],1,sum)
        area_lu$total <- apply(area_lu[c("vale","encosta","cabeceira")],1,sum)
        
        theta_s_lu$total <- apply(area_lu[c("vale","encosta","cabeceira")]*theta_s_lu[c("vale","encosta","cabeceira")],1,sum)/area_lu$total
        return(area_lu)
    }

get_tc_volume <- function(l,sbi=1)
    {
        wasa_input_dir <- l$wasa_input_dir
        LUMP_tables <- l$lump_out
        cat(LUMP_tables,"\n")
        lu <- getLUinfo(l)
        tc <- get_tc_properties(l)
        tclist <- unique(tc["tc_id"])
        lulist <- unique(tc["lu_id"])

        horizons <- read.table(paste0(LUMP_tables,"/horizons.dat"),header=TRUE) ### read soil properties and obtain porosity of soilto calculate available water holding capacity of LU
        theta_s_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        volume_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        area_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        
#        sbi <- 1
        
        for(tci in tclist[,1])
            {
                soil_in_tc <- tc$soil_id[tc$tc_id %in% tci]
                TCarea <- tc$TC_area[tc$tc_id %in% tci]
                soil_area <- tc$frac[tc$tc_id %in% tci]*TCarea
                lui <- unique(tc$lu_id[tc$tc_id %in% tci])
                theta_s_in_soil <- numeric()
                depth_in_soil <- numeric()
                for(soili in soil_in_tc)
                    {
                        hor_in_soil <- horizons[horizons$soil_id %in% soil_in_tc,]
                        theta_s_in_soil <- c(theta_s_in_soil,sum(hor_in_soil$theta_s*hor_in_soil$thickness)/sum(hor_in_soil$thickness))
                        depth_in_soil <- c(depth_in_soil,sum(hor_in_soil$thickness))
                    }
                theta_s_in_tc <- sum(tc$frac[tc$tc_id %in% tci]*theta_s_in_soil)
                depth_in_tc <- sum(tc$frac[tc$tc_id %in% tci]*depth_in_soil)
                
                position <- tc$TCposition[tc$tc_id %in% tci]

                switch(
                    position[1],
                    theta_s_lu$cabeceira[theta_s_lu$pid==lui] <- theta_s_in_tc,
                    theta_s_lu$encosta[theta_s_lu$pid==lui] <- theta_s_in_tc,
                    theta_s_lu$vale[theta_s_lu$pid==lui] <- theta_s_in_tc
                    )
                switch(
                    position[1],
                    volume_lu$cabeceira[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000,
                    volume_lu$encosta[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000,
                    volume_lu$vale[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000 #### in m3
                    )
                switch(
                    position[1],
                    area_lu$cabeceira[area_lu$pid==lui] <- TCarea[1]*1000000,
                    area_lu$encosta[area_lu$pid==lui] <- TCarea[1]*1000000,
                    area_lu$vale[area_lu$pid==lui] <- TCarea[1]*1000000 ### in m2
                    )
                
            }
        
        volume_lu$total <- apply(volume_lu[c("vale","encosta","cabeceira")],1,sum)
        area_lu$total <- apply(area_lu[c("vale","encosta","cabeceira")],1,sum)
        
        theta_s_lu$total <- apply(area_lu[c("vale","encosta","cabeceira")]*theta_s_lu[c("vale","encosta","cabeceira")],1,sum)/area_lu$total
        return(volume_lu)
    }

get_tc_theta_s <- function(l,sbi=1)
    {
        wasa_input_dir <- l$wasa_input_dir
        LUMP_tables <- l$lump_out
        lu <- getLUinfo(l)
        tc <- get_tc_properties(l)
        tclist <- unique(tc["tc_id"])
        lulist <- unique(tc["lu_id"])

        horizons <- read.table(paste0(LUMP_tables,"/horizons.dat"),header=TRUE) ### read soil properties and obtain porosity of soilto calculate available water holding capacity of LU
        theta_s_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        volume_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        area_lu <- data.frame(pid=lu$id,total=rep(NA,nrow(lu)),vale=rep(NA,nrow(lu)),encosta=rep(NA,nrow(lu)),cabeceira=rep(NA,nrow(lu))) #### in LUMP TC count from bottom to top of hillslope
        
#        sbi <- 1
        
        for(tci in tclist[,1])
            {
                soil_in_tc <- tc$soil_id[tc$tc_id %in% tci]
                TCarea <- tc$TC_area[tc$tc_id %in% tci]
                soil_area <- tc$frac[tc$tc_id %in% tci]*TCarea
                lui <- unique(tc$lu_id[tc$tc_id %in% tci])
                theta_s_in_soil <- numeric()
                depth_in_soil <- numeric()
                for(soili in soil_in_tc)
                    {
                        hor_in_soil <- horizons[horizons$soil_id %in% soil_in_tc,]
                        theta_s_in_soil <- c(theta_s_in_soil,sum(hor_in_soil$theta_s*hor_in_soil$thickness)/sum(hor_in_soil$thickness))
                        depth_in_soil <- c(depth_in_soil,sum(hor_in_soil$thickness))
                    }
                theta_s_in_tc <- sum(tc$frac[tc$tc_id %in% tci]*theta_s_in_soil)
                depth_in_tc <- sum(tc$frac[tc$tc_id %in% tci]*depth_in_soil)
                
                position <- tc$TCposition[tc$tc_id %in% tci]

                switch(
                    position[1],
                    theta_s_lu$cabeceira[theta_s_lu$pid==lui] <- theta_s_in_tc,
                    theta_s_lu$encosta[theta_s_lu$pid==lui] <- theta_s_in_tc,
                    theta_s_lu$vale[theta_s_lu$pid==lui] <- theta_s_in_tc
                    )
                switch(
                    position[1],
                    volume_lu$cabeceira[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000,
                    volume_lu$encosta[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000,
                    volume_lu$vale[volume_lu$pid==lui] <- 0.001*depth_in_tc*TCarea[1]*1000000 #### in m3
                    )
                switch(
                    position[1],
                    area_lu$cabeceira[area_lu$pid==lui] <- TCarea[1]*1000000,
                    area_lu$encosta[area_lu$pid==lui] <- TCarea[1]*1000000,
                    area_lu$vale[area_lu$pid==lui] <- TCarea[1]*1000000 ### in m2
                    )
                
            }
        
        volume_lu$total <- apply(volume_lu[c("vale","encosta","cabeceira")],1,sum)
        area_lu$total <- apply(area_lu[c("vale","encosta","cabeceira")],1,sum)
        
        theta_s_lu$total <- apply(area_lu[c("vale","encosta","cabeceira")]*theta_s_lu[c("vale","encosta","cabeceira")],1,sum)/area_lu$total
        return(theta_s_lu)
    }

check_soil_veg_dat <- function(wasa_input_dir)
    {
        f <- file(paste0(wasa_input_dir,"/Hillslope/soil_vegetation.dat"))
        lines <- readLines(f)
        
        f.wr <- file(paste0(wasa_input_dir,"/Hillslope/soil_vegetation.dat.new"))
        writeLines(readLines(f,n=3),f.wr)
        f.wr <- file(paste0(wasa_input_dir,"/Hillslope/soil_vegetation.dat.new"),open="a")
        
        count <- 1
        for(i in seq(4,length(lines),1))
            {
                count <- count-1
                
                linei <- scan(f,skip=i-1,nlines=1)
                SVCs <- linei[(5+1):(5+linei[5])]
                if(count==0)
                    {
                        nbrSVCs <- length(which(!is.na(SVCs)))
                        count <- 3
                    }
                cat(nbrSVCs,"\n")
##### update table with corrected values
                updtdline <- c(linei[1:4],nbrSVCs,SVCs[1:nbrSVCs])
                
##### append to new soil_vegetation.dat
                write(updtdline,f.wr,ncolumns=length(updtdline),append=TRUE,sep="\t")                
            }
##### replace soil_vegetation with soil_vegetation.new
        system(paste0("rm ",wasa_input_dir,"/Hillslope/soil_vegetation.dat"))
        system(paste0("mv ",wasa_input_dir,"/Hillslope/soil_vegetation.dat.new ",wasa_input_dir,"/Hillslope/soil_vegetation.dat"))
    }

getLUinfo <- function(l)
    {
        wasa_input_dir <- l$wasa_input_dir
        flu <- file(paste0(wasa_input_dir,"/Hillslope/hymo.dat"))
        lulines <- readLines(flu)
        lus <- data.frame(subbas_id=numeric(),subbas_area=numeric(),id=numeric(),areakm2=numeric(),LUfrac=numeric())
        for(i in seq(3,length(lulines)))
            {
                linei <- scan(flu,skip=i-1,nlines=1)
                nLU <- linei[3]
                area <- linei[2]
                subbas_id <- data.frame(subbas_id=rep(linei[1],nLU))
                subbas_area <- data.frame(subbas_area=rep(area,nLU))
                LUarea <- data.frame(areakm2=area*linei[(3+nLU+1):(3+2*nLU)])
                LUid <- data.frame(id=linei[(3+1):(3+nLU)])
                LUfrac <- data.frame(frac=linei[(3+nLU+1):(3+2*nLU)])
                lus <- rbind(lus,cbind(subbas_id,subbas_area,LUid,LUarea,LUfrac))            
            }
        return(lus)
    }


getTCfrac <- function(l)
    {
        wasa_input_dir <- l$wasa_input_dir
        ftc <- file(paste0(wasa_input_dir,"/Hillslope/terrain.dat"))
        tclines <- readLines(ftc)
        tcs <- data.frame(id=numeric(),frac=numeric(),position=numeric())

        for(i in seq(3,length(tclines)))
            {
                linei <- scan(ftc,skip=i-1,nlines=1)
                id <- linei[1]
                frac <- linei[2]
                position <- linei[4]
                tcs <- rbind(tcs,cbind(id,frac,position))            
            }
        return(tcs)
        
    }
#

get_tc_properties <- function(l)
{
    wasa_input_dir <- l$wasa_input_dir
    
    f <- file(paste0(wasa_input_dir,"/Hillslope/soil_vegetation.dat"))
    svclines <- readLines(f)


    TCfrac <- getTCfrac(l)
    LUinfo <- getLUinfo(l)
    
            
    getsoils <- function(i)
        {
            linei <- scan(f,skip=i-1,nlines=1)
            soils <- linei[(5+1):(5+linei[5])]
            return(soils)
        }
    getvegs <- function(i)
        {
            linei <- scan(f,skip=i-1,nlines=1)
            vegs <- linei[(5+1):(5+linei[5])]
            return(vegs)
        }
    getfracs <- function(i) 
        {
            linei <- scan(f,skip=i-1,nlines=1)
            fracs <- linei[(5+1):(5+linei[5])]
            return(fracs)
        }
    getnbrSVCs <- function(i) 
        {
            linei <- scan(f,skip=i-1,nlines=1)
            nbrSVCs <- linei[5]
            return(nbrSVCs)
        }
    getsubbasID <- function(i) 
        {
            linei <- scan(f,skip=i-1,nlines=1)
            subbasID <- linei[1]
            return(subbasID)
        }
    getLUid <- function(i) 
        {
            linei <- scan(f,skip=i-1,nlines=1)
            LUid <- linei[2]
            return(LUid)
        }
    getLUarea <- function(i)
        {
            id <- getLUid(i)
            area <- LUinfo$areakm2[LUinfo$id==id]
            return(area)
        }
    getTCid <- function(i) 
        {
            linei <- scan(f,skip=i-1,nlines=1)
            TCid <- linei[3]
            return(TCid)
        }
    getTCarea <- function(i)
        {
            id <- getTCid(i)
            area <- TCfrac$frac[TCfrac$id==id]*LUareai
            return(area)
        }
    getTCposition <- function(i)
        {
            id <- getTCid(i)
            position <- TCfrac$position[TCfrac$id==id]
            return(position)
        }
    
    tbl <- data.frame(subbas_id=numeric(),LU_id=numeric(),LUarea=numeric(),TC_id=numeric(),TCarea=numeric(),soils=numeric(),vegs=numeric(),fracs=numeric(),TCposition=numeric())

    for(i in seq(4,length(svclines),3))
        {
            nbrSVCs <- getnbrSVCs(i)
            soils <- data.frame(soil_id=getsoils(i))
            vegs <- data.frame(veg_id=getvegs(i+1))
            fracs <- data.frame(frac=getfracs(i+2))
            subbas <- data.frame(subbas_id=rep(getsubbasID(i),nbrSVCs))

            lui <- getLUid(i)
            LU <- data.frame(lu_id=rep(lui,nbrSVCs))
            
            LUareai <- getLUarea(i)
            LUarea <- data.frame(LU_area=rep(LUareai,nbrSVCs))

            tci <- getTCid(i)
            TC <- data.frame(tc_id=rep(tci,nbrSVCs))

            TCareai <- getTCarea(i)
            TCarea <- data.frame(TC_area=rep(TCareai,nbrSVCs))

            TCposition <- getTCposition(i)
            tbl <- rbind(tbl,cbind(subbas,LU,LUarea,TC,TCarea,soils,vegs,fracs,TCposition))           
        }
    return(tbl)
}

substrRight <- function(Str,start=NULL,Stop=0)
    {
        substr(Str, nchar(Str)-start+1, nchar(Str)-Stop)
    }


read_tc_theta_hourly<- function(l)
    {
        require("xts")
        output_dir <- l$wasa_output_dir
        df <- read.table(paste0(output_dir,"tc_theta.out"),skip=1,header=TRUE)        
        y <- df[,1]
        doy <- df[,2]
        h <- df[,3]
        df <- df[,-1:-3]
        xdate <- strptime(x=paste(y,doy,h,sep="-"),"%Y-%j-%H",tz="BRT")
        colnames(df) <- substr(colnames(df),2,nchar(colnames(df)))
        xtsObj <- xts(df,as.POSIXlt(xdate),tz="BRT")
        return(xtsObj)
    }

read_RiverFlow_daily <- function(l,spObj=NULL)
    {
        output_dir <- l$wasa_output_dir
        df <- read.table(paste0(output_dir,"/River_Flow.out"),skip=1,header=TRUE)
        y <- df[,1]
        doy <- df[,2]
        df <- df[,-1:-2]
        xdate <- strptime(x=paste(y,doy,sep="-"),"%Y-%j")
        if(is.null(colnames(df)))
            {
                xtsObj <- xts(df,as.POSIXlt(xdate))
                colnames(xtsObj) <- "1"
            } else
                {
                    colnames(df) <- substr(colnames(df),2,nchar(df))
                    xtsObj <- xts(df,as.POSIXlt(xdate))
                }
        return(xtsObj)
    }

read_qhorton_hourly <- function(l)
    {
        output_dir <- l$wasa_output_dir
        df <- read.table(paste0(output_dir,"/qhorton.out"),skip=1,header=TRUE)
        y <- df[,1]
        doy <- df[,2]
        h <- df[,3]
        df <- df[,-1:-3]
        xdate <- strptime(x=paste(y,doy,sep="-"),"%Y-%j",tz="BRT")
        colnames(df) <- substr(colnames(df),2,nchar(df))
        xtsObj <- xts(df,as.POSIXlt(xdate),tz="BRT")
        if(is.null(spObj))
        {
            return(xtsObj)
        } else
        {
            x <- coredata(xtsObj)
            xtime <- data.frame(time=time(xtsObj))
            xx <- cbind(xtime,x)
            xtsmelt <- melt(xx,id.vars="time")
            names(xtsmelt) <- c("time","ID","values")
            xtsmelt$ID <- as.numeric(xtsmelt$ID)
            xxx <- join(xtsmelt,as.data.frame(spObj),by="ID")
            xxxx <- xxx[with(xxx,order(time,ID)), ]
            stObj <- STFDF(sp=spObj,time=unique(xxxx$time),data=xxxx[c("values")])
            return(stObj)
        }
    }

read_rainfall_daily <- function(l)
    {
        input_dir <- l$wasa_input_dir
        cols <- scan(paste0(input_dir,"/rain_daily.dat"),nlines=1,skip=2)
        cols <- cols[-1:-2]
        x <- read.table(paste0(input_dir,"/rain_daily.dat"),header=FALSE,skip=3)
        xdate <- x[,1]
        xyear <- substrRight(xdate,4)
        xmonth <- substrRight(xdate,6,4)
        xday <- substrRight(xdate,8,6)
        xdatetime <- as.POSIXlt(ISOdate(year=xyear,month=xmonth,day=xday,tz="GMT"))
        dfObj <- as.data.frame(x[,-1:-2])
        colnames(dfObj) <- paste0("X",as.character(cols))
        xtsObj <- xts(x=dfObj,order.by=xdatetime)

        return(xtsObj)
    }


read_rainfall_hourly <- function(l)
    {
        input_dir <- l$wasa_input_dir
        cols <- scan(paste0(input_dir,"/Time_series/rain_hourly.dat"),nlines=1,skip=2)
        cols <- cols[-1:-2]
        x <- read.table(paste0(input_dir,"/Time_series/rain_hourly.dat"),header=FALSE,skip=3)
        xdate <- x[,1]
        xhour <- x[,2]
        xyear <- substrRight(xdate,4)
        xmonth <- substrRight(xdate,6,4)
        xday <- substrRight(xdate,8,6)
        xdatetime <- as.POSIXlt(ISOdate(year=xyear,month=xmonth,day=xday,hour=xhour,tz="GMT"))
        dfObj <- as.data.frame(x[,-1:-2])
        xtsObj <- xts(x=dfObj,order.by=xdatetime)
        colnames(xtsObj) <- as.character(cols)
        return(xtsObj)
    }



spacetime2WASA_T <- function(stObj,address)
{
    df <- as(stObj,"data.frame")
    mtx <- t(matrix(df$values,ncol=length(stObj@time),nrow=length(stObj@sp@data$ids)))

    try(system(paste0("mkdir ",address)))
    try(system(paste0("rm ",address,"/temperature.dat")))
    fileConn <- file(paste0(address,"/temperature.dat"),"a")
    cat("Daily average temperature (in degree Celsius)","Date\tNo. of days\tSubbasin-ID.", file = fileConn, sep = "\n")

    thetime <- list()
    thetime$posix <- time(stObj@time)
    thetime$str <- strftime(thetime$posix,"%d%m%Y")
    thetime$index <- seq(1,length(thetime$str),1)
    dfObj1 <- data.frame(thetime$str,thetime$index)
    colnames(dfObj1) <- c("0","0")
    dfObj2 <- as.data.frame(mtx)
    colnames(dfObj2) <- as.numeric(stObj@sp@data$ids)
    dfObj <- data.frame(dfObj1,format(dfObj2,digits=0,nsmall=1, scientific=FALSE))
    HEADER <- c(colnames(dfObj1),colnames(dfObj2))
    cat(HEADER, file = fileConn, sep = "\t")
    cat("\n",file = fileConn, sep = "")
    write.table(dfObj,file = fileConn, quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE, fileEncoding = "UTF-8")
    close(fileConn)            
}

spacetime2WASA_R <- function(stObj)
{
    df <- as(stObj,"data.frame")
    mtx <- t(matrix(df$values,ncol=length(stObj@time),nrow=length(stObj@sp@data$ids)))

    fileConn <- file("./radiation.dat","a")
    cat("Daily average radiation [W/m2]","Date\tNo. of days\tSubbasin-ID.", file = fileConn, sep = "\n")

    thetime <- list()
    thetime$posix <- time(stObj@time)
    thetime$str <- strftime(thetime$posix,"%d%m%Y")
    thetime$index <- seq(1,length(thetime$str),1)
    dfObj1 <- data.frame(thetime$str,thetime$index)
    colnames(dfObj1) <- c("0","0")
    dfObj2 <- as.data.frame(mtx)
    colnames(dfObj2) <- as.numeric(stObj@sp@data$ids)
    dfObj <- data.frame(dfObj1,format(dfObj2,digits=0,nsmall=1, scientific=FALSE))
    HEADER <- c(colnames(dfObj1),colnames(dfObj2))
    cat(HEADER, file = fileConn, sep = "\t")
    cat("\n",file = fileConn, sep = "")
    write.table(dfObj,file = fileConn, quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE, fileEncoding = "UTF-8")
    close(fileConn)            
}


spacetime2WASA_H <- function(stObj,address)
{
    df <- as(stObj,"data.frame")
    mtx <- t(matrix(df$values,ncol=length(stObj@time),nrow=length(stObj@sp@data$ids)))

    try(system(paste0("mkdir ",address)))
    try(system(paste0("rm ",address,"/humidity.dat")))
    fileConn <- file(paste0(address,"/humidity.dat"),"a")
    cat("Daily average humidity [%]","Date\tNo. of days\tSubbasin-ID.", file = fileConn, sep = "\n")

    thetime <- list()
    thetime$posix <- time(stObj@time)
    thetime$str <- strftime(thetime$posix,"%d%m%Y")
    thetime$index <- seq(1,length(thetime$str),1)
    dfObj1 <- data.frame(thetime$str,thetime$index)
    colnames(dfObj1) <- c("0","0")
    dfObj2 <- as.data.frame(mtx)
    colnames(dfObj2) <- as.numeric(stObj@sp@data$ids)
    dfObj <- data.frame(dfObj1,format(dfObj2,digits=0,nsmall=1, scientific=FALSE))
    HEADER <- c(colnames(dfObj1),colnames(dfObj2))
    cat(HEADER, file = fileConn, sep = "\t")
    cat("\n",file = fileConn, sep = "")
    write.table(dfObj,file = fileConn, quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE, fileEncoding = "UTF-8")
    close(fileConn)            
}

#' writes WASA input file
#' @param stObj is a spacetime object whose points are the centroids of the subbasins
#' @param address is the path where the file should be written
#' @export
spacetime2WASA_P <- function(stObj,address)
{
    thetime <- list()
    
    df <- as(stObj,"data.frame")

    if("cat" %in% colnames(df)) colnames(df)[colnames(df)=="cat"] <- "id"
    if("ID" %in% colnames(df)) colnames(df)[colnames(df)=="ID"] <- "id"
    
    df$id <- as.factor(df$id)
    mtx <- dcast(df[c("time","cat","value")],time ~ id)
    mtx <- mtx[,colnames(mtx)!=c("time")]
    thetime$posix <- time(stObj@time)
    try(system(paste0("mkdir ",address)))
    try(system(paste0("rm ",address,"/rain_daily.dat")))
    fileConn <- file(paste0(address,"/rain_daily.dat"),"a")
    cat("Daily total precipitation [mm] for each subasin, ordered according to Map-IDs","Date\t\tSubbasin-ID.", file = fileConn, sep = "\n")

    thetime$str <- strftime(thetime$posix,"%d%m%Y")
    thetime$index <- seq(1,length(thetime$str),1)
    dfObj1 <- data.frame(thetime$str,thetime$index)
    colnames(dfObj1) <- c("0","0")
    dfObj2 <- mtx

#    colnames(dfObj2) <- as.numeric(stObj@sp@data$ids)
    dfObj <- data.frame(dfObj1,format(dfObj2,digits=0,nsmall=1, scientific=FALSE))
    HEADER <- c(colnames(dfObj1),colnames(dfObj2))
    cat(HEADER, file = fileConn, sep = "\t")
    cat("\n",file = fileConn, sep = "")
    write.table(dfObj,file = fileConn, quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE, fileEncoding = "UTF-8")
    close(fileConn)
    return(address)
}



writePrecDaily <- function(tsObj)
    {
        require(xts)
        tsObj <- fillGapsDaily(tsObj)
        fileConn <- file("./rain_daily.dat","a")
        cat("Daily average precipitation [mm/d]","Date\tNo. of days\tSubbasin-ID.", file = fileConn, sep = "\n")
        
        thetime <- list()
        thetime$posix <- time(tsObj)
        thetime$str <- strftime(thetime$posix,"%d%m%Y")
        thetime$index <- seq(1,length(thetime$str),1)
        dfObj1 <- data.frame(thetime$str,thetime$index)
        colnames(dfObj1) <- c("0","0")
        dfObj2 <- as.data.frame(tsObj)
        colnames(dfObj2) <- colnames(tsObj)
        dfObj <- data.frame(dfObj1,dfObj2)
        HEADER <- c(colnames(dfObj1),colnames(dfObj2))
        cat(HEADER, file = fileConn, sep = "\t")
        cat("\n",file = fileConn, sep = "")
        write.table(dfObj,file = fileConn, quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE, fileEncoding = "UTF-8")
        close(fileConn)             
    }



