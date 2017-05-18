
### colector pertence a bacia em que a caixa com o mesmo número está integrada. Exemplo: caixa 4 está na subbacia 2 e colector 4 considera-se que pertence à subbacia 2. id do colector é o mesmo id da caixa a montante



load_colectores <- function(cxs)
{
## cxs must be of class caixa
exclude_nas <- which(!is.na(cxs@data$flows_to))

i_to <- cxs@data$flows_to


dists <- spDists(cxs[exclude_nas,],cxs[i_to[exclude_nas],])
sl <- (cxs@data$z[exclude_nas]-cxs@data$z[i_to[exclude_nas]])/diag(dists)

col <- data.frame(colector=cxs@data$caixa[exclude_nas],
                  subbacia=cxs@data$subbacia[exclude_nas],
                  bacia=cxs@data$bacia[exclude_nas],
                  length=diag(dists),
                  slope_colector=sl,
                  flows_to_cx=i_to[exclude_nas])
return(col)

}


load_subbacia <- function(col)
{
    subbacia <- ddply(col,"subbacia",function(df)c(mean(df$slope),sum(df$length),max(df$flows_to)))
    colnames(subbacia) <- c("subbacia","slope_colector","length","flows_to_cx")
    subbacia$bacia <- rep(col$bacia[1],nrow(subbacia))
    return(subbacia)
}


load_cxs <- function(dem,spObj_cxs)
{
require("rgdal")
cxs <- readOGR(dsn="/home/delgado/Projectos/SUDS_Famalicao/CAPITULO_2/CAD/cxs/",layer="cxs")
}

load_cxs_old <- function(dem,file_profundidades,files_cxs)
{

    prof_cxs <- read.table(file_profundidades,header=T,sep="\t")

    subbacias <- unique(prof_cxs$subbacia)

    fname <- as.list(paste0(files_cxs,subbacias))



    read_f_cad <- function(fname)
    {
        
                                        #    file.pipe <- pipe(paste0("awk 'BEGIN{i=0}{i++;if (i%3==0) print $1}' < ",fname))
        f <- try(file(fname,open="r"))
        if(class(f)=="try-error"){
        } else {
            res <- readLines(f)
            close(f)
            raw <- unlist(strsplit(res[seq(3,length(res),4)],"\\  |\\="))
            raw <- as.numeric(raw[seq(2,length(raw),2)])
            cooked <- t(matrix(raw,nrow=3))
            df <- as.data.frame(cooked)
            colnames(df) <- c("x","y","z")
            spdf <- SpatialPointsDataFrame(coords=df[c("x","y")],data=df["z"],proj4string=dem@proj4string)
            return(spdf)
        }
    }
    
    cxs_list <- lapply(fname,read_f_cad)
    isnull <- !unlist(lapply(cxs_list,is.null))
    cxs <- do.call("rbind", cxs_list[isnull]) 
    
    cxs@data$caixa <- as.numeric(rownames(cxs@data))
    
    cxs@data["z"] <- over(cxs,dem)
    
    
    cxs@data <- merge(prof_cxs,cxs@data,by="caixa")
    cxs@data$cota_soleira <- cxs@data$z-cxs@data$profundidade
    ix <- is.na(cxs@data$cota_soleira)
    cxs@data$cota_soleira[ix] <- cxs@data$z[ix]-1.5
    return(cxs)
}

make_sb <- function(nomebacia) ### make subbasins
{
    sb_sp <- readVECT(paste0("bacia_",nomebacia,"_centroids"))  #import dxf to grass

    sb_sp <- sb_sp[c("cat","nome_subbacia")]
    sb_sp@data$bacia <- rep(nomebacia,length(sb_sp@data$cat))
    return(sb_sp)

}
