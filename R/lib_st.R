#' @export
applyThiessenWeights <- function(X,thiessen_wgs,centroids)
{
    #' @param X is a spacetime obj
    #' @param thiessen_wgs is a data.frame with  $SubbasinID is the subbasin id corresponding to the ids given in centroids and   $thiessenID is the id of the observation station nearby which must be present in the spacetime dataframe X and      $weights are the weights given by thiessen_weights.R
    
    spObj <- centroids[0,]
    st0 <- matrix()
    st0 <- st0[0]
    
    for(j in subbas@data$SubbasinID)
    {
        weights <- thiessen_wgs[thiessen_wgs$SubbasinID==j,]
        Xi <- matrix()
        Xi <- Xi[0]
        for(k in seq(1,length(weights$weights)))
        {
            Xi <- cbind(Xi,as.matrix(X[as.character(weights$thiessenID[k])]$values*weights$weights[k],ncol=1))
        }
        spObj <- spRbind(spObj,centroids[centroids@data$ids==j,])
        st0 <- cbind(st0,as.matrix(rowSums(Xi),ncol=1))
    }
    
    Xsubbas <- STFDF(sp=spObj,time=X@time,data=data.frame(values=matrix(t(st0),ncol=1)))
    return(Xsubbas)

}



#' @export
fillGapsHourly <- function(tsObj)
    {
        require(xts)
        tsFull <- seq(time(first(tsObj)),time(last(tsObj)),by="hour")
        tsFilled <- merge(tsObj,tsFull,join="right")
        tsFilled <- na.approx(tsFilled)
        colnames(tsFilled) <- colnames(tsObj)
        return(tsFilled)
    }

#' @export
fillGapsDaily <- function(tsObj)
    {
        require(xts)
        tsFull <- seq(time(first(tsObj)),time(last(tsObj)),by="day")
        
        tsFilled <- merge(tsObj,tsFull,join="right")
        tsFilled <- na.approx(tsFilled)
        colnames(tsFilled) <- colnames(tsObj)
        return(tsFilled)
    }

#' @export
apply.months <- function(xtsObj,months,aggrFun="mean",na.rm=TRUE)
{
        #xtsObj is the time series to aggregate
    #require(quantmod)
    return(apply.yearly(xtsObj[.indexmon(xtsObj) %in% (months-1)],aggrFun,na.rm=TRUE))
}


#' @export
aggrMonth <- function(stObjmon,month,aggrFun="sum")
    {
        ###mean precipitation for a given month in all stations.
        SP <- stObjmon[,.indexmon(stObjmon@time) %in% (month-1)]
        x <- as(SP, "xts")
        x <- apply(x, 2, aggrFun, na.rm=TRUE)
        return(SpatialPointsDataFrame(coords=coordinates(stObjmon@sp),data=data.frame(values=x),proj4string=CRS(proj4string(stObjmon@sp))))
    }        

#' @export
aggrMonths.yearly <- function(stObj,months,aggrFun="mean")
{
    x <- aggregate(stObj[,.indexmon(stObj@time) %in% (months-1)],by="year",FUN=aggrFun)
    return(x)
}

#' @export
interp_grid <- function(gridObj,spPointsObj)
    {
        interp <- idw(values ~ 1, locations=spPointsObj, newdata=gridObj,idp=10)
        colnames(interp@data) <- c("values","nothing")
        interp@data <- data.frame(values=interp@data$values)
        return(interp)
    }



