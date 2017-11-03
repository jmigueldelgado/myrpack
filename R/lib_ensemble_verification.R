#' @export
apply_rel <- function(spdf,xfc,xobs,model,lead,plower,pupper)
{
    xfc <- aggr_month_after_lead(xfc,spdf,lead)
    xobs <- aggr_month_after_lead(xobs,spdf,lead)

    p <- c(plower,pupper)
    thresh <- quantile(xobs,probs=p)
    if(p[1]==0) thresh[1] <- 0
    if(p[2]==1) thresh[2] <- 2*max(xobs)

    event_tbl <- obs_fc_table(xobs,xfc,thresh)
    xfc_bins <- get_fc_prob_bins(event_tbl)

    fpi <- getfp(xfc_bins,event_tbl)

    
    Oi <- getOi(xfc_bins,event_tbl)
    NOi <- getNOi(xfc_bins,event_tbl)


    HRi <- getHRi_rel(xfc_bins,Oi,NOi)

    binhi <- xfc_bins$binhi
    binlo <- xfc_bins$binlo
    nbins <- length(binhi)

    climatology <- get_sample_climatology(xfc_bins,event_tbl)

    three.month.abb <- c("FMA","MAM","AMJ")
       
    rel <- data.frame(HR=HRi,fp=seq(1/(2*nbins),1,1/nbins),climatology=rep(climatology,nbins),binhi=binhi,binlo=binlo,model=model,forecast_month=month.abb[2+lead],region=spdf@data[1,1],event_threshold=paste0("[",paste0(p,collapse="-"),"]"))
    return(rel)
}


#' roc of ensemble forecast
#' @param xfc is a dataframe with the forecast, each column corresponds to one ensemble member
#' @param xobs is a dataframe with the observations, the row order must be the same as in xfc
#' @param prctile_upper is the value of the upper percentile for defining an event
#' @param prctile_lower is the value of the lower percentile defining an event
#' @param thresh_upper is optional. In case no percentiles are given, a threshold value must be given
#' @param thresh_lower is optional. In case no percentiles are given, a threshold value must be given 
#' @export
roc <- function(xfc,xobs,type="cont",prctile_upper=0.33,prctile_lower=NULL,thresh_upper=NULL,thresh_lower=NULL)
{
    if(type=="cont") # compute thresholds based on input parameters
    {
        if(is.null(thresh_upper) & is.null(thresh_lower))
        {
            p <- c(prctile_lower,prctile_pupper)
            thresh <- quantile(xobs,probs=p)
            thresh_lower <- thresh[1]
            thresh_upper <- thresh[2]
            if(p[1]==0) thresh_lower <- min(xobs) - 2*max(xobs) ## a low value
            if(p[2]==1) thresh_upper <- 10*max(xobs) ## a large value
        }
        
        if(!is.null(thresh_upper))
        {
            thresh_lower <- min(xobs) - 2*max(xobs) ## a low value
        }
        
        if(!is.null(thresh_lower))
        {
            thresh_upper <- 10*max(xobs) ## a large value
        }
        
        thresh <- c(thresh_lower,thresh_upper)
        
        event_tbl <- obs_fc_table(xobs,xfc,thresh)
    }

    if(type=="prob")
    {
        if(class(xobs)=="numeric") xobs <- data.frame(xobs)
        if(class(xfc)=="numeric") xfc <- data.frame(xfc)
        event_tbl <- bind_cols(xfc,xobs)
        colnames(event_tbl) <- c("forecasted","observed")
    }
    
    xfc_bins <- get_fc_prob_bins(event_tbl)
    
    fpi <- getfp(xfc_bins,event_tbl)
    
    
    Oi <- getOi(xfc_bins,event_tbl)
    NOi <- getNOi(xfc_bins,event_tbl)
    
    HRi <- getHRi_ROC(xfc_bins,Oi)
    FARi <- getFARi_ROC(xfc_bins,NOi)
    
   # three.month.abb <- c("FMA","MAM","AMJ")
    
    roc <- data.frame(HR=HRi,FAR=FARi,fp=fpi,binlo=xfc_bins$binlo,binhi=xfc_bins$binhi,binct=xfc_bins$binct)
 
    if(sum(roc$HR==0 & roc$FAR==0)>0)
    {
        return(roc)
    }
    else
    {
        return(rbind(roc,data.frame(HR=0,FAR=0,fp=NA,binlo=NA,binhi=NA,binct=NA)))
    }

}



#' aggregate forecast/observation for January+lead by polygon
#' @export
aggr_lead_plus_3months <- function(x,pol,lead)
{
    x <- x[pol,]
    issue_month <- 2
    xx <- aggrMonths.yearly(x,seq(issue_month+lead,issue_month+lead+2),"sum")
    xx <- aggregate(xx,by=pol,mean,na.rm=TRUE)
    xx <- xx[,!(names(xx) %in% c("timeIndex"))]
   return(xx)
}

#' aggregate forecast/observation for January+lead by polygon
#' @export
aggr_month_after_lead <- function(x,pol,lead)
{
    x <- x[pol,]
    issue_month <- 2
    xx <- aggrMonths.yearly(x,issue_month+lead,"sum")
    xx <- aggregate(xx,by=pol,mean,na.rm=TRUE)
    xx <- xx[,!(names(xx) %in% c("timeIndex"))]
   return(xx)
}



#' table of observed vs forecasted events
#' @export
obs_fc_table <- function(xobs,xfc,thresh)
{
    xfc_probs <- 1*(xfc>=thresh[1] & xfc<=thresh[2])
    xfc_probs <- apply(xfc_probs,1,mean)
    xfc_probs <- melt(xfc_probs)

    xobs_event <- as.data.frame(1*(xobs>=thresh[1] & xobs<=thresh[2]))
    
    event_tbl <- cbind(xobs_event,xfc_probs)
    colnames(event_tbl) <- c("observed","forecasted")

    return(event_tbl)
}


#' get forecast probability bins from set of forecasts probabilities
#' @export
get_fc_prob_bins <- function(event_tbl)
{

    step <- 1/6
    xfc_bins <- data.frame(binlo=seq(0,1-step,step),binhi=seq(step,1,step))
    xfc_bins$binlo[1] <- -0.1
    xfc_bins$binhi[length(xfc_bins$binhi)] <- 1.1
    binct <- rep(NA,length(xfc_bins$binlo))
    binct[1] <- sum(event_tbl$forecasted>=xfc_bins$binlo[1] & event_tbl$forecasted<=xfc_bins$binhi[1])
    for(i in seq(2,length(xfc_bins$binlo)))
    {
        binct[i] <- sum(event_tbl$forecasted>xfc_bins$binlo[i] & event_tbl$forecasted<=xfc_bins$binhi[i])
    }
    xfc_bins$binct <- binct
    return(xfc_bins)
}

#' marginal distribution of the forecasts as in wilks page 285 eq 7.36 and 7.37
#' @export
get_marginal_fc_dist <- function(xfc_bins,event_tbl)
{

    Ni <- xfc_bins$binct
    nbins <- length(Ni)
    n <- sum(Ni)
    return(as.numeric(Ni/n))
}

#' conditional average observation as in wilks page 286 eq 7.38
#' @export
get_cond_av_obs <- function(xfc_bins,event_tbl)
{

    Ni <- xfc_bins$binct
    nbins <- length(Ni)
    n <- sum(Ni)
    oi <- getOi(xfc_bins,event_tbl)/Ni # Oi is occurrences. notation comes from WMO http://www.wmo.int/pages/prog/www/DPS/LRF/ATTACHII-8SVSfrom%20WMO_485_Vol_I.pdf
    oi[is.na(oi)] <- 0
    return(oi)
}

#' @export
get_sample_climatology <- function(xfc_bins,event_tbl)
{

    o <- sum(event_tbl$observed)/sum(xfc_bins$binct)

    return(o)
}

#' determine occurrences in bins
#' @export
getOi <- function(xfc_bins,event_tbl)
{
    nbins <- length(xfc_bins$binhi)
    
    Oi <- rep(NA,nbins)
    Oi[1] <- sum(event_tbl$observed[event_tbl$forecasted>=xfc_bins$binlo[1] & event_tbl$forecasted<=xfc_bins$binhi[1]])
    for(i in seq(2,nbins))
    {
        Oi[i] <- sum(event_tbl$observed[event_tbl$forecasted>xfc_bins$binlo[i] & event_tbl$forecasted<=xfc_bins$binhi[i]])
    }

    return(Oi)
}

#' determine forecast probability of each bin
#' @export
getfp <- function(xfc_bins,event_tbl)
{
    nbins <- length(xfc_bins$binhi)
    
    fp <- rep(NA,nbins)
    fp[1] <- sum(event_tbl$forecasted[event_tbl$forecasted>xfc_bins$binlo[1] & event_tbl$forecasted<=xfc_bins$binhi[1]])/xfc_bins$binct[1]
    for(i in seq(2,nbins))
    {
        fp[i] <- sum(event_tbl$forecasted[event_tbl$forecasted>xfc_bins$binlo[i] & event_tbl$forecasted<=xfc_bins$binhi[i]])/xfc_bins$binct[i]
    }
    fp[xfc_bins$binct==0] <- 0

    return(fp)
}




#' determine non-occurrences in bins
#' @export
getNOi <- function(xfc_bins,event_tbl)
{
    nbins <- length(xfc_bins$binhi)
    
    NOi <- rep(NA,nbins)

    NOi[1] <- -sum(event_tbl$observed[event_tbl$forecasted>xfc_bins$binlo[1] & event_tbl$forecasted<=xfc_bins$binhi[1]]-1)
    for(i in seq(2,nbins))
    {
        NOi[i] <- -sum(event_tbl$observed[event_tbl$forecasted>xfc_bins$binlo[i] & event_tbl$forecasted<=xfc_bins$binhi[i]]-1)
    }

    return(NOi)

}


    #' notation and expression comes from http://www.wmo.int/pages/prog/www/DPS/LRF/ATTACHII-8SVSfrom%20WMO_485_Vol_I.pdf
getHRi_rel <- function(xfc_bins,Oi,NOi)
{
    nbins <- length(xfc_bins$binhi)
    HR <- rep(NA,nbins)
    HR <- Oi/(NOi+Oi)

    HR[is.na(HR)] <- 0
    
    return(HR)
}

#' @export
getHRi_ROC <- function(xfc_bins,Oi)
{
    nbins <- length(xfc_bins$binhi)
    HR <- rep(NA,nbins)
    for(n in seq(1,nbins))
    {
        HR[n] <- sum(Oi[n:nbins])/sum(Oi)
    }

    return(HR)
}

#' @export
getFARi_ROC <- function(xfc_bins,NOi)
{
    nbins <- length(xfc_bins$binhi)
    FAR <- rep(NA,nbins)
    for(n in seq(1,nbins))
    {
        FAR[n] <- sum(NOi[n:nbins])/sum(NOi)
    }

    return(FAR)
}
