### use alpha equal to 45 for average wet condition
### use kinematic wave aproximation with average B for alpha between 45 and 135
### calculate ck with page 287 of Chow (manning is necessary)
### calculate K as deltax/ck just like in page 302 chow
### calculate X as in page 302 chow



#' @export
circular_flow_area <- function(alpha,D) ## D is diameter of section, alpha is angle of wetted perimeter
{
    A <- ((D^2)/4)*(alpha-sin(2*alpha)/2)
    return(A)
}

#' @export
circular_wetted_angle <- function(h,D) # h is height and  D is diameter
{
    alpha <- acos(1-h/(D/2))
    return(alpha)
}

#' @export
circular_depth <- function(alpha,D)
{
    y <- D/2*(1-cos(alpha))
    return(y)
}

#' @export
circular_hydraulic_radius <- function(D,alpha) ## D is diameter, alpha is wetted angle
{
    R <- (D/4)*(1-sin(2*alpha)/(2*alpha))
    return(R)
}

#' @export
mannings_eq_flow <- function(Ks,slope,A,R) ## Ks is strickler coefficient, slope is slope of the energy curve, A is section area, R is hydraulic radius
{
    Q <- A*Ks*R^(2/3)*slope^(1/2)
    return(Q)
}

#' @export
dQdy_circular <- function(Ks,slope,A,R,y) ### for alpha between 45 and 135 I say dQdy can be approximated by a rectangular section
{
    dQdy <- Ks*slope^(1/2)*(5/3)*y^(2/3) ### as in page 287 of chow
    return(dQdy)
}

#' @export
kinematic_wave_celerity_constant_B <- function(B,dQdy) ### page 284 chow
### dQdy is the variation of flow with depth given by manning for each particular section. for the kinematic approximation the section should be close to a rectangular section
{
    ck <- (1/B)*dQdy
    return(ck)
}


#' @export
inlet_time <- function(L,n,i,S) ### time of concentration upstream of sewers, based on kinematic wave, chow 501, in minutes
{
    L <- L*3.281 ### meter to feet
    i <- i*0.0394 ### mm to inch
    inlettime <- (0.94*(L^0.6)*(n^0.6))/((i^0.4)*(S^0.3))
    return(inlettime) ### in minutes
}

#' @export
IDF <- function(D,T)
{
    require("dplyr")
    IDF <- read.table("./IDF",header=T)
    IDF <- IDF[IDF$T %in% T,]
    Daux <- D
    if(D<5)
    {
        D <- 5
        cat("\n Warning: D is lower than the defined domain for the IDF \n")
    }
    if(D>2880)
    {
        D <- 2880
        cat("\n Warning: D is greater than the defined domain for the IDF \n")
    }
    tmp <- IDF %>% filter(D0 <= D, D <= D1)
    i <- tmp$a*D^tmp$b
    P <- i*Daux/60
    return(P)
}

#' @param parm.df is SpatialPointsDataFrame pointing to the centroids of the catchment and containing C factor, n, initial loss, permanent loss, dt in s and the id of the upstream subbasin
#' @param I is the rainfall intensity in mm/h
#' @export
loss_model <- function(I,parm.df)
{
    x <- parm.df
    x@data$he <- ((I*parm.df@data$dt/3600)-parm.df@data$hi)*parm.df@data$c.factor
    x@data$he[x@data$he<0] <- 0
    
    return(x)
}


#' @param parm.df is SpatialPointsDataFrame pointing to the centroids of the catchment and containing C factor, n, initial loss, permanent loss and the id of the upstream subbasin
#' @export
Qin <- function(P,parm.df)
{

    
}

#' Q effluent from the documentation of citydrain 2
#' @param Qin is a SpatialPointsDataFrame inflows to reach in step i
#' @param V is a SpatialPointsDataFrame of volumes in reach in step i-1
#' @param Qout is a SpatialPointsDataFrame of outflows from reach in step i
#' @param parm.df is SpatialPointsDataFrame pointing to the centroids of the catchment and containing dt, K, X and the id of the upstream subbasin
#' @param i is the time step
#' @export
Qe_muskingum <- function(Qin,V,parm.df,i)
{
    Cx <- (dt/2-K*X)/(dt/2+K*(1-X))
    Cy <- 1/(dt/2+K*(1-X))
    Qout <- Cx*Qin+Cy*V
    return(Qe)
}

#' Stored volume in step i. from the documentation of citydrain 2
#' @param Qin is a SpatialPointsDataFrame inflows to reach in step i
#' @param V is a SpatialPointsDataFrame of volumes in reach in step i-1
#' @param Qout is a SpatialPointsDataFrame of outflows from reach in step i
#' @param parm.df is SpatialPointsDataFrame pointing to the centroids of the catchment and containing dt, K, X and the id of the upstream subbasin
#' @param i is the time step
#' @export
V_muskingum <- function(Qin,Qe,V,parm.df)
{ 
    Vpresent <- (Qin-Qout)*parm.df$dt[1]+V
    return(Vpresent)
}

#' Virtual volume in retention structure
#' @param Qin is a SpatialPointsDataFrame inflow to structure in step i
#' @param Qoverflow is a SpatialPointsDataFrame unregulated overflow from structure in step i
#' @param Vvirtual is the virtual volume in step i, including overflow volume
#' @param Vprevious is the volume in step i-1
#' @param param.df is a SpatialPointsDataFrame of parameters of each subbasin, including df (time step), Vmax (maximum allowable volume of structure) and Qoutmax (maximum regulated allowable outflow from structure)
#' @export
Virtual_retention <- function(Qin,Qoverflow,Vprevious,parm.df)
{
    Vvirtual <- (Qin-parm.df@data$Qoutmax)*parm.df@data$dt + Vprevious
    return(Vvirtual)
}

#' Actual volume in retention structure
#' @param Qin is a SpatialPointsDataFrame inflow to structure in step i
#' @param Qoverflow is a SpatialPointsDataFrame unregulated overflow from structure in step i
#' @param Vactual is the virtual volume in step i, including overflow volume
#' @param Vprevious is the volume in step i-1
#' @param param.df is a SpatialPointsDataFrame of parameters of each subbasin, including df (time step), Vmax (maximum allowable volume of structure) and Qoutmax (maximum regulated allowable outflow from structure)
#' @export
Actual_retention <- function(Qin,Qoverflow,Vprevious,parm.df)
{
    #'case 1 Vvirtual==0
    V <- 0

    #' case 2 Vvirtual>parm.df$Vmax
    V <- parm.df$Vmax
    #' case 3 Vvirtual >0 & Vvirtual< parm.df$Vmax
    V <- (Qin-parm.df$Qoutmax)*dt+Vprevious
    return(V)
}

#' Outflow from retention structure
#' @param Qin is a SpatialPointsDataFrame inflow to structure in step i
#' @param Qout is a SpatialPointsDataFrame regulated outflow from structure in step i
#' @param Vvirtual is the virtual volume in step i, including overflow volume
#' @param Vprevious is the volume in step i-1
#' @param param.df is a SpatialPointsDataFrame of parameters of each subbasin, including df (time step), Vmax (maximum allowable volume of structure) and Qoutmax (maximum regulated allowable outflow from structure)
#' @export
Qoutflow_ret_str <- function(Qin,Qout,Vprevious,Vvirtual,parm.df)
{
    #'case 1 Vvirtual==0
    V <- 0
    Qout <- Vprevious/dt + Qin

    #' case 2 Vvirtual>parm.df$Vmax
    Qout <- parm.df$Qoutmax
    
    #' case 3 Vvirtual >0 & Vvirtual< parm.df$Vmax
    Qout <- parm.df$Qoutmax
    
    return(Qout)
}

#' Overflow from retention structure
#' @param Qin is a SpatialPointsDataFrame inflow to structure in step i
#' @param Qout is a SpatialPointsDataFrame regulated outflow from structure in step i
#' @param Vvirtual is the virtual volume in step i, including overflow volume
#' @param Vprevious is the volume in step i-1
#' @param param.df is a SpatialPointsDataFrame of parameters of each subbasin, including df (time step), Vmax (maximum allowable volume of structure) and Qoutmax (maximum regulated allowable outflow from structure)
#' @export
Qoverflow_ret_str <- function(Qin,Qout,Vprevious,Vvirtual,parm.df)
{
    #'case 1 Vvirtual<0 Vvirtual==0
    Qover <- 0

    #' case 2 Vvirtual>parm.df$Vmax
    Qover <- Qin-parm.df$Qoutmax-(parm.df$Vmax-Vprevious)/dt
    
    #' case 3 Vvirtual >0 & Vvirtual< parm.df$Vmax
    Qover <- 0
    
    return(Qover)
}


#' define objects
#' @param Qsb
#' @export
define_obj_Qsb<- function(spObj,timeObj)
{
    Qsb <- STFDF(spObj,time=timeObj,data=data.frame(Qin=NA,Qout=NA,V=NA))    
    return(Qsb)
}

#' define objects
#' @param Qpipe
#' @export
define_obj_Qpipe<- function(spObj)
{
    Qpipe <- SpatialPointsDataFrame(coords=spObj@coords,data.frame(Qin=NA,Qout=NA,V=NA),proj4string=spObj@proj4string)
    return(Qpipe)
}

#' define objects
#' @param Qret_str
#' @export
define_obj_Qret_str<- function(spObj)
{
    Qret_str <- SpatialPointsDataFrame(coords=spObj@coords,data.frame(Qin=NA,Qout=NA,Qover=NA,V=NA),proj4string=spObj@proj4string)
    return(Qret_str)
}
