#### this script:
##        for use change point estimation for phenological transition dates in Xie&Wilson RSE 2020
##

library(segmented)
library(foreach)
library(parallel)
library(doParallel)
library(dplyr)
library(raster)
library(ncdf4)

# set working directory
fmosaic=""          # folder of mosaic MODIS nc files
fevi=""             # folder of evi nc files
fdata=""            # relevant datasets
fseg=""             # output folder

# functions----------------------------------------------------------------------
# get twice daily netCDF files info
filesinfo <-function(fmosaic)
{
  fs=data.frame(path=list.files(fmosaic,full=T,pattern=".nc"),file=list.files(fmosaic,pattern=".nc"),stringsAsFactors=F) # build file list
  fs[,c("sensor","datecode")]=do.call("rbind",strsplit(fs$file,".",fixed=T))[,c(1,2)]
  fs$sensordate=paste(fs$sensor,fs$datecode,sep=".")
  fs$date=as.Date(substr(fs$file,10,16),"%Y%j")
  fs$month=format(fs$date,"%m")
  fs$year=format(fs$date,"%Y")
  fs$time=ifelse(fs$sensor=="MOD09GA","10:30","13:30")
  fs$datetime=as.POSIXct(paste(fs$date,fs$time))
  fs$dateid=format(fs$date,"%Y%m%d")
  return(fs)
}

# make input dataframe
makedf=function(EVI,doyt,...){
  data=cbind.data.frame(EVI=EVI,doyt=doyt)
  data=na.omit(data)
  return(data)
}

# change point estimation functions
seg_fun <- function(data,na.rm=TRUE,...)
{
  #profvis({
  #data=makedf(EVI,doyt)
  #return(data$EVI[1:5])
  evi.lm<-tryCatch(lm(EVI~doyt,data=data), error= function(e) return(0))    # catch error from all NAs
  if (class(evi.lm)[1]=="numeric") {
    return(rep(NA,23)) } else
    {
      #data=cbind.data.frame(EVI=EVI,doyt=doyt)
      #data=na.omit(data)
      spr=log_spr(data)       # logistic curve fitting for spring
      aut=log_aut(data)       # logistic curve fitting for autumn
      p=c(spr[5],spr[6],aut[5],aut[6])      # set logistic curve fitting results as initial estimation values
      evi.seg<-tryCatch(segmented(evi.lm,seg.Z=~doyt,psi=p,
                                  control=seg.control(toll = 1e-08, it.max=100, h=1, stop.if.error=T, n.boot=300)),
                        error= function(e) return(0))      # catch error from failed segmented function
      if (class(evi.seg)[1] == "numeric") {
        #print(paste("Segmented estimation Fail!",sep=""))
        return(rep(-9999,23))} else
        {
          #print(paste("Segmented estimation Succeed!",sep=""))
          df=numeric(23)
          df=c(c(evi.seg$psi[,c(2,3)]), c(intercept(evi.seg)$doyt), c(slope(evi.seg)$doyt[,c(1,2)]))
          #varnames=c("greenup","maturity","senescence","dormancy","greenup_se","maturity_se","senescence_se","dormancy_se",
          #           "intercept1","intercept2","intercept3","intercept4","intercept5",
          #           "slope1","slope2","slope3","slope4","slope5","slope1_se","slope2_se","slope3_se","slope4_se","slope5_se")
          #structure(df, names=varnames)
          return(df)      # return a vector with 23 numbers
        }
    }
  # })
}

# run code --------------------------------------------------------------------------
# read location information for 20000 pixels
grd=read.table(file=paste(fdata,"/loc_sample20000.csv",sep=""),header=T,sep= ",",stringsAsFactors=FALSE)

# read evi nc files (one nc file stores evi values of all pixels at one time point)
fs2=filesinfo(fmosaic)        # get files info
fs2=fs2[order(fs2$datetime),]

filenames=paste(fevi,"/",fs2$sensordate,"_evi_QAls.nc",sep="")
d=stack(x=c(filenames))            # read all files

# output 23 variable names
varnames=c("greenup","maturity","senescence","dormancy","greenup_se","maturity_se","senescence_se","dormancy_se",
           "intercept1","intercept2","intercept3","intercept4","intercept5",
           "slope1","slope2","slope3","slope4","slope5","slope1_se","slope2_se","slope3_se","slope4_se","slope5_se")

year=2000:2019

grd=modloc_grdob   # for ground observation sites
for (j in 1:19)
{
  # get time series of day and time including current year and January in next year
  datetime=c(fs2[which(fs2$year==year[j]),10], fs2[which(fs2$year==year[j+1] & fs2$month=="01"),10])
  if (j>1) doyt=as.numeric(datetime-as.POSIXct(paste(year[j],"-01-01",sep="")))/24 else
    doyt=as.numeric(datetime-as.POSIXct(paste(year[j],"-01-01",sep="")))
  bands=c(which(fs2$year==year[j]),which(fs2$year==year[j+1] & fs2$month=="01"))

  # extract evi for 20000 pixels
  df=extract(d,grd$cellid,df=T,layer=bands[1],n=length(bands))

  # parallel computation
  registerDoParallel(30)
  getDoParWorkers()

  # change point estimation
  seg <- foreach(i=1:length(grd$cellid), .combine='rbind.data.frame') %dopar%
    {
      EVI=as.numeric(df[i,-1])
      data=makedf(EVI,doyt)
      seg_fun(data)
    }
   
  # stop parallel workers
  stopCluster()
  
  # assign variable names
  colnames(seg)=varnames

  # save output for each year
  write.table(seg,file=paste(fseg,"/seg_",year[j],"_smp20000.csv",sep=""),row.names = F, sep=",")
}

