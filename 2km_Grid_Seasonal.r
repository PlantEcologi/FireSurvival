# In this script are 3 primary sections
#  1  Create grid and collect geographic data
#  2  Prepare the temperature, precipitation, and soil moisture data
#  3  Prepare the fire data

library(foreign);library(zoo);library(PBSmapping); library(maptools); library(RNetCDF);library(rgdal)
source("/media/Data/Work/Stat/Scripts/Functions.r")

# Location of files
griddir <- function(x) paste("/media/Data/Work/Regional/ZA/CFR/FireAnalysis/Grid/",x,sep="")
climatedir <- function(x) paste("/media/Data/Work/Regional/ZA/Climate/",x,sep="")
datadir <- function(x) paste("/media/Data/Work/Regional/ZA/CFR/FireAnalysis/ModelData_2km/",x,sep="")


#Period of interest
    start=as.Date("6/1/1951", "%m/%d/%Y")  #must be a year later to include the first years' moving average
    startearly=as.Date("6/1/1950", "%m/%d/%Y")  #this is for the calculation of the previous 12 months of weather
    stop=as.Date("5/31/1999", "%m/%d/%Y")

# Some parameters
    newyears=6                  #Start of the new fire year.  This puts it just before winter.

### Set up grid
      domainx=c(17,26)  #define domain to speed things up for testing
      domainy=c(-35,-32)
      resolution=0.02 #this is the size of the grid cells 0.02 ~ 2km
      lon=seq(domainx[1],domainx[2],resolution)  
      lat=seq(domainy[1],domainy[2],resolution)
      grid=makeGrid(lon,lat,addSID=F,projection="LL")
      attr(grid,"PolyData")=calcArea(grid,rollup=1) #add area as attribute - sizes vary slightly due to geographic (WGS1984) projection

### Identify which cells to include
### First merge with Reserves
### Import Reserve map
      reserves=importShapefile("/media/Data/Work/Regional/ZA/CFR/CFR_reserves.shp",projection="LL")  #read in reserves
	x=subset(attr(reserves,"PolyData"),select=c("PID","Region","Include")) #get attribute table
        x=x[x$Include==1,]  # identify PIDs with with the "Include" field equal to "which"
        x$Region=as.factor(as.character(x$Region))
        reserves=reserves[reserves$PID%in%x$PID,] #subset polygons
        reserves=SpatialPolygons2PolySet(unionSpatialPolygons(PolySet2SpatialPolygons(reserves),x$Region))  #convert to sp, union by region, convert back to PolySet
        attr(reserves,"PolyData")=x #update attributes
        attr(reserves,"PolyData")$regionid=as.numeric(attr(reserves,"PolyData")$Region)
        attr(reserves,"PolyData")$PID=1:nrow(attr(reserves,"PolyData"))
        regiontable=unique(subset(attr(reserves,"PolyData"),select=c("Region","regionid")));regiontable=regiontable[order(regiontable$regionid),];rownames(regiontable)=1:nrow(regiontable)
        write.csv(regiontable,datadir("regiontable.csv"),row.names=F)
        writeOGR(SpatialPolygonsDataFrame(PolySet2SpatialPolygons(reserves),data=as.data.frame(regiontable)),datadir(""), "reserves", driver="ESRI Shapefile")

### Intersect grid with regions
intersectgrid <- function(grid,reserves,percentkeep){
### intersects the grid with the regions and keeps cells with at least 'percentkeep' of their cells in the region
	gridid=expand.grid(gid=unique(grid$PID),regionid=unique(as.numeric(reserves$PID)))  #identify all combinations to create new PIDs
	gridid$PID=1:nrow(gridid) #set up ids to result from intersect
	gridint=joinPolys(grid,reserves,"INT",maxVert=1e+8);gc()  #intersect 
	areas=calcArea(gridint,rollup=1) #calculate area of intersection
	areas=merge(areas,gridid,by="PID",all.x=T) #merges area with gid and regionid
	for(i in as.numeric(names(table(areas$gid)[table(areas$gid)>1]))){ #identify which cells are on a boundary between regions and assign it to region with maximum area
		drop=which.min(areas[areas$gid==i,2]) #identify the minimum
		keep=which.max(areas[areas$gid==i,2]) #identify the maximum
		pidremove=areas[areas$gid==i,1][drop] #identify which to remove
		pidkeep=areas[areas$gid==i,1][keep] #identify which to keep
		add=areas[areas$gid==i,2][drop] #get area to assign to pidkeep 
		areas=areas[!areas$PID==pidremove,] #remove value from areas table
		areas[areas$PID==pidkeep,2]=areas[areas$PID==pidkeep,2]+add #add area to pidkeep
	}
	areas=merge(areas,attr(grid,"PolyData"),by.x="gid",by.y="PID",all.x=T) #add original area attributes
	areas=areas[(areas$area.x/areas$area.y)>=percentkeep,] #drop grids with less than 'percentkeep' area in region
	grid2=grid[grid$PID%in%areas$gid,]  #subset original grid to selected polygons
	attr(grid2,"PolyData")=subset(areas,select=c("gid","regionid","area.y")) #define new attributes
	#attr(grid2,"PolyData")=subset(areas,select=c(PID="PID",regionid="regionid",area="area.x")) #define new attributes
	grid2$PID=as.numeric(as.factor(grid2$PID))# create new PID 1:n
	attr(grid2,"PolyData")$gid=as.numeric(as.factor(attr(grid2,"PolyData")$gid)) #update PIDs in attribute table
	attr(grid2,"PolyData")=attr(grid2,"PolyData")[order(attr(grid2,"PolyData")$gid),]
	rownames(attr(grid2,"PolyData"))=attr(grid2,"PolyData")$gid
	centers=calcCentroid(grid2,rollup=1)
	attr(grid2,"PolyData")=merge(attr(grid2,"PolyData"),centers,by.x="gid",by.y="PID") #update attribute table with new things
#	attr(grid2,"PolyData")=merge(attr(grid2,"PolyData"),attr(reserves,"PolyData"),by.x="regionid",by.y="regionid") #add reserve name
#	attr(grid2,"PolyData")=merge(attr(grid2,"PolyData"),regiontable,by="regionid") #add reserve name
	return(grid2)
}
grid=intersectgrid(grid,reserves,0.5); 	paste("Total cells created: ",length(unique(grid$PID)))  #run the function and replace the original grid
plotPolys(grid,col=attr(grid,"PolyData")$regionid)
writeOGR(SpatialPolygonsDataFrame(PolySet2SpatialPolygons(grid),data=as.data.frame(attr(grid,"PolyData"))),datadir(""), "grid", driver="ESRI Shapefile")

grid=importShapefile(datadir("grid.shp"),projection="LL")

					#plotPolys(grid,col=attr(grid,"PolyData")$regionid) 
static=attr(grid,"PolyData")[,-1]; colnames(static)[match(c("area.y","X","Y"),colnames(static))]=c("area","lon","lat")


#Intersected grid and vegrtype in grass, here we read in dbf file from vectors and compute the maximum
#veg=read.dbf("/media/Data/Work/grassdata/CFR/CFR/dbf/veggridfire.dbf")
#veg=subset(veg,select=c("a_GID","b_MAPCODE","area","b_NAME"));colnames(veg)=c("gid","mapcode","area","name")
#x=reshape(veg,timevar="mapcode",idvar="gid",direction="wide",drop="name");colnames(x)=gsub("area.","",colnames(x))
#x[is.na(x)]=0
x$eco=colnames(x)[-1][max.col(x[,2:ncol(x)]) ]
x=subset(x,select=c("gid","eco"))
	# eliminate small categories
      t=table(x$eco)  #make table
      t=names(t[t<10]) #identify which have less than 10 cells to lump as "other"
      x[x$eco%in%t,2]="Other"
ecotable=unique(cbind(as.character(veg$name),as.character(veg$mapcode)));colnames(ecotable)=c("econame","ecocode")
ecotable=rbind(ecotable,c("Other","Other"))

x=merge(static,x,by.x="gidold",by.y="gid",all.x=T);x=subset(x,select=c("gid","eco"))
x[is.na(x$eco),2]="Other"
x=merge(x,ecotable,by.x="eco",by.y="ecocode",all.x=T);x=x[order(x$gid),]
x$econame=as.factor(as.character(x$econame))
vegcodes=x; rm(x); 
vegcodes$ecoid=as.numeric(as.factor((vegcodes$eco)))
vegcodes=subset(vegcodes,select=c("eco","gid","econame","ecoid"))
write.csv(vegcodes,datadir("vegcodes.csv"))  #export to combine with other data

### chose majority and "mix" categories
### link to gids from static


#######################################################################################################
# Section two: Generate climate indices for points of interest

#get lat/lon for each point
#Get climate data out for each point
#summarize for months of interest

#write seasonal table to cfr
      x=expand.grid(year=as.numeric(format(start,"%Y")):as.numeric(format(stop,"%Y")),month=1:12)
      x$decyear=x[,1]+x[,2]/12-15/365
      x=x[!(x$year==as.numeric(format(start,"%Y"))&x$month<as.numeric(format(start,"%m"))),]  #get rid of missing months at beginning
      x=x[!(x$year==as.numeric(format(stop,"%Y"))&x$month>as.numeric(format(stop,"%m"))),]  #get rid of missing months at end
      season=as.data.frame(cbind(month=1:12,season=c(rep(4,2),rep(1,3),rep(2,3),rep(3,3),4),
          seasonname=c(rep("summer",2),rep("fall",3),rep("winter",3),rep("spring",3),"summer")))
      season=merge(x,season,by="month")
      season=season[order(season$decyear),]
      rownames(season)=1:nrow(season)
      season$seasondecyear=season$year+(as.numeric(season$season)/4)-(45/365)
      season$seasondecyear=ifelse(season$month%in%1:2,season$seasondecyear-1,season$seasondecyear)
      season$fireyear=ifelse(season$month<newyears,season$year-1,season$year)
      season$sid=match(season$seasondecyear,unique(season$seasondecyear))
      season$date=as.Date(paste(season$year,season$month,15,sep="/"),"%Y/%m/%d")
      write.csv(season,datadir("season.csv"),row.names=F)


tmax=open.nc(climatedir("tmax.nc"))  # Opens connection with netcdf file
tmin=open.nc(climatedir("tmin.nc"))  # Opens connection with netcdf file
precip=open.nc(climatedir("precip.nc"))  # Opens connection with netcdf file

#create matrix of date (year,month,day)
      dates=as.data.frame(utcal.nc("days since 1950-01-01", var.get.nc(ncfile=tmax,"time"), type="n")[,1:3])
      #dates2=as.character(as.Date(paste(dates[,2],"/",dates[,3],"/",dates[,1],sep=""), "%m/%d/%Y"))
      daystart1=c(as.numeric(format(start,"%Y")),as.numeric(format(start,"%m")),as.numeric(format(start,"%d")))
      daystart=which(dates[,1]==daystart1[1]&dates[,2]==daystart1[2]&dates[,3]==daystart1[3])
      daystartearly1=c(as.numeric(format(startearly,"%Y")),as.numeric(format(startearly,"%m")),as.numeric(format(startearly,"%d")))
      daystartearly=which(dates[,1]==daystartearly1[1]&dates[,2]==daystartearly1[2]&dates[,3]==daystartearly1[3])
      daystop1=c(as.numeric(format(stop,"%Y")),as.numeric(format(stop,"%m")),as.numeric(format(stop,"%d")))
      daystop=which(dates[,1]==daystop1[1]&dates[,2]==daystop1[2]&dates[,3]==daystop1[3])
      numdays=as.numeric(stop-start)
  #add fire year to meteo dates
      dates$fireyear=ifelse(dates$month<newyears,dates$year-1,dates$year)

#subset to dates of interest
    dates2=dates[daystart:(daystart+numdays),]
    dates2=merge(subset(season,select=c("year","fireyear","month","sid")),dates2)
    dates2=dates2[order(dates2$year,dates2$month),]
    dates=dates2 ; rm(dates2)
#save dates to disk
    write.csv(dates,datadir("dates.csv"),row.names=F)

##############
#Identify which meteo grids go with which 1km grids
#Temperature
    tlon=var.get.nc(tmax,"longitude")
    tlat=var.get.nc(tmax,"latitude")
    plon=var.get.nc(precip,"longitude")
    plat=var.get.nc(precip,"latitude")

#Identify which cells in the tmax, tmin, and precip datasets (netcdf) are completely null
# Add lat / lon to tid and pid tables to use for spatial random effects.

  #Temp
    # Get first days worth of data and identify which cells are null
        miss=var.get.nc(tmax,"t",start=c(1,1,daystart),count=c(length(tlon),length(tlat),1))
        tind=which(!is.na(miss),arr.ind=T)
        tpos=cbind(lonind=tind[,1],lon=tlon[tind[,1]],latind=tind[,2],lat=tlat[tind[,2]])
            # go through each 
              tlatid=matrix(ncol=5,nrow=nrow(static))
                for(i in 1:nrow(static)) {
                tlatid[i,1]=static$gid[i]
                tlatid[i,2:5]=tpos[which.min((abs(round(tpos[,4],3)-static$lat[i]))+(abs(round(tpos[,2],3)-static$lon[i]))),] }
                colnames(tlatid)=c("gid","lonind","lon","latind","lat")
                tid=unique(tlatid[,2:5]);tid=cbind(tid=1:nrow(tid),tid);colnames(tid)=c("tid","lonind","lon","latind","lat") #get unique TID values
                write.csv(as.data.frame(tid),datadir("tid.csv"),row.names=F) #save unique values
                tidgid=merge(tid,tlatid) #merge to get tid-gid connection
                static=merge(static,subset(tidgid,select=c("tid","gid")))

					#Precip
    # Get first days worth of data and identify which cells are null
    miss=var.get.nc(precip,"ppt",start=c(1,1,daystart),count=c(length(plon),length(plat),1))
        pind=which(!is.na(miss),arr.ind=T)
        ppos=cbind(lonind=pind[,1],lon=plon[pind[,1]],latind=pind[,2],lat=plat[pind[,2]])
        platid=matrix(ncol=5,nrow=nrow(static))
          for(i in 1:nrow(static)) {
          platid[i,1]=static$gid[i]
          platid[i,2:5]=ppos[which.min((abs(round(ppos[,4],3)-static$lat[i]))+(abs(round(ppos[,2],3)-static$lon[i]))),] }
          colnames(platid)=c("gid","lonind","lon","latind","lat")
          pid=unique(platid[,2:5]);pid=cbind(pid=1:nrow(pid),pid); colnames(tid)=c("tid","lonind","lon","latind","lat")
        write.csv(as.data.frame(pid),datadir("pid.csv"),row.names=F)
                pidgid=merge(pid,platid) #merge to get tid-gid connection
                static=merge(static,subset(pidgid,select=c("pid","gid")))

      #Clean up
      rm(miss);rm(tpos);rm(tind);rm(ppos);rm(pind);rm(tlatid);rm(platid);rm(pidgid);rm(tidgid)
  #save static to disk
      write.csv(as.data.frame(static),datadir("static.csv"),row.names=F)


###########################
# Extract data from the NetCDF files
# Get seasonal meteo data for each point


#Load Data
      library(foreign);library(RNetCDF);library(zoo)
      tmax=open.nc(climatedir("tmax.nc"))  # Opens connection with netcdf file
      tmin=open.nc(climatedir("tmin.nc"))  # Opens connection with netcdf file
      precip=open.nc(climatedir("precip.nc"))  # Opens connection with netcdf file
      static=read.csv(datadir("static.csv"))
      dates=read.csv(datadir("dates.csv"))
      tid=read.csv(datadir("tid.csv"))
      pid=read.csv(datadir("pid.csv"))
      gc() #prepare for some large datasets
     #dimensions are [lon,lat,time]
     tmaxdaily=var.get.nc(tmax,"t")
     tmaxdata=matrix(ncol=nrow(tid),nrow=length(daystartearly:(daystop)))
     for(i in 1:nrow(tid)){
	     tmaxdata[,i]=tmaxdaily[tid[i,2],tid[i,4],c(daystartearly:daystop)]
     }
     rm(tmaxdaily);gc();tmaxdaily=tmaxdata ; rm(tmaxdata); gc()
### Minimum Temperature
     tmindaily=var.get.nc(tmin,"t")
     tmindata=matrix(ncol=nrow(tid),nrow=length(daystartearly:(daystop)))
     for(i in 1:nrow(tid)){
	     tmindata[,i]=tmindaily[tid[i,2],tid[i,4],c(daystartearly:daystop)]
     }
     rm(tmindaily);gc();tmindaily=tmindata ; rm(tmindata); gc()
### Combine max, min, & ave in one tempdaily
     tmindaily[is.na(tmindaily)]=mean(tmindaily,na.rm=T)
     tmaxdaily[is.na(tmaxdaily)]=mean(tmaxdaily,na.rm=T)
     tempave=(tmaxdaily+tmindaily)/2
     tempannualave=as.matrix(rollapply(zoo(tempave),width=365,mean,na.pad=T,align="right")) #right aligned rolling mean - 10minutes
     temp6monthave=as.matrix(rollapply(zoo(tempave),width=180,mean,na.pad=T,align="right")) #right aligned rolling mean - 10minutes
     temphotweek=as.matrix(rollapply(zoo(tmaxdaily),width=7,mean,na.pad=T,align="right")) #right aligned rolling mean - 10minutes
     tempave=tempave[(daystart-daystartearly+1):nrow(tempave),]#cut to final sid values
     tmindaily=tmindaily[(daystart-daystartearly+1):nrow(tmindaily),]
     tmaxdaily=tmaxdaily[(daystart-daystartearly+1):nrow(tmaxdaily),]
     tempannualave=tempannualave[(daystart-daystartearly+1):nrow(tempannualave),]
     temp6monthave=temp6monthave[(daystart-daystartearly+1):nrow(temp6monthave),]
     temphotweek=temphotweek[(daystart-daystartearly+1):nrow(temphotweek),]
     tempdaily=stack(as.data.frame(tmaxdaily)); colnames(tempdaily)=c("tmax","tid")
     tempdaily$tid=as.numeric(gsub("V","",tempdaily$tid)) # get rid of "V"s in tid field
     tempdaily$tmin=stack(as.data.frame(tmindaily))[,1] ; rm(tmindaily);gc()
     tempdaily$sid=dates$sid
     tempdaily$annualave=stack(as.data.frame(tempannualave))[,1] ; rm(tempannualave);gc()
     tempdaily$sixmonthave=stack(as.data.frame(temp6monthave))[,1] ; rm(temp6monthave);gc()
     tempdaily$hotweek=stack(as.data.frame(temphotweek))[,1] ; rm(temphotweek);gc()
     tempdaily$tave=(tempdaily$tmax+tempdaily$tmin)/2
     write.csv(tempdaily,datadir("tempdaily.csv"),row.names=F)
     tempannualaveold=t(as.matrix(apply(tempave,2,function(x) tapply(x,list(as.factor(dates$fireyear)),mean,na.rm=T)))) #annual mean by year
     write.csv(tempannualaveold,datadir("tempannualaveold.csv"),row.names=F)

### Precipitation
     gc()
     precipdaily=var.get.nc(precip,"ppt")*.1 # load data and to convert to mm
     precipdata=matrix(ncol=nrow(pid),nrow=length(daystartearly:(daystop)))
     for(i in 1:nrow(pid)){
	     precipdata[,i]=precipdaily[pid[i,2],pid[i,4],c(daystartearly:daystop)]
     }
     precipdaily=precipdata ; rm(precipdata); gc()
     precipannualtotal=as.matrix(rollapply(zoo(precipdaily),width=365,sum,na.pad=T,align="right")) #right aligned rolling mean - 10minutes
     precip6monthtotal=as.matrix(rollapply(zoo(precipdaily),width=180,sum,na.pad=T,align="right")) #right aligned rolling mean - 10minutes
     precipdaily=precipdaily[(daystart-daystartearly+1):nrow(precipdaily),]
     precipannualtotal=precipannualtotal[(daystart-daystartearly+1):nrow(precipannualtotal),]
     precip6monthtotal=precip6monthtotal[(daystart-daystartearly+1):nrow(precip6monthtotal),]
     precipdaily=stack(as.data.frame(precipdaily)); colnames(precipdaily)=c("precip","pid")
     precipdaily$annualtot=stack(as.data.frame(precipannualtotal))[,1] ; rm(precipannualtotal);gc()
     precipdaily$sixmonthtot=stack(as.data.frame(precip6monthtotal))[,1] ; rm(precip6monthtotal);gc()
     precipdaily$pid=gsub("V","",precipdaily$pid) # get rid of "V"s in pid field
     precipdaily$sid=dates$sid #add sid
     write.csv(precipdaily,datadir("precipdaily.csv"),row.names=F)

    ##################################
    # Generate indices from daily data
    #Temp
tempdaily=read.csv(datadir("tempdaily.csv"))
precipdaily=read.csv(datadir("precipdaily.csv"))
dates=read.csv(datadir("dates.csv"))

conc=tapply(precipdaily$precip,list(precipdaily$pid),pptconc,months=dates$month,years=dates$year)
conc=cbind(pptconc=conc,pid=as.numeric(names(conc)))
write.csv(conc,datadir("conc.csv"))

## Plots of the daily data to see what the seasonal averages do to the data...
           smoothfunctionsplot <- function(tempdata,precipdata, id,xlim=c(1,400),lines=T) {
		   par(mfrow=c(2,1))
		   x=tempdata[tempdata$tid==id,];x$time=1:nrow(x)
		   plot(x$tave~x$time,type="n",xlim=xlim,ylim=c(0,40),xlab="Days",ylab="Temperature")
		   if(lines==T) abline(v=seq(xlim[1],xlim[2],90),col="grey",lty="dashed")
		   polygon(c(x$time,rev(x$time)),c(x$tmin,rev(x$tmax)),col="grey",border=NA)
		   lines(x$tave,lwd=1.5)
		   lines(x$annualave,col="red")
		   lines(x$sixmonthave,col="green")
		   x=precipdata[precipdata$pid==id,];x$time=1:nrow(x)
		   plot(scale(x$precip)~x$time,type="h",xlim=xlim,ylim=c(-2,12),xlab="Days",ylab="Precipitation",col="grey")
		   if(lines==T) abline(v=seq(xlim[1],xlim[2],90),col="grey",lty="dashed")
		   lines(scale(x$annualtot),col="blue")
		   lines(scale(x$sixmonthtot),col="brown")
	   }

           smoothfunctionsplot(tempdaily,precipdaily,id=100,xlim=c(1,17532),lines=F)

### Generate Seasonal Data
              TAAve=as.data.frame(tapply(tempdaily$annualave,list(tempdaily$tid,tempdaily$sid), tempave))
                write.table(TAAve,datadir("TAAve.csv"),row.names=F,col.names=F,sep=",")
              T6Ave=as.data.frame(tapply(tempdaily$sixmonthave,list(tempdaily$tid,tempdaily$sid), tempave))
                write.table(T6Ave,datadir("T6Ave.csv"),row.names=F,col.names=F,sep=",")
              TSAve=as.data.frame(tapply(tempdaily$tave,list(tempdaily$tid,tempdaily$sid),tempave))
                write.table(TSAve,datadir("TSAve.csv"),row.names=F,col.names=F,sep=",")
              TSHotweek=as.data.frame(tapply(tempdaily$tave,list(tempdaily$tid,tempdaily$sid),max, na.rm=T))
                write.table(TSHotweek,datadir("TSHotweek.csv"),row.names=F,col.names=F,sep=",")          

      #Precipitation
#          precipdaily=read.csv(datadir("precipdaily.csv"))
          precipdaily$ppt=ifelse(precipdaily$ppt<0,0,precipdaily$ppt)       #set negative (?) precip to zero.
        #Annual
          PATot= as.data.frame(tapply(precipdaily$annualtot,list(precipdaily$pid,precipdaily$sid),tempave)) #this takes the average of the annual sums for the season
                write.table(PATot,datadir("PATot.csv"),row.names=F,col.names=F,sep=",")
          P6Tot= as.data.frame(tapply(precipdaily$sixmonthtot,list(precipdaily$pid,precipdaily$sid),tempave)) #this takes the average of the annual sums for the season
                write.table(P6Tot,datadir("P6Tot.csv"),row.names=F,col.names=F,sep=",")
					#Seasonal
          PSTot= as.data.frame(tapply(precipdaily$precip,list(precipdaily$pid,precipdaily$sid),preciptot))
                write.table(PSTot,datadir("PSTot.csv"),row.names=F,col.names=F,sep=",")


#Other time varying variables
  #Climate Indicies
      #import data
      season=read.csv(datadir("season.csv")); season$date=as.Date(season$date)
      aao=climreshape("/media/Data/Work/Climate/indices/AAO_Marshall.txt","aao"); aao=aao[aao$date>=start&aao$date<=stop,]
      aao_old=climreshape("/media/Data/Work/Climate/indices/AAO.txt","aao_old"); aao_old=aao_old[aao_old$date>=start&aao_old$date<=stop,]
      soi=climreshape("/media/Data/Work/Climate/indices/soi.txt","soi"); soi=soi[soi$date>=start&soi$date<=stop,]
      best=climreshape("/media/Data/Work/Climate/indices/BESTenso.txt","best"); best=best[best$date>=start&best$date<=stop,]
      nino3=climreshape("/media/Data/Work/Climate/indices/nina3.txt","nino3"); nino3=nino3[nino3$date>=start&nino3$date<=stop,]
      solar=climreshape("/media/Data/Work/Climate/indices/SolarFlux.txt","solar"); solar=solar[solar$date>=start&solar$date<=stop,]
      tad=climreshape("/media/Data/Work/Climate/indices/TAD.txt","tad"); tad=tad[tad$date>=start&tad$date<=stop,]
      #merge
      indexdata=merge(season,subset(aao,select=c("aao","date")),by="date",all.x=T)
      indexdata=merge(indexdata,subset(aao_old,select=c("aao_old","date")),by="date",all.x=T)
      indexdata=merge(indexdata,subset(soi,select=c("soi","date")),by="date",all.x=T)
      indexdata=merge(indexdata,subset(best,select=c("best","date")),by="date",all.x=T)
      indexdata=merge(indexdata,subset(nino3,select=c("nino3","date")),by="date",all.x=T)
      indexdata=merge(indexdata,subset(tad,select=c("tad","date")),by="date",all.x=T)

					#write to disk
      write.csv(indexdata,datadir("indexdata.csv"),row.names=F)

#close.nc(precip);close.nc(tmax);close.nc(tmin)

####################################################################################################
#Time Since Fire
#write fires table to cfr
      season=read.csv(datadir("season.csv")); season$date=as.Date(season$date)
      static=read.csv(datadir("static.csv"))

fire=importShapefile("/media/Data/Work/Regional/ZA/CFR/FireAnalysis/FireData/RevisedFires.shp",projection="LL")
      grid=importShapefile(datadir("grid.shp"),projection="LL")

### Intersect grid and vegetation 
     intersectAreas <- function(polygon, grid) {
	     gridid=expand.grid(gid=unique(grid$PID),fireid=unique(fire$PID))  #identify all combinations to create new PIDs
	     gridid$joinid=1:nrow(gridid) #set up ids to result from intersect
	     firegrid=joinPolys(grid,fire,"INT",maxVert=1e+8)  #intersect grid and fire
#	     areas=merge(calcArea(firegrid,rollup=2),calcCentroid(firegrid,rollup=1),by="PID")
	     areas=calcArea(firegrid,rollup=2)
	     areas=merge(areas,gridid,by.x="PID",by.y="joinid",all.x=T) #merge intersected data and ids to regain correct ids
	     firegrid=merge(areas, attr(grid,"PolyData"),by.x="gid",by.y="PID",all.x=T) #add original grid area attributes
	     firegrid$percburn=firegrid$area/firegrid$area.y
	     firegrid=merge(firegrid,attr(fire,"PolyData"),by.x="fireid",by.y="PID",all.x=T) #add fire attributes
	     firegrid=subset(firegrid,select=c("gid","YEAR","Month","X","Y","percburn"))
	     firegrid$Month=ifelse(firegrid$Month==0,1,firegrid$Month) #switch fires with "0" month to Jan (summer)
	     firegrid$date=as.Date(paste(firegrid$YEAR,firegrid$Month,15,sep="/"),"%Y/%m/%d")
	     firegrid=merge(season,subset(firegrid,select=c("gid","percburn","date")),by="date")
	     firegrid=as.data.frame(t(tapply(firegrid$percburn,list(firegrid$sid,firegrid$gid),function(x) min(1,sum(x))))) #sum by cell (remove values >1)
	     firegrid = reshape(firegrid,varying=list(colnames(firegrid)),ids=as.numeric(rownames(firegrid)),times=as.numeric(colnames(firegrid)),direction="long",v.names="fire",timevar="sid",idvar="gid")
	     return(firegrid)
     }
firegrid=intersectAreas(fire,grid); gc()
					# create grid table with all possible grid/month combinations        
      gids=unique(static$gid)
      sids=unique(season$sid)
      grid=expand.grid(gids,sids); colnames(grid)=c("gid","sid")
      firegrid=merge(firegrid,grid,by=c("gid","sid"),all.y=T)
      firegrid$fire[is.na(firegrid$fire)]=0                                       #convert NAs to 0s
#extract fire occurence table    
    firegrid$fire=ifelse(firegrid$fire<0.25,0,1)                      #fires with less than x (0.25) burned are called 0s
    write.csv(firegrid,datadir("fire.csv"),row.names=F)

#Test and compare to original fire data
x=merge(firegrid,subset(static,select=c("gid","lon","lat")),by="gid")
x=x[x$fire==1,]
x=merge(x,season,by="sid")
x=x[x$sid>112,]
par(mfrow=c(7,3),mai=c(0.1,0.1,0.2,0.2));for(i in sort(unique(x$fireyear))) plot(x$lat[x$fireyear==i]~x$lon[x$fireyear==i],pch=15,cex=.5,main=i,ylim=c(-34.5,-32),xlim=c(18.5,23.5),col="red")

fire=importShapefile("/media/Data/Work/Regional/ZA/CFR/FireAnalysis/FireData/RevisedFires.shp",projection="LL")
attr(fire,"PolyData")$date=as.Date(paste(attr(fire,"PolyData")$YEAR,attr(fire,"PolyData")$Month,15,sep="/"),"%Y/%m/%d")
attr(fire,"PolyData")=merge(attr(fire,"PolyData"),season,by="date")
attr(fire,"PolyData")=attr(fire,"PolyData")[sort(attr(fire,"PolyData")$PID),]
par(mfrow=c(7,3),mai=c(.1,.1,0.2,0.2));for(i in unique(x$fireyear)){
	index=attr(fire,"PolyData")[attr(fire,"PolyData")$fireyear==i,"PID"]
	plotPolys(fire[fire$PID%in%index,],main=i,col="red",ylim=c(-34.5,-32),xlim=c(18.5,23.5))
}


#####
# Fill in the observed time since fires
fire=firegrid
firew=reshape(fire, timevar = "sid", idvar = "gid",direction="wide")[,-1]; colnames(firew)=1:length(sids)      #reshape to wide format

gid=sort(unique(fire$gid))
sid=sort(unique(fire$sid))
tsfw=matrix(nrow=length(gid),ncol=length(sid));colnames(tsfw)=1:length(sid)
event=data.frame()

for (g in gid){
    fired=as.numeric(firew[g,])
#    event=event[order(event$sid),]
    event=as.data.frame(cbind(fired,time=NA))
    for (s in 2:max(sid))
        event$time[s]= ifelse(event$fired[s-1]>0,1,event$time[s-1]+1)
        tsfw[g,]=event[,2] 
    print(paste(g,"out of",length(gid),sep=" "))
    }
    tsf=reshape(as.data.frame(tsfw), varying=list(colnames(tsfw)),ids=1:nrow(tsfw),times=1:ncol(tsfw),v.names="tsf",timevar = "sid", idvar = "gid",direction="long")
fire=    merge(fire,tsf,by=c("sid","gid"),all=T);fire=fire[order(fire$gid,fire$sid),]
write.csv(fire,datadir("fire.csv"),row.names=F)

################################################################################
################################################################################
# Bring it all together
tolong <- function(x,name,g=1:nrow(x),s=1:ncol(x)){
	x = x[g,s]
	x = as.data.frame(x)
	rownames(x) = 1:nrow(x)
	x = reshape(x,varying=list(names(x)),ids=rownames(x),times=1:ncol(x),direction="long",v.names=name,timevar="sid",idvar="gid")
	return(x)
}
season=read.csv(datadir("season.csv"))

yind=as.numeric(as.factor(unique(cbind(season$sid,season$fireyear))[,2]))
tind=static$tid
pind=static$pid
sind=as.numeric(as.factor(unique(cbind(season$sid,season$season))[,2]))
nGrids=nrow(static)
nEco=length(unique(static$regionid))

#Import data matrices
TAAveold=read.csv(datadir("tempannualaveold.csv"),header=T)
  TAAveold=scalem(TAAveold)
  TAAveold=tolong(TAAveold,"TAAveold",tind,yind)
TAAve=read.csv(datadir("TAAve.csv"),header=F)
  TAAve=scalem(TAAve)
  TAAve=tolong(TAAve,"TAAve",tind,sind)
T6Ave=read.csv(datadir("T6Ave.csv"),header=F)
  T6Ave=scalem(T6Ave)
  T6Ave=tolong(T6Ave,"T6Ave",tind,sind)
TSHotweek=read.csv(datadir("TSHotweek.csv"),header=F)
  TSHotweek=t(seasscale(TSHotweek,sind))
#  TSHotweek=scalem(TSHotweek)
  TSHotweek=tolong(TSHotweek,"TSHotweek",tind)
PATotold=read.csv(datadir("PATotold.csv"),header=F)
  PATotold=scalem(PATotold)
  PATotold=tolong(PATotold,"PATotold",pind,yind)
PATot=read.csv(datadir("PATot.csv"),header=F)
  PATot=scalem(PATot)
  PATot=tolong(PATot,"PATot",pind,sind)
P6Tot=read.csv(datadir("P6Tot.csv"),header=F)
  P6Tot=scalem(P6Tot)
  P6Tot=tolong(P6Tot,"P6Tot",pind,sind)
fire=read.csv(datadir("fire.csv"))
PSTot=read.csv(datadir("PSTot.csv"),header=F)
  PSTot=t(seasscale(PSTot,sind))
#  PSTot=scalem(PSTot)
  PSTot=tolong(PSTot,"PSTot",pind)
TSAve=read.csv(datadir("TSAve.csv"),header=F)
  TSAve=t(seasscale(TSAve,sind))
#  TSAve=scalem(TSAve)
  TSAve=tolong(TSAve,"TSAve",tind)
#lon=tolong(matrix(rep(lon,nSeasons),nrow=nGrids),"lon")
nYears=ncol(TAAve)
nSeasons=ncol(TSAve)

summary(unlist(TSAve))
summary(unlist(TSHotweek))
summary(unlist(TAAve))
summary(unlist(PSTot))
summary(unlist(PATot))

list=c("TAAve","TSAve","TSHotweek","PATot","PSTot") #,"MADroughtdays","MSDroughtdays","MADroughtmax","MSDroughtmax")
data=get(list[1])
for( i in 2:length(list)) {
	data=merge(data,get(list[i]),by=c("gid","sid"))
	print(paste(i," in ",length(list)))
}


#attach static items from static table to the monster
data=merge(cbind(gid=static$gid,tid=static$tid,pid=static$pid,eco=static$regionid),data,by="gid",all.y=T)   #rain=as.numeric(scale(static$rain))

#temporal data
  seastable=unique(cbind(season$season,season$sid,as.numeric(as.factor(season$fireyear))));colnames(seastable)=c("season","sid","yind")
  data=merge(seastable,data,by="sid",all.y=T)             #add season and yind to data
  indexdata=read.csv(datadir("indexdata.csv"))  # climate indices
  seastable=merge(seastable,subset(indexdata,select=c("sid","aao","soi","best","nino3","tad")),by="sid")
  nSeason=length(unique(season$sid))
  seasdata=data.frame(sid=1:nSeason)
  seasdata$aao=as.vector(tapply(seastable$aao,list(seastable$sid),mean))
  seasdata$soi=as.vector(tapply(seastable$soi,list(seastable$sid),mean))
  seasdata$best=as.vector(tapply(seastable$best,list(seastable$sid),mean))
  seasdata$nino3=as.vector(tapply(seastable$nino3,list(seastable$sid),mean))
  seasdata$tad=as.vector(tapply(seastable$tad,list(seastable$sid),mean))

rm(seastable)
  data=merge(seasdata,data,by="sid",all.y=T)             #add season id

#fire data
  fire= read.csv(datadir("fire.csv"))
  data=merge(fire,data,by=c("gid","sid"))
#Order and write it
  data=data[order(data$gid,data$sid),]
  write.csv(data,datadir("data.csv"),row.names=F)


       Fireind =        as.numeric(rownames(fire)[apply(fire[,-nSeasons],1,max)==1]) 
       Nfiregrids =     length(Fireind)
       leftlim=         vector(); for(i in 1:Nfiregrids) leftlim[i]=min((tsfleftobs[Fireind[i]]+1),NSeasons)
       sind =           as.numeric(season[,1])
       east =           ifelse(static$longitude>20.5,1,0)
       lon =            static$longitude


data2=data[data$sid>112,]
glm=glm(fire~aao+as.factor(season)+TAAve+TSAve+TSHotweek+PATot+PSTot,data=data2,family=binomial(link="probit")); summary(glm)

coef=cbind(glm$coefficients,confint(glm))
errbar(rownames(coef)[-1],coef[-1,1],yplus=coef[-1,2],yminus=coef[-1,3],main="Regression Coefficients   (+/- 95% CI)",cex.lab=2,cex.main=2,cex.axis=2);abline(v=0,lty="dashed")
text(.5,10,"Seasonal Anomolies")
eps(file=datadir("plots_seasonalanoms.eps"))

