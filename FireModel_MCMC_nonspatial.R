
setwd("/media/Data/Work/Regional/ZA/CFR/FireAnalysis/ModelData_2km/")

# Location of files
datadir <- function(x) paste("/media/Data/Work/Regional/ZA/CFR/FireAnalysis/ModelData_2km/",x,sep="")
statdir <- function(x) paste("/media/Data/Work/Stat/Scripts/FireModel/",x,sep="")

                                        #Run functions
source("/media/Data/Work/Stat/Scripts/Functions.r")
loadnecessaries <- function() {
    library(Umacs); library(snow); library(rpvm) ; 
    library(MCMCpack); library(coda); library(rv); library(MASS); library("geneplotter")
    source(statdir("FireModel_Functions.R"))
}; loadnecessaries(); library(grid); library(gridBase); library(Hmisc); library(maptools)


## Start PVM and spawn slaves
.PVM.start.pvmd()       # Starts the background PVM service
cl=makeCluster(3)       # Starts "cluster" cl with 3 nodes  (it is possible to do more nodes, but with three chains only three are needed)
# .PVM.halt()            # closes all pvm activities and kills all R processes!    file.remove("/tmp/pvmd.500")


### load the data
### NOTE NOW WE NEED ALL CLIMATE DATA FROM 1950 - 1999, so times are 1 through 200 or 201
static=read.csv(datadir("static.csv"))
static=static[order(static$gid),]
nGrid=nrow(static)
data=read.csv(datadir("data.csv"))
data=data[order(data$gid,data$sid),] #order it

### Extract precipitation seasonality from static table and merge it with the data
conc=read.csv(datadir("conc.csv"))[,-1]
static=merge(static,conc,by="pid")
static=static[order(static$gid),]
data=merge(data,cbind(gid=static$gid,pptconc=static$pptconc),by="gid")
data$pptconc=scale(data$pptconc)
data=data[order(data$gid,data$sid),] #order it

eco = data$regionid
nEco = length(unique(eco))
season = data$season
sid=data$sid
N = length(season)
nSeas = N/nGrid
gid = data$gid
gidunique = unique(gid)
gidrow = apply(t(gidunique), 1, match, gid) # vector indexing the first row of data for each gid
fire = data$fire
select=c("TAAve","TSAve","TSHotweek","PATot","PSTot","aao","pptconc") 
datamatrix = as.matrix(cbind(intercept=1,subset(data,select=select)))
datamatrix[,"aao"][is.na(datamatrix[,"aao"])]=mean(datamatrix[,"aao"],na.rm=T)

### sensor lengths
     leftcensorlength = 112 # how many seasons before our fire record starts (i.e. between 1950 and 1979)
     leftcensorindex = rep(F, N) # make a full-length index of which data to exclude from fitting Beta and random effects because outside of obs period
     for (i in 1:nGrid) leftcensorindex[gidrow[i]:(gidrow[i]+leftcensorlength-1)] = T
     fire[data$sid<=leftcensorlength]=0  # get rid of pre-1980 fires so the lastnoburn index is correct
# add indicators for seasons
     spring=ifelse(season==3,1,0)
     summer=ifelse(season==4,1,0)
     fall=ifelse(season==1,1,0)
     datamatrix=cbind(datamatrix,spring,summer,fall)
     nBeta=ncol(datamatrix)
# Make index for times and locations of observed fires
     firerows = grep(1, fire)
     fireTF = fire==1
     lastnoburn = rep(NA, nGrid) # vector to hold last interval for each cell before the first fire was observed
     for (i in 1:nGrid) {
	     temp = data$gid==gidunique[i] # select rows of data corresponding to grid cell i
	     lastnoburn[i] = match(1, fire[temp]) - 1
     }
     lastnoburn[is.na(lastnoburn)] = nSeas


##Build spatial Data - build adjancency matrix
     adjm=AdjMat(as.data.frame(cbind(gid=static$gid,y=static$lat,x=static$lon)), 0.021)
     adj=adjm$adj;num=adjm$num;SumNumNeigh=adjm$SumNumNeigh

rm(data)  # get rid of the original matrix to save space

################
# Initial Values
 beta={}
 log={}
setup <- function(iterations){
    ## Data Structures
    nBeta <<- ncol(datamatrix)
    ## Intital values
    beta.init=function() c(3.45,.18,-.3,-.05,-.15,.03,-.05,-.55,-.2,-1.1,-.8) #rnorm(nBeta,0,0.5)
    ecoRF.init=function() rnorm(nEco,0,0.5)
    fire0.init=function() sample(1:40, nGrid, replace=T)
#    propsd_W <<- 0.005
#    propsd_Beta <<- 0.005
    propsd_Beta <<- c(0.018,rep(0.005,7),rep(0.018,3))

    
### Data Structure
    beta <<- matrix(ncol=nBeta,nrow=iterations[length(iterations)]+1) ; beta[1,] <<- beta.init()
    fire0 <<- matrix(ncol=nGrid,nrow=iterations[length(iterations)]+1) ; fire0[1,] <<- fire0.init()
    PLF <<- matrix(ncol=2,nrow=iterations[length(iterations)]+1);colnames(PLF)=c("Gm","Pm");PLF[1,]=c(0,0)
    Dbar <<- vector(length=iterations[length(iterations)]+1); Dbar[1]=0
    sigma.beta <<- vector(length=iterations[length(iterations)]+1); sigma.beta[1]=0.1

}

sampler <- function(iterations) {
    ## MCMC Loop
    for(i in iterations) {
	XBeta = datamatrix%*%beta[i,]
        fire0[i+1,] <<- fire0.update(pnorm(XBeta), lastnoburn, nGrid, gidrow, leftcensorlength) # get times of unobserved initial fires        
        leftexclude = UpdateLeftExclude(fire0[i,], nGrid, nSeas)
        params.nonspatial = list(Y=fire,X=datamatrix,Beta=beta[i,], leftexclude=leftexclude, W=rep(0,nGrid), gid=gid, firerows=firerows)
        beta[i+1,] <<- metrop.norm.subset(beta[i,], propsd_Beta, CalcPostNonspatial, params.nonspatial, "Beta",which=list(c(1,9:11),c(2:8)))
	Dbar[i] <<- -2*CalcLogLik(Y=fire, X=datamatrix, Beta=beta[i,], leftexclude=leftexclude, W=rep(0,nGrid), gid=gid, firerows=firerows)
        print(i)
    }
}

setup2 <- function(iterations){
    beta <<- rbind(beta,matrix(ncol=nBeta,nrow=length(iterations)))
    fire0 <<- rbind(fire0,matrix(ncol=nGrid,nrow=length(iterations)))
}
setup(1:5) #set up for 5 iterations
system.time(sampler(1:5))[3]/5  #this runs 5 iterations and reports the time/iteration

setup2(5:15)  #this sets it up for x iterations, you have to start where setup() ended.
sampler(5:15) #this runs it for x iterations, reporting each iteration.  

## collate and plot
par(mfrow=c(2,5),mai=c(.2,.3,.2,.2));traceplot(as.mcmc(beta),smooth=F,ask=T)
par(mfrow=c(2,5),mai=c(.2,.3,.2,.2));traceplot(as.mcmc(W[,1:10]),smooth=F,ask=T)

## Send to other processes
clusterEvalQ(cl,ls())
t=clusterEvalQ(cl, library(MASS))  #Set up the sampler on each node
t=clusterEvalQ(cl, library(MCMCpack))  #Set up the sampler on each node
t=clusterExport(cl,ls())   #Export data and Sampler function to nodes
t=clusterEvalQ(cl, setup(1:5))  #Set up the sampler on each node
t=clusterEvalQ(cl, system.time(sampler(1:5)));t  #Test sampler on each node

clusterExport(cl,c("sampler","propsd_Beta","CalcPost","CalcPostNonspatial" ))
clusterEvalQ(cl, gc())  #Set up the sampler on each node
paste(round(12000*max(unlist(t))/5/3600,1)," Hours for this run",sep="");
t=clusterEvalQ(cl, setup2(6:12000))  #Setup data structure on each node
t=clusterEvalQ(cl, system.time(sampler(2001:12000)));t  #Run the sampler on each node

## restart from previous run



## set burnin
r = c(1,2000)
r2 = c(2001,2000)

#pdf(outputdir("ModelOutput.pdf"),width=11,height=8,paper="USr",family="Times")

############
## Convert to MCMCpack format
# Bring the data back
outputdir <- function(x) paste("/media/Data/Work/Regional/ZA/CFR/FireAnalysis/ModelData_2km/Results_NoSpaceFinal/",x,sep="")
clusterExport(cl,c("r","r2","outputdir"))   #Export data and Sampler function to nodes

log=clusterEvalQ(cl,log)                      #bring back to master
#Z=clusterEvalQ(cl,Z)                      #bring back to master


savedate = "14june08" # label for chain file -- gives the date of run

#betas
beta=as.mcmc.list(clusterEvalQ(cl,as.mcmc(beta)))                      #bring back to master
varnames(beta)=c(colnames(datamatrix))   #Names
save(beta,file=outputdir(paste("beta",savedate, ".R", sep="")))
write(beta[[1]], "betachain1.txt")
write(beta[[2]], "betachain2.txt")
write(beta[[3]], "betachain3.txt")

fire0=as.mcmc.list(clusterEvalQ(cl,window(as.mcmc(fire0),start=r2[1],end=r2[2])))                      #bring back to master
x=save(fire0,file=outputdir(paste("fire0",savedate, ".R", sep="")))
rm(fire0); gc()
write.table(fire0[[1]], "fire0chain1.txt", row.names=F, col.names=F)
write.table(fire0[[2]], "fire0chain2.txt", row.names=F, col.names=F)
write.table(fire0[[3]], "fire0chain3.txt", row.names=F, col.names=F)



Dbar=as.mcmc.list(clusterEvalQ(cl,as.mcmc(Dbar)))
save(Dbar,file=outputdir(paste("Dbar",savedate, ".R", sep="")))
write(Dbar[[1]], "Dbarchain1.txt")
write(Dbar[[2]], "Dbarchain2.txt")
write(Dbar[[3]], "Dbarchain3.txt")



x=clusterEvalQ(cl,gc()) #clean up

########################################################
#  Analysis
## Betas
load(outputdir("beta14june08.R"))
acceptrate(window(beta,r[1],r[2]))

par(mfrow=c(3,4),mai=c(.2,.3,.2,.2));traceplot(beta,smooth=F,ask=T,xlim=r)  #plot the betas
eps(outputdir(paste("betatrace", savedate, ".eps", sep="")))
gelman.diag(window(beta,r2[1],r2[2]),autoburnin=F)

betastats=as.data.frame(summary(window(beta,r2[1],r2[2]),quantiles=c(0.025,.5,0.975))[2])
betastats$names=c("intercept","Annual Average Temperature","Seasonal Average Temperature","Seasonal Hottest Week","Total Annual Precipitation","Total Seasonal Precipitation","AAO","Precipitation Concentration","Spring","Summer","Fall") #,"Precipitation Concentration",
write.csv(betastats,outputdir(paste("betastats", savedate, ".csv", sep="")),row.names=T)

betastats2=betastats[order(betastats$quantiles.50.),] #order it

# Create Seasonal Table
seasonstats=betastats2[c("intercept","spring","summer","fall"),];seasonstats=seasonstats[,-4]
seasonstats$sum=seasonstats$quantiles.50.[1]+seasonstats$quantiles.50.;seasonstats$sum[1]=seasonstats$quantiles.50.[1]
seasonstats$sum025=seasonstats$quantiles.2.5.[1]+seasonstats$quantiles.2.5.;seasonstats$sum025[1]=seasonstats$quantiles.2.5.[1]
seasonstats$sum975=seasonstats$quantiles.97.5.[1]+seasonstats$quantiles.97.5.;seasonstats$sum975[1]=seasonstats$quantiles.97.5.[1]
seasonstats$prob=round(pnorm(seasonstats$sum),3)
seasonstats
write.csv(round(seasonstats,2),outputdir(paste("seasonstats",savedate, ".csv", sep="")))

list=sort(match(c("TAAve","TSAve","TSHotweek","PATot","PSTot","pptconc","aao"),rownames(betastats2)))
par(mfrow=c(1,1),mai=c(1,3,0,.2));errbar(betastats2$names[list],betastats2[list,2],yplus=betastats2[list,1],yminus=betastats2[list,3],main="",xlab="Coefficient Value +/- 95% Credible Intervals",cex.lab=2,cex.main=2,cex.axis=2)
arrows(0,0,0,7.5,lty="dashed",length=0)
mtext("Coefficient Value +/- 95% Credible Intervals",1,padj=4,adj=0.5)
text(-0.125,7.5,"Increasing Fire Risk",cex=1.5)
arrows(-0.08,7.2,-0.17,7.2,length=0.1,lwd=2)
eps(outputdir(paste("betas", savedate,".eps", sep="")))

#gelman.diag(window(beta,r2[1],r2[2]))
clusterEvalQ(cl, gc())  

##Unobserved fires
load(outputdir("fire0.R"))
fire0stats=cbind(quicksummary(fire0),gid=1:nGrid);save(fire0stats,file=outputdir("fire0stats.R"))
fire0stats=merge(cbind(lat=static$lat,lon=static$lon,gid=static$gid),fire0stats,by="gid")
write.csv(fire0stats,outputdir(paste("fire0stats", savedate, ".csv", sep="")))
hist(fire0stats[,6]/4+1950,col="grey",main="mean of unobserved Fires (fire0)",xlab="Year")
rm(fire0);gc()

### Model Comparison
## DIC
load(outputdir("Dbar14june08.R"))
fire0summary=read.csv(outputdir("fire0stats14june08.csv"))
betasstats=read.csv(outputdir("betastats14june08.csv"))
#Wmean=read.csv(outputdir("Wmean.csv"))
leftexclude = UpdateLeftExclude(fire0summary[,"mean"], nrow(fire0summary), nSeas)  #get estimate from mean simulated fire
Dbar=summary(window(Dbar,r2[1],r2[2]))$statistics["Mean"] #mean sampled Dbar
DThetabar= -2*CalcLogLik(Y=fire, X=datamatrix, Beta=betasstats[,"quantiles.50."], leftexclude=leftexclude, W=rep(0,nGrid), gid=gid, firerows=firerows) 
pD=Dbar-DThetabar
DIC=pD+Dbar
DICsummary=rbind("Dbar = "=Dbar,"D(Thetabar) = "=DThetabar,"pD = "=pD,"DIC = "=DIC);DICsummary
write.csv(DICsummary,outputdir(paste("DICsummary", savedate, ".csv", sep="")))


#Write results to disk
save.image(file=outputdir(paste("Rrun2000", savedate,".RData", sep="")))

#Stop the Cluster
stopCluster(cl)
# .PVM.halt()            # closes all pvm activities and kills all R processes!



