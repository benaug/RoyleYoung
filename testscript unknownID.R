#"appears to work". Will remove this when tested.

library(nimble) #data simulator uses nimble
library(truncnorm) #required for data simulator
source("sim.RY.R")
source("Nimble Functions unknownID.R") #nimble functions used in data simulator
#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

#state space. Must start at (0,0)
xlim <- c(0,100)
ylim <- c(0,100)
res <- 2.5 #resolution, cell width/height
if(xlim[1]!=0|ylim[1]!=0)stop("xlim and ylim must start at 0.")
if((diff(range(xlim))/res)%%1!=0)stop("The range of xlim must be divisible by 'res'")
if((diff(range(ylim))/res)%%1!=0)stop("The range of ylim must be divisible by 'res'")

#make discrete state space objects
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res)
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res)
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov, higher cov.pars for large scale cov
set.seed(1320562)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(25,25),messages=FALSE)[[2]]
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#simulate an rsf cov, lower cov.pars for finer scale cov
set.seed(24674345)
rsf.cov <- grf(n.cells,grid=dSS,cov.pars=c(5,5),messages=FALSE)[[2]]
rsf.cov <- as.numeric(scale(rsf.cov)) #scale
image(x.vals,y.vals,matrix(rsf.cov,n.cells.x,n.cells.y),main="rsf.cov",xlab="X",ylab="Y",col=cols1)

#make state space mask - just making a circle
dists <- sqrt((dSS[,1]-mean(dSS[,1]))^2+(dSS[,2]-mean(dSS[,2]))^2)
rad <- min(diff(xlim),diff(ylim))/2
InSS <- 1*(dists<rad)
image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

#simulate effort - I admit, I am being lazy here, but this gives us
#effort that varies over space with replication. You need some minimal
#level of spatial replication in survey effort
K <- 12
#Let's implement a state space buffer and only sample inside that.
#buffer should be at least 3sigma (which we haven't defined, yet, I use 7.8 below)
target.sigma <- 2
search.buff <- 1*(dists<(rad-3*target.sigma))

#but let's consider we don't survey the entire state space on each occasion
#let's say there are 4 sections and we do 1 per occasion
split.x <- x.vals[round(length(x.vals)/2)]
split.y <- y.vals[round(length(y.vals)/2)]
idx.list <- vector("list",K)
idx.list[[1]] <- which(dSS[,1]<split.x&dSS[,2]<split.y&search.buff)
idx.list[[2]] <- which(dSS[,1]<split.x&dSS[,2]>=split.y&search.buff)
idx.list[[3]] <- which(dSS[,1]>=split.x&dSS[,2]>=split.y&search.buff)
idx.list[[4]] <- which(dSS[,1]>=split.x&dSS[,2]<split.y&search.buff)
#must have replication, using same survey cells as above
idx.list[[5]] <- idx.list[[1]]
idx.list[[6]] <- idx.list[[2]]
idx.list[[7]] <- idx.list[[3]]
idx.list[[8]] <- idx.list[[4]]
#adding 4 more
idx.list[[9]] <- idx.list[[1]]
idx.list[[10]] <- idx.list[[2]]
idx.list[[11]] <- idx.list[[3]]
idx.list[[12]] <- idx.list[[4]]

set.seed(3356735)
effort <- survey <- matrix(0,n.cells,K)
for(k in 1:K){
  tmp <- grf(n.cells,grid=dSS,cov.pars=c(15,15),nugget=1,messages=FALSE)[[2]]
  effort[idx.list[[k]],k] <- tmp[idx.list[[k]]]
  survey[idx.list[[k]],k] <- 1
}

#visualize effort summed over all occasions
image(x.vals,y.vals,matrix(rowSums(effort),n.cells.x,n.cells.y),main="Search Effort",xlab="X",ylab="Y",col=cols2)

# visualize effort on each occasion, scroll through k,
# can see different areas surveyed on each occasions 1-4 and 5-8
# k <- 1
# k <- k+1
# image(x.vals,y.vals,matrix(effort[,k],n.cells.x,n.cells.y),main="Search Effort",col=cols2)

#scale effort, but only in surveyed cell/occasions
effort[survey==1] <- as.numeric(scale(effort[survey==1]))


D.beta0 <- -4.5 #baseline D
D.beta1 <- 1.0 #density coefficient 
rsf.beta <- 1.5 #rsf coefficient
beta.p.int <- -0.5 #baseline detection prob
beta.p.effort <- 2.5 #effort effect on detection prob
sigma <- 2 #spatial scale of availability distribution
n.tel.inds <- 10 #number of telemetry individuals
K.tel <- 15 #number of telemetry locations per individual

#visualize distribution of p across cells given search effort and detection parameters
p.test <- plogis(beta.p.int + beta.p.effort*effort)
par(mfrow=c(1,1),ask=FALSE)
hist(p.test[survey==1])

# set.seed(32348)
data <- sim.RY(D.beta0=D.beta0,D.beta1=D.beta1,rsf.beta=rsf.beta,
               sigma=sigma,beta.p.int=beta.p.int,beta.p.effort=beta.p.effort,
               xlim=xlim,ylim=ylim,res=res,InSS=InSS,D.cov=D.cov,rsf.cov=rsf.cov,
               effort=effort,survey=survey,n.tel.inds=n.tel.inds,K=K,K.tel=K.tel)

data$truth$lambda #expected abundance from D cov inputs
data$truth$N #simulated realized abundance
data$truth$n #number of inds captured
table(rowSums(data$capture$y)) #number of inds captures X times

#given search effort and use distributions, here are individuals cumulative detection probabilities
#spatial het in space use and effort lead to lots of het in p
#this is relevant to doing nonspatial capture-recapture
hist(data$summaries$p.marg.i,breaks=50,main="Individual Cumulative Detection Probabilities|Search Effort",
     xlim=c(0,1),xlab="Detection Prob")

#can inspect every individual's availability and use distributions
#useful to check if simulated behavior is realistic
# par(mfrow=c(2,1)) #if you want to plot availability over use
# i <- 1
# i <- i + 1
# image(x.vals,y.vals,matrix(data$truth$avail.dist[i,],n.cells.x,n.cells.y),main="Availibility Distribution",xlab="X",ylab="Y")
# image(x.vals,y.vals,matrix(data$truth$use.dist[i,],n.cells.x,n.cells.y),main="Use Distribution",xlab="X",ylab="Y")
# par(mfrow=c(1,1)) #set back

##Fit Model
library(nimble)
library(coda)
nimbleOptions(determinePredictiveNodesInModel = FALSE)
source("init.data.RY.unknownID.R")
source("NimbleModel unknownID.R")
source("sSampler Dcov RSF Marginal.R")
M <- 200 #data augmentation limit. Must be larger than simulated N. If N posterior hits M, need to raise M and try again.
if(M<=data$truth$N)stop("Raise M to be larger than simulated N.")

inits <- list(sigma=1) #needs to be set somewhere in the ballpark of truth
nimbuild <- init.data.RY.unknownID(data=data,inits=inits,M=M)
n.samples <- nimbuild$n.samples

#plot initialized data
plot(nimbuild$s,pch=16,xlim=data$xlim,ylim=data$ylim,col="grey")
points(nimbuild$s[nimbuild$z==1,1],nimbuild$s[nimbuild$z==1,2],pch=16,cex=1.25)
for(l in 1:n.samples){
  points(data$u.obs[l,1],data$u.obs[l,2],pch=16,col="lightblue",cex=0.75)
  lines(x=c(data$u.obs[l,1],nimbuild$s[nimbuild$ID[l],1]),
        y=c(data$u.obs[l,2],nimbuild$s[nimbuild$ID[l],2]))
}


n.surveyed.cells <- colSums(survey)
max.surveyed.cells <- max(n.surveyed.cells)
surveyed.cells <- matrix(0,max.surveyed.cells,K)
surveyed.cells.effort <- matrix(0,max.surveyed.cells,K)
for(k in 1:K){
  surveyed.cells[1:n.surveyed.cells[k],k] <- which(survey[,k]==1)
  surveyed.cells.effort[1:n.surveyed.cells[k],k] <- effort[surveyed.cells[1:n.surveyed.cells[k],k],k]
}

#GPS locs
u.tel <- data$telemetry$u.tel
n.tel.inds <- nrow(u.tel)
n.locs.ind <- rowSums(!is.na(u.tel[,,1]))
n.locs.max <- max(n.locs.ind)
xlim.GPS <- range(u.tel[,,1],na.rm=TRUE) + c(-1,1) #add small arbitrary buffer for edge points
ylim.GPS <- range(u.tel[,,2],na.rm=TRUE) + c(-1,1)

points(u.tel[,,1],u.tel[,,2],pch=".",col="darkred",cex=3)

#Need u.cell.tel data for rsf
u.cell.tel <- matrix(NA,n.tel.inds,n.locs.max)
for(i in 1:n.tel.inds){
  for(k in 1:n.locs.ind[i]){
    u.cell.tel[i,k] <- getCellR(u.tel[i,k,],data$constants$res,data$constants$cells,xlim.GPS,ylim.GPS)
  }
}
any(u.cell.tel==0,na.rm=TRUE) #should be false

#these are observed, derived from u.cell.tel
dSS <- data$constants$dSS
u.xlim.tel <- u.ylim.tel <- array(0,dim=c(n.tel.inds,n.locs.max,2))
for(i in 1:n.tel.inds){
  for(k in 1:n.locs.ind[i]){
    u.xlim.tel[i,k,1] <- dSS[u.cell.tel[i,k],1]-res/2
    u.xlim.tel[i,k,2] <- dSS[u.cell.tel[i,k],1]+res/2
    u.ylim.tel[i,k,1] <- dSS[u.cell.tel[i,k],2]-res/2
    u.ylim.tel[i,k,2] <- dSS[u.cell.tel[i,k],2]+res/2
  }
}

Niminits <- list(z=nimbuild$z,s=nimbuild$s,
                 ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true),
                 y.true=nimbuild$y.true,N=nimbuild$N,
                 u=nimbuild$u,u.cell=nimbuild$u.cell,
                 D.beta0=-4,D.beta1=0,
                 sigma=inits$sigma,
                 beta.p.int=c(0),beta.p.effort=c(2),
                 rsf.beta=2,
                 s.tel=apply(data$telemetry$u.tel,c(1,3),mean,na.rm=TRUE))

#constants for Nimble
constants <- list(M=M,K=K,n.samples=n.samples,
                  n.cells=data$constants$n.cells,n.cells.x=data$constants$n.cells.x,
                  n.cells.y=data$constants$n.cells.y,res=data$constants$res,
                  x.vals=data$constants$x.vals,y.vals=data$constants$y.vals,
                  xlim=data$constants$xlim,ylim=data$constants$ylim,
                  n.locs.ind=data$constants$n.locs.ind,n.tel.inds=data$constants$n.tel.inds,
                  D.cov=D.cov,
                  cellArea=data$constants$res^2,surveyed.cells=surveyed.cells,n.surveyed.cells=n.surveyed.cells,
                  surveyed.cells.effort=surveyed.cells.effort,
                  u.xlim.tel=data$telemetry$u.xlim.tel,
                  u.ylim.tel=data$telemetry$u.ylim.tel)

#supply data to nimble
Nimdata <- list(y.true=matrix(NA,nrow=M,ncol=K),u=array(NA,dim=c(M,K,2)),
                u.tel=data$telemetry$u.tel,u.cell.tel=data$telemetry$u.cell.tel,
                cells=data$constants$cells,InSS=data$constants$InSS,
                z=nimbuild$z,
                dummy.data=rep(1,M),
                rsf.cov=rsf.cov)

# set parameters to monitor
parameters <- c('beta.p.int','beta.p.effort','rsf.beta','D.beta1',
              'sigma','N','D.beta0','lambda','n')

#can also monitor a different set of parameters with a different thinning rate
nt <- 2 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#tell nimble which nodes to configure so we don't waste time for samplers we will replace below
config.nodes <- c("beta.p.int","sigma",'D.beta0','D.beta1','rsf.beta','beta.p.effort')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = FALSE,
                      nodes=config.nodes) 

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
# conf$removeSampler("y.true")
calcNodes.y.true <- Rmodel$getDependencies("y.true")
calcNodes.u <- Rmodel$getDependencies("u")
calcNodes.all <- c(calcNodes.y.true,calcNodes.u)
map <- t(matrix(1:(M*K),K,M)) #map 2D node reference to 1D node vector calcNodes.y.true
conf$addSampler(target = paste0("y.true[1:",M,",1:",K,"]"),
                type = 'IDSampler',control = list(M=M,K=K,n.samples=nimbuild$n.samples,
                                                  u.obs=data$capture$u.obs2D,map=map,
                                                  this.k=nimbuild$this.k,
                                                  calcNodes.y.true=calcNodes.y.true,
                                                  calcNodes.u=calcNodes.u,
                                                  calcNodes.all=calcNodes.all),
                silent = TRUE)

###*required* sampler replacement for "alternative data augmentation" N/z update
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
# conf$removeSampler("N")
#nodes used for update, calcNodes + z nodes
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",K,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,K=K,M=M,
                                                 y.nodes=y.nodes,N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),
                silent = TRUE)

# replace default activity center sampler
# conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  calcNodes <- Rmodel$getDependencies(paste("s[",i,",1:2]"))
  conf$addSampler(target = paste("s[",i,",1:2]", sep=""),
                  type = 'sSamplerDcovRSF',control = list(i=i,K=K,xlim=data$constants$xlim,
                                                          ylim=data$constants$ylim,
                                                          n.cells.x=data$constants$n.cells.x,
                                                          n.cells.y=data$constants$n.cells.y,
                                                          res=data$constants$res,
                                                          calcNodes=calcNodes), silent = TRUE)
}

# more efficient s sampler for telemetry inds
# conf$removeSampler(paste("s.tel[1:",n.tel.inds,", 1:2]", sep=""))
for(i in 1:data$constants$n.tel.inds){
  calcNodes.s.tel <- Rmodel$getDependencies(paste("s.tel[",i,",1:2]"))
  conf$addSampler(target = paste("s.tel[",i,",1:2]", sep=""),
                  type = 'sSamplerDcovRSF.tel',control = list(i=i,xlim=data$constants$xlim,
                                                              ylim=data$constants$ylim,
                                                              calcNodes.s.tel=calcNodes.s.tel), silent = TRUE)
}

conf$addSampler(target = c("D.beta0","D.beta1"),
                type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(1000,reset=FALSE) #can keep running this line to extend sampler
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[10:nrow(mvSamples),])) #discarding some burnin here. Can't plot 1st sample which is all NA

data$truth$lambda #target expected abundance
data$truth$N #target realized abundance
data$truth$n #target number detected (n)

par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(Cmodel$lambda.cell,n.cells.x,n.cells.y),main="Final Iteration Spatial D and ACs")
points(Cmodel$s[Cmodel$z==1,],pch=16,cex=0.75,col="grey30")
#truth
image(x.vals,y.vals,matrix(data$truth$lambda.cell,length(x.vals),length(y.vals)),
      main="Spatially Explicit D and Realized ACs",xlab="X",ylab="Y")
points(data$truth$s,pch=16,cex=0.75,col="grey30")
