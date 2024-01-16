load("grid.objects.RData") #load some objects I'll use to make the grid system
library(nimble) #data simulator uses nimble
library(truncnorm) #required for data simulator
setwd("D:/Sync/Cornell/van Manen/Simulate")
source("sim.RY.R")
source("Nimble Functions.R") #nimble functions used in data simulator

#state space. Must start at (0,0)
xlim <- round(grid.objects$xlim) #rounding here bc upper bounds not exactly integers
ylim <- round(grid.objects$ylim)
dSS <- grid.objects$dSS
res <- grid.objects$res
x.vals <- grid.objects$x.vals
y.vals <- grid.objects$y.vals
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)
# rsf.cov <- 1*(grid.objects$rsf.cov>0)

#simulate a D.cov
set.seed(132402)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(10,150))[[2]]
D.cov <- as.numeric(scale(D.cov))
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov")

#make state space mask - just making a circle
dists <- sqrt((dSS[,1]-mean(dSS[,1]))^2+(dSS[,2]-mean(dSS[,2]))^2)
InSS <- 1*(dists<175)
image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="D.cov")

#simulate an RSF cov. Randomly distribute some cells animals really like
rsf.cov <- rbinom(n.cells,1,p=0.1)
image(x.vals,y.vals,matrix(rsf.cov,n.cells.x,n.cells.y),main="rsf.cov")

#simulate effort - I admit, I am being lazy here, but this gives us
#effort that varies over space with replication. You need some minimal
#level of spatial replication in survey effort
K <- 8
#Let's implement a state space buffer and only sample inside that.
#buffer should be at least 3sigma (which we haven't defined, yet, I use 7.8 below)
search.buff <- 1*(dists<150)

#but let's consider we don't survey the entire state space on each occasion
#let's say there are 4 sections and we do 1 per occasion
idx.list <- vector("list",K)
idx.list[[1]] <- which(dSS[,1]<175&dSS[,2]<200&search.buff)
idx.list[[2]] <- which(dSS[,1]<175&dSS[,2]>=200&search.buff)
idx.list[[3]] <- which(dSS[,1]>=175&dSS[,2]>=200&search.buff)
idx.list[[4]] <- which(dSS[,1]>=175&dSS[,2]<200&search.buff)
#must have replication, using same survey cells as above
idx.list[[5]] <- which(dSS[,1]<175&dSS[,2]<200&search.buff)
idx.list[[6]] <- which(dSS[,1]<175&dSS[,2]>=200&search.buff)
idx.list[[7]] <- which(dSS[,1]>=175&dSS[,2]>=200&search.buff)
idx.list[[8]] <- which(dSS[,1]>=175&dSS[,2]<200&search.buff)


effort <- survey <- matrix(0,n.cells,K)
for(k in 1:K){
  tmp <- grf(n.cells,grid=dSS,cov.pars=c(10,50))[[2]]
  effort[idx.list[[k]],k] <- tmp[idx.list[[k]]]
  survey[idx.list[[k]],k] <- 1
}

#visualize effort summed over all occasions
image(x.vals,y.vals,matrix(rowSums(effort),n.cells.x,n.cells.y),main="Search Effort")

# visualize effort on each occasion, scroll through k,
# can see different areas surveyed on each occasions 1-4 and 5-8
# k <- 1
# k <- k+1
# image(x.vals,y.vals,matrix(effort[,k],n.cells.x,n.cells.y),main="Search Effort")

#scale effort, but only in surveyed cell/occasions
effort[survey==1] <- as.numeric(scale(effort[survey==1]))


D.beta0 <- -7.0 #baseline D
D.beta1 <- 1.25 #density coefficient 
rsf.beta <- 4 #rsf coefficient
beta.p.int <- -0.5 #baseline detection prob
beta.p.effort <- 2.5 #effort effect no detection prob
sigma <- 7.8 #spatial scale of availability distribution
n.tel.inds <- 30 #number of telemetry individuals
K.tel <- 20 #number of telemetry locations per individual

#visualize distribution of p across cells given search effort and detection parameters
p.test <- plogis(beta.p.int + beta.p.effort*effort)
par(mfrow=c(1,1),ask=FALSE)
hist(p.test[survey==1])

set.seed(3234) #change this for new data set

data <- sim.RY(D.beta0=D.beta0,D.beta1=D.beta1,rsf.beta=rsf.beta,
               sigma=sigma,beta.p.int=beta.p.int,beta.p.effort=beta.p.effort,
               xlim=xlim,ylim=ylim,res=res,InSS=InSS,D.cov=D.cov,rsf.cov=rsf.cov,
               effort=effort,survey=survey,n.tel.inds=n.tel.inds,K=K,K.tel=K.tel)

data$truth$lambda #expected abundance from D cov inputs
data$truth$N #simulated realized abundance
data$truth$n.cap #number of inds captured
table(rowSums(data$capture$y)) #number of inds captures X times

#given search effort and use distributions, here are individuals cumulative detection probabilities
#spatial het in space use and effort lead to lots of het in p
#this is relevant to doing nonspatial capture-recapture
hist(data$summaries$p.marg.i,breaks=50,main="Individual Cumulative Detection Probabilities|Search Effort",
     xlim=c(0,1),xlab="Detection Prob")


##Fit Model
library(nimble)
library(coda)
nimbleOptions(determinePredictiveNodesInModel = FALSE)
source("init.data.RY.knownID.R")
source("NimbleModel knownID.R")
source("Nimble Functions knownID.R")
source("sSampler Dcov RSF Marginal.R")
M <- 300 #data augmentation limit. Must be larger than simulate N. If N posterior hits M, need to raise M and try again.
if(M<=data$truth$N)stop("Raise M to be larger than simulated N.")

inits <- list(sigma=7.8)
nimbuild <- init.data.RY.knownID(data=data,inits=inits,M=M)

n.surveyed.cells <- colSums(survey)
max.surveyed.cells <- max(n.surveyed.cells)
surveyed.cells <- matrix(0,max.surveyed.cells,K)
surveyed.cells.effort <- matrix(0,max.surveyed.cells,K)
for(k in 1:K){
  surveyed.cells[1:n.surveyed.cells[k],k] <- which(survey[,k]==1)
  surveyed.cells.effort[1:n.surveyed.cells[k],k] <- effort[surveyed.cells[1:n.surveyed.cells[k],k],k]
}

Niminits <- list(z=nimbuild$z,s=nimbuild$s,
                 N=nimbuild$N,D.beta0=-4,D.beta1=0,
                 sigma=inits$sigma,
                 beta.p.int=c(0),beta.p.effort=c(2),
                 rsf.beta=2,
                 s.tel=apply(data$telemetry$u.tel,c(1,3),mean,na.rm=TRUE))

#constants for Nimble
constants <- list(M=M,K=K,
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
Nimdata <- list(y=nimbuild$y,u=nimbuild$u,u.cell=nimbuild$u.cell,
                u.tel=data$telemetry$u.tel,u.cell.tel=data$telemetry$u.cell.tel,
                cells=data$constants$cells,InSS=data$constants$InSS,
                z=nimbuild$z,
                dummy.data=rep(1,M),
                rsf.cov=rsf.cov)

# set parameters to monitor
parameters<-c('beta.p.int','beta.p.effort','rsf.beta','D.beta1',
              'sigma','N','D.beta0','lambda')

#can also monitor a different set of parameters with a different thinning rate
nt <- 2 #thinning rate

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#tell nimble which nodes to configure so we don't waste time for samplers we will replace below
config.nodes <- c("beta.p.int","sigma",'D.beta0','D.beta1','rsf.beta','beta.p.effort')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = FALSE,
                      nodes=config.nodes) 

###*required* sampler replacement for "alternative data augmentation" N/z update
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
# conf$removeSampler("N")
#nodes used for update, calcNodes + z nodes
y.nodes <- Rmodel$expandNodeNames(paste("y[1:",M,",1:",K,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,y.nodes)
ind.detected <- 1*(rowSums(nimbuild$y)>0)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,K=K,M=M,ind.detected=ind.detected,
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
                  type = 'sSamplerDcovRSF.tel',control = list(i=i,xlim=xlim,ylim=ylim,
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
start.time2<-Sys.time()
Cmcmc$run(1000,reset=FALSE) #can keep running this line to extend sampler
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[20:nrow(mvSamples),])) #discarding some burnin here. Can't plot 1's sample which is all NA

data$truth$lambda #target expected abundance
data$truth$N #target realized abundance

par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(Cmodel$lambda.cell,n.cells.x,n.cells.y),main="Final Iteration Spatial D and ACs")
points(Cmodel$s[Cmodel$z==1,],pch=16,cex=0.75,col="grey30")
