getCellR = function(u,res,cells,xlim,ylim){
  inout <- 1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
  if(inout==1){
    this.cell.init <- cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
  }else{
    this.cell.init <- 0
  }
  return(this.cell.init)
}

init.data.RY.unknownID <- function(data=NA,M=NA,inits=NA){
  data <- c(data$constants,data$capture,data$telemetry)
  K <- data$K
  xlim <- data$xlim
  ylim <- data$ylim
  cells <- data$cells
  dSS <- data$dSS
  InSS <- data$InSS
  n.samples <- length(data$this.k)
  this.k <- data$this.k
  u.obs <- data$u.obs2D

  if(n.samples>M)stop("M must be at least as large as the number of observed samples for this initializer.")

  # init 1 observed sample per individual.
  ID <- seq_len(n.samples)
  y.true <- matrix(0,nrow=M,ncol=K)
  u <- array(0,dim=c(M,K,2))
  u.cell <- matrix(0,nrow=M,ncol=K)
  u.cell.survey <- matrix(0,nrow=M,ncol=K)
  obs.cell <- rep(0,n.samples)

  for(l in 1:n.samples){
    i <- ID[l]
    k <- this.k[l]
    y.true[i,k] <- 1
    u[i,k,1:2] <- u.obs[l,1:2]
    obs.cell[l] <- getCellR(u.obs[l,1:2],data$res,cells,xlim,ylim)
    if(obs.cell[l]==0)stop("Observed location outside state space.")
    u.cell[i,k] <- obs.cell[l]
    surveyed.cells.k <- which(data$survey[,k]==1)
    idx <- which(surveyed.cells.k==obs.cell[l])
    if(length(idx)!=1)stop("Observed location is not in a surveyed cell on its detection occasion.")
    u.cell.survey[i,k] <- idx
  }

  #Initialize z from the current latent-ID allocation.
  z.init <- c(rep(1,n.samples),rep(0,M-n.samples))
  s.init <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
  for(i in 1:M){
    if(sum(y.true[i,])>0){
      k.idx <- which(y.true[i,]>0)
      s.init[i,1] <- mean(u[i,k.idx,1])
      s.init[i,2] <- mean(u[i,k.idx,2])
    }
  }

  #move any initialized outside state space
  for(i in 1:M){
    s.cell.init <- getCellR(s.init[i,],data$res,cells,xlim,ylim)
    if(InSS[s.cell.init]==0){
      dists <- sqrt((dSS[s.cell.init,1]-dSS[,1])^2+(dSS[s.cell.init,2]-dSS[,2])^2)
      dists[InSS==0] <- Inf
      pick <- which(dists==min(dists))[1]
      s.init[i,] <- dSS[pick,]
    }
  }

  #z=0 individuals have activity centers fixed at 0
  s.init[z.init==0,] <- 0

  #initialize telemetry ACs from observed locations and move any outside InSS to the nearest valid cell.
  s.tel.init <- apply(data$u.tel,c(1,3),mean,na.rm=TRUE)
  for(i in 1:data$n.tel.inds){
    s.cell.init <- getCellR(s.tel.init[i,],data$res,cells,xlim,ylim)
    if(InSS[s.cell.init]==0){
      dists <- sqrt((dSS[s.cell.init,1]-dSS[,1])^2+(dSS[s.cell.init,2]-dSS[,2])^2)
      dists[InSS==0] <- Inf
      pick <- which(dists==min(dists))[1]
      s.tel.init[i,] <- dSS[pick,]
    }
  }

  capcounts <- rowSums(y.true)
  N.init <- sum(z.init)

  return(list(y.true=y.true,z=z.init,s=s.init,s.tel=s.tel.init,u=u,ID=ID,n.samples=n.samples,
              this.k=this.k,u.cell=u.cell,u.cell.survey=u.cell.survey,obs.cell=obs.cell,
              capcounts=capcounts,N=N.init))
}
