getCellR = function(u,res,cells,xlim,ylim){
  inout <- 1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
  if(inout==1){
    this.cell.init <- cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
  }else{
    this.cell.init <- 0
  }
  return(this.cell.init)
}

init.data.RY.knownID <- function(data=NA,M=NA,inits=NA){
  data <- c(data$constants,data$capture,data$telemetry)
  u.obs <- data$u.obs
  K <- data$K
  xlim <- data$xlim
  ylim <- data$ylim
  cells <- data$cells
  dSS <- data$dSS
  InSS <- data$InSS

  #get some inits, actually sigma is all we need
  sigma <- inits$sigma
  
  n <- nrow(data$y)
  y <- matrix(0,nrow=M,ncol=K)
  y[1:n,] <- data$y
  u <- array(NA,dim=c(M,K,2))
  u[1:n,,] <- data$u.obs
  u.cell <- matrix(0,nrow=M,ncol=K) #must set unobserved u.cells to 0
  u.cell[1:n,] <- data$u.cell
  
  #Initialize z, just using observed z's
  z.init <- c(rep(1,n),rep(0,M-n))
  s.init <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
  for(i in 1:M){
    if(sum(y[i,])>0){#if captured
      #set s init
      s.init[i,] <- colMeans(u[i,,],na.rm=TRUE)
    }
  }
  #move any initialized outside state space
  for(i in 1:M){
    s.cell.init <- getCellR(s.init[i,],data$res,cells,xlim,ylim)
    if(InSS[s.cell.init]==0){#not in SS, move to nearest cell
      dists <- sqrt((dSS[s.cell.init,1]-dSS[,1])^2+(dSS[s.cell.init,2]-dSS[,2])^2)
      dists[InSS==0] <- Inf
      pick <- which(dists==min(dists))[1] #if more than 1, just use first
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

  #flatten observed continuous locations so the fitted model only creates duInCell nodes for detections
  obs.idx <- which(y==1,arr.ind=TRUE)
  obs.ID <- obs.idx[,1]
  obs.k <- obs.idx[,2]
  n.samples <- nrow(obs.idx)
  obs.cell <- rep(0,n.samples)
  u.obs2D <- matrix(NA,nrow=n.samples,ncol=2)
  for(l in 1:n.samples){
    obs.cell[l] <- u.cell[obs.ID[l],obs.k[l]]
    u.obs2D[l,] <- u[obs.ID[l],obs.k[l],]
  }
  return(list(y=y,u=u,u.cell=u.cell,z=z.init,s=s.init,s.tel=s.tel.init,N=sum(z.init),
              obs.ID=obs.ID,obs.k=obs.k,obs.cell=obs.cell,u.obs2D=u.obs2D,n.samples=n.samples))
}