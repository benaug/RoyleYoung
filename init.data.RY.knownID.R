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
  data <- c(data$constants,data$capture) #restructure data list
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
    s.cell.init <- getCellR(s.init[i,],res,cells,xlim,ylim)
    if(InSS[s.cell.init]==0){#not in SS, move to nearest cell
      dists <- sqrt((dSS[s.cell.init,1]-dSS[,1])^2+(dSS[s.cell.init,2]-dSS[,2])^2)
      dists[InSS==0] <- Inf
      pick <- which(dists==min(dists))[1] #if more than 1, just use first
      s.init[i,] <- dSS[pick,]
    }
  }
  return(list(y=y,u=u,u.cell=u.cell,z=z.init,s=s.init,N=sum(z.init)))
}