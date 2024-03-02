getCellR = function(u,res,cells,xlim,ylim){
  inout=1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
  if(inout==1){
    this.cell=cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
  }else{
    this.cell=0
  }
  return(this.cell)
}

init.data.RY.unknownID=function(data=NA,M=NA,inits=NA){
  data <- c(data$constants,data$capture) #restructure data list
  library(abind)
  this.k <- data$this.k
  u.obs <- data$u.obs2D
  K <- data$K
  xlim <- data$xlim
  ylim <- data$ylim
  n.samples <- length(this.k)
  #get some inits
  sigma <- inits$sigma
  cells <- data$cells
  
  #init 1 obs per bear
  ID <- seq(1:n.samples)
  
  #Build y.true, u
  y.true <- matrix(0,nrow=M,ncol=K)
  u.init <- array(NA,dim=c(M,K,2))
  for(l in 1:n.samples){
    y.true[ID[l],this.k[l]] <- 1
    u.init[ID[l],this.k[l],] <- u.obs[l,]
  }
  
  known.vector <- c(rep(1,max(ID)),rep(0,M-max(ID)))
  
  #Initialize z
  z <- known.vector
  
  library(truncnorm)
  
  s.init <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
  for(i in 1:M){
    if(sum(y.true[i,])>0){#if captured
      #set s init
      s.init[i,] <- colMeans(u.init[i,,],na.rm=TRUE)
    }
  }
  #move any initialized outside state space
  for(i in 1:M){
    s.cell <- getCellR(s.init[i,],data$res,data$cells,data$xlim,data$ylim)
    if(data$InSS[s.cell]==0){#not in SS, move to nearest cell
      dists <- sqrt((data$dSS[s.cell,1]-data$dSS[,1])^2+(data$dSS[s.cell,2]-data$dSS[,2])^2)
      dists[data$InSS==0] <- Inf
      pick <- which(dists==min(dists))[1] #if more than 1, just use first
      s.init[i,] <- data$dSS[pick,]
    }
  }
  
  #Need u cell inits for rsf
  u.cell.init <- matrix(0,M,K)
  for(i in 1:M){
    for(k in 1:K){
      if(!is.na(u.init[i,k,1])){
        u.cell.init[i,k] <- getCellR(u.init[i,k,],data$res,data$cells,data$xlim,data$ylim)
      }
    }
  }
  
  any(u.cell.init==0,na.rm=TRUE) #should be false
  
  z.init <- z
  z.init[is.na(z.init)] <- 0
  N.init <- sum(z.init) #must use these z and N inits
  
  
  return(list(y.true=y.true,z=z.init,s=s.init,u=u.init,ID=ID,n.samples=n.samples,xlim=xlim,ylim=ylim,
              this.k=this.k,u.cell=u.cell.init,N=N.init))
}