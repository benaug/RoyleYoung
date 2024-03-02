
sim.RY <-
  function(D.beta0=NA,D.beta1=NA,rsf.beta=NA,sigma=NA,beta.p.int=NA,beta.p.effort=NA,
           xlim=NA,ylim=NA,res=NA,InSS=NA,D.cov=NA,rsf.cov=NA,effort=NA,survey=NA,
           K=NA,K.tel=NA,n.tel.inds=NA){
    if(K==1)stop("K must be >1")
    if(xlim[1]!=0|ylim[1]!=0)stop("xlim and ylim must start at 0.")
    if((diff(range(xlim))/res)%%1!=0)stop("The range of xlim must be divisible by 'res'")
    if((diff(range(ylim))/res)%%1!=0)stop("The range of ylim must be divisible by 'res'")
    # make discrete state space
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res)
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res)
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    cellArea <- res^2
    
    #Density Model
    D.intercept <- exp(D.beta0)*cellArea
    lambda.cell <- InSS*exp(D.beta1*D.cov)
    pi.denom <- sum(lambda.cell)
    pi.cell <- lambda.cell/pi.denom
    lambda <- D.intercept*pi.denom
    lambda.cell.plot <- lambda.cell
    lambda.cell.plot[lambda.cell==0] <- -Inf
    image(x.vals,y.vals,matrix(lambda.cell.plot,length(x.vals),length(y.vals)),
          main="Spatially Explicit D and Realized ACs",xlab="X",ylab="Y")
    N <- rpois(1,lambda)

    #Activity centers
    # s.cell <- rcat(N,pi.cell)
    s.cell <- sample(1:n.cells,N,replace=TRUE,prob=pi.cell)
    s <- matrix(NA,N,2)
    for(i in 1:N){
      s.xlim <- dSS[s.cell[i],1] + c(-res,res)/2
      s.ylim <- dSS[s.cell[i],2] + c(-res,res)/2
      s[i,1] <- runif(1,s.xlim[1],s.xlim[2])
      s[i,2] <- runif(1,s.ylim[1],s.ylim[2])
    }
    points(s[,1],s[,2],pch=16)
    
    #RSF Model
    rsf <- exp(rsf.beta*rsf.cov)
    #get BVN availability and use distributions
    avail.dist <- use.dist <- matrix(NA,N,n.cells)
    for(i in 1:N){
      avail.dist[i,] <- getAvail(s=s[i,1:2],sigma=sigma,res=res,x.vals=x.vals,
                                 y.vals=y.vals,n.cells.x=n.cells.x,n.cells.y=n.cells.y)
      use.dist[i,] <- rsf*avail.dist[i,]
      use.dist[i,] <- use.dist[i,]/sum(use.dist[i,])
    }
    
    #simulate u's
    u <- array(NA,dim=c(N,K,2))
    u.cell <- matrix(NA,N,K)
    for(i in 1:N){
      for(k in 1:K){
        # u.cell[i,k] <- rcat(1,use.dist[i,])
        u.cell[i,k] <- sample(1:n.cells,1,replace=TRUE,prob=use.dist[i,])
        u.xlim <- dSS[u.cell[i,k],1] + c(-res,res)/2
        u.ylim <- dSS[u.cell[i,k],2] + c(-res,res)/2
        u[i,k,1] <- rtruncnorm(1,a=u.xlim[1],b=u.xlim[2],
                                   mean=s[i,1],sd=sigma)
        u[i,k,2] <- rtruncnorm(1,a=u.ylim[1],b=u.ylim[2],
                                   mean=s[i,2],sd=sigma)
      }
    }
    rsf.cov.plot <- rsf.cov
    # rsf.cov.plot[InSS==0] <- -Inf
    image(x.vals,y.vals,matrix(rsf.cov.plot,length(x.vals),length(y.vals)),
          main="Individual ACs and Site Use",xlab="X",ylab="Y")
    grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
    for(i in 1:N){
      for(k in 1:K){
        lines(x=c(s[i,1],u[i,k,1]),y=c(s[i,2],u[i,k,2]),col="black")
      }
    }
    # points(u[,,1],u[,,2],pch=16,cex=1,col="royalblue")
    points(u[,,1],u[,,2],pch=16,cex=1,col=rgb(65,105,225,100,maxColorValue=225))
    points(s[,1],s[,2],pch=16,cex=1.25,col="grey30")
    # points(s[,1],s[,2],pch=16,cex=1.5,col=rgb(77,77,77,100,maxColorValue=225))
    
    #simulate detection
    p <- plogis(beta.p.int+beta.p.effort*effort)
    p[survey==0] <- 0 #zero out unsurveyed areas
    y <- matrix(NA,N,K)
    p.ik <- matrix(NA,N,K) #p conditional on selecting cell
    p.i <- rep(NA,N)
    p.marg.ik <- matrix(NA,N,K) #p unconditional. conditional on search effort
    p.marg.i <- rep(NA,N)
    for(i in 1:N){
      for(k in 1:K){
        y[i,k] <- rbinom(1,1,p[u.cell[i,k],k])
        p.ik[i,k] <- p[u.cell[i,k],k]
        p.marg.ik[i,k] <- sum(use.dist[i,]*p[,k])
      }
      p.i[i] <- 1-prod(1-p.ik[i,])
      p.marg.i[i] <- 1-prod(1-p.marg.ik[i,])
    }
    cum.p.c <- 1-apply(1-p,1,prod)
    cum.p.c.plot <- cum.p.c
    cum.p.c.plot[cum.p.c==0] <- -Inf
    image(x.vals,y.vals,matrix(cum.p.c.plot,length(x.vals),length(y.vals)),
          main="Cum Detect Prob over K|u.cell",xlab="X",ylab="Y")
    grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
    detected <- which(rowSums(y)>0)
    points(s[detected,1],s[detected,2],pch=16,col="grey30",cex=1.25)
    points(u[,,1][y==1],u[,,2][y==1],pch=16,col="royalblue")
    # points(u[,,1][y==1],u[,,2][y==1],pch=16,cex=1,col=rgb(65,105,225,100,maxColorValue=225))
    for(i in detected){
      for(k in 1:K){
        if(y[i,k]==1){
          lines(x=c(s[i,1],u[i,k,1]),y=c(s[i,2],u[i,k,2]))
        }
      }
    }
    
    ####Telemetry data#####
    #simulating from same D model here
    s.tel.cell <- rcat(n.tel.inds,pi.cell)
    s.tel <- matrix(NA,n.tel.inds,2)
    for(i in 1:n.tel.inds){
      s.xlim <- dSS[s.tel.cell[i],1] + c(-res,res)/2
      s.ylim <- dSS[s.tel.cell[i],2] + c(-res,res)/2
      s.tel[i,1] <- runif(1,s.xlim[1],s.xlim[2])
      s.tel[i,2] <- runif(1,s.ylim[1],s.ylim[2])
    }
    #get BVN availability and use distributions
    avail.dist.tel <- use.dist.tel <- matrix(NA,N,n.cells)
    for(i in 1:n.tel.inds){
      avail.dist.tel[i,] <- getAvail(s=s.tel[i,1:2],sigma=sigma,res=res,x.vals=x.vals,
                                 y.vals=y.vals,n.cells.x=n.cells.x,n.cells.y=n.cells.y)
      use.dist.tel[i,] <- rsf*avail.dist.tel[i,]
      use.dist.tel[i,] <- use.dist.tel[i,]/sum(use.dist.tel[i,])
    }
    #simulate u's
    u.tel <- array(NA,dim=c(n.tel.inds,K.tel,2))
    u.cell.tel <- matrix(NA,n.tel.inds,K.tel)
    u.xlim.tel <- u.ylim.tel <- array(NA,dim=c(n.tel.inds,K.tel,2))
    for(i in 1:n.tel.inds){
      for(k in 1:K.tel){
        # u.cell.tel[i,k] <- rcat(1,use.dist.tel[i,])
        u.cell.tel[i,k] <- sample(1:n.cells,1,replace=TRUE,prob=use.dist.tel[i,])
        u.xlim.tel[i,k,] <- dSS[u.cell.tel[i,k],1] + c(-res,res)/2
        u.ylim.tel[i,k,] <- dSS[u.cell.tel[i,k],2] + c(-res,res)/2
        u.tel[i,k,1] <- rtruncnorm(1,a=u.xlim.tel[i,k,1],b=u.xlim.tel[i,k,2],
                                   mean=s.tel[i,1],sd=sigma)
        u.tel[i,k,2] <- rtruncnorm(1,a=u.ylim.tel[i,k,1],b=u.ylim.tel[i,k,2],
                                   mean=s.tel[i,2],sd=sigma)
      }
    }
    image(x.vals,y.vals,matrix(rsf.cov.plot,length(x.vals),length(y.vals)),
          main="Telemetry Individual ACs and Site Use",xlab="X",ylab="Y")
    grid(n.cells.x,n.cells.y,lwd=1,lty=1,col="grey80")
    for(i in 1:n.tel.inds){
      for(k in 1:K.tel){
        lines(x=c(s.tel[i,1],u.tel[i,k,1]),y=c(s.tel[i,2],u.tel[i,k,2]),col="black")
      }
    }
    # points(u.tel[,,1],u.tel[,,2],pch=16,cex=1,col="royalblue")
    points(u.tel[,,1],u.tel[,,2],pch=16,cex=1,col=rgb(65,105,225,100,maxColorValue=225))
    points(s.tel[,1],s.tel[,2],pch=16,cex=1.25,col="grey30")
    
    n.locs.ind <- rep(K.tel,n.tel.inds)
    n.locs.max <- max(n.locs.ind)

    #discard uncaptured inds and disaggregate data
    caught <- which(rowSums(y)>0)
    n <- length(caught)
    y <- y[caught,]
    # u.full <- u
    u <- u.obs <- u[caught,,]
    u.cell <- u.cell[caught,]
    #discard unobserved u's
    u.obs[,,1][y==0] <- NA
    u.obs[,,2][y==0] <- NA
    u.cell.obs <- u.cell
    u.cell.obs[y==0] <- 0 #set unobserved to 0
    # s.full <- s
    s <- s[caught,]
    s.cell <- s.cell[caught]
    
    tmp <- which(y==1,arr.ind=TRUE)
    ID <- tmp[,1]
    this.k <- tmp[,2]
    n.samples <- length(this.k)
    
    u.obs2D <- matrix(NA,nrow=n.samples,ncol=2)
    for(l in 1:n.samples){
      u.obs2D[l,] <- u.obs[ID[l],this.k[l],]
    }

    #check data disaggregation
    y.check <- matrix(0,nrow=n,ncol=K)
    u.check <- u.obs*NA
    for(i in 1:length(ID)){
      y.check[ID[i],this.k[i]] <- 1
      u.check[ID[i],this.k[i],] <- u.obs2D[i,]
    }
    if(!all(y==y.check))stop("Error rebuilding data. Report Bug.")
    if(!all(u.check==u.obs,na.rm=TRUE))stop("Error rebuilding data. Report Bug.")
    
    constants <- list(K=K,K.tel=K.tel,xlim=xlim,ylim=ylim,dSS=dSS,res=res,cells=cells,x.vals=x.vals,y.vals=y.vals,
                      n.tel.inds=n.tel.inds,n.locs.ind=n.locs.ind,n.samples=n.samples,
                      InSS=InSS,n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    truth <- list(lambda=lambda,lambda.cell=lambda.cell,rsf=rsf,N=N,s=s,u=u,s.cell=s.cell,u.cell=u.cell,
                  ID=ID,n=n,s.tel=s.tel,s.tel.cell=s.tel.cell,
                  use.dist=use.dist,avail.dist=avail.dist)
    capture <- list(y=y,u.obs=u.obs, u.cell=u.cell.obs, #known ID
                    this.k=this.k,u.obs2D=u.obs2D) #unknown ID
    telemetry <- list(u.tel=u.tel,u.cell.tel=u.cell.tel,u.xlim.tel=u.xlim.tel,
                      u.ylim.tel=u.ylim.tel) #observed telemetry data
    summaries <- list(p=p,p.ik=p.ik,p.i=p.i,p.marg.ik=p.marg.ik,p.marg.i=p.marg.i,cum.p.c=cum.p.c)
    out <-list(constants=constants,truth=truth,capture=capture,telemetry=telemetry,summaries=summaries)
    return(out)
  }
