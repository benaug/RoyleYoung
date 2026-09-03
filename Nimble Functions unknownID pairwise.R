#Within-cell continuous likelihood using stored 1-D availability probabilities.
duInCell <- nimbleFunction(
  run = function(x = double(1),s = double(1),u.cell = double(0),sigma = double(0),
                 n.cells.x = integer(0),res = double(0),avail.x = double(1),avail.y = double(1),
                 log = integer(0)) {
    returnType(double(0))
    if(u.cell>0){
      u.cell.x <- u.cell%%n.cells.x
      u.cell.y <- floor(u.cell/n.cells.x)+1
      if(u.cell.x==0){
        u.cell.x <- n.cells.x
        u.cell.y <- u.cell.y-1
      }
      if(avail.x[u.cell.x]<=0 | avail.y[u.cell.y]<=0){
        return(-Inf)
      }
      logProb <- dnorm(x[1],s[1],sigma,log=TRUE)-log(avail.x[u.cell.x]) +
        dnorm(x[2],s[2],sigma,log=TRUE)-log(avail.y[u.cell.y])
    }else{
      logProb <- 0
    }
    return(logProb)
  }
)

ruInCell <- nimbleFunction(
  run = function(n = integer(0),s = double(1),u.cell = double(0),sigma = double(0),
                 n.cells.x = integer(0),res = double(0),avail.x = double(1),avail.y = double(1)) {
    returnType(double(1))
    return(c(0,0))
  })

dRYmargFactored <- nimbleFunction(
  run = function(x = double(0),u.cell.survey = double(0),z = integer(0),
                 p.rsf = double(1),surveyed.cell.x = double(1),surveyed.cell.y = double(1),
                 n.surveyed.cells = integer(0),avail.x = double(1),avail.y = double(1),
                 use.denom = double(0),log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(x==1){
        return(-Inf)
      }
      return(0)
    }
    if(use.denom<=0){
      return(-Inf)
    }
    if(x==1){
      if(u.cell.survey<=0){
        return(-Inf)
      }
      use.detect <- p.rsf[u.cell.survey]*
        avail.x[surveyed.cell.x[u.cell.survey]]*
        avail.y[surveyed.cell.y[u.cell.survey]]
      if(use.detect<=0){
        return(-Inf)
      }
      logProb <- log(use.detect)-log(use.denom)
    }else{
      detect.mass <- 0
      for(c in 1:n.surveyed.cells){
        detect.mass <- detect.mass + p.rsf[c]*
          avail.x[surveyed.cell.x[c]]*avail.y[surveyed.cell.y[c]]
      }
      detect.prob <- detect.mass/use.denom
      if(detect.prob>1 | detect.prob<0){
        return(-Inf)
      }
      logProb <- log(1-detect.prob)
    }
    return(logProb)
  }
)

rRYmargFactored <- nimbleFunction(
  run = function(n = integer(0),u.cell.survey = double(0),z = integer(0),
                 p.rsf = double(1),surveyed.cell.x = double(1),surveyed.cell.y = double(1),
                 n.surveyed.cells = integer(0),avail.x = double(1),avail.y = double(1),
                 use.denom = double(0)) {
    returnType(double(0))
    return(1)
  })

#used in data simulator
getAvail <- nimbleFunction(
  run = function(s = double(1),sigma=double(0),res=double(0),x.vals=double(1),y.vals=double(1),n.cells.x=integer(0),n.cells.y=integer(0)) {
    returnType(double(1))
    avail.dist.x <- rep(0,n.cells.x)
    avail.dist.y <- rep(0,n.cells.y)
    delta <- 1e-12 #this sets the degree of trimming used to get individual availability distributions
    x.limits <- qnorm(c(delta,1-delta),mean=s[1],sd=sigma)
    y.limits <- qnorm(c(delta,1-delta),mean=s[2],sd=sigma)
    #convert to grid edges instead of centroids
    x.vals.edges <- c(x.vals - res/2, x.vals[n.cells.x]+0.5*res)
    y.vals.edges <- c(y.vals - res/2, y.vals[n.cells.y]+0.5*res)
    #trim in x and y direction
    if(x.vals.edges[1]<x.limits[1]){
      x.start <- floor((x.limits[1] - x.vals.edges[1]) / res) + 1
    }else{
      x.start <- 1
    }
    if(x.vals.edges[n.cells.x]>x.limits[2]){
      x.stop <- ceiling((x.limits[2] - x.vals.edges[1]) / res)
    }else{
      x.stop <- n.cells.x
    }
    if(y.vals.edges[1]<y.limits[1]){
      y.start <- floor((y.limits[1] - y.vals.edges[1]) / res) + 1
    }else{
      y.start <- 1
    }
    if(y.vals.edges[n.cells.y]>y.limits[2]){
      y.stop <- ceiling((y.limits[2] - y.vals.edges[1]) / res)
    }else{
      y.stop <- n.cells.y
    }
    pnorm.x <- rep(0,n.cells.x+1)
    pnorm.y <- rep(0,n.cells.y+1)
    #get pnorms
    for(l in x.start:(x.stop+1)){
      pnorm.x[l] <- pnorm(x.vals.edges[l],mean=s[1],sd=sigma)
    }
    for(l in y.start:(y.stop+1)){
      pnorm.y[l] <- pnorm(y.vals.edges[l],mean=s[2],sd=sigma)
    }
    for(l in (x.start):(x.stop)){
      avail.dist.x[l] <- pnorm.x[l+1] - pnorm.x[l]
    }
    for(l in (y.start):(y.stop)){
      avail.dist.y[l] <- pnorm.y[l+1] - pnorm.y[l]
    }
    avail.dist.tmp <- matrix(0,n.cells.x,n.cells.y)
    sum.dist <- 0
    for(i in x.start:x.stop){
      for(j in y.start:y.stop){
        avail.dist.tmp[i,j] <- avail.dist.x[i]*avail.dist.y[j]
        sum.dist <- sum.dist + avail.dist.tmp[i,j]
      }
    }
    avail.dist <- c(avail.dist.tmp)
    #if any probability mass is outside state space, normalize
    if(sum.dist<1){
      avail.dist <- avail.dist/sum.dist
    }
    return(avail.dist)
  }
)

#Factored replacement for getAvail in the fitted model. Cell availability is
#avail.x[cell.x]*avail.y[cell.y], so only the two 1-D vectors are stored.
getAvail1D <- nimbleFunction(
  run = function(s = double(0),sigma = double(0),res = double(0),vals.edges = double(1),
                 n.cells = integer(0),avail.z = double(0),z = integer(0)) {
    returnType(double(1))
    avail <- rep(0,n.cells)
    if(z==0){
      return(avail)
    }
    lower <- s-sigma*avail.z
    upper <- s+sigma*avail.z
    if(vals.edges[1]<lower){
      idx.start <- floor((lower-vals.edges[1])/res)+1
    }else{
      idx.start <- 1
    }
    if(vals.edges[n.cells]>upper){
      idx.stop <- ceiling((upper-vals.edges[1])/res)
    }else{
      idx.stop <- n.cells
    }
    pnorm.vals <- rep(0,n.cells+1)
    for(l in idx.start:(idx.stop+1)){
      pnorm.vals[l] <- pnorm(vals.edges[l],mean=s,sd=sigma)
    }
    for(l in idx.start:idx.stop){
      avail[l] <- pnorm.vals[l+1]-pnorm.vals[l]
    }
    return(avail)
  }
)

#normalizing constant for RSF-weighted use. The scan identifies the
#nonzero x/y ranges so the nested loop does not use the full grid when sigma is small relative to grid extent.
getUseDenom <- nimbleFunction(
  run = function(rsf = double(1),avail.x = double(1),avail.y = double(1),
                 n.cells.x = integer(0),n.cells.y = integer(0),z = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }
    x.start <- n.cells.x
    x.stop <- 0
    y.start <- n.cells.y
    y.stop <- 0
    for(i in 1:n.cells.x){
      if(avail.x[i]>0){
        if(i<x.start){
          x.start <- i
        }
        if(i>x.stop){
          x.stop <- i
        }
      }
    }
    for(j in 1:n.cells.y){
      if(avail.y[j]>0){
        if(j<y.start){
          y.start <- j
        }
        if(j>y.stop){
          y.stop <- j
        }
      }
    }
    use.denom <- 0
    if(x.stop>0 & y.stop>0){
      for(j in y.start:y.stop){
        for(i in x.start:x.stop){
          cell <- i+(j-1)*n.cells.x
          use.denom <- use.denom + rsf[cell]*avail.x[i]*avail.y[j]
        }
      }
    }
    return(use.denom)
  }
)

#z-gated activity-center distribution combining pi.cell and the within-cell uniform density.
#When z=1, choose a cell from pi.cell and use a uniform density within that cell.
#When z=0, s is fixed at c(0,0), so inactive individuals do not retain latent ACs.
dAC <- nimbleFunction(
  run = function(x = double(1),pi.cell = double(1),res = double(0),
                 n.cells.x = integer(0),n.cells.y = integer(0),z = integer(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(x[1]==0 & x[2]==0){
        logProb <- 0
      }else{
        logProb <- -Inf
      }
    }else{
      cell.x <- trunc(x[1]/res)+1
      cell.y <- trunc(x[2]/res)+1
      cell <- cell.x+(cell.y-1)*n.cells.x
      if(pi.cell[cell]<=0){
        logProb <- -Inf
      }else{
        logProb <- log(pi.cell[cell])-2*log(res)
      }
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rAC <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(1),res = double(0),
                 n.cells.x = integer(0),n.cells.y = integer(0),z = integer(0)) {
    returnType(double(1))
    if(z==0){
      return(c(0,0))
    }else{
      n.cells <- n.cells.x*n.cells.y
      cell <- rcat(1,pi.cell[1:n.cells])
      cell.x <- cell%%n.cells.x
      cell.y <- floor(cell/n.cells.x)+1
      if(cell.x==0){
        cell.x <- n.cells.x
        cell.y <- cell.y-1
      }
      x.sim <- runif(1,(cell.x-1)*res,cell.x*res)
      y.sim <- runif(1,(cell.y-1)*res,cell.y*res)
      return(c(x.sim,y.sim))
    }
  }
)

getP <- nimbleFunction(
  run = function(surveyed.cells.effort=double(1),n.surveyed.cells=integer(0),
                 beta.p.int=double(0),beta.p.effort=double(0)) {
    returnType(double(1))
    p <- rep(0,n.surveyed.cells)
    for(l in 1:n.surveyed.cells){
      p[l] <- plogis(beta.p.int + beta.p.effort*surveyed.cells.effort[l])
    }#otherwise p=0, not surveyed
    return(p)
  }
)

#Telemetry likelihood with the categorical cell probability and within-cell
#truncated Normal combined. The availability cell masses cancel
dTelemetryFactored <- nimbleFunction(
  run = function(x = double(2),u.cell = double(1),s = double(1),sigma = double(0),
                 rsf = double(1),use.denom = double(0),
                 n.locs.ind = double(0),log = integer(0)) {
    returnType(double(0))
    if(use.denom<=0){
      return(-Inf)
    }
    logProb <- 0
    for(l in 1:n.locs.ind){
      this.cell <- u.cell[l]
      if(rsf[this.cell]<=0){
        return(-Inf)
      }
      logProb <- logProb + log(rsf[this.cell])-log(use.denom) +
        dnorm(x[l,1],s[1],sigma,log=TRUE) + dnorm(x[l,2],s[2],sigma,log=TRUE)
    }
    return(logProb)
  }
)

rTelemetryFactored <- nimbleFunction(
  run = function(n = integer(0),u.cell = double(1),s = double(1),sigma = double(0),
                 rsf = double(1),use.denom = double(0),
                 n.locs.ind = double(0)) {
    returnType(double(2))
    return(matrix(0,n.locs.ind,2))
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    K <- nimDim(y.true)[2]
    capcounts <- numeric(M,value=0)
    for(i in 1:M){
      capcounts[i] <- sum(y.true[i,1:K])
    }
    return(capcounts)
  }
)

Getncap <- nimbleFunction( 
  run = function(capcounts=double(1),ID=double(1)){
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    detected <- numeric(M,value=0)
    for(i in 1:M){
      if(capcounts[i]>0){
        detected[i] <- 1
      }
    }
    n.cap <- sum(detected)
    return(n.cap)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    K <- control$K
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    s.nodes <- control$s.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
    res <- control$res
    x.vals.edges <- control$x.vals.edges
    y.vals.edges <- control$y.vals.edges
    avail.z <- control$avail.z
    n.cells <- control$n.cells
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
  },
  run = function() {
    #Build undetected on/off lists once, then update after accepted proposals.
    #Detected status is dynamic under latent ID, so use the current capcounts.
    #Individuals with currently assigned detections are never included in z.on.
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in 1:M){
      if(model$capcounts[i]==0){
        if(model$z[i]==1){
          non.curr <- non.curr+1
          z.on[non.curr] <- i
        }else{
          noff.curr <- noff.curr+1
          z.off[noff.curr] <- i
        }
      }
    }

    for(up in 1:z.ups){
      updown <- rbinom(1,1,0.5)
      if(updown==0){ #subtract
        non.init <- non.curr
        if(non.init>0){
          pick.pos <- rcat(1,rep(1/non.init,non.init))
          pick <- z.on[pick.pos]
          N.init <- model$N[1]

          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- 0
          for(k in 1:K){
            y.idx <- pick+(k-1)*M
            lp.initial.y <- lp.initial.y + model$getLogProb(y.nodes[y.idx])
          }
          # lp.initial.s <- model$getLogProb(s.nodes[pick]) #cancels with reverse s proposal density

          #propose new N/z and unique inactive AC state
          model$N[1] <<- model$N[1]-1
          model$z[pick] <<- 0
          model$s[pick,1:2] <<- c(0,0)
          model$avail.x[pick,1:n.cells.x] <<- rep(0,n.cells.x)
          model$avail.y[pick,1:n.cells.y] <<- rep(0,n.cells.y)
          model$use.denom[pick] <<- 0

          #get proposed logprob for N
          lp.proposed.N <- model$calculate(N.node)
          #the y logProb is 0 when z=0, so skip calculating it until acceptance 
          lp.proposed.y <- 0
          # lp.proposed.s <- model$calculate(s.nodes[pick]) #is 0 for z=0; calculate after acceptance

          #s target and proposal terms cancel
          # lp.initial.s <- model$getLogProb(s.nodes[pick])
          # log.prop.back.s <- lp.initial.s
          # lp.proposed.s <- 0
          # log.prop.for.s <- 0

          #MH step
          #combinatorial/proposal correction is non.init/N.init
          #full ratio before s cancellation
          # log_MH_ratio <- (lp.proposed.N+lp.proposed.y+lp.proposed.s) -
          #   (lp.initial.N+lp.initial.y+lp.initial.s) + log(non.init/N.init) +
          #   log.prop.back.s-log.prop.for.s
          log_MH_ratio <- (lp.proposed.N+lp.proposed.y)-(lp.initial.N+lp.initial.y)+log(non.init/N.init)
          accept <- decide(log_MH_ratio)

          if(accept){
            #y and s were intentionally not calculated before the MH decision where their values are known/cancel;
            #calculate them now to synchronize accepted model logProbs
            for(k in 1:K){
              y.idx <- pick+(k-1)*M
              model$calculate(y.nodes[y.idx])
            }
            model$calculate(s.nodes[pick])
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["s",1][pick,1:2] <<- model[["s"]][pick,1:2]
            mvSaved["avail.x",1][pick,1:n.cells.x] <<- model[["avail.x"]][pick,1:n.cells.x]
            mvSaved["avail.y",1][pick,1:n.cells.y] <<- model[["avail.y"]][pick,1:n.cells.y]
            mvSaved["use.denom",1][pick] <<- model[["use.denom"]][pick]
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr-1
            noff.curr <- noff.curr+1
            z.off[noff.curr] <- pick
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["s"]][pick,1:2] <<- mvSaved["s",1][pick,1:2]
            model[["avail.x"]][pick,1:n.cells.x] <<- mvSaved["avail.x",1][pick,1:n.cells.x]
            model[["avail.y"]][pick,1:n.cells.y] <<- mvSaved["avail.y",1][pick,1:n.cells.y]
            model[["use.denom"]][pick] <<- mvSaved["use.denom",1][pick]
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1]<M){
          noff.init <- noff.curr
          if(noff.init>0){
            pick.pos <- rcat(1,rep(1/noff.init,noff.init))
            pick <- z.off[pick.pos]
            N.init <- model$N[1]

            #get initial logprobs for N and y
            lp.initial.N <- model$getLogProb(N.node)
            lp.initial.y <- 0
            # lp.initial.s <- model$getLogProb(s.nodes[pick]) #gated z=0 s logProb is 0

            #propose new N/z
            model$N[1] <<-  model$N[1] + 1
            model$z[pick] <<- 1
            
            #propose s from its dAC model prior when the individual is turned on
            model$s[pick,1:2] <<- rAC(1,pi.cell=model$pi.cell[1:n.cells],res=res,
                                      n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
            model$avail.x[pick,1:n.cells.x] <<- getAvail1D(s=model$s[pick,1],sigma=model$sigma[1],res=res,
                                                           vals.edges=x.vals.edges,n.cells=n.cells.x,
                                                           avail.z=avail.z,z=1)
            model$avail.y[pick,1:n.cells.y] <<- getAvail1D(s=model$s[pick,2],sigma=model$sigma[1],res=res,
                                                           vals.edges=y.vals.edges,n.cells=n.cells.y,
                                                           avail.z=avail.z,z=1)
            model$use.denom[pick] <<- getUseDenom(rsf=model$rsf[1:n.cells],
                                                  avail.x=model$avail.x[pick,1:n.cells.x],
                                                  avail.y=model$avail.y[pick,1:n.cells.y],
                                                  n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.node)
            lp.proposed.y <- 0
            for(k in 1:K){
              y.idx <- pick+(k-1)*M
              lp.proposed.y <- lp.proposed.y + model$calculate(y.nodes[y.idx])
            }
            # lp.proposed.s <- model$calculate(s.nodes[pick]) #cancels with proposal, calculate after acceptance
            
            #s target and proposal terms cancel 
            # lp.proposed.s <- model$calculate(s.nodes[pick])
            # log.prop.for.s <- lp.proposed.s
            # lp.initial.s <- 0
            # log.prop.back.s <- 0

            #MH step
            #combinatorial/proposal correction is (N.init+1)/(non.curr+1)
            #full ratio before s cancellation
            # log_MH_ratio <- (lp.proposed.N+lp.proposed.y+lp.proposed.s) -
            #   (lp.initial.N+lp.initial.y+lp.initial.s) + log((N.init+1)/(non.curr+1)) +
            #   log.prop.back.s-log.prop.for.s
            log_MH_ratio <- (lp.proposed.N+lp.proposed.y)-(lp.initial.N+lp.initial.y)+log((N.init+1)/(non.curr+1))
            accept <- decide(log_MH_ratio)

            if(accept){
              model$calculate(s.nodes[pick])
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["z",1][pick] <<- model[["z"]][pick]
              mvSaved["s",1][pick,1:2] <<- model[["s"]][pick,1:2]
              mvSaved["avail.x",1][pick,1:n.cells.x] <<- model[["avail.x"]][pick,1:n.cells.x]
              mvSaved["avail.y",1][pick,1:n.cells.y] <<- model[["avail.y"]][pick,1:n.cells.y]
              mvSaved["use.denom",1][pick] <<- model[["use.denom"]][pick]
              z.off[pick.pos] <- z.off[noff.curr]
              z.off[noff.curr] <- 0
              noff.curr <- noff.curr-1
              non.curr <- non.curr+1
              z.on[non.curr] <- pick
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["z"]][pick] <<- mvSaved["z",1][pick]
              model[["s"]][pick,1:2] <<- mvSaved["s",1][pick,1:2]
              model[["avail.x"]][pick,1:n.cells.x] <<- mvSaved["avail.x",1][pick,1:n.cells.x]
              model[["avail.y"]][pick,1:n.cells.y] <<- mvSaved["avail.y",1][pick,1:n.cells.y]
              model[["use.denom"]][pick] <<- mvSaved["use.denom",1][pick]
              for(k in 1:K){
                y.idx <- pick+(k-1)*M
                model$calculate(y.nodes[y.idx])
              }
              model$calculate(N.node)
            }
          }
        }
      }
    }
    #copy back to mvSaved once to synchronize all calcNode values/logProbs
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods = list(reset=function(){})
)

#Latent-ID update.
#proposal weights are computed from the factored use representation
#Pairwise match-score likelihoods are included in each ID proposal
IDSampler <- nimbleFunction( 
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    K <- control$K
    n.samples <- control$n.samples
    this.k <- control$this.k
    u.obs.cell <- control$u.obs.cell
    y.nodes <- control$y.nodes
    u.nodes <- control$u.nodes
    map.score <- control$map.score
    calcNodes.match.type <- control$calcNodes.match.type
    calcNodes.scores <- control$calcNodes.scores
    calcNodes.all <- control$calcNodes.all
    n.cells.x <- control$n.cells.x
  },
  run = function() {
    z.on <- model$z==1
    propprobs <- rep(0,M)
    for(l in 1:n.samples){
      this.cell <- u.obs.cell[l]
      cell.x <- this.cell%%n.cells.x
      cell.y <- floor(this.cell/n.cells.x)+1
      if(cell.x==0){
        cell.x <- n.cells.x
        cell.y <- cell.y-1
      }
      #For a fixed sample cell, the RSF term is the same for all candidate IDs, so it cancels
      #when proposal weights are normalized. Reuse stored availability and use.denom.
      for(i in 1:M){
        propprobs[i] <- 0
        if(z.on[i] & model$use.denom[i]>0){
          propprobs[i] <- model$avail.x[i,cell.x]*model$avail.y[i,cell.y]/model$use.denom[i]
        }
      }
      prop.sum <- sum(propprobs)
      if(prop.sum>0){
        propprobs <- propprobs/prop.sum
        pick <- rcat(1,prob=propprobs)
        if(model$ID[l]!=pick){
          swapped <- c(mvSaved["ID",1][l],pick) #order swap.out then swap.in
          propprob <- propprobs[swapped[2]]
          backprob <- propprobs[swapped[1]]
          ID.curr <- swapped[1]
          ID.cand <- swapped[2]
          k.curr <- this.k[l]
          y.idx.curr <- ID.curr+(k.curr-1)*M
          y.idx.cand <- ID.cand+(k.curr-1)*M
          these.nodes <- c(y.idx.curr,y.idx.cand)
          these.scores <- map.score[l,1:(n.samples-1)] #all pairwise scores involving focal sample l

          lp.initial.y <- model$getLogProb(y.nodes[these.nodes])
          lp.initial.u <- model$getLogProb(u.nodes[these.nodes])
          lp.initial.scores <- model$getLogProb(calcNodes.scores[these.scores])
          lp.initial.total <- lp.initial.y+lp.initial.u+lp.initial.scores

          model[["ID"]][l] <<- pick

          #select-sample proposal probs
          focalprob <- sum(mvSaved["ID",1][1:n.samples]==mvSaved["ID",1][l] & this.k==k.curr)/n.samples
          focalbackprob <- sum(model[["ID"]][1:n.samples]==model[["ID"]][l] & this.k==k.curr)/n.samples

          #If the proposed individual already has a sample on this occasion, swap the two IDs
          #so each individual can have at most one Bernoulli detection per occasion.
          check <- sum(mvSaved["ID",1][1:n.samples]==pick & this.k==k.curr)
          if(check==1){
            check.l <- which(mvSaved["ID",1][1:n.samples]==pick & this.k==k.curr)[1]
            these.scores.check.l <- map.score[check.l,1:(n.samples-1)]
            #The focal/check.l score node occurs in both score-node sets.
            #Its match status remains "different ID" before and after the ID swap,
            #so its duplicated likelihood contribution is identical on both sides of the MH ratio and cancels.
            lp.initial.total <- lp.initial.total + model$getLogProb(calcNodes.scores[these.scores.check.l])
            model[["ID"]][check.l] <<- mvSaved["ID",1][l]
          }

          #Swap the individual-occasion detection states.
          model[["y.true"]][ID.curr,k.curr] <<- mvSaved["y.true",1][ID.cand,k.curr]
          model[["y.true"]][ID.cand,k.curr] <<- mvSaved["y.true",1][ID.curr,k.curr]

          #Swap absolute and sparse surveyed-cell indices.
          model[["u.cell"]][ID.curr,k.curr] <<- mvSaved["u.cell",1][ID.cand,k.curr]
          model[["u.cell"]][ID.cand,k.curr] <<- mvSaved["u.cell",1][ID.curr,k.curr]
          model[["u.cell.survey"]][ID.curr,k.curr] <<- mvSaved["u.cell.survey",1][ID.cand,k.curr]
          model[["u.cell.survey"]][ID.cand,k.curr] <<- mvSaved["u.cell.survey",1][ID.curr,k.curr]

          #Swap continuous locations. Nondetection slots contain c(0,0).
          model[["u"]][ID.curr,k.curr,1] <<- mvSaved["u",1][ID.cand,k.curr,1]
          model[["u"]][ID.cand,k.curr,1] <<- mvSaved["u",1][ID.curr,k.curr,1]
          model[["u"]][ID.curr,k.curr,2] <<- mvSaved["u",1][ID.cand,k.curr,2]
          model[["u"]][ID.cand,k.curr,2] <<- mvSaved["u",1][ID.curr,k.curr,2]

          # match.type must be updated after changing ID and before evaluating score likelihoods.
          model$calculate(calcNodes.match.type[these.scores])
          lp.proposed.scores <- model$calculate(calcNodes.scores[these.scores])
          if(check==1){
            model$calculate(calcNodes.match.type[these.scores.check.l])
            lp.proposed.scores <- lp.proposed.scores + model$calculate(calcNodes.scores[these.scores.check.l])
          }

          lp.proposed.y <- model$calculate(y.nodes[these.nodes])
          lp.proposed.u <- model$calculate(u.nodes[these.nodes])
          lp.proposed.total <- lp.proposed.y+lp.proposed.u+lp.proposed.scores

          log_MH_ratio <- (lp.proposed.total+log(backprob)+log(focalbackprob)) -
            (lp.initial.total+log(propprob)+log(focalprob))
          accept <- decide(log_MH_ratio)

          if(accept){
            mvSaved["y.true",1][ID.curr,k.curr] <<- model[["y.true"]][ID.curr,k.curr]
            mvSaved["y.true",1][ID.cand,k.curr] <<- model[["y.true"]][ID.cand,k.curr]
            mvSaved["u.cell",1][ID.curr,k.curr] <<- model[["u.cell"]][ID.curr,k.curr]
            mvSaved["u.cell",1][ID.cand,k.curr] <<- model[["u.cell"]][ID.cand,k.curr]
            mvSaved["u.cell.survey",1][ID.curr,k.curr] <<- model[["u.cell.survey"]][ID.curr,k.curr]
            mvSaved["u.cell.survey",1][ID.cand,k.curr] <<- model[["u.cell.survey"]][ID.cand,k.curr]
            mvSaved["u",1][ID.curr,k.curr,1:2] <<- model[["u"]][ID.curr,k.curr,1:2]
            mvSaved["u",1][ID.cand,k.curr,1:2] <<- model[["u"]][ID.cand,k.curr,1:2]
            mvSaved["ID",1][l] <<- model[["ID"]][l]
            if(check==1){
              mvSaved["ID",1][check.l] <<- model[["ID"]][check.l]
            }
          }else{
            model[["y.true"]][ID.curr,k.curr] <<- mvSaved["y.true",1][ID.curr,k.curr]
            model[["y.true"]][ID.cand,k.curr] <<- mvSaved["y.true",1][ID.cand,k.curr]
            model[["u.cell"]][ID.curr,k.curr] <<- mvSaved["u.cell",1][ID.curr,k.curr]
            model[["u.cell"]][ID.cand,k.curr] <<- mvSaved["u.cell",1][ID.cand,k.curr]
            model[["u.cell.survey"]][ID.curr,k.curr] <<- mvSaved["u.cell.survey",1][ID.curr,k.curr]
            model[["u.cell.survey"]][ID.cand,k.curr] <<- mvSaved["u.cell.survey",1][ID.cand,k.curr]
            model[["u"]][ID.curr,k.curr,1:2] <<- mvSaved["u",1][ID.curr,k.curr,1:2]
            model[["u"]][ID.cand,k.curr,1:2] <<- mvSaved["u",1][ID.cand,k.curr,1:2]
            model$calculate(y.nodes[these.nodes])
            model$calculate(u.nodes[these.nodes])
            model[["ID"]][l] <<- mvSaved["ID",1][l]
            if(check==1){
              model[["ID"]][check.l] <<- mvSaved["ID",1][check.l]
            }
            #restore deterministic match types and score logProbs after restoring IDs.
            model$calculate(calcNodes.match.type[these.scores])
            model$calculate(calcNodes.scores[these.scores])
            if(check==1){
              model$calculate(calcNodes.match.type[these.scores.check.l])
              model$calculate(calcNodes.scores[these.scores.check.l])
            }
          }
        }
      }
    }

    #Detected status must reflect the current latent-ID allocation before zSampler runs.
    capcounts <- Getcapcounts(y.true=model$y.true[1:M,1:K])
    n <- Getncap(capcounts=capcounts,ID=model$ID[1:n.samples])
    model$capcounts[1:M] <<- capcounts
    model$n[1] <<- n
    # mvSaved capcounts/n are synchronized by the final copy through y.true dependencies.
    #Synchronize y/u values and logProbs once after all ID proposals.
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes.all,logProb=TRUE)
  },
  methods = list(reset=function(){})
)

