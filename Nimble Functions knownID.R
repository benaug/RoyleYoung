duInCell <- nimbleFunction(
  run = function(x = double(1), s = double(1), u.cell = double(0),
                 sigma = double(0),n.cells.x = integer(0),res=double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(u.cell>0){
      u.cell.x <- u.cell%%n.cells.x
      u.cell.y <- floor(u.cell/n.cells.x)+1
      if(u.cell.x==0){
        u.cell.x <- n.cells.x
        u.cell.y <- u.cell.y-1
      }
      xlim.cell <- c(u.cell.x-1,u.cell.x)*res
      ylim.cell <- c(u.cell.y-1,u.cell.y)*res
      logProb <- log(dnorm(x[1],s[1],sigma,log=FALSE)/
                       (pnorm(xlim.cell[2],s[1],sd=sigma) - pnorm(xlim.cell[1],s[1],sd=sigma))) + #x continuous likelihood
        log(dnorm(x[2],s[2],sigma,log=FALSE)/
              (pnorm(ylim.cell[2],s[2],sd=sigma) - pnorm(ylim.cell[1],s[2],sd=sigma))) #y continuous likelihood
    }else{
      logProb <- 0
    }
    return(logProb)
  }
)

#dummy RNG to make nimble happy
ruInCell <- nimbleFunction(
  run = function(n = integer(0), s = double(1), u.cell = double(0),
                 sigma = double(0),n.cells.x = integer(0),res=double(0)) {
    returnType(double(1))
    return(c(0,0))
  })

dRYmarg <- nimbleFunction(
  run = function(x = double(0),u.cell=double(0),z = integer(0), p = double(1),
                 sigma = double(0),use.dist = double(1), survey = double(1),
                 surveyed.cells = double(1), n.surveyed.cells = integer(0),
                 pos.cells = double(1),n.pos.cells = integer(0),n.cells = integer(0),
                 res=double(0), log = integer(0)) {
    returnType(double(0))
    if(z==1){
      #calculate occasion-level likelihood marginalized over u
      if(x==1){#If we observed a u...
        idx <- which(surveyed.cells==u.cell)[1]
        this.p <- p[idx]
        #obsmod
        # logProb <- dbinom(1,1,this.p,log=TRUE) + 
        logProb <- log(this.p) + #log(p) faster than dbinom
          log(use.dist[u.cell]) #RSF cell use likelihood
      }else{#unobserved u's
        overlap <- sum(survey[pos.cells[1:n.pos.cells]])>0
        if(overlap){ #if there is non-negligible overlap between this individuals use cells and survey effort on this occasion
          #only sum over relevant pos cells
          logProb.tmp <- rep(-Inf,n.cells)
          for(c in 1:n.pos.cells){ #start with use component
            this.cell <- pos.cells[c]
            logProb.tmp[this.cell] <- log(use.dist[this.cell]) #inefficient to be logging the use distribution for every k
          }
          for(c in 1:n.surveyed.cells){ #add detection component
            this.cell <- surveyed.cells[c]
            logProb.tmp[this.cell] <- logProb.tmp[this.cell] + log(1-p[c])
          }
          #keep use pos.cells only for remaining calculations
          logProb.tmp.reduced <- rep(0,n.pos.cells)
          for(c in 1:n.pos.cells){
            logProb.tmp.reduced[c] <- logProb.tmp[pos.cells[c]]
          }
          maxlp <- max(logProb.tmp.reduced)
          logProb <- log(sum(exp(logProb.tmp.reduced-maxlp)))+maxlp
        }else{ #if no overlap of individual site use and survey effort on this occasion
          logProb <- 0
        }
      }
    }else{#if z=0, we just integrate the use distribution giving us a logProb of 0
      logProb <- 0
    }
    return(logProb)
  }
)

#dummy RNG to make nimble happy
rRYmarg <- nimbleFunction(
  run = function(n = integer(0),u.cell=double(0),z = integer(0), p = double(1),
                 sigma = double(0), use.dist = double(1), survey = double(1),
                 surveyed.cells = double(1), n.surveyed.cells = integer(0),
                 pos.cells = double(1),n.pos.cells = integer(0),n.cells = integer(0),
                 res=double(0)) {
    returnType(double(0))
    return(1)
  })

getAvail <- nimbleFunction(
  run = function(s = double(1),sigma=double(0),res=double(0),x.vals=double(1),y.vals=double(1),n.cells.x=integer(0),n.cells.y=integer(0)) {
    returnType(double(1))
    avail.dist.x <- rep(0,n.cells.x)
    avail.dist.y <- rep(0,n.cells.y)
    delta <- 1e-8 #this sets the degree of trimming used to get individual availability distributions
    x.limits <- qnorm(c(delta,1-delta),mean=s[1],sd=sigma)
    y.limits <- qnorm(c(delta,1-delta),mean=s[2],sd=sigma)
    #convert to grid edges instead of centroids
    x.vals.edges <- c(x.vals - res/2, x.vals[n.cells.x]+0.5*res)
    y.vals.edges <- c(y.vals - res/2, y.vals[n.cells.y]+0.5*res)
    #trim in x and y direction
    if(x.vals.edges[1]<x.limits[1]){
      x.start <- sum(x.vals.edges<x.limits[1])
    }else{
      x.start <- 1
    }
    if(x.vals.edges[n.cells.x]>x.limits[2]){
      x.stop <- which(x.vals.edges>x.limits[2])[1]
    }else{
      x.stop <- n.cells.x
    }
    if(y.vals.edges[1]<y.limits[1]){
      y.start <- sum(y.vals.edges<y.limits[1])
    }else{
      y.start <- 1
    }
    if(y.vals.edges[n.cells.y]>y.limits[2]){
      y.stop <- which(y.vals.edges>y.limits[2])[1]
    }else{
      y.stop <- n.cells.y
    }
    #get pnorms
    pnorm.x <- rep(0,n.cells.x+1)
    pnorm.y <- rep(0,n.cells.y+1)
    pnorm.x[x.start:(x.stop+1)] <- pnorm(x.vals.edges[x.start:(x.stop+1)], mean=s[1], sd=sigma)
    pnorm.y[y.start:(y.stop+1)] <- pnorm(y.vals.edges[y.start:(y.stop+1)], mean=s[2], sd=sigma)
    # Compute availability distributions
    avail.dist.x[x.start:x.stop] <- pnorm.x[(x.start+1):(x.stop+1)] - pnorm.x[x.start:x.stop]
    avail.dist.y[y.start:y.stop] <- pnorm.y[(y.start+1):(y.stop+1)] - pnorm.y[y.start:y.stop]
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

getUse <- nimbleFunction(
  run = function(rsf = double(1),avail.dist=double(1),pos.cells=double(1),n.pos.cells=double(0),n.cells=double(0)){
    returnType(double(1))
    use.dist <- rep(0,n.cells)
    sum.dist <- 0
    for(c in 1:n.pos.cells){
      this.cell <- pos.cells[c]
      use.dist[this.cell] <- rsf[this.cell]*avail.dist[this.cell]
      sum.dist <- sum.dist + use.dist[this.cell]
    }
    for(c in 1:n.pos.cells){
      this.cell <- pos.cells[c]
      use.dist[this.cell] <- use.dist[this.cell]/sum.dist
    }
    return(use.dist)
  }
)

getUseTel <- nimbleFunction(
  run = function(rsf = double(1),avail.dist=double(1)) {
    returnType(double(1))
    use.dist <- rsf*avail.dist
    use.dist <- use.dist/sum(use.dist)
    return(use.dist)
  }
)

getPosCells <- nimbleFunction(
  run = function(avail.dist=double(1),n.cells=integer(0)) {
    returnType(double(1))
    pos.cells <- rep(0,n.cells)
    idx <- 1
    for(c in 1:n.cells){
      #sets level of trimming used to get use distribution and calculate y marginal logprobs
      if(avail.dist[c]>1e-5){ #not effectively zero
        pos.cells[idx] <- c
        idx <- idx + 1
      }
    }
    return(pos.cells)
  }
)

getNPosCells <- nimbleFunction(
  run = function(pos.cells=double(1)) {
    returnType(integer(0))
    n.pos.cells <- sum(pos.cells>0)
    return(n.pos.cells)
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

dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    K <- control$K
    ind.detected <- control$ind.detected
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    MK <- M*K #can precalculate, or supply via control.
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(ind.detected[pick]==1){#is this an individual with captured?
          reject <- TRUE
        }
        if(!reject){
          pick.idx <- seq(pick,MK,M) #used to reference correct y nodes
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick.idx])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          pick.idx <- seq(pick,MK,M) #used to reference correct y nodes
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick.idx])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)