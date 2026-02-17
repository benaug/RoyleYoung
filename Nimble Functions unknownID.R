dTruncNormVector <- nimbleFunction(
  run = function(x = double(2), s = double(1), sigma = double(0), u.xlim = double(2), u.ylim = double(2),
                 n.locs.ind = double(0), log = integer(0)) {
    returnType(double(0))
    if(n.locs.ind>0){
      logProb <- 0
      for(i in 1:n.locs.ind){
        # log density at point
        logpx <- dnorm(x[i,1], mean = s[1], sd = sigma, log = TRUE)
        logpy <- dnorm(x[i,2], mean = s[2], sd = sigma, log = TRUE)
        # log CDFs - need to deal with potential underflow
        lFxU <- pnorm(u.xlim[i,2], mean = s[1], sd = sigma, log.p = TRUE)
        lFxL <- pnorm(u.xlim[i,1], mean = s[1], sd = sigma, log.p = TRUE)
        lFyU <- pnorm(u.ylim[i,2], mean = s[2], sd = sigma, log.p = TRUE)
        lFyL <- pnorm(u.ylim[i,1], mean = s[2], sd = sigma, log.p = TRUE)
        # if U <= L numerically, density is undefined -> -Inf
        if (lFxU <= lFxL | lFyU <= lFyL){
          logProb <- -Inf
        }else{
          logDenX <- lFxU + log(1 - exp(lFxL - lFxU))
          logDenY <- lFyU + log(1 - exp(lFyL - lFyU))
          logProb <- logProb + (logpx - logDenX) + (logpy - logDenY)
        }
      }
    }else{
      logProb <- 0
    }
    return(logProb)
  }
)

rTruncNormVector <- nimbleFunction(
  run = function(n = integer(0), s = double(1), sigma = double(0), u.xlim = double(2), u.ylim = double(2), n.locs.ind = double(0)) {
    returnType(double(2))
    return(matrix(0,n.locs.ind,2))
  }
)

dCatVector <- nimbleFunction(
  run = function(x = double(1), use.dist = double(1), n.locs.ind = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    for(i in 1:n.locs.ind){
      logProb <- logProb + log(use.dist[x[i]])
    }
    return(logProb)
  }
)
rCatVector <- nimbleFunction(
  run = function(n = integer(0),use.dist = double(1), n.locs.ind = double(0)) {
    returnType(double(1))
    return(rep(0,n.locs.ind))
  }
)

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
                 survey.map=double(1),
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
            if(survey.map[this.cell]==0){ #cell not surveyed
              logProb.tmp[this.cell] <- log(use.dist[this.cell])
            }else{ #cell surveyed
              logProb.tmp[this.cell] <- log(use.dist[this.cell]) + log(1-p[survey.map[this.cell]])
            }
          }
          #use pos.cells only for remaining calculations
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
                 survey.map=double(1),
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
  run = function(avail.dist=double(1),InSS.cells=double(1)){
    returnType(double(1))
    n.InSS.cells <- nimDim(InSS.cells)[1]
    pos.cells <- rep(0,n.InSS.cells)
    idx <- 1
    for(c in 1:n.InSS.cells){
      #sets level of trimming used to get use distribution and calculate y marginal logprobs
      if(avail.dist[InSS.cells[c]]>1e-5){ #not effectively zero
        pos.cells[idx] <- InSS.cells[c]
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

getCell <- nimbleFunction(#cell 0 not allowed in this model, but leaving in as an error check
  run = function(u = double(1),res=double(0),cells=integer(2),xlim=double(1),ylim=double(1)) {
    returnType(double(0))
    inout <- 1*(u[1]>xlim[1]&u[1]<xlim[2]&u[2]>ylim[1]&u[2]<ylim[2])
    if(inout==1){
      u.cell <- cells[trunc(u[1]/res)+1,trunc(u[2]/res)+1]
    }else{
      u.cell <- 0
    }
    return(u.cell)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    K <- nimDim(y.true)[2]
    capcounts=numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i]=sum(y.true[i,1:K])
    }
    return(capcounts)
  }
)

Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
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
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    MK <- M*K #can precalculate, or supply via control.
    n.det <- model$n[1] #need for proposal probs
    ind.detected <- model$capcounts>0
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){ #subtract
        #find all z's currently on
        # z.on <- which(model$z==1)
        z.on <- which(model$z==1 & ind.detected==0) #exclude guys with detections
        n.z.on <- length(z.on)
        if(n.z.on==0){
          reject <- TRUE
        }
        if(!reject){
          pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
          pick <- z.on[pick]
          pick.idx <- seq(pick,MK,M) #used to reference correct y nodes
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
          N.initial <- model$N[1]
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) +
            log(N.initial - n.det) - log(N.initial)
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
          N.initial <- model$N[1]
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y) + 
            log(N.initial + 1) - log(N.initial + 1 - n.det)
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

#original, slightly less efficient zSampler
# zSampler <- nimbleFunction(
#   contains = sampler_BASE,
#   setup = function(model, mvSaved, target, control) {
#     M <- control$M
#     K <- control$K
#     z.ups <- control$z.ups
#     y.nodes <- control$y.nodes
#     N.node <- control$N.node
#     z.nodes <- control$z.nodes
#     calcNodes <- control$calcNodes
#   },
#   run = function() {
#     MK <- M*K #can precalculate, or supply via control.
#     for(up in 1:z.ups){ #how many updates per iteration?
#       #propose to add/subtract 1
#       updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
#       reject <- FALSE #we auto reject if you select a detected call
#       if(updown==0){#subtract
#         #find all z's currently on
#         z.on <- which(model$z==1)
#         n.z.on <- length(z.on)
#         pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
#         pick <- z.on[pick]
#         #prereject turning off individuals currently allocated samples
#         if(model$capcounts[pick]>0){#is this an individual with samples?
#           reject <- TRUE
#         }
#         if(!reject){
#           pick.idx <- seq(pick,MK,M) #used to reference correct y nodes
#           #get initial logprobs for N and y
#           lp.initial.N <- model$getLogProb(N.node)
#           lp.initial.y <- model$getLogProb(y.nodes[pick.idx])
#           
#           #propose new N/z
#           model$N[1] <<-  model$N[1] - 1
#           model$z[pick] <<- 0
#           
#           #get proposed logprobs for N and y
#           lp.proposed.N <- model$calculate(N.node)
#           lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0
#           
#           #MH step
#           log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
#           accept <- decide(log_MH_ratio)
#           
#           if(accept) {
#             mvSaved["N",1][1] <<- model[["N"]]
#             mvSaved["z",1][pick] <<- model[["z"]][pick]
#           }else{
#             model[["N"]] <<- mvSaved["N",1][1]
#             model[["z"]][pick] <<- mvSaved["z",1][pick]
#             model$calculate(y.nodes[pick.idx])
#             model$calculate(N.node)
#           }
#         }
#       }else{#add
#         if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
#           z.off <- which(model$z==0)
#           n.z.off <- length(z.off)
#           pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
#           pick <- z.off[pick]
#           pick.idx <- seq(pick,MK,M) #used to reference correct y nodes
#           #get initial logprobs for N and y
#           lp.initial.N <- model$getLogProb(N.node)
#           lp.initial.y <- model$getLogProb(y.nodes[pick.idx]) #will always be 0
#           
#           #propose new N/z
#           model$N[1] <<-  model$N[1] + 1
#           model$z[pick] <<- 1
#           
#           #get proposed logprobs for N and y
#           lp.proposed.N <- model$calculate(N.node)
#           lp.proposed.y <- model$calculate(y.nodes[pick.idx])
#           
#           #MH step
#           log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
#           accept <- decide(log_MH_ratio)
#           if(accept) {
#             mvSaved["N",1][1] <<- model[["N"]]
#             mvSaved["z",1][pick] <<- model[["z"]][pick]
#           }else{
#             model[["N"]] <<- mvSaved["N",1][1]
#             model[["z"]][pick] <<- mvSaved["z",1][pick]
#             model$calculate(y.nodes[pick.idx])
#             model$calculate(N.node)
#           }
#         }
#       }
#     }
#     #copy back to mySaved to update logProbs which was not done above
#     copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
#   },
#   methods = list( reset = function () {} )
# )

IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <-control$M
    K <- control$K
    n.samples <- control$n.samples
    this.k <- control$this.k
    # u.obs <- control$u.obs
    u.obs.cell <- control$u.obs.cell
    map <- control$map
    calcNodes.y.true <- control$calcNodes.y.true
    calcNodes.u <- control$calcNodes.u
    calcNodes.all <- control$calcNodes.all
  },
  run = function() {
    s <- model$s
    z <- model$z
    sigma <- model$sigma[1]
    z.on <- z==1 #precalculate
    for(l in 1:n.samples){ #propose for each sample l one at a time
      # propprobs <- rep(0,M)
      # for(i in 1:M){
      #   if(z.on[i]){#only need to consider inds in population
      #     dist <- sqrt((s[i,1]-u.obs[l,1])^2+(s[i,2]-u.obs[l,2])^2)
      #     propprobs[i] <- exp(-dist^2/(2*sigma^2))
      #   }
      # }
      # propprobs <- propprobs/sum(propprobs)
      #this is faster
      propprobs <- model$use.dist[,u.obs.cell[l]]*z.on
      propprobs <- propprobs/sum(propprobs)
      pick <- rcat(1,prob=propprobs)
      if(model$ID[l]!=pick){ #skip if propose same ID.
        swapped <- c(mvSaved["ID",1][l],pick) #order swap.out then swap.in
        propprob <- propprobs[swapped[2]]
        backprob <- propprobs[swapped[1]]
        these.nodes <- map[swapped,this.k[l]]
        lp.initial.y.true <- model$getLogProb(calcNodes.y.true[these.nodes])
        lp.initial.u <- model$getLogProb(calcNodes.u[these.nodes])
        lp.initial.total <- lp.initial.y.true + lp.initial.u #+ lp.initial.G.obs
        model[["ID"]][l] <<- pick
        #p(select sample from this individual on this occasion)
        focalprob <- sum(mvSaved["ID",1]==mvSaved["ID",1][l]&this.k==this.k[l])/n.samples
        focalbackprob <- sum(model[["ID"]]==model[["ID"]][l]&this.k==this.k[l])/n.samples
        #check if this guy already has sample assigned on this occasion
        #it is associated with another u.obs
        #don't do this until after computing proposal probs
        check <- sum(mvSaved["ID",1]==pick&this.k==this.k[l]) #either 0 or 1
        if(check==1){#if so, swap IDs
          check.l <- which(mvSaved["ID",1]==pick&this.k==this.k[l])[1]
          model[["ID"]][check.l] <<- mvSaved["ID",1][l]
        }
        #forgive me, but this is just easier for me to track here
        ID.curr <- swapped[1]
        ID.cand <- swapped[2]
        
        #update y.true - swap em
        model[["y.true"]][ID.curr,this.k[l]] <<- mvSaved["y.true",1][ID.cand,this.k[l]]
        model[["y.true"]][ID.cand,this.k[l]] <<- mvSaved["y.true",1][ID.curr,this.k[l]]
        
        #update u.cell - swap em
        model[["u.cell"]][ID.curr,this.k[l]] <<- mvSaved["u.cell",1][ID.cand,this.k[l]]
        model[["u.cell"]][ID.cand,this.k[l]] <<- mvSaved["u.cell",1][ID.curr,this.k[l]]
        
        #update u - swap em
        model[["u"]][ID.curr,this.k[l],1] <<- mvSaved["u",1][ID.cand,this.k[l],1]
        model[["u"]][ID.cand,this.k[l],1] <<- mvSaved["u",1][ID.curr,this.k[l],1]
        model[["u"]][ID.curr,this.k[l],2] <<- mvSaved["u",1][ID.cand,this.k[l],2]
        model[["u"]][ID.cand,this.k[l],2] <<- mvSaved["u",1][ID.curr,this.k[l],2]
        
        lp.proposed.y.true <- model$calculate(calcNodes.y.true[these.nodes])
        lp.proposed.u <- model$calculate(calcNodes.u[these.nodes])
        lp.proposed.total <- lp.proposed.y.true + lp.proposed.u
        
        #MH step
        log_MH_ratio <- (lp.proposed.total + log(backprob) + log(focalbackprob)) -
          (lp.initial.total + log(propprob) + log(focalprob))
        accept <- decide(log_MH_ratio)
        if(accept){
          #copy model to mvSaved - dont need to keep up with logprobs
          mvSaved["y.true",1][ID.curr,this.k[l]] <<- model[["y.true"]][ID.curr,this.k[l]]
          mvSaved["y.true",1][ID.cand,this.k[l]] <<- model[["y.true"]][ID.cand,this.k[l]]
          mvSaved["u.cell",1][ID.curr,this.k[l]] <<- model[["u.cell"]][ID.curr,this.k[l]]
          mvSaved["u.cell",1][ID.cand,this.k[l]] <<- model[["u.cell"]][ID.cand,this.k[l]]
          mvSaved["u",1][ID.curr,this.k[l],1] <<- model[["u"]][ID.curr,this.k[l],1]
          mvSaved["u",1][ID.cand,this.k[l],1] <<- model[["u"]][ID.cand,this.k[l],1]
          mvSaved["u",1][ID.curr,this.k[l],2] <<- model[["u"]][ID.curr,this.k[l],2]
          mvSaved["u",1][ID.cand,this.k[l],2] <<- model[["u"]][ID.cand,this.k[l],2]
          mvSaved["ID",1][l] <<- model[["ID"]][l]
          if(check==1){
            mvSaved["ID",1][check.l] <<- model[["ID"]][check.l]
          }
        }else{
          #set model back to mvSaved states
          model[["y.true"]][ID.curr,this.k[l]] <<- mvSaved["y.true",1][ID.curr,this.k[l]]
          model[["y.true"]][ID.cand,this.k[l]] <<- mvSaved["y.true",1][ID.cand,this.k[l]]
          model[["u.cell"]][ID.curr,this.k[l]] <<- mvSaved["u.cell",1][ID.curr,this.k[l]]
          model[["u.cell"]][ID.cand,this.k[l]] <<- mvSaved["u.cell",1][ID.cand,this.k[l]]
          model[["u"]][ID.curr,this.k[l],1] <<- mvSaved["u",1][ID.curr,this.k[l],1]
          model[["u"]][ID.cand,this.k[l],1] <<- mvSaved["u",1][ID.cand,this.k[l],1]
          model[["u"]][ID.curr,this.k[l],2] <<- mvSaved["u",1][ID.curr,this.k[l],2]
          model[["u"]][ID.cand,this.k[l],2] <<- mvSaved["u",1][ID.cand,this.k[l],2]
          #need to set logProbs back, not updated in mvSaved, so recalculate. not most efficient
          model$calculate(calcNodes.y.true[these.nodes])
          model$calculate(calcNodes.u[these.nodes])
          model[["ID"]][l] <<- mvSaved["ID",1][l]
          if(check==1){
            model[["ID"]][check.l] <<- mvSaved["ID",1][check.l]
          }
        }
      }
    }
    capcounts <- Getcapcounts(y.true=model$y.true[1:M,1:K])
    n <- Getncap(capcounts=capcounts,ID=model$ID)
    model$capcounts <<- capcounts
    model$n[1] <<- n
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes.all, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
