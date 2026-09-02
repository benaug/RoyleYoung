NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  #RSF coefficients
  rsf.beta ~ dnorm(0,sd=10) 
  #availability distribution spatial scale
  sigma ~ dunif(0,20)
  #detection coefficients
  beta.p.int ~ dlogis(0,1)
  beta.p.effort ~ dnorm(0,sd=10)
  #pairwise match-score means
  for(q in 1:2){
    lambda.match[q] ~ dunif(0,50) # q=1 same ID, q=2 different ID
  }
  
  #Density model
  D.intercept <- D0*cellArea
  # D.intercept <- exp(D.beta0)*cellArea
  #multiplying by InSS=0 prevents activity centers from living there, InSS=1 otherwise
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells]) #separate this component so s's do not depend on D.intercept
  # lambda.cell[1:n.cells] <- InSS[1:n.cells] #if no Dcov, this gives homogeneous D
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda)

  #Resource selection function evaluated across all cells
  rsf[1:n.cells] <- InSS[1:n.cells]*exp(rsf.beta*rsf.cov[1:n.cells])

  #Detection model
  for(k in 1:K){
    #p in each surveyed cell on occasion k - sparse matrix representation
    p[1:n.surveyed.cells[k],k] <- getP(surveyed.cells.effort=surveyed.cells.effort[1:n.surveyed.cells[k],k],
                                       n.surveyed.cells=n.surveyed.cells[k],
                                       beta.p.int=beta.p.int,beta.p.effort=beta.p.effort)
    #save p*RSF for surveyed cells because it is shared across all individuals.
    for(c in 1:n.surveyed.cells[k]){
      p.rsf[c,k] <- p[c,k]*rsf[surveyed.cells[c,k]]
    }
  }
  for(i in 1:M){
    # z-gated AC distribution. When z=0, s is fixed at c(0,0); when z=1, s is drawn
    # from pi.cell and then uniformly within the selected cell.
    s[i,1:2] ~ dAC(pi.cell=pi.cell[1:n.cells],res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=z[i])
    #Factored individual availability distribution. The BVN cell probability is
    #avail.x[cell.x]*avail.y[cell.y], no n.cells-length availability vector is stored.
    avail.x[i,1:n.cells.x] <- getAvail1D(s=s[i,1],sigma=sigma,res=res,vals.edges=x.vals.edges[1:(n.cells.x+1)],
                                         n.cells=n.cells.x,avail.z=avail.z,z=z[i])
    avail.y[i,1:n.cells.y] <- getAvail1D(s=s[i,2],sigma=sigma,res=res,vals.edges=y.vals.edges[1:(n.cells.y+1)],
                                         n.cells=n.cells.y,avail.z=avail.z,z=z[i])
    #RSF-weighted normalizing constant
    use.denom[i] <- getUseDenom(rsf=rsf[1:n.cells],avail.x=avail.x[i,1:n.cells.x],
                                avail.y=avail.y[i,1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=z[i])

    for(k in 1:K){
      #Latent-ID detection state; IDSampler moves detections among individual-occasion slots.
      y.true[i,k] ~ dRYmargFactored(u.cell.survey=u.cell.survey[i,k],z=z[i],
                                    p.rsf=p.rsf[1:n.surveyed.cells[k],k],
                                    surveyed.cell.x=surveyed.cell.x[1:n.surveyed.cells[k],k], 
                                    surveyed.cell.y=surveyed.cell.y[1:n.surveyed.cells[k],k], 
                                    n.surveyed.cells=n.surveyed.cells[k], 
                                    avail.x=avail.x[i,1:n.cells.x],avail.y=avail.y[i,1:n.cells.y], 
                                    use.denom=use.denom[i]) 
      #Latent-ID continuous location likelihood. Nondetection slots have u.cell=0 and contribute 0
      u[i,k,1:2] ~ duInCell(s=s[i,1:2],u.cell=u.cell[i,k],sigma=sigma,n.cells.x=n.cells.x,res=res,
                            avail.x=avail.x[i,1:n.cells.x],avail.y=avail.y[i,1:n.cells.y])
    }
  }

  #telemetry - could use these for mark-resight if random sample w.r.t. space or marking process is modeled
  for(i in 1:n.tel.inds){
    #use same density process as detected individuals
    # s.tel[i,1:2] ~ dAC(pi.cell=pi.cell[1:n.cells],res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
    #uniform over InSS
    s.tel[i,1:2] ~ dAC(pi.cell=pi.cell.tel[1:n.cells],res=res,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
    #Factored telemetry availability and normalizing constants.
    avail.x.tel[i,1:n.cells.x] <- getAvail1D(s=s.tel[i,1],sigma=sigma,res=res,vals.edges=x.vals.edges[1:(n.cells.x+1)],
                                             n.cells=n.cells.x,avail.z=avail.z,z=1)
    avail.y.tel[i,1:n.cells.y] <- getAvail1D(s=s.tel[i,2],sigma=sigma,res=res,vals.edges=y.vals.edges[1:(n.cells.y+1)],
                                             n.cells=n.cells.y,avail.z=avail.z,z=1)
    use.denom.tel[i] <- getUseDenom(rsf=rsf[1:n.cells],avail.x=avail.x.tel[i,1:n.cells.x],
                                    avail.y=avail.y.tel[i,1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y,z=1)
    #Cell-use probability and within-cell truncation cancel for telemetry locations.
    u.tel[i,1:n.locs.ind[i],1:2] ~ dTelemetryFactored(u.cell=u.cell.tel[i,1:n.locs.ind[i]],s=s.tel[i,1:2],
                                                      sigma=sigma,rsf=rsf[1:n.cells],use.denom=use.denom.tel[i],
                                                      n.cells.x=n.cells.x,n.locs.ind=n.locs.ind[i])
  }

  #Pairwise sample-level match scores
  #Only the upper triangle matrix is modeled because each unordered sample pair appears once.
  for(l1 in 1:(n.samples-1)){
    for(l2 in (l1+1):n.samples){
      match.type[l1,l2] <- 1*equals(ID[l1],ID[l2]) + 2*(1-equals(ID[l1],ID[l2]))
      scores[l1,l2] ~ dpois(lambda.match[match.type[l1,l2]])
    }
  }

  #Latent-ID derived variables.
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:K])
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples])
})
