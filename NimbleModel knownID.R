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
  #--------------------------------------------------------------
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
  rsf[1:n.cells] <- exp(rsf.beta*rsf.cov[1:n.cells])
  #Detection model
  for(k in 1:K){
    #p in each surveyed cell on occasion k - sparse matrix representation
    p[1:n.surveyed.cells[k],k] <- getP(surveyed.cells.effort=surveyed.cells.effort[1:n.surveyed.cells[k],k],
                                      n.surveyed.cells=n.surveyed.cells[k],
                                      beta.p.int=beta.p.int,beta.p.effort=beta.p.effort)
  }
  for(i in 1:M){
    #continuous activity center likelihood inside cell
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1] #extract activity center cell
    #categorical activity center likelihood for this cell, equivalent to zero's trick
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
    #Individual availability distributions conditioned on cells, bivariate Normal centered on activity center
    avail.dist[i,1:n.cells] <- getAvail(s=s[i,1:2],sigma=sigma,res=res,x.vals=x.vals[1:n.cells.x],
                                        y.vals=y.vals[1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    #dynamic sparse matrix representation of cells with nonzero prob of availability to trim 
    #use dist and marginalization calculations below, only needs to be recomputed when s_i or sigma update
    pos.cells[i,1:n.cells] <- getPosCells(avail.dist=avail.dist[i,1:n.cells],n.cells=n.cells)
    n.pos.cells[i] <- getNPosCells(pos.cells=pos.cells[i,1:n.cells])
    #Individual use distributions - multiply rsf and available distribution, normalize. trimmed
    use.dist[i,1:n.cells] <- getUse(rsf=rsf[1:n.cells],avail=avail.dist[i,1:n.cells],
                                    pos.cells=pos.cells[i,1:n.cells],n.pos.cells=n.pos.cells[i],n.cells=n.cells)
    for(k in 1:K){
      #observation likelihood
      #for detections, use detection likelihood conditional on cell of detection
      #for nondetections, marginalize detection over unobserved site use, u.cell, marginal prob(not detected). trimmed
      y[i,k] ~ dRYmarg(u.cell=u.cell[i,k],sigma=sigma,z=z[i],p=p[1:n.surveyed.cells[k],k],survey=survey[1:n.cells,k],
                            surveyed.cells=surveyed.cells[1:n.surveyed.cells[k],k],n.surveyed.cells=n.surveyed.cells[k],
                            use.dist=use.dist[i,1:n.cells],
                            pos.cells=pos.cells[i,1:n.cells],n.pos.cells=n.pos.cells[i],n.cells=n.cells,
                            res=res)
      #continuous use location likelihood conditioned on the cell of detection
      #split out of likelihood above bc not used when updating rsf beta, all parameters on p
      #logprob is 0 for all nondetects with u.cell[i,k]=0, which includes all k for z[i]=0 inds.
      u[i,k,1:2] ~ duInCell(s=s[i,1:2],u.cell=u.cell[i,k],sigma=sigma,n.cells.x=n.cells.x,res=res)
    }
  }
  #telemetry - could use these for mark-resight if random sample w.r.t. space or marking process is modeled
  for(i in 1:n.tel.inds){
    s.tel[i,1] ~ dunif(xlim[1],xlim[2])
    s.tel[i,2] ~ dunif(ylim[1],ylim[2])
    #can use telemetry data to inform D cov estimation, assumes inds captured at random wrt to expected D
    s.cell.tel[i] <- cells[trunc(s.tel[i,1]/res)+1,trunc(s.tel[i,2]/res)+1] #extract activity center cell
    dummy.data.tel[i] ~ dCell(pi.cell[s.cell.tel[i]])
    #Individual available distributions - bivariate Normal centered on activity center.
    avail.dist.tel[i,1:n.cells] <- getAvail(s=s.tel[i,1:2],sigma=sigma,res=res,x.vals=x.vals[1:n.cells.x],
                                            y.vals=y.vals[1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y)
    #Individual use distributions - multiply rsf and available distribution, normalize. not trimmed
    use.dist.tel[i,1:n.cells] <- getUseTel(rsf=rsf[1:n.cells],avail=avail.dist.tel[i,1:n.cells])
    for(k in 1:n.locs.ind[i]){
      u.cell.tel[i,k] ~ dcat(use.dist.tel[i,1:n.cells]) #likelihood of this cell being used
      #continuous use location likelihood conditioned on the cell
      u.tel[i,k,1] ~ T(dnorm(s.tel[i,1],sd=sigma),u.xlim.tel[i,k,1],u.xlim.tel[i,k,2])
      u.tel[i,k,2] ~ T(dnorm(s.tel[i,2],sd=sigma),u.ylim.tel[i,k,1],u.ylim.tel[i,k,2])
    }
  }
})# end model
