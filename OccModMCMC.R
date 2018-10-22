OccModMCMC <- function(Y,J,X,y.oos,J.oos,X.oos,n.mcmc,s2.beta=1.5^2,null.model=FALSE,
                       no.print=FALSE,get.D=FALSE,no.pred=FALSE){

#  MCMC to fit Bayesian hierarchical site-occupancy model as described in Hooten & Hobbs (2015) A guide to 
#  Bayesian model selection for ecologists. Ecol. App. 85(1): 3-28. MCMC code was modified slightly for use in
#  Parallel R Workshop at the Research Centre for Ecological Change Annual Meeting, October 25, 2018 in Lammi, 
#  Finland, but is otherwise reproduced from the Supplemental Material provided by Hooten and Hobbs (2015).
#
#  Please see the original disclaimer for use of MCMC code below.
# 
# 
#  (20140406)
#  This R function implements an MCMC algorithm for fitting the Bayesian hierarchical probit occupancy model.  
#  There are various flags in the inputs that allow for faster computation by limiting extra calculations
#  for things like, status updates, prediction, and criteria calculations.
#  
#  y, J, X are input data to the model (within-sample response and covariate data).
#  y.oos, J.oos, X.oos are out-of-sample data to obtain predictions for (something always needs to be specified for these).  
#  s2.beta is the prior variance for beta (i.e., the regularlization parameter). 
#  If an intercept-only model fit is desired use the 'null.model=TRUE' flag; in this case, you'll still specify X, X.oos
#  as placeholders (they don't affect output).  
#
#  This MCMC algorithm is for demonstration only; it is by no means a general algorithm for all occupancy models.  If
#  you feel compelled to modify the code provided here, please do so while following formal rules for MCMC and 
#  statistical modeling and use at your own risk.  
#

####
####  Libraries and Subroutines
####

rtn1 <- function(mu){
  mu + qnorm(log(runif(length(mu)))+pnorm(mu,log=TRUE),lower.tail=FALSE,log.p=TRUE)
}

rtn0 <- function(mu){
  mu - qnorm(log(runif(length(mu)))+pnorm(-mu,log=TRUE),lower.tail=FALSE,log.p=TRUE)
}

####
####  Setup Variables 
####

X=as.matrix(X)
X.oos=as.matrix(X.oos,,dim(X)[2])
n.oos=length(y.oos)

n=dim(Y)[1]
q.beta=dim(X)[2]
n.burn=round(0.25*n.mcmc)
# n.save=ceiling(0.75*n.mcmc)

z.mean=rep(0,n)
beta.0.save=rep(0,n.mcmc)
beta.save=matrix(0,n.mcmc,q.beta)
p.save=rep(0,n.mcmc)
# psi.save=matrix(0,n.mcmc,nrow(X))

if(get.D){
  y.save=matrix(0,n.mcmc,n)
}

ppd.save=matrix(0,n.mcmc,n)
ppd.oos.save=matrix(0,n.mcmc,n.oos)
lppd.save=matrix(0,n.mcmc,n)

####
####  Hyperparameters and Starting Values 
####

beta.0=0
beta=rep(0,q.beta)
Xbeta=X%*%beta
psi=pnorm(beta.0+Xbeta)
p=.5
y=apply(Y,1,sum,na.rm=TRUE)
z=rep(0,n)
z[y==0]=0
n0=sum(z==0)
n1=sum(z==1)

v=rep(0,n)

s2.0=1.5^2
Sig.beta=s2.beta*diag(q.beta)
Sig.beta.inv=solve(Sig.beta)
mu.beta=rep(0,q.beta)

alpha.1=1
alpha.2=1

XprimeX=t(X)%*%X
A.beta=XprimeX+Sig.beta.inv
A.beta.chol=chol(A.beta)
Sig.beta.inv.times.mu.beta=Sig.beta.inv%*%mu.beta

####
####  Begin Gibbs Loop 
####

if(!no.print){
  bar = txtProgressBar(min=1,max=(n.mcmc+1),initial=0,style=3,char="*",width=50,title=m)
}
for(k in 1:n.mcmc){
  # if(!no.print){
  #   if(k%%100==0) cat(k," ");flush.console()
  # }

  ####
  ####  Sample v 
  ####

  v[z==0]=rtn0(beta.0+Xbeta[z==0])
  v[z==1]=rtn1(beta.0+Xbeta[z==1])
 
  ####
  ####  Sample beta 
  ####
  ####  Fancy way to sample correlated MVN rv; super fast.
  ####

  if(!null.model){
    b.beta=t(X)%*%(v-beta.0)+Sig.beta.inv.times.mu.beta
    beta=backsolve(A.beta.chol, backsolve(A.beta.chol, b.beta, transpose = TRUE) + rnorm(q.beta))
    Xbeta=X%*%beta
  }
  if(null.model){
    beta=rep(0,q.beta)
    Xbeta=rep(0,n)
  }

  ####
  ####  Sample beta.0
  ####

  tmp.var=1/(n+1/s2.0)
  tmp.mn=tmp.var*(sum(v-Xbeta))
  beta.0=rnorm(1,tmp.mn,sqrt(tmp.var))
  psi=pnorm(beta.0+Xbeta)

  ####
  ####  Sample p 
  ####
 
  p=rbeta(1,sum(y)+alpha.1,sum(J-y)+alpha.2)

  ####
  ####  Sample z 
  ####

  psi.numer=psi*(1-p)^J
  psi.tmp=psi.numer/(psi.numer+1-psi)
  psi.tmp[psi.tmp>.9999]=.9999
  psi.tmp[psi.tmp<.0001]=.0001
  z[y==0]=rbinom(sum(y==0),1,psi.tmp[y==0])
  z[y>0]=1
  n0=sum(z==0)
  n1=sum(z==1)

  ####
  ####  Save Samples 
  ####

  beta.0.save[k]=beta.0
  beta.save[k,]=beta
  p.save[k]=p
  # psi.save[k-n.burn,] = psi
  if(k > n.burn){
    z.mean=z.mean+z/(n.mcmc-n.burn)
  }
  
  ppd.save[k,]=dbinom(y,J,p*z)
  lppd.save[k,]=dbinom(y,J,p*z,log=TRUE)

  if(get.D){
    y.save[k,]=rbinom(n,J,p*z)
  }

  if(!no.pred){
    z.oos=rbinom(n.oos,1,pnorm(beta.0+X.oos%*%beta)) 
    ppd.oos.save[k,]=dbinom(y.oos,J.oos,p*z.oos)
  }
  if(!no.print){
    setTxtProgressBar(bar,k)
  }
}
# if(!no.print){
#   cat("\n");flush.console()
# }

####
####  Compute Model Comparison Metrics
####

D.hat=-2*sum(dbinom(y,J,mean(p.save[n.burn:n.mcmc])*z.mean,log=TRUE))
D.bar=mean(-2*apply(lppd.save[n.burn:n.mcmc,],1,sum))
pD=D.bar-D.hat
DIC=D.hat+2*pD
tmp.log=log(apply(ppd.save[n.burn:n.mcmc,],2,mean))
tmp.sum=-2*sum(tmp.log)
pD.1=2*sum(tmp.log-apply(lppd.save[n.burn:n.mcmc,],2,mean))
pD.2=sum(apply(lppd.save[n.burn:n.mcmc,],2,var))
WAIC.1=tmp.sum+2*pD.1
WAIC.2=tmp.sum+2*pD.2
CPO.vec=(n.mcmc-n.burn)/apply(1/ppd.save[n.burn:n.mcmc,],2,sum)
CPO=-sum(log(CPO.vec))

if(get.D){
  y.mean=apply(y.save[n.burn:n.mcmc,],2,mean)
  sum.1=sum((y-y.mean)^2)
  sum.2=sum(apply(y.save[n.burn:n.mcmc,],2,var))
  D=sum.1+sum.2
}

score=0
if(!no.pred){
  score=sum(log(apply(ppd.oos.save[n.burn:n.mcmc,],2,mean)))
}

####
####  Write Output 
####

if(!get.D){
  list(p.save=p.save[(n.burn+1):n.mcmc],beta.0.save=beta.0.save[(n.burn+1):n.mcmc],beta.save=beta.save[(n.burn+1):n.mcmc,],
       z.mean=z.mean,n.mcmc=n.mcmc,WAIC.2=WAIC.2,pD.2=pD.2,WAIC.1=WAIC.1,pD.1=pD.1,score=score,CPO.vec=CPO.vec,CPO=CPO,pD=pD,DIC=DIC,X=X,s2.beta=s2.beta,y=y,J=J)
}
else{
  list(p.save=p.save[(n.burn+1):n.mcmc],beta.0.save=beta.0.save[(n.burn+1):n.mcmc],beta.save=beta.save[(n.burn+1):n.mcmc,],
       z.mean=z.mean,n.mcmc=n.mcmc,WAIC.2=WAIC.2,pD.2=pD.2,WAIC.1=WAIC.1,pD.1=pD.1,score=score,CPO.vec=CPO.vec,CPO=CPO,pD=pD,DIC=DIC,D=D,X=X,s2.beta=s2.beta,y=y,J=J)
}

}
