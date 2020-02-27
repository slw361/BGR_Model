get.TL.distmats <- function(L){
    ## L is a list, with each entry C1, C2, ... being a log-connectivity matrix
    ## 
    ## C[i,j] = directional covariate from i -> j
    ## Note that C[i,j] = 0 typically if i and j are not neighbors
    p=length(L)
    M=ncell(rast.stack[[1]])
    B=list()
    for(i in 1:p){
        B[[i]]=Matrix(0,nrow=M,ncol=M)
    }
    n.x=ncol(rast.stack[[1]])
    n.y=nrow(rast.stack[[1]])
    N=n.x*n.y

    xy=rowColFromCell(rast.stack[[1]],1:M)

    ## left cell
    idx=which(xy[,2]>1)
    left.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1],xy[idx,2]-1)
    B.idx=idx+N*(left.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],
            left.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],
             j=left.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }

    ## right cell
    idx=which(xy[,2]<n.x)
    right.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1],xy[idx,2]+1)
    B.idx=idx+N*(right.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],right.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],j=right.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }

    ## up cell
    idx=which(xy[,1]>1)
    up.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1]-1,xy[idx,2])
    B.idx=idx+N*(up.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],up.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],j=up.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }

    ## down cell
    idx=which(xy[,1]<n.y)
    down.cells=cellFromRowCol(rast.stack[[1]],xy[idx,1]+1,xy[idx,2])
    B.idx=idx+N*(down.cells-1)
    for(i in 1:p){
        vals=.5*(extract(rast.stack[[i]],idx)+extract(rast.stack[[i]],down.cells))
        idx.nonzero=which(vals!=0)
        B[[i]]=B[[i]]+sparseMatrix(i=idx[idx.nonzero],j=down.cells[idx.nonzero],x=vals[idx.nonzero],dims=c(N,N))
    }
    B
}


get.Q <-
function(TL,beta=0,model="ICAR"){
    p=length(TL)
    M=nrow(TL[[1]])
    A=Matrix(0,nrow=M,ncol=M,sparse=T)
    for(i in 1:p){
        A=A+beta[i]*TL[[i]]
    }
    A@x=exp(A@x)
    m=rowSums(A)
    if(model=="ICAR"){
        Q=Diagonal(M,m)-A
        Q=1/2*(Q+t(Q))
    }
    if(model=="SAR"){
        Q=Diagonal(M,m)-A
        Q=Q%*%t(Q)
    }
    rm(A)
    Q
}


mcmc.wish <- function(Dobs, TL, obs.idx, df=1,model="ICAR",
                           beta.start = rep(0, length(TL)),
                           beta.prior.mean = rep(0, length(TL)),
                           beta.prior.cov = diag(10, length(TL)),
                           tau.start = 0.1, tau.prior.var = 1,
                           theta.tune = diag(10^-4,length(TL)+1),
                           n.mcmc=100, adapt.max=10000, adapt.int=100,
                           print.iter=FALSE, output.trace.plot=FALSE){


    ## ## wrapper function to create covariance matrix of observations
    get.Psi=function(TL,beta,tau,cells=cells,K=K,L=L,cov.model=model){
        ## get precision matrix of whole graph
        Q=get.Q(TL,beta,model=cov.model)
        ## get precision matrix of observed nodes
        max.diag=max(diag(Q))
        Q=Q/max.diag
        Phi=get.Phi(Q,cells)
        ## get covariance matrix of observations
        Sigma.nodes=ginv(as.matrix(Phi))
        Sigma.nodes=Sigma.nodes/max.diag
        Psi=K%*%Sigma.nodes%*%t(K)+tau*diag(nrow(K))
        Psi
    }

    ## Preliminaries----------------------------------------
    ## make L and W=L'(-Dobs)L ~ Wishart(df,L'(2Sigma)L)
    n=nrow(Dobs)
    L=diag(n-1)
    L=cbind(L,-1)
    L=t(L)
    W=t(L)%*%(-Dobs)%*%L

    ## number of params
    p <- length(TL)
    ## number of observations
    n.obs <- length(obs.idx)
    ## unique indexes of observations
    cells.idx <- unique(obs.idx)
    n.cells <- length(cells.idx)
    ## matrix for creating psi covariance matrix (above equation 11)
    K <- matrix(0, nrow = n.obs, ncol = n.cells)
    for (i in 1:n.obs){
        K[i, which(cells.idx == obs.idx[i])] <- 1
    }

    ## starting values
    beta=beta.start
    tau=tau.start
    logtau=log(tau)

    theta=c(logtau,beta)

    ## log-adaptive proposal stuff
    Sig.hat=theta.tune
    s2.tune=2.4^2/length(theta)
    theta.tune=s2.tune*Sig.hat

    ## Starting value for ll
    Psi=get.Psi(TL,beta,tau,cells.idx,K,L)
    ll=dGenWish(Dobs,Psi,df,log=TRUE)
    ## MCMC ---------------------------------------------------
    beta.save=matrix(NA,n.mcmc,p)
    tau.save=rep(NA,n.mcmc)
    ll.save=rep(NA,n.mcmc)
    logprior.save=rep(NA,n.mcmc)
    dic.sum=0
    accept=0

    for(iter in 1:n.mcmc){
        if(print.iter)    cat(iter," ",beta,tau,"\n")
        if(iter%%100==0) cat (iter," ",beta,tau,"\n")
        if(iter%%1000==0 & output.trace.plot){
            pdf("traceOut.pdf",width=12,height=5)
            matplot(beta.save[1:(iter),],type="l")
            legend("topleft",legend=1:length(beta),lwd=1,col=1:length(beta),bg="white")
            abline(h=0,col="yellow",lwd=2)
            dev.off()
        }

        ## propose
        pd=0
        while(pd==0){
            theta.star=rmvnorm(1,theta,theta.tune)
            beta.star=theta.star[-1]
            tau.star=exp(theta.star[1])

            ## get ll for gen wish
            Psi.star=try(get.Psi(TL,beta.star,tau.star,cells.idx,K,L),silent=TRUE)
            ll.star=try(dGenWish(Dobs,Psi.star,df,log=TRUE),silent=TRUE)
            if(is.numeric(ll.star)){
                if(!is.na(ll.star)){
                    pd=1
                }
            }
        }
        ## MH step
        mh1=ll.star+dmvnorm(beta.star,beta.prior.mean,beta.prior.cov,log=TRUE)+dnorm(tau.star,0,sd=sqrt(tau.prior.var),log=TRUE)+log(tau.star) ## log(tau.star) is from jacobian for change of variables
        mh2=ll+dmvnorm(beta,beta.prior.mean,beta.prior.cov,log=TRUE)+dnorm(tau,0,sd=sqrt(tau.prior.var),log=TRUE)+log(tau) ## log(tau) is from jacobian for change of variables

        if(runif(1)<exp(as.numeric(mh1-mh2))){
            theta=theta.star
            tau=tau.star
            beta=beta.star
            ll=ll.star
            accept=accept+1
        }

        ## DIC update: (at end of chain, dic.sum=Dbar)
        ##   done by taking second half of chain
        if(iter>(n.mcmc/2)){
            dic.sum=dic.sum+1/(n.mcmc/2)*(-2*ll)
        }

        ## Save values---------------------------------------------

        beta.save[iter,] <- beta
        tau.save[iter] <- tau
        ll.save[iter] <- ll
        logprior.save[iter] <- dmvnorm(beta,beta.prior.mean,beta.prior.cov,log=TRUE)+dnorm(tau,0,sd=sqrt(tau.prior.var),log=TRUE)+log(tau)

        ## Adaptive MCMC for beta ---------------------------------
        if (iter < adapt.max + 1 & iter / adapt.int == round(iter/adapt.int)){
            #########################################
            #### Roberts and Rosenthal adaptive here
            #########################################
            ## theta.tune <- 2.4^2 /length(theta) *
            ##     var(cbind(log(tau.save[1:iter]),beta.save[1:iter,]))+
            ##     diag(10^-10,length(theta))
            #########################################
            #### Shaby Wells Log Adaptive below here
            #########################################
            gamma.1=1/(iter/adapt.int)^.8
            gamma.2=1*gamma.1
            Sig.hat.current=var(cbind(log(tau.save[(iter-adapt.int+1):iter]),beta.save[(iter-adapt.int+1):iter,]))
            Sig.hat=Sig.hat+gamma.1*(Sig.hat.current-Sig.hat)
            s2.tune=exp(log(s2.tune+gamma.2*(accept/adapt.int-.234)))
            theta.tune=s2.tune*Sig.hat
            accept=0
        }
    }

    ## DIC --------------------------------------------------------
    ##browser()
    Dbar=dic.sum
    if(ncol(beta.save)==1) {
        beta.hat=mean(beta.save[-c(1:n.mcmc/2),])
    } else {
        beta.hat=apply(beta.save[-c(1:n.mcmc/2),],2,mean)
    }

    tau.hat=mean(tau.save[-c(1:n.mcmc/2)])

    Psi=get.Psi(TL,beta.hat,tau.hat,cells.idx,K,L)
    ll.hat=dGenWish(Dobs,Psi,df,log=TRUE)
    
    Dhat=-2*ll.hat

    DIC=2*Dbar-Dhat
    DIC2=Dhat + var(-2*ll.save[-c(1:n.mcmc/2)])

    list(beta=beta.save,tau=tau.save,n.mcmc=n.mcmc,ll=ll.save,logprior=logprior.save,
         accept=accept,DIC=DIC,DIC2=DIC2,Dbar=Dbar,Dhat=Dhat)
}
