simdistr <- function (events,edge,edge.length,nspecies,N,rootu,sd,roota,rootb,rootx,rootdistr,rootdistro,rootniche, rootback,d,area,PAloc,nBK) {
upd <- function (events,edget,dendistrt,dendistrto,nichedistrt,backdistrt,a,b,u,x,sd,N,d) {
  if (nrow(events)!=0 && sum(dendistrt)>0) {
  if (events$time[1]>0) {
	x <- events$xmin[1]
	if (x>min(which(dendistrt>0))) {
		events$xmin[1] <- x <- min(which(dendistrt>0))
	}
	if (x>1) {backdistrt[1:(x-1+N)] <- 0}
	if (sum(backdistrt)==0) {
		stop(paste(paste(events[z,],collapse=","),"no distribution"))
	}
    backdistrt <- backdistrt/sum(backdistrt)
	dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*events$time[1]
	if (sum(dendistrto)==0) {
		stop(paste(paste(events[z,],collapse=","),"no distribution"))
	}
	dendistrto <- dendistrto/sum(dendistrto)
	dendistrt <- dendistrto[(N+1):(2*N)]
	dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
  }
  for (z in 1:nrow(events)) {
  	dendistrt2 <- dendistrt
  	if (events$type[z]==1) {
  		pt <- events$n[z]
  		tmp1 <- cumsum(dendistrto)
  		tmp <- min(which(tmp1>=abs(pt)))
  		if (pt<0) {
  			if (tmp<(2*N)) {dendistrto[(tmp+1):(2*N)] <- 0}
  			dendistrto[tmp] <- dendistrto[tmp]-tmp1[tmp]+abs(pt)
  		} else {
  			if (tmp>1) {dendistrto[1:(tmp-1)] <- 0}
  			dendistrto[tmp] <- tmp1[tmp]-abs(pt)
    	}
    	if (sum(dendistrto)==0) {
    		stop(paste(paste(events[z,],collapse=","),"no distribution"))
    	}
    	dendistrto <- dendistrto/sum(dendistrto)
    	dendistrt <- dendistrto[(N+1):(2*N)]
		dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
    	if (sum(dendistrt==dendistrt2)==N) {
    		stop(paste(paste(events[z,],collapse=","),"speciation leads to no change in distribution"))
    	}
    	tmpo <- dendistrto*c(rep(1,N),nichedistrt)
    	tmpt <- tmpo[(N+1):(2*N)]
		tmpt[1] <- tmpt[1]+sum(tmpo[1:N])
    	a <- a+events$a[z]*(sum(dendistrt*c(1:N))-sum(tmpt*c(1:N)))
    	b <- b+events$b[z]*(sum(tmpt*c(1:N))-sum(dendistrt*c(1:N)))
    	nichedistrt <- exp(-((0:(N-1))/(exp(a)+1e-6))^(exp(b)+1e-6))
		x <- events$xmin[z]
		if (x>min(which(dendistrt>0))) {
			events$xmin[z] <- x <- min(which(dendistrt>0))
		}
		if (x>1) {backdistrt[1:(x-1+N)] <- 0}
		if (sum(backdistrt)==0) {
			stop(paste(paste(events[z,],collapse=","),"no distribution"))
		}
        backdistrt <- backdistrt/sum(backdistrt)
		dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
		if (sum(dendistrto)==0) {
			stop(paste(paste(events[z,],collapse=","),"no distribution"))
		}
		dendistrto <- dendistrto/sum(dendistrto)
		dendistrt <- dendistrto[(N+1):(2*N)]
		dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
  	}
  	if (events$type[z]==2) {
		na <- events$n[z]
		tmp <- numeric(2*N)
        if (na>0) {
		tmp[(na+1):(2*N)] <- dendistrto[1:(2*N-na)]
		tmp[2*N] <- tmp[2*N]+sum(dendistrto[(2*N-na+1):(2*N)])
		dendistrto <- tmp*c(rep(1,N),nichedistrt)
        }
        if (na<0){
        tmp[1:(2*N+na)] <- dendistrto[(-na+1):(2*N)]
        tmp[1] <- tmp[1]+sum(dendistrto[1:(-na)])
        dendistrto <- tmp
        }
		if (sum(dendistrto)==0) {
			stop(paste(paste(events[z,],collapse=","),"no distribution"))
		}
		dendistrto <- dendistrto/sum(dendistrto)
		dendistrt <- dendistrto[(N+1):(2*N)]
		dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
		u <- u+na
        tmp <- sum(dendistrt2*c(1:N))+na-sum(dendistrt*c(1:N))
		if (na>0 && tmp>0) {
        a <- a+events$a[z]*tmp
        b <- b-events$b[z]*tmp
        nichedistrt <- exp(-((0:(N-1))/(exp(a)+1e-6))^(exp(b)+1e-6))
        }
		backdistrt <- pnorm(seq(-(N-1),N,1),mean=u,sd=sd)-pnorm(seq(-N,N-1,1),mean=u,sd=sd)
		x <- events$xmin[z]
		if (x>min(which(dendistrt>0))) {
			events$xmin[z] <- x <- min(which(dendistrt>0))
		}
		if (x>1) {backdistrt[1:(x-1+N)] <- 0}
		if (sum(backdistrt)==0) {
			stop(paste(paste(events[z,],collapse=","),"no distribution"))
		}
        backdistrt <- backdistrt/sum(backdistrt)
		dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
		if (sum(dendistrto)==0) {
			stop(paste(paste(events[z,],collapse=","),"no distribution"))
		}
		dendistrto <- dendistrto/sum(dendistrto)
		dendistrt <- dendistrto[(N+1):(2*N)]
		dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
	}
  }
  }
  result <- list(events=events,dendistrt=dendistrt,dendistrto=dendistrto,nichedistrt=nichedistrt,backdistrt=backdistrt,a=a,b=b,u=u,x=x)
  result
}
##evolve niches along phylogeny
	dendistr <- matrix(NA,dim(edge)[1],N)
	rownames(dendistr) <- edge[,2]
	dendistro <- matrix(NA,dim(edge)[1],2*N)
	rownames(dendistro) <- edge[,2]
    nichedistr <- matrix(NA,dim(edge)[1],N)
    rownames(nichedistr) <- edge[,2]
    backdistr <- matrix(NA,dim(edge)[1],2*N)
    rownames(backdistr) <- edge[,2]
	avec <- numeric(dim(edge)[1])
	names(avec) <- edge[,2]
	bvec <- numeric(dim(edge)[1])
	names(bvec) <- edge[,2]
	uvec <- numeric(dim(edge)[1])
	names(uvec) <- edge[,2]
	xvec <- numeric(dim(edge)[1])
	names(xvec) <- edge[,2]
	des <- which(edge[,1]==edge[1,1])
    j <- des[1]
    res <- upd(events[events$node==j,],edget=edge.length[j],dendistrt=rootdistr,dendistrto=rootdistro,nichedistrt=rootniche,backdistr=rootback,a=roota,b=rootb,u=rootu,x=rootx,sd,N,d)
    dendistr[j,] <- res$dendistrt
    dendistro[j,] <- res$dendistrto
    nichedistr[j,] <- res$nichedistrt
    backdistr[j,] <- res$backdistrt
    avec[j] <- res$a
    bvec[j] <- res$b
    uvec[j] <- res$u
    xvec[j] <- res$x
    events[events$node==j,] <- res$events
    j<-j+1
    while(j<=dim(edge)[1] && edge[j,1]>edge[des[1],1]) {
    		tmp <- as.character(edge[j,1])
    		res <- upd(events[events$node==j,],edget=edge.length[j],dendistrt=dendistr[tmp,],dendistrto=dendistro[tmp,],nichedistrt=nichedistr[tmp,],backdistrt=backdistr[tmp,],a=avec[tmp],b=bvec[tmp],u=uvec[tmp],x=xvec[tmp],sd,N,d)
    		dendistr[j,] <- res$dendistrt
    		dendistro[j,] <- res$dendistrto
            nichedistr[j,] <- res$nichedistrt
            backdistr[j,] <- res$backdistrt
    		avec[j] <- res$a
    		bvec[j] <- res$b
    		uvec[j] <- res$u
    		xvec[j] <- res$x
    		events[events$node==j,] <- res$events
    		j<-j+1
    		if (j>dim(edge)[1]) {break}
    	}
    	j <- des[2]
    	res <- upd(events[events$node==j,],edget=edge.length[j],dendistrt=rootdistr,dendistrto=rootdistro,nichedistrt=rootniche,backdistr=rootback,a=roota,b=rootb,u=rootu,x=rootx,sd,N,d)
    	dendistr[j,] <- res$dendistrt
    	dendistro[j,] <- res$dendistrto
        nichedistr[j,] <- res$nichedistrt
        backdistr[j,] <- res$backdistrt
    	avec[j] <- res$a
    	bvec[j] <- res$b
    	uvec[j] <- res$u
    	xvec[j] <- res$x
    	events[events$node==j,] <- res$events
    	j<-j+1
    	while(j<=dim(edge)[1] && edge[j,1]>edge[des[2],1]) {
    		tmp <- as.character(edge[j,1])
    		res <- upd(events[events$node==j,],edget=edge.length[j],dendistrt=dendistr[tmp,],dendistrto=dendistro[tmp,],nichedistrt=nichedistr[tmp,],backdistrt=backdistr[tmp,],a=avec[tmp],b=bvec[tmp],u=uvec[tmp],x=xvec[tmp],sd,N,d)
    		dendistr[j,] <- res$dendistrt
    		dendistro[j,] <- res$dendistrto
            nichedistr[j,] <- res$nichedistrt
            backdistr[j,] <- res$backdistrt
    		avec[j] <- res$a
    		bvec[j] <- res$b
    		uvec[j] <- res$u
    		xvec[j] <- res$x
    		events[events$node==j,] <- res$events
    		j<-j+1
    		if (j>dim(edge)[1]) {break}
       }
##extract tip species distribution data
backdistro <- backdistr[,(N+1):(2*N)]
backdistro[,1] <- rowSums(backdistr[,1:(N+1)])
    BK <- lapply(1:nspecies,function (i) sample.int(n=N,size=nBK,replace=T,prob=backdistro[as.character(i),]))
    rateall <- area*rowSums(nichedistr[as.character(1:nspecies),]*backdistro[as.character(1:nspecies),])
    nloc <- sapply(rateall,rpois,n=1)+1
    PO <- lapply(1:nspecies,function (i) sample.int(n=N,size=nloc[i],replace=T,prob=dendistr[as.character(i),]))
    PA <- lapply(1:nspecies,function (i) as.numeric(runif(length(PAloc))>=exp(-dendistr[as.character(i),PAloc])))
    out <- list(PO=PO,PA=PA,BK=BK,events=events,dendistr=dendistr,backdistr=backdistr,nichedistr=nichedistr,avec=avec,bvec=bvec,uvec=uvec,xvec=xvec)
    out
}
