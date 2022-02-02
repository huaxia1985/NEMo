simcode <- function (sp,mu,nspecies,A1,A2,rootu,sd,roota,rootb,rootx,d,N,r,pois.x,area,nBK) {
	eff_seed <- sample(1:2^15, 1)
	set.seed(eff_seed)
	samplet <- function(t1,t2) {
		integrand <- function (tt) {mu*(1-exp(-(sp-mu)*tt))/(sp-mu*exp(-(sp-mu)*tt))}
		tt <- runif(1)*(t1-t2)
		while (runif(1)>integrand(t1-tt)) {
			tt <- runif(1)*(t1-t2)
		}
		tt
	}
	rootniche <- exp(-((0:(N-1))/(exp(roota)+1e-6))^(exp(rootb)+1e-6))
	rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
    rootback <- rootback/sum(rootback)
	rootdistro <- c(rep(1,N),rootniche)*rootback
	rootdistro <- rootdistro/sum(rootdistro)
	rootdistr <- rootdistro[(N+1):(2*N)]
	rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
	library(diversitree)
	library(phytools)
	tree <- tree.bd(c(sp,mu),max.taxa=nspecies,include.extinct=T)
	tree <- try(prune(tree),silent=T)
	while(inherits(tree,"try-error")||is.null(tree)) {
		tree <- tree.bd(c(sp,mu),max.taxa=nspecies,include.extinct=T)
		tree <- try(prune(tree),silent=T)
	}
	t <- nodeHeights(tree)
	t <- max(t)-t
    integrand <- function (tt) {mu*(1-exp(-(sp-mu)*tt))/(sp-mu*exp(-(sp-mu)*tt))}
	edge.prob <- sapply(1:length(tree$edge.length),function (i) integrate(integrand,lower=t[i,2],upper=t[i,1])$value)
	tT1 <- sum(edge.prob)
	edge.prob1 <- edge.prob/tT1
	tT2 <- sum(tree$edge.length)
	edge.prob2 <- tree$edge.length/tT2
	edge <- tree$edge
	edge.length <- tree$edge.length
	k1 <- rpois(1,A1)
	k2 <- rpois(1,A2)
	events <- matrix(NA,k1+k2,7)
	events <- data.frame(events)
	names(events) <- c("node","type","time","a","b","n","xmin")
	events$type <- c(rep(1,k1),rep(2,k2))
	tmp <- which(runif(k1)<((nspecies-1)/(nspecies-1+sp*tT1)))
	events$node[tmp] <- sample(unique(edge[,1]),size=length(tmp),replace=F)
	sis <- sapply(events$node[tmp], function (i) which(edge[,1]==i)[2])
	events$node[tmp] <- sapply(events$node[tmp], function (i) which(edge[,1]==i)[1])
	events$time[tmp] <- 0
	events$node[c(1:k1)[-tmp]] <- sample.int(n=dim(tree$edge)[1],size=k1-length(tmp),replace=T,prob=edge.prob1)
	events$node[(k1+1):(k1+k2)] <- sample.int(n=dim(tree$edge)[1],size=k2,replace=T,prob=edge.prob2)
	events <- events[with(events,order(type,time,node)),]
	if (k1>length(tmp)) {
		events$time[(length(tmp)+1):k1] <- unlist(sapply(events$node[(length(tmp)+1):k1],function (z) samplet(t1=t[z,1],t2=t[z,2])))
	}
	events$time[(k1+1):(k1+k2)] <- unlist(sapply(unique(events$node[(k1+1):(k1+k2)]),function (z) runif(n=sum(events$node[(k1+1):(k1+k2)]==z),min=0,max=t[z,1]-t[z,2])))
	events$n[1:k1] <- sample(seq(-9,9,2)/10,k1,prob=sapply(seq(-9,9,2)/10,prior.n),replace=T)
    events$n[(k1+1):(k1+k2)] <- sample(c(-(N-1):(N-1)),size=k2,replace=T)
	events[(k1+k2+1):(k1+k2+length(tmp)),] <- cbind(sis,1,0,NA,NA,-events$n[1:length(tmp)],NA)
	events$a <- r
	events$b <- r
	events$xmin <- rpois(nrow(events),lambda=pois.x)+1
	events <- events[with(events,order(node,time)),]
	PAloc <- sample.int(10,1000,replace=T)
	sim <- try(simdistr(events=events,edge=tree$edge,edge.length=tree$edge.length,nspecies=nspecies,N=N,rootu=rootu,sd=sd,roota=roota,rootb=rootb,rootx=rootx,rootdistr=rootdistr,rootdistro=rootdistro,rootniche=rootniche,rootback=rootback,d=d,area=area,PAloc=PAloc,nBK=nBK),silent=T)
	if (!inherits(sim,"try-error")) {
		PA <- vector("list",nspecies)
		PO <- vector("list",nspecies)
		BK <- vector("list",nspecies)
		for (i in 1:length(PA)) {
			PA[[i]]$y <- sim$PA[[i]]
			PA[[i]]$x <- PAloc
			PO[[i]]$x <- sim$PO[[i]]
			BK[[i]]$x <- sim$BK[[i]]
		}
		out <- list(sim=sim,PA=PA,PO=PO,BK=BK,tree=tree,tT1=tT1,tT2=tT2,edge.prob1=edge.prob1,edge.prob2=edge.prob2,t=t,eff_seed=eff_seed)
	} else {
		out <- sim
	}
}
