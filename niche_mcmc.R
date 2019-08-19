niche_mcmc <- function (ngen, ...)
{
sample.vec <- function(x, ...) x[sample(length(x), ...)]

#### start MCMC
for (tt in 1:ngen) {
# update the location of an event
if (nrow(events)!=0) {
	tmp <- which(events$time==0)
	if (length(tmp)>0) {
	i <- 1
	tmp2 <- which(events$node[tmp]==setdiff(which(edge[,1]==edge[events$node[tmp[i]],1]),events$node[tmp[i]]))
	tmp <- tmp[-tmp2]
	i <- i+1
	while(i<=length(tmp)) {
		tmp2 <- which(events$node[tmp]==setdiff(which(edge[,1]==edge[events$node[tmp[i]],1]),events$node[tmp[i]]))
		tmp <- tmp[-tmp2]
		i <- i+1 
	}
	}
	tmp2 <- which(events$time!=0)
	idx <- sample.vec(c(tmp2,tmp),1)
	node <- events[idx,]$node
	type <- events[idx,]$type
	n <- events[idx,]$n
	xmin <- events[idx,]$xmin
	time <- events[idx,]$time
	if (type==1 && time==0) {
		node2 <- setdiff(which(edge[,1]==edge[node,1]),node)
		idx2 <-  which((events$type==1)*(events$node==node2)*(events$time==0)==1)
		tmp <- events[-c(idx,idx2),]
		node.new <- sample.int(dim(edge)[1],size=1)
		idx2 <- which((tmp$node==node.new)*(tmp$type==1)*(tmp$time==0)==1)
		if (length(idx2)>0) {
			tmp$node[idx2] <- node
			node.new2 <- setdiff(which(edge[,1]==edge[node.new,1]),node.new)
			idx2 <-  which((tmp$node==node.new2)*(tmp$type==1)*(tmp$time==0)==1)
			tmp$node[idx2] <- node2
		}
		out <- addspeciation(n,xmin,time=0,node.new,is.des=NULL,expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
		lik.new <- out$lik.new
		accept <- min(1,exp(sum(lik.old$res)-sum(lik.new$res)))
		if (runif(1)<accept) {
			lik.old <- lik.new
			events <- out$tmp
		}
	} else {
		if (runif(1)<0.1) {
			anc <- which(edge[,2]==edge[node,1])
			node.new <- sample.vec(c(anc,which(edge[,1]==edge[node,2])),1)
			if (length(node.new)>0) {
			time.new <- runif(n=1,min=0,max=t[node,1]-t[node,2])
			tmp <- events[-idx,]
			if (type==1) {
				out <- addspeciation(n,xmin,time.new,node.new,is.des=ifelse((length(anc)>0 && node.new==anc),F,T),expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
				r <- integrand(t[node.new,1]-time.new)/integrand(t[node,1]-time)
			}
			if (type==2) {
				out <- addadaptation(n,xmin,time.new,node.new,is.des=ifelse((length(anc)>0 && node.new==anc),F,T),expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
				r <- 1
			}
			lik.new <- out$lik.new
			accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*r)
			if (runif(1)<accept) {
				lik.old$res[lik.new$idx] <- lik.new$res
				lik.old$likPO[lik.new$idx] <- lik.new$likPO
				ii <- which(!is.na(lik.new$dendistr[,1]))
				lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
				lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
				lik.old$avec[ii] <- lik.new$avec[ii]
				lik.old$bvec[ii] <- lik.new$bvec[ii]
				lik.old$uvec[ii] <- lik.new$uvec[ii]
				lik.old$xvec[ii] <- lik.new$xvec[ii]
				events <- out$tmp
			}
			}
		} else {
			u <- runif(n=1,min=-tuning.t,max=tuning.t)
			time.new <- time+u
			while (time.new==0 || time.new==(t[node,1]-t[node,2])) {
				u <- runif(n=1,min=-tuning.t,max=tuning.t)
				time.new <- time+u
			}
			is.des <- F
			if (time.new < 0) {
				tmp <- which(edge[,2]==edge[node,1])
				if (length(tmp)>0) {
					node.new <- tmp
					time.new <- t[node.new,1]-t[node.new,2]+time.new
				} else {
					node.new <- node
					time.new <- time-u
				}
			} else if (time.new > (t[node,1]-t[node,2])) {
				tmp <- which(edge[,1]==edge[node,2])
				if (length(tmp)>0) {
					node.new <- sample(tmp,size=1)
					time.new <- time.new-t[node,1]+t[node,2]
					is.des <- T
				} else {
					node.new <- node
					time.new <- time-u
				}
			} else {
				node.new <- node
			}
			tmp <- events[-idx,]
			if (type==1) {
				out <- addspeciation(n,xmin,time.new,node.new,is.des,expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
				r <- integrand(t[node.new,1]-time.new)/integrand(t[node,1]-time)
			}
			if (type==2) {
				out <- addadaptation(n,xmin,time.new,node.new,is.des,expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
				r <- 1
			}
			lik.new <- out$lik.new
			accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*r)
			if (runif(1)<accept) {
			lik.old$res[lik.new$idx] <- lik.new$res
			lik.old$likPO[lik.new$idx] <- lik.new$likPO
			ii <- which(!is.na(lik.new$dendistr[,1]))
			lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
			lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
			lik.old$avec[ii] <- lik.new$avec[ii]
			lik.old$bvec[ii] <- lik.new$bvec[ii]
			lik.old$uvec[ii] <- lik.new$uvec[ii]
			lik.old$xvec[ii] <- lik.new$xvec[ii]
			events <- out$tmp
			}
		}
	}
}
# updating the parameters of an event
if (nrow(events)!=0) {
	tmp <- which(events$time==0)
	if (length(tmp)>0) {
	i <- 1
	tmp2 <- which(events$node[tmp]==setdiff(which(edge[,1]==edge[events$node[tmp[i]],1]),events$node[tmp[i]]))
	tmp <- tmp[-tmp2]
	i <- i+1
	while(i<=length(tmp)) {
		tmp2 <- which(events$node[tmp]==setdiff(which(edge[,1]==edge[events$node[tmp[i]],1]),events$node[tmp[i]]))
		tmp <- tmp[-tmp2]
		i <- i+1 
	}
	}
	tmp2 <- which(events$time!=0)
	idx <- sample.vec(c(tmp2,tmp),1)
	node <- events[idx,]$node
	type <- events[idx,]$type
	n <- events[idx,]$n
	xmin <- events[idx,]$xmin
	time <- events[idx,]$time
	if (type==1 && time==0) {
			node2 <- setdiff(which(edge[,1]==edge[node,1]),node)
			idx2 <-  which((events$node==node2)*(events$time==0)==1)
			tmp <- events[-c(idx,idx2),]
		} else {
			tmp <- events[-idx,]
		}
	if (type==1) {
		n <- round(n,1)
		if (n==-0.9) {
			n.new <- n+0.2
			r <- 0.5*prior.n(n.new)/prior.n(n)
		} else if (n==0.9) {
			n.new <- n-0.2
			r <- 0.5*prior.n(n.new)/prior.n(n)
		} else {
			n.new <- n+ifelse(runif(1)<=0.5,0.2,-0.2)
			r <- 1*prior.n(n.new)/prior.n(n)
		}
		out <- addspeciation(n.new,xmin,time,node,is.des=F,expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
	}
	if (type==2) {
		if (n==1) {
			n.new <- n+1
			r <- 0.5
		} else if (n==(N-1)) {
			n.new <- n-1
			r <- 0.5
		} else {
			n.new <- n+ifelse(runif(1)<=0.5,1,-1)
			r <- 1
		}
		out <- addadaptation(n.new,xmin,time,node,is.des=F,expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
	}
	lik.new <- out$lik.new
	accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*r)
	if (runif(1)<accept) {
		lik.old$res[lik.new$idx] <- lik.new$res
		lik.old$likPO[lik.new$idx] <- lik.new$likPO
		ii <- which(!is.na(lik.new$dendistr[,1]))
		lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
		lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
		lik.old$avec[ii] <- lik.new$avec[ii]
		lik.old$bvec[ii] <- lik.new$bvec[ii]
		lik.old$uvec[ii] <- lik.new$uvec[ii]
		lik.old$xvec[ii] <- lik.new$xvec[ii]
		events <- out$tmp
	}
	node <- events[idx,]$node
	type <- events[idx,]$type
	n <- events[idx,]$n
	xmin <- events[idx,]$xmin
	time <- events[idx,]$time
	tmp <- events[-idx,]
	if (xmin==1) {
		xmin.new <- xmin+1
		r <- 0.5
	} else if (xmin==(N-1)) {
		xmin.new <- xmin-1
		r <- 0.5
	} else {
		xmin.new <- xmin+ifelse(runif(1)<=0.5,1,-1)
		r <- 1
	}
	if (type==1) {
		out <- addspeciation(n,-xmin.new,time,node,is.des=F,expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
	}
	if (type==2) {
		out <- addadaptation(n,xmin.new,time,node,is.des=F,expinv.a,expinv.b,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
	}
	lik.new <- out$lik.new
	accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*r)
	if (runif(1)<accept) {
		lik.old$res[lik.new$idx] <- lik.new$res
		lik.old$likPO[lik.new$idx] <- lik.new$likPO
		ii <- which(!is.na(lik.new$dendistr[,1]))
		lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
		lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
		lik.old$avec[ii] <- lik.new$avec[ii]
		lik.old$bvec[ii] <- lik.new$bvec[ii]
		lik.old$uvec[ii] <- lik.new$uvec[ii]
		lik.old$xvec[ii] <- lik.new$xvec[ii]
		events <- out$tmp
	}
	}
#adding a speciation event
	if (runif(1)<((nspecies-1)/(nspecies-1+sp*tT1))) {
		tmp <-  events$node[which((events$type==1)*(events$time==0)==1)]
		if (length(tmp)>0 && length(tmp)<dim(edge)[1]) {
			node <- sample.vec(c(1:dim(edge)[1])[-tmp],size=1)
			time <- 0
		} else if (length(tmp)==0) {
			node <- sample.vec(c(1:dim(edge)[1]),size=1)
			time <- 0
		} else {
			node <- sample.int(n=dim(edge)[1],size=1,prob=edge.prob1)
			time <- runif(n=1,min=0,max=t[node,1]-t[node,2])
			while(time==0) {time <- runif(n=1,min=0,max=t[node,1]-t[node,2])}
		}
	} else {
		node <- sample.int(n=dim(edge)[1],size=1,prob=edge.prob1)
		time <- runif(n=1,min=0,max=t[node,1]-t[node,2])
		while(time==0) {time <- runif(n=1,min=0,max=t[node,1]-t[node,2])}
	}
	n <- sample(seq(-9,9,2)/10,1,prob=sapply(seq(-9,9,2)/10,prior.n))
	out <- addspeciation(n=n,xmin=sample.int(N-1,1),time,node,is.des=F,expinv.a,expinv.b,N,edge,edge.length,events,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
	lik.new <- out$lik.new
	if (time==0) {
		accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*A[1]/(k[1]+1)*ifelse(k[1]==0,0.5,1)*prior.n(n)/(N-1)/exp.a/exp.b)
	} else {
		accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*A[1]/(k[1]+1)*ifelse(k[1]==0,0.5,1)*sp*integrand(t[node,1]-time)*prior.n(n)/(N-1)/exp.a/exp.b)
	}
	if (runif(1)<accept) {
		lik.old$res[lik.new$idx] <- lik.new$res
		lik.old$likPO[lik.new$idx] <- lik.new$likPO
		ii <- which(!is.na(lik.new$dendistr[,1]))
			lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
			lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
			lik.old$avec[ii] <- lik.new$avec[ii]
			lik.old$bvec[ii] <- lik.new$bvec[ii]
			lik.old$uvec[ii] <- lik.new$uvec[ii]
			lik.old$xvec[ii] <- lik.new$xvec[ii]
		events <- out$tmp
		k[1] <- k[1]+1
	}
#adding an adaptation event
    node <- sample.int(n=dim(edge)[1],size=1,prob=edge.prob2)
    time <- runif(n=1,min=0,max=t[node,1]-t[node,2])
    out <- addadaptation(n=sample.int(N-1,1),xmin=sample.int(N-1,1),time,node,is.des=F,expinv.a,expinv.b,N,edge,edge.length,events,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
    lik.new <- out$lik.new	
	accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*A[2]/(k[2]+1)*ifelse(k[2]==0,0.5,1)/tT2/(N-1)/(N-1)/exp.a/exp.b)
	if (runif(1)<accept) {
		lik.old$res[lik.new$idx] <- lik.new$res
		lik.old$likPO[lik.new$idx] <- lik.new$likPO
		ii <- which(!is.na(lik.new$dendistr[,1]))
			lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
			lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
			lik.old$avec[ii] <- lik.new$avec[ii]
			lik.old$bvec[ii] <- lik.new$bvec[ii]
			lik.old$uvec[ii] <- lik.new$uvec[ii]
			lik.old$xvec[ii] <- lik.new$xvec[ii]
		events <- out$tmp
		k[2] <- k[2]+1
	}
#delete a speciation event
	if (k[1]>0) {
		idx <- sample.vec(which(events$type==1),1)
		node <- events$node[idx]
		anc <- edge[node,1]
		if (anc==edge[1,1]) {
			nodedistr <- rootdistr
			nodedistro <- rootdistro
			nodeu <- rootu
			nodea <- roota
			nodeb <- rootb
			nodex <- rootx
			anc <- 0
		} else {
			anc <- as.character(anc)
			nodedistr <- lik.old$dendistr[anc,]
			nodedistro <- lik.old$dendistro[anc,]
			nodeu <- lik.old$uvec[anc]
			nodea <- lik.old$avec[anc]
			nodeb <- lik.old$bvec[anc]
			nodex <- lik.old$xvec[anc]
			anc <- which(edge[,2]==anc)
		}
		if (events$time[idx]==0) {
			node2 <- setdiff(which(edge[,1]==edge[node,1]),node)
			idx2 <-  which((events$node==node2)*(events$time==0)==1)
			tmp <- events[-c(idx,idx2),]
			lik.new <- likcal(node=anc,root.included=F,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,d=d,bias=bias,area=area)
			accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*k[1]/A[1]*ifelse(k[1]==1,2,1)/prior.n(n)*(N-1)*exp.a*exp.b)
		} else {
			tmp <- events[-idx,]
			lik.new <- likcal(node=node,root.included=T,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,d=d,bias=bias,area=area)
			accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*k[1]/A[1]*ifelse(k[1]==1,2,1)/sp/integrand(t[node,1]-time)/prior.n(n)*(N-1)*exp.a*exp.b)
		}
		if (runif(1)<accept) {
			lik.old$res[lik.new$idx] <- lik.new$res
			lik.old$likPO[lik.new$idx] <- lik.new$likPO
			ii <- which(!is.na(lik.new$dendistr[,1]))
			lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
			lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
			lik.old$avec[ii] <- lik.new$avec[ii]
			lik.old$bvec[ii] <- lik.new$bvec[ii]
			lik.old$uvec[ii] <- lik.new$uvec[ii]
			lik.old$xvec[ii] <- lik.new$xvec[ii]
			events <- tmp
			k[1] <- k[1]-1
		}
	}
#delete an adaptation event
	if (k[2]>0) {
		idx <- sample.vec(which(events$type==2),1)
		node <- events$node[idx]
		anc <- edge[node,1]
		if (anc==edge[1,1]) {
			nodedistr <- rootdistr
			nodedistro <- rootdistro
			nodeu <- rootu
			nodea <- roota
			nodeb <- rootb
			nodex <- rootx
		} else {
			anc <- as.character(anc)
			nodedistr <- lik.old$dendistr[anc,]
			nodedistro <- lik.old$dendistro[anc,]
			nodeu <- lik.old$uvec[anc]
			nodea <- lik.old$avec[anc]
			nodeb <- lik.old$bvec[anc]
			nodex <- lik.old$xvec[anc]
		}
		tmp <- events[-idx,]
		lik.new <- likcal(node=node,root.included=T,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,d=d,bias=bias,area=area)
		accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*k[2]/A[2]*ifelse(k[2]==1,2,1)*tT2*(N-1)*(N-1)*exp.a*exp.b)
		if (runif(1)<accept) {
			lik.old$res[lik.new$idx] <- lik.new$res
			lik.old$likPO[lik.new$idx] <- lik.new$likPO
			ii <- which(!is.na(lik.new$dendistr[,1]))
			lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
			lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
			lik.old$avec[ii] <- lik.new$avec[ii]
			lik.old$bvec[ii] <- lik.new$bvec[ii]
			lik.old$uvec[ii] <- lik.new$uvec[ii]
			lik.old$xvec[ii] <- lik.new$xvec[ii]
			events <- tmp
			k[2] <- k[2]-1
		}
	}
#updating events occurence rate
	u <- runif(1)
	A1 <- A[1]*(nspecies-1+sp*tT1)
	tmp <- A1*exp(tuning.A[1]*(u-0.5))
	accept <- min(1,(tmp/A1)^k[1]*exp(-tmp+A1)*dexp(tmp,rate=exp.A[1])/dexp(A1,rate=exp.A[1])*exp(tuning.A[1]*(u-0.5)))
	if (runif(1)<accept) {
		A[1] <- tmp/(nspecies-1+sp*tT1)
	}
	u <- runif(1)
	tmp <- A[2]*exp(tuning.A[2]*(u-0.5))
	accept <- min(1,(tmp/A[2])^k[2]*exp(-tmp+A[2])*dexp(tmp,rate=exp.A[2])/dexp(A[2],rate=exp.A[2])*exp(tuning.A[2]*(u-0.5)))
	if (runif(1)<accept) {
		A[2] <- tmp
	}
#updating ancestral niches
	if (tt<=1000) {
	ua <- runif(1)
	ub <- runif(1)
	usd <- runif(1)
    uu <- runif(1)
	tmpa <- roota*exp((ua-0.5)*tuning.roota)
	tmpb <- rootb*exp((ub-0.5)*tuning.rootb)
	tmpsd <- sd*exp((usd-0.5)*tuning.sd)
    tmpu <- rootu+uu*tuning.rootu
	rootniche <- sapply((seq(1:N)/(exp(tmpa)+0.1))^(exp(tmpb)+0.1),gammainc,a=1/(exp(tmpb)+0.1))
	rootniche <- (rootniche-c(0,rootniche[-length(rootniche)]))/2/gamma(1/(exp(tmpb)+0.1))
	rootniche <- rootniche/rootniche[1]
	rootback <- pnorm(seq(-(N-1),N,1),mean=tmpu,sd=tmpsd)-pnorm(seq(-N,N-1,1),mean=tmpu,sd=tmpsd)
	rootdistrotmp <- c(rep(1,N),rootniche)*rootback
	rootdistrotmp <- rootdistrotmp/sum(rootdistrotmp)
	rootdistrtmp <- rootdistrotmp[(N+1):(2*N)]
	rootdistrtmp[1] <- rootdistrtmp[1]+sum(rootdistrotmp[1:N])
	lik.new <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=tmpu,sd=tmpsd,nodea=tmpa,nodeb=tmpb,nodex=rootx,nodedistr=rootdistrtmp,nodedistro=rootdistrotmp,d=d,bias=bias,area=area)
	r <- exp(sum(lik.old$res)-sum(lik.new$res))*dlnorm(tmpa,meanlog=lnorm.roota.m,sdlog=lnorm.roota)/dlnorm(roota,meanlog=lnorm.roota.m,sdlog=lnorm.roota)*exp(tuning.roota*(ua-0.5))*dlnorm(tmpb,meanlog=lnorm.rootb.m,sdlog=lnorm.rootb)/dlnorm(rootb,meanlog=lnorm.rootb.m,sdlog=lnorm.rootb)*exp(tuning.rootb*(ub-0.5))*dlnorm(tmpsd,meanlog=lnorm.sd.m,sdlog=lnorm.sd)/dlnorm(sd,meanlog=lnorm.sd.m,sdlog=lnorm.sd)*exp(tuning.sd*(usd-0.5))*dnorm(tmpu,mean=norm.rootu.m,sd=norm.rootu)/dnorm(rootu,mean=norm.rootu.m,sd=norm.rootu)
	accept <- min(1,ifelse(is.na(r),0,r))
	if (runif(1)<accept) {
		lik.old <- lik.new
		roota <- tmpa
		rootb <- tmpb
		sd <- tmpsd
        rootu <- tmpu
		rootdistr <- rootdistrtmp
        rootdistro <- rootdistrotmp
	}
#updating dispersal rate
	u <- runif(1)
	tmp <- d*exp((u-0.5)*tuning.d)
	lik.new <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,d=tmp,bias=bias,area=area)
	#print(11)
	r <- exp(sum(lik.old$res)-sum(lik.new$res))*dexp(tmp,rate=exp.d)/dexp(d,rate=exp.d)*exp(tuning.d*(u-0.5))
	accept <- min(1,ifelse(is.na(r),0,r))
	if (runif(1)<accept) {
		lik.old <- lik.new
		d <- tmp
	}
#updating adaptation rate
	if (dim(events)[1]>0) {
	u <- runif(1)
	tmp <- expinv.a*exp((u-0.5)*tuning.a)
	eventstmp <- events
	eventstmp$a <- tmp
	eventstmp$b <- tmp
	lik.new <- likcal(node=0,root.included=F,events=eventstmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,d=tmp,bias=bias,area=area)
	r <- exp(sum(lik.old$res)-sum(lik.new$res))*dexp(tmp,rate=exp.a)/dexp(expinv.a,rate=exp.a)*exp(tuning.a*(u-0.5))
	accept <- min(1,ifelse(is.na(r),0,r))
	if (runif(1)<accept) {
		lik.old <- lik.new
		expinv.a <- tmp
		expinv.b <- tmp
		events <- eventstmp
	}
	}
}
#updating coefficients for sampling bias
if (!is.null(bias)) {
if (runif(1)<0.01) {
	u <- runif(1)
	tmp <- bias
	tmp[1] <- bias[1]*exp(tuning.bias1*(u-0.5))
	lik.new <- biascal(dendistr=lik.old$dendistr,PO=PO,BK=BK,bias=tmp,avec=lik.old$avec,bvec=lik.old$bvec,N=N,uvec=lik.old$uvec,sd=sd,xvec=lik.old$xvec,edge=edge,area=area)
	#print(13)
	r <- exp(sum(lik.new)-sum(lik.old$likPO))*dexp(-tmp[1],rate=exp.bias1)/dexp(-bias[1],rate=exp.bias1)*exp(tuning.bias1*(u-0.5))
	accept <- min(1,ifelse(is.na(r),0,r))
	if (runif(1)<accept) {
		lik.old$res <- lik.old$res+lik.old$likPO-lik.new
		lik.old$likPO <- lik.new
		bias <- tmp
	}
	u <- runif(1)
	tmp <- bias
	tmp[2] <- bias[2]*exp(tuning.bias2*(u-0.5))
	lik.new <- biascal(dendistr=lik.old$dendistr,PO=PO,BK=BK,bias=tmp,avec=lik.old$avec,bvec=lik.old$bvec,N=N,uvec=lik.old$uvec,sd=sd,xvec=lik.old$xvec,edge=edge,area=area)
	#print(14)
	r <- exp(sum(lik.new)-sum(lik.old$likPO))*dexp(-tmp[2],rate=exp.bias2)/dexp(-bias[2],rate=exp.bias2)*exp(tuning.bias2*(u-0.5))
	accept <- min(1,ifelse(is.na(r),0,r))
	if (runif(1)<accept) {
		lik.old$res <- lik.old$res+lik.old$likPO-lik.new
		lik.old$likPO <- lik.new
		bias <- tmp
	}
}
# write output
if (tt %% 1000 ==0) {
	prior <- dlnorm(roota,meanlog=lnorm.roota.m,sdlog=lnorm.roota,log=T)+dlnorm(rootb,meanlog=lnorm.rootb.m,sdlog=lnorm.rootb,log=T)+dlnorm(sd,meanlog=lnorm.sd.m,sdlog=lnorm.sd,log=T)+dexp(d,rate=exp.d,log=T)+dnorm(rootu,mean=norm.rootu.m,sd=norm.rootu,log=T)+dexp(expinv.a,rate=exp.a,log=T)+dexp(expinv.b,rate=exp.b,log=T)+dexp(A[1],rate=exp.A[1],log=T)+dexp(A[2],rate=exp.A[2],log=T)+dexp(-bias[1],rate=exp.bias1,log=T)+dexp(-bias[2],rate=exp.bias2,log=T)
	suppressWarnings(write.table(events,eventsfilename,append=T,col.names=T))
	write.table(t(as.numeric(c(tt,sum(lik.old$res),prior,roota,rootb,rootu,sd,d,expinv.a,expinv.b,A[1],A[2],bias))),postfilename,append=T,col.names=F)
	print(tt)
}
} else {
if (tt %% 1000 ==0) {
	prior <- dlnorm(roota,meanlog=lnorm.roota.m,sdlog=lnorm.roota,log=T)+dlnorm(rootb,meanlog=lnorm.rootb.m,sdlog=lnorm.rootb,log=T)+dlnorm(sd,meanlog=lnorm.sd.m,sdlog=lnorm.sd,log=T)+dexp(d,rate=exp.d,log=T)+dnorm(rootu,mean=norm.rootu.m,sd=norm.rootu,log=T)+dexp(expinv.a,rate=exp.a,log=T)+dexp(expinv.b,rate=exp.b,log=T)+dexp(A[1],rate=exp.A[1],log=T)+dexp(A[2],rate=exp.A[2],log=T)
	suppressWarnings(write.table(events,eventsfilename,append=T,col.names=T))
	write.table(t(as.numeric(c(tt,sum(lik.old$res),prior,roota,rootb,rootu,sd,d,expinv.a,expinv.b,A[1],A[2]))),postfilename,append=T,col.names=F)
	print(tt)
}
}
}
}