addspeciation <- 
function (n,xmin,time,node,is.des,expinv.a,expinv.b,N,edge,edge.length,events,rootdistr,rootdistro,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
{
	type <- 1
	a <- sample(expinv.a,1,replace=T)
	a2 <- sample(expinv.a,1,replace=T)
	b <- sample(expinv.b,1,replace=T)
	b2 <- sample(expinv.b,1,replace=T)
	anc <- edge[node,1]
	if (anc==edge[1,1] || is.null(is.des)) {
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
	if (time==0) {
		if (xmin>0) {
			node2 <- setdiff(which(edge[,1]==edge[node,1]),node)
			tmp <- events
			tmp[nrow(tmp)+(1:2),] <-rbind(c(node,type,time,a,b,n,xmin),c(node2,type,time,a2,b2,-n,xmin))
		} else {
			tmp <- events
			tmp[nrow(tmp)+1,] <- c(node,type,time,a,b,n,-xmin)
		}
		tmp <- tmp[with(tmp,order(node,time)),]
		lik.new <- likcal(node=anc,root.included=F,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,d=d,bias=bias,area=area)
	} else {
		tmp <- events
		tmp[nrow(tmp)+1,] <- c(node,type,time,a,b,n,abs(xmin))
		tmp <- tmp[with(tmp,order(node,time)),]
		if (!is.null(is.des)) {
			lik.new <- likcal(node=ifelse(is.des,anc,node),root.included=ifelse(anc==0,F,T),events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,d=d,bias=bias,area=area)
		} else {
			lik.new <- likcal(node=0,root.included=F,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,d=d,bias=bias,area=area)
		}
	}
	out <- list(lik.new=lik.new,tmp=tmp)
	out
}
