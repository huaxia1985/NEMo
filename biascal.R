biascal <- function (dendistr,PO,BK,bias,avec,bvec,N,uvec,sd,xvec,edge,area) 
{
	rescal <- function (i,edge,PA,PO,BK,dendistr,avec,bvec,uvec,xvec,bias,area,N,sd) {
		j <-which(edge[,2]==i)
		a <- avec[j]
		b <- bvec[j]
		tmp <- sapply((seq(1:N)/(exp(a)+0.1))^(exp(b)+0.1),gammainc,a=1/(exp(b)+0.1))
		niche <- (tmp-c(0,tmp[-length(tmp)]))/2/gamma(1/(exp(b)+0.1))
		niche <- niche/niche[1]
		back <- pnorm(seq(2,N-1,1),mean=uvec[j],sd=sd)-pnorm(seq(1,N-2,1),mean=uvec[j],sd=sd)
		back <- c(pnorm(1,mean=uvec[j],sd=sd),back)
		back <- c(back,1-pnorm(N-1,mean=uvec[j],sd=sd))
		back[0:(xvec[j]-1)] <- 0
		back <- back*niche
		likPO <- sum(log(dendistr[j,PO[[i]]$x]))+sum(colSums(bias*t(PO[[i]]$bias)))-sum(exp(colSums(bias*t(BK[[i]]$bias)))*back[BK[[i]]$x])*area[i]
		likPO
	}
	likPO <- sapply(1:nspecies,rescal,edge,PA,PO,BK,dendistr,avec,bvec,uvec,xvec,bias,area,N,sd)
	likPO[is.na(likPO)] <- -10^7
	likPO[is.infinite(likPO)] <- -10^7
	likPO
}