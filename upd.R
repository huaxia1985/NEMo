upd <- function (events,edget,dendistrt,dendistrto,a,b,u,x,sd,N,d) {
  		if (nrow(events)!=0 && sum(dendistrt)>0) {
 			if (events$time[1]>0) {
				tmp <- sapply((seq(1:N)/(exp(a)+0.1))^(exp(b)+0.1),gammainc,a=1/(exp(b)+0.1))
				if ((!is.list(tmp)&&is.na(sum(tmp)))||is.list(tmp)) {dendistrt <- numeric(N)}
				nichedistrt <- (tmp-c(0,tmp[-length(tmp)]))/2/gamma(1/(exp(b)+0.1))
				nichedistrt <- nichedistrt/nichedistrt[1]
				backdistrt <- pnorm(seq(-(N-1),N,1),mean=u,sd=sd)-pnorm(seq(-N,N-1,1),mean=u,sd=sd)
				x <- events$xmin[1]
				if (x>min(which(dendistrt>0))) {dendistrt <- numeric(N)}
				if (x>1) {backdistrt[1:(x-1+N)] <- 0}
				dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*events$time[1]
				if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {dendistrt <- numeric(N)}
				dendistrto <- dendistrto/sum(dendistrto)
				dendistrt <- dendistrto[(N+1):(2*N)]
				dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
  			}
			for (z in 1:nrow(events)) {
				dendistrt2 <- dendistrt
				if (events$type[z]==1 && sum(dendistrto)>0) {
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
					if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {
						dendistrt <- numeric(N)
						break
					}
					dendistrto <- dendistrto/sum(dendistrto)
					dendistrt <- dendistrto[(N+1):(2*N)]
					dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
					if (sum(dendistrt==dendistrt2)==N) {
						dendistrt <- numeric(N)
						break
					}
					a <- a+events$a[z]*(-sum(dendistrt2*c(1:N))+sum(dendistrt*c(1:N)))
					b <- b+events$b[z]*(sum(dendistrt2*c(1:N))-sum(dendistrt*c(1:N)))
					tmp <- sapply((seq(1:N)/(exp(a)+0.1))^(exp(b)+0.1),gammainc,a=1/(exp(b)+0.1))
					if (((!is.list(tmp)&&is.na(sum(tmp)))||is.list(tmp))) {
						dendistrt <- numeric(N)
						break
					}
					nichedistrt <- (tmp-c(0,tmp[-length(tmp)]))/2/gamma(1/(exp(b)+0.1))
					nichedistrt <- nichedistrt/nichedistrt[1]
					backdistrt <- pnorm(seq(-(N-1),N,1),mean=u,sd=sd)-pnorm(seq(-N,N-1,1),mean=u,sd=sd)
					x <- events$xmin[z]
					if (x>min(which(dendistrt>0))) {
						dendistrt <- numeric(N)
						break
					}
					if (x>1) {backdistrt[1:(x-1+N)] <- 0}
					dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
					if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {
						dendistrt <- numeric(N)
						break
					}
					dendistrto <- dendistrto/sum(dendistrto)
					dendistrt <- dendistrto[(N+1):(2*N)]
					dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
					if (sum(dendistrt==dendistrt2)==N) {
						dendistrt <- numeric(N)
						break
					}
				}
				if (events$type[z]==2 && sum(dendistrt)>0) {
					tmp <- sapply((seq(1:N)/(exp(a)+0.1))^(exp(b)+0.1),gammainc,a=1/(exp(b)+0.1))
					if (((!is.list(tmp)&&is.na(sum(tmp)))||is.list(tmp))) {
					dendistrt <- numeric(N)
					break
				}
				nichedistrt <- (tmp-c(0,tmp[-length(tmp)]))/2/gamma(1/(exp(b)+0.1))
				nichedistrt <- nichedistrt/nichedistrt[1]
				na <- events$n[z]
				tmp <- numeric(2*N)
				tmp[(na+1):(2*N)] <- dendistrto[1:(2*N-na)]
				tmp[2*N] <- tmp[2*N]+sum(dendistrto[(2*N-na+1):(2*N)])
				dendistrto <- tmp*c(rep(1,N),nichedistrt)
				if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {
					dendistrt <- numeric(N)
					break
				}
				dendistrto <- dendistrto/sum(dendistrto)
				dendistrt <- dendistrto[(N+1):(2*N)]
				dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
				a <- a+events$a[z]*(-sum(dendistrt2*c(1:N))+sum(dendistrt*c(1:N)))
				b <- b+events$b[z]*(sum(dendistrt2*c(1:N))-sum(dendistrt*c(1:N)))
				u <- u+na
				tmp <- sapply((seq(1:N)/(exp(a)+0.1))^(exp(b)+0.1),gammainc,a=1/(exp(b)+0.1))
				if (((!is.list(tmp)&&is.na(sum(tmp)))||is.list(tmp))) {
					dendistrt <- numeric(N)
					break
				}
				nichedistrt <- (tmp-c(0,tmp[-length(tmp)]))/2/gamma(1/(exp(b)+0.1))
				nichedistrt <- nichedistrt/nichedistrt[1]
				backdistrt <- pnorm(seq(-(N-1),N,1),mean=u,sd=sd)-pnorm(seq(-N,N-1,1),mean=u,sd=sd)
				x <- events$xmin[z]
				if (x>min(which(dendistrt>0))) {
					dendistrt <- numeric(N)
					break
				}
				if (x>1) {backdistrt[1:(x-1+N)] <- 0}
				dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
				if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {
					dendistrt <- numeric(N)
					break
				}
				dendistrto <- dendistrto/sum(dendistrto)
				dendistrt <- dendistrto[(N+1):(2*N)]
				dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
			}
		}
	}
	result <- list(dendistrt=dendistrt,dendistrto=dendistrto,a=a,b=b,u=u,x=x)
	result
	}
