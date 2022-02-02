upd <- function (events,edget,dendistrt,dendistrto,nichedistrt,backdistrt,a,b,u,x,sd,N,d) {
  if (nrow(events)!=0 && sum(dendistrt)>0) {
  if (events$time[1]>0) {
	x <- events$xmin[1]
	if (x>min(which(dendistrt>0))) {
		dendistrt <- numeric(N)
    } else {
	if (x>1) {backdistrt[1:(x-1+N)] <- 0}
	if (sum(backdistrt)==0) {
		dendistrt <- numeric(N)
    } else {
    backdistrt <- backdistrt/sum(backdistrt)
	dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*events$time[1]
	if (sum(dendistrto)==0) {
		dendistrt <- numeric(N)
	}
    dendistrto <- dendistrto/sum(dendistrto)
	dendistrt <- dendistrto[(N+1):(2*N)]
	dendistrt[1] <- dendistrt[1]+sum(dendistrto[1:N])
    }
    }
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
    tmpo <- dendistrto*c(rep(1,N),nichedistrt)
    tmpt <- tmpo[(N+1):(2*N)]
	tmpt[1] <- tmpt[1]+sum(tmpo[1:N])
    a <- a+events$a[z]*(sum(dendistrt*c(1:N))-sum(tmpt*c(1:N)))
    b <- b+events$b[z]*(sum(tmpt*c(1:N))-sum(dendistrt*c(1:N)))
    nichedistrt <- exp(-((0:(N-1))/(exp(a)+1e-6))^(exp(b)+1e-6))
    x <- events$xmin[z]
    if (x>min(which(dendistrt>0))) {
      dendistrt <- numeric(N)
      break
    }
    if (x>1) {backdistrt[1:(x-1+N)] <- 0}
    if (sum(backdistrt)==0) {
      dendistrt <- numeric(N)
      break
    }
    backdistrt <- backdistrt/sum(backdistrt)
	dendistrto <- dendistrto+c(rep(1,N),nichedistrt)*backdistrt*d*(ifelse((z+1)<=nrow(events),events$time[z+1],edget)-events$time[z])
	if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {
		dendistrt <- numeric(N)
		break
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
	if (sum(dendistrto)==0 || sum(dendistrto<0)>0) {
		dendistrt <- numeric(N)
		break
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
		dendistrt <- numeric(N)
		break
	}
	if (x>1) {backdistrt[1:(x-1+N)] <- 0}
	if (sum(backdistrt)==0) {
      dendistrt <- numeric(N)
      break
    }
	backdistrt <- backdistrt/sum(backdistrt)
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
  result <- list(dendistrt=dendistrt,dendistrto=dendistrto,nichedistrt=nichedistrt,backdistrt=backdistrt,a=a,b=b,u=u,x=x)
  result
}
