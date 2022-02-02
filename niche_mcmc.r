niche_mcmc <- function (ngen, ...)
{
sample.vec <- function(x, ...) x[sample(length(x), ...)]

# start MCMC
for (tt in 1:ngen) {
    if (nrow(events)!=0) {
    # update the location of the event
    tmp <- which(events$time==0)
    if (length(tmp)>0) {
    i <- 1
    tmp2 <- which(events$node[tmp]==setdiff(which(edge[,1]==edge[events$node[tmp[i]],1]),events$node[tmp[i]]))
    if (length(tmp2)>0) {tmp <- tmp[-tmp2]}
    i <- i+1
    while(i<=length(tmp)) {
        tmp2 <- which(events$node[tmp]==setdiff(which(edge[,1]==edge[events$node[tmp[i]],1]),events$node[tmp[i]]))
        if (length(tmp2)>0) {tmp <- tmp[-tmp2]}
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
        out <- addspeciation(n,xmin,time=0,node.new,node=1,r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
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
            is.des <- T
            if (length(anc)>0) {
                if (node.new==anc) {is.des <- F}
            }
            if (type==1) {
                out <- addspeciation(n,xmin,time.new,node.new,ifelse(is.des,node,node.new),r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
                a <- integrand(t[node.new,1]-time.new)/integrand(t[node,1]-time)
            }
            if (type==2) {
                out <- addadaptation(n,xmin,time.new,node.new,ifelse(is.des,node,node.new),r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
                a <- 1
            }
            lik.new <- out$lik.new
            accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*a)
            if (runif(1)<accept) {
                lik.old$res[lik.new$idx] <- lik.new$res
                lik.old$likPO[lik.new$idx] <- lik.new$likPO
                ii <- which(!is.na(lik.new$dendistr[,1]))
                lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
                lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
                lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
                lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
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
            is.des <- T
            if (time.new < 0) {
                tmp <- which(edge[,2]==edge[node,1])
                if (length(tmp)>0) {
                    node.new <- tmp
                    time.new <- t[node.new,1]-t[node.new,2]+time.new
                    is.des <- F
                } else {
                    node.new <- node
                    time.new <- time-u
                }
            } else if (time.new > (t[node,1]-t[node,2])) {
                tmp <- which(edge[,1]==edge[node,2])
                if (length(tmp)>0) {
                    node.new <- sample(tmp,size=1)
                    time.new <- 0
                } else {
                    node.new <- node
                    time.new <- time-u
                }
            } else {
                node.new <- node
            }
            tmp <- events[-idx,]
            if (type==1) {
                out <- addspeciation(n,xmin,time.new,node.new,ifelse(is.des,node,node.new),r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
                a <- integrand(t[node.new,1]-time.new)/integrand(t[node,1]-time)
            }
            if (type==2) {
                out <- addadaptation(n,xmin,time.new,node.new,ifelse(is.des,node,node.new),r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
                a <- 1
            }
            lik.new <- out$lik.new
            accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*a)
            if (runif(1)<accept) {
            lik.old$res[lik.new$idx] <- lik.new$res
            lik.old$likPO[lik.new$idx] <- lik.new$likPO
            ii <- which(!is.na(lik.new$dendistr[,1]))
            lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
            lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
            lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
            lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
            lik.old$avec[ii] <- lik.new$avec[ii]
            lik.old$bvec[ii] <- lik.new$bvec[ii]
            lik.old$uvec[ii] <- lik.new$uvec[ii]
            lik.old$xvec[ii] <- lik.new$xvec[ii]
            events <- out$tmp
            }
        }
}
    }
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
    #updating n for the event
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
            a <- 0.5*prior.n(n.new)/prior.n(n)
        } else if (n==0.9) {
            n.new <- n-0.2
            a <- 0.5*prior.n(n.new)/prior.n(n)
        } else {
            n.new <- n+ifelse(runif(1)<=0.5,0.2,-0.2)
            a <- 1*prior.n(n.new)/prior.n(n)
        }
        out <- addspeciation(n.new,xmin,time,node,node,r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
    }
    if (type==2) {
        if (n==-(N-1)) {
            n.new <- n+1
            a <- 0.5
        } else if (n==(N-1)) {
            n.new <- n-1
            a <- 0.5
        } else {
            n.new <- n+ifelse(runif(1)<=0.5,ifelse(n==-1,2,1),ifelse(n==1,-2,-1))
            a <- 1
        }
        out <- addadaptation(n.new,xmin,time,node,node,r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
    }
    lik.new <- out$lik.new
    accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*a)
    if (runif(1)<accept) {
        lik.old$res[lik.new$idx] <- lik.new$res
        lik.old$likPO[lik.new$idx] <- lik.new$likPO
        ii <- which(!is.na(lik.new$dendistr[,1]))
        lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
        lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
        lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
        lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
        lik.old$avec[ii] <- lik.new$avec[ii]
        lik.old$bvec[ii] <- lik.new$bvec[ii]
        lik.old$uvec[ii] <- lik.new$uvec[ii]
        lik.old$xvec[ii] <- lik.new$xvec[ii]
        events <- out$tmp
    }
    #updating xmin for the event
    node <- events[idx,]$node
    type <- events[idx,]$type
    n <- events[idx,]$n
    xmin <- events[idx,]$xmin
    time <- events[idx,]$time
    if (xmin==1) {
        xmin.new <- xmin+1
        a <- 0.5
    } else if (xmin==(N-1)) {
        xmin.new <- xmin-1
        a <- 0.5
    } else {
        xmin.new <- xmin+ifelse(runif(1)<=0.5,1,-1)
        a <- 1
    }
    if (type==1) {
            tmp <- events[-idx,]
            out <- addspeciation(n,-xmin.new,time,node,node,r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
    }
    if (type==2) {
        out <- addadaptation(n,xmin.new,time,node,node,r,N,edge,edge.length,events=tmp,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
    }
    lik.new <- out$lik.new
    accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*a)
    if (runif(1)<accept) {
        lik.old$res[lik.new$idx] <- lik.new$res
        lik.old$likPO[lik.new$idx] <- lik.new$likPO
        ii <- which(!is.na(lik.new$dendistr[,1]))
        lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
        lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
        lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
        lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
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
    out <- addspeciation(n=n,xmin=sample.int(N-1,1),time,node,node,r,N,edge,edge.length,events,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
    lik.new <- out$lik.new
    if (time==0) {
        accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*A[1]/(k[1]+1)*ifelse(k[1]==0,0.5,1)*prior.n(n)/(N-1)/exp.r^2)
    } else {
        accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*A[1]/(k[1]+1)*ifelse(k[1]==0,0.5,1)*sp*integrand(t[node,1]-time)*prior.n(n)/(N-1)/exp.r^2)
    }
    if (runif(1)<accept) {
        lik.old$res[lik.new$idx] <- lik.new$res
        lik.old$likPO[lik.new$idx] <- lik.new$likPO
        ii <- which(!is.na(lik.new$dendistr[,1]))
            lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
            lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
            lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
            lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
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
    out <- addadaptation(n=sample.int(N-1,1)*sample(c(-1,1),1),xmin=sample.int(N-1,1),time,node,node,r,N,edge,edge.length,events,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old)
    lik.new <- out$lik.new
    accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*A[2]/(k[2]+1)*ifelse(k[2]==0,0.5,1)/tT2/(N-1)/(N-1)/exp.r^2)
    if (runif(1)<accept) {
        lik.old$res[lik.new$idx] <- lik.new$res
        lik.old$likPO[lik.new$idx] <- lik.new$likPO
        ii <- which(!is.na(lik.new$dendistr[,1]))
            lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
            lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
            lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
            lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
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
            nodeniche <- rootniche
            nodeback <- rootback
            nodeu <- rootu
            nodea <- roota
            nodeb <- rootb
            nodex <- rootx
            anc <- 0
        } else {
            anc <- as.character(anc)
            nodedistr <- lik.old$dendistr[anc,]
            nodedistro <- lik.old$dendistro[anc,]
            nodeniche <- lik.old$nichedistr[anc,]
            nodeback <- lik.old$backdistr[anc,]
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
            lik.new <- likcal(node=anc,root.included=F,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,nodeniche=nodeniche,nodeback=nodeback,d=d,bias=bias,area=area)
            accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*k[1]/A[1]*ifelse(k[1]==1,2,1)/prior.n(n)*(N-1)*exp.r^2)
        } else {
            tmp <- events[-idx,]
            lik.new <- likcal(node=node,root.included=T,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,nodeniche=nodeniche,nodeback=nodeback,d=d,bias=bias,area=area)
            accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*k[1]/A[1]*ifelse(k[1]==1,2,1)/sp/integrand(t[node,1]-time)/prior.n(n)*(N-1)*exp.r^2)
        }
        if (runif(1)<accept) {
            lik.old$res[lik.new$idx] <- lik.new$res
            lik.old$likPO[lik.new$idx] <- lik.new$likPO
            ii <- which(!is.na(lik.new$dendistr[,1]))
            lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
            lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
            lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
            lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
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
            nodeniche <- rootniche
            nodeback <- rootback
            nodeu <- rootu
            nodea <- roota
            nodeb <- rootb
            nodex <- rootx
        } else {
            anc <- as.character(anc)
            nodedistr <- lik.old$dendistr[anc,]
            nodedistro <- lik.old$dendistro[anc,]
            nodeniche <- lik.old$nichedistr[anc,]
            nodeback <- lik.old$backdistr[anc,]
            nodeu <- lik.old$uvec[anc]
            nodea <- lik.old$avec[anc]
            nodeb <- lik.old$bvec[anc]
            nodex <- lik.old$xvec[anc]
        }
        tmp <- events[-idx,]
        lik.new <- likcal(node=node,root.included=T,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,nodeniche=nodeniche,nodeback=nodeback,d=d,bias=bias,area=area)
        accept <- min(1,exp(sum(lik.old$res[lik.new$idx])-sum(lik.new$res))*k[2]/A[2]*ifelse(k[2]==1,2,1)*tT2*(N-1)*(N-1)*exp.r^2)
        if (runif(1)<accept) {
            lik.old$res[lik.new$idx] <- lik.new$res
            lik.old$likPO[lik.new$idx] <- lik.new$likPO
            ii <- which(!is.na(lik.new$dendistr[,1]))
            lik.old$dendistr[ii,] <- lik.new$dendistr[ii,]
            lik.old$dendistro[ii,] <- lik.new$dendistro[ii,]
            lik.old$nichedistr[ii,] <- lik.new$nichedistr[ii,]
            lik.old$backdistr[ii,] <- lik.new$backdistr[ii,]
            lik.old$avec[ii] <- lik.new$avec[ii]
            lik.old$bvec[ii] <- lik.new$bvec[ii]
            lik.old$uvec[ii] <- lik.new$uvec[ii]
            lik.old$xvec[ii] <- lik.new$xvec[ii]
            events <- tmp
            k[2] <- k[2]-1
        }
    }
    #updating As
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
    #update roota
    if (runif(1)<0.05) {
    u <- rnorm(1)
    tmpa <- roota+u*tuning.roota
    rootnichetmp <- exp(-((0:(N-1))/(exp(tmpa)+1e-6))^(exp(rootb)+1e-6))
    rootdistrotmp <- c(rep(1,N),rootnichetmp)*rootback
    rootdistrotmp <- rootdistrotmp/sum(rootdistrotmp)
    rootdistrtmp <- rootdistrotmp[(N+1):(2*N)]
    rootdistrtmp[1] <- rootdistrtmp[1]+sum(rootdistrotmp[1:N])
    lik.new <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=tmpa,nodeb=rootb,nodex=rootx,nodedistr=rootdistrtmp,nodedistro=rootdistrotmp,nodeniche=rootnichetmp,nodeback=rootback,d=d,bias=bias,area=area)
    a <- exp(sum(lik.old$res)-sum(lik.new$res))*dnorm(tmpa,mean=norm.roota.m,sd=norm.roota.sd)/dnorm(roota,mean=norm.roota.m,sd=norm.roota.sd)
    accept <- min(1,ifelse(is.na(a),0,a))
    if (runif(1)<accept) {
        lik.old <- lik.new
        roota <- tmpa
        rootdistr <- rootdistrtmp
        rootdistro <- rootdistrotmp
        rootniche <- rootnichetmp
    }
    #updating rootb
    u <- rnorm(1)
    tmpb <- rootb+u*tuning.rootb
    rootnichetmp <- exp(-((0:(N-1))/(exp(roota)+1e-6))^(exp(tmpb)+1e-6))
    rootdistrotmp <- c(rep(1,N),rootnichetmp)*rootback
    rootdistrotmp <- rootdistrotmp/sum(rootdistrotmp)
    rootdistrtmp <- rootdistrotmp[(N+1):(2*N)]
    rootdistrtmp[1] <- rootdistrtmp[1]+sum(rootdistrotmp[1:N])
    lik.new <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=tmpb,nodex=rootx,nodedistr=rootdistrtmp,nodedistro=rootdistrotmp,nodeniche=rootnichetmp,nodeback=rootback,d=d,bias=bias,area=area)
    a <- exp(sum(lik.old$res)-sum(lik.new$res))*dnorm(tmpb,mean=norm.rootb.m,sd=norm.rootb.sd)/dnorm(rootb,mean=norm.rootb.m,sd=norm.rootb.sd)
    accept <- min(1,ifelse(is.na(a),0,a))
    if (runif(1)<accept) {
        lik.old <- lik.new
        rootb <- tmpb
        rootdistr <- rootdistrtmp
        rootdistro <- rootdistrotmp
        rootniche <- rootnichetmp
    }
    #udating rootsd
    u <- runif(1)
    tmpsd <- sd*exp((u-0.5)*tuning.sd)
    rootbacktmp <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=tmpsd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=tmpsd)
    rootbacktmp <- rootbacktmp/sum(rootbacktmp)
    rootdistrotmp <- c(rep(1,N),rootniche)*rootbacktmp
    rootdistrotmp <- rootdistrotmp/sum(rootdistrotmp)
    rootdistrtmp <- rootdistrotmp[(N+1):(2*N)]
    rootdistrtmp[1] <- rootdistrtmp[1]+sum(rootdistrotmp[1:N])
    lik.new <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=tmpsd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistrtmp,nodedistro=rootdistrotmp,nodeniche=rootniche,nodeback=rootbacktmp,d=d,bias=bias,area=area)
    a <- exp(sum(lik.old$res)-sum(lik.new$res))*dlnorm(tmpsd,meanlog=lnorm.sd.m,sdlog=lnorm.sd.sd)/dlnorm(sd,meanlog=lnorm.sd.m,sdlog=lnorm.sd.sd)*exp(tuning.sd*(u-0.5))
    accept <- min(1,ifelse(is.na(a),0,a))
    if (runif(1)<accept) {
        lik.old <- lik.new
        sd <- tmpsd
        rootdistr <- rootdistrtmp
        rootdistro <- rootdistrotmp
        rootback <- rootbacktmp
    }
    #updating rootu
    u <- rnorm(1)
    tmpu <- rootu+u*tuning.rootu
    rootbacktmp <- pnorm(seq(-(N-1),N,1),mean=tmpu,sd=sd)-pnorm(seq(-N,N-1,1),mean=tmpu,sd=sd)
    rootbacktmp <- rootbacktmp/sum(rootbacktmp)
    rootdistrotmp <- c(rep(1,N),rootniche)*rootbacktmp
    rootdistrotmp <- rootdistrotmp/sum(rootdistrotmp)
    rootdistrtmp <- rootdistrotmp[(N+1):(2*N)]
    rootdistrtmp[1] <- rootdistrtmp[1]+sum(rootdistrotmp[1:N])
    lik.new <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=tmpu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistrtmp,nodedistro=rootdistrotmp,nodeniche=rootniche,nodeback=rootbacktmp,d=d,bias=bias,area=area)
    a <- exp(sum(lik.old$res)-sum(lik.new$res))*dnorm(tmpu,mean=norm.rootu.m,sd=norm.rootu.sd)/dnorm(rootu,mean=norm.rootu.m,sd=norm.rootu.sd)
    accept <- min(1,ifelse(is.na(a),0,a))
    if (runif(1)<accept) {
        lik.old <- lik.new
        rootu <- tmpu
        rootdistr <- rootdistrtmp
        rootdistro <- rootdistrotmp
        rootback <- rootbacktmp
    }
    #updating d
    u <- runif(1)
    tmp <- d*exp((u-0.5)*tuning.d)
    lik.new <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,nodeniche=rootniche,nodeback=rootback,d=tmp,bias=bias,area=area)
    a <- exp(sum(lik.old$res)-sum(lik.new$res))*dexp(tmp,rate=exp.d)/dexp(d,rate=exp.d)*exp(tuning.d*(u-0.5))
    accept <- min(1,ifelse(is.na(a),0,a))
    if (runif(1)<accept) {
        lik.old <- lik.new
        d <- tmp
    }
    #updating r
    if (dim(events)[1]>0) {
    u <- runif(1)
    tmp <- r*exp((u-0.5)*tuning.r)
    eventstmp <- events
    eventstmp$a <- tmp
    eventstmp$b <- tmp
    lik.new <- likcal(node=0,root.included=F,events=eventstmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,nodeniche=rootniche,nodeback=rootback,d=d,bias=bias,area=area)
    a <- exp(sum(lik.old$res)-sum(lik.new$res))*dexp(tmp,rate=exp.r)/dexp(r,rate=exp.r)*exp(tuning.r*(u-0.5))
    accept <- min(1,ifelse(is.na(a),0,a))
    if (runif(1)<accept) {
        lik.old <- lik.new
        r <- tmp
        events <- eventstmp
    }
    }
    }
if (tt %% 1000 ==0) {
    prior <- dnorm(roota,mean=norm.roota.m,sd=norm.roota.sd,log=T)+dnorm(rootb,mean=norm.rootb.m,sd=norm.rootb.sd,log=T)+dlnorm(sd,meanlog=lnorm.sd.m,sdlog=lnorm.sd.sd,log=T)+dexp(d,rate=exp.d,log=T)+dnorm(rootu,mean=norm.rootu.m,sd=norm.rootu.sd,log=T)+dexp(r,rate=exp.r,log=T)+dexp(A[1],rate=exp.A[1],log=T)+dexp(A[2],rate=exp.A[2],log=T)
    suppressWarnings(write.table(events,eventsfilename,append=T,col.names=T))
    write.table(t(as.numeric(c(tt,sum(lik.old$res),prior,roota,rootb,rootu,sd,d,r,A[1],A[2],bias))),postfilename,append=T,col.names=F)
    print(tt)
    save.image(filename)
}
}
}
