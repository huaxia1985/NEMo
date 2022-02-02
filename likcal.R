likcal <- function (node,root.included,events,edge,edge.length,nspecies,PO,PA,N=10,nodeu,sd,nodea,nodeb,nodex,nodedistr,nodedistro,nodeniche,nodeback,d,bias,area) {
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
if (node==0) {
des <- which(edge[,1]==edge[1,1])
} else {
des <- which(edge[,1]==edge[node,2])
}
if (root.included) {
res <- upd(events=events[events$node==node,],edget=edge.length[node],dendistrt=nodedistr,dendistrto=nodedistro,nichedistrt=nodeniche,backdistrt=nodeback,a=nodea,b=nodeb,u=nodeu,x=nodex,sd,N,d)
    nodedistr <- dendistr[node,] <- res$dendistrt
    nodedistro <- dendistro[node,]  <- res$dendistrto
    nodeniche <- nichedistr[node,] <- res$nichedistrt
    nodeback <- backdistr[node,] <- res$backdistrt
    nodea <- avec[node] <- res$a
    nodeb <- bvec[node] <- res$b
    nodeu <- uvec[node] <- res$u
    nodex <- xvec[node] <- res$x
}
if (length(des)>0) {
      j <- des[1]
      res <- upd(events[events$node==j,],edget=edge.length[j],dendistrt=nodedistr,dendistrto=nodedistro,nichedistrt=nodeniche,backdistrt=nodeback,a=nodea,b=nodeb,u=nodeu,x=nodex,sd,N,d)
      dendistr[j,] <- res$dendistrt
      dendistro[j,] <- res$dendistrto
      nichedistr[j,] <- res$nichedistrt
      backdistr[j,] <- res$backdistrt
      avec[j] <- res$a
      bvec[j] <- res$b
      uvec[j] <- res$u
      xvec[j] <- res$x
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
            j<-j+1
                if (j>dim(edge)[1]) {break}
      }
      j <- des[2]
      res <- upd(events[events$node==j,],edget=edge.length[j],dendistrt=nodedistr,dendistrto=nodedistro,nichedistrt=nodeniche,backdistrt=nodeback,a=nodea,b=nodeb,u=nodeu,x=nodex,sd,N,d)
      dendistr[j,] <- res$dendistrt
      dendistro[j,] <- res$dendistrto
      nichedistr[j,] <- res$nichedistrt
      backdistr[j,] <- res$backdistrt
      avec[j] <- res$a
      bvec[j] <- res$b
      uvec[j] <- res$u
      xvec[j] <- res$x
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
            j<-j+1
                if (j>dim(edge)[1]) {break}
       }
       }
###
rescal <- function (i,edge,PA,PO,BK,dendistr,nichedistr,backdistr,bias,area,N) {
j <-which(edge[,2]==i)
       likPA <- 0
       if (!is.null(PA[[i]])) {
           likPA <- sum(log(1-exp(-dendistr[j,PA[[i]]$x]))[PA[[i]]$y==1])-sum((1-PA[[i]]$y)*dendistr[j,PA[[i]]$x])
       }
       likPO <- 0
       if (!is.null(PO[[i]])) {
           backdistro <- backdistr[,(N+1):(2*N)]
           backdistro[,1] <- rowSums(backdistr[,1:(N+1)])
           likBK <- sum(log(backdistro[j,BK[[i]]$x]))
           back <- backdistro[j,]*nichedistr[j,]
           if (!is.null(bias)) {
               likPO <- sum(log(dendistr[j,PO[[i]]$x]))+sum(colSums(bias*t(PO[[i]]$bias)))-sum(exp(colSums(bias*t(BK[[i]]$bias)))*back[BK[[i]]$x])*area[i]
           } else {
               likPO <- sum(log(dendistr[j,PO[[i]]$x]))-area*sum(back)
           }
       res <- -likPA-likPO-likBK
       c(res,likPO)
       }
}
idx <- which(!is.na(dendistr[as.character(c(1:nspecies)),1]))
res <- sapply(idx,rescal,edge,PA,PO,BK,dendistr,nichedistr,backdistr,bias,area,N)
if (length(res)>2) {
likPO <- res[2,]
res <- res[1,]
} else {
likPO <- res[2]
res <- res[1]
}
       res[is.na(res)] <- 10^7
       res[is.infinite(res)] <- 10^7
       likPO[is.na(likPO)] <- -10^7
likPO[is.infinite(likPO)] <- -10^7
       out <- list(res=res,idx=idx,dendistr=dendistr,dendistro=dendistro,backdistr=backdistr,nichedistr=nichedistr,avec=avec,bvec=bvec,uvec=uvec,xvec=xvec,likPO=likPO)
       out
    }
