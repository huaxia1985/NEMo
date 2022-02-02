addadaptation <- function (n,xmin,time,node.new,node,r,N,edge,edge.length,events,rootdistr,rootdistro,rootniche,rootback,rootu,roota,rootb,rootx,nspecies,PO,PA,sd,d,bias,area,lik.old) {
type <- 2
a <- r
b <- r
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
tmp <- events
tmp[nrow(tmp)+1,] <- c(node.new,type,time,a,b,n,xmin)
tmp <- tmp[with(tmp,order(node,time)),]
lik.new <- likcal(node=node,root.included=T,events=tmp,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=nodeu,sd=sd,nodea=nodea,nodeb=nodeb,nodex=nodex,nodedistr=nodedistr,nodedistro=nodedistro,nodeniche=nodeniche,nodeback=nodeback,d=d,bias=bias,area=area)
out <- list(lik.new=lik.new,tmp=tmp)
out
}
