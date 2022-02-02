##### input: PO, PA, BK, timetree, area (D), N (number of bins).
#PO is a list of presence-only data, where each element is each tip species in the order corresponding to timetree$tip.label. In each element, PO[[i]]$x is the value of the environmental variale of interest at each location. PO[[i]]$bias is the value of the varialbes included to correct for sampling bias at each location
#PA is a list of presence-absence data, where each element is each tip species in the order corresponding to timetree$tip.label. In each element, PA[[i]]$x is the value of the environmental variale of interest at each location. PA[[i]]$y is a binary code recording whether each location is a presence or an absence
#BK is a list of background data, where each element is each tip species in the order corresponding to timetree$tip.label. In each element, BK[[i]]$x is the value of the environmental variale of interest at each location. BK[[i]]$bias is the value of the varialbes included to correct for sampling bias at each location
load(paste0(getwd(),"/arid"))
#or load(paste0(getwd(),"/salt")) for salt tolerance

##### output file name
# for one MCMC chain
eventsfilename <- paste0(getwd(),"/case study/events1")
postfilename <- paste0(getwd(),"/case study/post1")
workfilename <- paste0(getwd(),"/case study/work1")

##### set values
r <- 0.03
d <- 1
rootx <- 1

norm.roota.sd <- norm.rootb.sd <- 1
norm.rootu.sd <- 3
lnorm.sd.sd <- 0.5
exp.d <- 0.5
exp.r <- 10
exp.A <- c(10,10)
exp.bias1 <- exp.bias2 <- 0.1

tuning.roota <- tuning.rootb <- 0.1
tuning.rootu <- 1
tuning.sd <- 0.1
tuning.r <- 0.1
tuning.d <- 0.1
tuning.A <- c(0.1,0.1)
tuning.bias1 <- tuning.bias2 <- 0.1

#event occurrence rate and time prior
library(ape)
library(diversitree)
library(phytools)
lik <- make.bd(timetree)
fit <- find.mle(lik,c(0.1,0.03),method="subplex")
sp <- coef(fit)[1]
mu <- coef(fit)[2]
t <- nodeHeights(timetree)
t <- max(t)-t
integrand <- function (tt) {mu*(1-exp(-(sp-mu)*tt))/(sp-mu*exp(-(sp-mu)*tt))}
edge.prob <- sapply(1:length(timetree$edge.length),function (i) integrate(integrand,lower=t[i,2],upper=t[i,1])$value)
tT1 <- sum(edge.prob)
edge.prob1 <- edge.prob/tT1
tT2 <- sum(timetree$edge.length)
edge.prob2 <- timetree$edge.length/tT2
edge <- timetree$edge
edge.length <- timetree$edge.length
tuning.t <- min(edge.length)
nspecies <- length(timetree$tip.label)
A <- c(0.1/(nspecies-1+sp*tT1),0.1)

#get ML estiamtes of ancestral niches and coefficient for sampling bias
events <- matrix(NA,0,7)
events <- data.frame(events)
names(events) <- c("node","type","time","a","b","n","xmin")
fn <- function (pars) {
	bias <- c(-exp(pars[1]),-exp(pars[2]))
	roota <- pars[3]
	rootb <- pars[4]
	rootu <- exp(pars[5])
	sd <- exp(pars[6])
	rootniche <- exp(-((0:(N-1))/(exp(roota)+1e-6))^(exp(rootb)+1e-6))
	rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
	rootback <- rootback/sum(rootback)
	rootdistro <- c(rep(1,N),rootniche)*rootback
	rootdistro <- rootdistro/sum(rootdistro)
	rootdistr <- rootdistro[(N+1):(2*N)]
	rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
	tmp <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,nodeniche=rootniche,nodeback=rootback,d=d,bias=bias,area=area)
	sum(tmp$res)
}
bias <- optim(par=c(0,0,log(3),log(2),log(1),log(2)),fn)
roota <- norm.roota.m <- bias$par[3]
rootb <- norm.rootb.m <- bias$par[4]
rootu <- norm.rootu.m <- exp(bias$par[5])
sd <- lnorm.sd.m <- exp(bias$par[6])
bias <- bias.m <- -exp(bias$par[c(1,2)])

#calculate initial likelihood
rootniche <- exp(-((0:(N-1))/(exp(roota)+1e-6))^(exp(rootb)+1e-6))
rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
rootback <- rootback/sum(rootback)
rootdistro <- c(rep(1,N),rootniche)*rootback
rootdistro <- rootdistro/sum(rootdistro)
rootdistr <- rootdistro[(N+1):(2*N)]
rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
lik.old <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,nodeniche=rootniche,nodeback=rootback,d=d,bias=bias,area=area)
ngen <- 3*10^6
k <- c(0,0)

#### start MCMC
niche_mcmc(ngen)

#### summarize results
#read files
output <- read.table(postfilename,header=F,sep=" ")
output <- output[,-1]
names(output) <- c("iter","nll","prior","roota","rootb","rootu","sd","d","r","A1","A2","bias1","bias2")
out$posterior <- -out$nll+out$prior
plot(c(1:dim(output)[1]),output$posterior,type="l")
burnin <- 600 #picked by the plot
output <- output[-c(1:burnin),]
tmp <- read.table(eventsfilename,header=F,sep=" ",fill=T,stringsAsFactors=F)
idx <- which(is.na(tmp[,8]))
nsamples <- length(idx)-burnin
events.list <- vector("list",nsamples)
for (i in (burnin+1):length(idx)) {
    tmp2 <- tmp[(idx[i]+1):(ifelse(i<length(idx),idx[i+1],dim(tmp)[1]+1)-1),-1]
    names(tmp2) <- c("node","type","time","a","b","n","xmin")
    events.list[[i-burnin]] <- data.frame(apply(tmp2,2,function (x) as.numeric(x)))
}

#calculate niches
lik.list <- vector("list",nsample)
for (i in 1:nsamples) {
    roota <- output$roota[i]
    rootb <- output$rootb[i]
    rootu <- output$rootu[i]
    sd <- output$sd[i]
    d <- output$d[i]
    bias <- c(output$bias1[i],output$bias2[i])
    rootniche <- exp(-((0:(N-1))/(exp(roota)+1e-6))^(exp(rootb)+1e-6))
    rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
    rootback <- rootback/sum(rootback)
    rootdistro <- c(rep(1,N),rootniche)*rootback
    rootdistro <- rootdistro/sum(rootdistro)
    rootdistr <- rootdistro[(N+1):(2*N)]
    rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
    lik.list[[i]] <- likcal(node=0,root.included=F,events=events.list[[i]],edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,nodeniche=rootniche,nodeback=rootback,d=d,bias=bias,area=area)
}

#plot history of niche evolution
#this example uses the sample with the highst posterior probability
library(phytools)
idx <- which(out$posterior==max(out$posterior))
events <- events.list[[idx]]
best.avec <- lik.list[[idx]]$avec
best.avec <- c(out$roota[idx],best.avec)
names(best.avec)[1] <- "506"
#changes in parameter alpha in fundamental niche along each branch
best.avec <- best.avec[as.character(timetree$edge[,2])]-best.avec[as.character(timetree$edge[,1])]
colfunc <- colorRampPalette(c("blue","lightgrey"))
colfunc2 <- colorRampPalette(c("lightgrey","red"))
tmp <- hist(best.avec,breaks=20)
n <- sum(tmp$mids<0)
idx <- findInterval(best.avec[best.avec<0],tmp$breaks[1:(n+1)])
n2 <- sum(tmp$mids>=0)
idx2 <-findInterval(best.avec[best.avec>=0],tmp$breaks[(n+1):(n+n2+1)])
col <- colfunc(n)[idx]
col2 <- colfunc2(n2+1)[idx2]
colr <- numeric(length(timetree$edge.length))
colr[best.avec<0] <- col
colr[best.avec>=0] <- col2
plot(timetree,type="phylogram",show.tip.label=F,edge.color=colr,edge.width=6,no.margin=T,direction="downwards")
#plot color scale
legend_image <- as.raster(matrix(c(colfunc2(n2+1)[(n2+1):1],colfunc(n)[n:1]), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.5, y = seq(0,1,l=length(tmp$breaks)), labels = tmp$breaks)
rasterImage(legend_image, 0, 0, 1,1)

#plot adaptation event
u.edge <- events[events$type==2,"node"]
u.time <- events[events$type==2,"time"]
u.n <- events[events$type==2,"n"]
tmp <- seq(-4,4,1)
u.n[u.n<(-3)] <- -4
u.n[u.n>3] <- 4
colfunc <- colorRampPalette(c("lightgrey","blue"))
colfunc2 <- colorRampPalette(c("lightgrey","red"))
col <- colfunc(5)[abs(u.n[u.n<0])+1]
col2 <- colfunc2(5)[u.n[u.n>0]+1]
colr <- numeric(length(u.edge))
colr[u.n<0] <- col
colr[u.n>0] <- col2
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xx <- lastPP$xx[lastPP$edge[u.edge,2]]
yy <- lastPP$yy[lastPP$edge[u.edge,1]]
points(xx,yy-u.time,pch=24,bg=colr,cex=1.5)

#plot speciation event
u.edge <- events[events$type==1,"node"]
u.time <- events[events$type==1,"time"]
u.n <- events[events$type==1,"n"]
tmp <- seq(-1,1,0.2)
colfunc <- colorRampPalette(c("blue","lightgrey"))
colfunc2 <- colorRampPalette(c("lightgrey","red"))
col <- colfunc(6)[findInterval(u.n[u.n<0],tmp[1:6])]
col2 <- colfunc2(6)[findInterval(u.n[u.n>0],tmp[6:11])]
colr <- numeric(length(u.edge))
colr[u.n<0] <- col
colr[u.n>0] <- col2
xx <- lastPP$xx[lastPP$edge[u.edge,2]]
yy <- lastPP$yy[lastPP$edge[u.edge,1]]
points(xx,yy-u.time,pch=21,bg=colr,cex=1.5)

#plot ancestral occupied niches, to plot the fundamental niches, replace dendistr with nichedistr, to plot the available niche, replace dendistr with backdistr
heatmap.2(1-lik.list[[idx]]$dendistr[as.character(timetree$edge[timetree$edge[,2]>505,2]),],dendrogram='none',Rowv=F,Colv=F, trace='none', scale='none', key=T, key.title=NA, key.ylab=NA, key.xlab=NA, density.info='none', keysize=0.1,
lmat=rbind( c(3, 4), c(2,1)), lhei=c(0.3, 4), lwid=c(0.1,2))
#plot tip occupied niches 
heatmap.2(1-lik.list[[idx]]$dendistr[as.character(timetree$edge[timetree$edge[,2]<=505,2]),],dendrogram='none',Rowv=F,Colv=F, trace='none', scale='none', key=T, key.title=NA, key.ylab=NA, key.xlab=NA, density.info='none', keysize=0.1, margins=c(1,1),
lmat=rbind( c(3, 4), c(2,1)), lhei=c(0.3, 4), lwid=c(0.1,2))

#function to calculate confidence interval of the occurence density of events
GetConfidenceInterval <- function(x, alpha=0.05) {
        if (x == 0) {
            L <- 0
            U <- - log(alpha)
        } else {
        fl <- function(lambda) {
                ppois(x, lambda) - alpha/2
        }
        lambda.upper <- x + (1 + sqrt(2*x*alpha + 1)) / alpha
        U <- uniroot(fl, c(x, lambda.upper))$root
        fu <- function(lambda) {
                ppois(x - 1, lambda, lower.tail=FALSE) - alpha/2
        }
        L <- uniroot(fu, c(0, x))$root
        }
return(c(L,U))
}

#function to calculate confidence interval of difference in the occurence density of events
GetConfidenceIntervaldiff <- function (x1,x2,alpha=0.05) {
    CI1 <- GetConfidenceInterval(x1,alpha)
    CI2 <- GetConfidenceInterval(x2,alpha)
    L <- x1-x2-sqrt((x1-CI1[1])^2+(CI2[2]-x2)^2)
    U <- x1-x2+sqrt((CI1[2]-x1)^2+(x2-CI2[1])^2)
    return(c(L,U))
}

#calculate the occurence density of events in different drying periods
idx1 <- which(events$type==1)
idx2p <- which((events$type==2)*(events$n>0)==1)
idx2n <- which((events$type==2)*(events$n<0)==1)
time <- t[events$node[idx1],1]
time2p <- t[events$node[idx2p],1]-events$time[idx2p]
time2n <- t[events$node[idx2n],1]-events$time[idx2n]

idx1 <- timetree$edge[events$node[idx1],1]
time1 <- idx1[which(time<2.6)]
time2 <- idx1[which((time>2.6)*(time<5)==1)]
time3 <- idx1[which((time>5)*(time<14)==1)]
time4 <- idx1[which((time>14)*(time<23)==1)]

idx2p <- timetree$edge[events$node[idx2p],1]
time2p1 <- idx2p[which(time2p<2.6)]
time2p2 <- idx2p[which((time2p>2.6)*(time2p<5)==1)]
time2p3 <- idx2p[which((time2p>5)*(time2p<14)==1)]
time2p4 <- idx2p[which((time2p>14)*(time2p<23)==1)]

idx2n <- timetree$edge[events$node[idx2n],1]
time2n1 <- idx2n[which(time2n<2.6)]
time2n2 <- idx2n[which((time2n>2.6)*(time2n<5)==1)]
time2n3 <- idx2n[which((time2n>5)*(time2n<14)==1)]
time2n4 <- idx2n[which((time2n>14)*(time2n<23)==1)]

data <- rbind(c(length(time1)/c[1], GetConfidenceInterval(length(time1),alpha=0.05)/c[1],length(time2p1)/c[1], GetConfidenceInterval(length(time2p1),alpha=0.05)/c[1],length(time2n1)/c[1], GetConfidenceInterval(length(time2n1),alpha=0.05)/c[1]),
    c(length(time2)/c[2], GetConfidenceInterval(length(time2),alpha=0.05)/c[2],length(time2p2)/c[2], GetConfidenceInterval(length(time2p2),alpha=0.05)/c[2],length(time2n2)/c[2], GetConfidenceInterval(length(time2n2),alpha=0.05)/c[2]),
    c(length(time3)/c[3], GetConfidenceInterval(length(time3),alpha=0.05)/c[3],length(time2p3)/c[3], GetConfidenceInterval(length(time2p3),alpha=0.05)/c[3],length(time2n3)/c[3], GetConfidenceInterval(length(time2n3),alpha=0.05)/c[3]),
    c(length(time4)/c[4], GetConfidenceInterval(length(time4),alpha=0.05)/c[4],length(time2p4)/c[4], GetConfidenceInterval(length(time2p4),alpha=0.05)/c[4],length(time2n4)/c[4], GetConfidenceInterval(length(time2n4),alpha=0.05)/c[4]),
    c(length(time5)/c[5], GetConfidenceInterval(length(time5),alpha=0.05)/c[5],length(time2p5)/c[5], GetConfidenceInterval(length(time2p5),alpha=0.05)/c[5],length(time2n5)/c[5], GetConfidenceInterval(length(time2n5),alpha=0.05)/c[5]))

GetConfidenceIntervaldiff(length(time2p1),length(time2n1))/c[1]
GetConfidenceIntervaldiff(length(time2p2),length(time2n2))/c[2]
GetConfidenceIntervaldiff(length(time2p3),length(time2n3))/c[3]
GetConfidenceIntervaldiff(length(time2p4),length(time2n4))/c[4]

par(mfrow=c(1,2))
plot(x=1:5,y=data[,1],ylim=c(0,0.3),xaxt="n",xlab="Time intervals",ylab="Occurrence rate of speciation event")
axis(side=1,at=1:5,labels=c("0-2.6Ma", "2.6-5Ma","5-14Ma","14-23Ma","23Ma-root"))
segments(x0=c(1:5),x1=c(1:5),y0=data[,2],y1=data[,3])
plot(x=1:5+0.1,y=data[,4],xlim=c(0.8,5),ylim=c(0,0.2),xaxt="n",xlab="Time intervals",ylab="Occurrence rate of adaptation events",col="red")
axis(side=1,at=1:5,labels=c("0-2.6Ma", "2.6-5Ma","5-14Ma","14-23Ma","23Ma-root"))
segments(x0=c(1:5)+0.1,x1=c(1:5)+0.1,y0=data[,5],y1=data[,6],col="red")
points(x=1:5-0.1,y=data[,7],col="blue")
segments(x0=c(1:5)-0.1,x1=c(1:5)-0.1,y0=data[,8],y1=data[,9],col="blue")
legend(x="topright",legend=c("with environmetnal change to higher aridity", "with environmetnal change to lower aridity"),pch=1,col=c("red","blue"))
