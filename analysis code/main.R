##### input: PO, PA, BK, timetree, area (D), N (number of bins).
#PO is a list of presence-only data, where each element is each tip species in the order corresponding to timetree$tip.label. In each element, PO[[i]]$x is the value of the environmental variale of interest at each location. PO[[i]]$bias is the value of the varialbes included to correct for sampling bias at each location
#PA is a list of presence-absence data, where each element is each tip species in the order corresponding to timetree$tip.label. In each element, PA[[i]]$x is the value of the environmental variale of interest at each location. PA[[i]]$y is a binary code recording whether each location is a presence or an absence
#BK is a list of background data, where each element is each tip species in the order corresponding to timetree$tip.label. In each element, BK[[i]]$x is the value of the environmental variale of interest at each location. BK[[i]]$bias is the value of the varialbes included to correct for sampling bias at each location
load("arid")
#or load("salt")

##### output file name
# for one MCMC chain
eventsfilename <- "events1"
postfilename <- "post1"

##### set initial values
#adaptation rate
expinv.a <- expinv.b <- 0.03
exp.a <- exp.b <- 10
tuning.a <- tuning.b <- 0.1
#dispersal rate
d <- 1
exp.d <- 0.1
tuning.d <- 0.1

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
exp.A <- c(10,10)
tuning.A <- c(0.1,0.1)

#get ML estiamtes of ancestral niches and coefficient for sampling bias
lnorm.roota <- lnorm.rootb <- 0.5
tuning.roota <- tuning.rootb <- 0.1
norm.rootu <- 3
tuning.rootu <- 0.1
lnorm.sd <- 0.5
tuning.sd <- 0.1
rootx <- 1
exp.bias1 <- exp.bias2 <- 0.1
tuning.bias1 <- tuning.bias2 <- 0.1
events <- matrix(NA,0,7)
events <- data.frame(events)
names(events) <- c("node","type","time","a","b","n","xmin")
fn <- function (pars) {
	bias <- c(-exp(pars[1]),-exp(pars[2]))
	roota <- pars[3]
	rootb <- pars[4]
	rootu <- exp(pars[5])
	sd <- exp(pars[6])
	rootniche <- sapply((seq(1:N)/(exp(roota)+0.1))^(exp(rootb)+0.1),gammainc,a=1/(exp(rootb)+0.1))
	rootniche <- (rootniche-c(0,rootniche[-length(rootniche)]))/2/gamma(1/(exp(rootb)+0.1))
	rootniche <- rootniche/rootniche[1]
	rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
	rootdistro <- c(rep(1,N),rootniche)*rootback
	rootdistro <- rootdistro/sum(rootdistro)
	rootdistr <- rootdistro[(N+1):(2*N)]
	rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
	tmp <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,d=d,bias=bias,area=area)
	sum(tmp$res)
}
bias <- optim(par=c(0,0,log(3),log(2),log(1),log(2)),fn)
roota <- lnorm.roota.m <- bias$par[3]
rootb <- <- lnorm.rootb.m <- bias$par[4]
rootu <- norm.rootu.m <- exp(bias$par[5])
sd <- lnorm.sd.m <- exp(bias$par[6])
bias <- bias.m <- -exp(bias$par[c(1,2)])

#calculate initial likelihood
rootniche <- sapply((seq(1:N)/(exp(roota)+0.1))^(exp(rootb)+0.1),gammainc,a=1/(exp(rootb)+0.1))
rootniche <- (rootniche-c(0,rootniche[-length(rootniche)]))/2/gamma(1/(exp(rootb)+0.1))
rootniche <- rootniche/rootniche[1]
rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
rootdistro <- c(rep(1,N),rootniche)*rootback
rootdistro <- rootdistro/sum(rootdistro)
rootdistr <- rootdistro[(N+1):(2*N)]
rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
lik.old <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,d=d,bias=bias,area=area)
ngen <- 10^6
k <- c(0,0)

#### start MCMC
niche_mcmc(ngen)

#### summarize results
#read files
output <- read.table(postfilename,header=F,sep=" ")
output <- output[,-1]
names(output) <- c("iter","nll","prior","roota","rootb","rootu","sd","d","expinv.a","expinv.b","A1","A2","bias1","bias2")
plot(c(1:dim(output)[1]),output[1:dim(output)[1],2]-output[1:dim(output)[1],3],type="l")
burnin <- 600 #picked by the plot
output <- output[-c(1:burnin),]
tmp <- read.table(eventsfilename,header=F,sep=" ",fill=T,stringsAsFactors=F)
idx <- which(is.na(tmp[,8]))
nsamples <- length(idx)-burnin
events <- vector("list",nsamples)
for (i in (burnin+1):length(idx)) {
    tmp2 <- tmp[(idx[i]+1):(ifelse(i<length(idx),idx[i+1],dim(tmp)[1]+1)-1),-1]
    names(tmp2) <- c("node","type","time","a","b","n","xmin")
    events[[i-burnin]] <- data.frame(apply(tmp2,2,function (x) as.numeric(x)))
}
#amat records the fundamental niche of each tip species; event1mat and event2mat records the number of speciation event and adaptation event inferered on each edge; event1matn records the sum of p values of the inferred speciation event on each edge; event2matn records the sum of n values of the inferred adaptation event on each edge; t1 records the time of each inferred speciation event; t2 records the time of each inferred adaptation event; N2 records the n value of each inferred adaptation event; 
amat <- matrix(NA,nsamples,nspecies)
event1mat <- matrix(0,nsamples,dim(edge)[1])
event2mat <- matrix(0,nsamples,dim(edge)[1])
event1matn <- matrix(0,nsamples,dim(edge)[1])
event2matn <- matrix(0,nsamples,dim(edge)[1])
t1 <- NULL
t2 <- NULL
N2 <- NULL
for (i in 1:nsamples) {
    roota <- output[i,4]
    rootb <- output[i,5]
    rootu <- output[i,6]
    sd <- output[i,7]
    d <- output[i,8]
    bias <- output[i,c(13,14)]
    rootniche <- sapply((seq(1:N)/(exp(roota)+0.1))^(exp(rootb)+0.1),gammainc,a=1/(exp(rootb)+0.1))
    rootniche <- (rootniche-c(0,rootniche[-length(rootniche)]))/2/gamma(1/(exp(rootb)+0.1))
    rootniche <- rootniche/rootniche[1]
    rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
    rootdistro <- c(rep(1,N),rootniche)*rootback
    rootdistro <- rootdistro/sum(rootdistro)
    rootdistr <- rootdistro[(N+1):(2*N)]
    rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
    lik <- likcal(node=0,root.included=F,events=events[[i]],edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,d=d,bias=bias,area=area)
    amat[i,] <- lik$avec[as.character(1:nspecies)]
    tmp <- which(events[[i]]$type==1)
    t1 <- c(t1,t[events[[i]]$node[tmp],1]-events[[i]]$time[tmp])
    tmp <- table(events[[i]]$node[tmp])
    tmp <- data.frame(tmp)
    event1mat[i,as.numeric(as.character(tmp$Var1))] <- tmp$Freq
    event1matn[i,as.numeric(as.character(tmp$Var1))] <- sapply(as.numeric(as.character(tmp$Var1)),function (z) sum(events[[i]]$n[which((events[[i]]$node==z)*(events[[i]]$type==1)==1)]))
    tmp <- which(events[[i]]$type==2)
    t2 <- c(t2,t[events[[i]]$node[tmp],1]-events[[i]]$time[tmp])
    N2 <- c(N2,events[[i]]$n[tmp])
    tmp <- table(events[[i]]$node[tmp])
    tmp <- data.frame(tmp)
    event2mat[i,as.numeric(as.character(tmp$Var1))] <- tmp$Freq
    event2matn[i,as.numeric(as.character(tmp$Var1))] <- sapply(as.numeric(as.character(tmp$Var1)),function (z) sum(events[[i]]$n[which((events[[i]]$node==z)*(events[[i]]$type==2)==1)]))
}

#####repeat over all MCMC runs and plot
#test if known salt tolerant species has higher fitness under salt environment than the other species
tlist <-c(44, 487, 319, 323, 51, 293, 321,  45,  18, 339, 318,  12)
t.test(x=colMeans(amat)[tlist],y=colMeans(amat)[-tlist],alternative="greater",var.equal=F)

#plot histogram of occurence density of adaptation events
library(ape)
aa <- ltt.plot.coords(timetree)
a <- list()
a$time<- -aa[,1]+max(aa[,1])
a$ltt <- aa[,2]
t2new <- sapply(1:length(t2),function (i) rep(t2[i],N2[i]))
t2new <- unlist(t2new)
b <- hist(t2new,breaks=8,plot=F)
c <- numeric(length(b$breaks[-1]))
for (j in 1:length(b$breaks[-1])) {
    i <- b$breaks[-1][j]
    tmp <- min(which(a$time<=i))
    if (tmp!=1) {
        c[j] <- (i-a$time[tmp])*a$ltt[tmp-1]
    }
    c[j] <- c[j]+(a$time[tmp]-a$time[tmp+1])*a$ltt[tmp]
    tmp <- tmp+1
    while (tmp<length(a$ltt)) {
        c[j] <- c[j]+(a$time[tmp]-a$time[tmp+1])*a$ltt[tmp]
        tmp <- tmp+1
    }
}
c <- c(c[1],c[-1]-c[-length(c)])
plot(b$mids,b$counts/c/sum(b$counts/c),type="l")

#plot the amount of niche evolution due to adaptation events on tree
colfunc <- colorRampPalette(c("white","black"))
n <- ceiling(max(colMeans(event2matn)))+2
plot(timetree,type="fan",show.tip.label=F,edge.color=colfunc(n)[round(colMeans(event2matn))+2],edge.width=2,no.margin=T,direction="downwards")
legend_image <- as.raster(matrix(colfunc(n)[n:2], ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.5, y = seq(0,1,l=n), labels = seq(0,n-2,l=n))
rasterImage(legend_image, 0, 0, 1,1)

#plot the amount of niche evolution due to speciation events on tree
timetree.new <- timetree
timetree.new$edge.length[which(timetree.new$edge.length<=0.5)] <- 0.5
colfunc <- colorRampPalette(c("blue","lightgrey"))
colfunc2 <- colorRampPalette(c("lightgrey","red"))
tmp <- hist(colMeans(event1matn))
n <- sum(tmp$mids<0)
idx <- findInterval(colMeans(event1matn)[colMeans(event1matn)<0],tmp$breaks[1:(n+1)],all.inside=T)
n2 <- sum(tmp$mids>=0)
idx2 <-findInterval(colMeans(event1matn)[colMeans(event1matn)>=0],tmp$breaks[(n+1):(n+n2+1)],all.inside=T)
col <- colfunc(n)[idx]
col2 <- colfunc2(n2)[idx2]
colr <- length(length(timetree$edge.length))
colr[colMeans(event1matn)<0] <- col
colr[colMeans(event1matn)>=0] <- col2
plot(timetree.new,type="fan",show.tip.label=F,edge.color=colr,edge.width=2,no.margin=T)
legend_image <- as.raster(matrix(c(colfunc2(n2)[n2:1],colfunc(n)[n:1]), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.5, y = seq(0,1,l=length(tmp$breaks)), labels = tmp$breaks)
rasterImage(legend_image, 0, 0, 1,1)

