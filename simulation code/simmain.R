##### simulate a dataset
# set parameter value
expinv.a <- expinv.b <- 0.03
d <- 1
rootu <- 1
sd <- 2
roota <- log(3)
rootb <- log(2)
rootx <- 1
N <- 10
pois.x <- 1
sp <- 0.1
mu <- 0.03
nspecies <- 100
bias <- NULL
area <- 1000
nBK <- 1000
A1 <- 20
A2 <- 20
sim <- simcode(sp,mu,nspecies,A1,A2,rootu,sd,roota,rootb,rootx,d,N,expinv.a,expinv.b,pois.x,area,nBK)
while(inherits(sim,"try-error")) {
	sim <- simcode(sp,mu,nspecies,A1,A2,rootu,sd,roota,rootb,rootx,d,N,expinv.a,expinv.b,pois.x,area,nBK)
}

##### output file name
eventsfilename <- "events1"
postfilename <- "post1"

##### start analysis
# set initial values
tT1 <- sim$tT1
tT2 <- sim$tT2
t <- sim$t
PA <- sim$PA
PO <- sim$PO
BK <- sim$BK
edge <- sim$tree$edge
edge.length <- sim$tree$edge.length
edge.prob1 <- sim$edge.prob1
edge.prob2 <- sim$edge.prob2
exp.A <- c(10,10)
tuning.A <- c(0.1,0.1)
lnorm.roota.m <- log(3)
lnorm.rootb.m <- log(2)
lnorm.roota <- lnorm.rootb <-  0.5
tuning.roota <- tuning.rootb <- 0.1
norm.rootu.m <- 1
norm.rootu <- 3
tuning.rootu <- 0.1
lnorm.sd.m <- 2
lnorm.sd <- 0.5
tuning.sd <- 0.1
exp.d <- 0.5
tuning.d <- 0.1
tuning.t <- 0.1
exp.a <- exp.b <- 10
tuning.a <- tuning.b <- 0.1
A <- c(0.1/(nspecies-1+sp*tT1),0.1)

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

events <- matrix(NA,0,7)
events <- data.frame(events)
names(events) <- c("node","type","time","a","b","n","xmin")
ngen <- 10^6
k <- c(0,0)

# start MCMC
niche_mcmc(ngen)

#### summarize results
#summarize events and tip fundamental niche in the simulation
tmp <- which(sim$sim$events$type==1)
tmp <- table(sim$sim$events$node[tmp])
tmp <- data.frame(tmp)
event1 <- numeric(dim(edge)[1])
event1[as.numeric(as.character(tmp$Var1))] <- tmp$Freq
tmp <- which(sim$sim$events$type==2)
tmp <- table(sim$sim$events$node[tmp])
tmp <- data.frame(tmp)
event2 <- numeric(dim(edge)[1])
event2[as.numeric(as.character(tmp$Var1))] <- tmp$Freq
expected <- as.table(rbind(event1,event2)) #number of true speciation and adaptation events on each edge
avec <- sim$sim$avec[as.character(1:nspecies)] #fundamental niche of each tip species
#read files
output <- read.table(postfilenamne,header=F,sep=" ")
output <- output[,-1]
names(output) <- c("iter","nll","prior","roota","rootb","rootu","sd","d","expinv.a","expinv.b","A1","A2")
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
#amat records the inferred fundamental niche of each tip species; event1mat and event2mat records the inferred number of speciation event and adaptation event inferered on each edge
amat <- matrix(NA,nsamples,nspecies)
bmat <- matrix(NA,nsamples,nspecies)
event1mat <- matrix(0,nsamples,dim(edge)[1])
event2mat <- matrix(0,nsamples,dim(edge)[1])
for (i in 1:nsamples) {
    roota <- output[i,4]
    rootb <- output[i,5]
    rootu <- output[i,6]
    sd <- output[i,7]
    d <- output[i,8]
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
    tmp <- table(events[[i]]$node[tmp])
    tmp <- data.frame(tmp)
    event1mat[i,as.numeric(as.character(tmp$Var1))] <- tmp$Freq
    tmp <- which(events[[i]]$type==2)
    tmp <- table(events[[i]]$node[tmp])
    tmp <- data.frame(tmp)
    event2mat[i,as.numeric(as.character(tmp$Var1))] <- tmp$Freq
}
#regress the inferred fundamental niche of each tip species against the true niche in the simulation
tmp <- summary(lm(colMeans(amat)~avec))
nichear2 <- tmp$r.squared #R2
nicheacoef <- tmp$coefficients[2] #regression coefficient
#test whether an edge has significantly more occurences of adaptation events or speciation events than expected
Nedge <- dim(edge)[1]
NumberOfShifts1 <- NULL
NumberOfShifts2 <- NULL
for (i in 1:nsamples) {
  NumberOfShifts1 <- c(NumberOfShifts1,sum(events[[i]]$type==1))
  NumberOfShifts2 <- c(NumberOfShifts2,sum(events[[i]]$type==2))
}
expectedNumberOfShifts1 <- mean(NumberOfShifts1)
expectedNumberOfShifts2 <- mean(NumberOfShifts2)
threshold <- 5
Nmax <- 1000
geom_p1 <- 1 / (expectedNumberOfShifts1 + 1)
geom_p2 <- 1 / (expectedNumberOfShifts2 + 1)
prior1 <- dgeom(1:Nmax, geom_p1)
prior2 <- dgeom(1:Nmax, geom_p2)
pp1 <- numeric(Nedge)
pp2 <- numeric(Nedge)
for (i in 1:Nmax) {
	pp1 <- pp1 + prior1[i] * (1-dbinom(0,i,prob=1/Nedge))
	pp2 <- pp2 + prior2[i] * (1-dbinom(0,i,prob=1/Nedge))
}
prior1 <- pp1
prior2 <- pp2
post1 <- colMeans(event1mat)
post2 <- colMeans(event2mat)
oddsratio1 <- post1/prior1
oddsratio2 <- post2/prior2
core1 <- which(oddsratio1>=threshold)
core2 <- which(oddsratio2>=threshold)
#calculate the power and false positive rate to identifiy adaptation events and speciation events
event1corr <- sum(is.element(which(event1>0),core1))/length(which(event1>0)) #power to correctly infer true speciation events
event1wron <- sum(!is.element(core1,which(event1>0)))/length(core1) #false positive rate to infer wrong speciation events
event2corr <- sum(is.element(which(event2>0),core2))/length(which(event2>0)) #power to correctly infer true adaptation events
event2wron <- sum(!is.element(core2,which(event2>0)))/length(core2) #false positive rate to infer wrong adaptation events

