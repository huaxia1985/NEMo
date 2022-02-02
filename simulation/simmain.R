##### section 1: simulate a dataset
iii <- 1
set.seed(iii)

r <- 0.03
d <- 1
roota <- -1
rootb <- -1
rootu <- 1
sd <- 2
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
sim <- simcode(sp,mu,nspecies,A1,A2,rootu,sd,roota,rootb,rootx,d,N,r,pois.x,area,nBK)
while(inherits(sim,"try-error")) {
	sim <- simcode(sp,mu,nspecies,A1,A2,rootu,sd,roota,rootb,rootx,d,N,r,pois.x,area,nBK)
}
integrand <- function (tt) {mu*(1-exp(-(sp-mu)*tt))/(sp-mu*exp(-(sp-mu)*tt))}

##### section 2: setting up the output files for results from NEMo

#two chains per simulation
for (aaa in 1:2) {
eventsfilename <- paste0(getwd(),"/simulation/events",iii,"_",aaa)
postfilename <- paste0(getwd(),"/simulation/post",iii,"_",aaa)
filename <- paste0(getwd(),"/simulation/work",iii,"_",aaa)

##### section 3: setting up the initial values NEMo

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
tuning.t <- min(edge.length)
A <- c(0.1/(nspecies-1+sp*tT1),0.1)

norm.roota.m <- norm.rootb.m <- -1
norm.rootu.m <- 1
sd <- lnorm.sd.m <- 2
norm.roota.sd <- norm.rootb.sd <- 1
norm.rootu.sd <- 3
lnorm.sd.sd <- 0.5
exp.d <- 0.5
exp.r <- 10
exp.A <- c(10,10)

tuning.roota <- tuning.rootb <- 0.1
tuning.rootu <- 1
tuning.sd <- 0.1
tuning.r <- 0.1
tuning.d <- 0.1
tuning.A <- c(0.1,0.1)

events <- matrix(NA,0,7)
events <- data.frame(events)
names(events) <- c("node","type","time","a","b","n","xmin")

#calculate initial likelihood
rootniche <- exp(-((0:(N-1))/(exp(roota)+1e-6))^(exp(rootb)+1e-6))
rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
rootback <- rootback/sum(rootback)
rootdistro <- c(rep(1,N),rootniche)*rootback
rootdistro <- rootdistro/sum(rootdistro)
rootdistr <- rootdistro[(N+1):(2*N)]
rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
lik.old <- likcal(node=0,root.included=F,events=events,edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,nodeniche=rootniche,nodeback=rootback,d=d,bias=bias,area=area)

#number of MCMC generations to run
ngen <- 3*10^6
k <- c(0,0)

##### section 4: running the NEMo analysis on simulate data

niche_mcmc(ngen)

##### section 5: summarizing results from NEMo

nichear2 <- matrix(NA,20,2)
nicheur2 <- matrix(NA,20,2)
event1corr <- matrix(NA,20,2)
event1wron <- matrix(NA,20,2)
event2corr <- matrix(NA,20,2)
event2wron <- matrix(NA,20,2)

for (iii in 1:20) {
for (aaa in 1:2) {
filename <- paste0(getwd(),"/simulation/work",iii,"_",aaa)
load(filename)
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
avec <- sim$sim$avec[as.character(1:nspecies)] #true fundamental niche of each tip species
uvec <- sim$sim$uvec[as.character(1:nspecies)] #true available niche of each tip species
#read in results
eventsfilename <- paste0("events",iii,"_",aaa)
postfilename <- paste0("post",iii,"_",aaa)
output <- read.table(postfilename,header=F,sep=" ")
output <- output[,-1]
names(output) <- c("iter","nll","prior","roota","rootb","rootu","sd","d","r","A1","A2")
if (aaa==1) {
plot(c(1:dim(output)[1]),output[1:dim(output)[1],2]-output[1:dim(output)[1],3],type="l")
} else {
points(c(1:dim(output)[1]),output[1:dim(output)[1],2]-output[1:dim(output)[1],3],type="l")
}
burnin <- 600 #picked by the plot
output <- output[-c(1:burnin),]
tmp <- read.table(eventsfilename,header=F,sep=" ",fill=T,stringsAsFactors=F)
idx <- which(is.na(tmp[,8]))
nsamples <- 3000-burnin
events <- vector("list",nsamples)
for (i in (burnin+1):3000) {
    tmp2 <- tmp[(idx[i]+1):(ifelse(i<length(idx),idx[i+1],dim(tmp)[1]+1)-1),-1]
    names(tmp2) <- c("node","type","time","a","b","n","xmin")
    events[[i-burnin]] <- data.frame(apply(tmp2,2,function (x) as.numeric(x)))
}
#extracting posterior samples
amat <- matrix(NA,nsamples,nspecies) #alpha for fundamental niche
umat <- matrix(NA,nsamples,nspecies) #mu for available niche
event1mat <- matrix(0,nsamples,dim(edge)[1]) #number of speciation events per branch
event2mat <- matrix(0,nsamples,dim(edge)[1]) #number of adaptation events per branch
for (i in 1:nsamples) {
    roota <- output[i,4]
    rootb <- output[i,5]
    rootu <- output[i,6]
    sd <- output[i,7]
    d <- output[i,8]
    rootniche <- exp(-((0:(N-1))/(exp(roota)+1e-6))^(exp(rootb)+1e-6))
    rootback <- pnorm(seq(-(N-1),N,1),mean=rootu,sd=sd)-pnorm(seq(-N,N-1,1),mean=rootu,sd=sd)
    rootback <- rootback/sum(rootback)
    rootdistro <- c(rep(1,N),rootniche)*rootback
    rootdistro <- rootdistro/sum(rootdistro)
    rootdistr <- rootdistro[(N+1):(2*N)]
    rootdistr[1] <- rootdistr[1]+sum(rootdistro[1:N])
    lik <- likcal(node=0,root.included=F,events=events[[i]],edge=edge,edge.length=edge.length,nspecies=nspecies,PO=PO,PA=PA,N=N,nodeu=rootu,sd=sd,nodea=roota,nodeb=rootb,nodex=rootx,nodedistr=rootdistr,nodedistro=rootdistro,nodeniche=rootniche,nodeback=rootback,d=d,bias=bias,area=area)
    amat[i,] <- lik$avec[as.character(1:nspecies)]
    umat[i,] <- lik$uvec[as.character(1:nspecies)]
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
nichear2[iii,aaa] <- tmp$r.squared
#regress the inferred available niche of each tip species against the true niche in the simulation
tmp <- summary(lm(colMeans(umat)~uvec))
nicheur2[iii,aaa] <- tmp$r.squared
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
event1corr[iii,aaa] <- sum(is.element(which(event1>0),core1))/length(which(event1>0)) #power to correctly infer true speciation events
event1wron[iii,aaa] <- sum(!is.element(core1,which(event1>0)))/length(core1) #false positive rate to infer wrong speciation events
event2corr[iii,aaa] <- sum(is.element(which(event2>0),core2))/length(which(event2>0)) #power to correctly infer true adaptation events
event2wron[iii,aaa] <- sum(!is.element(core2,which(event2>0)))/length(core2) #false positive rate to infer wrong adaptation events

#generate figure 3
par(mfrow=c(2,6))
boxplot(as.numeric(nichear2),ylim=c(0,1))
abline(h=0.8)
boxplot(as.numeric(nicheur2),ylim=c(0,1))
abline(h=0.8)
boxplot(as.numeric(event1corr),ylim=c(0,1))
abline(h=0.8)
boxplot(as.numeric(event1wron),ylim=c(0,1))
abline(h=0.2)
boxplot(as.numeric(event2corr),ylim=c(0,1))
abline(h=0.8)
boxplot(as.numeric(event2wron),ylim=c(0,1))
abline(h=0.2)

##### section 6: performing ENM on simulated data

library(splines)
#install.packages("multispeciesPP_1.0.tar.gz"). Downloaded from https://github.com/wfithian/multispeciesPP
library(multispeciesPP)

for (iii in 1:20) {
    filename <- paste0(getwd(),"/simulation/work",iii,"_",1)
    load(filename)
    background <- unlist(sim$BK)
    (x.knots <- quantile(background,(1:5)/6,na.rm=TRUE))
    x.basis <- ns(background,knots=x.knots)
    sdm.formula <- formula(~predict(x.basis,newx=x))
    bias.formula <- formula(~y)
    PA <- cbind(sim$PA[[1]]$x,sim$PA[[1]]$y)
    for (sss in 2:100) {
        PA <- cbind(PA,sim$PA[[sss]]$y)
    }
    colnames(PA) <- c("x",sim$tree$tip.label)
    PA <- as.data.frame(PA)
    PA$isPO <- 0
    PA$y <- rnorm(dim(PA)[1])
    PO <- vector("list",100)
    for (sss in 1:100) {
        PO[[sss]] <- data.frame(x=sim$PO[[sss]],y=rnorm(length(sim$PO[[sss]])))
    }
    names(PO) <- sim$tree$tip.label
    mod <-
        multispeciesPP(sdm.formula,bias.formula,
        PA=PA,PO=PO,BG=data.frame(x=background,y=rnorm(length(background))),species=sim$tree$tip.label,
                       control=list(trace=TRUE,maxit=100,epsilon=1e-6))
    ry=model.matrix(~ns(c(1:10),knots=x.knots))
    species.niche <- matrix(NA,100,10)
    for (i in 1:100) {
      prd <- exp(mod$species.coef[1:7,i]%*%t(ry))
      prd <- prd/sum(prd)
      species.niche[i,] <- prd
    }
    filename <- paste0(getwd(),"/simulation/niche",iii)
    save.image(filename)
}

##### section 7: plotting results on phylogeny

#the example here uses the last posterior sample of chain 1 of each simulated data
library(phytools)
par(mfrow=c(4,5))
KS_nemo <- matrix(NA,20,100)
KS_enm <- matrix(NA,20,100)
for (iii in 1:20) {
filename <- paste0(getwd(),"/simulation/work",iii,"_",1)
load(filename)
speciesniche <- sim$sim$dendistr[as.character(c(1:100)),]
speciesniche <- cbind(speciesniche,NA,lik.old$dendistr[as.character(c(1:100)),])
sdmfilename <- paste0(getwd(),"/simulation/niche",iii)
load(sdmfilename)
speciesniche <- cbind(speciesniche,NA,species.niche)
PO.size <- sapply(sim$sim$PO,function (i) length(i))
KS_nemo[iii,] <- sapply(c(1:100), function (i) max(abs(cumsum(speciesniche[i,1:10])-cumsum(speciesniche[i,12:21]))))
KS_enm[iii,] <- sapply(1:100, function (i) max(abs(cumsum(speciesniche[i,1:10])-cumsum(speciesniche[i,23:32]))))
#plot(PO.size,KS_enm[iii,],pch=16)
#points(PO.size,KS_nemo[iii,],col="red",pch=16)
PO <- sapply(sim$sim$PO,function (i) table(i)[as.character(c(1:10))])
PO[is.na(PO)] <- 0
PO <- apply(PO,2,function (i) i/sum(i))
speciesniche <- cbind(speciesniche,NA,t(PO))
speciesniche <- cbind(speciesniche,NA,1,sim$sim$nichedistr[as.character(c(1:100)),-1]/rowSums(sim$sim$nichedistr[as.character(c(1:100)),-1]),NA,1,lik.old$nichedistr[as.character(c(1:100)),-1]/rowSums(lik.old$nichedistr[as.character(c(1:100)),-1]))
speciesniche <- cbind(speciesniche,NA,rowSums(sim$sim$backdistr[as.character(c(1:100)),1:11]),sim$sim$backdistr[as.character(c(1:100)),12:20],NA,rowSums(lik.old$backdistr[as.character(c(1:100)),1:11]),lik.old$backdistr[as.character(c(1:100)),12:20])
rownames(speciesniche) <- sim$tree$tip.label

#generate figure 4
#plot tree with tip niches
figfilename <- paste0(getwd(),"/simulation/fig",iii,"_",aaa,".pdf")
pdf(figfilename)
phylo.heatmap(tree=sim$tree,X=speciesniche,labels=F,tip.labels=F)
#plot events on tree
u.edge <- sim$sim$events[sim$sim$events$type==2,"node"]
u.time <- sim$sim$events[sim$sim$events$type==2,"time"]
u.n <- sim$sim$events[sim$sim$events$type==2,"n"]
u.edge <- u.edge[u.n!=0]
u.time <- u.time[u.n!=0]
u.n <- u.n[u.n!=0]
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
yy <- lastPP$yy[lastPP$edge[u.edge,2]]
tt <- max(lastPP$xx)-min(lastPP$xx)
points(xx-tt*(sim$tree$edge.length[u.edge]-u.time)/max(t),yy,pch=24,bg=colr,cex=2)

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

xx <- lastPP$xx[lastPP$edge[u.edge,2]]
yy <- lastPP$yy[lastPP$edge[u.edge,2]]
points(xx-tt*(sim$tree$edge.length[u.edge]-u.time)/max(t),yy,pch=24,bg=colr,col="white",cex=1.5)
dev.off()
}
