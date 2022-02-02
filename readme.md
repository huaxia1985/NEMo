# NEMo: Niche Evolution Model

This contains the source code of NEMo introduced in paper: 
"NEMo: a Niche Evolution Model for the simultaneous inference of contemporary niches and the history of niche evolution"
By Xia Hua, Marcel Cardillo, and Lindell Bromham

Please cite the paper if you use the code

For example of using NEMo on empirical data, please see
"Adapting to extremes: reconstructing evolution in response to changing climate over time and space in the diverse Australian plant genus Acacia"
By Xia Hua, Marcel Cardillo, and Lindell Bromham


The code implements our proposed Bayesian approach to infer the environmental niches of tip species and their history of niche evolution along the phylogeny.

This repo gives information on:
1. replicating the results in the papers
2. using the method on your own data

## Replicating the results in the papers

#### Step 1

Download this github repository to your computer. From here on in, we'll assume that this is now located at `~/Desktop/NEMo/`

#### Step 2

Source all the functions to R and load the related libraries

```r
  file.sources = list.files(path = "~/Desktop/NEMo/", pattern="*.R", full.names = T)
  sapply(file.sources,source,.GlobalEnv)
```

#### Step 3

Replicating the simulation study

The main code is in `~/Desktop/NEMo/simulation/simmain`

The code leads you through each step to simulate a phylogeny and a set of species presence/absence locations, as well as run the analysis on the simulatd data.
To use the code, simply run each line of the code.

#### Section 1 sets up parameter value to simulate a set of data, where

r = heritability of individual tolerance range;

d = dispersal rate;

rootu, sd = the mean and the standard deviation of the root available niche, which is assumed to be a normal distribution;

roota, rootb = parameter alpha and beta of the root fundamental niche;

rootx = location of dispersal barrier in the root available niche;

N = the number of bin of the environmental axis;

pois.x = the mean of a poisson distribution that simulates parameter x in a dispersal event;

sp, mu = speciation and extinction rate to simulate the phylogeny;

nspecies = the number of tip species to simulate;

bias = NULL, do not simulate sampling bias;

area = parameter D, the area of available niche used in the ENM part of the method;

nBK = the number of background locations to sample;

A1, A2 = the number of gradients speciation events and in-situ adaptation events to simulate on the phylogeny.

Then call function simcode to simulate the data.

#### Section 2 sets up the output files for the analysis, where

1) eventsfilename = the file to write the event matrix that records each inferred events on the phylogeny in every 1000 steps in a MCMC, for example:
```r
   node type      time          a          b   n xmin
1     2    1 5.8020000 0.04359012 0.04359012 0.3    2
2     3    2 4.1170970 0.04359012 0.04359012 6.0    3
```
where each row is an event,

'node' = the branch where the event occurred; 

'type' = the type of the events, 1 being a gradient speciation event, 2 being an in-situ adaptation event;

'time' = the time interval from the begnning of the branch to the occurence of the event;

'a' and 'b' = the heritability of parameter alpha and beta of the fundamental niche, these are assumed the same for now;

'n' = the parameter n or p of the event;

'xmin' = the location of barrier in the dispersal event following the event.

2) postfilename = the file to write -lnL, prior, and inferred global parameter values in every 1000 steps in a MCMC, for example:
```r
iter     nll     prior    roota     rootb    rootu       sd        d
1 1000 4962319 -28.14734 1.064112 0.6291913 2.186328 3.539255 1.765835
           r          A1       A2      bias1      bias2
1 0.04359012 7.00663e-05 2.356269 -0.1545956 -0.1095036
```
where,'iter' = the step of the MCMC;

'nll' = the negative log likelihood of the model;

'prior' = the log prior probability of the model;

'roota' and 'rootb' = the parameter alpha and beta of the root fundamental niche;

'rootu' and 'sd' = the mean and the standard deviation of the root avaiable niche;

'd' = dispersal rate;

'r' = heritability;

'A1' = occurence rate of gradient speciation events per time step on the phylogeny;

'A2' = overall occurence rate of in-situ adaptation events on the phylogeny;

'bias1' and 'bias2' = the coefficient for the included factors that may influence sampling efforts (two factors are included here).

#### Section 3 sets up the initial values for the MCMC analysis, where

edge.prob1, edge.prob2 and tuning.t = the prior and the tuning of the time of a gradient speciation event and an in-situ adaptation event;

exp.A and tuning.A = the exponential prior and the tuning of the overall occurence rate of gradient speciation events and an in-situ adaptation events on the phylogeny;

norm.roota.m and norm.rootb.m, norm.roota and norm.rootb, tuning.roota and tuning.rootb = the mean and sd of the normal prior and the tuning of parameter alpha and beta on log scale of the root fundamental niche;

norm.rootu.m, norm.rootu, and tuning.rootu = the mean and sd of the normal prior and the tuning of the mean of the root available niche;

lnorm.sd.m, lnorm.sd, and tuning.sd  = the mean and sd of the lognormal prior and the tuning of the sd of the root available niche;

exp.d and tuning.d = the exponential prior and the tuning of the dispersal rate;

exp.r and tuning.r = the exponential prior and the tuning of the heritability;

events = and empty matrix to initialize the event matrix;

k = a vector of c(0,0) to initialize the number of inferred gradient speciation events and in-situ adaptation events.

Then call likcal to calculate the initial likelihood of the model

#### Section 4 calls niche_mcmc to run the actual MCMC.

You may want to repeat all the above steps to replicate simulations and analyze each simulated data.

#### Section 5 summarizes the results stored in the output files. 

This generates figure 3 in the paper.
1. nichear2 = the estimated R2 in regressing parameter a and b inferred for each tip fundamental niche to the simulated niches.
2. nicheur2 = the estimated R2 in regressing parameter u inferred for each tip available niche to the simulated niches.
2. event1corr and event1wrong = the true positive rate and the false positive rate in inferring gradient speciatin events
3. event2corr and event2wrong = the true positive rate and the false positive rate in inferring in-situ adaptation events
4. posterior samples of global parameters are in output, where roota, rootb, rootu, sd are parameters for the root niches

#### Section 6 performs ENM_only analysis on simulated data

This requires install the "multispeciesPP" package for Poisson point process with correction for bias.
The results are stored in file called "niche".

#### Section 7 plots results on phylogeny

This generates figure 4 in the paper.
The simulated niches at each node and tip are stored in the results of the simdistr function.
The inferred niches at each node and tip are stored in the results of the likcal function.
1. backdistr = available niche
2. nichedistr = fundamental niche
3. dendistr = occupied niche

#### Step 4

Replicating the case study

The data is stored in `~/Desktop/NEMo/case study/`.

'arid' for Acacia drought tolerance and 'salt' for Acacia salt tolerance

Each file includes:
1. area = parameter D of each tip species used in the ENM part of the method
2. BK = a set of background locations and the related variables extracted from each location
3. N = the number of bins of the environmental axis
4. PA = a set of presence-absence locations and the related variables extracted from each location
5. PO = a set of presence-only locations and the related variables extracted from each location
6. timetree = dated phylogeny of Acacia

To use the data, load the file into R, using drought tolerance for example,
```r
load("~/Desktop/NEMo/case study/arid")
```

The main code is in `~/Desktop/NEMo/case study/main`.

The code leads you through each step to run one MCMC on Acacia drought tolerance and Acacia salt tolerance.
To use the code, simply run each line of the code.

Same as the simulation code, the analysis code starts with setting initial parameter values, priors, and the output file names.
Then call niche_mcmc to run one MCMC. 

You may want to repeat call niche_mcmc multiple times to run multiple MCMCs.

The last section summarize results. Same as the simulation code, the summary takes several steps:

1. reading in posterior samples from eventsfile and postfile 
2. using likcal function to calculate niches at each node and tip, given each posterior sample
3. Plotting changes in parameter alpha, events of niche evolution, and estimated niches on the phylogeny. This generates Figure S1 in the Acacia case study paper.
4. Calculating the confidence intervals of the occurence density of events of niche evolution within diffrent time periods. This generates Figure 3 in the Acacia case study paper.

## Using the method on your own data

#### Step 1

Generating your data in the same format as the data in '~/Desktop/NEMo/case study/arid'

The data includes:
1) a dated phylogeny in a "phylo" class
2) presence-only data: a list, where each element corresponds to each tip species, in the same order of the tip.label in the phylogeny. 
In each element, there is a vector named 'x' that records the environmental condition (the bin index on the environmental axis) of each presence-only location of the species; and a matrix, with each row corresponding to a location and recording the value of included factors at that location that may affect sampling effort.
3) background data: same as the presence-only data, except that each location is randomly sampled around the presence locations. These locations are sampled to approximate the integral term in the ENM model.
4) presence-ansence data: similar to the presence-only data, this is a list with each element corresponding to each tip species. In each element, there is the vector 'x' that records the environmental condition of each presence-absence location; and a vector named 'y' that records whether the species is present (1) or absent (0) in each presence-absence location.
5) area: the total area around each presence locations where background data is sampled.

To prepare the species distribution data, you can use the following code for the type of environmental variable you are interested in:
1) load a fine resolution map that records the environmental conditions of your study area, for example,
```r
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(ENMTools) 
salt <- raster("where you save the map")
proj4string(salt) <- CRS(the projection of the map)
```

2) read in species location coordinates, for example if you save the coordinates of each presence-only location in a csv file, with column name "SpeciesName","Longitude" and "Latitude":
```r
sp.pts <- read.csv("where you save the coordinates")
species <- sp.pts$SpeciesName
coordinates(sp.pts) <- c("Longitude", "Latitude")
proj4string(sp.pts) <- CRS(the projection of the map)
```

3) extract the environmental condition of each presence-only location and construct the presence-only data:
```r
po <- numeric()
po$long <- sp.pts$Longitude
po$lat <- sp.pts$Latitude
po$salt <- extract(salt,sp.pts)
po <- as.data.frame(po)
PO <- lapply(as.list(timetree$tip.label),function (sp) po[which(species==sp),])
names(PO) <- timetree$tip.label
```

3) sample background locations from each species presence locations and calculate area, for example, if you want to sample the same number of background locations as the presence-only locations in the 20km radius around each presence locations:
```r
mask <- raster("where you save a mask map of the study area")
nspecies <- length(length(timetree$tip.label)
bk <- vector("list",nspecies)
area <- numeric(nspecies)
for (i in 1:nspecies) {
	bk[[i]] <- background.points(points=cbind(PO[[i]]$long,PO[[i]]$lat),radius=20000,n=length(PO[[i]]$long),mask=mask)
	area[i] <- background.area(points=cbind(PO[[i]]$long,PO[[i]]$lat),radius=20000)
}
```

4) extract the environmental condition of each background location:
```r
BK <- vector("list",nspecies)
for (i in 1:nspecies) {
    BK[[i]]$salt <- extract(salt,bg[[i]])
    names(BK) <- timetree$tip.label
}
```

5) use the distribution of the environmental condition of background locations to determine the number of bins to discretise the enironmental axis:
```r
tmp <- hist(log(unlist(BK),plot=F)
N <- length(tmp$mids)
x.axis <- tmp$breaks
```

6) translate the environmental conditions in PO and BK to discrete units:
```r
for (i in 1:nspecies) {
    PO[[i]]$x <- findInterval(PO[[i]]$salt,x.axis,all.inside=T)
    BK[[i]]$x <- findInterval(BK[[i]]$salt,x.axis,all.inside=T)
}
```

#### Step 2

Analysing your data using the code in '~/Desktop/NEMo/case study/main'.

You may want to modify the initial parameter values and the prior hyperparameters, and run niche_mcmc multiple times for multiple MCMCs.

Then use the code to summarize the results.
