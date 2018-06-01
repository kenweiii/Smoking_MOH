#MCMC for the 5x5 transitional matrix
library(dplyr)
library(gtools)
library(truncnorm)
library(Rcpp)
setwd("C:/Users/ephtkw/Desktop/PHMA/Smoking_MOH")
ecigdata <- read.csv("ecig_prevalence.csv") %>% .[15:17,]
sourceCpp("simforward.cpp")
iterations <- 10000

#simulate starting population of 100,000 people with 2014 ratio
#N:1, C:2, Ex:3, EC:4, D:5
counts <- round(ecigdata[1,7:11], 3) * 1000
simpop <- NULL
for(i in 1:length(counts))
{
    simpop <- c(simpop, rep(i, counts[i]))
}

################metropolis hastings algorithm to estimate the transitional probability matrix##################
#transitional matrix structure
#           N       C       Ex      EC      D   FUTURE STATE            parameters
#       N   a       b       0       c       0                           N1, N2, N3
#       C   0       d       e       f       g                           C1, C2, C3, C4
#       Ex  0       h       i       j       0                           Ex1, Ex2, Ex3
#       EC  0       k       l       m       n                           EC1, EC2, EC3, EC4
#       D   0       o       p       q       r                           D1, D2, D3, D4
#CURRENT STATE

#initial transitional matrix
tm <- matrix(0,5,5)

N1 <- round(runif(1, 1, 10000))
N2 <- round(runif(1, 1, 10000))
N3 <- round(runif(1, 1, 10000))
abc <- rdirichlet(1, c(N1, N2, N3))
tm[1,c(1,2,4)] <- abc

C1 <- round(runif(1, 1, 10000))
C2 <- round(runif(1, 1, 10000))
C3 <- round(runif(1, 1, 10000))
C4 <- round(runif(1, 1, 10000))
defg <- rdirichlet(1, c(C1, C2, C3, C4))
tm[2,2:5] <- defg

Ex1 <- round(runif(1, 1, 10000))
Ex2 <- round(runif(1, 1, 10000))
Ex3 <- round(runif(1, 1, 10000))
hij <- rdirichlet(1, c(Ex1, Ex2, Ex3))
tm[3,2:4] <- hij

EC1 <- round(runif(1, 1, 10000))
EC2 <- round(runif(1, 1, 10000))
EC3 <- round(runif(1, 1, 10000))
EC4 <- round(runif(1, 1, 10000))
klmn <- rdirichlet(1, c(EC1, EC2, EC3, EC4))
tm[4,2:5] <- klmn

D1 <- round(runif(1, 1, 10000))
D2 <- round(runif(1, 1, 10000))
D3 <- round(runif(1, 1, 10000))
D4 <- round(runif(1, 1, 10000))
opqr <- rdirichlet(1, c(D1, D2, D3, D4))
tm[5,2:5] <- opqr

#Trace plot of each parameter 18 in total
parameterstore <- matrix(0, nrow = 18, ncol = iterations)
rownames(parameterstore) <- c("N1", "N2", "N3", "C1", "C2", "C3", "C4", "Ex1", "Ex2", "Ex3", "EC1", "EC2", "EC3", "EC4", "D1", "D2", "D3", "D4")
logstore <- rep(0, iterations)

#proposal distribution std dev
propose <- 10

###################################for loop for each iteration of the algorithm##############################################################
logpost <- -Inf
for(i in 1:iterations)
{
    #make proposal distributions
    N1 <- round(rtruncnorm(1, 1, Inf, N1, propose))
    N2 <- round(rtruncnorm(1, 1, Inf, N2, propose))
    N3 <- round(rtruncnorm(1, 1, Inf, N3, propose))
    abc <- rdirichlet(1, c(N1, N2, N3))
    tm[1,c(1,2,4)] <- abc
    
    C1 <- round(rtruncnorm(1, 1, Inf, C1, propose))
    C2 <- round(rtruncnorm(1, 1, Inf, C2, propose))
    C3 <- round(rtruncnorm(1, 1, Inf, C3, propose))
    C4 <- round(rtruncnorm(1, 1, Inf, C4, propose))
    defg <- rdirichlet(1, c(C1, C2, C3, C4))
    tm[2,2:5] <- defg
    
    Ex1 <- round(rtruncnorm(1, 1, Inf, Ex1, propose))
    Ex2 <- round(rtruncnorm(1, 1, Inf, Ex2, propose))
    Ex3 <- round(rtruncnorm(1, 1, Inf, Ex3, propose))
    hij <- rdirichlet(1, c(Ex1, Ex2, Ex3))
    tm[3,2:4] <- hij
    
    EC1 <- round(rtruncnorm(1, 1, Inf, EC1, propose))
    EC2 <- round(rtruncnorm(1, 1, Inf, EC2, propose))
    EC3 <- round(rtruncnorm(1, 1, Inf, EC3, propose))
    EC4 <- round(rtruncnorm(1, 1, Inf, EC4, propose))
    klmn <- rdirichlet(1, c(EC1, EC2, EC3, EC4))
    tm[4,2:5] <- klmn
    
    D1 <- round(rtruncnorm(1, 1, Inf, D1, propose))
    D2 <- round(rtruncnorm(1, 1, Inf, D2, propose))
    D3 <- round(rtruncnorm(1, 1, Inf, D3, propose))
    D4 <- round(rtruncnorm(1, 1, Inf, D4, propose))
    opqr <- rdirichlet(1, c(D1, D2, D3, D4))
    tm[5,2:5] <- opqr
    
    #for loop to simulate through each individual, written in C++
    finalpop <- simforward(tm = tm, population = simpop, years = 1)
    
    #NOTE:
    #the log likelihood should be across all datasets or multiple years
    #use finalpop to calculate the log likelihood
    
    #evaluate log posterior this iteration and save or reject values
    thislogpost <- 0 #INSERT LOG POSTERIOR FUNCTION VALUE HERE
    
    if(log(runif(1)) < (thislogpost - logstore[i-1])) #save values
    {
        parameterstore[1,i] <- N1
        parameterstore[2,i] <- N2
        parameterstore[3,i] <- N3
        parameterstore[4,i] <- C1
        parameterstore[5,i] <- C2
        parameterstore[6,i] <- C3
        parameterstore[7,i] <- C4
        parameterstore[8,i] <- Ex1
        parameterstore[9,i] <- Ex2
        parameterstore[10,i] <- Ex3
        parameterstore[11,i] <- EC1
        parameterstore[12,i] <- EC2
        parameterstore[13,i] <- EC3
        parameterstore[14,i] <- EC4
        parameterstore[15,i] <- D1
        parameterstore[16,i] <- D2
        parameterstore[17,i] <- D3
        parameterstore[18,i] <- D4
        logstore[i] <- thislogpost
    }
    else #reject values, use values from previous iteration
    {
        N1 <- parameterstore[1,i] <- parameterstore[1,i-1]
        N2 <- parameterstore[2,i] <- parameterstore[2,i-1]
        N3 <- parameterstore[3,i] <- parameterstore[3,i-1]
        C1 <- parameterstore[4,i] <- parameterstore[4,i-1]
        C2 <- parameterstore[5,i] <- parameterstore[5,i-1]
        C3 <- parameterstore[6,i] <- parameterstore[6,i-1]
        C4 <- parameterstore[7,i] <- parameterstore[7,i-1]
        Ex1 <- parameterstore[8,i] <- parameterstore[8,i-1]
        Ex2 <- parameterstore[9,i] <- parameterstore[9,i-1]
        Ex3 <- parameterstore[10,i] <- parameterstore[10,i-1]
        EC1 <- parameterstore[11,i] <- parameterstore[11,i-1]
        EC2 <- parameterstore[12,i] <- parameterstore[12,i-1]
        EC3 <- parameterstore[13,i] <- parameterstore[13,i-1]
        EC4 <- parameterstore[14,i] <- parameterstore[14,i-1]
        D1 <- parameterstore[15,i] <- parameterstore[15,i-1]
        D2 <- parameterstore[16,i] <- parameterstore[16,i-1]
        D3 <- parameterstore[17,i] <- parameterstore[17,i-1]
        D4 <- parameterstore[18,i] <- parameterstore[18,i-1]
        logstore[i] <- logstore[i-1]
    }
}

#at this point, check the trace plots of parameterstore and logstore













