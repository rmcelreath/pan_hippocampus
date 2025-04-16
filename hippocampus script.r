# chimpanzee hippocampus script
# data provided by Maëlig Chauvel as of March 2025 as part of EBC consortium
# goal: estimate contrast between captive and wild reared chimpanzees in relative hippocampus volume
# background: hippocampus functionally associated with memory and navigation - wild chimps learn many trees and navigate complex environments - captive chimps sit on their asses and eat monkey chow

######################################################
# libraries

library(rethinking)
library("readxl") # to read the provided sheet

######################################################
# set up data

d <- read_excel("Hippocampus_volumes_Global_CorrectedValues.xlsx",sheet="ForR")
d <- as.data.frame(d)

dat <- list(
    Y = d[,'hippocampus volume(mm3)'],
    B = d[,'Brain Volume(mm3)'],
    A = d$Age,
    R = ifelse( d$Residence=="captive" , 1 , 2 )
)
# R = 1 indicates captive, R = 2 indicates wild

dat$logY <- log(dat$Y)
dat$logB <- log(dat$B)
dat$hl <- d[,'LEFT hippocampus volume(mm3)']
dat$hr <- d[,'RIGHT hippocampus volume(mm3)']
dat$logHL <- log(dat$hl)
dat$logHR <- log(dat$hr)

# anterior/posterior right and left
dat$logAHR <- log(d[,'aHPC_Right_volume(mm3)'])
dat$logAHL <- log(d[,'aHPC_Left_volume(mm3)'])

dat$logPHR <- log(d[,'pHPC_Right_volume(mm3)'])
dat$logPHL <- log(d[,'pHPC_Left_volume(mm3)'])

######################################################
# visualize sample

# pairs plot for Y B A, color by R
pairs( dat$Y ~ dat$B + dat$A , col=dat$R , cex=2 , lwd=2 , labels=c("Hippocampus mm3","Brain mm3","Age") )

plot( dat$A , dat$Y/dat$B , col=dat$R , cex=2 , lwd=2 , xlab="Age" , ylab="Hippocampus/Brain" )
mtext("red=wild, black=captive")

# pairs showing left/right
d2 <- data.frame(
    HRL=c(dat$hr,dat$hl),
    id=rep(1:length(dat$hl),times=2),
    side=rep(c("R","L"),each=length(dat$hl)),
    B=rep(dat$B,times=2),
    A=rep(dat$A,times=2),
    R=rep(dat$R,times=2)
)

pairs( d2$HRL ~ d2$B + d2$A , col=d2$R , pch=ifelse(d2$side=="R",1,16) , cex=2 , lwd=2 , labels=c("Hippocampus mm3 (left & right)","Brain mm3","Age") )

plot( d2$A , d2$HRL/d2$B , col=d2$R , cex=2 , xlab="Age" , ylab="Hippocampus/Brain" )

######################################################
# baseline model using average of left and right

m0 <- ulam(
    alist(
        logY ~ normal(mu,sigma),
        mu <- a[R] + bB*logB + bA*A,
        a[R] ~ normal(0,10),
        bB ~ normal(1,1),
        bA ~ normal(1,1),
        sigma ~ exponential(1)
    ), data=dat , chains=8 , cores=8 )

precis(m0,2)

post <- extract.samples(m0)
diff <- post$a[,2]-post$a[,1]
quantile(diff)
dens(post$a[,2]-post$a[,1],show.zero=TRUE)

# mass below zero
sum(diff<0)/length(diff)

######################################################
# baseline model with left/right instead of average

# without covariance, just as scaffold to develop the covariance version
m1_dev <- ulam(
    alist(
        logHL ~ normal(mu,sigma),
        logHR ~ normal(mu,sigma),
        mu <- a[R] + bB*logB + bA*A,
        a[R] ~ normal(-5,2),
        bB ~ normal(1,1),
        bA ~ normal(1,1),
        sigma ~ exponential(1)
    ), data=dat , chains=4 , cores=4 )

# covariance version
# priors here from prior predictive logic in doc
m1 <- ulam(
    alist(
        c(logHL,logHR) ~ multi_normal(c(muL,muR),Rho,Sigma),
        muR <- muL + delta,
        muL <- a[R] + bB*logB + g*A,
        a[R] ~ normal(-5.7,2),
        delta ~ normal(0,0.5),
        bB ~ normal(1,0.5),
        g ~ normal(0,0.5),
        Rho ~ lkj_corr(4),
        vector[2]:Sigma <<- rep_vector(sigma,2),
        sigma ~ exponential(1)
    ), data=dat , chains=8 , cores=8 , refresh=1000 )

precis(m1,3)

post <- extract.samples(m1)
diff <- post$a[,2]-post$a[,1]
quantile(diff)
blank2()
dens(post$a[,2]-post$a[,1],show.zero=TRUE,xlab="contrast [wild-captive]")

# mass below zero
mbz <- sum(diff<0)/length(diff)
text(0.15,2,round(1-mbz,3))

# prior predictive for covariance model

N <- 50
Hcaptive <- matrix(NA,N,2)
Hwild <- matrix(NA,N,2)
logB <- rep(NA,N)
A <- rep(NA,N)
for ( i in 1:N ) {
    a <- rnorm(2,-5.7,2)
    delta <- rnorm(1,0,0.1)
    bB <- rnorm(1,1,0.25)
    Rho <- rlkjcorr(1,K=2,eta=4)
    sigma <- rnorm(1,0.1,0.1)
    idx <- sample(1:16,size=1)
    logB[i] <- dat$logB[idx]
    A[i] <- dat$A[idx]
    mucaptive <- a[1] + bB*logB[i]
    Hcaptive[i,] <- rmvnorm2(1,Mu=c(mucaptive,mucaptive+delta),Rho=Rho,sigma=rep(sigma,2))
    muwild <- a[2] + bB*logB[i]
    Hwild[i,] <- rmvnorm2(1,Mu=c(muwild,muwild+delta),Rho=Rho,sigma=rep(sigma,2))
}

# blank2()
plot( (logB) , (Hcaptive[,1]) , xlab="log Brain" , ylab="log Hippocampus" )
points( (logB) , (Hwild[,1]) , col=2 )

######################################################
# add anterior/posterior volumes

m2 <- ulam(
    alist(
        c(logHL,logHR,logAHL,logAHR,logPHL,logPHR) ~ multi_normal(c(muL,muR,muAL,muAR,muPL,muPR),Rho,Sigma),
        muR <- muL + delta,
        muL <- a[R] + bB*logB,
        muAR <- muAL + delta,
        muAL <- aA[R] + bB*logB,
        muPR <- muPL + delta,
        muPL <- aP[R] + bB*logB,
        a[R] ~ normal(-5.7,1),
        aA[R] ~ normal(-6,1),
        aP[R] ~ normal(-6,1),
        delta ~ normal(0,0.5),
        bB ~ normal(1,0.5),
        Rho ~ lkj_corr(4),
        Sigma ~ exponential(1)
    ), data=dat , chains=8 , cores=8 , refresh=1000 )

precis(m2,3)

######################################################
# model hippocampus as age-residence-specific proportion of brain B
# Y = p[A,R]*B   =>   log(Y) = log(B) + log(p[A,R])
# so like log-normal with age-residence-specific intercepts, but constant exposure with log(B)

# covariance version
m3 <- ulam(
    alist(
        # hippocampus model
        c(logHL,logHR) ~ multi_normal(c(muL,muR),Rho,Sigma),
        muR <- muL + delta,
        muL <- a[R] + p*A + logB,
        p ~ normal(0,0.5),
        delta ~ normal(0,0.5),
        # brain size model
        logB ~ normal(nu,tau),
        nu <- b + log(1-exp(-phi*A)),
        # priors
        a[R] ~ normal(-5,5),
        Rho ~ lkj_corr(4),
        vector[2]:Sigma <<- rep_vector(sigma,2),
        sigma ~ exponential(1),
        b ~ normal(10,10),
        phi ~ exponential(10),
        tau ~ exponential(1)
    ), data=dat , chains=8 , cores=8 , sample=TRUE )

precis(m3,3)

# plot brain growth against age
rep_vector <- function(x,n) rep(x,times=n)
pred_dat <- data.frame(A=1:50,R=1,logB=10)
nu_pred <- link(m2,data=pred_dat)$nu
plot( dat$A , dat$B , xlab="age" , ylab="brain volume" , col=dat$R )
lines( apply(exp(nu_pred),2,mean) , lwd=2 )
ci <- apply(exp(nu_pred),2,PI)
lines( 1:50 , ci[1,] )
lines( 1:50 , ci[2,] )
