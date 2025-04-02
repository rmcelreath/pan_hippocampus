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
dens(post$a[,2]-post$a[,1])

# mass below zero
sum(diff<0)/length(diff)

######################################################
# baseline model with left/right instead of average

m1 <- ulam(
    alist(
        logHL ~ normal(mu,sigma),
        logHR ~ normal(mu,sigma),
        mu <- a[R] + bB*logB + bA*A,
        a[R] ~ normal(0,5),
        bB ~ normal(1,1),
        bA ~ normal(1,1),
        sigma ~ exponential(1)
    ), data=dat , chains=4 , cores=4 )

# covariance version
m1 <- ulam(
    alist(
        c(logHL,logHR) ~ multi_normal(c(muL,muR),Rho,Sigma),
        muL <- a[R] + bB*logB,
        muR <- a[R] + bB*logB + delta,
        a[R] ~ normal(0,2),
        delta ~ normal(0,0.5),
        bB ~ normal(1,0.5),
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

######################################################
# model hippocampus as age-residence-specific proportion of brain B
# Y = p[A,R]*B   =>   log(Y) = log(B) + log(p[A,R])
# so like log-normal with age-residence-specific intercepts, but constant exposure with log(B)

# covariance version
m2 <- ulam(
    alist(
        # hippocampus model
        c(logHL,logHR) ~ multi_normal(c(muL,muR),Rho,Sigma),
        muL <- a[R] + p*A + logB,
        muR <- a[R] + p*A + logB + delta,
        p ~ normal(0,0.5),
        delta ~ normal(0,0.5),
        # brain size model
        logB ~ normal(nu,tau),
        nu <- b + log(1-exp(-phi*A)),
        # priors
        a[R] ~ normal(0,10),
        Rho ~ lkj_corr(4),
        vector[2]:Sigma <<- rep_vector(sigma,2),
        sigma ~ exponential(1),
        b ~ normal(10,10),
        phi ~ exponential(10),
        tau ~ exponential(1)
    ), data=dat , chains=8 , cores=8 , sample=TRUE )

precis(m2,3)

# plot brain growth against age
rep_vector <- function(x,n) rep(x,times=n)
pred_dat <- data.frame(A=1:50,R=1,logB=10)
nu_pred <- link(m2,data=pred_dat)$nu
plot( dat$A , dat$B , xlab="age" , ylab="brain volume" , col=dat$R )
lines( apply(exp(nu_pred),2,mean) , lwd=2 )
ci <- apply(exp(nu_pred),2,PI)
lines( 1:50 , ci[1,] )
lines( 1:50 , ci[2,] )
