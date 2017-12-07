spin = scan("input_R.txt")

library(cutoff)
hist(spin,10,F,xlab="Spin",ylab="Density", main=NULL,col="light blue")


# A kernel density estimation of the distribution:
lines(density(spin),lwd=1.5,col="blue")


# Estimating the parameters of the finite mixture model:
(spin_out <- em(spin,"normal","normal"))

# The confidence interval of the parameter estimates (nb is the number of MonteCarlo simulations):
confint(spin_out,nb=300,level=.95) 

# The plot:
hist(spin,10,F,xlab="Spin Frequency",ylab="Density",xlim=c(100,700),ylim=c(0,0.006),main=NULL,col="light blue")
lines(spin_out,lwd=1.5,col="black")

# Estimating a cutoff value from this fitted finite mixture model:
(cut_off <- cutoff(spin_out))

# Plotting it:
polygon(c(cut_off[-1],rev(cut_off[-1])),c(0,0,.55,.55),col=rgb(0.2,0.5,0.5,.5),border=NA)

#abline(v=cut_off[-1],lty=1,col="red", lwd=2)
abline(v=cut_off[1],col="red", lwd=3)


#### Plot the KDE ########
density(spin, kernel=c("gaussian"))
plot(density(spin, kernel=c("rectangular")), xlab="Spin Frequency", lwd=3, cex.lab=1.4, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
plot(density(spin, kernel=c("gaussian")), xlab="Spin Frequency", lwd=3, cex.lab=1.4, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)



######## Show qq-plot and theoretical line
# Here you can show that spin is non-normal and spiders is normal
# Use also a shapiro.wilk test to show the non-normal vs. normal data
qqnorm(spin)
qqline(spin)
qqnorm(spiders)
qqline(spiders)
qqnorm(spin, lwd=2, cex.lab=1.4, cex.axis=1.3, cex.main=1.3, cex.sub=1.3, ylab = 'Spin Frequency (Hz)')
qqline(spin, lwd=3)



############ Test 1: Effect Size

## 1. Effect Size RMSPs vs. AMXPs

library(compute.es)
mes(mean(amxps),mean(rmsps), sd(amxps), sd(rmsps), 19, 399, level=95)
### Effect size d = 0.74
### Then we calculate the power (assuming both distributions are normal, so we use a t-test. If they are not normal then also the MC simulations of Papitto are wrong).
library(pwr)
pwr.t2n.test(399,10, 0.74, sig.level=0.05, alternative=c("greater"))


## 1. Effect Size RMSPs vs. SPIN (AMXPs+NXPs)
mes(mean(spin),mean(rmsps), sd(spin), sd(rmsps), 29, 399, level=95)



### Plot the Cullen and Frey graph
library(fitdistrplus)
descdist(rmsps, discrete=FALSE, boot=500)



### Fit a distribution and plot results
# Fit
fit_w  <- fitdist(my_data, "weibull")
#fit_g  <- fitdist(my_data, "gamma")
fit_n  <- fitdist(my_data, "norm")
#fit_u  <- fitdist(my_data, "unif")
fit_ln <- fitdist(my_data, "lnorm")
summary(fit_ln)

# Plot 
par(mfrow=c(2,2))
plot.legend <- c("Weibull", "Lognormal", "Normal")
denscomp(list(fit_w, fit_ln, fit_n), legendtext = plot.legend,xlab= "Spin Frequency (Hz)", lwd=3)
cdfcomp (list(fit_w, fit_ln, fit_n), legendtext = plot.legend,xlab= "Spin Frequency (Hz)", lwd=3)
qqcomp  (list(fit_w, fit_ln, fit_n), legendtext = plot.legend, lwd=1)
ppcomp  (list(fit_w, fit_ln, fit_n), legendtext = plot.legend, lwd=1)

## Check goodness of the three fits
## Look at the BIC. If Delta_BIC ~0-2, not much to say
### If Delta_Bic between 2-6 positive
### If Delta_Bic between 6-10 strong
### If Delta_Bic >10 very strong (from https://en.wikipedia.org/wiki/Bayesian_information_criterion)
library(fitdistrplus)
gofstat(list(fit_w, fit_n,fit_ln))




### Create EPS
# setEPS()
# postscript("dist.eps")
# Plot stuff and then:
# dev.off()
####


###### Simulation of distribution with spin-down ###
nudot <-scan("nudot.txt")
meansamp <- vector("numeric", 10L)
shapvec <- vector("numeric", 10L)
maximum <- vector("numeric", 10L)
###dipvec <- vector("numeric", 10L)
for(i in 1:1000) {
#rage <- rnorm(100000, mean=mean(log10(age)), sd=sd(log10(age)))
#rrage <- 10**rage
rrage <- runif(100000, min=1e8, max=1e10)
rnudot <- rnorm(100000, mean=mean(log10(nudot)), sd=sd(log10(nudot)))
rrnudot <- 10**rnudot
delta <- rrage*rrnudot*86400*365
x<-rweibull(100000, shape=2, scale=292)
y <- x+delta
###hist(y, xlim=c(100,1000), breaks=1000)
## Select all elements in y with value <1000
y <- y[y<700]
mean(y)
sd(y)
median(y)
### Randomly sample 19 elements from y (19 because of the slow-population)
sample <- sample(y, size=19)
maximum[i] <- max(sample)
meansamp[i] <- mean(sample)
shapvec[i] <- shapiro.test(sample)$p.value
###dipvec[i] <-dip.test(sample)$p.value
}
shapvec[shapvec<0.05]
##########################################################


### Spin up backward radiomsp population

meansamp <- vector("numeric", 10L)
shapvec <- vector("numeric", 10L)
sdsamp <- vector("numeric", 10L)
delta <- vector("numeric",10L)
y <- vector("numeric",10L)

rmsp <- scan("../RMSP_NOSPIDERS_EXCLUDED.txt")


pulsnudot <- scan("nudot.txt")
pulsspin <- sample(rmsp, size=100)

for(j in 1:1000)
{     
for(i in 1:100) {
###### Simulation of distribution with spin-down
rrage <- runif(1, min=1e8, max=1e10)
delta[i] <- rrage*pulsnudot[i]*86400*365
y[i] = pulsspin[i] - delta[i]
}
y <- y[y<1500]
y <- y[y>0]
### Randomly sample 29 elements from y
sample<- sample(y, size=29)
meansamp[j] <- mean(sample)
sdsamp[j] <- sd(sample)
shapvec[j] <- shapiro.test(sample)$p.value
}


##############

### Bootstrap
library(mixtools)
growth <- spin
growth.boot <- boot.comp(growth,max.comp=5,mix.type="normalmix",maxit=400,epsilon=1e-2, B=100000)

#############

### Calculate BIC distribution from mixture model (1 to 3 components in the example below)

library(mclust)
bicn <- mclustBIC(spin, G=1:3)
plot(bicn)



###### Simulation of distribution with GW spin-down ###

meansamp <- vector("numeric", 10L)
shapvec <- vector("numeric", 10L)
maximum <- vector("numeric", 10L)
###dipvec <- vector("numeric", 10L)
## uniform dist dM between 1e20 and 1e22 g
rdM <- runif(100000,min=0.01, max=1) 
dT = 100*1e6*rdM*0.5
dTq = 0.03*dT
q22 =3e35*dTq/1e5
for(i in 1:1000) {
x<- rnorm(10000, mean=mean(spin), sd=sd(spin))
x <- x[x>100]
rrage <- runif(10000, min=1e8, max=1e9)
rrnudot <- -1.4e-13*(x/500.0)**5*(q22/1e37)**2
delta <- rrage*rrnudot*86400*365
y <- x+delta
y <- y[y>100]
mean(y)
sd(y)
median(y)
### Randomly sample 29 elements from y
sample <- sample(y, size=29)
maximum[i] <- max(sample)
meansamp[i] <- mean(sample)
shapvec[i] <- shapiro.test(sample)$p.value
}

##########################################################
