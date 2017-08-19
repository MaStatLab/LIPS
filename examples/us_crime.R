### US crime data
library(MASS)
data(UScrime)
attach(UScrime)

X = cbind(log(M),So,log(Ed),log(Po1),log(Po2),log(LF),log(M.F), log(Pop),log(NW),log(U1),log(U2),log(GDP),log(Ineq),log(Prob),log(Time))
Y = log(y)


set.seed(12345)

### Run one of the below

## k = 3 and 50000 particles
ans = lips(Y,X,T=ncol(X),kstep = 3, alpha = length(Y), nparticle = 50000,resample.param=0)

## k = 4 and 5000 particles
ans = lips(Y,X,T=ncol(X),kstep = 4, alpha = length(Y), nparticle = 5000,resample.param=0)


models=ans$models
weights = exp(ans$logweights)
incl.probs = apply(models*weights,2,sum)/sum(weights)

incl.probs.true = c(0.7636530,0.2093355,0.9211455,0.6661690,0.4224025,0.1542015,0.1699795,0.3013900,0.5702085,0.1897870,0.5065560,0.2856305,0.9924710,0.7852880,0.2865580)
barplot(rbind(incl.probs.true,incl.probs),xlab="Predictors",col=c("darkblue","red"),legend=c("Exact","LIPS"), beside=TRUE, main="Posterior marginal inclusion probabilities")
axis(1,at=seq(from=2,by=3,length=15),labels=paste("X",1:15,sep=""))

