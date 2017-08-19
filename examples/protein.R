library(BAS)
data(protein)

protein.lm <- lm(prot.act4 ~
        buf + buf*pH + buf*NaCl + buf*con + buf*ra +
         buf*det + buf*MgCl2 + buf*temp +
        pH + pH*NaCl + pH*con + pH*ra + pH*det +
         pH*MgCl2 +pH*temp +
        NaCl + NaCl*con + NaCl*ra + NaCl*det +
         NaCl*MgCl2 + NaCl*temp +
        con + con*ra + con*det +
         con*MgCl2 +con*temp +
        ra + ra*det + ra*MgCl2 + ra*temp +
        det + det*MgCl2 + det*temp +
        MgCl2 + MgCl2*temp + I(NaCl^2) + I(pH^2) + I(con^2) + I(temp^2),
        data=protein,x=T)

protein.designmat <- protein.lm$x
proteindata <-  data.frame(cbind(protein.designmat[,-1],protein$prot.act4));
names(proteindata)[89] <- "prot.act4"


Y=proteindata$prot.act4;
X=protein.designmat[,-1]
colnames(X)=paste("X",1:88,sep="")

### Run one of the below

## k = 3 and 10000 particles with resampling
ans = lips(Y,X,kstep = 3, nparticle = 1000, resample.param = 1)

models=ans$models
weights = exp(ans$logweights)
incl.probs = apply(models*weights,2,sum)/sum(weights)
barplot(incl.probs)

