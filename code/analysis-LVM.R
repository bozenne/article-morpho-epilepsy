## * 1- load packages
library(lava)
library(lavaSearch2)
library(qqtest)
library(multcomp)
library(data.table)

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
}else{
    path <- ""
}
dfMDD <- readRDS(file.path(path,"data","dfMDD.rds"))

## * 3- exploratory (univariate)

## linear models
## WARNING: due to missing data, model are fitted on different dataset

## ipsilateral
ls.lm <- list(ipsi.ofc = lm(ipsi.ofc~Group2+Age, data = dfMDD),
             ipsi.fusi = lm(ipsi.fusi~Group2+Age, data = dfMDD),
             ipsi.insula.c = lm(ipsi.insula.c~Group2+Age, data = dfMDD),
             ipsi.ros.acc=lm(ipsi.ros.acc~Group2+Age, data=dfMDD),
             ipsi.pos.cc = lm(ipsi.pos.cc~Group2+Age, data = dfMDD)
             )
lapply(ls.lm, function(iLM){summary(iLM)$coef})

## contralateral 
ls.lmc <- list(con.ofc = lm(con.ofc~Group2+Age, data = dfMDD),
              con.fusi = lm(con.fusi~Group2+Age, data = dfMDD),
              con.insula.c = lm(con.insula.c~Group2+Age, data = dfMDD),
              con.ros.acc=lm(con.ros.acc~Group2+Age, data=dfMDD),
              con.pos.cc = lm(con.pos.cc~Group2+Age, data = dfMDD)
              )
lapply(ls.lmc, function(iLM){summary(iLM)$coef})

## library(LMMstar)
## summary(lmm(con.insula.c~Group2+Age, data = dfMDD,
##             structure = IND(~Group2)))

## * 4- blue print (LVM)
## WARNING: we decided to mix thickness and volume - may not make sense

## ** a) define the simplified LVM (to solve convergence issues)
## which is just a list of linear regressions
## WARNING: this is a wrong model just used for initialization of an algorithm
mLVM0.blueIP <- lvm(ipsi.ofc      ~ Age + 1*eta,
                    ipsi.fusi     ~ Age + 1*eta,
                    ipsi.insula.c ~ Age + 1*eta,
                    ipsi.ros.acc  ~ Age + 1*eta,
                    ipsi.pos.cc   ~ Age + 1*eta,
                    ipsi.khippo   ~ Age + kETIV + 1*eta)
latent(mLVM0.blueIP) <- ~eta

mLVM0.blueCO <- lvm(con.ofc       ~ Age + 1*eta,
                    con.fusi      ~ Age + 1*eta,
                    con.insula.c  ~ Age + 1*eta,
                    con.ros.acc   ~ Age + 1*eta,
                    con.pos.cc    ~ Age + 1*eta,
                    con.khippo    ~ Age + kETIV + 1*eta)
latent(mLVM0.blueCO) <- ~eta

## ** b) define the LVM
## ipsilateral hemisphere

## with group (under the alternative i.e. model of interest) SHOULD ALL THEE GROUPS BE INCLUDED HERE ? 
mLVM.blueIP <- lvm(ipsi.ofc      ~ Age + eta,
                   ipsi.fusi     ~ Age + eta,
                   ipsi.insula.c ~ Age + eta,
                   ipsi.ros.acc  ~ Age + eta,
                   ipsi.pos.cc   ~ Age + eta,
                   ipsi.khippo   ~ Age + eta + kETIV, ## ETIV should only affect volumes not thickness
                   eta           ~ Group2)
latent(mLVM.blueIP) <- ~eta
#covariance(mLVM.blueIP) <- ipsi.sup.fron~ipsi.mid.temp

## controlateral hemisphere
mLVM.blueCO <- lvm(con.ofc      ~ Age + eta,
                   con.fusi     ~ Age + eta,
                   con.insula.c ~ Age + eta,
                   con.ros.acc  ~ Age + eta,
                   con.pos.cc   ~ Age + eta,
                   con.khippo   ~ Age + eta + kETIV, ## ETIV should only affect volumes not thickness
                   eta          ~ Group2)
latent(mLVM.blueCO)<-~eta
#covariance(mLVM.blueCO) <- con.sup.fron~con.ros.mid.fron

## sanity check: spelling mistakes in variables?
## all(manifest(mLVM.blueIP) %in% names(dfMDD))
## all(manifest(mLVM.blueCO) %in% names(dfMDD))

## ** c) estimate the LVM
dfMDD$Group2 <- relevel(dfMDD$Group2, "C")

## NROW(dfMDD)
eLVM0.blueIP <- estimate(mLVM0.blueIP, data=dfMDD)

eLVM.blueIP <- estimate(mLVM.blueIP, data = dfMDD,
                        control = list(constrain = TRUE, start = coef(eLVM0.blueIP)))
## logLik(eLVM.blueIP)
## 'log Lik.' 159.1133 (df=27)

eLVM0.blueCO <- estimate(mLVM0.blueCO, data = dfMDD)

eLVM.blueCO <- estimate(mLVM.blueCO, data = dfMDD,
                        control=list(constrain=TRUE, start = coef(eLVM0.blueCO)))

## logLik(eLVM.blueCO)
## 'log Lik.' 236.1798 (df=27)

## ** d) model diagnostic
## LVM assumes a somehow similar group effect across regions
##             a flexible but not fully unstructured correlation pattern         
## test missing links (e.g. region-specific group effect or complex correlation pattern)
modelsearch(eLVM.blueIP)
 ## 7.852    0.005076  ipsi.ofc~~Group2MDDP        0.4873  0.05753
 ## 7.852    0.005076  ipsi.ofc~Group2MDDP         0.4873  0.05753
 ## 7.852    0.005076  Group2MDDP~ipsi.ofc         0.4873  0.05753
 ## 7.93     0.004862  ipsi.pos.cc~~kETIV          0.4813  0.05753
 ## 7.93     0.004862  ipsi.pos.cc~kETIV           0.4813  0.05753
 ## 7.93     0.004862  kETIV~ipsi.pos.cc           0.4813  0.05753
 ## 11.11    0.0008597 ipsi.insula.c~~kETIV        0.08769 0.02923
 ## 11.11    0.0008597 ipsi.insula.c~kETIV         0.08769 0.02923
 ## 11.11    0.0008597 kETIV~ipsi.insula.c         0.08769 0.02923
##set.seed(10) ## fix initial conditions for random draw

par(mfrow = c(2,3))
xx <- lapply(endogenous(eLVM.blueIP), function(iName){ ## iName <- endogenous(eLVM.blueIP)[1]
    qqtest(na.omit(residuals(eLVM.blueIP)[,iName]), main = iName)
})

modelsearch(eLVM.blueCO)
 ## 6.875    0.008739 con.ofc~~con.insula.c     0.8652 0.1486
 ## 6.875    0.008739 con.ofc~con.insula.c      0.8652 0.1486
 ## 6.875    0.008739 con.insula.c~con.ofc      0.8652 0.1486
 ## 7.413    0.006475 con.insula.c~~con.khippo  0.6604 0.1486
 ## 7.413    0.006475 con.insula.c~con.khippo   0.6604 0.1486
 ## 7.413    0.006475 con.khippo~con.insula.c   0.6604 0.1486

par(mfrow = c(2,3))
xx <- lapply(endogenous(eLVM.blueCO), function(iName){ ## iName <- endogenous(eLVM.blueIP)[1]
    qqtest(na.omit(residuals(eLVM.blueCO)[,iName]), main = iName)
})


## ** e) Post hoc LVM

## *** Ipsilateral
mLVM.blueIP2 <- mLVM.blueIP
regression(mLVM.blueIP2) <- ipsi.ofc~Group2

eLVM.blueIP2 <- estimate(mLVM.blueIP2, data = dfMDD,
                         control = list(constrain = TRUE, start = coef(eLVM.blueIP)))
modelsearch(eLVM.blueIP2)
## 4.286    0.03843   ipsi.ros.acc~~kETIV           1       0.4099 
## 4.286    0.03843   ipsi.ros.acc~kETIV            1       0.4099 
## 4.286    0.03843   kETIV~ipsi.ros.acc            1       0.4099 
## 8.337    0.003885  ipsi.pos.cc~~kETIV            0.3613  0.06216
## 8.337    0.003885  ipsi.pos.cc~kETIV             0.3613  0.06216
## 8.337    0.003885  kETIV~ipsi.pos.cc             0.3613  0.06216
## 11.6     0.0006586 ipsi.insula.c~~kETIV          0.06322 0.02107
## 11.6     0.0006586 ipsi.insula.c~kETIV           0.06322 0.02107
## 11.6     0.0006586 kETIV~ipsi.insula.c           0.06322 0.02107

compare(eLVM.blueIP,eLVM.blueIP2)
## chisq = 17.276, df = 2, p-value = 0.0001773
## sample estimates:
## log likelihood (model 1) log likelihood (model 2) 
##                 159.1133                 167.7511 

## *** Contralateral
mLVM.blueCO2 <- mLVM.blueCO
covariance(mLVM.blueCO2) <- con.khippo~con.insula.c

eLVM.blueCO2 <- estimate(mLVM.blueCO2, data = dfMDD,
                         control = list(constrain = TRUE, start = coef(eLVM.blueCO)))

modelsearch(eLVM.blueCO2)
## 5.447    0.01961 con.ofc~~con.insula.c        1    0.2282
## 5.447    0.01961 con.ofc~con.insula.c         1    0.2282
## 5.447    0.01961 con.insula.c~con.ofc         1    0.2282


compare(eLVM.blueCO,eLVM.blueCO2)
## chisq = 7.256, df = 1, p-value = 0.007066
## sample estimates:
## log likelihood (model 1) log likelihood (model 2) 
##                 236.1798                 239.8079 


## ** f) export
saveRDS(eLVM.blueIP, file = file.path(path,"data","lvm-blueIP.rds"))
saveRDS(eLVM.blueCO, file = file.path(path,"data","lvm-blueCO.rds"))

saveRDS(eLVM.blueIP2, file = file.path(path,"data","lvm-blueIP2.rds"))
saveRDS(eLVM.blueCO2, file = file.path(path,"data","lvm-blueCO2.rds"))

## * 6- Connectivity
## testing whether there is at all a difference in connectivity or global binding between depression groups


dfMDD.NNA <- dfMDD[rowSums(is.na(dfMDD))==0,]

dfMDD.NNA$ipsi.ofc <- scale(dfMDD.NNA$ipsi.ofc)
dfMDD.NNA$ipsi.fusi <- scale(dfMDD.NNA$ipsi.fusi)
dfMDD.NNA$ipsi.insula.c <- scale(dfMDD.NNA$ipsi.insula.c)
dfMDD.NNA$ipsi.ros.acc <- scale(dfMDD.NNA$ipsi.ros.acc)
dfMDD.NNA$ipsi.pos.cc <- scale(dfMDD.NNA$ipsi.pos.cc)
dfMDD.NNA$ipsi.khippo <- scale(dfMDD.NNA$ipsi.khippo)

dfMDD.NNA$con.ofc <- scale(dfMDD.NNA$con.ofc)
dfMDD.NNA$con.fusi <- scale(dfMDD.NNA$con.fusi)
dfMDD.NNA$con.insula.c <- scale(dfMDD.NNA$con.insula.c)
dfMDD.NNA$con.ros.acc <- scale(dfMDD.NNA$con.ros.acc)
dfMDD.NNA$con.pos.cc <- scale(dfMDD.NNA$con.pos.cc)
dfMDD.NNA$con.khippo <- scale(dfMDD.NNA$con.khippo)

## ** a) define the LVM
## *** ipsi
## 0 model
mLVM.meanIP0.HC <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                       ipsi.fusi[mu2:sigma2]     ~ beta2*Age + 1*eta,
                       ipsi.insula.c[mu3:sigma3] ~ beta3*Age + 1*eta,
                       ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + 1*eta,
                       ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + 1*eta,
                       ipsi.khippo[mu6:sigma6]  ~ beta6*Age + 1*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                       eta[alpha1:tau1] ~ 1
                       )
latent(mLVM.meanIP0.HC)<-~eta

mLVM.meanIP0.Dpost <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                          ipsi.fusi[mu2:sigma2]     ~ beta2*Age + 1*eta,
                          ipsi.insula.c[mu3:sigma3] ~ beta3*Age + 1*eta,
                          ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + 1*eta,
                          ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + 1*eta,
                          ipsi.khippo[mu6:sigma6]  ~ beta6*Age + 1*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                          eta[alpha1:tau1] ~ 1
                          )
latent(mLVM.meanIP0.Dpost)<-~eta

mLVM.meanIP0.Dpre <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                         ipsi.fusi[mu2:sigma2]     ~ beta2*Age + 1*eta,
                         ipsi.insula.c[mu3:sigma3] ~ beta3*Age + 1*eta,
                         ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + 1*eta,
                         ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + 1*eta,
                         ipsi.khippo[mu6:sigma6]  ~ beta6*Age + 1*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                         eta[alpha1:tau1] ~ 1
                         )
latent(mLVM.meanIP0.Dpre)<-~eta

## mean model
mLVM.meanIP.HC <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                       ipsi.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                       ipsi.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                       ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                       ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                       ipsi.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                       eta[alpha1:tau1] ~ 1
                       )
latent(mLVM.meanIP.HC)<-~eta

mLVM.meanIP.Dpost <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                          ipsi.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                          ipsi.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                          ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                          ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                          ipsi.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                          eta[alpha1.bis:tau1] ~ 1
                          )
latent(mLVM.meanIP.Dpost)<-~eta

mLVM.meanIP.Dpre <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                         ipsi.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                         ipsi.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                         ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                         ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                         ipsi.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                         eta[alpha1.ter:tau1] ~ 1
                         )
latent(mLVM.meanIP.Dpre)<-~eta

## connectivity model
mLVM.strataIP.HC <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                        ipsi.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                        ipsi.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                        ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                        ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                        ipsi.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                        eta[alpha1:tau1] ~ 1
                        )
latent(mLVM.strataIP.HC)<-~eta

mLVM.strataIP.Dpost <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                           ipsi.fusi[mu2:sigma2]     ~ beta2*Age + lambda2.bis*eta,
                           ipsi.insula.c[mu3:sigma3] ~ beta3*Age + lambda3.bis*eta,
                           ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4.bis*eta,
                           ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5.bis*eta,
                           ipsi.khippo[mu6:sigma6]   ~ beta6*Age + lambda6.bis*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                           eta[alpha1.bis:tau1.bis] ~ 1
                           )
latent(mLVM.strataIP.Dpost)<-~eta

mLVM.strataIP.Dpre <- lvm(ipsi.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                          ipsi.fusi[mu2:sigma2]     ~ beta2*Age + lambda2.ter*eta,
                          ipsi.insula.c[mu3:sigma3] ~ beta3*Age + lambda3.ter*eta,
                          ipsi.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4.ter*eta,
                          ipsi.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5.ter*eta,
                          ipsi.khippo[mu6:sigma6]   ~ beta6*Age + lambda6.ter*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                          eta[alpha1.ter:tau1.ter] ~ 1
                          )
latent(mLVM.strataIP.Dpre)<-~eta

## *** contralateral
## 0 model
mLVM.meanCO0.HC <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                       con.fusi[mu2:sigma2]     ~ beta2*Age + 1*eta,
                       con.insula.c[mu3:sigma3] ~ beta3*Age + 1*eta,
                       con.ros.acc[mu4:sigma4]  ~ beta4*Age + 1*eta,
                       con.pos.cc[mu5:sigma5]   ~ beta5*Age + 1*eta,
                       con.khippo[mu6:sigma6]  ~ beta6*Age + 1*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                       eta[alpha1:tau1] ~ 1
                       )
latent(mLVM.meanCO0.HC)<-~eta

mLVM.meanCO0.Dpost <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                          con.fusi[mu2:sigma2]     ~ beta2*Age + 1*eta,
                          con.insula.c[mu3:sigma3] ~ beta3*Age + 1*eta,
                          con.ros.acc[mu4:sigma4]  ~ beta4*Age + 1*eta,
                          con.pos.cc[mu5:sigma5]   ~ beta5*Age + 1*eta,
                          con.khippo[mu6:sigma6]  ~ beta6*Age + 1*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                          eta[alpha1:tau1] ~ 1
                          )
latent(mLVM.meanCO0.Dpost)<-~eta

mLVM.meanCO0.Dpre <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                         con.fusi[mu2:sigma2]     ~ beta2*Age + 1*eta,
                         con.insula.c[mu3:sigma3] ~ beta3*Age + 1*eta,
                         con.ros.acc[mu4:sigma4]  ~ beta4*Age + 1*eta,
                         con.pos.cc[mu5:sigma5]   ~ beta5*Age + 1*eta,
                         con.khippo[mu6:sigma6]  ~ beta6*Age + 1*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                         eta[alpha1:tau1] ~ 1
                         )
latent(mLVM.meanCO0.Dpre)<-~eta

## mean model
mLVM.meanCO.HC <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                       con.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                       con.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                       con.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                       con.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                       con.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                       eta[alpha1:tau1] ~ 1
                       )
latent(mLVM.meanCO.HC)<-~eta

mLVM.meanCO.Dpost <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                          con.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                          con.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                          con.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                          con.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                          con.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                          eta[alpha1.bis:tau1] ~ 1
                          )
latent(mLVM.meanCO.Dpost)<-~eta

mLVM.meanCO.Dpre <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                         con.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                         con.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                         con.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                         con.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                         con.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                         eta[alpha1.ter:tau1] ~ 1
                         )
latent(mLVM.meanCO.Dpre)<-~eta

## connectivity model
mLVM.strataCO.HC <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                        con.fusi[mu2:sigma2]     ~ beta2*Age + lambda2*eta,
                        con.insula.c[mu3:sigma3] ~ beta3*Age + lambda3*eta,
                        con.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4*eta,
                        con.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5*eta,
                        con.khippo[mu6:sigma6]   ~ beta6*Age + lambda6*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                        eta[alpha1:tau1] ~ 1
                        )
latent(mLVM.strataCO.HC)<-~eta

mLVM.strataCO.Dpost <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                           con.fusi[mu2:sigma2]     ~ beta2*Age + lambda2.bis*eta,
                           con.insula.c[mu3:sigma3] ~ beta3*Age + lambda3.bis*eta,
                           con.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4.bis*eta,
                           con.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5.bis*eta,
                           con.khippo[mu6:sigma6]   ~ beta6*Age + lambda6.bis*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                           eta[alpha1.bis:tau1.bis] ~ 1
                           )
latent(mLVM.strataCO.Dpost)<-~eta

mLVM.strataCO.Dpre <- lvm(con.ofc[0:sigma1]        ~ beta1*Age + 1*eta,
                          con.fusi[mu2:sigma2]     ~ beta2*Age + lambda2.ter*eta,
                          con.insula.c[mu3:sigma3] ~ beta3*Age + lambda3.ter*eta,
                          con.ros.acc[mu4:sigma4]  ~ beta4*Age + lambda4.ter*eta,
                          con.pos.cc[mu5:sigma5]   ~ beta5*Age + lambda5.ter*eta,
                          con.khippo[mu6:sigma6]   ~ beta6*Age + lambda6.ter*eta + gamma*kETIV, ## ETIV should only affect volumes not thickness
                          eta[alpha1.ter:tau1.ter] ~ 1
                          )
latent(mLVM.strataCO.Dpre)<-~eta

## ** b) estimate the LVM
## *** ipsi
## sanity check
eLVM.meanIP00 <- estimate(mLVM.meanIP0.HC, data = dfMDD.NNA, control = list(constrain = TRUE, trace = FALSE))
## summary(eLVM.meanIP00)
logLik(eLVM.meanIP00)
## 'log Lik.' -668.0951 (df=20)

eLVM.meanIP0 <- estimate(list(mLVM.meanIP0.HC,mLVM.meanIP0.Dpost,mLVM.meanIP0.Dpre), 
                        data = split(dfMDD.NNA, dfMDD.NNA$Group),
                        control = list(constrain = TRUE, trace = FALSE))
## summary(eLVM.meanIP0)
logLik(eLVM.meanIP0)
## 'log Lik.' -668.0951 (df=20)


## mean model
mLVM.meanIP.bis <- mLVM.meanIP.HC
regression(mLVM.meanIP.bis) <- eta ~ Group
logLik(estimate(mLVM.meanIP.bis, data = dfMDD.NNA, control = list(constrain = TRUE, trace = FALSE)))
## 'log Lik.' -655.7329 (df=27)

eLVM.meanIP <- estimate(list(mLVM.meanIP.HC,mLVM.meanIP.Dpost,mLVM.meanIP.Dpre), 
                        data = split(dfMDD.NNA, dfMDD.NNA$Group),
                        control = list(constrain = TRUE, trace = FALSE, start = coef(eLVM.meanIP0)))
## summary(eLVM.meanIP)
logLik(eLVM.meanIP)
## 'log Lik.' -655.7329 (df=27)

## connectivity model
eLVM.strataIP <- estimate(list(mLVM.strataIP.HC,mLVM.strataIP.Dpost,mLVM.strataIP.Dpre), 
                        data = split(dfMDD.NNA, dfMDD.NNA$Group),
                        control = list(constrain = TRUE, trace = FALSE, start = coef(eLVM.meanIP)))
## summary(eLVM.meanIP)
logLik(eLVM.strataIP)
## 'log Lik.' -655.7329 (df=27)

lava::compare(eLVM.strataIP, eLVM.meanIP)
## data:  
## chisq = 22.353, df = 12, p-value = 0.03375
## sample estimates:
## log likelihood (model 1) log likelihood (model 2) 
##                -644.5564                -655.7329 


## *** contralateral

## sanity check
eLVM.meanCO00 <- estimate(mLVM.meanCO0.HC, data = dfMDD.NNA, control = list(constrain = TRUE, trace = FALSE))
## summary(eLVM.meanCO00)
logLik(eLVM.meanCO00)
## 'log Lik.' -655.9556 (df=20)

eLVM.meanCO0 <- estimate(list(mLVM.meanCO0.HC,mLVM.meanCO0.Dpost,mLVM.meanCO0.Dpre), 
                        data = split(dfMDD.NNA, dfMDD.NNA$Group),
                        control = list(constrain = TRUE, trace = FALSE))
## summary(eLVM.meanCO0)
logLik(eLVM.meanCO0)
## ''log Lik.' -655.9556 (df=20)


## mean model
mLVM.meanCO.bis <- mLVM.meanCO.HC
regression(mLVM.meanCO.bis) <- eta ~ Group
logLik(estimate(mLVM.meanCO.bis, data = dfMDD.NNA, control = list(constrain = TRUE, trace = FALSE)))
## 'log Lik.' -633.0835 (df=27)

eLVM.meanCO <- estimate(list(mLVM.meanCO.HC,mLVM.meanCO.Dpost,mLVM.meanCO.Dpre), 
                        data = split(dfMDD.NNA, dfMDD.NNA$Group),
                        control = list(constrain = TRUE, trace = FALSE, start = coef(eLVM.meanCO0)))
## summary(eLVM.meanCO)
logLik(eLVM.meanCO)
## 'log Lik.' -633.0835 (df=27)

## connectivity model
eLVM.strataCO <- estimate(list(mLVM.strataCO.HC,mLVM.strataCO.Dpost,mLVM.strataCO.Dpre), 
                        data = split(dfMDD.NNA, dfMDD.NNA$Group),
                        control = list(constrain = TRUE, trace = FALSE, start = coef(eLVM.meanCO)))
## summary(eLVM.meanCO)
logLik(eLVM.strataCO)
## 'log Lik.' -624.094 (df=39)



## ** c) compare models
## depression effect only on the mean vs. on both the mean and correlation
lava::compare(eLVM.strataIP, eLVM.meanIP)
## data:  
## chisq = 22.353, df = 12, p-value = 0.03375
## sample estimates:
## log likelihood (model 1) log likelihood (model 2) 
##                -644.5564                -655.7329 

lava::compare(eLVM.strataCO, eLVM.meanCO)
## 	- Likelihood ratio test -

## data:  
## chisq = 17.979, df = 12, p-value = 0.1163
## sample estimates:
## log likelihood (model 1) log likelihood (model 2) 
##                -624.0940                -633.0835 


## ** d) export

con.newname <- setdiff(names(coef(eLVM.strataIP)),names(coef(eLVM.meanIP)))
con.newname
##  [1] "ipsi.fusi~eta@2"     "ipsi.insula.c~eta@2" "ipsi.ros.acc~eta@2" 
##  [4] "ipsi.pos.cc~eta@2"   "ipsi.khippo~eta@2"   "eta~~eta@2"         
##  [7] "ipsi.fusi~eta@3"     "ipsi.insula.c~eta@3" "ipsi.ros.acc~eta@3" 
## [10] "ipsi.pos.cc~eta@3"   "ipsi.khippo~eta@3"   "eta~~eta@3"         

saveRDS(eLVM.meanIP, file.path(path,"data","lvm-meanIP.rds"))
saveRDS(eLVM.meanCO, file.path(path,"data","lvm-meanCO.rds"))
saveRDS(eLVM.strataIP, file.path(path,"data","lvm-strataIP.rds"))
saveRDS(eLVM.strataCO, file.path(path,"data","lvm-strataCO.rds"))



