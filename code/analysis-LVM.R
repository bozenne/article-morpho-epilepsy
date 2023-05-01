## * 1- load packages
library(lava)
library(lavaSearch2)
library(qqtest)
library(multcomp)
library(data.table)
library(ggpubr) ## install.packages("ggpubr")
library(officer)

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
}else{
    path <- ""
}
dfMDD.rds <- readRDS(file.path(path,"data","dfMDD.rds"))

## * 3- exploratory (univariate)

## graphical display
ggplot(dfMDD.long.ipsi, aes(x=value)) + facet_wrap(~variable, scales = "free") + geom_histogram()

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


## ** e) Statistical inference

#### test whether hippocampus volume correlated with thickness
coef(eLVM.blueIP, type = 9)["ipsi.khippo~eta",]
#### note from the sumary all thickness correlates together
## to get the correlation between region AFTER adjusting for covariates effect (i.e. remove group, age, ETIV effects)
M.cor <- cov2cor(attr(predict(eLVM.blueIP),"cond.var")) ##ipsilateral 
M.cor
M.cor <- cov2cor(attr(predict(eLVM.blueCO),"cond.var")) ##contralateral
M.cor

##is it possible to get the above, but  one for each group ? 

#### over all regions ipsilateral
coef(eLVM.blueIP, type = 9)["eta~Group2MDD",] ##is this for all three groups or ?##
#Estimate  Std. Error     Z-value     P-value 
#-0.07939275  0.03098253 -2.56250058  0.01039214
coef(eLVM.blueIP, type = 9)["eta~Group2MDDP",]
#Estimate   Std. Error      Z-value      P-value 
#-0.087706566  0.030999073 -2.829328708  0.004664576 

#### over all regions contralaterallateral
coef(eLVM.blueCO, type = 9)["eta~Group2MDD",]
#Estimate   Std. Error      Z-value      P-value 
#-0.106465486  0.037980340 -2.803173595  0.005060242 
coef(eLVM.blueCO, type = 9)
coef(eLVM.blueCO, type = 9)["eta~Group2MDDP",]
#Estimate  Std. Error     Z-value     P-value 
#-0.07552045  0.03744400 -2.01689050  0.04370694

#### region specific (multiply by loading to make it region specific and return to original scale) GROUP MDD
## effect in ipsi.mid.temp
effects(eLVM.blueIP, ipsi.ofc~Group2MDD)
effects(eLVM.blueIP, ipsi.fusi~Group2MDD)
effects(eLVM.blueIP, ipsi.insula.c~Group2MDD)
effects(eLVM.blueIP, ipsi.ros.acc~Group2MDD)
effects(eLVM.blueIP, ipsi.pos.cc~Group2MDD)
effects(eLVM.blueIP, ipsi.khippo~Group2MDD)

#### region specific (multiply by loading to make it region specific and return to original scale) GROUP2MDD
## effect in con.mid.temp 
effects(eLVM.blueIP, con.ofc~Group2MDD)
effects(eLVM.blueIP, con.fusi~Group2MDD)
effects(eLVM.blueIP, con.insula.c~Group2MDD)
effects(eLVM.blueIP, con.ros.acc~Group2MDD)
effects(eLVM.blueIP, con.pos.cc~Group2MDD)
effects(eLVM.blueIP, con.khippo~Group2MDD)

## ** f) export
saveRDS(eLVM.blueIP, file = file.path(path,"data","lvm-blueIP.rds"))
saveRDS(eLVM.blueCO, file = file.path(path,"data","lvm-blueCO.rds"))

## * 6- Connectivity
## testing whether there is at all a difference in connectivity or global binding between depression groups


if(FALSE){ ## TO BE UPDATED

## a) define the simplified LVM (to solve convergence issues)

## not needed here instead we used rescaled data

## b) define the LVM
## ipsilateral hemisphere

## mLVM.blueIP2 is equivalent to mLVM.blueIP just with the rescaled data and ipsi.hippo insteand of ipsi.khippo
## this is because this dataset is easier to handle for the optimizer in complex model like the connectivity one (mLVM.conIP.D)
mLVM.blueIP2 <- lvm(ipsi.mid.temp~Age+eta,
                    ipsi.hippo~Age+eta+ETIV, 
                    ipsi.ofc~Age+Group2MDD+eta,
                    ipsi.pos.cc~Age+eta,
                    ipsi.ros.mid.fron~Age+eta,
                    eta~Group2)
latent(mLVM.blueIP2)<-~eta

mLVM.blueCO2 <- lvm(con.mid.temp~Age+eta,
                    con.hippo~Age+eta+ETIV, 
                    con.ofc~Age+Group2MDD+eta,
                    con.pos.cc~Age+eta,
                    con.ros.mid.fron~Age+eta,
                    eta~Group2)
latent(mLVM.blueCO2)<-~eta

## mLVM.blueIP2.HC and mLVM.blueIP2.D are together equivalent to mLVM.blueIP2
## this parametrisation enable to have group dependent loading thus connectivity model
mLVM.meanIP.HC <- lvm(ipsi.mid.temp[0:sigma1]~beta1*Age+1*eta,
                      ipsi.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+lambda2*eta,
                      ipsi.ofc[mu3:sigma3]~beta4*Age+lambda3*eta,
                      ipsi.pos.cc[mu4:sigma4]~beta5*Age+lambda4*eta,
                      ipsi.ros.mid.fron[mu5:sigma5]~beta8*Age+lambda5*eta,
                      eta[alpha1:tau1]~1)
latent(mLVM.meanIP.HC)<-~eta
mLVM.meanIP.D <- lvm(ipsi.mid.temp[0:sigma1]~beta1*Age+1*eta,
                     ipsi.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+lambda2*eta,
                     ipsi.ofc[mu3.bis:sigma3]~beta4*Age+lambda3*eta,
                     ipsi.pos.cc[mu4:sigma4]~beta5*Age+lambda4*eta,
                     ipsi.ros.mid.fron[mu5:sigma5]~beta8*Age+lambda5*eta,
                     eta[alpha1.bis:tau1]~1)
latent(mLVM.meanIP.D)<-~eta

mLVM.meanIP.DP <- lvm(ipsi.mid.temp[0:sigma1]~beta1*Age+1*eta,
                     ipsi.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+lambda2*eta,
                     ipsi.ofc[mu3.biz:sigma3]~beta4*Age+lambda3*eta,
                     ipsi.pos.cc[mu4:sigma4]~beta5*Age+lambda4*eta,
                     ipsi.ros.mid.fron[mu5:sigma5]~beta8*Age+lambda5*eta,
                     eta[alpha1.biz:tau1]~1)
latent(mLVM.meanIP.DP)<-~eta

mLVM.meanCO.HC <- lvm(con.mid.temp[0:sigma1]~beta1*Age+1*eta,
                      con.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+lambda2*eta,
                      con.ofc[mu3:sigma3]~beta4*Age+lambda3*eta,
                      con.pos.cc[mu4:sigma4]~beta5*Age+lambda4*eta,
                      con.ros.mid.fron[mu5:sigma5]~beta8*Age+lambda5*eta,
                      eta[alpha1:tau1]~1)
latent(mLVM.meanCO.HC)<-~eta
mLVM.meanCO.D <- lvm(con.mid.temp[0:sigma1]~beta1*Age+1*eta,
                     con.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+lambda2*eta,
                     con.ofc[mu3.bis:sigma3]~beta4*Age+lambda3*eta,
                     con.pos.cc[mu4:sigma4]~beta5*Age+lambda4*eta,
                     con.ros.mid.fron[mu5:sigma5]~beta8*Age+lambda5*eta,
                     eta[alpha1.bis:tau1]~1)
latent(mLVM.meanCO.D)<-~eta
##new bis=>biz##
mLVM.meanCO.DP <- lvm(con.mid.temp[0:sigma1]~beta1*Age+1*eta,
                      con.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+lambda2*eta,
                      con.ofc[mu3.biz:sigma3]~beta4*Age+lambda3*eta,
                      con.pos.cc[mu4:sigma4]~beta5*Age+lambda4*eta,
                      con.ros.mid.fron[mu5:sigma5]~beta8*Age+lambda5*eta,
                      eta[alpha1.biz:tau1]~1)
latent(mLVM.meanCO.DP)<-~eta

## connectivity model
mLVM.conIP.HC <- lvm(ipsi.mid.temp[0:sigma1]~beta1*Age+1*eta,
                     ipsi.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+eta,
                     ipsi.ofc[mu3:sigma3]~beta4*Age+eta,
                     ipsi.pos.cc[mu4:sigma4]~beta5*Age+eta,
                     ipsi.ros.mid.fron[mu5:sigma5]~beta8*Age+eta,
                     eta[alpha1:tau1]~1)
latent(mLVM.conIP.HC)<-~eta
mLVM.conIP.D <- lvm(ipsi.mid.temp[0:sigma1]~beta1*Age+1*eta,
                    ipsi.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+eta,
                    ipsi.ofc[mu3.bis:sigma3]~beta4*Age+eta,
                    ipsi.pos.cc[mu4:sigma4]~beta5*Age+eta,
                    ipsi.ros.mid.fron[mu5:sigma5]~beta8*Age+eta,
                    eta[alpha1.bis:tau1.bis]~1)
latent(mLVM.conIP.D)<-~eta
##new##
mLVM.conIP.DP <- lvm(ipsi.mid.temp[0:sigma1]~beta1*Age+1*eta,   
                    ipsi.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+eta,
                    ipsi.ofc[mu3.biz:sigma3]~beta4*Age+eta,
                    ipsi.pos.cc[mu4:sigma4]~beta5*Age+eta,
                    ipsi.ros.mid.fron[mu5:sigma5]~beta8*Age+eta,
                    eta[alpha1.biz:tau1.biz]~1)

mLVM.conCO.HC <- lvm(con.mid.temp[0:sigma1]~beta1*Age+1*eta,
                     con.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+eta,
                     con.ofc[mu3:sigma3]~beta4*Age+eta,
                     con.pos.cc[mu4:sigma4]~beta5*Age+eta,
                     con.ros.mid.fron[mu5:sigma5]~beta8*Age+eta,
                     eta[alpha1:tau1]~1)
latent(mLVM.conCO.HC)<-~eta

mLVM.conCO.D <- lvm(con.mid.temp[0:sigma1]~beta1*Age+1*eta,
                    con.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+eta,
                    con.ofc[mu3.bis:sigma3]~beta4*Age+eta,
                    con.pos.cc[mu4:sigma4]~beta5*Age+eta,
                    con.ros.mid.fron[mu5:sigma5]~beta8*Age+eta,
                    eta[alpha1.bis:tau1.bis]~1)
latent(mLVM.conCO.D)<-~eta

mLVM.conCO.DP <- lvm(con.mid.temp[0:sigma1]~beta1*Age+1*eta,
                    con.hippo[mu2:sigma2]~beta2*Age+gamma*ETIV+eta,
                    con.ofc[mu3.biz:sigma3]~beta4*Age+eta,
                    con.pos.cc[mu4:sigma4]~beta5*Age+eta,
                    con.ros.mid.fron[mu5:sigma5]~beta8*Age+eta,
                    eta[alpha1.biz:tau1.biz]~1)
latent(mLVM.conCO.DP)<-~eta

## b) estimate the LVM - added LVm.DP model 
eLVM.H0IP2 <- estimate(mLVM.conIP.HC,
                       data = dfMDD.rescale,
                       control = list(constrain = TRUE))

eLVM.blueIP2 <- estimate(mLVM.blueIP2,
                         data = dfMDD.rescale,
                         control = list(constrain = TRUE))

eLVM.meanIP <- estimate(list(mLVM.meanIP.HC,mLVM.meanIP.D,mLVM.meanIP.DP), #!!!!!#
                       data = split(dfMDD.rescale, dfMDD.rescale$Group),
                       control = list(constrain = TRUE, trace = FALSE))

eLVM.conIP <- estimate(list(mLVM.conIP.HC,mLVM.conIP.D,mLVM.conIP.DP), 
                       data = split(dfMDD.rescale, dfMDD.rescale$Group),
                       control = list(constrain = TRUE, trace = FALSE))

eLVM.blueCO2 <- estimate(mLVM.blueCO2,
                         data = dfMDD.rescale,
                         control = list(constrain = TRUE))

eLVM.meanCO <- estimate(list(mLVM.meanCO.HC,mLVM.meanCO.D,mLVM.meanCO.DP), 
                        data = split(dfMDD.rescale, dfMDD.rescale$Group),
                        control = list(constrain = TRUE, trace = FALSE))

eLVM.conCO <- estimate(list(mLVM.conCO.HC,mLVM.conCO.D,mLVM.conCO.DP), 
                       data = split(dfMDD.rescale, dfMDD.rescale$Group),
                       control = list(constrain = TRUE, trace = FALSE))


logLik(eLVM.meanIP) ## 'log Lik.' -399.6732 (df=23)
logLik(eLVM.blueIP2) ## 'log Lik.' -399.6732 (df=23)
logLik(eLVM.conIP) ## 'log Lik.' -397.2052 (df=28)
range(eigen(vcov(eLVM.blueIP2))$value) ## [1] 0.004396119 0.097991250 (ok - should be strictly positive)
range(eigen(vcov(eLVM.meanIP))$value) ## [1] -0.1533372  0.8620473 (not ok - should be strictly positive)
range(eigen(vcov(eLVM.conIP))$value) ## [1] -0.4327652  0.7246956 (not ok - should be strictly positive)


logLik(eLVM.meanCO) ## 'log Lik.' -385.6391 (df=23)
logLik(eLVM.blueCO2) ## 'log Lik.' -385.6391 (df=23)
logLik(eLVM.conCO) ## 'log Lik.' -381.4767 (df=28)

## c) compare models
## depression effect only on the mean vs. on both the mean and correlation
lava::compare(eLVM.conIP,eLVM.meanIP)
## 	- Likelihood ratio test -

## data:  
## chisq = 4.9359, df = 5, p-value = 0.4238
## sample estimates:
## log likelihood (model 1) log likelihood (model 2) 
##                -397.2052                -399.6732 

lava::compare(eLVM.conCO,eLVM.meanCO)
# - Likelihood ratio test -
#   
#   data:  
#   chisq = 8.3249, df = 5, p-value = 0.1392
# sample estimates:
#   log likelihood (model 1) log likelihood (model 2) 
# -381.4767                -385.6391 

con.newname <- setdiff(names(coef(eLVM.conIP)),names(coef(eLVM.meanIP)))
con.newname
## [1] "con.khippo~eta@2"       "ipsi.ofc~eta@2"         !!!! ###unsure here####
## [3] "ipsi.pos.cc~eta@2"       "ipsi.ros.mid.fron~eta@2"
## [5] "eta~~eta@2"

## extract coefficients
sigmaIP.mean <- coef(eLVM.meanIP)[paste0(endogenous(eLVM.meanIP),"~~",endogenous(eLVM.meanIP),"@1")]
sigmaIP.con <- coef(eLVM.conIP)[paste0(endogenous(eLVM.conIP),"~~",endogenous(eLVM.conIP),"@1")]
lambdaIP.mean2 <- c(setNames(1,paste0(endogenous(eLVM.blueIP2)[1],"~",latent(eLVM.blueIP2))),
                  coef(eLVM.blueIP2)[paste0(endogenous(eLVM.blueIP2)[-1],"~",latent(eLVM.blueIP2))])
lambdaIP.mean <- c(setNames(1,paste0(endogenous(eLVM.meanIP)[1],"~",latent(eLVM.meanIP),"@1")),
                 coef(eLVM.meanIP)[paste0(endogenous(eLVM.meanIP)[-1],"~",latent(eLVM.meanIP),"@1")])
lambdaIP.con.HC <- c(setNames(1,paste0(endogenous(eLVM.conIP)[1],"~",latent(eLVM.conIP),"@1")),
                   coef(eLVM.conIP)[paste0(endogenous(eLVM.conIP)[-1],"~",latent(eLVM.conIP),"@1")])
lambdaIP.con.D <- c(setNames(1,paste0(endogenous(eLVM.conIP)[1],"~",latent(eLVM.conIP),"@2")),
                  coef(eLVM.conIP)[paste0(endogenous(eLVM.conIP)[-1],"~",latent(eLVM.conIP),"@2")])
lambdaIP.con.DP <- c(setNames(1,paste0(endogenous(eLVM.conIP)[1],"~",latent(eLVM.conIP),"@3")),
                    coef(eLVM.conIP)[paste0(endogenous(eLVM.conIP)[-1],"~",latent(eLVM.conIP),"@3")]) ##neww#
tauIP.mean <- coef(eLVM.meanIP)[paste0(latent(eLVM.meanIP),"~~",latent(eLVM.meanIP),"@1")]
tauIP.con.HC <- coef(eLVM.conIP)[paste0(latent(eLVM.conIP),"~~",latent(eLVM.conIP),"@1")]
tauIP.con.D <- coef(eLVM.conIP)[paste0(latent(eLVM.conIP),"~~",latent(eLVM.conIP),"@2")]
tauIP.con.DP <- coef(eLVM.conIP)[paste0(latent(eLVM.conIP),"~~",latent(eLVM.conIP),"@3")]##new##
range(lambdaIP.mean2 - lambdaIP.mean) ## [1] -1.021931e-06  3.350759e-06


cbind("common" = c(lambdaIP.mean,tauIP.mean),
      "HC"=c(lambdaIP.con.HC,tauIP.con.HC),
      "D"=c(lambdaIP.con.D,tauIP.con.D),
      "DP"=c(lambdaIP.con.DP,tauIP.con.DP))

sigmaCO.mean <- coef(eLVM.meanCO)[paste0(endogenous(eLVM.meanCO),"~~",endogenous(eLVM.meanCO),"@1")]
sigmaCO.con <- coef(eLVM.conCO)[paste0(endogenous(eLVM.conCO),"~~",endogenous(eLVM.conCO),"@1")]
lambdaCO.mean2 <- c(setNames(1,paste0(endogenous(eLVM.blueCO2)[1],"~",latent(eLVM.blueCO2))),
                    coef(eLVM.blueCO2)[paste0(endogenous(eLVM.blueCO2)[-1],"~",latent(eLVM.blueCO2))])
lambdaCO.mean <- c(setNames(1,paste0(endogenous(eLVM.meanCO)[1],"~",latent(eLVM.meanCO),"@1")),
                   coef(eLVM.meanCO)[paste0(endogenous(eLVM.meanCO)[-1],"~",latent(eLVM.meanCO),"@1")])
lambdaCO.con.HC <- c(setNames(1,paste0(endogenous(eLVM.conCO)[1],"~",latent(eLVM.conCO),"@1")),
                     coef(eLVM.conCO)[paste0(endogenous(eLVM.conCO)[-1],"~",latent(eLVM.conCO),"@1")])
lambdaCO.con.D <- c(setNames(1,paste0(endogenous(eLVM.conCO)[1],"~",latent(eLVM.conCO),"@2")),
                    coef(eLVM.conCO)[paste0(endogenous(eLVM.conCO)[-1],"~",latent(eLVM.conCO),"@2")])
lambdaCO.con.DP <- c(setNames(1,paste0(endogenous(eLVM.conCO)[1],"~",latent(eLVM.conCO),"@3")),
                    coef(eLVM.conCO)[paste0(endogenous(eLVM.conCO)[-1],"~",latent(eLVM.conCO),"@3")])

tauCO.mean <- coef(eLVM.meanCO)[paste0(latent(eLVM.meanCO),"~~",latent(eLVM.meanCO),"@1")]
tauCO.con.HC <- coef(eLVM.conCO)[paste0(latent(eLVM.conCO),"~~",latent(eLVM.conCO),"@1")]
tauCO.con.D <- coef(eLVM.conCO)[paste0(latent(eLVM.conCO),"~~",latent(eLVM.conCO),"@2")]
tauCO.con.DP <- coef(eLVM.conCO)[paste0(latent(eLVM.conCO),"~~",latent(eLVM.conCO),"@3")]

cbind("common" = c(lambdaCO.mean,tauCO.mean),
      "HC"=c(lambdaCO.con.HC,tauCO.con.HC),
      "D"=c(lambdaCO.con.D,tauCO.con.D),
      "DP"=c(lambdaCO.con.DP,tauCO.con.DP))


## d) rebuild correlation matrix

## estimate residual variance ##added DP 
resvarIP.con.HC <- lambdaIP.con.HC^2*tauIP.con.HC + sigmaIP.con^2
resvarIP.con.D <- lambdaIP.con.D^2*tauIP.con.D + sigmaIP.con^2
resvarIP.con.DP <- lambdaIP.con.DP^2*tauIP.con.DP + sigmaIP.con^2 ##

resvarCO.con.HC <- lambdaCO.con.HC^2*tauCO.con.HC + sigmaCO.con^2
resvarCO.con.D <- lambdaCO.con.D^2*tauCO.con.D + sigmaCO.con^2
resvarCO.con.DP <- lambdaCO.con.DP^2*tauCO.con.DP + sigmaCO.con^2 ##
## estimate residual correlation
rescovIP.con.HC <- tcrossprod(lambdaIP.con.HC)*tauIP.con.HC + diag(sigmaIP.con^2)
rescovIP.con.D <- tcrossprod(lambdaIP.con.D)*tauIP.con.D + diag(sigmaIP.con^2)
rescovIP.con.DP <- tcrossprod(lambdaIP.con.DP)*tauIP.con.DP + diag(sigmaIP.con^2)
rescorIP.con.HC <- cov2cor(rescovIP.con.HC)
dimnames(rescorIP.con.HC) <- list(endogenous(eLVM.conIP),endogenous(eLVM.conIP))
rescorIP.con.D <- cov2cor(rescovIP.con.D)
dimnames(rescorIP.con.D) <- list(endogenous(eLVM.conIP),endogenous(eLVM.conIP))
rescorIP.con.DP <- cov2cor(rescovIP.con.DP)
dimnames(rescorIP.con.DP) <- list(endogenous(eLVM.conIP),endogenous(eLVM.conIP))


rescovCO.con.HC <- tcrossprod(lambdaCO.con.HC)*tauCO.con.HC + diag(sigmaCO.con^2)
rescovCO.con.D <- tcrossprod(lambdaCO.con.D)*tauCO.con.D + diag(sigmaCO.con^2)
rescovCO.con.DP <- tcrossprod(lambdaCO.con.DP)*tauCO.con.DP + diag(sigmaCO.con^2)
rescorCO.con.HC <- cov2cor(rescovCO.con.HC)
dimnames(rescorCO.con.HC) <- list(endogenous(eLVM.conCO),endogenous(eLVM.conCO))
rescorCO.con.D <- cov2cor(rescovCO.con.D)
dimnames(rescorCO.con.D) <- list(endogenous(eLVM.conCO),endogenous(eLVM.conCO))
rescorCO.con.DP <- cov2cor(rescovCO.con.DP)
dimnames(rescorCO.con.DP) <- list(endogenous(eLVM.conCO),endogenous(eLVM.conCO))
library(ggpubr)

gg.figExtraIP <- ggarrange(ggHeatmap(rescorIP.con.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Control"),
                         ggHeatmap(rescorIP.con.D, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Depressed"),
                         ggHeatmap(rescorIP.con.DP, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Pre-Depressed"),
                         ggHeatmap(rescorIP.con.D-rescorIP.con.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Difference"),
                         common.legend = TRUE, legend = "bottom", nrow = 1)
gg.figExtraIP

gg.figExtraCO <- ggarrange(ggHeatmap(rescorCO.con.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Control"),
                           ggHeatmap(rescorCO.con.D, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Depressed"),
                           ggHeatmap(rescorCO.con.DP, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Pre-Depressed"),
                           ggHeatmap(rescorCO.con.D-rescorCO.con.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Difference"),
                           common.legend = TRUE, legend = "bottom", nrow = 1)
gg.figExtraCO
}
