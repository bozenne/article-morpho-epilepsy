### figure3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  9 2023 (18:57) 
## Version: 
## Last-Updated: maj  9 2023 (19:08) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * 1- load packages
library(ggplot2)
library(data.table)
library(lava)
library(ggpubr)
library(butils)

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
}else{
    path <- ""
}
eLVM.meanIP <- readRDS(file.path(path,"data","lvm-meanIP.rds"))
eLVM.meanCO <- readRDS(file.path(path,"data","lvm-meanCO.rds"))
eLVM.strataIP <- readRDS(file.path(path,"data","lvm-strataIP.rds"))
eLVM.strataCO <- readRDS(file.path(path,"data","lvm-strataCO.rds"))

## * 3- rebuild correlation matrix
## ** extract coefficients

sigmaIP.mean <- coef(eLVM.meanIP)[paste0(endogenous(eLVM.meanIP),"~~",endogenous(eLVM.meanIP),"@1")]
sigmaIP.strata <- coef(eLVM.strataIP)[paste0(endogenous(eLVM.strataIP),"~~",endogenous(eLVM.strataIP),"@1")]

lambdaIP.mean <- c(setNames(1,paste0(endogenous(eLVM.meanIP)[1],"~",latent(eLVM.meanIP),"@1")),
                   coef(eLVM.meanIP)[paste0(endogenous(eLVM.meanIP)[-1],"~",latent(eLVM.meanIP),"@1")])
lambdaIP.strata.HC <- c(setNames(1,paste0(endogenous(eLVM.strataIP)[1],"~",latent(eLVM.strataIP),"@1")),
                        coef(eLVM.strataIP)[paste0(endogenous(eLVM.strataIP)[-1],"~",latent(eLVM.strataIP),"@1")])
lambdaIP.strata.D <- c(setNames(1,paste0(endogenous(eLVM.strataIP)[1],"~",latent(eLVM.strataIP),"@2")),
                       coef(eLVM.strataIP)[paste0(endogenous(eLVM.strataIP)[-1],"~",latent(eLVM.strataIP),"@2")])
lambdaIP.strata.DP <- c(setNames(1,paste0(endogenous(eLVM.strataIP)[1],"~",latent(eLVM.strataIP),"@3")),
                        coef(eLVM.strataIP)[paste0(endogenous(eLVM.strataIP)[-1],"~",latent(eLVM.strataIP),"@3")]) 

tauIP.mean <- coef(eLVM.meanIP)[paste0(latent(eLVM.meanIP),"~~",latent(eLVM.meanIP),"@1")]
tauIP.strata.HC <- coef(eLVM.strataIP)[paste0(latent(eLVM.strataIP),"~~",latent(eLVM.strataIP),"@1")]
tauIP.strata.D <- coef(eLVM.strataIP)[paste0(latent(eLVM.strataIP),"~~",latent(eLVM.strataIP),"@2")]
tauIP.strata.DP <- coef(eLVM.strataIP)[paste0(latent(eLVM.strataIP),"~~",latent(eLVM.strataIP),"@3")]

## cbind("common" = c(lambdaIP.mean,tauIP.mean),
##       "HC"=c(lambdaIP.strata.HC,tauIP.strata.HC),
##       "D"=c(lambdaIP.strata.D,tauIP.strata.D),
##       "DP"=c(lambdaIP.strata.DP,tauIP.strata.DP))

sigmaCO.mean <- coef(eLVM.meanCO)[paste0(endogenous(eLVM.meanCO),"~~",endogenous(eLVM.meanCO),"@1")]
sigmaCO.strata <- coef(eLVM.strataCO)[paste0(endogenous(eLVM.strataCO),"~~",endogenous(eLVM.strataCO),"@1")]

lambdaCO.mean <- c(setNames(1,paste0(endogenous(eLVM.meanCO)[1],"~",latent(eLVM.meanCO),"@1")),
                   coef(eLVM.meanCO)[paste0(endogenous(eLVM.meanCO)[-1],"~",latent(eLVM.meanCO),"@1")])
lambdaCO.strata.HC <- c(setNames(1,paste0(endogenous(eLVM.strataCO)[1],"~",latent(eLVM.strataCO),"@1")),
                     coef(eLVM.strataCO)[paste0(endogenous(eLVM.strataCO)[-1],"~",latent(eLVM.strataCO),"@1")])
lambdaCO.strata.D <- c(setNames(1,paste0(endogenous(eLVM.strataCO)[1],"~",latent(eLVM.strataCO),"@2")),
                    coef(eLVM.strataCO)[paste0(endogenous(eLVM.strataCO)[-1],"~",latent(eLVM.strataCO),"@2")])
lambdaCO.strata.DP <- c(setNames(1,paste0(endogenous(eLVM.strataCO)[1],"~",latent(eLVM.strataCO),"@3")),
                    coef(eLVM.strataCO)[paste0(endogenous(eLVM.strataCO)[-1],"~",latent(eLVM.strataCO),"@3")])

tauCO.mean <- coef(eLVM.meanCO)[paste0(latent(eLVM.meanCO),"~~",latent(eLVM.meanCO),"@1")]
tauCO.strata.HC <- coef(eLVM.strataCO)[paste0(latent(eLVM.strataCO),"~~",latent(eLVM.strataCO),"@1")]
tauCO.strata.D <- coef(eLVM.strataCO)[paste0(latent(eLVM.strataCO),"~~",latent(eLVM.strataCO),"@2")]
tauCO.strata.DP <- coef(eLVM.strataCO)[paste0(latent(eLVM.strataCO),"~~",latent(eLVM.strataCO),"@3")]

## cbind("common" = c(lambdaCO.mean,tauCO.mean),
##       "HC"=c(lambdaCO.con.HC,tauCO.con.HC),
##       "D"=c(lambdaCO.con.D,tauCO.con.D),
##       "DP"=c(lambdaCO.con.DP,tauCO.con.DP))

## ** coefficients to covariance

## resvarIP.strata.HC <- lambdaIP.strata.HC^2*tauIP.strata.HC + sigmaIP.strata^2
## resvarIP.strata.D <- lambdaIP.strata.D^2*tauIP.strata.D + sigmaIP.strata^2
## resvarIP.strata.DP <- lambdaIP.strata.DP^2*tauIP.strata.DP + sigmaIP.strata^2 ##

## resvarCO.strata.HC <- lambdaCO.strata.HC^2*tauCO.strata.HC + sigmaCO.strata^2
## resvarCO.strata.D <- lambdaCO.strata.D^2*tauCO.strata.D + sigmaCO.strata^2
## resvarCO.strata.DP <- lambdaCO.strata.DP^2*tauCO.strata.DP + sigmaCO.strata^2 ##

rescovIP.strata.HC <- tcrossprod(lambdaIP.strata.HC)*tauIP.strata.HC + diag(sigmaIP.strata^2)
rescovIP.strata.D <- tcrossprod(lambdaIP.strata.D)*tauIP.strata.D + diag(sigmaIP.strata^2)
rescovIP.strata.DP <- tcrossprod(lambdaIP.strata.DP)*tauIP.strata.DP + diag(sigmaIP.strata^2)

rescovCO.strata.HC <- tcrossprod(lambdaCO.strata.HC)*tauCO.strata.HC + diag(sigmaCO.strata^2)
rescovCO.strata.D <- tcrossprod(lambdaCO.strata.D)*tauCO.strata.D + diag(sigmaCO.strata^2)
rescovCO.strata.DP <- tcrossprod(lambdaCO.strata.DP)*tauCO.strata.DP + diag(sigmaCO.strata^2)

## ** covariance to correlation
rescorIP.strata.HC <- cov2cor(rescovIP.strata.HC)
dimnames(rescorIP.strata.HC) <- list(endogenous(eLVM.strataIP),endogenous(eLVM.strataIP))

rescorIP.strata.D <- cov2cor(rescovIP.strata.D)
dimnames(rescorIP.strata.D) <- list(endogenous(eLVM.strataIP),endogenous(eLVM.strataIP))

rescorIP.strata.DP <- cov2cor(rescovIP.strata.DP)
dimnames(rescorIP.strata.DP) <- list(endogenous(eLVM.strataIP),endogenous(eLVM.strataIP))

rescorCO.strata.HC <- cov2cor(rescovCO.strata.HC)
dimnames(rescorCO.strata.HC) <- list(endogenous(eLVM.strataCO),endogenous(eLVM.strataCO))

rescorCO.strata.D <- cov2cor(rescovCO.strata.D)
dimnames(rescorCO.strata.D) <- list(endogenous(eLVM.strataCO),endogenous(eLVM.strataCO))

rescorCO.strata.DP <- cov2cor(rescovCO.strata.DP)
dimnames(rescorCO.strata.DP) <- list(endogenous(eLVM.strataCO),endogenous(eLVM.strataCO))

## * 4- graphical display
figure3.1 <- ggarrange(ggHeatmap(rescorIP.strata.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1), legend_title = "correlation\n (ipsilateral)")$plot + ggtitle("Control"),
                       ggHeatmap(rescorIP.strata.D, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("De novo MDD"),
                       ggHeatmap(rescorIP.strata.DP, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Prevalent MDD"),
                       NULL,
                       ggHeatmap(rescorIP.strata.D-rescorIP.strata.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("De novo MDD vs. control"),
                       ggHeatmap(rescorIP.strata.DP-rescorIP.strata.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Prevalent MDD vs. control"),
                       common.legend = TRUE, legend = "bottom", nrow = 2, ncol = 3)
## figure3.1

figure3.2 <- ggarrange(ggHeatmap(rescorCO.strata.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1), legend_title = "correlation\n (contralateral)")$plot + ggtitle("Control"),
                       ggHeatmap(rescorCO.strata.D, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("De novo MDD"),
                       ggHeatmap(rescorCO.strata.DP, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Prevalent MDD"),
                       NULL,
                       ggHeatmap(rescorCO.strata.D-rescorCO.strata.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("De novo MDD vs. control"),
                       ggHeatmap(rescorCO.strata.DP-rescorCO.strata.HC, plot = FALSE, add.text = "fill", round = 2, limits = c(-1,1))$plot + ggtitle("Prevalent MDD vs. control"),
                       common.legend = TRUE, legend = "bottom", nrow = 2, ncol = 3)
## figure3.2
 
## * 5- export
ggsave(figure3.1, filename = file.path(path,"figures","figure3.1.pdf"), width = 15, height = 10)
ggsave(figure3.2, filename = file.path(path,"figures","figure3.1.png"), width = 15, height = 10)

ggsave(figure3.1, filename = file.path(path,"figures","figure3.2.pdf"), width = 15, height = 10)
ggsave(figure3.2, filename = file.path(path,"figures","figure3.2.png"), width = 15, height = 10)

##----------------------------------------------------------------------
### figure3.R ends here
