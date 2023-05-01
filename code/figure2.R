### figure2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2023 (18:14) 
## Version: 
## Last-Updated: maj  1 2023 (18:45) 
##           By: Brice Ozenne
##     Update #: 3
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

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
}else{
    path <- ""
}
dfMDD <- readRDS(file.path(path,"data","dfMDD.rds"))
eLVM.blueIP <- readRDS(file.path(path,"data","lvm-blueIP.rds"))
eLVM.blueCO <- readRDS(file.path(path,"data","lvm-blueCO.rds"))

## * 3- prepare data for graphical display

## evaluate partial residuals
adj.dfMDD.ipsi <- data.table(ID = dfMDD$ID, Group = dfMDD$Group,
                             ipsi.ofc = dfMDD$ipsi.ofc - coef(eLVM.blueIP)["ipsi.ofc~Age"]*dfMDD$Age,
                             ipsi.fusi = dfMDD$ipsi.fusi - coef(eLVM.blueIP)["ipsi.fusi~Age"]*dfMDD$Age,
                             ipsi.insula.c = dfMDD$ipsi.insula.c - coef(eLVM.blueIP)["ipsi.insula.c~Age"]*dfMDD$Age,
                             ipsi.ros.acc = dfMDD$ipsi.ros.acc - coef(eLVM.blueIP)["ipsi.ros.acc~Age"]*dfMDD$Age,
                             ipsi.pos.cc = dfMDD$ipsi.pos.cc - coef(eLVM.blueIP)["ipsi.pos.cc~Age"]*dfMDD$Age,
                             ipsi.khippo = dfMDD$ipsi.khippo - coef(eLVM.blueIP)["ipsi.khippo~Age"]*dfMDD$Age - coef(eLVM.blueIP)["ipsi.khippo~kETIV"]*dfMDD$kETIV,
                             ipsi.lv = predict(eLVM.blueIP, x = manifest(eLVM.blueIP), y = latent(eLVM.blueIP))[,"eta"])

adj.dfMDD.long.ipsi <- melt(as.data.table(adj.dfMDD.ipsi), id.vars = c("ID", "Group"), 
                            measure.vars = c("ipsi.ofc","ipsi.fusi","ipsi.insula.c","ipsi.ros.acc","ipsi.pos.cc","ipsi.khippo","ipsi.lv"))
adj.dfMDD.long.ipsi$hemisphere <- "Ipsilateral hemisphere"
adj.dfMDD.long.ipsi$region <- factor(adj.dfMDD.long.ipsi$variable, 
                                     levels = c("ipsi.ofc","ipsi.fusi","ipsi.insula.c","ipsi.ros.acc","ipsi.pos.cc","ipsi.khippo","ipsi.lv"),
                                     labels = c("orbital frontal cortex","fusiform","insula","rostral anterior cingulate","posterior cingulate cortex","hippocampus","latent variable"))

data.figure2 <- adj.dfMDD.long.ipsi
data.figure2$Depression <- factor(data.figure2$Group, levels = c("MTLE-control","MTLE-MDD-Post","MTLE-MDD-Pre"), labels = c("Control","De Novo\nMDD","Prevalent\nMDD"))

## * 4- graphical display
figure2a <- ggplot(data.figure2[data.figure2$region %in% c("hippocampus","latent variable")==FALSE], aes(y=value, x = Depression, fill = Depression)) + geom_boxplot()
figure2a <- figure2a + facet_grid(hemisphere~region)
figure2a <- figure2a + theme(legend.position = "bottom", 
                             text = element_text(size = 13),
                             axis.line = element_line(size = 1.25),
                             axis.ticks = element_line(size = 2),
                             axis.ticks.length=unit(.25, "cm")) 
figure2a <- figure2a + guides(fill = "none") + ylab("cortical thickness [mm]") + xlab(NULL)

figure2b <- ggplot(data.figure2[data.figure2$region=="hippocampus"], aes(y=value, x = Depression, fill = Depression)) + geom_boxplot()
figure2b <- figure2b + facet_grid(hemisphere~region)
figure2b <- figure2b + theme(legend.position = "bottom", 
                             text = element_text(size = 13),
                             axis.line = element_line(size = 1.25),
                             axis.ticks = element_line(size = 2),
                             axis.ticks.length=unit(.25, "cm")) 
figure2b <- figure2b + guides(fill = "none") + ylab("volume [mm3 TOFIX]") + xlab(NULL)

figure2c <- ggplot(data.figure2[data.figure2$region=="latent variable"], aes(y=value, x = Depression, fill = Depression)) + geom_boxplot()
figure2c <- figure2c + facet_grid(hemisphere~region)
figure2c <- figure2c + theme(legend.position = "bottom", 
                             text = element_text(size = 13),
                             axis.line = element_line(size = 1.25),
                             axis.ticks = element_line(size = 2),
                             axis.ticks.length=unit(.25, "cm")) 
figure2c <- figure2c + guides(fill = "none") + ylab("") + xlab(NULL)

figure2 <- ggarrange(figure2a,figure2b,figure2c, common.legend = TRUE, widths = c(4,1,1), nrow = 1)
figure2

## * 5- export
ggsave(figure2, filename = file.path(path,"figures","figure2.pdf"), width = 15)
ggsave(figure2, filename = file.path(path,"figures","figure2.png"), width = 15)


##----------------------------------------------------------------------
### figure2.R ends here
