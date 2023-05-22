### figure1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2023 (17:55) 
## Version: 
## Last-Updated: maj 22 2023 (11:41) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## Patients with a presurgical history of depression were coined TLE-DP,
## whereas those without a presurgical history of depression who developed de novo depression after surgery were coined TLE-DN.
## Finally, patients without any lifetime psychiatric history (pre – and postsurgical) were defined as TLE-HC.
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * 1- load packages
library(ggplot2)
library(data.table)
library(ggpubr)

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
}else{
    path <- ""
}
dfMDD <- readRDS(file.path(path,"data","dfMDD-NNA.rds"))

## * 3- reshape data in the log format for graphical display
## IPSI
dfMDD.long.ipsi <- melt(as.data.table(dfMDD), id.vars = c("ID", "Group"), 
                        measure.vars = c("ipsi.ofc", "ipsi.fusi", "ipsi.insula.c", "ipsi.ros.acc", "ipsi.pos.cc", "ipsi.hippo")
                        )
dfMDD.long.ipsi$hemisphere <- "Ipsilateral hemisphere"
dfMDD.long.ipsi$region <- factor(dfMDD.long.ipsi$variable,
                                 levels = c("ipsi.ofc", "ipsi.fusi", "ipsi.insula.c", "ipsi.ros.acc", "ipsi.pos.cc", "ipsi.hippo"),
                                 labels = c("orbital frontal cortex","fusiform","insula","rostral anterior cingulate","posterior cingulate cortex","hippocampus"))

## CON
dfMDD.long.con <- melt(as.data.table(dfMDD), id.vars = c("ID", "Group"), 
                       measure.vars = c("con.ofc", "con.fusi", "con.insula.c", "con.ros.acc", "con.pos.cc", "con.hippo"))
dfMDD.long.con$hemisphere <- "Contralateral hemisphere"
dfMDD.long.con$region <- factor(dfMDD.long.con$variable,
                                levels = c("con.ofc", "con.fusi", "con.insula.c", "con.ros.acc", "con.pos.cc", "con.hippo"),
                                labels = c("orbital frontal cortex","fusiform","insula","rostral anterior cingulate","posterior cingulate cortex","hippocampus"))

## Assemble
data.figure1 <- rbind(dfMDD.long.ipsi,dfMDD.long.con)
data.figure1$Depression <- factor(data.figure1$Group,
                                  levels = c("MTLE-control","MTLE-MDD-Post","MTLE-MDD-Pre"),
                                  labels = c("TLE-C","TLE-DN","TLE-DP"))


## * 4- graphical display
figure1a <- ggplot(data.figure1[data.figure1$region!="hippocampus"], aes(y=value, x = Depression, fill = Depression)) + geom_boxplot()
figure1a <- figure1a + facet_grid(hemisphere~region)
figure1a <- figure1a + theme(legend.position = "bottom", 
                             text = element_text(size = 13),
                             axis.line = element_line(size = 1.25),
                             axis.ticks = element_line(size = 2),
                             axis.ticks.length=unit(.25, "cm")) 
figure1a <- figure1a + guides(fill = FALSE) + ylab("cortical thickness [mm]") + xlab(NULL)
figure1a

figure1b <- ggplot(data.figure1[data.figure1$region=="hippocampus"], aes(y=value, x = Depression, fill = Depression)) + geom_boxplot()
figure1b <- figure1b + facet_grid(hemisphere~region)
figure1b <- figure1b + theme(legend.position = "bottom", 
                             text = element_text(size = 13),
                             axis.line = element_line(size = 1.25),
                             axis.ticks = element_line(size = 2),
                             axis.ticks.length=unit(.25, "cm")) 
figure1b <- figure1b + guides(fill = "none") + ylab("volume [mm3]") + xlab(NULL)
figure1b

figure1 <- ggarrange(figure1a,figure1b, common.legend = TRUE, widths = c(4,1.15))
figure1

## * 5- export
ggsave(figure1, filename = file.path(path,"figures","figure1.pdf"), width = 14)
ggsave(figure1, filename = file.path(path,"figures","figure1.png"), width = 14)


##----------------------------------------------------------------------
### figure1.R ends here
