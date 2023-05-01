### table2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2023 (18:37) 
## Version: 
## Last-Updated: maj  1 2023 (18:46) 
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
library(officer)

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
}else{
    path <- ""
}
dfMDD <- readRDS(file.path(path,"data","dfMDD.rds"))
eLVM.blueIP <- readRDS(file.path(path,"data","lvm-blueIP.rds"))
eLVM.blueCO <- readRDS(file.path(path,"data","lvm-blueCO.rds"))

## * 3- prepare data for table
digit.table <- 4

ls.effect <- list(ofc = effects(eLVM.blueIP, ipsi.ofc~Group2MDD),
                  fusi = effects(eLVM.blueIP, ipsi.fusi~Group2MDD),
                  insula = effects(eLVM.blueIP, ipsi.insula.c~Group2MDD),
                  ros.acc = effects(eLVM.blueIP, ipsi.ros.acc~Group2MDD),
                  pos.cc = effects(eLVM.blueIP, ipsi.pos.cc~Group2MDD),
                  khippo = effects(eLVM.blueIP, ipsi.khippo~Group2MDD)
                  )

data.table2 <- do.call(rbind,lapply(ls.effect, function(iL){ ## iL <- ls.effect[[1]]
    data.frame(as.list(c(summary(iL)$coef["Total",], confint(iL)["Total",])))
}))

colnames(data.table2)[colnames(data.table2)=="Pr...z.."] <- "p.value"
colnames(data.table2)[colnames(data.table2)=="X2.5."] <- "lower"
colnames(data.table2)[colnames(data.table2)=="X97.5."] <- "upper"

## * table
table2 <- data.table(
  "Region" = rownames(data.table2),
  "Estimate" = data.table2$Estimate,
  "Unit" = c("mm","mm","mm","mm","mm","??"),
  "95% CI" = paste0("[", round(data.table2$lower, digits = digit.table), " ;", round(data.table2$upper, digits = digit.table),"]"),
  "p-value" = format.pval(data.table2$p.value,digits = digit.table)
)
table2 #### WARNING not adjusted for multiple comparisons over regions (5 comparisons)

## * 5- export
myTable2.doc <- read_docx()
myTable2.doc <- body_add_table(x = myTable2.doc, value = table2)
print(myTable2.doc, target = file.path(path,"tables","table2.docx"))

##----------------------------------------------------------------------
### table2.R ends here
