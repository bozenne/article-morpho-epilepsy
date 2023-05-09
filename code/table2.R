### table2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2023 (18:37) 
## Version: 
## Last-Updated: maj  9 2023 (17:58) 
##           By: Brice Ozenne
##     Update #: 9
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
library(crosstable)

lvmfit2effect <- function(object, name.coef){
    ## calls effects function to output the total effect
    name.endo <- endogenous(object)
    name.latent <- latent(object)

    split.endo <- strsplit(name.endo, split = ".", fixed = TRUE)
    side <- sapply(split.endo,"[",1)
    region <- sapply(split.endo,function(iE){paste(iE[-1],collapse = ".")})
    ls.formula <- lapply(paste(name.endo,"~",name.coef), as.formula)
    
    if(all(grepl(name.coef, names(coef(object)), fixed = TRUE)==FALSE)){
        stop("Could not find \"",name.coef,"\" in the model coefficients. \n")
    }

    ls.effects <- lapply(ls.formula, FUN = function(iF){
        iE <- effects(object, iF)
        data.frame(as.list(c(summary(iE)$coef["Total",], confint(iE)["Total",])))
    })
    df.effects <- do.call(rbind,ls.effects)
    colnames(df.effects)[colnames(df.effects)=="Pr...z.."] <- "p.value"
    colnames(df.effects)[colnames(df.effects)=="X2.5."] <- "lower"
    colnames(df.effects)[colnames(df.effects)=="X97.5."] <- "upper"

    return(cbind(side = side, region = region, df.effects))
}

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
}else{
    path <- ""
}
dfMDD <- readRDS(file.path(path,"data","dfMDD.rds"))
eLVM.blueIP <- readRDS(file.path(path,"data","lvm-blueIP.rds"))
eLVM.blueCO <- readRDS(file.path(path,"data","lvm-blueCO.rds"))

eLVM.blueIP2 <- readRDS(file.path(path,"data","lvm-blueIP2.rds"))
eLVM.blueCO2 <- readRDS(file.path(path,"data","lvm-blueCO2.rds"))

## * 3- prepare data for table
digit.table <- 4

df.region <- rbind(cbind(Group = "MDD_Post", lvmfit2effect(eLVM.blueIP, name.coef = "Group2MDD_Post")),
                   cbind(Group = "MDD_Pre", lvmfit2effect(eLVM.blueIP, name.coef = "Group2MDD_Pre")),
                   cbind(Group = "MDD_Post", lvmfit2effect(eLVM.blueCO, name.coef = "Group2MDD_Post")),
                   cbind(Group = "MDD_Pre", lvmfit2effect(eLVM.blueCO, name.coef = "Group2MDD_Pre")))

df.region2 <- rbind(cbind(Group = "MDD_Post", lvmfit2effect(eLVM.blueIP2, name.coef = "Group2MDD_Post")),
                    cbind(Group = "MDD_Pre", lvmfit2effect(eLVM.blueIP2, name.coef = "Group2MDD_Pre")),
                    cbind(Group = "MDD_Post", lvmfit2effect(eLVM.blueCO2, name.coef = "Group2MDD_Post")),
                    cbind(Group = "MDD_Pre", lvmfit2effect(eLVM.blueCO2, name.coef = "Group2MDD_Pre")))
eLVM.blueIP2


## * table
table2 <- data.table(
  "Group" = df.region$Group,
  "Side" = df.region$side,
  "Region" = df.region$region,
  "Estimate" = round(df.region$Estimate, digits = digit.table),
  "Unit" = c("mm","mm","mm","mm","mm","mm3"),
  "95% CI" = paste0("[", round(df.region$lower, digits = digit.table), " ;", round(df.region$upper, digits = digit.table),"]"),
  "p-value" = format.pval(df.region$p.value,digits = digit.table)
)
table2 #### WARNING not adjusted for multiple comparisons over regions (6 comparisons)
##        Group Side   Region Estimate Unit             95% CI   p-value
##  1: MDD_Post ipsi      ofc  -0.0660   mm [-0.1077 ;-0.0242] 0.0019480
##  2: MDD_Post ipsi     fusi  -0.0991   mm [-0.1645 ;-0.0338] 0.0029497
##  3: MDD_Post ipsi insula.c  -0.0978   mm [-0.1572 ;-0.0384] 0.0012507
##  4: MDD_Post ipsi  ros.acc  -0.1270   mm [-0.2107 ;-0.0432] 0.0029574
##  5: MDD_Post ipsi   pos.cc  -0.0838   mm [-0.1409 ;-0.0267] 0.0040183
##  6: MDD_Post ipsi   khippo  -0.0886  mm3  [-0.2439 ;0.0667] 0.2634554
##  7:  MDD_Pre ipsi      ofc  -0.0431   mm [-0.0814 ;-0.0047] 0.0276365
##  8:  MDD_Pre ipsi     fusi  -0.0647   mm [-0.1236 ;-0.0059] 0.0310773
##  9:  MDD_Pre ipsi insula.c  -0.0639   mm [-0.1196 ;-0.0082] 0.0246174
## 10:  MDD_Pre ipsi  ros.acc  -0.0829   mm [-0.1583 ;-0.0075] 0.0311012
## 11:  MDD_Pre ipsi   pos.cc  -0.0547   mm [-0.1053 ;-0.0041] 0.0341307
## 12:  MDD_Pre ipsi   khippo  -0.0579  mm3  [-0.1655 ;0.0498] 0.2922748
## 13: MDD_Post  con      ofc  -0.0985   mm [-0.1456 ;-0.0514] 4.134e-05
## 14: MDD_Post  con     fusi  -0.1033   mm  [-0.157 ;-0.0496] 0.0001630
## 15: MDD_Post  con insula.c  -0.1314   mm [-0.1982 ;-0.0647] 0.0001135
## 16: MDD_Post  con  ros.acc  -0.1355   mm  [-0.2121 ;-0.059] 0.0005225
## 17: MDD_Post  con   pos.cc  -0.0342   mm  [-0.0819 ;0.0136] 0.1605361
## 18: MDD_Post  con   khippo  -0.0882  mm3   [-0.1894 ;0.013] 0.0875774
## 19:  MDD_Pre  con      ofc  -0.0463   mm [-0.0893 ;-0.0033] 0.0348571
## 20:  MDD_Pre  con     fusi  -0.0486   mm [-0.0948 ;-0.0024] 0.0393191
## 21:  MDD_Pre  con insula.c  -0.0618   mm [-0.1202 ;-0.0034] 0.0379807
## 22:  MDD_Pre  con  ros.acc  -0.0637   mm  [-0.126 ;-0.0015] 0.0447320
## 23:  MDD_Pre  con   pos.cc  -0.0161   mm  [-0.0419 ;0.0098] 0.2228374
## 24:  MDD_Pre  con   khippo  -0.0415  mm3  [-0.0994 ;0.0165] 0.1605005
##        Group Side   Region Estimate Unit             95% CI   p-value

table22 <- data.table(
  "Group" = df.region2$Group,
  "Side" = df.region2$side,
  "Region" = df.region2$region,
  "Estimate" = round(df.region2$Estimate, digits = digit.table),
  "Unit" = c("mm","mm","mm","mm","mm","mm3"),
  "95% CI" = paste0("[", round(df.region2$lower, digits = digit.table), " ;", round(df.region2$upper, digits = digit.table),"]"),
  "p-value" = format.pval(df.region2$p.value,digits = digit.table)
)
table22 #### WARNING not adjusted for multiple comparisons over regions (6 comparisons)
##        Group Side   Region Estimate Unit             95% CI   p-value
##  1: MDD_Post ipsi      ofc  -0.1144   mm [-0.1614 ;-0.0674] 1.846e-06
##  2: MDD_Post ipsi     fusi  -0.0772   mm [-0.1428 ;-0.0115] 0.0212200
##  3: MDD_Post ipsi insula.c  -0.0762   mm  [-0.138 ;-0.0145] 0.0155611
##  4: MDD_Post ipsi  ros.acc  -0.1022   mm [-0.1885 ;-0.0159] 0.0202498
##  5: MDD_Post ipsi   pos.cc  -0.0653   mm [-0.1219 ;-0.0087] 0.0238301
##  6: MDD_Post ipsi   khippo  -0.0852  mm3  [-0.2157 ;0.0452] 0.2004778
##  7:  MDD_Pre ipsi      ofc  -0.1065   mm [-0.1522 ;-0.0609] 4.864e-06
##  8:  MDD_Pre ipsi     fusi  -0.0368   mm   [-0.097 ;0.0233] 0.2302241
##  9:  MDD_Pre ipsi insula.c  -0.0364   mm   [-0.095 ;0.0223] 0.2243008
## 10:  MDD_Pre ipsi  ros.acc  -0.0488   mm  [-0.1283 ;0.0307] 0.2292660
## 11:  MDD_Pre ipsi   pos.cc  -0.0311   mm    [-0.0823 ;0.02] 0.2327092
## 12:  MDD_Pre ipsi   khippo  -0.0407  mm3  [-0.1248 ;0.0435] 0.3439385
## 13: MDD_Post  con      ofc  -0.1045   mm [-0.1523 ;-0.0566] 1.889e-05
## 14: MDD_Post  con     fusi  -0.1038   mm [-0.1572 ;-0.0503] 0.0001426
## 15: MDD_Post  con insula.c  -0.1275   mm [-0.1928 ;-0.0622] 0.0001294
## 16: MDD_Post  con  ros.acc  -0.1391   mm  [-0.2163 ;-0.062] 0.0004068
## 17: MDD_Post  con   pos.cc  -0.0315   mm  [-0.0796 ;0.0167] 0.2006804
## 18: MDD_Post  con   khippo  -0.0563  mm3  [-0.1583 ;0.0457] 0.2790676
## 19:  MDD_Pre  con      ofc  -0.0493   mm [-0.0935 ;-0.0051] 0.0288483
## 20:  MDD_Pre  con     fusi  -0.0490   mm [-0.0944 ;-0.0036] 0.0345376
## 21:  MDD_Pre  con insula.c  -0.0602   mm [-0.1159 ;-0.0045] 0.0341894
## 22:  MDD_Pre  con  ros.acc  -0.0656   mm  [-0.128 ;-0.0033] 0.0390008
## 23:  MDD_Pre  con   pos.cc  -0.0148   mm  [-0.0403 ;0.0106] 0.2530342
## 24:  MDD_Pre  con   khippo  -0.0266  mm3  [-0.0789 ;0.0257] 0.3192979
##        Group Side   Region Estimate Unit             95% CI   p-value

## * 5- export
keep.cols <- c("Region","Estimate","Unit","95% CI","p-value")

myTable2.doc <- read_docx()

myTable2.doc <- body_add_table(x = myTable2.doc, value = table2[Group=="MDD_Post" & Side == "ipsi",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "De Novo MDD vs. Control (ipsilateral side, a priori LVM)")
body_add_par(value = "\n")
myTable2.doc <- body_add_table(x = myTable2.doc, value = table2[Group=="MDD_Pre" & Side == "ipsi",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "Prevalent MDD vs. Control (ipsilateral side, a priori LVM)")
body_add_par(value = "\n")
myTable2.doc <- body_add_table(x = myTable2.doc, value = table2[Group=="MDD_Post" & Side == "con",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "De Novo MDD vs. Control (contralateral side, a priori LVM)")
body_add_par(value = "\n")
myTable2.doc <- body_add_table(x = myTable2.doc, value = table2[Group=="MDD_Pre" & Side == "con",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "Prevalent MDD vs. Control (contralateral side, a priori LVM)")

myTable2.doc <- body_add_break(myTable2.doc)

myTable2.doc <- body_add_table(x = myTable2.doc, value = table22[Group=="MDD_Post" & Side == "ipsi",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "De Novo MDD vs. Control (ipsilateral side, data-driven LVM)")
body_add_par(value = "\n")
myTable2.doc <- body_add_table(x = myTable2.doc, value = table22[Group=="MDD_Pre" & Side == "ipsi",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "Prevalent MDD vs. Control (ipsilateral side, data-driven LVM)")
body_add_par(value = "\n")
myTable2.doc <- body_add_table(x = myTable2.doc, value = table22[Group=="MDD_Post" & Side == "con",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "De Novo MDD vs. Control (contralateral side, data-driven LVM)")
body_add_par(value = "\n")
myTable2.doc <- body_add_table(x = myTable2.doc, value = table22[Group=="MDD_Pre" & Side == "con",.SD,.SDcols = keep.cols])
myTable2.doc <- body_add_table_legend(myTable2.doc, "Prevalent MDD vs. Control (contralateral side, data-driven LVM)")


print(myTable2.doc, target = file.path(path,"tables","table2.docx"))

##----------------------------------------------------------------------
### table2.R ends here
