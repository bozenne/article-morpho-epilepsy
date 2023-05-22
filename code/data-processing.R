### data-processing.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2023 (17:49) 
## Version: 
## Last-Updated: maj 22 2023 (12:10) 
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
library(readxl)
library(data.table)

## * 2- load data
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){  
    path <- "c:/Users/hpl802/Documents/Github/article-morpho-epilepsy"
    dfMDD0 <- read_excel(file.path(path,"source","LVM _12_01_2022.xlsx"))
    dfFree1 <- read_excel(file.path(path,"source","TLE_MDD_ALL_DATA.xlsx"), sheet = 1)
    dfFreeLH <- read_excel(file.path(path,"source","TLE_MDD_ALL_DATA.xlsx"), sheet = 4)
    dfFreeRH <- read_excel(file.path(path,"source","TLE_MDD_ALL_DATA.xlsx"), sheet = 5)
}else{
    dfMDD0 <- read_excel("/old_data1/Philip_Fink/EXCEL_DATA_SHEETS/LVM-MDD/LVM _12_01_2022.xlsx")
    dfFree1 <- read_excel("/old_data1/Philip_Fink/EXCEL_DATA_SHEETS/TLE_MDD_ALL_DATA.xlsx", sheet = 1)
    dfFreeLH <- read_excel("/old_data1/Philip_Fink/EXCEL_DATA_SHEETS/TLE_MDD_ALL_DATA.xlsx", sheet = 4)
    dfFreeRH <- read_excel("/old_data1/Philip_Fink/EXCEL_DATA_SHEETS/TLE_MDD_ALL_DATA.xlsx", sheet = 5)
}

## * 3- data processing

## merge data files
dfFree <- data.frame(dfFree1[,c("ID","Side")],
                     dfFreeLH[,"lh_insula_thickness",drop=FALSE],
                     dfFreeRH[,"rh_insula_thickness",drop=FALSE])
dfFree$ipsi.insula.c <- NA
dfFree$ipsi.insula.c[dfFree$Side == "right"] <- dfFree$rh_insula_thickness[dfFree$Side == "right"]
dfFree$ipsi.insula.c[dfFree$Side == "left"] <- dfFree$lh_insula_thickness[dfFree$Side == "left"]
    
dfFree$con.insula.c <- NA
dfFree$con.insula.c[dfFree$Side == "right"] <- dfFree$lh_insula_thickness[dfFree$Side == "right"]
dfFree$con.insula.c[dfFree$Side == "left"] <- dfFree$rh_insula_thickness[dfFree$Side == "left"]

dfMDD <- merge(x = as.data.frame(dfMDD0), 
               y = as.data.frame(dfFree)[,c("ID","ipsi.insula.c","con.insula.c")],
               by = c("ID"), all.y = FALSE)

## convert group to factor with shorter names
dfMDD$Group2 <- factor(dfMDD$Group,
                       levels = c( "MTLE-control", "MTLE-MDD-Post", "MTLE-MDD-Pre" ),
                       labels = c("C","MDD_Post","MDD_Pre"))

## rescale volume
dfMDD$ipsi.khippo <- dfMDD$ipsi.hippo/1000
dfMDD$con.khippo <- dfMDD$con.hippo/1000
dfMDD$kETIV <- dfMDD$ETIV/1000

## remove observations with missing data
dfMDD.NNA <- dfMDD[rowSums(is.na(dfMDD))==0,]

dfMDD[rowSums(is.na(dfMDD))>0,]
##          ID Age         Group ipsi.ofc ipsi.prec ipsi.fusi ipsi.mid.temp
## 86 TLE_0244  31 MTLE-MDD-Post 2.688995     2.594     3.024         3.094
## 88 TLE_0248  36 MTLE-MDD-Post 2.535657     2.434     2.698         2.810
##    ipsi.ros.acc ipsi.pos.cc ipsi.sup.fron ipsi.ros.mid.fron ipsi.hippo
## 86        3.256       2.885         2.933             2.479     4460.1
## 88        2.646       2.734         2.823             2.445     4872.1
##    ipsi.lat.ofc ipsi.med.ofc ipsi.cingulate ipsi.acc  con.ofc con.prec con.fusi
## 86        2.775        2.566       2.880359 3.127730 2.613928    2.483    2.947
## 88        2.646        2.364       2.735694 2.811378 2.433714    2.442    2.527
##    con.mid.temp con.ros.acc con.pos.cc con.sup.fron con.ros.mid.fron con.hippo
## 86        3.156       3.119      2.630        2.949            2.337    4364.0
## 88        2.862       2.814      2.603        2.823            2.485    4866.3
##    con.lat.ofc con.med.ofc con.cingulate  con.acc    ETIV ipsi.insula.c
## 86        2.67       2.536      2.606611 2.762682 1278563            NA
## 88        2.52       2.320      2.767578 2.827332 1513632            NA
##    con.insula.c   Group2 ipsi.khippo con.khippo    kETIV
## 86           NA MDD_Post      4.4601     4.3640 1278.563
## 88           NA MDD_Post      4.8721     4.8663 1513.632

## * quality check
name.thickness <- c("ipsi.ofc", "ipsi.fusi", "ipsi.insula.c", "ipsi.ros.acc","ipsi.pos.cc")
any(is.na(dfMDD[,name.thickness]))
## [1] TRUE
range(dfMDD[,name.thickness], na.rm = TRUE)
## [1] 2.000 3.517

name.volume <- c("ipsi.khippo")
any(is.na(dfMDD[,name.volume]))
## [1] FALSE
range(dfMDD[,name.volume])
## [1] 2.4353 5.2023

table(dfMDD$Group2, useNA = "always")
      ##  C MDD-Post  MDD-Pre     <NA> 
      ## 42       24       23        0 

## * Export
if(system("whoami",intern=TRUE)=="unicph\\hpl802"){
    saveRDS(dfMDD, file = file.path(path,"data","dfMDD.rds"))
    saveRDS(dfMDD.NNA, file = file.path(path,"data","dfMDD-NNA.rds"))
}
##----------------------------------------------------------------------
### data-processing.R ends here
