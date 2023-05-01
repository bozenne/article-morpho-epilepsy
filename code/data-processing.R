### data-processing.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  1 2023 (17:49) 
## Version: 
## Last-Updated: maj  1 2023 (17:54) 
##           By: Brice Ozenne
##     Update #: 4
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
                       labels = c("C","MDD-Post","MDD-Pre"))

## rescale volume
dfMDD$ipsi.khippo <- dfMDD$ipsi.hippo/1000
dfMDD$con.khippo <- dfMDD$con.hippo/1000
dfMDD$kETIV <- dfMDD$ETIV/1000

## quality check
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
}
##----------------------------------------------------------------------
### data-processing.R ends here
