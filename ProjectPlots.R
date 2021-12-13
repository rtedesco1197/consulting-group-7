##Joelle Strom
##STAT 488-002 Consulting
##Visualizations for project report
##Dec. 13, 2021

let <- read.csv("D:/Documents/Applied Stats MS/Fall 2021/STAT 488_002/LetrozolelDat.csv")
let <- let[1:432,]
let$TNK <- as.character(let$TNK)
let <- data.frame(Day=let$d,
                  Tank=factor(let$TNK),
                  Trt=factor(let$TRT, levels=c("control","0.075ug/kg","0.75ug/kg")),
                  Eggs_gFem=let$eggs_gfem,
                  PctFert=let$PctFert,
                  PctVBL=let$PctVBL,
                  Compound=rep("Letrozole", times=nrow(let)))
levels(let$Trt) <- c("Control", "Low Treatment", "High Treatment")

anast <- read.csv("D:/Documents/Applied Stats MS/Fall 2021/STAT 488_002/AnastrozolelDat.csv")
anast <- anast[1:432,]
anast$TNK <- as.character(anast$TNK)
anast <- data.frame(Day=anast$d,
                    Tank=factor(anast$TNK),
                    Trt=factor(anast$TRT, levels=c("control", "0.3mg/kg", "3mg/kg")),
                    Eggs_gFem=anast$Eggs_gFem,
                    PctFert=anast$PctFert,
                    PctVBL=anast$PctVBL,
                    Compound=rep("Anastrozole", times=nrow(anast)))
levels(anast$Trt) <- c("Control", "Low Treatment", "High Treatment")

exe <- read.csv("D:/Documents/Applied Stats MS/Fall 2021/STAT 488_002/ExemestaneDat.csv")
exe <- exe[1:432,]
exe$TNK <- as.character(exe$TNK)
exe <- data.frame(Day=exe$d,
                  Tank=factor(exe$TNK),
                  Trt=factor(exe$TRT, levels=c("control", "0.075ug/kg", "0.75ug/kg")),
                  Eggs_gFem=exe$Eggs_gFem,
                  PctFert=exe$PctFert,
                  PctVBL=exe$PctVBL,
                  Compound=rep("Exemestane", times=nrow(exe)))
levels(exe$Trt) <- c("Control", "Low Treatment", "High Treatment")

fulldat <- rbind(let, anast, exe)
fulldat$Compound <- factor(fulldat$Compound, levels=c("Letrozole", "Anastrozole", "Exemestane"))

library(ggplot2)

EggsPlot <- ggplot(data=fulldat,aes(x=Day,y=Eggs_gFem, group=Trt, colour=Trt)) + stat_summary(geom = "line", fun = mean) + labs(y="# Eggs per g female") + facet_wrap(~Compound)

FertPlot <- ggplot(data=fulldat,aes(x=Day,y=PctFert, group=Trt, colour=Trt)) + stat_summary(geom = "line", fun = mean) + labs(y="% Fertility") + facet_wrap(~Compound)

VBLPlot <- ggplot(data=fulldat,aes(x=Day,y=PctVBL, group=Trt, colour=Trt)) + stat_summary(geom = "line", fun = mean) + labs(y="% Viability") + facet_wrap(~Compound)

png(filename="D:/Documents/Applied Stats MS/Fall 2021/STAT 488_002/EggsPlot.png", height=544, width=788)
EggsPlot
dev.off()

png(filename="D:/Documents/Applied Stats MS/Fall 2021/STAT 488_002/FertPlot.png", height=544, width=788)
FertPlot
dev.off()

png(filename="D:/Documents/Applied Stats MS/Fall 2021/STAT 488_002/VBLPlot.png", height=544, width=788)
VBLPlot
dev.off()