########################################################
# analysis.R
# This file will run the log-linear analysis for all models.
# Many different models are run here, but we use a certain syntax
# to identify them. 
#
# First, any model that has a WB in the name is run on
# the collapsed data where browns and blacks are treated as 
# one group. 
#
# Second, each model has a number that indicates how we treated
# the issue of different patterns of HE*WE by race, as discussed
# in the article and the appendix. These are:
# (1) treat interracial couples like the endogamous lighter group
# (2) treat interracial couples like the endogamous darker group
# (3) use the geometric mean of endogamous groups
# (4) just use a pooled HE*WE - this is the model we report in the paper.
# (5) Hou and Myles approach, with light as ref
# (6) Hou and Myles approach, with darks as ref
#
# Third, there are four different model types:
# basemodel - baseline model with no exchange parameters
# semodel - just dyadic exchange terms
# ebomodel - just market exchange terms
# ebmodel - both dyadic and market exchange terms
########################################################

#load the data
data <- read.csv("output/table_marriages.csv")
dataWB <- read.csv("output/tableWB_marriages.csv")

###Calculate overall odds ratios of marriage###
tableWB <- tapply(dataWB$Freq,dataWB[,c("HR","WR")],sum)
lOR.wb.collapse <- log(tableWB["W","Bl"]*tableWB["Bl","W"]/(tableWB["W","W"]*tableWB["Bl","Bl"]))

table <- tapply(data$Freq,data[,c("HR","WR")],sum)
lOR.wb <- log(table["Bl","W"]*table["W","Bl"]/(table["W","W"]*table["Bl","Bl"]))
lOR.wr <- log(table["Br","W"]*table["W","Br"]/(table["W","W"]*table["Br","Br"]))
lOR.br <- log(table["Br","Bl"]*table["Bl","Br"]/(table["Br","Br"]*table["Bl","Bl"]))

####White/Black collapsed first###

#code SE terms
dataWB$SEBMWFUp <- dataWB$HE>dataWB$WE & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$SEBMWFDown <- dataWB$HE<dataWB$WE & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$SEWMBFUp <- dataWB$WE>dataWB$HE & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$SEWMBFDown <- dataWB$WE<dataWB$HE & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$SEBMBFUp <- dataWB$HE>dataWB$WE & dataWB$HR=="Bl" & dataWB$WR=="Bl"
dataWB$SEBMBFDown <- dataWB$HE<dataWB$WE & dataWB$HR=="Bl" & dataWB$WR=="Bl"
dataWB$SEWMWFUp <- dataWB$HE>dataWB$WE & dataWB$HR=="W" & dataWB$WR=="W"
dataWB$SEWMWFDown <- dataWB$HE<dataWB$WE & dataWB$HR=="W" & dataWB$WR=="W"

#symmetric terms
dataWB$SEBMWF <- 0
dataWB$SEBMWF[dataWB$HE>dataWB$WE & dataWB$HR=="Bl" & dataWB$WR=="W"] <- 1
dataWB$SEBMWF[dataWB$HE<dataWB$WE & dataWB$HR=="Bl" & dataWB$WR=="W"] <- -1
dataWB$SEWMBF <- 0
dataWB$SEWMBF[dataWB$HE<dataWB$WE & dataWB$HR=="W" & dataWB$WR=="Bl"] <- 1
dataWB$SEWMBF[dataWB$HE>dataWB$WE & dataWB$HR=="W" & dataWB$WR=="Bl"] <- -1
dataWB$SE <- dataWB$SEBMWF+dataWB$SEWMBF

#code boundary terms
dataWB$EBBM1 <- dataWB$HE>1 & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$EBBM2 <- dataWB$HE>2 & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$EBBM3 <- dataWB$HE>3 & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$EBBM4 <- dataWB$HE>4 & dataWB$HR=="Bl" & dataWB$WR=="W"

dataWB$EBWF1 <- dataWB$WE>1 & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$EBWF2 <- dataWB$WE>2 & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$EBWF3 <- dataWB$WE>3 & dataWB$HR=="Bl" & dataWB$WR=="W"
dataWB$EBWF4 <- dataWB$WE>4 & dataWB$HR=="Bl" & dataWB$WR=="W"

dataWB$EBWM1 <- dataWB$HE>1 & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$EBWM2 <- dataWB$HE>2 & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$EBWM3 <- dataWB$HE>3 & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$EBWM4 <- dataWB$HE>4 & dataWB$HR=="W" & dataWB$WR=="Bl"

dataWB$EBBF1 <- dataWB$WE>1 & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$EBBF2 <- dataWB$WE>2 & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$EBBF3 <- dataWB$WE>3 & dataWB$HR=="W" & dataWB$WR=="Bl"
dataWB$EBBF4 <- dataWB$WE>4 & dataWB$HR=="W" & dataWB$WR=="Bl"

#Various EAM terms

dataWB$WWEAM1.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==2
dataWB$WWEAM1.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==2
dataWB$WWEAM1.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==2)
dataWB$BBEAM1.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==2
dataWB$BBEAM1.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==2
dataWB$BBEAM1.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==2)

dataWB$WWEAM2.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==3
dataWB$WWEAM2.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==3
dataWB$WWEAM2.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==3)
dataWB$BBEAM2.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==3
dataWB$BBEAM2.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==3
dataWB$BBEAM2.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==3)

dataWB$WWEAM3.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==4
dataWB$WWEAM3.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==4
dataWB$WWEAM3.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==4)
dataWB$BBEAM3.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==4
dataWB$BBEAM3.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==4
dataWB$BBEAM3.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==4)

dataWB$WWEAM4.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==5
dataWB$WWEAM4.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==5
dataWB$WWEAM4.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==2 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==5)
dataWB$BBEAM4.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==5
dataWB$BBEAM4.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==5
dataWB$BBEAM4.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==2 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==2 & dataWB$WE==5)

dataWB$WWEAM5.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==2
dataWB$WWEAM5.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==2
dataWB$WWEAM5.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==2)
dataWB$BBEAM5.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==2
dataWB$BBEAM5.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==2
dataWB$BBEAM5.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==2)

dataWB$WWEAM6.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==3
dataWB$WWEAM6.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==3
dataWB$WWEAM6.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==3)
dataWB$BBEAM6.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==3
dataWB$BBEAM6.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==3
dataWB$BBEAM6.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==3)

dataWB$WWEAM7.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==4
dataWB$WWEAM7.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==4
dataWB$WWEAM7.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==4)
dataWB$BBEAM7.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==4
dataWB$BBEAM7.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==4
dataWB$BBEAM7.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==4)

dataWB$WWEAM8.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==5
dataWB$WWEAM8.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==5
dataWB$WWEAM8.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==3 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==5)
dataWB$BBEAM8.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==5
dataWB$BBEAM8.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==5
dataWB$BBEAM8.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==3 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==3 & dataWB$WE==5)

dataWB$WWEAM9.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==2
dataWB$WWEAM9.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==2
dataWB$WWEAM9.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==2)
dataWB$BBEAM9.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==2
dataWB$BBEAM9.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==2
dataWB$BBEAM9.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==2)

dataWB$WWEAM10.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==3
dataWB$WWEAM10.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==3
dataWB$WWEAM10.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==3)
dataWB$BBEAM10.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==3
dataWB$BBEAM10.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==3
dataWB$BBEAM10.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==3)

dataWB$WWEAM11.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==4
dataWB$WWEAM11.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==4
dataWB$WWEAM11.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==4)
dataWB$BBEAM11.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==4
dataWB$BBEAM11.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==4
dataWB$BBEAM11.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==4)

dataWB$WWEAM12.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==5
dataWB$WWEAM12.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==5
dataWB$WWEAM12.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==4 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==5)
dataWB$BBEAM12.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==5
dataWB$BBEAM12.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==5
dataWB$BBEAM12.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==4 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==4 & dataWB$WE==5)

dataWB$WWEAM13.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==2
dataWB$WWEAM13.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==2
dataWB$WWEAM13.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==2)
dataWB$BBEAM13.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==2
dataWB$BBEAM13.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==2
dataWB$BBEAM13.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==2) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==2)

dataWB$WWEAM14.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==3
dataWB$WWEAM14.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==3
dataWB$WWEAM14.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==3)
dataWB$BBEAM14.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==3
dataWB$BBEAM14.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==3
dataWB$BBEAM14.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==3) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==3)

dataWB$WWEAM15.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==4
dataWB$WWEAM15.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==4
dataWB$WWEAM15.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==4)
dataWB$BBEAM15.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==4
dataWB$BBEAM15.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==4
dataWB$BBEAM15.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==4) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==4)

dataWB$WWEAM16.1 <- ((dataWB$HR=="W" & dataWB$WR=="W") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==5
dataWB$WWEAM16.2 <- dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==5
dataWB$WWEAM16.3 <- (dataWB$HR=="W" & dataWB$WR=="W" & dataWB$HE==5 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==5)
dataWB$BBEAM16.1 <- dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==5
dataWB$BBEAM16.2 <- ((dataWB$HR=="Bl" & dataWB$WR=="Bl") | (dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==5
dataWB$BBEAM16.3 <- (dataWB$HR=="Bl" & dataWB$WR=="Bl" & dataWB$HE==5 & dataWB$WE==5) + 0.5 * (((dataWB$HR=="W" & dataWB$WR=="Bl") | (dataWB$HR=="Bl" & dataWB$WR=="W"))  & dataWB$HE==5 & dataWB$WE==5)

dataWB2 <- dataWB 

dataWB2$HE <- as.factor(dataWB2$HE)
dataWB2$WE <- as.factor(dataWB2$WE)

##Run WB models

basemodel.wb <- glm(I(round(Freq,0))~HR*WR+HE*HR+WE*WR, family=poisson, data=dataWB2)
basemodel1.wb <- update(basemodel.wb,.~.+WWEAM1.1+WWEAM2.1+WWEAM3.1+WWEAM4.1+WWEAM5.1+WWEAM6.1+WWEAM7.1+WWEAM8.1+
									+WWEAM9.1+WWEAM10.1+WWEAM11.1+WWEAM12.1+WWEAM13.1+WWEAM14.1+WWEAM15.1+WWEAM16.1
									+BBEAM1.1+BBEAM2.1+BBEAM3.1+BBEAM4.1+BBEAM5.1+BBEAM6.1+BBEAM7.1+BBEAM8.1+
									+BBEAM9.1+BBEAM10.1+BBEAM11.1+BBEAM12.1+BBEAM13.1+BBEAM14.1+BBEAM15.1+BBEAM16.1)
semodel1.wb <- update(basemodel1.wb, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown)
ebmodel1.wb <- update(semodel1.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)
ebomodel1.wb <- update(basemodel1.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)

basemodel2.wb <- update(basemodel.wb,.~.+WWEAM1.2+WWEAM2.2+WWEAM3.2+WWEAM4.2+WWEAM5.2+WWEAM6.2+WWEAM7.2+WWEAM8.2+
									+WWEAM9.2+WWEAM10.2+WWEAM11.2+WWEAM12.2+WWEAM13.2+WWEAM14.2+WWEAM15.2+WWEAM16.2
									+BBEAM1.2+BBEAM2.2+BBEAM3.2+BBEAM4.2+BBEAM5.2+BBEAM6.2+BBEAM7.2+BBEAM8.2+
									+BBEAM9.2+BBEAM10.2+BBEAM11.2+BBEAM12.2+BBEAM13.2+BBEAM14.2+BBEAM15.2+BBEAM16.2)
semodel2.wb <- update(basemodel2.wb, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown)
ebmodel2.wb <- update(semodel2.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)
ebomodel2.wb <- update(basemodel2.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)

basemodel3.wb <- update(basemodel.wb,.~.+WWEAM1.3+WWEAM2.3+WWEAM3.3+WWEAM4.3+WWEAM5.3+WWEAM6.3+WWEAM7.3+WWEAM8.3+
									+WWEAM9.3+WWEAM10.3+WWEAM11.3+WWEAM12.3+WWEAM13.3+WWEAM14.3+WWEAM15.3+WWEAM16.3
									+BBEAM1.3+BBEAM2.3+BBEAM3.3+BBEAM4.3+BBEAM5.3+BBEAM6.3+BBEAM7.3+BBEAM8.3+
									+BBEAM9.3+BBEAM10.3+BBEAM11.3+BBEAM12.3+BBEAM13.3+BBEAM14.3+BBEAM15.3+BBEAM16.3)
semodel3.wb <- update(basemodel3.wb, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown)
ebmodel3.wb <- update(semodel3.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)
ebomodel3.wb <- update(basemodel3.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)

basemodel4.wb <- update(basemodel.wb,.~.+HE*WE)
semodel4.wb <- update(basemodel4.wb, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown)
ebmodel4.wb <- update(semodel4.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)
ebomodel4.wb <- update(basemodel4.wb, .~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)

semodel5.wb <- update(basemodel4.wb,.~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown+SEBMBFUp+SEBMBFDown)
ebmodel5.wb <- update(semodel5.wb,.~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)
semodel6.wb <- update(basemodel4.wb,.~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown+SEWMWFUp+SEWMWFDown)
ebmodel6.wb <- update(semodel6.wb,.~.+EBBM1+EBBM2+EBBM3+EBBM4+EBWF1+EBWF2+EBWF3+EBWF4
									+EBWM1+EBWM2+EBWM3+EBWM4+EBBF1+EBBF2+EBBF3+EBBF4)

####White/Brown/Black######

# R - brown

#code SE terms
data$SEBMWFUp <- data$HE>data$WE & data$HR=="Bl" & data$WR=="W"
data$SEBMWFDown <- data$HE<data$WE & data$HR=="Bl" & data$WR=="W"
data$SEWMBFUp <- data$WE>data$HE & data$HR=="W" & data$WR=="Bl"
data$SEWMBFDown <- data$WE<data$HE & data$HR=="W" & data$WR=="Bl"

data$SERMWFUp <- data$HE>data$WE & data$HR=="Br" & data$WR=="W"
data$SERMWFDown <- data$HE<data$WE & data$HR=="Br" & data$WR=="W"
data$SEWMRFUp <- data$WE>data$HE & data$HR=="W" & data$WR=="Br"
data$SEWMRFDown <- data$WE<data$HE & data$HR=="W" & data$WR=="Br"

data$SEBMRFUp <- data$HE>data$WE & data$HR=="Bl" & data$WR=="Br"
data$SEBMRFDown <- data$HE<data$WE & data$HR=="Bl" & data$WR=="Br"
data$SERMBFUp <- data$WE>data$HE & data$HR=="Br" & data$WR=="Bl"
data$SERMBFDown <- data$WE<data$HE & data$HR=="Br" & data$WR=="Bl"

#test a single term with scalar
data$SEDMup <- data$SEBMWFUp | data$SERMWFUp | data$SEBMRFUp
data$SEDMdown <- data$SEBMWFDown | data$SERMWFDown | data$SEBMRFDown
data$SEDMup2 <- data$SEBMWFUp | data$SERMWFUp
data$SEDMdown2 <- data$SEBMWFDown | data$SERMWFDown
data$SEDMup3 <- data$SEBMWFUp
data$SEDMdown3 <- data$SEBMWFDown
data$SEDWup <- data$SEWMBFUp | data$SEWMRFUp | data$SERMBFUp
data$SEDWdown <- data$SEWMBFDown | data$SEWMRFDown | data$SERMBFDown
data$SEDWup2 <- data$SEWMBFUp | data$SEWMRFUp
data$SEDWdown2 <- data$SEWMBFDown | data$SEWMRFDown
data$SEDWup3 <- data$SEWMBFUp
data$SEDWdown3 <- data$SEWMBFDown

#scalar
data$SEDMUpScale <- 0
data$SEDMUpScale[data$SEBMRFUp] <- 1
data$SEDMUpScale[data$SERMWFUp] <- 2
data$SEDMUpScale[data$SEBMWFUp] <- 3
data$SEDMDownScale <- 0
data$SEDMDownScale[data$SEBMRFDown] <- 1
data$SEDMDownScale[data$SERMWFDown] <- 2
data$SEDMDownScale[data$SEBMWFDown] <- 3

data$SEDFUpScale <- 0
data$SEDFUpScale[data$SERMBFUp] <- 1
data$SEDFUpScale[data$SEWMRFUp] <- 2
data$SEDFUpScale[data$SEWMBFUp] <- 3
data$SEDFDownScale <- 0
data$SEDFDownScale[data$SERMBFDown] <- 1
data$SEDFDownScale[data$SEWMRFDown] <- 2
data$SEDFDownScale[data$SEWMBFDown] <- 3

data$SEDMUpScale2 <- 0
data$SEDMUpScale2[data$SEBMRFUp] <- 1
data$SEDMUpScale2[data$SERMWFUp] <- 2
data$SEDMUpScale2[data$SEBMWFUp] <- 2
data$SEDMDownScale2 <- 0
data$SEDMDownScale2[data$SEBMRFDown] <- 1
data$SEDMDownScale2[data$SERMWFDown] <- 2
data$SEDMDownScale2[data$SEBMWFDown] <- 2

data$SEDFUpScale2 <- 0
data$SEDFUpScale2[data$SERMBFUp] <- 1
data$SEDFUpScale2[data$SEWMRFUp] <- 2
data$SEDFUpScale2[data$SEWMBFUp] <- 2
data$SEDFDownScale2 <- 0
data$SEDFDownScale2[data$SERMBFDown] <- 1
data$SEDFDownScale2[data$SEWMRFDown] <- 2
data$SEDFDownScale2[data$SEWMBFDown] <- 2

#educational boundaries
data$EBBMW1 <- data$HE>1 & data$HR=="Bl" & data$WR=="W"
data$EBBMW2 <- data$HE>2 & data$HR=="Bl" & data$WR=="W"
data$EBBMW3 <- data$HE>3 & data$HR=="Bl" & data$WR=="W"
data$EBBMW4 <- data$HE>4 & data$HR=="Bl" & data$WR=="W"

data$EBWFB1 <- data$WE>1 & data$HR=="Bl" & data$WR=="W"
data$EBWFB2 <- data$WE>2 & data$HR=="Bl" & data$WR=="W"
data$EBWFB3 <- data$WE>3 & data$HR=="Bl" & data$WR=="W"
data$EBWFB4 <- data$WE>4 & data$HR=="Bl" & data$WR=="W"

data$EBWMB1 <- data$HE>1 & data$HR=="W" & data$WR=="Bl"
data$EBWMB2 <- data$HE>2 & data$HR=="W" & data$WR=="Bl"
data$EBWMB3 <- data$HE>3 & data$HR=="W" & data$WR=="Bl"
data$EBWMB4 <- data$HE>4 & data$HR=="W" & data$WR=="Bl"

data$EBBFW1 <- data$WE>1 & data$HR=="W" & data$WR=="Bl"
data$EBBFW2 <- data$WE>2 & data$HR=="W" & data$WR=="Bl"
data$EBBFW3 <- data$WE>3 & data$HR=="W" & data$WR=="Bl"
data$EBBFW4 <- data$WE>4 & data$HR=="W" & data$WR=="Bl"

data$EBRMW1 <- data$HE>1 & data$HR=="Br" & data$WR=="W"
data$EBRMW2 <- data$HE>2 & data$HR=="Br" & data$WR=="W"
data$EBRMW3 <- data$HE>3 & data$HR=="Br" & data$WR=="W"
data$EBRMW4 <- data$HE>4 & data$HR=="Br" & data$WR=="W"

data$EBWFR1 <- data$WE>1 & data$HR=="Br" & data$WR=="W"
data$EBWFR2 <- data$WE>2 & data$HR=="Br" & data$WR=="W"
data$EBWFR3 <- data$WE>3 & data$HR=="Br" & data$WR=="W"
data$EBWFR4 <- data$WE>4 & data$HR=="Br" & data$WR=="W"

data$EBWMR1 <- data$HE>1 & data$HR=="W" & data$WR=="Br"
data$EBWMR2 <- data$HE>2 & data$HR=="W" & data$WR=="Br"
data$EBWMR3 <- data$HE>3 & data$HR=="W" & data$WR=="Br"
data$EBWMR4 <- data$HE>4 & data$HR=="W" & data$WR=="Br"

data$EBRFW1 <- data$WE>1 & data$HR=="W" & data$WR=="Br"
data$EBRFW2 <- data$WE>2 & data$HR=="W" & data$WR=="Br"
data$EBRFW3 <- data$WE>3 & data$HR=="W" & data$WR=="Br"
data$EBRFW4 <- data$WE>4 & data$HR=="W" & data$WR=="Br"

data$EBBMR1 <- data$HE>1 & data$HR=="Bl" & data$WR=="Br"
data$EBBMR2 <- data$HE>2 & data$HR=="Bl" & data$WR=="Br"
data$EBBMR3 <- data$HE>3 & data$HR=="Bl" & data$WR=="Br"
data$EBBMR4 <- data$HE>4 & data$HR=="Bl" & data$WR=="Br"

data$EBRFB1 <- data$WE>1 & data$HR=="Bl" & data$WR=="Br"
data$EBRFB2 <- data$WE>2 & data$HR=="Bl" & data$WR=="Br"
data$EBRFB3 <- data$WE>3 & data$HR=="Bl" & data$WR=="Br"
data$EBRFB4 <- data$WE>4 & data$HR=="Bl" & data$WR=="Br"

data$EBRMB1 <- data$HE>1 & data$HR=="Br" & data$WR=="Bl"
data$EBRMB2 <- data$HE>2 & data$HR=="Br" & data$WR=="Bl"
data$EBRMB3 <- data$HE>3 & data$HR=="Br" & data$WR=="Bl"
data$EBRMB4 <- data$HE>4 & data$HR=="Br" & data$WR=="Bl"

data$EBBFR1 <- data$WE>1 & data$HR=="Br" & data$WR=="Bl"
data$EBBFR2 <- data$WE>2 & data$HR=="Br" & data$WR=="Bl"
data$EBBFR3 <- data$WE>3 & data$HR=="Br" & data$WR=="Bl"
data$EBBFR4 <- data$WE>4 & data$HR=="Br" & data$WR=="Bl"

#scalar approach (leave out black/brown)

data$EBDM1 <- 0
data$EBDM1[data$EBRMW1] <- 1
data$EBDM1[data$EBBMW1] <- 2
data$EBDM2 <- 0
data$EBDM2[data$EBRMW2] <- 1
data$EBDM2[data$EBBMW2] <- 2
data$EBDM3 <- 0
data$EBDM3[data$EBRMW3] <- 1
data$EBDM3[data$EBBMW3] <- 2
data$EBDM4 <- 0
data$EBDM4[data$EBRMW4] <- 1
data$EBDM4[data$EBBMW4] <- 2
data$EBDF1 <- 0
data$EBDF1[data$EBRFW1] <- 1
data$EBDF1[data$EBBFW1] <- 2
data$EBDF2 <- 0
data$EBDF2[data$EBRFW2] <- 1
data$EBDF2[data$EBBFW2] <- 2
data$EBDF3 <- 0
data$EBDF3[data$EBRFW3] <- 1
data$EBDF3[data$EBBFW3] <- 2
data$EBDF4 <- 0
data$EBDF4[data$EBRFW4] <- 1
data$EBDF4[data$EBBFW4] <- 2

data$EBLM1 <- 0
data$EBLM1[data$EBWMR1] <- 1
data$EBLM1[data$EBWMB1] <- 2
data$EBLM2 <- 0
data$EBLM2[data$EBWMR2] <- 1
data$EBLM2[data$EBWMB2] <- 2
data$EBLM3 <- 0
data$EBLM3[data$EBWMR3] <- 1
data$EBLM3[data$EBWMB3] <- 2
data$EBLM4 <- 0
data$EBLM4[data$EBWMR4] <- 1
data$EBLM4[data$EBWMB4] <- 2
data$EBLF1 <- 0
data$EBLF1[data$EBWFR1] <- 1
data$EBLF1[data$EBWFB1] <- 2
data$EBLF2 <- 0
data$EBLF2[data$EBWFR2] <- 1
data$EBLF2[data$EBWFB2] <- 2
data$EBLF3 <- 0
data$EBLF3[data$EBWFR3] <- 1
data$EBLF3[data$EBWFB3] <- 2
data$EBLF4 <- 0
data$EBLF4[data$EBWFR4] <- 1
data$EBLF4[data$EBWFB4] <- 2

data$EBDM1.2 <- data$EBRMW1 | data$EBBMW1
data$EBDM2.2 <- data$EBRMW2 | data$EBBMW2
data$EBDM3.2 <- data$EBRMW3 | data$EBBMW3
data$EBDM4.2 <- data$EBRMW4 | data$EBBMW4
data$EBDF1.2 <- data$EBRFW1 | data$EBBFW1
data$EBDF2.2 <- data$EBRFW2 | data$EBBFW2
data$EBDF3.2 <- data$EBRFW3 | data$EBBFW3
data$EBDF4.2 <- data$EBRFW4 | data$EBBFW4

data$EBLM1.2 <- data$EBWMB1 | data$EBWMR1
data$EBLM2.2 <- data$EBWMB2 | data$EBWMR2
data$EBLM3.2 <- data$EBWMB3 | data$EBWMR3
data$EBLM4.2 <- data$EBWMB4 | data$EBWMR4
data$EBLF1.2 <- data$EBWFB1 | data$EBWFR1
data$EBLF2.2 <- data$EBWFB2 | data$EBWFR2
data$EBLF3.2 <- data$EBWFB3 | data$EBWFR3
data$EBLF4.2 <- data$EBWFB4 | data$EBWFR4

#accounting for EAM is a little more complicated because we have more options than like white or like black
# I will do this so its like the lighter group or like the darker group

data$WWEAM1.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==2
data$WWEAM1.2 <- data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==2
data$WWEAM1.3 <- (data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==2)
data$RREAM1.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==2 & data$WE==2
data$RREAM1.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==2 & data$WE==2
data$RREAM1.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==2 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==2)
data$BBEAM1.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==2
data$BBEAM1.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==2 & data$WE==2
data$BBEAM1.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==2)

data$WWEAM2.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==3
data$WWEAM2.2 <- data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==3
data$WWEAM2.3 <- (data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==3)
data$RREAM2.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==2 & data$WE==3
data$RREAM2.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==2 & data$WE==3
data$RREAM2.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==2 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==3)
data$BBEAM2.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==3
data$BBEAM2.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==2 & data$WE==3
data$BBEAM2.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==3)

data$WWEAM3.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==4
data$WWEAM3.2 <- data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==4
data$WWEAM3.3 <- (data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==4)
data$RREAM3.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==2 & data$WE==4
data$RREAM3.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==2 & data$WE==4
data$RREAM3.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==2 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==4)
data$BBEAM3.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==4
data$BBEAM3.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==2 & data$WE==4
data$BBEAM3.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==4)

data$WWEAM4.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==5
data$WWEAM4.2 <- data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==5
data$WWEAM4.3 <- (data$HR=="W" & data$WR=="W" & data$HE==2 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==2 & data$WE==5)
data$RREAM4.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==2 & data$WE==5
data$RREAM4.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==2 & data$WE==5
data$RREAM4.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==2 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==5)
data$BBEAM4.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==5
data$BBEAM4.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==2 & data$WE==5
data$BBEAM4.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==2 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==2 & data$WE==5)

data$WWEAM5.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==2
data$WWEAM5.2 <- data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==2
data$WWEAM5.3 <- (data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==2)
data$RREAM5.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==3 & data$WE==2
data$RREAM5.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==3 & data$WE==2
data$RREAM5.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==3 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==2)
data$BBEAM5.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==2
data$BBEAM5.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==3 & data$WE==2
data$BBEAM5.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==2)

data$WWEAM6.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==3
data$WWEAM6.2 <- data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==3
data$WWEAM6.3 <- (data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==3)
data$RREAM6.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==3 & data$WE==3
data$RREAM6.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==3 & data$WE==3
data$RREAM6.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==3 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==3)
data$BBEAM6.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==3
data$BBEAM6.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==3 & data$WE==3
data$BBEAM6.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==3)

data$WWEAM7.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==4
data$WWEAM7.2 <- data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==4
data$WWEAM7.3 <- (data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==4)
data$RREAM7.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==3 & data$WE==4
data$RREAM7.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==3 & data$WE==4
data$RREAM7.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==3 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==4)
data$BBEAM7.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==4
data$BBEAM7.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==3 & data$WE==4
data$BBEAM7.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==4)

data$WWEAM8.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==5
data$WWEAM8.2 <- data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==5
data$WWEAM8.3 <- (data$HR=="W" & data$WR=="W" & data$HE==3 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==3 & data$WE==5)
data$RREAM8.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==3 & data$WE==5
data$RREAM8.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==3 & data$WE==5
data$RREAM8.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==3 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==5)
data$BBEAM8.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==5
data$BBEAM8.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==3 & data$WE==5
data$BBEAM8.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==3 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==3 & data$WE==5)

data$WWEAM9.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==2
data$WWEAM9.2 <- data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==2
data$WWEAM9.3 <- (data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==2)
data$RREAM9.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==4 & data$WE==2
data$RREAM9.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==4 & data$WE==2
data$RREAM9.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==4 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==2)
data$BBEAM9.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==2
data$BBEAM9.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==4 & data$WE==2
data$BBEAM9.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==2)

data$WWEAM10.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==3
data$WWEAM10.2 <- data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==3
data$WWEAM10.3 <- (data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==3)
data$RREAM10.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==4 & data$WE==3
data$RREAM10.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==4 & data$WE==3
data$RREAM10.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==4 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==3)
data$BBEAM10.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==3
data$BBEAM10.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==4 & data$WE==3
data$BBEAM10.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==3)

data$WWEAM11.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==4
data$WWEAM11.2 <- data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==4
data$WWEAM11.3 <- (data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==4)
data$RREAM11.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==4 & data$WE==4
data$RREAM11.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==4 & data$WE==4
data$RREAM11.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==4 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==4)
data$BBEAM11.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==4
data$BBEAM11.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==4 & data$WE==4
data$BBEAM11.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==4)

data$WWEAM12.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==5
data$WWEAM12.2 <- data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==5
data$WWEAM12.3 <- (data$HR=="W" & data$WR=="W" & data$HE==4 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==4 & data$WE==5)
data$RREAM12.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==4 & data$WE==5
data$RREAM12.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==4 & data$WE==5
data$RREAM12.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==4 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==5)
data$BBEAM12.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==5
data$BBEAM12.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==4 & data$WE==5
data$BBEAM12.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==4 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==4 & data$WE==5)

data$WWEAM13.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==2
data$WWEAM13.2 <- data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==2
data$WWEAM13.3 <- (data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==2)
data$RREAM13.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==5 & data$WE==2
data$RREAM13.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==5 & data$WE==2
data$RREAM13.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==5 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==2)
data$BBEAM13.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==2
data$BBEAM13.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==5 & data$WE==2
data$BBEAM13.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==2) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==2)

data$WWEAM14.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==3
data$WWEAM14.2 <- data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==3
data$WWEAM14.3 <- (data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==3)
data$RREAM14.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==5 & data$WE==3
data$RREAM14.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==5 & data$WE==3
data$RREAM14.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==5 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==3)
data$BBEAM14.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==3
data$BBEAM14.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==5 & data$WE==3
data$BBEAM14.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==3) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==3)

data$WWEAM15.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==4
data$WWEAM15.2 <- data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==4
data$WWEAM15.3 <- (data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==4)
data$RREAM15.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==5 & data$WE==4
data$RREAM15.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==5 & data$WE==4
data$RREAM15.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==5 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==4)
data$BBEAM15.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==4
data$BBEAM15.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==5 & data$WE==4
data$BBEAM15.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==4) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==4)

data$WWEAM16.1 <- ((data$HR=="W" & data$WR=="W") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==5
data$WWEAM16.2 <- data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==5
data$WWEAM16.3 <- (data$HR=="W" & data$WR=="W" & data$HE==5 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W")  | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W"))  & data$HE==5 & data$WE==5)
data$RREAM16.1 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br")) & data$HE==5 & data$WE==5
data$RREAM16.2 <- ((data$HR=="Br" & data$WR=="Br") | (data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W")) & data$HE==5 & data$WE==5
data$RREAM16.3 <- (data$HR=="Br" & data$WR=="Br" & data$HE==5 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Br") | (data$HR=="Br" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==5)
data$BBEAM16.1 <- data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==5
data$BBEAM16.2 <- ((data$HR=="Bl" & data$WR=="Bl") | (data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Br" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="Br"))  & data$HE==5 & data$WE==5
data$BBEAM16.3 <- (data$HR=="Bl" & data$WR=="Bl" & data$HE==5 & data$WE==5) + 0.5 * (((data$HR=="W" & data$WR=="Bl") | (data$HR=="Bl" & data$WR=="W") | (data$HR=="Bl" & data$WR=="Br") | (data$HR=="Br" & data$WR=="Bl"))  & data$HE==5 & data$WE==5)


data2 <- data

data2$HE <- as.factor(data2$HE)
data2$WE <- as.factor(data2$WE)

basemodel.full <- glm(I(round(Freq,0))~HR*WR+HE*HR+WE*WR, family=poisson, data=data2)
basemodel1.full <- update(basemodel.full,.~.+WWEAM1.1+WWEAM2.1+WWEAM3.1+WWEAM4.1+WWEAM5.1+WWEAM6.1+WWEAM7.1+WWEAM8.1+
									+WWEAM9.1+WWEAM10.1+WWEAM11.1+WWEAM12.1+WWEAM13.1+WWEAM14.1+WWEAM15.1+WWEAM16.1
									+RREAM1.1+RREAM2.1+RREAM3.1+RREAM4.1+RREAM5.1+RREAM6.1+RREAM7.1+RREAM8.1+
									+RREAM9.1+RREAM10.1+RREAM11.1+RREAM12.1+RREAM13.1+RREAM14.1+RREAM15.1+RREAM16.1
									+BBEAM1.1+BBEAM2.1+BBEAM3.1+BBEAM4.1+BBEAM5.1+BBEAM6.1+BBEAM7.1+BBEAM8.1+
									+BBEAM9.1+BBEAM10.1+BBEAM11.1+BBEAM12.1+BBEAM13.1+BBEAM14.1+BBEAM15.1+BBEAM16.1)
semodel1.full <- update(basemodel1.full, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown+SERMWFUp+SERMWFDown+SEWMRFUp+SEWMRFDown+SEBMRFUp+SEBMRFDown+SERMBFUp+SERMBFDown)
ebmodel1.full <- update(semodel1.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)
ebomodel1.full <- update(basemodel1.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)

basemodel2.full <- update(basemodel.full,.~.+WWEAM1.2+WWEAM2.2+WWEAM3.2+WWEAM4.2+WWEAM5.2+WWEAM6.2+WWEAM7.2+WWEAM8.2+
									+WWEAM9.2+WWEAM10.2+WWEAM11.2+WWEAM12.2+WWEAM13.2+WWEAM14.2+WWEAM15.2+WWEAM16.2
									+RREAM1.2+RREAM2.2+RREAM3.2+RREAM4.2+RREAM5.2+RREAM6.2+RREAM7.2+RREAM8.2+
									+RREAM9.2+RREAM10.2+RREAM11.2+RREAM12.2+RREAM13.2+RREAM14.2+RREAM15.2+RREAM16.2
									+BBEAM1.2+BBEAM2.2+BBEAM3.2+BBEAM4.2+BBEAM5.2+BBEAM6.2+BBEAM7.2+BBEAM8.2+
									+BBEAM9.2+BBEAM10.2+BBEAM11.2+BBEAM12.2+BBEAM13.2+BBEAM14.2+BBEAM15.2+BBEAM16.2)
semodel2.full <- update(basemodel2.full, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown+SERMWFUp+SERMWFDown+SEWMRFUp+SEWMRFDown+SEBMRFUp+SEBMRFDown+SERMBFUp+SERMBFDown)
ebmodel2.full <- update(semodel2.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)
ebomodel2.full <- update(basemodel2.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)
basemodel3.full <- update(basemodel.full,.~.+WWEAM1.3+WWEAM2.3+WWEAM3.3+WWEAM4.3+WWEAM5.3+WWEAM6.3+WWEAM7.3+WWEAM8.3+
									+WWEAM9.3+WWEAM10.3+WWEAM11.3+WWEAM12.3+WWEAM13.3+WWEAM14.3+WWEAM15.3+WWEAM16.3
									+RREAM1.3+RREAM2.3+RREAM3.3+RREAM4.3+RREAM5.3+RREAM6.3+RREAM7.3+RREAM8.3+
									+RREAM9.3+RREAM10.3+RREAM11.3+RREAM12.3+RREAM13.3+RREAM14.3+RREAM15.3+RREAM16.3
									+BBEAM1.3+BBEAM2.3+BBEAM3.3+BBEAM4.3+BBEAM5.3+BBEAM6.3+BBEAM7.3+BBEAM8.3+
									+BBEAM9.3+BBEAM10.3+BBEAM11.3+BBEAM12.3+BBEAM13.3+BBEAM14.3+BBEAM15.3+BBEAM16.3)
semodel3.full <- update(basemodel3.full, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown+SERMWFUp+SERMWFDown+SEWMRFUp+SEWMRFDown+SEBMRFUp+SEBMRFDown+SERMBFUp+SERMBFDown)
ebmodel3.full <- update(semodel3.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)
ebomodel3.full <- update(basemodel3.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)
basemodel4.full <- update(basemodel.full,.~.+HE*WE)
semodel4.full <- update(basemodel4.full, .~.+SEBMWFUp+SEBMWFDown+SEWMBFUp+SEWMBFDown+SERMWFUp+SERMWFDown+SEWMRFUp+SEWMRFDown+SEBMRFUp+SEBMRFDown+SERMBFUp+SERMBFDown)
ebmodel4.full <- update(semodel4.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)
ebomodel4.full <- update(basemodel4.full, .~.+EBBMW1+EBBMW2+EBBMW3+EBBMW4+EBWFB1+EBWFB2+EBWFB3+EBWFB4
									+EBWMB1+EBWMB2+EBWMB3+EBWMB4+EBBFW1+EBBFW2+EBBFW3+EBBFW4
									+EBBMR1+EBBMR2+EBBMR3+EBBMR4+EBRFB1+EBRFB2+EBRFB3+EBRFB4
									+EBBFR1+EBBFR2+EBBFR3+EBBFR4+EBRMB1+EBRMB2+EBRMB3+EBRMB4
									+EBRMW1+EBRMW2+EBRMW3+EBRMW4+EBWFR1+EBWFR2+EBWFR3+EBWFR4
									+EBRFW1+EBRFW2+EBRFW3+EBRFW4+EBWMR1+EBWMR2+EBWMR3+EBWMR4)

semodel4_scale.full <- update(basemodel4.full,.~.+SEDMUpScale+SEDMDownScale+SEDFUpScale+SEDFDownScale)
semodel4_scale2.full <- update(basemodel4.full,.~.+SEDMUpScale2+SEDMDownScale2+SEDFUpScale2+SEDFDownScale2)

ebmodel4_scale.full <- update(semodel4.full, .~.+EBDM1+EBDM2+EBDM3+EBDM4+EBDF1+EBDF2+EBDF3+EBDF4+EBLM1+EBLM2+EBLM3+EBLM4+EBLF1+EBLF2+EBLF3+EBLF4)
ebmodel4_scale2.full <- update(semodel4.full, .~.+EBDM1.2+EBDM2.2+EBDM3.2+EBDM4.2+EBDF1.2+EBDF2.2+EBDF3.2+EBDF4.2
											+EBLM1.2+EBLM2.2+EBLM3.2+EBLM4.2+EBLF1.2+EBLF2.2+EBLF3.2+EBLF4.2)
ebmodel4_scale3.full <- update(semodel4.full, .~.+EBDM1.2+EBDM2.2+EBDM3.2+EBDM4.2+EBDF1.2+EBDF2.2+EBDF3.2+EBDF4.2
											+EBLM1+EBLM2+EBLM3+EBLM4+EBLF1+EBLF2+EBLF3+EBLF4)
											