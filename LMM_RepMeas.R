## Run LMM Models on Time-Series Data ##
## Janssen Data using Multiple-Measures ##
## July 29, 2014 ##
## Kristopher Standish ##

library(nlme)

##############################################################
## LOAD DATA #################################################
##############################################################

## Set Date
DATE <- "20150226"

# ## TSCC Paths
# PathToData <- "/projects/janssen/clinical/"
# PathToSave <- "/projects/janssen/clinical/Plots/20140522/"

## Mac Paths
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,sep="" )

## Previously Compiled Data
TAB <- read.table( PathToData, sep="\t",header=T )

##############################################################
## LINEAR MODELS #############################################
##############################################################

## Basic Linear Model w/ only DRUG effect
LM.1 <- lm( DAS ~ DRUG, data=TAB )
summary(LM.1)

## Add DRUG*SEX Interaction
LM.2 <- lm( DAS ~ DRUG*SEX, data=TAB )
summary(LM.2)
## Change to DRUG+SEX
LM.3 <- lm( DAS ~ DRUG+SEX, data=TAB )
summary(LM.3)

## Remove SEX, add DRUG*AGE
LM.4 <- lm( DAS ~ DRUG*AGE, data=TAB )
summary(LM.4)
## Change to DRUG+AGE
LM.5 <- lm( DAS ~ DRUG+AGE, data=TAB )
summary(LM.5)

## Remove SEX, add DRUG*DIS_DUR
LM.6 <- lm( DAS ~ DRUG*DIS_DUR, data=TAB )
summary(LM.6)
## Change to DRUG+DIS_DUR
LM.7 <- lm( DAS ~ DRUG+DIS_DUR, data=TAB )
summary(LM.7)

##############################################################
## LMM MODELS ################################################
##############################################################

## When using ANOVA to compare models:
 # Greater (less negative) logLik is better model
 # Smaller (less positive) AIC/BIC is better model

## Linear Mixed Model w/ only DRUG effect
LME.1 <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB )
summary(LME.1)

## Add Drug as Random Effect
LME.1b <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB )
summary(LME.1b)
anova(LME.1,LME.1b)
 # Choose LME.2

## Add DRUG*SEX Interaction
LME.2 <- lme( fixed = DAS ~ DRUG*SEX, random = ~ DRUG | IID, data=TAB )
summary(LME.2)
## Change to DRUG+SEX
LME.3 <- lme( fixed = DAS ~ DRUG+SEX, random = ~ DRUG | IID, data=TAB )
summary(LME.3)

## Remove SEX, add DRUG*AGE
LME.4 <- lme( fixed = DAS ~ DRUG*AGE, random = ~ DRUG | IID, data=TAB )
summary(LME.4)
## Change to DRUG+AGE
LME.5 <- lme( fixed = DAS ~ DRUG+AGE, random = ~ DRUG | IID, data=TAB )
summary(LME.5)

## Remove AGE, add DRUG*DIS_DUR
LME.6 <- lme( fixed = DAS ~ DRUG*DIS_DUR, random = ~ DRUG | IID, data=TAB )
summary(LME.6)
## Remove AGE, add DRUG+DIS_DUR
LME.7 <- lme( fixed = DAS ~ DRUG+DIS_DUR, random = ~ DRUG | IID, data=TAB )
summary(LME.7)

##############################################################
## LIST OF LINEAR MODELS #####################################
##############################################################

## Linear Mixed Model w/ only DRUG effect
TAB.G1 <- groupedData( DAS ~ DRUG | IID, data=TAB) # , outer=~SEX )
TAB.G1 <- TAB.G1[,-c("CRP","SJC","TJC","SJC28","TJC28")]
LS.1 <- lmList( DAS ~ DRUG, data=TAB.G1 )
summary(LS.1)

## Add
LME.2 <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB )
summary(LME.2)















##############################################################
## END OF DOC ################################################
##############################################################
