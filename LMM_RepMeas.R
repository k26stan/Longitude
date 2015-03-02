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
# Get unique Weeks for plotting purposes
WKS <- unique( TAB$WK )

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
LM.2a <- lm( DAS ~ DRUG+SEX, data=TAB )
summary(LM.2a)

## Remove SEX, add DRUG*AGE
LM.3 <- lm( DAS ~ DRUG*AGE, data=TAB )
summary(LM.3)
## Change to DRUG+AGE
LM.3a <- lm( DAS ~ DRUG+AGE, data=TAB )
summary(LM.3a)

## Remove SEX, add DRUG*DIS_DUR
LM.4 <- lm( DAS ~ DRUG*DIS_DUR, data=TAB )
summary(LM.4)
## Change to DRUG+DIS_DUR
LM.4a <- lm( DAS ~ DRUG+DIS_DUR, data=TAB )
summary(LM.4a)

## Remove DIS_DUR, add DRUG*HT
LM.5 <- lm( DAS ~ DRUG*HT, data=TAB )
summary(LM.5)
## Change to DRUG+HT
LM.5a <- lm( DAS ~ DRUG+HT, data=TAB )
summary(LM.5a)

## Change to DRUG*WT
LM.6 <- lm( DAS ~ DRUG*WT, data=TAB )
summary(LM.6)
## Change to DRUG+WT
LM.6a <- lm( DAS ~ DRUG+WT, data=TAB )
summary(LM.6a)

## Change to DRUG*BMI
LM.7 <- lm( DAS ~ DRUG*BMI, data=TAB )
summary(LM.7)
## Change to DRUG+BMI
LM.7a <- lm( DAS ~ DRUG+BMI, data=TAB )
summary(LM.7a)

## Change to DRUG*ACPA
LM.8 <- lm( DAS ~ DRUG*ACPA, data=TAB, subset=ACPA!="" )
summary(LM.8)
## Change to DRUG+ACPA
LM.8a <- lm( DAS ~ DRUG+ACPA, data=TAB, subset=ACPA!="" )
summary(LM.8a)

## Change to DRUG*RF
LM.9 <- lm( DAS ~ DRUG*RF, data=TAB )
summary(LM.9)
## Change to DRUG+RF
LM.9a <- lm( DAS ~ DRUG+RF, data=TAB )
summary(LM.9a)

## Change to DRUG*RF_ACPA
LM.10 <- lm( DAS ~ DRUG*RF_ACPA, data=TAB )
summary(LM.10)
## Change to DRUG+RF_ACPA
LM.10a <- lm( DAS ~ DRUG+RF_ACPA, data=TAB )
summary(LM.10a)

## Change to DRUG*WK
LM.11 <- lm( DAS ~ DRUG*WK, data=TAB )
summary(LM.11)
## Change to DRUG+WK
LM.11a <- lm( DAS ~ DRUG+WK, data=TAB )
summary(LM.11a)

## Include PLAC (all 2-way interactions)
LM.12 <- lm( DAS ~ (DRUG+WK+PLAC)^2, data=TAB )
summary(LM.12)
## Include PLAC (2-way interactions w/o PLAC:DRUG)
LM.12a <- lm( DAS ~ DRUG*WK+PLAC+PLAC:WK, data=TAB )
summary(LM.12a)
## Include PLAC (w/o PLAC interactions)
LM.12b <- lm( DAS ~ DRUG*WK+PLAC, data=TAB )
summary(LM.12b)

## Use 12b, add SEX interactions
LM.13 <- lm( DAS ~ SEX*(DRUG*WK+PLAC), data=TAB )
summary(LM.13)
## Remove SEX:DRUG:WK
LM.13a <- lm( DAS ~ SEX*(DRUG*WK+PLAC)-SEX:DRUG:WK, data=TAB )
summary(LM.13a)
## Remove SEX:WK
LM.13b <- lm( DAS ~ SEX*(DRUG*WK+PLAC)-SEX:DRUG:WK-SEX:WK, data=TAB )
summary(LM.13b)
## Remove SEX:PLAC
LM.13c <- lm( DAS ~ (DRUG*WK+PLAC)+SEX*DRUG, data=TAB )
summary(LM.13c)
## Remove SEX:DRUG
LM.13d <- lm( DAS ~ (DRUG*WK+PLAC)+SEX, data=TAB )
summary(LM.13d)

## Use 12b, add RF_ACPA interactions
LM.14 <- lm( DAS ~ RF_ACPA*(DRUG*WK+PLAC), data=TAB )
summary(LM.14)
## Remove RF_ACPA:DRUG:WK
LM.14a <- lm( DAS ~ RF_ACPA*(DRUG*WK+PLAC)-RF_ACPA:DRUG:WK, data=TAB )
summary(LM.14a)
## Remove RF_ACPA:WK
LM.14b <- lm( DAS ~ RF_ACPA*(DRUG*WK+PLAC)-RF_ACPA:DRUG:WK-RF_ACPA:WK, data=TAB )
summary(LM.14b)
## Remove RF_ACPA:PLAC
LM.14c <- lm( DAS ~ (DRUG*WK+PLAC)+RF_ACPA*DRUG, data=TAB )
summary(LM.14c)
## Remove RF_ACPA:DRUG
LM.14d <- lm( DAS ~ (DRUG*WK+PLAC)+RF_ACPA, data=TAB )
summary(LM.14d)

## Use 12b, add DIS_DUR interactions
LM.15 <- lm( DAS ~ DIS_DUR*(DRUG*WK+PLAC), data=TAB )
summary(LM.15)
## Remove DIS_DUR:DRUG:WK
LM.15a <- lm( DAS ~ DIS_DUR*(DRUG*WK+PLAC)-DIS_DUR:DRUG:WK, data=TAB )
summary(LM.15a)
## Remove DIS_DUR:WK
LM.15b <- lm( DAS ~ DIS_DUR*(DRUG*WK+PLAC)-DIS_DUR:DRUG:WK-DIS_DUR:WK, data=TAB )
summary(LM.15b)
## Remove DIS_DUR:PLAC
LM.15c <- lm( DAS ~ (DRUG*WK+PLAC)+DIS_DUR*DRUG, data=TAB )
summary(LM.15c)
## Remove DIS_DUR:DRUG
LM.15d <- lm( DAS ~ (DRUG*WK+PLAC)+DIS_DUR, data=TAB )
summary(LM.15d)

## Use 12b, add BMI interactions
LM.16 <- lm( DAS ~ BMI*(DRUG*WK+PLAC), data=TAB )
summary(LM.16)
## Remove BMI:DRUG:WK
LM.16a <- lm( DAS ~ BMI*(DRUG*WK+PLAC)-BMI:DRUG:WK, data=TAB )
summary(LM.16a)
## Remove BMI:WK
LM.16b <- lm( DAS ~ BMI*(DRUG*WK+PLAC)-BMI:DRUG:WK-BMI:WK, data=TAB )
summary(LM.16b)
## Remove BMI:PLAC
LM.16c <- lm( DAS ~ (DRUG*WK+PLAC)+BMI*DRUG, data=TAB )
summary(LM.16c)
## Remove BMI:DRUG
LM.16d <- lm( DAS ~ (DRUG*WK+PLAC)+BMI, data=TAB )
summary(LM.16d)


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
LME.2a <- lme( fixed = DAS ~ DRUG+SEX, random = ~ DRUG | IID, data=TAB )
summary(LME.2a)

## Change to DRUG*AGE
LME.3 <- lme( fixed = DAS ~ DRUG*AGE, random = ~ DRUG | IID, data=TAB )
summary(LME.3)
## Change to DRUG+AGE
LME.3a <- lme( fixed = DAS ~ DRUG+AGE, random = ~ DRUG | IID, data=TAB )
summary(LME.3a)

## Change to DRUG*DIS_DUR
LME.4 <- lme( fixed = DAS ~ DRUG*DIS_DUR, random = ~ DRUG | IID, data=TAB )
summary(LME.4)
## Change to DRUG+DIS_DUR
LME.4a <- lme( fixed = DAS ~ DRUG+DIS_DUR, random = ~ DRUG | IID, data=TAB )
summary(LME.4a)

## Remove DIS_DUR, add DRUG*HT
LME.5 <- lme( fixed = DAS ~ DRUG*HT, random = ~ DRUG | IID, data=TAB )
summary(LME.5)
## Change to DRUG+HT
LME.5a <- lme( fixed = DAS ~ DRUG+HT, random = ~ DRUG | IID, data=TAB )
summary(LME.5a)

## Change to DRUG*WT
LME.6 <- lme( fixed = DAS ~ DRUG*WT, random = ~ DRUG | IID, data=TAB )
summary(LME.6)
## Change to DRUG+WT
LME.6a <- lme( fixed = DAS ~ DRUG+WT, random = ~ DRUG | IID, data=TAB )
summary(LME.6a)

## Change to DRUG*BMI
LME.7 <- lme( fixed = DAS ~ DRUG*BMI, random = ~ DRUG | IID, data=TAB )
summary(LME.7)
## Change to DRUG+BMI
LME.7a <- lme( fixed = DAS ~ DRUG+BMI, random = ~ DRUG | IID, data=TAB )
summary(LME.7a)

## Change to DRUG*ACPA
LME.8 <- lme( fixed = DAS ~ DRUG*ACPA, random = ~ DRUG | IID, data=TAB, subset=ACPA!="" )
summary(LME.8)
## Change to DRUG+ACPA
LME.8a <- lme( fixed = DAS ~ DRUG+ACPA, random = ~ DRUG | IID, data=TAB, subset=ACPA!="" )
summary(LME.8a)

## Change to DRUG*RF
LME.9 <- lme( fixed = DAS ~ DRUG*RF, random = ~ DRUG | IID, data=TAB )
summary(LME.9)
## Change to DRUG+RF
LME.9a <- lme( fixed = DAS ~ DRUG+RF, random = ~ DRUG | IID, data=TAB )
summary(LME.9a)

## Change to DRUG*RF_ACPA
LME.10 <- lme( fixed = DAS ~ DRUG*RF_ACPA, random = ~ DRUG | IID, data=TAB )
summary(LME.10)
## Change to DRUG+RF_ACPA
LME.10a <- lme( fixed = DAS ~ DRUG+RF_ACPA, random = ~ DRUG | IID, data=TAB )
summary(LME.10a)

## Change to DRUG*WK
LME.11 <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB )
summary(LME.11)
## Change to DRUG+WK
LME.11a <- lme( fixed = DAS ~ DRUG+WK, random = ~ DRUG | IID, data=TAB )
summary(LME.11a)
## Back to DRUG*WK, add WK & DRUG*WK as Random Effect
LME.11b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.11b)
anova(LME.11,LME.11b)

## Include PLAC (all 2-way interactions)
LME.12 <- lme( fixed = DAS ~ (DRUG+WK+PLAC)^2, random = ~ DRUG | IID, data=TAB )
summary(LME.12)
## Include PLAC (2-way interactions w/o PLAC:DRUG)
LME.12a <- lme( fixed = DAS ~ DRUG*WK+PLAC+PLAC:WK, random = ~ DRUG | IID, data=TAB )
summary(LME.12a)
## Include PLAC (w/o PLAC interactions)
LME.12b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.12b)
## Add PLAC as Random Effect
LME.12c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK+PLAC | IID, data=TAB )
summary(LME.12c) # !!! DID NOT CONVERGE !!!
## Try w/ only DRUG as Random Effect
LME.12d <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | IID, data=TAB )
summary(LME.12d)
anova(LME.12b,LME.12d) # Stick w/ WK in Random Effects (b)

## Use 12b, add SEX interactions
LME.13 <- lme( fixed = DAS ~ SEX*(DRUG*WK+PLAC), random = ~ DRUG+WK | IID, data=TAB )
summary(LME.13)
## Remove SEX:DRUG:WK
LME.13a <- lme( fixed = DAS ~ SEX*(DRUG*WK+PLAC)-SEX:DRUG:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.13a)
## Remove SEX:WK
LME.13b <- lme( fixed = DAS ~ SEX*(DRUG*WK+PLAC)-SEX:DRUG:WK-SEX:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.13b)
## Remove SEX:PLAC
LME.13c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+SEX*DRUG, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.13c)
## Remove SEX:DRUG
LME.13d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+SEX, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.13d)
 # !!! No effect of SEX at all !!!

## Use 12b, add RF_ACPA interactions
LME.14 <- lme( fixed = DAS ~ RF_ACPA*(DRUG*WK+PLAC), random = ~ DRUG+WK | IID, data=TAB )
summary(LME.14)
## Remove RF_ACPA:DRUG:WK
LME.14a <- lme( fixed = DAS ~ RF_ACPA*(DRUG*WK+PLAC)-RF_ACPA:DRUG:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.14a)
## Remove RF_ACPA:WK
LME.14b <- lme( fixed = DAS ~ RF_ACPA*(DRUG*WK+PLAC)-RF_ACPA:DRUG:WK-RF_ACPA:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.14b)
## Remove RF_ACPA:PLAC
LME.14c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+RF_ACPA*DRUG, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.14c)
## Remove RF_ACPA:DRUG
LME.14d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+RF_ACPA, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.14d)
## Include RF_ACPA as Grouping
LME.14c1 <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | RF_ACPA/IID, data=TAB )
summary(LME.14c1) ## Didn't Converge...first try...try again...
## Include RF_ACPA as Random Effect
LME.14da <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+RF_ACPA, random = ~ DRUG+WK | RF_ACPA/IID, data=TAB )
summary(LME.14da)
## Decide: w/ or w/o RF_ACPA:DRUG interaction??
cor(data.frame( TAB$DAS, predict(LME.14c), predict(LME.14d) ))
sd(resid(LME.14c)) ; sd(resid(LME.14d)) ; sd(resid(LME.12b))
 # !!! Include RF_ACPA, but not RF_ACPA:DRUG !!!

## Use 12b, add DIS_DUR interactions
LME.15 <- lme( fixed = DAS ~ DIS_DUR*(DRUG*WK+PLAC), random = ~ DRUG+WK | IID, data=TAB )
summary(LME.15)
## Remove DIS_DUR:DRUG:WK
LME.15a <- lme( fixed = DAS ~ DIS_DUR*(DRUG*WK+PLAC)-DIS_DUR:DRUG:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.15a)
## Remove DIS_DUR:WK
LME.15b <- lme( fixed = DAS ~ DIS_DUR*(DRUG*WK+PLAC)-DIS_DUR:DRUG:WK-DIS_DUR:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.15b)
## Remove DIS_DUR:PLAC
LME.15c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DIS_DUR*DRUG, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.15c)
## Remove DIS_DUR:DRUG
LME.15d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DIS_DUR, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.15d)
## Decide: w/ or w/o RF_ACPA:DRUG interaction??
cor(data.frame( TAB$DAS, predict(LME.15b), predict(LME.15c), predict(LME.15d) ))
sd(resid(LME.15b)) ; sd(resid(LME.15c)) ; sd(resid(LME.15d)) ; sd(resid(LME.12b))
 # !!! Include ... !!!

## Use 12b, add BMI interactions
LME.16 <- lme( fixed = DAS ~ BMI*(DRUG*WK+PLAC), random = ~ DRUG+WK | IID, data=TAB )
summary(LME.16)
## Remove BMI:DRUG:WK
LME.16a <- lme( fixed = DAS ~ BMI*(DRUG*WK+PLAC)-BMI:DRUG:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.16a)
## Remove BMI:WK
LME.16b <- lme( fixed = DAS ~ BMI*(DRUG*WK+PLAC)-BMI:DRUG:WK-BMI:WK, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.16b)
## Remove BMI:PLAC
LME.16c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+BMI*DRUG, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.16c)
## Remove BMI:DRUG
LME.16d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+BMI, random = ~ DRUG+WK | IID, data=TAB )
summary(LME.16d)
## Decide: w/ or w/o RF_ACPA:DRUG interaction??
cor(data.frame( TAB$DAS, predict(LME.16b), predict(LME.16c), predict(LME.16d) ))
sd(resid(LME.16b)) ; sd(resid(LME.16c)) ; sd(resid(LME.16d)) ; sd(resid(LME.12b))
 # !!! Include ... !!!

head(coef(LME.12b))
head(coef(LME.14b))
##############################################################
## PLOT SOME PREDICTED VALUES ################################
##############################################################

## Histogram of Coefficients (random Effects)
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
COLS <- sample(COLS.list)
par(mfrow=c(1,4))
hist( coef(LME.12b)[,"(Intercept)"], col=COLS[1] )
hist( coef(LME.12b)[,"DRUG"], col=COLS[2] )
hist( coef(LME.12b)[,"WK"], col=COLS[3] )
hist( coef(LME.12b)[,"WK"]+coef(LME.12b)[,"DRUG:WK"], col=COLS[4] )

## All Data Points
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
COLS <- COLS.list
LIM <- c(0,10)
png( paste(PathToPlot,"/Pl_Scatter_ObsExp.png",sep=""), height=2400,width=1800, pointsize=26 )
par(mfrow=c(4,3))
plot( data.frame( TAB$DAS, predict(LME.1) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LME.1b) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LM.1) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LME.11) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LME.11b) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LM.11) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( 0,0,type="n", col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LME.12b) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LM.12) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( 0,0,type="n", col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LME.14c) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
plot( data.frame( TAB$DAS, predict(LM.14) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
abline( 0,1, lty=2 )
dev.off()

## All Data Points
pairs( data.frame( TAB$DAS, predict(LME.1b),predict(LME.11b),predict(LME.12b) ), col=COLS[factor(TAB$PLAC)], pch="." )

## A few specific Patients
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
N.patients <- 4
First.patient <- sample(1:436,1) # 3
COLS <- colorRampPalette(sample(COLS.list))(N.patients)
plot( 0,0,type="n", xlim=c(0,100),ylim=c(0,10),xlab="Weeks",ylab="DAS", main="Predicted vs Actual DAS Over Time" )
abline( h=seq(0,10,1), lty=2,col="grey50",lwd=1 )
for ( n in 1:N.patients ) {
	samp <- unique(TAB$IID)[First.patient+n-1]
	ind <- which(TAB$IID==samp)
	PCHS <- c(1,16)[factor(TAB$DRUG[ind])]
	OBS <- TAB[ ind, "DAS" ]
	EXP.LME.1 <- predict(LME.1b)[ind]
	EXP.LME.12 <- predict(LME.12b)[ind]
	points( TAB[ind,"WK"], OBS, col=COLS[n], type="o",pch=PCHS,lty=1,lwd=2 )
	points( TAB[ind,"WK"], EXP.LME.1, col=COLS[n], type="o",pch=PCHS,lty=2,lwd=2 )
	points( TAB[ind,"WK"], EXP.LME.12, col=COLS[n], type="o",pch=PCHS,lty=3,lwd=2 )
}
coef(LME.12b)[ First.patient:(First.patient+N.patients-1) , ]

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
