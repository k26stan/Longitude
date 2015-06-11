## Do Bayesian Analysis w/ Mixed Model & Longitudinal Data ##
## Janssen Data using Multiple-Measures ##
## May 22, 2015 ##
## Kristopher Standish ##

## Game Plan: Follow N-of-1 Template ##
 # Use only Patients who were on Placebo also
   # And Remove Baseline Measure??
 # Fill out Table 1A
   # Calculate Mean Difference b/n Drug & Placebo (per person)
   # Calculate Significance (based on 1-sided unpaired t-test)
 # Calculate Basic Mixed Effects Model
   # w/ Drug as Random Effect
   # w/ IID as Grouping Factor
 # Calculate Prior Distributions for:
   # Population Mean effect size (from Mixed Model)
   # b/n Patient Variance (from Mixed Model)
   # w/i Patient Variance (from Mixed Model)


##############################################################
## LOAD DATA #################################################
##############################################################
library(nlme)
library(gplots)

## Set Date
DATE <- "20150610"

## Mac Paths
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150506_Resp_v_Time.txt"
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150512_Resp_v_Time.txt"
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,sep="" )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( paste(PathToFT,"20141229_Full_Table.txt",sep=""),sep="\t",header=T)
TAB.2 <- merge( TAB.l, FUL[,c("ID_2",paste("PC",1:3,sep=""))], by.x="IID",by.y="ID_2")
# Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Load Candidate Genotype Files
COMP.l <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", sep="\t",header=T)
COMP <- COMP.l[ which(!duplicated(COMP.l$SNP)), ]
RAW <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/TEST_CND.raw", sep="",header=T)

##############################################################
## FILTER DATA ###############################################
##############################################################

TAB <- TAB.l

## Take out Patients who left before getting DRUG
RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
RM.exit <- which( TAB$IID %in% RM.exit.id )
TAB <- TAB[-RM.exit,c(1:15,17)]

## Take out Patients from Golimumab Arm
RM.gol.id <- as.character( FUL$ID_2[which(FUL$GRP=="G")] )
RM.gol <- which( TAB$IID %in% RM.gol.id )
TAB <- TAB[ -RM.gol, ]

## Take out WK==0
TAB <- TAB[ which(TAB$WK!=0), ]

## Remove NA Values for DAS
TAB <- TAB[ which(!is.na(TAB$DAS)), ]
dim(TAB)

##############################################################
## LMM MODELS ################################################
##############################################################
	## When using ANOVA to compare models:
	 # Greater (less negative) logLik is better model
	 # Smaller (less positive) AIC/BIC is better model

## Mixed Effects Model
LME.1 <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB )
LME.2 <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB )
LME.3 <- lme( fixed = DAS ~ DRUG, random = ~ DRUG - 1 | IID, data=TAB )
LME.4 <- lme( fixed = DAS ~ DRUG - 1, random = ~ DRUG - 1 | IID, data=TAB )
 # Got with LME.2
LME <- LME.2

## Compile Stats for Bayesian Inference
 # mu(i) estimates
 # mu(0) estimate
 # tau estimate
 # sigma(i) estimates
MU.0 <- fixef(LME)
MU.i <- coef(LME)
TAU <- as.numeric( VarCorr(LME)[1:2,"StdDev"] ) ; names(TAU) <- names(MU.0)
# TAU <- intervals( LME, which="var-cov" )
SIG.i <- aggregate( x=data.frame(resid(LME)), by=list(TAB$IID), sd )
# SIG.i <- LME$sigma

POST <- merge( MU.i, SIG.i, by.x="row.names",by.y="Group.1", all=T )
colnames(POST) <- c("ID","INT","PMD","PSD")

 # Plot Posterior Distributions for Individuals
XX <- seq( -10,10,.05 )
COLS.list <- c("black","slateblue2","steelblue2","springgreen2","gold1","chocolate2","firebrick2")
COLS <- colorRampPalette(COLS.list)(nrow(POST))
PROBS <- numeric( nrow(POST) )
plot( 0,0,type="n",xlim=c(-5,5),ylim=c(0,2),xlab="Posterior Mean Difference",ylab="" )
abline( v=0,lty=2,lwd=2 )
N <- 20 ; WHICH <- sample(1:nrow(POST),N)
for ( i in 1:nrow(POST) ) {
# for ( w in 1:N ) {
   # i <- WHICH[w]
   points( XX, dnorm( XX,POST[i,"PMD"],POST[i,"PSD"]), type="l",col=COLS[i] )
   PROBS[i] <- pnorm( 0, POST[i,"PMD"],POST[i,"PSD"], lower.tail=T )
   text( 3,2-quantile(c(0,2),(w-1)/N), label=round(PROBS[i],2),col=COLS[i] )
}

##############################################################
## LIS MODELS ################################################
##############################################################

## Independent Linear Regression Models
LIS <- lmList( DAS ~ DRUG | IID, data=TAB )
LIS.I <- summary(LIS)$coefficients[,,"(Intercept)"]
LIS.D <- summary(LIS)$coefficients[,,"DRUG"]

## Compile Stats for Bayesian Inference
MU.i.lis <- coef(LIS)
MU.0.lis <- colMeans( MU.i.lis )
# TAU <- as.numeric( VarCorr(LIS)[1:2,"StdDev"] ) ; names(TAU) <- names(MU.0)
SIG.i.lis <- aggregate( x=data.frame(resid(LIS)), by=list(TAB$IID), sd )

POST.lis <- merge( MU.i.lis, SIG.i.lis, by.x="row.names",by.y="Group.1", all=T )
colnames(POST.lis) <- c("ID","INT","PMD","PSD")

 # Plot Posterior Distributions for Individuals
XX <- seq( -10,10,.05 )
COLS.list <- c("black","slateblue2","steelblue2","springgreen2","gold1","chocolate2","firebrick2")
COLS <- colorRampPalette(COLS.list)(nrow(POST.lis))
PROBS.lis <- numeric( nrow(POST.lis) )
plot( 0,0,type="n",xlim=c(-5,5),ylim=c(0,2),xlab="Posterior Mean Difference",ylab="" )
abline( v=0,lty=2,lwd=2 )
N <- 20 ; WHICH <- sample(1:nrow(POST.lis),N)
for ( i in 1:nrow(POST.lis) ) {
# for ( w in 1:N ) {
   # i <- WHICH[w]
   points( XX, dnorm( XX,POST.lis[i,"PMD"],POST.lis[i,"PSD"]), type="l",col=COLS[i] )
   PROBS.lis[i] <- pnorm( 0, POST.lis[i,"PMD"],POST.lis[i,"PSD"], lower.tail=T )
   text( 3,2-quantile(c(0,2),(w-1)/N), label=round(PROBS.lis[i],2),col=COLS[i] )
}

##############################################################
## LM MODELS #################################################
##############################################################

## Single Linear Regression
LM <- lm( DAS ~ DRUG, data=TAB )

## Compile Stats for Bayesian Inference
MU.i.lm <- coef(LM)
# TAU <- as.numeric( VarCorr(LIS)[1:2,"StdDev"] ) ; names(TAU) <- names(MU.0)
SIG.i.lm <- confint(LM)

POST.lm <- merge( MU.i.lm, SIG.i.lm, by.x="row.names",by.y="Group.1", all=T )
colnames(POST.lm) <- c("ID","INT","PMD","PSD")

 # Plot Posterior Distributions for Individuals
XX <- seq( -10,10,.05 )
COLS.lmt <- c("black","slateblue2","steelblue2","springgreen2","gold1","chocolate2","firebrick2")
COLS <- colorRampPalette(COLS.lmt)(nrow(POST.lm))
PROBS.lm <- numeric( nrow(POST.lm) )
plot( 0,0,type="n",xlim=c(-5,5),ylim=c(0,2),xlab="Posterior Mean Difference",ylab="" )
abline( v=0,lty=2,lwd=2 )
N <- 20 ; WHICH <- sample(1:nrow(POST.lm),N)
for ( i in 1:nrow(POST.lm) ) {
# for ( w in 1:N ) {
   # i <- WHICH[w]
   points( XX, dnorm( XX,POST.lm[i,"PMD"],POST.lm[i,"PSD"]), type="l",col=COLS[i] )
   PROBS.lm[i] <- pnorm( 0, POST.lm[i,"PMD"],POST.lm[i,"PSD"], lower.tail=T )
   text( 3,2-quantile(c(0,2),(w-1)/N), label=round(PROBS.lm[i],2),col=COLS[i] )
}



##############################################################
## COMPARE MODELS ############################################
##############################################################








COLS <- c("deepskyblue2","chocolate2","chartreuse3")
XLIM <- range( MU.i[,1],MU.i.lis[,1] )
YLIM <- range( MU.i[,2],MU.i.lis[,2] )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="DRUG" )
arrows( MU.i.lis[,1], MU.i.lis[,2], MU.i[,1], MU.i[,2], col="grey50",length=.1, lwd=rowMeans(cbind(SIG.i.lis[,2],SIG.i[,2])) )
points( MU.i.lis[,1],MU.i.lis[,2], col=COLS[1], pch=20, cex=SIG.i.lis[,2] )
points( MU.i[,1],MU.i[,2], col=COLS[2], pch=20, cex=SIG.i[,2] )
abline( v=mean(MU.i.lis[,1]),col=COLS[1],lty=2,lwd=2 ) ; abline( h=mean(MU.i.lis[,2]),col=COLS[1],lty=2,lwd=2 )
abline( v=mean(MU.i[,1]),col=COLS[2],lty=2,lwd=1 ) ; abline( h=mean(MU.i[,2]),col=COLS[2],lty=2,lwd=1 )
points( MU.i.lm[1],MU.i.lm[2], pch=10,col=COLS[3],cex=2,lwd=2 )
legend( "bottomleft",legend=c("LIS","LME","LM","Hi_Var","Lo_Var"),col=c(COLS,"black","black"),pch=c(20,20,10,20,20),pt.cex=c(1,1,2,2,.5))
# hist( MU.i.lis[,1], freq=F,add=T, density=20, angle=45, col=COLS[1])
# hist( MU.i[,1], freq=F,add=T, density=20, angle=-45, col=COLS[2])






COLS <- c("deepskyblue2","chocolate2")
boxplot( DAS ~ DRUG + IID, data=TAB, col=COLS,las=2 )
MN.id.D0 <- aggregate( DAS ~ IID, data=TAB, subset=DRUG==0, mean )
MN.id.D1 <- aggregate( DAS ~ IID, data=TAB, subset=DRUG==1, mean )
hist( MN.id.D0[,2], xlim=c(1,9),ylim=c(0,30), density=20, angle=45, col=COLS[1] )
hist( MN.id.D1[,2], xlim=c(1,9),ylim=c(0,30), density=20, angle=-45, col=COLS[2], add=T )



summary(lm(DAS~DRUG*(WK+RF_ACPA),data=TAB))$coefficients
summary(lme(DAS~DRUG*(WK+RF_ACPA),random=~DRUG|IID,data=TAB))$tTable
summary(lme(DAS~DRUG*(WK+RF_ACPA),random=~DRUG+WK|IID,data=TAB))$tTable









# LIM.I <- range( MU.i[,1],MU.i.lis[,1] )
# LIM.D <- range( MU.i[,2],MU.i.lis[,2] )
# par(mfrow=c(1,2))
# plot( MU.i.lis[,1], MU.i[,1], xlim=LIM.I, ylim=LIM.I, main="Intercept" ) ; abline(0,1)
# plot( MU.i.lis[,2], MU.i[,2], xlim=LIM.D, ylim=LIM.D, main="Drug" ) ; abline(0,1)
