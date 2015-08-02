## Do Bayesian Analysis w/ Mixed Model & Longitudinal Data ##
 # Use .../1997_Nof1_Trial_Analysis.pdf (Zucker, et al.) as Reference
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

library(nlme)
library(gplots)

##############################################################
## LOAD DATA #################################################
##############################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Mac Paths
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150520_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_LME_Bayes/",sep="" )
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Load Candidate Genotype Files
# TAB.PC <- merge( TAB.l, FUL[,c("ID_2",paste("PC",1:3,sep=""))], by.x="IID",by.y="ID_2")
# COMP.l <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", sep="\t",header=T)
# COMP <- COMP.l[ which(!duplicated(COMP.l$SNP)), ]
# RAW <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/TEST_CND.raw", sep="",header=T)

##############################################################
## FILTER DATA ###############################################
##############################################################

## Take out Patients who left within 4 weeks of getting DRUG
RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
RM.exit <- which( TAB.l$IID %in% RM.exit.id )
TAB <- TAB.l[-RM.exit,c(1:15,17)]

# ## Take out Patients from Golimumab Arm
# RM.gol.id <- as.character( FUL$ID_2[which(FUL$GRP=="G")] )
# RM.gol <- which( TAB$IID %in% RM.gol.id )
# TAB <- TAB[ -RM.gol, ]

# ## Take out WK==0
# TAB <- TAB[ which(TAB$WK!=0), ]

## Remove NA Values for DAS
TAB <- TAB[ which(!is.na(TAB$DAS)), ]
dim(TAB)

##############################################################
## FCT: MODEL, ANALYSE, PLOT #################################
##############################################################

## Create Function to Build Models, Analyse, & Plot
DO_IT <- function( TAB, file_tag ) {

   ##############################################################
   ## LMM MODELS ################################################
   	## When using ANOVA to compare models:
   	 # Greater (less negative) logLik is better model
   	 # Smaller (less positive) AIC/BIC is better model

   ## Mixed Effects Model
   LME <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB )

   ## Compile Stats for Bayesian Inference
    # mu(i) estimates
    # mu(0) estimate
    # tau estimate
    # sigma(i) estimates
   MU.0.lme <- fixef(LME)
   MU.i.lme <- coef(LME)
   TAU.lme <- as.numeric( VarCorr(LME)[1:2,"StdDev"] ) ; names(TAU.lme) <- names(MU.0.lme)
   SIG.i.lme <- aggregate( x=data.frame(resid(LME)), by=list(TAB$IID), sd )

   POST.lme <- merge( MU.i.lme, SIG.i.lme, by.x="row.names",by.y="Group.1", all=T )
   colnames(POST.lme) <- c("ID","INT","PMD","PSD")

   ##############################################################
   ## LIS MODELS ################################################

   ## Independent Linear Regression Models
   LIS <- lmList( DAS ~ DRUG | IID, data=TAB )
   LIS.I <- summary(LIS)$coefficients[,,"(Intercept)"]
   LIS.D <- summary(LIS)$coefficients[,,"DRUG"]

   ## Compile Stats for Bayesian Inference
   MU.i.lis <- coef(LIS)
   MU.0.lis <- colMeans( MU.i.lis )
   # TAU.lis <- as.numeric( VarCorr(LIS)[1:2,"Variance"] ) # as.numeric( VarCorr(LIS)[1:2,"StdDev"] ) ; names(TAU.lis) <- names(MU.0)
   # SIG.i.lis <- aggregate( x=data.frame(resid(LIS)), by=list(TAB$IID), sd )
   SIG.i.lis <- cbind( LIS.I[,"Std. Error"], LIS.D[,"Std. Error"] ) ; colnames(SIG.i.lis) <- c("(Intercept)","DRUG")

   POST.lis <- merge( MU.i.lis, SIG.i.lis[,"DRUG"], by="row.names", all=T ) # ,by.y="Group.1"
   colnames(POST.lis) <- c("ID","INT","PMD","PSD")

   ##### MANUAL LIS MODE #####
   MAN <- list()
   for ( samp in unique(TAB$IID) ) { MAN[[samp]] <- lm( DAS ~ DRUG, data=TAB, subset=IID==samp ) }
   MU.i.man <- matrix( unlist(lapply( MAN, coef )),byrow=T,ncol=2 ) ; rownames(MU.i.man) <- names(MAN) ; colnames(MU.i.man) <- colnames(MU.i.lme)

   ##############################################################
   ## LM MODELS #################################################

   ## Single Linear Regression
   LM <- lm( DAS ~ DRUG, data=TAB )

   ## Compile Stats for Bayesian Inference
   MU.0.lm <- coef(LM)
   SIG.0.lm <- summary(LM)$coefficients[,"Std. Error"]

   POST.lm <- c( MU.0.lm, sqrt(SIG.0.lm) )
   names(POST.lm) <- c("INT","PMD","INT_SD","PSD")

   ##############################################################
   ## LMI MODELS ################################################

   ## Linear Regression w/ DRUG*IID Interaction
   # LMI <- lm( DAS ~ DRUG*IID, data=TAB )

   # MU.i.lmi <- coef(LMI)[c("(Intercept)","DRUG")]
   # MU.i.lmi <- rbind( MU.i.lmi, t(MU.i.lmi+t(cbind( coef(LMI)[grep("^IID",names(coef(LMI)))], coef(LMI)[grep("DRUG:",names(coef(LMI)))] ))) )
   # colnames(MU.i.lmi) <- c("(Intercept)","DRUG") ; rownames(MU.i.lmi) <- gsub("IID","",rownames(MU.i.lmi)) ; rownames(MU.i.lmi)[1] <- setdiff( rownames(MU.i.lis),rownames(MU.i.lmi) )
   # TAU.lmi <- apply( MU.i.lmi, 2, var) ; # names(TAU.lmi) <- names(MU.0.lmi)
   # # SIG.i.lmi <- confint(LMI)
   # SIG.i.lmi <- summary(LMI)$coefficients[,"Std. Error"]
   # SIG.0.lmi <- summary(LMI)$coefficients["DRUG","Std. Error"]

   ##############################################################
   ## CALCULATE POSTERIOR PROBABILITIES #########################
   PROBS.0.lis <- PROBS.0.lme <- PROBS.1.lis <- PROBS.1.lme <- numeric( nrow(POST.lis) )
   for ( i in 1:nrow(POST.lis) ) {
      PROBS.0.lme[i] <- pnorm( 0, POST.lme[i,"PMD"],POST.lme[i,"PSD"], lower.tail=T )
      PROBS.1.lme[i] <- pnorm( -1, POST.lme[i,"PMD"],POST.lme[i,"PSD"], lower.tail=T )
      PROBS.0.lis[i] <- pnorm( 0, POST.lis[i,"PMD"],POST.lis[i,"PSD"], lower.tail=T )
      PROBS.1.lis[i] <- pnorm( -1, POST.lis[i,"PMD"],POST.lis[i,"PSD"], lower.tail=T )
   }
   PROBS.0.lm <- pnorm( 0, POST.lm["PMD"],POST.lm["PSD"], lower.tail=T )
   PROBS.1.lm <- pnorm( -1, POST.lm["PMD"],POST.lm["PSD"], lower.tail=T )

   ##############################################################
   ## PLOT MODELS ###############################################

   png( paste(PathToPlot,"SHRINK_",file_tag,".png",sep=""), height=1000,width=2000, pointsize=30 )
   par(mfrow=c(1,2))

   ## Plot Shrinkage of Coefficients w/ Mixed Effects Model
   COLS <- c("deepskyblue2","chocolate2","slateblue3")
   XLIM <- c(1,9) # range( MU.i.lm[,1],MU.i.lis[,1] )
   YLIM <- c(-5,3) # range( MU.i.lm[,2],MU.i.lis[,2] )
   plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="DRUG", main="Individual vs Mixed Effects Model Coefficients" )
   arrows( MU.i.lis[,1], MU.i.lis[,2], MU.i.lme[,1], MU.i.lme[,2], col="grey50",length=.2, lwd=2*rowMeans(cbind(SIG.i.lis[,2],SIG.i.lme[,2])) )
   points( MU.i.lis[,1],MU.i.lis[,2], col=COLS[1], pch=20, cex=SIG.i.lis[,2] )
   points( MU.i.lme[,1],MU.i.lme[,2], col=COLS[2], pch=20, cex=SIG.i.lme[,2] )
   abline( v=mean(MU.i.lis[,1]),col=COLS[1],lty=2,lwd=2 ) ; abline( h=mean(MU.i.lis[,2]),col=COLS[1],lty=2,lwd=2 )
   abline( v=mean(MU.i.lme[,1]),col=COLS[2],lty=2,lwd=1 ) ; abline( h=mean(MU.i.lme[,2]),col=COLS[2],lty=2,lwd=1 )
   points( MU.0.lm[1],MU.0.lm[2], pch=10,col=COLS[3],cex=2,lwd=2 )
   legend( "bottomleft",legend=c("LIS","LME","LM","Hi_Var","Lo_Var"),col=c(COLS,"black","black"),pch=c(20,20,10,20,20),pt.cex=c(1,1,2,2,.5))

   ## Plot Probabilities for LME vs LIS
   plot( 0,0,type="n",xlim=c(0,1),ylim=c(0,1), main="Probability of Treatment Response",xlab="Individual Model",ylab="Mixed Effects Model" )
   abline(0,1,lty=1,lwd=1) ; abline(h=seq(0,1,.2),lty=3,col="grey50",lwd=1) ; abline(v=seq(0,1,.2),lty=3,col="grey50",lwd=1)
    # Pr<0
   points( PROBS.0.lis, PROBS.0.lme, col="firebrick1",pch=1,lwd=2,cex=rowMeans(cbind(SIG.i.lis[,2],SIG.i.lme[,2])) )
   abline(lm( PROBS.0.lme ~ PROBS.0.lis ), col="firebrick4",lwd=2,lty=2 )
   # hist( PROBS.0.lis, freq=F,col="firebrick1",density=20,angle=45, add=T )
   points( PROBS.1.lis, PROBS.1.lme, col="chartreuse1",pch=1,lwd=2,cex=rowMeans(cbind(SIG.i.lis[,2],SIG.i.lme[,2])) )
   abline(lm( PROBS.1.lme ~ PROBS.1.lis ), col="chartreuse4",lwd=2,lty=2 )
   # hist( PROBS.1.lis, freq=F,col="chartreuse1",density=20,angle=-45, add=T )
   dev.off() # Close png Output

   # COLS <- c("deepskyblue2","chocolate2")
   # boxplot( DAS ~ DRUG + IID, data=TAB, col=COLS,las=2 )
   # MN.id.D0 <- aggregate( DAS ~ IID, data=TAB, subset=DRUG==0, mean )
   # MN.id.D1 <- aggregate( DAS ~ IID, data=TAB, subset=DRUG==1, mean )
   # hist( MN.id.D0[,2], xlim=c(1,9),ylim=c(0,30), density=20, angle=45, col=COLS[1] )
   # hist( MN.id.D1[,2], xlim=c(1,9),ylim=c(0,30), density=20, angle=-45, col=COLS[2], add=T )

   # summary(lm(DAS~DRUG*(WK+RF_ACPA),data=TAB))$coefficients
   # summary(lme(DAS~DRUG*(WK+RF_ACPA),random=~DRUG|IID,data=TAB))$tTable
   # summary(lme(DAS~DRUG*(WK+RF_ACPA),random=~DRUG+WK|IID,data=TAB))$tTable
      
   ## Calculate Shrinkage
   SHRINK.int <- MU.i.lis[,1] - MU.i.lme[,1]
   SHRINK.drug <- MU.i.lis[,2] - MU.i.lme[,2]

   ## Output Compiled Results
   COMPILE <- list( SHRINK.int, SHRINK.drug, PROBS.0.lme, PROBS.0.lis, PROBS.1.lme, PROBS.1.lis, POST.lme, POST.lis )
   names(COMPILE) <- c("SH.i","SH.d","PROBS.0.lme","PROBS.0.lis","PROBS.1.lme","PROBS.1.lis","POST.lme","POST.lis")
   return(COMPILE)

} # Close Function

TEST <- DO_IT(TAB, "Full")

##############################################################
## OPTIONAL: DOWN SAMPLE TO ASSESS SHRINKAGE #################
##############################################################
 # Down sample to fraction of observations per person
 # Output Amount of Shrinkage (Int & Drug) at each level of down-sampling
Obs_per_Person <- table( TAB$IID )
barplot( table( Obs_per_Person ))
RM.exit.2.id <- names(Obs_per_Person)[which(Obs_per_Person<15)]
RM.exit.2 <- which(TAB$IID %in% RM.exit.2.id)
TAB.2 <- TAB[ -RM.exit.2, ] # Should be 120 patients, each with 15 measurements

## Number of Observations to Down Sample to
N.pats <- length(unique( TAB.2$IID ))
N.rep <- 20
N.obs <- c( 15, rep(12,N.rep), rep(9,N.rep), rep(6,N.rep), rep(3,N.rep) )
SHRINK <- list()
for ( o in 1:length(N.obs) ) {
   obs <- N.obs[o]
   tag <- paste( "It",o,"_",obs,"obs",sep="")
   print(tag)
   ## Down Sample Table
   which_rows <- sort(c(apply( matrix(1:nrow(TAB.2),byrow=T,ncol=15), 1, function(x) c(sample(x[1:5],obs/3),sample(x[6:10],obs/3),sample(x[11:15],obs/3)) )))
   TEMP_TAB <- TAB.2[ which_rows, ]
   SHRINK[[tag]] <- DO_IT(TEMP_TAB, tag )
}

## Compile Shrinkage Results
SHRINK.int <- unlist(lapply( SHRINK, function(x) mean(abs(x$SH.i)) ))
SHRINK.drug <- unlist(lapply( SHRINK, function(x) mean(abs(x$SH.d)) ))
SHRINK.pyth <- sqrt( SHRINK.int^2 + SHRINK.drug^2 )
COLS.box <- c("slateblue3","cadetblue2","tomato2")
YLIM <- c(0,1)
png( paste(PathToPlot,"SHRINK_2-Shrink_by_Nobs.png",sep=""), height=700,width=2100, pointsize=26 )
par(mfrow=c(1,3))
boxplot( SHRINK.int ~ factor(N.obs), ylim=YLIM,col=COLS.box[1], main="Intercept Coef Shrinkage",xlab="# Observations per Patient",ylab="Mean( Coef[LIS] - Coef[LME] )" )
points( SHRINK.int ~ factor(N.obs), pch=20 )
boxplot( SHRINK.drug ~ factor(N.obs), ylim=YLIM,col=COLS.box[2], main="Drug Coef Shrinkage",xlab="# Observations per Patient",ylab="Mean( Coef[LIS] - Coef[LME] )" )
points( SHRINK.drug ~ factor(N.obs), pch=20 )
boxplot( SHRINK.pyth ~ factor(N.obs), ylim=YLIM,col=COLS.box[3], main="Pythagorean Coef Shrinkage",xlab="# Observations per Patient",ylab="Mean( Coef[LIS] - Coef[LME] )" )
points( SHRINK.pyth ~ factor(N.obs), pch=20 )
dev.off()

## See if within Patient Variance correlates w/ size of Correction
SD.lis <- lapply( SHRINK, function(x) x$POST.lis[,"PSD"] )
SD.lme <- lapply( SHRINK, function(x) x$POST.lme[,"PSD"] )
SHR.drug <- lapply( SHRINK, function(x) abs(x$SH.d) )
SHR.int <- lapply( SHRINK, function(x) abs(x$SH.i) )

VARR <- data.frame( SD.lis=unlist(SD.lis), SD.lme=unlist(SD.lme), SHR.drug=unlist(SHR.drug), SHR.int=unlist(SHR.int) )
VARR.obs <- as.numeric(sapply(strsplit( sapply(strsplit( rownames(VARR),"_" ),"[",2), "obs"),"[",1))
VARR <- data.frame( VARR, OBS=VARR.obs )
VARR <- data.frame( VARR, SD=(VARR$SD.lis+VARR$SD.lme)/2 )
VARR <- data.frame( VARR, SHR.pyth=sqrt(VARR$SHR.drug^2+VARR$SHR.drug^2) )
MOD.drug <- lm( SHR.drug ~ SD*factor(OBS), data=VARR )
MOD.int <- lm( SHR.int ~ SD*factor(OBS), data=VARR )
summary(MOD.drug) ; summary(MOD.int)
COLS.list <- c("slateblue1","steelblue1","springgreen1","gold1","chocolate1","firebrick1")[c(1:3,5,6)]
COLS.3.list <- c("slateblue3","steelblue3","springgreen3","gold3","chocolate3","firebrick3")[c(1:3,5,6)]
ORDER <- sample( 1:nrow(VARR) )
png( paste(PathToPlot,"SHRINK_3-Shrink_vs_Var.png",sep=""), height=800,width=2400, pointsize=30 )
par(mfrow=c(1,3))
 # Intercept
plot( VARR$SHR.int[ORDER] ~ VARR$SD[ORDER], col=COLS.list[factor(VARR$OBS[ORDER])], main="Intercept Coef Shrinkage vs Variance",xlab="Within Patient Variance",ylab="Coef[LIS] - Coef[LME]", )
abline(lm(SHR.int~SD,data=VARR,subset=OBS==3), col=COLS.3.list[1],lwd=2 )
abline(lm(SHR.int~SD,data=VARR,subset=OBS==6), col=COLS.3.list[2],lwd=2 )
abline(lm(SHR.int~SD,data=VARR,subset=OBS==9), col=COLS.3.list[3],lwd=2 )
abline(lm(SHR.int~SD,data=VARR,subset=OBS==12), col=COLS.3.list[4],lwd=2 )
abline(lm(SHR.int~SD,data=VARR,subset=OBS==15), col=COLS.3.list[5],lwd=2 )
 # Drug
plot( VARR$SHR.drug[ORDER] ~ VARR$SD[ORDER], col=COLS.list[factor(VARR$OBS[ORDER])], main="Drug Coef Shrinkage vs Variance",xlab="Within Patient Variance",ylab="Coef[LIS] - Coef[LME]", )
abline(lm(SHR.drug~SD,data=VARR,subset=OBS==3), col=COLS.3.list[1],lwd=2 )
abline(lm(SHR.drug~SD,data=VARR,subset=OBS==6), col=COLS.3.list[2],lwd=2 )
abline(lm(SHR.drug~SD,data=VARR,subset=OBS==9), col=COLS.3.list[3],lwd=2 )
abline(lm(SHR.drug~SD,data=VARR,subset=OBS==12), col=COLS.3.list[4],lwd=2 )
abline(lm(SHR.drug~SD,data=VARR,subset=OBS==15), col=COLS.3.list[5],lwd=2 )
 # Pyth
plot( VARR$SHR.pyth[ORDER] ~ VARR$SD[ORDER], col=COLS.list[factor(VARR$OBS[ORDER])], main="Pyth Coef Shrinkage vs Variance",xlab="Within Patient Variance",ylab="Coef[LIS] - Coef[LME]", )
abline(lm(SHR.pyth~SD,data=VARR,subset=OBS==3), col=COLS.3.list[1],lwd=2 )
abline(lm(SHR.pyth~SD,data=VARR,subset=OBS==6), col=COLS.3.list[2],lwd=2 )
abline(lm(SHR.pyth~SD,data=VARR,subset=OBS==9), col=COLS.3.list[3],lwd=2 )
abline(lm(SHR.pyth~SD,data=VARR,subset=OBS==12), col=COLS.3.list[4],lwd=2 )
abline(lm(SHR.pyth~SD,data=VARR,subset=OBS==15), col=COLS.3.list[5],lwd=2 )
legend( "topleft", fill=COLS.list,legend=rev(unique(VARR$OBS)),title="# Obs/Patient" )
dev.off()

##############################################################
## PLOT POSTERIORS of FULL MODEL #############################
##############################################################

OUT <- DO_IT(TAB, "Full")
POST.lis <- OUT$POST.lis
POST.lme <- OUT$POST.lme
# PROBS.lis <- OUT$PR.lis
# PROBS.lme <- OUT$PR.lme

## Plot Posterior Distributions for Individuals
 # Parameters
N <- 8 ; WHICH <- sample(1:nrow(POST.lis),N) ; rand <- sample(1:100,1)
XX <- seq( -10,10,.01 )
COLS <- c("deepskyblue2","chocolate2","slateblue3")
COLS.lis.list <- c("black",COLS[1])
COLS.lis <- colorRampPalette(COLS.lis.list)(2+N)[1:N+2]
COLS.lme.list <- c("black",COLS[2])
COLS.lme <- colorRampPalette(COLS.lme.list)(2+N)[1:N+2]
PROBS.0.lis <- PROBS.0.lme <- PROBS.1.lis <- PROBS.1.lme <- numeric( nrow(POST.lis) )
 # Open Plot
png( paste(PathToPlot,"SHRINK_4-PostProb_Distribs.",rand,".png",sep=""), height=1000,width=1500, pointsize=28 )
plot( 0,0,type="n",xlim=c(-5,5),ylim=c(0,2),xlab="Posterior Mean Difference",ylab="" )
abline( v=c(0,-1),lty=2,lwd=2,col=c("firebrick1","chartreuse1") )
abline( v=3.8 )
text( 3, 2, label="Pr < 0",col="firebrick1" )
text( 4.5, 2, label="Pr < -1",col="chartreuse1" )
w <- 0
for ( i in 1:nrow(POST.lis) ) {
   PROBS.0.lme[i] <- pnorm( 0, POST.lme[i,"PMD"],POST.lme[i,"PSD"], lower.tail=T )
   PROBS.1.lme[i] <- pnorm( -1, POST.lme[i,"PMD"],POST.lme[i,"PSD"], lower.tail=T )
   PROBS.0.lis[i] <- pnorm( 0, POST.lis[i,"PMD"],POST.lis[i,"PSD"], lower.tail=T )
   PROBS.1.lis[i] <- pnorm( -1, POST.lis[i,"PMD"],POST.lis[i,"PSD"], lower.tail=T )
   if ( i %in% WHICH ) { w <- w + 1
      # LIS
      color <- COLS.lis[w]
      points( XX, dnorm( XX,POST.lis[i,"PMD"],POST.lis[i,"PSD"]), type="l",col=color,lwd=2 )
      arrows( POST.lis[i,"PMD"],1.85,POST.lis[i,"PMD"],1.75, col=color,length=.1,lwd=2 )
      text( 2.7,quantile(c(.1,1.8),(w)/N), label=round(PROBS.0.lis[i],2),col=color )
      text( 4.2,quantile(c(.1,1.8),(w)/N), label=round(PROBS.1.lis[i],2),col=color )
      # LME
      color <- COLS.lme[w]
      points( XX, dnorm( XX,POST.lme[i,"PMD"],POST.lme[i,"PSD"]), type="l",col=color,lwd=2 )
      arrows( POST.lme[i,"PMD"],1.95,POST.lme[i,"PMD"],1.85, col=color,length=.1,lwd=2 )
      text( 3.3,quantile(c(.1,1.8),(w)/N), label=round(PROBS.0.lme[i],2),col=color )
      text( 4.8,quantile(c(.1,1.8),(w)/N), label=round(PROBS.1.lme[i],2),col=color )
   }
}
legend( "topleft", fill=COLS[1:2], legend=c("LIS","LME") )
dev.off()















# o <- sample( 1:length(SHRINK), 1, prob=1:length(SHRINK) ) ; print(o)
# plot( colMeans(rbind(SD.lis[[o]],SD.lme[[o]])), SHR.int[[o]], pch="+", col=COLS.box[1], main=names(SHRINK)[o] ) ; abline(lm( SHR.int[[o]]~colMeans(rbind(SD.lis[[o]],SD.lme[[o]])) ) ,lty=2,col=COLS.box[1] )
# points( colMeans(rbind(SD.lis[[o]],SD.lme[[o]])), SHR.drug[[o]], col=COLS.box[2] ) ; abline(lm( SHR.drug[[o]]~colMeans(rbind(SD.lis[[o]],SD.lme[[o]])) ) ,lty=2,col=COLS.box[2] )

# LIM.I <- range( MU.i.lm[,1],MU.i.lis[,1] )
# LIM.D <- range( MU.i.lm[,2],MU.i.lis[,2] )
# par(mfrow=c(1,2))
# plot( MU.i.lis[,1], MU.i.lm[,1], xlim=LIM.I, ylim=LIM.I, main="Intercept" ) ; abline(0,1)
# plot( MU.i.lis[,2], MU.i.lm[,2], xlim=LIM.D, ylim=LIM.D, main="Drug" ) ; abline(0,1)






