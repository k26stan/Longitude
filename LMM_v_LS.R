## Run LMM Models on Time-Series Data ##
## Janssen Data using Multiple-Measures ##
## July 29, 2014 ##
## Kristopher Standish ##

library(nlme)
library(gplots)
##############################################################
## LOAD DATA #################################################
##############################################################

## Set Date
DATE <- "20150513"

# ## TSCC Paths
# PathToData <- "/projects/janssen/clinical/"
# PathToSave <- "/projects/janssen/clinical/Plots/20140522/"

## Mac Paths
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150506_Resp_v_Time.txt"
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150512_Resp_v_Time.txt"
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,sep="" )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T )
FUL <- read.table( paste(PathToFT,"20141229_Full_Table.txt",sep=""),sep="\t",header=T)
TAB.l.2 <- merge( TAB.l, FUL[,c("ID_2",paste("PC",1:3,sep=""))], by.x="IID",by.y="ID_2")
# Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Load Candidate Genotype Files
COMP.l <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", sep="\t",header=T)
COMP <- COMP.l[ which(!duplicated(COMP.l$SNP)), ]
RAW <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/TEST_CND.raw", sep="",header=T)

##############################################################
## FILTER DATA ###############################################
##############################################################

## Remove some patients & timepoints
 # Patients who left before getting DRUG
RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
RM.exit <- which( TAB.l$IID %in% RM.exit.id )
 # Timepoints where DAS==NA
RM.na <- which(is.na( TAB.l$DAS ))

## Remove
RM <- union( RM.exit, RM.na )
TAB <- TAB.l[-RM,c(1:15,17)]

##############################################################
## CREATE PLOTTING FUNCTIONS #################################
##############################################################

PLOT_MODS <- function( MODELS ) {
	N.mod <- length(MODELS)
	Names.mod <- names(MODELS)
	Nums.mod <- gsub("[^0-9]", "",Names.mod)
	Swhich <- which(!duplicated(Nums.mod))
	COMP <- array( ,c(4,N.mod))
	rownames(COMP) <- c("logLik","AIC","BIC","SD")
	colnames(COMP) <- Names.mod
	for ( i in 1:N.mod ) {
		COMP["logLik",] <- unlist(lapply( MODELS, logLik ))
		COMP["AIC",] <- unlist(lapply( MODELS, AIC ))
		COMP["BIC",] <- unlist(lapply( MODELS, BIC ))
		COMP["SD",] <- unlist(lapply( MODELS, function(x) sd(resid(x)) ))
	}
	## Plot
	COLS <- c("dodgerblue2","tomato2","slateblue3","chartreuse2")
	png( paste(PathToPlot,"1-LME_ModAssess.png",sep="/"), height=800,width=2400,pointsize=30 )
	par(mfrow=c(1,3))
	 # logLik
	plot( 1:N.mod, COMP["logLik",], type="o", col=COLS[1], pch=20, main="log(Likelihood) of Models",ylab="log(Lik)",xlab="Model",xaxt="n" )
	axis( 1, at=1:N.mod, label=colnames(COMP), las=2 )
	abline( v=which.max(COMP["logLik",]),lty=2,col=COLS[1],lwd=2 )
	abline( v=Swhich-.5,lty=1,col="grey50",lwd=1 )
	 # AIC/BIC
	YLIM <- range( COMP[c("AIC","BIC"),] )
	plot( 1:N.mod, COMP["AIC",], type="o", ylim=YLIM, col=COLS[2], pch=20, main="Information Criteria of Models",ylab="Information Criteria",xlab="Model",xaxt="n" )
	points( 1:N.mod, COMP["BIC",], type="o", col=COLS[3], pch=20 )
	axis( 1, at=1:N.mod, label=colnames(COMP), las=2 )
	abline( v=which.min(COMP["AIC",]),lty=2,col=COLS[2],lwd=2 )
	abline( v=which.min(COMP["BIC",]),lty=2,col=COLS[3],lwd=2 )
	abline( v=Swhich-.5,lty=1,col="grey50",lwd=1 )
	 # SD Res
	plot( 1:N.mod, COMP["SD",], type="o", col=COLS[4], pch=20, main="SD of Residuals of Models",ylab="SD Resids",xlab="Model",xaxt="n" )
	axis( 1, at=1:N.mod, label=colnames(COMP), las=2 )
	abline( v=which.min(COMP["SD",]),lty=2,col=COLS[4],lwd=2 )
	abline( v=Swhich-.5,lty=1,col="grey50",lwd=1 )
	dev.off()
	## Return Stats
	return(COMP)
}
# TEMP <- PLOT_MODS(LME)

HEAT_MODS <- function( MODELS ) {

	N.mod <- length(MODELS)
	Names.mod <- names(MODELS)
	Nums.mod <- gsub("[^0-9]", "",Names.mod)
	Swhich <- which(!duplicated(Nums.mod))

	## Specify Number and Names of Models
	MODELS.N.models <- length(MODELS)
	MODELS.Names.models <- names(MODELS)

	## Determine which models have DRUG interactions
	MODELS.drugint <- as.numeric(unlist( lapply( MODELS, function(x) any(grepl( "DRUG:|:DRUG", names(coef(x)) ))) ))

	################################################
	## Compile Terms and Significance for each model

	## Get all Terms used in any model
	MODELS.terms <- Reduce( union, lapply( MODELS, function(x) colnames(coef(x)) ) )
	# MODELS.terms <- sort( MODELS.terms )
	MODELS.N.terms <- length(MODELS.terms)

	## Determine which Terms are in each model
	MODELS.ARR.terms <- array( 0, c(MODELS.N.terms,MODELS.N.models) )
	MODELS.ARR.terms.p <- array( 1, c(MODELS.N.terms,MODELS.N.models) )
	colnames(MODELS.ARR.terms) <- colnames(MODELS.ARR.terms.p) <- MODELS.Names.models
	rownames(MODELS.ARR.terms) <- rownames(MODELS.ARR.terms.p) <- MODELS.terms
	for ( m in 1:length(MODELS) ) {
		mod <- names(MODELS)[m]
		which.in <- which( MODELS.terms %in% colnames(coef(MODELS[[mod]])) )
		which.terms <- MODELS.terms[which.in]
		MODELS.ARR.terms[which.in,mod] <- 1
		which.terms.2 <- which.terms[ which( which.terms %in% rownames(summary(MODELS[[mod]])$tTable) ) ]
		which.in.2 <- which( MODELS.terms %in% which.terms.2 )
		MODELS.ARR.terms.p[which.in.2,mod] <- summary(MODELS[[mod]])$tTable[which.terms.2,"p-value"]
	}
	 # Plot P-Values for each Term
	COLS.list <- c("firebrick1","chocolate1","gold1","springgreen2","steelblue2","slateblue3","black")
	COLS <- c( "firebrick1",colorRampPalette(COLS.list[3:7])(19) ) # colorRampPalette(COLS.list)(20)
	ROW.COLS <- rep("deepskyblue2",MODELS.N.terms) ; ROW.COLS[grep("DRUG",MODELS.terms)] <- "chartreuse2"
	png( paste(PathToPlot,"1-LME_ModBuild.png",sep="/"), height=1200,width=1600,pointsize=30 )
	heatmap.2( MODELS.ARR.terms.p, main="Inclusion/Significance of Terms in Models",xlab="Model",ylab="Term",trace="none",scale="none",Colv=F,Rowv=F,dendrogram="none",col=COLS,margin=c(6,11),lhei=c(2,9),lwid=c(1,6),RowSideColors=ROW.COLS, colsep=Swhich-1 )
	dev.off()
	## Return P-Vals
	return( MODELS.ARR.terms.p )
}
# TEMP <- HEAT_MODS(LME)

##############################################################
## LMM MODELS ################################################
##############################################################

## When using ANOVA to compare models:
 # Greater (less negative) logLik is better model
 # Smaller (less positive) AIC/BIC is better model

LME_MOD <- function( TAB ) {
	start_time <- proc.time()
	LME <- list()

	## Put together some Basic Models
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 1"))
	# Linear Mixed Model w/ only DRUG effect
	LME$M1a <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB )
	# Add Drug as Random Effect
	LME$M1b <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB )
	# Add Nested Random Effect DRUG %in% IID (IID/DRUG)
	LME$M1c <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID/DRUG, data=TAB )
	# Try removing Random Effects for Intercept
	LME$M1d <- lme( fixed = DAS ~ DRUG, random = ~ DRUG-1 | IID, data=TAB )
	# Use Nested Random Effect IID %in% COUN (COUN/IID)
	LME$M1e <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | COUN/IID, data=TAB )
	# Add Nested Random Effect DRUG %in% IID (IID/DRUG)
	LME$M1f <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID/DRUG, data=TAB )
	 # Choose Amongst Models
	anova(LME$M1a,LME$M1b,LME$M1c,LME$M1d,LME$M1e,LME$M1f)
	COMP <- PLOT_MODS( LME )

	## Include some Auto-correlation within Patient
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 2"))
	LME$M2a <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M2b <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M2c <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(form = ~1 | IID/DRUG) )
	LME$M2d <- lme( fixed = DAS ~ DRUG, random = ~ DRUG-1 | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M2e <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M2a,LME$M2b,LME$M2c,LME$M2d,LME$M2e)
	COMP <- PLOT_MODS( LME )

	## Include Auto-correlation for WK within Patient
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 3"))
	LME$M3a <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M3b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M3c <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(form = ~1 | IID/DRUG) )
	LME$M3d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG-1 | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M3e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M3a,LME$M3b,LME$M3c,LME$M3d,LME$M3e)
	COMP <- PLOT_MODS( LME )

	## Include Auto-correlation for WK within Patient
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 4"))
	LME$M4a <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID, data=TAB, correlation=corCAR1(form = ~WK | IID) )
	LME$M4b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(form = ~WK | IID) )
	LME$M4c <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(form = ~WK | IID/DRUG) )
	LME$M4d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG-1 | IID, data=TAB, correlation=corCAR1(form = ~WK | IID) )
	LME$M4e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(form = ~WK | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M4a,LME$M4b,LME$M4c,LME$M4d,LME$M4e)
	COMP <- PLOT_MODS( LME )

	## Include Auto-correlation for WK (but remove Intercept) within Patient
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 5"))
	LME$M5a <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID, data=TAB, correlation=corCAR1(form = ~WK-1 | IID) )
	LME$M5b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(form = ~WK-1 | IID) )
	LME$M5c <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(form = ~WK-1 | IID/DRUG) )
	LME$M5e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(form = ~WK-1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M5a,LME$M5b,LME$M5c,LME$M5e)
	COMP <- PLOT_MODS( LME )

	## Revert to Intercept only in Correlation
	## Include WK in Random Effects
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 6"))
	LME$M6a <- lme( fixed = DAS ~ DRUG*WK, random = ~ WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M6b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M6c <- lme( fixed = DAS ~ DRUG*WK, random = ~ WK | IID/DRUG, data=TAB, correlation=corCAR1(form = ~1 | IID/DRUG) )
	LME$M6e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M6a,LME$M6b,LME$M6c,LME$M6e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )

	## Include Fixed PLAC Effect
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 7"))
	LME$M7a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ 1 | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M7b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M7c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(form = ~1 | IID/DRUG) )
	LME$M7e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M7a,LME$M7b,LME$M7c,LME$M7e)
	COMP <- PLOT_MODS( LME )

	## Include Fixed RF_ACPA*DRUG Interaction
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 8"))
	LME$M8b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	# LME$M8c <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ WK | IID/DRUG, data=TAB, correlation=corCAR1(form = ~1 | IID/DRUG) )
	LME$M8e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M8b,LME$M8e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )

	## Try Fixed ACPA*DRUG Interaction Instead of RF_ACPA*DRUG
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 9"))
	LME$M9b <- lme( fixed = DAS ~ DRUG*(WK+I(ACPA=="Positive"))+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	# LME$M9c <- lme( fixed = DAS ~ DRUG*(WK+I(ACPA=="Positive"))+PLAC, random = ~ WK | IID/DRUG, data=TAB, correlation=corCAR1(form = ~1 | IID/DRUG) )
	LME$M9e <- lme( fixed = DAS ~ DRUG*(WK+I(ACPA=="Positive"))+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M9b,LME$M9e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )

	## Try Fixed RF*DRUG Interaction Instead of RF_ACPA*DRUG
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 10"))
	LME$M10b <- lme( fixed = DAS ~ DRUG*(WK+RF)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M10e <- lme( fixed = DAS ~ DRUG*(WK+RF)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M10b,LME$M10e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )

	## Add SEX*DRUG Interaction
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 11"))
	LME$M11b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+SEX)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M11e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+SEX)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M11b,LME$M11e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )
	HEAT <- HEAT_MODS( LME )

	## Add BMI*DRUG Interaction
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 12"))
	LME$M12b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+BMI)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M12e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+BMI)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M12b,LME$M12e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )
	HEAT <- HEAT_MODS( LME )

	## Add AGE*DRUG Interaction
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 13"))
	LME$M13b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+AGE)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M13e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+AGE)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M13b,LME$M13e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )
	HEAT <- HEAT_MODS( LME )

	## Add DIS_DUR*DRUG Interaction
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 14"))
	LME$M14b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+DIS_DUR)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M14e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+DIS_DUR)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M14b,LME$M14e)
	COMP <- PLOT_MODS( LME[6:length(LME)] )
	COMP <- PLOT_MODS( LME[25:length(LME)] )
	HEAT <- HEAT_MODS( LME )

	# ## Add DIS_DUR*DRUG Interaction
	# print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 15"))
	# LME$M15b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+DIS_DUR)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	# LME$M15e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA+DIS_DUR)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	#  # Choose Amongst Models
	# anova(LME$M15b,LME$M15e)
	# COMP <- PLOT_MODS( LME[6:length(LME)] )
	# COMP <- PLOT_MODS( LME[25:length(LME)] )
	# HEAT <- HEAT_MODS( LME )
	return(LME)

}
LME <- LME_MOD(TAB)

##############################################################
## LMM MODELS ################################################
##############################################################

## Estimate Top Models using lsList
LIS_MOD <- function( TAB ) {

	LIS <- list()

	## Put together some Basic Models
	# Linear Mixed Model w/ only DRUG effect
	LIS$M1b <- lmList( DAS ~ DRUG | IID, data=TAB )
	# Add Nested Random Effect DRUG %in% IID (IID/DRUG)
	LIS$M1c <- lmList( DAS ~ DRUG | IID/DRUG, data=TAB )
	# Use Nested Random Effect IID %in% COUN (COUN/IID)
	LIS$M1e <- lmList( DAS ~ DRUG | COUN/IID, data=TAB )

	## Include DRUG*WK Interaction
	LIS$M2b <- lmList( DAS ~ DRUG+WK | IID, data=TAB )
	LIS$M2c <- lmList( DAS ~ DRUG+WK | IID/DRUG, data=TAB )
	LIS$M2e <- lmList( DAS ~ DRUG+WK | COUN/IID, data=TAB )

	## Include DRUG*WK Interaction
	LIS$M3b <- lmList( DAS ~ DRUG*WK | IID, data=TAB )
	LIS$M3c <- lmList( DAS ~ DRUG*WK | IID/DRUG, data=TAB )
	LIS$M3e <- lmList( DAS ~ DRUG*WK | COUN/IID, data=TAB )

	## Include Fixed PLAC Effect
	LIS$M6b <- lmList( DAS ~ DRUG*WK+PLAC | IID, data=TAB )
	LIS$M6c <- lmList( DAS ~ DRUG*WK+PLAC | IID/DRUG, data=TAB )
	LIS$M6e <- lmList( DAS ~ DRUG*WK+PLAC | COUN/IID, data=TAB )

	## Include Fixed PLAC*WK Interaction
	LIS$M7b <- lmList( DAS ~ DRUG*WK+PLAC | IID, data=TAB )
	LIS$M7c <- lmList( DAS ~ DRUG*WK+PLAC | IID/DRUG, data=TAB )
	LIS$M7e <- lmList( DAS ~ DRUG*WK+PLAC | COUN/IID, data=TAB )

	## Return Models
	return(LIS)
}
LIS <- LIS_MOD(TAB)

##############################################################
## ESTIMATE COEFFICIENTS #####################################
##############################################################

DIST_PLOT <- function( CF.lis, CF.lme ) {
	N.vars <- ncol(CF.lis)
	VARS <- colnames(CF.lis)
	# COLS.list <- c("firebrick2","chocolate2","gold1","springgreen2","steelblue2","slateblue3")
	# COLS <- sample( colorRampPalette(COLS.list)(N.vars) )
	# COL_SPLIT <- function(COL1,COL2) { colorRampPalette(c(COL1,COL2))(3)[2] }
	COLS <- c("chocolate1","dodgerblue2","chartreuse2")
	COLS.4 <- c("chocolate4","dodgerblue4","chartreuse4")
	par(mfrow=c(N.vars,N.vars+1))
	for ( row in 1:N.vars ) {
		for ( col in 1:(N.vars+1) ) {
			if ( col <= N.vars ) { 
				## Plot
				if ( row > col ) {
					XLIM <- range( CF.lis[,col],CF.lme[,col], na.rm=T )
					YLIM <- range( CF.lis[,row],CF.lme[,row], na.rm=T )
					plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xlab=VARS[col],ylab=VARS[row],main=paste(VARS[col],"vs",VARS[row]) )
					abline( v=seq(floor(XLIM[1]),ceiling(XLIM[2]),10^floor(log10(diff(XLIM)))),lty=2,col="grey50" )
					abline( h=seq(floor(YLIM[1]),ceiling(YLIM[2]),10^floor(log10(diff(YLIM)))),lty=2,col="grey50" )
					points( CF.lis[,col],CF.lis[,row], col=COLS[1] )
					points( CF.lme[,col],CF.lme[,row], col=COLS[2] )
					abline(lm( CF.lis[,row]~CF.lis[,col] ),col=COLS.4[1],lty=2,lwd=2 )
					abline(lm( CF.lme[,row]~CF.lme[,col] ),col=COLS.4[2],lty=2,lwd=2 )
					points( mean(CF.lis[,col],na.rm=T),mean(CF.lis[,row],na.rm=T),pch=10,col=COLS.4[1],cex=2.5,lwd=3 )
					points( mean(CF.lme[,col],na.rm=T),mean(CF.lme[,row],na.rm=T),pch=10,col=COLS.4[2],cex=2.5,lwd=3 )
				}
				## Hist
				if ( row == col ) {
					XLIM <- range( CF.lis[,row],CF.lme[,row], na.rm=T )
					BRKS <- seq( floor(XLIM[1]),ceiling(XLIM[2]),round(diff(XLIM)/10,ceiling(-log10(diff(XLIM)/10))) )
					hist( CF.lme[,row], col=COLS[2],density=20,angle=-45 ,xlim=XLIM,breaks=BRKS,main=paste("Beta Estimates:",VARS[row]),xlab=VARS[row] )
					hist( CF.lis[,row], col=COLS[1],density=20,angle=45, add=T, breaks=BRKS )
				}
				## Corr
				if ( row < col ) {
					plot( 0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",main=paste(VARS[row],"vs",VARS[col]) )
					COR.lis <- cor( CF.lis[,row],CF.lis[,col], use="pairwise.complete.obs" )
					COR.lme <- cor( CF.lme[,row],CF.lme[,col], use="pairwise.complete.obs" )
					text( .5,.6, label=paste("LIS:",round(COR.lis,3)), cex=1.5, col=COLS[1] )
					text( .5,.4, label=paste("LME:",round(COR.lme,3)), cex=1.5, col=COLS[2] )
				}
			}else{
				# XLIM <- range( CF.lis[,row], na.rm=T )
				# YLIM <- range( CF.lme[,row], na.rm=T )
				LIM <- range( CF.lis[,row],CF.lme[,row], na.rm=T )
				plot( 0,0,type="n", xlim=LIM,ylim=LIM,xlab="LIS",ylab="LME",main=paste(VARS[row],"- LIS vs LME") )
				abline( v=seq(floor(LIM[1]),ceiling(LIM[2]),10^floor(log10(diff(LIM)))),lty=2,col="grey50" )
				abline( h=seq(floor(LIM[1]),ceiling(LIM[2]),10^floor(log10(diff(LIM)))),lty=2,col="grey50" )
				abline( 0,1,col="black" )
				points( CF.lis[,row],CF.lme[,row], col=COLS[3] )
				abline(lm( CF.lme[,row]~CF.lis[,row] ),col=COLS.4[3],lty=2,lwd=2 )
				points( mean(CF.lis[,row],na.rm=T),mean(CF.lme[,row],na.rm=T),pch=10,col=COLS.4[3],cex=2.5,lwd=3 )			
			}
			
		}
	}
} # Close DIST_PLOT function
DIST_PLOT( CF.lis, CF.lme )

## Notes
 # lmList will estimate separate coefs for each variable
 # Can only compare to LME models where all variables are random

## Best Direct Comparisons
 # DRUG: LIS$M1b vs LME$M2b
 # DRUG+WK: LIS$M2B vs XXXX
 # DRUG*WK: LIS$M3b vs XXXX

## Basic Model w/ Drug Effect
CF.lis.1 <- coef(LIS$M1b)
CF.lme.1 <- coef(LME$M2b)
png( paste(PathToPlot,"/VS_1.png",sep=""), height=100+400*ncol(CF.lis.3),width=400*(1+ncol(CF.lis.3)), pointsize=30 )
DIST_PLOT( CF.lis.1, CF.lme.1 )
dev.off()

## Model w/ Drug+WK Effects
# TEMP.2 <- lme( fixed = DAS ~ DRUG+WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
# TEMP.2 <- lme( fixed = DAS ~ DRUG+WK, random = ~ DRUG+WK | IID, data=TAB )
CF.lis.2 <- coef(LIS$M2b)
CF.lme.2 <- coef(TEMP.2)
png( paste(PathToPlot,"/VS_2.png",sep=""), height=100+400*ncol(CF.lis.3),width=400*(1+ncol(CF.lis.3)), pointsize=30 )
DIST_PLOT( CF.lis.2, CF.lme.2 )
dev.off()

## Model w/ Drug*WK Effects
# TEMP.3 <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
# TEMP.3 <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB )
CF.lis.3 <- coef(LIS$M3b)
CF.lme.3 <- coef(TEMP.3) # coef(LME$M6b)
png( paste(PathToPlot,"/VS_3.png",sep=""), height=100+400*ncol(CF.lis.3),width=400*(1+ncol(CF.lis.3)), pointsize=30 )
DIST_PLOT( CF.lis.3, CF.lme.3 )
dev.off()






##############################################################
## TESTING CorStruct #########################################
##############################################################


COR <- list()

for ( v in 1:10 ) {
	COR[[paste("M",v,sep=".")]] <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = v/10, form = ~1 | IID) )	
}



cor.1 <- corCAR1(value = .1, form = ~1 | IID)
cor.1 <- Initialize( cor.1, data=TAB )
cor.1.mat <- corMatrix( cor.1 )
heatmap.2( cor.1.mat$`F320833-27`, scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none" )

cor.2 <- corCAR1(value = .2, form = ~1 | IID)
cor.2 <- Initialize( cor.2, data=TAB )
cor.2.mat <- corMatrix( cor.2 )
heatmap.2( cor.2.mat$`F320833-27`, scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none" )

cor.9 <- corCAR1(value = .9, form = ~1 | IID)
cor.9 <- Initialize( cor.9, data=TAB )
cor.9.mat <- corMatrix( cor.9 )
heatmap.2( cor.9.mat$`F320833-27`, scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none" )

cor.9 <- corCAR1(value = .9, form = ~1 | DRUG/IID)
cor.9 <- Initialize( cor.9, data=TAB )
cor.9.mat <- corMatrix( cor.9 )
heatmap.2( cor.9.mat$`1/B012328-28`, scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none" )
quartz() ; heatmap.2( cor.9.mat$`0/B012328-28`, scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none" )


cor.9 <- corCAR1(value = .9, form = ~WK | DRUG/IID)
cor.9 <- Initialize( cor.9, data=TAB )
cor.9.mat <- corMatrix( cor.9 )
heatmap.2( cor.9.mat$`1/B012328-28`, scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none" )
quartz() ; heatmap.2( cor.9.mat$`0/B012328-28`, scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none" )





LME$M8b$call
TEST <- TAB[ which(TAB$IID==TAB$IID[40]),]
TEST.pred.0 <- predict( LME$M8b, newdata=TEST, level=0 )
TEST.pred.1 <- predict( LME$M8b, newdata=TEST, level=1 )
plot( TEST$WK, TEST$DAS, type="l",lty=2 )
points( TEST$WK, TEST.pred.1, type="l",lty=1 )
points( TEST$WK, TEST.pred.0, type="l",lty=3 )

TEST.2 <- TEST
TEST.2$WK <- sort( sample( 1:100, nrow(TEST.2) ) )
TEST.2$DRUG <- 1 ; SWITCH <- sample(20:60,1)
TEST.2$DRUG[ which(TEST.2$WK<SWITCH) ] <- 0
TEST.2$PLAC <- 0
TEST.2$PLAC[ which(TEST.2$DRUG!=1 & TEST.2$WK>0) ] <- 1
TEST.pred.2 <- predict( LME$M8b, newdata=TEST.2, level=1 )
points( TEST.2$WK, TEST.pred.2, type="o",lty=1,pch=20 )
























##############################################################
## SCRAPS!!!! ################################################
##############################################################