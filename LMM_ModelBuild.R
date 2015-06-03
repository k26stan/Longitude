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
DATE <- "20150530"

# ## TSCC Paths
# PathToData <- "/projects/janssen/clinical/"
# PathToSave <- "/projects/janssen/clinical/Plots/20140522/"

## Mac Paths
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150506_Resp_v_Time.txt"
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150512_Resp_v_Time.txt"
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_LME/",sep="" )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T )
FUL <- read.table( paste(PathToFT,"20141229_Full_Table.txt",sep=""),sep="\t",header=T)
TAB.pc <- merge( TAB.l, FUL[,c("ID_2",paste("PC",1:3,sep=""))], by.x="IID",by.y="ID_2")
# Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Load Candidate Genotype Files
COMP.l <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", sep="\t",header=T)
COMP <- COMP.l[ which(!duplicated(COMP.l$SNP)), ]
RAW <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/TEST_CND.raw", sep="",header=T)

##############################################################
## FILTER DATA ###############################################
##############################################################

## Take out Patients who left before getting DRUG
RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
RM.exit <- which( TAB.l$IID %in% RM.exit.id )
TAB <- TAB.l[-RM.exit,c(1:15,17)]

## Remove DAS values that are NA
RM.na <- which(is.na( TAB$DAS ))
if (length(RM.na) > 0 ) { TAB <- TAB[-RM.na, ] }

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
	LME$M2a <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M2b <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M2c <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	LME$M2d <- lme( fixed = DAS ~ DRUG, random = ~ DRUG-1 | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M2e <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	# LME$M2f <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M2a,LME$M2b,LME$M2c,LME$M2d,LME$M2e)
	COMP <- PLOT_MODS( LME )

	## Include DRUG*WK Interaction as Fixed Effect
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 3"))
	LME$M3a <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M3b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M3c <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	LME$M3d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG-1 | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M3e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	LME$M3f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M3a,LME$M3b,LME$M3c,LME$M3d,LME$M3e,LME$M3f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )

	## Include Auto-correlation for WK within Patient
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 4"))
	LME$M4a <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M4b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M4c <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	LME$M4d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG-1 | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M4e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | COUN/IID) )
	LME$M4f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M4a,LME$M4b,LME$M4c,LME$M4d,LME$M4e,LME$M4f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )

	## Revert to Order only in Correlation
	 # (Include WK in Random Effects)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 5"))
	LME$M5a <- lme( fixed = DAS ~ DRUG*WK, random = ~ WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M5b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M5c <- lme( fixed = DAS ~ DRUG*WK, random = ~ WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	LME$M5e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	# LME$M5f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M5a,LME$M5b,LME$M5c,LME$M5e,LME$M5f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )

	## Test WK in Correlation again
	 # (Keep WK in Random Effects)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 6"))
	LME$M6a <- lme( fixed = DAS ~ DRUG*WK, random = ~ WK | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M6b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M6c <- lme( fixed = DAS ~ DRUG*WK, random = ~ WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	LME$M6e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | COUN/IID) )
	# LME$M6f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$Ms6a,LME$M6b,LME$M6c,LME$M6e,LME$M6f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )

	## Revert to Order only in Correlation
	 # (Include Fixed PLAC Effect)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 7"))
	LME$M7a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M7b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	# LME$M7c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	LME$M7e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	# LME$M7f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M7a,LME$M7b,LME$M7c,LME$M7e,LME$M7f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )

	## Test WK in Correlation again
	 # (Keep Fixed PLAC Effect)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 8"))
	LME$M8a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ WK | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M8b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M8c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	LME$M8e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | COUN/IID) )
	# LME$M8f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M8a,LME$M8b,LME$M8c,LME$M8e,LME$M8f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )

	## Switch to ARMA Correlation Structure
	 # Q = How many adjacent positions (pos vector) are used in Moving Average Correlation
	 # P = How many adjacent positions (pos vector) are used in Autoregression Correlation
	 # e.g., Moving Average w/ 2 adjacent positions: Q=2, P=0
	 # e.g., Autoregression w/ 3 adjacent positions: Q=0, P=3
	 # e.g., Hybrid using 3 for Moving Average & 2 for Autoregression: Q=3, P=2

	## Revert to Order only in Correlation
	 # (Keep Fixed PLAC Effect)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 9"))
	Q <- 2 ; P <- 0 
	LME$M9a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID) )
	LME$M9b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID) )
	LME$M9c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ WK | IID/DRUG, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG) )
	LME$M9e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | COUN/IID) )
	LME$M9f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID/DRUG, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M9b,LME$M9e,LME$M9f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )

	## Test WK in Correlation again
	 # (Keep Fixed PLAC Effect)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 10"))
	LME$M10b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = c(seq(.8,.1,length.out=Q),seq(.8,.1,length.out=P)), q=Q,p=P, form = ~WK | IID) )
	LME$M10e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = c(seq(.8,.1,length.out=Q),seq(.8,.1,length.out=P)), q=Q,p=P, form = ~WK | COUN/IID) )
	LME$M10f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID/DRUG, data=TAB, correlation=corCAR1(value = c(seq(.8,.1,length.out=Q),seq(.8,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M10b,LME$M10e,LME$M10f)
	COMP <- PLOT_MODS( LME )
	HEAT <- HEAT_MODS( LME )







	## Include Fixed PLAC Effect
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 7"))
	LME$M7a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ 1 | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M7b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$M7c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(form = ~1 | IID/DRUG) )
	LME$M7e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	LME$M7f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
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

LME_MOD <- function( TAB ) {

	LME <- list()
	## When using ANOVA to compare models:
	 # Greater (less negative) logLik is better model
	 # Smaller (less positive) AIC/BIC is better model

	## Linear Mixed Model w/ only DRUG effect
	LME$LME.1 <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB )
	# summary(LME$LME.1)
	## Add Drug as Random Effect
	LME$LME.1b <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.1b)
	## Add Nested Random Effect w/ IID/DRUG
	LME$LME.1c <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID/DRUG, data=TAB )
	# summary(LME$LME.1c)
	# anova(LME$LME.1,LME$LME.1b,LME$LME.1c)
	 # Choose LME$LME.1b b/c has lowest AIC/BIC and less negative logLik
	## Add Nested Random Effect w/ IID/DRUG
	LME$LME.1d <- lme( fixed = DAS ~ DRUG, random = ~ DRUG-1 | IID, data=TAB )
	# summary(LME$LME.1d)
	# anova(LME$LME.1b,LME$LME.1d)
	 # Choose LME$LME.1b b/c has lowest AIC/BIC and less negative logLik

	## Add DRUG*SEX Interaction
	LME$LME.2 <- lme( fixed = DAS ~ DRUG*SEX, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.2)
	## Change to DRUG+SEX
	LME$LME.2a <- lme( fixed = DAS ~ DRUG+SEX, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.2a)

	## Change to DRUG*AGE
	LME$LME.3 <- lme( fixed = DAS ~ DRUG*AGE, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.3)
	## Change to DRUG+AGE
	LME$LME.3a <- lme( fixed = DAS ~ DRUG+AGE, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.3a)

	## Change to DRUG*DIS_DUR
	LME$LME.4 <- lme( fixed = DAS ~ DRUG*DIS_DUR, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.4)
	## Change to DRUG+DIS_DUR
	LME$LME.4a <- lme( fixed = DAS ~ DRUG+DIS_DUR, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.4a)

	## Remove DIS_DUR, add DRUG*HT
	LME$LME.5 <- lme( fixed = DAS ~ DRUG*HT, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.5)
	## Change to DRUG+HT
	LME$LME.5a <- lme( fixed = DAS ~ DRUG+HT, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.5a)

	## Change to DRUG*WT
	LME$LME.6 <- lme( fixed = DAS ~ DRUG*WT, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.6)
	## Change to DRUG+WT
	LME$LME.6a <- lme( fixed = DAS ~ DRUG+WT, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.6a)

	## Change to DRUG*BMI
	LME$LME.7 <- lme( fixed = DAS ~ DRUG*BMI, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.7)
	## Change to DRUG+BMI
	LME$LME.7a <- lme( fixed = DAS ~ DRUG+BMI, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.7a)

	## Change to DRUG*ACPA
	LME$LME.8 <- lme( fixed = DAS ~ DRUG*ACPA, random = ~ DRUG | IID, data=TAB, subset=ACPA!="" )
	# summary(LME$LME.8)
	## Change to DRUG+ACPA
	LME$LME.8a <- lme( fixed = DAS ~ DRUG+ACPA, random = ~ DRUG | IID, data=TAB, subset=ACPA!="" )
	# summary(LME$LME.8a)

	## Change to DRUG*RF
	LME$LME.9 <- lme( fixed = DAS ~ DRUG*RF, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.9)
	## Change to DRUG+RF
	LME$LME.9a <- lme( fixed = DAS ~ DRUG+RF, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.9a)

	## Change to DRUG*RF_ACPA
	LME$LME.10 <- lme( fixed = DAS ~ DRUG*RF_ACPA, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.10)
	## Change to DRUG+RF_ACPA
	LME$LME.10a <- lme( fixed = DAS ~ DRUG+RF_ACPA, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.10a)

	## Change to DRUG*WK
	LME$LME.11 <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.11) # *** WK & DRUG:WK terms are 99% correlated
	## Change to DRUG+WK
	LME$LME.11a <- lme( fixed = DAS ~ DRUG+WK, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.11a)
	## Back to DRUG*WK, add WK as Random Effect
	LME$LME.11b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.11b)
	## Back to DRUG+WK, keep WK as Random Effect
	LME$LME.11c <- lme( fixed = DAS ~ DRUG+WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.11c)
	# anova(LME$LME.11a,LME$LME.11c)
	## Back to DRUG*WK, add DRUG*WK as Random Effect
	# LME$LME.11d <- lme( fixed = DAS ~ DRUG+WK, random = ~ DRUG*WK | IID, data=TAB )
	# summary(LME$LME.11d)
	# anova(LME$LME.11a,LME$LME.11d)

	## Include PLAC (all 2-way interactions)
	# LME$LME.12 <- lme( fixed = DAS ~ (DRUG+WK+PLAC)^2, random = ~ DRUG | IID, data=TAB )
	# # summary(LME$LME.12)
	## Include PLAC (2-way interactions w/o PLAC:DRUG)
	# LME$LME.12a <- lme( fixed = DAS ~ DRUG*WK+PLAC+PLAC:WK, random = ~ DRUG | IID, data=TAB )
	# # summary(LME$LME.12a)
	## Include PLAC (w/o PLAC interactions)
	LME$LME.12b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.12b)
	## Add PLAC as Random Effect
	# LME$LME.12c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK+PLAC | IID, data=TAB )
	# # summary(LME$LME.12c) # !!! DID NOT CONVERGE !!!
	## Try w/ only DRUG as Random Effect
	LME$LME.12d <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | IID, data=TAB )
	# summary(LME$LME.12d)
	# anova(LME$LME.12b,LME$LME.12d) # Stick w/ WK in Random Effects (b)

	## Use 12b, add SEX interactions
	LME$LME.13 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*SEX, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.13)
	## Remove SEX:DRUG:WK
	LME$LME.13a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*SEX-DRUG:SEX:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.13a)
	## Remove SEX:WK
	LME$LME.13b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*SEX-DRUG:SEX:WK-SEX:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.13b)
	## Remove SEX:PLAC
	LME$LME.13c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*SEX, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.13c)
	## Remove SEX:DRUG
	LME$LME.13d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+SEX, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.13d)
	 # !!! No effect of SEX at all !!!

	## Use 12b, add RF_ACPA interactions
	LME$LME.14 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*RF_ACPA, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.14)
	## Remove RF_ACPA:DRUG:WK
	LME$LME.14a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*RF_ACPA-DRUG:RF_ACPA:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.14a)
	## Remove RF_ACPA:WK
	LME$LME.14b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*RF_ACPA-DRUG:RF_ACPA:WK-RF_ACPA:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.14b)
	## Remove RF_ACPA:PLAC
	LME$LME.14c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.14c)
	## Remove RF_ACPA:DRUG
	LME$LME.14d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+RF_ACPA, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.14d)
	## Include RF_ACPA as Grouping
	# LME$LME.14ca <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | RF_ACPA/IID, data=TAB )
	# # summary(LME$LME.14ca) ## Didn't Converge...first try...try again...
	## Include RF_ACPA as Random Effect
	# LME$LME.14da <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+RF_ACPA, random = ~ DRUG+WK | RF_ACPA/IID, data=TAB )
	# # summary(LME$LME.14da)
	## Decide: w/ or w/o RF_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.14c), predict(LME$LME.14d) ))
	# sd(resid(LME$LME.14c)) ; sd(resid(LME$LME.14d)) ; sd(resid(LME$LME.12b))

	## Use 12b, add DIS_DUR interactions
	LME$LME.15 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*DIS_DUR, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.15)
	## Remove DIS_DUR:DRUG:WK
	LME$LME.15a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*DIS_DUR-DRUG:DIS_DUR:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.15a)
	## Remove DIS_DUR:WK
	LME$LME.15b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*DIS_DUR-DRUG:DIS_DUR:WK-DIS_DUR:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.15b)
	## Remove DIS_DUR:PLAC
	LME$LME.15c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*DIS_DUR, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.15c)
	## Remove DIS_DUR:DRUG
	LME$LME.15d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DIS_DUR, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.15d)
	## Decide: w/ or w/o RF_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.15b), predict(LME$LME.15c), predict(LME$LME.15d) ))
	# sd(resid(LME$LME.15b)) ; sd(resid(LME$LME.15c)) ; sd(resid(LME$LME.15d)) ; sd(resid(LME$LME.12b))

	## Use 12b, add BMI interactions
	LME$LME.16 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*BMI, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.16)
	## Remove BMI:DRUG:WK
	LME$LME.16a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*BMI-DRUG:BMI:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.16a)
	## Remove BMI:WK
	LME$LME.16b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*BMI-DRUG:BMI:WK-BMI:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.16b)
	## Remove BMI:PLAC
	LME$LME.16c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*BMI, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.16c)
	## Remove BMI:DRUG
	LME$LME.16d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+BMI, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.16d)
	## Decide: w/ or w/o RF_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.16b), predict(LME$LME.16c), predict(LME$LME.16d) ))
	# sd(resid(LME$LME.16b)) ; sd(resid(LME$LME.16c)) ; sd(resid(LME$LME.16d)) ; sd(resid(LME$LME.12b))

	## Use 12b, add AGE interactions
	LME$LME.17 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*AGE, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.17)
	## Remove AGE:DRUG:WK
	LME$LME.17a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*AGE-DRUG:AGE:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.17a)
	## Remove AGE:WK
	LME$LME.17b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*AGE-DRUG:AGE:WK-AGE:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.17b)
	## Remove AGE:PLAC
	LME$LME.17c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*AGE, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.17c)
	## Remove AGE:DRUG
	LME$LME.17d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+AGE, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.17d)
	## Decide: w/ or w/o RF_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.17b), predict(LME$LME.17c), predict(LME$LME.17d) ))
	# sd(resid(LME$LME.17b)) ; sd(resid(LME$LME.17c)) ; sd(resid(LME$LME.17d)) ; sd(resid(LME$LME.12b))

	## Use 12b, add RF interactions
	LME$LME.18 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*RF, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.18)
	## Remove RF:DRUG:WK
	LME$LME.18a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*RF-DRUG:RF:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.18a)
	## Remove RF:WK
	LME$LME.18b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*RF-DRUG:RF:WK-RF:WK, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.18b)
	## Remove RF:PLAC
	LME$LME.18c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.18c)
	## Remove RF:DRUG
	LME$LME.18d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+RF, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.18d)
	## Decide: w/ or w/o RF_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.18b), predict(LME$LME.18c), predict(LME$LME.18d) ))
	# sd(resid(LME$LME.18b)) ; sd(resid(LME$LME.18c)) ; sd(resid(LME$LME.18d)) ; sd(resid(LME$LME.12b))

	## Use 12b, add ACPA interactions
	LME$LME.19 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*ACPA, random = ~ DRUG+WK | IID, data=TAB, subset=ACPA!="" )
	# summary(LME$LME.19)
	## Remove ACPA:DRUG:WK
	LME$LME.19a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*ACPA-DRUG:ACPA:WK, random = ~ DRUG+WK | IID, data=TAB, subset=ACPA!="" )
	# summary(LME$LME.19a)
	## Remove ACPA:WK
	LME$LME.19b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)*ACPA-DRUG:ACPA:WK-ACPA:WK, random = ~ DRUG+WK | IID, data=TAB, subset=ACPA!="" )
	# summary(LME$LME.19b)
	## Remove ACPA:PLAC
	LME$LME.19c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*ACPA, random = ~ DRUG+WK | IID, data=TAB, subset=ACPA!="" )
	# summary(LME$LME.19c)
	## Remove ACPA:DRUG
	LME$LME.19d <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+ACPA, random = ~ DRUG+WK | IID, data=TAB, subset=ACPA!="" )
	# summary(LME$LME.19d)
	## Decide: w/ or w/o ACPA_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.19b), predict(LME$LME.19c), predict(LME$LME.19d) ))
	# sd(resid(LME$LME.19b)) ; sd(resid(LME$LME.19c)) ; sd(resid(LME$LME.19d)) ; sd(resid(LME$LME.12b))

	###### SO FAR, 14c & 14d are only models with ALL significant terms ######
	############## Use 14c (includes RF_ACPA*DRUG) and go from here ##########

	## Use 14c, add DIS_DUR*DRUG interactions
	LME$LME.20 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA+DRUG*DIS_DUR, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.20)
	## Remove RF_ACPA:PLAC
	LME$LME.20a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA+DIS_DUR, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.20a)
	## Decide: w/ or w/o ACPA_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.12b), predict(LME$LME.14c), predict(LME$LME.14d), predict(LME$LME.20), predict(LME$LME.20a) ))
	# sd(resid(LME$LME.12b)); sd(resid(LME$LME.14c)) ; sd(resid(LME$LME.14d)) ; sd(resid(LME$LME.20)) ; sd(resid(LME$LME.20a))

	## Use 14c, add AGE*DRUG interactions
	LME$LME.21 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA+DRUG*AGE, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.21)
	## Remove RF_ACPA:PLAC
	LME$LME.21a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA+AGE, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.21a)
	## Decide: w/ or w/o ACPA_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.12b), predict(LME$LME.14c), predict(LME$LME.14d), predict(LME$LME.21), predict(LME$LME.21a) ))
	# sd(resid(LME$LME.12b)); sd(resid(LME$LME.14c)) ; sd(resid(LME$LME.14d)) ; sd(resid(LME$LME.21)) ; sd(resid(LME$LME.21a))

	## Use 14c, add BMI*DRUG interactions
	LME$LME.22 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA+DRUG*BMI, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.22)
	## Remove RF_ACPA:PLAC
	LME$LME.22a <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA+BMI, random = ~ DRUG+WK | IID, data=TAB )
	# summary(LME$LME.22a)
	## Decide: w/ or w/o ACPA_ACPA:DRUG interaction??
	# cor(data.frame( TAB$DAS, predict(LME$LME.12b), predict(LME$LME.14c), predict(LME$LME.14d), predict(LME$LME.22), predict(LME$LME.22a) ))
	# sd(resid(LME$LME.12b)); sd(resid(LME$LME.14c)) ; sd(resid(LME$LME.14d)) ; sd(resid(LME$LME.22)) ; sd(resid(LME$LME.22a))

	###### SO FAR, 14c & 14d are only models with ALL significant terms ######
	############## Use 14c (includes RF_ACPA*DRUG) and go from here ##########
	## PROBLEM: Trajectory during Placebo may not match Trajectory during Treatment

	## Use 14c, include DRUG*WK interactions in Random Effects
	# LME$LME.23 <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA, random = ~ DRUG*WK | IID, data=TAB )
	# # summary(LME$LME.23) # Did not Converge, Try Simpler Fixed Effects
	## Get Rid of RF_ACPA, just keep WK, DRUG, PLAC
	# LME$LME.23a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG*WK | IID, data=TAB )
	# # summary(LME$LME.23a) # Did not Converge!!!!
	## Try including PLAC*WK interaction as well, instead of DRUG*WK in Random Effects
	# LME$LME.23b <- lme( fixed = DAS ~ DRUG*WK+PLAC*WK, random = ~ DRUG+WK | IID, data=TAB )
	# # summary(LME$LME.23b) # Did not Converge!!!!
	## Try including PLAC*WK interaction as well, instead of DRUG*WK in Random Effects
	# LME$LME.23b <- lme( fixed = DAS ~ WK*factor(TRT), random = ~ 1 | IID, data=TAB )
	# # summary(LME$LME.23b) # Did not Converge!!!!

	## Use 12b, include COUN as Hierarchical Grouping Factor as well
	print("Running Mods: 23")
	LME$LME.23 <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB )
	 # Switch to COUN/IID (IID %in% COUN) in simple model
	LME$LME.23a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | COUN/IID, data=TAB )
	 # Switch to COUN/IID (IID %in% COUN) in simple model
	LME$LME.23b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG | COUN/IID, data=TAB )
	 # Switch to COUN/IID (IID %in% COUN) in simple model
	LME$LME.23c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | COUN/IID, data=TAB )
	anova( LME$LME.12b, LME$LME.23, LME$LME.12d, LME$LME.23b )

	## Use 23, include CorStruct
	LME$LME.24 <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	LME$LME.24a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$LME.24b <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(form = ~1 | IID) )
	LME$LME.24c <- lme( fixed = DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(form = ~1 | COUN/IID) )
	anova( LME$LME.24, LME$LME.24a )

	###### SO FAR, 14c & 14d are only models with ALL significant terms ######
	## Return Compiled Models
	return(LME)

}

LME <- LME_MOD( TAB )

head(coef(LME$LME.12b))
head(coef(LME$LME.14c))
head(coef(LME$LME.20))

##############################################################
## COMPILE INFO ABOUT MODELS #################################
##############################################################

## Specify Number and Names of Models
LME.N.models <- length(LME)
LME.Names.models <- names(LME)

## Determine which models have DRUG interactions
LME.drugint <- as.numeric(unlist( lapply( LME, function(x) any(grepl( "DRUG:|:DRUG", names(coef(x)) ))) ))

################################################
## Compile Terms and Significance for each model

## Get all Terms used in any model
LME.terms <- Reduce( union, lapply( LME, function(x) colnames(coef(x)) ) )
# LME.terms <- sort( LME.terms )
LME.N.terms <- length(LME.terms)

## Determine which Terms are in each model
LME.ARR.terms <- array( 0, c(LME.N.terms,LME.N.models) )
LME.ARR.terms.p <- array( 1, c(LME.N.terms,LME.N.models) )
colnames(LME.ARR.terms) <- colnames(LME.ARR.terms.p) <- LME.Names.models
rownames(LME.ARR.terms) <- rownames(LME.ARR.terms.p) <- LME.terms
for ( m in 1:length(LME) ) {
	mod <- names(LME)[m]
	which.in <- which( LME.terms %in% colnames(coef(LME[[mod]])) )
	which.terms <- LME.terms[which.in]
	LME.ARR.terms[which.in,mod] <- 1
	which.terms.2 <- which.terms[ which( which.terms %in% rownames(summary(LME[[mod]])$tTable) ) ]
	which.in.2 <- which( LME.terms %in% which.terms.2 )
	LME.ARR.terms.p[which.in.2,mod] <- summary(LME[[mod]])$tTable[which.terms.2,"p-value"]
}
 # Plot P-Values for each Term
COLS.list <- c("firebrick1","chocolate1","gold1","springgreen2","steelblue2","slateblue3","black")
COLS <- c( "firebrick1",colorRampPalette(COLS.list[3:7])(19) ) # colorRampPalette(COLS.list)(20)
ROW.COLS <- rep("deepskyblue2",LME.N.terms) ; ROW.COLS[grep("DRUG",LME.terms)] <- "chartreuse2"
png( paste(PathToPlot,"1-LME_ModBuild.png",sep="/"), height=1200,width=1600,pointsize=30 )
heatmap.2( LME.ARR.terms.p, main="Inclusion/Significance of Terms in Models",xlab="Model",ylab="Term",trace="none",scale="none",Colv=F,Rowv=F,dendrogram="none",col=COLS,margin=c(6,11),lhei=c(2,9),lwid=c(1,6),RowSideColors=ROW.COLS )
dev.off()

################################################
## Compile Fit Stats for each Model
 # AIC/BIC/logLik
LME.BIC <- unlist(lapply( LME, BIC ))
LME.AIC <- unlist(lapply( LME, AIC ))
LME.LL <- unlist(lapply( LME, logLik ))
 # Residual Vairance
LME.RES <- unlist(lapply( LME, function(x) sd(resid(x)) ))

## Plot Stats
COLS <- c("cadetblue2","tomato2","mediumpurple3","seagreen2")
png( paste(PathToPlot,"1-LME_ModBuild_Assess.png",sep="/"), height=1000,width=2000,pointsize=30 )
par( mfrow=c(1,3) )
 # LL
plot( 1:LME.N.models, LME.LL, col=COLS[1], main="Model Assessment: LogLik",ylab="LogLik",xlab="Model", pch=20,type="o" )
abline( v=which.max(LME.LL), col="firebrick1", lty=2 )
text( which.max(LME.LL)-1,quantile(range(LME.LL),.8), label=names(LME.LL)[which.max(LME.LL)],srt=90,col="firebrick1")
 # IC
plot( 1:LME.N.models, LME.AIC, col=COLS[2], main="Model Assessment: Info Criteria",ylab="AIC & BIC",xlab="Model", pch=20,type="o" )
points( 1:LME.N.models, LME.BIC, col=COLS[3], pch=20,type="o" )
abline( v=which.min(LME.AIC), col="firebrick1", lty=2 )
abline( v=which.min(LME.BIC), col="firebrick1", lty=2 )
legend( "topright", col=COLS[2:3],legend=c("AIC","BIC"), pch=20,lwd=1 )
text( which.min(LME.BIC)-1,quantile(range(LME.BIC),.2), label=names(LME.BIC)[which.min(LME.BIC)],srt=90,col="firebrick1")
 # RES
plot( 1:LME.N.models, LME.RES, col=COLS[4], main="Model Assessment: sd(Residuals)",ylab="sd(Residuals)",xlab="Model", pch=20,type="o" )
abline( v=which.min(LME.RES), col="firebrick1", lty=2 )
text( which.min(LME.RES)-1,quantile(range(LME.RES),.2), label=names(LME.RES)[which.min(LME.RES)],srt=90,col="firebrick1")
dev.off()
## More Refined Models
WHICH_MODS <- which(apply( LME.ARR.terms, 2, function(x) any(grep( "WK", rownames( LME.ARR.terms )[which(x==1)] )) ))
COLS <- c("cadetblue2","tomato2","mediumpurple3","seagreen2")
png( paste(PathToPlot,"1-LME_ModBuild_AssessWK.png",sep="/"), height=1000,width=2000,pointsize=30 )
par( mfrow=c(1,3) )
 # LL
plot( WHICH_MODS, LME.LL[WHICH_MODS], col=COLS[1], main="Model Assessment: LogLik",ylab="LogLik",xlab="Model", pch=20,type="o" )
abline( v=which.max(LME.LL), col="firebrick1", lty=2 )
text( which.max(LME.LL)-1,quantile(range(LME.LL),.8), label=names(LME.LL)[which.max(LME.LL)],srt=90,col="firebrick1")
 # IC
plot( WHICH_MODS, LME.AIC[WHICH_MODS], col=COLS[2], main="Model Assessment: Info Criteria",ylab="AIC & BIC",xlab="Model", pch=20,type="o" )
points( WHICH_MODS, LME.BIC[WHICH_MODS], col=COLS[3], pch=20,type="o" )
abline( v=which.min(LME.AIC), col="firebrick1", lty=2 )
abline( v=which.min(LME.BIC), col="firebrick1", lty=2 )
legend( "topright", col=COLS[2:3],legend=c("AIC","BIC"), pch=20,lwd=1 )
text( which.min(LME.BIC)-1,quantile(range(LME.BIC),.2), label=names(LME.BIC)[which.min(LME.BIC)],srt=90,col="firebrick1")
 # RES
plot( WHICH_MODS, LME.RES[WHICH_MODS], col=COLS[4], main="Model Assessment: sd(Residuals)",ylab="sd(Residuals)",xlab="Model", pch=20,type="o" )
abline( v=which.min(LME.RES), col="firebrick1", lty=2 )
text( which.min(LME.RES)-1,quantile(range(LME.RES),.2), label=names(LME.RES)[which.min(LME.RES)],srt=90,col="firebrick1")
dev.off()

##############################################################
## PLOT SOME PREDICTED VALUES ################################
##############################################################

## A few random Patients
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
N.patients <- 4
First.patient <- sample(1:436,1)
COLS <- colorRampPalette(sample(COLS.list))(N.patients)
COLS <- c("#6959CD", "#A75268", "#9ED527", "#EE7621")
png( paste(PathToPlot,"1-LME_Mod_Pred.png",sep="/"), height=1000,width=1400,pointsize=30 )
plot( 0,0,type="n", xlim=c(0,100),ylim=c(1,9),xlab="Weeks",ylab="DAS", main="Predicted vs Actual DAS Over Time" )
abline( h=seq(0,10,1), lty=2,col="grey50",lwd=1 )
for ( n in 1:N.patients ) {
	samp <- unique(TAB$IID)[First.patient+n-1]
	ind <- which(TAB$IID==samp)
	PCHS <- c(1,16)[factor(TAB$DRUG[ind])]
	SWITCH <- TAB$WK[ind][which(TAB$DRUG[ind]==1)][1]
	OBS <- TAB[ ind, "DAS" ]
	EXP.LME.1 <- predict(LME$LME.1b)[ind]
	EXP.LME.12 <- predict(LME$LME.12b)[ind]
	EXP.LME.14 <- predict(LME$LME.14c)[ind]
	points( TAB[ind,"WK"], OBS, col=COLS[n], type="o",pch=PCHS,lty=3,lwd=2 )
	points( TAB[ind,"WK"], EXP.LME.1, col=COLS[n], type="l",pch=PCHS,lty=2,lwd=3 )
	points( TAB[ind,"WK"], EXP.LME.14, col=COLS[n], type="l",pch=PCHS,lty=1,lwd=3 )
	abline( v=jitter(SWITCH,amount=.5), col=COLS[n],lty=2,lwd=1 )
	# points( TAB[ind,"WK"], EXP.LME.12, col=COLS[n], type="o",pch=PCHS,lty=3,lwd=2 )
}
dev.off()
coef(LME$LME.12b)[ First.patient:(First.patient+N.patients-1) , ]

## Histogram of Coefficients (random Effects)
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
COLS <- sample(COLS.list)
COLS <- COLS.list[c(2,5,3,6)]
png( paste(PathToPlot,"1-LME_Coef_Dist.png",sep="/"), height=1000,width=2000,pointsize=30 )
par(mfrow=c(1,4))
hist( coef(LME$LME.14)[,"(Intercept)"], col=COLS[1], main="Coefficient Distribution", xlab="Intercept" )
hist( coef(LME$LME.14)[,"DRUG"], col=COLS[2], main="Coefficient Distribution", xlab="Drug" )
hist( coef(LME$LME.14)[,"WK"], col=COLS[3], main="Coefficient Distribution", xlab="Week" )
hist( coef(LME$LME.14)[,"WK"]+coef(LME$LME.14)[,"DRUG:WK"], col=COLS[4], main="Coefficient Distribution", xlab="Week (on Drug)" )
dev.off()

## Pairs/Distribution/Correlation b/n Beta Estimates
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
COLS <- sample(COLS.list)
COLS <- COLS.list[c(2,5,3,6)]
COLS <- COLS.list[c(1,3,4,5)]
B.1 <- coef(LME$LME.14)[,"(Intercept)"]
B.2 <- coef(LME$LME.14)[,"DRUG"]
B.3 <- coef(LME$LME.14)[,"WK"]
B.4 <- coef(LME$LME.14)[,"WK"] + coef(LME$LME.14)[,"DRUG:WK"]
png( paste(PathToPlot,"1-LME_Coef_Dist_Cor.png",sep="/"), height=2000,width=2000,pointsize=30 )
par(mfrow=c(4,4))
 # 1
hist( B.1, col=COLS[1] )
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1)) ; text(.5,.5, round(cor(B.1,B.2),3),col=colorRampPalette(COLS[c(1,2)])(3)[2], cex=2.5 )
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1)) ; text(.5,.5, round(cor(B.1,B.3),3),col=colorRampPalette(COLS[c(1,3)])(3)[2], cex=2.5 )
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1)) ; text(.5,.5, round(cor(B.1,B.4),3),col=colorRampPalette(COLS[c(1,4)])(3)[2], cex=2.5 )
 # 2
plot( B.1, B.2, pch="+", col=colorRampPalette(COLS[c(1,2)])(3)[2] )
hist( B.2, col=COLS[2] )
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1)) ; text(.5,.5, round(cor(B.2,B.3),3),col=colorRampPalette(COLS[c(2,3)])(3)[2], cex=2.5 )
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1)) ; text(.5,.5, round(cor(B.2,B.4),3),col=colorRampPalette(COLS[c(2,4)])(3)[2], cex=2.5 )
 # 3
plot( B.1, B.3, pch="+", col=colorRampPalette(COLS[c(1,3)])(3)[2] )
plot( B.2, B.3, pch="+", col=colorRampPalette(COLS[c(2,3)])(3)[2] )
hist( B.3, col=COLS[3] )
plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1)) ; text(.5,.5, round(cor(B.3,B.4),3),col=colorRampPalette(COLS[c(3,4)])(3)[2], cex=2.5 )
 # 4
plot( B.1, B.4, pch="+", col=colorRampPalette(COLS[c(1,4)])(3)[2] )
plot( B.2, B.4, pch="+", col=colorRampPalette(COLS[c(2,4)])(3)[2] )
plot( B.3, B.4, pch="+", col=colorRampPalette(COLS[c(3,4)])(3)[2] )
hist( B.4, col=COLS[4] )
dev.off()

## All Data Points
# COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
# COLS <- COLS.list
# LIM <- c(0,10)
# png( paste(PathToPlot,"/Pl_Scatter_ObsExp.png",sep=""), height=2400,width=1800, pointsize=26 )
# par(mfrow=c(4,3))
# plot( data.frame( TAB$DAS, predict(LME.1) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LME.1b) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LM.1) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LME.11) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LME.11b) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LM.11) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( 0,0,type="n", col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LME.12b) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LM.12) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( 0,0,type="n", col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LME.14c) ), col=COLS[c(2,5)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# plot( data.frame( TAB$DAS, predict(LM.14) ), col=COLS[c(1,4)][factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM )
# abline( 0,1, lty=2 )
# dev.off()

# ## All Data Points
# pairs( data.frame( TAB$DAS, predict(LME.1b),predict(LME.11b),predict(LME.12b) ), col=COLS[factor(TAB$PLAC)], pch="." )

# png( "Dropbox/Schork/JNJ11/Slides/20150416_PA_Visit/MEM_Ex_2_hist2.png",height=1400,width=1400,pointsize=30 )
# par(mfrow=c(2,2))
# hist( coef(LME$LME.1b)[,1], main="Distribution of Intercept Terms",xlab="Beta(0)",ylab="# Patients",col="dodgerblue1" )
# plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(-1,1),ylim=c(-1,1) )
# text( 0,0, round( cor( coef(LME$LME.1b)[,1],coef(LME$LME.1b)[,2] ) ,3), col="purple2",cex=3 )
# plot( coef(LME$LME.1b)[,1],coef(LME$LME.1b)[,2], main="Correlation b/n Random Effects Estimates",xlab="Beta(0)",ylab="Beta(dr)",col="purple2" )
# hist( coef(LME$LME.1b)[,2], main="Distribution of Drug Terms",xlab="Beta(dr)",ylab="# Patients",col="firebrick1" )
# dev.off()

##############################################################
## INCORPORATE VARIANTS ######################################
##############################################################

############ I FORGOT WHAT I WAS DOING!!!!!!!!!! ############
## CHANGING COLNAMES & STUFF SO I KNOW WHERE VARIANTS ARE ON THE PLOT #####

## Pull out Candidate SNPs
THRSH <- 5e-6 # 5e-5
CAND.which <- which( COMP$P_Assoc < THRSH )
# CAND.which <- which( COMP$P_Assoc >= THRSH )
CAND.p <- COMP$P_Assoc[CAND.which]
CAND.id <- COMP$SNP[CAND.which] # paste( COMP$SNP[CAND.which], COMP$Allele[CAND.which], sep="_" )
CAND.loc <- paste( COMP$CHR[CAND.which], COMP$BP[CAND.which], sep="_" )
CAND.tag <- data.frame( ID=CAND.id, LOC=CAND.loc, P=CAND.p )

RAW.colnames <- unlist(lapply(strsplit( colnames(RAW), "_" ), function(x) paste(x[-length(x)],collapse="_") ))
RAW.colnames <- gsub(".",":",RAW.colnames,fixed=T )
RAW.colnames.which.X <- which( substr(RAW.colnames,1,2) %in% paste("X",1:22,sep="") )
RAW.colnames[RAW.colnames.which.X] <- gsub("X","",RAW.colnames[RAW.colnames.which.X],fixed=T )
# RAW.colnames <- RAW.colnames[ which(RAW.colnames!="") ]

CAND.intersect <- intersect( RAW.colnames, CAND.tag$ID )
CAND.rawcols <- which( RAW.colnames %in% CAND.intersect )
CAND.tagrows <- which( CAND.tag$ID %in% CAND.intersect )
RAW.cand <- RAW[,c(1,CAND.rawcols)]
CAND.tag.2 <- CAND.tag[ CAND.tagrows, ]
CAND.colnames <- colnames(RAW.cand)[2:ncol(RAW.cand)]

## Merge Candidate SNPs w/ Clinical Variables
MG <- merge( x=TAB.2, y=RAW.cand, by="FID" )

## Test Model that includes each Variant Individually
LME.VAR <- list()
GWAS.p <- numeric( length(CAND.colnames) )
GWAS.id <- rep("", length(CAND.colnames) )
start_time <- proc.time()
for ( c in 1:length(CAND.colnames) ) {
# for ( c in 1:5 ) {
	cand <- CAND.colnames[c]
	cand.loc <- as.character(CAND.tag.2$LOC[c])
	TEMP_DAT <- MG[,c("IID","DAS","DRUG",cand,"WK","PLAC","PC1","PC2")]
	colnames(TEMP_DAT)[4] <- "cand"
	# LME.VAR[[cand.loc]] <- lme( DAS ~ DRUG*cand, random= ~ DRUG | IID, data=TEMP_DAT )
	# LME.VAR[[cand.loc]] <- lme( DAS ~ DRUG*cand+PC1+PC2, random= ~ DRUG | IID, data=TEMP_DAT )
	LME.VAR[[cand.loc]] <- lme( fixed = DAS ~ (DRUG*WK)*cand+PLAC+PC1+PC2, random = ~ DRUG+WK | IID, data=TEMP_DAT )
	GWAS.p[c] <- CAND.tag.2$P[c]
	GWAS.id[c] <- as.character(CAND.tag.2$ID[c])
	names(GWAS.p)[c] <- cand.loc
	run_time <- round( proc.time()-start_time, 2)[3]
	if ( c%%2==0 ) { print(paste( "Done with",c,"of",length(CAND.colnames),"-",run_time )) }
}

## Pull out p-value for Variants
LME.VAR.p.cand <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable["cand","p-value"] ))
LME.VAR.p.cand.dr <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable["DRUG:cand","p-value"] ))
LME.VAR.p.cand.wk <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable["WK:cand","p-value"] ))
LME.VAR.p.cand.dr.wk <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable["DRUG:WK:cand","p-value"] ))
P.comp <- data.frame( ID=GWAS.id, GWAS=GWAS.p, DIS=LME.VAR.p.cand, DR=LME.VAR.p.cand.dr, WK=LME.VAR.p.cand.wk, DR_WK=LME.VAR.p.cand.dr.wk )

## Plot P-Values for Variants vs Drug Response
COLS <- c("firebrick2","dodgerblue2","chartreuse2","purple2","gold1","chocolate2")
png( paste(PathToPlot,"V-LME_P-Vals.png",sep="/"), height=1000+nrow(P.comp),width=1600+5*nrow(P.comp),pointsize=30 )
plot( -log10(P.comp$GWAS), main="Variant Effect in Mixed Model", col=COLS[1],pch="+", ylim=c(0,9),xaxt="n",xlab="Variant",ylab="-log10(p)",cex=1.5 )
points( -log10(P.comp$DIS), col=COLS[2],pch="+",cex=1.5 )
points( -log10(P.comp$DR), col=COLS[3],pch="+",cex=1.5 )
points( -log10(P.comp$WK), col=COLS[4],pch="+",cex=1.5 )
points( -log10(P.comp$DR_WK), col=COLS[5],pch="+",cex=1.5 )
axis( 1, at=1:nrow(P.comp), label=rownames(P.comp),las=2,cex.axis=.7 )
abline( h=seq(0,15,1),lty=2,col="grey50")
abline( h=-log10(5e-8),lty=2,col="firebrick2")
abline( h=-log10(.05),lty=2,col="firebrick2")
legend( "bottomleft",legend=colnames(P.comp)[-1],col=COLS,pch="+",pt.cex=1.5,ncol=2)
dev.off()
 # Only variants w/ P<1e-6
P.comp.2 <- P.comp[ which(apply( P.comp[c("GWAS","DIS","DR")], 1, function(x) any(x<1e-6) )) , ]
png( paste(PathToPlot,"V-LME_P-Vals_CAND.png",sep="/"), height=1000+nrow(P.comp.2),width=1600+5*nrow(P.comp.2),pointsize=30 )
plot( -log10(P.comp.2$GWAS), main="Variant Effect in Mixed Model", col=COLS[1],pch="+", ylim=c(0,9),xaxt="n",xlab="Variant",ylab="-log10(p)",cex=1.5 )
points( -log10(P.comp.2$DIS), col=COLS[2],pch="+",cex=1.5 )
points( -log10(P.comp.2$DR), col=COLS[3],pch="+",cex=1.5 )
points( -log10(P.comp.2$WK), col=COLS[4],pch="+",cex=1.5 )
points( -log10(P.comp.2$DR_WK), col=COLS[5],pch="+",cex=1.5 )
axis( 1, at=1:nrow(P.comp.2), label=rownames(P.comp.2),las=2,cex.axis=.7 )
abline( h=seq(0,15,1),lty=2,col="grey50")
abline( h=-log10(5e-8),lty=2,col="firebrick2")
abline( h=-log10(.05),lty=2,col="firebrick2")
legend( "bottomleft",legend=colnames(P.comp.2)[-1],col=COLS,pch="+",pt.cex=1.5,ncol=2)
dev.off()
P.comp.3 <- P.comp[ which(apply( P.comp[c("GWAS","DIS","DR")], 1, function(x) any(x<1e-7) )) , ]

## Plot GWAS results vs LMM Results
quartz()
LIM <- c(0, -log10(min(P.comp[,-1],na.rm=T)) )
plot( -log10(P.comp$GWAS),-log10(P.comp$DR),xlim=LIM,ylim=LIM,xlab="GWAS",ylab="LMM",main="GWAS vs LMM for Drug Response",pch="+",cex=1.5,col=COLS[3]) # ,col=colorRampPalette(COLS[c(1,3)])(3)[2] )
abline( v=seq(0,10,1),lty=2,col="grey50") ; abline( h=seq(0,10,1),lty=2,col="grey50") 
abline(0,1,lty=1,col="black")
abline( h=-log10(5e-8),lty=2,col="firebrick2") ; abline( v=-log10(5e-8),lty=2,col="firebrick2")
cor( -log10(P.comp$GWAS),-log10(P.comp$DR) )

##############################################################
## END OF DOC ################################################
##############################################################



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
## LINEAR MODELS #############################################
##############################################################
LM <- list()

## Basic Linear Model w/ only DRUG effect
LM$LM.1 <- lm( DAS ~ DRUG, data=TAB )
summary(LM$LM.1)

## Add DRUG*SEX Interaction
LM$LM.2 <- lm( DAS ~ DRUG*SEX, data=TAB )
summary(LM$LM.2)
## Change to DRUG+SEX
LM$LM.2a <- lm( DAS ~ DRUG+SEX, data=TAB )
summary(LM$LM.2a)

## Remove SEX, add DRUG*AGE
LM$LM.3 <- lm( DAS ~ DRUG*AGE, data=TAB )
summary(LM$LM.3)
## Change to DRUG+AGE
LM$LM.3a <- lm( DAS ~ DRUG+AGE, data=TAB )
summary(LM$LM.3a)

## Remove SEX, add DRUG*DIS_DUR
LM$LM.4 <- lm( DAS ~ DRUG*DIS_DUR, data=TAB )
summary(LM$LM.4)
## Change to DRUG+DIS_DUR
LM$LM.4a <- lm( DAS ~ DRUG+DIS_DUR, data=TAB )
summary(LM$LM.4a)

## Remove DIS_DUR, add DRUG*HT
LM$LM.5 <- lm( DAS ~ DRUG*HT, data=TAB )
summary(LM$LM.5)
## Change to DRUG+HT
LM$LM.5a <- lm( DAS ~ DRUG+HT, data=TAB )
summary(LM$LM.5a)

## Change to DRUG*WT
LM$LM.6 <- lm( DAS ~ DRUG*WT, data=TAB )
summary(LM$LM.6)
## Change to DRUG+WT
LM$LM.6a <- lm( DAS ~ DRUG+WT, data=TAB )
summary(LM$LM.6a)

## Change to DRUG*BMI
LM$LM.7 <- lm( DAS ~ DRUG*BMI, data=TAB )
summary(LM$LM.7)
## Change to DRUG+BMI
LM$LM.7a <- lm( DAS ~ DRUG+BMI, data=TAB )
summary(LM$LM.7a)

## Change to DRUG*ACPA
LM$LM.8 <- lm( DAS ~ DRUG*ACPA, data=TAB, subset=ACPA!="" )
summary(LM$LM.8)
## Change to DRUG+ACPA
LM$LM.8a <- lm( DAS ~ DRUG+ACPA, data=TAB, subset=ACPA!="" )
summary(LM$LM.8a)

## Change to DRUG*RF
LM$LM.9 <- lm( DAS ~ DRUG*RF, data=TAB )
summary(LM$LM.9)
## Change to DRUG+RF
LM$LM.9a <- lm( DAS ~ DRUG+RF, data=TAB )
summary(LM$LM.9a)

## Change to DRUG*RF_ACPA
LM$LM.10 <- lm( DAS ~ DRUG*RF_ACPA, data=TAB )
summary(LM$LM.10)
## Change to DRUG+RF_ACPA
LM$LM.10a <- lm( DAS ~ DRUG+RF_ACPA, data=TAB )
summary(LM$LM.10a)

## Change to DRUG*WK
LM$LM.11 <- lm( DAS ~ DRUG*WK, data=TAB )
summary(LM$LM.11)
## Change to DRUG+WK
LM$LM.11a <- lm( DAS ~ DRUG+WK, data=TAB )
summary(LM$LM.11a)

## Include PLAC (all 2-way interactions)
LM$LM.12 <- lm( DAS ~ (DRUG+WK+PLAC)^2, data=TAB )
summary(LM$LM.12)
## Include PLAC (2-way interactions w/o PLAC:DRUG)
LM$LM.12a <- lm( DAS ~ DRUG*WK+PLAC+PLAC:WK, data=TAB )
summary(LM$LM.12a)
## Include PLAC (w/o PLAC interactions)
LM$LM.12b <- lm( DAS ~ DRUG*WK+PLAC, data=TAB )
summary(LM$LM.12b)

## Use 12b, add SEX interactions
LM$LM.13 <- lm( DAS ~ (DRUG*WK+PLAC)*SEX, data=TAB )
summary(LM$LM.13)
## Remove SEX:DRUG:WK
LM$LM.13a <- lm( DAS ~ (DRUG*WK+PLAC)*SEX-DRUG:SEX:WK, data=TAB )
summary(LM$LM.13a)
## Remove SEX:WK
LM$LM.13b <- lm( DAS ~ (DRUG*WK+PLAC)*SEX-DRUG:SEX:WK-SEX:WK, data=TAB )
summary(LM$LM.13b)
## Remove SEX:PLAC
LM$LM.13c <- lm( DAS ~ (DRUG*WK+PLAC)+DRUG*SEX, data=TAB )
summary(LM$LM.13c)
## Remove SEX:DRUG
LM$LM.13d <- lm( DAS ~ (DRUG*WK+PLAC)+SEX, data=TAB )
summary(LM$LM.13d)

## Use 12b, add RF_ACPA interactions
LM$LM.14 <- lm( DAS ~ (DRUG*WK+PLAC)*RF_ACPA, data=TAB )
summary(LM$LM.14)
## Remove RF_ACPA:DRUG:WK
LM$LM.14a <- lm( DAS ~ (DRUG*WK+PLAC)*RF_ACPA-DRUG:RF_ACPA:WK, data=TAB )
summary(LM$LM.14a)
## Remove RF_ACPA:WK
LM$LM.14b <- lm( DAS ~ (DRUG*WK+PLAC)*RF_ACPA-DRUG:RF_ACPA:WK-RF_ACPA:WK, data=TAB )
summary(LM$LM.14b)
## Remove RF_ACPA:PLAC
LM$LM.14c <- lm( DAS ~ (DRUG*WK+PLAC)+DRUG*RF_ACPA, data=TAB )
summary(LM$LM.14c)
## Remove RF_ACPA:DRUG
LM$LM.14d <- lm( DAS ~ (DRUG*WK+PLAC)+RF_ACPA, data=TAB )
summary(LM$LM.14d)

## Use 12b, add DIS_DUR interactions
LM$LM.15 <- lm( DAS ~ (DRUG*WK+PLAC)*DIS_DUR, data=TAB )
summary(LM$LM.15)
## Remove DIS_DUR:DRUG:WK
LM$LM.15a <- lm( DAS ~ (DRUG*WK+PLAC)*DIS_DUR-DRUG:DIS_DUR:WK, data=TAB )
summary(LM$LM.15a)
## Remove DIS_DUR:WK
LM$LM.15b <- lm( DAS ~ (DRUG*WK+PLAC)*DIS_DUR-DRUG:DIS_DUR:WK-DIS_DUR:WK, data=TAB )
summary(LM$LM.15b)
## Remove DIS_DUR:PLAC
LM$LM.15c <- lm( DAS ~ (DRUG*WK+PLAC)+DRUG*DIS_DUR, data=TAB )
summary(LM$LM.15c)
## Remove DIS_DUR:DRUG
LM$LM.15d <- lm( DAS ~ (DRUG*WK+PLAC)+DIS_DUR, data=TAB )
summary(LM$LM.15d)

## Use 12b, add BMI interactions
LM$LM.16 <- lm( DAS ~ (DRUG*WK+PLAC)*BMI, data=TAB )
summary(LM$LM.16)
## Remove BMI:DRUG:WK
LM$LM.16a <- lm( DAS ~ (DRUG*WK+PLAC)*BMI-DRUG:BMI:WK, data=TAB )
summary(LM$LM.16a)
## Remove BMI:WK
LM$LM.16b <- lm( DAS ~ (DRUG*WK+PLAC)*BMI-DRUG:BMI:WK-BMI:WK, data=TAB )
summary(LM$LM.16b)
## Remove BMI:PLAC
LM$LM.16c <- lm( DAS ~ (DRUG*WK+PLAC)+DRUG*BMI, data=TAB )
summary(LM$LM.16c)
## Remove BMI:DRUG
LM$LM.16d <- lm( DAS ~ (DRUG*WK+PLAC)+BMI, data=TAB )
summary(LM$LM.16d)


################################################
## Compile Info about all Models

## Specify Number and Names of Models
LM.N.models <- length(LM)
LM.Names.models <- names(LM)

## Get all Terms used in any model
LM.terms <- Reduce( union, lapply( LM, function(x) names(coef(x)) ) )
# LM.terms <- sort( LM.terms )
LM.N.terms <- length(LM.terms)

## Determine which Terms are in each model
LM.ARR.terms <- array( 0, c(LM.N.terms,LM.N.models) )
LM.ARR.terms.p <- array( 1, c(LM.N.terms,LM.N.models) )
colnames(LM.ARR.terms) <- colnames(LM.ARR.terms.p) <- LM.Names.models
rownames(LM.ARR.terms) <- rownames(LM.ARR.terms.p) <- LM.terms
for ( m in 1:length(LM) ) {
	mod <- names(LM)[m]
	which.in <- which( LM.terms %in% names(coef(LM[[mod]])) )
	which.terms <- LM.terms[which.in]
	LM.ARR.terms[which.in,mod] <- 1
	which.terms.2 <- which.terms[ which( which.terms %in% rownames(summary(LM[[mod]])$coefficients) ) ]
	which.in.2 <- which( LM.terms %in% which.terms.2 )
	LM.ARR.terms.p[which.in.2,mod] <- summary(LM[[mod]])$coefficients[which.terms.2,"Pr(>|t|)"]
}
 # Presence of each Term
# COLS <- c("black","green")
# heatmap.2( LM.ARR.terms, trace="none",scale="none",Colv=F,Rowv=F,dendrogram="none",col=COLS,margin=c(5,10) )
 # P-Values for each Term
COLS.list <- c("firebrick1","chocolate1","gold1","springgreen2","steelblue2","slateblue3","black")
COLS <- c( "firebrick1",colorRampPalette(COLS.list[3:7])(19) ) # colorRampPalette(COLS.list)(20)
ROW.COLS <- rep("deepskyblue2",LM.N.terms) ; ROW.COLS[grep("DRUG",LM.terms)] <- "chartreuse2"
png( paste(PathToPlot,"1-LM_ModBuild.png",sep="/"), height=1200,width=1600,pointsize=30 )
# heatmap.2( LM.ARR.terms.p, trace="none",scale="none",Colv=F,Rowv=F,dendrogram="none",col=COLS,margin=c(5,10),RowSideColors=ROW.COLS )
heatmap.2( LM.ARR.terms.p, main="Inclusion/Significance of Terms in Models",xlab="Model",ylab="Term",trace="none",scale="none",Colv=F,Rowv=F,dendrogram="none",col=COLS,margin=c(6,11),lhei=c(2,9),lwid=c(1,6),RowSideColors=ROW.COLS )
dev.off()

## Determine which models have DRUG interactions
LM.drugint <- as.numeric(unlist( lapply( LM, function(x) any(grepl( "DRUG:|:DRUG", names(coef(x)) ))) ))

