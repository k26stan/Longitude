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
DATE <- gsub("-","",Sys.Date())

## Mac Paths
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150520_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_LME_ModSel/",sep="" )
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
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

## Add PC's onto Table for later use
TAB.PC <- merge( TAB, FUL[,c("ID_2",paste("PC",1:3,sep=""))], by.x="IID",by.y="ID_2")

##############################################################
## CREATE PLOTTING FUNCTIONS #################################
##############################################################

## Plot Summary Statistics about Various LME Models
PLOT_MODS <- function( MODELS, tag ) {
	N.mod <- length(MODELS)
	Names.mod <- names(MODELS)
	Nums.mod <- gsub("[^0-9]", "",Names.mod)
	Swhich <- which(!duplicated(Nums.mod))
	COMP <- array( ,c(4,N.mod))
	rownames(COMP) <- c("logLik","AIC","BIC","SD")
	colnames(COMP) <- Names.mod
	## Extract Model Stats
	COMP["logLik",] <- unlist(lapply( MODELS, logLik ))
	COMP["AIC",] <- unlist(lapply( MODELS, AIC ))
	COMP["BIC",] <- unlist(lapply( MODELS, BIC ))
	COMP["SD",] <- unlist(lapply( MODELS, function(x) sd(resid(x)) ))
	## Plot
	COLS <- c("dodgerblue2","tomato2","slateblue3","chartreuse2")
	png( paste(PathToPlot,"/1-LME_ModAssess.",tag,".png",sep=""), height=700,width=2100,pointsize=32 )
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
# TEST <- PLOT_MODS( LME, "test" )

## Plot Presence & Significance of Covariates in LME Models
HEAT_MODS <- function( MODELS, tag ) {

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
	MODELS.terms.2 <- gsub( "I(ACPA == \"Positive\")TRUE","ACPAPositive",MODELS.terms, fixed=T )
	# MODELS.terms <- sort( MODELS.terms )
	MODELS.N.terms <- length(MODELS.terms)

	## Determine which Terms are in each model
	 # And extract P-Value for each Term
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
	rownames(MODELS.ARR.terms) <- rownames(MODELS.ARR.terms.p) <- MODELS.terms.2
	COLS.list <- c("firebrick1","chocolate1","gold1","springgreen2","steelblue2","slateblue3","black")
	COLS <- c( "firebrick1",colorRampPalette(COLS.list[3:7])(19) ) # colorRampPalette(COLS.list)(20)
	ROW.COLS <- rep("deepskyblue2",MODELS.N.terms) ; ROW.COLS[grep("DRUG",MODELS.terms)] <- "chartreuse2"
	png( paste(PathToPlot,"/1-LME_ModBuild.",tag,".png",sep=""), height=1200,width=1600,pointsize=30 )
	heatmap.2( MODELS.ARR.terms.p, main="Inclusion/Significance of Terms in Models",xlab="Model",ylab="Term",trace="none",scale="none",Colv=F,Rowv=F,dendrogram="none",col=COLS,margin=c(6,11),lhei=c(2,9),lwid=c(1,6),RowSideColors=ROW.COLS, colsep=Swhich-1 )
	dev.off()
	## Return P-Vals
	return( MODELS.ARR.terms.p )
}
# TEST <- HEAT_MODS( LME, "test" )

## Plot Predicted Values for Subset of Patients in specified LME Models
PLOT_PRED <- function( MODELS, DATA, N_SAMPS, tag ) {

	## Plotting Parameters
	 # Colors
	COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
	COLS <- colorRampPalette(sample(COLS.list))(N_SAMPS)
	if ( N_SAMPS==4 ) { COLS <- c("#6959CD", "#A75268", "#9ED527", "#EE7621") }
	COLS.1 <- sapply( COLS, function(x) colorRampPalette(c("white",x))(4)[2] )
	 # Which Samples
	Samps.which <- sample(unique(DATA$IID),N_SAMPS) # First.patient <- sample( 1:nrow(DATA), 1 )
	## Open Plot
	png( paste(PathToPlot,"2-LME_Mod_Pred.",tag,".png",sep=""), height=1000,width=1400,pointsize=30 )
	plot( 0,0,type="n", xlim=c(0,110),ylim=c(1,9),xlab="Weeks",ylab="DAS", main="Predicted vs Actual DAS Over Time" )
	abline( h=seq(0,10,1), lty=2,col="grey50",lwd=1 )
	## Plot Data
	for ( s in 1:N_SAMPS ) {
		# Specify Sample
		samp <- Samps.which[s]
		ind <- which(DATA$IID==samp)
		# Plotting Parameters
		PCHS <- c(21,16)[factor(DATA$DRUG[ind])]
		SWITCH <- DATA$WK[ind][which(DATA$DRUG[ind]==1)][1]
		# Plot Actual Data
		OBS <- DATA[ ind, "DAS" ]
		points( DATA[ind,"WK"], OBS, col=COLS[s], bg=COLS.1[s], type="o",pch=PCHS,lty=3,lwd=1 )
		PRED <- list()
		for ( m in 1:length(MODELS) ) {
			model <- names(MODELS)[m]
			PRED[[model]] <- predict(MODELS[[model]])[ind]
			points( DATA[ind,"WK"], PRED[[model]], col=COLS[s], type="l",pch=PCHS,lty=m,lwd=3 )
		}
		abline( v=jitter(SWITCH,amount=.5), col=COLS[s],lty=2,lwd=1 )
		# Add Sample Names
		text( 100, OBS[length(OBS)], label=samp, col=COLS[s],pos=4 )
	}
	dev.off()
	## Return List of Patients Plotted
	return(Samps.which)
}
# TEST <- PLOT_PRED( LME[c("M1e","M20e")], TAB, 4, "test" )

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
	COMP <- PLOT_MODS( LME, "Full" )

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
	COMP <- PLOT_MODS( LME, "Full" )

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
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Include Auto-correlation for WK within Patient
	 # Remove "d" as viable option
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 4"))
	LME$M4a <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M4b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M4c <- lme( fixed = DAS ~ DRUG*WK, random = ~ 1 | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	LME$M4e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | COUN/IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | COUN/IID) )
	LME$M4f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID/DRUG, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID/DRUG) )
	 # Choose Amongst Models
	anova(LME$M4a,LME$M4b,LME$M4c,LME$M4e,LME$M4f)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Include WK in Random Effects
	 # Revert to Order only in Correlation
	 # Remove "a", "c", and "f" as viable options
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 5"))
	LME$M5b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M5e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M5b,LME$M5e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Test WK in Correlation again
	 # (Keep WK in Random Effects)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 6"))
	LME$M6b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M6e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M6b,LME$M6e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Include Fixed PLAC Effect
	 # Revert to Order only in Correlation
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 7"))
	LME$M7b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M7e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M7b,LME$M7e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Test WK in Correlation again
	 # (Keep Fixed PLAC Effect)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 8"))
	LME$M8b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | IID) )
	LME$M8e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .8, form = ~WK | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M8b,LME$M8e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Switch to ARMA Correlation Structure
	 # Q = How many adjacent positions (pos vector) are used in Moving Average Correlation
	 # P = How many adjacent positions (pos vector) are used in Autoregression Correlation
	 # e.g., Moving Average w/ 2 adjacent positions: Q=2, P=0
	 # e.g., Autoregression w/ 3 adjacent positions: Q=0, P=3
	 # e.g., Hybrid using 3 for Moving Average & 2 for Autoregression: Q=3, P=2

	## Revert to Order only in Correlation, using ARMA (MA=2)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 9"))
	Q <- 2 ; P <- 0 
	LME$M9b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID) )
	LME$M9e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M9b,LME$M9e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Test WK in Correlation again, using ARMA (MA=2)
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 10"))
	LME$M10b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.8,.1,length.out=Q),seq(.8,.1,length.out=P)), q=Q,p=P, form = ~WK | IID) )
	LME$M10e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corARMA(value = c(seq(.8,.1,length.out=Q),seq(.8,.1,length.out=P)), q=Q,p=P, form = ~WK | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M10b,LME$M10e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## THIS FAILED!!!
	## Try hybrid AR & MA...using ARMA (MA=2,AR=1)
	 # STOP USING "WK" in Correlation...it doesn't work...
	# print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 11"))
	# Q <- 2 ; P <- 1
	# LME$M11b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID) )
	# LME$M11e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | COUN/IID) )
	#  # Choose Amongst Models
	# anova(LME$M11b,LME$M11e)
	# COMP <- PLOT_MODS( LME, "Full" )
	# HEAT <- HEAT_MODS( LME )

	## Test Clinical Covaraites in Model
	 # SEX, AGE, HT, WT, BMI, DIS_DUR, RF_ACPA, ACPA=Positive, RF

	## Include SEX as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 12"))
	LME$M12b <- lme( fixed = DAS ~ DRUG*(WK+SEX)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M12e <- lme( fixed = DAS ~ DRUG*(WK+SEX)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M12b,LME$M12e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Include AGE as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 13"))
	LME$M13b <- lme( fixed = DAS ~ DRUG*(WK+AGE)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M13e <- lme( fixed = DAS ~ DRUG*(WK+AGE)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M13b,LME$M13e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Include DIS_DUR as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 14"))
	LME$M14b <- lme( fixed = DAS ~ DRUG*(WK+DIS_DUR)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M14e <- lme( fixed = DAS ~ DRUG*(WK+DIS_DUR)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M14b,LME$M14e)
	COMP <- PLOT_MODS( LME, "Full" )
	HEAT <- HEAT_MODS( LME )

	## Include HT as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 15"))
	LME$M15b <- lme( fixed = DAS ~ DRUG*(WK+HT)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M15e <- lme( fixed = DAS ~ DRUG*(WK+HT)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M15b,LME$M15e)
	COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))>=5)], "Latter" )
	HEAT <- HEAT_MODS( LME )

	## Include WT as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 16"))
	LME$M16b <- lme( fixed = DAS ~ DRUG*(WK+WT)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M16e <- lme( fixed = DAS ~ DRUG*(WK+WT)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M16b,LME$M16e)
	COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))>=5)], "Latter" )
	HEAT <- HEAT_MODS( LME )

	## Include BMI as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 17"))
	LME$M17b <- lme( fixed = DAS ~ DRUG*(WK+BMI)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M17e <- lme( fixed = DAS ~ DRUG*(WK+BMI)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M17b,LME$M17e)
	COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))>=7)], "Latter" )
	HEAT <- HEAT_MODS( LME )

	## Include ACPA as Covariate (use I(ACPA=="Positive") )
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 18"))
	LME$M18b <- lme( fixed = DAS ~ DRUG*(WK+I(ACPA=="Positive"))+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M18e <- lme( fixed = DAS ~ DRUG*(WK+I(ACPA=="Positive"))+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M18b,LME$M18e)
	COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))>=7)], "Latter" )
	HEAT <- HEAT_MODS( LME )

	## Include RF as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 19"))
	LME$M19b <- lme( fixed = DAS ~ DRUG*(WK+RF)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M19e <- lme( fixed = DAS ~ DRUG*(WK+RF)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M19b,LME$M19e)
	COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))>=7)], "Latter" )
	HEAT <- HEAT_MODS( LME )

	## Include RF_ACPA as Covariate
	 # Revert to Order only in Correlation (using Autocorrlation, AR=1..."CAR1")
	print(paste(round(proc.time()-start_time,1)[3],"- Running Mod 20"))
	LME$M20b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME$M20e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | COUN/IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | COUN/IID) )
	 # Choose Amongst Models
	anova(LME$M20b,LME$M20e)
	COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))>=7)], "Latter" )
	HEAT <- HEAT_MODS( LME )

	return(LME)

}
LME <- LME_MOD(TAB)
COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))<=10)], "Start" )
COMP <- PLOT_MODS( LME[which(as.numeric(gsub("[a-z]","",names(LME),ignore.case=T))>=7)], "Latter" )
COMP <- PLOT_MODS( LME, "Full" )
HEAT <- HEAT_MODS( LME, "Full" )

## Plot Some Predicted Values of Various Models
WHICH_MODS <- c("M20e","M1e","M3e")
PRED <- PLOT_PRED( LME[WHICH_MODS], TAB, 4, "test" )
PRED <- PLOT_PRED( LME[WHICH_MODS], TAB, 5, "test" )

## Histogram of Coefficients (random Effects)
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
COLS <- sample(COLS.list)
COLS <- COLS.list[c(2,5,3,6)]
png( paste(PathToPlot,"3-LME_Coef_Dist.png",sep="/"), height=800,width=2000,pointsize=30 )
par(mfrow=c(1,4))
hist( coef(LME$M20b)[,"(Intercept)"], col=COLS[1], main="Coefficient Distribution", xlab="Intercept" )
hist( coef(LME$M20b)[,"DRUG"], col=COLS[2], main="Coefficient Distribution", xlab="Drug" )
hist( coef(LME$M20b)[,"WK"], col=COLS[3], main="Coefficient Distribution", xlab="Week" )
hist( coef(LME$M20b)[,"WK"]+coef(LME$M20b)[,"DRUG:WK"], col=COLS[4], main="Coefficient Distribution", xlab="Week (on Drug)" )
dev.off()

## Pairs/Distribution/Correlation b/n Beta Estimates
COLS.list <- c("firebrick3","chocolate2","gold2","springgreen2","steelblue2","slateblue3")
COLS <- sample(COLS.list)
COLS <- COLS.list[c(2,5,3,6)]
COLS <- COLS.list[c(1,3,4,5)]
B.1 <- coef(LME$M20b)[,"(Intercept)"]
B.2 <- coef(LME$M20b)[,"DRUG"]
B.3 <- coef(LME$M20b)[,"WK"]
B.4 <- coef(LME$M20b)[,"WK"] + coef(LME$M20b)[,"DRUG:WK"]
png( paste(PathToPlot,"3-LME_Coef_Dist_Cor.png",sep="/"), height=2000,width=2000,pointsize=36 )
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

## Correlation b/n DRUG & Intercept Random Effects Estimates
CORR.di <- unlist(lapply( LME[grep("b",names(LME))], function(x) cor(data.frame(coef(x)[,c("(Intercept)","DRUG")]))[1,2] ))
COLS <- colorRampPalette(COLS.list)(length(grep("b",names(LME))))
png( paste(PathToPlot,"/3-LME_Coef_Corr.png",sep=""), height=1200,width=1200, pointsize=26 )
par(mfrow=c(2,2))
plot( coef(LME$M1b)[,1], coef(LME$M1b)[,2], col=COLS[1],main="Beta Correlation: Model 1b",xlab="(Intercept)",ylab="DRUG" )
plot( coef(LME$M2b)[,1], coef(LME$M2b)[,2], col=COLS[2],main="Beta Correlation: Model 2b",xlab="(Intercept)",ylab="DRUG" )
plot( coef(LME$M20b)[,1], coef(LME$M20b)[,2], col=COLS[length(COLS)],main="Beta Correlation: Model 20b",xlab="(Intercept)",ylab="DRUG" )
barplot( CORR.di, main="Correlation: Intercept & Drug Betas",ylab="Correlation Coefficient", col=COLS,las=2)
dev.off()

## Correlation b/n Observed and Predicted Data using 3 different Models
CORR.op <- unlist(lapply( LME[grep("b",names(LME))], function(x) cor(TAB$DAS,predict(x)) ))
LIM <- c(0,10)
png( paste(PathToPlot,"/4-LME_Pred_v_Obs.png",sep=""), height=1200,width=1200, pointsize=26 )
par(mfrow=c(2,2))
 # Model 1
plot( data.frame( TAB$DAS, predict(LME$M1b) ), col=c(COLS[1],"grey20")[factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM, xlab="Observed",ylab="Predicted",main="Observed vs Predicted DAS - Model 1b" )
abline(lm( predict(LME$M1b) ~ TAB$DAS ), col=COLS[1],lty=2,lwd=3 )
abline( 0,1,lty=2 ) ; abline( h=seq(0,10,1),lty=3,col="grey50" ) ; abline( v=seq(0,10,1),lty=3,col="grey50" )
 # Model 1
plot( data.frame( TAB$DAS, predict(LME$M2b) ), col=c(COLS[2],"grey20")[factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM, xlab="Observed",ylab="Predicted",main="Observed vs Predicted DAS - Model 2b" )
abline(lm( predict(LME$M2b) ~ TAB$DAS ), col=COLS[2],lty=2,lwd=3 )
abline( 0,1,lty=2 ) ; abline( h=seq(0,10,1),lty=3,col="grey50" ) ; abline( v=seq(0,10,1),lty=3,col="grey50" )
 # Model 20
plot( data.frame( TAB$DAS, predict(LME$M20b) ), col=c(COLS[length(COLS)],"grey20")[factor(TAB$DRUG)],pch=3,xlim=LIM,ylim=LIM, xlab="Observed",ylab="Predicted",main="Observed vs Predicted DAS - Model 20b" )
abline(lm( predict(LME$M20b) ~ TAB$DAS ), col=COLS[length(COLS)],lty=2,lwd=3 )
abline( 0,1,lty=2 ) ; abline( h=seq(0,10,1),lty=3,col="grey50" ) ; abline( v=seq(0,10,1),lty=3,col="grey50" )
 # Correlation
barplot( CORR.op, main="Correlation: Observed & Predicted DAS",ylab="Correlation Coefficient", col=COLS,las=2)
dev.off()

##############################################################
## INCORPORATE VARIANTS ######################################
##############################################################

############ I FORGOT WHAT I WAS DOING!!!!!!!!!! ############
## CHANGING COLNAMES & STUFF SO I KNOW WHERE VARIANTS ARE ON THE PLOT #####

## Pull out Candidate SNPs
THRSH <- 1e-5 # 5e-5 # 5e-7 #
CAND.which <- which( COMP$P_Assoc < THRSH )
CAND.p <- COMP$P_Assoc[CAND.which]
CAND.id <- COMP$SNP[CAND.which] # paste( COMP$SNP[CAND.which], COMP$Allele[CAND.which], sep="_" )
CAND.loc <- paste( COMP$CHR[CAND.which], COMP$BP[CAND.which], sep="_" )
CAND.tag <- data.frame( ID=CAND.id, LOC=CAND.loc, P=CAND.p )

RAW.colnames <- unlist(lapply(strsplit( colnames(RAW), "_" ), function(x) paste(x[-length(x)],collapse="_") ))
RAW.colnames <- gsub(".",":",RAW.colnames,fixed=T )
RAW.colnames.which.X <- which( substr(RAW.colnames,1,2) %in% paste("X",1:22,sep="") )
RAW.colnames[RAW.colnames.which.X] <- gsub("X","",RAW.colnames[RAW.colnames.which.X],fixed=T )

CAND.intersect <- intersect( RAW.colnames, CAND.tag$ID )
CAND.rawcols <- which( RAW.colnames %in% CAND.intersect )
CAND.tagrows <- which( CAND.tag$ID %in% CAND.intersect )
RAW.cand <- RAW[,c(1,CAND.rawcols)]
CAND.tag.2 <- CAND.tag[ CAND.tagrows, ]
CAND.colnames <- colnames(RAW.cand)[2:ncol(RAW.cand)]

## Merge Candidate SNPs w/ Clinical Variables
MG <- merge( x=TAB.PC, y=RAW.cand, by="FID" )

## Test Model that includes each Variant Individually
LME.VAR <- list()
GWAS.p <- numeric( length(CAND.colnames) )
GWAS.id <- rep("", length(CAND.colnames) )
start_time <- proc.time()
for ( c in 1:length(CAND.colnames) ) {
# for ( c in 1:5 ) {
	cand <- CAND.colnames[c]
	cand.loc <- as.character(CAND.tag.2$LOC[c])
	TEMP_DAT <- MG[,c("IID","DAS","DRUG",cand,"WK","PLAC","RF_ACPA","PC1","PC2")]
	FORMULA <- as.formula( paste("DAS ~ DRUG*(",cand,"+WK+RF_ACPA)+PLAC+PC1+PC2",sep="" ) )
	# colnames(TEMP_DAT)[4] <- "cand"
	# LME.VAR[[cand.loc]] <- lme( DAS ~ DRUG*cand, random= ~ DRUG | IID, data=TEMP_DAT )
	# LME.VAR[[cand.loc]] <- lme( DAS ~ DRUG*cand+PC1+PC2, random= ~ DRUG | IID, data=TEMP_DAT )
	# LME.VAR[[cand.loc]] <- lme( fixed = DAS ~ (DRUG*WK)*cand+PLAC+PC1+PC2, random = ~ DRUG+WK | IID, data=TEMP_DAT )
	# LME.VAR[[cand.loc]] <- lme( fixed = DAS ~ DRUG*(cand+WK+RF_ACPA)+PLAC+PC1+PC2, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )
	LME.VAR[[cand.loc]] <- lme( fixed = FORMULA, random = ~ DRUG+WK | IID, data=MG, correlation=corCAR1(value = .5, form = ~1 | IID) )
	GWAS.p[c] <- CAND.tag.2$P[c]
	GWAS.id[c] <- as.character(CAND.tag.2$ID[c])
	names(GWAS.p)[c] <- cand.loc
	run_time <- round( proc.time()-start_time, 2)[3]
	if ( c%%2==0 ) { print(paste( "Done with",c,"of",length(CAND.colnames),"-",run_time )) }
}

## Pull out p-value for Variants
LME.VAR.p.cand <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable[grep("^rs",rownames(summary(x)$tTable)),"p-value"] ))
LME.VAR.p.cand.dr <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable[grep("^DRUG:rs",rownames(summary(x)$tTable)),"p-value"] ))
# LME.VAR.p.cand.wk <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable["^WK:rs","p-value"] ))
# LME.VAR.p.cand.dr.wk <- unlist(lapply( LME.VAR, function(x) summary(x)$tTable["DRUG:WK:cand","p-value"] ))
# P.comp <- data.frame( ID=GWAS.id, GWAS=GWAS.p, DIS=LME.VAR.p.cand, DR=LME.VAR.p.cand.dr, WK=LME.VAR.p.cand.wk, DR_WK=LME.VAR.p.cand.dr.wk )
P.comp <- data.frame( ID=GWAS.id, GWAS=GWAS.p, DIS=LME.VAR.p.cand, DR=LME.VAR.p.cand.dr )

## Plot P-Values for Variants vs Drug Response
COLS <- c("firebrick2","dodgerblue2","chartreuse2","purple2","gold1","chocolate2")
png( paste(PathToPlot,"V-LME_P-Vals.png",sep="/"), height=1000+nrow(P.comp),width=1600+5*nrow(P.comp),pointsize=30 )
plot( -log10(P.comp$GWAS), main="Variant Effect in Mixed Model", col=COLS[1],pch="+", ylim=c(0,9),xaxt="n",xlab="Variant",ylab="-log10(p)",cex=1.5 )
points( -log10(P.comp$DIS), col=COLS[2],pch="+",cex=1.5 )
points( -log10(P.comp$DR), col=COLS[3],pch="+",cex=1.5 )
# points( -log10(P.comp$WK), col=COLS[4],pch="+",cex=1.5 )
# points( -log10(P.comp$DR_WK), col=COLS[5],pch="+",cex=1.5 )
axis( 1, at=1:nrow(P.comp), label=rownames(P.comp),las=2,cex.axis=.7 )
abline( h=seq(0,15,1),lty=2,col="grey50")
abline( h=-log10(5e-8),lty=2,col="firebrick2")
abline( h=-log10(.05),lty=2,col="firebrick2")
legend( "bottomleft",legend=colnames(P.comp)[-1],col=COLS,pch="+",pt.cex=1.5,ncol=2,bg="white")
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