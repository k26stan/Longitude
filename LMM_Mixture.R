## Run Mixtures of LMM Models on Time-Series Data ##
## Janssen Data using Multiple-Measures ##
## July 14, 2015 ##
## Kristopher Standish ##

## How many groups are Ideal?
## Which parameters should be estimated separately for each group?
## Which parameters are fixed and which are sampled from a random population?

## Should also providing estimates for initial values to achieve global maximum

## Load Packages
library(lcmm)
library(lattice)
library(gplots)
library(ggplot2)
library(nlme)

##############################################################
## LOAD DATA #################################################
##############################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Mac Paths
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150520_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_LME_Mixtures/",sep="" )
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Load Candidate Genotype Files
# COMP.l <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", sep="\t",header=T)
# COMP <- COMP.l[ which(!duplicated(COMP.l$SNP)), ]
# RAW <- read.table( "/Users/kstandis/Data/Burn/Results/20150413_GWAS/TEST_CND.raw", sep="",header=T)

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





## OLD FILTERING APPROACH!!! ##

# ## For Response Vectors, identify Missing Data
# RESP_COLS <- c(grep("CRP",names(TAB.l)),grep("DAS",names(TAB.l)),grep("JC",names(TAB.l)))
# RESP_COLS <- grep("DAS|CRP|JC",names(TAB.l)) ; names(RESP_COLS) <- grep("DAS|CRP|JC",names(TAB.l),value=T)
# MISSING <- list()
# len_MISS <- numeric(length(RESP_COLS))
# for ( resp in 1:length(RESP_COLS) ) {
# 	MISSING[[resp]] <- which(is.na(TAB.l[,RESP_COLS[resp]]))
# 	names(MISSING)[resp] <- colnames(TAB.l)[RESP_COLS[resp]]
# 	len_MISS[resp] <- length(MISSING[[resp]])
# }
# len_MISS
# length(intersect(MISSING[["DAS"]],MISSING[["SJC"]]))
#   # 492, 0, 452, 452, 452, 452
#   # Last 4 have all the same missing values
#   # 450 overlap between missing CRP and missing Joint data
# TEMP <- setdiff( MISSING[["DAS"]], MISSING[["SJC"]] )
# TEMP2 <- head(data.frame( TAB.l$DAS[TEMP], TAB.l$DAS[TEMP+1], TAB.l$DAS[TEMP+2] ))
# which( TEMP2[,1]==TEMP2[,2] & TEMP2[,2]==TEMP2[,3] )
#   # only the 3rd one.
#   # Which means, the DAS was changing after some missing values in the CRP data

# ## Remove Values with missing Joint Data
# RM <- MISSING[["SJC"]]
# TAB.2 <- TAB.l[-RM,]

# ## Check out Patients with only a FEW measurements or NONE after taking the drug
# LT3m <- which( table(TAB.2$FID)<3 )
# NO_DRUG <- which( table(TAB.2$FID,TAB.2$DRUG)[,2]==0 )
# RM2_SAMPS <- union( names(LT3m), names(NO_DRUG) )

# RM2 <- which(TAB.2$FID %in% RM2_SAMPS)
# TAB.3 <- TAB.2[-RM2,]


# dim(TAB.l)
# dim(TAB.2)
# dim(TAB.3)

# TAB <- TAB.3
# i <- 1


# hlme( DAS ~ WK+DRUG+SEX+RF_ACPA+DRUG:SEX+DRUG:RF_ACPA, random=~WK+AGE+BMI+DRUG:WK+DRUG:AGE+DRUG:BMI+AGE:SEX, subject="IID", ng=i, data=TAB, idiag=T)

##############################################################
## FCT: USE LCMM FOR MIXTURE MODELS ##########################
##############################################################

## Create Function to Run LCMM on Data
 # g[[i]] <- hlme( DAS ~ WK+DRUG+SEX+RF_ACPA+DRUG:SEX+DRUG:RF_ACPA, mixture=~WK+DRUG+RF_ACPA+DRUG:RF_ACPA, random=~WK+AGE+BMI+DRUG:WK+DRUG:AGE+DRUG:BMI+AGE:SEX, subject="IID", ng=i, data=DAT, idiag=T)
RUN_LCMM <- function( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS ) {
	MOD <- list()
	RUNTIME <- numeric( MAX_GRPS )
	for ( N_GRPS in 1:MAX_GRPS ) {
		start_time <- proc.time()
		tag <- paste( "Grp",N_GRPS,sep="_" )
		print(paste( "Running",N_GRPS,"of",MAX_GRPS) )
		if ( N_GRPS==1 ) {
			if ( COR==F ) {
				MOD[[tag]] <- hlme( fixed=FIXED, random=RANDOM, subject=SUBJECT, ng=N_GRPS, data=DATA, idiag=DIAG)
			}else{
				MOD[[tag]] <- hlme( fixed=FIXED, random=RANDOM, subject=SUBJECT, cor=AR(WK), ng=N_GRPS, data=DATA, idiag=DIAG)
			}
		}else{ # >1 Group
			if ( COR==F ) {
				MOD[[tag]] <- hlme( fixed=FIXED, mixture=MIXTURE, random=RANDOM, subject=SUBJECT, ng=N_GRPS, data=DATA, idiag=DIAG)
			}else{
				MOD[[tag]] <- hlme( fixed=FIXED, mixture=MIXTURE, random=RANDOM, subject=SUBJECT, cor=AR(WK), ng=N_GRPS, data=DATA, idiag=DIAG)
			}
		}
		RUNTIME[N_GRPS] <- (proc.time()-start_time)[3]
	}
	# Return List of Models
	COMPILE <- list( MOD=MOD, TIME=RUNTIME, FIXED=FIXED, MIXTURE=MIXTURE, RANDOM=RANDOM, COR=COR, DIAG=DIAG )
	return(COMPILE)
}
# # Test Function
# TAG <- "Test"
# DATA <- TAB
# FIXED <- as.formula( "DAS ~ DRUG" )
# MIXTURE <- as.formula( "~ 1" )
# RANDOM <- as.formula( "~ DRUG" )
# SUBJECT <- "IID"
# COR <- F
# DIAG <- F
# MAX_GRPS <- 3
# TEST <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

## Pull out Model Informations
PAR <- MIX[[1]]
PULL_MOD_INFO <- function(PAR) {
	MOD <- PAR$MOD
	n_grps <- length(MOD)
	## Get Convergence Status of Models
	CONV <- unlist(lapply( MOD, function(x) x$conv ))
	## Pull Information Criterion Info
	IC <- array(, dim=c(n_grps,4)) ; colnames(IC) <- c("nParam","logLik","AIC","BIC") ; rownames(IC) <- paste(1:n_grps,"Grps",sep="_")
	IC[,"nParam"] <- unlist(lapply( MOD, function(x) length(x$best) ))
	IC[,"logLik"] <- unlist(lapply( MOD, function(x) x$loglik ))
	IC[,"AIC"] <- 2*IC[,"nParam"] - 2*IC[,"logLik"]
	IC[,"BIC"] <- unlist(lapply( MOD, function(x) x$ns ))*IC[,"nParam"] - 2*IC[,"logLik"]
	## Pull Coefficients
	COEF <- lapply( MOD, function(x) summary(x)[,"coef"] )
	PVALS <- lapply( MOD, function(x) summary(x)[,"p-value"] )
	ALL_COEFS <- lapply( MOD, function(x) estimates(x,F) )
	# for ( i in 1:length(ALL_COEFS) ) {
	# 	names(ALL_COEFS[[i]])[grep("varcov",names(ALL_COEFS[[i]]))] <- c("intercept",sapply(strsplit(as.character(PAR$RANDOM[i])," + ",fixed=T),"[",2:length(grep("varcov",names(ALL_COEFS[[i]])))-1) )
	# }
	## Get Groups
	GRPS <- lapply( MOD, function(x) x$pprob )
	## Compile Output
	COMPILE <- list( CONV=CONV, IC=IC, COEF=COEF, PVALS=PVALS, ALL_COEFS=ALL_COEFS, GRPS=GRPS )
	return(COMPILE)
}
# SUMM <- lapply( MIX, function(x) PULL_MOD_INFO(x) )

##############################################################
## TEST DIFFERENT MIXTURE MODELS #############################
##############################################################

## Potential Covariates
 # WK - Random+Fixed
 # SEX - Fixed
 # AGE - Random
 # HT - Random
 # WT - Random
 # BMI - Random
 # DIS_DUR - Random
 # RF_ACPA - Fixed
 # ACPA - Fixed
 # RF - Fixed
 # DRUG - Fixed

## How many groups are Ideal?
## Which parameters should be estimated separately for each group?
## Which parameters are fixed and which are sampled from a random population?

TEST_MODELS <- function() {

	MIX <- list()

	############################################
	## Basic Models ##

	## Basic Model using DRUG as FIXED, MIXTURE, RANDOM
	print("Running M1")
	TAG <- "M1"
	DATA <- TAB
	FIXED <- as.formula( "DAS ~ DRUG" )
	MIXTURE <- as.formula( "~ DRUG" )
	RANDOM <- as.formula( "~ DRUG" )
	SUBJECT <- "IID"
	COR <- F
	DIAG <- F
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## Add PLAC as FIXED, MIXTURE, RANDOM
	print("Running M2")
	TAG <- "M2"
	FIXED <- as.formula( "DAS ~ DRUG+PLAC" )
	MIXTURE <- as.formula( "~ DRUG+PLAC" )
	RANDOM <- as.formula( "~ DRUG" )
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## Add DRUG*WK Interaction as FIXED, MIXTURE
	print("Running M3")
	TAG <- "M3"
	FIXED <- as.formula( "DAS ~ DRUG*WK+PLAC" )
	MIXTURE <- as.formula( "~ DRUG*WK+PLAC" )
	RANDOM <- as.formula( "~ DRUG" )
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## Add WK as RANDOM
	print("Running M4")
	TAG <- "M4"
	FIXED <- as.formula( "DAS ~ DRUG*WK+PLAC" )
	MIXTURE <- as.formula( "~ DRUG*WK+PLAC" )
	RANDOM <- as.formula( "~ DRUG+WK" )
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	############################################
	## BEST Models ##
	 # LME$M20e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )

	## COR=F,Diag=F
	print("Running B1a")
	TAG <- "B1a"
	DATA <- TAB
	FIXED <- as.formula( "DAS ~ DRUG*(WK+RF_ACPA)+PLAC" )
	MIXTURE <- as.formula( "~ DRUG" )
	RANDOM <- as.formula( "~ DRUG+WK" )
	SUBJECT <- "IID"
	COR <- F
	DIAG <- F
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## COR=T,Diag=F
	print("Running B1b")
	TAG <- "B1b"
	COR <- T
	DIAG <- F
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## COR=T,Diag=T
	print("Running B1c")
	TAG <- "B1c"
	COR <- T
	DIAG <- T
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## COR=F,Diag=T
	print("Running B1d")
	TAG <- "B1d"
	COR <- F
	DIAG <- T
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	############################################
	## More MIXTURE Parameters
	 # COR=T,Diag=T
	print("Running B2a")
	TAG <- "B2a"
	DATA <- TAB
	FIXED <- as.formula( "DAS ~ DRUG*(WK+RF_ACPA)+PLAC" )
	MIXTURE <- as.formula( "~ DRUG*WK" )
	RANDOM <- as.formula( "~ DRUG+WK" )
	COR <- T
	DIAG <- T
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## MIXTURE=DRUG*WK
	print("Running B2b")
	TAG <- "B2b"
	FIXED <- as.formula( "DAS ~ DRUG*(WK+RF_ACPA)+PLAC" )
	MIXTURE <- as.formula( "~ DRUG*(WK+RF_ACPA)+PLAC" )
	RANDOM <- as.formula( "~ DRUG+WK" )
	COR <- T
	DIAG <- T
	MAX_GRPS <- 3
	MIX[[TAG]] <- RUN_LCMM( DATA, FIXED, MIXTURE, RANDOM, SUBJECT, COR, DIAG, MAX_GRPS )

	## Return "MIX" Object
	return(MIX)
}
MIX <- TEST_MODELS()

##############################################################
## TEST DIFFERENT MIXTURE MODELS #############################
##############################################################


SUMM <- lapply( MIX, function(x) PULL_MOD_INFO(x) )

## Further Compile Data for Plotting
O.CONV <- unlist(lapply( SUMM, function(x) x$CONV ))
O.IC <- Reduce( rbind, lapply(SUMM,function(x)x$IC) )
O.GRPS <- as.numeric(gsub("_Grps","",rownames(O.IC) ))
O.MODS.1 <- unlist(lapply( SUMM, function(x)nrow(x$IC) ))
O.MODS <- rep(names(O.MODS.1),O.MODS.1)
O.TAGS <- paste( O.MODS,O.GRPS,sep="_" )
# rownames(O.IC) <- paste( rep(names(TEMP),TEMP), rownames(O.IC), sep="_" )

O <- data.frame( TAG=O.TAGS, MOD=O.MODS, GRPS=O.GRPS, CONV=O.CONV, O.IC )
# O <- O[ order(O$TAG), ]

GRPS.unq <- unique(O$GRP)
MODS.unq <- 1:length(unique(O$MOD)) ; names(MODS.unq) <- unique(O$MOD)
COLS.list <- c("firebrick1","chocolate1","gold1","chartreuse1","dodgerblue1","slateblue1")
COLS.4.list <- gsub("1","4",COLS.list)
COLS <- sample( COLS.list, length(GRPS.unq) )
COLS.4 <- gsub("1","4",COLS)
par(mfrow=c(2,2))
## Fit & Information
for ( met in c("logLik","AIC","BIC","nParam") ) {
	YLIM <- extendrange( O[,met] )
	XLIM <- c( 1,length(MODS.unq) )
	plot( 0,0,type="n", xaxt="n",xlab="",ylab=met,ylim=YLIM,xlim=XLIM )
	axis( 1, at=1:length(MODS.unq), label=names(MODS.unq), las=2 )
	# abline( v=which(!duplicated(O$MOD))-.5, col="black", lty=1)
	# text( 1:nrow(O), quantile(YLIM,1), label=O$GRPS )
	# text( 1:nrow(O), quantile(YLIM,.95), label=O$nParam )	
	for ( g in 1:length(GRPS.unq) ) {
		grp <- GRPS.unq[g]
		WHICH <- which(O$GRP==grp)
		points( MODS.unq[ as.character(O$MOD[WHICH]) ], type="o", O[WHICH,met], pch=c(21,1,4,13)[factor(O$CONV[WHICH])], bg=COLS[g],col=COLS.4[g] )
		# points( (1:nrow(O))[WHICH], O$AIC[WHICH], type="o", col="chartreuse2",pch=c(16,1,4,13)[factor(O$CONV[WHICH])] )
		if ( met=="logLik" ) {
			abline( v=which.max(O[WHICH,met]), col=COLS[g],lty=2 )	
		}else{
			abline( v=which.min(O[WHICH,met]), col=COLS[g],lty=2 )	
		}
	}
}









##############################################################
## END OF DOC ################################################
##############################################################
