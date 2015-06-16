## Do Longitudinal Analysis looking for Deviation from Baseline after Treatment ##
 # Use .../2013_Longitud_Cancer.pdf (Drescher, et al.) as Reference
## Janssen Data using Multiple-Measures ##
## June 11, 2015 ##
## Kristopher Standish ##

## Game Plan: Longitudinal Analysis as Template, but w/ Changes ##
 # Use only Patients who were on Placebo also
   # And Remove Baseline Measure??
 # Calculate Relevant Metrics for each Patient
   # Calculate Mean Baseline and Post-Treatment Scores (per person)
   # Calculate B Values
   # Calculate Variance Metrics
 # Calculate Z-Scores
   # Using Modified Formula

##############################################################
## LOAD DATA #################################################
##############################################################
library(nlme)
library(gplots)

## Set Date
DATE <- gsub("-","",Sys.Date())

## Mac Paths
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150226_Resp_v_Time.txt"
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150506_Resp_v_Time.txt"
# PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150512_Resp_v_Time.txt"
PathToData <- "/Users/kstandis/Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt"
PathToRawFiles <- "/Users/kstandis/Data/Burn/Data/Phenos/Raw_Files/"
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,sep="" )
dir.create( PathToPlot )

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
## FCT: CALCULATE Z-SCORES FOR EACH INDVIDUAL ################
##############################################################
# TAB <- NULL_TAB
Z_SCORE <- function( TAB ) {
	# TAB <- NULL_TAB

	##############################################################
	## CALCULATE POPULATION METRICS ##############################
	## Calculate Relevant Metrics for each Patient
	   # Calculate Mean Baseline and Post-Treatment Scores (per person)
	   # Calculate B Values
	   # Calculate Variance Metrics

	## Calculate Population Mean Values
	MN.0.P <- mean( TAB$DAS[ which(TAB$DRUG==0) ] )
	MN.0.D <- mean( TAB$DAS[ which(TAB$DRUG==1) ] )
	N.0.P <- length( which(TAB$DRUG==0) )
	N.0.D <- length( which(TAB$DRUG==1) )

	## Calculate Individual Mean/Variance Values
	MN.i.P <- aggregate( DAS ~ IID, data=TAB, mean, subset=DRUG==0 )
	# VAR.i.P <- aggregate( DAS ~ IID, data=TAB, var, subset=DRUG==0 ) ; VAR.i.P[which(is.na(VAR.i.P[,2])),2] <- 0
	VAR.i.P <- var( TAB$DAS[which(TAB$DRUG==0)] )
	MN.i.D <- aggregate( DAS ~ IID, data=TAB, mean, subset=DRUG==1 )
	# VAR.i.D <- aggregate( DAS ~ IID, data=TAB, var, subset=DRUG==1 ) ; VAR.i.D[which(is.na(VAR.i.D[,2])),2] <- 0
	VAR.i.D <- var( TAB$DAS[which(TAB$DRUG==1)] )

	## Calculate Between Patient (Population) Variance
	TAU.0.P <- var( MN.i.P[,2] )
	TAU.0.D <- var( MN.i.D[,2] )

	## Calculate # Observations & B Values for Individuals
	N.i.P <- aggregate( DAS ~ IID, data=TAB, length, subset=DRUG==0 )
	# B.i.P <- TAU.0.P / ( TAU.0.P + VAR.i.P[,2]/N.i.P[,2] )
	B.i.P <- TAU.0.P / ( TAU.0.P + VAR.i.P/N.i.P[,2] )
	N.i.D <- aggregate( DAS ~ IID, data=TAB, length, subset=DRUG==1 )
	# B.i.D <- TAU.0.D / ( TAU.0.D + VAR.i.D[,2]/N.i.D[,2] )
	B.i.D <- TAU.0.D / ( TAU.0.D + VAR.i.D/N.i.D[,2] )

	# B.i.1.P <- TAU.0.P / ( TAU.0.P + VAR.i.P[,2] )
	B.i.1.P <- TAU.0.P / ( TAU.0.P + VAR.i.P )
	# B.i.1.D <- TAU.0.D / ( TAU.0.D + VAR.i.D[,2] )
	B.i.1.D <- TAU.0.D / ( TAU.0.D + VAR.i.D )

	# pairs( data.frame( B.i.P, B.i.1.P, B.i.1, B.i.1.D, B.i.D, B.i.1*(B.i.D+B.i.P)/2 ) )

	##############################################################
	## CALCULATE Z SCORES ########################################

	## Numerator
	NUM.i.D <- (1-B.i.D)*MN.0.D + B.i.D*MN.i.D[,2]
	NUM.i.P <- (1-B.i.P)*MN.0.P + B.i.P*MN.i.P[,2]
	NUM <- NUM.i.D - NUM.i.P
	NUM.frame <- data.frame( NUM.i.D, NUM.i.P, NUM )

	## Denominator
	N.i.sum <- (N.i.P[,2]+N.i.D[,2])
	TAU.0.wtd <- ( N.0.P*TAU.0.P + N.0.D*TAU.0.D ) / (N.0.P+N.0.D)
	# VAR.i.wtd <- (N.i.P[,2]*VAR.i.P[,2]+N.i.D[,2]*VAR.i.D[,2]) / N.i.sum
	VAR.i.wtd <- (N.i.P[,2]*VAR.i.P+N.i.D[,2]*VAR.i.D) / N.i.sum
	# VAR.i.wtd.2 <- (N.i.P[,2]*VAR.i.P[,2]+N.i.D[,2]*VAR.i.D[,2]) / N.i.sum^2
	V.i <- VAR.i.wtd + TAU.0.wtd
	B.i.1 <- TAU.0.wtd / V.i
	B.i.1.wtd <- (N.i.P[,2]*B.i.1.P+N.i.D[,2]*B.i.1.D) / N.i.sum
	# B.i.1.2 <- TAU.0.wtd / ( VAR.i.wtd + TAU.0.wtd )
	B.i.n.wtd <- (N.i.P[,2]*B.i.P+N.i.D[,2]*B.i.D) / N.i.sum
	# B.i.n <- TAU.0.wtd / ( TAU.0.wtd + VAR.i.wtd/N.i.sum )
	# B.i.n.2 <- ( N.i.sum*B.i.1 ) / ( N.i.sum*B.i.1 + (1-B.i.1) )
	# DENOM <- sqrt( V.i * ( 1 - B.i.1 * B.i.n ) )
	DENOM <- sqrt( V.i * ( 1 - B.i.1.wtd * B.i.n.wtd ) )
	DENOM.2 <- sqrt( V.i * ( 1 - B.i.1 * (B.i.D+B.i.P)/2 ) )
	DENOM.3 <- sqrt( V.i * ( 1 - B.i.1 * sqrt(B.i.D*B.i.P) ) )
	# DENOM <- sqrt( V.i )
	# DENOM <- sqrt( VAR.i.wtd )

	## Z-Score
	Z.i <- NUM / DENOM
	Z.i.2 <- NUM / DENOM.2
	Z.i.3 <- NUM / DENOM.3

	## Histogram of Z-Scores
	LIM <- range( Z.i )
	XLIM <- c( floor(LIM[1]),ceiling(LIM[2]) )
	BRKS <- seq( XLIM[1],XLIM[2], .25 )
	COLS <- colorRampPalette(c("cadetblue2","aquamarine4"))(4)[3]
	hist( Z.i, xlim=XLIM,breaks=BRKS, col=COLS,main="Z-Score Distribution",xlab="Patient Z-Scores" )
	MN.Z <- mean(Z.i) ; SD.Z <- sd(Z.i)
	SHAP.P <- shapiro.test(Z.i)$p.value
	text( quantile(XLIM,.05),length(Z.i)/20, label=paste("Mean =",round(MN.Z,2)),pos=4 )
	text( quantile(XLIM,.05),length(Z.i)/25, label=paste("SD =",round(SD.Z,2)),pos=4 )
	text( quantile(XLIM,.05),length(Z.i)/35, label=paste("Shap.P =",round(SHAP.P,2)),pos=4 )

	## Return Z Scores
	return(Z.i)

} # Close Z-Score Function
TEST <- Z_SCORE(TAB[1:100,])


##############################################################
## CHECK VARIANCE OF Z-SCORES UNDER VARIOUS SCENARIOS ########
##############################################################

N.patients <- 1000
N.measures <- 50
MN.DAS <- mean(TAB$DAS) ; SD.DAS <- sd(TAB$DAS)
NULL_TAU <- 2 ; NULL_VAR <- 1
NULL_IDS <- rep(1:N.patients,rep(N.measures,N.patients))
NULL_DRUG <- c(replicate( N.patients, sort(c(0,1,sample(0:1,N.measures-2,replace=T))) ))
NULL_DAS <- c(replicate( N.patients, rnorm(N.measures,rnorm(1,MN.DAS,NULL_TAU),NULL_VAR) ))
NULL_TAB <- data.frame( IID=NULL_IDS, DRUG=NULL_DRUG, DAS=NULL_DAS )
NULL_TEST <- Z_SCORE(NULL_TAB)

## Check
N.patients <- 100
CHECK_VAR <- function(N.patients) {
	NULL_TEST <- list()
	N.measures.list <- sort( rep(seq(5,50,5),5) )
	for ( n in 1:length(N.measures.list) ) {
		N.measures <- N.measures.list[n]
		tag <- paste(N.measures,"meas",n,sep="_")
		## Create Data
		MN.DAS <- mean(TAB$DAS) ; SD.DAS <- sd(TAB$DAS)
		NULL_TAU <- 2 ; NULL_VAR <- 1
		NULL_IDS <- rep(1:N.patients,rep(N.measures,N.patients))
		NULL_DRUG <- c(replicate( N.patients, sort(c(0,1,sample(0:1,N.measures-2,replace=T))) ))
		NULL_DAS <- c(replicate( N.patients, rnorm(N.measures,rnorm(1,MN.DAS,NULL_TAU),NULL_VAR) ))
		NULL_TAB <- data.frame( IID=NULL_IDS, DRUG=NULL_DRUG, DAS=NULL_DAS )
		## Run Test to Get Distribution
		NULL_TEST[[tag]] <- Z_SCORE(NULL_TAB)
		if ( n%%10==0 ) { print(paste( "Done with",n,"of",length(N.measures.list) )) }
	}
	return(NULL_TEST)
}
NULL_TEST <- CHECK_VAR(N.patients) ; plot( N.measures.list, unlist(lapply( NULL_TEST_PATS[[p]], sd )) )


## Collect Z-Score Distribution over different # of Measurements for Different Population Sizes
NULL_TEST_PATS <- list()
N.patients.list <- seq(20,200,20)
for ( p in 1:length(N.patients.list) ) {
	N.patients <- N.patients.list[p]
	tag <- paste(N.patients,"pats",p,sep="_")
	NULL_TEST_PATS[[tag]] <- CHECK_VAR(N.patients)
}


par(mfrow=c(2,5))
COLS.list <- c("slateblue1","steelblue1","springgreen1","gold1","chocolate1","firebrick1")
COLS <- colorRampPalette(COLS.list)(length(N.patients.list))
plot( 0,0,type="n",xlim=c(0,50),ylim=c(0,.4) )
for ( p in 1:length(N.patients.list) ) {
	# par(ask=T)
	points( N.measures.list, unlist(lapply( NULL_TEST_PATS[[p]], sd )), col=COLS[p] )
}
##### !!!!!!!!!!!! SD OF Z-SCORE DISTRIBUTION DEPENDS ON NUMBER OF MEASUREMENTS PER PERSON !!!!!!!!!! ######

##############################################################
## CALCULATE Z SCORES FOR INDIVIDUALS ########################
##############################################################

## Only Patients who were in Placebo Arm of Study to Begin
Z.plac <- Z_SCORE(TAB)























##############################################################
## END OF DOC ################################################
##############################################################
