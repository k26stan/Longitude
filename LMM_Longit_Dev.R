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
DATE <- "20150611"

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
## CALCULATE POPULATION METRICS ##############################
##############################################################
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
VAR.i.P <- aggregate( DAS ~ IID, data=TAB, var, subset=DRUG==0 )
MN.i.D <- aggregate( DAS ~ IID, data=TAB, mean, subset=DRUG==1 )
VAR.i.D <- aggregate( DAS ~ IID, data=TAB, var, subset=DRUG==1 )

## Calculate Between Patient (Population) Variance
TAU.0.P <- var( MN.i.P[,2] )
TAU.0.D <- var( MN.i.D[,2] )

## Calculate # Observations & B Values for Individuals
N.i.P <- aggregate( DAS ~ IID, data=TAB, length, subset=DRUG==0 )
B.i.P <- TAU.0.P / ( TAU.0.P + VAR.i.P[,2]/N.i.P[,2] )
N.i.D <- aggregate( DAS ~ IID, data=TAB, length, subset=DRUG==0 )
B.i.D <- TAU.0.D / ( TAU.0.D + VAR.i.D[,2]/N.i.D[,2] )

B.i.1.P <- TAU.0.P / ( TAU.0.P + VAR.i.P[,2] )
B.i.1.D <- TAU.0.D / ( TAU.0.D + VAR.i.D[,2] )

##############################################################
## CALCULATE Z SCORES ########################################
##############################################################

## Numerator
NUM.i.D <- (1-B.i.D)*MN.0.D + B.i.D*MN.i.D[,2]
NUM.i.P <- (1-B.i.P)*MN.0.P + B.i.P*MN.i.P[,2]
NUM <- NUM.i.D - NUM.i.P
NUM.frame <- data.frame( NUM.i.D, NUM.i.P, NUM )

## Denominator
V.i <- (1/(N.i.P[,2]+N.i.D[,2]))*(N.i.P[,2]*VAR.i.P[,2]+N.i.D[,2]*VAR.i.D[,2]) + (1/(N.0.P+N.0.D))*(N.0.P*TAU.0.P+N.0.D*TAU.0.D)
B.i.1 <- (1/(N.0.P+N.0.D))*(N.0.P*TAU.0.P+N.0.D*TAU.0.D) / V.i
B.i.n <- 


























##############################################################
## END OF DOC ################################################
##############################################################
