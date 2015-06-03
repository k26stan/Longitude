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
DATE <- "20150522"

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

## Take out Patients who left before getting DRUG
RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
RM.exit <- which( TAB.l$IID %in% RM.exit.id )
TAB <- TAB.l[-RM.exit,c(1:15,17)]

## Take out Patients from Golimumab Arm
RM.gol.id <- as.character( FUL$ID_2[which(FUL$GRP=="G")] )
RM.gol <- which( TAB$IID %in% RM.gol.id )
TAB <- TAB[ -RM.gol, ]

## Take out WK==0
TAB <- TAB[ which(TAB$WK!=0), ]

##############################################################
## LMM MODELS ################################################
##############################################################
	## When using ANOVA to compare models:
	 # Greater (less negative) logLik is better model
	 # Smaller (less positive) AIC/BIC is better model

LME.1 <- lme( fixed = DAS ~ DRUG, random = ~ 1 | IID, data=TAB )
LME.2 <- lme( fixed = DAS ~ DRUG, random = ~ DRUG | IID, data=TAB )
LME.3 <- lme( fixed = DAS ~ DRUG, random = ~ DRUG - 1 | IID, data=TAB )
LME.4 <- lme( fixed = DAS ~ DRUG - 1, random = ~ DRUG - 1 | IID, data=TAB )





























