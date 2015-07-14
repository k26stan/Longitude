## Do Bayesian Analysis w/ Mixed Model & Longitudinal Data ##
## Janssen Data using Multiple-Measures ##
## May 22, 2015 ##
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
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_LME_CorStruct/",sep="" )
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Load Candidate Genotype Files
TAB.PC <- merge( TAB.l, FUL[,c("ID_2",paste("PC",1:3,sep=""))], by.x="IID",by.y="ID_2")
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
## TESTING CorStruct #########################################
##############################################################

## Which Correlation Structure?
 # corCAR1, corAR1, corARMA
## What Value? (strength of correlation b/n adjacent measurements)
 # 0 < value < 1
## Which variables to include in Form?

TEST <- list()
TEST$M1 <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB )	
TEST$M2a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~1 | IID) )
TEST$M2b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~1 | IID/DRUG) )
TEST$M2c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~WK | IID) )
TEST$M2d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~WK | IID/DRUG) )
TEST$M3a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )	
TEST$M3b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )	
TEST$M3c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~WK | IID) )	
TEST$M3d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~WK | IID/DRUG) )	
Q <- 1 ; P <- 0
TEST$M4a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST$M4b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 2 ; P <- 0
TEST$M4c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST$M4d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 3 ; P <- 0
TEST$M4e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST$M4f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 1 ; P <- 0
TEST$M5a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST$M5b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 2 ; P <- 0
TEST$M5c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST$M5d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 4 ; P <- 0
TEST$M5e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST$M5f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 6 ; P <- 0
TEST$M5g <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST$M5h <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 0 ; P <- 1
TEST$M6a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~1 | IID, q=Q,p=P) )
TEST$M6b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~1 | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 2
TEST$M6c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~1 | IID, q=Q,p=P) )
TEST$M6d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~1 | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 3
TEST$M6e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~1 | IID, q=Q,p=P) )
TEST$M6f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~1 | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 1
TEST$M7a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST$M7b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 2
TEST$M7c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST$M7d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 3
TEST$M7e <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST$M7f <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 1 ; P <- 1
TEST$M8a <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST$M8b <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 1 ; P <- 1
TEST$M8c <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST$M8d <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
 # Test Models vs One Another
anova( TEST$M1, TEST$M2a, TEST$M2b, TEST$M3a, TEST$M3b, TEST$M4a, TEST$M4b, TEST$M5a, TEST$M5b )
 # Plot logLik, BIC
par(mfrow=c(1,2))
plot( unlist(lapply(TEST,logLik)), type="o",xaxt="n" ) ; axis(1, at=1:length(TEST), label=names(TEST),las=2) ; abline( v=which.max(unlist(lapply(TEST,logLik))), lty=2 )
plot( unlist(lapply(TEST,BIC)), type="o",xaxt="n" ) ; axis(1, at=1:length(TEST), label=names(TEST),las=2) ; abline( v=which.min(unlist(lapply(TEST,BIC))), lty=2 )
points( unlist(lapply(TEST,AIC)), type="o",lty=2, ) ; axis(1, at=1:length(TEST), label=names(TEST),las=2) ; abline( v=which.min(unlist(lapply(TEST,BIC))), lty=2 )

GRK <- lapply( TEST, function(x) x$modelStruct$corStruct)
plot( 0,0,type="n",xlim=c(1,length(GRK)),ylim=c(-1,1),xaxt="n" ) ; axis(1, at=1:length(TEST), label=names(TEST),las=2)
for ( i in 1:length(GRK) ) {
	Num <- length(GRK[[i]])
	points( rep(i,Num), coef(GRK[[i]],unconstrained=F) )
}
# Diagnostic Plots
plot( ACF(TEST$M2a, maxLag = 10, resType = "n"), alpha = 0.01 )
# plot( Variogram( TEST$M5c, form = ~ WK ) )
# plot( Variogram( TEST$M5c, form = ~ WK, maxDist=50 ) )
# plot( Variogram( TEST$M5c, form = ~ WK, maxDist=50, resType="n", robust=T ) )



TEST.p <- list()
TEST.p$M1 <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB )	
TEST.p$M2a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~1 | IID) )
TEST.p$M2b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~1 | IID/DRUG) )
TEST.p$M2c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~WK | IID) )
TEST.p$M2d <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~WK | IID/DRUG) )
TEST.p$M3a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )	
TEST.p$M3b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )	
TEST.p$M3c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~WK | IID) )	
TEST.p$M3d <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~WK | IID/DRUG) )	
Q <- 1 ; P <- 0
TEST.p$M4a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST.p$M4b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 2 ; P <- 0
TEST.p$M4c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST.p$M4d <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 3 ; P <- 0
TEST.p$M4e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST.p$M4f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 1 ; P <- 0
TEST.p$M5a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.p$M5b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 2 ; P <- 0
TEST.p$M5c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.p$M5d <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 4 ; P <- 0
TEST.p$M5e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.p$M5f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 6 ; P <- 0
TEST.p$M5g <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.p$M5h <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 0 ; P <- 1
TEST.p$M6a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST.p$M6b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 2
TEST.p$M6c <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST.p$M6d <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 3
TEST.p$M6e <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST.p$M6f <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 1 ; P <- 1
TEST.p$M7a <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.p$M7b <- lme( fixed = DAS ~ DRUG*WK+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
 # Test Models vs One Another
anova( TEST.p$M1, TEST.p$M2a, TEST.p$M2b, TEST.p$M3a, TEST.p$M3b, TEST.p$M4a, TEST.p$M4b, TEST.p$M5a, TEST.p$M5b )
 # Plot logLik, BIC
par(mfrow=c(1,2))
plot( unlist(lapply(TEST.p,logLik)), type="o",xaxt="n" ) ; axis(1, at=1:length(TEST.p), label=names(TEST.p),las=2) ; abline( v=which.max(unlist(lapply(TEST.p,logLik))), lty=2 )
plot( unlist(lapply(TEST.p,BIC)), type="o",xaxt="n" ) ; axis(1, at=1:length(TEST.p), label=names(TEST.p),las=2) ; abline( v=which.min(unlist(lapply(TEST.p,BIC))), lty=2 )
points( unlist(lapply(TEST.p,AIC)), type="o",lty=2, ) ; axis(1, at=1:length(TEST.p), label=names(TEST.p),las=2) ; abline( v=which.min(unlist(lapply(TEST.p,BIC))), lty=2 )

GRK <- lapply( TEST.p, function(x) x$modelStruct$corStruct)
plot( 0,0,type="n",xlim=c(1,length(GRK)),ylim=c(-1,1),xaxt="n" ) ; axis(1, at=1:length(TEST), label=names(TEST),las=2)
for ( i in 1:length(GRK) ) {
	Num <- length(GRK[[i]])
	points( rep(i,Num), coef(GRK[[i]],unconstrained=F) )
}
# Diagnostic Plots
plot( ACF(TEST.p$M2a, maxLag = 10, resType = "n"), alpha = 0.01 )
plot( ACF(TEST.p$M4e, maxLag = 10, resType = "n"), alpha = 0.01 )
# plot( Variogram( TEST.p$M5c, form = ~ WK ) )
# plot( Variogram( TEST.p$M5c, form = ~ WK, maxDist=50 ) )
# plot( Variogram( TEST.p$M5c, form = ~ WK, maxDist=50, resType="n", robust=T ) )

N.plot <- 5
IDS <- as.character(sample(unique( TAB$IID), N.plot ))
COLS.resid <- COLS[round(seq(1,length(COLS),length.out=N.plot),0)]
plot( 0,0,type="n",xlim=c(0,100),ylim=c(-2,2),xlab="WK",ylab="Resids")
abline( h=seq(-6,6,1),lty=rep(c(2,1,2),c(6,1,6)) )
for ( i in 1:N.plot ) {
	id <- IDS[i]
	which.id <- which( TAB$IID==id )
	points( TAB$WK[which.id], resid(TEST.p$M4e)[which.id], type="o",lty=2, col=COLS.resid[i],lwd=2 )
}

BRKS <- seq( -10,10,.5 ) ; XLIM <- c(-2.5,2.5)
par(mfrow=c(4,4))
for ( i in 1:16 ) {
	which.wk <- which( TAB$WK==WKS[i] )
	MEAN <- round( mean( resid(TEST.p$M4e)[which.wk] ), 3)
	SD <- round( sd( resid(TEST.p$M4e)[which.wk] ), 3)
	hist( resid(TEST.p$M4e)[which.wk], breaks=BRKS,xlim=XLIM, col=COLS[30],main=paste("MN=",MEAN,": SD=",SD,sep="") )
}
quartz()
par(mfrow=c(4,4))
for ( i in 1:16 ) {
	which.wk <- which( TAB$WK==WKS[i] )
	MEAN <- round( mean( resid(TEST.p$M2a)[which.wk] ), 3)
	SD <- round( sd( resid(TEST.p$M2a)[which.wk] ), 3)
	hist( resid(TEST.p$M2a)[which.wk], breaks=BRKS,xlim=XLIM, col=COLS[60],main=paste("MN=",MEAN,": SD=",SD,sep="") )
}
quartz()
par(mfrow=c(4,4))
for ( i in 1:16 ) {
	which.wk <- which( TAB$WK==WKS[i] )
	MEAN <- round( mean( resid(TEST.p$M1)[which.wk] ), 3)
	SD <- round( sd( resid(TEST.p$M1)[which.wk] ), 3)
	hist( resid(TEST.p$M1)[which.wk], breaks=BRKS,xlim=XLIM, col=COLS[90],main=paste("MN=",MEAN,": SD=",SD,sep="") )
}


####################### ADD RF_ACPA ################################
TEST.pr <- list()
TEST.pr$M1 <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB )	
TEST.pr$M2a <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~1 | IID) )
TEST.pr$M2b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~1 | IID/DRUG) )
TEST.pr$M2c <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~WK | IID) )
TEST.pr$M2d <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corAR1(value = .5, form = ~WK | IID/DRUG) )
TEST.pr$M3a <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID) )	
TEST.pr$M3b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~1 | IID/DRUG) )	
TEST.pr$M3c <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~WK | IID) )	
TEST.pr$M3d <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corCAR1(value = .5, form = ~WK | IID/DRUG) )	
Q <- 1 ; P <- 0
TEST.pr$M4a <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST.pr$M4b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 2 ; P <- 0
TEST.pr$M4c <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST.pr$M4d <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 3 ; P <- 0
TEST.pr$M4e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID ) )
TEST.pr$M4f <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~1 | IID/DRUG ) )
Q <- 1 ; P <- 0
TEST.pr$M5a <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.pr$M5b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 2 ; P <- 0
TEST.pr$M5c <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.pr$M5d <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 4 ; P <- 0
TEST.pr$M5e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.pr$M5f <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 6 ; P <- 0
TEST.pr$M5g <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.pr$M5h <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
Q <- 0 ; P <- 1
TEST.pr$M6a <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST.pr$M6b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 2
TEST.pr$M6c <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST.pr$M6d <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 0 ; P <- 3
TEST.pr$M6e <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID, q=Q,p=P) )
TEST.pr$M6f <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA( value=c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P) )
Q <- 1 ; P <- 1
TEST.pr$M7a <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID ) )
TEST.pr$M7b <- lme( fixed = DAS ~ DRUG*(WK+RF_ACPA)+PLAC, random = ~ DRUG+WK | IID, data=TAB, correlation=corARMA(value = c(seq(.5,.1,length.out=Q),seq(.5,.1,length.out=P)), q=Q,p=P, form = ~WK | IID/DRUG ) )
 # Test Models vs One Another
anova( TEST.pr$M1, TEST.pr$M2a, TEST.pr$M2b, TEST.pr$M3a, TEST.pr$M3b, TEST.pr$M4a, TEST.pr$M4b, TEST.pr$M5a, TEST.pr$M5b )
 # Plot logLik, BIC
par(mfrow=c(1,2))
plot( unlist(lapply(TEST.pr,logLik)), type="o",xaxt="n" ) ; axis(1, at=1:length(TEST.pr), label=names(TEST.pr),las=2) ; abline( v=which.max(unlist(lapply(TEST.pr,logLik))), lty=2 )
plot( unlist(lapply(TEST.pr,BIC)), type="o",xaxt="n" ) ; axis(1, at=1:length(TEST.pr), label=names(TEST.pr),las=2) ; abline( v=which.min(unlist(lapply(TEST.pr,BIC))), lty=2 )
points( unlist(lapply(TEST.pr,AIC)), type="o",lty=2, ) ; axis(1, at=1:length(TEST.pr), label=names(TEST.pr),las=2) ; abline( v=which.min(unlist(lapply(TEST.pr,BIC))), lty=2 )

GRK <- lapply( TEST.pr, function(x) x$modelStruct$corStruct)
plot( 0,0,type="n",xlim=c(1,length(GRK)),ylim=c(-1,1),xaxt="n" ) ; axis(1, at=1:length(TEST), label=names(TEST),las=2)
for ( i in 1:length(GRK) ) {
	Num <- length(GRK[[i]])
	points( rep(i,Num), coef(GRK[[i]],unconstrained=F) )
}
# Diagnostic Plots
plot( ACF(TEST.pr$M2a, maxLag = 10, resType = "n"), alpha = 0.01 )
plot( ACF(TEST.pr$M4e, maxLag = 10, resType = "n"), alpha = 0.01 )
# plot( Variogram( TEST.pr$M5c, form = ~ WK ) )
# plot( Variogram( TEST.pr$M5c, form = ~ WK, maxDist=50 ) )
# plot( Variogram( TEST.pr$M5c, form = ~ WK, maxDist=50, resType="n", robust=T ) )

N.plot <- 5
IDS <- as.character(sample(unique( TAB$IID), N.plot ))
COLS.resid <- COLS[round(seq(1,length(COLS),length.out=N.plot),0)]
plot( 0,0,type="n",xlim=c(0,100),ylim=c(-2,2),xlab="WK",ylab="Resids")
abline( h=seq(-6,6,1),lty=rep(c(2,1,2),c(6,1,6)) )
for ( i in 1:N.plot ) {
	id <- IDS[i]
	which.id <- which( TAB$IID==id )
	points( TAB$WK[which.id], resid(TEST.pr$M2a)[which.id], type="o",lty=2, col=COLS.resid[i],lwd=2 )
}

BRKS <- seq( -10,10,.5 ) ; XLIM <- c(-2.5,2.5)
par(mfrow=c(4,4))
for ( i in 1:16 ) {
	which.wk <- which( TAB$WK==WKS[i] )
	MEAN <- round( mean( resid(TEST.pr$M4c)[which.wk] ), 3)
	SD <- round( sd( resid(TEST.pr$M4c)[which.wk] ), 3)
	hist( resid(TEST.pr$M4c)[which.wk], breaks=BRKS,xlim=XLIM, col=COLS[30],main=paste("MN=",MEAN,": SD=",SD,sep="") )
}
quartz()
par(mfrow=c(4,4))
for ( i in 1:16 ) {
	which.wk <- which( TAB$WK==WKS[i] )
	MEAN <- round( mean( resid(TEST.pr$M2a)[which.wk] ), 3)
	SD <- round( sd( resid(TEST.pr$M2a)[which.wk] ), 3)
	hist( resid(TEST.pr$M2a)[which.wk], breaks=BRKS,xlim=XLIM, col=COLS[60],main=paste("MN=",MEAN,": SD=",SD,sep="") )
}
quartz()
par(mfrow=c(4,4))
for ( i in 1:16 ) {
	which.wk <- which( TAB$WK==WKS[i] )
	MEAN <- round( mean( resid(TEST.pr$M1)[which.wk] ), 3)
	SD <- round( sd( resid(TEST.pr$M1)[which.wk] ), 3)
	hist( resid(TEST.pr$M1)[which.wk], breaks=BRKS,xlim=XLIM, col=COLS[90],main=paste("MN=",MEAN,": SD=",SD,sep="") )
}

















COR.2 <- list()
for ( v in seq(1,9,2) ) {
	COR.2[[paste("M",v,sep=".")]] <- lme( fixed = DAS ~ DRUG*WK, random = ~ DRUG | IID/DRUG, data=TAB, correlation=corCAR1(value = v/10, form = ~WK | IID/DRUG, fixed=T) )	
	print(paste( "Done with",v ))
}
unlist(lapply( COR.2, logLik ))
unlist(lapply( COR.2, BIC ))


COLS.list <- c("firebrick1","chocolate1","gold1","springgreen2","steelblue2","slateblue3","black")
COLS <- colorRampPalette(rev(COLS.list))(100)

ar.6a <- corAR1(value = .6, form = ~1 | IID/DRUG)
ar.6a <- Initialize( ar.6a, data=TAB )
ar.6a.mat <- corMatrix( ar.6a )
quartz() ; heatmap.2( car.6a.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( car.6a.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

ar.6b <- corAR1(value = .6, form = ~WK | IID/DRUG)
ar.6b <- Initialize( ar.6b, data=TAB )
ar.6b.mat <- corMatrix( ar.6b )
quartz() ; heatmap.2( car.6a.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( car.6a.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )


car.6a <- corCAR1(value = .6, form = ~1 | IID/DRUG)
car.6a <- Initialize( car.6a, data=TAB )
car.6a.mat <- corMatrix( car.6a )
quartz() ; heatmap.2( car.6a.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( car.6a.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

car.6b <- corCAR1(value = .6, form = ~WK | IID/DRUG)
car.6b <- Initialize( car.6b, data=TAB )
car.6b.mat <- corMatrix( car.6b )
quartz() ; heatmap.2( car.6b.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( car.6b.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

## Q = How many adjacent positions (pos vector) are used in Moving Average Correlation
## P = How many adjacent positions (pos vector) are used in Autoregression Correlation
 # e.g., Moving Average w/ 2 adjacent positions: Q=2, P=0
 # e.g., Autoregression w/ 3 adjacent positions: Q=0, P=3
 # e.g., Hybrid using 3 for Moving Average & 2 for Autoregression: Q=3, P=2
       # for Q<=3, it's hybrid ; for Q>3, it's Autoregressive
Q <- 2 ; P <- 0 # No autoregressive
arma.6b <- corARMA( value=rep(c(.6,.6),c(Q,P)), form = ~1 | IID/DRUG, q=Q,p=P, fixed=T)
arma.6b <- Initialize( arma.6b, data=TAB )
arma.6b.mat <- corMatrix( arma.6b )
quartz() ; heatmap.2( arma.6b.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( arma.6b.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

Q <- 8 ; P <- 1
arma.6c <- corARMA( value=c(seq(.6,.1,length.out=Q),seq(.6,.1,length.out=P)), form = ~WK | IID/DRUG, q=Q,p=P, fixed=T)
arma.6c <- Initialize( arma.6c, data=TAB )
arma.6c.mat <- corMatrix( arma.6c )
quartz() ; heatmap.2( arma.6c.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( arma.6c.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

Q <- 8 ; P <- 1
arma.6d <- corARMA( value=rep(c(.6,.6),c(Q,P)), form = ~WK | IID/DRUG, q=Q,p=P, fixed=T)
arma.6d <- Initialize( arma.6d, data=TAB )
arma.6d.mat <- corMatrix( arma.6d )
quartz() ; heatmap.2( arma.6d.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( arma.6d.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

Q <- 0 ; P <- 1
arma.6e <- corARMA( value=rep(c(.9,.9),c(Q,P)), form = ~WK | IID/DRUG, q=Q,p=P, fixed=T)
arma.6e <- Initialize( arma.6e, data=TAB )
arma.6e.mat <- corMatrix( arma.6e )
quartz() ; heatmap.2( arma.6e.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( arma.6e.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

Q <- 0 ; P <- 3
arma.6f <- corARMA( value=rep(c(.6,.2),c(Q,P)), form = ~WK | IID/DRUG, q=Q,p=P, fixed=T)
arma.6f <- Initialize( arma.6f, data=TAB )
arma.6f.mat <- corMatrix( arma.6f )
quartz() ; heatmap.2( arma.6f.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( arma.6f.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )

Q <- 0 ; P <- 2
arma.6f <- corARMA( value=c(.3,.2), form = ~WK | IID/DRUG, q=Q,p=P, fixed=T)
arma.6f <- Initialize( arma.6f, data=TAB )
arma.6f.mat <- corMatrix( arma.6f )
quartz() ; heatmap.2( arma.6f.mat$`B012328-28/1`, main="DRUG==1",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )
quartz() ; heatmap.2( arma.6f.mat$`B012328-28/0`, main="DRUG==0",scale="none",Colv=F,Rowv=F, trace="none",dendrogram="none", col=COLS )




