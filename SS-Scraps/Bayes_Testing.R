TAB.l <- read.table( "Data/Burn/Data/Phenos/Time_Series/20150530_Resp_v_Time.txt",sep="\t",header=T)

AA <- TAB.l[which(TAB.l$IID %in% unique(TAB.l$IID)[1:20]), ]
N.meas <- table( as.character(AA$IID) )
DF <- length(N.meas) / sum( 1/N.meas)
XX <- seq( -10,10,.1)
plot( XX, dinvchisq(XX,DF) )
Sigma.const <- dinvchisq( N.meas, DF )
S.ind <- aggregate( DAS ~ ID, data=AA, var )

S.ind <- numeric( length(N.meas))
for ( i in 1:length(N.meas)) { S.ind[i] <- var( AA$DAS[which(AA$IID==names(N.meas)[i])], na.rm=T) }

Sigma.ind <- S.ind * Sigma.const
Sigma.ind

## Build LME Model
MOD <- lme( DAS ~ )








##################################################
## Recreate Analyses from Zucker, et al. (1997) ##
##################################################

## Input Data
MEAS <- c(3,3,3,3,3,3,3,3,3,3,4,6,3,4,3,4,4,3,4,3,3,3,3)
N.pats <- length(MEAS)
ID.ind <- c(1,22,19,20,2,3,4,10,13,7,9,18,21,12,16,23,17,5,15,14,6,11,8)
ID <- rep( ID.ind,MEAS)
DOS.ind <- c(10,10,10,10,50,10,50,10,10,10,10,25,5,20,20,10,25,10,50,20,10,10,25)
DOS <- rep( DOS.ind,MEAS)
DIFF <- c(-.43,-.64,0,
	.35,-2,-.92,
	.07,-.07,.15,
	.29,-1.42,.14,
	.57,-2.57,-1.07,
	.08,.43,-1.07,
	-.36,-.14,.64,
	-.21,1.28,-1.5,
	0.68,-.25,-.09,
	-.64,1,-.29,
	.05,-.22,.57,.36,
	.64,1.08,-.36,.79,-.64,1.5,
	-.29,1.08,1.5,
	4.29,3.15,0.78,4.49,
	0.64,0,.5,
	1.22,1.07,-.08,.5,
	-.08,.86,1.07,1.15,
	1.85,.78,.54,
	.86,1.43,.65,1.86,
	.86,1.39,1.85,
	1.35,1.36,0.79,
	.47,.93,.86,
	.42,.71,.66 )
TAB.1 <- data.frame( ID, DOS, DIFF )
TAB <- TAB.1[ order(TAB.1$ID),]
MEAS <- MEAS[ order(ID.ind) ]

## Derived Data
PMD <- c(-.05,.25,.09,.28,.3,.27,.28,.37,.3,.34,.29,.46,.51,.67,.4,.53,.56,.63,.74,.85,.88,.66,.57)
PSD <- c(.32,.51,.12,.48,.53,.45,.37,.51,.36,.45,.29,.44,.46,.57,.28,.39,.38,.43,.4,.4,.32,.23,.15)
DER <- data.frame( ID.ind, PMD, PSD )
 # Plot Posterior Distributions for Individuals
XX <- seq( -10,10,.05 )
COLS.list <- c("black","slateblue2","steelblue2","springgreen2","gold1","chocolate2","firebrick2")
COLS <- colorRampPalette(COLS.list)(nrow(DER))
PROBS.der <- numeric( nrow(DER) )
plot( 0,0,type="n",xlim=c(-2,4),ylim=c(0,2),xlab="Posterior Mean Difference",ylab="" )
abline( v=0,lty=2,lwd=2 )
for ( i in 1:nrow(DER) ) {
	points( XX, dnorm( XX,DER[i,"PMD"],DER[i,"PSD"]), type="l",col=COLS[i] )
	PROBS.der[i] <- pnorm( 0, DER[i,"PMD"],DER[i,"PSD"], lower.tail=F )
	text( 3,2-quantile(c(0,2),i/nrow(DER)), label=round(PROBS.der[i],2),col=COLS[i] )
}

## Calculate Summary Stats
MEANS <- aggregate( DIFF ~ ID, data=TAB, mean )
MEANS <- MEANS[ unique( TAB$ID ), ]
round( MEANS, 2)
SDS <- aggregate( DIFF ~ ID, data=TAB, sd )
SDS <- SDS[ unique( TAB$ID ), ]
round( SDS, 2)
RAW <- data.frame( ID=ID.ind, MEANS=MEANS[,"DIFF"], SDS=SDS[,"DIFF"] )
 # Plot Posterior Distributions for Individuals (Using Raw Data)
XX <- seq( -10,10,.05 )
COLS.list <- c("black","slateblue2","steelblue2","springgreen2","gold1","chocolate2","firebrick2")
COLS <- colorRampPalette(COLS.list)(nrow(RAW))
PROBS.raw <- numeric( nrow(RAW) )
plot( 0,0,type="n",xlim=c(-2,4),ylim=c(0,2),xlab="Posterior Mean Difference",ylab="" )
abline( v=0,lty=2,lwd=2 )
for ( i in 1:nrow(RAW) ) {
	points( XX, dnorm( XX,RAW[i,"MEANS"],RAW[i,"SDS"]), type="l",col=COLS[i] )
	PROBS.raw[i] <- pnorm( 0, RAW[i,"MEANS"],RAW[i,"SDS"], lower.tail=F )
	text( 3,2-quantile(c(0,2),i/nrow(RAW)), label=round(PROBS.raw[i],2),col=COLS[i] )
}

## Calculate Priors for Sigma(i)
DF <- N.pats / sum( 1/MEAS)
XX <- seq( 0,50,.2)
plot( XX, dinvchisq(XX,DF), type="l" )
Sigma.const <- dinvchisq( MEAS, DF )
S.ind <- aggregate( DIFF ~ ID, data=TAB, var )
Sigma.ind <- S.ind[,"DIFF"] * Sigma.const
TEMP <- data.frame( MEAS=MEAS, SD=SDS[,"DIFF"], S.i=S.ind[,"DIFF"], Sig.c=Sigma.const, Sig.i=Sigma.ind )
pairs(TEMP)

## LME w/ their Data
MOD <- lme( DIFF ~ 1, data=TAB, random= ~1 | ID )
RAW.means <- MEANS[order(MEANS$ID),"DIFF"]
FIT.means <- ranef(MOD)[,1] + fixef(MOD)
RAW.sds <- SDS[order(MEANS$ID),"DIFF"]
LIM <- range(TAB$DIFF) # range( RAW.means, FIT.means )
plot( 0,0,type="n", xlim=LIM,ylim=LIM,xlab="RAW",ylab="FIT",main="Individual Mean Estimates" )
abline( lm(FIT.means~RAW.means), lty=2 )
abline( v=mean(RAW.means) )
abline( h=fixef(MOD) )
abline( 0,1 )
points( TAB$DIFF, TAB$DIFF, pch=1, col=rep(COLS,MEAS),cex=.7 )
points( RAW.means, FIT.means, cex=2*RAW.sds, col=COLS,pch=20 )



TEMP <- data.frame( TAB, RES=MOD$residuals[,"ID"], RAW.MN=rep(RAW.means,MEAS), FIT.MNS=rep(FIT.means,MEAS) )
TEMP.sd <- aggregate( RES ~ ID, data=TEMP, sd )
TEMP.var <- aggregate( RES ~ ID, data=TEMP, sd )














