## Simulate Clinical Trial by Staggering Start Time, then Assess Drug Efficacy over Time ##
 # Do Longitudinal Analysis looking for Deviation from Baseline after Treatment ##
 # Use .../2013_Longitud_Cancer.pdf (Drescher, et al.) as Reference
## Janssen Data using Multiple-Measures ##
## June 15, 2015 ##
## Kristopher Standish ##

## Rough Game Plan: Longitudinal Analysis as Template ##
 # Calculate DRUG efficacy w/ all patients, Placebo efficacy w/ subset
 # Simulate Start Dates for Patients (uniformly distributed over 100 weeks)
   # Add "Absolute" week as variable (for tracking over time)
 # Starting at week 0, calculate DRUG & PLAC effect size over time
   # Calculate Effect Size for each Individual and Overall Effect
## Specific Game Plan:
 # Simulate Start Dates and Start Moving over Time
 # Calculate Relevant Metrics for each Patient (in trial at absoluate time point)
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
PathToData <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_LongTrial/",sep="" )
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

##############################################################
## FILTER DATA ###############################################
##############################################################

TAB <- TAB.l

## Add ID Column w/ Short Hand
TAB.id.short <- sapply( strsplit( TAB$IID, "-" ),"[",1 )
TAB <- data.frame( ID=TAB.id.short, TAB[,-1] )

## Change ACPA to 1/0
TAB$ACPA <- as.numeric(TAB$ACPA=="Positive")

## Take out Patients who left before getting DRUG
RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
RM.exit <- which( TAB$FID %in% RM.exit.id )
TAB <- TAB[ -RM.exit, ]
SAMPS <- unique(as.character(TAB$FID))
SAMPS.short <- unique(as.character(TAB$ID))
N.SAMPS <- length(SAMPS)

## Take out Patients from Golimumab Arm
RM.gol.id <- as.character( FUL$ID_2[which(FUL$GRP=="G")] )
RM.gol <- which( TAB$FID %in% RM.gol.id )
TAB <- TAB[ -RM.gol, ]

## Make Clinical Table for each Phenotype
RESP_PHENOS <- c("DAS","lCRP","rSJC","rTJC")

## FCT: Create Clin Table for each Response Phenotype
MAKE_CLIN_TAB <- function( RESP_PHENO ) {
	## Pull Response Phenotype of Interest
	TEMP <- TAB[ , c(1:grep("PLAC",colnames(TAB)),which(colnames(TAB)==RESP_PHENO)) ]
	## Remove DAS values that are NA
	RM.na <- which(is.na( TEMP[,RESP_PHENO] ))
	if (length(RM.na) > 0 ) { TEMP <- TEMP[-RM.na, ] }
	## Return Table
	return(TEMP)
}

## Make Clinical Table for each Phenotype
RESP_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
CLIN_TABS <- lapply( RESP_PHENOS, function(x)MAKE_CLIN_TAB(x) )
names(CLIN_TABS) <- RESP_PHENOS
lapply(CLIN_TABS,dim)
TAB <- CLIN_TABS$DAS

# ## Remove NA Values for DAS
# TAB <- TAB[ which(!is.na(TAB$DAS)), ]
# dim(TAB)

##############################################################
## FCT: CALCULATE DRUG EFFECT AT GIVEN TIMEPOINT #############
##############################################################

## Create Function to Assess Drug Efficacy at given Timepoint
BY_WK <- function( TAB.a, wk ) {
	## Remove Measurements after "wk"
	MEAS <- TAB.a[ which(TAB.a$WK.ab<=wk), ]
	n.DRUG <- length(which(MEAS$DRUG==1))
	n.PLAC <- length(which(MEAS$PLAC==1))
	if ( n.DRUG>1 & n.PLAC>1 ) {
		LME <- try( lme( fixed= DAS ~ DRUG+PLAC, data=MEAS, random= ~DRUG+PLAC|FID ), silent=T)
		if ( class(LME)=="try-error" ) {
			print(paste( "P+D: Fail.1 -",wk ))
			LME <- try( lme( fixed= DAS ~ DRUG+PLAC, data=MEAS, random= ~DRUG|FID ), silent=T)
		}
		if ( class(LME)=="try-error" ) {
			print(paste( "P+D: Fail.2 -",wk ))
			LME <- try( lme( fixed= DAS ~ DRUG+PLAC, data=MEAS, random= ~PLAC|FID ), silent=T)
		}
		if ( class(LME)=="try-error" ) {
			print(paste( "P+D: Fail.3 -",wk ))
			LME <- try( lme( fixed= DAS ~ DRUG+PLAC, data=MEAS, random= ~1|FID ), silent=T)
		}
		if ( class(LME)=="try-error" ) { print(paste( "P+D: Fail.All -",wk )) }
	}else{
		if ( n.DRUG>1 ) {
			LME <- try( lme( fixed= DAS ~ DRUG, data=MEAS, random= ~DRUG|FID ), silent=T)
			if ( class(LME)=="try-error" ) {
				print(paste( "D: Fail -",wk ))
				LME <- lme( fixed= DAS ~ DRUG+PLAC, data=MEAS, random= ~1|FID )	
			}
		}
		if ( n.PLAC>1 ) {
			LME <- try( lme( fixed= DAS ~ PLAC, data=MEAS, random= ~PLAC|FID ), silent=T)
			if ( class(LME)=="try-error" ) {
				print(paste( "P: Fail -",wk ))
				LME <- lme( fixed= DAS ~ PLAC, data=MEAS, random= ~1|FID )	
			}
		}
	}

	## Compile Output Stats
	if ( (n.DRUG>1 | n.PLAC>1) & class(LME)!="try-error" ) {
		COEF.0 <- summary(LME)$tTable
		COEF.i <- coef(LME)
		VAR.i <- LME$sigma
		VAR.0 <- VarCorr(LME)
		COMPILE <- list( COEF.0, COEF.i, VAR.i, VAR.0, LME )
		names(COMPILE) <- c( "COEF.0","COEF.i","VAR.i","VAR.0","LME_Model" )
		return(COMPILE)		
	}else{
		return("No Model")
	}
}

##############################################################
## SIMULATE CLINICAL TRIAL ENROLLMENT ########################
##############################################################

## Get Number of Patients
Samps <- as.character(unique( TAB$FID ))
N.samps <- length(Samps)

## Simulate Start Dates for Patients
Enroll.countby <- 2
Enroll.start <- 0
Enroll.end <- 100
Enroll.samps <- sample( seq(Enroll.start, Enroll.end, Enroll.countby), N.samps, replace=T )
names(Enroll.samps) <- Samps

## Add Absolute Week to Data Table
TAB.a <- data.frame( TAB, WK.ab=0 )
for ( s in 1:N.samps ) {
	samp <- Samps[s]
	TAB.a[ which(TAB.a$FID==samp), "WK.ab" ] <- TAB.a[ which(TAB.a$FID==samp), "WK" ] + Enroll.samps[s]
} ; hist( TAB.a$WK.ab, breaks=seq(0,100+Enroll.end,Enroll.countby),col="dodgerblue2" )
 # Get Range of Weeks
WKS <- min(TAB.a$WK.ab):max(TAB.a$WK.ab)
WKS <- seq( min(TAB.a$WK.ab), max(TAB.a$WK.ab), Enroll.countby )

##############################################################
## COMPILE DATA OVER COURSE OF TRIAL #########################
##############################################################

## Compile data over Trial
OUT <- list()
start_time <- proc.time()
for ( w in 1:length(WKS) ) {
# for ( w in seq(10,50,5) ) {
# for ( w in 1:100 ) {
	wk <- WKS[w]
	tag <- paste("Wk",wk,sep="_")
	OUT[[tag]] <- BY_WK( TAB.a, wk )
	if ( wk%%10==0 ) { print(paste( "Done with week",wk,". (",w,"of",length(WKS),") -",round(proc.time()-start_time,3)[3] )) }
}
if ( exists("OUT.save") ) { OUT.save[[length(OUT.save)+1]] <- OUT }else{ OUT.save <- list() ; OUT.save[[1]] <- OUT }

##############################################################
## PARSE TRIAL UPDATES OVER TIME #############################
##############################################################

## Which Models failed Entirely
FAIL.which <- which(summary(OUT)[,"Length"]==1) # which(OUT=="No Model")
if ( length(FAIL.which)>0 ) {
	OUT.2 <- OUT[-FAIL.which]
}else{ OUT.2 <- OUT }

## Which Weeks yielded actual results?
WKS.plot <- as.numeric( gsub("Wk_","",names(OUT.2)) )

##########################
## POPULATION STATS ######

## Which Weeks have Valid PLAC/DRUG Estimates
 # Fixed
DRUG.which <- which( unlist(lapply( OUT.2, function(x) "DRUG"%in%rownames(x$COEF.0) )) )
PLAC.which <- which( unlist(lapply( OUT.2, function(x) "PLAC"%in%rownames(x$COEF.0) )) )
 # Random
DRUG.r.which <- which( unlist(lapply( OUT.2, function(x) "DRUG"%in%rownames(x$VAR.0) )) )
PLAC.r.which <- which( unlist(lapply( OUT.2, function(x) "PLAC"%in%rownames(x$VAR.0) )) )

## Compile Population Estimates for each Week
ARR.colnames <- c( "DRUG.f","PLAC.f","DRUG.r","PLAC.r","DRUG.f.se","PLAC.f.se" )
ARR <- array( ,c(length(OUT.2),length(ARR.colnames)) ) ; colnames(ARR) <- ARR.colnames ; rownames(ARR) <- names(OUT.2)
ARR[ DRUG.which, "DRUG.f" ] <- as.numeric( unlist(lapply( OUT.2[DRUG.which], function(x) x$COEF.0["DRUG","Value"] )) )
ARR[ PLAC.which, "PLAC.f" ] <- as.numeric( unlist(lapply( OUT.2[PLAC.which], function(x) x$COEF.0["PLAC","Value"] )) )
ARR[ DRUG.r.which, "DRUG.r" ] <- as.numeric( unlist(lapply( OUT.2[DRUG.r.which], function(x) x$VAR.0["DRUG","Variance"] )) )
ARR[ PLAC.r.which, "PLAC.r" ] <- as.numeric( unlist(lapply( OUT.2[PLAC.r.which], function(x) x$VAR.0["PLAC","Variance"] )) )
ARR[ DRUG.which, "DRUG.f.se" ] <- as.numeric( unlist(lapply( OUT.2[DRUG.which], function(x) x$COEF.0["DRUG","Std.Error"] )) )
ARR[ PLAC.which, "PLAC.f.se" ] <- as.numeric( unlist(lapply( OUT.2[PLAC.which], function(x) x$COEF.0["PLAC","Std.Error"] )) )

## Calculate Posterior Probability of Drug Efficacy > Placebo Efficacy by X
 # Use very naive approach: Normal distribution around mean
EFF_SIZES <- seq( -2,1,.25 )
POST.pop.drug <- POST.pop.plac <- list()
for ( e in 1:length(EFF_SIZES) ) {
	eff_size <- EFF_SIZES[e]
	tag <- paste("Eff",eff_size,sep="_")
	POST.pop.drug[[tag]] <- pnorm( eff_size, ARR[DRUG.which,"DRUG.f"],ARR[DRUG.which,"DRUG.f.se"], lower.tail=T )
	POST.pop.plac[[tag]] <- pnorm( eff_size, ARR[PLAC.which,"PLAC.f"],ARR[PLAC.which,"PLAC.f.se"], lower.tail=T )
}

##########################
## INDIVIDUAL STATS ######

## Create Array w/ Samples by Week
SAMP.i.which <- lapply( OUT.2, function(x) rownames(x$COEF.i) )

## Compile Individual Estimates for each Week
 # Set up Array
IND.colnames <- c("DRUG","PLAC")
IND.dimnames.3 <- names(OUT.2)
IND <- array( ,c(N.samps,length(IND.colnames),length(IND.dimnames.3)) )
rownames(IND) <- Samps ; colnames(IND) <- IND.colnames ; dimnames(IND)[[3]] <- IND.dimnames.3
 # Fill Array w/ Numbers
for ( w in 1:length(OUT.2) ) {
	tag <- names(OUT.2)[w]
	if ( w %in% DRUG.which) { IND[SAMP.i.which[[tag]],"DRUG",tag] <- OUT.2[[tag]]$COEF.i[SAMP.i.which[[tag]],"DRUG"] }
	if ( w %in% PLAC.which) { IND[SAMP.i.which[[tag]],"PLAC",tag] <- OUT.2[[tag]]$COEF.i[SAMP.i.which[[tag]],"PLAC"] }
}

## Model Final Coefficient Estimates as Function of Previous Coefficient Estimates (over time)
COR.mods.drug <- COR.mods.plac <- list()
CORS.drug <- numeric( length(WKS.plot) ) ; names(CORS.drug) <- dimnames(IND)[[3]]
CORS.plac <- CORS.drug
for ( w in 1:length(WKS.plot) ) {
	wk <- WKS.plot[w]
	tag <- dimnames(IND)[[3]][w]
	temp_drug <- IND[,"DRUG",tag]
	if ( length(which(!is.na(temp_drug)))>2 ) {
		COR.mods.drug[[tag]] <- lm( IND[,"DRUG",dim(IND)[3]] ~ temp_drug )
	}else{ COR.mods.drug[[tag]] <- "No Model" }
	temp_plac <- IND[,"PLAC",tag]
	if ( length(which(!is.na(temp_plac)))>2 ) {
		COR.mods.plac[[tag]] <- lm( IND[,"PLAC",dim(IND)[3]] ~ temp_plac )
	}else{ COR.mods.plac[[tag]] <- "No Model" }
	CORS.drug[w] <- cor( IND[,"DRUG",dim(IND)[3]], IND[,"DRUG",tag], use="pairwise.complete.obs" )
	CORS.plac[w] <- cor( IND[,"PLAC",dim(IND)[3]], IND[,"PLAC",tag], use="pairwise.complete.obs" )
}

## Calculate Posterior Probability of Drug/Plac having Effect greater than E at each week
 # Use Bayesian Approach
   # Flat Prior at week 0 (alpha=0,beta=0)
   # Update Prior each week of trial
   # Calculate Prob of Delta Greater than E for each week
   # Do this for DRUG and PLAC
EFF_SIZES <- seq( -2,1,.5 )
BAYES.pop.drug <- BAYES.pop.plac <- list()
for ( e in 1:length(EFF_SIZES) ) {
	eff_size <- EFF_SIZES[e]
	tag <- paste("Eff",eff_size,sep="_")
	xvals <- seq(0,1,.001)
	o.pos.drug <- apply( IND, 3, function(x) length(which(x[,"DRUG"]<eff_size)) )
	o.neg.drug <- apply( IND, 3, function(x) length(which(x[,"DRUG"]>=eff_size)) )
	BAYES.pop.drug[[tag]] <- apply( data.frame(cumsum(o.pos.drug)+1,cumsum(o.neg.drug)+1), 1, function(x) dbeta(xvals,x[1],x[2]) )
	o.pos.plac <- apply( IND, 3, function(x) length(which(x[,"PLAC"]<eff_size)) )
	o.neg.plac <- apply( IND, 3, function(x) length(which(x[,"PLAC"]>=eff_size)) )
	BAYES.pop.plac[[tag]] <- apply( data.frame(cumsum(o.pos.plac)+1,cumsum(o.neg.plac)+1), 1, function(x) dbeta(xvals,x[1],x[2]) )
	rownames(BAYES.pop.plac[[tag]]) <- rownames(BAYES.pop.drug[[tag]]) <- xvals
}

##############################################################
## PLOT TRIAL UPDATES OVER TIME ##############################
##############################################################

## Plot Trial Enrollment Distribution over Time
png( paste(PathToPlot,"0-Observation_Hist.png",sep=""), height=1200,width=1600,pointsize=24 )
par(mfrow=c(2,1))
 # Individual Patient Enrollment over Time
XLIM <- range(WKS)
YLIM <- c(0,N.samps)
COLS.list <- c("slateblue1","steelblue1","springgreen1","gold1","chocolate1","firebrick1")
COLS <- colorRampPalette(COLS.list)(N.samps)
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Weeks",ylab="Patient",main="Patient Enrollment over Time")
for ( s in 1:N.samps ) {
	samp.plot <- names(Enroll.samps)[order(Enroll.samps)][s]
	which_rows <- which(TAB.a$FID==samp.plot)
	points( TAB.a$WK.ab[which_rows], rep(s,length(which_rows)), col=COLS[s],type="o",pch=20,lwd=2 )
}
 # Time of Measurements
# barplot( table(TAB.a$WK.ab), las=2,col="dodgerblue2" )
hist( TAB.a$WK.ab, breaks=seq(0,100+Enroll.end,Enroll.countby),col="dodgerblue2" )
dev.off()

##########################
## POPULATION STATS ######

## Plot Polulation-Level Fixed Coefficients over Time (w/ Error)
 # Parameters
XLIM <- c( 0, max(WKS.plot) )
YLIM <- range( ARR[,1]+ARR[,3], ARR[,2]+ARR[,4], ARR[,1]-ARR[,3], ARR[,2]-ARR[,4] ,na.rm=T) ; YLIM <- c( max(-4,YLIM[1]),min(2,YLIM[2]) )
COLS.which <- c(1,3)
COLS <- COLS.list[COLS.which]
COLS.pol <- sapply( COLS, function(x) colorRampPalette(c("white",x))(3)[2] )
 # Create Plot
png( paste(PathToPlot,"1-Population_Plot.png",sep=""), height=1200,width=1600,pointsize=24 )
par(mfrow=c(2,1))
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Weeks",ylab="Effect Size Estimate",main="Effect Size Coefficient over Time" )
abline( h=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
abline( v=seq(0,XLIM[2]+10,10),lty=3,col="grey50",lwd=1 )
 # Plot Data
polygon( c(WKS.plot[PLAC.r.which],rev(WKS.plot[PLAC.r.which])), c(ARR[PLAC.r.which,"PLAC.f"]+ARR[PLAC.r.which,"PLAC.r"],ARR[rev(PLAC.r.which),"PLAC.f"]-ARR[rev(PLAC.r.which),"PLAC.r"]), col=COLS.pol[2],border=COLS[2],density=50,angle=45 )
polygon( c(WKS.plot[PLAC.which],rev(WKS.plot[PLAC.which])), c(ARR[PLAC.which,"PLAC.f"]+ARR[PLAC.which,"PLAC.f.se"],ARR[rev(PLAC.which),"PLAC.f"]-ARR[rev(PLAC.which),"PLAC.f.se"]), col=COLS.pol[2],border=COLS[2],density=80,angle=45 )
points( WKS.plot, ARR[,"PLAC.f"], type="l",col=COLS[2], lwd=3 )
polygon( c(WKS.plot[DRUG.r.which],rev(WKS.plot[DRUG.r.which])), c(ARR[DRUG.r.which,"DRUG.f"]+ARR[DRUG.r.which,"DRUG.r"],ARR[rev(DRUG.r.which),"DRUG.f"]-ARR[rev(DRUG.r.which),"DRUG.r"]), col=COLS.pol[1],border=COLS[1],density=50,angle=-45 )
polygon( c(WKS.plot[DRUG.which],rev(WKS.plot[DRUG.which])), c(ARR[DRUG.which,"DRUG.f"]+ARR[DRUG.which,"DRUG.f.se"],ARR[rev(DRUG.which),"DRUG.f"]-ARR[rev(DRUG.which),"DRUG.f.se"]), col=COLS.pol[1],border=COLS[1],density=80,angle=-45 )
points( WKS.plot, ARR[,"DRUG.f"], type="l",col=COLS[1], lwd=3 )
 # Legend
legend( "bottomleft", fill=COLS, legend=c("DRUG","PLAC") )
# dev.off()

## Plot Polulation-Level Posterior Probabilities (w/ Error)
 # Parameters
XLIM <- c( 0, max(WKS.plot) )
YLIM <- c(0,1)
COLS <- COLS.list[COLS.which]
COLS.ramp <- sapply( COLS, function(x) colorRampPalette(c("white",x,"black"))(length(EFF_SIZES)+4)[1:length(EFF_SIZES)+2] )
COLS.ramp.leg <- colorRampPalette(c("white","grey50","black"))(length(EFF_SIZES)+4)[1:length(EFF_SIZES)+2]
 # Create Plot
# png( paste(PathToPlot,"1-Population_Posterior.png",sep=""), height=800,width=1600,pointsize=24 )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Weeks",ylab="Posterior Probability Estimate",main="Posterior Probability of Effect > |k|" )
abline( h=seq(0,1,.2),lty=3,col="grey50",lwd=1 )
abline( v=seq(0,XLIM[2]+10,10),lty=3,col="grey50",lwd=1 )
 # Plot Data
for ( e in 1:length(EFF_SIZES) ) {
	eff_size <- EFF_SIZES[e]
	tag <- paste("Eff",eff_size,sep="_")
	points( WKS.plot[DRUG.which], POST.pop.drug[[tag]], type="l",col=COLS.ramp[e,1], lwd=3,lty=2 )
	text( max(WKS.plot), POST.pop.drug[[tag]][length(POST.pop.drug[[tag]])], label=eff_size,col=COLS.ramp[e,1],pos=4 )
	points( WKS.plot[PLAC.which], POST.pop.plac[[tag]], type="l",col=COLS.ramp[e,2], lwd=3,lty=2 )
	text( max(WKS.plot), POST.pop.plac[[tag]][length(POST.pop.plac[[tag]])], label=eff_size,col=COLS.ramp[e,2],pos=4 )
}
 # Legend
legend( quantile(XLIM,.9),quantile(YLIM,.8), fill=COLS, legend=c("DRUG","PLAC") )
legend( quantile(XLIM,.8),quantile(YLIM,.8), col=COLS.ramp.leg, legend=EFF_SIZES,title="Effect Sizes",lwd=3,lty=2 )
dev.off()


##########################
## INDIVIDUAL STATS ######

## Plot Individual Coefficients over Time
 # Parameters
File_Tag <- sample(1:100,1)
N.samps.plot <- 5
Samps.id.plot <- sample( Samps, N.samps.plot, replace=F )
XLIM <- c( 0, max(WKS.plot) ) + c(0,20)
YLIM <- range(IND,na.rm=T)
YLIM <- range(ARR[,1:2],IND[Samps.id.plot,,],na.rm=T)
COLS.4.list <- gsub("1","4",COLS.list)
COLS <- colorRampPalette(COLS.list)(N.samps.plot)
COLS.4 <- colorRampPalette(COLS.4.list)(N.samps.plot)
 # Create Plot
png( paste(PathToPlot,"2-Indiv_Plot_",File_Tag,".png",sep=""), height=800,width=1600,pointsize=24 )
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Weeks",ylab="Effect Size Estimate",main="Effect Size Coefficient over Time" )
abline( h=seq(-10,10,.5),lty=3,col="grey50",lwd=1 )
abline( v=seq(0,XLIM[2]+10,10),lty=3,col="grey50",lwd=1 )
 # Population-Level Estimates
points( WKS.plot, ARR[,"DRUG.f"], type="l",col="black", lwd=5,lty=1 )
points( WKS.plot, ARR[,"PLAC.f"], type="l",col="black", lwd=5,lty=2 )
 # Individual Data
for ( s in 1:N.samps.plot ) {
	samp.plot <- Samps.id.plot[s]
	TEMP <- TAB.a[which(TAB.a$FID==samp.plot),]
	 # Lines
	points( WKS.plot, IND[samp.plot,"DRUG",], type="l",col=COLS[s],lty=1,lwd=3 )
	points( WKS.plot, IND[samp.plot,"PLAC",], type="l",col=COLS[s],lty=2,lwd=3 )
	 # Points (Measured)
	TEMP.wks <- TEMP$WK.ab[ which(TEMP$WK.ab %in% WKS.plot) ]
	points( TEMP.wks, IND[samp.plot,"DRUG",paste("Wk",TEMP.wks,sep="_")],col=COLS.4[s],pch=c(1,16)[factor(TEMP$DRUG)] )
	points( TEMP.wks, IND[samp.plot,"PLAC",paste("Wk",TEMP.wks,sep="_")],col=COLS.4[s],pch=c(1,16)[factor(TEMP$DRUG)] )
	 # Sample Names & Range of Trial Involvement
	text( 200, IND[samp.plot,"DRUG",dim(IND)[3]], label=samp.plot, col=COLS[s], pos=4 )
	# abline( v=TEMP.wks[c(1,length(TEMP.wks))], lty=4,col=COLS[s],lwd=2 )
	arrows( TEMP.wks[c(1,length(TEMP.wks))], quantile(YLIM,.05), TEMP.wks[c(1,length(TEMP.wks))], YLIM[1], lwd=3,col=COLS[s],length=.2 )
}
legend( "topright",lty=c(2,1),pch=c(1,16),lwd=3,legend=c("PLAC","DRUG"), bg="white" )
 # Show which weeks had Random Effects Estimates
# points( WKS.plot[DRUG.r.which],rep(YLIM[2],length(DRUG.r.which)),pch=16 )
# points( WKS.plot[PLAC.r.which],rep(YLIM[2]+.05,length(PLAC.r.which)),pch=1 )
dev.off()

## Plot Correlation b/n Coefficient Estimates over Time (w/ Final Estimates)
drug.which <- which( COR.mods.drug!="No Model" )
drug.which <- which(unlist(lapply( COR.mods.drug[drug.which], function(x) "temp_drug" %in% rownames(summary(x)$coefficients) ))) # 
COR.p.drug <- unlist(lapply( COR.mods.drug[names(drug.which)], function(x) summary(x)$coefficients["temp_drug",'Pr(>|t|)'] ))
COR.p.drug[which(COR.p.drug==0)] <- min(COR.p.drug[which(COR.p.drug>0)])
plac.which <- which( COR.mods.plac!="No Model" )
plac.which <- which(unlist(lapply( COR.mods.plac[plac.which], function(x) "temp_plac" %in% rownames(summary(x)$coefficients) ))) # 
COR.p.plac <- unlist(lapply( COR.mods.plac[names(plac.which)], function(x) summary(x)$coefficients["temp_plac",'Pr(>|t|)'] ))
COR.p.plac[which(COR.p.plac==0)] <- min(COR.p.plac[which(COR.p.plac>0)])
 # Write File
png( paste(PathToPlot,"3-Coef_Estim_Corr.png",sep=""), height=1200,width=1600,pointsize=24 )
par(mfrow=c(2,1))
 # Parameters
XLIM <- c(0, max(WKS.plot) )
YLIM <- range( CORS.drug,CORS.plac,na.rm=T)
COLS <- COLS.list[COLS.which]
 # Create Plot
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Weeks",ylab="Correlation Coefficient",main="Correlation of Coefficient Estimates\nOver Time w/ Final Estimates" )
abline( h=seq(-1,1,.25),lty=3,col="grey50",lwd=1 )
abline( h=c(-1,0,1),lty=1,col="black",lwd=1 )
abline( v=seq(0,XLIM[2]+10,10),lty=3,col="grey50",lwd=1 )
 # Plot Results
points( as.numeric(gsub("Wk_","",names(CORS.drug))),CORS.drug, type="l", col=COLS[1],lwd=3 )
points( as.numeric(gsub("Wk_","",names(CORS.plac))),CORS.plac, type="l", col=COLS[2],lwd=3 )
 # Parameters
YLIM <- c(0, -log10(min(COR.p.drug,COR.p.plac)) )
 # Create Plot
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Weeks",ylab="-log10(p)",main="Correlation of Coefficient Estimates\nOver Time w/ Final Estimates" )
abline( h=seq(0,YLIM[2]+50,50),lty=3,col="grey50",lwd=1 )
abline( v=seq(0,XLIM[2]+10,10),lty=3,col="grey50",lwd=1 )
 # Plot Results
points( as.numeric(gsub("Wk_","",names(COR.p.drug))),-log10(COR.p.drug), type="l", col=COLS[1],lwd=3 )
points( as.numeric(gsub("Wk_","",names(COR.p.plac))),-log10(COR.p.plac), type="l", col=COLS[2],lwd=3 )
dev.off()

## Plot Posterior Probability of an Individual having Effect Size greater than E
COLS <- COLS.list[COLS.which]
for ( e in 1:length(EFF_SIZES) ) {
	eff_size <- EFF_SIZES[e]
	tag <- paste("Eff",eff_size,sep="_")
	png( paste(PathToPlot,"4-Bayes_Posterior_D_",tag,".png",sep=""), height=1200,width=1600,pointsize=24 )
	heatmap.2( sqrt(BAYES.pop.drug[[tag]][nrow(BAYES.pop.drug[[tag]]):1,]), col=colorRampPalette(c("black",COLS[1],"white"))(100), scale="none",Rowv=F,Colv=F,dendrogram="none",trace="none" )
	dev.off()
	png( paste(PathToPlot,"4-Bayes_Posterior_P_",tag,".png",sep=""), height=1200,width=1600,pointsize=24 )
	heatmap.2( sqrt(BAYES.pop.plac[[tag]][nrow(BAYES.pop.plac[[tag]]):1,]), col=colorRampPalette(c("black",COLS[2],"white"))(100), scale="none",Rowv=F,Colv=F,dendrogram="none",trace="none" )
	dev.off()
	# plot( 0,0,type="n",xlim=c(0,1),ylim=c(0,100) )
	# for ( col in seq(1,100,5) ) { points( xvals, BAYES.pop.drug[[3]][,col],type="l",col=COLS[col] ) }
	# for ( col in seq(1,100,5) ) { points( xvals, BAYES.pop.plac[[3]][,col],type="l",col=COLS[col] ) }
}

BAYES.pop.both <- list()
for ( e in 1:length(EFF_SIZES) ) {
	eff_size <- EFF_SIZES[e]
	tag <- paste("Eff",eff_size,sep="_")
	BAYES.pop.both[[tag]] <- sqrt(BAYES.pop.drug[[tag]]) - sqrt(BAYES.pop.plac[[tag]])
	TEST.plot <- BAYES.pop.both[[tag]][nrow(BAYES.pop.both[[tag]]):1,]
	png( paste(PathToPlot,"4-Bayes_Posterior_Both_",tag,".png",sep=""), height=1200,width=1600,pointsize=24 )
	heatmap.2( TEST.plot, col=colorRampPalette(c("white",COLS[2],"black",COLS[1],"white"))(200), scale="none",Rowv=F,Colv=F,dendrogram="none",trace="none",breaks=c( seq(min(TEST.plot),0,length.out=101)[1:100],0,seq(0,max(TEST.plot),length.out=101)[1:100] ) )
	dev.off()
}

##############################################################
## END OF DOC ################################################
##############################################################

























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

