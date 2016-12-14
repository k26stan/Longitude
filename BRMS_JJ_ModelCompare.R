## Compare Time-Series Results ##
## Compare "Mean", Individual, and LMM (w/ & w/o Cor) Models
## Using "brms" Package
## Janssen Data using Repeated-Measures ##
## July 01, 2016 ##
## Kristopher Standish ##

###################################################
## Brainstorming Ideas ############################

## Probe potential advantages of this approach
 # How accurate is this approach at forecasting?
   # In clinical trial setting, is there an advantage to this type of approach?
     # Can we determine the efficacy of a drug earlier?
   # Given a set of new patients w/ a couple preliminary measurements, how well can we predict their future response profile?
     # e.g., given 1 baseline measurements & 1 treatment measurement, how well does the model predict subsequent measurements?
 # How accurate is this approach at imputing missing measurements?
   # If we downsample a subset of patients
 # Does this approach give us a better idea of what a responder looks like?
   # We get credible intervals on DRUG effect sizes
   # 

# Rscript <.R> Prior 50 ValPrior_50 XXX
#############################################################
## GLOBALS ##################################################
#############################################################

# ## Parse Command Line
# LINE <- commandArgs(trailingOnly = TRUE)
# # LINE <- c("Null",421,"Null_All","m1,m2,m3,m4,m5,m6,m7,m8,m9,m10")
# # LINE <- c("HLA",421,"HLA_TEST_m1","XXXX")
# Goal <- LINE[1]
# N.Samps <- LINE[2]
# Dir.Tag <- LINE[3]
# Mod.Names.1 <- LINE[4]
# Mod.Names <- unlist(strsplit(Mod.Names.1,",")[[1]])
#  # Print Inputs
# print(paste("Goal:",Goal))
# print(N.Samps)
# print(Dir.Tag)
# print(Mod.Names)

## COLORS: Set Color Scheme
 # FCT: Blend 2 Colors
BLEND <- function( colors ) { colorRampPalette(colors)(3)[2] }
SHOW_COLORS <- function( colors ) { barplot( 1:length(colors),col=colors,names.arg=1:length(colors)) }
 # Set Color Palette
COLS.list <- c("firebrick2","chocolate2","gold1","chartreuse2","cadetblue2","dodgerblue2","slateblue2","magenta2")
COLS.list.2 <- c("aquamarine3","deeppink2",BLEND(c("gold1","chartreuse2")),BLEND(c("steelblue3","slateblue3")),"tomato2","deepskyblue1",BLEND(c("sienna2","goldenrod1")) )
COLS.list.ord <- COLS.list.2[c(5,7,3,1,6,4,2)]
 # Model Colors
COLS.mods <- COLS.list.ord[c(4,2,7,5)]
names(COLS.mods) <- c("Null_Tn","Null_Tnc","Ind","LeadUp")
 # Trial Arm Colors
COLS.arm <- COLS.list.ord[c(1,3,6)]
names(COLS.arm) <- c("G","P","PE")

##############################################################
## LOAD DATA #################################################
##############################################################

## Load Packages
library(nlme)
library(gplots)
library(brms)
library(vioplot)
library(xtable)

## Mac or TSCC?
Mac_Root <- "/Users/kstandis/Data/"
TSCC_Root <- "/projects/janssen/Mac_Data/"
if ( file.exists(Mac_Root) ) {
	Root <- Mac_Root
	PathToNullMod.T <- paste(Root,"Janssen/Data/Longitude/20160702_Null_Trt/",sep="")
	PathToNullMod.Tc <- paste(Root,"Janssen/Data/Longitude/20160702_Null_TrtCor/",sep="")
	PathToNullMod.Tn <- paste(Root,"Janssen/Data/Longitude/20160719_Null_Trt_NoInt/",sep="")
	PathToNullMod.Tnc <- paste(Root,"Janssen/Data/Longitude/20160719_Null_TrtCor_NoInt/",sep="")
	PathToDownMod.1 <- paste(Root,"Janssen/Data/Longitude/20160614_DownSample/",sep="")
	PathToIndMod.1 <- paste(Root,"Janssen/Data/Longitude/20160630_Indiv_Bayes/Rdata.Model.Summarizeind4.Rdata",sep="")
	PathToIndMod.1 <- paste(Root,"Janssen/Data/Longitude/20160719_Indiv_Bayes/Rdata.Model.Summarize.ind3.Rdata",sep="")
	PathToLeadUp <- paste(Root,"Janssen/Data/Longitude/20160706_LeadUp_m91011/",sep="")
	PathToLeadUp <- paste(Root,"Janssen/Data/Longitude/20160726_LeadUp_m91011/",sep="")
	# SNP Model
	PathToSNPMod.1 <- paste(Root,"Janssen/Data/Longitude/20160813_Cand_SNP/",sep="")	
	# DWAI
	PathToNullMod.Tcdi <- paste(Root,"Janssen/Data/Longitude/20160702_Null_TrtCorDrI/",sep="")
	PathToNullMod.Tcid <- paste(Root,"Janssen/Data/Longitude/20160702_Null_TrtCorIDr/",sep="")
}else{
	Root <- TSCC_Root
	PathToNullMod.T <- paste(Root,"Janssen/Plots_Mac/20160702_Null_Trt/",sep="")
	PathToNullMod.Tc <- paste(Root,"Janssen/Plots_Mac/20160702_Null_TrtCor/",sep="")
	PathToNullMod.Tn <- paste(Root,"Janssen/Plots_Mac/20160719_Null_Trt_NoInt/",sep="")
	PathToNullMod.Tnc <- paste(Root,"Janssen/Plots_Mac/20160719_Null_TrtCor_NoInt/",sep="")
	PathToDownMod.1 <- paste(Root,"Janssen/Plots_Mac/20160614_DownSample/",sep="")
	PathToIndMod.1 <- paste(Root,"Janssen/Plots_Mac/20160630_Indiv_Bayes/Rdata.Model.Summarizeind4.Rdata",sep="")
	PathToIndMod.1 <- paste(Root,"Janssen/Plots_Mac/20160719_BayesLMM/Rdata.Model.Summarize.ind3.Rdata",sep="")
	PathToLeadUp <- paste(Root,"Janssen/Plots_Mac/20160706_LeadUp_m91011/",sep="")
	PathToLeadUp <- paste(Root,"Janssen/Plots_Mac/20160726_LeadUp_m91011/",sep="")
	# SNP Model
	PathToSNPMod.1 <- paste(Root,"Janssen/Plots_Mac/20160813_Cand_SNP/",sep="")	
	# DWAI
	PathToNullMod.Tcdi <- paste(Root,"Janssen/Plots_Mac/20160702_Null_TrtCorDrI/",sep="")
	PathToNullMod.Tcid <- paste(Root,"Janssen/Plots_Mac/20160702_Null_TrtCorIDr/",sep="")
}

## Set Date
DATE <- gsub("-","",Sys.Date())

## New Mac Paths
PathToTypes <- paste(Root,"Janssen/Data/HLA/SOAP_HLA_Types/20151211_HLA_Types.Rdata",sep="")
PathToAA <- paste(Root,"Janssen/Data/HLA/Amino_Acids/20160126_HLA_AA.Rdata",sep="")
PathToRawFiles <- paste(Root,"Janssen/Data/Pheno/Raw_Files/",sep="")
PathToAssoc <- paste(Root,"Janssen/Data/HLA/Association/20160614_HLA_Assoc_",sep="")
PathToData <- paste(Root,"Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt",sep="")
PathToFT <- paste(Root,"Janssen/Data/Pheno/Derived/20151015_Full_Table.txt",sep="")
PathToDer <- paste(Root,"Janssen/Data/Pheno/Derived/20150619_Derived_Pheno_Table.txt",sep="")
PathToPlot <- paste(Root,"Janssen/Plots_Mac/",DATE,"_Long_Compare/",sep="")
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Write Update to File
write(paste(date(),"- Data Loaded ######"), paste(PathToPlot,"Update.txt",sep=""),append=F)
##############################################################
## FILTER DATA ###############################################
##############################################################
write(paste(date(),"- Filtering Data ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

X <- 1:10

#########################################
## CLINICAL DATA ########################

## Add Column for "TRT" (Treatment: All but Week 0)
TRT <- rep(1,nrow(TAB.l))
TRT[which(TAB.l$WK==0)] <- 0
TAB <- data.frame( TAB.l, TRT )

## Add ID Column w/ Short Hand
TAB.id.short <- sapply( strsplit( TAB$IID, "-" ),"[",1 )
TAB <- data.frame( ID=TAB.id.short, TAB[,-1] )

## Change ACPA to 1/0
TAB$ACPA <- as.numeric(TAB$ACPA=="Positive")
TAB$RF <- as.numeric(TAB$RF=="Positive")
TAB$RF_ACPA <- as.numeric(TAB$RF_ACPA=="Positive")

# ## Add PC's onto Table for later use
# TAB <- merge( TAB, FUL[,c("ID_2",paste("PC",1:3,sep=""))], by.x="FID",by.y="ID_2")

## Take out Patients who left before getting DRUG
RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
RM.exit <- which( TAB$FID %in% RM.exit.id )
TAB <- TAB[ -RM.exit, ]
SAMPS <- unique(as.character(TAB$FID))
SAMPS.short <- unique(as.character(TAB$ID))
N.SAMPS <- length(SAMPS)
SAMPS.list <- list()
SAMPS.list$g <- as.character(FUL$ID_2[ FUL$GRP=="G" & FUL$ID_2%in%SAMPS ])
SAMPS.list$p <- as.character(FUL$ID_2[ FUL$GRP!="G" & FUL$ID_2%in%SAMPS ])
SAMPS.list$ne <- as.character(FUL$ID_2[ FUL$GRP=="P" & FUL$ID_2%in%SAMPS ])
SAMPS.list$ee <- as.character(FUL$ID_2[ FUL$GRP=="PE" & FUL$ID_2%in%SAMPS ])
SAMPS.list$sg <- as.character(FUL$ID[ FUL$GRP=="G" & FUL$ID_2%in%SAMPS ])
SAMPS.list$sp <- as.character(FUL$ID[ FUL$GRP!="G" & FUL$ID_2%in%SAMPS ])
SAMPS.list$sne <- as.character(FUL$ID[ FUL$GRP=="P" & FUL$ID_2%in%SAMPS ])
SAMPS.list$see <- as.character(FUL$ID[ FUL$GRP=="PE" & FUL$ID_2%in%SAMPS ])

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

#############################################################
## SET UP MODELS ############################################
#############################################################
write(paste(date(),"- Setting up Stan Parameters ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

####################################
## Get Set Up ######################

## Package Parameters
 # From Tutorial: http://www.r-bloggers.com/r-users-will-now-inevitably-become-bayesians/
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

## Model Parameters
Mod.iter <- 3000
Mod.warm <- 600
Mod.chain <- 3
Mod.fam <- "gaussian"

####################################
## Get Models Set Up ###############
JJ.priors <- JJ.form <- JJ.cor <- list()
if ( !exists("JJ") ) { JJ <- JJ.time <- JJ.samps <- list() }

####################################
## Set Priors for Models ###########

## Set Priors for each Variable
JJ.priors.list <- list()
# JJ.priors.list$Intercept <- set_prior("normal(5,1)", class="b",coef="Intercept")
JJ.priors.list$Intercept <- set_prior("normal(5,1)", class="Intercept")
JJ.priors.list$DRUG <- set_prior("normal(-1,1)", class="b",coef="DRUG")
JJ.priors.list$WK <- set_prior("normal(0,.1)", class="b",coef="WK")
JJ.priors.list$PLAC <- set_prior("normal(0,1)", class="b",coef="PLAC")
JJ.priors.list$TRT <- set_prior("normal(0,1)", class="b",coef="TRT")
JJ.priors.list$ACPA <- set_prior("normal(0,1)", class="b",coef="ACPA")
JJ.priors.list$DRUG_WK <- set_prior("normal(0,.1)", class="b",coef="DRUG:WK")
JJ.priors.list$DRUG_ACPA <- set_prior("normal(0,1)", class="b",coef="DRUG:ACPA")
JJ.priors.list$ID.sd <- set_prior("cauchy(0,2)", class="sd",group="ID")
JJ.priors.list$cor <- set_prior("lkj(1.5)", class="cor")
JJ.priors.list$ar <- set_prior("normal(0,.75)", class="ar")

#############################################################
## LOAD PREVIOUSLY RUN MODELS ###############################
#############################################################

## FCT: Load Previously Compiled Models & Data
LOAD_PREV_MODS <- function( Path, Tag ) {
	Temp.files <- grep( "Rdata$",list.files( Path ), value=T )
	Temp.files.2 <- grep( "Model_Summary", Temp.files, value=T,invert=T )
	Temp.files.mods <- grep("Model",Temp.files.2, value=T,invert=F )
	Temp.files.meta <- grep("Model",Temp.files.2, value=T,invert=T )
	# for ( file in Temp.files.meta ) { load(paste(Path,file,sep="")) }
	OUT <- list()
	for ( file in Temp.files.mods ) {
		print(file)
		mod.tag <- gsub( "Rdata.Model.","", gsub(".Rdata","",file, fixed=T),fixed=T)
		print(mod.tag)
		load(paste(Path,file,sep=""))
		OUT[[mod.tag]] <- temp.model
	}
	return(OUT)
}

# ## Load Null Models (TRT)
# if ( !exists("JJ") ) { JJ <- list() }
# load.tag <- "Null_T"
# JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.T, load.tag )

# ## Load Null Models (TRT, w/ Corr)
# load.tag <- "Null_Tc"
# if ( !exists("JJ") ) { JJ <- list() }
# JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.Tc, load.tag )

## Load Null Models (TRT, no DRUG:WK)
load.tag <- "Null_Tn"
if ( !exists("JJ") ) { JJ <- list() }
JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.Tn, load.tag )

## Load Null Models (TRT, w/ Corr, no DRUG:WK)
load.tag <- "Null_Tnc"
if ( !exists("JJ") ) { JJ <- list() }
JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.Tnc, load.tag )

# ## Load SNP Models
# load.tag <- "SNP"
# if ( !exists("JJ") ) { JJ <- list() }
# JJ[[load.tag]] <- LOAD_PREV_MODS( PathToSNPMod.1, load.tag )

# ## Load Null Models (TRT, w/ Corr (DRUG:ID))
# load.tag <- "Null_Tcdi"
# if ( !exists("JJ") ) { JJ <- list() }
# JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.Tcdi, load.tag )
# ## Load Null Models (TRT, w/ Corr)
# load.tag <- "Null_Tcid"
# if ( !exists("JJ") ) { JJ <- list() }
# JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.Tcid, load.tag )

# if ( !( "m11i1"%in%names(JJ$Null_T) ) ) {
	## Null Model 11 
	# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160703_Null_Trt_m11/Rdata.Model.m11i1.Rdata" )
	# JJ$Null_T$m11i1 <- temp.model
	# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160703_Null_TrtCor_m11/Rdata.Model.m11i1.Rdata" )
	# JJ$Null_Tc$m11i1 <- temp.model
	# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160703_Null_TrtCorDrI_m11/Rdata.Model.m11i1.Rdata" )
	# JJ$Null_Tcdi$m11i1 <- temp.model
	# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160703_Null_TrtCorIDr_m11/Rdata.Model.m11i1.Rdata" )
	# JJ$Null_Tcid$m11i1 <- temp.model	
# }
# ## Load Down Sampled Models
# load.tag <- "down9"
# if ( !exists("JJ") ) { JJ <- list() }
# JJ[[load.tag]] <- LOAD_PREV_MODS( PathToDownMod.1, load.tag )
# load( paste(PathToDownMod.1,"Rdata.ModObs.",load.tag,".Rdata",sep="") )

## Load Individual Model Summary
load.tag <- "Ind4"
load.tag <- "Ind3"
load.tag.ind <- load.tag
load( PathToIndMod.1 )
JJ[[load.tag]] <- COMPILE

## Load Individual Model Summary
load.tag <- "LeadUp"
if ( !exists("JJ") ) { JJ <- list() }
JJ[[load.tag]] <- LOAD_PREV_MODS( PathToLeadUp, load.tag )

#############################################################
## CALCULATE MODEL FITS (WAIC) ##############################
#############################################################

## MODEL KEY!!!
 # Set Colors
KEY.cols <- COLS.mods
 # Set Tags
KEY.tag <- c("Null_Tn","Null_Tnc",load.tag.ind,"LeadUp")
KEY.tag.2 <- c("Tn","Tnc","BI","LeadUp")
 # Compile to Data Frame
KEY <- data.frame( COL=KEY.cols, TAG=KEY.tag, TAG.2=KEY.tag.2, stringsAsFactors=F )
rownames(KEY) <- KEY$TAG

## Calculate WAIC of Null Models
# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160702_Long_Compare/Rdata.Null_WAIC.Rdata" )
# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160703_Long_Compare/Rdata.Null_WAIC.Rdata" )
# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160704_Long_Compare/Rdata.Null_WAIC.Rdata" )
# load( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160705_Long_Compare/Rdata.Null_WAIC.Rdata" )
# load( "/Users/kstandis/Data/Janssen/Data/Longitude/20160705_Null_WAIC.Rdata" )
# load( "/Users/kstandis/Data/Janssen/Plots_Mac/20160719_Long_Compare/Rdata.Null_WAIC.Rdata" )
# load( "/Users/kstandis/Data/Janssen/Plots_Mac/20160720_Long_Compare/Rdata.Null_WAIC.Rdata" )
load( "/Users/kstandis/Data/Janssen/Plots_Mac/20160727_Long_Compare/Rdata.Null_WAIC.Rdata" )
if ( !exists("JJ.waic") ) { JJ.waic <- list() }
# if ( !exists("JJ.loo") ) { JJ.loo <- list() }
goals.temp <- grep("Null",names(JJ),value=T)
goals.temp <- "LeadUp"
for ( load.tag in goals.temp ) {
	if ( !(load.tag%in%names(JJ.waic)) ) { JJ.waic[[load.tag]] <- list() }
	for ( m.tag in names(JJ[[load.tag]]) ) {
		print(paste(load.tag,m.tag))
		if ( !(m.tag %in% names(JJ.waic[[load.tag]])) ) {
			JJ.waic[[load.tag]][[m.tag]] <- WAIC( JJ[[load.tag]][[m.tag]] )	
			# JJ.loo[[load.tag]][[m.tag]] <- LOO( JJ[[load.tag]][[m.tag]] )	
		}
	}
}
# save( JJ.waic, file=paste(PathToPlot,"Rdata.Null_WAIC.Rdata",sep="") )

## Plot Model Fits (WAIC)
 # Which Models/Cors to Plot
goals.temp <- c("Null_T","Null_Tc","Null_Tn","Null_Tnc")
goals.temp <- grep("Null",names(JJ),value=T)
mods.temp <- paste("m",1:11,"i1",sep="")
mods.temp <- mods.temp[-4]
 # Plotting Parameters
XLIM <- c(1,length(mods.temp)) # c(1,Reduce(max,lapply(JJ.waic,length)))
YLIM <- Reduce( range, lapply(JJ.waic,function(x)lapply(x,function(y)y$waic)) )
 # Plot
png( paste(PathToPlot,"Model_Selection.WAIC.png",sep=""), height=1000,width=1600,pointsize=32 )
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Model",ylab="WAIC",main="Model Selection: WAIC",xaxt="n" )
axis( 1, at=XLIM[1]:XLIM[2] )
abline(h=seq(0,100e3,2e3),v=0:20, lty=3,col="grey50",lwd=1 )
for ( load.tag in goals.temp ) {
	mods.temp.2 <- intersect( mods.temp, names(JJ.waic[[load.tag]]) )
	xvals <- which(mods.temp %in% mods.temp.2 ) # 1:length(mods.temp.2) # length(JJ.waic[[load.tag]])
	if ( grepl("_Tn",load.tag) ) { LTY=1 }else{ LTY=2 }
	waic.temp <- unlist(lapply(mods.temp.2,function(y)JJ.waic[[load.tag]][[y]]$waic))
	se.temp <- unlist(lapply(mods.temp.2,function(y)JJ.waic[[load.tag]][[y]]$se_waic))
	points( xvals, waic.temp, type="o",pch=16,col=KEY[load.tag,"COL"],lwd=4,lty=LTY )
	arrows( xvals, waic.temp+se.temp, xvals, waic.temp-se.temp, col=KEY[load.tag,"COL"],lwd=2,code=3,angle=90,length=.1 )
}
legend( "topright",pch=16,lty=1,lwd=4,col=KEY[goals.temp,"COL"],legend=goals.temp )
dev.off()

## Plot WAIC for LeadUp Models
Cols.temp <- COLS.list.ord[5]
mods.temp <- paste("m",9:11,sep="")
mods.temp.lab <- paste("m",9:11-1,sep="")
waic.temp <- unlist(lapply( JJ.waic$LeadUp, function(x)x$waic ))
se.temp <- unlist(lapply( JJ.waic$LeadUp, function(x)x$se_waic ))
YLIM <- range( c(waic.temp+se.temp,waic.temp-se.temp) )
XLIM <- c(1,length(mods.temp)) + c(-1,1)*.5
png( paste(PathToPlot,"LeadUp.WAIC.png",sep=""), height=1000,width=600,pointsize=32 )
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Model",ylab="WAIC",main="Model Selection:\nSimulated Baseline",xaxt="n" )
axis( 1, at=1:3, label=mods.temp.lab )
abline( h=seq(16e3,20e3,2e2),lty=2,col="grey50",lwd=1 )
points( 1:length(mods.temp), waic.temp[mods.temp], type="p",col=Cols.temp,pch=16, cex=1.5,lwd=4 )
arrows( 1:length(mods.temp), waic.temp[mods.temp]-se.temp[mods.temp], 1:length(mods.temp), waic.temp[mods.temp]+se.temp[mods.temp], angle=90,code=3,length=.1,col=Cols.temp,lwd=4 )
dev.off()

##############################################################
## COMPILE SUMMARY INFO FROM MODELS ##########################
##############################################################

######################################
## MODEL SUMMARIES ###################

## Get some info from models 9 to 11
# load("/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160705_Long_Compare/Rdata.M911.Rdata" )
# load( "/Users/kstandis/Data/Janssen/Data/Longitude/20160705_Null_M911.Rdata" )
# load( "/Users/kstandis/Data/Janssen/Plots_Mac/20160719_Long_Compare/Rdata.M911.Rdata" )
# load( "/Users/kstandis/Data/Janssen/Plots_Mac/20160727_Long_Compare/Rdata.M911.Rdata" )
load( "/Users/kstandis/Data/Janssen/Plots_Mac/20160823_Long_Compare/Rdata.M911.Rdata" )

# goals.temp <- grep("Null",names(JJ),value=T)
# mods.temp <- paste("m",9:11,"i1",sep="")
# goals.temp <- "LeadUp"
# mods.temp <- paste("m",9:11,sep="")
# JJ.sum <- JJ.res <- JJ.coef <- JJ.pred <- list()
# for ( load.tag in goals.temp ) {
# 	for ( m.tag in mods.temp ) {
# 	# for ( m.tag in "m11i1" ) {
# 		print(paste(load.tag,m.tag))
# 		# JJ.sum[[load.tag]][[m.tag]] <- summary( JJ[[load.tag]][[m.tag]] )
# 		JJ.res[[load.tag]][[m.tag]] <- resid( JJ[[load.tag]][[m.tag]] )
# 		JJ.coef[[load.tag]][[m.tag]] <- coef( JJ[[load.tag]][[m.tag]] )$ID
# 		JJ.pred[[load.tag]][[m.tag]] <- predict( JJ[[load.tag]][[m.tag]] )
# 	}
# }

## Compile & Save Model Summaries, Predictions, Residuals, & Coefficients
# M.9_11 <- list( Res=JJ.res, Pred=JJ.pred, Coef=JJ.coef )
# M.9_11$Res$Null_Tn <- JJ.res$Null_Tn
# M.9_11$Coef$Null_Tn <- JJ.coef$Null_Tn
# M.9_11$Pred$Null_Tn <- JJ.pred$Null_Tn
# M.9_11$Res$Null_Tnc <- JJ.res$Null_Tnc
# M.9_11$Coef$Null_Tnc <- JJ.coef$Null_Tnc
# M.9_11$Pred$Null_Tnc <- JJ.pred$Null_Tnc
# M.9_11$Res$LeadUp <- JJ.res$LeadUp
# M.9_11$Coef$LeadUp <- JJ.coef$LeadUp
# M.9_11$Pred$LeadUp <- JJ.pred$LeadUp
# save( M.9_11, file=paste(PathToPlot,"Rdata.M911.Rdata",sep="") )

## FCT: Pull Model Data
 # model <- JJ$Null_T$m9i1 ; tag <- "T_m9"
SUMM_MODEL <- function( model ) {
	write(paste(date(),"Model:"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Collect General Model Info
	m <- d <- list()
	summ <- summary(model,waic=F)
	m$obs <- summ$nobs
	m$iter <- summ$iter
	m$warm <- summ$warmup
	m$chain <- summ$chains
	# m$waic <- summ$WAIC
	m$pheno <- as.character(model$formula)[2]
	m$covs <- as.character(model$formula)[3]
	m$form <- paste( m$pheno, "~", m$covs )

	## Pull Model Outputs
	 # Distributions
	d$prior <- model$prior
	d$post <- posterior_samples(model)
	 # Fixed Effects
	f.eff <- summ$fixed
	m$covs.f <- rownames(f.eff)
	m$covs.f.tag <- paste("b_",m$covs.f,sep="")
	 # Correlation Effects
	c.eff <- summ$cor_pars
	m$covs.c <- m$covs.c.tag <- rownames(c.eff)
	 # Sigma Pheno
	s.eff <- summ$spec_pars
	m$covs.s <- rownames(s.eff)
	m$covs.s.tag <- gsub(")","", gsub("(","_", m$covs.s, fixed=T),fixed=T )
	 # Random Effects
	m$rand_TF <- length(summ$random)>0
	if ( m$rand_TF ) {
		m$r.grp <- summ$group
		m$r.ngrps <- summ$ngrps
		 # Individual Estimates
		r.eff.i <- ranef(model)[[1]]
		m$covs.r <- colnames(r.eff.i)[1:ncol(r.eff.i)]
		# m$covs.r.tag <- unlist(lapply( m$covs.r, function(x)paste( "sd_",m$r.grp,"_",x,sep="") ))
		m$covs.r.i.tag <- unlist(lapply( m$covs.r, function(x)paste( "r_",m$r.grp,"[",m$samps,",",x,"]",sep="") ))
		 # Population Distribution Metrics
		r.eff.p <- summ$random
		m$covs.r.p <- rownames(r.eff.p$ID)
		m$covs.r.p.tag <- gsub(")","", gsub("(","_ID_", gsub(",","_", m$covs.r.p, fixed=T),fixed=T),fixed=T )
		 # Samples
		m$samps <- rownames(r.eff.i)
		m$n.samps <- length(m$samps)
		## Compile Effects
		all.eff <- list( f.eff, c.eff, Reduce( rbind, r.eff.p ), s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","R","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		mod.fit$Eff[grep("cor(",rownames(mod.fit),fixed=T)] <- "RC"
		mod.covs.tags <- c(m$covs.f.tag, m$covs.c.tag, m$covs.r.p.tag, m$covs.s.tag )
		mod.fit <- data.frame( mod.fit, Var_Tag=mod.covs.tags )
	}else{
		all.eff <- list( f.eff, c.eff, s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		mod.covs.tags <- c(m$covs.f.tag, m$covs.c.tag, m$covs.s.tag )
		mod.fit <- data.frame( mod.fit, Var_Tag=mod.covs.tags )
	}

	## Compile Outputs
	COMPILE <- list( Summ=summ, Fit=mod.fit, D=d, M=m )
}

## Get Model Summaries
 # Null Models
goals.temp <- grep("Null",names(JJ),value=T)
mods.temp <- paste("m",9:11,"i1",sep="")
mods.temp <- paste("m",10:11,"i1",sep="")
 # Lead Up Models
goals.temp <- "LeadUp"
mods.temp <- paste("m",9:11,sep="")
 # Null & Leadup Models
goals.temp <- grep("Null|Lead",names(JJ),value=T)
mods.temp <- paste("m",9:11,"i1",sep="")
if ( !exists("MOD.DAT") ) { MOD.DAT <- list() }
for ( goal in goals.temp ) {
	for ( mod in mods.temp ) {
		tag <- paste( gsub("Null_","",goal), gsub("i1","",mod), sep="_" )
		print( tag )
		if ( grepl("Lead",goal ) ) { mod <- gsub("i1","",mod ) }
		MOD.DAT[[tag]] <- SUMM_MODEL( JJ[[goal]][[mod]] )
		MOD.DAT[[tag]]$Resid <- M.9_11$Res[[goal]][[mod]]
		MOD.DAT[[tag]]$Coef <- M.9_11$Coef[[goal]][[mod]]
		MOD.DAT[[tag]]$Pred <- M.9_11$Pred[[goal]][[mod]]
	}
}

## Write out Model Fit Tables
tags <- names(MOD.DAT)
for ( tag in tags ) {
	writeLines( print(xtable(MOD.DAT[[tag]]$Fit[,c(7,1:4)],digits=3), include.rownames=T ), con=paste(PathToPlot,"TAB-",tag,".LaTeX.txt",sep="") )
	write.table( MOD.DAT[[tag]]$Fit, paste(PathToPlot,"TAB-",tag,".csv",sep=""), sep=",",row.names=T,col.names=T,quote=F )
}

######################################
## COMPILE MODEL ESTIMATES ###########

## From 4 models:
 # "Mean" model (MN)
 # BayesLMM model (BL)
 # BayesLMMcor model (BLC)
 # BayesIndiv model (BI)

## Get Sample Names
Samps <- sort(as.character(unique( TAB$ID )))
N.Samps <- length(Samps)

## Pull Drug Estimates
goals.temp <- grep("Null",names(JJ),value=T)
mods.temp <- paste("m",9:11,"i1",sep="")
Which_Vars <- data.frame( Var=c("DRUG","Intercept","TRT","WK"), Var.tag=c("Drug","Int","Trt","Wk"), stringsAsFactors=F )
COEF <- list()
for ( var in Which_Vars$Var ) {
	var.tag <- Which_Vars$Var.tag[ Which_Vars$Var==var ]
	COEF[[var.tag]] <- data.frame( row.names=Samps )
	## Mean Model
	if ( var.tag=="Drug" ) { COEF[[var.tag]]$MN <- FUL$DEL_MNe_MN[ match(Samps,FUL$ID) ] }
	if ( var.tag=="Int" ) { COEF[[var.tag]]$MN <- FUL$DAS_BL_MN[ match(Samps,FUL$ID) ] }
	## Mixed Models
	 # Has "TRT" as Random Effect?
	if ( var.tag %in% c("Drug","Int","Wk") ) { mods.temp.2 <- mods.temp
	}else{ mods.temp.2 <- grep("m11|m8",mods.temp,invert=T,value=T) }
	for ( goal in goals.temp ) {
		for ( mod in mods.temp.2 ) {
			tag <- paste( gsub("Null_","",goal), gsub("i1","",mod), sep="_" )
			COEF[[var.tag]][[tag]] <- MOD.DAT[[tag]]$Coef[Samps,var]
		}
	}
	# Individual Models
	ind.tag <- grep("Ind",names(JJ),value=T)
	COEF[[var.tag]]$BI <- JJ[[ind.tag]]$EST[paste(tolower(ind.tag),Samps,sep="_"),var]
}

## Pull Sample Info
Samp.inf <- list()
Samp.inf$Samps <- Samps
Samp.inf$GRP <- as.character( FUL$GRP[ match(Samps,FUL$ID) ] )
Samp.inf$ACPA <- as.numeric( FUL$ACPA[ match(Samps,FUL$ID) ]=="Positive" )
Samp.inf$RF <- as.numeric( FUL$RF[ match(Samps,FUL$ID) ]=="Positive" )
Samp.inf$ACR50_28 <- as.character( FUL$ACR50_28wk[ match(Samps,FUL$ID) ] )# as.numeric( FUL$ACR50_28wk[ match(Samps,FUL$ID) ]=="Y" )
Samp.inf$ACR50_52 <- as.character( FUL$ACR50_52wk[ match(Samps,FUL$ID) ] )# as.numeric( FUL$ACR50_52wk[ match(Samps,FUL$ID) ]=="Y" )
Samp.inf$EULAR_28 <- as.character( FUL$EUL_28_BL[ match(Samps,FUL$ID) ] )
Samp.inf$EULAR_52 <- as.character( FUL$EUL_52_BL[ match(Samps,FUL$ID) ] )
Samp.inf$COUN <- as.character( FUL$COUN[ match(Samps,FUL$ID) ] )
Samp.inf$RACE <- as.character( FUL$RACE[ match(Samps,FUL$ID) ] )
Samp.inf$ETHN <- as.character( FUL$ETHN[ match(Samps,FUL$ID) ] )
Samp.inf$JNJ_ANC <- as.character( FUL$JNJ_ANC[ match(Samps,FUL$ID) ] )
Samp.inf$SEX <- as.character( FUL$SEX[ match(Samps,FUL$ID) ] )
for ( i in 1:length(Samp.inf) ) { names(Samp.inf[[i]]) <- Samps }

##############################################################
## PLOTTING MODEL ESTIMATES ##################################
##############################################################


######################################
## MODEL KEY #########################

## Pull out Tags for Candidate Models
tags.cand <- grep("_m1", grep("^T",colnames(COEF$Int),value=T), value=T )

## Specify (Random) Variable Info
VAR.names <- c("DRUG","Intercept","TRT","WK")
VAR.tags <- c("Drug","Int","Trt","Wk")
VAR.labels <- c("GOL Effect","Baseline DAS","Treatment Effect","Trajectory")
Which_Vars <- data.frame( VAR=VAR.names, TAG=VAR.tags, LAB=VAR.labels, stringsAsFactors=F )
rownames(Which_Vars) <- Which_Vars$VAR

## Create Key for Storing Model Meta Info
 # Combined Goal/Model Tag
KEY.tags.comb <- colnames(COEF$Int)
 # Goal Tag
KEY.tags.goal.1 <- sapply( strsplit( KEY.tags.comb, "_" ),"[",1 )
KEY.tags.goal <- paste( "Null",KEY.tags.goal.1,sep="_" )
KEY.tags.goal[grep("^T",KEY.tags.goal.1,invert=T)] <- "Ind"
 # Model Tag
KEY.tags.mod.1 <- sapply( strsplit( KEY.tags.comb, "_" ),"[",2 )
KEY.tags.mod.1[grep("^T",KEY.tags.goal.1,invert=T)] <- "Ind"
KEY.tags.mod <- paste(KEY.tags.mod.1,"i1",sep="")
 # Combine
KEY.2 <- data.frame( TAG.c=KEY.tags.comb, TAG.g=KEY.tags.goal, TAG.g1=KEY.tags.goal.1, TAG.m=KEY.tags.mod, TAG.m1=KEY.tags.mod.1, stringsAsFactors=F )
KEY.2 <- merge( KEY.2, KEY[,c("COL","TAG.2")], by.x="TAG.g1", by.y="TAG.2" )
rownames(KEY.2) <- KEY.2$TAG.c

## Add LeadUp Meta Data to KEY.2
KEY.2l <- c( "LeadUp","LeadUp_m9","LeadUp","m9","m9",COLS.mods["LeadUp"] )
KEY.2 <- Reduce( rbind, list( KEY.2, KEY.2l, gsub("m9","m10",KEY.2l), gsub("m9","m11",KEY.2l) ) )
rownames(KEY.2) <- KEY.2[,"TAG.c"]

######################################
## SCRAP PLOT IDEAS ##################

## Distributions of Random Effects across models & Variables
par(mfrow=c(2,2))
boxplot( COEF$Int, las=2, col=COLS.list.2[factor(sapply(strsplit( colnames(COEF$Int), "_" ),"[",1))] )
boxplot( COEF$Drug, las=2, col=COLS.list.2[factor(sapply(strsplit( colnames(COEF$Drug), "_" ),"[",1))] )
boxplot( COEF$Trt, las=2, col=COLS.list.2[factor(sapply(strsplit( colnames(COEF$Trt), "_" ),"[",1))] )
boxplot( COEF$Wk, las=2, col=COLS.list.2[factor(sapply(strsplit( colnames(COEF$Wk), "_" ),"[",1))] )

## Actual DAS values across ARM before and after treament
TEMP <- merge( TAB, FUL[,c("ID","GRP")] )
boxplot( DAS ~ DRUG + GRP, data=TEMP )

## Standard Deviation of the Residuals across Models
barplot( unlist(lapply( MOD.DAT, function(x)sd(x$Resid[,1]) )) )

## Correlation b/n Predicted and Actual Values
TAB.das <- CLIN_TABS$DAS
lapply( MOD.DAT,function(x)cor(x$Pred[,1],TAB.das$DAS) )
lapply( MOD.DAT,function(x)cor(x$Resid[,1],TAB.das$DAS) )
lapply( MOD.DAT,function(x)cor(x$Resid[,1],x$Pred[,1]) )

## Residuals of Candidate Models over the Weeks
par(mfrow=c(2,2))
Scrap <- lapply(tags.cand,function(x)boxplot( MOD.DAT[[x]]$Resid[,1] ~ TAB.das$WK, main=x) )

## Correlation b/n Residuals of Candidate Models
cor( Reduce( cbind, lapply( tags.cand,function(x) MOD.DAT[[x]]$Resid[,1] )) )

## Standard Deviation of Within-Patient Residuals b/n Models
TAB.das <- CLIN_TABS$DAS
MN.wi.pat.Tn_m10 <- unlist(lapply( Samp.inf$Samps, function(x)mean(MOD.DAT$Tn_m10$Resid[TAB.das$ID==x,1]) ))
MN.wi.pat.BI <- unlist(lapply( Samp.inf$Samps, function(x)mean(TAB.das$DAS[TAB.das$ID==x]-JJ$Ind3$PRED[[paste("ind3",x,sep="_")]][,1]) ))
SD.wi.pat.Tn_m10 <- unlist(lapply( Samp.inf$Samps, function(x)sd(MOD.DAT$Tn_m10$Resid[TAB.das$ID==x,1]) ))
SD.wi.pat.BI <- unlist(lapply( Samp.inf$Samps, function(x)sd(TAB.das$DAS[TAB.das$ID==x]-JJ$Ind3$PRED[[paste("ind3",x,sep="_")]][,1]) ))
par(mfrow=c(1,2))
plot( MN.wi.pat.Tn_m10 ~ MN.wi.pat.BI, col=adjustcolor(Cols.a[Samp.inf$GRP],.6),pch=16, main="MN" ) ; abline(0,1)
plot( SD.wi.pat.Tn_m10 ~ SD.wi.pat.BI, col=adjustcolor(Cols.a[Samp.inf$GRP],.6),pch=16, main="SD" ) ; abline(0,1)
head(sort( -abs(SD.wi.pat.Tn_m10 - SD.wi.pat.BI) ))

## Within Patient Error of Estimate (using Posterior Distributions)
POST.wi.sd <- list()
POST.wi.sd$int <- lapply( tags.cand,function(x) apply(MOD.DAT[[x]]$D$post[,paste("r_ID[",Samp.inf$Samps,",Intercept]",sep="")], 2, sd) )
POST.wi.sd$drug <- lapply( tags.cand,function(x) apply(MOD.DAT[[x]]$D$post[,paste("r_ID[",Samp.inf$Samps,",DRUG]",sep="")], 2, sd) )
POST.wi.sd$wk <- lapply( tags.cand,function(x) apply(MOD.DAT[[x]]$D$post[,paste("r_ID[",Samp.inf$Samps,",WK]",sep="")], 2, sd) )
for ( v in names(POST.wi.sd) ) { names(POST.wi.sd[[v]]) <- tags.cand }
POST.wi.sd$int$BI <- unlist(lapply( JJ$Ind3$POST_RAW, function(x)sd(x[,"b_Intercept"]) ))
POST.wi.sd$drug$BI <- unlist(lapply( JJ$Ind3$POST_RAW, function(x)sd(x[,"b_DRUG"]) ))
POST.wi.sd$wk$BI <- unlist(lapply( JJ$Ind3$POST_RAW, function(x)sd(x[,"b_WK"]) ))
# par(mfrow=c(3,5))
# Scrap <- lapply( names(POST.wi.sd),
# 	function(x)lapply( names(POST.wi.sd[[x]]),
# 		function(y)hist(POST.wi.sd[[x]][[y]],main=paste(x,y) )
# 	)
# )

temp.tags <- grep("Tnc_",names(POST.wi.sd$drug),value=T,invert=T)
temp.tags <- grep("Tnc_|11",names(POST.wi.sd$drug),value=T,invert=T)
par(mfrow=c(length(temp.tags),length(temp.tags)))
par(mar=c(4,4,3,.5))
lapply( temp.tags,
	function(x)lapply( temp.tags,
		function(y)plot(POST.wi.sd$drug[[x]] ~ POST.wi.sd$int[[y]], col=Cols.a[Samp.inf$GRP],main=paste(x,y),ylab=paste("Drug",x),xlab=paste("Int",y) )
		# function(y)plot(POST.wi.sd$drug[[x]] ~ POST.wi.sd$wk[[y]], col=Cols.a[Samp.inf$GRP],main=paste(x,y),ylab=paste("Drug",x),xlab=paste("Wk",y) )
		# function(y)plot(POST.wi.sd$int[[x]] ~ POST.wi.sd$wk[[y]], col=Cols.a[Samp.inf$GRP],main=paste(x,y),ylab=paste("Int",x),xlab=paste("Wk",y) )
	)
)


## Boxplot Random Effects vs Arms
par(mfrow=c(2,4))
layout(matrix(1:8,ncol=4))
boxplot( tag.var ~ GRP, data=mg.var, col=Cols.a,border=Cols.m[1],lwd=2 )
boxplot( vs.var ~ GRP, data=mg.var, col=Cols.a,border=Cols.m[2],lwd=2 )
boxplot( tag.var ~ GRP, data=mg.var, border=Cols.a,col=Cols.m[1],lwd=2 )
boxplot( vs.var ~ GRP, data=mg.var, border=Cols.a,col=Cols.m[2],lwd=2 )
boxplot( tag.var ~ GRP, data=mg.var, col=Cols.a,lwd=2 )
points( tag.var ~ GRP, data=mg.var, col=Cols.m.7[1],pch=16 )
boxplot( vs.var ~ GRP, data=mg.var, col=Cols.a,lwd=2 )
points( vs.var ~ GRP, data=mg.var, col=Cols.m.7[2],pch=16 )
boxplot( tag.var ~ GRP, data=mg.var, col=Cols.a,lwd=2 ) ; arrows(1,min(mg.var$tag.var),3,min(mg.var$tag.var),col=Cols.m[1],lwd=2,length=0 )
boxplot( vs.var ~ GRP, data=mg.var, col=Cols.a,lwd=2 ) ; arrows(1,min(mg.var$vs.var),3,min(mg.var$vs.var),col=Cols.m[2],lwd=2,length=0 )
apply( MOD.DAT$Tn_m10$Coef[,1:4], 2, function(x)anova(lm(x ~ Samp.inf$GRP))[1,"Pr(>F)"] )
apply( JJ$Ind3$EST[,1:4], 2, function(x)anova(lm(x ~ Samp.inf$GRP))[1,"Pr(>F)"] )
apply( MOD.DAT$Tn_m10$Coef[,1:4], 2, function(x)boxplot(x ~ Samp.inf$GRP,col=Cols.a) )
apply( JJ$Ind3$EST[,1:4], 2, function(x)boxplot(x ~ Samp.inf$GRP,col=Cols.a) )

######################################
## FUNCTIONS #########################

## FCT: More Model Selection Parameters
PLOT_MOD_SELECT <- function( mods.cand ) {
	mod.dat <- MOD.DAT[mods.cand]

	## Standard Deviation of the Residuals across Models
	mod.res.sd <- unlist(lapply( mod.dat, function(x)sd(x$Resid[,1]) ))
	goal.tags <- unique(sapply(strsplit(mods.cand,"_"),"[",1))
	mod.tags <- unique(sapply(strsplit(mods.cand,"_"),"[",2))
	dim(mod.res.sd) <- c(length(mod.tags),length(goal.tags))
	colnames(mod.res.sd) <- goal.tags
	rownames(mod.res.sd) <- mod.tags
	## Correlation b/n Predicted and Actual Values
	TAB.das <- CLIN_TABS$DAS
	mod.pred.cor <- unlist(lapply( mod.dat,function(x)cor(x$Pred[,1],TAB.das$DAS,method="pearson") ))
	dim(mod.pred.cor) <- c(length(mod.pred.cor)/2,2)
	colnames(mod.pred.cor) <- unique(sapply(strsplit(mods.cand,"_"),"[",1))
	rownames(mod.pred.cor) <- unique(sapply(strsplit(mods.cand,"_"),"[",2))
	
	## Plot it
	png( paste(PathToPlot,"Model_Selection.StDev.Cor.png",sep=""),height=1200,width=1200,pointsize=28)
	par(mfrow=c(1,2))
	 # S.D. of Residuals
	MAIN.1 <- "Model Residuals"
	YLAB.1 <- "St.Dev. of Residuals"
	barplot( t(mod.res.sd), ylim=c(0,1),ylab=YLAB.1,xlab="Model",main=MAIN.1, beside=T, col=KEY[match(colnames(mod.res.sd),KEY$TAG.2),"COL"],border=NA )
	legend("topright", title="Correlation Structure",legend=colnames(mod.res.sd),fill=KEY[match(colnames(mod.res.sd),KEY$TAG.2),"COL"],border=NA,ncol=2,bg="white")
	 # Correlation of Predicted vs Observed
	MAIN.2 <- "Pearson Correlation of\nPredicted & Observed Values"
	YLAB.2 <- "R"
	barplot( t(mod.pred.cor), ylim=c(0,1),ylab=YLAB.2,xlab="Model",main=MAIN.2, beside=T, col=KEY[match(colnames(mod.pred.cor),KEY$TAG.2),"COL"],border=NA )
	dev.off()
}
mods.cand <- grep("^T",names(MOD.DAT),value=T)
mods.cand <- grep("m1",mods.cand,value=T)
PLOT_MOD_SELECT( mods.cand )

## FCT: Plot Model Fixed Effects
PLOT_FIXED <- function( tag, which_plots ) {

	## Specify Model
	tag.goal <- KEY.2[tag,"TAG.g"]
	tag.mod <- KEY.2[tag,"TAG.m"]
	model <- JJ[[tag.goal]][[tag.mod]]
	## Pull out Model Summary
	if ( tag %in% names(MOD.DAT) ) {
		mod.dat <- MOD.DAT[[tag]]
	}else{ mod.dat <- SUMM_MODEL( model, tag ) }

	## Set Color Palette
	COLS.list.heat <- c("firebrick3","chocolate2","gold1","springgreen1","steelblue2","slateblue3")
	COLS <- adjustcolor(COLS.list.2,alpha=.6)
	COLS.eff <- COLS[c(1:4,6)]
	names(COLS.eff) <- c("F","C","R","RC","S")

	## 1 - Plot Model Summary
	if ( 1 %in% which_plots ) {	
		print("Plotting #1 - Model Summary")
		png( paste(PathToPlot,"ModSumm_",tag,".F1-Converged.png",sep=""),height=200*nrow(mod.dat$Fit),width=800,pointsize=30)
		par(mfrow=c(nrow(mod.dat$Fit),2))
		plot( model, N=nrow(mod.dat$Fit) )
		dev.off()
	}

	# ## 2 - Plot Chains
	#  # (aka) How to Pull all Sampling Data
	# print("Plotting #2 - Chains")
	# if ( 2 %in% which_plots ) {	
	# 	which_vars <- mod.dat$M$covs.f
	# 	which_vars <- grep( "sd(",rownames(mod.dat$Fit),fixed=T,value=T )

	# 	png( paste(PathToPlot,"ModSumm_",tag,".1-Chains.png",sep=""),height=200*nrow(mod.dat$Fit),width=800,pointsize=30)
	# 	par( mfrow=c(length(which_vars),1) )
	# 	for ( v in which_vars ) {
	# 		v.tag <- grep( substr(gsub("(",paste("_",r.grp,"_",sep=""),v,fixed=T),1,10), model$fit@sim$fnames_oi )
	# 		YLIM <- Reduce(range,lapply(model$fit@sim$samples,function(x)range(x[[paste("b",v,sep="_")]])))
	# 		plot( 0,0,type="n",xlim=c(0,m.iter),ylim=YLIM,xlab="Iteration",ylab=v,main="Chains" )
	# 		abline( h=-100:100,lty=3,col="grey50",lwd=1 )
	# 		# SCRAP <- lapply( seq(YLIM[1],YLIM[2],.025),function(x)abline(h=x,col=adjustcolor(COLS.list.2[4],alpha=2*dnorm(x,5,1)),lwd=2 ))
	# 		for ( c in 1:m.chain ) {
	# 			# TEMP <- model$fit@sim$samples[[c]]$b_Intercept
	# 			points( 1:m.iter, model$fit@sim$samples[[c]][[paste("b",v,sep="_")]], type="l",col=adjustcolor(COLS.list.2[c],alpha=.7),lwd=2 )
	# 		}
	# 	}
	# dev.off()
	
	## 3 - Fixed Effect Sizes
	 # including prior distributions
	print("Plotting #3 - Effect Sizes")

	if ( 3 %in% which_plots ) {	
		 # Plotting Parameters
		YLIM <- extendrange(mod.dat$Fit[,"Estimate"], f=.2)
		png( paste(PathToPlot,"ModSumm_",tag,".F3-EffSize.png",sep=""),height=1000,width=400+100*nrow(mod.dat$Fit),pointsize=26)
		par(mar=c( 7,5,5,3 ))
		TEMP <- 1:nrow(mod.dat$Fit)
		plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=range(TEMP),ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
		axis( 1, at=TEMP,label=rownames(mod.dat$Fit), las=2 )
		axis( 2, at=seq(-10,10,2), las=2 )
		abline( h=-10:10, lty=3,col="grey50",lwd=1 )
		abline( h=0, lty=1,col="black",lwd=1 )
		 # Plot Prior/Posterior Distributions
		for ( v in 1:nrow(mod.dat$Fit) ) {
			var <- rownames(mod.dat$Fit)[v]
			var.tag <- as.character( mod.dat$Fit$Var_Tag[v] )
			# Priors
			if ( var %in% mod.dat$M$covs.f ) {
				if ( var=="Intercept" ) {
					temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",mod.dat$D$prior[mod.dat$D$prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
				}else{
					temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",mod.dat$D$prior[mod.dat$D$prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
				}
				if ( length(temp.priors)==2 ) {
					temp.y <- seq(floor(YLIM[1]),ceiling(YLIM[2]),.01)
					temp.distr <- dnorm( temp.y, temp.priors[1], temp.priors[2] )
					temp.distr <- temp.distr / max(temp.distr) / 2.5
					polygon( c(TEMP[v]-temp.distr,rev(TEMP[v]+temp.distr)), c(temp.y,rev(temp.y)), col=adjustcolor("black",alpha=.2),border=NA )

				}
			}else{
				if ( grepl("^sd",var.tag) ) {
					temp.y <- seq(0,ceiling(YLIM[2]),.01)
					temp.distr <- dcauchy( temp.y, -1, 1 )
					temp.distr <- temp.distr / max(temp.distr) / 2.5
					polygon( c(TEMP[v]-temp.distr,rev(TEMP[v]+temp.distr)), c(temp.y,rev(temp.y)), col=adjustcolor("black",alpha=.2),border=NA )
				}
				if ( grepl("^ar",var.tag) ) {
					temp.y <- seq(-1,1,.01)
					temp.distr <- dunif( temp.y, -1, 1 )
					temp.distr <- temp.distr / max(temp.distr) / 2.5
					polygon( c(TEMP[v]-temp.distr,rev(TEMP[v]+temp.distr)), c(temp.y,rev(temp.y)), col=adjustcolor("black",alpha=.2),border=NA )
				}
				if ( grepl("^sigma",var.tag) ) {
					temp.y <- seq(0,ceiling(YLIM[2]),.01)
					temp.distr <- dt( temp.y, 3, 0 )
					temp.distr <- temp.distr / max(temp.distr) / 2.5
					polygon( c(TEMP[v]-temp.distr,rev(TEMP[v]+temp.distr)), c(temp.y,rev(temp.y)), col=adjustcolor("black",alpha=.2),border=NA )
				}
			}
			# Posteriors
			if ( var.tag %in% colnames(mod.dat$D$post) ) {
				vioplot( mod.dat$D$post[,var.tag], at=TEMP[v],col=COLS.eff[mod.dat$Fit[v,"Eff"]],add=T,drawRect=F )
			}
		}
		arrows( TEMP,mod.dat$Fit[,"l.95..CI"],TEMP,mod.dat$Fit[,"u.95..CI"],lwd=3,length=0 )
		arrows( TEMP-diff(TEMP)[1]/2*.6,mod.dat$Fit[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.dat$Fit[,"Estimate"],lwd=3,length=0 )
		legend("topright",fill=c(adjustcolor("black",alpha=.2),COLS.eff),border=NA,legend=c("Prior",names(COLS.eff)),ncol=(1+length(COLS.eff))/2,title="Effect Type",bg="white")
		dev.off()
	}
}

## FCT: Plot Shrinkage vs "Mean/BI" Model (1 Variable)
PLOT_SHRINK.1 <- function( tag, vs_tag, var, which_plots ) {

	## General Info
	 # Variable Name
	var.tag <- Which_Vars[ var, "TAG" ]
	var.lab <- Which_Vars[ var, "LAB" ]
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.m <- KEY.2[ match( c(tag,vs_tag), KEY.2$TAG.c ), "COL" ]
	Cols.m.7 <- adjustcolor( Cols.m, .7 )
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .7 )

	## Specify Model Info (tag)
	tag.goal <- KEY.2[tag,"TAG.g"]
	tag.mod <- KEY.2[tag,"TAG.m"]
	model <- JJ[[tag.goal]][[tag.mod]]
	mod.dat <- MOD.DAT[[tag]]
	 # Random Effects
	if ( var=="Intercept" ) { tag.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
	if ( var=="DRUG" ) { tag.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
	if ( !(var%in%c("DRUG","Intercept")) ) { tag.var <- mod.dat$Coef[,var] }
	 # Fixed Effect & Variance
	tag.fx <- mean( tag.var )
	tag.fx.sd <- sd( tag.var )

	## Specify Model Info (vs_tag)
	if ( vs_tag %in% c("MN","BI") ) {
		vs.var <- COEF[[var.tag]][,vs_tag]
	}else{
		vs.goal <- KEY.2[vs_tag,"TAG.g"]
		vs.mod <- KEY.2[vs_tag,"TAG.m"]
		model <- JJ[[vs.goal]][[vs.mod]]
		mod.dat <- MOD.DAT[[tag]]
		 # Random Effects
		if ( var=="Intercept" ) { vs.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
		if ( var=="DRUG" ) { vs.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
		if ( !(var%in%c("DRUG","Intercept")) ) { vs.var <- mod.dat$Coef[,var] }
	}
	 # Fixed Effect & Variance
	vs.fx <- mean( vs.var )
	vs.fx.sd <- sd( vs.var )

	## Merge Estimates from 2 Models
	mg.var <- data.frame( ID=Samp.inf$Samp, tag.var, vs.var, GRP=Samp.inf$GRP )
	mg.var.ord <- mg.var[ order(mg.var$tag.var), ]

	## Plot 1: Distributions
	if ( 1 %in% which_plots ) {
		XLIM <- range( c(mg.var$vs.var,mg.var$tag.var) )
		if ( var=="WK" ) { X_BIN <- .01 }else{ X_BIN <- .25 }
		BRKS <- seq( floor(XLIM[1]), ceiling(XLIM[2])+X_BIN, X_BIN)
		YLIM <- c( 0, Reduce( max, lapply(list(mg.var$tag.var,mg.var$vs.var),function(x)hist(x,breaks=BRKS,plot=F)$counts) ) )
		hist( mg.var$tag.var, col=Cols.m.7[1],border=NA,breaks=BRKS, xlim=XLIM,ylim=YLIM, main=paste(var.lab,"Distributions"),xlab=var.lab,ylab="# Patients" )
		hist( mg.var$vs.var, col=Cols.m.7[2],border=NA,breaks=BRKS, add=T )
		abline( v=0 )
		abline( v=c(tag.fx,vs.fx), col=Cols.m,lwd=4,lty=2 )
		arrows( tag.fx-tag.fx.sd, .0*YLIM[2], tag.fx+tag.fx.sd, .0*YLIM[2], col=Cols.m[1],lwd=4,angle=90,code=3,length=.1 )
		arrows( vs.fx-vs.fx.sd, .05*YLIM[2], vs.fx+vs.fx.sd, .05*YLIM[2], col=Cols.m[2],lwd=4,angle=90,code=3,length=.1 )
		legend( "topleft", legend=c(tag,vs_tag),title="Model",fill=Cols.m,border=NA)	
		text( XLIM[1], .74*YLIM[2], paste(round(tag.fx,2),"\u00B1",round(tag.fx.sd,2)), col=Cols.m[1], pos=4 )
		text( XLIM[1], .68*YLIM[2], paste(round(vs.fx,2),"\u00B1",round(vs.fx.sd,2)), col=Cols.m[2], pos=4 )
	}

	## Plot 2: Boxplot by Arm
	if ( 2 %in% which_plots ) {
		XLIM <- c(0,7)
		YLIM <- range( c(mg.var$vs.var,mg.var$tag.var) )
		plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",xlab="ARM",ylab=var,main=paste(var,"Distribution by Arm") )
		abline( h=-10:10,lty=3,col="grey50",lwd=1 )
		abline( h=0 )
		axis( 1,at=seq(1.5,7,2),label=ARMS )
		SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$tag.var[mg.var$GRP==ARMS[x]],at=2*x-1,col=Cols.m.7[1],border=NA,add=T) )
		SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$vs.var[mg.var$GRP==ARMS[x]],at=2*x,col=Cols.m.7[2],border=NA,add=T) )
		abline( h=c(tag.fx,vs.fx), col=Cols.m,lwd=4,lty=2 )
		# legend( "topleft", legend=c(tag,vs_tag),title="Model",fill=COLS.a[1:2],border=NA)	
		arrows( c(1,3,5),rep(YLIM[1],3),c(2,4,6),rep(YLIM[1],3), length=0,lwd=4,col=Cols.a )
	}

	## Plot 1 & 2: Distributions
	if ( 12 %in% which_plots ) {
		XLIM <- range( c(mg.var$vs.var,mg.var$tag.var) )
		if ( var=="WK" ) { X_BIN <- .01 }else{ X_BIN <- .25 }
		BRKS <- seq( floor(XLIM[1]), ceiling(XLIM[2])+X_BIN, X_BIN)
		YLIM <- c( 0, Reduce( max, lapply(list(mg.var$tag.var,mg.var$vs.var),function(x)hist(x,breaks=BRKS,plot=F)$counts) ) )
		YLIM[1] <- -.3*YLIM[2]
		 # Plot Distributions
		hist( mg.var$tag.var, col=Cols.m.7[1],border=NA,breaks=BRKS, xlim=XLIM,ylim=YLIM, main=paste(var,"Distributions"),xlab=var,ylab="# Patients",yaxt="n" )
		axis( 2, at=seq(0,YLIM[2],20), las=2 )
		hist( mg.var$vs.var, col=Cols.m.7[2],border=NA,breaks=BRKS, add=T )
		 # Add Lines for 0, means, sds
		abline( v=0 )
		abline( v=c(tag.fx,vs.fx), col=Cols.m,lwd=4,lty=2 )
		arrows( tag.fx-tag.fx.sd, .0*YLIM[2], tag.fx+tag.fx.sd, .0*YLIM[2], col=Cols.m[1],lwd=4,angle=90,code=3,length=.1 )
		arrows( vs.fx-vs.fx.sd, .05*YLIM[2], vs.fx+vs.fx.sd, .05*YLIM[2], col=Cols.m[2],lwd=4,angle=90,code=3,length=.1 )
		 # Legend
		legend( "topleft", legend=c(tag,vs_tag),title="Model",fill=Cols.m,border=NA)	
		 # Write Mean +/- SD
		text( XLIM[1], .7*YLIM[2], paste(round(tag.fx,2),"\u00B1",round(tag.fx.sd,2)), col=Cols.m[1], pos=4 )
		text( XLIM[1], .65*YLIM[2], paste(round(vs.fx,2),"\u00B1",round(vs.fx.sd,2)), col=Cols.m[2], pos=4 )
		 # Arm Distributions
		SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$tag.var[mg.var$GRP==ARMS[x]],at=YLIM[1]/6*(2*x-1),col=Cols.m.7[1],colMed=Cols.a[x],border=NA,add=T,horizontal=T,wex=6) )
		SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$vs.var[mg.var$GRP==ARMS[x]],at=YLIM[1]/6*(2*x),col=Cols.m.7[2],colMed=Cols.a[x],border=NA,add=T,horizontal=T,wex=6) )
		# legend( "topleft", legend=c(tag,vs_tag),title="Model",fill=COLS.a[1:2],border=NA)	
		# arrows( rep(XLIM[1],3),YLIM[1]/6*((2*1:3)-1),rep(XLIM[1],3),YLIM[1]/6*((2*1:3)), length=0,lwd=4,col=COLS[3:5] )
	}

	## Plot 3: MN vs BRMS Model Estimates
	if ( 3 %in% which_plots ) {
		XLIM <- range( mg.var$vs.var )
		YLIM <- range( mg.var$tag.var )
		X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5 )
		Y_LNS <- seq( floor(YLIM[1]-5), ceiling(YLIM[2]+5), .5 )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=paste("Individual Estimate:",vs_tag),ylab=paste("Mixed Model Estimate:",tag),main=paste(var.lab,"Correlations") )
		abline( v=X_LNS,h=Y_LNS, lty=3,col="grey50",lwd=1 )
		abline( h=0,v=0 )
		abline( 0,1 )
		points( mg.var$tag.var ~ mg.var$vs.var, data=mg.var, pch=16,cex=1,col=Cols.a.7[factor(GRP)] )
		abline( h=tag.fx, col=Cols.m[1],lwd=4,lty=2 )
		abline( v=vs.fx, col=Cols.m[2],lwd=4,lty=2 )
		legend( "topleft",pch=16,col=Cols.a,legend=levels(factor(mg.var$GRP)),title="Trial Arm",bg="white")
		temp.cor <- round( cor(mg.var$tag.var,mg.var$vs.var,method="pearson"), 3 )
		text( XLIM[1], quantile(YLIM,.65), paste("R =",temp.cor), col="black", pos=4 )
	}

	## Plot 4: Shrinkage by Person
	if ( 4 %in% which_plots ) {
		XLIM <- range( c(mg.var.ord$vs.var,mg.var.ord$tag.var) )
		X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5)
		YLIM <- c( 1, nrow(mg.var) )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=var,ylab="Patient",main=paste(var,"- Estimate Shrinkage") )
		abline( v=0 )
		abline( v=X_LNS, lty=3,col="grey50",lwd=1 )
		arrows( mg.var.ord$vs.var,1:nrow(mg.var.ord),mg.var.ord$tag.var,1:nrow(mg.var.ord), angle=30,length=.2,lwd=2,col=Cols.a[factor(mg.var.ord$GRP)] )
		points( mg.var.ord$tag.var,1:nrow(mg.var.ord), pch=16,cex=1,col=Cols.m.7[1] )
		points( mg.var.ord$vs.var,1:nrow(mg.var.ord), pch=16,cex=1,col=Cols.m.7[2] )
		abline( v=c(tag.fx,vs.fx), col=Cols.m,lwd=4,lty=2 )
		legend("topleft",col=c(Cols.m,Cols.a),legend=c(paste("Mod:",c(tag,vs_tag)),paste("Arm:",ARMS)),pch=c(16,16,NA,NA,NA),lty=c(0,0,1,1,1),lwd=2 )
	}
}

## FCT: Plot Shrinkage vs "Mean/BI" Model (2 Variables)
PLOT_SHRINK.2 <- function( tag, vs_tag, var_1, var_2, which_plots ) {

	## General Info
	 # Variable Name
	var.tag_1 <- Which_Vars[ var_1, "TAG" ]
	var.lab_1 <- Which_Vars[ var_1, "LAB" ]
	var.tag_2 <- Which_Vars[ var_2, "TAG" ]
	var.lab_2 <- Which_Vars[ var_2, "LAB" ]
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.m <- KEY.2[ match( c(tag,vs_tag), KEY.2$TAG.c ), "COL" ]
	Cols.m.7 <- adjustcolor( Cols.m, .7 )
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .8 )
	Cols.a.2 <- adjustcolor( Cols.a, .2 )

	## Specify Model Info (tag)
	tag.goal <- KEY.2[tag,"TAG.g"]
	tag.mod <- KEY.2[tag,"TAG.m"]
	model <- JJ[[tag.goal]][[tag.mod]]
	mod.dat <- MOD.DAT[[tag]]
	 # Random Effects
	if ( var_1=="Intercept" ) { tag.var_1 <- mod.dat$Coef[,var_1] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
	if ( var_1=="DRUG" ) { tag.var_1 <- mod.dat$Coef[,var_1] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
	if ( !(var_1%in%c("DRUG","Intercept")) ) { tag.var_1 <- mod.dat$Coef[,var_1] }
	if ( var_2=="Intercept" ) { tag.var_2 <- mod.dat$Coef[,var_2] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
	if ( var_2=="DRUG" ) { tag.var_2 <- mod.dat$Coef[,var_2] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
	if ( !(var_2%in%c("DRUG","Intercept")) ) { tag.var_2 <- mod.dat$Coef[,var_2] }
	 # Fixed Effect & Variance
	tag.fx_1 <- mean( tag.var_1 )
	tag.fx_1.sd <- sd( tag.var_1 )
	tag.fx_2 <- mean( tag.var_2 )
	tag.fx_2.sd <- sd( tag.var_2 )

	## Specify Model Info (vs_tag)
	if ( vs_tag %in% c("MN","BI") ) {
		vs.var_1 <- COEF[[var.tag_1]][,vs_tag]
		vs.var_2 <- COEF[[var.tag_2]][,vs_tag]
	}else{
		vs.goal <- KEY.2[vs_tag,"TAG.g"]
		vs.mod <- KEY.2[vs_tag,"TAG.m"]
		model <- JJ[[vs.goal]][[vs.mod]]
		mod.dat <- MOD.DAT[[tag]]
		 # Random Effects
		if ( var_1=="Intercept" ) { vs.var_1 <- mod.dat$Coef[,var_1] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
		if ( var_1=="DRUG" ) { vs.var_1 <- mod.dat$Coef[,var_1] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
		if ( !(var_1%in%c("DRUG","Intercept")) ) { vs.var_1 <- mod.dat$Coef[,var_1] }
		if ( var_2=="Intercept" ) { vs.var_2 <- mod.dat$Coef[,var_2] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
		if ( var_2=="DRUG" ) { vs.var_2 <- mod.dat$Coef[,var_2] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
		if ( !(var_2%in%c("DRUG","Intercept")) ) { vs.var_2 <- mod.dat$Coef[,var_2] }
	}
	 # Fixed Effect & Variance
	vs.fx_1 <- mean( vs.var_1 )
	vs.fx_1.sd <- sd( vs.var_1 )
	vs.fx_2 <- mean( vs.var_2 )
	vs.fx_2.sd <- sd( vs.var_2 )

	## Merge Estimates from 2 Models
	mg.var <- data.frame( ID=Samp.inf$Samp, tag.var_1, vs.var_1, tag.var_2, vs.var_2, GRP=Samp.inf$GRP )

	## Plot Shrinkage in 2 Dimensions
	if ( 1 %in% which_plots ) {
		png( paste(PathToPlot,"ModSumm_",tag,".Shrink2.",vs_tag,".",var_1,"_",var_2,".png",sep=""),height=1200,width=1200,pointsize=32)
		XLIM <- range( mg.var[,c("tag.var_1","vs.var_1")] )
		YLIM <- range( mg.var[,c("tag.var_2","vs.var_2")] )
		X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5 )
		Y_LNS <- seq( floor(YLIM[1]-5), ceiling(YLIM[2]+5), .5 )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=var_1,ylab=var_2,main=paste(tag,"v",vs_tag) )
		abline( v=X_LNS,h=Y_LNS, lty=3,col="grey50",lwd=1 )
		abline( 0,1 )
		abline( h=0,v=0 )
		# points( mg.var$vs.var_2 ~ mg.var$vs.var_1, data=mg.var, pch=16,cex=1,col=COLS[2] )
		# arrows( mg.var$vs.var_1, mg.var$vs.var_2, mg.var$tag.var_1, mg.var$tag.var_2, angle=30,length=_1,lwd=1,col=adjustcolor("grey50",.7) )
		# points( mg.var$tag.var_2 ~ mg.var$tag.var_1, data=mg.var, pch=16,cex=1,col=COLS[1] )
		# legend( "topleft", legend=c(tag,vs_tag),title="Model",pch=16,col=COLS[1:2])	
		points( mg.var$vs.var_2 ~ mg.var$vs.var_1, data=mg.var, pch=16,cex=1,col=Cols.a.2[factor(Samp.inf$GRP)] ) # col="grey80"
		arrows( mg.var$vs.var_1, mg.var$vs.var_2, mg.var$tag.var_1, mg.var$tag.var_2, angle=30,length=.1,lwd=1,col=adjustcolor("grey70",.7) )
		points( mg.var$tag.var_2 ~ mg.var$tag.var_1, data=mg.var, pch=16,cex=1,col=Cols.a.7[factor(Samp.inf$GRP)] )
		legend( "topleft", legend=ARMS,title="Arm",pch=16,col=Cols.a)	
		legend( "topright", legend=c(tag,vs_tag),title="Model",pch=NA,lty=2,lwd=4,col=Cols.m)
		 # Population Means
		abline( v=c(tag.fx_1,vs.fx_1), col=Cols.m,lwd=4,lty=2 )
		abline( h=c(tag.fx_2,vs.fx_2), col=Cols.m,lwd=4,lty=2 )
		dev.off()
	}
}

## FCT: Plot Random Effects (Individual Estimates)
PLOT_RAND <- function( tag, which_plots ) {

	## General Info
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .8 )
	Cols.a.2 <- adjustcolor( Cols.a, .2 )

	## Specify Model Info (tag)
	## Specify Model Info (vs_tag)
	if ( tag %in% c("MN","BI") ) {
		ind.tag <- grep("Ind",names(JJ),value=T)
		temp.reff <- JJ[[ind.tag]]$EST
		rownames(temp.reff) <- Samp.inf$Samps
		covs.r <- colnames(temp.reff)
		n.covs.r <- length(covs.r)
	}else{
		tag.goal <- KEY.2[tag,"TAG.g"]
		tag.mod <- KEY.2[tag,"TAG.m"]
		model <- JJ[[tag.goal]][[tag.mod]]
		mod.dat <- MOD.DAT[[tag]]
		covs.r <- mod.dat$M$covs.r
		n.covs.r <- length(covs.r)
		 # Random Effects
		temp.reff <- data.frame( row.names=rownames(mod.dat$Coef) )
		for ( var in covs.r ) {
			if ( var=="Intercept" ) { tag.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
			if ( var=="DRUG" ) { tag.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
			if ( !(var%in%c("DRUG","Intercept")) ) { tag.var <- mod.dat$Coef[,var] }	
			temp.reff <- data.frame( temp.reff, tag.var )
			colnames(temp.reff)[ncol(temp.reff)] <- var
		}
	}
	 # Add Arm
	temp.reff <- data.frame( temp.reff, GRP=Samp.inf$GRP )


	## Plot Pairs of Random Effects
	if ( 1 %in% which_plots ) {
		panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...) {
			usr <- par("usr"); on.exit(par(usr))
			par(usr = c(0, 1, 0, 1))
			r <- abs(cor(x, y))
			r <- cor(x, y)
			txt <- format(c(r, 0.123456789), digits = digits)[1]
			txt <- paste0(prefix, txt)
			# if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
			cex.cor <- 5
			text(0.5, 0.5, txt, cex = cex.cor * sqrt(abs(r)) )
		}

		print("Plotting #R1 - Pairs")
		if ( n.covs.r > 1 ) {
			# temp.reff <- mod.dat$Coef[ samps ,covs.r ]
			# temp.grp <- FUL$GRP[ match(rownames(temp.reff),FUL$ID) ] # col=adjustcolor(COLS.list.2[factor(temp.grp)],alpha=.4)
			# temp.reff <- merge( FUL[,c("ID","GRP")], temp.reff, by.x="ID",by.y="row.names" )
			# temp.reff <- data.frame( temp.reff, GRP=Samp.inf$GRP[samps] )
			png( paste(PathToPlot,"ModSumm_",tag,".R1-PairRand.png",sep=""),height=400*n.covs.r,width=400*n.covs.r,pointsize=24+2*n.covs.r)
			# pairs( temp.reff, pch=16,col=adjustcolor(COLS.eff["RC"],alpha=.5),cex=2,upper.panel=panel.cor )
			pairs( temp.reff[,covs.r], pch=16,col=Cols.a.7[factor(temp.reff$GRP)],cex=2,upper.panel=panel.cor )
			dev.off()
		}
	}

	## Plot Distributions w/ Variances of Random Effects of each Variable
}

## FCT: Plot Posterior Probabilities of Clinically Meaningful Response
PLOT_HEAT_COR <- function( tag, vs_tag, vars, subset ) {

	## General Info
	 # Variant Info
	var.tags <- Which_Vars$TAG[ match(vars,Which_Vars$VAR) ]
	var.labs <- Which_Vars$LAB[ match(vars,Which_Vars$VAR) ]
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.heat.list <- c(COLS.list.ord[5],COLS.arm)
	names(Cols.heat.list) <- c("All",ARMS)
	 # Patient Info
	if ( subset %in% Samp.inf$GRP ) {
		Samps.temp <- Samps[ Samp.inf$GRP==subset ]
		Cols.heat <- colorRampPalette(c("black",Cols.heat.list[subset]))(100)
	}else{
		Samps.temp <- Samps
		Cols.heat <- colorRampPalette(c("black",Cols.heat.list["All"]))(100)
	}

	## Specify Model Info (tag)
	tag.goal <- KEY.2[tag,"TAG.g"]
	tag.mod <- KEY.2[tag,"TAG.m"]
	model <- JJ[[tag.goal]][[tag.mod]]
	mod.dat <- MOD.DAT[[tag]]
	covs.r <- intersect( vars, mod.dat$M$covs.r )
	n.covs.r <- length(covs.r)
	var.tags.1 <- Which_Vars$TAG[ match(covs.r,Which_Vars$VAR) ]
	## Pull Model Coefficients
	temp.reff <- FUL[ FUL$ID%in%Samps.temp, c("ID","DAS_0wk") ]
	for ( v in 1:n.covs.r ) {
		var <- covs.r[v]
		var.tag <- var.tags.1[v]
		if ( var=="Intercept" ) { tag.var <- mod.dat$Coef[Samps.temp,var] + mod.dat$Coef[Samps.temp,"ACPA"]*Samp.inf$ACPA[Samps.temp] }
		if ( var=="DRUG" ) { tag.var <- mod.dat$Coef[Samps.temp,var] + mod.dat$Coef[Samps.temp,"DRUG:ACPA"]*Samp.inf$ACPA[Samps.temp] }
		if ( !(var%in%c("DRUG","Intercept")) ) { tag.var <- mod.dat$Coef[Samps.temp,var] }	
		temp.reff <- data.frame( temp.reff, tag.var )
		colnames(temp.reff)[ncol(temp.reff)] <- paste(var.tag,tag,sep=".")
	}

	## Specify Model Info (vs_tag)
	if ( vs_tag %in% c("MN","BI") ) {
		for ( var in var.tags ) {
			temp.tab <- data.frame( ID=rownames(COEF[[var]]), COEF[[var]][,vs_tag] )
			temp.reff <- merge( temp.reff, temp.tab, by="ID" )
			colnames(temp.reff)[ncol(temp.reff)] <- paste(var,vs_tag,sep=".")
		}
		# vs.var_1 <- COEF[[var.tag_1]][,vs_tag]
	}else{
		vs.goal <- KEY.2[vs_tag,"TAG.g"]
		vs.mod <- KEY.2[vs_tag,"TAG.m"]
		model <- JJ[[vs.goal]][[vs.mod]]
		mod.dat <- MOD.DAT[[vs_tag]]
		covs.r <- intersect( vars, mod.dat$M$covs.r )
		n.covs.r <- length(covs.r)
		var.tags.2 <- Which_Vars$TAG[ match(covs.r,Which_Vars$VAR) ]
		## Pull Model Coefficients
		for ( v in 1:n.covs.r ) {
			var <- covs.r[v]
			var.tag <- var.tags.1[v]
			if ( var=="Intercept" ) { vs.var <- mod.dat$Coef[Samps.temp,var] + mod.dat$Coef[Samps.temp,"ACPA"]*Samp.inf$ACPA[Samps.temp] }
			if ( var=="DRUG" ) { vs.var <- mod.dat$Coef[Samps.temp,var] + mod.dat$Coef[Samps.temp,"DRUG:ACPA"]*Samp.inf$ACPA[Samps.temp] }
			if ( !(var%in%c("DRUG","Intercept")) ) { vs.var <- mod.dat$Coef[Samps.temp,var] }	
			temp.reff <- data.frame( temp.reff, vs.var )
			colnames(temp.reff)[ncol(temp.reff)] <- paste(var.tag,tag,sep=".")
		}
	}

	## Pull out Model Coefficients
	# MG.coefs <- FUL[ FUL$ID%in%Samps.temp, c("ID","DAS_0wk") ]
	# for ( var in var.tags ) {
	# 	MG.coefs <- merge( MG.coefs, COEF[[var]][,c(tag,vs_tag)], by.x="ID",by.y="row.names" )
	# 	colnames(MG.coefs)[match(c(tag,vs_tag),colnames(MG.coefs))] <- paste(var,c(tag,vs_tag),sep=".")
	# }

	## Calculate Correlation b/n Variables
	# temp.reff <- temp.reff[ ,c(1:2,order(colnames(temp.reff)[-(1:2)])) ]
	temp.reff.ord <- c("ID","DAS_0wk", sort(grep("Int",colnames(temp.reff),value=T)) )
	temp.reff.ord <- c(temp.reff.ord, sort(setdiff(colnames(temp.reff),temp.reff.ord)) )
	# temp.reff <- temp.reff[ ,c("ID","DAS_0wk",sort(colnames(temp.reff)[-(1:2)])) ]
	temp.reff <- temp.reff[ , temp.reff.ord ]
	MG.cor <- cor( temp.reff[,-1], method="spearman" )

	## Plot Heatmap of Correlations
	 # Side Colors
 	Cols.side <- rep( KEY.2[tag,"COL"],ncol(MG.cor) )
	Cols.side[1] <- "white"
	Cols.side[grep("BI",colnames(MG.cor))] <- KEY.2["BI","COL"] # Cols.mod["Ind"]
	 # Labels/Breaks
	MAIN <- paste("Estimate Correlation:",subset)
	BRKS <- seq(0,1,length.out=101)
	 # Heatmap
	heatmap.2( abs(MG.cor), main=MAIN, col=Cols.heat,breaks=BRKS,
		Colv=NA,Rowv=NA,dendrogram="none",scale="none",trace="none",symm=T,
		lhei=c(1,5),lwid=c(1,6),key=F,margins=c(8,8),
		cellnote=round(MG.cor,2),notecol="white",notecex=1.7,
		RowSideColors=Cols.side,ColSideColors=Cols.side )

	## Return Correlation Matrix
	return( MG.cor )
}

## FCT: Plot Within Patient Error of Estimates
PLOT_WI_PAT_ERR <- function( tag, vs_tag, var_1, var_2 ) {
	## General Info
	 # Variable Name
	var.tag_1 <- Which_Vars[ var_1, "TAG" ]
	var.lab_1 <- Which_Vars[ var_1, "LAB" ]
	var.tag_2 <- Which_Vars[ var_2, "TAG" ]
	var.lab_2 <- Which_Vars[ var_2, "LAB" ]
	 # Model Tags
	tags.cand <- c(tag,vs_tag)
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .8 )
	Cols.a.2 <- adjustcolor( Cols.a, .2 )

	## Calculate Within Patient Error of Estimate (using Posterior Distributions)
	POST.wi.sd <- list()
	POST.wi.sd[[paste(var_1,tag,sep=".")]] <- apply(MOD.DAT[[tag]]$D$post[,paste("r_ID[",Samp.inf$Samps,",",var_1,"]",sep="")], 2, sd)
	POST.wi.sd[[paste(var_2,tag,sep=".")]] <- apply(MOD.DAT[[tag]]$D$post[,paste("r_ID[",Samp.inf$Samps,",",var_2,"]",sep="")], 2, sd)
	POST.wi.sd[[paste(var_1,vs_tag,sep=".")]] <- unlist(lapply( JJ$Ind3$POST_RAW, function(x)sd(x[,paste("b_",var_1,sep="")]) ))
	POST.wi.sd[[paste(var_2,vs_tag,sep=".")]] <- unlist(lapply( JJ$Ind3$POST_RAW, function(x)sd(x[,paste("b_",var_2,sep="")]) ))

	## Plot it
	png( paste(PathToPlot,"ModSumm_",tag,".WithinVar.vs.",vs_tag,"_",var_1,".",var_2,".png",sep=""),height=800,width=1600,pointsize=28)
	par(mfrow=c(1,2))
	 # First Model
	XLIM <- c(0,max( POST.wi.sd[[paste(var_1,tag,sep=".")]] ))
	YLIM <- c(0,max( POST.wi.sd[[paste(var_2,tag,sep=".")]] ))
	XLIM <- YLIM <- range(XLIM,YLIM)
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main=paste("Within-Patient Error:",tag),ylab=paste("St.Dev",var.lab_2),xlab=paste("St.Dev",var.lab_1) )
	abline( h=seq(0,10,.1),v=seq(0,10,.1), lty=3,col="grey50",lwd=1 )
	abline( 0,1 )
	abline( h=0,v=0 )
	points(POST.wi.sd[[paste(var_2,tag,sep=".")]] ~ POST.wi.sd[[paste(var_1,tag,sep=".")]], pch=16,col=Cols.a.7[factor(Samp.inf$GRP)] )
	legend( "topleft", legend=ARMS,col=Cols.a,pch=16, title="Trial Arm" )
	 # Second Model
	XLIM <- c(0,max( POST.wi.sd[[paste(var_1,vs_tag,sep=".")]] ))
	YLIM <- c(0,max( POST.wi.sd[[paste(var_2,vs_tag,sep=".")]] ))
	XLIM <- YLIM <- range(XLIM,YLIM)
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main=paste("Within-Patient Error:",vs_tag),ylab=paste("St.Dev",var.lab_2),xlab=paste("St.Dev",var.lab_1) )
	abline( h=seq(0,10,.1),v=seq(0,10,.1), lty=3,col="grey50",lwd=1 )
	abline( 0,1 )
	abline( h=0,v=0 )
	points(POST.wi.sd[[paste(var_2,vs_tag,sep=".")]] ~ POST.wi.sd[[paste(var_1,vs_tag,sep=".")]], pch=16,col=Cols.a.7[factor(Samp.inf$GRP)],main=vs_tag,ylab=paste(var_2,vs_tag),xlab=paste(var_1,vs_tag) )
	dev.off()
}
PLOT_WI_PAT_ERR( tag, vs_tag, var_1, var_2 )

## FCT: Plot Posterior Probabilities of Clinically Meaningful Response
PLOT_POSTERIORS <- function( tag, vs_tag, var, thresh ) {
	## General Info
	 # Variable Name
	var.tag <- Which_Vars[ var, "TAG" ]
	var.lab <- Which_Vars[ var, "LAB" ]
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .8 )
	Cols.m <- KEY.2[c(tag,vs_tag),"COL"] # COLS.mods
	Cols.m.7 <- adjustcolor( Cols.m, .8 )

	## Pull out Model Summary
	if ( tag %in% names(MOD.DAT) ) {
		mod.dat.1 <- MOD.DAT[[tag]]
		post.p.1 <- mod.dat.1$D$post[ ,paste("b_",var,sep="") ]
		post.i.1 <- mod.dat.1$D$post[ ,paste( "r_ID[",Samps,",",var,"]",sep="" ) ]
		which_acpa <- grep(paste(Samp.inf$Samps[Samp.inf$ACPA==1],collapse="|"),colnames(post.i.1))
		if ( var=="Intercept" ) {
			post.acpa <- mod.dat.1$D$post[ ,"b_ACPA"]
			post.i.1[,which_acpa] <- post.i.1[,which_acpa] + post.acpa
		}
		if ( var=="DRUG" ) {
			post.acpa <- mod.dat.1$D$post[ ,"b_DRUG:ACPA"]
			post.i.1[,which_acpa] <- post.i.1[,which_acpa] + post.acpa
		}
		post.1 <- post.p.1 + post.i.1
		prob.1 <- apply( post.1, 2, function(x)length(which(x<thresh)) ) / nrow(post.1)
	}
	if ( vs_tag %in% names(MOD.DAT) ) {
		mod.dat.2 <- MOD.DAT[[tag]]
		post.p.2 <- mod.dat.2$D$post[ ,paste("b_",var,sep="") ]
		post.i.2 <- mod.dat.2$D$post[ ,paste( "r_ID[",Samps,",",var,"]",sep="" ) ]
		which_acpa <- grep(paste(Samp.inf$Samps[Samp.inf$ACPA==1],collapse="|"),colnames(post.i.2))
		if ( var=="Intercept" ) {
			post.acpa <- mod.dat.2$D$post[ ,"b_ACPA"]
			post.i.2[,which_acpa] <- post.i.2[,which_acpa] + post.acpa
		}
		if ( var=="DRUG" ) {
			post.acpa <- mod.dat.2$D$post[ ,"b_DRUG:ACPA"]
			post.i.2[,which_acpa] <- post.i.2[,which_acpa] + post.acpa
		}
		post.2 <- post.p.2 + post.i.2
		prob.2 <- apply( post.2, 2, function(x)length(which(x<thresh)) ) / nrow(post.2)
	}else{
		if( vs_tag=="BI" ) {
			# JJ$Ind4$EFF$DRUG[,1]
			ind.tag <- grep("ind",names(JJ),ignore.case=T,value=T)
			prob.2 <- JJ[[ind.tag]]$POST[[paste("b_",var,sep="")]][,paste("lt_",thresh,sep="")]


			# ind.tag <- grep("ind",names(JJ),ignore.case=T,value=T)
			# ind.tag.low <- tolower(ind.tag)
			# # prob.2 <- JJ[[ind.tag]]$POST[[paste("b_",var,sep="")]][ paste(ind.tag.low,samp,sep="_"), paste("lt_",thresh,sep="") ]
			# post.2 <- JJ[[ind.tag]]$POST_RAW[[paste(ind.tag.low,samp,sep="_")]][, paste("b_",var,sep="") ]
			# prob.2 <- lapply( thresh, function(x) length(which(post.2<x))/7200 )

		}
	}

	## Plot Distributions & Scatter of Posterior Probabilities
	 # Plot Parameters
	# COLS <- COLS.list.2[1:5]
	# COLS.a <- adjustcolor( COLS, .7 )
	MAIN <- paste("Posterior Probability:",var,"<",thresh)
	 # Densities
	dens.1 <- density( prob.1 )
	dens.1.grp <- lapply( ARMS,function(x) density(prob.1[Samp.inf$GRP==x]) ) ; names(dens.1.grp) <- ARMS
	dens.2 <- density( prob.2 )
	dens.2.grp <- lapply( ARMS,function(x) density(prob.2[Samp.inf$GRP==x]) ) ; names(dens.2.grp) <- ARMS
	 # Plot it
	png( paste(PathToPlot,"ModPost_",var,"_",thresh,".",tag,".v.",vs_tag,".png",sep=""),height=1600,width=1600,pointsize=32)
	layout( matrix(c(1,4,2,3),ncol=2), widths=c(1,4),height=c(4,1) )
	 # Distribution 1
	par( mar=c(3,4,4,1) )
	plot( dens.1$y[dens.1$x<=1 & dens.1$x>=0], dens.1$x[dens.1$x<=1 & dens.1$x>=0], type="l", col=Cols.m[1],lwd=6, ylim=c(0,1), xlab="",ylab=tag,xaxt="n" )
	lapply( ARMS, function(a)points( .5*dens.1.grp[[a]]$y[dens.1.grp[[a]]$x<=1 & dens.1.grp[[a]]$x>=0], dens.1.grp[[a]]$x[dens.1.grp[[a]]$x<=1 & dens.1.grp[[a]]$x>=0], type="l", col=Cols.a[a],lwd=3,lty=2 ) )
	abline( h=seq(0,1,.2), lty=3,col="grey50",lwd=1 )
	abline( h=seq(0,1,.5),v=0, lty=1,col="black",lwd=1 )
	 # Scatter
	par( mar=c(3,3,4,1) )
	plot( 0,0,type="n", xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main=MAIN )
	abline( h=seq(0,1,.2),v=seq(0,1,.2), lty=3,col="grey50",lwd=1 )
	abline( h=seq(0,1,.5),v=seq(0,1,.5), lty=1,col="black",lwd=1 )
	abline( 0,1, lty=1,col="black",lwd=1 )
	points(prob.1 ~ prob.2, col=Cols.a.7[factor(Samp.inf$GRP)],pch=16 )
	 # Distribution 2
	par( mar=c(4,3,1,1) )
	plot( dens.2$x[dens.2$x<=1 & dens.2$x>=0], dens.2$y[dens.2$x<=1 & dens.2$x>=0], type="l", col=Cols.m[2],lwd=4, xlim=c(0,1), ylab="",xlab=vs_tag,yaxt="n" )
	lapply( ARMS, function(a)points( dens.2.grp[[a]]$x[dens.2.grp[[a]]$x<=1 & dens.2.grp[[a]]$x>=0], .5*dens.2.grp[[a]]$y[dens.2.grp[[a]]$x<=1 & dens.2.grp[[a]]$x>=0], type="l", col=Cols.a[a],lwd=3,lty=2 ) )
	abline( v=seq(0,1,.2), lty=3,col="grey50",lwd=1 )
	abline( v=seq(0,1,.5),h=1, lty=1,col="black",lwd=1 )
	 # Legend
	plot( 0,0,type="n",xlim=c(0,1),ylim=c(0,1),yaxt="n",xaxt="n",xlab="",ylab="")
	legend( 0,1, col=c(Cols.m,Cols.a),legend=c(paste("Mod:",c(tag,vs_tag)),paste("Arm:",ARMS)),pch=c(16,16,NA,NA,NA),lty=c(0,0,1,1,1),lwd=6,box.lwd=0,cex=.7)
	# legend("topleft",col=c(Cols.m,Cols.a),legend=c(paste("Mod:",c(tag,vs_tag)),paste("Arm:",ARMS)),pch=c(16,16,NA,NA,NA),lty=c(0,0,1,1,1),lwd=2 )
	dev.off()

	## Return Posterior Probabilities
	OUT <- data.frame( prob.1, prob.2 )
	colnames(OUT) <- paste( var,thresh,c(tag,vs_tag), sep="." )
	return(OUT)
}

## FCT: Plot Posterior Probabilities of Clinically Meaningful Response
 # Color by EULAR/ACR50 outcomes
PLOT_POSTERIORS.2 <- function( tag, vs_tag, var, thresh, category ) {
	## General Info
	 # Variable Name
	var.tag <- Which_Vars[ var, "TAG" ]
	var.lab <- Which_Vars[ var, "LAB" ]
	if ( tag=="Tn_m10" ) { tag.print <- "m9" }
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Categorical Outcome Variables
	if ( grepl("acr",category,ignore.case=T) ) {
		CATS <- as.character(unique( Samp.inf$ACR50_28 ))
		CAT.tag.28 <- paste("ACR50_28")
		CAT.tag.52 <- paste("ACR50_52")
		Cols.cat <- COLS.list.2[c(6,7)]
	}else{
		CATS <- as.character(unique( Samp.inf$EULAR_28 ))
		CAT.tag.28 <- paste("EULAR_28")
		CAT.tag.52 <- paste("EULAR_52")
		Cols.cat <- COLS.list.2[c(6,3,7)]
	}
	 # Specify Colors
	Cols.cat.7 <- adjustcolor( Cols.cat, .8 )
	Cols.cat.3 <- adjustcolor( Cols.cat, .6 )
	names(Cols.cat) <- names(Cols.cat.7) <- names(Cols.cat.3) <- CATS
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .8 )
	Cols.m <- KEY.2[c(tag,vs_tag),"COL"] # COLS.mods
	Cols.m.7 <- adjustcolor( Cols.m, .8 )

	## Pull out Model Summary
	if ( tag %in% names(MOD.DAT) ) {
		mod.dat.1 <- MOD.DAT[[tag]]
		post.p.1 <- mod.dat.1$D$post[ ,paste("b_",var,sep="") ]
		post.i.1 <- mod.dat.1$D$post[ ,paste( "r_ID[",Samps,",",var,"]",sep="" ) ]
		which_acpa <- grep(paste(Samp.inf$Samps[Samp.inf$ACPA==1],collapse="|"),colnames(post.i.1))
		if ( var=="Intercept" ) {
			post.acpa <- mod.dat.1$D$post[ ,"b_ACPA"]
			post.i.1[,which_acpa] <- post.i.1[,which_acpa] + post.acpa
		}
		if ( var=="DRUG" ) {
			post.acpa <- mod.dat.1$D$post[ ,"b_DRUG:ACPA"]
			post.i.1[,which_acpa] <- post.i.1[,which_acpa] + post.acpa
		}
		post.1 <- post.p.1 + post.i.1
		prob.1 <- apply( post.1, 2, function(x)length(which(x<thresh)) ) / nrow(post.1)
	}
	if ( vs_tag %in% names(MOD.DAT) ) {
		mod.dat.2 <- MOD.DAT[[tag]]
		post.p.2 <- mod.dat.2$D$post[ ,paste("b_",var,sep="") ]
		post.i.2 <- mod.dat.2$D$post[ ,paste( "r_ID[",Samps,",",var,"]",sep="" ) ]
		which_acpa <- grep(paste(Samp.inf$Samps[Samp.inf$ACPA==1],collapse="|"),colnames(post.i.2))
		if ( var=="Intercept" ) {
			post.acpa <- mod.dat.2$D$post[ ,"b_ACPA"]
			post.i.2[,which_acpa] <- post.i.2[,which_acpa] + post.acpa
		}
		if ( var=="DRUG" ) {
			post.acpa <- mod.dat.2$D$post[ ,"b_DRUG:ACPA"]
			post.i.2[,which_acpa] <- post.i.2[,which_acpa] + post.acpa
		}
		post.2 <- post.p.2 + post.i.2
		prob.2 <- apply( post.2, 2, function(x)length(which(x<thresh)) ) / nrow(post.2)
	}else{
		if( vs_tag=="BI" ) {
			# JJ$Ind4$EFF$DRUG[,1]
			ind.tag <- grep("ind",names(JJ),ignore.case=T,value=T)
			prob.2 <- JJ[[ind.tag]]$POST[[paste("b_",var,sep="")]][,paste("lt_",thresh,sep="")]
		}
	}

	## Plot Distributions & Scatter of Posterior Probabilities
	 # Plot Parameters
	MAIN <- paste("Posterior Probability:",var,"<",thresh)
	 # Densities
	dens.1 <- density( prob.1 )
	dens.1.cat.28 <- lapply( CATS,function(x) density(prob.1[Samp.inf[[CAT.tag.28]]==x]) )
	dens.1.cat.52 <- lapply( CATS,function(x) density(prob.1[Samp.inf[[CAT.tag.52]]==x]) )
	names(dens.1.cat.28) <- names(dens.1.cat.52) <- CATS
	dens.2 <- density( prob.2 )
	dens.2.cat.28 <- lapply( CATS,function(x) density(prob.2[Samp.inf[[CAT.tag.28]]==x]) )
	dens.2.cat.52 <- lapply( CATS,function(x) density(prob.2[Samp.inf[[CAT.tag.52]]==x]) )
	names(dens.2.cat.28) <- names(dens.2.cat.52) <- CATS
	 # Plot it
	png( paste(PathToPlot,"ModPost.",category,"_",var,"_",thresh,".",tag,".v.",vs_tag,".png",sep=""),height=1600,width=1600,pointsize=32)
	layout( matrix(c(1,4,2,3),ncol=2), widths=c(1,4),height=c(4,1) )
	 # Distribution 1
	par( mar=c(3,4,4,1) )
	plot( dens.1$y[dens.1$x<=1 & dens.1$x>=0], dens.1$x[dens.1$x<=1 & dens.1$x>=0], type="l", col=Cols.m[1],lwd=6, ylim=c(0,1), xlab="",ylab=tag.print,xaxt="n" )
	lapply( CATS, function(a)points( .5*dens.1.cat.28[[a]]$y[dens.1.cat.28[[a]]$x<=1 & dens.1.cat.28[[a]]$x>=0], dens.1.cat.28[[a]]$x[dens.1.cat.28[[a]]$x<=1 & dens.1.cat.28[[a]]$x>=0], type="l", col=Cols.cat.7[a],lwd=3,lty=2 ) )
	lapply( CATS, function(a)points( .5*dens.1.cat.52[[a]]$y[dens.1.cat.52[[a]]$x<=1 & dens.1.cat.52[[a]]$x>=0], dens.1.cat.52[[a]]$x[dens.1.cat.52[[a]]$x<=1 & dens.1.cat.52[[a]]$x>=0], type="l", col=Cols.cat.7[a],lwd=3,lty=3 ) )
	# lapply( ARMS, function(a)points( .5*dens.1.grp[[a]]$y[dens.1.grp[[a]]$x<=1 & dens.1.grp[[a]]$x>=0], dens.1.grp[[a]]$x[dens.1.grp[[a]]$x<=1 & dens.1.grp[[a]]$x>=0], type="l", col=Cols.a[a],lwd=3,lty=2 ) )
	abline( h=seq(0,1,.2), lty=3,col="grey50",lwd=1 )
	abline( h=seq(0,1,.5),v=0, lty=1,col="black",lwd=1 )
	 # Scatter
	par( mar=c(3,3,4,1) )
	plot( 0,0,type="n", xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main=MAIN )
	abline( h=seq(0,1,.2),v=seq(0,1,.2), lty=3,col="grey50",lwd=1 )
	abline( h=seq(0,1,.5),v=seq(0,1,.5), lty=1,col="black",lwd=1 )
	abline( 0,1, lty=1,col="black",lwd=1 )
	# points(prob.1 ~ prob.2, col=Cols.cat.7[factor(Samp.inf[[CAT.tag.28]])],bg=Cols.cat.3[factor(Samp.inf[[CAT.tag.52]])],pch=21,lwd=5, cex=1.8 )
	points(prob.1 ~ prob.2, col=Cols.cat.7[Samp.inf[[CAT.tag.28]]],bg=Cols.cat.3[Samp.inf[[CAT.tag.52]]],pch=21,lwd=5, cex=1.8 )
	   # (Scatter Legend)
	Legend.1 <- c( paste(category,": ",CATS,sep=""), paste("Week",c(28,52)) )
	Legend.1.pch <- c( rep(15,length(CATS)), 16, 1 )
	Legend.1.col <- c( Cols.cat, rep("grey20",2) )
	Legend.1.cex <- c( rep(2.4,length(CATS)), rep(1.8,2) )
	legend( "topleft", legend=Legend.1, pch=Legend.1.pch, col=Legend.1.col, pt.lwd=5, lty=0, pt.cex=Legend.1.cex,bg="white",ncol=2 )
	 # Distribution 2
	par( mar=c(4,3,1,1) )
	plot( dens.2$x[dens.2$x<=1 & dens.2$x>=0], dens.2$y[dens.2$x<=1 & dens.2$x>=0], type="l", col=Cols.m[2],lwd=4, xlim=c(0,1), ylab="",xlab=vs_tag,yaxt="n" )
	# lapply( CATS, function(a)points( .5*dens.1.cat.28[[a]]$y[dens.1.cat.28[[a]]$x<=1 & dens.1.cat.28[[a]]$x>=0], dens.1.cat.28[[a]]$x[dens.1.cat.28[[a]]$x<=1 & dens.1.cat.28[[a]]$x>=0], type="l", col=Cols.cat[a],lwd=3,lty=2 ) )
	lapply( CATS, function(a)points( dens.2.cat.28[[a]]$x[dens.2.cat.28[[a]]$x<=1 & dens.2.cat.28[[a]]$x>=0], .5*dens.2.cat.28[[a]]$y[dens.2.cat.28[[a]]$x<=1 & dens.2.cat.28[[a]]$x>=0], type="l", col=Cols.cat.7[a],lwd=3,lty=2 ) )
	lapply( CATS, function(a)points( dens.2.cat.52[[a]]$x[dens.2.cat.52[[a]]$x<=1 & dens.2.cat.52[[a]]$x>=0], .5*dens.2.cat.52[[a]]$y[dens.2.cat.52[[a]]$x<=1 & dens.2.cat.52[[a]]$x>=0], type="l", col=Cols.cat.7[a],lwd=3,lty=3 ) )
	abline( v=seq(0,1,.2), lty=3,col="grey50",lwd=1 )
	abline( v=seq(0,1,.5),h=1, lty=1,col="black",lwd=1 )
	 # Legend
	Legend <- c( paste("Mod:",c(tag.print,vs_tag)), paste("Week",c(28,52)) )
	Legend.lty <- c( 1,1, 2,3 )
	Legend.lwd <- c( 6,6, 3,3 )
	Legend.col <- c( Cols.m, rep("grey20",2) )
	plot( 0,0,type="n",xlim=c(0,1),ylim=c(0,1),yaxt="n",xaxt="n",xlab="",ylab="")
	legend( "topleft", legend=Legend, lty=Legend.lty, col=Legend.col, lwd=Legend.lwd, box.lwd=0,cex=.9 )
	dev.off()

	## Return Posterior Probabilities
	OUT <- data.frame( prob.1, prob.2 )
	colnames(OUT) <- paste( var,thresh,c(tag,vs_tag), sep="." )
	return(OUT)
}

## FCT: Plot Posterior Probabilities of Clinically Meaningful Response
PLOT_POST_DENS.1 <- function( tag, vs_tag, var, samp ) {
	## General Info
	 # Variable Name
	var.tag <- Which_Vars[ var, "TAG" ]
	var.lab <- Which_Vars[ var, "LAB" ]
	 # Print Tag
	if ( tag=="Tn_m10" ) { tag.print <- "m9" }else{ tag.print <- tag }
	 # Thresholds
	thresh <- c(0,-.6,-1.2,-2)
	thresh <- c(0,-.6,-1.2)
	thresh <- c(0,-1.2)
	# thresh <- c(0,-.6,-1.2,-1.8)
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .6 )
	Cols.a.2 <- adjustcolor( Cols.a, .2 )
	Cols.m <- KEY.2[c(tag,vs_tag),"COL"] # COLS.mods
	Cols.m.7 <- adjustcolor( Cols.m, .6 )
	Cols.m.2 <- adjustcolor( Cols.m, .2 )

	## Pull out Model Summary
	if ( tag %in% names(MOD.DAT) ) {
		mod.dat.1 <- MOD.DAT[[tag]]
		post.p.1 <- mod.dat.1$D$post[ ,paste("b_",var,sep="") ]
		post.i.1 <- mod.dat.1$D$post[ ,paste( "r_ID[",samp,",",var,"]",sep="" ) ]
		if ( samp %in% Samp.inf$Samps[Samp.inf$ACPA==1] ) {
			if ( var=="Intercept" ) {
				post.acpa <- mod.dat.1$D$post[ ,"b_ACPA"]
				post.i.1 <- post.i.1 + post.acpa
			}
			if ( var=="DRUG" ) {
				post.acpa <- mod.dat.1$D$post[ ,"b_DRUG:ACPA"]
				post.i.1 <- post.i.1 + post.acpa
			}
		}
		post.1 <- post.p.1 + post.i.1
		prob.1 <- lapply( thresh, function(x) length(which(post.1<x))/7200 )
	}
	if ( vs_tag %in% names(MOD.DAT) ) {
		mod.dat.2 <- MOD.DAT[[tag]]
		post.p.2 <- mod.dat.2$D$post[ ,paste("b_",var,sep="") ]
		post.i.2 <- mod.dat.2$D$post[ ,paste( "r_ID[",Samps,",",var,"]",sep="" ) ]
		which_acpa <- grep(paste(Samp.inf$Samps[Samp.inf$ACPA==1],collapse="|"),colnames(post.i.2))
		if ( samp %in% Samp.inf$Samps[Samp.inf$ACPA==1] ) {
			if ( var=="Intercept" ) {
				post.acpa <- mod.dat.2$D$post[ ,"b_ACPA"]
				post.i.2 <- post.i.2 + post.acpa
			}
			if ( var=="DRUG" ) {
				post.acpa <- mod.dat.2$D$post[ ,"b_DRUG:ACPA"]
				post.i.2 <- post.i.2 + post.acpa
			}
		}
		post.2 <- post.p.2 + post.i.2
		prob.2 <- lapply( thresh, function(x) length(which(post.2<x))/7200 )
	}else{
		if( vs_tag=="BI" ) {
			# JJ$Ind4$EFF$DRUG[,1]
			ind.tag <- grep("ind",names(JJ),ignore.case=T,value=T)
			ind.tag.low <- tolower(ind.tag)
			# prob.2 <- JJ[[ind.tag]]$POST[[paste("b_",var,sep="")]][ paste(ind.tag.low,samp,sep="_"), paste("lt_",thresh,sep="") ]
			post.2 <- JJ[[ind.tag]]$POST_RAW[[paste(ind.tag.low,samp,sep="_")]][, paste("b_",var,sep="") ]
			prob.2 <- lapply( thresh, function(x) length(which(post.2<x))/7200 )
		}
	}

	## Plot Posterior Distributions for Patient
	dens.1 <- density( post.1 )
	dens.2 <- density( post.2 )

	# par(mfrow=c(2,1))
	# layout( matrix(1:2,ncol=1), heights=c(2,2) )
	par(mar=c(1,4,4,2))
	XLIM <- range( c(post.1,post.2) )
	YLIM <- c(0,max( c(dens.1$y,dens.2$y) ))
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=var.lab,ylab="Density",main=paste(var.lab,"-",samp),xaxt="n" )
	abline( v=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
	points( dens.1$x, dens.1$y, type="l",lwd=8,col=Cols.m[1] )
	points( dens.2$x, dens.2$y, type="l",lwd=8,col=Cols.m[2] )
	polygon( c(dens.1$x[dens.1$x<=0],rev(dens.1$x[dens.1$x<=0])), c(dens.1$y[dens.1$x<=0],rep(0,sum(dens.1$x<=0))), border=Cols.m.2[1],col=Cols.m.2[1] )
	polygon( c(dens.2$x[dens.2$x<=0],rev(dens.2$x[dens.2$x<=0])), c(dens.2$y[dens.2$x<=0],rep(0,sum(dens.2$x<=0))), border=Cols.m.2[2],col=Cols.m.2[2] )
	abline( h=0,v=thresh,lwd=2 ) # , lty=2:(1+length(thresh)),lwd=2 )
	# polygon( c(dens.1$x[dens.1$x< -1.2],rev(dens.1$x[dens.1$x< -1.2])), c(dens.1$y[dens.1$x< -1.2],rep(0,sum(dens.1$x< -1.2))), border=NA,col=Cols.m.7[1] )
	# polygon( c(dens.2$x[dens.2$x< -1.2],rev(dens.2$x[dens.2$x< -1.2])), c(dens.2$y[dens.2$x< -1.2],rep(0,sum(dens.2$x< -1.2))), border=NA,col=Cols.m.7[2] )
	legend( "topleft", col=Cols.m,lty=1,lwd=6, title="Model",legend=c(tag.print,vs_tag),bg="white" )

	## Plot Posterior Distributions for Patient
	par(mar=c(4,4,1,2))
	dens.1 <- density( post.1 )
	dens.2 <- density( post.2 )
	XLIM <- range( c(post.1,post.2) )
	YLIM <- c(0,1)
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=var.lab,ylab="Pr(X<x)" ) # main="Cummulative Density"
	abline( v=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
	points( dens.1$x, cumsum(dens.1$y)/sum(dens.1$y), type="l",lwd=6,col=Cols.m[1] )
	points( dens.2$x, cumsum(dens.2$y)/sum(dens.2$y), type="l",lwd=6,col=Cols.m[2] )
	points( thresh, unlist(prob.1), pch=16,cex=1.5,col=Cols.m[1] )
	points( thresh, unlist(prob.2), pch=16,cex=1.5,col=Cols.m[2] )
	# polygon( c(dens.1$x[dens.1$x<=0],rev(dens.1$x[dens.1$x<=0])), c((cumsum(dens.1$y)/sum(dens.1$y))[dens.1$x<=0],rep(0,sum(dens.1$x<=0))), border=NA,col=Cols.m.2[1] )
	# polygon( c(dens.2$x[dens.2$x<=0],rev(dens.2$x[dens.2$x<=0])), c((cumsum(dens.2$y)/sum(dens.2$y))[dens.2$x<=0],rep(0,sum(dens.2$x<=0))), border=NA,col=Cols.m.2[2] )
	 # Lines w/ Probabilities
	# abline( v=0 )
	abline( v=thresh,lwd=2 ) # , lty=2:(1+length(thresh)),lwd=2 )
	arrows( rep(-10,length(thresh)), unlist(prob.1), thresh, unlist(prob.1), length=0,lwd=4,lty=2:6,col=Cols.m[1] )
	arrows( rep(-10,length(thresh)), unlist(prob.2), thresh, unlist(prob.2), length=0,lwd=4,lty=2:6,col=Cols.m[2] )
	legend( "bottomright", lty=2:(1+length(thresh)),lwd=4, title="Pr(X<x)",legend=thresh,bg="white" )
	# Write Posterior Probs
	XPOS <- .9
	CEX <- .9
	text( quantile(XLIM,XPOS),quantile(YLIM,.9), "Pr(GOL Eff.) < 0", cex=CEX )
	text( quantile(XLIM,XPOS),quantile(YLIM,.82), round(prob.1[[1]],2), col=Cols.m[1],pos=2 )
	text( quantile(XLIM,XPOS),quantile(YLIM,.82), round(prob.2[[1]],2), col=Cols.m[2],pos=4 )
	text( quantile(XLIM,XPOS),quantile(YLIM,.7), "Pr(GOL Eff.) < -1.2", cex=CEX )
	text( quantile(XLIM,XPOS),quantile(YLIM,.62), round(prob.1[[2]],2), col=Cols.m[1],pos=2 )
	text( quantile(XLIM,XPOS),quantile(YLIM,.62), round(prob.2[[2]],2), col=Cols.m[2],pos=4 )
}

## FCT: Plot Posterior Probabilities of Clinically Meaningful Response
PLOT_POST_DENS.2 <- function( tag, vs_tag, var ) {
	## General Info
	 # Variable Name
	var.tag <- Which_Vars[ var, "TAG" ]
	var.lab <- Which_Vars[ var, "LAB" ]
	 # Print Tag
	if ( tag=="Tn_m10" ) { tag.print <- "m9" }else{ tag.print <- tag }
	 # Thresholds
	thresh <- c(0,-.6,-1.2,-2)
	thresh <- c(0,-.6,-1.2)
	thresh <- c(0,-1.2)
	## General Info
	 # Arm Names
	ARMS <- unique(as.character(Samp.inf$GRP))
	Samps <- Samp.inf$Samps
	N.Samps <- length(Samps)
	 # Specify Colors
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .8 )
	Cols.m <- KEY.2[c(tag,vs_tag),"COL"] # COLS.mods
	Cols.m.7 <- adjustcolor( Cols.m, .8 )

	## Pull out Model Summary
	if ( tag %in% names(MOD.DAT) ) {
		mod.dat.1 <- MOD.DAT[[tag]]
		post.p.1 <- mod.dat.1$D$post[ ,paste("b_",var,sep="") ]
		post.i.1 <- mod.dat.1$D$post[ ,paste( "r_ID[",Samps,",",var,"]",sep="" ) ]
		which_acpa <- grep(paste(Samp.inf$Samps[Samp.inf$ACPA==1],collapse="|"),colnames(post.i.1))
		if ( var=="Intercept" ) {
			post.acpa <- mod.dat.1$D$post[ ,"b_ACPA"]
			post.i.1[,which_acpa] <- post.i.1[,which_acpa] + post.acpa
		}
		if ( var=="DRUG" ) {
			post.acpa <- mod.dat.1$D$post[ ,"b_DRUG:ACPA"]
			post.i.1[,which_acpa] <- post.i.1[,which_acpa] + post.acpa
		}
		post.1 <- post.p.1 + post.i.1
		quant.1 <- apply( post.1, 2, function(x)quantile(x,c(.025,.25,.5,.75,.975)) )
		re.order.1 <- order(quant.1["50%",])
		quant.1 <- quant.1[ ,re.order.1 ]
	}
	if ( vs_tag %in% names(MOD.DAT) ) {
		mod.dat.2 <- MOD.DAT[[tag]]
		post.p.2 <- mod.dat.2$D$post[ ,paste("b_",var,sep="") ]
		post.i.2 <- mod.dat.2$D$post[ ,paste( "r_ID[",Samps,",",var,"]",sep="" ) ]
		which_acpa <- grep(paste(Samp.inf$Samps[Samp.inf$ACPA==1],collapse="|"),colnames(post.i.2))
		if ( var=="Intercept" ) {
			post.acpa <- mod.dat.2$D$post[ ,"b_ACPA"]
			post.i.2[,which_acpa] <- post.i.2[,which_acpa] + post.acpa
		}
		if ( var=="DRUG" ) {
			post.acpa <- mod.dat.2$D$post[ ,"b_DRUG:ACPA"]
			post.i.2[,which_acpa] <- post.i.2[,which_acpa] + post.acpa
		}
		post.2 <- post.p.2 + post.i.2
		prob.2 <- apply( post.2, 2, function(x)length(which(x<thresh)) ) / nrow(post.2)
	}else{
		if( vs_tag=="BI" ) {
			# JJ$Ind4$EFF$DRUG[,1]
			# ind.tag <- grep("ind",names(JJ),ignore.case=T,value=T)
			# prob.2 <- JJ[[ind.tag]]$POST[[paste("b_",var,sep="")]][,paste("lt_",thresh,sep="")]


			ind.tag <- grep("ind",names(JJ),ignore.case=T,value=T)
			ind.tag.low <- tolower(ind.tag)
			# prob.2 <- JJ[[ind.tag]]$POST[[paste("b_",var,sep="")]][ paste(ind.tag.low,samp,sep="_"), paste("lt_",thresh,sep="") ]
			post.2 <- lapply( Samps, function(samp)JJ[[ind.tag]]$POST_RAW[[paste(ind.tag.low,samp,sep="_")]][, paste("b_",var,sep="") ])
			# post.2 <- JJ[[ind.tag]]$POST_RAW[[paste(ind.tag.low,Samps,sep="_")]][, paste("b_",var,sep="") ]
			# prob.2 <- lapply( thresh, function(x) length(which(post.2<x))/7200 )

			quant.2 <- Reduce( cbind, lapply( post.2, function(x)quantile(x,c(.025,.25,.5,.75,.975)) ) )
			colnames(quant.2) <- Samps
			re.order.2 <- order(quant.2["50%",])
			quant.2 <- quant.2[ ,re.order.2 ]

		}
	}

	## Plot Posterior Distributions for Patient
	XLIM <- c(5,N.Samps-5) # c(0,N.Samps+1)
	YLIM <- range( quant.1 )
	png( paste(PathToPlot,"ModPost_Quants_",var,".",tag,".v.",vs_tag,".png",sep=""),height=2000,width=3500,pointsize=32)
	par(mfrow=c(2,1))
	 # Model 1
	par(mar=c(4,4,4,1))
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Patient",ylab=paste(var.lab,"- DAS"),main=paste("Posterior Distributions:",var.lab,"-",tag),xaxt="n" )
	abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
	# abline( h=thresh )
	abline( h=c(0,5) )
	arrows( 1:N.Samps, quant.1["2.5%",], 1:N.Samps, quant.1["97.5%",], col=Cols.a[Samp.inf$GRP[re.order.1]], lwd=2, length=0 )
	arrows( 1:N.Samps, quant.1["25%",], 1:N.Samps, quant.1["75%",], col=Cols.a[Samp.inf$GRP[re.order.1]], lwd=4, length=0 )
	points( 1:N.Samps, quant.1["50%",], pch=21, col="black",bg=Cols.m[1],cex=.5 )
	 # Legend
	arrows( 0, quantile(YLIM,.7), 0, quantile(YLIM,.9), col="grey20",lwd=2, length=0 )
	arrows( 0, quantile(YLIM,.75), 0, quantile(YLIM,.85), col="grey20",lwd=4, length=0 )
	points( 0, quantile(YLIM,.8), pch=21,col="black",bg="grey70",cex=0.5 )
	LEG.labs <- c(.88,.84,.8)
	text( 5,quantile(YLIM,LEG.labs), c("95% CI","Quantiles","Median"),pos=4,cex=.9 )
	arrows( 7,quantile(YLIM,LEG.labs),1,quantile(YLIM,LEG.labs) )
	 # Model 2
	YLIM <- range( quant.2 )
	par(mar=c(4,4,4,1))
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Patient",ylab=paste(var.lab,"- DAS"),main=paste("Posterior Distributions:",var.lab,"-",vs_tag),xaxt="n" )
	abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
	# abline( h=thresh )
	abline( h=c(0,5) )
	arrows( 1:N.Samps, quant.2["2.5%",], 1:N.Samps, quant.2["97.5%",], col=Cols.a[Samp.inf$GRP[re.order.2]], lwd=2, length=0 )
	arrows( 1:N.Samps, quant.2["25%",], 1:N.Samps, quant.2["75%",], col=Cols.a[Samp.inf$GRP[re.order.2]], lwd=4, length=0 )
	points( 1:N.Samps, quant.2["50%",], pch=21, col="black",bg=Cols.m[2],cex=.5 )
	dev.off()

	## Compile & Return Quantiles
	OUT <- list( Q1=quant.1, Q2=quant.2 )
	return(OUT)
}

## FCT: Plot Random Effects (Individual Estimates)
PLOT_PRED <- function( tag, vs_tag, samp ) {

	## General Info
	 # Variables to Print
	which_paste <- c("Intercept","DRUG","TRT")
	label_paste <- c("BL","GOL","TRT")
	if ( tag=="Tn_m10" ) { tag.print <- "m9" }else{ tag.print <- tag }
	 # Specify Colors
	COL.tru <- adjustcolor("black",.8)
	Cols.m <- KEY.2[c(tag,vs_tag),"COL"] # COLS.mods
	Cols.m.7 <- adjustcolor( Cols.m, .8 )
	Cols.m.2 <- adjustcolor( Cols.m, .25 )

	## Pull out Actual Data
	TAB.temp <- TAB[ !is.na(TAB$DAS), ]
	which_rows <- which( TAB.temp$ID==samp )
	tru.dat <- TAB.temp[ which_rows, ]
	tru.pch <- paste( tru.dat$DRUG, tru.dat$TRT, sep="" )
	PCHS <- c(4,1,16) ; names(PCHS) <- c("00","01","11")

	## Pull out Model Summary
	if ( tag %in% names(MOD.DAT) ) {
		mod.dat.1 <- MOD.DAT[[tag]]
		mod.pred.1 <- mod.dat.1$Pred[ which_rows, ]
		coef.1 <- mod.dat.1$Coef[ samp, ]
		coef.1.print <- coef.1[which_paste]
		coef.1.print["Intercept"] <- coef.1["Intercept"] + coef.1["ACPA"]*tru.dat$ACPA[1]
		coef.1.print["DRUG"] <- coef.1["DRUG"] + coef.1["DRUG:ACPA"]*tru.dat$ACPA[1]
	}
	if ( vs_tag %in% names(MOD.DAT) ) {
		mod.dat.2 <- MOD.DAT[[vs_tag]]
		mod.pred.2 <- mod.dat.2$Pred[ which_rows, ]
		coef.2 <- mod.dat.2$Coef[ samp, ]
		coef.2.print <- coef.2[which_paste]
		coef.2.print["Intercept"] <- coef.2["Intercept"] + coef.2["ACPA"]*tru.dat$ACPA[1]
		coef.2.print["DRUG"] <- coef.2["DRUG"] + coef.2["DRUG:ACPA"]*tru.dat$ACPA[1]
	}else{
		if( vs_tag=="BI" ) {
			ind.tag <- grep("Ind",names(JJ),value=T)
			ind.tag.samp <- paste( tolower(ind.tag),samp,sep="_" )
			coef.2 <- JJ[[ind.tag]]$EST[ind.tag.samp,]
			coef.2.print <- coef.2
			mod.pred.2 <- JJ[[ind.tag]]$PRED[[ind.tag.samp]]
		}
	}

	## Plot Predicted vs Real DAS over Time
	XLIM <- c(0,100)
	YLIM <- c(1,10)
	MAIN <- paste(samp,"-",tag.print,"vs",vs_tag)
	par(mar=c(4,4,4,2))
	plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Week",ylab="DAS",main=MAIN)
	# abline( v=WKS,h=0:10,lty=3,col="grey50",lwd=1 )
	SCRAP <- lapply( Samps, function(x)points(DAS~WK,TAB.temp,subset=ID==x,type="l",col="grey90",lwd=1) )
	abline( v=seq(0,100,10),h=0:10,lty=3,col="grey50",lwd=1 )
	polygon( c(tru.dat$WK,rev(tru.dat$WK)),c(mod.pred.2[,"2.5%ile"],rev(mod.pred.2[,"97.5%ile"])), col=Cols.m.2[2],border=Cols.m.7[2] )
	polygon( c(tru.dat$WK,rev(tru.dat$WK)),c(mod.pred.1[,"2.5%ile"],rev(mod.pred.1[,"97.5%ile"])), col=Cols.m.2[1],border=Cols.m.7[1] )
	points( DAS ~ WK, data=tru.dat, pch=PCHS[tru.pch],type="p",col="black",lwd=4,cex=1.5 )
	points( mod.pred.2[,"Estimate"] ~ tru.dat$WK, pch=16,type="l",col=Cols.m.7[2],lwd=8,lty=1 )
	points( mod.pred.1[,"Estimate"] ~ tru.dat$WK, pch=16,type="l",col=Cols.m.7[1],lwd=8,lty=1 )
	 # Legend
	legend("topright",col=Cols.m,legend=c(tag.print,vs_tag),lwd=8,lty=1,ncol=1, bg="white", title="Model" )
	legend("topleft",col=COL.tru,legend=c("Baseline","PBO","GOL"),pch=PCHS,ncol=1,pt.lwd=4,pt.cex=1.5, bg="white", title="Visit" )
	 # Arrows: Effect Sizes
	DAS.bl <- c(coef.1.print["Intercept"],coef.2.print["Intercept"])
	DAS.post.trt <- c(coef.1.print["Intercept"],coef.2.print["Intercept"])+c(coef.1.print["TRT"],coef.2.print["TRT"])
	DAS.post.drug <- c(coef.1.print["Intercept"],coef.2.print["Intercept"])+c(coef.1.print["TRT"],coef.2.print["TRT"])+c(coef.1.print["DRUG"],coef.2.print["DRUG"])
	points( c(98,100), DAS.bl, col=Cols.m, cex=1.2,pch="-" )
	arrows( c(98,100), DAS.bl, c(98,100), DAS.post.trt, col=Cols.m,lwd=4,angle=25,length=.12 )
	arrows( c(98,100), DAS.bl, c(98,100), DAS.post.drug, col=Cols.m,lwd=4,angle=25,length=.24 )
	 # Text: Effect Sizes
	paste.eff.1 <- paste(label_paste,round(coef.1.print[which_paste],3),sep="=",collapse="\n")
	paste.eff.2 <- paste(label_paste,round(coef.2.print[which_paste],3),sep="=",collapse="\n")
	text( 16, 9.5, label=paste.eff.1, col=Cols.m[1], pos=4 )
	text( 32, 9.5, label=paste.eff.2, col=Cols.m[2], pos=4 )
}

## FCT: Violin Plot & Anova of Random Effects across Trial Arm
PLOT_ARM_ANOVA <- function( tag, vs_tag, var ) {

	## General Info
	 # Variable Name
	var.tag <- Which_Vars[ var, "TAG" ]
	var.lab <- Which_Vars[ var, "LAB" ]
	 # Arm Names
	ARMS <- unique(as.character(FUL$GRP))
	 # Specify Colors
	Cols.m <- KEY.2[ match( c(tag,vs_tag), KEY.2$TAG.c ), "COL" ]
	Cols.m.7 <- adjustcolor( Cols.m, .7 )
	Cols.a <- COLS.arm
	Cols.a.7 <- adjustcolor( Cols.a, .7 )

	## Specify Model Info (tag)
	tag.goal <- KEY.2[tag,"TAG.g"]
	tag.mod <- KEY.2[tag,"TAG.m"]
	model <- JJ[[tag.goal]][[tag.mod]]
	mod.dat <- MOD.DAT[[tag]]
	 # Random Effects
	if ( var=="Intercept" ) { tag.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
	if ( var=="DRUG" ) { tag.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
	if ( !(var%in%c("DRUG","Intercept")) ) { tag.var <- mod.dat$Coef[,var] }
	 # Fixed Effect & Variance
	tag.fx <- mean( tag.var )
	tag.fx.sd <- sd( tag.var )

	## Specify Model Info (vs_tag)
	if ( vs_tag %in% c("MN","BI") ) {
		vs.var <- COEF[[var.tag]][,vs_tag]
	}else{
		vs.goal <- KEY.2[vs_tag,"TAG.g"]
		vs.mod <- KEY.2[vs_tag,"TAG.m"]
		model <- JJ[[vs.goal]][[vs.mod]]
		mod.dat <- MOD.DAT[[tag]]
		 # Random Effects
		if ( var=="Intercept" ) { vs.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"ACPA"]*Samp.inf$ACPA }
		if ( var=="DRUG" ) { vs.var <- mod.dat$Coef[,var] + mod.dat$Coef[,"DRUG:ACPA"]*Samp.inf$ACPA }
		if ( !(var%in%c("DRUG","Intercept")) ) { vs.var <- mod.dat$Coef[,var] }
	}
	 # Fixed Effect & Variance
	vs.fx <- mean( vs.var )
	vs.fx.sd <- sd( vs.var )

	## Merge Estimates from 2 Models
	mg.var <- data.frame( ID=Samp.inf$Samp, tag.var, vs.var, GRP=Samp.inf$GRP )

	## ANOVA across Arms
	tag.anv <- anova(lm( tag.var ~ GRP, data=mg.var ))
	vs.anv <- anova(lm( vs.var ~ GRP, data=mg.var ))
	tag.p <- tag.anv["GRP","Pr(>F)"]
	vs.p <- vs.anv["GRP","Pr(>F)"]

	XLIM <- c(.5,6.5)
	YLIM <- extendrange( mg.var[,c("tag.var","vs.var")], f=.2 )
	YBAR <- extendrange( mg.var[,c("tag.var","vs.var")], f=.05 )[1]
	if ( var=="WK" ) { BRK <- .05 }else{ BRK <- .5 }
	YLIN <- seq( floor(YLIM[1]), ceiling(YLIM[2]), BRK )
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, ylab=var.lab, xlab="Model/Arm",main=paste(var.lab,"by Clinical Trial Arm"), xaxt="n" )
	axis( 1, at=c(2,5),label=c(tag,vs_tag) )
	abline( h=YLIN,lty=3,col="grey50",lwd=1 )
	abline( h=0 )
	for ( i in 1:3 ) { vioplot(mg.var$tag.var[mg.var$GRP==ARMS[i]], col=Cols.a.7[i], at=i, border=NA, add=T ) }
	for ( i in 1:3 ) { vioplot(mg.var$vs.var[mg.var$GRP==ARMS[i]], col=Cols.a.7[i], at=i+3, border=NA, add=T ) }
	arrows( c(.75,3.75),rep(YBAR,2),c(3.25,6.25),rep(YBAR,2), lwd=4,col=Cols.m,length=0 )
	text( c(2,5), rep(YBAR,2), paste("p =",formatC(c(tag.p,vs.p),2,format="e")), col=Cols.m, pos=1 )
	legend( "topleft", fill=c(Cols.m,Cols.a.7),legend=c(tag,vs_tag,ARMS), bg="white",ncol=5,border=NA )
}

## FCT: Plot Individual Estimates across Models
Scatter_Mod_Pairs <- function( vars, mods ) {
	N.mods <- length( mods )
	par( mfrow=c(N.mods,N.mods) )
	# Plot Colors
	COLS.temp <- COLS.list.2[c(2,3,6)] ; names(COLS.temp) <- c("G","P","PE")
	COLS.temp <- COLS.list.2[3:5] ; names(COLS.temp) <- c("G","P","PE")
	for ( r in 1:N.mods ) {
		for ( c in 1:N.mods ) {
			if ( c==r ) {
				plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
				text( 0,0, mods[c],cex=3 )
				next
			}
			# Set Variable & Parameters
			if ( c > r ) { var <- vars[1] }else{ var <- vars[2] }
			main.key <- data.frame( Var=c("Int","Drug","Trt"),Main=c("Baseline DAS","Drug Effect","Treatment Effect") )
			main <- main.key$Main[ main.key$Var==var ]
			lim <- range( COEF[[var]][,mods] )
			brks <- seq( floor(lim[1]), ceiling(lim[2]), 1 )
			## Create Plot
			plot( 0,0,type="n", xlim=lim,ylim=lim, xlab=mods[c],ylab=mods[r],main=main )
			abline(h=brks,v=brks,lty=3,col="grey50",lwd=1 )
			abline(0,1) ; abline(h=0,v=0)
			points( COEF[[var]][,mods[c]], COEF[[var]][,mods[r]], pch=16,col=adjustcolor(COLS.temp[Samp.inf$GRP],.5) )
		}
	}
}

## FCT: Additional Rough Drafts of Plots for Random Effects
OTHER_RAND_PLOTS <- function( dwai ) {
		## Boxplot of Individual Posterior Distributions
		print("Plotting #R2 - Boxplot")
		# png( paste(PathToPlot,"ModSumm_",tag,".R2-BoxRand.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
		# par(mfrow=c(n.covs.r,1))
		# par(mar=c(6,4,4,1))
		# for ( z in 1:n.covs.r ) {
		# 	z.name <- m.covs.r[z]
		# 	z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
		# 	z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
		# 	z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
		# 	if ( z>0 ) { z.ord <- order( colMeans(z.temp) ) }
		# 	YLIM <- range( z.temp )
		# 	boxplot( z.temp[,z.ord], col=COLS.eff["R"],main=paste("Random",z.name,"by Patient"),ylim=YLIM,ylab=z.name,names=m.samps[z.ord],las=2,pch=16 )
		# 	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
		# 	abline( h=0, lty=1,col="black",lwd=1 )
		# 	boxplot( z.temp[,z.ord], col=COLS.eff["R"],names=m.samps[z.ord],las=2,pch=16,add=T )
		# 	if ( any(z.name %in% c("DRUG","PLAC","TRT")) ) {
		# 		abline(h=c(0,-1),lwd=2,col=COLS.list.2[6:7] )
		# 		text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.5 )
		# 		text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[7],cex=.5 )
		# 	}else{
		# 		legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
		# 	}			
		# }
		# dev.off()
		 # Plot Confidence Intervals of Individual Posterior Distributions (quicker)
		png( paste(PathToPlot,"ModSumm_",tag,".R2-BoxRand.2.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
		par(mfrow=c(n.covs.r,1))
		par(mar=c(6,4,4,1))
		for ( z in 1:n.covs.r ) {
			# Pull Data
			z.name <- m.covs.r[z]
			z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
			z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
			z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
			# z.samps.temp <- gsub("r_ID[","", sapply(strsplit(colnames(z.temp),","),"[",1) ,fixed=T)
			z.grp <- FUL$GRP[ match(m.samps,as.character(FUL$ID)) ]
			which_quants <- c(.025,.25,.5,.75,.975)
			z.quants <- apply( z.temp, 2, function(x)quantile(x,which_quants) )
			if ( z>0 ) { z.ord <- order( colMeans(z.temp) ) }
			# Parameters
			XLIM <- c( 1,ncol(z.temp) )
			YLIM <- range( z.temp )
			# Build Plot
			plot( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp[z.ord]], main=paste("Random",z.name,"by Patient"),xlab="",ylab=z.name,xaxt="n",xlim=XLIM,ylim=YLIM )
			axis( 1, at=1:ncol(z.temp),label=m.samps[z.ord],las=2,cex=.6 )
			abline( h=-10:10,v=1:ncol(z.temp), lty=3,col="grey50",lwd=1 )
			abline( h=0, lty=1,col="black",lwd=1 )
			# Text/Legend
			if ( any(z.name %in% c("DRUG","PLAC","TRT")) ) {
				abline(h=c(0,-1),lwd=2,col=COLS.list.2[6:7] )
				Post.eff0 <- round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1)
				Post.eff1 <- round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1)
				text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=Post.eff0,col=COLS.list.2[6],cex=.7,srt=90 )
				text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=Post.eff1,col=COLS.list.2[7],cex=.7,srt=90 )
			}else{
				# legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
				legend( "topleft",pch=c(rep(16,3),16,1),cex=1.5,col=c(COLS.list.2[3:5],rep("black",2)),legend=c("Arm: G","Arm: P","Arm: PE","Median","CI: 95%"),lwd=3 )
			}
			# Points/Arrows
			arrows( 1:ncol(z.temp),z.quants["25%",z.ord],1:ncol(z.temp),z.quants["75%",z.ord],col=COLS[1],lwd=5,code=3,angle=90,length=.1 )
			points( rep(1:ncol(z.temp),2),c(z.quants["2.5%",z.ord],z.quants["97.5%",z.ord]),pch=1,lwd=3,col=COLS.list.2[3:5][z.grp[z.ord]],cex=1.5 )
			points( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp[z.ord]],cex=1.5 )
		}
		dev.off()

		# ## Violin Plot: Prior vs Posterior (Compiled)
		# print("Plotting #R2 - Violin Plot")
		# png( paste(PathToPlot,"ModSumm_",tag,".R2b-ViolRand.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
		# par(mfrow=c(n.covs.r,1))
		# par(mar=c(6,4,4,1))
		# for ( z in 1:n.covs.r ) {
		# 	z.name <- m.covs.r[z]
		# 	z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
		# 	z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
		# 	z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
		# 	if ( z==1 ) { z.ord <- order( colMeans(z.temp) ) }
		# 	XLIM <- c(1,ncol(z.temp))
		# 	YLIM <- extendrange(z.temp,f=.2)
		# 	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="",main="Posterior Population Estimate (by Patient)",xaxt="n" )
		# 	axis( 1, at=1:ncol(z.temp), label=m.samps[z.ord], las=2)
		# 	abline(h=-10:10,lty=3,col="grey50",lwd=1 )
		# 	abline(h=0)
		# 	# Get Prior
		# 	if ( z.name=="Intercept" ) { temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
		# 	}else{ temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==z.name,"prior"], fixed=T),fixed=T),"," )[[1]]) }
		# 	names(temp.priors) <- c("Mean","SD")
		# 	temp.prior.distr <- rnorm(10000,temp.priors["Mean"],temp.priors["SD"])
		# 	# Plot Prior/Posterior
		# 	SCRAP <- lapply( 1:ncol(z.temp), function(x)vioplot(temp.prior.distr,at=x, col=adjustcolor("black",alpha=.2),border=NA,add=T ))
		# 	SCRAP <- lapply( 1:ncol(z.temp), function(x)vioplot(z.temp[,z.ord[x]],at=x, col=COLS.eff["R"],add=T ))
		# 	if ( any(z.name %in% c("DRUG","PLAC")) ) {
		# 		abline(h=c(0,-1),lwd=2,col=COLS.list.2[5:6] )
		# 		text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[5],cex=.5 )
		# 		text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.5 )
		# 	}else{
		# 		legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
		# 	}			
		# }
		# dev.off()

		## Plot Posterior Probabilities for Random Effects
		print("Plotting #R3 - Heatmaps of Posteriors")
		 # Plotting Parameters
		COLS.heat <- colorRampPalette(c(COLS.list.heat[1:3],"white",COLS.list.heat[4:6]))(100)
		BRKS.heat <- seq(0,1,length.out=101)
		n.iter.tot <- m.chain * ( m.iter - m.warm )
		post.prob.brks <- list( Intercept=seq(3.5,7.5,.125),
			DRUG=seq(-3,0,.125),
			PLAC=seq(-2,1,.125),
			TRT=seq(-2,1,.125) )
		 # Calculate/Plot Posteriors at Various Thresholds
		post.prob <- list()
		for ( z in 1:n.covs.r ) {
			z.name <- m.covs.r[z]
			z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
			z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
			z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
			z.grp <- FUL$GRP[ match(m.samps,as.character(FUL$ID)) ]
			# Calculate Posteriors
			post.prob[[z.name]] <- Reduce( cbind, lapply( post.prob.brks[[z.name]],function(x)apply(z.temp,2,function(y)length(which(y<x))) ) ) / n.iter.tot
			colnames(post.prob[[z.name]]) <- paste("LT",post.prob.brks[[z.name]],sep="_")
			rownames(post.prob[[z.name]]) <- colnames(z.temp)
			if ( z.name=="Intercept" ) { post.prob[[z.name]] <- 1 - post.prob[[z.name]] ; colnames(post.prob[[z.name]]) <- gsub("LT","GT",colnames(post.prob[[z.name]])) }
			# Plot Heatmap
			png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.",z.name,".png",sep=""),height=1200,width=2000+10*n.samps,pointsize=28)
			heatmap.2( t(post.prob[[z.name]]), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5+n.samps/100),main=z.name,ColSideColors=COLS.list.2[3:5][z.grp] )
			dev.off()
		}
		 # Heatmap of All Random Effect Posteriors
		z.grp <- FUL$GRP[ match(m.samps,as.character(FUL$ID)) ]
		post.prob.temp <- Reduce( cbind, post.prob )
		png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.All.png",sep=""),height=800+400*n.covs.r,width=2000+10*n.samps,pointsize=28)
		heatmap.2( t(post.prob.temp), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5+n.samps/100),rowsep=cumsum(lapply(post.prob,ncol)),ColSideColors=COLS.list.2[3:5][z.grp] )
		dev.off()

		# png( paste(PathToPlot,"ModSumm_",tag,".R3-PostScatter.png",sep=""),height=800+400*n.covs.r,width=2000+10*n.samps,pointsize=28)
		# par(mfrow=c(n.covs.r,1))
		# for ( z in 1:n.covs.r ) {
		# 	z.name <- m.covs.r[z]
		# 	temp.split <- as.numeric(sapply(strsplit(colnames(post.prob[[z.name]]),"_"),"[",2))
		# 	which_cols <- which( temp.split %in% unique(round(2*temp.split)/2) )
		# 	z.ord <- order( post.prob[[z.name]][,which.max(apply(post.prob[[z.name]],2,sd))] )
		# 	COLS.temp <- colorRampPalette(COLS.list.ord)(length(which_cols))
		# 	plot( 0,0,type="n", xlim=c(1,n.samps),ylim=c(0,1) )
		# 	lapply( 1:length(which_cols),function(x)points(1:n.samps,post.prob[[z.name]][z.ord,which_cols[x]],col=adjustcolor(COLS.temp[x],alpha=.7),pch=16,type="o",lwd=2) )
		# }
		# dev.off()
		
		## Plot Distribution (across patients) of Standard Deviations (across iteration) of Random Effects
		# TEST <- lapply(m.covs.r,function(r)unlist(lapply( m.samps, function(s) sd(d.post[,paste("r_ID[",s,",",r,"]",sep="")]) )) )
		# names(TEST) <- m.covs.r

		## Plot DRUG vs Intercept for each Simulation/Patient
		# if ( all(c("DRUG","Intercept")%in%m.covs.r) ) {
		# 	print("Plotting #R4 - Drug v Int per Patient")
		# 	x.name <- "Intercept"
		# 	x.col.1 <- d.post.f.which[ d.post.f.which==paste("b",x.name,sep="_") ] # grep(x.name,d.post.f.which,value=T)
		# 	x.col.2 <- grep( paste(",",x.name,"]",sep=""),d.post.r2.which,value=T)
		# 	x.temp <- d.post[,x.col.1] + d.post[,x.col.2]
		# 	y.name <- "DRUG"
		# 	y.col.1 <- d.post.f.which[ d.post.f.which==paste("b",y.name,sep="_") ] # grep(y.name,d.post.f.which,value=T)
		# 	y.col.2 <- grep( paste(",",y.name,"]",sep=""),d.post.r2.which,value=T)
		# 	y.temp <- d.post[,y.col.1] + d.post[,y.col.2]
		# 	XLIM <- range(x.temp)
		# 	YLIM <- range(y.temp)
		# 	COLS.ind <- colorRampPalette(COLS.list)(ncol(x.temp))
		# 	COLS.ind.2 <- adjustcolor(COLS.ind,alpha=.01)
		# 	png( paste(PathToPlot,"ModSumm_",tag,".R4-Sims.DRvINT.png",sep=""),height=1000,width=1600,pointsize=30)
		# 	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="Drug",main="Individual Simulation Estimates")
		# 	abline(h=-10:10,v=-10:10,lty=3,col="grey50",lwd=1 )
		# 	abline(h=0)
		# 	SCRAP <- lapply( 1:ncol(x.temp),function(i)points(x.temp[,i],y.temp[,i],col=COLS.ind.2[i],pch=16) )
		# 	mod.temp <- lm(unlist(y.temp)~unlist(x.temp))
		# 	abline( mod.temp )
		# 	points( colMeans(x.temp),colMeans(y.temp),pch="+",col=adjustcolor(COLS.ind,alpha=.7),cex=2,lwd=2 )
		# 	dev.off()
		# }

		## Plot a Few Individual Patients' Profiles
		if ( all(c("DRUG","Intercept")%in%m.covs.r) & any(c("PLAC","TRT")%in%m.covs.r) ) {
		# if ( all(c("DRUG","PLAC","Intercept")%in%m.covs.r) ) {
			print("Plotting Individual Patients")
			## PLAC or TRT?
			PBO <- grep("PLAC|TRT",m.covs.r,value=T )
			## Which Patients?
			 # How Many?
			plot.rows <- 4
			plot.cols <- 5
			z.samps.n <- plot.rows * plot.cols
			 # Specify Samples
			# z.samps <- sample(m.samps, z.samps.n )
			z.samps.temp.ord.int <- order( colMeans(d.post[,grep( paste(",Intercept]",sep=""),d.post.r2.which,value=T)]) )
			z.samps.temp.ord.dr <- order( colMeans(d.post[,grep( paste(",DRUG]",sep=""),d.post.r2.which,value=T)]) )
			z.samps.temp <- m.samps[ c( head(z.samps.temp.ord.int,plot.rows),tail(z.samps.temp.ord.int,plot.rows), head(z.samps.temp.ord.dr,plot.rows),tail(z.samps.temp.ord.dr,plot.rows) ) ]
			z.samps <- sample( setdiff(m.samps,z.samps.temp), z.samps.n-length(z.samps.temp) )
			z.samps <- c( head(z.samps.temp,2*plot.rows), z.samps, tail(z.samps.temp,2*plot.rows) )
			# Pull out Sample Data
			z.data <- model$data[model$data[,r.grp]%in%z.samps,]
			z.temp <- predict( model, newdata=z.data )
			z.pred <- cbind( z.data, z.temp )
			z.ranef <- t(f.eff[m.covs.r,"Estimate"] + t(r.eff[z.samps,]))
			if ( PBO=="TRT" ) { z.ranef[,"DRUG"] <- z.ranef[,"DRUG"] + z.ranef[,"TRT"]}

			# print("Plotting #R5 - Predicted Values")
			## FCT: Plot Patient Profile
			sample <- z.samps[1]
			PLOT_IND <- function( sample, which_plots ) {
				## Plot 1: Real vs Predicted Values
				if ( 1 %in% which_plots ) {
					COLS.tru <- COLS.list.2[c(6,2,1)]
					COLS.conf <- COLS.list.2[4]
					sub.tab <- z.pred[ z.pred[,r.grp]==sample, ]
					plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main=sample )
					abline( h=0:10,lty=3,col="grey50",lwd=1)
					points( Estimate ~ WK, data=sub.tab,type="l",lwd=4 )
					polygon( c(sub.tab$WK,rev(sub.tab$WK)), c(sub.tab[,"2.5%ile"],rev(sub.tab[,"97.5%ile"])), col=adjustcolor(COLS.conf,alpha=.2),border=NA )
					points( DAS ~ WK, data=sub.tab,type="p",lwd=3,cex=1.5,pch=c(4,16)[factor(DRUG)],col=COLS.tru[1+rowSums(sub.tab[,c(PBO,"DRUG")])] )
					# pre.post <- cumsum(unlist( FUL[FUL$ID==sample,c("DAS_BL_MN","DEL_MNe_MN")] ))
					# arrows( c(0,24),pre.post,c(24,100),pre.post, col=COLS.pred[5],lwd=4,lty=2,length=0 )
					pre.post.2 <- z.ranef[sample,c(1,3,2)] + c(0,rep(z.ranef[sample,1],2))
					arrows( c(-5,2,24),pre.post.2,c(100,24,100),pre.post.2, col=COLS.tru,lwd=6,lty=2,length=0 )
					text( rep(100,3),c(9,8.5,8), paste(colnames(z.ranef),round(z.ranef[sample,],2),sep="=")[c(1,3,2)], pos=2,col=COLS.tru )
				}

				## Plot 2: Scatter of Individual Monte Carlo Estimates
				if ( 2 %in% which_plots ) {
					x.name <- "Intercept"
					x.col.1 <- d.post.f.which[ d.post.f.which==paste("b",x.name,sep="_") ] # grep(x.name,d.post.f.which,value=T)
					x.col.2 <- grep( paste(sample,",",x.name,"]",sep=""),d.post.r2.which,value=T)
					x.temp <- d.post[,x.col.1] + d.post[,x.col.2]
					y.name <- "DRUG"
					y.col.1 <- d.post.f.which[ d.post.f.which==paste("b",y.name,sep="_") ] # grep(y.name,d.post.f.which,value=T)
					y.col.2 <- grep( paste(sample,",",y.name,"]",sep=""),d.post.r2.which,value=T)
					y.temp <- d.post[,y.col.1] + d.post[,y.col.2]
					w.name <- PBO
					w.col.1 <- d.post.f.which[ d.post.f.which==paste("b",w.name,sep="_") ] # grep(w.name,d.post.f.which,value=T)
					w.col.2 <- grep( paste(sample,",",w.name,"]",sep=""),d.post.r2.which,value=T)
					w.temp <- d.post[,w.col.1] + d.post[,w.col.2]
					if ( PBO=="TRT" ) { y.temp <- w.temp + y.temp }
					XLIM <- c(1,8) # range(x.temp)
					YLIM <- c(-3,2) # range(y.temp)
					# Plot it
					plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="Drug",main="Individual Simulation Estimates")
					abline(h=-10:10,v=-10:10,lty=3,col="grey50",lwd=1 )
					abline(h=0,v=5)
					SCRAP <- points(x.temp,y.temp,col=adjustcolor(COLS.tru[3],.03),pch=16)
					SCRAP <- points(x.temp,w.temp,col=adjustcolor(COLS.tru[2],.03),pch=16)
					points( mean(x.temp),mean(y.temp),pch="+",col=COLS.tru[3],cex=2,lwd=2 )
					points( mean(x.temp),mean(w.temp),pch="+",col=COLS.tru[2],cex=2,lwd=2 )
					legend( "topleft",col=COLS.tru[2:3],pch=16,legend=c("TRT","DRUG+TRT"),title="Effect Sizes" )
				}					
			}

			for ( samp in z.samps ) {
			# for ( samp in z.samps[seq(1,20,4)] ) {
				png( paste(PathToPlot,"ModSumm_",tag,".Ri.",samp,"-Pred.png",sep=""),height=1000,width=2000,pointsize=30)
				par(mfrow=c(1,2))
				SCRAP <- PLOT_IND( samp, 1:2 )
				dev.off()
			}

		}
}

######################################
## MAKE PLOTS ########################
tags.cand <- "Tn_m10"
vars.rand <- Which_Vars$VAR

## PLOT_FIXED
 # Plot Fixed Effects
for ( tag in tags.cand ) {
	PLOT_FIXED( tag, 3 )
}
 # w/ Simulated Baseline Measurements
PLOT_FIXED( "LeadUp_m10", 3 )

# ## PLOT_SHRINK.1
#  # Plot Shrinkage vs BI Model for 1 Variable
# vs_tag <- "BI"
# for ( tag in tags.cand ) {
# 	for ( var in vars.rand ) {
# 		if ( !(var %in% MOD.DAT[[tag]]$M$covs.r) ) { next }
# 		png( paste(PathToPlot,"ModSumm_",tag,".Shrink1.",vs_tag,".",var,".png",sep=""),height=2000,width=2000,pointsize=32)
# 		layout( matrix(c(1:4,4,4),byrow=F,ncol=2), heights=c(2,2,3),widths=c(3,4) )
# 		PLOT_SHRINK.1( tag, vs_tag, var, 1:4 )
# 		dev.off()
# 	}
# }

# ## PLOT_SHRINK.1 (Combined Plots)
#  # Plot Shrinkage (combined plot) vs BI Model for 1 Variable
# vs_tag <- "BI"
# for ( tag in tags.cand ) {
# 	for ( var in vars.rand ) {
# 		if ( !(var %in% MOD.DAT[[tag]]$M$covs.r) ) { next }
# 		png( paste(PathToPlot,"ModSumm_",tag,".Shrink1b.",vs_tag,".",var,".png",sep=""),height=2000,width=2000,pointsize=32)
# 		layout( matrix(c(1:2,3,3),byrow=F,ncol=2), heights=c(4,3),widths=c(3,4) )
# 		PLOT_SHRINK.1( tag, vs_tag, var, c(12,3,4) )
# 		dev.off()
# 	}
# }

## PLOT_SHRINK.1 (Hist & Cor Only)
 # Plot Shrinkage (Hist & Cor Only) vs BI Model for 1 Variable
vs_tag <- "BI"
for ( tag in tags.cand ) {
	for ( var in vars.rand ) {
		if ( !(var %in% MOD.DAT[[tag]]$M$covs.r) ) { next }
		png( paste(PathToPlot,"ModSumm_",tag,".Shrink1c.",vs_tag,".",var,".png",sep=""),height=2000,width=1000,pointsize=32)
		par(mfrow=c(2,1))
		PLOT_SHRINK.1( tag, vs_tag, var, c(1,3) )
		dev.off()
	}
}

# ## PLOT_SHRINK.2:
#  # Plot Shrinkage vs BI Model for 2 Variables
# vs_tag <- "BI"
#  # DRUG vs Intercept
# var_1 <- "Intercept"
# var_2 <- "DRUG"
# for ( tag in tags.cand ) {
# 	PLOT_SHRINK.2( tag, vs_tag, var_1, var_2, 1:10 )
# }
#  # w/ TRT
# for ( tag in grep("10",tags.cand,value=T ) ) {
# 	 # DRUG vs TRT
# 	var_1 <- "TRT"
# 	var_2 <- "DRUG"
# 	PLOT_SHRINK.2( tag, vs_tag, var_1, var_2, 1:10 )
# 	 # DRUG vs Intercept
# 	var_1 <- "Intercept"
# 	var_2 <- "TRT"
# 	PLOT_SHRINK.2( tag, vs_tag, var_1, var_2, 1:10 )
# }

## PLOT_RAND:
 # Plot Pairwise Random Effects for a Model
for ( tag in tags.cand ) {
	PLOT_RAND( tag, 1:10 )
}
 # w/ Simulated Baseline Measurements
PLOT_RAND( "LeadUp_m10", 1:10 )

## PLOT_HEAT_COR:
 # Heatmaps of Correlations b/n Estimates
ARMS <- as.character(unique( FUL$GRP ))
MG.cor <- list()
vs_tag <- "BI"
for ( tag in tags.cand ) {
	if ( grepl("10",tag) ) { vars <- vars.rand }else{ vars <- setdiff(vars.rand,"TRT") }
	for ( subset in c("All",ARMS) ) {
		png( paste(PathToPlot,"ModHeat_",tag,".",vs_tag,".",subset,".png",sep=""),height=1000,width=1000,pointsize=26)
		MG.cor[[subset]] <- PLOT_HEAT_COR( tag, vs_tag, vars, subset )
		dev.off()
	}
}
# MG.diff.GP <- round(MG.cor$G - MG.cor$P, 3)
# MG.diff.GPE <- round(MG.cor$G - MG.cor$PE, 3)
# round( ( MG.diff.GP + MG.diff.GPE ) / 2 ,2 )

## PLOT_ARM_ANOVA:
 # Violin Plot & ANOVA of Variant Estimates across Arms
vs_tag <- "BI"
png( paste(PathToPlot,"ModArm_",tag,".vs.",vs_tag,".png",sep=""),height=1600,width=2000,pointsize=30)
par(mfrow=c(2,2))
PLOT_ARM_ANOVA( tag, vs_tag, "Intercept" )
PLOT_ARM_ANOVA( tag, vs_tag, "DRUG" )
PLOT_ARM_ANOVA( tag, vs_tag, "TRT" )
PLOT_ARM_ANOVA( tag, vs_tag, "WK" )
dev.off()


## PLOT_POST_DENS.2:
 # Plot Posterior Densities for all Patients (Boxplot-ish)
vs_tag <- "BI"
QUANTS <- list()
for ( tag in "Tn_m10" ) {
# for ( tag in tags.cand ) {
# for ( tag in "LeadUp_m10" ) {
	for ( var in c("Intercept","DRUG","TRT") ) {
		if ( grepl("m11",tag) & var=="TRT" ) { next }
		temp.tag <- paste(tag,var,sep=".")
		QUANTS[[temp.tag]] <- PLOT_POST_DENS.2( tag, vs_tag, var )
	}
}
Which_Samps.dr <- gsub(",DRUG]","", head( colnames(QUANTS$Tn_m10.DRUG$Q1), 6 ), fixed=T )
Which_Samps.int <- gsub(",Intercept]","", head( colnames(QUANTS$Tn_m10.Intercept$Q1), 6 ), fixed=T )
Which_Samps <- union( Which_Samps.dr,Which_Samps.int )
Which_Samps.quant <- gsub( "r_ID[","", Which_Samps,fixed=T)

## PLOT_POSTERIORS:
 # Plot Posterior Probabilities of Response
vs_tag <- "BI"
POSTS <- list()
for ( tag in tags.cand ) {
	# for ( var in c("DRUG","TRT") ) {
	for ( var in c("DRUG") ) {
		for ( thresh in c(-1.8,-1.2,-.6,0) ) {
			if ( grepl("m11",tag) & var=="TRT" ) { next }
			temp.tag <- paste(tag,var,thresh,sep=".")
			# POSTS[[temp.tag]] <- PLOT_POSTERIORS( tag, vs_tag, var, thresh )
			scrap <- PLOT_POSTERIORS.2( tag, vs_tag, var, thresh, "ACR50" )
			scrap <- PLOT_POSTERIORS.2( tag, vs_tag, var, thresh, "EULAR" )
		}
	}
}
 # Which Patients have Biggest Disagreement
DIFF <- Which_Samps <- list()
DIFF$DR12 <- apply( POSTS$`Tn_m10.DRUG.-1.2`, 1, diff )
Which_Samps$DR12 <- head( sort( abs(DIFF$DR12),decreasing=T ), 6 )
DIFF$DR0 <- apply( POSTS$`Tn_m10.DRUG.0`, 1, diff )
Which_Samps$DR0 <- head( sort( abs(DIFF$DR0),decreasing=T ), 6 )
Which_Samps.post <- gsub( "r_ID[","", gsub(",DRUG]","", Reduce(union,lapply(Which_Samps,names)), fixed=T),fixed=T)

## PLOT_PRED & PLOT_POST_DENS.1:
 # Plot a few samples (predicted vs actual over time)
Which_Samps.all <- union( Which_Samps.quant, Which_Samps.post)
Which_Samps.all <- sample(Samp.inf$Samps[Samp.inf$GRP!="G"], 40 )
Which_Samps.all <- "B012328"
vs_tag <- "BI"
# for ( tag in tags.cand ) {
for ( tag in "Tn_m10" ) {
	# DIFF <- Which_Samps <- list()
	# DIFF$DR12 <- apply( POSTS$`Tn_m10.DRUG.-1.2`, 1, diff )
	# Which_Samps$DR12 <- head( sort( abs(DIFF$DR12),decreasing=T ), 5 )
	# DIFF$DR0 <- apply( POSTS$`Tn_m10.DRUG.0`, 1, diff )
	# Which_Samps$DR0 <- head( sort( abs(DIFF$DR0),decreasing=T ), 5 )
	# Which_Samps.all <- gsub( "r_ID[","", gsub(",DRUG]","", Reduce(union,lapply(Which_Samps,names)), fixed=T),fixed=T)
	for ( samp in Which_Samps.all ) {	
		# # Predicted & Actual Values of Individuals (over time)
		# png( paste(PathToPlot,"ModSamp_",samp,"_",tag,".vs.",vs_tag,".png",sep=""),height=1600,width=2200,pointsize=36)
		# PLOT_PRED( tag, vs_tag, samp )
		# dev.off()
		# # Posterior Probabilities of Individuals
		# png( paste(PathToPlot,"ModSamp_",samp,"_",tag,".vs.",vs_tag,".",var,".png",sep=""),height=1600,width=1600,pointsize=30)
		# par(mfrow=c(2,1))
		# PLOT_POST_DENS.1( tag, vs_tag, var, samp )
		# dev.off()

		# All Together
		png( paste(PathToPlot,"ModSamp.2_",samp,"_",tag,".vs.",vs_tag,".png",sep=""),height=1200,width=2400,pointsize=30)
		layout( matrix(c(1,1,2,3),ncol=2,byrow=F), widths=c(4,3) )
		PLOT_PRED( tag, vs_tag, samp )
		PLOT_POST_DENS.1( tag, vs_tag, var, samp )
		dev.off()
	}
}

## Plot Within-Patient Error of Estimates
tag <- "Tn_m10"
vs_tag <- "BI"
var_1 <- "Intercept"
var_2 <- "DRUG"
PLOT_WI_PAT_ERR( tag, vs_tag, var_1, var_2 )


######################################
## PLOT LEADUP MODEL SUMMARIES #######

PLOT_FIXED( "LeadUp_m10", 1:10 )
PLOT_RAND( "LeadUp_m10", 1:10 )


PLOT_RAND( "BI", 1:10 )


######################################
## PLOT POSTERIOR RESP PROBS #########
## ...for each patient, ranked #######
temp.ord <- order( POSTS$Tn_m10.DRUG.0$DRUG.0.Tn_m10 )
temp.ord <- order( POSTS$`Tn_m10.DRUG.-1.8`$`DRUG.-1.8.Tn_m10` )
plot( 0,0,type="n", xlim=c(1,421),ylim=c(0,1.2), xlab="Patient",ylab="Pr[x<X]", main="Posterior Response Probability", xaxt="n" )
axis( 2, at=seq(0,1,.2) )
abline( h=seq(0,1,.2), lty=3,col="grey50",lwd=1 )
abline( h=c(0,.5,1) )
points( 1:421, POSTS$Tn_m10.DRUG.0$DRUG.0.Tn_m10[temp.ord], col=adjustcolor(COLS.list.ord[1],.5),pch=16 )
points( 1:421, POSTS$`Tn_m10.DRUG.-0.6`$`DRUG.-0.6.Tn_m10`[temp.ord], col=adjustcolor(COLS.list.ord[2],.5),pch=16, ylim=c(0,1) )
points( 1:421, POSTS$`Tn_m10.DRUG.-1.2`$`DRUG.-1.2.Tn_m10`[temp.ord], col=adjustcolor(COLS.list.ord[3],.5),pch=16, ylim=c(0,1) )
points( 1:421, POSTS$`Tn_m10.DRUG.-1.8`$`DRUG.-1.8.Tn_m10`[temp.ord], col=adjustcolor(COLS.list.ord[4],.5),pch=16, ylim=c(0,1) )
legend( "topleft", legend=paste("X=",rev(c(-1.8,-1.2,-.6,0)),sep=""), col=COLS.list.ord[1:4], pch=16, ncol=4,bg="white" )

######################################
## PLOT ACR/EULAR OUTCOMES ###########

## Plot it
png( paste(PathToPlot,"Categorical_Outcomes.png",sep=""),height=1200,width=1800,pointsize=30)
layout( matrix(1:2,ncol=2), widths=c(5,3) )
 # EULAR
CATS <- as.character(unique( Samp.inf$EULAR_28 ))
Cols.cat <- COLS.list.2[c(6,3,7)]
COMP.eul <- table( Samp.inf$EULAR_28, Samp.inf$EULAR_52 )
temp <- barplot( COMP.eul, beside=T, col=Cols.cat,border=NA, main="EULAR Response at Week 28 vs 52",xlab="Week 52",ylab="# Patients",ylim=c(0,1.1*max(COMP.eul)) )
abline( h=seq(0,500,20), lty=3,col="grey50",lwd=1 )
abline( h=0 )
temp <- barplot( COMP.eul, beside=T, col=Cols.cat,border=NA, main="",xlab="",ylab="", add=T )
legend( "topleft", fill=Cols.cat, legend=CATS, border=NA, title="Week 28")
text( c(temp), c(COMP.eul)+8, c(COMP.eul), srt=90 )
 # ACR50
CATS <- as.character(unique( Samp.inf$ACR50_28 ))
Cols.cat <- COLS.list.2[c(6,7)]
Cols.cat <- COLS.list.2[c(7,6)]
COMP.acr <- table( Samp.inf$ACR50_28, Samp.inf$ACR50_52 )
temp <- barplot( COMP.acr, beside=T, col=Cols.cat,border=NA, main="ACR50 Response",xlab="Week 52",ylab="# Patients",ylim=c(0,1.1*max(COMP.acr)) )
abline( h=seq(0,500,50), lty=3,col="grey50",lwd=1 )
abline( h=0 )
temp <- barplot( COMP.acr, beside=T, col=Cols.cat,border=NA, main="",xlab="",ylab="", add=T )
legend( "topright", fill=Cols.cat, legend=CATS, border=NA, title="Week 28")
text( c(temp), c(COMP.acr)+15, c(COMP.acr), srt=90 )
dev.off()

## (MESSING AROUND)
CATS <- as.character(unique( Samp.inf$EULAR_28 ))
Cols.cat <- COLS.list.2[c(6,3,7)]

png( paste(PathToPlot,"Categorical_Outcomes.vs.Posterior.Boxplots.png",sep=""),height=1000,width=2000,pointsize=30)
par(mfrow=c(1,4))
# Week 28, Mixed Model
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ Samp.inf$EULAR_28, col=Cols.cat, las=2, main="Mixed Model vs EULAR (week 28)",ylab="Pr(x < -1.2)",xlab="" )
abline(h=seq(0,1,.1),lty=3) ; abline(h=seq(0,1,.5))
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ Samp.inf$EULAR_28, col=Cols.cat, las=2, main="Mixed Model vs EULAR (week 28)",ylab="Pr(x < -1.2)",xlab="", add=T )
points( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ as.factor(Samp.inf$EULAR_28), pch=16,col=adjustcolor("black",.6) )
# Week 28, N-of-1 Model
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ Samp.inf$EULAR_28, col=Cols.cat, las=2, main="N-of-1 vs EULAR (week 28)",ylab="Pr(x < -1.2)",xlab="" )
abline(h=seq(0,1,.1),lty=3) ; abline(h=seq(0,1,.5))
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ Samp.inf$EULAR_28, col=Cols.cat, las=2, main="N-of-1 vs EULAR (week 28)",ylab="Pr(x < -1.2)",xlab="", add=T )
points( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ as.factor(Samp.inf$EULAR_28), pch=16,col=adjustcolor("black",.6) )
# Week 52, Mixed Model
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ Samp.inf$EULAR_52, col=Cols.cat, las=2, main="Mixed Model vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="" )
abline(h=seq(0,1,.1),lty=3) ; abline(h=seq(0,1,.5))
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ Samp.inf$EULAR_52, col=Cols.cat, las=2, main="Mixed Model vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="", add=T )
points( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ as.factor(Samp.inf$EULAR_52), pch=1,lwd=2,col=adjustcolor("black",.6) )
# Week 52, N-of-1 Model
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ Samp.inf$EULAR_52, col=Cols.cat, las=2, main="N-of-1 vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="" )
abline(h=seq(0,1,.1),lty=3) ; abline(h=seq(0,1,.5))
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ Samp.inf$EULAR_52, col=Cols.cat, las=2, main="N-of-1 vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="", add=T )
points( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ as.factor(Samp.inf$EULAR_52), pch=1,lwd=2,col=adjustcolor("black",.6) )
dev.off()

png( paste(PathToPlot,"Categorical_Outcomes.vs.Posterior.Boxplots.2.png",sep=""),height=1000,width=1000,pointsize=22)
par(mfrow=c(1,2))
# Agree, Mixed Model
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ Samp.inf$EULAR_52, subset=Samp.inf$EULAR_52==Samp.inf$EULAR_28, col=Cols.cat, las=2, main="Mixed Model vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="" )
abline(h=seq(0,1,.1),lty=3) ; abline(h=seq(0,1,.5))
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ Samp.inf$EULAR_52, subset=Samp.inf$EULAR_52==Samp.inf$EULAR_28, col=Cols.cat, las=2, main="Mixed Model vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="", add=T )
points( POSTS$`Tn_m10.DRUG.-1.2`[,1] ~ as.factor(Samp.inf$EULAR_52), subset=Samp.inf$EULAR_52==Samp.inf$EULAR_28, pch=16,col=adjustcolor("black",.6) )
# Agree, N-of-1 Model
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ Samp.inf$EULAR_52, subset=Samp.inf$EULAR_52==Samp.inf$EULAR_28, col=Cols.cat, las=2, main="N-of-1 vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="" )
abline(h=seq(0,1,.1),lty=3) ; abline(h=seq(0,1,.5))
boxplot( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ Samp.inf$EULAR_52, subset=Samp.inf$EULAR_52==Samp.inf$EULAR_28, col=Cols.cat, las=2, main="N-of-1 vs EULAR (week 52)",ylab="Pr(x < -1.2)",xlab="", add=T )
points( POSTS$`Tn_m10.DRUG.-1.2`[,2] ~ as.factor(Samp.inf$EULAR_52), subset=Samp.inf$EULAR_52==Samp.inf$EULAR_28, pch=16,col=adjustcolor("black",.6) )
dev.off()


#############################################################
## PLOT MODEL COMPARISONS ###################################
##### b/n models ############################################

######################################
## COMPARE DRUG/INT/TRT ESTIAMTES ####


## Drug vs Intercept
vars <- c("Int","Drug")
 # Model 9
# mods <- c("MN","T_m9","Tc_m9","Tcdi_m9","Tcid_m9","BI")
# png( paste(PathToPlot,"1-Mod_Scatter.INTvDR.9.png",sep=""), height=500*length(mods),width=500*length(mods),pointsize=32 )
# Scatter_Mod_Pairs( COEF, vars, mods )
# dev.off()
 # Model 10
mods <- c("MN","T_m10","Tc_m10","Tcdi_m10","Tcid_m10","BI")
mods <- c("T_m10","Tc_m10","BI")
png( paste(PathToPlot,"Mod_Scatter.INTvDR.10.png",sep=""), height=500*length(mods),width=500*length(mods),pointsize=32 )
Scatter_Mod_Pairs( COEF, vars, mods )
dev.off()
 # Model 11
mods <- c("MN","T_m11","Tc_m11","Tcdi_m11","Tcid_m11","BI")
mods <- c("T_m11","Tc_m11","BI")
png( paste(PathToPlot,"Mod_Scatter.INTvDR.11.png",sep=""), height=500*length(mods),width=500*length(mods),pointsize=32 )
Scatter_Mod_Pairs( COEF, vars, mods )
dev.off()

 # Model w/o Cor
mods <- c("MN","T_m9","T_m10","T_m11","BI")
mods <- c("T_m10","T_m11","BI")
png( paste(PathToPlot,"Mod_Scatter.INTvDR.T.png",sep=""), height=500*length(mods),width=500*length(mods),pointsize=32 )
Scatter_Mod_Pairs( COEF, vars, mods )
dev.off()
 # Model w/ Cor | ID
mods <- c("MN","Tc_m9","Tc_m10","Tc_m11","BI")
mods <- c("Tc_m10","Tc_m11","BI")
png( paste(PathToPlot,"Mod_Scatter.INTvDR.Tc.png",sep=""), height=500*length(mods),width=500*length(mods),pointsize=32 )
Scatter_Mod_Pairs( COEF, vars, mods )
dev.off()

## Treatment vs Drug
vars <- c("Trt","Drug")
 # Model 9
mods <- c("T_m9","Tc_m9","Tcdi_m9","Tcid_m9","BI")
png( paste(PathToPlot,"1-Mod_Scatter.TRTvDR.9.png",sep=""), height=500*length(mods),width=500*length(mods),pointsize=32 )
Scatter_Mod_Pairs( COEF, vars, mods )
dev.off()
 # Model 10
mods <- c("T_m10","Tc_m10","Tcdi_m10","Tcid_m10","BI")
png( paste(PathToPlot,"1-Mod_Scatter.TRTvDR.10.png",sep=""), height=500*length(mods),width=500*length(mods),pointsize=32 )
Scatter_Mod_Pairs( COEF, vars, mods )
dev.off()





COLS <- adjustcolor( COLS.list.2[3:5], .5 )
par(mfrow=c(1,3))
for ( var in names(COEF)[c(1,2,4)] ) {
	plot( T_m11 ~ T_m10, data=COEF[[var]], pch=16,col=COLS[factor(Samp.inf$GRP)], main=var )
	abline(0,1)
}

#############################################################
## END OF DOC ###############################################
#############################################################


## FCT: Plot Model
PLOT_MOD <- function( model, tag, plot_rand ) {
	write(paste(date(),"Model:",tag), paste(PathToPlot,"Update.txt",sep=""),append=T)

	# ## Collect General Model Info
	# summ <- summary(model,waic=F)
	# m.obs <- summ$nobs
	# m.iter <- summ$iter
	# m.warm <- summ$warmup
	# m.chain <- summ$chains
	# # m.waic <- summ$WAIC
	# m.pheno <- as.character(model$formula)[2]
	# m.covs <- as.character(model$formula)[3]
	# m.form <- paste( m.pheno, "~", m.covs )
	# RAND <- length(summ$random)>0

	# ## Collect Model Outputs
	# d.prior <- model$prior
	# d.post <- posterior_samples(model)
	# f.eff <- summ$fixed
	# m.covs.f <- rownames(f.eff)
	# c.eff <- summ$cor_pars
	# s.eff <- summ$spec_pars
	# d.post.f.which <- paste("b_",m.covs.f,sep="")
	# d.post.c.which <- rownames(c.eff)
	# if ( RAND==T ) {
	# 	r.grp <- summ$group
	# 	r.ngrps <- summ$ngrps
	# 	r.eff <- ranef(model)[[1]]
	# 	r.eff.2 <- summ$random
	# 	m.covs.r <- colnames(r.eff)[1:ncol(r.eff)]
	# 	m.samps <- rownames(r.eff)
	# 	n.samps <- length(m.samps)
	# 	all.eff <- list( f.eff, c.eff, Reduce( rbind, r.eff.2 ), s.eff )
	# 	mod.fit <- Reduce( rbind, all.eff )
	# 	mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","R","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
	# 	mod.fit$Eff[grep("cor(",rownames(mod.fit),fixed=T)] <- "RC"
	# 	d.post.r.which <- unlist(lapply( m.covs.r, function(x)paste( "sd_",r.grp,"_",x,sep="") ))
	# 	d.post.r2.which <- unlist(lapply( m.covs.r, function(x)paste( "r_",r.grp,"[",m.samps,",",x,"]",sep="") ))
	# }else{
	# 	all.eff <- list( f.eff, c.eff, s.eff )
	# 	mod.fit <- Reduce( rbind, all.eff )
	# 	mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
	# }
	# print("Model Parsed")

	# ###################################
	# ## FIXED EFFECTS PLOTS ############
	#  # Set Color Palette
	# COLS.list.heat <- c("firebrick3","chocolate2","gold1","springgreen1","steelblue2","slateblue3")
	# COLS <- adjustcolor(COLS.list.2,alpha=.6)
	# COLS.eff <- COLS[c(1:4,6)]
	# names(COLS.eff) <- c("F","C","R","RC","S")

	# ## Plot Model Summary
	# print("Plotting #1 - Model Summary")
	# png( paste(PathToPlot,"ModSumm_",tag,".1-PlotFct.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	# par(mfrow=c(nrow(mod.fit),2))
	# plot( model, N=nrow(mod.fit) )
	# dev.off()

	# # ## Plot Chains
	# #  # (aka) How to Pull all Sampling Data
	# # which_vars <- m.covs.f
	# # which_vars <- grep( "sd(",rownames(mod.fit),fixed=T,value=T )
	# # png( paste(PathToPlot,"ModSumm_",tag,".1-Chains.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	# # par( mfrow=c(length(which_vars),1) )
	# # for ( v in which_vars ) {
	# # 	v.tag <- grep( substr(gsub("(",paste("_",r.grp,"_",sep=""),v,fixed=T),1,10), model$fit@sim$fnames_oi )
	# # 	YLIM <- Reduce(range,lapply(model$fit@sim$samples,function(x)range(x[[paste("b",v,sep="_")]])))
	# # 	plot( 0,0,type="n",xlim=c(0,m.iter),ylim=YLIM,xlab="Iteration",ylab=v,main="Chains" )
	# # 	abline( h=-100:100,lty=3,col="grey50",lwd=1 )
	# # 	# SCRAP <- lapply( seq(YLIM[1],YLIM[2],.025),function(x)abline(h=x,col=adjustcolor(COLS.list.2[4],alpha=2*dnorm(x,5,1)),lwd=2 ))
	# # 	for ( c in 1:m.chain ) {
	# # 		# TEMP <- model$fit@sim$samples[[c]]$b_Intercept
	# # 		points( 1:m.iter, model$fit@sim$samples[[c]][[paste("b",v,sep="_")]], type="l",col=adjustcolor(COLS.list.2[c],alpha=.7),lwd=2 )
	# # 	}
	# # }
	# # dev.off()
	
	# ## Fixed Effect Sizes
	# print("Plotting #2 - Effect Sizes")
	# YLIM <- extendrange(mod.fit[,"Estimate"], f=.2)
	# png( paste(PathToPlot,"ModSumm_",tag,".2-EffSize.png",sep=""),height=1000,width=400+100*nrow(mod.fit),pointsize=26)
	# par(mar=c( 7,5,5,3 ))
	# # TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n",ylim=YLIM,ylab="Effect Size",main="Effect Size Estimates" )
	# TEMP <- 1:nrow(mod.fit)
	# plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=range(TEMP),ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
	# axis( 1, at=TEMP,label=rownames(mod.fit), las=2 )
	# axis( 2, at=seq(-10,10,2), las=2 )
	# abline( h=-10:10, lty=3,col="grey50",lwd=1 )
	# abline( h=0, lty=1,col="black",lwd=1 )
	#  # Plot Prior Distributions
	# for ( v in 1:nrow(mod.fit) ) {
	# 	var <- rownames(mod.fit)[v]
	# 	if ( var %in% m.covs.f ) {
	# 		if ( var=="Intercept" ) {
	# 			temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
	# 		}else{
	# 			temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
	# 		}
	# 		if ( length(temp.priors)==2 ) {
	# 			names(temp.priors) <- c("Mean","SD")
	# 			vioplot( rnorm(1e5,temp.priors["Mean"],temp.priors["SD"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
	# 		}
	# 	}
	# 	# if ( var %in% paste("sd(",m.covs.r,")",sep="") ) {
	# 	# 	temp.priors <- as.numeric(strsplit( gsub(")","",gsub("cauchy(","",d.prior[which(d.prior$class=="sd"&d.prior$group==r.grp)[1],"prior"], fixed=T),fixed=T),"," )[[1]])
	# 	# 	names(temp.priors) <- c("Loc","Scale")
	# 	# 	vioplot( rnorm(1e5,temp.priors["Loc"],temp.priors["Scale"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )
	# 	# }
	# }
	#  # Plot Posterior Distributions
	# for ( v in 1:nrow(mod.fit) ) {
	# 	var <- rownames(mod.fit)[v]
	# 	var.tag <- var
	# 	if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
	# 	if ( grepl("sd(",var,fixed=T) ) { var.tag <- gsub(")","", gsub("sd(",paste("sd_",r.grp,"_",sep=""),var,fixed=T),fixed=T) }
	# 	if ( grepl("cor(",var,fixed=T) ) { var.tag <- gsub(",","_", gsub(")","", gsub("cor(",paste("cor_",r.grp,"_",sep=""),var,fixed=T),fixed=T),fixed=T) }
	# 	if ( var==paste("sigma(",m.pheno,")",sep="") ) { var.tag <- paste("sigma",m.pheno,sep="_") }
	# 	if ( var.tag %in% colnames(d.post) ) {
	# 		vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.eff[mod.fit[v,"Eff"]],add=T,drawRect=F )
	# 	}
	# }
	# # TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n", add=T)
	# arrows( TEMP,mod.fit[,"l.95..CI"],TEMP,mod.fit[,"u.95..CI"],lwd=3,length=0 )
	# arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],lwd=3,length=0 )
	# # legend("topright",fill=COLS.eff,border=NA,legend=names(COLS.eff),ncol=length(COLS.eff),title="Effect Type",bg="white")
	# legend("topright",fill=c(adjustcolor("black",alpha=.2),COLS.eff),border=NA,legend=c("Prior",names(COLS.eff)),ncol=(1+length(COLS.eff))/2,title="Effect Type",bg="white")
	# dev.off()

	# ## RANDOM EFFECTS PLOTS ###########
	# if ( RAND==T & plot_rand!=F ) {
	# 	print("Plotting Random Effects")
	# 	n.covs.r <- length(m.covs.r)

	# 	## Plot Distributions & Shrinkage of Random Effects
	# 	print("Plotting #R1 - Distributions")
	# 	PLOT_SHRINK <- function( var, which_plots ) {
	# 		## Pull Data
	# 		brm.var <- r.eff[,var] + f.eff[var,"Estimate"]
	# 		if ( var %in% c("Intercept","DRUG") ) {
	# 			if ( var=="Intercept" ) {
	# 				mn.var <- FUL[,"DAS_BL_MN"]
	# 				brm.fx <- f.eff[var,"Estimate"]
	# 			}
	# 			if ( var=="DRUG" ) {
	# 				mn.var <- FUL[,"DEL_MNe_MN"]
	# 				brm.fx <- f.eff[var,"Estimate"] + f.eff["TRT","Estimate"]
	# 				if ( "TRT"%in%m.covs.f ) { brm.var <- brm.var + f.eff["TRT","Estimate"] }
	# 			}
	# 			mn.mod <- T
	# 			names(mn.var) <- as.character(FUL$ID)
	# 		}else{
	# 			mn.var <- brm.var
	# 			mn.mod <- F
	# 			brm.fx <- f.eff[var,"Estimate"]
	# 		}
	# 		mg.var <- merge( data.frame(brm.var),data.frame(mn.var), by="row.names" )
	# 		mg.var <- merge( FUL[,c("ID","GRP")], mg.var, by.x="ID",by.y="Row.names" )
	# 		mg.var.ord <- mg.var[ order(mg.var$mn.var), ]
	# 		mn.fx <- mean(mg.var$mn.var)

	# 		## Plot 1: Distributions
	# 		if ( 1 %in% which_plots ) {
	# 			XLIM <- range( c(mg.var$mn.var,mg.var$brm.var) )
	# 			X_BIN <- .25
	# 			BRKS <- seq( floor(XLIM[1]), ceiling(XLIM[2])+X_BIN, X_BIN)
	# 			YLIM <- c( 0, Reduce( max, lapply(list(mg.var$brm.var,mg.var$mn.var),function(x)hist(x,breaks=BRKS,plot=F)$counts) ) )
	# 			hist( mg.var$brm.var, col=COLS[1],border=NA,breaks=BRKS, xlim=XLIM,ylim=YLIM, main=paste(var,"- Mean vs BRMS"),xlab=var,ylab="# Patients" )
	# 			abline( v=brm.fx, col=COLS.list.2[1],lwd=4,lty=1 )
	# 			if ( mn.mod==T ) {
	# 				hist( mg.var$mn.var, col=COLS[2],border=NA,breaks=BRKS, add=T )
	# 				abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=1 )
	# 				legend( "topleft", legend=c("BRMS","Mean"),title="Model",fill=COLS[1:2],border=NA)	
	# 			}			
	# 		}
	
	# 		## Plot 2: Boxplot by Arm
	# 		if ( 2 %in% which_plots ) {
	# 			ARMS <- as.character( unique( mg.var$GRP ) )
	# 			XLIM <- c(0,7)
	# 			YLIM <- range( c(mg.var$mn.var,mg.var$brm.var) )
	# 			plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",xlab="ARM",ylab=var,main=paste(var,"Distribution by Arm") )
	# 			abline( h=-10:10,lty=3,col="grey50",lwd=1 )
	# 			axis( 1,at=seq(1.5,7,2),label=ARMS )
	# 			SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$brm.var[mg.var$GRP==ARMS[x]],at=2*x-1,col=COLS[1],border=NA,add=T) )
	# 			abline( h=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
	# 			if ( mn.mod==T ) {
	# 				SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$mn.var[mg.var$GRP==ARMS[x]],at=2*x,col=COLS[2],border=NA,add=T) )
	# 				abline( h=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
	# 				legend( "topleft", legend=c("BRMS","Mean"),title="Model",fill=COLS[1:2],border=NA)	
	# 				arrows( c(1,3,5),rep(YLIM[1],3),c(2,4,6),rep(YLIM[1],3), length=0,lwd=4,col=COLS.list.2[3:5] )
	# 			}		
	# 		}

	# 		## Plot 3: Shrinkage by Person
	# 		if ( 3 %in% which_plots ) {
	# 			XLIM <- range( c(mg.var.ord$mn.var,mg.var.ord$brm.var) )
	# 			X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5)
	# 			YLIM <- c( 1, n.samps )
	# 			plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=var,ylab="Patient",main=paste(var,"- Estimate Shrinkage") )
	# 			abline( v=X_LNS, lty=3,col="grey50",lwd=1 )
	# 			points( mg.var.ord$brm.var,1:nrow(mg.var.ord), pch=16,cex=1,col=COLS[1] )
	# 			if ( mn.mod==T ) {
	# 				points( mg.var.ord$mn.var,1:nrow(mg.var.ord), pch=16,cex=1,col=COLS[2] )
	# 				arrows( mg.var.ord$mn.var,1:nrow(mg.var.ord),mg.var.ord$brm.var,1:nrow(mg.var.ord), angle=30,length=.2,lwd=2,col=COLS.list.2[3:5][factor(mg.var.ord$GRP)] )
	# 				abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
	# 			}
	# 			abline( v=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
	# 			legend("topleft",col=COLS.list.2[1:5],legend=c("Mod: BRMS","Mod: Mean","Arm: G","Arm: P","Arm: PE"),pch=c(16,16,NA,NA,NA),lty=c(0,0,1,1,1),lwd=2 )
	# 		}

	# 		## Plot 4: MN vs BRMS Model Estimates
	# 		if ( 4 %in% which_plots ) {
	# 			XLIM <- range( mg.var$mn.var )
	# 			YLIM <- range( mg.var$brm.var )
	# 			X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5 )
	# 			Y_LNS <- seq( floor(YLIM[1]-5), ceiling(YLIM[2]+5), .5 )
	# 			plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Mean Estimate",ylab="BRMS Estimate",main=var )
	# 			abline( v=X_LNS,h=Y_LNS, lty=3,col="grey50",lwd=1 )
	# 			abline( 0,1 )
	# 			abline( h=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
	# 			abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
	# 			points( mg.var$brm.var ~ mg.var$mn.var, data=mg.var, pch=16,cex=1,col=COLS[3:5][factor(GRP)] )
	# 			legend( "topleft",pch=16,col=COLS[3:5],legend=levels(factor(mg.var$GRP)))
	# 		}
	# 	}
	# 	 # Plot it
	# 	for ( var in m.covs.r ) {
	# 		if ( var %in% c("Intercept","DRUG") ) {
	# 			png( paste(PathToPlot,"ModSumm_",tag,".R1-Dist.",var,".png",sep=""),height=2000,width=2000,pointsize=32)
	# 			layout( matrix(c(1:4,4,4),byrow=F,ncol=2) )
	# 			PLOT_SHRINK( var, c(1,2,4) )
	# 			PLOT_SHRINK( var, c(3) )
	# 			dev.off()
	# 		}else{
	# 			png( paste(PathToPlot,"ModSumm_",tag,".R1-Dist.",var,".png",sep=""),height=2000,width=2000,pointsize=32)
	# 			layout( matrix(c(1:2,3,3),byrow=F,ncol=2) )
	# 			PLOT_SHRINK( var, c(1,2) )
	# 			PLOT_SHRINK( var, c(3) )
	# 			dev.off()
	# 		}

	# 	}
	
	# 	## Plot Pairs of Random Effects
	# 	panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...) {
	# 		usr <- par("usr"); on.exit(par(usr))
	# 		par(usr = c(0, 1, 0, 1))
	# 		r <- abs(cor(x, y))
	# 		txt <- format(c(r, 0.123456789), digits = digits)[1]
	# 		txt <- paste0(prefix, txt)
	# 		if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	# 		text(0.5, 0.5, txt, cex = cex.cor * r)
	# 	}
	# 	print("Plotting #R1 - Pairs")
	# 	if ( n.covs.r > 1 ) {
	# 		temp.reff <- t( f.eff[m.covs.r,"Estimate"] + t(r.eff) ) 
	# 		# temp.grp <- FUL$GRP[ match(rownames(temp.reff),FUL$ID) ] # col=adjustcolor(COLS.list.2[factor(temp.grp)],alpha=.4)
	# 		temp.reff <- merge( FUL[,c("ID","GRP")], temp.reff, by.x="ID",by.y="row.names" )
	# 		png( paste(PathToPlot,"ModSumm_",tag,".R1-PairRand.png",sep=""),height=400*n.covs.r,width=400*n.covs.r,pointsize=24+2*n.covs.r)
	# 		# pairs( temp.reff, pch=16,col=adjustcolor(COLS.eff["RC"],alpha=.5),cex=2,upper.panel=panel.cor )
	# 		pairs( temp.reff[,m.covs.r], pch=16,col=COLS[3:5][factor(temp.reff$GRP)],cex=2,upper.panel=panel.cor )
	# 		dev.off()
	# 	}

	# 	## Boxplot of Individual Posterior Distributions
	# 	print("Plotting #R2 - Boxplot")
	# 	# png( paste(PathToPlot,"ModSumm_",tag,".R2-BoxRand.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
	# 	# par(mfrow=c(n.covs.r,1))
	# 	# par(mar=c(6,4,4,1))
	# 	# for ( z in 1:n.covs.r ) {
	# 	# 	z.name <- m.covs.r[z]
	# 	# 	z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
	# 	# 	z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
	# 	# 	z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
	# 	# 	if ( z>0 ) { z.ord <- order( colMeans(z.temp) ) }
	# 	# 	YLIM <- range( z.temp )
	# 	# 	boxplot( z.temp[,z.ord], col=COLS.eff["R"],main=paste("Random",z.name,"by Patient"),ylim=YLIM,ylab=z.name,names=m.samps[z.ord],las=2,pch=16 )
	# 	# 	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
	# 	# 	abline( h=0, lty=1,col="black",lwd=1 )
	# 	# 	boxplot( z.temp[,z.ord], col=COLS.eff["R"],names=m.samps[z.ord],las=2,pch=16,add=T )
	# 	# 	if ( any(z.name %in% c("DRUG","PLAC","TRT")) ) {
	# 	# 		abline(h=c(0,-1),lwd=2,col=COLS.list.2[6:7] )
	# 	# 		text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.5 )
	# 	# 		text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[7],cex=.5 )
	# 	# 	}else{
	# 	# 		legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
	# 	# 	}			
	# 	# }
	# 	# dev.off()
	# 	 # Plot Confidence Intervals of Individual Posterior Distributions (quicker)
	# 	png( paste(PathToPlot,"ModSumm_",tag,".R2-BoxRand.2.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
	# 	par(mfrow=c(n.covs.r,1))
	# 	par(mar=c(6,4,4,1))
	# 	for ( z in 1:n.covs.r ) {
	# 		# Pull Data
	# 		z.name <- m.covs.r[z]
	# 		z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
	# 		z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
	# 		z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
	# 		# z.samps.temp <- gsub("r_ID[","", sapply(strsplit(colnames(z.temp),","),"[",1) ,fixed=T)
	# 		z.grp <- FUL$GRP[ match(m.samps,as.character(FUL$ID)) ]
	# 		which_quants <- c(.025,.25,.5,.75,.975)
	# 		z.quants <- apply( z.temp, 2, function(x)quantile(x,which_quants) )
	# 		if ( z>0 ) { z.ord <- order( colMeans(z.temp) ) }
	# 		# Parameters
	# 		XLIM <- c( 1,ncol(z.temp) )
	# 		YLIM <- range( z.temp )
	# 		# Build Plot
	# 		plot( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp[z.ord]], main=paste("Random",z.name,"by Patient"),xlab="",ylab=z.name,xaxt="n",xlim=XLIM,ylim=YLIM )
	# 		axis( 1, at=1:ncol(z.temp),label=m.samps[z.ord],las=2,cex=.6 )
	# 		abline( h=-10:10,v=1:ncol(z.temp), lty=3,col="grey50",lwd=1 )
	# 		abline( h=0, lty=1,col="black",lwd=1 )
	# 		# Text/Legend
	# 		if ( any(z.name %in% c("DRUG","PLAC","TRT")) ) {
	# 			abline(h=c(0,-1),lwd=2,col=COLS.list.2[6:7] )
	# 			Post.eff0 <- round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1)
	# 			Post.eff1 <- round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1)
	# 			text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=Post.eff0,col=COLS.list.2[6],cex=.7,srt=90 )
	# 			text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=Post.eff1,col=COLS.list.2[7],cex=.7,srt=90 )
	# 		}else{
	# 			# legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
	# 			legend( "topleft",pch=c(rep(16,3),16,1),cex=1.5,col=c(COLS.list.2[3:5],rep("black",2)),legend=c("Arm: G","Arm: P","Arm: PE","Median","CI: 95%"),lwd=3 )
	# 		}
	# 		# Points/Arrows
	# 		arrows( 1:ncol(z.temp),z.quants["25%",z.ord],1:ncol(z.temp),z.quants["75%",z.ord],col=COLS[1],lwd=5,code=3,angle=90,length=.1 )
	# 		points( rep(1:ncol(z.temp),2),c(z.quants["2.5%",z.ord],z.quants["97.5%",z.ord]),pch=1,lwd=3,col=COLS.list.2[3:5][z.grp[z.ord]],cex=1.5 )
	# 		points( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp[z.ord]],cex=1.5 )
	# 	}
	# 	dev.off()

	# 	# ## Violin Plot: Prior vs Posterior (Compiled)
	# 	# print("Plotting #R2 - Violin Plot")
	# 	# png( paste(PathToPlot,"ModSumm_",tag,".R2b-ViolRand.png",sep=""),height=700*n.covs.r,width=800+22*n.samps,pointsize=24+2*n.covs.r)
	# 	# par(mfrow=c(n.covs.r,1))
	# 	# par(mar=c(6,4,4,1))
	# 	# for ( z in 1:n.covs.r ) {
	# 	# 	z.name <- m.covs.r[z]
	# 	# 	z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
	# 	# 	z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
	# 	# 	z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
	# 	# 	if ( z==1 ) { z.ord <- order( colMeans(z.temp) ) }
	# 	# 	XLIM <- c(1,ncol(z.temp))
	# 	# 	YLIM <- extendrange(z.temp,f=.2)
	# 	# 	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="",main="Posterior Population Estimate (by Patient)",xaxt="n" )
	# 	# 	axis( 1, at=1:ncol(z.temp), label=m.samps[z.ord], las=2)
	# 	# 	abline(h=-10:10,lty=3,col="grey50",lwd=1 )
	# 	# 	abline(h=0)
	# 	# 	# Get Prior
	# 	# 	if ( z.name=="Intercept" ) { temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
	# 	# 	}else{ temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==z.name,"prior"], fixed=T),fixed=T),"," )[[1]]) }
	# 	# 	names(temp.priors) <- c("Mean","SD")
	# 	# 	temp.prior.distr <- rnorm(10000,temp.priors["Mean"],temp.priors["SD"])
	# 	# 	# Plot Prior/Posterior
	# 	# 	SCRAP <- lapply( 1:ncol(z.temp), function(x)vioplot(temp.prior.distr,at=x, col=adjustcolor("black",alpha=.2),border=NA,add=T ))
	# 	# 	SCRAP <- lapply( 1:ncol(z.temp), function(x)vioplot(z.temp[,z.ord[x]],at=x, col=COLS.eff["R"],add=T ))
	# 	# 	if ( any(z.name %in% c("DRUG","PLAC")) ) {
	# 	# 		abline(h=c(0,-1),lwd=2,col=COLS.list.2[5:6] )
	# 	# 		text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[5],cex=.5 )
	# 	# 		text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.5 )
	# 	# 	}else{
	# 	# 		legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
	# 	# 	}			
	# 	# }
	# 	# dev.off()

	# 	## Plot Posterior Probabilities for Random Effects
	# 	print("Plotting #R3 - Heatmaps of Posteriors")
	# 	 # Plotting Parameters
	# 	COLS.heat <- colorRampPalette(c(COLS.list.heat[1:3],"white",COLS.list.heat[4:6]))(100)
	# 	BRKS.heat <- seq(0,1,length.out=101)
	# 	n.iter.tot <- m.chain * ( m.iter - m.warm )
	# 	post.prob.brks <- list( Intercept=seq(3.5,7.5,.125),
	# 		DRUG=seq(-3,0,.125),
	# 		PLAC=seq(-2,1,.125),
	# 		TRT=seq(-2,1,.125) )
	# 	 # Calculate/Plot Posteriors at Various Thresholds
	# 	post.prob <- list()
	# 	for ( z in 1:n.covs.r ) {
	# 		z.name <- m.covs.r[z]
	# 		z.col.1 <- d.post.f.which[ d.post.f.which==paste("b",z.name,sep="_") ] # grep(z.name,d.post.f.which,value=T)
	# 		z.col.2 <- grep( paste(",",z.name,"]",sep=""),d.post.r2.which,value=T)
	# 		z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
	# 		z.grp <- FUL$GRP[ match(m.samps,as.character(FUL$ID)) ]
	# 		# Calculate Posteriors
	# 		post.prob[[z.name]] <- Reduce( cbind, lapply( post.prob.brks[[z.name]],function(x)apply(z.temp,2,function(y)length(which(y<x))) ) ) / n.iter.tot
	# 		colnames(post.prob[[z.name]]) <- paste("LT",post.prob.brks[[z.name]],sep="_")
	# 		rownames(post.prob[[z.name]]) <- colnames(z.temp)
	# 		if ( z.name=="Intercept" ) { post.prob[[z.name]] <- 1 - post.prob[[z.name]] ; colnames(post.prob[[z.name]]) <- gsub("LT","GT",colnames(post.prob[[z.name]])) }
	# 		# Plot Heatmap
	# 		png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.",z.name,".png",sep=""),height=1200,width=2000+10*n.samps,pointsize=28)
	# 		heatmap.2( t(post.prob[[z.name]]), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5+n.samps/100),main=z.name,ColSideColors=COLS.list.2[3:5][z.grp] )
	# 		dev.off()
	# 	}
	# 	 # Heatmap of All Random Effect Posteriors
	# 	z.grp <- FUL$GRP[ match(m.samps,as.character(FUL$ID)) ]
	# 	post.prob.temp <- Reduce( cbind, post.prob )
	# 	png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.All.png",sep=""),height=800+400*n.covs.r,width=2000+10*n.samps,pointsize=28)
	# 	heatmap.2( t(post.prob.temp), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5+n.samps/100),rowsep=cumsum(lapply(post.prob,ncol)),ColSideColors=COLS.list.2[3:5][z.grp] )
	# 	dev.off()

	# 	# png( paste(PathToPlot,"ModSumm_",tag,".R3-PostScatter.png",sep=""),height=800+400*n.covs.r,width=2000+10*n.samps,pointsize=28)
	# 	# par(mfrow=c(n.covs.r,1))
	# 	# for ( z in 1:n.covs.r ) {
	# 	# 	z.name <- m.covs.r[z]
	# 	# 	temp.split <- as.numeric(sapply(strsplit(colnames(post.prob[[z.name]]),"_"),"[",2))
	# 	# 	which_cols <- which( temp.split %in% unique(round(2*temp.split)/2) )
	# 	# 	z.ord <- order( post.prob[[z.name]][,which.max(apply(post.prob[[z.name]],2,sd))] )
	# 	# 	COLS.temp <- colorRampPalette(COLS.list.ord)(length(which_cols))
	# 	# 	plot( 0,0,type="n", xlim=c(1,n.samps),ylim=c(0,1) )
	# 	# 	lapply( 1:length(which_cols),function(x)points(1:n.samps,post.prob[[z.name]][z.ord,which_cols[x]],col=adjustcolor(COLS.temp[x],alpha=.7),pch=16,type="o",lwd=2) )
	# 	# }
	# 	# dev.off()
		
	# 	## Plot Distribution (across patients) of Standard Deviations (across iteration) of Random Effects
	# 	# TEST <- lapply(m.covs.r,function(r)unlist(lapply( m.samps, function(s) sd(d.post[,paste("r_ID[",s,",",r,"]",sep="")]) )) )
	# 	# names(TEST) <- m.covs.r

	# 	## Plot DRUG vs Intercept for each Simulation/Patient
	# 	# if ( all(c("DRUG","Intercept")%in%m.covs.r) ) {
	# 	# 	print("Plotting #R4 - Drug v Int per Patient")
	# 	# 	x.name <- "Intercept"
	# 	# 	x.col.1 <- d.post.f.which[ d.post.f.which==paste("b",x.name,sep="_") ] # grep(x.name,d.post.f.which,value=T)
	# 	# 	x.col.2 <- grep( paste(",",x.name,"]",sep=""),d.post.r2.which,value=T)
	# 	# 	x.temp <- d.post[,x.col.1] + d.post[,x.col.2]
	# 	# 	y.name <- "DRUG"
	# 	# 	y.col.1 <- d.post.f.which[ d.post.f.which==paste("b",y.name,sep="_") ] # grep(y.name,d.post.f.which,value=T)
	# 	# 	y.col.2 <- grep( paste(",",y.name,"]",sep=""),d.post.r2.which,value=T)
	# 	# 	y.temp <- d.post[,y.col.1] + d.post[,y.col.2]
	# 	# 	XLIM <- range(x.temp)
	# 	# 	YLIM <- range(y.temp)
	# 	# 	COLS.ind <- colorRampPalette(COLS.list)(ncol(x.temp))
	# 	# 	COLS.ind.2 <- adjustcolor(COLS.ind,alpha=.01)
	# 	# 	png( paste(PathToPlot,"ModSumm_",tag,".R4-Sims.DRvINT.png",sep=""),height=1000,width=1600,pointsize=30)
	# 	# 	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="Drug",main="Individual Simulation Estimates")
	# 	# 	abline(h=-10:10,v=-10:10,lty=3,col="grey50",lwd=1 )
	# 	# 	abline(h=0)
	# 	# 	SCRAP <- lapply( 1:ncol(x.temp),function(i)points(x.temp[,i],y.temp[,i],col=COLS.ind.2[i],pch=16) )
	# 	# 	mod.temp <- lm(unlist(y.temp)~unlist(x.temp))
	# 	# 	abline( mod.temp )
	# 	# 	points( colMeans(x.temp),colMeans(y.temp),pch="+",col=adjustcolor(COLS.ind,alpha=.7),cex=2,lwd=2 )
	# 	# 	dev.off()
	# 	# }

	# 	## Plot a Few Individual Patients' Profiles
	# 	if ( all(c("DRUG","Intercept")%in%m.covs.r) & any(c("PLAC","TRT")%in%m.covs.r) ) {
	# 	# if ( all(c("DRUG","PLAC","Intercept")%in%m.covs.r) ) {
	# 		print("Plotting Individual Patients")
	# 		## PLAC or TRT?
	# 		PBO <- grep("PLAC|TRT",m.covs.r,value=T )
	# 		## Which Patients?
	# 		 # How Many?
	# 		plot.rows <- 4
	# 		plot.cols <- 5
	# 		z.samps.n <- plot.rows * plot.cols
	# 		 # Specify Samples
	# 		# z.samps <- sample(m.samps, z.samps.n )
	# 		z.samps.temp.ord.int <- order( colMeans(d.post[,grep( paste(",Intercept]",sep=""),d.post.r2.which,value=T)]) )
	# 		z.samps.temp.ord.dr <- order( colMeans(d.post[,grep( paste(",DRUG]",sep=""),d.post.r2.which,value=T)]) )
	# 		z.samps.temp <- m.samps[ c( head(z.samps.temp.ord.int,plot.rows),tail(z.samps.temp.ord.int,plot.rows), head(z.samps.temp.ord.dr,plot.rows),tail(z.samps.temp.ord.dr,plot.rows) ) ]
	# 		z.samps <- sample( setdiff(m.samps,z.samps.temp), z.samps.n-length(z.samps.temp) )
	# 		z.samps <- c( head(z.samps.temp,2*plot.rows), z.samps, tail(z.samps.temp,2*plot.rows) )
	# 		# Pull out Sample Data
	# 		z.data <- model$data[model$data[,r.grp]%in%z.samps,]
	# 		z.temp <- predict( model, newdata=z.data )
	# 		z.pred <- cbind( z.data, z.temp )
	# 		z.ranef <- t(f.eff[m.covs.r,"Estimate"] + t(r.eff[z.samps,]))
	# 		if ( PBO=="TRT" ) { z.ranef[,"DRUG"] <- z.ranef[,"DRUG"] + z.ranef[,"TRT"]}

	# 		# print("Plotting #R5 - Predicted Values")
	# 		## FCT: Plot Patient Profile
	# 		sample <- z.samps[1]
	# 		PLOT_IND <- function( sample, which_plots ) {
	# 			## Plot 1: Real vs Predicted Values
	# 			if ( 1 %in% which_plots ) {
	# 				COLS.tru <- COLS.list.2[c(6,2,1)]
	# 				COLS.conf <- COLS.list.2[4]
	# 				sub.tab <- z.pred[ z.pred[,r.grp]==sample, ]
	# 				plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main=sample )
	# 				abline( h=0:10,lty=3,col="grey50",lwd=1)
	# 				points( Estimate ~ WK, data=sub.tab,type="l",lwd=4 )
	# 				polygon( c(sub.tab$WK,rev(sub.tab$WK)), c(sub.tab[,"2.5%ile"],rev(sub.tab[,"97.5%ile"])), col=adjustcolor(COLS.conf,alpha=.2),border=NA )
	# 				points( DAS ~ WK, data=sub.tab,type="p",lwd=3,cex=1.5,pch=c(4,16)[factor(DRUG)],col=COLS.tru[1+rowSums(sub.tab[,c(PBO,"DRUG")])] )
	# 				# pre.post <- cumsum(unlist( FUL[FUL$ID==sample,c("DAS_BL_MN","DEL_MNe_MN")] ))
	# 				# arrows( c(0,24),pre.post,c(24,100),pre.post, col=COLS.pred[5],lwd=4,lty=2,length=0 )
	# 				pre.post.2 <- z.ranef[sample,c(1,3,2)] + c(0,rep(z.ranef[sample,1],2))
	# 				arrows( c(-5,2,24),pre.post.2,c(100,24,100),pre.post.2, col=COLS.tru,lwd=6,lty=2,length=0 )
	# 				text( rep(100,3),c(9,8.5,8), paste(colnames(z.ranef),round(z.ranef[sample,],2),sep="=")[c(1,3,2)], pos=2,col=COLS.tru )
	# 			}

	# 			## Plot 2: Scatter of Individual Monte Carlo Estimates
	# 			if ( 2 %in% which_plots ) {
	# 				x.name <- "Intercept"
	# 				x.col.1 <- d.post.f.which[ d.post.f.which==paste("b",x.name,sep="_") ] # grep(x.name,d.post.f.which,value=T)
	# 				x.col.2 <- grep( paste(sample,",",x.name,"]",sep=""),d.post.r2.which,value=T)
	# 				x.temp <- d.post[,x.col.1] + d.post[,x.col.2]
	# 				y.name <- "DRUG"
	# 				y.col.1 <- d.post.f.which[ d.post.f.which==paste("b",y.name,sep="_") ] # grep(y.name,d.post.f.which,value=T)
	# 				y.col.2 <- grep( paste(sample,",",y.name,"]",sep=""),d.post.r2.which,value=T)
	# 				y.temp <- d.post[,y.col.1] + d.post[,y.col.2]
	# 				w.name <- PBO
	# 				w.col.1 <- d.post.f.which[ d.post.f.which==paste("b",w.name,sep="_") ] # grep(w.name,d.post.f.which,value=T)
	# 				w.col.2 <- grep( paste(sample,",",w.name,"]",sep=""),d.post.r2.which,value=T)
	# 				w.temp <- d.post[,w.col.1] + d.post[,w.col.2]
	# 				if ( PBO=="TRT" ) { y.temp <- w.temp + y.temp }
	# 				XLIM <- c(1,8) # range(x.temp)
	# 				YLIM <- c(-3,2) # range(y.temp)
	# 				# Plot it
	# 				plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Intercept",ylab="Drug",main="Individual Simulation Estimates")
	# 				abline(h=-10:10,v=-10:10,lty=3,col="grey50",lwd=1 )
	# 				abline(h=0,v=5)
	# 				SCRAP <- points(x.temp,y.temp,col=adjustcolor(COLS.tru[3],.03),pch=16)
	# 				SCRAP <- points(x.temp,w.temp,col=adjustcolor(COLS.tru[2],.03),pch=16)
	# 				points( mean(x.temp),mean(y.temp),pch="+",col=COLS.tru[3],cex=2,lwd=2 )
	# 				points( mean(x.temp),mean(w.temp),pch="+",col=COLS.tru[2],cex=2,lwd=2 )
	# 				legend( "topleft",col=COLS.tru[2:3],pch=16,legend=c("TRT","DRUG+TRT"),title="Effect Sizes" )
	# 			}					
	# 		}

	# 		for ( samp in z.samps ) {
	# 		# for ( samp in z.samps[seq(1,20,4)] ) {
	# 			png( paste(PathToPlot,"ModSumm_",tag,".Ri.",samp,"-Pred.png",sep=""),height=1000,width=2000,pointsize=30)
	# 			par(mfrow=c(1,2))
	# 			SCRAP <- PLOT_IND( samp, 1:2 )
	# 			dev.off()
	# 		}

	# 	}
	# } ## End Random Effects Plots

	## Return some info about Models
	# COMPILE <- list( )
	# return(COMPILE)
} ## Close "PLOT_MOD" Function













## FCT: Plot Distributions of Individual Estimates
Plot_Viol <- function( var, params ) {
	## Plotting Params
	MAIN <- paste(params$ylab,"by Trial Arm")
	XLAB <- "Arm"
	YLAB <- params$ylab
	XLIM <- c(.5,3.5)
	YLIM <- params$ylim
	# Plot/Lines
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab=XLAB,ylab=YLAB,main=MAIN,xaxt="n" )
	axis( 1, at=1:3,label=c("G","P","PE") )
	abline( h=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
	abline( h=0 )
	vioplot( EST.2[EST.2$GRP=="G",var], at=1,col=adjustcolor(COLS.temp["G"],.8),border=NA,add=T )
	vioplot( EST.2[EST.2$GRP=="P",var], at=2,col=adjustcolor(COLS.temp["P"],.8),border=NA,add=T )
	vioplot( EST.2[EST.2$GRP=="PE",var], at=3,col=adjustcolor(COLS.temp["PE"],.8),border=NA,add=T )
}
## Save Plot of Distributions
png( paste(PathToPlot,"Ind_3-VioPlots.png",sep=""),height=800,width=1600,pointsize=30 )
par(mfrow=c(1,3))
SCRAP <- lapply( Vars, function(x)Plot_Viol( x, PARAMS[[x]] ) )
dev.off()












RES.9c <- resid( JJ$Null_Tc$m9i1 )
RES.9 <- resid( JJ$Null_T$m9i1 )

plot( RES.9[,1] ~ TAB$WK[!is.na(TAB$DAS)],pch=16,col=adjustcolor("blue2",.4) )
par(mfrow=c(1,2))
boxplot( RES.9[,1] ~ TAB$WK[!is.na(TAB$DAS)],col=adjustcolor("blue2",.4) ) ; abline(h=-10:10,lty=3)
boxplot( RES.9c[,1] ~ TAB$WK[!is.na(TAB$DAS)],col=adjustcolor("blue2",.4) ) ; abline(h=-10:10,lty=3)

