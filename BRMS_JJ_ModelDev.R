## Validate use of Bayes LMM Models on Time-Series Data ##
## Using "brms" Package
## Janssen Data using Repeated-Measures ##
## May 13, 2016 ##
## Kristopher Standish ##

###################################################
## Brainstorming Ideas ############################

## Validate Robustness of use of Priors
 # Vary sample size to see how Estimates are affected
   # Test increasing sample sizes with convservative priors (assume no drug effect)
 # Use simple/small models (n=50?)
   # Test many different means and sds for priors
   # time(m1)=30sec
 # Use complex/small model (n=50)
   # Test a few different means for priors
   # time(m8)=216sec
   # time(m9)=252sec
   # time(m10)=264sec
 # Use simple/large model (n=250+)
   # Test a few different means for priors
   # time(m1)=98sec
   # time(m2)=156sec

## Determine Best Null Model for Clinical Data
 # Test different fixed effects
 # Test different random effects
 # Test correlation structure(s)

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

## Parse Command Line
LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("Null",421,"Null_All","m1,m2,m3,m4,m5,m6,m7,m8,m9,m10")
# LINE <- c("HLA",421,"HLA_TEST_m1","XXXX")
Goal <- LINE[1]
N.Samps <- LINE[2]
Dir.Tag <- LINE[3]
Mod.Names.1 <- LINE[4]
Mod.Names <- unlist(strsplit(Mod.Names.1,",")[[1]])
 # Print Inputs
print(paste("Goal:",Goal))
print(N.Samps)
print(Dir.Tag)
print(Mod.Names)

## Specify Phenos/Covs
PHENOS <- c("DEL_MNe_MN","DEL_lCRP_MNe_MN","DEL_rSJC_MNe_MN","DEL_rTJC_MNe_MN",
	"DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA","RF" )
COVS <- c("RF_ACPA","RF_ACPA","RF_ACPA","RF_ACPA",
	"RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)",
	"","")
COVS <- c("","","","",
	"RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)",
	"","")
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)",
	"","")
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)","ACPA+log(DIS_DUR)",
	"","")
PH.COV.1 <- data.frame(PHENOS,COVS)
 # Look at Fewer Phenotypes
PH.COV <- PH.COV.1[c(1,5,9),]
PHENOS <- as.character(PH.COV[,1])
COVS <- as.character(PH.COV[,2])

## COLORS: Set Color Scheme
 # FCT: Blend 2 Colors
BLEND <- function( colors ) { colorRampPalette(colors)(3)[2] }
 # Set Color Palette
COLS.list <- c("firebrick2","chocolate2","gold1","chartreuse2","cadetblue2","dodgerblue2","slateblue2","magenta2")
COLS.list.2 <- c("aquamarine3","deeppink2",BLEND(c("gold1","chartreuse2")),BLEND(c("steelblue3","slateblue3")),"tomato2","deepskyblue1",BLEND(c("sienna2","goldenrod1")) )
COLS.list.ord <- COLS.list.2[c(5,7,3,1,6,4,2)]
## COLORS: Pick Specific Colors
 # Multiple Hypothis Correction Color
COLS.cor <- COLS.list[1]
 # Platform Colors
PLATS <- c("SOP","CHP","SEQ","LAB")
COLS.plat <- COLS.list.2[1:length(PLATS)]
names(COLS.plat) <- PLATS
 # Haplotype/AA Frequency Colors
COLS.fr <- COLS.list.2[7]
COLS.AA <- c(colorRampPalette(COLS.list)(26),"black","grey90","grey50") ; names(COLS.AA) <- c(LETTERS,"*",".","?")
 # Phenotype Colors
COLS.ph <- COLS.list.2[c(6,3,2)]
names(COLS.ph) <- PHENOS
 # Null LMM Model Colors
COLS.null <- COLS.list.2[c(1,4,5)]
COLS.beta <- COLS.list.2[5]

##############################################################
## LOAD DATA #################################################
##############################################################

## Load Packages
library(nlme)
library(gplots)
library(brms)
library(vioplot)

## Mac or TSCC?
Mac_Root <- "/Users/kstandis/Data/"
TSCC_Root <- "/projects/janssen/Mac_Data/"
if ( file.exists(Mac_Root) ) {
	Root <- Mac_Root
}else{
	Root <- TSCC_Root
}

## Set Date
DATE <- gsub("-","",Sys.Date())
if ( is.na(Dir.Tag) ) { Dir.Tag <- "BayesLMM" }

## New Mac Paths
PathToTypes <- paste(Root,"Janssen/Data/HLA/SOAP_HLA_Types/20151211_HLA_Types.Rdata",sep="")
PathToAA <- paste(Root,"Janssen/Data/HLA/Amino_Acids/20160126_HLA_AA.Rdata",sep="")
PathToRawFiles <- paste(Root,"Janssen/Data/Pheno/Raw_Files/",sep="")
PathToAssoc <- paste(Root,"Janssen/Data/HLA/Association/20160614_HLA_Assoc_",sep="")
PathToData <- paste(Root,"Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt",sep="")
PathToFT <- paste(Root,"Janssen/Data/Pheno/Derived/20151015_Full_Table.txt",sep="")
PathToDer <- paste(Root,"Janssen/Data/Pheno/Derived/20150619_Derived_Pheno_Table.txt",sep="")
PathToNullMod <- paste(Root,"Janssen/Plots_Mac/20160517_Null_All/",sep="")
# PathToNullMod.T <- paste(Root,"Janssen/Plots_Mac/20160523_Null_TRT_421/",sep="")
PathToNullMod.T <- paste(Root,"Janssen/Plots_Mac/20160601_Null_TRT_421/",sep="")
PathToNullMod.Tc <- paste(Root,"Janssen/Plots_Mac/20160601_Null_TRTcor_421/",sep="")
PathToDownMod.1 <- paste(Root,"Janssen/Plots_Mac/20160614_DownSample/",sep="")
# PathToHLAMod.1 <- paste(Root,"Janssen/Plots_Mac/20160519_HLA_Fr20/",sep="")
PathToHLAMod.1 <- paste(Root,"Janssen/Plots_Mac/20160614_HLA_HapTyp/",sep="")
PathToSNPMod.1 <- paste(Root,"Janssen/Plots_Mac/20160813_Cand_SNP/",sep="")
PathToPlot <- paste(Root,"Janssen/Plots_Mac/",DATE,"_",Dir.Tag,"/",sep="")
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

## Write Update to File
write(paste(date(),"#### Goal:",Goal,"####"), paste(PathToPlot,"Update.txt",sep=""),append=F)
write(paste(date(),"- Data Loaded ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)
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

##############################################################
## CREATE FUNCTIONS FOR LATER USE ############################
##############################################################
write(paste(date(),"- Creating Functions ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

PLOT_FIXED <- function( model, tag ) {

	## Collect General Model Info
	summ <- summary(model,waic=F)
	m.obs <- summ$nobs
	m.iter <- summ$iter
	m.warm <- summ$warmup
	m.chain <- summ$chains
	# m.waic <- summ$WAIC
	m.pheno <- as.character(model$formula)[2]
	m.covs <- as.character(model$formula)[3]
	m.form <- paste( m.pheno, "~", m.covs )
	RAND <- length(summ$random)>0

	## Collect Model Outputs
	d.prior <- model$prior
	d.post <- posterior_samples(model)
	f.eff <- summ$fixed
	m.covs.f <- rownames(f.eff)
	c.eff <- summ$cor_pars
	s.eff <- summ$spec_pars
	d.post.f.which <- paste("b_",m.covs.f,sep="")
	d.post.c.which <- rownames(c.eff)
	if ( RAND==T ) {
		r.grp <- summ$group
		r.ngrps <- summ$ngrps
		r.eff <- ranef(model)[[1]]
		r.eff.2 <- summ$random
		m.covs.r <- colnames(r.eff)[1:ncol(r.eff)]
		m.samps <- rownames(r.eff)
		n.samps <- length(m.samps)
		all.eff <- list( f.eff, c.eff, Reduce( rbind, r.eff.2 ), s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","R","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		mod.fit$Eff[grep("cor(",rownames(mod.fit),fixed=T)] <- "RC"
		d.post.r.which <- unlist(lapply( m.covs.r, function(x)paste( "sd_",r.grp,"_",x,sep="") ))
		d.post.r2.which <- unlist(lapply( m.covs.r, function(x)paste( "r_",r.grp,"[",m.samps,",",x,"]",sep="") ))
	}else{
		all.eff <- list( f.eff, c.eff, s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
	}
	print("Model Parsed")
	## Compile List & Save
	# out.m <- list( N.Obs=m.obs, N.iter=m.iter, N.warm=m.warm, N.chain=m.chain,
	# 	Pheno=m.pheno, Covs=m.covs, Formula=m.form, Covs.F=m.covs.f )
	# if ( RAND==T ) { out.r <- list( Grp=r.grp, N.grps=r.ngrps, Effects=r.eff, Effects.2=r.eff.2, Samps=m.samps ) }else{ out.r <- "No Random Effects"}
	# out.compile <- list( Model=out.m, Effects=mod.fit, Posterior=d.post, Prior=d.prior, Random=out.r )
	# save( out.compile, file=paste(PathToPlot,"Rdata.Model_Summary.",tag,".Rdata",sep="") )
	# print("Model Summary Saved")

	###################################
	## FIXED EFFECTS PLOTS ############
	 # Set Color Palette
	COLS.list.heat <- c("firebrick3","chocolate2","gold1","springgreen1","steelblue2","slateblue3")
	COLS <- adjustcolor(COLS.list.2,alpha=.6)
	COLS.eff <- COLS[c(1:4,6)]
	names(COLS.eff) <- c("F","C","R","RC","S")

	## Plot Model Summary
	print("Plotting #1 - Model Summary")
	png( paste(PathToPlot,"ModSumm_",tag,".1-PlotFct.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	par(mfrow=c(nrow(mod.fit),2))
	plot( model, N=nrow(mod.fit) )
	dev.off()

	# ## Plot Chains
	#  # (aka) How to Pull all Sampling Data
	# which_vars <- m.covs.f
	# which_vars <- grep( "sd(",rownames(mod.fit),fixed=T,value=T )
	# png( paste(PathToPlot,"ModSumm_",tag,".1-Chains.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	# par( mfrow=c(length(which_vars),1) )
	# for ( v in which_vars ) {
	# 	v.tag <- grep( substr(gsub("(",paste("_",r.grp,"_",sep=""),v,fixed=T),1,10), model$fit@sim$fnames_oi )
	# 	YLIM <- Reduce(range,lapply(model$fit@sim$samples,function(x)range(x[[paste("b",v,sep="_")]])))
	# 	plot( 0,0,type="n",xlim=c(0,m.iter),ylim=YLIM,xlab="Iteration",ylab=v,main="Chains" )
	# 	abline( h=-100:100,lty=3,col="grey50",lwd=1 )
	# 	# SCRAP <- lapply( seq(YLIM[1],YLIM[2],.025),function(x)abline(h=x,col=adjustcolor(COLS.list.2[4],alpha=2*dnorm(x,5,1)),lwd=2 ))
	# 	for ( c in 1:m.chain ) {
	# 		# TEMP <- model$fit@sim$samples[[c]]$b_Intercept
	# 		points( 1:m.iter, model$fit@sim$samples[[c]][[paste("b",v,sep="_")]], type="l",col=adjustcolor(COLS.list.2[c],alpha=.7),lwd=2 )
	# 	}
	# }
	# dev.off()
	
	## Fixed Effect Sizes
	print("Plotting #2 - Effect Sizes")
	YLIM <- extendrange(mod.fit[,"Estimate"], f=.2)
	png( paste(PathToPlot,"ModSumm_",tag,".2-EffSize.png",sep=""),height=1000,width=400+100*nrow(mod.fit),pointsize=26)
	par(mar=c( 7,5,5,3 ))
	# TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n",ylim=YLIM,ylab="Effect Size",main="Effect Size Estimates" )
	TEMP <- 1:nrow(mod.fit)
	plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=range(TEMP),ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
	axis( 1, at=TEMP,label=rownames(mod.fit), las=2 )
	axis( 2, at=seq(-10,10,2), las=2 )
	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
	abline( h=0, lty=1,col="black",lwd=1 )
	 # Plot Prior Distributions
	for ( v in 1:nrow(mod.fit) ) {
		var <- rownames(mod.fit)[v]
		if ( var %in% m.covs.f ) {
			if ( var=="Intercept" ) {
				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
			}else{
				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
			}
			if ( length(temp.priors)==2 ) {
				names(temp.priors) <- c("Mean","SD")
				vioplot( rnorm(1e5,temp.priors["Mean"],temp.priors["SD"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
			}
		}
		# if ( var %in% paste("sd(",m.covs.r,")",sep="") ) {
		# 	temp.priors <- as.numeric(strsplit( gsub(")","",gsub("cauchy(","",d.prior[which(d.prior$class=="sd"&d.prior$group==r.grp)[1],"prior"], fixed=T),fixed=T),"," )[[1]])
		# 	names(temp.priors) <- c("Loc","Scale")
		# 	vioplot( rnorm(1e5,temp.priors["Loc"],temp.priors["Scale"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )
		# }
	}
	 # Plot Posterior Distributions
	for ( v in 1:nrow(mod.fit) ) {
		var <- rownames(mod.fit)[v]
		var.tag <- var
		if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
		if ( grepl("sd(",var,fixed=T) ) { var.tag <- gsub(")","", gsub("sd(",paste("sd_",r.grp,"_",sep=""),var,fixed=T),fixed=T) }
		if ( grepl("cor(",var,fixed=T) ) { var.tag <- gsub(",","_", gsub(")","", gsub("cor(",paste("cor_",r.grp,"_",sep=""),var,fixed=T),fixed=T),fixed=T) }
		if ( var==paste("sigma(",m.pheno,")",sep="") ) { var.tag <- paste("sigma",m.pheno,sep="_") }
		if ( var.tag %in% colnames(d.post) ) {
			vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.eff[mod.fit[v,"Eff"]],add=T,drawRect=F )
		}
	}
	# TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n", add=T)
	arrows( TEMP,mod.fit[,"l.95..CI"],TEMP,mod.fit[,"u.95..CI"],lwd=3,length=0 )
	arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],lwd=3,length=0 )
	# legend("topright",fill=COLS.eff,border=NA,legend=names(COLS.eff),ncol=length(COLS.eff),title="Effect Type",bg="white")
	legend("topright",fill=c(adjustcolor("black",alpha=.2),COLS.eff),border=NA,legend=c("Prior",names(COLS.eff)),ncol=(1+length(COLS.eff))/2,title="Effect Type",bg="white")
	dev.off()
}

## FCT: Plot Model
PLOT_MOD <- function( model, tag, plot_rand ) {
	write(paste(date(),"Model:",tag), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Collect General Model Info
	summ <- summary(model,waic=F)
	m.obs <- summ$nobs
	m.iter <- summ$iter
	m.warm <- summ$warmup
	m.chain <- summ$chains
	# m.waic <- summ$WAIC
	m.pheno <- as.character(model$formula)[2]
	m.covs <- as.character(model$formula)[3]
	m.form <- paste( m.pheno, "~", m.covs )
	RAND <- length(summ$random)>0

	## Collect Model Outputs
	d.prior <- model$prior
	d.post <- posterior_samples(model)
	f.eff <- summ$fixed
	m.covs.f <- rownames(f.eff)
	c.eff <- summ$cor_pars
	s.eff <- summ$spec_pars
	d.post.f.which <- paste("b_",m.covs.f,sep="")
	d.post.c.which <- rownames(c.eff)
	if ( RAND==T ) {
		r.grp <- summ$group
		r.ngrps <- summ$ngrps
		r.eff <- ranef(model)[[1]]
		r.eff.2 <- summ$random
		m.covs.r <- colnames(r.eff)[1:ncol(r.eff)]
		m.samps <- rownames(r.eff)
		n.samps <- length(m.samps)
		all.eff <- list( f.eff, c.eff, Reduce( rbind, r.eff.2 ), s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","R","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		mod.fit$Eff[grep("cor(",rownames(mod.fit),fixed=T)] <- "RC"
		d.post.r.which <- unlist(lapply( m.covs.r, function(x)paste( "sd_",r.grp,"_",x,sep="") ))
		d.post.r2.which <- unlist(lapply( m.covs.r, function(x)paste( "r_",r.grp,"[",m.samps,",",x,"]",sep="") ))
	}else{
		all.eff <- list( f.eff, c.eff, s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
	}
	print("Model Parsed")
	## Compile List & Save
	# out.m <- list( N.Obs=m.obs, N.iter=m.iter, N.warm=m.warm, N.chain=m.chain,
	# 	Pheno=m.pheno, Covs=m.covs, Formula=m.form, Covs.F=m.covs.f )
	# if ( RAND==T ) { out.r <- list( Grp=r.grp, N.grps=r.ngrps, Effects=r.eff, Effects.2=r.eff.2, Samps=m.samps ) }else{ out.r <- "No Random Effects"}
	# out.compile <- list( Model=out.m, Effects=mod.fit, Posterior=d.post, Prior=d.prior, Random=out.r )
	# save( out.compile, file=paste(PathToPlot,"Rdata.Model_Summary.",tag,".Rdata",sep="") )
	# print("Model Summary Saved")

	###################################
	## FIXED EFFECTS PLOTS ############
	 # Set Color Palette
	COLS.list.heat <- c("firebrick3","chocolate2","gold1","springgreen1","steelblue2","slateblue3")
	COLS <- adjustcolor(COLS.list.2,alpha=.6)
	COLS.eff <- COLS[c(1:4,6)]
	names(COLS.eff) <- c("F","C","R","RC","S")

	## Plot Model Summary
	print("Plotting #1 - Model Summary")
	png( paste(PathToPlot,"ModSumm_",tag,".1-PlotFct.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	par(mfrow=c(nrow(mod.fit),2))
	plot( model, N=nrow(mod.fit) )
	dev.off()

	# ## Plot Chains
	#  # (aka) How to Pull all Sampling Data
	# which_vars <- m.covs.f
	# which_vars <- grep( "sd(",rownames(mod.fit),fixed=T,value=T )
	# png( paste(PathToPlot,"ModSumm_",tag,".1-Chains.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	# par( mfrow=c(length(which_vars),1) )
	# for ( v in which_vars ) {
	# 	v.tag <- grep( substr(gsub("(",paste("_",r.grp,"_",sep=""),v,fixed=T),1,10), model$fit@sim$fnames_oi )
	# 	YLIM <- Reduce(range,lapply(model$fit@sim$samples,function(x)range(x[[paste("b",v,sep="_")]])))
	# 	plot( 0,0,type="n",xlim=c(0,m.iter),ylim=YLIM,xlab="Iteration",ylab=v,main="Chains" )
	# 	abline( h=-100:100,lty=3,col="grey50",lwd=1 )
	# 	# SCRAP <- lapply( seq(YLIM[1],YLIM[2],.025),function(x)abline(h=x,col=adjustcolor(COLS.list.2[4],alpha=2*dnorm(x,5,1)),lwd=2 ))
	# 	for ( c in 1:m.chain ) {
	# 		# TEMP <- model$fit@sim$samples[[c]]$b_Intercept
	# 		points( 1:m.iter, model$fit@sim$samples[[c]][[paste("b",v,sep="_")]], type="l",col=adjustcolor(COLS.list.2[c],alpha=.7),lwd=2 )
	# 	}
	# }
	# dev.off()
	
	## Fixed Effect Sizes
	print("Plotting #2 - Effect Sizes")
	YLIM <- extendrange(mod.fit[,"Estimate"], f=.2)
	png( paste(PathToPlot,"ModSumm_",tag,".2-EffSize.png",sep=""),height=1000,width=400+100*nrow(mod.fit),pointsize=26)
	par(mar=c( 7,5,5,3 ))
	# TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n",ylim=YLIM,ylab="Effect Size",main="Effect Size Estimates" )
	TEMP <- 1:nrow(mod.fit)
	plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=range(TEMP),ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
	axis( 1, at=TEMP,label=rownames(mod.fit), las=2 )
	axis( 2, at=seq(-10,10,2), las=2 )
	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
	abline( h=0, lty=1,col="black",lwd=1 )
	 # Plot Prior Distributions
	for ( v in 1:nrow(mod.fit) ) {
		var <- rownames(mod.fit)[v]
		if ( var %in% m.covs.f ) {
			if ( var=="Intercept" ) {
				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
			}else{
				temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
			}
			if ( length(temp.priors)==2 ) {
				names(temp.priors) <- c("Mean","SD")
				vioplot( rnorm(1e5,temp.priors["Mean"],temp.priors["SD"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
			}
		}
		# if ( var %in% paste("sd(",m.covs.r,")",sep="") ) {
		# 	temp.priors <- as.numeric(strsplit( gsub(")","",gsub("cauchy(","",d.prior[which(d.prior$class=="sd"&d.prior$group==r.grp)[1],"prior"], fixed=T),fixed=T),"," )[[1]])
		# 	names(temp.priors) <- c("Loc","Scale")
		# 	vioplot( rnorm(1e5,temp.priors["Loc"],temp.priors["Scale"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )
		# }
	}
	 # Plot Posterior Distributions
	for ( v in 1:nrow(mod.fit) ) {
		var <- rownames(mod.fit)[v]
		var.tag <- var
		if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
		if ( grepl("sd(",var,fixed=T) ) { var.tag <- gsub(")","", gsub("sd(",paste("sd_",r.grp,"_",sep=""),var,fixed=T),fixed=T) }
		if ( grepl("cor(",var,fixed=T) ) { var.tag <- gsub(",","_", gsub(")","", gsub("cor(",paste("cor_",r.grp,"_",sep=""),var,fixed=T),fixed=T),fixed=T) }
		if ( var==paste("sigma(",m.pheno,")",sep="") ) { var.tag <- paste("sigma",m.pheno,sep="_") }
		if ( var.tag %in% colnames(d.post) ) {
			vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.eff[mod.fit[v,"Eff"]],add=T,drawRect=F )
		}
	}
	# TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n", add=T)
	arrows( TEMP,mod.fit[,"l.95..CI"],TEMP,mod.fit[,"u.95..CI"],lwd=3,length=0 )
	arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit[,"Estimate"],lwd=3,length=0 )
	# legend("topright",fill=COLS.eff,border=NA,legend=names(COLS.eff),ncol=length(COLS.eff),title="Effect Type",bg="white")
	legend("topright",fill=c(adjustcolor("black",alpha=.2),COLS.eff),border=NA,legend=c("Prior",names(COLS.eff)),ncol=(1+length(COLS.eff))/2,title="Effect Type",bg="white")
	dev.off()

	## RANDOM EFFECTS PLOTS ###########
	if ( RAND==T & plot_rand!=F ) {
		print("Plotting Random Effects")
		n.covs.r <- length(m.covs.r)

		## Plot Distributions & Shrinkage of Random Effects
		print("Plotting #R1 - Distributions")
		PLOT_SHRINK <- function( var, which_plots ) {
			## Pull Data
			brm.var <- r.eff[,var] + f.eff[var,"Estimate"]
			if ( var %in% c("Intercept","DRUG") ) {
				if ( var=="Intercept" ) {
					mn.var <- FUL[,"DAS_BL_MN"]
					brm.fx <- f.eff[var,"Estimate"]
				}
				if ( var=="DRUG" ) {
					mn.var <- FUL[,"DEL_MNe_MN"]
					brm.fx <- f.eff[var,"Estimate"] + f.eff["TRT","Estimate"]
					if ( "TRT"%in%m.covs.f ) { brm.var <- brm.var + f.eff["TRT","Estimate"] }
				}
				mn.mod <- T
				names(mn.var) <- as.character(FUL$ID)
			}else{
				mn.var <- brm.var
				mn.mod <- F
				brm.fx <- f.eff[var,"Estimate"]
			}
			mg.var <- merge( data.frame(brm.var),data.frame(mn.var), by="row.names" )
			mg.var <- merge( FUL[,c("ID","GRP")], mg.var, by.x="ID",by.y="Row.names" )
			mg.var.ord <- mg.var[ order(mg.var$mn.var), ]
			mn.fx <- mean(mg.var$mn.var)

			## Plot 1: Distributions
			if ( 1 %in% which_plots ) {
				XLIM <- range( c(mg.var$mn.var,mg.var$brm.var) )
				X_BIN <- .25
				BRKS <- seq( floor(XLIM[1]), ceiling(XLIM[2])+X_BIN, X_BIN)
				YLIM <- c( 0, Reduce( max, lapply(list(mg.var$brm.var,mg.var$mn.var),function(x)hist(x,breaks=BRKS,plot=F)$counts) ) )
				hist( mg.var$brm.var, col=COLS[1],border=NA,breaks=BRKS, xlim=XLIM,ylim=YLIM, main=paste(var,"- Mean vs BRMS"),xlab=var,ylab="# Patients" )
				abline( v=brm.fx, col=COLS.list.2[1],lwd=4,lty=1 )
				if ( mn.mod==T ) {
					hist( mg.var$mn.var, col=COLS[2],border=NA,breaks=BRKS, add=T )
					abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=1 )
					legend( "topleft", legend=c("BRMS","Mean"),title="Model",fill=COLS[1:2],border=NA)	
				}			
			}
	
			## Plot 2: Boxplot by Arm
			if ( 2 %in% which_plots ) {
				ARMS <- as.character( unique( mg.var$GRP ) )
				XLIM <- c(0,7)
				YLIM <- range( c(mg.var$mn.var,mg.var$brm.var) )
				plot( 0,0,type="n", xlim=XLIM,ylim=YLIM,xaxt="n",xlab="ARM",ylab=var,main=paste(var,"Distribution by Arm") )
				abline( h=-10:10,lty=3,col="grey50",lwd=1 )
				axis( 1,at=seq(1.5,7,2),label=ARMS )
				SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$brm.var[mg.var$GRP==ARMS[x]],at=2*x-1,col=COLS[1],border=NA,add=T) )
				abline( h=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
				if ( mn.mod==T ) {
					SCRAP <- lapply( 1:3, function(x)vioplot(mg.var$mn.var[mg.var$GRP==ARMS[x]],at=2*x,col=COLS[2],border=NA,add=T) )
					abline( h=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
					legend( "topleft", legend=c("BRMS","Mean"),title="Model",fill=COLS[1:2],border=NA)	
					arrows( c(1,3,5),rep(YLIM[1],3),c(2,4,6),rep(YLIM[1],3), length=0,lwd=4,col=COLS.list.2[3:5] )
				}		
			}

			## Plot 3: Shrinkage by Person
			if ( 3 %in% which_plots ) {
				XLIM <- range( c(mg.var.ord$mn.var,mg.var.ord$brm.var) )
				X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5)
				YLIM <- c( 1, n.samps )
				plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab=var,ylab="Patient",main=paste(var,"- Estimate Shrinkage") )
				abline( v=X_LNS, lty=3,col="grey50",lwd=1 )
				points( mg.var.ord$brm.var,1:nrow(mg.var.ord), pch=16,cex=1,col=COLS[1] )
				if ( mn.mod==T ) {
					points( mg.var.ord$mn.var,1:nrow(mg.var.ord), pch=16,cex=1,col=COLS[2] )
					arrows( mg.var.ord$mn.var,1:nrow(mg.var.ord),mg.var.ord$brm.var,1:nrow(mg.var.ord), angle=30,length=.2,lwd=2,col=COLS.list.2[3:5][factor(mg.var.ord$GRP)] )
					abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
				}
				abline( v=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
				legend("topleft",col=COLS.list.2[1:5],legend=c("Mod: BRMS","Mod: Mean","Arm: G","Arm: P","Arm: PE"),pch=c(16,16,NA,NA,NA),lty=c(0,0,1,1,1),lwd=2 )
			}

			## Plot 4: MN vs BRMS Model Estimates
			if ( 4 %in% which_plots ) {
				XLIM <- range( mg.var$mn.var )
				YLIM <- range( mg.var$brm.var )
				X_LNS <- seq( floor(XLIM[1]-5), ceiling(XLIM[2]+5), .5 )
				Y_LNS <- seq( floor(YLIM[1]-5), ceiling(YLIM[2]+5), .5 )
				plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Mean Estimate",ylab="BRMS Estimate",main=var )
				abline( v=X_LNS,h=Y_LNS, lty=3,col="grey50",lwd=1 )
				abline( 0,1 )
				abline( h=brm.fx, col=COLS.list.2[1],lwd=4,lty=2 )
				abline( v=mn.fx, col=COLS.list.2[2],lwd=4,lty=2 )
				points( mg.var$brm.var ~ mg.var$mn.var, data=mg.var, pch=16,cex=1,col=COLS[3:5][factor(GRP)] )
				legend( "topleft",pch=16,col=COLS[3:5],legend=levels(factor(mg.var$GRP)))
			}
		}
		 # Plot it
		for ( var in m.covs.r ) {
			if ( var %in% c("Intercept","DRUG") ) {
				png( paste(PathToPlot,"ModSumm_",tag,".R1-Dist.",var,".png",sep=""),height=2000,width=2000,pointsize=32)
				layout( matrix(c(1:4,4,4),byrow=F,ncol=2) )
				PLOT_SHRINK( var, c(1,2,4) )
				PLOT_SHRINK( var, c(3) )
				dev.off()
			}else{
				png( paste(PathToPlot,"ModSumm_",tag,".R1-Dist.",var,".png",sep=""),height=2000,width=2000,pointsize=32)
				layout( matrix(c(1:2,3,3),byrow=F,ncol=2) )
				PLOT_SHRINK( var, c(1,2) )
				PLOT_SHRINK( var, c(3) )
				dev.off()
			}

		}
	
		## Plot Pairs of Random Effects
		panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...) {
			usr <- par("usr"); on.exit(par(usr))
			par(usr = c(0, 1, 0, 1))
			r <- abs(cor(x, y))
			txt <- format(c(r, 0.123456789), digits = digits)[1]
			txt <- paste0(prefix, txt)
			if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
			text(0.5, 0.5, txt, cex = cex.cor * r)
		}
		print("Plotting #R1 - Pairs")
		if ( n.covs.r > 1 ) {
			temp.reff <- t( f.eff[m.covs.r,"Estimate"] + t(r.eff) ) 
			# temp.grp <- FUL$GRP[ match(rownames(temp.reff),FUL$ID) ] # col=adjustcolor(COLS.list.2[factor(temp.grp)],alpha=.4)
			temp.reff <- merge( FUL[,c("ID","GRP")], temp.reff, by.x="ID",by.y="row.names" )
			png( paste(PathToPlot,"ModSumm_",tag,".R1-PairRand.png",sep=""),height=400*n.covs.r,width=400*n.covs.r,pointsize=24+2*n.covs.r)
			# pairs( temp.reff, pch=16,col=adjustcolor(COLS.eff["RC"],alpha=.5),cex=2,upper.panel=panel.cor )
			pairs( temp.reff[,m.covs.r], pch=16,col=COLS[3:5][factor(temp.reff$GRP)],cex=2,upper.panel=panel.cor )
			dev.off()
		}

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
			plot( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp], main=paste("Random",z.name,"by Patient"),xlab="",ylab=z.name,xaxt="n",xlim=XLIM,ylim=YLIM )
			axis( 1, at=1:ncol(z.temp),label=m.samps[z.ord],las=2,cex=.6 )
			abline( h=-10:10,v=1:ncol(z.temp), lty=3,col="grey50",lwd=1 )
			abline( h=0, lty=1,col="black",lwd=1 )
			# Text/Legend
			if ( any(z.name %in% c("DRUG","PLAC","TRT")) ) {
				abline(h=c(0,-1),lwd=2,col=COLS.list.2[6:7] )
				text( 1:ncol(z.temp), rep(YLIM[1]+.5,ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < 0)))/nrow(z.temp),1),col=COLS.list.2[6],cex=.7,srt=90 )
				text( 1:ncol(z.temp), rep(YLIM[1],ncol(z.temp)), label=round(100*apply(z.temp[,z.ord],2,function(x)length(which(x < -1)))/nrow(z.temp),1),col=COLS.list.2[7],cex=.7,srt=90 )
			}else{
				# legend( "topleft",fill=c(COLS.eff["R"],adjustcolor("black",alpha=.2)),border=NA,legend=c("Posterior","Prior"),ncol=2 )
				legend( "topleft",pch=c(rep(16,3),16,1),cex=1.5,col=c(COLS.list.2[3:5],rep("black",2)),legend=c("Arm: G","Arm: P","Arm: PE","Median","CI: 95%"),lwd=3 )
			}
			# Points/Arrows
			arrows( 1:ncol(z.temp),z.quants["25%",z.ord],1:ncol(z.temp),z.quants["75%",z.ord],col=COLS[1],lwd=5,code=3,angle=90,length=.1 )
			points( rep(1:ncol(z.temp),2),c(z.quants["2.5%",z.ord],z.quants["97.5%",z.ord]),pch=1,lwd=3,col=COLS.list.2[3:5][z.grp],cex=1.5 )
			points( 1:ncol(z.temp),z.quants["50%",z.ord],pch=16,lwd=3,col=COLS.list.2[3:5][z.grp],cex=1.5 )
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
	} ## End Random Effects Plots
} ## Close "PLOT_MOD" Function

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

## Take random sample of Patients
if (is.na(N.Samps)) { N.Samps <- 50 }
Samps <- sample( SAMPS.short, N.Samps )

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
## COMPARE PREVIOUSLY RUN MODELS ############################
#############################################################
if ( is.na(Mod.Names.1) ) { 

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

	## Load Null Models (TRT)
	if ( !exists("JJ") ) { JJ <- list() }
	load.tag <- "Null_T"
	JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.T, load.tag )
	## Load Null Models (TRT, w/ Corr)
	load.tag <- "Null_Tc"
	if ( !exists("JJ") ) { JJ <- list() }
	JJ[[load.tag]] <- LOAD_PREV_MODS( PathToNullMod.Tc, load.tag )
	## Load Down Sampled Models
	load.tag <- "down9"
	if ( !exists("JJ") ) { JJ <- list() }
	JJ[[load.tag]] <- LOAD_PREV_MODS( PathToDownMod.1, load.tag )
	load( paste(PathToDownMod.1,"Rdata.ModObs.",load.tag,".Rdata",sep="") )
	## Load HLA Models (1)
	load.tag <- "HLA"
	if ( !exists("JJ") ) { JJ <- list() }
	JJ[[load.tag]] <- LOAD_PREV_MODS( PathToHLAMod.1, load.tag )
	## Load SNP Models (1)
	load.tag <- "SNP"
	if ( !exists("JJ") ) { JJ <- list() }
	JJ[[load.tag]] <- LOAD_PREV_MODS( PathToSNPMod.1, load.tag )

	## Calculate WAIC of Models
	JJ.waic <- list()
	for ( load.tag in names(JJ) ) {
		JJ.waic[[load.tag]] <- list()
		for ( m.tag in names(JJ[[load.tag]]) ) {
			print(paste(load.tag,m.tag))
			if ( !(m.tag %in% names(JJ.waic[[load.tag]])) ) {
				JJ.waic[[load.tag]][[m.tag]] <- WAIC( JJ[[load.tag]][[m.tag]] )	
			}
		}
	}
	save( JJ.waic, file=paste(PathToPlot,"Rdata.Null_WAIC.Rdata",sep="") )

	## Plot Model Comparisons
	XLIM <- c(1,Reduce(max,lapply(JJ.waic,length)))
	YLIM <- Reduce( range, lapply(JJ.waic,function(x)lapply(x,function(y)y$waic)) )
	COLS.mod <- adjustcolor(COLS.list.2,.8)[1:length(JJ)]
	names(COLS.mod) <- names(JJ)
	png( paste(PathToPlot,"Model_Selection.WAIC.png",sep=""), height=1200,width=1600,pointsize=36 )
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Model",ylab="WAIC",main="Model Selection: WAIC" )
	abline(h=seq(0,100e3,2e3),v=0:20, lty=3,col="grey50",lwd=1 )
	for ( load.tag in names(JJ) ) {
		points( 1:length(JJ.waic[[load.tag]]), unlist(lapply(JJ.waic[[load.tag]],function(y)y$waic)), type="o",pch=16,col=COLS.mod[load.tag],lwd=2 )
	}
	legend( "topright",pch=16,lty=1,lwd=2,col=COLS.mod,legend=names(COLS.mod) )
	dev.off()

	## Plot Model Summaries
	for ( load.tag in names(JJ) ) {
		for ( m.tag in names(JJ[[load.tag]]) ) {
			file.tag <- paste(load.tag,m.tag,sep="-")
			print(file.tag)
			PLOT_MOD( JJ[[load.tag]][[m.tag]], file.tag, 1 )
		}
	}

}
write(paste(date(),"- Which Goal to Pursue??"), paste(PathToPlot,"Update.txt",sep=""),append=T)
write(paste(date(),Goal), paste(PathToPlot,"Update.txt",sep=""),append=T)
print(paste("Goal:",Goal))
##############################################################
## NULL LMM MODEL ############################################
##############################################################
if ( grepl("Null",Goal,ignore.case=T) ) {

	write(paste(date(),"- Building Null Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Set Default Tag
	tag <- "Null"
	MOD.priors <- MOD.form <- list()
	JJ.priors[[tag]] <- JJ.form[[tag]] <- list()

	## Null Model Formulae & Priors
	 # M1: Drug Only
	m.tag <- "m1"
	MOD.form[[m.tag]] <- "DAS ~ DRUG"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG )
	 # M2: Drug + Week
	m.tag <- "m2"
	MOD.form[[m.tag]] <- "DAS ~ DRUG+WK"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK )
	 # M3: Drug + Week + Placebo
	m.tag <- "m3"
	MOD.form[[m.tag]] <- "DAS ~ DRUG+WK+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC )
	 # M4: Drug*Week + Placebo
	m.tag <- "m4"
	MOD.form[[m.tag]] <- "DAS ~ DRUG*WK+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK )
	 # M5: Drug*Week + ACPA + Placebo
	m.tag <- "m5"
	MOD.form[[m.tag]] <- "DAS ~ DRUG*WK+ACPA+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA )
	 # M6: Drug*(Week+ACPA) + Placebo
	m.tag <- "m6"
	MOD.form[[m.tag]] <- "DAS ~ DRUG*(WK+ACPA)+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA )
	 # M7: Drug*(Week+ACPA) + Placebo + Random Intercept
	m.tag <- "m7"
	MOD.form[[m.tag]] <- "DAS ~ (1|ID)+DRUG*(WK+ACPA)+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd )
	 # M8: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug
	m.tag <- "m8"
	MOD.form[[m.tag]] <- "DAS ~ (1+DRUG|ID)+DRUG*(WK+ACPA)+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )
	 # M9: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac
	m.tag <- "m9"
	MOD.form[[m.tag]] <- "DAS ~ (1+DRUG+PLAC|ID)+DRUG*(WK+ACPA)+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )
	 # M10: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac & Week
	m.tag <- "m10"
	MOD.form[[m.tag]] <- "DAS ~ (1+DRUG+PLAC+WK|ID)+DRUG*(WK+ACPA)+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )
	 # M11: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Week
	m.tag <- "m11"
	MOD.form[[m.tag]] <- "DAS ~ (1+DRUG+WK|ID)+DRUG*(WK+ACPA)+PLAC"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )

	####################################
	## Specify Alternative Null Models #

	## Models w/ Remove DRUG*WK Interactions
	if ( grepl("noint",Goal,ignore.case=T) ) {
		tag <- paste(tag,"n",sep="_")
		MOD.form$m5 <- "DAS ~ DRUG+WK+ACPA+PLAC"
		MOD.form$m6 <- "DAS ~ DRUG*ACPA+WK+PLAC"
		MOD.form$m7 <- "DAS ~ (1|ID)+DRUG*ACPA+WK+PLAC"
		MOD.form$m8 <- "DAS ~ (1+DRUG|ID)+DRUG*ACPA+WK+PLAC"
		MOD.form$m9 <- "DAS ~ (1+DRUG+PLAC|ID)+DRUG*ACPA+WK+PLAC"
		MOD.form$m10 <- "DAS ~ (1+DRUG+PLAC+WK|ID)+DRUG*ACPA+WK+PLAC"
		MOD.form$m11 <- "DAS ~ (1+DRUG+WK|ID)+DRUG*ACPA+WK+PLAC"
	}

	## Models w/ TRT instead of PLAC
	if ( grepl("trt",Goal,ignore.case=T) ) {
		tag <- paste(tag,"T",sep="_")
		MOD.form <- lapply( MOD.form, function(x)gsub("PLAC","TRT",x) )
		for ( m in 1:length(MOD.priors) ) {
			MOD.priors[[m]]$coef <- gsub("PLAC","TRT",MOD.priors[[m]]$coef)
		}
	}

	## Models w/ Correlation Stucture
	if ( grepl("cor",Goal,ignore.case=T) ) {
		# Default Correlation Structure b/n Weeks
		tag <- paste(tag,"c",sep="")
		JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
		# Week within Drug Status
		if ( grepl("idr",Goal,ignore.case=T) ) {
			tag <- paste(tag,"id",sep="")
			JJ.cor[[tag]] <- cor_ar( ~WK | ID:DRUG, p=1, cov=F )
		}
		# Week within Drug Status
		if ( grepl("dri",Goal,ignore.case=T) ) {
			tag <- paste(tag,"di",sep="")
			JJ.cor[[tag]] <- cor_ar( ~WK | DRUG:ID, p=1, cov=F )
		}
			# # TEST.1 <- brm( mod.form,data=TAB.temp,, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,prior=mod.prior, autocor=JJ.cor$Null_Tcd )
			# # summary(TEST.1)$cor_pars
			# # apply( coef(TEST.1)$ID, 2, sd )
			# JJ.cor[[tag]] <- cor_ar( ~WK | DRUG:ID, p=1, cov=F )
			# # TEST.2 <- brm( mod.form,data=TAB.temp,, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,prior=mod.prior, autocor=JJ.cor$Null_Tcd )
			# # summary(TEST.2)$cor_pars
			# # apply( coef(TEST.2)$ID, 2, sd )
			# JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
			# # TEST.3 <- brm( mod.form,data=TAB.temp,, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,prior=mod.prior, autocor=JJ.cor$Null_Tcd )
			# # summary(TEST.3)$cor_pars
			# # apply( coef(TEST.3)$ID, 2, sd )
			# JJ.cor[[tag]] <- cor_ar( ~WK | DRUG, p=1, cov=F )
			# # TEST.4 <- brm( mod.form,data=TAB.temp,, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,prior=mod.prior, autocor=JJ.cor$Null_Tcd )
			# # summary(TEST.4)$cor_pars
			# # apply( coef(TEST.4)$ID, 2, sd )
		print( JJ.cor[[tag]] )
	}

	####################################
	## Run Models ######################

	## Pull/Specify Model Names/Info/Etc...
	 # Mod Names
	if ( is.na(Mod.Names.1) ) { Mod.Names <- grep("^m",names(MOD.form),value=T) }
	Mod.N <- length(Mod.Names)
	 # Priors & Formulae
	JJ.priors[[tag]] <- MOD.priors
	JJ.form[[tag]] <- MOD.form

	## Specify Iterations & Samples
	# N.Samps <- 100
	Pr.n.iter <- 1
	Pr.n.iter.tot <- Mod.N * Pr.n.iter
	JJ.samps[[tag]] <- replicate(Pr.n.iter, sample(SAMPS.short,N.Samps) )
	JJ.samps[[tag]] <- array( SAMPS.short, c(length(SAMPS.short),1) )
	save( JJ.samps, file=paste(PathToPlot,"Rdata.Samps.Rdata",sep="") )

	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	if ( tag %in% names(JJ.cor) ) { write(paste("Cor:",as.character(JJ.cor[[tag]])[[1]]), paste(PathToPlot,"Update.txt",sep=""),append=T) }
	## Run Models
	print(paste("Starting:",tag))
	print(paste("...including models:",paste(Mod.Names,collapse=",")))
	JJ[[tag]] <- list()
	for ( m in 1:Mod.N ) {
		m.tag <- Mod.Names[m]
		mod.form <- as.formula( JJ.form[[tag]][[m.tag]] )
		mod.prior <- JJ.priors[[tag]][[m.tag]]
		JJ[[tag]][[m.tag]] <- JJ.time[[tag]][[m.tag]] <- list()
		for ( i in 1:Pr.n.iter ) {
			i.tag <- paste(m.tag,"i",i,sep="")
			print(paste("Running Model:",i.tag))
			write(paste(date(),"- Running Model:",i.tag,"######"), paste(PathToPlot,"Update.txt",sep=""),append=T)
			print(mod.form)
			write(paste(date(),"-",mod.form), paste(PathToPlot,"Update.txt",sep=""),append=T)
			TAB.i <- TAB[ TAB$ID%in%JJ.samps[[tag]][,i], ]
			if ( i==1 ) {
				if ( tag %in% names(JJ.cor) ) {
					mod.cor <- JJ.cor[[tag]] ; print(mod.cor)
					JJ.time[[tag]][[m.tag]][[i.tag]] <- system.time(
						JJ[[tag]][[m.tag]][[i.tag]] <- brm( mod.form,
						data=TAB.i, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
			            prior=mod.prior, autocor=mod.cor )
					)
					i.tag.1 <- i.tag
				}else{
					JJ.time[[tag]][[m.tag]][[i.tag]] <- system.time(
						JJ[[tag]][[m.tag]][[i.tag]] <- brm( mod.form,
						data=TAB.i, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
			            prior=mod.prior )
					)
					i.tag.1 <- i.tag
				}
			}else{
				JJ.time[[tag]][[m.tag]][[i.tag]] <- system.time(
					JJ[[tag]][[m.tag]][[i.tag]] <- update( JJ[[tag]][[m.tag]][[i.tag.1]],
					newdata=TAB.i, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm ) # autocor=mod.cor
				)
			}
			write(paste(date(),"Model",i.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
			## Print Some Results
			temp.print <- paste(rownames(fixef(JJ[[tag]][[m.tag]][[i.tag]])),round(fixef(JJ[[tag]][[m.tag]][[i.tag]]),6),sep="=")
			print(paste("Model:",i.tag,"- FIXEF:", temp.print ))
			print(paste("Model:",i.tag,"- TIME:", round( JJ.time[[tag]][[m.tag]][[i.tag]][3], 3 ) ))
			## Print Some Results
			temp.model <- JJ[[tag]][[m.tag]][[i.tag]]
			save( temp.model, file=paste(PathToPlot,"Rdata.Model.",i.tag,".Rdata",sep="") )
		}
	}
	write(paste(date(),"Done Building Models!!"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Write Time Results
	save( JJ.time, file=paste(PathToPlot,"Rdata.Time.Rdata",sep="") )

}
##############################################################
## HLA LMM MODEL #############################################
##############################################################
if ( grepl("HLA",Goal,ignore.case=T) ) {
	write(paste(date(),"- Building HLA Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Set Default Tag
	tag <- "HLA"
	# MOD.priors <- MOD.form <- list()
	JJ.priors[[tag]] <- JJ.form[[tag]] <- list()

	## Set Priors for HLA Variables
	JJ.priors.list$HLA <- set_prior("normal(0,.3)",class="b",coef="HLA")
	JJ.priors.list$DRUG_HLA <- set_prior("normal(0,.3)",class="b",coef="DRUG:HLA")

	#########################################
	## HLA DATA #############################

	## Load HLA Types
	HLA_AA.l <- read.table(paste(PathToAssoc,"HLA_AA_Table.txt",sep=""),header=T,sep="\t" )
	HLA_TYP.l <- read.table(paste(PathToAssoc,"HLA_Types_Table.txt",sep=""),header=T,sep="\t" )
	HLA_HAP.l <- read.table(paste(PathToAssoc,"HLA_Hap_Table.txt",sep=""),header=T,sep="\t" )
	colnames(HLA_AA.l) <- gsub(".","_",colnames(HLA_AA.l),fixed=T)
	colnames(HLA_AA.l) <- gsub("-","_",colnames(HLA_AA.l),fixed=T)
	colnames(HLA_TYP.l) <- gsub(".","_",colnames(HLA_TYP.l),fixed=T)
	colnames(HLA_TYP.l) <- gsub("-","_",colnames(HLA_TYP.l),fixed=T)
	colnames(HLA_HAP.l) <- gsub(".","_",colnames(HLA_HAP.l),fixed=T)
	colnames(HLA_HAP.l) <- gsub("-","_",colnames(HLA_HAP.l),fixed=T)

	# ## Load Janssen HLA Results
	#  # Types
	# load( PathToTypes )
	# TYPES.l <- COMPILE
	#  # Amino Acids
	# load( PathToAA )
	# AA.l <- COMPILE

	# ## Load Betas from Mean Models
	#  # HLA Types
	# load(paste(PathToAssoc,"0-B_Precise.Rdata",sep=""))
	# summary(B.out)
	#  # Collapsed Haplotypes
	# B.out[["pSE"]] <- read.table( paste(PathToAssoc,"TAB_DRB1_pSE_Table.txt",sep=""),sep="\t",header=T )
	# B.out[["p117174"]] <- read.table( paste(PathToAssoc,"TAB_DRB1_p117174_Table.txt",sep=""),sep="\t",header=T )

	# ## Pull out Patient Data
	# PAT_DOS <- TYPES.l$GENES.2.list
	# PAT_TYP <- TYPES.l$GENES
	# PAT_AA <- AA.l$PAT_AA
	#  # Print Example
	# # PAT_AA$DRB1[X,X]

	## HLA Tables
	HLA_TYP <- HLA_TYP.l[ which(HLA_TYP.l$ID %in% SAMPS.short), ]
	HLA_AA <- HLA_AA.l[ which(HLA_AA.l$ID %in% SAMPS.short), ]
	HLA_HAP <- HLA_HAP.l[ which(HLA_HAP.l$ID %in% SAMPS.short), ]

	####################################
	## Get Models Set Up ###############
	write(paste(date(),"Setting Model Parameters"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Which HLA predictors to use?
	Mods.hla.dat <- list()
	 # DRB1 Type
	temp.types <- grep( "Pr4_DRB",colnames(HLA_TYP),value=T )
	temp.types <- grep( "__",temp.types,value=T,invert=T )	
	Mods.hla.dat$Type <- HLA_TYP[, c("ID",temp.types) ]
	 # AA
	temp.which_pos <- c(11,13,70:74)
	Mods.hla.dat$AA <- HLA_AA[, c("ID",grep(paste(paste( "Pr4_DRB1_Pos_",temp.which_pos,sep="" ),collapse="|"),colnames(HLA_AA),value=T )) ]
	 # Haps
	temp.haps.p117174 <- c("DRE","GRQ","LEA","LRA","PAA","PRA","RRA","SEA","SKR","SRA","SRE","SRL","VEA","VKA","VRA","VRE")
	Mods.hla.dat$p117174 <- HLA_HAP[, c("ID",temp.haps.p117174) ]
	temp.haps.p11137174 <- c("DFRE","GYRQ","LFEA","LFRA","PRAA","PRRA","RSRA","SGRA","SGRE","SGRL","SSEA","SSKR","SSRA","SSRE","VFRA","VHEA","VHKA","VHRA","VHRE")
	Mods.hla.dat$p11137174 <- HLA_HAP[, c("ID",temp.haps.p11137174) ]
	temp.haps.pSE <- c("DERAA","DRRAA","DRRAL","DRRGQ","QARAA","QKRAA","QKRGR","QRRAA","QRRAE","RRRAA","RRRAE")
	Mods.hla.dat$pSE <- HLA_HAP[, c("ID",temp.haps.pSE) ]
	 # HLA-A Type
	temp.types.a <- grep( "Pr2_A",colnames(HLA_TYP),value=T )
	temp.types.a <- grep( "__",temp.types.a,value=T,invert=T )	
	Mods.hla.dat$HLA_A <- HLA_TYP[, c("ID",temp.types.a) ]

	## HLA9: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac
	MOD.form <- "DAS ~ (1+DRUG+PLAC|ID)+DRUG*(WK+ACPA+HLA)+PLAC"
	MOD.priors <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )
	JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )

	## HLA1: Drug Only (Simple Model for Testing Purposes)
	# tag <- "hla1"
	# JJ.priors[[tag]] <- c(JJ.priors.list$Intercept,
	# 	JJ.priors.list$DRUG )
	# # JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
	# JJ.form[[tag]] <- "DAS ~ DRUG"

	## FCT: Modify Priors & Formula for various HLA Predictiors
	Modify_HLA_Inputs <- function( mod.name ) {
		## Pull Predictor Info
		hla.data <- Mods.hla.dat[[mod.name]]
		hla.predictors <- colnames(hla.data)[-1]
		## Modify Formula
		hla.form.preds <- paste(hla.predictors,collapse="+")
		OUT.form <- gsub( "HLA",hla.form.preds,MOD.form )
		## Modify Priors
		OUT.priors <- MOD.priors
		for ( hla.pred in hla.predictors ) {
			OUT.priors <- rbind( OUT.priors, JJ.priors.list$HLA, JJ.priors.list$DRUG_HLA )
			OUT.priors[which(OUT.priors$coef=="HLA"),"coef"] <- hla.pred
			OUT.priors[which(OUT.priors$coef=="DRUG:HLA"),"coef"] <- paste("DRUG:",hla.pred,sep="")
		}
		## Return Modified Formula/Priors
		OUT <- list( Formula=OUT.form, Priors=OUT.priors, Data=hla.data )
		return(OUT)
	}

	####################################
	## Run Models ######################

	# ## Specific HLA Alleles to Test
	# HLA.drb.all <- grep("Pr4_DRB1",colnames(HLA_TYP),value=T)
	# HLA.drb.all.freq <- apply( HLA_TYP[,HLA.drb.all], 2, sum )
	# HLA.drb <- HLA.drb.all[ HLA.drb.all.freq > 20 ]
	# HLA.N <- length(HLA.drb)
	# ## Merge Clinical Table w/ HLA Table
	# TAB.hla <- merge( TAB, HLA_TYP[,c("ID",HLA.drb)], by="ID" )

	## Specify Iterations & Samples
	# N.Samps <- 421
	JJ.samps[[tag]] <- # sample(SAMPS.short,N.Samps)
	save( JJ.samps, file=paste(PathToPlot,"Rdata.Samps.Rdata",sep="") )

	## Run Models
	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	if ( !exists("JJ") ) { JJ[[tag]] <- JJ.time[[tag]] <- list() }
	for ( m in 1:length(Mod.Names) ) {
		mod <- Mod.Names[m]
		m.tag <- paste("hla",mod,sep="_")
		mod.inputs <- Modify_HLA_Inputs( mod )
		# Specify Model Parameters/Data
		JJ.priors[[tag]][[m.tag]] <- mod.inputs$Priors
		JJ.form[[tag]][[m.tag]] <- mod.inputs$Formula
		HLA.dat <- mod.inputs$Data
		TAB.hla <- merge( TAB, HLA.dat, by="ID" )
		TAB.h <- TAB.hla[ TAB.hla$ID%in%JJ.samps[[tag]], ]

		# Run Model
		mod.form <- as.formula( JJ.form[[tag]][[m.tag]] )
		mod.prior <- JJ.priors[[tag]][[m.tag]]
		mod.cor <- JJ.cor[[tag]]
		print(paste("Running Model:",m.tag)) ; print(mod.form)
		JJ.time[[tag]][[m.tag]] <- system.time(
			JJ[[tag]][[m.tag]] <- brm( mod.form,
			data=TAB.h, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
            prior=mod.prior, autocor=mod.cor )
		)

		## Save Model Output
		write(paste(date(),"Model",m.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		temp.model <- JJ[[tag]][[m.tag]]
		save( temp.model, file=paste(PathToPlot,"Rdata.Model.",m.tag,".Rdata",sep="") )
		## Print Some Results
		temp.print <- paste(rownames(fixef(JJ[[tag]][[m.tag]])),round(fixef(JJ[[tag]][[m.tag]]),6),sep="=")
		print(paste("Model:",m.tag,"- FIXEF:", temp.print ))
		print(paste("Model:",m.tag,"- TIME:", round( JJ.time[[tag]][[m.tag]][3], 3 ) ))

	}

	## FCT: Plot Certain Fixed Effects from HLA Models
	PLOT_FIXED <- function( mod.fit.table, eff.tag ) {
		mod.fit.table.names <- gsub("_",":",gsub( "Pr4_DRB1_","DRB1*",rownames(mod.fit.table) ))
		YLIM <- extendrange(mod.fit.table[,"Estimate"], f=.2)
		YLIM <- YLIM - c(.1*diff(YLIM), 0)
		png( paste(PathToPlot,"HLA-",m.tag,".2-EffSize.",eff.tag,".png",sep=""), height=1000,width=600+45*nrow(mod.fit.table),pointsize=26)
		par(mar=c( 8,5,5,3 ))
		TEMP <- 1:nrow(mod.fit.table)
		plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=range(TEMP)+c(-.5,.5),ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
		axis( 1, at=TEMP,label=mod.fit.table.names, las=2 )
		axis( 2, at=seq(-10,10,2), las=2 )
		abline( h=-10:10, lty=3,col="grey50",lwd=1 )
		abline( h=0, lty=1,col="black",lwd=1 )
		 # Plot Prior Distributions
		for ( v in 1:nrow(mod.fit.table) ) {
			var <- rownames(mod.fit.table)[v]
			if ( var %in% m.covs.f ) {
				if ( var=="Intercept" ) {
					temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
				}else{
					temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
				}
				if ( length(temp.priors)==2 ) {
					names(temp.priors) <- c("Mean","SD")
					vioplot( rnorm(1e4,temp.priors["Mean"],temp.priors["SD"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
				}
			}
		}
		 # Plot Posterior Distributions
		COLS.temp.list <- COLS.list.2[c(5,4,1)]
		COLS.temp <- adjustcolor(COLS.temp.list[1:2][factor(grepl("DRUG",mod.fit.table.names))],.8)
		COLS.temp[grep(paste(hla.preds,collapse="|"),rownames(mod.fit.table),invert=T)] <- adjustcolor(COLS.temp.list[3],.8)
		for ( v in 1:nrow(mod.fit.table) ) {
			var <- rownames(mod.fit.table)[v]
			var.tag <- var
			if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
			if ( var.tag %in% colnames(d.post) ) {
				vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.temp[v],add=T,drawRect=F )
			}
		}
		arrows( TEMP,mod.fit.table[,"l.95..CI"],TEMP,mod.fit.table[,"u.95..CI"],lwd=3,length=0 )
		arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],lwd=3,length=0 )
		legend("topright",fill=c(adjustcolor("black",.2),adjustcolor(COLS.temp.list,.8)),border=NA,legend=c("Prior","HLA Disease Severity","HLA Drug Response","Clinical Covariate"),ncol=1,title="Effect Type",bg="white")
		temp.text <- abs(.5 - round(d.post.prob[rownames(mod.fit.table),"Pr_l0"],2)) + .5
		text( TEMP, YLIM[1], temp.text, cex=.9 )
		dev.off()
	}

	## Compile & Plot Results
	HLA.mod.fits <- list()
	# Mod.Names <- sapply(strsplit( names(JJ[[tag]]),"_"),"[",2 )
	for ( m in 1:length(Mod.Names) ) {
		mod <- Mod.Names[m]
		m.tag <- paste("hla",mod,sep="_")
		model <- JJ[[tag]][[m.tag]]
		summ <- summary(model,waic=F)
		
		## Collect HLA Info
		mod.inputs <- Modify_HLA_Inputs( mod )
		hla.preds <- setdiff( colnames(mod.inputs$Data), "ID" )
		hla.freqs <- colSums( mod.inputs$Data[,-1] )
		if ( mod=="Type" ) {
			gene.tag <- "DRB1"
			hla.tags <- gsub("_",":",gsub("Pr4_DRB1_","",hla.preds ))
		}else{ hla.tags <- hla.preds }
		
		## Collect General Model Info
		m.obs <- summ$nobs
		m.iter <- summ$iter
		m.warm <- summ$warmup
		m.chain <- summ$chains
		# m.waic <- summ$WAIC
		m.pheno <- as.character(model$formula)[2]
		m.covs <- as.character(model$formula)[3]
		m.form <- paste( m.pheno, "~", m.covs )
		RAND <- length(summ$random)>0

		## Collect Model Outputs
		f.eff <- summ$fixed
		c.eff <- summ$cor_pars
		s.eff <- summ$spec_pars
		m.covs.f <- rownames(f.eff)
		all.eff <- list( f.eff, c.eff, s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		 # Posterior Probabilities
		d.prior <- model$prior
		d.post <- posterior_samples(model)
		d.post.f.which <- paste("b_",m.covs.f,sep="")
		d.post.c.which <- rownames(c.eff)
		d.post.prob <- data.frame( F.eff=d.post.f.which, Pr_l0=round(unlist(lapply( d.post.f.which,function(x)length(which(d.post[,x]>0)) )) / 7200, 2) )
		rownames(d.post.prob) <- m.covs.f
		d.post.prob.rank <- order( abs(d.post.prob$Pr_l0-.5), decreasing=T )
		HLA.mod.fits[[m.tag]] <- mod.fit

		## PLOTS ########

		## Plot Severity vs Response
		# COLS.temp <- COLS.list.2[c(1,4)]
		COLS.temp <- COLS.list.2[c(6,1)]
		png( paste(PathToPlot,"HLA-",m.tag,".1-RespVsSev.png",sep=""), height=1200,width=1200,pointsize=32 )
		# png( paste(PathToPlot,"HLA-",m.tag,".1-RespVsSev.png",sep=""), height=1200,width=2400,pointsize=32 ) ; par(mfrow=c(1,2))
		## Plot Betas
		XVALS <- f.eff[hla.preds,"Estimate"]
		YVALS <- f.eff[paste("DRUG:",hla.preds,sep=""),"Estimate"]
		XLIM <- extendrange( XVALS )
		YLIM <- extendrange( YVALS )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Disease Severity",ylab="Drug Response",main=paste("HLA Allele Effect Sizes:",m.tag) )
		abline( h=seq(-5,5,.1),v=seq(-5,5,.1),lty=3,col="grey50",lwd=1 )
		abline( h=0,v=0 )
		points( XVALS, YVALS, pch=16,col=adjustcolor(COLS.temp[1],.5),cex=2*log10(1+hla.freqs) )
		text( XVALS, YVALS, hla.tags, pos=3,cex=.8 )
		 # Trend Line
		MOD <- lm( YVALS ~ XVALS, weights=hla.freqs )
		PVAL <- summary(MOD)$coefficients["XVALS","Pr(>|t|)"]
		abline( MOD, lty=2,lwd=6,col=COLS.temp[1] )
		text( XLIM[1],YLIM[1],paste("p=",formatC(PVAL,3,format="e"),sep=""), pos=4 )
		# ## Plot Posterior Probs
		# XVALS <- d.post.prob[d.post.prob$F.eff%in%paste("b_",hla.preds,sep=""),"Pr_l0"]
		# YVALS <- d.post.prob[d.post.prob$F.eff%in%paste("b_DRUG:",hla.preds,sep=""),"Pr_l0"]
		# XLIM <- c(0,1) # extendrange( XVALS )
		# YLIM <- c(0,1) # extendrange( YVALS )
		# plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Disease Severity",ylab="Drug Response",main=paste("HLA Allele Posterior Probabilities (B>0):",m.tag),xaxt="n",yaxt="n" )
		# axis( 1, at=seq(0,1,.1) )
		# axis( 2, at=seq(0,1,.1), las=2 )
		# abline( h=seq(0,1,.1),v=seq(0,1,.1),lty=3,col="grey50",lwd=1 )
		# abline( h=.5,v=.5 )
		# points( XVALS, YVALS, pch=16,col=adjustcolor(COLS.temp[2],.5),cex=2*log10(1+hla.freqs) )
		# text( XVALS, YVALS, hla.tags, pos=3,cex=.8 )
		#  # Trend Line
		# MOD <- lm( YVALS ~ XVALS )
		# PVAL <- summary(MOD)$coefficients["XVALS","Pr(>|t|)"]
		# abline( MOD, lty=2,lwd=6,col=COLS.temp[2] )
		# text( XLIM[1],YLIM[1],paste("p=",formatC(PVAL,3,format="e"),sep=""), pos=4 )
		dev.off()

		## Fixed Effect Sizes
		 # Plot top 25 predictors
		temp.n <- 25
		mod.fit.f <- mod.fit[m.covs.f,][d.post.prob.rank[1:temp.n],]
		PLOT_FIXED( mod.fit.f, "top25" )
		 # Plot Clinical Covariates
		mod.fit.c <- mod.fit[grep(paste(hla.preds,collapse="|"),rownames(mod.fit),invert=T),]
		mod.fit.c <- mod.fit.c[intersect(rownames(mod.fit.c),m.covs.f),]
		PLOT_FIXED( mod.fit.c, "clin" )
	}
	save( HLA.mod.fits, file=paste(PathToPlot,"Rdata.HLA_ModFits.Rdata",sep="") )

}
##############################################################
## SNP LMM MODEL #############################################
##############################################################
if ( grepl("SNP",Goal,ignore.case=T) ) {
	write(paste(date(),"- Building SNP Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Set Default Tag
	tag <- "SNP"
	# MOD.priors <- MOD.form <- list()
	JJ.priors[[tag]] <- JJ.form[[tag]] <- list()

	## Set Priors for SNP Variables
	JJ.priors.list$SNP <- set_prior("normal(0,.3)",class="b",coef="SNP")
	JJ.priors.list$DRUG_SNP <- set_prior("normal(0,.3)",class="b",coef="DRUG:SNP")

	#########################################
	## LOAD SNP DATA ########################

	## Load SNP Types
	SNP.raw.l <- read.table( "/projects/janssen/ASSOCIATION/20151201_Previous_CandGWAS/Previous_GWAS_SNPs.tpA.traw", header=T )	
	SNP.cand <- readLines( "/projects/janssen/ASSOCIATION/Previous_GWAS_rsIDs.Uniq.txt" )
	 # Pull out Candidate SNPs
	SNP.raw <- SNP.raw.l[ SNP.raw.l$SNP %in% SNP.cand, ]
	colnames(SNP.raw) <- sapply( strsplit(colnames(SNP.raw),".",fixed=T), "[",1 )
	SNPs <- t( SNP.raw[,-(1:6)] )
	colnames(SNPs) <- SNP.raw$SNP

	####################################
	## Get Models Set Up ###############
	write(paste(date(),"Setting Model Parameters"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Which SNP predictors to use?
	Mods.snp.dat <- list()
	 # DRB1 Type
	Mods.snp.dat$Cand <- data.frame( ID=rownames(SNPs), SNPs ) # SNP_TYP[, c("ID",temp.types) ]

	Mod.Names <- "Cand"

	## SNP9: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac
	# MOD.form <- "DAS ~ (1+DRUG+PLAC|ID)+DRUG*(WK+ACPA+SNP)+PLAC"
	MOD.form <- "DAS ~ (1+DRUG+TRT+WK|ID)+DRUG*(ACPA+SNP)+WK+TRT"
	MOD.priors <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$TRT,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )
	JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )

	####################################
	## Specify Alternative Null Models #

	## FCT: Modify Priors & Formula for various SNP Predictiors
	Modify_SNP_Inputs <- function( mod.name ) {
		## Pull Predictor Info
		snp.data <- Mods.snp.dat[[mod.name]]
		snp.predictors <- colnames(snp.data)[-1]
		## Modify Formula
		snp.form.preds <- paste(snp.predictors,collapse="+")
		OUT.form <- gsub( "SNP",snp.form.preds,MOD.form )
		## Modify Priors
		OUT.priors <- MOD.priors
		for ( snp.pred in snp.predictors ) {
			OUT.priors <- rbind( OUT.priors, JJ.priors.list$SNP, JJ.priors.list$DRUG_SNP )
			OUT.priors[which(OUT.priors$coef=="SNP"),"coef"] <- snp.pred
			OUT.priors[which(OUT.priors$coef=="DRUG:SNP"),"coef"] <- paste("DRUG:",snp.pred,sep="")
		}
		## Return Modified Formula/Priors
		OUT <- list( Formula=OUT.form, Priors=OUT.priors, Data=snp.data )
		return(OUT)
	}

	####################################
	## Run Models ######################

	## Specify Iterations & Samples
	# N.Samps <- 421
	JJ.samps[[tag]] <- SAMPS.short # sample(SAMPS.short,N.Samps)

	## Run Models
	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	if ( !exists("JJ") ) { JJ[[tag]] <- JJ.time[[tag]] <- list() }
	for ( m in 1:length(Mod.Names) ) {
		mod <- Mod.Names[m]
		m.tag <- paste("snp",mod,sep="_")
		mod.inputs <- Modify_SNP_Inputs( mod )
		# Specify Model Parameters/Data
		JJ.priors[[tag]][[m.tag]] <- mod.inputs$Priors
		JJ.form[[tag]][[m.tag]] <- mod.inputs$Formula
		SNP.dat <- mod.inputs$Data
		TAB.snp <- merge( TAB, SNP.dat, by="ID" )
		TAB.h <- TAB.snp[ TAB.snp$ID%in%JJ.samps[[tag]], ]

		# Run Model
		mod.form <- as.formula( JJ.form[[tag]][[m.tag]] )
		mod.prior <- JJ.priors[[tag]][[m.tag]]
		mod.cor <- JJ.cor[[tag]]
		print(paste("Running Model:",m.tag)) ; print(mod.form)
		JJ.time[[tag]][[m.tag]] <- system.time(
			JJ[[tag]][[m.tag]] <- brm( mod.form,
			data=TAB.h, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
            prior=mod.prior ) # , autocor=mod.cor )
		)

		## Save Model Output
		write(paste(date(),"Model",m.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		temp.model <- JJ[[tag]][[m.tag]]
		save( temp.model, file=paste(PathToPlot,"Rdata.Model.",m.tag,".Rdata",sep="") )
		## Print Some Results
		temp.print <- paste(rownames(fixef(JJ[[tag]][[m.tag]])),round(fixef(JJ[[tag]][[m.tag]]),6),sep="=")
		print(paste("Model:",m.tag,"- FIXEF:", temp.print ))
		print(paste("Model:",m.tag,"- TIME:", round( JJ.time[[tag]][[m.tag]][3], 3 ) ))

	}

	####################################
	## Plot Models #####################

	## FCT: Plot Certain Fixed Effects from SNP Models
	PLOT_FIXED <- function( mod.fit.table, eff.tag ) {
		mod.fit.table.names <- gsub("_",":",gsub( "Pr4_DRB1_","DRB1*",rownames(mod.fit.table) ))
		YLIM <- extendrange(mod.fit.table[,"Estimate"], f=.2)
		YLIM <- YLIM - c(.1*diff(YLIM), 0)
		if ( grepl("drug_snp",eff.tag) ) { YLIM <- c(-2,1.5) }
		png( paste(PathToPlot,"SNP-",m.tag,".2-EffSize.",eff.tag,".png",sep=""), height=1000,width=600+45*nrow(mod.fit.table),pointsize=26)
		par(mar=c( 8,5,5,3 ))
		TEMP <- 1:nrow(mod.fit.table)
		plot( 0,0,type="n", xaxt="n",yaxt="n",xlim=range(TEMP)+c(-.5,.5),ylim=YLIM,xlab="",ylab="Effect Size",main="Effect Size Estimates" )
		axis( 1, at=TEMP,label=mod.fit.table.names, las=2 )
		axis( 2, at=seq(-10,10,2), las=2 )
		abline( h=-10:10, lty=3,col="grey50",lwd=1 )
		abline( h=0, lty=1,col="black",lwd=1 )
		 # Plot Prior Distributions
		for ( v in 1:nrow(mod.fit.table) ) {
			var <- rownames(mod.fit.table)[v]
			if ( var %in% m.covs.f ) {
				if ( var=="Intercept" ) {
					temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$class=="temp_Intercept","prior"], fixed=T),fixed=T),"," )[[1]])
				}else{
					temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",d.prior[d.prior$coef==var,"prior"], fixed=T),fixed=T),"," )[[1]])
				}
				if ( length(temp.priors)==2 ) {
					names(temp.priors) <- c("Mean","SD")
					vioplot( rnorm(1e4,temp.priors["Mean"],temp.priors["SD"]), at=TEMP[v],col=adjustcolor("black",alpha=.2),border=NA,add=T,drawRect=F )	
				}
			}
		}
		 # Plot Posterior Distributions
		COLS.temp.list <- COLS.list.2[c(5,4,1)]
		COLS.temp <- adjustcolor(COLS.temp.list[1:2][factor(grepl("DRUG",mod.fit.table.names))],.8)
		COLS.temp[grep(paste(snp.preds,collapse="|"),rownames(mod.fit.table),invert=T)] <- adjustcolor(COLS.temp.list[3],.8)
		if ( grepl("drug_snp",eff.tag) ) { COLS.temp <- rep(adjustcolor(COLS.temp.list[2],.8),length(COLS.temp)) }
		if ( grepl("drug_snp",eff.tag) ) { COLS.temp <- rep(adjustcolor(COLS.list.2[3],.8),length(COLS.temp)) }
		for ( v in 1:nrow(mod.fit.table) ) {
			var <- rownames(mod.fit.table)[v]
			var.tag <- var
			if ( var %in% m.covs.f ) { var.tag <- paste("b",var,sep="_") }
			if ( var.tag %in% colnames(d.post) ) {
				vioplot( d.post[,var.tag], at=TEMP[v],col=COLS.temp[v],add=T,drawRect=F )
			}
		}
		arrows( TEMP,mod.fit.table[,"l.95..CI"],TEMP,mod.fit.table[,"u.95..CI"],lwd=3,length=0 )
		arrows( TEMP-diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],TEMP+diff(TEMP)[1]/2*.6,mod.fit.table[,"Estimate"],lwd=3,length=0 )
		if ( !grepl("drug_snp",eff.tag) ) { legend("topright",fill=c(adjustcolor("black",.2),adjustcolor(COLS.temp.list,.8)),border=NA,legend=c("Prior","SNP Disease Severity","SNP Drug Response","Clinical Covariate"),ncol=2,title="Effect Type",bg="white") }
		temp.text <- abs(.5 - round(d.post.prob[rownames(mod.fit.table),"Pr_l0"],3)) + .5
		text( TEMP, quantile(YLIM,c(0,0.025)), temp.text, cex=.9 )
		dev.off()
	}

	## Compile & Plot Results
	SNP.mod.fits <- list()
	# Mod.Names <- sapply(strsplit( names(JJ[[tag]]),"_"),"[",2 )
	for ( m in 1:length(Mod.Names) ) {
		mod <- Mod.Names[m]
		m.tag <- paste("snp",mod,sep="_")
		model <- JJ[[tag]][[m.tag]]
		summ <- summary(model,waic=F)
		
		## Collect SNP Info
		mod.inputs <- Modify_SNP_Inputs( mod )
		snp.preds <- setdiff( colnames(mod.inputs$Data), "ID" )
		snp.freqs <- colSums( mod.inputs$Data[,-1] )
		if ( mod=="Type" ) {
			gene.tag <- "DRB1"
			snp.tags <- gsub("_",":",gsub("Pr4_DRB1_","",snp.preds ))
		}else{ snp.tags <- snp.preds }
		
		## Collect General Model Info
		m.obs <- summ$nobs
		m.iter <- summ$iter
		m.warm <- summ$warmup
		m.chain <- summ$chains
		# m.waic <- summ$WAIC
		m.pheno <- as.character(model$formula)[2]
		m.covs <- as.character(model$formula)[3]
		m.form <- paste( m.pheno, "~", m.covs )
		RAND <- length(summ$random)>0

		## Collect Model Outputs
		f.eff <- summ$fixed
		c.eff <- summ$cor_pars
		s.eff <- summ$spec_pars
		m.covs.f <- rownames(f.eff)
		all.eff <- list( f.eff, c.eff, s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		 # Posterior Probabilities
		d.prior <- model$prior
		d.post <- posterior_samples(model)
		d.post.f.which <- paste("b_",m.covs.f,sep="")
		d.post.c.which <- rownames(c.eff)
		d.post.prob <- data.frame( F.eff=d.post.f.which, Pr_l0=round(unlist(lapply( d.post.f.which,function(x)length(which(d.post[,x]>0)) )) / 7200, 4) )
		rownames(d.post.prob) <- m.covs.f
		d.post.prob.rank <- order( abs(d.post.prob$Pr_l0-.5), decreasing=T )
		SNP.mod.fits[[m.tag]] <- mod.fit

		## PLOTS ########

		## Plot Severity vs Response
		# COLS.temp <- COLS.list.2[c(1,4)]
		COLS.temp <- COLS.list.2[c(6,1)]
		png( paste(PathToPlot,"SNP-",m.tag,".1-RespVsSev.png",sep=""), height=1200,width=1200,pointsize=32 )
		# png( paste(PathToPlot,"SNP-",m.tag,".1-RespVsSev.png",sep=""), height=1200,width=2400,pointsize=32 ) ; par(mfrow=c(1,2))
		## Plot Betas
		XVALS <- f.eff[snp.preds,"Estimate"]
		YVALS <- f.eff[paste("DRUG:",snp.preds,sep=""),"Estimate"]
		XLIM <- extendrange( XVALS )
		YLIM <- extendrange( YVALS )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Disease Severity",ylab="Drug Response",main=paste("SNP Allele Effect Sizes:",m.tag) )
		abline( h=seq(-5,5,.1),v=seq(-5,5,.1),lty=3,col="grey50",lwd=1 )
		abline( h=0,v=0 )
		points( XVALS, YVALS, pch=16,col=adjustcolor(COLS.temp[1],.5),cex=2*log10(1+snp.freqs) )
		text( XVALS, YVALS, snp.tags, pos=3,cex=.8 )
		 # Trend Line
		MOD <- lm( YVALS ~ XVALS, weights=snp.freqs )
		PVAL <- summary(MOD)$coefficients["XVALS","Pr(>|t|)"]
		abline( MOD, lty=2,lwd=6,col=COLS.temp[1] )
		text( XLIM[1],YLIM[1],paste("p=",formatC(PVAL,3,format="e"),sep=""), pos=4 )
		# ## Plot Posterior Probs
		# XVALS <- d.post.prob[d.post.prob$F.eff%in%paste("b_",snp.preds,sep=""),"Pr_l0"]
		# YVALS <- d.post.prob[d.post.prob$F.eff%in%paste("b_DRUG:",snp.preds,sep=""),"Pr_l0"]
		# XLIM <- c(0,1) # extendrange( XVALS )
		# YLIM <- c(0,1) # extendrange( YVALS )
		# plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="Disease Severity",ylab="Drug Response",main=paste("SNP Allele Posterior Probabilities (B>0):",m.tag),xaxt="n",yaxt="n" )
		# axis( 1, at=seq(0,1,.1) )
		# axis( 2, at=seq(0,1,.1), las=2 )
		# abline( h=seq(0,1,.1),v=seq(0,1,.1),lty=3,col="grey50",lwd=1 )
		# abline( h=.5,v=.5 )
		# points( XVALS, YVALS, pch=16,col=adjustcolor(COLS.temp[2],.5),cex=2*log10(1+snp.freqs) )
		# text( XVALS, YVALS, snp.tags, pos=3,cex=.8 )
		#  # Trend Line
		# MOD <- lm( YVALS ~ XVALS )
		# PVAL <- summary(MOD)$coefficients["XVALS","Pr(>|t|)"]
		# abline( MOD, lty=2,lwd=6,col=COLS.temp[2] )
		# text( XLIM[1],YLIM[1],paste("p=",formatC(PVAL,3,format="e"),sep=""), pos=4 )
		dev.off()

		## Fixed Effect Sizes
		 # Plot top 25 predictors
		temp.n <- 25
		mod.fit.f <- mod.fit[m.covs.f,][d.post.prob.rank[1:temp.n],]
		PLOT_FIXED( mod.fit.f, "top25" )
		 # Plot top 10 DRUG:SNP interactions
		temp.n <- 15
		dr_snp.which <- grep("DRUG:rs",rownames(mod.fit))
		dr_snp.name <- grep("DRUG:rs",rownames(mod.fit),value=T)
		dr_snp.order <- order( abs(.5-d.post.prob[ dr_snp.name, "Pr_l0"]), decreasing=T )
		mod.fit.dr_snp <- mod.fit[dr_snp.name, ]
		mod.fit.dr_snp.n <- mod.fit.dr_snp[ dr_snp.order[1:temp.n], ]
		PLOT_FIXED( mod.fit.dr_snp.n, "15_drug_snp" )
		 # Plot Clinical Covariates
		mod.fit.c <- mod.fit[grep(paste(snp.preds,collapse="|"),rownames(mod.fit),invert=T),]
		mod.fit.c <- mod.fit.c[intersect(rownames(mod.fit.c),m.covs.f),]
		PLOT_FIXED( mod.fit.c, "clin" )
	}
	save( SNP.mod.fits, file=paste(PathToPlot,"Rdata.SNP_ModFits.Rdata",sep="") )

}
##############################################################
## DOWNSAMPLE MODEL ##########################################
##############################################################
if ( grepl("Down",Goal,ignore.case=T) ) {
	write(paste(date(),"- Downsampling Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	####################################
	## Get Models Set Up ###############
	write(paste(date(),"Setting Model Parameters"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Down9: Drug*Week + Placebo + ACPA + Correlation Structure + Random Intercept & Drug & Placebo
	tag <- "down9"
	JJ.form[[tag]] <- "DAS ~ (1+DRUG+TRT|ID)+DRUG*(WK+ACPA)+TRT"
	JJ.priors[[tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$TRT,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )

	####################################
	## Run Models ######################

	## Specify Iterations & Samples
	N.Obs.list <- c(100,250,500,1000,2000)
	N.iters <- 10
	N.iters <- c( 10,8,6,4,2 ) # / 2
	N.iter.tot <- sum( N.iters )
	j <- 1
	which_obs <- list()
	 # Samples
	which_obs <- lapply( 1:length(N.iters), function(x)replicate(N.iters[x],sample(1:nrow(TAB),N.Obs.list[x],replace=F)) )
	# names(which_obs) <- paste("o",rep(N.Obs.list,N.iters),"i",unlist(lapply(N.iters,function(x)1:x)),sep="")
	names(which_obs) <- obs.tag <- paste("o",N.Obs.list,sep="")

	for ( o in 1:length(N.Obs.list) ) {
	# for ( o in 2:length(N.Obs.list) ) {
		obs <- N.Obs.list[o]
		obs.tag <- paste("o",obs,sep="")
		for ( i in 1:N.iters[o] ) {
			## Set Tag
			i.tag <- paste(obs.tag,"i",i,sep="")
			print(paste("Running:",i.tag))
			## Sample from Clinical Data
			# TAB.i <- TAB[ sample(1:nrow(TAB),obs,replace=F), ]
			# TAB.i <- TAB.i[ order(TAB.i$ID,TAB.i$WK), ]
			# which_obs[[i.tag]] <- paste(TAB.i$ID,"_wk",TAB.i$WK,sep="")
			TAB.i <- TAB[ sort(which_obs[[obs.tag]][,i]), ]
			## Run Model
			if ( j==1 ) {
				JJ.time[[tag]][[i.tag]] <- system.time(
					JJ[[tag]][[i.tag]] <- brm( as.formula(JJ.form[[tag]]),
					data=TAB.i, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
					prior=JJ.priors[[tag]] ) # , autocor=mod.cor
				)
				i.tag.1 <- i.tag
			}else{
				JJ.time[[tag]][[i.tag]] <- system.time(
					JJ[[tag]][[i.tag]] <- update( JJ[[tag]][[i.tag.1]],
					newdata=TAB.i, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
					)
				)
			}
			j <- j+1
			write(paste(date(),"Model",i.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
			## Print Some Results
			temp.print <- paste(rownames(fixef(JJ[[tag]][[i.tag]])),round(fixef(JJ[[tag]][[i.tag]]),6),sep="=")
			print(paste("Model:",i.tag,"- FIXEF:", temp.print ))
			print(paste("Model:",i.tag,"- TIME:", round( JJ.time[[tag]][[i.tag]][3], 3 ) ))
			## Print Some Results
			temp.model <- JJ[[tag]][[i.tag]]
			save( temp.model, file=paste(PathToPlot,"Rdata.Model.",i.tag,".Rdata",sep="") )
			save( which_obs, file=paste(PathToPlot,"Rdata.ModObs.",tag,".Rdata",sep="") )
		}
	}

	# tag <- "down9"
	AGGREGATE_MODS <- function( i.tag ) {
		tag.split <- strsplit(i.tag,"i")[[1]]
		temp.obs.ind <- sort( which_obs[[tag.split[1]]][,as.numeric(tag.split[2])] )
		tab.obs <- TAB[ temp.obs.ind, ]
		temp.obs.samp <- as.character( tab.obs$ID )
		temp.obs.dr <- tab.obs$DRUG
		temp.obs.tr <- tab.obs$TRT
		# temp.obs <- paste( TAB$ID[temp.obs.ind], TAB$WK[temp.obs.ind], sep="_")
		# print(head(temp.obs))
		# temp.obs.samp <- sapply(strsplit(temp.obs,"_"),"[",1)
		# tab.obs <- TAB[ match(temp.obs,paste(TAB$ID,"_wk",TAB$WK,sep="")), ]
		# temp.obs.dr <- TAB$DRUG[ match(temp.obs,paste(TAB$ID,"_wk",TAB$WK,sep="")) ]
		# temp.obs.id <- TAB$ID[ match(temp.obs,paste(TAB$ID,"_wk",TAB$WK,sep="")) ]
		# temp.obs.trt <- TAB$TRT[ match(temp.obs,paste(TAB$ID,"_wk",TAB$WK,sep="")) ]
		temp.model <- JJ[[tag]][[i.tag]]
		f.eff <- fixef(temp.model)
		covs.f <- rownames(f.eff)
		r.eff <- coef(temp.model)$ID
		r.eff <- r.eff[ ,apply(r.eff,2,sd)>0 ]
		covs.r <- colnames(r.eff)
		r.diff <- t(t(r.eff)-f.eff[covs.r,])
		r.dist <- apply( r.diff, 1, function(x)sqrt(sum(x^2)) )
		m.obs <- table(temp.obs.samp) # unlist(lapply(m.samps,function(x)length(which(x==temp.obs.samp)) ))
		m.obs.dr <- aggregate( temp.obs.dr, by=list(samp=temp.obs.samp), sum )
		m.obs.tr <- aggregate( temp.obs.tr, by=list(samp=temp.obs.samp), sum )
		m.samps <- names(m.obs)
		mod.tab <- merge( data.frame(Diff=abs(r.diff),Dist=r.dist), data.frame(m.obs,Obs.dr=m.obs.dr[,-1],Obs.tr=m.obs.tr[,-1]), by.x="row.names",by.y="temp.obs.samp" )
		mods.out <- apply( mod.tab[,2:5], 2, function(x) lm(x ~ mod.tab$Freq ) )
		mods.out$Dr <- lm( mod.tab$Diff.DRUG ~ mod.tab$Obs.dr )
		# print( summary(mods.out$Dist) )
		# lapply( mods.out, summary )

		## Compile Outputs
		OUT <- list( Table=mod.tab, Models=mods.out )
	}

	## Aggregate results from each model
	 # Pull & Split model tags
	i.tags <- names(JJ[[tag]])
	i.tag.o <- unlist(lapply( strsplit(i.tags,"i"),function(x)as.numeric(gsub("o","",x[1])) ))
	i.tag.i <- unlist(lapply( strsplit(i.tags,"i"),function(x)as.numeric(x[2]) ))
	 # Get model results
	MODS.agg <- lapply( i.tags, function(x)AGGREGATE_MODS(x) )
	names(MODS.agg) <- i.tags
	 # Pull model coefficients (etc...) for "Dist" models
	MODS.beta.dist <- Reduce( rbind, lapply( MODS.agg, function(x)summary(x$Models$Dist)$coefficients[2,] ))
	MODS.beta.dist <- data.frame( N_obs=i.tag.o, Iter=i.tag.i, MODS.beta.dist )
	rownames(MODS.beta.dist) <- i.tags
	MODS.beta.dist <- MODS.beta.dist[ order(MODS.beta.dist$N_obs,MODS.beta.dist$Iter), ]
	 # Pull model coefficients (etc...) for "Dr" models
	MODS.beta.dr <- Reduce( rbind, lapply( MODS.agg, function(x)summary(x$Models$Dr)$coefficients[2,] ))
	MODS.beta.dr <- data.frame( N_obs=i.tag.o, Iter=i.tag.i, MODS.beta.dr )
	rownames(MODS.beta.dr) <- i.tags
	MODS.beta.dr <- MODS.beta.dr[ order(MODS.beta.dr$N_obs,MODS.beta.dr$Iter), ]
	## FCT: Plot "Dist" & "Dr" Models
	PLOT_MODS <- function( i.tag ) {
		table <- MODS.agg[[i.tag]]$Table
		png( paste(PathToPlot,"Mod-1_DistDr.",i.tag,".png",sep=""),height=1000,width=2000,pointsize=30 )
		par(mfrow=c(1,2))
		# Dist Model
		XLIM <- c( 0, max(table$Freq) )
		YLIM <- c( 0, max(table$Dist) )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="# Observations/Person",ylab="Pythagorean Distance",main="Pythagorean Distance from Pop. Mean vs No. Obs" )
		abline( v=seq(0,20,2),h=seq(0,3,.5),lty=3,col="grey50",lwd=1 )
		points( table$Dist ~ table$Freq, pch=16,col=adjustcolor(COLS.fr,.5) )
		abline(lm( table$Dist ~ table$Freq ), lty=2,col=COLS.fr,lwd=4 )
		text( XLIM[2],YLIM[1], paste("p=",formatC(MODS.beta.dist[i.tag,6],2,format="e"),sep=""),pos=2 )
		text( XLIM[1],YLIM[1], i.tag,pos=4 )
		# Dr Model
		XLIM <- c( 0, max(table$Freq) )
		YLIM <- c( 0, max(table$Diff.DRUG) )
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, xlab="# Observations (on Drug)/Person",ylab="Distance (Drug Coefficient)",main="Distance (Dr) from Pop. vs No. Obs" )
		abline( v=seq(0,20,2),h=seq(0,3,.5),lty=3,col="grey50",lwd=1 )
		points( table$Diff.DRUG ~ table$Freq, pch=16,col=adjustcolor(COLS.fr,.5) )
		abline(lm( table$Diff.DRUG ~ table$Freq ), lty=2,col=COLS.fr,lwd=4 )
		text( XLIM[2],YLIM[1], paste("p=",formatC(MODS.beta.dr[i.tag,6],2,format="e"),sep=""),pos=2 )
		dev.off()
	}
	lapply( i.tags, PLOT_MODS )

}
##############################################################
## MULTIPLE TRIALS MODEL #####################################
##############################################################
if ( grepl("Mult",Goal,ignore.case=T) ) {
	write(paste(date(),"- Building Clinical Trial Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	####################################
	## Get Models Set Up ###############
	write(paste(date(),"Setting Model Parameters"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Load Data w/ Simulated LeadUp Observations
	TAB.l <- read.table( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160726_LeadUp_m91011/Sim.Data.txt",sep="\t",header=T )
	 # Use this data set if specified in "Goal"
	if ( grepl("Lead|Ld",Goal,ignore.case=T) ) {
		TAB.mt <- TAB.l
	}else{
		TAB.mt <- TAB
	}

	## Set Priors

	 # Mult9: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac
	# tag <- "mult9"
	# JJ.priors[[tag]] <- list()
	# JJ.priors[[tag]][[paste(tag,"_1",sep="")]] <- c(JJ.priors.list$Intercept,
	# 	JJ.priors.list$DRUG,
	# 	JJ.priors.list$WK,
	# 	JJ.priors.list$TRT, # TRT [or] PLAC???
	# 	JJ.priors.list$cor,
	# 	JJ.priors.list$ID.sd,
	# 	JJ.priors.list$ACPA,
	# 	JJ.priors.list$DRUG_WK,
	# 	JJ.priors.list$DRUG_ACPA )
	# JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
	# JJ.form[[tag]] <- "DAS ~ (1+DRUG+TRT|ID)+DRUG*(WK+ACPA)+TRT"

	 # Mult10: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac & Week
	tag <- "mult10"
	JJ.form[[tag]] <- "DAS ~ (1+DRUG+TRT+WK|ID)+DRUG*(WK+ACPA)+TRT"
	DRUG.prior.0 <- set_prior("normal(0,2)", class="b",coef="DRUG")
	JJ.priors[[tag]] <- list()
	JJ.priors[[tag]]$Init <- c(JJ.priors.list$Intercept,
		DRUG.prior.0,
		JJ.priors.list$WK,
		JJ.priors.list$TRT,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )

	## Remove DRUG*WK Interactions from Formulae
	# JJ.form$m9 <- "DAS ~ (1+DRUG+TRT|ID)+DRUG*ACPA+WK+TRT"
	JJ.form[[tag]] <- "DAS ~ (1+DRUG+TRT+WK|ID)+DRUG*ACPA+WK+TRT"
	# JJ.form$m11 <- "DAS ~ (1+DRUG+WK|ID)+DRUG*ACPA+WK+TRT"

	## Specify Iterations & Samples
	if ( !exists("JJ.samps") ) { JJ.samps <- list() }
	JJ.samps[[tag]] <- list()
	JJ.samps[[tag]]$P <- as.character( FUL$ID[FUL$GRP=="P"] )
	JJ.samps[[tag]]$G <- as.character( FUL$ID[FUL$GRP=="G"] )
	JJ.samps[[tag]]$PE <- as.character( FUL$ID[FUL$GRP=="PE"] )
	Mod.Names <- c("PE","P","G") # names(JJ.samps[[tag]])
	Mod.N <- length(Mod.Names)

	####################################
	## Run Models ######################

	## Run Models
	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	JJ[[tag]] <- JJ.time[[tag]] <- list()
	for ( m in 1:Mod.N ) {
		## Set Params & Data
		mod <- Mod.Names[m]
		m.tag <- paste(m,mod,sep="_")
		mod.form <- as.formula( JJ.form[[tag]] )
		if ( m==1 ) { JJ.priors[[tag]][[paste(tag,m,sep="_")]] <- JJ.priors[[tag]]$Init }
		mod.prior <- JJ.priors[[tag]][[paste(tag,m,sep="_")]]
		TAB.m <- TAB.mt[ TAB.mt$ID %in% JJ.samps[[tag]][[mod]], ]
		## Run Model
		print(paste("Running Model:",m.tag))
		write(paste(date(),"Running Model:",m.tag), paste(PathToPlot,"Update.txt",sep=""),append=T)
		if ( tag %in% names(JJ.cor) ) {
			mod.cor <- JJ.cor[[tag]]
			JJ.time[[tag]][[m.tag]] <- system.time(
				JJ[[tag]][[m.tag]] <- brm( mod.form,
				data=TAB.m, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior, autocor=mod.cor )
			)
		}else{
			JJ.time[[tag]][[m.tag]] <- system.time(
				JJ[[tag]][[m.tag]] <- brm( mod.form,
				data=TAB.m, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior )
			)
		}
		write(paste(date(),"Model",m.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		## Save Results
		temp.model <- JJ[[tag]][[m.tag]]
		save( temp.model, file=paste(PathToPlot,"Rdata.Model.",m.tag,".Rdata",sep="") )
		## Set Posteriors as New Priors for Next Model
		summ <- summary( temp.model )
		temp.post <- posterior_samples(temp.model)
		new.prior <- temp.model$prior
		 # Fixed Effects
		temp.cov.f <- rownames(summ$fixed) # grep("^b_",colnames(temp.post),value=T)
		temp.cov.f.tag <- paste( "b_",temp.cov.f,sep="" )
		 # Random Effects
		temp.cov.r <- rownames(summ$random$ID)
		temp.cov.r.tag <- gsub(")","", gsub(",","_", gsub("(","_ID_",temp.cov.r,fixed=T) ))
		 # Sigma Pheno
		temp.cov.s <- rownames(summ$spec_pars)
		temp.cov.s.tag <- gsub(")","", gsub("(","_",temp.cov.s,fixed=T) )
		 # Go through all variables & set new prior based on posterior
		temp.cov <- c( temp.cov.f, temp.cov.r, temp.cov.s )
		temp.cov.tag <- c( temp.cov.f.tag, temp.cov.r.tag, temp.cov.s.tag )
		for ( c in 1:length(temp.cov) ) {
		# for ( c in 1:6 ) {
			cov <- temp.cov[c]
			cov.tag <- temp.cov.tag[c]
			cov.tag.1 <- strsplit(cov.tag,"_")[[1]][1]
			if ( cov=="Intercept" ) {
				which_row <- grep("Intercept",new.prior[,"class"])
				new.mn <- round( mean( temp.post[,cov.tag] ), 4)
				new.sd <- round( sd( temp.post[,cov.tag] ), 4)
				new.text <- paste("normal(",new.mn,",",new.sd,")",sep="")
				new.prior[ which_row, "prior" ] <- new.text
				next
			}
			if ( cov.tag.1=="b" ) {
				which_row <- which( new.prior[,"class"]==cov.tag.1 & new.prior[,"coef"]==cov )
				new.mn <- round( mean( temp.post[,cov.tag] ), 4)
				new.sd <- round( sd( temp.post[,cov.tag] ), 4)
				new.sd <- 1.5 * new.sd
				new.text <- paste("normal(",new.mn,",",new.sd,")",sep="")
				new.prior[ which_row, "prior" ] <- new.text
				next
			}
			if ( cov.tag.1=="sd" ) {
				cov.tag.2 <- tail( strsplit(cov.tag,"_")[[1]], 1 )
				which_row <- which( new.prior[,"class"]==cov.tag.1 & new.prior[,"coef"]==cov.tag.2 )
				new.mn <- round( mean( temp.post[,cov.tag] ), 4)
				new.sd <- round( sd( temp.post[,cov.tag] ), 4)
				new.sd <- 1.5 * new.sd
				new.text <- paste("cauchy(",new.mn,",",new.sd,")",sep="")
				new.prior[ which_row, "prior" ] <- new.text
				next
			}
			if ( cov.tag.1=="sigma" ) {
				cov.tag.2 <- tail( strsplit(cov.tag,"_")[[1]], 1 )
				which_row <- which( new.prior[,"class"]==cov.tag.1 & new.prior[,"coef"]==cov.tag.2 )
				new.mn <- round( mean( temp.post[,cov.tag] ), 4)
				new.sd <- round( sd( temp.post[,cov.tag] ), 4)
				new.sd <- 1.5 * new.sd
				new.text <- paste("normal(",new.mn,",",new.sd,")",sep="")
				new.prior[ which_row, "prior" ] <- new.text
				next
			}

		}
		 # Compile as New Priors
		JJ.priors[[tag]][[paste(tag,m+1,sep="_")]] <- new.prior
		## Print Some Results
		temp.print <- paste(rownames(fixef(JJ[[tag]][[m.tag]])),round(fixef(JJ[[tag]][[m.tag]]),6),sep="=")
		print(paste("Model:",m.tag,"- FIXEF:", temp.print ))
		print(paste("Model:",m.tag,"- TIME:", round( JJ.time[[tag]][[m.tag]][3], 3 ) ))
	}

}
##############################################################
## INDIVIDUAL MODELS #########################################
##############################################################
if ( grepl("Ind",Goal,ignore.case=T) ) {
	write(paste(date(),"- Building Clinical Trial Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	####################################
	## Get Models Set Up ###############
	write(paste(date(),"Setting Model Parameters"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	## Set Priors for each Variable

	## IND4: Drug*Week + Placebo
	if ( grepl("ind4",Mod.Names) ) {
		tag <- "ind4"
		JJ.priors[[tag]] <- c(JJ.priors.list$Intercept,
			JJ.priors.list$DRUG,
			JJ.priors.list$WK,
			JJ.priors.list$TRT,
			JJ.priors.list$DRUG_WK )
		JJ.form[[tag]] <- "DAS ~ DRUG*WK+TRT"	
	}
	
	## IND3: Drug + Week + Placebo
	if ( grepl("ind3",Mod.Names) ) {
		tag <- "ind3"
		JJ.priors[[tag]] <- c(JJ.priors.list$Intercept,
			JJ.priors.list$DRUG,
			JJ.priors.list$WK,
			JJ.priors.list$TRT )
		JJ.form[[tag]] <- "DAS ~ DRUG+WK+TRT"
	}

	####################################
	## Run Models ######################

	## Pull out Sample Names
	Samps <- as.character(unique( TAB$ID )) # as.character( FUL$ID )
	N.Samps <- length(Samps)

	## Run Models
	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	JJ[[tag]] <- JJ.time[[tag]] <- list()
	s.tag.1g <- s.tag.1p <- ""
	for ( s in 1:N.Samps ) {
	# for ( s in 1:30 ) {
		## Set Params & Data
		samp <- Samps[s]
		s.tag <- paste(tag,samp,sep="_")
		TAB.s <- TAB[ TAB$ID==samp, ]
		mod.form <- as.formula( JJ.form[[tag]] )
		mod.prior <- JJ.priors[[tag]]
		## Run Model
		print(paste("Running Model",s,":",s.tag))
		## Run Model
		if ( s==1 ) {
			JJ.time[[tag]][[s.tag]] <- system.time(
				JJ[[tag]][[s.tag]] <- brm( mod.form,
				data=TAB.s, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior )
			)
			s.tag.1 <- s.tag
		}else{
			JJ.time[[tag]][[s.tag]] <- system.time(
				JJ[[tag]][[s.tag]] <- update( JJ[[tag]][[s.tag.1]],
				newdata=TAB.s, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
				)
			)
		}
		write(paste(date(),"Model",s.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		print(paste( "Done with:",s,"-",s.tag ))
	}
	
	## Save Results (***TOO BIG TO SAVE***)
	# temp.model <- JJ[[tag]]
	# save( temp.model, file=paste(PathToPlot,"Rdata.Model.",tag,".Rdata",sep="") )
	

	##########################################
	## Pull/Compile Model Results ############
	model <- JJ[[tag]]
	COMPILE <- list()

	## Calculate WAIC for each Patient's model
	WAIC <- Reduce( rbind, lapply( model, WAIC ) )
	rownames(WAIC) <- names(model)
	COMPILE$WAIC <- WAIC

	## Calculate Predicted Values for each Patient's model
	PRED <- lapply( model, predict )
	COMPILE$PRED <- PRED

	## Pull out Estimates for Fixed Effects for each Patient
	covs.f <- rownames( fixef(model[[1]]) )
	model.fixef <- lapply( model, function(x)summary(x)$fixed )
	 # Compile into a table for each variable (including intervals)
	EFF <- lapply( covs.f, function(x)Reduce(rbind,lapply(model.fixef,function(y)y[x,] )) )
	names(EFF) <- covs.f
	for ( cov in covs.f ) { rownames(EFF[[cov]]) <- names(model.fixef) }
	COMPILE$EFF <- EFF
	 # Compile Betas for all variables/patients into 1 table
	EST <- Reduce( cbind, lapply(EFF,function(x)x[,"Estimate"]) )
	colnames(EST) <- covs.f
	COMPILE$EST <- EST

	## Calculate Prob of effect size beyond various thresholds
	which_vars <- paste("b_",c("Intercept","DRUG","WK","TRT"), sep="" )
	POST.1 <- lapply( model,function(x)posterior_samples(x)[,which_vars] )
	Thresholds <- seq(-2,0,.1) # list( DRUG=seq(-2,0,.2), TRT=seq(-1,0,.2) )
	POST.2 <- lapply(POST.1,
		function(x)Reduce( rbind, lapply(Thresholds,
			function(t)apply(x[,c("b_DRUG","b_TRT")],2,
				function(y)length(which(y<t)))/7200
			)
		)
	)
	POST <- lapply(colnames(POST.2[[1]]),function(v)Reduce( rbind,lapply(POST.2,function(x)t(x[,v])) ) )
	names(POST) <- colnames(POST.2[[1]])
	for ( i in 1:length(POST) ) {
		colnames(POST[[i]]) <- paste("lt",Thresholds,sep="_")
		rownames(POST[[i]]) <- names(POST.1)
	}
	COMPILE$POST <- POST
	COMPILE$POST_RAW <- POST.1
	save( COMPILE, file=paste(PathToPlot,"Rdata.Model.Summarize.",tag,".Rdata",sep="") )
	# save( POST.1, file=paste(PathToPlot,"Rdata.Model.Summarize.PostRaw.",tag,".Rdata",sep="") )

	##########################################
	## Plot Model Summaries ##################

	EST.2 <- COMPILE$EST
	rownames(EST.2) <- gsub(paste(tag,"_",sep=""),"",rownames(EST.2))
	EST.2 <- merge( EST.2, FUL[,c("ID","DAS_BL_MN","DEL_MNe_MN","GRP")], by.x="row.names",by.y="ID" )

	## 1 - Distributions of Intercept, Drug, & Treatment Effects
	Vars <- c("Intercept","DRUG","TRT")
	 # Set Plotting Parameters
	COLS.temp <- COLS.list.2[c(1,7,4)]
	PARAMS <- list( Intercept=list(breaks=seq(3,8,.5),color=COLS.temp[1],xlab="Baseline DAS",main="Individual Baseline DAS"),
		DRUG=list(breaks=seq(-3,1,.25),color=COLS.temp[2],xlab="Drug Effect",main="Individual Drug Effect"),
		TRT=list(breaks=seq(-2,2,.25),color=COLS.temp[3],xlab="Treatment Effect",main="Individual Treatment Effect")
	)
	## FCT: Plot Distributions of Individual Estimates
	Plot_Distrib <- function( var, params ) {
		## Plotting Params
		BRKS <- params$breaks
		Color <- params$color
		MAIN <- params$main
		XLAB <- params$xlab
		hist( EST.2[,var], breaks=BRKS,col=Color, xlab=XLAB,main=MAIN,ylab="# Patients" )
		# hist( EST.2[,var], breaks=BRKS,col="grey80", xlab=XLAB,main=MAIN,ylab="# Patients" )
		# hist( EST.2[EST.2$GRP=="G",var],breaks=BRKS,col=adjustcolor(COLS.temp["G"],.5),add=T )
		# hist( EST.2[EST.2$GRP=="P",var],breaks=BRKS,col=adjustcolor(COLS.temp["P"],.5),add=T )
		# hist( EST.2[EST.2$GRP=="PE",var],breaks=BRKS,col=adjustcolor(COLS.temp["PE"],.5),add=T )
	}
	## Save Plot of Distributions
	png( paste(PathToPlot,"Ind_1-Distributions.png",sep=""),height=800,width=1600,pointsize=30 )
	par(mfrow=c(1,3))
	SCRAP <- lapply( Vars, function(x)Plot_Distrib( x, PARAMS[[x]] ) )
	dev.off()

	## 2 - Scatter plots of Intercept, Drug, & Treatment Effects
	 # Set Plotting Parameters
	COLS.temp <- COLS.list.2[c(2,3,6)] ; names(COLS.temp) <- c("G","P","PE")
	COLS.scat <- COLS.temp[ EST.2$GRP ]
	PARAMS <- list( Int_v_Drug=list(main="Individual Estimates",xlab="Baseline DAS",ylab="Drug Effect",xlim=c(3,8),ylim=c(-3,1),colors=COLS.scat ),
		Int_v_Trt=list(main="Individual Estimates",xlab="Baseline DAS",ylab="Treatment Effect",xlim=c(3,8),ylim=c(-2,2),colors=COLS.scat ),
		Drug_v_Trt=list(main="Individual Estimates",ylab="Treatment Effect",xlab="Drug Effect",xlim=c(-3,1),ylim=c(-2,2),colors=COLS.scat )
	)
	## FCT: Plot Distributions of Individual Estimates
	Plot_Scatter <- function( var_1, var_2, params ) {
		## Plotting Params
		MAIN <- params$main
		XLAB <- params$xlab
		YLAB <- params$ylab
		XLIM <- params$xlim
		YLIM <- params$ylim
		Colors <- params$colors
		# Plot/Lines
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab=XLAB,ylab=YLAB,main=MAIN )
		abline( h=seq(-10,10,1),v=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
		# Data/Fit
		points( EST.2[,var_1], EST.2[,var_2], col=adjustcolor(Colors,.5),pch=16 )
		MOD <- lm( EST.2[,var_2] ~ EST.2[,var_1] )
		abline( MOD, lty=2,lwd=4,col="grey20" )
		P <- summary(MOD)$coefficients[2,4]
		text( XLIM[1],YLIM[1], paste("p=",formatC(P,2,format="e")),pos=4 )
	}
	## Save Plot of Distributions
	png( paste(PathToPlot,"Ind_2-Scatter.png",sep=""),height=600,width=1800,pointsize=30 )
	par(mfrow=c(1,3))
	Plot_Scatter( "Intercept","DRUG", PARAMS$Int_v_Drug )
	Plot_Scatter( "Intercept","TRT", PARAMS$Int_v_Trt )
	Plot_Scatter( "DRUG","TRT", PARAMS$Drug_v_Trt )
	dev.off()

	## 3 - Violin Plots of Estimates by Arms
	COLS.temp <- COLS.list.2[c(2,3,6)] ; names(COLS.temp) <- c("G","P","PE")
	PARAMS <- list( Intercept=list(ylim=c(3,8),ylab="Baseline DAS"),
		DRUG=list(ylim=c(-3,1),ylab="Drug Effect"),
		TRT=list(ylim=c(-2,2),ylab="Treatment Effect")
	)
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


	# plot( EST.2$DEL_MNe_MN ~ POST$b_DRUG[,"lt_-1"],pch=16,col=adjustcolor(COLS.temp[factor(EST.2$GRP)],.6) )

}
##############################################################
## LEAD-UP BASELINE MEASUREMENT MODEL ########################
##############################################################
if ( grepl("Lead",Goal,ignore.case=T) ) {
	write(paste(date(),"- Building Lead-Up Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Set Default Tag
	tag <- "Lead"
	MOD.priors <- MOD.form <- list()
	JJ.priors[[tag]] <- JJ.form[[tag]] <- list()

	 # M9: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac
	m.tag <- "m9"
	MOD.form[[m.tag]] <- "DAS ~ (1+DRUG+TRT|ID)+DRUG*(WK+ACPA)+TRT"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$TRT,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )
	 # M10: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Plac & Week
	m.tag <- "m10"
	MOD.form[[m.tag]] <- "DAS ~ (1+DRUG+TRT+WK|ID)+DRUG*(WK+ACPA)+TRT"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$TRT,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )
	 # M11: Drug*(Week+ACPA) + Placebo + Random Intercept & Drug & Week
	m.tag <- "m11"
	MOD.form[[m.tag]] <- "DAS ~ (1+DRUG+WK|ID)+DRUG*(WK+ACPA)+TRT"
	MOD.priors[[m.tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$TRT,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$ID.sd,
		JJ.priors.list$cor )


	## Models w/ Remove DRUG*WK Interactions
	MOD.form$m9 <- "DAS ~ (1+DRUG+TRT|ID)+DRUG*ACPA+WK+TRT"
	MOD.form$m10 <- "DAS ~ (1+DRUG+TRT+WK|ID)+DRUG*ACPA+WK+TRT"
	MOD.form$m11 <- "DAS ~ (1+DRUG+WK|ID)+DRUG*ACPA+WK+TRT"

	####################################
	## Specify Alternative Null Models #

	## Models w/ Correlation Stucture
	if ( grepl("cor",Goal,ignore.case=T) ) {
		# Default Correlation Structure b/n Weeks
		tag <- paste(tag,"c",sep="")
		JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
		print( JJ.cor[[tag]] )
	}

	####################################
	## Run Models ######################

	## Pull/Specify Model Names/Info/Etc...
	 # Mod Names
	if ( is.na(Mod.Names.1) ) { Mod.Names <- grep("^m",names(MOD.form),value=T) }
	Mod.N <- length(Mod.Names)
	 # Priors & Formulae
	JJ.priors[[tag]] <- MOD.priors
	JJ.form[[tag]] <- MOD.form

	## Specify Iterations & Samples
	# N.Samps <- 100
	JJ.samps[[tag]] <- sample(SAMPS.short,N.Samps)
	JJ.samps[[tag]] <- SAMPS.short
	# save( JJ.samps, file=paste(PathToPlot,"Rdata.Samps.Rdata",sep="") )

	## Simulate Lead-Up Measurements for Patients
	 # Based on the lone baseline measurement
	TAB.temp <- TAB[ which( TAB$ID%in%JJ.samps[[tag]] & TAB$WK==0), ]
	TAB.temp$WK <- -2
	TAB.temp$DAS <- TAB.temp$DAS + rnorm(N.Samps,0,1)
	TAB.l <- rbind( TAB[TAB$ID%in%JJ.samps[[tag]],], TAB.temp )
	TAB.l <- TAB.l[ order(TAB.l$ID,as.numeric(TAB.l$WK)), ]
	TAB.l <- read.table( "/projects/janssen/Mac_Data/Janssen/Plots_Mac/20160726_LeadUp_m91011/Sim.Data.txt",sep="\t",header=T )
	write.table( TAB.l, paste(PathToPlot,"Sim.Data.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	if ( tag %in% names(JJ.cor) ) { write(paste("Cor:",as.character(JJ.cor[[tag]])[[1]]), paste(PathToPlot,"Update.txt",sep=""),append=T) }
	## Run Models
	print(paste("Starting:",tag))
	print(paste("...including models:",paste(Mod.Names,collapse=",")))
	JJ[[tag]] <- list()
	for ( m in 1:Mod.N ) {
		m.tag <- Mod.Names[m]
		mod.form <- as.formula( JJ.form[[tag]][[m.tag]] )
		mod.prior <- JJ.priors[[tag]][[m.tag]]
		JJ[[tag]][[m.tag]] <- JJ.time[[tag]][[m.tag]] <- list()
		print(paste("Running Model:",m.tag))
		write(paste(date(),"- Running Model:",m.tag,"######"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		print(mod.form)
		if ( tag %in% names(JJ.cor) ) {
			mod.cor <- JJ.cor[[tag]] ; print(mod.cor)
			JJ.time[[tag]][[m.tag]][[m.tag]] <- system.time(
				JJ[[tag]][[m.tag]][[m.tag]] <- brm( mod.form,
				data=TAB.l, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior, autocor=mod.cor )
			)
		}else{
			JJ.time[[tag]][[m.tag]][[m.tag]] <- system.time(
				JJ[[tag]][[m.tag]][[m.tag]] <- brm( mod.form,
				data=TAB.l, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior )
			)
			m.tag.1 <- m.tag
		}
		write(paste(date(),"Model",m.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		## Print Some Results
		temp.print <- paste(rownames(fixef(JJ[[tag]][[m.tag]][[m.tag]])),round(fixef(JJ[[tag]][[m.tag]][[m.tag]]),6),sep="=")
		print(paste("Model:",m.tag,"- FIXEF:", temp.print ))
		print(paste("Model:",m.tag,"- TIME:", round( JJ.time[[tag]][[m.tag]][[m.tag]][3], 3 ) ))
		## Print Some Results
		temp.model <- JJ[[tag]][[m.tag]][[m.tag]]
		save( temp.model, file=paste(PathToPlot,"Rdata.Model.",m.tag,".Rdata",sep="") )
	}
	write(paste(date(),"Done Building Models!!"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Write Time Results
	save( JJ.time, file=paste(PathToPlot,"Rdata.Time.Rdata",sep="") )

}
##############################################################
## PLOT MODEL SUMMARIES ######################################
##############################################################
write(paste(date(),"- Plotting Model Summaries ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

# # Plot/Save Models
# PLOT_MOD( JJ$m1$m1i1, "m1Test" )
# PLOT_MOD( JJ$m5$m5i1, "m5Test" )
# SCRAP <- lapply( grep("^m",names(JJ),value=T),function(x)PLOT_MOD(JJ[[x]][[1]],names(JJ[[x]])[1]) )
# SCRAP <- lapply( grep("^m",names(JJ),value=T),function(x)lapply(names(JJ[[x]]),function(y)PLOT_MOD(JJ[[x]][[y]],y) ))

# for ( x in grep("^m",names(JJ),value=T) ) {
# 	PLOT_MOD(JJ[[x]][[1]],names(JJ[[x]])[1])
# }

# for ( x in names(JJ$hla10) ) {
# 	PLOT_MOD(JJ$hla10[[x]],x)
# }

#############################################################
#############################################################
# END OF DOC ################################################
#############################################################
#############################################################

#############################################################
## VALIDATE USE OF PRIORS ###################################
#############################################################
# # Validate Robustness w/ use of Priors
#  Use simple/small models (n=50?)
#    Test many different means and sds for priors
#  Use complex/small model (n=50)
#    Test a few different means for priors
#  Use simple/large model (n=250+)
#    Test a few different means for priors
#  How does Sample Size influence priors?
#    Larger sample size means smaller influence of priors...
#    Show this?

####################################
## Sample Size Model ###############
## Vary Sample Size ################
if ( grepl("Prior",Goal,ignore.case=T) & grepl("Size",Mod.Names.1,ignore.case=T) ) {
	write(paste(date(),"- Validate: Sample Size ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Simp1: Simple Model
	tag <- "Size1"
	JJ.time[[tag]] <- JJ[[tag]] <- list()
	 # Specify Simple Model Formula
	JJ.form[[tag]] <- "DAS ~ DRUG"
	 # Specify Priors (Mean & SD)
	Pr.int.mn <- 5
	Pr.int.sd <- 2
	Pr.int <- paste("normal(",Pr.int.mn,",",Pr.int.sd,")",sep="")
	Pr.dr.mn <- 0 # seq(-3,1,.5)
	Pr.dr.sd <- 1 # seq(.5,4,.5)
	Pr.dr <- paste("normal(",Pr.dr.mn,",",Pr.dr.sd,")",sep="")
	Pr.prior <- c( set_prior(Pr.int,class="b",coef="Intercept"),
		set_prior(Pr.dr,class="b",coef="DRUG") )
	JJ.priors[[tag]] <- Pr.prior
	 # Specify Sample Sizes to use & Number of Iterations
	Pr.sizes <- c( 5,10,15,20,30,40,50, 75,100, 150,200,300,length(SAMPS.short) )
	Pr.n.iters <- c( rep(10,7), 5,5, 3,3,3,3 )
	# Pr.n.iters <- rep( 2, length(Pr.sizes) )
	Pr.n.iter.tot <- sum( Pr.n.iters )
	JJ.samps[[tag]] <- lapply( 1:length(Pr.sizes), function(x)replicate(Pr.sizes[x],sample(SAMPS.short,Pr.n.iters[x])) )
	names(JJ.samps[[tag]]) <- paste("s",Pr.sizes,sep="")
	j <- 1
	for ( s in 1:length(Pr.sizes) ) {
		sz <- Pr.sizes[s]
		sz.tag <- paste("s",sz,sep="")
		for ( i in 1:Pr.n.iters[s] ) {
			## Set Tag
			pr.tag <- paste("s",sz,"i",i,sep="")
			print(paste("Running:",pr.tag))
			## Sample from Clinical Data
			# print(JJ.samps[[tag]][[sz.tag]][i,])
			TAB.pr <- TAB[ TAB$ID%in%JJ.samps[[tag]][[sz.tag]][i,], ]
			# print(dim(TAB.pr))
			## Run Model
			if ( j==1 ) {
				JJ.time[[tag]][[pr.tag]] <- system.time(
					JJ[[tag]][[pr.tag]] <- brm( as.formula(JJ.form[[tag]]),
					data=TAB.pr, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
					prior=Pr.prior ) # , autocor=mod.cor
				)
				pr.tag.1 <- pr.tag
			}else{
				JJ.time[[tag]][[pr.tag]] <- system.time(
					JJ[[tag]][[pr.tag]] <- update( JJ[[tag]][[pr.tag.1]],
					newdata=TAB.pr, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
					) # prior=Pr.prior, autocor=mod.cor
				)
			}
			## Print Some Results
			temp.print <- paste(rownames(fixef(JJ[[tag]][[pr.tag]])),round(fixef(JJ[[tag]][[pr.tag]]),6),sep="=")
			print(paste("Model:",pr.tag,"- FIXEF:", temp.print ))
			print(paste("Model:",pr.tag,"- TIME:", round( JJ.time[[tag]][[pr.tag]][3], 3 ) ))
			# Next J
			j <- j + 1
		}
	}

	## Pull out Fixed Effects & Distributions
	Sz.fixef <- data.frame( MOD=names(JJ.time[[tag]]),gsub("[a-z]","",t(sapply(strsplit(names(JJ.time[[tag]]),"i"),"[",1:2))), t(Reduce( cbind, lapply(JJ[[tag]],fixef) )), stringsAsFactors=F )
	colnames(Sz.fixef)[2:3] <- c("Size","Iter")
	Sz.time <- data.frame( MOD=names(JJ.time[[tag]]), Reduce( rbind, JJ.time[[tag]] ), stringsAsFactors=F )
	Sz.mg <- merge( Sz.fixef,Sz.time, by="MOD" )
	 # Write Tables for Later Use
	write.table( Sz.mg, paste(PathToPlot,"TAB.",tag,"-Fixef.Time.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	Sz.mg <- read.table( paste(PathToPlot,"TAB.",tag,"-Fixef.Time.txt",sep=""),sep="\t",header=T )

	## Summarize Fixed Effects by Sample Size
	Sz.fef.mn <- aggregate( Sz.mg[,c("Intercept","DRUG")], by=list(Size=Sz.mg$Size), mean )
	Sz.fef.sd <- aggregate( Sz.mg[,c("Intercept","DRUG")], by=list(Size=Sz.mg$Size), sd )
	Sz.fef.n <- table(Sz.mg$Size)
	Sz.comp <- data.frame( N=Sz.fef.n, MN=Sz.fef.mn[,-1], SD=Sz.fef.sd[,-1], stringsAsFactors=F )
	colnames(Sz.comp)[1:2] <- c("Size","N.Iter")
	Sz.comp$Size <- as.numeric(as.character(Sz.comp$Size))
	Sz.comp <- Sz.comp[order(as.numeric(Sz.comp$Size)), ]

	###################
	## Plot ###########
	COLS <- COLS.list.2

	## Plot Posterior Estimates
	YLIM <- extendrange( Sz.comp[,c("MN.Intercept","MN.DRUG")], f=.2)
	XLIM <- range( Sz.comp$Size )
	png( paste(PathToPlot,tag,".1-BetaVsSize.png",sep=""),height=800,width=1200,pointsize=24)
	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Posterior Effect Size",xlab="Sample Size",main="Posterior Estimate vs Sample Size" )
	abline(h=-10:10,lty=3,col="grey50",lwd=1 )
	abline(h=0)
	points( Intercept ~ Size, data=Sz.mg, pch=16,col=adjustcolor(COLS[1],alpha=.5) )
	points( DRUG ~ Size, data=Sz.mg, pch=16,col=adjustcolor(COLS[2],alpha=.5) )
	points( MN.Intercept ~ Size, data=Sz.comp, type="o",pch=1,col=COLS[1],cex=1.5,lwd=3 )
	points( MN.DRUG ~ Size, data=Sz.comp, type="o",pch=1,col=COLS[2],cex=1.5,lwd=3 )
	legend("topright",pch=16,col=COLS[1:2],legend=c("Intercept","DRUG"),ncol=2 )
	dev.off()

	## Violin Plot: Prior vs Posterior
	COLS.sz <- rep(COLS.list.2[6],nrow(Sz.comp)) # colorRampPalette(COLS.list.ord)(nrow(Sz.comp))
	temp.priors <- as.numeric(strsplit( gsub(")","",gsub("normal(","",JJ.priors[[tag]][JJ.priors[[tag]]$coef=="DRUG","prior"], fixed=T),fixed=T),"," )[[1]])
	names(temp.priors) <- c("Mean","SD")
	# XLIM <- c(1,nrow(Sz.mg))
	# YLIM <- c(-6,4)
	# png( paste(PathToPlot,tag,".2-PriorVsPost.png",sep=""),height=800,width=1600,pointsize=24)
	# plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="",main="Posterior Estimate vs Prior",xaxt="n" )
	# axis( 1, at=1:nrow(Sz.mg), label=Sz.mg$MOD, las=2)
	# abline(h=-10:10,lty=3,col="grey50",lwd=1 )
	# abline(h=0)
	# SCRAP <- lapply( 1:nrow(Sz.mg), function(x)vioplot(rnorm(10000,temp.priors["Mean"],temp.priors["SD"]),at=x, col=adjustcolor(COLS.sz[factor(Sz.mg$Size)[x]],alpha=.7),border=NA,add=T ))
	# SCRAP <- lapply( 1:nrow(Sz.mg), function(x)vioplot(posterior_samples(JJ[[tag]][[x]])$b_DRUG,at=x, col=adjustcolor("black",alpha=.7),border=NA,add=T ))
	# legend( "topleft",fill=adjustcolor(c("black",COLS.sz[1:length(unique(Sz.mg$Size))]),alpha=.7),border=NA,legend=c("Posterior",paste("Prior: Size=",unique(Sz.mg$Size),sep="")),ncol=length(unique(Sz.mg$Size))+1 )
	# dev.off()
	## Violin Plot: Prior vs Posterior (Compiled)
	XLIM <- c(1,nrow(Sz.comp))
	YLIM <- c(-6,4)
	png( paste(PathToPlot,tag,".2b-PriorVsPost.png",sep=""),height=800,width=1600,pointsize=24)
	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="# Patients",main="Posterior Population Estimate vs Sample Size",xaxt="n" )
	axis( 1, at=1:nrow(Sz.comp), label=paste("s",Sz.comp$Size,sep=""), las=2)
	abline(h=-10:10,lty=3,col="grey50",lwd=1 )
	abline(h=0)
	SCRAP <- lapply( 1:nrow(Sz.comp), function(x)vioplot(rnorm(10000,temp.priors["Mean"],temp.priors["SD"]),at=x, col=adjustcolor(COLS.sz[factor(Sz.comp$Size)[x]],alpha=.7),border=NA,add=T ))
	SCRAP <- lapply( 1:nrow(Sz.comp), function(x)vioplot(unlist(lapply( JJ[[tag]][grep(paste("s",Sz.comp[x,"Size"],sep=""),names(JJ[[tag]]))], function(y)posterior_samples(y)$b_DRUG)),at=x, col=adjustcolor("black",alpha=.7),border=NA,add=T ))
	legend( "topleft",fill=adjustcolor(c("black",COLS.sz),alpha=.7),border=NA,legend=c("Posterior","Prior"),ncol=2 )
	text( 1:nrow(Sz.comp),YLIM[1]+.6,label=paste("N.iter=",Sz.comp$N.Iter,sep=""),srt=90)
	dev.off()

	## Plot Time vs Size
	YLIM <- c(0,max(Sz.mg$elapsed[Sz.mg$Iter!=1]) )
	png( paste(PathToPlot,tag,".3-TimeVsSize.png",sep=""),height=800,width=800,pointsize=24)
	plot( elapsed ~ Size, data=Sz.mg,subset=Iter!=1, pch=16,ylim=YLIM,col=adjustcolor(COLS[1],alpha=.7),xlab="# Patients",ylab="Time (s)",main="Computational Burden vs Sample Size" )
	abline( h=seq(0,YLIM[2]+10,10),v=seq(0,500,20),lty=3,col="grey50",lwd=1 )
	dev.off()
}
####################################
## Small/Simple Model ##############
## Vary Mean & SD ##################
if ( grepl("Prior",Goal,ignore.case=T) & grepl("Simp1",Mod.Names.1,ignore.case=T) ) {
	write(paste(date(),"- Validate: Mean/SD on Simple Mods ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Simp1: Simple Model
	tag <- "Simp1"
	JJ.time[[tag]] <- JJ[[tag]] <- list()
	 # Specify Simple Model Formula
	JJ.form[[tag]] <- "DAS ~ DRUG"
	 # Specify Levels of Mean & SD
	Pr.int.mn <- 5
	Pr.int.sd <- 2
	Pr.int <- paste("normal(",Pr.int.mn,",",Pr.int.sd,")",sep="")
	Pr.dr.mns <- -3:1 # seq(-3,1,.5)
	Pr.dr.sds <- c(.03,.1,.3,1,2) # seq(.5,4,.5)
	 # Run Model a bunch of times
	Pr.n.iter <- 5
	Pr.n.iter.tot <- length(Pr.dr.mns) * length(Pr.dr.sds) * Pr.n.iter
	JJ.samps[[tag]] <- replicate(Pr.n.iter.tot, sample(SAMPS.short,N.Samps) )
	j <- 1
	for ( mn in Pr.dr.mns ) {
		for ( sd in Pr.dr.sds ) {
			for ( i in 1:Pr.n.iter ) {
				# Set Tag
				pr.tag <- paste("m",mn,"v",sd,"i",i,sep="")
				print(paste("Running:",pr.tag))
				# Sample from Clinical Data
				print(JJ.samps[[tag]][,j])
				TAB.pr <- TAB[ TAB$ID%in%JJ.samps[[tag]][,j], ]
				print(dim(TAB.pr))
				# Set Priors
				Pr.dr <- paste("normal(",mn,",",sd,")",sep="")
				print(Pr.dr)
				temp.prior <- c( set_prior(Pr.int,class="b",coef="Intercept"),
					set_prior(Pr.dr,class="b",coef="DRUG") )
				## Run Model
				if ( i==1 ) {
					JJ.time[[tag]][[pr.tag]] <- system.time(
						JJ[[tag]][[pr.tag]] <- brm( as.formula(JJ.form[[tag]]),
						data=TAB.pr, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
						prior=temp.prior ) # , autocor=mod.cor
					)
					pr.tag.1 <- pr.tag
				}else{
					JJ.time[[tag]][[pr.tag]] <- system.time(
						JJ[[tag]][[pr.tag]] <- update( JJ[[tag]][[pr.tag.1]],
						newdata=TAB.pr, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm ) # autocor=mod.cor
					)
				}
				## Print Some Results
				temp.print <- paste(rownames(fixef(JJ[[tag]][[pr.tag]])),round(fixef(JJ[[tag]][[pr.tag]]),6),sep="=")
				print(paste("Model:",pr.tag,"- FIXEF:", temp.print ))
				print(paste("Model:",pr.tag,"- TIME:", round( JJ.time[[tag]][[pr.tag]][3], 3 ) ))
				# Next J
				j <- j + 1
			}
		}
	}

	## Pull out Fixed Effects & Distributions
	Si1.fixef.1 <- t(Reduce( cbind, lapply(JJ[[tag]],fixef) ))
	Si1.fixef.2 <- t(sapply(strsplit( gsub("[a-z]","_",names(JJ[[tag]])),"_" ),"[",2:4 ))
	Si1.fixef <- data.frame( MOD=names(JJ[[tag]]), Si1.fixef.2, Si1.fixef.1, stringsAsFactors=F )
	colnames(Si1.fixef)[2:4] <- c("Pr.Mean","Pr.SD","Iter")
	Si1.time <- data.frame( MOD=names(JJ.time[[tag]]), Reduce( rbind, JJ.time[[tag]] ), stringsAsFactors=F )
	Si1.mg <- merge( Si1.fixef,Si1.time, by="MOD" )
	Si1.mg <- Si1.mg[ order(as.numeric(Si1.mg$Pr.Mean),as.numeric(Si1.mg$Pr.SD)), ]
	 # Write Tables for Later Use
	write.table( Si1.mg, paste(PathToPlot,"TAB.",tag,"-Fixef.Time.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	# Si1.mg <- read.table( paste(PathToPlot,"TAB.",tag,"-Fixef.Time.txt",sep=""),sep="\t",header=T )

	## Summarize Fixed Effects by Mean & SD
	Si1.fef.mn <- aggregate( Si1.mg[,c("Intercept","DRUG")], by=list(Mean=Si1.mg$Pr.Mean,SD=Si1.mg$Pr.SD), mean )
	Si1.fef.sd <- aggregate( Si1.mg[,c("Intercept","DRUG")], by=list(Mean=Si1.mg$Pr.Mean,SD=Si1.mg$Pr.SD), sd )
	Si1.fef.n <- aggregate( Si1.mg[,"Intercept"], by=list(Mean=Si1.mg$Pr.Mean,SD=Si1.mg$Pr.SD), length )
	Si1.comp <- data.frame( Pr=Si1.fef.n, MN=Si1.fef.mn[,-(1:2)], SD=Si1.fef.sd[,-(1:2)], stringsAsFactors=F )
	colnames(Si1.comp)[3] <- "N.Iter"
	Si1.comp$Pr.Mean <- as.numeric(Si1.comp$Pr.Mean)
	Si1.comp$Pr.SD <- as.numeric(Si1.comp$Pr.SD)
	Si1.comp <- Si1.comp[ order(Si1.comp$Pr.Mean,Si1.comp$Pr.SD), ]

	###################
	## Plot ###########

	## Plot Posterior vs Mean/S.D.
	COLS <- COLS.list.2
	YLIM <- range( Si1.mg[,"DRUG"] ) + c(0,.3*diff(range(Si1.mg$DRUG)))
	XLIM.1 <- range( Si1.comp$Pr.Mean )
	XLIM.2 <- range( Si1.comp$Pr.SD )
	png( paste(PathToPlot,tag,".1-BetaVsSi1.png",sep=""),height=1600,width=2000,pointsize=32)
	layout( matrix(c(1:3,3),byrow=T,ncol=2) )
	 # vs Mean
	plot( 0,0,type="n", xlim=XLIM.1, ylim=YLIM,ylab="Posterior Effect Size",xlab="Prior Mean",main="Drug: Posterior Estimate vs Prior Mean" )
	abline(h=seq(-10,10,.2),lty=3,col="grey50",lwd=1 )
	abline(h=0)
	points( DRUG ~ Pr.Mean, data=Si1.mg, pch=16,col=adjustcolor(COLS,alpha=.5)[factor(Pr.SD)] )
	SCRAP <- lapply( Pr.dr.sds,function(x)points( MN.DRUG ~ Pr.Mean, data=Si1.comp[Si1.comp$Pr.SD==x,], type="o",pch=1,col=COLS[which(Pr.dr.sds==x)],cex=1.5,lwd=4 ) )
	legend("topright",pch=16,col=COLS[1:length(unique(Si1.comp$Pr.SD))],legend=levels(factor(Si1.comp$Pr.SD)),ncol=length(unique(Si1.comp$Pr.SD)),title="Prior St.Dev." )
	 # vs S.D.
	plot( 0,0,type="n", xlim=XLIM.2, ylim=YLIM,ylab="Posterior Effect Size",xlab="Prior St.Dev.",main="Drug: Posterior Estimate vs Prior St.Dev." )
	abline(h=seq(-10,10,.2),lty=3,col="grey50",lwd=1 )
	abline(h=0)
	points( DRUG ~ Pr.SD, data=Si1.mg, pch=16,col=adjustcolor(COLS,alpha=.5)[factor(Pr.Mean)] )
	SCRAP <- lapply( Pr.dr.mns,function(x)points( MN.DRUG ~ Pr.SD, data=Si1.comp[Si1.comp$Pr.Mean==x,], type="o",pch=1,col=COLS[which(Pr.dr.mns==x)],cex=1.5,lwd=4 ) )
	legend("topright",pch=16,col=COLS[1:length(unique(Si1.comp$Pr.Mean))],legend=levels(factor(Si1.comp$Pr.Mean)),ncol=length(unique(Si1.comp$Pr.Mean)),title="Prior Mean" )
	 # Boxplot
	boxplot( DRUG ~ as.numeric(Pr.Mean)+Pr.SD, data=Si1.mg, pch=16,col=rep(COLS[1:length(unique(Si1.comp$Pr.SD))],each=length(unique(Si1.comp$Pr.Mean))),xlab="Prior_Mean.Prior_SD",main="Drug: Posterior vs Prior Mean & St.Dev.",ylab="Posterior Effect Size",las=2,ylim=YLIM )
	abline( h=seq(-10,10,.2),lty=3,col="grey50",lwd=1 )
	abline( v=seq(.5,100,length(unique(Si1.comp$Pr.Mean))),lty=1,col="grey50",lwd=1)
	boxplot( DRUG ~ as.numeric(Pr.Mean)+Pr.SD, data=Si1.mg, pch=16,col=rep(COLS[1:length(unique(Si1.comp$Pr.SD))],each=length(unique(Si1.comp$Pr.Mean))),xlab="Prior_Mean.Prior_SD",main="Drug: Posterior vs Prior Mean & St.Dev.",ylab="Posterior Effect Size",las=2,ylim=YLIM,add=T )
	legend("topright",fill=COLS[1:length(unique(Si1.comp$Pr.SD))],legend=levels(factor(Si1.comp$Pr.SD)),ncol=length(unique(Si1.comp$Pr.SD)),title="Prior St.Dev." )
	TEMP.xvals <- as.numeric(as.factor(as.numeric(Si1.mg$Pr.Mean)))+length(unique(Si1.comp$Pr.Mean))*(as.numeric(as.factor(as.numeric(Si1.mg$Pr.SD)))-1)
	points( Si1.mg$DRUG ~ TEMP.xvals, pch=16,col=adjustcolor("black",.5) )
	dev.off()

	# ## Violin Plot: Prior vs Posterior
	# XLIM <- c(1,nrow(Si1.mg))
	# YLIM <- c(-6,4)
	# png( paste(PathToPlot,tag,".2-PriorVsPost.png",sep=""),height=800,width=2400,pointsize=24)
	# plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="",main="Posterior Estimate vs Prior",xaxt="n" )
	# axis( 1, at=1:nrow(Si1.mg), label=Si1.mg$MOD, las=2)
	# abline(h=-10:10,lty=3,col="grey50",lwd=1 )
	# abline(h=0)
	# SCRAP <- lapply( 1:nrow(Si1.mg), function(x)vioplot(rnorm(10000,as.numeric(Si1.mg[x,"Pr.Mean"]),as.numeric(Si1.mg[x,"Pr.SD"])),at=x, col=adjustcolor(COLS[factor(Si1.mg$Pr.SD)[x]],alpha=.6),border=NA,add=T ))
	# SCRAP <- lapply( 1:nrow(Si1.mg), function(x)vioplot(posterior_samples(JJ$Simp1[[x]])$b_DRUG,at=x, col=adjustcolor("black",alpha=.6),border=NA,add=T ))
	# legend( "topleft",fill=adjustcolor(c("black",COLS[1:length(unique(Si1.mg$Pr.SD))]),alpha=.6),border=NA,legend=c("Posterior",paste("Prior: St.Dev.=",unique(Si1.mg$Pr.SD),sep="")),ncol=length(unique(Si1.mg$Pr.SD))+1 )
	# dev.off()
	# ## Violin Plot: Prior vs Posterior (Compiled)
	# XLIM <- c(1,nrow(Si1.comp))
	# YLIM <- c(-6,4)
	# png( paste(PathToPlot,tag,".2b-PriorVsPost.png",sep=""),height=1000,width=2000,pointsize=24)
	# plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="",main="Posterior Estimate vs Prior",xaxt="n" )
	# axis( 1, at=1:nrow(Si1.comp), label=paste("m",Si1.comp$Pr.Mean,"v",Si1.comp$Pr.SD,sep=""), las=2)
	# abline(h=-10:10,lty=3,col="grey50",lwd=1 )
	# abline(h=0)
	# SCRAP <- lapply( 1:nrow(Si1.comp), function(x)vioplot(rnorm(10000,as.numeric(Si1.comp[x,"Pr.Mean"]),as.numeric(Si1.comp[x,"Pr.SD"])),at=x, col=adjustcolor(COLS[factor(Si1.comp$Pr.SD)[x]],alpha=.6),border=NA,add=T ))
	# SCRAP <- lapply( 1:nrow(Si1.comp), function(x)vioplot(unlist(lapply( JJ$Simp1[grep(paste("m",Si1.comp[x,"Pr.Mean"],"v",Si1.comp[x,"Pr.SD"],sep=""),names(JJ$Simp1))], function(y)posterior_samples(y)$b_DRUG)),at=x, col=adjustcolor("black",alpha=.6),border=NA,add=T ))
	# legend( "topleft",fill=adjustcolor(c("black",COLS[1:length(unique(Si1.comp$Pr.SD))]),alpha=.6),border=NA,legend=c("Posterior",paste("Prior: St.Dev.=",unique(Si1.comp$Pr.SD),sep="")),ncol=length(unique(Si1.comp$Pr.SD))+1 )
	# dev.off()
	 # Sort by S.D.
	Si1.comp.2 <- Si1.comp[order(Si1.comp$Pr.SD),]
	XLIM <- c(1,nrow(Si1.comp.2))
	YLIM <- c(-6,4)
	png( paste(PathToPlot,tag,".2c-PriorVsPost.png",sep=""),height=1000,width=2000,pointsize=24)
	plot( 0,0,type="n", xlim=XLIM, ylim=YLIM,ylab="Prior/Posterior Effect Size",xlab="",main="Posterior Estimate vs Prior",xaxt="n" )
	axis( 1, at=1:nrow(Si1.comp.2), label=paste("m",Si1.comp.2$Pr.Mean,"v",Si1.comp.2$Pr.SD,sep=""), las=2)
	abline(h=-10:10,lty=3,col="grey50",lwd=1 )
	abline(h=0)
	SCRAP <- lapply( 1:nrow(Si1.comp.2), function(x)vioplot(rnorm(10000,as.numeric(Si1.comp.2[x,"Pr.Mean"]),as.numeric(Si1.comp.2[x,"Pr.SD"])),at=x, col=adjustcolor(COLS[factor(Si1.comp.2$Pr.SD)[x]],alpha=.6),border=NA,add=T ))
	SCRAP <- lapply( 1:nrow(Si1.comp.2), function(x)vioplot(unlist(lapply( JJ$Simp1[grep(paste("m",Si1.comp.2[x,"Pr.Mean"],"v",Si1.comp.2[x,"Pr.SD"],sep=""),names(JJ$Simp1))], function(y)posterior_samples(y)$b_DRUG)),at=x, col=adjustcolor("black",alpha=.6),border=NA,add=T ))
	legend( "topleft",fill=adjustcolor(c("black",COLS[1:length(unique(Si1.comp.2$Pr.SD))]),alpha=.6),border=NA,legend=c("Posterior",paste("Prior: St.Dev.=",unique(Si1.comp.2$Pr.SD),sep="")),ncol=length(unique(Si1.comp.2$Pr.SD))+1 )
	dev.off()

	# ## Plot Time vs Size
	# YLIM <- range(Si1.mg$elapsed[Si1.mg$Iter!=1])
	# png( paste(PathToPlot,tag,".3-TimeVsSi1.png",sep=""),height=800,width=1600,pointsize=24)
	# par(mfrow=c(1,2))
	# plot( elapsed ~ Pr.Mean, data=Si1.mg,subset=Iter!=1, pch=16,ylim=YLIM,col=adjustcolor(COLS[factor(Pr.SD)],alpha=.7),xlab="Prior Mean",ylab="Time (s)",main="Computational Burden vs Prior Mean" )
	# abline( h=seq(0,YLIM[2]+10,10),v=seq(0,500,20),lty=3,col="grey50",lwd=1 )
	# # plot( elapsed ~ Pr.SD, data=Si1.mg,subset=Iter!=1, pch=16,ylim=YLIM,col=adjustcolor(COLS[factor(Pr.SD)],alpha=.7),xlab="Prior St.Dev.",ylab="Time (s)",main="Computational Burden vs Prior St.Dev." )
	# # abline( h=seq(0,YLIM[2]+10,10),v=seq(0,500,20),lty=3,col="grey50",lwd=1 )
	# plot( Si1.mg$elapsed[Si1.mg$Iter!=1] ~ (1:nrow(Si1.mg))[Si1.mg$Iter!=1], pch=16,ylim=YLIM,col=adjustcolor(COLS[factor(Si1.mg$Pr.SD[Si1.mg$Iter!=1])],alpha=.7),xlab="Prior St.Dev.",ylab="Time (s)",main="Computational Burden vs Prior St.Dev." )
	# abline( h=seq(0,YLIM[2]+10,10),v=seq(0,500,20),lty=3,col="grey50",lwd=1 )
	# dev.off()
}
####################################
## Large/Simple Model ##############
## Vary Mean & SD ##################
if ( grepl("Prior",Goal,ignore.case=T) & grepl("Simp2",Mod.Names.1,ignore.case=T) ) {
	write(paste(date(),"- Validate: Sample Size ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Simp2: Simple Model
	# DWAI
}
####################################
## Small/Complex Model #############
## Vary Mean #######################
if ( grepl("Prior",Goal,ignore.case=T) & grepl("Comp",Mod.Names.1,ignore.case=T) ) {
	write(paste(date(),"- Validate: Sample Size ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Comp1: Complex Model
	# DWAI
}
##############################################################
## CLINICAL TRIAL MODEL ######################################
##############################################################
if ( grepl("Clin|Trial",Goal,ignore.case=T) ) {
	write(paste(date(),"- Building Clinical Trial Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	####################################
	## Get Models Set Up ###############
	write(paste(date(),"Setting Model Parameters"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Set Priors for each Variable
	 # Build Fixed Linear Model
	MOD.fixed <- lm( DAS ~ DRUG*(WK+ACPA)+PLAC, data=TAB )
	 # Specify Prior Distribution for each Variable
	JJ.priors.list <- list()
	JJ.priors.list$Intercept <- set_prior("normal(5,1)",class="b",coef="Intercept")
	JJ.priors.list$DRUG <- set_prior("normal(-1,1)",class="b",coef="DRUG")
	JJ.priors.list$WK <- set_prior("normal(0,.1)",class="b",coef="WK")
	JJ.priors.list$PLAC <- set_prior("normal(0,1)",class="b",coef="PLAC")
	JJ.priors.list$TRT <- set_prior("normal(0,1)",class="b",coef="TRT")
	JJ.priors.list$ACPA <- set_prior("normal(0,1)",class="b",coef="ACPA")
	JJ.priors.list$DRUG_WK <- set_prior("normal(0,.1)",class="b",coef="DRUG:WK")
	JJ.priors.list$DRUG_ACPA <- set_prior("normal(0,1)",class="b",coef="DRUG:ACPA")
	JJ.priors.list$ID.sd <- set_prior("cauchy(0,2)",class="sd",group ="ID")
	JJ.priors.list$cor <- set_prior("lkj(1.5)", class = "cor")
	JJ.priors.list$HLA <- set_prior("normal(0,.2)",class="b",coef="HLA")
	JJ.priors.list$DRUG_HLA <- set_prior("normal(0,.2)",class="b",coef="DRUG:HLA")

	## CLIN10: Drug*(Week+ACPA) + Placebo + Correlation Structure + Random Intercept & Drug & Placebo
	tag <- "clin10"
	JJ.priors[[tag]] <- list()
	JJ.priors[[tag]][[paste(tag,"wk0",sep="")]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$TRT, # TRT [or] PLAC???
		JJ.priors.list$cor,
		JJ.priors.list$ID.sd,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$DRUG_ACPA )
	JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
	JJ.form[[tag]] <- "DAS ~ (1+DRUG+PLAC|ID)+DRUG*(WK+ACPA)+PLAC"

	## clin1: Drug Only (Simple Model for Testing Purposes)
	# tag <- "clin1"
	# JJ.priors[[tag]] <- c(JJ.priors.list$Intercept,
	# 	JJ.priors.list$DRUG )
	# # JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
	# JJ.form[[tag]] <- "DAS ~ DRUG"

	####################################
	## Run Models ######################

	## Specify Iterations & Samples
	# N.Samps <- 421
	# if ( !exists("JJ.samps") ) {
	JJ.samps[[tag]] <- sample(SAMPS.short,N.Samps)
	save( JJ.samps, file=paste(PathToPlot,"Rdata.Samps.Rdata",sep="") )
	# }

	## Simulate Start Dates for Patients
	Enroll.countby <- 10
	Enroll.start <- 0
	Enroll.end <- 100
	Enroll.samps <- round( runif(N.Samps,0,100) ) # sample( seq(Enroll.start, Enroll.end, Enroll.countby), N.Samps, replace=T )
	names(Enroll.samps) <- JJ.samps[[tag]]

	## Add Absolute Week to Data Table
	TAB.a <- data.frame( TAB[TAB$ID%in%JJ.samps[[tag]],], WK.ab=0 )
	for ( s in 1:N.Samps ) {
		samp <- JJ.samps[[tag]][s]
		TAB.a[ which(TAB.a$ID==samp), "WK.ab" ] <- TAB.a[ which(TAB.a$ID==samp), "WK" ] + Enroll.samps[s]
	}
	# hist( TAB.a$WK.ab, breaks=seq(0,100+Enroll.end,Enroll.countby),col="dodgerblue2" )
	write.table( TAB.a, paste(PathToPlot,"TAB.",tag,"-WKab_DataTable.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	 # Get Range of Weeks
	WKS.ab <- min(TAB.a$WK.ab):max(TAB.a$WK.ab)
	WKS.ab <- seq( min(TAB.a$WK.ab)+Enroll.countby/2, max(TAB.a$WK.ab), Enroll.countby )
	WKS.ab <- sort(unique( c(WKS.ab,max(TAB.a$WK.ab)) ))
	WKS.N <- length(WKS.ab)

	## Run Models
	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	JJ[[tag]] <- JJ.time[[tag]] <- list()
	for ( w in 1:WKS.N ) {
		## Set Params & Data
		wk <- WKS.ab[w]
		w.tag <- paste(tag,"wk",wk,sep="")
		mod.form <- as.formula( JJ.form[[tag]] )
		if ( w==1 ) {
			prev.wk <- min(TAB.a$WK.ab)-1
			mod.prior <- JJ.priors[[tag]][[paste(tag,"wk0",sep="")]]
		}else{
			prev.wk <- WKS.ab[w-1]
			mod.prior <- JJ.priors[[tag]][[paste(tag,"wk",prev.wk,sep="")]]
		}
		TAB.w <- TAB.a[ which( TAB.a$WK.ab<=wk & TAB.a$WK.ab>prev.wk ) , ]
		## Run Model
		print(paste("Running Model:",w.tag))
		if ( tag %in% names(JJ.cor) ) {
			mod.cor <- JJ.cor[[tag]]
			JJ.time[[tag]][[w.tag]] <- system.time(
				JJ[[tag]][[w.tag]] <- brm( mod.form,
				data=TAB.w, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior, autocor=mod.cor )
			)
		}else{
			JJ.time[[tag]][[w.tag]] <- system.time(
				JJ[[tag]][[w.tag]] <- brm( mod.form,
				data=TAB.w, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior )
			)
		}
		write(paste(date(),"Model",w.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		## Save Results
		temp.model <- JJ[[tag]][[w.tag]]
		save( temp.model, file=paste(PathToPlot,"Rdata.Model.",w.tag,".Rdata",sep="") )
		## Set Posteriors as New Priors for Next Model
		new.prior <- mod.prior
		temp.post <- posterior_samples(temp.model)
		temp.cov.f <- rownames(fixef(temp.model)) # grep("^b_",colnames(temp.post),value=T)
		for ( cov in temp.cov.f ) {
			cov.tag <- paste("b",cov,sep="_")
			post.mn <- round(mean( temp.post[,cov.tag] ),5)
			post.sd <- round(sd( temp.post[,cov.tag] ),5)
			cov.prior <- paste("normal(",post.mn,",",post.sd,")",sep="")
			new.prior[ which(new.prior$class=="b" & new.prior$coef==cov), "prior" ] <- cov.prior
		}
		JJ.priors[[tag]][[w.tag]] <- new.prior
		## Print Some Results
		temp.print <- paste(rownames(fixef(JJ[[tag]][[w.tag]])),round(fixef(JJ[[tag]][[w.tag]]),6),sep="=")
		print(paste("Model:",w.tag,"- FIXEF:", temp.print ))
		print(paste("Model:",w.tag,"- TIME:", round( JJ.time[[tag]][[w.tag]][3], 3 ) ))
	}
}










































# ##############################################################
# ## ASSIGN FUNCTIONS FOR LMM/HLA ANALYSES #####################
# ##############################################################

# #########################################
# ## FCT: Collapse Amino Acid Haplotypes ##

# COLLAPSE <- function( Gene, Positions ) {
# 	## Pull together Haplotypes from Specified Positions
# 	HAP <- apply( PAT_AA[[Gene]][,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
# 	HAP <- gsub("NA","-",HAP)
# 	HAP.uniq <- sort(unique(HAP))
# 	N.HAP <- length(HAP.uniq)
# 	 # Get Haplotype Frequencies
# 	HAP.freq <- table(HAP)
# 	HAP.rare <- names(HAP.freq)[which(HAP.freq < 15 )]
# 	 # Convert to Array for Additive Analysis
# 	HAP.arr <- array( 0,c(N.SAMPS,length(HAP.uniq)) )
# 	colnames(HAP.arr) <- HAP.uniq
# 	rownames(HAP.arr) <- SAMPS.short
# 	for ( pat in SAMPS.short ) {
# 		HAP.pat <- HAP[ grep(pat,names(HAP)) ]
# 		if ( HAP.pat[1]==HAP.pat[2] ) { HAP.arr[pat,HAP.pat[1]] <- 2
# 		}else{ HAP.arr[pat,HAP.pat] <- 1 }
# 	}
# 	## Compile Outputs
# 	HAP.arr <- data.frame( ID=rownames(HAP.arr), HAP.arr )
# 	OUT <- list( ARR=HAP.arr, FREQ=HAP.freq )
# 	return(OUT)
# }

# #########################################
# ## SPECIFY FUNCTION INPUTS ##############

# ## Haplotype Table
#  # Type Level
# HAP_TAB <- HLA_TYP.l
#  # Amino-Acid Level
# # HAP_TAB <- HLA_AA.l
# ## Clinical Table
# CLIN_TAB <- CLIN_TABS[["DAS"]]
# Pr_Gene <- "Pr4_DRB1"
# ## Set Formula for LMM
#  # Fixed
# FORMULA.fixed <- paste("DAS ~ DRUG*(cand+WK+RF_ACPA)+PLAC",sep="" )
# FORMULA.fixed <- paste("DAS ~ DRUG*(cand+WK+ACPA)+PLAC",sep="" )
#  # Random
# FORMULA.rand <- as.formula( "~ DRUG+WK | ID" )

# #########################################
# ## FUNCTIONS TO RUN/PLOT MODELS #########

# ## FCT: LMM vs HLA for Specific Gene
# LMM_HLA <- function( HAP_TAB, CLIN_TAB, Pr_Gene, FORMULA.fixed, FORMULA.rand ) {
# 	## Create List for Output
# 	LMM.OUT <- list()
# 	 # Add Formulae to Output
# 	LMM.OUT$Form <- list( FORM.fx=FORMULA.fixed, FORM.rn=FORMULA.rand )
# 	## Parse Input
# 	if ( length(Pr_Gene)==1 ) {
# 		# Pull out relevant columns from HAP_TAB
# 		HAP.which <- grep(Pr_Gene,colnames(HAP_TAB))
# 		HAP.names <- colnames(HAP_TAB)[HAP.which]
# 		# Add Tag to Output for Plotting Purposes
# 		LMM.OUT$Tag <- strsplit(Pr_Gene,"_")[[1]]
# 		names(LMM.OUT$Tag) <- c("Pr","Gene")
# 	}else{
# 		if ( substr(Pr_Gene["tag"],1,2)=="AA" ) {
# 			# Pull out relevant columns from HAP_TAB
# 			Pr_Gene.grep <- paste(Pr_Gene,"_",sep="")
# 			HAP.which <- Reduce( union, lapply( Pr_Gene.grep, function(x)grep(x,colnames(HAP_TAB)) ))
# 			HAP.names <- colnames(HAP_TAB)[HAP.which]
# 			# Add Tag to Output for Plotting Purposes
# 			LMM.OUT$Tag <- strsplit(Pr_Gene[["tag"]],"_")[[1]]
# 			names(LMM.OUT$Tag) <- c("Pr","Gene")
# 		}
# 		if ( substr(Pr_Gene["tag"],1,3)=="HAP" ) {
# 			# Pull out relevant columns from HAP_TAB
# 			HAP.which <- Reduce( union, lapply( Pr_Gene, function(x)grep(x,colnames(HAP_TAB)) ))
# 			HAP.names <- colnames(HAP_TAB)[HAP.which]
# 			# Add Tag to Output for Plotting Purposes
# 			LMM.OUT$Tag <- strsplit(Pr_Gene[["tag"]],"_")[[1]]
# 			names(LMM.OUT$Tag) <- c("Pr","Gene","Pos")
# 		}
# 		if ( substr(Pr_Gene["tag"],1,6)=="HAPGRP" ) {
# 			# Pull out relevant columns from HAP_TAB
# 			HAP.which <- Reduce( union, lapply( Pr_Gene, function(x)grep(x,colnames(HAP_TAB)) ))
# 			HAP.names <- colnames(HAP_TAB)[HAP.which]
# 			# Add Tag to Output for Plotting Purposes
# 			LMM.OUT$Tag <- strsplit(Pr_Gene[["tag"]],"_")[[1]]
# 			names(LMM.OUT$Tag) <- c("Pr","Gene","Pos")
# 		}	
# 	}

# 	## Pull out Haplotypes of Interest
# 	 # Calculate Frequency of each Hap
# 	FREQ.hap <- apply(HAP_TAB[,HAP.which], 2, sum)
# 	LMM.OUT$Freq <- FREQ.hap
# 	 # Common Haplotypes
# 	HAP.which.com <- HAP.which[ which(FREQ.hap>=10 & !(HAP.names%in%paste(Pr_Gene,"__",sep=""))) ]
# 	HAP.names.com <- colnames(HAP_TAB)[HAP.which.com]
# 	HAP_TAB.2 <- HAP_TAB[ ,c(1,HAP.which.com) ]

# 	## Merge Candidate SNPs w/ Clinical Variables
# 	MG <- merge( x=CLIN_TAB, y=HAP_TAB.2, by="ID" )

# 	## Test Model that includes each Variant Individually
# 	start_time <- proc.time()
# 	for ( c in 1:length(HAP.names.com) ) {
# 		cand <- HAP.names.com[c]
# 		if ( length(unique(MG[,cand]))==1 ) { print(paste( "Skipped",c,"-",cand )) ; next }
# 		formula <- as.formula( gsub("cand",cand,FORMULA.fixed) )
# 		LMM.OUT[[cand]] <- lme( fixed = formula, random = FORMULA.rand, data=MG, correlation=corCAR1(value = .5, form = ~WK | ID) )
# 		run_time <- round( proc.time()-start_time, 2)[3]
# 		if ( c%%2==0 ) { print(paste( "Done with",c,"of",length(HAP.names.com),"-",run_time )) }
# 	}
# 	return( LMM.OUT )
# }













