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
if ( is.na(Dir.Tag) ) { Dir.Tag <- "BayesLMM_HLA" }

## New Mac Paths
PathToTypes <- paste(Root,"Janssen/Data/HLA/SOAP_HLA_Types/20151211_HLA_Types.Rdata",sep="")
PathToAA <- paste(Root,"Janssen/Data/HLA/Amino_Acids/20160126_HLA_AA.Rdata",sep="")
PathToRawFiles <- paste(Root,"Janssen/Data/Pheno/Raw_Files/",sep="")
PathToAssoc <- paste(Root,"Janssen/Data/HLA/Association/20160503_HLA_Assoc_",sep="")
PathToData <- paste(Root,"Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt",sep="")
PathToFT <- paste(Root,"Janssen/Data/Pheno/Derived/20151015_Full_Table.txt",sep="")
PathToDer <- paste(Root,"Janssen/Data/Pheno/Derived/20150619_Derived_Pheno_Table.txt",sep="")
PathToNullMod <- paste(Root,"Janssen/Plots_Mac/20160517_Null_All/",sep="")
# PathToNullMod.T <- paste(Root,"/Janssen/Plots_Mac/20160523_Null_TRT_421/",sep="")
PathToNullMod.T <- paste(Root,"/Janssen/Plots_Mac/20160601_Null_TRT_421/",sep="")
PathToNullMod.Tc <- paste(Root,"/Janssen/Plots_Mac/20160601_Null_TRT_421/",sep="")
PathToHLAMod.1 <- paste(Root,"Janssen/Plots_Mac/20160519_HLA_Fr20/",sep="")
PathToPlot <- paste(Root,"Janssen/Plots_Mac/",DATE,"_",Dir.Tag,"/",sep="")
dir.create( PathToPlot )

## Previously Compiled Data
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FUL <- read.table( PathToFT,sep="\t",header=T)

## Load Previously Compiled Models & Data
# JJ <- list()
#  # Null Models
# Temp.files <- grep( "Rdata$",list.files( PathToNullMod ), value=T )
# Temp.files.2 <- grep( "Model_Summary", Temp.files, value=T,invert=T )
# Temp.files.mods <- grep("Model",Temp.files.2, value=T,invert=F )
# Temp.files.meta <- grep("Model",Temp.files.2, value=T,invert=T )
# JJ[["Null"]] <- list()
# for ( file in Temp.files.meta ) { load(paste(PathToNullMod,file,sep="")) }
# for ( file in Temp.files.mods ) {
# 	print(file)
# 	mod.tag <- gsub( "Rdata.Model.","", gsub("i1.Rdata","",file, fixed=T),fixed=T)
# 	print(mod.tag)
# 	load(paste(PathToNullMod,file,sep=""))
# 	JJ[["Null"]][[mod.tag]] <- model
# }
#  # HLA Models (1)
# if ( !exists("JJ") ) { JJ <- list() }
# Temp.files <- grep( "Rdata$",list.files( PathToHLAMod.1 ), value=T )
# Temp.files.2 <- grep( "Model_Summary", Temp.files, value=T,invert=T )
# Temp.files.mods <- grep("Model",Temp.files.2, value=T,invert=F )
# Temp.files.meta <- grep("Model",Temp.files.2, value=T,invert=T )
# for ( file in Temp.files.meta ) { load(paste(PathToHLAMod.1,file,sep="")) }
# JJ[["hla10"]] <- list()
# for ( file in Temp.files.mods ) {
# 	print(file)
# 	mod.tag <- gsub( "Rdata.Model.","", gsub(".Rdata","",file, fixed=T),fixed=T)
# 	print(mod.tag)
# 	load(paste(PathToHLAMod.1,file,sep=""))
# 	JJ[["hla10"]][[mod.tag]] <- temp.model
# }
 # Null TRT Models
if ( !exists("JJ") ) { JJ <- list() }
Temp.files <- grep( "Rdata$",list.files( PathToNullMod.T ), value=T )
Temp.files.2 <- grep( "Model_Summary", Temp.files, value=T,invert=T )
Temp.files.mods <- grep("Model",Temp.files.2, value=T,invert=F )
Temp.files.meta <- grep("Model",Temp.files.2, value=T,invert=T )
for ( file in Temp.files.meta ) { load(paste(PathToMod.T,file,sep="")) }
JJ[["Null_T"]] <- list()
for ( file in Temp.files.mods ) {
	print(file)
	mod.tag <- gsub( "Rdata.Model.","", gsub(".Rdata","",file, fixed=T),fixed=T)
	print(mod.tag)
	load(paste(PathToNullMod.T,file,sep=""))
	JJ[["Null_T"]][[mod.tag]] <- temp.model
}
 # Null TRT Models (w/ Correlation Structure)
if ( !exists("JJ") ) { JJ <- list() }
Temp.files <- grep( "Rdata$",list.files( PathToNullMod.Tc ), value=T )
Temp.files.2 <- grep( "Model_Summary", Temp.files, value=T,invert=T )
Temp.files.mods <- grep("Model",Temp.files.2, value=T,invert=F )
Temp.files.meta <- grep("Model",Temp.files.2, value=T,invert=T )
for ( file in Temp.files.meta ) { load(paste(PathToMod.Tc,file,sep="")) }
JJ[["Null_Tc"]] <- list()
for ( file in Temp.files.mods ) {
	print(file)
	mod.tag <- gsub( "Rdata.Model.","", gsub(".Rdata","",file, fixed=T),fixed=T)
	print(mod.tag)
	load(paste(PathToNullMod.Tc,file,sep=""))
	JJ[["Null_Tc"]][[mod.tag]] <- temp.model
}

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
lapply(CLIN_TABS,dim)

##############################################################
## CREATE FUNCTIONS FOR LATER USE ############################
##############################################################
write(paste(date(),"- Creating Functions ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

## FCT: Plot Model
PLOT_MOD <- function( model, tag, plot_rand ) {
	write(paste(date(),"Model:",tag), paste(PathToPlot,"Update.txt",sep=""),append=T)
	## Save Model
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
		 # Custom (quicker) "Boxplot" of Individual Posterior Distributions
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

#############################################################
## VALIDATE USE OF PRIORS ###################################
#############################################################
write(paste(date(),"- Which Goal to Pursue??"), paste(PathToPlot,"Update.txt",sep=""),append=T)
write(paste(date(),Goal), paste(PathToPlot,"Update.txt",sep=""),append=T)
print(paste("Goal:",Goal))

## Validate Robustness w/ use of Priors
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
 # How does Sample Size influence priors?
   # Larger sample size means smaller influence of priors...
   # Show this?

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
## NULL LMM MODEL ############################################
##############################################################
if ( grepl("Null",Goal,ignore.case=T) ) {
	write(paste(date(),"- Building Null Models ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

	## Set Default Tag
	tag <- "Null"
	MOD.priors <- MOD.form <- list()
	JJ.priors[[tag]] <- JJ.form[[tag]] <- list()

	####################################
	## Set Priors for Models ###########

	## Set Priors for each Variable
	JJ.priors.list <- list()
	JJ.priors.list$Intercept <- set_prior("normal(5,1)",class="b",coef="Intercept")
	JJ.priors.list$DRUG <- set_prior("normal(-1,1)",class="b",coef="DRUG")
	JJ.priors.list$WK <- set_prior("normal(0,.1)",class="b",coef="WK")
	JJ.priors.list$PLAC <- set_prior("normal(0,1)",class="b",coef="PLAC")
	JJ.priors.list$TRT <- set_prior("normal(0,1)",class="b",coef="TRT")
	JJ.priors.list$ACPA <- set_prior("normal(0,1)",class="b",coef="ACPA")
	JJ.priors.list$DRUG_WK <- set_prior("normal(0,.1)",class="b",coef="DRUG:WK")
	JJ.priors.list$DRUG_ACPA <- set_prior("normal(0,1)",class="b",coef="DRUG:ACPA")
	JJ.priors.list$ID.sd <- set_prior("cauchy(0,2)", class = "sd", group = "ID")
	JJ.priors.list$cor <- set_prior("lkj(1.5)", class = "cor")

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
	 # M6: Drug*(Week+ACPA) + Placebo + Random Intercept
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

	## M8: Drug + Week + Placebo + ACPA + Correlation Structure + Random Intercept & Drug & Placebo
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

	## M9: Drug*Week + Placebo + ACPA + Correlation Structure + Random Intercept & Drug & Placebo
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

	####################################
	## Specify Alternative Null Models #

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
		if ( grepl("dr",Goal,ignore.case=T) ) {
			tag <- paste(tag,"d",sep="")
			JJ.cor[[tag]] <- cor_ar( ~WK | ID:DRUG, p=1, cov=F )
		}
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
	save( JJ.samps, file=paste(PathToPlot,"Rdata.Samps.Rdata",sep="") )

	## Run Models
	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
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
			print(mod.form)
			# TAB.i <- TAB[ TAB$ID%in%JJ.samps[["Null"]][,j], ]
			TAB.i <- TAB[ TAB$ID%in%JJ.samps[[tag]][,i], ]
			if ( i==1 ) {
				if ( tag %in% names(JJ.cor) ) {
					mod.cor <- JJ.cor[[tag]]
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

	#########################################
	## HLA DATA #############################

	## Load HLA Types
	HLA_AA.l <- read.table(paste(PathToAssoc,"HLA_AA_Table.txt",sep=""),header=T,sep="\t" )
	HLA_TYP.l <- read.table(paste(PathToAssoc,"HLA_Types_Table.txt",sep=""),header=T,sep="\t" )
	colnames(HLA_AA.l) <- gsub(".","_",colnames(HLA_AA.l),fixed=T)
	colnames(HLA_AA.l) <- gsub("-","_",colnames(HLA_AA.l),fixed=T)
	colnames(HLA_TYP.l) <- gsub(".","_",colnames(HLA_TYP.l),fixed=T)
	colnames(HLA_TYP.l) <- gsub("-","_",colnames(HLA_TYP.l),fixed=T)

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
	JJ.priors.list$ACPA <- set_prior("normal(0,1)",class="b",coef="ACPA")
	JJ.priors.list$DRUG_WK <- set_prior("normal(0,.1)",class="b",coef="DRUG:WK")
	JJ.priors.list$DRUG_ACPA <- set_prior("normal(0,1)",class="b",coef="DRUG:ACPA")
	JJ.priors.list$ID.sd <- set_prior("cauchy(0,2)",class="sd",group ="ID")
	JJ.priors.list$cor <- set_prior("lkj(1.5)", class = "cor")
	JJ.priors.list$HLA <- set_prior("normal(0,.2)",class="b",coef="HLA")
	JJ.priors.list$DRUG_HLA <- set_prior("normal(0,.2)",class="b",coef="DRUG:HLA")

	## HLA10: Drug*(Week+ACPA) + Placebo + Correlation Structure + Random Intercept & Drug & Placebo
	tag <- "hla10"
	JJ.priors[[tag]] <- c(JJ.priors.list$Intercept,
		JJ.priors.list$DRUG,
		JJ.priors.list$WK,
		JJ.priors.list$PLAC,
		JJ.priors.list$cor,
		JJ.priors.list$ID.sd,
		JJ.priors.list$ACPA,
		JJ.priors.list$DRUG_WK,
		JJ.priors.list$DRUG_ACPA,
		JJ.priors.list$HLA,
		JJ.priors.list$DRUG_HLA )
	JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
	JJ.form[[tag]] <- "DAS ~ (1+DRUG+PLAC|ID)+DRUG*(WK+ACPA+HLA)+PLAC"

	## HLA1: Drug Only (Simple Model for Testing Purposes)
	# tag <- "hla1"
	# JJ.priors[[tag]] <- c(JJ.priors.list$Intercept,
	# 	JJ.priors.list$DRUG )
	# # JJ.cor[[tag]] <- cor_ar( ~WK | ID, p=1, cov=F )
	# JJ.form[[tag]] <- "DAS ~ DRUG"

	####################################
	## Run Models ######################

	## Specific HLA Alleles to Test
	HLA.drb.all <- grep("Pr4_DRB1",colnames(HLA_TYP),value=T)
	HLA.drb.all.freq <- apply( HLA_TYP[,HLA.drb.all], 2, sum )
	HLA.drb <- HLA.drb.all[ HLA.drb.all.freq > 20 ]
	HLA.N <- length(HLA.drb)

	## Merge Clinical Table w/ HLA Table
	TAB.hla <- merge( TAB, HLA_TYP[,c("ID",HLA.drb)], by="ID" )

	## Specify Iterations & Samples
	# N.Samps <- 421
	# if ( !exists("JJ.samps") ) {
	JJ.samps[[tag]] <- sample(SAMPS.short,N.Samps)
	save( JJ.samps, file=paste(PathToPlot,"Rdata.Samps.Rdata",sep="") )
	# }

	## Run Models
	write(paste(date(),"Model Parameters Set"), paste(PathToPlot,"Update.txt",sep=""),append=T)
	JJ[[tag]] <- JJ.time[[tag]] <- list()
	for ( h in 1:HLA.N ) {
		# h.tag <- paste(tag,"h",h,sep="")
		h.tag <- HLA.drb[h]
		mod.form <- as.formula( JJ.form[[tag]] )
		mod.prior <- JJ.priors[[tag]]
		print(paste("Running Model:",h.tag))

		TAB.h <- TAB.hla[ TAB.hla$ID%in%JJ.samps[[tag]], ]
		colnames(TAB.h)[grep(h.tag,colnames(TAB.h))] <- "HLA"
		if ( tag %in% names(JJ.cor) ) {
			mod.cor <- JJ.cor[[tag]]
			JJ.time[[tag]][[h.tag]] <- system.time(
				JJ[[tag]][[h.tag]] <- brm( mod.form,
				data=TAB.h, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior, autocor=mod.cor )
			)
		}else{
			JJ.time[[tag]][[h.tag]] <- system.time(
				JJ[[tag]][[h.tag]] <- brm( mod.form,
				data=TAB.h, family=Mod.fam, chains=Mod.chain, iter=Mod.iter, warmup=Mod.warm,
	            prior=mod.prior )
			)
		}
		write(paste(date(),"Model",h.tag,"Done"), paste(PathToPlot,"Update.txt",sep=""),append=T)
		temp.model <- JJ[[tag]][[h.tag]]
		save( temp.model, file=paste(PathToPlot,"Rdata.Model.",h.tag,".Rdata",sep="") )
		## Print Some Results
		temp.print <- paste(rownames(fixef(JJ[[tag]][[h.tag]])),round(fixef(JJ[[tag]][[h.tag]]),6),sep="=")
		print(paste("Model:",h.tag,"- FIXEF:", temp.print ))
		print(paste("Model:",h.tag,"- TIME:", round( JJ.time[[tag]][[h.tag]][3], 3 ) ))
	}

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
##############################################################
## PLOT MODEL SUMMARIES ######################################
##############################################################
write(paste(date(),"- Plotting Model Summaries ######"), paste(PathToPlot,"Update.txt",sep=""),append=T)

## Plot/Save Models
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

##############################################################
## END OF DOC ################################################
##############################################################



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













