## Bayesian LMM Analysis of Longitudinal Janssen Data ##
## Use: brsm - Bayesian Regression Models w/ Stan Package ##
## May 6, 2016 ##
## Kristopher Standish ##

## Game Plan:
 # 

###############################################
## LOAD DATA ##################################
###############################################

## Load Packages
# library(bsts)
library(brms)
library(gplots)
# library(BLR)
# library(rstanarm)
# library(lattice)

## Set Date
DATE <- gsub("-","",Sys.Date())

## New Mac Paths
PathToData <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_LongitBRSM/",sep="" )
dir.create( PathToPlot )

## Load Janssen Data Sets
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FT.l <- read.table( PathToFT,sep="\t",header=T)

## Get unique Weeks for plotting purposes
WKS <- unique( TAB.l$WK )

##############################################################
## FILTER DATA ###############################################
##############################################################

#########################################
## CLINICAL DATA ########################

## Add Column for "TRT" (Treatment: All but Week 0)
TRT <- rep(1,nrow(TAB.l))
TRT[which(TAB.l$WK==0)] <- 0
TAB <- data.frame( TAB.l, TRT )

## Add ID Column w/ Short Hand
TAB.id.short <- sapply( strsplit( TAB$IID, "-" ),"[",1 )
TAB <- data.frame( ID=TAB.id.short, TAB[,-1] )
SAMPS <- as.character( unique(TAB$ID) )
SAMPS.G <- as.character( FT.l$ID[FT.l$GRP=="G"] )
SAMPS.P <- as.character( FT.l$ID[FT.l$GRP!="G"] )
SAMPS.EE <- as.character( FT.l$ID[FT.l$GRP=="PE"] )
SAMPS.NE <- as.character( FT.l$ID[FT.l$GRP=="P"] )

## Change ACPA to 1/0
TAB$ACPA <- as.numeric(TAB$ACPA=="Positive")
TAB$RF <- as.numeric(TAB$RF=="Positive")
TAB$RF_ACPA <- as.numeric(TAB$RF_ACPA=="Positive")

## Take out Patients who left before getting DRUG
# RM.exit.id <- as.character( FUL$ID_2[which(FUL$IN<=4)] )
# RM.exit <- which( TAB$FID %in% RM.exit.id )
# TAB <- TAB[ -RM.exit, ]
# SAMPS <- unique(as.character(TAB$FID))
# SAMPS.short <- unique(as.character(TAB$ID))
# N.SAMPS <- length(SAMPS)

# ## FCT: Create Clin Table for each Response Phenotype
# MAKE_CLIN_TAB <- function( RESP_PHENO ) {
# 	## Pull Response Phenotype of Interest
# 	TEMP <- TAB[ , c(1:grep("PLAC",colnames(TAB)),which(colnames(TAB)==RESP_PHENO)) ]
# 	## Remove DAS values that are NA
# 	RM.na <- which(is.na( TEMP[,RESP_PHENO] ))
# 	if (length(RM.na) > 0 ) { TEMP <- TEMP[-RM.na, ] }
# 	## Return Table
# 	return(TEMP)
# }
#  # Make Clinical Table for each Phenotype
# RESP_PHENOS <- c("DAS","lCRP","rSJC","rTJC")
# CLIN_TABS <- lapply( RESP_PHENOS, function(x)MAKE_CLIN_TAB(x) )
# names(CLIN_TABS) <- RESP_PHENOS
# lapply(CLIN_TABS,dim)

##############################################################
## BUILD/COMPARE MODELS ######################################
##############################################################

## Package Parameters
 # From Tutorial: http://www.r-bloggers.com/r-users-will-now-inevitably-become-bayesians/
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores


## Specify Table of Samples
TAB.j <- TAB[ which(TAB$ID%in%SAMPS.P), ]
TAB.j <- TAB[ which(TAB$ID%in%SAMPS.P), ]

JJ <- JJ.priors <- list()
JJ.time <- c()

## M1: Drug Only
tag <- "m1"
JJ.priors[[tag]] <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG") )
JJ.time[[tag]] <- system.time( JJ[[tag]] <- brm(DAS ~ DRUG, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600, prior=JJ.priors[[tag]] ) )

## M2: Drug + Week
tag <- "m2"
JJ.priors[[tag]] <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK") )
JJ.time[[tag]] <- system.time( JJ[[tag]] <- brm(DAS ~ DRUG+WK, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors[[tag]] ) )

## M3: Drug + Week + Placebo
tag <- "m3"
JJ.priors[[tag]] <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.time[[tag]] <- system.time( JJ[[tag]] <- brm(DAS ~ DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors[[tag]] ) )

## M4: Drug + Week + Placebo + Correlation Structure
tag <- "m4"
JJ.priors[[tag]] <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | ID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.time[[tag]] <- system.time( JJ[[tag]] <- brm(DAS ~ DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors[[tag]],autocor=JJ.cor  ) )

## M5: Drug + Week + Placebo + Correlation Structure + Random Intercept
tag <- "m5"
JJ.priors[[tag]] <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | ID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.time[[tag]] <- system.time( JJ[[tag]] <- brm(DAS ~ (1|ID)+DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors[[tag]],autocor=JJ.cor  ) )

## M6: Drug + Week + Placebo + Correlation Structure + Random Intercept & Drug
tag <- "m6"
JJ.priors[[tag]] <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | ID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.time[[tag]] <- system.time( JJ[[tag]] <- brm(DAS ~ (1+DRUG|ID)+DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors[[tag]],autocor=JJ.cor  ) )

## M7: Drug + Week + Placebo + Correlation Structure + Random Intercept & Drug
tag <- "m7"
JJ.priors[[tag]] <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | ID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.time[[tag]] <- system.time( JJ[[tag]] <- brm(DAS ~ (1+DRUG+PLAC|ID)+DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors,autocor=JJ.cor  )







## FCT: Plot Model
PLOT_MOD <- function( model, tag ) {	
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
	m.form
	RAND <- length(summ$random)>0

	## Collect Model Outputs
	d.prior <- model$prior
	d.post <- posterior_samples(model)
	f.eff <- summ$fixed
	m.covs.f <- rownames(f.eff)
	c.eff <- summ$cor_pars
	s.eff <- summ$spec_pars
	d.post.f.which <- paste("b_",m.covs.f,sep="")
	d.post.c.which <- rownames(c.eff) # grep(rownames(c.eff),colnames(d.post),value=T,fixed=T)
	if ( RAND==T ) {
		r.grp <- summ$group
		r.ngrps <- summ$ngrps
		r.eff <- ranef(model)[[1]]
		r.eff.2 <- summ$random
		m.covs.r <- colnames(r.eff)[1:ncol(r.eff)]
		m.samps <- rownames(r.eff)
		all.eff <- list( f.eff, c.eff, Reduce( rbind, r.eff.2 ), s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","R","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
		d.post.r.which <- unlist(lapply( m.covs.r, function(x)paste( "sd_",r.grp,"_",x,sep="") ))
		d.post.r2.which <- unlist(lapply( m.covs.r, function(x)paste( "r_",r.grp,"[",m.samps,",",x,"]",sep="") ))
	}else{
		all.eff <- list( f.eff, c.eff, s.eff )
		mod.fit <- Reduce( rbind, all.eff )
		mod.fit <- data.frame( mod.fit, Eff=rep(c("F","C","S"),lapply(all.eff,nrow)), stringsAsFactors=F )
	}

	## FIXED PLOTS ############
	 # Set Color Palette
	BLEND <- function( colors ) { colorRampPalette(colors)(3)[2] }
	COLS.list <- c("firebrick2","chocolate2","gold1","chartreuse2","cadetblue2","dodgerblue2","slateblue2","magenta2")
	COLS.list.heat <- c("firebrick3","chocolate2","gold1","springgreen1","steelblue2","slateblue3")
	COLS.list.2 <- c("aquamarine3","deeppink2",BLEND(c("gold1","chartreuse2")),BLEND(c("steelblue3","slateblue3")),"tomato2","deepskyblue1",BLEND(c("sienna2","goldenrod1")) )
	COLS <- adjustcolor(COLS.list.2,alpha=.6)
	COLS.eff <- COLS[1:length(unique(mod.fit[,"Eff"]))]
	names(COLS.eff) <- unique(mod.fit[,"Eff"])

	## Plot Model Summary
	png( paste(PathToPlot,"ModSumm_",tag,".1-PlotFct.png",sep=""),height=200*nrow(mod.fit),width=800,pointsize=30)
	par(mfrow=c(nrow(mod.fit),2))
	plot( model, N=nrow(mod.fit) )
	dev.off()

	## Fixed Effect Sizes
	YLIM <- extendrange(mod.fit[,"Estimate"], f=.2)
	png( paste(PathToPlot,"ModSumm_",tag,".2-EffSize.png",sep=""),height=1000,width=400+100*nrow(mod.fit),pointsize=30)
	TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n",ylim=YLIM,ylab="Effect Size" )
	axis( 1, at=TEMP,label=rownames(mod.fit), las=2 )
	axis( 2, at=seq(-10,10,2), las=2 )
	abline( h=-10:10, lty=3,col="grey50",lwd=1 )
	abline( h=0, lty=1,col="black",lwd=1 )
	TEMP <- barplot( mod.fit[,"Estimate"],col=COLS.eff[mod.fit[,"Eff"]],border=NA,yaxt="n", add=T)
	arrows( TEMP,mod.fit[,"l.95..CI"],TEMP,mod.fit[,"u.95..CI"],lwd=3,length=0 )
	legend("topright",fill=COLS.eff,border=NA,legend=names(COLS.eff),ncol=length(COLS.eff),title="Effect Type",bg="white")
	dev.off()

	## RANDOM PLOTS ###########
	if ( RAND==T ) {
		n.covs.r <- length(m.covs.r)
	
		## Plot Pairs of Random Effects
		if ( n.covs.r > 1 ) {
			temp.reff <- t( f.eff[m.covs.r,"Estimate"] + t(r.eff) ) 
			png( paste(PathToPlot,"ModSumm_",tag,".R1-PairRand.png",sep=""),height=400*n.covs.r,width=400*n.covs.r,pointsize=24+2*n.covs.r)
			pairs( temp.reff, pch=16,col=COLS.eff["R"],cex=2 )
			dev.off()
		}

		## Boxplot of Individual Estimates across Sims (d.post)
		png( paste(PathToPlot,"ModSumm_",tag,".R2-BoxRand.png",sep=""),height=600*n.covs.r,width=500+15*length(m.samps),pointsize=24+2*n.covs.r)
		par(mfrow=c(n.covs.r,1))
		par(mar=c(6,4,4,1))
		for ( z in 1:n.covs.r ) {
			z.name <- m.covs.r[z]
			z.col.1 <- grep(z.name,d.post.f.which,value=T)
			z.col.2 <- grep(z.name,d.post.r2.which,value=T)
			z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
			if ( z>0 ) { z.ord <- order( colMeans(z.temp) ) }
			boxplot( z.temp[,z.ord], col=COLS.eff["R"],main=paste("Random",z.name,"by Patient"),ylab=z.name,names=m.samps[z.ord],las=2,pch=16 )
			abline( h=-10:10, lty=3,col="grey50",lwd=1 )
			abline( h=0, lty=1,col="black",lwd=1 )
			boxplot( z.temp[,z.ord], col=COLS.eff["R"],names=m.samps[z.ord],las=2,pch=16,add=T )
		}
		dev.off()

		## Plot Posterior Probabilities for Random Effects
		COLS.heat <- colorRampPalette(c(COLS.list.heat[1:3],"white",COLS.list.heat[4:6]))(100)
		BRKS.heat <- seq(0,1,length.out=101)
		n.iter.tot <- m.chain * ( m.iter - m.warm )
		post.prob <- list()
		post.prob.brks <- list( Intercept=seq(3,8,.25),
			DRUG=seq(-3,0,.125),
			PLAC=seq(-2,1,.125) )
		for ( z in 1:n.covs.r ) {
			z.name <- m.covs.r[z]
			z.col.1 <- grep(z.name,d.post.f.which,value=T)
			z.col.2 <- grep(z.name,d.post.r2.which,value=T)
			z.temp <- d.post[,z.col.1] + d.post[,z.col.2]
			post.prob[[z.name]] <- Reduce( cbind, lapply( post.prob.brks[[z.name]],function(x)apply(z.temp,2,function(y)length(which(y<x))) ) ) / n.iter.tot
			colnames(post.prob[[z.name]]) <- paste("LT",post.prob.brks[[z.name]],sep="_")
			rownames(post.prob[[z.name]]) <- colnames(z.temp)
			if ( z.name=="Intercept" ) { post.prob[[z.name]] <- 1 - post.prob[[z.name]] ; colnames(post.prob[[z.name]]) <- gsub("LT","GT",colnames(post.prob[[z.name]])) }
			## Plot Heatmap
			png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.",z.name,".png",sep=""),height=1200,width=2000+10*length(m.samps),pointsize=28)
			heatmap.2( t(post.prob[[z.name]]), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5) )
			dev.off()
		}
		post.prob.temp <- Reduce( cbind, post.prob )
		png( paste(PathToPlot,"ModSumm_",tag,".R3-PostHeat.All.png",sep=""),height=1200,width=2000+10*length(m.samps),pointsize=28)
		heatmap.2( t(post.prob.temp), col=COLS.heat,breaks=BRKS.heat,scale="none",trace="none",Rowv=NA,dendrogram="column",lwid=c(1,5),rowsep=cumsum(lapply(post.prob,ncol)) )
		dev.off()

		## Plot Distribution (across patients) of Standard Deviations (across iteration) of Random Effects
		TEST <- lapply(m.covs.r,function(r)unlist(lapply( m.samps, function(s) sd(d.post[,paste("r_ID[",s,",",r,"]",sep="")]) )) )
		names(TEST) <- m.covs.r

		# ## Plot a Few Individual Patients over Time
		# z.samps.n <- 12
		# z.samps <- sample(m.samps, z.samps.n )
		# z.data <- model$data[model$data[,r.grp]%in%z.samps,]
		# z.temp <- predict( model, newdata=z.data )
		# # z.temp <- predict( model, newdata=z.data,re_formula=NA )
		# z.pred <- cbind( z.data, z.temp )
		# z.ranef <- t(f.eff[m.covs.r,"Estimate"] + t(r.eff[z.samps,]))
		# PRED_v_REAL <- function( sample, pred.tab ) {
		# 	COLS.pred <- COLS.list.2[c(1,3,4,6,7,3,5)]
		# 	sub.tab <- pred.tab[ pred.tab[,r.grp]==sample, ]
		# 	plot( DAS ~ WK, data=sub.tab,type="p",pch=c(16,4)[factor(PLAC)],col=COLS.pred[factor(DRUG)],xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main=sample )
		# 	abline( h=0:10,lty=3,col="grey50",lwd=1)
		# 	points( Estimate ~ WK, data=sub.tab,type="l",lwd=3 )
		# 	polygon( c(sub.tab$WK,rev(sub.tab$WK)), c(sub.tab[,"2.5%ile"],rev(sub.tab[,"97.5%ile"])), col=adjustcolor(COLS.pred[3],alpha=.2),border=NA )
		# 	points( DAS ~ WK, data=sub.tab,type="p",pch=c(16,4)[factor(PLAC)],col=COLS.pred[factor(DRUG)],lwd=2 )
		# 	points( DAS ~ WK, data=sub.tab[1,],type="p",pch=16,col=COLS.pred[4] )
		# 	pre.post <- cumsum(unlist( FT.l[FT.l$ID==sample,c("DAS_BL_MN","DEL_MNe_MN")] ))
		# 	arrows( c(0,24),pre.post,c(24,100),pre.post, col=COLS.pred[5],lwd=4,lty=2,length=0 )
		# 	pre.post.2 <- z.ranef[sample,c(1,3,2)] + c(0,rep(z.ranef[sample,1],2))
		# 	arrows( c(0,2,24),pre.post.2,c(2,24,100),pre.post.2, col=COLS.pred[c(4,1,2)],lwd=4,lty=2,length=0 )
		# 	text( 0,8, paste(round(z.ranef[sample,1],2),collapse="\n"), pos=4,col=COLS.pred[4] )
		# 	text( 0,2, paste(round(z.ranef[sample,3],2),collapse="\n"), pos=4,col=COLS.pred[1] )
		# 	text( 100,8, paste(round(z.ranef[sample,2],2),collapse="\n"), pos=2,col=COLS.pred[2] )
		# }
		# png( paste(PathToPlot,"ModSumm_",tag,".R4-Pred.png",sep=""),height=1000,width=1600,pointsize=30)
		# par(mfrow=c(3,4))
		# par(mar=c(5,4,1,1))
		# PTS <- lapply( z.samps, function(x)PRED_v_REAL(x,z.pred) )
		# dev.off()
	}
}

SCRAP <- lapply( names(JJ),function(x)PLOT_MOD(JJ[[x]],x) )
SCRAP <- lapply( names(JJ)[5:7],function(x)PLOT_MOD(JJ[[x]],x) )



PLOT_MOD( JJ$m6, "m6" )
PLOT_MOD( JJ$m7, "m7" )
LOO( JJ$m1, JJ$m2 )
LOO( JJ$m6, JJ$m7 )



NEW <- data.frame( DRUG=rep(0:1,c(21,80)), PLAC=rep(c(0,1,0),c(1,20,80)), WK=0:100 )
NEW.pred <- predict( JJ$m3, newdata=NEW )
plot( NEW$WK, NEW.pred[,1], col=COLS.pred[1+NEW$DRUG+2*NEW$PLAC],ylim=c(1,9),type="p" )
polygon( c(NEW$WK,rev(NEW$WK)), c(NEW.pred[,"2.5%ile"],rev(NEW.pred[,"97.5%ile"])), col=adjustcolor(COLS.pred[3],alpha=.2),border=NA )

PRED_NEW <- function( model, new_data ) {
	temp.wks <- sort(unique(new_data$WK))
	temp.in <- new_data[ new_data$WK%in%temp.wks[1:2], ]
	temp.in <- new_data[ new_data$WK==temp.wks[1], ]
	temp.in <- new_data
	# temp.comp <- list()
	for ( wk in temp.wks[-1] ) {
		whick <- which(temp.in$WK < wk)
		temp.out <- predict( model, newdata=temp.in[whick,] )
		temp.in$DAS[whick] <- temp.out[,"Estimate"]
	}
	temp.out <- predict( model, newdata=temp.in )
	temp.in$DAS <- temp.out[,"Estimate"]
	if ( wk==tail(temp.wks,1) ) { temp.in <- cbind( temp.in,temp.out) }
	return(temp.in)
}

NEW <- data.frame( DRUG=rep(0:1,c(21,80)), PLAC=rep(c(0,1,0),c(1,20,80)), WK=0:100, DAS=rep(c(0,NA),c(1,100)), ID=z.samps[3] )
NEW <- Reduce( rbind, lapply(z.samps,function(x) data.frame( DRUG=rep(0:1,c(21,80)), PLAC=rep(c(0,1,0),c(1,20,80)), WK=0:100, DAS=rep(c(0,NA),c(1,100)), ID=x ) ))
NEW.pred <- PRED_NEW( JJ$m7, NEW[seq(1,101,1),] )
NEW.pred <- PRED_NEW( JJ$m7, NEW[NEW$WK%in%seq(0,100,10),] )
# plot( NEW.pred$IN$WK, NEW.pred$IN$Estimate, col=COLS.pred[1+NEW.pred$IN$DRUG+2*NEW.pred$IN$PLAC],ylim=c(1,9),type="p" )
# polygon( c(NEW.pred$IN$WK,rev(NEW.pred$IN$WK)), c(NEW.pred[,"2.5%ile"],rev(NEW.pred[,"97.5%ile"])), col=adjustcolor(COLS.pred[3],alpha=.2),border=NA )
plot( NEW.pred$WK, NEW.pred$Estimate, col=COLS.pred[1+NEW.pred$DRUG+2*NEW.pred$PLAC],type="o",pch=16 )
polygon( c(NEW.pred$WK,rev(NEW.pred$WK)), c(NEW.pred[,"2.5%ile"],rev(NEW.pred[,"97.5%ile"])), col=adjustcolor(COLS.pred[3],alpha=.2),border=NA )

PLOT.samp <- function(sample) {
	plot( NEW.pred[NEW.pred$ID==sample,]$WK, NEW.pred[NEW.pred$ID==sample,]$Estimate, col=COLS.pred[1+NEW.pred[NEW.pred$ID==sample,]$DRUG+2*NEW.pred[NEW.pred$ID==sample,]$PLAC],type="o",pch=16,ylim=c(1,9) )
	polygon( c(NEW.pred[NEW.pred$ID==sample,]$WK,rev(NEW.pred[NEW.pred$ID==sample,]$WK)), c(NEW.pred[NEW.pred$ID==sample,"2.5%ile"],rev(NEW.pred[NEW.pred$ID==sample,"97.5%ile"])), col=adjustcolor(COLS.pred[3],alpha=.2),border=NA )
}
# plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9) )
par(mfrow=c(3,4))
lapply( z.samps, PLOT.samp )

Reduce(rbind, lapply(NEW.pred$COMP,function(x)head(x,1)) )



Reduce(rbind, lapply(temp.comp,function(x)head(x,1)) )


NEW.pred <- replicate( 10, PRED_NEW( JJ$m4, NEW[seq(1,101,5),] ) )

