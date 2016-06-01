## Learn Bayesian Regression Packages ##
## R Packages: bsts, brms, BLR ##
## May 5, 2016 ##

## Load Packages
library(bsts)
library(brms)
library(BLR)
library(rstanarm)
library(lattice)

###############################################
## LOAD DATA ##################################
###############################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## New Mac Paths
PathToData <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20160112_Resp_v_Time.txt"
PathToFT <- "/Users/kstandis/Data/Janssen/Data/Pheno/Derived/20151015_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Janssen/Plots_Mac/",DATE,"_BayesLong_Test/",sep="" )
# dir.create( PathToPlot )

## Load Janssen Data Sets
TAB.l <- read.table( PathToData, sep="\t",header=T, stringsAsFactors=F )
FT.l <- read.table( PathToFT,sep="\t",header=T)

## Pull out only PBO arms
Samps.PBO <- as.character( FT.l$ID_2[FT.l$GRP=="P"] )
TAB <- TAB.l[ TAB.l$IID%in%Samps.PBO, ]
TAB.2 <- TAB.l[ TAB.l$IID%in%Samps.PBO[1:20], ]
TAB.3 <- TAB.l[ TAB.l$IID%in%Samps.PBO[1:40], ]
TAB.test <- TAB.l[ TAB.l$IID%in%tail(Samps.PBO,10), ]

###############################################
## brms: Bayesian Regression Models w/ Stan ###
###############################################

############################
## Tutorial from: http://www.r-bloggers.com/r-users-will-now-inevitably-become-bayesians/
 # Package Parameters
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores
 # Data & Seed
set.seed (3875)
ir <- data.frame (scale (iris[, -5]), Species=iris[, 5])

## With improper prior it takes about 12 minutes, with about 40% CPU utilization and fans running,
 # so you probably don't want to casually run the next line...
 # NOTES on COMMAND:
   # family="categorical" - multiple pairwise logistic regression to determine class (baseline class vs each)
   # chains=3 - Number of Markov chains (Default=4)
   # iter=3000 - Number of iterations per chain
   # warmup=600 - Number of "burn in" iterations (not to be used for analysis)
# system.time (b1 <- brm (Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width, data=ir, family="categorical", n.chains=3, n.iter=3000, n.warmup=600)) # deprecated
# system.time (b1 <- brm (Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width, data=ir, family="categorical", chains=3, iter=3000, warmup=600))

## Try again, but with informative Priors
system.time (b2 <- brm (Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width, data=ir,
                        family="categorical", chains=3, iter=3000, warmup=600,
                        prior=c(set_prior ("normal (0, 8)"))))
plot( b2 )
pairs( b2 )

############################
## Tutorial: http://www.r-bloggers.com/bayesian-regression-models-using-stan-in-r/
 # Specify Data
temp <- c(11.9,14.2,15.2,16.4,17.2,18.1,18.5,19.4,22.1,22.6,23.4,25.1)
units <- c(185L,215L,332L,325L,408L,421L,406L,412L,522L,445L,544L,614L)
log_units <- log(units)
n <- length(units)
market.size <- rep(800, n)

## Create Bayesian Regression Model
 # Linear
lin.mod <- brm(units ~ temp, family="gaussian")
 # Log-transformed Linear Gaussian model
log.lin.mod <- brm(log_units ~ temp, family="gaussian")
 # Poisson model
pois.mod <- brm(units ~ temp, family="poisson")
 # Binomial model
bin.mod <- brm(units | trials(market.size) ~ temp, family="binomial")

## Plot Summaries
plot( lin.mod )
plot( log.lin.mod )
plot( pois.mod )
plot( bin.mod )

## Compile Model Results
modelData <- data.frame(
  Model=factor(c(rep("Linear model", n), 
                 rep("Log-transformed LM", n),
                 rep("Poisson (log)",n),
                 rep("Binomial (logit)",n)),  
               levels=c("Linear model", 
                        "Log-transformed LM",
                        "Poisson (log)",
                        "Binomial (logit)"), 
               ordered = TRUE),
  Temperature=rep(temp, 4),
  Units_sold=rep(units, 4),
  rbind(predict(lin.mod),
        exp(predict(log.lin.mod) + 
              0.5 * mean(extract(log.lin.mod$fit)[["sigma_log_units"]])),
        predict(pois.mod),
        predict(bin.mod)
  )
)
head(modelData)
colnames(modelData)[6:7] <- c("l.95..CI","u.95..CI")

## Plot Predictions from each Model
key <- list(
  rep=FALSE, 
  lines=list(col=c("black", "firebrick2"), type=c("p","l"), pch=1),
  text=list(lab=c("Observation","Estimate")),
  rectangles = list(col=adjustcolor("chartreuse2", alpha.f=0.5), border="grey"),
  text=list(lab="95% Prediction credible interval"))
xyplot(l.95..CI + u.95..CI + Estimate + Units_sold ~ Temperature | Model, 
       data=modelData, as.table=TRUE, main="Ice cream model comparision",
       xlab="Temperatures (C)", ylab="Units sold", 
       scales=list(alternating=1), key=key,
       panel=function(x, y){
         n <- length(x)
         k <- n/2
         upper <- y[(k/2+1):k]
         lower <- y[1:(k/2)]
         x <- x[1:(k/2)]
         panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                       col = adjustcolor("chartreuse2", alpha.f = 0.5), 
                       border = "grey")
         panel.lines(x, y[(k+1):(k+n/4)], col="firebrick2")
         panel.points(x, y[(n*3/4+1):n], lwd=2, col="black")
       })

## Predict Max (97.5%ile) Sales for Hot Day
A <- function(samples){
  as.matrix(samples[,c("b_Intercept" ,"b_temp")])
}
x <- c(1, 35)
prob <- 0.975 
lin.samples <- posterior_samples(lin.mod)
n <- nrow(lin.samples)
mu <- A(lin.samples) %*% x
sigma <- lin.samples[,"sigma_units"]
LIN <- rnorm(n,mu,sigma)
(lin.q <- quantile(LIN, prob))
#    97.5% 
# 1025.077
log.lin.samples <- posterior_samples(log.lin.mod)
mu <- A(log.lin.samples) %*% x
sigma <- log.lin.samples[,"sigma_log_units"]
LOG.LIN <- exp(rnorm(n, mu +  0.5*sigma^2, sigma))
(log.lin.q <- quantile(LOG.LIN, prob))
#    97.5% 
# 2460.108
pois.samples <- posterior_samples(pois.mod)
mu <- exp(A(pois.samples) %*% x)
POIS <- rpois(n, mu)
(pois.q <- quantile(POIS , prob))
#  97.5% 
#   1515 
bin.samples <- posterior_samples(bin.mod)
inv.logit <- function(u) exp(u)/(1+exp(u))
mu <- inv.logit( A(bin.samples) %*% x)
BIN <- rbinom(n, size = 800, mu)
(bin.q <- quantile(BIN, prob))
#   97.5% 
# 761.025
percentiles <- c(lin.q, log.lin.q, pois.q, bin.q)
b <- barplot(percentiles, 
             names.arg = c("Linear", "Log-transformed",
                           "Poisson", "Binomial"),
             ylab="Predicted ice cream units",
             main="Predicted 97.5%ile at 35ºC",
             ylim=c(0,max(percentiles)*1.2),col="royalblue2")
text(b, percentiles-75, round(percentiles))


boxplot( LIN, LOG.LIN, POIS, BIN, xaxt="n",at=b, add=T,col="chartreuse2")
points( rep(b,each=n), c(LIN,LOG.LIN,POIS,BIN),pch=16,col=adjustcolor("black",alpha=.1) )

############################
## Example from Publication

## Time to recurrence of infection
data("kidney")
fit1 <- brm(formula = time | cens(censored) ~ age + sex + disease
            + (1 + age|patient),
            data = kidney, family = gaussian("log"),
            prior = c(set_prior("normal(0,5)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, iter = 2000, chains = 4,
            control = list(adapt_delta = 0.95))

SUMM <- summary(fit1,waic=T)
SUMM$fixed
SUMM$random
## To test whether SD of "age" is smaller than the SD of Intercept, use hypothesis method:
hypothesis(fit1, "Intercept - age > 0", class = "sd", group = "patient")
 # Hypothesis being tested (Int-age>0) is ~90x more likely than alternative (Int-age<0)
 # Should "age" be included as Random Effect?

fit2 <- brm(formula = time | cens(censored) ~ age + sex + disease
	+ (1|patient),
    data = kidney, family = gaussian("log"),
    prior = c(set_prior("normal(0,5)", class = "b"),
              set_prior("cauchy(0,2)", class = "sd")),
    warmup = 1000, iter = 2000, chains = 4,
    control = list(adapt_delta = 0.95))
## Compare with "Leave One Out" cross-validation (LOO)
LOO(fit1, fit2)

############################
## Play w/ Janssen Data
TAB.j <- TAB.3
## Fixed Models (No Random Effects)
 # Drug Only
JJ.priors <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG") )
JJ.1 <- brm(DAS ~ DRUG, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600, prior=JJ.priors )
 # Drug + Week
JJ.priors <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK") )
JJ.2 <- brm(DAS ~ DRUG+WK, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors )
 # Drug + Week + Placebo
JJ.priors <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.3 <- brm(DAS ~ DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors )
 # Drug + Week + Placebo + Correlation Structure
JJ.priors <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | IID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.4 <- brm(DAS ~ DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors,autocor=JJ.cor  )
 # Drug + Week + Placebo + Correlation Structure + Random Intercept
JJ.priors <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | IID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.5 <- brm(DAS ~ (1|IID)+DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors,autocor=JJ.cor  )
 # Drug + Week + Placebo + Correlation Structure + Random Intercept & Drug
JJ.priors <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | IID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.6 <- brm(DAS ~ (1+DRUG|IID)+DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors,autocor=JJ.cor  )
 # Drug + Week + Placebo + Correlation Structure + Random Intercept & Drug
JJ.priors <- c( set_prior("normal(5,4)",class="b",coef="Intercept"),
	set_prior("normal(0,4)",class="b",coef="DRUG"),
	set_prior("normal(0,.5)",class="b",coef="WK"),
	set_prior("normal(0,4)",class="b",coef="PLAC") )
JJ.cor <- cor_ar( ~WK | IID, p=1, cov=F ) # corCAR1(value = .5, form = ~WK | ID)
JJ.7 <- brm(DAS ~ (1+DRUG+PLAC|IID)+DRUG+WK+PLAC, data=TAB.j, family="gaussian", chains=3, iter=3000, warmup=600,
                        prior=JJ.priors,autocor=JJ.cor  )


# plot( JJ.6, theme=scale_color_brewer() )
JJ.p <- JJ.6
PRED.tab <- predict( JJ.p, newdata=TAB.j )
PRED.tab <- fitted.values( JJ.p )
PRED <- cbind( TAB.j[,c("IID","DAS","WK","DRUG","PLAC")], PRED.tab )
RANEF <- t(fixef(JJ.p)[1:2,] + t(ranef(JJ.p)[[1]]))

# plot( JJ.6, theme=scale_color_brewer() )
JJ.p <- JJ.7
PRED.tab <- predict( JJ.p, newdata=TAB.j )
PRED <- cbind( TAB.j[,c("IID","DAS","WK","DRUG","PLAC")], PRED.tab )
RANEF <- t(fixef(JJ.p)[1:3,] + t(ranef(JJ.p)[[1]]))





# plot( JJ.4 )
# pairs( JJ.4 )
# PRED <- cbind( TAB.2[,c("IID","DAS","WK","DRUG")], predict( JJ.4, newdata=TAB.2 ) )

plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9) )
PTS <- lapply( Samps.PBO[1:20], function(x)points( DAS ~ WK, data=TAB.2,subset=IID==x,type="o",pch=16 ,col=adjustcolor(col="black",alpha=.5) ))

PRED_v_REAL <- function( sample ) {
	COLS <- c("black","chartreuse2")
	x <- sample
	plot( DAS ~ WK, data=PRED,subset=IID==x,type="p",pch=c(16,1)[factor(PLAC)] ,col=adjustcolor(col=COLS,alpha=.5)[factor(DRUG)],xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main=x )
	points( Estimate ~ WK, data=PRED,subset=IID==x,type="l" )
	text( 0,2, paste(round(RANEF[x,1],2),collapse="\n"), pos=4,col=COLS[1] )
	text( 100,8, paste(round(RANEF[x,2],2),collapse="\n"), pos=2,col=COLS[2] )
}
par(mfrow=c(4,5))
par(mar=c(5,4,1,1))
PTS <- lapply( Samps.PBO[1:20], PRED_v_REAL )




TAB.null <- data.frame( IID=rep("A",100),DAS=rep(5,100),WK=1:100,DRUG=rep(0:1,c(30,70)),PLAC=rep(c(0,1,0),c(1,29,70)) )
TAB.null <- Reduce( rbind, lapply(unique(TAB.j$IID),function(x)data.frame( IID=rep(x,101),DAS=rep(NA,101),WK=0:100,DRUG=rep(0:1,c(30,71)),PLAC=rep(c(0,1,0),c(1,29,71)) ) ))
PRED.100 <- predict(JJ.6,newdata=TAB.null)
PRED <- cbind( TAB.null, PRED.100 )



TEMP <- posterior_samples(JJ.p)
RAND.drug.which <- grep( ",DRUG]$",colnames(TEMP),value=T )
RAND.drug <- TEMP[,"b_DRUG"] + TEMP[,RAND.drug.which]
RAND.drug.order <- order( apply(RAND.drug,2,median) )
par(mar=c(6,3,3,1))
boxplot( RAND.drug[,RAND.drug.order],pch=16,xaxt="n" )
abline(h=seq(-10,10,1),lty=3,col="grey50",lwd=1 )
abline(h=0)
boxplot( RAND.drug[,RAND.drug.order],pch=16,xaxt="n",add=T )
axis( 1, at=1:ncol(RAND.drug),label=gsub("r_IID\\[|,DRUG\\]","",colnames(RAND.drug)[RAND.drug.order]),las=2 )
# lapply( 1:ncol(RAND.drug),function(x)points(rep(x,nrow(RAND.drug)),RAND.drug[,x],col=adjustcolor("black",alpha=.2)) )


head(predict(JJ.6,newdata=TAB.null,re_formula=NA,allow_new_levels=T))
###############################################
## bsts: Bayesian Structural Time Series ######
###############################################




