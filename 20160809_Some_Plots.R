## Some Plots for Lab Meeting
 # Also for Committee Meeting...


SAMPS <- sample( unique(TAB.l$IID), 10 )
COLS <- colorRampPalette(COLS.list)(length(SAMPS))
names(COLS) <- SAMPS
png("~/Dropbox/Schork/JNJ11/Slides/20160803_LabMtg_Longitude/Variable_Response.png",height=1200,width=1600,pointsize=30 )
plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main="DAS over Time" )
abline( h=0:10,lty=3,col="grey50" )
Scrap <- lapply( SAMPS, function(x)points( DAS ~ WK, data=TAB.l,subset=IID==x,type="o",lwd=2,col=COLS[x],pch=c(1,16)[factor(DRUG)] ) )
dev.off()

SAMP <- "B012327-28"
png("~/Dropbox/Schork/JNJ11/Slides/20160803_LabMtg_Longitude/Variable_Response.1.png",height=1200,width=1600,pointsize=30 )
plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main="DAS over Time" )
abline( h=0:10,lty=3,col="grey50" )
Scrap <- lapply( SAMP, function(x)points( DAS ~ WK, data=TAB.l,subset=IID==x,type="o",lwd=2,col=COLS[x],pch=c(1,16)[factor(DRUG)] ) )
dev.off()


SAMP <- "M600319-27"
png("~/Dropbox/Schork/JNJ11/Slides/20160803_LabMtg_Longitude/Variable_Response.2.png",height=1200,width=1600,pointsize=30 )
plot( 0,0,type="n",xlim=c(0,100),ylim=c(1,9),xlab="Week",ylab="DAS",main="DAS over Time" )
abline( h=0:10,lty=3,col="grey50" )
Scrap <- lapply( SAMP, function(x)points( DAS ~ WK, data=TAB.l,subset=IID==x,type="o",lwd=2,col=COLS[x],pch=c(1,16)[factor(DRUG)] ) )
dev.off()



COLS <- c("chartreuse2","gold1","firebrick2")
YLIM <- c(0,250)
png("~/Dropbox/Schork/JNJ11/Slides/20160803_LabMtg_Longitude/Variable_Response.EULAR.png",height=1200,width=1600,pointsize=30 )
layout( matrix(c(1:3,3),ncol=2,byrow=T) )
barplot( table(FUL$EUL_28_BL), col=COLS,ylim=YLIM, main="EULAR Response: 28 weeks",ylab="# Patients" )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
barplot( table(FUL$EUL_52_BL), col=COLS,ylim=YLIM, main="EULAR Response: 52 weeks" )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
YLIM <- c(0,150)
barplot( table(FUL$EUL_28_BL,FUL$EUL_52_BL), beside=T,col=COLS,ylim=YLIM, main="EULAR Response: 28 & 52 weeks",xlab="Week 52 Status",yaxt="n" )
axis( 2, at=seq(0,500,50) )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
legend( "topleft", fill=COLS, legend=unique(FUL$EUL_28_BL), title="Week 28 Status",ncol=3 )
dev.off()


COLS <- rev(c("chartreuse2","firebrick2"))
YLIM <- c(0,400)
png("~/Dropbox/Schork/JNJ11/Slides/20160803_LabMtg_Longitude/Variable_Response.ACR50.png",height=1200,width=1000,pointsize=30 )
layout( matrix(c(1:3,3),ncol=2,byrow=T) )
barplot( table(FUL$ACR50_28wk), col=COLS,ylim=YLIM, main="ACR50 Response: 28 weeks",ylab="# Patients" )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
barplot( table(FUL$ACR50_52wk), col=COLS,ylim=YLIM, main="ACR50 Response: 52 weeks" )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
YLIM <- c(0,350)
barplot( table(FUL$ACR50_28wk,FUL$ACR50_52wk), beside=T,col=COLS,ylim=YLIM, main="ACR50 Response: 28 & 52 weeks",xlab="Week 52 Status",yaxt="n" )
axis( 2, at=seq(0,500,50) )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
legend( "topright", fill=COLS, legend=unique(FUL$ACR50_28wk), title="Week 28 Status",ncol=3 )
dev.off()

COLS <- rev(c("chartreuse2","firebrick2"))
YLIM <- c(0,400)
png("~/Dropbox/Schork/JNJ11/Slides/20160803_LabMtg_Longitude/Variable_Response.ACR20.png",height=1200,width=1000,pointsize=30 )
layout( matrix(c(1:3,3),ncol=2,byrow=T) )
barplot( table(FUL$ACR20_28wk), col=COLS,ylim=YLIM, main="ACR20 Response: 28 weeks",ylab="# Patients" )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
barplot( table(FUL$ACR20_52wk), col=COLS,ylim=YLIM, main="ACR20 Response: 52 weeks" )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
YLIM <- c(0,350)
barplot( table(FUL$ACR20_28wk,FUL$ACR20_52wk), beside=T,col=COLS,ylim=YLIM, main="ACR20 Response: 28 & 52 weeks",xlab="Week 52 Status",yaxt="n" )
axis( 2, at=seq(0,500,50) )
abline( h=seq(0,500,50),lty=3,col="grey50",lwd=1 )
legend( "topright", fill=COLS, legend=unique(FUL$ACR20_28wk), title="Week 28 Status",ncol=3 )
dev.off()








## POPULATION ##

LM.pop <- lm( DAS ~ WK+DRUG+TRT, data=TAB )
new.data <- data.frame( WK=0:100,DRUG=rep(0:1,c(25,76)),TRT=c(0,rep(1,100) ) )
PRED.pop <- predict( LM.pop, newdata=new.data)

png( "/Users/kstandis/Dropbox/UCSD/2016_Committee_Meeting/Figures/Long-1_Population.png", height=1000,width=1600, pointsize=30 )
plot( 0,0,type="n", xlim=c(0,100),ylim=c(0,9), xlab="Week",ylab="DAS", main="DAS over 100-week Trial" )
abline( v=seq(0,100,10),h=0:10, lty=3,col="grey50",lwd=1 )
Scrap <- lapply( SAMPS.short, function(x)points( DAS ~ WK, data=TAB, subset=ID==x, col="grey80",type="l") )
points( 0:100, PRED.pop, type="l", lwd=4, col="black", lty=2 )
dev.off()



## INDIVIDUALS ##
Samps <- sample( FUL$ID[FUL$GRP!="G"], 3 )
# Samps <- sample( FUL$ID[FUL$GRP!="G"], 2 )
LM.ind <- lapply( Samps, function(x)lm( DAS ~ WK+DRUG+TRT, data=TAB, subset=ID==x ) )
new.data <- data.frame( WK=0:100,DRUG=rep(0:1,c(25,76)),TRT=c(0,rep(1,100) ) )
PRED.ind <- lapply( LM.ind, function(x)predict( x ) )

COLS <- c("firebrick2","dodgerblue2","chartreuse2")
png( "/Users/kstandis/Dropbox/UCSD/2016_Committee_Meeting/Figures/Long-2_Individual.png", height=1000,width=1600, pointsize=30 )
plot( 0,0,type="n", xlim=c(0,100),ylim=c(0,9), xlab="Week",ylab="DAS", main="DAS over 100-week Trial" )
abline( v=seq(0,100,10),h=0:10, lty=3,col="grey50",lwd=1 )
Scrap <- lapply( SAMPS.short, function(x)points( DAS ~ WK, data=TAB, subset=ID==x, col="grey80",type="l") )
points( 0:100, PRED.pop, type="l", lwd=4, col="black", lty=2 )
Scrap <- lapply( 1:length(Samps), function(x)points( DAS ~ WK, data=TAB, subset=ID==Samps[x], col=COLS[x], pch=16,cex=1.2 ) )
Scrap <- lapply( 1:length(Samps), function(x)points( TAB$WK[TAB$ID==Samps[x]], PRED.ind[[x]], type="l", lwd=4, col=COLS[x] ) )
dev.off()
