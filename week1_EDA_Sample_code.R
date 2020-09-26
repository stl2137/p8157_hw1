#------------------------------------------------------------------------------------------------------------------------------------------------------
# exploring pig growth data
pigs<-read.table("PIGS.txt",sep="\t",header=TRUE)
colnames(pigs)<-1:9
matplot(t(pigs),type="l",xlab="Time in Weeks",ylab="Weight",lwd=1,main="Pigs Weight Growth")
sd1<-apply(pigs,2,sd)
pigs2<-data.frame(pigs)
colnames(pigs2)<-1:9
par(mfcol=c(1,2))
boxplot(pigs2,xlab="Time in Weeks",ylab="Weight",col=5,main="Pigs Weight Growth",cex=0.75)
plot(1:9,sd1,xlab="Time in Weeks",ylab="Std. Dev of Weight",main="Pigs Weight Growth",col=5,pch=19)
abline(lm(sd1~c(1:9)),lty=2)
pigs3<-sweep(pigs,2,apply(pigs,2,mean))
pigs3<-sweep(pigs3,2,sd1,FUN="/")
matplot(t(pigs3),type="l",xlab="Time in Weeks",ylab="Standardized Weight",lwd=1,main="Pigs Weight Growth")


#------------------------------------------------------------------------------------------------------------------------------------------------------
# CD4+ data exploration
#------------------------------------------------------------------------------------------------------------------------------------------------------
library(lattice)
macs<-read.table("MACS.txt",sep="\t",header=TRUE)
# SMoothing kernel
x<-macs$time
y<-macs$cd4
sid<-macs$id

# choose 100 with at least 7 observations
len<-sapply(split(x,sid),length)
id<-which(len>=7)
id<-sample(id,100,replace=FALSE)
pid<-unique(sid)
pid<-pid[id]
pid<-sid%in%pid
x<-x[pid]
y<-y[pid]
sid<-sid[pid]
macs<-macs[pid,]


k1<-ksmooth(x,y,kernel="normal",bandwidth=0.1)
k2<-ksmooth(x,y,kernel="normal",bandwidth=0.5)
k3<-ksmooth(x,y,kernel="normal",bandwidth=1.0)
k4<-ksmooth(x,y,kernel="normal",bandwidth=2.0)
ylims<-c(400,1200)
xlims<-range(x)
plot(x,y,pch=19,cex=0.5,xlim=xlims,ylim=ylims,col="gray",xlab="Years since Seroconversion",ylab="CD4+ Count",main="Kernel Smoothers")
par(new=T)
plot(k1,type="l",col=5,xlim=xlims,ylim=ylims,lwd=3,lty=4,xlab="Years since Seroconversion",ylab="CD4+ Count",main="Kernel Smoothers")
par(new=T)
plot(k2,type="l",col=4,xlim=xlims,ylim=ylims,lwd=3,lty=3,xlab="Years since Seroconversion",ylab="CD4+ Count",main="Kernel Smoothers")
par(new=T)
plot(k3,type="l",col=3,xlim=xlims,ylim=ylims,lwd=3,lty=2,xlab="Years since Seroconversion",ylab="CD4+ Count",main="Kernel Smoothers")
par(new=T)
plot(k4,type="l",col=2,xlim=xlims,ylim=ylims,lwd=3,lty=1,xlab="Years since Seroconversion",ylab="CD4+ Count",main="Kernel Smoothers")
temp <- legend("topright", legend = c(" ", " "," "," "),
               text.width = strwidth("Bandwidth"),
               lty = c(4,3,2,1), col=c(5,4,3,2), lwd=3,xjust = 1, yjust = 1,
               title = "Bandwidth")
text(temp$rect$left + temp$rect$w, temp$text$y,
     c("0.10", "0.50","1.0","2.0"), pos=2)
# comparing smoothers
k1<-ksmooth(x,y,kernel="normal",bandwidth=1.0)
s1<-smooth.spline(x,y)
l1<-loess.smooth(x,y,span=0.5,family="gaussian",evaluation=length(x))
ylims<-c(400,1200)
xlims<-range(x)
plot(x,y,pch=19,cex=0.5,xlim=xlims,ylim=ylims,col="gray",xlab="Years since Seroconversion",ylab="CD4+ Count",main="Comparing Smoothers")
par(new=T)
plot(k1,type="l",col=5,xlim=xlims,ylim=ylims,lwd=3,lty=4,xlab="Years since Seroconversion",ylab="CD4+ Count",main="Comparing Smoothers")
par(new=T)
plot(s1,type="l",col=4,xlim=xlims,ylim=ylims,lwd=3,lty=3,xlab="Years since Seroconversion",ylab="CD4+ Count",main="Comparing Smoothers")
par(new=T)
plot(l1,type="l",col=3,xlim=xlims,ylim=ylims,lwd=3,lty=2,xlab="Years since Seroconversion",ylab="CD4+ Count",main="Comparing Smoothers")
temp <- legend("topright", legend = c(" ", " "," "),
               text.width = strwidth("Smoothers"),
               lty = c(4,3,2), col=c(5,4,3), lwd=3,xjust = 1, yjust = 1,
               title = "Smoothers")
text(temp$rect$left + temp$rect$w, temp$text$y,
     c("Kernel", "Spline","LOWESS"), pos=2)
# Time plot / profile plot
xyplot(cd4~time,groups=id,data=macs,type="l",xlab="Years since Seroconversion",ylab="CD4+ Count")
xyplot(cd4~time,groups=id,data=macs,type="l",xlab="Years since Seroconversion",ylab="CD4+ Count",col="lightgray",
       panel=function(...){
         panel.xyplot(...)
         panel.loess(x,y)
       }
)
len<-length(unique(macs$id))
id2<-rbinom(len,1,p=0.05)
id2[id2==0]<-"lightgray"
id2[id2==1]<-"red"
xyplot(cd4~time,groups=id,data=macs,type="l",xlab="Years since Seroconversion",ylab="CD4+ Count",col=id2)
