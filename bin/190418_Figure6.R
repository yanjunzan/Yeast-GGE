save(trait_ordered,file = "./results/trait_ordered.RData")
load("./results/trait_ordered.RData")
## plot the VA VG V
##########################################################################
rh <- read.table("doc/repeatibility_heritibility.txt",header=T,stringsAsFactors = F)
Vp <- c()
for(i in 1:nrow(rh)){
  Vp <- c(Vp,var(phdata(data)[,get.Gen.id(t.name = rh$trait[i])],na.rm = T))
}

V_G <- rh$repeatibility*Vp
V_A <- rh$heritibility*Vp
V <- cbind.data.frame(Vp,V_G,V_A,rh$repeatibility,rh$heritibility)#[c(1,2,4,3),]
colnames(V) <- c("Vp","Vg","Va","r","h2")
rownames(V) <- rh$trait

# Run 181002_Figure
#source("~/Dropbox/home_bin/yeast/181002_Figure_phenotypic_cor.R")
V <- V[trait_ordered,]

#################### Set YPD an dYNB to 0
V[19,1:3] <- V[19,1:3]/5.7566746
V[20,1:3] <- V[20,1:3]/1.1428124
# pal <- brewer.pal(7, "Set3")[c(1,3,4,5,7)]
# dev.off()
# plot(0,col="white",ylim = c(0,1),xlim = c(0,25),xaxt="n",frame.plot = F,main="",xlab = "")
# #lines(x=c(1:4),y=H2.vc,col=pal3[1])
# lines(x=c(1:20),y=V$Vp,col=pal[1],lwd=2)
# lines(x=c(1:20),y=V$Vg,col=pal[2],lwd=2)#,lty="dashed")
# lines(x=c(1:20),y=V$Va,col=pal[3],lwd=2)#,lty="twodash")
# 
# lines(x=c(1:20),y=V$r,col=pal[4],lwd=2,lty="dotted")#,lty="twodash")
# lines(x=c(1:20),y=V$h2,col=pal[5],lwd=2,lty="dotted")#,lty="twodash")
# #points(x = 1:20,out,col=pal3[1],pch=19)
# leg <- c(expression('V'[P]),expression('V'[G]),expression('V'[A]),expression(italic(r)^2),expression(italic(h)^2))
# legend("topright",legend = c(leg),pch=19,col = pal,border = "white",bty = "n",cex=1.5)
# par(xpd=T)
# text(x = 1:20,y = -0.1,labels = rownames(V),srt=90)
#barplot(rbind(dat$Number[1:12],c(dat$Number[13:16],rep(NA,8))),
#        beside = TRUE,names.arg = y)
dev.off()
par(mfrow=c(5,1),mar=c(0,4,1,4))
#plot(0,col="white",ylim = c(0,20),xlim = c(0,20),xaxt="n",frame.plot = F,main="",xlab = "hub",yaxt="n",)
library("RColorBrewer")
pal <- brewer.pal(7, "Set3")[c(1,3,4,5,6,7)]
barplot(rbind(V$Vp,V$Vg,V$Va),beside = TRUE,col=rev(pal[1:3]),ylim = c(0,1),ylab = "Variance", border=NA)
legend(x = 39,y=0.95,legend = c("Phenotypic variance","Genetic variance","Additive variance"),pch = 19,col=pal[1:3],bg="white",bty="n")

grid(nx = 0,ny=5)
par(mar=c(1,4,1,4))
barplot(rbind(V$r,V$h2),beside = TRUE,col=rev(pal[4:5]),ylim = c(0,1),ylab="Heritability", border=NA)
legend(x = 28,y=0.95,legend = c("Broad sense","Narrow sense"),pch = 19,col=pal[4:5],bg="white",bty="n")

grid(nx = 0,ny=5)
par(xpd=T)
axis(1,at = seq(2,60,by=3),labels = F,pos=-0.02,tck=-0.09)
plot(0,col="white",ylim = c(0,20),xlim = c(0,59),xaxt="n",frame.plot = F,main="",xlab = "",ylab="",yaxt="n")
text(x = seq(1,60,by=3),y = 18,labels = rownames(V),srt=25)
#save(trait_ordered,file="./results/181205-trait_ordered.RData")

################################################################
load("./results/181205-V-num_epi.RData")
num_epi_uni <- num_epi[rownames(V)]
diff_H_h <- V$r -V$h2
dev.off()
# get the nuber of interactors and variance of additive effects
load("./results/181205-nun_interactors_and-var.RData")
#par(mfrow=c(2,1),mar=c(4,4,0,2))

par(mfrow=c(2,1),mar=c(4,5,1,2))
# plot(nun_interactors,var,frame.plot = F,cex=0.8,pch=19,col="blue",xlab="Number of interactors",cex.lab=0.9,ylab="Variance of additive effects")
# sum(nun_interactors > 4)
# hubs_plot <- names(var)[nun_interactors > 4]
# lm <- lm(var~nun_interactors)
# abline(lm,lty="dashed",col="red")
# mtext("A",side=3,line=-1.5, at=-1,cex=1.2)
# plot(num_epi_uni,diff_H_h,pch=19,col="blue",frame.plot = F,xlab = "Number of pairwise epistatic interactions",cex.lab=0.9,ylab = expression(paste("Difference between H and ",italic("h"),""^{2}*"" )))
# te_lm <-lm(diff_H_h~num_epi_uni)
# summary(te_lm)
# abline(te_lm,col="red",lty="dashed")
# mtext("B",side=3,line=-0.5, at=7.2,cex=1.2)
# diff
diff_H <- outer(X = V$r,Y = V$r,FUN = "-")
num_epi_uni2 <- num_epi_uni
num_epi_uni2[2] <- 0
names(num_epi_uni2) <- rownames(V)
diff_epi_num <- outer(X = num_epi_uni2,Y = num_epi_uni2,FUN = "-")
plot(x=diff_epi_num[upper.tri(diff_epi_num)],y=diff_H[upper.tri(diff_H)],pch=19,col="blue",frame.plot = F,cex.lab=0.9,cex=0.7,ylab = "Pairwise difference in H",xlab="Pairwise difference in the number of pairwise epistatic intearctions")
te_lm2 <-lm(diff_H[upper.tri(diff_H)]~diff_epi_num[upper.tri(diff_epi_num)])
a<- summary(te_lm2)
a$coefficients
abline(te_lm2,col="red",lty="dashed")
mtext("A",side=3,line=-1.2, at=-35,cex=1.2)


# diff_H_h_p <- outer(X = diff_H_h -diff_H_h,Y = V$r,FUN = "-")
# num_epi_uni2 <- num_epi_uni
# num_epi_uni2[2] <- 0
# names(num_epi_uni2) <- rownames(V)
# diff_epi_num <- outer(X = num_epi_uni2,Y = num_epi_uni2,FUN = "-")
# plot(x=diff_epi_num[upper.tri(diff_epi_num)],y=diff_H_h_p[upper.tri(diff_H_h_p)],pch=19,col="blue",frame.plot = F,cex.lab=0.9,cex=0.7,ylab = "Pairwise difference in H",xlab="Pairwise difference in the number of pairwise epistatic intearctions")
# te_lm2 <-lm(diff_H[upper.tri(diff_H)]~diff_epi_num[upper.tri(diff_epi_num)])
# a<- summary(te_lm2)
# a$coefficients
# abline(te_lm2,col="red",lty="dashed")
# mtext("A",side=3,line=-1.2, at=-35,cex=1.2)


