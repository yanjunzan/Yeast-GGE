#setwd("/Users/yanjunzan/Library/Mobile Documents/com~apple~CloudDocs/Projects/yeast/")
#source(file = "~/Dropbox/home_bin/yeast/Functions.yeast.R")
load("./data/yeast.GenABEL.Data")
ph <- data@phdata
col.meana <- grep("mean",colnames(ph))
ph.mn <- ph[,col.meana]
rownames(ph.mn) <- ph[,1]
load("./data/Trait2name.RData")
Trait_name <- c("Cobalt Chloride","Copper Sulfate","Diamide","E6-Berbamine","Ethanol","Formamide",
                "Hydroxyurea","Indol acetic acid" ,"Lactate","Lactose","Magnesium Chloride",
                "Manganese Sulfate","Menadione","Neomycin","Raffinose","Trehalose","Xylose","YNB",
                "YPD","Zeocin")
# col name 
colnames(ph.mn) <- Trait_name
na <- complete.cases(ph.mn)

# phenotype correlation
# require(corrplot)
# par(mar=c(2.3,2,4,0))
# layout(matrix(c(1,2),nrow=1),widths = c(1,3))
# cor.freq <- cor(data.matrix(ph.mn[na,]),method = c("spearman"))
# d <- hclust(as.dist(1-cor.freq),"complete")
# 
# 
# #plot(as.dendrogram(d),horiz=T,xaxt = "n", yaxt = "n",leaflab="none")
# plot(as.dendrogram(d),horiz=T,xaxt = "n", yaxt = "n",leaflab="none")
# par(mar=c(2,10,4.7,2),xpd=T)
# plate <- c(colorRampPalette(c("skyblue","black"))(100),colorRampPalette(c("black","red"))(100))
# tl.col <- ifelse(Trait_name[rev(d$order)] %in% c("Manganese Sulfate","Raffinose","Lactate","Indol acetic acid"),"tomato","black")
# corrplot(cor.freq[rev(d$order),rev(d$order)],type="lower",tl.col=tl.col,method = "number",col=plate,tl.srt=45,is.corr = FALSE,number.cex = 0.4,tl.cex=0.8)
# Trait_name2 <- c("CobaltChloride","CopperSulfate","Diamide","E6-Berbamine","Ethanol","Formamide",
#                 "Hydroxyurea","IndolaceticAcid" ,"Lactate","Lactose","MagnesiumChloride",
#                 "ManganeseSulfate","Menadione","Neomycin","Raffinose","Trehalose","Xylose","YNB",
#                 "YPD","Zeocin")
# trait_ordered <- Trait_name2[rev(d$order)]
# 



######################### Test differnet way of cluster

# phenotype correlation

#heatmap(mat)
mat <- data.matrix(ph.mn[na,])
mat <- scale(mat)
hcr <- hclust(dist(t(mat)))
ddr <- as.dendrogram(hcr)
#dev.off()
par(mfrow=c(3,1))
plot(rev(ddr))

library(ggfortify)
require(ggfortify)
#autoplot(prcomp(t(mat),scale. = T))
p <- prcomp(t(mat),scale. = T)
autoplot(p, data = t(mat), label = TRUE, label.size = 4)
# 
# pc <- prcomp(mat)
# plot(pc$rotation[,1],pc$rotation[,2])

cor.freq <- cor(data.matrix(ph.mn[na,]),method = c("spearman"))
#d <- hclust(as.dist(1-cor.freq),"complete")
#plot(as.dendrogram(d),horiz=T,xaxt = "n", yaxt = "n",leaflab="none")
require(corrplot)
dev.off()
par(mar=c(2.3,2,4,0))
layout(matrix(c(1,2),nrow=1),widths = c(1,3))
plot(as.dendrogram(ddr),horiz=T,xaxt = "n", yaxt = "n",leaflab="none")
par(mar=c(2,10,4.7,6),xpd=T)
plate <- c(colorRampPalette(c("skyblue","black"))(100),colorRampPalette(c("black","red"))(100))
tl.col <- ifelse(Trait_name[rev(ddr$order)] %in% c("Manganese Sulfate","Raffinose","Lactate","Indol acetic acid"),"tomato","black")
a <- cor.freq[rev(hcr$order),rev(hcr$order)]*100
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
#round_df(x = a,digits = 0)
plate <- c(colorRampPalette(c("skyblue","black"))(100),colorRampPalette(c("black","red"))(100))
tl.col <- ifelse(Trait_name[rev(hcr$order)] %in% c("Manganese Sulfate","Raffinose","Lactate","Indol acetic acid"),"tomato","black")
par(mar=c(2,15,4.7,6),xpd=T)
corrplot(round_df(a,0),type="lower",tl.col=tl.col,method = "number",col=plate,tl.srt=45,is.corr = FALSE,number.cex = 0.6,tl.cex=0.8)

## Testing ggplot
library(reshape2)
library(ggplot2)
get_lower_tri = function(cormat) {
  cormat[upper.tri(cormat)] = NA
  return(cormat)
}

get_upper_tri = function(cormat) {
  cormat[lower.tri(cormat)] = NA
  return(cormat)
}

upper_tri <- get_upper_tri(round_df(a,0))
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-100,100), space = "Lab", 
                       name="Spearman Correlation\n times 100") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


Trait_name2 <- c("CobaltChloride","CopperSulfate","Diamide","E6-Berbamine","Ethanol","Formamide",
                 "Hydroxyurea","IndolaceticAcid" ,"Lactate","Lactose","MagnesiumChloride",
                 "ManganeseSulfate","Menadione","Neomycin","Raffinose","Trehalose","Xylose","YNB",
                 "YPD","Zeocin")
trait_ordered <- Trait_name2[rev(hcr$order)]



################################## Generating the 3D plot
#library("plot3D")
# p <- prcomp(mat,scale. = T,center = T)
# #dev.off()
# par(mfrow=c(2,1))
# scatter3D(p$rotation[,1], p$rotation[,2], p$rotation[,3])
# text3D(p$rotation[,1], p$rotation[,2], p$rotation[,3],labels = rownames(p$rotation),cex=0.8,phi=40,bty="g")

p <- princomp(mat)
#dev.off()
scatter3D(p$loadings[,1], p$loadings[,2], p$loadings[,3],bty="g")
cl1 <- c("Neomycin","Manganese Sulfate","E6-Berbamine","Hydroxyurea","Cobalt Chloride")
cl2 <- c("Menadione","Formamide","Zeocin","Diamide","Indol acetic acid")
cl3 <- c("Copper Sulfate","Ethanol","Magnesium Chloride")
cl4 <-c("YPD","YNB","Raffinose","Xylose","Trehalose","Lactose")
library("RColorBrewer")
pal <- c("tomato","blue","purple")
#display.brewer.pal(5,"Set3")
col <- ifelse(rownames(p$loadings) %in% cl1,pal[1],ifelse(rownames(p$loadings) %in% cl2,pal[2],ifelse(rownames(p$loadings) %in% cl3,pal[3],"black")))
text3D(p$loadings[,1], p$loadings[,2], p$loadings[,3],labels = rownames(p$loadings),cex=0.8,phi=40,bty="g",col=col)
legend(x = -0.39,y=0.30,legend=c("Carbon sources","Oxidative stress","Unclear","Ca2+ Signaling releated"),fill = c("black","blue","purple","tomato"),cex=.8,bty = "n")

Col_tarit <- col
names(Col_tarit) <- rownames(p$loadings)
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
plot(0,col="white",ylim = c(0,20),xlim = c(0,20),xaxt="n",frame.plot = F,main="",xlab = "hub",yaxt="n")
library("RColorBrewer")
pal <- brewer.pal(7, "Set3")[c(1,3,4,5,6,7)]
barplot(rbind(V$Vp,V$Vg,V$Va),beside = TRUE,col=pal[1:3],ylim = c(0,1),ylab = "Variance", border=NA)
legend(x = 39,y=0.95,legend = c("Phenotypic variance","Genetic variance","Additive variance"),pch = 19,col=pal[1:3],bg="white",bty="n")

grid(nx = 0,ny=5)
par(mar=c(1,4,1,4))
barplot(rbind(V$r,V$h2),beside = TRUE,col=pal[4:5],ylim = c(0,1),ylab="Heritability", border=NA)
legend(x = 28,y=0.95,legend = c("Broad sense","Narrow sense"),pch = 19,col=pal[4:5],bg="white",bty="n")

grid(nx = 0,ny=5)
par(xpd=T)
axis(1,at = seq(2,60,by=3),labels = F,pos=-0.02,tck=-0.09)
plot(0,col="white",ylim = c(0,20),xlim = c(0,59),xaxt="n",frame.plot = F,main="",xlab = "",ylab="",yaxt="n")
text(x = seq(1,60,by=3),y = 18,labels = rownames(V),srt=25)
#save(trait_ordered,file="./results/181205-trait_ordered.RData")
