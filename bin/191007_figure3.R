load("./results/figure3.RData")
#dev.off()
par(xpd=T,mar=c(4,4,2,1),mfrow=c(4,4),cex.lab=1.52,cex.axis=1.1)
plot(p$loadings[,1], p$loadings[,2],col="black",frame.plot = F,xlab="PC1",ylab = "PC2",pch= ifelse(rownames(p$loadings) %in% index_selected,19,0) )
text(x =p$loadings[,1],y=p$loadings[,2]+0.04,labels =rownames(p$loadings),cex=cex_now)
#legend(x=-0.1,y=0.1,legend=c("Carbon sources","Oxidative stress","Unclear","Ca2+ Signaling releated"),fill = c("black","blue","purple","tomato"),cex=.6,bty = "n")
mtext(text = "A",side = 3,line = 0,at = -0.5,cex=2)
hist(as.numeric(num_add_3),col="lightblue",xpd=T,xlab = "Number of environments \ndetected as additive loci",main="",ylab = "Count")
mtext(text = "B",side = 3,line = -1,at = 0,cex=2)
hist(num_epi_3,col="lightblue",xlab = "Number of environments \ndetected as epistatic loci",main="",ylab = "Count",xpd=T)
mtext(text = "C",side = 3,line = -1,at = 0,cex=2)

hist(as.numeric(ed_rank_3),col="lightblue",xpd=T,xlab = "Number of environments detected for\n a particular pairwise interaction",main="",ylab = "Count")
mtext(text = "D",side = 3,line = -1,at = 0,cex=2)

dev.off()
par(mfrow=c(3,3),mar=c(2,2,2,2))
#lab <- c("D","E","F")
lab <- c("D","E","F")
require(RColorBrewer)
col <- brewer.pal(9,"Reds")[c(2,4,9)] #colorpanel(100,"lightblue","darkblue")

lab_main <- c("IAA","Formamide","Raffinose")
for ( i in 1:3){
  if(i ==1){
    #a <- layout.fruchterman.reingold(Netall_update,weights=rep(2,39))
    #a <- layout.graphopt(Netall_update,charge =  0.005,spring.length=1)
    #save(a,file="./results/190207_layout.RData")
    load("./results/190207_layout2.RData")
  }
  trait_now <- trait_selected[i]
  add_now <- Add_for_select[[trait_now]]
  epi_now <- epi[[trait_now]]
  hub_now <- hubs[[trait_now]]
  ######################## Shape 
  V(Netall_update)$shape <- ifelse(names(V(Netall_update)) %in% add_now,"circle",ifelse(names(V(Netall_update)) %in% epi_now ,"square","sphere"))
  V(Netall_update)$label <- "" # test_name[names(V(Netall))]
  V(Netall_update)$size <- ifelse(names(V(Netall_update)) %in% "4928552_chrVIII_98622_C_G",12,10) #"4928552_chrVIII_98622_C_G" "9680784_chrXIV_433148_G_A" 
  #   V(network)$label <- network_nrNeighbors
  #   V(network)$label.cex <- 1.6
  set1 <- edge.attributes(Netall_update)$trait %in% trait_now
  set2 <- rep(F,length(set1))
  set2[grep(pattern = paste0(trait_selected[i],"-updated"),x =edge.attributes(Netall_update)$trait)] <- T
  set3 <- edge.attributes(Netall_update)$index_IAA & edge.attributes(Netall_update)$trait %in% trait_now 
  
  if(any(set2)){
    E(Netall_update)$color <- ifelse(set1,"tomato",ifelse(set2,"tomato",rgb(0,0,0,alpha = 0.5)))
    ###E(Netall_update)$color <- ifelse(set2 & shared_edgeall,"red",ifelse(set2 & edge2,"darkorange2","darkgreen"))
    E(Netall_update)$width <-  ifelse(set3,1,ifelse(set1,1,ifelse(set2,1,0.25)))
  }else{
    set_new <- edge.attributes(Netall_update)$trait %in% trait_now
    E(Netall_update)$color <- ifelse(edge.attributes(Netall_update)$trait %in% trait_now ,"tomato",rgb(0,0,0,alpha = 0.5))
    ###E(Netall_update)$color <- ifelse(set_new & shared_edgeall,"red",ifelse(set_new & edge2,"darkorange2","darkgreen"))
    
    E(Netall_update)$width <-  ifelse(set3,1,ifelse(set1,1,ifelse(set2,1,0.25)))
  }
  loci_now <- names(V(Netall_update)) %in% epi_now
  epi1_in <- names(V(Netall_update)) %in% epi1
  epi2_in <- names(V(Netall_update)) %in% epi2
  epi3_in <- names(V(Netall_update)) %in% epi3
  
  #V(Netall_update)$color <- ifelse(loci_now & Unique_loci[[i]],col[1],ifelse(loci_now & Shared_loci_all,col[3],ifelse(loci_now & Shared_loci[[i]],col[2],"white")))
  V(Netall_update)$color <- ifelse(loci_now & epi3_in,col[3],ifelse(loci_now & epi2_in,col[2],ifelse(loci_now & epi1_in,col[1],"white")))
  #V(Netall_update)$color <- ifelse(num_epi_3
  #ifelse(names(V(Netall)) %in% inter_all,"lightpink2",ifelse(names(V(Netall)) %in% unique_loci[[i]],"grey","lightblue"))
  #plot(Netall_update)
  plot(Netall_update,layout=b,main=lab_main[i])
  #legend("topleft",cex=2,legend = c("Shared with 3","Shared with 2","Unique"),col=c("lightpink2","lightblue","grey"),pch = 19,bty="n")
  ### plot network evolution
  if(i==1){
    mtext(lab[i],side=3,line=0, at=-1,cex=2)
  }else{
    mtext(lab[i],side=3,line=0, at=-1.4,cex=2)
  }
}
par(mar = c(4.5, 4.5, 3, 2) + .1,xpd=F)
#plot(0,type="n",col="white",yaxt="n",xaxt="n",frame.plot = F,main="",xlab = "",ylab = "")
#layout(matrix(c(1,2:4,1,5:7),ncol =4,byrow = T))
## panel B first plot G-P on Formamide
trait_now <- "IndolaceticAcid"
t <- unname(get.Gen.id(t.name = trait_now))
# box.pheno <- Sort_bxp(data = data,trait_now = t,snps = snps)$bxp
# tmp <- Sort_bxp(data = data,trait_now = t,snps = snps)$tmp
# bxp(box.pheno, ylim = c(-4, 1.9), xaxt = "n", boxfill = c(rep(c3, 32), rep(c1, 32)), frame = F, yaxt = "n", outcex = .3,pch=19,outcol="grey",whiskcol="grey",staplecol="grey")
geno_IAA <- Sort_bxp(data = data,trait_now = t,snps = snps)$geno
pheno_IAA <- Sort_bxp(data = data,trait_now = t,snps = snps)$pheno
na <- Sort_bxp(data = data,trait_now = t,snps = snps)$NNA
num <- apply(geno_IAA[,2:6], 1, FUN = function(x)sum(x=="H"))

bplot <- function(pheno=pheno_IAA,geno_mat = geno_IAA,count = num,text1="Number of growth increasing alleles\n at six-locus IAA network",text2 = "Growth on IAA"){
  group1.count <- num[geno_IAA[,1] == "RM"]
  group2.count <- num[geno_IAA[,1] == "BY"]
  group1 <- geno_IAA[,1] == "RM"
  group2 <- geno_IAA[,1] == "BY"
  ylim = range(pheno_IAA)
  box.group1 <- boxplot(pheno[group1] ~ group1.count, boxwex = .15, at = seq(.9, 5.9),ylim=ylim, col = "blue", xaxt = "n", frame = F, yaxt = "n", outcex = .3,pch=19,outcol="grey",whiskcol="grey",staplecol="grey")
  box.group2 <- boxplot(pheno[group2] ~ group2.count, boxwex = .15, at = seq(1.1, 6.1),ylim=ylim, add = T, col = "tomato", xaxt = "n", frame = F, yaxt = "n", outcex = .3,pch=19,outcol="grey",whiskcol="grey",staplecol="grey")
  axis(side = 1, at = 1:6, labels = 0:5, cex.axis = 1.6, padj = .5, line = .8)
  axis(side = 2, cex.axis = 1.6, line = -2)
  mtext(text1, side = 1.5, cex = 1, line = 4.5)
  mtext(text2, side = 2, cex = 1, line = .5)
  #legend(x = .6, y = -2.5, c("BY hub-QTL allele", "RM hub-QTL allele"), col = c(cols[2], cols[1]), pch = 15, cex = .75, bty = "n")
  #legend(x = .6, y = -2, c("Additive model fit", "Exponential model fit"), col = c("black", "blue"), lty = "solid", lwd = 3, cex = .75, bty = "n")
}


# boxplot(pheno_IAA ~ geno_IAA[,1] + num ,names=c(rep(0:5,each=2)),col=c("tomato","blue"),ylab="Growth on IAA",xlab="Number of growth increasing allele\n at six-locus IAA network",frame=F, outcex = .3,pch=19,outcol="grey",whiskcol="grey",staplecol="grey")
bplot(pheno=pheno_IAA,geno_mat = geno_IAA,count = num,text1="Number of growth increasing alleles\n at six-locus IAA network",text2 = "Growth on IAA" )
legend("bottomright",legend = c("BY allele","RM allele"),fill= c("tomato","blue"),bty="n")
mtext("H",side=3,line=1.1, at=-1,cex=2)

## continue to For and raf

trait_now <- "Formamide"
t <- unname(get.Gen.id(t.name = trait_now))
pheno_For <- phdata(data)[na,t]
#boxplot(pheno_For ~ geno_IAA[,1] + num ,names=c(rep(0:5,each=2)),col=c("tomato","blue"),ylab="Growth on Formamide",xlab="Number of growth increasing allele\n at six-locus IAA network",frame=F, outcex = .3,pch=19,outcol="grey",whiskcol="grey",staplecol="grey")
bplot(pheno=pheno_For,geno_mat = geno_IAA,count = num,text1="Number of growth increasing alleles\n at six-locus IAA network",text2 = "Growth on Formamide" )

mtext("I",side=3,line=1.1, at=-1.4,cex=2)



trait_now <- "Raffinose"
t <- unname(get.Gen.id(t.name = trait_now))
pheno_raf <- phdata(data)[na,t]
#boxplot(pheno_raf ~ geno_IAA[,1] + num ,names=c(rep(0:5,each=2)),col=c("tomato","blue"),ylab="Growth on Raffinose",xlab="Number of growth increasing allele\n at six-locus IAA network",frame=F, outcex = .3,pch=19,outcol="grey",whiskcol="grey",staplecol="grey")
bplot(pheno=pheno_raf,geno_mat = geno_IAA,count = num,text1="Number of growth increasing alleles\n at six-locus IAA network",text2 = "Growth on Raffinose" )

mtext("J",side=3,line=1.1, at=-1.4,cex=2)



col = c("tomato","blue")
plot(phdata(data)[,t1],phdata(data)[,t2],col=col,frame.plot=F,pch=19,cex=0.3,xlab="Growth on IAA ",ylab="Growth on Formamide")
a <- lm(phdata(data)[,t2][g==0] ~phdata(data)[,t1][g==0])
abline(a,col="tomato")
b <- lm(phdata(data)[t2][g==2] ~phdata(data)[,t1][g==2])
abline(b,col="blue")
mtext("K",side=3,line=1.1, at=-3.4,cex=2)
## label correlation 
r1 <- cor(phdata(data)[,t2][g==0] ,phdata(data)[,t1][g==0],use="pairwise.complete")
r2 <- cor(phdata(data)[,t2][g==2] ,phdata(data)[,t1][g==2],use="pairwise.complete")
r1.f <- as.numeric(format(r1,digits = 2))
r2.f <- format(r2,digits = 2)

legend("topleft",legend =c(paste0("Pearson r^2 = ",r1.f),paste0("Pearson r^2 = ",r2.f)) ,fill= c("tomato","blue"),bty="n")

#eval(paste(expression(paste("Pearson r"^"2")),"=",r1.f))
#my_string <- "Pearson r"
#bquote(.(my_string)^2~"big")
#plot(1,1, main=)


plot(phdata(data)[,t1],phdata(data)[,t3],col=col,frame.plot=F,pch=19,cex=0.3,ylim=c(-4,4),xlab=" Growth on IAA ",ylab="Growth on Raffinose")
a <- lm(phdata(data)[,t3][g==0] ~phdata(data)[,t1][g==0])
abline(a,col="tomato")
b <- lm(phdata(data)[,t3][g==2] ~phdata(data)[,t1][g==2])
abline(b,col="blue")
mtext("L",side=3,line=1.1, at=-3.4,cex=2)
r1 <- cor(phdata(data)[,t3][g==0] ,phdata(data)[,t1][g==0],use="pairwise.complete")
r2 <- cor(phdata(data)[,t3][g==2] ,phdata(data)[,t1][g==2],use="pairwise.complete")
r1.f <- round(r1,digits = 2)
r2.f <- round(r2,digits = 2)

legend("topleft",legend =c(paste0("Pearson r^2 = ",r1.f),paste0("Pearson r^2 = ",r2.f)) ,fill= c("tomato","blue"),bty="n")

# plot(phdata(data)[,t2],phdata(data)[,t3],col=col,frame.plot=F,pch=19,cex=0.3,xlab=" Growth on Formamide ",ylab="Growth on Raffinose")
# a <- lm(phdata(data)[,t3][g==0] ~phdata(data)[,t2][g==0])
# abline(a,col="tomato")
# b <- lm(phdata(data)[,t3][g==2] ~phdata(data)[,t2][g==2])
# abline(b,col="blue")

# ## genetic robustness
phe <- apply(phdata(data)[,c(t1,t2,t3)],MARGIN = 2,FUN=scale)
var <- apply(phe,MARGIN = 1,FUN=var)
geno_hub <- as.double.gwaa.data(data[,"4944074_chrVIII_114144_A_G"])
geno_hub[geno_hub[,1]==0,] <- "BY"
geno_hub[geno_hub[,1]==2,] <- "RM"
boxplot(var~geno_hub,cex=0.2,xlab="Hub genotype",ylab="Within strain growth variance",ylim=c(0,6),col=c("tomato","blue"),frame=F,border="black")
arrows(x0 = 1.05,y0 = 5,x1 = 1.95,y1 = 5,col = "black",angle = 90,length = 0.04)
arrows(x1 = 1.05,y1 = 5,x0 = 1.95,y0 = 5,col = "black",angle = 90,length = 0.04)
text(x = 1.5,y=4.5,labels = "P = 1.03e-18")
text(x = 1.5,y=5.5,labels = "***",cex=1.5)
mtext("M",side=3,line=1.1, at=-0,cex=2)


