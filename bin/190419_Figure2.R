setwd("~/Documents/ye/")
require(igraph)
require(RColorBrewer)
load(file = "./results/190306_intermedidate_files.RData")
load(file = "./results/190314_big_manully.RData")

bloom_add_epi_epi_SNP_uni_dup_index <- paste0(bloom_add_epi_epi_SNP_uni$trait,bloom_add_epi_epi_SNP_uni$epiSNP1,bloom_add_epi_epi_SNP_uni$epiSNP2)
bloom_add_epi_epi_SNP_uni <- bloom_add_epi_epi_SNP_uni[!duplicated(bloom_add_epi_epi_SNP_uni_dup_index),]

Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
                                           LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)
bloom_add_epi_epi_SNP_uni[grep("IndolaceticAcid",bloom_add_epi_epi_SNP_uni$trait),]
#################################

## Shared vertice
##################################
loci_all <- names(V(Netall))

num_loci_all <- loci_all
names(num_loci_all) <- loci_all
for(i in 1:length(num_loci_all)){
  num_now <- 0
  
  for( j in 1:length(networks)){
    if(any(names(V(networks[[j]])) %in% num_loci_all[i] ))
      num_now <- num_now +1
  }
  num_loci_all[i] <- num_now
  cat(i, "\n")
}

num_loci_all <- as.numeric(num_loci_all)
names(num_loci_all) <- loci_all
loci_5 <- names(num_loci_all[num_loci_all > 5])
loci_10 <- names(num_loci_all[num_loci_all > 10])
index <- findInterval(num_loci_all,seq(0,16,length.out = 9))
require(gplots)
#col_pal2 <- colorpanel(length(unique(index)),"blue","red")[index]
col_pal2 <- brewer.pal(9,"YlOrRd")[1:9][index]
#plot(1:8,col=col_pal2)
#display.brewer.pal(name = "YlOrRd",n=9)
#################################

## Shared edge
##################################
ed_all <- attr(E(Netall),"vnames")
ed_rank <- sort(table(attr(E(Netall),"vnames")))
length(match(ed_all,names(ed_rank)))
col_ed <- match(ed_all,names(ed_rank))

#col_ed2 <- findInterval(col_ed,c(unname(quantile(col_ed))))
require(plot3D)
col_pal <- ramp.col(col=c("grey","lightblue"),n=length(col_ed),alpha = 1)[col_ed]
#col_pal2 <- ramp.col(col=c("grey","tomato"),n=5,alpha = 1)[col_ed2]

#range(match(names(ed_rank),ed_all))
E(Netall)$color <- col_pal#ifelse(edge.attributes(Netall)$trait %in% trait_now ,"tomato","grey")
E(Netall)$width <- seq(0.35,1.5,length.out = length(col_pal))[col_ed]
V(Netall)$label <- ""
V(Netall)$size <- 5
V(Netall)$frame.width <- 0.05
V(Netall)$color <- col_pal2
V(Netall)$frame.color <- "grey"
V(Netall)$shape <- "circle"
#dev.off()
par(mar=c(4,4,4,0),xpd=T)

layout(matrix(1:2,ncol=2), width = c(10,1),height = c(1,1))
plot(Netall,layout=big_manully)
#plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))
par(mar=c(0,0,0,0),xpd=T)
plot(0, col="white",axes = F,xlab = '', ylab = '', main = '',xlim = c(-0.6,0.9),ylim = c(-1.5,1.3))
#legend_image <- as.raster(matrix(rep(rev(colorpanel(9,"grey","red")),each=100)),col=1)
legend_image <- as.raster(matrix(rep(rev(brewer.pal(9,"YlOrRd")[2:9]),each=100)),col=1)
rasterImage(legend_image, -0.3, -0.5, 0.8,0.8)
text(x = c(0.8,0.8,0.8,0.8,0.8)-1.3,seq(-0.5,0.8,length.out = 5),labels = c(0,4,8,12,16))
#axis(4,at=seq(-0.5,0.8,length.out = 5),labels = seq(0,16,length.out = 5))
#dev.off()
#par(mfrow=c(1,3))
layout(matrix(c(1:2,3,4,5,5),ncol=3), width = c(5,5,10),height = c(1,1))
par(mar=c(4,4,2,1))
hist(num_loci_all,col="lightblue",xlab = "Number of detected \nenviroments (as epistatic loci)",main="",ylab = "Count",xpd=T)
mtext(text = "A",side = 3,line = -1,at = -5)
hist(as.numeric(ed_rank),col="lightblue",xpd=T,xlab = "Number of detected \n enviroments (pairwise interactions)",main="",ylab = "Count")
mtext(text = "B",side = 3,line = -1,at = 0)
#save(Num_add_for_epi,r2,file = "./results/19043_Num_add_for_epi_r2.RData") # 190403_epi_add_overlap
load(file = "./results/19043_Num_add_for_epi_r2.RData")
hist(as.numeric(Num_add_for_epi),col="lightblue",xpd=T,xlab = "Number of detected enviroments\n (as additive loci)",main="",ylab = "Count")
mtext(text = "C",side = 3,line = -1,at = -5)
hist(as.numeric(r2),col="lightblue",xpd=T,xlab = " Addtive variance \nexplained (%)",main="",ylab = "Count")
mtext(text = "D",side = 3,line = -1,at = 0.55)
#save(Netall,num_loci_all,ed_rank,Num_add_for_epi,r2,file = "results/190412_Figure2.RData")

