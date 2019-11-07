setwd("/Users/yanjunzan/Library/Mobile Documents/com~apple~CloudDocs/Projects/yeast/")
source(file = "~/Dropbox/home_bin/yeast/Functions.yeast.R")
source("./bin/Functions.yeast2.R")
load("./data/yeast.GenABEL.Data")
require(GenABEL)
load("./data/Trait2name.RData")
require(RColorBrewer)
#save(networks,hubs,Netall,add,epi,V,qtl,qtl_select,bloom_add_epi_epi_SNP_uni,trait_ordered,file = "./results/190306_intermedidate_files.RData")
load("./results/190306_intermedidate_files.RData")
trait_selected <- c("IndolaceticAcid","Formamide","Raffinose")
require(plot3D)
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
mat <- data.matrix(ph.mn[na,])
mat <- scale(mat)
p <- princomp(mat)
#dev.off()
#scatter3D(p$loadings[,1], p$loadings[,2], p$loadings[,3],bty="g")
cl1 <- c("Neomycin","Manganese Sulfate","E6-Berbamine","Hydroxyurea","Cobalt Chloride")
cl2 <- c("Menadione","Formamide","Zeocin","Diamide","Indol acetic acid")
cl3 <- c("Copper Sulfate","Ethanol","Magnesium Chloride")
cl4 <-c("YPD","YNB","Raffinose","Xylose","Trehalose","Lactose")
library("RColorBrewer")
pal <- c("tomato","blue","purple")
#display.brewer.pal(5,"Set3")
col <- ifelse(rownames(p$loadings) %in% cl1,pal[1],ifelse(rownames(p$loadings) %in% cl2,pal[2],ifelse(rownames(p$loadings) %in% cl3,pal[3],"black")))
#dev.off()
index_selected <- c("Indol acetic acid","Raffinose","Formamide")
cex_now <- ifelse(rownames(p$loadings) %in% index_selected,0.8,0.01)
col2 <- col
col2[!(rownames(p$loadings) %in% index_selected)] <- "white"
#text3D(p$loadings[,1], p$loadings[,2], p$loadings[,3],labels = rownames(p$loadings),cex=0.8,phi=40,bty="g",col=col2)
#legend(x = -0.39,y=0.30,legend=c("Carbon sources","Oxidative stress","Unclear","Ca2+ Signaling releated"),fill = c("black","blue","purple","tomato"),cex=.8,bty = "n")
# par(xpd=T,mar=c(4,4,4,4))
# plot(p$loadings[,1], p$loadings[,2],col="black",frame.plot = F,xlab="Pc1",ylab = "Pc2",pch= ifelse(rownames(p$loadings) %in% index_selected,19,0) )
# text(x =p$loadings[,1],y=p$loadings[,2]+0.04,labels =rownames(p$loadings),cex=cex_now)
# #legend(x=-0.1,y=0.1,legend=c("Carbon sources","Oxidative stress","Unclear","Ca2+ Signaling releated"),fill = c("black","blue","purple","tomato"),cex=.6,bty = "n")
# mtext(text = "A",side = 3,line = 2,at = -0.6)


require(igraph)
bloom_add_epi_epi_SNP_uni_dup_index <- paste0(bloom_add_epi_epi_SNP_uni$trait,bloom_add_epi_epi_SNP_uni$epiSNP1,bloom_add_epi_epi_SNP_uni$epiSNP2)
bloom_add_epi_epi_SNP_uni <- bloom_add_epi_epi_SNP_uni[!duplicated(bloom_add_epi_epi_SNP_uni_dup_index),]

inte_IAA_For <- intersect(epi$IndolaceticAcid,epi$Formamide)
inte_IAA_Raf <- intersect(epi$IndolaceticAcid,epi$Raffinose)
inte_Raf_For <- intersect(epi$Raffinose,epi$Formamide)
inter_all <- intersect(inte_Raf_For,epi$IndolaceticAcid)
# unique to 1 
unique_IAA <- setdiff(epi$IndolaceticAcid,unique(epi$Formamide,epi$Raffinose))
# unique to 2
unique_For <- setdiff(epi$Formamide,unique(epi$IndolaceticAcid,epi$Raffinose)) 
unique_raf <- setdiff(epi$Raffinose,unique(epi$IndolaceticAcid,epi$Formamide)) 
unique_loci <- list(unique_IAA,unique_For,unique_raf)
names(unique_loci) <- trait_selected
#######################
#######################
# color
#######################
#######################
col <- brewer.pal(9,"Reds")[c(2,4,9)] #colorpanel(100,"lightblue","darkblue")

#plot(1:19,col=col)

#######################
#######################
# subset 
#######################
#######################
epi <- epi[trait_selected]
add <- add[trait_selected]
networks <- networks[trait_selected]
hubs <- hubs[trait_selected]
bloom_add_epi_epi_SNP_uni <- subset(bloom_add_epi_epi_SNP_uni,bloom_add_epi_epi_SNP_uni$trait %in% c(trait_selected,paste0(trait_selected,"-updated")))

a1 <- bloom_add_epi_epi_SNP_uni[grep("IndolaceticAcid",bloom_add_epi_epi_SNP_uni$trait),]$epiSNP1
a2 <-  bloom_add_epi_epi_SNP_uni[ grep("IndolaceticAcid",bloom_add_epi_epi_SNP_uni$trait),]$epiSNP2
dup <- duplicated(paste0(a1,a2))

#dev.off()
# Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
#                                            LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)
# 
# test_name <- rep(1:length(names(V(Netall))))
# names(test_name) <- names(V(Netall))
# dev.off()
# #par(mfrow=c(3,4),mar=c(1,1,1,1))
# layout(matrix(c(1,2,3,4,4,4),nrow = 2,byrow = T),heights = c(1,2))
# layout(matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,byrow = T),heights = c(1,1,2))
# hubs_4 <- list()
# par(mar=c(1,1,1,1))
# for(i in 1:length(networks)){
#   network <- networks[[i]]
#   trait_now <- names(networks)[i]
#   trait <- names(networks)[i]
#   network_nrNeighbors <- sapply(X = V(network), FUN = function(x){length(neighbors(graph = network, v = x))})
#   V(network)$label <- test_name[names(V(network))]
#   V(network)$size <- ifelse(names(V(network)) %in% hubs[[trait_now]],30,20)
#   #   V(network)$label <- network_nrNeighbors
#   #   V(network)$label.cex <- 1.6
#   E(network)$width <- 3
#   if(any(network_nrNeighbors >= 4))
#     hubs_4[[trait_now]] <- names(V(network))[network_nrNeighbors >= 4]
#   
#   #loci_other <- unique(unlist(epi[!(names(networks) %in% trait_now)]))
#   V(network)$color <- ifelse(names(V(network)) %in% inter_all,"lightpink2",ifelse(names(V(network)) %in% unique_loci[[i]],"grey","lightblue"))
#   #V(network)$color <- ifelse(names(V(network)) %in% unique_loci[[i]],"black","lightblue")
#   #V(network)$color[names(V(network)) %in% hubs[[trait_now]]] <- "red"
#   #pdf(file = paste("../figures/fig1/", trait, ".pdf", sep = ""), width = 9, height = 6)
#   plot(network,main=trait)
#   #dev.off()
#   #if(i ==1)
#   #legend("topleft",legend = c("Shared with 3","Shared with 2","Unique"),col=c("lightpink2","lightblue","grey"),pch = 19,bty="n")
# }
# #plot(c(0.5,1,1.5),c(1,1,1),xlim = c(0,2),col="black",pch=19,frame.plot = F,xaxt='n',yaxt='n',ann=FALSE)
# ## build another net trait 1 and 2
# #bloom_add_epi_epi_SNP_uni_sub1 <- filter(bloom_add_epi_epi_SNP_uni,trait %in% names(networks)[c(1,2)])
# #a <- duplicated(paste0(bloom_add_epi_epi_SNP_uni$epiSNP1,bloom_add_epi_epi_SNP_uni$epiSNP2))
# #sum(a)

Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
                                           LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)


###########################################
###########################################
### Shared connections
###########################################
###########################################
## Get to all network

########################################### ########################################### ########################################### 

########################################### Level of  and shared loci in a network: kick out the residuals

########################################### ########################################### ########################################### 
index_IAA <- rep(F,nrow(bloom_add_epi_epi_SNP_uni))

for( i in 1:nrow(bloom_add_epi_epi_SNP_uni)){
  V_now <- unique(c(bloom_add_epi_epi_SNP_uni$epiSNP1[i],bloom_add_epi_epi_SNP_uni$epiSNP2[i]))
  if(sum(V_now %in% "4928552_chrVIII_98622_C_G")==1)
    index_IAA[i] <- T
}

Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
                                           LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait,index_IAA=index_IAA),directed = F)
network_nrNeighbors <- sapply(X = V(Netall), FUN = function(x){length(neighbors(graph = Netall, v = x))})

cl <- clusters(Netall)
dg <- decompose.graph(Netall)
Netall_update <- dg[[1]]

###########################################
###########################################
# All loci for epistatic
###########################################
###########################################
all_loci_3 <- names(V(Netall_update))

num_epi_3 <- rep(0,length(all_loci_3))
names(num_epi_3) <- all_loci_3
for( i in 1:length(num_epi_3)){
  loci_now <- all_loci_3[i]
  if( loci_now %in% epi[[1]])
    num_epi_3[i] = num_epi_3[i] +1
  if( loci_now %in% epi[[2]])
    num_epi_3[i] = num_epi_3[i] +1
  if( loci_now %in% epi[[3]])
    num_epi_3[i] = num_epi_3[i] +1
  cat(i,"\n")
}
par(mfrow=c(3,3))
hist(num_epi_3,breaks = seq(0,3))
sum(num_epi_3 ==2)

epi2 <- names(num_epi_3)[num_epi_3 ==2] 
epi1 <- names(num_epi_3)[num_epi_3 ==1] 
epi3 <- names(num_epi_3)[num_epi_3 ==3] 

#sum(num_epi_3)
###########################################
###########################################
### All loci for additive
###########################################
###########################################
find_best_modify <- function(snp_now2,genable,snp_list,r2_cut=0.9){
  #snp_now="9714855_chrXIV_467219_A_G";genable=data;snp_list=qtl_select
  # select the chromsome
  if(snp_now2 %in% snp_list){
    return(T)
  }else{
    chr <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_now2)
    chrs <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_list)
    index  <- c(chr,chrs) == chr
    if(length(unique(c(snp_now2,snp_list)[index]))==1){
      return(F)
    }else{
      r2_all <- r2fast(data = genable,snpsubset = unique(c(snp_now2,snp_list)[index]))
      if(max(r2_all[snp_now2,],na.rm = T) > r2_cut){
        return(T)
      }else{
        #warning(i," No hits with r2 above ",r2_cut,", the max is: ",max(r2_all[snp_now,],na.rm = T))
        return(F)
      }
    }
  }
}
Add_for_select <- list()
My_if_add <- function(epi_loci,add_loci,gwaadata,cut){
  out <- rep(F,length(epi_loci))
  for( k in 1:length(epi_loci)){
    out[k] <- find_best_modify(snp_now2 = epi_loci[k],genable = gwaadata,snp_list = add_loci,r2_cut = cut)
  }
  return(out)
}

loci_epi <- unique(c(bloom_add_epi_epi_SNP_uni$epiSNP1,bloom_add_epi_epi_SNP_uni$epiSNP2))
for(i in 1:length(trait_selected)){
  trait_now <- trait_selected[i]
  add_now <- qtl[[trait_now]]
  tf <- My_if_add(epi_loci = loci_epi,add_loci = add_now,gwaadata = data,cut=0.6)
  Add_for_select[[trait_now]] <- loci_epi[tf]
}

num_add_3 <- rep(0,length(all_loci_3))
names(num_add_3) <- all_loci_3
all_epi <- unlist(Add_for_select)
a <- table(all_epi[all_epi %in% all_loci_3])
num_add_3[names(a)] <- a
hist(num_add_3)
###########################################
###########################################
### Shared edge
###########################################
###########################################
ed_all_3 <- attr(E(Netall_update),"vnames")
ed_rank_3 <- sort(table(attr(E(Netall_update),"vnames")))


#Netall_update <- Netall
#plot(Netall,layout=layout.grid(Netall))
#plot(Netall_update,layout=layout.grid(Netall_update))
## Get all the loci in IAA, network
IAA_V_all <- c("4928552_chrVIII_98622_C_G","2357646_chrIV_997624_A_G", "1242017_chrIII_198615_T_G",
               "2975503_chrV_83548_C_T","9680784_chrXIV_433148_G_A","8733525_chrXIII_410320_T_C",
               "10185286_chrXV_153317_T_C","7455567_chrXII_210539_A_G","9623949_chrXIV_376313_C_T")
#IAA_V_GP <- c("4928552_chrVIII_98622_C_G", "9680784_chrXIV_433148_G_A", "1242017_chrIII_198615_T_G", "2357646_chrIV_997624_A_G", "8733525_chrXIII_410320_T_C", "7455567_chrXII_210539_A_G")
#IAA_V_GP2 <- neighbors(networks[["IndolaceticAcid"]],v = "4928552_chrVIII_98622_C_G")(
IAA_V_GP <- c("2371958_chrIV_1011936_A_T" , "1242017_chrIII_198615_T_G", "9728717_chrXIV_481081_C_T", "8732658_chrXIII_409453_T_C", "7473402_chrXII_228374_G_C")
##################################################### 

###################################

############################################# a intersection of the network
Netinter12 <- intersection(as.undirected(networks[[trait_selected[1]]]),as.undirected(networks[[trait_selected[2]]]),keep.all.vertices=F)
Netinter13 <- intersection(as.undirected(networks[[trait_selected[1]]]),as.undirected(networks[[trait_selected[3]]]))
Netinter23 <- intersection(as.undirected(networks[[trait_selected[3]]]),as.undirected(networks[[trait_selected[2]]]))
Netinterall <- intersection(Netinter12,as.undirected(networks[[trait_selected[3]]]))
shared_edge12 <- E(Netall_update) %in% E(Netinter12)
shared_edge13 <- E(Netall_update) %in% E(Netinter13)
shared_edge23 <- E(Netall_update) %in% E(Netinter23)
shared_edgeall <- E(Netinterall)

edge2 <- (shared_edge12 | shared_edge13 | shared_edge23)
#edge3 <- (shared_edge12 & shared_edge13 & shared_edge23)
########################################### Shared loci
Shared_loci <- list()
Unique_loci <- list()
for( i in 1:3){
  other <- c(1:3)[!(c(1:3 %in% i))]
  Shared_loci[[i]] <- names(V(Netall_update)) %in% epi[[i]][epi[[i]] %in% unique(epi[[other[1]]],epi[[other[2]]])]
  Unique_loci[[i]] <- names(V(Netall_update)) %in% epi[[i]][!(epi[[i]] %in% unique(epi[[other[1]]],epi[[other[2]]]))]

  #Shared_loci[[i]] <- epi[[i]] %in% unique(epi[[other[1]]],epi[[other[2]]])
  cat(length(Shared_loci[[i]])," ","\n")
  cat(length(V(networks[[i]]))," ","\n")
}
names(Shared_loci) <- names(networks)
Shared_loci_all <- names(V(Netall_update)) %in% intersect(epi[[1]],intersect(epi[[2]],epi[[3]]))

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
    mtext(lab[i],side=3,line=0, at=-1,cex=1.2)
  }else{
    mtext(lab[i],side=3,line=0, at=-1.4,cex=1.2)
  }
}
save.image("./results/figure3.RData")
