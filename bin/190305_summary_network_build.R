rm(list=ls())
setwd("/Users/yanjunzan/Library/Mobile Documents/com~apple~CloudDocs/Projects/yeast/")
source(file = "~/Dropbox/home_bin/yeast/Functions.yeast.R")
source(file = "./bin/Functions.yeast2.R")
load("./data/yeast.GenABEL.Data")
require(GenABEL)
load("./data/Trait2name.RData")

mrk <- read.table("./data/bloom2015_int.txt",header = T,sep = "\t",stringsAsFactors = F)
phe.pick <- unique(mrk$trait)
map.all <- map(data)
snp.all <- snpnames(data)


######################################################################################################
## read in the data from Bloom et al 
######################################################################################################

require(gdata)
bloom.int <- read.table(file = "./data/bloom2015_int.txt", sep = "\t", header = T,stringsAsFactors = F) #The pairwise interactions mapped in Bloom2015

## update SNP index to SNP name
epiSNP1 <- snp.all[bloom.int$Q1_index]
epiSNP2 <- snp.all[bloom.int$Q2_index]
addSNP1 <-  snp.all[bloom.int$closest.additive.Q1]
addSNP2 <- snp.all[bloom.int$closest.additive.Q2]
bloom_add_SNP <- data.frame(bloom.int,epiSNP1,epiSNP2,addSNP1,addSNP2,stringsAsFactors = F) # add SNP names


######################################################################################################
## Unifing the additive QTL names in this supplements with the joint set detected in our study
## Keep the epistatic QTL name as the same in Forsberg 2017 Nat Gen
######################################################################################################
load("results/181031_all_additive_qtl.RData") # # this is all the additive QTL saved in four list
qtl_select <- unique(unlist(qtl)) # all QTL
all_add <- unique(c(bloom_add_SNP$addSNP1,bloom_add_SNP$addSNP2)) # all the additive QTL in Bloom et al
all_add_diff <- setdiff(all_add,qtl_select) # the additive QTL in Bloom et al that are not overlapping with ours, due to LD
counter_diff <- numeric(length(all_add_diff))
names(counter_diff) <- all_add_diff

############# Here the distance and r2 will be checked 
out_summary <- data.frame(array(NA,c(length(all_add_diff),4)))
for( i in 1:length(all_add_diff)){
  snp_now <- all_add_diff[i]
  # find best hit
  counter_diff[i] <-  find_best(snp_now = snp_now,genable = data,snp_list = qtl_select,r2_cut = 0.9)
  out_summary[i,1] <- snp_now
  out_summary[i,2] <- counter_diff[i]
  out_summary[i,3] <- abs(map(data[,snp_now]) - map(data[,counter_diff[i]]))
  out_summary[i,4] <- r2fast(data = data,snpsubset = c(counter_diff[i],snp_now))[1,2] 
}
sum(out_summary$X3 <2e4)
sum(out_summary$X3 <2e4 | out_summary$X4 > 0.8 )/nrow(out_summary)
#all these additive QTL(within 20kb and r2 > 0.9,) in this supplements of Bloom et al will be replaced 

for( i in 1:nrow(bloom_add_SNP)){
  check1 <- bloom_add_SNP$addSNP1[i] %in% qtl_select
  check2 <- bloom_add_SNP$addSNP2[i] %in% qtl_select
  if(!check1)
    bloom_add_SNP$addSNP1[i] <- counter_diff[bloom_add_SNP$addSNP1[i]]
  if(!check2)
    bloom_add_SNP$addSNP2[i] <- counter_diff[bloom_add_SNP$addSNP2[i]]
}

bloom_add_epi_SNP_uni <-  bloom_add_SNP
trait_selected <- unique(bloom_add_epi_SNP_uni$trait)
dim(bloom_add_epi_SNP_uni)
######################################################################################################
## If there are QTLs detected for two or more media located within 50 index(defined in Forsberg 2017 Nat Ge)
## Keep the epistatic QTL name as the same in Forsberg 2017 Nat Gen
######################################################################################################
index <- unique(c(bloom_add_epi_SNP_uni$Q1_index,bloom_add_epi_SNP_uni$Q2_index))
test <- find_neighbours(SNPs = unique(c(bloom_add_epi_SNP_uni$epiSNP1,bloom_add_epi_SNP_uni$epiSNP2)),index = index,genable = data,n=50)

#test <- remove_tight.loci(SNPs = unique(c(bloom_add_epi_SNP_uni$epiSNP1,bloom_add_epi_SNP_uni$epiSNP2)),genable = data,cut_dis = 200e3,cut_r2 = 0.6)
for( i in 1:length(test$re)){
  Index1 <- bloom_add_epi_SNP_uni$epiSNP1 %in% test$re[[i]]
  if(any(Index1))
    bloom_add_epi_SNP_uni$epiSNP1[Index1] <- names(test$re[i])
  Index2 <- bloom_add_epi_SNP_uni$epiSNP2 %in% test$re[[i]]
  if(any(Index2))
    bloom_add_epi_SNP_uni$epiSNP2[Index2] <- names(test$re[i])
}
######################################################################################################
## Testing if there are epistatic interaction detected in other media, also significant in current media
##  under a relexed BF threshold.
######################################################################################################

###################### load residuals from Bloom et al 2015 Nat Comm
#save(nbr,file="./data/181128_nbr_josh_all_traits_residuals.RData")
load("./data/181128_nbr_josh_all_traits_residuals.RData")
#p <- numeric(nrow(bloom_add_epi_epi_SNP_uni))
#nbr=matrix(NA,4390,20)
bloom_add_epi_epi_SNP_uni <- bloom_add_epi_SNP_uni

p <- matrix(data = NA,nrow = nrow(bloom_add_epi_epi_SNP_uni),ncol = length(trait_selected))
colnames(p) <- trait_selected
for ( i in 1:nrow(bloom_add_epi_epi_SNP_uni)){
  for(ii in 1:length(trait_selected)){
    mrks <- c(bloom_add_epi_epi_SNP_uni$epiSNP1[i],bloom_add_epi_epi_SNP_uni$epiSNP2[i])
    T_now <-  trait_selected[ii]
    p[i,ii] <-  Test_pairwise_epi_josh(mrks = mrks,data = data,trait = nbr[,names(get.Gen.id(t.name = T_now))])
    cat(i,"_",ii,"\n")
  }
}

sum(p < 0.05/(nrow(bloom_add_epi_epi_SNP_uni)*ncol(bloom_add_epi_epi_SNP_uni)))
bloom_add_epi_epi_SNP_uni_update <- bloom_add_epi_epi_SNP_uni
## Get list to append
for( i in 1:length(trait_selected)){
  index <- bloom_add_epi_epi_SNP_uni$trait != trait_selected[i]
  index_p <- p[index,trait_selected[i]] < 0.05/(nrow(bloom_add_epi_epi_SNP_uni)*3)
  if(any(index_p)){
    update_now <- bloom_add_epi_epi_SNP_uni[index,][index_p,]
    update_now$trait <- paste0(trait_selected[i],"-updated")
    
    ## compare if this interaction is already there
    all <- paste0(bloom_add_epi_epi_SNP_uni$epiSNP1,bloom_add_epi_epi_SNP_uni$epiSNP2)[bloom_add_epi_epi_SNP_uni$trait %in% trait_selected[i] ]
    now <- paste0(update_now$epiSNP1,update_now$epiSNP2)
    if(any(!(now %in% all)))
      bloom_add_epi_epi_SNP_uni_update <- rbind.data.frame(bloom_add_epi_epi_SNP_uni_update,update_now[!(now %in% all),])
  }
  cat(i,nrow(update_now),"\n")
}

#save(bloom_add_epi_epi_SNP_uni_update,file = "./results/181205_update_all_net-bloom_add_epi_epi_SNP_uni_update.RData")
#load(file = "./results/181205_update_all_net-bloom_add_epi_epi_SNP_uni_update.RData")

bloom_add_epi_epi_SNP_uni <- bloom_add_epi_epi_SNP_uni_update

## build networks
networks <- as.list(1:length(trait_selected))
names(networks) <- trait_selected
require(igraph)
for(trait in trait_selected){
  tmp <- bloom_add_epi_epi_SNP_uni[grep(pattern = trait,x = bloom_add_epi_epi_SNP_uni$trait), ]
  # filter out the duplicated rows
  dup <- !duplicated(paste(tmp$epiSNP1,"-",tmp$epiSNP2))
  tmp2 <- tmp[dup,]
  networks[[trait]] <- graph_from_data_frame(data.frame(tmp2$epiSNP1, tmp2$epiSNP2, LOD = tmp2$LOD,trait=tmp2$trait), directed = F)
}
epi <- list()
add <- list()
for( i in 1:length(trait_selected)){
  epi[[i]] <- get_loci_add_epi(tmp = bloom_add_epi_epi_SNP_uni,trait =trait_selected[i])$epi
  add[[i]] <- get_loci_add_epi(tmp = bloom_add_epi_epi_SNP_uni,trait =trait_selected[i])$add
}
names(epi) <- trait_selected
names(add) <- trait_selected
num_epi <- sapply(epi, length)
#num_epi <- sapply(epi, length)

save(V,num_epi,file="./results/181205-V-num_epi.RData")

############################################ pariwise overlap
# #
# load("./results/181205-trait_ordered.RData")
# inter_mat <- data.frame(array(NA,c(20,20)))
# dimnames(inter_mat) <- list(trait_ordered,trait_ordered)
# 
# for( i in 1:20){
#   t1 <- rownames(inter_mat)[i]
#   for( j in 1:20){
#     t2 <- colnames(inter_mat)[j]
#     inter_mat[i,j] <- length(intersect(epi[[t2]],epi[[t1]]))
#   }
# }
# plate <- c(colorRampPalette(c("skyblue","black"))(100),colorRampPalette(c("black","red"))(100))
# tl.col <- "black"#ifelse(trait_ordered %in% c("Manganese Sulfate","Raffinose","Lactate","Indol acetic acid"),"tomato","black")
# par(mar=c(2,15,4.7,6),xpd=T)
# require(corrplot)
# corrplot(data.matrix(inter_mat[trait_ordered,trait_ordered]),type="lower",tl.col=tl.col,method = "number",col=plate,tl.srt=45,is.corr = FALSE,number.cex = 0.6,tl.cex=0.8)
# #require(ggplot2)
# heatmap(data.matrix(inter_mat))


Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
                                           LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)
test_name <- rep(1:length(names(V(Netall))))
names(test_name) <- names(V(Netall))

hubs <- list()
par(mfrow=c(4,5),mar=c(1,1,1,1))
for(i in 1:length(networks)){
  network <- networks[[i]]
  trait_now <- names(networks)[i]
  trait <- names(networks)[i]
  set1 <- edge.attributes(network)$trait %in% trait_now
  set2 <- rep(F,length(set1))
  set2[grep(pattern = paste0(trait_now,"-updated"),x =edge.attributes(network)$trait)] <- T
  if(any(set2)){
    E(network)$color <- ifelse(set1,"tomato",ifelse(set2,"purple","black"))
    E(network)$width <-  ifelse(set1,1,ifelse(set2,1,0.5))
    
  }else{
    E(network)$color <- ifelse(edge.attributes(network)$trait %in% trait_now ,"tomato","black")
    E(network)$width <-  ifelse(set1,1,ifelse(set2,1,0.5))
  }
  network_nrNeighbors <- sapply(X = V(network), FUN = function(x){length(neighbors(graph = network, v = x))})
  V(network)$label <- test_name[names(V(network))]
  #V(network)$label <- ""
  #   V(network)$label <- network_nrNeighbors
  #   V(network)$label.cex <- 1.6
  E(network)$width <- 1
  if(any(network_nrNeighbors > 4))
    hubs[[trait_now]] <- names(V(network))[network_nrNeighbors > 4]
  V(network)$color <- "lightblue"
  V(network)$color[names(V(network)) %in% unlist(hubs)] <- "lightpink2"#V(network)$color[names(V(network)) %in% hubs[[trait_now]]] <- "red"
  #pdf(file = paste("../figures/fig1/", trait, ".pdf", sep = ""), width = 9, height = 6)
  plot(network,main=trait)
  #dev.off()
}


## build a big, whole network
Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
                                           LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)
#dup <- !duplicated(paste(bloom_add_epi_epi_SNP_uni$epiSNP1,"-",bloom_add_epi_epi_SNP_uni$epiSNP2))
#tmp2 <- bloom_add_epi_epi_SNP_uni[dup,]
#hub_all <- list()
dev.off()
# order the trait
#trait_ordered
par(mfrow=c(4,5),mar=c(1,1,1,1))

#trait_ordered <-
for(i in 1:length(trait_ordered)){
  trait_now <- trait_ordered[i]
  ## Need to replace the duplicated
  Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
                                             LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)
  network_nrNeighbors_all <- sapply(X = V(Netall), FUN = function(x){length(neighbors(graph = Netall, v = x))})
  #hub_all[[trait_now]] <- names(V(Netall))[network_nrNeighbors_all > 4]
  if(i==1){
    a <- layout.auto(Netall)
  }
  # 
  #save(a,file = "./results/181203_layout.RData")
  
  V(Netall)$label <- ""
  V(Netall)$size <- 5
  V(Netall)$frame.width <- 0.05
  ## 
  #E(Netall)$color <- ifelse(edge.attributes(Netall)$trait %in% trait_now ,"tomato","grey")
  # set1
  set1 <- edge.attributes(Netall)$trait %in% trait_now
  set2 <- grepl(pattern = paste0(trait_now,"-updated"),x =edge.attributes(Netall)$trait)
  if(any(set2)){
    E(Netall)$color <- ifelse(set1,"tomato",ifelse(set2,"purple","grey"))
    E(Netall)$width <-  ifelse(set1,1.5,ifelse(set2,1.5,0.5))
    
  }else{
    
    E(Netall)$color <- ifelse(edge.attributes(Netall)$trait %in% trait_now ,"tomato","grey")
    E(Netall)$width <-  ifelse(set1,1.5,ifelse(set2,1.5,0.5))
  }
  if( any(names(hubs) %in% trait_now)){
    hub_now <- hubs[[trait_now]]#hub_all_20$hub_snp[hub_all_20$trait %in% trait_now]
    epi_now<- epi[[trait_now]][!(epi[[trait_now]] %in% hub_now)]
    V(Netall)$color <- ifelse(names(V(Netall)) %in% hub_now,"red",ifelse(names(V(Netall)) %in% epi_now,"blue","grey"))
    #V(Netall)$frame.color <- ifelse(names(V(Netall)) %in% hub_now,"red",ifelse(names(V(Netall)) %in% epi_now,"black","grey"))
  }else{
    
    V(Netall)$color <- ifelse(names(V(Netall)) %in% epi[[trait_now]],"blue","grey")
    #V(Netall)$frame.color <- ifelse(names(V(Netall)) %in% epi[[trait_now]],"grey","grey")
  }
  V(Netall)$frame.color <- "grey"
  V(Netall)$shape <- "circle"#ifelse(names(V(Netall)) %in% unique(unlist(qtl)),"circle", "square")
  plot(Netall,layout=a,main=trait_now)
  #mtext(lab[i],side=3,line=1.1, at=-1.4,cex=1.2)
}
save(networks,hubs,Netall,add,epi,V,qtl,qtl_select,bloom_add_epi_epi_SNP_uni,trait_ordered,file = "./results/190306_intermedidate_files.RData")
# 
# tkplot(Netall,layout=layout_with_dh(Netall))
# library("animation")
# library("igraph")
# 
# ani.options("convert") # Check that the package knows where to find ImageMagick
# ani.options(convert="/opt/ImageMagick/bin/convert") 
# 
# saveGIF(for(i in 1:length(trait_ordered)){
#   trait_now <- trait_ordered[i]
#   ## Need to replace the duplicated
#   Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
#                                              LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)
#   network_nrNeighbors_all <- sapply(X = V(Netall), FUN = function(x){length(neighbors(graph = Netall, v = x))})
#   #hub_all[[trait_now]] <- names(V(Netall))[network_nrNeighbors_all > 4]
#   if(i==1){
#     a <- layout.auto(Netall)
#   }
#   # 
#   #save(a,file = "./results/181203_layout.RData")
#   
#   V(Netall)$label <- ""
#   V(Netall)$size <- 5
#   V(Netall)$frame.width <- 0.05
#   ## 
#   #E(Netall)$color <- ifelse(edge.attributes(Netall)$trait %in% trait_now ,"tomato","grey")
#   # set1
#   set1 <- edge.attributes(Netall)$trait %in% trait_now
#   set2 <- grepl(pattern = paste0(trait_now,"-updated"),x =edge.attributes(Netall)$trait)
#   if(any(set2)){
#     E(Netall)$color <- ifelse(set1,"tomato",ifelse(set2,"purple","grey"))
#     E(Netall)$width <-  ifelse(set1,1.5,ifelse(set2,1.5,0.5))
#     
#   }else{
#     
#     E(Netall)$color <- ifelse(edge.attributes(Netall)$trait %in% trait_now ,"tomato","grey")
#     E(Netall)$width <-  ifelse(set1,1.5,ifelse(set2,1.5,0.5))
#   }
#   if( any(names(hubs) %in% trait_now)){
#     hub_now <- hubs[[trait_now]]#hub_all_20$hub_snp[hub_all_20$trait %in% trait_now]
#     epi_now<- epi[[trait_now]][!(epi[[trait_now]] %in% hub_now)]
#     V(Netall)$color <- ifelse(names(V(Netall)) %in% hub_now,"red",ifelse(names(V(Netall)) %in% epi_now,"blue","grey"))
#     #V(Netall)$frame.color <- ifelse(names(V(Netall)) %in% hub_now,"red",ifelse(names(V(Netall)) %in% epi_now,"black","grey"))
#   }else{
#     
#     V(Netall)$color <- ifelse(names(V(Netall)) %in% epi[[trait_now]],"blue","grey")
#     #V(Netall)$frame.color <- ifelse(names(V(Netall)) %in% epi[[trait_now]],"grey","grey")
#   }
#   V(Netall)$frame.color <- "grey"
#   V(Netall)$shape <- "circle"#ifelse(names(V(Netall)) %in% unique(unlist(qtl)),"circle", "square")
#   plot(Netall,layout=a,main=trait_now)
#   #mtext(lab[i],side=3,line=1.1, at=-1.4,cex=1.2)
# } ,
# interval = 1.2, movie.name="./doc/fig/181104_network_animation1.gif" )
# 
# ############################################################################################

