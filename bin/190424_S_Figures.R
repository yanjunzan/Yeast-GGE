setwd("~/Documents/ye/")
require(igraph)
load(file = "./results/190306_intermedidate_files.RData")
bloom_add_epi_epi_SNP_uni_dup_index <- paste0(bloom_add_epi_epi_SNP_uni$trait,bloom_add_epi_epi_SNP_uni$epiSNP1,bloom_add_epi_epi_SNP_uni$epiSNP2)
bloom_add_epi_epi_SNP_uni <- bloom_add_epi_epi_SNP_uni[!duplicated(bloom_add_epi_epi_SNP_uni_dup_index),]

Netall <- graph_from_data_frame(data.frame(bloom_add_epi_epi_SNP_uni$epiSNP1, bloom_add_epi_epi_SNP_uni$epiSNP2, 
                                           LOD = bloom_add_epi_epi_SNP_uni$LOD,trait=bloom_add_epi_epi_SNP_uni$trait),directed = F)
bloom_add_epi_epi_SNP_uni[grep("IndolaceticAcid",bloom_add_epi_epi_SNP_uni$trait),]
V_names <- names(V(Netall))
#Netall <- simplify(Netall,remove.multiple = T)
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
    #a <- layout.auto(Netall)
    load(file="./results/190314_big_manully.RData")
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
    E(Netall)$color <- ifelse(set1,"tomato",ifelse(set2,"tomato","grey"))
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
  # if(i ==19){
  #   plot(0,yaxt="n",xaxt="n",frame.plot = F,col="white")
  #   plot(0,yaxt="n",xaxt="n",frame.plot = F,col="white")
  # }
  plot(Netall,layout=big_manully,main=trait_now)
  #mtext(lab[i],side=3,line=1.1, at=-1.4,cex=1.2)
}
##########################################################
##########################################################
##### Figure S3
##########################################################
##########################################################
# #################################
# 
# ## Shared vertice
# ##################################
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
#shared <- loci_all >2

load(file="./results/190314_big_manully.RData")

update_names <- function(input_name){
  row_names <- input_name#rownames(Num_connectors)
  # read in my annotation file
  require(gdata)
  anno <- read.xls("./doc/181218_network.vertice.id.xls")
  gene_snp <- as.character(anno$possible.candidate[match(row_names,anno$v_id)])
  #gene_snp <- c("Unknown","HAP1","CUP1","VPS41","KRE33","MKT1","Unknown","MAT-alpha","GPA1","IRA2","Unknown","PMR1","YPT7")
  chr <- gsub(pattern = "\\d+_chr(.*)_(\\d+)_.+_.+",replacement = "\\1",x = row_names)
  pos <- gsub(pattern = "\\d+_chr(.*)_(\\d+)_.+_.+",replacement = "\\2",x = row_names)
  return(gene_snp)
  #output_name <- gene_snp#paste0("Chr",chr,":",pos,":",gene_snp)
}

#V_names[match(hub_snps,V_names)] <- 
  

V_label <- rep("",length(names(V(Netall))))#ifelse(!(V_names %in% hub_snps),"",update_names(V_names[match(hub_snps,V_names)]))
Genes <- update_names(input_name = V_names[match(hub_snps,V_names)])
Genes[Genes == ""] <- "Unknown"
V_label[match(hub_snps,V_names)] <- Genes


snps_connected <- c()
for ( i in 1:13){
  snp_now <- neighbors(graph = Netall,v = hub_snps[i])
  snps_connected <- c(snps_connected,names(snp_now))
  
}
snps_all <- unique(snps_connected)
length(snps_all)


#range(match(names(ed_rank),ed_all))
E(Netall)$color <- "grey"#ifelse(edge.attributes(Netall)$trait %in% trait_now ,"tomato","grey")
E(Netall)$width <- 1 #seq(0.35,1.5,length.out = length(col_pal))[col_ed]
V_names <- names(V(Netall))
hub_snps <- unique(unlist(hubs))
V(Netall)$label <- V_label
V(Netall)$size <- 5
V(Netall)$frame.width <- 0.05
V(Netall)$color <- ifelse(names(V(Netall)) %in% hub_snps,"red",ifelse(names(V(Netall)) %in% snps_all,"yellow","lightblue"))
V(Netall)$frame.color <- "grey"
V(Netall)$shape <- "circle"
dev.off()
par(mar=c(4,4,4,0),xpd=T)

layout(matrix(1:2,ncol=2), width = c(10,1),height = c(1,1))
plot(Netall,layout=big_manully)

tkplot(Netall,layout=big_manully)
#big_manully <- tk_coords(tkp.id = 31)
#save(big_manully,file="./results/190314_big_manully.RData")
### highlight loci shared by more than 5 


# loci_5 <- names(num_loci_all[num_loci_all > 5])
# loci_10 <- names(num_loci_all[num_loci_all > 10])
# index <- findInterval(num_loci_all,seq(0,16,length.out = 9))
# require(gplots)
# #col_pal2 <- colorpanel(length(unique(index)),"blue","red")[index]
# col_pal2 <- brewer.pal(9,"YlOrRd")[1:9][index]
# #plot(1:8,col=col_pal2)
# #display.brewer.pal(name = "YlOrRd",n=9)
# #################################
# 
# ## Shared edge
# ##################################
# ed_all <- attr(E(Netall),"vnames")
# ed_rank <- sort(table(attr(E(Netall),"vnames")))
# length(match(ed_all,names(ed_rank)))
# col_ed <- match(ed_all,names(ed_rank))
# 
# #col_ed2 <- findInterval(col_ed,c(unname(quantile(col_ed))))
# require(plot3D)
# col_pal <- ramp.col(col=c("grey","lightblue"),n=length(col_ed),alpha = 1)[col_ed]
# #col_pal2 <- ramp.col(col=c("grey","tomato"),n=5,alpha = 1)[col_ed2]

