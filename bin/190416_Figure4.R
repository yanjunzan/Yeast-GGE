setwd("~/Documents/ye")
load("./data/yeast.GenABEL.Data")
load("./data/Trait2name.RData")
# load previous results
#save(networks,hubs,Netall,add,epi,V,qtl,qtl_select,bloom_add_epi_epi_SNP_uni,trait_ordered,file = "./results/190306_intermedidate_files.RData")
load(file = "./results/190306_intermedidate_files.RData")
require(igraph)
length(qtl)
length(qtl_select)
#load("./results/190306_intermedidate_files.RData")
hub_snps <- unique(unlist(hubs))
snps_connected <- c()
for ( i in 1:13){
  snp_now <- neighbors(graph = Netall,v = hub_snps[i])
  snps_connected <- c(snps_connected,names(snp_now))
  
}
snps_all <- unique(snps_connected)
length(snps_all)

#range(match(names(ed_rank),ed_all))
E(Netall)$color <- "grey"#ifelse(edge.attributes(Netall)$trait %in% trait_now ,"tomato","grey")
E(Netall)$width <- 1#seq(0.35,1.5,length.out = length(col_pal))[col_ed]
V(Netall)$label <- ""
V(Netall)$size <- 5
V(Netall)$frame.width <- 0.05
V(Netall)$color <- ifelse(names(V(Netall)) %in% hub_snps,"red",ifelse(names(V(Netall)) %in% snps_all,"yellow","lightblue"))
V(Netall)$frame.color <- "grey"
V(Netall)$shape <- "circle"
dev.off()
par(mar=c(4,4,4,0),xpd=T)
load("./results/190314_big_manully.RData")
#layout(matrix(1:2,ncol=2), width = c(10,1),height = c(1,1))
plot(Netall,layout=big_manully)


#####

###################################
###################################
# full set of interactors in the big network in local networks
###################################
###################################


hub_snps <- unique(unlist(hubs))
# connectors across Es
Num_connectors <- data.frame(array(0,dim = c(length(hub_snps),length(networks))),row.names = hub_snps)
colnames(Num_connectors) <- names(networks)

for (i in 1:length(networks)){
  trait_now <- names(networks)[i]
  network_now <- networks[[trait_now]]
  for( ii in 1:length(hub_snps)){
    snp_now <- hub_snps[ii]
    total_now <- length(neighbors(graph = Netall,v = snp_now))
    if(snp_now %in% names(V(network_now))){
      Num_connectors[snp_now,trait_now] <- length(neighbors(graph = network_now,v = snp_now))#/total_now
    }
  }
}




# sort the phentoype in colnames
#save(trait_ordered,file="./results/181205-trait_ordered.RData")
load("./results/181205-trait_ordered.RData")
trait_ordered_update <- trait_ordered[!(trait_ordered %in% "Hydroxyurea")]
Num_connectors <- Num_connectors[,trait_ordered_update]
#source("~/Dropbox/home_bin/yeast/181002_Figure_phenotypic_cor.R")
my_palette <- colorRampPalette(c("lightblue", "blue","black"))(n = 299)
mean_test <- apply(Num_connectors,1,FUN = function(x) return(sum(x>4)))
require(gplots)
## creat dendrogrem
ph <- data@phdata
col.meana <- grep("mean",colnames(ph))
ph.mn <- ph[,col.meana]
Trait_name <- c("Cobalt Chloride","Copper Sulfate","Diamide","E6-Berbamine","Ethanol","Formamide",
                "Hydroxyurea","Indol acetic acid" ,"Lactate","Lactose","Magnesium Chloride",
                "Manganese Sulfate","Menadione","Neomycin","Raffinose","Trehalose","Xylose","YNB",
                "YPD","Zeocin")
# col name 
colnames(ph.mn) <- Trait_name
na <- complete.cases(ph.mn)
mat <- data.matrix(ph.mn[na,])
mat <- mat[,!(colnames(mat) %in% "Hydroxyurea")] 
mat <- scale(mat)
hcr <- hclust(dist(t(mat)))
ddr <- as.dendrogram(hcr)
# dev.off()
# par(mfrow=c(3,1))
# plot(rev(ddr))
# 


# update SNP names to gene names
update_names <- function(input_name){
  row_names <- input_name#rownames(Num_connectors)
  # read in my annotation file
  require(gdata)
  anno <- read.xls("./doc/181218_network.vertice.id.xls")
  gene_snp <- as.character(anno$possible.candidate[match(row_names,anno$v_id)])
  gene_snp <- c("Unknown","HAP1","CUP1","VPS41","KRE33","MKT1","Unknown","MAT-alpha","GPA1","IRA2","Unknown","PMR1","YPT7")
  chr <- gsub(pattern = "\\d+_chr(.*)_(\\d+)_.+_.+",replacement = "\\1",x = row_names)
  pos <- gsub(pattern = "\\d+_chr(.*)_(\\d+)_.+_.+",replacement = "\\2",x = row_names)
  
  output_name <- paste0("Chr",chr,":",pos,":",gene_snp)
}
Num_connectors_plot <- Num_connectors[order(mean_test,decreasing = T),]
rownames(Num_connectors_plot) <- update_names(input_name = rownames(Num_connectors_plot))
######### Making heatmap
heatmap.2(as.matrix(Num_connectors_plot),
          density.info="none", key.title ="",cexRow=1,srtCol=30,
          col = my_palette, Rowv=FALSE, Colv= FALSE,dendrogram="none",
          trace='none',cexCol=1,cex=0.5,lwid = c(1,5),lhei = c(1,4.5),
          margins = c(5.5, 11.5),
          sepwidth=c(0.01,0.01),
          sepcolor="white",
          colsep=1:ncol(Num_connectors),
          rowsep=1:nrow(Num_connectors))


#### Updating the connector to only include the same connection
rownames(Num_connectors) <- hub_snps
num_now <- as.data.frame(Num_connectors[order(mean_test,decreasing = T),])
num_now2 <- num_now
for( i in 1:nrow(num_now2)){
  snp_now <-  rownames(num_now2)[i]
  trait_now <- colnames(num_now2)[which.max(num_now2[i,])]
  radial_all <- names(neighbors(graph = networks[[trait_now]],v = snp_now))
  for( ii in 1:ncol(num_now2)){
    trait_now_col <- colnames(num_now2)[ii]
    if(snp_now %in% names(V(networks[[trait_now_col]]))){
      radial_now <- names(neighbors(graph = networks[[trait_now_col]],v = snp_now))
      if(any(radial_now %in% radial_all)){
        num_now2[i,ii] <- sum(radial_now %in% radial_all)
      }else{
        num_now2[i,ii] <- 0
      }
      
    }else{
      num_now2[i,ii] <- 0
    }
    
  }
  
}
#dev.off()
heatmap.2(as.matrix(num_now2),
          density.info="none", key.title ="",cexRow=1,srtCol=30,
          col = my_palette, Rowv=FALSE, Colv= FALSE,dendrogram="none",
          trace='none',cexCol=1,cex=0.5,lwid = c(1,5),lhei = c(1,4.5),
          margins = c(5.5, 11.5),
          sepwidth=c(0.01,0.01),
          sepcolor="white",
          colsep=1:ncol(Num_connectors),
          rowsep=1:nrow(Num_connectors))


############################################
##### NEW connection formed
############################################

num_new_connection <- as.matrix(Num_connectors_plot) - as.matrix(num_now2)

heatmap.2(as.matrix(num_new_connection),
          density.info="none", key.title ="",cexRow=1,srtCol=30,
          col = my_palette, Rowv=FALSE, Colv= FALSE,dendrogram="none",
          trace='none',cexCol=1,cex=0.5,lwid = c(1,5),lhei = c(1,4.5),
          margins = c(5.5, 11.5),
          sepwidth=c(0.01,0.01),
          sepcolor="white",
          colsep=1:ncol(Num_connectors),
          rowsep=1:nrow(Num_connectors))
####################################################################################################
############ For these that are gone gone, how many form new epistatic intearction
####################################################################################################
####################################################################################################
############################# Loci that are completely gone ########################################

My_if_add <- function(epi_loci,add_loci,gwaadata,cut){
  out <- rep(F,length(epi_loci))
  
  for( k in 1:length(epi_loci)){
    out[k] <- find_best_modify(snp_now2 = epi_loci[k],genable = gwaadata,snp_list = add_loci,r2_cut = cut)
  }
  return(out)
}

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

### Number of SNP that are gone

num_epi_epi <- data.frame(array(0,dim = c(length(hub_snps),length(networks))),row.names = hub_snps)
colnames(num_epi_epi) <- colnames(num_now2)
num_epi_add <- data.frame(array(0,dim = c(length(hub_snps),length(networks))),row.names = hub_snps)
colnames(num_epi_add) <-  colnames(num_now2)
num_addAndepi <- data.frame(array(0,dim = c(length(hub_snps),length(networks))),row.names = hub_snps)
colnames(num_addAndepi) <-  colnames(num_now2)

for( i in 1:nrow(num_now2)){
  snp_now <-  rownames(num_now2)[i]
  trait_now <- colnames(num_now2)[which.max(num_now2[i,])]
  radial_all <- names(neighbors(graph = networks[[trait_now]],v = snp_now))
  epi_epi <- rep(0,ncol(num_now2))
  epi_add <- rep(0,ncol(num_now2))
  epi_addAndepi <- rep(0,ncol(num_now2))
  #epi_gone <- rep(0,ncol(num_now2))
  for( ii in 1:ncol(num_now2)){
    trait_now_col <- colnames(num_now2)[ii]
    ######################
    ######################
    ## This is how much remained as epi
    
    if(snp_now %in% names(V(networks[[trait_now_col]]))){
      radial_now <- names(neighbors(graph = networks[[trait_now_col]],v = snp_now))
      if(any(radial_now %in% radial_all)){
        epi_epi_loci <- radial_now[radial_now %in% radial_all]
        epi_epi[ii] <- sum(radial_now %in% radial_all)
      }else{
        epi_epi_loci <- NULL
        epi_epi[ii] <- 0
      }
      
    }else{
      epi_epi_loci <- NULL
      epi_epi[ii] <- 0
    }
    ######################
    ######################
    # This is ow much remained as add
    add_now <- qtl[[trait_now_col]]
    ifAdd <- My_if_add(epi_loci = radial_all,add_loci = add_now,gwaadata = data,cut = 0.8)
    epi_add[ii] <- length(radial_all[ifAdd])
    ######################
    ######################
    # This is ow much remained as add and epi
    if(!is.null(epi_epi_loci) & epi_add[ii] > 0 & any(epi_epi_loci %in% radial_all[ifAdd])){
      
      epi_addAndepi[ii] <- length(intersect(radial_all[ifAdd],epi_epi_loci))
    }else{
      epi_addAndepi[ii] <-0
    }
    cat(i, "_", ii, "\n")
  }
  ### return values
  num_epi_epi[i,] <- epi_epi
  num_epi_add[i,] <- epi_add
  num_addAndepi[i,] <-epi_addAndepi
}

