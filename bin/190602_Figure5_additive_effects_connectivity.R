setwd("~/Documents/ye/")
require(igraph)
require(GenABEL)
load(file = "./data/yeast.GenABEL.Data")
load(file = "./results/190306_intermedidate_files.RData")
load("./results/190415_trait.name.RData")
source("~/Dropbox/home_bin/yeast/Functions.yeast.R")
trait_picked <- trait.name[-c(18,19)]#get.Gen.id(t.name = phe.pick[!(phe.pick %in% c("YNB","YPD"))])#c(get.Gen.id(t.name = phe.pick[1]),get.Gen.id(t.name = phe.pick[2]),get.Gen.id(t.name = phe.pick[3]),get.Gen.id(t.name = phe.pick[4]))

qtl_select <- unique(unlist(qtl)) #unique(c(BE.yeast.1_v3$out5, BE.yeast.2_v3$out5, BE.yeast.3_v3$out5, BE.yeast.4_v3$out5))
# out_inter_qtl <- fit.GbyE.allele(trait = trait_picked,genable = data,snp_test =qtl_select )
# sum(out_inter_qtl$GxE=="Yes")/nrow(out_inter_qtl) # 
# write.table(out_inter_qtl,file = "./results/181031_Table_S2_P_value_QTLxE.txt", sep = "\t", row.names = F,col.names = T,quote = F) #The pairwise interactions mapped in Bloom2015
out_inter_qtl <- read.table("./results/181031_Table_S2_P_value_QTLxE.txt",sep="\t",header=T,stringsAsFactors = F)
rownames(out_inter_qtl) <- out_inter_qtl$SNP_name
dim(out_inter_qtl)
table(out_inter_qtl$GxE)
################################################################################################
############### Estimate additive effects
################################################################################################
hub_snps <- unique(unlist(hubs))

#hub_snp1 <- networkNeighbours.hubs$snpName #unlist(hubs)
# remove high ld stuff
rm <- remove_tight.loci(SNPs = hub_snps,cut_dis = 50e3,genable = data)
hub_snp <- unique(rm$snp)
length(hub_snp)

eff_picked <- est.genetic.effect(trait= trait_picked,genable=data,snp_test=qtl_select)
var_eff <- function(input){
  return(var(input[seq(1,length(input),3)]))
}
eff_no_inter <- eff_picked[,2:ncol(eff_picked)]
var <- apply(eff_no_inter,MARGIN = 2,FUN = var_eff)
out <- rbind.data.frame(eff_no_inter,var)
out_sorted <- out[,order(-out[nrow(out),])]
n <- ncol(eff_no_inter)

eff_size_only <- out[seq(1,54,3),order(-out[nrow(out),])]
rownames(eff_size_only) <- gsub(pattern = "1_(.*)_1",replacement = "\\1",names(trait_picked))

#dev.off()
########################################################################################################################
#####################################################################################################################


#### Number of interactors and variance in additive effects


####################################################################################################################
# load simon's hub stuff
load("yeast-epistasis-paper/data/networkInfo.RData")
networkNeighbours.hubs <- networkNeighbours[networkNeighbours$nrNeighbours > 0,]
rownames(eff_size_only)[!(rownames(eff_size_only) %in% unique(networkNeighbours.hubs$trait))]
#
  
hub_mt_max <- aggregate(networkNeighbours.hubs,list(networkNeighbours.hubs$snpName),FUN=max)
# match add and epi in this list
nun_interactors <- numeric(length(var))
nun_hits <- numeric(length(var))
hit_snps_all <- list()
for(i in 1:length(var)){
  # is there a hit or not
  snp_now <- names(var)[i]
  snp_list <- hub_mt_max$snpName
  chr <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_now)
  chrs <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_list)
  index  <- c(chr,chrs) == chr
  r2_all <- r2fast(data = data,snpsubset = c(snp_now,snp_list)[index])
  diag(r2_all) <- 0
  r2_all[lower.tri(r2_all)] <- r2_all[upper.tri(r2_all)]
  if(any(r2_all[snp_now,] > 0.9)){ # r2 above 0.95 was regarded as the same loci
    hit_snps <- colnames(r2_all)[r2_all[snp_now,] > 0.9]
    hit_snps_all[[snp_now]] <- hit_snps
    nun_hits[i] <- length(hit_snps)
    cat(i,length(hit_snps),"\n")
    num_now <- max(hub_mt_max$nrNeighbours[match(hit_snps,table = hub_mt_max$snpName)])
  }else{
    num_now <- 0
  }
  nun_interactors[i] <- num_now
}
#save(nun_interactors,var,file="./results/181205-nun_interactors_and-var.RData")
hubs_plot <- names(var)[nun_interactors > 4]


###################################################################################
#dev.off()
########################################################################################################################
#####################################################################################################################
################### additive effects 

#layout(rbind(c(1,1,1),c(2,3,4),c(5,6,7)), widths=c(1,1,1))
par(mfrow=c(3,1),mar=c(4,4,2,2))
plot(1,ylim = c(-0.5,0.5),xlim = c(0,n*7),frame.plot = F,col="white",xaxt="n",xlab="",ylab="Additive effects")
text(x =1000,y=-0.5,labels = "All 311 additive QTL" )
library("RColorBrewer")
#display.brewer.all()
pal1 <- brewer.pal(9, "Set3")
pal2 <- brewer.pal(9, "Paired")

col2 <- c(pal1,pal2)#c("blue","tomato","darkolivegreen","purple")
grid(nx = 80,ny = 20)
for( i in 1:length(seq(1,nrow(eff_no_inter),3))){
  j <- seq(1,nrow(eff_no_inter),3)[i]
  points(7*(1:n-1)+1,out_sorted[j,],cex=0.4,pch=19,col=col2[i])
}
# points(7*(1:n-1)+1,out_sorted[1,],cex=0.5,pch=19,col=col2[1])
# points(7*(1:n-1)+1,out_sorted[4,],cex=0.5,pch=19,col=col2[2])
# points(7*(1:n-1)+1,out_sorted[7,],cex=0.5,pch=19,col=col2[3])
# points(7*(1:n-1)+1,out_sorted[10,],cex=0.5,pch=19,col=col2[4])

sig.snp <- which(as.character(out_inter_qtl[colnames(out_sorted),"GxE"]) =="No")
sig.x <- c(7*(1:n-1)+1)[sig.snp]
arrows(x0 = sig.x,x1 = sig.x,y0 = -0.43,y1 = -0.4,length = 0.02)
#legend("topright",legend = c("Manganese Sulfate","Lactate","Indol acetic acid","Raffinose"),fill=col2,bty = "n")
abline(h=0,lty="dashed",col="black",lwd=0.5)
#find_best(r2_cut = ,snp_list = ,genable = )
hubs_plot <- names(var)[nun_interactors > 4]
SNPs_hub_out_sorted <- match(hubs_plot,colnames(out_sorted))
#col3 <- c("blue","blue","blue","tomato","tomato","darkolivegreen","purple")
#abline(v=c(7*(1:n-1)+1)[SNPs_hub_out_sorted],lty="dashed",col="red",lwd=0.4)
arrows(x0 = c(7*(1:n-1)+1)[SNPs_hub_out_sorted],x1 = c(7*(1:n-1)+1)[SNPs_hub_out_sorted],y0 =0.43,y1 = 0.4,col="red",length = 0.03)
###################################
###################################
####
###################################

#layout(rbind(c(1,1,1),c(2,2,2)), widths=c(1,1,1))
par(mar=c(4,4,2,2))
plot(nun_interactors,var,frame.plot = F,cex=0.8,pch=19,col="blue",xlab="Number of interactors",ylab="Variance of estimated additive effects")
sum(nun_interactors > 4)

lm <- lm(var~nun_interactors)
abline(lm,lty="dashed",col="red")



###################################################################################

find_best2 <- function(snp_now,genable,snp_list,r2_cut=0.9){
  #snp_now="9714855_chrXIV_467219_A_G";genable=data;snp_list=qtl_select
  # select the chromsome
  chr <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_now)
  chrs <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_list)
  index  <- c(chr,chrs) == chr
  r2_all <- r2fast(data = genable,snpsubset = c(snp_now,snp_list)[index])
  if(max(r2_all[snp_now,],na.rm = T) > r2_cut){
    return(rownames(r2_all)[which.max(r2_all[snp_now,])])
    
  }else{
    warning(i," No hits with r2 above ",r2_cut,", the max is: ",max(r2_all[snp_now,],na.rm = T))
    return(NULL)
    
  }
}

#return number of interactors

get_ineractors <- function(snp_input = hubs_plot[5],trait_list= rownames(eff_size_only),net_list = networkNeighbours.hubs){
  num_list <- numeric(length(trait_list))
  for ( i in 1:length(trait_list)){
    trait_now <- trait_list[i]
    if(any(net_list$trait == trait_now)){
      net_list_now <- net_list[net_list$trait == trait_now,]
      hit <- find_best(snp_now = snp_input,genable = data,snp_list =net_list_now$snpName,r2_cut = 0.8 )
      if(length(hit) != 0){
        num <- net_list_now$nrNeighbours[net_list_now$snpName == hit]
      }else{
        num <- 0
      }
      
    }else{
      num <- 0
    }
    num_list[i] <- num
  }
  return(num_list)
}

## get interactor for all snps 

out <- list()

for(i in 1:length(hubs_plot)){
  out[[i]] <- get_ineractors(snp_input = hubs_plot[i],trait_list= rownames(eff_size_only),net_list = networkNeighbours.hubs)
  cat(i,"\n")
}
# par(mfrow=c(4,4))
# for(i in 1:11){
#   plot(out[[i]],eff_size_only[,hubs_plot[i]])
#   lm <- lm(eff_size_only[,hubs_plot[i]]~out[[i]])
#   abline(lm)
# }
#p_list <- numeric(11)
#par(mfrow=c(4,4))
x <- c()
y <- c()
for(i in 1:11){
  # lm <- lm(abs(eff_size_only[,hubs_plot[i]])~out[[i]])
  # lm_s <- summary(lm)
  # abline(lm,lty="dashed",col="red")
  # p <- format(lm_s$coefficients[2,4],scientific = T,digits = 2)
  # lab <- paste(hubs_plot[i]," (" ,p,")")
  # plot(out[[i]],abs(eff_size_only[,hubs_plot[i]]),frame.plot = F,xlab = "Num of connectors",ylab = "Effect size",main =lab)
  # p_list[i] <- lm_s$coefficients[2,4]
  x <- c(x,out[[i]])
  y <- c(y,abs(eff_size_only[,hubs_plot[i]]))
  #y <- c(y,eff_size_only[,hubs_plot[i]])
}

lm <- lm(y~x)
lm_s <- summary(lm)
abline(lm,lty="dashed",col="red")
p <- format(lm_s$coefficients[2,4],scientific = T,digits = 2)
lab <- paste("P value = ",p)
#plot(x,y,frame.plot = F,xlab = "Num of connectors",ylab = "abs(Effect size)",main =lab,pch=19,cex=0.8,col="blue")
#abline(lm,lty="dashed",col="red")


