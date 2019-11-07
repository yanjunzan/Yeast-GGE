##########################################################################################

# This is the analysis code for manuscript  "Dynamic genetic architecture of yeast growth response
#     to environmental perturbation shed light on origin of hidden genetic variation" Zan and Carlborg
#######################################################################################




#######################################################################################
#>>>>>>>>>>> This study is a continuation of two previous studies Bloom et al Nat comm and Forsberg et al 2017 Nat Genetic. 
#Raw and intermediate data from previous study are used in our study, which are avaiable at:
# 1. Forsberg, S. K. G., Bloom, J. S., Sadhu, M. J., Kruglyak, L. & Carlborg, Ö. Accounting for genetic interactions improves modeling of individual quantitative trait phenotypes in yeast. Nat. Genet. 49, 497–503 (2017).
# 2. Bloom, J. S. et al. Genetic interactions contribute less than additive effects to quantitative trait variation in yeast. Nat. Commun. (2015). doi:10.1038/ncomms9712
#######################################################################################



#######################################################################################
#>>>>>>>>>>> To be able to run the code please make sure the folder structure are:
##

#submission
# -bin
# -results
# -doc
# -data
# Please note some of the analysis are splited as external R script, which generate variable that are necessary for success running of
# downstream scrpts. Therefore, please run these line by line
#######################################################################################


#######################################################################################
#### analysis starts here
setwd("~/Dropbox/home_bin/yeast/submission/") # please modify this accordingly
# a few functions need to be sourced first.

source("./bin/Functions.yeast.R")

#######################################################################################
#>>>>>>>>>>> First download and format the raw data <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######################################################################################

#1 The raw data is avaiable as supplementary file in a previous paper -Forsberg et al 2017 Nature Genetics

#nevagate to project folder -submission- and in bash run:
git clone https://github.com/simfor/yeast-epistasis-paper.git

#2 some raw data is avaiable as supplementary file in a previous paper -Bloom et al Nature communication
git clone https://github.com/joshsbloom/4000BYxRM.git


#######################################################################################
# format the data into a genable format

require(GenABEL)
load("./yeast-epistasis-paper/data/genotypes.RData")
load("./yeast-epistasis-paper/data/phenotypes.RData")

#ped file
gdata.ped <- gdata
gdata.ped[gdata.ped == 1] <- 2
gdata.ped[gdata.ped == -1] <- 1
gdata.ped <- gdata.ped[, rep(x = colnames(gdata.ped), times = rep(x = 2, times = ncol(gdata.ped)))] #duplicate every column in the genotype matrix
ids <- rownames(gdata.ped)
gdata.ped <- data.frame(ids, ids, "0 0 0 -9", gdata.ped)
write.table(x = gdata.ped, file = "./yeast-epistasis-paper/data/gdata.ped", sep = "\t", quote = F, row.names = F, col.names = F)

#map file
markers <- colnames(gdata)
tmp <- strsplit(markers, split = "_")
map_SacCer3 <- as.numeric(sapply(X = tmp, FUN = function(x){x[3]}))
chr <- sapply(X = tmp, FUN = function(x){x[2]})
tmp <- unique(chr)
for(i in 1:length(tmp)){
  chr[chr == tmp[i]] <- i
}
chr <- as.numeric(chr)
gdata.map <- data.frame(chrom = chr, name = markers, position = map_SacCer3)
write.table(x = gdata.map, file = "./yeast-epistasis-paper/data/gdata.map", sep = "\t", quote = F, row.names = F, col.names = T)

#phenotype file
pheno <- data.frame(id = ids, sex = 1, matrix(nrow = length(pheno_raw[[1]]), ncol = 60))
names <- c() #the names of the traits
names[seq(from = 1, to = 58, by = 3)] <- paste(paste("trait", 1:20, sep = "_"), 1, sep = ".")
names[seq(from = 2, to = 59, by = 3)] <- paste(paste("trait", 1:20, sep = "_"), 2, sep = ".")
names[seq(from = 3, to = 60, by = 3)] <- paste(paste("trait", 1:20, sep = "_"), "mean", sep = ".")
names(pheno)[3:62] <- names
j <- 3
for(i in 1:length(pheno_raw)){
  trait <- pheno_raw[[i]]
  trait1 <- sapply(X = trait, FUN = function(x){x[1]})
  trait2 <- sapply(X = trait, FUN = function(x){x[2]})
  pheno[,j] <- trait1
  j <- j+1
  pheno[,j] <- trait2
  j <- j+1
  pheno[,j] <- sapply(X = trait, FUN = function(x){ mean(c(x[1], x[2]), na.rm = T) })
  j <- j+1
}
write.table(x = pheno, file = "./yeast-epistasis-paper/data/phenotypes_meanAndRep", quote = F, sep = "\t", row.names = F)

#Create raw file and load it as a GenABEL object
convert.snp.ped(pedfile = "./yeast-epistasis-paper/data/gdata.ped", mapfile = "./yeast-epistasis-paper/data/gdata.map", outfile = "./yeast-epistasis-paper/data/gdata.raw")
data <- load.gwaa.data(phenofile = "./yeast-epistasis-paper/data/phenotypes_meanAndRep", genofile = "./yeast-epistasis-paper/data/gdata.raw")
save(list = "data", file = "./data/yeast.GenABEL.Data")
load("./data/yeast.GenABEL.Data")

#######################################################################################
#>>>>>>>>>>> For each enviroment compile a list of additive loci
#######################################################################################
# Bloom et al performed QTL mapping for each medium independently, and provided a list of 
# QTL for each medium(avaialbe as Supplementary information ncomms9712-s4.txt). As Loci mapped
# for these medium might be in LD with each other, we first run a multi-locus analysis to select a
# independent set of loci for each medium and compile a list of additive loci contributing to growth in at least 1 enviroment.

### in GenABLE data the phenotype name need to match the medium
trait.name <- paste(paste("trait", 1:20, sep = "_"), "mean", sep = ".")
names(trait.name) <- names(pheno_raw)
save(trait.name,file="./data/Trait2name.RData")

load("./data/yeast.GenABEL.Data")
load("./data/Trait2name.RData")
require(GenABEL)
# pick 
map.all <- GenABEL::map(data)
# Extract all additive snp for picked phenotype
snp.all <- snpnames(data)
mrk <- read.table("./doc/ncomms9712-s4.txt",header = T,sep = "\t",stringsAsFactors = F)
phe.pick <- unique(mrk$trait) ## select a few 
mrk.sub <- subset(mrk,subset=mrk$trait %in% phe.pick)
chr.loca <- paste(mrk.sub$chr,mrk.sub$pos,sep="_")
get.index <- function(x) return(snp.all[grep(pattern = paste0(x,"_"),x = snp.all)])
chr.loca.full <- unique(unlist(lapply(chr.loca, FUN = get.index)))
mrk.idp4 <- remove_tight.loci(SNPs =chr.loca.full,genable = data,cut_dis = 20e3,cut_r2 = 0.9) # loci located with 20kb and r2 >0.9 are viewed as the same loci
#length(unique(mrk.idp3$snp))
length(unique(mrk.idp4$snp))
geno <- as.double.gwaa.data(data[,unique(mrk.idp4$snp)])
##################################################
#>>> Please note that part is paralleled across 7 threads, and will run very long time
##################################################

library(doSNOW)
library(foreach)
ptm <- proc.time()
#BE_all_traits <- 
cl <- makeCluster(7) # update this accroding to your computer
registerDoSNOW(cl)
foreach(i = 1:length(phe.pick)) %dopar% {
  source(file = "./bin/Functions.yeast.R")
  #source(file = "~/Dropbox/home_bin/yeast/Functions.yeast.R")
  require(GenABEL)
  triat_now <- phe.pick[i]
  out <- BE_analysis_no_cov(phe =phdata(data)[,get.Gen.id(t.name = triat_now)],data = geno)
  #return(list("BE"=out,"trait"=triat_now))
  output <- list("BE"=out,"trait"=triat_now)
  output_name <- paste0("./results/181019_Be_all_",triat_now,".RData") 
  save(output,file = output_name)
}
stopCluster(cl)
####################################### Summarise the result ##############################l
local <- list.files(path ="./results/",pattern = "181019_Be_all_trait")
dip <- list.files(path ="./results/",pattern = "181019_Be_all" )
inter <- list.files(path ="./results/",pattern = "181019_Be_all_inter" )
dip <- dip[!(dip %in% c(local,inter))]

trait_local <- gsub(pattern = "181019_Be_all_trait(.*)\\.RData",replacement = "\\1",local)
trait_dip <- gsub(pattern = "181019_Be_all_(.*)\\.RData","\\1",dip)
trait_dip <- trait_dip[!(trait_dip %in% trait_local )]

link_local <- paste0("181019_Be_all_trait",trait_local,".RData")
link_dip <-  paste0("181019_Be_all_",trait_dip,".RData")
length(unique(c(link_dip,link_local)))
# generate a list
Trait_all  <- c(trait_local,trait_dip)
link_all <- c(link_local,link_dip)
trait_name <- c()
qtl <- list() # this is a list object of length 20, storing the QTL mapped for each medium
for( i in 1:20){
  link_now  <- paste0("./results/",link_all[i])
  qtl_now <- get(load(link_now))
  if(qtl_now$trait !=Trait_all[i])
    stop("trait does not match")
  qtl[[i]] <- qtl_now$BE$out5
  cat(i,"\n")
  
}
qtl_all <- unique(unlist(qtl)) # this is a list of QTL detected for at least 1 medium
names(qtl) <- Trait_all
save(qtl_all,qtl,file = "./results/181031_all_additive_qtl.RData")
load("./results/181031_all_additive_qtl.RData")


##################################################################################################################################
##>>>>>>>>>>>>>>>>>>>>testing if QTLs mapped for other medium are making significant contribution to growth at focal medium
##################################################################################################################################
p_value <- numeric(20)
#qtl_select <- unique(c(BE.yeast.1$out5, BE.yeast.2$out5, BE.yeast.3$out5, BE.yeast.4$out5))
qtl_select <- unique(unlist(qtl))
geno_full <- as.double.gwaa.data(data[,qtl_select])
require(lmtest)
phe.pick <- names(qtl)#c("ManganeseSulfate", "Lactate","Raffinose","IndolaceticAcid")
for(i in 1:20){
  ph <- phdata(data)[,get.Gen.id(t.name = phe.pick[i])]
  lm_full  <- lm(ph~geno_full)
  geno_now <- as.double.gwaa.data(data[,qtl[[i]]])
  lm_now <- lm(ph~geno_now)
  lmtest_now  <- lmtest::lrtest(lm_full,lm_now)
  p_value[i] <- lmtest_now$`Pr(>Chisq)`[2]
}
sum(p_value < 0.05/20)
pass <- rep("Yes",20)
pass[p_value > 0.05/20] <-rep("No",sum(p_value > 0.05/20)) 

num <- c() # number of QTL for each medium
for( i in 1:length(qtl)){
  num <- c(num,length(unique(qtl[[i]])))
}

out <- data.frame("Trait"=phe.pick,"Number.QTL"=num,"P.values"=format(p_value,scientific = T,digits = 3),"Significant"=pass)
write.table(out,file = "./results/181030_test_for_polygenetic.txt",row.names = F,col.names = T,quote = F,sep="\t")


##################################################################################################################################
#>>>>>>>>>>>testing for QTL by medium interaction
##################################################################################################################################
trait_picked <- trait.name[-c(18,19)]
qtl_select <- unique(unlist(qtl)) 
out_inter_qtl <- fit.GbyE.allele(trait = trait_picked,genable = data,snp_test =qtl_select ) # a function to do the test
#sum(out_inter_qtl$GxE=="Yes")/nrow(out_inter_qtl)
write.table(out_inter_qtl,file = "./results/181031_Table_S2_P_value_QTLxE.txt", sep = "\t", row.names = F,col.names = T,quote = F) #The pairwise interactions mapped in Bloom2015

##################################################################################################################################
#>>>>>>>>>>>Visulize phenotypic PCA plot across medium
##################################################################################################################################
source("./bin/181002_Figure_phenotypic_cor.R")


#######################################################################################
#>>>>>>>>>>>  Build a across enviroment epistatic network
#######################################################################################
# Bloom et al 2015 nature comm mapped QTL-QTL interactions and Forsberg et al nature genetics

source("./bin/190305_summary_network_build.R")

#######################################################################################
#>>>>>>>>>>>>>>>>>>>>>>>>>>>> Making figures<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######################################################################################

## Please note that some figures are ploted in R and pieced together using adobe illustrator

##################################
# >>>>>>>>>>>>>>>>Figure 1<<<<<<<<
##################################
source("./bin/190412_Figure1.R")

##################################
# >>>>>>>>>>>>>>>>Figure 2<<<<<<<<
##################################
setwd("~/Dropbox/home_bin/yeast/submission/")
source("./bin/190419_Figure2.R")

##################################
# >>>>>>>>>>>>>>>>Figure 3<<<<<<<<
##################################
source("./bin/190507_figure3_update.R")
source("./bin/191007_figure3.R")

##################################
# >>>>>>>>>>>>>>>>Figure 4<<<<<<<<
##################################
source("./bin/190416_Figure4.R")


##################################
# >>>>>>>>>>>>>>>>Figure 5<<<<<<<<
##################################
source("./bin/190602_Figure5_additive_effects_connectivity.R")


##################################
# >>>>>>>>>>>>>>>>Figure 6<<<<<<<<
##################################
source("./bin/190418_Figure6.R")


##################################
# >>>>>>>>>>>>>>>>S Figure s<<<<<<<<
##################################
source("./bin/190424_S_Figures.R")

# some s Figures are generate by source("./bin/190416_Figure4.R")
