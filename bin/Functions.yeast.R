
LDheatmap.GenABEL.snp.my <- function(data, nth = 5, SNP.name = NULL, returnGeno = F,LDmeasure="r",col,...){
  if(!require(LDheatmap)){
    stop("Could not load package LDheatmap")
  }
  #id <- length(SNP.name) %% nth ==0
  gtps <- as.genotype.gwaa.data(data[,SNP.name])
  dist <- rank(as.vector(data@gtdata@map[SNP.name]))
  
  #gtps <- as.genotype.gwaa.data(data)
  #dist <- rank(as.vector(data@gtdata@map))
  
  
  if(returnGeno){
    return(gtps)
  }
  #   names <- names(gtps)
  #   names <- names[seq(1,length(names),by=2)]
  #   coords <- region[seq(1,length(region),by=2)]
  #   labels <- paste(coords," bp (",names,")",sep="")
  #mycols <- colorRampPalette(colors=c("red","orange","white","slateblue"))
  #mycols <- colorRampPalette(colors=c("red","orange","white","green"))
  #plot.LD(gtps, dist, type="physical", col=mycols(1000), LD.measure="r", SNPnames=names)
  MyLD <- LDheatmap(gtps, dist, "physical", color=col, LDmeasure=LDmeasure, title=substitute(paste("Pairwise LD measured by ", LDmeasure)), 
                    name="myLDgrob", SNP.name=SNP.name, ...)
  return(MyLD)
}
#
get.Gen.id <- function(t.name,input=get(load("./data/Trait2name.RData")) ){
  return(input[paste("1_",t.name,"_1",sep="")])
}


###############
BE_analysis_no_cov <- function(phe,data,fdr=c(0.05,0.2),mfactor=1){ 
  # data is a data frame contain cov genotype, cov is the name of covariants
  # build full model 
  nna <- complete.cases(cbind(phe,data))
  mrk <- colnames(data)
  id.full<- seq(1:length(mrk))
  geno_add <- paste("as.numeric(data[,",id.full,"])",sep="",collapse="+")
  phe<-phe[nna]
  data <- data[nna,]
  #test_fx1<-fx1[nna]
  fm <- as.formula(paste("phe ~ ",geno_add,sep=""))
  reg.lm <- lm(fm, y=TRUE)
  
  # min model
  #id.min<-c(1:ncol(geno_reg))[colnames(geno_reg) %in% sigMrk_sub]
  #geno_bc.add<-paste("as.numeric(geno_reg[,",id.min,"])",sep="",collapse="+")
  #fm.min<- as.formula(paste("test_phe ~ test_fx1 +",geno_bc.add,sep=""))
  #min.lm <- lm(fm.min,y=TRUE)
  #min.lm <- lm(test_phe ~ test_fx1)
  #Perform Backward-Elimination in the original data at different adaptive FDR thresholds
  fitFDR5<- BEFDR_nofix( maximal.lm = reg.lm, FDR.q = fdr[1],mfactor=mfactor)
  
  if(sum(grepl("data",rownames(summary(fitFDR5)$coefficients)))>0){
    terms5<-rownames(summary(fitFDR5)$coefficients)[grep("data",rownames(summary(fitFDR5)$coefficients))]
    idx5<-unlist(strsplit(terms5,",|[]]"))
    name5<-colnames(data)[as.numeric(idx5[seq(2,length(idx5),3)])]
    #id.in5<- name5[(name5%in% mrksT)]
    #mrks5<-c(mrks5,id.in5)
    #out5[[as.numeric(i)]]<- id.in5
  }
  
  fitFDR20<-BEFDR_nofix( maximal.lm = reg.lm, FDR.q = fdr[2],mfactor=mfactor)
  if(sum(grepl("data",rownames(summary(fitFDR20)$coefficients)))>0){
    terms20 <-rownames(summary(fitFDR20)$coefficients)[grep("data",rownames(summary(fitFDR20)$coefficients))]
    idx20 <- unlist(strsplit(terms20,",|[]]"))
    name20 <-colnames(data)[as.numeric(idx20[seq(2,length(idx20),3)])]
    #id.in20 <- gsub(pattern = "X(.*)",replacement = "\\1",x = name20)
    #mrks20<-c(mrks20,id.in20)
    #out20[[as.numeric(i)]]<- id.in20
  }
  return(list("out5"=name5,"out20"=name20))
}



BEFDR_nofix <- function ( maximal.lm, FDR.q, mfactor = 1) {
  compute.Lambda <- function(k, m, Q) {
    i <- c(1:k)
    return((1/(k + 1)) * sum(qnorm((Q/2) * (i/(m + 1 - i * 
                                                 (1 - Q))))^2))
  }
  get.model.size <- function(a.lm) {
    require(MASS)
    return(extractAIC(a.lm)[1] - 1)
  }
  require(MASS)
  #the.scope <- list(lower = minimal.lm, upper = maximal.lm)
  m <- mfactor * get.model.size(maximal.lm)
  new.model.size <- get.model.size(maximal.lm)
  for (i in 1:m) {
    old.model.size <- new.model.size
    Lambda <- compute.Lambda(k = old.model.size - 1, m, Q = FDR.q)
    ## no min model
    new.model <- stepAIC(maximal.lm, direction = "backward", k = Lambda, trace = FALSE)
    new.model.size <- get.model.size(new.model)
    if (new.model.size >= old.model.size) 
      break
  }
  new.lm <- lm(new.model)
  return(new.lm)
}

################################### With cov

BEFDR <- function (minimal.lm, maximal.lm, FDR.q, mfactor = 1) {
  compute.Lambda <- function(k, m, Q) {
    i <- c(1:k)
    return((1/(k + 1)) * sum(qnorm((Q/2) * (i/(m + 1 - i * 
                                                 (1 - Q))))^2))
  }
  get.model.size <- function(a.lm) {
    require(MASS)
    return(extractAIC(a.lm)[1] - 1)
  }
  require(MASS)
  the.scope <- list(lower = minimal.lm, upper = maximal.lm)
  m <- mfactor * get.model.size(maximal.lm)
  new.model.size <- get.model.size(maximal.lm)
  for (i in 1:m) {
    old.model.size <- new.model.size
    Lambda <- compute.Lambda(k = old.model.size - 1, m, Q = FDR.q)
    new.model <- stepAIC(maximal.lm, direction = "backward", 
                         scope = the.scope, k = Lambda, trace = FALSE)
    new.model.size <- get.model.size(new.model)
    if (new.model.size >= old.model.size) 
      break
  }
  new.lm <- lm(new.model)
  return(new.lm)
}

BE_analysis_cov<- function(cov,phe,data,fdr=c(0.05,0.2),mfactor=1){ 
  # data is a data frame contain cov genotype, cov is the name of covariants
  # build full model 
  nna <- complete.cases(cbind(phe,geno.mat))
  mrk <- colnames(data)
  id.full<- seq(1:length(mrk))
  geno_add <- paste("as.numeric(data[,",id.full,"])",sep="",collapse="+")
  phe<-phe[nna]
  data <- data[nna,]
  #test_fx1<-fx1[nna]
  fm <- as.formula(paste("phe ~ ",geno_add,sep=""))
  reg.lm <- lm(fm, y=TRUE)
  
  # min model
  id.min<-c(1:ncol(data))[mrk %in% cov]
  geno_bc.add<-paste("as.numeric(data[,",id.min,"])",sep="",collapse="+")
  fm.min<- as.formula(paste("phe ~",geno_bc.add,sep=""))
  min.lm <- lm(fm.min,y=TRUE)
  
  #Perform Backward-Elimination in the original data at different adaptive FDR thresholds
  fitFDR5<-BEFDR( minimal.lm = min.lm, maximal.lm = reg.lm, FDR.q = fdr[1],mfactor=mfactor)
  
  if(sum(grepl("data",rownames(summary(fitFDR5)$coefficients)))>0){
    terms5<-rownames(summary(fitFDR5)$coefficients)[grep("data",rownames(summary(fitFDR5)$coefficients))]
    idx5<-unlist(strsplit(terms5,",|[]]"))
    name5<-colnames(data)[as.numeric(idx5[seq(2,length(idx5),3)])]
    #id.in5<- name5[!(name5%in% sigMrk_sub)]
    #mrks5<-c(mrks5,id.in5)
    #out5[[as.numeric(i)]]<- id.in5
  }
  
  fitFDR20<-BEFDR( minimal.lm = min.lm, maximal.lm = reg.lm, FDR.q = fdr[2],mfactor=mfactor)
  
  if(sum(grepl("data",rownames(summary(fitFDR20)$coefficients)))>0){
    terms20 <-rownames(summary(fitFDR20)$coefficients)[grep("data",rownames(summary(fitFDR20)$coefficients))]
    idx20 <- unlist(strsplit(terms20,",|[]]"))
    name20 <-colnames(data)[as.numeric(idx20[seq(2,length(idx20),3)])]
    #id.in20 <- gsub(pattern = "X(.*)",replacement = "\\1",x = name20)
    #mrks20<-c(mrks20,id.in20)
    #out20[[as.numeric(i)]]<- id.in20
  }
  
  return(list("out5"=name5,"out20"=name20))
}
# function 
# selecting indpendent mrk
get.indpen.mrk <- function(r2.mat=r$LDmatrix,cut=0.9){
  diag(r2.mat) <-0
  r2.mat[lower.tri(r2.mat)] <- r2.mat[upper.tri(r2.mat)]
  
  all.marker <- colnames(r2.mat)
  all.chr <- gsub(pattern = ".*_(chr.*)_\\d+_.*",replacement = "\\1",x = all.marker)
  chr.level <- unique(all.chr)
  snp <-c()
  snp.rm <-c()
  for( chr in chr.level){
    all.marker.sub <- all.marker[all.chr %in% chr]
    r2.mat.sub <- r2.mat[all.marker.sub,all.marker.sub]
    # all comb
    if(length(all.marker.sub)>2){
      all.com <- combn(x = all.marker.sub,m = 2)
      yes <- rep(F,ncol(all.com))
      for( j in 1:ncol(all.com)){
        if(r2.mat.sub[all.com[1,j],all.com[2,j]] > cut) # add distance and save these two marker, which is saved and which is removed
          yes[j] <- T
      }
      mrk.rm <- unique(all.com[1,yes])
      snp.rm <- c(snp.rm,mrk.rm)
      snp <- c(snp,all.marker.sub[!(all.marker.sub %in% mrk.rm)])
    }
  }
  ### make sure removed is all right
  #r2.mat.rm <- r2.mat[snp.rm,]
  #if(any(r2.mat.rm > cut))
  return(snp)
}

remove_tight.loci <- function(SNPs=full_name_g,genable=data,cut_dis=10e3,cut_r2=0.9){
  if(!require(igraph))
    require(igraph)
  if(!require(GenABEL))
    require(GenABEL)
  # rule1 physical distance
  pos.all <- map(genable[,unique(SNPs)])
  names(pos.all) <- unique(SNPs)
  pos <- pos.all[SNPs]
  diffs <- abs(outer(pos, pos, FUN = "-")) #all pairwise differences in physical position
  diffs[lower.tri(diffs)] <- 1000000
  diag(diffs) <- 1000000
  #identical(rownames(diffs),rownames(r2))
  ## rule 2 linkage
  r2 <- r2fast(data = genable,snpsubset = SNPs)
  diag(r2) <- 0
  r2[lower.tri(r2)] <- 0
  ## rule3 chromsome
  chrs <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = SNPs)
  chrs.order <- match(chrs,unique(chrs))
  diffs.chr <- abs(outer(chrs.order, chrs.order, FUN = "-")) #all pairwise differences in index
  diffs.chr[lower.tri(diffs.chr)] <- 1
  diag(diffs.chr) <- 1
  
  clump <- which(diffs < cut_dis & r2 > cut_r2 & diffs.chr==0) #differences smaller than n
  if(length(clump) <=1)
    stop("nothing")
  ## diffs.col here is wrong should be row and should be clump %% nrow(diffs)
  ## diffs.row should be ceiling(clump/ncol(diffs))
  diffs.col <- ceiling(clump/ncol(diffs)) #differences smaller than n, columns in the distance matrix
  diffs.row <- clump - (diffs.col - 1)*nrow(diffs) #differences smaller than n, rows in the distance matrix
  #mycol <- clump %% nrow(diffs)
  #myrow <- ceiling(clump/ncol(diffs))
  
  #require(igraph)
  diffs.graph <- graph_from_data_frame(data.frame(diffs.col, diffs.row))
  diffs.clust <- clusters(diffs.graph)
  
  snps <- SNPs
  relation <- list()
  name <- c()
  for(i in 1:diffs.clust$no){
    nodesInCluster <- as.numeric(names(diffs.clust$membership[diffs.clust$membership == i]))
    #chrs[nodesInCluster] <- chrs[nodesInCluster[1]]
    relation[[i]] <- snps[nodesInCluster]
    snps[nodesInCluster] <- snps[nodesInCluster[1]]
    name <- c(name, snps[nodesInCluster[1]])
    #pos[nodesInCluster] <- pos[nodesInCluster[1]]
  }
  names(relation) <- name
  if(identical(snps,SNPs)){
    warning("no change ")
    return(list("snp"=snps,"re"=relation))
  }else{
    return(list("snp"=snps,"re"=relation))
  }
}

get.snpnames <- function(chr.loca,allsnp=snp.all){
  #chr.loca <- paste(data.sub$chr,data.sub$pos,sep="_")
  get.index <- function(x) return(allsnp[grep(pattern = paste0(x,"_"),x = allsnp)])
  return(c(unlist(lapply(chr.loca, FUN = get.index))))
}

find_best <- function(snp_now,genable,snp_list,r2_cut=0.9){
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
    return(rownames(r2_all)[which.max(r2_all[snp_now,])])
    
  }
}


cmp_snps <- function(add_snp,epi_snp,genable,r2_cut=0.9){
  #add_snp <- tmp.before$addSNP1[i]
  #epi_snp <-tmp.before$epiSNP1[i]
  add_chr <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = add_snp)
  epi_chr <-gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = epi_snp)
  r2 <- r2fast(data = genable,snpsubset = c(add_snp,epi_snp))[1,2]
  if(add_chr==epi_chr & r2 >r2_cut){
    return(T)
  }else{
    warning(" No hits with r2 above ",r2_cut,  " r2 is" ,r2, "or they are in different chr")
    return(F)
  }
}

fit.GbyE.allele <- function(trait= trait_picked,genable=data,snp_test=qtl_select){
  if(!require(reshape2))
    require(reshape2)
  if(!require(GenABEL))
    require(GenABEL)
  warning("Phenotypes must have been scaled","\n")
  trait2 <- paste0(gsub(pattern = "(.*)\\.mean","\\1",trait_picked),c(rep(".1",length(trait)),rep(".2",length(trait))))
  num_trait <- length(trait2)
  geno <- as.double.gwaa.data(data[,snp_test])
  ph <- phdata(genable)[,c("id",trait2)]
  na <- complete.cases(cbind.data.frame(ph,geno))  
  ph.v <- melt(data = ph,id.vars="id",measure.vars=2:ncol(ph))
  #update variable
  ph.v$variable <- gsub(pattern = "(.*)\\.\\d{1}",replacement = "\\1",x = ph.v$variable )
  geno.all <- do.call("rbind", replicate(num_trait, geno, simplify = FALSE))
  if(!identical(rownames(geno.all),ph.v$id))
    stop("id in the genotype matrix and phentoype matrix do not match")
  
  # build full model
  id.full<- 1:ncol(geno.all)
  geno_add <- paste("as.numeric(geno.all[nna,",id.full,"])",sep="",collapse="+")
  nna <- complete.cases(cbind(ph.v$value,ph.v$variable,geno.all))
  test_phe <- ph.v$value[nna]
  test_fx1 <- ph.v$variable[nna]
  fm <- as.formula(paste("test_phe ~ test_fx1 +",geno_add,sep=""))
  lm <- lm(fm)
  allele.all <- levels(as.factor(colnames(geno.all)))
  num <- length(allele.all)
  p<- numeric(num)
  
  for( i in 1:num){
    fm2 <- as.formula(paste("`test_phe` ~ `test_fx1` +",geno_add,"+ `test_fx1`:as.numeric(geno.all[nna,",i,"])",sep=""))
    #fm2 <- as.formula(paste("`test_phe` ~ `test_fx1` +",geno_add,sep=""))
    #lm(fm2)
    if(!require(lmtest))
      require(lmtest)
    lm.inter2 <- try(silent = T,lm(fm2))
    if(!inherits(lm.inter2,"try-error")){
      an <- lrtest(lm,lm.inter2)
      p[i] <- an$`Pr(>Chisq)`[2]
    }
    cat(i,"\n")
  }
  pos <- as.numeric(gsub(pattern = "\\d+_chr.*_(\\d+)_.*_.*",replacement = "\\1",x = snp_test))
  Inter <- ifelse(p<0.05/length(snp_test),"Yes","No")
  out <- cbind.data.frame(colnames(geno),chromosome(data[,snp_test]),pos,format(p,scientific = T,digits = 3),Inter)
  colnames(out) <- c("SNP_name","Chromosome","position","P_vlaue","GxE")
  return(out)
  #save(p,file="results/180727_lm_model_allele.RData") ## 89 allele
  #load("results/180727_lm_model_allele.RData")
}
est.genetic.effect <- function(trait= trait_picked,genable=data,snp_test=qtl_select){
  if(!require(reshape2))
    require(reshape2)
  warning("Phenotypes must have been scaled","\n")
  num_trait <- length(trait)
  geno <- as.double.gwaa.data(data[,snp_test])
  ph <- phdata(genable)[,c("id",trait_picked)]
  na <- complete.cases(cbind.data.frame(ph,geno))
  # ph.v <- melt(data = ph[na,],id.vars="id",measure.vars=2:5)
  # geno.all <- do.call("rbind", replicate(num_trait, geno[na,], simplify = FALSE))
  # if(!identical(rownames(geno.all),ph.v$id))
  #   stop("id in the gentype matrix and phentoype matrix do not match")
  # 
  # # only additive no interaction fitted
  # lm <- lm(ph.v)
  
  ## iteratively fit lm and extract the output
  out_put <- data.frame(array(NA,c(length(trait)*3,length(qtl_select)+1)))
  for( i in 1:length(trait)){
    out <- lm(ph[na,i+1]~geno[na,])
    out.s <- summary(out)
    index <- seq((3*i-2),3*i)
    out_put[index,] <- t(out.s$coefficients[,c(1,2,4)])
  }
  colnames(out_put) <- gsub(pattern = "geno\\[.*\\](.*)",replacement = "\\1",row.names(out.s$coefficients))
  return(out_put)
}

get_loci_add_epi <- function(tmp,trait){
  #tmp <- tmp.clumped;trait="IndolaceticAcid"
  epi_loci <- unique(c(tmp[tmp$trait==trait,"epiSNP1"],tmp[tmp$trait==trait,"epiSNP2"]))
  add_loci <- unique(c(tmp[tmp$trait==trait,"addSNP1"],tmp[tmp$trait==trait,"addSNP2"]))
  return(list("epi"=epi_loci,"add"=add_loci))
}
find_closest <- function(x,snp.all=snp.all){
  #x <- "chrVII_122421"
  chr_x <- gsub("(.*)_.*",replacement = "\\1",x)
  pos_x <- as.numeric(gsub(".*_(.*)",replacement = "\\1",x))
  chr_all <- gsub(".*_(.*)_.*_.*_.*",replacement = "\\1",snp.all)
  pos_all <-  as.numeric(gsub(".*_.*_(.*)_.*_.*",replacement = "\\1",snp.all))
  
  dis <- abs(pos_all[chr_all==chr_x] - pos_x)
  pos_min <- pos_all[chr_all==chr_x][which.min(dis)[1]]
  return(paste0(chr_x,"_",pos_min))
}


find_all_above <- function(snp_now,genable,snp_list,r2_cut=0.9){
  #snp_now="9714855_chrXIV_467219_A_G";genable=data;snp_list=qtl_select
  # select the chromsome
  chr <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_now)
  chrs <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = snp_list)
  index  <- c(chr,chrs) == chr
  r2_all <- r2fast(data = genable,snpsubset = c(snp_now,snp_list)[index])
  if(max(r2_all[snp_now,],na.rm = T) > r2_cut){
    return(rownames(r2_all)[which(r2_all[snp_now,] > r2_cut)])
    
  }else{
    warning(i," No hits with r2 above ",r2_cut,", the max is: ",max(r2_all[snp_now,],na.rm = T))
    return(rownames(r2_all)[which.max(r2_all[snp_now,])])
    
  }
}

c.z.hglm <- function(kin){
  relmat <- kin
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d))
  return(Z)
}



LDheatmap.GenABEL.snp.my <- function(data, nth = 5, SNP.name = NULL, returnGeno = F,LDmeasure="r",col,...){
  if(!require(LDheatmap)){
    stop("Could not load package LDheatmap")
  }
  #id <- length(SNP.name) %% nth ==0
  gtps <- as.genotype.gwaa.data(data[,SNP.name])
  dist <- rank(as.vector(data@gtdata@map[SNP.name]))
  
  #gtps <- as.genotype.gwaa.data(data)
  #dist <- rank(as.vector(data@gtdata@map))
  
  
  if(returnGeno){
    return(gtps)
  }
  #   names <- names(gtps)
  #   names <- names[seq(1,length(names),by=2)]
  #   coords <- region[seq(1,length(region),by=2)]
  #   labels <- paste(coords," bp (",names,")",sep="")
  #mycols <- colorRampPalette(colors=c("red","orange","white","slateblue"))
  #mycols <- colorRampPalette(colors=c("red","orange","white","green"))
  #plot.LD(gtps, dist, type="physical", col=mycols(1000), LD.measure="r", SNPnames=names)
  MyLD <- LDheatmap(gtps, dist, "physical", color=col, LDmeasure=LDmeasure, title=substitute(paste("Pairwise LD measured by ", LDmeasure)), 
                    name="myLDgrob", SNP.name=SNP.name, ...)
  return(MyLD)
}
#
get.Gen.id <- function(t.name,input=get(load("./data/Trait2name.RData")) ){
  return(input[paste("1_",t.name,"_1",sep="")])
}


###############
BE_analysis_no_cov <- function(phe,data,fdr=c(0.05,0.2),mfactor=1){ 
  # data is a data frame contain cov genotype, cov is the name of covariants
  # build full model 
  nna <- complete.cases(cbind(phe,data))
  mrk <- colnames(data)
  id.full<- seq(1:length(mrk))
  geno_add <- paste("as.numeric(data[,",id.full,"])",sep="",collapse="+")
  phe<-phe[nna]
  data <- data[nna,]
  #test_fx1<-fx1[nna]
  fm <- as.formula(paste("phe ~ ",geno_add,sep=""))
  reg.lm <- lm(fm, y=TRUE)
  
  # min model
  #id.min<-c(1:ncol(geno_reg))[colnames(geno_reg) %in% sigMrk_sub]
  #geno_bc.add<-paste("as.numeric(geno_reg[,",id.min,"])",sep="",collapse="+")
  #fm.min<- as.formula(paste("test_phe ~ test_fx1 +",geno_bc.add,sep=""))
  #min.lm <- lm(fm.min,y=TRUE)
  #min.lm <- lm(test_phe ~ test_fx1)
  #Perform Backward-Elimination in the original data at different adaptive FDR thresholds
  fitFDR5<- BEFDR_nofix( maximal.lm = reg.lm, FDR.q = fdr[1],mfactor=mfactor)
  
  if(sum(grepl("data",rownames(summary(fitFDR5)$coefficients)))>0){
    terms5<-rownames(summary(fitFDR5)$coefficients)[grep("data",rownames(summary(fitFDR5)$coefficients))]
    idx5<-unlist(strsplit(terms5,",|[]]"))
    name5<-colnames(data)[as.numeric(idx5[seq(2,length(idx5),3)])]
    #id.in5<- name5[(name5%in% mrksT)]
    #mrks5<-c(mrks5,id.in5)
    #out5[[as.numeric(i)]]<- id.in5
  }
  
  fitFDR20<-BEFDR_nofix( maximal.lm = reg.lm, FDR.q = fdr[2],mfactor=mfactor)
  if(sum(grepl("data",rownames(summary(fitFDR20)$coefficients)))>0){
    terms20 <-rownames(summary(fitFDR20)$coefficients)[grep("data",rownames(summary(fitFDR20)$coefficients))]
    idx20 <- unlist(strsplit(terms20,",|[]]"))
    name20 <-colnames(data)[as.numeric(idx20[seq(2,length(idx20),3)])]
    #id.in20 <- gsub(pattern = "X(.*)",replacement = "\\1",x = name20)
    #mrks20<-c(mrks20,id.in20)
    #out20[[as.numeric(i)]]<- id.in20
  }
  return(list("out5"=name5,"out20"=name20))
}



BEFDR_nofix <- function ( maximal.lm, FDR.q, mfactor = 1) {
  compute.Lambda <- function(k, m, Q) {
    i <- c(1:k)
    return((1/(k + 1)) * sum(qnorm((Q/2) * (i/(m + 1 - i * 
                                                 (1 - Q))))^2))
  }
  get.model.size <- function(a.lm) {
    require(MASS)
    return(extractAIC(a.lm)[1] - 1)
  }
  require(MASS)
  #the.scope <- list(lower = minimal.lm, upper = maximal.lm)
  m <- mfactor * get.model.size(maximal.lm)
  new.model.size <- get.model.size(maximal.lm)
  for (i in 1:m) {
    old.model.size <- new.model.size
    Lambda <- compute.Lambda(k = old.model.size - 1, m, Q = FDR.q)
    ## no min model
    new.model <- stepAIC(maximal.lm, direction = "backward", k = Lambda, trace = FALSE)
    new.model.size <- get.model.size(new.model)
    if (new.model.size >= old.model.size) 
      break
  }
  new.lm <- lm(new.model)
  return(new.lm)
}

################################### With cov

BEFDR <- function (minimal.lm, maximal.lm, FDR.q, mfactor = 1) {
  compute.Lambda <- function(k, m, Q) {
    i <- c(1:k)
    return((1/(k + 1)) * sum(qnorm((Q/2) * (i/(m + 1 - i * 
                                                 (1 - Q))))^2))
  }
  get.model.size <- function(a.lm) {
    require(MASS)
    return(extractAIC(a.lm)[1] - 1)
  }
  require(MASS)
  the.scope <- list(lower = minimal.lm, upper = maximal.lm)
  m <- mfactor * get.model.size(maximal.lm)
  new.model.size <- get.model.size(maximal.lm)
  for (i in 1:m) {
    old.model.size <- new.model.size
    Lambda <- compute.Lambda(k = old.model.size - 1, m, Q = FDR.q)
    new.model <- stepAIC(maximal.lm, direction = "backward", 
                         scope = the.scope, k = Lambda, trace = FALSE)
    new.model.size <- get.model.size(new.model)
    if (new.model.size >= old.model.size) 
      break
  }
  new.lm <- lm(new.model)
  return(new.lm)
}

BE_analysis_cov<- function(cov,phe,data,fdr=c(0.05,0.2),mfactor=1){ 
  # data is a data frame contain cov genotype, cov is the name of covariants
  # build full model 
  nna <- complete.cases(cbind(phe,geno.mat))
  mrk <- colnames(data)
  id.full<- seq(1:length(mrk))
  geno_add <- paste("as.numeric(data[,",id.full,"])",sep="",collapse="+")
  phe<-phe[nna]
  data <- data[nna,]
  #test_fx1<-fx1[nna]
  fm <- as.formula(paste("phe ~ ",geno_add,sep=""))
  reg.lm <- lm(fm, y=TRUE)
  
  # min model
  id.min<-c(1:ncol(data))[mrk %in% cov]
  geno_bc.add<-paste("as.numeric(data[,",id.min,"])",sep="",collapse="+")
  fm.min<- as.formula(paste("phe ~",geno_bc.add,sep=""))
  min.lm <- lm(fm.min,y=TRUE)
  
  #Perform Backward-Elimination in the original data at different adaptive FDR thresholds
  fitFDR5<-BEFDR( minimal.lm = min.lm, maximal.lm = reg.lm, FDR.q = fdr[1],mfactor=mfactor)
  
  if(sum(grepl("data",rownames(summary(fitFDR5)$coefficients)))>0){
    terms5<-rownames(summary(fitFDR5)$coefficients)[grep("data",rownames(summary(fitFDR5)$coefficients))]
    idx5<-unlist(strsplit(terms5,",|[]]"))
    name5<-colnames(data)[as.numeric(idx5[seq(2,length(idx5),3)])]
    #id.in5<- name5[!(name5%in% sigMrk_sub)]
    #mrks5<-c(mrks5,id.in5)
    #out5[[as.numeric(i)]]<- id.in5
  }
  
  fitFDR20<-BEFDR( minimal.lm = min.lm, maximal.lm = reg.lm, FDR.q = fdr[2],mfactor=mfactor)
  
  if(sum(grepl("data",rownames(summary(fitFDR20)$coefficients)))>0){
    terms20 <-rownames(summary(fitFDR20)$coefficients)[grep("data",rownames(summary(fitFDR20)$coefficients))]
    idx20 <- unlist(strsplit(terms20,",|[]]"))
    name20 <-colnames(data)[as.numeric(idx20[seq(2,length(idx20),3)])]
    #id.in20 <- gsub(pattern = "X(.*)",replacement = "\\1",x = name20)
    #mrks20<-c(mrks20,id.in20)
    #out20[[as.numeric(i)]]<- id.in20
  }
  
  return(list("out5"=name5,"out20"=name20))
}
# function 
# selecting indpendent mrk
get.indpen.mrk <- function(r2.mat=r$LDmatrix,cut=0.9){
  diag(r2.mat) <-0
  r2.mat[lower.tri(r2.mat)] <- r2.mat[upper.tri(r2.mat)]
  
  all.marker <- colnames(r2.mat)
  all.chr <- gsub(pattern = ".*_(chr.*)_\\d+_.*",replacement = "\\1",x = all.marker)
  chr.level <- unique(all.chr)
  snp <-c()
  snp.rm <-c()
  for( chr in chr.level){
    all.marker.sub <- all.marker[all.chr %in% chr]
    r2.mat.sub <- r2.mat[all.marker.sub,all.marker.sub]
    # all comb
    if(length(all.marker.sub)>2){
      all.com <- combn(x = all.marker.sub,m = 2)
      yes <- rep(F,ncol(all.com))
      for( j in 1:ncol(all.com)){
        if(r2.mat.sub[all.com[1,j],all.com[2,j]] > cut) # add distance and save these two marker, which is saved and which is removed
          yes[j] <- T
      }
      mrk.rm <- unique(all.com[1,yes])
      snp.rm <- c(snp.rm,mrk.rm)
      snp <- c(snp,all.marker.sub[!(all.marker.sub %in% mrk.rm)])
    }
  }
  ### make sure removed is all right
  #r2.mat.rm <- r2.mat[snp.rm,]
  #if(any(r2.mat.rm > cut))
  return(snp)
}

remove_tight.loci <- function(SNPs=full_name_g,genable=data,cut_dis=10e3,cut_r2=0.9){
  if(!require(igraph))
    require(igraph)
  if(!require(GenABEL))
    require(GenABEL)
  # rule1 physical distance
  pos.all <- map(genable[,unique(SNPs)])
  names(pos.all) <- unique(SNPs)
  pos <- pos.all[SNPs]
  diffs <- abs(outer(pos, pos, FUN = "-")) #all pairwise differences in physical position
  diffs[lower.tri(diffs)] <- 1000000
  diag(diffs) <- 1000000
  #identical(rownames(diffs),rownames(r2))
  ## rule 2 linkage
  r2 <- r2fast(data = genable,snpsubset = SNPs)
  diag(r2) <- 0
  r2[lower.tri(r2)] <- 0
  ## rule3 chromsome
  chrs <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = SNPs)
  chrs.order <- match(chrs,unique(chrs))
  diffs.chr <- abs(outer(chrs.order, chrs.order, FUN = "-")) #all pairwise differences in index
  diffs.chr[lower.tri(diffs.chr)] <- 1
  diag(diffs.chr) <- 1
  
  clump <- which(diffs < cut_dis & r2 > cut_r2 & diffs.chr==0) #differences smaller than n
  if(length(clump) <=1)
    stop("nothing")
  ## diffs.col here is wrong should be row and should be clump %% nrow(diffs)
  ## diffs.row should be ceiling(clump/ncol(diffs))
  diffs.col <- ceiling(clump/ncol(diffs)) #differences smaller than n, columns in the distance matrix
  diffs.row <- clump - (diffs.col - 1)*nrow(diffs) #differences smaller than n, rows in the distance matrix
  #mycol <- clump %% nrow(diffs)
  #myrow <- ceiling(clump/ncol(diffs))
  
  #require(igraph)
  diffs.graph <- graph_from_data_frame(data.frame(diffs.col, diffs.row))
  diffs.clust <- clusters(diffs.graph)
  
  snps <- SNPs
  relation <- list()
  name <- c()
  for(i in 1:diffs.clust$no){
    nodesInCluster <- as.numeric(names(diffs.clust$membership[diffs.clust$membership == i]))
    #chrs[nodesInCluster] <- chrs[nodesInCluster[1]]
    relation[[i]] <- snps[nodesInCluster]
    snps[nodesInCluster] <- snps[nodesInCluster[1]]
    name <- c(name, snps[nodesInCluster[1]])
    #pos[nodesInCluster] <- pos[nodesInCluster[1]]
  }
  names(relation) <- name
  if(identical(snps,SNPs)){
    warning("no change ")
    return(list("snp"=snps,"re"=relation))
  }else{
    return(list("snp"=snps,"re"=relation))
  }
}

get.snpnames <- function(chr.loca,allsnp=snp.all){
  #chr.loca <- paste(data.sub$chr,data.sub$pos,sep="_")
  get.index <- function(x) return(allsnp[grep(pattern = paste0(x,"_"),x = allsnp)])
  return(c(unlist(lapply(chr.loca, FUN = get.index))))
}

find_best <- function(snp_now,genable,snp_list,r2_cut=0.9){
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
    return(rownames(r2_all)[which.max(r2_all[snp_now,])])
    
  }
}


cmp_snps <- function(add_snp,epi_snp,genable,r2_cut=0.9){
  #add_snp <- tmp.before$addSNP1[i]
  #epi_snp <-tmp.before$epiSNP1[i]
  add_chr <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = add_snp)
  epi_chr <-gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = epi_snp)
  r2 <- r2fast(data = genable,snpsubset = c(add_snp,epi_snp))[1,2]
  if(add_chr==epi_chr & r2 >r2_cut){
    return(T)
  }else{
    warning(" No hits with r2 above ",r2_cut,  " r2 is" ,r2, "or they are in different chr")
    return(F)
  }
}
# This is for 4 traits no use anymore
# fit.GbyE.allele <- function(trait= trait_picked,genable=data,snp_test=qtl_select){
#   if(!require(reshape2))
#     require(reshape2)
#   if(!require(GenABEL))
#     require(GenABEL)
#   warning("Phenotypes must have been scaled","\n")
#   trait2 <- paste0(gsub(pattern = "(.*)\\.mean","\\1",trait_picked),c(rep(".1",4),rep(".2",4)))
#   num_trait <- length(trait2)
#   geno <- as.double.gwaa.data(data[,snp_test])
#   ph <- phdata(genable)[,c("id",trait2)]
#   na <- complete.cases(cbind.data.frame(ph,geno))  
#   ph.v <- melt(data = ph,id.vars="id",measure.vars=2:9)
#   #update variable
#   ph.v$variable <- gsub(pattern = "(.*)\\.\\d{1}",replacement = "\\1",x = ph.v$variable )
#   geno.all <- do.call("rbind", replicate(num_trait, geno, simplify = FALSE))
#   if(!identical(rownames(geno.all),ph.v$id))
#     stop("id in the genotype matrix and phentoype matrix do not match")
#   
#   # build full model
#   id.full<- 1:ncol(geno.all)
#   geno_add <- paste("as.numeric(geno.all[nna,",id.full,"])",sep="",collapse="+")
#   nna <- complete.cases(cbind(ph.v$value,ph.v$variable,geno.all))
#   test_phe <- ph.v$value[nna]
#   test_fx1 <- ph.v$variable[nna]
#   fm <- as.formula(paste("test_phe ~ test_fx1 +",geno_add,sep=""))
#   lm <- lm(fm)
#   allele.all <- levels(as.factor(colnames(geno.all)))
#   num <- length(allele.all)
#   p<- numeric(num)
#   
#   for( i in 1:num){
#     fm2 <- as.formula(paste("`test_phe` ~ `test_fx1` +",geno_add,"+ `test_fx1`:as.numeric(geno.all[nna,",i,"])",sep=""))
#     #fm2 <- as.formula(paste("`test_phe` ~ `test_fx1` +",geno_add,sep=""))
#     #lm(fm2)
#     if(!require(lmtest))
#       require(lmtest)
#     lm.inter2 <- try(silent = T,lm(fm2))
#     if(!inherits(lm.inter2,"try-error")){
#       an <- lrtest(lm,lm.inter2)
#       p[i] <- an$`Pr(>Chisq)`[2]
#     }
#     cat(i,"\n")
#   }
#   pos <- as.numeric(gsub(pattern = "\\d+_chr.*_(\\d+)_.*_.*",replacement = "\\1",x = snp_test))
#   Inter <- ifelse(p<0.05/length(snp_test),"Yes","No")
#   out <- cbind.data.frame(colnames(geno),chromosome(data[,snp_test]),pos,format(p,scientific = T,digits = 3),Inter)
#   colnames(out) <- c("SNP_name","Chromosome","position","P_vlaue","GxE")
#   return(out)
#   #save(p,file="results/180727_lm_model_allele.RData") ## 89 allele
#   #load("results/180727_lm_model_allele.RData")
# }
est.genetic.effect <- function(trait= trait_picked,genable=data,snp_test=qtl_select){
  if(!require(reshape2))
    require(reshape2)
  warning("Phenotypes must have been scaled","\n")
  num_trait <- length(trait)
  geno <- as.double.gwaa.data(data[,snp_test])
  ph <- phdata(genable)[,c("id",trait_picked)]
  na <- complete.cases(cbind.data.frame(ph,geno))
  # ph.v <- melt(data = ph[na,],id.vars="id",measure.vars=2:5)
  # geno.all <- do.call("rbind", replicate(num_trait, geno[na,], simplify = FALSE))
  # if(!identical(rownames(geno.all),ph.v$id))
  #   stop("id in the gentype matrix and phentoype matrix do not match")
  # 
  # # only additive no interaction fitted
  # lm <- lm(ph.v)
  
  ## iteratively fit lm and extract the output
  out_put <- data.frame(array(NA,c(length(trait)*3,length(qtl_select)+1)))
  for( i in 1:length(trait)){
    out <- lm(ph[na,i+1]~geno[na,])
    out.s <- summary(out)
    index <- seq((3*i-2),3*i)
    out_put[index,] <- t(out.s$coefficients[,c(1,2,4)])
  }
  colnames(out_put) <- gsub(pattern = "geno\\[.*\\](.*)",replacement = "\\1",row.names(out.s$coefficients))
  return(out_put)
}

get_loci_add_epi <- function(tmp,trait){
  #tmp <- tmp.clumped;trait="IndolaceticAcid"
  index <- grepl(pattern = trait,tmp$trait)
  epi_loci <- unique(c(tmp[index,"epiSNP1"],tmp[index,"epiSNP2"]))
  add_loci <- unique(c(tmp[index,"addSNP1"],tmp[index,"addSNP2"]))
  return(list("epi"=epi_loci,"add"=add_loci))
}
Test_pairwise_epi_josh <- function(mrks,data=data,trait){
  if(!require(lmtest))
    require(lmtest)
  havePheno <- !is.na(trait)
  geno <- as.double.gwaa.data(data[havePheno, mrks])
  #geno[geno == 2] <- 1
  pheno <- trait[havePheno]
  if(!require(lmtest))
    require(lmtest)
  # build full model
  lm_full <- lm(pheno~ geno[,1]*geno[,2])
  lm_null <- lm(pheno~ geno[,1]+geno[,2])
  p <- lrtest(lm_full,lm_null)$`Pr(>Chisq)`
  return(p[2])
}
find_neighbours <- function(SNPs=full_name_g,index ,genable=data,n=50){
  if(!require(igraph))
    require(igraph)
  if(!require(GenABEL))
    require(GenABEL)
  ##
  #snps <- as.numeric(gsub(pattern = "(.*)_chr.*",replacement = "\\1",x=SNPs))
  #pos <- as.numeric(gsub(pattern = ".*_chr.*_(.*)_.*_.*",replacement = "\\1",x=SNPs))
  diffs <- abs(outer(index, index, FUN = "-")) #all pairwise differences in index
  diffs[lower.tri(diffs)] <- 1000
  diag(diffs) <- 1000
  clump <- which(diffs < n) #differences smaller than n
  if(length(clump) <=1)
    stop("nothing")
  ## diffs.col here is wrong should be row and should be clump %% nrow(diffs)
  ## diffs.row should be ceiling(clump/ncol(diffs))
  diffs.col <- ceiling(clump/ncol(diffs)) #differences smaller than n, columns in the distance matrix
  diffs.row <- clump - (diffs.col - 1)*nrow(diffs) #differences smaller than n, rows in the distance matrix
  #mycol <- clump %% nrow(diffs)
  #myrow <- ceiling(clump/ncol(diffs))
  
  #require(igraph)
  diffs.graph <- graph_from_data_frame(data.frame(diffs.col, diffs.row))
  diffs.clust <- clusters(diffs.graph)
  
  snps <- SNPs
  relation <- list()
  name <- c()
  for(i in 1:diffs.clust$no){
    nodesInCluster <- as.numeric(names(diffs.clust$membership[diffs.clust$membership == i]))
    #chrs[nodesInCluster] <- chrs[nodesInCluster[1]]
    relation[[i]] <- snps[nodesInCluster]
    snps[nodesInCluster] <- snps[nodesInCluster[1]]
    name <- c(name, snps[nodesInCluster[1]])
    #pos[nodesInCluster] <- pos[nodesInCluster[1]]
  }
  names(relation) <- name
  if(identical(snps,SNPs)){
    warning("no change ")
    return(list("snp"=snps,"re"=relation))
  }else{
    return(list("snp"=snps,"re"=relation))
  }
}

