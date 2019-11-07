# Prep data
require(igraph)
require(RColorBrewer)
#dev.off()
#display.brewer.all()
##################
##################
##define color
##################
##################

c1 <- brewer.pal(9,"Blues")[3]
c2 <-  brewer.pal(9,"YlOrRd")[3]
c3 <-  brewer.pal(9,"Reds")[4]

##################
##################
##simulate  G-P map
##################
##################
g1 <- sample(x = c(0,1),size = 200,replace = T,prob = c(0.5,0.5))# locus A
g2 <- sample(x = c(0,1),size = 200,replace = T,prob = c(0.5,0.5)) # locus B
g3 <- sample(x = c(0,1),size = 200,replace = T,prob = c(0.5,0.5)) #
g4 <- sample(x = c(0,1),size = 200,replace = T,prob = c(0.5,0.5)) # locus D

y1 <- g1*g2+g1*g3+ g1*g4*0.05+rnorm(200,mean = 0,sd = 0.3)
y2 <- g2*1+g1*g4+rnorm(200,mean = 0,sd = 0.3)


##################
##################
##simulate G-G Network
##################
##################

e1 <- data.frame("locus1"=c("A","B"),"locus"=c("C","A"))
e1_graph <-graph_from_data_frame(e1, directed = F, vertices = NULL)
V(e1_graph)$color <- c1

e2 <- data.frame("locus1"=c("A","B"),"locus"=c("D","B"))
e2_graph <-graph_from_data_frame(e2, directed = F, vertices = NULL)
V(e2_graph)$color <- c2

e3 <- rbind(e1,e2)
e3_graph <-graph_from_data_frame(e3, directed = F, vertices = NULL)
V1 <- names(V(e1_graph))
V2 <- names(V(e2_graph))
inter <- intersect(V1,V2)
V(e3_graph)$color <- ifelse(names(V(e3_graph)) %in% as.character(inter),c3,ifelse(names(V(e3_graph)) %in% V1,c1,c2))
a <- layout.auto(e3_graph)
a_index12 <- names(V(e3_graph))
a_index1 <- names(V(e1_graph))
a_index2<- names(V(e2_graph))
#dev.off()
#par(mfrow=c(3,1))
plot(e1_graph)
plot(e2_graph)
plot(e3_graph)
legend("topleft",legend = c("Shared","Unique to Env1","Unique to Env2"),col = c(c3,c1,c2),fill = c(c3,c1,c2),border =c(c3,c1,c2),bty = "n" )

##################
##################
##plot
##################
##################
par(mfrow=c(3,3),mar=c(4,3,2,2))
x = 0.2
boxplot(y1~g1*g2,frame.plot=F,col=c1,las=1,names=c("AABB","AAbb","aaBB","aabb"),boxwex=x)

plot(0,col="white",frame.plot = F,xaxt="n",yaxt="n",xlab = "",ylab = "")
boxplot(y2~g1*g4,frame.plot=F,col=c2,las=1,names=c("AADD","AAdd","aaDD","aadd"),boxwex=x)

#plot(0,col="white",frame.plot = F,xaxt="n",yaxt="n")

boxplot(y1~g1*g3,frame.plot=F,col=c1,las=1,names=c("AACC","AAcc","aaCC","aacc"),boxwex=x)

plot(0,col="white",frame.plot = F,xaxt="n",yaxt="n",xlab = "",ylab = "")

boxplot(y2~g2,frame.plot=F,col=c2,las=1,names=c("BB","bb"),boxwex=0.4*x)
my_box <- function(y,col){
  #boxplot(y1 ~ g2*g1*g3*g4,frame.plot=F,col=c1,las=2,plot=T,names=nm)
  #nm <- paste0(rep(c("AA","aa"),8),rep(c("BB","BB","bb","bb"),4),rep(c(rep("CC",4),rep("cc",4)),2),c(rep("DD",8),rep("dd",8)))
  nm <- paste0(rep(c("A","a"),8),rep(c("B","B","b","b"),4),rep(c(rep("C",4),rep("c",4)),2),c(rep("D",8),rep("d",8)))
  
  b<-boxplot(y ~ g1*g2*g3*g4,frame.plot=F,col=c1,las=2,plot=F,names=nm) #boxplot(y1 ~ g2*g1*g3*g4,frame.plot=F,col=c1,las=2,plot=T,names=c("AABBCC","aaBBCC","AAbbCC","aabbCC","AABBcc","aaBBcc","AAbbcc","aabbcc"))
  b <- b
  or <- order(b$names)#order(b$stats[3,])
  b$stats <- b$stats[,or]
  b$n <- b$n[or]
  for(i in 1:length(or)){
    b$group[b$group == or[i]] <- i
  }
  b$conf <- b$conf[or]
  b$names <-b$names[or]
  bxp(b,frame.plot=F,las=2,boxfill=col)
}

## joint network
par(mar=c(6,4,0,2))
my_box(y = y1,col=c(rep(c3,8),rep(c1,8)))
plot(0,col="white",frame.plot = F,xaxt="n",yaxt="n",xlab = "",ylab = "")

my_box(y = y2,col=c(rep(c3,8),rep(c2,8)))


# b<- boxplot(y1 ~ g2*g1*g3,frame.plot=F,col=c1,las=2,plot=T,names=c("AABBCC","aaBBCC","AAbbCC","aabbCC","AABBcc","aaBBcc","AAbbcc","aabbcc"))
# b <- b
# or <- order(b$stats[3,])
# b$stats <- b$stats[,or]
# b$n <- b$n[or]
# for(i in 1:length(or)){
#   b$group[b$group == or[i]] <- i
# }
# b$conf <- b$conf[or]
# b$names <-b$names[or]
# bxp(b,frame.plot=F,col=c1,las=2,boxfill=c1)

############### Second enviroment
# 
# #c(rnorm(100,mean = 2,sd = 0.3),rnorm(100,mean=3,sd=0.35))
# geno <- rep(c("AABB","AAbb","aaBB","aabb"),each=50)
# boxplot(y~geno,frame.plot=F,col=c1,las=2)
# 
# y <- c(rnorm(150,mean = 2,sd = 0.3),rnorm(50,mean=3,sd=0.35))
# geno <- rep(c("AACC","AAcc","aaCC","aacc"),each=50)




## Pannel A