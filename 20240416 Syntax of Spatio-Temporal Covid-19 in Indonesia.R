#Article: "Spatio-temporal Clustering using Generalized Lasso to Identify the Spread of Covid-19 in Indonesia according to Provincial Flight Route-based Connections"
#Authors: Septian RAHARDIANTORO*, Sachnaz Desta OKTARINA, Anang KURNIA, Nickyta Shavira MAHARANI, & Alfidhia Rahman Nasa JUHANDA.
#Corresponding Email*: septianrahardiantoro@apps.ipb.ac.id 
#Last Modified: 2024-04-16

#Spatio-temporal for Covid-19 Indonesia
######Data Preparation
dt <- read.csv("20240416 Weekly Covid-19 Cases in Indonesia.csv",sep=";",header=T)
dim(dt) #34 135

##Weekly Covid-19 Cases
dt1 <- dt[,2:134]
dim(dt1) #34 133

##Response variable y_it = ln(Weekly Covid-19 Cases*100000/Population)
dt3 <- log((dt1+1)*100000/dt$pop)

######Trend-Filtering
##Creating D for temporal effect (D2 for trend filtering)
row <- ncol(dt3)-2
D1 <- diag(rep(1,row))
dim(D1) #131 131
D2 <- (-2)*D1
D3 <- D1
nol <- matrix(0,row,1)
D11 <- cbind(D1,nol,nol)
D21 <- cbind(nol,D2,nol)
D31 <- cbind(nol,nol,D3)
D.period <- D11+D21+D31 #D2
dim(D.period) #131 133

####Temporal effect detection
##Response variable preparation
y2.tbar <- apply(dt3,2,mean)
y.tbar <- y2.tbar

##Predictor matrix preparation
Xa <- diag(rep(1,length(y2.tbar)))

##Generalized Lasso: Temporal Effects Detection
library(genlasso);library(sp);library(spdep);library(igraph)
library(rgdal);library(pracma)

##Setting maximum degree of freedom
df1.max <- ncol(dt3)*3/4

##Trend filtering
beta.tbar0 <- 0

zt <- y.tbar-beta.tbar0
mod1 <- genlasso(zt,Xa,D=D.period)

### Start of selecting lambda_T ###
df.T <- mod1$df[mod1$df<=df1.max]
lam1 <- mod1$lambda[mod1$df<=df1.max]
lam1 <- exp(seq(log(min(lam1)), log(max(lam1)), length=101))  
alpha <- coef(mod1, lambda=lam1, type="both") 
beta1 <- alpha$beta 
u1 <- alpha$u
df1 <- alpha$df

####ALOCV & GCV
alot <- c()
gt <- c()
for(j in 1:length(lam1)){
  E <- which(abs( abs(u1[,j])-lam1[j] )<0.00001)
  #print(E)
  if (length(E) == 0) {Dbaru <- D.period
  } else {Dbaru <- D.period[-E,]}
  if (length(E) == nrow(D.period)-1) Dbaru <- matrix(Dbaru,1,ncol(D.period)) 
  A <- Xa %*% nullspace(Dbaru)
  Aplus <- MASS::ginv(A)
  H <- A%*%Aplus
  hj <- diag(H)
  trH <- sum(hj)/ncol(dt3)
  hj <- ifelse(hj >=0.999,0.999,hj)
  alo <- mean(((zt-beta1[,j])/(1-hj))^2)
  alot <- c(alot,alo)
  go <- mean(((zt-beta1[,j])/(1-trH))^2)
  gt <- c(gt,go)
}
#ALOCV
best.lam1 <- lam1[which.min(alot)]
best.df1 <- df1[which.min(alot)]
#GCV
best.lam2 <- lam1[which.min(gt)]
best.df2 <- df1[which.min(gt)]

##Table 1. Table summary ALOCV & GCV for trend filtering
Method <- c("GenLasso with ALOCV","GenLasso with GCV")
Lambda_T <- c(best.lam1,best.lam2)
df <- c(best.df1,best.df2)
table1 <- data.frame(Method,Lambda_T,df)
table1
#               Method   Lambda_T df
#1 GenLasso with ALOCV 0.60756068 43
#2   GenLasso with GCV 0.01466067 97

##Temporal effects estimates: aplha_t hat
#ALOCV  
alpha.hat <- coef(mod1, lambda=best.lam1)$beta
#GCV
alpha.hat2 <- coef(mod1, lambda=best.lam2)$beta

##Figure trend filtering estimation
#Figure 5
Methods <- rep(c("Unpenalized","GenLasso with ALOCV","GenLasso with GCV"),each=133)
Weeks <- rep(1:133,3)
est <- c(zt,alpha.hat,alpha.hat2)
dtemp <- data.frame(Methods,Weeks,est)
library(ggplot2)
cols <- c("#D43F3A", "#EEA236", "#5CB85C")
#install.packages("hrbrthemes")
library(hrbrthemes)

ggplot(dtemp, aes(x = Weeks, y = est, color = Methods)) +
  geom_line(linewidth=0.7,alpha=0.7) +
  theme_ipsum()+
  ylab("Target Variable")

######Spatial effect detection
####D Penalty Matrices Preparations
####Function for convertion M (neighboorhood) to D (penalty matrix)
matriksD <- function(w3) {
  nk <- ncol(w3)
  indc <- which(w3 == 1,arr.ind = TRUE)
  indc.sort <- indc[order(indc[,1],decreasing=FALSE),]
  
  for(j in 1:nrow(indc.sort)){
    if(indc.sort[j,1]>indc.sort[j,2]) {
      a <- indc.sort[j,1]
      indc.sort[j,1] <- indc.sort[j,2]
      indc.sort[j,2] <- a
    }
  }
  indc.sort
  
  indc.sort.uniq <- indc.sort[!duplicated(indc.sort), ]
  
  D2 <- NULL
  for(i in 1:nrow(indc.sort.uniq)){
    dd <- rep(0,nk)
    dd[indc.sort.uniq[i,1]] <- c(-1)
    dd[indc.sort.uniq[i,2]] <- c(1)
    D2 <- rbind(D2,dd)
  }
  return(D2)
}

##D based on contiguity (border based neighborhood)
m1 <- read.csv("20240416 provinces neighborhood.csv",sep=";",header=T)
dim(m1)
mt1 <- m1[,-1]
dim(mt1)
D1 <- matriksD(mt1) #D based on contiguity
dim(D1) #32 34

##D based on flight connection
m2 <- read.csv("20240416 flight connection.csv",sep=";",header=T)
dim(m2)
mt2 <- m2[,-1]
dim(mt2)
D2 <- matriksD(mt2) #D based on flight connection
dim(D2) #276  34

#####Generalized Lasso: Spatial Effects Detection
library(genlasso);library(sp);library(spdep);library(igraph)
library(rgdal);library(pracma)

##Setting maximum degree of freedom
df2.max <- nrow(dt3)*3/4 #for spatial effects detection

#####Fused Lasso for a Common lambda_S
##Variable preparation
y2.tbar <- apply(dt3,2,mean)
y.tbar <- y2.tbar
m.tbar <- matrix(rep(y.tbar,each=nrow(dt3)),nrow(dt3),ncol(dt3))

rit <- dt3-m.tbar
l.time <- ncol(dt3)

Xb1 <- diag(rep(1,nrow(dt3)))

##Checking lambda_S and DF for each t
####D based on contiguity (border based neighborhood)
c.lam <- c()
c.df <- c()
for(z in 1:l.time){
  md <- genlasso(rit[,z],Xb1,D=D1)
  mdl <- md$lambda[md$df<=df2.max]
  mdd <- md$df[md$df<=df2.max]
  c.lam <- c(c.lam,mdl)
  c.df	<- c(c.df,mdd)
}
range(c.lam)
#[1] 0.1764182 3.9841617
range(c.df)
#[1] 12 25

##Start of selecting lambda_S using ALOCV & GCV
lam2 <- seq(min(c.lam),max(c.lam),length=101)

beta.k.all <- NULL
alot2.all <- NULL
gt2.all <- NULL
df2 <- NULL

for(t in 1:l.time){
  mod2 <- genlasso(rit[,t],Xb1,D=D1)
  beta <- coef(mod2, lambda=lam2, type="both") 
  beta2 <- beta$beta
  u2 <- beta$u
  df2 <- rbind(df2,beta$df)
  
  alot2 <- c()
  gt2 <- c()
  
  for(k in 1:length(lam2)){
    E2 <- which(abs( abs(u2[,k])-lam2[k] )<0.00001)
    if (length(E2) == 0) {Dbaru2 <- D1
    } else {Dbaru2 <- D1[-E2,]}
    if (length(E2) == nrow(D1)-1) Dbaru2 <- matrix(Dbaru2,1,ncol(D1)) 
    A2 <- Xb1 %*% nullspace(Dbaru2)
    Aplus2 <- MASS::ginv(A2)
    H2 <- A2%*%Aplus2
    hj2 <- diag(H2)
    trH2 <- sum(hj2)/nrow(dt3)
    hj2 <- ifelse(hj2 >=0.999,0.999,hj2)
    alo2 <- mean(((rit[,t]-beta2[,k])/(1-hj2))^2)
    alot2 <- c(alot2,alo2)
    go2 <- mean(((rit[,t]-beta2[,k])/(1-trH2))^2)
    gt2 <- c(gt2,go2)
  }
  beta.k.all <- rbind(beta.k.all,beta2)
  alot2.all <- rbind(alot2.all,alot2)
  gt2.all <- rbind(gt2.all,gt2)
}
dfS <- apply(df2,2,mean)

alocv.error <- apply(alot2.all,2,mean)
ind.lamb <- which.min(alocv.error)
best.lamS <- lam2[ind.lamb]
beta.hatS <- beta.k.all[,ind.lamb]
best.dfS <- dfS[ind.lamb]

gcv.error <- apply(gt2.all,2,mean)
ind.lamb2 <- which.min(gcv.error)
best.lamS2 <- lam2[ind.lamb2]
beta.hatS2 <- beta.k.all[,ind.lamb2]
best.dfS2 <- dfS[ind.lamb2]
## End of selecting lambda_S ###

##Table 2
Method <- c("GenLasso with ALOCV","GenLasso with GCV")
Lambda_S <- c(best.lamS,best.lamS2)
df_Spatial <- c(best.dfS,best.dfS2)
table2 <- data.frame(Method,Lambda_S,df_Spatial)
table2
#               Method Lambda_S df_Spatial
#1 GenLasso with ALOCV 3.984162         12
#2   GenLasso with GCV 3.984162         12

##plotting heatmap
beta.hat.D1 <- matrix(beta.hatS2,nrow(dt3),ncol(dt3))
prov.names <- dt$PROVINSI
prov.names[5] <- "DI YOGYAKARTA"
rownames(beta.hat.D1) <- prov.names
dd<-seq(as.Date("2020/3/07"), as.Date("2022/9/17"),length.out=133)
colnames(beta.hat.D1)<- dd

##Figure 7
library(RColorBrewer);library(gplots);library(heatmaply)
mypalette<-brewer.pal(9,"OrRd")
pheat.D1 <- heatmaply(beta.hat.D1,dendrogram = "none",colors=mypalette,fontsize_row = 7, fontsize_col = 6.5,plot_method="ggplot",
                      labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat.D1


#############
##Checking lambda_S and DF for each t
####D based on flight connection
c.lam <- c()
c.df <- c()
for(z in 1:l.time){
  md <- genlasso(rit[,z],Xb1,D=D2)
  mdl <- md$lambda[md$df<=df2.max]
  mdd <- md$df[md$df<=df2.max]
  c.lam <- c(c.lam,mdl)
  c.df	<- c(c.df,mdd)
}
range(c.lam)
#[1] 0.01892362 0.67903089
range(c.df)
#[1] 0 25

##Start of selecting lambda_S using ALOCV & GCV
lam2 <- seq(min(c.lam),max(c.lam),length=101)

beta.k.all <- NULL
alot2.all <- NULL
gt2.all <- NULL
df2 <- NULL

for(t in 1:l.time){
  mod2 <- genlasso(rit[,t],Xb1,D=D2)
  beta <- coef(mod2, lambda=lam2, type="both") 
  beta2 <- beta$beta
  u2 <- beta$u
  df2 <- rbind(df2,beta$df)
  
  alot2 <- c()
  gt2 <- c()
  
  for(k in 1:length(lam2)){
    E2 <- which(abs( abs(u2[,k])-lam2[k] )<0.00001)
    if (length(E2) == 0) {Dbaru2 <- D2
    } else {Dbaru2 <- D2[-E2,]}
    if (length(E2) == nrow(D2)-1) Dbaru2 <- matrix(Dbaru2,1,ncol(D2)) 
    if(is.null(nullspace(Dbaru2))) {nsm <- Xb1
    } else {nsm <- nullspace(Dbaru2)}
    A2 <- Xb1 %*% nsm
    Aplus2 <- MASS::ginv(A2)
    H2 <- A2%*%Aplus2
    hj2 <- diag(H2)
    trH2 <- sum(hj2)/nrow(dt3)
    hj2 <- ifelse(hj2 >=0.999,0.999,hj2)
    alo2 <- mean(((rit[,t]-beta2[,k])/(1-hj2))^2)
    alot2 <- c(alot2,alo2)
    go2 <- mean(((rit[,t]-beta2[,k])/(1-trH2))^2)
    gt2 <- c(gt2,go2)
  }
  beta.k.all <- rbind(beta.k.all,beta2)
  alot2.all <- rbind(alot2.all,alot2)
  gt2.all <- rbind(gt2.all,gt2)
}
dfS <- apply(df2,2,mean)

alocv.error.2 <- apply(alot2.all,2,mean)
ind.lamb.2 <- which.min(alocv.error.2)
best.lamS.2 <- lam2[ind.lamb.2]
beta.hatS.2 <- beta.k.all[,ind.lamb.2]
best.dfS.2 <- dfS[ind.lamb.2]

gcv.error.2 <- apply(gt2.all,2,mean)
ind.lamb2.2 <- which.min(gcv.error.2)
best.lamS2.2 <- lam2[ind.lamb2.2]
beta.hatS2.2 <- beta.k.all[,ind.lamb2.2]
best.dfS2.2 <- dfS[ind.lamb2.2]
## End of selecting lambda_S ###

##Table 3
Method <- c("GenLasso with ALOCV","GenLasso with GCV")
Lambda_S2 <- c(best.lamS.2,best.lamS2.2)
df_Spatial2 <- c(best.dfS.2,best.dfS2.2)
table3 <- data.frame(Method,Lambda_S2,df_Spatial2)
table3
#               Method Lambda_S df_Spatial
#1 GenLasso with ALOCV 0.01892362    31.79699
#2   GenLasso with GCV 0.11133864    13.08271


###########################
#####Figure heatmap
##Figure 6 Unpenalized beta
beta.hat.unp <- rit
prov.names <- dt$PROVINSI
prov.names[5] <- "DI YOGYAKARTA"
rownames(beta.hat.unp) <- prov.names
dd<-seq(as.Date("2020/3/07"), as.Date("2022/9/17"),length.out=133)
colnames(beta.hat.unp)<- dd

library(RColorBrewer);library(gplots);library(heatmaply)
mypalette<-brewer.pal(9,"OrRd")
pheat.unp <- heatmaply(beta.hat.unp,dendrogram = "none",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                       labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat.unp 

##Figure 8 ALOCV
beta.hat.D2alo <- matrix(beta.hatS.2,nrow(dt3),ncol(dt3))
rownames(beta.hat.D2alo) <- prov.names
colnames(beta.hat.D2alo)<- dd

pheat.D2alo <- heatmaply(beta.hat.D2alo,dendrogram = "none",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                         labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat.D2alo

##Figure 9 GCV
beta.hat.D2g <- matrix(beta.hatS2.2,nrow(dt3),ncol(dt3))
rownames(beta.hat.D2g) <- prov.names
colnames(beta.hat.D2g)<- dd

pheat.D2g <- heatmaply(beta.hat.D2g,dendrogram = "none",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                       labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat.D2g

##Figure 10: selected GCV result
pheat4 <- heatmaply(beta.hat.D2g,dendrogram = "row",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                    labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat4

#######################################################################################
###End of analysis



