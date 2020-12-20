## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp20023)

## -----------------------------------------------------------------------------
set.seed(0)
s <- 7; n<- 35; p <- 100; r <- 0.1; beta <- c(1,-0.5,0.7,-1.2,-0.9,3,0.55);
data <- simu_sirs_generation(s,n,p,r,beta)
X <- data[[1]]
y <- data[[2]]
a <- 0.5
beta_hat <- sirs(X,y,a)
sum(beta_hat!=0)
beta_hat[1:s]

## ----eval=FALSE---------------------------------------------------------------
#  as <- c(0,0.05,0.1,0.2,1,2,5)
#  rs <- c(0,0.2,0.5)
#  beta <- c(1,-0.5,0.7,-1.2,-0.9,3,0.55)
#  s <- 7
#  n <- 35
#  p <- 1000
#  p_as <- length(as)
#  p_rs <- length(rs)
#  n_data_set <- 100;
#  success_percent <- matrix(0,p_as,p_rs)
#  
#  for (i in 1:p_as){
#    for (j in 1:p_r){
#      a <- as[i]
#      r <- rs[j]
#      record <- 0
#      for (k in 1:n_data_set){
#        res <- simu_sirs_generation(s,n,p,r,beta)
#        X <- res[[1]]
#        y <- res[[2]]
#        beta_sirs <- sirs(X,y,a)
#        if (identical(beta_sirs[beta_sirs!=0],1:s)){
#          record <- record + 1
#        }
#      }
#      success_percent[i,j] <- record / n_data_set
#    }
#  }

## -----------------------------------------------------------------------------
data(PMdata_original)

## ----eval=FALSE---------------------------------------------------------------
#  K <- 100; b_threshold <- 3; Rmse_threshold <- 4;
#  l <- geo_anomaly(PMdata_original,K,b_threshold,Rmse_threshold)

## -----------------------------------------------------------------------------
library(ggplot2)
data(l)

## -----------------------------------------------------------------------------
PMdata_type <- cbind(PMdata_original,type = rep(1,nrow(PMdata_original)))
for (i in 1:nrow(PMdata_original)){
  if (PMdata_type[i,1] %in% l[[1]])
    PMdata_type$type[i] <- 2
  else{
    if (PMdata_type[i,1] %in% l[[2]])
      PMdata_type$type[i] <- 3
    else{
      if(PMdata_type[i,1] %in% l[[3]])
        PMdata_type$type[i] <- 4
    }
  }
}
PMdata_type$type <- as.factor(PMdata_type$type)

## -----------------------------------------------------------------------------
ggplot(PMdata_type,aes(x=lon,y=lat,color=type,shape=type))+
  scale_color_manual(values = c("green","black","orange","red"))+
  geom_point(size=1)+
  coord_equal(ratio=1)+
  labs(title="AirBox locations in Taiwan",x="longitude",y="latitude")+
  theme(plot.title = element_text(hjust=0.5))

## -----------------------------------------------------------------------------
ggplot(PMdata_type,aes(x=lon,y=lat))+
  geom_point(size=1)+
  facet_grid(.~type)+
  labs(title="AirBox locations in Taiwan",x="longitude",y="latitude")+
  theme(plot.title = element_text(hjust=0.5))

## -----------------------------------------------------------------------------
PMdata <- PMdata_original[,-1]
plot(c(1,57),c(0,480),type = "n",xlab="Time",ylab="PM2.5",main = "Normal AirBoxes")
ccolor <- c("green","blue","black","red")
for (i in (1:(nrow(PMdata)))){
  if (PMdata_type[i,]$type == 1)
    lines(1:(ncol(PMdata)-2),PMdata[i,-c(1,2)],col=ccolor[PMdata_type[i,]$type])
}

plot(c(1,57),c(0,480),type = "n",xlab="Time",ylab="PM2.5",main = "AirBoxes with relatively high PM2.5")
ccolor <- c("green","blue","black","red")
for (i in (1:(nrow(PMdata)))){
  if (PMdata_type[i,]$type == 2)
    lines(1:(ncol(PMdata)-2),PMdata[i,-c(1,2)],col=ccolor[PMdata_type[i,]$type])
}

plot(c(1,57),c(0,480),type = "n",xlab="Time",ylab="PM2.5",main = "AirBoxes with relatively low PM2.5")
ccolor <- c("green","blue","black","red")
for (i in (1:(nrow(PMdata)))){
  if (PMdata_type[i,]$type == 3)
    lines(1:(ncol(PMdata)-2),PMdata[i,-c(1,2)],col=ccolor[PMdata_type[i,]$type])
}

plot(c(1,57),c(0,480),type = "n",xlab="Time",ylab="PM2.5",main = "Malfunctioed AirBoxes")
ccolor <- c("green","blue","black","red")
for (i in (1:(nrow(PMdata)))){
  if (PMdata_type[i,]$type == 4)
    lines(1:(ncol(PMdata)-2),PMdata[i,-c(1,2)],col=ccolor[PMdata_type[i,]$type])
}

