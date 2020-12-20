## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp20023)

## ---- echo=FALSE, out.width="50%", fig.cap="I would like some French fries"----
knitr::include_graphics("example1.jpeg")

## -----------------------------------------------------------------------------
knitr::kable(head(iris), "simple")

## -----------------------------------------------------------------------------
library(ggplot2)
library(gridExtra)

## -----------------------------------------------------------------------------
pareto_pdf <- function(x, lambda = 1, k = 1){
    density <- (k*(lambda^k)) / (x^(k + 1))
    return(density)
}
set.seed(33)
a <- 2
b <- 2
n <- 1000
u <- runif(n)
x <- b / (1-u)^{1/a}
P1 <- ggplot(data = data.frame(x),aes(x)) +
  geom_histogram(aes(y = ..density.. ),binwidth = 0.5,boundary = 2)+
  labs(title = "Density Histogram of Samples ", x = "X", Y ="Density")
P2 <- ggplot(data = data.frame(x[x<15]),aes(x[x<15])) +
  xlim(2,15) +
  geom_histogram(aes(y = ..density.. ),binwidth = 0.5,boundary = 2,color="#e9ecef")+
  stat_function(fun = pareto_pdf, args = list(lambda = 2, k =2), colour = "deeppink")+
    annotate("text", x = 4, y = 0.8, parse = TRUE, size = 3,
           label = "f(x) == frac(8, x^3)")+
  labs(title = "Density Histogram of Samples [<15]", x = "X", Y ="Density")
grid.arrange(P1,P2,ncol=2)

## -----------------------------------------------------------------------------
Ep_pdf <- function(x)(0.75*(1-x^2))
set.seed(339)
n <- 100000
EP_generator <- function(n){#function that generates random samples
  u <- runif(n * 3)
  U <- matrix(u,ncol = 3)
  x <- U[,3]*(1-I(U[,2]<= U[,3])*I(U[,1] <= U[,3])) + U[,2]*I(U[,2]<= U[,3])*I(U[,1] <= U[,3])
  x <- x*sample(c(-1,1),prob= c(0.5,0.5),size = n, replace = TRUE)
  return(x)
}
x <- EP_generator(n)
ggplot(data = data.frame(as.numeric(x)),aes(as.numeric(x))) +
  xlim(-1,1) +
  geom_histogram(aes(y = ..density.. ),binwidth = 0.05,boundary = -1,fill="#69b3a2", 
                 color="#e9ecef", alpha=0.9)+
  stat_function(fun=Ep_pdf,colour="deeppink")+annotate("text", x = 0.75, y = 0.6, parse = TRUE, size = 3,
           label = "f(x) == frac(3, 4) (1-x^2)")+
  labs(title ="Density Histogram of the Samples",x="X")

## -----------------------------------------------------------------------------
#Create the Pareto density function
pareto_pdf2 <- function(x,beta = 2,r = 4){
  density <- beta^r * r / (beta+x)^(r+1)
}
#Simulate random obserbations from the mixture
beta <- 2
r <- 4
n <- 1000
set.seed(313)
lambda <- rgamma(n,r,beta)
x <- rexp(n,lambda)
ggplot(data = data.frame(x),aes(x)) +
  xlim(0,15) +
  geom_histogram(aes(y = ..density..),binwidth = 0.2,boundary = 0,color="#e9ecef")+annotate("text", x = 2, y = 1.5, parse = TRUE, size = 3,
           label = "f(x) == frac(64, {(2+x)}^5)")+
  stat_function(fun = pareto_pdf2,args = list(beta = 2, r=4),colour = "deeppink")+
  labs(title = "Density Hsitogram of the Samples",x = "X")

## -----------------------------------------------------------------------------
set.seed(51)
num <- 1e4
x_simu <- runif(num, min=0, max=pi/3)
results <- sin(x_simu) * pi / 3
result_hat <- mean(sin(x_simu))*pi/3
se <- sd(results) / sqrt(num)
print(c(result_hat, 0.5)) #print my estimate along with the real value of the intergral
round(c(result_hat - 1.96*se,result_hat + 1.96*se),4) #print the 95% CI

## -----------------------------------------------------------------------------
MC.anti <- function(iter = 1e4, anti = TRUE){
  u <- runif(iter/2)
  if(!anti) v <- runif(iter/2) else v <- 1-u
  u <- c(u,v)
  return(mean(exp(u)))
}#The function that generates random samples with simple MC or antithetic variable approach

## -----------------------------------------------------------------------------
iter1 <- 1e4
iter2 <- 1e4
MC_simple <- sapply(rep(iter2,iter1), FUN = MC.anti, anti = FALSE) #Simple Monte Carlo
MC_anti <- sapply(rep(iter2,iter1), FUN = MC.anti, anti = TRUE) #Antithetic approach
print(c(mean(MC_simple),mean(MC_anti))) #print estimates by simple Monte Carlo and Antithetic approach
print((var(MC_simple) - var(MC_anti))/var(MC_simple)) #print the estimate of percent reduction in variance

## -----------------------------------------------------------------------------
library(kableExtra)

## -----------------------------------------------------------------------------
g <- function(x) x^2 * exp(-x^2 / 2) / sqrt(2*pi)
f1 <- function(x) 2 * exp(-(x-1)^2 / 2) / sqrt(2*pi)
f2 <- function(x) exp(-(x-1))
f1_ratio <- function(x) f1(x) / g(x)
f2_ratio <- function(x) f2(x) / g(x)

## -----------------------------------------------------------------------------
p1 <- ggplot(data.frame(x=c(1,5)), aes(x=x)) + #plot g(x) as well as two importance functions
  geom_line(stat='function',fun = g,aes(color = 'g')) +
  geom_line(stat = 'function',fun = f1,aes(color = 'f1')) +
  geom_line(stat='function',fun = f2,aes(color = 'f2')) +
  scale_colour_manual("Function",values = c('g' = 'black','f1'='red','f2' = 'blue'))
p2 <- ggplot(data.frame(x=c(1,3.5)), aes(x=x)) + #plot the ratio of two importance functions versus g(x)
  geom_line(stat='function',fun = f1_ratio,aes(colour = 'f1/g')) +
  geom_line(stat='function',fun = f2_ratio,aes(colour = 'f2/g')) + 
  scale_colour_manual('Ratio',values = c('f1/g'='red','f2/g'='blue'))
grid.arrange(p1,p2,ncol = 2,top = 'Figure1: Importance functions and their ratios')

## -----------------------------------------------------------------------------
m <- 1e4
theta.hat <- se <- numeric(2)
x <- abs(rnorm(m)) + 1 #using f1
fg <- g(x) / f1(x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
x <- rexp(m) + 1 #using f2
fg <- g(x) / f2(x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)
res <- rbind(theta = round(theta.hat,4),se = round(se,4))
colnames(res) <- paste("f",1:2)
kbl(res) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
#Divide [0,1] into 5 intervals
F_inv <- function(j,k=5){ #j th end point in case where we divide [0,1] into k intervals
  -log(1-(1-exp(-1))*j/k)
}
end_points <- sapply(0:5,FUN = F_inv)
g <- function(x) exp(-x) / (1+x^2)
f <- function(x) 5*exp(-x) / (1-exp(-1))
Fj_inv <- function(x,j) -log(exp(-end_points[j]) - (1-exp(-1))*x/5)

## -----------------------------------------------------------------------------
set.seed(515)
m <- 1e4 #number of replicates
k <- 5 #number of strata
r <- m/k #replicates per stratum
T <- matrix(0,r,k)
for (j in 1:k){
  x <- Fj_inv(runif(r),j)
  T[,j] <- g(x) / f(x)
}

## -----------------------------------------------------------------------------
sum(apply(T,2,mean))

## -----------------------------------------------------------------------------
sqrt(k*sum(apply(T,2,var)))

## -----------------------------------------------------------------------------
set.seed(64)
alpha = 0.05
mu = 0
sigma = 1
m <- 1e6
n <- 20
CI <- matrix(0,m,2)
for (i in 1:m){
  x <- rlnorm(n,meanlog = mu,sdlog = sigma)
  y <- log(x)
  z <- sd(y)*qt(1-alpha/2, df= n-1)/sqrt(n)
  CI[i,] <- c(mean(y) - z, mean(y) + z)
}
cl <- mean((CI[,1]<= mu)*(CI[,2] >= mu))
print(cl)

## -----------------------------------------------------------------------------
set.seed(65)
alpha <- 0.05
m <- 1e6
n <- 20
mu <- 2
CI <- matrix(0,m,2)
for (i in 1:m){
  x <- rchisq(n,df = 2)
  z <- sd(x) * qt(1-alpha/2,df = n-1) / sqrt(n)
  CI[i,] <- c(mean(x)-z, mean(x)+z)
}
cl <- mean((CI[,1]<= mu)*(CI[,2] >= mu))
print(cl)

## -----------------------------------------------------------------------------
library(MASS)

## -----------------------------------------------------------------------------
alpha_test <- 0.05
n <- 30
m <- 2000
cv <- qnorm(1-alpha_test/2,0,sqrt(6*(n-2)/((n+1)*(n+3))))
sk <- function(x){
  #compute the sample skewness coefficient
  x_bar <- mean(x)
  m3 <- mean((x-x_bar)^3)
  m2 <- mean((x-x_bar)^2)
  return(m3/m2^1.5)
}

## -----------------------------------------------------------------------------
set.seed(671)
#against symmetric beta distribution
alpha_beta <- beta_beta <- seq(2,100,2)
N <- length(alpha_beta)
pwr_beta <- numeric(N)
for (k in 1:N){
  sktests_beta <- numeric(m)
  for (i in 1:m){
    x <- rbeta(n,alpha_beta[k],beta_beta[k])
    sktests_beta[i] <- as.integer(abs(sk(x))>= cv)
  }
  pwr_beta[k] <- mean(sktests_beta)
}

## -----------------------------------------------------------------------------
data <- data.frame(alpha_beta,pwr_beta)
colnames(data) <- c("alpha","power")
ggplot(data,aes(x = alpha,y=power)) + geom_point() + ggtitle("power v.s. alpha in beta distribution")

## ----echo = FALSE-------------------------------------------------------------
ggplot(data.frame(x =c(0,1)),aes(x = x))+
  geom_line(stat="function",fun = dbeta,args = list(shape1 = 2,shape2 = 2),aes(color="alpha = 2"))+
  geom_line(stat="function",fun = dbeta,args = list(shape1 = 10,shape2 = 10),aes(color="alpha = 10"))+
    geom_line(stat="function",fun = dbeta,args = list(shape1 = 20,shape2 = 20),aes(color="alpha = 20"))+
  scale_colour_manual("Function",values = c('alpha = 2' = 'green','alpha = 10'='red','alpha = 20' = 'blue'))+
  ggtitle("Symmetric Beta probability density distribution function")

## -----------------------------------------------------------------------------
set.seed(67)
#against heavy-tailed t distribution
df <- seq(1,50,1)
N <- length(df)
pwr_t <- numeric(N)
for (j in 1:N){
  sktests_t <- numeric(m)
  for (i in 1:m){
    x <- rt(n,df[j])
    sktests_t[i] <- as.integer(abs(sk(x))>= cv) 
  }
  pwr_t[j] <- mean(sktests_t)
}

## -----------------------------------------------------------------------------
data <- data.frame(df,pwr_t)
colnames(data) <- c("df","power")
ggplot(data,aes(x = df,y=power)) + geom_point() + ggtitle("power v.s. df in t distribution")

## -----------------------------------------------------------------------------
count5test <- function(x,y){
  X <- x - mean(x)
  Y <- y - mean(y)
  out_x <- sum(X > max(Y)) + sum(X < min(Y))
  out_y <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(out_x,out_y))>5))
}
F_test <- function(x,y,alpha){
  X <- x - mean(x)
  Y <- y - mean(y)
  F_res <- var.test(X,Y,alternative = "two.sided")
  return(as.integer(F_res$p.value <= alpha))
}
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
n <- c(20,100,2000)
m <- 2000
alpha_F <- 0.055

## -----------------------------------------------------------------------------
set.seed(68)
power <- matrix(0,2,length(n))
rownames(power) <- c("Count Five","F test")
colnames(power) <- paste("n=",n,sep="")
for (i in 1:length(n)){
  power_CF <- power_F <- numeric(m)
  for (j in 1:m){
        x <- rnorm(n[i],mu1,sigma1)
        y <- rnorm(n[i],mu2,sigma2)
        power_CF[j] <- count5test(x,y)
        power_F[j] <- F_test(x,y,alpha_F)
  }
  power[,i] <- c(mean(power_CF),mean(power_F))
}
kbl(power) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
sk_multi <- function(x){ #a function that computes b_{1,d}
  d <- ncol(x)
  n <- nrow(x)
  x <- scale(x,center = TRUE,scale = FALSE)
  sigma_es <- (n-1)*cov(x)/n
  sigma_es_inv <- chol2inv(chol(sigma_es))
  outer_matrix <- x %*% sigma_es_inv %*% t(x)
  return(sum(outer_matrix^3)/n^2)
}

## -----------------------------------------------------------------------------
set.seed(63120)
alpha <- 0.05
n <- c(10,20,30,50,100,500)
sigma <- matrix(c(1,0.5,0.5,1),nrow = 2)
mu <- c(0,1)
d <- ncol(sigma)
df <- d*(d+1)*(d+2)/6
c <- (n+1)*(n+3)*(d+1)/(n*(n+1)*(d+1) - 6)
cv <- qchisq(1-alpha,df)
p.reject <- numeric(length(n)) #to store simulation results
p.reject_modi <- numeric(length(n))
m <- 1e4 #number of replicates 
for (i in 1:length(n)){
  sktests <- numeric(m)
  sktests_modi <- numeric(m)
  for (j in 1:m){
    x <- mvrnorm(n[i],mu,sigma)
    z <- sk_multi(x)
    sktests[j] <- as.integer(z >= 6*cv/n[i])
    sktests_modi[j] <- as.integer(z >= 6*cv/(n[i]*c[i]))
  }
  p.reject[i] <- mean(sktests)
  p.reject_modi[i] <- mean(sktests_modi)
}
print(p.reject)

## ----echo=FALSE---------------------------------------------------------------
result_simu <- rbind(paste(n),p.reject)
rownames(result_simu) <- c("sample size","estimate")
kbl(result_simu) %>% kable_classic_2(full_width = F,position = 'left')

## ----echo = FALSE-------------------------------------------------------------
result_modi_simu <- rbind(c(10,20,30,50,100,500),p.reject_modi)
rownames(result_modi_simu) <- c("sample size","estimate")
kbl(result_modi_simu) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
set.seed(63121)
alpha <- 0.1
n <- 30
m <- 2500
epsilon <- c(seq(0,0.15,0.01),seq(0.15,1,0.05))
N <- length(epsilon)
pwr <- numeric(N)
d <- 2
df <- d*(d+1)*(d+2)/6
cv <- qchisq(1-alpha,df)
c <- (n+1)*(n+3)*(d+1)/(n*(n+1)*(d+1) - 6)
for (j in 1:N){
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m){
    n2 <- sum(sample(c(0,1),replace = TRUE,size = n,prob = c(1-e,e)))
    if (n2 == 0) {
      x <- mvrnorm(n,c(0,0),diag(2))
      }else {
      if (n2 == n){
        x <- mvrnorm(n,c(0,0),100*diag(2))
      }else {
        x <- rbind(mvrnorm(n-n2,c(0,0),diag(2)),mvrnorm(n2,c(0,0),100*diag(2))) 
      }
      }
    sktests[i] <- as.integer(sk_multi(x) >= 6*cv/(c*n))
  }
  pwr[j] <- mean(sktests)
}

## -----------------------------------------------------------------------------
#plot power vs epsilon
plot(epsilon,pwr,type = "b",xlab = bquote(epsilon),ylim = c(0,1))
abline(h=0.1,lty=3)
se <- sqrt(pwr*(1-pwr)/m) #add standard errors
lines(epsilon,pwr+se,lty=3)
lines(epsilon,pwr-se,lty=3)

## -----------------------------------------------------------------------------
library(bootstrap)
library(boot,warning(.call = FALSE))

## -----------------------------------------------------------------------------
theta.hat <- cor(law$LSAT,law$GPA)
n <- nrow(law)
theta.jack <- numeric(n)
for (i in 1:n){
  theta.jack[i] <- cor(law$LSAT[-i],law$GPA[-i])
}
bias <- (n-1)*(mean(theta.jack) - theta.hat)
se <- sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2))
cat("A jackknife estimate of the bias of the correlation statistic is: ",bias,"\n")
cat("A jackknife estimate of the standard error of the correlation statistic is: ",se)

## -----------------------------------------------------------------------------
exp_mean <- function(x,i){
  mean(x[i,1])
}
n <- nrow(aircondit)
obj <- boot(data = aircondit,statistic = exp_mean,R = 2000)
boot_cis <- boot.ci(obj,conf=0.95,type = c("norm","basic","perc","bca"))
ci.norm <- paste("(",paste(as.character(boot_cis$normal[2]),as.character(boot_cis$normal[3]),sep = ","),")",sep ="")
ci.basic <- paste("(",paste(as.character(boot_cis$basic[4]),as.character(boot_cis$basic[5]),sep = ","),")",sep ="")
ci.perc <- paste("(",paste(as.character(boot_cis$percent[4]),as.character(boot_cis$percent[5]),sep = ","),")",sep ="")
ci.bca <- paste("(",paste(as.character(boot_cis$bca[4]),as.character(boot_cis$bca[5]),sep = ","),")",sep="")
ci <- rbind(ci.norm,ci.basic,ci.perc,ci.bca)
colnames(ci) <- c("Confidence Interval")
rownames(ci) <- c('Standard Normal','Basic','Percentil','BCa')
kbl(ci) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
result <- c(boot_cis$normal[3] - boot_cis$normal[2],boot_cis$basic[5] - boot_cis$basic[4],boot_cis$percent[5] - boot_cis$percent[4],boot_cis$bca[5]-boot_cis$bca[4])
result <- data.frame(matrix(result,nrow=1))
colnames(result) <- c("Standard Normal","Basic","Percentile","BCa")
rownames(result) <- c("Length of interval")
kbl(result) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
n <- nrow(scor)
lambda.hat <- eigen(cov(scor))$values
theta.hat <- lambda.hat[1] / sum(lambda.hat)
# Jacknife estimate
theta.jack <- numeric(n)
for (i in 1:n){
  lambda.jackhat <- eigen(cov(scor[-i,]))$values
  theta.jack[i] <- lambda.jackhat[1] / sum(lambda.jackhat)
}
bias.jack <- (n-1) * (mean(theta.jack) - theta.hat) #jackknife estimate of bias
se.jack <- sqrt((n-1)^2 * var(theta.jack) / n) #jackknife estimate of se
result <- cbind(bias.jack,se.jack)
colnames(result) <- c('Bias',"SE")
rownames(result) <- c('Jackknife')
kbl(result) %>% kable_classic_2(full_width = F,position = 'left')

## ----message=FALSE------------------------------------------------------------
library(DAAG,warning(.call = FALSE))
attach(ironslag)

## -----------------------------------------------------------------------------
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- matrix(0,n*(n-1)/2,2)
# for leave-two-out cross validation
for (j in 2:n){
  for (i in 1:(j-1)){
    y <- magnetic[c(-i,-j)]
    x <- chemical[c(-i,-j)]
    #Linear model
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i,j)]
    e1[(j-1)*(j-2)/2 + i,] <- magnetic[c(i,j)] - yhat1
    #Quadratic model
    J2 <- lm(y~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2]*chemical[c(i,j)] + J2$coef[3]*chemical[c(i,j)]^2
    e2[(j-1)*(j-2)/2 + i,] <- magnetic[c(i,j)] - yhat2
    #Exponential model
    J3 <- lm(log(y)~x)
    logyhat3 <- J3$coef[1] + J3$coef[2]*chemical[c(i,j)]
    yhat3 <- exp(logyhat3)
    e3[(j-1)*(j-2)/2 + i,] <- magnetic[c(i,j)] - yhat3
    #Log-Log model
    J4 <- lm(log(y)~log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2]*log(chemical[c(i,j)])
    yhat4 <- exp(logyhat4)
    e4[(j-1)*(j-2)/2 + i,] <- magnetic[c(i,j)] - yhat4
  }
}

## -----------------------------------------------------------------------------
result <- data.frame(matrix(c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2)),nrow = 1))
colnames(result) <- c("Linear","Quadratic","Exponential","Log-Log")
rownames(result) <- c("Estimate for prediction error")
kbl(result) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
L2 <- lm(magnetic ~ chemical + I(chemical^2))
L2
detach(ironslag)

## ---- message=FALSE-----------------------------------------------------------
library(doParallel)
n_cores <- 2
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

## -----------------------------------------------------------------------------
counttest <- function(x,y){
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx,outy))))
}

## -----------------------------------------------------------------------------
set.seed(83)
n1 <- 20
n2 <- 30
K <- 1:(n1+n2)
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
R <- 999
x <- rnorm(n1,mu1,sigma1)
y <- rnorm(n2,mu2,sigma2)
z <- c(x,y)
D <- numeric(R)
D0 <- counttest(x,y)
for (i in 1:R){
  k <- sample(K,size = n1, replace = FALSE)
  x1 <- z[k]
  x2 <- z[-k]
  D[i] <- counttest(x1,x2)
}

## -----------------------------------------------------------------------------
p <- mean(c(D0,D) >= D0)
print(p)

## -----------------------------------------------------------------------------
hist(D,freq = FALSE,xlab = "replicates of maxout statistic",main = "")
points(D0,0,cex=1,pch=16)
abline(v =D0,col = 'red')

## ----warning=FALSE------------------------------------------------------------
library(RANN)
library(energy)
library(Ball)
library(boot)

## -----------------------------------------------------------------------------
Tn <- function(z,ix,sizes,k){
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix,]
  NN <- nn2(data = z, k= k+1)
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5)
  i2 <- sum(block2 > n1 + .5)
  (i1 + i2) / (k*n)
}
eqdist.nn <- function(z,sizes,k,R){
  boot.obj <- boot(data = z,statistic = Tn, R = R,sim = "permutation",sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts >= ts[1])
  list(statistic = ts[1], p.value = p.value)
}
rMix <- function(n,p){ #function that produces bimodel distribution
  ind <- sample(c(1,2),n,prob = c(0.5,0.5),replace = TRUE)
  y <- matrix(0,n,p)
  y[which(ind == 1),] <- matrix(rnorm(p*length(which(ind == 1)),-0.1,2),ncol=p)
  y[which(ind == 2),] <- matrix(rnorm(p*length(which(ind == 2)),0.1,2),ncol=p)
  return(y)
}

## ----message=FALSE------------------------------------------------------------
mess <- clusterEvalQ(cl,{library("boot")
  library("RANN")
  library("Ball")
  library("energy")})
clusterExport(cl,c("Tn","eqdist.nn","rMix")) 

## -----------------------------------------------------------------------------
test1 <- function(i){
k <- 3; p <- 2; sigma <- 2;
n1 <- n2 <- 50; R <- 999; n <- n1 + n2; N <- c(n1,n2)
p.values <- numeric(3);
x <- matrix(rnorm(n1*p),ncol = p)
y <- cbind(rnorm(n2),rnorm(n2,sd = sigma))
z <- rbind(x,y)
p.values[1] <- eqdist.nn(z,N,k,R)$p.value
p.values[2] <- eqdist.etest(z,sizes = N,R=R)$p.value
p.values[3] <- bd.test(x=x,y=y,R = 999,seed = i*520)$p.value
return(p.values)
}

## -----------------------------------------------------------------------------
test.1 <- parLapply(cl,1:100,test1)
test.1 <- matrix(unlist(test.1),ncol =3,byrow = TRUE)

## -----------------------------------------------------------------------------
alpha <- 0.1
pow1 <- colMeans(test.1 < alpha)
pow1 <- data.frame(matrix(pow1,nrow = 1))
rownames(pow1) <- c("power")
colnames(pow1) <- c("NN","Energy","Ball")
kbl(pow1) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
test2 <- function(i){
k <- 3; p <- 2; mu <- 0.5; sigma <- 2;
n1 <- n2 <- 50; R <- 999; n <- n1 + n2; N <- c(n1,n2)
p.values <- numeric(3);
x <- matrix(rnorm(n1*p),ncol = p)
y <- cbind(rnorm(n2),rnorm(n2, mean = mu, sd = sigma))
z <- rbind(x,y)
p.values[1] <- eqdist.nn(z,N,k,R)$p.value
p.values[2] <- eqdist.etest(z,sizes = N,R=R)$p.value
p.values[3] <- bd.test(x=x,y=y,R = 999,seed = i*520)$p.value
return(p.values)
}

## -----------------------------------------------------------------------------
test.2 <- parLapply(cl,1:100,test2)
test.2 <- matrix(unlist(test.2),ncol =3,byrow = TRUE)

## -----------------------------------------------------------------------------
alpha <- 0.1
pow2 <- colMeans(test.2 < alpha)
pow2 <- data.frame(matrix(pow2,nrow = 1))
rownames(pow2) <- c("power")
colnames(pow2) <- c("NN","Energy","Ball")
kbl(pow2) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
test3 <- function(i){
k <- 3; p <- 2; 
n1 <- n2 <- 50; R <- 999; n <- n1 + n2; N <- c(n1,n2)
p.values <- numeric(3);
x <- matrix(rt(n1*p,df=1),ncol = p)
y <- rMix(n2,p)
z <- rbind(x,y)
p.values[1] <- eqdist.nn(z,N,k,R)$p.value
p.values[2] <- eqdist.etest(z,sizes = N,R=R)$p.value
p.values[3] <- bd.test(x=x,y=y,R = 999,seed = i*520)$p.value
return(p.values)
}

## -----------------------------------------------------------------------------
test.3 <- parLapply(cl,1:100,test3)
test.3 <- matrix(unlist(test.3),ncol =3,byrow = TRUE)

## -----------------------------------------------------------------------------
alpha <- 0.1
pow3 <- colMeans(test.3 < alpha)
pow3 <- data.frame(matrix(pow3,nrow = 1))
rownames(pow3) <- c("power")
colnames(pow3) <- c("NN","Energy","Ball")
kbl(pow3) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
test4 <- function(i){
k <- 3; p <- 2; mu <- 0.5; sigma <- 2;
n1 <-20; n2 <- 200; R <- 999; n <- n1 + n2; N <- c(n1,n2)
p.values <- numeric(3);
x <- matrix(rnorm(n1*p),ncol = p)
y <- cbind(rnorm(n2),rnorm(n2, mean = mu,sd = sigma))
z <- rbind(x,y)
p.values[1] <- eqdist.nn(z,N,k,R)$p.value
p.values[2] <- eqdist.etest(z,sizes = N,R=R)$p.value
p.values[3] <- bd.test(x=x,y=y,R = 999,seed = i*520)$p.value
return(p.values)
}

## -----------------------------------------------------------------------------
test.4 <- parLapply(cl,1:100,test4)
test.4 <- matrix(unlist(test.4),ncol =3,byrow = TRUE)

## -----------------------------------------------------------------------------
alpha <- 0.1
pow4 <- colMeans(test.4 < alpha)
pow4 <- data.frame(matrix(pow4,nrow = 1))
rownames(pow4) <- c("power")
colnames(pow4) <- c("NN","Energy","Ball")
kbl(pow4) %>% kable_classic_2(full_width = F,position = 'left')

## -----------------------------------------------------------------------------
stopCluster(cl)

## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma,x0,N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in  2:N){
    y <- rnorm(1,x[i-1],sigma)
    if (u[i] <= exp(abs(x[i-1]) - abs(y))){
      x[i] <- y
    }else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x=x,k=k))
}
d_lap <- function(x){
  return(exp(-abs(x))/2)
}

## -----------------------------------------------------------------------------
set.seed(94)
sigma <- c(0.3,1,10,30)
n <- length(sigma)
b <- 6001 # discard the burn-in sample
N <- 10000
x0 = 30
rw <- vector(mode="list",length = n)
for (i in 1:4){
  rw[[i]] <- rw.Metropolis(sigma[i],x0,N)
}

## ----echo=FALSE---------------------------------------------------------------
rw2 <- vector(mode="list",length = n)
for (i in 1:n){
  rw2[[i]] <- (rw[[i]]$x)[b:N]
  #rw2[[i]] <- rw2[[i]][which(abs(rw2[[i]]) < 8)]
}
#plot the histogram
p1 <- ggplot(data=data.frame(rw2[[1]]),aes(rw2[[1]]))+
  geom_histogram(aes(y = ..density.. ),binwidth = 0.05)+
  #xlim(-8,8)+
  geom_line(stat="function",fun = d_lap,color = "deeppink")+
  ggtitle(paste("sigma =",sigma[1],sep = " "))+
  xlab("x")
p2 <- ggplot(data=data.frame(rw2[[2]]),aes(rw2[[2]]))+
  geom_histogram(aes(y = ..density.. ),binwidth = 0.05)+
  #xlim(-8,8)+
  geom_line(stat="function",fun = d_lap,color = "deeppink")+
  ggtitle(paste("sigma =",sigma[2],sep = " "))+
  xlab("x")
p3 <- ggplot(data=data.frame(rw2[[3]]),aes(rw2[[3]]))+
  geom_histogram(aes(y = ..density.. ),binwidth = 0.05)+
  #xlim(-8,8)+
  geom_line(stat="function",fun = d_lap,color = "deeppink")+
  ggtitle(paste("sigma =",sigma[3],sep = " "))+
  xlab("x")
p4 <- ggplot(data=data.frame(rw2[[4]]),aes(rw2[[4]]))+
  geom_histogram(aes(y = ..density.. ),binwidth = 0.05)+
  #xlim(-8,8)+
  geom_line(stat="function",fun = d_lap,color = "deeppink")+
  ggtitle(paste("sigma =",sigma[4],sep = " "))+
  xlab("x")
grid.arrange(p1,p2,p3,p4,ncol=2)

## -----------------------------------------------------------------------------
rw3 <- cbind(rw[[1]]$x,rw[[2]]$x,rw[[3]]$x,rw[[4]]$x)
refline <- c(-qexp(0.95),qexp(0.95))
for (j in 1:n){
  plot(rw3[,j],type="l",xlab = bquote(sigma == .(round(sigma[j],3))),ylab="X",ylim=range(rw3[,j]))
  abline(h = refline)
}

## -----------------------------------------------------------------------------
#print out the reject rate
rej_rate <- numeric(n)
acep_rate <- numeric(n)
for (i in 1:n){
  rej <- rw[[i]]$k
  rej_rate[i] <- rej / N
  acep_rate[i] <- 1 - rej_rate[i]
}
acep_rate <- data.frame(matrix(acep_rate,nrow = 1))
rownames(acep_rate) <- c("aceptance rate")
colnames(acep_rate) <- paste("sigma = ",sigma,sep = "")
kbl(acep_rate) %>% kable_classic_2(full_width=FALSE, position = "left")

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi){
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
sigma <- 1
k <- 4
n <- 15000
b <- 1000

#choose overdispersed inital values
x0 <- c(-15,-5,5,15)

#generate the chains
X <- matrix(0,nrow= k,ncol = n)
for (i in 1:k)
  X[i,] <- rw.Metropolis(sigma,x0[i],n)$x

#compute diagnostic statistics
psi <- t(apply(X,1,cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] /(1:ncol(psi))

## -----------------------------------------------------------------------------
#plot psi for the four chains
for (i in 1:k)
  plot(psi[i,(b+1):n],type="l",xlab=i,ylab=bquote(psi))

## -----------------------------------------------------------------------------
#plot the sequence of R-hat statistics
rhat <- rep(0,n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n],type="l",xlab ="",ylab="R")
abline(h = 1.2, lty = 2)

## -----------------------------------------------------------------------------
f <- function(x,df){ #we convert the problem into finding the root 
  q2 <- sqrt(x^2 * df / (df+1-x^2))
  s2 <- 1 - pt(q2,df)
  df <- df - 1
  if (x^2 < (df + 1)){
  q1 <- sqrt(x^2 * df / (df+1-x^2))
  s1 <- 1 - pt(q1,df)
  }else{
    return(NULL)
  }
  return(s1-s2)
}
f_vec <- Vectorize(f,vectorize.args = "x")
f_vec2 <- Vectorize(f,vectorize.args = "df")
k = c(4:25,100,500,1000)
lower <- 0
uppers <- sqrt(k)

## -----------------------------------------------------------------------------

plot(c(0,2),c(-0.01,0.01),type = "n",xlab = "X",ylab = "f",main = "k = 4")
curve(f_vec(x,4),from = 0,to = 1.99,add = TRUE)
abline(h= 0 ,lty = 2)

plot(c(0,26),c(-0.0003,0.0003),type = "n",xlab = "X",ylab = "f",main = "k = 25")
curve(f_vec(x,25),from = 0,to = 4.99,add = TRUE)
abline(h= 0 ,lty = 2)

plot(c(0,32),c(-0.0000003,0.0000003),type = "n",xlab = "X",ylab = "f",main = "k = 1000")
curve(f_vec(x,1000),from = 0,to = 31.5,add = TRUE)
abline(h= 0 ,lty = 2)

plot(c(0,2),c(-0.01,0.01),type = "n",xlab = "X",ylab = "f")
for (i in 1:length(k)){
  curve(f_vec(x,k[i]),from = 0,to = 1.99,add = TRUE)
}

## -----------------------------------------------------------------------------
test <- cbind(f_vec2(1,df = k),f_vec2(1.99,df = k))
c(min(test[,1]),max(test[,2]))

## -----------------------------------------------------------------------------
outs <- vector(mode="list",length = length(k))
for (i in 1:length(k)){
  outs[[i]] <- uniroot(f,lower = 1,upper = 1.99,df = k[i])
}

## -----------------------------------------------------------------------------
roots <- numeric(length(k))
for (i in 1:length(k)){
  roots[i] <- outs[[i]]$root
}
roots <- data.frame(matrix(roots,nrow =1))
colnames(roots) <- paste("k = ",k,sep ="")
rownames(roots) <- c("root")
kbl(roots) %>%  kable_paper() %>%
  scroll_box(width = "900px", height = "100%")

## -----------------------------------------------------------------------------
ABO_type <- function(na,nb,noo,nab,tol,N){
  #initialize the parameter
  ps <- qs <- numeric(N)
  xs <- ys <- numeric(N)
  log_record <- numeric(N)
  ps[1] <- qs[1] <- 0.3
  
  for (i in 2:N){
    #E step
    xs[i-1] <- ps[i-1] * na / (2-ps[i-1]-2*qs[i-1]) 
    ys[i-1] <- qs[i-1] * nb / (2-qs[i-1]-2*ps[i-1])
    
    #M step
    a11 <- 2*noo + 2*na + nab + nb - ys[i-1]
    a1 <- a12 <- xs[i-1] + na + nab
    a2 <- a21 <- ys[i-1] + nb + nab
    a22 <- 2*noo + 2*nb + nab + na - xs[i-1]
    a <- matrix(c(a11,a21,a12,a22),ncol = 2)
    b <- matrix(c(a1,a2),ncol = 1)
    result <- solve(a,b)
    ps[i] <- result[1,1]
    qs[i] <- result[2,1]
    
    if (max(abs(ps[i]-ps[i-1]), abs(qs[i] - qs[i-1])) < tol)
      break
  }
  
  #compute the log-maximum likelihood values in each E-M step
  ps <- ps[1:i]
  qs <- qs[1:i]
  p_est <- ps[i]
  q_est <- qs[i]
  rs <- 1 - ps - qs
  log_record <- na*log(ps^2 + 2*ps*rs) + nb * log(qs^2 + 2*qs*2*rs) + 2*noo*log(rs) + nab*log(2*ps*qs)
  return(list(p_est = p_est, q_est = q_est, ps = ps, qs = qs , log_record = log_record))
}

## -----------------------------------------------------------------------------
na <- 444; nb <- 132; noo <- 361; nab <- 63
tol <- 0.0000001
N <- 1000
pq_em <- ABO_type(na,nb,noo,nab,tol,N)

## -----------------------------------------------------------------------------
cat("MLE of p is: ", pq_em$p_est,"\n")
cat("MLE of q is: ", pq_em$q_est,"\n")

## -----------------------------------------------------------------------------
cat("The values of p in each E-M step are:", pq_em$ps,"\n")
cat("The values of q in each E-M step are:", pq_em$qs,"\n")
cat("The corresponding log-maximum likelihood values are: ",pq_em$log_record)

## -----------------------------------------------------------------------------
#use for loops to fit linear models
formulas <- list(mpg~disp,mpg~I(1/disp),mpg~disp+wt,mpg~I(1/disp)+wt)
out <- vector("list",length(formulas))
for (i in seq_along(formulas)){
  out[[i]] <- lm(formula = formulas[[i]],data = mtcars)
}
print(out)

## -----------------------------------------------------------------------------
#using lapply to fit linear models
lapply(formulas, function(formula) lm(formula = formula,data = mtcars))

## -----------------------------------------------------------------------------
trials <- replicate( 100,t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)

## -----------------------------------------------------------------------------
#Extract the p-value from every trial using sapply and an anonymous function
p_values1 <- sapply(trials,function(x) x$p.value)
print(p_values1)

## -----------------------------------------------------------------------------
#get rid of the anonymous function by using [[ directly
p_values2 <- sapply(trials, "[[",3)
print(p_values2)

## -----------------------------------------------------------------------------
example_list <- list(mtcars, faithful)
lapply(example_list, function(x) vapply(x, mean, numeric(1)))

## -----------------------------------------------------------------------------
lapply_mv <- function(X, FUN, FUN.VALUE,simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE){
    return(out)
    }
  unlist(out)
}
lapply_mv(example_list, mean, numeric(1))

## -----------------------------------------------------------------------------
library(microbenchmark)
library(rmutil,.Call(warning = FALSE))

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
#  List rw_MetropolisC(double sigma,double x_0, int N) {
#      NumericVector x(N);
#      x[0] = x_0;
#      NumericVector u = runif(N);
#      double y;
#      int k = 0;
#      for(int i = 1; i < N; ++i){
#          y = R::rnorm(x[i-1],sigma);
#          if (u[i-1] <= exp(abs(x[i-1]) - abs(y))){
#              x[i] = y;
#          } else{
#              x[i] = x[i-1];
#              ++k;
#          }
#      }
#      return List::create(
#                          _["x"] = x,
#                          _["k"] = k
#                          );
#  }

## -----------------------------------------------------------------------------
set.seed(1)
N <- 20000
sigma <- 1
x_0 <- 0
b <- 2001 #to discard burn-in samples
rwC_pre <- rw_MetropolisC(sigma,x_0,N)
rwC <- rwC_pre[[1]][b:N]

## -----------------------------------------------------------------------------
p1 <- ggplot(data=data.frame(rwC),aes(rwC))+
  geom_histogram(aes(y = ..density.. ),binwidth = 0.05)+
  #xlim(-8,8)+
  geom_line(stat="function",fun = dlaplace,color = "deeppink")+
  ggtitle(paste("sigma =",sigma,"with a Rcpp function",sep = " "))+
  xlab("x")
p1

## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma,x0,N){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in  2:N){
    y <- rnorm(1,x[i-1],sigma)
    if (u[i] <= exp(abs(x[i-1]) - abs(y))){
      x[i] <- y
    }else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x=x,k=k))
}

## -----------------------------------------------------------------------------
set.seed(2)
rwR_pre <- rw.Metropolis(sigma,x_0,N)
rwR <- rwR_pre[[1]][b:N]

## -----------------------------------------------------------------------------
p2 <- ggplot(data=data.frame(rwR),aes(rwR))+
  geom_histogram(aes(y = ..density.. ),binwidth = 0.05)+
  #xlim(-8,8)+
  geom_line(stat="function",fun = dlaplace,color = "deeppink")+
  ggtitle(paste("sigma =",sigma,"with a R function",sep = " "))+
  xlab("x")
p2

## -----------------------------------------------------------------------------
par(mfrow =c(1,2))
qqplot(qlaplace(ppoints((N-b+1))),rwC,main = "Q-Q plot for Cpp")
lines(c(-7,7),c(-7,7),col="red")
qqplot(qlaplace(ppoints((N-b+1))),rwR,main = "Q-Q plot for R")
lines(c(-7,7),c(-7,7),col="red")

## -----------------------------------------------------------------------------
N <- 10000
x_0 <- 0
sigma <- 1
ts <- microbenchmark(rw.MetropolisR = rw.Metropolis(sigma,x_0,N), rw.MetropolisCpp = rw_MetropolisC(sigma,x_0,N))
summary(ts)[,c(1,3,5,6)]

