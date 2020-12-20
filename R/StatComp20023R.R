#' @title The original PM 2.5 data by airboxes in Taiwan in consecutive 57 hours.
#' @name PMdata_original
#' @description This data is used to provide an example which the function geo-anomaly can handle.
#' @examples
#' \dontrun{
#' data(PMdata_original)
#' }
NULL

#' @title Anomaly found by applying geo_anomaly() to PMdata_original
#' @name l
#' @description Anomaly found by applying geo_anomaly() to PMdata_original with K = 100, b_threshold = 3, Rmse_threshold = 4
#' @examples
#' \dontrun{
#' data(l)
#' }
NULL


#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of Cpp and R.
#' @examples
#' \dontrun{
#' N <- 10000
#' x_0 <- 0
#' sigma <- 1
#' ts <- microbenchmark(rw.MetropolisR = rw.Metropolis(sigma,x_0,N), 
#' rw.MetropolisCpp = rw_MetropolisC(sigma,x_0,N))
#' summary(ts)[,c(1,3,5,6)]
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rgamma
#' @useDynLib StatComp20023
NULL


#' @title Sequentially and iteratively reweighted squares algorithm for sparse recovery
#' @description SIRS algorithm was proposed by Jinchi Lv and Yingying Fan in 
#' @param A the design matrix
#' @param y the vector of response
#' @param a the parameter in sica penalty
#' @return a estimated vector of coefficients
#' @references Lv, Jinchi; Fan, Yingying. A unified approach to model selection and sparse recovery using regularized least squares. Ann. Statist. 37 (2009), no. 6A, 3498--3528. doi:10.1214/09-AOS683. 
#' @examples
#' \dontrun{
#' s <- 7; n<- 35; p <- 100; r <- 0.1; beta <- c(1,-0.5,0.7,-1.2,-0.9,3,0.55);
#' data <- simu_sirs_generation(s,n,p,r,beta)
#' X <- data[[1]]
#' y <- data[[2]]
#' a <- 0.5
#' beta_hat <- sirs(X,y,a)
#' }
#' @useDynLib StatComp20023
#' @export
sirs <- function(A, y, a = 0.1){
  n <- nrow(A); p <- ncol(A);
  delta <- 1e-6; x0 <- matrix(1,p,1); maxsize <- min(ceiling(n/2),p); eps <- 1/p; thresh <- 1e-6; maxiter <- 50; maxseq <- min(maxiter,maxsize); tol <- 1e-6;
  
  x_ini <- x0;
  rep <- 1;
  move <- 1;
  xp0 <- matrix(0,p,1);
  
  while(move>tol & rep <= maxseq){
    xp <- sirscore(A,y,a,delta,x_ini,maxsize,thresh,maxiter,tol)
    
    num <- sum(xp!= 0)
    if (num < maxsize){
      break
    }else{
      estmod <- (1:p)[xp!=0]
      xd <- xp[estmod,]
      index <- order(abs(xd),decreasing = TRUE)
      
      x_ini <- matrix(eps,p,1)
      x_ini[estmod[index[1:rep]],] <- matrix(1,rep,1)
      
      move <- sqrt(sum((xp - xp0)^2))
      xp0 <- xp
      rep <- rep + 1
    }
  }
  
  if (sum(xp!= 0) > n/2){
    cat("Warning: The solution found is non-sparse. Try a smaller value of a!")
  }
  return(xp)
}


#' @title Function generating data for simulation 1 in SICA paper
#' @description Function that generates data for simulation 1 in the paper of Lv and Fan
#' @seealso \code{\link{sirs}}
#' @param s the number of support
#' @param n nrow in design matrix
#' @param p ncol in design matrix
#' @param r the correlation coefficients
#' @param beta the true coefficients of the model
#' @return the design matrix as well as its response vector as a list
#' @references Lv, Jinchi; Fan, Yingying. A unified approach to model selection and sparse recovery using regularized least squares. Ann. Statist. 37 (2009), no. 6A, 3498--3528. doi:10.1214/09-AOS683. 
#' @examples
#' \dontrun{
#' s <- 7; n<- 35; p <- 100; r <- 0.1; beta <- c(1,-0.5,0.7,-1.2,-0.9,3,0.55);
#' data <- simu_sirs_generation(s,n,p,r,beta)
#' }
#' @import MASS
#' @export
simu_sirs_generation <- function(s,n,p,r,beta){
  support <- 1:s
  mu <- rep(0,p)
  Sig <- diag(p)*(1-r) + matrix(r,p,p)
  A <- mvrnorm(n,mu,Sig)
  A <- sweep(A,2,apply(A,2,mean));
  A <- sweep(A,2,sqrt(apply(A^2,2,sum)),"/")
  y <- A[,support] %*% beta
  return(list(A,y))
}

sirscore <- function(A,y,a,delta,x0,maxsize,thresh,maxiter,tol){
  n <- nrow(A); p <- ncol(A);
  x <- x0
  
  D <- computing_D(x,a)
  
  k <- 1
  update <- 1
  
  while(update > tol & k <= maxiter){
    k <- k+1
    xold <- x
    
    if(p <= ceiling(log(n))*n){
      D1 <- sqrt(D)
      x <- D1 %*% solve(delta*diag(p) + D1 %*% t(A) %*% A %*% D1) %*% D1 %*% t(A) %*% y
    }else{
      x <- D%*%t(A)%*%solve(delta*diag(n) + A%*%D%*%t(A))%*%y
    }
    update <- sqrt(sum((x - xold)^2)) 
    D <- computing_D(x,a)
  }
  
  xthre <- (abs(x) > thresh)
  x <- x*xthre
  num <- sum(xthre)
  
  if (num < maxsize){
    xp <- x
    estmod <- (1:p)[xp!=0]
    A_mod <- A[,estmod]
    xp[estmod] <- solve(t(A_mod) %*% A_mod)%*%t(A_mod)%*%y
    xpthre <- (abs(xp) > thresh)
    xp <- xp*xpthre
  }else{
    xp <- x
  }
}

computing_D <- function(t,a){
  b <- (abs(t) + a)/(a+1)
  return(diag(as.numeric(abs(t)*b)))
}


#' @title An anomaly detection function
#' @description An anomaly detection function with method proposed by Huang G. et al.
#' @param data the data matrix, which should provides device-id for its first column, latitude and longitude with colnames "lat" and "lon" for the second and third column and data in order of time in following columns
#' @param K the number of multi-resolution spline basis functions used
#' @param b_threshold the threshold used in controling the number of devices found in terms of relatively high and low measurements
#' @param Rmse_threshold the threshold used in controling the number of devices found in terms of malfunctioned devices
#' @return A list consisting of devices with relatively high measurements, devices with relatively low measurements and malfunctioned devices, in order.
#' @references Huang G, Chen L-J, Hwang W-H, Tzeng S, Huang H-C. Real-time PM2.5 mapping and anomaly detection from AirBoxes in Taiwan. Environmetrics. 2018;e2537. https://doi.org/10.1002/env.2537
#' @examples
#' \dontrun{
#' data(PMdata_original)
#' K <- 100; b <- 3; Rmse <- 4
#' l <- geo_anomaly(PMdata_original,K,b,Rmse)
#' }
#' @useDynLib StatComp20023
#' @import autoFRK MASS gstat
#' @importFrom stats mahalanobis residuals var
#' @export
geo_anomaly <- function(data,K,b_threshold,Rmse_threshold){
  Device_ID <- data[,1]
  Device_ID <- as.character(Device_ID)
  data1 <- data[,-1]
  dat_list <- dat_creator(data1,K)
  data2 <- dat_list[[1]]
  data3 <- dat_list[[2]]
  dat <- (data1[,-c(1,2)] - data2[,-c(1,2)]) / data3[,-c(1,2)]
  Rmse <- apply(dat^2,1,mean,na.rm = TRUE)
  Rmse <- sqrt(Rmse)
  b <- apply(dat,1,mean,na.rm = TRUE)
  V <- apply(dat,1,var,na.rm = TRUE)
  group1 <- Device_ID[b >= b_threshold]
  group2 <- Device_ID[b <= -b_threshold]
  group3 <- Device_ID[-b_threshold< b & b < b_threshold & Rmse > Rmse_threshold ]
  return(list(group1,group2,group3))
}

dat_creator <- function(data1,K){
  times <- ncol(data1) - 2
  num <- nrow(data1)
  
  knot <- data1[,c(1,2)]
  mm <- matrix(rep(-1,num*times),ncol = times)
  mm <- as.data.frame(mm)
  data2 <- cbind(knot,mm)
  data3 <- cbind(knot,mm)
  
  phia <- mrts(knot,K)
  
  for(i in 1:times){
    y <- data1[,(i+2)]
    fit <- rlm(x = phia, y = y,k2 = 1.345)
    beta_hat <- numeric(K)
    for (j in 1:K){
      beta_hat[j] <- fit$coefficients[[j]]
    }
    data_res <- residuals(fit)
    v <- variogram(object = data_res~1,data = data1,locations = ~lat+lon,cressie = TRUE,width= 0.1)
    v.fit <- suppressWarnings(fit.variogram(v,vgm(model = "Exp",Err = 1)))
    var_y <- v.fit$psill[2]
    lambda <- v.fit$range[2]
    var_epi <- v.fit$psill[1]
    As <- B_compute(knot,var_y,lambda)
    B <- As
    for (j in 1:num){
      B[j,j] <- B[j,j] + var_epi
    }
    C <- solve(B)
    data2[,(i+2)] <- phia%*%beta_hat + As%*%C%*%data_res
    data3[,(i+2)] <- sqrt(var_epi + var_y - mahalanobis(x = As,center = FALSE,cov = C,inverted = TRUE))
  }
  return(list(data2,data3))
}

B_compute <- function(locations,var_y,lambda){
  num <- nrow(locations)
  B <- outer(1:num,1:num,function(a,b) var_y*exp(-sqrt((locations[a,1]-locations[b,1])^2 + (locations[a,2]-locations[b,2])^2)/lambda))
  return(B)
}
