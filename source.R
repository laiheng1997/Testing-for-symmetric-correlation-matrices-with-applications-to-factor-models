##packages

##install.packages("gtools")
library("gtools")

#------------------------------------------------------------------------------

## Construct the B matrix for contrasting c_{ij} and c_{ji} in CCF

############## Input ###############
## k: the dimension of the series ##
####################################

symB <- function(k){
  l <- k-1
  B <- matrix(0,nrow=k*(k-1)/2,ncol=k^2)
  c <- 1
  for (i in 1:l){
    s <- i+1
    for (j in s:k){
      B[c,(i-1)*k+j] <- 1
      B[c,(j-1)*k+i] <- -1
      c <- c+1
    }
  }
  return(B)
}

########### Output ############
## the contrasting matrix B  ##
###############################

#-------------------------------------------------------------------------------

## Use data find vectorized C_hat or vectorized R_hat

################################## Input ###################################
## data: the time series data (in matrix form, row: time, col:covariates) ##
## type: "covariance" -- output vectorized auto-covariance matrix c_l     ##
##       "correlation" -- output vectorized auto-correlation matrix r_l   ##
############################################################################

vec_c <- function(data, type=c("covariance","correlation")){
  
  Y_t <- t(data)
  
  ## Find X_bar
  s<-0
  for(i in 1:ncol(Y_t)){
    s<-s+Y_t[,i]
  }
  mean<-s/ncol(Y_t)
  
  ## Find k and T
  k<-nrow(Y_t)
  tt<-ncol(Y_t)
  
  type <- match.arg(type)
  
  ## Find C_hat (vecterized autocovariance matrix)
  
  C_l<-matrix(NA,nrow=k^2,ncol=2*tt+1)
  T_bound<-tt-1
  for(j in 0:T_bound){
    down<-j+1
    C<-matrix(0,nrow=k,ncol=k)
    for(i in down:tt){
      C<-C+(Y_t[,i]-mean)%*%t(Y_t[,i-j]-mean)} ## future*past
    C_l[,tt+j+1]<-as.vector(C/tt)
    C_l[,tt-j+1]<-as.vector(t(C)/tt)
  }
  
  C_l<-C_l[,c(-1,-2*tt-1)]
  
  if( k == 1){
    C_l<-t(as.matrix(C_l,nrow=1,ncol=tt))
  }
  
  if (type == "covariance"){
    return(C_l)
  }
  
  if (type == "correlation"){
    corrr<-matrix(NA,nrow=k^2,ncol=2*tt-1)
    C_0<-matrix(C_l[,tt],ncol=k,nrow=k)
    o<-tt-1
    for (l in -o:o){
      covacf1<-matrix(C_l[,tt+l],ncol=k,nrow=k)
      corracf1<-matrix(NA,ncol=k,nrow=k)
      for (i in 1:k){
        for (j in 1:k){
          corracf1[i,j]<-covacf1[i,j]/(sqrt(C_0[i,i])*sqrt(C_0[j,j]))
        }
        corrr[,l+tt]<-as.vector(corracf1)
      }
    }
    return(corrr)
  }
}

################################## Output ###################################
## the vectorized auto-covariance or auto-correlation function c_l or r_l  ##
## (Note: here the vectorized lag-l sample covariance matrix is            ##
##        the (tt+l)-th column of c_l)                                     ##
#############################################################################

#-------------------------------------------------------------------------------

## Given l, find lag-l autocovariance matrix and correlation matrix

################################## Input ###################################
## data: the time series data (in matrix form, row: time, col:covariates) ##
## type: "covariance" -- output vectorized auto-covariance matrix C_l     ##
##       "correlation" -- output vectorized auto-correlation matrix R_l   ##
## lagl: the lag-l                                                        ##
############################################################################

laglcov <- function(data, type=c("covariance","correlation"), lagl){
  
  k<-ncol(data)
  tt<-nrow(data)
  
  type <- match.arg(type)
  
  ans<-vec_c(data,type=type)
  
  return(matrix(ans[,tt+lagl],ncol=k,nrow=k))
  
}

########################### Output #############################
## the auto-covariance or auto-correlation matrix C_l or R_l  ##
################################################################

#------------------------------------------------------------------------------

## Modified Bartlett window (Melard et. al, 1986)

########## Input ############
## x: input variable       ##
#############################

barwin<-function(x){
  if(abs(x)<=1){
    return(1-abs(x))
  }
  else return(0)
}

###### Output #######
## output variable ##
#####################

#------------------------------------------------------------------------------

## Theta and Delta function for aysmptotic variance

# Theta-function (According to the equation in Roy, 1989)

############################## Input ################################
## c_hat: the vectorized acf c_ell constructed by function "vec_c" ##
## i,j,l,m: the index of c_{ij} and c_{lm}                         ##
## p: time lag p                                                   ##
#####################################################################

THETA <- function(c_hat,i,j,l,m,p){
  k <- sqrt(dim(c_hat)[1])
  if(p>=0){
    a<-c(rep(0,p),(c_hat[(j-1)*k+i,]))%*%c((c_hat[(m-1)*k+l,]),rep(0,p))  
  }
  else{
    a<-c((c_hat[(j-1)*k+i,]),rep(0,-p))%*%c(rep(0,-p),(c_hat[(m-1)*k+l,]))
  }
  return(a)
}

######## Output ########
## THETA_p (i,j,l,m)  ##
########################

# Delta-function (According to the equation in Roy, 1989)

############################## Input ################################
## c_hat: the vectorized acf c_ell constructed by function "vec_c" ##
## i,j,l,m: the index of c_{ij} and c_{lm}                         ##
## p: time lag p                                                   ##
#####################################################################

DELTA <- function(c_hat,i,j,l,m,p){
  k <- sqrt(dim(c_hat)[1])
  tt <- (dim(c_hat)[2]+1)/2
  a <- THETA(c_hat,i,j,l,m,p)/sqrt(c_hat[(i-1)*k+i,tt]*c_hat[(j-1)*k+j,tt]*c_hat[(l-1)*k+l,tt]*c_hat[(m-1)*k+m,tt])
  return(a)
}

######## Output ########
## DELTA_p (i,j,l,m)  ##
########################

#-------------------------------------------------------------------------------

## Given two different lags and their vectorized auto-covariance and auto-coorelation function, find its corresponding W and V

## input the vecterized acvf C_l, vecterized acf and time lag-p,q, p,q>=0
## output the asymptotically variance matrix W_l and V_l

############################## Input ################################
## c_hat: the vectorized acf c_ell constructed by function "vec_c" ##
## r_hat: the vectorized acf r_ell constructed by function "vec_c" ##
## lag.k: the first lag-k                                          ##
## lag.h: the second lag-h                                         ##
## THETA: the Theta function                                       ##
## DELTA: the Delta function                                       ##
## window: the modified window for consistent estimator            ##
## output: covariance -- output for W_{k,h}                        ##
##         correlation -- output for V_{k,h}                       ##
#####################################################################

V_twolag <- function(c_hat, r_hat, lag.k, lag.h, THETA, DELTA, window = NA, output = c("covariance","correlation")){ # Chat must put matrix
  
  k <- sqrt(dim(c_hat)[1])
  tt <- (dim(c_hat)[2]+1)/2
  
  p <- lag.k
  q <- lag.h
  
  output <- match.arg(output)
  
  if (is.function(window)){
    for (i in 1:ncol(c_hat)){
      c_hat[,i]<-c_hat[,i]*window((i-tt)*1/sqrt(tt))
    }  
  }
  
  perm <- permutations(n = k, r = 4, v=seq(1:k), repeats.allowed = TRUE)
  
  if(output == "covariance"){
    W_hat<-matrix(NA,nrow=k^2,ncol=k^2)

    for (entry in 1:nrow(perm)){
      i<-perm[entry,1];j<-perm[entry,2];l<-perm[entry,3];m<-perm[entry,4];
      W_hat[(j-1)*k+i,(m-1)*k+l]<-
        THETA(c_hat,i,l,j,m,q-p)+THETA(c_hat,j,l,i,m,p+q) # (According to the equation (1) in Roy, 1989)
    }
    return(list("What"=W_hat))
  }

  if(output == "correlation"){
    V_hat<-matrix(NA,nrow=k^2,ncol=k^2)
    
    for (entry in 1:nrow(perm)){
      a<-perm[entry,1];b<-perm[entry,2];d<-perm[entry,3];e<-perm[entry,4];
      V_hat[(b-1)*k+a,(e-1)*k+d]<-
        0.5*r_hat[(b-1)*k+a,tt+p]*r_hat[(e-1)*k+d,tt+q]*
        (DELTA(c_hat,a,d,a,d,0)+DELTA(c_hat,a,e,a,e,0)+DELTA(c_hat,b,d,b,d,0)+DELTA(c_hat,b,e,b,e,0))-
        r_hat[(b-1)*k+a,tt+p]*(DELTA(c_hat,a,d,a,e,q)+DELTA(c_hat,b,d,b,e,q))-
        r_hat[(e-1)*k+d,tt+q]*(DELTA(c_hat,b,d,a,d,p)+DELTA(c_hat,b,e,a,e,p))+
        DELTA(c_hat,a,d,b,e,q-p)+DELTA(c_hat,b,d,a,e,p+q) # (According to the equation (5) in Roy, 1989)
    }
    return(list("Vhat"=V_hat))
  }
}

####### Output #######
## What: W_{k,h}    ##
## Vhat: V_{k,h}    ##
######################

#-------------------------------------------------------------------------------


## the ij-entry of V_hat by given covariance, correlation, i,j,l,m-entry and lag-k,lag-h

############################## Input ################################
## c_hat: the vectorized acf c_ell constructed by function "vec_c" ##
## r_hat: the vectorized acf r_ell constructed by function "vec_c" ##
## lag.k: the first lag-k                                          ##
## lag.h: the second lag-h                                         ##
## THETA: the Theta function                                       ##
## DELTA: the Delta function                                       ##
## i,j,l,m: the index of c_{ij}(lag.k) and c_{lm}(lag.h)           ##
#####################################################################

V_entry = function(c_hat, r_hat, lag.k, lag.h, THETA, DELTA, i, j, l, m){
  
  k<-sqrt(dim(c_hat)[1])
  tt<-(dim(c_hat)[2]+1)/2
  
  p <- lag.k
  q <- lag.h
  
  a <- i; b <- j; d <- l; e <- m;
  
  return(0.5*r_hat[(b-1)*k+a,tt+p]*r_hat[(e-1)*k+d,tt+q]*
           (DELTA(c_hat,a,d,a,d,0)+DELTA(c_hat,a,e,a,e,0)+DELTA(c_hat,b,d,b,d,0)+DELTA(c_hat,b,e,b,e,0))-
           r_hat[(b-1)*k+a,tt+p]*(DELTA(c_hat,a,d,a,e,q)+DELTA(c_hat,b,d,b,e,q))-
           r_hat[(e-1)*k+d,tt+q]*(DELTA(c_hat,b,d,a,d,p)+DELTA(c_hat,b,e,a,e,p))+
           DELTA(c_hat,a,d,b,e,p-q)+DELTA(c_hat,b,d,a,e,p+q))
}

############### Output ##############
## the V_hat (i,j,l,m)-th entry    ##
#####################################
