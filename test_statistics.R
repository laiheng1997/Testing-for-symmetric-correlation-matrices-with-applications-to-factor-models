library("MTS")
library("pracma")
library("mvtnorm")
library("gtools")
library("magic")

# -------------------------------------------------------------------------

## The data generation

#################### Input #########################
## tt: time length of the series                  ##
## ar: the order of ar                            ##
## p1: the matrix phi of VAR                      ##
## sig_eta: the covariance matrix Sigma_eta       ##
## sig_xi: the covariance matrix Sigma_xi         ##
## L: the loading matrix L                        ##
####################################################

data_generation <- function(tt, ar, p1, sig_eta, L, sig_xi){
  m1 = VARMAsim(tt, arlags=ar, phi=p1, sigma=sig_eta)
  zt = m1$series
  k = dim(L)[1]
  if (dim(L)[2] != dim(zt)[2])
    stop("The dimension of f_t is different with the number of columns of L")
  simyt = matrix(NA,tt,k)
  for (j in 1:tt){
    simyt[j,] = t(L%*%(zt[j,])) + rmvnorm(1, mean=rep(0,k), sigma = sig_xi)
  }
  return (simyt)
}

##### Output ######
## the series    ##
###################

# -------------------------------------------------------------------------

## The Wald test statistics based on sample CCF

################################ Input #####################################
## data: data with dimension T*k (T:time length, k = dimension of series) ##
## lagmax: the order of lags                                              ##
## window: the modified window                                            ##
## port: evaluate the Portmanteau test or not                             ##
############################################################################

symmetry_CCF_test <- function(data, lagmax, window, port = T){
  
  kk = ncol(data)
  nT = nrow(data)
  B = symB(kk)
  
  cc = vec_c(data, type="covariance")
  rr = vec_c(data, type="correlation")
  
  if (port == T){
    
    if (lagmax == 1){
      stop("Not to do Portmanteau test for lag = 1!")
    }
    
    Voutput = array(NA, dim=c(kk^2, kk^2, lagmax, lagmax))
    
    for (i in 1:lagmax){
      for (j in i:lagmax){
        Voutput[,,i,j] = V_twolag(cc, rr, i, j, THETA, DELTA, window = window, output="correlation")$Vhat
      }
    }
    
    test = c()
    cumtest = c()
    
    for (i in 1:lagmax){
      temp = nT*t(B%*%rr[,nT+i])%*%solve(B%*%Voutput[,,i,i]%*%t(B))%*%B%*%rr[,nT+i]
      test = c(test,temp)
    }
    
    cumtest = c(cumtest,test[1])
    V = Voutput[,,1,1]
    BigB = B
    rall = rr[,nT+1]
    
    for (i in 2:lagmax){
      V = adiag(V,Voutput[,,i,i])
      BigB = adiag(BigB,B)
      rall = c(rall,rr[,nT+i])
      for (j in 1:(i-1)){
        colindex = ((i-1)*kk^2+1):(i*kk^2)
        rowindex = ((j-1)*kk^2+1):(j*kk^2)
        V[rowindex,colindex] = Voutput[,,j,i]
        V[colindex,rowindex] = t(Voutput[,,j,i])
      }
      temp = nT*t(BigB%*%rall)%*%solve(BigB%*%V%*%t(BigB))%*%BigB%*%rall
      cumtest = c(cumtest,temp)
    }
    
    cumtest = as.data.frame(cbind(test,cumtest))
    colnames(cumtest) = c("Individual Lag Test", "Portmanteau Test")
    rownames(cumtest) = paste0("Lag", 1:lagmax)
    cumtest$p_val_ind = 1 - pchisq(cumtest[,1], df = kk*(kk-1)/2)
    cumtest$p_val_port = 1 - pchisq(cumtest[,2], df = (1:lagmax)*kk*(kk-1)/2)
    
  } else {
    
    Voutput = array(NA, dim=c(kk^2, kk^2, lagmax, lagmax))
    
    for (i in 1:lagmax){
      Voutput[,,i,i] = V_twolag(cc, rr, i, i, THETA, DELTA, window = window, output="correlation")$Vhat
    }
    
    test = c()
    
    for (i in 1:lagmax){
      temp = nT*t(B%*%rr[,nT+i])%*%solve(B%*%Voutput[,,i,i]%*%t(B))%*%B%*%rr[,nT+i]
      test = c(test,temp)
    }
    
    cumtest = as.data.frame(test)
    colnames(cumtest) = c("Individual Lag Test")
    rownames(cumtest) = paste0("Lag", 1:lagmax)
    cumtest$p_val_ind = 1 - pchisq(cumtest[,1], df = kk*(kk-1)/2)
  }
  
  return(cumtest)
}

###################### Output ######################
## the testing result (matrix)                    ##
## Individual Lag Test: test statistics Q_ell     ##
## Portmanteau Test   : test statistics Q_{1:ell} ##
## p_val_ind          : the p-values of Q_ell     ##
## p_val_port         : the p-values of Q_{1:ell} ##
####################################################

# -------------------------------------------------------------------------

## The test statistics based on BH procedures

################################ Input #####################################
## data: data with dimension T*k (T:time length, k = dimension of series) ##
## lagmax: the order of lags                                              ##
## window: the modified window                                            ##
############################################################################

symmetry_BH_test <- function(data, lagmax, window){
  
  kk = ncol(data)
  nT = nrow(data)
  B = symB(kk)
  dimm = kk*(kk-1)/2
  
  cc = vec_c(data, type = "covariance")
  rr = vec_c(data, type = "correlation")
  
  Z_j = matrix(NA, nrow = kk*(kk-1)/2, ncol = lagmax)
  
  if (is.function(window)){
    for (m in 1:ncol(cc)){
      cc[,m]<-cc[,m]*window((m-nT)*1/sqrt(nT))
    }  
  }
  
  if (is.null(colnames(data))){
    idx = 1:kk
  } else {
    idx = colnames(data)
  }
  
  tmp = c()
  
  for (j in 1:dimm){
    aa = which(diag(dimm)[,j]%*%B!=0)
    a = ifelse(aa[1]%%kk==0, kk, aa[1]%%kk)
    b = ceiling(aa[1]/kk)
    d = ifelse(aa[2]%%kk==0, kk, aa[2]%%kk)
    e = ceiling(aa[2]/kk)
    
    tmp = c(tmp, paste0(idx[a],":",idx[b]))
    
    for (i in 1:lagmax){
      Z_j[j,i] = sqrt(nT)*diag(kk*(kk-1)/2)[,j]%*%B%*%rr[,nT+i]/
        sqrt(
            V_entry(cc,rr,i,i,THETA,DELTA,a,b,a,b)+
            V_entry(cc,rr,i,i,THETA,DELTA,d,e,d,e)-
            V_entry(cc,rr,i,i,THETA,DELTA,a,b,d,e)-
            V_entry(cc,rr,i,i,THETA,DELTA,d,e,a,b)
        )
    }
  }

  pval_ind = matrix(NA, nrow = dimm, ncol = lagmax)
  rownames(pval_ind) = tmp
  pval_adj_ind = list()
  pval_adj_port = list()
  tmp = c()
  
  for (i in 1:lagmax){
    Z_j[,i] = abs(Z_j[,i])
    pval_ind[,i] = 2*(1-pnorm(Z_j[,i],0,1))
    tmp = c(tmp, pval_ind[,i])
    pval_adj_ind[[i]] = p.adjust(pval_ind[,i], method="fdr")
    pval_adj_port[[i]] = p.adjust(tmp, method="fdr")
  }

  return(list(p.value = pval_ind, p.adj_ind = pval_adj_ind, p.adj_port = pval_adj_port))
    
}

################################# Output ###############################
## the testing result (list)                                          ##
## p.value    : the p-value of Z_j^{ell}                              ##
## p.adj_ind  : the adjusted p-values of Z_j^{ell} for individual lag ##
## p.adj_port : the adjusted p-values of Z_j^{ell} for multiple lags  ##
########################################################################
