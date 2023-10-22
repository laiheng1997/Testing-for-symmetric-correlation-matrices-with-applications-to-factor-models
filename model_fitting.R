## packages

library(MARSS)

### Input the dimension of the series and common factors,
### output the matrix structure of matrix Z and matrix B for the model list in function MARSS()

############### Input ######################
## p: the order of AR                     ##
## k: the dimension of the series         ##
## r: the dimension of the common factors ##
## symmetry: CCF is symmetry or not       ##
## seed: random seed                      ##
############################################

MARSS_DFM_Z_B_AR <- function(p, k, r, symmetry = TRUE, seed = NULL){
  
  if (k < r)
    stop("The dimension of the common factors (r) is larger than the dimension of the series (k)")
  
  if (p <= 0)
    stop("The order of AR must greater than zero")
  
  tmp = c()
  for (i in 1:r){
    tmp = c(tmp, paste0("L_",1:k,"_",i))
  }
  
  ### Construct the loading matrix L (in "MARSS" called Z)
  
  Z = matrix(tmp, k, r)
  Z_star = matrix(list(0), k, r*(p-1))
  Z = cbind(Z, Z_star)
  
  ### Set initial state (x0)
  
  set.seed(seed)
  x0 = matrix(runif(r * p, -0.5, 0.5), r * p, 1)
  
  ### Construct the AR matrix Phi (in "MARSS" called B)
  
  B <- list()
  
  if (symmetry == TRUE){
    for (j in 1:p){
      B[[j]] <- matrix(list(0), r, r)
      for (i in 1:r){
        B[[j]][i,i] <- paste0("b_",i,"_",i,"_",j)
      }
    }
  } else {
    for (j in 1:p){
      tmp = c()
      for (i in 1:r){
        tmp = c(tmp, paste0("b_",1:r,"_",i,"_",j))
      }
      B[[j]] = matrix(tmp, r, r)
    }
  }
  
  if (p == 1){
    B = B[[1]]
    Q = "identity"
    return(list(Z=Z, B=B, x0=x0, Q=Q))
  } else {
    B_star <- matrix(list(0), r*p, r*p)
    B_star[1:r,] <- do.call(cbind,B)
    B_star[(r+1):(r*p),1:(r*(p-1))] <- diag(rep(list(1),r*(p-1)))
    Q = matrix(rep(0,r^2*p^2), nrow = r*p, ncol = r*p)
    Q[1:r,1:r] = diag(r)
    return(list(Z = Z, B = B_star, x0=x0, Q=Q))
  }
}

####################### Output: ############################
## Z: the matrix structure of matrix L (dimension: k*r)   ##
## B: the matrix structure of matrix Phi (dimension: r*r) ##
## x0: the initial state x(0)                             ##
############################################################


# -------------------------------------------------------------------------

## Model fitting

################################ Input ######################################
## data: data with dimension T*k (T:time length, k = dimension of series)  ##
## p: the order of AR                                                      ##
## r: the dimension of the common factors                                  ##
## symmetry: CCF is symmetry or not                                        ##
## seed: random seed                                                       ##
## ...: other arguments in function "MARSS"                                ##
#############################################################################

MARSS_model_fit <- function(data, p, r, symmetry = FALSE, seed = NULL, ...){
  
  kk <- ncol(data)
  
  param = MARSS_DFM_Z_B_AR(p, kk, r, symmetry, seed)
  
  model.list = list(
    B = param$B,
    U = "zero",
    Q = param$Q,
    Z = param$Z,
    A = "zero",
    R = "diagonal and equal",
    x0 = param$x0,
    tinitx = 0
  )
  
  fit = MARSS::MARSS(t(data), model = model.list, ...)
  
  return(fit)
}

######## Output ###########
## model fitting result  ##
###########################
