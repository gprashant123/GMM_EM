gmm_em <- function(recons_data,clus){

library(pracma)

#Preprocess
#Mean-center normalization
norm_recons_dat = scale(recons_data)

#EM Algorithm
size = dim(norm_recons_dat)


Parameters = list() #initiating parameter list

#Multivariate Gaussian Distribution 

multivar_gauss <- function(X,mu,sigma){
  EPS = 2^-1023
  size = dim(X)
  n = size[2]
  pX = exp((-0.5)*(X - pracma::repmat(mu,size[1],1))%*%inv(sigma)%*%t(X - pracma::repmat(mu,size[1],1))) *(1/sqrt(det(sigma)))
  p = diag(pX)
  p = replace(p,which(p %in% 0),EPS)
  return(p)
}

#Initialize Parameters


#library(mvtnorm)

Parameters$means = list()
Parameters$covariances = list()
Parameters$weights = runif(clus,0,1)
Parameters$weights = Parameters$weights/sum(Parameters$weights)

for (i in 1:clus){
  #Parameters$means[[i]] = t(runif(size[2],-1,1))
  #Parameters$means[[i]] = t(runif(size[2],min(norm_recons_dat),max(norm_recons_dat)))
  #Parameters$means[[i]] = array(0,c(1,size[2]))
  #Parameters$means[[i]] = colSums(norm_recons_dat)/size[1]
  Parameters$means[[i]] = norm_recons_dat[sample(c(1:size[1]),1),]
  #CoV = lower.tri(matrix(runif(size[2]*size[2]), ncol=size[2]))
  #CoV = CoV + t(CoV)
  #diag(CoV) = diag(CoV)/2
  Parameters$covariances[[i]] = diag(size[2])
  
}

p_k = 0
for (q in 1:clus){
  
  p_k = p_k + Parameters$weights[q]*multivar_gauss(norm_recons_dat,Parameters$means[[q]],Parameters$covariances[[q]])    
}
p_X_i = sum(log(p_k))
p_X_init = p_X_i


p_X_final = Inf
tol = 10^-4 #tolerance for convergence
iter = 1
R = array(0, c(size[1],clus))

#Iteration
 
while (is.nan(abs(p_X_final-p_X_init)/p_X_init) || abs(p_X_final-p_X_init)/abs(p_X_init) > tol){
  
  print(iter)
  iter = iter + 1
  p_X_init = p_X_final
  
  for (i in 1:clus){
    print(i)
    R[,i] = Parameters$weights[i]*t(t(multivar_gauss(norm_recons_dat,Parameters$means[[i]],Parameters$covariances[[i]])))
  }

  R = R/rowSums(R)
  mc = colSums(R)
  Parameters$weights = mc/size[1]
  
  p_k = 0
  for (j in 1:clus){
    print(j)
    Parameters$means[[j]] = (1/mc[j])*(t(R[,j]) %*% norm_recons_dat)
    Parameters$covariances[[j]] = (1/mc[j])*(t(norm_recons_dat - pracma::repmat(Parameters$means[[j]],size[1],1)) %*% diag(R[,j])) %*% (norm_recons_dat - pracma::repmat(Parameters$means[[j]],size[1],1))
    p_k = p_k + Parameters$weights[j]*multivar_gauss(norm_recons_dat,Parameters$means[[j]],Parameters$covariances[[j]])    
  }
  
  p_X_final = sum(log(p_k))
  print(p_X_final)
}

return(Parameters)

}

X = load("SRA635094_SRS2725547.sparse-RPKM.RData")
recons_data = PCA(X)
#Number of clusters
clus = 3
Parameters = gmm_em(recons_data,clus)

for (i in 1:clus){
  R[,i] = Parameters$weights[i]*t(t(multivar_gauss(norm_recons_dat,Parameters$means[[i]],Parameters$covariances[[i]])))
}
R = R/rowSums(R)

#Clustering result
classification = cbind(1:nrow(R), max.col(R, 'first'))
View(classification)
