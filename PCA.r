PCA <- function(X){
  #Performing Singular Value Decomposition
  SVD_decomp = svd(t(X)) 
  #getting the diagonal elements of the Sigma Matrix
  D = SVD_decomp$d 
  #Largest Singular Value
  sigma_max = D[1]  
  #Threshold value of ratio of singular value with largest singular value
  k = sum(D/sigma_max > 0.0005)   
  
  #Reconstructed data matrix
  recons_data = SVD_decomp$u[,1:k]%*%diag(SVD_decomp$d)[1:k,1:k]
  
  return(recons_data)
}
