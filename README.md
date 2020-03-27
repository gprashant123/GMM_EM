# GMM_EM
Estimation of parameters of Gaussian Mixture Models using Expectation maximization algorithm

## PCA.r
To perform dimensionality reduction using PCA
#### Arguments

X - Data matrix

#### Returns 
recons_data - Reconstructed Data Matrix with reduced number of dimensions

## gmm_em.r
To implement Gaussian Mixture Models using Expectation-Maximization Algorithm
#### Arguments

recons_data - Data matrix after dimensionality reduction

clus - Number of Gaussian clusters

#### Returns 
Parameters - A list containing the parameters (mean, covariance matrix and weights) corresponding to each cluster.
