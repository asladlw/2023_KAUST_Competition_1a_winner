library(GpGp)
library(MASS)

method="percentile"
setting="(1, 0.02108, 1.5, 0.27)"
size=500          # number of points

m=60              # number of the nearest neighbors
B=100             # number of bootstrap samples
repea=100         # to compute coverage probability

sigma2=1          # partial sill
beta=0.02108      # range parameter
nu= 1.5           # smoothness parameter
tau2=0.27         # nugget effect


x1=runif(size);x2=runif(size)
locs=cbind(x1,x2)
sig1=matern_isotropic(c(sigma2,beta,nu,tau2/sigma2),locs) #GpGp的nugget value為tau2/sigma2



h=t(chol(sig1))
q=rnorm(repea*size)
k=matrix(q, ncol = size)
g=h %*% t(k) 
sam =t(g)




sam_mle=matrix(NA, repea, 4)
for(w in 1:repea){
  sam_mle[w,]=fit_model(sam[w,], locs ,X=NULL, "matern_isotropic",m_seq=m)$covparms
}

sigma2_CI=c()
beta_CI=c()
nu_CI=c()
tau2_CI=c()
total_bootstrap_mle=c()
for(i in 1:repea){
  sam_mle_a=sam_mle[i,]
  sig2=matern_isotropic(sam_mle_a,locs)
  L=t(chol(sig2))
  x=rnorm(B*size)
  Z=matrix(x, ncol = size)
  Y=L %*% t(Z) 
  bootstrap_sam =t(Y)
  bootstrap_mle=c()
  for(j in 1:B){
    bootstrap_mle=rbind(bootstrap_mle,fit_model(bootstrap_sam[j,], locs ,X=NULL, "matern_isotropic",m_seq=m)$covparms)
  } #每一列為一個拔靴法樣本的最大概似估計值
  #-------------------------------------------------#GpGp的nugget value為tau2/sigma2
  sam_mle_a[4]=sam_mle_a[4]*sam_mle_a[1]
  bootstrap_mle[,4]=bootstrap_mle[,4]*bootstrap_mle[,1]
  #---------------------------------------------------------------------------------
  delta_star=matrix(NA,B,4)
  for(k in 1:B){
    delta_star[k,]=sam_mle_a-bootstrap_mle[k,]  ##bootstrap_mle[k,]-sam_mle_a是difference method,sam_mle_a-bootstrap_mle[k,]是percentile method
  }
  
  delta_star_quantile=rbind(quantile(delta_star[,1],c(0.975,0.025)),
                            quantile(delta_star[,2],c(0.975,0.025)),
                            quantile(delta_star[,3],c(0.975,0.025)),
                            quantile(delta_star[,4],c(0.975,0.025)))
  
  CI=sam_mle_a-delta_star_quantile
  sigma2_CI=rbind(sigma2_CI,CI[1,])
  beta_CI=rbind(beta_CI,CI[2,])
  nu_CI=rbind(nu_CI,CI[3,])
  tau2_CI=rbind(tau2_CI,CI[4,])
  total_bootstrap_mle=rbind(total_bootstrap_mle,bootstrap_mle)
  
}


a=0;b=0;c=0;d=0
for(i in 1:repea){
  if(sigma2>sigma2_CI[i,][1] & sigma2<sigma2_CI[i,][2]){a=a+1}
  if(beta  >beta_CI[i,][1] & beta  <beta_CI[i,][2]){b=b+1}
  if(nu>nu_CI[i,][1] & nu<nu_CI[i,][2]){c=c+1}
  if(tau2>tau2_CI[i,][1] & tau2<tau2_CI[i,][2]){d=d+1}
}
c(a,b,c,d)/repea


