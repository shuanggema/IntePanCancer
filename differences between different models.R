########################################### differences among models with and without shrinkage ############################

## load data
load("data.RData")
source("shrinkage-based functions.R")
## each term of X is normalized and represents for the gene expression 
## each term of Y is the log-transformed overall survival time
## get the gene names
genelist

## calculate the Kaplan-Meier weights for AFT model
W<-list()
p=ncol(X[[1]])
m=length(lengths(X))

for(k in 1:m)
{
  w<-vector()
  y<-Y[[k]]
  x<-X[[k]]
  del<-delta[[k]]
  n=nrow(x)
  comb<-cbind(y,del)
  oo=order(y)
  ocomb<-comb[oo,]
  y<-ocomb[,1]
  del<-ocomb[,2]
  for(i in 1:nrow(x))
  {
    if(i==1){w[i]<-del[1]/n}
    else if(del[i]*1==1)
    {
      multiple=1
      for(j in 1: (i-1))
      {
        multiple=multiple*((n-j)/(n-j+1))^(del[j])
      }
      w[i]<-multiple/(n-i+1)
    }
    else {w[i]=0}
  }
  
  w[oo]=w
  W[[k]]=w*nrow(X[[k]])/sum(w)
}


## center Y & X for each cancer using their weighted means
for(k in 1:m)
{
  y<-Y[[k]]
  x<-X[[k]]
  w<-W[[k]]
  y_<-t(w)%*%y/sum(w)
  x_<-rep(0,p)
  for(i in 1:nrow(x))
  {
    x_=x_+w[i]*x[i,]/sum(w)
  }
  for(i in 1:nrow(x)){y[i]<-sqrt(w[i])*(y[i]-y_)}        ## y
  for(i in 1:nrow(x)){x[i,]<-sqrt(w[i])*(x[i,]-x_)}      ## x
  X[[k]]<-x
  Y[[k]]<-y
}



########################################################### marginal analysis############################################

##### coefficient matrices estimation

  ## marginal analysis without shrinkage
  library(ncvreg)
  set.seed(123)
  coef_marginal_mcp<-matrix(0,p,m)
  for(k in 1: m )
    { for(specific in 1:p)
      { x<-as.matrix(X[[k]][,specific])
        y<-Y[[k]]
        fit= ncvreg(x, y,penalty="MCP", gamma=3)
        lam <- fit$lambda[which.min(BIC(fit))]
        coef_marginal_mcp[specific,k]<-coef(fit,lambda=lam)[-1]   
      }
    }

  ## marginal analysis with two shrinkage-based models
  set.seed(123)
  S<-vector();for(k in 1:m){S<-cbind(S,sign(cor(X[[k]], Y[[k]])))}
  N<-0;for(k in 1:m){N<-N+nrow(X[[k]])}
  gamma=3
  xi=0.01
  
    ## magnitude-based shrinkage
    lambda1= c(0.06);lambda2= c(0.25) ;coef_marginal_magnitude=matrix(0,p,m);for(specific in 1:p){ coef_marginal_magnitude[specific,]=magnitude_marginal(X,Y,lambda1,lambda2,gamma,specific)};coef_marginal_magnitude[coef_marginal_magnitude<0.05]=0
    ## sign-based shrinkage
    lambda1= c(0.25);lambda2= c(0.012);coef_marginal_sign=matrix(0,p,m);for(specific in 1:p){ coef_marginal_sign[specific,]=sign_marginal(X,Y,lambda1,lambda2,gamma,specific)};coef_marginal_sign[coef_marginal_sign<0.05]=0                 
 
    colnames(coef_marginal_magnitude)=colnames(coef_marginal_sign)=colnames(coef_marginal_mcp)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
    rownames(coef_marginal_magnitude)=rownames(coef_marginal_sign)=rownames(coef_marginal_mcp)=gene_list

    
        
##### differences among coefficient matrices of different models 
    
  ## differences in overlapping genes 
    ## without shrinkage
    marginal_overlap_mcp<-vector();for(i in 1:ncol(X[[1]])){marginal_overlap_mcp[i]<-sum(coef_marginal_mcp[i,]!=0)}              
    cbind(as.character(gene_list[order(-marginal_overlap_mcp)]),marginal_overlap_mcp[order(-marginal_overlap_mcp)])
    ## magnitude-based shrinkage
    marginal_overlap_magnitude<-vector();for(i in 1:ncol(X[[1]])){marginal_overlap_magnitude[i]<-sum(coef_marginal_magnitude[i,]!=0)}              
    cbind(as.character(gene_list[order(-marginal_overlap_magnitude)]),marginal_overlap_magnitude[order(-marginal_overlap_magnitude)])
    ## sign-based shrinkage
    marginal_overlap_sign<-vector();for(i in 1:ncol(X[[1]])){marginal_overlap_sign[i]<-sum(coef_marginal_sign[i,]!=0)}              
    cbind(as.character(gene_list[order(-marginal_overlap_sign)]),marginal_overlap_sign[order(-marginal_overlap_sign)])

    
    
  ## differences in relative overlapping 
    ## without shrinkage
    ROL_marginal_mcp<-matrix(0,m,m); for(i in 1:m){for(j in 1:m){ROL_marginal_mcp[i,j]<-length(intersect(which(coef_marginal_mcp[,i]!=0),which(coef_marginal_mcp[,j]!=0)))/length(union(which(coef_marginal_mcp[,i]!=0),which(coef_marginal_mcp[,j]!=0)))}}
    rownames(ROL_marginal_mcp)=colnames(ROL_marginal_mcp)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
    ## magnitude-based shrinkage
    ROL_marginal_magnitude<-matrix(0,m,m); for(i in 1:m){for(j in 1:m){ROL_marginal_magnitude[i,j]<-length(intersect(which(coef_marginal_magnitude[,i]!=0),which(coef_marginal_magnitude[,j]!=0)))/length(union(which(coef_marginal_magnitude[,i]!=0),which(coef_marginal_magnitude[,j]!=0)))}}
    rownames(ROL_marginal_magnitude)=colnames(ROL_marginal_magnitude)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
    ## sign-based shrinkage
    ROL_marginal_sign<-matrix(0,m,m); for(i in 1:m){for(j in 1:m){ROL_marginal_sign[i,j]<-length(intersect(which(coef_marginal_sign[,i]!=0),which(coef_marginal_sign[,j]!=0)))/length(union(which(coef_marginal_sign[,i]!=0),which(coef_marginal_sign[,j]!=0)))}}
    rownames(ROL_marginal_sign)=colnames(ROL_marginal_sign)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
  
    
  ## differences in relative Euclidean distances between estimated coefficient matrices
    ## without shrinkage
    dis_marginal_mcp<-matrix(0,m,m);for(i in 1:m){for(j in 1:m){dis_marginal_mcp[i,j]<-sum((coef_marginal_mcp[,i]-coef_marginal_mcp[,j])^2)/(sqrt(sum((coef_marginal_mcp[,i])^2)*sum((coef_marginal_mcp[,j])^2)))}}
    rownames(dis_marginal_mcp)=colnames(dis_marginal_mcp)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
    ## magnitude-based shrinkage
    dis_marginal_magnitude<-matrix(0,m,m);for(i in 1:m){for(j in 1:m){dis_marginal_magnitude[i,j]<-sum((coef_marginal_magnitude[,i]-coef_marginal_magnitude[,j])^2)/(sqrt(sum((coef_marginal_magnitude[,i])^2)*sum((coef_marginal_magnitude[,j])^2)))}}
    rownames(dis_marginal_magnitude)=colnames(dis_marginal_magnitude)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
    ## sign-based shrinkage
    dis_marginal_sign<-matrix(0,m,m);for(i in 1:m){for(j in 1:m){dis_marginal_sign[i,j]<-sum((coef_marginal_sign[,i]-coef_marginal_sign[,j])^2)/(sqrt(sum((coef_marginal_sign[,i])^2)*sum((coef_marginal_sign[,j])^2)))}}
    rownames(dis_marginal_sign)=colnames(dis_marginal_sign)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
    

    
    
    
    
    

    
    ########################################################### joint analysis############################################
    
    ##### coefficient matrices estimation
    
    ## joint analysis without shrinkage
    library(ncvreg)
    set.seed(123)
    coef_joint_mcp<-matrix(0,p,m)
    for(k in 1:m)
      {
          x<-X[[k]]
          y<-Y[[k]]
          fit= ncvreg(x, y,penalty="MCP", gamma=3)
          lam <- fit$lambda[which.min(BIC(fit))]
          coef_joint_mcp[,k] <- coef(fit, lambda=lam)[-1]
      }
    
    
    ## joint analysis with two shrinkage-based models
    set.seed(123)
    S<-vector();for(k in 1:m){S<-cbind(S,sign(cor(X[[k]], Y[[k]])))}
    N<-0;for(k in 1:m){N<-N+nrow(X[[k]])}
    gamma=3
    xi=0.01

    ## magnitude-based shrinkage
    lambda1_seq = c(0.04);lambda2_seq = c(0,0.065);coef_joint_magnitude=matrix(0,p,m);coef_joint_magnitude=magnitude_joint(X,Y,lambda1,lambda2,gamma)
    ## sign-based shrinkage
    lambda1=0.007;lambda2=0.00001;coef_joint_sign=matrix(0,p,m); coef_joint_sign=sign_joint(X,Y,lambda1,lambda2,gamma)              
    
    colnames(coef_joint_magnitude)=colnames(coef_joint_sign)=colnames(coef_joint_mcp)=c("BRCA","LUAD","BLCA","LUSC","GBM","HNSC","LAML","OV","PAAD")
    rownames(coef_joint_magnitude)=rownames(coef_joint_sign)=rownames(coef_joint_mcp)=gene_list
    
