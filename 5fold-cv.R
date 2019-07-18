############################################ cross validation to select the best parameters################################

## load data
load("data.RData")
source("shrinkage-based functions.R")
## each term of X is normalized and represents for the gene expression 
## each term of Y is the log-transformed overall survival time

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




######################################################### for marginal analysis###########################################
lambda1_seq = c(0.04,0.06,0.8);lambda2_seq = c(0.25)                                       ##maginitude-based shrinkage
#lambda1_seq = c(0.1,0.25,0.03);lambda2_seq = c(0.012)                                     ##sign-based shrinkage

S<-vector();for(k in 1:m){S<-cbind(S,sign(cor(X[[k]], Y[[k]])))}
N<-0;for(k in 1:m){N<-N+nrow(X[[k]])}
gamma=3 
xi=0.01
result<-matrix(0,length(lambda1_seq)*length(lambda2_seq),14)
RSS<-rep(0,length(lambda1_seq)*length(lambda2_seq))
time=1

set<-list();for (k in 1: m){set[[k]]=sample(nrow(X[[k]]),nrow(X[[k]]),replace = F)}

for(qq in 1:length(lambda1_seq))
{for(ww in 1:length(lambda2_seq))
{
  lambda1<-lambda1_seq[qq] 
  lambda2<-lambda2_seq[ww]
  
  
  rss<-0
  for (kk in 1: 5)
    {  set.seed(kk*8)
       ## divide the data set into training set and test set
       X_train=list();X_test=list();Y_train=list();Y_test=list()
       for(k in 1:m)
        {  start=(floor(length(set[[k]])/5)*(kk-1)+1)
           end=floor(length(set[[k]])/5)*kk  
           train<-set[[k]][-c(start:end)]
           X_train[[k]]<-X[[k]][train,];Y_train[[k]]<-as.matrix(Y[[k]][train,])
           X_test[[k]]<-X[[k]][-train,];Y_test[[k]]<-as.matrix(Y[[k]][-train,])
         }
    
       ##count coefficients
       BETA=vector()
       for(specific in 1:p)
         { 
            beta=magnitude_marginal(X_train,Y_train,lambda1,lambda2,gamma,specific)                 ## magnitude-based shrinkage
            #beta<-sign_marginal(X_train,Y_train,lambda1,lambda2,gamma,specific)                    ## sign-based shrinkage
            BETA<-rbind(BETA,beta)
          }
        BETA[BETA<0.05]=0

        ##count rss 
        for(k in 1:m)
          {
            Rss<-rep(0,nrow(X_test[[k]]))
            for(row in 1:nrow(X_test[[k]])){Rss[row]=(Y_test[[k]][row,1]-X_test[[k]][row,]%*%BETA[,k])^2}
            Rss[which(Rss==0)]=min(Rss[which(Rss!=0)])
            rss=rss+sum(Rss)
          }
    }
  
  df=vector()
  result[time,1]=lambda1;result[time,2]=lambda2;result[time,3]=sum(BETA!=0);result[time,4]=log(rss)/N
  sum=0;for(j in 1:p){if(sum(BETA[j,]==0)==0){sum=sum+1}};result[time,5]<-sum
  for(k in 1:m){result[time,5+k]<-sum(BETA[,k]!=0)}
  time=time+1
  
  print(k)
}
  
}


colnames(result)<-c("lambda1","lambda2","total-df","log_rss/N","overlap","data1_nonzero_coef","data2_nonzero_coef","data3_nonzero_coef","data4_nonzero_coef","data5_nonzero_coef","data6_nonzero_coef","data7_nonzero_coef","data8_nonzero_coef","data9_nonzero_coef")
result












######################################################### for joint analysis###########################################
lambda1_seq = c(0.04,0.05);lambda2_seq = c(0.065)                                       ## maginitude-based shrinkage
#lambda1_seq = c(0.005,0.007,0.009);lambda2_seq = c(0.00001)                                  ## sign-based shrinkage

S<-vector();for(k in 1:m){S<-cbind(S,sign(cor(X[[k]], Y[[k]])))}
N<-0;for(k in 1:m){N<-N+nrow(X[[k]])}

gamma=3
xi=0.01
time=1
RSS<-rep(0,length(lambda1_seq)*length(lambda2_seq))
loglik<-vector()
result<-matrix(0,length(lambda1_seq)*length(lambda2_seq),14)


###############################


for(qq in 1:length(lambda1_seq))
  { for(ww in 1:length(lambda2_seq))
    {
      lambda1<-lambda1_seq[qq] 
      lambda2<-lambda2_seq[ww]
    
      set<-list()
      for (k in 1: m) {set[[k]]=sample(nrow(X[[k]]),nrow(X[[k]]),replace = F)}
      rss<-0
     
      for (kk in 1: 5)
        {
          set.seed(kk)
          ##divide the data into training and test sets
          X_train=list();X_test=list();Y_train=list();Y_test=list()
          set[[k]]=sample(nrow(X[[k]]),nrow(X[[k]]),replace = F)
          for(k in 1:m)
            {
              start=(floor(length(set[[k]])/5)*(kk-1)+1)
              end=floor(length(set[[k]])/5)*kk  
              train<-set[[k]][-c(start:end)]
              X_train[[k]]<-X[[k]][train,];Y_train[[k]]<-as.matrix(Y[[k]][train,])
              X_test[[k]]<-X[[k]][-train,];Y_test[[k]]<-as.matrix(Y[[k]][-train,])
            }
    
          beta=magnitude_joint(X_train,Y_train,lambda1,lambda2,gamma)                                          ## for magnitude-based shrinkage
          #beta=magnitude_joint(X_train,Y_train,lambda1,lambda2,gamma)                                          ## for sign-based shrinkage
          for(k in 1:m)
            {
              Rss<-rep(0,nrow(X_test[[k]]))
              for(row in 1:nrow(X_test[[k]])){Rss[row]=(Y_test[[k]][row,1]-X_test[[k]][row,]%*%beta[,k])^2}
              Rss[which(Rss==0)]=min(Rss[which(Rss!=0)])
              rss=rss+sum(Rss)
            }
        }
    
      df=vector()
      result[time,1]=lambda1;result[time,2]=lambda2;result[time,3]=sum(beta!=0);result[time,4]=log(rss/N)
      sum=0;for(j in 1:p){if(sum(beta[j,]==0)==0){sum=sum+1}};result[time,5]<-sum
      for(k in 1:m){result[time,5+k]<-sum(beta[,k]!=0)}
      time=time+1
      print(k) 
    }
}


colnames(result)<-c("lambda1","lambda2","total-df","log_rss/N","overlap","data1_nonzero_coef","data2_nonzero_coef","data3_nonzero_coef","data4_nonzero_coef","data5_nonzero_coef","data6_nonzero_coef","data7_nonzero_coef","data8_nonzero_coef","data9_nonzero_coef")
result
