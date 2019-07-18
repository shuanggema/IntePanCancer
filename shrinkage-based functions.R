########################################### magrginal analysis #####################################

  ## define magnitude-based shrinkage for marginal analysis
  magnitude_marginal<-function(X,Y,lambda1,lambda2,gamma,specific)
    {  change<-vector()
       it = 0
       maxit = 200
       beta<-matrix(0,1,m); beta<-apply(beta,2,as.numeric)
       a = 0
       b = 0
       eps=10
       while ((it<=200) & (eps >10^(-3)) )
         { it<-it+1
           Lold = beta
           for ( k in 1:m) 
             {  x = X[[k]][,specific]
                y = Y[[k]]
                b3=0;b4=0 
                for (r in 1:m){ if(r != k & (sign(S[specific,k])==sign(S[specific,r]))){b3 <-b3+ beta[r];b4<-b4+1} }
                a = t(y)%*%x /nrow(X[[k]])+ 2.0*lambda2*b3
                b = t(x)%*%x/nrow(X[[k]]) + lambda2* (m-1+b4)
                temp=a/b
                if (abs(temp)>gamma*lambda1){ beta[k]=temp } 
                else 
                  { if (abs(a)>lambda1){ beta[k]=(a-sign(a)*lambda1)/(b-1/gamma)} 
                    else { beta[k]=0 }
                  }
             }
           Lnew<-beta
           change [it]<- sum(abs(Lold-Lnew)) 
           eps<-change[it]
         } 
       
       return(beta=beta)
   }



  ## define sign-based shrinkage for marginal analysis
  sign_marginal<-function(X,Y,lambda1,lambda2,gamma,specific)
    { change<-vector()
      it = 0
      maxit = 200
      beta<-matrix(0,1,m);beta<-apply(beta,2,as.numeric)
      a = 0
      b = 0
      eps=10
      while ((it<=200) & (eps >10^(-3)) )
        { it<-it+1
          Lold = beta
          for ( k in 1:m) 
            { x = X[[k]][,specific]
              y = Y[[k]]
              a5=0
              for (r in 1:m) {if(r != k )  {a5 <-a5+ lambda2*sign(beta[r])/((abs(beta[k])+xi)) } }
              a =t(y)%*%x /nrow(X[[k]])+ a5
              b = t(x)%*%x/nrow(X[[k]])+lambda2*(m-1)/(abs(beta[k])+xi)^2
              temp=a/b
              if (abs(temp)>gamma*lambda1){ beta[k]=temp } 
              else 
                { if (abs(a)>lambda1){ beta[k]=(a-sign(a)*lambda1)/(b-1/gamma)} 
                  else { beta[k]=0 }
                 }
             }
           Lnew<-beta
           change [it]<- sum(abs(Lold-Lnew)) 
           eps<-change[it]
        }
      return(beta=beta)
   }

  
  
########################################### joint analysis #####################################
  
  ## define magnitude-based shrinkage for joint analysis
  magnitude_joint<-function(X,Y,lambda1,lambda2,gamma)
    {
      change<-vector()
      it = 0
      maxit = 200
      a = 0;b = 0;eps=10
      beta<-matrix(0,p,m);beta<-apply(beta,2,as.numeric)
      while ((it<=200) & (eps >10^(-3)) )
        { it<-it+1
          Lold = beta
          for( k in 1:m)
             {  x = X[[k]];  y = Y[[k]]
                res_l=y - x %*% beta[,k]        
                for( l in 1:p)
                   { 
                      res_l = res_l+x[,l]*beta[l,k]
                      b3=0;b4=0 
                      for (r in 1:m){ if(r != k & (sign(S[l,k])==sign(S[l,r]))){b3 <-b3+ beta[l, r];b4<-b4+1} }
                      a = t(res_l)%*%x[,l] /nrow(x)+ 2.0*lambda2*b3
                      b = t(x[,l])%*%x[,l ]/nrow(x) + lambda2* (m-1+b4)
                      temp=a/b 
                      if (abs(temp)>gamma*lambda1){  beta[l,k]=temp } 
                      else{
                            if (abs(a)>lambda1){ beta[l,k]=(a-sign(a)*lambda1)/(b-1/gamma)} 
                            else{beta[l,k]=0}
                          }
                      res_l = res_l - x[,l]*beta[l,k]
                    }   
              }
          Lnew<-beta
          change [it]<- sum(abs(Lold-Lnew)) 
          eps<-change[it]
        } 
      return(beta=beta)
    }
  
  
  ## define sign-based shrinkage for joint analysis 
  sign_joint<-function(X,Y,lambda1,lambda2,gamma)
    {
      change<-vector()
      it = 0
      maxit = 200
      m = length(X)
      beta<-matrix(0,p,m)
      beta<-apply(beta,2,as.numeric)
      a = 0
      b = 0
      eps=10
      while ((it<=200) & (eps >10^(-3)))
        { 
          it<-it+1
          Lold = beta
          for ( k in 1:m)
            {
              x = X[[k]]
              y = Y[[k]]
              res_l=y - x%*%beta[,k]     
              for( l in 1:p)
                { res_l = res_l+x[,l]*beta[l,k]
                  a5=0
                  for (r in 1:m) {if(r != k )  {a5 <-a5+ lambda2*sign(beta[l,r])/((abs(beta[l,k])+xi)) } }
                  a = t(res_l)%*%x[,l] /nrow(x)+ a5
                  b = t(x[,l])%*%x[,l]/nrow(x) +lambda2*(m-1)/(abs(beta[l,k])+xi)^2
                  temp=a/b
                  if (abs(temp)>gamma*lambda1){ beta[l,k]=temp } 
                  else
                    {
                      if (abs(a)>lambda1){ beta[l,k]=(a-sign(a)*lambda1)/(b-1/gamma)} 
                      else{beta[l,k]=0}
                    }
                  res_l = res_l - x[,l]*beta[l,k]
                } 
              }
          Lnew<-beta
          change [it]<- sum(abs(Lold-Lnew)) 
          eps<-change[it]
        }
      return(beta=beta)
    }
  
  