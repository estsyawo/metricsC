# Set working directory to source file
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
list.files()
#set directory to source file



MGS<- function(A,m,n)
{ Q = matrix(NA,m,n)
  R = matrix(0,n,n)
  for (k in 1:n) {
    R[k,k] = sqrt(sum(A[,k]^2))
    Q[,k] = A[,k]/R[k,k]
    
    if(k<n)
    {for (j in (k+1):n) {
      R[k,j] = t(Q[,k])%*%A[,j]
      A[,j] = A[,j] - Q[,k]*R[k,j]
    }
      }
  }
  list(Q=Q,R=R)
}

# Test:
dat<- read.table("data.txt",header = F)

X = dat[,-11]; Y = dat[,11]

(msObj<-MGS(X,nrow(X),ncol(X)))
