## Function to transform Gamma from matrix in vector form (I think lexicographic)

Gamma2vec <- function(Gamma){
  m = ncol(Gamma)
  G <- vector(length = choose(m,2))
  G<- Gamma[lower.tri(Gamma, diag=FALSE)]  
  return(G)
}

## Function to transform a Gamma vector back to matrix form

Gvec2Gamma <-function(Gvec){
  n=length(Gvec)
  m=sqrt(2*n+1/4)+1/2
  G = matrix(0,nrow = m,ncol = m)
  G[lower.tri(G, diag=FALSE)] <- Gvec
  return(G+t(G))
}

## Function to convert Theta matrix to Q matrix

Theta2Q <- function(Theta){
  Q=-Theta
  diag(Q)=rep(0,nrow(Theta))
  return(Q)
}

## Function to convert Q matrix to Theta matrix

Q2Theta <- function(Q){
  Theta=-Q
  diag(Theta)=rowSums(Q)
  return(Theta)
}

## Function to compute Cayley--Menger matrix for a given Gamma

CM = function(Gam){
  m=ncol(Gam)
  CM = rbind(cbind(-Gam/2,rep(1,m)),c(rep(-1,m),0))
  return(CM)}

## Score function for a given Gamma, sample variogram in a graphical model

scores = function(Gam,vario,graph){
  edges = ends(graph,E(graph))
  scores = (Gam[edges]-vario[edges])/2
  return(scores)
}

## Function to compute the Information of a given Gamma for an uncolored graphical model

Information = function(Gam,graph){
  edges = ends(graph,E(graph))
  dim = nrow(edges)
  Information = matrix(0,ncol = dim,nrow = dim)
  for (i in 1:dim) {
    for (j in 1:dim) {
      Information[i,j]=(1/8)*(Gam[edges[i,1],edges[j,1]]-Gam[edges[i,1],edges[j,2]]
                              -Gam[edges[i,2],edges[j,1]]+Gam[edges[j,2],edges[i,2]])^2
    }
  }
  return(Information)
}

## Score function for a colored model

scoresCol = function(omega,vario,graph,partition){
  Q=omega2Q(omega=omega,graph = graph,partition = partition)
  Gamma=Theta2Gamma(Q2Theta(Q))  
  scores = scores(Gam=Gamma,vario = vario,graph = graph)
  dim2=length(omega)
  scoresCol=rep(0,dim2)
  for (i in 1:dim2) {
    scoresCol[i]=sum(scores[partition[[i]]])
  }
  return(scoresCol)
}

## Information matrix in a colored graphical model with parameter omega. Coloring is determined by a partition.

InformationCol = function(omega,graph,partition){
  Q=omega2Q(omega=omega,graph = graph,partition = partition)
  Gamma=Theta2Gamma(Q2Theta(Q))
  Information=Information(Gam=Gamma,graph = graph)
  dim2=length(omega)
  InformationCol = matrix(0,ncol=dim2,nrow=dim2)
  for (i in 1:(dim2)) {
    for (j in 1:(dim2)) {
      InformationCol[i,j]= sum(Information[partition[[i]],partition[[j]]])
    }
  }
  return(InformationCol)
}

## Function converts the parameter omega of a colored model to a matrix Q

omega2Q = function(omega,graph,partition){
  edges=ends(graph,E(graph))
  dim1=length(V(graph))
  dim2=length(omega)
  Q=matrix(0,nrow = dim1,ncol=dim1)
  for (i in 1:dim2) {
    Q[edges][partition[[i]]]=rep(omega[i],length(partition[[i]]))
  }
  Q[lower.tri(Q,diag = FALSE)]=t(Q)[lower.tri(Q,diag = FALSE)]
  return(Q)
}

## Function to compute the starting point for the scoring algorithm. It averages over the empirical weights in each color class.

start = function(vario,graph,partition){
  edges=ends(graph,E(graph))
  Gam=complete_Gamma(Gamma = vario,graph = graph)
  Q=Theta2Q(Gamma2Theta(Gam))
  dim2=length(partition)
  omega=rep(0,dim2)
  for (i in 1:dim2) {
    if(length(edges[partition[[i]],])==2){
      omega[i]=mean(Q[edges[partition[[i]],1],edges[partition[[i]],2]])}
    else{
      omega[i]=mean(Q[edges[partition[[i]],]])}
    }
  return(omega)
}

## Scoring algorithm for our setting.

Scoring_alg_col = function(vario,graph,partition,start,tol = 1e-10,maxiter=1000,verbose=FALSE){
  edges=ends(graph,E(graph))
  omega=start
  S=scoresCol(omega=omega, vario = vario, graph = graph,partition=partition)  
  gap=Inf
  iter=0
  while (gap>tol&&iter<maxiter) {
    omega=omega+ solve(InformationCol(omega=omega,graph = graph,partition=partition)+S%*%t(S))%*%S
    S=scoresCol(omega=omega, vario = vario, graph = graph,partition=partition)  
    gap=sum(abs(S))
    iter=iter+1
    if (verbose == TRUE) {
      cat(iter, "\t  | ", gap, "\n")
    }
  }  
  return(list(omega=omega,gap=gap,iter=iter))
}

## Reciprocal Score function for a given Q (sample adjacency in a graphical model)

scores_recipr = function(Q,Qhat,graph){
  edges = ends(graph,E(graph))
  scores_recipr = (Q[edges]-Qhat[edges])/2
  return(scores_recipr)
}

## Reciprocal Information

Information_recipr = function(Q,graph){
  edges = ends(graph,E(graph))
  dim = nrow(edges)
  Theta= Q2Theta(Q)
  Information_recipr = matrix(0,ncol = dim,nrow = dim)
  for (i in 1:dim) {
    for (j in 1:dim) {
      Information_recipr[i,j]=(1/4)*(Theta[edges[i,1],edges[j,2]]*Theta[edges[i,2],edges[j,1]]+Theta[edges[i,1],edges[j,1]]*Theta[edges[i,2],edges[j,2]])
    }
  }
  return(Information_recipr)
}

## Function converts the parameter omega of a colored model to a matrix Q

nu2Gamma = function(nu,graph,partition,start=FALSE,start_Gamma=NA,N=20000){
  edges=ends(graph,E(graph))
  dim1=length(V(graph))
  dim2=length(nu)
  if(start==TRUE){Gamma=start_Gamma}
  else{Gamma=matrix(NA,nrow = dim1,ncol=dim1)}
  for (i in 1:dim2) {
    Gamma[edges][partition[[i]]]=rep(nu[i],length(partition[[i]]))
  }
  Gamma[lower.tri(Gamma,diag = FALSE)]=t(Gamma)[lower.tri(Gamma,diag = FALSE)]
  diag(Gamma)=rep(0,dim1) 
  Gamma=complete_Gamma(Gamma,graph = graph,N=N)
  return(Gamma)
}

## Reciprocal score function for a colored model

scoresCol_recipr = function(Gamma,Qhat,graph,partition){
  Q=Theta2Q(Gamma2Theta(Gamma)) 
  scores_recipr = scores_recipr(Q=Q,Qhat = Qhat,graph = graph)
  dim2=length(partition)
  scoresCol_recipr=rep(0,dim2)
  for (i in 1:dim2) {
    scoresCol_recipr[i]=sum(scores_recipr[partition[[i]]])
  }
  return(scoresCol_recipr)
}

## Reciprocal information matrix in a colored graphical model. Coloring is determined by a partition.

InformationCol_recipr = function(Gamma,graph,partition){
  Q=Theta2Q(Gamma2Theta(Gamma)) 
  Information_recipr=Information_recipr(Q=Q,graph = graph)
  dim2=length(partition)
  InformationCol_recipr = matrix(0,ncol=dim2,nrow=dim2)
  for (i in 1:(dim2)) {
    for (j in 1:(dim2)) {
      InformationCol_recipr[i,j]= sum(Information_recipr[partition[[i]],partition[[j]]])
    }
  }
  return(InformationCol_recipr)
}

## Function to compute the starting point for the reciorocal scoring algorithm. It averages over the empirical weights in each color class.

start_recipr = function(Qhat,graph,partition){
  edges=ends(graph,E(graph))
  Gamma= Theta2Gamma(Q2Theta(Qhat))
  dim2=length(partition)
  nu=rep(0,dim2)
  for (i in 1:dim2) {
    if(length(edges[partition[[i]],])==2){
      nu[i]=mean(Gamma[edges[partition[[i]],1],edges[partition[[i]],2]])}
    else{
      nu[i]=mean(Gamma[edges[partition[[i]],]])}  }
  return(nu)
}

## Reciprocal scoring algorithm.

Scoring_alg_col_recipr = function(Qhat,graph,partition,start_recipr,tol = 1e-10,maxiter=1000,verbose=FALSE){
  edges=ends(graph,E(graph))
  nu=start_recipr
  Gamma= nu2Gamma(nu,graph,partition,N=1000)
  S=scoresCol_recipr(Gamma=Gamma, Qhat = Qhat, graph = graph,partition=partition)  
  gap=Inf
  iter=0
  while (gap>tol&&iter<maxiter) {
    nu=nu+solve(InformationCol_recipr(Gamma=Gamma,graph = graph,partition=partition)+S%*%t(S))%*%S
    Gamma= nu2Gamma(nu,graph,partition,start = TRUE,start_Gamma=Gamma)
    S=scoresCol_recipr(Gamma=Gamma, Qhat = Qhat, graph = graph,partition=partition)  
    gap=sum(abs(S))
    iter=iter+1
    if (verbose == TRUE) {
      cat(iter, "\t  | ", gap, "\n")
    } 
  } 
  return(list(nu=nu,Gamma=Gamma,scores=S,gap=gap,iter=iter))
}
