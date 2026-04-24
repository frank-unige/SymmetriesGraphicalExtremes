##load necessary packages
# install the latest version of graphical extremes with
## install.packages("devtools")
#devtools::install_github("sebastian-engelke/graphicalExtremes")
# Note that only the gitHub version contains the full flights data, without the code below will lead to errors

library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(latex2exp)
library(xtable)
library(cluster)
library(Matrix)
source("scoring_algorithm_functions.R")

# Filter and select observations without NAs

airports=which(rowSums(apply(flights$flightCounts[,,1:16], c(3),colSums)>2000)==16&rowSums(apply(flights$flightCounts[,,1:16], c(3),rowSums)>2000)==16)
d=length(airports)
graph_training_data=na.omit(flights$delays[1:1826,airports,"arrivals"]+flights$delays[1:1826,airports,"departures"])
color_training_data=na.omit(flights$delays[1827:3439,airports,"arrivals"]+flights$delays[1827:3439,airports,"departures"])
validation_data=na.omit(flights$delays[3440:5235,airports,"arrivals"]+flights$delays[3440:5235,airports,"departures"])

#Estimate full empirical variogram
p=0.85
vario_emp_full=emp_vario(graph_training_data,p=p)

#Clustering

nClusters <- 4
diss <- vario_emp_full

kmed <- pam(
  diss,
  diss = TRUE,
  k = nClusters,
  nstart = 10,
  variant = 'original'
)

############################################################
## Structure learning choice:emtp2 and eglearn comparison ##
############################################################

loglik_eg_best=list()
graph_eg_list=list()
loglik_emp_graph_training_data_list=list()
Gam_eg_best_list=list()

for (clust in 1:nClusters) {
  cluster=which(kmed$clustering==clust)
  airports_cluster = airports[cluster]
  Y= data2mpareto(graph_training_data[,cluster],p=p)
  n=nrow(Y)
  vario_emp=emp_vario(Y)
  loglik_emp_graph_training_data_list[[clust]]=loglik_HR(data=validation_data[,cluster],p=p, graph = make_full_graph(length(cluster)), Gamma=vario_emp, cens=FALSE)[1]
  
  rholist=seq(0, 0.2, length.out = 21)
  lasso.est = eglearn(graph_training_data[,which(kmed$clustering==clust)],p=p,rholist=rholist)
  loglik_eg <- list()
  Gam_eg_list=list()
  for (i in 1:length(rholist)){
    Gam_eg_list[[i]] = complete_Gamma(vario_emp,graph=lasso.est$graph[[i]],final_tol=1e-6)
    loglik_eg[[i]] <- loglik_HR(data=validation_data[,which(kmed$clustering==clust)],p=p, graph = lasso.est$graph[[i]],
                                Gamma = Gam_eg_list[[i]], cens = FALSE)[1]
  }
  
  lasso.best=which(unlist(loglik_eg)==max(unlist(loglik_eg)))
  loglik_eg_best[[clust]]=loglik_eg[[lasso.best]]
  graph_eg_list[[clust]]=lasso.est$graph[[lasso.best]]
  Gam_eg_best_list[[clust]]=Gam_eg_list[[lasso.best]]
  
}

############################
### Learning RVAR models ###
############################

Gam_colListList=list()
colors_list = list()
partiListList=list()
results_col_list = list()
results_col_start_list = list()

for (clust in 1:nClusters) {
  ## Build cluster and variogram
  cluster=which(kmed$clustering==clust)
  airports_cluster = airports[cluster]
  
  Y= data2mpareto(color_training_data[,cluster],p=p)
  n=nrow(Y)
  vario_emp=emp_vario(Y)

  
  # Choose coloring graph
  
  graph_coloring = graph_eg_list[[clust]]
  

  
  ###########################
  ### Clustering approach ###
  ###########################
  
  # Assign coloring for each k
  
  edges = ends(graph_coloring,E(graph_coloring))
  kList = 1:floor(ecount(graph_coloring)/2)
  colors = vector(mode = "list", length = length(kList))
  for (i in 1:length(kList)) {
    colors[[i]] <- pam(
      vario_emp[edges],
      diss = FALSE,
      k = kList[i],
      nstart = 10,
      variant = 'original'
    )$clustering
  }
  colors_list[[clust]]=colors
  partiList = vector(mode = "list", length = length(kList))
  for (i in 1:length(kList)) {
    parti=vector(mode = "list", length = kList[i])
    for (j in 1:kList[i]) {
      parti[[j]]=which(colors[[i]]==j)
    }
    partiList[[i]]=parti
  }
  partiListList[[clust]]=partiList
  
  # Apply reciprocal scoring algorithm
  # Computation can be sped up using the edmcr package when calculating the starting Gamma, by setting start= FALSE in nu2Gamma. As this is not supported in R 4.5, we provide here an alternative solution.
  
  Gam_colList= vector(mode = "list", length = length(kList))
  Gam_colstartList= vector(mode = "list", length = length(kList))
  for (i in length(kList):1) {
    start_value=start_recipr(vario=Gam_eg_best_list[[clust]],graph = graph_coloring,partition = partiList[[i]])
    if(i==length(kList)){Gam_colstartList[[i]]=nu2Gamma(nu=start_value,graph = graph_coloring,start = TRUE ,start_Gamma=Gam_eg_best_list[[clust]] ,partition = partiList[[i]], N=1000)}
    if(i==1){Gam_colstartList[[i]]=nu2Gamma(nu=start_value,graph = graph_coloring,start = TRUE ,start_Gamma=matrix(start_value,nrow = nrow(vario_emp),ncol = ncol(vario_emp)) ,partition = partiList[[i]], N=1000)}
    if(i<length(kList) & i>1){Gam_colstartList[[i]]=nu2Gamma(nu=start_value,graph = graph_coloring,start = TRUE ,start_Gamma= Gam_colstartList[[i+1]] ,partition = partiList[[i]], N=1000)}
    if(!is_valid_Gamma(Gam_colstartList[[i]])){print("Invalid starting point in cluster") 
      print(clust) 
      print("and color class") 
      print(i)}
    res=Scoring_alg_col_recipr(Qhat=Theta2Q(Gamma2Theta(Gam_eg_best_list[[clust]])),graph = graph_coloring,partition = partiList[[i]],start_recipr=start_value,start_Gamma=Gam_colstartList[[i]],verbose = FALSE,tol=1e-6)
    Gam_colList[[i]]=res$Gamma
    if(!is_valid_Gamma(Gam_colList[[i]])){print("Invalid estimate in cluster") 
      print(clust) 
      print("and color class") 
      print(i)}
    }
  
  # Evaluate validation log-likelihood
  
  results_col = vector(mode = "list", length = length(kList))
  results_col_start = vector(mode = "list", length = length(kList))
  for (i in 1:length(kList)) {
    results_col[[i]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_coloring,
                               Gamma = Gam_colList[[i]], cens = FALSE)[1]
    results_col_start[[i]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_coloring,
                               Gamma =  Gam_colstartList[[i]], cens = FALSE)[1]
  }
  
  results_col_list[[clust]]=results_col
  results_col_start_list[[clust]]=results_col_start
  Gam_colListList[[clust]]=Gam_colList
  
  print(clust)
}




###################################################
## Plotting log-likelihoods vs. number of colors ##
###################################################

pointshapes = c(16,17,18,3)
for (clust in 1:nClusters) {
  kList = 1:floor(ecount(graph_eg_list[[clust]])/2)
  pdf(paste0("plot_lik_RVAR_", clust,".pdf"),width = 7,height = 4.5)
  par(cex = 1.25, cex.lab = 1.3, cex.axis = 1, cex.main = 1.5,
      mar = c(4,4,3,2) +.1)
  matplot(cbind(kList,kList), cbind(unlist(results_col_list[[clust]]),unlist(results_col_start_list[[clust]])), type = "b",
          xlab = expression(paste("number of colors ", k )),
          ylab = "log-likelihood",
          ylim = c(min(min(unlist(results_col_list[[clust]])),min(unlist(results_col_start_list[[clust]])), loglik_eg_best[[clust]],loglik_emp_graph_training_data_list[[clust]]),
                   max(max(unlist(results_col_list[[clust]])),max(unlist(results_col_start_list[[clust]])), loglik_eg_best[[clust]],loglik_emp_graph_training_data_list[[clust]])),
          col=c('black'), pch=c(16,17), lwd=2)
  abline(loglik_eg_best[[clust]], b=0, lty=2, col="red", lwd=2,)
  abline(loglik_emp_graph_training_data_list[[clust]], b=0, lty=1, col="blue", lwd=2,)
  grid()
  dev.off()
}


#######################
## Result comparison ##
#######################

kbest_list = rep(0,nClusters)
for (clust in 1:nClusters) {
  kbest_list[clust]=which(unlist(results_col_list[[clust]])==max(unlist(results_col_list[[clust]])))
}

## Compare results of different estimators via test data
results =  matrix(nrow = 4, ncol = 3)
colnames(results) = c("Airports","nb par eglearn", "nb par best k")
rownames(results) = c("Southern","Western",  "Central","Eastern")
results[,1] = c(ncol(Gam_colListList[[1]][[kbest_list[1]]]),ncol(Gam_colListList[[2]][[kbest_list[2]]]),ncol(Gam_colListList[[3]][[kbest_list[3]]]),ncol(Gam_colListList[[4]][[kbest_list[4]]]))
results[,2] = c(ecount(graph_eg_list[[1]]),ecount(graph_eg_list[[2]]),ecount(graph_eg_list[[3]]),ecount(graph_eg_list[[4]]))
results[,3] = kbest_list

table = xtable(results,digits = 0)
print(table,type = "latex")  


##################
## Result plots ##
##################
maptype=c("state","world","state","world")
for (clust in 1:nClusters) {
  pdf(paste0("colored_graph_RVAR_", clust,".pdf"),width = 7,height = 4.5)
  colored_plot=plotFlights(names(airports[which(kmed$clustering==clust)]),graph=graph_eg_list[[clust]],map=maptype[clust],clipMap=1.2,edgeColors = as.character(colors_list[[clust]][[kbest_list[clust]]]),edgeAlpha = 1,returnGGPlot = TRUE)+ theme(legend.position = "none")
  plot(colored_plot)
  dev.off()
}



############################################
## Empirical correlation comparison plots ##
############################################

library(ggm)
source("plot_function.R")

for (clust in 1:nClusters) {
  pdf(paste0("correlation_plot_RVAR_", clust,".pdf"),width = 7,height = 4.5)
  G0 <- emp_vario(validation_data[,names(airports[which(kmed$clustering==clust)])],p=p)
  chi0 <- Gamma2chi(G0)
  chi1 <- Gamma2chi(Gam_colListList[[clust]][[kbest_list[clust]]])
  ggp <- plot_fitted_params(chi0, chi1, colors = makeColorMat(graph=graph_eg_list[[clust]], parti = partiListList[[clust]][[kbest_list[clust]]]), isFirst = TRUE)
  plot(ggp)
  dev.off()
}

