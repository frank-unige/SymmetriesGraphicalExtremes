##load necessary packages
# install the latest version of graphical extremes with
## install.packages("devtools")
#devtools::install_github("sebastian-engelke/graphicalExtremes")
# Note that only the gitHub version contains the full flights data

library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(latex2exp)
library(xtable)
library(cluster)
source("scoring_algorithm_functions.R")

# Filter and select observations without NAs

airports=which(rowSums(apply(flights$flightCounts[,,1:16], c(3),colSums)>2000)==16&rowSums(apply(flights$flightCounts[,,1:16], c(3),rowSums)>2000)==16)
d=length(airports)
training_data=na.omit(flights$delays[1:2191,airports,"arrivals"]+flights$delays[1:2191,airports,"departures"])
validation_data=na.omit(flights$delays[2192:5235,airports,"arrivals"]+flights$delays[2192:5235,airports,"departures"])

#Estimate full empirical variogram
plist=c(0.8,0.85,0.9)
kmedlist=list()
for (i in 1:3) {
vario_emp_full=emp_vario(training_data,p=plist[i])

#Clustering

nClusters <- 4
diss <- vario_emp_full

kmedlist[[i]] <- pam(
  diss,
  diss = TRUE,
  k = nClusters,
  nstart = 10,
  variant = 'original'
)
}

# Cluster comparison

kmedlist[[1]]$clustering-kmedlist[[2]]$clustering
kmedlist[[3]]$clustering-kmedlist[[2]]$clustering

## clustering is stable

kmed=kmedlist[[1]]

############################################
### RCON comparison for southern cluster ###
############################################

Gam_colListList=list()
colors_list = list()
partiListList=list()
results_col_list = list()
results_col_start_list = list()
loglik_emp_list=list()
graph_emtp2_list = list()
loglik_EMTP2_list=list()
for (i in 1:3) {
  p=plist[i]
  ## Build cluster and variogram
  cluster=which(kmed$clustering==1)
  airports_cluster = airports[cluster]
  
  Y= data2mpareto(training_data[,cluster],p=p)
  n=nrow(Y)
  vario_emp=emp_vario(Y)
  Theta_emp=Gamma2Theta(vario_emp)
  loglik_emp_list[[i]]=loglik_HR(data=validation_data[,cluster],p=p, graph = make_full_graph(length(cluster)), Gamma=vario_emp, cens=FALSE)[1]
  
  
  # EMTP2 estimate
  Gam_emtp2=emtp2(vario_emp,verbose = FALSE)$G_emtp2
  graph_emtp2=Gamma2graph(Gam_emtp2)
  graph_emtp2_list[[i]]=graph_emtp2
  loglik_EMTP2_list[[i]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_emtp2, Gamma=Gam_emtp2, cens=FALSE)[1]
  
  # Choose coloring graph
  
  graph_coloring = graph_emtp2
  
  ###########################
  ### Clustering approach ###
  ###########################
  
  # Assign coloring for each k 
  
  edges = ends(graph_coloring,E(graph_coloring))
  kList = 1:25
  colors = vector(mode = "list", length = length(kList))
  for (l in 1:length(kList)) {
    colors[[l]] <- pam(
      Theta_emp[edges],
      diss = FALSE,
      k = kList[l],
      nstart = 10,
      variant = 'original'
    )$clustering
  }
  colors_list[[i]]=colors
  partiList = vector(mode = "list", length = length(kList))
  for (l in 1:length(kList)) {
    parti=vector(mode = "list", length = kList[l])
    for (j in 1:kList[l]) {
      parti[[j]]=which(colors[[l]]==j)
    }
    partiList[[l]]=parti
  }
  partiListList[[i]]=partiList
  
  # Apply scoring algorithm
  
  Gam_colList= vector(mode = "list", length = length(kList))
  Gam_colstartList= vector(mode = "list", length = length(kList))
  for (l in 1:length(kList)) {
    start_value=start(vario = vario_emp,graph = graph_coloring,partition = partiList[[l]])
    res=Scoring_alg_col(vario=vario_emp,graph = graph_coloring,partition = partiList[[l]],start=start_value)
    Gam_colList[[l]]= Theta2Gamma(Q2Theta(omega2Q(res$omega,graph = graph_coloring,partition = partiList[[l]]))) 
    Gam_colstartList[[l]]= Theta2Gamma(Q2Theta(omega2Q(start_value,graph = graph_coloring,partition = partiList[[l]])))
  }
  
  # Evaluate validation log-likelihood
  
  results_col = vector(mode = "list", length = length(kList))
  results_col_start = vector(mode = "list", length = length(kList))
  for (l in 1:length(kList)) {
    results_col[[l]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_coloring,
                               Gamma = Gam_colList[[l]], cens = FALSE)[1]
    results_col_start[[l]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_coloring,
                                     Gamma = Gam_colstartList[[l]], cens = FALSE)[1]
  }
  
  results_col_list[[i]]=results_col
  results_col_start_list[[i]]=results_col_start
  Gam_colListList[[i]]=Gam_colList
}

###################################################
## Plotting log-likelihoods vs. number of colors ##
###################################################

for (i in 1:3) {
  pdf(paste0("plot_lik_sensitivity_", i,".pdf"),width = 7,height = 4.5)
  par(cex = 1.25, cex.lab = 1.3, cex.axis = 1, cex.main = 1.5,
      mar = c(4,4,3,2) +.1)
  matplot(cbind(kList,kList), cbind(unlist(results_col_list[[i]]),unlist(results_col_start_list[[i]])), type = "b",
          xlab = expression(paste("number of colors ", k )),
          ylab = "log-likelihood", #main = expression(paste("Log-likelihood vs. ", k )),
          ylim = c(min(min(unlist(results_col_list[[i]])),min(unlist(results_col_start_list[[i]])), loglik_EMTP2_list[[i]],loglik_emp_list[[i]]),
                   max(max(unlist(results_col_list[[i]])),max(unlist(results_col_start_list[[i]])), loglik_EMTP2_list[[i]],loglik_emp_list[[i]])),
          col=c('black'), pch=c(16,17), lwd=2)
  abline(loglik_EMTP2_list[[i]], b=0, lty=2, col="orange", lwd=2,)
  abline(loglik_emp_list[[i]], b=0, lty=1, col="blue", lwd=2,)
  grid()
  dev.off()
}

kbest_list = rep(0,3)
for (l in 1:3) {
  kbest_list[l]=which(unlist(results_col_list[[l]])==max(unlist(results_col_list[[l]])))
}

for (l in 1:3) {
  pdf(paste0("colored_graph_sensitivity", l,".pdf"),width = 7,height = 4.5)
  colored_plot=plotFlights(names(airports[which(kmed$clustering==1)]),graph=graph_emtp2_list[[l]],map="state",clipMap=1.2,edgeColors = as.character(colors_list[[l]][[kbest_list[l]]]),edgeAlpha = 1,returnGGPlot = TRUE)+ theme(legend.position = "none")
  plot(colored_plot)
  dev.off()
}




###################################
### RVAR sensitivity comparison ###
###################################

Gam_colListList=list()
colors_list = list()
partiListList=list()
results_col_list = list()
results_col_start_list = list()
loglik_emp_list=list()
graph_emtp2_list = list()
loglik_EMTP2_list=list()
for (l in 1:3) {
  p=plist[l]
  ## Build cluster and variogram
  cluster=which(kmed$clustering==1)
  airports_cluster = airports[cluster]
  
  Y= data2mpareto(training_data[,cluster],p=p)
  n=nrow(Y)
  vario_emp=emp_vario(Y)
  loglik_emp_list[[l]]=loglik_HR(data=validation_data[,cluster],p=p, graph = make_full_graph(length(cluster)), Gamma=vario_emp, cens=FALSE)[1]
  
  
  # EMTP2 estimate
  Gam_emtp2=emtp2(vario_emp,verbose = FALSE)$G_emtp2
  Q_emtp2= Theta2Q(Gamma2Theta(Gam_emtp2))
  graph_emtp2=Gamma2graph(Gam_emtp2)
  graph_emtp2_list[[l]]=graph_emtp2
  loglik_EMTP2_list[[l]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_emtp2, Gamma=Gam_emtp2, cens=FALSE)[1]
  
  # Choose coloring graph
  
  graph_coloring = graph_emtp2
  
  ###########################
  ### Clustering approach ###
  ###########################
  
  # Assign coloring for each k
  
  edges = ends(graph_coloring,E(graph_coloring))
  kList = 1:25
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
  colors_list[[l]]=colors
  partiList = vector(mode = "list", length = length(kList))
  for (i in 1:length(kList)) {
    parti=vector(mode = "list", length = kList[i])
    for (j in 1:kList[i]) {
      parti[[j]]=which(colors[[i]]==j)
    }
    partiList[[i]]=parti
  }
  partiListList[[l]]=partiList
  
  # Apply reciprocal scoring algorithm
  
  Gam_colList= vector(mode = "list", length = length(kList))
  Gam_colstartList= vector(mode = "list", length = length(kList))
  for (i in 1:length(kList)) {
    start_value=start_recipr(vario=Gam_emtp2,graph = graph_coloring,partition = partiList[[i]])
    res=Scoring_alg_col_recipr(Qhat=Q_emtp2,graph = graph_coloring,partition = partiList[[i]],start_recipr=start_value,verbose = FALSE,tol=1e-6)
    Gam_colList[[i]]=res$Gamma   
    Gam_colstartList[[i]]=nu2Gamma(nu=start_value,graph = graph_coloring,partition = partiList[[i]], N=1000)
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
  
  results_col_list[[l]]=results_col
  results_col_start_list[[l]]=results_col_start
  Gam_colListList[[l]]=Gam_colList
}

###############
## Plotting  ##
###############

for (l in 1:3) {
  pdf(paste0("plot_lik_RVAR_sensitivity_", l,".pdf"),width = 7,height = 4.5)
  par(cex = 1.25, cex.lab = 1.3, cex.axis = 1, cex.main = 1.5,
      mar = c(4,4,3,2) +.1)
  matplot(cbind(kList,kList), cbind(unlist(results_col_list[[l]]),unlist(results_col_start_list[[l]])), type = "b",
          xlab = expression(paste("number of colors ", k )),
          ylab = "log-likelihood", #main = expression(paste("Log-likelihood vs. ", k )),
          ylim = c(min(min(unlist(results_col_list[[l]])),min(unlist(results_col_start_list[[l]])), loglik_EMTP2_list[[l]],loglik_emp_list[[l]]),
                   max(max(unlist(results_col_list[[l]])),max(unlist(results_col_start_list[[l]])), loglik_EMTP2_list[[l]],loglik_emp_list[[l]])),
          col=c('black'), pch=c(16,17), lwd=2)
  abline(loglik_EMTP2_list[[l]], b=0, lty=2, col="orange", lwd=2,)
  abline(loglik_emp_list[[l]], b=0, lty=1, col="blue", lwd=2,)
  grid()
  dev.off()
}


kbest_list = rep(0,3)
for (l in 1:3) {
  kbest_list[l]=which(unlist(results_col_list[[l]])==max(unlist(results_col_list[[l]])))
}

for (l in 1:3) {
  pdf(paste0("colored_graph_RVAR_sensitivity", l,".pdf"),width = 7,height = 4.5)
  colored_plot=plotFlights(names(airports[which(kmed$clustering==1)]),graph=graph_emtp2_list[[l]],map="state",clipMap=1.2,edgeColors = as.character(colors_list[[l]][[kbest_list[l]]]),edgeAlpha = 1,returnGGPlot = TRUE)+ theme(legend.position = "none")
  plot(colored_plot)
  dev.off()
}


######################
## graph comparison ##
######################

graph_emtp2_list[[1]]-graph_emtp2_list[[2]]
graph_emtp2_list[[2]]-graph_emtp2_list[[1]]
graph_emtp2_list[[3]]-graph_emtp2_list[[2]]
graph_emtp2_list[[2]]-graph_emtp2_list[[3]]


