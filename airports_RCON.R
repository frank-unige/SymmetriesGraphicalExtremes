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
p=0.85
vario_emp_full=emp_vario(training_data,p=p)
#chi_emp=emp_chi(training_data,p=p)
#vario_emp_full=chi2Gamma(chi_emp)

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


flightsPerConnection <- apply(flights$flightCounts[, , 1:6], c(1, 2), sum)

# Compute undirected flights per connection
flightsPerConnectionUD <- flightsPerConnection + t(flightsPerConnection)
# Consider only connections between selected airports
flightsPerConnectionUD1 <- flightsPerConnectionUD[names(airports[which(kmed$clustering==1)]), names(airports[which(kmed$clustering==1)])]
flightsPerConnectionUD2 <- flightsPerConnectionUD[names(airports[which(kmed$clustering==2)]), names(airports[which(kmed$clustering==2)])]
flightsPerConnectionUD3 <- flightsPerConnectionUD[names(airports[which(kmed$clustering==3)]), names(airports[which(kmed$clustering==3)])]
flightsPerConnectionUD4 <- flightsPerConnectionUD[names(airports[which(kmed$clustering==4)]), names(airports[which(kmed$clustering==4)])]


# Make flight graph for connections with on average at least monthly flights
A= matrix(0,nrow = d, ncol = d)
minNConnections <- 72
A[which(kmed$clustering==1),which(kmed$clustering==1)] <- 1 * (flightsPerConnectionUD1 > minNConnections)
A[which(kmed$clustering==2),which(kmed$clustering==2)] <- 1 * (flightsPerConnectionUD2 > minNConnections)
A[which(kmed$clustering==3),which(kmed$clustering==3)] <- 1 * (flightsPerConnectionUD3 > minNConnections)
A[which(kmed$clustering==4),which(kmed$clustering==4)] <- 1 * (flightsPerConnectionUD4 > minNConnections)



flight_graph <- graph_from_adjacency_matrix(A, diag = FALSE, mode = "undirected",add.colnames = NA)
vertexShapes = kmed$clustering
names(vertexShapes)=names(airports)
clusterplot=plotFlights(names(airports), map="world",graph = flight_graph,clipMap = 1.2,vertexShapes = vertexShapes,returnGGPlot = TRUE)+ theme(legend.position = "none")
pdf("airports_cluster_map.pdf",width = 7,height = 4.5)
plot(clusterplot)
dev.off()

############################################################
## Structure learning choice:emtp2 and eglearn comparison ##
############################################################

loglik_eg_best=list()
graph_emtp2_list = list()
loglik_EMTP2_list=list()

for (clust in 1:nClusters) {
  cluster=which(kmed$clustering==clust)
  airports_cluster = airports[cluster]
  Y= data2mpareto(training_data[,cluster],p=p)
  n=nrow(Y)
  vario_emp=emp_vario(Y)
  rholist=seq(0, 0.2, length.out = 21)
  lasso.est = eglearn(training_data[,which(kmed$clustering==clust)],p=p,rholist=rholist)
  loglik_eg <- list()
  for (i in 1:length(rholist)){
    Gamma_eg = complete_Gamma(vario_emp,graph=lasso.est$graph[[i]],final_tol=1e-6)
    loglik_eg[[i]] <- loglik_HR(data=validation_data[,which(kmed$clustering==clust)],p=p, graph = lasso.est$graph[[i]],
                                Gamma = Gamma_eg, cens = FALSE)[1]
  }
  
  lasso.best=which(unlist(loglik_eg)==max(unlist(loglik_eg)))
  loglik_eg_best[[clust]]=loglik_eg[[lasso.best]]
  
  # EMTP2 estimate
  Gam_emtp2=emtp2(vario_emp,verbose = FALSE)$G_emtp2
  graph_emtp2=Gamma2graph(Gam_emtp2)
  graph_emtp2_list[[clust]]=graph_emtp2
  loglik_EMTP2_list[[clust]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_emtp2, Gamma=Gam_emtp2, cens=FALSE)[1]
  
}

#Table comparison

results =  matrix(nrow = 2, ncol = 4)
rownames(results) = c("EMTP2", "eglearn")
colnames(results) = c("Southern","Western",  "Central","Eastern")
results[1,] = unlist(loglik_EMTP2_list)
results[2,] = unlist(loglik_eg_best)

table = xtable(results,digits = 0)
print(table,type = "latex")  


############################
### Learning RCON models ###
############################

Gam_colListList=list()
colors_list = list()
partiListList=list()
results_col_list = list()
results_col_start_list = list()
loglik_emp_list=list()
for (clust in 1:nClusters) {
## Build cluster and variogram
cluster=which(kmed$clustering==clust)
airports_cluster = airports[cluster]

Y= data2mpareto(training_data[,cluster],p=p)
n=nrow(Y)
vario_emp=emp_vario(Y)
Theta_emp=Gamma2Theta(vario_emp)
loglik_emp_list[[clust]]=loglik_HR(data=validation_data[,cluster],p=p, graph = make_full_graph(length(cluster)), Gamma=vario_emp, cens=FALSE)[1]


# Choose coloring graph

graph_coloring = graph_emtp2_list[[clust]]

###########################
### Clustering approach ###
###########################

# Assign coloring for each k 

edges = ends(graph_coloring,E(graph_coloring))
kList = 1:25
colors = vector(mode = "list", length = length(kList))
for (i in 1:length(kList)) {
  colors[[i]] <- pam(
    Theta_emp[edges],
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

# Apply scoring algorithm

Gam_colList= vector(mode = "list", length = length(kList))
Gam_colstartList= vector(mode = "list", length = length(kList))
for (i in 1:length(kList)) {
  start_value=start(vario = vario_emp,graph = graph_coloring,partition = partiList[[i]])
  res=Scoring_alg_col(vario=vario_emp,graph = graph_coloring,partition = partiList[[i]],start=start_value)
  Gam_colList[[i]]= Theta2Gamma(Q2Theta(omega2Q(res$omega,graph = graph_coloring,partition = partiList[[i]]))) 
  Gam_colstartList[[i]]= Theta2Gamma(Q2Theta(omega2Q(start_value,graph = graph_coloring,partition = partiList[[i]])))
}

# Evaluate validation log-likelihood

results_col = vector(mode = "list", length = length(kList))
results_col_start = vector(mode = "list", length = length(kList))
for (i in 1:length(kList)) {
  results_col[[i]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_coloring,
                             Gamma = Gam_colList[[i]], cens = FALSE)[1]
  results_col_start[[i]]=loglik_HR(data=validation_data[,cluster],p=p, graph = graph_coloring,
                             Gamma = Gam_colstartList[[i]], cens = FALSE)[1]
}

results_col_list[[clust]]=results_col
results_col_start_list[[clust]]=results_col_start
Gam_colListList[[clust]]=Gam_colList
}



###################################################
## Plotting log-likelihoods vs. number of colors ##
###################################################

pointshapes = c(16,17,18,3)
for (clust in 1:nClusters) {
pdf(paste0("plot_lik_", clust,".pdf"),width = 7,height = 4.5)
par(cex = 1.25, cex.lab = 1.3, cex.axis = 1, cex.main = 1.5,
    mar = c(4,4,3,2) +.1)
matplot(cbind(kList,kList), cbind(unlist(results_col_list[[clust]]),unlist(results_col_start_list[[clust]])), type = "b",
        xlab = expression(paste("number of colors ", k )),
        ylab = "log-likelihood", #main = expression(paste("Log-likelihood vs. ", k )),
        ylim = c(min(min(unlist(results_col_list[[clust]])),min(unlist(results_col_start_list[[clust]])), loglik_EMTP2_list[[clust]],loglik_emp_list[[clust]]),
                 max(max(unlist(results_col_list[[clust]])),max(unlist(results_col_start_list[[clust]])), loglik_EMTP2_list[[clust]],loglik_emp_list[[clust]])),
        col=c('black'), pch=c(16,17), lwd=2)
abline(loglik_EMTP2_list[[clust]], b=0, lty=2, col="orange", lwd=2,)
abline(loglik_emp_list[[clust]], b=0, lty=1, col="blue", lwd=2,)
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
colnames(results) = c("Airports","nb par EMTP2", "nb par best k")
rownames(results) = c("Southern","Western",  "Central","Eastern")
results[,1] = c(ncol(Gam_colListList[[1]][[kbest_list[1]]]),ncol(Gam_colListList[[2]][[kbest_list[2]]]),ncol(Gam_colListList[[3]][[kbest_list[3]]]),ncol(Gam_colListList[[4]][[kbest_list[4]]]))
results[,2] = c(ecount(graph_emtp2_list[[1]]),ecount(graph_emtp2_list[[2]]),ecount(graph_emtp2_list[[3]]),ecount(graph_emtp2_list[[4]]))
results[,3] = kbest_list

table = xtable(results,digits = 0)
print(table,type = "latex")  

##################
## Result plots ##
##################
maptype=c("state","world","state","world")
for (clust in 1:nClusters) {
pdf(paste0("colored_graph_", clust,".pdf"),width = 7,height = 4.5)
colored_plot=plotFlights(names(airports[which(kmed$clustering==clust)]),graph=graph_emtp2_list[[clust]],map=maptype[clust],clipMap=1.2,edgeColors = as.character(colors_list[[clust]][[kbest_list[clust]]]),edgeAlpha = 1,returnGGPlot = TRUE)+ theme(legend.position = "none")
plot(colored_plot)
dev.off()
}

## Western cluster for k=4
pdf(paste0("colored_graph_", 5,".pdf"),width = 7,height = 4.5)
colored_plot=plotFlights(names(airports[which(kmed$clustering==2)]),graph=graph_emtp2_list[[2]],map=maptype[2],clipMap=1.2,edgeColors = as.character(colors_list[[2]][[4]]),edgeAlpha = 1,returnGGPlot = TRUE)+ theme(legend.position = "none")
plot(colored_plot)
dev.off()

## Central cluster for k=3
pdf(paste0("colored_graph_", 6,".pdf"),width = 7,height = 4.5)
colored_plot=plotFlights(names(airports[which(kmed$clustering==3)]),graph=graph_emtp2_list[[3]],map=maptype[3],clipMap=1.2,edgeColors = as.character(colors_list[[3]][[3]]),edgeAlpha = 1,returnGGPlot = TRUE)+ theme(legend.position = "none")
plot(colored_plot)
dev.off()

############################################
## Empirical correlation comparison plots ##
############################################

library(ggm)
source("plot_function.R")

for (clust in 1:nClusters) {
  pdf(paste0("correlation_plot_RCON_", clust,".pdf"),width = 7,height = 4.5)
  G0 <- emp_vario(validation_data[,names(airports[which(kmed$clustering==clust)])],p=p)
  chi0 <- Gamma2chi(G0)
  chi1 <- Gamma2chi(Gam_colListList[[clust]][[kbest_list[clust]]])
  ggp <- plot_fitted_params(chi0, chi1, colors = makeColorMat(graph=graph_emtp2_list[[clust]], parti = partiListList[[clust]][[kbest_list[clust]]]), isFirst = TRUE)
  plot(ggp)
  dev.off()
}




