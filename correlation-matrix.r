# library(igraph)
# library(Hmisc)

# data <- t(fifty.d0)
# N <- NCOL(data)

# spcor <- rcorr(data,type="spearman")

# adja <- matrix(0,N,N)
# pthresh <- 0.05
# effectthresh <- 0.5
# for (i in 1:N){
# 	for(j in 1:N){
# 		if( i!= j){
# 			if( spcor$P[i,j] <= pthresh){
# 			 	sp <- spcor$r[i,j]
# 			 	if(sp < -1*effectthresh){
# 			 		adja[i,j] <- -1
# 			 	} else if(sp > effectthresh){
# 			 		adja[i,j] <- 1
# 			 	}
# 			}
# 		}
# 	}
# }

# net <- graph.adjacency(adja,mode="undirected", weighted=TRUE,diag=FALSE)
# E(net)[weight < 0]$color <- "RED"
# E(net)[weight > 0]$color <- "BLUE"

# print("Plotting graph...")
# plot(net, vertex.color="SKYBLUE")
