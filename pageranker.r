
# source("microarray.r")
# m <- readInMicroarray("data/sample_probe_profile.matrix")
# pg <- loadProbeToGeneMap("data/probe2gene.map")
# print("Finished loading files.")
# r <- splitMicroarrayMatrixByRat(m)
# rA <- r$R1+r$R2+r$R3+r$R4
# rA <- rA / 4
# delta    <- cbind(abs(rA[1]-rA[5]),sign(rA[1]-rA[5]))
# delta    <- delta[ order(delta[,1]),,drop=FALSE]
# delta    <- delta[ dim(delta)[1]:1,,drop=FALSE]
# k    <- delta[1:1000,,drop=FALSE]

# k.gene    <- sapply(rownames(k),probeToGene,pg)

# print("Finished making deltas.")

# r5 <- m[1:1000,7,drop=FALSE]
# r1 <- r$R1[,1,drop=FALSE]
# r2 <- r$R2[1:1000,1,drop=FALSE]
# r3 <- r$R3[1:1000,1,drop=FALSE]
# r4 <- r$R4[1:1000,1,drop=FALSE]

# r1.gene <- convertProbeDataToGeneData(r1,pg)
# print("Finished converting 1/5 probe data.")
# r2.gene <- convertProbeDataToGeneData(r2,pg)
# print("Finished converting 2/5 probe data.")
# r3.gene <- convertProbeDataToGeneData(r3,pg)
# print("Finished converting 3/5 probe data.")
# r4.gene <- convertProbeDataToGeneData(r4,pg)
# print("Finished converting 4/5 probe data.")
# r5.gene <- convertProbeDataToGeneData(r5,pg)
# print("Finished converting all probe data.")

# k1 <- r1.gene[names(k.gene),1,drop=FALSE]
# k2 <- r2.gene[rownames(k.gene),1,drop=FALSE]
# k3 <- r3.gene[rownames(k.gene),1,drop=FALSE]
# k4 <- r4.gene[rownames(k.gene),1,drop=FALSE]
# k5 <- r5.gene[rownames(k.gene),1,drop=FALSE]

# k.d0.gene <- cbind(k1,k2,k3,k4,k5)
# data <- t(k.d0.gene)
# source("correlation-matrix.r")
# anet <- graph.adjacency(abs(adja),mode="undirected",weighted=TRUE,diag=FALSE)
# catch <- page.rank(anet,vids=V(anet))
# top <- catch$vector[order(catch$vector)][850:1]
# print("Top 50 implicated:")
# top[1:50]

# ns <- names(top[1:20])
# vs <- V(net)[name %in% ns]
# sub <- induced.subgraph(net,vs)
# E(sub)[weight < 0]$color <- "RED"
# E(sub)[weight > 0]$color <- "BLUE"

# # for (i in V(sub)){
# #   g <- V(sub)[i]$name
# #   if(g %in% rownames(dirk.gene)){
# #     if( dirk.gene[g,1] < 0){
# #       V(sub)[i]$color <- "GREEN"
# #     } else {
# #       V(sub)[i]$color <- "SKYBLUE"
# #     }
# #   } else {
# #     V(sub)[i]$color <- "RED"
# #   }
# # }

# print("Plotting graph...")
# plot(sub)