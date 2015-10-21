
# source("microarray.r")
# pg <- loadProbeToGeneMap("data/probe2gene.map")
# m <- readInMicroarray("data/sample_probe_profile.matrix")
# r <- splitMicroarrayMatrixByRat(m)
# r1 <- r$R1
# r1.gene <- convertProbeDataToGeneData(r1,pg)
# r1.gene <- r1.gene[,1:7]
# r1.base <- r1.gene[,1]
# r1.base.m <- cbind(r1.base,r1.base,r1.base)
# r1.base.m <- cbind(r1.base.m,r1.base.m,r1.base)
# r1.delta <- r1.gene - r1.base

# data <- t(r1.delta)
# data <- t(r1.delta[1:5000,])
# source("correlation-matrix.r")

# rowinterp <- function(r){
#   #                   d0   d1   d2 d3   d4 d5 d6  
#   infr <- na.spline(c(r[1],r[2],NA,r[3],NA,NA,r[4],
#   #                   d7 d8 d9 d10 d11 d12 d13 d14
#                       NA,NA,NA,NA, NA, NA, NA, r[5],
#   #                   d15 d16 d17 d18 d19 d20 d22
#                       NA, NA, NA, NA, NA, NA, r[6],
#   #                   d23 d24 d25 d26 d27 d28                    
#                       NA, NA, NA, NA, NA, r[7]))
#   return(infr)
# }

#  c <- matrix(apply(r1.delta[1:2,],1,FUN=rowinterp),byrow=TRUE,nrow=2)


# > source("clustering.r")
# > d.r1.infdel <- makeDistMatrix(r1.delta.inf
# r1.delta.inf
# > d.r1.infdel <- makeDistMatrix(r1.delta.inf, 'dtw')
# ^C
# ^C  
# ^C^C> 
# > d.r1.infdel <- makeDistMatrix(r1.delta.inf[1:100,], 'dtw')
# > d.r1.5 <- kmeans(d.r1.infdel,5)
# > clustering <- d.r1.5
# > series <-r1.
# -r1.base       -r1.base.m     -r1.delta      -r1.delta.inf  -r1.gene
# > series <- r1.delta.inf
# > source('workflows.r')
# Error: object 'cnum' not found
# > cnum <- -1
# > source('workflows.r')
# Error in colors[iclust] : invalid subscript type 'list'
# > colors
# [1] "#FF0000FF" "#CCFF00FF" "#00FF66FF" "#0066FFFF" "#CC00FFFF"
# > colors[2]
# [1] "#CCFF00FF"
# > nseries
# [1] 17330
# > clustering <- d.r1.5[1:100,]
# Error in d.r1.5[1:100, ] : incorrect number of dimensions
# > clustering <- d.r1.5
# > series <- r1.delta.inf[1:100,]
# > d.r1.5$cluster[3]
# 0610005H09RIK 
#             2 
# > clustering <- d.r1.5$cluster
# > source('workflows.r')
# > 


# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 10, to=100, by=5))
# # Call the network topology analysis function
# sft = pickSoftThreshold(d0, powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# net = blockwiseModules(d0, power = 16,
# TOMType = "unsigned", minModuleSize = 30,
# reassignThreshold = 0, mergeCutHeight = 0.25,
# numericLabels = TRUE, pamRespectsDendro = FALSE,
# saveTOMs = TRUE,
# saveTOMFileBase = "mouse-toms",
# verbose = 3)

# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
# "Module colors",
# dendroLabels = FALSE, hang = 0.03,
# addGuide = TRUE, guideHang = 0.05)

# MEs = moduleEigengenes(d0, moduleColors)$eigengenes 
# MET = orderMEs(MEs)
# sizeGrWindow(5,7.5);
# par(cex = 0.9)
# plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
# = 90)

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory.  On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Here are input parameters of the simulation model
# number of samples or microarrays in the training data
no.obs=50
# now we specify the true measures of eigengene significance
# recall that ESturquoise=cor(y,MEturquoise)
ESturquoise=0;   ESbrown= -.6;
ESgreen=.6;ESyellow=0
# Note that we don’t specify the eigengene significance of the blue module
# since it is highly correlated with the turquoise module.
ESvector=c(ESturquoise,ESbrown,ESgreen,ESyellow)
# number of genes 
nGenes1=3000
# proportion of genes in the turquoise, blue, brown, green, and yellow module #respectively.
simulateProportions1=c(0.2,0.15, 0.08, 0.06, 0.04)
# Note that the proportions don’t add up to 1. The remaining genes will be colored grey,
# ie the grey genes are non-module genes.
# set the seed of the random number generator. As a homework exercise change this seed.
set.seed(1)
#Step 1: simulate a module eigengene network.
# Training Data Set I
MEgreen=rnorm(no.obs)
scaledy=MEgreen*ESgreen+sqrt(1-ESgreen^2)*rnorm(no.obs)
y=ifelse( scaledy>median(scaledy),2,1)
MEturquoise= ESturquoise*scaledy+sqrt(1-ESturquoise^2)*rnorm(no.obs)
# we simulate a strong dependence between MEblue and MEturquoise
MEblue= .6*MEturquoise+ sqrt(1-.6^2) *rnorm(no.obs)
MEbrown= ESbrown*scaledy+sqrt(1-ESbrown^2)*rnorm(no.obs)
MEyellow= ESyellow*scaledy+sqrt(1-ESyellow^2)*rnorm(no.obs)
ModuleEigengeneNetwork1=data.frame(y,MEturquoise,MEblue,MEbrown,MEgreen, MEyellow)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


dat1=simulateDatExpr5Modules(MEturquoise=ModuleEigengeneNetwork1$MEturquoise,
          MEblue=ModuleEigengeneNetwork1$MEblue,
          MEbrown=ModuleEigengeneNetwork1$MEbrown,
          MEyellow=ModuleEigengeneNetwork1$MEyellow,
          MEgreen=ModuleEigengeneNetwork1$MEgreen, 
          nGenes=nGenes1, 
          simulateProportions=simulateProportions1)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


datExpr = dat1$datExpr;
truemodule = dat1$truemodule;
datME = dat1$datME;
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


table(truemodule)
dim(datExpr)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


datExpr=data.frame(datExpr)
ArrayName=paste("Sample",1:dim(datExpr)[[1]], sep="" )   
# The following code is useful for outputting the simulated data 
GeneName=paste("Gene",1:dim(datExpr)[[2]], sep="" )   
dimnames(datExpr)[[1]]=ArrayName
dimnames(datExpr)[[2]]=GeneName


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


rm(dat1); collectGarbage();
# The following command will save all variables defined in the current session.
save.image("Simulated-dataSimulation.RData");