rat1.setDefaultNProbes <- function(){
	nprobes <<- 5000
}

rat1.plotTOMHeatmap <- function(){
	library(WGCNA)
	if(!exists('rat1.net')){
		rat1.getNet()
	}
	if(!exists('rat1.dissTOM')){
		rat1.dissTOM <<- 1 - TOMsimilarityFromExpr( t(rat1), power=16)
	}
	plotTOM = rat1.dissTOM^7 # Makes the resulting plot easier to see
	# Set diagonal to NA for a nicer plot
	diag(plotTOM) = NA;
	# Call the plot function
	sizeGrWindow(9,9)
	TOMplot(plotTOM, rat1.geneTree, rat1.moduleColors, main="Network heatmap plot, all genes")
}

rat1.plotSoftThreshPowers <- function(){
	if(!exists('rat1')){
		if(!exists('nprobes')){
			rat1.setDefaultNProbes()
		}
		source('load-data.r')
		getRatOneVariableGenes(nprobes)
	}
	library(WGCNA)
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to=40, by=2))
	# Call the network topology analysis function
	sft = pickSoftThreshold(t(rat1), powerVector = powers, verbose = 5)
	# Plot the results:
	sizeGrWindow(9, 5)
	par(mfrow = c(1,2));
	cex1 = 0.9;
	# Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	     main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     labels=powers,cex=cex1,col="red");
	# this line corresponds to using an R^2 cut-off of h
	abline(h=0.80,col="red")
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
	     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	     main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

rat1.getNet <- function(p=16){
	if(!exists('rat1')){
		source('load-data.r')
		if(!exists('nprobes')){
			rat1.setDefaultNProbes()
		}
		getRatOneVariableGenes(nprobes)
	}
	library(WGCNA)
	rat1.net <<- blockwiseModules(t(rat1), power = p,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)

	rat1.moduleLabels <<- rat1.net$colors
	rat1.moduleColors <<- labels2colors(rat1.net$colors)
	rat1.MEs <<- rat1.net$MEs;
	rat1.geneTree <<- rat1.net$dendrograms[[1]];
	return(rat1.net)
}


rat1.plotDendo <- function(){
	if(!exists('rat1.geneTree')){
		rat1.getNet()
	}

	staticHeight=0.99
#	staticColoring <- as.character( cutreeStaticColor(rat1.geneTree, cutHeight=staticHeight, minSize=20))
#	dynamicColoring <- labels2colors( cutreeDynamic(rat1.geneTree, method='tree'))
	sizeGrWindow(10,5)
	plotDendroAndColors(rat1.geneTree, colors=rat1.moduleColors, #,staticColoring,dynamicColoring), 
						abHeight=staticHeight,dendroLabels=FALSE)
}

rat1.plotTOMDendo <- function(){
	library(WGCNA)
	if(!exists('rat1.net')){
		rat1.getNet()
	}
	if(!exists('rat1.dissTOM')){
		rat1.dissTOM <<- 1 - TOMsimilarityFromExpr( t(rat1), power=16)
	}
	# Calculate the dendrogram
	rat1.hierTOM <<- hclust(as.dist(rat1.dissTOM),method="average");
	staticHeight=0.95
	colorStaticTOM <- as.character(cutreeStaticColor(rat1.hierTOM, cutHeight=staticHeight, minSize=20))
	colorDynamicTOM <- labels2colors (cutreeDynamic(rat1.hierTOM,method="tree"))
	colorDynamicHybridTOM <- labels2colors(cutreeDynamic(rat1.hierTOM, distM=rat1.dissTOM , cutHeight=staticHeight,
	                       deepSplit=2, pamRespectsDendro = FALSE))
	# Now we plot the results
	sizeGrWindow(10,5)
	plotDendroAndColors(rat1.hierTOM, 
       	colors=data.frame(colorStaticTOM,colorDynamicTOM,colorDynamicHybridTOM), 
       	dendroLabels = FALSE,
       	main = "Gene dendrogram and module colors, TOM dissimilarity",
       	abHeight=staticHeight)
	selectedTOMColoring <<- colorStaticTOM
}

rat1.tableIntramodularConnectivity <- function(){
	if(!exists('selectedTOMColoring')){
		rat1.plotTOMDendo()
	}
	library(WGCNA)

	ADJ1 <- abs( cor(t(rat1),use='p'))^6
	Alldegrees1 <- intramodularConnectivity(ADJ1, rat1.moduleColors)

	# Makes a table of graphs comparing clusters to eigengenes.
	# Appears that most clusters (if the clustering is valid)
	# have a corresponding eigengene


	rat1.KME <<- signedKME(t(rat1), rat1.MEs, outputColumnName='MM.')
	allColors <- names(table(rat1.moduleColors)) # 4
	allEigengenes <- colnames(rat1.KME) # 9 

	sizeGrWindow(16,12)
	par(mar=rep(2,4), mfrow=c(length(allColors),length(allEigengenes)))

	for (color in allColors){
		which.color <- color 
		restrictGenes <- rat1.moduleColors==which.color 
		intramodularConnectivityForGenesInCluster <- Alldegrees1$kWithin[ restrictGenes]
		for (eigengeneMembership in allEigengenes)
		verboseScatterplot( intramodularConnectivityForGenesInCluster, 
			(rat1.KME[restrictGenes, eigengeneMembership])^6,
			col=which.color)
			# xlab="Intramodular Connectivity", 
			# ylab="(Module Membership)^6")
	}

}