d0.plotTOMHeatmap <- function(ngenes=1000){
	library(WGCNA)
	if(!exists('d0.net')){
		d0.getNet(ngenes=1000)
	}
	if(!exists('d0.dissTOM')){
		d0.dissTOM <<- 1 - TOMsimilarityFromExpr( t(d0), power=16)
	}
	plotTOM = d0.dissTOM^7 # Makes the resulting plot easier to see
	# Set diagonal to NA for a nicer plot
	diag(plotTOM) = NA;
	# Call the plot function
	sizeGrWindow(9,9)
	TOMplot(plotTOM, d0.geneTree, d0.moduleColors, main="Network heatmap plot, all genes")
}

d0.plotSoftThreshPowers <- function(ngenes=1000){
	if(!exists('d0')){
		source('load-data.r')
		getDayZeroVariableGenes(ngenes)
	}
	library(WGCNA)
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to=40, by=2))
	# Call the network topology analysis function
	sft = pickSoftThreshold(t(d0), powerVector = powers, verbose = 5)
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

d0.getNet <- function(p=16,ngenes=1000){
	if(!exists('d0')){
		source('load-data.r')
		getDayZeroVariableGenes(ngenes)
	}
	library(WGCNA)
	d0.net <<- blockwiseModules(t(d0), power = p,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)

	d0.moduleLabels <<- d0.net$colors
	d0.moduleColors <<- labels2colors(d0.net$colors)
	d0.MEs <<- d0.net$MEs;
	d0.geneTree <<- d0.net$dendrograms[[1]];
	return(d0.net)
}


d0.plotDendo <- function(){
	if(!exists('d0.geneTree')){
		d0.getNet()
	}
	staticHeight=0.99
	staticColoring <- as.character( cutreeStaticColor(d0.geneTree, cutHeight=staticHeight, minSize=20))
	dyanmicColoring <- labels2colors( cutreeDynamic(d0.geneTree, method='tree'))
	sizeGrWindow(10,5)
	plotDendroAndColors(d0.geneTree, colors=data.frame(d0.moduleColors,staticColoring,dyanmicColoring), 
						abHeight=staticHeight,dendroLabels=FALSE)
}

d0.plotTOMDendo <- function(){
	library(WGCNA)
	if(!exists('d0.net')){
		d0.getNet(ngenes=1000)
	}
	if(!exists('d0.dissTOM')){
		d0.dissTOM <<- 1 - TOMsimilarityFromExpr( t(d0), power=16)
	}
	# Calculate the dendrogram
	d0.hierTOM <<- hclust(as.dist(d0.dissTOM),method="average");
	staticHeight=0.95
	colorStaticTOM <- as.character(cutreeStaticColor(d0.hierTOM, cutHeight=staticHeight, minSize=20))
	colorDynamicTOM <- labels2colors (cutreeDynamic(d0.hierTOM,method="tree"))
	colorDynamicHybridTOM <- labels2colors(cutreeDynamic(d0.hierTOM, distM=d0.dissTOM , cutHeight=staticHeight,
	                       deepSplit=2, pamRespectsDendro = FALSE))
	# Now we plot the results
	sizeGrWindow(10,5)
	plotDendroAndColors(d0.hierTOM, 
       	colors=data.frame(colorStaticTOM,colorDynamicTOM,colorDynamicHybridTOM), 
       	dendroLabels = FALSE,
       	main = "Gene dendrogram and module colors, TOM dissimilarity",
       	abHeight=staticHeight)
	selectedTOMColoring <<- colorStaticTOM
}

d0.tableIntramodularConnectivity <- function(){
	if(!exists('selectedTOMColoring')){
		d0.plotTOMDendo()
	}
	library(WGCNA)

	ADJ1 <- abs( cor(t(d0),use='p'))^6
	Alldegrees1 <- intramodularConnectivity(ADJ1, d0.moduleColors)

	# Makes a table of graphs comparing clusters to eigengenes.
	# Appears that most clusters (if the clustering is valid)
	# have a corresponding eigengene


	d0.KME <<- signedKME(t(d0), d0.MEs, outputColumnName='MM.')
	allColors <- names(table(d0.moduleColors)) # 4
	allEigengenes <- colnames(d0.KME) # 9 

	sizeGrWindow(16,12)
	par(mar=rep(2,4), mfrow=c(length(allColors),length(allEigengenes)))

	for (color in allColors){
		which.color <- color 
		restrictGenes <- d0.moduleColors==which.color 
		intramodularConnectivityForGenesInCluster <- Alldegrees1$kWithin[ restrictGenes]
		for (eigengeneMembership in allEigengenes)
		verboseScatterplot( intramodularConnectivityForGenesInCluster, 
			(d0.KME[restrictGenes, eigengeneMembership])^6,
			col=which.color)
			# xlab="Intramodular Connectivity", 
			# ylab="(Module Membership)^6")
	}

}
