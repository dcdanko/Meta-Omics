
plotTOMHeatmap <- function(day='peak',ngenes=1000){
	library(WGCNA)
	if(!exists('dnet')){
		getNet(day=day, ngenes=ngenes)
	}
	data <- days[[day]]
	dissTOM <<- 1 - TOMsimilarityFromExpr( t(data), power=16)
	plotTOM = dissTOM^7 # Makes the resulting plot easier to see
	# Set diagonal to NA for a nicer plot
	diag(plotTOM) = NA;
	# Call the plot function
	sizeGrWindow(9,9)
	TOMplot(plotTOM, dnet$day$dendrograms[[1]], labels2colors(dnet$day$colors), main="Network heatmap plot, all genes")
}

plotSoftThreshPowers <- function(day='peak', ngenes=1000){
	library(WGCNA)
	if(!exists('days')){
		source('load-data.r')
		getEveryStageVariableGenes(ngenes)
	}
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to=40, by=2))
	# Call the network topology analysis function
	data <- days[[day]]
	sft = pickSoftThreshold(t(data), powerVector = powers, verbose = 5)
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

getNet <- function(day='peak',p=16,ngenes=1000){
	if(!exists('days')){
		source('load-data.r')
		getEveryStageVariableGenes(ngenes=ngenes)
	}
	library(WGCNA)
	dnet <<- list()
	dnet$day <<- blockwiseModules(t(days[[day]]), power = p,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)

	# d0.moduleLabels <<- dnet$day$colors
	# d0.moduleColors <<- labels2colors(dnet$day$colors)
	# d0.MEs <<- dnet$day$MEs;
	# d0.geneTree <<- dnet$day$dendrograms[[1]];
}


plotTOMDendo <- function(day='peak',ngenes=1000){
	library(WGCNA)
	if(!exists('dnet')){
		getNet(day=day, ngenes=ngenes)
	}
	dissTOM <<- 1 - TOMsimilarityFromExpr( t(days[[day]]), power=16)
	# Calculate the dendrogram
	hierTOM <<- hclust(as.dist(dissTOM),method="average");
	staticHeight=0.95
	staticcut <<- cutreeStaticColor(hierTOM, cutHeight=staticHeight, minSize=20)
	colorStaticTOM <- as.character(staticcut)
	dynamiccut <<- cutreeDynamic(hierTOM,method="tree")
	colorDynamicTOM <- labels2colors (dynamiccut)
	colorDynamicHybridTOM <- labels2colors(cutreeDynamic(hierTOM, distM=dissTOM , cutHeight=staticHeight,
	                       deepSplit=2, pamRespectsDendro = FALSE))
	# Now we plot the results
	sizeGrWindow(10,5)
	plotDendroAndColors(hierTOM, 
       	colors=data.frame(colorStaticTOM,colorDynamicTOM,colorDynamicHybridTOM), 
       	dendroLabels = FALSE,
       	main = "Gene dendrogram and module colors, TOM dissimilarity",
       	abHeight=staticHeight)
	selectedTOMColoring <<- colorStaticTOM
}