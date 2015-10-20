makeTOMHeatmap <- function(){
	library(WGCNA)
	if(!exists('d0.net')){
		getNet()
	}
	dissTOM = 1 - TOMsimilarityFromExpr( t(d0), power=16)
	plotTOM = dissTOM^7 # Makes the resulting plot easier to see
	# Set diagonal to NA for a nicer plot
	diag(plotTOM) = NA;
	# Call the plot function
	sizeGrWindow(9,9)
	TOMplot(plotTOM, d0.geneTree, d0.moduleColors, main="Network heatmap plot, all genes")
}

plotSoftThreshPowers <- function(ngenes=1000){
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

getNet <- function(p=16,ngenes=1000){
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

