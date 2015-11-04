
getNormNOGs <- function(){
	if(!exists('normnogs')){
		normnogs <<- getClipNormNOGs(8000)
	}
	return(normnogs)
}

getClipNormNOGs <- function(ngenes){
	source('microarray.r')
	print("Reading in normalized NOG counts")
	normnogs <- readInMicroarray('hh_ail10r_timecourse_metagenomics/data/normalised_counts/gene_counts.norm.matrix')

	splits <- strsplit(colnames(normnogs), split='[.]')
	rats <- sapply(sapply(splits,head,n=3), tail,n=1)
	urats <- unique(rats)
	r <- list(mode="matrix")
	for (arat in urats){
		r[[arat]] <- normnogs[,grepl(arat,colnames(normnogs))]
	}

	print("Reducing to the 8,000 most varying NOG counts")
	rA <- r$R1+r$R2+r$R3+r$R4+r$R5+r$R6+r$R7+r$R8
	rA <- rA / 8
	rA.d0 <- rA[,1,drop=FALSE]
	baseline <- cbind(rA.d0,rA.d0,rA.d0,rA.d0,rA.d0)
	rA.logfold <- log(rA / baseline)

	rowf <- function(row){
		return( sum(abs(row)))
	}

	deltavec <- apply(rA.logfold,1,rowf)
	delta <- matrix(deltavec,ncol=1)
	rownames(delta) <- names(deltavec)
	delta <- delta[ order(delta[,1]),,drop=FALSE]
	delta <- delta[ dim(delta)[1]:1,,drop=FALSE]
	# plot(1:dim(delta)[1],delta[,1])
	delta <- delta[1:ngenes,,drop=FALSE]


	normnogs <<- normnogs[rownames(delta),]

}


getNormNOGsByDay <- function(ngenes=-1){
	if(!exists('normnogs.byday')){
		if(ngenes == -1){
			normnogs <- getNormNOGs()
		}
		source('microarray.r')
		# names are of the form: 'stool.day6.R1.diamond_count'
		# make a list of days
		splits <- strsplit(colnames(normnogs), split='[.]')
		days <- sapply(sapply(splits,head,n=2), tail,n=1)
		udays <- unique(days)
		normnogs.byday <<- list(mode="matrix")
		for (aDay in udays){
			normnogs.byday[[aDay]] <<- normnogs[,grepl(aDay,colnames(normnogs))]
		}
	}
	return(normnogs.byday)
}

getNormNOGsByRat <- function(){
	if(!exists('normnogs.byrat')){
		normnogs <- getNormNOGs()

		# names are of the form: 'stool.day6.R1.diamond_count'
		# make a list of days
		splits <- strsplit(colnames(normnogs), split='[.]')
		rats <- sapply(sapply(splits,head,n=3), tail,n=1)
		urats <- unique(rats)
		normnogs.byrat <<- list(mode="matrix")
		for (arat in urats){
			normnogs.byrat[[arat]] <<- normnogs[,grepl(arat,colnames(normnogs))]
		}
	}
	return(normnogs.byrat)
}

getNOGCorrelationNet <- function(day='day6',p=16){
	daynog <- getNormNOGsByDay()[[day]]
	library(WGCNA)
	nognet <- blockwiseModules(t(daynog), power = p,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)
	return(nognet)
}

plotSoftThreshPowers <- function(day='day6'){
	library(WGCNA)
	data <- getNormNOGsByDay()[[day]]
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to=30, by=2))
	# Call the network topology analysis function
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


plotNOGTomHeatmap <- function(day='day6',p=16){
	daynog <- getNormNOGsByDay()[[day]]
	nognet <- getNOGCorrelationNet(day=day,p=p)
	daynog <- daynog[nognet$goodGenes == TRUE,]
	netcols <- nognet$colors[nognet$goodGenes == TRUE]

	library(WGCNA)
	dissTOM <- 1 - TOMsimilarityFromExpr(t(daynog),power=p)
	plotTOM <- dissTOM^7
	diag(plotTOM) <- NA

	fname <- paste(day,'-tom-heatmap')
	#fname <- paste('tom-heatmaps/',fname)
	png(fname)

	sizeGrWindow(9,9)
	TOMplot(plotTOM, nognet$dendrograms[[1]], labels2colors(netcols), main="NOG TOM Network")
	dev.off()
}