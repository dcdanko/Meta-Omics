
readInMicroarray <- function(filename,idsfirst=FALSE){

	microarrayFile <- read.table(filename, header = TRUE)

	if(idsfirst){
		M <- microarrayFile[,-1]
		ids <- microarrayFile[,1]
	} else {
		M <- microarrayFile[,-dim(microarrayFile)[2]]
		ids <- microarrayFile[,dim(microarrayFile)[2]]
	}

	M <- sapply(M,as.numeric)
	rownames(M) <- ids
	return(M)
}

getNorms <- function(fname){
	if(!exists('norms')){
		norms <<- getClipNorms(fname, 12000) # 300 for genus, 4,000 for nogs
	}
	return(norms)
}

getClipNorms <- function(fname, ngenes){
	norms <- readInMicroarray(fname, idsfirst=FALSE)

	splits <- strsplit(colnames(norms), split='[.]')
	rats <- sapply(sapply(splits,head,n=3), tail,n=1)
	urats <- unique(rats)
	r <- list(mode="matrix")
	for (arat in urats){
		r[[arat]] <- norms[,grepl(arat,colnames(norms))]
	}

	infostr <- paste("Reducing to the ", ngenes)
	infostr <- paste(infostr," most varying counts that are not NAN")
	print(infostr)

	rA <- r$R1+r$R2+r$R3+r$R4+r$R5+r$R6+r$R7+r$R8
	rA <- rA / 8
	rA.d0 <- rA[,1,drop=FALSE]
	baseline <- cbind(rA.d0,rA.d0,rA.d0,rA.d0,rA.d0)
	rA.logfold <- log(rA / baseline)

	rowf <- function(row){
		return( sum(abs(row)))
	}

	deltavec <- apply(rA.logfold,1,rowf)
	deltavec <- deltavec[!is.na(deltavec)]
	delta <- matrix(deltavec,ncol=1)
	rownames(delta) <- names(deltavec)
	delta <- delta[ order(delta[,1]),,drop=FALSE]
	delta <- delta[ nrow(delta):1,,drop=FALSE]

	plot(1:nrow(delta),delta[,1])
	delta <- delta[1:ngenes,,drop=FALSE]


	return(norms[rownames(delta),])

}


getNormsByDay <- function(fname){
	if(!exists('norms.byday')){
		norms <- getNorms(fname)
		# names are of the form: 'stool.day6.R1.diamond_count'
		# make a list of days
		splits <- strsplit(colnames(norms), split='[.]')
		days <- sapply(sapply(splits,head,n=2), tail,n=1)
		udays <- unique(days)
		norms.byday <<- list(mode="matrix")
		for (aDay in udays){
			norms.byday[[aDay]] <<- norms[,grepl(aDay,colnames(norms))]
		}
	}
	return(norms.byday)
}

getNormsByRat <- function(fname){
	if(!exists('norms.byrat')){
		norms <- getNorms(fname)

		# names are of the form: 'stool.day6.R1.diamond_count'
		# make a list of days
		splits <- strsplit(colnames(norms), split='[.|-]')
		rats <- sapply(sapply(splits,head,n=3), tail,n=1)
		urats <- unique(rats)
		norms.byrat <<- list(mode="matrix")
		for (arat in urats){
			norms.byrat[[arat]] <<- norms[,grepl(arat,colnames(norms))]
		}
	}
	return(norms.byrat)
}

getCorrelationNet <- function(fname,day='day6',p=16){
	day <- getNormsByDay(fname)[[day]]
	library(WGCNA)
	net <- blockwiseModules(t(day), power = p,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)
	return(net)
}

plotSoftThreshPowers <- function(fname,day='day6'){
	library(WGCNA)
	data <- getNormsByDay(fname)[[day]]
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


plotTomHeatmap <- function(fname, day='day6',p=16){
	nday <- getNormsByDay(fname)[[day]]
	net <- getCorrelationNet(day=day,p=p)
	nday <- nday[net$goodGenes == TRUE,]
	netcols <- net$colors[net$goodGenes == TRUE]

	library(WGCNA)
	dissTOM <- 1 - TOMsimilarityFromExpr(t(nday),power=p)
	plotTOM <- dissTOM^7
	diag(plotTOM) <- NA

	sizeGrWindow(9,9)
	TOMplot(plotTOM, net$dendrograms[[1]], labels2colors(netcols), main="Normalized count TOM Network")

}


analyze <- function(){
	# print( paste('plotting normalized nogs on',d))
	# plotTomHeatmap('hh_ail10r_timecourse_metagenomics/data/normalised_counts/gene_counts.norm.matrix',day=d)
	# print( paste('plotting normalized genus on',d))
	

	# print( paste('plotting raw nogs on',d))
	# plotTomHeatmap('gene_counts.tsv',day=d)

	# dontprint <- getNormsByRat('hh_ail10r_timecourse_metagenomics/data/normalised_counts/gene_counts.norm.matrix')
	r1 <- norms.byrat$R1
	delr1 <- r1[,2:5] - r1[,1:4]
	baser1 <- r1[,1:4]
	r2 <- norms.byrat$R2
	delr2 <- r2[,2:5] - r2[,1:4]
	baser2 <- r2[,1:4]
	r3 <- norms.byrat$R3
	delr3 <- r3[,2:5] - r3[,1:4]
	baser3 <- r3[,1:4]
	r4 <- norms.byrat$R4
	delr4 <- r4[,2:5] - r4[,1:4]
	baser4 <- r4[,1:4]
	r5 <- norms.byrat$R5
	delr5 <- r5[,2:5] - r5[,1:4]
	baser5 <- r5[,1:4]
	r6 <- norms.byrat$R6
	delr6 <- r6[,2:5] - r6[,1:4]
	baser6 <- r6[,1:4]
	r7 <- norms.byrat$R7
	delr7 <- r7[,2:5] - r7[,1:4]
	baser7 <- r7[,1:4]
	r8 <- norms.byrat$R8
	delr8 <- r8[,2:5] - r8[,1:4]
	baser8 <- r8[,1:4]

	delall <- cbind(delr1,delr2,delr3,delr4,delr5,delr6,delr7,delr8)
	baseall <- cbind(baser1,baser2,baser3,baser4,baser5,baser6,baser7,baser8)

	colnames(delall) <- c('r1.del.d0d3','r1.del.d3d6','r1.del.d6d14','r1.del.d14d28',
		'r2.del.d0d3','r2.del.d3d6','r2.del.d6d14','r2.del.d14d28',
		'r3.del.d0d3','r3.del.d3d6','r3.del.d6d14','r3.del.d14d28',
		'r4.del.d0d3','r4.del.d3d6','r4.del.d6d14','r4.del.d14d28',
		'r5.del.d0d3','r5.del.d3d6','r5.del.d6d14','r5.del.d14d28',
		'r6.del.d0d3','r6.del.d3d6','r6.del.d6d14','r6.del.d14d28',
		'r7.del.d0d3','r7.del.d3d6','r7.del.d6d14','r7.del.d14d28',
		'r8.del.d0d3','r8.del.d3d6','r8.del.d6d14','r8.del.d14d28')

	colnames(baseall) <- c('r1.base.d0','r1.base.d3','r1.base.d6','r1.base.d14',
		'r2.base.d0','r2.base.d3','r2.base.d6','r2.base.d14',
		'r3.base.d0','r3.base.d3','r3.base.d6','r3.base.d14',
		'r4.base.d0','r4.base.d3','r4.base.d6','r4.base.d14',
		'r5.base.d0','r5.base.d3','r5.base.d6','r5.base.d14',
		'r6.base.d0','r6.base.d3','r6.base.d6','r6.base.d14',
		'r7.base.d0','r7.base.d3','r7.base.d6','r7.base.d14',
		'r8.base.d0','r8.base.d3','r8.base.d6','r8.base.d14')

	intersects <- c()
	slopes <- c()
	sps <- c()
	ips <- c()
	inds <- c()
	deps <- c()


	for (i in 1:nrow(r1)){
		for (j in 1:nrow(r1)){
			ind <- baseall[j,]
			dep <- delall[i,]
			if( length(ind[ind==0]) < 9 && length(dep[dep==0]) < 9){ # limit regression to data with enough non-zero values
				fit <- lm(ind ~ dep)
				slopepval <- summary(fit)$coefficients[2,4] # easier way to get this?
				interceptpval <- summary(fit)$coefficients[1,4] # easier way to get this?
				if(!is.na(slopepval) && slopepval < 0.05){ 
					intersects <- append(intersects, summary(fit)$coefficients[1,1])
					slopes <- append(slopes, summary(fit)$coefficients[2,1])
					sps <- append(sps, slopepval)
					ips <- append(ips, interceptpval)
					inds <- append(inds, rownames(r1)[j])
					deps <- append(deps, rownames(r1)[i])
				}
			}
		}
	}

	timedependencies <<- data.frame(inds, deps, slopes, intersects, sps, ips)

	# days, dis and norms are already in system
	# library(Hmisc)

	# disease <- c( rep(0,8), rep(1,8), rep(2,8), rep(4,8), rep(3,8))
	# nogs <- c()
	# rhos <- c()
	# ps <- c()
	# for (i in 1:nrow(norms)){
	# 	val <- norms[i,]
	# 	c <- rcorr(disease, val, type='spearman')
	# 	nog <- rownames(norms)[i]
	# 	rho <- c$r[1,2]
	# 	p <- c$P[1,2]
	# 	if(p < 0.05){
	# 		nogs <- append(nogs, nog)
	# 		rhos <- append(rhos, rho)
	# 		ps <- append(ps, p)
	# 	} 
	# }

	# spearmancorrs <<- data.frame(nogs,rhos,ps)
}