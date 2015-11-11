reload <- function(){
	print('reloading... ')
	source('hh_analysis.r')
}

readInMicroarray <- function(filename,idsfirst=FALSE){

	extension <- tail( strsplit(filename,'[.]')[[1]], 1)
	if(extension == "csv"){
		microarrayFile <- read.table(filename, sep=",", header = TRUE)
	} else {
		microarrayFile <- read.table(filename, header = TRUE)
	}
	# print(colnames(microarrayFile))
	# print(microarrayFile[1:5,])
	# microarrayFile <- as.matrix( microarrayFile)
	# print(microarrayFile[1:5,])

	# if(idsfirst){
		M <- microarrayFile[,-which(colnames(microarrayFile) == 'taxa')]
		ids <- microarrayFile[,'taxa']
	# } else {
	# 	M <- microarrayFile[,-dim(microarrayFile)[2]]
	# 	ids <- microarrayFile[,dim(microarrayFile)[2]]
	# }

	M <- sapply(M,as.numeric)
	rownames(M) <- ids
	return(M)
}

getNorms <- function(fname,sortbyday=FALSE){
	norms <- getClipNorms(fname, -1) 

	# sort norms by day then rat
	if(sortbyday){
		m <- do.call(rbind, strsplit(colnames(norms), "[.]"))[,-1]
		ratsort <- order(do.call(rbind,strsplit(m[,2],'_count'))[,1])
		norms <- norms[,ratsort]

		m <- do.call(rbind, strsplit(colnames(norms), "[.]"))[,-1]
		daysort <- order(as.numeric(do.call(rbind,strsplit(m[,1],'day'))[,2]))
		norms <- norms[,daysort]
	} else {
		# sort norms by day then rat
		m <- do.call(rbind, strsplit(colnames(norms), "[.]"))[,-1]
		daysort <- order(as.numeric(do.call(rbind,strsplit(m[,1],'day'))[,2]))
		norms <- norms[,daysort]

		m <- do.call(rbind, strsplit(colnames(norms), "[.]"))[,-1]
		ratsort <- order(do.call(rbind,strsplit(m[,2],'_count'))[,1])
		norms <- norms[,ratsort]
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

	if(ngenes > 1){
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
		norms <- norms[rownames(delta),]
	}


	return(norms)

}


getNormsByDay <- function(fname){
	norms <- getNorms(fname,sortbyday=TRUE)

	norms.byday <<- list(mode="matrix")

	norms.byday[['d0']] <<- norms[,1:8]
	norms.byday[['d3']] <<- norms[,9:16]
	norms.byday[['d6']] <<- norms[,17:24]
	norms.byday[['d14']] <<- norms[,25:32]
	norms.byday[['d28']] <<- norms[,33:40]

	return(norms.byday)
}

getNormsByRat <- function(fname){
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

}

analyze_drivercorr_rise <- function(){
	dnames <- getDataNames()
	onames <- getOutNames('hh_results/driver_corr/','.norm.d0d3d6_rise_only.driver_corr_pearsons.csv')

	for (i in 1:length(dnames)){
		dname <- dnames[i]
		oname <- onames[i]

		dontprint <- getNormsByRat(dname)
		r1 <- norms.byrat$R1
		delr1 <- r1[,2:3] - r1[,1:2]
		baser1 <- r1[,1:2]
		r2 <- norms.byrat$R2
		delr2 <- r2[,2:3] - r2[,1:2]
		baser2 <- r2[,1:2]
		r3 <- norms.byrat$R3
		delr3 <- r3[,2:3] - r3[,1:2]
		baser3 <- r3[,1:2]
		r4 <- norms.byrat$R4
		delr4 <- r4[,2:3] - r4[,1:2]
		baser4 <- r4[,1:2]
		r5 <- norms.byrat$R5
		delr5 <- r5[,2:3] - r5[,1:2]
		baser5 <- r5[,1:2]
		r6 <- norms.byrat$R6
		delr6 <- r6[,2:3] - r6[,1:2]
		baser6 <- r6[,1:2]
		r7 <- norms.byrat$R7
		delr7 <- r7[,2:3] - r7[,1:2]
		baser7 <- r7[,1:2]
		r8 <- norms.byrat$R8
		delr8 <- r8[,2:3] - r8[,1:2]
		baser8 <- r8[,1:2]

		delall <- cbind(delr1,delr2,delr3,delr4,delr5,delr6,delr7,delr8)
		baseall <- cbind(baser1,baser2,baser3,baser4,baser5,baser6,baser7,baser8)

		colnames(delall) <- c('r1.del.d0d3','r1.del.d3d6',
			'r2.del.d0d3','r2.del.d3d6',
			'r3.del.d0d3','r3.del.d3d6',
			'r4.del.d0d3','r4.del.d3d6',
			'r5.del.d0d3','r5.del.d3d6',
			'r6.del.d0d3','r6.del.d3d6',
			'r7.del.d0d3','r7.del.d3d6',
			'r8.del.d0d3','r8.del.d3d6')

		colnames(baseall) <- c('r1.base.d0','r1.base.d3',
			'r2.base.d0','r2.base.d3',
			'r3.base.d0','r3.base.d3',
			'r4.base.d0','r4.base.d3',
			'r5.base.d0','r5.base.d3',
			'r6.base.d0','r6.base.d3',
			'r7.base.d0','r7.base.d3',
			'r8.base.d0','r8.base.d3')

		analyze_drivercorr(dname,oname,delall,baseall)
	}


}

analyze_drivercorr_all <- function(){
	dnames <- getDataNames()
	onames <- getOutNames('hh_results/driver_corr/','.norm.d0d3d6_rise_only.driver_corr.csv')

	for (i in 1:length(dnames)){
		dname <- dnames[i]
		oname <- onames[i]

		dontprint <- getNormsByRat(dname)
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

		analyze_drivercorr(dname,oname,delall,baseall)
	}
}

analyze_drivercorr <- function(dname, oname, delall, baseall){
	library(Hmisc)

	intersects <- c()
	slopes <- c()
	slopepvals <- c()
	interceptpvals <- c()
	inds <- c()
	deps <- c()
	pearsonsrs <- c()
	pearsonsrpvals <- c()

	baseall <- cbind(rownames(baseall),baseall)
	delall <- cbind(rownames(delall),delall)

	depfunc <- function(dep){
		apply(baseall,1,basefunc,dep)
		
	}

	basefunc <- function(indy,depy){
		indname <- indy[1]
		ind <- indy[-1]
		depname <- depy[1]
		dep <- depy[-1]

		zerolim <- length(ind) / 4

		if( length(ind[ind==0]) <= zerolim && length(dep[dep==0]) <= zerolim){ # limit regression to data with enough non-zero values

			c <- rcorr(dep, ind, type='pearson')
			pr <- c$r[1,2]
			prpval <- c$P[1,2]
			if(prpval < 0.05){

				fit <- lm(dep ~ ind)
				slopepval <- summary(fit)$coefficients[2,4] # easier way to get this?
				interceptpval <- summary(fit)$coefficients[1,4] # easier way to get this?
				# if(!is.na(slopepval) && slopepval < 0.05){ 
					intersects <<- append(intersects, summary(fit)$coefficients[1,1])
					slopes <<- append(slopes, summary(fit)$coefficients[2,1])
					slopepvals <<- append(slopepvals, slopepval)
					interceptpvals <<- append(interceptpvals, interceptpval)
					inds <<- append(inds, indname)
					deps <<- append(deps, depname)
					pearsonsrs <<- append(pearsonsrs, pr)
					pearsonsrpvals <<- append(pearsonsrpvals, prpval)
			}
		}
	}

	apply(delall,1,depfunc)

	timedependencies <- data.frame(inds, deps, slopes, intersects, slopepvals, interceptpvals, pearsonsrs, pearsonsrpvals)
	write.csv(timedependencies,oname)
}

analyze_daycorr <- function(){
	# find taxa with high sp rho to day and regress out the effects of day

	dnames <- getDataNames()
	dnames <- dnames
	onames <- getOutNames('hh_results/diseaseday_corr/','.norm.diseaseday_corr.csv')
	onames <- onames
	adjnames <- getOutNames('hh_results/diseaseday_corr/','.norm.diseaseday_adjusted.matrix.csv')

	for (i in 1:length(dnames)){
		dname <- dnames[i]
		oname <- onames[i]
		adjname <- adjnames[i]

		norms <- getNorms(dname)
		disease <- rep(c(0,1,2,4,3),8)

		# delall <- cbind(delr1,delr2,delr3,delr4,delr5,delr6,delr7,delr8)
		# baseall <- cbind(baser1,baser2,baser3,baser4,baser5,baser6,baser7,baser8)


		intersects <- c()
		slopes <- c()
		sps <- c()
		ips <- c()
		names <- c()
		adjnorms <- matrix(ncol=40,nrow=0)

		basefunc <- function(normy){
			name <- normy[1]
			data <- as.numeric(normy[-1])
			adjusted <- FALSE

			if( length(data[data==0]) <=8 ){ # limit regression to data with enough non-zero values
				fit <- lm(data ~ disease)
				slopepval <- summary(fit)$coefficients[2,4] # easier way to get this?
				interceptpval <- summary(fit)$coefficients[1,4] # easier way to get this?
				if(!is.na(slopepval) && slopepval < 0.05){ 
					intersects <<- append(intersects, summary(fit)$coefficients[1,1])
					slopes <<- append(slopes, summary(fit)$coefficients[2,1])
					sps <<- append(sps, slopepval)
					ips <<- append(ips, interceptpval)
					names <<- append(names, name)
					adjusted <- TRUE

					predictedvals <- predict(fit, data.frame(disease=c(0,1,2,4,3)))
					adjdata <- data - predictedvals

					oldnames <- rownames(adjnorms)
					adjnorms <<- rbind( adjnorms, adjdata)
					rownames(adjnorms) <<- c(oldnames,name)

				}
			}

			if( !adjusted){
				oldnames <- rownames(adjnorms)
				adjnorms <<- rbind( adjnorms, data)
				rownames(adjnorms) <<- c(oldnames,name)
			}
		}
		normswithnames <- cbind(rownames(norms),norms)
		apply(normswithnames,1,basefunc)
		colnames(adjnorms) <- colnames(norms)


		timedependencies <- data.frame(names, slopes, intersects, sps,ips)
		write.csv(timedependencies,oname)

		write.csv(adjnorms,adjname)
	}
}

analyze_spcorr <- function(){
	norms <- getNorms('hh_data_norm/class.diamond.aggregated.counts.norm.matrix')
	library(Hmisc)

	disease <- c( rep(0,8), rep(1,8), rep(2,8), rep(4,8), rep(3,8))
	nogs <- c()
	rhos <- c()
	ps <- c()
	for (i in 1:nrow(norms)){
		val <- norms[i,]
		c <- rcorr(disease, val, type='spearman')
		nog <- rownames(norms)[i]
		rho <- c$r[1,2]
		p <- c$P[1,2]
		if(p < 0.05){
			nogs <- append(nogs, nog)
			rhos <- append(rhos, rho)
			ps <- append(ps, p)
		} 
	}

	spearmancorrs <<- data.frame(nogs,rhos,ps)
	write.csv(spearmancorrs, 'hh_results/class.norms.diseasebyday.spearman_corr.csv')
}



analyze_abundances <- function(){
	dnames <- getDataNames()
	onames <- getOutNames('hh_results/abundances/','.norm.abundances.csv')
	for (i in 1:length(dnames)){
		dname <- dnames[i]
		oname <- onames[i]

		norms <- getNorms(dname)
		sc <- scale(norms, center=FALSE, scale=colSums(norms))
		write.csv(sc, oname)
	}
}

getDataNames <- function(){
	# names <- c('hh_data_norm/species.diamond.aggregated.counts.norm.matrix',
	# 	'hh_data_norm/genus.diamond.aggregated.counts.norm.matrix',
	# 	'hh_data_norm/family.diamond.aggregated.counts.norm.matrix',
	# 	'hh_data_norm/order.diamond.aggregated.counts.norm.matrix',
	# 	'hh_data_norm/class.diamond.aggregated.counts.norm.matrix',
	# 	'hh_data_norm/phylum.diamond.aggregated.counts.norm.matrix',
	# 	'hh_data_norm/gene_counts.norm.matrix')

	names <- c('hh_data_norm/genus.diamond.aggregated.counts.norm.matrix')

	# names <- c('/Users/dcdanko/Science/Meta-Omics/hh_results/diseaseday_corr/genus.norm.diseaseday_adjusted.matrix.csv')
	return(names)


}

getOutNames <- function(prefix, suffix){
	# names <- c('species','genus','family','order','class','phylum','gene_counts')
	names <- c('genus')
	onames <- paste(prefix,names,sep='')
	onames <- paste(onames,suffix,sep='')
	return(onames)
}

getBaseNames <- function(){
	base <- c('species','genus','family','order','class','phylum','gene_counts')
	return(base)
}

analyze_movers <- function(){
	dnames <- getDataNames()
	onames <- getOutNames('hh_results/early_movers/','.norm.8.early_movers.csv')
	bnames <- getBaseNames()

	for(startday in 1:4){
		for(endday in (startday+1):5){

			dnames <- getDataNames()
			suffix <- paste('.norm.8.',startday,sep='')
			suffix <- paste(suffix,'to',sep='')
			suffix <- paste(suffix,endday,sep='')
			suffix <- paste(suffix,'.movers.csv',sep='')
			onames <- getOutNames('hh_results/movers/',suffix)
			bnames <- getBaseNames()


	for (i in 1:length(dnames)){


		dname <- dnames[i]
		oname <- onames[i]
		print(oname)
		bname <- bnames[i]

		dontprint <- getNormsByRat(dname)
		r1 <- norms.byrat$R1
		br1 <- (r1[,endday,drop=FALSE] - r1[,startday,drop=FALSE])
		r2 <- norms.byrat$R2
		br2 <- (r2[,endday,drop=FALSE] - r2[,startday,drop=FALSE])
		r3 <- norms.byrat$R3
		br3 <- (r3[,endday,drop=FALSE] - r3[,startday,drop=FALSE])
		r4 <- norms.byrat$R4
		br4 <- (r4[,endday,drop=FALSE] - r4[,startday,drop=FALSE])

		r5 <- norms.byrat$R5
		br5 <- (r5[,endday,drop=FALSE] - r5[,startday,drop=FALSE])
		r6 <- norms.byrat$R6
		br6 <- (r6[,endday,drop=FALSE] - r6[,startday,drop=FALSE])
		r7 <- norms.byrat$R7
		br7 <- (r7[,endday,drop=FALSE] - r7[,startday,drop=FALSE])
		r8 <- norms.byrat$R8
		br8 <- (r8[,endday,drop=FALSE] - r8[,startday,drop=FALSE])

		sumr <- sign(br1)+sign(br2)+sign(br3)+sign(br4)+sign(br5)+sign(br6)+sign(br7)+sign(br8)
		aver <- (br1+br2+br3+br4+br5+br6+br7+br8)/8
		early_risers <- rownames(sumr[sumr >= 8,,drop=FALSE])
		early_droppers <- rownames(sumr[sumr <= -8,,drop=FALSE])
		aver_risers <- aver[early_risers,,drop=FALSE]
		aver_droppers <- aver[early_droppers,,drop=FALSE]
		aver <- rbind(aver_risers,aver_droppers)
		aver <- cbind(rownames(aver),aver[,1])
		rownames(aver) <- NULL
		title <- paste('ave_change_',startday,sep='')
		title <- paste(title,'_to_',sep='')
		title <- paste(title,endday,sep='')
		title <- paste(title,'_direction_conserved_in_8_of_8',sep='')
		colnames(aver) <- c(bname,title)
		write.csv(aver,oname)
	}
	}}

}

get_driver_matrix <- function(fname, pcut=0.001){
	drives <- read.csv(fname)[,-1]

	m <- matrix(nrow=0,ncol=0)
	matrixmaker <- function(driver){
		#print(driver)
		ind <- driver[1]
		dep <- driver[2]
		slope <- as.numeric(driver[3])
		if(length(driver)>=7){ # indicates we have a p value for pearson's r
			p <- as.numeric(driver[7])
		} else {
			p <- as.numeric(driver[5])
		}
		if(p<=pcut){
			#print(m)
			if( ind %in% colnames(m)){
				if( !(dep %in% rownames(m))) {
					newrow <- rep(NA,ncol(m))
					oldnames <- rownames(m)
					m <<- rbind(m,newrow)
					rownames(m) <<- c(oldnames,dep)
				}
				m[dep,ind] <<- slope
			} else if(dep %in% rownames(m)){
				# we already know ind is not in colnames(m)
				newcol <- rep(NA,nrow(m))
				oldnames <- colnames(m)
				m <<- cbind(m,newcol)
				colnames(m) <<- c(oldnames,ind)
				m[dep,ind] <<- slope
			} else {
				newrow <- rep(NA,ncol(m))
				oldnames <- rownames(m)
				m <<- rbind(m,newrow)
				rownames(m) <<- c(oldnames,dep)
				newcol <- rep(NA,nrow(m))
				oldnames <- colnames(m)
				m <<- cbind(m,newcol)
				colnames(m) <<- c(oldnames,ind)
				m[dep,ind] <<- slope
			}
		}
	}

	apply(drives,1,matrixmaker)
	return(m)
}



plot_driver_heatmap <- function(fname){
	library(gplots)
	
	m <- get_driver_matrix(fname)
	m <- m[order(rownames(m)),order(colnames(m))]
	heatmap.2(m,Rowv='manhattan',Colv='manhattan',dendrogram='none',col=cm.colors(256),na.color='gray',margins=c(10,10),trace='none')
}

plot_driver_heatmap_abundance_normal <- function(fname,sname){
	library(gplots)
	
	m <- get_driver_matrix(fname)
	m <- m[order(rownames(m)),order(colnames(m))]
	g <- get_ave_abundances_by_day_d0d3(sname)

	for (c in 1:ncol(m)){
		n <- g[colnames(m)[c]]
		m[,c] <- m[,c] / n
	}

	heatmap.2(m,Rowv=NA,Colv=NA,col=cm.colors(256),na.color='gray',margins=c(10,10),trace='none')
}

plot_driver_barplot_abundance_normal <- function(fname,sname){
	m <- get_driver_matrix(fname)
	g <- get_ave_abundances_by_day_d0d3(sname)

	for (c in 1:ncol(m)){
		n <- g[colnames(m)[c]]
		m[,c] <- m[,c] / n
	}

	# dependent plot
	depm <- matrix(nrow=nrow(m),ncol=5)
	rownames(depm) <- rownames(m)
	colnames(depm) <- c('Num. Pos.', 'Num. Neg.', 'Total Pos.', 'Total Neg.', 'Total')

	rowfunc <- function(r){
		name <- r[1]
		data <- as.numeric(r[-1])
		npos <- 0
		nneg <- 0
		tpos <- 0
		tneg <- 0
		for (val in data){
			if(!is.na(val)){
				if(val > 0){
					npos <- npos + 1
					tpos <- tpos + val
				} else {
					nneg <- nneg + 1
					tneg <- tneg + abs(val)
				}
			}
		}
		t <- tpos + tneg
		depm[name,] <<- c(npos, nneg, tpos, tneg, t)
	}

	apply(cbind(rownames(m),m), 1, rowfunc)

	sums <- depm[,3] + abs(depm[,4])
	depm <- depm[rev(order(sums)),]
	

	# independent plot

	indm <<- matrix(nrow=ncol(m),ncol=5)
	rownames(indm) <<- colnames(m)
	colnames(indm) <<- c('Num. Pos.', 'Num. Neg.', 'Total Pos.', 'Total Neg.', 'Total')

	rowfunc2 <- function(r){
		name <- r[1]
		data <- as.numeric(r[-1])
		npos <- 0
		nneg <- 0
		tpos <- 0
		tneg <- 0
		for (val in data){
			if(!is.na(val)){
				if(val > 0){
					npos <- npos + 1
					tpos <- tpos + val
				} else {
					nneg <- nneg + 1
					tneg <- tneg + abs(val)
				}
			}
		}
		t <- tpos + tneg
		indm[name,] <<- c(npos, nneg, tpos, tneg, t)
	}

	apply(cbind(colnames(m),t(m)), 1, rowfunc2)

	sums <- indm[,3] + indm[,4]
	sums <- indm[,5]
	indm <<- indm[rev(order(sums)),]

	# nf <- layout(matrix(c(1,2),nrow=2))
	# layout.show(nf)

	par(mfrow=c(3,2))
	barplot(t(depm[1:25,]), beside=TRUE, col=cm.colors(5)[c(1,4,2,5,3)], las=2)
	barplot(t(indm[1:25,]), beside=TRUE, legend=colnames(indm), col=cm.colors(5)[c(1,4,2,5,3)], las=2)

	barplot(t(depm[1:25,1:2]), col=cm.colors(2), las=2, axisnames=FALSE)
	barplot(t(indm[1:25,1:2]), col=cm.colors(2), las=2, axisnames=FALSE)

	barplot(t(depm[1:25,3:4]), col=cm.colors(2), las=2, axisnames=FALSE)
	barplot(t(indm[1:25,3:4]), col=cm.colors(2), las=2, axisnames=FALSE)

}

plot_driver_barplot <- function(fname){
	m <- get_driver_matrix(fname)

	# dependent plot
	depm <- matrix(nrow=nrow(m),ncol=5)
	rownames(depm) <- rownames(m)
	colnames(depm) <- c('Num. Pos.', 'Num. Neg.', 'Total Pos.', 'Total Neg.', 'Total')

	rowfunc <- function(r){
		name <- r[1]
		data <- as.numeric(r[-1])
		npos <- 0
		nneg <- 0
		tpos <- 0
		tneg <- 0
		for (val in data){
			if(!is.na(val)){
				if(val > 0){
					npos <- npos + 1
					tpos <- tpos + val
				} else {
					nneg <- nneg + 1
					tneg <- tneg + abs(val)
				}
			}
		}
		t <- tpos + tneg
		depm[name,] <<- c(npos, nneg, tpos, tneg, t)
	}

	apply(cbind(rownames(m),m), 1, rowfunc)

	sums <- depm[,3] + abs(depm[,4])
	depm <- depm[rev(order(sums)),]
	

	# independent plot

	indm <- matrix(nrow=ncol(m),ncol=5)
	rownames(indm) <- colnames(m)
	colnames(indm) <- c('Num. Pos.', 'Num. Neg.', 'Total Pos.', 'Total Neg.', 'Total')

	rowfunc2 <- function(r){
		name <- r[1]
		data <- as.numeric(r[-1])
		npos <- 0
		nneg <- 0
		tpos <- 0
		tneg <- 0
		for (val in data){
			if(!is.na(val)){
				if(val > 0){
					npos <- npos + 1
					tpos <- tpos + val
				} else {
					nneg <- nneg + 1
					tneg <- tneg + abs(val)
				}
			}
		}
		t <- tpos + tneg
		indm[name,] <<- c(npos, nneg, tpos, tneg, t)
	}

	apply(cbind(colnames(m),t(m)), 1, rowfunc2)

	sums <- indm[,3] + indm[,4]

	indm <- indm[rev(order(sums)),]

	# nf <- layout(matrix(c(1,2),nrow=2))
	# layout.show(nf)

	par(mfrow=c(3,2))
	barplot(t(depm[1:25,]), beside=TRUE, col=cm.colors(5)[c(1,4,2,5,3)], las=2)
	barplot(t(indm[1:25,]), beside=TRUE, legend=colnames(indm), col=cm.colors(5)[c(1,4,2,5,3)], las=2)

	barplot(t(depm[1:25,1:2]), col=cm.colors(2), las=2, axisnames=FALSE)
	barplot(t(indm[1:25,1:2]), col=cm.colors(2), las=2, axisnames=FALSE)

	barplot(t(depm[1:25,3:4]), col=cm.colors(2), las=2, axisnames=FALSE)
	barplot(t(indm[1:25,3:4]), col=cm.colors(2), las=2, axisnames=FALSE)



}

get_ave_abundances_by_day_d0d3 <- function(fname){
	normsbyday <- getNormsByDay(fname)
	d0 <- normsbyday$d0 
	d3 <- normsbyday$d3

	dtotal <- d0 + d3
	dave <- dtotal / 2
	rowsums <- apply(dave, 1, sum)
	rowave <- rowsums / 8

	return(rowave)
}

get_driver_adjacency <- function(fname,pcut=0.001,power=4,edgecut=0.5){
	m <- get_driver_matrix(fname,pcut=pcut)
	allnodes <- union(colnames(m),rownames(m))
	N <- length(allnodes)
	adj <- matrix(0, ncol=N, nrow=N)
	rownames(adj) <- allnodes
	colnames(adj) <- allnodes

	adjfunc <- function(x,pow){
		val <- abs(x)^pow 
		if(val > edgecut){
			return(sign(x) * val)
		}
		return(0)
	}

	for (c in 1:ncol(m)){
		ind <- colnames(m)[c]
		for (r in 1:nrow(m)){
			dep <- rownames(m)[r]
			val <- m[r,c]
			if(!is.na(val)){
				adj[ind,dep] <- adjfunc( val, power)
			}
		}
	}

	zerorows <- apply(abs(adj),1,sum)
	zerorows <- zerorows[zerorows==0]
	zerocols <- apply(abs(adj),2,sum)
	zerocols <- zerocols[zerocols==0]
	zeronames <- intersect(names(zerorows),names(zerocols))

	adj <- adj[!rownames(adj) %in% zeronames, !colnames(adj) %in% zeronames]

	return(adj)
}

plot_driver_adjacency_graph <- function(driverf){
	library(igraph)
	adj <- get_driver_adjacency(driverf)
	net <- graph.adjacency(adj,weighted=TRUE,mode='directed')
	V(net)$size <- degree(net) / 2
	V(net)$color <- 'plum'
	E(net)[weight < 0]$color <- 'lightblue'
	E(net)[weight >= 0]$color <- 'lightcoral'
	med <- median( abs(E(net)$weight))
	E(net)[weight > med]$color <- 'coral'
	E(net)[weight < (-1*med)]$color <- 'lightslategrey'
	par(mai=c(0,0,1,0)) 
	plot(net, vertex.label.cex=0.7,layout=layout.fruchterman.reingold,vertex.label.dist=0.2,vertex.label.color='black',vertex.frame.color='plum')

}

plot_driver_adjacency_graph_with_abundance <- function(driverf, normf,pcut=0.001,power=4,edgecut=0.5,sizecut=1){
	library(igraph)
	abund <- get_ave_abundances_by_day_d0d3(normf)

	adj <- get_driver_adjacency(driverf,pcut=pcut,power=power,edgecut=edgecut)
	net <- graph.adjacency(adj,weighted=TRUE,mode='directed')

	# edit the graph
	V(net)$size <- abund[V(net)$name] 
	vertstoremove <- V(net)[V(net)$size <= sizecut]
	net <- delete.vertices(net, vertstoremove)
	vertstoremove <- V(net)[which( degree(net) <= 0)]
	net <- delete.vertices(net, vertstoremove)

	# plot a pretty graph
	V(net)$size <- V(net)$size / 2
	V(net)$color <- 'plum'
	E(net)[weight < 0]$color <- 'lightblue'
	E(net)[weight >= 0]$color <- 'lightcoral'
	med <- median( abs(E(net)$weight))
	E(net)[weight > med]$color <- 'coral'
	E(net)[weight < (-1*med)]$color <- 'lightslategrey'
	par(mai=c(0,0,1,0)) 
	plot(net, vertex.label.cex=0.7,layout=layout.fruchterman.reingold,vertex.label.dist=0.2,vertex.label.color='black',vertex.frame.color='plum')

}

get_driver_overlap_matrix <- function(driverf){
	adj <- get_driver_adjacency(driverf)


	weightedOverlap <- function(v1,v2){
		a1 <- abs(v1)
		a2 <- abs(v2)

		denominator <- min( sum(a1), sum(a2))
		if(denominator == 0){
			return(0)
		}
		denominator <- max(1, denominator)
		numerator <- sum( mapply(min, a1, a2))
		wo <- numerator / denominator

		return(wo)
	}

	vectorWO <- function(v2, adjmatrix){
		return( apply(adjmatrix,1,weightedOverlap, v2))
	}

	return(apply(adj, 1, vectorWO, adj))
}

plot_driver_overlap_heatmap <- function(driverf){
	om <- get_driver_overlap_matrix(driverf)
	omdist <- 1 - om 
	tree <-  hclust( as.dist(omdist), method='average')
	dynamicMods <- cutreeDynamic(dendro=tree, distM=omdist,deepSplit=2, pamRespectsDendro=FALSE)
	dyncols <- labels2colors(dynamicMods)

	melist <- moduleEigengenes(om, colors=dyncols)
	mediss <- 1 - cor(melist$eigengenes)
	# metree <- hclust( as.dist(mediss), method='average')
	mergeheight <- 0.4
	# plot(metree)
	# abline(h=mergeheight,col='red')
	merge <- mergeCloseModules(om, dyncols, cutHeight=mergeheight)

	# plotDendroAndColors(tree,cbind(dyncols,merge$colors), c('Dynamic Tree Cut', 'Merged Dynamic'), dendroLabels=FALSE, hang=0.03)
	TOMplot(omdist, tree, merge$colors)


}