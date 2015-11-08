
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

getNorms <- function(fname){
	norms <- getClipNorms(fname, -1) # 300 for genus, 4,000 for nogs

	# sort norms by day then rat
	m <- do.call(rbind, strsplit(colnames(norms), "[.]"))[,-1]
	daysort <- order(as.numeric(do.call(rbind,strsplit(m[,1],'day'))[,2]))
	norms <- norms[,daysort]

	m <- do.call(rbind, strsplit(colnames(norms), "[.]"))[,-1]
	ratsort <- order(do.call(rbind,strsplit(m[,2],'_count'))[,1])
	norms <- norms[,ratsort]

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

analyze_drivercorr <- function(){
	dnames <- getDataNames()
	dnames <- dnames
	onames <- getOutNames('hh_results/driver_corr/','.norm.diseaseday_adjusted.driver_corr.csv')
	onames <- onames

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

		intersects <- c()
		slopes <- c()
		sps <- c()
		ips <- c()
		inds <- c()
		deps <- c()

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

			if( length(ind[ind==0]) <= 8 && length(dep[dep==0]) <= 8){ # limit regression to data with enough non-zero values
					fit <- lm(dep ~ ind)
					slopepval <- summary(fit)$coefficients[2,4] # easier way to get this?
					interceptpval <- summary(fit)$coefficients[1,4] # easier way to get this?
					if(!is.na(slopepval) && slopepval < 0.05){ 
						intersects <<- append(intersects, summary(fit)$coefficients[1,1])
						slopes <<- append(slopes, summary(fit)$coefficients[2,1])
						sps <<- append(sps, slopepval)
						ips <<- append(ips, interceptpval)
						inds <<- append(inds, indname)
						deps <<- append(deps, depname)
					}
				}
		}

		apply(delall,1,depfunc)

		# for (i in 1:nrow(r1)){
		# 	if( i %% 100 == 0){
		# 		print( paste('row:',i))
		# 	}
		# 	for (j in 1:nrow(r1)){
		# 		ind <- baseall[j,]
		# 		dep <- delall[i,]
		# 		if( length(ind[ind==0]) < 9 && length(dep[dep==0]) < 9){ # limit regression to data with enough non-zero values
		# 			fit <- lm(ind ~ dep)
		# 			slopepval <- summary(fit)$coefficients[2,4] # easier way to get this?
		# 			interceptpval <- summary(fit)$coefficients[1,4] # easier way to get this?
		# 			if(!is.na(slopepval) && slopepval < 0.05){ 
		# 				intersects <- append(intersects, summary(fit)$coefficients[1,1])
		# 				slopes <- append(slopes, summary(fit)$coefficients[2,1])
		# 				sps <- append(sps, slopepval)
		# 				ips <- append(ips, interceptpval)
		# 				inds <- append(inds, rownames(r1)[j])
		# 				deps <- append(deps, rownames(r1)[i])
		# 			}
		# 		}
		# 	}
		# }

		timedependencies <- data.frame(inds, deps, slopes, intersects, sps, ips)
		write.csv(timedependencies,oname)
	}
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

	names <- c('~/Dropbox/hh_data/genus.diamond.aggregated.counts.norm.matrix')

	names <- c('/Users/dcdanko/Science/Meta-Omics/hh_results/diseaseday_corr/genus.norm.diseaseday_adjusted.matrix.csv')
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

get_driver_matrix <- function(fname){
	drives <- read.csv(fname)[,-1]

	m <- matrix(nrow=0,ncol=0)
	matrixmaker <- function(driver){
		#print(driver)
		ind <- driver[1]
		dep <- driver[2]
		slope <- as.numeric(driver[3])
		p <- as.numeric(driver[5])
		if(p<0.01){
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
	heatmap.2(m,Rowv=NA,Colv=NA,col=cm.colors(256),na.color='gray',margins=c(10,10),trace='none')
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
		t <- tpos - tneg
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
		t <- tpos - tneg
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

	barplot(t(indm[1:25,3:4]), col=cm.colors(2), las=2, axisnames=FALSE)
	barplot(t(depm[1:25,3:4]), col=cm.colors(2), las=2, axisnames=FALSE)



}