# cluster 'data' into 'nclust' clusters

plotRatOne <- function(ngenes=1000,nclust=5){
	if(!exists('inf.rat1')){
		source('load-data.r')
		getRatOneVariableGenesWithInferredDays(ngenes)
	}
	cnum <- -1
	baseline <- inf.rat1[,1]
	baseline <- matrix(baseline,byrow=FALSE,nrow=nrow(inf.rat1),ncol=ncol(inf.rat1))
	fold <- inf.rat1 / baseline
	logfold <- log(fold,2)
	data <- logfold
	clusterAndPlotMedoidsAndSeries(data,nclust)
}

clusterAndPlotAllSeries <- function(series,nclust=5){
	library(cluster)
	library(dtw)
	c.series <- pam(series,nclust,metric='dtw')

	time <- 1:dim(series)[2]
	nseries <- dim(series)[1]

	plot(range(time), range(series), type='n', xlab="Day", ylab="Exp. Level")
	colors <- rainbow( nclust)

	for (i in 1:nseries){
		iclust <- c.series$cluster[i] # clustering$cluster[i]
		if(cnum == -1 || is.na(cnum) || cnum == iclust ){
			lines(time, series[i,], col=colors[iclust])

		}
	}

	title("Expression level by day")
}

clusterAndPlotMedoids <- function(series,nclust=5){
	library(cluster)
	library(dtw)
	clustering <- pam(series,nclust,metric='dtw')

	time <- 1:dim(series)[2]

	plot(range(time), range(series), type='n', xlab="Day", ylab="Exp. Level")
	colors <- rainbow( nclust)

	for (i in 1:nclust){
		lines(time, clustering$medoids[i,], lty=2, col=colors[i])
	}

	seriesByClust <- list()

	for (c in 1:nclust){
		seriesByClust[[c]] <- matrix(0.0, 0, dim(series)[2])
	}
	augseries <- cbind(rownames(series),clustering$clustering,series)
	rf <- function(r) {
		genename <- r[1]
		clustNum <- r[2]
		oldmatrix <- seriesByClust[[clustNum]]
		newrow <- as.numeric(r[3:length(r)])
		newmatrix <- rbind(oldmatrix, newrow)
		rownames(newmatrix) <- c(rownames(oldmatrix), genename)
		seriesByClust[[clustNum]] <<- newmatrix

	} 
	apply(augseries,1,rf)

	extremesByClust <- list()
	for (c in 1:nclust){
		extremesByClust[[c]] <- matrix(0.0, 2, dim(series)[2])
	}

	for (c in 1:nclust){
		clust <- seriesByClust[[as.character(c)]]
		for (t in time){
			clust_t <- clust[,t]
			extremesByClust[[c]][1,t] <- maxval <- max(clust_t)
			extremesByClust[[c]][2,t] <- minval <- min(clust_t)
			rand <- runif(1,-0.1,0.1)
			if(FALSE && t %in% c(0,1,3,6,14,21,28)){
				lines(c(t,t),c(minval,maxval),lty=2,col=colors[c])
			}
		}
		lines(time,extremesByClust[[c]][1,],col=colors[c])
		lines(time,extremesByClust[[c]][2,],col=colors[c])
	}


	title("Central expression level by day")

	return(seriesByClust)

}

clusterAndPlotMedoidsAndSeries <- function(series,nclust=5){
	library(cluster)
	library(dtw)
	c.series <- pam(series,nclust,metric='dtw')

	time <- 1:dim(series)[2]
	nseries <- dim(series)[1]

	plot(range(time), range(series), type='n', xlab="Day", ylab="Exp. Level")
	colors <- rainbow( nclust)

	for (i in 1:nseries){
		iclust <- c.series$cluster[i] # clustering$cluster[i]
		if(cnum == -1 || is.na(cnum) || cnum == iclust ){
			lines(time, series[i,], col=colors[iclust])

		}
	}


	for (i in 1:nclust){
		lines(time, c.series$medoids[i,], col="BLACK")
	}

	title("Expression level by day")


}

plotNumClusters <- function(data){
	wssplot(data)                                                #2
	library(NbClust)
	set.seed(1234)
	nc <- NbClust(data, min.nc=2, max.nc=15, method="kmeans")
 	barplot(table(nc$Best.n[1,]),
          xlab="Numer of Clusters", ylab="Number of Criteria",
          main="Number of Clusters Chosen by 26 Criteria")
}



wssplot <- function(data, nc=15, seed=1234){
               wss <- (nrow(data)-1)*sum(apply(data,2,var))
               for (i in 2:nc){
                    set.seed(seed)
                    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
                plot(1:nc, wss, type="b", xlab="Number of Clusters",
                     ylab="Within groups sum of squares")}

inferPoints <- function(data){
	library(zoo)
	rowinterp <- function(r){
		#                   d0   d1   d2 d3   d4 d5 d6  
		infr <- na.spline(c(r[1],r[2],NA,r[3],NA,NA,r[4],
		#                   d7 d8 d9 d10 d11 d12 d13 d14
		                    NA,NA,NA,NA, NA, NA, NA, r[5],
		#                   d15 d16 d17 d18 d19 d20 d22
		                    NA, NA, NA, NA, NA, NA, r[6],
		#                   d23 d24 d25 d26 d27 d28                    
		                    NA, NA, NA, NA, NA, r[7]))
		return(infr)
	}

	apped <- apply(data,1,FUN=rowinterp)
	return( matrix(apped,byrow=TRUE,nrow=dim(data)[1]))
}