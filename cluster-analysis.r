source('network-analyzer.r')

clusterMembers <- function(clustname){
	clusteredgenes <- cbind(dnet$day$colors, rownames(days$d0))
	onecluster <- clusteredgenes[clusteredgenes[,1]==clustname,]
	return(onecluster[,2])
}

clusterMembers2 <- function(clusteredgenes, clustname){
	onecluster <- clusteredgenes[clusteredgenes[,1]==clustname,]
	return(onecluster[,2])
}


allClustsMetric <- function(aclust,bclust){
	anames <- unique(aclust[,1])
	bnames <- unique(bclust[,1])

	tani <<- matrix(0,ncol=length(anames),nrow=length(bnames))
	rownames(tani) <<- bnames
	colnames(tani) <<- anames
	over <<- matrix(0,ncol=length(anames),nrow=length(bnames))
	rownames(over) <<- bnames
	colnames(over) <<- anames
	for (aname in anames){
		for (bname in bnames){
			a <-  unique(aclust[aclust[,1]==aname,])
			b <-  unique(bclust[bclust[,1]==bname,])
			alen <- nrow(a)
			blen <- nrow(b)

			int <- length( intersect(a,b))
			uni <- length( union(a,b))
			min <- min(alen,blen)

			tanimoto <- int/uni
			overlap <- int/min


			if(tanimoto > 0.5 || overlap > 0.5){
				print( paste(aname,bname))
				print( paste('tanimoto: ',tanimoto))
				print( paste('overlap: ',overlap))
				ia <- int/alen
				ib <- int/blen
				print( paste('i / a: ', ia))
				print( paste('i / b: ', ib))
				cat('\n')
				# print( paste('intersect: ',int))
				# print( paste('union: ',uni))
				# print( paste('alen: ',alen))
				# print( paste('blen: ', blen))

			}
		}
	}

}

# if(!exists('d0.gclust.stat')){
# 	plotTOMDendo(day='d0')
# 	d0.gclust.stat <- cbind(staticcut, rownames(days$d0))
# 	d0.gclust.dyna <- cbind(dynamiccut, rownames(days$d0))

# 	plotTOMDendo(day='d1')
# 	d1.gclust.stat <- cbind(staticcut, rownames(days$d0))
# 	d1.gclust.dyna <- cbind(dynamiccut, rownames(days$d0))

# 	plotTOMDendo(day='d3')
# 	d3.gclust.stat <- cbind(staticcut, rownames(days$d0))
# 	d3.gclust.dyna <- cbind(dynamiccut, rownames(days$d0))

# 	plotTOMDendo(day='d6')
# 	d6.gclust.stat <- cbind(staticcut, rownames(days$d0))
# 	d6.gclust.dyna <- cbind(dynamiccut, rownames(days$d0))

# 	plotTOMDendo(day='d14')
# 	d14.gclust.stat <- cbind(staticcut, rownames(days$d0))
# 	d14.gclust.dyna <- cbind(dynamiccut, rownames(days$d0))

# 	plotTOMDendo(day='d21')
# 	d21.gclust.stat <- cbind(staticcut, rownames(days$d0))
# 	d21.gclust.dyna <- cbind(dynamiccut, rownames(days$d0))

# 	plotTOMDendo(day='d28')
# 	d28.gclust.stat <- cbind(staticcut, rownames(days$d0))
# 	d28.gclust.dyna <- cbind(dynamiccut, rownames(days$d0))
# }

# print('static d0 -> d1')
# allClustsMetric(d0.gclust.stat, d1.gclust.stat)
# cat('\n')
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)
# # print('dynamic d0 -> d1')
# # allClustsMetric(d0.gclust.dyna, d1.gclust.dyna)
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)

# print('static d1 -> d3')
# allClustsMetric(d1.gclust.stat, d3.gclust.stat)
# cat('\n')
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)
# # print('dynamic d1 -> d3')
# # allClustsMetric(d1.gclust.dyna, d3.gclust.dyna)
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)

# print('static d3 -> d6')
# allClustsMetric(d3.gclust.stat, d6.gclust.stat)
# cat('\n')
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)
# # print('dynamic d3 -> d6')
# # allClustsMetric(d3.gclust.dyna, d6.gclust.dyna)
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)

# print('static d6 -> d14')
# allClustsMetric(d6.gclust.stat, d14.gclust.stat)
# cat('\n')
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)
# # print('dynamic d6 -> d14')
# # allClustsMetric(d6.gclust.dyna, d14.gclust.dyna)
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)

# print('static d14 -> d21')
# allClustsMetric(d14.gclust.stat, d21.gclust.stat)
# cat('\n')
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)
# # print('dynamic d14 -> d21')
# # allClustsMetric(d14.gclust.dyna, d21.gclust.dyna)
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)

# print('static d21 -> d28')
# allClustsMetric(d21.gclust.stat, d28.gclust.stat)
# cat('\n')
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)
# # allClustsMetric(d21.gclust.dyna, d28.gclust.dyna)
# # print('dynamic d21 -> d28')
# # print('tanimoto')
# # print(tani)
# # print('over')
# # print(over)

# print('static d0 -> d21')
# allClustsMetric(d0.gclust.stat, d21.gclust.stat)
# cat('\n')

# print('static d0 -> d28')
# allClustsMetric(d0.gclust.stat, d28.gclust.stat)
# cat('\n')



