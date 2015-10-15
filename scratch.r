

for (i in 1:50){
	for(j in 1:50){
		if( i!= j){
			if( catcher$P[i,j] <= 0.05){
			 	sp <- catcher$r[i,j]
			 	if(sp < 0.4){
			 		adja[i,j] <- -1
			 	} else if(sp > 0.4){
			 		adja[i,j] <- 1
			 	}
			}
		}
	}
}