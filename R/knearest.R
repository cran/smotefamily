knearest <-
function(D,P,n_clust)
{
    if(!requireNamespace("FNN", quietly = TRUE))
    {
       install.packages("FNN", quiet = TRUE)
    }
	  knD<-FNN::knnx.index(D,P,k=(n_clust+1), algo="kd_tree")
	  knD=knD*(knD!=row(knD))
	  que=which(knD[,1]>0)
	  for (i in que)
	  {    knD[i,which(knD[i,]==0)]=knD[i,1]
	       knD[i,1]=0
		   #print(i)
	  }
      return(knD[,2:(n_clust+1)])

}
