DBSMOTE =function(X,target,dupSize=0,MinPts=NULL,eps=NULL)
{
#rm(SynP)
  
	ncD=ncol(X) #The number of attributes
	n_target=table(target)
	classP=names(which.min(n_target))
	P_set=subset(X,target==names(which.min(n_target)))[sample(min(n_target)),]     #Extract a set of positive instances
	N_set=subset(X,target!=names(which.min(n_target)))
	P_class=rep(names(which.min(n_target)),nrow(P_set))
	N_class=target[target!=names(which.min(n_target))]
    sizeP=nrow(P_set) #The number of positive instances
	sizeN=nrow(N_set) #The number of negative instances
	set.seed(ceiling(runif(1)*10000))
	if(is.null(MinPts))
	{
	MinPts=2+ceiling(log(sizeP,10))-1
	}
	if(is.null(eps))
	{
		FNN::get.knn(P_set, k=MinPts+1, algorithm="kd_tree")->Pknn 
		quantile(Pknn$nn.dist[,MinPts+1],0.75)->eps
	}
	
	ds <- dbscan::dbscan(P_set,eps,MinPts)
	size_Pgen = nrow(P_set[ds$cluster!=0,])
	Outcast = P_set[ds$cluster==0,]

    sum_dup=n_dup_max(sizeP+sizeN,size_Pgen,sizeN,dupSize)
    
	clus=sort(unique(ds$cluster)[unique(ds$cluster)>0])

syn_dat=NULL	
for(i in clus)
{
   Vert_Cl=which(ds$cluster==i)
   Adj_m_Cl=matrix(0,length(Vert_Cl),length(Vert_Cl))    # a matrix i x j -> ith positive is directly density reachable to jith
   colnames(Adj_m_Cl)=Vert_Cl
   rownames(Adj_m_Cl)=Vert_Cl
   as.matrix(dist(P_set[Vert_Cl,]))->distP
   NNq=rowSums(distP<=eps)
   for(j in 1:length(Vert_Cl))
   {
		if(length(Vert_Cl)>=MinPts+1)
		{ 
			Adj_m_Cl[j,]=(distP[j,]<=eps&distP[j,]>0)*(NNq>=MinPts)
		}
		else
		{
			Adj_m_Cl[j,]=(distP[j,]<=eps&distP[j,]>0)
		}
	}
   Centroid=colMeans(P_set[Vert_Cl,])
   pDCinxD=order(as.matrix(dist(rbind(P_set[Vert_Cl,],Centroid)))[,length(Vert_Cl)+1])[2]
   pseudo_Cen=P_set[Vert_Cl[pDCinxD],]

#install.packages("igraph")
	#suppressWarnings(suppressMessages(library("igraph")))
	igraph::graph.adjacency(Adj_m_Cl)->graph_clus

	igraph::get.shortest.paths(graph_clus, pDCinxD, mode = "all", weights = subset(as.vector(t(distP*Adj_m_Cl)),as.vector(t(distP*Adj_m_Cl))>0))->short_Path

	
	syn_dat_cl=NULL
	for(j in 1:length(Vert_Cl))
	{   
        path=short_Path$vpath[[j]]
		if(length(path)>1)
		{   
			print(length(path))
			rand_edge=ceiling(runif(sum_dup)*(length(path)-1))
			g = runif(sum_dup)
			P_i=P_set[Vert_Cl[path[rand_edge]],]
			Q_i= P_set[Vert_Cl[path[rand_edge+1]],]
            syn_i = P_i + g*(Q_i - P_i)
			syn_dat_cl = rbind(syn_dat_cl,syn_i)
		}
	}	
	syn_dat=rbind(syn_dat,syn_dat_cl)
}

P_set[,ncD+1] = P_class		
colnames(P_set)=c(colnames(X),"class")
N_set[,ncD+1] = N_class
colnames(N_set)=c(colnames(X),"class")
rownames(syn_dat)= NULL
syn_dat = data.frame(syn_dat)
syn_dat[,ncD+1] = rep(names(which.min(n_target)),nrow(syn_dat))
colnames(syn_dat)=c(colnames(X),"class")		
NewD=rbind(P_set,syn_dat,N_set)
rownames(NewD) = NULL
rownames(Outcast) <- c()
new_target_out=rep(classP,nrow(as.matrix(Outcast)))
Outcast_df=data.frame(Outcast)
Outcast_df[,ncD+1]=new_target_out
colnames(Outcast_df)=c(colnames(X),"class")

D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set, orig_P = P_set, K = NULL, K_all = NULL, dup_size = sum_dup, outcast = Outcast_df, eps=eps, method = "DBSMOTE")
class(D_result) = "gen_data"		
print("DBSMOTE is Done")
return(D_result)
}
