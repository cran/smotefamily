RSLS = function(X,target,K=5,C=5,dupSize=0)
{
    ncD=ncol(X) #The number of attributes
	n_target=table(target)
	classP=names(which.min(n_target))
	P_set=subset(X,target==names(which.min(n_target)))[sample(min(n_target)),]     #Extract a set of positive instances
	N_set=subset(X,target!=names(which.min(n_target)))
	P_class=rep(classP,nrow(P_set))
	N_class=target[target!=names(which.min(n_target))]
    sizeP=nrow(P_set) #The number of positive instances
	sizeN=nrow(N_set) #The number of negative instances
	Darr=rbind(P_set,N_set)  #The re-arranged data with positive first
#main

	knear_D=knearest(Darr,P_set,C)
	knear_P=knearest(P_set,P_set,K)
	knct=kncount(knear_D,sizeP)
    safe_level=knct[,1]
    sum_dup = n_dup_max(sizeP+sizeN,length(safe_level>0),sizeN,dupSize)
	Outcast=P_set[knct[,1]==0,]  
	syn_dat=NULL
	for (i in 1:sizeP)
	{ 
		slp=safe_level[i] #which k neighbour of positive i is positive
	    if(slp>0)
		{
			      pair_idx=knear_P[i,ceiling(runif(sum_dup)*K)]
				  rbind(knear_D[i,],knear_D[pair_idx,])->all_neig
				  neg_to_check = unique(all_neig[all_neig>sizeP])
	              sln=safe_level[pair_idx] #which k neighbor of random one of positive i is positive
			
			g = mapply(gap,sl_p=rep(slp,sum_dup),sl_n=sln)
			P_i = matrix(unlist(P_set[i,]),sum_dup,ncD,byrow=TRUE)
			Q_i = as.matrix(P_set[pair_idx,])
			temp_i = P_i + g*(Q_i - P_i)
			if(length(neg_to_check)>0)
			{
				for(j in 1:sum_dup)
				{
					dist_temp_toP = min(dist(rbind(temp_i[j,],unlist(P_set[i,]),unlist(P_set[pair_idx[j],]))))			
					dist_chk=matrix(dist(rbind(temp_i[j,],Darr[neg_to_check,])))[1:length(neg_to_check)] #check if 
					while(sum(dist_chk<dist_temp_toP)>0)
					{
			            if(slp/sln[j]>=1)
						{	
							start_j = unlist(P_set[i,])
							end_j = temp_i[j,]				
						} else
						{
							start_j = temp_i[j,]
							end_j = unlist(P_set[pair_idx[j],])
						}
						temp_i[j,] = start_j + gap()*(end_j-start_j)
						dist_temp_toP = min(dist(rbind(temp_i[j,],unlist(P_set[i,]),unlist(P_set[pair_idx[j],]))))			
						dist_chk=matrix(dist(rbind(temp_i[j,],Darr[neg_to_check,])))[1:length(neg_to_check)]
						#cat("*") #improve activation
					}
				}
			}
			syn_dat = rbind(syn_dat,temp_i)
		}
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

D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set, orig_P = P_set, K = K, K_all = C, dup_size = sum_dup, outcast = Outcast_df, eps=NULL, method = "RSLS")
class(D_result) = "gen_data"		
		
print("RSLS is done")
return(D_result)
}
