SMOTE <-
function(X,target,K=5,dup_size=0)   
{
    ncD=ncol(X) #The number of attributes
	n_target=table(target)
	classP=names(which.min(n_target))
	P_set=subset(X,target==names(which.min(n_target)))[sample(min(n_target)),]     #Extract a set of positive instances
	N_set=subset(X,target!=names(which.min(n_target)))
	P_class=rep(names(which.min(n_target)),nrow(P_set))
	N_class=target[target!=names(which.min(n_target))]
    sizeP=nrow(P_set) #The number of positive instances
	sizeN=nrow(N_set) #The number of negative instances
    knear=knearest(P_set,P_set,K)
	sum_dup = n_dup_max(sizeP+sizeN,sizeP,sizeN,dup_size)
	syn_dat=NULL
    for(i in 1:sizeP)
	 {       if(is.matrix(knear))
                {    pair_idx = knear[i,ceiling(runif(sum_dup)*K)]
			    }else 
				{    pair_idx = rep(knear[i],sum_dup)
				}				
                g = runif(sum_dup)
				P_i = matrix(unlist(P_set[i,]),sum_dup,ncD,byrow=TRUE)
				Q_i = as.matrix(P_set[pair_idx,])
				syn_i = P_i + g*(Q_i - P_i)
				syn_dat = rbind(syn_dat,syn_i)

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
D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set, orig_P = P_set, K = K, K_all = NULL, dup_size = sum_dup, outcast = NULL, eps = NULL, method="SMOTE")
class(D_result) = "gen_data"
return(D_result)
}
