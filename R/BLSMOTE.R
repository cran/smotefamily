BLSMOTE <-
function(X,target,K=5,C=5,dupSize=0,method =c("type1","type2"))
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
	Darr=rbind(P_set,N_set)  #The re-arranged data with positive first


	knear_D=knearest(Darr,P_set,C)
	knct=kncount(knear_D,sizeP)
    safe_level=knct[,1]/C
	Danger_idx = which(safe_level<=0.5&safe_level>0)
	Danger = P_set[Danger_idx,]
	SAFE = P_set[safe_level>0.5,]
	Outcast = P_set[safe_level==0,]
	sum_dup = n_dup_max(sizeP+sizeN,nrow(Danger),sizeN,dupSize)
	
	syn_dat=NULL
	if(match.arg(method) == "type1")
    {	
		if(K>=nrow(Danger))
		{   stop("K cannot exceed the size of DANGER")
		}	
		knear_P=knearest(Danger,Danger,K)
		for (i in 1:nrow(Danger))
		{	    	       
			    pair_idx=knear_P[i,ceiling(runif(sum_dup)*K)]		
				g = mapply(gap,sl_p=rep(1,sum_dup),sl_n=rep(1,sum_dup))
				P_i = matrix(unlist(Danger[i,]),sum_dup,ncD,byrow=TRUE)
				Q_i = as.matrix(Danger[pair_idx,])
				syn_i = P_i + g*(Q_i - P_i)
				syn_dat = rbind(syn_dat,syn_i)
			#print(i)
		
		}
	}
	else if(match.arg(method) == "type2") 
		{
			for (i in Danger_idx)
			{	    
					pair_idx=knear_D[i,ceiling(runif(sum_dup)*C)]		
					g = rep(0,sum_dup)
					for(j in 1:sum_dup)
					{   if(pair_idx[j]<=sizeP)
						{  g[j] = runif(1) }
						else
						{   g[j] = runif(1,min = 0, max = 0.5) 
						}
					}
					P_i = matrix(unlist(P_set[i,]),sum_dup,ncD,byrow=TRUE)
					Q_i = as.matrix(P_set[pair_idx,])
					syn_i = P_i + g*(Q_i - P_i)
					syn_dat = rbind(syn_dat,syn_i)
			#print(i)
				

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
D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set, orig_P = P_set, K = K, K_all = C, dup_size= sum_dup, outcast = NULL, eps=NULL, method = paste("Borderline-SMOTE",match.arg(method),sep=" "))
class(D_result) = "gen_data"

print("Borderline-SMOTE done")
return(D_result)
}
