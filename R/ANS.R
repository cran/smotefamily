ANS <-
function(X,target,dupSize=0)
{          
    nrD=nrow(X)
    ncD=ncol(X) #The number of attributes
	n_target=table(target)
	classP=names(which.min(n_target))
	P_set=subset(X,target==classP)[sample(min(n_target)),]     #Extract a set of positive instances
	N_set=subset(X,target!=classP)			#Extract a set of negative instances  
	P_class=rep(classP,nrow(P_set))
	N_class=target[target!=classP]
    sizeP=nrow(P_set) #The number of positive instances
	sizeN=nrow(N_set) #The number of negative instances
	Darr=rbind(P_set,N_set)  #The re-arranged data with positive first
#main
   
    OClevel=rep(0,sizeP)

	if (requireNamespace("FNN", quietly = TRUE)) {
	knearD<-FNN::knnx.index(Darr,P_set,k=ceiling(0.25*nrD)+1, algo="kd_tree")[,-1]
	}
    P_idx=1:sizeP
	i=1
	O_Cl=0
	while(length(P_idx)>0&i<30)
	{  
	   temp=which(knearD[,i]<=sizeP&OClevel==0)
	   #print(temp)
	   OClevel[temp]=rep(i,length(temp))
	   #print(OClevel)
	   O_Cl[i]=sum(OClevel==0)/sizeP
	   P_idx=which(OClevel==0)
	   i=i+1
	   #print(i)
	}
	
	tol = 0.01
	
	#print(OCl)
	if(length(which(abs(diff(O_Cl[-1]))<tol))>=1)
	{
	Cvalue=which(abs(diff(O_Cl))<tol)[1]+1
	}else{   Cvalue=5
	}


	k_for_all=rep(NA,sizeP)
	#print(paste(" C is",Cvalue))
	#print(Cvalue)
    Pgen=P_set[OClevel<Cvalue,]
	Pgenidx=which(OClevel<Cvalue)
    Outcast=P_set[OClevel>=Cvalue,]   
    #print(Pgen)
	if (requireNamespace("FNN", quietly = TRUE)) {
	FNN::get.knn(Pgen, k=(nrow(Pgen)-1), algorithm="kd_tree")->knn_ofP 
    }
    eps=knn_ofP$nn.dist[which.max(knn_ofP$nn.dist[,1]),2]
    k_for_pgen=rowSums(knn_ofP$nn.dist<=eps)   
	k_for_all[Pgenidx]=k_for_pgen
    knearP=knn_ofP$nn.index
	
##############Set the number of rounds to generate, 0 for generating till balanced###########

	sum_dup =n_dup_max(nrD,length(Pgenidx),sizeN,dupSize)
###############################################################################################

	syn_dat=c()
	for(i in 1:length(Pgenidx))
	 {      pair_idx=NULL 
			if(is.matrix(knearP))
                {    pair_idx = knearP[i,ceiling(runif(sum_dup)*k_for_pgen[i])]
			    }else 
				{    pair_idx = rep(knearP[i],sum_dup)
				}				
                g = runif(sum_dup)
				P_i = matrix(unlist(Pgen[i,]),sum_dup,ncD,byrow=TRUE)
				Q_i = as.matrix(Pgen[pair_idx,])
				syn_i = P_i + g*(Q_i - P_i)
				syn_dat = rbind(syn_dat,syn_i)
			#  print(i)
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

D_result = list(data = NewD, syn_data = syn_dat, orig_N = N_set, orig_P = P_set, K = k_for_all, K_all = Cvalue, dup_size=sum_dup , outcast = Outcast_df, eps = eps,method = "ANS")
class(D_result) = "gen_data"	

print("ANS is done")
return(D_result)



}
