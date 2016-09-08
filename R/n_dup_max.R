n_dup_max <-
function(size_input,size_P,size_N,dup_size=0)
{   #Size_P is the number of positive used for generating not actual size of P
	if(is.vector(dup_size)&&length(dup_size)>1)
	{
		if(length(which(dup_size==0))>0)
			{
				sizeM=floor((2*size_N-size_input)/size_P)
			}
		if(length(which(dup_size==0))==0)
			{ 
				sizeM=max(dup_size)
			}
	}
	if(!is.vector(dup_size)||length(dup_size)==1)
	{
		if(dup_size==0)
			{
				sizeM=floor((2*size_N-size_input)/size_P)
			}
	if(dup_size!=0)
			{
				sizeM=dup_size
			}
	}
    return(sizeM)	
}
