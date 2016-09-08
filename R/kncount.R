kncount <-
function(knidex, classArray)
{       
#if classArray is numeric, it represent the index of last element of 1st class (binary class with already separated), if it is a vector, it represents a class vector corresponding to knidex
if(is.matrix(knidex))
		{
        if(length(classArray)<=1)
		{
		knc=matrix(0,nrow(knidex),2)
		knc[,1]=rowSums(matrix(1,nrow(knidex),ncol(knidex))*(knidex<=classArray))
        knc[,2]=rowSums(matrix(1,nrow(knidex),ncol(knidex))*(knidex>classArray))
		}else { nClass=length(unique(classArray))
		knc=matrix(0,nrow(knidex),nClass)
		colnames(knc)=unique(classArray)
		for(i in 1:nClass)
		{
		knc[,i]= rowSums(matrix(as.numeric(classArray[knidex]==unique(classArray)[i]),nrow(knidex),ncol(knidex)))
		}
		
		}
		}
		return(knc)   
}
