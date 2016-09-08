gap <-
function(sl_p=1,sl_n=1)
{ 
    if(sl_n==0)
	{	 
		return(0)
	}
	else if(sl_p==sl_n)
		{
			return(runif(1,0,1))
		}
		else if(sl_p>sl_n)
			{
			    return(runif(1,0,sl_n/sl_p))
			}
			else
			{
				return(runif(1,1-(sl_p/sl_n),1))
			}
}
