sample_generator =function(n, ratio = 0.8, xlim = c(0,1), ylim = c(0,1),radius = 0.25,overlap = -0.05, outcast_ratio=0.01)
{

seed=runif(n)
center = c(xlim[1] + 0.5*diff(xlim),ylim[1] + 0.5*diff(ylim))
x=matrix(0,n,2)
cx=matrix(1,n,1)
for(i in 1:length(seed))
{
if(seed[i]<(1-ratio))
{   r = rep(0,2)
	r[1] = runif(1,min = 0, max = radius*diff(xlim))
	r[2] = runif(1,min = 0, max = radius*diff(ylim))
	theta = runif(2,max = 2)*acos(-1)
	x[i,1]= center[1] + r[1]*cos(theta[1])
    x[i,2]= center[2] + r[2]*cos(theta[2])
    cx[i,1]=2
}else if(seed[i]>(1-ratio))
{
    x[i,1:2]=runif(2,0,1)
	while(dist(rbind(x[i,1:2],center))<(radius+overlap)|dist(rbind(x[i,1:2],center))<(radius+overlap))
	{
	 x[i,1:2]=runif(2,0,1)
	}
	if(runif(1,0,1)<outcast_ratio)
	{    cx[i,1]=2
	}else{
	cx[i,1]=1
	}
}
	}
x=cbind(x,cx) 
result=matrix("n",n,1)
result[x[,3]==2]="p"
data.frame(x[,1:2])->dx
dx=cbind(dx,result)
dx[,3]=as.character(dx[,3])
return(dx)
}
