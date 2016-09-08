\name{knearest}
\alias{knearest}

\title{
   The function to find n_clust nearest neighbors of each instance, always removing the index of that instance if it is reported.
} 
\description{
   The function will find n_clust nearest neighbors of each instance using Fast nearest neighbors (through KD-tree method) but will correct the result if it reports the index of that instance as its neighbors.
}
\usage{
knearest(D, P, n_clust)
}

\arguments{
  \item{D}{
        a query data matrix.
}
  \item{P}{
		an input data matrix
}
  \item{n_clust}{
		the maximum number of nearest neighbors to search
}
}
\details{
	This function will perform K-nearest neighbor of instances in P on instances in P based on FNN. Then, it will verify if one of neighbors of each instance is itself then removes if it is. 
}
\value{
	The index matrix of K nearest neighbour of each instance 
}

\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}


\examples{
	data_example = sample_generator(10000,ratio = 0.80)
	P = data_example[data_example[,3]=="p",-3]
	N = data_example[data_example[,3]=="n",-3]
	D = rbind(P,N)
	knear = knearest(D,P,n_clust = 5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
