\name{kncount}
\alias{kncount}

\title{
	Counting the number of each class in K nearest neighbor
}
\description{
	The function to count how many neighbor of each instance belong to each class. 
}
\usage{
kncount(knidex, classArray)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{knidex}{
	The matrix of K nearest neighbor of dataset
}
  \item{classArray}{
	The index of last instance of the first class in the dataset or the vector containing indices of last instances of each class.
}
}
\details{
   The dataset is expected to be sorted as all m1 instances in the first class are in the first m1 instances of the dataset following with all m2 instances in the next m2 instances etc. before performing k-nearest neighbor with the knearest function. 
}
\value{
	The matrix with the number of columns equal to the number of classes. Each a[i][j] represents the number of K-nearest neighbors of i th instance belonging to the class j th
	}

\author{
 Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}

\examples{
	 D = sample_generator(1000,ratio = 0.8)
	 P = D[D[,3]=="p",]
	 N = D[D[,3]=="n",]
	 D_arr=rbind(P,N)
     knear=knearest(D_arr[,-3],P[,-3],5)
	 kncount_result = kncount(knear,nrow(P))
}