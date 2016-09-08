\name{DBSMOTE}
\alias{DBSMOTE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Density-based SMOTE
}
\description{
	Generate a oversampling dataset from imbalance dataset using Density-based SMOTE. Using density reachability concept to cluster minority instances and generate synthetic instances.
}
\usage{
DBSMOTE(X, target, dupSize = 0, MinPts = NULL, eps = NULL)
}

\arguments{
  \item{X}{
	A data frame or matrix of numeric-attributed dataset 
}
  \item{target}{
	A vector of a target class attribute 
}
  \item{dupSize}{
	A number of vector representing the desired times of synthetic minority instances over the original number of majority instances
}
  \item{MinPts}{
	The minimum instance parameter to decide whether each instance inside eps is reachable, the automatic algorithm is used to find the value instead if there is no positive integer value given for it.
}
  \item{eps}{
	The radius to consider neighbor.
}

}

\value{
	\item{data}{A resulting dataset consists of original minority instances, synthetic minority instances and original majority instances with a vector of their respective target class appended at the last column}
	\item{syn_data}{A set of synthetic minority instances with a vector of minority target class appended at the last column}
	\item{orig_N}{A set of original instances whose class is not oversampled with a vector of their target class appended at the last column}
	\item{orig_P}{A set of original instances whose class is oversampled with a vector of their target class appended at the last column}
	\item{K}{Unavailable for this method}
	\item{K_all}{Unavailable for this method}
	\item{dup_size}{The maximum times of synthetic minority instances over original majority instances in the oversampling}
	\item{outcast}{A set of original minority instances which is defined as NOISE/minority outcast}
	\item{eps}{The value of parameter eps}
	\item{method}{The name of oversampling method used for this generated dataset}
}
\references{
Bunkhumpornpat, C., Sinapiromsaran, K. and Lursinsap, C. 2012. DBSMOTE: Density-based synthetic minority oversampling technique. Applied Intelligence. 36, 664-684.
}
\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}




\examples{
data_example = sample_generator(5000,ratio = 0.90)
genData = DBSMOTE(data_example[,-3],data_example[,3])
}

