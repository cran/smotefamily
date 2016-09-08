\name{ADASYN}
\alias{ADAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Adaptive Synthetic Sampling Approach for Imbalanced Learning
}
\description{
	Generate synthetic positive instances using ADASYN algorithm. The number of majority neighbors of each minority instance determines the number of synthetic instances generated from the minority instance.
}
\usage{
ADAS(X,target,K=5)  
}
\arguments{
  \item{X}{
	A data frame or matrix of numeric-attributed dataset 
}
  \item{target}{
	A vector of a target class attribute corresponding to a dataset X.
}
  \item{K}{
	The number of nearest neighbors during sampling process
}

}
\value{
	\item{data}{A resulting dataset consists of original minority instances, synthetic minority instances and original majority instances with a vector of their respective target class appended at the last column}
	\item{syn_data}{A set of synthetic minority instances with a vector of minority target class appended at the last column}
	\item{orig_N}{A set of original instances whose class is not oversampled with a vector of their target class appended at the last column}
	\item{orig_P}{A set of original instances whose class is oversampled with a vector of their target class appended at the last column}
	\item{K}{The value of parameter K for nearest neighbor process used for generating data}
	\item{K_all}{Unavailable for this method}
	\item{dup_size}{A vector of times of synthetic minority instances over original majority instances in the oversampling in each instances}
	\item{outcast}{Unavailable for this method}
	\item{eps}{Unavailable for this method}
	\item{method}{The name of oversampling method used for this generated dataset (ADASYN)}
}
\references{
	He, H., Bai, Y., Garcia, E. and Li, S. 2008. ADASYN: Adaptive synthetic sampling approach for imbalanced learning. Proceedings of IJCNN 2008. (IEEE World Congress on Computational Intelligence). IEEE International Joint Conference. pp.1322-1328.
}
\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}

\examples{
data_example = sample_generator(10000,ratio = 0.80)
genData = ADAS(data_example[,-3],data_example[,3])
genData_2 = ADAS(data_example[,-3],data_example[,3],K=7)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

