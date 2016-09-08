\name{SLS}
\alias{SLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Safe-level SMOTE
}
\description{
	Generate synthetic positive instances using Safe-level SMOTE algorithm. Using the parameter "Safe-level" to determine the possible location of synthetic instances.
}
\usage{
SLS(X, target, K = 5, C = 5, dupSize = 0)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{C}{
	The number of nearest neighbors during calculating safe-level process 
}
  \item{dupSize}{
	The number or vector representing the desired times of synthetic minority instances over the original number of majority instances
}
 
}

\value{
	\item{data}{A resulting dataset consists of original minority instances, synthetic minority instances and original majority instances with a vector of their respective target class appended at the last column}
	\item{syn_data}{A set of synthetic minority instances with a vector of minority target class appended at the last column}
	\item{orig_N}{A set of original instances whose class is not oversampled with a vector of their target class appended at the last column}
	\item{orig_P}{A set of original instances whose class is oversampled with a vector of their target class appended at the last column}
	\item{K}{The value of parameter K for nearest neighbor process used for generating data}
	\item{K_all}{The value of parameter C for nearest neighbor process used for calculating safe-level}
	\item{dup_size}{The maximum times of synthetic minority instances over original majority instances in the oversampling}
	\item{outcast}{A set of original minority instances which has safe-level equal to zero and is defined as the minority outcast}
	\item{eps}{Unavailable for this method}
	\item{method}{The name of oversampling method used for this generated dataset (SLS)}
}
\references{
	Bunkhumpornpat, C., Sinapiromsaran, K. and Lursinsap, C. 2009. Safe-level-SMOTE: Safe-level-synthetic minority oversampling technique for handling the class imbalanced problem. Proceedings of the 13th Pacific-Asia Conference on Advances in Knowledge Discovery and Data Mining. 2009, 475-482.
}
\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}

\examples{
    data_example = sample_generator(5000,ratio = 0.80)
	genData = SLS(data_example[,-3],data_example[,3])
	genData_2 = SLS(data_example[,-3],data_example[,3],K=7, C=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
