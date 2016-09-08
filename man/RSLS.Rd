\name{RSLS}
\alias{RSLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Relocating Safe-level SMOTE
}
\description{
	Generate synthetic positive instances using Relocating Safe-level SMOTE algorithm. Using the parameter "Safe-Level" to determine the possible location and relocating synthetic instances if there is too close to majority instances. 
}
\usage{
RSLS(X, target, K = 5, C = 5, dupSize = 0)
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
	\item{method}{The name of oversampling method used for this generated dataset (RSLS)}
}
\references{
	Siriseriwan, W. and Sinapiromsaran, K. The Effective Redistribution for Imbalance Dataset : Relocating Safe-Level SMOTE with Minority Outcast Handling. Chiang Mai Journal of Science. 43(1), 234 - 246.
}
\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}

\examples{
	library(smotefamily)
	data_example = sample_generator(5000,ratio = 0.80)
    genData = RSLS(data_example[,-3],data_example[,3])
	genData_2 = RSLS(data_example[,-3],data_example[,3],K=7, C=5)


}
