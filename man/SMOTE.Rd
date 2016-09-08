\name{SMOTE}
\alias{SMOTE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Synthetic Minority Oversampling TEchnique
}
\description{
	Generate synthetic positive instances using SMOTE algorithm
}
\usage{
SMOTE(X, target, K = 5, dup_size = 0)
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
  \item{dup_size}{
	The number or vector representing the desired times of synthetic minority instances over the original number of majority instances
}
  
}

\value{
	\item{data}{A resulting dataset consists of original minority instances, synthetic minority instances and original majority instances with a vector of their respective target class appended at the last column}
	\item{syn_data}{A set of synthetic minority instances with a vector of minority target class appended at the last column}
	\item{orig_N}{A set of original instances whose class is not oversampled with a vector of their target class appended at the last column}
	\item{orig_P}{A set of original instances whose class is oversampled with a vector of their target class appended at the last column}
	\item{K}{The value of parameter K for nearest neighbor process used for generating data}
	\item{K_all}{Unavailable for this method}
	\item{dup_size}{The maximum times of synthetic minority instances over original majority instances in the oversampling}
	\item{outcast}{Unavailable for this method}
	\item{eps}{Unavailable for this method}
	\item{method}{The name of oversampling method used for this generated dataset (SMOTE)}
}
\references{
	Chawla, N., Bowyer, K., Hall, L. and Kegelmeyer, W. 2002. SMOTE: Synthetic minority oversampling technique. Journal of Artificial Intelligence Research. 16, 321-357.
}
\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}

\examples{
    data_example = sample_generator(10000,ratio = 0.80)
	genData = SMOTE(data_example[,-3],data_example[,3])
	genData_2 = SMOTE(data_example[,-3],data_example[,3],K=7)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
