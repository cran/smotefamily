\name{Borderline-SMOTE}
\alias{BLSMOTE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Borderline-SMOTE
}
\description{
	Generate synthetic positive instances using Borderline-SMOTE algorithm. The number of majority neighbor of each minority instance is used to divide minority instances into 3 groups; SAFE/DANGER/NOISE, only the DANGER are used to generate synthetic instances. 
}
\usage{
BLSMOTE(X,target,K=5,C=5,dupSize=0,method =c("type1","type2"))
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
	The number or vector representing the desired times of synthetic minority instances over the original number of majority instances, 0 for duplicating until balanced
}
  \item{method}{
	A parameter to indicate which type of Borderline-SMOTE presented in the paper is used
}

}

\value{
	\item{data}{A resulting dataset consists of original minority instances, synthetic minority instances and original majority instances with a vector of their respective target class appended at the last column}
	\item{syn_data}{A set of synthetic minority instances with a vector of minority target class appended at the last column}
	\item{orig_N}{A set of original instances whose class is not oversampled with a vector of their target class appended at the last column}
	\item{orig_P}{A set of original instances whose class is oversampled with a vector of their target class appended at the last column}
	\item{K}{The value of parameter K for nearest neighbor process used for generating data}
	\item{K_all}{The value of parameter C for nearest neighbor process used for determining SAFE/DANGER/NOISE}
	\item{dup_size}{The maximum times of synthetic minority instances over original majority instances in the oversampling}
	\item{outcast}{Unavailable for this method}
	\item{eps}{Unavailable for this method}
	\item{method}{The name of oversampling method and type used for this generated dataset (BLSMOTE type1/2)}
}
\references{
	Han, H., Wang, W.Y. and Mao, B.H.  Borderline-SMOTE: a new over-sampling method in imbalanced data sets learning. In Proceedings of the 2005 international conference on Advances in Intelligent Computing - Volume Part I (ICIC'05), De-Shuang Huang, Xiao-Ping Zhang, and Guang-Bin Huang (Eds.), Vol. Part I. Springer-Verlag, Berlin, Heidelberg, 2005. 878-887.  DOI=http://dx.doi.org/10.1007/11538059_91
}
\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}

\examples{
    data_example = sample_generator(5000,ratio = 0.80)
	genData = BLSMOTE(data_example[,-3],data_example[,3])
	genData_2 = BLSMOTE(data_example[,-3],data_example[,3],K=7, C=5, method = "type2")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
