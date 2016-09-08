\name{ANS}
\alias{ANS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Adaptive Neighbor Synthetic Majority Oversampling TEchnique
}
\description{
	Generate a oversampling dataset from imbalanced dataset using Adaptive Neighbor SMOTE which provides the parameter K to each minority instance automatically
}
\usage{
ANS(X, target, dupSize = 0)
}

\arguments{
  \item{X}{
	A data frame or matrix of numeric-attributed dataset 
}
  \item{target}{
	A vector of a target class attribute corresponding to a dataset X.
}
  \item{dupSize}{
	A number of vector representing the desired times of synthetic minority instances over the original number of majority instances, 0 for balanced dataset.
}
  
}

\value{
	\item{data}{A resulting dataset consists of original minority instances, synthetic minority instances and original majority instances with a vector of their respective target class appended at the last column}
	\item{syn_data}{A set of synthetic minority instances with a vector of minority target class appended at the last column}
	\item{orig_N}{A set of original instances whose class is not oversampled with a vector of their target class appended at the last column}
	\item{orig_P}{A set of original instances whose class is oversampled with a vector of their target class appended at the last column}
	\item{K}{A vector of parameter K for each minority instance }
	\item{K_all}{The value of parameter C for nearest neighbor process used for identifying outcasts}
	\item{dup_size}{The maximum times of synthetic minority instances over original majority instances in the oversampling}
	\item{outcast}{A set of original minority instances which is defined as minority outcast}
	\item{eps}{The value of eps which determines automatic K}
	\item{method}{The name of oversampling method used for this generated dataset (ANS)}

	}
\references{
   Siriseriwan, W. and Sinapiromsaran, K. Adaptive neighbor Synthetic Minority Oversampling TEchnique under 1NN outcast handling.Songklanakarin Journal of Science and Technology.

}
\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}



\examples{
	data_example = sample_generator(5000,ratio = 0.80)
	genData = ANS(data_example[,-3],data_example[,3])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
