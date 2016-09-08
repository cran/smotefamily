\name{n_dup_max}
\alias{n_dup_max}

\title{
	The function to calculate the maximum round each sampling is repeated
}
\description{
	The function to calculate the maximum round each sampling is repeated, if dup_size is given as 0 then, it calculates the maximum round the number of positive instances to be duplicated to nearly match the number of negative instances 
}
\usage{
n_dup_max(size_input, size_P, size_N, dup_size = 0)
}

\arguments{
  \item{size_input}{
	The size of overall dataset
}
  \item{size_P}{
	The number of positive instances
}
  \item{size_N}{
	The number of negative instances
}
  \item{dup_size}{
	A number or vector of the number of times to be duplicated. The default is zero which means duplicating until nearly balanced.
}
}

\value{
	If dup_size is zero or contains zero, the number of rounds to duplicate positive to nearly equal to the number of negative instances
	If dup_size is not zero or contains no zero, the maximum value in dup_size
}

\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}


\examples{
	data_example = sample_generator(10000,ratio = 0.80)
	P = data_example[data_example[,3]=="p",-3]
	N = data_example[data_example[,3]=="n",-3]
	D = rbind(P,N)
	max_round =n_dup_max(nrow(D),nrow(P),nrow(N),dup_size= 0)

}

