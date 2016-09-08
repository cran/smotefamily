\name{gap}
\alias{gap}

\title{
    The function to provide a random number which is used as a location of synthetic instance
}
\description{
  The function to provide a random number which uses to identify the location of each synthetic instance. The interval of possible values depends from safe-level values of instances in a pair.
}
\usage{
gap(sl_p = 1, sl_n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sl_p}{
     The safe-level value of the first instance
}
  \item{sl_n}{
    The safe-level value of the second instance
}
}

\value{
  A value between 0 to 1 which is used to identify the location of synthetic instance
    If sl_p >= sl_n, it gives the random number between 0 to sl_n/sl_p
	If sl_p < sl_n, it gives the random number between 1-sl_p/sl_n to 1
}

\author{
	Wacharasak Siriseriwan <wacharasak.s@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
  
	r_num = gap()
	r_num_2 = gap(sl_p = 4, sl_n = 2)	
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
