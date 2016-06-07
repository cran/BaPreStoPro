
\name{Virkler}
\alias{Virkler}
\title{
Crack Growth
}
\description{
  68 measurement series in 164 time points
}
\usage{
Virkler
}
\format{
Sixty-eight replicate constant amplitude tests in aluminum alloy were carried out to investigate the fatigue crack propagation. In each of these tests, the number of cycles that leads to fixed crack lengths was observed. Against the natural assumption that something is observed at fixed times, here the time is the dependent variable and the crack length is the independent variable. Therefore, from the mathematical viewpoint the crack length will here be treated as time vector \eqn{t}.

The Virkler data comes as a data-frame of 164 rows and 69 columns where the first column contains the crack lengths (in mm) and the 68 following the series of observed times in load cycles up to a fixed crack length.

We want to thank Eric J. Tuegel for providing us the data that were collected by Prof. B. M. Hillberry, published in Virkler et al. (1979).
}

\examples{
data(Virkler)

Y <- t(Virkler[,-1]/10000)
t <- Virkler[,1]

plot(t, Y[1,], type = 'l', ylim = range(Y), xlab = "crack length in mm", 
   ylab = "time in load cycles / 10000")
for (i in 2:nrow(Y)){
 lines(t, Y[i,])
}

}
\source{
  \emph{Eric J. Tuegel}
}
\references{
  \emph{Virkler, D. A., Hillberry, B. M. and Goel, P. K. (1979).
 The Statistical Nature of Fatigue Crack Propagation.
Journal of Engineering Materials and Technology 101, 148--153.}
}





\keyword{datasets}
