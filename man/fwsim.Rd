\name{fwsim}
\alias{fwsim}
\title{Fisher-Wright Population Simulation}
\description{Simulates Fisher-Wright populations with growth.

Let M be the population size at generation i and N the population size at generation i + 1.
Then the model is that N is Poisson(alpha*M) distributed.

For each type x occuring m times in the i'th generation, the number
of children n is Poisson(alpha*m) distributed. 
It then follows that the sum of types follows a Poisson(alpha*N) distribution 
and that n | N follows a Binomial(N, m/M) as expected.

The mutation model is the stepwise mutation model with locuswise mutation of mu.
It is implemented by dividing the number of children, following Poisson(alpha*m), 
into Poisson(dbinom(x = d, size = r, prob = mu)*alpha*m) for d = 0, 1, ..., r
number of simultaneous mutations.
}
\usage{
fwsim(g = 10, k = 20000, r = 7, alpha = 1, mu = 0.003, mutmodel = 1, trace = TRUE, trace.loc.mut = 4, ...)
}
\arguments{
\item{g}{ number of generations to evolve. }
\item{k}{ size of the initial population. }
\item{r}{ number of loci. }
\item{alpha}{ growth rate. }
\item{mu}{ mutation rate (per locus). }
\item{mutmodel}{ not used (for future use). }
\item{trace}{ whether to print trace of the evolution. }
\item{trace.loc.mut}{ threshold for number of simultaneous mutations before printed. }
\item{\dots}{ not used. }
}
\value{
A list with elements
\item{haplotypes}{ haplotypes in end population after g generations }
\item{sizes}{ the population size for each generation }
}
\author{Mikkel Meyer Andersen <mikl@math.aau.dk>}
\examples{
pop <- fwsim(r = 2)
xtabs(N ~ Locus1 + Locus2, pop$haplotypes)
}
\keyword{Fisher-Wright}

