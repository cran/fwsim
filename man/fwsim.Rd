\name{fwsim}
\alias{fwsim}
\title{Fisher-Wright Population Simulation}
\description{This package provides tools to simulate a population under the Fisher-Wright model with a stepwise neutral mutation process on \eqn{r} loci. The population sizes are random Poisson distributed and exponential growth is supported.

Let \eqn{M} be the population size at generation \eqn{i} and \eqn{N} the population size at generation \eqn{i + 1}.
Then we assume that \eqn{N} conditionally on \eqn{M} is \eqn{\mbox{Poisson}(\alpha M)}{Poisson(\alpha*M)} distributed for \eqn{\alpha > 0} (\eqn{\alpha > 1} gives growth and \eqn{0 < \alpha < 1} gives decrease).

For each haplotype \eqn{x} occuring \eqn{m} times in the \eqn{i}'th generation, the number
of children \eqn{n} is \eqn{\mbox{Poisson}(\alpha m)}{Poisson(\alpha*m)} distributed independently of other haplotypes. 
It then follows that the sum of the number of haplotypes follows a \eqn{\mbox{Poisson}(\alpha M)}{Poisson(\alpha*M)} distribution (as just stated in the previous paragraph) and that \eqn{n} conditionally on \eqn{N} follows a \eqn{\mbox{Binomial}(N, m/M)}{Binomial(N, m/M)} as expected.

The mutation model is the stepwise neutral mutation model. There is possibility for a per locus per direction mutation rate (different up and down rates for each locus). The mutation model is implemented by dividing the number of children, following a \eqn{\mbox{Poisson}(\alpha m)}{Poisson(\alpha*m)} distribution, into categories depending on how many times they mutate. There are \eqn{r+1} categories, namely for \eqn{d = 0, 1, \ldots, r}{d = 0, 1, ..., r} mutations. The count for \eqn{d} mutations follows a \eqn{\mbox{Poisson}(\eta_d \alpha m)}{Poisson(\eta_d*\alpha*m)} distribution, where \eqn{\eta_d} is the probability for \eqn{d} mutations with \eqn{\sum_d \eta_d = 1}. See details for futher explanation.
}
\usage{
fwsim(g = 10, k = 20000, r = 7, alpha = 1, mu = 0.003, trace = TRUE, ...)
}
\arguments{
\item{g}{ number of generations to evolve. }
\item{k}{ size of the initial population. }
\item{r}{ number of loci. }
\item{alpha}{ growth rate, either a vector of size 1 or g, the latter to use a different growth rate each generation. }
\item{mu}{ mutation rate - se details. }
\item{trace}{ whether to print trace of the evolution. }
\item{\dots}{ not used. }
}
\details{
If \eqn{\mu} is a number, this is used across all loci with \eqn{\mu/2} as upwards mutation rate and \eqn{\mu/2} as downwards mutation rate for each locus.

If \eqn{\mu} is a vector of length \eqn{r}, \eqn{\mu_1, \ldots, \mu_r}, say, this is interpreted as locus specific mutation rates with equal upwards and downwards mutation rate of \eqn{\mu_i/2} for the \eqn{i}'th locus.

If \eqn{\mu} is a \eqn{2 \times r}{2 x r} matrix, then the first row in the \eqn{i}'th column is upwards mutation rate for the \eqn{i}'th locus and second row is the downwards mutation rate for the \eqn{i}'th locus.
}
\value{
A list with elements
\item{haplotypes}{ haplotypes in end population after g generations }
\item{sizes}{ the population size for each generation }
}
\author{Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen}
\examples{
pop <- fwsim(r = 2)
xtabs(N ~ Locus1 + Locus2, pop$haplotypes)

mut <- matrix(c(seq(0.001, 0.003, length.out = 13), 0.4), nrow = 2, byrow = TRUE)
mut
pop <- fwsim(g = 10, r = 7, alpha = 1.1, mu = mut, trace = TRUE)
pop
xtabs(N ~ Locus1 + Locus2, pop$haplotypes)
xtabs(N ~ Locus1 + Locus7, pop$haplotypes)
plot(pop$sizes)
}
\keyword{Fisher-Wright}

