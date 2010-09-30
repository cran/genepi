\name{haplotypeOddsRatio}
\alias{haplotypeOddsRatio}
\alias{print.haploOR}
\title{
  Calculate haplotype disease risk.
}
\description{
  Haplotype disease risk is calculated resolving haplotype ambiguity and
  adjusting for covariates and population stratification.
}
\usage{
haplotypeOddsRatio(formula, gtypevar, data, stratvar=NULL, nsim=100, tol=1e-8)
\method{print}{haploOR}(x, \dots)
}
\arguments{
  \item{formula}{
    The formula for logistic regression without the haplotype variable.
  }
  \item{gtypevar}{
    The variable names in the data frame corresponding to the loci of
    interest.  Each variables counts the number of mutant genotypes a
    subject has at that locus.  Legal values are 0, 1, 2 & NA.
  }
  \item{data}{
    The name of the dataframe being analyzed.  It should have all the
    variables in the formula as well as those in genotype and stratvar.
  }
  \item{stratvar}{
    Name of the stratification variable.  This is used to account for
    population stratification.  The haplotype frequencies are estimated
    within each stratum.
  }
  \item{nsim}{
    Variance should be inflated to account for inferred ambiguous
    haplotypes.  The estimates are recalculated by simulating the
    disease haplotype copy number and variance added to average.
  }
  \item{tol}{Tolerance limit for the EM algorithm convergence.}
  \item{x}{Object of class haploOR.}
  \item{...}{Other print options.}
}
\details{
  This implements the method in the reference below.
}
\value{
  It is a list of class haploOR
  \item{call}{The function call that produced this output.}
  \item{coef}{Table with estimated coefficients, standard error,
    Z-statistic and p-value.}
  \item{var}{Covariance matrix of the estimated log odds-ratiios.}
  \item{deviance}{Average of the simulated deviances.  Its theoretical
    properties are unknown.}
  \item{aic}{Average of the simulated aic.}
  \item{null.deviance}{Deviance of null model.}
  \item{df.null}{Degrees of freedom of null model.}
  \item{df.residual}{Degrees of freedom of full model.}

  The "print" method formats the results into a user-friendly table.
}
\references{
   Venkatraman ES, Mitra N, Begg CB. (2004) A method for evaluating the
   impact of individual haplotypes on disease incidence in molecular
   epidemiology studies. \emph{Stat Appl Genet Mol Biol.} v3:Article27.
}
\author{
  Venkatraman E. Seshan
}

\examples{
# simulated data with 2 loci haplotypes 1=00, 2=01, 3=10, 4=11
# control haplotype probabilities p[i]  i=1,2,3,4
# haplotype pairs (i<=j) i=j: probs = p[i]^2 ; i<j: p[i]*p[j]
p <- c(0.25, 0.2, 0.2, 0.35)
p0 <- rep(0, 10)
l <- 0
for(i in 1:4) {for(j in i:4) {l <- l+1; p0[l] <- 2*p[i]*p[j]/(1+1*(i==j))}}
controls <- as.numeric(cut(runif(1000), cumsum(c(0,p0)), labels=1:10))
# case probabilities disease haplotype is 11
or <- c(2, 5)
p1 <- p0*c(1,1,1,2,1,1,2,1,2,8); p1 <- p1/sum(p1)
cases <- as.numeric(cut(runif(1000), cumsum(c(0,p1)), labels=1:10))
# now pool them together and set up the data frame
dat <- data.frame(status=rep(0:1, c(1000, 1000)))
# number of copies of mutant variant for locus 1
dat$gtype1 <- c(0,0,1,1,0,1,1,2,2,2)[c(controls, cases)]
# number of copies of mutant variant for locus 2
dat$gtype2 <- c(0,1,0,1,2,1,2,0,1,2)[c(controls, cases)]
# true number of copies of disease haplotype
dat$hcn <- c(0,0,0,1,0,0,1,0,1,2)[c(controls, cases)]
# model with genotypes only
haplotypeOddsRatio(status ~ 1, c("gtype1","gtype2"), dat)
# model from the logistic fit using the number of copies of disease haplotype
glm(status ~ factor(hcn), dat, family=binomial)
}
% library(genepi, lib.loc="genepi.Rcheck/")
% xpd <- read.table("~/Gdrive/SHGchapter/gem-xpd.txt", header=T, sep="\t")
% haplotypeOddsRatio(status ~ age*sex, gtypevar=c("xpd312","xpd751"), stratvar="centre", data=xpd)

\keyword{regression}
\keyword{odds-ratio}
