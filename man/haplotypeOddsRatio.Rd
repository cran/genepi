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
    Variables should be as is.  Functions including I() can give error.
    E.g:  "disease ~ age*sex" is correct but "disease ~ I(age^2)" is not.
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

% \examples{
% library(genepi, lib.loc="genepi.Rcheck/")
% xpd <- read.table("~/Gdrive/SHGchapter/gem-xpd.txt", header=T, sep="\t")
% haplotypeOddsRatio(status ~ age*sex, gtypevar=c("xpd312","xpd751"), stratvar="centre", data=xpd)
% }
\keyword{regression}
\keyword{odds-ratio}
