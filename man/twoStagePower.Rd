\name{twoStagePower}
\alias{twoStagePower}
\title{
  Calculate the power of a two stage design for GWAS
}
\description{
  Calculate the power of a two stage design for GWAS under sample size
  or cost constraints.  Implements methods in the refereces below.
}
\usage{
twoStagePower(n=NULL, Cost=NULL, m=5000, mu=0.045, mu.loc=0.5, p=0.10,
            f=NULL, relcost=1, true.needed=1, rho=0, rho0=0, nsim=2000)
}
\arguments{
  \item{n}{The maximum sample size for the study.}
  \item{Cost}{
    Maximum available resource.  One of \code{Cost} or \code{n} must be
    specified.
  }
  \item{m}{
    The number of markers in the study. Default is 5000.  It will take a
    a long time to compute power for very large numbers e.g. 100000
  }
  \item{mu}{
    The mean vector for the markers that are associated with endpoint.
  }
  \item{mu.loc}{
    The locations of the true markers.  Since the chromosome is mapped
    to the unit interval (0,1) the numbers should be between 0 and 1.
  }
  \item{p}{
    The proportion of markers taken to the second stage.  The default is
    0.1 which is found to be optimal.
  } 
  \item{f}{
    The fraction of Cost or sample size allocated to the first stage.
    If not specified it uses 0.75 for the Cost constraint scenario and
    0.5 for the sample size contraint scenario.
  }
  \item{relcost}{
    Specifies how expensive it is to genotype in the second stage
    compared to the first stage.
  }
  \item{true.needed}{
    The number of markers selected in the end.  Can be a maximum of
    length of mu.loc (or mu).
  }
  \item{rho, rho0}{correlation between markers}
  \item{nsim}{Number of Monte Carlo replications to compute power.}
}
\details{
  This implements the method in the reference below.
}
\value{
  It returns the power as a single numeric value
}
\references{
  Satagopan JM, Venkatraman ES, Begg CB. (2004) Two-stage designs for
  gene-disease association studies with sample size constraints.
  \emph{Biometrics} 
}
\author{
  Jaya M. Satagopan & Venkatraman E. Seshan
}

\examples{
twoStagePower(n=1000)
twoStagePower(Cost=1000)
% library(genepi, lib.loc="genepi.Rcheck/")
% xpd <- read.table("~/Gdrive/SHGchapter/gem-xpd.txt", header=T, sep="\t")
% haplotypeOddsRatio(status ~ age*sex, gtypevar=c("xpd312","xpd751"), stratvar="centre", data=xpd)
}
\keyword{design}

