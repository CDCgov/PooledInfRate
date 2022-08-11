\name{vectorIndex}
\alias{vectorIndex}
\alias{vectorIndex.default}
\alias{vectorIndex.formula}
\alias{VI}
\alias{VI.default}
\alias{VI.formula}
%\alias{vectorIndexList}
%\alias{VIList}
\alias{vectorIndex.default}
\alias{VI.default}
\alias{print.vectorIndex}
\alias{print.VI}
%\alias{print.vectorIndexList}
%\alias{print.VIList}
\alias{[.vectorIndex}
\alias{[.VI}

\title{
Compute the Vector Index from pooled/group testing data
}
\description{
Calculates the Vector Index based on pooled or group testing data containing various different pool sizes and accomodating both pool-based and exogenous sources for population density computations.
}
\usage{
vectorIndex(x, \dots)

VI(x, \dots)


\method{vectorIndex}{formula}(x, data,
 n.use.traps = TRUE, n.use.na = FALSE,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)

\method{vectorIndex}{default}(x, m, n=rep(1, length(x)), vector, trap.time=rep(1,length(x)), group,
 n.use.traps = TRUE, n.use.na = FALSE,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)

\method{VI}{formula}(x, data,
 n.use.traps = TRUE, n.use.na = FALSE,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)

\method{VI}{default}(x, m, n = rep(1, length(x)), vector, trap.time = rep(1,length(x)), group,
 n.use.traps = TRUE, n.use.na = FALSE,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)

\method{print}{vectorIndex}(x, \dots)
%\method{print}{vectorIndexList}(x, \dots)
\method{as.data.frame}{vectorIndex}(x)
\method{print}{VI}(x, \dots)
%\method{print}{VIList}(x, \dots)
\method{as.data.frame}{VI}(x)

\method{[}{vectorIndex}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)
\method{[}{VI}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)

}
\arguments{
  \item{x}{If an object of class formula (or one that can be coerced to that class): a symbolic representation of the variables used to identify the number of positive pools, corresponding pool size, and number of pools of corresponding pool size, along with an optional grouping variable and required indication of the variable identifying the vector subgroups. A variable (trap.time) may be included for collection effort measured in collection time units for traps (e.g., mosquito trapping efforts as 'trap-nights') differ; the default is 1 for each record. See the 'Details' section below.  Otherwise a vector, specifying the observed number of positive pools, among the number of pools tested (\code{n}).  For methods, objects of class \code{vectorIndex}/\code{VI}.}

  \item{data}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables specified in the formula.}

  \item{n.use.traps}{logical value to indicate whether to use \emph{only} the  pool (\code{m}) and sample (\code{n}) sizes for records with missing numbers of positives (\code{x}), so including those with positivity results (\code{x}), when computing the measure of poplation size.  Default is TRUE.}

  \item{n.use.na}{logical value to indicate whether the pool (\code{m}) and sample (\code{n}) sizes for records with missing numbers of positives (\code{x}) should be included (TRUE) or not (FALSE) when computing the measure of poplation size.  Default is FALSE.}


  \item{m}{ a vector of pool sizes, must have the same length as \code{x}}

  \item{n}{ a vector of the corresponding number of pools of sizes \code{m}}

  \item{vector}{a vector indicating the "vector" to which the pool or count belongs}

  \item{trap.time}{a numeric vector recording the time-based collection effort (with standard time unit, e.g., "trap nights"); the defafult is 1 for each record}

  \item{group}{ a vector of group identifiers corresponding to the number of positive pools \code{x} and pool size \code{m}}

  \item{pt.method}{ a character string, specifying the point estimate to compute, with the following options: \code{"firth"}, bias-corrected maximum likelihood estimate (MLE)  using Firth's correction (the default);\code{"gart"} and \code{"bc-mle"}, the bias-corrected MLE; \code{"mle"}, MLE; and \code{"mir"}, minimum infection rate (MIR). }

  \item{ci.method}{ a character string, specifying the confidence interval to compute, with options: \code{"skew-score"}, skewness-corrected score interval (the default); \code{"score"} the score interval; \code{"bc-skew-score"} bias- and skewness-corrected score interval;  \code{"lrt"} likelihood ratio test interval; \code{"wald"} Wald interval; and \code{"mir"} Wald binomial interval based on the MIR.}

  \item{scale}{ a single numeric, coefficient to scale the point estimates and intervals bounds in the print and summary method (\code{\link{print.pooledBin}}, \code{\link{summary.pooledBin}}) }

  \item{alpha}{ a single numeric, specifying the type-I-error level; confidence level is 100(1-alpha)\%}

  \item{tol}{accuracy required for iterations in internal functions}

  \item{i,j}{elements to extract or replace. For [ and [[, these are numeric or character or, for [ only, empty or logical. Numeric values are coerced to integer as if by as.integer.}

  \item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left.}


  \item{\dots}{future arguments}

}


\details{

The vector index (VI) is a measure used to quantify the infectious burden in a population, such as the West Nile virus burden in populations of \emph{Culex} vector mosquito species. The language here is adopted from that for mosquito collections and testing for viral pathogens. The VI is computed as the sum over vector species (\eqn{i}{i}) of the average number of indivdiuals collected per collection effort ('trap night') \eqn{\overline{N}_i}{N_i} times an estimate of the infection prevalence or 'rate' \eqn{\widehat{p}_i}{p_i}, generally computed using pooled samples.  It is possible that not all indivdiuals are tested in a collection when computing the VI, so that the data set recording the numbers of individuals collected may be different or larger than that used in estimation of \code{p}.  To allow for a single data set to be used in computation of the VI, we permit data configuration and function parameter specifications (\code{n.use.na} and \code{n.no.traps}) to allow different data subsets to be used in computation of the components \eqn{\overline{N}_i}{N_i} and \eqn{\widehat{p}_i}{p_i}.


\code{formula} specification:  the formula is an extension of that used in the \code{\link{pooledBin}}/\code{\link{pIR}} functions; see the 'Details' section in their documentation. The extension in the formula specification to identify the variable representing the vector identifying variable is to include the variable after a 'forward slash' character.  Letting the variable \code{V} represent the vector variable, such as specification is then generically \code{X ~ m(M) + n(N) / V} when there is no grouping variable, and \code{X ~ m(M) + n(N) | G / V} when there is a grouping variable G, and multiple grouping variables are allowed, as with \code{\link{pooledBin}}/\code{\link{pIR}}.

Note that the vector variable \code{V} \emph{\bold{must}} come after the complete formula specification for the underlying \code{\link{pooledBin}}/\code{\link{pIR}} call, which includes, as needed, the grouping variable \code{G}.

A further extension to include differences in collection effort underlying the pools (e.g., the number of 'trap nights' a trap is set for mosquito trapping and all individuals in the pool are from traps set for the same amount of time), typically expressed in time of collection.  Let \code{T} represent this time variable.  Then this is included in the formula interface after a colon (:) separator \emph{after} the 'vector' variable, so extending the expression in the above paragraph:  \code{X ~ m(M) + n(N) / V:T} and \code{X ~ m(M) + n(N) | G / V:T} in the absence or presence of a grouping variable \code{G}, respectively.  By default, if such \code{T} is not specified, all collection efforts are set to the default unit 1.


Point estimate refers to underlying probability (infection rate) parameters used in the computation of the Vector Index, and they may be viewed in using a call to \code{summary}.  For point estimation of the probabilties (infection rates), the same options available for \code{pooledBin}/\code{pIR} are available here:  bias preventative ("firth") and bias corrected ("gart") estimators are recommended, with details described in Hepworth G, Biggerstaff BJ (2017). Use of MLE (\code{"mle"}), and MIR (\code{"mir"}) estimators is not recommended, but their computation is provided for historical reasons.

Note that confidence intervals are \emph{not reported} at this time for the Vector Index (even though there is an option for the method to use), as a direct population measure for inference is unclear; theoretical work on confidence intervals in this context therefore remains undeveloped.  Use of a confidence interval to accompany the Vector Index to aid characterization of uncertainty rather than for formal population inference may be justified, and inclusion of such in future releases of the PooledInfRate package is under consideration.

}

\value{
A object of class 'vectorIndex' ('VI') or, if more than one group, of class 'vectorIndexList' ('VIList').  These have a list structure with elements

\item{p}{the estimated proportion}
\item{lcl}{the lower confidence limit}
\item{ucl}{the upper confidence limit}
\item{pt.method}{the method used for point estimation }
\item{ci.method}{the method used for interval estimation }
\item{alpha}{the type-I-error level }
\item{x}{the numbers of postive pools }
\item{m}{the size of the pools }
\item{n}{the numbers of pools with corresponding pool sizes m}
\item{scale}{Scaling coefficient for the output}

along with attributes \code{class}, \code{vector}, \code{vectors.var}, \code{traptime.var}, \code{p}, \code{n}, \code{n.use.na}, \code{n.use.traps}, \code{pt.method}, \code{ci.method}, and \code{call}.
}


\references{

\emph{Walter SD, Hildreth SW, Beaty BJ:} Estimation of infection rates in population of organisms using pools of variable size. \emph{Am J Epidemiol 1980, 112(1):124-128}

\emph{Hepworth G: Estimation of proportions by group testing. PhD Dissertation.
   Melbourne, Australia: The University of Melbourne; 1999.}

\emph{Biggerstaff BJ (2008)}: Confidence interval for the difference of proportions estmimated from pooled samples. \emph{Journal of Agricultural Biological and Environmental Statistics 2008, 13(4):478-496.}

\emph{Hepworth G, Biggerstaff BJ:} Bias correction in estimating Proportions by pooled testing. \emph{JABES 2017, to appear.}
}

\author{Brad Biggerstaff}


\examples{

#######################################################################
# Consider an imaginary example for a single vector species, where pools of size
# 1, 5, 10 and 50 are tested, 5 pools of each size
# among the 5 pools with size 1 and 5, no pool is positive,
# while among the 5 pools of size 10 and 50, 1 and 2 positive
# pools are identified, respectively.
#######################################################################
# For another vector species, tested with the same design, one might find:
# a pool of size 1 individual is negative,
# 1 positive result among the pools pooling 5 elements,
# no positive result among the pools pooling 10 elements,
# 4 positive results among the pools pooling 50 elements,
#######################################################################

x1 <- c(0,0,1,2)
m1 <- c(1,5,10,50)
n1 <- c(5,5,5,5)

x2<-c(0,1,0,4)
m2 <- c(1,5,10,50)
n2 <- c(5,5,5,5)

ex.dat <- data.frame(NumPos = c(x1,x2),
                     PoolSize = c(m1,m2),
                     NumPools = c(n1,n2),
                     Species = rep(1:2, each=4),
                     TrapNights = c(2,1,1,2,2,3,2,1))

# the Vector Index is thus computed using

vectorIndex(NumPos ~ PoolSize + n(NumPools) / Species:TrapNights , data = ex.dat)

ex2.dat <- rbind(ex.dat,ex.dat)
ex2.dat$Group <- rep(LETTERS[1:2],each=8)
ex2.dat$NumPos[16] <- 0 # just to make them different

vectorIndex(NumPos ~ PoolSize + n(NumPools) |Group / Species:TrapNights , data = ex2.dat)
summary(vectorIndex(NumPos ~ PoolSize + n(NumPools) |Group / Species:TrapNights , data = ex2.dat))


}

