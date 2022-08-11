\name{pooledBin}
\alias{pooledBin}
\alias{pIR}
\alias{pooledBin.default}
\alias{pooledBin.formula}
\alias{pIR.default}
\alias{pIR.formula}
\alias{print.pooledBin}
\alias{print.pIR}
%\alias{print.pooledBinList}
%\alias{print.pIRList}
\alias{plot.pooledBin}
\alias{plot.pIR}
%\alias{plot.pooledBinList}
%\alias{plot.pIRList}
\alias{[.pooledBin}
\alias{[.pIR}

\title{
Estimates for a proportion from pooled/group testing
}

\description{
Calculates point estimates and confidence intervals for a single proportion based on pooled or group testing data with equal or different pool sizes.
}

\usage{
pooledBin(x, \dots)

pIR(x, \dots)

\method{pooledBin}{formula}(x, data,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)

\method{pooledBin}{default}(x, m, n = rep(1, length(x)), group,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)

\method{pIR}{formula}(x, data,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)

\method{pIR}{default}(x, m, n = rep(1, length(x)), group,
 pt.method = c("firth", "gart", "bc-mle", "mle", "mir"),
 ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "wald", "mir"),
 scale = 1, alpha = 0.05, tol = .Machine$double.eps^0.5, \dots)


\method{print}{pooledBin}(x, \dots)
%\method{print}{pooledBinList}(x, \dots)
\method{plot}{pooledBin}(x, pch=16, refline=TRUE, printR2=TRUE,layout=NULL,\dots)
%\method{plot}{pooledBinList}(x, pch=16, refline=TRUE, printR2=TRUE,layout=NULL,\dots)
\method{as.data.frame}{pooledBin}(x)
\method{[}{pooledBin}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)

\method{print}{pIR}(x, \dots)
%\method{print}{pIRList}(x, \dots)
\method{plot}{pIR}(x, pch=16, refline=TRUE, printR2=TRUE,layout=NULL,\dots)
%\method{plot}{pIRList}(x, pch=16, refline=TRUE, printR2=TRUE, layout=NULL,\dots)
\method{as.data.frame}{pIR}(x)
\method{[}{pIR}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)
%\method{rame}{pIRList}(x, row.names = NULL, optional = FALSE, \dots)



}

\arguments{
  \item{x}{If an object of class formula (or one that can be coerced to that class): a symbolic representation of the variables used to identify the number of positive pools, corresponding pool size, and number of pools of corresponding pool size, along with an optional grouping variable.  See the 'Details' section below. Otherwise a vector, specifying the observed number of positive pools, among the number of pools tested (\code{n}). Missing data are omitted. For methods, objects of class \code{pooledBin}/\code{pIR}.}

  \item{data}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables specified in the formula. Records with mising data are omitted.}

  \item{m}{a vector of pool sizes, must have the same length as \code{x}.  Missing data are omitted.}

  \item{n}{a vector of the corresponding number of pools of sizes \code{m}. Missing data are omitted.}

  \item{group}{ a vector of group identifiers corresponding to the number of positive pools \code{x} and pool size \code{m}. Missing data are omitted.}

  \item{pt.method}{a character string, specifying the point estimate to compute, with the following options: \code{"firth"}, bias-corrected maximum likelihood estimate (MLE)  using Firth's correction (the default);\code{"gart"} and \code{"bc-mle"}, the bias-corrected MLE; \code{"mle"}, MLE; and \code{"mir"}, minimum infection rate (MIR). }

  \item{ci.method}{a character string, specifying the confidence interval to compute, with options: \code{"skew-score"}, skewness-corrected score interval (the default); \code{"score"} the score interval; \code{"bc-skew-score"} bias- and skewness-corrected score interval;  \code{"lrt"} likelihood ratio test interval; \code{"wald"} Wald interval; and \code{"mir"} Wald binomial interval based on the MIR.}

  \item{scale}{a single numeric, coefficient to scale the point estimates and intervals bounds in the print and summary method (\code{\link{print.pooledBin}}, \code{\link{summary.pooledBin}}) }

  \item{alpha}{a single numeric, specifying the type-I-error level; confidence level is 100(1-alpha)\%}

  \item{tol}{accuracy required for iterations in internal functions}

  \item{\dots}{future arguments}

  \item{pch}{lotting character}

  \item{refline}{logical valute to indicate whether to inlcude the reference line for graphic interpretation}

  \item{printR2}{include printing of the \eqn{R^2} value on the graphic}

  \item{layout}{2-element, \code{c(row, col)}, numeric vector for the layout of the plots for multiple groups.  Default is NULL, for which a layout is computed}

  \item{i,j}{elements to extract or replace. For [ and [[, these are numeric or character or, for [ only, empty or logical. Numeric values are coerced to integer as if by as.integer.}

  \item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left.}

}


\details{

The model for pooled binary data is: Assume there are independent binary (\code{0/1, negative/positive}) observations \eqn{Y_j}, \eqn{j = 1, 2, \dots, N}, positive with probability \eqn{p}, and that they are pooled/grouped into batches/pools of possbily varying sizes, \eqn{m_i}, with \eqn{n_i} pools of size \eqn{m_i}.  The observed values \code{X_i} are binary results of some test, assay, or other ascertainment of positivity of the pool.  For \eqn{i = 1, 2, \dots, d}, the number of positive pools \eqn{X_i} of size \eqn{m_i} is distributed \eqn{X_i ~ Binomial(n_i, 1-(1-p)^{m_i})}, with \eqn{n_i} pools of size \eqn{m_i}, independently.  The number of individuals (having the unobserved \code{Y_j} values) is \code{N} \eqn{= \sum_{i=1}^d m_i n_i}.


\code{formula} specification:    The basic structure of the formula interface echoes the model formula structure used in standard R model functions like \code{lm} and \code{glm}: 'number of positive pools' ~ 'pool size'.  For commonly used binary (0/1) variable \code{X} = 'number of positive pools' for pools of sizes \code{M} = 'pool size', the formula is \code{X ~ M}.  As a generalization, \code{X} = 'number of positive pools' may be a number > 1 representing the number of positive pools of \code{N} = 'number of pools of size \code{M}' (so 0  \eqn{\le} \code{X} \eqn{\le} \code{N}). The formula representation now requires identification of both the pool size variable (\code{M}) and the number of pools variable (\code{N}).  This is done using functional notation in the formula, using \code{m()} and \code{n()} to identify the variables for 'pool size' and 'number of pools', respectively, so that the basic formula is extended to \code{X ~ m(M) + n(N)}.  Because the pool size variable identified by \code{m()} is \emph{required} for use of these functions to make sense, specificaiton by \code{m()} is optional to avoid the annoyance of having to type \code{m()} for each call; examples are given below.  Note that if the 'number of pools' variable is needed, use of \code{n()} to identify this variable is required. The final extension for the formula is to indicate a grouping variable, so that estiamtes are produced for each group separately.  This is indicated in the formula using a 'conditioning' indicator '\code{|}' separating the part of the formula above from the grouping variable(s), \code{Group}. The resulting basic formula specification is \code{X ~ m(M) + n(N) | Group}, and multiple grouping variables may be specified (as of package verison 1.2) and indicated by separating them with a \code{*}, as \code{X ~ m(M) + n(N) | Group1 * Group2}.  Since the \code{m()} indication is optional, the following, like-identified forms (a--d) are equivalent formula specifications:
\tabular{rrcl}{
a) \tab \code{X} \tab \code{~} \tab \code{m(M)} \cr
a) \tab \code{X} \tab \code{~} \tab \code{M} \cr
\cr
\cr
b) \tab \code{X} \tab \code{~} \tab \code{m(M) + n(N)}\cr
b) \tab \code{X} \tab \code{~} \tab \code{n(N) + m(M))}\cr
b) \tab \code{X} \tab \code{~} \tab \code{M  + n(N)}\cr
b) \tab \code{X} \tab \code{~} \tab \code{n(N) + M}\cr
\cr
\cr
c) \tab \code{X} \tab \code{~} \tab \code{m(M) | Group}\cr
c) \tab \code{X} \tab \code{~} \tab \code{M | Group}\cr
\cr
\cr
d) \tab \code{X} \tab \code{~} \tab \code{m(M) + n(N) | Group}\cr
d) \tab \code{X} \tab \code{~} \tab \code{M  + n(N) | Group}\cr
d) \tab \code{X} \tab \code{~} \tab \code{n(N) + m(M) | Group}\cr
d) \tab \code{X} \tab \code{~} \tab \code{n(N)  + M | Group}
\cr
\cr
e) \tab \code{X} \tab \code{~} \tab \code{m(M) | Group1 * Group2}\cr
e) \tab \code{X} \tab \code{~} \tab \code{M | Group1*Group2}\cr
}
Calls with \code{formula} with other than the detailed formula structures, \emph{e.g.}, with more than two variables in the RHS, will result in error or incorrect results.  Similarly, numerical input (such as \code{c(1,0,0,0) ~ c(50,25,10,5)}) in the formula will result in an error; such input will work directly with a comma instead of a \code{~} symbol as \code{pooledBin(c(1,0,0,0), c(50,25,10,5))} via an internal call to \code{pooledBin.default}. Finally, we just note that the use of symbols 'm' and 'n' for \code{m()} and \code{n()} comes from the mathematical development in the key references.

Point estimation: bias preventative ("firth") and bias corrected ("gart") estimators are recommended, with details described in Hepworth G, Biggerstaff BJ (2017). Use of MLE (\code{"mle"}), and MIR (\code{"mir"}) estimators are not recommended, but their computation is provided for historical reasons.

Confidence intervals: All confidence intervals are based on asymptotic methods.  Score-based intervals are recommended, with the skewness-corrected score interval showing geenerally best performance and so is the default (see Hepworth (2005) or Biggerstaff (2008)).  The likelihood ratio test-based confidence interval is also included, as are the Wald interval and an interval based in the MIR.  Note that the methods "mir" and "wald" cannot be recommended, because they return too narrow intervals in relevant situations, "mir" because it ignores the pooling, and "wald" because it relies on simple large sample methods. For computational details and simulation results of the methods, see Hepworth (2005) and Biggerstaff (2008).

The subsetting or extractor functions \code{[.pooledBin} and \code{[.pIR} mimic the same \code{[} behavior as with data frames.

}


\value{
A object of class 'pooledBin' ('pIR').  These have a list structure with elements

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

and attributes \code{class}, \code{group.names}, \code{group.var}, and \code{class}.


}

\references{

Walter SD, Hildreth SW, Beaty BJ. Estimation of infection rates in population of organisms using pools of variable size. \emph{American Journal of  Epidemiology}, 112(1):124-128, 1980.

Hepworth G: Estimation of proportions by group testing. PhD Dissertation. Melbourne, Australia: The University of Melbourne; 1999.

Hepworth G. Confidence intervals for proportions estimated by group testing with groups of unequal size, \emph{Journal of Agricultural Biological and Environmental Statistics}, 10(4):478-497, 2005.

Biggerstaff BJ. Confidence interval for the difference of proportions estmimated from pooled samples. \emph{Journal of Agricultural Biological and Environmental Statistics}, 13(4):478-496, 2008.

Hepworth G, Biggerstaff BJ. Bias correction in estimating proportions by pooled testing. \emph{Journal of Agricultural Biological and Environmental Statistics}, 22(4):602-614, 2017.
}


\author{Brad Biggerstaff}


\examples{

# Consider an imaginary example, where pools of size
# 1, 5, 10 and 50 are tested, 5 pools of each size.
# Among the 5 pools with sizes 1 and 5, no pool is positive,
# while among the 5 pools of sizes 10 and 50, 1 and 2 positive
# pools are identified, respectively.

x1 <- c(0,0,1,2)
m1 <- c(1,5,10,50)
n1 <- c(5,5,5,5)

pooledBin(x=x1, m=m1, n=n1)
pooledBin(x=x1, m=m1, n=n1, scale=1000)

pooledBin(x=x1, m=m1, n=n1)

summary(pooledBin(x=x1, m=m1, n=n1), scale=1000)

### to use the formula interface, store the data in a data frame
ex.dat <- data.frame(NumPos = x1, PoolSize = m1, NumPools = n1)

pooledBin(NumPos ~ PoolSize + n(NumPools), data = ex.dat)

# without the NumPools variable, just as an example
pooledBin(NumPos ~ m(PoolSize), data = subset(ex.dat,NumPos<2))
pooledBin(NumPos ~ PoolSize, data = subset(ex.dat,NumPos<2))

#######################################################################
# For another population, tested with the same design, one might find:
# 1 positive result among the pools pooling 5 elements,
# no positive result among the pools pooling 10 elements,
# 4 positive results among the pools pooling 50 elements,
#######################################################################

x2<-c(0,1,0,4)
m2 <- c(1,5,10,50)
n2 <- c(5,5,5,5)

pooledBin(x=x2, m=m2, n=n2)

# Some other methods for the confidence bounds:

pooledBin(x=x2, m=m2, n=n2, ci.method="lrt")

### These two populations may be analyzed in one call if they are included in a single data frame

ex2.dat <- data.frame(NumPos = c(x1,x2),
                      PoolSize = c(m1,m2),
                      NumPools = c(n1,n2),
                      Population = rep(1:2, each=4))

pooledBin(NumPos ~ PoolSize + n(NumPools) | Population, data = ex2.dat)

####################################################
# Reproducing some of the estimates from Table 1 in
# Hepworth & Biggerstaff (2017):
####################################################

pooledBin(x=c(1,2), m=c(20,5), n=c(8,8), pt.method="firth", ci.method="lrt")

pooledBin(x=c(7,8), m=c(20,5), n=c(8,8), pt.method="firth", ci.method="lrt")

}

