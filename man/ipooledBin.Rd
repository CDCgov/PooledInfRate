\name{ipooledBin}
\alias{ipooledBin}
\alias{ipooledBin.default}
\alias{ipooledBin.formula}
\alias{print.ipooledBin}
%\alias{print.ipooledBinList}
\alias{[.ipooledBin}

\title{
Estimates for a proportion from pooled/group testing using an imperfect test
}

\description{
Calculates point estimates for a single proportion based on pooled or group testing data with equal or different pool sizes using an imperfect test.
}

\usage{
ipooledBin(x, \dots)

\method{ipooledBin}{formula}(x, data, pt.method = c("firth","mle"),
sens = 1, spec = 1,
scale = 1, tol=.Machine$double.eps^0.5, max.iter=10000, p.start=NULL, \dots)

\method{ipooledBin}{default}(x,m,n=rep(1,length(x)), group, pt.method = c("firth","mle"),
sens=rep(1,length(x)), spec=rep(1,length(x)), scale=1, tol=.Machine$double.eps^0.5,
max.iter=10000, p.start=NULL, \dots)

\method{print}{ipooledBin}(x, \dots)
\method{as.data.frame}{ipooledBin}(x)
\method{[}{ipooledBin}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)

}

\arguments{
  \item{x}{If an object of class formula (or one that can be coerced to that class): a symbolic representation of the variables used to identify the number of positive pools, corresponding pool size, and number of pools of corresponding pool size, along with an optoinal grouping variable.  See the 'Details' section below. Otherwise a vector, specifying the observed number of positive pools, among the number of pools tested (\code{n}). Missing data are omitted For methods, objects of class \code{ipooledBin}.}

  \item{data}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables specified in the formula. Records with mising data are omitted.}

  \item{m}{a vector of pool sizes, must have the same length as \code{x}.  Missing data are omitted.}

  \item{n}{a vector of the corresponding number of pools of sizes \code{m}. Missing data are omitted.}

  \item{group}{ a vector of group identifiers corresponding to the number of positive pools \code{x} and pool size \code{m}. Missing data are omitted.}

  \item{sens}{a vector of test sensitivities, either a single value applied to each pool, or a vector the length of x to reflect the test sensitivity for the corresponding pool.  By default set to 1.}

  \item{spec}{a vector of test specificities, either a single value applied to each pool, or a vector the length of x to reflect the test specificity for the corresponding pool.  By default set to 1.}

  \item{pt.method}{a character string, specifying the point estimate to compute, with the following options: \code{"firth"}, \code{"mle"}, MLE. }

  \item{scale}{a single numeric, coefficient to scale the point estimates and intervals bounds in the print and summary method (\code{\link{print.pooledBin}}, \code{\link{summary.pooledBin}}) }

  \item{max.iter}{maximum number of iterations for the Newton-Raphson algorithm before halting execution}

  \item{p.start}{starting value for the iterations for point estimation using Firth's correction.  A default is computed of NULL.}

  \item{tol}{accuracy required for iterations in internal functions}


  \item{i,j}{elements to extract or replace. For [ and [[, these are numeric or character or, for [ only, empty or logical. Numeric values are coerced to integer as if by as.integer.}

  \item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left.}

  \item{\dots}{future arguments}

}


\details{

The model for pooled binary data with a (possibly) imperfect test with sensitivity \eqn{a_i} and specificity \eqn{b_i} is (notation follows Hepworth and Biggerstaff (2021)): Assume there are independent binary (\code{0/1, negative/positive}) observations \eqn{Y_j}, \eqn{j = 1, 2, \dots, N}, positive with probability \eqn{p}, and that they are pooled/grouped into batches/pools of possbily varying sizes, \eqn{m_i}, with \eqn{n_i} pools of size \eqn{m_i}.  The observed values \code{X_i} are binary results of some test, assay, or other ascertainment of positivity of the pool.  For \eqn{i = 1, 2, \dots, d}, the number of positive pools \eqn{X_i} of size \eqn{m_i} is distributed \eqn{X_i ~ Binomial(n_i, a_i - r_i(1-p)^{m_i})}, with \eqn{r = a_i + b_i - 1} and with \eqn{n_i} pools of size \eqn{m_i}, independently.  The number of individuals (having the unobserved \code{Y_j} values) is \code{N} \eqn{= \sum_{i=1}^d m_i n_i}.


\code{formula} specification:    The specification of the formula interface matches \code{\link{pooledBin}}.  The basic structure of the formula interface echoes the model formula structure used in standard R model functions like \code{lm} and \code{glm}: 'number of positive pools' ~ 'pool size'.  For commonly used binary (0/1) variable \code{X} = 'number of positive pools' for pools of sizes \code{M} = 'pool size', the formula is \code{X ~ M}.  As a generalization, \code{X} = 'number of positive pools' may be a number > 1 representing the number of positive pools of \code{N} = 'number of pools of size \code{M}' (so 0  \eqn{\le} \code{X} \eqn{\le} \code{N}). The formula representation now requires identification of both the pool size variable (\code{M}) and the number of pools variable (\code{N}).  This is done using functional notation in the formula, using \code{m()} and \code{n()} to identify the variables for 'pool size' and 'number of pools', respectively, so that the basic formula is extended to \code{X ~ m(M) + n(N)}.  Because the pool size variable identified by \code{m()} is \emph{required} for use of these functions to make sense, specificaiton by \code{m()} is optional to avoid the annoyance of having to type \code{m()} for each call; examples are given below.  Note that if the 'number of pools' variable is needed, use of \code{n()} to identify this variable is required. The final extension for the formula is to indicate a grouping variable, so that estiamtes are produced for each group separately.  This is indicated in the formula using a 'conditioning' indicator '\code{|}' separating the part of the formula above from the grouping variable, say \code{Group}. The resulting formula specification is \code{X ~ m(M) + n(N) | Group}, and multiple grouping variables may be specified using a \code{*}, as \code{X ~ m(M) + n(N) | Group1 * Group2}.  Since the \code{m()} indication is optional, the following, like-identified forms (a--d) are equivalent formula specifications:
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

Point estimation: the bias preventative ("firth") estimate is recommended, with details described in Hepworth G, Biggerstaff BJ (2017). Use of MLE (\code{"mle"}) estimate is not recommended, but its computation is provided for users who need it.

Note that not all sensitivity (\code{sens}) and specificity (\code{spec}) values are compatible with all data (\code{x}, \code{m}, \code{n}) configiurations.  If convergence issues arise, consider revising \code{sens} or \code{spec}; see Hepworth and Biggerstaff (2021).

No confidence intervals are computed or reported for imperfect tests, as clear recommendations await further research.

The subsetting or extractor functions \code{[.pooledBin} and \code{[.pIR} mimic the same \code{[} behavior as with data frames.

}

\value{
A object of class 'ipooledBin' or, if more than one group, of class 'ipooledBinList'.  These have a list structure with elements

\item{p}{the estimated proportion(s)}

with attributes \code{class}, \code{group.names} (List), \code{group.var} (List), \code{x}, \code{m}, \code{n}, \code{sens}, \code{spec}, \code{scale}, \code{pt.method}, and \code{call}.

}

\references{

Hepworth G, Biggerstaff BJ. Bias correction in estimating proportions by pooled testing. \emph{Journal of Agricultural Biological and Environmental Statistics}, 22(4):602-614, 2017.

Hepworth G, Biggerstaff BJ. Bias correction in estimating proportions by imperfect pooled testing. \emph{Journal of Agricultural Biological and Environmental Statistics}, 26(1):90-104, 2021.

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

ipooledBin(x=x1, m=m1, n=n1, sens=0.95, spec=0.99)
ipooledBin(x=x1, m=m1, n=n1, sens=0.95, spec=0.99, scale = 1000)

summary(ipooledBin(x=x1, m=m1, n=n1, sens=0.95, spec=0.99, scale=1000))

### to use the formula interface, store the data in a data frame
ex.dat <- data.frame(NumPos = x1, PoolSize = m1, NumPools = n1)

ipooledBin(NumPos ~ PoolSize + n(NumPools), data = ex.dat, sens=0.95, spec=0.99)

# without the NumPools variable, just as an example
ipooledBin(NumPos ~ m(PoolSize), data = subset(ex.dat,NumPos<2), sens=0.95, spec=0.99)
summary(ipooledBin(NumPos ~ PoolSize, data = subset(ex.dat,NumPos<2), sens=0.95, spec=0.99))

}

