
################################################################################
# One-sample functions
################################################################################


"pooledbinom.loglike" <-
  function(p,x,m,n=rep(1,length(m)))
  {
    if(p > 0 & p < 1)
      sum(x*log((1-(1-p)^m))) + log(1-p)*sum(m*(n-x))
    else 0
  }

"pooledbinom.loglike.vec" <-
  function(p,x,m,n=rep(1,length(m)))
  {
    np <- length(p)
    lik <- vector(length=np)
    for(i in 1:np) lik[i] <- pooledbinom.loglike(p[i],x,m,n)
    lik
  }

"score.p" <-
  function(p, x, m, n = rep(1,length(m)))
  {
    if(sum(x) == 0 & p >= 0 & p < 1) return(-sum(m*n)/(1-p))
    if(p < 0 | p >= 1) return(0)
    sum(m*x/(1-(1-p)^m) - m*n)/(1-p)
  }

"pooledbinom.mir" <-
  function(x,m,n=rep(1,length(x)))
  {
    sum(x)/sum(m*n)
  }

"pooledbinom.mle" <-
  function(x, m, n = rep(1., length(x)), tol = 1e-008)
  {
    #
    # This is the implementation using Newton-Raphson, as given
    # in the Walter, Hildreth, Beaty paper, Am. J. Epi., 1980
    #
    if(length(m) == 1.) m <- rep(m, length(x)) else if(length(m) != length(x))
      stop("\n ... x and m must have same length if length(m) > 1")
    if(any(x > n))
      stop("x elements must be <= n elements")
    if(all(x == 0.))
      return(0.)
    if(sum(x) == sum(n)) return(NA)
    p.new <- 1 - (1 - sum(x)/sum(n))^(1/mean(m)) # starting value
    done <- 0
    N <- sum(n * m)
    while(!done) {
      p.old <- p.new
      p.new <- p.old - (N - sum((m * x)/(1 - (1 - p.old)^m)))/
        sum((m^2 * x * (1 - p.old)^(m - 1))/(1 - (1 -
                                                    p.old)^m)^2)
      if(abs(p.new - p.old) < tol)
        done <- 1
    }
    p.new
  }


"pooledbinom.bias" <-
  function(p,m,n=rep(1,length(m)))
  {
    if(p > 0)
      sum((m-1)*m^2*n*(1-p)^(m-3)/(1-(1-p)^m)) *
      pooledbinom.mle.var(p,m,n)^2 / 2
    else 0
  }

# c is bias-Corrected
"pooledbinom.cmle" <-
  function(x, m, n = rep(1, length(m)), tol= 1e-8)
  {
    if(sum(x)==sum(n)) return(NA)
    phat <- pooledbinom.mle(x,m,n)
    bias <- pooledbinom.bias(phat,m,n)
    phat - bias
  }

"pooledbinom.mle.var" <-
  function(p, m, n = rep(1, length(m)))
  {
    if(p > 0 & p < 1)
      1/sum((m^2 * n * (1 - p)^(m - 2))/(1 - (1 - p)^m))
    else 0 #1/sum(m*n)
  }

"pooledbinom.I" <-
  function(p,m,n=rep(1,length(m)))
  {
    if(p > 0)
      sum((m^2 * n * (1 - p)^(m - 2))/(1 - (1 - p)^m))
    else 0 # Not sure whether to set this to 0 or not
  }

"pooledbinom.mu3" <-
  function(p, m, n = rep(1,length(m)))
  {
    if(p > 0)
      sum((m^3 * n * (1-p)^(m-3) * (2 * (1-p)^m - 1))/(1-(1-p)^m)^2)
    else 0
  }

"pooledbinom.firth" <-
  function(x, m, n = rep(1., length(x)), tol = 1e-008, rel.tol = FALSE)
  {
    #
    # This is the implementation using Newton-Raphson
    #
    if(length(m) == 1.) m <- rep(m, length(x)) else if(length(m) != length(x))
      stop("\n ... x and m must have same length if length(m) > 1")
    if(any(x > n))
      stop("x elements must be <= n elements")
    if(all(x == 0)) return(0)
    if(sum(x) == sum(n)) return(NA)

    # assign a convergence criterion function to avoid repeated
    # checks of rel.tol during iteration
    if(rel.tol) "converge.f" <- function(old,new) abs((old-new)/old)
    else "converge.f" <- function(old,new) abs(old-new)

    N <- sum(n * m)
    # use proportion of positive pools, sum(x)/sum(n)
    # and inverse of average pool size, sum(n)/N
    # ...but this doesn't work when all pools are positive, while the
    # ...MIR does, so use MIR in that case
    # 06/12/2019: nevermind--return NA for all pools positive, as above
    #if(sum(x) == sum(n)){
    #   p.new <- sum(x)/N
    #} else {
    p.new <- 1 - (1 - sum(x)/sum(n))^(sum(n)/N)
    #}

    done <- FALSE
    while(!done) {
      p.old <- p.new
      vi.p <- m^2 * n * (1-p.old)^(m-2) / (1 - (1-p.old)^m)
      wi.p <- vi.p / sum(vi.p)
      vi.p.prime <- m^2 * n * (1-p.old)^(m-3) * (2 * (1-(1-p.old)^m) - m) / (1-(1-p.old)^m)^2
      wi.p.prime <- (vi.p.prime * sum(vi.p) - vi.p * sum(vi.p.prime)) / sum(vi.p)^2
      p.new <- p.old +  (sum(m*x/(1-(1-p.old)^m) - m*wi.p/2) - (N-1/2)) /
        sum(m^2*x*(1-p.old)^(m-1) / (1-(1-p.old)^m)^2 + m * wi.p.prime / 2)
      if(p.new <= 0) p.new <- sum(x)/sum(m*n)
        if(p.new > 1) p.new <- pooledbinom.mle(x,m,n)
      if(converge.f(p.old, p.new) < tol) done <- TRUE
    }
    p.new
  }

"pooledbinom.score.ci" <-
  function(x, m, n = rep(1., length(x)), tol = .Machine$double.eps^0.75, alpha =
             0.05)
  {
    f <- function(p0, x, mm, n, alpha)
    {
      # NOTE: I use the square of the score statistic and chisq
      score.p(p0, x, mm, n)^2 * pooledbinom.mle.var(p0, mm, n) -
        qchisq(1 - alpha, 1)
    }
    # The MLE is just used for a sensible cut-value for the search ...
    # it's not needed in the computation
    p.hat <- pooledbinom.mle(x, m, n, tol)
    if(sum(x) == 0) {
      lower.limit <- 0
      z <- qnorm(1-alpha/2)
      N <- sum(m*n)
      # Upper limit should agree with unpooled case when all pools negative when using a perfect test
      # this is the Wilson score upper limit
      upper.limit <-  (z^2/2/N +  z * sqrt((z^2/4/N)/N))/(1 + z^2/N)
    }
    else {
      root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd =
                                          0., upper = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n,
                                        alpha = alpha)
      lower.limit <- uniroot(f, lower = root.brak[1], upper =
                               root.brak[2], tol = .Machine$double.eps^0.75, x = x,
                             mm = m, n = n, alpha = alpha)$root
    }
    if(sum(x) == sum(n)) {
      upper.limit <- 1
    }
    else {
      root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10,
                                        lbnd = p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x,
                                        mm = m, n = n, alpha = alpha)
      upper.limit <- uniroot(f, lower = root.brak[1], upper =
                               root.brak[2], tol = .Machine$double.eps^0.75, x = x,
                             mm = m, n = n, alpha = alpha)$root
    }
    c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
    )
  }


"pooledbinom.cscore.ci"<-
  function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
  {
    # use gamm not gam b/c of the function gam()
    f <- function(p, x, mm, n, alpha)
    {
      gamm <- function(p0, mm, n = rep(1, length(mm)))
      {
        pooledbinom.mu3(p0, mm, n) * pooledbinom.mle.var(p0, mm, n)^(3/2)
      }
      # NOTE: I use the square of the score statistic and chisq
      (score.p(p, x, mm, n) * sqrt(pooledbinom.mle.var(p, mm, n)) -
          (gamm(p, mm, n) * (qchisq(1 - alpha, 1) - 1))/6)^2 -
        qchisq(1. - alpha, 1)
    }
    # The MLE is just used for a sensible cut-value for the search ...
    # it's not needed in the computation
    p.hat <- pooledbinom.firth(x, m, n, tol) # used to be cmle, then mle
    if(sum(x) == 0)
      return(pooledbinom.score.ci(x, m, n, tol, alpha))
    # Effectively "else"
    root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
                                      = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, alpha = alpha)
    lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
      2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
      = alpha)$root
    # use mm because m is confused with maxiter
    root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd =
                                        p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
                                      alpha = alpha)
    upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
      2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
      = alpha)$root
    c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
    )
  }


"pooledbinom.bcscore.ci"<-
  function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
  {
    f <- function(p, x, mm, n, alpha)
    {
      gamm <- function(p0, mm, n = rep(1, length(mm)))
      {
        pooledbinom.mu3(p0, mm, n) * pooledbinom.mle.var(p0,
                                                         mm, n)^(3/2)
      }
      # NOTE: I use the square of the score statistic and chisq
      (score.p(p, x, mm, n) * sqrt(pooledbinom.mle.var(p, mm, n)) -
          pooledbinom.bias(p, mm, n) - (gamm(p, mm, n) * (qchisq(
            1 - alpha, 1) - 1))/6)^2 - qchisq(1. - alpha, 1)
    }
    # The MLE is just used for a sensible cut-value for the search ...
    # it's not needed in the computation
    p.hat <- pooledbinom.firth(x, m, n, tol) # used to be cmle
    p.hat <- pooledbinom.firth(x, m, n, tol) # used to be cmle
    if(sum(x) == 0)
      return(pooledbinom.score.ci(x, m, n, tol, alpha))
    # Effectively "else"
    root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
                                      = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, alpha = alpha)
    lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
      2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
      = alpha)$root
    # use mm because m is confused with maxiter
    root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd =
                                        p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
                                      alpha = alpha)
    upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
      2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
      = alpha)$root
    c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
    )
  }


"pooledbinom.lrt.ci"<-
  function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
  {
    p.hat <- pooledbinom.mle(x, m, n, tol)
    f <- function(p, x, mm, n, f.tol, alpha, p.hat)
    {
      if(p.hat == 0.)
        return( - Inf)
      else -2. * (pooledbinom.loglike(p, x, mm, n) -
                    pooledbinom.loglike(p.hat, x, mm, n)) - qchisq(
                      1. - alpha, 1.)
    }
    if(sum(x) == 0)
      return(pooledbinom.score.ci(x, m, n, tol, alpha))
    # Effectively "else"
    root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
                                      = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, f.tol = tol,
                                      alpha = alpha, p.hat = p.hat)
    lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
      2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, f.tol
      = tol, alpha = alpha, p.hat = p.hat)$root
    # use mm because m is confused with maxiter
    #print(lower.limit)
    root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd =
                                        p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
                                      f.tol = tol, alpha = alpha, p.hat = p.hat)
    upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
      2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, f.tol
      = tol, alpha = alpha, p.hat = p.hat)$root
    c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha)
  }



"pooledbinom.wald.ci"<-
  function(x, m, n = rep(1, length(x)), tol = 1e-008, alpha = 0.05)
  {
    if(length(m) == 1)
      m <- rep(m, length(x))
    p.hat <- pooledbinom.mle(x, m, n, tol)
    p.stderr <- sqrt(pooledbinom.mle.var(p.hat, m, n))
    z <- qnorm(1 - alpha/2)
    c(p = p.hat, lower = max(0, p.hat - z * p.stderr), upper = min(1, p.hat +
                                                                     z * p.stderr), alpha = alpha)
  }

"pooledbinom.mir.ci" <-
  function(x, m, n = rep(1, length(x)), tol = 1e-008, alpha = 0.05)
  {
    N <- sum(m*n)
    mir <- sum(x)/N
    mir.stderr <- sqrt(mir*(1-mir)/N)
    z <- qnorm(1-alpha/2)
    c(p=mir,lower=max(0, mir - z * mir.stderr), upper = min(1, mir + z * mir.stderr),alpha=alpha)
  }

