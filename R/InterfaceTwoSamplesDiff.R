
"pooledBinDiff" <- function(x,...){
  UseMethod("pooledBinDiff")
}

"pooledBinDiff.default" <-
  function(x,m,x2,m2,n=rep(1,length(x1)),n2=rep(1,length(x2)),
           group.names = c("1","2"),
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
    # had to remove the 1s to match the generic definition
    x1 <- x
    m1 <- m
    n1 <- n

    xmn1.nomiss <- (!is.na(x1) & !is.na(m1) & !is.na(n1))
    x1 <- x1[xmn1.nomiss]
    m1 <- m1[xmn1.nomiss]
    n1 <- n1[xmn1.nomiss]

    m1n1.nozero <- (m1>0 & n1>0)
    x1 <- x1[m1n1.nozero]
    m1 <- m1[m1n1.nozero]
    n1 <- n1[m1n1.nozero]

    xmn2.nomiss <- (!is.na(x2) & !is.na(m2) & !is.na(n2))
    x2 <- x2[xmn2.nomiss]
    m2 <- m2[xmn2.nomiss]
    n2 <- n2[xmn2.nomiss]

    m2n2.nozero <- (m2>0 & n2>0)
    x2 <- x2[m2n2.nozero]
    m2 <- m2[m2n2.nozero]
    n2 <- n2[m2n2.nozero]

    call <- match.call()
    call[[1]] <- as.name("pooledBinDiff")
    pt.method <- match.arg(pt.method)
    ci.method <- match.arg(ci.method)
    if(ci.method=="mir" | pt.method=="mir"){
      ci.method <- "mir"
      pt.methid <- "mir"
    }
    if(pt.method == "gart") pt.method <- "bc-mle" # backward compatibility


    switch(pt.method,
           "firth" = {d <- pooledbinom.firth(x1,m1,n1,tol) - pooledbinom.firth(x2,m2,n2,tol)},
           "mle" = { d <- pooledbinom.mle(x1,m1,n1,tol) - pooledbinom.mle(x2,m2,n2,tol)},
           "bc-mle" ={ d <- pooledbinom.cmle(x1,m1,n1,tol) - pooledbinom.cmle(x2,m2,n2,tol)},
           "mir" = {d <- pooledbinom.mir(x1,m1,n1) - pooledbinom.mir(x2,m2,n2)}
    )
    if( (pooledbinom.cmle(x1,m1,n1,tol) < 0 | pooledbinom.cmle(x2,m2,n2,tol) < 0 ) & pt.method=="bc-mle"){
      pt.method <- "mle"
      warning("Bias-correction results in a negative point estimate; using MLE\n")
      d <- pooledbinom.mle(x1,m1,n1,tol) - pooledbinom.mle(x2,m2,n2,tol)
    }

    switch(ci.method,
           "skew-score" = { ci.d <- pooledbinom.diff.cscore.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "bc-skew-score" ={ ci.d <- pooledbinom.diff.bcscore.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "score" = {ci.d <- pooledbinom.diff.score.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "lrt" = {ci.d <- pooledbinom.diff.lrt.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "wald" = {ci.d <- pooledbinom.diff.wald.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "mir" = {ci.d <- pooledbinom.diff.mir.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]}
    )
    if(is.na(d)) ci.d <- c(NA,NA)
    ans.grp <- list(pooledBin.fit(x1,m1,n1,pt.method=pt.method,ci.method=ci.method,scale=scale,alpha=alpha,tol=tol),
                    pooledBin.fit(x2,m2,n2,pt.method=pt.method,ci.method=ci.method,scale=scale,alpha=alpha,tol=tol))

    ans.lst <- list(d=d,lcl=ci.d[1],ucl=ci.d[2],pt.method=pt.method,ci.method=ci.method,alpha=alpha,
                   scale=scale,x1=x1,m1=m1,n1=n1,x2=x2,m2=m2,n2=n2,call=call)

  ans <- data.frame(Diff = d,
                    Lower = ci.d[1],
                    Upper = ci.d[2],
                    Scale = scale)
  if(scale == 1) ans$Scale <- NULL

  structure(ans,class="pooledBinDiff",fullList = ans.lst,
            pt.method=pt.method,ci.method=ci.method,alpha=alpha,nComparisons=1,
            comparisonNames = paste0(group.names,collapse=" - "),
            grp.pooledBin = ans.grp, group.var="Group",group.names=group.names,scale=scale,call=call)

  }


"pooledBinDiff.formula" <-
  function(x, data,
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
    call <- match.call()
    call[[1]] <- as.name("pooledBinDiff")

    pt.method <- match.arg(pt.method)
    ci.method <- match.arg(ci.method)
    if(ci.method=="mir" | pt.method=="mir"){
      ci.method <- "mir"
      pt.methid <- "mir"
    }
    if(pt.method == "gart") pt.method <- "bc-mle" # backward compatability

    have.df <- (!missing(data))
    if(!have.df)
      data <- environment(x)
    vars <- pooledBinParseFormula(x, data)
    #if(any(sapply(vars,length)>1)) stop("only variable names permitted in formula; perhaps use the default call")

    # omit records with missing data -- note, if data contains records missing
    # anywhere (even not in X, M, N, Group, they are omitted), so care should be
    # used in subsetting before the call to be assured of the desired analysis
    # restrict to only those variables needed before subsetting to avoid that issue.
    if(have.df){
      # restrict the data to the variables needed before removing records with any NA
      # --this doesn't help when no data are specified, so that's dealt with below
      data <- data[,unlist(vars)]
      data <- na.omit(data)
    }

    # retrieve values from data using the character name
    # use the eval(parse(text=XX), data) construct in case data is left unspecified,
    # and because we have the names of the variables available
    # from pooledBin.deparseFormula
    x <- eval(parse(text=vars$x), data)
    m <- eval(parse(text=vars$m), data)
    if(!is.null(vars$n))
      n <- eval(parse(text=vars$n), data)
    else n <- rep(1,length(x)) # default n

    if(is.null(vars$group)) stop("Group variable must be specified for computation of group differences.")

    group <- eval(parse(text=vars$group), data)

    # restrict to data with no missing x, m, n, group
    #if(missing(data)){
      xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n) & !is.na(group))
      x <- x[xmng.nomiss]
      m <- m[xmng.nomiss]
      n <- n[xmng.nomiss]
      group <- group[xmng.nomiss]
      groups <- unique(group)
      nGroups <- length(groups)
    #}

    mn.nozero <- (m>0 & n>0)
    x <- x[mn.nozero]
    m <- m[mn.nozero]
    n <- n[mn.nozero]
    group <- group[mn.nozero]


    groups <- unique(group)
    nGroups <- length(groups)
    if(nGroups < 2) stop("Must have at least 2 groups for computation of group differences.")
    nComparisons <- nGroups * (nGroups - 1) / 2
    ans <- vector(mode="list",length=nComparisons)
    ans.grp <- vector(mode="list",length=nGroups)
    iter <- 0
    comparisonNames <- vector(length=nComparisons)
    # this is really just to have the individual results to pass on to the summary.pooledBin() function
    for(i in 1:nGroups){
      grp <- group == groups[i]
      ans.grp[[i]] <- pooledBin.fit(x[grp],m[grp],n[grp],
                                        pt.method=pt.method,ci.method=ci.method,
                                        scale=scale,alpha=alpha,tol=tol)
    }
    for(i in 1:(nGroups-1)){
      for(j in (i+1):nGroups){
        iter <- iter + 1
        g1 <- group == groups[i]
        g2 <- group == groups[j]
        x1 <- x[g1]
        m1 <- m[g1]
        x2 <- x[g2]
        m2 <- m[g2]
        n1 <- n[g1]
        n2 <- n[g2]

        comparisonNames[iter] <- paste0(as.character(groups[i])," - ",as.character(groups[j]))
        ans[[iter]] <- attr(pooledBinDiff.default(x1,m1,x2,m2,n1,n2,
                                             pt.method=pt.method,
                                             ci.method=ci.method,
                                             scale=scale,alpha=alpha,tol=tol,
                                             #comparisonNames = comparisonNames[iter])
                                             group.names = trimws(unlist(strsplit(comparisonNames[iter],"-")))),
                            "fullList")

      }
    }
    names(ans) <- comparisonNames
    if(length(ans) == 1) ans <- ans[[1]]

    ans.lst <- ans

    # annoying thing that with 2 groups the list doesn't quite work right
    # as the elements aren't accessible using $
    if(nGroups == 2){
        ans <- data.frame(Comparison = comparisonNames,
                     Diff = scale * ans.lst$d,
                     Lower = scale * ans.lst$lcl,
                     Upper = scale * ans.lst$ucl,
                     Scale = scale)

    } else {
        ans <- data.frame(Comparison = comparisonNames,
                     Diff = scale * sapply(ans.lst, function(u) u$d),
                     Lower = scale * sapply(ans.lst, function(u) u$lcl),
                     Upper = scale * sapply(ans.lst, function(u) u$ucl),
                     Scale = scale)
    }


    if(scale == 1) ans$Scale <- NULL

    structure(ans, class = "pooledBinDiff",fullList = ans.lst, pt.method=pt.method,ci.method=ci.method,alpha=alpha,
              scale=scale, nComparisons = nComparisons, comparisonNames = comparisonNames, grp.pooledBin = ans.grp, group.var = vars$group,
              group.names=groups, scale=scale, call=call)

  }



"print.pooledBinDiff" <- function(x, ...){
  print(as.data.frame(unclass(x)),...)
  invisible(x)
}


# to build the package, had to use a version of "[.data.frame" without the .Internal(copyDFattr(xx,x)) call
# this just copied x as a list along with the attributes of the data frame
#"[.pooledBinDiff" <- subsetDF #get("[.data.frame")
#"[.summary.pooledBinDiff" <- subsetDF #get("[.data.frame")

"[.pooledBinDiff" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
  mdrop <- missing(drop)
  Narg <- nargs() - !mdrop
  has.j <- !missing(j)
  if (!all(names(sys.call()) %in% c("", "drop")) &&
      !isS4(x))
    warning("named arguments other than 'drop' are discouraged")
  if (Narg < 3L) {
    if (!mdrop)
      warning("'drop' argument will be ignored")
    if (missing(i))
      return(x)
    if (is.matrix(i))
      return(as.matrix(x)[i])
    nm <- names(x)
    if (is.null(nm))
      nm <- character()
    if (!is.character(i) && anyNA(nm)) {
      names(nm) <- names(x) <- seq_along(x)
      y <- NextMethod("[")
      cols <- names(y)
      if (anyNA(cols))
        stop("undefined columns selected")
      cols <- names(y) <- nm[cols]
    }
    else {
      y <- NextMethod("[")
      cols <- names(y)
      if (!is.null(cols) && anyNA(cols))
        stop("undefined columns selected")
    }
    if (anyDuplicated(cols))
      names(y) <- make.unique(cols)
    attr(y, "row.names") <- .row_names_info(x, 0L)
    attr(y, "class") <- oldClass(x)
    return(y)
  }
  if (missing(i)) {
    if (drop && !has.j && length(x) == 1L)
      return(.subset2(x, 1L))
    nm <- names(x)
    if (is.null(nm))
      nm <- character()
    if (has.j && !is.character(j) && anyNA(nm)) {
      names(nm) <- names(x) <- seq_along(x)
      y <- .subset(x, j)
      cols <- names(y)
      if (anyNA(cols))
        stop("undefined columns selected")
      cols <- names(y) <- nm[cols]
    }
    else {
      y <- if (has.j)
        .subset(x, j)
      else x
      cols <- names(y)
      if (anyNA(cols))
        stop("undefined columns selected")
    }
    if (drop && length(y) == 1L)
      return(.subset2(y, 1L))
    if (anyDuplicated(cols))
      names(y) <- make.unique(cols)
    nrow <- .row_names_info(x, 2L)
    if (drop && !mdrop && nrow == 1L)
      return(structure(y, class = NULL, row.names = NULL))
    else {
      attr(y, "class") <- oldClass(x)
      attr(y, "row.names") <- .row_names_info(x,
                                              0L)
      return(y)
    }
  }
  xx <- x
  cols <- names(xx)
  x <- vector("list", length(x))
  #x <- .Internal(copyDFattr(xx, x))
  x <- as.list(xx)
  attributes(x) <- attributes(xx)
  oldClass(x) <- attr(x, "row.names") <- NULL
  if (has.j) {
    nm <- names(x)
    if (is.null(nm))
      nm <- character()
    if (!is.character(j) && anyNA(nm))
      names(nm) <- names(x) <- seq_along(x)
    x <- x[j]
    cols <- names(x)
    if (drop && length(x) == 1L) {
      if (is.character(i)) {
        rows <- attr(xx, "row.names")
        i <- pmatch(i, rows, duplicates.ok = TRUE)
      }
      xj <- .subset2(.subset(xx, j), 1L)
      return(if (length(dim(xj)) != 2L) xj[i] else xj[i,
                                                      , drop = FALSE])
    }
    if (anyNA(cols))
      stop("undefined columns selected")
    if (!is.null(names(nm)))
      cols <- names(x) <- nm[cols]
    nxx <- structure(seq_along(xx), names = names(xx))
    sxx <- match(nxx[j], seq_along(xx))
  }
  else sxx <- seq_along(x)
  rows <- NULL
  if (is.character(i)) {
    rows <- attr(xx, "row.names")
    i <- pmatch(i, rows, duplicates.ok = TRUE)
  }
  for (j in seq_along(x)) {
    xj <- xx[[sxx[j]]]
    x[[j]] <- if (length(dim(xj)) != 2L)
      xj[i]
    else xj[i, , drop = FALSE]
  }
  if (drop) {
    n <- length(x)
    if (n == 1L)
      return(x[[1L]])
    if (n > 1L) {
      xj <- x[[1L]]
      nrow <- if (length(dim(xj)) == 2L)
        dim(xj)[1L]
      else length(xj)
      drop <- !mdrop && nrow == 1L
    }
    else drop <- FALSE
  }
  if (!drop) {
    if (is.null(rows))
      rows <- attr(xx, "row.names")
    rows <- rows[i]
    if ((ina <- anyNA(rows)) | (dup <- anyDuplicated(rows))) {
      if (!dup && is.character(rows))
        dup <- "NA" %in% rows
      if (ina)
        rows[is.na(rows)] <- "NA"
      if (dup)
        rows <- make.unique(as.character(rows))
    }
    if (has.j && anyDuplicated(nm <- names(x)))
      names(x) <- make.unique(nm)
    if (is.null(rows))
      rows <- attr(xx, "row.names")[i]
    attr(x, "row.names") <- rows
    oldClass(x) <- oldClass(xx)
  }
  x
}

"[.summary.pooledBinDiff" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
  mdrop <- missing(drop)
  Narg <- nargs() - !mdrop
  has.j <- !missing(j)
  if (!all(names(sys.call()) %in% c("", "drop")) &&
      !isS4(x))
    warning("named arguments other than 'drop' are discouraged")
  if (Narg < 3L) {
    if (!mdrop)
      warning("'drop' argument will be ignored")
    if (missing(i))
      return(x)
    if (is.matrix(i))
      return(as.matrix(x)[i])
    nm <- names(x)
    if (is.null(nm))
      nm <- character()
    if (!is.character(i) && anyNA(nm)) {
      names(nm) <- names(x) <- seq_along(x)
      y <- NextMethod("[")
      cols <- names(y)
      if (anyNA(cols))
        stop("undefined columns selected")
      cols <- names(y) <- nm[cols]
    }
    else {
      y <- NextMethod("[")
      cols <- names(y)
      if (!is.null(cols) && anyNA(cols))
        stop("undefined columns selected")
    }
    if (anyDuplicated(cols))
      names(y) <- make.unique(cols)
    attr(y, "row.names") <- .row_names_info(x, 0L)
    attr(y, "class") <- oldClass(x)
    return(y)
  }
  if (missing(i)) {
    if (drop && !has.j && length(x) == 1L)
      return(.subset2(x, 1L))
    nm <- names(x)
    if (is.null(nm))
      nm <- character()
    if (has.j && !is.character(j) && anyNA(nm)) {
      names(nm) <- names(x) <- seq_along(x)
      y <- .subset(x, j)
      cols <- names(y)
      if (anyNA(cols))
        stop("undefined columns selected")
      cols <- names(y) <- nm[cols]
    }
    else {
      y <- if (has.j)
        .subset(x, j)
      else x
      cols <- names(y)
      if (anyNA(cols))
        stop("undefined columns selected")
    }
    if (drop && length(y) == 1L)
      return(.subset2(y, 1L))
    if (anyDuplicated(cols))
      names(y) <- make.unique(cols)
    nrow <- .row_names_info(x, 2L)
    if (drop && !mdrop && nrow == 1L)
      return(structure(y, class = NULL, row.names = NULL))
    else {
      attr(y, "class") <- oldClass(x)
      attr(y, "row.names") <- .row_names_info(x,
                                              0L)
      return(y)
    }
  }
  xx <- x
  cols <- names(xx)
  x <- vector("list", length(x))
  #x <- .Internal(copyDFattr(xx, x))
  x <- as.list(xx)
  attributes(x) <- attributes(xx)
  oldClass(x) <- attr(x, "row.names") <- NULL
  if (has.j) {
    nm <- names(x)
    if (is.null(nm))
      nm <- character()
    if (!is.character(j) && anyNA(nm))
      names(nm) <- names(x) <- seq_along(x)
    x <- x[j]
    cols <- names(x)
    if (drop && length(x) == 1L) {
      if (is.character(i)) {
        rows <- attr(xx, "row.names")
        i <- pmatch(i, rows, duplicates.ok = TRUE)
      }
      xj <- .subset2(.subset(xx, j), 1L)
      return(if (length(dim(xj)) != 2L) xj[i] else xj[i,
                                                      , drop = FALSE])
    }
    if (anyNA(cols))
      stop("undefined columns selected")
    if (!is.null(names(nm)))
      cols <- names(x) <- nm[cols]
    nxx <- structure(seq_along(xx), names = names(xx))
    sxx <- match(nxx[j], seq_along(xx))
  }
  else sxx <- seq_along(x)
  rows <- NULL
  if (is.character(i)) {
    rows <- attr(xx, "row.names")
    i <- pmatch(i, rows, duplicates.ok = TRUE)
  }
  for (j in seq_along(x)) {
    xj <- xx[[sxx[j]]]
    x[[j]] <- if (length(dim(xj)) != 2L)
      xj[i]
    else xj[i, , drop = FALSE]
  }
  if (drop) {
    n <- length(x)
    if (n == 1L)
      return(x[[1L]])
    if (n > 1L) {
      xj <- x[[1L]]
      nrow <- if (length(dim(xj)) == 2L)
        dim(xj)[1L]
      else length(xj)
      drop <- !mdrop && nrow == 1L
    }
    else drop <- FALSE
  }
  if (!drop) {
    if (is.null(rows))
      rows <- attr(xx, "row.names")
    rows <- rows[i]
    if ((ina <- anyNA(rows)) | (dup <- anyDuplicated(rows))) {
      if (!dup && is.character(rows))
        dup <- "NA" %in% rows
      if (ina)
        rows[is.na(rows)] <- "NA"
      if (dup)
        rows <- make.unique(as.character(rows))
    }
    if (has.j && anyDuplicated(nm <- names(x)))
      names(x) <- make.unique(nm)
    if (is.null(rows))
      rows <- attr(xx, "row.names")[i]
    attr(x, "row.names") <- rows
    oldClass(x) <- oldClass(xx)
  }
  x
}







"summary.pooledBinDiff" <- function(object, simple=FALSE, ...){
    x <- attr(object,"fullList")
    args <- list(...)
    if(is.null(args$digits)) digits <- 4
    else digits <- args$digits
    #print(x, scale=scale, ...)
    cat("\n")
    switch(attr(object,"pt.method"),
           "firth" = PtEstName <- "Firth's Correction",
           "gart" = PtEstName <- "Gart's Correction",
           "bc-mle" = PtEstName <- "Gart's Correction",
           "mle" = PtEstName <- "Maximum Likelihood",
           "mir" = PtEstName <- "Minimum Infection Rate"
    )
    switch(attr(object,"ci.method"),
           "skew-score" = CIEstName <- "Skew-Corrected Score (Gart)",
           "bc-skew-score" = CIEstName <- "Bias- & Skew-Corrected Score (Gart)",
           "score" = CIEstName <- "Score",
           "lrt" = CIEstName <- "Likelihood Ratio Test Inversion",
           "wald" = CIEstName <- "Wald",
           "mir" = CIEstName <- "Minimum Infection Rate"
    )
    grp.pooledBin <- attr(object, "grp.pooledBin")
    grp.pb <- t(sapply(grp.pooledBin,
                      function(x) c(x$scale*x$p, x$scale*x$lcl, x$scale*x$ucl, x$scale, sum(x$n * x$m), sum(x$n), sum(x$x))))


    grp.pb <- as.data.frame(grp.pb)
    names(grp.pb) <- c("P","Lower","Upper","Scale","Individuals","Pools","Positive Pools")
    grp.pb$Group <- attr(object,"group.names")
    grp.pb <- grp.pb[,c(8,1:7)]
    names(grp.pb)[1] <- attr(object,"group.var")
    rownames(grp.pb) <- 1:nrow(grp.pb)



    if(simple) return(grp.pb)

    #cat("\nCall: ", deparse(x$call), "\n\n")
    structure(object, class="summary.pooledBinDiff", df=object, grp.pooledBin = grp.pb, PtEstName = PtEstName, CIEstName = CIEstName,
              scale=scale, call = attr(object,"call"))
  }



"print.summary.pooledBinDiff" <-
  function(x, ...){
    args <- list(...)
    if(is.null(attr(x,"scale"))) scale <- 1
    else scale <- attr(x,"scale")
    if(is.null(args$digits)) digits <- 4
    else digits <- args$digits
    cat("Estimation of Difference of Binomial Proportions for Pooled Data\n\n")
    cat(paste0("Call: ", deparse(attr(x,"call"),width.cutoff = 120),"\n\n"))
    cat(paste("Point estimator:",attr(x,"PtEstName"),"\n"))
    cat(paste("CI method:",attr(x,"CIEstName"),"\n\n"))
    nComparisons <- attr(x, "nComparisons")
    if(is.null(nComparisons)) nComparisons <- 1
    comparisonNames <- attr(x,"comparisonNames")

    print(attr(x,"df"))

    cat("\n")
    cat("Group summaries:\n\n")
    print(attr(x,"grp.pooledBin"))
    invisible(x)
  }

"as.data.frame.pooledBinDiff" <- function(x){
  as.data.frame(unclass(x))
}

