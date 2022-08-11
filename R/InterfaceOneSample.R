
"pooledBin" <- function(x,...){
  UseMethod("pooledBin")
}

"pooledBin.default" <- function(x,m,n=rep(1,length(x)), group,
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
    call <- match.call()
    call[[1]] <- as.name("pooledBin")

    xmn.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n))
    x <- x[xmn.nomiss]
    m <- m[xmn.nomiss]
    n <- n[xmn.nomiss]

    mn.nozero <- (m>0 & n>0)
    x <- x[mn.nozero]
    m <- m[mn.nozero]
    n <- n[mn.nozero]

    if(!missing(group)){
      group.var <- "Group" #deparse(substitute(group))
      xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n) & !is.na(group))
      group <- group[xmng.nomiss]
      groups <- unique(group)
      nGroups <- length(groups)
      group.dat <- data.frame(Group=group)
    } else {
      xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n))
      group <- rep("SINGLE",length(x[xmng.nomiss]))
      groups <- unique(group)
      nGroups <- length(groups)
      group.var <- ""
      group.dat <- data.frame(Group=rep(1,length(x[xmng.nomiss])))
    }
      ans <- vector(mode="list",length=nGroups)

      for(i in 1:nGroups){
        ans[[i]] <- pooledBin.fit(x[group==groups[i]],
                                                    m[group==groups[i]],
                                                    n[group==groups[i]],
                                                    pt.method=pt.method,
                                                    ci.method=ci.method,
                                                    scale=scale,alpha=alpha,tol=tol)
      }

      ans.lst <- ans

      ans <- cbind(group.dat[!duplicated(group),],
                   data.frame(
                        P = rep(0,  nGroups),
                        Lower = rep(0, nGroups),
                        Upper = rep(0, nGroups),
                        Scale = rep(1,  nGroups))
              )
      #ans <- data.frame(Group = groups,
      #                  P = rep(0,  nGroups),
      #                  Lower = rep(0, nGroups),
      #                  Upper = rep(0, nGroups),
      #                  Scale = rep(1,  nGroups))


      #print(group.var)
      #print(ans)
      #print(names(ans))
      names(ans)[1] <- ifelse(group.var != "", group.var, "Group")
      #if(group.var != "") names(ans)[1] <- group.var

      for (i in 1:nGroups) ans[i, 2:5] <- ans.lst[[i]]$scale * c(ans.lst[[i]]$p,  ans.lst[[i]]$lcl, ans.lst[[i]]$ucl, 1)
      if (all(ans$Scale == 1)) ans$Scale <- NULL

      if(nGroups == 1) ans$Group <- NULL

      # fix row names after subsetting
      group.dat <- as.data.frame(group.dat[!duplicated(group),])
      rownames(group.dat) <- 1:nrow(group.dat)

      ans <- structure(ans, class = "pooledBin",fullList = ans.lst,
                       group.names = groups, group.var = group.var,n.groups = nGroups,
                       group.dat = group.dat,
                       scale=scale,call=call)

    ans
  }


"pooledBin.formula" <- function(x, data,
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...){
    call <- match.call()
    call[[1]] <- as.name("pooledBin")

    if(missing(data))
      data <- environment(x)

    vars <- pooledBinParseFormula(x, data)
    if(any(sapply(vars[1:3],length)>1)) stop("only variable names permitted in formula; perhaps use the default call")
    # omit records with missing data -- note, if data contains records missing
    # anywhere (even not in X, M, N, Group, they are omitted), so care should be
    # used in subsetting before the call to be assured of the desired analysis
    # restrict to only those variables needed before subsetting to avoid that issue.
    if(!missing(data)){
      # restrict the data to the variables needed before removing records with any NA
      # --this doesn't help when no data are specified, so that's dealt with below
      data <- data[,unlist(vars)]
      data <- na.omit(data)
    }

    # retrieve values from data using the character name
    # use the eval(parse(text=XX), data) construct in case data is left unspecified,
    # and because we have the names of the variables available
    # from pooledBinParseFormula
    x <- eval(parse(text=vars$x), data)
    m <- eval(parse(text=vars$m), data)
    if(!is.null(vars$n))
      n <- eval(parse(text=vars$n), data)
    else n <- rep(1,length(x)) # default n

    if(!is.null(vars$group[1])){
      # NEW as of 5/24/2022 to allow multiple grouping variables
      n.group.var <- length(vars$group)
      group.dat <- as.data.frame(matrix(nrow=length(x),ncol=n.group.var))
      names(group.dat) <- vars$group

      for(i in 1:n.group.var) group.dat[,i] <- eval(parse(text=vars$group[i]),data)

      #print(interaction(group.dat,drop=TRUE)) # drop unused levels
      #return(group.dat)

      # end of NEW

      #group <- eval(parse(text=vars$group), data)
      group <- interaction(group.dat,drop=TRUE)
      groups <- unique(group)
      nGroups <- length(groups)
      group.var <- vars$group
    }  else {
      group <- rep("SINGLE",length(x))
      groups <- unique(group)
      nGroups <- 1
      group.var <- ""
      group.dat <- data.frame()
    }

      # restrict to data with no missing x, m, n, group
    if(missing(data)){
      mn.nozero <- (m>0 & n>0)
      x <- x[mn.nozero]
      m <- m[mn.nozero]
      n <- n[mn.nozero]

      if(nGroups == 1){
        xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n))
      } else {
        xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n) & !is.na(group))
      }
      x <- x[xmng.nomiss]
      m <- m[xmng.nomiss]
      n <- n[xmng.nomiss]

      group <- group[xmng.nomiss]
      groups <- unique(group)
      nGroups <- length(groups)
    }

    mn.nozero <- (m>0 & n>0)
    x <- x[mn.nozero]
    m <- m[mn.nozero]
    n <- n[mn.nozero]
    group <- group[mn.nozero]
    groups <- unique(group)
    nGroups <- length(groups)

    #if(nGroups > 1){
      ans <- vector(mode="list",length=nGroups)
      for(i in 1:nGroups){
        ans[[i]] <- pooledBin.fit(x[group==groups[i]],
                                  m[group==groups[i]],
                                  n[group==groups[i]],
                                  pt.method=pt.method,
                                  ci.method=ci.method,
                                  scale=scale,alpha=alpha,tol=tol)
      }
      names(ans) <- groups
      ans.lst <- ans

      # ans <- data.frame(Group = groups,
      #                   P = rep(0,  nGroups),
      #                   Lower = rep(0, nGroups),
      #                   Upper = rep(0, nGroups),
      #                   Scale = rep(1,  nGroups))
      ans <- cbind(group.dat[!duplicated(group),],
                   data.frame(
                        P = rep(0,  nGroups),
                        Lower = rep(0, nGroups),
                        Upper = rep(0, nGroups),
                        Scale = rep(1,  nGroups))
      )
      if (!is.null(vars$group)){
        #names(ans)[1] <- paste0(vars$group,collapse=".") # new as of 5/26/2022 for multiple grouping variables
        names(ans)[1:length(vars$group)] <- vars$group
      }
      for (i in 1:nGroups) ans[i,(ncol(group.dat)-1)+ (2:5)] <- ans.lst[[i]]$scale * c(ans.lst[[i]]$p,  ans.lst[[i]]$lcl, ans.lst[[i]]$ucl, 1)
      if (all(ans$Scale == 1)) ans$Scale <- NULL

      if(nGroups == 1) ans$Group <- NULL

      # fix row names after subsetting
      group.dat <- as.data.frame(group.dat[!duplicated(group),])
      rownames(group.dat) <- 1:nrow(group.dat)

      ans <- structure(ans, class = "pooledBin", fullList = ans.lst,
                       x.var = vars$x, m.var = vars$m, n.var = vars$n,
                       group.names = groups, group.var = group.var,
                       group.dat = group.dat,
                       n.groups = nGroups, scale=scale,call=call)
    ans

  }

"pooledBinParseFormula" <- function(formula,data=NULL){
  # this function reads a formula of the form
  # x ~ m(m.var) + n(n.var) | group
  # for a call to pooledBin and returns
  # a list with elements (x, m, n, group), which are character strings
  # with the respective names for those variables;
  # missing n() or group variables are set to NULL;
  # if a single variable is used it is assumed to be for m(), in which case
  # the m() construct may be excluded;
  # a grouping variable is optional;
  # errors for ambiguous calls are issued;

  # first determine if there's a grouping variable, indicated by "|" in the formula
  # note: grepl returns logicals for inclusion, whereas grep does not
  has.group <- any(grepl("\\|",formula))

  if(has.group){
    tmp.tms <- terms(formula,specials=c("m","n"),data=data)
    if(attr(tmp.tms,"response")==0) stop("must have a response variable representing the number of positives in the pool")
    group.var<- as.character(formula[[3]])[3]
    lhs <- as.character(formula)[2]
    rhs <- as.character(formula[[3]])[2]
    formula <- as.formula(paste0(lhs," ~ ",rhs))
    # NEW as of 5/24/2022 to allow multiple grouping variables
    group.var <- unlist(strsplit(group.var," * ",fixed=TRUE))
    # end of NEW
  } else group.var <- NULL

  tms <- terms(formula, specials = c("m","n"),data=data)
  spec <- attr(tms,"specials")
  if(attr(tms,"response")==0) stop("must have a response variable representing the number of positives in the pool")
  if(length(attr(tms,"variables")) > 4) stop("at most pool size [m()] and number of pools [n()] accepted in RHS of formula...see help page for specification")
  x.var <- as.character(attr(tms,"variables")[[2]])
  m.ind <- spec$m
  if(!is.null(m.ind)){
    m.var <- as.character(attr(tms,"variables")[[m.ind+1]][[2]])
  } else m.var <- NULL
  n.ind <- spec$n
  # ifelse didn't want to play nicely with a NULL return
  if(!is.null(n.ind)){
    n.var <- as.character(attr(tms,"variables")[[n.ind+1]][[2]])
  } else n.var <- NULL

  if(is.null(m.var) & !is.null(n.var)){
    #stop("formula must contain an m() term if it contains an n() term\n")
    m.ind <- ifelse(n.ind==2,3,2)
    m.var <- as.character(attr(tms,"variables")[[m.ind+1]])
  }

  if(is.null(m.var) & is.null(n.var)) {
    if(length(attr(tms,"variables")) > 3) stop("if no m() or n() term specified, only one variable is permitted -- interpreted as the m() term")
    m.var <- as.character(attr(tms,"variables")[[3]])
  }

  list(x = x.var, m = m.var, n = n.var, group = group.var)
}


"pooledBin.fit" <- function(x,m,n=rep(1,length(x)),
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5,...)
  {
  # check for x > n
  if(any(x>n)) stop("must have x <= n for each element")
  # if all m==1, there is no pooling, so just use the regular, Wilson score interval
  if(all(m==1)){
    z <- qnorm(1-alpha/2)
    p <- sum(x)/sum(n)
    # the usual score interval, when all pools are size 1
    cl <- (p + z^2/2/sum(n) + c(-1, 1) * z * sqrt((p * (1 -  p) + z^2/4/sum(n))/sum(n)))/(1 + z^2/sum(n))

    return(list(p=p,lcl=cl[1],ucl=cl[2],pt.method="mle",ci.method="score",alpha=alpha,x=x,m=m,n=n,scale=scale))
  }
  # check for all x == n
  if(all(x == n)){
    #warning("estimation from pooled samples with all pools positive not available--estimates set to NA")
    return(list(p=NA,lcl=NA,ucl=NA,pt.method=pt.method,ci.method=ci.method,alpha=alpha,x=x,m=m,n=n,scale=scale))
  }

    pt.method <- match.arg(pt.method)
    ci.method <- match.arg(ci.method)
    if(ci.method=="mir" | pt.method=="mir"){
      ci.method <- "mir"
      pt.method <- "mir"
    }
    if(pt.method == "gart") pt.method <- "bc-mle" # backward compatability
    switch(pt.method,
           "firth" = {p <- pooledbinom.firth(x,m,n,tol)},
           "mle" = { p <- pooledbinom.mle(x,m,n,tol)},
           "bc-mle" ={ p <- pooledbinom.cmle(x,m,n,tol)},
           "mir" = {p <- pooledbinom.mir(x,m,n)}
    )
    if(p < 0 & pt.method=="bc-mle"){
      pt.method <- "mle"
      warning("Bias-correction results in negative point estimate; using MLE\n")
      p <- pooledbinom.mle(x,m,n,tol)
    }
    switch(ci.method,
           "skew-score" = { ci.p <- pooledbinom.cscore.ci(x,m,n,tol,alpha)[2:3]},
           "bc-skew-score" ={ ci.p <- pooledbinom.bcscore.ci(x,m,n,tol,alpha)[2:3]},
           "score" = {ci.p <- pooledbinom.score.ci(x,m,n,tol,alpha)[2:3]},
           "lrt" = {ci.p <- pooledbinom.lrt.ci(x,m,n,tol,alpha)[2:3]},
           "wald" = {ci.p <- pooledbinom.wald.ci(x,m,n,tol,alpha)[2:3]},
           "mir" = {ci.p <- pooledbinom.mir.ci(x,m,n,tol,alpha)[2:3]}
    )
    list(p=p,lcl=ci.p[1],ucl=ci.p[2],pt.method=pt.method,ci.method=ci.method,alpha=alpha,x=x,m=m,n=n,scale=scale)
  }


"print.pooledBin" <- function(x, ...){
    print(as.data.frame(unclass(x)),...)
    invisible(x)
  }


# to build the package, had to use a version of "[.data.frame" without the .Internal(copyDFattr(xx,x)) call
# this just copied x as a list along with the attributes of the data frame
#"[.pooledBin" <- subsetDF #get("[.data.frame")
#"[.summary.pooledBin" <- subsetDF # get("[.data.frame")
#"[.pooledBinList" <- get("[.data.frame")

"[.pooledBin" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
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

"[.summary.pooledBin" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
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








"summary.pooledBin" <- function(object, simple = FALSE, ...){
  grp.names <- as.character(attr(object,"group.names"))
  "sumf" <- function(x) c(x, N = sum(x$n * x$m), NumPools = sum(x$n), NumPosPools = sum(x$x),
                          PtEstName = x$pt.method,
                          CIEstName = x$ci.method,
                          Alpha = x$alpha) # c() just adds to the list x
  out <- lapply(attr(object,"fullList"), sumf)
  #attributes(out) <- list(class = "summary.pooledBinList", names = attr(object,"names"),group.var = attr(object,"group.var"))
  #attributes(out) <- list(class = "summary.pooledBinList",
  #                        names = grp.names,group.var = attr(object,"group.var"))

  n <- length(out)
  group.dat <- as.data.frame(attr(object,"group.dat"))
  if(ncol(group.dat)==0){
    group.dat <- data.frame(Group=0)
  } else {
    names(group.dat) <- names(object)[1:length(attr(object,"group.var"))]
  }

  ans <- cbind(group.dat,
               data.frame(
                    P = rep(0,n),
                    Lower = rep(0,n),
                    Upper = rep(0,n),
                    Scale = rep(1,n),
                    N = rep(0,n),
                    NumPools = rep(0,n),
                    NumPosPools = rep(0,n))
  )

  # ans <- data.frame(Group = grp.names,
  #                   P = rep(0,n),
  #                   Lower = rep(0,n),
  #                   Upper = rep(0,n),
  #                   Scale = rep(1,n),
  #                   N = rep(0,n),
  #                   NumPools = rep(0,n),
  #                   NumPosPools = rep(0,n))

  #if(!is.null(attr(object,"group.var"))) names(ans)[1] <- attr(object,"group.var")
  #if(attr(object,"group.var") != "") names(ans)[1:length(attr(object,"group.var"))] <- attr(object,"group.var")

  for(i in 1:n)
    ans[i,(length(attr(object,"group.var"))-1) + (2:8)] <- c(out[[i]]$scale * c(out[[i]]$p, out[[i]]$lcl,out[[i]]$ucl), out[[i]]$scale,
                    out[[i]]$N, out[[i]]$NumPools, out[[i]]$NumPosPools)

  if(attr(object,"n.groups") == 1) ans$Group <- NULL

  if(simple) return(ans)

  structure(ans, class = "summary.pooledBin",
            df = as.data.frame(unclass(ans)),
            fullList = attr(object,"fullList"),
            #names = grp.names,
            group.var = attr(object,"group.var"),
            call = attr(object,"call"))
}

"print.summary.pooledBin" <- function(x, ...){

  cat("\nEstimation of Binomial Proportion for Pooled Data\n")
  cat(paste0("\nCall: ", deparse(attr(x,"call"),width.cutoff = 100),"\n\n"))
  switch(attr(x,"fullList")[[1]]$pt.method,
         "firth"  = PtEstName <- "Firth's Correction",
         "gart"   = PtEstName <- "Gart's Correction",
         "bc-mle" = PtEstName <- "Gart's Correction",
         "mle"    = PtEstName <- "Maximum Likelihood",
         "mir"    = PtEstName <- "Minimum Infection Rate"
  )
  switch(attr(x,"fullList")[[1]]$ci.method,
         "skew-score"    = CIEstName <- "Skew-Corrected Score (Gart)",
         "bc-skew-score" = CIEstName <- "Bias- & Skew-Corrected Score (Gart)",
         "score"         = CIEstName <- "Score",
         "lrt"           = CIEstName <- "Likelihood Ratio Test Inversion",
         "wald"          = CIEstName <- "Wald",
         "mir"           = CIEstName <- "Minimum Infection Rate"
  )
  cat(paste("Point estimator        :",PtEstName,"\n"))
  cat(paste("CI method              :",CIEstName,"\n"))
  cat(paste("Confidence coefficient : ",100*(1-attr(x,"fullList")[[1]]$alpha),"%\n\n",sep=""))

  print(attr(x,"df"))
  invisible(x)
}


"plot.pooledBin" <-
  function(x,pch=16,refline=TRUE,printR2=TRUE,layout=NULL,...){
    x.lst <- attr(x,"fullList")
    n.groups <- length(x.lst)
    n.groups.root <- ceiling(sqrt(n.groups))
    if(is.null(layout)) layout <- c(n.groups.root, n.groups.root)
    par(mfrow = layout)
    for(i in 1:n.groups){
      # reference: Chen & Swallow
      if(all(x.lst[[i]]$n==1)) {
        xmn <- as.list(by(x.lst[[i]]$x, x.lst[[i]]$m, function(x) c(sum(x),length(x))))
        m <- as.numeric(names(xmn))
        xx <- sapply(xmn,function(x) x[1])
        n <- sapply(xmn,function(x) x[2])
      } else {
        xx <- x.lst[[i]]$x
        m <- x.lst[[i]]$m
        n <- x.lst[[i]]$n
      }
      y <- log((xx+0.5)/(n+0.5))
      cc <- lm(y ~ m)
      plot(m,y,xlab="Pool Size",ylab="log((x+0.5)/(n+0.5))",...)
      if(refline) abline(cc)
      #if(printR2) cat(paste("R-squared for diagnostic line fit =",round(summary(cc)$r.squared,4),"\n"))
      #if(printR2) title(expression(paste(names(x)[[i]],"\n(",plain(R)^2,"=",round(summary(cc)$r.squared,2),")")),cex=0.75)
      #if(printR2) title(bquote(atop(.(names(x)[[i]]),"(" ~ R^2 ~ "=" ~ .(round(summary(cc)$r.squared,2)) ~ ")")), cex = 0.75)
      #if(printR2) title(bquote(.(names(x)[[i]]) ~ ":" ~ R^2 ~ "=" ~ .(round(summary(cc)$r.squared,2))), cex = 0.75)
      if(n.groups > 1){
        title.group <- as.character(attr(x,"group.names"))[[i]]
        if(printR2) title(bquote(.(title.group) ~ ":" ~ R^2 ~ "=" ~ .(round(summary(cc)$r.squared,2))), cex = 0.75)
          else title(bquote(.(title.group)), cex = 0.75)
      } else {
        if(printR2) title(bquote(R^2 ~ "=" ~ .(round(summary(cc)$r.squared,2))), cex = 0.75)
          else title("Model Diagnostic",cex=0.75)
      }
    }
    invisible(x)
  }

"as.data.frame.pooledBin" <- function(x){
  as.data.frame(unclass(x))
}
