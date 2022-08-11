
"ipooledBin" <- function(x,...){
  UseMethod("ipooledBin")
}

"ipooledBin.default" <- function(x,m,n=rep(1,length(x)), group,
           pt.method = c("firth","mle"),
           sens=rep(1,length(x)),spec=rep(1,length(x)),  scale=1,tol=.Machine$double.eps^0.5,
           max.iter=10000,p.start=NULL,...) {
    call <- match.call()
    call[[1]] <- as.name("ipooledBin")
    pt.method <- match.arg(pt.method)

      if(!missing(group)){
        group.var <- deparse(substitute(group))
        xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n) & !is.na(group))
        group <- group[xmng.nomiss]
        groups <- unique(group)
        nGroups <- length(groups)
      } else {
        xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n))
        group <- rep("SINGLE",length(x[xmng.nomiss]))
        groups <- unique(group)
        nGroups <- length(groups)
        group.var <- ""
      }

      ans <- vector(mode="list",length=nGroups)

      for(i in 1:nGroups){
        switch(pt.method,
               firth = {
                tmp.p <- ipooledbinom.firth(x[group==groups[i]],
                                             m[group==groups[i]],
                                             n[group==groups[i]],
                                             sens=sens,spec=spec,tol=tol,
                                             max.iter=max.iter,p.start=p.start)
                },
                mle = {
                  tmp.p <- ipooledbinom.firth(x[group==groups[i]],
                                               m[group==groups[i]],
                                               n[group==groups[i]],
                                               sens=sens,spec=spec,tol=tol,
                                               max.iter=max.iter)
                  }
                  )

        ans[[i]] <- list(p=tmp.p,x=x[group==groups[i]],
                              m=m[group==groups[i]],
                              n=n[group==groups[i]],class="ipooledBin",
                         scale = scale,
                              call=call)
      }

      ans.lst <- ans

      ans <- data.frame(Group = groups,
                        P = rep(0,  nGroups),
                        Scale = rep(scale, nGroups))

      if(group.var != "") names(ans)[1] <- group.var

      for (i in 1:nGroups) ans[i, 2:3] <- ans.lst[[i]]$scale * c(ans.lst[[i]]$p, 1)
      if (all(ans$Scale == 1)) ans$Scale <- NULL


      if(nGroups == 1) ans$Group <- NULL

      ans <- structure(ans, class = "ipooledBin", fullList = ans.lst, groups= groups, group.var = group.var,
                       n.groups = nGroups,
                       x=x,m=m,n=n,
                       sens=sens,spec=spec,scale = scale, pt.method=pt.method,call=call)

    ans
  }


"ipooledBin.formula" <- function(x, data,
           pt.method = c("firth","mle"),
           sens = 1, spec = 1,
           scale = 1, tol=.Machine$double.eps^0.5,
           max.iter=10000,p.start=NULL,...){
    call <- match.call()
    call[[1]] <- as.name("ipooledBin")

    pt.method <- match.arg(pt.method)

    if(missing(data))
      data <- environment(x)

    # omit records with missing data -- note, if data contains records missing
    # anywhere (even not in X, M, N, Group, they are omitted), so care should be
    # used in subsetting before the call to be asured of the desired analysis
    #data <- na.omit(data)

    vars <- pooledBinParseFormula(x, data)
    if(any(sapply(vars,length)[1:3]>1)) stop("only variable names permitted in formula; perhaps use the default call")

    # retrieve values from data using the character name
    # use the eval(parse(text=XX), data) construct in case data is left unspecified,
    # and because we have the names of the variables available
    # from pooledBin.deparseFormula
    x <- eval(parse(text=vars$x), data)
    m <- eval(parse(text=vars$m), data)
    if(!is.null(vars$n))
      n <- eval(parse(text=vars$n), data)
    else n <- rep(1,length(x)) # default n

    if(!is.null(vars$group)){
      # NEW as of 5/27/2022 to allow multiple grouping variables
      n.group.var <- length(vars$group)
      group.dat <- as.data.frame(matrix(nrow=length(x),ncol=n.group.var))
      names(group.dat) <- vars$group

      for(i in 1:n.group.var) group.dat[,i] <- eval(parse(text=vars$group[i]),data)

      group.var <- vars$group
      #group <- eval(parse(text=vars$group), data)
      group <- interaction(group.dat, drop=TRUE)
      groups <- unique(group)
      nGroups <- length(groups)
    }  else{
      group <- rep("SINGLE",length(x))
      groups <- unique(group)
      nGroups <- 1
      group.var <- ""
      group.dat <- data.frame()
    }

      ans <- vector(mode="list",length=nGroups)
      for(i in 1:nGroups){
        switch(pt.method,
               firth = {
                 tmp.p <- ipooledbinom.firth(x[group==groups[i]],
                                              m[group==groups[i]],
                                              n[group==groups[i]],
                                              sens=sens,spec=spec,tol=tol,
                                              max.iter=max.iter,p.start=p.start)
               },
               mle = {
                 tmp.p <- ipooledbinom.mle(x[group==groups[i]],
                                              m[group==groups[i]],
                                              n[group==groups[i]],
                                              sens=sens,spec=spec,tol=tol,
                                              max.iter=max.iter)
               }
        )
        ans[[i]] <- list(p=tmp.p,x=x[group==groups[i]],
                              m=m[group==groups[i]],
                              n=n[group==groups[i]],scale=scale,class="ipooledBIn",call=call)

      }


      ans.lst <- ans

      ans <- cbind(group.dat[!duplicated(group),],
                   data.frame(
                        P = rep(0,  nGroups),
                        Scale = rep(scale, nGroups))
      )
      #ans <- data.frame(Group = groups,
      #                  P = rep(0,  nGroups),
      #                  Scale = rep(scale, nGroups))

      #if(group.var != "") names(ans)[1] <- group.var

      if (!is.null(vars$group)){
        #names(ans)[1] <- paste0(vars$group,collapse=".") # new as of 5/26/2022 for multiple grouping variables
        names(ans)[1:length(vars$group)] <- vars$group
      }

      for (i in 1:nGroups) ans[i, (ncol(group.dat)-1) + (2:3)] <- ans.lst[[i]]$scale * c(ans.lst[[i]]$p, 1)
      if (all(ans$Scale == 1)) ans$Scale <- NULL

      if(nGroups == 1) ans$Group <- NULL

      # fix row names after subsetting
      group.dat <- as.data.frame(group.dat[!duplicated(group),])
      rownames(group.dat) <- 1:nrow(group.dat)


      ans <- structure(ans, class = "ipooledBin", fullList = ans.lst, groups = groups, group.var = group.var,
                       n.groups = nGroups,group.dat = group.dat,
                       x=x,m=m,n=n,
                       sens=sens,spec=spec,scale = scale, pt.method=pt.method,call=call)

      ans

  }



"print.ipooledBin" <- function(x, ...){
    print(as.data.frame(unclass(x)),...)
    invisible(x)
  }

"summary.ipooledBin" <- function(object, ...){
  groups  <- attr(object,"groups")
  "sumf" <- function(x) c(P=x$p,N = sum(x$n * x$m), NumPools = sum(x$n), NumPosPools = sum(x$x))
  out <- lapply(attr(object,"fullList"), sumf)
  #attributes(out) <- list(class = "summary.pooledBinList", names = attr(object,"names"),group.var = attr(object,"group.var"))
  #attributes(out) <- list(class = "summary.pooledBinList",
  #                        names = grp.names,group.var = attr(object,"group.var"))

  group.dat <- as.data.frame(attr(object,"group.dat"))

  if(ncol(group.dat)==0){
    group.dat <- data.frame(Group=0)
  } else {
    names(group.dat) <- names(object)[1:length(attr(object,"group.var"))]
  }

  structure(out, class = "summary.ipooledBin",
            n.groups = attr(object,"n.groups"),
            groups = groups,
            group.var = attr(object,"group.var"),
            group.dat = group.dat,
            scale = attr(object,"scale"),
                       x=attr(object,"x"),m=attr(object,"m"),n=attr(object,"n"),
            sens = attr(object,"sens"),
            spec = attr(object,"spec"),
            pt.method = attr(object,"pt.method"),
            call = attr(object,"call"))
}

"print.summary.ipooledBin" <- function(x, ...){
  n <- length(x)
  out <- cbind(attr(x,"group.dat"),
               data.frame(
                    PointEst = rep(0,n),
                    N = rep(0,n),
                    NumPools = rep(0,n),
                    NumPosPools = rep(0,n),
                    Scale = rep(1,n))
  )
  #out <- data.frame(Group = attr(x,"groups"),
  #                  PointEst = rep(0,n),
  #                  N = rep(0,n),
  #                  NumPools = rep(0,n),
  #                  NumPosPools = rep(0,n),
  #                  Scale = rep(1,n))
  print(length(attr(x,"group.var")))
  if(attr(x,"n.groups") != 1) names(out)[1:length(attr(x,"group.var"))] <- attr(x,"group.var")

  for(i in 1:n)
    out[i,(length(attr(x,"group.var"))-1) + (2:6)] <- c(attr(x,"scale") * x[[i]][["P"]], x[[i]][["N"]], x[[i]][["NumPools"]], x[[i]][["NumPosPools"]], attr(x,"scale"))

  if(attr(x,"n.groups") == 1) out$Group <- NULL

  #if(is.null(digits)) digits <- 4
  #else digits <- args$digits
  cat("\n")
  cat("Estimation of Binomial Proportion for Pooled Data using an Imperfect Test\n")
  cat(paste0("\nCall: ", deparse(attr(x,"call"),width.cutoff=120),"\n\n"))
  if(all(attr(x,"sens")==1)) sens <- 1
  else sens <- attr(x,"sens")
  if(all(attr(x,"spec")==1)) spec <- 1
  else spec <- attr(x,"spec")
  cat(paste0("Sensitivity : ",paste0(" ",sens,collapse=","),"\n"))
  cat(paste0("Specificity : ",paste0(" ",spec,collapse=","),"\n"))
  cat(paste("Point estimator :",attr(x,"pt.method"),"\n"))
  cat("\n")
  print(out,...)
  invisible(x)
}


"[.ipooledBin" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
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

"[.summary.ipooledBin" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
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

"as.data.frame.ipooledBin" <- function(x){
  as.data.frame(unclass(x))
}


