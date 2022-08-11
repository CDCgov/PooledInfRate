
"vectorIndex" <- function(x,...){
  UseMethod("vectorIndex")
}
"VI" <- function(x,...){
  UseMethod("VI")
}


"vectorIndex.default" <- function(x,m,n=rep(1,length(x)), vector, trap.time=rep(1,length(x)), group,
                         n.use.traps = TRUE,
                         n.use.na = FALSE,
                         pt.method = c("firth","gart","bc-mle","mle","mir"),
                         ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                         scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
  call <- match.call()
  call[[1]] <- as.name("VI")

  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  if(!missing(group)) group.var <- deparse(substitute(group))
  else group.var = ""
  if(missing(vector)) stop("Must have vector specified to use vectorIndex().\n")
  vectors.var <- deparse(substitute(vector))
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  if(missing(trap.time)){
    traptime.var <- "VectorTime"
  } else {
    traptime.var <- deparse(substitute(trap.time))
  }
  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

  if(!missing(group)){
    group.var <- deparse(substitute(group))
    groups <- unique(group)
    nGroups <- length(groups)
    group.dat <- data.frame(Group=group)
  } else {
    group.var <= ""
    group <- rep("SINGLE",length(x))
    groups <- unique(group)
    nGroups <- 1
    group.dat <- data.frame()
  }


  ans <- vector(mode="list",length=nGroups)
  ans.p <- vector(mode="list",length=nGroups)
  names(ans.p) <- groups
  ans.nbar <- vector(mode="list",length=nGroups)
  names(ans.nbar) <- groups
  for(i in 1:nGroups){
    ans[[i]] <- 0
    ans.p[[i]] <- vector(mode="list",length=nVectors)
    names(ans.p[[i]]) <- Vectors
    ans.nbar[[i]] <- vector(mode="list",length=nVectors)
    names(ans.nbar[[i]]) <- Vectors
    for(j in 1:nVectors){
      sub <- (group==groups[i]) & (vector == Vectors[j]) & (m>0) & (n>0) # can't do this above because of NAs
      tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                              pt.method=pt.method,
                              ci.method=ci.method,
                              scale=scale,alpha=alpha,tol=tol)
      ans.p[[i]][[j]] <- tmp.pb$p
      # include pools not tested (x=NA) in the estimate of abundance
      # but only if there are missing values in x
      if(x.with.na & n.use.na){
        if(n.use.traps){
          ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                   sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
        } else {
          ans.nbar[[i]][[j]] <-  (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
        }
      } else {
        ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
      }

      ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
    }
  }
  ans.lst <- ans

  #ans <- data.frame(Group = groups,
  #                  VI = rep(0,nGroups))
  ans <- cbind(group.dat[!duplicated(group),],
               data.frame(
                    VI = rep(0,nGroups))
  )

  for(i in 1:nGroups) ans$VI[i] <- ans.lst[[i]]

  if(group.var != "") names(ans)[1] <- group.var

  if(nGroups == 1) ans$Group <- NULL
  # fix row names after subsetting
  group.dat <- as.data.frame(group.dat[!duplicated(group),])
  rownames(group.dat) <- 1:nrow(group.dat)

  ans <- structure(ans,class = "vectorIndex", fullList = ans.lst, group.names = groups, group.var = group.var,
                   n.groups = nGroups, group.dat = group.dat,
                   vector = Vectors, vectors.var = vectors.var, traptime.var = traptime.var, #deparse(substitute(vectors)),
                   p = ans.p, n = ans.nbar,
                   n.use.na = n.use.na, n.use.traps = n.use.traps,
                   pt.method = pt.method, ci.method = ci.method,
                   call=call)

  ans
}


"vectorIndex.formula" <- function(x, data,
                                  n.use.traps = TRUE,
                                  n.use.na = FALSE,
                                                  pt.method = c("firth","gart","bc-mle","mle","mir"),
                                                  ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                                                  scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
  call <- match.call()
  call[[1]] <- as.name("vectorIndex")
  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  #print(call)
  #call[[1]] <- as.name(strsplit(as.character(call[[1]]),"\\.")[[1]][1])
  if(missing(data))
    data <- environment(x)


  vars <- VIParseFormula(x,data)
  # -4 is the grouping variable, which can (now) have length > 1
  if(any(sapply(vars,length)[-4] >1)) stop("only variable names permitted in formula; perhaps use the default call")

  if(is.na(vars$vector) | is.null(vars$vector) | vars$vector=="")
    stop("Must have vectors specified to use vectorIndex().\n")
  vector <- eval(parse(text=vars$vector),data)
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  x <- eval(parse(text=vars$x), data)
  m <- eval(parse(text=vars$m), data)
  if(!is.null(vars$n))
    n <- eval(parse(text=vars$n), data)
  else n <- rep(1,length(x)) # default n
  if(vars$traptime == ""){
    trap.time <- rep(1,length(x))
  } else{
    trap.time <- eval(parse(text=vars$traptime),data)
  }
  if(!is.null(vars$group[1])){
    # NEW as of 5/27/2022 to allow multiple grouping variables
    n.group.var <- length(vars$group)
    group.dat <- as.data.frame(matrix(nrow=length(x),ncol=n.group.var))
    names(group.dat) <- vars$group

    for(i in 1:n.group.var) group.dat[,i] <- eval(parse(text=vars$group[i]),data)

    #print(interaction(group.dat,drop=TRUE)) # drop unused levels
    #return(group.dat)

    # end of NEW



    group.var <- vars$group
    #group <- eval(parse(text=vars$group), data)
    group <- interaction(group.dat, drop = TRUE)
    groups <- unique(group)
    nGroups <- length(groups)
  }  else{
    group.var <- ""
    group <- rep("SINGLE",length(x))
    groups <- unique(group)
    nGroups <- 1
    group.dat <- data.frame()
  }

  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

    ans <- vector(mode="list",length=nGroups)
    ans.p <- vector(mode="list",length=nGroups)
    names(ans.p) <- groups
    ans.nbar <- vector(mode="list",length=nGroups)
    names(ans.nbar) <- groups


    for(i in 1:nGroups){
      ans[[i]] <- 0
      ans.p[[i]] <- vector(mode="list",length=nVectors)
      names(ans.p[[i]]) <- Vectors
      ans.nbar[[i]] <- vector(mode="list",length=nVectors)
      names(ans.nbar[[i]]) <- Vectors
      for(j in 1:nVectors){
        sub <- (group == groups[i]) & (vector == Vectors[j]) & (m>0) & (n>0) # can't do this above because of NAs
        tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                               pt.method=pt.method,
                               ci.method=ci.method,
                               scale=scale,alpha=alpha,tol=tol)
        ans.p[[i]][[j]] <- tmp.pb$p
        ans.nbar[[i]][[j]] <- sum(m[sub] * n[sub]) / sum(trap.time[sub])# number of each vector
        # include pools not tested (x=NA) in the estimate of abundance
        # but only if there are missing values in x
        if(x.with.na & n.use.na){
          if(n.use.traps){
          ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                   sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          } else {
            ans.nbar[[i]][[j]] <- (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
                                   (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          }
        } else {
          ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
        }

        # number of each vector scaled by collection effort ("trap nights")
        ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
      }
    }

    ans.lst <- ans

    ans <- cbind(group.dat[!duplicated(group),],
                 data.frame(VI = rep(0, nGroups))
    )
    #ans <- data.frame(Group = groups,
    #                  VI = rep(0,nGroups))

    for(i in 1:nGroups) ans$VI[i] <- ans.lst[[i]]

    if(group.var[1] != "") names(ans)[1:length(vars$group)] <- group.var

    if(nGroups == 1) ans$Group <- NULL

    # fix row names after subsetting
    group.dat <- as.data.frame(group.dat[!duplicated(group),])
    rownames(group.dat) <- 1:nrow(group.dat)

    ans <- structure(ans,class = "vectorIndex", fullList = ans.lst, group.names = groups,
                     group.var = vars$group, n.groups = nGroups,
                     group.dat = group.dat,
                     vector = Vectors, vectors.var = vars$vector, traptime.var=vars$traptime,
                     p = ans.p, n = ans.nbar,
                     n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method=ci.method,
                     call=call)

  ans
}



"VI.default" <- function(x,m,n=rep(1,length(x)), vector, trap.time=rep(1,length(x)), group,
                                  n.use.traps = TRUE,
                                  n.use.na = FALSE,
                                  pt.method = c("firth","gart","bc-mle","mle","mir"),
                                  ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                                  scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
  call <- match.call()
  call[[1]] <- as.name("VI")

  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  if(!missing(group)) group.var <- deparse(substitute(group))
  else group.var = ""
  if(missing(vector)) stop("Must have vector specified to use vectorIndex().\n")
  vectors.var <- deparse(substitute(vector))
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  if(missing(trap.time)){
    traptime.var <- "VectorTime"
  } else {
    traptime.var <- deparse(substitute(trap.time))
  }
  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

  if(!missing(group)){
    group.var <- deparse(substitute(group))
    groups <- unique(group)
    nGroups <- length(groups)
    group.dat <- data.frame(Group=group)
  } else {
    group.var <= ""
    group <- rep("SINGLE",length(x))
    groups <- unique(group)
    nGroups <- 1
    group.dat <- data.frame()
  }


  ans <- vector(mode="list",length=nGroups)
  ans.p <- vector(mode="list",length=nGroups)
  names(ans.p) <- groups
  ans.nbar <- vector(mode="list",length=nGroups)
  names(ans.nbar) <- groups
  for(i in 1:nGroups){
    ans[[i]] <- 0
    ans.p[[i]] <- vector(mode="list",length=nVectors)
    names(ans.p[[i]]) <- Vectors
    ans.nbar[[i]] <- vector(mode="list",length=nVectors)
    names(ans.nbar[[i]]) <- Vectors
    for(j in 1:nVectors){
      sub <- (group==groups[i]) & (vector == Vectors[j]) & (m>0) & (n>0) # can't do this above because of NAs
      tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                              pt.method=pt.method,
                              ci.method=ci.method,
                              scale=scale,alpha=alpha,tol=tol)
      ans.p[[i]][[j]] <- tmp.pb$p
      # include pools not tested (x=NA) in the estimate of abundance
      # but only if there are missing values in x
      if(x.with.na & n.use.na){
        if(n.use.traps){
          ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                   sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
        } else {
          ans.nbar[[i]][[j]] <-  (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
        }
      } else {
        ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
      }

      ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
    }
  }
  ans.lst <- ans

  #ans <- data.frame(Group = groups,
  #                  VI = rep(0,nGroups))
  ans <- cbind(group.dat[!duplicated(group),],
               data.frame(
                 VI = rep(0,nGroups))
  )

  for(i in 1:nGroups) ans$VI[i] <- ans.lst[[i]]

  if(group.var != "") names(ans)[1] <- group.var

  if(nGroups == 1) ans$Group <- NULL
  # fix row names after subsetting
  group.dat <- as.data.frame(group.dat[!duplicated(group),])
  rownames(group.dat) <- 1:nrow(group.dat)

  ans <- structure(ans,class = "VI", fullList = ans.lst, group.names = groups, group.var = group.var,
                   n.groups = nGroups, group.dat = group.dat,
                   vector = Vectors, vectors.var = vectors.var, traptime.var = traptime.var, #deparse(substitute(vectors)),
                   p = ans.p, n = ans.nbar,
                   n.use.na = n.use.na, n.use.traps = n.use.traps,
                   pt.method = pt.method, ci.method = ci.method,
                   call=call)

  ans





}




"VI.formula" <- function(x, data,
                                  n.use.traps = TRUE,
                                  n.use.na = FALSE,
                                  pt.method = c("firth","gart","bc-mle","mle","mir"),
                                  ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                                  scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {

  call <- match.call()
  call[[1]] <- as.name("vectorIndex")
  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  #print(call)
  #call[[1]] <- as.name(strsplit(as.character(call[[1]]),"\\.")[[1]][1])
  if(missing(data))
    data <- environment(x)


  vars <- VIParseFormula(x,data)
  # -4 is the grouping variable, which can (now) have length > 1
  if(any(sapply(vars,length)[-4] >1)) stop("only variable names permitted in formula; perhaps use the default call")

  if(is.na(vars$vector) | is.null(vars$vector) | vars$vector=="")
    stop("Must have vectors specified to use vectorIndex().\n")
  vector <- eval(parse(text=vars$vector),data)
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  x <- eval(parse(text=vars$x), data)
  m <- eval(parse(text=vars$m), data)
  if(!is.null(vars$n))
    n <- eval(parse(text=vars$n), data)
  else n <- rep(1,length(x)) # default n
  if(vars$traptime == ""){
    trap.time <- rep(1,length(x))
  } else{
    trap.time <- eval(parse(text=vars$traptime),data)
  }
  if(!is.null(vars$group[1])){
    # NEW as of 5/27/2022 to allow multiple grouping variables
    n.group.var <- length(vars$group)
    group.dat <- as.data.frame(matrix(nrow=length(x),ncol=n.group.var))
    names(group.dat) <- vars$group

    for(i in 1:n.group.var) group.dat[,i] <- eval(parse(text=vars$group[i]),data)

    #print(interaction(group.dat,drop=TRUE)) # drop unused levels
    #return(group.dat)

    # end of NEW



    group.var <- vars$group
    #group <- eval(parse(text=vars$group), data)
    group <- interaction(group.dat, drop = TRUE)
    groups <- unique(group)
    nGroups <- length(groups)
  }  else{
    group.var <- ""
    group <- rep("SINGLE",length(x))
    groups <- unique(group)
    nGroups <- 1
    group.dat <- data.frame()
  }

  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

  ans <- vector(mode="list",length=nGroups)
  ans.p <- vector(mode="list",length=nGroups)
  names(ans.p) <- groups
  ans.nbar <- vector(mode="list",length=nGroups)
  names(ans.nbar) <- groups


  for(i in 1:nGroups){
    ans[[i]] <- 0
    ans.p[[i]] <- vector(mode="list",length=nVectors)
    names(ans.p[[i]]) <- Vectors
    ans.nbar[[i]] <- vector(mode="list",length=nVectors)
    names(ans.nbar[[i]]) <- Vectors
    for(j in 1:nVectors){
      sub <- (group == groups[i]) & (vector == Vectors[j]) & (m>0) & (n>0) # can't do this above because of NAs
      tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                              pt.method=pt.method,
                              ci.method=ci.method,
                              scale=scale,alpha=alpha,tol=tol)
      ans.p[[i]][[j]] <- tmp.pb$p
      ans.nbar[[i]][[j]] <- sum(m[sub] * n[sub]) / sum(trap.time[sub])# number of each vector
      # include pools not tested (x=NA) in the estimate of abundance
      # but only if there are missing values in x
      if(x.with.na & n.use.na){
        if(n.use.traps){
          ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                   sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
        } else {
          ans.nbar[[i]][[j]] <- (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
        }
      } else {
        ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
      }

      # number of each vector scaled by collection effort ("trap nights")
      ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
    }
  }

  ans.lst <- ans

  ans <- cbind(group.dat[!duplicated(group),],
               data.frame(VI = rep(0, nGroups))
  )
  #ans <- data.frame(Group = groups,
  #                  VI = rep(0,nGroups))

  for(i in 1:nGroups) ans$VI[i] <- ans.lst[[i]]

  if(group.var[1] != "") names(ans)[1:length(vars$group)] <- group.var

  if(nGroups == 1) ans$Group <- NULL

  # fix row names after subsetting
  group.dat <- as.data.frame(group.dat[!duplicated(group),])
  rownames(group.dat) <- 1:nrow(group.dat)

  ans <- structure(ans,class = "VI", fullList = ans.lst, group.names = groups,
                   group.var = vars$group, n.groups = nGroups,
                   group.dat = group.dat,
                   vector = Vectors, vectors.var = vars$vector, traptime.var=vars$traptime,
                   p = ans.p, n = ans.nbar,
                   n.use.na = n.use.na, n.use.traps = n.use.traps,
                   pt.method = pt.method, ci.method=ci.method,
                   call=call)

  ans



}






"VIParseFormula" <- function(formula,data=NULL){
  # this function reads a formula of the form
  # x ~ m(m.var) + n(n.var) | group / vector
  # or
  # x ~ m(m.var) + n(n.var) / vector
  # for a call to VI and returns
  # a list with elements (x, m, n, group, vector), which are character strings
  # with the respective names for those variables;
  # missing n() or group variables are set to NULL;
  # if a single variable is used it is assumed to be for m(), in which case
  # the m() construct may be excluded;
  # a grouping variable is optional;
  # errors for ambiguous calls are issued;

  # first determine if there's a grouping variable, indicated by "|" in the formula
  # note: grepl returns logicals for inclusion, whereas grep does not
  has.group <- any(grepl("\\|",formula))
  has.vectors <- any(grepl("/",formula))
  has.traptime <- any(grepl(":",formula))

  if(!has.vectors) stop("formula must specify the variable for the 'vectors' using / after the base formula")

  # pick off the bits before the / to get the vector & traptime variables
  vectors.var <- unlist(strsplit(as.character(formula),"/"))[[4]]
  if(has.traptime) {
    splt <- strsplit(vectors.var,":")
    vectors.var <- splt[[1]][[1]]
    traptime.var <- splt[[1]][[2]]
  } else {
    traptime.var <- ""
  }


  # strip off the bits after the / symbol
  # ...this leaves a formula that has the same structure as the pooledBin (or pIR) ones
  # ...so we can use that formula deparser to get the correct names
  formula <- as.formula(paste0(unlist(strsplit(as.character(formula), "/"))[-4][c(2,1,3)],collapse=" "))
  #
  lst <- pooledBinParseFormula(formula,data)
  # add the vectors variable names
  lst$vector <- vectors.var
  lst$traptime <- traptime.var
  lst
}

"print.VI" <- function(x, ...){
  print(as.data.frame(unclass(x)),...)

  invisible(x)
}

"print.vectorIndex" <- function(x, ...){
  print(as.data.frame(unclass(x)),...)
  invisible(x)
}


# to build the package, had to use a version of "[.data.frame" without the .Internal(copyDFattr(xx,x)) call
# this just copied x as a list along with the attributes of the data frame
#"[.pooledBin" <- subsetDF #get("[.data.frame")
#"[.summary.pooledBin" <- subsetDF # get("[.data.frame")
#"[.pooledBinList" <- get("[.data.frame")

"[.VI" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
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


"[.vectorIndex" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
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



"summary.vectorIndex" <- function(object, simple=FALSE, ...){
  x <- object
  nGroups <- attr(x,"n.groups")
  nVectors <- length(attr(x,"vector"))
  out.lst <-  list(vi=x$VI,vector=attr(x,"vectors"),
                 p=attr(x,"p"),n=attr(x,"n"),nGroups = attr(x,"n.groups"), nVectors =length(attr(x,"vector")),
                 vi.obj=object)

  if(nGroups > 1){
    #out.dat <- data.frame(rep(attr(out.lst$vi.obj,"group.names"),each=nVectors),
    #                      rep(attr(out.lst$vi.obj, "vector"),nGroups),
    #                      as.vector(unlist(out.lst$n)),
    #                      as.vector(unlist(out.lst$p)),
    #                      as.vector(unlist(out.lst$n) * unlist(out.lst$p)))
    group.dat <- attr(object,"group.dat")
    out.dat <- data.frame(group.dat[rep(1:nrow(group.dat),each=nVectors),],
                          rep(attr(out.lst$vi.obj, "vector"),nGroups),
                          as.vector(unlist(out.lst$n)),
                          as.vector(unlist(out.lst$p)),
                          as.vector(unlist(out.lst$n) * unlist(out.lst$p)))
    names(out.dat) <- c(attr(out.lst$vi.obj,"group.var"),attr(out.lst$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")
    rownames(out.dat) <- 1:nrow(out.dat)
  } else {
    out.dat <- data.frame(rep(attr(out.lst$vi.obj, "vector"),nGroups),
                          as.vector(unlist(out.lst$n)),
                          as.vector(unlist(out.lst$p)),
                          as.vector(unlist(out.lst$n) * unlist(out.lst$p)))
    names(out.dat) <- c(attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")

  }

  if(simple) return(out.dat)

  structure(out.lst, class="summary.vectorIndex",out.dat=out.dat)

}

"print.summary.vectorIndex" <- function(x, ...){
  nGroups <- x$nGroups
  nVectors <- x$nVectors
  out.dat <- attr(x,"out.dat")
  # if(nGroups > 1){
  #   out.dat <- data.frame(rep(attr(x$vi.obj,"group.names"),each=nVectors),
  #                         rep(attr(x$vi.obj, "vector"),nGroups),
  #                         as.vector(unlist(x$n)),
  #                         as.vector(unlist(x$p)),
  #                         as.vector(unlist(x$n) * unlist(x$p)))
  #   names(out.dat) <- c(attr(x$vi.obj,"group.var"),attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")
  # } else {
  #   out.dat <- data.frame(rep(attr(x$vi.obj, "vector"),nGroups),
  #                         as.vector(unlist(x$n)),
  #                         as.vector(unlist(x$p)),
  #                         as.vector(unlist(x$n) * unlist(x$p)))
  #   names(out.dat) <- c(attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")
  #
  # }

  cat(paste("\nCall: ", deparse(attr(x$vi.obj,"call"),width.cutoff=100,nlines=1),"\n\n"))

  cat(paste0("Prevalence estimate method : ",attr(x$vi.obj,"pt.method"),"\n"))
  cat(paste0("Use non-tested individuals in abundance estimate : ",attr(x$vi.obj,"n.use.na"),"\n"))
  cat(paste0("Use tested pools' pool sizes in abundance estimate : ",attr(x$vi.obj,"n.use.traps"),"\n"))
  cat("\n")

  print(x$vi.obj,...) # print the VIs by group

  cat("\n")
  cat("Detail by group and vector:\n\n")
  print(out.dat,...)
  invisible(x)
}


"summary.VI" <- function(object, simple=FALSE, ...){
  x <- object
  nGroups <- attr(x,"n.groups")
  nVectors <- length(attr(x,"vector"))
  out.lst <-  list(vi=x$VI,vector=attr(x,"vectors"),
                   p=attr(x,"p"),n=attr(x,"n"),nGroups = attr(x,"n.groups"), nVectors =length(attr(x,"vector")),
                   vi.obj=object)

  if(nGroups > 1){
    #out.dat <- data.frame(rep(attr(out.lst$vi.obj,"group.names"),each=nVectors),
    #                      rep(attr(out.lst$vi.obj, "vector"),nGroups),
    #                      as.vector(unlist(out.lst$n)),
    #                      as.vector(unlist(out.lst$p)),
    #                      as.vector(unlist(out.lst$n) * unlist(out.lst$p)))
    group.dat <- attr(object,"group.dat")
    out.dat <- data.frame(group.dat[rep(1:nrow(group.dat),each=nVectors),],
                          rep(attr(out.lst$vi.obj, "vector"),nGroups),
                          as.vector(unlist(out.lst$n)),
                          as.vector(unlist(out.lst$p)),
                          as.vector(unlist(out.lst$n) * unlist(out.lst$p)))
    names(out.dat) <- c(attr(out.lst$vi.obj,"group.var"),attr(out.lst$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")
    rownames(out.dat) <- 1:nrow(out.dat)
  } else {
    out.dat <- data.frame(rep(attr(out.lst$vi.obj, "vector"),nGroups),
                          as.vector(unlist(out.lst$n)),
                          as.vector(unlist(out.lst$p)),
                          as.vector(unlist(out.lst$n) * unlist(out.lst$p)))
    names(out.dat) <- c(attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")

  }

  if(simple) return(out.dat)

  structure(out.lst, class="summary.VI",out.dat=out.dat)

}

"print.summary.VI" <- function(x, ...){
  nGroups <- x$nGroups
  nVectors <- x$nVectors
  out.dat <- attr(x, "out.dat")
  # if(nGroups > 1){
  #   out.dat <- data.frame(rep(attr(x$vi.obj,"group.names"),each=nVectors),
  #                         rep(attr(x$vi.obj, "vector"),nGroups),
  #                         as.vector(unlist(x$n)),
  #                         as.vector(unlist(x$p)),
  #                         as.vector(unlist(x$n) * unlist(x$p)))
  #   names(out.dat) <- c(attr(x$vi.obj,"group.var"),attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")
  # } else {
  #   out.dat <- data.frame(rep(attr(x$vi.obj, "vector"),nGroups),
  #                         as.vector(unlist(x$n)),
  #                         as.vector(unlist(x$p)),
  #                         as.vector(unlist(x$n) * unlist(x$p)))
  #   names(out.dat) <- c(attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")
  #
  # }

  cat(paste("\nCall: ", deparse(attr(x$vi.obj,"call"),width.cutoff=100,nlines=1),"\n\n"))

  cat(paste0("Prevalence estimate method : ",attr(x$vi.obj,"pt.method"),"\n"))
  cat(paste0("Use non-tested individuals in abundance estimate : ",attr(x$vi.obj,"n.use.na"),"\n"))
  cat(paste0("Use tested pools' pool sizes in abundance estimate : ",attr(x$vi.obj,"n.use.traps"),"\n"))
  cat("\n")

  print(x$vi.obj,...) # print the VIs by group

  cat("\n")
  cat("Detail by group and vector:\n\n")
  print(out.dat,...)
  invisible(x)
}

"as.data.frame.vectorIndex" <- function(x){
  as.data.frame(unclass(x))
}

"as.data.frame.VI" <- function(x){
  as.data.frame(unclass(x))
}

# "[.summary.vectorIndex" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
#   mdrop <- missing(drop)
#   Narg <- nargs() - !mdrop
#   has.j <- !missing(j)
#   if (!all(names(sys.call()) %in% c("", "drop")) &&
#       !isS4(x))
#     warning("named arguments other than 'drop' are discouraged")
#   if (Narg < 3L) {
#     if (!mdrop)
#       warning("'drop' argument will be ignored")
#     if (missing(i))
#       return(x)
#     if (is.matrix(i))
#       return(as.matrix(x)[i])
#     nm <- names(x)
#     if (is.null(nm))
#       nm <- character()
#     if (!is.character(i) && anyNA(nm)) {
#       names(nm) <- names(x) <- seq_along(x)
#       y <- NextMethod("[")
#       cols <- names(y)
#       if (anyNA(cols))
#         stop("undefined columns selected")
#       cols <- names(y) <- nm[cols]
#     }
#     else {
#       y <- NextMethod("[")
#       cols <- names(y)
#       if (!is.null(cols) && anyNA(cols))
#         stop("undefined columns selected")
#     }
#     if (anyDuplicated(cols))
#       names(y) <- make.unique(cols)
#     attr(y, "row.names") <- .row_names_info(x, 0L)
#     attr(y, "class") <- oldClass(x)
#     return(y)
#   }
#   if (missing(i)) {
#     if (drop && !has.j && length(x) == 1L)
#       return(.subset2(x, 1L))
#     nm <- names(x)
#     if (is.null(nm))
#       nm <- character()
#     if (has.j && !is.character(j) && anyNA(nm)) {
#       names(nm) <- names(x) <- seq_along(x)
#       y <- .subset(x, j)
#       cols <- names(y)
#       if (anyNA(cols))
#         stop("undefined columns selected")
#       cols <- names(y) <- nm[cols]
#     }
#     else {
#       y <- if (has.j)
#         .subset(x, j)
#       else x
#       cols <- names(y)
#       if (anyNA(cols))
#         stop("undefined columns selected")
#     }
#     if (drop && length(y) == 1L)
#       return(.subset2(y, 1L))
#     if (anyDuplicated(cols))
#       names(y) <- make.unique(cols)
#     nrow <- .row_names_info(x, 2L)
#     if (drop && !mdrop && nrow == 1L)
#       return(structure(y, class = NULL, row.names = NULL))
#     else {
#       attr(y, "class") <- oldClass(x)
#       attr(y, "row.names") <- .row_names_info(x,
#                                               0L)
#       return(y)
#     }
#   }
#   xx <- x
#   cols <- names(xx)
#   x <- vector("list", length(x))
#   #x <- .Internal(copyDFattr(xx, x))
#   x <- as.list(xx)
#   attributes(x) <- attributes(xx)
#   oldClass(x) <- attr(x, "row.names") <- NULL
#   if (has.j) {
#     nm <- names(x)
#     if (is.null(nm))
#       nm <- character()
#     if (!is.character(j) && anyNA(nm))
#       names(nm) <- names(x) <- seq_along(x)
#     x <- x[j]
#     cols <- names(x)
#     if (drop && length(x) == 1L) {
#       if (is.character(i)) {
#         rows <- attr(xx, "row.names")
#         i <- pmatch(i, rows, duplicates.ok = TRUE)
#       }
#       xj <- .subset2(.subset(xx, j), 1L)
#       return(if (length(dim(xj)) != 2L) xj[i] else xj[i,
#                                                       , drop = FALSE])
#     }
#     if (anyNA(cols))
#       stop("undefined columns selected")
#     if (!is.null(names(nm)))
#       cols <- names(x) <- nm[cols]
#     nxx <- structure(seq_along(xx), names = names(xx))
#     sxx <- match(nxx[j], seq_along(xx))
#   }
#   else sxx <- seq_along(x)
#   rows <- NULL
#   if (is.character(i)) {
#     rows <- attr(xx, "row.names")
#     i <- pmatch(i, rows, duplicates.ok = TRUE)
#   }
#   for (j in seq_along(x)) {
#     xj <- xx[[sxx[j]]]
#     x[[j]] <- if (length(dim(xj)) != 2L)
#       xj[i]
#     else xj[i, , drop = FALSE]
#   }
#   if (drop) {
#     n <- length(x)
#     if (n == 1L)
#       return(x[[1L]])
#     if (n > 1L) {
#       xj <- x[[1L]]
#       nrow <- if (length(dim(xj)) == 2L)
#         dim(xj)[1L]
#       else length(xj)
#       drop <- !mdrop && nrow == 1L
#     }
#     else drop <- FALSE
#   }
#   if (!drop) {
#     if (is.null(rows))
#       rows <- attr(xx, "row.names")
#     rows <- rows[i]
#     if ((ina <- anyNA(rows)) | (dup <- anyDuplicated(rows))) {
#       if (!dup && is.character(rows))
#         dup <- "NA" %in% rows
#       if (ina)
#         rows[is.na(rows)] <- "NA"
#       if (dup)
#         rows <- make.unique(as.character(rows))
#     }
#     if (has.j && anyDuplicated(nm <- names(x)))
#       names(x) <- make.unique(nm)
#     if (is.null(rows))
#       rows <- attr(xx, "row.names")[i]
#     attr(x, "row.names") <- rows
#     oldClass(x) <- oldClass(xx)
#   }
#   x
# }
#
#
# "[.summary.VI" <-  function (x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)  {
#   mdrop <- missing(drop)
#   Narg <- nargs() - !mdrop
#   has.j <- !missing(j)
#   if (!all(names(sys.call()) %in% c("", "drop")) &&
#       !isS4(x))
#     warning("named arguments other than 'drop' are discouraged")
#   if (Narg < 3L) {
#     if (!mdrop)
#       warning("'drop' argument will be ignored")
#     if (missing(i))
#       return(x)
#     if (is.matrix(i))
#       return(as.matrix(x)[i])
#     nm <- names(x)
#     if (is.null(nm))
#       nm <- character()
#     if (!is.character(i) && anyNA(nm)) {
#       names(nm) <- names(x) <- seq_along(x)
#       y <- NextMethod("[")
#       cols <- names(y)
#       if (anyNA(cols))
#         stop("undefined columns selected")
#       cols <- names(y) <- nm[cols]
#     }
#     else {
#       y <- NextMethod("[")
#       cols <- names(y)
#       if (!is.null(cols) && anyNA(cols))
#         stop("undefined columns selected")
#     }
#     if (anyDuplicated(cols))
#       names(y) <- make.unique(cols)
#     attr(y, "row.names") <- .row_names_info(x, 0L)
#     attr(y, "class") <- oldClass(x)
#     return(y)
#   }
#   if (missing(i)) {
#     if (drop && !has.j && length(x) == 1L)
#       return(.subset2(x, 1L))
#     nm <- names(x)
#     if (is.null(nm))
#       nm <- character()
#     if (has.j && !is.character(j) && anyNA(nm)) {
#       names(nm) <- names(x) <- seq_along(x)
#       y <- .subset(x, j)
#       cols <- names(y)
#       if (anyNA(cols))
#         stop("undefined columns selected")
#       cols <- names(y) <- nm[cols]
#     }
#     else {
#       y <- if (has.j)
#         .subset(x, j)
#       else x
#       cols <- names(y)
#       if (anyNA(cols))
#         stop("undefined columns selected")
#     }
#     if (drop && length(y) == 1L)
#       return(.subset2(y, 1L))
#     if (anyDuplicated(cols))
#       names(y) <- make.unique(cols)
#     nrow <- .row_names_info(x, 2L)
#     if (drop && !mdrop && nrow == 1L)
#       return(structure(y, class = NULL, row.names = NULL))
#     else {
#       attr(y, "class") <- oldClass(x)
#       attr(y, "row.names") <- .row_names_info(x,
#                                               0L)
#       return(y)
#     }
#   }
#   xx <- x
#   cols <- names(xx)
#   x <- vector("list", length(x))
#   #x <- .Internal(copyDFattr(xx, x))
#   x <- as.list(xx)
#   attributes(x) <- attributes(xx)
#   oldClass(x) <- attr(x, "row.names") <- NULL
#   if (has.j) {
#     nm <- names(x)
#     if (is.null(nm))
#       nm <- character()
#     if (!is.character(j) && anyNA(nm))
#       names(nm) <- names(x) <- seq_along(x)
#     x <- x[j]
#     cols <- names(x)
#     if (drop && length(x) == 1L) {
#       if (is.character(i)) {
#         rows <- attr(xx, "row.names")
#         i <- pmatch(i, rows, duplicates.ok = TRUE)
#       }
#       xj <- .subset2(.subset(xx, j), 1L)
#       return(if (length(dim(xj)) != 2L) xj[i] else xj[i,
#                                                       , drop = FALSE])
#     }
#     if (anyNA(cols))
#       stop("undefined columns selected")
#     if (!is.null(names(nm)))
#       cols <- names(x) <- nm[cols]
#     nxx <- structure(seq_along(xx), names = names(xx))
#     sxx <- match(nxx[j], seq_along(xx))
#   }
#   else sxx <- seq_along(x)
#   rows <- NULL
#   if (is.character(i)) {
#     rows <- attr(xx, "row.names")
#     i <- pmatch(i, rows, duplicates.ok = TRUE)
#   }
#   for (j in seq_along(x)) {
#     xj <- xx[[sxx[j]]]
#     x[[j]] <- if (length(dim(xj)) != 2L)
#       xj[i]
#     else xj[i, , drop = FALSE]
#   }
#   if (drop) {
#     n <- length(x)
#     if (n == 1L)
#       return(x[[1L]])
#     if (n > 1L) {
#       xj <- x[[1L]]
#       nrow <- if (length(dim(xj)) == 2L)
#         dim(xj)[1L]
#       else length(xj)
#       drop <- !mdrop && nrow == 1L
#     }
#     else drop <- FALSE
#   }
#   if (!drop) {
#     if (is.null(rows))
#       rows <- attr(xx, "row.names")
#     rows <- rows[i]
#     if ((ina <- anyNA(rows)) | (dup <- anyDuplicated(rows))) {
#       if (!dup && is.character(rows))
#         dup <- "NA" %in% rows
#       if (ina)
#         rows[is.na(rows)] <- "NA"
#       if (dup)
#         rows <- make.unique(as.character(rows))
#     }
#     if (has.j && anyDuplicated(nm <- names(x)))
#       names(x) <- make.unique(nm)
#     if (is.null(rows))
#       rows <- attr(xx, "row.names")[i]
#     attr(x, "row.names") <- rows
#     oldClass(x) <- oldClass(xx)
#   }
#   x
# }
