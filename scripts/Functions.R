#### Parameters ####

.figures = "%MYFOLDER%/figures/"
.data = "%MYFOLDER%/data/"
.plot_par = list(
  cex.axis = 0.8,
  cex.lab = 1.1,
  cex.main = 1.1,
  cex.sub = 1,
  family = "serif",
  pch = 19,
  las = 2
)
.ticks_par = list(xaxp = c(-50, 50, 50),
                  yaxp = c(-50, 50, 50))
options(scipen = 5)


#### Packages ####

.pkg = c(
  "RColorBrewer",
  "latex2exp",
  "mixtools",
  "mvnfast",
  "ellipse",
  "ICS",
  "ICSOutlier",
  "compositions",
  "pracma",
  "expm",
  "xtable"
)
lapply(.pkg, require, character.only = TRUE)
require(knitr)
write_bib(.pkg, "%MYFOLDER%/Rpackages.bib")


#### Data ####

load(paste0(.data, "BDDSegX.RData"))
BDDSegXvol = BDDSegX[c("V_A", "V_B", "V_C", "V_D", "V_E")]
BDDSegXcomp = acomp(BDDSegX[c("S_A", "S_B", "S_C", "S_D", "S_E")])


#### Utilities ####

.ones = function(n) {
  return(rep(1 / n, n) %*% t(rep(1, n)))
}

identity.acomp = function(D) {
  return(diag(D) - .ones(D))
}

.rotation = function(theta) {
  rbind(c(cos(theta),-sin(theta)), c(sin(theta), cos(theta)))
}

cov.acomp = function (x,
                      y = NULL,
                      ...,
                      robust = getOption("robust"),
                      use = "all.obs",
                      giveCenter = FALSE) {
  control <- attr(robust, "control")
  if (is.logical(robust))
    robust <- if (robust)
      "mcd"
  else
    "pearson"
  if (has.missings(x) || has.missings(y)) {
    if (!(is.character(robust) && robust == "pearson"))
      warning("var.*: Robust estimation with losts not yet implemented")
    if (is.null(y)) {
      if (use == "pairwise.complete.obs") {
        return(gsi.varwithlosts(cdt(x), giveCenter = giveCenter))
      }
      else {
        tk = as.logical(gsi.geometricmeanRow(is.NMV(x)))
        xaux = x[tk,]
        class(xaux) = class(x)
        return(var(cdt(xaux), giveCenter = giveCenter))
      }
    }
    else {
      warning("Covariance with losts not yet implemented. Omitting lost values.")
      tk = as.logical(gsi.geometricmeanRow(is.NMV(cbind(x,
                                                        y))))
      xaux = x[tk,]
      class(xaux) = class(x)
      yaux = y[tk,]
      class(yaux) = class(y)
      return(var(unclass(cdt(xaux)), unclass(cdt(yaux)),
                 ...))
    }
  }
  else {
    switch(robust, pearson = {
      erg <- var(unclass(cdt(x)), unclass(cdt(y)), ...)
      if (giveCenter)
        attr(erg, "center") <- mean(x, robust = FALSE)
      class(erg) = "acompmatrix"
      erg
    }, mcd = {
      if (is.null(y)) {
        erg <- ilrvar2clr(
          if (is.null(control))
            covMcd(unclass(idt(x)),
                   ...)$cov
          else
            covMcd(unclass(idt(x)), ...,
                   control = control)$cov,
          x = x
        )
        if (giveCenter)
          attr(erg, "center") <- mean(x,
                                      robust = FALSE)
        class(erg) = "acompmatrix"
        erg
      } else {
        Dx <- ncol(x)
        Dy <- ncol(y)
        x1 <- idt(x)
        y1 <- idt(y)
        erg <- if (is.null(control))
          covMcd(cbind(x1,
                       y1), ...)
        else
          covMcd(cbind(x1, y1), ...,
                 control = control)
        m <- erg$center
        erg <- erg$cov
        erg <- t(ilrBase(D = Dx)) %*% erg[1:(Dx - 1),
                                          Dx:ncol(erg)] %*% ilrBase(D = Dy)
        row.names(erg) <- colnames(x)
        colnames(erg) <- colnames(y)
        if (giveCenter)
          attr(erg, "center") <- idtInv(m,
                                        x)
        class(erg) = "acompmatrix"
        erg
      }
    }, stop("var.*: unkown robust method", robust))
  }
}

ilr.default = ilr

ilr = function(x) {
  UseMethod("ilr", x)
}

ilr.acompmatrix = function (A, V = ilrBase(D = ncol(A)), ...)
{
  t(V) %*% A %*% V
}

ilrInv.default = ilrInv

ilrInv = function(x) {
  UseMethod("ilrInv", x)
}

ilrInv.matrix = function (A, V = ilrBase(D = ncol(A) + 1), ...)
{
  res = V %*% A %*% t(V)
  class(res) = "acompmatrix"
  res
}

solve.acompmatrix = function(A) {
  D = nrow(A)
  B = solve.default(A + .ones(D)) - .ones(D)
  class(B) = "acompmatrix"
  B
}

sqrtm.matrix = sqrtm

sqrtm = function(x) {
  UseMethod("sqrtm", x)
}

sqrtm.acompmatrix = function(A) {
  class(A) = "matrix"
  B = sqrtm(A)
  class(B) = "acompmatrix"
  B
}

MeanCov.default = MeanCov

MeanCov = function(x)
  UseMethod("MeanCov", x)

MeanCov.acomp = function(x) {
  list(location = mean(x), scatter = cov(x))
}

Mean3Cov4.default = Mean3Cov4

Mean3Cov4 = function(x)
  UseMethod("Mean3Cov4", x)

Mean3Cov4.acomp = function(x) {
  V = ilrBase(D = ncol(x))
  list(location = ilrInv(mean3(ilr(x))),
       scatter = V %*% cov4(data.frame(ilr(x))) %*% t(V))
}

#### Plotting & printing ####

print_table = function(tab, file) {
  text = function(x) {
    paste0("\\mathtt{", x, "}")
    #x
  }
  print(
    xtable(tab),
    file = paste0(.figures, file, ".tex"),
    floating = FALSE,
    tabular.environment = "array",
    hline.after = NULL,
    include.rownames = T,
    include.colnames = T,
    sanitize.rownames.function = text,
    sanitize.colnames.function = text
    #latex.environments = ""
  )
}

grid.acomp = function(ilr_max = 6,
                      eps = 1,
                      A = diag(2),
                      b = c(0, 0),
                      data = NULL) {
  if (A == "cov") {
    A = sqrtm(cov(ilr(data)))
    b = mean(ilr(data))
  }
  for (i in (-ilr_max / eps):(ilr_max / eps)) {
    lines(ilrInv(A %*% rbind(c(-ilr_max, i * eps), c(ilr_max, i * eps)) + b),
          col = "lightgray",
          lty = "dotted")
    lines(ilrInv(A %*% rbind(c(i * eps,-ilr_max), c(i * eps, ilr_max)) +
                   b),
          col = "lightgray",
          lty = "dotted")
  }
}

.plot_data_make = function(data,
                           grid,
                           cov_ellipse,
                           cov4_ellipse,
                           col,
                           asp) {
  UseMethod(".plot_data_make", data)
}

.plot_data_make.default = function(data,
                                   grid,
                                   cov_ellipse,
                                   cov4_ellipse,
                                   col,
                                   asp) {
  stop("plot_data only works on objects whose class is 'data.frame' or 'acomp'.")
}

.plot_data_make.acomp = function(data,
                                 grid,
                                 cov_ellipse,
                                 cov4_ellipse,
                                 col,
                                 asp) {
  par(mar = c(0, 1, 0, 1))
  plot.acomp(data,
             labels = TeX(names(data)),
             col = col)
  if (cov_ellipse) {
    ellipses.acomp(mean(data),
                   var(data),
                   r = 3.5,
                   col = "blue")
    plot.acomp(t(mean(data)),
               col = "blue",
               pch = 1,
               add = TRUE)
  }
  if (grid) {
    grid.acomp()
  }
}

.plot_data_make.data.frame = function(data,
                                      grid,
                                      cov_ellipse,
                                      cov4_ellipse,
                                      col,
                                      asp) {
  if (ncol(data) == 2) {
    plot(
      data,
      asp = asp,
      xlab = TeX(names(data))[1],
      ylab = TeX(names(data))[2],
      col = col
    )
  }
  else {
    plot(data,
         labels = TeX(names(data)),
         col = col)
  }
  if (cov_ellipse) {
    lines(
      ellipse::ellipse(cov(data),
                       centre = colMeans(data),
                       level = 0.98),
      type = "l",
      lty = 2,
      col = "blue"
    )
    points(t(colMeans(data)), col = "blue", pch = 3)
  }
  if (cov4_ellipse) {
    lines(
      ellipse::ellipse(cov4(data),
                       centre = colMeans(data),
                       level = 0.98),
      type = "l",
      lty = 2,
      col = "red"
    )
  }
  if (grid) {
    #par(.ticks_par)
    #grid()
    if (class(data[, 1]) == "Date") {
      abline(
        h = pretty(extendrange(data[, 2])),
        v = axis.Date(1, data[, 1]),
        col = "lightgray",
        lty = "dotted"
      )
    }
    else{
      abline(
        h = pretty(extendrange(data[, 2])),
        v = pretty(extendrange(data[, 1])),
        col = "lightgray",
        lty = "dotted"
      )
    }
  }
}

plot_data = function(data_list,
                     grid = FALSE,
                     cov_ellipse = FALSE,
                     cov4_ellipse = FALSE,
                     col = rgb(0, 0, 0, 0.5),
                     asp = 1,
                     file = NULL,
                     width = 7,
                     height = 5) {
  if (class(data_list) != "list")
    data_list = list(data_list)
  #data_list = lapply(data_list, data.frame)
  
  if (!is.null(file)) {
    pdf(
      file = paste0(.figures, file, ".pdf"),
      width = width,
      height = height
    )
  }
  
  if (col == "Blues") {
    col = colorRampPalette(rev(brewer.pal(9, "Blues")))(nrow(data_list[[1]]))
  }
  
  par(.plot_par)
  par(mar = c(5, 4, 2, 2))
  
  n = length(data_list)
  par(mfrow = c(1, n))
  for (i in 1:n) {
    data = data_list[[i]]
    .plot_data_make(data, grid, cov_ellipse, cov4_ellipse, col, asp)
  }
  
  if (!is.null(file)) {
    dev.off()
  }
}


#### Whitening ####

whiten = function(x)
  UseMethod("whiten", x)

whiten.data.frame = function (data, center = TRUE) {
  X = as.matrix(data)
  W = Re(sqrtm(solve(cov(X))))
  Z = data.frame(X %*% W)
  if (center)
    Z = sweep(Z, 2, colMeans(Z))
  colnames(Z) = paste(rep("Z", ncol(Z)), 1:ncol(Z), sep = "_")
  Z
}


whiten.acomp = function(data, center = TRUE) {
  X = data
  W = Re(sqrtm(as.matrix(solve(cov(
    X
  )))))
  class(W) = "acompmatrix"
  Z = clrInv(clr(X) %*% W)
  if (center)
    Z = Z - mean(Z)
  names(Z) = paste(rep("Z", ncol(Z)), 1:ncol(Z), sep = "_")
  Z
}


#### Mixtures CoDa ####

rmvnormmix.acomp = function (n,
                             lambda = 1,
                             mu = 0,
                             sigma = diag(length)) {
  m = length(lambda)
  if (is.null(D <- ncol(mu))) {
    D = length(mu)
  }
  if (is.null(nrow(mu)) || nrow(mu) != m) {
    mu = matrix(mu, m, D, TRUE)
  }
  sigma = matrix(sigma, m * D, D, TRUE)
  if (D != ncol(sigma)) {
    stop("mu and sigma must have the same number of columns",
         call. = FALSE)
  }
  sigma = lapply(split(sigma, rep(1:m, each = D, D)), function(x) {
    matrix(x, D, D, T)
  })
  z = sort(sample(m, n, replace = TRUE, prob = lambda), decreasing = TRUE)
  X = data.frame()
  for (k in 1:n) {
    X = rbind(X, as.vector(clrInv(
      ilr2clr(rmvn(1, rep(0, D - 1), diag(D - 1))) %*% Re(sqrtm(sigma[z][[k]])) + clr(mu[z, ][k, ])
    )))
  }
  names(X) = paste(rep("X", D), 1:D, sep = "_")
  list(data = X, index = z)
}


#### ICS CoDa ####

# ics.acomp = function (X,
#                       V = NULL,
#                       S1 = MeanCov,
#                       S2 = Mean3Cov4,
#                       S1args = list(),
#                       S2args = list(),
#                       na.action = na.fail)
# {
#   X <- na.action(X)
#   D <- ncol(X)
#   data.matrix <- acomp(X)
#   S1name <- deparse(substitute(S1))
#   S2name <- deparse(substitute(S2))
#   if (!is.list(S1))
#     S1 <- do.call("S1", c(list(data.matrix), S1args))
#   if (!is.list(S2))
#     S2 <- do.call("S2", c(list(data.matrix), S2args))
#   T1.X <- S1[[1]]
#   S1.X <- S1[[2]]
#   T2.X <- S2[[1]]
#   S2.X <- S2[[2]]
#   S1.X.eigen <- eigen(S1.X, symmetric = TRUE)
#   B1 <-
#     S1.X.eigen$vectors %*% tcrossprod(diag(S1.X.eigen$values ^ (-0.5)),
#                                       S1.X.eigen$vectors)
#   S2.Y <- B1 %*% S2.X %*% B1
#   S2.Y.eigen <- eigen(S2.Y, symmetric = TRUE)
#   U2 <- S2.Y.eigen$vectors
#   gKurt <- S2.Y.eigen$values
#   B <- crossprod(U2, B1)
#   T1.Z <- T1.X %*% B
#   T2.Z <- T2.X %*% B
#   gSkew <- T1.Z - T2.Z
#   skew.signs <- ifelse(gSkew >= 0, 1,-1)
#   B.res <- sweep(B, 1, skew.signs, "*")
#   data.matrix.C <- sweep(data.matrix, 2, T1.X, "-")
#   Z <- as.data.frame(tcrossprod(data.matrix.C, B.res))
#   names(Z) <- paste(rep("IC", p), 1:p, sep = ".")
#   if (is.null(colnames(X)) == TRUE)
#     names.X <- paste(rep("X", p), 1:p, sep = ".")
#   else
#     names.X <- colnames(X)
#   stdB <- "Z"
#   stdKurt <- FALSE
#   res <- new(
#     "ics2",
#     gKurt = gKurt,
#     UnMix = B.res,
#     S1 = S1.X,
#     S2 = S2.X,
#     S1name = S1name,
#     S2name = S2name,
#     S1args = S1args,
#     S2args = S2args,
#     Scores = Z,
#     T1 = T1.X,
#     T2 = T2.X,
#     gSkew = as.vector(skew.signs *
#                         gSkew),
#     DataNames = names.X,
#     StandardizeB = stdB,
#     StandardizegKurt = stdKurt
#   )
#   return(res)
# }

.V0 = function(D){
  V0 = matrix(1:(D * (D - 1)), D - 1, D)
  V0 = t(apply(V0,  c(1, 2), function(x) {
    i = which(V0 == x, arr.ind = T)[, 1]
    j = which(V0 == x, arr.ind = T)[, 2]
    if (j >= i + 1) {
      1 / sqrt((D - i) * (D - i + 1))
    }
    else if (j == i) {
      -sqrt((D - i) / (D - i + 1))
    }
    else{
      0
    }
  }))
  V0
}

ics.acomp = function(X,
                     S1 = "cov",
                     S2 = "cov4",
                     V = .V0(ncol(X))) {
  D = ncol(X)
  X = acomp(X)
  V0 = .V0(D)
  
  if (S1 == "cov") {
    S1 = cov(X)
  }
  if (S2 == "cov4") {
    S2 = V %*% cov4(data.frame(clr(X) %*% V)) %*% t(V)
  }
  
  W = Re(sqrtm(solve(cov(X))))
  A = W %*% S2 %*% W
  Q_s = eigen(t(V) %*% A %*% V, symmetric = TRUE)$vectors
  Q_s = t(V) %*% V0 %*% apply(t(V0) %*% V %*% Q_s, 1, function(x)
    if (x[1] < 0) {
      -x
    } else{
      x
    })
  Q = V %*% Q_s %*% t(V0)
  H = W %*% t(Q)
  Z = clrInv(clr(X - mean(X)) %*% H)
  colnames(Z) = paste(rep("IC", D), 1:D, sep = ".")
  class(H) = "acompmatrix"
  cov_clr = solve(t(H))
  cov_log = cov(log(data.frame(X)), log(data.frame(Z)))
  rownames(cov_clr) = colnames(X)
  colnames(cov_clr) = colnames(Z)
  rownames(cov_log) = colnames(X)
  list(
    data = Z,
    H = H,
    Q = Q,
    S2 = S2,
    kurtosis = eigen(t(V) %*% A %*% V, symmetric = TRUE)$values,
    cov_clr = cov_clr,
    cov_log = cov_log
  )
}