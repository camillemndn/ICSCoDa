source("%MYFOLDER%/scripts/Functions.R")

## ---- 1.1.0_mixtures

.q = 3
.eps_0 = 0.9
.theta = runif(.q, 0, pi / 3)

mixture = data.frame(rmvnormmix(
  n = 300,
  lambda = c(.eps_0, rep((1 - .eps_0) / .q, .q)),
  mu = t(matrix(c(
    0, 0, 5 * rbind(cos(.theta), sin(.theta))
  ), 2, .q + 1)),
  sigma = t(matrix(c(1, 1, rep(
    0.2, 2 * .q
  )), 2, .q + 1))
))
names(mixture) = paste(rep("X", 2), 1:2, sep = "_")

plot_data(mixture,
          grid = TRUE,
          file = "1.1.0_mixture")


## ---- 1.1.1_elliptical

elliptical = data.frame(rmvn(n = 300, c(5, 1), rbind(c(2, 1), c(1, 2))))
names(elliptical) = paste(rep("X", 2), 1:2, sep = "_")

plot_data(elliptical,
          grid = TRUE,
          cov_ellipse = TRUE,
          file = "1.1.1_elliptical")


## ---- 1.1.2_cov4

plot_data(
  list(elliptical, mixture),
  grid = TRUE,
  cov_ellipse = TRUE,
  cov4_ellipse = TRUE,
  file = "1.1.2_cov4",
  width = 8,
  height = 3
)


## ---- 1.1.3_BDDSegXvol_whitening

plot_data(BDDSegXvol, col = "Blues", file = "1.1.3_BDDSegXvol")
plot_data(whiten(BDDSegXvol), col = "Blues", file = "1.1.3_BDDSegXvol_whitening")
print_table(cor(BDDSegXvol, whiten(BDDSegXvol)), "1.1.3_BDDSegXvol_whitening_cor")
#solve(diag(sqrt(diag(cov(BDDSegXvol))))) %*% sqrtm(cov(BDDSegXvol))


## ---- 1.1.4_BDDSegXvol_ICS

.BDDSegXvol_ICS = ics2(BDDSegXvol)

# IC
BDDSegXvol_ICS_scores = .BDDSegXvol_ICS@Scores

plot_data(BDDSegXvol_ICS_scores, col = "Blues", file = "1.1.4_BDDSegXvol_ICS_scores")

print_table(cor(BDDSegXvol, BDDSegXvol_ICS_scores),
            "1.1.4_BDDSegXvol_ICS_scores_cor")

# Select
BDDSegXvol_ICS_select = BDDSegXvol_ICS_scores[c(1, 2)]
agostino.test(BDDSegXvol_ICS_scores)

plot_data(BDDSegXvol_ICS_select,
          grid = TRUE,
          col = "Blues",
          file = "1.1.4_BDDSegXvol_ICS_select")

# Distances
.BDDSegXvol_ICS_outliers = ics.outlier(.BDDSegXvol_ICS)@ics.distances
BDDSegXvol_ICS_distances = data.frame(cbind(BDDSegX$Date, .BDDSegXvol_ICS_outliers))
colnames(BDDSegXvol_ICS_distances) = c("Date", "ICS distances")
BDDSegXvol_ICS_distances$Date = as.Date(BDDSegX$Date)

plot_data(
  BDDSegXvol_ICS_distances,
  grid = TRUE,
  col = "Blues",
  file = "1.1.4_BDDSegXvol_ICS_distances",
  asp = 0
)
