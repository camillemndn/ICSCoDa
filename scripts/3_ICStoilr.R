source("%MYFOLDER%/scripts/Functions.R")

## ---- 3.2_mixturecomp

# IC
.mixturecomp = rmvnormmix.acomp(150,
                                c(0.92, 0.08),
                                matrix(c(rep(1, 5), 4, 1, 1, 1, 1), 2, 5, TRUE),
                                0.1 * identity.acomp(5))
.mixturecomp_index = .mixturecomp$index

.B = matrix(rnorm(25), 5, 5)
#.B = identity.acomp(5)
mixturecomp = clrInv(clr(.mixturecomp$data) %*% .B)
.mixturecomp_ICS = ics.acomp(mixturecomp)
mixturecomp_ICS_scores = .mixturecomp_ICS$data

plot_data(
  mixturecomp,
  col = .mixturecomp_index,
  width = 14,
  height = 10,
  file = "3.2_mixturecomp"
)

plot_data(
  mixturecomp_ICS_scores,
  col = .mixturecomp_index,
  width = 14,
  height = 10,
  file = "3.2_mixturecomp_ICS_scores"
)

print_table(.mixturecomp_ICS$cov_log, file = "3.2_mixturecomp_ICS_scores_cov_log")

# Select
.mixturecomp_ICS_outliers = norm.acomp(groupparts(.mixturecomp_ICS$data, 1, 2, 3:5))
#.mixturecomp_ICS_outliers = apply(matrix(clr(.mixturecomp_ICS$data)[, 2]), 1, function(x) {
#norm(matrix(x), "2")
#})
plot_data(groupparts(.mixturecomp_ICS$data, 1, 2, 3:5),
          col = .mixturecomp_index,
          file = "3.2_mixturecomp_ICS_select")


# Distances
mixturecomp_ICS_distances = data.frame(cbind(
  1:length(.mixturecomp_ICS_outliers),
  .mixturecomp_ICS_outliers
))
colnames(mixturecomp_ICS_distances) = c("Index", "ICS distances")

plot_data(
  mixturecomp_ICS_distances,
  grid = TRUE,
  col = .mixturecomp_index,
  asp = 0,
  file = "3.2_mixturecomp_ICS_distances"
)


## ---- 3.3_ICS_indep_ilr

.D = 5
.V1 = ilrBase(D = .D)
.V2 = t(matrix(ilr2clr(randortho(.D - 1)), .D - 1, .D))

plot_data(
  ics.acomp(mixturecomp, V = .V1)$data,
  col = .mixturecomp_index,
  width = 14,
  height = 10,
  file = "3.3_mixturecomp_ICS_scores_V1"
)

plot_data(
  ics.acomp(mixturecomp, V = .V2)$data,
  col = .mixturecomp_index,
  width = 14,
  height = 10,
  file = "3.3_mixturecomp_ICS_scores_V2"
)

# Non, à cause de la non-unicité des vecteurs propres et de l'algorithme "eigen"


## ---- 3.4_BDDSegXcomp_ICS


# IC
.BDDSegXcomp_ICS = ics.acomp(BDDSegXcomp)
BDDSegXcomp_ICS_scores = .BDDSegXcomp_ICS$data

plot_data(
  BDDSegXcomp_ICS_scores,
  col = "Blues",
  width = 14,
  height = 10,
  file = "3.4_BDDSegXcomp_ICS_scores"
)

print_table(.BDDSegXcomp_ICS$cov_log, file = "3.4_BDDSegXcomp_ICS_scores_cov_log")

# Select
.BDDSegXcomp_ICS_outliers = norm.acomp(groupparts(.BDDSegXcomp_ICS$data, 4, c(1:3, 5)))
#.BDDSegXcomp_ICS_outliers = apply(matrix(clr(.BDDSegXcomp_ICS$data)[, 2]), 1, function(x) {
#norm(matrix(x), "2")
#})
plot_data(groupparts(.BDDSegXcomp_ICS$data, 4, 5, 1:3),
          col = "Blues",
          file = "3.4_BDDSegXcomp_ICS_select")


# Distances
BDDSegXcomp_ICS_distances = data.frame(cbind(
  BDDSegX$Date,
  .BDDSegXcomp_ICS_outliers
))
colnames(BDDSegXcomp_ICS_distances) = c("Date", "ICS distances")
BDDSegXcomp_ICS_distances$Date = BDDSegX$Date

plot_data(
  BDDSegXcomp_ICS_distances,
  grid = TRUE,
  col = "Blues",
  asp = 0,
  file = "3.4_BDDSegXcomp_ICS_distances"
)
