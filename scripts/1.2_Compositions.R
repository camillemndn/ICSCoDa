source("%MYFOLDER%/scripts/Functions.R")

## ---- 1.2.1_simplex_operations

# Simplex
pdf(
  file = paste0(.figures, "1.2.1_simplex.pdf"),
  width = 2.5,
  height = 2.5
)
plot_data.acomp(matrix(1, 1, 3))
dev.off()

# Perturbations
t21 = acomp(read.csv(paste0(.data, "data.txt")))
# dist(clo(t21[1:2,]))

a = acomp(c(0.7, 0.4, 0.8))
# b = acomp(c(0.2, 0.8, 0.1))
# e1 = acomp(c(exp(1), 1, 1))
# e2 = acomp(c(1, exp(1), 1))

plot(t21)
plot(t21 + a)

# a %*% b
