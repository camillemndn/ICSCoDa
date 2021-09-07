source("%MYFOLDER%/scripts/Functions.R")

## ---- 2.2_BDDSegXcomp_whitening

# BDDSegXcomp
pdf(
  file = paste0(.figures, "2.0_BDDSegXcomp.pdf"),
  width = 14,
  height = 10
)
plot_data.acomp(BDDSegXcomp, grid = FALSE, col = TRUE)
dev.off()

# Whitening
doubleplot.acomp(BDDSegXcomp[c(1, 2, 3)],
                 "2.2_BDDSegXcomp_whitening",
                 width = 6,
                 height = 3)


## ---- 2.3_gaussiancomp

gaussiancomp = rmvnormmix.acomp(100, mu= c(1, 6, 4), sigma = 0.1 * identity.acomp(3))$data
doubleplot.acomp(gaussiancomp, "2.3_gaussiancomp", width = 7, height = 4, grid = TRUE)


## ---- student

stud_lin = rmvt(100, c(1, 0), diag(c(0.05, 0.05)), 1)
outlierplot(ilrInv(stud_lin))
outlierplot(Aar)


## ---- dirichlet


## ---- simplex_matrices

V = ilrBase(D = 3)
I_s = diag(c(1, 1))
G = V %*% I_s %*% t(V)
one = diag(c(1, 1, 1)) - G
A_s = rbind(c(1, 1), c(1, 0))
A = V %*% A_s %*% t(V)
B_s = inv(A_s)
B = V %*% B_s %*% t(V)
D = A %*% B
C = V %*% A_s %*% B_s %*% t(V)

Sigma_s = diag(c(1, 2))
Sigma = V %*% Sigma_s %*% t(V)
e_1 = acomp(clrInv(V[, 1]))
e_2 = acomp(clrInv(V[, 2]))
eigen(Sigma)
