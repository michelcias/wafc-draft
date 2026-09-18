# E1.4. Conferência numérica da Gram do desenho de produtos
# (derivations/03-desenho-produtos.tex). Roda ANTES da prova; imprime OK ou
# falha com stop(). Notação: docs/notacao.md (congelada em E1.1).
#
# O que confere, com n = 200, J = 3, p = 2 (X_1 = 1), q = 2, Daublets com
# filtro 8, U ~ Unif[0,1]^2 (logo c_U = C_U = 1 e Sigma_Psi = I):
#   1. a base periodizada e a constante formam um sistema ortonormal em [0,1]^q
#      e vale a cota pontual sum_{jk} psi_{jk}(u)^2 <= C_psi 2^J;
#   2. (i)  X indep. de U: Sigma = P (E[XX'] (x) Sigma_Psi) P' exatamente, e
#           lambda_min(Sigma) = lambda_min(E[XX']) lambda_min(Sigma_Psi);
#   3. (ii) X dependente de U com lambda_min(E[XX'|U]) >= kappa_1 > 0:
#           kappa_1 lambda_min(Sigma_Psi) <= lambda_min(Sigma) e
#           lambda_max(Sigma) <= kappa_2 lambda_max(Sigma_Psi);
#   4. contra-exemplo com kappa_1 = 0 (X_2 função de U_1): lambda_min(Sigma)
#      cai com J;
#   5. (iii) versão empírica: posto cheio em n = 200; ||Sigma_hat - Sigma||_op
#           cai com n; lambda_min(Sigma_hat) >= lambda_min(Sigma)/2 em n grande;
#           cota da norma das linhas ||Z_i||^2 <= p B_X^2 (1 + q C_psi 2^J);
#   6. constante de compatibilidade phi^2(S) (cone ||d_{S^c}||_1 <= 3||d_S||_1,
#      S com os p não penalizados) por QP sobre padrões de sinal, e a
#      desigualdade trivial phi^2(S) >= lambda_min(Sigma) usada na derivação.
# A Gram populacional é calculada por quadratura (ponto médio) em [0,1]^2 com
# E[XX' | U = u] em forma fechada, não por Monte Carlo.
#
# Dependências: WaveBased (wbasis, wtable), quadprog (solve.QP).

suppressPackageStartupMessages({
  library(WaveBased)
  library(quadprog)
})
set.seed(20260918)

n  <- 200L
J  <- 3L
p  <- 2L
q  <- 2L
fs <- 8L                       # filter.size (Daublets, 4 momentos nulos)
NJ <- 2L^J - 1L                # wavelets por bloco (j0 = 0, constante descartada)
K  <- 1L + q * NJ              # colunas de Psi(u): constante e q blocos
D  <- p * K                    # colunas de Z: p não penalizadas + p q N_J
tb <- wtable(filter.size = fs)
ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (cond) "ok" else "FALHA", msg))
  if (!cond) ok <<- FALSE
  invisible(cond)
}

## ---- base, desenho e permutação ------------------------------------------

# Psi(U): (1, psi_{jk}(U_1), ..., psi_{jk}(U_q)); a primeira coluna do wbasis()
# com j0 = 0 periódico é a constante, descartada de cada bloco.
psi_block <- function(u, J = J) {
  wbasis(u, j0 = 0, J = J, filter.size = fs, wavelet.table = tb)[, -1, drop = FALSE]
}
Psi <- function(U, J = J) {
  cbind(1, do.call(cbind, lapply(seq_len(ncol(U)), function(m) psi_block(U[, m], J))))
}
# Z = X (x) Psi(U) linha a linha (ordem l-major: bloco de X_l = X_l * Psi(U)),
# depois permutado para a ordem congelada (D12): X_1..X_p, blocos (l, m)
# lexicográficos, dentro do bloco j e k crescentes.
kron_rows <- function(X, P) {
  X[, rep(seq_len(ncol(X)), each = ncol(P)), drop = FALSE] *
    P[, rep(seq_len(ncol(P)), times = ncol(X)), drop = FALSE]
}
perm_D12 <- function(p, q, NJ) {
  K <- 1L + q * NJ
  kron_idx <- function(l, col) as.integer((l - 1L) * K + col)
  c(vapply(seq_len(p), function(l) kron_idx(l, 1L), 1L),
    unlist(lapply(seq_len(p), function(l)
      lapply(seq_len(q), function(m) kron_idx(l, 1L + (m - 1L) * NJ + seq_len(NJ))))))
}
perm <- perm_D12(p, q, NJ)
stopifnot(length(perm) == D, !anyDuplicated(perm))
design <- function(X, U, J = J) kron_rows(X, Psi(U, J))[, perm_D12(ncol(X), ncol(U), 2L^J - 1L), drop = FALSE]

## ---- cenários: E[XX' | U = u] em forma fechada ---------------------------

v_eta <- 1 / 3                  # Var(Unif[-1, 1])
scen <- list(
  # (i) X indep. de U: X_2 = 1/2 + eta; E[XX'] = [[1, 1/2], [1/2, 1/4 + 1/3]]
  A = list(name = "X indep. de U",
           draw_X = function(U) cbind(1, 0.5 + runif(nrow(U), -1, 1)),
           M      = function(u) list(m11 = rep(1, nrow(u)), m12 = rep(0.5, nrow(u)),
                                      m22 = rep(0.25 + v_eta, nrow(u)))),
  # (ii) X dependente de U: X_2 = (U_1 - 1/2) + eta; kappa_1 > 0 porque Var(eta) > 0
  B = list(name = "X dependente de U, kappa_1 > 0",
           draw_X = function(U) cbind(1, (U[, 1] - 0.5) + runif(nrow(U), -1, 1)),
           M      = function(u) list(m11 = rep(1, nrow(u)), m12 = u[, 1] - 0.5,
                                      m22 = (u[, 1] - 0.5)^2 + v_eta)),
  # contra-exemplo: X_2 função de U_1, E[XX' | U] de posto 1, kappa_1 = 0
  C = list(name = "X_2 = U_1 - 1/2 (kappa_1 = 0)",
           draw_X = function(U) cbind(1, U[, 1] - 0.5),
           M      = function(u) list(m11 = rep(1, nrow(u)), m12 = u[, 1] - 0.5,
                                      m22 = (u[, 1] - 0.5)^2))
)
B_X <- 1.5                      # sup |X_l| nos cenários A e B

# Sigma = E[ E[XX'|U] (x) Psi(U) Psi(U)' ] por quadratura no grid de ponto
# médio G x G em [0,1]^2 (U uniforme), montada bloco a bloco e permutada.
G <- 400L
g1 <- (seq_len(G) - 0.5) / G
Ug <- as.matrix(expand.grid(u1 = g1, u2 = g1))
pop_gram <- function(sc, J = J) {
  P  <- Psi(Ug, J)
  Kj <- ncol(P)
  M  <- sc[["M"]](Ug)
  S  <- matrix(0, p * Kj, p * Kj)
  blk <- function(w) crossprod(P, w * P) / nrow(P)
  S[1:Kj, 1:Kj]                 <- blk(M[["m11"]])
  S[1:Kj, Kj + 1:Kj]            <- blk(M[["m12"]])
  S[Kj + 1:Kj, 1:Kj]            <- t(S[1:Kj, Kj + 1:Kj])
  S[Kj + 1:Kj, Kj + 1:Kj]       <- blk(M[["m22"]])
  pm <- perm_D12(p, q, 2L^J - 1L)
  S[pm, pm]
}
lam <- function(S) { e <- eigen(S, symmetric = TRUE, only.values = TRUE)[["values"]]; c(min = min(e), max = max(e)) }
opnorm <- function(A) max(abs(eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)[["values"]]))

## ---- 1. ortonormalidade e cota pontual -----------------------------------

cat("1. Sistema {1, psi_jk(u_m)} em [0,1]^q e cota pontual\n")
Pg <- Psi(Ug, J)
Sigma_Psi <- crossprod(Pg) / nrow(Pg)
dev_orth <- max(abs(Sigma_Psi - diag(K)))
chk(dev_orth < 1e-4, sprintf("Gram de Psi no grid = I (desvio máximo %.1e, erro da quadratura com G = %d)", dev_orth, G))
lam_Psi <- lam(Sigma_Psi)
gf <- (seq_len(2^14) - 0.5) / 2^14
C_psi <- vapply(2:6, function(Jj) max(rowSums(psi_block(gf, Jj)^2)) / 2^Jj, 1)
names(C_psi) <- paste0("J", 2:6)
cat("  C_psi(J) = max_u sum_jk psi_jk(u)^2 / 2^J:", paste(sprintf("%.3f", C_psi), collapse = " "), "\n")
chk(all(diff(C_psi) >= 0) && diff(range(C_psi[3:5])) < 0.1,
    "C_psi cresce e estabiliza em J (série geométrica em 2^{-j})")
Cpsi <- C_psi[["J3"]]

## ---- 2. (i) fatoração sob independência ----------------------------------

cat("2. (i) X indep. de U: Sigma = P (E[XX'] (x) Sigma_Psi) P'\n")
S_A  <- pop_gram(scen[["A"]], J)
EXX  <- matrix(c(1, 0.5, 0.5, 0.25 + v_eta), 2, 2)
S_Ak <- (EXX %x% Sigma_Psi)[perm, perm]
dev_kron <- max(abs(S_A - S_Ak))
chk(dev_kron < 1e-10, sprintf("Sigma coincide com o Kronecker permutado (desvio %.1e)", dev_kron))
lA <- lam(S_A); lX <- lam(EXX)
chk(abs(lA[["min"]] - lX[["min"]] * lam_Psi[["min"]]) < 1e-6 &&
    abs(lA[["max"]] - lX[["max"]] * lam_Psi[["max"]]) < 1e-6,
    sprintf("lambda_min(Sigma) = %.4f = %.4f x %.4f; lambda_max = %.4f = %.4f x %.4f",
            lA[["min"]], lX[["min"]], lam_Psi[["min"]], lA[["max"]], lX[["max"]], lam_Psi[["max"]]))

## ---- 3. (ii) caso geral sob lambda_min(E[XX'|U]) >= kappa_1 --------------

cat("3. (ii) X dependente de U: kappa_1 lambda_min(Sigma_Psi) <= lambda_min(Sigma)\n")
S_B <- pop_gram(scen[["B"]], J)
MB  <- scen[["B"]][["M"]](Ug)
lamM <- function(M) {                 # autovalores de [[m11, m12], [m12, m22]] ponto a ponto
  tr <- M[["m11"]] + M[["m22"]]; dt <- M[["m11"]] * M[["m22"]] - M[["m12"]]^2
  cbind(min = (tr - sqrt(pmax(tr^2 - 4 * dt, 0))) / 2, max = (tr + sqrt(pmax(tr^2 - 4 * dt, 0))) / 2)
}
lM <- lamM(MB)
kappa_1 <- min(lM[, "min"]); kappa_2 <- max(lM[, "max"])
lB <- lam(S_B)
chk(lB[["min"]] >= kappa_1 * lam_Psi[["min"]] - 1e-10,
    sprintf("lambda_min(Sigma) = %.4f >= kappa_1 lambda_min(Sigma_Psi) = %.4f x %.4f = %.4f",
            lB[["min"]], kappa_1, lam_Psi[["min"]], kappa_1 * lam_Psi[["min"]]))
chk(lB[["max"]] <= kappa_2 * lam_Psi[["max"]] + 1e-10,
    sprintf("lambda_max(Sigma) = %.4f <= kappa_2 lambda_max(Sigma_Psi) = %.4f x %.4f = %.4f",
            lB[["max"]], kappa_2, lam_Psi[["max"]], kappa_2 * lam_Psi[["max"]]))
# E[XX'] = E[E[XX'|U]] tem lambda_min >= kappa_1 (concavidade de lambda_min)
EXX_B <- matrix(c(mean(MB[["m11"]]), mean(MB[["m12"]]), mean(MB[["m12"]]), mean(MB[["m22"]])), 2, 2)
chk(lam(EXX_B)[["min"]] >= kappa_1 - 1e-10,
    sprintf("lambda_min(E[XX']) = %.4f >= kappa_1 = %.4f (a hipótese condicional implica a marginal)",
            lam(EXX_B)[["min"]], kappa_1))

## ---- 4. contra-exemplo: kappa_1 = 0 ---------------------------------------

cat("4. Contra-exemplo X_2 = U_1 - 1/2 (E[XX'|U] singular): lambda_min(Sigma) por J\n")
lC <- vapply(2:4, function(Jj) lam(pop_gram(scen[["C"]], Jj))[["min"]], 1)
names(lC) <- paste0("J", 2:4)
cat("  lambda_min(Sigma):", paste(sprintf("%s = %.4f", names(lC), lC), collapse = ", "), "\n")
chk(all(diff(lC) < 0) && lC[["J3"]] < 0.5 * lB[["min"]],
    "cai com J e, em J = 3, fica abaixo da metade do caso (ii)")

## ---- 5. (iii) versão empírica ---------------------------------------------

cat("5. (iii) Gram empírica: posto, concentração em n e cota das linhas\n")
R <- 200L
emp <- function(sc, S_pop, n) {
  l_pop <- lam(S_pop)[["min"]]
  out <- t(replicate(R, {
    U <- matrix(runif(n * q), n, q)
    X <- sc[["draw_X"]](U)
    Z <- design(X, U, J)
    Sh <- crossprod(Z) / n
    c(lmin = lam(Sh)[["min"]], dev = opnorm(Sh - S_pop), rown = max(rowSums(Z^2)))
  }))
  c(frac_half = mean(out[, "lmin"] >= l_pop / 2), min_lmin = min(out[, "lmin"]),
    med_ratio = median(out[, "lmin"]) / l_pop, med_dev = median(out[, "dev"]),
    max_rown = max(out[, "rown"]))
}
row_bound <- p * B_X^2 * (1 + q * Cpsi * 2^J)
for (nm in c("A", "B")) {
  S_pop <- if (nm == "A") S_A else S_B
  e200  <- emp(scen[[nm]], S_pop, 200L)
  e3200 <- emp(scen[[nm]], S_pop, 3200L)
  cat(sprintf("  cenário %s (%s), D = %d colunas:\n", nm, scen[[nm]][["name"]], D))
  cat(sprintf("    n = 200 : min lambda_min(Sigma_hat) = %.4f; mediana lambda_min(Sigma_hat)/lambda_min(Sigma) = %.3f; P(>= 1/2) = %.3f; mediana ||Sigma_hat - Sigma||_op = %.3f\n",
              e200[["min_lmin"]], e200[["med_ratio"]], e200[["frac_half"]], e200[["med_dev"]]))
  cat(sprintf("    n = 3200: min lambda_min(Sigma_hat) = %.4f; mediana lambda_min(Sigma_hat)/lambda_min(Sigma) = %.3f; P(>= 1/2) = %.3f; mediana ||Sigma_hat - Sigma||_op = %.3f\n",
              e3200[["min_lmin"]], e3200[["med_ratio"]], e3200[["frac_half"]], e3200[["med_dev"]]))
  chk(e200[["min_lmin"]] > 1e-3, sprintf("posto cheio em todas as %d réplicas com n = 200", R))
  chk(e3200[["frac_half"]] >= 0.95, "lambda_min(Sigma_hat) >= lambda_min(Sigma)/2 em >= 95% das réplicas com n = 3200")
  chk(e200[["med_dev"]] / e3200[["med_dev"]] > 2.5,
      sprintf("||Sigma_hat - Sigma||_op cai de n = 200 para 3200 por fator %.2f (1/sqrt(n) daria 4)",
              e200[["med_dev"]] / e3200[["med_dev"]]))
  chk(max(e200[["max_rown"]], e3200[["max_rown"]]) <= row_bound,
      sprintf("max_i ||Z_i||^2 = %.1f <= p B_X^2 (1 + q C_psi 2^J) = %.1f",
              max(e200[["max_rown"]], e3200[["max_rown"]]), row_bound))
}

## ---- 6. compatibilidade no cone por QP ------------------------------------

cat("6. Constante de compatibilidade phi^2(S) no cone ||d_Sc||_1 <= 3 ||d_S||_1\n")
# phi^2(S) = s_0 min { d' Sigma d : ||d_S||_1 = 1, ||d_Sc||_1 <= 3 }. Para cada
# padrão de sinal sigma de d_S o conjunto {sigma_i d_i >= 0, sigma' d_S = 1} é
# uma face convexa; d_Sc = dplus - dminus com dplus, dminus >= 0 e
# sum(dplus + dminus) <= 3. Um QP por padrão (2^{s_0} QPs), o mínimo é exato.
compat <- function(S, Sset) {
  Dd <- ncol(S); s0 <- length(Sset); Sc <- setdiff(seq_len(Dd), Sset); r <- length(Sc)
  A <- matrix(0, Dd, s0 + 2 * r)                       # d = A z, z = (d_S, dplus, dminus)
  A[cbind(Sset, seq_len(s0))] <- 1
  A[cbind(Sc, s0 + seq_len(r))] <- 1
  A[cbind(Sc, s0 + r + seq_len(r))] <- -1
  Dmat <- crossprod(A, S %*% A); Dmat <- (Dmat + t(Dmat)) / 2 + 1e-9 * diag(ncol(A))
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), s0)))
  best <- Inf
  for (i in seq_len(nrow(signs))) {
    sg <- signs[i, ]
    # restrições: sg' d_S = 1 (igualdade); sg_i d_i >= 0; dplus, dminus >= 0; 3 - sum(dplus + dminus) >= 0
    Amat <- cbind(c(sg, rep(0, 2 * r)),
                  rbind(diag(sg, s0), matrix(0, 2 * r, s0)),
                  rbind(matrix(0, s0, 2 * r), diag(2 * r)),
                  c(rep(0, s0), rep(-1, 2 * r)))
    bvec <- c(1, rep(0, s0), rep(0, 2 * r), -3)
    sol <- solve.QP(Dmat, rep(0, ncol(A)), Amat, bvec, meq = 1)
    best <- min(best, sol[["value"]] * 2)               # solve.QP minimiza (1/2) z' D z
  }
  s0 * best
}
s_pen <- 5L
Sset <- c(seq_len(p), p + sort(sample(D - p, s_pen)))    # os p não penalizados sempre em S
s0 <- length(Sset)
R2 <- 20L
for (nm in c("A", "B")) {
  S_pop <- if (nm == "A") S_A else S_B
  ph_pop <- compat(S_pop, Sset)
  emp2 <- t(replicate(R2, {
    U <- matrix(runif(n * q), n, q); X <- scen[[nm]][["draw_X"]](U)
    Sh <- crossprod(design(X, U, J)) / n
    c(phi = compat(Sh, Sset), lmin = lam(Sh)[["min"]])
  }))
  cat(sprintf("  cenário %s, |S| = %d: phi^2(S) populacional = %.4f (lambda_min = %.4f); empírica n = 200 em %d réplicas: mediana %.4f, mínimo %.4f (mediana lambda_min = %.4f)\n",
              nm, s0, ph_pop, lam(S_pop)[["min"]], R2, median(emp2[, "phi"]), min(emp2[, "phi"]), median(emp2[, "lmin"])))
  # phi^2(S) = s_0 min d'Sigma d sobre ||d_S||_1 = 1 >= s_0 lambda_min ||d_S||_2^2 >= lambda_min
  chk(ph_pop >= lam(S_pop)[["min"]] - 1e-6 && all(emp2[, "phi"] >= emp2[, "lmin"] - 1e-6),
      "phi^2(S) >= lambda_min(Sigma), populacional e empírica (a cota usada na derivação)")
  chk(median(emp2[, "phi"]) >= 0.5 * ph_pop,
      "mediana empírica de phi^2(S) em n = 200 é pelo menos metade da populacional")
}

cat("\n")
if (ok) cat("OK\n") else stop("E1.4: conferência numérica FALHOU (ver linhas acima)")
