# derivations/check/01-identificabilidade.R
#
# Conferência numérica de E1.2 (Proposição 1 e Lema 1 de
# derivations/01-identificabilidade.md). Forma densa, n = 200, J = 3 no caso
# periódico; o caso boundary = "interval" exige j0 >= 4 com filter.size = 8 e
# usa J = 5. Imprime OK no fim ou para no primeiro stopifnot() que falhar.
#
# O que se confere:
#   1. Lema 1 (base): no caso periódico com j0 = 0, a única função de escala
#      é a constante e cada psi_{jk} tem integral de Lebesgue zero.
#   2. Posto cheio da matriz de desenho Z (ordem de colunas D12) em n = 200,
#      com U dependente entre coordenadas e marginais não uniformes.
#   3. Recuperação exata de (c*, theta*) a partir de y = sum_l beta_l(U) X_l,
#      sem ruído, por mínimos quadrados.
#   4. O que as hipóteses excluem, cada caso com a deficiência de posto
#      prevista: manter phi_00 (deficiência p q); moduladora repetida
#      (deficiência p N_J); X_l função de U dentro do sieve (deficiência >= 1).
#   5. boundary = "interval": a constante está no espaço das funções de
#      escala e não é uma coluna isolada; manter as 2^{j0} funções de escala
#      dá deficiência p q; a reparametrização Phi Q, com Q base ortonormal de
#      {alpha : mu' alpha = 0} e mu_k = int phi_{j0 k}, restaura o posto cheio
#      e a recuperação exata.
#
# Regra do projeto (instrucoes.md, §5): elemento de lista com [[ ]] e nome
# completo, nunca com $.

suppressMessages(library(WaveBased))
set.seed(20260918)

n  <- 200L
p  <- 2L      # X_1 == 1, X_2 ~ N(0, 1)
q  <- 2L      # U_1 ~ Unif(0,1), U_2 ~ Beta(2,3), copula gaussiana rho = 0.6
J  <- 3L
L  <- 8L      # Daublets, filter.size = 8 (N = 4 momentos nulos)
tol_rank <- 1e-8   # razão sigma_min/sigma_max abaixo disto conta como singular

rank_of <- function(Z, tol = tol_rank) {
  d <- svd(Z, nu = 0, nv = 0)[["d"]]
  sum(d > tol * d[1])
}

# ---- bases -------------------------------------------------------------------

# Bloco periódico, j0 = 0: wbasis() devolve [phi_00 | psi_00 | psi_10 psi_11 | ...].
# A coluna phi_00 é a constante; Lema 1 descarta-a.
per_basis <- function(u, J, drop_phi = TRUE) {
  W <- wbasis(u, j0 = 0, J = J, family = "Daublets", filter.size = L,
              boundary = "periodic")
  colnames(W) <- c("phi0.0", jk_names(0L, J))
  if (drop_phi) W[, -1, drop = FALSE] else W
}

jk_names <- function(j0, J) {
  unlist(lapply(j0:(J - 1L), function(j) paste0("psi", j, ".", 0:(2^j - 1L))))
}

# Matriz de desenho Z na ordem D12: os p termos não penalizados, depois os
# blocos (l, m) em ordem lexicográfica, dentro do bloco j e k crescentes.
design <- function(X, U, block_fun) {
  p <- ncol(X); q <- ncol(U)
  blocks <- list(X)
  nm <- paste0("X", seq_len(p))
  for (l in seq_len(p)) for (m in seq_len(q)) {
    B <- block_fun(U[, m])
    blocks[[length(blocks) + 1L]] <- X[, l] * B
    nm <- c(nm, paste0("X", l, ".U", m, ".", colnames(B)))
  }
  Z <- do.call(cbind, blocks)
  colnames(Z) <- nm
  Z
}

# ---- 1. Lema 1: a base periódica impõe a centralização ----------------------

grid <- (seq_len(2^14) - 0.5) / 2^14                 # regra do ponto médio
Wg <- per_basis(grid, J, drop_phi = FALSE)
NJ <- 2^J - 1L
stopifnot(ncol(Wg) == 2^J)
stopifnot(max(abs(Wg[, 1] - 1)) < 1e-12)              # phi_00 == 1
leb_mean <- colMeans(Wg[, -1, drop = FALSE])           # int psi_{jk} du
stopifnot(max(abs(leb_mean)) < 1e-6)
gram_err <- max(abs(crossprod(Wg) / length(grid) - diag(2^J)))
stopifnot(gram_err < 1e-6)                             # {1} U {psi_jk} ortonormal
cat(sprintf("1. periódico j0=0, J=%d: phi_00 == 1; max |int psi_jk| = %.1e; max |Gram - I| = %.1e\n",
            J, max(abs(leb_mean)), gram_err))

# ---- dados -------------------------------------------------------------------

rho <- 0.6
G <- matrix(rnorm(2L * n), n, 2L)
G[, 2] <- rho * G[, 1] + sqrt(1 - rho^2) * G[, 2]
U <- cbind(pnorm(G[, 1]), qbeta(pnorm(G[, 2]), 2, 3))  # dependentes, marginais distintas
X <- cbind(1, rnorm(n))                                 # X_1 == 1 permitido (D1)
stopifnot(all(U > 0 & U < 1))

# ---- 2. posto cheio em n = 200 -------------------------------------------------

Z <- design(X, U, function(u) per_basis(u, J))
d_full <- p + p * q * NJ
stopifnot(ncol(Z) == d_full)
sv <- svd(Z, nu = 0, nv = 0)[["d"]]
stopifnot(rank_of(Z) == d_full)
lmin <- min(eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]])
cat(sprintf("2. Z: n=%d, %d colunas, posto %d; sigma_min/sigma_max = %.2e; lambda_min(Z'Z/n) = %.3e\n",
            n, ncol(Z), rank_of(Z), min(sv) / max(sv), lmin))

# A centralização da base é de Lebesgue, não de P_U: com U_2 ~ Beta(2,3),
# mean(psi_jk(U_2)) fica longe de zero (só informação; ver o handoff).
pu_mean <- colMeans(per_basis(U[, 2], J))
cat(sprintf("   (info) max |mean psi_jk(U_2)| com U_2 ~ Beta(2,3): %.3f (Lebesgue: %.1e)\n",
            max(abs(pu_mean)), max(abs(leb_mean))))

# ---- 3. recuperação exata de (c*, theta*) -----------------------------------------

c_star <- c(1.5, -0.7)
theta_star <- array(rnorm(p * q * NJ), dim = c(NJ, q, p))   # theta[jk, m, l]
g <- function(l, m, u) drop(per_basis(u, J) %*% theta_star[, m, l])
beta_fun <- function(l, Umat) {
  c_star[l] + Reduce(`+`, lapply(seq_len(q), function(m) g(l, m, Umat[, m])))
}
y <- Reduce(`+`, lapply(seq_len(p), function(l) beta_fun(l, U) * X[, l]))
est <- qr.solve(Z, y)
truth <- c(c_star, as.vector(theta_star))   # ordem D12: (l, m) lexicográfico, jk dentro
rec_err <- max(abs(est - truth))
stopifnot(rec_err < 1e-8)
cat(sprintf("3. recuperação de (c*, theta*) sem ruído: max |erro| = %.1e\n", rec_err))

# ---- 4. o que as hipóteses excluem -----------------------------------------------

# (a) manter phi_00: a coluna X_l phi_00(U_m) = X_l repete a não penalizada.
Za <- design(X, U, function(u) per_basis(u, J, drop_phi = FALSE))
def_a <- ncol(Za) - rank_of(Za)
stopifnot(def_a == p * q)
# (b) moduladora repetida (U_2 = U_1): os blocos (l,1) e (l,2) coincidem.
Zb <- design(X, cbind(U[, 1], U[, 1]), function(u) per_basis(u, J))
def_b <- ncol(Zb) - rank_of(Zb)
stopifnot(def_b == p * NJ)
# (c) X_2 função de U_1 dentro do sieve (E[XX'|U] singular): X_2 = 1 + psi_10(U_1).
Xc <- cbind(1, 1 + per_basis(U[, 1], J)[, "psi1.0"])
Zc <- design(Xc, U, function(u) per_basis(u, J))
def_c <- ncol(Zc) - rank_of(Zc)
stopifnot(def_c >= 1L)
cat(sprintf("4. deficiências de posto: phi_00 mantida %d (= pq); U_2 = U_1 %d (= p N_J); X_2 = 1 + psi_10(U_1) %d\n",
            def_a, def_b, def_c))

# ---- 5. boundary = "interval" ------------------------------------------------------

j0i <- 4L; Ji <- 5L    # com filter.size = 8 o wbasis() exige j0 >= 4
int_basis <- function(u, j0 = j0i, J = Ji) {
  W <- wbasis(u, j0 = j0, J = J, family = "Daublets", filter.size = L,
              boundary = "interval")
  colnames(W) <- c(paste0("phi", j0, ".", 0:(2^j0 - 1L)), jk_names(j0, J))
  W
}
Wi <- int_basis(grid)
stopifnot(ncol(Wi) == 2^Ji)
# As funções de borda são menos regulares e a quadratura no grid converge em
# O(h^2): 6e-6 com 2^14 pontos, 4e-7 com 2^16. A tolerância é a da quadratura.
gram_i <- max(abs(crossprod(Wi) / length(grid) - diag(2^Ji)))
stopifnot(gram_i < 1e-4)
Phi <- Wi[, seq_len(2^j0i)]
Psi <- Wi[, -seq_len(2^j0i)]
# a constante pertence a V_{j0}: 1 = sum_k mu_k phi_{j0 k}, mu_k = int phi_{j0 k}
mu <- qr.solve(Phi, rep(1, length(grid)))
res_const <- max(abs(Phi %*% mu - 1))
stopifnot(res_const < 1e-6)
stopifnot(max(abs(mu - colMeans(Phi))) < 1e-5)
stopifnot(max(abs(colMeans(Psi))) < 1e-5)             # wavelets de intervalo: integral zero
stopifnot(max(abs(colMeans(Phi))) > 0.1)              # funções de escala: não
cat(sprintf("5. intervalo j0=%d, J=%d: max |Gram - I| = %.1e; 1 em V_j0 (resíduo %.1e); max |int psi| = %.1e; max |int phi| = %.2f\n",
            j0i, Ji, gram_i, res_const, max(abs(colMeans(Psi))), max(abs(colMeans(Phi)))))

# (a) manter as 2^{j0} funções de escala por bloco: deficiência p q. A
# identidade 1 = sum_k mu_k phi_{j0 k} vale a menos da precisão da avaliação
# CDV (~3e-9), logo essas p q direções aparecem em sigma ~ 1e-9 sigma_max, não
# em 1e-15; tol_rank = 1e-8 as conta como nulas.
#
# Com o U populacional (U_2 ~ Beta(2,3), max U_2 < 1 - 1/16 em n = 200) a
# última célula diádica fica vazia e as 2 x 4 funções de borda direita do
# bloco U_2 (phi_{4,12..15}, psi_{4,12..15}) ficam linearmente dependentes na
# amostra: deficiência extra, de amostra finita, não populacional. O wall()
# evita isso com boundary = "interval" reescalando cada U_m para [0, 1] pela
# amplitude amostral (eps = 0), que é o que se faz abaixo.
Zi_raw <- design(X, U, int_basis)
def_raw <- ncol(Zi_raw) - rank_of(Zi_raw)
stopifnot(def_raw >= p * q)
U01 <- apply(U, 2L, function(u) (u - min(u)) / (max(u) - min(u)))
Zi_full <- design(X, U01, int_basis)
sv_i <- svd(Zi_full, nu = 0, nv = 0)[["d"]]
def_i <- ncol(Zi_full) - rank_of(Zi_full)
stopifnot(ncol(Zi_full) == p + p * q * 2^Ji)
stopifnot(def_i == p * q)
cat(sprintf("   escala mantida (%d colunas): deficiência %d com U populacional (max U_2 = %.3f), %d (= pq) com U reescalado à amplitude amostral; sigma/sigma_max nas 6 últimas: %s\n",
            ncol(Zi_full), def_raw, max(U[, 2]), def_i, paste(signif(tail(sv_i, 6) / sv_i[1], 2), collapse = " ")))

# (b) reparametrização: Phi Q com Q base ortonormal do complemento de mu.
Q <- qr.Q(qr(mu), complete = TRUE)[, -1, drop = FALSE]    # 2^{j0} x (2^{j0} - 1)
stopifnot(max(abs(crossprod(Q, mu))) < 1e-12)
stopifnot(max(abs(colMeans(Phi %*% Q))) < 1e-5)           # as novas colunas têm integral zero
int_basis_c <- function(u) {
  W <- int_basis(u)
  out <- cbind(W[, seq_len(2^j0i)] %*% Q, W[, -seq_len(2^j0i)])
  colnames(out) <- c(paste0("phiQ", j0i, ".", seq_len(ncol(Q))), colnames(W)[-seq_len(2^j0i)])
  out
}
Zi <- design(X, U01, int_basis_c)
NJi <- 2^Ji - 1L
stopifnot(ncol(Zi) == p + p * q * NJi)
stopifnot(rank_of(Zi) == ncol(Zi))
lmin_i <- min(eigen(crossprod(Zi) / n, symmetric = TRUE, only.values = TRUE)[["values"]])
# recuperação exata com a verdade no sieve reparametrizado
coef_i <- array(rnorm(p * q * NJi), dim = c(NJi, q, p))
gi <- function(l, m, u) drop(int_basis_c(u) %*% coef_i[, m, l])
beta_i <- function(l, Umat) {
  c_star[l] + Reduce(`+`, lapply(seq_len(q), function(m) gi(l, m, Umat[, m])))
}
yi <- Reduce(`+`, lapply(seq_len(p), function(l) beta_i(l, U01) * X[, l]))
rec_i <- max(abs(qr.solve(Zi, yi) - c(c_star, as.vector(coef_i))))
stopifnot(rec_i < 1e-8)
cat(sprintf("   Phi Q: %d colunas, posto cheio, lambda_min(Z'Z/n) = %.3e, recuperação max |erro| = %.1e\n",
            ncol(Zi), lmin_i, rec_i))

cat("OK\n")
