# E1.9. Conferência numérica da cota inferior minimax no nível da taxa
# (derivations/09-cota-inferior.tex: Lema 17, Teorema 5, Corolário 16). Roda
# ANTES da prova; imprime OK ou falha com stop(). Notação: docs/notacao.md
# (congelada em E1.1); símbolos locais do 09: o conjunto de blocos K, o
# conjunto de coordenadas C (blocos de K no nível j), M = |C| = |K| 2^j, o
# hipercubo Omega = {0,1}^C, a altura zeta, a distância de Hamming d_H e
# gamma = kappa_1 c_U.
#
# Desenho (X dependente de U, q = 2, marginal de U_1 não uniforme): Daublets
# com filtro 8; p = 2 com X_1 = 1 e X_2 = h(U) + V, h(u) = 1/2 + sin(2 pi u_1)
# u_2 / 2 em [0, 1], V ~ Unif(-1, 1) independente de U; U com densidade
# conjunta p_U(u) = 1 + 0.6 (u_1 - 1/2) + 0.3 cos(2 pi (u_1 + u_2)) em
# [0,1]^2, logo c_U = 0.4, C_U = 1.6, B_X = 2, e E(XX' | U = u) =
# [[1, h], [h, h^2 + 1/3]] em forma fechada. As integrais em (X, U) são
# quadratura de ponto médio em [0,1]^2; Monte Carlo só onde se confere uma
# probabilidade (o log da razão de verossimilhança, o risco de Bayes, o evento
# da Gram empírica).
#
# Partes:
#   I    a base: ortonormalidade de Lebesgue e média zero no nível j;
#   II   as constantes do desenho: c_U, C_U, kappa_1, kappa_2, B_X;
#   III  o Kullback-Leibler (Lema 17(ii)): E(X_l^2 psi_jk(U_m)^2) <=
#        kappa_2 C_U <= B_X^2 C_U para todo (l, m, j, k); e a fórmula
#        KL = n E(Delta^2) / (2 sigma^2) com X dependente de U, por Monte
#        Carlo do log da razão de verossimilhança (vizinhos e um par qualquer);
#   IV   a separação (Lema 17(iii)): lambda_min da Gram do hipercubo em L_2(P)
#        >= kappa_1 c_U para todo K, com os termos cruzados entre moduladoras e
#        entre covariáveis; a distância das componentes é zeta^2 d_H;
#   V    a redução (Lema 17(iv)): para estimadores arbitrários, fora do
#        hipercubo e não aditivos, L_comp >= zeta^2 d_H(w_hat, w) / 4,
#        L_P >= gamma zeta^2 d_H / 4 e, no evento lambda_min(Sigma_hat_CC) >=
#        gamma / 2, L_n >= gamma zeta^2 d_H / 8;
#   VI   Assouad: o risco de Bayes de Hamming sob a priori uniforme (com o
#        estimador de Bayes exato, a maioria das marginais a posteriori) fica
#        acima de (M/2)(1 - sqrt(alpha/2)) com o KL do enunciado;
#   VII  o pertencimento (Lema 17(i)): a norma de Besov em sequência de
#        g_omega, calculada da função por quadratura, é zeta 2^{j(s+1/2)} <=
#        C_g para todo pi, inclusive pi = infinito;
#   VIII o expoente de n e a constante do Teorema 5, com j_n do enunciado, ao
#        lado das cotas superiores do Corolário 11, do Teorema 4 e do
#        Corolário 5;
#   IX   a construção esparsa contra a densa, na escala de expoentes;
#   X    sem ruído gaussiano: a constante p_* da condição (2.29) de Tsybakov
#        (2009) para uma mistura gaussiana.
#
# Nos scripts de conferência, elemento de lista se acessa com [[ ]] e nome
# completo (instrucoes.md, §5).
#
# Dependências: WaveBased (wbasis, wtable).

suppressPackageStartupMessages(library(WaveBased))
set.seed(20261005)
t_start <- proc.time()[["elapsed"]]

fs <- 8L
tb <- wtable(filter.size = fs)
ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "ok" else "FALHA", msg))
  if (!isTRUE(cond)) ok <<- FALSE
  invisible(cond)
}
c_A <- (1 - 2^(-1/2)) / 2          # Assouad com KL <= 1: (M/2)(1 - sqrt(1/2)) = c_A M

## ---- o desenho -----------------------------------------------------------

p <- 2L; q <- 2L
B_X <- 2
a_dens <- 0.6; b_dens <- 0.3
dens_U <- function(u1, u2) 1 + a_dens * (u1 - 0.5) + b_dens * cos(2 * pi * (u1 + u2))
h_fun <- function(u1, u2) 0.5 + 0.5 * sin(2 * pi * u1) * u2
vV <- 1 / 3                          # Var(Unif(-1, 1))
c_U_th <- 1 - 0.5 * a_dens - b_dens  # 0.4
C_U_th <- 1 + 0.5 * a_dens + b_dens  # 1.6
lam_min_h <- function(h) { tr <- 1 + h^2 + vV; (tr - sqrt(tr^2 - 4 * vV)) / 2 }
lam_max_h <- function(h) { tr <- 1 + h^2 + vV; (tr + sqrt(tr^2 - 4 * vV)) / 2 }
kappa1_th <- lam_min_h(1)            # decrescente em h, h <= 1
kappa2_th <- lam_max_h(1)
gam <- kappa1_th * c_U_th            # gamma = kappa_1 c_U

# amostra do desenho: U por rejeição sob a uniforme, X_2 = h(U) + V
r_design <- function(n) {
  U <- matrix(numeric(0), 0, 2)
  while (nrow(U) < n) {
    m <- 2L * (n - nrow(U)) + 10L
    cand <- matrix(runif(2 * m), m, 2)
    acc <- runif(m) * C_U_th <= dens_U(cand[, 1], cand[, 2])
    U <- rbind(U, cand[acc, , drop = FALSE])
  }
  U <- U[seq_len(n), , drop = FALSE]
  X <- cbind(1, h_fun(U[, 1], U[, 2]) + runif(n, -1, 1))
  list(X = X, U = U)
}

# colunas do nível j na saída do wbasis() sem a constante: 2^j + (0:(2^j - 1))
lev_cols <- function(j) 2^j + seq_len(2^j) - 1L
psi_lev <- function(u, j, table = NULL) {
  B <- if (is.null(table)) wbasis(u, j0 = 0, J = j + 1L, filter.size = fs)
       else wbasis(u, j0 = 0, J = j + 1L, filter.size = fs, wavelet.table = table)
  B[, -1, drop = FALSE][, lev_cols(j), drop = FALSE]
}

# coordenadas do hipercubo: blocos K (matriz de pares (l, m)) no nível j
coords <- function(K, j) {
  do.call(rbind, lapply(seq_len(nrow(K)), function(r)
    cbind(l = K[r, 1], m = K[r, 2], k = 0:(2^j - 1))))
}
all_blocks <- as.matrix(expand.grid(l = seq_len(p), m = seq_len(q)))

# colunas do hipercubo numa amostra: z_a(X_i, U_i) = X_il psi_jk(U_im), a em C
# (tabela do wbasis, D31: as partes que usam isto são de Monte Carlo)
Zmat <- function(d, C, j) {
  Pm <- lapply(seq_len(q), function(m)
    if (any(C[, "m"] == m)) psi_lev(d[["U"]][, m], j, tb) else NULL)
  do.call(cbind, lapply(seq_len(nrow(C)), function(a)
    d[["X"]][, C[a, "l"]] * Pm[[C[a, "m"]]][, C[a, "k"] + 1L]))
}

## ---- quadratura em [0,1]^2 -----------------------------------------------

N2 <- 1024L
ug <- (seq_len(N2) - 0.5) / N2
U1g <- matrix(ug, N2, N2)                 # linhas: u_1
U2g <- matrix(ug, N2, N2, byrow = TRUE)   # colunas: u_2
Pg <- dens_U(U1g, U2g)
Hg <- h_fun(U1g, U2g)
wq <- Pg / N2^2
# W[[l, l']] = E(X_l X_l' | u) p_U(u) du na grade
Wq <- list("11" = wq, "12" = wq * Hg, "22" = wq * (Hg^2 + vV))
Wll <- function(l, lp) Wq[[paste0(min(l, lp), max(l, lp))]]
cat(sprintf("quadratura: massa de p_U = %.10f\n", sum(wq)))

Jg <- 6L                                   # níveis 0..5 na grade (estimadores vão até j + 2)
Bg <- wbasis(ug, j0 = 0, J = Jg, filter.size = fs)[, -1, drop = FALSE]
psi_grid <- function(j) Bg[, lev_cols(j), drop = FALSE]

# Gram em L_2(P) das colunas z_a(x, u) = x_l psi_jk(u_m), a em C, montada por
# pares de blocos: E{X_l X_l' psi(U_m) psi(U_m')} com a marginal em (u_m, u_m')
RSq <- lapply(Wq, rowSums); CSq <- lapply(Wq, colSums)
key <- function(l, lp) paste0(min(l, lp), max(l, lp))
gram_P <- function(C, j) {
  Ps <- psi_grid(j)
  M <- nrow(C); G <- matrix(0, M, M)
  blk <- unique(C[, c("l", "m"), drop = FALSE])
  for (r in seq_len(nrow(blk))) for (t in seq_len(nrow(blk))) {
    l <- blk[r, "l"]; m <- blk[r, "m"]; lp <- blk[t, "l"]; mp <- blk[t, "m"]
    ia <- which(C[, "l"] == l & C[, "m"] == m); ib <- which(C[, "l"] == lp & C[, "m"] == mp)
    Pa <- Ps[, C[ia, "k"] + 1L, drop = FALSE]; Pb <- Ps[, C[ib, "k"] + 1L, drop = FALSE]
    W <- Wq[[key(l, lp)]]
    G[ia, ib] <- if (m == 1 && mp == 1) crossprod(Pa, RSq[[key(l, lp)]] * Pb)
                 else if (m == 2 && mp == 2) crossprod(Pa, CSq[[key(l, lp)]] * Pb)
                 else if (m == 1 && mp == 2) crossprod(Pa, W %*% Pb)
                 else crossprod(Pa, crossprod(W, Pb))
  }
  (G + t(G)) / 2
}

# funções da forma f(x, u) = sum_l x_l R_l(u), R_l uma matriz N2 x N2
fun_from_coords <- function(C, j, theta) {
  Ps <- psi_grid(j)
  lapply(seq_len(p), function(l) {
    R <- matrix(0, N2, N2)
    for (a in which(C[, "l"] == l)) {
      v <- theta[a] * Ps[, C[a, "k"] + 1L]
      if (C[a, "m"] == 1) R <- R + v else R <- R + matrix(v, N2, N2, byrow = TRUE)
    }
    R
  })
}
norm2_P <- function(R) {                   # E{(sum_l X_l R_l(U))^2}
  s <- 0
  for (l in seq_len(p)) for (lp in seq_len(p)) s <- s + sum(Wll(l, lp) * R[[l]] * R[[lp]])
  s
}
proj_rhs <- function(C, j, R) {            # E{z_a(X, U) f(X, U)}, a em C
  Ps <- psi_grid(j); out <- numeric(nrow(C))
  for (l in seq_len(p)) {
    acc1 <- numeric(N2); acc2 <- numeric(N2)
    for (lp in seq_len(p)) {
      WR <- Wll(l, lp) * R[[lp]]
      acc1 <- acc1 + rowSums(WR); acc2 <- acc2 + colSums(WR)
    }
    for (a in which(C[, "l"] == l))
      out[a] <- sum((if (C[a, "m"] == 1) acc1 else acc2) * Ps[, C[a, "k"] + 1L])
  }
  out
}

## ---- I. a base -----------------------------------------------------------

cat("\nI. A base periodizada (ortonormalidade de Lebesgue e média zero)\n")
N1 <- 8192L
u1d <- (seq_len(N1) - 0.5) / N1
B1d <- wbasis(u1d, j0 = 0, J = Jg, filter.size = fs)[, -1, drop = FALSE]
G1d <- crossprod(B1d) / N1
chk(max(abs(G1d - diag(ncol(B1d)))) < 1e-6 && max(abs(colMeans(B1d))) < 1e-10,
    sprintf("{psi_jk}, j < %d: |Gram de Lebesgue - I| = %.1e, |média| = %.1e",
            Jg, max(abs(G1d - diag(ncol(B1d)))), max(abs(colMeans(B1d)))))
# na grade da quadratura 2D, nos níveis que as Partes III a V usam (j <= 4); o
# erro fica abaixo das folgas medidas na Parte V
cols_used <- seq_len(2^5 - 1)
G2d <- crossprod(Bg[, cols_used]) / N2
err2d <- max(abs(G2d - diag(length(cols_used))))
chk(err2d < 1e-4,
    sprintf("na grade da quadratura 2D (N = %d), níveis j <= 4: |Gram - I| = %.1e", N2, err2d))

## ---- II. constantes do desenho -------------------------------------------

cat("\nII. Constantes do desenho\n")
cU <- min(Pg); CU <- max(Pg)
lmin_g <- lam_min_h(Hg); lmax_g <- lam_max_h(Hg)
cat(sprintf("  c_U = %.4f, C_U = %.4f, kappa_1 = %.4f, kappa_2 = %.4f, gamma = %.4f, B_X = %g\n",
            c_U_th, C_U_th, kappa1_th, kappa2_th, gam, B_X))
chk(cU >= c_U_th - 1e-9 && CU <= C_U_th + 1e-9 && abs(cU - c_U_th) < 1e-3 && abs(CU - C_U_th) < 1e-3,
    sprintf("densidade conjunta em [c_U, C_U] (na grade: %.4f a %.4f)", cU, CU))
chk(min(lmin_g) >= kappa1_th - 1e-12 && max(lmax_g) <= kappa2_th + 1e-12,
    sprintf("autovalores de E(XX'|U) em [kappa_1, kappa_2] (na grade: %.4f a %.4f)", min(lmin_g), max(lmax_g)))
marg1 <- rowSums(Pg) / N2
chk(diff(range(marg1)) > 0.5,
    sprintf("a marginal de U_1 não é uniforme (varia de %.3f a %.3f)", min(marg1), max(marg1)))

## ---- III. Kullback-Leibler ----------------------------------------------

cat("\nIII. Kullback-Leibler dos vizinhos: E(X_l^2 psi_jk(U_m)^2) contra kappa_2 C_U e B_X^2 C_U\n")
tab3 <- do.call(rbind, lapply(0:4, function(j) {
  Ps <- psi_grid(j)
  do.call(rbind, lapply(seq_len(p), function(l) do.call(rbind, lapply(seq_len(q), function(m) {
    W <- Wll(l, l)
    marg <- if (m == 1) rowSums(W) else colSums(W)
    e <- colSums(marg * Ps^2)
    data.frame(j = j, l = l, m = m, e_max = max(e), e_min = min(e))
  }))))
}))
print(tab3, digits = 4, row.names = FALSE)
chk(all(tab3[["e_max"]] <= kappa2_th * C_U_th) && kappa2_th <= B_X^2,
    sprintf("max E(X_l^2 psi_jk(U_m)^2) = %.4f <= kappa_2 C_U = %.4f <= B_X^2 C_U = %.4f, em todo (l, m, j, k)",
            max(tab3[["e_max"]]), kappa2_th * C_U_th, B_X^2 * C_U_th))
chk(all(tab3[["e_min"]] >= kappa1_th * c_U_th),
    sprintf("e o mínimo %.4f >= kappa_1 c_U = %.4f (a separação de uma coordenada)",
            min(tab3[["e_min"]]), gam))

# a fórmula do KL com desenho aleatório, por Monte Carlo do log da razão de
# verossimilhança: E_omega log(dP_omega / dP_omega') = n zeta^2 v' G v / (2 sigma^2)
sigma <- 1
kl_mc <- function(C, j, w, wp, n, zeta, R = 4000L) {
  llr <- replicate(R, {
    d <- r_design(n)
    Z <- Zmat(d, C, j)
    f0 <- drop(Z %*% (zeta * w)); f1 <- drop(Z %*% (zeta * wp))
    Y <- f0 + sigma * rnorm(n)
    sum((Y - f1)^2 - (Y - f0)^2) / (2 * sigma^2)
  })
  c(media = mean(llr), ep = sd(llr) / sqrt(R))
}
C3 <- coords(all_blocks, 2L); G3 <- gram_P(C3, 2L)
n3 <- 200L; zeta3 <- sqrt(2 * sigma^2 / (n3 * B_X^2 * C_U_th))   # KL de vizinho <= 1
w0 <- rbinom(nrow(C3), 1, 0.5)
a_nb <- which(C3[, "l"] == 2 & C3[, "m"] == 2)[2]                  # coordenada de X_2 (dependente de U)
w1 <- w0; w1[a_nb] <- 1 - w1[a_nb]
w2 <- rbinom(nrow(C3), 1, 0.5)
kl_nb_th <- n3 * zeta3^2 * G3[a_nb, a_nb] / (2 * sigma^2)
kl_any_th <- n3 * zeta3^2 * drop(crossprod(w2 - w0, G3 %*% (w2 - w0))) / (2 * sigma^2)
mc_nb <- kl_mc(C3, 2L, w0, w1, n3, zeta3)
mc_any <- kl_mc(C3, 2L, w0, w2, n3, zeta3)
cat(sprintf("  vizinho (l, m) = (2, 2): KL = %.4f (quadratura), %.4f +- %.4f (Monte Carlo); cota do Lema 17 = 1\n",
            kl_nb_th, mc_nb[["media"]], mc_nb[["ep"]]))
cat(sprintf("  par qualquer (d_H = %d): KL = %.4f (quadratura), %.4f +- %.4f (Monte Carlo); cota d_H = %d\n",
            sum(w0 != w2), kl_any_th, mc_any[["media"]], mc_any[["ep"]], sum(w0 != w2)))
chk(abs(mc_nb[["media"]] - kl_nb_th) < 4 * mc_nb[["ep"]] && abs(mc_any[["media"]] - kl_any_th) < 4 * mc_any[["ep"]],
    "KL = n E(Delta^2)/(2 sigma^2) com X dependente de U e q = 2, a menos de 4 erros-padrão, nos dois pares")
chk(kl_nb_th <= 1 && kl_any_th <= sum(w0 != w2) * kappa2_th / B_X^2,
    "o KL do vizinho fica sob 1 com zeta^2 = 2 sigma^2/(n B_X^2 C_U); o de um par qualquer, sob d_H kappa_2/B_X^2")

## ---- IV. separação ------------------------------------------------------

cat("\nIV. Separação: lambda_min da Gram do hipercubo em L_2(P) contra gamma = kappa_1 c_U\n")
Ksets <- list("(1,1)" = all_blocks[1, , drop = FALSE],
              "(2,2)" = all_blocks[4, , drop = FALSE],
              "(2,1),(2,2)" = all_blocks[c(2, 4), , drop = FALSE],
              "(1,1),(1,2)" = all_blocks[c(1, 3), , drop = FALSE],
              "todos (pq = 4)" = all_blocks)
tab4 <- do.call(rbind, lapply(names(Ksets), function(nm) do.call(rbind, lapply(0:3, function(j) {
  C <- coords(Ksets[[nm]], j); G <- gram_P(C, j)
  ev <- eigen(G, symmetric = TRUE, only.values = TRUE)[["values"]]
  # termos cruzados: entradas entre blocos distintos
  same <- outer(C[, "l"], C[, "l"], "==") & outer(C[, "m"], C[, "m"], "==")
  cross <- if (any(!same)) max(abs(G[!same])) else 0
  data.frame(K = nm, j = j, M = nrow(C), lmin = min(ev), lmax = max(ev),
             razao = min(ev) / gam, cruzado = cross)
}))))
print(tab4, digits = 4, row.names = FALSE)
chk(all(tab4[["lmin"]] >= gam) && all(tab4[["lmax"]] <= kappa2_th * C_U_th),
    sprintf("gamma <= lambda_min <= lambda_max <= kappa_2 C_U em todo K e j <= 3 (razão mínima %.2f)",
            min(tab4[["razao"]])))
chk(max(tab4[["cruzado"]]) > 0.05,
    sprintf("os termos cruzados entre blocos estão presentes (até %.3f): a Gram não é bloco-diagonal",
            max(tab4[["cruzado"]])))
# pares quaisquer: ||f_w - f_w'||_P^2 >= gamma zeta^2 d_H, ||g_w - g_w'||^2 = zeta^2 d_H
C4 <- coords(all_blocks, 2L); G4 <- gram_P(C4, 2L)
viol4 <- 0; razao4 <- Inf; comp4 <- 0
Ps1 <- B1d[, lev_cols(2L), drop = FALSE]
for (r in 1:500) {
  w <- rbinom(nrow(C4), 1, 0.5); wp <- rbinom(nrow(C4), 1, 0.5); dH <- sum(w != wp)
  if (dH == 0) next
  dist <- drop(crossprod(w - wp, G4 %*% (w - wp)))
  if (dist < gam * dH) viol4 <- viol4 + 1
  razao4 <- min(razao4, dist / (gam * dH))
  # componentes: soma sobre os blocos da norma de Lebesgue, por quadratura em 1D
  dc <- 0
  for (b in seq_len(nrow(all_blocks))) {
    idx <- which(C4[, "l"] == all_blocks[b, 1] & C4[, "m"] == all_blocks[b, 2])
    dc <- dc + mean((Ps1 %*% (w[idx] - wp[idx]))^2)
  }
  comp4 <- max(comp4, abs(dc - dH))
}
chk(viol4 == 0, sprintf("||f_w - f_w'||_P^2 >= gamma zeta^2 d_H(w, w') em 500 pares (razão mínima %.2f)", razao4))
chk(comp4 < 1e-6, sprintf("sum ||g_w - g_w'||^2 = zeta^2 d_H(w, w') (erro %.1e)", comp4))

## ---- V. a redução -------------------------------------------------------

cat("\nV. A redução de um estimador qualquer a w_hat (Lema 17(iv))\n")
jV <- 2L; CV <- coords(all_blocks, jV); GV <- gram_P(CV, jV); MV <- nrow(CV)
lminV <- min(eigen(GV, symmetric = TRUE, only.values = TRUE)[["values"]])
zetaV <- 0.3
f_omega <- function(w) fun_from_coords(CV, jV, zetaV * w)
# estimador: f_{w'} + ruído em níveis 0..4 nas duas moduladoras + nível
# constante + termo não aditivo x_l (u_1 u_2 - 1/4) + termo em x_l fora da base
rand_est <- function(wp, amp) {
  R <- f_omega(wp)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      v <- drop(Bg[, seq_len(2^5 - 1)] %*% (amp * zetaV * rnorm(2^5 - 1) / sqrt(2^5)))
      R[[l]] <- R[[l]] + if (m == 1) v else matrix(v, N2, N2, byrow = TRUE)
    }
    R[[l]] <- R[[l]] + amp * zetaV * (rnorm(1) + rnorm(1) * (U1g * U2g - 0.25) +
                                        rnorm(1) * cos(6 * pi * U1g) * U2g)
  }
  R
}
resV <- t(vapply(1:200, function(r) {
  w <- rbinom(MV, 1, 0.5)
  wp <- if (r %% 3 == 0) w else rbinom(MV, 1, 0.5)
  amp <- c(0.05, 0.3, 1, 3)[(r %% 4) + 1]
  Rhat <- rand_est(wp, amp)
  Rw <- f_omega(w)
  D <- lapply(seq_len(p), function(l) Rhat[[l]] - Rw[[l]])
  LP <- norm2_P(D)
  a <- solve(GV, proj_rhs(CV, jV, Rhat))
  wh <- as.numeric(a >= zetaV / 2)
  dH <- sum(wh != w)
  quad <- drop(crossprod(a - zetaV * w, GV %*% (a - zetaV * w)))
  c(pit = LP - quad, eig = quad - lminV * sum((a - zetaV * w)^2),
    thr = sum((a - zetaV * w)^2) - zetaV^2 * dH / 4,
    fim = LP - gam * zetaV^2 * dH / 4, dH = dH)
}, numeric(5)))
chk(min(resV[, "pit"]) > -1e-8 && min(resV[, "eig"]) > -1e-10 && min(resV[, "thr"]) > -1e-12,
    sprintf("população: L_P >= (a - zeta w)'G(a - zeta w) >= lambda_min ||a - zeta w||^2 >= lambda_min zeta^2 d_H/4 (folgas mínimas %.1e, %.1e, %.1e)",
            min(resV[, "pit"]), min(resV[, "eig"]), min(resV[, "thr"])))
chk(min(resV[, "fim"]) >= 0 && max(resV[, "dH"]) > 0 && mean(resV[, "dH"] > 0) > 0.5,
    sprintf("L_P >= gamma zeta^2 d_H(w_hat, w)/4 em 200 estimadores (d_H até %d; %.0f%% com d_H > 0)",
            max(resV[, "dH"]), 100 * mean(resV[, "dH"] > 0)))
# componentes, em L_2[0,1] de Lebesgue, com estimadores fora da base
PsV1 <- B1d[, lev_cols(jV), drop = FALSE]
resC <- t(vapply(1:300, function(r) {
  w <- rbinom(2^jV, 1, 0.5); wp <- if (r %% 3 == 0) w else rbinom(2^jV, 1, 0.5)
  amp <- c(0.05, 0.3, 1, 3)[(r %% 4) + 1]
  g <- drop(PsV1 %*% (zetaV * w))
  gh <- drop(PsV1 %*% (zetaV * wp)) + amp * zetaV *
    (drop(B1d %*% (rnorm(ncol(B1d)) / sqrt(ncol(B1d)))) + rnorm(1) * cos(6 * pi * u1d) + rnorm(1) * u1d^2)
  Lc <- mean((gh - g)^2)
  a <- drop(crossprod(PsV1, gh)) / N1
  wh <- as.numeric(a >= zetaV / 2)
  c(bessel = Lc - sum((a - zetaV * w)^2), fim = Lc - zetaV^2 * sum(wh != w) / 4)
}, numeric(2)))
chk(min(resC[, "bessel"]) > -1e-9 && min(resC[, "fim"]) >= 0,
    sprintf("componentes: ||g_hat - g_w||^2 >= sum_k (<g_hat, psi_jk> - zeta w_k)^2 >= zeta^2 d_H/4 (folgas %.1e, %.1e)",
            min(resC[, "bessel"]), min(resC[, "fim"])))
# empírico, no evento lambda_min(Sigma_hat_CC) >= gamma/2
emp_once <- function(n, j, K) {
  C <- coords(K, j); d <- r_design(n)
  Z <- Zmat(d, C, j)
  S <- crossprod(Z) / n
  list(Z = Z, d = d, C = C, lmin = min(eigen(S, symmetric = TRUE, only.values = TRUE)[["values"]]))
}
viola_n <- 0; nA <- 0
for (r in 1:200) {
  e <- emp_once(200L, jV, all_blocks)
  if (e[["lmin"]] < gam / 2) next
  nA <- nA + 1
  w <- rbinom(MV, 1, 0.5); wp <- rbinom(MV, 1, 0.5)
  # estimador avaliado na amostra: f_{w'} + nível + termo não aditivo + ruído de nível 4
  X <- e[["d"]][["X"]]; U <- e[["d"]][["U"]]
  B4 <- wbasis(U[, 1], j0 = 0, J = 5L, filter.size = fs, wavelet.table = tb)[, -1]
  Fh <- drop(e[["Z"]] %*% (zetaV * wp)) + 0.3 * zetaV *
    (X[, 2] * (U[, 1] * U[, 2] - 0.25) + rnorm(1) + X[, 1] * drop(B4 %*% rnorm(ncol(B4))) / 4)
  Fw <- drop(e[["Z"]] %*% (zetaV * w))
  Ln <- mean((Fh - Fw)^2)
  a <- drop(solve(crossprod(e[["Z"]]), crossprod(e[["Z"]], Fh)))
  wh <- as.numeric(a >= zetaV / 2)
  if (Ln < gam * zetaV^2 * sum(wh != w) / 8 - 1e-12) viola_n <- viola_n + 1
}
chk(viola_n == 0 && nA > 0,
    sprintf("empírico (n = 200, M = %d): L_n >= gamma zeta^2 d_H/8 no evento, em %d amostras do evento", MV, nA))
freqA <- vapply(c(200L, 800L, 3200L), function(n)
  mean(replicate(200, emp_once(n, jV, all_blocks)[["lmin"]] >= gam / 2)), 0)
cat(sprintf("  P(lambda_min(Sigma_hat_CC) >= gamma/2), M = %d: n = 200: %.3f, 800: %.3f, 3200: %.3f\n",
            MV, freqA[1], freqA[2], freqA[3]))
chk(freqA[3] >= 0.95 && freqA[3] >= freqA[1],
    "a probabilidade do evento sobe com n e passa de 0.95 em n = 3200 (Proposição 3(iv) numa submatriz)")

## ---- VI. Assouad --------------------------------------------------------

cat("\nVI. Assouad: risco de Bayes de Hamming (priori uniforme) contra (M/2) max(e^{-alpha}/2, 1 - sqrt(alpha/2))\n")
bayes_risk <- function(K, j, n, zeta, R = 3000L) {
  C <- coords(K, j); M <- nrow(C)
  Om <- as.matrix(expand.grid(rep(list(0:1), M)))
  loss <- replicate(R, {
    d <- r_design(n)
    Z <- Zmat(d, C, j)
    w <- rbinom(M, 1, 0.5)
    Y <- drop(Z %*% (zeta * w)) + sigma * rnorm(n)
    Fit <- Z %*% t(zeta * Om)                      # n x 2^M
    ll <- -colSums((Y - Fit)^2) / (2 * sigma^2)
    post <- exp(ll - max(ll)); post <- post / sum(post)
    wh <- as.numeric(drop(crossprod(Om, post)) > 0.5)
    sum(wh != w)
  })
  c(risco = mean(loss), ep = sd(loss) / sqrt(R))
}
assouad_rhs <- function(M, alpha) (M / 2) * max(exp(-alpha) / 2, 1 - sqrt(alpha / 2))
tab6 <- do.call(rbind, lapply(list(list(K = all_blocks[c(2, 4), , drop = FALSE], j = 1L),
                                   list(K = all_blocks, j = 1L)), function(cf) {
  C <- coords(cf[["K"]], cf[["j"]]); M <- nrow(C); n <- 100L
  G <- gram_P(C, cf[["j"]]); kmax <- max(diag(G))
  out <- list()
  # (a) zeta do enunciado: KL de vizinho <= 1 pela cota B_X^2 C_U
  z_a <- sqrt(2 * sigma^2 / (n * B_X^2 * C_U_th))
  br <- bayes_risk(cf[["K"]], cf[["j"]], n, z_a)
  out[[1]] <- data.frame(M = M, zeta = "enunciado", KLmax = n * z_a^2 * kmax / (2 * sigma^2),
                         alpha = 1, rhs = assouad_rhs(M, 1), risco = br[["risco"]], ep = br[["ep"]])
  # (b) zeta com o KL exato de vizinho igual a 1 (a cota mais apertada)
  z_b <- sqrt(2 * sigma^2 / (n * kmax))
  br <- bayes_risk(cf[["K"]], cf[["j"]], n, z_b)
  out[[2]] <- data.frame(M = M, zeta = "KL exato = 1", KLmax = 1, alpha = 1,
                         rhs = assouad_rhs(M, 1), risco = br[["risco"]], ep = br[["ep"]])
  do.call(rbind, out)
}))
print(tab6, digits = 4, row.names = FALSE)
chk(all(tab6[["risco"]] - 3 * tab6[["ep"]] >= tab6[["rhs"]]),
    "o risco de Bayes de Hamming fica acima da cota de Assouad, com X dependente de U e os termos cruzados, nos dois zeta")

## ---- VII. pertencimento -------------------------------------------------

cat("\nVII. Pertencimento: a norma de Besov em sequência de g_omega, da função\n")
C_g <- 1
jn_fun <- function(n, s, Cg = C_g, sig = sigma, BX = B_X, CU = C_U_th) {
  x <- Cg^2 * n * BX^2 * CU / (2 * sig^2)
  if (x < 1) return(NA_integer_)
  floor(log2(x) / (2 * s + 1) + 1e-12)
}
zeta_fun <- function(n, sig = sigma, BX = B_X, CU = C_U_th) sqrt(2 * sig^2 / (n * BX^2 * CU))
seqnorm <- function(theta_by_level, s, pi_) {
  max(vapply(seq_along(theta_by_level), function(i) {
    j <- i - 1L; th <- theta_by_level[[i]]
    nrm <- if (is.infinite(pi_)) max(abs(th)) else sum(abs(th)^pi_)^(1 / pi_)
    2^(j * (s + 0.5 - if (is.infinite(pi_)) 0 else 1 / pi_)) * nrm
  }, 0))
}
NB <- 2^14
uB <- (seq_len(NB) - 0.5) / NB
tab7 <- do.call(rbind, lapply(c(0.3, 1, 2.5), function(s) {
  n <- 500; jn <- jn_fun(n, s); z <- zeta_fun(n)
  BB <- wbasis(uB, j0 = 0, J = jn + 3L, filter.size = fs)[, -1, drop = FALSE]
  do.call(rbind, lapply(c("um", "aleatorio"), function(tipo) {
    w <- if (tipo == "um") rep(1, 2^jn) else rbinom(2^jn, 1, 0.5)
    g <- drop(BB[, lev_cols(jn), drop = FALSE] %*% (z * w))
    th <- drop(crossprod(BB, g)) / NB                  # coeficientes da função, por quadratura
    by_lev <- lapply(0:(jn + 2L), function(j) th[lev_cols(j)])
    err <- max(abs(th[lev_cols(jn)] - z * w), abs(th[-lev_cols(jn)]))
    alvo <- if (tipo == "um") z * 2^(jn * (s + 0.5)) else NA_real_
    do.call(rbind, lapply(c(1, 1.5, 2, 4, Inf), function(pi_)
      data.frame(s = s, j_n = jn, omega = tipo, pi = pi_, norma = seqnorm(by_lev, s, pi_),
                 alvo = alvo, erro_coef = err)))
  }))
}))
print(tab7, digits = 4, row.names = FALSE)
um <- tab7[tab7[["omega"]] == "um", ]
chk(all(tab7[["norma"]] <= C_g * (1 + 1e-6)) && max(tab7[["erro_coef"]]) < 1e-6 &&
      max(abs(um[["norma"]] - um[["alvo"]])) < 1e-6,
    "||g_omega||_{b^s_{pi,inf}} = zeta 2^{j(s+1/2)} <= C_g para omega = 1, e <= C_g para omega qualquer, em todo pi (os coeficientes da função são os do hipercubo)")
chk(all(um[["norma"]] > C_g * 2^(-(um[["s"]] + 0.5))),
    "e j_n é o maior nível admissível: zeta 2^{j_n(s+1/2)} > C_g 2^{-(s+1/2)}")

## ---- VIII. o expoente de n e a constante --------------------------------

cat("\nVIII. O expoente de n e a constante do Teorema 5\n")
K_all <- p * q
lb_exact <- function(n, s) {                 # (c_A/4) K 2^{j_n} zeta_n^2: a cota das componentes
  jn <- jn_fun(n, s); (c_A / 4) * K_all * 2^jn * zeta_fun(n)^2
}
rho2 <- function(n, s) C_g^(2 / (2 * s + 1)) * (sigma^2 / (n * B_X^2 * C_U_th))^(2 * s / (2 * s + 1))
n0 <- ceiling(2 * sigma^2 / (C_g^2 * B_X^2 * C_U_th))
ngrid <- unique(round(10^seq(log10(max(n0, 10)), 12, length.out = 600)))
tab8 <- do.call(rbind, lapply(c(0.25, 0.5, 1, 2, 4), function(s) {
  lb <- vapply(ngrid, lb_exact, 0, s = s)
  base <- (c_A / 8) * K_all * vapply(ngrid, rho2, 0, s = s)
  rat <- lb / base
  slope <- unname(coef(lm(log(lb) ~ log(ngrid)))[2])
  # na sequência x_n = 2^{k(2s+1)} o piso não oscila: inclinação exata
  k <- 1:8; nk <- 2^(k * (2 * s + 1)) * 2 * sigma^2 / (C_g^2 * B_X^2 * C_U_th)
  lbk <- vapply(nk, lb_exact, 0, s = s)
  slope_k <- unname(coef(lm(log(lbk) ~ log(nk)))[2])
  data.frame(s = s, alvo = -2 * s / (2 * s + 1), incl = slope, incl_exata = slope_k,
             razao_min = min(rat), razao_max = max(rat),
             teto = 2^((4 * s + 1) / (2 * s + 1)))
}))
print(tab8, digits = 5, row.names = FALSE)
chk(all(tab8[["razao_min"]] >= 1) && all(tab8[["razao_max"]] <= tab8[["teto"]] * (1 + 1e-9)),
    "(c_A/4) K 2^{j_n} zeta_n^2 >= (c_A/8) K rho_n^2 em todo n >= n_0, e fica a um fator 2^{(4s+1)/(2s+1)} dela")
chk(all(abs(tab8[["incl_exata"]] - tab8[["alvo"]]) < 1e-9) && all(abs(tab8[["incl"]] - tab8[["alvo"]]) < 0.01),
    "o expoente de n é -2s/(2s+1): exato na sequência diádica, a 0.01 na grade de 600 valores")
# ao lado das cotas superiores: a razão cota superior / cota inferior em n
cat("  razão das taxas superiores à inferior (constantes à parte):\n")
tab8b <- do.call(rbind, lapply(list(c(1, 2), c(1, 1), c(2, 1.5), c(0.5, 4)), function(sp) {
  s <- sp[1]; pi_ <- sp[2]; se <- s - max(1 / pi_ - 0.5, 0)
  do.call(rbind, lapply(c(1e3, 1e6, 1e9), function(n) data.frame(
    s = s, pi = pi_, n = n,
    cor11 = (log(n))^(max(2 / pi_ - 1, 0) / (2 * s + 1)),
    teo4 = n^(-2 * se / (2 * se + 1)) / n^(-2 * s / (2 * s + 1)),
    cor5 = (log(n))^(2 * s / (2 * s + 1)))))
}))
print(tab8b, digits = 4, row.names = FALSE)
chk(all(tab8b[["cor11"]][tab8b[["pi"]] >= 2] == 1) && all(diff(tab8b[["cor11"]][tab8b[["pi"]] == 1]) > 0),
    "Corolário 11: razão 1 em pi >= 2 (ótimo) e crescente como (log n)^{(2/pi-1)/(2s+1)} em pi < 2")

## ---- IX. esparso contra denso -------------------------------------------

cat("\nIX. A construção esparsa (d picos em 2^j) contra a densa, em expoentes\n")
# 2^j = n^a, d = n^b, b = a - delta; Besov: b/pi + a(s + 1/2 - 1/pi) <= 1/2;
# risco ~ n^{b - 1} (com logaritmos); o denso é delta = 0.
best_b <- function(s, pi_, amax = 4) {
  agr <- seq(0, amax, length.out = 801); dgr <- seq(0, amax, length.out = 801)
  best <- c(b = -Inf, delta = NA)
  for (dl in dgr) {
    a <- agr[agr >= dl]
    ip <- if (is.infinite(pi_)) 0 else 1 / pi_
    okc <- (a - dl) * ip + a * (s + 0.5 - ip) <= 0.5 + 1e-12
    if (!any(okc)) next
    b <- max(a[okc] - dl)
    if (b > best[["b"]] + 1e-12) best <- c(b = b, delta = dl)
  }
  best
}
tab9 <- do.call(rbind, lapply(c(0.1, 0.3, 0.7, 1.5, 3), function(s)
  do.call(rbind, lapply(c(1, 1.3, 1.8, 2, 4, Inf), function(pi_) {
    se <- s - max((if (is.infinite(pi_)) 0 else 1 / pi_) - 0.5, 0)
    bb <- best_b(s, pi_)
    data.frame(s = s, pi = pi_, s_ef = se, b_max = bb[["b"]], delta = bb[["delta"]],
               denso = 1 / (2 * s + 1))
  }))))
# a fronteira s' = 0 fica fora das duas leituras (nela um nível não basta para
# ver que a classe deixa de ser limitada em L_2: é a soma sobre os níveis)
dentro <- tab9[tab9[["s_ef"]] > 1e-9, ]
fora <- tab9[tab9[["s_ef"]] < -1e-9, ]
print(tab9, digits = 4, row.names = FALSE)
chk(all(dentro[["b_max"]] <= dentro[["denso"]] + 1e-9) && all(dentro[["delta"]] == 0),
    sprintf("em s' > 0 (%d pares) nenhuma construção esparsa supera a densa: o máximo é delta = 0, b = 1/(2s+1)", nrow(dentro)))
if (nrow(fora) > 0)
  chk(all(fora[["b_max"]] >= 1),
      sprintf("em s' < 0 (%d pares, fora da hipótese) o expoente esparso passa de 1: a classe não é limitada em L_2", nrow(fora)))

## ---- X. sem ruído gaussiano ---------------------------------------------

cat("\nX. A condição (2.29) de Tsybakov: K(p, p(. + v)) <= p_* v^2 para |v| <= v_0\n")
# em log-densidades, numa janela onde as duas leis têm massa 1 a menos de 1e-20
kl_shift <- function(ldens, v) integrate(function(u) {
  la <- ldens(u); lb <- ldens(u + v)
  exp(la) * (la - lb)
}, -15, 15, rel.tol = 1e-12, subdivisions = 5000L)[["value"]]
lgauss <- function(u) dnorm(u, 0, sigma, log = TRUE)
lmix <- function(u) {
  a <- dnorm(u, -1, 0.5, log = TRUE); b <- dnorm(u, 1, 0.5, log = TRUE)
  mx <- pmax(a, b); log(0.5) + mx + log(exp(a - mx) + exp(b - mx))
}
vv <- c(1e-3, 0.01, 0.1, 0.3, 0.6, 1)
rg <- vapply(vv, function(v) kl_shift(lgauss, v) / v^2, 0)
rm <- vapply(vv, function(v) kl_shift(lmix, v) / v^2, 0)
fisher_mix <- integrate(function(u) {
  d1 <- 0.5 * dnorm(u, -1, 0.5) * (-(u + 1) / 0.25) + 0.5 * dnorm(u, 1, 0.5) * (-(u - 1) / 0.25)
  d1^2 / exp(lmix(u))
}, -12, 12, rel.tol = 1e-12, subdivisions = 5000L)[["value"]]
cat(sprintf("  gaussiana: K/v^2 = %s (1/(2 sigma^2) = %.4f)\n", paste(sprintf("%.4f", rg), collapse = " "), 1 / (2 * sigma^2)))
cat(sprintf("  mistura 0.5 N(-1, 1/4) + 0.5 N(1, 1/4): K/v^2 = %s; I/2 = %.4f; p_* (v_0 = 1) = %.4f\n",
            paste(sprintf("%.4f", rm), collapse = " "), fisher_mix / 2, max(rm)))
chk(max(abs(rg - 1 / (2 * sigma^2))) < 1e-6,
    "gaussiana: K(p, p(. + v)) = v^2/(2 sigma^2), isto é, p_* = 1/(2 sigma^2) para todo v_0")
chk(is.finite(max(rm)) && abs(rm[1] - fisher_mix / 2) / (fisher_mix / 2) < 1e-2,
    "mistura gaussiana (sub-gaussiana, não gaussiana): p_* finito, e K/v^2 -> I/2 quando v -> 0")

cat(sprintf("\n  tempo: %.0f s\n", proc.time()[["elapsed"]] - t_start))
if (ok) cat("OK\n") else stop("E1.9: conferência numérica FALHOU (ver linhas acima)")
