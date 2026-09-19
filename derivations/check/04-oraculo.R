# E1.5. Conferência numérica da desigualdade oráculo do LASSO no desenho de
# produtos (derivations/04-oraculo.tex). Roda ANTES da prova; imprime OK ou
# falha com stop(). Notação: docs/notacao.md (congelada em E1.1).
#
# O estimador conferido é
#     (c_hat, theta_hat) = argmin  ||Y - Z (c, theta)||_n^2 + 2 lambda ||theta||_1,
# com os p primeiros (os niveis c_l) NÃO penalizados e a ordem de colunas de
# D12. O ajuste é do glmnet, e a seção 0 verifica pelas condições KKT que o
# que o glmnet devolve é o minimizador DESTE objetivo (o glmnet reescala
# internamente penalty.factor para somar nvars, logo o lambda efetivo é
# lambda_glmnet * D / d; a seção 0 confere isso em vez de acreditar).
#
# O que confere, com J = 3, p = 2 (X_1 = 1), q = 2, Daublets com filtro 8,
# U ~ Unif[0,1]^2 (logo c_U = C_U = 1) e erro gaussiano:
#   0. KKT: A' r = 0 e |B' r / n|_inf <= lambda, com igualdade e sinal no
#      suporte ativo (o glmnet resolve o objetivo acima);
#   1. Lema do complemento de Schur: lambda_min(B~'B~/n) >= lambda_min(Sigma_hat),
#      com B~ = (I - P) B o bloco penalizado residualizado nos não penalizados.
#      É o que dispensa a condição de cone com coordenadas não penalizadas;
#   2. Lema das normas das colunas: sigma_max_hat^2 = max_a Sigma_hat_aa é
#      O_p(1) em n e NÃO cresce como 2^J, contra ||Z||_max^2 que cresce;
#   3. Lema de calibração: P(||B~' eps / n||_inf <= lambda_0) >= 1 - alpha com
#      lambda_0 = sigma sigma_max_hat sqrt(2 log(2d/alpha)/n);
#   4. Lema do ell_1 do oráculo em Besov: ||theta*||_1 <= p q C_g sum_j 2^{j(1/2-s)};
#   5. taxa lenta: ||B~(theta_hat - theta*)||_n^2 <= 6 lambda ||theta*||_1 + 4||b||_n^2
#      em TODA réplica, inclusive com theta* denso (sem condição de desenho);
#   6. taxa rápida: ||B~(theta_hat - theta*)||_n^2 <= 64 lambda^2 s_0 / gamma_til
#      + 16 ||b||_n^2, e o erro de predição ||f_hat - f||_n^2 escala como
#      lambda^2 s_0 ao variar n (a razão estabiliza);
#   7. termo das coordenadas não penalizadas: ||P eps||_n^2 ~ sigma^2 p / n;
#   8. cenário COM viés (g_11 = sin(2 pi u) fora do sieve): a desigualdade
#      completa vale e o viés empírico ||b||_n fica na ordem de ||f - f_J||_{L2(P)}.
#
# Dependências: WaveBased (wbasis, wtable), glmnet.
# Tempo nesta máquina: cerca de 40 s.

suppressPackageStartupMessages({
  library(WaveBased)
  library(glmnet)
})
set.seed(20260919)

J  <- 3L
p  <- 2L
q  <- 2L
fs <- 8L                        # filter.size (Daublets, 4 momentos nulos)
NJ <- 2L^J - 1L                 # wavelets por bloco (j0 = 0, constante descartada)
d  <- p * q * NJ                # colunas penalizadas
D  <- p + d                     # colunas de Z
sig <- 0.5                      # desvio do erro gaussiano
alpha <- 0.05
tb <- wtable(filter.size = fs)

ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "ok" else "FALHA", msg))
  if (!isTRUE(cond)) ok <<- FALSE
  invisible(cond)
}

## ---- desenho (mesma construção de check/03-desenho-produtos.R) ------------

psi_block <- function(u, Jl = J) {
  wbasis(u, j0 = 0, J = Jl, filter.size = fs, wavelet.table = tb)[, -1, drop = FALSE]
}
Psi <- function(U, Jl = J) {
  cbind(1, do.call(cbind, lapply(seq_len(ncol(U)), function(m) psi_block(U[, m], Jl))))
}
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
design <- function(X, U, Jl = J) {
  kron_rows(X, Psi(U, Jl))[, perm_D12(ncol(X), ncol(U), 2L^Jl - 1L), drop = FALSE]
}
draw_X <- function(U) cbind(1, 0.5 + runif(nrow(U), -1, 1))   # X indep. de U (cenário A de E1.4)
B_X <- 1.5

# residualização nos p não penalizados: B~ = (I - P) B, e P y
resid_pen <- function(Z) {
  A <- Z[, seq_len(p), drop = FALSE]
  Bm <- Z[, -seq_len(p), drop = FALSE]
  qa <- qr(A)
  list(A = A, B = Bm, Bt = Bm - qr.fitted(qa, Bm), proj = function(v) qr.fitted(qa, v))
}

## ---- ajuste do LASSO com os p primeiros não penalizados -------------------

# glmnet minimiza (1/2n)||y - Z beta||^2 + lambda_g sum_a pf_a |beta_a|, com pf
# reescalado para somar nvars. Nosso objetivo dividido por 2 é
# (1/2n)||y - Z beta||^2 + lambda ||theta||_1.
#
# ARMADILHA: o glmnet descarta colunas de variância zero, e a coluna não
# penalizada X_1 = 1 é uma delas; passada dentro de Z com penalty.factor = 0
# ela sai do ajuste (coeficiente identicamente nulo) e A' r não zera. As duas
# saídas, conferidas na seção 0:
#   (a) resolver o problema RESIDUALIZADO, theta_hat = argmin ||Y~ - B~ theta||_n^2
#       + 2 lambda ||theta||_1, e recuperar c_hat por mínimos quadrados. É a
#       perfilagem que a prova usa, e não tem coluna não penalizada nenhuma,
#       logo lambda_g = lambda sem reescalonamento;
#   (b) tirar a coluna constante de Z e deixá-la com o intercepto do glmnet
#       (que é não penalizado), com penalty.factor = 0 nas demais não
#       penalizadas; aí lambda_g = lambda * sum(pf) / ncol.
# (a) é o primário; (b) é a receita prática, conferida contra (a).
fit_resid <- function(rp, y, lambda) {
  yt <- y - rp[["proj"]](y)
  fit <- glmnet(rp[["Bt"]], yt, family = "gaussian", lambda = c(4, 2, 1) * lambda,
                standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-14))
  th <- as.numeric(fit[["beta"]][, 3])
  c(qr.coef(qr(rp[["A"]]), y - rp[["B"]] %*% th), th)
}
fit_joint <- function(Z, y, lambda) {           # rota (b), só para conferência
  Z2 <- Z[, -1, drop = FALSE]
  pf <- c(rep(0, p - 1), rep(1, d))
  lg <- lambda * sum(pf) / ncol(Z2)
  fit <- glmnet(Z2, y, family = "gaussian", lambda = c(4, 2, 1) * lg,
                penalty.factor = pf, standardize = FALSE, intercept = TRUE,
                control = list(thresh = 1e-14))
  c(fit[["a0"]][3], as.numeric(fit[["beta"]][, 3]))
}
kkt_dev <- function(Z, y, beta, lambda) {
  r <- as.numeric(y - Z %*% beta)
  gA <- crossprod(Z[, seq_len(p), drop = FALSE], r) / length(y)
  gB <- as.numeric(crossprod(Z[, -seq_len(p), drop = FALSE], r)) / length(y)
  th <- beta[-seq_len(p)]
  act <- which(abs(th) > 1e-10)
  c(unpen = max(abs(gA)),
    dual  = max(abs(gB)) - lambda,                                  # <= 0
    activ = if (length(act)) max(abs(gB[act] - lambda * sign(th[act]))) else 0)
}

## ---- oráculo e simulação ---------------------------------------------------

c_true <- c(1, -0.5)
make_theta <- function(s0, seed) {           # theta* esparso: s_0 entradas +-1
  set.seed(seed); th <- numeric(d)
  th[sample(d, s0)] <- sample(c(-1, 1), s0, TRUE)
  th
}
# theta* "denso", decaimento de Besov por nível dentro de cada bloco
lev_of <- rep(0:(J - 1), 2^(0:(J - 1)))                    # nível de cada wavelet do bloco
theta_besov <- function(s, Cg, seed) {
  set.seed(seed)
  unlist(lapply(seq_len(p * q), function(b) {
    out <- numeric(NJ)
    for (jj in 0:(J - 1)) {
      idx <- which(lev_of == jj); a <- rnorm(length(idx))
      out[idx] <- Cg * 2^(-jj * (s + 0.5 - 0.5)) * a / sqrt(sum(a^2))   # ||theta_j.||_2 = Cg 2^{-js}
    }
    out
  }))
}

sim_one <- function(n, theta, lambda_fun, g_extra = NULL) {
  U <- matrix(runif(n * q), n, q)
  X <- draw_X(U)
  Z <- design(X, U)
  fJ <- as.numeric(Z %*% c(c_true, theta))
  f  <- if (is.null(g_extra)) fJ else fJ + X[, 1] * g_extra(U[, 1])
  y  <- f + rnorm(n, sd = sig)
  rp <- resid_pen(Z)
  smax <- sqrt(max(colSums(Z[, -seq_len(p), drop = FALSE]^2) / n))
  lam <- lambda_fun(n, smax)
  beta <- fit_resid(rp, y, lam)
  v <- beta[-seq_len(p)] - theta
  b <- f - fJ
  list(n = n, lambda = lam, smax = smax, zmax = max(abs(Z[, -seq_len(p)])),
       kkt = kkt_dev(Z, y, beta, lam),
       dev_joint = max(abs(beta - fit_joint(Z, y, lam))),
       pred_Bt = mean((rp[["Bt"]] %*% v)^2),                 # ||B~ v||_n^2
       pred_f  = mean((as.numeric(Z %*% beta) - f)^2),       # ||f_hat - f||_n^2
       l1 = sum(abs(v)), l2sq = sum(v^2),
       bias_n = mean(b^2),                                   # ||b||_n^2
       score = max(abs(crossprod(rp[["Bt"]], y - f) / n)),   # ||B~' eps / n||_inf
       Peps = mean(rp[["proj"]](y - f)^2),                   # ||P eps||_n^2
       lmin_sig = min(eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]]),
       lmin_schur = min(eigen(crossprod(rp[["Bt"]]) / n, symmetric = TRUE,
                              only.values = TRUE)[["values"]]))
}
lam_theory <- function(n, smax) 2 * sig * smax * sqrt(2 * log(2 * d / alpha) / n)

## ===========================================================================
cat("0. KKT: o glmnet resolve ||y - Z(c,theta)||_n^2 + 2 lambda ||theta||_1\n")
## ===========================================================================

th4 <- make_theta(4L, 11L)
s0 <- sum(th4 != 0)
set.seed(101)
r0 <- lapply(1:20, function(i) sim_one(400L, th4, lam_theory))
kk <- do.call(rbind, lapply(r0, function(z) z[["kkt"]]))
cat(sprintf("  max |A' r / n| = %.2e; max (|B' r / n|_inf - lambda) = %.2e; max desvio no suporte ativo = %.2e\n",
            max(kk[, "unpen"]), max(kk[, "dual"]), max(kk[, "activ"])))
chk(max(kk[, "unpen"]) < 1e-8, "estacionariedade nas p colunas não penalizadas (A' r = 0)")
chk(max(kk[, "dual"]) < 1e-7, "viabilidade dual: |B' r / n|_inf <= lambda")
chk(max(kk[, "activ"]) < 1e-6, "igualdade com sinal no suporte ativo")
dj <- vapply(r0, function(z) z[["dev_joint"]], 1)
cat(sprintf("  rota residualizada vs rota do intercepto do glmnet: desvio máximo %.2e\n", max(dj)))
chk(max(dj) < 1e-7, "a perfilagem dos não penalizados devolve o mesmo minimizador que o ajuste conjunto")

## ===========================================================================
cat("\n1. Lema do complemento de Schur: lambda_min(B~'B~/n) >= lambda_min(Sigma_hat)\n")
## ===========================================================================

ls_ <- vapply(r0, function(z) z[["lmin_schur"]], 1)
lg_ <- vapply(r0, function(z) z[["lmin_sig"]], 1)
gam <- 0.25                     # kappa_1 c_U no cenário A de E1.4: lambda_min(E[XX']) x 1
cat(sprintf("  n = 400: mediana lambda_min(Sigma_hat) = %.4f; mediana lambda_min do Schur = %.4f; razão mínima = %.3f\n",
            median(lg_), median(ls_), min(ls_ / lg_)))
chk(all(ls_ >= lg_ - 1e-10), "vale em todas as réplicas (é a desigualdade usada na prova)")
chk(min(ls_) > 0, "o bloco penalizado residualizado tem posto cheio")
cat(sprintf("  cadeia consumida pelo teorema: gamma_til >= lambda_min(Sigma_hat) >= gamma/2 = kappa_1 c_U / 2 = %.4f em %.0f%% das réplicas\n",
            gam / 2, 100 * mean(lg_ >= gam / 2)))
chk(mean(lg_ >= gam / 2) >= 0.9,
    "a forma consumível de E1.4 (lambda_min(Sigma_hat) >= kappa_1 c_U / 2) vale em n = 400")

## ===========================================================================
cat("\n2. Lema das normas das colunas: sigma_max_hat = O_p(1), contra ||Z||_max ~ 2^{J/2}\n")
## ===========================================================================

for (nn in c(200L, 800L, 3200L)) {
  rr <- lapply(1:20, function(i) sim_one(nn, th4, lam_theory))
  cat(sprintf("  n = %5d: mediana sigma_max_hat = %.3f (cota populacional sqrt(B_X^2 C_U) = %.3f); mediana ||Z||_max = %.3f\n",
              nn, median(vapply(rr, function(z) z[["smax"]], 1)), B_X,
              median(vapply(rr, function(z) z[["zmax"]], 1))))
  if (nn == 3200L) {
    chk(median(vapply(rr, function(z) z[["smax"]], 1)) < 1.6 * B_X,
        "sigma_max_hat fica na ordem de B_X sqrt(C_U), sem fator 2^{J/2}")
  }
}
smax_J <- vapply(2:6, function(Jl) {
  U <- matrix(runif(4000 * q), 4000, q); X <- draw_X(U); Z <- design(X, U, Jl)
  c(smax = sqrt(max(colSums(Z[, -seq_len(p), drop = FALSE]^2) / 4000)),
    zmax = max(abs(Z[, -seq_len(p)])))
}, c(0, 0))
colnames(smax_J) <- paste0("J", 2:6)
cat("  em J = 2..6 (n = 4000): sigma_max_hat =", paste(sprintf("%.2f", smax_J["smax", ]), collapse = " "),
    "; ||Z||_max =", paste(sprintf("%.2f", smax_J["zmax", ]), collapse = " "), "\n")
chk(max(smax_J["smax", ]) < 2 * B_X && smax_J["zmax", "J6"] / smax_J["zmax", "J2"] > 3,
    "sigma_max_hat é estável em J e ||Z||_max cresce como 2^{J/2}: a calibração usa o primeiro")

## ===========================================================================
cat("\n3. Lema de calibração: P(||B~' eps / n||_inf <= lambda_0) >= 1 - alpha\n")
## ===========================================================================

R <- 400L
for (nn in c(200L, 1600L)) {
  rr <- lapply(seq_len(R), function(i) sim_one(nn, th4, lam_theory))
  sc <- vapply(rr, function(z) z[["score"]], 1)
  l0 <- vapply(rr, function(z) z[["lambda"]] / 2, 1)
  cov_ <- mean(sc <= l0)
  cat(sprintf("  n = %5d: lambda_0 mediano = %.4f; ||B~' eps/n||_inf mediano = %.4f; cobertura = %.3f (nominal >= %.2f)\n",
              nn, median(l0), median(sc), cov_, 1 - alpha))
  chk(cov_ >= 1 - alpha, "o evento de calibração tem a probabilidade nominal")
}

## ===========================================================================
cat("\n4. Lema do ell_1 do oráculo em Besov: ||theta_j.||_1 <= 2^{j(1-1/pi)} ||theta_j.||_pi\n")
## ===========================================================================

set.seed(7)
dev4 <- c()
for (pii in c(1, 1.5, 2, 4, Inf)) for (jj in 2:8) {
  a <- rnorm(2^jj) * (runif(2^jj) < 0.4)
  npi <- if (is.infinite(pii)) max(abs(a)) else sum(abs(a)^pii)^(1 / pii)
  dev4 <- c(dev4, sum(abs(a)) / (2^(jj * (1 - 1 / pii)) * npi))
}
cat(sprintf("  razão ||a||_1 / (2^{j(1-1/pi)} ||a||_pi) sobre pi em {1,1.5,2,4,inf} e j = 2..8: máximo %.4f\n",
            max(dev4)))
chk(max(dev4) <= 1 + 1e-12, "a desigualdade de Hölder por nível vale (é o Passo 1 do lema)")
# a soma: sob ||theta_{lm,j.}||_pi <= Cg 2^{-j(s+1/2-1/pi)}, ||theta*||_1 <= p q Cg sum_{j<J} 2^{j(1/2-s)}
for (s in c(0.3, 0.5, 1.5)) {
  Cg <- 1
  th <- theta_besov(s, Cg, 42L)
  bound <- p * q * Cg * sum(2^((0:(J - 1)) * (0.5 - s)))
  cat(sprintf("  s = %.1f: ||theta*||_1 = %.3f <= p q Cg sum_{j<J} 2^{j(1/2-s)} = %.3f\n",
              s, sum(abs(th)), bound))
  chk(sum(abs(th)) <= bound + 1e-9, "a cota do lema vale para a sequência de Besov construída")
}

## ===========================================================================
cat("\n5. Taxa lenta: ||B~ v||_n^2 <= 6 lambda ||theta*||_1 + 4 ||b||_n^2, sem condição de desenho\n")
## ===========================================================================

for (cs in list(list(nm = "theta* esparso (s_0 = 4)", th = th4),
                list(nm = "theta* denso de Besov (s = 0.3)", th = theta_besov(0.3, 1, 42L)))) {
  rr <- lapply(1:60, function(i) sim_one(800L, cs[["th"]], lam_theory))
  lhs <- vapply(rr, function(z) z[["pred_Bt"]], 1)
  rhs <- vapply(rr, function(z) 6 * z[["lambda"]] * sum(abs(cs[["th"]])) + 4 * z[["bias_n"]], 1)
  cat(sprintf("  %s: n = 800, ||theta*||_1 = %.2f; mediana LHS = %.5f; mediana RHS = %.5f; folga mínima = %.2fx\n",
              cs[["nm"]], sum(abs(cs[["th"]])), median(lhs), median(rhs), min(rhs / lhs)))
  chk(all(lhs <= rhs), "a desigualdade lenta vale em todas as réplicas")
}

## ===========================================================================
cat("\n6. Taxa rápida: ||B~ v||_n^2 <= 64 lambda^2 s_0 / gamma_til, e escala de lambda^2 s_0 em n\n")
## ===========================================================================

ns <- c(200L, 400L, 800L, 1600L, 3200L, 6400L)
tabF <- t(vapply(ns, function(nn) {
  rr <- lapply(1:60, function(i) sim_one(nn, th4, lam_theory))
  lhs <- vapply(rr, function(z) z[["pred_Bt"]], 1)
  gt  <- vapply(rr, function(z) z[["lmin_schur"]], 1)
  lam <- vapply(rr, function(z) z[["lambda"]], 1)
  pf  <- vapply(rr, function(z) z[["pred_f"]], 1)
  c(lambda = median(lam), pred = median(pf), Btv = median(lhs),
    razao = median(pf / (lam^2 * s0)), n_pred = nn * median(pf),
    viola = mean(lhs > 64 * lam^2 * s0 / gt),
    folga = median(64 * lam^2 * s0 / gt / lhs),
    l1 = median(vapply(rr, function(z) z[["l1"]], 1)),
    l1b = median(40 * lam * s0 / gt),
    l2 = median(vapply(rr, function(z) z[["l2sq"]], 1)),
    l2b = median(64 * lam^2 * s0 / gt^2))
}, numeric(11)))
rownames(tabF) <- paste0("n=", ns)
print(round(tabF[, c("lambda", "pred", "Btv", "razao", "n_pred", "viola", "folga")], 4))
chk(all(tabF[, "viola"] == 0), "a cota rápida 64 lambda^2 s_0 / gamma_til vale em todas as réplicas e todos os n")
rz <- tabF[, "razao"]
cat(sprintf("  razão ||f_hat - f||_n^2 / (lambda^2 s_0): %s (máx/mín = %.2f)\n",
            paste(sprintf("%.3f", rz), collapse = " "), max(rz) / min(rz)))
chk(max(rz) / min(rz) < 2, "a razão ao alvo lambda^2 s_0 estabiliza ao variar n por um fator 32")
np <- tabF[, "n_pred"]
cat(sprintf("  n * ||f_hat - f||_n^2: %s (a taxa é 1/n a J fixo)\n", paste(sprintf("%.2f", np), collapse = " ")))
chk(max(np) / min(np) < 2.5, "n * erro de predição é estável: o erro cai como 1/n")
cat(sprintf("  ell_1: mediana ||theta_hat - theta*||_1 = %s contra a cota 40 lambda s_0 / gamma_til = %s\n",
            paste(sprintf("%.3f", tabF[, "l1"]), collapse = " "),
            paste(sprintf("%.3f", tabF[, "l1b"]), collapse = " ")))
cat(sprintf("  ell_2: mediana ||theta_hat - theta*||_2^2 = %s contra a cota 64 lambda^2 s_0 / gamma_til^2 = %s\n",
            paste(sprintf("%.4f", tabF[, "l2"]), collapse = " "),
            paste(sprintf("%.4f", tabF[, "l2b"]), collapse = " ")))
chk(all(tabF[, "l1"] <= tabF[, "l1b"]) && all(tabF[, "l2"] <= tabF[, "l2b"]),
    "as cotas ell_1 e ell_2 do corolário valem nas medianas")

## ===========================================================================
cat("\n7. Termo das coordenadas não penalizadas: ||P eps||_n^2 ~ sigma^2 p / n\n")
## ===========================================================================

tab7 <- t(vapply(c(200L, 800L, 3200L), function(nn) {
  rr <- lapply(1:100, function(i) sim_one(nn, th4, lam_theory))
  pe <- vapply(rr, function(z) z[["Peps"]], 1)
  c(media = mean(pe), alvo = sig^2 * p / nn, razao = mean(pe) / (sig^2 * p / nn))
}, numeric(3)))
rownames(tab7) <- paste0("n=", c(200, 800, 3200))
print(signif(tab7, 4))
chk(all(abs(tab7[, "razao"] - 1) < 0.2),
    "E||P eps||_n^2 = sigma^2 p / n (a cota de Markov do teorema é exata em média)")

## ===========================================================================
cat("\n8. Cenário com viés: g_11(u) = sin(2 pi u) fora do sieve (J = 3)\n")
## ===========================================================================

gsin <- function(u) sin(2 * pi * u)
# projeção de sin(2 pi u) em W_3 por quadratura, para separar theta* do viés
ug <- (seq_len(2^14) - 0.5) / 2^14
Wg <- psi_block(ug, J)
th_sin <- drop(crossprod(Wg, gsin(ug))) / length(ug)
bias_L2 <- sqrt(mean((gsin(ug) - Wg %*% th_sin)^2))      # ||g - Pi_J g||_{L2[0,1]}
th_b <- numeric(d); th_b[seq_len(NJ)] <- th_sin          # bloco (l,m) = (1,1), primeiro do vetor
g_res <- function(u) gsin(u) - drop(psi_block(u, J) %*% th_sin)
cat(sprintf("  ||g - Pi_J g||_{L2[0,1]} = %.4f (C_U = 1, X_1 = 1, logo ||f - f_J||_{L2(P)} = o mesmo)\n", bias_L2))
for (nn in c(400L, 1600L)) {
  rr <- lapply(1:60, function(i) sim_one(nn, th_b, lam_theory, g_extra = g_res))
  bn <- vapply(rr, function(z) z[["bias_n"]], 1)
  lhs <- vapply(rr, function(z) z[["pred_Bt"]], 1)
  gt <- vapply(rr, function(z) z[["lmin_schur"]], 1)
  lam <- vapply(rr, function(z) z[["lambda"]], 1)
  s0b <- sum(th_b != 0)
  rhs <- 64 * lam^2 * s0b / gt + 16 * bn
  pf <- vapply(rr, function(z) z[["pred_f"]], 1)
  pe <- vapply(rr, function(z) z[["Peps"]], 1)
  rhs_pred <- 3 * (lhs + 4 * bn + pe)
  cat(sprintf("  n = %4d: mediana ||b||_n^2 = %.5f (populacional %.5f); mediana ||B~v||_n^2 = %.5f <= cota %.5f; violações %d\n",
              nn, median(bn), bias_L2^2, median(lhs), median(rhs), sum(lhs > rhs)))
  chk(all(lhs <= rhs), "a desigualdade rápida com viés vale em todas as réplicas")
  chk(all(pf <= rhs_pred), "||f_hat - f||_n^2 <= 3(||B~v||_n^2 + 4||b||_n^2 + ||P eps||_n^2)")
  chk(abs(median(bn) / bias_L2^2 - 1) < 0.25,
      "||b||_n^2 concentra em ||f - f_J||_{L2(P)}^2 (a passagem por Markov do teorema)")
}

cat("\n")
if (ok) cat("OK\n") else stop("E1.5: conferência numérica FALHOU (ver linhas acima)")
