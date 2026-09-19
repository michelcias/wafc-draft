# E1.6. Conferência numérica das taxas do LASSO no desenho de produtos
# (derivations/05-taxas.tex). Roda ANTES da prova; imprime OK ou falha com
# stop(). Notação: docs/notacao.md (congelada em E1.1).
#
# O que confere, em cinco partes:
#
#   I.  (determinístico, sem simulação)
#       1. Lema 9(i): a Hipótese de Besov de E1.3 implica weak-ell_tau com
#          tau = (s + 1/2)^{-1} e constante C_tau = Cg (p q A)^{s+1/2},
#          A = 1 + (1 - 2^{1 - pi(s+1/2)})^{-1}; confere também a contagem
#          N(eps) <= p q A (Cg/eps)^tau, que é o passo da prova;
#       2. Lema 9(ii): a cauda da melhor aproximação com s termos,
#          sum_{k>s} theta_(k)^2 <= C_tau^2 tau/(2-tau) s^{1-2/tau};
#       3. a álgebra do balanço do Corolário 5: o minimizador em s de
#          h(s) = a s + b s^{1-2/tau} e o valor h(s*) = 2 a s*/(2-tau);
#       4. os expoentes do Teorema 2 (os TRÊS termos -- estimação, viés e a
#          cota lenta lambda ||theta*||_1 -- têm a mesma ordem em J_n) e a
#          janela de c do Corolário 5: max{tau/2, (2-tau)/(4 s')} <= c < 1
#          é não vazia sse s' > (2-tau)/4, que com tau = (s+1/2)^{-1} é
#          s' > s/(2s+1).
#
#   II. o regime J_n = min{J : 2^J >= (n/log n)^{1/(2s'+1)}}: a condição
#       empírica de E1.4 (p q 2^J log(p q 2^J)/n -> 0) e a cadeia
#       gamma_til >= lambda_min(Sigma_hat) >= kappa_1 c_U / 2 ao longo dela.
#
#   III+IV. simulação com pi = 2 e s = s' = 1/4 (caso denso, s_0 = d): o erro
#       de predição em J_n contra (J_n/n)^{2s'/(2s'+1)} (Teorema 2), o erro
#       L_2 das componentes (Corolário 4), e a comparação de J_n com o J que
#       minimiza o erro realizado.
#
#   V.  simulação com pi = 1 (caso comprimível): o Lema 8 (comparador
#       arbitrário) para TODO s, a dimensão efetiva s*_n = argmin da cota, a
#       sua escala n^{tau/2}, a escala lambda^{2-tau} do valor mínimo e o
#       ganho sobre a leitura densa s = d.
#
# Dependências: WaveBased (wbasis, wtable), glmnet.
# Tempo nesta máquina: cerca de 7 min.

suppressPackageStartupMessages({
  library(WaveBased)
  library(glmnet)
})
set.seed(20260919)

p  <- 2L                        # X_1 = 1, X_2 = 1/2 + Unif(-1,1)  (cenário A de E1.4)
q  <- 2L
fs <- 8L                        # filter.size (Daublets, 4 momentos nulos)
sig <- 0.5
alpha <- 0.05
B_X <- 1.5
gam <- 0.25                     # kappa_1 c_U no cenário A (U uniforme: c_U = C_U = 1)
tb <- wtable(filter.size = fs)

ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "ok" else "FALHA", msg))
  if (!isTRUE(cond)) ok <<- FALSE
  invisible(cond)
}

## ---- desenho (mesma construção de check/03 e check/04) --------------------

psi_block <- function(u, Jl) {
  wbasis(u, j0 = 0, J = Jl, filter.size = fs, wavelet.table = tb)[, -1, drop = FALSE]
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
design_from <- function(X, Wlist, Jl) {           # Wlist[[m]] já avaliada em J >= Jl
  NJ <- 2L^Jl - 1L
  P <- cbind(1, do.call(cbind, lapply(Wlist, function(W) W[, seq_len(NJ), drop = FALSE])))
  kron_rows(X, P)[, perm_D12(ncol(X), length(Wlist), NJ), drop = FALSE]
}
resid_pen <- function(Z) {
  A <- Z[, seq_len(p), drop = FALSE]
  Bm <- Z[, -seq_len(p), drop = FALSE]
  qa <- qr(A)
  list(A = A, B = Bm, Bt = Bm - qr.fitted(qa, Bm), proj = function(v) qr.fitted(qa, v))
}
# LASSO com os p níveis fora da penalidade, pela rota residualizada (a
# perfilagem do Lema 4 de E1.5; é a rota conferida por KKT em check/04).
fit_resid <- function(rp, y, lambda) {
  yt <- y - rp[["proj"]](y)
  fit <- glmnet(rp[["Bt"]], yt, family = "gaussian", lambda = c(4, 2, 1) * lambda,
                standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-14))
  th <- as.numeric(fit[["beta"]][, 3])
  list(c = drop(qr.coef(qr(rp[["A"]]), y - rp[["B"]] %*% th)), theta = th)
}

## ---- sequências de Besov ---------------------------------------------------

lev_of <- function(Jl) rep(0:(Jl - 1L), 2^(0:(Jl - 1L)))     # nível de cada wavelet
# theta com ||theta_{j.}||_pi = Cg 2^{-j(s + 1/2 - 1/pi)} em TODO nível (satura
# a Hipótese de Besov de E1.3). spread = "uniform" dá o perfil weak-ell_tau
# justo; spread = "random" dá uma direção aleatória dentro do nível.
theta_besov <- function(s, pii, Cg, Jl, nblocks, spread = c("random", "uniform"), seed = 1L) {
  spread <- match.arg(spread)
  set.seed(seed)
  lv <- lev_of(Jl)
  unlist(lapply(seq_len(nblocks), function(bb) {
    out <- numeric(length(lv))
    for (jj in 0:(Jl - 1L)) {
      idx <- which(lv == jj)
      a <- if (spread == "uniform") rep(1, length(idx)) * sample(c(-1, 1), length(idx), TRUE)
           else rnorm(length(idx))
      npi <- if (is.infinite(pii)) max(abs(a)) else sum(abs(a)^pii)^(1 / pii)
      out[idx] <- Cg * 2^(-jj * (s + 0.5 - 1 / pii)) * a / npi
    }
    out
  }))
}
tau_of <- function(s) 1 / (s + 0.5)
A_of <- function(s, pii) 1 + 1 / (1 - 2^(1 - pii * (s + 0.5)))
Ctau_of <- function(s, pii, Cg, nblocks) Cg * (nblocks * A_of(s, pii))^(s + 0.5)

## ===========================================================================
cat("PARTE I. Lema 9, a álgebra do balanço e os expoentes (determinístico)\n")
cat("\n1. Lema 9(i): Besov  =>  weak-ell_tau com tau = (s + 1/2)^{-1}\n")
## ===========================================================================

Jl <- 12L
grid1 <- expand.grid(s = c(0.4, 0.8, 1.5, 3.0), pii = c(1, 1.5, 2, 4),
                     KEEP.OUT.ATTRS = FALSE)
grid1 <- grid1[grid1[["s"]] > 1 / grid1[["pii"]] - 0.5, ]      # s' > 0 (E1.3)
row1 <- t(vapply(seq_len(nrow(grid1)), function(i) {
  s <- grid1[["s"]][i]; pii <- grid1[["pii"]][i]; Cg <- 1
  tau <- tau_of(s); Ct <- Ctau_of(s, pii, Cg, p * q)
  res <- vapply(c("random", "uniform"), function(sp) {
    th <- theta_besov(s, pii, Cg, Jl, p * q, spread = sp, seed = 17L)
    srt <- sort(abs(th), decreasing = TRUE)
    kk <- seq_along(srt)
    # (a) envelope ordenado; (b) contagem N(eps) <= p q A (Cg/eps)^tau
    eps <- srt[srt > 0]
    Ne <- vapply(eps, function(e) sum(abs(th) > e * (1 - 1e-12)), 1)
    c(max(srt / (Ct * kk^(-1 / tau))),
      max(Ne / (p * q * A_of(s, pii) * (Cg / eps)^tau)))
  }, c(0, 0))
  c(s = s, pi = pii, tau = tau, Ctau = Ct, env = max(res[1, ]), cont = max(res[2, ]))
}, numeric(6)))
print(round(row1, 4))
chk(max(row1[, "env"]) <= 1, "theta*_(k) <= C_tau k^{-1/tau} em todas as sequências (envelope do lema)")
chk(max(row1[, "cont"]) <= 1, "N(eps) <= p q A (Cg/eps)^tau (o passo de contagem da prova)")

cat("\n2. Lema 9(ii): cauda da melhor aproximação com s termos\n")
tail2 <- t(vapply(seq_len(nrow(grid1)), function(i) {
  s <- grid1[["s"]][i]; pii <- grid1[["pii"]][i]; Cg <- 1
  tau <- tau_of(s); Ct <- Ctau_of(s, pii, Cg, p * q)
  th <- theta_besov(s, pii, Cg, Jl, p * q, spread = "uniform", seed = 17L)
  srt <- sort(abs(th), decreasing = TRUE)
  ss <- c(1L, 4L, 16L, 64L, 256L, 1024L)
  r <- vapply(ss, function(sv) {
    lhs <- sum(srt[-seq_len(sv)]^2)
    rhs <- Ct^2 * tau / (2 - tau) * sv^(1 - 2 / tau)
    lhs / rhs
  }, 1)
  c(s = s, pi = pii, razao_max = max(r))
}, numeric(3)))
print(round(tail2, 5))
chk(max(tail2[, "razao_max"]) <= 1,
    "sum_{k>s} theta_(k)^2 <= C_tau^2 tau/(2-tau) s^{1-2/tau} em s = 1..1024")

cat("\n3. Álgebra do balanço (Corolário 5): s* e h(s*) = 2 a s*/(2 - tau)\n")
dev3 <- c()
for (tau in c(0.3, 0.588, 0.8, 1.2, 1.8)) for (a in c(1e-4, 1e-2)) for (b in c(1e-3, 1, 10)) {
  h <- function(sv) a * sv + b * sv^(1 - 2 / tau)
  sstar <- (b * (2 - tau) / (a * tau))^(tau / 2)
  num <- optimize(function(lg) h(exp(lg)), log(c(1e-8, 1e12)), tol = 1e-12)
  dev3 <- c(dev3, abs(sstar / exp(num[["minimum"]]) - 1),
            abs(h(sstar) / (2 * a * sstar / (2 - tau)) - 1),
            abs(h(sstar) / num[["objective"]] - 1))
  # arredondar para o inteiro acima custa no máximo um fator 2 quando s* >= 1
  if (sstar >= 1) dev3 <- c(dev3, max(0, h(ceiling(sstar)) / (2 * h(sstar)) - 1))
}
cat(sprintf("  desvio máximo entre a forma fechada e a minimização numérica: %.2e\n", max(dev3)))
chk(max(dev3) < 1e-6,
    "s* = [b(2-tau)/(a tau)]^{tau/2}, h(s*) = 2 a s*/(2-tau) e h(ceil(s*)) <= 2 h(s*)")

cat("\n4. Expoentes do Teorema 2 e janela de c do Corolário 5\n")
# com 2^{J} = (n/L)^{1/(2s'+1)}, os três termos -- estimação 2^J L/n, viés
# 2^{-2Js'} e a cota lenta sqrt(L/n) 2^{J(1/2-s)} quando s < 1/2 -- têm a
# MESMA ordem (L/n)^{2s'/(2s'+1)}
dev4 <- c()
for (sp in c(0.25, 0.5, 1, 1.5, 3)) for (n in 10^c(3, 5, 8, 12)) {
  L <- log(n); tJ <- (n / L)^(1 / (2 * sp + 1))
  est <- tJ * L / n; vies <- tJ^(-2 * sp); alvo <- (L / n)^(2 * sp / (2 * sp + 1))
  dev4 <- c(dev4, abs(est / alvo - 1), abs(vies / alvo - 1))
  if (sp < 0.5) {                             # pi >= 2, s = s' < 1/2: ||theta*||_1 ~ 2^{J(1/2-s)}
    lento <- sqrt(L / n) * tJ^(0.5 - sp)
    dev4 <- c(dev4, abs(lento / alvo - 1))
  }
}
cat(sprintf("  desvio máximo dos termos ao alvo (L/n)^{2s'/(2s'+1)}: %.2e\n", max(dev4)))
chk(max(dev4) < 1e-10,
    "o balanço 2^{J_n} = (n/log n)^{1/(2s'+1)} iguala estimação, viés e (em s < 1/2) a cota lenta")
win <- t(vapply(seq_len(nrow(grid1)), function(i) {
  s <- grid1[["s"]][i]; pii <- grid1[["pii"]][i]
  sp <- s - max(0, 1 / pii - 0.5); tau <- tau_of(s)
  c(s = s, pi = pii, sprime = sp, tau = tau,
    c_min = max(tau / 2, (2 - tau) / (4 * sp)), piso = (2 - tau) / 4, piso2 = s / (2 * s + 1))
}, numeric(7)))
print(round(win, 4))
chk(max(abs(win[, "piso"] - win[, "piso2"])) < 1e-12,
    "com tau = (s+1/2)^{-1}, o piso (2-tau)/4 é s/(2s+1)")
chk(all((win[, "c_min"] < 1) == (win[, "sprime"] > win[, "piso"])),
    "a janela max{tau/2,(2-tau)/(4s')} <= c < 1 é não vazia exatamente quando s' > (2-tau)/4")

## ===========================================================================
cat("\nPARTE II. O regime J_n e a condição empírica de E1.4\n")
## ===========================================================================

J_rule <- function(n, sp) max(1L, as.integer(ceiling(log2((n / log(n))^(1 / (2 * sp + 1))))))
ns <- c(200L, 400L, 800L, 1600L, 3200L, 6400L, 12800L)
sp_dense <- 0.25                                  # pi = 2, s = 1/4 (logo s' = s)
E14 <- function(n, tJ) p * q * tJ * log(p * q * tJ) / n
tabII <- t(vapply(ns, function(n) {
  Jn <- J_rule(n, sp_dense); NJ <- 2L^Jn - 1L
  tJc <- (n / log(n))^(1 / (2 * sp_dense + 1))    # versão contínua da regra
  c(n = n, J_n = Jn, d = p * q * NJ, E14_cont = E14(n, tJc), E14_Jn = E14(n, 2^Jn))
}, numeric(5)))
print(round(tabII, 4))
# a regra contínua é monótona; o arredondamento de J_n ao inteiro não é, e é a
# única razão de E14_Jn oscilar
big <- 10^c(4, 6, 9, 12, 16)
E14_big <- vapply(big, function(n) E14(n, (n / log(n))^(1 / (2 * sp_dense + 1))), 1)
cat(sprintf("  extrapolação da regra contínua (n = 1e4, 1e6, 1e9, 1e12, 1e16): %s\n",
            paste(sprintf("%.4f", E14_big), collapse = " ")))
chk(all(diff(tabII[, "E14_cont"]) < 0) && all(diff(E14_big) < 0) && min(E14_big) < 0.01,
    "p q 2^J log(p q 2^J)/n decresce e tende a zero na regra contínua (condição de E1.4)")
chk(max(tabII[, "E14_Jn"] / tabII[, "E14_cont"]) < 2.5,
    "o arredondamento de J_n ao inteiro custa no máximo um fator 2 nessa condição")

draw_X <- function(n) cbind(1, 0.5 + runif(n, -1, 1))
lmin_at <- function(n, Jn, R = 20L) {
  vapply(seq_len(R), function(i) {
    U <- matrix(runif(n * q), n, q); X <- draw_X(n)
    Wl <- lapply(seq_len(q), function(m) psi_block(U[, m], Jn))
    Z <- design_from(X, Wl, Jn)
    rp <- resid_pen(Z)
    c(lmin = min(eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]]),
      gtil = min(eigen(crossprod(rp[["Bt"]]) / n, symmetric = TRUE, only.values = TRUE)[["values"]]))
  }, c(0, 0))
}
tabII2 <- t(vapply(ns, function(n) {
  e <- lmin_at(n, J_rule(n, sp_dense))
  c(n = n, J_n = J_rule(n, sp_dense), lmin = median(e["lmin", ]), gtil = median(e["gtil", ]),
    frac = mean(e["lmin", ] >= gam / 2), schur = mean(e["gtil", ] >= e["lmin", ] - 1e-10))
}, numeric(6)))
print(round(tabII2, 4))
chk(all(tabII2[, "schur"] == 1), "gamma_til >= lambda_min(Sigma_hat) em todas as réplicas e todos os n")
cat(sprintf("  fração das réplicas com lambda_min(Sigma_hat) >= gamma/2 = %.3f: %s\n",
            gam / 2, paste(sprintf("%.2f", tabII2[, "frac"]), collapse = " ")))
# a cadeia de E1.4 é assintótica: a cota gamma/2 só é atingida no maior n desta
# grade. O que se confere aqui é a MONOTONIA em n a J fixo (é 2^J/n que governa)
# e a convergência a lambda_min(Sigma) = kappa_1 c_U com J fixo.
mono <- all(vapply(unique(tabII2[, "J_n"]), function(Jv) {
  ii <- which(tabII2[, "J_n"] == Jv)
  length(ii) < 2 || all(diff(tabII2[ii, "lmin"]) > 0)
}, TRUE))
chk(mono, "a J_n fixo, lambda_min(Sigma_hat) cresce com n (é 2^{J}/n que governa)")
chk(tabII2[nrow(tabII2), "lmin"] >= gam / 2,
    "no maior n da grade a mediana de lambda_min(Sigma_hat) já passa de kappa_1 c_U / 2")
fix4 <- t(vapply(c(400L, 1600L, 6400L, 12800L), function(n) {
  e <- lmin_at(n, 4L, R = 10L)
  c(n = n, J = 4, lmin = median(e["lmin", ]), alvo = gam)
}, numeric(4)))
print(round(fix4, 4))
chk(all(diff(fix4[, "lmin"]) > 0) && fix4[nrow(fix4), "lmin"] > 0.8 * gam,
    "com J fixo, lambda_min(Sigma_hat) converge para lambda_min(Sigma) = kappa_1 c_U")

## ===========================================================================
cat("\nPARTE III+IV. pi = 2, s = s' = 1/4: predição (Teorema 2) e componentes (Corolário 4)\n")
## ===========================================================================

Jmax <- 10L                     # g é a soma finita dos níveis j < Jmax: satisfaz a
                                # Hipótese de Besov para todo j e dá viés exato
s_d <- 0.25; pi_d <- 2; Cg_d <- 0.3
theta_full_d <- theta_besov(s_d, pi_d, Cg_d, Jmax, p * q, spread = "random", seed = 23L)
c_true <- c(1, -0.5)
lvJ <- lev_of(Jmax)

sim_data <- function(n, theta_full) {
  U <- matrix(runif(n * q), n, q); X <- draw_X(n)
  Wl <- lapply(seq_len(q), function(m) psi_block(U[, m], Jmax))
  Zf <- design_from(X, Wl, Jmax)
  f <- as.numeric(Zf %*% c(c_true, theta_full))
  list(X = X, Wl = Wl, f = f, y = f + rnorm(n, sd = sig))
}
keep_idx <- function(Jl) {                       # recorta o vetor completo no sieve J
  NJ <- 2L^Jl - 1L
  as.integer(unlist(lapply(seq_len(p * q), function(bb)
    (bb - 1L) * length(lvJ) + seq_len(NJ))))
}
fit_at <- function(dat, Jl, theta_full) {
  n <- length(dat[["f"]]); NJ <- 2L^Jl - 1L; d <- p * q * NJ
  Z <- design_from(dat[["X"]], dat[["Wl"]], Jl)
  rp <- resid_pen(Z)
  smax <- sqrt(max(colSums(Z[, -seq_len(p), drop = FALSE]^2) / n))
  lam <- 2 * sig * smax * sqrt(2 * log(2 * d / alpha) / n)
  ki <- keep_idx(Jl)
  th_or <- theta_full[ki]
  fJ <- as.numeric(Z %*% c(c_true, th_or))
  fit <- fit_resid(rp, dat[["y"]], lam)
  v <- fit[["theta"]] - th_or
  gt <- min(eigen(crossprod(rp[["Bt"]]) / n, symmetric = TRUE, only.values = TRUE)[["values"]])
  list(J = Jl, d = d, lambda = lam, gtil = gt, bias_n = mean((dat[["f"]] - fJ)^2),
       pred = mean((as.numeric(Z %*% c(fit[["c"]], fit[["theta"]])) - dat[["f"]])^2),
       l2sq = sum(v^2), Peps = mean(rp[["proj"]](dat[["y"]] - dat[["f"]])^2),
       tail2 = sum(theta_full[-ki]^2), l1_or = sum(abs(th_or)),
       theta = fit[["theta"]], th_or = th_or, rp = rp)
}

R3 <- 15L
tabIII <- t(vapply(ns, function(n) {
  Jn <- J_rule(n, sp_dense)
  r <- lapply(seq_len(R3), function(i) fit_at(sim_data(n, theta_full_d), Jn, theta_full_d))
  pr <- vapply(r, function(z) z[["pred"]], 1)
  bd <- vapply(r, function(z) 3 * (64 * z[["lambda"]]^2 * sum(z[["th_or"]] != 0) / z[["gtil"]] +
                                   20 * z[["bias_n"]] + z[["Peps"]]), 1)
  l2 <- vapply(r, function(z) z[["l2sq"]] + z[["tail2"]], 1)
  bl2 <- vapply(r, function(z) 2 * (64 * z[["lambda"]]^2 * sum(z[["th_or"]] != 0) / z[["gtil"]]^2 +
                                    16 * z[["bias_n"]] / z[["gtil"]]) + 2 * z[["tail2"]], 1)
  alvo <- (Jn / n)^(2 * sp_dense / (2 * sp_dense + 1))
  c(n = n, J_n = Jn, d = p * q * (2^Jn - 1L), pred = median(pr),
    vies = median(vapply(r, function(z) z[["bias_n"]], 1)),
    lento = median(vapply(r, function(z) 6 * z[["lambda"]] * z[["l1_or"]], 1)),
    alvo = alvo, razao = median(pr) / alvo, viola = mean(pr > bd),
    L2 = median(l2), razao_L2 = median(l2) / alvo, viola_L2 = mean(l2 > bl2))
}, numeric(12)))
print(signif(tabIII, 4))
chk(all(tabIII[, "viola"] == 0),
    "o erro de predição fica sob a cota do Teorema 2 em todas as réplicas e todos os n")
chk(all(tabIII[, "viola_L2"] == 0), "a cota do Corolário 4 vale em todas as réplicas")
for (nm in c("razao", "razao_L2")) {
  rz <- tabIII[, nm]
  cat(sprintf("  %s / (J_n/n)^{2s'/(2s'+1)}: %s (máx/mín = %.2f)\n",
              if (nm == "razao") "||f_hat - f||_n^2" else "sum ||g_hat - g||^2_{L2}",
              paste(sprintf("%.3f", rz), collapse = " "), max(rz) / min(rz)))
  chk(max(rz) / min(rz) < 2, "a razão ao alvo é estável ao variar n por um fator 64")
}
sl <- coef(lm(log(tabIII[, "pred"]) ~ log(tabIII[, "alvo"])))[2]
cat(sprintf("  inclinação de log(erro) contra log do alvo: %.3f (alvo 1)\n", sl))
chk(abs(sl - 1) < 0.2, "o erro realizado segue a taxa (J_n/n)^{2s'/(2s'+1)} do Teorema 2")

cat("\n  J_n contra o J que minimiza o erro realizado (diagnóstico do pré-assintótico)\n")
tabIIIb <- t(vapply(ns[ns <= 6400L], function(n) {
  Jn <- J_rule(n, sp_dense)
  Js <- seq.int(max(1L, Jn - 2L), min(Jmax - 1L, Jn + 1L))
  dat <- lapply(seq_len(8L), function(i) sim_data(n, theta_full_d))
  pr <- vapply(Js, function(Jl)
    median(vapply(dat, function(dd) fit_at(dd, Jl, theta_full_d)[["pred"]], 1)), 1)
  c(n = n, J_n = Jn, J_opt = Js[which.min(pr)], erro_Jn = pr[Js == Jn], erro_opt = min(pr))
}, numeric(5)))
print(signif(tabIIIb, 4))
cat(sprintf("  J_opt - J_n: %s\n",
            paste(tabIIIb[, "J_opt"] - tabIIIb[, "J_n"], collapse = " ")))
chk(all(tabIIIb[, "erro_Jn"] / tabIIIb[, "erro_opt"] < 1.5),
    "o erro em J_n fica a menos de 50% acima do melhor J do intervalo varrido")

## ===========================================================================
cat("\nPARTE V. Comprimível (pi = 1): Lema 8 para todo s e a dimensão efetiva\n")
## ===========================================================================

s_c <- 1.2; pi_c <- 1; Cg_c <- 0.6
sp_c <- s_c - (1 / pi_c - 0.5)                 # s' = 0.7
tau_c <- tau_of(s_c)                           # 1/1.7 = 0.588
theta_full_c <- theta_besov(s_c, pi_c, Cg_c, Jmax, p * q, spread = "uniform", seed = 31L)
cat(sprintf("  s = %.1f, pi = %d: s' = %.3f, tau = %.3f; expoente denso 2s'/(2s'+1) = %.3f,\n",
            s_c, pi_c, sp_c, tau_c, 2 * sp_c / (2 * sp_c + 1)))
cat(sprintf("  expoente comprimível (2-tau)/2 = %.3f; piso s' > (2-tau)/4 = %.3f: %s\n",
            (2 - tau_c) / 2, (2 - tau_c) / 4, if (sp_c > (2 - tau_c) / 4) "vale" else "NÃO vale"))
chk(sp_c > (2 - tau_c) / 4, "o cenário comprimível satisfaz o piso de regularidade do Corolário 5")

c_comp <- max(tau_c / 2, (2 - tau_c) / (4 * sp_c))     # a ponta inferior da janela
J_comp <- function(n) max(1L, as.integer(ceiling(c_comp * log2(n))))
R5 <- 12L
tabV <- t(vapply(ns, function(n) {
  Jn <- min(J_comp(n), Jmax - 1L)
  out <- t(vapply(seq_len(R5), function(i) {
    z <- fit_at(sim_data(n, theta_full_c), Jn, theta_full_c)
    th_or <- z[["th_or"]]; srt <- order(abs(th_or), decreasing = TRUE)
    Bt <- z[["rp"]][["Bt"]]; Bm <- z[["rp"]][["B"]]
    ss <- unique(pmin(z[["d"]], c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L, 256L, z[["d"]])))
    # Lema 8 com o comparador theta-barra = melhor aproximação com s termos
    viola <- 0; hs <- numeric(length(ss))
    for (ii in seq_along(ss)) {
      tb_ <- numeric(length(th_or)); idx <- srt[seq_len(ss[ii])]; tb_[idx] <- th_or[idx]
      bb <- mean((as.numeric(Bm %*% (th_or - tb_)))^2)
      bnbar <- (sqrt(z[["bias_n"]]) + sqrt(bb))^2     # (||b||_n + ||B(theta*-tb)||_n)^2
      lhs <- mean((Bt %*% (z[["theta"]] - tb_))^2)
      rhs <- 64 * z[["lambda"]]^2 * ss[ii] / z[["gtil"]] + 16 * bnbar
      if (lhs > rhs) viola <- viola + 1
      hs[ii] <- rhs
    }
    c(viola = viola, sstar = ss[which.min(hs)], hmin = min(hs), hd = hs[length(hs)],
      pred = z[["pred"]], lambda = z[["lambda"]])
  }, numeric(6)))
  c(n = n, J_n = Jn, d = p * q * (2^Jn - 1L), viola = sum(out[, "viola"]),
    sstar = median(out[, "sstar"]), hmin = median(out[, "hmin"]),
    ganho = median(out[, "hd"] / out[, "hmin"]), lambda = median(out[, "lambda"]),
    pred = median(out[, "pred"]),
    razao_c = median(out[, "pred"]) / (Jn / n)^((2 - tau_c) / 2),
    razao_d = median(out[, "pred"]) / (Jn / n)^(2 * sp_c / (2 * sp_c + 1)))
}, numeric(11)))
print(signif(tabV, 4))
chk(all(tabV[, "viola"] == 0),
    "o Lema 8 (comparador arbitrário) vale para TODO s de aproximação, em todas as réplicas")
chk(all(tabV[, "sstar"] < tabV[, "d"]) && min(tabV[, "ganho"]) > 5,
    "o minimizador em s da cota é interior e bate por larga margem a leitura densa s = d")
sl5 <- coef(lm(log(tabV[, "sstar"]) ~ log(tabV[, "n"])))[2]
cat(sprintf("  inclinação log-log de s*_n em n: %.3f (alvo tau/2 = %.3f)\n", sl5, tau_c / 2))
chk(abs(sl5 - tau_c / 2) < 0.25, "a dimensão efetiva s*_n escala como n^{tau/2}")
sl5b <- coef(lm(log(tabV[, "hmin"]) ~ log(tabV[, "lambda"])))[2]
cat(sprintf("  inclinação log-log de h(s*) em lambda: %.3f (alvo 2 - tau = %.3f)\n", sl5b, 2 - tau_c))
chk(abs(sl5b - (2 - tau_c)) < 0.3, "o valor mínimo da cota escala como lambda^{2-tau}")
cat(sprintf("  erro realizado / (J_n/n)^{(2-tau)/2}: %s (máx/mín = %.2f)\n",
            paste(sprintf("%.2f", tabV[, "razao_c"]), collapse = " "),
            max(tabV[, "razao_c"]) / min(tabV[, "razao_c"])))
cat(sprintf("  erro realizado / (J_n/n)^{2s'/(2s'+1)}: %s (máx/mín = %.2f)\n",
            paste(sprintf("%.2f", tabV[, "razao_d"]), collapse = " "),
            max(tabV[, "razao_d"]) / min(tabV[, "razao_d"])))
cat("  (os dois expoentes distam 0.12; com n variando por um fator 64 eles ainda não se\n")
cat("   separam empiricamente, e o que a conferência decide é a cota, não o realizado)\n")

cat("\n")
if (ok) cat("OK\n") else stop("E1.6: conferência numérica FALHOU (ver linhas acima)")
