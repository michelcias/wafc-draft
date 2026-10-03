# E1.12. Conferência numérica de derivations/08-blocos.tex (a teoria em
# blocos do estimador de D44: block LASSO na forma balanceada, níveis c_l
# livres, pesos de razão limitada) e da cota do estimador limiarizado
# acrescentada a derivations/06-selecao-limiar.tex (D45). Roda ANTES da prova;
# imprime OK ou falha com stop(). Notação: docs/notacao.md (congelada em
# E1.1); os símbolos da teoria em blocos são os de E1.11 (08a, §1.4).
#
# Estende as Partes B, C e E de check/08a-blocos.R a J_n e à taxa, como
# check/05-taxas.R faz para o LASSO. Não chama nada de wafc/R/: o ajuste em
# blocos é o FISTA de check/08a-blocos.R e o LASSO é o glmnet na rota
# perfilada de check/05-taxas.R.
#
# Partes:
#   I.   determinístico: a Proposição 7 (a forma balanceada: tamanhos dos
#        pedaços, |G| <= pq(1 + 2^J/b), rho^2 <= (2b - 1)/(b - 1) <= 3 com
#        sqrt(|G|)); o Lema 15 (risco ideal por pedaços sob Besov, com um eta
#        pelo pedaço grosso), em três formas de sequência e pi até infinito;
#        a aritmética do Teorema 4 e do Corolário 11 (balanço, janela, ganho,
#        termos de ordem menor) e o que o Teorema 2 de K&P não cobre.
#   II.  o regime: J_n das duas regras, a condição empírica de E1.4, os
#        eventos sigma_max^2 <= 2 B_X^2 C_U e lambda_max(Sigma_hat) <= Lambda,
#        a forma fechada do Lema 14 (exata <= fechada <= determinística) e a
#        cobertura do evento T_{G,w}.
#   III. denso (pi = 2, s = s' = 1/2): o Teorema 3 com o comparador theta*,
#        o Teorema 4 e o Corolário 10 contra n^{-2s'/(2s'+1)}.
#   IV.  comprimível (pi = 1, s = 1,2): o Teorema 3 com comparadores
#        truncados, o Corolário 9 (a forma de risco ideal), o Lema 15 no
#        theta* do desenho e o Corolário 11 contra
#        n^{-2s/(2s+1)} (log n)^{(2/pi - 1)_+/(2s+1)}.
#   V.   a variante branca (Observação do 08): norma ||Q_G^{1/2} theta_G|| com
#        Q_G a Gram centrada do pedaço, a do grpreg; calibração pivotal e o
#        Teorema 3 com q_max W no lugar de W.
#   VI.  o limiarizado (06, Lema 16, Corolários 12 e 13): o lema em
#        configurações aleatórias, a passagem pela norma de predição com o
#        autovalor, e os dois ajustes (LASSO e blocos) em simulação.
#
# Nos scripts de conferência, elemento de lista se acessa com [[ ]] e nome
# completo (instrucoes.md, §5).
#
# Dependências: WaveBased (wbasis, wtable), glmnet.
# Tempo nesta máquina: ver o parágrafo da conferência no 08-blocos.tex.

suppressPackageStartupMessages({
  library(WaveBased)
  library(glmnet)
})
set.seed(20261003)
t_start <- proc.time()[["elapsed"]]

p  <- 2L                        # X_1 = 1, X_2 = 1/2 + Unif(-1,1) (cenário A de E1.4)
q  <- 2L
fs <- 8L                        # filter.size (Daublets, 4 momentos nulos)
sig <- 0.5
alpha <- 0.05
B_X <- 1.5
C_U <- 1
gam <- 0.25                     # kappa_1 c_U (U uniforme: c_U = C_U = 1)
kap2 <- 4 / 3                   # lambda_max(E XX'): autovalores de [[1, 1/2], [1/2, 7/12]]
Lam <- kap2 * C_U + gam / 2     # Lambda de E1.6 (Lema 9(iii))
Jmax <- 10L                     # a verdade é a soma finita dos níveis j < Jmax
tb <- wtable(filter.size = fs)

ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "ok" else "FALHA", msg))
  if (!isTRUE(cond)) ok <<- FALSE
  invisible(cond)
}

## ---- desenho (construção de check/05-taxas.R) ------------------------------

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
design_from <- function(X, Wlist, Jl) {            # Wlist[[m]] avaliada em Jmax
  NJ <- 2L^Jl - 1L
  P <- cbind(1, do.call(cbind, lapply(Wlist, function(W) W[, seq_len(NJ), drop = FALSE])))
  kron_rows(X, P)[, perm_D12(ncol(X), length(Wlist), NJ), drop = FALSE]
}
draw_X <- function(n) cbind(1, 0.5 + runif(n, -1, 1))
resid_pen <- function(Z) {
  A <- Z[, seq_len(p), drop = FALSE]
  Bm <- Z[, -seq_len(p), drop = FALSE]
  qa <- qr(A)
  list(A = A, B = Bm, Bt = Bm - qr.fitted(qa, Bm), proj = function(v) qr.fitted(qa, v))
}
lev_of <- function(Jl) rep(0:(Jl - 1L), 2^(0:(Jl - 1L)))
keep_idx <- function(Jl) {                          # recorta o vetor de Jmax níveis em J
  NJ <- 2L^Jl - 1L
  NJm <- 2L^Jmax - 1L
  as.integer(unlist(lapply(seq_len(p * q), function(bb) (bb - 1L) * NJm + seq_len(NJ))))
}
blk_full <- rep(seq_len(p * q), each = 2L^Jmax - 1L)

## ---- sequências de Besov (check/05-taxas.R) --------------------------------

theta_besov_blk <- function(s, pii, Cg, Jl, spread = c("random", "uniform"), seed = 1L) {
  spread <- match.arg(spread)
  set.seed(seed)
  lv <- lev_of(Jl)
  out <- numeric(length(lv))
  for (jj in 0:(Jl - 1L)) {
    idx <- which(lv == jj)
    a <- if (spread == "uniform") sample(c(-1, 1), length(idx), TRUE) else rnorm(length(idx))
    npi <- if (is.infinite(pii)) max(abs(a)) else sum(abs(a)^pii)^(1 / pii)
    out[idx] <- Cg * 2^(-jj * (s + 0.5 - 1 / pii)) * a / npi
  }
  out
}
theta_besov <- function(s, pii, Cg, spread, seed) {
  unlist(lapply(seq_len(p * q), function(bb) theta_besov_blk(s, pii, Cg, Jmax, spread, seed + bb)))
}

## ---- a forma balanceada (wafc_kp_groups(balanced = TRUE), E2.5e) ----------

bal_sizes <- function(Jl, b) {
  levs <- 0:(Jl - 1L)
  co <- levs[2^levs < b]
  out <- if (length(co) > 0L) sum(2L^co) else integer(0)
  for (jj in setdiff(levs, co)) {
    nj <- 2L^jj
    k <- nj %/% b
    sz <- rep(b, k)
    sz[k] <- sz[k] + nj - k * b
    out <- c(out, sz)
  }
  as.integer(out)
}
kp_groups_bal <- function(Jl, b, nblocks = p * q) {
  sz <- bal_sizes(Jl, b)
  one <- rep(seq_along(sz), sz)
  as.integer(unlist(lapply(seq_len(nblocks), function(bb) one + (bb - 1L) * length(sz))))
}

## ---- solver: group LASSO no problema perfilado (check/08a-blocos.R) --------

kkt_gl <- function(G, h, grp, w, lam, x) {
  r <- as.numeric(h - G %*% x)
  nx <- sqrt(as.numeric(rowsum(x^2, grp)))
  nr <- sqrt(as.numeric(rowsum(r^2, grp)))
  act <- nx > 0
  dev_act <- 0
  if (any(act)) {
    u <- x / pmax(nx, 1e-300)[grp]
    dd <- r - lam * w[grp] * u
    dev_act <- max(sqrt(as.numeric(rowsum(dd^2, grp)))[act]) / lam
  }
  dev_ina <- if (any(!act)) max(0, max(nr[!act] - lam * w[!act])) / lam else 0
  c(ativo = dev_act, inativo = dev_ina)
}
gl_fit <- function(G, h, grp, w, lam, Lmax, x0 = NULL, tol = 1e-9, maxit = 100000L) {
  d <- length(h)
  L <- 2 * Lmax
  thr <- 2 * lam * w / L
  x <- if (is.null(x0)) numeric(d) else x0
  y <- x
  tk <- 1
  for (it in seq_len(maxit)) {
    z <- y + 2 * as.numeric(h - G %*% y) / L
    nz <- sqrt(as.numeric(rowsum(z^2, grp)))
    xn <- z * pmax(0, 1 - thr / pmax(nz, 1e-300))[grp]
    if (sum((y - xn) * (xn - x)) > 0) {
      tk <- 1
      y <- xn
    } else {
      tn <- (1 + sqrt(1 + 4 * tk^2)) / 2
      y <- xn + ((tk - 1) / tn) * (xn - x)
      tk <- tn
    }
    x <- xn
    if (it %% 20L == 0L && max(kkt_gl(G, h, grp, w, lam, x)) < tol) break
  }
  x
}
coef_c <- function(rp, y, th) drop(qr.coef(qr(rp[["A"]]), y - rp[["B"]] %*% th))

## ---- quantidades da teoria --------------------------------------------------

# Lema 14: lambda_{0,G} exato (HKZ), forma fechada (sigma_max, lambda_max(Sigma_hat))
# e determinística (B_X sqrt(2 C_U), Lambda).
lam0_exact <- function(n, M, trv, opv) sig / sqrt(n) * (sqrt(trv) + sqrt(2 * opv * log(M / alpha)))
lam0_cf <- function(n, M, gsize, smax, lmaxS) sig / sqrt(n) * (smax * sqrt(gsize) + sqrt(2 * lmaxS * log(M / alpha)))
lam0_det <- function(n, M, gsize) sig / sqrt(n) * (B_X * sqrt(2 * C_U * gsize) + sqrt(2 * Lam * log(M / alpha)))
lam_w <- function(l0, w) 2 * max(l0 / w)
Wset <- function(th, grp, w) {
  act <- which(sqrt(as.numeric(rowsum(th^2, grp))) > 0)
  sum(w[act]^2)
}
Rideal <- function(th, grp, eta) sum(pmin(as.numeric(rowsum(th^2, grp)), eta))
tau_of <- function(s) 1 / (s + 0.5)
sp_of <- function(s, pii) s - max(0, 1 / pii - 0.5)
A15 <- function(s, pii) if (pii <= 2) 2 + 1 / (1 - 2^(1 - pii * (s + 0.5))) else 2 + 1 / (1 - 2^(-2 * s))
bound15 <- function(s, pii, Cg, b, eta, nblocks = 1L) {
  e <- if (is.infinite(pii)) 0 else 2 / pii
  nblocks * (A15(s, pii) * Cg^tau_of(s) * (eta / b)^(2 * s / (2 * s + 1)) *
               b^(max(0, e - 1) / (2 * s + 1)) + eta)
}

# Um conjunto de dados com o desenho avaliado em Jmax e a verdade em theta_full.
c_true <- c(1, -0.5)
sim_data <- function(n, theta_full) {
  U <- matrix(runif(n * q), n, q)
  X <- draw_X(n)
  Wl <- lapply(seq_len(q), function(m) psi_block(U[, m], Jmax))
  NJm <- 2L^Jmax - 1L
  f <- numeric(n)                                   # f = sum_l X_l {c_l + sum_m g_lm(U_m)}, blocos em D12
  for (l in seq_len(p)) {
    bl <- rep(c_true[l], n)
    for (m in seq_len(q)) {
      bb <- (l - 1L) * q + m
      bl <- bl + as.numeric(Wl[[m]] %*% theta_full[(bb - 1L) * NJm + seq_len(NJm)])
    }
    f <- f + X[, l] * bl
  }
  eps <- rnorm(n, sd = sig)
  list(X = X, Wl = Wl, f = f, eps = eps, y = f + eps)
}
# Tudo o que um ajuste em J precisa e que não depende de lambda.
prep_J <- function(dat, Jl, theta_full) {
  n <- length(dat[["f"]])
  Z <- design_from(dat[["X"]], dat[["Wl"]], Jl)
  rp <- resid_pen(Z)
  Bt <- rp[["Bt"]]
  G <- crossprod(Bt) / n
  evG <- eigen(G, symmetric = TRUE, only.values = TRUE)[["values"]]
  evS <- eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]]
  ki <- keep_idx(Jl)
  th_or <- theta_full[ki]
  fJ <- as.numeric(Z %*% c(c_true, th_or))
  yt <- as.numeric(dat[["y"]] - rp[["proj"]](dat[["y"]]))
  tail_blk <- as.numeric(rowsum(theta_full[-ki]^2, blk_full[-ki]))
  list(n = n, Z = Z, rp = rp, G = G, h = as.numeric(crossprod(Bt, yt)) / n,
       gtil = min(evG), Lmax = max(evG), lmin = min(evS), lmax = max(evS),
       smax = sqrt(max(colSums(rp[["B"]]^2) / n)), th_or = th_or, fJ = fJ,
       bias = mean((dat[["f"]] - fJ)^2), Peps = mean(rp[["proj"]](dat[["eps"]])^2),
       tail_blk = tail_blk, blk = rep(seq_len(p * q), each = 2L^Jl - 1L))
}
group_stats <- function(Bt, grp) {
  n <- nrow(Bt)
  M <- max(grp)
  tr <- numeric(M)
  op <- numeric(M)
  for (g in seq_len(M)) {
    ii <- which(grp == g)
    ev <- eigen(crossprod(Bt[, ii, drop = FALSE]) / n, symmetric = TRUE, only.values = TRUE)[["values"]]
    tr[g] <- sum(ev)
    op[g] <- max(ev)
  }
  list(tr = tr, op = op)
}

## ===========================================================================
cat("PARTE I. Determinístico: Proposição 7, Lema 15 e a aritmética das taxas\n")
## ===========================================================================

## ---- I.1 Proposição 7 ---------------------------------------------------------
P7 <- do.call(rbind, lapply(2:200, function(b) {
  do.call(rbind, lapply(1:14, function(Jl) {
    sz <- bal_sizes(Jl, b)
    levs <- 0:(Jl - 1L)
    jst <- max(levs[2^levs < b])                    # j*: o nível mais fino com 2^j < b
    has_fine <- any(2^levs >= b)
    coarse <- sz[1L]
    fine <- if (has_fine) sz[-1L] else integer(0)
    c(b = b, J = Jl,
      soma = sum(sz) == 2L^Jl - 1L,
      fino = if (has_fine) all(fine >= b & fine <= 2L * b - 1L) else TRUE,
      grosso = if (has_fine) coarse == 2L^(jst + 1L) - 1L && coarse >= b - 1L && coarse <= 2L * b - 3L
               else coarse == 2L^Jl - 1L,
      card = length(sz) <= 1 + 2^Jl / b,
      rho2 = max(sz) / min(sz) <= (2 * b - 1) / (b - 1) + 1e-12 && (2 * b - 1) / (b - 1) <= 3)
  }))
}))
chk(all(P7[, c("soma", "fino", "grosso", "card", "rho2")] == 1),
    "Proposição 7, b = 2..200, J = 1..14: pedaços finos entre b e 2b - 1; grosso com 2^(j*+1) - 1 colunas entre b - 1 e 2b - 3; |G| por bloco <= 1 + 2^J/b; max|G|/min|G| <= (2b - 1)/(b - 1) <= 3")
cat(sprintf("  rho^2 dos pesos sqrt(|G|) nos b_n de n = 250, 500, 1000, 4000 (J = 8): %s\n",
            paste(vapply(c(250, 500, 1000, 4000), function(n) {
              sz <- bal_sizes(8L, as.integer(ceiling(log(n))))
              sprintf("%.3f", max(sz) / min(sz))
            }, ""), collapse = " ")))

## ---- I.2 Lema 15 ----------------------------------------------------------------
# Energias dos pedaços balanceados de um bloco (formas "uniforme" e
# "aleatoria", check/08a-blocos.R, Parte E4), e a "espiga": um pico de altura
# sqrt(eta) por pedaço fino, em tantos pedaços quanto o orçamento ell_pi do
# nível permite, com o pedaço grosso contado como eta.
budget_of <- function(s, pii, Cg, jj) Cg * 2^(-jj * (s + 0.5 - (if (is.infinite(pii)) 0 else 1 / pii)))
energies_bal <- function(s, pii, Cg, Jl, b, shape, lev_vals = NULL) {
  levs <- 0:(Jl - 1L)
  co <- levs[2^levs < b]
  e_lev <- function(jj) {
    if (shape == "uniforme") {
      cj <- budget_of(s, pii, Cg, jj) / (if (is.infinite(pii)) 1 else (2^jj)^(1 / pii))
      return(rep(cj^2, 2^jj))
    }
    lev_vals[[jj + 1L]]^2
  }
  fine <- unlist(lapply(setdiff(levs, co), function(jj) {
    nj <- 2L^jj
    k <- nj %/% b
    sz <- rep(b, k)
    sz[k] <- sz[k] + nj - k * b
    as.numeric(rowsum(e_lev(jj), rep(seq_len(k), sz)))
  }))
  list(coarse = if (length(co) > 0L) sum(unlist(lapply(co, e_lev))) else 0, fine = fine)
}
R_bal <- function(s, pii, Cg, Jl, b, eta, shape, en = NULL) {
  levs <- 0:(Jl - 1L)
  co <- levs[2^levs < b]
  if (shape == "espiga") {
    rf <- vapply(setdiff(levs, co), function(jj) {
      bud <- budget_of(s, pii, Cg, jj)
      kmax <- if (is.infinite(pii)) (if (bud >= sqrt(eta)) Inf else 0) else floor((bud / sqrt(eta))^pii)
      min((2^jj) %/% b, kmax) * eta
    }, 1)
    return((length(co) > 0L) * eta + sum(rf))
  }
  min(en[["coarse"]], eta) + sum(pmin(en[["fine"]], eta))
}
grid15 <- expand.grid(s = c(0.3, 0.8, 1.5, 3), pi = c(1, 1.5, 2, 4, Inf))
grid15 <- grid15[grid15[["pi"]] * (grid15[["s"]] + 0.5) > 1, ]
JC <- 20L
viol15 <- 0
maxr15 <- 0
n15 <- 0
for (i in seq_len(nrow(grid15))) {
  s <- grid15[i, "s"]
  pii <- grid15[i, "pi"]
  set.seed(200L + i)
  lv <- lapply(0:(JC - 1L), function(jj) {
    a <- rnorm(2^jj)
    npi <- if (is.infinite(pii)) max(abs(a)) else sum(abs(a)^pii)^(1 / pii)
    budget_of(s, pii, 1, jj) * a / npi
  })
  for (b in c(2L, 3L, 4L, 6L, 7L, 9L, 12L, 32L)) {
    en <- list(uniforme = energies_bal(s, pii, 1, JC, b, "uniforme"),
               aleatoria = energies_bal(s, pii, 1, JC, b, "aleatoria", lev_vals = lv))
    for (eta in 10^c(-9, -7, -5, -3, -1)) {
      B <- bound15(s, pii, 1, b, eta)
      for (shape in c("uniforme", "aleatoria", "espiga")) {
        R <- R_bal(s, pii, 1, JC, b, eta, shape, en = en[[shape]])
        maxr15 <- max(maxr15, R / B)
        n15 <- n15 + 1
        if (R > B * (1 + 1e-12)) viol15 <- viol15 + 1
      }
    }
  }
}
cat(sprintf("  Lema 15: %d sequências (%d pares (s, pi) com pi até infinito, b em {2,...,32}, eta em 1e-9..1e-1, J = %d): maior R/cota = %.3f\n",
            n15, nrow(grid15), JC, maxr15))
chk(viol15 == 0, "Lema 15: R_G(theta*; eta) <= A C_g^tau (eta/b)^{2s/(2s+1)} b^{(2/pi-1)_+/(2s+1)} + eta em todas as sequências")
chk(maxr15 > 0.5, "Lema 15: a cota não é frouxa por mais de um fator 2 na pior sequência construída")

## ---- I.3 a aritmética das taxas -------------------------------------------------
gridD <- expand.grid(s = c(0.3, 0.6, 1, 1.5, 3), pi = c(1, 1.5, 2, 4, Inf))
gridD <- gridD[gridD[["pi"]] * (gridD[["s"]] + 0.5) > 1, ]
ar <- t(vapply(seq_len(nrow(gridD)), function(i) {
  s <- gridD[i, "s"]
  pii <- gridD[i, "pi"]
  sp <- sp_of(s, pii)
  tau <- tau_of(s)
  e2 <- if (is.infinite(pii)) 0 else 2 / pii
  cmin <- s / ((2 * s + 1) * sp)
  # Corolário 11: com eta/b = kappa/n, a parcela principal é n^{-2s/(2s+1)} b^{(2/pi-1)_+/(2s+1)}
  bn <- log(1e6)
  lhs <- (bn * 1e-6 / bn)^(2 * s / (2 * s + 1)) * bn^(max(0, e2 - 1) / (2 * s + 1))
  rhs <- 1e6^(-2 * s / (2 * s + 1)) * log(1e6)^(max(0, e2 - 1) / (2 * s + 1))
  # termos de ordem menor contra a taxa do Corolário 11: (log n)/n e 1/n
  rate <- function(n) n^(-2 * s / (2 * s + 1)) * log(n)^(max(0, e2 - 1) / (2 * s + 1))
  lo <- vapply(10^c(4, 8, 12, 16, 24), function(n) (log(n) / n) / rate(n), 1)
  # Teorema 4: 2^J = n^{1/(2s'+1)} iguala 2^J/n e 2^{-2Js'}
  n <- 1e6
  tJ <- n^(1 / (2 * sp + 1))
  c(janela_tau = cmin >= tau / 2 - 1e-12,
    janela_vazia = (cmin < 1) == (sp > s / (2 * s + 1)),
    ident = abs(lhs - rhs) / rhs,
    ganho = abs((2 * s / (2 * s + 1) - max(0, e2 - 1) / (2 * s + 1)) - 2 * sp / (2 * s + 1)),
    menor = all(diff(lo) < 0) && lo[length(lo)] < 0.05,
    thm4 = abs(tJ / n - tJ^(-2 * sp)) / (tJ / n),
    kp = sp > 0.5, nosso = sp > s / (2 * s + 1))
}, numeric(8)))
chk(all(ar[, "janela_tau"] == 1),
    "a janela do Corolário 11 é a do Corolário 5: s/((2s+1)s') >= tau/2 em todo (s, pi)")
chk(all(ar[, "janela_vazia"] == 1), "a janela s/((2s+1)s') <= c < 1 é não vazia sse s' > s/(2s+1)")
chk(max(ar[, c("ident", "ganho", "thm4")]) < 1e-10,
    "Corolário 11: (eta/b)^{2s/(2s+1)} b^{(2/pi-1)_+/(2s+1)} = n^{-2s/(2s+1)} b^{...} com eta/b = 1/n; ganho sobre o Corolário 5 = (log n)^{2s'/(2s+1)}; Teorema 4: 2^J/n = 2^{-2Js'} em 2^J = n^{1/(2s'+1)}")
chk(all(ar[, "menor"] == 1), "(log n)/n é de ordem menor que a taxa do Corolário 11 em todo (s, pi) (razão decrescente, < 0,05 em n = 1e24)")
fora_kp <- gridD[ar[, "nosso"] == 1 & ar[, "kp"] == 0, ]
cat(sprintf("  pares (s, pi) da grade cobertos pelo Corolário 11 e não pelo Teorema 2 de K&P (que pede r* = s' > 1/(2 varsigma) > 1/2): %s\n",
            paste(sprintf("(%.1f, %s)", fora_kp[["s"]], format(fora_kp[["pi"]])), collapse = " ")))
chk(sp_of(1, 1) == 0.5 && 0.5 > 1 / 3,
    "o cenário não homogêneo de D27 (s' = 1/2: pi = 1, s = 1) está no Corolário 11 (s' > s/(2s+1) = 1/3) e fora do Teorema 2 de K&P (r* > 1/2)")

## ---- I.4 os dois expoentes do Corolário 11 no risco ideal --------------------
# A taxa do Corolário 11 é (eta/b)^{2s/(2s+1)} b^{(2/pi-1)_+/(2s+1)} com
# eta/b ~ 1/n e b = b_n ~ log n. Os dois expoentes são lidos separadamente nos
# pedaços finos da forma balanceada, num bloco com 60 níveis (sem truncamento,
# que só diminuiria R_G, e sem o pedaço grosso, que é o termo eta de ordem
# menor): o de n, pela inclinação de log R contra log(eta/b) com b = 16 fixo; o
# do logaritmo, pela inclinação de log R contra log b com eta/b fixo e a
# transição entre cabeça e cauda perto do nível 30 (como na Parte C de
# check/08a-blocos.R). A "espiga" (desfavorável aos pedaços) deve pagar
# b^{(2/pi-1)_+/(2s+1)}; a "uniforme" (mesmo módulo no nível) não paga nada.
R_fine <- function(s, pii, b, eta, shape, Jl = 60L) {
  levs <- 0:(Jl - 1L)
  fine <- levs[2^levs >= b]
  sum(vapply(fine, function(jj) {
    bud <- budget_of(s, pii, 1, jj)
    k <- floor(2^jj / b)
    if (shape == "espiga") {
      kmax <- if (is.infinite(pii)) (if (bud >= sqrt(eta)) Inf else 0) else floor((bud / sqrt(eta))^pii)
      return(min(k, kmax) * eta)
    }
    cj2 <- (bud / (if (is.infinite(pii)) 1 else (2^jj)^(1 / pii)))^2
    r <- 2^jj - k * b
    (k - 1) * min(b * cj2, eta) + min((b + r) * cj2, eta)
  }, 1))
}
tab14 <- do.call(rbind, lapply(list(c(1.2, 1), c(1, 1), c(0.8, 1.5), c(1.5, 2), c(1, 4)), function(sp_) {
  s <- sp_[1]
  pii <- sp_[2]
  e2 <- max(0, 2 / pii - 1) / (2 * s + 1)
  # expoente de n: eta/b = 2^{-k(2s+1)}, k de 20 a 36 em quartos (a transição percorre os níveis 20 a 36)
  epsn <- 2^(-seq(20, 36, by = 0.25) * (2 * s + 1))
  incl_n <- vapply(c("espiga", "uniforme"), function(sh) {
    R <- vapply(epsn, function(e) R_fine(s, pii, 16, 16 * e, sh), 1)
    unname(coef(lm(log(R) ~ log(epsn)))[2])
  }, 1)
  # expoente do logaritmo: eta = b eps com eps fixo, transição perto do nível 30 em b = 32
  a <- pii * (s + 0.5)
  epsb <- if (pii <= 2) (32 * 2^(-30 * a))^(2 / pii) / 32 else 2^(-30 * (2 * s + 1))
  bs <- c(8, 16, 32, 64, 128)
  incl_b <- vapply(c("espiga", "uniforme"), function(sh) {
    R <- vapply(bs, function(bb) R_fine(s, pii, bb, bb * epsb, sh), 1)
    unname(coef(lm(log(R) ~ log(bs)))[2])
  }, 1)
  c(s = s, pi = pii, expo_n = 2 * s / (2 * s + 1), n_espiga = incl_n[["espiga"]], n_uniforme = incl_n[["uniforme"]],
    expo_log = e2, log_espiga = incl_b[["espiga"]], log_uniforme = incl_b[["uniforme"]])
}))
print(round(tab14, 3))
chk(max(abs(tab14[, c("n_espiga", "n_uniforme")] - tab14[, "expo_n"])) < 0.05,
    "o expoente de n do risco ideal é 2s/(2s+1) nas duas formas (inclinação em log(eta/b) a 0,05)")
lt2 <- tab14[, "pi"] < 2
chk(all(abs(tab14[lt2, "log_espiga"] - tab14[lt2, "expo_log"]) < 0.1) && all(tab14[!lt2, "log_espiga"] < 0.05),
    "o expoente do logaritmo é (2/pi-1)/(2s+1) na espiga em pi < 2 (a 0,1) e não positivo em pi >= 2")
chk(all(abs(tab14[, "log_uniforme"]) < 0.1),
    "a sequência espalhada no nível não paga o logaritmo residual (expoente de b abaixo de 0,1)")

## ===========================================================================
cat("\nPARTE II. O regime J_n, os eventos e a calibração ao longo dele\n")
## ===========================================================================

J_thm4 <- function(n, sp) max(1L, as.integer(ceiling(log2(n^(1 / (2 * sp + 1))))))
J_cor11 <- function(n, s, sp) max(1L, as.integer(ceiling(s / ((2 * s + 1) * sp) * log2(n))))
E14 <- function(n, Jl) p * q * 2^Jl * log(p * q * 2^Jl) / n
nsII <- c(250L, 500L, 1000L, 2000L, 4000L)
RII <- 40L
tabII <- t(vapply(nsII, function(n) {
  Jl <- J_thm4(n, 0.5)
  b <- as.integer(ceiling(log(n)))
  grp <- kp_groups_bal(Jl, b)
  M <- max(grp)
  gsize <- as.numeric(table(grp))
  ws <- sqrt(gsize)
  l0d <- lam0_det(n, M, gsize)
  lamD <- lam_w(l0d, ws)
  lam1D <- 2 * max(l0d)
  rr <- t(vapply(seq_len(RII), function(r) {
    U <- matrix(runif(n * q), n, q)
    X <- draw_X(n)
    Wl <- lapply(seq_len(q), function(m) psi_block(U[, m], Jl))
    Z <- design_from(X, Wl, Jl)
    rp <- resid_pen(Z)
    Bt <- rp[["Bt"]]
    evS <- eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]]
    smax <- sqrt(max(colSums(rp[["B"]]^2) / n))
    gs <- group_stats(Bt, grp)
    eps <- rnorm(n, sd = sig)
    nrm <- sqrt(as.numeric(rowsum((as.numeric(crossprod(Bt, eps)) / n)^2, grp)))
    l0e <- lam0_exact(n, M, gs[["tr"]], gs[["op"]])
    l0c <- lam0_cf(n, M, gsize, smax, max(evS))
    evs <- smax^2 <= 2 * B_X^2 * C_U
    evl <- max(evS) <= Lam
    c(evs = evs, evl = evl, evmin = min(evS) >= gam / 2,
      e_c = all(l0e <= l0c + 1e-12), lmax = max(evS),
      cobC = max(nrm / ws) <= lam_w(l0c, ws) / 2, cobD = max(nrm / ws) <= lamD / 2,
      gt_ls = min(eigen(crossprod(Bt) / n, symmetric = TRUE, only.values = TRUE)[["values"]]) >= min(evS) - 1e-10)
  }, numeric(8)))
  c(n = n, J = Jl, b = b, d = p * q * (2^Jl - 1), M = M, E14 = E14(n, Jl),
    E14c = p * q * n^(1 / 2) * log(p * q * n^(1 / 2)) / n,
    f_smax = mean(rr[, "evs"]), f_lmax = mean(rr[, "evl"]), f_lmin = mean(rr[, "evmin"]),
    e_c = mean(rr[, "e_c"]), lmax = median(rr[, "lmax"]),
    cobC = mean(rr[, "cobC"]), cobD = mean(rr[, "cobD"]), schur = mean(rr[, "gt_ls"]),
    nlam1_b = n * lam1D^2 / b, custo = max((lamD * ws)^2) / lam1D^2, rho2 = max(gsize) / min(gsize))
}, numeric(18)))
print(round(tabII, 4))
E14_big <- vapply(10^c(4, 6, 9, 12, 16), function(n) p * q * n^(1 / 2) * log(p * q * n^(1 / 2)) / n, 1)
chk(all(diff(tabII[, "E14c"]) < 0) && all(diff(E14_big) < 0) && min(E14_big) < 0.01 &&
      max(tabII[, "E14"] / tabII[, "E14c"]) < 2.5,
    "pq 2^J log(pq 2^J)/n decresce e vai a zero na regra contínua do Teorema 4 (condição de E1.4); o arredondamento de J_n custa menos que 2,5")
chk(all(tabII[, "f_smax"] == 1), "sigma_max^2 <= 2 B_X^2 C_U (Lema 6 de E1.5) em todas as réplicas: a primeira metade da forma determinística")
chk(all(tabII[, "e_c"] == 1), "Lema 14(ii): lambda_{0,G} exato <= forma fechada (sigma_max, lambda_max(Sigma_hat)) em toda réplica")
chk(all(tabII[, c("cobC", "cobD")] >= 1 - alpha), "a calibração do Lema 14 tem a probabilidade nominal com lambda fechado e com lambda_n^G determinístico")
chk(all(tabII[, "schur"] == 1), "gamma_til >= lambda_min(Sigma_hat) em todas as réplicas")
chk(all(tabII[, "custo"] <= tabII[, "rho2"] + 1e-12), "(lambda_n^G w_G)^2 <= rho^2 (lambda_{n,1}^G)^2 com os pesos sqrt(|G|)")
chk(max(tabII[, "nlam1_b"]) / min(tabII[, "nlam1_b"]) < 1.5, "n (lambda_{n,1}^G)^2 / b_n é limitado ao longo de n: lambda^2 ~ sigma^2 log n / n")
cat(sprintf("  lambda_max(Sigma_hat) mediano ao longo de J_n: %s, contra Lambda = %.3f e lambda_max(Sigma) = %.3f: o evento de E1.4 ainda não vale nesses n\n",
            paste(sprintf("%.2f", tabII[, "lmax"]), collapse = " "), Lam, kap2 * C_U))
# O evento {||Sigma_hat - Sigma|| < gamma/2} é assintótico (como em check/05-taxas.R,
# Parte II): a J fixo ele chega. J = 4 (d = 60), n de 1000 a 64000.
fixJ <- t(vapply(c(1000L, 4000L, 16000L, 64000L), function(n) {
  e <- vapply(seq_len(5L), function(r) {
    U <- matrix(runif(n * q), n, q)
    Z <- design_from(draw_X(n), lapply(seq_len(q), function(m) psi_block(U[, m], 4L)), 4L)
    ev <- eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]]
    c(max(ev), min(ev))
  }, c(0, 0))
  c(n = n, lmax = median(e[1, ]), lmin = median(e[2, ]))
}, numeric(3)))
print(round(fixJ, 4))
chk(all(diff(fixJ[, "lmax"]) < 0) && all(diff(fixJ[, "lmin"]) > 0) &&
      fixJ[nrow(fixJ), "lmax"] <= Lam && fixJ[nrow(fixJ), "lmin"] >= gam / 2,
    "a J fixo, lambda_max(Sigma_hat) desce e lambda_min(Sigma_hat) sobe com n, e em n = 64000 os dois estão no evento de E1.4 (<= Lambda, >= gamma/2)")

## ===========================================================================
cat("\nPARTE III. Denso, pi = 2, s = s' = 1/2: Teorema 3, Teorema 4 e Corolário 10\n")
## ===========================================================================

s_d <- 0.5
pi_d <- 2
Cg_d <- 0.3
th_full_d <- theta_besov(s_d, pi_d, Cg_d, "random", 40L)
nsIII <- c(250L, 500L, 1000L, 2000L, 4000L, 8000L)
RIII <- 12L
one_III <- function(n) {
  dat <- sim_data(n, th_full_d)
  Jl <- J_thm4(n, sp_of(s_d, pi_d))
  b <- as.integer(ceiling(log(n)))
  grp <- kp_groups_bal(Jl, b)
  M <- max(grp)
  gsize <- as.numeric(table(grp))
  ws <- sqrt(gsize)
  pj <- prep_J(dat, Jl, th_full_d)
  lam <- lam_w(lam0_cf(n, M, gsize, pj[["smax"]], pj[["lmax"]]), ws)
  th <- gl_fit(pj[["G"]], pj[["h"]], grp, ws, lam, pj[["Lmax"]])
  kk <- max(kkt_gl(pj[["G"]], pj[["h"]], grp, ws, lam, th))
  ch <- coef_c(pj[["rp"]], dat[["y"]], th)
  v <- th - pj[["th_or"]]
  W0 <- Wset(pj[["th_or"]], grp, ws)
  gt <- pj[["gtil"]]
  pred <- mean((as.numeric(pj[["Z"]] %*% c(ch, th)) - dat[["f"]])^2)
  Btv <- mean((pj[["rp"]][["Bt"]] %*% v)^2)
  comp <- sum(v^2) + sum(pj[["tail_blk"]])
  c(n = n, J = Jl, kkt = kk, pred = pred, Btv = Btv,
    bBtv = 64 * lam^2 * W0 / gt + 16 * pj[["bias"]],
    bpred = 3 * (64 * lam^2 * W0 / gt + 20 * pj[["bias"]]) + 3 * pj[["Peps"]],
    l2 = sum(v^2), bl2 = 64 * lam^2 * W0 / gt^2 + 16 * pj[["bias"]] / gt,
    pen = sum(ws * sqrt(as.numeric(rowsum(v^2, grp)))), bpen = 40 * lam * W0 / gt + 10 * pj[["bias"]] / lam,
    comp = comp, bcomp = 2 * (64 * lam^2 * W0 / gt^2 + 16 * pj[["bias"]] / gt) + 2 * sum(pj[["tail_blk"]]),
    schur = gt >= pj[["lmin"]] - 1e-10, bias = pj[["bias"]])
}
resIII <- lapply(nsIII, function(n) t(vapply(seq_len(RIII), function(r) one_III(n), numeric(15))))
tabIII <- t(vapply(seq_along(nsIII), function(i) {
  rr <- resIII[[i]]
  n <- nsIII[i]
  alvo <- n^(-2 * s_d / (2 * s_d + 1))
  c(n = n, J = rr[1, "J"], kkt = max(rr[, "kkt"]), pred = median(rr[, "pred"]),
    alvo = alvo, razao = median(rr[, "pred"]) / alvo,
    razao_log = median(rr[, "pred"]) / (log(n) / n)^(2 * s_d / (2 * s_d + 1)),
    comp = median(rr[, "comp"]), razao_comp = median(rr[, "comp"]) / alvo,
    viola = sum(rr[, "Btv"] > rr[, "bBtv"]) + sum(rr[, "pred"] > rr[, "bpred"]) +
      sum(rr[, "l2"] > rr[, "bl2"]) + sum(rr[, "pen"] > rr[, "bpen"]) + sum(rr[, "comp"] > rr[, "bcomp"]),
    folga = min(rr[, "bpred"] / rr[, "pred"]), schur = mean(rr[, "schur"]))
}, numeric(12)))
print(signif(tabIII, 4))
chk(max(tabIII[, "kkt"]) < 1e-7, "o FISTA resolve o objetivo em blocos com pesos sqrt(|G|) (KKT a 1e-7 de lambda)")
chk(all(tabIII[, "viola"] == 0), "as cotas do Teorema 3(ii)-(iii) com o comparador theta* e a do Corolário 10 valem em todas as réplicas")
chk(all(tabIII[, "schur"] == 1), "gamma_til >= lambda_min(Sigma_hat): sem cone também ao longo de J_n")
for (nm in c("razao", "razao_comp")) {
  rz <- tabIII[, nm]
  cat(sprintf("  %s / n^{-2s'/(2s'+1)}: %s (máx/mín %.2f)\n",
              if (nm == "razao") "||f_hat - f||_n^2" else "sum ||g_hat - g||^2", paste(sprintf("%.3f", rz), collapse = " "),
              max(rz) / min(rz)))
  chk(max(rz) / min(rz) < 2, "a razão ao alvo do Teorema 4 é estável ao variar n por um fator 32")
}
slIII <- unname(coef(lm(log(tabIII[, "pred"]) ~ log(tabIII[, "alvo"])))[2])
cat(sprintf("  inclinação de log(erro) contra log n^{-2s'/(2s'+1)}: %.3f (alvo 1); razão a (log n/n)^{2s'/(2s'+1)}: %s\n",
            slIII, paste(sprintf("%.3f", tabIII[, "razao_log"]), collapse = " ")))
chk(abs(slIII - 1) < 0.25, "o erro realizado segue n^{-2s'/(2s'+1)}, a taxa do Teorema 4")
slIIIlog <- unname(coef(lm(log(tabIII[, "pred"]) ~ log((log(tabIII[, "n"]) / tabIII[, "n"])^(2 * s_d / (2 * s_d + 1)))))[2])
cat(sprintf("  inclinação contra log (log n/n)^{2s'/(2s'+1)} (a forma do Teorema 2, com logaritmo): %.3f\n", slIIIlog))
chk(abs(slIII - 1) < abs(slIIIlog - 1),
    "o erro realizado fica mais perto da taxa sem logaritmo do Teorema 4 que da forma com logaritmo do Teorema 2")

## ===========================================================================
cat("\nPARTE IV. Comprimível, pi = 1, s = 1,2: Corolário 9, Lema 15 e Corolário 11\n")
## ===========================================================================

s_c <- 1.2
pi_c <- 1
Cg_c <- 0.6
sp_c <- sp_of(s_c, pi_c)                            # 0,7
th_full_c <- theta_besov(s_c, pi_c, Cg_c, "uniform", 60L)
cat(sprintf("  s = %.1f, pi = %d: s' = %.2f, tau = %.3f; janela do Corolário 11: c em [%.3f, 1)\n",
            s_c, pi_c, sp_c, tau_of(s_c), s_c / ((2 * s_c + 1) * sp_c)))
nsIV <- c(250L, 500L, 1000L, 2000L, 4000L, 8000L)
RIV <- 10L
one_IV <- function(n) {
  dat <- sim_data(n, th_full_c)
  Jl <- J_cor11(n, s_c, sp_c)
  b <- as.integer(ceiling(log(n)))
  grp <- kp_groups_bal(Jl, b)
  M <- max(grp)
  gsize <- as.numeric(table(grp))
  ws <- sqrt(gsize)
  pj <- prep_J(dat, Jl, th_full_c)
  lam <- lam_w(lam0_cf(n, M, gsize, pj[["smax"]], pj[["lmax"]]), ws)
  th <- gl_fit(pj[["G"]], pj[["h"]], grp, ws, lam, pj[["Lmax"]])
  ch <- coef_c(pj[["rp"]], dat[["y"]], th)
  gt <- pj[["gtil"]]
  Lh <- pj[["lmax"]]                                  # Lambda >= lambda_max(Sigma_hat): o menor válido
  th_or <- pj[["th_or"]]
  Bm <- pj[["rp"]][["B"]]
  Bt <- pj[["rp"]][["Bt"]]
  pred <- mean((as.numeric(pj[["Z"]] %*% c(ch, th)) - dat[["f"]])^2)
  # Teorema 3 com comparadores truncados theta_bar = theta* 1{G em T(eta')}
  en <- as.numeric(rowsum(th_or^2, grp))
  viola_cmp <- 0
  for (et in c(0, en[en > 0] * (1 - 1e-9), Inf)) {
    keepG <- en > et
    thb <- th_or * keepG[grp]
    vb <- th - thb
    bbar <- mean((dat[["f"]] - as.numeric(pj[["Z"]] %*% c(c_true, thb)))^2)
    Wb <- Wset(thb, grp, ws)
    lhs1 <- mean((Bt %*% vb)^2)
    rhs1 <- 64 * lam^2 * Wb / gt + 16 * bbar
    lhs3 <- pred
    rhs3 <- 3 * (64 * lam^2 * Wb / gt + 20 * bbar) + 3 * pj[["Peps"]]
    if (lhs1 > rhs1 || lhs3 > rhs3 || sum(vb^2) > 64 * lam^2 * Wb / gt^2 + 16 * bbar / gt) viola_cmp <- viola_cmp + 1
  }
  # Corolário 9: a forma de risco ideal, com eta = 1,6 max_G (lambda w_G)^2 / (Lambda gamma_til)
  eta <- 1.6 * max((lam * ws)^2) / (Lh * gt)
  R <- Rideal(th_or, grp, eta)
  b9 <- 120 * Lh * R + 120 * pj[["bias"]] + 3 * pj[["Peps"]]
  c9 <- (82 * Lh * R + 64 * pj[["bias"]]) / gt
  # Lema 15 no theta* do desenho, com b = b_n e os pq blocos
  B15 <- bound15(s_c, pi_c, Cg_c, b, eta, nblocks = p * q)
  # o eta determinístico do Corolário 11: lambda_n^G, Lambda e gamma/2 no lugar de lambda, Lambda_hat e gamma_til
  etaD <- 1.6 * max((lam_w(lam0_det(n, M, gsize), ws) * ws)^2) / (Lam * gam / 2)
  alvo <- n^(-2 * s_c / (2 * s_c + 1)) * log(n)^(max(0, 2 / pi_c - 1) / (2 * s_c + 1))
  c(n = n, J = Jl, viola_cmp = viola_cmp, pred = pred, b9 = b9,
    l2 = sum((th - th_or)^2), c9 = c9, R = R, B15 = B15, gt = gt,
    alvo = alvo, comp = sum((th - th_or)^2) + sum(pj[["tail_blk"]]),
    etaD_b = etaD * n / b, RD = Rideal(th_or, grp, etaD), B15D = bound15(s_c, pi_c, Cg_c, b, etaD, nblocks = p * q))
}
resIV <- lapply(nsIV, function(n) t(vapply(seq_len(RIV), function(r) one_IV(n), numeric(15))))
tabIV <- t(vapply(seq_along(nsIV), function(i) {
  rr <- resIV[[i]]
  al <- unname(rr[1, "alvo"])
  c(n = nsIV[i], J = unname(rr[1, "J"]), viola_cmp = sum(rr[, "viola_cmp"]),
    viola9 = sum(rr[, "pred"] > rr[, "b9"]) + sum(rr[, "l2"] > rr[, "c9"]),
    viola15 = sum(rr[, "R"] > rr[, "B15"]) + sum(rr[, "RD"] > rr[, "B15D"]),
    gtil = median(rr[, "gt"]), R = median(rr[, "R"]), etaD_b = median(rr[, "etaD_b"]),
    RD = median(rr[, "RD"]), pred = median(rr[, "pred"]), alvo = al,
    razaoRD = median(rr[, "RD"]) / al, razao = median(rr[, "pred"]) / al,
    razao_comp = median(rr[, "comp"]) / al)
}, numeric(14)))
print(signif(tabIV, 4))
chk(all(tabIV[, "viola_cmp"] == 0), "Teorema 3 com comparador arbitrário: as cotas valem para todo truncamento T(eta') de theta*, em todas as réplicas")
chk(all(tabIV[, "viola9"] == 0), "Corolário 9: ||f_hat - f||_n^2 <= 120 Lambda R_G(theta*; eta) + 120 ||b||^2 + 3 ||P_A eps||^2 e a cota de componentes valem em todas as réplicas")
chk(all(tabIV[, "viola15"] == 0), "Lema 15 no theta* do desenho, com b = b_n e a forma balanceada, em todas as réplicas")
cat(sprintf("  gamma_til mediano: %s (gamma/2 = %.3f): com o gamma_til realizado, o eta do Corolário 9 é grande e R_G satura nesses n\n",
            paste(sprintf("%.3f", tabIV[, "gtil"]), collapse = " "), gam / 2))
kappaIV <- 3.2 * 3 * 16 * (2 * B_X^2 * C_U + Lam * (1 + 0.5 * log(2 * p * q / alpha))) * sig^2 / (Lam * gam)
chk(max(tabIV[, "etaD_b"]) <= kappaIV,
    sprintf("n eta_n / b_n <= kappa = 3,2 rho^2 C_lambda sigma^2/(Lambda gamma) = %.0f no eta determinístico (prova do Corolário 11): eta_n ~ log n / n", kappaIV))
cat(sprintf("  R_G(theta*; eta_n determinístico): %s, contra a energia total %.3f: a constante da calibração satura R_G nesses n; os expoentes estão na Parte I.4\n",
            paste(sprintf("%.3f", tabIV[, "RD"]), collapse = " "), sum(th_full_c[keep_idx(7L)]^2)))
rzP <- tabIV[, "razao"]
cat(sprintf("  ||f_hat - f||_n^2 / alvo: %s (máx/mín %.2f)\n", paste(sprintf("%.3f", rzP), collapse = " "), max(rzP) / min(rzP)))
chk(max(rzP) / min(rzP) < 3, "o erro realizado fica a uma razão limitada do alvo do Corolário 11 ao variar n por um fator 32")

## ===========================================================================
cat("\nPARTE V. A variante branca (o grpreg): Q_G a Gram centrada do pedaço\n")
## ===========================================================================

# Penalidade sum_G w_G ||Q_G^{1/2} theta_G||, Q_G = B^c_G'B^c_G/n (colunas
# centradas). Como X_1 = 1 está em A, M_A <= M_1 e Psi~_G <= Q_G, de modo que
# tr(Q^{-1} Psi~) <= |G| e ||Q^{-1/2} Psi~ Q^{-1/2}|| <= 1: lambda_{0,G} pivotal,
# sigma n^{-1/2}(sqrt|G| + sqrt(2 log(|G|/alpha))). O Teorema 3 vale com
# q_max W no lugar de W, q_max = max_G lambda_max(Q_G).
nV <- 500L
JV <- 5L
bV <- as.integer(ceiling(log(nV)))
grpV <- kp_groups_bal(JV, bV)
MV <- max(grpV)
gsV <- as.numeric(table(grpV))
wV <- sqrt(gsV)
th_full_s <- numeric(length(th_full_d))
ug <- (seq_len(2^14) - 0.5) / 2^14
th_full_s[blk_full == 1L] <- drop(crossprod(psi_block(ug, Jmax), sin(2 * pi * ug))) / length(ug)
truthsV <- list(denso = th_full_d, seno = th_full_s)
one_V <- function(th_full) {
  dat <- sim_data(nV, th_full)
  pj <- prep_J(dat, JV, th_full)
  Bm <- pj[["rp"]][["B"]]
  Bt <- pj[["rp"]][["Bt"]]
  Bc <- sweep(Bm, 2, colMeans(Bm))
  D <- Bt
  Rinv <- vector("list", MV)
  stat <- numeric(MV)
  trr <- numeric(MV)
  opr <- numeric(MV)
  qmax <- 0
  sc <- as.numeric(crossprod(Bt, dat[["eps"]])) / nV
  for (g in seq_len(MV)) {
    ii <- which(grpV == g)
    Qg <- crossprod(Bc[, ii, drop = FALSE]) / nV
    qmax <- max(qmax, max(eigen(Qg, symmetric = TRUE, only.values = TRUE)[["values"]]))
    Rg <- chol(Qg)
    Rinv[[g]] <- backsolve(Rg, diag(length(ii)))
    D[, ii] <- Bt[, ii, drop = FALSE] %*% Rinv[[g]]
    Pg <- crossprod(D[, ii, drop = FALSE]) / nV          # Q^{-1/2} Psi~ Q^{-1/2} (forma de Cholesky)
    evP <- eigen(Pg, symmetric = TRUE, only.values = TRUE)[["values"]]
    trr[g] <- sum(evP) / length(ii)
    opr[g] <- max(evP)
    stat[g] <- sqrt(sum((t(Rinv[[g]]) %*% sc[ii])^2))  # ||R^{-T} z_G|| = ||Q^{-1/2} z_G|| em norma
  }
  l0 <- sig / sqrt(nV) * (sqrt(gsV) + sqrt(2 * log(MV / alpha)))
  lam <- lam_w(l0, wV)
  GD <- crossprod(D) / nV
  hD <- as.numeric(crossprod(D, dat[["y"]] - pj[["rp"]][["proj"]](dat[["y"]]))) / nV
  phi <- gl_fit(GD, hD, grpV, wV, lam, max(eigen(GD, symmetric = TRUE, only.values = TRUE)[["values"]]))
  th <- phi
  for (g in seq_len(MV)) {
    ii <- which(grpV == g)
    th[ii] <- Rinv[[g]] %*% phi[ii]
  }
  ch <- coef_c(pj[["rp"]], dat[["y"]], th)
  v <- th - pj[["th_or"]]
  W0 <- Wset(pj[["th_or"]], grpV, wV)
  gt <- pj[["gtil"]]
  pred <- mean((as.numeric(pj[["Z"]] %*% c(ch, th)) - dat[["f"]])^2)
  Btv <- mean((Bt %*% v)^2)
  c(cob = max(stat / wV) <= lam / 2, tr = max(trr), op = max(opr),
    viola = (Btv > 64 * lam^2 * qmax * W0 / gt + 16 * pj[["bias"]]) +
      (sum(v^2) > 64 * lam^2 * qmax * W0 / gt^2 + 16 * pj[["bias"]] / gt) +
      (pred > 3 * (64 * lam^2 * qmax * W0 / gt + 20 * pj[["bias"]]) + 3 * pj[["Peps"]]),
    qmax = qmax, lmax = pj[["lmax"]], folga = (64 * lam^2 * qmax * W0 / gt + 16 * pj[["bias"]]) / Btv)
}
tabV <- t(vapply(names(truthsV), function(nm) {
  rr <- t(vapply(seq_len(40L), function(r) one_V(truthsV[[nm]]), numeric(7)))
  c(cob = mean(rr[, "cob"]), tr_max = max(rr[, "tr"]), op_max = max(rr[, "op"]), viola = sum(rr[, "viola"]),
    qmax_lmax = max(rr[, "qmax"] / rr[, "lmax"]), folga = min(rr[, "folga"]))
}, numeric(6)))
print(round(tabV, 4))
chk(all(tabV[, "tr_max"] <= 1 + 1e-10) && all(tabV[, "op_max"] <= 1 + 1e-10),
    "Psi~_G <= Q_G: tr(Q^{-1} Psi~)/|G| <= 1 e ||Q^{-1/2} Psi~ Q^{-1/2}|| <= 1 em todos os pedaços (X_1 = 1 em A)")
chk(all(tabV[, "cob"] >= 1 - alpha), "a calibração pivotal da variante branca tem a probabilidade nominal")
chk(all(tabV[, "viola"] == 0), "Teorema 3 na variante branca, com q_max W(G_0) no lugar de W(G_0): as três cotas valem em todas as réplicas")
chk(all(tabV[, "qmax_lmax"] <= 1 + 1e-10), "q_max = max_G lambda_max(Q_G) <= lambda_max(Sigma_hat)")

## ===========================================================================
cat("\nPARTE VI. O estimador limiarizado (06: Lema 16, Corolários 12 e 13)\n")
## ===========================================================================

## ---- VI.1 Lema 16 em configurações aleatórias -------------------------------
set.seed(77)
viol16 <- c(igual = 0, sem_hip = 0, com_hip = 0, pred_sem = 0, pred_com = 0)
n16 <- 0
for (rep in seq_len(2000L)) {
  K <- sample(2:12, 1)
  szb <- sample(3:15, K, TRUE)
  blk <- rep(seq_len(K), szb)
  dd <- length(blk)
  act <- runif(K) < 0.5
  ths <- rnorm(dd) * runif(1, 0.05, 2) * act[blk]           # theta* por bloco, alguns nulos
  that <- ths + rnorm(dd) * runif(1, 0.01, 1)
  cs <- rnorm(2)
  chh <- cs + rnorm(2) * 0.1
  Zr <- matrix(rnorm((dd + 2) * 3 * (dd + 2)), 3 * (dd + 2), dd + 2) %*% diag(runif(dd + 2, 0.3, 2))
  Sh <- crossprod(Zr) / nrow(Zr)
  ev <- eigen(Sh, symmetric = TRUE, only.values = TRUE)[["values"]]
  nu <- sqrt(as.numeric(rowsum(ths^2, blk)))
  nuh <- sqrt(as.numeric(rowsum(that^2, blk)))
  Del <- sqrt(as.numeric(rowsum((that - ths)^2, blk)))
  S <- nu > 0
  for (t in c(0, quantile(nuh, c(0.1, 0.3, 0.5, 0.7, 0.9)), max(nuh) + 1)) {
    n16 <- n16 + 1
    Sh_t <- nuh > t
    tht <- that * Sh_t[blk]
    err_t <- sum((tht - ths)^2)
    lhs_eq <- sum(Del[Sh_t]^2) + sum(nu[S & !Sh_t]^2)
    if (abs(err_t - lhs_eq) > 1e-10 * max(1, err_t)) viol16["igual"] <- viol16["igual"] + 1
    if (err_t > 2 * sum(Del^2) + 2 * t^2 * sum(S & !Sh_t) + 1e-12) viol16["sem_hip"] <- viol16["sem_hip"] + 1
    if (all(Sh_t[S]) && err_t > sum(Del^2) + 1e-12) viol16["com_hip"] <- viol16["com_hip"] + 1
    dlt <- c(chh - cs, that - ths)
    dlt_t <- c(chh - cs, tht - ths)
    pz <- drop(t(dlt) %*% Sh %*% dlt)
    pz_t <- drop(t(dlt_t) %*% Sh %*% dlt_t)
    if (pz_t > max(ev) / min(ev) * 2 * pz + 2 * max(ev) * t^2 * sum(S & !Sh_t) + 1e-10) viol16["pred_sem"] <- viol16["pred_sem"] + 1
    if (all(Sh_t[S]) && pz_t > max(ev) / min(ev) * pz + 1e-10) viol16["pred_com"] <- viol16["pred_com"] + 1
  }
}
cat(sprintf("  %d pares (configuração, limiar): violações %s\n", n16, paste(names(viol16), viol16, sep = " = ", collapse = ", ")))
chk(viol16[["igual"]] == 0, "Lema 16(i): erro do limiarizado = soma dos Delta^2 mantidos + soma dos nu^2 dos ativos zerados")
chk(viol16[["sem_hip"]] == 0, "Lema 16(i): erro do limiarizado <= 2 sum Delta^2 + 2 t^2 |S minus S_hat(t)|, sem hipótese")
chk(viol16[["com_hip"]] == 0, "Lema 16(ii): no evento S contido em S_hat(t), o erro do limiarizado <= o do ajuste")
chk(viol16[["pred_sem"]] == 0 && viol16[["pred_com"]] == 0,
    "Lema 16(iii): a passagem à norma de predição com lambda_max/lambda_min de Sigma_hat, com e sem a hipótese")
# O fator 2 é atingido: um bloco ativo com theta* = (t + D) e, theta_hat = t e.
tt <- 0.3
DD <- 0.3
err_ex <- (tt + DD)^2
chk(abs(err_ex - (2 * DD^2 + 2 * tt^2)) < 1e-15,
    "a cota 2 Delta^2 + 2 t^2 é atingida (bloco ativo com norma estimada t, alinhado, t = Delta)")

## ---- VI.2 Corolários 12 e 13 nos dois ajustes ------------------------------
# Verdade com estrutura: (1,1) Besov pi = 2, s = 1/2; (2,2) Besov pi = 1,
# s = 1,2; (1,2) e (2,1) nulos. Os dois ajustes: o LASSO perfilado (glmnet)
# no lambda do Corolário 2 de E1.5 com sigma_max e o block LASSO balanceado
# com sqrt(|G|) no lambda fechado do Lema 14; cada um no lambda da teoria e
# em um décimo dele (onde os blocos nulos ficam com energia), e o limiar t em grade.
th_full_st <- numeric(length(th_full_d))
th_full_st[blk_full == 1L] <- theta_besov_blk(0.5, 2, 0.6, Jmax, "random", 91L)
th_full_st[blk_full == 4L] <- theta_besov_blk(1.2, 1, 0.6, Jmax, "uniform", 92L)
nu_true <- sqrt(as.numeric(rowsum(th_full_st^2, blk_full)))
S_true <- nu_true > 0
cat(sprintf("  normas dos blocos verdadeiros: %s\n", paste(sprintf("%.3f", nu_true), collapse = " ")))
fit_lasso <- function(pj, y, lam) {
  yt <- as.numeric(y - pj[["rp"]][["proj"]](y))
  fit <- glmnet(pj[["rp"]][["Bt"]], yt, family = "gaussian", lambda = c(4, 2, 1) * lam,
                standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-14))
  as.numeric(fit[["beta"]][, 3])
}
tgrid <- seq(0, 1.2, by = 0.02)
one_VI <- function(n) {
  dat <- sim_data(n, th_full_st)
  Jl <- J_thm4(n, 0.5)
  b <- as.integer(ceiling(log(n)))
  grp <- kp_groups_bal(Jl, b)
  M <- max(grp)
  gsize <- as.numeric(table(grp))
  ws <- sqrt(gsize)
  pj <- prep_J(dat, Jl, th_full_st)
  d <- length(pj[["th_or"]])
  lamL <- 2 * sig * pj[["smax"]] * sqrt(2 * log(2 * d / alpha) / n)
  lamB <- lam_w(lam0_cf(n, M, gsize, pj[["smax"]], pj[["lmax"]]), ws)
  fits <- list(lasso = fit_lasso(pj, dat[["y"]], lamL), lasso10 = fit_lasso(pj, dat[["y"]], lamL / 10),
               blocos = gl_fit(pj[["G"]], pj[["h"]], grp, ws, lamB, pj[["Lmax"]]),
               blocos10 = gl_fit(pj[["G"]], pj[["h"]], grp, ws, lamB / 10, pj[["Lmax"]]))
  blk <- pj[["blk"]]
  S_J <- sqrt(as.numeric(rowsum(pj[["th_or"]]^2, blk))) > 0
  out <- lapply(names(fits), function(nm) {
    th <- fits[[nm]]
    ch <- coef_c(pj[["rp"]], dat[["y"]], th)
    nuh <- sqrt(as.numeric(rowsum(th^2, blk)))
    Del2 <- as.numeric(rowsum((th - pj[["th_or"]])^2, blk)) + pj[["tail_blk"]]   # ||g_hat - g||^2 por bloco
    dlt <- c(ch - c_true, th - pj[["th_or"]])
    pz <- mean((pj[["Z"]] %*% dlt)^2)                                            # ||f_hat - f_J||_n^2
    # Corolário 12 (blocos) e Corolário 8 (LASSO): Delta^2 <= 2||v||^2 + 2 max tail, com ||v||^2 pela cota realizada
    v2b <- if (grepl("blocos", nm)) {
      lam <- if (nm == "blocos") lamB else lamB / 10
      64 * lam^2 * Wset(pj[["th_or"]], grp, ws) / pj[["gtil"]]^2 + 16 * pj[["bias"]] / pj[["gtil"]]
    } else {
      lam <- if (nm == "lasso") lamL else lamL / 10
      64 * lam^2 * sum(pj[["th_or"]] != 0) / pj[["gtil"]]^2 + 16 * pj[["bias"]] / pj[["gtil"]]
    }
    Dbar2 <- 2 * v2b + 2 * max(pj[["tail_blk"]])
    teoria <- nm %in% c("lasso", "blocos")           # só no lambda da teoria o evento T vale
    viola12 <- if (teoria) max(Del2) > Dbar2 else NA
    tri <- if (teoria) all(!(nuh > sqrt(Dbar2))[!S_true]) else NA               # triagem em t = Delta_bar
    rr <- t(vapply(tgrid, function(t) {
      kp <- nuh > t
      tht <- th * kp[blk]
      err_t <- sum(Del2[kp]) + sum(nu_true[!kp]^2)
      dlt_t <- c(ch - c_true, tht - pj[["th_or"]])
      pz_t <- mean((pj[["Z"]] %*% dlt_t)^2)
      nfn <- sum(S_true & !kp)
      nfnJ <- sum(S_J & !kp)
      c(sem = err_t <= 2 * sum(Del2) + 2 * t^2 * nfn + 1e-12,
        com = if (nfn == 0) err_t <= sum(Del2) + 1e-12 else NA,
        psem = pz_t <= pj[["lmax"]] / pj[["lmin"]] * 2 * pz + 2 * pj[["lmax"]] * t^2 * nfnJ + 1e-12,
        pcom = if (nfnJ == 0) pz_t <= pj[["lmax"]] / pj[["lmin"]] * pz + 1e-12 else NA,
        exato = all(kp == S_true), ganho = err_t / sum(Del2),
        ganho_pred = (mean((pj[["Z"]] %*% c(ch, tht) - dat[["f"]])^2)) / mean((pj[["Z"]] %*% c(ch, th) - dat[["f"]])^2))
    }, numeric(7)))
    list(viola12 = viola12, tri = tri, sem = all(rr[, "sem"] == 1), com = all(rr[, "com"] == 1, na.rm = TRUE),
         psem = all(rr[, "psem"] == 1), pcom = all(rr[, "pcom"] == 1, na.rm = TRUE),
         frac_exato = mean(rr[, "exato"]),
         ganho = if (any(rr[, "exato"] == 1)) median(rr[rr[, "exato"] == 1, "ganho"]) else NA,
         ganho_pred = if (any(rr[, "exato"] == 1)) median(rr[rr[, "exato"] == 1, "ganho_pred"]) else NA,
         nulos = sum(Del2[!S_true]) / sum(Del2))
  })
  names(out) <- names(fits)
  out
}
nsVI <- c(500L, 2000L)
RVI <- 15L
tabVI <- do.call(rbind, lapply(nsVI, function(n) {
  res <- lapply(seq_len(RVI), function(r) one_VI(n))
  do.call(rbind, lapply(c("lasso", "lasso10", "blocos", "blocos10"), function(nm) {
    g <- function(k) vapply(res, function(z) as.numeric(z[[nm]][[k]]), 1)
    data.frame(n = n, ajuste = nm, viola12 = sum(g("viola12")), triagem = mean(g("tri")),
               sem = all(g("sem") == 1), com = all(g("com") == 1), psem = all(g("psem") == 1),
               pcom = all(g("pcom") == 1), frac_exato = median(g("frac_exato")),
               ganho = median(g("ganho"), na.rm = TRUE), ganho_pred = median(g("ganho_pred"), na.rm = TRUE),
               nulos = median(g("nulos")))
  }))
}))
print(tabVI, digits = 3, row.names = FALSE)
teo <- tabVI[["ajuste"]] %in% c("lasso", "blocos")
chk(all(tabVI[teo, "viola12"] == 0), "Corolários 8 e 12: o erro máximo por bloco fica sob Delta_bar (forma realizada) nos dois ajustes no lambda da teoria, em todas as réplicas")
chk(all(tabVI[teo, "triagem"] == 1), "triagem: com t = Delta_bar nenhum bloco nulo sobrevive (Lema 13(i)), nos dois ajustes")
chk(all(tabVI[["sem"]]) && all(tabVI[["com"]]),
    "Corolário 13, componentes: <= 2 sum ||g_hat - g||^2 + 2 t^2 |S minus S_hat(t)| em todo t, e <= sum ||g_hat - g||^2 quando S está em S_hat(t)")
chk(all(tabVI[["psem"]]) && all(tabVI[["pcom"]]),
    "Corolário 13, predição: a passagem com lambda_max/lambda_min de Sigma_hat vale nos dois ajustes em todo t")
g10 <- tabVI[tabVI[["ajuste"]] %in% c("lasso10", "blocos10"), ]
chk(all(g10[["ganho"]] <= 1) && all(g10[["frac_exato"]] > 0),
    "com lambda/10, o platô de acerto existe e nele o erro das componentes do limiarizado não passa o do ajuste (Lema 16(ii)); a coluna nulos é a parte do erro que ele tira")
cat("  (ganho_pred: razão do erro de predição limiarizado / ajuste no platô; pode passar de 1, dentro da cota de Lema 16(iii))\n")

cat(sprintf("\n  tempo: %.0f s\n", proc.time()[["elapsed"]] - t_start))
if (ok) cat("OK\n") else stop("E1.12: conferência numérica FALHOU (ver linhas acima)")
