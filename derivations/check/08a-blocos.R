# E1.11. Conferência numérica da sondagem derivations/08a-sondagem-blocos.md:
# o que a teoria de E1.4 a E1.7c vira quando a penalidade ||theta||_1 dá lugar
# ao block LASSO de Klopp & Pensky (2015) com os níveis c_l livres (D3), com
# pedaços de até b translações consecutivas dentro de cada bloco (l,m), sem
# atravessar nível, como em wafc_kp_groups() (wafc/R/competitors.R).
#
# É conferência de SONDAGEM: confere os números que o documento afirma, não
# uma prova. Imprime OK ou falha com stop(). Notação: docs/notacao.md.
#
# O estimador conferido é
#     (c_hat, theta_hat) = argmin ||Y - A c - B theta||_n^2 + 2 lambda sum_g w_g N_g(theta_g),
# com N_g a norma euclidiana do pedaço ("euclid", a de K&P e de Lounici et al.)
# ou a norma do ajuste do pedaço, ||B~_g theta_g||_n ("white", a que o grpreg
# penaliza depois de ortonormalizar cada grupo), e w_g = 1 (K&P) ou
# sqrt(|G_g|) (o padrão do grpreg, que é o que o método klopp do piloto usa).
# O ajuste é por FISTA com reinício no problema perfilado (Lema 4 de E1.5);
# a Parte B confere as condições KKT e, no caso de grupos unitários, o glmnet.
#
# Partes:
#   A. calibração: cobertura do evento max_g ||(B~' eps / n)_g|| / w_g <= lambda/2
#      com o lambda de blocos (Hsu, Kakade & Zhang 2012), para os dois pesos e
#      as duas normas, e a escala do lambda de blocos contra o do LASSO;
#   B. o Teorema 1 de E1.5 em blocos: a cota rápida com W_S = sum_{g em S} w_g^2
#      no lugar de s_0, sem cone (gamma_til >= lambda_min(Sigma_hat)), a razão
#      ||f_hat - f||_n^2 / (lambda^2 W_S) contra n, e blocos contra LASSO no
#      lambda da teoria e no melhor lambda de uma grade;
#   C. Besov implica a cota do risco ideal por blocos (a forma do Lema 4 de
#      K&P, arXiv v2), sem hipótese nova, e o expoente de b que ela prevê;
#   D. a aritmética dos expoentes das taxas (LASSO contra blocos).
#
# Nos scripts de conferência, elemento de lista se acessa com [[ ]] e nome
# completo (instrucoes.md, §5).
#
# Dependências: WaveBased (wbasis, wtable), glmnet.
# Tempo nesta máquina: cerca de 2 min.

suppressPackageStartupMessages({
  library(WaveBased)
  library(glmnet)
})
set.seed(20260930)

p  <- 2L                        # X_1 = 1, X_2 = 1/2 + Unif(-1,1) (cenário A de E1.4)
q  <- 2L
fs <- 8L                        # filter.size (Daublets, 4 momentos nulos)
sig <- 0.5
alpha <- 0.05
gam <- 0.25                     # kappa_1 c_U no cenário A (U uniforme: c_U = C_U = 1)
tb <- wtable(filter.size = fs)

ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "ok" else "FALHA", msg))
  if (!isTRUE(cond)) ok <<- FALSE
  invisible(cond)
}

## ---- desenho (mesma construção de check/04-oraculo.R) ---------------------

psi_block <- function(u, Jl) {
  wbasis(u, j0 = 0, J = Jl, filter.size = fs, wavelet.table = tb)[, -1, drop = FALSE]
}
Psi <- function(U, Jl) {
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
design <- function(X, U, Jl) {
  kron_rows(X, Psi(U, Jl))[, perm_D12(ncol(X), ncol(U), 2L^Jl - 1L), drop = FALSE]
}
draw_X <- function(U) cbind(1, 0.5 + runif(nrow(U), -1, 1))
resid_pen <- function(Z) {
  A <- Z[, seq_len(p), drop = FALSE]
  Bm <- Z[, -seq_len(p), drop = FALSE]
  qa <- qr(A)
  list(A = A, B = Bm, Bt = Bm - qr.fitted(qa, Bm), proj = function(v) qr.fitted(qa, v))
}

## ---- os pedaços de K&P ------------------------------------------------------

# Grupo de cada coluna penalizada, na ordem de D12 (blocos (l,m) lexicográficos,
# dentro do bloco j e k crescentes). "nivel": pedaços de até b translações
# consecutivas sem atravessar nível, como wafc_kp_groups(); "consecutivo":
# pedaços de b índices consecutivos da base, atravessando níveis, como a
# norma (3.1) de K&P.
kp_groups <- function(Jl, b, nblocks = p * q, kind = c("nivel", "consecutivo")) {
  kind <- match.arg(kind)
  NJ <- 2L^Jl - 1L
  out <- integer(0)
  cur <- 0L
  for (bb in seq_len(nblocks)) {
    if (kind == "consecutivo") {
      k <- ceiling(NJ / b)
      out <- c(out, cur + rep(seq_len(k), each = b)[seq_len(NJ)])
      cur <- cur + k
    } else {
      for (jj in 0:(Jl - 1L)) {
        nj <- 2L^jj
        k <- ceiling(nj / b)
        out <- c(out, cur + rep(seq_len(k), each = b)[seq_len(nj)])
        cur <- cur + k
      }
    }
  }
  as.integer(out)
}

## ---- solver: group LASSO no problema perfilado -----------------------------

# minimiza (1/n)||yt - D x||^2 + 2 lam sum_g w_g ||x_g||_2, com G = D'D/n e
# h = D'yt/n; FISTA com reinício por gradiente (O'Donoghue e Candès).
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
gl_fit <- function(G, h, grp, w, lam, x0 = NULL, tol = 1e-9, maxit = 50000L) {
  d <- length(h)
  L <- 2 * max(eigen(G, symmetric = TRUE, only.values = TRUE)[["values"]])
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

# Um ajuste completo: devolve (c_hat, theta_hat) no parâmetro original. Para
# a norma "white" o problema é resolvido nas colunas ortonormalizadas por
# pedaço, C_g = B~_g R_g^{-1} com R_g'R_g = B~_g'B~_g/n, e theta_g = R_g^{-1} phi_g.
prep_fit <- function(rp, y, grp, norm = c("euclid", "white")) {
  norm <- match.arg(norm)
  Bt <- rp[["Bt"]]
  n <- nrow(Bt)
  yt <- as.numeric(y - rp[["proj"]](y))
  if (norm == "euclid") {
    D <- Bt
    back <- function(x) x
  } else {
    Rinv <- vector("list", max(grp))
    D <- Bt
    for (g in seq_len(max(grp))) {
      ii <- which(grp == g)
      Rg <- chol(crossprod(Bt[, ii, drop = FALSE]) / n)
      Rinv[[g]] <- backsolve(Rg, diag(length(ii)))
      D[, ii] <- Bt[, ii, drop = FALSE] %*% Rinv[[g]]
    }
    back <- function(x) {
      out <- x
      for (g in seq_len(max(grp))) {
        ii <- which(grp == g)
        out[ii] <- Rinv[[g]] %*% x[ii]
      }
      out
    }
  }
  list(G = crossprod(D) / n, h = as.numeric(crossprod(D, yt)) / n, back = back)
}
coef_full <- function(rp, y, th) c(qr.coef(qr(rp[["A"]]), y - rp[["B"]] %*% th), th)

## ---- lambdas da teoria --------------------------------------------------------

# Lema 5 de E1.5 (LASSO): lambda = 2 lambda_0, lambda_0 = sigma smax sqrt(2 log(2d/alpha)/n).
lam_lasso <- function(n, d, smax) 2 * sig * smax * sqrt(2 * log(2 * d / alpha) / n)
# Calibração em blocos (HKZ 2012 com t = log(M/alpha) e união sobre os M pedaços):
#   ||(B~'eps/n)_g|| <= (sigma/sqrt n)(sqrt(tr Psi_g) + sqrt(2 ||Psi_g|| log(M/alpha))),
# Psi_g = B~_g'B~_g/n. Forma fechada usada no documento, com tr Psi_g <= |G_g| smax^2
# e ||Psi_g|| <= Lhat = max_g ||Psi_g||.
lam0_g <- function(n, M, gsize, smax, Lhat) {
  sig / sqrt(n) * (smax * sqrt(gsize) + sqrt(2 * Lhat * log(M / alpha)))
}
# Na norma "white" o escore do pedaço é ||P_g eps||/sqrt(n): pivotal, qui-quadrado.
lam0_white <- function(n, M, gsize) sig / sqrt(n) * (sqrt(gsize) + sqrt(2 * log(M / alpha)))

group_stats <- function(Bt, grp) {
  n <- nrow(Bt)
  M <- max(grp)
  tr <- numeric(M)
  op <- numeric(M)
  Q <- vector("list", M)
  for (g in seq_len(M)) {
    ii <- which(grp == g)
    Pg <- crossprod(Bt[, ii, drop = FALSE]) / n
    ev <- eigen(Pg, symmetric = TRUE, only.values = TRUE)[["values"]]
    tr[g] <- sum(ev)
    op[g] <- max(ev)
    Q[[g]] <- qr.Q(qr(Bt[, ii, drop = FALSE]))
  }
  list(tr = tr, op = op, Q = Q)
}

## ===========================================================================
cat("PARTE A. Calibração em blocos e a escala do lambda\n")
## ===========================================================================

cenA <- list(list(n = 100L, J = 4L), list(n = 200L, J = 4L),
             list(n = 500L, J = 4L), list(n = 500L, J = 5L))
RA <- 300L
tabA <- t(vapply(cenA, function(cc) {
  n <- cc[["n"]]
  Jl <- cc[["J"]]
  b <- as.integer(ceiling(log(n)))
  grp <- kp_groups(Jl, b)
  M <- max(grp)
  gsize <- as.numeric(table(grp))
  d <- length(grp)
  res <- t(vapply(seq_len(RA), function(r) {
    U <- matrix(runif(n * q), n, q)
    X <- draw_X(U)
    Z <- design(X, U, Jl)
    rp <- resid_pen(Z)
    Bt <- rp[["Bt"]]
    eps <- rnorm(n, sd = sig)
    sc <- as.numeric(crossprod(Bt, eps)) / n
    smax <- sqrt(max(colSums(rp[["B"]]^2) / n))
    gs <- group_stats(Bt, grp)
    Lhat <- max(gs[["op"]])
    nrm <- sqrt(as.numeric(rowsum(sc^2, grp)))
    white <- vapply(seq_len(M), function(g) sqrt(sum(crossprod(gs[["Q"]][[g]], eps)^2)), 1) / sqrt(n)
    lL <- lam_lasso(n, d, smax)
    l0 <- lam0_g(n, M, gsize, smax, Lhat)
    lE1 <- 2 * max(l0)                      # euclid, w = 1 (K&P), forma fechada
    lEs <- 2 * max(l0 / sqrt(gsize))        # euclid, w = sqrt(|G|)
    w0 <- lam0_white(n, M, gsize)
    lW1 <- 2 * max(w0)                      # white, w = 1
    lWs <- 2 * max(w0 / sqrt(gsize))        # white, w = sqrt(|G|): o padrão do grpreg
    c(cL = max(abs(sc)) <= lL / 2, cE1 = max(nrm) <= lE1 / 2,
      cEs = max(nrm / sqrt(gsize)) <= lEs / 2, cW1 = max(white) <= lW1 / 2,
      cWs = max(white / sqrt(gsize)) <= lWs / 2,
      lL = lL, lE1 = lE1, lEs = lEs, lW1 = lW1, lWs = lWs, smax = smax, Lhat = Lhat,
      lmin = min(eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]]),
      fL = (lL / 2) / max(abs(sc)), fE1 = (lE1 / 2) / max(nrm))
  }, numeric(15)))
  bmax <- max(gsize)
  c(n = n, J = Jl, b = b, d = d, M = M, cobL = mean(res[, "cL"]), cobE1 = mean(res[, "cE1"]),
    cobEs = mean(res[, "cEs"]), cobW1 = mean(res[, "cW1"]), cobWs = mean(res[, "cWs"]),
    lamL = median(res[, "lL"]), lamE1 = median(res[, "lE1"]), lamEs = median(res[, "lEs"]),
    # custo por coordenada de um pedaço cheio (W_S / coordenadas = 1/b) e de
    # um coeficiente isolado (um pedaço por coeficiente), contra o LASSO
    cheio = median(res[, "lE1"]^2 / bmax / res[, "lL"]^2),
    isolado = median(res[, "lE1"]^2 / res[, "lL"]^2),
    grpreg = median(res[, "lEs"]^2 / res[, "lL"]^2),
    whiteW1 = median(res[, "lW1"]^2 / bmax / res[, "lL"]^2),
    smax = median(res[, "smax"]), Lhat = median(res[, "Lhat"]), lmin = median(res[, "lmin"]),
    folgaL = median(res[, "fL"]), folgaE1 = median(res[, "fE1"]))
}, numeric(22)))
print(round(tabA[, c("n", "J", "b", "d", "M", "cobL", "cobE1", "cobEs", "cobW1", "cobWs", "folgaL", "folgaE1")], 3))
print(round(tabA[, c("n", "J", "lamL", "lamE1", "lamEs", "cheio", "isolado", "grpreg", "whiteW1", "smax", "Lhat")], 3))
chk(all(tabA[, c("cobL", "cobE1", "cobEs", "cobW1", "cobWs")] >= 1 - alpha),
    "o evento de calibração em blocos tem a probabilidade nominal nos quatro esquemas (e o do LASSO)")
chk(all(tabA[, "cheio"] < 0.5),
    "num pedaço cheio, o custo por coordenada do lambda de blocos (pesos 1) é menos da metade do LASSO")
chk(all(tabA[, "isolado"] > 1),
    "num coeficiente isolado (um pedaço por coeficiente), o custo é MAIOR que o do LASSO")
chk(all(tabA[, "grpreg"] > 0.8),
    "com pesos sqrt(|G|) o lambda da teoria é ditado pelos pedaços unitários e o custo por coordenada volta à ordem do LASSO")

# A escala em n, só por fórmula (smax = Lhat = 1): a razão do custo por
# coordenada de um pedaço cheio decai como 1/log n.
formA <- t(vapply(c(1e2, 1e3, 1e4, 1e5, 1e6), function(n) {
  vapply(c(1 / 3, 0.6), function(cc) {
    Jl <- max(2L, as.integer(ceiling(cc * log2(n))))
    b <- as.integer(ceiling(log(n)))
    grp <- kp_groups(Jl, b)
    M <- max(grp)
    d <- length(grp)
    (sqrt(max(table(grp))) + sqrt(2 * log(M / alpha)))^2 / max(table(grp)) /
      (2 * log(2 * d / alpha))
  }, 1)
}, numeric(2)))
dimnames(formA) <- list(paste0("n=1e", 2:6), c("c=1/3", "c=0.6"))
cat("  razão (lambda_blocos^2 / b) / lambda_LASSO^2, fórmula com smax = Lhat = 1, J = ceiling(c log2 n):\n")
print(round(formA, 3))
chk(all(diff(formA[, 1]) < 0) && all(diff(formA[, 2]) < 0),
    "a vantagem por coordenada cresce com n (a razão decai), como 1/log n")

## ===========================================================================
cat("\nPARTE B. O Teorema 1 de E1.5 com a norma de blocos\n")
## ===========================================================================

JB <- 4L
bB <- 4L
grpB <- kp_groups(JB, bB)                   # por bloco: {1}, {2,3}, {4..7}, {8..11}, {12..15}
MB <- max(grpB)
dB <- length(grpB)
gsB <- as.numeric(table(grpB))
NJB <- 2L^JB - 1L
c_true <- c(1, -0.5)
pos <- function(bb, k) (bb - 1L) * NJB + k  # coluna penalizada k do bloco bb
set.seed(11)
th_cheio <- numeric(dB)                     # T1: três pedaços cheios
th_cheio[pos(1L, 4:11)] <- sample(c(-1, 1), 8, TRUE)
th_cheio[pos(4L, 12:15)] <- sample(c(-1, 1), 4, TRUE)
# T2: os mesmos 12 coeficientes, um por pedaço. As posições de (1,m) e (2,m)
# são distintas: neste desenho X_1 psi_jk(U_m) e X_2 psi_jk(U_m) têm correlação
# 0,65 (é o kappa_1 = 0,25), e repetir a wavelet nos dois blocos misturaria o
# efeito dos pedaços com o do mau condicionamento.
th_esp <- numeric(dB)
posE <- list(c(1L, 4L, 8L), c(2L, 5L, 12L), c(3L, 6L, 9L), c(1L, 7L, 14L))
for (bb in 1:4) th_esp[pos(bb, posE[[bb]])] <- sample(c(-1, 1), 3, TRUE)
gsin <- function(u) sin(2 * pi * u)         # T3: g_11 = sin(2 pi u), com viés
ug <- (seq_len(2^14) - 0.5) / 2^14
th_sin_blk <- drop(crossprod(psi_block(ug, JB), gsin(ug))) / length(ug)
th_sin <- numeric(dB)
th_sin[pos(1L, 1:NJB)] <- th_sin_blk
g_res <- function(u) gsin(u) - drop(psi_block(u, JB) %*% th_sin_blk)
truths <- list(cheio = list(th = th_cheio, extra = NULL),
               espalhado = list(th = th_esp, extra = NULL),
               seno = list(th = th_sin, extra = g_res))
Wset <- function(th, grp, w) {
  act <- which(sqrt(as.numeric(rowsum(th^2, grp))) > 0)
  sum(w[act]^2)
}
cat(sprintf("  J = %d, b = %d, d = %d, M = %d pedaços; W_S (pesos 1): cheio %d, espalhado %d, seno %d; coeficientes: 12, 12, %d\n",
            JB, bB, dB, MB, Wset(th_cheio, grpB, rep(1, MB)), Wset(th_esp, grpB, rep(1, MB)),
            Wset(th_sin, grpB, rep(1, MB)), sum(th_sin != 0)))

one_B <- function(n, tr) {
  th <- tr[["th"]]
  U <- matrix(runif(n * q), n, q)
  X <- draw_X(U)
  Z <- design(X, U, JB)
  fJ <- as.numeric(Z %*% c(c_true, th))
  f <- if (is.null(tr[["extra"]])) fJ else fJ + X[, 1] * tr[["extra"]](U[, 1])
  y <- f + rnorm(n, sd = sig)
  rp <- resid_pen(Z)
  Bt <- rp[["Bt"]]
  smax <- sqrt(max(colSums(rp[["B"]]^2) / n))
  gs <- group_stats(Bt, grpB)
  lamB <- 2 * max(lam0_g(n, MB, gsB, smax, max(gs[["op"]])))
  lamL <- lam_lasso(n, dB, smax)
  pe <- prep_fit(rp, y, grpB, "euclid")
  thB <- gl_fit(pe[["G"]], pe[["h"]], grpB, rep(1, MB), lamB)
  kk <- kkt_gl(pe[["G"]], pe[["h"]], grpB, rep(1, MB), lamB, thB)
  thL <- gl_fit(pe[["G"]], pe[["h"]], seq_len(dB), rep(1, dB), lamL)
  vB <- thB - th
  gt <- min(eigen(crossprod(Bt) / n, symmetric = TRUE, only.values = TRUE)[["values"]])
  lsig <- min(eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]])
  WS <- Wset(th, grpB, rep(1, MB))
  bn <- mean((f - fJ)^2)
  predB <- mean((as.numeric(Z %*% coef_full(rp, y, thB)) - f)^2)
  predL <- mean((as.numeric(Z %*% coef_full(rp, y, thL)) - f)^2)
  c(lam = lamB, lamL = lamL, Btv = mean((Bt %*% vB)^2), bound = 64 * lamB^2 * WS / gt + 16 * bn,
    l2 = sum(vB^2), l2b = 64 * lamB^2 * WS / gt^2 + 16 * bn / gt,
    penv = sum(sqrt(as.numeric(rowsum(vB^2, grpB)))), penb = 40 * lamB * WS / gt + 10 * bn / lamB,
    gt = gt, lsig = lsig, kkt = max(kk), predB = predB, predL = predL,
    razao = predB / (lamB^2 * WS), bias = bn)
}
nsB <- c(125L, 250L, 500L)
RB <- 40L
resB <- lapply(names(truths), function(nm) {
  lapply(nsB, function(n) t(vapply(seq_len(RB), function(r) one_B(n, truths[[nm]]), numeric(15))))
})
names(resB) <- names(truths)
tabB <- do.call(rbind, lapply(names(truths), function(nm) {
  t(vapply(seq_along(nsB), function(i) {
    rr <- resB[[nm]][[i]]
    c(n = nsB[i], lambda = median(rr[, "lam"]), viola = sum(rr[, "Btv"] > rr[, "bound"]),
      folga = min(rr[, "bound"] / rr[, "Btv"]), viola_l2 = sum(rr[, "l2"] > rr[, "l2b"]),
      viola_pen = sum(rr[, "penv"] > rr[, "penb"]), schur = sum(rr[, "gt"] < rr[, "lsig"] - 1e-10),
      kkt = max(rr[, "kkt"]), razao = median(rr[, "razao"]),
      blk_lasso = median(rr[, "predB"] / rr[, "predL"]))
  }, numeric(10)))
}))
rownames(tabB) <- paste(rep(names(truths), each = length(nsB)), paste0("n=", nsB))
print(signif(tabB, 4))
chk(max(tabB[, "kkt"]) < 1e-7, "o FISTA resolve o objetivo em blocos (KKT a 1e-7 de lambda)")
chk(all(tabB[, "schur"] == 0), "gamma_til >= lambda_min(Sigma_hat) em todas as réplicas: a perfilagem dispensa o cone também em blocos")
chk(all(tabB[, c("viola", "viola_l2", "viola_pen")] == 0),
    "as três cotas do Teorema 1(ii) em blocos (W_S no lugar de s_0) valem em todas as réplicas, com e sem viés")
rzc <- tabB[grep("^cheio", rownames(tabB)), "razao"]
rze <- tabB[grep("^espalhado", rownames(tabB)), "razao"]
cat(sprintf("  razão ||f_hat - f||_n^2 / (lambda^2 W_S): cheio %s (máx/mín %.2f); espalhado %s (máx/mín %.2f)\n",
            paste(sprintf("%.3f", rzc), collapse = " "), max(rzc) / min(rzc),
            paste(sprintf("%.3f", rze), collapse = " "), max(rze) / min(rze)))
chk(max(rzc) / min(rzc) < 2 && max(rze) / min(rze) < 2,
    "o erro de predição escala como lambda^2 W_S ao variar n por um fator 4")
chk(all(tabB[grep("^cheio", rownames(tabB)), "blk_lasso"] < 1),
    "no lambda da teoria, blocos vence o LASSO quando o suporte é de pedaços cheios")
chk(all(tabB[grep("^espalhado", rownames(tabB)), "blk_lasso"] > 1),
    "e perde quando o mesmo número de coeficientes vem um por pedaço")

# Solver conferido contra o glmnet no caso de grupos unitários (LASSO).
set.seed(5)
U <- matrix(runif(300 * q), 300, q)
X <- draw_X(U)
Z <- design(X, U, JB)
y <- as.numeric(Z %*% c(c_true, th_cheio)) + rnorm(300, sd = sig)
rp <- resid_pen(Z)
pe <- prep_fit(rp, y, seq_len(dB), "euclid")
lam_t <- 0.05
x_f <- gl_fit(pe[["G"]], pe[["h"]], seq_len(dB), rep(1, dB), lam_t)
yt <- as.numeric(y - rp[["proj"]](y))
fit_g <- glmnet(rp[["Bt"]], yt, family = "gaussian", lambda = c(4, 2, 1) * lam_t,
                standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-14))
x_g <- as.numeric(fit_g[["beta"]][, 3])
cat(sprintf("  FISTA contra glmnet no LASSO perfilado (lambda = %.2f): desvio máximo %.1e\n", lam_t, max(abs(x_f - x_g))))
chk(max(abs(x_f - x_g)) < 1e-6, "o solver de blocos reproduz o glmnet quando os pedaços são unitários")

# Os pesos e a norma, no melhor lambda de uma grade (n = 500): o que o código
# de fato roda contra o estimador analisado.
cat("  no melhor lambda de uma grade (n = 500, 20 réplicas), razão do erro de predição ao do LASSO:\n")
grpS <- kp_groups(JB, bB, kind = "consecutivo")
metodos <- list(
  lasso      = list(grp = seq_len(dB), w = function(gs) rep(1, length(gs)), norm = "euclid"),
  kp_w1      = list(grp = grpB, w = function(gs) rep(1, length(gs)), norm = "euclid"),
  kp_sqrt    = list(grp = grpB, w = function(gs) sqrt(gs), norm = "euclid"),
  grpreg     = list(grp = grpB, w = function(gs) sqrt(gs), norm = "white"),
  kp_consec  = list(grp = grpS, w = function(gs) rep(1, length(gs)), norm = "euclid"))
mult <- 2^seq(2, -10, by = -1)
best_err <- function(rp, y, f, Z, me, lam_ref) {
  pe <- prep_fit(rp, y, me[["grp"]], me[["norm"]])
  w <- me[["w"]](as.numeric(table(me[["grp"]])))
  x <- NULL
  errs <- vapply(mult, function(mm) {
    x <<- gl_fit(pe[["G"]], pe[["h"]], me[["grp"]], w, mm * lam_ref, x0 = x, tol = 1e-7)
    mean((as.numeric(Z %*% coef_full(rp, y, pe[["back"]](x))) - f)^2)
  }, 1)
  c(err = min(errs), pos = which.min(errs))
}
nmM <- names(metodos)
tabO <- t(vapply(names(truths), function(nm) {
  tr <- truths[[nm]]
  rr <- t(vapply(seq_len(20L), function(r) {
    n <- 500L
    U <- matrix(runif(n * q), n, q)
    X <- draw_X(U)
    Z <- design(X, U, JB)
    fJ <- as.numeric(Z %*% c(c_true, tr[["th"]]))
    f <- if (is.null(tr[["extra"]])) fJ else fJ + X[, 1] * tr[["extra"]](U[, 1])
    y <- f + rnorm(n, sd = sig)
    rp <- resid_pen(Z)
    lam_ref <- lam_lasso(n, dB, sqrt(max(colSums(rp[["B"]]^2) / n)))
    out <- vapply(metodos, function(me) best_err(rp, y, f, Z, me, lam_ref), numeric(2))
    pe <- prep_fit(rp, y, seq_len(dB), "euclid")
    ols <- mean((as.numeric(Z %*% coef_full(rp, y, solve(pe[["G"]], pe[["h"]]))) - f)^2)
    c(out["err", ], ols = ols, topo = sum(out["pos", ] == 1),
      setNames(out["pos", ] == length(mult), paste0("fundo_", nmM)))
  }, numeric(2L * length(metodos) + 2L)))
  c(vapply(c(nmM[-1], "ols"), function(m) median(rr[, m] / rr[, "lasso"]), 1),
    topo = sum(rr[, "topo"]),
    fundo = vapply(nmM, function(m) sum(rr[, paste0("fundo_", m)]), 1))
}, numeric(2L * length(metodos) + 1L)))
print(round(tabO, 3))
chk(all(tabO[, "topo"] == 0), "nenhum método tem o melhor lambda no topo da grade (o ajuste nulo)")
chk(tabO["cheio", "kp_w1"] < 1 && tabO["seno", "kp_w1"] < 1.05,
    "no melhor lambda, blocos com pesos 1 vence o LASSO no suporte de pedaços cheios e não perde no seno")

## ===========================================================================
cat("\nPARTE C. Besov implica a cota do risco ideal por blocos\n")
## ===========================================================================

# Risco ideal por blocos, R(eta) = sum_g min(||theta_g||^2, eta), num bloco
# (l,m) com os níveis j = 0..J-1 e pedaços de até b sem atravessar nível. A
# cota proposta no documento (§3), com a = pi(s + 1/2) > 1 e tau = 1/(s + 1/2):
#   R(eta) <= A_b Cg^tau (eta/b)^{2s/(2s+1)} b^{(2/pi - 1)_+/(2s+1)} + (x_+ + 1) eta,
#   A_b = 2 + 1/(1 - 2^{1-a}) (pi <= 2) ou 2 + 1/(1 - 2^{-2s}) (pi >= 2),
#   x = log2(b^{tau/pi} Cg^tau eta^{-tau/2}) (pi <= 2) ou log2((Cg^2 b/eta)^{1/(2s+1)}) (pi >= 2).
# Três formas de sequência que saturam a Hipótese de Besov por nível,
# ||theta_j.||_pi = Cg 2^{-j(s+1/2-1/pi)}: "uniforme" (mesmo módulo no nível),
# "aleatoria" (direção gaussiana) e "espiga" (a desfavorável aos blocos: um
# coeficiente de altura sqrt(eta) por pedaço, em tantos pedaços quanto o
# orçamento ell_pi do nível permite).
# Contribuição de cada nível ao risco ideal, nas formas "espiga" e "uniforme"
# (analíticas, sem construir o vetor) e "aleatoria" (energias dos pedaços
# calculadas uma vez por (s, pi, b) a partir de uma sequência sorteada).
budget_of <- function(s, pii, Cg, jj) Cg * 2^(-jj * (s + 0.5 - 1 / pii))   # ||theta_j.||_pi
R_levels <- function(s, pii, Cg, Jl, b, eta, shape) {
  vapply(0:(Jl - 1L), function(jj) {
    nj <- 2^jj
    k <- ceiling(nj / b)
    bud <- budget_of(s, pii, Cg, jj)
    if (shape == "espiga") {
      # picos de altura sqrt(eta), um por pedaço: k' picos custam k'^{1/pi} sqrt(eta) em ell_pi
      return(min(k, floor((bud / sqrt(eta))^pii)) * eta)
    }
    cj2 <- (bud / nj^(1 / pii))^2                     # "uniforme": mesmo módulo no nível
    full <- floor(nj / b)
    rest <- nj - full * b
    full * min(b * cj2, eta) + (rest > 0) * min(rest * cj2, eta)
  }, 1)
}
chunk_energy_random <- function(s, pii, Cg, Jl, b, seed) {
  set.seed(seed)
  lapply(0:(Jl - 1L), function(jj) {
    nj <- 2^jj
    a <- rnorm(nj)
    vals <- budget_of(s, pii, Cg, jj) * a / sum(abs(a)^pii)^(1 / pii)
    k <- ceiling(nj / b)
    as.numeric(rowsum(vals^2, rep(seq_len(k), each = b)[seq_len(nj)]))
  })
}
R_ideal <- function(s, pii, Cg, Jl, b, eta, shape, energies = NULL) {
  if (shape == "aleatoria") return(sum(vapply(energies, function(e) sum(pmin(e, eta)), 1)))
  sum(R_levels(s, pii, Cg, Jl, b, eta, shape))
}
bound_C <- function(s, pii, Cg, b, eta) {
  tau <- 1 / (s + 0.5)
  if (pii <= 2) {
    a <- pii * (s + 0.5)
    A <- 2 + 1 / (1 - 2^(1 - a))
    x <- log2(b^(tau / pii) * Cg^tau * eta^(-tau / 2))
  } else {
    A <- 2 + 1 / (1 - 2^(-2 * s))
    x <- log2((Cg^2 * b / eta)^(1 / (2 * s + 1)))
  }
  A * Cg^tau * (eta / b)^(2 * s / (2 * s + 1)) * b^(max(0, 2 / pii - 1) / (2 * s + 1)) +
    (max(x, 0) + 1) * eta
}
set.seed(3)
grid_sp <- expand.grid(s = c(0.8, 1.5, 3), pi = c(1, 1.5, 2, 4))
JC <- 20L
viol <- 0
maxr <- 0
for (i in seq_len(nrow(grid_sp))) {
  s <- grid_sp[i, "s"]
  pii <- grid_sp[i, "pi"]
  for (b in c(1L, 2L, 4L, 8L, 16L, 32L)) {
    en <- chunk_energy_random(s, pii, 1, JC, b, seed = 100L + i)
    for (eta in 10^c(-9, -7, -5, -3)) {
      for (shape in c("uniforme", "aleatoria", "espiga")) {
        R <- R_ideal(s, pii, 1, JC, b, eta, shape, energies = en)
        B <- bound_C(s, pii, 1, b, eta)
        maxr <- max(maxr, R / B)
        if (R > B) viol <- viol + 1
      }
    }
  }
}
cat(sprintf("  %d pares (s, pi), b em {1,...,32}, eta em 1e-9..1e-3, J = %d, três formas: maior razão R/cota = %.3f\n",
            nrow(grid_sp), JC, maxr))
chk(viol == 0, "a cota do risco ideal por blocos vale em todas as sequências de Besov construídas")

# O expoente de b: com eta = b * eps (o lambda^2 de blocos é proporcional a b),
# a forma "espiga" deve crescer como b^{(2/pi - 1)/(2s+1)} quando pi < 2, e a
# "uniforme" não deve crescer (é aí que o pedaço não custa nada). Os níveis
# com 2^j < b, um pedaço cada, ficam de fora: são o termo (x_+ + 1) eta, de
# ordem menor. eps é escolhido por (s, pi) para que a transição entre cabeça e
# cauda fique no nível ~20 com b = 32, longe das duas pontas (J = 40).
JCe <- 40L
bs <- c(8L, 16L, 32L, 64L, 128L)
tabC <- t(vapply(seq_len(nrow(grid_sp)), function(i) {
  s <- grid_sp[i, "s"]
  pii <- grid_sp[i, "pi"]
  a <- pii * (s + 0.5)
  epsC <- 2^(-(20 * a - (1 - pii / 2) * log2(32)) * 2 / pii)
  fine <- function(b) (2^(0:(JCe - 1L))) >= b
  Rs <- vapply(bs, function(b) sum(R_levels(s, pii, 1, JCe, b, b * epsC, "espiga")[fine(b)]), 1)
  Ru <- vapply(bs, function(b) sum(R_levels(s, pii, 1, JCe, b, b * epsC, "uniforme")[fine(b)]), 1)
  c(s = s, pi = pii, previsto = max(0, 2 / pii - 1) / (2 * s + 1),
    espiga = unname(coef(lm(log(Rs) ~ log(bs)))[2]),
    uniforme = unname(coef(lm(log(Ru) ~ log(bs)))[2]))
}, numeric(5)))
print(round(tabC, 3))
lt2 <- tabC[, "pi"] < 2
chk(all(abs(tabC[lt2, "espiga"] - tabC[lt2, "previsto"]) < 0.1),
    "pi < 2: o expoente de b na forma desfavorável é (2/pi - 1)/(2s + 1), a 0.1")
chk(all(tabC[, "uniforme"] < 0.1),
    "sinal espalhado dentro do nível não paga nada pelos pedaços (expoente de b ~ 0)")

## ===========================================================================
cat("\nPARTE D. Os expoentes das taxas\n")
## ===========================================================================

# LASSO (Corolário 5): (log n / n)^{2s/(2s+1)}. Blocos com b ~ log n:
# n^{-2s/(2s+1)} (log n)^{(2/pi - 1)_+/(2s+1)}. O ganho é (log n)^{2s'/(2s+1)},
# s' = s - (1/pi - 1/2)_+, e a janela de c é a mesma do Corolário 5.
gridD <- expand.grid(s = c(0.3, 0.8, 1.5, 3), pi = c(1, 1.5, 2, 4, 10))
gridD <- gridD[gridD[["pi"]] * (gridD[["s"]] + 0.5) > 1, ]
eD <- t(vapply(seq_len(nrow(gridD)), function(i) {
  s <- gridD[i, "s"]
  pii <- gridD[i, "pi"]
  tau <- 1 / (s + 0.5)
  sp <- s - max(0, 1 / pii - 0.5)
  lasso <- 2 * s / (2 * s + 1)                 # expoente de log n no LASSO
  blocos <- max(0, 2 / pii - 1) / (2 * s + 1)  # expoente de log n em blocos
  c(e1 = abs((1 - tau / 2) - 2 * s / (2 * s + 1)),
    e2 = if (pii <= 2) abs(tau * (1 / pii - 0.5) - (2 / pii - 1) / (2 * s + 1)) else 0,
    e3 = abs((lasso - blocos) - 2 * sp / (2 * s + 1)),
    ganho = 2 * sp / (2 * s + 1))
}, numeric(4)))
cat(sprintf("  %d pares (s, pi): identidades a menos de %.1e; ganho (log n)^{2s'/(2s+1)} com expoente entre %.3f e %.3f\n",
            nrow(gridD), max(eD[, c("e1", "e2", "e3")]), min(eD[, "ganho"]), max(eD[, "ganho"])))
chk(max(eD[, c("e1", "e2", "e3")]) < 1e-14,
    "1 - tau/2 = 2s/(2s+1); tau(1/pi - 1/2) = (2/pi - 1)/(2s+1); diferença dos expoentes = 2s'/(2s+1)")
chk(all(eD[, "ganho"] > 0), "o expoente do ganho é positivo sempre que s' > 0")

cat("\n")
if (ok) cat("OK\n") else stop("E1.11: conferência numérica FALHOU (ver linhas acima)")
