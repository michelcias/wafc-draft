# E1.7a. Sondagem numérica da condição de irrepresentabilidade em grupos no
# desenho de produtos (derivations/06a-sondagem-irrepresentabilidade.md).
# É SONDAGEM, não conferência de prova: parte do que sai aqui é medida, não
# asserção. Os testes que são asserção imprimem [ok] ou [FALHA]; o resto sai
# em tabela. Notação: docs/notacao.md (congelada em E1.1).
#
# A condição é a de Bach (2008, JMLR 9, 1179-1225), condição (4), na forma
# forte (suficiente para consistência de padrão); com grupos de tamanho igual
# N_J os pesos d_b se cancelam. Escrita aqui já com as p colunas não
# penalizadas X_l perfiladas por complemento de Schur:
#
#   IRR(b) = || Sig~_{b,S} Sig~_{S,S}^{-1} zeta_S ||_2 ,   b = (l0, m0) fora de S,
#   Sig~_{A,B} = Sig_{A,B} - Sig_{A,C} Sig_{C,C}^{-1} Sig_{C,B},  C = colunas X_l,
#   zeta_{l m} = theta*_{lm} / ||theta*_{lm}||_2  para (l, m) em S.
#
# O que confere e mede:
#   1. sob X ind. de U e U ~ Unif[0,1]^q: Sig = P (Omega (x) I_K) P' e as
#      colunas não penalizadas são ORTOGONAIS às penalizadas (Schur é nulo);
#   2. a redução conjecturada pelo chat principal: IRR(l0, m0) = sqrt(a' V a),
#      com a = Omega_{l0, S_{m0}} Omega_{S_{m0} S_{m0}}^{-1} (só Omega) e
#      V = Gram dos zeta do modulador m0 (só as componentes verdadeiras);
#      logo IRR NÃO DEPENDE DA BASE nem de J, a não ser por V;
#   3. as cotas ||a||_2 sqrt(lmin(V)) <= sqrt(a'Va) <= ||a||_1, e o ganho do
#      agrupamento sobre a condição de Zhao & Yu (2006) para o LASSO, que é
#      ||a||_1: exemplo com ||a||_1 > 1 e sqrt(a'Va) < 1;
#   4. a margem nos cenários de wafc/R/dgp.R e numa extensão densa deles,
#      contra J, contra p q e contra a correlação de X;
#   5. os três obstáculos: U não uniforme, U com coordenadas dependentes, e
#      X dependente de U (cenário B de check/03-desenho-produtos.R). Mede-se
#      o erro da predição por Omega e o VAZAMENTO do oráculo populacional
#      para blocos estruturalmente nulos;
#   6. a versão empírica em n = 250, 500, 1000, 50 réplicas.
#
# A Gram populacional sai por quadratura de ponto médio, com E[XX' | U = u]
# em forma fechada, como em check/03-desenho-produtos.R. A ordem das colunas
# usada aqui é a de Kronecker (l-major), não a de D12: IRR é invariante a
# permutação de colunas, e os índices são montados por conjunto.
#
# Dependências: WaveBased (wbasis, wtable) e wafc/R/ (wafc_component,
# wafc_scenario), carregado por wafc/R/load.R.

suppressPackageStartupMessages(library(WaveBased))
root <- if (file.exists("wafc/R/load.R")) "." else "../.."
source(file.path(root, "wafc/R/load.R"))
set.seed(20260920)

fs <- 8L                                   # Daublets, filtro 8 (4 momentos nulos)
tb <- wtable(filter.size = fs)
ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "ok" else "FALHA", msg))
  if (!isTRUE(cond)) ok <<- FALSE
  invisible(cond)
}
hdr <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")

## ---- base e bookkeeping --------------------------------------------------

# psi_{jk}, j = 0..J-1: wbasis() com j0 = 0 periódico devolve a constante na
# primeira coluna, descartada de cada bloco (D2 / Lema 1 de E1.2).
psi_block <- function(u, J) {
  wbasis(u, j0 = 0, J = J, filter.size = fs, wavelet.table = tb)[, -1, drop = FALSE]
}
# Psi(u) = (1, psi(u_1), ..., psi(u_q)) em R^K, K = 1 + q N_J.
Psi <- function(U, J) {
  cbind(1, do.call(cbind, lapply(seq_len(ncol(U)), function(m) psi_block(U[, m], J))))
}
# Índices na ordem de Kronecker: bloco de X_l ocupa (l-1)K + 1:K, com a
# constante em 1 e o modulador m em 1 + (m-1)N_J + 1:N_J.
ix <- function(p, q, NJ) {
  K <- 1L + q * NJ
  K <- as.integer(K)
  list(K = K, D = p * K,
       C = as.integer(vapply(seq_len(p), function(l) (l - 1L) * K + 1L, 1)),
       B = function(l, m) (l - 1L) * K + 1L + (m - 1L) * NJ + seq_len(NJ))
}

## ---- a estatística IRR ---------------------------------------------------

# Complemento de Schur que perfila as colunas não penalizadas C.
schur_fun <- function(Sig, C) {
  Wcc <- solve(Sig[C, C, drop = FALSE])
  function(A, B) {
    Sig[A, B, drop = FALSE] -
      Sig[A, C, drop = FALSE] %*% Wcc %*% Sig[C, B, drop = FALSE]
  }
}
# IRR(b) para cada bloco candidato, dado o suporte S (lista de vetores de
# índices) e as direções zeta (lista de vetores unitários de R^{N_J}).
irr_stat <- function(Sig, C, S_idx, zeta, cand_idx) {
  sc <- schur_fun(Sig, C)
  Sidx <- unlist(S_idx, use.names = FALSE)
  w <- solve(sc(Sidx, Sidx), unlist(zeta, use.names = FALSE))
  vapply(cand_idx, function(bi) sqrt(sum((sc(bi, Sidx) %*% w)^2)), 0)
}

## ---- coeficientes verdadeiros e direções zeta ----------------------------

# theta^Leb_{lm} = (<g_{lm}, psi_{jk}>)_{jk} por quadratura 1-D fina. Sob U
# uniforme e X ind. de U este é o oráculo populacional (verificado na §1).
G1 <- 2L^14L
g1grid <- (seq_len(G1) - 0.5) / G1
theta_leb <- function(g, J) {
  P <- psi_block(g1grid, J)
  as.numeric(crossprod(P, g(g1grid)) / G1)
}
# Estrutura de um cenário: matriz p x q de nomes de componentes.
#   "dgp"      a de wafc_scenario(): (1,1), (1,2) e (2,1) ativos;
#   "densa"    quatro covariáveis lineares dividem o modulador 1, que é o
#              regime em que |S_m| > 2 e a condição pode morder;
#   "alinhada" a mesma, com as quatro componentes IGUAIS, que é o pior caso
#              de V (todas as direções colineares, V = 1 1').
struct_of <- function(scenario, p, q, mode = c("dgp", "densa", "alinhada"),
                      na = NULL) {
  mode <- match.arg(mode)
  st <- wafc_scenario(scenario, p, q)
  if (mode == "dgp") return(st)
  nm <- if (scenario == "smooth") c("sine", "cubic", "cosine")
        else c("bumps", "blocks", "heavisine")
  st[] <- ""
  # ativos no modulador 1; o padrão deixa ao menos um candidato fora de S
  if (is.null(na)) na <- min(4L, max(1L, p - 2L))
  for (l in seq_len(na)) {
    st[l, 1L] <- if (mode == "alinhada") nm[1L] else nm[(l - 1L) %% 3L + 1L]
  }
  if (q >= 2L) st[1L, 2L] <- nm[2L]
  st
}
# which(..., arr.ind = TRUE) precisa da matriz: nzchar() devolve vetor sem dim.
pairs_of <- function(st, active = TRUE) {
  m <- matrix(nzchar(st), nrow(st), ncol(st))
  which(if (active) m else !m, arr.ind = TRUE)
}
zeta_of <- function(st, J, amp = NULL) {
  p <- nrow(st); q <- ncol(st)
  Z <- vector("list", p * q); dim(Z) <- c(p, q)
  TH <- vector("list", p * q); dim(TH) <- c(p, q)
  for (l in seq_len(p)) for (m in seq_len(q)) {
    if (!nzchar(st[l, m])) next
    a <- if (is.null(amp)) 1 else amp[l, m]
    th <- a * theta_leb(wafc_component(st[l, m]), J)
    TH[[l, m]] <- th
    Z[[l, m]] <- th / sqrt(sum(th^2))
  }
  list(zeta = Z, theta = TH)
}

## ---- três famílias de Omega = E(XX') -------------------------------------

omega_ar <- function(p, rho) {
  Om <- diag(p); Om[1, 1] <- 1
  if (p >= 2L) { k <- p - 1L
    Om[2:p, 2:p] <- outer(seq_len(k), seq_len(k), function(i, j) rho^abs(i - j)) }
  Om
}
omega_equi <- function(p, rho) {
  Om <- diag(p); Om[1, 1] <- 1
  if (p >= 2L) { k <- p - 1L
    E <- matrix(rho, k, k); diag(E) <- 1
    Om[2:p, 2:p] <- E }
  Om
}
omega_adv <- function(p, rho, act, r12 = 0.2) {
  # act: índices (em 1..p) das covariáveis cujo bloco no modulador 1 é ativo.
  Om <- diag(p); Om[1, 1] <- 1
  A <- setdiff(act, 1L); B <- setdiff(2:p, A)
  if (length(A) > 1L) {
    E <- matrix(r12, length(A), length(A)); diag(E) <- 1
    Om[A, A] <- E
  }
  if (length(A) && length(B)) { Om[B, A] <- rho; Om[A, B] <- rho }
  Om
}
## =========================================================================
## 1. A estrutura da Gram sob X ind. de U e U uniforme
## =========================================================================
hdr("1. Gram sob X independente de U e U ~ Unif[0,1]^2 (quadratura)")

p <- 4L; q <- 2L; J <- 3L
NJ <- 2L^J - 1L
I4 <- ix(p, q, NJ)
# Omega = E(XX'): X_1 = 1, (X_2, X_3, X_4) gaussianas centradas com
# correlação AR(1) rho = 0.6.
# Omega adversária: X_2 e X_3 (blocos ativos em U_1) correlacionados 0.2, e
# X_4 (candidato) correlacionado 0.7 com os dois. É o desenho em que ||a||_1
# passa de 1 e mais de uma direção zeta contribui, que é o caso que interessa.
Omega <- omega_adv(p, 0.7, act = 1:3)

G <- 400L
gg <- (seq_len(G) - 0.5) / G
Ug <- as.matrix(expand.grid(u1 = gg, u2 = gg))
Pg <- Psi(Ug, J)
# Sigma = E[ E(XX'|U) (x) Psi Psi' ]; aqui E(XX'|U) = Omega (independência).
pop_gram_indep <- function(Omega, Pg, wts = NULL) {
  K <- ncol(Pg); p <- nrow(Omega)
  A <- if (is.null(wts)) crossprod(Pg) / nrow(Pg) else crossprod(Pg, wts * Pg) / sum(wts) * 1
  kronecker(Omega, A)
}
Sig_quad <- pop_gram_indep(Omega, Pg)
Sig_kron <- kronecker(Omega, diag(I4$K))
chk(max(abs(Sig_quad - Sig_kron)) < 1e-4,
    sprintf("Sigma = Omega (x) I_K a %.1e (erro de quadratura)",
            max(abs(Sig_quad - Sig_kron))))
cross <- Sig_quad[I4$C, setdiff(seq_len(I4$D), I4$C), drop = FALSE]
chk(max(abs(cross)) < 1e-4,
    sprintf("colunas não penalizadas ortogonais às penalizadas (max %.1e): o complemento de Schur é nulo",
            max(abs(cross))))

## =========================================================================
## 2. A redução conjecturada: IRR = sqrt(a' V a), com a vindo só de Omega
## =========================================================================
hdr("2. Redução a Omega: IRR(l0,m0) = sqrt(a' V a)")

# Predição por Omega e pelas direções, sem tocar na base.
irr_pred <- function(Omega, Zt, l0, m0, Sm) {
  a <- as.numeric(Omega[l0, Sm, drop = FALSE] %*% solve(Omega[Sm, Sm, drop = FALSE]))
  Vm <- outer(Sm, Sm, Vectorize(function(i, j) sum(Zt[[i, m0]] * Zt[[j, m0]])))
  c(irr = sqrt(max(0, drop(a %*% Vm %*% a))), l1 = sum(abs(a)), l2 = sqrt(sum(a^2)),
    lminV = min(eigen(Vm, symmetric = TRUE, only.values = TRUE)$values))
}
# Suporte: cenário suave denso com p = 4 -> S_1 = {1,2,3,4} menos o candidato.
st <- struct_of("smooth", p, q, mode = "densa", na = 3L)  # (4,1) fica fora de S
zz <- zeta_of(st, J)
S_pairs <- pairs_of(st)
S_idx <- lapply(seq_len(nrow(S_pairs)),
                function(r) I4$B(S_pairs[r, 1], S_pairs[r, 2]))
zeta_S <- lapply(seq_len(nrow(S_pairs)),
                 function(r) zz$zeta[[S_pairs[r, 1], S_pairs[r, 2]]])
cand <- pairs_of(st, FALSE)
cand_idx <- lapply(seq_len(nrow(cand)), function(r) I4$B(cand[r, 1], cand[r, 2]))
irr_num <- irr_stat(Sig_quad, I4$C, S_idx, zeta_S, cand_idx)
irr_prd <- vapply(seq_len(nrow(cand)), function(r) {
  l0 <- cand[r, 1]; m0 <- cand[r, 2]
  Sm <- which(nzchar(st[, m0]))
  if (length(Sm) == 0L) return(0)
  irr_pred(Omega, zz$zeta, l0, m0, Sm)[["irr"]]
}, 0)
tab <- data.frame(bloco = sprintf("(%d,%d)", cand[, 1], cand[, 2]),
                  IRR_Sigma = round(irr_num, 6), IRR_Omega = round(irr_prd, 6))
print(tab, row.names = FALSE)
chk(max(abs(irr_num - irr_prd)) < 1e-4,
    sprintf("a redução vale exatamente (max |diferença| = %.1e)", max(abs(irr_num - irr_prd))))

# Invariância em J: só V muda, e V converge para as correlações L_2 de g.
hdr("2b. Dependência em J (U uniforme): só através de V")
Jv <- 2:7
rowsJ <- do.call(rbind, lapply(Jv, function(Jx) {
  zx <- zeta_of(st, Jx)
  r <- vapply(seq_len(nrow(cand)), function(rr) {
    l0 <- cand[rr, 1]; m0 <- cand[rr, 2]; Sm <- which(nzchar(st[, m0]))
    if (!length(Sm)) return(0)
    irr_pred(Omega, zx$zeta, l0, m0, Sm)[["irr"]]
  }, 0)
  data.frame(J = Jx, NJ = 2^Jx - 1, IRR_max = max(r))
}))
print(rowsJ, row.names = FALSE)
chk(diff(range(rowsJ$IRR_max)) < 0.05,
    sprintf("IRR estável em J = 2..7 (amplitude %.4f), como a redução prevê",
            diff(range(rowsJ$IRR_max))))

## ---- 2c. o núcleo algébrico, com matrizes aleatórias ---------------------
hdr("2c. Onde a redução vale: Sigma~ = Omega (x) T com T bloco-diagonal em m")

# Sigma~ (já perfilada) é Omega (x) T, com T = Cov(psi_{jk}(U_m), psi_{j'k'}(U_m')),
# de tamanho q N_J. Se T for bloco-diagonal no modulador, o sistema desacopla em
# m e T_m CANCELA: a redução vale para T_m arbitrário, não só para T_m = I.
rspd <- function(k, df = k + 3L) { A <- matrix(rnorm(df * k), df, k); crossprod(A) / df + 0.3 * diag(k) }
red_test <- function(Tm) {
  pz <- 5L; qz <- 3L; NJz <- 4L
  Om <- rspd(pz)
  Bi <- function(l, m) (l - 1L) * qz * NJz + (m - 1L) * NJz + seq_len(NJz)
  Sg <- kronecker(Om, Tm)
  Sl <- list(c(1, 1), c(2, 1), c(3, 1), c(1, 2), c(4, 3)); b <- c(4, 1)
  zt <- lapply(Sl, function(s) { v <- rnorm(NJz); v / sqrt(sum(v^2)) })
  Si <- unlist(lapply(Sl, function(s) Bi(s[1], s[2])))
  w <- solve(Sg[Si, Si], unlist(zt))
  lhs <- sqrt(sum((Sg[Bi(b[1], b[2]), Si, drop = FALSE] %*% w)^2))
  Sm <- 1:3
  a <- as.numeric(Om[4, Sm, drop = FALSE] %*% solve(Om[Sm, Sm, drop = FALSE]))
  Vm <- outer(Sm, Sm, Vectorize(function(i, j) sum(zt[[i]] * zt[[j]])))
  c(exato = lhs, predito = sqrt(drop(a %*% Vm %*% a)))
}
set.seed(11)
tt <- rbind(
  `T = I`                   = red_test(diag(12)),
  `T = D (x) I`             = red_test(kronecker(diag(c(1, 2, 0.5)), diag(4))),
  `T bloco-diagonal em m`   = red_test(as.matrix(Matrix::bdiag(rspd(4), rspd(4), rspd(4)))),
  `T cheia (m acoplados)`   = red_test(rspd(12)))
print(round(cbind(tt, erro = tt[, 1] - tt[, 2]), 6))
chk(max(abs(tt[1:3, 1] - tt[1:3, 2])) < 1e-10,
    "a redução é EXATA sempre que T é bloco-diagonal no modulador, com T_m qualquer")
chk(abs(tt[4, 1] - tt[4, 2]) > 1e-3,
    "e falha quando T acopla moduladores diferentes: é aí que a base volta a entrar")

# Tamanho do erro com T cheia, sobre suportes aleatórios: a falha não é
# patológica, mas também não é uniforme (há suportes em que ela não aparece).
red_rand <- function(bloco_diag, nrep = 300L) {
  pz <- 6L; qz <- 3L; NJz <- 4L
  vapply(seq_len(nrep), function(r) {
    Om <- rspd(pz)
    Tm <- if (bloco_diag) as.matrix(Matrix::bdiag(rspd(NJz), rspd(NJz), rspd(NJz)))
          else rspd(qz * NJz)
    Bi <- function(l, m) (l - 1L) * qz * NJz + (m - 1L) * NJz + seq_len(NJz)
    act <- which(matrix(runif(pz * qz) < 0.35, pz, qz), arr.ind = TRUE)
    if (nrow(act) < 2L) return(NA_real_)
    cnd <- which(matrix(!(seq_len(pz * qz) %in% ((act[, 2] - 1L) * pz + act[, 1])),
                        pz, qz), arr.ind = TRUE)
    m0s <- cnd[cnd[, 2] %in% act[, 2], , drop = FALSE]
    if (!nrow(m0s)) return(NA_real_)
    b <- m0s[sample.int(nrow(m0s), 1L), ]
    Sl <- lapply(seq_len(nrow(act)), function(i) act[i, ])
    zt <- lapply(Sl, function(s) { v <- rnorm(NJz); v / sqrt(sum(v^2)) })
    Si <- unlist(lapply(Sl, function(s) Bi(s[1], s[2])))
    Sg <- kronecker(Om, Tm)
    w <- tryCatch(solve(Sg[Si, Si], unlist(zt)), error = function(e) NULL)
    if (is.null(w)) return(NA_real_)
    lhs <- sqrt(sum((Sg[Bi(b[1], b[2]), Si, drop = FALSE] %*% w)^2))
    Sm <- act[act[, 2] == b[2], 1]
    a <- as.numeric(Om[b[1], Sm, drop = FALSE] %*% solve(Om[Sm, Sm, drop = FALSE]))
    zm <- zt[act[, 2] == b[2]]
    Vm <- outer(seq_along(zm), seq_along(zm), Vectorize(function(i, j) sum(zm[[i]] * zm[[j]])))
    lhs - sqrt(max(0, drop(a %*% Vm %*% a)))
  }, 0)
}
set.seed(12)
e_bd <- red_rand(TRUE); e_full <- red_rand(FALSE)
cat(sprintf("\n  suportes aleatórios (p=6, q=3, N_J=4), |erro da redução|:\n"))
cat(sprintf("    T bloco-diagonal: max %.2e sobre %d sorteios\n",
            max(abs(e_bd), na.rm = TRUE), sum(!is.na(e_bd))))
cat(sprintf("    T cheia:          mediana %.3f, q90 %.3f, max %.3f; |erro| < 1e-8 em %.0f%% dos sorteios\n",
            median(abs(e_full), na.rm = TRUE), quantile(abs(e_full), .9, na.rm = TRUE),
            max(abs(e_full), na.rm = TRUE), 100 * mean(abs(e_full) < 1e-8, na.rm = TRUE)))

## =========================================================================
## 3. O ganho do agrupamento sobre a condição do LASSO (Zhao & Yu 2006)
## =========================================================================
hdr("3. sqrt(a'Va) contra ||a||_1 (LASSO) e ||a||_2 (melhor caso)")

# Desenho adversário: X_1 = 1; X_2, X_3 ativos com correlação r12; X_4 inativo,
# correlacionado r com os dois. Para o LASSO, ||a||_1 = 2r/(1+r12) > 1 já com
# r = 0.7, r12 = 0.2; para o grupo, entra o cosseno entre as duas componentes.
adv <- function(r, r12, cosv) {
  Om <- diag(4); Om[1, 1] <- 1; Om[1, 2:4] <- 0; Om[2:4, 1] <- 0
  Om[2, 3] <- Om[3, 2] <- r12
  Om[2, 4] <- Om[4, 2] <- r
  Om[3, 4] <- Om[4, 3] <- r
  a <- as.numeric(Om[4, 2:3, drop = FALSE] %*% solve(Om[2:3, 2:3, drop = FALSE]))
  Vm <- matrix(c(1, cosv, cosv, 1), 2, 2)
  c(l1 = sum(abs(a)), grupo = sqrt(drop(a %*% Vm %*% a)), l2 = sqrt(sum(a^2)))
}
grid_adv <- expand.grid(r = c(0.5, 0.6, 0.7, 0.8), cosv = c(0, 0.3, 0.6, 0.9))
res_adv <- t(apply(grid_adv, 1, function(z) adv(z[["r"]], 0.2, z[["cosv"]])))
print(cbind(grid_adv, round(res_adv, 4)), row.names = FALSE)
i0 <- which(grid_adv$r == 0.7 & grid_adv$cosv == 0)
chk(res_adv[i0, "l1"] > 1 && res_adv[i0, "grupo"] < 1,
    sprintf("existe desenho com ||a||_1 = %.3f > 1 (LASSO falha) e sqrt(a'Va) = %.3f < 1 (grupos vale)",
            res_adv[i0, "l1"], res_adv[i0, "grupo"]))

# O cosseno é a correlação L_2 entre as componentes que dividem o modulador.
hdr("3b. Cossenos <zeta_{l m}, zeta_{l' m}> das componentes de dgp.R (J = 5)")
for (sc in c("smooth", "inhomogeneous")) {
  nms <- if (sc == "smooth") c("sine", "cubic", "cosine") else c("bumps", "blocks", "heavisine")
  th <- lapply(nms, function(nm) { t <- theta_leb(wafc_component(nm), 5L); t / sqrt(sum(t^2)) })
  Vm <- outer(1:3, 1:3, Vectorize(function(i, j) sum(th[[i]] * th[[j]])))
  cat(" ", sc, ":", paste(sprintf("<%s,%s> = %+.4f", nms[c(1, 1, 2)], nms[c(2, 3, 3)],
                                  Vm[cbind(c(1, 1, 2), c(2, 3, 3))]), collapse = "; "), "\n")
}


## =========================================================================
## 4. A margem nos cenários, e como ela degrada em J, em p q e em Omega
## =========================================================================
hdr("4. Margem eta = 1 - max_b IRR(b) (U uniforme, X ind. de U, J = 4)")

# Três famílias de Omega = E(XX'), com X_1 = 1 e X_2..X_p de variância 1.
#   ar    AR(1) de parâmetro rho: Markov, a regressão de um inativo sobre os
#         ativos carrega um coeficiente só, e a condição quase não morde;
#   equi  equicorrelacionada rho: ||a||_1 = s rho/(1 + (s-1) rho) < 1 sempre;
#   adv   ativos equicorrelacionados r12 = 0.2 e cada inativo correlacionado
#         rho com TODOS os ativos: é a construção de Zhao & Yu (2006) em que
#         ||a||_1 > 1 e a condição do LASSO falha.
margin_of <- function(st, Omega, J, zx = NULL) {
  if (is.null(zx)) zx <- zeta_of(st, J)
  cd <- pairs_of(st, FALSE)
  if (!nrow(cd)) return(c(irr = 0, l1 = 0, s = 0))
  r <- vapply(seq_len(nrow(cd)), function(rr) {
    l0 <- cd[rr, 1]; m0 <- cd[rr, 2]; Sm <- which(nzchar(st[, m0]))
    if (!length(Sm)) return(c(0, 0))
    v <- irr_pred(Omega, zx$zeta, l0, m0, Sm)
    c(v[["irr"]], v[["l1"]])
  }, c(0, 0))
  c(irr = max(r[1, ]), l1 = max(r[2, ]),
    s = max(colSums(matrix(nzchar(st), nrow(st), ncol(st)))))
}
grid4 <- expand.grid(cenario = c("smooth", "inhomogeneous"),
                     forma = c("dgp", "densa", "alinhada"),
                     familia = c("ar", "equi", "adv"),
                     rho = c(0.3, 0.6, 0.8), stringsAsFactors = FALSE)
res4 <- do.call(rbind, lapply(seq_len(nrow(grid4)), function(r) {
  z <- grid4[r, ]; pz <- 8L
  stx <- struct_of(z$cenario, pz, 2L, mode = z$forma)
  act <- which(nzchar(stx[, 1L]))
  Om <- switch(z$familia, ar = omega_ar(pz, z$rho), equi = omega_equi(pz, z$rho),
               adv = omega_adv(pz, z$rho, act))
  m <- margin_of(stx, Om, 4L)
  data.frame(z, S1 = length(act), IRR = round(m[["irr"]], 4),
             LASSO_l1 = round(m[["l1"]], 4), eta = round(1 - m[["irr"]], 4))
}))
print(res4[order(res4$familia, res4$cenario, res4$forma, res4$rho), ], row.names = FALSE)
cat(sprintf("\n  falha em grupos (IRR > 1): %d de %d; falha no LASSO (||a||_1 > 1): %d de %d\n",
            sum(res4$IRR > 1), nrow(res4), sum(res4$LASSO_l1 > 1), nrow(res4)))

hdr("4b. Degradação em p q (familia adv, rho = 0.7, forma densa, J = 4)")
res4b <- do.call(rbind, lapply(c(4L, 6L, 8L, 12L, 20L, 40L), function(pz) {
  do.call(rbind, lapply(c("smooth", "inhomogeneous"), function(sc) {
    stx <- struct_of(sc, pz, 2L, mode = "densa")
    Om <- omega_adv(pz, 0.7, which(nzchar(stx[, 1L])))
    m <- margin_of(stx, Om, 4L)
    data.frame(cenario = sc, p = pz, pq = 2L * pz, S1 = sum(nzchar(stx[, 1L])),
               IRR = round(m[["irr"]], 4), LASSO_l1 = round(m[["l1"]], 4),
               eta = round(1 - m[["irr"]], 4))
  }))
}))
print(res4b[order(res4b$cenario, res4b$p), ], row.names = FALSE)
chk(diff(range(res4b$IRR[res4b$cenario == "smooth" & res4b$p >= 6])) < 1e-8,
    "IRR não depende de p q uma vez fixado o suporte: só o número de candidatos cresce")

hdr("4b'. Degradação em p q com Omega aleatória (modelo de um fator, 200 sorteios)")

# Omega adversária é simétrica nos candidatos, então o max sobre b não cresce
# com p q. Com Omega aleatória cresce, porque o max é sobre p q - |S| blocos.
omega_fator <- function(p, carga) {
  k <- p - 1L
  L <- carga * matrix(rnorm(k), k, 1L)
  S <- tcrossprod(L) + (1 - carga^2) * diag(k)
  d <- sqrt(diag(S)); S <- S / outer(d, d)
  Om <- diag(p); Om[1, 1] <- 1; Om[2:p, 2:p] <- S
  Om
}
set.seed(4)
res4br <- do.call(rbind, lapply(c(4L, 6L, 10L, 20L, 40L, 80L), function(pz) {
  stx <- struct_of("inhomogeneous", pz, 2L, mode = "densa")
  zx <- zeta_of(stx, 4L)
  do.call(rbind, lapply(c(0.9, 0.97), function(cg) {
    v <- replicate(200L, margin_of(stx, omega_fator(pz, cg), 4L, zx)[["irr"]])
    data.frame(carga = cg, p = pz, pq = 2L * pz, cand = sum(!nzchar(stx)),
               IRR_mediana = round(median(v), 4),
               IRR_q90 = round(unname(quantile(v, .9)), 4),
               falhas_pct = round(100 * mean(v > 1), 1))
  }))
}))
print(res4br[order(res4br$carga, res4br$p), ], row.names = FALSE)

hdr("4c. Degradação em J (forma densa e alinhada, familia adv, rho = 0.7, p = 8)")
res4c <- do.call(rbind, lapply(2:7, function(Jx) {
  do.call(rbind, lapply(c("densa", "alinhada"), function(fm) {
    do.call(rbind, lapply(c("smooth", "inhomogeneous"), function(sc) {
      stx <- struct_of(sc, 8L, 2L, mode = fm)
      Om <- omega_adv(8L, 0.7, which(nzchar(stx[, 1L])))
      m <- margin_of(stx, Om, Jx)
      data.frame(J = Jx, forma = fm, cenario = sc, IRR = round(m[["irr"]], 4))
    }))
  }))
}))
print(reshape(res4c, idvar = c("forma", "cenario"), timevar = "J", direction = "wide"),
      row.names = FALSE)

## =========================================================================
## 5. Os três obstáculos à redução
## =========================================================================
hdr("5. Obstáculos: U não uniforme, U dependente, X dependente de U")

# Gram e oráculo populacionais na forma geral, por quadratura em [0,1]^2.
# Mfun(u) devolve a lista p x p de E(X_l X_l' | U = u); dens(u) a densidade.
pop_general <- function(Mfun, dens, Ug, Pg, p, beta_fun) {
  K <- ncol(Pg)
  w <- dens(Ug); w <- w / sum(w)
  M <- Mfun(Ug)
  Sig <- matrix(0, p * K, p * K)
  Om <- matrix(0, p, p)
  for (l in seq_len(p)) for (lp in l:p) {
    A <- crossprod(Pg, (w * M[[l, lp]]) * Pg)
    Sig[(l - 1) * K + 1:K, (lp - 1) * K + 1:K] <- A
    if (lp > l) Sig[(lp - 1) * K + 1:K, (l - 1) * K + 1:K] <- t(A)
    Om[l, lp] <- Om[lp, l] <- sum(w * M[[l, lp]])
  }
  Bt <- beta_fun(Ug)                                      # n x p: beta_l(u)
  b <- numeric(p * K)
  for (l in seq_len(p)) {
    s <- numeric(nrow(Ug))
    for (lp in seq_len(p)) s <- s + M[[l, lp]] * Bt[, lp]
    b[(l - 1) * K + 1:K] <- crossprod(Pg, w * s)
  }
  list(Sigma = Sig, Omega = Om, theta = solve(Sig, b))
}
# X_1 = 1; X_l = mu + eta_l para l >= 3 e X_2 = mu + tau sqrt(12)(u_1 - 1/2) +
# sqrt(1 - tau^2) eta_2, com (eta_2, ..., eta_p) ~ N(0, R). Var(X_l) = 1 para
# todo tau. A média mu != 0 é o que deixa Omega_{1,l} != 0 e, com isso, faz o
# acoplamento entre moduladores chegar ao candidato: com X centrado o bloco
# (1, m') não alcança (l_0, m_0) e o obstáculo fica invisível.
mk_M <- function(R, tau, mu = 0.5) {
  p <- nrow(R) + 1L
  function(u) {
    n <- nrow(u)
    mu2 <- mu + tau * sqrt(12) * (u[, 1] - 0.5)
    M <- vector("list", p * p); dim(M) <- c(p, p)
    S <- (1 - tau^2) * R
    for (l in seq_len(p)) for (lp in seq_len(p)) {
      M[[l, lp]] <- rep(
        if (l == 1L && lp == 1L) 1
        else if (l == 1L || lp == 1L) mu
        else S[l - 1L, lp - 1L] + mu^2, n)
    }
    M[[1, 2]] <- M[[2, 1]] <- mu2
    M[[2, 2]] <- mu2^2 + S[1, 1]
    for (lp in seq_len(p)[-(1:2)]) M[[2, lp]] <- M[[lp, 2]] <- mu2 * mu + S[1, lp - 1L]
    M
  }
}

p5 <- 4L; J5 <- 3L; NJ5 <- 2L^J5 - 1L; I5 <- ix(p5, 2L, NJ5)
G5 <- 256L; g5 <- (seq_len(G5) - 0.5) / G5
U5 <- as.matrix(expand.grid(u1 = g5, u2 = g5)); P5 <- Psi(U5, J5)
st5 <- struct_of("smooth", p5, 2L, mode = "densa", na = 3L)   # candidato: (4,1)
z5 <- zeta_of(st5, J5)
S5p <- pairs_of(st5)
S5 <- lapply(seq_len(nrow(S5p)), function(r) I5$B(S5p[r, 1], S5p[r, 2]))
zt5 <- lapply(seq_len(nrow(S5p)), function(r) z5$zeta[[S5p[r, 1], S5p[r, 2]]])
cd5 <- pairs_of(st5, FALSE)
ci5 <- lapply(seq_len(nrow(cd5)), function(r) I5$B(cd5[r, 1], cd5[r, 2]))
R5 <- omega_adv(p5, 0.7, which(nzchar(st5[, 1L])))[2:p5, 2:p5, drop = FALSE]
cc5 <- rep_len(c(1, 2, -1.5, 0.5), p5)
beta5 <- function(u) {
  B <- matrix(rep(cc5, each = nrow(u)), nrow(u), p5)
  for (l in seq_len(p5)) for (m in 1:2) if (nzchar(st5[l, m]))
    B[, l] <- B[, l] + wafc_component(st5[l, m])(u[, m])
  B
}
dens_unif <- function(u) rep(1, nrow(u))
dens_beta <- function(a, b) function(u) dbeta(u[, 2], a, b)
dens_cop <- function(r) function(u) {
  z <- qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
  exp(-(r^2 * (z[, 1]^2 + z[, 2]^2) - 2 * r * z[, 1] * z[, 2]) / (2 * (1 - r^2))) /
    sqrt(1 - r^2)
}
# T = Cov(psi_{jk}(U_m), psi_{j'k'}(U_m')): o bloco l = 1 de Sigma é Sigma_Psi
# (X_1 = 1), e T é o seu complemento de Schur em relação à constante.
cross_mod <- function(Sig, K, q, NJ) {
  Sp <- Sig[1:K, 1:K, drop = FALSE]
  T <- Sp[-1, -1, drop = FALSE] - outer(Sp[-1, 1], Sp[1, -1]) / Sp[1, 1]
  if (q < 2L) return(0)
  mx <- 0
  for (m in 1:(q - 1L)) for (mp in (m + 1L):q) {
    mx <- max(mx, max(abs(T[(m - 1L) * NJ + seq_len(NJ), (mp - 1L) * NJ + seq_len(NJ)])))
  }
  mx / max(abs(diag(T)))
}
probe5 <- function(nome, Mfun, dens) {
  pg <- pop_general(Mfun, dens, U5, P5, p5, beta5)
  nrm <- function(idx) sqrt(sum(pg$theta[idx]^2))
  # zeta VEM DO ORÁCULO POPULACIONAL deste caso, não da projeção em Lebesgue.
  zt <- lapply(S5, function(idx) pg$theta[idx] / nrm(idx))
  zl <- vector("list", p5 * 2L); dim(zl) <- c(p5, 2L)
  for (r in seq_len(nrow(S5p))) zl[[S5p[r, 1], S5p[r, 2]]] <- zt[[r]]
  irr <- max(irr_stat(pg$Sigma, I5$C, S5, zt, ci5))
  prd <- max(vapply(seq_len(nrow(cd5)), function(rr) {
    l0 <- cd5[rr, 1]; m0 <- cd5[rr, 2]; Sm <- which(nzchar(st5[, m0]))
    if (!length(Sm)) return(0)
    irr_pred(pg$Omega, zl, l0, m0, Sm)[["irr"]]
  }, 0))
  data.frame(caso = nome, IRR_Sigma = round(irr, 4), IRR_Omega = round(prd, 4),
             erro = signif(irr - prd, 3),
             cruzM = signif(cross_mod(pg$Sigma, I5$K, 2L, NJ5), 3),
             cruz_CS = round(max(abs(pg$Sigma[I5$C, unlist(S5)])), 4),
             vazamento = signif(max(vapply(ci5, nrm, 0)) / min(vapply(S5, nrm, 0)), 3))
}
res5 <- rbind(
  probe5("referência: U unif., X ind. de U", mk_M(R5, 0), dens_unif),
  probe5("U_2 ~ Beta(2,3)", mk_M(R5, 0), dens_beta(2, 3)),
  probe5("U_2 ~ Beta(0.6,0.6) (densidade em U)", mk_M(R5, 0), dens_beta(0.6, 0.6)),
  probe5("U_1, U_2 cópula gaussiana 0.3", mk_M(R5, 0), dens_cop(0.3)),
  probe5("U_1, U_2 cópula gaussiana 0.6", mk_M(R5, 0), dens_cop(0.6)),
  probe5("X_2 depende de U_1, tau = 0.3", mk_M(R5, 0.3), dens_unif),
  probe5("X_2 depende de U_1, tau = 0.6", mk_M(R5, 0.6), dens_unif),
  probe5("X_2 depende de U_1, tau = 0.9", mk_M(R5, 0.9), dens_unif))
print(res5, row.names = FALSE)
cat("\n  cruzM:   maior |T| entre moduladores distintos, relativo à diagonal de T\n",
    " cruz_CS: maior |Sigma| entre coluna não penalizada e coluna de S (0 = Schur trivial)\n",
    " vazamento: ||theta*|| do maior bloco estruturalmente nulo / do menor bloco ativo\n", sep = "")

## =========================================================================
## 6. Versão empírica
## =========================================================================
hdr("6. IRR na Gram empírica (U uniforme, X ind. de U, Omega adversária 0.7)")

emp <- function(n, nrep = 50L) {
  Rk <- chol(R5)
  v <- numeric(nrep); s <- numeric(nrep)
  for (r in seq_len(nrep)) {
    U <- matrix(runif(n * 2L), n, 2L)
    X <- cbind(1, 0.5 + matrix(rnorm(n * (p5 - 1L)), n, p5 - 1L) %*% Rk)
    Pn <- Psi(U, J5)
    Zn <- X[, rep(seq_len(p5), each = ncol(Pn))] * Pn[, rep(seq_len(ncol(Pn)), p5)]
    Sh <- crossprod(Zn) / n
    v[r] <- max(irr_stat(Sh, I5$C, S5, zt5, ci5))
    s[r] <- min(eigen(Sh, symmetric = TRUE, only.values = TRUE)$values)
  }
  list(v = v, s = s)
}
pop6 <- max(irr_stat(pop_general(mk_M(R5, 0), dens_unif, U5, P5, p5, beta5)$Sigma,
                     I5$C, S5, zt5, ci5))
res6 <- do.call(rbind, lapply(c(250L, 500L, 1000L, 4000L), function(n) {
  e <- emp(n)
  data.frame(n = n, D = I5$D, mediana = round(median(e$v), 4),
             q90 = round(unname(quantile(e$v, 0.9)), 4), maximo = round(max(e$v), 4),
             falhas = sum(e$v > 1), lmin_hat = round(median(e$s), 4))
}))
res6$populacional <- round(pop6, 4)
print(res6, row.names = FALSE)
chk(abs(res6$mediana[nrow(res6)] - pop6) < 0.05,
    sprintf("a mediana empírica converge para o valor populacional (%.4f contra %.4f em n = 4000)",
            res6$mediana[nrow(res6)], pop6))

cat("\n", if (ok) "OK" else "FALHOU", "\n", sep = "")
if (!ok) stop("Alguma asserção da sondagem falhou.")
