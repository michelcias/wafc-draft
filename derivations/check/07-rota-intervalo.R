## E1.8 -- conferência numérica da rota do intervalo, do caso eps > 0 e da
## análise de s' (derivations/07-rota-intervalo.tex).
##
## Roda com
##   Rscript derivations/check/07-rota-intervalo.R
## e termina imprimindo OK. Leva cerca de 3 min (a parte C constrói a base
## até o nível 12 em três regimes).
##
## O que se confere, na ordem do .tex:
##
## Parte A (a rota do intervalo, base de Cohen-Daubechies-Vial):
##   (A1) a base CDV é ortonormal em L_2[0,1] com a medida de Lebesgue, que
##        é a única propriedade da base que a Proposição 3 de E1.4 usa;
##   (A2) a reparametrização [Phi Q | Psi] da §5 de 01-identificabilidade.md:
##        Q'Q = I, as colunas de Phi Q têm integral zero e são ortogonais às
##        wavelets, e o bloco fica com 2^J - 1 colunas, como no periódico;
##   (A3) a cota pontual sum_{jk} psi^2_{jk}(u) <= C_psi 2^J, que E1.4 supõe
##        (Hipótese 3 de lá) e conferiu só na base periódica: aqui se mede
##        C_psi na CDV e se compara com o 1.35 da periódica;
##   (A4) o preço dos p q (2^{j0} - 1) coeficientes de escala não penalizados
##        no termo sigma^2 (não penalizados)/n do Teorema 1: contagem e
##        medição de E||P_A eps||_n^2 nos dois desenhos.
##
## Parte B (o caso eps > 0):
##   (B1) o cruzamento exato das duas exigências na família
##        eps = a (L-1) 2^{-J}: lambda_min(G_eps) > 0 pede a < 1/2, e o ganho
##        de taxa pede a >= 1; as duas faixas são disjuntas, e onde a Gram
##        sobrevive a queda por nível é a mesma 1/2 bit da periodização sem
##        margem nenhuma (ganho assintótico zero, D26);
##   (B2) a faixa admissível de J com eps fixo: o maior J com
##        lambda_min(G_eps) > 0, contra o limiar 2^J >= (L-1)/(2 eps).
##
## Parte C (a análise de s' dos cenários de wafc/R/dgp.R):
##   (C1) valor e derivada de cada componente nos dois extremos, que é o que
##        decide a ordem r da descontinuidade na emenda 0 = 1;
##   (C2) a queda por nível de ||theta_{j.}||_2 de cada componente nos três
##        regimes (periódico com eps = 0, periódico com margem e extensão,
##        intervalo), contra a tabela do .tex;
##   (C3) o s' de cada cenário, que é o mínimo sobre as suas componentes:
##        3/2 e 1/2 no regime da teoria (D27).

suppressPackageStartupMessages(library(WaveBased))

fam  <- "Daublets"
fs   <- 8L                  # filter.size = 8: N = 4 momentos nulos, L - 1 = 7
Nvm  <- fs / 2
Lm1  <- fs - 1
j0i  <- 4L                  # o wbasis() exige j0 >= 4 na CDV com fs = 8
wtab <- wtable(family = fam, filter.size = fs)

slope <- function(e) -diff(log2(e))
ok_all <- TRUE
chk <- function(cond, msg){
  ok_all <<- ok_all && isTRUE(cond)
  cat(sprintf("    %-72s %s\n", msg, if (isTRUE(cond)) "ok" else "FALHA"))
  invisible(cond)
}

## ===========================================================================
## Parte A. A rota do intervalo
## ===========================================================================
cat("=== Parte A: a base do intervalo (Cohen-Daubechies-Vial) ===\n")

ngA <- 2^14
uA  <- (seq_len(ngA) - 0.5) / ngA

## (A1) ortonormalidade em L_2[0,1] ----------------------------------------
cat("\n(A1) ortonormalidade da base CDV em L_2[0,1] (Lebesgue), j0 = 4:\n")
devA1 <- sapply(4:8, function(J){
  W <- wbasis(uA, j0 = j0i, J = J, family = fam, filter.size = fs,
              boundary = "interval")
  G <- crossprod(W) / ngA
  max(abs(G - diag(ncol(G))))
})
names(devA1) <- paste0("J=", 4:8)
print(signif(devA1, 2))
chk(all(devA1 < 5e-3), "max |G - I| abaixo do erro de quadratura da grade 2^14")

## (A2) a reparametrização [Phi Q | Psi] ------------------------------------
cat("\n(A2) reparametrizacao [Phi Q | Psi] (01-identificabilidade.md, §5):\n")
JA2 <- 6L
W   <- wbasis(uA, j0 = j0i, J = JA2, family = fam, filter.size = fs,
              boundary = "interval")
nsc <- 2^j0i                                  # colunas de escala de V_{j0}
Phi <- W[, 1:nsc, drop = FALSE]
Psi <- W[, -(1:nsc), drop = FALSE]
mu  <- colMeans(Phi)                          # mu_k = int_0^1 phi^int_{j0 k}
Q   <- qr.Q(qr(cbind(mu, diag(nsc))))[, -1, drop = FALSE]  # base de {mu'a = 0}
PhiQ <- Phi %*% Q
cat(sprintf("    1 = sum_k mu_k phi_k? desvio maximo na grade: %.1e\n",
            max(abs(drop(Phi %*% mu) - 1))))
cat(sprintf("    ||mu||^2 = %.6f (deve ser 1, pois ||1||_{L_2[0,1]} = 1)\n",
            sum(mu^2)))
cat(sprintf("    colunas: escala reparametrizada %d + wavelets %d = %d = 2^J - 1 = %d\n",
            ncol(PhiQ), ncol(Psi), ncol(PhiQ) + ncol(Psi), 2^JA2 - 1))
gPhiQ <- crossprod(PhiQ) / ngA
chk(max(abs(drop(Phi %*% mu) - 1)) < 5e-3, "1 pertence a V_{j0} (reproducao da constante)")
chk(abs(sum(mu^2) - 1) < 5e-3, "mu' mu = 1")
chk(max(abs(drop(mu %*% Q))) < 1e-10, "mu' Q = 0: as colunas Phi Q tem integral zero")
chk(max(abs(gPhiQ - diag(ncol(PhiQ)))) < 5e-3, "(Phi Q)'(Phi Q) = I: ortonormais em L_2[0,1]")
chk(max(abs(crossprod(PhiQ, Psi) / ngA)) < 5e-3, "Phi Q ortogonal as wavelets")
chk(ncol(PhiQ) + ncol(Psi) == 2^JA2 - 1, "o bloco fica com 2^J - 1 colunas, como no periodico")
chk(max(abs(colMeans(PhiQ))) < 5e-3 && max(abs(colMeans(Psi))) < 5e-3,
    "toda coluna do bloco tem integral zero (parte (i) de E1.4 transfere)")

## (A3) a cota pontual sum psi^2 <= C_psi 2^J -------------------------------
cat("\n(A3) cota pontual: sup_u sum_col W(u)^2 / 2^J, nas duas bases:\n")
Cpsi <- t(sapply(4:9, function(J){
  Wi  <- wbasis(uA, j0 = j0i, J = J, family = fam, filter.size = fs,
                boundary = "interval")
  mui <- colMeans(Wi[, 1:nsc, drop = FALSE])
  Qi  <- qr.Q(qr(cbind(mui, diag(nsc))))[, -1, drop = FALSE]
  Wr  <- cbind(Wi[, 1:nsc, drop = FALSE] %*% Qi, Wi[, -(1:nsc), drop = FALSE])
  Wp  <- wbasis(uA, j0 = 0L, J = J, family = fam, filter.size = fs,
                wavelet.table = wtab)[, -1, drop = FALSE]
  c(cdv_VJ = max(rowSums(Wi^2)) / 2^J,
    cdv_rep = max(rowSums(Wr^2)) / 2^J,
    periodica = max(rowSums(Wp^2)) / 2^J)
}))
rownames(Cpsi) <- paste0("J=", 4:9)
print(round(Cpsi, 3))
cat(sprintf("    C_psi(CDV)/C_psi(periodica) = %.1f\n",
            max(Cpsi[, "cdv_rep"]) / max(Cpsi[, "periodica"])))
chk(max(Cpsi[, "cdv_rep"]) < 10 &&
    max(Cpsi[, "cdv_rep"]) / min(Cpsi[, "cdv_rep"]) < 1.5,
    "C_psi da CDV e limitado e estavel em J (a Hipotese 3 de E1.4 vale na CDV)")
chk(max(Cpsi[, "cdv_rep"]) > 3 * max(Cpsi[, "periodica"]),
    "C_psi da CDV e varias vezes o da periodica (as funcoes de borda sao mais altas)")

## (A4) o preco dos coeficientes de escala nao penalizados ------------------
cat("\n(A4) coeficientes nao penalizados e o termo sigma^2 (nao pen.)/n:\n")
pA <- 3L; qA <- 2L; sigA <- 1
npen_per <- pA                                   # periodica, j0 = 0 (D2)
npen_cdv <- pA * (1L + qA * (2L^j0i - 1L))       # CDV: p + p q (2^{j0} - 1)
cat(sprintf("    p = %d, q = %d, j0 = %d: nao penalizados %d (periodica) e %d (CDV)\n",
            pA, qA, j0i, npen_per, npen_cdv))
set.seed(20260919)
nA  <- 1000L
nrep <- 200L
uS  <- matrix(runif(nA * qA), nA, qA)
xS  <- cbind(1, matrix(rnorm(nA * (pA - 1L)), nA, pA - 1L))
JS  <- 6L
Bcdv <- do.call(cbind, lapply(seq_len(pA), function(l){
  do.call(cbind, lapply(seq_len(qA), function(m){
    Wi  <- wbasis(uS[, m], j0 = j0i, J = JS, family = fam, filter.size = fs,
                  boundary = "interval")[, 1:nsc, drop = FALSE]
    Wg  <- wbasis(uA,      j0 = j0i, J = JS, family = fam, filter.size = fs,
                  boundary = "interval")[, 1:nsc, drop = FALSE]
    Qm  <- qr.Q(qr(cbind(colMeans(Wg), diag(nsc))))[, -1, drop = FALSE]
    xS[, l] * (Wi %*% Qm)
  }))
}))
Aper <- xS
Acdv <- cbind(xS, Bcdv)
cat(sprintf("    colunas de A: %d (periodica) e %d (CDV); posto %d e %d\n",
            ncol(Aper), ncol(Acdv), qr(Aper)$rank, qr(Acdv)$rank))
proj_sq <- function(A){
  qrA <- qr(A)
  mean(replicate(nrep, { e <- rnorm(nA, sd = sigA); mean(qr.fitted(qrA, e)^2) }))
}
mper <- proj_sq(Aper); mcdv <- proj_sq(Acdv)
cat(sprintf("    E||P_A eps||_n^2 medido: %.5f (periodica) e %.5f (CDV)\n", mper, mcdv))
cat(sprintf("    sigma^2 (nao pen.)/n previsto: %.5f e %.5f\n",
            sigA^2 * npen_per / nA, sigA^2 * npen_cdv / nA))
chk(abs(mper / (sigA^2 * npen_per / nA) - 1) < 0.15 &&
    abs(mcdv / (sigA^2 * npen_cdv / nA) - 1) < 0.15,
    "o termo e sigma^2 (nao penalizados)/n nos dois desenhos")
chk(qr(Acdv)$rank == npen_cdv, "A da CDV tem posto cheio (o Lema 4 de E1.5 exige)")

## ===========================================================================
## Parte B. O caso eps > 0: o cruzamento exato
## ===========================================================================
cat("\n=== Parte B: a margem eps ===\n")

ngB <- 2^14
uB  <- (seq_len(ngB) - 0.5) / ngB
lmin_G <- function(J, eps){
  W   <- wbasis(uB, j0 = 0L, J = J, family = fam, filter.size = fs,
                wavelet.table = wtab)          # constante + as 2^J - 1 wavelets
  idx <- which(uB >= eps & uB <= 1 - eps)
  min(eigen(crossprod(W[idx, , drop = FALSE]) / ngB,
            symmetric = TRUE, only.values = TRUE)$values)
}
## Funcao-teste: C^inf no interior mas g(0) != g(1), que e o caso em que a
## periodizacao cobra o termo 2^{-J} da Proposicao 2. Nao se usa um polinomio
## de grau < N (o u - 1/2 de E1.3b), que os momentos nulos anulam no interior
## e cujo erro restrito satura no zero de maquina.
g_lin <- exp(uB) - mean(exp(uB))
Wc <- list()
best_restr <- function(J, eps){
  key <- as.character(J)
  if (is.null(Wc[[key]]))
    Wc[[key]] <<- wbasis(uB, j0 = 0L, J = J, family = fam, filter.size = fs,
                         wavelet.table = wtab)[, -1, drop = FALSE]
  idx <- which(uB >= eps & uB <= 1 - eps)
  r   <- g_lin[idx] - qr.fitted(qr(Wc[[key]][idx, , drop = FALSE]), g_lin[idx])
  sqrt(sum(r^2) / ngB)
}

## (B1) a família eps = a (L-1) 2^{-J} --------------------------------------
cat("\n(B1) familia eps = a (L-1) 2^{-J}, L - 1 = 7, J = 5..9:\n")
JB <- 5:9
aa <- c(0, 0.25, 0.5, 1)
LMB <- sapply(aa, function(a) sapply(JB, function(J) lmin_G(J, a * Lm1 * 2^(-J))))
BRB <- sapply(aa, function(a) sapply(JB, function(J) best_restr(J, a * Lm1 * 2^(-J))))
dimnames(LMB) <- dimnames(BRB) <- list(paste0("J=", JB), paste0("a=", aa))
cat("lambda_min(G_eps):\n");                      print(signif(LMB, 3))
cat("melhor aproximacao de exp(u) medida em [eps, 1-eps]:\n"); print(signif(BRB, 3))
SB <- apply(BRB, 2, slope); rownames(SB) <- paste0("J=", JB[-1])
cat("queda por nivel (bits):\n");                 print(round(SB, 2))
cat(sprintf("    media das quedas: %s\n",
            paste(sprintf("a=%.2f: %.2f", aa, colMeans(SB)), collapse = "; ")))
gram_ok <- apply(LMB, 2, function(z) all(z > 1e-8))
razao <- BRB[nrow(BRB), ] / BRB[nrow(BRB), 1]
cat(sprintf("    erro restrito em J = 9, contra o de eps = 0: %s\n",
            paste(sprintf("a=%.2f: %.0e (%.0e x)", aa, BRB[nrow(BRB), ], razao),
                  collapse = "; ")))
cat(sprintf("    Gram nao degenerada em a = %s\n",
            paste(aa[gram_ok], collapse = ", ")))
chk(!any(gram_ok & (colMeans(SB) > 1)),
    "nenhum a da ao mesmo tempo Gram viva e queda acima de 1/2 bit")
chk(all(abs(colMeans(SB)[gram_ok] - 0.5) < 0.15),
    "onde a Gram sobrevive a queda e 1/2 bit: o mesmo da periodizacao sem margem")
chk(gram_ok[["a=0.25"]] && !gram_ok[["a=0.5"]] && !gram_ok[["a=1"]],
    "a < 1/2 preserva a Gram e a >= 1/2 a degenera (limiar 2 eps < (L-1) 2^{-J})")
chk(all(razao[!gram_ok] < 1e-4) && razao[["a=0.25"]] > 1e-3,
    "onde a Gram degenera o erro restrito desaba ao piso: e o mesmo fenomeno")

## (B2) a faixa admissível de J com eps fixo --------------------------------
cat("\n(B2) eps fixo: maior J com lambda_min(G_eps) > 0, contra o limiar\n")
epsF <- c(0.02, 0.05, 0.10)
JB2  <- 2:9
LMF  <- sapply(epsF, function(e) sapply(JB2, function(J) lmin_G(J, e)))
dimnames(LMF) <- list(paste0("J=", JB2), paste0("eps=", epsF))
print(signif(LMF, 3))
Jobs <- sapply(seq_along(epsF), function(i){
  z <- which(LMF[, i] > 1e-8); if (length(z)) max(JB2[z]) else NA_integer_ })
Jsuf <- floor(log2(Lm1 / (2 * epsF)))            # 2^J < (L-1)/(2 eps)
Jbia <- ceiling(log2(Lm1 / epsF))                # 2^J >= (L-1)/eps
tabB2 <- rbind(`maior J com Gram nao degenerada (medido)` = Jobs,
               `limiar suficiente de degenerescencia` = Jsuf,
               `menor J que exclui a faixa contaminada` = Jbia)
colnames(tabB2) <- paste0("eps=", epsF)
print(tabB2)
chk(all(Jobs <= Jsuf) && all(Jobs >= Jsuf - 1),
    "o J medido fica no limiar ou um nivel abaixo dele")
chk(all(Jbia > Jobs), "a faixa da Gram e a da exclusao nao se encontram para nenhum eps")

## ===========================================================================
## Parte C. A análise de s' dos cenários de wafc/R/dgp.R
## ===========================================================================
cat("\n=== Parte C: a regularidade efetiva das componentes de dgp.R ===\n")

## As seis componentes, na forma crua de wafc/R/dgp.R (a centralização e a
## normalização de lá são afins e não mudam a ordem da emenda nem a taxa).
raw <- list(
  sine   = function(u) sin(2 * pi * u),
  cosine = function(u) cos(4 * pi * u),
  cubic  = function(u) u^3 - 1.4 * u^2 + 0.4 * u,
  bumps  = function(u){
    t <- c(0.1,0.13,0.15,0.23,0.25,0.40,0.44,0.65,0.76,0.78,0.81)
    h <- c(4,5,3,4,5,4.2,2.1,4.3,3.1,5.1,4.2)
    w <- c(0.005,0.005,0.006,0.01,0.01,0.03,0.01,0.01,0.005,0.008,0.005)
    out <- numeric(length(u))
    for (i in seq_along(t)) out <- out + h[i] * (1 + abs(u - t[i]) / w[i])^(-4)
    out },
  blocks = function(u){
    t <- c(0.1,0.13,0.15,0.23,0.25,0.40,0.44,0.65,0.76,0.78,0.81)
    h <- c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
    out <- numeric(length(u))
    for (i in seq_along(t)) out <- out + h[i] * (1 + sign(u - t[i])) / 2
    out },
  heavisine = function(u) 4 * sin(4 * pi * u) - sign(u - 0.3) - sign(0.72 - u)
)
comps <- names(raw)

## (C1) a emenda: valor e derivada nos dois extremos ------------------------
cat("\n(C1) a emenda 0 = 1: valor e derivada de cada componente nos extremos\n")
hh <- 1e-6
ugC <- (seq_len(2^13) - 0.5) / 2^13
C1 <- t(sapply(comps, function(nm){
  g  <- raw[[nm]]
  gv <- g(ugC)
  c(`g(0)` = g(0), `g(1)` = g(1), `salto` = g(1) - g(0),
    `g'(0)` = (g(hh) - g(0)) / hh, `g'(1)` = (g(1) - g(1 - hh)) / hh,
    `sup|g|` = max(abs(gv)), `sup|g'|` = max(abs(diff(gv))) * 2^13)
}))
C1 <- cbind(C1, `salto de g'` = C1[, "g'(1)"] - C1[, "g'(0)"])
print(signif(C1[, c("g(0)", "g(1)", "salto", "sup|g|")], 4))
print(signif(C1[, c("g'(0)", "g'(1)", "salto de g'", "sup|g'|")], 4))
## A ordem r da descontinuidade na emenda, medida contra a escala da propria
## funcao: o que importa e se o salto e visivel, nao se e exatamente nulo.
rel_v <- abs(C1[, "salto"]) / C1[, "sup|g|"]
rel_d <- abs(C1[, "salto de g'"]) / C1[, "sup|g'|"]
cat("\n    salto relativo de g e de g' (contra sup|g| e sup|g'|):\n")
print(signif(cbind(`g` = rel_v, `g'` = rel_d), 3))
r_seam <- ifelse(rel_v > 1e-3, 0, ifelse(rel_d > 1e-3, 1, Inf))
cat("    ordem r da descontinuidade na emenda (0 salto, 1 quina, Inf casa):\n")
print(r_seam)
chk(all(rel_v < 1e-3),
    "nenhuma das seis salta na emenda: as tres de Donoho-Johnstone somam zero")
chk(r_seam[["cubic"]] == 1, "so a cubica tem quina na emenda (g'(0) = 0.4, g'(1) = 0.6)")
chk(all(is.infinite(r_seam[c("sine","cosine","bumps","blocks","heavisine")])),
    "as outras cinco casam valor e derivada na emenda")

## (C2) queda por nivel nos tres regimes ------------------------------------
ngC  <- 2^15
uC   <- (seq_len(ngC) - 0.5) / ngC
jmax <- 12L
epsC <- 0.10                              # margem fixa do regime (b)

smoothstep <- function(t){
  a <- exp(-1 / pmax(t, 1e-300)); b <- exp(-1 / pmax(1 - t, 1e-300))
  ifelse(t <= 0, 0, ifelse(t >= 1, 1, a / (a + b)))
}
chi <- smoothstep(2 * uC / epsC - 1) * smoothstep(2 * (1 - uC) / epsC - 1)
## Regime (b): a moduladora vive em [eps, 1-eps] (wafc_rescale), a componente
## e lida nessa escala, e a extensao admissivel e a propria formula fora do
## suporte, cortada pelo perfil C^inf chi_eps (Passo 2 do Lema 10).
arg <- (uC - epsC) / (1 - 2 * epsC)

## Coeficientes de todos os niveis numa passagem por regime. A coluna RUIDO
## carrega a constante: como toda wavelet das duas bases tem integral zero
## por construcao, ||(<psi_{jk}, 1>)_k||_2 mede o erro de avaliacao da base
## no nivel j, e e o piso abaixo do qual nada pode ser lido.
lvlnorms <- function(fmat, boundary, j0){
  acc <- NULL; chunk <- 2^11
  for (st in seq(1, ngC, by = chunk)){
    idx <- st:min(st + chunk - 1, ngC)
    W <- if (boundary == "periodic")
      wbasis(uC[idx], j0 = j0, J = jmax + 1L, family = fam, filter.size = fs,
             wavelet.table = wtab)
    else
      wbasis(uC[idx], j0 = j0, J = jmax + 1L, family = fam, filter.size = fs,
             boundary = "interval")
    if (is.null(acc)) acc <- matrix(0, ncol(W), ncol(fmat) + 1L)
    acc <- acc + crossprod(W, cbind(fmat[idx, , drop = FALSE], 1))
  }
  acc <- acc / ngC
  lv  <- c(rep(-1L, 2^j0), rep(j0:jmax, 2^(j0:jmax)))
  out <- t(sapply(j0:jmax, function(j) sqrt(colSums(acc[lv == j, , drop = FALSE]^2))))
  rownames(out) <- paste0("j=", j0:jmax)
  colnames(out) <- c(colnames(fmat), "RUIDO")
  out
}

Fa <- sapply(comps, function(nm) raw[[nm]](uC))
Fb <- sapply(comps, function(nm) chi * raw[[nm]](arg))
cat(sprintf("\n(C2) ||theta_{j.}||_2 por nivel, grade 2^15, j <= %d, eps = %.2f no regime (b)\n",
            jmax, epsC))
Lper <- lvlnorms(Fa, "periodic", 0L)
Lmar <- lvlnorms(Fb, "periodic", 0L)
Lint <- lvlnorms(Fa, "interval", j0i)
regs <- list(`(a) periodico eps=0` = Lper, `(b) margem + extensao` = Lmar,
             `(c) intervalo (CDV)` = Lint)
for (tag in names(regs)){
  cat(sprintf("\n  %s:\n", tag)); print(signif(regs[[tag]], 3))
}

## A janela de medicao comeca no primeiro nivel em que 2^{-j} fica abaixo da
## menor escala presente na componente: o periodo (sine, cosine), a escala em
## que a singularidade passa a dominar o termo suave (a cubica), a largura
## minima das protuberancias (0.005 em bumps, isto e 2^{-7.6}) e a menor
## distancia entre saltos (0.03 em blocks, 0.42 em heavisine). No regime (b)
## entra mais uma escala, a largura da transicao do corte chi_eps, que e
## eps/2 = 0.05 ~ 2^{-4.3} e pede tres niveis de folga. A janela termina no
## ultimo nivel acima de dez vezes o piso de avaliacao da base (coluna RUIDO).
jstart   <- c(sine = 3, cosine = 3, cubic = 5, bumps = 9, blocks = 4, heavisine = 4)
jstart_b <- pmax(jstart, 8)
win_slope <- function(M, j0, jini){
  js <- as.integer(sub("j=", "", rownames(M)))
  sapply(comps, function(nm){
    use <- which(js >= max(jini[[nm]], j0) & M[, nm] > 10 * M[, "RUIDO"])
    if (length(use) < 2) return(NA_real_)
    use <- min(use):max(use)
    (log2(M[min(use), nm]) - log2(M[max(use), nm])) / (length(use) - 1)
  })
}
janela <- function(M, j0, jini){
  js <- as.integer(sub("j=", "", rownames(M)))
  sapply(comps, function(nm){
    use <- which(js >= max(jini[[nm]], j0) & M[, nm] > 10 * M[, "RUIDO"])
    if (length(use) < 2) "no piso" else sprintf("%d..%d", js[min(use)], js[max(use)])
  })
}
SLO <- rbind(`(a) periodico eps=0`   = win_slope(Lper, 0L, jstart),
             `(b) margem + extensao` = win_slope(Lmar, 0L, jstart_b),
             `(c) intervalo (CDV)`   = win_slope(Lint, j0i, jstart))
JAN <- rbind(`(a) periodico eps=0`   = janela(Lper, 0L, jstart),
             `(b) margem + extensao` = janela(Lmar, 0L, jstart_b),
             `(c) intervalo (CDV)`   = janela(Lint, j0i, jstart))
cat("\n  janela de niveis efetivamente usada:\n"); print(JAN)
cat("  queda media por nivel na janela:\n");       print(round(SLO, 2))

esperado <- rbind(`(a) periodico eps=0`   = c(sine = 4, cosine = 4, cubic = 1.5,
                                              bumps = 1.5, blocks = 0.5, heavisine = 0.5),
                  `(b) margem + extensao` = c(4, 4, 4, 1.5, 0.5, 0.5),
                  `(c) intervalo (CDV)`   = c(4, 4, 4, 1.5, 0.5, 0.5))
colnames(esperado) <- comps
cat("\n  s' declarado pela tabela do .tex:\n"); print(esperado)
cat("\n")
for (i in seq_len(nrow(SLO))) for (nm in comps){
  obs <- SLO[i, nm]; exp_ <- esperado[i, nm]; tag <- rownames(SLO)[i]
  if (is.na(obs)){
    ## unico caso previsto: a cubica na CDV, de grau 3 < N = 4, reproduzida
    ## exatamente pelas funcoes de escala de V_{j0}.
    chk(nm == "cubic" && tag == "(c) intervalo (CDV)",
        sprintf("%s, %s: coeficientes abaixo do piso (reproducao exata)", tag, nm))
  } else if (exp_ >= Nvm){
    chk(obs >= Nvm - 0.5,
        sprintf("%s, %s: queda %.2f >= N - 1/2 = %.1f", tag, nm, obs, Nvm - 0.5))
  } else if (nm == "bumps"){
    chk(obs >= exp_ - 0.25,
        sprintf("%s, %s: queda %.2f >= s' - 1/4 = %.2f", tag, nm, obs, exp_ - 0.25))
  } else {
    chk(abs(obs - exp_) < 0.15,
        sprintf("%s, %s: queda %.2f contra s' = %.2f", tag, nm, obs, exp_))
  }
}

## (C3) o s' de cada cenário ------------------------------------------------
cat("\n(C3) s' de cada cenario = minimo sobre as suas componentes (D27):\n")
cen <- list(smooth = c("sine", "cubic", "cosine"),
            inhomogeneous = c("bumps", "blocks", "heavisine"))
tabC3 <- sapply(names(cen), function(cn) apply(esperado[, cen[[cn]], drop = FALSE], 1, min))
print(tabC3)
chk(tabC3["(a) periodico eps=0", "smooth"] == 1.5 &&
    tabC3["(a) periodico eps=0", "inhomogeneous"] == 0.5,
    "D27 confirmada: 3/2 no suave e 1/2 no nao homogeneo, no regime da teoria")
chk(tabC3["(b) margem + extensao", "smooth"] == 4 &&
    tabC3["(c) intervalo (CDV)", "smooth"] == 4,
    "com margem ou na CDV o cenario suave leria 4, e nao 3/2")

cat("\n")
if (ok_all) cat("OK\n") else stop("E1.8: conferencia numerica FALHOU (ver linhas acima)")
