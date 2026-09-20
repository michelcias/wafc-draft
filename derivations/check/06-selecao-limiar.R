# E1.7c. Conferência numérica da seleção de estrutura por limiarização
# (derivations/06-selecao-limiar.tex). Roda ANTES da prova; imprime OK ou
# falha com stop(). Notação: docs/notacao.md (congelada em E1.1).
#
# O que confere, em cinco partes:
#
#   I.   (determinístico, sem simulação) o Lema 13: o sanduíche
#        {||g_lm|| > t + D} subset S_hat subset {||g_lm|| > t - D} para todo
#        t >= 0 e todo D que domine o erro máximo por bloco; a janela
#        D <= t < delta - D dando S_hat = S; e a exatidão da janela, isto é,
#        que em delta = 2D existe configuração em que nenhum t acerta.
#
#   II.  a ponte com o código: com rescale = FALSE (o regime eps = 0 da
#        teoria, D26) a base é ortonormal em L_2[0,1] e a estatística do
#        limiar ||g_hat_lm||_{L_2} é EXATAMENTE a norma euclidiana do bloco
#        de coeficientes que wafc_blocks() devolve. Avaliação exata da base
#        (use.table = "never"), pela exceção de D31.
#
#   III. o Corolário 8 nos cenários de wafc/R/dgp.R, em duas regras de
#        sintonia (a da teoria, com sigma conhecido, e cv.min, que é o
#        padrão de D20): o erro máximo por bloco Delta_n DECOMPOSTO em viés
#        do sieve e estimação, a separação efetiva min ||Pi_J g_lm||, a
#        curva de acerto de estrutura contra o limiar e a probabilidade de
#        que ALGUM limiar acerte.
#
#   IV.  a hipótese de separação, medida: um bloco ativo com amplitude
#        delta em {1, 1/2, 1/4, 1/8}, a n e J fixos num regime em que o
#        viés do sieve é pequeno, para ver a janela fechar quando delta
#        deixa de dominar Delta_n.
#
#   V.   o LASSO limiarizado contra o sparse group LASSO no mesmo desenho e
#        nas mesmas dobras: acerto exato de estrutura e RMSE fora da
#        amostra.
#
# Uso: Rscript derivations/check/06-selecao-limiar.R [n_rep]
# Dependências: WaveBased, glmnet, Matrix (wafc/R/load.R), sparsegl (parte V).
# Tempo nesta máquina, com o padrão n_rep = 20: cerca de 7 min.

suppressPackageStartupMessages(source("wafc/R/load.R"))

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1L) as.integer(args[[1L]]) else 20L
set.seed(20260920)

ok <- TRUE
chk <- function(cond, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "ok" else "FALHA", msg))
  if (!isTRUE(cond)) ok <<- FALSE
  invisible(cond)
}

## Base fixada e avaliada por tabela nas partes de simulação (D31); as
## partes I e II não usam tabela (I não usa base; II mede precisão fina).
FS <- 8L
TB <- WaveBased::wtable(filter.size = FS)

## s' declarado de cada cenário (D27, e a Parte C de 07-rota-intervalo.tex).
S_PRIME <- c(smooth = 3 / 2, inhomogeneous = 1 / 2, null = 3 / 2)

## ---------------------------------------------------------------------------
## Parte I. O Lema 13, que é determinístico
## ---------------------------------------------------------------------------

cat("\nI. Lema 13 (o sanduiche e a janela), deterministico\n")

s_hat <- function(nh, t) nh > t           # a regra do enunciado, com ">"

sandwich_fail <- 0L
window_fail <- 0L
for (it in seq_len(2000L)) {
  pq <- sample(2:12, 1L)
  nt <- numeric(pq)                       # ||g_lm||, com alguns zeros
  act <- sample(c(TRUE, FALSE), pq, replace = TRUE)
  if (!any(act)) act[[1L]] <- TRUE
  nt[act] <- runif(sum(act), 0.05, 3)
  sgn <- sample(c(-1, 1), pq, replace = TRUE)
  nh <- abs(nt + sgn * runif(pq, 0, 0.6))
  err <- abs(nh - nt)                     # |  ||g_hat|| - ||g||  | <= ||g_hat - g||
  D <- max(err)
  delta <- min(nt[act])
  for (t in c(0, sort(runif(20L, 0, 3.5)))) {
    Sh <- s_hat(nh, t)
    if (any((nt > t + D) & !Sh)) sandwich_fail <- sandwich_fail + 1L
    if (any(Sh & !(nt > t - D))) sandwich_fail <- sandwich_fail + 1L
    if (t >= D && t < delta - D && !identical(Sh, act)) {
      window_fail <- window_fail + 1L
    }
  }
}
chk(sandwich_fail == 0L,
    sprintf("sanduiche {||g|| > t+D} <= S_hat <= {||g|| > t-D} em 2000 configuracoes (%d falhas)",
            sandwich_fail))
chk(window_fail == 0L,
    sprintf("janela D <= t < delta - D da S_hat = S (%d falhas)", window_fail))

## Exatidão da janela: em delta = 2D ela é vazia, e há configuração em que
## nenhum t acerta -- um bloco ativo de norma 2D e um nulo, com o erro
## saturando D nos dois, nas direções que os confundem.
D0 <- 0.3
nt0 <- c(2 * D0, 0)
nh0 <- c(2 * D0 - D0, 0 + D0)
chk(all(abs(nh0 - nt0) <= D0 + 1e-12) && nh0[[1L]] <= nh0[[2L]],
    "em delta = 2D as duas estatisticas empatam: nenhum limiar separa")
worst <- vapply(c(0, sort(runif(200L, 0, 1))),
                function(t) identical(s_hat(nh0, t), c(TRUE, FALSE)), TRUE)
chk(!any(worst),
    "nenhum t acerta a estrutura nessa configuracao (a janela e exata)")

## ---------------------------------------------------------------------------
## Parte II. A estatística do limiar é a norma do bloco de coeficientes
## ---------------------------------------------------------------------------

cat("\nII. ||g_hat_lm||_{L_2[0,1]} = ||theta_hat_lm||_2 com rescale = FALSE\n")

d2 <- simulate_wafc(400L, p = 3L, q = 2L, scenario = "smooth", seed = 11L)
fit2 <- wafc(d2[["x"]], d2[["u"]], d2[["y"]], J = 4L, rescale = FALSE,
             use.table = "never", filter.size = FS)
s2 <- fit2[["lambda"]][[60L]]
ng <- 2L^14L
grid2 <- matrix(rep((seq_len(ng) - 0.5) / ng, 2L), ng, 2L)
fn2 <- wafc_functions(fit2, s = s2, grid = grid2)
bl2 <- wafc_blocks(fit2, s = s2)
num2 <- matrix(0, 3L, 2L)
for (l in seq_len(3L)) {
  for (m in seq_len(2L)) num2[l, m] <- sqrt(mean(fn2[["g"]][[l, m]]^2))
}
rel2 <- max(abs(num2 - bl2[["norm"]]) / pmax(bl2[["norm"]], 1e-12))
chk(rel2 < 1e-5,
    sprintf("norma numerica em L_2[0,1] = ||theta_hat||_2 nos 6 blocos (erro relativo maximo %.2e)",
            rel2))

## Com a margem do código (rescale = TRUE, eps = 0.05) a identidade vale na
## variável reescalada: ||theta_hat_lm||_2 é a norma de g_hat no intervalo
## ALARGADO, de comprimento scale, e não no intervalo observado. A parcela
## que fica de fora é a energia da margem, e é o que separa a estatística do
## código da do enunciado. Mede-se aqui, porque é o que E2.5 vai usar.
fit2b <- wafc(d2[["x"]], d2[["u"]], d2[["y"]], J = 4L, use.table = "never",
              filter.size = FS)
s2b <- fit2b[["lambda"]][[60L]]
sc <- fit2b[["design"]][["scale"]]
fn2b <- wafc_functions(fit2b, s = s2b, n_grid = ng)   # grade = amplitude observada
bl2b <- wafc_blocks(fit2b, s = s2b)
raz <- matrix(0, 3L, 2L)
for (l in seq_len(3L)) {
  for (m in seq_len(2L)) {
    obs <- sqrt(mean(fn2b[["g"]][[l, m]]^2) * (1 - 2 * fit2b[["design"]][["eps"]][[m]]))
    raz[l, m] <- obs / max(bl2b[["norm"]][l, m], 1e-12)
  }
}
act2 <- d2[["structure"]] != ""
cat(sprintf("     com rescale = TRUE (eps = %.2f, scale = %.4f): razao entre a norma no\n",
            fit2b[["design"]][["eps"]][[1L]], sc[[1L]]))
cat(sprintf("     intervalo observado e ||theta_hat||_2 sqrt(scale) = %.4f a %.4f nos blocos ativos\n",
            min(raz[act2]), max(raz[act2])))
chk(all(raz[act2] > 0.8) && all(raz[act2] <= 1 + 1e-8),
    "a margem retem parte da energia, e a estatistica do codigo domina a do intervalo observado")

## ---------------------------------------------------------------------------
## Auxiliares das simulações
## ---------------------------------------------------------------------------

GRID_N <- 1024L
grid_mid <- (seq_len(GRID_N) - 0.5) / GRID_N

## Normas por bloco e erro por bloco, no regime eps = 0 da teoria, com o
## erro decomposto nas duas parcelas que compõem o D do Lema 13: o viés do
## sieve ||Pi_J g - g|| (Lema 1 de E1.3) e a estimação ||g_hat - Pi_J g||
## (Corolário 3 de E1.5). A projeção Pi_J g sai por quadratura na mesma
## grade, que é a mesma conta que check/02 usa.
##
## A estatística do limiar devolvida em 'norm' é a norma de g_hat calculada
## por quadratura, e não wafc_blocks()$norm; a Parte II é o que autoriza
## tratá-las como o mesmo número (elas coincidem a 1e-09 com rescale =
## FALSE). A Parte V usa wafc_blocks() diretamente, que é o caminho do
## código.
block_stats <- function(fit, s, dgp) {
  design <- fit[["design"]]
  p <- design[["p"]]
  q <- design[["q"]]
  gr <- matrix(rep(grid_mid, q), GRID_N, q)
  fn <- wafc_functions(fit, s = s, grid = gr)
  dg <- wafc_design(matrix(1, GRID_N, p), gr, spec = design)
  nh <- nt <- npj <- er <- bias <- est <- matrix(0, p, q)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[paste0(design[["xnames"]][[l]], ":",
                                        design[["unames"]][[m]])]]
      Psi <- as.matrix(dg[["Z"]][, idx, drop = FALSE])
      gh <- fn[["g"]][[l, m]]
      g0 <- dgp[["g"]][[l, m]]
      gt <- if (is.null(g0)) rep_len(0, GRID_N) else g0(grid_mid)
      gp <- as.numeric(Psi %*% (crossprod(Psi, gt) / GRID_N))
      nh[l, m] <- sqrt(mean(gh^2))
      nt[l, m] <- sqrt(mean(gt^2))
      npj[l, m] <- sqrt(mean(gp^2))
      er[l, m] <- sqrt(mean((gh - gt)^2))
      bias[l, m] <- sqrt(mean((gp - gt)^2))
      est[l, m] <- sqrt(mean((gh - gp)^2))
    }
  }
  list(norm = nh, truth = nt, proj = npj, err = er, bias = bias, est = est)
}

## Uma réplica. A regra "theory" usa o par (J_n, lambda_n) do Teorema 2 e do
## Corolário 2 com sigma conhecido; a regra "cv.min" usa a grade padrão de
## J de cv.wafc() e lambda.min (D20). As duas leem o mesmo dado.
one_rep <- function(n, scenario, seed, amplitude = 1, Jfix = NULL) {
  sp <- S_PRIME[[scenario]]
  d <- simulate_wafc(n, p = 3L, q = 2L, scenario = scenario, seed = seed,
                     amplitude = amplitude)
  Jth <- if (is.null(Jfix)) wafc_J_theory(n, s = sp) else Jfix
  des <- wafc_design(d[["x"]], d[["u"]], J = Jth, rescale = FALSE,
                     wavelet.table = TB, filter.size = FS)
  lam_th <- wafc_lambda_theory(des, sigma = d[["sigma"]])
  fit_th <- wafc(design = des, y = d[["y"]],
                 lambda = wafc_path_to(lam_th, des, d[["y"]]))
  Jgrid <- if (is.null(Jfix)) NULL else Jfix
  cv <- cv.wafc(d[["x"]], d[["u"]], d[["y"]], J = Jgrid, nfolds = 5L,
                rescale = FALSE, wavelet.table = TB, filter.size = FS)
  list(J.theory = Jth, J.cv = cv[["J.min"]],
       theory = block_stats(fit_th, lam_th, d),
       cv.min = block_stats(cv[["wafc.fit"]], cv[["lambda.min"]], d))
}

T_GRID <- exp(seq(log(0.002), log(2.5), length.out = 80L))

summarise <- function(stats) {
  R <- length(stats)
  act1 <- stats[[1L]][["truth"]] > 0
  Delta <- vapply(stats, function(z) max(z[["err"]]), 0)
  Dbias <- vapply(stats, function(z) max(z[["bias"]]), 0)
  Dest <- vapply(stats, function(z) max(z[["est"]]), 0)
  delta <- vapply(stats, function(z) min(z[["truth"]][z[["truth"]] > 0]), 0)
  deltaJ <- vapply(stats, function(z) min(z[["proj"]][z[["truth"]] > 0]), 0)
  sep <- vapply(stats, function(z) {
    a <- z[["truth"]] > 0
    (if (any(!a)) max(z[["norm"]][!a]) else -Inf) <
      (if (any(a)) min(z[["norm"]][a]) else Inf)
  }, TRUE)
  hit <- vapply(stats, function(z) {
    a <- as.vector(z[["truth"]] > 0)
    vapply(T_GRID, function(t) identical(as.vector(z[["norm"]] > t), a), TRUE)
  }, logical(length(T_GRID)))
  fp <- vapply(stats, function(z) {
    a <- z[["truth"]] > 0
    vapply(T_GRID, function(t) sum((z[["norm"]] > t) & !a), 0)
  }, numeric(length(T_GRID)))
  fn_ <- vapply(stats, function(z) {
    a <- z[["truth"]] > 0
    vapply(T_GRID, function(t) sum(!(z[["norm"]] > t) & a), 0)
  }, numeric(length(T_GRID)))
  curve <- rowMeans(matrix(hit, length(T_GRID), R))
  list(Delta = mean(Delta), bias = mean(Dbias), est = mean(Dest),
       delta = mean(delta), deltaJ = mean(deltaJ), sep = mean(sep),
       curve = curve, fp = rowMeans(matrix(fp, length(T_GRID), R)),
       fn = rowMeans(matrix(fn_, length(T_GRID), R)),
       best = T_GRID[[which.max(curve)]], best_rate = max(curve),
       active = act1)
}

## ---------------------------------------------------------------------------
## Parte III. O Corolário 8 nos cenários de dgp.R
## ---------------------------------------------------------------------------

cat(sprintf("\nIII. Corolario 8 nos cenarios (%d replicas por celula)\n", n_rep))

NS <- c(250L, 500L, 1000L, 2000L)
SC <- c("smooth", "inhomogeneous")
tab3 <- NULL
curves <- list()
for (sc_ in SC) {
  for (n in NS) {
    reps <- lapply(seq_len(n_rep), function(r)
      one_rep(n, sc_, seed = 1000L * match(sc_, SC) + 10L * round(log2(n)) + r))
    for (rule in c("theory", "cv.min")) {
      z <- summarise(lapply(reps, function(w) w[[rule]]))
      key <- paste(sc_, n, rule, sep = "/")
      curves[[key]] <- z
      sp <- S_PRIME[[sc_]]
      tab3 <- rbind(tab3, data.frame(
        scenario = sc_, n = n, rule = rule,
        J = if (rule == "theory") reps[[1L]][["J.theory"]]
            else mean(vapply(reps, function(w) w[["J.cv"]], 0)),
        Delta = z[["Delta"]], vies = z[["bias"]], estim = z[["est"]],
        delta = z[["delta"]], deltaJ = z[["deltaJ"]],
        razao = (log(n) / n)^(-sp / (2 * sp + 1)) * z[["est"]],
        sep = z[["sep"]], best_t = z[["best"]], acerto = z[["best_rate"]],
        stringsAsFactors = FALSE))
    }
  }
}
print(format(tab3, digits = 3))

for (sc_ in SC) {
  sp <- S_PRIME[[sc_]]
  for (rule in c("theory", "cv.min")) {
    sub <- tab3[tab3[["scenario"]] == sc_ & tab3[["rule"]] == rule, ]
    target <- (log(sub[["n"]]) / sub[["n"]])^(sp / (2 * sp + 1))
    slope <- unname(coef(lm(log(sub[["estim"]]) ~ log(target)))[[2L]])
    cat(sprintf("     %s/%s: parte de estimacao de %.4f a %.4f, inclinacao contra o alvo = %.2f\n",
                sc_, rule, sub[["estim"]][[1L]], sub[["estim"]][[nrow(sub)]],
                slope))
    chk(sub[["estim"]][[nrow(sub)]] < sub[["estim"]][[1L]],
        sprintf("%s/%s: a parte de estimacao de Delta_n decresce de n = %d a n = %d",
                sc_, rule, sub[["n"]][[1L]], sub[["n"]][[nrow(sub)]]))
    chk(slope > 0,
        sprintf("%s/%s: e acompanha a taxa do Corolario 4 (inclinacao > 0)",
                sc_, rule))
  }
  sub <- tab3[tab3[["scenario"]] == sc_ & tab3[["rule"]] == "cv.min", ]
  chk(sub[["acerto"]][[nrow(sub)]] >= sub[["acerto"]][[1L]] - 1e-12,
      sprintf("%s/cv.min: o acerto no melhor limiar nao piora com n (%.2f -> %.2f)",
              sc_, sub[["acerto"]][[1L]], sub[["acerto"]][[nrow(sub)]]))
  chk(sub[["sep"]][[nrow(sub)]] >= 0.8,
      sprintf("%s/cv.min: em n = %d algum limiar acerta em %.0f%% das replicas",
              sc_, sub[["n"]][[nrow(sub)]], 100 * sub[["sep"]][[nrow(sub)]]))
}

cat("\n     curva de acerto contra o limiar (cv.min, n = 2000)\n")
idx <- round(seq(6L, 72L, length.out = 9L))
for (sc_ in SC) {
  z <- curves[[paste(sc_, 2000L, "cv.min", sep = "/")]]
  cat(sprintf("     %-14s t     = %s\n", sc_,
              paste(sprintf("%6.3f", T_GRID[idx]), collapse = " ")))
  cat(sprintf("     %-14s P(=S) = %s\n", "",
              paste(sprintf("%6.2f", z[["curve"]][idx]), collapse = " ")))
  cat(sprintf("     %-14s falso+= %s\n", "",
              paste(sprintf("%6.2f", z[["fp"]][idx]), collapse = " ")))
  cat(sprintf("     %-14s falso-= %s\n", "",
              paste(sprintf("%6.2f", z[["fn"]][idx]), collapse = " ")))
}

## ---------------------------------------------------------------------------
## Parte IV. A hipótese de separação, medida
## ---------------------------------------------------------------------------

cat("\nIV. A hipotese de separacao: a janela fecha quando delta deixa de dominar\n")

DELTAS <- c(1, 1 / 2, 1 / 4, 1 / 8)
J4 <- 4L
tab4 <- NULL
for (dl in DELTAS) {
  amp <- matrix(1, 3L, 2L)
  amp[2L, 1L] <- dl                       # o bloco (2,1) do cenario "smooth"
  reps <- lapply(seq_len(n_rep), function(r)
    one_rep(1000L, "smooth", seed = 5000L + 100L * round(8 * dl) + r,
            amplitude = as.vector(amp), Jfix = J4))
  z <- summarise(lapply(reps, function(w) w[["cv.min"]]))
  tab4 <- rbind(tab4, data.frame(
    delta = dl, Delta = z[["Delta"]], estim = z[["est"]],
    razao = dl / z[["Delta"]], janela = as.integer(dl > 2 * z[["Delta"]]),
    sep = z[["sep"]], best_t = z[["best"]], acerto = z[["best_rate"]],
    stringsAsFactors = FALSE))
}
print(format(tab4, digits = 3))
chk(tab4[["acerto"]][[1L]] > tab4[["acerto"]][[nrow(tab4)]],
    sprintf("o acerto cai de %.2f (delta = 1) a %.2f (delta = 1/8)",
            tab4[["acerto"]][[1L]], tab4[["acerto"]][[nrow(tab4)]]))
chk(all(diff(tab4[["razao"]]) < 0),
    "a razao delta/Delta_n decresce monotonamente com delta")

## ---------------------------------------------------------------------------
## Parte V. LASSO limiarizado contra o sparse group LASSO
## ---------------------------------------------------------------------------

cat("\nV. LASSO limiarizado contra o sparse group LASSO, no mesmo desenho\n")

if (!requireNamespace("sparsegl", quietly = TRUE)) {
  cat("     sparsegl nao instalado; parte V pulada\n")
} else {
  n5 <- 1000L
  nrep5 <- max(10L, n_rep %/% 2L)
  tab5 <- NULL
  for (sc_ in SC) {
    J <- 4L
    acc_l <- acc_r <- acc_s <- rmse_l <- rmse_s <- numeric(nrep5)
    for (r in seq_len(nrep5)) {
      d <- simulate_wafc(n5, p = 3L, q = 2L, scenario = sc_,
                         seed = 9000L + 100L * match(sc_, SC) + r)
      dte <- simulate_wafc(1000L, p = 3L, q = 2L, scenario = sc_,
                           seed = 90000L + 100L * match(sc_, SC) + r)
      fid <- sample(rep_len(seq_len(10L), n5))
      act <- as.vector(d[["structure"]] != "")
      cvl <- cv.wafc(d[["x"]], d[["u"]], d[["y"]], J = J, foldid = fid,
                     penalty = "lasso", rescale = FALSE, wavelet.table = TB,
                     filter.size = FS)
      cvs <- cv.wafc(d[["x"]], d[["u"]], d[["y"]], J = J, foldid = fid,
                     penalty = "sglasso", rescale = FALSE, wavelet.table = TB,
                     filter.size = FS)
      bl <- wafc_blocks(cvl[["wafc.fit"]], s = cvl[["lambda.min"]])
      bs <- wafc_blocks(cvs[["wafc.fit"]], s = cvs[["lambda.min"]])
      ## dois limiares: o melhor da grade, que é o teto do que a
      ## limiarização pode dar, e a regra invariante de escala
      ## t = 0.15 max_lm ||g_hat_lm||, que é escolha empírica e não teoria
      ## (a calibração é de E2.4/E2.5).
      acc_l[[r]] <- as.numeric(any(vapply(T_GRID, function(t)
        identical(as.vector(bl[["norm"]] > t), act), TRUE)))
      acc_r[[r]] <- as.numeric(identical(
        as.vector(bl[["norm"]] > 0.15 * max(bl[["norm"]])), act))
      acc_s[[r]] <- as.numeric(identical(as.vector(bs[["nonzero"]] > 0L), act))
      rmse_l[[r]] <- sqrt(mean((dte[["y"]] -
        predict(cvl[["wafc.fit"]], dte[["x"]], dte[["u"]],
                s = cvl[["lambda.min"]]))^2))
      rmse_s[[r]] <- sqrt(mean((dte[["y"]] -
        predict(cvs[["wafc.fit"]], dte[["x"]], dte[["u"]],
                s = cvs[["lambda.min"]]))^2))
    }
    tab5 <- rbind(tab5, data.frame(
      scenario = sc_, n = n5, J = J, replicas = nrep5,
      lasso_melhor_t = mean(acc_l), lasso_regra_015 = mean(acc_r),
      sglasso = mean(acc_s), rmse_lasso = mean(rmse_l),
      rmse_sglasso = mean(rmse_s), sigma = round(d[["sigma"]], 3),
      stringsAsFactors = FALSE))
  }
  print(format(tab5, digits = 3))
  chk(all(tab5[["lasso_melhor_t"]] >= tab5[["sglasso"]]),
      "o LASSO limiarizado nao fica atras do sparse group LASSO em estrutura")
  chk(all(tab5[["rmse_lasso"]] <= 1.02 * tab5[["rmse_sglasso"]]),
      "e nao paga predicao por isso (RMSE dentro de 2%)")
}

## ---------------------------------------------------------------------------

cat("\n")
if (ok) cat("OK\n") else stop("Alguma verificacao falhou.", call. = FALSE)
