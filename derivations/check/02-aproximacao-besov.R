## E1.3 e E1.3b -- conferência numérica dos lemas de aproximação em Besov
## (derivations/02-aproximacao-besov.tex).
##
## Mede o erro de projeção || g - Pi_J g || em L_2[0,1] e em L_inf, para a
## base periodizada com j0 = 0 e a constante descartada (W_J = V_J menos as
## constantes, dimensão 2^J - 1), e compara a queda por nível com 2^{-J s'},
## s' = s - (1/pi - 1/2)_+. As bases são avaliadas pelo WaveBased (wbasis(),
## wtable()), nunca reimplementadas. Roda com
##   Rscript derivations/check/02-aproximacao-besov.R
## e termina imprimindo OK. Leva cerca de 2 min (partes B e C).
##
## O que se confere:
##   (1) Parseval: o erro L_2 da projeção coincide com a soma de cauda dos
##       coeficientes ao quadrado, que é o que a prova usa;
##   (2) sin(2 pi u), periódica e C^inf: queda por nível >= N (momentos
##       nulos do filtro), em L_2 e em L_inf;
##   (3) bumps (Donoho & Johnstone), com bumps(0) ~ bumps(1) ~ 0: as normas
##       por nível ||theta_{j.}||_2 caem como 2^{-j s'} com s' = 3/2 (a função
##       tem quinas |u - t_i|, logo está em B^s_{2,inf} só para s <= 3/2); a
##       fase da quina em relação à grade diádica faz a queda oscilar de nível
##       para nível, e o que é estável é a média sobre quatro níveis. A
##       largura mínima das protuberâncias é 0.005 ~ 2^{-7.6}, e por isso a
##       taxa só aparece de j ~ 9 em diante: uma quina periódica isolada,
##       |sin(pi (u - 0.4))|, serve de referência exata da taxa 3/2;
##   (4) função não periódica (g(0) != g(1)) na base periodizada: a queda em
##       L_2 trava em 1/2, qualquer que seja a regularidade interior, e o
##       erro L_inf não vai a zero (o termo de borda da Proposição do .tex);
##   (5) a mesma função na base do intervalo (Cohen-Daubechies-Vial) volta a
##       ter a taxa cheia (a linear é reproduzida exatamente).
## A parte C acrescenta, para E1.3b, (6) a taxa restrita ao suporte, (7)
## lambda_min da Gram restrita e (8) a degradação de C(eps); o cabeçalho da
## parte diz o que cada uma mede.
##
## Grades: 2^12 pontos (projeção por mínimos quadrados, J <= 9), 2^15 pontos
## (coeficientes por quadratura, um nível de cada vez, j <= 12) e 2^14 pontos
## (partes C7 e C8). É mais que o "J <= 3" da regra geral de conferência
## porque medir uma taxa exige vários níveis assintóticos; o custo é de
## segundos.

suppressPackageStartupMessages(library(WaveBased))

fam  <- "Daublets"
fs   <- 8L                 # filter.size = 8: N = 4 momentos nulos
Nvm  <- fs / 2
wtab <- wtable(family = fam, filter.size = fs)

## --- funções-teste ---------------------------------------------------------
bumps <- function(t){
  pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
  hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
  wth <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
  out <- numeric(length(t))
  for (i in seq_along(pos))
    out <- out + hgt[i] * (1 + abs((t - pos[i]) / wth[i]))^(-4)
  out
}
centre <- function(g) g - mean(g)   # E[g(U)] = 0 com U uniforme na grade
slope  <- function(e) -diff(log2(e))  # queda por nível, em bits

## ===========================================================================
## Parte A: projeção por mínimos quadrados em W_J, grade de 2^12, J <= 9
## ===========================================================================
ngrid <- 2^12
u     <- (seq_len(ngrid) - 0.5) / ngrid
Jmax  <- 9L
Js    <- 1:Jmax

funs <- list(
  sin    = centre(sin(2 * pi * u)),
  bumps  = centre(bumps(u)),
  linear = centre(u),                 # g(0) != g(1): não periódica
  expo   = centre(exp(u))             # idem, C^inf no interior
)

proj_err <- function(g, J, boundary = "periodic", j0 = 0L){
  W <- wbasis(u, j0 = j0, J = J, family = fam, filter.size = fs,
              boundary = boundary,
              wavelet.table = if (boundary == "periodic") wtab else NULL)
  if (boundary == "periodic" && j0 == 0L){
    stopifnot(all(abs(W[, 1] - 1) < 1e-10))   # a primeira coluna é phi_00 = 1
    W <- W[, -1, drop = FALSE]                 # W_J: as 2^J - 1 wavelets
  }
  r <- g - W %*% qr.coef(qr(W), g)
  c(L2 = sqrt(mean(r^2)), Linf = max(abs(r)))
}

## (1) Parseval: erro de projeção = soma de cauda dos coeficientes
Wfull <- wbasis(u, j0 = 0L, J = Jmax, family = fam, filter.size = fs,
                wavelet.table = wtab)[, -1]
G     <- crossprod(Wfull) / ngrid
cat(sprintf("ortonormalidade na grade de 2^12 (max |G - I|): %.1e\n",
            max(abs(G - diag(ncol(G))))))
lev   <- rep(0:(Jmax - 1), 2^(0:(Jmax - 1)))
g     <- funs[["bumps"]]
theta <- drop(crossprod(Wfull, g)) / ngrid
lvl2  <- tapply(theta^2, lev, sum)
rem   <- proj_err(g, Jmax)[["L2"]]^2          # resto além de Jmax
tail_pred <- sapply(1:(Jmax - 1), function(J)
  sqrt(sum(lvl2[as.character(J:(Jmax - 1))]) + rem))
tail_obs  <- sapply(1:(Jmax - 1), function(J) proj_err(g, J)[["L2"]])
dev1 <- max(abs(tail_obs - tail_pred) / tail_obs)
ok1  <- dev1 < 1e-2
cat("\n(1) bumps: erro L_2 da projeção em W_J vs soma de cauda dos coeficientes\n")
print(round(rbind(J = 1:(Jmax - 1), projecao = tail_obs, cauda = tail_pred), 5))
cat(sprintf("    desvio relativo máximo %.1e -> %s\n", dev1, if (ok1) "ok" else "FALHA"))

## erros por nível para as quatro funções
err <- lapply(funs, function(g) t(sapply(Js, function(J) proj_err(g, J))))
tab  <- sapply(err, function(e) e[, "L2"]);   rownames(tab)  <- paste0("J=", Js)
tabi <- sapply(err, function(e) e[, "Linf"]); rownames(tabi) <- paste0("J=", Js)
sl   <- sapply(err, function(e) slope(e[, "L2"]));   rownames(sl)  <- paste0("J=", Js[-1])
sli  <- sapply(err, function(e) slope(e[, "Linf"])); rownames(sli) <- paste0("J=", Js[-1])
cat("\nerro L_2 de projeção em W_J (base periodizada):\n");  print(signif(tab, 3))
cat("queda por nível em L_2:\n");                              print(round(sl, 2))
cat("\nerro L_inf de projeção em W_J:\n");                     print(signif(tabi, 3))
cat("queda por nível em L_inf:\n");                            print(round(sli, 2))

## (2) sin: queda por nível >= N acima do piso numérico (J >= 4)
e_sin <- err[["sin"]][, "L2"]
use   <- which(e_sin[-1] > 1e-9); use <- use[use >= 3]     # quedas J-1 -> J, J >= 4
ok2   <- length(use) >= 3 && all(sl[use, "sin"] >= Nvm - 0.25) &&
         all(sli[use, "sin"] >= Nvm - 0.25)
cat(sprintf("\n(2) sin: quedas L_2 em J=%s: %s; L_inf: %s; esperado >= N - 1/4 = %.2f -> %s\n",
            paste(use + 1, collapse = ","),
            paste(round(sl[use, "sin"], 2), collapse = ", "),
            paste(round(sli[use, "sin"], 2), collapse = ", "),
            Nvm - 0.25, if (ok2) "ok" else "FALHA"))

## (4) não periódica na base periodizada: L_2 trava em 1/2; L_inf não decai
ok4 <- TRUE
for (nm in c("linear", "expo")){
  s_nm <- sl[(Jmax - 3):(Jmax - 1), nm]
  jump <- abs(funs[[nm]][ngrid] - funs[[nm]][1])
  okk  <- all(abs(s_nm - 0.5) < 0.1) && tabi[Jmax, nm] > 0.25 * jump
  ok4  <- ok4 && okk
  cat(sprintf("(4) %s (|g(1)-g(0)| = %.2f), periodizada: quedas L_2 em J=%d..%d: %s (esperado 1/2); erro L_inf em J=%d: %.2f -> %s\n",
              nm, jump, Jmax - 2, Jmax, paste(round(s_nm, 2), collapse = ", "),
              Jmax, tabi[Jmax, nm], if (okk) "ok" else "FALHA"))
}

## (5) base do intervalo (CDV): o wbasis() exige j0 >= 4 para filter.size = 8
j0i   <- 4L
Ji    <- j0i:Jmax
err_i <- sapply(c("linear", "expo"), function(nm)
  sapply(Ji, function(J) proj_err(funs[[nm]], J, boundary = "interval", j0 = j0i)[["L2"]]))
rownames(err_i) <- paste0("J=", Ji)
cat("\n(5) erro L_2 na base do intervalo (V_J inteiro, j0 = 4) vs periodizada:\n")
print(signif(cbind(err_i, tab[paste0("J=", Ji), c("linear", "expo")]), 3))
ok5 <- all(err_i[, "linear"] < 1e-8) &&
       all(err_i[, "expo"] < 1e-3 * tab[paste0("J=", Ji), "expo"])
cat(sprintf("    linear reproduzida (< 1e-8) e expo com erro < 1e-3 do periodizado -> %s\n",
            if (ok5) "ok" else "FALHA"))

## ===========================================================================
## Parte B: normas por nível por quadratura, grade de 2^15, níveis j <= 12
## ===========================================================================
ngridB <- 2^15
uB     <- (seq_len(ngridB) - 0.5) / ngridB
jB     <- 6:12
lvlnorm <- function(g, j, chunk = 2^11){
  acc <- numeric(2^j)
  for (st in seq(1, ngridB, by = chunk)){
    idx <- st:min(st + chunk - 1, ngridB)
    W   <- wbasis(uB[idx], j0 = j, J = j + 1, family = fam, filter.size = fs,
                  wavelet.table = wtab)[, -(1:2^j), drop = FALSE]   # só psi_{j.}
    acc <- acc + drop(crossprod(W, g[idx]))
  }
  sqrt(sum((acc / ngridB)^2))
}
funsB <- list(bumps = centre(bumps(uB)),
              kink  = centre(abs(sin(pi * (uB - 0.4)))))   # quina periódica isolada
lvlB  <- sapply(funsB, function(g) sapply(jB, lvlnorm, g = g))
rownames(lvlB) <- paste0("j=", jB)
slB   <- apply(lvlB, 2, slope); rownames(slB) <- paste0("j=", jB[-1])
cat("\n(3) normas por nível ||theta_{j.}||_2 (grade de 2^15):\n"); print(signif(lvlB, 3))
cat("queda por nível:\n"); print(round(slB, 2))
sp     <- 1.5
last   <- (length(jB) - 3):(length(jB) - 1)           # as três últimas quedas
m_b    <- mean(slB[last, "bumps"])
m_k    <- mean(slB[(length(jB) - 4):(length(jB) - 1), "kink"])  # quatro: uma fase
ok3    <- m_b >= sp - 0.25 && abs(m_k - sp) < 0.15
cat(sprintf("    bumps: queda média em j=%d..%d = %.2f (esperado ~ s' = %.1f, tolerância 1/4); quina periódica: média de quatro quedas = %.2f -> %s\n",
            jB[length(jB) - 3], jB[length(jB)], m_b, sp, m_k, if (ok3) "ok" else "FALHA"))
## erro de projeção por soma de cauda (resto além de 12 pela taxa 3/2)
tailB <- sapply(seq_along(jB), function(i)
  sqrt(sum(lvlB[i:length(jB), "bumps"]^2) + lvlB[length(jB), "bumps"]^2 / (2^(2 * sp) - 1)))
cat("    bumps: erro L_2 por soma de cauda, J = 6..12, e erro * 2^{J s'}:\n")
print(round(rbind(J = jB, erro = tailB, `erro*2^(J s')` = tailB * 2^(jB * sp)), 4))

## ===========================================================================
## Parte C (E1.3b): a margem eps, o lema de extensão e a Gram restrita
## ===========================================================================
## Sob D23 a moduladora não vive em [0,1] e sim no suporte [eps, 1-eps], e a
## base é periodizada no intervalo maior. O que se confere aqui:
##   (6) o viés que E1.5 consome é o inf sobre W_J da norma L_2(P_U), isto é a
##       melhor aproximação medida SÓ em [eps, 1-eps]; o Lema 10 diz que ela
##       volta à taxa cheia 2^{-J s'} quando eps > 0 é fixo, contra os 1/2 bit
##       por nível da Proposição do custo da periodização (que é o caso
##       eps = 0). Com eps acoplado a J a constante C(eps) degrada junto e o
##       ganho é parcial;
##   (7) lambda_min(G_eps), com G_eps a Gram da base restrita ao suporte, que
##       é a constante que substitui c_U em E1.4: eps > 0 fixo derruba-a a
##       zero assim que V_J contém uma função de escala que cabe na margem
##       periodizada, de comprimento 2 eps, isto é 2^J >= (L-1)/(2 eps);
##   (8) a degradação de C(eps) do Lema 10: a norma de sequência da extensão
##       explícita de u - 1/2 cresce como eps^{-(s - 1/pi)}.

eps_of <- list("0" = function(J) 0, "2^-(J+1)" = function(J) 2^(-J - 1),
               "1.9^-J" = function(J) 1.9^(-J),
               "0.05" = function(J) 0.05, "0.10" = function(J) 0.10)

## (6) melhor aproximação em W_J medida só no suporte -----------------------
JC <- 4:Jmax
Wc <- lapply(JC, function(J)
  wbasis(u, j0 = 0L, J = J, family = fam, filter.size = fs,
         wavelet.table = wtab)[, -1, drop = FALSE])
names(Wc) <- as.character(JC)
best_restr <- function(g, J, eps){
  W   <- Wc[[as.character(J)]]
  idx <- which(u >= eps & u <= 1 - eps)
  r   <- g[idx] - qr.fitted(qr(W[idx, , drop = FALSE]), g[idx])   # qr.fitted
  sqrt(sum(r^2) / ngrid)                     # || . ||_{L_2[eps, 1-eps]}
}
ok6 <- TRUE
for (nm in c("linear", "expo")){
  g <- funs[[nm]]
  M <- sapply(names(eps_of), function(en)
              sapply(JC, function(J) best_restr(g, J, eps_of[[en]](J))))
  rownames(M) <- paste0("J=", JC)
  S <- apply(M, 2, slope); rownames(S) <- paste0("J=", JC[-1])
  cat(sprintf("\n(6) %s: melhor aproximação em W_J medida em [eps, 1-eps]:\n", nm))
  print(signif(M, 3))
  cat("queda por nível:\n"); print(round(S, 2))
  m19 <- mean(S[, "1.9^-J"])
  okk <- all(abs(S[, "0"] - 0.5) < 0.05) &&                       # Proposição
         all(abs(S[(nrow(S) - 1):nrow(S), "2^-(J+1)"] - 0.5) < 0.1) &&
         m19 > 0.7 && m19 < 1.3 &&
         all(S[1:2, "0.05"] > 1.5) && all(S[1:2, "0.10"] > 1.5)
  ok6 <- ok6 && okk
  cat(sprintf("    eps = 0: 1/2 bit; eps = 2^-(J+1): 1/2 bit; eps = 1.9^-J: %.2f em média; margem fixa: acima de 3/2 até o piso numérico -> %s\n",
              m19, if (okk) "ok" else "FALHA"))
}

## (7) lambda_min da Gram restrita ao suporte -------------------------------
ngridG <- 2^14
uG     <- (seq_len(ngridG) - 0.5) / ngridG
JG     <- 2:8
specG  <- function(J, eps){
  W   <- wbasis(uG, j0 = 0L, J = J, family = fam, filter.size = fs,
                wavelet.table = wtab)          # constante + as 2^J - 1 wavelets
  idx <- which(uG >= eps & uG <= 1 - eps)
  ev  <- sort(eigen(crossprod(W[idx, , drop = FALSE]) / ngridG,
                    symmetric = TRUE, only.values = TRUE)$values)
  c(lmin = ev[1], nzero = sum(ev < 1e-8))
}
SP <- lapply(names(eps_of), function(en)
             sapply(JG, function(J) specG(J, eps_of[[en]](J))))
names(SP) <- names(eps_of)
LM <- sapply(SP, function(z) z["lmin", ]);  rownames(LM) <- paste0("J=", JG)
NZ <- sapply(SP, function(z) z["nzero", ]); rownames(NZ) <- paste0("J=", JG)
## limiar previsto: uma função de escala de V_J cabe na margem periodizada,
## que na circunferência é um único intervalo de comprimento 2 eps
Jstar <- sapply(names(eps_of), function(en){
  Jt <- JG[sapply(JG, function(J){ e <- eps_of[[en]](J)
        e > 0 && (fs - 1) / 2^J <= 2 * e })]
  if (length(Jt)) min(Jt) else NA_integer_
})
Jobs <- sapply(colnames(NZ), function(cn){
  Jt <- JG[NZ[, cn] > 0]; if (length(Jt)) min(Jt) else NA_integer_
})
cat("\n(7) lambda_min de G_eps (constante + wavelets, q = 1, Lebesgue no suporte):\n")
print(signif(LM, 3))
cat("direções nulas (autovalor < 1e-8):\n"); print(NZ)
cat("primeiro J com direção nula, observado e previsto por 2^J >= (L-1)/(2 eps):\n")
print(rbind(observado = Jobs, previsto = Jstar))
ok7 <- all(abs(LM[, "0"] - 1) < 1e-6) &&
       all(abs(LM[, "2^-(J+1)"] - LM[1, "2^-(J+1)"]) < 0.05) &&
       all(abs(LM[, "1.9^-J"] / LM[1, "1.9^-J"] - 1) < 0.35) &&
       is.na(Jobs[["0"]]) && is.na(Jobs[["2^-(J+1)"]]) && is.na(Jobs[["1.9^-J"]]) &&
       all(abs(Jobs[c("0.05", "0.10")] - Jstar[c("0.05", "0.10")]) <= 1)
cat(sprintf("    eps = 0 dá 1 (ortonormal); margem acoplada a J estabiliza em %.3f (2^-(J+1)) e %.4f (1.9^-J); margem fixa degenera a um nível do limiar -> %s\n",
            LM[1, "2^-(J+1)"], LM[1, "1.9^-J"], if (ok7) "ok" else "FALHA"))

## (8) como a constante C(eps) do Lema 10 degrada ---------------------------
ngridE <- 2^14
uE     <- (seq_len(ngridE) - 0.5) / ngridE
smoothstep <- function(t){                    # 0 em t <= 0, 1 em t >= 1, C^inf
  a <- exp(-1 / pmax(t, 1e-300)); b <- exp(-1 / pmax(1 - t, 1e-300))
  ifelse(t <= 0, 0, ifelse(t >= 1, 1, a / (a + b)))
}
glue <- function(g, eps)                      # g~ = chi_eps * g, emenda em 0 = 1
  smoothstep(2 * uE / eps - 1) * smoothstep(2 * (1 - uE) / eps - 1) * g
lvlnormE <- function(g, j, chunk = 2^11){
  acc <- numeric(2^j)
  for (st in seq(1, ngridE, by = chunk)){
    idx <- st:min(st + chunk - 1, ngridE)
    W   <- wbasis(uE[idx], j0 = j, J = j + 1, family = fam, filter.size = fs,
                  wavelet.table = wtab)[, -(1:2^j), drop = FALSE]
    acc <- acc + drop(crossprod(W, g[idx]))
  }
  sqrt(sum((acc / ngridE)^2))
}
epsE <- c(0.4, 0.2, 0.1, 0.05); jE <- 0:10; sE <- c(1, 1.5, 2)
g0   <- uE - 0.5                              # g(0) != g(1): o caso da Proposição
Bn   <- sapply(epsE, function(e){
  gt <- glue(g0, e); gt <- gt - mean(gt)
  ln <- sapply(jE, lvlnormE, g = gt)
  sapply(sE, function(s) max(2^(jE * s) * ln))   # || . ||_{b^s_{2,inf}}
})
dimnames(Bn) <- list(paste0("s=", sE), paste0("eps=", epsE))
expE <- t(apply(Bn, 1, function(b) -diff(log2(b)) / diff(log2(epsE))))
colnames(expE) <- paste0("eps=", epsE[-1])
cat("\n(8) norma de sequência da extensão explícita de u - 1/2:\n"); print(signif(Bn, 3))
cat("expoente medido de 1/eps:\n"); print(round(expE, 2))
ok8 <- all(abs(expE[, ncol(expE)] - (sE - 0.5)) < 0.15)
cat(sprintf("    esperado s - 1/pi = %s; medido no último par %s -> %s\n",
            paste(sE - 0.5, collapse = ", "),
            paste(round(expE[, ncol(expE)], 2), collapse = ", "),
            if (ok8) "ok" else "FALHA"))

## --- veredito --------------------------------------------------------------
cat("\n")
if (ok1 && ok2 && ok3 && ok4 && ok5 && ok6 && ok7 && ok8) cat("OK\n") else
  stop("FALHA na conferência de E1.3 / E1.3b")
