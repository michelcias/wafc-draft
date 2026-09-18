# Literatura

Uma linha por trabalho: o que ele resolve, o que **não** resolve para nós, e o
status de verificação. Nenhuma entrada vai ao `.bib` sem status `verificado`
(título, autores, ano, veículo conferidos na fonte, de preferência Crossref).

Legenda de status: `web` (dados bibliográficos vistos em busca na web em
2026-09-18, sem Crossref; entrada provisória em
[`referencias-verificadas.bib`](referencias-verificadas.bib) marcada com
`note = {[L1: confirmar]}`), `resumo` (só o resumo foi lido), `lido` (o PDF
foi lido neste projeto), `no wall` (já está no `references_theo_1.bib` do
`wall-manuscript`, verificado lá; copiar de lá), `verificado` (conferido em
L1), `[VERIFICAR]` (citado de memória; nada de `.bib` até L1).

---

## Coeficientes funcionais e coeficientes aditivos

| Referência | O que resolve | O que não cobre | Status |
|---|---|---|---|
| Hastie & Tibshirani (1993), *JRSS B* 55(4) 757–796, "Varying-coefficient models" | o modelo; backfitting | sem estrutura aditiva nos coeficientes; sem penalização | `[VERIFICAR]` |
| Cai, Fan & Yao (2000), *JASA* 95(451) 941–956, "Functional-coefficient regression models for nonlinear time series" | polinômios locais para coeficientes funcionais; a nomenclatura "functional coefficient" | uma moduladora; sem seleção | `[VERIFICAR]` |
| Fan & Zhang (2008), *Stat. Interface* 1, 179–195, "Statistical methods with varying coefficient models" | revisão | | `[VERIFICAR]` |
| **Xue & Yang (2006)**, *Statistica Sinica* 16(4) 1423–1446, "Additive coefficient modeling via polynomial spline" | o additive coefficient model; splines; taxa univariada; BIC para seleção | suavidade global; sem wavelets; sem LASSO | `web` |
| Xue & Liang (2010), *Scand. J. Statist.*, "Polynomial spline estimation for a generalized additive coefficient model" | versão GLM | idem | `web` (DOI 10.1111/j.1467-9469.2009.00655.x) |
| Liu & Yang (2010), *Econometric Theory*, "Spline-backfitted kernel smoothing of additive coefficient model" | oráculo por componente | idem | `[VERIFICAR]` |
| Antoniadis, Gijbels & Verhasselt (2014?), *Statistical Papers*, "Penalized estimation in additive varying coefficient models using grouped regularization" | penalização em grupos no modelo aditivo de coeficientes | splines/P-splines; não wavelets | `web` (DOI 10.1007/s00362-013-0522-1; autores `[VERIFICAR]`) |
| **Wei, Huang & Li (2011)**, *Statistica Sinica* 21(4) 1515–1540, "Variable selection and estimation in high-dimensional varying-coefficient models" | group LASSO adaptativo com B-splines; propriedade oráculo de seleção | uma moduladora; splines | `web` |
| **Xue & Qu (2012)**, *JMLR* 13, "Variable selection in high-dimensional varying-coefficient models with global optimality" | seleção com otimalidade global | idem | `web` |
| Wang & Xia (2009), *JASA* 104(486) 747–757, "Shrinkage estimation of the varying coefficient model" (KLASSO) | kernel + LASSO | idem | `[VERIFICAR]` |
| Wang, Li & Huang (2008), *JASA* 103(484) 1556–1569 | SCAD em grupos, medidas repetidas | idem | `[VERIFICAR]` |

## Wavelets em regressão e em coeficientes variáveis

| Referência | O que resolve | O que não cobre | Status |
|---|---|---|---|
| Donoho & Johnstone (1994), *Biometrika* 81(3) 425–455; (1998), *Ann. Statist.* 26(3) 879–921 | limiarização; adaptação minimax em Besov | desenho equiespaçado | `no wall` |
| Antoniadis & Fan (2001), *JASA* 96(455) 939–967, "Regularization of wavelet approximations" | regularização (inclusive `ℓ_1`) com desenho não equiespaçado; oráculo; minimax adaptativo | regressão univariada | `web` |
| Cohen, Daubechies & Vial (1993), *ACHA* 1, 54–81 | wavelets no intervalo | | `no wall` (o `WaveBased` implementa) |
| Daubechies & Lagarias (1992) | avaliação de `φ`, `ψ` em pontos arbitrários | | `no wall` |
| **Zhou & You (2004)**, *Statist. Probab. Lett.* 68(1) 91–104, "Wavelet estimation in varying-coefficient partially linear regression models" | wavelets em coeficientes variáveis parcialmente lineares; normalidade assintótica sem suavidade forte | sem estrutura aditiva; sem LASSO; uma moduladora | `web` (DOI 10.1016/j.spl.2004.01.018) |
| Kovac & Silverman (2000), *JASA* 95(449) 172–183 | wavelets com desenho irregular por interpolação | | `[VERIFICAR]` |
| Amato, Antoniadis, De Feis & Gijbels (2019?), "Wavelet regression and additive models for irregularly spaced data" (arXiv 1903.04631) | aditivo com wavelets e desenho irregular | sem `X_j` | `web`; `[VERIFICAR]` autores e veículo |
| Zhang & Wong (2003), *Ann. Statist.* 31(1) 152–173, "Wavelet threshold estimation for additive regression models" | limiarização em modelos aditivos | | `[VERIFICAR]` |

## Modelos aditivos esparsos e LASSO com bases

| Referência | O que resolve | O que não cobre | Status |
|---|---|---|---|
| **Sardy & Tseng (2004)**, *JCGS* 13(2) 283–309, "AMlet, RAMlet, and GAMlet: automatic nonlinear fitting of additive models, robust and generalized, with wavelets" | aditivo com wavelets e `ℓ_1` | sem coeficientes funcionais | `web` |
| **Sardy & Ma (2024)**, *Scand. J. Statist.*, "Sparse additive models in high dimensions with wavelets" (arXiv 2207.04083, 2022) | aditivo esparso em alta dimensão com wavelets; convexo; sem CV | sem `X_j`: o mais próximo do nosso | `web`; volume e páginas `[VERIFICAR]` |
| Meier, van de Geer & Bühlmann (2009), *Ann. Statist.* 37(6B) 3779–3821, "High-dimensional additive modeling" | oráculo para aditivos com penalidade de suavidade e esparsidade | | `no wall`? conferir |
| Ravikumar, Lafferty, Liu & Wasserman (2009), *JRSS B* 71(5) 1009–1030 (SpAM) | aditivo esparso por backfitting | | `[VERIFICAR]` |
| Huang, Horowitz & Wei (2010), *Ann. Statist.* 38(4) 2282–2313, "Variable selection in nonparametric additive models" | group LASSO adaptativo com B-splines; seleção consistente | | `web` (arXiv 1010.4115) |
| Bühlmann & van de Geer (2011), *Statistics for High-Dimensional Data*, Springer | compatibilidade, oráculo (cap. 6), aditivos (cap. 8) | | `no wall` |
| Simon, Friedman, Hastie & Tibshirani (2013), *JCGS* 22(2) 231–245, "A sparse-group lasso" | a penalidade em grupos com esparsidade interna | | `[VERIFICAR]` |
| Zou (2006), *JASA* 101(476) 1418–1429, "The adaptive lasso and its oracle properties" | pesos adaptativos | | `[VERIFICAR]` |
| Bickel, Ritov & Tsybakov (2009), *Ann. Statist.* 37(4) 1705–1732 | autovalor restrito; oráculo | | `no wall`? conferir |
| Lounici, Pontil, van de Geer & Tsybakov (2011), *Ann. Statist.* 39(4) 2164–2204, "Oracle inequalities and optimal inference under group sparsity" | oráculo para group LASSO | | `web` (arXiv 1007.1771) |
| Klopp & Pensky (2013?), "Sparse high-dimensional varying coefficient model: non-asymptotic minimax study" (arXiv 1312.4087) | minimax não assintótico em coeficientes variáveis esparsos; block thresholding LASSO | uma moduladora; sem aditividade | `web`; `[VERIFICAR]` veículo |

## O próprio grupo

| Referência | O que resolve | Status |
|---|---|---|
| WALL teórico (`wall-manuscript/manuscript/theo/`) | sieve LASSO com wavelets no modelo aditivo logístico: aproximação em Besov, oráculo de taxa lenta, trilha rápida com compatibilidade e corolário de compressibilidade | `lido` (resumo e estrutura, 2026-09-18); é o molde |
| WALL aplicado (`manuscript/app/`) | benchmark do `wall()` em 14 bases | `lido` (estrutura) |
| Montoril, Pinheiro & Vidakovic (2019), *Scand. J. Statist.* 46(1) 215–234, "Wavelet-based estimators for mixture regression" | wavelets em regressão de mistura | `verificado` (`inst/CITATION` do `WaveBased`) |
| Motta & Montoril (2026), *Comm. Statist. Simul. Comput.* 55(6) 2426–2434 | mistura bayesiana com wavelets | `verificado` (idem) |

## Buscas pendentes (L2)

1. "additive coefficient model" + wavelet; "varying coefficient" + wavelet +
   lasso; "functional coefficient" + wavelet.
2. Varredura da *Statistica Sinica* e da EJS desde 2015 por título com
   "varying coefficient", "additive coefficient", "wavelet".
3. Sardy & Ma (2024): ler o PDF inteiro; registrar o que a teoria deles cobre
   e o que a condição de desenho de produtos acrescenta.
4. Antoniadis, Gijbels & Verhasselt (Statistical Papers): ler; confirmar se a
   "grouped regularization" ali é em B-splines e se há teoria.
