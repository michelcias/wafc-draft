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
| Hastie & Tibshirani (1993), *JRSS B* 55(4) 757–779 (discussão 779–796), "Varying-coefficient models" | o modelo; backfitting | sem estrutura aditiva nos coeficientes; sem penalização | `verificado` (`Hastie-Tibshirani-1993`) |
| Cai, Fan & Yao (2000), *JASA* 95(451) 941–956, "Functional-coefficient regression models for nonlinear time series" | polinômios locais para coeficientes funcionais; a nomenclatura "functional coefficient" | uma moduladora; sem seleção | `verificado` (`Cai-Fan-Yao-2000`) |
| Fan & Zhang (2008), *Stat. Interface* 1(1) 179–195, "Statistical methods with varying coefficient models" | revisão | | `verificado` (`Fan-Zhang-2008`) |
| **Xue & Yang (2006)**, *Statistica Sinica* 16(4) 1423–1446, "Additive coefficient modeling via polynomial spline" | o additive coefficient model; splines; taxa univariada; BIC para seleção | suavidade global; sem wavelets; sem LASSO | `verificado` (`Xue-Yang-2006`; sem DOI: a SS não registrou 2006 no Crossref; conferido no PDF da editora) |
| Xue & Liang (2010), *Scand. J. Statist.* 37(1) 26–46, "Polynomial spline estimation for a generalized additive coefficient model" | versão GLM | idem | `verificado` (`Xue-Liang-2010`) |
| Liu & Yang (2010), *Econometric Theory* 26(1) 29–59, "Spline-backfitted kernel smoothing of additive coefficient model" | oráculo por componente | idem | `verificado` (`Liu-Yang-2010`) |
| Antoniadis, Gijbels & Lambert-Lacroix (2014), *Statistical Papers* 55(3) 727–750, "Penalized estimation in additive varying coefficient models using grouped regularization" | penalização em grupos no modelo aditivo de coeficientes | splines/P-splines; não wavelets | `verificado` (`Antoniadis-Gijbels-LambertLacroix-2014`; o terceiro autor é Lambert-Lacroix, não Verhasselt) |
| **Wei, Huang & Li (2011)**, *Statistica Sinica* 21(4) 1515–1540, "Variable selection and estimation in high-dimensional varying-coefficient models" | group LASSO adaptativo com B-splines; propriedade oráculo de seleção | uma moduladora; splines | `verificado` (`Wei-Huang-Li-2011`) |
| **Xue & Qu (2012)**, *JMLR* 13(63) 1973–1998, "Variable selection in high-dimensional varying-coefficient models with global optimality" | seleção com otimalidade global | idem | `verificado` (`Xue-Qu-2012`; sem DOI, conferido na página do JMLR) |
| Wang & Xia (2009), *JASA* 104(486) 747–757, "Shrinkage estimation of the varying coefficient model" (KLASSO) | kernel + LASSO | idem | `verificado` (`Wang-Xia-2009`) |
| Wang, Li & Huang (2008), *JASA* 103(484) 1556–1569, "Variable selection in nonparametric varying-coefficient models for analysis of repeated measurements" | SCAD em grupos, medidas repetidas | idem | `verificado` (`Wang-Li-Huang-2008`) |

## Wavelets em regressão e em coeficientes variáveis

| Referência | O que resolve | O que não cobre | Status |
|---|---|---|---|
| Donoho & Johnstone (1994), *Biometrika* 81(3) 425–455, "Ideal spatial adaptation by wavelet shrinkage"; (1998), *Ann. Statist.* 26(3) 879–921, "Minimax estimation via wavelet shrinkage" | limiarização; adaptação minimax em Besov | desenho equiespaçado | `verificado` (`Donoho-Johnstone-1994`, `Donoho-Johnstone-1998`; nenhum dos dois está no WALL, que cita só o de *PTRF* 1994) |
| Antoniadis & Fan (2001), *JASA* 96(455) 939–967, "Regularization of wavelet approximations" | regularização (inclusive `ℓ_1`) com desenho não equiespaçado; oráculo; minimax adaptativo | regressão univariada | `verificado` (`Antoniadis-Fan-2001`) |
| Cohen, Daubechies & Vial (1993), *ACHA* 1(1) 54–81 | wavelets no intervalo | | `verificado` (`cohen1993wavelets`, copiada do WALL; o `WaveBased` implementa) |
| Daubechies & Lagarias (1991), *SIAM J. Math. Anal.* 22(5) 1388–1410; (1992), 23(4) 1031–1079 | avaliação de `φ`, `ψ` em pontos arbitrários | | `verificado` (`Daubechies-Lagarias-1991`, `Daubechies-Lagarias-1992`; não está no WALL; o algoritmo de produtos de matrizes é da parte I) |
| **Zhou & You (2004)**, *Statist. Probab. Lett.* 68(1) 91–104, "Wavelet estimation in varying-coefficient partially linear regression models" | wavelets em coeficientes variáveis parcialmente lineares; normalidade assintótica sem suavidade forte | sem estrutura aditiva; sem LASSO; uma moduladora | `verificado` (`Zhou-You-2004`) |
| Kovac & Silverman (2000), *JASA* 95(449) 172–183, "Extending the scope of wavelet regression methods by coefficient-dependent thresholding" | wavelets com desenho irregular por interpolação | | `verificado` (`Kovac-Silverman-2000`) |
| Haris, Simon & Shojaie (2018), *NeurIPS* 31, 8987–8997, "Wavelet regression and additive models for irregularly spaced data" (arXiv 1903.04631); o candidato de Amato, Antoniadis, De Feis & Gijbels é (2022), *Stat. Comput.* 32(1) art. 11, "Wavelet-based robust estimation and variable selection in nonparametric additive models" | aditivo com wavelets e desenho irregular | sem `X_j` | `verificado` (`Haris-Simon-Shojaie-2018`, `Amato-Antoniadis-DeFeis-Gijbels-2022`; o arXiv citado é de Haris et al., não de Amato et al.) |
| Zhang & Wong (2003), *Ann. Statist.* 31(1) 152–173, "Wavelet threshold estimation for additive regression models" | limiarização em modelos aditivos | | `verificado` (`Zhang-Wong-2003`) |

## Modelos aditivos esparsos e LASSO com bases

| Referência | O que resolve | O que não cobre | Status |
|---|---|---|---|
| **Sardy & Tseng (2004)**, *JCGS* 13(2) 283–309, "AMlet, RAMlet, and GAMlet: automatic nonlinear fitting of additive models, robust and generalized, with wavelets" | aditivo com wavelets e `ℓ_1` | sem coeficientes funcionais | `verificado` (`Sardy-Tseng-2004`) |
| **Sardy & Ma (2024)**, *Scand. J. Statist.* 51(1) 89–108, "Sparse additive models in high dimensions with wavelets" (arXiv 2207.04083, 2022) | aditivo esparso em alta dimensão com wavelets; convexo; sem CV | sem `X_j`: o mais próximo do nosso | `verificado` (`Sardy-Ma-2024`; o WALL cita só o preprint) |
| Meier, van de Geer & Bühlmann (2009), *Ann. Statist.* 37(6B) 3779–3821, "High-dimensional additive modeling" | oráculo para aditivos com penalidade de suavidade e esparsidade | | `verificado` (`meier2009high`, copiada do WALL) |
| Ravikumar, Lafferty, Liu & Wasserman (2009), *JRSS B* 71(5) 1009–1030, "Sparse additive models" (SpAM) | aditivo esparso por backfitting | | `verificado` (`ravikumar2009sparse`, copiada do WALL) |
| Huang, Horowitz & Wei (2010), *Ann. Statist.* 38(4) 2282–2313, "Variable selection in nonparametric additive models" | group LASSO adaptativo com B-splines; seleção consistente | | `verificado` (`Huang-Horowitz-Wei-2010`) |
| Bühlmann & van de Geer (2011), *Statistics for High-Dimensional Data*, Springer | compatibilidade, oráculo (cap. 6), aditivos (cap. 8) | | `verificado` (`buhlmann2011statistics`, copiada do WALL) |
| Simon, Friedman, Hastie & Tibshirani (2013), *JCGS* 22(2) 231–245, "A sparse-group lasso" | a penalidade em grupos com esparsidade interna | | `verificado` (`Simon-Friedman-Hastie-Tibshirani-2013`) |
| Zou (2006), *JASA* 101(476) 1418–1429, "The adaptive lasso and its oracle properties" | pesos adaptativos | | `verificado` (`Zou-2006`) |
| Bickel, Ritov & Tsybakov (2009), *Ann. Statist.* 37(4) 1705–1732, "Simultaneous analysis of lasso and Dantzig selector" | autovalor restrito; oráculo | | `verificado` (`Bickel-Ritov-Tsybakov-2009`; não está no WALL) |
| Lounici, Pontil, van de Geer & Tsybakov (2011), *Ann. Statist.* 39(4) 2164–2204, "Oracle inequalities and optimal inference under group sparsity" | oráculo para group LASSO | | `verificado` (`Lounici-Pontil-vandeGeer-Tsybakov-2011`) |
| Klopp & Pensky (2015), *Ann. Statist.* 43(3) 1273–1299, "Sparse high-dimensional varying coefficient model: nonasymptotic minimax study" (arXiv 1312.4087) | minimax não assintótico em coeficientes variáveis esparsos; block thresholding LASSO | uma moduladora; sem aditividade | `verificado` (`Klopp-Pensky-2015`) |

## O próprio grupo

| Referência | O que resolve | Status |
|---|---|---|
| WALL teórico (`wall-manuscript/manuscript/theo/`) | sieve LASSO com wavelets no modelo aditivo logístico: aproximação em Besov, oráculo de taxa lenta, trilha rápida com compatibilidade e corolário de compressibilidade | `lido` (resumo e estrutura, 2026-09-18); é o molde |
| WALL aplicado (`manuscript/app/`) | benchmark do `wall()` em 14 bases | `lido` (estrutura) |
| Montoril, Pinheiro & Vidakovic (2019), *Scand. J. Statist.* 46(1) 215–234, "Wavelet-based estimators for mixture regression" | wavelets em regressão de mistura | `verificado` (`Montoril-Pinheiro-Vidakovic-2019`; Crossref em L1) |
| Motta & Montoril (2026), *Comm. Statist. Simul. Comput.* 55(6) 2426–2434 | mistura bayesiana com wavelets | `verificado` (`Motta-Montoril-2026`; Crossref em L1) |

## Buscas pendentes (L2)

1. "additive coefficient model" + wavelet; "varying coefficient" + wavelet +
   lasso; "functional coefficient" + wavelet.
2. Varredura da *Statistica Sinica* e da EJS desde 2015 por título com
   "varying coefficient", "additive coefficient", "wavelet".
3. Sardy & Ma (2024): ler o PDF inteiro; registrar o que a teoria deles cobre
   e o que a condição de desenho de produtos acrescenta.
4. Antoniadis, Gijbels & Verhasselt (Statistical Papers): ler; confirmar se a
   "grouped regularization" ali é em B-splines e se há teoria.
