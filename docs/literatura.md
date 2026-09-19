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
L1), `[VERIFICAR]` (citado de memória; nada de `.bib` até L1). **A segunda rodada (L3, 2026-09-19) fechou os pendentes**: nenhuma linha abaixo continua com status `resumo` ou `[VERIFICAR]`, e onde a conferência mudou um dado (ano, autor, páginas) a mudança está anotada na própria linha.

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

| Ma, Carroll, Liang & Xu (2015), *Ann. Statist.* 43(5) 2102–2131, "Estimation and inference in generalized additive coefficient models for nonlinear interactions with high-dimensional covariates" (DOI 10.1214/15-AOS1344) | GACM em alta dimensão: penalização em grupos, identificação de estrutura consistente, bandas de confiança simultâneas, splines em dois passos | splines; sem wavelets; sem adaptatividade | `verificado` (`Ma-Carroll-Liang-Xu-2015`; páginas 2102–2131 conferidas no Project Euclid) |
| Wang, Jiang & Liu (2024), *JCGS* 33(2), "Varying coefficient model via adaptive spline fitting" (DOI 10.1080/10618600.2023.2267616; arXiv 2201.10063) | splines com nós adaptativos por coeficiente; programação dinâmica; EQM menor que nós equiespaçados | uma moduladora; sem teoria; sem esparsidade; é o concorrente adaptativo para E4 | `verificado` (`Wang-Jiang-Liu-2024`; páginas 614–624; on-line em 2023, fascículo 33(2) em 2024) |
| Tibshirani & Friedman (2020), *JCGS* 29(1) 215–225, "A pliable lasso" (DOI 10.1080/10618600.2019.1648271) | LASSO com coeficientes modificados linearmente por variáveis modificadoras | modificação linear; sem wavelets | `verificado` (`Tibshirani-Friedman-2020`) |
| Deshpande, Bai, Balocchi, Starling & Weiss (2026), *Bayesian Anal.* 21(1) 281–308, "VCBART: Bayesian trees for varying coefficients" (DOI 10.1214/24-BA1470) | coeficientes variáveis em várias moduladoras por árvores bayesianas; adaptativo sem aditividade | sem wavelets; sem teoria de taxa; concorrente não aditivo para E4 | `verificado` (`Deshpande-Bai-Balocchi-Starling-Weiss-2026`; **o ano é 2026, não 2024**: o `24` do DOI é o ano de aceitação e o fascículo saiu em março de 2026) |

## Wavelets em regressão e em coeficientes variáveis

| Referência | O que resolve | O que não cobre | Status |
|---|---|---|---|
| Donoho & Johnstone (1994), *Biometrika* 81(3) 425–455, "Ideal spatial adaptation by wavelet shrinkage"; (1998), *Ann. Statist.* 26(3) 879–921, "Minimax estimation via wavelet shrinkage" | limiarização; adaptação minimax em Besov | desenho equiespaçado | `verificado` (`Donoho-Johnstone-1994`, `Donoho-Johnstone-1998`; nenhum dos dois está no WALL, que cita só o de *PTRF* 1994) |
| Johnstone (2019), *Gaussian Estimation: Sequence and Wavelet Models*, versão de rascunho de 16 de setembro de 2019 | o livro do modelo de sequência gaussiano: minimax sobre bolas `ℓ_τ` fracas, limiarização e a tradução entre sequência e função; é a referência de apoio da observação de taxa ao lado de Donoho & Johnstone (1998) | não é regressão com coeficientes; sem desenho aleatório | `verificado` (`johnstone2019gaussian`, copiada do WALL; **não é livro publicado**: sem registro no Crossref e zero resultados na busca em livros da Cambridge University Press em 2026-09-19, ao contrário do `publisher` que a entrada do WALL traz; a fonte é a página do autor, que oferece o rascunho de 2019-09-16) |
| Antoniadis & Fan (2001), *JASA* 96(455) 939–967, "Regularization of wavelet approximations" | regularização (inclusive `ℓ_1`) com desenho não equiespaçado; oráculo; minimax adaptativo | regressão univariada | `verificado` (`Antoniadis-Fan-2001`) |
| Cohen, Daubechies & Vial (1993), *ACHA* 1(1) 54–81 | wavelets no intervalo | | `verificado` (`cohen1993wavelets`, copiada do WALL; o `WaveBased` implementa) |
| Daubechies & Lagarias (1991), *SIAM J. Math. Anal.* 22(5) 1388–1410; (1992), 23(4) 1031–1079 | avaliação de `φ`, `ψ` em pontos arbitrários | | `verificado` (`Daubechies-Lagarias-1991`, `Daubechies-Lagarias-1992`; não está no WALL; o algoritmo de produtos de matrizes é da parte I) |
| **Zhou & You (2004)**, *Statist. Probab. Lett.* 68(1) 91–104, "Wavelet estimation in varying-coefficient partially linear regression models" | wavelets em coeficientes variáveis parcialmente lineares; normalidade assintótica sem suavidade forte | sem estrutura aditiva; sem LASSO; uma moduladora | `verificado` (`Zhou-You-2004`) |
| Kovac & Silverman (2000), *JASA* 95(449) 172–183, "Extending the scope of wavelet regression methods by coefficient-dependent thresholding" | wavelets com desenho irregular por interpolação | | `verificado` (`Kovac-Silverman-2000`) |
| Haris, Simon & Shojaie (2018), *NeurIPS* 31, 8987–8997, "Wavelet regression and additive models for irregularly spaced data" (arXiv 1903.04631); o candidato de Amato, Antoniadis, De Feis & Gijbels é (2022), *Stat. Comput.* 32(1) art. 11, "Wavelet-based robust estimation and variable selection in nonparametric additive models" | aditivo com wavelets e desenho irregular | sem `X_ℓ` | `verificado` (`Haris-Simon-Shojaie-2018`, `Amato-Antoniadis-DeFeis-Gijbels-2022`; o arXiv citado é de Haris et al., não de Amato et al.) |
| Zhang & Wong (2003), *Ann. Statist.* 31(1) 152–173, "Wavelet threshold estimation for additive regression models" | limiarização em modelos aditivos | | `verificado` (`Zhang-Wong-2003`) |

| **Montoril, Morettin & Chiann (2018)**, *Int. J. Wavelets Multiresolut. Inf. Process.* 16(1) 1850004, "Wavelet estimation of functional coefficient regression models" (DOI 10.1142/S0219691318500042) | wavelets clássicas e *warped* em coeficientes funcionais de séries temporais; taxas; AIC/BIC para os níveis; previsão | uma moduladora; mínimos quadrados sem penalização; sem aditividade; precursor do grupo, obrigatório citar | `verificado` (`Montoril-Morettin-Chiann-2018`) |
| Zhou, Xu & Lin (2017), *Statist. Probab. Lett.* 122, 179–189 (DOI 10.1016/j.spl.2016.11.009); Zhou, Xu & Lin (2018), *Comm. Statist. Theory Methods* 47(10) 2504–2519 (DOI 10.1080/03610926.2017.1339801); Zhou, Ni & Zhu (2019), *Lith. Math. J.* 59(2) 276–293 (DOI 10.1007/s10986-019-09440-1); Zhou, Yang & Xiang (2022), *Mathematics* 10(13) 2321 (DOI 10.3390/math10132321) | estimadores lineares por núcleo de wavelet em coeficientes variáveis no tempo: censura, erro de medida, α-mixing, quantis; normalidade assintótica | uma moduladora; sem limiarização nem penalização; sem aditividade | `verificado` (`Zhou-Xu-Lin-2017`, `Zhou-Xu-Lin-2018`, `Zhou-Ni-Zhu-2019`, `Zhou-Yang-Xiang-2022`; **o primeiro é de 2017**, não 2016, porque o fascículo 122 é de março de 2017; **o terceiro autor do de 2022 é Yu Xiang**, não "Yu"); escolher uma ou duas para citar |
| Zhao, Ogden & Reiss (2012), *JCGS* 21(3) 600–617, "Wavelet-based LASSO in functional linear regression" (DOI 10.1080/10618600.2012.679241); Yu, Zhang, Mizera, Jiang & Kong (2019), *Comput. Statist. Data Anal.* 136, 12–29, "Sparse wavelet estimation in quantile regression with multiple functional predictors" (DOI 10.1016/j.csda.2018.12.002) | função-coeficiente de regressão funcional em wavelets com LASSO (2012) e sparse group LASSO por preditor (2019) | outro modelo (escalar sobre curva); sem coeficientes variáveis | `verificado` (`Zhao-Ogden-Reiss-2012`, `Yu-Zhang-Mizera-Jiang-Kong-2019`) |
| Schnaidt Grez & Vidakovic (2018), arXiv 1803.04558 e 1804.03015 | mínimos quadrados com wavelets periódicas em aditivos com desenho aleatório; consistência e taxas | sem `X_ℓ`; sem LASSO; sem veículo | `verificado` (`SchnaidtGrez-Vidakovic-2018a`, `SchnaidtGrez-Vidakovic-2018b`; **continuam sem veículo** em 2026-09-19: a API do arXiv não registra `journal_ref` nem DOI para nenhum dos dois); só se a lista de aditivos com wavelets tiver de ser completa |

## Modelos aditivos esparsos e LASSO com bases

| Referência | O que resolve | O que não cobre | Status |
|---|---|---|---|
| **Tibshirani (1996)**, *JRSS B* 58(1) 267–288, "Regression shrinkage and selection via the lasso" | o lasso: a penalidade `‖θ‖_1`, a seleção que vem junto com o encolhimento | desenho fixo e de dimensão pequena; sem bases; sem teoria de alta dimensão | `verificado` (`tibshirani1996regression`, copiada do WALL; DOI acrescentado em L4) |
| **Sardy & Tseng (2004)**, *JCGS* 13(2) 283–309, "AMlet, RAMlet, and GAMlet: automatic nonlinear fitting of additive models, robust and generalized, with wavelets" | aditivo com wavelets e `ℓ_1` | sem coeficientes funcionais | `verificado` (`Sardy-Tseng-2004`) |
| **Sardy & Ma (2024)**, *Scand. J. Statist.* 51(1) 89–108, "Sparse additive models in high dimensions with wavelets" (arXiv 2207.04083, 2022) | aditivo esparso em alta dimensão com wavelets; convexo; sem CV | sem `X_ℓ`: o mais próximo no método | `verificado` (`Sardy-Ma-2024`; o WALL cita só o preprint) |
| Meier, van de Geer & Bühlmann (2009), *Ann. Statist.* 37(6B) 3779–3821, "High-dimensional additive modeling" | oráculo para aditivos com penalidade de suavidade e esparsidade | | `verificado` (`meier2009high`, copiada do WALL) |
| Ravikumar, Lafferty, Liu & Wasserman (2009), *JRSS B* 71(5) 1009–1030, "Sparse additive models" (SpAM) | aditivo esparso por backfitting | | `verificado` (`ravikumar2009sparse`, copiada do WALL) |
| Huang, Horowitz & Wei (2010), *Ann. Statist.* 38(4) 2282–2313, "Variable selection in nonparametric additive models" | group LASSO adaptativo com B-splines; seleção consistente | | `verificado` (`Huang-Horowitz-Wei-2010`) |
| Bühlmann & van de Geer (2011), *Statistics for High-Dimensional Data*, Springer | compatibilidade, oráculo (cap. 6), aditivos (cap. 8) | | `verificado` (`buhlmann2011statistics`, copiada do WALL) |
| Simon, Friedman, Hastie & Tibshirani (2013), *JCGS* 22(2) 231–245, "A sparse-group lasso" | a penalidade em grupos com esparsidade interna | | `verificado` (`Simon-Friedman-Hastie-Tibshirani-2013`) |
| Zou (2006), *JASA* 101(476) 1418–1429, "The adaptive lasso and its oracle properties" | pesos adaptativos | | `verificado` (`Zou-2006`) |
| Bickel, Ritov & Tsybakov (2009), *Ann. Statist.* 37(4) 1705–1732, "Simultaneous analysis of lasso and Dantzig selector" | autovalor restrito; oráculo | | `verificado` (`Bickel-Ritov-Tsybakov-2009`; não está no WALL) |
| Lounici, Pontil, van de Geer & Tsybakov (2011), *Ann. Statist.* 39(4) 2164–2204, "Oracle inequalities and optimal inference under group sparsity" | oráculo para group LASSO | | `verificado` (`Lounici-Pontil-vandeGeer-Tsybakov-2011`) |
| **Klopp & Pensky (2015)**, *Ann. Statist.* 43(3) 1273–1299, "Sparse high-dimensional varying coefficient model: nonasymptotic minimax study" (arXiv 1312.4087) | desenho de produtos com Gram `Ω ⊗ Φ`, concentração da Gram empírica restrita, block LASSO, oráculo não assintótico, taxa adaptativa em Besov com `ν < 2` e cota inferior minimax | uma moduladora; `W ⊥ t`; penalidade em blocos; sem aditividade, software, simulação ou aplicação | `lido` (L2, arXiv 1312.4087v2); `verificado` (`Klopp-Pensky-2015`). **É o trabalho mais próximo na teoria**; detalhe em [`busca-novidade.md`](busca-novidade.md) §1 |

| Dalalyan, Ingster & Tsybakov (2014), *Probab. Theory Related Fields* 158(3–4) 513–532, "Statistical inference in compound functional models" (arXiv 1208.6402) | minimax não assintótico em modelos compostos, inclusive aditivos esparsos | sem coeficientes variáveis | `verificado` (`Dalalyan-Ingster-Tsybakov-2014`; DOI 10.1007/s00440-013-0487-y) |
| Giacobino, Sardy, Diaz-Rodriguez & Hengartner (2017), *Electron. J. Statist.* 11(2) 4701–4722, "Quantile universal threshold" | a regra QUT de `λ` (pivotal, sem `σ`) que Sardy & Ma usam; candidata para E2.3 | | `verificado` (`Giacobino-Sardy-DiazRodriguez-Hengartner-2017`; DOI 10.1214/17-EJS1366) |

## O próprio grupo

| Referência | O que resolve | Status |
|---|---|---|
| WALL teórico (`wall-manuscript/manuscript/theo/`) | sieve LASSO com wavelets no modelo aditivo logístico: aproximação em Besov, oráculo de taxa lenta, trilha rápida com compatibilidade e corolário de compressibilidade | `lido` (resumo e estrutura, 2026-09-18); é o molde |
| WALL aplicado (`manuscript/app/`) | benchmark do `wall()` em 14 bases | `lido` (estrutura) |
| Montoril, Pinheiro & Vidakovic (2019), *Scand. J. Statist.* 46(1) 215–234, "Wavelet-based estimators for mixture regression" | wavelets em regressão de mistura | `verificado` (`Montoril-Pinheiro-Vidakovic-2019`; Crossref em L1) |
| Motta & Montoril (2026), *Comm. Statist. Simul. Comput.* 55(6) 2426–2434 | mistura bayesiana com wavelets | `verificado` (`Motta-Montoril-2026`; Crossref em L1) |

## Buscas (L2, fechada em 2026-09-18)

As quatro buscas foram feitas e estão em
[`busca-novidade.md`](busca-novidade.md), com uma tabela por busca, a
varredura da *Statistica Sinica* e da EJS desde 2015 (interseção
"coeficientes variáveis × wavelet": zero nas duas) e o veredito por
contribuição. Resultado curto: ninguém combina coeficientes aditivos em
várias moduladoras, wavelets e LASSO; o trabalho mais próximo na teoria é
Klopp & Pensky (2015), e no método são Sardy & Ma (2024) e Amato et al.
(2022).

As linhas que L2 acrescentou acima foram conferidas em L3 (2026-09-19),
como as de L1, e estão em
[`referencias-verificadas.bib`](referencias-verificadas.bib). A rodada
corrigiu três dados que L2 trazia errados: o ano de Deshpande et al. (2026,
não 2024), o ano de Zhou, Xu & Lin (2017, não 2016) e o terceiro autor de
Zhou, Yang & Xiang (2022).
