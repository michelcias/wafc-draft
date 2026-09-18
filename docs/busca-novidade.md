# Busca de novidade (L2)

Chat de tarefa, 2026-09-18. Entregável do plano (`plano-projeto.md`, L2):
as quatro buscas de `literatura.md` ("Buscas pendentes"), uma tabela por
busca, a leitura de Sardy & Ma (2024), a varredura da *Statistica Sinica* e
da EJS desde 2015 e um veredito por contribuição de `alvo-revista.md` §4.
Notação de `notacao.md` (congelada em E1.1): linear `X_ℓ`, moduladora `U_m`,
wavelet `ψ_{jk}`, coeficiente `θ_{ℓm,jk}`.

**Como foi feito.** Busca na web (títulos e resumos); Crossref por ISSN para
a varredura das revistas; OpenAlex para o grafo de citações dos trabalhos-
âncora (quem cita Zhou & You 2004, Klopp & Pensky 2015, Sardy & Tseng 2004,
Sardy & Ma 2024, Antoniadis et al. 2014, Wei, Huang & Li 2011, Montoril et
al. 2018), filtrando título e resumo por "wavelet", "lasso", "sparse",
"Besov" e "additive coefficient"; arXiv para os textos integrais lidos. Os
metadados citados abaixo foram vistos no Crossref ou no OpenAlex hoje; a
verificação formal é de L1 (chaves do `.bib` entre parênteses quando já
existem).

**Limites.** O texto publicado de Sardy & Ma (Wiley) está atrás de uma
verificação anti-robô e o de Antoniadis et al. (Springer) atrás de paywall,
sem cópia aberta em Unpaywall, CORE, HAL nem no site da autora (404). A
varredura das revistas é por título, não por resumo. Google Scholar não
está acessível pelas ferramentas; periódicos em chinês não foram cobertos.

Legenda da coluna "ameaça": **alta** (cobre parte do que o artigo diz que é
seu), **média** (o referee vai pedir a comparação ou a citação), **baixa**
(citar em uma linha ou nem isso).

---

## 1. Busca 1: "additive coefficient" + wavelet; "varying coefficient" + wavelet + lasso; "functional coefficient" + wavelet

Resultado seco: **nenhum trabalho combina coeficientes aditivos em várias
moduladoras, base de wavelets e LASSO.** "Additive coefficient model" +
wavelet devolve só modelos aditivos (sem `X_ℓ`). O que existe está em três
famílias: (a) LASSO em blocos sobre o desenho de produtos com **uma**
moduladora e teoria minimax (Klopp & Pensky); (b) wavelets lineares (sem
limiarização) em coeficientes variáveis no tempo (Zhou & You e a escola de
Zhou Xingcai; o próprio autor com Morettin & Chiann); (c) wavelets + LASSO
em regressão **funcional** (escalar sobre curva), que é outro modelo.

| Trabalho | O que faz | O que não faz | Ameaça à novidade | Status |
|---|---|---|---|---|
| **Klopp & Pensky (2015)**, *Ann. Statist.* 43(3) 1273–1299, DOI 10.1214/15-AOS1309 (`Klopp-Pensky-2015`) | modelo `Y = W'f(t) + σξ`, `W ⊥ t`, `t ∈ [0,1]` com densidade limitada; expande cada `f_ℓ` numa base ortonormal de `L_2[0,1]` (a wavelet periódica é o exemplo (b) da hipótese A1); desenho de produtos `B_i = Vec(φ(t_i) W_i')`, Gram populacional `Σ = Ω ⊗ Φ` com `Ω = E[WW']`, `Φ = E[φφ']` (eq. 1.8); **block LASSO** com blocos de tamanho `≈ log n` dentro de cada `f_ℓ` e o coeficiente da constante como bloco próprio (norma 3.1); hipótese A5 "restricted isometry in expectation" sobre `Ω` restrita a `|Λ| ≤ ℵ`; **Lema 1**: concentração da Gram empírica restrita em torno de `Ω_Λ ⊗ Φ` com `n ≥ N(ℵ) ≍ ℵ (L+1) log(p+L)`; esparsidade/regularidade por `Σ_k |a_{ℓk}|^ν (k+1)^{ν r'} ≤ C^ν` (Besov quando a base é wavelet, `ν < 2` = não homogênea); **cotas inferiores minimax** (Teorema 1) e **cota superior adaptativa** não assintótica com `δ ≍ σ sqrt(log p / n)` (Teorema 2), ótima a menos de `log p` (Corolário 1); erro sub-gaussiano; risco em esperança com restrição extra (Teorema 3) | **uma moduladora** (`t` escalar; sem soma em `m`); assume `W ⊥ t` (o caso "X dependente de U" fica fora); penaliza em blocos, não coeficiente a coeficiente; a constante entra na penalidade; exige `L + 1 ≥ n^{1/2}` e `r ≥ 2`; sem `L_2(P_U)` com densidade (usa `Φ` com autovalores limitados); sem software, sem simulação, sem aplicação; não trata `j_0`, periodização nem identificabilidade no nível da base | **alta** para as contribuições 1 e 2 no caso `q = 1`, `X ⊥ U`: a fatoração `E[XX'] ⊗ E[ψψ']`, a concentração da Gram restrita e a taxa adaptativa em Besov com desigualdade oráculo **já estão lá**, inclusive com cota inferior. O WAFC precisa se posicionar como extensão (aditividade em `q ≥ 2` moduladoras, `X` dependente de `U`, LASSO puro com corolário de compressibilidade) e citar K&P como base do caso `q = 1` | `lido` (arXiv 1312.4087v2, 2014-05-15; é a versão do AoS a menos de revisão) |
| **Montoril, Morettin & Chiann (2018)**, *Int. J. Wavelets Multiresolut. Inf. Process.* 16(1) 1850004, DOI 10.1142/S0219691318500042 | estimação por wavelets (clássicas e *warped*) do modelo de coeficientes funcionais para séries temporais não lineares; taxas de convergência; AIC/AICc/BIC para escolher os níveis grosso e fino; previsão multi-passos | uma moduladora; mínimos quadrados sem penalização; sem aditividade; sem desigualdade oráculo | **baixa** como concorrente; **obrigatório citar**: é o precursor do próprio grupo e o referee o encontra pelo nome do autor | `resumo` (OpenAlex); metadados no Crossref |
| Zhou, Xu & Lin (2016), *Statist. Probab. Lett.* 122, 179–189, DOI 10.1016/j.spl.2016.11.009; Zhou, Xu & Lin (2018), *Comm. Statist. Theory Methods* 47(10) 2504–2519, DOI 10.1080/03610926.2017.1339801; Zhou, Ni & Zhu (2019), *Lith. Math. J.* 59(2) 276–293, DOI 10.1007/s10986-019-09440-1; Zhou & Zhu (2020), *Discrete Dyn. Nat. Soc.* 2020, 1025452; Zhou & Lv (2021), *Anal. Appl.*, DOI 10.1142/S0219530521500032; Zhou, Yang & Yu (2022), *Mathematics* 10(13) 2321, DOI 10.3390/math10132321 | a escola que cita Zhou & You (2004): estimadores lineares por "núcleo de wavelet" (projeção em `V_J`) em modelos de coeficiente variável no tempo, com erros de medida, censura, dependência α-mixing, M-estimação e quantis; normalidade assintótica e cotas de Berry–Esseen | uma moduladora (tempo); sem penalização nem limiarização (logo sem adaptatividade espacial nem Besov com `π < 2`); sem aditividade; sem alta dimensão | **baixa**; uma frase na introdução ("wavelet estimators without thresholding in time-varying coefficient models") com duas ou três citações | `resumo` (2018, 2019, 2022); os demais só título e metadados |
| Zhao, Ogden & Reiss (2012), *JCGS* 21(3) 600–617, DOI 10.1080/10618600.2012.679241; Zhao, Chen & Ogden (2015), *JCGS* 24(3) 655–675; Yu, Zhang, Mizera, Jiang & Kong (2019), *Comput. Statist. Data Anal.* 136, 12–29, DOI 10.1016/j.csda.2018.12.002 | regressão linear funcional (escalar sobre curva) com a função-coeficiente `β(t)` em wavelets e LASSO (2012, 2015); versão quantílica com vários preditores funcionais e **sparse group LASSO** (grupo = preditor funcional; LASSO dentro do grupo) (2019); argumentam exatamente a adaptatividade a traços locais | outro modelo: o `t` indexa a curva-preditora, não uma moduladora; sem coeficientes variáveis | **baixa** para o modelo; **média** para a retórica "wavelet + LASSO adapta a traços locais", que já é deles; a combinação SGL + wavelets bi-nível existe (2019) e a variante de E2.2 deve citá-la | `resumo` |
| Tibshirani & Friedman (2020), *JCGS* 29(1) 215–225, DOI 10.1080/10618600.2019.1648271, "A pliable lasso" | LASSO cujos coeficientes são modificados **linearmente** por variáveis modificadoras; conexão explícita com coeficientes variáveis | modificação linear, não não paramétrica; sem wavelets | **baixa**; citar como o caso paramétrico do mesmo problema | `resumo` |
| Wang, Jiang & Liu (2024), *JCGS* 33(2), DOI 10.1080/10618600.2023.2267616, "Varying coefficient model via adaptive spline fitting" (arXiv 2201.10063) | splines com nós **adaptativos e específicos por coeficiente**, escolhidos por programação dinâmica; EQM menor que splines equiespaçados | uma moduladora; sem teoria de taxa; sem esparsidade | **média** para a contribuição 3: é o spline que adapta localmente; se E4 comparar só com `mgcv`, o referee pergunta por este | `resumo` |
| Haris, Simon & Shojaie (2018), *NeurIPS* 31, 8987–8997 (`Haris-Simon-Shojaie-2018`), waveMesh | wavelets em desenho irregular por malha de interpolação; gradiente proximal; **taxas minimax adaptativas**; extensão a aditivos esparsos com taxa ótima **sob condição de compatibilidade fraca** | sem `X_ℓ`; sem coeficientes variáveis | **média** para a arquitetura de prova: "aditivo com wavelets + compatibilidade + taxa adaptativa" já existe; o que sobra para o WAFC é o desenho de produtos | `resumo` |
| Amato, Antoniadis, De Feis & Gijbels (2022), *Stat. Comput.* 32(1) art. 11 (`Amato-Antoniadis-DeFeis-Gijbels-2022`) | aditivos com wavelets em desenho não equiespaçado e `n` qualquer; componentes com regularidades distintas; **seleção bi-nível** (componentes e coeficientes dentro da componente) com penalidades não convexas em grupo; M-estimação robusta; taxas ótimas e consistência de seleção sob compatibilidade fraca | sem `X_ℓ`; penalidade não convexa | **média** para E2.2/E1.7: a seleção bi-nível com wavelets em aditivos existe; a novidade só pode ser "no desenho de produtos" | `resumo` (OpenAlex) |
| Schnaidt Grez & Vidakovic (2018), arXiv 1803.04558 e 1804.03015 | mínimos quadrados com wavelets periódicas em `[0,1]` para aditivos com desenho aleatório; consistência forte e taxas; versão com coeficientes empíricos limiarizados | sem `X_ℓ`; sem LASSO; preprints sem veículo | **baixa**; citar se a seção de aditivos com wavelets pedir a lista completa (Vidakovic é coautor do autor) | `resumo` (arXiv) |
| Dalalyan, Ingster & Tsybakov (2014), *Probab. Theory Related Fields* 158, "Statistical inference in compound functional models" (arXiv 1208.6402) | taxas minimax não assintóticas em modelos compostos (inclui aditivos esparsos); K&P citam como comparação | sem coeficientes variáveis | **baixa** | título e resumo; volume e páginas `[VERIFICAR]` |
| Benhaddou, Chokri & Pinschenat (2026), arXiv 2603.08538, "Minimax estimation for varying coefficient model via Laguerre series" | séries de Laguerre e mínimos quadrados em coeficientes variáveis; taxas minimax em Laguerre–Sobolev; normalidade e intervalos | uma moduladora; sem esparsidade; sem wavelets | **baixa** | `resumo` |
| Chokri & Bouzebda (2024), *Comm. Statist. Theory Methods* 53(23) 8376–8411, DOI 10.1080/03610926.2023.2286905 | normalidade assintótica das componentes do modelo parcialmente linear aditivo estimado por wavelets | sem coeficientes variáveis; sem penalização | **baixa** | título e metadados |
| Deshpande, Bai, Balocchi et al. (2024), *Bayesian Anal.*, DOI 10.1214/24-BA1470, VCBART; Franco-Villoria, Ventrucci & Rue (2019), *EJS* 13, DOI 10.1214/19-EJS1653 | coeficientes variáveis em **várias** moduladoras por árvores bayesianas (adaptativo, sem aditividade imposta); revisão unificada dos modelos bayesianos de coeficiente variável | sem wavelets; sem teoria de taxa | **média** para E4: o VCBART é o concorrente não aditivo com várias moduladoras; **baixa** para a teoria | `resumo` (VCBART); título (EJS 2019); metadados `[VERIFICAR]` |

Também visto e descartado como ameaça (uma linha cada): Bai, Boland & Chen
(2019, arXiv 1907.06477, VC bayesiano escalável, sem wavelets); Lee &
Mammen (2016, *EJS* 10, local linear em VC esparso de alta dimensão, sem
wavelets); Cheng, Honda & Zhang (2016, *JASA*, seleção *forward* em VC,
B-splines); Honda (2019, *AISM*, group LASSO desviesado em VC, B-splines);
Petersen, Witten & Simon (2016, *JCGS*, fused lasso additive model: aditivo
adaptativo por variação total, sem `X_ℓ`); Brooks, Zhu & Lu (arXiv
1411.5725, seleção local em VC espacial por adaptive group lasso, local
linear).

## 2. Busca 2: *Statistica Sinica* e EJS desde 2015

Método: API do Crossref, filtro `issn` (SS 1017-0405; EJS 1935-7524) e
`from-pub-date:2015-01-01`, `query.bibliographic` com cada termo, e depois
filtro exato pelas palavras no título. Os anos são os do campo `issued` do
Crossref; a SS atribui a artigos aceitos o ano do volume futuro (por isso
aparecem 2027 e 2029).

| Termo no título | SS (n) | EJS (n) |
|---|---|---|
| varying coefficient | 28 | 7 |
| additive coefficient | 1 (Liu, Tu, Bao & Jiang, "varying-coefficient additive model with functional response", 2027) | 0 |
| functional coefficient | 5 (inclui os dois VC acima) | 0 |
| wavelet | 2 (Penev & Hall 2018, médias erráticas com erro de medida; Karamikabir & Afshari 2021, limiares SURE) | 6 (Autin, Claeskens & Freyermuth 2015; Luo, Qi & Wang 2016; Chau & von Sachs 2016; Chichignoud et al. 2017; Wishart 2019; Loosveldt & Tudor 2025) |
| **interseção VC/aditivo × wavelet** | **0** | **0** |

**Conclusão:** nenhum artigo nas duas revistas, desde 2015, junta coeficiente
variável ou aditivo com wavelets. Na SS o tema "coeficiente variável" está
vivo (28 títulos, cerca de 2,5 por ano) e os vizinhos temáticos são:

| Trabalho (SS) | Por que importa |
|---|---|
| Tu, Park & Wang (2020), "Estimation of functional sparsity in nonparametric varying coefficient models for longitudinal data analysis", DOI 10.5705/ss.202017.0246 | esparsidade **local** (regiões nulas) em VC: é a versão "spline + penalidade funcional" da adaptatividade; comparação natural em E4 |
| Zhong, Zhang & Zhang (2024), "Locally sparse estimator of generalized varying coefficient model for asynchronous longitudinal data", DOI 10.5705/ss.202022.0196 | idem, GLM |
| Zhang, Zhou, He & Wong (2024), "Multivariate varying-coefficient models via tensor decomposition", DOI 10.5705/ss.202022.0103 | várias moduladoras por tensor, sem aditividade |
| Chen & He (2018), "Inference of high-dimensional linear models with time-varying coefficients", DOI 10.5705/ss.202015.0202 | alta dimensão em VC no tempo; cita Klopp & Pensky |
| Yang, Yang & Li (2020), "Feature screening in ultrahigh dimensional generalized varying-coefficient models", DOI 10.5705/ss.202017.0362 | triagem; a linhagem de Wei, Huang & Li continua na revista |
| Park, Oh & Lee (2022), "Lévy adaptive B-spline regression via overcomplete systems", DOI 10.5705/ss.202021.0288 | não é VC, mas é o artigo recente da SS sobre **adaptatividade espacial** com Besov; mostra que o tema entra na revista |

Na EJS os sete títulos de VC são local linear esparso (Lee & Mammen 2016),
bayesiano (Franco-Villoria et al. 2019), inferência quantílica em alta
dimensão (Dai & Kolar 2021), discriminante (Bao & Liu 2022), longitudinal
assíncrono (Liu, Sun & Cao 2023), índice único com interações G×E (Guan,
Zhao & Cui 2023) e espacial (Jin & Mu 2026); os seis de wavelets são teoria
de adaptação em outros modelos. Leitura recomendada para calibrar o
referee provável: Lee & Mammen (2016, EJS), Tu, Park & Wang (2020, SS) e
Zhang et al. (2024, SS).

## 3. Sardy & Ma (2024): o que a teoria cobre

Lido o arXiv 2207.04083v1 (2022-07-08, única versão; 28 páginas). A versão
publicada é *Scand. J. Statist.* 51(1) 89–108 (online 2023-08-31, DOI
10.1111/sjos.12680; `Sardy-Ma-2024`); o resumo publicado (Semantic
Scholar) coincide com o do arXiv a menos de "the best" → "a good FDR–TPR
trade-off"; o texto integral está atrás da verificação anti-robô da Wiley.

**O modelo e o desenho.** Aditivo puro, `µ(x) = c + Σ_ℓ µ_ℓ(x_ℓ)`, erro
gaussiano. Cada `µ_ℓ` é escrita numa matriz `Φ` **`n × n` ortonormal** de
wavelets "isométricas para amostras não equiespaçadas" (Sardy, Percival,
Bruce, Gao & Stuetzle 1999; Kerkyacharian & Picard 2004): uma função de
escala (a constante) e `n − 1` wavelets, aplicada às observações
**ordenadas** pela permutação `P_ℓ`. A matriz de regressão é `W = [P_1'Φ …
P_p'Φ]`, com `pn` coeficientes; não há nível de sieve `J`, não há avaliação
de `ψ` em pontos arbitrários e não há Gram a controlar dentro de cada bloco
(é ortonormal por construção).

**O estimador.** SRAMlet: square-root LASSO, `min ‖y − c1 − Wβ‖_2 + λ‖β‖_1`.
O Teorema 1 mostra que a versão em **grupos** com raiz quadrada é
degenerada em bloco ortonormal (solução de mínimos quadrados se `λ < 1`,
zero se `λ > 1`), o que os leva ao LASSO puro. O Teorema 2 dá a forma
fechada do bloco univariado (limiarização suave com limiar implícito, que
embute um estimador de `σ`, eq. 11); o Teorema 3 dá o algoritmo (relaxação
por blocos de coordenadas, um bloco ortonormal por vez); o Teorema 4 dá a
SURE. `λ` é escolhido pelo **quantile universal threshold** (QUT):
quantil `1 − α` da estatística `λ_0 = ‖W'(y − ȳ1)‖_∞ / ‖y − ȳ1‖_2` sob o
modelo nulo, que é **pivotal** (não depende de `σ` nem de `c`), ou pela
SURE (que precisa de `σ`, MAD).

**O que a teoria cobre.** Só a otimização e a seleção de `λ`: forma da
solução, degenerescência do grupo, convergência do algoritmo, SURE. **Não
há** desigualdade oráculo, taxa de convergência, espaço de Besov, condição
de desenho (compatibilidade, autovalor restrito) nem consistência de
seleção. A evidência é Monte Carlo (`blocks`, `bumps`, `heavisine`,
`doppler` com SNR 3, `U` uniforme, Daubechies "extremal phase" com 4
filtros, `n = 2^10` e `p ∈ {10, 100, 1000}`; depois `(n, p) = (2^j,
2^{j+1})`, `j = 8, …, 11`, onde só os métodos com wavelets rodam) e três
bases reais (`n ∈ {215, 166, 166}`, `p ∈ {100, 235, 235}`), onde o LASSO
linear vence em erro preditivo e o SRAMlet vence em parcimônia. Código em
`github.com/StatisticsL/SRAMlet`.

**O que a condição de desenho de produtos acrescenta (e o que não).** Em
relação a Sardy & Ma, tudo de E1.4 e E1.5 é novo: o desenho deles é uma
concatenação de blocos ortonormais e a interação entre blocos
(`E[ψ(U_m) ψ(U_{m'})']`, `m ≠ m'`) nem é analisada; no WAFC o bloco é
`X_ℓ ψ_{jk}(U_m)`, avaliado em pontos (Daubechies–Lagarias) num sieve de
nível `J_n`, e a Gram não é ortonormal nem dentro do bloco. Mas a
comparação relevante **não é com Sardy & Ma**: é com Klopp & Pensky (2015),
que já têm a fatoração `E[XX'] ⊗ E[ψψ']`, a concentração da Gram restrita e
a taxa adaptativa para `q = 1`, `X ⊥ U`, e com Haris et al. (2018) e Amato
et al. (2022), que já têm compatibilidade e taxa em aditivos esparsos com
wavelets (sem `X_ℓ`). O que resta genuinamente sem dono é a Gram do desenho
**aditivo em `q ≥ 2` moduladoras multiplicado por `X`**, isto é, os blocos
cruzados `E[X_ℓ X_{ℓ'} ψ_{jk}(U_m) ψ_{j'k'}(U_{m'})]`, e o caso `X`
dependente de `U`.

**Duas coisas a aproveitar deles.** (i) A escolha de `λ` por QUT é pivotal
e dispensa `σ` e validação cruzada; a estatística `λ_0` está definida para
qualquer matriz de desenho, logo é simulável no desenho do WAFC e cabe como
alternativa em E2.3 ao lado de BIC/EBIC. (ii) O cenário de simulação deles
(as quatro funções de Donoho–Johnstone com SNR 3, `U` uniforme) é o molde
natural do cenário não homogêneo de E2.4/E4, e usá-lo facilita a comparação
com a literatura de aditivos com wavelets.

## 4. Antoniadis, Gijbels & Lambert-Lacroix (2014)

Correção de atribuição (também feita por L1): o terceiro autor é
**Lambert-Lacroix**, não Verhasselt; os artigos com Verhasselt são os de
seleção com P-splines em aditivos (*Technometrics* 2012) e em coeficientes
variáveis (*JCGS* 2012), citados por este. *Statistical Papers* 55(3)
727–750, recebido 2011-11-30, online 2013-04-26, número de agosto de 2014,
DOI 10.1007/s00362-013-0522-1 (`Antoniadis-Gijbels-LambertLacroix-2014`).

Lidos o resumo e a lista de referências (página da Springer); o texto
integral não está acessível. O que se conclui com segurança:

- **Base: B-splines.** As referências de base são de Boor (1978),
  Nürnberger (1989), Bhatti & Bracken (2006, "the calculation of integrals
  involving B-splines"), Huang, Wu & Zhou (2002) e Qingguo & Longsheng
  (2012, "componentwise B-spline estimation for varying coefficient
  models"); a única referência a wavelets é Donoho & Johnstone (1995),
  genérica. Não há "grouped regularization" em wavelets ali.
- **Penalidades em grupo, em quadro unificado:** group LASSO (Yuan & Lin
  2006), group SCAD (Wang, Chen & Li 2007; Kim, Choi & Oh 2008), group
  bridge (Huang, Ma, Xie & Zhang 2009), penalidades bi-nível e MCP
  (Breheny & Huang 2009, 2011; Zhang 2010) e `ℓ_1`–`ℓ_q` (Liu & Zhang
  2008); algoritmos de gradiente projetado (Birgin et al. 2000; Figueiredo
  et al. 2007; van den Berg et al. 2008).
- **Há teoria:** o resumo promete "the variable and estimation consistency
  of the methods" para cada penalidade, com Bickel, Ritov & Tsybakov (2009),
  Huang & Zhang (2010) e Bach (2008) nas referências; o regime (`p` fixo ou
  crescente) e a forma das taxas ficam `[VERIFICAR]` no texto integral.
- **O modelo:** "allowing the regression coefficients to be functions of
  other variables"; pelo título e pelas referências (Hastie & Tibshirani
  1993; Wei, Huang & Li 2011; Wang, Li & Huang 2008) é o additive
  varying-coefficient model de Xue & Yang, com um grupo por componente.
  `[VERIFICAR]` se a aditividade é nas moduladoras (como no WAFC) e quantas
  moduladoras entram nos exemplos (Boston housing; MACS).

Ameaça: **média** para a contribuição 4 (seleção de "quais moduladoras
afetam quais coeficientes" por grupos já existe com B-splines, com teoria) e
para a variante em grupos de E1.7; **baixa** para as contribuições 1 a 3.
Deve ser citado na introdução ao lado de Xue & Yang (2006) e de Ma, Carroll,
Liang & Xu (2015).

## 5. Veredito por contribuição (`alvo-revista.md` §4)

| # | Contribuição como está em `alvo-revista.md` | Veredito | Razão |
|---|---|---|---|
| 1 | Estimador e teoria: oráculo no desenho `X_ℓ ψ_{jk}(U_m)`, taxas em Besov, corolário de compressibilidade | **novidade parcial** | Klopp & Pensky (2015) já têm, para `q = 1` e `X ⊥ U`, LASSO em blocos com desigualdade oráculo não assintótica, taxa adaptativa em Besov (inclusive `π < 2`) e cota inferior minimax. O que fica: `q ≥ 2` com aditividade nas moduladoras; `X` dependente de `U`; LASSO puro com o corolário weak-`ℓ_τ` do WALL; constantes não penalizadas e identificabilidade no nível da base. A frase-tese precisa dizer "extends the single-index result of Klopp and Pensky to additive multi-modulator coefficients and dependent designs" |
| 2 | Condição de desenho para produtos (E1.4) | **nova no que importa; conhecida na parte (i)** | E1.4 (i), fatoração sob `X ⊥ U`, é a eq. (1.8) e o Lema 1 de K&P; deve ser citada, não provada como nova. O conteúdo novo é a parte (ii) (`λ_min(E[XX' \| U]) ≥ κ_1`) e, sobretudo, o termo cruzado entre moduladoras, `E[ψ_{jk}(U_m) ψ_{j'k'}(U_{m'})]` com `m ≠ m'`, que não é produto de Kronecker e que nenhum dos trabalhos acima trata com `X_ℓ` multiplicando |
| 3 | Evidência numérica de adaptatividade contra splines | **fica nova, com o conjunto de concorrentes ampliado** | ninguém mediu isso em coeficientes aditivos; mas o referee vai pedir, além de `mgcv`/Xue & Yang: o spline adaptativo de Wang, Jiang & Liu (2024, JCGS), o LASSO em blocos de K&P sobre o mesmo desenho (trivial com `grpreg`/`gglasso`), e um não aditivo em várias moduladoras (VCBART). O cenário não homogêneo deve usar as funções de Donoho–Johnstone como em Sardy & Ma |
| 4 | Software com a interface do `wall()` | **nova, de peso baixo** | não há código público para coeficientes variáveis com wavelets e LASSO (K&P não têm; Sardy & Ma têm SRAMlet só para aditivos). Na SS conta como reprodutibilidade, não como contribuição |
| (E1.7) | Seleção de estrutura pela variante em grupos | **não é nova** como ideia | em splines: Antoniadis et al. (2014), Ma et al. (2015), Wei, Huang & Li (2011); em wavelets sem `X_ℓ`: Amato et al. (2022, bi-nível), Yu et al. (2019, SGL). Só é nova como "bi-nível com wavelets no desenho de produtos". Recomendação: manter opcional ou apresentar como variante sem teorema de seleção |

**Veredito geral.** O nicho existe e a busca não achou o artigo do WAFC
escrito por outros; mas o enquadramento "primeira teoria de LASSO com
wavelets em coeficientes variáveis" é falso por causa de Klopp & Pensky
(2015), e o enquadramento "Sardy & Ma com um `X_ℓ` multiplicando" subestima
o que já existe em aditivos com wavelets (Haris et al. 2018; Amato et al.
2022). A linha de risco do `plano-projeto.md` ("Sardy & Ma cobre mais do que
parece → o peso vai para E1.4 e para a aplicação") se confirma, com a
ameaça vindo de Klopp & Pensky e não de Sardy & Ma. O peso do artigo deve
ir para: (a) a estrutura aditiva em várias moduladoras e o desenho
dependente (E1.4 ii), (b) a comparação numérica ampla, (c) a aplicação com
mais de uma moduladora, onde nenhum dos trabalhos acima chega.

## 6. O que muda nos documentos (propostas; decisão do chat principal)

1. `alvo-revista.md` §4: reescrever a frase-tese e a contribuição 1 como
   extensão de Klopp & Pensky; acrescentar à lista "o que o referee vai
   perguntar": "Why not block LASSO as in Klopp and Pensky?" (resposta:
   `glmnet`, LASSO puro com compressibilidade, e a comparação em E2).
2. `plano-projeto.md` E1.4: a parte (i) vira "recordar K&P (1.8) e Lema 1 e
   estender ao desenho aditivo"; o entregável central passa a ser o termo
   cruzado entre moduladoras e a parte (ii).
3. `plano-projeto.md` E2.3: incluir o QUT de Sardy & Ma (Giacobino et al.
   2017) como regra de `λ` sem `σ`, ao lado de BIC/EBIC.
4. `plano-projeto.md` E2.4/E4: concorrentes mínimos: `mgcv`, spline
   adaptativo (Wang, Jiang & Liu 2024), block LASSO de K&P no mesmo desenho,
   VCBART; funções de teste de Donoho–Johnstone no cenário não homogêneo.
5. `proposta-metodo.md` §5: trocar "o mais próximo é Sardy & Ma" por "o
   mais próximo na teoria é Klopp & Pensky (2015); o mais próximo no método
   é Sardy & Ma (2024) e Amato et al. (2022)"; citar Montoril, Morettin &
   Chiann (2018) como precursor do grupo.
