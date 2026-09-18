# Inventário de código existente (E0.5)

**Data:** 2026-09-18. O que o `WaveBased`, o `wall` e o `wall-manuscript` têm,
o que serve ao WAFC e o que não serve. Conclusão primeiro: **o motor inteiro
já existe.** A matriz de desenho do WAFC é a do `wall()` com cada bloco
multiplicado por uma coluna de `X`; o ajuste é o mesmo `glmnet` com família
gaussiana em vez de binomial. O que é novo é (i) o produto por `X_j`, (ii) os
`c_j` não penalizados no lugar do intercepto único, (iii) a variante em
grupos, e (iv) os métodos que reconstroem `β̂_j(u)` e `ĝ_{jk}(u)`.

**Como isso se usa depois de D4 (código em `wafc/`, não no `WaveBased`):** as
funções exportadas do `WaveBased` (`wbasis()`, `wtable()`, `PHI()`, `PSI()`)
são chamadas como dependência; as internas do `wall()` (`.wall_design` e
companhia) **não são chamadas** (são privadas e mudariam sem aviso): servem
de leitura para reimplementar o desenho por blocos em `wafc/R/design.R`
com nomes próprios.

---

## 1. `michelcias/WaveBased` 2.6-0 (`../../WaveBased`)

| Camada | Arquivos | Serve ao WAFC? |
|---|---|---|
| avaliação de `φ`, `ψ` em pontos arbitrários (Daubechies–Lagarias), famílias Daublets/Symmlets/Coiflets, filtros próprios | `R/PHI.R`, `R/PSI.R`, `src/phi_psi_vec.c`, `src/wav_filters_*.c` | **sim, sem mudança** |
| tabelas de interpolação (`wtable()`, `wtable_cache()`), 30 a 500× mais rápido | `R/wtable.R`, `src/wav_table.c`, `src/phi_psi_interp.c` | sim, sem mudança; obrigatório em `n` grande |
| base decomposta `wbasis()` com `boundary = "periodic" / "none" / "interval"` (CDV) | `R/wbasis.R`, `src/wav_basis.c`, `src/cdv_edge.c`, `src/wav_basis_sparse.c` | sim, sem mudança |
| `wall()`: entrada `x`, resposta, `J` por covariável, `j0`, reescalonamento para `[ε, 1−ε]` (`.wall_rescale_pars`, `.wall_eps`), tabela automática (`.wall_table`), desenho por blocos denso ou esparso (`.wall_design`, `.wall_wbasis_sparse`), descarte da função de escala constante quando `j0 = 0` periódico (`drop.phi`), `penalty.factor` (`.wall_penalty`), chamada ao `glmnet` | `R/wall.R` (~900 linhas) | **sim, como referência de leitura**: `wafc/R/design.R` reimplementa `.wall_rescale_pars`, `.wall_eps`, `.wall_design` (com o produto por `X_j`), `.wall_penalty` e `.wall_colnames` com nomes próprios, citando a origem em comentário; `.wall_wbasis_sparse` chama C interno e a versão esparsa em `wafc/` monta o `dgCMatrix` a partir de `wbasis()` denso por bloco |
| `cv.wall()`: CV sobre `(J, λ)`, `share.design`, paralelo | `R/cv.wall.R` | sim, como molde de `wafc/R/tune.R`; muda `type.measure` (`"mse"`, `"mae"`) |
| métodos `predict`, `coef`, `plot` (path, componentes, rede) | `R/wall.R` | sim, como molde; `plot.wafc` ganha painel por `(j, k)` |
| testes do `wall` | `tests/testthat/test-wall.R` | molde dos testes em `wafc/tests/` |
| `glmnet` como dependência (`Imports`) | `DESCRIPTION` | `wafc/` declara `glmnet` em `R/load.R`; sparse group LASSO é dependência própria de `wafc/` (ver §4) |
| CI: `R-CMD-check` em 5 plataformas | `.github/workflows/R-CMD-check.yaml` | verde; nada a fazer |
| `renv.lock` do pacote | | existe; o compêndio fixa o `WaveBased` por commit |

**O que não serve:** `wmixreg`, `bwmixreg`, `bwregime`, `wdensity`,
`bayesthresh` (outros problemas). A camada bayesiana pode inspirar uma
extensão (spike-and-slab nos coeficientes de wavelet do WAFC), fora deste
artigo.

**O `WaveBased` não recebe código deste projeto** (D4). Se E3.3 decidir
empacotar dentro dele, esta tabela diz o que refatorar.

**Padrões a reproduzir do `wall()` sem discutir:** `J` escalar ou por
moduladora; `j0 = 0` periódico como padrão com `drop.phi`; `eps` padrão
`1.9^{-J}`; `use.table = "auto"` por carga de trabalho; `sparse = "auto"`;
`standardize = FALSE` por padrão (base ortonormal); erro informativo em
argumento desconhecido (lição do 2.6-0).

## 2. `wall` (compêndio, `../../wall`)

`R/` funções, `scripts/` orquestração numerada, `config/*.yaml` dados,
`renv.lock` fixando o `WaveBased` por commit, cache por unidade `dataset ×
método` retomável, semente mestra única, `INSTRUCTIONS.md` rastreado. É o
molde do `wafc-studies` (E4.1). O que muda: unidades são `cenário × método ×
réplica` (simulação) em vez de `dataset × método × fold`; métricas são ISE
das funções, RMSE de predição, seleção de suporte e tempo; e a aplicação é
um script próprio.

Lições registradas lá que valem aqui: `renv::restore()` e não
`renv::install()` (o `WaveBased` não está no CRAN); scripts que refazem ajuste
dependem de BLAS e ficam numa máquina só; tabelas agregadas se movem em toda
rodada que acrescente unidade ao cache. Diferença: o compêndio do WAFC fixa
também o código de `wafc/` por commit do `wafc-draft`.

## 3. `wall-manuscript` (`../../wall-manuscript`)

- `manuscript/theo/ms_theo_1.tex`: Seções 2 (base, vetor de parâmetros,
  desenho), 3 (estimador), 5 (hipóteses A1 a A7 e o teorema), 7 (trilha
  lenta em quatro passos) e 8 (trilha rápida, compatibilidade, corolário de
  compressibilidade). **É o molde de E1 e de E5a §3.** As hipóteses de
  densidade das covariáveis (A5), de base (A6) e de nível e regularização
  (A7) transferem quase literalmente; a de compatibilidade (A8) é a que muda
  (desenho de produtos).
- `manuscript/theo/references_theo_1.bib`: fonte das entradas `no wall` de
  `literatura.md`; copiar, não redigitar.
- `docs/instrucoes-revisao.md`: origem das convenções de marcação.

## 4. Ferramentas presentes nesta máquina (2026-09-18)

R 4.6.1; `glmnet`, `grpreg`, `gglasso` (group LASSO), `Matrix`, `bench`,
`testthat`, `devtools`, `roxygen2`, `pkgdown`, `renv`, `wavethresh`,
`foreach`, `doParallel`, `future`, `ggplot2`, `xtable`, `knitr`,
`rmarkdown`, `mgcv` (competidor: `gam(y ~ x1 + s(u1, by = x1) + ...)` ajusta
o modelo de coeficientes aditivos por splines penalizadas), `splines2`,
`posterior`; `latexmk`, `pdflatex`. **Não instalados:** `SGL` (sparse group
LASSO, a referência original) e `targets`. O `WaveBased` 2.6-0 foi instalado
de `../../WaveBased` (`R CMD INSTALL .`) em 2026-09-18.

Na máquina de trabalho do autor (conferida em 2026-09-18), `grpreg`,
`gglasso` e `sparsegl` vieram do repositório do Ubuntu (`r-cran-*`, em
`/usr/lib/R/site-library`) e o `gh` está autenticado como `michelcias`.

Para a variante em grupos (E2.2): `sparsegl` (Liang, Cohen, Sólon Heinsfeld,
Pestilli & McDonald, 2024, *JSS*) é a opção com interface tipo `glmnet` e
matriz esparsa; `SGL` é a referência original (Simon et al. 2013), mais lenta.
Decidir em E2.2 e registrar em `CONTINUAR.md` e em `wafc/R/load.R`.
