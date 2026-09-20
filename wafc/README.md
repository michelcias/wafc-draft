# wafc

O código do método WAFC (decisão D4: pasta dedicada dentro do `wafc-draft`,
não o `WaveBased`; **privado por enquanto**, sem pacote público nem cópia
para outro repositório até E3.3). Organizado como pacote sem ser pacote, para que
empacotar (E3.3) seja mover arquivos:

```
wafc/
├── R/          # funções, um arquivo por tema, cabeçalhos roxygen; load.R carrega tudo
├── tests/      # testthat: Rscript -e 'testthat::test_dir("wafc/tests")'
├── scripts/    # numerados: fumaça, piloto, comparações
└── cache/      # resultados intermediários, não versionado
```

Dependências: `WaveBased` instalado (só `wbasis()`, `wtable()`, para as
bases), `glmnet`, `Matrix`; `sparsegl` para a variante em grupos (E2.2);
`mgcv`, `grpreg`, `splines` e `VCBART` para os competidores (E2.4). Tudo
declarado em `R/load.R`.

Carregar tudo numa sessão, da raiz do repositório:

```r
source("wafc/R/load.R")
```

| Arquivo | Etapa | O que tem |
|---|---|---|
| `R/load.R` ✓ | E2.1 | `wafc_depends`, `wafc_attach()`, `wafc_check_suggests()` e o `source()` de `R/*.R`; carrega com `source("wafc/R/load.R")` da raiz |
| `R/dgp.R` ✓ | E2.1 | `wafc_component()` (seno, cosseno, cúbica, bumps, blocks, heavisine de Donoho–Johnstone, todas centradas em Lebesgue e de norma `L_2[0,1]` igual a 1), `wafc_scenario()` (a estrutura ativa: `β_1` aditivo em duas moduladoras, `β_2` em uma, `β_ℓ` constante para `ℓ ≥ 3`), `simulate_wafc(n, p, q, scenario, seed, snr, sigma, x_dist, u_dist, u_rho, cc, amplitude)` e `wafc_beta()` |
| `R/design.R` ✓ | E2.1 | `wafc_rescale()` e `wafc_design(x, u, J, j0, family, filter.size, boundary, rescale, eps, use.table, wavelet.table, sparse, spec)`: blocos `X_ℓ ⊙ wbasis(U_m)` na ordem de D12, `φ_{00}` descartada, colunas nomeadas `x2:u1:psi3.5`, `penalty.factor`, `blocks`, `constant`, versão esparsa, e `spec` para reconstruir o desenho em dados novos |
| `R/fit.R` ✓ | E2.2 | `wafc(x, u, y, J, penalty, lambda, nlambda, lambda.min.ratio, asparse, intercept, standardize, thresh, maxit, design, ...)`: LASSO por `glmnet` e sparse group LASSO por `sparsegl` (grupo = bloco `(ℓ, m)`), `intercept = TRUE` quando há covariável constante e o intercepto somado ao nível dela; `coef.wafc` (intercepto já dobrado no nível), `predict.wafc` (`response`, `beta`, `coefficients`, `nonzero`), `print.wafc` e `wafc_kkt()`, que confere as condições de otimalidade e é o que fixa a escala de `λ` |
| `R/reconstruct.R` ✓ | E2.2 | `wafc_functions(fit, s, grid, n_grid)`: `ĉ_ℓ` e `ĝ_{ℓm}(grid)` pelo `spec` do desenho (grade padrão = amplitude amostral de cada `U_m`); `wafc_blocks(fit, s)`: não nulos e norma `ℓ_2` por bloco; `β̂_ℓ(u)` por `predict(type = "beta")` |
| `R/tune.R` ✓ | E2.3 | `cv.wafc(x, u, y, J, penalty, nfolds, foldid, lambda, nlambda, lambda.min.ratio, type.measure, trace, ...)`: validação cruzada sobre `(J, λ)` com dobras fixas para toda a grade, um desenho e um caminho de `λ` por `J`, `lambda.min` e `lambda.1se` lidos dentro do `J` selecionado, mais `print`, `coef` e `predict`; `wafc_bic()` e `wafc_ebic(gamma)`, com graus de liberdade = não nulos + os `p` níveis; `wafc_J_theory(n, s)` e `wafc_lambda_theory(design, sigma, alpha)`, a regra `J_n = min{J : 2^J ≥ (n/log n)^{1/(2s'+1)}}` do Teorema 2 de E1.6 com o `λ_n` do Corolário 2 de E1.5; `wafc_sigma()` para o `σ` que essa regra exige; e `wafc_tune(rule)`, a entrada única das cinco regras (`cv.min`, `cv.1se`, `bic`, `ebic`, `theory`) |
| `R/competitors.R` ✓ | E2.4 | `wafc_competitor(method, x, u, y, active, ...)`, entrada única dos sete concorrentes de L2d, todos com a mesma interface (`beta(u)`, `g(grid)`, `blocks`, `predict`): `gam` (`mgcv`, `s(u_m, by = x_ℓ)`, `select = TRUE`), `bsgl` (B-splines + group LASSO por bloco, `grpreg`), `klopp` (o block LASSO de Klopp & Pensky no desenho do WAFC: pedaços de `⌈log n⌉` dentro do bloco, sem atravessar nível, e os `c_ℓ` penalizados como bloco próprio), `aspline` (nós adaptativos por bloco, inserção progressiva com o tamanho escolhido por BIC, no espírito de Wang, Jiang & Liu 2024), `vcbart` (`VCBART::VCBART_ind`, sem decomposição aditiva), `linear` (MQO em `x`) e `oracle` (o WAFC restrito aos blocos ativos); mais `wafc_grid_components()`, a métrica comum, e `wafc_lambda_qut()`, a regra QUT de Giacobino et al. (2017), pivotal em `σ` |
| `R/plot.R` | E3.2 | `plot.wafc` (painel por bloco `(ℓ,m)`; caminho das normas), `plot.cv.wafc` |
| `tests/test-design.R` ✓ | E2.1 | com `θ*` na base e sem ruído, `glmnet` com `λ → 0` recupera `θ*`; posto cheio; nomes e ordem de D12 contra a construção de referência de `derivations/check/03-desenho-produtos.R`; esparso igual a denso; tabela igual à avaliação exata; os cenários do `dgp.R` |
| `tests/test-fit.R` ✓ | E2.2 | recuperação exata por `wafc()` nas duas penalidades; KKT ao longo de todo o caminho, com a evidência de que `λ` do objetivo é `λ` do `glmnet` vezes `nvars/npen`; intercepto e nível da covariável constante; `predict`, `coef`, reconstrução contra `θ*`; blocos zerados pela variante em grupos |
| `tests/test-tune.R` ✓ | E2.3 | `wafc_kkt()` nas cinco regras, que é o que confere que o `λ` devolvido é o do objetivo (D17); `wafc_J_theory()` contra a busca exaustiva em quatro `s'` e seis `n`, inclusive `J = 7` em `n = 12800` com `s' = 1/4`, o valor da conferência de E1.6; `wafc_lambda_theory()` contra a fórmula do Corolário 2 e a calibração de `σ̂_max` contra `max|Z|` em `J = 2..5`; BIC e EBIC recalculados à mão; dobras fixas reproduzindo a tabela; recuperação exata pela validação cruzada sem ruído |
| `tests/test-competitors.R` ✓ | E2.4 | recuperação exata por cada concorrente num modelo da forma dele e pelo oráculo de estrutura com `θ*` na base; a interface comum conferida nas coordenadas do modelo (`predict = rowSums(x · beta(u))`, componentes centradas somando a `β_ℓ` menos o nível); os pedaços do block LASSO contra um vetor de grupos montado à mão, e a esparsidade por pedaço; o QUT contra a própria definição (invariância de escala, taxa de morte sob o nulo acima de `1 − α`, e mais conservador que `cv.min`) |
| `scripts/01-smoke.R` ✓ | E2.1 | `Rscript wafc/scripts/01-smoke.R [n] [J] [cenário]` (padrão `500 4 smooth`): desenho, `cv.glmnet`, ISE por bloco em `lambda.min` e `lambda.1se`, erro de predição fora da amostra e `wafc/cache/01-smoke.png` com `ĝ_{ℓm}` contra a verdade |
| `scripts/04-pilot.R` ✓ | E2.4 | `Rscript wafc/scripts/04-pilot.R [n_rep] [partes] [ncores]` (padrão `50 all`), em quatro partes: `competitors` (cenários × `n` × métodos × réplicas; ISE por componente, RMSE fora da amostra, estrutura por bloco, segundos), `lambda` (as cinco regras de E2.3 mais o QUT, contra o oráculo da grade), `margin` (`λ_min(G_eps)` determinístico e o erro contra `eps`, na grade de D23 mais a família `a(L−1)2^{−J}`) e `j1` (se a grade de `cv.wafc()` deve começar em `J = 1`); grava um `.rds` por parte em `WAFC_OUT` |

Regras (`docs/instrucoes.md`, §6): toda função pública tem teste; o `wall()`
é referência de leitura, não de cópia, e as internas dele não são chamadas;
dependência nova entra em `R/load.R` e em `docs/CONTINUAR.md` na mesma
rodada.
