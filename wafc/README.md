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
bases), `glmnet`, `Matrix`; `sparsegl` se E2.2 adotar a variante em grupos;
`mgcv`, `grpreg` para os competidores. Tudo declarado em `R/load.R`.

Carregar tudo numa sessão, da raiz do repositório:

```r
source("wafc/R/load.R")
```

| Arquivo (a criar) | Etapa | O que tem |
|---|---|---|
| `R/load.R` | E2.1 | `library()` das dependências e `source()` de `R/*.R` |
| `R/dgp.R` | E2.1 | cenários: `g_{jk}` suaves (seno, polinômio), não homogêneas (bumps, blocks, heavisine de Donoho–Johnstone reescaladas), nulas; `X` gaussiano ou uniforme, com `X_1 ≡ 1`; `U` uniforme ou beta; erro gaussiano; `simulate_wafc(n, p, q, scenario, seed)` |
| `R/design.R` | E2.1 | `wafc_rescale()`, `wafc_design(X, U, J, j0, family, filter.size, boundary, eps, table)`: blocos `X_j ⊙ wbasis(U_k)`, colunas nomeadas `x{j}:u{k}:l{l}m{m}`, `penalty.factor`, versão esparsa |
| `R/fit.R` | E2.2 | `wafc(x, u, y, J, ..., penalty = c("lasso", "sglasso"), lambda)`; `predict.wafc`, `coef.wafc`, `print.wafc` |
| `R/reconstruct.R` | E2.2 | `wafc_functions(fit, s, grid)`: `ĉ_j`, `ĝ_{jk}(grid)`, `β̂_j(u)` |
| `R/tune.R` | E2.3 | `cv.wafc()` sobre `(J, λ)`; `wafc_bic()`, `wafc_ebic()` |
| `R/competitors.R` | E2.4 | `mgcv::gam` com `s(u_k, by = x_j)`; B-splines + group LASSO (`grpreg`); linear oráculo; oráculo de suporte |
| `R/plot.R` | E3.2 | `plot.wafc` (painel por `(j,k)`; caminho das normas), `plot.cv.wafc` |
| `tests/test-design.R` | E2.1 | com `θ*` na base e sem ruído, `glmnet` com `λ → 0` recupera `θ*`; posto cheio; nomes |
| `tests/test-fit.R`, `test-tune.R` | E2.2, E2.3 | idem por função |
| `scripts/01-smoke.R` | E2.1 | um cenário, `n = 500`: recupera as funções; gráfico de `ĝ_{jk}` contra a verdade |
| `scripts/02-pilot.R` | E2.4 | cenários × `n` × métodos × réplicas; ISE por função, RMSE de predição, suporte, tempo; `Rscript wafc/scripts/02-pilot.R [n_rep] [--quick]` |

Regras (`docs/instrucoes.md`, §6): toda função pública tem teste; o `wall()`
é referência de leitura, não de cópia, e as internas dele não são chamadas;
dependência nova entra em `R/load.R` e em `docs/CONTINUAR.md` na mesma
rodada.
