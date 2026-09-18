# prototype

Protótipo do estimador WAFC em R (etapa E2 do plano). Chama o `WaveBased`
instalado para as bases (`wbasis()`, `wtable()`) e o `glmnet` (LASSO) ou o
pacote de sparse group LASSO escolhido em E2.2 para o ajuste. Um harness
comum (mesmos dados, mesma semente) para comparar variantes do estimador e
competidores. Não é pacote; não tem `DESCRIPTION`.

Serve para (i) fixar a matriz de desenho e a reconstrução das funções,
(ii) escolher a variante do estimador (LASSO puro ou em grupos) e a regra de
sintonia `(J, λ)`, (iii) medir contra `mgcv` e splines com group LASSO no
piloto (E2.4) e decidir o go/no-go. Depois vira oráculo de teste de `wafc()`
no `WaveBased`. Nada daqui entra no pacote por cópia.

| Arquivo (a criar) | O que tem |
|---|---|
| `00-dgp.R` | cenários: funções `g_{jk}` suaves (seno, polinômio), não homogêneas (bumps, blocks, heavisine de Donoho–Johnstone reescaladas), nulas; `X` gaussiano ou uniforme, com `X_1 ≡ 1`; `U` uniforme ou beta; erro gaussiano; `simulate_wafc(n, p, q, scenario, seed)` |
| `01-design.R` | `wafc_design(X, U, J, j0, family, filter.size, boundary, eps, table)`: blocos `X_j ⊙ wbasis(U_k)`, colunas nomeadas `x{j}:u{k}:l{l}m{m}`, `penalty.factor`, versão esparsa |
| `02-fit.R` | `wafc_fit(design, y, penalty = c("lasso", "sglasso"), lambda)`, `wafc_reconstruct(fit, s, grid)` devolvendo `ĉ_j`, `ĝ_{jk}(grid)`, `β̂_j(u)` |
| `03-tune.R` | CV sobre `(J, λ)`; BIC e EBIC; `wafc_tune()` |
| `04-competitors.R` | `mgcv::gam` com `s(u_k, by = x_j)`; B-splines + group LASSO (`grpreg`/`gglasso`); regressão linear oráculo; oráculo de suporte |
| `05-smoke.R` | um cenário, `n = 500`: recupera as funções; gráfico de `ĝ_{jk}` contra a verdade |
| `06-pilot.R` | piloto E2.4: cenários × `n` × métodos × réplicas; ISE por função, RMSE de predição, suporte, tempo; `Rscript prototype/06-pilot.R [n_rep] [--quick]` |

Rodar da raiz do repositório: `Rscript prototype/05-smoke.R`. Resultados
(`.rds`) não são versionados; o handoff e o `ESTADO.md` registram os números.
