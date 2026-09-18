# Notação

Fonte única dos símbolos. O `.tex` e os `derivations/` obedecem a este
arquivo; símbolo novo entra aqui primeiro. **Estado: esboço** (2026-09-18);
congela em E1.1 por ratificação do autor. Herda do WALL teórico
(`wall-manuscript/manuscript/theo/ms_theo_1.tex`, Seções 2 e 3) o que não
conflita: `ψ_{jk}` para as wavelets, `J_n` para o nível do sieve, `‖·‖_n`
para a norma empírica, `B^s_{p,q}` para Besov, e a convenção de coeficientes
em `β`. Muda o que os coeficientes funcionais exigem.

---

## 1. Dimensões e índices

| Símbolo | Significado | Nota |
|---|---|---|
| `n` | tamanho da amostra | `i = 1, …, n` |
| `p` | número de covariáveis lineares `X_j` | `j = 1, …, p`; `X_1 ≡ 1` permitido |
| `q` | número de covariáveis moduladoras `U_k` | `k = 1, …, q` |
| `j_0` | nível mais grosso da base | `j_0 = 0` no caso periódico |
| `J`, `J_n` | nível de resolução do sieve | cresce com `n` |
| `l, m` | nível e translação da wavelet `ψ_{lm}` | `l = j_0, …, J − 1`; `m = 0, …, 2^l − 1` |
| `N_J` | número de wavelets por par `(j, k)` | `2^J − 2^{j_0}` (`= 2^J − 1` se `j_0 = 0`) |
| `d` | número total de colunas penalizadas | `p q N_J` |
| `s_0` | esparsidade do oráculo | número de coeficientes não nulos de `θ*` |

**Conflito a resolver em E1.1:** `j` indexa a covariável linear e também é o
índice de nível usual das wavelets. Proposta: nível de wavelet é `l`, e a
covariável linear é `j`. O WALL usa `j` para o nível; a mudança tem de ser
consciente.

## 2. O modelo

```
Y_i = Σ_{j=1}^{p} β_j(U_i) X_{ij} + ε_i,           i = 1, …, n
β_j(u) = c_j + Σ_{k=1}^{q} g_{jk}(u_k),             E[g_{jk}(U_k)] = 0
```

| Símbolo | Significado | Nota |
|---|---|---|
| `Y` | resposta escalar | |
| `X = (X_1, …, X_p)'` | covariáveis lineares | `X ∈ ℝ^p`, limitadas (hipótese) |
| `U = (U_1, …, U_q)'` | covariáveis moduladoras | `U ∈ [0,1]^q` após reescalonamento |
| `β_j` | coeficiente funcional de `X_j` | `β_j : [0,1]^q → ℝ` |
| `c_j` | nível do coeficiente `j` | não penalizado |
| `g_{jk}` | componente aditiva de `β_j` em `U_k` | integral zero |
| `ε` | erro | `E[ε | X, U] = 0`, sub-gaussiano de parâmetro `σ` |
| `f(x, u)` | função de regressão `Σ_j β_j(u) x_j` | |

## 3. Wavelets e sieve

| Símbolo | Significado |
|---|---|
| `φ`, `ψ` | função de escala e wavelet-mãe (suporte compacto, `N` momentos nulos) |
| `ψ_{lm}(u) = 2^{l/2} ψ(2^l u − m)` | wavelet periodizada em `[0, 1]` |
| `V_J` | espaço de multirresolução de nível `J` |
| `Π_J g` | projeção `L_2` de `g` em `V_J` |
| `θ_{jk,lm}` | coeficiente de `ψ_{lm}` em `g_{jk}` |
| `θ` | vetor de todos os `θ_{jk,lm}`, `θ ∈ ℝ^d` |
| `θ*` | coeficientes do oráculo (projeção da função verdadeira) |
| `Z_i` | linha `i` da matriz de desenho: `(X_{ij}) ∪ (X_{ij} ψ_{lm}(U_{ik}))` |
| `Z` | matriz `n × (p + d)` |
| `Σ = E[Z Z']` | Gram populacional; `Σ̂ = Z'Z/n` a empírica |

## 4. Estimador e penalidade

| Símbolo | Significado |
|---|---|
| `λ`, `λ_n` | parâmetro de penalidade |
| `w_{jk,lm}` | pesos da penalidade (1 no LASSO; adaptativos em E1.7) |
| `(ĉ, θ̂)` | estimador LASSO |
| `ĝ_{jk}`, `β̂_j`, `f̂` | as funções reconstruídas |
| `‖h‖_n² = (1/n) Σ_i h(X_i, U_i)²` | norma empírica |
| `‖h‖²_{L_2(P)}` | norma populacional |
| `φ_0` | constante de compatibilidade (Bühlmann & van de Geer) |
| `S_0`, `S` | suporte do oráculo; um suporte genérico |

## 5. Espaços de funções

| Símbolo | Significado |
|---|---|
| `B^s_{π,r}[0,1]` | espaço de Besov de regularidade `s`, integrabilidade `π`, índice fino `r` |
| `s' = s − (1/π − 1/2)_+` | regularidade efetiva em `L_2` |
| `wℓ_τ` | bola weak-`ℓ_τ` (compressibilidade), `τ < 2` |

## 6. A congelar em E1.1

1. Índice de nível `l` (proposta) contra `j` (WALL). Ver §1.
2. Se o vetor de coeficientes se chama `θ` (proposta, distinto de `β`, que é
   função) ou `β` como no WALL.
3. Como escrever o par `(j, k)`: `θ_{jk,lm}` (proposta) ou `θ^{(j,k)}_{lm}`.
4. Reescalonamento das `U_k`: assumido em `[0, 1]` já na população (mais
   simples), ou tratado como transformação empírica (mais honesto).
5. Nome do método: WAFC (*wavelet additive functional coefficients*),
   provisório.
