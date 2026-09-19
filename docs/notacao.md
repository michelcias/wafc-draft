# Notação

Fonte única dos símbolos. O `.tex` e os `derivations/` obedecem a este
arquivo; símbolo novo entra aqui primeiro. **Estado: congelada em
2026-09-18** (E1.1), pelas cinco decisões da §6. Alterar exige passar pelo
chat principal e registrar em `ESTADO.md`.

Herda do WALL teórico (`wall-manuscript/manuscript/theo/ms_theo_1.tex`,
Seções 2 e 3) o que não conflita: `ψ_{jk}` para as wavelets, `J_n` para o
nível do sieve, `‖·‖_n` para a norma empírica, `B^s_{π,r}` para Besov. O que
muda é o que os coeficientes funcionais exigem: as duas covariáveis saem de
`j` e `k`, que ficam com a base, e o vetor de coeficientes é `θ`, porque `β`
já é função. O WALL achata o par `(j,k)` num índice único `b_m`; aqui o par
fica visível, porque o bloco `(ℓ,m)` é a unidade do desenho e do grupo em
E2.2.

---

## 1. Dimensões e índices

| Símbolo | Significado | Nota |
|---|---|---|
| `n` | tamanho da amostra | `i = 1, …, n` |
| `p` | número de covariáveis lineares `X_ℓ` | `ℓ = 1, …, p`; `X_1 ≡ 1` permitido |
| `q` | número de covariáveis moduladoras `U_m` | `m = 1, …, q` |
| `j` | nível de resolução da wavelet | `j = j_0, …, J − 1` |
| `k` | translação da wavelet | `k = 0, …, 2^j − 1` |
| `j_0` | nível mais grosso da base | `j_0 = 0` no caso periódico |
| `J`, `J_n` | nível de resolução do sieve | cresce com `n` |
| `N_J` | número de wavelets por par `(ℓ, m)` | `2^J − 2^{j_0}` (`= 2^J − 1` se `j_0 = 0`) |
| `d` | número total de colunas penalizadas | `p q N_J` |
| `s_0` | esparsidade do oráculo | número de coeficientes não nulos de `θ*` |

O par `(ℓ, m)` identifica um bloco: a componente de `β_ℓ` que depende de
`U_m`. O par `(j, k)` identifica uma wavelet dentro do bloco.

**Colisão a vigiar.** `ℓ` é índice de covariável e também a letra usual da
penalidade `ℓ_1` e da bola weak-`ℓ_τ`. Convenção: como índice, `ℓ` aparece
sempre em subscrito de `X`, `β`, `c`, `g` ou `θ`; a penalidade é escrita
`‖θ‖_1`, nunca "`ℓ_1`", e a bola de compressibilidade é escrita
`weak-ℓ_τ` só em texto corrido, nunca colada a um índice. O `m` de
moduladora não colide: o WALL usava `m` para o índice achatado, que aqui não
existe.

## 2. O modelo

```
Y_i = Σ_{ℓ=1}^{p} β_ℓ(U_i) X_{iℓ} + ε_i,            i = 1, …, n
β_ℓ(u) = c_ℓ + Σ_{m=1}^{q} g_{ℓm}(u_m),             ∫_0^1 g_{ℓm}(u) du = 0
```

| Símbolo | Significado | Nota |
|---|---|---|
| `Y` | resposta escalar | |
| `X = (X_1, …, X_p)'` | covariáveis lineares | `X ∈ ℝ^p`, limitadas (hipótese) |
| `U = (U_1, …, U_q)'` | covariáveis moduladoras | `U ∈ [0,1]^q` por hipótese; ver §6.4 |
| `β_ℓ` | coeficiente funcional de `X_ℓ` | `β_ℓ : [0,1]^q → ℝ` |
| `c_ℓ` | nível do coeficiente `ℓ` | não penalizado |
| `g_{ℓm}` | componente aditiva de `β_ℓ` em `U_m` | integral zero em Lebesgue (D22) |
| `ε` | erro | `E[ε | X, U] = 0`, sub-gaussiano de parâmetro `σ` |
| `f(x, u)` | função de regressão `Σ_ℓ β_ℓ(u) x_ℓ` | |

## 3. Wavelets e sieve

| Símbolo | Significado |
|---|---|
| `φ`, `ψ` | função de escala e wavelet-mãe (suporte compacto, `N` momentos nulos) |
| `ψ_{jk}(u) = 2^{j/2} ψ(2^j u − k)` | wavelet periodizada em `[0, 1]` |
| `V_J` | espaço de multirresolução de nível `J` |
| `W_J` | complemento das constantes dentro de `V_J`; dimensão `N_J` |
| `Π_J g` | projeção `L_2` de `g` em `W_J` |
| `θ_{ℓm,jk}` | coeficiente de `ψ_{jk}` em `g_{ℓm}` |
| `θ` | vetor de todos os `θ_{ℓm,jk}`, `θ ∈ ℝ^d` |
| `θ*` | coeficientes do oráculo (projeção da função verdadeira) |
| `Z_i` | linha `i` da matriz de desenho: `(X_{iℓ}) ∪ (X_{iℓ} ψ_{jk}(U_{im}))` |
| `Z` | matriz `n × (p + d)` |
| `Σ = E[Z Z']` | Gram populacional; `Σ̂ = Z'Z/n` a empírica |

A expansão de um bloco é

```
g_{ℓm}(u) = Σ_{j=j_0}^{J−1} Σ_{k=0}^{2^j−1} θ_{ℓm,jk} ψ_{jk}(u).
```

A ordem das colunas de `Z` é: primeiro os `p` termos não penalizados
`X_{iℓ}`, depois os blocos `(ℓ, m)` em ordem lexicográfica, e dentro de cada
bloco as wavelets por `j` crescente e, dentro do nível, por `k` crescente.

## 4. Estimador e penalidade

| Símbolo | Significado |
|---|---|
| `λ`, `λ_n` | parâmetro de penalidade |
| `w_{ℓm,jk}` | pesos da penalidade (1 no LASSO; adaptativos em E1.7) |
| `(ĉ, θ̂)` | estimador LASSO |
| `ĝ_{ℓm}`, `β̂_ℓ`, `f̂` | as funções reconstruídas |
| `‖h‖_n² = (1/n) Σ_i h(X_i, U_i)²` | norma empírica |
| `‖h‖²_{L_2(P)}` | norma populacional |
| `φ_0` | constante de compatibilidade (Bühlmann & van de Geer) |
| `S_0`, `S` | suporte do oráculo; um suporte genérico |

O grupo da variante em grupos (E2.2) é o bloco `(ℓ, m)`: o subvetor
`θ_{ℓm,·} ∈ ℝ^{N_J}`.

## 5. Espaços de funções

| Símbolo | Significado |
|---|---|
| `B^s_{π,r}[0,1]` | espaço de Besov de regularidade `s`, integrabilidade `π`, índice fino `r` |
| `s' = s − (1/π − 1/2)_+` | regularidade efetiva em `L_2` |
| `wℓ_τ` | bola weak-`ℓ_τ` (compressibilidade), `τ < 2` |

## 6. Decisões congeladas em E1.1 (2026-09-18)

Ratificadas pelo autor; entram na tabela de decisões do `ESTADO.md` como D9
a D12.

1. **Índices.** `ψ_{jk}` mantém o padrão da literatura de wavelets: `j` é o
   nível e `k` a translação. As covariáveis é que mudam de letra: a linear é
   `X_ℓ`, `ℓ = 1, …, p`, e a moduladora é `U_m`, `m = 1, …, q`. (A proposta
   do esboço era o contrário, `ψ_{lm}` com a covariável em `j`; o autor
   preferiu não mexer na base.)
2. **Coeficientes em `θ`.** `β_ℓ` é função e `θ` é o vetor estimado, contra o
   `β` do WALL, que lá não tem coeficiente funcional disputando a letra.
3. **Subscrito duplo `θ_{ℓm,jk}`**, com a vírgula separando o bloco da
   wavelet; `ℓm` não é produto. A forma `θ^{(ℓm)}_{jk}` fica descartada.
4. **`U ∈ [0,1]^q` por hipótese populacional**, com densidade limitada longe
   de `0` e de `∞`. A transformação monótona que a prática exige é discutida
   na seção de computação e não entra na teoria de E1.3 a E1.6.
5. **Nome do método:** WAFC (D8, ratificada em 2026-09-19).

## 7. Emenda de 2026-09-19 (D22): centralização de Lebesgue

A §2 escrevia `E[g_{ℓm}(U_m)] = 0`. A restrição oficial passa a ser
`∫_0^1 g_{ℓm}(u) du = 0`, que é o que a base periódica com `j_0 = 0` impõe
sem restrição numérica (Lema 1 de `01-identificabilidade.md`) e o que o
estimador de fato estima. As duas coincidem quando `U_m` é uniforme; sob
D11, que só pede densidade limitada longe de `0` e de `∞`, elas diferem por
um nível, e a passagem é `c_ℓ ↦ c_ℓ + Σ_m E{g_{ℓm}(U_m)}`, que deixa a forma
das componentes intacta.

O que as fontes fazem, conferido nos originais: **Klopp & Pensky (2015)**,
hipóteses (A1) e (A2), tomam a base ortonormal em `L_2([0,1])`, isto é em
Lebesgue, com `φ_0 ≡ 1`, e não impõem centralização nenhuma, porque sem
decomposição aditiva não há competição entre a constante e as componentes;
a densidade entra só por `Φ = E(φφ')`. **Xue & Yang (2006)**, pp. 1425-1426,
seguindo Stone (1985), centralizam na distribuição, `E{α_ls(X_s)} = 0`, e a
equação (3.4) deles recentraliza as componentes ajustadas **empiricamente**,
jogando a diferença no nível. O WALL teórico também adota centralização de
Lebesgue (`ms_theo_1.tex`, §2.1).

A versão de Xue & Yang continua disponível como observação: recentralizar
`ĝ_{ℓm}` subtraindo a média empírica e somá-la a `ĉ_ℓ`. O que impede
adotá-la como hipótese é o alvo do Corolário 4 de E1.6, que mede
`‖ĝ_{ℓm} − g_{ℓm}‖_{L_2}`: com alvo centrado em `P` e estimador centrado em
Lebesgue sobra uma constante que não desaparece sem um passo de
recentralização que não está em prova nenhuma hoje.
