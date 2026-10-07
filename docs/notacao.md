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
| `J`, `J_n` | nível de resolução do espaço de aproximação | cresce com `n` |
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

## 3. Wavelets e espaço de aproximação

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
| `B_X` | cota de `‖X‖_∞` (D24; era `C_X` em E1.3 e `B_X` em E1.4) |
| `c_U`, `C_U` | cotas da densidade conjunta de `U` no seu suporte (D14, D23) |
| `κ_1`, `κ_2` | cotas dos autovalores de `E(XX' \| U)` |
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

## 8. Emenda de 2026-09-19 (D24): as exigências de estilo da revista

A §5 das instruções da *Statistica Sinica* manda escrever `E(X)` para
esperança, em romano e com parênteses, e não `\mathbb{E}[X]` nem um `E`
caligráfico; o mesmo para `Var(X)` e para a probabilidade. O manuscrito já
nasceu assim, e a notação passa a ser essa **também nas derivações**, para
que os dois não imprimam símbolos diferentes para a mesma coisa.

Na mesma emenda, `C_X` (E1.3) e `B_X` (E1.4) eram dois nomes para a cota de
`‖X‖_∞`: fica **`B_X`**, como na tabela da §4.

A troca no `derivations/macros.tex` está feita: `\E` e `\Prob` imprimem
`E` e `P` em romano. Os `.tex` das derivações recompilaram sem alteração de
fonte, porque `\E` continua sendo `\E`; o que muda é como ele imprime.

(D33, 2026-09-21: nos documentos de trabalho, "sieve" também cede lugar a
"espaço de aproximação"; os símbolos não mudam.)

## 9. Emenda de 2026-10-01 (D40): a seleção de estrutura

A seleção por limiarização (E1.7c, Corolário 8; §3.6 e S7 do manuscrito)
usava quatro símbolos que colidiam com outros já fixados. Ficam:

| Símbolo | Significado | Era | Colidia com |
|---|---|---|---|
| `ν_{ℓm}`, `ν̂_{ℓm}` | norma `L_2[0,1]` da componente `g_{ℓm}` e da estimada; a segunda é a norma euclidiana do bloco `θ̂_{ℓm,·}` | `N_{ℓm}`, `N̂_{ℓm}` | `N` (momentos nulos) e `N_J` |
| `𝒮`, `𝒮̂(t)` | a estrutura, `{(ℓ,m) : ν_{ℓm} > 0}`, e a estimada pelo limiar `t` | `S`, `Ŝ(t)` | `S` (suporte genérico) e `S_0` |
| `ν_min`, `ν_{min,n}` | a separação, `min_{(ℓ,m)∈𝒮} ν_{ℓm}`, fixa ou na sequência triangular | `δ`, `δ_n` | o vetor `δ` da prova do oráculo |
| `Δ̄`, `Δ̄_n` | a cota do erro uniforme por bloco `Δ` | `D`, `D_n` | `D = p + d` (proposição do desenho) |

`Δ`, `ρ_n` e `t_n` não mudam. `ρ̃_n = (log n/n)^{s/(2s+1)}` é a raiz da taxa do
Corolário 5 de E1.6 (Theorem 1 do manuscrito), usada na versão do Corolário 8
nessa resolução (pergunta 32(c)); `ρ̃_n ≤ ρ_n`. No LaTeX: `\nu`, `\widehat{\nu}`,
`\mathcal{S}`, `\widehat{\mathcal{S}}`, `\nu_{\min}`, `\nu_{\min,n}`,
`\bar{\Delta}`.


## 10. Emenda de 2026-10-03 (D47): a penalidade em blocos

Com o block LASSO balanceado como o WAFC (D44), os símbolos da teoria em
blocos (E1.11, E1.12; `derivations/08-blocos.tex` e o adendo do
`06-selecao-limiar.tex`) entram aqui. Dois foram trocados antes de entrar,
por colisão: a razão dos pesos é `ϱ` (`\varrho`), e não `ρ`, que é a taxa
por componente `ρ_n` (§9); a soma dos pesos ao quadrado é `𝒲(·)`
(`\mathcal{W}`), e não `W(·)`, que encosta no espaço `W_J`. As macros
continuam locais no preâmbulo do `08` e do `06`, como no `07`; o manuscrito
as copia em `k = 4`.

| Símbolo | Significado | Nota |
|---|---|---|
| `𝒢`, `G`, `\|𝒢\|`, `\|G\|` | partição das colunas penalizadas em pedaços dentro de cada bloco `(ℓ,m)`; um pedaço; o número de pedaços; o número de colunas de um pedaço | o grupo da penalidade não é o bloco `(ℓ,m)`; a norma que se limiariza continua sendo a do bloco, `ν̂_{ℓm}` |
| `b_n = ⌈log n⌉`, `j*` | tamanho mínimo de um pedaço fino; o nível mais fino com `2^j < b_n` | escalar; não confundir com o viés `𝐛`, em negrito |
| `w_G`, `ϱ`, `ϱ̄` | peso do pedaço; razão `max w_G / min w_G`; sua cota em `n` | forma balanceada com `w_G = sqrt(\|G\|)`: `ϱ² ≤ 3` (Proposição 7) |
| `‖θ‖_{𝒢,w}`, `θ_G` | norma de grupos `Σ_G w_G ‖θ_G‖_2`; subvetor do pedaço | |
| `𝒢_0(θ)`, `𝒲(ℋ)` | pedaços ativos de `θ`; `Σ_{G∈ℋ} w_G²` | com `w_G = 1`, `𝒲(𝒢_0)` é o número de pedaços ativos, o `s_0` de E1.5 |
| `Ψ̃_G` | Gram do pedaço residualizada nas colunas não penalizadas | |
| `R_𝒢(θ; η)` | risco ideal por pedaços, `Σ_G min(‖θ_G‖², η)` | Corolário 9, Lema 15 |
| `λ_{0,G}`, `λ̄_{0,G}`, `λ_w`, `λ_n^𝒢`, `λ^𝒢_{n,1}` | calibração exata do pedaço; sua cota determinística; o nível com pesos; o da teoria; o dos pesos 1 | Lema 14 |
| `𝒯_{𝒢,w}`, `ℰ`, `ℰ_Σ`, `𝒜_n` | eventos da calibração, do ruído, da Gram (`‖Σ̂ − Σ‖ < γ/2`) e do regime | `ℰ_Σ` dá `λ_min` e `λ_max` juntos |
| `A^𝒢_{s,π}`, `C_λ` | constante do Lema 15; constante da ordem de `λ` | `A^𝒢 = 2 + (1 − 2^{1−π(s+1/2)})^{−1}` em `π ≤ 2` |
| `ρ_n^𝒢`, `ρ̃_n^𝒢` | taxas por componente do Teorema 4 (`n^{−s'/(2s'+1)}`) e do Corolário 11 | par de `ρ_n` e `ρ̃_n` da §9 |
| `Q_G`, `q_max` | Gram da variante branca (o `grpreg` ortonormaliza o pedaço); seu maior autovalor | Observação 2 do `08` |
| `θ̂^{(t)}`, `ĝ^{(t)}`, `f̂^{(t)}`, `Δ_{ℓm}`, `𝒮_J`, `Δ̄_n^𝒢` | o ajuste limiarizado em `t`; erro por bloco; estrutura no nível `J`; cota do Corolário 12 | Lema 16, Corolários 12 e 13 |
| `‖θ‖_{𝒢,1}` | norma de grupos com pesos 1, `Σ_G ‖θ_G‖_2` | Corolário 14 (E1.13; D49); macro local `\nGone` |
| `μ_J`, `λ_n^+`, `c_ψ` | fator de inflação da calibração do LASSO acima do teto de concentração, `μ_J = 1 + c_ψ 2^J log d/n`; o `λ` inflado, `μ_{J_n}^{1/2}λ_n`; a constante, `7‖ψ‖²_∞/(12C_U)` | Proposição 4(iii) e (iv) do `05` (E1.14; D52) |
| `μ^𝒢_J`, `λ_n^{𝒢,+}`, `c_μ`, `ξ` | o mesmo nos blocos, cobrindo também `λ_max(Σ̂)`; `ξ = Λ − κ_2C_U` | Corolário 14(iii) e (iv) do `08` (E1.14; D52); `ξ` no lugar do `Δ` da tarefa, que colidia com `Δ` e `Δ̄_n` (D40); macro local `\lamGp` |
| `β̄ = (c̄, θ̄)`, `b̄`, `v̄`, `𝒢̄ = 𝒢_0(θ̄)`, `𝒲̄ = 𝒲(𝒢̄)` | o comparador do oráculo em blocos, o seu viés, o seu erro e os seus pedaços ativos | Teorema 3 do `08`, Theorem 2 do `ms_4` (D52) |
| `M_1`, `Σ̂_{GG}` | o centrador da variante branca; a Gram empírica do pedaço | Observação 2 do `08`, Remark S5.2 (D52) |
| `ℛ_1(θ; t)`, `ℛ_{𝒢,w}(θ; λ)`, `ℛ_{𝒢,1}(θ; t)` | `Σ_a min(t\|θ_a\|, θ_a²)` e os pares em norma de grupos, `Σ_G min(λw_G‖θ_G‖_2, ‖θ_G‖_2²)` e com pesos 1 | Proposição 8 do `05` e Corolário 15 do `08` (E1.15; D54); `\mathcal{R}`, sem uso anterior; não confundir com o risco ideal `R_𝒢(θ; η)` do Lema 15 |
| `ρ̲_n`, `𝒦`, `K = \|𝒦\|`, `c_A` | a taxa da cota inferior, `C_g^{2/(2s+1)}{σ²/(nB_X²C_U)}^{2s/(2s+1)}`, da família de `ρ_n` (§9); o conjunto de blocos do hipercubo e o seu tamanho; a constante de Assouad, `(1 − 2^{−1/2})/2` | Teorema 5 do `09` (E1.9; D57); os demais símbolos do `09` ficam locais. **No supp** (E5g, D58), `K = \|𝒦\|` fica escrito `\|𝒦\|`, porque `K` já é a dimensão de `Ψ(u)` em S4, e o `M = \|𝒞\|` do `09` fica `\|𝒞\|` (colide com `M_ε` e com o `M` da S8) |
| `ϑ_π` | expoente da contagem, `min(1 − 1/π, 1/2)`: nos níveis finos, `‖θ*‖_{𝒢,1}` fica abaixo de `‖θ*‖_1` por `b_n^{−ϑ_π}` | Corolário 14 (E1.13; D49); `\vartheta_{\pi}`, sem uso no `ms_3` nem no `supp_3` |

No LaTeX: `\mathcal{G}`, `b_{n}`, `w_{G}`, `\varrho`, `\bar\varrho`,
`\lVert\theta\rVert_{\mathcal{G},w}`, `\mathcal{G}_{0}`, `\mathcal{W}`,
`\widetilde{\Psi}_{G}`, `R_{\mathcal{G}}`, `\lambda^{\mathcal{G}}_{n}`,
`\rho^{\mathcal{G}}_{n}`.

## 11. Emenda de 2026-10-06 (D68(d), E5i): os símbolos locais do §5.1

Locais à frase do terceiro braço do estudo de simulação (`ms_5`, §5.1),
sem uso fora dela:

| Símbolo | O quê | Onde |
|---|---|---|
| `ζ_ℓ` (`\zeta_\cv`) | o sorteio normal padrão, independente de `U`, que forma `X_ℓ = (1 − ρ²)^{1/2} ζ_ℓ + ρ h(U_{m(ℓ)})` | §5.1 do `ms_5`; trocado de `Z_ℓ`, que colidia com a matriz de desenho `\bZ` e as entradas `Z_{ia}` (D68(d)) |
| `ρ` | a correlação de `X_ℓ` com o seu par, `0.5` no braço | §5.1 do `ms_5` |
| `h(u) = √12 (u − 1/2)` | o escore padronizado de uma moduladora uniforme | §5.1 do `ms_5` |

Colisões conhecidas, aceitas por serem locais e distantes: `ζ` e `ζ_n` são
a altura do hipercubo da cota inferior (S6 do supp); `ρ` sem índice não é a
taxa `ρ_n` (§9) nem `\rhoG`, `\rhoGt` (S7).

