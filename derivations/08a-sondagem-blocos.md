# E1.11. Sondagem: a teoria de E1.4 a E1.7c com a penalidade em blocos de Klopp & Pensky

Documento de **sondagem**, não de prova. A pergunta é a 33(b) do
[`../docs/ESTADO.md`](../docs/ESTADO.md): se o block LASSO de Klopp & Pensky
(2015) com os níveis `c_ℓ` livres (o `klopp.free` de E2.5b) passar a ser a
variante do WAFC, o que acontece com a teoria de E1.4 a E1.7c, que é do LASSO
coordenado. A entrega é um veredito com números. Não há numeração global de
resultados aqui: os enunciados são nomeados, não numerados, e as provas são
esboços que dizem o que muda em cada passo de E1.5 e E1.6.

Conferência numérica: [`check/08a-blocos.R`](check/08a-blocos.R) (imprime
`OK` em 2026-09-30, 1 min 43 s nesta máquina; dependências `WaveBased` e
`glmnet`). Notação: `../docs/notacao.md` (congelada em E1.1); os símbolos
novos da §1.4 são proposta e estão listados no handoff.

---

## 0. Veredito

**A teoria transfere, e a taxa melhora por um fator logarítmico. O que não
fecha é o estimador que o código roda, não a teoria.** Um item por linha do
catálogo:

1. **(i) O oráculo sem cone transfere verbatim.** A perfilagem dos níveis
   (Lema 4 de E1.5) não depende da penalidade, e o passo que dispensa o cone
   (autovalor cheio no lugar da compatibilidade) só usa Cauchy-Schwarz. Muda
   um lema, a calibração: com pedaços de peso 1, como na norma (3.1) de K&P,
   `λ ≍ σ sqrt((b_n + log|𝒢|)/n)`, que é da mesma ordem do `λ` do LASSO mas
   penaliza a norma de `b_n ≈ log n` coeficientes. Por coordenada de um
   pedaço cheio, o custo é **0,35 do LASSO em `n = 500`**, e cai como
   `1/log n`.
2. **(ii) O corolário de compressibilidade vira uma cota do risco ideal por
   pedaços**, `Σ_G min(‖θ*_G‖², η)`, que é a forma do Lema 4 de K&P; a
   hipótese de Besov de E1.3 **ainda o implica sem hipótese nova**, pelo mesmo
   argumento de contagem por nível do Lema 9, mais um Hölder dentro do
   pedaço.
3. **(iii) A taxa ganha um fator logarítmico:** de `(log n/n)^{2s/(2s+1)}`
   (Corolário 5) para `n^{−2s/(2s+1)} (log n)^{(2/π−1)_+/(2s+1)}`, isto é, um
   ganho de `(log n)^{2s'/(2s+1)}`. Em `π ≥ 2` o logaritmo some inteiro, e a
   taxa atinge, para `q = 1` com `X ⊥ U`, a cota inferior de K&P.
4. **(iv) O Corolário 8 fica como está.** O Lema 13 não sabe qual é o
   estimador; o Corolário 8 recebe a taxa por componente nova,
   `ρ_n = n^{−s'/(2s'+1)}`, sem logaritmo, e a hipótese de separação
   enfraquece na mesma medida. A estatística continua sendo a norma do bloco
   `(ℓ, m)`, não a do pedaço.
5. **(v) E1.4 não muda.** O custo estimado é de **~5 páginas de derivação e 2 a
   3 dias**, sem matemática nova além de dois lemas curtos.

**O que não fecha:** o `klopp` do piloto penaliza com os pesos padrão do
`grpreg`, `sqrt(|G|)`, sobre colunas ortonormalizadas por pedaço. Com esses
pesos o `λ` da teoria é ditado pelos pedaços unitários do nível 0, e a cota
**volta à ordem do LASSO** (1,26 a 1,30 do custo por coordenada do LASSO na
Parte A). A teoria cobre o estimador com pesos 1 (os de K&P) ou com pesos por
pedaço (os de Lounici et al.), não o padrão do `grpreg`. No melhor `λ` de uma
grade a diferença entre os pesos existe e muda de sinal com a verdade (§7).

**Recomendação.** A objeção teórica da pergunta 33(b) cai: adotar a
penalidade em blocos não custa a teoria e ainda a melhora. Se a variante for
adotada, com pesos 1 e não os do `grpreg`, o que pede uma medição a mais
antes de E2.5 decidir (§9).

---

## 1. O estimador

### 1.1 Klopp & Pensky (2015)

Lido no arXiv 1312.4087v2 (a versão que L2 leu); a numeração abaixo é a do
arXiv, e a do *Annals* não foi conferida.

- **A norma (3.1)** é `‖A‖_block = Σ_j Σ_{l=0}^{M} ‖a_{jl}‖_2`, sem pesos: para
  cada função `f_j`, o bloco 0 é só o coeficiente da constante e os demais são
  `M = L/d` blocos de `d ≈ log n` coeficientes **consecutivos no índice da
  base**, isto é, atravessando níveis (a Observação 1 deles ajusta os tamanhos
  para `d` ou `d + 1` quando `L/d` não é inteiro). O estimador (3.2) minimiza
  `n^{−1}‖Y − Bα‖² + δ‖α‖_block`; na nossa normalização, `δ = 2λ`.
- **A escolha de `δ`** que precede o Teorema 2 é
  `δ̂ = 2(σ C_ω K sqrt(μ) + 1) sqrt((1 + h) φ_max ω_max(1) log p / n)`, com
  `p^μ ≥ 2n` entre as condições do Teorema 2. Com `p` fixo,
  `μ log p ≥ log(2n)`, e `δ̂/2 ≍ σ sqrt(log n/n)`: o `log p` do enunciado
  esconde um `log n`.
- **O Teorema 2** exige `min_k(r_k ∧ r'_k) ≥ 2` e `L + 1 ≥ n^{1/2}`, e a
  prova (pelo Lema 3 deles) usa a (A6), `(s + s_0)(1 + log n) ≤ p`. As duas
  primeiras existem porque o termo de viés é controlado pela norma dual, como
  o de ruído (Lema 2(ii) deles); a (A6) é de alta dimensão e **falha para `p`
  fixo e `n` grande**, que é o regime de E1.6. A taxa do teorema, com `μ log p ≍ log n`, é
  `(σ²/n)^{2r/(2r+1)} (log n)^{(2−ν)_+/(ν(2r+1))}`, que é exatamente a de §4
  com `(r, ν) = (s, π)`.
- **O Lema 4 deles** (risco ideal por blocos sob a classe (A3)) é
  `Σ_l min(‖a_l‖², εd) ≤ C_a^{2/(2r+1)} ε^{2r/(2r+1)} d^{(2−ν)_+/(ν(2r+1))}`,
  e é o que a §3 reproduz sob a hipótese de E1.3.

### 1.2 Lounici, Pontil, van de Geer & Tsybakov (2011)

Lido no arXiv 1007.1771v3 (numeração do arXiv). O group LASSO deles é
`min N^{−1}‖Xβ − y‖² + 2 Σ_j λ_j ‖β^j‖`, com um `λ_j` por grupo, que é a
nossa normalização. O **Lema 3.1** calibra, com ruído gaussiano,
`λ_j ≥ (2σ/sqrt N) sqrt(tr Ψ_j + 2‖Ψ_j‖(2q log M + sqrt(K_j q log M)))`,
`Ψ_j = X_{G_j}'X_{G_j}/N`, com probabilidade `1 − 2M^{1−q}`; o **Teorema 3.1**
dá `N^{−1}‖X(β̂ − β*)‖² ≤ (16/κ²) Σ_{j ∈ J(β*)} λ_j²` sob a condição RE em
grupos, e o **Teorema 3.2**, a versão com viés. O **Teorema 7.1** é a cota
inferior para o LASSO com o `λ` universal: ele paga `σ² log K / N` por
coeficiente selecionado, que é o logaritmo que o group LASSO evita quando o
tamanho do grupo passa de `log M`. A seção 8 estende ao ruído não gaussiano
sob uma hipótese técnica; aqui o ruído sub-gaussiano entra pela desigualdade
de Hsu, Kakade & Zhang (2012), §2.3.

### 1.3 O código

`wafc_kp_groups()` (`wafc/R/competitors.R`) corta cada bloco `(ℓ, m)` em
pedaços de até `⌈log n⌉` translações consecutivas **sem atravessar nível**:
os níveis com `2^j ≤ b` são um pedaço cada, e um nível mais fino é cortado
em pedaços de tamanho `b`, com o resto no fim (`n = 500`, `b = 7`: o nível 3
vira `7 + 1`, e o nível 6, `9 × 7 + 1`). O `wafc_fit_klopp()` chama
`grpreg::cv.grpreg(..., penalty = "grLasso")`, e o `grpreg` 3.6.0 (função
interna `newXG`) **padroniza as colunas, ortonormaliza cada grupo e usa
`group.multiplier = sqrt(|G|)` por padrão**. Na escala original a penalidade
é, portanto, `λ Σ_G sqrt(|G|) ‖B^c_G θ_G‖_n`, com `B^c` as colunas centradas,
e `λ` é escolhido por validação cruzada. Com `penalize.levels = FALSE` os
níveis são o grupo 0, não penalizado: é D3.

### 1.4 O objeto da sondagem

Uma partição `𝒢` das `d = pq N_J` coordenadas penalizadas em pedaços `G`,
cada um dentro de um bloco `(ℓ, m)` e de um nível `j`, com `|G| ≤ b_n`; no
nível `j` há `⌈2^j/b_n⌉` pedaços, e `|𝒢| ≤ pq(2^J/b_n + J)`. Pesos
`w_G > 0` e

```
(ĉ, θ̂) ∈ argmin_{c, θ}  ‖Y − A c − B θ‖_n² + 2λ ‖θ‖_{𝒢,w},      ‖θ‖_{𝒢,w} = Σ_{G ∈ 𝒢} w_G ‖θ_G‖_2 .
```

Três escolhas de peso: **`w_G = 1`** (K&P), **`w_G = sqrt(|G|)`** (o padrão do
`grpreg`) e **`w_G` proporcional a `λ_{0,G}`** da §2.3 (Lounici et al.). A
variante "branca", `Σ_G w_G ‖B̃_G θ_G‖_n`, é a que o `grpreg` penaliza depois
de ortonormalizar; a §2.4 diz o que muda nela. Escrevo `𝒢_0 = {G : θ*_G ≠ 0}`
para os pedaços ativos de um comparador, `W(𝒢_0) = Σ_{G ∈ 𝒢_0} w_G²`,
`Ψ̃_G = B̃_G'B̃_G/n`, e `R_𝒢(θ; η) = Σ_G min(‖θ_G‖_2², η)` para o risco
ideal por pedaços. `b_n` é escalar; o vetor de viés continua `𝐛 = f − f_J`,
em negrito, como em E1.5.

---

## 2. Frente (i): o oráculo sem cone, na norma de blocos

### 2.1 O que passa sem mudança

A prova do Teorema 1 de E1.5 tem seis passos. Cinco passam trocando `‖·‖_1`
por `‖·‖_{𝒢,w}` e `s_0` por `W(𝒢_0)`:

- **Lema 4 (perfilagem e Schur).** A prova só usa que a penalidade não toca
  `c`: `θ̂` minimiza `‖Ỹ − B̃θ‖_n² + 2λ‖θ‖_{𝒢,w}`, a decomposição
  `f̂ − f_J = B̃v + P_A(𝐛 + ε)` vale, `γ̃ ≥ λ_min(Σ̂)` vale, e a unicidade vem da
  convexidade estrita.
- **Passo 1 (básica)** e **Passo 6 (predição, com `σ²p/n`)**: sem mudança.
- **Passo 2 (os dois ruídos).** A dualidade usa a norma dual de
  `‖·‖_{𝒢,w}`, que é `max_G ‖z_G‖_2 / w_G`; o viés continua por Young.
- **Passos 3 e 4.** Desigualdade triangular e decomponibilidade de
  `‖·‖_{𝒢,w}` sobre `𝒢_0`, as duas propriedades de qualquer norma de grupos.
- **Passo 5, o que dispensa o cone.** `‖v_{𝒢_0}‖_{𝒢,w} ≤ sqrt(W(𝒢_0)) ‖v‖_2
  ≤ sqrt(W(𝒢_0)/γ̃) ‖B̃v‖_n`: Cauchy-Schwarz e o Lema 4(iii), sobre o espaço
  inteiro. Lounici et al. precisam da condição RE em grupos no cone
  ponderado (a Hipótese 3.1 deles); aqui ela é dominada pelo autovalor cheio
  de E1.4, exatamente como em E1.5.

O Lema 8 de E1.6 (comparador arbitrário) passa pela mesma lista.

### 2.2 O enunciado

**Teorema 1 em blocos (proposto).** Sob as Hipóteses 1 e 2 de E1.5, com `A`
de posto `p`, no evento `𝒯_𝒢 = {max_G ‖(B̃'ε/n)_G‖_2 / w_G ≤ λ/2}` e com
`v = θ̂ − θ*`:

```
(i)   ½‖B̃v‖_n² + λ‖θ̂‖_{𝒢,w}  ≤  3λ‖θ*‖_{𝒢,w} + 2‖𝐛‖_n²
(ii)  se γ̃ > 0:   ‖B̃v‖_n² ≤ 64 λ² W(𝒢_0)/γ̃ + 16‖𝐛‖_n²
                   ‖v‖_{𝒢,w} ≤ 40 λ W(𝒢_0)/γ̃ + 10‖𝐛‖_n²/λ
                   ‖v‖_2²   ≤ 64 λ² W(𝒢_0)/γ̃² + 16‖𝐛‖_n²/γ̃
(iii) ‖f̂ − f‖_n ≤ ‖B̃v‖_n + 2‖𝐛‖_n + ‖P_A ε‖_n     (sem mudança)
```

As constantes são as de E1.5. Com `w_G = 1`, `W(𝒢_0)` é o número de pedaços
ativos; com `w_G = sqrt(|G|)`, é o número de coordenadas nos pedaços ativos.
Para a taxa lenta (i), `‖θ*‖_{𝒢,1} ≤ ‖θ*‖_1`, de modo que o Lema 7 serve
como está, como cota superior.

### 2.3 O que muda: a calibração

**Lema de calibração em blocos (proposto).** Sob as Hipóteses 1 e 2 de E1.5,
para todo `α ∈ (0, 1)`,

```
P( existe G:  ‖(B̃'ε/n)_G‖_2 > λ_{0,G} )  ≤  α,
λ_{0,G} = σ n^{−1/2} ( sqrt(tr Ψ̃_G) + sqrt(2 ‖Ψ̃_G‖_op log(|𝒢|/α)) ),     Ψ̃_G = B̃_G'B̃_G/n .
```

*Esboço.* Condicione em `(X, U)`. Pela Hipótese 2 e pela independência da
amostra, `E[exp(a'ε) | X, U] ≤ exp(σ²‖a‖²/2)` para todo `a ∈ ℝ^n`, que é a
hipótese do Teorema 1 de Hsu, Kakade & Zhang (2012) com `μ = 0`. Tome a
matriz deles como `A = B̃_G'/n`, de modo que `(B̃'ε/n)_G = Aε` e
`Σ = A'A = B̃_G B̃_G'/n²` tem `tr Σ = tr Ψ̃_G / n` e `‖Σ‖ = ‖Ψ̃_G‖/n`. O
teorema dá `‖Aε‖² ≤ σ²(tr Σ + 2 sqrt(tr(Σ²) t) + 2‖Σ‖t)` com probabilidade
`1 − e^{−t}`.
Como `tr(Σ²) ≤ ‖Σ‖ tr Σ`, o lado direito é no máximo
`σ²(sqrt(tr Σ) + sqrt(2‖Σ‖t))²`. Tome `t = log(|𝒢|/α)` e some sobre os
`|𝒢|` pedaços. ∎

É o Lema 5 de E1.5 com o máximo de coordenadas trocado pelo máximo de
normas de pedaço, e é a versão sub-gaussiana do Lema 3.1 de Lounici et al. O
evento `𝒯_𝒢` tem probabilidade ao menos `1 − α` sempre que
`λ w_G ≥ 2λ_{0,G}` para todo `G`.

**Forma fechada.** `M_A` é projeção, logo `Ψ̃_G ⪯ Σ̂_GG`; daí
`tr Ψ̃_G ≤ |G| σ̂²_max` (Lema 6 de E1.5) e `‖Ψ̃_G‖ ≤ λ_max(Σ̂) ≤ Λ = κ_2 C_U + γ/2`
no evento da Proposição 3(iv) de E1.4 (é o Lema 9(iii) de E1.6). No mesmo
evento do Corolário 2 de E1.5,

```
λ_n^𝒢 = 2σ n^{−1/2} ( B_X sqrt(2 C_U b_n) + sqrt(2Λ log(|𝒢|/α)) )   ≍   σ sqrt((b_n + log|𝒢|)/n)
```

é admissível para `w_G = 1`. É o análogo do `λ_n = 2σ B_X sqrt(2C_U) sqrt(2 log(2d/α)/n)`
do Corolário 2.

### 2.4 Com qual `λ`, e o que cada peso compra

| pesos | `λ` admissível | `λ² W(𝒢_0)` para `s` coordenadas em `k` pedaços |
|---|---|---|
| LASSO (E1.5) | `2σ σ̂_max sqrt(2 log(2d/α)/n)` | `≍ σ² s log d / n` |
| `w_G = 1` (K&P) | `λ_n^𝒢` acima, dominado pelo pedaço maior | `≍ σ² k (b_n + log\|𝒢\|)/n` |
| `w_G = sqrt(\|G\|)` (`grpreg`) | `2σ n^{−1/2}(σ̂_max + sqrt(2Λ log(\|𝒢\|/α)))`, **ditado pelo pedaço unitário do nível 0** | `≍ σ² s log\|𝒢\| / n`: a ordem do LASSO |
| `w_G ∝ λ_{0,G}` (Lounici) | `λ w_G = 2λ_{0,G}` | `≍ σ² (s + k log\|𝒢\|)/n` |

Três leituras:

- **Com pesos 1 e pedaços cheios** (`s = k b_n`), a cota é
  `σ² s (1 + log|𝒢|/b_n)/n`: com `b_n ≍ log n` o logaritmo sai. Com um
  coeficiente por pedaço (`k = s`), ela é `σ² s (b_n + log|𝒢|)/n`, **pior**
  que a do LASSO. É o compromisso que Lounici et al. discutem na seção 7, e
  é o que a Parte B mede.
- **O `δ̂` de K&P é o `λ_n^𝒢` com pesos 1**, a menos de constantes: os dois
  são `≍ σ sqrt(log n / n)` com `p` fixo. A diferença de método é que K&P
  controlam o viés pela norma dual (daí `r* ≥ 2` e `L + 1 ≥ n^{1/2}`), e
  E1.5 o controla por Young, sem essas restrições.
- **Os pesos do `grpreg` desfazem o ganho na teoria.** Com um único `λ` e
  pesos `sqrt(|G|)`, o pedaço de uma coordenada do nível 0, que sempre existe,
  exige `λ` na escala do LASSO, e os pedaços cheios ficam superpenalizados
  por um fator
  `sqrt(b_n)(1 + sqrt(2 log(|𝒢|/α)))/(sqrt(b_n) + sqrt(2 log(|𝒢|/α)))`, 1,9 em
  `n = 500`. A cota resultante é a do LASSO. Isso é sobre a cota com o
  `λ` da teoria; com `λ` por validação cruzada nenhuma das duas variantes tem
  teorema (como o LASSO de E1.5).

**A variante branca.** Na norma `‖B̃_G θ_G‖_n` a estatística do pedaço é
`‖P_G ε‖_2/sqrt(n)`, com `P_G` a projeção no espaço-coluna de `B̃_G`:
qui-quadrado, **pivotal**, e
`λ_{0,G} = σ n^{−1/2}(sqrt|G| + sqrt(2 log(|𝒢|/α)))` sem `σ̂_max` nem `Λ`. O
preço aparece no Passo 5: `‖B̃_G v_G‖_n ≤ ‖Ψ̃_G‖^{1/2}‖v_G‖_2`, e a constante
ganha um fator `Λ̂ = max_G ‖Ψ̃_G‖`. As duas versões ficam na mesma ordem
(Parte A: 0,24 a 0,29 por coordenada na branca, contra 0,33 a 0,42 na
euclidiana, antes do fator `Λ̂` de 1,24 a 1,51).

**Números (Parte A**, `p = 2` com `X_1 ≡ 1`, `q = 2`, Daublets de filtro 8,
`U` uniforme, `σ = 0,5`, `α = 0,05`, 300 réplicas por célula,
`b_n = ⌈log n⌉`**):**

| `n` | `J` | `b_n` | `\|𝒢\|` | cobertura (5 esquemas) | folga LASSO / blocos | `λ` LASSO | `λ_n^𝒢` | custo por coordenada, pedaço cheio | coeficiente isolado | pesos `sqrt(\|G\|)` |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 4 | 5 | 20 | 1,000 | 2,12 / 2,47 | 0,481 | 0,699 | **0,421** | 2,11 | 1,29 |
| 200 | 4 | 6 | 20 | 1,000 | 2,00 / 2,35 | 0,322 | 0,487 | **0,379** | 2,27 | 1,30 |
| 500 | 4 | 7 | 20 | 1,000 | 1,85 / 2,19 | 0,193 | 0,302 | **0,346** | 2,42 | 1,30 |
| 500 | 5 | 7 | 32 | 1,000 | 1,87 / 2,15 | 0,216 | 0,329 | **0,331** | 2,32 | 1,26 |

As três últimas colunas são razões ao LASSO: `(λ_n^𝒢)²/b_max` sobre `λ²`
do LASSO (pedaço cheio), `(λ_n^𝒢)²` sobre `λ²` (um coeficiente por pedaço),
e `λ²` com pesos `sqrt(|G|)` sobre `λ²` do LASSO. A calibração é
conservadora por um fator 2,2 em blocos e 1,9 no LASSO, como E1.5 já tinha
medido (fator 2). **O `λ` de blocos é 1,45 a 1,56 vezes o do LASSO**, mas
penaliza a norma de até `b_n` coeficientes. Por fórmula, com
`σ̂_max = Λ = 1` e `J = ⌈log_2 n / 3⌉`, a razão do pedaço cheio é 0,50, 0,34,
0,26, 0,23 e 0,20 em `n = 10^2` a `10^6`: decai devagar, como `1/log n`.

---

## 3. Frente (ii): o que substitui a compressibilidade

### 3.1 O oráculo na forma de risco ideal

O Lema 8 de E1.6 com o comparador `θ̄ = θ* 1{G ∈ T}`, `T = {G : ‖θ*_G‖² > η}`,
e o Lema 9(iii) (`‖B(θ* − θ̄)‖_n² ≤ Λ‖θ* − θ̄‖_2²`) dão, com `w_G = 1` e no
evento `𝒯_𝒢 ∩ {λ_min(Σ̂) ≥ γ/2}`,

```
‖f̂ − f‖_n²  ≤  120 Λ R_𝒢(θ*; η) + 120 ‖𝐛‖_n² + 3 ‖P_A ε‖_n²,        η = 1,6 λ² / (Λ γ̃),
```

porque `3 · 64 λ² |T|/γ̃ = 120 Λ η |T|` e `3 · 20 · 2Λ Σ_{G ∉ T} ‖θ*_G‖²` é o
resto. É o risco ideal da projeção por pedaços no limiar `η ≍ λ²`, o
análogo em blocos do que o Corolário 5 de E1.6 controla por weak-`ℓ_τ` com
`b_n = 1`. Esta é a forma que substitui o corolário de compressibilidade:
ela cobre `π ≥ 2` e os pedaços grossos numa linha, onde a versão weak-`ℓ_τ`
por pedaços precisa de dois casos.

### 3.2 Besov implica a cota, sem hipótese nova

**Lema do risco ideal por pedaços (proposto).** Sob a hipótese de Besov de
E1.3 por nível, `‖θ*_{ℓm,j·}‖_π ≤ C_g 2^{−j(s+1/2−1/π)}`, com
`π(s + 1/2) > 1` e pedaços dentro de níveis de tamanho até `b`, para todo
`η > 0`:

```
R_𝒢(θ*; η)  ≤  pq { A_𝒢 C_g^τ (η/b)^{2s/(2s+1)} b^{(2/π − 1)_+/(2s+1)}  +  (x_+ + 1) η },       τ = (s + 1/2)^{−1},
A_𝒢 = 2 + (1 − 2^{1 − π(s+1/2)})^{−1}  (π ≤ 2),      2 + (1 − 2^{−2s})^{−1}  (π ≥ 2),
x = log_2( b^{τ/π} C_g^τ η^{−τ/2} )  (π ≤ 2),        log_2( (C_g² b / η)^{1/(2s+1)} )  (π ≥ 2).
```

*Esboço (`π ≤ 2`).* Separe os níveis em `j_1 = ⌈x⌉`. Na cabeça, cada pedaço
custa no máximo `η`, e há no máximo `2^j/b + 1` pedaços no nível `j`: soma
`≤ (2^{j_1}/b + j_1) η`. Na cauda, `min(u, η) ≤ η^{1−π/2} u^{π/2}`, e dentro de
um pedaço `‖·‖_2 ≤ ‖·‖_π`, logo cada nível contribui no máximo
`η^{1−π/2} ‖θ*_{ℓm,j·}‖_π^π ≤ η^{1−π/2} C_g^π 2^{−j(π(s+1/2) − 1)}`, uma série
geométrica. Com `2^x = (b C_g^π η^{−π/2})^{τ/π}` as duas partes são
`b^{τ/π − 1} C_g^τ η^{1−τ/2}`, e `τ/π − τ/2 = (2/π − 1)/(2s + 1)`.
Para `π ≥ 2`, a cauda vai por `‖θ*_{ℓm,j·}‖_2 ≤ C_g 2^{−js}` e o `b` some. ∎

É o Lema 9(i) de E1.6 relido: a mesma contagem por nível, mais um Hölder
dentro do pedaço. **A hipótese de Besov de E1.3 continua bastando**, como no
Lema 9, e é mais fraca que a classe (A3) de K&P: a deles é uma bola `ℓ_ν`
ponderada sobre todos os coeficientes, `B^r_{ν,ν}`, e a nossa é por nível,
`B^s_{π,∞}`, que a contém. O argumento por nível não precisa da soma sobre
níveis, e o expoente sai o mesmo do Lema 4 de K&P, com
`(d, ε, r, ν) = (b, η/b, s, π)`.

### 3.3 A leitura em weak-`ℓ_τ`

Para `π ≤ 2`, a mesma contagem dá, para o número de pedaços com norma acima
de `ε`, `N_𝒢(ε) ≤ pq{(log_2 b + 1) + 2A b^{τ/π − 1}(C_g/ε)^τ}`, com o `A` do
Lema 9 (conta feita, não conferida à parte): há `b^{1 − τ/π}` vezes menos
pedaços grandes que coeficientes grandes, mas cada um custa `η ≍ σ² b_n/n`
em vez de `σ² log d/n`. A forma de risco ideal da §3.1
faz essa conta sozinha, e é a que eu levaria à derivação.

### 3.4 Números (Parte C)

Três formas de sequência que saturam a hipótese nível a nível (módulo
constante no nível; direção gaussiana; e a "espiga", desfavorável aos
pedaços: um coeficiente de altura `sqrt(η)` por pedaço, em tantos pedaços
quanto o orçamento `ℓ_π` do nível permite), 12 pares `(s, π)` com
`s ∈ {0,8; 1,5; 3}` e `π ∈ {1; 1,5; 2; 4}`, `b ∈ {1, …, 32}`,
`η ∈ {10^{−9}, …, 10^{−3}}`, `J = 20`: **a cota vale nas 864 sequências**, e
a maior razão entre o risco e a cota é **0,786**, isto é, a cota está a um
fator 1,3 da pior sequência construída.

O expoente de `b` com `η = b ε` (o `λ²` de blocos é proporcional a `b_n`),
sem os níveis com `2^j < b`, que são o termo `(x_+ + 1)η`:

| `s` | `π` | previsto `(2/π − 1)_+/(2s+1)` | espiga | módulo constante |
|---|---|---|---|---|
| 0,8 | 1 | 0,385 | 0,383 | 0 |
| 1,5 | 1 | 0,250 | 0,250 | 0 |
| 3 | 1 | 0,143 | 0,150 | 0 |
| 0,8 | 1,5 | 0,128 | 0,128 | 0 |
| 1,5 | 1,5 | 0,083 | 0,087 | 0 |
| 3 | 1,5 | 0,048 | 0,063 | 0 |
| qualquer | 2 | 0 | 0,000 | 0 |
| 0,8 a 3 | 4 | 0 | −0,21 a −0,18 | 0 |

O expoente previsto é atingido pela espiga em `π < 2`, e **um sinal
espalhado dentro do nível não paga nada pelos pedaços**. O logaritmo que
sobra em `π < 2` é o preço dos picos isolados, que é justamente o que a
penalidade em blocos não agrupa.

---

## 4. Frente (iii): a taxa

### 4.1 O análogo do Corolário 5

Com `w_G = 1`, `b_n = ⌈log n⌉`, `λ = λ_n^𝒢` e `J_n = ⌈c log_2 n⌉`, a §3.1 e a
§3.2 dão, com `η_n ≍ σ² b_n/n`, isto é, `η_n/b_n ≍ σ²/n`:

```
‖f̂ − f‖_n²  =  O_p( n^{−2s/(2s+1)} (log n)^{(2/π − 1)_+/(2s+1)} ),     s/((2s+1)s') ≤ c < 1,
```

contra `O_p((log n/n)^{2s/(2s+1)})` do Corolário 5. As peças:

- `(λ_n^𝒢)² ≍ σ²(b_n + log|𝒢|)/n ≍ σ² b_n/n`, porque
  `log|𝒢| ≤ c log n + O(1)` e `c < 1`: o pedaço tem de ser **pelo menos do
  tamanho de `log|𝒢|`**, e `⌈log n⌉` é.
- O viés `2^{−2 J_n s'} ≍ n^{−2cs'}` não domina se `c ≥ s/((2s+1)s')`; a
  condição empírica de E1.4 pede `c < 1`. A janela é **a mesma do
  Corolário 5**, porque `(2 − τ)/(4s') = s/((2s+1)s') ≥ τ/2`; o piso
  `s' > s/(2s+1)` não muda.
- O termo `pq(x_+ + 1)η_n` dos pedaços grossos é `O((log n)²/n)`, e
  `σ²p/n` continua de ordem menor.

**O ganho é `(log n)^{2s'/(2s+1)}` sempre** (Parte D, identidade a `1e−16`
em 19 pares `(s, π)`): em `π ≥ 2` é o logaritmo inteiro,
`(log n)^{2s/(2s+1)}`; em `π < 2` é parte dele, e sobra
`(log n)^{(2/π−1)/(2s+1)}`. Com `π = 1` e `s = 1`, isto é o `s' = 1/2` que
D27 declara para o cenário não homogêneo, o ganho é `(log n)^{1/3}`, 1,9 em
`n = 1000`.

### 4.2 O análogo do Teorema 2 e do Corolário 4

Com todos os pedaços até `J_n` ativos (a leitura densa, `s_0 = d_n`),
`(λ_n^𝒢)² |𝒢_{J_n}| ≍ σ² pq 2^{J_n}/n + O(σ²(log n)²/n)`, e
`2^{J_n} ≍ n^{1/(2s'+1)}` dá `O_p(n^{−2s'/(2s'+1)})`, **sem o logaritmo** do
Teorema 2. O Corolário 4 (componentes) acompanha, pelo `‖v‖_2²` do Teorema 1
em blocos.

### 4.3 Diante de K&P e da cota inferior

Para `q = 1` e `X ⊥ U`, a cota inferior do Teorema 1 de K&P é
`(σ²/n)^{2r/(2r+1)}` sem logaritmo: **em `π ≥ 2` a variante em blocos é
ótima a menos de constante**, e em `π < 2` fica a
`(log n)^{(2/π−1)/(2s+1)}` dela, que é o Corolário 1 de K&P lido com `p`
fixo, e é também o expoente da limiarização em blocos de wavelets de Cai
(1999) **[VERIFICAR: não lido]**. O LASSO de E1.6 fica a
`(log n)^{2s/(2s+1)}`, e o Teorema 7.1 de
Lounici et al. diz que esse logaritmo não é artefato da cota superior, ao
menos no modelo multitarefa deles (não se transpôs aquele argumento ao
desenho de produtos). Para `q ≥ 2` a cota inferior continua não existindo
(E1.9).

Há um ponto em que a transferência cobre mais que o Teorema 2 de K&P: a (A6)
deles exige `p ≥ (s + s_0)(1 + log n)`, que falha com `p` fixo, e o controle
do viés pela norma dual exige `r* ≥ 2`. A arquitetura de E1.5 (perfilagem e
Young) não precisa de nenhuma das duas.

### 4.4 O que o ganho vale nos `n` do piloto

Pouco, e é por isso que ele é de taxa e não de tabela. Com `b_n ≈ 6` a 7 e
`log(|𝒢|/α) ≈ 6`, o `λ` de blocos por coordenada é 0,35 do LASSO em `λ²` em
`n = 500` (Parte A), e a razão cai só como `1/log n`. O ganho realizado no
diagnóstico de E2.5a (`klopp.free` contra `wafc.lasso`, `n = 1000`, 10
réplicas), de 3,5% a 8,7% em `rmse_f` e de 7% a 16% em ISE, é da mesma
direção e bem menor, porque o erro nesses `n` é dominado pelo viés e pelos
coeficientes grossos, onde os pedaços são pequenos e não agrupam nada.

---

## 5. Frente (iv): o Corolário 8

**O Lema 13 fica como está.** Ele é determinístico e cego ao estimador, e a
Observação "o que muda se a variante em grupos ganhar taxa" de
`06-selecao-limiar.tex` já dizia que o Corolário 8 vale para qualquer
estimador com taxa por componente, trocando `ρ_n`.

**O Corolário 8 em blocos** é o mesmo enunciado com três trocas, todas
mecânicas:

```
D_n² = 512 (λ_n^𝒢)² |𝒢_{J_n}| / γ²  +  64 𝓑²_{J_n}/(α'γ)  +  2 C_g² (1 − 2^{−2s'})^{−1} 2^{−2 J_n s'},
ρ_n^𝒢 = n^{−s'/(2s'+1)}   (com 2^{J_n} ≍ n^{1/(2s'+1)}),
Hipótese S:  δ_n / ρ_n^𝒢 → ∞ .
```

A primeira parcela é a do Teorema 1(ii) em blocos com `γ̃ ≥ γ/2`, na leitura
densa que o Corolário 8 já usa. A triagem (parte (ii)) continua **sem
hipótese nenhuma**, e a separação exigida pela parte (iii) fica mais fraca
por um fator `(log n)^{s'/(2s'+1)}`. Numericamente, isso é coberto pela cota
de `‖v‖_2²` da Parte B, que não falha em réplica nenhuma (§7), porque
`Δ² ≤ 2‖v‖_2² + ` viés, como na prova do Corolário 8(i).

**A norma que se limiariza continua sendo a do bloco `(ℓ, m)`**,
`N̂_{ℓm} = ‖θ̂_{ℓm,·}‖_2`, que agrega todos os pedaços do bloco. Limiarizar a
norma do pedaço responderia a outra pergunta (em que escala e em que região de
`U_m` o efeito está), e sob Besov as normas dos pedaços decaem com `j`, de
modo que não há separação a supor; é o "suporte dentro do bloco" que
E1.7c já exclui.

**O estimador em blocos continua não selecionando.** Ele zera pedaços, não
blocos `(ℓ, m)`, e zerar todos os pedaços de um bloco nulo pediria uma
condição de irrepresentabilidade (E1.7a), que os pedaços não mudam: a
redução de E1.7a é sobre os blocos `(ℓ, m)`. Continua estimação seguida de
limiar.

---

## 6. Frente (v): E1.4, e o custo

**E1.4 não muda.** A condição de desenho é de autovalor e não sabe qual é a
penalidade. A única quantidade nova que a calibração consome,
`‖Ψ̃_G‖_op`, é controlada pelo mesmo evento da Proposição 3(iv) que E1.5 e
E1.6 já usam (`λ_max(Σ̂) ≤ Λ`, Lema 9(iii) de E1.6). D13 (`X` dependente de
`U`) e o termo cruzado entre moduladoras também não entram: tudo passa pelo
autovalor. A condição empírica ao longo de `J_n` é
`pq 2^{J_n} log(pq 2^{J_n})/n ≍ n^{c−1} log n → 0`, a mesma.

**Custo, se a variante for adotada:**

| arquivo | o que entra | páginas | dias |
|---|---|---|---|
| `04-oraculo.tex` | lema de calibração em blocos (HKZ e união), Teorema 1 em blocos como lista de substituições, a forma de risco ideal pelo Lema 8 | ~2 | 1 |
| `05-taxas.tex` | lema do risco ideal por pedaços sob Besov (dois regimes de `π`), análogos do Teorema 2, do Corolário 4 e do Corolário 5, com a janela | ~2,5 | 1 |
| `06-selecao-limiar.tex` | uma observação: o Corolário 8 com `ρ_n^𝒢` | ~0,3 | 0,25 |
| `03-desenho-produtos.tex` | nada | 0 | 0 |
| `check/` | as Partes B e C estendidas a `J_n` e à taxa, como em `05-taxas.R` | — | 0,5 |

**Total: ~5 páginas e 2 a 3 dias.** Não há matemática nova além dos dois
lemas propostos, e os dois já estão conferidos aqui. No manuscrito (mapa no
cabeçalho do `ms_2.tex`): o Theorem 2 (o Teorema 1) passa à norma de blocos,
o Theorem 1 (o Corolário 5) muda o expoente do logaritmo, o Theorem 3 (o
Teorema 2) perde o logaritmo, o Lemma 3 (o Lema 9) dá lugar ao lema do risco
ideal por pedaços, e o Corollary 2 recebe o `ρ_n` novo; no suplemento, as
Seções S5 a S7. Mais 1 a 2 dias de E5, e a lista de contribuições deixa de
dizer "LASSO puro".

---

## 7. Os números (Parte B)

Estimador de §1.4 com `w_G = 1`, resolvido por FISTA com reinício no problema
perfilado; KKT a menos de `10^{−9}` de `λ` em todas as réplicas, e o mesmo
solver com pedaços unitários reproduz o `glmnet` a `1,9e−07`. `J = 4`
(`d = 60`), `b = 4` (por bloco, pedaços `{1}, {2,3}, {4..7}, {8..11},
{12..15}`, `|𝒢| = 20`), `n ∈ {125, 250, 500}`, 40 réplicas, `λ` da teoria
(`2 max_G λ_{0,G}`). Três verdades: **cheio**, 12 coeficientes `±1` em três
pedaços inteiros; **espalhado**, os mesmos 12 coeficientes, um por pedaço, em
posições distintas nos blocos `(1, m)` e `(2, m)` (neste desenho
`X_1 ψ_{jk}(U_m)` e `X_2 ψ_{jk}(U_m)` têm correlação 0,65, que é o
`κ_1 = 0,25`, e repetir a wavelet misturaria o efeito dos pedaços com o do
mau condicionamento); **seno**, `g_{11} = sin(2πu)` com o viés fora de `W_4`.

| verdade | `n` | violações das três cotas | folga mínima | `γ̃ < λ_min(Σ̂)` | razão erro / `(λ² W(𝒢_0))` | blocos / LASSO, `λ` da teoria |
|---|---|---|---|---|---|---|
| cheio | 125 | 0 | 1547× | 0 | 1,22 | **0,49** |
| cheio | 250 | 0 | 491× | 0 | 1,40 | **0,48** |
| cheio | 500 | 0 | 323× | 0 | 1,38 | **0,48** |
| espalhado | 125 | 0 | 1871× | 0 | 0,94 | 1,31 |
| espalhado | 250 | 0 | 600× | 0 | 1,28 | 1,70 |
| espalhado | 500 | 0 | 316× | 0 | 1,50 | 1,98 |
| seno | 125 a 500 | 0 | ≥ 1001× | 0 | — | 1,16 a 1,27 |

- **As três cotas do Teorema 1(ii) em blocos valem em todas as 360
  réplicas**, com e sem viés, e `γ̃ ≥ λ_min(Σ̂)` também: o cone não volta.
  As folgas são da ordem das de E1.5 (lá, medianas de 187× a 371×): as
  constantes de Bühlmann e van de Geer são grosseiras também em blocos.
- **O erro escala como `λ² W(𝒢_0)`**: a razão varia 1,15 vezes (cheio) e 1,60
  vezes (espalhado) quando `n` varia por 4, contra 1,05 do LASSO em E1.5 com
  `n` variando por 32.
- **No `λ` da teoria o compromisso aparece inteiro:** blocos tem metade do
  erro do LASSO no suporte cheio e o dobro no espalhado. O seno perde 16% a
  27% porque a energia de `sin(2πu)` está nos níveis grossos, onde os pedaços
  têm 1, 2 e 4 coeficientes e não agrupam nada, e o `λ` de blocos é 1,5 vez o
  do LASSO.

**No melhor `λ` de uma grade** (13 valores de `4λ` a `λ/1024`, com `λ` o do
LASSO da teoria, comum a todos os métodos; `n = 500`, 20 réplicas; mediana da
razão do erro de predição ao do LASSO na mesma réplica; o topo da grade nunca
é escolhido):

| verdade | `w_G = 1` | `w_G = sqrt(\|G\|)` | `grpreg` (branca, `sqrt(\|G\|)`) | K&P literal (pedaços atravessando níveis, `w_G = 1`) | mínimos quadrados |
|---|---|---|---|---|---|
| cheio | **0,773** | 0,836 | 0,850 | 1,086 | 2,081 |
| espalhado | 1,629 | 1,596 | 1,591 | 1,665 | 1,847 |
| seno | 0,952 | **0,773** | 0,853 | 0,851 | 3,258 |

- **Os pesos não são neutros na prática, e a ordem muda com a verdade:**
  pesos 1 ganham no suporte de pedaços cheios (0,77 contra 0,85 do padrão do
  `grpreg`), e pesos `sqrt(|G|)` ganham no seno (0,77 contra 0,95), onde
  penalizar menos os pedaços pequenos protege os coeficientes grossos.
- **O corte dos pedaços também não é neutro,** mas aqui a comparação é
  viciada: as verdades cheio e espalhado foram construídas com os pedaços que
  não atravessam nível, de modo que a coluna de K&P literal perde por
  construção no cheio. No seno, que não sabe dos pedaços, ela fica no meio.
- No espalhado, três réplicas (uma em `w_G = 1`, uma em `sqrt(|G|)`, uma em
  K&P literal) têm o melhor `λ` no fundo da grade, isto é, preferem mínimos
  quadrados; o LASSO nunca.

---

## 8. O que a sondagem não cobre

- **Não é prova.** Os dois lemas propostos (§2.3 e §3.2) têm esboço e
  conferência, e o Teorema 1 em blocos é a lista de substituições da §2.1;
  nada disso foi escrito em `.tex` nem revisado linha a linha.
- **O `λ` por validação cruzada**, que é o que o piloto usa, não tem teoria
  aqui, como não tem em E1.5. A afirmação "a teoria cobre o `klopp.free`" só
  vale para `w_G = 1` (ou pesos por pedaço) e `λ` determinístico.
- **O estimador exato do `grpreg`**: centra as colunas em vez de
  residualizar em `A`, padroniza, e usa `sqrt(|G|)`. A centragem é a
  perfilagem na coluna constante, e a métrica branca é a da §2.4; o peso é o
  que a teoria não cobre.
- **Pedaços que atravessam níveis** (a norma (3.1) literal): o esboço da §3.2
  usa pedaços dentro de níveis, que é o que o código faz; com pedaços
  atravessando níveis cada pedaço toca no máximo dois níveis e a contagem
  muda por constante, mas isso não foi escrito.
- **`b_n` diferente de `⌈log n⌉`**: a taxa da §4.1 pede `b_n ≳ log|𝒢|`; com
  `b_n` menor volta parte do logaritmo, com `b_n` maior cresce o preço dos
  picos isolados. O ótimo de ordem é `log n`, como em K&P; a constante não
  foi estudada.
- **A cota inferior para `q ≥ 2`** (E1.9), **`σ` desconhecido**, **`p`
  crescente** e **a base do intervalo**: pelas mesmas razões de E1.5 e E1.6.
- **O efeito sobre a observação de E1.7a** (irrepresentabilidade em grupos
  `(ℓ, m)`): não muda, porque os grupos daquela condição são os blocos, não os
  pedaços; não foi recalculada.

---

## 9. Veredito detalhado

**O que fecha.** (i) O Teorema 1 de E1.5 em norma de blocos, sem cone e com
as mesmas constantes, com `s_0` trocado por `W(𝒢_0)`; a calibração por Hsu,
Kakade & Zhang, com `λ_n^𝒢 ≍ σ sqrt((b_n + log|𝒢|)/n)` para pesos 1, que é o
`δ̂` de K&P a menos de constante e o `λ_j` de Lounici et al. no caso
sub-gaussiano; cobertura 1,000 em 1200 réplicas. (ii) O risco ideal por
pedaços no lugar do weak-`ℓ_τ`, implicado pela hipótese de Besov de E1.3 sem
hipótese nova, com o expoente do Lema 4 de K&P; a cota vale em 864 sequências
e está a 1,3 da pior. (iii) A taxa `n^{−2s/(2s+1)}(log n)^{(2/π−1)_+/(2s+1)}`,
ganho de `(log n)^{2s'/(2s+1)}` sobre o Corolário 5, sem logaritmo em
`π ≥ 2` e na mesma janela de `J_n`. (iv) O Corolário 8 com `ρ_n` sem
logaritmo. (v) E1.4 intacta.

**O que não fecha.** O estimador do piloto: com os pesos `sqrt(|G|)` do
`grpreg`, o `λ` da teoria é ditado pelos pedaços unitários e a cota volta à
do LASSO, de modo que a teoria desta sondagem não se aplica ao `klopp.free`
como ele está. Com pesos 1 se aplica.

**O que é conjectura.** Que o `klopp.free` com `λ` por validação cruzada
realize, em `n` grande, a vantagem da taxa (o mesmo estatuto do LASSO de
E1.5). Que o logaritmo residual `(log n)^{(2/π−1)/(2s+1)}` em `π < 2` seja
necessário para este estimador (K&P não o provam; nem esta sondagem). Que o
ganho de 3,5% a 8,7% em `rmse_f` medido em E2.5a seja este mecanismo: ele tem
a direção prevista, mas nos `n` do piloto a vantagem da cota é de constante
(§4.4).

**Recomendação.**

1. **A pergunta 33(b) perde a objeção teórica.** Adotar a penalidade em
   blocos com os níveis livres não obriga a trocar a arquitetura de E1.5 pela
   de Lounici et al.: a perfilagem e o autovalor cheio de E1.4 continuam
   fazendo o trabalho, e o enunciado principal melhora. O que muda no artigo
   é o posicionamento diante de K&P (D18): de "LASSO puro, com
   compressibilidade" para "o estimador de K&P, estendido a `q ≥ 2`, a `X`
   dependente de `U` e a níveis livres, com a taxa deles e sem as hipóteses
   `r* ≥ 2` e (A6)".
2. **Se adotar, com pesos 1**, e medir antes. O `klopp.free` de E2.5b usa os
   pesos do `grpreg`; a variante que a teoria cobre é
   `group.multiplier = rep(1, número de pedaços)`, um argumento. A Parte B
   mostra que a troca move o erro em 10% a 20% nos dois sentidos, conforme a
   verdade; só as células de E2.5b dizem para que lado nos cenários do
   artigo. Alternativa que acomoda os dois: juntar num pedaço só os níveis
   com `2^j < b_n` de cada bloco (o que a norma (3.1) de K&P faz ao atravessar
   níveis), que tira os pedaços unitários sem mudar os pesos dos cheios; não
   medida aqui.
3. **Não abrir a prova antes da decisão de E2.5.** Ela custa 2 a 3 dias e não
   tem risco conhecido; se a decisão for pelo LASSO, esta sondagem vira uma
   observação no `05-taxas.tex` (a variante em blocos ganha o logaritmo) e
   uma resposta pronta ao referee que perguntar "why not block LASSO?".
4. **O Corolário 8 e a limiarização seguem como estão**, com qualquer das duas
   penalidades.

---

## 10. Referências

- **Klopp, O. and Pensky, M. (2015).** Sparse high-dimensional varying
  coefficient model: nonasymptotic minimax study. *The Annals of
  Statistics* 43(3), 1273--1299. `Klopp-Pensky-2015`, verificada em L1. Lido
  aqui no arXiv 1312.4087v2 (código-fonte da versão reenviada): norma (3.1),
  estimador (3.2), a escolha de `δ̂` e o Teorema 2 com as condições
  `min(r ∧ r') ≥ 2`, `L + 1 ≥ n^{1/2}` e (A6), o Lema 1 (Gram restrita), o
  Lema 2 (termos aleatório e de viés pela norma dual, com Hanson-Wright), o
  Lema 4 (risco ideal por blocos) e o Corolário 1. A numeração do *Annals*
  não foi conferida.
- **Lounici, K., Pontil, M., van de Geer, S. and Tsybakov, A. B. (2011).**
  Oracle inequalities and optimal inference under group sparsity. *The
  Annals of Statistics* 39(4), 2164--2204.
  `Lounici-Pontil-vandeGeer-Tsybakov-2011`, verificada em L3. Lido aqui no
  arXiv 1007.1771v3: o estimador com `λ_j` por grupo, o Lema 3.1, os
  Teoremas 3.1 e 3.2, o Teorema 7.1 (cota inferior para o LASSO) e a seção 8
  (ruído não gaussiano). A numeração do *Annals* não foi conferida.
- **Hsu, D., Kakade, S. M. and Zhang, T. (2012).** A tail inequality for
  quadratic forms of subgaussian random vectors. *Electronic Communications
  in Probability* 17. DOI `10.1214/ECP.v17-2079`. Conferida no Crossref em
  2026-09-30 (título, autores, veículo, volume, ano, DOI); o Crossref e o
  OpenAlex não registram páginas nem número do artigo
  **[VERIFICAR: número do artigo e páginas]**. O enunciado usado (Teorema 1,
  com `μ = 0`) foi lido no arXiv 1110.2842. Não está no `.bib`.
  Chave sugerida: `Hsu-Kakade-Zhang-2012`.
- **Cai, T. T. (1999).** Adaptive wavelet estimation: a block thresholding
  and oracle inequality approach. *The Annals of Statistics* 27(3). DOI
  `10.1214/aos/1018031262`. Conferida no Crossref em 2026-09-30 (título,
  autor, veículo, volume, número, ano, DOI), sem páginas. **Não lida**: é a
  origem clássica de blocos de tamanho `log n` em limiarização de wavelets e
  do mesmo expoente residual em `π < 2`, citada aqui de memória
  **[VERIFICAR: páginas e o enunciado]**. K&P não a citam. Não está no `.bib`.
  Chave sugerida: `Cai-1999`.
- **Bühlmann, P. and van de Geer, S. (2011)** e **Donoho, D. L. and
  Johnstone, I. M. (1998)**: usados só por meio de E1.5 e E1.6.
- **`grpreg` 3.6.0**, instalado nesta máquina: a padronização, a
  ortonormalização por grupo e o `group.multiplier = sqrt(|G|)` padrão foram
  lidos no código da função interna `newXG`.
