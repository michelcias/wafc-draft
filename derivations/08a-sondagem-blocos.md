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

> **Nota de 2026-10-01 (E2.5e):** onde esta seção e a §9 dizem que a
> teoria não cobre os pesos `sqrt(|G|)` do `grpreg`, vale a ressalva da
> §11: com os pedaços balanceados (todo pedaço fino entre `b_n` e
> `2b_n − 1`) ela os cobre, ao preço de uma constante `ρ² ≤ 3` (1,57 a
> 1,67 nos `n` do piloto) no termo de estimação.

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

Lido no arXiv 1312.4087v2 (a versão que L2 leu) e conferido no *Annals*
(L7); a numeração abaixo é a do *Annals*. A norma (3.1) e o estimador (3.2)
estão na p. 1280, a Observação 1 na p. 1281, a (A3) e a (A6) na §2.2
(p. 1278), o Lema 1 na p. 1282 e o Teorema 2 na p. 1283. Os Lemas 2, 3 e 4 não estão no
artigo, e sim no suplemento (DOI `10.1214/15-AOS1309SUPP`, a referência
[21] deles); a prova do Teorema 2 os cita com esses números, e o suplemento
não foi relido.

- **A norma (3.1)** é `‖A‖_block = Σ_j Σ_{l=0}^{M} ‖a_{jl}‖_2`, sem pesos: para
  cada função `f_j`, o bloco 0 é só o coeficiente da constante e os demais são
  `M = L/d` blocos de `d ≈ log n` coeficientes **consecutivos no índice da
  base**, isto é, atravessando níveis (a Observação 1 deles ajusta os tamanhos
  para `d` ou `d + 1` quando `L/d` não é inteiro). O estimador (3.2) minimiza
  `n^{−1}‖Y − Bα‖² + δ‖α‖_block`; na nossa normalização, `δ = 2λ`.
- **A escolha de `δ`** é a eq. (3.14),
  `δ̂ = (C_K K sqrt(μ) + 1) sqrt((1 + h) φ_max ω_max(1) log p / n)`, com
  `C_K` dependendo só da lei do ruído e `K` o parâmetro sub-gaussiano da
  (A4) (`K = σ` no caso gaussiano); o Teorema 2 usa `δ = 2δ̂`, isto é,
  `λ = δ̂` na nossa normalização. A (3.16) põe `p^μ ≥ 2n` entre as condições
  do Teorema 2. Com `p` fixo, `μ log p ≥ log(2n)`, e `λ = δ̂ ≍ σ sqrt(log n/n)`:
  o `log p` do enunciado esconde um `log n`.
- **O Teorema 2** exige `L + 1 = n^ς` com `1/2 ≤ ς < 1` e
  `r* = min_k(r_k ∧ r'_k) ≥ r_0* > (2ς)^{−1}` (eq. 3.13, p. 1283; o par
  `r* ≥ 2`, `L + 1 = n^{1/2}` do arXiv é o caso `ς = 1/2` da Observação 2),
  e a prova (pelo Lema 3 do suplemento) usa a (A6), `(s + s_0)(1 + log n) ≤ p`.
  As duas primeiras controlam o viés: o de truncamento cai de
  `C_a² s n^{−2r*ς}` a `C_a² s/n` (p. 1294), e o termo de viés no ajuste é
  controlado pela norma dual, como o de ruído (eq. (5.25) do Lema 2 do
  suplemento); a (A6) é de alta dimensão e **falha para `p`
  fixo e `n` grande**, que é o regime de E1.6. A taxa do teorema, com `μ log p ≍ log n`, é
  `(σ²/n)^{2r/(2r+1)} (log n)^{(2−ν)_+/(ν(2r+1))}`, que é exatamente a de §4
  com `(r, ν) = (s, π)`.
- **O Lema 4 do suplemento** (risco ideal por blocos sob a classe (A3);
  eq. (5.31), conferido em `refs/klopp2015-supp.pdf` em 2026-10-01) é
  `Σ_l min(‖a_l‖², εd) ≤ C_a^{2/(2r+1)} ε^{2r/(2r+1)} d^{(2−ν)_+/(ν(2r+1))}`,
  e é o que a §3 reproduz sob a hipótese de E1.3. O suplemento rotula a
  condição como "(A4)", mas a que ele escreve, eq. (5.30), é a (A3) do
  artigo (p. 1278), e a (A4) de lá é a do ruído: é rótulo trocado no
  suplemento, não neste arquivo. O lema pede `1 ≤ ν < ∞` e
  `r > min(1/2, 1/ν)`, de modo que a correspondência `(r, ν) = (s, π)` com
  o lema deles só vale com `π ≥ 1`.

### 1.2 Lounici, Pontil, van de Geer & Tsybakov (2011)

Lido no arXiv 1007.1771v3 e conferido no *Annals* (L7): a numeração é a
mesma nas duas versões. O group LASSO deles é
`min N^{−1}‖Xβ − y‖² + 2 Σ_j λ_j ‖β^j‖`, com um `λ_j` por grupo, que é a
nossa normalização. O **Lema 3.1** (eq. 3.1, p. 2173) calibra, com ruído gaussiano,
`λ_j ≥ (2σ/sqrt N) sqrt(tr Ψ_j + 2‖Ψ_j‖(2q log M + sqrt(K_j q log M)))`,
`Ψ_j = X_{G_j}'X_{G_j}/N`, com probabilidade `1 − 2M^{1−q}`; o **Teorema 3.1**
dá `N^{−1}‖X(β̂ − β*)‖² ≤ (16/κ²) Σ_{j ∈ J(β*)} λ_j²` sob a condição RE em
grupos (eq. 3.10, pp. 2176–2177), e o **Teorema 3.2** (p. 2177), a versão
com viés, cuja prova o *Annals* omite e remete ao arXiv. O **Teorema 7.1**
(p. 2187) é a cota
inferior para o LASSO com o `λ` universal: ele paga `σ² log K / N` por
coeficiente selecionado, que é o logaritmo que o group LASSO evita quando o
tamanho do grupo passa de `log M`. A seção 8 (p. 2189) estende ao ruído não
gaussiano, no modelo multitarefa, sob a Hipótese 8.1; aqui o ruído sub-gaussiano entra pela desigualdade
de Hsu, Kakade & Zhang (2012), Teorema 2.1.

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
hipótese (2.1) do Teorema 2.1 de Hsu, Kakade & Zhang (2012) com `μ = 0`.
Tome a matriz deles como `A = B̃_G'/n`, de modo que `(B̃'ε/n)_G = Aε` e
`Σ = A'A = B̃_G B̃_G'/n²` tem `tr Σ = tr Ψ̃_G / n` e `‖Σ‖ = ‖Ψ̃_G‖/n`. O
teorema dá `‖Aε‖² ≤ σ²(tr Σ + 2 sqrt(tr(Σ²) t) + 2‖Σ‖t)` com probabilidade
`1 − e^{−t}`. O enunciado publicado toma `A` quadrada (`n × n`), e esta é
`|G| × n`; não muda nada, porque só `Σ = A'A` entra: `‖Aε‖ = ‖Σ^{1/2}ε‖`, e
`Σ^{1/2}` é quadrada com o mesmo `Σ`.
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
`(log n)^{(2/π−1)/(2s+1)}` dela, que é o que as eqs. (3.19) e (3.20) de
K&P dão com `p` fixo (o enunciado do Corolário 1, eq. (3.21), traz
`2(2−ν)_+/(ν(2r+1))`, com um fator 2 que (3.6), (3.19) e (3.20) não
sustentam: a razão entre (3.19) e (3.6) é `(log n)^{(2−ν)_+/(ν(2r+1))}`
vezes `(log p/log n)^{2r/(2r+1)}`, limitado sob `n^β ≥ p`; o arXiv v2 tinha
outro corolário, com `(log p)^{2r/(2r+1)}`, e não há errata no Crossref;
conferido em 2026-10-02), e é também o expoente da limiarização em blocos de wavelets de Cai
(1999, Teorema 4, eq. 5.4, p. 908), que vale com `α ≥ 1/p` (no nosso
símbolo, `s ≥ 1/π`) e no modelo de sequência, não no desenho de produtos.
O LASSO de E1.6 fica a `(log n)^{2s/(2s+1)}`, e o Teorema 7.1 de
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
  `min(r ∧ r') ≥ 2`, `L + 1 ≥ n^{1/2}` e (A6) (no *Annals*, `L + 1 = n^ς` e
  `r* > (2ς)^{−1}`, §1.1), o Lema 1 (Gram restrita), o
  Lema 2 (termos aleatório e de viés pela norma dual, com Hanson-Wright), o
  Lema 4 (risco ideal por blocos) e o Corolário 1. Numeração conferida no
  *Annals* em L7: a mesma, com `δ̂` na eq. (3.14) e os Lemas 2, 3 e 4 no
  suplemento (`10.1214/15-AOS1309SUPP`), que não foi relido.
- **Lounici, K., Pontil, M., van de Geer, S. and Tsybakov, A. B. (2011).**
  Oracle inequalities and optimal inference under group sparsity. *The
  Annals of Statistics* 39(4), 2164--2204.
  `Lounici-Pontil-vandeGeer-Tsybakov-2011`, verificada em L3. Lido aqui no
  arXiv 1007.1771v3: o estimador com `λ_j` por grupo, o Lema 3.1, os
  Teoremas 3.1 e 3.2, o Teorema 7.1 (cota inferior para o LASSO) e a seção 8
  (ruído não gaussiano). Numeração conferida no *Annals* em L7: a mesma,
  com as páginas da §1.2.
- **Hsu, D., Kakade, S. M. and Zhang, T. (2012).** A tail inequality for
  quadratic forms of subgaussian random vectors. *Electronic Communications
  in Probability* 17, no. 52, 1--6. DOI `10.1214/ECP.v17-2079`.
  `Hsu-Kakade-Zhang-2012`, verificada em L5 (o número do artigo e as
  páginas vêm do cabeçalho do PDF, que o Crossref não guarda). O enunciado
  usado é o **Teorema 2.1**, com `μ = 0` (Observação 2.2); este arquivo o
  leu primeiro no arXiv 1110.2842, onde é o "Theorem 1" e toma `A` de
  `m × n`. O publicado toma `A` quadrada, o que não muda a §2.1: só
  `Σ = A'A` entra.
- **Cai, T. T. (1999).** Adaptive wavelet estimation: a block thresholding
  and oracle inequality approach. *The Annals of Statistics* 27(3),
  898--924. DOI `10.1214/aos/1018031262`. `Cai-1999`, verificada em L5 e
  lida no PDF da editora. É a origem dos blocos de tamanho `log n` em
  limiarização de wavelets (o BlockJS das seções 3 e 4), e o **Teorema 4**
  (eq. 5.4, p. 908) dá, em Besov com `1 ≤ p < 2` e `α ≥ 1/p`, a taxa
  minimax a menos de `(log n)^{(2/p−1)/(1+2α)}`, o mesmo expoente da §4.3;
  é um resultado no modelo de sequência, não no desenho de produtos. K&P
  não a citam.
- **Bühlmann, P. and van de Geer, S. (2011)** e **Donoho, D. L. and
  Johnstone, I. M. (1998)**: usados só por meio de E1.5 e E1.6.
- **Donoho, D. L. and Johnstone, I. M. (1994).** Ideal spatial adaptation by
  wavelet shrinkage. *Biometrika* 81(3), 425--455. DOI
  `10.1093/biomet/81.3.425`. `Donoho-Johnstone-1994`, verificada em L1 e
  lida no PDF da editora em L7. Usada na §12: a VisuShrink é a Definição 2
  (§4.2, p. 445), e os níveis abaixo de `j_0` deixados intocados, em número
  fixo, com o custo `2^{j_0}σ²/n`, vêm da §2.4 (p. 440).
- **`grpreg` 3.6.0**, instalado nesta máquina: a padronização, a
  ortonormalização por grupo e o `group.multiplier = sqrt(|G|)` padrão foram
  lidos no código da função interna `newXG`.

---

## 11. Adendo de E2.5e: pesos de razão limitada

Seção acrescentada por E2.5e (2026-09-30), sem mudar as anteriores. A
pergunta é o item (iv) de E2.5e no `docs/TAREFA.md`: o argumento das §2 a §5,
escrito para pesos 1, aceita pesos `w_G` com `ρ = max_G w_G / min_G w_G`
limitado, e com que constante? A razão de perguntar é E2.5c: os pesos 1 são
a forma do block LASSO que pior prediz, e a forma balanceada do código
(`wafc_kp_groups(balanced = TRUE)`: os níveis com `2^j < b_n` de cada bloco
num pedaço só, e a sobra de cada nível mais fino absorvida no pedaço
anterior do mesmo nível) põe todo pedaço fino entre `b_n` e `2b_n − 1`
colunas, de modo que os pesos `sqrt(|G|)` do `grpreg` ficam a razão
limitada. Conferência: Parte E de [`check/08a-blocos.R`](check/08a-blocos.R).

**Veredito: aceita, e a constante é `ρ²`, só no termo de estimação.** A
calibração não vê os pesos, o Teorema 1 em blocos já estava escrito para
pesos quaisquer, e o custo `λ² W(𝒢_0)` com pesos `w` e o `λ` que a
calibração pede é no máximo `ρ²` vezes o dos pesos 1. A taxa da §4 fica,
com a constante multiplicada por no máximo `ρ²`, desde que `ρ` seja limitado
em `n`. **Os pedaços balanceados com os pesos do `grpreg` têm
`ρ² = max|G|/min|G| ≤ (2b_n − 1)/(b_n − 1) ≤ 3`** (1,67 e 1,57 nos `n` do
piloto) **e passam a ser cobertos**; o `klopp.free` e o `klopp.merged` têm
`ρ²` da ordem de `b_n`, que come o logaritmo que os blocos ganham. É o "volta à ordem do LASSO" da §2.4, dito em termos de `ρ`.

### 11.1 O argumento

Pesos `w_G > 0` quaisquer, `ρ = max_G w_G / min_G w_G`, e `λ_{0,G}` o do
lema de calibração da §2.3.

1. **O Teorema 1 em blocos (§2.2) já é para pesos quaisquer.** Nenhum dos
   seis passos da §2.1 usa `w_G = 1`: o peso entra só na norma dual,
   `max_G ‖z_G‖_2 / w_G`, e em `W(𝒢_0) = Σ_{G ∈ 𝒢_0} w_G²` (Passo 5).
2. **A calibração não depende dos pesos.** O lema da §2.3 controla o evento
   por pedaço `ℰ = {‖(B̃'ε/n)_G‖_2 ≤ λ_{0,G} para todo G}`, com
   `P(ℰ) ≥ 1 − α`, e `ℰ` não sabe de peso. Com

   ```
   λ_w = 2 max_G λ_{0,G} / w_G ,
   ```

   vale `ℰ ⊂ 𝒯_{𝒢,w}` para todo `w`, porque em `ℰ` cada
   `‖(B̃'ε/n)_G‖_2 / w_G ≤ λ_{0,G}/w_G ≤ λ_w/2`. O peso só escolhe qual
   pedaço dita o `λ`: o de maior `λ_{0,G}/w_G`.
3. **O custo.** As cotas consomem `λ_w² W(𝒢_0) = Σ_{G ∈ 𝒢_0} (λ_w w_G)²`, e
   para todo `G`

   ```
   λ_w w_G  =  2 max_{G'} λ_{0,G'} w_G / w_{G'}  ≤  ρ · 2 max_{G'} λ_{0,G'}  =  ρ λ_1 ,
   ```

   com `λ_1 = 2 max_G λ_{0,G}` o `λ` dos pesos 1. Logo
   `λ_w² W(𝒢_0) ≤ ρ² λ_1² |𝒢_0|`.

**O que isso dá em cada frente**, com `λ = λ_w`, no evento `ℰ`:

- **(i) Teorema 1(ii) em blocos:** `‖B̃v‖_n² ≤ 64 ρ² λ_1² |𝒢_0|/γ̃ + 16‖𝐛‖_n²`,
  e o mesmo fator na cota de `‖v‖_2²`. A cota de `‖v‖_{𝒢,w}` fica na norma
  ponderada, que está entre `w_min ‖v‖_{𝒢,1}` e `w_max ‖v‖_{𝒢,1}`.
- **(ii) Risco ideal (§3.1):** com o comparador em
  `T = {G : ‖θ*_G‖² > ρ² η_1}`, `η_1 = 1,6 λ_1²/(Λγ̃)` o limiar dos pesos 1,
  a parcela de estimação é `192 λ_w² W(T)/γ̃ ≤ 120 Λ ρ² η_1 |T|`, e o resto,
  `120 Λ Σ_{G ∉ T} ‖θ*_G‖²`, vem do viés do comparador no Lema 8 de E1.6,
  que não vê a penalidade. Como `min(u, ρ²η) ≤ ρ² min(u, η)` para `ρ ≥ 1`,

  ```
  ‖f̂ − f‖_n²  ≤  120 Λ ρ² R_𝒢(θ*; η_1) + 120 ‖𝐛‖_n² + 3 ‖P_A ε‖_n² .
  ```

- **(iii) Taxa (§4):** `η_n` vira `ρ² η_n`. Com `ρ` limitado em `n`, a ordem
  `n^{−2s/(2s+1)} (log n)^{(2/π−1)_+/(2s+1)}` fica, e a constante cresce no
  máximo `ρ²` (no termo principal, `ρ^{4s/(2s+1)}`, porque a cota da §3.2 é
  homogênea de grau `2s/(2s+1)` em `η`). A condição da §4.1, pedaço ao menos
  do tamanho de `log|𝒢|`, continua valendo, e o pedaço maior entra em
  `λ_1² ≍ σ²(b_max + log|𝒢|)/n`, de modo que `b_max ≤ 2b_n − 1` não muda a
  ordem.
- **(iv) Corolário 8 (§5):** a primeira parcela de `D_n²` ganha o fator
  `ρ²`; `ρ_n^𝒢` não muda de ordem.
- **A variante branca (§2.4)**, a que o `grpreg` resolve: o mesmo argumento,
  com os `λ_{0,G}` pivotais dela; o fator `Λ̂` do Passo 5 não muda.

**`ρ²` é o preço, não folga da prova.** Com pedaços de uma coluna ao lado de
pedaços de `b_n` e pesos `sqrt(|G|)`, `ρ² = b_n`; o `λ_w` é ditado pelo
pedaço unitário, `λ_w ≍ σ n^{−1/2}(σ̂_max + sqrt(2Λ log(|𝒢|/α)))`, e um
pedaço cheio custa `λ_w² b_n ≍ σ² b_n log|𝒢| / n`, contra
`λ_1² ≍ σ²(b_n + log|𝒢|)/n` com pesos 1. Como `b_n ≍ log|𝒢|`, a razão é da
ordem de `b_n`: a cota atinge `ρ²` a menos de constante, e é a conta da
§2.4.

**Um pouco melhor que `ρ²`.** Como `λ_w w_G ≤ ρ_* · 2λ_{0,G} ≤ ρ_* λ_1`, com
`ρ_* = max_G(λ_{0,G}/w_G) / min_G(λ_{0,G}/w_G)`, a constante pode ser
`min(ρ, ρ_*)²`, e com `sqrt(|G|)` nos pedaços balanceados o custo atinge
`ρ_*²` (§11.4). `ρ_* = 1` são os pesos de Lounici et al.
(`w_G ∝ λ_{0,G}`), os que igualam o custo de cada pedaço a `4λ_{0,G}²`. Na
forma fechada da §2.3, `λ_{0,G} ≤ σ n^{−1/2}(a sqrt|G| + c)`, com
`a = B_X sqrt(2C_U)` e `c = sqrt(2Λ log(|𝒢|/α))`, e os pesos `sqrt(|G|)`
têm `ρ_* = (a + c/sqrt(b_min)) / (a + c/sqrt(b_max)) < ρ`. Os pesos 1 têm
`ρ_* = (a sqrt(b_max) + c)/(a sqrt(b_min) + c)`, também maior que 1 quando os
pedaços têm tamanhos diferentes: nenhum dos dois é o de Lounici et al., e
qual fica mais perto depende de `a sqrt|G|` contra `c` (números na §11.4).

### 11.2 As quatro formas do piloto

`ρ² = max|G|/min|G|` entre os pedaços de wavelet de um bloco (os de `c_ℓ`
ficam livres, D3), com `J` a partir de 5 (`b_n = 6`) e 6 (`b_n = 7`); em
`J` menor os valores são menores ou iguais.

| forma | pedaços | pesos | `ρ²`, `n = 250` (`b_n = 6`) | `ρ²`, `n = 500` e `1000` (`b_n = 7`) | limitado em `n` | coberta |
|---|---|---|---|---|---|---|
| `klopp.unit` | dentro do nível, sobra no fim | 1 | 1 | 1 | sim | sim (§2 a §5) |
| `klopp.free` | dentro do nível, sobra no fim | `sqrt(\|G\|)` | 6 | 7 | não, `= b_n` | não |
| `klopp.merged` | níveis grossos juntos, sobra no fim | `sqrt(\|G\|)` | 3,5 | 7 | não | não |
| `klopp.balanced` | níveis grossos juntos, sobra absorvida | `sqrt(\|G\|)` | 1,67 | 1,57 | sim, `≤ (2b_n − 1)/(b_n − 1) ≤ 3` | **sim** |

- **`klopp.free`:** o nível 0 é um pedaço de uma coluna, e `ρ² = b_n`.
- **`klopp.merged`:** tira o pedaço unitário do nível 0, mas a sobra do fim
  de cada nível, `2^j mod b_n`, continua. Com `b_n` ímpar ela vale 1 sempre
  que `j` é múltiplo da ordem de 2 módulo `b_n` (em `b_n = 7`, nos níveis 3
  e 6), e então `ρ²` é o tamanho do pedaço grosso, até `2b_n − 3`.
- **`klopp.balanced`:** o pedaço fino maior tem até `2b_n − 1` colunas e o
  grosso, `2^{j*+1} − 1 ≥ b_n − 1`, com `j*` o nível mais fino com
  `2^{j*} < b_n`. Fica abaixo de `b_n` só com `b_n` potência de 2 (`b_n − 1`
  colunas) e no bloco sem nível fino (um pedaço só, `ρ = 1`). A cota
  `(2b_n − 1)/(b_n − 1)` é 3 em `b_n = 2`, 2,2 em `b_n = 6` e tende a 2.

### 11.3 A §3.2 com os pedaços balanceados

O lema do risco ideal por pedaços vale como está, com `b = b_n`, o tamanho
mínimo de um pedaço fino. A cabeça usa só que o nível `j` tem no máximo
`2^j/b_n` pedaços (agora `⌊2^j/b_n⌋`, menos que antes); a cauda usa só o
Hölder dentro de um pedaço contido num nível, e os pedaços finos estão; e o
pedaço grosso custa no máximo `η` por bloco, o que pode substituir o termo
`(x_+ + 1)η`. O tamanho máximo, até `2b_n − 1`, não entra no lema; entra só
em `λ_1`, pela §11.1(iii).

### 11.4 Números (Parte E)

Mesmo desenho da Parte A (`p = 2` com `X_1 ≡ 1`, `q = 2`, Daublets de
filtro 8, `U` uniforme, `σ = 0,5`, `α = 0,05`, `b_n = ⌈log n⌉`), com os
pedaços balanceados; o script inteiro imprime `OK` em 65 a 69 s nesta máquina.

**Calibração e custo** (300 réplicas por linha; `ρ²` dos pesos `sqrt(|G|)`
é `max|G|/min|G|`; os máximos são sobre as réplicas; o custo por coordenada
é `λ_w²` sobre o `λ²` do LASSO, a coluna que na Parte A dava 1,26 a 1,30
com os pedaços de E1.11):

| `n` | `J` | `b_n` | tamanhos | cobertura, 6 esquemas | `ρ²` | `max_G (λ_w w_G)²/λ_1²` | `ρ_*²`, `sqrt(\|G\|)` / pesos 1 | custo por coordenada, `sqrt(\|G\|)` | pesos aleatórios: `ρ²` / `max_G (λ_w w_G)²/λ_1²` |
|---|---|---|---|---|---|---|---|---|---|
| 50 | 4 | 4 | 3 a 4 | 1,000 | 1,33 | 1,21 | 1,21 / 1,11 | **0,57** | 2,72 / 2,72 |
| 100 | 5 | 5 | 5 a 8 | ≥ 0,997 | 1,60 | 1,33 | 1,33 / 1,23 | **0,39** | 3,89 / 3,89 |
| 200 | 5 | 6 | 6 a 10 | 1,000 | 1,67 | 1,34 | 1,34 / 1,26 | **0,34** | 2,93 / 2,64 |
| 500 | 6 | 7 | 7 a 11 | 1,000 | 1,57 | 1,30 | 1,30 / 1,23 | **0,31** | 3,78 / 3,24 |

- **A calibração não vê os pesos:** os seis esquemas (o evento por pedaço,
  sem peso; pesos 1; `sqrt(|G|)` com `λ_{0,G}` exato e na forma fechada;
  pesos aleatórios de razão até 2; e `sqrt(|G|)` na variante branca) têm
  cobertura de 0,997 a 1,000, e o evento por pedaço está contido nos quatro
  eventos ponderados euclidianos em todas as 1 200 réplicas.
- **O custo fica abaixo de `ρ²`, e atinge `ρ_*²`:** com `sqrt(|G|)`,
  `max_G (λ_w w_G)²/λ_1²` é 1,21 a 1,34, contra `ρ²` de 1,33 a 1,67, e é
  igual a `ρ_*²` réplica a réplica, porque o pedaço maior é ao mesmo tempo o
  de maior `λ_{0,G}` e o de menor `λ_{0,G}/w_G`. Com os pesos aleatórios a
  cota `ρ²` vale em toda réplica.
- **Com os pedaços balanceados, os pesos do `grpreg` recuperam o ganho dos
  blocos:** o custo por coordenada cai de 1,26 a 1,30 do LASSO (Parte A,
  pedaços de E1.11) para 0,31 a 0,57; com pesos 1 ele é 0,24 a 0,47 no
  pedaço maior e 0,38 a 0,63 no menor.
- **Nestes `n`, os pesos 1 ficam mais perto dos de Lounici et al. que os
  `sqrt(|G|)`** (`ρ_*²` de 1,11 a 1,26 contra 1,21 a 1,34). A vantagem de
  predição dos `sqrt(|G|)` medida em E2.5c, com `λ` por validação cruzada,
  não é explicada por esta cota.

**O Teorema 1 em blocos com `sqrt(|G|)`** (`J = 5`, `b = 6`, pedaços
`{7, 8, 6, 10}` por bloco, `ρ² = 1,67`; duas verdades, dois pedaços cheios
de 6 e 8 coeficientes, e `g_{11} = sin(2πu)` com viés fora de `W_5`;
`n ∈ {250, 375, 500}`, 40 réplicas; em `n = 125` o desenho tem 126 colunas
e `γ̃ < 0`, fora da hipótese): **as três cotas valem nos 240 ajustes**, com
folga mínima de 598 vezes, `γ̃ ≥ λ_min(Σ̂)` em todos e KKT a menos de
`10^{−9}`. A razão realizada `λ_w² W(𝒢_0)/(λ_1² |𝒢_0|)` é 0,93 (os pedaços
ativos têm 6 e 8 colunas, menos que o maior, de 10) e 1,03 (seno), bem
abaixo de `ρ²`.

**O risco ideal (§11.3)**, nas três formas de sequência da Parte C, com os
pedaços balanceados e `b ∈ {2, 3, 4, 6, 7, 12, 32}` (1 008 sequências): a
cota da §3.2 vale em todas, com razão máxima 0,783 (era 0,786 na Parte C), e
continua valendo com o termo `(x_+ + 1)η` trocado por `η` (razão máxima
0,794).

### 11.5 O que isto não cobre

- **O `λ` por validação cruzada**, que é o que o `klopp.balanced` do piloto
  usa. Como em toda a sondagem, a teoria é para `λ` determinístico:
  "coberta", na tabela da §11.2, quer dizer que o estimador com os mesmos
  pedaços e pesos e `λ = λ_w` tem a taxa da §4, a menos de `ρ²`.
- **O `grpreg` exato** padroniza as colunas e ortonormaliza cada pedaço
  antes de pesar por `sqrt(posto)`: é a variante branca, coberta pelo mesmo
  argumento, com `ρ² = max|G|/min|G|` quando todo pedaço tem posto cheio.
- **Pesos que dependem dos dados** (por exemplo `w_G` estimados de `Ψ̃_G`):
  o argumento pede `ρ` determinístico, ou limitado num evento de
  probabilidade alta; não foi escrito.
- **Não é prova**, como o resto do documento: a §11.1 é álgebra de três
  linhas sobre as §2.2 e §2.3, e a §11.3 é a releitura da §3.2; as duas
  estão conferidas na Parte E.
- **Se a forma balanceada prediz como o `klopp.free`** é a medição de E2.5e,
  relatada no `docs/handoff-E2.5e.md`, não nesta seção.

## 12. Adendo de E2.5f: os níveis grossos livres

Seção acrescentada por E2.5f (2026-10-01), sem mudar as anteriores. A
pergunta é o item (iv) de E2.5f no `docs/TAREFA.md` (pergunta 37 do
`docs/ESTADO.md`): se os níveis com `2^j < b_n` de cada bloco, que a forma
balanceada da §11 penaliza juntos num pedaço, saem da penalidade e vão para
o bloco não penalizado `A`, ao lado dos `c_ℓ`, o que muda no Teorema 1 em
blocos, e a taxa da §4 fica? No código é o `klopp.freecoarse`:
`wafc_kp_groups(balanced = TRUE, free.coarse = TRUE)`, com esses níveis no
grupo 0 do `grpreg`, os níveis finos nos pedaços balanceados e os pesos
`sqrt(|G|)`. Conferência: Parte F de [`check/08a-blocos.R`](check/08a-blocos.R).

Escrevo `j*` para o nível mais fino com `2^{j*} < b_n`, de modo que
`2^{j*} < b_n ≤ 2^{j*+1}`; cada bloco solta os níveis `0` a `j*`, isto é,
`2^{j*+1} − 1` colunas, e `𝒢_F` é a partição dos níveis `j*+1` a `J − 1` em
pedaços balanceados.

**Veredito: o Teorema 1 em blocos vale como está, com `A` trocado por
`A_F = (c_ℓ, níveis grossos)`, de posto**

```
p_0 = p + pq (2^{j*+1} − 1) = p {1 + q (2^{j*+1} − 1)},      b_n − 1 ≤ 2^{j*+1} − 1 ≤ 2b_n − 3,
```

**e o termo dos níveis passa a `σ² p_0/n`; os pedaços que sobram têm
`ρ² ≤ (2b_n − 1)/b_n < 2`; e a taxa da §4 fica, porque `σ² p_0/n` é
`O(log n/n)`.** É o `p_0` da Proposição 5 de E1.8 com `j_0 = j* + 1`: na
base periódica centrada os níveis `0` a `j*` geram o mesmo espaço que as
funções de escala do nível `j* + 1` menos a constante, e a forma livre é a
convenção clássica da limiarização de wavelets, que deixa intocados os
coeficientes abaixo de um nível inicial `j_0` (a VisuShrink de Donoho &
Johnstone, 1994, Definição 2, §4.2, p. 445, com `j_0` fixo; a razão, que os
coeficientes abaixo de `j_0` são "a fixed number, independent of n", e o
custo `2^{j_0}σ²/n` estão na §2.4, p. 440), aqui com
`j_0 = j* + 1 ≍ log_2 b_n`. A diferença para E1.8 é que lá `p_0`
não depende de `n` e aqui cresce como `b_n ≍ log n`. Nos `n` do piloto,
`2^{j*+1} − 1 = 7` (`b_n = 6` e `7`, `j* = 2`), `p_0 = 45` nas células
`q = 2` e `116` na `mixed` (`p = q = 4`).

### 12.1 O que passa, passo a passo

1. **A perfilagem (Lema 4 de E1.5).** A prova só usa que a penalidade não
   toca `A`; com `A_F` no lugar de `A`, `θ̂_F` minimiza
   `‖Ỹ − B̃_F θ‖_n² + 2λ‖θ‖_{𝒢_F,w}` com `B̃_F = M_{A_F} B_F`, e a
   decomposição `f̂ − f_J = B̃_F v + P_{A_F}(𝐛 + ε)` vale. A constante
   melhora em vez de piorar: `γ̃_F = λ_min(B̃_F'B̃_F/n)` é o complemento de
   Schur do bloco grosso dentro do complemento de Schur da forma balanceada,
   e a inversa de um complemento de Schur é submatriz principal da inversa,
   logo `γ̃_F ≥ γ̃ ≥ λ_min(Σ̂)` (Parte F2, em todas as réplicas).
2. **O posto.** `A_F` tem posto `p_0` sempre que `Σ̂` é definida positiva,
   porque as colunas de `A_F` são colunas de `Z`; é a hipótese da Proposição
   5 de E1.8, e E1.4 não muda, porque o desenho `Z` é o mesmo.
3. **A calibração (§2.3)** roda só sobre `𝒢_F`: `|𝒢_F| = |𝒢| − pq`, e
   `Ψ̃_G ⪯ Σ̂_GG` continua, porque `M_{A_F}` é projeção. A §11.1 vale com
   `ρ` dos pedaços finos.
4. **O termo dos níveis (Teorema 1(iii)).** `E[‖P_{A_F} ε‖_n² | X, U] =
   σ² p_0/n`, exato em média como na Proposição 5(iii) (Parte F2: razão de
   0,986 a 1,007 em 2 000 sorteios). A conta direta dá até um pouco mais que
   o (iii) de E1.5: `f̂ − f = B̃_F v − M_{A_F} 𝐛 + P_{A_F} ε`, logo
   `‖f̂ − f‖_n ≤ ‖B̃_F v‖_n + ‖𝐛‖_n + ‖P_{A_F} ε‖_n`.
5. **Os coeficientes grossos** saem da perfilagem, e não da cota de `v`:
   `A_F(ĉ − c*) = P_{A_F}(𝐛 + ε − B_F v)`, logo

   ```
   ‖ĉ − c*‖_2  ≤  (‖B_F v‖_n + ‖𝐛‖_n + ‖P_{A_F} ε‖_n) / sqrt(λ_min(A_F'A_F/n)),
   ```

   com `λ_min(A_F'A_F/n) ≥ λ_min(Σ̂)` (submatriz principal). O erro de uma
   componente é `‖v_{ℓm}‖² + ‖ĉ_{ℓm} − c*_{ℓm}‖²`, e a parte nova é
   `O_p(σ² p_0/(nγ))` mais os termos que o Teorema 1 já controla.

### 12.2 O `ρ²` dos pedaços que sobram

Sem o pedaço grosso, todo pedaço penalizado é fino e tem entre `b_n` e
`2b_n − 1` colunas: `ρ² = max|G|/min|G| ≤ (2b_n − 1)/b_n < 2`, contra
`(2b_n − 1)/(b_n − 1) ≤ 3` da forma balanceada (Parte F1, `b = 2` a `200`,
`J = 2` a `12`). Nos `n` do piloto o valor não muda (1,67 em `b_n = 6` e
1,57 em `b_n = 7`, a partir de `J = 5` e `J = 6`), porque o pedaço grosso de
7 colunas nunca era o extremo; em `J = 4` resta um pedaço de 8 por bloco
(`ρ = 1`), e em `J ≤ j* + 1 = 3` não resta nada para penalizar: o estimador
é o de mínimos quadrados na peneira de nível `J`, que o código ajusta com
`lm.fit` e o mesmo erro de validação cruzada do `grpreg`.

### 12.3 A taxa

**O risco ideal (§3.1 com a §11.1(ii)):** no evento `𝒯_{𝒢_F}`,

```
‖f̂ − f‖_n²  ≤  120 Λ ρ² R_{𝒢_F}(θ*_F; η_1) + 120 ‖𝐛‖_n² + 3 ‖P_{A_F} ε‖_n²,
```

com o risco ideal só sobre os pedaços finos. No lema da §3.2 os níveis
grossos saem da soma: o termo `(x_+ + 1)η`, ou o `η` por bloco da §11.3,
desaparece, e em troca entra `3σ² p_0/n` em média. A cabeça da §3.2 começa
em `j* + 1` e a cauda não muda.

**A ordem.** `σ² p_0/n ≤ σ²(p + pq(2b_n − 3))/n = O(log n/n)`, e a razão para
a taxa da §4 é da ordem de `log n · n^{−1/(2s+1)}` (dividida pelo fator
`(log n)^{(2/π−1)_+/(2s+1)}` quando `π < 2`), que vai a zero: **a taxa da §4
fica**, na mesma janela de `J_n` e sem condição nova. A convergência é lenta
quando `s` é grande. Sem as constantes da taxa, que a razão não vê, a razão
das ordens com `π = 1` e `p_0` do piloto é 2,8, 7,4 e 16 em `n = 250` para
`s = 0,8`, `1,5` e `3`, e 0,17, 1,5 e 8,9 em `n = 10^6` (Parte F4). O mesmo
vale para a leitura densa da §4.2, cuja taxa `n^{−2s'/(2s'+1)}` domina
`log n/n`. A §4.1 pedia o pedaço ao menos do tamanho de `log|𝒢|`; com
`|𝒢_F| < |𝒢|` isso só fica mais fácil.

**O Corolário 8 (§5)** continua com a estatística da norma do bloco
`(ℓ, m)`, agora com os coeficientes grossos dentro, e a cota do item 5 da
§12.1 acrescenta a `D_n²` um termo de ordem `σ² p_0/(nγ)`, menor que
`ρ_n² = n^{−2s'/(2s'+1)}`. O que se perde é a seleção por zeros: um bloco
com nível grosso nunca é zero, e a estrutura só se lê pela limiarização do
Lema 13 (ou pela parte penalizada, que é o que o piloto registra).

### 12.4 O que a troca custa e o que compra

Na cota, soltar é sempre melhor: por bloco, o pedaço grosso penalizado
custa até `120 Λ ρ² η_1 = 192 ρ² λ_1²/γ̃`, com
`λ_1² ≍ σ²(b_n + log|𝒢|)/n` (§2.3), e o livre custa
`3σ²(2^{j*+1} − 1)/n ≤ 3σ²(2b_n − 3)/n`: a mesma ordem, com a constante 3
no lugar de `192 ρ²/γ̃` vezes a de `λ_1²`. Mas a cota é cota: o custo realizado do pedaço penalizado é o
de `min(‖θ*_G‖², η)`, nulo quando os coeficientes grossos são nulos, e o do
livre é `σ²(2^{j*+1} − 1)/n` sempre. **No nulo a forma livre paga
`σ² p_0/n` que a balanceada não paga.**

Números (Parte F3; `J = 5`, `b = 6`, 7 colunas livres por bloco e pedaços
finos `{8, 6, 10}`, `p_0 = 30`, `ρ² = 1,67`; `λ` da teoria, `λ_w` sobre os
pedaços finos para a livre e o da §11 para a balanceada; 40 réplicas por
linha; média de `‖f̂ − f‖_n²`):

| verdade | `n` | forma livre | forma balanceada | livre / `σ²p_0/n` |
|---|---|---|---|---|
| dois pedaços finos cheios, grossos nulos | 250, 375, 500 | 0,888, 0,498, 0,340 | 0,827, 0,466, 0,322 | 30, 25, 23 |
| os mesmos, com os grossos de dois blocos ligados | 250, 375, 500 | 0,815, 0,489, 0,348 | 2,173, 1,639, 1,219 | 27, 24, 23 |
| `g_{11} = sin(2πu)`, com viés | 250, 375, 500 | 0,0295, 0,0215, 0,0143 | 0,265, 0,160, 0,111 | 0,98, 1,07, 0,96 |
| nulo | 250, 375, 500 | 0,0301, 0,0185, 0,0143 | 0,0027, 0,0011, 0,0010 | 1,00, 0,92, 0,96 |

- **As três cotas do Teorema 1(ii) em blocos valem nos 480 ajustes**, com
  `W` e `λ_w` só sobre os pedaços finos e folga mínima de 444 (no nulo e no
  seno o `λ` da teoria zera todos os pedaços finos); a (iii) com o `A_F` e a
  cota dos coeficientes grossos também; KKT a menos de `10^{−9}`.
- **Onde há energia nos níveis grossos a forma livre ganha muito**: 2,7 a
  3,5 vezes com os grossos ligados e 7,4 a 9,0 vezes no seno, cuja energia está
  quase toda nos níveis 0 a 2. **No nulo perde por 11 a 16 vezes**, e o seu
  erro é `σ² p_0/n` a menos de 8%, como a §12.1 prevê; o da balanceada é
  `σ² p/n`. Com os grossos nulos e pedaços finos ativos as duas empatam,
  porque o erro está no encolhimento dos pedaços finos.
- **Leitura para o piloto, como previsão e não medição:** com `σ = 0,62` no
  nulo e `p_0 = 45`, o termo `σ² p_0/n` vale 0,069 em `n = 250` e 0,017 em
  `n = 1000`, contra um `rmse_f²` do `wafc.lasso` da ordem de 0,007 e 0,002
  (E2.5a); no `smooth`, cuja energia está nos níveis grossos, é onde a forma
  livre pode recuperar o que a balanceada perde. O `λ` da teoria
  superencolhe (razão erro/`σ²p_0/n` de 23 a 30 com pedaços finos ativos), e
  com `λ` por validação cruzada os números são outros: a medição é a de
  E2.5f, no `docs/handoff-E2.5f.md`.
- **Medido em E2.5f** (validação cruzada, 50 réplicas por célula): a
  previsão do nulo se confirma (ISE 5 a 9 vezes o da forma balanceada; a
  validação cruzada escolhe `J = 2`, mínimos quadrados em 21 colunas, em 98%
  a 100% das réplicas), e o custo dos níveis livres nos blocos **inativos** domina
  também onde há componente: no `smooth` o ISE dos blocos ativos cai até 26%
  em `n = 1000`, mas o dos inativos sobe 3 a 6 vezes, e a forma livre é a
  pior das cinco em quase toda célula. A troca da cota (§12.4) favorece
  soltar porque conta só os blocos ativos; o termo `σ² p_0/n` é pago por
  todos os `pq` blocos.

### 12.5 O que isto não cobre

- **O `λ` por validação cruzada**, como em toda a sondagem, e a escolha de
  `J` entre um ajuste penalizado (`J ≥ 4`) e um de mínimos quadrados
  (`J ≤ 3`): o primeiro é o desta seção, o segundo é a peneira por mínimos
  quadrados, cuja cota é `σ²(p + pq(2^J − 1))/n` mais o viés de `V_J`, e a
  validação cruzada compara os dois sem teorema que a sustente.
- **`p_0` crescente.** A Proposição 5 de E1.8 tinha `p_0` fixo; aqui
  `p_0 ≍ pq log n`, e o argumento da taxa usa só que `σ² p_0/n = O(log n/n)`.
  Nenhum passo do Teorema 1 usa `p_0` fixo; a condição de posto de `A_F`
  vem do autovalor de E1.4, que já é sobre `Z` inteiro.
- **Não é prova**, como o resto do documento: a §12.1 é a releitura das §2.1,
  §2.3 e §11.1 com `A_F`, e a §12.3, a da §3.2 sem os níveis grossos; as
  duas estão conferidas na Parte F.
