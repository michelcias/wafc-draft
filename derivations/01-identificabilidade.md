# 01. Identificabilidade (E1.2)

Resultado da etapa E1.2 do plano. Enunciado, hipóteses, prova, o que muda com
`boundary = "interval"`, o que não cobre e a conferência numérica
(`check/01-identificabilidade.R`, que imprime `OK`). Notação de
[`../docs/notacao.md`](../docs/notacao.md), congelada em E1.1. Numeração
global dos resultados (`TAREFA.md`, §5): este arquivo registra a
**Proposição 1** e o **Lema 1**; as hipóteses são numeradas localmente como
(A1) a (A3) até o chat principal congelar a lista de E1.

**Estado: conferido numericamente (2026-09-18) e provado.** O que se afirma
aqui é qualitativo (injetividade, posto cheio); nada quantitativo sobre
`λ_min` da Gram, que é E1.4.

---

## 1. O problema

O modelo (D1) é

```
Y = Σ_{ℓ=1}^{p} β_ℓ(U) X_ℓ + ε,            E[ε | X, U] = 0,
β_ℓ(u) = c_ℓ + Σ_{m=1}^{q} g_{ℓm}(u_m),    u ∈ [0,1]^q,
```

com `f(x, u) = Σ_ℓ β_ℓ(u) x_ℓ = E[Y | X = x, U = u]`. A distribuição de
`(Y, X, U)` determina `f` como elemento de `L_2(P_{X,U})`, isto é, a menos
de conjuntos `P_{X,U}`-nulos. Identificar o modelo é mostrar que a aplicação

```
(c_1, …, c_p, (g_{ℓm})_{ℓ,m})  ↦  f
```

é injetiva. Há duas fontes de não unicidade que qualquer hipótese tem de
excluir:

1. **Entre coeficientes funcionais.** Se `X` for, dado `U`, confinado a um
   hiperplano (por exemplo `X_2 = h(U_1)` com `X_1 ≡ 1`), então
   `β_1(u) + β_2(u) h(u_1)` só determina uma combinação de `β_1` e `β_2`.
2. **Dentro de um coeficiente.** A decomposição aditiva `c + Σ_m g_m(u_m)`
   não é única sem centralização (uma constante passa de `g_m` para `c`),
   nem quando as moduladoras são funções uma da outra (`U_2 = h(U_1)`:
   `g_1(u_1) + g_2(h(u_1))` só determina a soma).

A primeira fonte é nova em relação ao WALL (o produto `X_ℓ ψ_{jk}(U_m)`); a
segunda é a dos modelos aditivos. A escolha da centralização é o ponto que
liga a hipótese ao desenho: no caso periódico com `j_0 = 0` a base
**impõe** `∫_0^1 g_{ℓm} = 0` (Lema 1), o que é o que D2 antecipava.

## 2. Hipóteses

- **(A1) Desenho linear não degenerado dado `U`.** `E‖X‖² < ∞` e a matriz
  `E[X X' | U]` é não singular `P_U`-quase certamente. Com `X_1 ≡ 1` isso
  equivale a `Var(X_{−1} | U)` não singular q.c., onde `X_{−1}` são as demais
  coordenadas: as covariáveis lineares têm de variar dentro de cada valor de
  `U`. Dependência entre `X` e `U` é permitida.
- **(A2) Moduladoras com suporte cheio.** `U` tem densidade `p_U` em relação
  à medida de Lebesgue em `[0,1]^q`, com `p_U > 0` Lebesgue-q.c. A hipótese
  D11 (`0 < c ≤ p_U ≤ C`) implica (A2); aqui só a positividade é usada.
  Dependência entre as `U_m` é permitida; o que fica excluído é uma
  moduladora ser função de outra (a lei conjunta seria singular).
- **(A3) Centralização de Lebesgue e integrabilidade.** `c_ℓ ∈ ℝ`,
  `g_{ℓm} ∈ L_2[0,1]` (medida de Lebesgue) e `∫_0^1 g_{ℓm}(u) du = 0` para
  todo `(ℓ, m)`.

O `notacao.md` (§2) escreve a restrição como `E[g_{ℓm}(U_m)] = 0`. A
Proposição 1 vale com qualquer das duas (Observação 1); a versão de
Lebesgue é a que a base impõe de graça (Lema 1) e a que o WALL adota
(`ms_theo_1.tex`, §2.1, "Lebesgue centering"). A escolha é proposta ao chat
principal no handoff.

## 3. Enunciados

**Proposição 1 (identificabilidade populacional).** Sob (A1) a (A3):

(a) Se `β = (β_1, …, β_p)` e `β̃` são funções mensuráveis `[0,1]^q → ℝ^p`
com `Σ_ℓ β_ℓ(U) X_ℓ = Σ_ℓ β̃_ℓ(U) X_ℓ` `P_{X,U}`-q.c., então `β_ℓ = β̃_ℓ`
`P_U`-q.c. para todo `ℓ`. Só (A1) é usada.

(b) Se `c_ℓ + Σ_m g_{ℓm}(u_m) = c̃_ℓ + Σ_m g̃_{ℓm}(u_m)` para `P_U`-q.t.
`u`, com `(c_ℓ, g_{ℓm})` e `(c̃_ℓ, g̃_{ℓm})` satisfazendo (A3), então
`c_ℓ = c̃_ℓ` e `g_{ℓm} = g̃_{ℓm}` Lebesgue-q.c. em `[0,1]`, para todo `m`.
Só (A2) e (A3) são usadas.

Em consequência, `(c, (g_{ℓm})) ∈ ℝ^p × Π_{ℓ,m} L_2[0,1]` é determinado por
`f ∈ L_2(P_{X,U})`.

**Lema 1 (a base periódica com `j_0 = 0` impõe (A3) e parametriza o sieve
injetivamente).** Seja `{ψ_{jk}}` a base de wavelets periodizada em `[0,1]`
gerada por um par `(φ, ψ)` ortonormal de suporte compacto com `N ≥ 1`
momentos nulos, e `V_J = span{φ_{00}} ⊕ W_J`, `W_J = span{ψ_{jk} : 0 ≤ j < J}`.

(i) `φ_{00} ≡ 1` e `∫_0^1 ψ_{jk}(u) du = 0` para todo `(j, k)`. Logo
`W_J = {g ∈ V_J : ∫_0^1 g = 0}` e toda função `g = Σ_{jk} θ_{jk} ψ_{jk}`
satisfaz (A3) sem restrição sobre `θ`.

(ii) Sob (A1) e (A2), para cada `J` fixo a aplicação
`(c, θ) ∈ ℝ^p × ℝ^d ↦ f_{c,θ}(x, u) = Σ_ℓ x_ℓ [c_ℓ + Σ_m Σ_{jk} θ_{ℓm,jk} ψ_{jk}(u_m)]`
é injetiva em `L_2(P_{X,U})`. Equivalentemente, a Gram populacional
`Σ = E[Z Z']`, com `Z` na ordem D12, é definida positiva.

(iii) Se a coluna `X_ℓ φ_{00}(U_m) = X_ℓ` for mantida em cada bloco, `Σ` é
singular, com nulidade exatamente `p q`.

## 4. Provas

**Prova da Proposição 1(a).** Seja `δ(u) = β(u) − β̃(u) ∈ ℝ^p`. A hipótese
diz que `δ(U)' X = 0` q.c., logo `(δ(U)' X)² = 0` q.c. e a esperança
condicional dessa variável não negativa dado `U` é zero q.c. Como
`E‖X‖² < ∞` e `δ(U)` é `σ(U)`-mensurável,

```
0 = E[(δ(U)' X)² | U] = δ(U)' E[X X' | U] δ(U)     P_U-q.c.
```

(nenhum momento de `β` é necessário: a igualdade vale em `[0, ∞]` e o lado
esquerdo é zero). Por (A1), `E[X X' | U]` é definida positiva q.c., e uma
forma quadrática definida positiva só se anula em zero; logo `δ(U) = 0`
`P_U`-q.c. ∎

**Prova da Proposição 1(b).** Por (A2), um conjunto `A ⊂ [0,1]^q` é
`P_U`-nulo se e só se é Lebesgue-nulo: `P_U(A) = ∫_A p_U = 0` força
`Leb(A ∩ {p_U > 0}) = 0`, e `{p_U = 0}` é Lebesgue-nulo. Assim a igualdade
da hipótese vale Lebesgue-q.c. em `[0,1]^q`. Escreva `d = c_ℓ − c̃_ℓ` e
`δ_m = g_{ℓm} − g̃_{ℓm} ∈ L_2[0,1] ⊂ L_1[0,1]`:

```
d + Σ_{m=1}^{q} δ_m(u_m) = 0     Lebesgue-q.c. em [0,1]^q.
```

Cada `δ_m(u_m)` é integrável no cubo (Tonelli e `δ_m ∈ L_1`). Integrando
em todas as coordenadas, `d + Σ_m ∫_0^1 δ_m = d = 0`, porque cada `δ_m` tem
integral zero por (A3). Fixado `m`, integrando nas demais coordenadas
(Fubini), para Lebesgue-q.t. `u_m`:
`δ_m(u_m) + Σ_{m' ≠ m} ∫_0^1 δ_{m'} = δ_m(u_m) = 0`. ∎

**Observação 1 (centralização em `P`).** Se (A3) for trocada por
`g_{ℓm} ∈ L_2(P_{U_m})` com `E[g_{ℓm}(U_m)] = 0`, a prova é a mesma sob
(A2) reforçada para D11 (`p_U ≤ C` dá `L_2(P_{U_m}) ⊃ L_2[0,1]` e
`p_U ≥ c` dá a inclusão inversa, logo os dois espaços coincidem): o passo de
Fubini dá `δ_m` constante q.c., e a constante é `E[δ_m(U_m)] = 0`. As duas
parametrizações estão ligadas por `c^{P}_ℓ = c_ℓ + Σ_m E[g_{ℓm}(U_m)]` e
`g^{P}_{ℓm} = g_{ℓm} − E[g_{ℓm}(U_m)]`. A diferença não é cosmética para o
estimador: com `U_2 ~ Beta(2,3)` a média `P` de `ψ_{jk}(U_2)` chega a
`0.48` (§7), e o que o desenho sem `φ_{00}` estima é a versão centrada em
Lebesgue.

**Prova do Lema 1(i).** A periodização é `φ^{per}_{00}(u) = Σ_{l∈ℤ} φ(u + l)`
e `ψ^{per}_{jk}(u) = Σ_{l∈ℤ} 2^{j/2} ψ(2^j (u + l) − k)`, somas finitas em
cada `u` porque os suportes são compactos. A primeira é a partição da
unidade `Σ_l φ(x − l) ≡ 1`, válida para toda função de escala ortonormal de
suporte compacto com `∫ φ = 1` (a periodização de `φ` é a §9.3 de Daubechies
1992, "Wavelets for `L^p([0,1])`", p. 304; a identidade em si não aparece
enunciada lá, e L3 não achou fonte verificada para ela — item do handoff);
o `wbasis()` devolve essa coluna igual a `1` até a última
casa decimal (§7). Para a segunda, trocando soma e integral (soma finita),
`∫_0^1 ψ^{per}_{jk} = ∫_ℝ 2^{j/2} ψ(2^j t − k) dt = 2^{−j/2} ∫_ℝ ψ = 0` pelo
momento nulo de ordem zero. Que `{1} ∪ {ψ^{per}_{jk}}` seja base ortonormal
de `L_2[0,1]` é a construção padrão da análise de multirresolução
periodizada (Daubechies 1992, §9.3; Härdle, Kerkyacharian, Picard &
Tsybakov 1998), e é o que o WALL assume (`ms_theo_1.tex`,
§2.2). Como `dim V_J = 2^J` e as `2^J − 1` wavelets de nível `< J` são
ortogonais a `1`, `W_J` é exatamente o complemento das constantes em `V_J`,
e toda `g ∈ W_J` tem integral zero. ∎

**Prova do Lema 1(ii).** Suponha `f_{c,θ} = f_{c̃,θ̃}` `P_{X,U}`-q.c. As
funções `β^J_ℓ(u) = c_ℓ + Σ_m g^J_{ℓm}(u_m)`, com
`g^J_{ℓm} = Σ_{jk} θ_{ℓm,jk} ψ_{jk}`, são mensuráveis e limitadas. Pela
Proposição 1(a), `β^J_ℓ = β̃^J_ℓ` `P_U`-q.c. Cada `g^J_{ℓm}` satisfaz (A3)
pelo item (i), logo, pela Proposição 1(b), `c_ℓ = c̃_ℓ` e
`Σ_{jk} (θ − θ̃)_{ℓm,jk} ψ_{jk} = 0` Lebesgue-q.c.; a ortonormalidade das
`ψ_{jk}` dá `θ = θ̃`. Para a equivalência com a Gram: `E‖Z‖² < ∞` porque
`E‖X‖² < ∞` e cada `ψ_{jk}` é limitada; para `a = (a_c, a_θ) ≠ 0`,
`a' Σ a = E[(a' Z)²] = ‖f_{a_c, a_θ}‖²_{L_2(P)} > 0` pela injetividade. ∎

**Prova do Lema 1(iii).** Com `φ_{00} ≡ 1`, a coluna `X_ℓ φ_{00}(U_m)`
coincide com a coluna não penalizada `X_ℓ`, para cada um dos `p q` pares
`(ℓ, m)`. Os `p q` vetores `e_{X_ℓ} − e_{(ℓm, φ)}` são linearmente
independentes (suportes disjuntos nas coordenadas `(ℓm, φ)`) e estão no
núcleo de `Σ`. Não há mais: um vetor do núcleo define
`Σ_ℓ X_ℓ [a_ℓ + Σ_m (α_{ℓm} + Σ_{jk} a_{ℓm,jk} ψ_{jk}(U_m))] = 0` q.c.; pela
Proposição 1(a) e pelo argumento de (b) aplicado à parte de integral zero,
`a_{ℓm,jk} = 0` e `a_ℓ + Σ_m α_{ℓm} = 0`, o que é o espaço de dimensão `p q`
gerado pelos vetores acima. ∎

## 5. O que muda com `boundary = "interval"`

A base de Cohen, Daubechies e Vial (1993; `cohen1993wavelets`), como
implementada no `wbasis()`, mantém em cada nível `j` os `2^j − L` transladados
interiores (`L` é o tamanho do filtro, `N = L/2` momentos nulos) e os
completa com `L/2` funções de borda em cada extremo, de modo que os
polinômios de grau `≤ N − 1` são reproduzidos em `[0,1]`. O `wbasis()` exige
`j_0 ≳ log_2(5L)`: `j_0 ≥ 3` com `L = 4`, `j_0 ≥ 4` com `L = 6` ou `8`,
`j_0 ≥ 6` com o `L = 20` padrão do `wall()`. Três coisas mudam em relação ao
Lema 1:

1. **A constante está em `V_{j_0}`, mas não é uma coluna.** Como `1` é
   polinômio de grau zero, `1 ∈ V_{j_0}`, e na base ortonormal
   `{φ^{int}_{j_0 k}}_{k=0}^{2^{j_0}−1}` de `V_{j_0}` ela se escreve
   `1 = Σ_k μ_k φ^{int}_{j_0 k}`, `μ_k = ⟨1, φ^{int}_{j_0 k}⟩ = ∫_0^1 φ^{int}_{j_0 k}`
   (`μ_k = 2^{−j_0/2}` nos transladados interiores; nas funções de borda,
   valores próprios, com alguns nulos). Não há uma única coluna a descartar.
2. **As wavelets continuam centradas; as funções de escala, não.** Cada
   `ψ^{int}_{jk}` com `j ≥ j_0` é ortogonal a `V_{j_0} ∋ 1`, logo
   `∫ ψ^{int}_{jk} = 0`, e a parte penalizada do bloco satisfaz (A3) sem
   restrição, como no caso periódico. As `2^{j_0}` funções de escala não:
   `{g ∈ V_J : ∫ g = 0} = {Σ_k α_k φ^{int}_{j_0 k} + Σ_{jk} θ_{jk} ψ^{int}_{jk} : μ' α = 0}`,
   de dimensão `2^J − 1`, com uma restrição linear explícita por bloco.
3. **Duas implementações equivalentes da restrição.** (a) Manter as `2^{j_0}`
   colunas de escala por bloco e impor `μ' α_{ℓm} = 0` (`p q` restrições
   lineares; o `glmnet` não as aceita diretamente). (b) Reparametrizar
   `α_{ℓm} = Q γ_{ℓm}`, com `Q ∈ ℝ^{2^{j_0} × (2^{j_0} − 1)}` base
   ortonormal de `{α : μ' α = 0}`; o bloco `(ℓ, m)` passa a ter colunas
   `[X_ℓ Φ(U_m) Q | X_ℓ Ψ(U_m)]`, com `2^{j_0} − 1` colunas de escala não
   penalizadas e `2^J − 2^{j_0}` wavelets penalizadas, total `2^J − 1` como
   no caso periódico. As colunas `Φ Q` são ortonormais em `L_2[0,1]` (`Q'Q = I`)
   e têm integral zero (`μ' Q = 0`), de modo que o Lema 1(ii) vale
   literalmente com `[Φ Q | Ψ]` no lugar de `Ψ`: Proposição 1(b) mais
   ortonormalidade. Manter as `2^{j_0}` colunas sem restrição repete a
   constante em cada bloco: `Σ` singular com nulidade `p q`, como em (iii).
   É o que o `wall()` faz hoje com `boundary = "interval"` (`drop.phi = FALSE`,
   todas as funções de escala não penalizadas, mais o intercepto do
   `glmnet`), o que ali não compromete a função ajustada, só a leitura dos
   coeficientes de escala; registrado no handoff.

O custo da opção é o número de parâmetros não penalizados:
`p q (2^{j_0} − 1)`, isto é, `15 p q` com `L = 8` e `63 p q` com `L = 20`,
contra zero no caso periódico. Se a parte de escala for penalizada também, a
reparametrização (b) deixa de ser inócua (`‖γ‖_1 ≠ ‖Q γ‖_1`); é decisão de
E2, não de identificabilidade.

## 6. O que isso não cobre

- **Nada quantitativo.** O Lema 1(ii) dá `λ_min(Σ) > 0` para cada `J`
  fixo; a taxa em que pode ir a zero com `J`, a versão restrita ao cone e a
  constante em função de `κ_1 = inf λ_min(E[XX' | U])` e de `c ≤ p_U` são
  E1.4. Toda taxa de E1.5 e E1.6 depende disso, não deste arquivo.
- **A Gram empírica.** `Σ̂ = Z'Z/n` definida positiva exige `n ≥ p + d` e
  observações em toda célula diádica onde alguma coluna vive. A conferência
  mostra a falha no caso de intervalo: com `U_2 ~ Beta(2,3)`, `n = 200` e
  `j_0 = 4`, a última célula `(15/16, 1]` fica vazia e as oito funções de
  borda direita do bloco `U_2` ficam linearmente dependentes na amostra
  (nulidade extra `2`). É fenômeno de amostra finita, coberto em E1.4(iii)
  por concentração sob `p_U ≥ c`; na prática, o reescalonamento de cada
  `U_m` para `[0,1]` pela amplitude amostral (o que o `wall()` faz com
  `boundary = "interval"`) povoa as células extremas por construção.
- **Só a menos de conjuntos nulos.** `β_ℓ` fica identificada `P_U`-q.c. e
  `g_{ℓm}` Lebesgue-q.c.; nada pontual. Continuidade das `g_{ℓm}` (que
  E1.3 pode dar sob Besov com `s > 1/π`) transforma isso em identificação
  em todo ponto.
- **Moduladoras discretas ou com atomos, ou suportadas em conjunto de
  dimensão menor.** Excluídas por (A2) e genuinamente não identificadas no
  caso funcional (`U_2 = h(U_1)`). Não é limitação da prova, é do modelo.
- **`X` função de `U`.** `X_2 = h(U_1)` viola (A1) e o modelo não identifica
  `(β_1, β_2)`; a conferência exibe a deficiência de posto com `h ∈ V_J`.
  (A1) permite dependência entre `X` e `U`, desde que reste variação
  condicional: com `X_1 ≡ 1` e `p = 2`, `det E[XX' | U] = Var(X_2 | U)`.
- **Unicidade da solução do LASSO** (`Σ̂` singular ou não) é outra
  questão, de posição geral das colunas (R. J. Tibshirani 2013);
  E1.5 não precisa dela.
- **O reescalonamento empírico de `U`** para `[ε, 1 − ε]` (computação, fora
  da teoria por D11) é uma bijeção monótona por coordenada; preserva a
  identificabilidade, mas troca a centralização de Lebesgue pela
  centralização na coordenada reescalada.
- **Viés de sieve.** Se `g_{ℓm} ∉ W_J`, o objeto identificado no sieve é a
  projeção `Π_J g_{ℓm}`, com o mesmo enunciado aplicado a ela. A distância
  entre as duas é E1.3.

## 7. Conferência numérica

`check/01-identificabilidade.R`, executado em 2026-09-18 com o `WaveBased`
2.6-0 (Daublets, `filter.size = 8`, `N = 4`), `n = 200`, `p = 2`
(`X_1 ≡ 1`, `X_2 ~ N(0,1)`), `q = 2` (`U_1 ~ Unif`, `U_2 ~ Beta(2,3)`,
cópula gaussiana `ρ = 0.6`: dependência entre moduladoras e marginal não
uniforme, de propósito), forma densa, tolerância de posto
`σ/σ_max > 10^{−8}`. Imprime `OK`.

| Conferência | Resultado |
|---|---|
| Lema 1(i), periódico `j_0 = 0`, `J = 3`, grade de `2^{14}` pontos | `φ_{00} ≡ 1` exato; `max |∫ψ_{jk}| = 5.4e−17`; `max |Gram − I| = 2.0e−09` |
| Lema 1(ii), posto de `Z` (`30 = 2 + 4 · 7` colunas) | posto `30`; `σ_min/σ_max = 5.3e−02`; `λ_min(Z'Z/n) = 7.8e−03` |
| Recuperação de `(c*, θ*)` sem ruído por mínimos quadrados, `y` calculado das funções `β_ℓ(U_i) X_iℓ` | `max |erro| = 3.0e−15` |
| Lema 1(iii), `φ_{00}` mantida (`38` colunas) | nulidade `4 = p q` |
| Fora de (A2): `U_2 = U_1` | nulidade `14 = p N_J` |
| Fora de (A1): `X_2 = 1 + ψ_{10}(U_1)` | nulidade `1` |
| Observação 1: `max |mean ψ_{jk}(U_2)|`, `U_2 ~ Beta(2,3)` | `0.481` (contra `5.4e−17` em Lebesgue) |
| Intervalo, `j_0 = 4`, `J = 5`, grade `2^{14}` | `max |Gram − I| = 6.1e−06` (quadratura `O(h²)`: `3.8e−07` com `2^{16}`); `1 = Σ μ_k φ_k` com resíduo `3.1e−09`; `max |∫ψ^{int}| = 2.9e−07`; `max |∫φ^{int}| = 0.64` |
| Intervalo, escala mantida (`130` colunas), `U` populacional | nulidade `6 = p q + 2` (`max U_2 = 0.883`, célula `(15/16, 1]` vazia) |
| Intervalo, escala mantida, `U` reescalado à amplitude amostral | nulidade `4 = p q`; as `4` direções nulas aparecem em `σ/σ_max ≈ 2e−11` a `7e−11` (precisão da avaliação CDV), as seguintes em `2.2e−03` |
| Intervalo, `[Φ Q | Ψ]` (`126 = 2 + 4 · 31` colunas) | posto cheio; `λ_min(Z'Z/n) = 3.9e−05`; recuperação `max |erro| = 1.4e−13` |

Sinal para E1.4 e E2: a `λ_min` empírica do desenho de intervalo no menor
`j_0` admissível é duas ordens de grandeza abaixo da do periódico com
`J = 3`, em parte por ter `126` colunas contra `30` com o mesmo `n`.

## 8. O que foi transposto do WALL e o que é novo

- **Transposto:** a convenção de centralização de Lebesgue e o argumento de
  que a ausência de constante em `W_J` a impõe (`ms_theo_1.tex`, §2.1 e
  §2.2); a leitura de que a periodização é o que faz `{1} ∪ {ψ_{jk}}` ser
  base ortonormal.
- **Novo:** a Proposição 1(a), que é o preço do produto `X_ℓ ψ_{jk}(U_m)` e
  não tem análogo no WALL; a Proposição 1(b) para `q` moduladoras dentro
  de cada coeficiente (o WALL tem uma soma aditiva só, e a Observação
  `rem:identif` diz que lá a identificação das componentes nem é
  necessária, porque só a log-odds agregada entra no risco). No WAFC as
  `β_ℓ` são o objeto de interpretação, e a identificação é parte do
  enunciado. A §5 (intervalo) não tem análogo no WALL, que só usa a base
  periódica.

## Referências citadas

Todas estão em `docs/referencias-verificadas.bib` desde L3 (2026-09-19):
`cohen1993wavelets` (Cohen, Daubechies & Vial 1993) e `Xue-Yang-2006`
(origem do modelo de coeficientes aditivos), já de L1, mais `Daubechies-1992`
(*Ten Lectures on Wavelets*, CBMS-NSF 61, SIAM), `hardle1998wavelets`
(Härdle, Kerkyacharian, Picard & Tsybakov 1998, copiada do WALL) e
`Tibshirani-2013` (Ryan J. Tibshirani, "The lasso problem and uniqueness",
*Electron. J. Statist.* 7, 1456–1490 — **Ryan**, não Robert).

O que L3 não conseguiu conferir: a partição da unidade `Σ_l φ(x − l) ≡ 1`
atribuída ao cap. 5 de Daubechies (1992). A expressão "partition of unity"
não ocorre no livro (busca de texto integral no Google Books), e a §9.3, que
constrói a base periodizada, não a enuncia. A âncora ficou provisória e está
no handoff de L3.
