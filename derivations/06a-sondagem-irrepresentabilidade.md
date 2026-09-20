# E1.7a. Sondagem: irrepresentabilidade em grupos no desenho de produtos

Documento de **sondagem**, não de prova (D28). A pergunta é se a saída (a) de
[`../docs/selecao-estrutura.md`](../docs/selecao-estrutura.md) — a propriedade
oráculo de seleção para a variante em grupos — é viável no desenho de produtos,
e a entrega é um veredito com números, não um teorema. Não há numeração global
de resultados aqui; o que se prova abaixo é uma identidade algébrica, e ela é
nomeada, não numerada.

Conferência numérica:
[`check/06a-irrepresentabilidade.R`](check/06a-irrepresentabilidade.R)
(rodou e imprimiu `OK` em 2026-09-20). Notação: `../docs/notacao.md`
(congelada em E1.1).

---

## 0. Veredito

**Vale a pena tentar, com o escopo reduzido, e não antes de E2.5.** Em três
linhas:

1. **A condição não falha por construção.** Sob `E(XX' | U) = Ω` constante e
   moduladoras não correlacionadas através da base, a condição de
   irrepresentabilidade em grupos do desenho de produtos **se reduz
   exatamente** a uma condição que não envolve a base de wavelets: nem `ψ`,
   nem `J`, nem `N_J`, nem as densidades marginais das `U_m`. Esse era o risco
   principal da saída (a), e ele não se materializou.
2. **Mas não se reduz a `Ω` sozinha**, como o chat principal conjecturou. Sobra
   um fator: a matriz de correlações `L_2` entre as componentes verdadeiras que
   dividem a mesma moduladora. A condição continua sendo conjunta sobre desenho
   e verdade, como em Bach (2008); o que a redução compra é que ela vive em
   dimensão `|S_m| × |S_m|` em vez de `p q N_J`.
3. **A redução resolve o risco, não o custo.** O que falta para a saída (a) é
   um teorema de seleção para o LASSO em grupos com `d = p q N_J → ∞`, o passo
   do viés de aproximação, e uma hipótese de independência entre moduladoras
   mais forte que D11/D13. Isso continua sendo da ordem de E1.4 mais E1.5, como
   `selecao-estrutura.md` estimou, e continua condicionado a E2.5 escolher a
   variante em grupos, o que os números de E2.2 e E2.3 não indicam.

O que se recomenda incorporar já: **a redução como observação**, meia página no
artigo ou uma observação no `03-desenho-produtos.tex`. Ela é um fato estrutural
sobre o desenho de produtos, barata, e é o único pedaço de (a) que não depende
de E2.5.

O que se recomenda **não** fazer agora: abrir a prova da propriedade oráculo em
grupos. A saída (c) (E1.7c) continua sendo a que entrega teorema por preço de
corolário, e continua valendo para o estimador base.

---

## 1. A condição, na forma da literatura

Bach (2008, JMLR 9, 1179--1225), Seção 2.5, condições (4) e (5): com grupos
`G_1, …, G_M`, pesos `d_b`, coeficiente populacional `w` e conjunto ativo `J`,
a condição forte, **suficiente** para consistência de padrão (o Teorema 2 dele),
é

```
max_{b fora de J}  (1/d_b) || Σ_{G_b, J} Σ_{J,J}^{-1} Diag(d_i/||w_i||) w_J ||_2  <  1,
```

e a versão com `≤ 1` é **necessária** (o Teorema 3 dele). Quando os grupos têm
dimensão um, as duas são a condição de Zhao e Yu (2006) para o LASSO.

No WAFC os grupos são os blocos `(ℓ, m)`, todos de tamanho `N_J`, de modo que
os pesos `d_b` se cancelam, e `Diag(d_i/‖w_i‖)w_J` é a concatenação das direções
unitárias

```
ζ_{ℓm} = θ*_{ℓm} / ||θ*_{ℓm}||_2  em  R^{N_J},   (ℓ, m) em S,
```

com `θ*` o minimizador populacional do risco quadrático sobre as colunas do
desenho (o `w` de Bach, definido sob má especificação: a hipótese (A3) dele não
exige que a regressão seja linear nas colunas).

Falta um ajuste que não está em Bach: as `p` colunas não penalizadas `X_ℓ`
(os níveis `c_ℓ`, D3). Elas entram por **complemento de Schur**. Escrevendo
`C` para o conjunto dessas colunas e

```
Σ~_{A,B} = Σ_{A,B} − Σ_{A,C} Σ_{C,C}^{-1} Σ_{C,B},
```

a estatística que decide é

```
IRR(ℓ0, m0) = || Σ~_{(ℓ0,m0), S} Σ~_{S,S}^{-1} ζ_S ||_2 ,   (ℓ0, m0) fora de S.
```

Vale a pena registrar o que Bach observa na página 1183, porque é exatamente a
comparação entre as saídas (a) e (c) de `selecao-estrutura.md`: a condição
contrasta com limiarizar o estimador de mínimos quadrados conjunto, que é
consistente **sem condição nenhuma** além da inversibilidade de `Σ`, desde que o
limiar seja escolhido abaixo da menor norma ativa. É o argumento de (c),
escrito pela própria fonte de (a).

---

## 2. Frente (i): a álgebra

### 2.1 A Gram, depois do complemento de Schur

Seja `Ψ(u) = (1, ψ_{jk}(u_1), …, ψ_{jk}(u_q))' ∈ R^K`, `K = 1 + q N_J`, como em
`03-desenho-produtos.tex`, e `Σ_Ψ = E[Ψ(U)Ψ(U)']`. A parte (ii) da Proposição 3
de E1.4 dá, quando `E(XX' | U) = Ω` é constante quase certamente (o que a
independência `X ⊥ U` implica, e que é estritamente mais fraco),

```
Σ = Π (Ω ⊗ Σ_Ψ) Π' .
```

As colunas não penalizadas são exatamente as coordenadas `(ℓ, 0)` desse produto
de Kronecker, porque `Ψ_0 ≡ 1` e a coluna `X_ℓ · 1` é a coluna não penalizada de
`X_ℓ`. Como `(Σ_Ψ)_{00} = 1`, o complemento de Schur sai em uma linha:

```
Σ~ = Ω ⊗ Σ_Ψ,rest − (Ω ⊗ E[ψ]) (Ω^{-1} ⊗ 1) (Ω ⊗ E[ψ]')
   = Ω ⊗ ( Σ_Ψ,rest − E[ψ] E[ψ]' )
   = Ω ⊗ T ,        T = Cov( ψ_{jk}(U_m), ψ_{j'k'}(U_{m'}) ) ,
```

uma matriz `q N_J × q N_J`. Ou seja: **perfilar os níveis `c_ℓ` troca o segundo
momento da base pela covariância da base, e preserva a forma de Kronecker.**
É por isso que a densidade marginal de `U_m` não sobrevive à redução, embora
ela faça `E[ψ_{jk}(U_m)] ≠ 0` e portanto as colunas não penalizadas deixarem de
ser ortogonais às penalizadas (medido na §3.4: até `0.44` no caso Beta).

### 2.2 A redução

Suponha, além disso, que `T` seja **bloco-diagonal na moduladora**, isto é

```
Cov( ψ_{jk}(U_m), ψ_{j'k'}(U_{m'}) ) = 0   para  m ≠ m'  e  todo (j,k), (j',k').   (BD)
```

(BD) vale, por exemplo, quando `U_1, …, U_q` são independentes, com densidades
marginais **quaisquer**, limitadas ou não. Escreva `T = ⊕_m T_m`, `S_m = {ℓ :
(ℓ, m) ∈ S}` e, para um candidato `(ℓ0, m0)` fora de `S`,

```
a = Ω_{ℓ0, S_{m0}} ( Ω_{S_{m0}, S_{m0}} )^{-1}  ∈ R^{|S_{m0}|} ,
```

o vetor de coeficientes da regressão populacional (em segundos momentos) de
`X_{ℓ0}` sobre `{X_ℓ : ℓ ∈ S_{m0}}`. Então

```
IRR(ℓ0, m0) = || Σ_{ℓ ∈ S_{m0}} a_ℓ ζ_{ℓ m0} ||_2 = sqrt( a' V^{(m0)} a ) ,
V^{(m0)}_{ℓℓ'} = < ζ_{ℓ m0} , ζ_{ℓ' m0} > .
```

**Prova.** Por (BD), `Σ~ = Ω ⊗ (⊕_m T_m)`, e o bloco `(ℓ,m) × (ℓ',m')` de `Σ~`
é `Ω_{ℓℓ'} T_m` quando `m = m'` e `0` caso contrário. Logo o sistema
`Σ~_{S,S} w = ζ_S` **desacopla em `m`**: para cada `m`,
`(Ω_{S_m S_m} ⊗ T_m) w_{S_m} = ζ_{S_m}`, isto é
`w_{S_m} = (Ω_{S_m S_m}^{-1} ⊗ T_m^{-1}) ζ_{S_m}`. E o bloco do candidato só
enxerga `m = m0`:

```
Σ~_{(ℓ0,m0), S} w = (Ω_{ℓ0, S_{m0}} ⊗ T_{m0}) w_{S_{m0}}
                  = (Ω_{ℓ0,S_{m0}} Ω_{S_{m0}S_{m0}}^{-1} ⊗ T_{m0} T_{m0}^{-1}) ζ_{S_{m0}}
                  = (a' ⊗ I_{N_J}) ζ_{S_{m0}} = Σ_ℓ a_ℓ ζ_{ℓ m0} .
```

`T_{m0}` cancela. ∎

O que cancela é tudo que vem da base: o nível `J`, o número `N_J` de wavelets,
a família de wavelets, a condição de número de `T_{m0}` e, com ela, as
densidades marginais das moduladoras. O que sobra é `Ω` e o **ângulo entre as
componentes verdadeiras que dividem a moduladora `m0`**.

### 2.3 Resposta à conjectura, e o que ela custa

A conjectura do chat principal era que a condição se reduzisse a `Ω` sozinha.
**Não se reduz.** Sobra `V^{(m0)}`, que é propriedade da verdade e não do
desenho. Mas a redução tem consequências que valem por si:

- **Cotas.** `sqrt(a' V a) ≤ ||a||_1`, com igualdade quando as direções são
  colineares e de mesmo sinal; `= ||a||_2` quando são ortonormais. Como
  `||a||_1` é exatamente a condição forte de Zhao e Yu (2006) para o LASSO
  ordinário no desenho `Ω`, **a condição em grupos é uniformemente mais fraca
  que a do LASSO**, com ganho de até `sqrt(|S_{m0}|)`.
- **Uma condição suficiente que é só sobre `Ω`:** `max_{(ℓ0,m0) ∉ S}
  ||Ω_{ℓ0,S_{m0}} Ω_{S_{m0}S_{m0}}^{-1}||_1 < 1`. É a condição de Zhao e Yu
  aplicada, moduladora por moduladora, ao suporte `S_{m0}`. A conjectura vale,
  portanto, como **suficiência**, não como equivalência.
- **O modo de falhar fica legível.** `IRR` só passa de `1` quando `||a||_1 > 1`
  **e** as componentes que dividem a moduladora são quase proporcionais. Em
  palavras: a condição falha quando duas covariáveis lineares correlacionadas
  têm o efeito modulado pela mesma `U_m` com a mesma forma. Isso é uma hipótese
  que um referee entende e que um aplicador pode checar.

### 2.4 Onde a redução quebra

Três lugares, os dois primeiros com consequência.

- **(BD) falha.** Se as moduladoras se correlacionam através da base, `T` tem
  blocos cruzados e o sistema não desacopla: o `T_{m0}` deixa de cancelar e a
  base volta a entrar. Isso é estritamente mais forte que a hipótese de E1.4:
  a Proposição 3 só pede `λ_min(Σ_Ψ) ≥ c_U`, que permite `U_m` dependentes
  (a Observação "Joint versus marginal density" de E1.4 exclui apenas
  concurvidade exata). A redução exige independência; a Gram positiva-definida,
  não.
- **`E(XX' | U)` não constante.** É o caso geral que D13 escolheu assumir. Sem
  a fatoração de Kronecker não há redução nenhuma; o que se mede na §3.4 é
  quanto ela erra, não se ela vale.
- **O suporte.** A prova usa só (BD); a redução não pede que `S` seja um
  produto `S_ℓ × S_m`.

---

## 3. Frente (ii): os números

Tudo em `check/06a-irrepresentabilidade.R`. Base: Daublets com filtro 8,
periodizada, `j_0 = 0`, `ε = 0` (a base em que o artigo enuncia, D23/D26). As
componentes verdadeiras e a estrutura ativa vêm de `wafc/R/dgp.R`.

### 3.1 A redução, conferida

| verificação | resultado |
|---|---|
| `Σ = Π(Ω ⊗ I_K)Π'` sob `U` uniforme, por quadratura `400²` | erro `3.0e-05` (quadratura) |
| colunas não penalizadas ortogonais às penalizadas, `U` uniforme | `4.8e-14` |
| `IRR` da `Σ` cheia contra `sqrt(a'Va)` (cenário suave denso, `p = 4`, `Ω` adversária) | diferença `1.0e-17` |
| `IRR` contra `J = 2, 3, 4, 5, 6, 7` (`N_J` de 3 a 127) | `0.8733`, `0.8455`, `0.8456`, `0.8456`, `0.8456`, `0.8456` |

A estabilidade em `J` é o ponto: `N_J` passa de 3 para 127 e `IRR` se move
`0.028`, todo esse movimento entre `J = 2` e `J = 3`, onde a projeção grosseira
ainda distorce os ângulos. A partir de `J = 4` a condição é numericamente
independente do nível.

Com matrizes aleatórias (`p = 6`, `q = 3`, `N_J = 4`, 300 sorteios de `Ω`, de
`T` e do suporte):

| `T` | `|IRR_Σ − sqrt(a'Va)|` |
|---|---|
| bloco-diagonal na moduladora (arbitrário dentro do bloco) | máximo `4.4e-16` |
| cheia (moduladoras acopladas) | mediana `0.111`, q90 `0.405`, máximo `1.020`; abaixo de `1e-8` em 8% dos sorteios |

Isto é (BD) sendo necessária e suficiente na prática, e não apenas suficiente:
com `T` cheia a redução erra por ordem `0.1`, o que é a ordem da própria
margem.

### 3.2 O ganho do agrupamento sobre o LASSO

Desenho adversário de Zhao e Yu: `X_1 ≡ 1`; `X_2, X_3` ativos, correlação
`0.2`; `X_4` inativo, correlação `r` com os dois. `cos` é o cosseno entre as
duas componentes ativas.

| `r` | `cos = 0` | `cos = 0.3` | `cos = 0.6` | `cos = 0.9` | `‖a‖_1` (LASSO) |
|---|---|---|---|---|---|
| 0.5 | 0.589 | 0.672 | 0.745 | 0.812 | 0.833 |
| 0.6 | 0.707 | 0.806 | 0.894 | 0.975 | 1.000 |
| 0.7 | **0.825** | 0.941 | **1.044** | 1.137 | **1.167** |
| 0.8 | 0.943 | 1.075 | 1.193 | 1.300 | 1.333 |

Em `r = 0.7` o LASSO falha (`1.167 > 1`) e a variante em grupos vale
(`0.825 < 1`) quando as duas componentes são ortogonais. O ganho é exatamente
`‖a‖_1 / sqrt(a'Va)`, no máximo `sqrt(|S_m|)`.

O cosseno não é hipótese abstrata: é a correlação `L_2` entre as componentes.
Nos cenários de `dgp.R` (`J = 5`):

| cenário | pares |
|---|---|
| `smooth` | `<sine, cubic> = +0.969`; `<sine, cosine> = 0.000`; `<cubic, cosine> = +0.051` |
| `inhomogeneous` | `<bumps, blocks> = +0.014`; `<bumps, heavisine> = −0.114`; `<blocks, heavisine> = +0.100` |

O `sine` e o `cubic` do cenário suave são 97% correlacionados: o cenário suave
é, sem querer, o pior caso para a condição, e o não homogêneo é quase o melhor.
Isso é exatamente o que se espera — componentes espacialmente não homogêneas
são quase ortogonais entre si.

### 3.3 A margem nos cenários, e como ela degrada

`η = 1 − max_b IRR(b)`, `U` uniforme, `E(XX'|U) = Ω`, `J = 4`, `p = 8`,
`q = 2`. Três famílias de `Ω`: AR(1) de parâmetro `ρ`; equicorrelacionada `ρ`;
adversária (ativos a `0.2`, cada inativo a `ρ` de todos os ativos). Três
formas do suporte: `dgp` (a de `wafc_scenario()`, `|S_1| = 2`), `densa`
(`|S_1| = 4`, componentes distintas), `alinhada` (`|S_1| = 4`, componentes
iguais, que é `V = 11'`).

| forma | `Ω` | `ρ = 0.3` | `ρ = 0.6` | `ρ = 0.8` |
|---|---|---|---|---|
| `dgp` | qualquer | `η = 0.70` | `0.40` | `0.20` |
| `densa`, suave | equi | `0.579` | `0.388` | `0.309` |
| `densa`, não homogêneo | equi | `0.665` | `0.513` | `0.451` |
| `densa`, suave | adversária | `0.519` | `0.038` | **`−0.283`** |
| `densa`, não homogêneo | adversária | `0.618` | `0.235` | **`−0.020`** |
| `alinhada` | adversária | `0.357` | **`−0.286`** | **`−0.714`** |

Em 54 configurações a condição em grupos falha em 6 e a do LASSO em 8; as 6
são as adversárias com `ρ ≥ 0.6`, e em todas elas a do LASSO já tinha falhado.
Nos cenários como estão em `dgp.R`, com `|S_1| = 2`, a condição vale com
margem larga em toda a grade: `IRR = ρ` exatamente.

**Em `p q`.** Com `Ω` simétrica nos candidatos, `IRR` não depende de `p q` uma
vez fixado o suporte: de `p = 4` a `p = 40` (`p q` de 8 a 80) o valor é
constante, e só o número de candidatos cresce. Com `Ω` aleatória (modelo de um
fator, 200 sorteios por linha), o máximo sobre candidatos cresce, mas devagar:

| carga | `p = 4` | `6` | `10` | `20` | `40` | `80` |
|---|---|---|---|---|---|---|
| `0.90`, mediana | 0.691 | 0.652 | 0.680 | 0.713 | 0.712 | 0.715 |
| `0.90`, q90 | 0.882 | 0.806 | 0.841 | 0.871 | 0.883 | 0.891 |
| `0.97`, mediana | 0.880 | 0.701 | 0.720 | 0.738 | 0.742 | 0.758 |
| `0.97`, q90 | 0.964 | 0.881 | 0.830 | 0.915 | 0.899 | 0.922 |

Nenhuma falha em 2400 sorteios. A degradação em `p q` é lenta: a mediana sobe de
`0.69` para `0.72` quando `p` vai de 4 a 80, e o q90 fica preso em torno de
`0.88`. É o comportamento de um máximo sobre `p q − |S|` termos de cauda leve, e
não é aí que a sondagem encontra dificuldade.

**Em `J`.** `densa` e adversária com `ρ = 0.7`, `p = 8`:

| forma, cenário | `J = 2` | `3` | `4` | `5` | `6` | `7` |
|---|---|---|---|---|---|---|
| densa, suave | 1.138 | 1.122 | 1.123 | 1.123 | 1.123 | 1.123 |
| densa, não homogêneo | 0.394 | 0.944 | 0.893 | 0.866 | 0.855 | 0.866 |
| alinhada (qualquer) | 1.500 | 1.500 | 1.500 | 1.500 | 1.500 | 1.500 |

**`IRR` não degrada com `J`.** No cenário não homogêneo ele sobe de `J = 2`
para `J = 3` (a projeção em três wavelets alinha `bumps`, `blocks` e
`heavisine` muito mais do que elas são de fato) e depois estabiliza. Esta é a
boa notícia numérica mais relevante da sondagem: a resolução do sieve, que é a
fonte da alta dimensão, **não é** a fonte da dificuldade de seleção.

### 3.4 Os obstáculos, medidos

`p = 4`, `q = 2`, `J = 3`, suporte `{(1,1),(2,1),(3,1),(1,2)}`, candidato
`(4,1)`, `Ω` adversária `r = 0.7`, `E X_ℓ = 0.5`; `ζ` vindo do **oráculo
populacional de cada caso**, não da projeção em Lebesgue. `cruzM` é o maior
`|T|` entre moduladoras distintas, relativo à diagonal; `cruz_CS` é o maior
`|Σ|` entre coluna não penalizada e coluna de `S`; `vazamento` é
`‖θ*‖` do maior bloco estruturalmente nulo dividido pela do menor bloco ativo.

| caso | `IRR_Σ` | `IRR_Ω` | erro | `cruzM` | `cruz_CS` | vazamento |
|---|---|---|---|---|---|---|
| referência: `U` unif., `X ⊥ U` | 0.7922 | 0.7922 | `−8e−14` | `1e−17` | 0.000 | `1e−12` |
| `U_2 ~ Beta(2,3)` | 0.7922 | 0.7922 | `1e−16` | `2e−17` | **0.444** | `7e−13` |
| `U_2 ~ Beta(0.6,0.6)` | 0.7922 | 0.7922 | `−4e−15` | `2e−17` | 0.246 | `5e−13` |
| `U_1,U_2` cópula 0.3 | 0.7924 | 0.7924 | `−6e−15` | 0.100 | 0.000 | 0.0029 |
| `U_1,U_2` cópula 0.6 | 0.7931 | 0.7931 | `−3e−15` | 0.295 | 0.000 | 0.0067 |
| `X_2` depende de `U_1`, `τ = 0.3` | 0.7456 | 0.7808 | **−0.035** | `1e−17` | 0.207 | 0.019 |
| `X_2` depende de `U_1`, `τ = 0.6` | 0.6724 | 0.7456 | **−0.073** | `1e−17` | 0.467 | 0.029 |
| `X_2` depende de `U_1`, `τ = 0.9` | 0.5875 | 0.7194 | **−0.132** | `1e−17` | 0.784 | 0.066 |

Leitura, linha a linha:

- **Densidade marginal não uniforme não faz nada.** `cruz_CS` sobe a `0.44` —
  as colunas não penalizadas deixam de ser ortogonais às penalizadas — e mesmo
  assim o erro da redução é `1e−16`. É o complemento de Schur da §2.1
  absorvendo exatamente a não ortogonalidade. Isso é mais forte do que a
  sondagem esperava: D11 pede densidade limitada, e a redução não precisa nem
  disso.
- **Moduladoras acopladas quebram (BD) mas não apareceram aqui.** `cruzM`
  chega a `0.295` e o erro continua em `1e−15`. Não é virtude do desenho: é
  que neste suporte, com `q = 2`, o acoplamento não alcança o candidato. Em
  suportes aleatórios com `q = 3` (§3.1) o erro mediano é `0.111`. **A
  conclusão de (BD) vem da §3.1, não desta linha.**
- **`X` dependente de `U` é o obstáculo real e mensurável.** O erro cresce
  monotonamente com a força da dependência, de `−0.035` a `−0.132`. O sinal é
  negativo: a predição por `Ω` é **conservadora** nestes cenários, isto é,
  prevê uma condição mais apertada do que a verdadeira. Não há razão para que
  o sinal seja esse em geral, e a sondagem não tentou estabelecê-lo.
- **Vazamento.** O oráculo populacional carrega massa em blocos
  estruturalmente nulos assim que as moduladoras são dependentes (até `0.7%`
  do menor bloco ativo) ou `X` depende de `U` (até `6.6%`). Com `X ⊥ U` e `U`
  de coordenadas independentes o vazamento é `1e−12`, isto é, o suporte em
  blocos do oráculo **é** o suporte estrutural. Isso importa porque a saída (a)
  promete recuperar o suporte de `θ*`, e só quando não há vazamento é que esse
  suporte é o que o artigo quer nomear.

### 3.5 Versão empírica

`Σ̂` no lugar de `Σ`, mesmo `ζ`, `D = 60`, 50 réplicas, valor populacional
`0.7922`:

| `n` | mediana | q90 | máximo | réplicas com `IRR > 1` | `λ_min(Σ̂)` mediano |
|---|---|---|---|---|---|
| 250 | 0.8135 | 0.883 | **1.050** | 1 | 0.034 |
| 500 | 0.7928 | 0.843 | 0.888 | 0 | 0.052 |
| 1000 | 0.7932 | 0.816 | 0.847 | 0 | 0.067 |
| 4000 | 0.7934 | 0.810 | 0.818 | 0 | 0.084 |

Com margem populacional `η = 0.21`, uma réplica em 50 já viola a condição em
`n = 250`. O viés do máximo empírico sobre o populacional é de `0.02` na
mediana e `0.26` no pior caso em `n = 250`, e some em `n = 1000`. Tradução
prática: **a condição é um enunciado de amostra grande nos tamanhos de E2**, e
uma margem populacional abaixo de `0.2` não se traduz em seleção correta em
`n = 250`.

---

## 4. O que a sondagem não cobre

- **Não é um teorema de seleção.** Bach (2008) trabalha com número de grupos
  fixo e `n → ∞` (as hipóteses (A1) a (A3) dele são de momento e de
  inversibilidade, sem regime assintótico em `M`). Aqui `d = p q N_J → ∞` com
  `n`, porque `N_J = 2^{J_n} − 1` cresce. O que falta é a versão em dimensão
  crescente; Nardi e Rinaldo (2008, *EJS* 2, 605--633) é o lugar para procurar
  [VERIFICAR: a forma exata da condição de seleção deles, que não foi lida].
- **Wei e Huang (2010) não usam irrepresentabilidade, e isso é informação.**
  A outra referência que o catálogo sugeria prova seleção consistente para o
  LASSO **adaptativo** em grupos, a partir da condição de Riesz esparsa (SRC)
  mais um estimador inicial consistente em zero mais uma condição sobre a menor
  norma ativa. Não há condição de irrepresentabilidade em nenhum lugar. A SRC
  deles é, essencialmente, o que a Proposição 3 de E1.4 já entrega. Ou seja:
  **existe uma rota para seleção em grupos que não passa por (a)** e que
  consome E1.4 como está, ao custo de trocar o estimador pelo adaptativo e de
  assumir uma condição de separação — que é a mesma hipótese de (c).
- **O viés de aproximação não foi tocado.** `θ*` é projeção, e o modelo é mal
  especificado dentro de cada bloco ativo por `g − Π_J g`. A testemunha
  primal-dual precisa controlar `Z'(g − Π_J g)/n`, o que é uma estimativa nova,
  não presente em E1.5 nem em E1.6. A sondagem calculou tudo no nível
  populacional e com `ζ` exato.
- **Falsos negativos.** A condição de irrepresentabilidade dá ausência de
  falsos positivos. Para `P(Ŝ = S) → 1` é preciso também uma condição de
  separação (`beta-min`), a mesma de (c). Portanto **(a) não domina (c)**: ela
  troca o limiar `t_n` por uma condição de desenho, e mantém a hipótese de
  separação.
- **`ε = 0`.** Tudo aqui é na base periodizada sem margem, que é onde o artigo
  enuncia (D23, D26). O efeito de `ε > 0` sobre a condição não foi medido; é
  uma pergunta para E2.4, que já mede `λ_min(G_ε)`.
- **A necessidade.** Mediu-se a condição forte (4) de Bach, suficiente. A fraca
  (5), necessária, difere apenas pelo `≤`; a sondagem não distinguiu as duas em
  nenhuma tabela, porque as margens medidas nunca ficaram no fio.

---

## 5. Veredito detalhado

**O que fechou.** A redução. Sob `E(XX'|U) = Ω` constante e (BD), a condição de
irrepresentabilidade em grupos do desenho de produtos é
`max sqrt(a' V^{(m)} a) < 1`, com `a` vindo só de `Ω` e `V^{(m)}` a matriz de
correlações das componentes verdadeiras que dividem `U_m`. A base, o nível `J` e
as densidades marginais cancelam. Conferido a `1e−17` no desenho do WAFC e a
`4e−16` em 300 sorteios aleatórios.

**O que isso vale.** Meia página de artigo, agora, sem depender de mais nada. É
a observação de que o desenho de produtos não paga nada, em irrepresentabilidade,
pelo fato de ser de alta dimensão em wavelets: a dificuldade de seleção é toda
do lado de `Ω` e do ângulo entre componentes. Num artigo cuja contribuição é o
desenho, isso é substantivo, e é o tipo de coisa que responde à pergunta do
referee sem gastar um teorema.

**O que não fechou, e o preço.** A saída (a) inteira. Falta o teorema em
dimensão crescente, o passo do viés, e a hipótese (BD), que é mais forte que
D11 e que D13 (que escolheu o caso geral `E(XX'|U)` não constante). Sob D13, a
redução não vale, e o erro medido cresce com a dependência (até `−0.132`). O
custo continua sendo o de `selecao-estrutura.md`: da ordem de E1.4 mais E1.5.
A sondagem retirou o **risco** (a condição não falha por construção no desenho
de produtos), não o **custo**.

**Recomendação.**

1. Incorporar a redução como observação, agora. Lugar natural: uma observação
   nova em `03-desenho-produtos.tex`, ao lado da parte (ii) da Proposição 3,
   que já tem a fatoração de Kronecker na mão. Custo: uma rodada curta.
2. Não abrir a prova de (a) antes de E2.5. Se E2.5 escolher o LASSO puro, como
   E2.2 e E2.3 apontam, (a) perde o objeto e a observação do item 1 é tudo o
   que sobra dela no artigo.
3. Se (a) for aberta depois de E2.5, abrir na forma de Wei e Huang e não na de
   Bach: LASSO adaptativo em grupos sobre a SRC que E1.4 já provou, mais
   separação. Consome E1.4 como está, e a hipótese de separação já vai estar
   escrita por E1.7c.
4. E1.7c segue como planejada e independente. Nada aqui a altera.

---

## 6. Referências

Conferidas antes de citar, como em L1, L3 e L4.

- **Bach, F. R. (2008).** Consistency of the group Lasso and multiple kernel
  learning. *Journal of Machine Learning Research* 9(40), 1179--1225.
  Conferido na página oficial do JMLR (`jmlr.org/papers/v9/bach08b.html`):
  título, autor, volume, número e páginas. **Sem DOI** (o JMLR não deposita no
  Crossref, como L1 já registrou para Xue e Qu 2012). As condições (4) e (5) e
  os Teoremas 2 e 3 estão na Seção 2.5; a observação sobre limiarizar os
  mínimos quadrados está na discussão que segue o Teorema 3.
  Chave sugerida: `bach2008consistency`.
- **Wei, F. and Huang, J. (2010).** Consistent group selection in
  high-dimensional linear regression. *Bernoulli* 16(4), 1369--1384.
  DOI `10.3150/10-BEJ252`. Conferido no Crossref (título, autores, veículo,
  volume, número, ano) e no OpenAlex (páginas); o texto foi lido no arXiv
  (1011.6161). A SRC é a equação (2.4) da Seção 2; as condições (C1) a (C4) e o
  Teorema 3.1 estão na Seção 3.
  Chave sugerida: `wei2010consistent`.
- **Zhao, P. and Yu, B. (2006).** On model selection consistency of Lasso.
  *Journal of Machine Learning Research* 7(90), 2541--2563. Conferido na página
  oficial do JMLR (`jmlr.org/papers/v7/zhao06a.html`). **Sem DOI.**
  Chave sugerida: `zhao2006model`.
- **Nardi, Y. and Rinaldo, A. (2008).** On the asymptotic properties of the
  group lasso estimator for linear models. *Electronic Journal of Statistics*
  2, 605--633. DOI `10.1214/08-EJS200`. Conferido no Crossref (título, autores,
  veículo, volume, ano) e no Semantic Scholar (páginas). **O texto não foi
  lido**: a afirmação de que eles tratam seleção em grupos em dimensão
  crescente é de memória e fica `[VERIFICAR]`.
  Chave sugerida: `nardi2008asymptotic`.

Nenhuma delas está em `docs/referencias-verificadas.bib`; a inclusão é do chat
principal (o catálogo de `TAREFA.md` não dá o `.bib` a esta tarefa).
