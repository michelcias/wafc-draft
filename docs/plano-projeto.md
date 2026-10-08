# Plano do projeto

**Objetivo.** Um estimador para modelos de regressão com coeficientes
funcionais aditivos, `Y = Σ_ℓ β_ℓ(U) X_ℓ + ε` com `β_ℓ(u) = c_ℓ + Σ_m
g_{ℓm}(u_m)`, em que cada `g_{ℓm}` é expandida numa base ortonormal de
wavelets e todos os coeficientes são estimados por LASSO num único problema
convexo. Três produtos: a teoria (desigualdade oráculo no desenho de
produtos, taxas em Besov, adaptação), o código do método na pasta `wafc/`
deste repositório (com testes; empacotamento decidido depois de E2.5),
e um artigo num periódico Q1 de Statistics and Probability (alvo ratificado
em D5: *Statistica Sinica*, confirmada Q1 no SCImago; reserva: EJS).

**Notação.** Este arquivo foi escrito antes de E1.1 e usava `j` e `k` para as
covariáveis. Desde 2026-09-18 vale `notacao.md`, congelada: covariável linear
`X_ℓ`, moduladora `U_m`, wavelet `ψ_{jk}`, coeficiente `θ_{ℓm,jk}` (D9 a D12),
e centralização `∫_0^1 g_{ℓm} = 0` (D22). Os trechos abaixo foram passados
para ela em 2026-09-19.

**Como ler.** O trabalho está dividido em etapas `E0` a `E7`, mais duas de
literatura (`L1`, `L2`), desenhadas para serem o mais independentes possível.
Cada etapa diz de quem depende, o que entrega, e como se verifica que acabou.
Etapas sem dependência entre si podem andar em paralelo ou em qualquer ordem,
em chats separados (`instrucoes.md`, §7). Os critérios de saída estão
escritos para serem verificáveis, não argumentados.

**Regra de ouro.** Nenhuma etapa fecha sem o critério de saída satisfeito; uma
etapa pode *começar* antes da anterior fechar quando a dependência é só
parcial (a tabela abaixo diz quando).

---

## Mapa de dependências

```
E0 infraestrutura ─────────────────────────────────────────────────────────┐
L1 verificação bibliográfica ──┐                                           │
L2 busca de novidade ──────────┼──► E5a manuscrito: seções 1 a 4 ──► E5b ──┤
                               │         ▲                                 │
E1 teoria ─────────────────────┘─────────┘                                 ├──► E7 fechamento
   (E1.1 notação primeiro)                                                 │
E2 código em wafc/ (go/no-go) ──► E3 consolidação ──┬──► E4 simulação ──┤
                                                    └──► E6 aplicação ──┘
```

| Etapa | Depende de | Pode andar junto com |
|---|---|---|
| E0 infraestrutura | nada | tudo |
| L1 verificação bibliográfica | nada | tudo |
| L2 busca de novidade | nada | tudo; deve fechar antes de E5a |
| E1 teoria | E1.1 (notação) | E0, E2, L1, L2 |
| E2 código em `wafc/` | E1.1 (notação) e o enunciado de E1.2 (identificabilidade) | E1.3 a E1.7 |
| E3 consolidação do código | E2.5 (variante e sintonia decididas) | E1 provas, E5a |
| E4 simulação | E3.1 (interface congelada) | E6, E5a |
| E5a manuscrito, seções 1 a 4 | enunciados de E1.2 a E1.6, L2 | E2, E3, E4, E6 |
| E5b manuscrito, resultados | E4 e E6 | |
| E6 aplicação | E3.1 e E6.1 | E4 |
| E7 fechamento | E4, E5, E6, E3.4 | |

---

## E0. Infraestrutura e decisões pequenas

Tudo aqui é barato e independente do resto; o valor está em não deixar essas
decisões para a hora errada.

### E0.1 Nomes

| O quê | Nome | Estado (2026-09-30) |
|---|---|---|
| método | WAFC, *wavelet additive functional coefficients* | ratificado (D8, 2026-09-19) |
| funções | `wafc()`, `cv.wafc()`, classe `"wafc"`, na pasta `wafc/` | D4 fechada: pasta dedicada, não o `WaveBased` |
| repositório de rascunho | `michelcias/wafc-draft` (este) | publicado (privado) em 2026-09-18 |
| compêndio | pasta `wafc-studies/` deste repositório; `michelcias/wafc-studies` no fim | nasce em E4.1 como pasta (D60); o repositório é criado pelo autor no fim, por cópia |

### E0.2 Onde o código vive em cada momento

| Código | Onde | Por quê |
|---|---|---|
| conferência numérica de um resultado (`n ≤ 200`, denso) | `derivations/check/` aqui | é parte da prova, viaja com ela |
| o método: `wafc()`, `cv.wafc()`, métodos, testes, scripts | `wafc/` aqui (`R/`, `tests/`, `scripts/`) | decisão do autor (D4): pasta dedicada dentro do `wafc-draft`; o `WaveBased` instalado é só dependência para as bases; empacotar (pacote próprio ou `WaveBased`) decide-se depois de E2.5 |
| estudo de simulação e aplicação | pasta `wafc-studies/` aqui; repositório próprio só no fim (D60) | pipeline pesado com cache; o `wall` é o molde; o repositório final reproduz só o que está no manuscrito e no supp |

### E0.3 Template e instruções da revista

- **Feito:** `manuscript/ejs-template/` (imsart 2025/03/18, `ejsv2`) compila.
- **Feito (2026-09-18/19):** instruções da *Statistica Sinica* transcritas
  em `docs/ss-instrucoes-autores.md` (teto de 40 páginas em espaço duplo,
  referências incluídas; tipo de revisão não declarado), template em
  `manuscript/ss-template/` compilando, e Q1 confirmado pelo autor no
  SCImago.

### E0.4 Ferramentas

Presentes: R 4.6.1, `glmnet`, `grpreg`, `gglasso`, `mgcv`, `Matrix`,
`bench`, `testthat`, `devtools`, `roxygen2`, `pkgdown`, `renv`, `latexmk`.
`WaveBased` 2.6-0 instalado de `../../WaveBased` em 2026-09-18; `sparsegl`,
`VCBART`, `quadprog` e `gh` instalados depois. A lista em dia está em
`CONTINUAR.md` §2.

### E0.5 Inventário

Feito em `inventario-codigo.md`: o desenho por blocos do `wall()` é a
referência de leitura para `wafc_design()`; o `WaveBased` instalado fornece
`wbasis()` e `wtable()`.

### E0.6 Convenções do código em `wafc/`, fixadas agora

- `wafc()` tem a mesma assinatura de `wall()` até onde fizer sentido (`J`,
  `j0`, `family`, `filter.size`, `boundary`, `rescale`, `eps`, `use.table`,
  `sparse`, `lambda`, `standardize`, `weights`), mais `x` (lineares) e `u`
  (moduladoras) separados, e `penalty = c("lasso", "sglasso")` se E2 mantiver
  a variante.
- Organização de pacote sem ser pacote: `wafc/R/` (um arquivo por tema,
  roxygen nos cabeçalhos desde já), `wafc/tests/` (`testthat`),
  `wafc/scripts/` (numerados), `wafc/R/load.R` que carrega tudo e declara as
  dependências. O compêndio fixa o código por commit deste repositório.
- Estilo do `WaveBased`: `snake_case` nas funções internas, `.` só nos
  métodos S3; erro informativo em argumento desconhecido.

**Critério de saída de E0:** D5 e D8 ratificados; instruções da SS
transcritas e template compilando; `WaveBased` instalado; repositório
publicado no GitHub.

---

## L1. Verificação bibliográfica

Conferir no Crossref (título, autores, ano, veículo, volume, páginas, DOI)
cada entrada de `docs/referencias-verificadas.bib` e cada `[VERIFICAR]` de
`docs/literatura.md`; copiar de `wall-manuscript/manuscript/theo/references_theo_1.bib`
as entradas `no wall`. Entregável: `.bib` sem notas `[L1: confirmar]`;
`literatura.md` sem `[VERIFICAR]` (ou com o motivo de não ter achado).

**Critério de saída:** todas as entradas com status `verificado`.

## L2. Busca de novidade

As quatro buscas de `literatura.md`, "Buscas pendentes". Entregável:
`docs/busca-novidade.md` com uma tabela por busca (trabalho, o que faz, o que
não faz, ameaça à novidade) e um veredito por contribuição de
`alvo-revista.md` §4. Ler o PDF de Sardy & Ma (2024) inteiro e registrar o
que a teoria deles cobre.

**Critério de saída:** veredito escrito; lista de contribuições reescrita se
preciso (decisão do autor).

---

## E1. Teoria

Cada resultado é um arquivo em `derivations/`, com conferência numérica em
`derivations/check/` **antes** da prova (`instrucoes.md`, §5). O molde é o
WALL teórico; cada arquivo diz o que foi transposto e o que é novo.

### E1.1 Congelar a notação

**Fechada em 2026-09-18.** `notacao.md` saiu do estado de esboço: os cinco
pontos da §6 decididos pelo autor (D9 a D12), `derivations/macros.tex`
reescrito. Emendas posteriores: D22 (centralização de Lebesgue) e D24
(estilo da revista), nas §7 e §8 do arquivo.

### E1.2 `01-identificabilidade.md`

**Fechada em 2026-09-18** (Proposição 1, Lema 1). Enunciado: sob
`E(XX' | U)` não singular q.c., densidade conjunta de `U` positiva em
`[0,1]^q` e `∫_0^1 g_{ℓm} = 0` (D22), o vetor `(c_ℓ, g_{ℓm})` é
identificado pela função de regressão. Mostrar que a base periódica com `j0
= 0` sem a função de escala impõe a restrição no nível da base, e o que
muda com `boundary = "interval"` (restrição explícita). Conferência: em `n`
pequeno, a matriz de desenho tem posto cheio e a projeção recupera
coeficientes conhecidos.

### E1.3 `02-aproximacao-besov.tex`

**Fechada em 2026-09-18** (Lemas 2 e 3, Corolário 1, Proposição 2).
Lema de aproximação: `g ∈ B^s_{π,r}[0,1]`, `s > 1/π`, base com `N > s`
momentos nulos; `‖g − Π_J g‖_{L_2(P_U)} ≤ C 2^{−J s'}` com densidade de `U`
limitada; e a versão em `L_∞` se precisar para a compatibilidade. Transposto
do WALL (Passo 1); registrar a diferença por causa da periodização (o WALL
já a trata?). Conferência: erro de projeção medido em `bumps` e `sin` contra
a curva `2^{−J s'}`.

### E1.4 `03-desenho-produtos.tex`

**Fechada em 2026-09-19; o texto abaixo é o que ela virou depois de L2
(proposta L2b, ratificada).** (i) O caso `X ⊥ U`, em que `Σ` fatora como
`E(XX') ⊗ Σ_ψ`, **não é novo**: é a eq. (1.8) e o Lema 1 de Klopp & Pensky
(2015), a recordar e citar, não a provar. (ii) O conteúdo novo é o termo
cruzado entre moduladoras, `E{X_ℓ X_{ℓ'} ψ_{jk}(U_m) ψ_{j'k'}(U_{m'})}` com
`m ≠ m'`, que não é produto de Kronecker, e o caso geral sob
`λ_min(E(XX' | U)) ≥ κ_1` q.c. com a densidade conjunta de `U` limitada por
baixo, que saiu **mais forte do que o previsto**: autovalor mínimo cheio,
`λ_min(Σ) ≥ κ_1 c_U`, sem condição de cone. (iii) Versão empírica por
concentração, com `p q 2^J log(p q 2^J)/n → 0`. Conferência feita em
`n = 200`, `J = 3`, nos dois casos.

Sob D23, a margem do reescalonamento muda a constante para
`c_U λ_min(G_eps)`, com `G_eps` a Gram da base restrita ao suporte; a
emenda é E1.3b e a medida é de E2.4.

### E1.5 `04-oraculo.tex`

Desigualdade oráculo para o LASSO com perda quadrática, erro sub-gaussiano,
**Fechada em 2026-09-19** (Lemas 4 a 7, Teorema 1, Corolários 2 e 3), e sem
condição de cone: a perfilagem dos níveis não penalizados dá
`λ_min(B̃'B̃/n) ≥ λ_min(Σ̂)`. `c_ℓ` não penalizados;
`λ ≍ σ σ̂_max sqrt(log(p q N_J)/n)`, com `σ̂_max = max_a sqrt(Σ̂_aa)`, que é a
leitura certa de `‖Z‖_max` (a literal é de ordem `2^{J/2}` e destruiria a
taxa). Segue Bühlmann & van de Geer (2011), **Corolário 6.1** (taxa lenta) e
**Teorema 6.2** (versão com viés), com as coordenadas não penalizadas na
§6.9 e não na §6.2.3 — numeração conferida em L3.
Conferência: em simulação com `s_0` pequeno, o erro de predição escala como
`λ² s_0` ao variar `n`.

### E1.6 `05-taxas.tex`

Taxas: com `J_n ≍ log_2 n/(2s'+1)`, erro de predição `n^{−2s'/(2s'+1)}` a
menos de logaritmos; erro `L_2` de cada `ĝ_{ℓm}` pela condição de desenho
(a norma de predição controla cada componente); corolário de
compressibilidade (weak-`ℓ_τ`): o LASSO adapta à dimensão efetiva, com a
taxa melhor que a do sieve cheio, transposto do WALL (trilha rápida).
**Fechada em 2026-09-19** (Lemas 8 e 9, Proposição 4, Teorema 2, Corolários
4 e 5), com um achado que não estava previsto: a hipótese de Besov **já
implica** weak-`ℓ_τ`, então a compressibilidade não custa hipótese. O
enunciado principal do artigo sai daqui e é o Corolário 5 (D16).

### E1.7 Seleção de estrutura

**Decidida (D28, D32) e fechada em 2026-09-20:** E1.7c entregou a seleção
por limiarização (`06-selecao-limiar.tex`, Lema 13 e Corolário 8), e E1.7a
deu veredito de escopo reduzido para a saída (a)
(`06a-sondagem-irrepresentabilidade.md`). O texto abaixo é o registro da
pergunta original. Era: **o escopo desta etapa está em decisão (L2f).** O plano original era
seleção consistente de `{(ℓ,m): g_{ℓm} ≢ 0}` com sparse group LASSO
adaptativo, seguindo Wei, Huang & Li (2011) e Huang, Horowitz & Wei (2010).
L2 mostrou que a ideia não é nova, e as três saídas estão escritas em
[`selecao-estrutura.md`](selecao-estrutura.md): a propriedade oráculo em
grupos, a variante sem teorema, e a seleção por limiarização como corolário
do Corolário 4. O autor decide.

### E1.3b (emenda) A margem e o lema de extensão

**Fechada em 2026-09-19.** Hipótese 5 e Lema 10 em
`02-aproximacao-besov.tex`: com margem fixa, a componente pode ser
substituída por uma extensão que emende, e o viés volta a `O(2^{−2Js'})`. O
preço é `C(eps) ≍ eps^{−(s−1/π)}`, exato, e a margem **não pode encolher com
`J`**. A tensão com a Gram restrita e a decisão de manter a teoria em
`eps = 0` estão em D26.

### E1.8 A rota do intervalo (fechada em 2026-09-19)

`07-rota-intervalo.tex`: a teoria sem periodicidade, pronta para o caso de
um referee pedir — o que transfere verbatim para a base CDV, o que precisa
ser conferido nela, o enunciado verdadeiro para `eps > 0`, e a análise de
`s'` por componente nos três regimes. **Não vai ao manuscrito agora** (D23,
D26).

### E1.13 A taxa lenta em blocos (fechada em 2026-10-03, D47(c))

O enunciado sem condição de desenho (a Proposição 4 de E1.6, Proposition 3
do manuscrito) na versão em blocos, como Corolário 14 do `08-blocos.tex`,
a partir do Teorema 3(i). Fica no corpo do artigo em `k = 4`.

### E1.14 A taxa lenta em todo `s' > 0` (fechada em 2026-10-04)

A Proposição 4 de E1.6 e o Corolário 14 de E1.12 completados nos regimes
que D49 excluiu (`s' ≤ 1/4` com `s > 1/2`, `s' ≤ s/2` com `s < 1/2`) e na
fronteira `s = 1/2`, parando `J_n` no maior nível admissível (pergunta 45
do `ESTADO.md`).

### E1.15 O comparador truncado na taxa lenta (fechada em 2026-10-05)

Proposição 8 no `05` e Corolário 15 no `08`: em `s < 1/2`, a taxa lenta pelo
comparador truncado alcança a taxa do Theorem 1 sem condição de desenho; em
`s ≥ 1/2` nenhum comparador melhora. Números no `ESTADO.md` §2.

### E1.9 Cota inferior para `q ≥ 2` (fechada em 2026-10-05)

**Fechada:** `09-cota-inferior.tex` (Lema 17, Teorema 5, Corolário 16), a
cota no nível da taxa com constante proporcional a `pq`, para qualquer `q`
e `X` dependente de `U`; números no `ESTADO.md` §2.

**Catalogada em 2026-10-05, com escopo reduzido:** no nível da taxa, pela
subclasse em que só uma componente é não nula (regressão não paramétrica
com desenho aleatório), com o lema de Assouad; a constante em `pq` é
opcional; o logaritmo em `π < 2` fica fora. Com `p` e `q` fixos a
dificuldade da cota de K&P, que é de alta dimensão, não se aplica. O texto
abaixo é o registro de quando ela não estava aberta.

Klopp & Pensky (2015) têm a cota inferior minimax para `q = 1` com
`X ⊥ U`; para o desenho aditivo com `q ≥ 2` e `X` dependente de `U` **não
há nenhuma**, e o `05-taxas.tex` diz isso explicitamente, comparando com a
referência do modelo de sequência. Construí-la é trabalho do porte de E1.4
mais E1.5 somadas e **não estava no plano**.

**A decisão de abrir ou não fica para depois dos resultados** (E2.5 e E4):
se a evidência numérica sustentar o artigo, a cota inferior é cortesia que
um referee da SS pode pedir e que se responde em revisão; se a contribuição
numérica ficar fraca, ela vira o peso que falta. Sob D18 o artigo se declara
extensão de quem tem a dele, o que torna a ausência mais visível — é o risco
assumido, registrado na pergunta 20 do `ESTADO.md`.

### E1.12 A teoria em blocos (fechada em 2026-10-03, D44)

**Fechada:** `08-blocos.tex` (Proposição 7 a Corolário 11) e o adendo do
`06` (Lema 16, Corolários 12 e 13), com `check/08-blocos.R` imprimindo `OK`;
números no `ESTADO.md` §2, pendências na pergunta 42. O texto abaixo é o
registro do catálogo.

Com o block LASSO balanceado como o WAFC (D44), a sondagem de E1.11
(`08a-sondagem-blocos.md`, com os adendos das §§ 11 e 12) vira resultado
numerado em `08-blocos.tex`: a calibração em blocos, o oráculo sem cone com
pesos de razão limitada, o risco ideal por pedaços sob Besov, a taxa do
espaço de aproximação e as taxas adaptativas (o análogo do Corolário 5, sem
logaritmo em `π ≥ 2`). O `06-selecao-limiar.tex` recebe o Corolário 8 em
blocos e a cota de risco do estimador limiarizado (D45). Os resultados do
LASSO (`04`, `05`) ficam como estão, para a opção coordenada (D43). Custo
estimado em E1.11: ~5 páginas e 2 a 3 dias. Entregável e arquivos no
catálogo do `TAREFA.md`.

**Critério de saída de E1:** E1.2 a E1.6 com prova e script de conferência
imprimindo `OK`; hipóteses numeradas e congeladas; o teorema principal
enunciado na forma que vai ao manuscrito. **Atingido em 2026-09-19**; E1.7 é
condicional e E1.8 e E1.9 são opcionais.

---

## E2. O código em `wafc/` e o go/no-go

Vive em `wafc/` (`R/`, `tests/`, `scripts/`; ver `wafc/README.md`). R sobre
`WaveBased` (bases) e `glmnet`. `n ≤ 2000`.

### E2.1 Cenários e desenho

`R/dgp.R`, `R/design.R`, `tests/test-design.R`, `scripts/01-smoke.R`. O
desenho é construído com `wbasis()` do `WaveBased` bloco a bloco e
multiplicado por `X_ℓ`; colunas nomeadas; `penalty.factor` zero nos `c_ℓ`.
Conferência (e primeiro teste): com `g_{ℓm}` na base (`θ*` conhecido) e sem
ruído, o LASSO com `λ → 0` recupera `θ*`. **Fechada em 2026-09-19**, com a
ordem de colunas de D12 e o descarte de `φ_{00}`; a emenda E2.1b desacoplou
a margem `eps` do nível `J`.

### E2.2 Variantes do estimador

`R/fit.R`, `R/reconstruct.R`, `tests/test-fit.R`: `wafc()` com LASSO
(`glmnet`, gaussiano) e sparse group LASSO (`sparsegl`; grupo por bloco
`(ℓ, m)`); `predict`, `coef`; reconstrução de `ĝ_{ℓm}` e `β̂_ℓ`. Registrar em
`CONTINUAR.md` o pacote escolhido. **Fechada em 2026-09-19** (D17); o
`sparsegl` ficou registrado como exigido só por `penalty = "sglasso"`.

### E2.3 Sintonia

`R/tune.R`, `tests/test-tune.R`: `cv.wafc()` sobre `(J, λ)` como o
`cv.wall()`; BIC e EBIC com graus de liberdade = número de não nulos.
Comparar as três regras nos cenários de E2.1 em `n ∈ {250, 1000}`.
**Fechada em 2026-09-19** (D19, D20: `cv.min` é o padrão).

### E2.4 Piloto

`R/competitors.R` e `scripts/04-pilot.R`: WAFC (LASSO; grupos) contra os
**concorrentes mínimos de L2d** (proposta ratificada): `mgcv::gam` com
`s(u_m, by = x_ℓ)`; B-splines mais group LASSO (`grpreg`); o **spline
adaptativo** de Wang, Jiang & Liu (2024), que é o concorrente que adapta
localmente; o **block LASSO de Klopp & Pensky** no mesmo desenho, que sai
com `grpreg`/`gglasso` e é a pergunta que o posicionamento D18 convida; o
**VCBART**, não aditivo e com várias moduladoras; e a regressão linear
oráculo. Cenários: (a) todas `g_{ℓm}` suaves; (b) não homogêneas, com as
funções de teste de **Donoho & Johnstone** (bumps, blocks, heavisine) e
SNR 3, como em Sardy & Ma; (c) mistura com metade das `g_{ℓm}` nulas
(`p = 4`, `q = 4`). `n ∈ {250, 500, 1000}`; `p q` até 16 no piloto.
Métricas: ISE de cada `ĝ_{ℓm}` no interior, RMSE de predição fora da
amostra, acerto de estrutura por bloco `(ℓ,m)`, segundos. 50 réplicas.

Além disso, e vindo das rodadas de 2026-09-19: medir **`λ_min(G_eps)` e o
ISE contra `eps`** (D23, com `eps ∈ {0, 2^{−J−1}, 1.9^{−J}}`), que é o que
decide o padrão da margem e conserta o `eps` que hoje varia com `J`; e
incluir o **QUT** de Giacobino et al. (2017) como regra de `λ` sem `σ`
(proposta L2c, ratificada), ao lado das cinco regras que E2.3 já compara.

**Fechada em 2026-09-20 (E2.4) e 2026-09-21 (emendas E2.4b):** números no
`ESTADO.md` §2. O que ela deixou para E2.5 está no `ESTADO.md` §5.

**E2.4c e E2.5a (fechadas em 2026-09-30):** a repetição da parte
`competitors` com a grade `2:8` (D34), o `bsgl` na mesma grade (D36), o
`gam` em dimensão casada por `bam` (`gam.matched`) e o cenário `uneven`
(D30). Resultado: **no-go para a variante LASSO** pelo critério abaixo; o
rumo está na pergunta 33 do `ESTADO.md`.

### E2.5 Go/no-go e variante principal

**Go** se o WAFC vence os competidores no cenário (b) por margem clara e
não perde no (a) por mais de um fator a fixar (proposta: 1,5 em ISE). Se
não: antes de mudar de rumo, testar `boundary = "interval"`, pesos
adaptativos e limiarização em blocos. Registrar em `ESTADO.md`: variante
(LASSO ou grupos), regra de sintonia, `J` padrão, e os números. Desde D32 e
E1.7c a escolha de variante é entre o **LASSO limiarizado** e o sparse group
LASSO, e E2.5 calibra também o limiar `t_n` (pergunta 11 do `ESTADO.md`).
A grade padrão (`2:8`, D34) e a margem (`0`, D35) já estão decididas.

**Medido de 2026-09-30 a 2026-10-02 (E2.5a a E2.5h, números no
`ESTADO.md` §2):** no-go literal para o LASSO puro (E2.5a); cinco formas do
block LASSO de K&P, das quais a balanceada é a que a teoria cobre e a melhor
fora do suave (E2.5b a E2.5f); a estimação seguida de limiar, que derruba o
fator do suave abaixo de 1,5 contra o `gam.matched` (E2.5g); e o `gam`
sintonizado por REML, validação cruzada e GCV numa grade única (D41), que
não muda o veredito (E2.5h). Ficaram por testar `boundary = "interval"` e
pesos adaptativos (o primeiro pede código antes: a reparametrização do
bloco de escala, pergunta 33(f)). A decisão de rumo é do autor (perguntas
33, 38 e 41). **E2.5j** (fechada em 2026-10-03) mediu a regra de `t`: a de
um erro-padrão (`cv1se`) chega ao teto do limiar oráculo no suave, no
`uneven`, no nulo e em `n = 1000`; o `c` relativo repete o `+cv`; a porta do
QUT resolve o nulo das outras regras. Nada medido em E2.5 é descartado
(D43).

**Decidido em 2026-10-03 (D44 a D46, `ESTADO.md`): go reposicionado.** A
variante principal é o block LASSO na forma balanceada seguido do limiar
`cv1se`, com o LASSO coordenado como opção; regra de sintonia `cv.min` (D20)
na grade `2:8` (D34); o alvo continua a *Statistica Sinica*, e a resposta
da tabela de riscos (JCGS) não foi necessária. **E2 fechada.**

**Critério de saída de E2:** testes de E2.1 a E2.3 passando; tabela do
piloto; decisão registrada.

---

## E3. Consolidação do código

Depois do go de E2.5, o código de `wafc/` deixa de ser exploratório.

### E3.1 Interface congelada

**Fechada em 2026-10-03** com o estimador de D44 e D45: o block LASSO
balanceado como `penalty = "block"`, padrão, e o limiar `cv1se` na
interface; a interface foi ratificada em D48. Assinaturas de `wafc()`, `cv.wafc()`, `predict`, `coef`, `plot` fixadas e
registradas em `wafc/README.md`; nenhuma variante de E2.5 é removida sem
nova decisão (D43): as que não forem a principal ficam como opção;
`tests/` cobrindo cada função pública; `Rscript -e 'testthat::test_dir("wafc/tests")'`
limpo em menos de 60 s.

### E3.2 Gráficos e documentação

**Fechada em 2026-10-03**, com a correção do `wafc_thr_object()` de D48; números no `ESTADO.md` §2.

`plot.wafc` com painel por `(j,k)` (função reconstruída e, se houver, a
verdade) e o caminho das normas por par; `plot.cv.wafc`; cabeçalhos roxygen
completos; exemplo reproduzível em `wafc/README.md`.

### E3.4 `J` por moduladora (fechada em 2026-10-05)

O nível de cada moduladora limitado pelo número de valores distintos, para
que um bloco de moduladora discreta fique identificado (lição da E6.1b,
D57).

### E3.3 Empacotamento (decisão do autor)

Com o código testado, decidir: pacote próprio `wafc`, função dentro do
`WaveBased`, ou permanecer em `wafc/` e ser citado pelo repositório. O que a
revista pede é código disponível e citável; a decisão afeta E7 (DOI), não
E4 e E6. A pasta já está organizada como `R/` + `tests/` para que qualquer
das três saídas seja mover arquivos.

**Decidida em 2026-10-06 (D61):** `wafc/` fica onde está enquanto o
trabalho durar e vai ao `WaveBased` só quando tudo estiver pronto, antes da
E7; o compêndio extraído passa então a depender dessa versão do
`WaveBased`, fixada pelo `renv`.

**Critério de saída de E3:** interface congelada; testes limpos; um usuário
que não é o autor ajusta o exemplo do `wafc/README.md` seguindo só ele;
decisão de E3.3 registrada.

---

## E4. Estudo de simulação: o compêndio `wafc-studies`

### E4.1 Nascimento

Nos moldes do `wall`: `R/`, `scripts/` numerados, `config/*.yaml`, `renv`
fixando o `WaveBased` por commit e o código de `wafc/` pelo commit do
`wafc-draft` (submódulo ou cópia com o hash registrado em `PROVENANCE.md`),
cache por unidade retomável, semente mestra única, `INSTRUCTIONS.md` e
`CLAUDE.md` desde o primeiro commit.

**Emenda de D60 (2026-10-06):** o compêndio nasce como a pasta
`wafc-studies/` deste repositório, não como repositório próprio. Enquanto
for pasta, lê o código de `../wafc/` diretamente; na extração, o código vai
junto como a versão do `WaveBased` que o contiver (D61), com o hash do `wafc-draft` no
`PROVENANCE.md`. O repositório é criado pelo autor no fim, **por cópia, sem
o histórico** (um `git subtree split` levaria as tentativas), e reproduz só
o que está no manuscrito e no supp. Por isso: a documentação da pasta é
escrita desde o início como a do estudo final, sem narrar o caminho; o
piloto (E4.3), as medições e as configurações descartadas ficam fora da
pasta (`wafc/scripts/` e os documentos de trabalho), ou saem dela antes da
cópia; o cache e os produtos pesados não são versionados.

### E4.2 Desenho

| Fator | Níveis |
|---|---|
| regularidade das `g_{ℓm}` | suaves (`s' = 3/2`); **suaves de curvatura desigual** (`C^∞`, dentro da hipótese de Xue & Yang, mas com escala variando: gaussiana estreita, `doppler` truncado — D30); não homogêneas (`s' = 1/2`); mistura (D27) |
| esparsidade | todas ativas; metade nulas; esparso (`p q = 50`, 6 ativas) |
| `n` | 250, 500, 1000, 2000 |
| `(p, q)` | (2, 2), (4, 4), (10, 5) |
| dependência `X`–`U` | independentes; correlacionados |
| ruído | dois níveis de razão sinal-ruído |
| métodos | WAFC (variante de E2.5); `mgcv`; B-splines + group LASSO; spline adaptativo (Wang, Jiang & Liu 2024); block LASSO de Klopp & Pensky no mesmo desenho; VCBART; linear oráculo (L2d) |
| réplicas | fixadas pelo piloto (E4.3) |

Métricas: ISE por função e total; RMSE de predição; suporte; tempo; **a
predição e o ISE relatados em separado** (D30: ganhar em um e perder no
outro é informação). Cada concorrente é sintonizado nos termos dele — o
spline com número de nós por BIC, como em Xue & Yang, e o `mgcv` com REML —,
senão um ganho do WAFC é sobre sintonia e não sobre método. A
tabela principal do artigo é regularidade × método em `n = 1000`; a figura
principal é `ĝ_{ℓm}` sobreposta à verdade em bumps e blocks, WAFC contra
`mgcv`.

**Proposta revista de 2026-10-06** (chat principal, depois de D44 a D59;
**ratificada pelo autor como D60**, com 100 réplicas, os três braços e a
escala sob a condição de custo abaixo). A tabela acima é o desenho
original e fica como registro.

- **Princípio:** o corpo do artigo leva uma tabela (regularidade × método em
  `n = 1000`, com uma coluna de estrutura) e uma figura (`ĝ` sobre a verdade
  em `bumps` e `blocks`, WAFC contra `gam.reml`); o resto vai ao supp (D21).
  A E2.5 serviu de piloto: grade de `J` `2:8` (D34), 10 dobras, tabela da
  base (D31).
- **Métodos:** WAFC `+cv1se` (o principal, D44, D45) e `+cv` (a predição,
  como na aplicação, D57); o LASSO coordenado só no supp (D43); `gam.reml` e
  `gam.gcv` (D46; o `gam.gcv` fora da `mixed`, D41); `bsgl`; `aspline`;
  `klopp` (o block LASSO de K&P na forma deles); o linear oráculo; e o
  **VCBART por curiosidade**, porque é bayesiano e a decisão de usá-lo no
  artigo vem depois (autor, 2026-10-06). Ficam no código e fora do estudo
  (D43): o sparse group LASSO, as outras formas do block LASSO, o
  `gam.matched`, o `gam.k128` e o `wafc.gcv`.
- **O núcleo:** as cinco células do piloto (`smooth`, `uneven`, não
  homogêneo, `mixed`, nulo) × `n ∈ {250, 500, 1000, 2000}` × 100 réplicas.
  Custo medido no piloto, por réplica em `n = 1000`: ~450 s somando os
  métodos (`gam.gcv` ~230 s, `bsgl` ~120 s, WAFC ~40 s, `gam.reml` ~8 s;
  `gam.reml` ~86 s na `mixed`); `n = 2000` estimado em ~2,5 vezes. Com 100
  réplicas, ~290 h de processador, ~37 h em 8 núcleos; com 50, a metade.
- **Braços**, um fator por vez em `n = 1000` nas células não homogênea e
  `mixed`, 100 réplicas, ~3 h cada: ruído alto (razão sinal-ruído 1 contra o
  padrão 3 das duas células, D63); moduladoras correlacionadas (`u_rho = 0,5`, que exercita o termo
  cruzado); `X` dependente de `U` (D13; pede uma opção nova no `dgp.R`); e,
  opcional, escala `(p, q) = (6, 4)` com 6 blocos ativos e a grade `2:7`, com
  o custo medido na E4.3 (o `(10, 5)` original tem 12 750 colunas em `J = 8`
  e não cabe).
- **~~O topo das grades em `n = 2000`~~ substituído por D69 (2026-10-06):** a grade de `n = 2000` é fixada de antemão (`J` em `2:9`, `k = 120`), e o piloto mede só tempo e custo. Texto original (D60(e)): a E4.3 mede a fração das
  réplicas no topo com as grades de D34 e D41; acima de 20% (o gatilho de
  D34), a grade daquele método sobe só em `n = 2000` (`J` em `2:9` no WAFC,
  `k = 120` acrescentado no `gam`, com o custo do `gam.gcv` medido). A
  regra é por `n`, não por célula (D41). A escala entra se custar menos de
  ~10 h.
- **Fora:** erro dependente (D56; a aplicação cobre o caso) e a dimensão
  casada (E6.1a).
- **Métricas:** erro de predição e ISE em separado, o ISE dividido em blocos
  ativos e nulos; `P(Ŝ = S)` com as duas regras; o tempo com a busca; o `J`
  escolhido; no supp, a inclinação do erro contra `n` diante da taxa da
  teoria.

### E4.3 Piloto

**Fechada em 2026-10-07 (E4.3a, E4.3, E4.3b; D69, D70, D77 a D80):** a
escala entra, a grade de `n = 2000` foi fixada de antemão, a do `gam.reml`
cresce até 240 com a regra de tamanho e o limite de iterações; a E4.4
custa ~592 h de processador.

Réplicas, grade de `J` e de `λ`, e `nfolds` fixados por piloto **antes** da
produção, registrados em `config/`. A **tabela de avaliação da base** também:
uma por `(family, filter.size)`, construída com `WaveBased::wtable()` no
início da corrida e passada em `wavelet.table`, porque é o caminho rápido e a
regra `auto` não dispara nos `n` do estudo (D31).

### E4.4 Produção e agregação

**Próximo passo (2026-10-08):** na desktop do autor, 10 trabalhadores,
depois da junção entre máquinas (D81; `CONTINUAR.md` §2b).

Figuras e tabelas com os nomes que o `.tex` vai referenciar; cópia para
`results/` deste repositório; um parágrafo por figura.

**Critério de saída de E4:** estudo completo sob controles fixos; números em
`results/tables/`; parágrafo por figura; `PROVENANCE.md` no compêndio.

---

## E5. Manuscrito

### E5a. Seções que dependem só da teoria (abre quando E1.2 a E1.6 tiverem enunciado e L2 fechar)

`manuscript/ms_1.tex`, `supp_1.tex`, `references_1.bib`. Seções 1 a 4 de
`alvo-revista.md` §5: Introduction; Model and wavelet approximation space;
Theory (enunciados; provas no suplementar); Computation and tuning. Nome do
método fixado aqui. **Fechada em 2026-09-19.**

### E5c. Versão `k = 2` (fechada em 2026-09-30)

As edições acumuladas desde E5a, marcadas em `colR1` (D38 para a troca de
"sieve"), e o Corolário 8 na §3.6 com a prova no `supp` (D32). `ms_2` em
33 páginas e `supp_2` em 30; as pendências estão na pergunta 32 do
`ESTADO.md`.

### E5d. Versão `k = 4` (fechada em 2026-10-04)

O estimador e a teoria de D44 a D51 no manuscrito: o posicionamento de D50,
a penalidade em blocos na §2.3, a teoria em blocos na §3 e no supp (o
Corolário 11 como teorema principal), a taxa lenta em blocos no corpo, a
interface de D48 na §4, e a teoria do LASSO numa seção do supp (D51). As
remoções inteiras ficam em blocos cinza sem tachado (D51). Sem as Seções 5
a 7, que são da E5b.

### E5g. A cota inferior na `k = 4` (fechada em 2026-10-05)

A introdução e a terceira leitura do Theorem 1 passam à otimalidade em
`π ≥ 2` para todo `q`, e o Teorema 5 entra no supp no lugar da Remark S6.1
(D57).

### E5b. Seções que dependem dos números (abre quando E4 e E6 fecharem)

Simulation; Application; Discussion. Verificação de teto; cada número
conferido contra `results/tables/`.

**Critério de saída de E5:** compila limpo; o autor leu inteiro; números
conferem; checklist de `alvo-revista.md` §6.

---

## E6. Aplicação

### E6.1 Escolher os dados

Critérios: efeito de `X_ℓ` que plausivelmente varia com duas ou mais
moduladoras, `n` na casa dos milhares, dados públicos e citáveis, e
interpretação que um leitor da revista reconheça. Candidatos a levantar
(decisão do autor): salários (efeito de escolaridade variando com idade e
experiência; CPS ou PNAD); consumo de energia ou poluição (efeito de
temperatura variando com hora e umidade); dados de saúde com efeito de
tratamento modulado por idade e IMC. Entregável: `docs/aplicacao.md` com a
escolha e a razão.

**E6.1a (2026-09-20)** sondou três bases com veredito negativo no critério
de adaptação (`aplicacao-candidatas.md`). **E6.1b (catalogada em
2026-10-03)** refaz a sondagem no critério de D44, estrutura estável com
predição competitiva, nas três e em até três candidatas novas com salto
documentado; a escolha continua do autor (pergunta 2 do `ESTADO.md`).

**E6.1c (fechada em 2026-10-06):** marylebone com dados de licença
declarada, que fica como a aplicação principal (UK-AIR e vento de fonte declarada); **D57 escolheu duas
aplicações**, marylebone (estrutura) e beijing.heat (predição).

### E6.2 Ajustar e comparar

**Fechada em 2026-10-07 (E6.2, E6.2b; D73, D75, D76):** as duas bases
rodaram no compêndio (168 unidades, 0 falhas), as tabelas e figuras em
`results/`, os números na `k = 5`.

WAFC contra `mgcv` e linear; diagnóstico das funções estimadas; figura com
`β̂_ℓ(u)` interpretável.

**Critério de saída de E6:** reproduzido pelo compêndio de ponta a ponta;
figura e parágrafo de interpretação prontos para E5b.

---

## E7. Fechamento e submissão

1. Código de `wafc/` movido ao `WaveBased` (D61), versão etiquetada; DOI Zenodo do
   código e do compêndio.
2. Suplementar montado: provas, tabelas extras, diagnósticos.
3. Checklist de `alvo-revista.md` §6 inteiro; instruções da revista relidas
   na semana da submissão.
4. Carta de apresentação com a frase-tese e editores associados sugeridos.
5. Submissão.

---

## Ordem sugerida de execução

Quatro trilhas que podem correr em paralelo depois de E1.1:

- **Trilha teórica:** E1.2 → E1.3 → E1.4 → E1.5 → E1.6 (E1.3 e E1.4 são
  independentes entre si; E1.5 precisa das duas), e daí E5a.
- **Trilha computacional:** E2.1 assim que E1.1 e o enunciado de E1.2
  existirem; E2.2 e E2.3 em paralelo depois de E2.1; E2.4 e E2.5; E3; E4 e
  E6 em paralelo. Tudo em `wafc/`.
- **Trilha de literatura:** L1 e L2 a qualquer momento, antes de E5a.
- **Trilha editorial:** E0.3 (autor); E5a quando a teoria tiver enunciados;
  E5b no fim.

O ponto de sincronização obrigatório é E2.5: antes dele não se congela a
interface nem se escreve a seção de computação em versão final.

---

## Riscos e o que se faz com cada um

| Risco | Sinal | Resposta |
|---|---|---|
| A condição de desenho dos produtos não fecha no caso geral | E1.4 (ii) não sai; conferência numérica mostra `λ_min` restrito indo a zero com `X` dependente de `U` | enunciar sob `X ⊥ U` (fatoração) e a versão geral como hipótese de alto nível, com a conferência numérica como evidência |
| A adaptatividade não aparece nos números | E2.4: WAFC não vence `mgcv` em bumps/blocks | `boundary = "interval"`, pesos adaptativos, limiarização em blocos; se persistir, o artigo é sobre seleção de estrutura (grupos) e o alvo muda para JCGS |
| Sardy & Ma (2024) cobre mais do que parece | L2 | o peso vai para E1.4 (produtos) e para a aplicação |
| Teto de 40 páginas, referências incluídas | `ms` passa | provas e tabelas secundárias ao suplementar; cenários da simulação reduzidos no corpo |
| Custo do desenho em `p q 2^J` grande | `glmnet` lento em E4 com `p q = 50` | desenho esparso; `wtable()`; grade de `λ` curta; a grade de `J` vai a 8 (D34), e o custo medido está no `ESTADO.md` (2026-09-28) |
| Sem aplicação convincente | E6.1 sem candidato | decidir cedo (pergunta 2 do `ESTADO.md`); a EJS tolera aplicação mais leve |
