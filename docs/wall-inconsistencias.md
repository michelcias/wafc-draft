# Inconsistências e pontos a conferir no material do WALL

Documento de apoio, escrito em 2026-09-21 a pedido do autor. Reúne o que o
trabalho do WAFC encontrou, **de passagem**, no material dos repositórios
vizinhos `../../wall-manuscript` e `../../WaveBased`. Nada aqui foi
investigado por si: são achados colaterais de tarefas que tinham outro
objetivo, e cada linha diz onde foi encontrado e com que evidência, para que
a verificação lá seja barata.

**O que este documento não é.** Não é auditoria do WALL, não foi escrito com
o material deles na mão, e não distingue o que é erro do que é escolha
deliberada cuja justificativa o WAFC passou a duvidar. Os dois tipos estão
marcados.

---

## 1. Bibliografia

### 1.1 `johnstone2019gaussian` atribuído à Cambridge University Press

**Tipo:** provável erro. **Onde:** `references_theo_1.bib` do WALL teórico.
**Encontrado em:** L4 (2026-09-19), ao copiar a entrada para cá.

A entrada registra o livro *Gaussian Estimation: Sequence and Wavelet Models*,
de Iain Johnstone, como livro da CUP. **Não se confirma:** não há registro no
Crossref; a busca em livros da Cambridge University Press por "Gaussian
estimation" + Johnstone devolveu **zero** resultados em 2026-09-19; a página
do autor oferece "Book Draft, version of September 16, 2019" e a folha de
rosto do PDF diz "Draft version, September 16, 2019", sem editora.

No `referencias-verificadas.bib` daqui a entrada ficou como `unpublished`,
com `note` e `url`. **Se o WALL for submetido citando como livro da CUP, a
referência sai errada.**

**Desfecho (2026-09-22):** corrigido no WALL (commit `4529c1d` do `wall-manuscript`): a entrada virou `unpublished`,
com `note = {Book draft, version of September 16}` e `url`, e a chave foi
mantida.

### 1.2 `hardle1998wavelets` grafa "Tsybakov, Alexandre"

**Tipo:** divergência de grafia. **Onde:** mesma fonte. **Encontrado em:** L3
(2026-09-19).

O Crossref e o Springer grafam **"Alexander"**. A entrada daqui foi copiada
sem redigitar, como manda a regra de cópia, e por isso herdou a grafia; a
divergência continua nos dois repositórios.

**Desfecho (2026-09-22):** corrigido no WALL (commit `4529c1d` do `wall-manuscript`) para "Alexander". A
cópia daqui, em `referencias-verificadas.bib`, ainda traz "Alexandre".

### 1.3 Sardy & Ma citado só como preprint

**Tipo:** desatualização. **Encontrado em:** L1 (2026-09-18).

O WALL cita `sardy2022sparse`, o preprint. O trabalho saiu em *Scandinavian
Journal of Statistics* **51**(1), 89-108, DOI `10.1111/sjos.12680`. Para o
WAFC isso importa porque é o vizinho mais próximo no método; para o WALL, é
só atualizar a entrada.

**Desfecho (2026-09-22):** corrigido no WALL (commit `4529c1d` do `wall-manuscript`). A entrada aponta para
a versão da *SJS*, com DOI, e a chave `sardy2022sparse` foi mantida para não
tocar nos `\cite`; o texto passa a mostrar "Sardy and Ma (2024)".

---

## 2. Código (`WaveBased`)

### 2.1 `wall()` com `boundary = "interval"`: a constante entra `d + 1` vezes

**Tipo:** defeito de identificabilidade, sem efeito no ajuste.
**Encontrado em:** E1.2 (2026-09-18), ao derivar a identificabilidade do
desenho do WAFC e comparar com o que o `wall()` faz.

Com `boundary = "interval"` o `wall()` mantém as `2^{j_0}` funções de escala
de **cada** covariável não penalizadas (`drop.phi = FALSE`) **e** o intercepto
do `glmnet`. Como a constante pertence ao espaço gerado pelas funções de
escala de cada bloco, ela entra `d + 1` vezes, e a matriz de Gram fica
**singular, com nulidade `d`**.

**Consequência:** a função ajustada não sofre — o `glmnet` resolve um problema
convexo com solução não única em direções que não afetam o ajuste —, mas **a
leitura dos coeficientes de escala não é interpretável**, porque o nível está
repartido arbitrariamente entre `d + 1` lugares. Quem ler `coef()` por bloco
naquele caso está lendo uma decomposição que o dado não identifica.

A correção, do lado de cá, foi a reparametrização `[Φ Q | Ψ]` descrita na §5
de `derivations/01-identificabilidade.md`, que impõe `μ' α = 0` por bloco e
devolve `2^J − 1` colunas, o mesmo número do caso periódico.

### 2.2 O padrão `eps = 1.9^{−J}` acopla a margem à resolução

**Tipo:** escolha cuja justificativa não se sustenta assintoticamente.
**Encontrado em:** E1.3b e E2.1b (2026-09-19), E2.4 (2026-09-20).

O `wall()` reescalona as covariáveis para `[eps, 1 − eps]` com
`eps = 1.9^{−J}`, para afastar os dados da faixa onde a base periodizada
carrega artefato. O WAFC herdou o padrão e depois o descartou, por três
medições:

- **Margem acoplada a `J` não compra taxa.** Com `eps ∝ 2^{−J}` o ganho é
  **zero** (queda medida de `0,52` a `0,58` bit por nível, contra `0,50` sem
  margem nenhuma); com `1,9^{−J}` é parcial (`0,96` bit). Só margem **fixa**
  compra a taxa cheia. A razão é que as duas exigências — excluir a faixa
  contaminada e não esconder wavelet inteira na margem — medem **a mesma
  largura**, `(L−1)2^{−J}`, e são disjuntas (Proposição 6 de
  `derivations/07-rota-intervalo.tex`).
- **A margem varia com `J`**, então cada candidato da validação cruzada
  estima um alvo com suporte reescalado diferente. O efeito no erro é
  pequeno, mas é confusão de **estimando**, não de estimador.
- `eps = 1,9^{−J}` **cai fora de `[0, 0.5)` em `J = 1`**, o que torna aquele
  candidato inconstruível.

Nada disso é erro do `wall()` no que ele se propõe a fazer; é aviso de que a
justificativa numérica da margem não sobrevive a `J → ∞`, e de que o valor
herdado é dominado (E2.4 mediu `eps = 0` melhor no cenário suave por 12% a
15% e pior por no máximo 2,4% no não homogêneo).

### 2.3 A grade de `J` de `cv.wall()` trunca o espaço de aproximação no ótimo

**Tipo:** possível defeito de sintonia, herdado e medido aqui.
**Encontrado em:** E2.4 (2026-09-20), e confirmado em dado real por E6.1a.

A grade `2:⌈log_2 n / 2⌉` foi herdada do `cv.wall()`. No piloto do WAFC, **o
oráculo da grade escolhe o `J` do topo em 100% das réplicas** de toda célula
com componente. Alargando para `2:8`, no cenário não homogêneo com
`n = 1000`, o erro de predição cai **17,6%** e o ISE **32%**, em **50 de 50
réplicas**; no cenário suave nada muda.

Em dado real (E6.1a) a grade foi *binding* em duas de três bases, com a curva
de validação ainda caindo no topo.

**Por que isso interessa ao WALL:** se os estudos de lá usam a mesma regra, a
resolução pode estar sendo truncada onde a função exige mais — e o efeito
aparece justamente nas funções não homogêneas, que são o argumento do método.
**Ressalva importante**, também medida aqui: em série com dependência, sob
partição aleatória a validação cruzada **pede mais resolução do que deveria**;
retendo semanas inteiras, o `J` pedido caiu de 8 para 4. Alargar a grade sem
arrumar a partição piora o diagnóstico.

---

## 3. Teoria (WALL teórico)

### 3.1 A trilha rápida precisa de localização, e isso custa hipótese

**Tipo:** diferença de arquitetura, não erro. **Encontrado em:** E1.6
(2026-09-19).

A trilha rápida do WALL, por localizar o estimador não modificado, precisa de
`(M*)² J 2^J = o(n)` — o que restringe `τ < 2/3` — e de um piso em `s`. Com
perda quadrática, o WAFC não precisa de nada disso: o Lema 8 de
`derivations/05-taxas.tex` é a prova da desigualdade oráculo relida com outro
comparador, e o intervalo de `τ` é `(0,2)` inteiro, com piso `s' > (2−τ)/4`,
que é o que o WALL só obtém para o **estimador restrito**.

Registrado aqui porque é material de comparação entre os dois artigos, e
porque, se o WALL quiser a mesma folga, o caminho está escrito.

### 3.2 O vocabulário: "sieve" é raro no público-alvo

**Tipo:** escolha de redação, não erro. **Encontrado em:** 2026-09-21, ao
preparar a terminologia do manuscrito do WAFC (D33).

O WALL teórico usa **"sieve" 38 vezes**, e é de lá que o termo entrou no
vocabulário do WAFC. Ele é padrão em estatística teórica e em econometria
semiparamétrica (Grenander 1981; Geman & Hwang 1982; Shen & Wong 1994;
Chen 2007), mas **incomum na literatura em que os dois artigos querem se
inserir**: conferido, **Klopp & Pensky (2015) e Xue & Yang (2006) não o usam
nenhuma vez**, e "Approximation space" é **palavra-chave** do artigo de Xue
& Yang.

O WAFC trocou (D33): nas derivações, por "espaço de aproximação", com uma
única menção retida entre parênteses no ponto em que o espaço é definido
(E1.10, 2026-09-21); no manuscrito, por "approximation space", na rodada de
`k = 2`. O dicionário usado está no handoff de E1.10 e é reaproveitável:
"sieve linear" virou **"aproximação linear"**, para fazer par com
"aproximação não linear", que já era o termo do texto; "viés do sieve" virou
**"viés de aproximação"**.

**Para o WALL isto é decisão de lá**, e há um argumento honesto dos dois
lados: o público da revista que o WALL mira pode ser outro, e o termo tem a
vantagem de nomear a construção de uma vez. O que este documento registra é
que, no público do WAFC, ele custa uma explicação que o texto não tem espaço
para dar.

### 3.3 Referências de adaptação ausentes do `.bib` do WALL

**Tipo:** lacuna. **Encontrado em:** L1 (2026-09-18).

Donoho & Johnstone (1994, *Biometrika*), Donoho & Johnstone (1998, *Ann.
Statist.*), Daubechies & Lagarias (1991, 1992) e Bickel, Ritov & Tsybakov
(2009) **não estão** em `references_theo_1.bib`; o WALL cita só o Donoho &
Johnstone de *PTRF* 1994. Todas foram verificadas aqui e estão em
`docs/referencias-verificadas.bib`, prontas para copiar sem redigitar.

**Desfecho (2026-09-22):** as cinco foram copiadas para `references_theo_1.bib`
(commit `4529c1d` do `wall-manuscript`), com as chaves `donoho1994ideal`, `donoho1998minimax`,
`daubechies1991two`, `daubechies1992two` e `bickel2009simultaneous`. O
manuscrito ainda não cita nenhuma delas, então nada muda no PDF até alguém
incluir um `\cite`.

---

## 4. Como usar este documento

Cada item traz onde foi encontrado e com que evidência; nenhum exige refazer
medição para ser conferido no WALL. A ordem de custo crescente para tratar
lá: as três da §1 são edição de `.bib`; a §2.1 é um argumento `drop.phi` e a
leitura de `coef()`; a §2.3 é uma linha na grade, mas pede repetir os estudos
que dependem dela; a §2.2 e a §3.2 são decisões de redação e de desenho; a
§3.1 é material de discussão.

Se algum item for tratado no WALL, vale registrar aqui o desfecho, para que
os dois projetos não redescubram o mesmo.

---

## 5. Conferência no WALL (2026-09-22)

Leitura deste documento contra o material de lá, isto é, o `wall-manuscript`,
os scripts do benchmark em `~/Documentos/wall` e o `WaveBased`. Os itens de
§1 e §3.3 foram tratados (ver os desfechos acima). O resto ficou assim:

- **Caminhos.** Os repositórios estão em `~/Documentos/wall-manuscript` e
  `~/Documentos/WaveBased`, ou seja, `../../../` a partir desta pasta, e não
  `../../`.
- **§2.1 tem alcance maior.** Em `wall.R`, `drop.phi` só é `TRUE` com
  `boundary = "periodic"` **e** `j0 == 0`. O caso periódico com `j0 > 0`
  mantém as funções de escala e tem a mesma singularidade.
- **§2.2, `J = 1`: o candidato não é inconstruível.** `.wall_eps()` devolve
  o padrão `1,9^{-J}` antes de validar o intervalo `[0, 0.5)`. Com
  `eps ≈ 0,526`, o mapa leva `[0,1]` para `[0,474; 0,526]`, invertido, sem
  nenhum aviso. O benchmark do artigo aplicado inclui `J = 1` na grade
  (protocolo de `ms_app_5.tex`), e `wall_J_profile.csv` dá
  `inner_score = 0,5` em `J = 1` nas 138 dobras. Esse nível nunca foi
  selecionado, então os resultados não mudam, mas a afirmação
  `x* ∈ [ε, 1−ε]` da equação `eq:rescale` é falsa em `J = 1`.
- **§2.3 não vale para o artigo aplicado.** O benchmark não usa a grade
  padrão do `cv.wall()`. Usa `J ∈ {1, …, J_max}` com
  `J_max = min{9, ⌊log₂(4n/d)⌋}`, e o `J` selecionado nunca passou de 7. O
  alerta continua valendo para quem usa o `cv.wall()` com a grade padrão.
- **§3.2.** A contagem confere: são 38 linhas, com 39 ocorrências.
