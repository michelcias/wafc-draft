# Estado do trabalho, handoff de continuidade

**Última atualização:** 2026-09-19.
**Etapa corrente:** **E1 fechada** (E1.7 é condicional a E2.5) e L1, L2 e
L3 fechadas; de E2 falta E2.3 em diante. E0 deve só a ratificação de D5.
Liberada para chat de tarefa: E2.3. Nenhum chat de tarefa em curso. Duas
decisões do autor adiadas (D5, D8), seis propostas de L2 e duas decisões
novas (D16, D17) a ratificar, e doze perguntas na §4. **E5a está a um passo:
depende só de D5 e D8.**
**Versão viva do manuscrito:** nenhuma (nasce em E5a como `k = 1`).
**Cor da rodada corrente:** `colR1` (entra em uso quando existir `k = 2`).

Este documento é o ponto de partida de cada sessão. Ele diz onde o trabalho
parou, o que já foi decidido (para não reabrir) e o que vem a seguir.

---

## 1. Para começar uma sessão

1. Ler `CLAUDE.md` (as duas regras de abertura: sem coautoria, respostas
   curtas).
2. Ler este arquivo.
3. Ler [`instrucoes.md`](instrucoes.md). É a regra, não a sugestão.
4. Ler a etapa em foco em [`plano-projeto.md`](plano-projeto.md). Em E1,
   também [`proposta-metodo.md`](proposta-metodo.md).
5. Se for escrever fórmula: [`notacao.md`](notacao.md) antes, e as Seções 2,
   3 e 5 de `../../wall-manuscript/manuscript/theo/ms_theo_1.tex`, que é o
   molde da prova.

### Mensagem de abertura sugerida

Chat principal (o que edita este arquivo e commita):

> Leia `CLAUDE.md`, `docs/ESTADO.md` e `docs/instrucoes.md` para se situar.
> Hoje vamos <foco do dia>.

Outra máquina (instalação e disposição das pastas em `CONTINUAR.md`):

> Leia `docs/CONTINUAR.md` e me diga o que falta instalar; depois leia
> `docs/ESTADO.md` e continuamos de onde parou.

Chat de tarefa (uma etapa, arquivos restritos, termina em handoff):

> Leia `docs/TAREFA.md` e execute a tarefa <etapa>.

### Git

O autor commita e faz push. O assistente lê o estado e lembra de commitar nos
momentos apropriados. **Sem coautoria nas mensagens**, sem exceção.

---

## 2. Onde o trabalho está

### 2026-09-18: o repositório nasceu

Sessão de avaliação de viabilidade e de planejamento. O autor descreveu a
ideia (coeficientes funcionais, cada um aditivo nas moduladoras, cada
componente em wavelets, estimação por LASSO); o assistente leu o `bdm-draft`
(molde das convenções), o `WaveBased` (`wall()`, `wbasis()`, `wtable()`), o
`wall` (compêndio) e o `wall-manuscript` (o teórico é o molde da prova), e
fez buscas de literatura. O que saiu:

- **Viável e com nicho.** A combinação "coeficientes aditivos + wavelets +
  LASSO" não apareceu na busca inicial; os três pilares existem em separado
  (Xue & Yang 2006; Zhou & You 2004; Sardy & Tseng 2004 e Sardy & Ma 2024;
  Wei, Huang & Li 2011). O trabalho mais próximo é Sardy & Ma (2024), sem o
  `X_j` multiplicando. Detalhe em [`proposta-metodo.md`](proposta-metodo.md)
  §5 e [`literatura.md`](literatura.md); a confirmação é L2.
- **O motor já existe.** A matriz de desenho do WAFC é a do `wall()` com cada
  bloco multiplicado por `X_j`; o ajuste é o mesmo `glmnet`. Inventário em
  [`inventario-codigo.md`](inventario-codigo.md). A arquitetura de prova do
  WALL teórico (sieve, Besov, oráculo, compatibilidade, compressibilidade)
  transfere; o que é novo na teoria é a condição de desenho para produtos
  (E1.4).
- **O código fica aqui.** O autor decidiu (D4) que o `WaveBased` não recebe
  código por enquanto: o método vive na pasta `wafc/` (`R/`, `tests/`,
  `scripts/`), organizada como pacote sem ser pacote, e o `WaveBased`
  instalado é só dependência para as bases. Empacotar decide-se em E3.3.
- **Revista.** Proposta: *Statistica Sinica* como alvo primário (linhagem do
  modelo, formato completo, 30 páginas), EJS como reserva (sem teto, open
  access, template já versionado e compilando em `manuscript/ejs-template/`).
  Comparação e requisitos em [`alvo-revista.md`](alvo-revista.md). O SCImago
  e o site da SS estavam inacessíveis desta máquina; a confirmação é do
  autor (E0.3).
- **Plano** em etapas E0 a E7 mais L1 e L2, com dependências e trilhas
  paralelas ([`plano-projeto.md`](plano-projeto.md)); catálogo de tarefas
  abertas em [`TAREFA.md`](TAREFA.md).
- **Ferramentas locais:** R 4.6.1, `glmnet`, `grpreg`, `gglasso`, `mgcv`,
  `latexmk`; o `WaveBased` 2.6-0 foi instalado de `../../WaveBased` e o
  `wall()` roda. Faltam `sparsegl` (E2.2) e `gh`.

### 2026-09-18: E0.3 fechada (instruções e template da Statistica Sinica)

O autor capturou as páginas oficiais da revista e baixou os templates; as
instruções estão transcritas em inglês em
[`ss-instrucoes-autores.md`](ss-instrucoes-autores.md) e o
[`alvo-revista.md`](alvo-revista.md) foi corrigido. Os números que fecham a
etapa:

- **Teto de 40 páginas em espaço duplo no template, referências e apêndice
  incluídos** (12pt, margens de 1 polegada), e não as 30 que a busca na web
  tinha dado. A estrutura-alvo soma ~29 páginas, com ~11 de folga não
  alocada.
- **Template obrigatório e conferido na recepção**: manuscrito fora do
  template volta antes de qualquer revisão. Versionado em
  `manuscript/ss-template/`: `SS-template.tex` (bibliografia manual),
  `SS-template-bib.tex` (BibTeX, `chicago`) e `supp-temp_20240820.tex`;
  os três compilam com o TeX Live 2023 local (avisos de `fancyheadings`
  obsoleto e de `\headheight`, sem erro).
- **Triagem:** ~40% das submissões passam dos co-editores; o AE manda a
  referee ~60% do que recebe. Alvo de prazo: 2 a 3 meses na primeira
  rodada.
- **Suplementar:** um único PDF de até 10 MB, mesmo título e mesmos autores,
  revisado junto; seção "Supplementary Material" como última seção do corpo,
  antes dos agradecimentos. A revista *prefere* provas, lemas técnicos e
  detalhes de simulação no suplementar, o que é o que o plano já fazia.
- **Reprodutibilidade é exigência editorial**, não cortesia: "software
  producing the evidence should be available for examination as well as
  pertinent datasets".
- **O que a página não diz**, listado na §7 de `ss-instrucoes-autores.md`: o
  **tipo de revisão** (nem cega simples nem dupla; sem exigência de
  anonimização), taxas, licença, política formal de dados e código, limite
  de resumo e de palavras-chave, e qual dos dois templates (Windows ou Mac)
  é o zip baixado. Um símbolo da seção de fórmulas veio como imagem quebrada
  na captura.

### 2026-09-18: E1.1 fechada (notação congelada)

O autor ratificou os pontos da §6 do [`notacao.md`](notacao.md), que sai do
estado de esboço; [`../derivations/macros.tex`](../derivations/macros.tex)
foi reescrito para implementá-lo e compila (documento de teste com o modelo,
a expansão do bloco e os ambientes). O que ficou:

- **A base fica com `j` e `k`.** `ψ_{jk}` mantém o padrão da literatura de
  wavelets; quem muda de letra são as covariáveis: linear `X_ℓ`,
  `ℓ = 1, …, p`; moduladora `U_m`, `m = 1, …, q`. Isso inverte a proposta do
  esboço, que punha a covariável em `j` e a wavelet em `ψ_{lm}`.
- **Correção de um erro do esboço:** o `notacao.md` dizia que o WALL usa `j`
  para o nível. Não usa: o WALL teórico escreve `ψ_{ℓ,k}` e depois achata o
  par num índice único `b_m` (`ms_theo_1.tex`, §2.2). A notação do WAFC
  mantém o par visível, porque o bloco `(ℓ, m)` é a unidade do desenho e o
  grupo da variante de E2.2.
- **Coeficientes em `θ`**, com `β_ℓ` reservado ao coeficiente funcional; o
  subscrito é `θ_{ℓm,jk}`, bloco antes da wavelet.
- **`U ∈ [0,1]^q` por hipótese populacional**; a transformação monótona da
  prática fica na seção de computação e fora da teoria de E1.3 a E1.6.
- **Colisão registrada:** `ℓ` é índice e é a letra da penalidade `ℓ_1`. A
  convenção que ficou é escrever a penalidade `‖θ‖_1` e nunca colar `ℓ_τ` a
  um índice (`notacao.md`, §1).
- **Ordem das colunas de `Z`** fixada: os `p` termos não penalizados, depois
  os blocos `(ℓ, m)` em ordem lexicográfica, e dentro do bloco `j` crescente
  e `k` crescente. É o que E2.1 tem de implementar.

### 2026-09-18: L1 fechada (verificação bibliográfica)

Chat de tarefa, integrado neste. `docs/referencias-verificadas.bib` tem
**35 entradas**, cada uma conferida no Crossref, e `literatura.md` não tem
mais `[VERIFICAR]` nem `[L1: confirmar]`. Conferido aqui: 35 `@`-entradas,
35 `\bibitem` com `bibtex` e `plain`, sem erro nem aviso; 32 com DOI, e as
três sem DOI são Xue & Yang (2006), Xue & Qu (2012, JMLR) e Haris, Simon &
Shojaie (2018, NeurIPS), que não têm registro.

O que a verificação corrigiu, e que muda o que o artigo cita:

- O arXiv 1903.04631 que `literatura.md` atribuía a Amato, Antoniadis, De
  Feis & Gijbels é de **Haris, Simon & Shojaie** (NeurIPS 2018): aditivos
  com wavelets em desenho irregular. O candidato ao que a linha queria
  citar é Amato et al. (2022, *Stat. Comput.*). As duas entradas ficaram no
  `.bib`, e **as duas entram na varredura de L2**; pelo resumo, nenhuma tem
  o `X_ℓ` multiplicando.
- Antoniadis & Gijbels (2014): o terceiro autor é **Lambert-Lacroix**, não
  Verhasselt.
- Klopp & Pensky é de **2015** (*Ann. Statist.* 43(3) 1273-1299).
- Sardy & Ma (2024) saiu em *Scand. J. Statist.* 51(1) 89-108
  (10.1111/sjos.12680); o WALL cita só o preprint. É o trabalho mais
  próximo, então a citação tem de ser a publicada.
- Donoho & Johnstone (1994 *Biometrika*; 1998 *Ann. Statist.*), Daubechies
  & Lagarias e Bickel, Ritov & Tsybakov estavam marcados `no wall` mas não
  estão no `references_theo_1.bib`; entraram com chave própria.

**Lição de ferramenta** (vale para L2 e para E5b): Crossref e OpenAlex
respondem desta máquina sem chave; o Semantic Scholar limita a ~1
consulta/s; o Project Euclid bloqueia `curl`; a *Statistica Sinica* serve os
PDFs antigos em `statistica/oldpdf/A{vol}n{iss}{art}.pdf`.

### 2026-09-18: L2 fechada (busca de novidade), e o vizinho é outro

Chat de tarefa, integrado neste. `docs/busca-novidade.md` traz as quatro
buscas, a varredura da *Statistica Sinica* e da EJS desde 2015 e o veredito
por contribuição; as linhas novas de `literatura.md` foram coladas aqui,
porque L1 tinha o arquivo na mão quando L2 rodou.

- **A novidade central se sustenta:** nenhum trabalho combina coeficientes
  aditivos em várias moduladoras, wavelets e LASSO. Na varredura por título,
  a interseção "coeficientes variáveis × wavelet" é **zero** na SS e **zero**
  na EJS desde 2015 (SS: 28 `varying coefficient`, 1 `additive coefficient`,
  2 `wavelet`; EJS: 7, 0, 6).
- **Mas o vizinho perigoso não é Sardy & Ma: é Klopp & Pensky (2015,
  *Ann. Statist.* 43(3) 1273-1299).** Para `q = 1` e `X ⊥ U`, eles já têm o
  desenho de produtos com Gram `Ω ⊗ Φ` (eq. 1.8), a concentração da Gram
  empírica restrita (Lema 1), o LASSO em blocos com desigualdade oráculo não
  assintótica, a taxa adaptativa em Besov com `ν < 2` e a cota inferior
  minimax. Ou seja: E1.4 (i), E1.5 e E1.6 estão publicados no caso de uma
  moduladora com desenho independente.
- **O que sobra para o WAFC**, e que passa a ser o eixo do artigo:
  aditividade em `q ≥ 2` moduladoras, com o termo cruzado
  `E[X_ℓ X_{ℓ'} ψ_{jk}(U_m) ψ_{j'k'}(U_{m'})]`, `m ≠ m'`, que não é produto
  de Kronecker; `X` dependente de `U`; LASSO puro com o corolário de
  compressibilidade; identificabilidade no nível da base com constantes não
  penalizadas; e o software, a simulação e a aplicação, que K&P não têm.
- **Sardy & Ma (2024) não é ameaça teórica:** os quatro teoremas deles são
  de otimização (degenerescência do group square-root LASSO em bloco
  ortonormal, limiarização suave em forma fechada, relaxação por blocos,
  SURE), sem teoria estatística e sem sieve.
- **Precursor do próprio grupo, obrigatório citar:** Montoril, Morettin &
  Chiann (2018, *IJWMIP* 16(1) 1850004), wavelets em coeficientes funcionais
  sem penalização.
- **Concorrentes que o referee vai cobrar em E4:** o spline adaptativo de
  Wang, Jiang & Liu (2024, *JCGS*), o block LASSO de K&P no mesmo desenho
  (sai com `grpreg`/`gglasso`) e um não aditivo em várias moduladoras
  (VCBART).
- **Lição de ferramenta:** o Crossref por ISSN com `from-pub-date` varre uma
  revista por título; o OpenAlex dá grafo de citações e resumos; a Wiley
  bloqueia leitura automatizada.

As linhas que L2 acrescentou ao `literatura.md` estão com status `resumo` ou
`[VERIFICAR]` e **não podem entrar no `.bib`** antes de uma verificação como
a de L1.

### 2026-09-18: E1.2, E1.3 e E1.4 fechadas (o núcleo da teoria)

Três chats de tarefa, integrados neste. Rodei as três conferências aqui
antes de integrar, e as três imprimem `OK` e reproduzem os números dos
handoffs: `01-identificabilidade.R` (7 s), `03-desenho-produtos.R` (9 s),
`02-aproximacao-besov.R` (53 s).

**E1.2, identificabilidade** (`derivations/01-identificabilidade.md`).
Proposição 1 (das `β_ℓ` a partir de `f` sob `E[XX' | U]` não singular q.c.;
de `(c_ℓ, g_{ℓm})` a partir de `β_ℓ` sob densidade positiva e centralização)
e Lema 1 (a base periódica com `j_0 = 0` impõe a centralização no nível da
base; `(c, θ)` é injetiva e `Σ` definida positiva a cada `J`). Números:
`Z` com 30 colunas tem posto 30 em `n = 200`, `λ_min(Z'Z/n) = 7.8e-03`,
recuperação sem ruído com erro `3e-15`. As deficiências previstas
aparecem exatas: manter `φ_{00}` dá nulidade `4 = p q`; `U_2 = U_1` dá
`14 = p N_J`; `X_2 = 1 + ψ_{10}(U_1)` dá `1`.

**E1.3, aproximação em Besov** (`derivations/02-aproximacao-besov.tex`,
7 páginas). Lema de aproximação em `L_2` e em `L_∞`, corolário para
`f - f_J`, e uma proposição sobre o custo da periodização. Números: em
`sin(2πu)` a queda por nível é `4.00` bits (isto é `2^{−JN}`, `N = 4`); em
`bumps` a queda média em `j = 9..12` é `1.37`, contra `s' = 3/2`; numa
função que não emenda em `0 ≡ 1`, a queda trava em `0.50` bits e o erro
uniforme **não decai** (`0.54` e `0.92` em `J = 9`).

**E1.4, Gram do desenho de produtos**
(`derivations/03-desenho-produtos.tex`, 6 páginas), que entregou mais do que
o catálogo pedia: no caso geral sai o **autovalor mínimo cheio**, não só a
compatibilidade no cone. `λ_min(Σ) ≥ κ_1 c_U` e `λ_max(Σ) ≤ κ_2 C_U`.
Números: sob `X ⊥ U` a fatoração `Π(E[XX'] ⊗ Σ_Ψ)Π'` bate a `7.6e-13`; com
`X` dependente de `U`, `λ_min(Σ) = 0.2747` contra a cota `κ_1 c_U = 0.2503`
(10% de folga); com `κ_1 = 0`, a Gram degenera com `J`, como previsto
(`4e-04`, `<1e-04`, `<1e-04` em `J = 2, 3, 4`).

O que isso muda no rumo:

- **A pergunta 4 da §4 está respondida: o caso geral fecha** (vira D13). A
  independência compra só a fatoração exata, não constante melhor.
- **A linha de risco "a condição de desenho dos produtos não fecha em
  geral" sai** da tabela de `proposta-metodo.md` §6, com número.
- **A taxa `2^{−Js'}` só aparece quando `2^{−J}` fica abaixo da escala da
  função** (em `bumps`, `J ≥ 9`). Os `n` de E4 (`≤ 2000`) com
  `J_n ≍ log_2 n/(2s'+1)` não chegam lá: a simulação vive no regime
  pré-assintótico, que é justamente onde o sieve linear perde e o LASSO
  ganha. Isso é matéria da Seção 5 do artigo, não defeito da conferência.
- **Periodização custa `1/2`:** para qualquer `g` que não emende, a taxa em
  `L_2` trava em `2^{−J/2}` e o erro uniforme não converge. Daí a pergunta 6
  da §4.
- **Em `n` pequeno a compatibilidade é mais saudável que o autovalor cheio:**
  `φ²(S)` fica em 65-75% do populacional em `n = 200`, contra 40-50% de
  `λ_min(Σ̂)`. É `φ²(S)` que E1.5 consome.
- **Achado sobre o `wall()`** (não tocado, D4): com `boundary = "interval"`
  ele mantém as `2^{j_0}` funções de escala de cada covariável e ainda o
  intercepto do `glmnet`, logo a constante entra `d + 1` vezes e a Gram fica
  singular com nulidade `d`. O ajuste não sofre; a leitura dos coeficientes
  de escala, sim. Fica registrado para o autor decidir se avisa lá.

**Numeração global dos resultados** (TAREFA.md §5 tem a tabela): E1.2 fica
com Proposição 1 e Lema 1; E1.3 com Lema 2, Lema 3, Corolário 1 e
Proposição 2; E1.4 com Proposição 3. Os arquivos mantêm os contadores
locais; o mapa é a autoridade quando E5a montar o manuscrito.

### 2026-09-19: E1.5 e E2.1 fechadas (o oráculo e o primeiro código)

Dois chats de tarefa, integrados aqui. Conferências rodadas nesta máquina
antes de integrar: `04-oraculo.R` imprime `OK` em 16 s;
`testthat::test_dir("wafc/tests")` dá **103 passam, 0 falham** em 4,6 s;
`wafc/scripts/01-smoke.R` imprime `OK` em 2,8 s; o `.pdf` de E1.5 tem
9 páginas.

**E1.5, desigualdade oráculo** (`derivations/04-oraculo.tex`). Lemas 4 a 7,
Teorema 1 e Corolários 2 e 3. A estrutura da prova evita a condição de cone
por **perfilagem**: com `B̃ = (I − Π_A)B`, o estimador resolve o problema
residualizado e `λ_min(B̃'B̃/n) ≥ λ_min(Σ̂)`, que é a forma consumível de
E1.4. Números: com `s_0 = 4` e sem viés, a razão `‖f̂ − f‖²_n/(λ² s_0)` fica
entre `1.541` e `1.617` ao variar `n` de `200` a `6400` (máximo sobre
mínimo `1.05`), que é exatamente a escala que o plano mandava medir; nenhuma
violação da cota rápida em 360 ajustes, com folga de `371×` a `187×`.

- **D13 pagou aqui:** é o autovalor cheio de E1.4 que dispensa o cone e dá
  constante que **não depende de `S_0`**.
- **Um termo que Klopp & Pensky não têm: `σ² p / n`**, preço de deixar os
  `c_ℓ` fora da penalidade (D3), já que K&P penalizam a constante. A
  diferença para eles não é só de generalidade, aparece no enunciado.
- **Calibração de `λ`:** a leitura certa de `‖Z‖_max` no plano é
  `σ̂_max = max_a sqrt(Σ̂_aa)`, que é `O_p(1)`; a leitura literal
  `max_{i,a}|Z_{ia}|` é de ordem `2^{J/2}` e **destruiria a taxa de E1.6**.
  Medido: `σ̂_max` vai de `1.02` a `1.14` em `J = 2..6` enquanto
  `max|Z_{ia}|` vai de `3.04` a `11.47`.
- **O Lema 7 dá `‖θ*‖_1` sem dependência de `π`**, o que permite enunciar a
  taxa lenta em `B^s_{π,r}` geral, e não só em `B^s_{2,∞}` como o WALL.
- **Sob Besov, em geral `s_0 = d`:** o enunciado incondicional é a taxa
  lenta, e trocar `s_0` por dimensão efetiva é o corolário de
  compressibilidade de E1.6.

**E2.1, desenho e cenários** (`wafc/R/load.R`, `dgp.R`, `design.R`,
`tests/test-design.R`, `scripts/01-smoke.R`). `wafc_design()` implementa a
ordem de D12 e descarta só o que E1.2 manda (a coluna `φ_{00}` de cada
bloco); a ordem bate com a construção de referência de
`check/03-desenho-produtos.R` a `1e-12`; mínimos quadrados sem ruído com
erro `1.1e-15` e `glmnet` com `λ → 0` a `4.3e-06`. No piloto de fumaça, o
cenário nulo zera as 90 colunas penalizadas (RMSE `0.045` contra
`σ = 0.62`): o desenho não inventa estrutura.

- **Achado que fixa E2.2:** o `glmnet` **descarta colunas de variância
  zero** mesmo com `standardize = FALSE`, logo `X_1 ≡ 1` sai do ajuste. Com
  `intercept = FALSE` perde-se `c_1` inteiro; com `intercept = TRUE` o
  intercepto carrega `c_1` exatamente (`1.00000025` contra `1`). `wafc()`
  ajusta com `intercept = TRUE` e soma o intercepto ao nível da covariável
  constante; `wafc_design()` já devolve `constant` com esses índices. E1.5
  chegou ao mesmo achado por outro caminho e acrescentou uma terceira rota
  conferida (resolver o problema residualizado), além de duas notas: o
  `glmnet` **reescala `penalty.factor` para somar `nvars`**, e `thresh` está
  obsoleto na versão 5.0 (vai em `control`).
- **`lambda.min` deixa lixo nos blocos nulos** (3 e 8 coeficientes nos dois
  blocos de `β_3`, ISE `0.0010` e `0.0080`), e `lambda.1se` reduz sem zerar.
  É o argumento numérico a favor da variante em grupos de E2.2.
- **O reescalonamento herdado do `wall()`** (`eps = 1.9^{−J}`, isto é
  `[0.077, 0.923]` em `J = 4`) faz a componente estimada integrar zero na
  amplitude reescalada, não em `[0,1]`. É deslocamento de nível, não de
  forma, e é matéria da seção de computação.
- **Lição de ferramenta:** `Matrix::cbind2` só aceita dois argumentos, e
  `do.call(Matrix::cbind2, lista)` devolve em silêncio o `cbind` dos dois
  primeiros; usar `Reduce`.

**Numeração global** (mapa na §5 do `TAREFA.md`): E1.5 ficou com Lemas 4 a
7, Teorema 1 e Corolários 2 e 3, e **passou a imprimir o número global** via
`\setcounter`. Adotei a prática; `02` e `03` continuam com contador local e
o mapa é a autoridade até alguém alinhá-los.

### 2026-09-19: E1.6, E2.2 e L3 fechadas; E1 inteira fechada

Três chats de tarefa, integrados aqui. Conferido nesta máquina antes de
integrar: `05-taxas.R` imprime `OK` (24 verificações); a bateria de `wafc/`
dá **231 passam, 0 falham** em 14 s; `01-smoke.R` imprime `OK`; o `.bib` com
57 entradas compila com `bibtex` sem erro nem aviso.

**E1.6, taxas e dimensão efetiva** (`derivations/05-taxas.tex`, 10 páginas).
Lemas 8 e 9, Proposição 4, Teorema 2, Corolários 4 e 5. O eixo é um achado
que muda o enunciado principal do artigo: **a hipótese de Besov de E1.3 já
contém a compressibilidade** (Lema 9(i): weak-`ℓ_τ` com
`τ = (s + 1/2)^{−1}`). O corolário de compressibilidade não é hipótese
extra, é leitura mais fina da mesma hipótese, e é ele que recupera a taxa de
aproximação **não linear** `n^{−2s/(2s+1)}` onde o sieve linear só chega a
`n^{−2s'/(2s'+1)}`; as duas coincidem em `π ≥ 2` e se separam em `π < 2`,
que é a diferença que E1.3 tinha registrado.

- **A perda quadrática dispensa a localização do WALL.** A trilha rápida de
  lá precisa de `(M*)² J 2^J = o(n)` (o que restringe `τ < 2/3`) e de um
  piso em `s`, porque localiza o estimador não modificado. Aqui o Lema 8 é a
  prova de E1.5 relida com outro comparador: `τ ∈ (0,2)` inteiro e piso
  `s' > (2−τ)/4`, que o WALL só obtém para o estimador restrito. É vantagem
  técnica citável na comparação com o companheiro.
- Números: a razão do erro de predição ao alvo `(J_n/n)^{2s'/(2s'+1)}` vai
  de `2.95` a `3.70` com `n` de 200 a 12800 (máx/mín `1.25`; nas
  componentes, `1.72`); a dimensão efetiva `M*` fica de `51×` a `103×`
  abaixo da leitura densa `M = d`.
- **O pré-assintótico, medido de novo:** com `s' = 1/4` a regra leva `J_n` a
  7 só em `n = 12800`, e é só aí que a condição empírica de E1.4 vale em
  95% das réplicas (contra ≤10% em `n ≤ 6400`); a regra fica **dois níveis
  acima** do `J` que minimiza o erro realizado, ao custo de 5% a 11%. A
  Seção 5 do artigo tem de dizer isso.
- **Não há cota inferior**, e construir a de `q ≥ 2` com desenho aditivo é
  trabalho do porte de E1.4 mais E1.5, fora do plano. Ver pergunta 12.

**E2.2, o estimador** (`wafc/R/fit.R`, `reconstruct.R`, `tests/test-fit.R`).
`wafc()` com LASSO e sparse group LASSO (grupo = bloco `(ℓ,m)`), `coef`,
`predict`, `wafc_functions()`, `wafc_blocks()` e `wafc_kkt()`. As KKT fecham
em `74/74` pontos do caminho no LASSO e `100/100` nos grupos.

- **A escala de `λ` ficou fixada por KKT:** o `λ` do objeto é o do objetivo
  de E1.5, e o do `glmnet` é esse dividido por `nvars/npen` (medido:
  `1.016129` contra razão observada `1.016136`). O `sparsegl` não reescala.
- **Resposta com número à pergunta 3 da §4** (a variante em grupos entra em
  `wafc()` ou fica no piloto): **entra**. No mesmo desenho e com ~99 não
  nulos, o LASSO zera `0` dos 6 blocos e o sparse group LASSO zera
  exatamente os `3` blocos nulos verdadeiros; o preço é predição (`RMSE`
  `0.2633` contra `0.2078`, com `σ = 0.7314`) e tempo (`3,7 s` contra
  `0,10 s`). A escolha entre as duas continua sendo de E2.5.
- **Correção do achado de E2.1:** o `glmnet` só descarta a coluna constante
  quando o desenho é **denso**; no esparso ele a mantém e reparte o nível
  com o próprio intercepto. A rota "somar o intercepto ao nível" é a única
  que funciona nas duas formas, e bate a `1e-13` entre elas.
- Dois defeitos de motor contornados: o `sparsegl` devolve o vetor nulo no
  primeiro ponto do caminho (inclusive os `c_ℓ`, que não é solução), e
  `asparse = 1` não é jeito confiável de obter o LASSO por ele (violações de
  2% a 5% de `λ`).

**L3, segunda rodada bibliográfica.** `referencias-verificadas.bib` vai de
35 a **57 entradas**; `literatura.md` não tem mais nenhuma linha em `resumo`
nem `[VERIFICAR]`; os `derivations/` não têm mais marca bibliográfica
pendente. O que ela corrigiu, e que teria ido para o manuscrito errado:

- **Três anos e um autor de L2 estavam errados:** Deshpande et al. é **2026**
  (*Bayesian Anal.* 21(1) 281-308), não 2024; Zhou, Xu & Lin é **2017**, não
  2016; e o terceiro autor de Zhou, Yang & Xiang (2022) é **Xiang**. A
  origem dos dois erros de ano é o prefixo do DOI (`10.1214/24-…` é o ano de
  **aceitação**), e a lição é conferir `published-print` no Crossref, não o
  campo `issued`.
- **Duas citações de teorema do `04-oraculo.tex` estavam trocadas:** em
  Bühlmann & van de Geer (2011) a taxa lenta é o **Corolário 6.1**, não o
  Teorema 6.1, e as coordenadas não penalizadas são a **§6.9** ("The
  weighted Lasso"), não a §6.2.3. A matemática de E1.5 não muda, os rótulos
  sim. O Corolário 6.8 que E1.4 cita está certo. O mapa completo ficou em
  comentário na entrada `buhlmann2011statistics`, para não se refazer a
  busca em E5a.
- Sobra uma única marca `verificar` nos `derivations/`, e é **matemática**,
  não bibliográfica: o esboço da prova da extensão na Proposição 2 de E1.3.

### Decisões tomadas

| # | Data | Decisão | Razão |
|---|---|---|---|
| D1 | 09-18 | O modelo é `Y = Σ_j β_j(U) X_j + ε`, `β_j(u) = c_j + Σ_k g_{jk}(u_k)`, com `X_1 ≡ 1` permitido (o aditivo puro e o parcialmente linear aditivo são casos particulares) | pedido do autor |
| D2 | 09-18 | Base: wavelets ortonormais de suporte compacto do `WaveBased`; padrão periódico com `j0 = 0` e a função de escala constante descartada (identificabilidade no nível da base, como no `wall()`); `boundary = "interval"` como opção | herda o `wall()`; é o que faz a restrição `∫ g_{jk} = 0` sair de graça |
| D3 | 09-18 | Estimador base: LASSO sobre todos os coeficientes de wavelet, `c_j` não penalizados; sparse group LASSO por par `(j,k)` é a variante a medir em E2 | pedido do autor (LASSO); a variante em grupos é a candidata natural à seleção de estrutura |
| D4 | 09-18 | **O código do método vive na pasta `wafc/` deste repositório** (`R/`, `tests/`, `scripts/`); o `WaveBased` não recebe código por enquanto e é usado só como dependência para as bases (`wbasis()`, `wtable()`), confirmado pelo autor; **as funções criadas ficam privadas por enquanto** (neste repositório privado, sem pacote público, sem `install_github`); empacotar e publicar decide-se em E3.3, com o código testado | decisão do autor, contra a proposta de implementar dentro do `WaveBased` |
| D9 | 09-18 | **Índices:** wavelet `ψ_{jk}` (nível `j`, translação `k`); covariável linear `X_ℓ`, `ℓ = 1, …, p`; moduladora `U_m`, `m = 1, …, q` | decisão do autor: não mexer no padrão da base; a colisão sai das covariáveis |
| D10 | 09-18 | **Coeficientes em `θ`**, com `θ_{ℓm,jk}` (bloco antes da wavelet); `β_ℓ` fica sendo só o coeficiente funcional | `β` não pode ser função e vetor na mesma seção |
| D11 | 09-18 | **`U ∈ [0,1]^q` por hipótese populacional**, com densidade limitada longe de `0` e de `∞`; reescalonamento empírico só na seção de computação | teoria limpa em E1.3 a E1.6; o termo extra não vale o custo agora |
| D12 | 09-18 | Ordem das colunas de `Z`: não penalizados, depois blocos `(ℓ, m)` lexicográficos, dentro do bloco `j` e `k` crescentes | fixa a interface de `wafc_design()` em E2.1 |
| D13 | 09-18 | **A teoria assume o caso geral** `λ_min(E[XX' \| U]) ≥ κ_1 > 0` q.c.; `X ⊥ U` vira observação (a fatoração de Kronecker) | E1.4 fechou o caso geral com `λ_min(Σ) ≥ κ_1 c_U`, conferido a 10% da verdade; assumir independência custaria generalidade sem comprar constante |
| D14 | 09-18 | **D11 é sobre a densidade conjunta** de `U` em `[0,1]^q`, não sobre as marginais | `U_2 = U_1` tem marginais uniformes e `Σ_Ψ` singular; a prova de E1.4 usa `c_U` da conjunta, e é daí que sai a "não colinearidade entre moduladoras" |
| D15 | 09-19 | Padrões do `wafc_design()`, ratificados do handoff de E2.1: `filter.size = 8` (contra os 20 do `wall()`), nomes de coluna `x2:u1:psi3.5` e de bloco `x2:u1`, erro informativo em `j0 != 0` e `boundary = "interval"`, e cenários com `β_1` aditivo em duas moduladoras, `β_2` em uma e `β_ℓ` constante para `ℓ ≥ 3` | o filtro 8 é o das três conferências de `derivations/check/`; o erro em vez da implementação mantém a pergunta 6 aberta sem código morto; o cenário tem de exibir o termo cruzado, que L2 apontou como o eixo do artigo |
| D16 | 09-19 | **O enunciado principal do artigo é o Corolário 5 de E1.6** (taxa `n^{−2s/(2s+1)}` a menos de logaritmos, sob a mesma hipótese de Besov de E1.3), com o Teorema 2 como a taxa do sieve e a Proposição 4 como enunciado incondicional | o Lema 9(i) mostra que a compressibilidade não custa hipótese, e é a única das três taxas ótima para `π < 2`; **a ratificar, junto com a pergunta 7** |
| D17 | 09-19 | **Interface de `wafc()`** (E2.2): o `λ` do objeto é o do objetivo, não o do motor; `intercept` resolvido por presença de covariável constante, com erro informativo nos casos ambíguos; `coef()` dobra o intercepto no nível e `predict()` usa os coeficientes crus; mínimo quadrado escrito no ponto nulo do caminho; `wafc_kkt()` e `wafc_blocks()` públicas | a escala de `λ` é o que liga o código à teoria de E1.5, e as outras quatro saem dos defeitos de motor medidos; **a ratificar** |
| D6 | 09-18 | Documentos de trabalho em português; manuscrito em inglês americano; convenções de git, marcação e continuidade herdadas do `bdm-draft` | pedido do autor ("em linha com o bdm-draft") |
| D7 | 09-18 | Compêndio de simulação e aplicação em repositório próprio, `wafc-studies`, nos moldes do `wall` | o `wall` já resolveu cache, `renv` por commit e proveniência |

**Propostas de L2, a ratificar (2026-09-18).** Nenhuma bloqueia E1.2, E1.3 ou
E1.4; todas mudam documento, não resultado.

| # | Proposta | Onde |
|---|---|---|
| L2a | Frase-tese e contribuição 1 reescritas como extensão de Klopp & Pensky a coeficientes aditivos em várias moduladoras e desenho dependente; "Why not block LASSO?" entra nas perguntas do referee | `alvo-revista.md` §4 |
| L2b | E1.4 (i) vira "recordar K&P (eq. 1.8, Lema 1) e estender ao desenho aditivo"; o entregável central passa a ser o termo cruzado entre moduladoras e a parte (ii) | `plano-projeto.md` E1.4 |
| L2c | O QUT (Giacobino et al. 2017), que Sardy & Ma usam, entra em E2.3 como regra de `λ` sem `σ`, ao lado de BIC/EBIC | `plano-projeto.md` E2.3 |
| L2d | Concorrentes mínimos de E2.4/E4: `mgcv`, spline adaptativo (Wang, Jiang & Liu 2024), block LASSO de K&P no mesmo desenho, VCBART; cenário não homogêneo com as funções de Donoho-Johnstone | `plano-projeto.md` E2.4, E4 |
| L2e | "O mais próximo na teoria é Klopp & Pensky; no método, Sardy & Ma e Amato et al."; citar Montoril, Morettin & Chiann (2018) | `proposta-metodo.md` §5 |
| L2f | E1.7 fica opcional ou vira variante sem teorema de seleção: a ideia já existe em splines (Antoniadis et al. 2014; Ma et al. 2015) e em wavelets sem `X_ℓ` (Amato et al. 2022) | `plano-projeto.md` E1.7 |

**Adiadas pelo autor ("cobre-me depois", 2026-09-18); propostas do assistente:**

| # | Proposta | Razão | Alternativa |
|---|---|---|---|
| D5 | Alvo primário **Statistica Sinica**; reserva EJS | `alvo-revista.md` §2 | EJS como primário se a teoria pesar mais que a aplicação |
| D8 | Nome do método **WAFC** (*wavelet additive functional coefficients*) | curto, ecoa o WALL, e a sigla é a do repositório | "WAVC" (varying coefficients); "wavelet additive coefficient LASSO" |

---

## 3. O que NÃO reabrir

- **Pasta dentro do `wafc-draft`, `WaveBased` ou repositório próprio para o
  código.** Respondido em D4 e D7, detalhe em `plano-projeto.md` E0.2: o
  método em `wafc/` aqui, compêndio próprio. O empacotamento volta como
  decisão só em E3.3, com o código testado.
- **Se a restrição de identificabilidade precisa de multiplicador ou
  centralização numérica.** Não precisa no caso periódico com `j0 = 0`: a
  constante é a única função de escala e é descartada (D2). Só volta se
  `boundary = "interval"` virar o padrão.

---

## 4. Perguntas em aberto

Ordenadas pelo que bloqueia mais.

1. **D5 e D8** (tabela acima), adiadas. E1.1 fechou, que era o marco para
   cobrá-las: D5 fixa o formato de E5a e D8 é o nome que vai no título.
2. **Aplicação (E6.1):** o autor tem uma base em mente? Os candidatos de
   `plano-projeto.md` E6.1 são genéricos. Decidir cedo evita desenhar a
   simulação longe do caso real.
3. **~~Sparse group LASSO em `wafc()`~~ respondida por E2.2:** entrou como
   opção, com o número de seleção de estrutura (3 de 3 blocos nulos contra
   0 do LASSO) e o preço em predição. Qual das duas é a principal continua
   sendo de E2.5.
4. **Centralização das componentes:** o `notacao.md` §2 escreve
   `E[g_{ℓm}(U_m)] = 0` e a base impõe `∫_0^1 g_{ℓm} = 0`; o estimador
   estima a versão de Lebesgue, e as duas ficam longe quando `U_m` não é
   uniforme (na conferência de E1.2, `max |mean ψ_{jk}(U_2)| = 0.48` com
   `U_2 ~ Beta(2,3)`). A Proposição 1 vale com as duas, mas o `.tex` tem de
   dizer uma. Proposta: **Lebesgue**, como no WALL, com uma linha no
   `notacao.md` §2. É a única mudança de hipótese que E1.2 pede.
5. **Periodicidade no enunciado:** (a) hipótese em forma de sequência na
   base periodizada, como o WALL, com a proposição do custo da periodização
   como justificativa do reescalonamento e da opção `interval`; (b) teoria
   na base do intervalo (CDV), sem periodicidade, pagando `p q (2^{j_0} − 1)`
   parâmetros não penalizados; (c) as duas, teoria em (a) e `boundary` como
   opção do código. Proposta de E1.3: **(c)**, que é o que o `wall()` faz.
6. **`boundary = "interval"` no `wafc()` de E2.1:** entra com a
   reparametrização `Φ Q` (colunas de escala não penalizadas) ou fica de
   fora da primeira versão? Com o filtro padrão do `wall()` a opção exige
   `j_0 ≥ 6`, o que são `63 p q` parâmetros não penalizados.
7. **Como o artigo se posiciona diante de Klopp & Pensky** (a mais
   importante desta rodada): o WAFC se apresenta como extensão deles
   (honesto, e o referee da SS reconhece) ou como modelo diferente, com K&P
   citado como caso particular? A primeira muda a frase-tese de
   `alvo-revista.md` §4; a segunda obriga a Seção 2 a mostrar que o desenho
   aditivo não é o deles. Decide-se antes de E5a; E1.4 não depende.
8. **Block LASSO de K&P:** entra em `wafc()` como opção de penalidade, ao
   lado do sparse group LASSO de E2.2, ou fica só como concorrente em E2/E4?
9. **Idioma do código e das derivações.** E2.1 escreveu o Roxygen e os
   comentários de `wafc/` em inglês (é o que vira pacote em E3.3, e o
   `wall.R` é todo em inglês); os `derivations/` estão misturados, `02` e
   `04` em português e `03` em inglês. D6 fala de documentos e manuscrito,
   não de código. Proposta: **inglês no código**, **português nas
   derivações** (traduzir só quando E5a montar o manuscrito), e alinhar o
   `03`. Se aceita, entra uma correção de uma linha por ambiente no
   `macros.tex`, que hoje imprime "Lemma 4" enquanto a prosa diz "Lema 4".
10. **`rescale = TRUE` como padrão do `wafc_design()`?** É o que o `wall()`
   faz e o que a prática exige, mas desloca a centralização (item acima).
   Alternativa proposta por E2.1: manter o padrão e avisar quando a
   amplitude amostral já estiver dentro de `[0,1]`.
11. **~~Quem escolhe `J`~~ não era pergunta:** o plano sempre disse
   `cv.wafc()` sobre `(J, λ)`, como o `cv.wall()` (`plano-projeto.md` E2.3);
   o handoff de E1.6 leu o plano como se ele só falasse de `λ`. O que fica
   de E1.6 é **matéria de medição para E2.3**, não bloqueio: a regra teórica
   `J_n` do Teorema 2 ficou dois níveis acima do `J` que minimiza o erro
   realizado em toda a varredura, ao custo de 5% a 11% de erro, e as
   escolhas de `J_n` e de `c` dependem de `s'` e `τ`, que ninguém conhece.
   E2.3 compara a regra teórica com a validação cruzada.
12. **Cota inferior.** Não existe aqui, e Klopp & Pensky têm a deles para
   `q = 1` com `X ⊥ U`. Três saídas: (a) citar K&P e dizer que a cota
   superior atinge a referência do modelo de sequência, que é o que o
   `05-taxas.tex` faz hoje; (b) abrir E1.8 e construir a cota para `q ≥ 2`,
   trabalho do porte de E1.4 mais E1.5; (c) restringir a afirmação de
   otimalidade a `q = 1`. O referee da SS pode cobrar a (b).
13. **`p` crescente com `n`.** A teoria fixa `p` e `q`. Se a aplicação de E6
   tiver `p` grande, o termo `σ² p / n` deixa de ser de ordem menor e a
   janela do Corolário 5 estreita; mudar isso mexe em E1.3 e E1.4, não só em
   E1.6.
14. **Bibliografia, pontos de L3** (nenhum bloqueia): a identidade
   `Σ_l φ(x − l) ≡ 1`, usada na prova do Lema 1(i) de E1.2, ficou **sem
   âncora** — a expressão não ocorre em Daubechies (1992), a quem estava
   atribuída, e os candidatos a conferir são Härdle et al. (1998, cap. 5) e
   o Mallat; a citação de Meyer (1992, cap. III §11) para a base
   periodizada não foi conferida (a CUP bloqueia acesso automatizado) e foi
   substituída no texto pela de Daubechies (1992, §9.3), que foi; e
   `hardle1998wavelets`, copiada do WALL, grafa "Tsybakov, Alexandre" onde
   o Crossref diz "Alexander" (uniformizar nos dois repositórios ou deixar?).
   Proposta: abrir uma frente curta só para a âncora da partição da unidade
   quando E5a precisar dela.
15. **Bibliografia, quatro pontos deixados por L1** (nenhum bloqueia; o
   `.bib` fica como está até a resposta):
   - Amato et al. (2022) ou Haris, Simon & Shojaie (2018) na linha que
     citava o arXiv 1903.04631? Proposta: **as duas**, que são trabalhos
     distintos e ambos interessam a L2.
   - Hastie & Tibshirani (1993): páginas 757-779 (só o artigo) ou 757-796
     (com a discussão)? Proposta: **757-779**, que é o que o Crossref
     registra.
   - Daubechies & Lagarias: as duas partes ou só a parte I, que é a do
     algoritmo do `WaveBased`? Proposta: **só a parte I**, e a parte II
     entra se alguma prova precisar dela.
   - Chaves do `.bib`: padrão único `Autor-Autor-Ano` ou manter as quatro
     herdadas do WALL (`cohen1993wavelets` etc.)? Proposta: **manter as do
     WALL**, porque o `ms_theo_1.tex` é o molde da prova e a citação
     cruzada fica direta.

---

## 5. Próximos passos

(a) e (b) correm em paralelo; (c) espera E1.2.

- (a) **Autor:** ratificar D5 e D8 (a notação já não depende disso);
  responder a pergunta 2 se já tiver a base da aplicação.
- (b) **Chat de tarefa, pode abrir já:** E2.3 (sintonia), catalogada com
  a comparação entre a regra teórica de `J_n` e a validação cruzada.
- (b') **L3 fechada.** Quando o autor ratificar as propostas de L2, entram
  as edições de `alvo-revista.md` e `plano-projeto.md` da §6 de
  `busca-novidade.md`, que são do chat principal.
- (c) **Chat de tarefa, depois do enunciado de E1.2:** E2.1, que precisa
  saber o que a identificabilidade descarta do desenho.
- (d) **Chat principal:** integrar os handoffs e catalogar E1.5, E1.6 e
  E2.2 a E2.5 quando as dependências fecharem.
- (e) **Depois de E1.3 e E1.4:** E1.5, E1.6; E2.2 a E2.5.

## 6. Histórico de sessões

| Data | O que aconteceu |
|---|---|
| 2026-09-18 | Avaliação de viabilidade; criação do repositório e dos documentos de trabalho; template da EJS; plano E0 a E7 |
| 2026-09-18 | D4 decidida pelo autor (código em `wafc/`, não no `WaveBased`); D5 e D8 adiadas; `prototype/` virou `wafc/`; plano E2 e E3 reescritos; repositório publicado; o autor confirmou o `WaveBased` como dependência e que as funções ficam privadas |
| 2026-09-19 | E1.6, E2.2 e L3 fechadas e integradas: **E1 inteira**, com a compressibilidade saindo de graça da hipótese de Besov (D16); `wafc()` com as duas penalidades e KKT fechando no caminho inteiro (D17); `.bib` com 57 entradas e duas citações de teorema corrigidas |
| 2026-09-19 | E1.5 e E2.1 fechadas em dois chats de tarefa e integradas: desigualdade oráculo sem cone, com a razão `‖f̂−f‖²_n/(λ²s_0)` estável a 5% em `n` de 200 a 6400; `wafc/` nasce com 103 testes passando; D15 |
| 2026-09-18 | E1.2, E1.3 e E1.4 fechadas em três chats de tarefa e integradas; as três conferências rodam e imprimem `OK` aqui; D13 e D14; a numeração global dos resultados fixada |
| 2026-09-18 | L2 fechada em chat de tarefa e integrada: novidade central confirmada (interseção zero na SS e na EJS), mas Klopp & Pensky (2015) cobre E1.4 (i), E1.5 e E1.6 para `q = 1` e `X ⊥ U`; seis propostas de mudança de rumo a ratificar |
| 2026-09-18 | L1 fechada em chat de tarefa e integrada: 35 referências verificadas, quatro correções de atribuição, dois trabalhos novos para L2 olhar |
| 2026-09-18 | E1.1 fechada: notação congelada (`ψ_{jk}`, `X_ℓ`, `U_m`, `θ_{ℓm,jk}`, `U ∈ [0,1]^q`), D9 a D12, `macros.tex` reescrito e compilando |
| 2026-09-18 | Máquina nova conferida (R 4.6.1, `WaveBased` 2.6-0, `grpreg`, `gglasso`, `sparsegl` por `apt`); E0.3 fechada: instruções da SS transcritas, templates versionados e compilando, teto corrigido de 30 para 40 páginas |
