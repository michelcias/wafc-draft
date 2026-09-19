# Executar uma tarefa

Arquivo de abertura para qualquer chat que vá executar **uma** etapa ou
subetapa do plano. A mensagem de abertura é uma linha:

> Leia `docs/TAREFA.md` e execute a tarefa L2.

Nada mais precisa ser dito; o que a tarefa entrega, em quais arquivos, e o
que ler antes está aqui e no plano.

---

## 1. Ler antes, nesta ordem

1. [`../CLAUDE.md`](../CLAUDE.md): as duas regras (sem coautoria em commit;
   respostas curtas).
2. [`instrucoes.md`](instrucoes.md): fluxo, escrita, matemática, código, git.
   A §7 (chats paralelos) vale para todo chat de tarefa.
3. [`ESTADO.md`](ESTADO.md): onde o trabalho está, decisões, o que não
   reabrir, perguntas em aberto.
4. A seção da tarefa em [`plano-projeto.md`](plano-projeto.md): entregável e
   critério de saída.
5. O que a tarefa toca: [`notacao.md`](notacao.md) se houver fórmula;
   [`proposta-metodo.md`](proposta-metodo.md) para a ideia; os
   `derivations/` de que ela depende; `wafc/` se for código;
   [`inventario-codigo.md`](inventario-codigo.md) e `wafc/README.md` se for código;
   [`alvo-revista.md`](alvo-revista.md) se for manuscrito. Se for prova, o
   passo correspondente do WALL teórico
   (`../../wall-manuscript/manuscript/theo/ms_theo_1.tex`).

## 2. O que já foi feito

| Etapa | Estado | Onde está |
|---|---|---|
| E0.1 nomes | D4 fechada (código em `wafc/`); D5, D8 adiadas | `plano-projeto.md` E0.1, `ESTADO.md` |
| E0.2 onde o código vive | fechada (D4, D7) | `plano-projeto.md` E0.2 |
| E0.3 template e instruções | fechada, menos o quartil da SS no SCImago | `manuscript/ejs-template/`, `manuscript/ss-template/`, `ss-instrucoes-autores.md` |
| E0.4 ferramentas | fechada; `WaveBased` 2.6-0, `grpreg`, `gglasso`, `sparsegl` e `gh` instalados | `inventario-codigo.md` §4, `CONTINUAR.md` |
| E0.5 inventário | fechada | `inventario-codigo.md` |
| E0.6 convenções do pacote | fixadas | `plano-projeto.md` E0.6 |
| E1.1 notação | **fechada** (D9 a D12): `ψ_{jk}`, `X_ℓ`, `U_m`, `θ_{ℓm,jk}` | `notacao.md`, `derivations/macros.tex` |
| E1.2 identificabilidade | **fechada** (2026-09-18): Proposição 1 e Lema 1; conferência `OK` | `derivations/01-identificabilidade.md` |
| E1.3 aproximação em Besov | **fechada** (2026-09-18): dois lemas, corolário e o custo da periodização; conferência `OK` | `derivations/02-aproximacao-besov.tex` |
| E1.4 desenho de produtos | **fechada** (2026-09-18): autovalor mínimo cheio no caso geral (D13); conferência `OK` | `derivations/03-desenho-produtos.tex` |
| E1.5 oráculo | **fechada** (2026-09-19): Lemas 4 a 7, Teorema 1, Corolários 2 e 3; conferência `OK` | `derivations/04-oraculo.tex` |
| E2.1 desenho e cenários | **fechada** (2026-09-19): `wafc_design()`, cenários, 103 testes passando | `wafc/R/`, `wafc/tests/`, `wafc/scripts/` |
| E1.6 taxas | **fechada** (2026-09-19): Lemas 8 e 9, Proposição 4, Teorema 2, Corolários 4 e 5; conferência `OK` | `derivations/05-taxas.tex` |
| E2.2 estimador | **fechada** (2026-09-19): `wafc()` com as duas penalidades; 231 testes passam | `wafc/R/fit.R`, `wafc/R/reconstruct.R` |
| L3 bibliografia, 2ª rodada | **fechada** (2026-09-19): 57 entradas verificadas; duas citações de teorema corrigidas | `referencias-verificadas.bib`, `literatura.md` |
| **E1** | **fechada** (E1.7 é condicional a E2.5) | `derivations/` |
| E2.3 sintonia | não aberta; **pode abrir já** | catálogo da §3 |
| E5a manuscrito, Seções 1 a 4 | não aberta; **pode abrir já** (D5, D8, D16 e D18 decididas) | catálogo da §3 |
| E1.7, E2.4, E2.5, E3, E4, E5b, E6, E7 | não abertas | |
| L1 verificação bibliográfica | **fechada** (2026-09-18): 35 entradas verificadas | `referencias-verificadas.bib`, `literatura.md` |
| L2 busca de novidade | **fechada** (2026-09-18): novidade confirmada, Klopp & Pensky (2015) é o vizinho | `busca-novidade.md`, `literatura.md` |

Esta tabela é atualizada pelo chat principal quando uma etapa fecha; em caso
de dúvida, o `ESTADO.md` manda.

## 3. Catálogo das tarefas abertas

Cada linha diz o que entregar, de que depende, e **os únicos arquivos que a
tarefa pode criar ou editar**.

| Tarefa | Entregável | Depende de | Arquivos permitidos |
|---|---|---|---|
| E0.3 (autor) | só falta o quartil 2024 da SS conferido no SCImago e anotado em `alvo-revista.md` §1; a transcrição e o template estão feitos | nada | `docs/alvo-revista.md` (§1), `docs/handoff-E0.3.md` |

| E2.3 | `wafc/R/tune.R`, `wafc/tests/test-tune.R`: `cv.wafc()` sobre `(J, λ)` nos moldes do `cv.wall()` (dobras fixas, caminho de `λ` por `J`, `lambda.min` e `lambda.1se`), mais BIC e EBIC com graus de liberdade = número de coeficientes não nulos; as três regras comparadas nos três cenários de E2.1 em `n ∈ {250, 1000}`, com erro de predição fora da amostra e ISE por bloco. **Acrescentar uma quarta coluna à comparação: a regra teórica `J_n = min{J : 2^J ≥ (n/log n)^{1/(2s'+1)}}` do Teorema 2 de E1.6**, que na conferência daquela etapa ficou dois níveis acima do `J` de menor erro realizado, ao custo de 5% a 11%; medir quanto isso custa nos cenários. O `λ` de todo objeto é o do objetivo, não o do motor (D17), e `wafc_kkt()` confere | E2.2, fechada; `wafc()` e `wafc_blocks()` como estão | `wafc/R/tune.R`, `wafc/tests/test-tune.R`, `wafc/README.md` (só a tabela), `docs/handoff-E2.3.md` |
| E5a | `manuscript/ms_1.tex`, `manuscript/supp_1.tex` e `manuscript/references_1.bib`: o manuscrito nasce em `k = 1`, no template da *Statistica Sinica* (D5; copiar `manuscript/ss-template/SS-template-bib.tex` e `supp-temp_20240820.tex`, que compilam), com as Seções 1 a 4 de `alvo-revista.md` §5, isto é Introduction, Model and wavelet sieve, Theory (só enunciados; as provas vão ao `supp_1.tex`) e Computation and tuning. O método se chama **WAFC** (D8). A Introduction usa **como está** o parágrafo de posicionamento aprovado em `alvo-revista.md` §4 (D18: o artigo é extensão de Klopp & Pensky), e a Seção 3 abre pelo **Corolário 5** de `05-taxas.tex` (D16), com o Teorema 2 como taxa do sieve e a Proposição 4 como enunciado incondicional. Os enunciados são os de `derivations/01` a `05`, na numeração global da §5 deste arquivo, renumerados por rótulo do LaTeX; a notação é a de `notacao.md`, congelada; `references_1.bib` sai de `docs/referencias-verificadas.bib` (57 entradas verificadas), copiando só o que for citado. Critério de saída: `latexmk -pdf` limpo nos dois, dentro das 40 páginas em espaço duplo com referências, e toda afirmação de taxa ou de condição com o resultado correspondente citado | E1 fechada; L2 e L3 fechadas; D5, D8, D16, D18 | `manuscript/ms_1.tex`, `manuscript/supp_1.tex`, `manuscript/references_1.bib`, os `.pdf` correspondentes, `docs/handoff-E5a.md` |
Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**Quem abrir E5a** lê antes a §4 e a §5 de
[`alvo-revista.md`](alvo-revista.md) e a §3 de
[`ss-instrucoes-autores.md`](ss-instrucoes-autores.md), e obedece às regras
de escrita da §4 de [`instrucoes.md`](instrucoes.md): inglês americano, sem
negrito no corpo, sem em-dash. Três coisas que a rodada de 2026-09-19 deixou
prontas e que não se reabrem: o posicionamento (D18) com o parágrafo já
escrito, o enunciado de entrada (D16) e o nome (D8). Quatro perguntas da §4
do `ESTADO.md` ainda tocam a Seção 2 (centralização, periodicidade,
`boundary`, reescalonamento): enquanto não forem respondidas, escrever o que
não depende delas e levar as alternativas ao handoff, sem decidir.

**Quem abrir E2.3** não toca `fit.R` nem `design.R`; trabalha com
`wafc()` como está. Duas coisas da rodada de 2026-09-19 a governam: o `λ` de
todo objeto é o do objetivo de E1.5, não o do motor (D17, e o `glmnet`
reescala o dele por `nvars/npen`), e `wafc_kkt()` é o que confere isso em
vez de se confiar no `penalty.factor`.

## 5. Numeração dos resultados

Os `derivations/` numeram os resultados em sequência global (Proposição 1,
Lema 1, ...), começando em E1.2. O manuscrito terá numeração própria, por
rótulo do LaTeX. Resultado novo numera a partir do último registrado aqui
pelo chat principal, com o mapa abaixo.

| Arquivo | Resultados, na numeração global |
|---|---|
| `01-identificabilidade.md` (E1.2) | Proposição 1, Lema 1 |
| `02-aproximacao-besov.tex` (E1.3) | Lema 2, Lema 3, Corolário 1, Proposição 2 |
| `03-desenho-produtos.tex` (E1.4) | Proposição 3 |
| `04-oraculo.tex` (E1.5) | Lema 4, Lema 5, Lema 6, Lema 7, Teorema 1, Corolário 2, Corolário 3 |
| `05-taxas.tex` (E1.6) | Lema 8, Lema 9, Proposição 4, Teorema 2, Corolário 4, Corolário 5 |

**Este mapa é a autoridade**, e é ele que E5a usa ao montar o manuscrito. O
`04-oraculo.tex` já imprime o número global (via `\setcounter` no
preâmbulo), prática adotada daqui em diante; `02` e `03` ainda imprimem o
contador local. O próximo resultado numera a partir de Proposição 4, Lema 9,
Teorema 2 e Corolário 5.

**O enunciado que vai ao resumo do artigo é o Corolário 5** (D16, a
ratificar).

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
