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
| E1.6, E2.2, L3 | não abertas; **podem abrir já** | catálogo da §3 |
| E1.7, E2.3 a E2.5, E3 a E7 | não abertas | |
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

| E1.6 | `derivations/05-taxas.tex` (+ `.pdf`): com `J_n ≍ log_2 n/(2s'+1)`, taxa de predição `n^{−2s'/(2s'+1)}` a menos de logaritmos; erro `L_2` de cada `ĝ_{ℓm}` pelo Corolário 3 de E1.5; corolário de compressibilidade em weak-`ℓ_τ`, que é o que troca `s_0` por dimensão efetiva (sob Besov genérico `s_0 = d`, e sem ele o enunciado incondicional é só a taxa lenta); transposto da trilha rápida do WALL. `check/05-taxas.R` medindo o erro de predição contra `n^{−2s'/(2s'+1)}` ao variar `n`, e a dimensão efetiva sob coeficientes compressíveis | E1.5, fechada (Corolários 2 e 3, Lema 7); o regime de `J_n` tem de satisfazer a condição empírica de E1.4 | `derivations/05-taxas.tex`, `derivations/05-taxas.pdf`, `derivations/check/05-taxas.R`, `docs/handoff-E1.6.md` |
| E2.2 | `wafc/R/fit.R`, `wafc/R/reconstruct.R`, `wafc/tests/test-fit.R`: `wafc()` com LASSO (`glmnet`, gaussiano) e a variante sparse group LASSO (`sparsegl`; grupo = bloco `(ℓ, m)`, que é o que D12 deixa contíguo); `predict`, `coef`, reconstrução de `ĝ_{ℓm}` e `β̂_ℓ` pelo `spec` de `wafc_design()`. **Ajustar com `intercept = TRUE` e somar o intercepto ao nível da covariável constante** (o `glmnet` descarta colunas de variância zero; `wafc_design()` devolve `constant` com os índices), e conferir por KKT em vez de confiar no `penalty.factor`, que o `glmnet` reescala para somar `nvars` | E2.1, fechada; o pacote da variante em grupos vai ao `CONTINUAR.md` na mesma rodada | `wafc/R/fit.R`, `wafc/R/reconstruct.R`, `wafc/tests/test-fit.R`, `wafc/README.md` (só a tabela), `docs/handoff-E2.2.md` |
| L3 | segunda rodada de verificação bibliográfica, com três frentes: (a) as linhas que L2 deixou em `literatura.md` com status `resumo` ou `[VERIFICAR]`, conferidas no Crossref como em L1 e promovidas a `verificado` com entrada no `.bib`; (b) as referências que E1.3, E1.4 e E1.5 citam e que não estão em `referencias-verificadas.bib` (`hardle1998wavelets` e `tropp2012user` são cópia de `../../wall-manuscript/manuscript/theo/references_theo_1.bib`, sem redigitar; Meyer 1992, Cohen 2003, Daubechies 1992, Triebel 1983 e Tibshirani 2013 são novas); (c) as marcas `[verificar]` de **numeração de teorema** deixadas nos `derivations/` (Bühlmann & van de Geer 2011, Teoremas 6.1 e 6.2, §6.2.3 e Corolário 6.8; Härdle et al. 1998 cap. 9; Cohen 2003 cap. 3; Triebel 1983), conferidas na edição citada e substituídas pelo número certo ou pelo enunciado correto. O que não for encontrado fica marcado, com o motivo | nada; L1 e L2 fecharam | `docs/referencias-verificadas.bib`, `docs/literatura.md` (só status e dados bibliográficos), `derivations/02-aproximacao-besov.tex`, `derivations/03-desenho-produtos.tex`, `derivations/04-oraculo.tex`, `derivations/01-identificabilidade.md` (só as marcas de verificação e as referências; **nada de enunciado, hipótese ou prova**), os `.pdf` correspondentes recompilados, `docs/handoff-L3.md` |
Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**Quem abrir E1.6** lê antes o Teorema 1 e os Corolários 2 e 3 de
`derivations/04-oraculo.tex`, e a §1 de
[`busca-novidade.md`](busca-novidade.md): a taxa adaptativa em Besov de
Klopp & Pensky (2015) cobre `q = 1` com `X ⊥ U`, inclusive com cota
inferior, e o que é novo aqui é o desenho aditivo. Atenção à calibração:
`‖Z‖_max` é `max_a sqrt(Σ̂_aa)`, que é `O_p(1)`, e não `max_{i,a}|Z_{ia}|`,
que é de ordem `2^{J/2}` e destruiria a taxa.

**Quem abrir E2.2** parte de `wafc_design()` como está e não mexe nele; os
dois achados que economizam meio dia estão no catálogo acima (coluna
constante e `penalty.factor`).

**L3 encosta nos `.tex` de E1.3 a E1.5**, que estão fechados: pode corrigir
número de teorema citado, dado bibliográfico e marca de verificação, e nada
mais. Enunciado, hipótese ou prova que pareçam errados viram item do
handoff. Se E1.6 estiver aberta em paralelo, ela não toca esses arquivos, e
L3 não toca o `05-taxas.tex`.

Tarefa que não está no catálogo: pedir ao chat principal para catalogá-la
antes de abrir. E1.5, E1.6, E2.2 a E2.5 e E3 entram no catálogo quando as
dependências fecharem. Todo código de tarefa vai para `wafc/`, nunca para o
`WaveBased` (D4).

## 4. Regras do chat de tarefa

- **Só o que foi pedido.** A tarefa é a nomeada na mensagem de abertura; o
  que se nota fora dela vira item no handoff, não trabalho.
- **Só os arquivos permitidos.** Não editar `ESTADO.md`, `plano-projeto.md`,
  `notacao.md`, este arquivo, nem arquivos de outra tarefa. Símbolo novo,
  decisão ou mudança fora do escopo: propor no handoff.
- **Conferência numérica antes da prova** (`instrucoes.md`, §5): o script em
  `derivations/check/` roda e imprime `OK` antes de o enunciado ser escrito.
- **Não commitar nem fazer push**, salvo pedido explícito no chat. Se pedido,
  `git add` por caminho explícito, nunca `-A`, e sem coautoria.
- **Terminar escrevendo `docs/handoff-<tarefa>.md`** com quatro seções: o que
  foi feito (arquivos); o que a conferência mostrou (números); o que deve
  entrar no `ESTADO.md` (decisões, sinais, lições); perguntas em aberto. O
  chat principal integra e apaga o handoff.
- Resposta final no chat: veredito, os números que importam, o comando de
  commit sugerido por caminho explícito.

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

**Este mapa é a autoridade**, e é ele que E5a usa ao montar o manuscrito. O
`04-oraculo.tex` já imprime o número global (via `\setcounter` no
preâmbulo), prática adotada daqui em diante; `02` e `03` ainda imprimem o
contador local. O próximo resultado numera a partir de Proposição 3, Lema 7,
Teorema 1 e Corolário 3.

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
