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
| E0.3 template e instruções | **fechada** (quartil Q1 confirmado em 2026-09-19) | `manuscript/ejs-template/`, `manuscript/ss-template/`, `ss-instrucoes-autores.md` |
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
| E2.3 sintonia | **fechada** (2026-09-19): `cv.wafc()`, BIC, EBIC e a regra da teoria; 402 testes passam | `wafc/R/tune.R`, `wafc/scripts/02-tune.R` |
| E5a manuscrito, Seções 1 a 4 | **fechada** (2026-09-19): `k = 1` compila limpo, 28 e 26 páginas | `manuscript/ms_1.tex`, `supp_1.tex` |
| E1.3b emenda de extensão | **fechada** (2026-09-19): Lema 10 e três observações; conferência `OK` | `derivations/02-aproximacao-besov.tex` |
| E2.1b conserto do `eps` | **fechada** (2026-09-19): margem fixa, `J = 1` acessível; 420 testes passam | `wafc/R/design.R` |
| L4 referências de software | **fechada** (2026-09-19): `.bib` com 63 entradas | `referencias-verificadas.bib` |
| E1.8 rota do intervalo | **fechada** (2026-09-19): Lemas 11 e 12, Proposições 5 e 6, Corolários 6 e 7; conferência `OK` | `derivations/07-rota-intervalo.tex` |
| E2.4 piloto | **fechada** (2026-09-20): sete concorrentes, seis regras de `λ`, varredura de `eps`; 544 testes passam | `wafc/R/competitors.R`, `wafc/scripts/04-pilot.R` |
| E2.4b emendas do piloto | não aberta; **pode abrir já** | catálogo da §3 |
| E1.7a sondagem de irrepresentabilidade | **fechada** (2026-09-20): veredito de escopo reduzido (D32) | `derivations/06a-sondagem-irrepresentabilidade.md` |
| E1.7c seleção por limiarização | **fechada** (2026-09-20): Lema 13 e Corolário 8; conferência `OK` | `derivations/06-selecao-limiar.tex` |
| E1.4c tradução do `03` | **fechada** (2026-09-20) | `derivations/03-desenho-produtos.tex` |
| E6.1a sondagem de bases | **fechada** (2026-09-20): veredito negativo nas três candidatas | `docs/aplicacao-candidatas.md` |
| E2.5, E3, E4, E5b, E6, E7 | não abertas | |
| L1 verificação bibliográfica | **fechada** (2026-09-18): 35 entradas verificadas | `referencias-verificadas.bib`, `literatura.md` |
| L2 busca de novidade | **fechada** (2026-09-18): novidade confirmada, Klopp & Pensky (2015) é o vizinho | `busca-novidade.md`, `literatura.md` |

Esta tabela é atualizada pelo chat principal quando uma etapa fecha; em caso
de dúvida, o `ESTADO.md` manda.

## 3. Catálogo das tarefas abertas

Cada linha diz o que entregar, de que depende, e **os únicos arquivos que a
tarefa pode criar ou editar**.

| Tarefa | Entregável | Depende de | Arquivos permitidos |
|---|---|---|---|

| E2.4b | as emendas que E2.4 não podia fazer, porque os arquivos estavam fora do catálogo dela. **(i)** `wafc/R/dgp.R`: gravar `s'` como atributo do cenário (D27: `3/2` no suave, `1/2` no não homogêneo, com o **regime** anotado, porque com margem o suave leria `4`), hoje só em `wafc_sprime` no alto de `04-pilot.R`; e **acrescentar a componente suave de curvatura desigual** que D30 pede (gaussiana estreita ou `doppler` truncado longe da singularidade), `C^∞` e portanto dentro da hipótese de Xue & Yang, mas com escala variando ao longo do domínio — é ela que testa a afirmação de D30, e é pré-requisito de E2.5. **(ii)** `wafc/R/tune.R`: mover `wafc_lambda_qut()` de `competitors.R` para cá e expô-la como sexta regra de `wafc_tune(rule = "qut")`, que é onde ela pertence (E2.4 a escreveu em `competitors.R` por restrição de catálogo). **(iii)** `wafc/R/load.R`: `VCBART` entra em `wafc_suggests`. **(iv)** O `gam` como concorrente: `k` por moduladora e ajuste por `bam`, porque E6.1a mostrou que a comparação em dimensão **não** casada é o que separa "o WAFC ganha 6 a 15%" de "empata", e o piloto rodou com o padrão `k = 10`. **(v)** Se o autor ratificar P1 e P2 (perguntas 8 e 9 da §4 do `ESTADO.md`), aplicar as duas linhas: grade de `J` até `⌈log_2 n⌉` em `wafc_J_grid()` e `wafc_eps_periodic` de `0.05` para `0`; se não ratificar, deixar como estão e dizer no handoff. **(vi) Dois defeitos de desempenho, medidos no chat principal em 2026-09-21, que são o que faz o WAFC parecer caro na tabela do piloto.** O padrão `thresh = 1e-10` de `wafc()` custa **28×** no candidato dominante (`J = 5`, `n = 1000`: validação cruzada de 5,95 s contra 0,21 s com `1e-7`) **sem mudar a sintonia** (`lambda.min` idêntico, `cvm` mínimo diferente na quinta casa); o `glmnet` sozinho, no mesmo desenho, leva `0,01 s`. Rever o padrão (a tolerância apertada pertence às conferências de KKT, não ao ajuste de produção) e **consertar o roteamento do `...` de `cv.wafc()`**, que hoje manda os argumentos para `wafc_design()` e faz `cv.wafc(..., thresh = 1e-7)` parar com "argumento não utilizado". **(vii) `wafc/scripts/06-timing.R`**, curto e só de tempo: um cenário, dois ou três `n`, todos os concorrentes, com o tempo **separado em sintonia e ajuste final** — hoje o tempo do WAFC inclui a validação cruzada inteira e o do `gam` é um ajuste só, com os parâmetros de suavização escolhidos por REML por dentro, e a tabela do piloto não diz isso. O `gam` medido em `k = 10` **e** em `k` casado à dimensão do sieve, que é a comparação honesta de E6.1a e tem custo diferente; com a opção de `bam` do item do `competitors.R`, medir as duas. A bateria inteira roda antes do handoff | E2.4, fechada; D27, D30, e P1/P2 conforme a resposta do autor | `wafc/R/dgp.R`, `wafc/R/tune.R`, `wafc/R/load.R`, `wafc/R/design.R` (só a linha da margem), `wafc/R/competitors.R` (a remoção do QUT **e** o `k` do `gam`: passar a aceitar `k` por moduladora, hoje um só para todos os suavizadores, e ganhar opção de ajuste por `mgcv::bam` com `fREML` e `discrete = TRUE`, que é duas ordens de grandeza mais rápido nestes `n` — E6.1a §3.2), `wafc/tests/test-dgp.R`, `wafc/tests/test-tune.R`, `wafc/tests/test-design.R`, `docs/handoff-E2.4b.md` |
Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**Nada em curso. Aberta no catálogo: E2.4b**, que recolhe as emendas que
E2.4 não podia fazer. **E2.5 (go/no-go) só deve abrir depois de E2.4b**, por
três razões medidas: o veredito contra os concorrentes muda com a grade de
`J` (P1); o cenário suave que D30 quer medir depende da componente nova em
`dgp.R`; e a comparação com o `gam` foi feita em dimensão não casada, que é
o que E6.1a mostrou separar empate de vitória.

**Quem abrir E1.7a** entrega um **veredito**, não um teorema: a tarefa
existe para decidir se vale gastar semanas na saída (a), e um "não fecha,
e eis o contra-exemplo" é resultado tão bom quanto um "fecha".

 Não compartilham arquivo: E1.8 está em
`derivations/`, E2.4 em `wafc/`. A numeração "E1.8" é desta tarefa; a cota
inferior, que o handoff de E1.6 chegou a chamar de E1.8, será **E1.9** se
algum dia for aberta.

**Quem abrir E1.8** escreve teoria que **não vai ao manuscrito agora** (D23,
D26): o artigo enuncia na base periodizada com `eps = 0`. O arquivo existe
para estar pronto se um referee pedir, e para guardar a análise de `s'` que
E2.4 e E4 vão citar.

**Quem abrir E2.4** herda um conserto de três linhas em `wafc/R/tune.R`, que
descrevem o padrão de `eps` anterior a E2.1b e ficaram falsas. A próxima depois dela é E2.4 (piloto), que espera a
ratificação de D19 e D20 e a resposta sobre `s'` de cada cenário (pergunta
16). Tarefa nova pede catalogação antes de abrir.

**Quem abrir E1.3b** lê antes a §7 do `notacao.md` (D22) e a entrada D23 da
tabela de decisões: a margem `eps` deixou de ser conveniência numérica e
passou a ter função declarada, que é comprar a extensão. O arquivo é de uma
etapa fechada; a emenda acrescenta, não reescreve, e a Proposição 2 continua
onde está, agora com o caso `eps = 0` explicitado.

**O manuscrito está vivo em `k = 1`.** Qualquer alteração nele passa pela
regra do `CLAUDE.md`: perguntar antes se deve nascer `k = 2`, e `ms`, `supp`
e `references` andam juntos. O mapa da numeração global para os rótulos do
LaTeX está no cabeçalho do `ms_1.tex`.

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
| `02-aproximacao-besov.tex` (emenda E1.3b) | Lema 10 |
| `07-rota-intervalo.tex` (E1.8) | Lema 11, Proposição 5, Corolário 6, Proposição 6, Corolário 7, Lema 12 |
| `06-selecao-limiar.tex` (E1.7c) | Lema 13, Corolário 8 |

**Este mapa é a autoridade**, e é ele que E5a usa ao montar o manuscrito. O
`04-oraculo.tex` já imprime o número global (via `\setcounter` no
preâmbulo), prática adotada daqui em diante; `02` e `03` ainda imprimem o
contador local. O próximo resultado numera a partir de Proposição 6, **Lema 13**, Teorema 2
e **Corolário 8**. O `06-selecao-limiar.tex` imprime `Hipótese S`, com letra
em vez de número, porque é citada lado a lado com a `Hipótese 1` de E1.6;
é desvio local e aceito.

**O enunciado que vai ao resumo do artigo é o Corolário 5** (D16, a
ratificar).

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
