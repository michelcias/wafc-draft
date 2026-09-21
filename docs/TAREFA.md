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
| E2.4b emendas do piloto | **fechada** (2026-09-21): cenário `uneven`, `s'` como atributo, QUT em `tune.R`, `k` casado e `bam` no `gam`, tolerância e `...` consertados; 642 testes passam | `wafc/R/dgp.R`, `wafc/R/tune.R`, `wafc/scripts/06-timing.R` |
| E1.10 terminologia | **fechada** (2026-09-21): "sieve" sai das derivações, uma menção retida | `derivations/*.tex` |
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

Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**O catálogo está vazio.** A próxima é E2.5 (go/no-go), e ela **não deve
abrir antes de três coisas**: o conserto do repasse de `wavelet.table` aos
concorrentes, que hoje faz a coluna do `vcbart` sumir em silêncio; a decisão
sobre P1 e P2, que E2.4b devolveu sem ratificação; e a repetição da parte
`competitors` com `k` casado no `gam` (`wafc_k_matched()` e motor `bam`, já
disponíveis).

**Ao catalogar, o arquivo vai junto do item.** Duas tarefas seguidas
esbarraram em coluna de arquivos que não cobria o que o próprio texto
mandava fazer: E2.4 pôs o QUT em `competitors.R` por isso, e E2.4b teve de
tocar `fit.R` e criar `06-timing.R` fora da coluna. **E2.5 (go/no-go) só deve abrir depois de E2.4b**, por
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
