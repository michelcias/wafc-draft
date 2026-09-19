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
| E1.8 rota do intervalo | não aberta; **pode abrir já** | catálogo da §3 |
| E2.4 piloto | não aberta; **pode abrir já** (D27 fechou o que faltava) | catálogo da §3 |
| E1.7, E2.5, E3, E4, E5b, E6, E7 | não abertas; E1.7 espera a decisão de L2f (`selecao-estrutura.md`) | |
| L1 verificação bibliográfica | **fechada** (2026-09-18): 35 entradas verificadas | `referencias-verificadas.bib`, `literatura.md` |
| L2 busca de novidade | **fechada** (2026-09-18): novidade confirmada, Klopp & Pensky (2015) é o vizinho | `busca-novidade.md`, `literatura.md` |

Esta tabela é atualizada pelo chat principal quando uma etapa fecha; em caso
de dúvida, o `ESTADO.md` manda.

## 3. Catálogo das tarefas abertas

Cada linha diz o que entregar, de que depende, e **os únicos arquivos que a
tarefa pode criar ou editar**.

| Tarefa | Entregável | Depende de | Arquivos permitidos |
|---|---|---|---|

| E1.8 | `derivations/07-rota-intervalo.tex` (+ `.pdf`) e `check/07-rota-intervalo.R`: deixar **pronta e escrita** a teoria sem periodicidade, para o dia em que for preciso, em três partes. **(A) O que transfere para a base do intervalo (CDV)**, com o que é verbatim e o que não é: a aproximação, porque a CDV caracteriza `B^s_{π,r}[0,1]` sem periodicidade (mesma prova do Lema 2, sem operador de extensão e sem `C(eps)`); a identificabilidade, que já está feita na §5 de `01-identificabilidade.md` com a reparametrização `[Φ Q \| Ψ]`, cujas colunas são ortonormais em `L_2[0,1]` e de integral zero; o desenho, porque a prova de E1.4 usa só ortonormalidade em `[0,1]` e as cotas da densidade; e o oráculo e as taxas, que consomem só essas peças. **Provar ou conferir o que não transfere de graça:** a cota pontual `Σ_{jk} ψ²_{jk}(u) ≤ C_ψ 2^J` na base CDV (E1.4 conferiu só na periódica) e o efeito dos `p q (2^{j_0} − 1)` coeficientes de escala não penalizados no termo `σ² · (não penalizados)/n` do Teorema 1. **(B) O caso `eps > 0`**, consolidando o Lema 10 e D26: o cruzamento exato das duas exigências (`eps ≳ (L−1)2^{−J}` para excluir a faixa contaminada, `2 eps < (L−1)2^{−J}` para não esconder wavelet), por que o ganho assintótico é exatamente zero, e **qual enunciado é verdadeiro**: um resultado de amostra finita, para `eps` fixo e `J` na faixa admissível, e não um teorema de taxa. **(C) A análise de `s'`**, que é o que o autor pediu por escrito: para cada componente dos cenários de `dgp.R` (seno, cosseno, cúbica, bumps, blocks, heavisine), a regularidade efetiva sob os três regimes — periódica com `eps = 0`, periódica com margem e extensão, e intervalo —, explicando por que a cúbica lê `3/2` no primeiro e `4` nos outros dois (quina da extensão periódica contra os `N = 4` momentos nulos) e por que `blocks` e `heavisine` leem `1/2` nos três (o salto é interior). A conferência mede a queda por nível de cada componente nos três regimes e bate com a tabela. Numeração global: continuar de Proposição 4, Lema 10, Teorema 2, Corolário 5 | E1.2 a E1.6 e E1.3b, fechadas; D26 e D27 | `derivations/07-rota-intervalo.tex`, `derivations/07-rota-intervalo.pdf`, `derivations/check/07-rota-intervalo.R`, `docs/handoff-E1.8.md` |
| E2.4 | `wafc/R/competitors.R` e `wafc/scripts/04-pilot.R`, mais `wafc/tests/test-competitors.R`: o piloto de `plano-projeto.md` E2.4, já com o que as rodadas de 2026-09-19 acrescentaram. **Concorrentes** (L2d): `mgcv::gam` com `s(u_m, by = x_ℓ)`; B-splines mais group LASSO (`grpreg`); o spline adaptativo de Wang, Jiang & Liu (2024); o **block LASSO de Klopp & Pensky** no mesmo desenho (`grpreg`/`gglasso`), que é a pergunta que D18 convida; o VCBART; e a regressão linear oráculo. **Cenários** de `dgp.R`, com `s'` declarado por D27 (`3/2` no suave, `1/2` no não homogêneo) e o número gravado como atributo do cenário. **Regras de `λ`**: as cinco de E2.3 mais o **QUT** de Giacobino et al. (2017) (L2c), e o oráculo da grade como coluna permanente, que é o que dá sentido a "custo". **Margem**: medir ISE e `λ_min(G_eps)` em `eps ∈ {0, 2^{−(J+1)}, 1.9^{−J}, 0.02, 0.05, 0.10}` e na família `a(L−1)2^{−J}` com `a` em torno de 1; o resultado fixa `wafc_eps_periodic`, hoje provisório em `0.05`. **Consertar de passagem** as três linhas de `wafc/R/tune.R` que descrevem o padrão antigo de `eps` (Roxygen de `cv.wafc(J = )`, a guarda da regra da teoria que para em erro quando `J_n = 1`, e o comentário de `wafc_J_grid()`) e decidir, com número, se a grade passa a incluir `J = 1`. 50 réplicas; `n ∈ {250, 500, 1000}` | E2.1b, E2.2 e E2.3, fechadas; D19, D20, D23, D26, D27 | `wafc/R/competitors.R`, `wafc/scripts/04-pilot.R`, `wafc/tests/test-competitors.R`, `wafc/R/tune.R` (só as três linhas), `wafc/README.md` (só a tabela), `docs/handoff-E2.4.md` |
Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**Abertas: E1.8 e E2.4.** Não compartilham arquivo: E1.8 está em
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

**Este mapa é a autoridade**, e é ele que E5a usa ao montar o manuscrito. O
`04-oraculo.tex` já imprime o número global (via `\setcounter` no
preâmbulo), prática adotada daqui em diante; `02` e `03` ainda imprimem o
contador local. O próximo resultado numera a partir de Proposição 4, Lema 10,
Teorema 2 e Corolário 5.

**O enunciado que vai ao resumo do artigo é o Corolário 5** (D16, a
ratificar).

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
