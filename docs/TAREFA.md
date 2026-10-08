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
| E0.1 nomes | fechada: D4 (código em `wafc/`), D5 e D8 ratificadas em 2026-09-19 | `plano-projeto.md` E0.1, `ESTADO.md` |
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
| Revisão de `wafc/` | **feita** (2026-09-28): cinco defeitos, `call` compacto, só o melhor ajuste no `cv.wafc`, busca de nós rápida; 672 testes passam, piloto idêntico | `wafc/R/`, `wafc/scripts/04-pilot.R`, `wafc/tests/` |
| P1, P2 e perguntas 29 a 31 | **decididas** (2026-09-28): grade `2:8` (D34), margem `0` (D35), `bsgl` na mesma grade (D36), nenhuma aceleração (D37); 674 testes passam | `wafc/R/tune.R`, `wafc/R/design.R`, `wafc/R/competitors.R`, `wafc/scripts/07-accel.R`, `wafc/scripts/08-sgl-null.R` |
| E6.1a sondagem de bases | **fechada** (2026-09-20): veredito negativo nas três candidatas | `docs/aplicacao-candidatas.md` |
| E2.4c preparar o piloto | **fechada** (2026-09-30): célula `uneven`, `gam.matched` por `bam`, `WAFC_TAG`, sementes independentes da restrição; 674 testes, fumaça sem falha | `wafc/scripts/04-pilot.R` |
| E2.5a repetir o piloto | **fechada** (2026-09-30): 6 450 linhas, 0 falhas; no-go para a variante LASSO pelo critério literal (pergunta 33 do `ESTADO.md`) | `wafc/cache/e25a/` (não versionado), `ESTADO.md` §2 |
| E5c manuscrito `k = 2` | **fechada** (2026-09-30): pergunta 15 aplicada e marcada, Corolário 8 no artigo (§3.6, S7); 33 e 30 páginas, `.bib` com 37 entradas | `manuscript/ms_2.tex`, `supp_2.tex`, `references_2.bib` |
| E2.5b medir o block LASSO com níveis livres | **fechada** (2026-09-30): o `klopp.free` vence o `wafc.lasso` em toda célula com componente e empata no nulo; a `mixed` em 50 réplicas; `WAFC_METHODS` e `WAFC_REPS_MIXED`; 678 testes | `wafc/cache/e25b/` (não versionado), `ESTADO.md` §2 |
| E1.11 sondagem da teoria com blocos | **fechada** (2026-09-30): transfere sem cone, ganha `(log n)^{2s'/(2s+1)}`, pesos do `grpreg` fora da teoria; conferência `OK` | `derivations/08a-sondagem-blocos.md` |
| E2.5c os pesos do block LASSO | **fechada** (2026-09-30): os pesos 1 da teoria predizem pior; `grpreg` melhor no suave, níveis grossos juntos no não homogêneo e na `mixed`; 696 testes | `wafc/R/competitors.R`, `wafc/cache/e25c/` |
| E2.5d o `gam` autônomo | **fechada** (2026-09-30): empata com o `gam.matched`, confirma E2.5a, fator do suave pior; reprodução exata | `wafc/scripts/09-gam-autonomo.R`, `wafc/cache/e25d/` |
| L5 bibliografia curta | **fechada** (2026-09-30): 69 entradas; a partição da unidade ancorada em Mallat (2009); D39 aplicada | `referencias-verificadas.bib`, `literatura.md` |
| E2.5e pedaços balanceados | **fechada** (2026-10-01): coberta pela teoria com os pesos do `grpreg` (`ρ² ≤ 1,67`), melhor das quatro formas fora do `smooth`; 730 testes | `wafc/R/competitors.R`, `derivations/08a-sondagem-blocos.md` §11, `wafc/cache/e25e/` |
| E2.5f níveis grossos livres | **fechada** (2026-10-01): não paga; a pior das cinco formas; teoria cobre com `σ²p_0/n`; 807 testes | `wafc/R/competitors.R`, `derivations/08a-sondagem-blocos.md` §12, `wafc/cache/e25f/` |
| E2.5g estimação seguida de limiar | **fechada** (2026-10-01): `wafc_threshold()`; o limiar passa o fator do suave (1,31 a 1,37 no `wafc.lasso+max`); `c = 0,4` acerta 0,88 a 1 em `n ≥ 500`; 921 testes | `wafc/R/threshold.R`, `wafc/scripts/04-pilot.R`, `wafc/cache/e25g/` |
| L6 correções bibliográficas nas derivações | **fechada** (2026-10-01): Mallat, Hsu et al., Cai e van de Geer et al. onde as provas os usam; `06-selecao-limiar.pdf` recompilado | `derivations/01-identificabilidade.md`, `08a-sondagem-blocos.md`, `06-selecao-limiar.tex` |
| L7 numeração conferida nas fontes publicadas | **fechada** (2026-10-01): K&P, Lounici et al., Huang et al., Donoho & Johnstone e Bühlmann & van de Geer conferidos; o Teorema 2 de K&P mudou do arXiv ao *Annals* | `derivations/08a-sondagem-blocos.md`, `06-selecao-limiar.tex` |
| E2.5h o `gam` por REML, GCV e validação cruzada, e o GCV do WAFC | **fechada** (2026-10-02): o veredito se mantém com o `gam` sintonizado; `gam.reml` = `gam.cv`; o GCV do WAFC não ganha; 1 006 testes | `wafc/R/competitors.R`, `wafc/R/tune.R`, `wafc/scripts/04-pilot.R`, `wafc/cache/e25h/` |
| L8 referências do `mgcv` e da escolha da dimensão | **fechada** (2026-10-01): 82 entradas; Kauermann & Opsomer usam ML, e o catálogo da E2.5h foi corrigido | `referencias-verificadas.bib`, `literatura.md` |
| L9 as duas frases que L7 deixou | **fechada** (2026-10-01): a entrada de 1994 saiu; o risco minimax de 1998 é sobre corpos de Besov (pergunta 35(h)) | `derivations/05-taxas.tex`, `docs/busca-novidade.md` |
| E2.5j a regra de `t` e a porta do nulo | **fechada** (2026-10-03): o `cv1se` chega ao teto no suave, no `uneven`, no nulo e em `n = 1000`; o `cvrel` repete o `+cv`; a porta do QUT com nível acima do nominal; os cortes dos motores não tocam o `λ` escolhido; 1 133 testes | `wafc/R/threshold.R`, `wafc/R/tune.R`, `wafc/R/fit.R`, `wafc/R/competitors.R`, `wafc/scripts/04-pilot.R`, `wafc/cache/e25j/` |
| E2.5 go/no-go e variante principal | **fechada** (2026-10-03, D44 a D46): go reposicionado na *Statistica Sinica*; o WAFC é o block LASSO balanceado com limiar `cv1se`, o LASSO fica como opção; `gam.reml` e `gam.gcv` em E4; E2 inteira fechada | `ESTADO.md` §2 e decisões |
| E1.12 a teoria em blocos | **fechada** (2026-10-03): Proposição 7, Lema 14, Teorema 3, Corolário 9, Lema 15, Teorema 4, Corolários 10 e 11 no `08`; Lema 16 e Corolários 12 e 13 no `06`; conferência `OK` | `derivations/08-blocos.tex`, `06-selecao-limiar.tex`, `check/08-blocos.R` |
| E3.1 a interface congelada | **fechada** (2026-10-03): `penalty = "block"` e `threshold = "cv1se"` como padrão do `cv.wafc()`; oitava junção exata; 1 075 testes em 46 s (1 300 com `WAFC_SLOW_TESTS=1`); interface ratificada (D48) | `wafc/R/`, `wafc/tests/`, `wafc/README.md`, `wafc/scripts/02` a `08`, `wafc/cache/e31/` |
| E1.13 a taxa lenta em blocos | **fechada** (2026-10-03): Corolário 14 no `08`, com um logaritmo a menos que a Proposição 4 e, em `π ≥ 2`, a taxa do Teorema 4 sem condição de desenho; conferência `OK` | `derivations/08-blocos.tex`, `check/08-blocos.R` |
| E3.2 gráficos e documentação | **fechada** (2026-10-03): `plot.wafc`, `plot.cv.wafc`, roxygen auditado, o `coef` limiarizado corrigido, exemplo no README; 1 202 testes em 48 s | `wafc/R/plot.R`, `wafc/tests/test-plot.R`, `wafc/README.md`, `wafc/man-figures/` |
| E1.14 a taxa lenta em todo `s' > 0` | **fechada** (2026-10-04): itens (iii) e (iv) da Proposição 4 e do Corolário 14, com o `λ` inflado acima do teto; conferências `OK` | `derivations/05-taxas.tex`, `08-blocos.tex`, `check/05-taxas.R`, `check/08-blocos.R` |
| E5d o manuscrito em `k = 4` | **fechada** (2026-10-04): `ms_4` e `supp_4` com a teoria em blocos e a S8 do lasso, 50 e 60 páginas (43 e 59 sem o removido); sem os itens de E1.14 | `manuscript/ms_4.tex`, `supp_4.tex`, `references_4.bib` |
| E1.15 o comparador truncado na taxa lenta | **fechada** (2026-10-05): Proposição 8 no `05` e Corolário 15 no `08`; em `s < 1/2`, a taxa do Theorem 1 sem condição de desenho; conferências `OK` | `derivations/05-taxas.tex`, `08-blocos.tex`, os checks |
| E5e a E1.14 na `k = 4` | **fechada** (2026-10-05): Proposition S6.1, os itens da S8.1, "program", "thresholding"; 50 e 68 páginas | `manuscript/ms_4.tex`, `supp_4.tex` |
| L10 a bibliografia do `grpreg` e do group lasso | **fechada** (2026-10-05): Yuan & Lin (2006), Breheny & Huang (2015), `grpreg` 3.6.0; 86 entradas | `docs/referencias-verificadas.bib`, `docs/literatura.md` |
| E5f a E1.15 e a L10 na `k = 4` | **fechada** (2026-10-05): Propositions S6.2 e S8.2, a frase do §3 qualificada, as citações do `grpreg` e do group lasso; 52 e 80 páginas | `manuscript/ms_4.tex`, `supp_4.tex`, `references_4.bib` |
| L11 Simon & Tibshirani e Breheny & Huang no periódico | **fechada** (2026-10-05): `Simon-Tibshirani-2012` verificado, 87 entradas; a paginação de Breheny & Huang inferida | `docs/referencias-verificadas.bib`, `docs/literatura.md` |
| E1.9 a cota inferior no nível da taxa | **fechada** (2026-10-05): Lema 17, Teorema 5, Corolário 16 em `09-cota-inferior.tex`; constante proporcional a `pq`; conferência `OK` | `derivations/09-cota-inferior.tex`, `check/09-cota-inferior.R` |
| E6.1b a aplicação depois de D44 | **fechada** (2026-10-05): seis bases, 20 partições por bloco; nenhuma passa contra o `gam.cv` em blocos; marylebone é a melhor estrutura | `wafc/scripts/10-sondagem-aplicacao-b.R`, `docs/aplicacao-candidatas.md` §9 a §16, `wafc/cache/e61b/` |
| E3.4 `J` por moduladora | **fechada** (2026-10-05): `cap.J = TRUE`, a contagem no círculo; nona junção exata; 1 304 testes em 55 s | `wafc/R/design.R`, `fit.R`, `tune.R`, `threshold.R`, `plot.R`, `wafc/tests/test-cap.R` |
| E5g a E1.9 na `k = 4` | **fechada** (2026-10-05): Theorem S6.1 e Lemma S6.3 no supp, a otimalidade na introdução e na terceira leitura; 52 e 90 páginas | `manuscript/ms_4.tex`, `supp_4.tex`, `references_4.bib` |
| E6.1c marylebone com dados de licença declarada | **fechada** (2026-10-06): UK-AIR e ERA5; o degrau de 2003 em 20 de 20 partições; o `+cv` empata com os splines de D46; marylebone fica como principal | `wafc/scripts/10-sondagem-aplicacao-b.R`, `docs/aplicacao-candidatas.md` §17, `wafc/cache/e61c/` |
| E4.2 desenho do estudo | **decidido** (2026-10-06, D60): 100 réplicas, três braços, a escala sob condição de custo, a regra do topo das grades em `n = 2000`; o compêndio nasce como a pasta `wafc-studies/` | `plano-projeto.md` E4.1 e E4.2 |
| L13 as duas fontes restantes da aplicação | **fechada** (2026-10-06): 99 entradas; Liang et al. (2015) e Opsomer, Wang & Yang (2001) lidos no PDF; as datas da temporada são nominais | `docs/referencias-verificadas.bib`, `docs/literatura.md` |
| L12 as fontes da aplicação | **fechada** (2026-10-06): 97 entradas; cinco artigos lidos no PDF; o degrau de Marylebone Road é de ~10 a ~23 vol% (Carslaw 2005, Fig. 3(a)), não a média de Londres | `docs/referencias-verificadas.bib`, `docs/literatura.md`, `refs/` |
| E6.2b o topo das grades da aplicação | **fechada** (2026-10-07): a grade fica pelo critério declarado; as leituras sem os buracos | `wafc/scripts/12-app-grid-probe.R`, `wafc-studies/R/application_report.R` |
| E5j marylebone na `k = 5` | **fechada** (2026-10-07): a figura com buracos interrompidos e rug; os números de marylebone no §6.3; D72; 73 e 97 páginas | `wafc-studies/R/application_report.R`, `results/`, `manuscript/` |
| E5i as decisões de texto na `k = 5` | **fechada** (2026-10-06): D68(a), (d), (i), D69 e D71 aplicadas; 72 e 97 páginas; 64 entradas | `manuscript/ms_5.tex`, `supp_5.tex`, `references_5.bib` |
| L14 a citação que a UCI pede | **fechada** (2026-10-06): 101 entradas; Zhang et al. (2017) e `Chen-UCI-2017`; a troca de estação do Dongsi em março de 2015 | `docs/referencias-verificadas.bib`, `docs/literatura.md`, `refs/` |
| E5k a aplicação fechada na `k = 5` | **fechada** (2026-10-08): os marcadores `[E6.2]` todos trocados; Tabelas S5 a S7; Figura S1; 73 e 101 páginas | `manuscript/` |
| E4.3 o piloto do estudo | **fechada** (2026-10-07): 884 unidades, 0 falhas, 2 h 26 min; a E4.4 projeta ~59 h em 8 trabalhadores; a escala entra; os splines no topo (pergunta 67) | `wafc/cache/e43/` |
| E4.3a preparar o piloto | **fechada** (2026-10-06): `wafc/scripts/11-pilot.yaml` e `11-pilot-read.R`, mestre `43000000`; fumaça de 38 unidades sem falha; a rodada (~3 h em 8, sozinha) espera o fim da E6.2 e o aviso | `wafc/scripts/11-*`, `wafc/cache/e43/` |
| E5h as Seções 5 e 6 sem números | **fechada** (2026-10-06): a `k = 5` aberta; Seções 5 e 6, S9 e S10, com marcadores `[E4.3]`, `[E4.4]`, `[E6.2]`; 72 e 97 páginas (64 e 96 sem o removido) | `manuscript/ms_5.tex`, `supp_5.tex`, `references_5.bib`, `tables/`, `figures/` |
| E6.2 a aplicação no compêndio | **rodada** (terminada em 2026-10-07 às 16h31, 168 unidades, 0 falhas; antes, código pronto em 2026-10-06): as duas bases em `wafc-studies/`, dados versionados idênticos ao script 10, junção exata da partição 1 fora o WAFC da `beijing.heat`; **as rodadas (~70 h de processador) esperam o aviso do autor** (pergunta 59) | `wafc-studies/R/application*.R`, `scripts/03_application.R`, `config/application.yaml`, `data-raw/`, `data/` |
| E4.1c o oráculo em blocos e o `structure` | **fechada** (2026-10-06): `penalty = "block"` no `wafc_fit_oracle()`, `structure` no `simulate_wafc()` sem sorteio novo; 1 730 testes; décima primeira junção exata | `wafc/R/competitors.R`, `wafc/R/dgp.R`, os testes, `wafc-studies/` |
| E4.1 o nascimento do compêndio | **fechada** (2026-10-06): a pasta `wafc-studies/` com `renv`, configuração por `n`, cache retomável e sementes independentes do piloto; fumaça sem falha; décima junção exata | `wafc-studies/`, `wafc/cache/e41/` (não versionado) |
| E4.1b `X` dependente de `U` | **fechada** (2026-10-06): `x_u_rho` no `simulate_wafc()`, `h` linear e limitada, o padrão reproduz o `HEAD` em 900 de 900; 1 422 testes em 54 s | `wafc/R/dgp.R`, `wafc/tests/test-dgp.R` |
| E3.3 empacotamento | **decidida** (2026-10-06, D61): `wafc/` fica até o fim e vai então ao `WaveBased` | `plano-projeto.md` E3.3 |
| E4.4, E5b, E7 | não abertas | |
| L1 verificação bibliográfica | **fechada** (2026-09-18): 35 entradas verificadas | `referencias-verificadas.bib`, `literatura.md` |
| L2 busca de novidade | **fechada** (2026-09-18): novidade confirmada, Klopp & Pensky (2015) é o vizinho | `busca-novidade.md`, `literatura.md` |

Esta tabela é atualizada pelo chat principal quando uma etapa fecha; em caso
de dúvida, o `ESTADO.md` manda.

## 3. Catálogo das tarefas abertas

Cada linha diz o que entregar, de que depende, e **os únicos arquivos que a
tarefa pode criar ou editar**.

| Tarefa | Entregável | Depende de | Arquivos permitidos |
|---|---|---|---|
| **E4.3b** a grade nova do `gam.reml` | (o) **Emenda de D80 (2026-10-08, depois da primeira rodada com a regra; o handoff dela fica em `docs/handoff-E4.3b.md` e é reescrito no fim):** (1) a regra de tamanho vale **só para os candidatos com `k` acima de 80**, a grade comum de D41, e a junção exata passa a cobrir também a `mixed` em `n = 250` e `500`; (2) **um limite determinístico de iterações** do otimizador do `mgcv` (o argumento de controle das iterações externas do REML, ou o equivalente no `bam`, com o valor escolhido e justificado no handoff) em cada candidato da busca, e o não convergir (aviso ou indicador de convergência do `mgcv`) conta como falha do candidato, que sai da busca com a razão em `error`; conferir em minutos nas réplicas 1 e 3 do `urho` em `n = 1000` (a 3 é a que passou de 2 h 50) que terminam e o que escolhem, e que nas células onde tudo converge os ajustes não mudam (junção); refazer a projeção da E4.4. **Se o limite não resolver o `urho`, parar e voltar como pergunta** (aceitar e medir na produção está descartado). (i) **a regra de tamanho de D79 (emenda de 2026-10-08):** antes de ajustar, um candidato `k` com `(número de suavizações) × k > 2n` sai da busca, com a razão na tabela `gam_k` (a regra vale para os dois splines; no `gam.gcv` a grade atual não a aciona, conferir); **e o descarte do candidato inviável** em `wafc/R/competitors.R`, que fica como rede de segurança: na busca de `k` do `gam` (`wafc_fit_gam()` e o `wafc_gam_reml()`), um candidato cujo ajuste ou cujo escore falha (o `stop` "the penalty of a smooth has rank below the one mgcv declares", ou mais coeficientes que observações) sai da busca com a razão registrada na tabela `gam_k` do `extra`, em vez de derrubar o ajuste; se nenhum candidato servir, o erro de hoje. Teste em `wafc/tests/test-competitors.R` (um caso pequeno em que o maior `k` é inviável e o ajuste escolhe entre os demais) e a **junção exata**: onde todos os candidatos são viáveis, os ajustes de hoje não mudam (as unidades `gam.reml` do piloto em `wafc/cache/e43/` refeitas idênticas fora o tempo, numa célula barata e na `mixed` em `n = 1000`); a suíte inteira passa. (ii) **A grade** no `wafc-studies/config/study.yaml`: nas entradas `"1000"` e `"2000"` do `tuning`, `methods: {gam.reml: {k: [5, 10, 20, 40, 80, 120, 160, 240]}}`, deixando o `gam.gcv` na grade de hoje (`80` em `n = 1000`, `120` em `n = 2000`); o comentário no estilo do estudo final (sem D-números). (iii) **Conferência na `mixed` em `n = 1000`** (refeita com a regra de D79; a primeira tentativa, sem ela, travou com `k = 160` e `240` e foi parada pelo chat principal; as unidades de `wafc/cache/e43b/probe/` e a de `inhomogeneous.urho` rep 3 em `timing/` não existem ou estão incompletas): as réplicas 1 a 5 da sondagem (`wafc/scripts/13-gam-grid-probe.yaml`, mestre do piloto) com o código novo: o `k` escolhido, o maior viável, o `rmse_f` contra o piloto (a razão, para completar a tabela de D78), o tempo e o pico de memória. (iv) **A reprojeção da E4.4** com os tempos novos do `gam.reml` (o leitor `wafc/scripts/11-pilot-read.R` ou uma conta à parte), em horas de processador e de relógio em 8. Rodadas de até ~2 h, sozinha ou ao lado de nada pesado. O `PROVENANCE.md` muda com o commit (o chat principal atualiza). No handoff: o diff de comportamento, os números das conferências e a projeção | D41, D69, D77, D78 | `wafc/R/competitors.R`, `wafc/tests/test-competitors.R`, `wafc/README.md` (só a linha do `competitors.R`, se mudar), `wafc-studies/config/study.yaml`, `wafc/cache/e43b/` (não versionado), `docs/handoff-E4.3b.md` |

Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**E4.3b emendada em 2026-10-08** (D79, a regra de tamanho; e D80, a regra só além de `k = 80` e o limite de iterações; o chat da tarefa relê a entrada). **E5k fechou em 2026-10-08** (pergunta 68). **E4.3b catalogada em 2026-10-08** (D78); a E4.4 espera por ela. **E5k catalogada em 2026-10-07** (a Figura S1 refeita depois do piloto já está em `results/`). **E6.2b fechou em 2026-10-07** (pergunta 66); a aplicação terminou e o piloto roda desde 16h31. **E6.2b catalogada em 2026-10-07** (D75): a sondagem tem de terminar antes do fim da `beijing.heat`, senão atrasa o piloto. **E5j fechou em 2026-10-07** (pergunta 65). **E5j catalogada em 2026-10-07** (D73; sem ajuste, só o relatório de marylebone, enquanto a cadeia roda a `beijing.heat`). **E5i fechou em 2026-10-06** (pergunta 63). **L14 fechou em 2026-10-06** (pergunta 62). **E5h fechou em 2026-10-06** (pergunta 60; o manuscrito vivo é a `k = 5`). **E4.3a e E5h catalogadas em 2026-10-06**, sem arquivo em comum e sem rodada: a E4.3a só prepara o piloto (a rodada espera o fim da E6.2 e o aviso, sozinha na máquina, D67), a E5h só escreve na `k = 5`. **E6.2 com o código pronto em 2026-10-06**; as rodadas esperam o aviso (pergunta 59). **E4.1c fechou em 2026-10-06**; a E4.3 pode ser catalogada. **L13 fechou em 2026-10-06** (pergunta 58). **L13 catalogada em 2026-10-06**, sem arquivo em comum com a E4.1c e a E6.2. **E6.2 catalogada em 2026-10-06**, em paralelo à E4.1c, sem arquivo em comum: a E6.2 só cria arquivos novos no compêndio e deixa no handoff o que muda no `README.md`, no `run_all.R` e no `renv.lock` dele; se precisar mexer nos `R/` comuns (`cli.R`, `run.R`), também vai ao handoff. **L12 fechou em 2026-10-06** (pergunta 57). **E4.1c catalogada em 2026-10-06** (D63), sem arquivo em comum com a L12; a E4.3 espera por ela. **E4.1 fechou em 2026-10-06** (pendências na pergunta 56; a próxima junção exata inclui o nulo em `n = 250`). **E4.1b fechou em 2026-10-06** (o braço passa `x_u_rho = 0.5` também à amostra de teste). **E4.1 e E4.1b catalogadas em 2026-10-06** (D60), sem arquivo em comum: a E4.1 só cria dentro de `wafc-studies/` e lê `wafc/`; a E4.1b só toca o `dgp.R`, o teste dele e talvez o `wafc/README.md`. A E4.1 pode escrever a configuração do braço com `x_u_rho` antes, e a fumaça desse braço espera a E4.1b; a junção exata da E4.1 não depende dela, porque o padrão `x_u_rho = 0` reproduz os sorteios. Nenhuma das duas faz rodada longa.

**L12 catalogada em 2026-10-06** (as fontes da aplicação, D59). Os PDFs dos artigos entraram em `refs/` no mesmo dia. E6.1c
fechou em 2026-10-06. **E3.4 e E5g fecharam em 2026-10-05.**
Catalogadas em 2026-10-05 (D57), sem arquivo em
comum: a E6.1c só toca o script 10 (a base nova), os documentos da aplicação
e os dados, e lê o código do retrato `b0ea096`; a E3.4 só `wafc/R/`,
`wafc/tests/` e o `wafc/README.md`; a E5g só a `k = 4`. A E6.1c e a E3.4
disputam núcleos só na rodada da E6.1c (até 3 processos) e na junção exata
da E3.4 (~10 min).

**E5f e L11 fecharam em 2026-10-05.**

**E5e, E1.15 e L10 fecharam em 2026-10-05.** A E6.1b segue rodando.

**E1.14 e E5d fecharam em 2026-10-04.** A E6.1b está rodando (4
processos, 20 partições); pendências de E1.14 e E5d nas perguntas 46 e 47
do `ESTADO.md`. **O manuscrito está vivo em `k = 4`.**

**E1.13 e E3.2 fecharam em 2026-10-03.** A E6.1b segue no catálogo,
pronta e à espera da ordem do autor para a rodada de 4 processos.

**E1.12 e E3.1 fecharam em 2026-10-03.** **A E6.1b segue no catálogo**; a
E6.1b parou preparada antes da rodada, à espera do aviso e da permissão de
download do autor (`ESTADO.md` §5). Catalogadas em 2026-10-03 para correr em
paralelo: a E1.12 só toca `derivations/`, a E3.1 só `wafc/R/`, `wafc/tests/`
e os scripts 01 a 09, e a E6.1b só o script 10 e os documentos da
aplicação, lendo o código de um retrato do `HEAD`. A E1.12 não usa núcleos;
a E3.1 usa até 6 por ~15 min na junção exata, e a E6.1b até 6, com rodada
acima de 1 h esperando o aviso do autor.

**Catálogo vazio até 2026-10-03**: E2.5j fechou. A tabela mais recente é
`wafc/cache/e25j/e25j-joined.rds` (37 350 linhas, sétima junção exata), e é
a ela que um método novo se junta. Desde D43 nenhuma variante sai do
código. E2.5 fechou com D44 a D46; o que vem a seguir está no
`ESTADO.md` §5 (a teoria em blocos nas derivações numeradas, a catalogar).

**Catálogo vazio até 2026-10-02**: E2.5h (com a E2.5i dentro), L8 e L9 fecharam; E2.5e, E2.5f, E2.5g, L6 e L7 fecharam. A repetição do piloto (E2.4c,
E2.5a) devolveu no-go para a variante LASSO, e as medições que a pergunta
33 do `ESTADO.md` pedia fecharam no mesmo dia: E2.5b (o block LASSO com
níveis livres), E2.5c (os pesos), E2.5d (o `gam` autônomo), E1.11 (a
teoria com blocos) e L5 (bibliografia). E2.5 é do chat principal e espera
a decisão de rumo do autor, que desde E2.5g tem todas as medições (perguntas 33 e 38). Para acrescentar um método às tabelas, o
procedimento validado seis vezes é rodá-lo com `WAFC_METHODS` e conferir um
método já existente numa célula barata antes de juntar ao
`wafc/cache/e25h/e25h-joined.rds` (28 350 linhas, a tabela mais recente; junção exata seis vezes).

**Ao catalogar, o arquivo vai junto do item.** Duas tarefas seguidas
esbarraram em coluna de arquivos que não cobria o que o próprio texto
mandava fazer: E2.4 pôs o QUT em `competitors.R` por isso, e E2.4b teve de
tocar `fit.R` e criar `06-timing.R` fora da coluna.

A cota inferior para `q ≥ 2`, que o handoff de E1.6 chegou a chamar de E1.8,
é a **E1.9**, catalogada em 2026-10-05 (a numeração E1.8 ficou com a rota do
intervalo).

**O manuscrito está vivo em `k = 5`** (aberto em 2026-10-06 pela E5h, com as Seções 5 e 6 sem números; a `k = 4` foi aberta em 2026-10-04 pela E5d com
a teoria em blocos; `k = 1` a `k = 3` ficam intactas). Toda alteração é marcada em `colR1`; `ms`, `supp` e
`references` andam juntos. O mapa da numeração global para os rótulos do
LaTeX está no cabeçalho do `ms_4.tex`.

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
| `06-selecao-limiar.tex` (E1.7c) | Lema 13, Corolário 8 (no manuscrito: Lemma S7.1 e Corollary 2, com a Hipótese S como Assumption 6; mapa no cabeçalho do `ms_4.tex`) |
| `08-blocos.tex` (E1.12) | Proposição 7, Lema 14, Teorema 3, Corolário 9, Lema 15, Teorema 4, Corolário 10, Corolário 11; Hipótese B (letra, como a S); reimprime o Lema 4 em blocos sem número novo |
| `06-selecao-limiar.tex` (adendo de E1.12) | Lema 16, Corolário 12, Corolário 13 |
| `08-blocos.tex` (E1.13) | Corolário 14 (`cor:lenta-blocos`) |
| `05-taxas.tex` (E1.15) | Proposição 8 (`prop:truncada`) |
| `08-blocos.tex` (E1.15) | Corolário 15 (`cor:truncada-blocos`) |
| `09-cota-inferior.tex` (E1.9) | Lema 17, Teorema 5, Corolário 16 (a Proposição 9 ficou livre) |

**Este mapa é a autoridade**, e é ele que E5a usa ao montar o manuscrito. O
`04-oraculo.tex` já imprime o número global (via `\setcounter` no
preâmbulo), prática adotada daqui em diante; `02` e `03` ainda imprimem o
contador local. O próximo resultado novo é **Proposição 9**, **Lema 18**, **Teorema 6**
ou **Corolário 17**, conforme o tipo. O `06-selecao-limiar.tex` imprime `Hipótese S`, com letra
em vez de número, porque é citada lado a lado com a `Hipótese 1` de E1.6;
é desvio local e aceito.

**O enunciado que vai ao resumo do artigo é o Corolário 11** (`08-blocos.tex`;
D47, 2026-10-03, que emenda D16, onde era o Corolário 5).

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
