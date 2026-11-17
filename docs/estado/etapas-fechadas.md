# TAREFA antigo: §2 O que já foi feito
Movido de docs/TAREFA.md em 2026-11-16, sem alteração de linhas.
Arquivo morto: não ler inteiro; consultar por grep pelo identificador (D56, E6.1b, pergunta 33) e ler o trecho.

<!-- conteudo movido abaixo -->
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
| E4.3b a grade nova do `gam.reml` | **fechada** (2026-10-08): a regra de tamanho acima de 80, o limite de 80 iterações, o descarte; `k` até 240 no `gam.reml` em `n = 1000` e `2000`; junção exata onde tudo converge; a E4.4 em ~74 h | `wafc/R/competitors.R`, `wafc/tests/test-competitors.R`, `wafc-studies/config/study.yaml` |
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

