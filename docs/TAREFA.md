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
| E3.3, E4, E5b, E6.2, E7 | não abertas | |
| L1 verificação bibliográfica | **fechada** (2026-09-18): 35 entradas verificadas | `referencias-verificadas.bib`, `literatura.md` |
| L2 busca de novidade | **fechada** (2026-09-18): novidade confirmada, Klopp & Pensky (2015) é o vizinho | `busca-novidade.md`, `literatura.md` |

Esta tabela é atualizada pelo chat principal quando uma etapa fecha; em caso
de dúvida, o `ESTADO.md` manda.

## 3. Catálogo das tarefas abertas

Cada linha diz o que entregar, de que depende, e **os únicos arquivos que a
tarefa pode criar ou editar**.

| Tarefa | Entregável | Depende de | Arquivos permitidos |
|---|---|---|---|
| **L12** as fontes da aplicação | verificar, como em L1 a L11 (D39 nos nomes, D25 para software e dados), as referências que a aplicação de D57 e D59 pede, hoje `[VERIFICAR]` no `docs/literatura.md`: Carslaw (2005, *Atmos. Environ.* 39, 4793–4802, a fração primária de NO2 em Londres), **lido no PDF** para conferir o degrau de 2003 (posição e tamanho, de ~5–6% a ~17%) e se o texto explica o mergulho do começo de 2002 e o pico do fim de 2003 que a E6.1b e a E6.1c viram; Carslaw & Beevers (2005); Clapp & Jenkin (2001, o oxidante `NO2 + O3`, com o segundo autor conferido); Hersbach et al. (2020, *QJRMS*, o ERA5) e o conjunto ERA5 do C3S (a forma de citação do DataCite); o Open-Meteo (Zippenfenig 2024); o UK-AIR (Defra) e o relatório de conversão de unidades; e, se o artigo usar o erro-padrão corrigido, Nadeau & Bengio (2003). Os PDFs que não abrirem por `curl` são pedidos ao autor, com a referência completa, para ele pôr em `refs/` no padrão sobrenome e ano. Entradas no `docs/referencias-verificadas.bib`; as linhas do `literatura.md` passam de `[VERIFICAR]` a verificadas. No handoff: as entradas, as localizações, o que Carslaw (2005) diz do degrau e dos ressaltos, e as frases de atribuição que as licenças pedem (UK-AIR, ERA5, Open-Meteo) | D57, D59 | `docs/referencias-verificadas.bib`, `docs/literatura.md`, `refs/`, `docs/handoff-L12.md` |

Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**L12 catalogada em 2026-10-06** (as fontes da aplicação, D59). E6.1c
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

**O manuscrito está vivo em `k = 4`** (aberto em 2026-10-04 pela E5d com
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
