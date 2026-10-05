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
| **E6.1c** marylebone com dados de licença declarada | refazer a base marylebone da E6.1b com dados de licença limpa, para ela ser a aplicação principal (D57). (i) **Os dados:** a série horária da estação Marylebone Road (MY1) do UK-AIR (Defra, Open Government Licence v3) com NO2, NOx e O3 de 1998 a 2005, e o vento (velocidade e direção) de uma fonte com licença declarada e próxima o bastante (estação meteorológica de Londres ou reanálise); **cada download é pedido ao autor antes**, com o nome do arquivo, a fonte, a licença e o tamanho. Registrar no `docs/aplicacao-candidatas.md` a origem, a licença e a soma SHA-256 de cada arquivo. (ii) **A rodada:** a mesma especificação da E6.1b para marylebone (`Y = NO2 + O3`, `X = (1, NOx/100)`, `U = (data, vento)`, 20 partições por bloco de semana, `J = 2:8`, `+cv1se` e `+cv`, `gam.reml` e `gam.gcv` de D46, o linear), com o mesmo retrato do código (`b0ea096`), para os números serem comparáveis aos da E6.1b; a nova base entra no `wafc/scripts/10-sondagem-aplicacao-b.R` como uma base a mais (por exemplo `marylebone.ukair`), sem mexer nas outras. Até 3 processos; **rodada acima de 1 h espera o aviso do autor**. (iii) **A leitura:** se o degrau de `NOx × data` em 2002–2003 se repete (posição, tamanho, partições); se `NOx × vento` é zerado também pelo `+cv` (D57: só se afirma ausência de efeito onde os dois zeram); a predição contra os dois splines de D46; e a comparação com a E6.1b, base a base. No handoff: o veredito (marylebone fica como principal ou beijing.heat sobe), as licenças e as linhas `[VERIFICAR]` da literatura que a escolha pede | E6.1b, D56, D57 | `wafc/scripts/10-sondagem-aplicacao-b.R` (só a base nova e a sua preparação), `docs/aplicacao-candidatas.md` (seções novas no fim), `docs/literatura.md` (só linhas novas), `docs/handoff-E6.1c.md`; dados e saídas em `wafc/cache/data/` e `wafc/cache/e61c/`, não versionados |
| **E3.4** `J` por moduladora limitado pelos valores distintos | a lição 3 da E6.1b (D57): com uma moduladora discreta de poucos valores (a hora, 24), um bloco com `2^J − 1` colunas acima do número de valores distintos menos um não é identificado, e a norma que o limiar lê passa a depender de direções sem dado. Limitar o nível por moduladora, `J_m = min(J, ⌊log_2(valores distintos de U_m)⌋)`, como o `wafc_k_matched()` faz com o `k` do `gam`, no `wafc_design()` e em tudo o que lê a estrutura do desenho (os pedaços do block LASSO, a validação cruzada sobre a grade de `J`, o limiar, a reconstrução, os gráficos), com o `J` comum mantido como o argumento do usuário e o `J_m` efetivo guardado no objeto e mostrado no `print`. Opção para desligar. **Prova de que nada se moveu onde não devia:** nas moduladoras contínuas o comportamento é idêntico; testes que comparam o ajuste com e sem a opção a `1e-12` num desenho contínuo; e o `04-pilot.R` na célula `smooth` com `WAFC_METHODS=wafc.block` reproduzindo as linhas do `wafc/cache/e31/` em toda coluna fora o tempo (junção exata). Testes novos com uma moduladora de 24 valores: o bloco fica identificado, e a norma não depende de direções sem dado. A suíte padrão segue abaixo de 60 s. No handoff: a interface (a ratificar) e o que muda na aplicação | E3.1, E3.2, D57 | `wafc/R/` (todos), `wafc/tests/` (todos), `wafc/README.md`, `docs/handoff-E3.4.md`; os `.rds` da conferência em `wafc/cache/e34/`, não versionado. **Não tocar** em `wafc/scripts/` |
| **E5g** a E1.9 na `k = 4` | levar à `k = 4` (sem abrir `k = 5`; `colR1`; blocos cinza de D51 para o que sair inteiro; comentário `% E5g` nas inserções em grupos já azuis) o que a pergunta 51 do `ESTADO.md` decidiu (D57): (a) na introdução, "so that with a single modulator independent of the covariates it is minimax optimal there" passa a "where it is minimax optimal", e a frase "A lower bound for several modulators remains open." (D50) passa a "The lower bound that makes it optimal is proved here for any number of modulators and for covariates that depend on them; Klopp and Pensky have it with a single modulator independent of the covariates."; (b) no supp, em S6, o Teorema 5 de `derivations/09-cota-inferior.tex` (com o Lema 17, nas três perdas, e a prova) no lugar da Remark S6.1, que vira uma frase (K&P é o caso `q = 1`, `X ⊥ U`); no corpo, a terceira leitura do Theorem 1 com o texto guardado na pergunta 51 (ajustar o número "Theorem S6.1" ao que o supp der); (c) uma frase na prova do Theorem 1 no supp sobre a uniformidade ("the event does not depend on `f` and the bound depends on `f` only through `C_g`"); (d) `Tsybakov-2009` copiado do `docs/referencias-verificadas.bib`, sem mudar campo, no `references_4.bib`, se o supp citar o lema de Assouad; (e) as macros de `ρ̲_n`, `𝒦`, `K` e `c_A` no preâmbulo do supp, com os nomes do `notacao.md` §10. Os dois compilam sem referência indefinida; páginas com a marcação e sem o removido no handoff | E1.9, D57 | `manuscript/ms_4.tex`, `manuscript/supp_4.tex`, `manuscript/references_4.bib`, os seus PDFs, `docs/handoff-E5g.md` |

Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**E6.1c, E3.4 e E5g catalogadas em 2026-10-05** (D57), sem arquivo em
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
