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
| E2.4 piloto | **em curso** (aberta em 2026-09-19) | catálogo da §3 |
| E1.7a sondagem de irrepresentabilidade | não aberta; **pode abrir já** | catálogo da §3 |
| E1.7c seleção por limiarização | não aberta; **pode abrir já** | catálogo da §3 |
| E1.4c tradução do `03` | não aberta; **pode abrir já**, prioridade baixa | catálogo da §3 |
| E6.1a sondagem de bases | não aberta; **pode abrir já** | catálogo da §3 |
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

| E2.4 | `wafc/R/competitors.R` e `wafc/scripts/04-pilot.R`, mais `wafc/tests/test-competitors.R`: o piloto de `plano-projeto.md` E2.4, já com o que as rodadas de 2026-09-19 acrescentaram. **Concorrentes** (L2d): `mgcv::gam` com `s(u_m, by = x_ℓ)`; B-splines mais group LASSO (`grpreg`); o spline adaptativo de Wang, Jiang & Liu (2024); o **block LASSO de Klopp & Pensky** no mesmo desenho (`grpreg`/`gglasso`), que é a pergunta que D18 convida; o VCBART; e a regressão linear oráculo. **Cenários** de `dgp.R`, com `s'` declarado por D27 (`3/2` no suave, `1/2` no não homogêneo) e o número gravado como atributo do cenário. **Regras de `λ`**: as cinco de E2.3 mais o **QUT** de Giacobino et al. (2017) (L2c), e o oráculo da grade como coluna permanente, que é o que dá sentido a "custo". **Margem**: medir ISE e `λ_min(G_eps)` em `eps ∈ {0, 2^{−(J+1)}, 1.9^{−J}, 0.02, 0.05, 0.10}` e na família `a(L−1)2^{−J}` com `a` em torno de 1; o resultado fixa `wafc_eps_periodic`, hoje provisório em `0.05`. **Consertar de passagem** as três linhas de `wafc/R/tune.R` que descrevem o padrão antigo de `eps` (Roxygen de `cv.wafc(J = )`, a guarda da regra da teoria que para em erro quando `J_n = 1`, e o comentário de `wafc_J_grid()`) e decidir, com número, se a grade passa a incluir `J = 1`. 50 réplicas; `n ∈ {250, 500, 1000}` | E2.1b, E2.2 e E2.3, fechadas; D19, D20, D23, D26, D27 | `wafc/R/competitors.R`, `wafc/scripts/04-pilot.R`, `wafc/tests/test-competitors.R`, `wafc/R/tune.R` (só as três linhas), `wafc/README.md` (só a tabela), `docs/handoff-E2.4.md` |
| E1.7a | **sondagem, não prova** (D28): decidir se a saída (a) de [`selecao-estrutura.md`](selecao-estrutura.md) é viável, em duas frentes. **(i) Álgebra:** sob `X ⊥ U` a Gram é `Σ = Π(Ω ⊗ Σ_Ψ)Π'` com `Ω = E(XX')` (eq. 1.8 de Klopp & Pensky, recordada em `03-desenho-produtos.tex`); com os grupos sendo os blocos `(ℓ,m)` e a base ortonormal dentro de cada bloco, **a condição de irrepresentabilidade em grupos se reduz a uma condição sobre `Ω` sozinha?** É conjectura do chat principal, não resultado; ou se prova a redução, ou se exibe o obstáculo. Seguir a forma da condição em Bach (2008) ou Wei & Huang (2010), conferindo a referência antes de citá-la. **(ii) Numérica:** com o suporte verdadeiro conhecido nos cenários de `dgp.R`, calcular a condição em grupos e reportar se vale, com que folga, e como ela degrada quando `J` cresce, quando `p q` cresce e quando `X` depende de `U` (cenário B de `check/03-desenho-produtos.R`). **Entregável: um veredito explícito**, "vale a pena tentar" ou "não fecha, e eis onde", com os números. Sem teorema e sem numeração global | E1.4 e E2.1, fechadas; D28 | `derivations/06a-sondagem-irrepresentabilidade.md`, `derivations/check/06a-irrepresentabilidade.R`, `docs/handoff-E1.7a.md` |
| E1.7c | `derivations/06-selecao-limiar.tex` (+ `.pdf`) e `check/06-selecao-limiar.R`: a saída (c) de [`selecao-estrutura.md`](selecao-estrutura.md), que entra no artigo de qualquer forma (D28). Enunciar e provar que, com `Ŝ = {(ℓ,m) : ‖ĝ_{ℓm}‖_{L_2} > t_n}` e uma **hipótese de separação em nível de bloco** (`min_{(ℓ,m) ∈ S} ‖g_{ℓm}‖ ≫ r_n`, com `r_n` a taxa do Corolário 4 de E1.6), vale `P(Ŝ = S) → 1` para `t_n` entre `r_n` e a separação. **Não precisa de irrepresentabilidade nem de condição de desenho nova**: consome o Corolário 4 como está, e vale para o **LASSO puro**, que é o estimador base (D3). O enunciado tem de dizer com todas as letras que é estimação seguida de limiar, e não seleção pelo estimador, e a hipótese de separação tem de aparecer como hipótese. Registrar em "o que isto não cobre" que `t_n` depende de constantes desconhecidas na prática, e que a escolha empírica é de E2.4/E2.5. Numeração global: **Corolário 8** (o 6 e o 7 foram usados por E1.8). A conferência mede, nos cenários de `dgp.R`, a curva de acerto de estrutura contra o limiar, e compara o LASSO limiarizado com o sparse group LASSO no mesmo ajuste (`wafc_blocks()` já devolve a norma por bloco) | E1.5 e E1.6, fechadas; D28 | `derivations/06-selecao-limiar.tex`, `derivations/06-selecao-limiar.pdf`, `derivations/check/06-selecao-limiar.R`, `docs/handoff-E1.7c.md` |
| E1.4c | traduzir `derivations/03-desenho-produtos.tex` do inglês para o português (D29), recompilando o `.pdf`. **Só o idioma da prosa**: enunciados, hipóteses, provas, constantes e numeração ficam como estão, e a tradução não pode alterar nenhuma afirmação matemática. Os termos técnicos seguem os outros arquivos em português (`02`, `04`, `05`, `07`), que são o padrão. Motivo: com D29 o `macros.tex` imprime os ambientes em português, e hoje o `03` sai com cabeçalho "Proposição" sobre prosa em inglês. Prioridade baixa: nada depende dela | D29 | `derivations/03-desenho-produtos.tex`, `derivations/03-desenho-produtos.pdf`, `docs/handoff-E1.4c.md` |
| E6.1a | **sondagem, com veredito**: qual base sustenta a aplicação de E6. O critério não é o usual — a base precisa ter um efeito cuja **modulação seja não homogênea**, porque se `β_ℓ(u)` for suave o spline ganha e o argumento do artigo cai (`proposta-metodo.md` §6; Xue & Yang 2006 provam a taxa ótima univariada **sob `α_ls ∈ C^{p+1}[0,1]`**, isto é, na classe em que splines são ótimos). Levantar **duas ou três** candidatas públicas e citáveis, com `n` na casa dos milhares, e para cada uma registrar: fonte, licença, entrada de citação conferida (como em L1, L3 e L4), `n`, e o mapeamento explícito `Y`, `X_ℓ`, `U_m`. Candidatas sugeridas, a confirmar ou substituir: demanda horária de energia ou aluguel de bicicletas (efeito da temperatura modulado por hora do dia e dia do ano, onde a quina do início e do fim da jornada é substantiva); preço de imóvel com efeito da área modulado por latitude e longitude (salto em fronteira de bairro; a aditividade em coordenadas é hipótese forte e tem de ser discutida); poluição horária (efeito do vento modulado por hora e umidade). **A sondagem é empírica:** ajustar o WAFC (`cv.wafc`, padrões atuais) e o `mgcv::gam` com `s(u_m, by = x_ℓ)` em cada candidata, comparar erro de predição fora da amostra e **olhar a estrutura das componentes estimadas** — concentração de energia nos níveis finos, quinas visíveis, e `wafc_blocks()` para ver quais blocos ficam ativos. Veredito por base: sustenta o argumento, é neutra, ou favorece o spline. Os dados **não são versionados**: vão para `wafc/cache/`, com URL e soma de verificação registradas no `.md` | nada; usa `wafc/` como está | `docs/aplicacao-candidatas.md`, `wafc/scripts/05-sondagem-aplicacao.R`, `docs/handoff-E6.1a.md` |
Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**Abertas: E2.4 (em curso), E1.7a, E1.7c, E1.4c e E6.1a.** Arquivos
disjuntos, então correm juntas; E2.4 e E6.1a só se cruzam em `wafc/`, e cada
uma escreve os seus (`04-pilot.R` e `05-sondagem-aplicacao.R`).

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

**Este mapa é a autoridade**, e é ele que E5a usa ao montar o manuscrito. O
`04-oraculo.tex` já imprime o número global (via `\setcounter` no
preâmbulo), prática adotada daqui em diante; `02` e `03` ainda imprimem o
contador local. O próximo resultado numera a partir de Proposição 6, Lema 12,
Teorema 2 e Corolário 7. **Atenção:** E1.7c foi catalogada com "Corolário 6",
número que E1.8 já usou; ela passa a ser **Corolário 8**.

**O enunciado que vai ao resumo do artigo é o Corolário 5** (D16, a
ratificar).

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
