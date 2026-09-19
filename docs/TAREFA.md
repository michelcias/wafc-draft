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
| E2.3 sintonia | **fechada** (2026-09-19): `cv.wafc()`, BIC, EBIC e a regra da teoria; 402 testes passam | `wafc/R/tune.R`, `wafc/scripts/02-tune.R` |
| E5a manuscrito, Seções 1 a 4 | **fechada** (2026-09-19): `k = 1` compila limpo, 28 e 26 páginas | `manuscript/ms_1.tex`, `supp_1.tex` |
| E1.3b emenda de extensão | **em curso** (aberta em 2026-09-19) | catálogo da §3 |
| L4 referências de software | não aberta; **pode abrir já** | catálogo da §3 |
| E1.7, E2.4, E2.5, E3, E4, E5b, E6, E7 | não abertas; E2.4 espera a ratificação de D19 e D20 | |
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

| E1.3b | emenda a `derivations/02-aproximacao-besov.tex` (+ `.pdf`): um **lema de extensão** e a observação de que a Proposição 2 (o custo da periodização) é o caso `eps = 0`. O enunciado a provar: se `U_m` tem suporte em `[eps, 1 − eps]` com `eps > 0` e `g_{ℓm}` pertence à classe de Besov **no suporte**, então existe extensão `g̃` a `[0,1]` que emenda em `0 ≡ 1`, está na mesma classe com norma `≤ C(eps)‖g‖`, e o viés em `L_2(P_U)` volta a `O(2^{−2Js'})`, sem o termo `2^{−J}` da Proposição 2. O operador de extensão é o que o próprio arquivo já cita (Triebel 1983, cap. 3; Cohen 2003, §3.9); registrar como a constante degrada quando `eps → 0`, que é o preço da margem. Registrar também, em observação, que a mesma margem muda a constante de E1.4: com suporte próprio, `E[(a'ψ(U))²] ≥ c_U λ_min(G_eps)‖a‖²`, com `G_eps` a Gram da base **restrita ao suporte**, e não `c_U‖a‖²` — quem mede `λ_min(G_eps)` contra `eps` e `J` é E2.4. Numeração global: o lema novo é o **Lema 10**. `check/02-aproximacao-besov.R` ganha uma seção que mede, numa `g` que não emenda (`u − 1/2`, `e^u`), o erro de projeção **restrito a `[eps, 1−eps]`** contra o erro em `[0,1]`, em `eps ∈ {0, 2^{−J−1}, 1.9^{−J}}`, mostrando a taxa de `1/2` bit por nível virar a taxa cheia; e `λ_min(G_eps)` nos mesmos `eps` | E1.3 e E1.4, fechadas; a origem é a discussão de 2026-09-19 registrada em D23 | `derivations/02-aproximacao-besov.tex`, `derivations/02-aproximacao-besov.pdf`, `derivations/check/02-aproximacao-besov.R`, `docs/handoff-E1.3b.md` |
| L4 | frente curta de verificação, nos moldes de L1 e L3, para as referências que o manuscrito precisa citar e que não estão em `docs/referencias-verificadas.bib`: **Tibshirani (1996)** para o LASSO; **Friedman, Hastie & Tibshirani (2010, JSS)** para o `glmnet`; a **citação do R** (`citation()`, com o ano da versão usada aqui, 4.6.1); o **`sparsegl`** (Liang, Cohen, Sólon Heinsfeld, Pestilli & McDonald, 2024, JSS, pela linha de `inventario-codigo.md` §4, a conferir); o **`WaveBased`**, que é pacote do autor e não tem veículo, e portanto entra como software com versão e URL do repositório, não como artigo; e **Johnstone**, *Gaussian Estimation: Sequence and Wavelet Models*, que o WALL cita como `johnstone2019gaussian` e que a observação sobre a referência minimax quer ao lado de Donoho & Johnstone (1998) — copiar do `references_theo_1.bib` do WALL sem redigitar, se estiver lá. Cada entrada conferida no Crossref (ou na fonte oficial do pacote, quando não houver DOI), com o comentário de conferência acima da entrada, como em L1 e L3. Registrar em `literatura.md` só o que for trabalho, não o que for software | nada; L1 e L3 fecharam | `docs/referencias-verificadas.bib`, `docs/literatura.md` (só linhas novas), `docs/handoff-L4.md` |
Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**Abertas agora: E1.3b (em curso) e L4.** Elas não compartilham arquivo. A próxima depois dela é E2.4 (piloto), que espera a
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
