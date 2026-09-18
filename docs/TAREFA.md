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
| E1.5, E2.1 | não abertas; **podem abrir já** | catálogo da §3 |
| E1.6, E1.7, E2.2 a E2.5, E3 a E7 | não abertas; E1.6 abre quando E1.5 fechar | |
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
| E1.5 | `derivations/04-oraculo.tex` (+ `.pdf`): desigualdade oráculo para o LASSO com perda quadrática e erro sub-gaussiano, `c_ℓ` não penalizados, `λ ≍ σ ‖Z‖_max sqrt(log(p q N_J)/n)`, taxa lenta (sem condição de desenho) e taxa rápida (com E1.4); o viés entra pelo Corolário 1 de E1.3, e a constante de compatibilidade pela forma consumível de E1.4 (`φ_0²(S; Σ̂) ≥ κ_1 c_U / 2`, sem cone); `check/04-oraculo.R` medindo o erro de predição contra `λ² s_0` ao variar `n`, com `s_0` pequeno | E1.3 e E1.4, fechadas; segue Bühlmann & van de Geer (2011), Teoremas 6.1 e 6.2 | `derivations/04-oraculo.tex`, `derivations/04-oraculo.pdf`, `derivations/check/04-oraculo.R`, `docs/handoff-E1.5.md` |
| E2.1 | `wafc/R/load.R` (carrega `R/*.R`, declara dependências), `wafc/R/dgp.R` (cenários suave, não homogêneo, nulo; `simulate_wafc()`), `wafc/R/design.R` (`wafc_design()` sobre `wbasis()` do `WaveBased`, blocos `X_j ⊙ ψ(U_k)`, colunas nomeadas, `penalty.factor`, versão esparsa), `wafc/tests/test-design.R` (com `θ*` na base e sem ruído, `glmnet` com `λ → 0` recupera `θ*`; posto cheio; nomes das colunas) e `wafc/scripts/01-smoke.R` | o enunciado de E1.2 (a hipótese de identificabilidade fixa o que o desenho descarta); a ordem das colunas é D12 | `wafc/R/load.R`, `wafc/R/dgp.R`, `wafc/R/design.R`, `wafc/tests/test-design.R`, `wafc/scripts/01-smoke.R`, `wafc/README.md` (só a tabela), `docs/handoff-E2.1.md` |

Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**Quem abrir E1.5** lê antes o Corolário 1 de
`derivations/02-aproximacao-besov.tex` e a Proposição de
`derivations/03-desenho-produtos.tex`, além da §1 de
[`busca-novidade.md`](busca-novidade.md): a desigualdade oráculo de Klopp &
Pensky (2015) cobre o caso `q = 1` com `X ⊥ U`, e o que é novo aqui é o
desenho aditivo.

**Quem abrir E2.1** implementa a ordem de colunas de D12 e descarta o que
E1.2 manda descartar: no periódico com `j_0 = 0`, a coluna `φ_{00}` do
`wbasis()` e nada mais. A construção `Z = (X ⊗ Ψ(U))[, perm]` do
`check/03-desenho-produtos.R` (`kron_rows()`, `perm_D12()`) serve de teste
de referência.

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

Os arquivos mantêm os contadores locais do LaTeX; **este mapa é a
autoridade**, e é ele que E5a usa ao montar o manuscrito. O próximo
resultado numera a partir de Proposição 3, Lema 3 e Corolário 1.

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
