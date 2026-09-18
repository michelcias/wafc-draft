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
| E1.2, E1.3, E1.4 | não abertas; **podem abrir já** (E1.1 fechou) | catálogo da §3 |
| E1.5 a E1.7, E2, E3, E4, E5, E6, E7 | não abertas; E2.1 abre com o enunciado de E1.2 | |
| L1, L2 | não abertas; **podem abrir já** | `literatura.md` |

Esta tabela é atualizada pelo chat principal quando uma etapa fecha; em caso
de dúvida, o `ESTADO.md` manda.

## 3. Catálogo das tarefas abertas

Cada linha diz o que entregar, de que depende, e **os únicos arquivos que a
tarefa pode criar ou editar**.

| Tarefa | Entregável | Depende de | Arquivos permitidos |
|---|---|---|---|
| L1 | `docs/referencias-verificadas.bib` com cada entrada conferida no Crossref (título, autores, ano, veículo, volume, páginas, DOI) e sem a nota `[L1: confirmar]`; cada `[VERIFICAR]` de `docs/literatura.md` resolvido (status `verificado` com a entrada no `.bib`, ou o motivo de não ter achado); entradas `no wall` copiadas de `../../wall-manuscript/manuscript/theo/references_theo_1.bib` sem redigitar | nada | `docs/referencias-verificadas.bib`, `docs/literatura.md` (só a coluna de status e as linhas de bibliografia), `docs/handoff-L1.md` |
| L2 | `docs/busca-novidade.md`: as quatro buscas de `literatura.md` ("Buscas pendentes"), uma tabela por busca (trabalho, o que faz, o que não faz, ameaça à novidade, status), leitura do PDF de Sardy & Ma (2024) com o que a teoria deles cobre, varredura da *Statistica Sinica* e da EJS desde 2015, e um veredito por contribuição de `alvo-revista.md` §4 | nada | `docs/busca-novidade.md`, `docs/literatura.md` (só linhas novas nas tabelas), `docs/handoff-L2.md` |
| E0.3 (autor) | só falta o quartil 2024 da SS conferido no SCImago e anotado em `alvo-revista.md` §1; a transcrição e o template estão feitos | nada | `docs/alvo-revista.md` (§1), `docs/handoff-E0.3.md` |
| E1.2 | `derivations/01-identificabilidade.md` (enunciado, hipóteses, prova, "o que não cobre") e `derivations/check/01-identificabilidade.R` imprimindo `OK`: posto cheio do desenho em `n` pequeno e recuperação de coeficientes conhecidos; o que muda com `boundary = "interval"` | E1.1, fechada: a notação de `notacao.md` é para seguir como está, e símbolo novo vai ao handoff, não ao arquivo | `derivations/01-identificabilidade.md`, `derivations/check/01-identificabilidade.R`, `docs/handoff-E1.2.md` |
| E1.3 | `derivations/02-aproximacao-besov.tex` (+ `.pdf`) com o lema de aproximação em `L_2(P_U)` e, se preciso, em `L_∞`, transposto do Passo 1 do WALL com a diferença por periodização registrada; `check/02-aproximacao-besov.R` medindo o erro de projeção em `bumps` e `sin` contra `2^{−J s'}` | E1.1, fechada (idem) | `derivations/02-aproximacao-besov.tex`, `derivations/02-aproximacao-besov.pdf`, `derivations/check/02-aproximacao-besov.R`, `docs/handoff-E1.3.md` |
| E1.4 | `derivations/03-desenho-produtos.tex` (+ `.pdf`): Gram populacional do desenho de produtos; (i) fatoração sob `X ⊥ U`; (ii) autovalor restrito no caso geral sob `λ_min(E[XX' \| U]) ≥ κ_1`; (iii) versão empírica por concentração; `check/03-desenho-produtos.R` com `λ_min` restrito (ou compatibilidade) em `n = 200`, `J = 3`, `X ⊥ U` e `X` dependente de `U` | E1.1, fechada (idem) | `derivations/03-desenho-produtos.tex`, `derivations/03-desenho-produtos.pdf`, `derivations/check/03-desenho-produtos.R`, `docs/handoff-E1.4.md` |
| E2.1 | `wafc/R/load.R` (carrega `R/*.R`, declara dependências), `wafc/R/dgp.R` (cenários suave, não homogêneo, nulo; `simulate_wafc()`), `wafc/R/design.R` (`wafc_design()` sobre `wbasis()` do `WaveBased`, blocos `X_j ⊙ ψ(U_k)`, colunas nomeadas, `penalty.factor`, versão esparsa), `wafc/tests/test-design.R` (com `θ*` na base e sem ruído, `glmnet` com `λ → 0` recupera `θ*`; posto cheio; nomes das colunas) e `wafc/scripts/01-smoke.R` | o enunciado de E1.2 (a hipótese de identificabilidade fixa o que o desenho descarta); a ordem das colunas é D12 | `wafc/R/load.R`, `wafc/R/dgp.R`, `wafc/R/design.R`, `wafc/tests/test-design.R`, `wafc/scripts/01-smoke.R`, `wafc/README.md` (só a tabela), `docs/handoff-E2.1.md` |

Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 tocam `literatura.md` em partes distintas (status de
linhas existentes; linhas novas); se colidirem, L2 deixa as linhas no
handoff.

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
pelo chat principal: **nenhum ainda**.

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
