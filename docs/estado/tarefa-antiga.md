# TAREFA antigo: início e §3 Catálogo
Movido de docs/TAREFA.md em 2026-11-16, sem alteração de linhas.
Arquivo morto: não ler inteiro; consultar por grep pelo identificador (D56, E6.1b, pergunta 33) e ler o trecho.

<!-- conteudo movido abaixo -->
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
   (`../wall-manuscript/manuscript/theo/ms_theo_1.tex`).

## 3. Catálogo das tarefas abertas

Cada linha diz o que entregar, de que depende, e **os únicos arquivos que a
tarefa pode criar ou editar**.

| Tarefa | Entregável | Depende de | Arquivos permitidos |
|---|---|---|---|

Duas tarefas não podem editar o mesmo arquivo ao mesmo tempo; se o
catálogo tiver duas que tocam o mesmo arquivo, a segunda deixa as linhas
no handoff. L1 e L2 fecharam, então nenhuma tarefa aberta encosta no
`literatura.md`.

**E4.3b fechou em 2026-10-08** (pergunta 69); a E4.4 está pronta e espera o aviso. **E4.3b emendada em 2026-10-08** (D79, a regra de tamanho; e D80, a regra só além de `k = 80` e o limite de iterações; o chat da tarefa relê a entrada). **E5k fechou em 2026-10-08** (pergunta 68). **E4.3b catalogada em 2026-10-08** (D78); a E4.4 espera por ela. **E5k catalogada em 2026-10-07** (a Figura S1 refeita depois do piloto já está em `results/`). **E6.2b fechou em 2026-10-07** (pergunta 66); a aplicação terminou e o piloto roda desde 16h31. **E6.2b catalogada em 2026-10-07** (D75): a sondagem tem de terminar antes do fim da `beijing.heat`, senão atrasa o piloto. **E5j fechou em 2026-10-07** (pergunta 65). **E5j catalogada em 2026-10-07** (D73; sem ajuste, só o relatório de marylebone, enquanto a cadeia roda a `beijing.heat`). **E5i fechou em 2026-10-06** (pergunta 63). **L14 fechou em 2026-10-06** (pergunta 62). **E5h fechou em 2026-10-06** (pergunta 60; o manuscrito vivo é a `k = 5`). **E4.3a e E5h catalogadas em 2026-10-06**, sem arquivo em comum e sem rodada: a E4.3a só prepara o piloto (a rodada espera o fim da E6.2 e o aviso, sozinha na máquina, D67), a E5h só escreve na `k = 5`. **E6.2 com o código pronto em 2026-10-06**; as rodadas esperam o aviso (pergunta 59). **E4.1c fechou em 2026-10-06**; a E4.3 pode ser catalogada. **L13 fechou em 2026-10-06** (pergunta 58). **L13 catalogada em 2026-10-06**, sem arquivo em comum com a E4.1c e a E6.2. **E6.2 catalogada em 2026-10-06**, em paralelo à E4.1c, sem arquivo em comum: a E6.2 só cria arquivos novos no compêndio e deixa no handoff o que muda no `README.md`, no `run_all.R` e no `renv.lock` dele; se precisar mexer nos `R/` comuns (`cli.R`, `run.R`), também vai ao handoff. **L12 fechou em 2026-10-06** (pergunta 57). **E4.1c catalogada em 2026-10-06** (D63), sem arquivo em comum com a L12; a E4.3 espera por ela. **E4.1 fechou em 2026-10-06** (pendências na pergunta 56; a próxima junção exata inclui o nulo em `n = 250`). **E4.1b fechou em 2026-10-06** (o braço passa `x_u_rho = 0.5` também à amostra de teste). **E4.1 e E4.1b catalogadas em 2026-10-06** (D60), sem arquivo em comum: a E4.1 só cria dentro de `wafc-studies/` e lê `wafc/`; a E4.1b só toca o `dgp.R`, o teste dele e talvez o `wafc/README.md`. A E4.1 pode escrever a configuração do braço com `x_u_rho` antes, e a fumaça desse braço espera a E4.1b; a junção exata da E4.1 não depende dela, porque o padrão `x_u_rho = 0` reproduz os sorteios. Nenhuma das duas faz rodada longa.

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

