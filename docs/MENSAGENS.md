# As mensagens que eu mando

Cola de uso do autor. Copiar e colar; nada aqui é normativo.
A regra por trás está em [`instrucoes.md`](instrucoes.md), §7 e §8.

---

## Abrir a conversa principal

Uma por vez. É a única que edita `ESTADO.md`, `TAREFA.md`, `notacao.md`,
`plano-projeto.md` e a única que commita.

> Leia `CLAUDE.md`, `docs/ESTADO.md` e `docs/instrucoes.md` para se situar.
> Hoje vamos <foco do dia>.

## Abrir em outra máquina

Depois de clonar `wafc-draft` (o resto o arquivo explica).

> Leia `docs/CONTINUAR.md` e me diga o que falta instalar; depois leia
> `docs/ESTADO.md` e continuamos de onde parou.

## Abrir uma conversa em paralelo

Uma por tarefa, na mesma pasta. A tarefa tem de estar no catálogo de
[`TAREFA.md`](TAREFA.md), que diz o que ela entrega e em quais arquivos pode
mexer.

> Leia `docs/TAREFA.md` e execute a tarefa L2

Trocando o nome da tarefa. Para saber quais estão abertas agora, é a tabela
da §2 e o catálogo da §3 do `TAREFA.md`, ou pergunte na conversa principal:

> Quais tarefas posso rodar em paralelo agora?

## Pedir uma tarefa que não está no catálogo

Na conversa principal, antes de abrir o chat da tarefa:

> Catalogue <o que você quer> como tarefa.

Ela volta com o nome, o que a tarefa entrega e a mensagem de abertura.

## Quando uma tarefa termina

A conversa de tarefa escreve `docs/handoff-<nome>.md` e para. Avise a
principal:

> <nome da tarefa> terminou

Ela integra no `ESTADO.md`, apaga o handoff e commita.

## Ratificar uma decisão

Na conversa principal:

> Ratifico D4 e D5. D8: prefiro <nome>.

Ela move a linha para a tabela de decisões e não reabre.

## Fechar o dia, ou trocar de máquina

> Commita e faz push. Confira se está tudo documentado para eu continuar de
> outra máquina.

---

## O que só eu posso fazer

- **Commit e push** são meus por padrão; a conversa só executa quando peço.
- **Decisões** ficam listadas na §4 do `ESTADO.md`. Quando eu respondo, elas
  viram uma linha na tabela de decisões e não se reabrem.
- **Páginas inacessíveis da máquina** (site da *Statistica Sinica*, SCImago)
  e PDFs atrás de login.
- **Criar repositório** e qualquer coisa fora desta árvore: a conversa
  pergunta antes.
