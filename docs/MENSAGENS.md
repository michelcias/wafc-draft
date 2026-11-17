# As mensagens que eu mando

Cola de uso do autor. Copiar e colar; nada aqui é normativo.
A regra por trás está em [`instrucoes.md`](instrucoes.md), §7 e §8.

---

## Abrir o orquestrador

Um por vez, em **Opus, esforço médio**. É o único que edita a cabeça do
`ESTADO.md`, os arquivos de `docs/estado/`, os cartões e o `plano-fase2.md`,
e o único que commita. Cada sessão é um lote e termina; a seguinte começa de
contexto novo.

> Leia `CLAUDE.md` e `docs/ESTADO.md`. Integre os handoffs pendentes e
> despache o próximo lote.

## Abrir uma tarefa

Uma por cartão, na mesma pasta, com o **modelo e o esforço que a segunda
linha do cartão diz**.

> Leia `docs/TAREFA.md` e execute o cartão `docs/tarefas/E9.1a.md`.

Para saber quais estão abertas: o quadro na cabeça do `ESTADO.md` e a pasta
`docs/tarefas/`.

## Pedir uma tarefa que não está no quadro

No orquestrador:

> Escreva o cartão de <o que você quer>.

Ele volta com o nome, o modelo, o esforço e a linha de abertura.

## Dar o aviso de uma rodada

A tarefa de preparação termina com o custo medido e o comando. A rodada só
começa com:

> Pode rodar a E9.1b.

## Quando uma tarefa termina

Ela escreve `docs/handoff-<ID>.md` e para. No orquestrador:

> E9.1a terminou.

Ele integra (diário do mês, decisões, quadro), apaga o handoff e o cartão.

## Ratificar uma decisão

No orquestrador:

> Ratifico a pergunta 78. Na 79 prefiro <...>.

Ela vira uma linha em `docs/estado/decisoes.md` e não se reabre.

## Abrir em outra máquina

Depois de clonar `wafc-draft` (o resto o arquivo explica).

> Leia `docs/CONTINUAR.md` e me diga o que falta instalar; depois leia
> `docs/ESTADO.md` e continuamos de onde parou.

## Fechar o dia, ou trocar de máquina

> Commita e faz push. Confira se está tudo documentado para eu continuar de
> outra máquina.

---

## O que só eu posso fazer

- **Commit e push** são meus por padrão; a conversa só executa quando peço.
- **Decisões:** as perguntas que me esperam ficam na cabeça do `ESTADO.md`.
  Quando eu respondo, viram linha em `docs/estado/decisoes.md`.
- **O aviso de rodada longa** e a escolha da máquina (o tempo que vai a
  tabela vem só da desktop, sozinha).
- **A rota** depois dos quatro gates (`plano-fase2.md` §3).
- **Páginas inacessíveis da máquina** (site da *Statistica Sinica*, SCImago)
  e PDFs atrás de login.
- **Criar repositório** e qualquer coisa fora desta árvore: a conversa
  pergunta antes.
- **A leitura linha a linha das provas** é de um coautor (D90).
