# Executar uma tarefa

Protocolo do executor. A mensagem de abertura é uma linha, e o modelo e o
esforço da conversa são os que o cartão diz:

> Leia `docs/TAREFA.md` e execute o cartão `docs/tarefas/E9.1a.md`.

## 1. O que ler

1. [`../CLAUDE.md`](../CLAUDE.md): as duas regras (sem coautoria; respostas curtas).
2. O cartão. Ele é a tarefa inteira: objetivo, o que ler, o que escrever,
   passos, critério de saída.
3. Só o que o cartão lista em "lê", nos intervalos de linha dados. Do
   [`instrucoes.md`](instrucoes.md), só as seções que o cartão citar (§3 e
   §4 para `.tex`; §5 para matemática; §6 para código).

**Não se lê:** `docs/estado/` inteiro (consulta-se por `grep` do
identificador, `D56` ou `E6.1b`, e lê-se o trecho); `docs/plano-fase2.md`
(é do orquestrador); versões `k` anteriores do manuscrito; PDFs compilados;
`wafc/cache/`; CSV de réplicas; `docs/parecer-wafc.html`.

## 2. Regras

- **Escreve só nos arquivos da lista "escreve".** O que cair fora dela
  (símbolo novo, decisão, arquivo comum a outra tarefa) vira proposta no
  handoff, não edição.
- **Não commita, não faz push, não faz `git add`.**
- **Rodada de mais de alguns minutos não começa aqui.** A tarefa prepara,
  mede o custo em pequena escala, escreve o comando no handoff e para. A
  rodada é um script com o PID gravado e um marcador de fim, e só começa com
  o aviso do autor; nenhum modelo fica à espera dela.
- **Saída de comando no contexto: até 50 linhas.** Log vai a arquivo e
  lê-se o fim. Tabela de resultados lê-se por
  `Rscript wafc/scripts/15-digest.R <csv> [--cells a,b] [--n 500,1000]
  [--methods x,y] [--cols c1,c2] [--ratio M]` (até 40 linhas).
- **Número não se digita:** vem da tabela, por script, e confere-se por script.
- **Dúvida de intenção vira pergunta**, no handoff ou ao autor.
- **Critério de saída falhou duas vezes:** para e relata a causa; o
  orquestrador reescreve o cartão um nível acima.
- **No manuscrito:** a versão viva e a cor da rodada estão na cabeça do
  [`ESTADO.md`](ESTADO.md); nenhuma edição é silenciosa (`instrucoes.md` §3).

## 3. Handoff

A tarefa termina escrevendo `docs/handoff-<ID>.md`, com até 60 linhas:

1. **O que foi feito:** os arquivos.
2. **O que a conferência mostrou:** os números numa tabela curta, o comando
   que os produz e o caminho do arquivo completo.
3. **O que entra no ESTADO:** uma a três linhas; decisão ou lição proposta.
4. **Perguntas em aberto.**

A mensagem final tem até 10 linhas: o veredito contra o critério de saída e
o caminho do handoff.

## 4. Onde está o resto

Decisões: `estado/decisoes.md`, por `grep '^| D56 '`. O que uma etapa mediu:
`estado/diario-AAAA-MM.md`. Numeração dos resultados:
[`numeracao-resultados.md`](numeracao-resultados.md). Notação:
[`notacao.md`](notacao.md). Máquinas e pacotes: [`CONTINUAR.md`](CONTINUAR.md).
