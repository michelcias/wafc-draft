# wafc-draft

Repositório de trabalho do projeto **WAFC: regressão com coeficientes
funcionais aditivos, estimados por LASSO em bases de wavelets** (nome
ratificado em D8), com alvo de publicação em periódico Q1 de Statistics and
Probability (SCImago). Alvo primário: *Statistica Sinica* (D5);
reserva: *Electronic Journal of Statistics*. Detalhe em
[`docs/alvo-revista.md`](docs/alvo-revista.md).

**As duas regras que valem desde a primeira mensagem de qualquer conversa:**

1. **Commit e push nunca levam coautoria.** Sem `Co-Authored-By:`, sem
   "Generated with", sem trailer que nomeie o assistente. Vale para as
   mensagens redigidas, para os comandos passados ao autor e para os commits
   executados a pedido dele. Se a ferramenta injetar um trailer, remova antes
   de commitar; se um trailer escapar, não reescreva o commit, avise.
2. **As respostas no chat são as mais curtas possível.** Veredito, a razão em
   uma ou duas linhas, e o bloco de código ou `.tex` quando couber. Nada de
   reimprimir contexto que o autor já tem, enumerar alternativas descartadas
   ou anunciar o que vai fazer antes de fazer.

Depois disso, depende do papel da conversa (`docs/instrucoes.md`, §7):

- **Orquestrador** (o chat principal): leia [`docs/ESTADO.md`](docs/ESTADO.md),
  que é uma cabeça curta (onde o trabalho está, o quadro da fase, o próximo
  passo), e do [`docs/instrucoes.md`](docs/instrucoes.md) as seções do que for
  fazer.
- **Executor de uma tarefa:** leia [`docs/TAREFA.md`](docs/TAREFA.md) e o
  cartão indicado (`docs/tarefas/<ID>.md`). O cartão é a tarefa inteira; o
  estado não se lê.

O arquivo morto (`docs/estado/`: diário, decisões, perguntas fechadas) se
consulta por `grep` do identificador, nunca inteiro. O `instrucoes.md` é o
documento normativo (fluxo, escrita, marcação, git); em caso de conflito,
ele manda.

## Resumo das demais regras

### Git

- **O commit e o push são do autor por padrão.** Ler o estado (`git status`,
  `git log`, `git diff`) é livre; alterar o repositório por iniciativa própria,
  não. Terminar uma rodada de edições não autoriza commitar. Quando o autor
  pedir, executar sem devolver a pergunta; o pedido vale para aquele commit.
- **Nunca reescrever histórico**: `rebase`, `commit --amend`, `push --force`,
  `reset --hard`, apagar branch. Se parecer necessário, dizer o comando e o
  porquê, e esperar.
- Mensagens no imperativo, em inglês, com prefixo de escopo (`docs:`, `math:`,
  `proto:`, `ms:`, `plan:`); corpo explica o porquê.

### Continuidade

- **Etapa fechada é etapa registrada**, com os números que a fecharam, antes
  de a rodada terminar: um bloco em `docs/estado/diario-AAAA-MM.md` e a
  linha do quadro na cabeça do `docs/ESTADO.md`. Dependência nova de script
  entra no `docs/CONTINUAR.md` na mesma rodada. Regra completa em
  [`docs/instrucoes.md`](docs/instrucoes.md), §8.
- Edição nesses arquivos é **verificada**: substituição por âncora falha em
  silêncio quando a âncora mudou.

### Execução e economia de contexto (fase 2, D83)

- **Orquestrador e executores são no máximo Opus.** O modelo e o esforço de
  cada tarefa são os do cartão; os papéis estão em `.claude/agents/`.
- **Lê-se por lista fechada:** o cartão diz os arquivos e os intervalos de
  linha. Não se leem versões `k` anteriores do manuscrito, PDFs compilados,
  `wafc/cache/` nem CSV de réplicas.
- **Saída de comando no contexto: até 50 linhas.** Tabela de resultados se
  lê por `Rscript wafc/scripts/15-digest.R`.
- **Rodada longa é script**, com o aviso do autor antes; nenhum modelo fica
  à espera.

### Manuscrito

- Antes de aplicar qualquer alteração ao `.tex`, **perguntar se deve criar a
  versão `{k+1}`**. `ms_{k}.tex`, `supp_{k}.tex` e `references_{k}.bib` andam
  juntos. A pergunta é uma vez por rodada.
- **Sem edição silenciosa** depois que um `.tex` existe: remoção em
  `\textcolor{gray}{\sout{...}}`, inclusão em `\textcolor{colR1}{...}`.
- **Inglês americano**, registro acadêmico formal, sem negrito e sem em-dash
  (—) como pausa.
- **Escopo amplo, mão conservadora.** Levantar tudo que puder elevar o trecho,
  inclusive substância; o que for duvidoso ou substantivo vira pergunta no
  chat antes, não reescrita.
- O texto sugerido e a versão `.tex` vão **na própria resposta**, em blocos de
  código, nunca em artefato ou arquivo à parte.

### Matemática e código

- Toda afirmação de taxa, de condição de desenho ou de propriedade do
  estimador vem com a razão ou a referência. Se não houver, é sinalizada
  como conjectura.
- Referência bibliográfica só entra no `.bib` depois de verificada (título,
  autores, ano, veículo). Uma citação de memória é marcada `[VERIFICAR]` até lá.
- O código do método vive na pasta dedicada `wafc/` deste repositório
  (decisão D4): funções em `wafc/R/`, testes em `wafc/tests/`, scripts em
  `wafc/scripts/`. O `WaveBased` **não recebe código** por enquanto; o
  pacote instalado é só uma dependência para avaliar as bases. O estudo de
  simulação e a aplicação vão ao compêndio `wafc-studies`. Aqui ficam
  também a conferência numérica (`derivations/check/`) e o que documenta a
  reprodutibilidade.

## Repositórios vizinhos

| Repositório | Papel |
|---|---|
| `../WaveBased` | pacote R/C do autor (`michelcias/WaveBased`); fornece as bases de wavelets (`wbasis()`, `wtable()`); o `wall()` é referência de leitura para o desenho por blocos. Não recebe código deste projeto por enquanto (D4) |
| `../wall-manuscript` | artigos do WALL (teórico e aplicado); o teórico é o molde da arquitetura de prova (sieve, Besov, desigualdade oráculo, compressibilidade) |
| `../wall` | compêndio de benchmark do WALL; molde do `wafc-studies` |
| `../bdm-draft` | projeto irmão; a origem destas convenções de trabalho |
| `wafc-draft` (este) | manuscrito, documentos de trabalho, derivações, o código do método (`wafc/`) e o plano |
| `wafc-studies/` (pasta deste repositório, D60) | compêndio de simulação e aplicação; vira repositório próprio só no fim, por cópia, com só o que reproduz o manuscrito e o supp |
