# Instruções de trabalho

Documento normativo do projeto. Em caso de conflito com qualquer outro arquivo
deste repositório, este manda. Adaptado do `instrucoes.md` do `bdm-draft`, com
o que mudou: o alvo é uma revista de estatística geral (não bayesiana), a
teoria é de regressão penalizada (desigualdades oráculo, taxas em Besov), e o
código do método vive na pasta `wafc/` deste repositório, com o `WaveBased`
só como dependência para as bases (D4); empacotar decide-se em E3.3.

---

## 0. As duas regras de abertura

Elas estão no `CLAUDE.md` e são repetidas aqui porque valem antes de qualquer
outra coisa.

**Commits e push sem coautoria.** Nenhuma mensagem de commit leva atribuição de
ferramenta: sem `Co-Authored-By:`, sem "Generated with", sem trailer que nomeie
o assistente. Vale para mensagens redigidas, para comandos passados ao autor e
para commits executados a pedido dele. A ferramenta pode injetar o trailer
sozinha; a verificação é do assistente, antes de commitar. Se um trailer
escapar, **não se reescreve o commit**; avisa-se. A razão é de autoria: o
manuscrito é dos autores, e a ferramenta usada não é um fato que o histórico
precise registrar.

**Respostas curtas.** O padrão de resposta é: veredito, a razão em uma ou duas
linhas, e o bloco de código ou `.tex` quando couber. Não reimprimir contexto
que o autor já tem. Não enumerar alternativas que não serão seguidas. Não
anunciar o que vai fazer antes de fazer. Quando a pergunta é aberta e pede
análise, a análise pode ser longa, mas cada parágrafo precisa carregar algo que
o autor não sabia.

---

## 1. Fluxo de trabalho

### Como a conversa começa

1. Ler `CLAUDE.md`.
2. **Orquestrador:** ler [`ESTADO.md`](ESTADO.md), que é uma cabeça curta
   (onde o trabalho está, o quadro da fase, o próximo passo), os
   `docs/handoff-*.md` pendentes e, do [`plano-fase2.md`](plano-fase2.md), só
   a seção da tarefa em foco. **Executor:** ler [`TAREFA.md`](TAREFA.md) e o
   cartão; o estado não se lê (§7).
3. O arquivo morto (`docs/estado/`: diário, decisões, perguntas fechadas) se
   consulta por `grep` do identificador; não se lê inteiro.
4. Se o foco for o manuscrito, ler o `.tex` da seção em questão **antes** de
   opinar sobre ela.
5. Se o foco for matemática, ler o arquivo correspondente em `derivations/` e
   [`notacao.md`](notacao.md) antes de escrever qualquer fórmula. A prova do
   WALL teórico (`../wall-manuscript/manuscript/theo/ms_theo_1.tex`,
   Seções 5, 7 e 8) é o molde de arquitetura; ler o passo correspondente
   antes de reescrevê-lo para o caso de regressão.

### Perguntas direcionadas e escopo

As perguntas costumam ser direcionadas a um resultado, uma seção ou um
parágrafo. A lista de coisas a verificar (clareza, concisão, correção,
coerência) descreve o padrão, não uma cerca: o alvo e a amplitude de cada
análise são definidos pelo pedido do momento.

### Escopo amplo, mão conservadora

Quanto mais aberta a pergunta, mais largo o trecho pode ser lido, inclusive em
**substância**: um argumento que não fecha, uma afirmação sem apoio, uma taxa
afirmada sem a condição de desenho que a sustenta, uma comparação numérica sem
o competidor óbvio. A amplitude está no que se **nota e levanta**; a
contenção, em **quanto se mexe sem aval**. O substantivo ou duvidoso vira
pergunta no chat antes, não reescrita por conta própria.

### Dúvidas antes de aplicar

Dúvida genuína sobre intenção não se resolve aplicando primeiro e explicando
depois: vira pergunta. Dúvidas estruturais sobre o modo de trabalhar são
resolvidas de uma vez, no início; no decorrer, só o pontual é perguntado.

---

## 2. Git

### Quem commita

**O commit e o push são do autor por padrão; do assistente só quando pedido.**

- Por iniciativa própria, não. Ler o estado (`git status`, `git log`,
  `git diff`) é livre; alterar o repositório, não. Terminar uma rodada não
  autoriza commitar.
- Quando pedido, sim, sem devolver a pergunta. "Commita", "commita e faz push"
  autorizam `git add`, `git commit` e `git push` de fast-forward. O pedido vale
  para aquele commit, não para a sessão.
- Reescrever histórico está fora, sempre: `rebase`, `commit --amend`,
  `push --force`, `reset --hard`, apagar branch. Se parecer necessário, dizer o
  comando e o porquê, e esperar.
- Apagar arquivo rastreado se avisa antes: o que sai e o que deixa de funcionar.
- O `WaveBased` não é tocado por este projeto (D4); se E3.3 mudar isso,
  commit lá segue as mesmas regras.

O que continua esperado é o **lembrete**, uma linha no fecho da resposta, nos
momentos úteis: fim de uma rodada que compila, criação de versão `{k+1}`,
atualização substantiva de documento de trabalho, fim de sessão com trabalho
pendente.

### Mensagens

Assunto no imperativo, em inglês, com prefixo de escopo (`docs:`, `math:`,
`proto:`, `ms:`, `supp:`, `plan:`; no pacote, `feat(wafc):`, `fix(wafc):`,
`test(wafc):`); corpo explica o *porquê*, não o *o quê*, que o diff já
mostra. Sem coautoria (regra 0).

---

## 3. Manuscrito

### Versionamento

`ms_{k}.tex`, `supp_{k}.tex` e `references_{k}.bib` andam juntos com o mesmo
índice `k`. **Antes de aplicar qualquer alteração ao `.tex`, perguntar se deve
criar a versão `{k+1}`:**

- **Sim**: copiar os três, atualizar `\bibliography{references_{k+1}}` nos dois
  `.tex` novos, trabalhar só neles. A versão `{k}` fica intacta.
- **Não**: aplicar na versão `{k}` corrente.

A pergunta é feita uma vez por rodada, não a cada edição. O primeiro `.tex`
nasce como `k = 1`, sem marcação, porque não há versão anterior contra a qual
marcar.

A partir de `k = 7` (D86), o texto vive em `manuscript/k{k}/`, uma seção por
arquivo, e `ms_{k}.tex` e `supp_{k}.tex` guardam o preâmbulo e os `\input`; a
versão `{k+1}` copia também a pasta. Uma tarefa lê só a seção em que mexe.
Quando o autor aceita uma rodada, `manuscript/tools/accept-round.py` tira a
marcação e `split-sections.py` reparte o resultado; a conferência é a do
bloco da E8.4 em `estado/diario-2026-11.md`.

### Marcação de alterações

A partir de `k = 2`, nenhuma edição é silenciosa:

| Finalidade | Comando | Cor |
|---|---|---|
| Remoção sugerida | `\textcolor{gray}{\sout{...}}` | cinza tachado |
| Inserção da rodada corrente | `\textcolor{colR{r}}{...}` | azul (`RGB 0, 90, 200`) |

O índice `r` da cor avança quando o autor aceita uma rodada (despigmenta), e o
`ESTADO.md` diz qual é o corrente. O preâmbulo define as cores; se uma
alteração precisar de pacote ou cor não definidos, sinalizar em vez de assumir.

Cuidados que já custaram compilações:

- Comandos frágeis dentro de `\sout` (`\cite`, `\ref`, `\eqref`) vão em
  `\mbox{}`.
- `\sout` não atravessa matemática deslocada nem quebra de parágrafo. Para uma
  remoção que cruze isso, marcar de outra forma ou perguntar.
- Granularidade adaptativa: reescrita ampla tacha a frase inteira e insere a
  nova; troca pontual marca palavra a palavra.
- `\textcolor` não atravessa parágrafo nem ambiente de teorema: inserção de
  vários parágrafos, ou com enunciado, vai num grupo `{\color{colR1} ... }`
  (E5c).
- `\textcolor` em título de seção quebra a compilação quando o cabeçalho
  corrente põe o título em maiúsculas (a cor vira `COLR1`): título curto em
  `\section[...]{...}` e `\texorpdfstring` para os marcadores (E5c).
- Cor dentro de `\citep{...}` não funciona; citação acrescentada ao lado de
  uma existente vai por `\citetext` com `\citealp` (E5c).
- Remoção de um enunciado inteiro, de uma prova ou de matemática deslocada,
  que o `\sout` não atravessa: o texto antigo fica no lugar num grupo
  `{\color{gray} ...}` aberto por uma linha "[removed in k = r]", sem
  tachado, e o novo vem logo depois em `{\color{colR1} ...}`; trocas
  pontuais seguem com `\sout` (D51, para a troca da teoria em `k = 4`).
- `\textcolor{colR1}{ texto}` perde o espaço inicial (o `\color` termina em
  `\ignorespaces`) e cola a palavra anterior: o espaço vai fora do
  `\textcolor` (E5d).
- O PDF com a marcação superestima o artigo: para o teto de páginas (D21),
  contar numa cópia sem o texto removido (E5d).

Exceção registrada: **uniformização de notação** aplicada por script sobre
lista fechada de padrões pode entrar sem marcação, desde que o `ESTADO.md`
registre o padrão e a contagem.

### Formato da resposta a uma sugestão

Até quatro partes: (1) o porquê, curto; (2) o texto sugerido; (3) a versão
`.tex` com marcação, em bloco de código na própria resposta; (4) o fecho com a
opção de incorporar ou continuar aprimorando. No vai-e-vem, só o pedaço que
muda; quando o autor pede o trecho inteiro como entrega, o trecho inteiro.

---

## 4. Escrita

- **Idioma:** inglês americano no manuscrito; português nos documentos de
  trabalho.
- **Registro:** acadêmico formal, terminologia precisa, tom objetivo. Sem
  negrito no corpo. Sem em-dash (—) como pausa; en-dash em intervalos
  (`pp. 10–20`) é tipografia e fica.
- **Consistência:** preservar a voz do texto; intervenções pontuais, não
  reformulação ampla.
- **Respaldo:** afirmação sobre desempenho, comparação ou propriedade de método
  vem ancorada em citação ou em evidência interna (tabela, prova,
  experimento). Sem apoio é lacuna a sinalizar.
- **Economia de paper:** justificativa proporcional ao peso do ponto; cortar
  andaime; não derrubar alternativas que ninguém proporia; não repetir o que
  outra seção já estabelece; **preservar o que sustenta o argumento** (a
  concessão que antecipa o referee, a razão técnica, a citação-âncora).

### Específico da revista-alvo

Ver [`alvo-revista.md`](alvo-revista.md). O que muda na prática enquanto o
alvo for a *Statistica Sinica*: teto de **40 páginas em espaço duplo no
template, referências e apêndice incluídos** (`ss-instrucoes-autores.md`),
provas no suplementar, o leitor é metodólogo de estatística geral (a exposição de
wavelets precisa ser autossuficiente em uma página; a de LASSO pode ser
econômica), e o padrão da revista para este tipo de artigo é **método +
teoria + simulação + aplicação real**, todos os quatro. Se o alvo mudar para
a EJS, cai o teto de páginas e sobe a expectativa de teoria.

---

## 5. Matemática

- **Um resultado por arquivo** em `derivations/`, com enunciado, hipóteses
  explícitas, prova e um parágrafo de "o que isso não cobre". O `.tex` do
  manuscrito cita o arquivo de origem em comentário (`% derivations/xxx.tex`).
- **A notação é a de [`notacao.md`](notacao.md).** Símbolo novo entra lá
  primeiro; símbolo que muda, muda lá e depois em todo lugar.
- **Toda taxa vem com as três peças**: o erro de aproximação (Besov, nível
  `J`), o erro de estimação (dimensão efetiva, `log` do número de colunas,
  `n`) e a condição de desenho que liga os dois (compatibilidade ou
  autovalor restrito da matriz de produtos `X_ℓ ψ_{jk}(U_m)`). Uma taxa sem a
  condição de desenho é conjectura.
- **Hipótese de desenho é hipótese sobre `(X, U)` conjuntos**, não sobre cada
  um: o que entra na Gram é `E(X_ℓ X_{ℓ'} ψ_{jk}(U_m) ψ_{j'k'}(U_{m'}))`. Escrever a
  hipótese nessa forma e só depois simplificá-la.
- **Conferência numérica antes de prova.** Um resultado novo é primeiro testado
  em R com `n` pequeno (`n ≤ 200`, `J ≤ 3`) na forma densa; só depois se
  escreve a prova. O script fica em `derivations/check/`.
- **Nos scripts de conferência, elemento de lista se acessa com `[[ ]]` e
  nome completo**, nunca com `$` (casamento parcial já produziu falhas
  espúrias no projeto irmão).

---

## 6. Código

O código do método vive **na pasta `wafc/` deste repositório** (decisão D4
em `ESTADO.md`; detalhe em `plano-projeto.md`, E0.2). O `WaveBased` não
recebe código por enquanto: o pacote instalado é uma dependência, como o
`glmnet`, usada só para avaliar as bases (`wbasis()`, `wtable()`).

| Código | Onde |
|---|---|
| funções do método: desenho, ajuste, sintonia, reconstrução, predição, gráficos | `wafc/R/`, um arquivo por tema, carregados por `source()` via `wafc/R/load.R` |
| testes | `wafc/tests/`, `testthat` rodado com `testthat::test_dir("wafc/tests")` |
| scripts de fumaça, piloto e comparação | `wafc/scripts/`, numerados |
| estudo de simulação e aplicações | a pasta `wafc-studies/` deste repositório, compêndio nos moldes do `wall`; vira o repositório `michelcias/wafc-studies` só no fim, por cópia sem histórico, com só o que reproduz o manuscrito e o supp (D60) |

Neste repositório ficam também:

- scripts curtos de conferência numérica de um resultado (`derivations/check/`);
- a cópia de referência de figuras e tabelas (`results/`), copiadas do
  compêndio, nunca geradas aqui;
- a documentação de reprodutibilidade que o artigo cita.

Regras do código em `wafc/`:

- **Toda função exportável tem teste** em `wafc/tests/`; o teste de
  recuperação exata (`θ*` na base, sem ruído, `λ → 0`) é o primeiro.
- **O `wall()` é referência de leitura, não de cópia.** O desenho por blocos,
  o reescalonamento para `[ε, 1 − ε]` e o descarte da função de escala são
  reimplementados aqui com nomes próprios (`wafc_design()`, `wafc_rescale()`),
  citando no comentário a função do `wall.R` que inspirou.
- **Dependências declaradas** em `wafc/R/load.R` e em `CONTINUAR.md` na
  mesma rodada em que entram.
- **Se um dia virar pacote** (próprio ou dentro do `WaveBased`), a decisão é
  tomada depois de E2.5, com o código já testado; a pasta é organizada desde
  já como `R/` + `tests/` para que a migração seja mover arquivos.

Quando o `results/` e o compêndio divergirem, o compêndio é a verdade.

---

## 7. Orquestrador e executores

A fase 2 ([`plano-fase2.md`](plano-fase2.md), D83) corre com um orquestrador
e executores, todos no máximo Opus, na **mesma árvore de trabalho**. O que
evita perda e gasto:

- **Um orquestrador por vez.** É o chat principal: o único que edita a
  cabeça do `ESTADO.md`, os arquivos de `docs/estado/`, o `plano-fase2.md`,
  o `notacao.md` e os cartões, e o único que commita (quando o autor pede).
  Uma sessão dele é um lote: integra os handoffs pendentes, compara cada
  medição com o critério declarado, escreve os cartões destravados,
  despacha, atualiza o quadro e termina. Ele não lê saída bruta de tarefa,
  só o handoff.
- **Uma tarefa é um cartão** (`docs/tarefas/<ID>.md`, até 60 linhas; o
  modelo está no `plano-fase2.md` §5): objetivo, modelo e esforço, a lista
  fechada do que ler (com intervalos de linha), a lista fechada do que
  escrever, os passos e o critério de saída verificável por comando. O
  executor lê o `CLAUDE.md`, o [`TAREFA.md`](TAREFA.md) e o cartão; não lê
  o estado.
- **Modelo e esforço são os do cartão.** Tarefa que cabe num turno e não
  pede aviso do autor vai a um subagente (os papéis estão em
  `.claude/agents/`); as outras abrem em chat próprio, com o modelo e o
  esforço escolhidos na abertura. Tarefa que falha o critério de saída duas
  vezes recebe cartão novo um nível acima, com a causa escrita.
- **Um executor só toca nos arquivos do cartão.** Não commita, não faz
  push. Símbolo novo, decisão ou mudança fora do escopo vira proposta no
  handoff, não edição.
- **Todo executor termina escrevendo `docs/handoff-<ID>.md`**, com até 60
  linhas e as quatro seções do `TAREFA.md` §3. O orquestrador integra (o
  bloco no diário do mês, a decisão em `estado/decisoes.md`, o quadro na
  cabeça), commita por caminho explícito quando o autor pedir e apaga o
  handoff e o cartão.
- **Commit por caminho explícito**, nunca `git add -A`, enquanto houver
  mais de uma tarefa aberta: o `git status` mostra o que é de quem.
- **Paralelo só com arquivos disjuntos.** Código novo vai em arquivo novo;
  o que é comum a duas tarefas corre em série. Etapa que depende de
  resultado ainda em curso não abre.
- **Rodada longa é script, não conversa:** PID gravado, marcador de fim, e
  o aviso do autor antes de começar. Quem lê o resultado é uma tarefa nova.
- **Verificação por script** sempre que existir (testes, compilação, junção
  exata, conferência numérica, texto contra tabela). Segundo modelo, em
  contexto novo, só para matemática e para as afirmações do resumo.

### Mensagens de abertura

Orquestrador:

> Leia `CLAUDE.md` e `docs/ESTADO.md`. Integre os handoffs pendentes e
> despache o próximo lote.

Executor, com o modelo e o esforço do cartão:

> Leia `docs/TAREFA.md` e execute o cartão `docs/tarefas/E9.1a.md`.

---

## 8. Continuidade entre máquinas

O trabalho troca de computador com frequência, e o que viaja é o que está
commitado. Estes arquivos são o contrato de continuidade, e mantê-los em dia
é parte de fechar uma etapa, não tarefa à parte:

| Arquivo | O que tem de estar em dia |
|---|---|
| [`ESTADO.md`](ESTADO.md) (a cabeça, até 150 linhas) | onde o trabalho está, o quadro da fase, as perguntas que pedem o autor, o próximo passo |
| `estado/diario-AAAA-MM.md` | um bloco por etapa fechada, com os números que a fecharam e a lição que custou horas |
| `estado/decisoes.md` | a decisão, numa linha da tabela; pergunta respondida vira linha aqui |
| `tarefas/<ID>.md` | um cartão por tarefa aberta, com os arquivos permitidos |
| [`CONTINUAR.md`](CONTINUAR.md) | a disposição das pastas, as ferramentas e **todo pacote novo que um script passou a usar**, mais os comandos de conferência |

As regras que vêm disso:

- **Etapa fechada é etapa registrada.** Não se encerra uma rodada com o
  resultado só no chat: o que não está no diário e no quadro não existe na
  máquina seguinte. O mesmo vale para uma lição que custou horas.
- **A cabeça do `ESTADO.md` não cresce.** O que deixa de ser corrente sai
  dela para o diário; ela diz onde o trabalho está, não como chegou lá.
- **Arquivo morto se consulta, não se lê:** `grep` pelo identificador (`D56`,
  `E6.1b`) e leitura do trecho.
- **Dependência nova entra no `CONTINUAR.md` na mesma rodada** em que o
  script passa a usá-la.
- **Edição nesses arquivos é verificada, não presumida.** Depois de
  editar, conferir que o trecho novo está lá.
- **Produto de compilação não é documentação.** PDF compilado viaja quando
  for a entrega; `.aux`, `.bbl` e afins ficam de fora até haver razão.
- **Antes de trocar de máquina**, `git status` limpo em `wafc-draft`; `HEAD`
  igual a `origin/main`; nenhum `docs/handoff-*.md` pendente.
