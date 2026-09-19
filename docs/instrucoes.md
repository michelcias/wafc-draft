# Instruções de trabalho

Documento normativo do projeto. Em caso de conflito com qualquer outro arquivo
deste repositório, este manda. Adaptado do `instrucoes.md` do `bdm-draft`, com
o que mudou: o alvo é uma revista de estatística geral (não bayesiana), a
teoria é de regressão penalizada (desigualdades oráculo, taxas em Besov), e o
código de produção entra num pacote que **já existe** (`WaveBased`).

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
2. Ler [`ESTADO.md`](ESTADO.md): onde o trabalho parou e o que já foi decidido,
   para não reabrir.
3. Ler a etapa em foco em `plano-projeto.md` (`ESTADO.md` diz qual é).
4. Se o foco for o manuscrito, ler o `.tex` da seção em questão **antes** de
   opinar sobre ela.
5. Se o foco for matemática, ler o arquivo correspondente em `derivations/` e
   [`notacao.md`](notacao.md) antes de escrever qualquer fórmula. A prova do
   WALL teórico (`../../wall-manuscript/manuscript/theo/ms_theo_1.tex`,
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
alvo for a *Statistica Sinica*: teto de **30 páginas de manuscrito** (a
confirmar na página oficial, que estava inacessível em 2026-09-18), provas
no suplementar, o leitor é metodólogo de estatística geral (a exposição de
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
| estudo de simulação e aplicações | `michelcias/wafc-studies`, compêndio nos moldes do `wall`, que fixa o código de `wafc/` por commit deste repositório |

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

## 7. Chats paralelos

Etapas sem dependência entre si (`plano-projeto.md`, tabela de dependências)
podem correr em chats separados na **mesma árvore de trabalho**. O que evita
perda:

- **Um chat principal por vez.** É o único que edita `ESTADO.md`,
  `plano-projeto.md` e `notacao.md`, e o único que commita. Os demais são
  chats de tarefa.
- **Um chat de tarefa faz uma etapa (ou subetapa) e só toca nos arquivos
  dela**, listados no catálogo de [`TAREFA.md`](TAREFA.md). Não commita, não
  faz push, não edita arquivo de outra etapa. Símbolo novo, decisão ou mudança
  fora do escopo vira proposta no handoff, não edição.
- **Todo chat de tarefa termina escrevendo `docs/handoff-<etapa>.md`** com
  quatro seções: o que foi feito (arquivos), o que a conferência numérica
  mostrou (números), o que deve entrar no `ESTADO.md` (decisões, sinais,
  lições), perguntas em aberto. O chat principal integra o handoff no
  `ESTADO.md`, commita por caminho explícito e apaga o arquivo de handoff.
- **Commit por caminho explícito**, nunca `git add -A`, enquanto houver mais
  de um chat aberto: o `git status` mostra o que é de quem.
- Etapa que depende de resultado ainda em curso não abre em paralelo.

### Mensagem de abertura de um chat de etapa

Uma linha, porque [`TAREFA.md`](TAREFA.md) carrega o resto:

> Leia `docs/TAREFA.md` e execute a tarefa E1.3.

---

## 8. Continuidade entre máquinas

O trabalho troca de computador com frequência, e o que viaja é o que está
commitado. Três arquivos são o contrato de continuidade, e mantê-los em dia
é parte de fechar uma etapa, não tarefa à parte:

| Arquivo | O que tem de estar em dia |
|---|---|
| [`ESTADO.md`](ESTADO.md) | um bloco na §2 por etapa fechada, com os números que a fecharam; a decisão na tabela; a pergunta na §4; o próximo passo na §5 |
| [`TAREFA.md`](TAREFA.md) | a tabela de estado da §2 e o catálogo da §3, com os arquivos permitidos de cada tarefa aberta |
| [`CONTINUAR.md`](CONTINUAR.md) | a disposição das pastas, as ferramentas e **todo pacote novo que um script passou a usar**, mais os comandos de conferência |

As regras que vêm disso:

- **Etapa fechada é etapa registrada.** Não se encerra uma rodada com o
  resultado só no chat: o que não está no `ESTADO.md` não existe na máquina
  seguinte. O mesmo vale para uma lição que custou horas.
- **Dependência nova entra no `CONTINUAR.md` na mesma rodada** em que o
  script passa a usá-la.
- **Edição nesses três arquivos é verificada, não presumida.** Depois de
  editar, conferir que o trecho novo está lá.
- **Produto de compilação não é documentação.** PDF compilado viaja quando
  for a entrega; `.aux`, `.bbl` e afins ficam de fora até haver razão.
- **Antes de trocar de máquina**, `git status` limpo em `wafc-draft`; `HEAD`
  igual a `origin/main`; nenhum `docs/handoff-*.md` pendente.
