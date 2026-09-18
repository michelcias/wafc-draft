# Estado do trabalho, handoff de continuidade

**Última atualização:** 2026-09-18.
**Etapa corrente:** E0 aberta (falta E0.3 pelo autor e a publicação no
GitHub); E1.1 é o próximo do chat principal. Nenhum chat de tarefa em curso.
Três decisões do autor pendentes (D4, D5, D8) e três perguntas na §4.
**Versão viva do manuscrito:** nenhuma (nasce em E5a como `k = 1`).
**Cor da rodada corrente:** `colR1` (entra em uso quando existir `k = 2`).

Este documento é o ponto de partida de cada sessão. Ele diz onde o trabalho
parou, o que já foi decidido (para não reabrir) e o que vem a seguir.

---

## 1. Para começar uma sessão

1. Ler `CLAUDE.md` (as duas regras de abertura: sem coautoria, respostas
   curtas).
2. Ler este arquivo.
3. Ler [`instrucoes.md`](instrucoes.md). É a regra, não a sugestão.
4. Ler a etapa em foco em [`plano-projeto.md`](plano-projeto.md). Em E1,
   também [`proposta-metodo.md`](proposta-metodo.md).
5. Se for escrever fórmula: [`notacao.md`](notacao.md) antes, e as Seções 2,
   3 e 5 de `../../wall-manuscript/manuscript/theo/ms_theo_1.tex`, que é o
   molde da prova.

### Mensagem de abertura sugerida

Chat principal (o que edita este arquivo e commita):

> Leia `CLAUDE.md`, `docs/ESTADO.md` e `docs/instrucoes.md` para se situar.
> Hoje vamos <foco do dia>.

Outra máquina (instalação e disposição das pastas em `CONTINUAR.md`):

> Leia `docs/CONTINUAR.md` e me diga o que falta instalar; depois leia
> `docs/ESTADO.md` e continuamos de onde parou.

Chat de tarefa (uma etapa, arquivos restritos, termina em handoff):

> Leia `docs/TAREFA.md` e execute a tarefa <etapa>.

### Git

O autor commita e faz push. O assistente lê o estado e lembra de commitar nos
momentos apropriados. **Sem coautoria nas mensagens**, sem exceção.

---

## 2. Onde o trabalho está

### 2026-09-18: o repositório nasceu

Sessão de avaliação de viabilidade e de planejamento. O autor descreveu a
ideia (coeficientes funcionais, cada um aditivo nas moduladoras, cada
componente em wavelets, estimação por LASSO); o assistente leu o `bdm-draft`
(molde das convenções), o `WaveBased` (`wall()`, `wbasis()`, `wtable()`), o
`wall` (compêndio) e o `wall-manuscript` (o teórico é o molde da prova), e
fez buscas de literatura. O que saiu:

- **Viável e com nicho.** A combinação "coeficientes aditivos + wavelets +
  LASSO" não apareceu na busca inicial; os três pilares existem em separado
  (Xue & Yang 2006; Zhou & You 2004; Sardy & Tseng 2004 e Sardy & Ma 2024;
  Wei, Huang & Li 2011). O trabalho mais próximo é Sardy & Ma (2024), sem o
  `X_j` multiplicando. Detalhe em [`proposta-metodo.md`](proposta-metodo.md)
  §5 e [`literatura.md`](literatura.md); a confirmação é L2.
- **O motor já existe.** A matriz de desenho do WAFC é a do `wall()` com cada
  bloco multiplicado por `X_j`; o ajuste é o mesmo `glmnet`. Inventário em
  [`inventario-codigo.md`](inventario-codigo.md). A arquitetura de prova do
  WALL teórico (sieve, Besov, oráculo, compatibilidade, compressibilidade)
  transfere; o que é novo na teoria é a condição de desenho para produtos
  (E1.4).
- **Revista.** Proposta: *Statistica Sinica* como alvo primário (linhagem do
  modelo, formato completo, 30 páginas), EJS como reserva (sem teto, open
  access, template já versionado e compilando em `manuscript/ejs-template/`).
  Comparação e requisitos em [`alvo-revista.md`](alvo-revista.md). O SCImago
  e o site da SS estavam inacessíveis desta máquina; a confirmação é do
  autor (E0.3).
- **Plano** em etapas E0 a E7 mais L1 e L2, com dependências e trilhas
  paralelas ([`plano-projeto.md`](plano-projeto.md)); catálogo de tarefas
  abertas em [`TAREFA.md`](TAREFA.md).
- **Ferramentas locais:** R 4.6.1, `glmnet`, `grpreg`, `gglasso`, `mgcv`,
  `latexmk`; o `WaveBased` 2.6-0 foi instalado de `../../WaveBased` e o
  `wall()` roda. Faltam `sparsegl` (E2.2) e `gh`.

### Decisões tomadas

| # | Data | Decisão | Razão |
|---|---|---|---|
| D1 | 09-18 | O modelo é `Y = Σ_j β_j(U) X_j + ε`, `β_j(u) = c_j + Σ_k g_{jk}(u_k)`, com `X_1 ≡ 1` permitido (o aditivo puro e o parcialmente linear aditivo são casos particulares) | pedido do autor |
| D2 | 09-18 | Base: wavelets ortonormais de suporte compacto do `WaveBased`; padrão periódico com `j0 = 0` e a função de escala constante descartada (identificabilidade no nível da base, como no `wall()`); `boundary = "interval"` como opção | herda o `wall()`; é o que faz a restrição `∫ g_{jk} = 0` sair de graça |
| D3 | 09-18 | Estimador base: LASSO sobre todos os coeficientes de wavelet, `c_j` não penalizados; sparse group LASSO por par `(j,k)` é a variante a medir em E2 | pedido do autor (LASSO); a variante em grupos é a candidata natural à seleção de estrutura |
| D6 | 09-18 | Documentos de trabalho em português; manuscrito em inglês americano; convenções de git, marcação e continuidade herdadas do `bdm-draft` | pedido do autor ("em linha com o bdm-draft") |
| D7 | 09-18 | Compêndio de simulação e aplicação em repositório próprio, `wafc-studies`, nos moldes do `wall` | o `wall` já resolveu cache, `renv` por commit e proveniência |

**A ratificar pelo autor (propostas do assistente):**

| # | Proposta | Razão | Alternativa |
|---|---|---|---|
| D4 | `wafc()` e `cv.wafc()` **dentro do `WaveBased`**, ao lado de `wall()`, reaproveitando as internas de desenho | 90 % do código já está lá; um pacote novo duplicaria; o `WaveBased` já é o pacote citado pelo grupo | pacote novo `wafc` que importa o `WaveBased` (mais limpo para citar, mais caro para manter) |
| D5 | Alvo primário **Statistica Sinica**; reserva EJS | `alvo-revista.md` §2 | EJS como primário se a teoria pesar mais que a aplicação |
| D8 | Nome do método **WAFC** (*wavelet additive functional coefficients*) | curto, ecoa o WALL, e a sigla é a do repositório | "WAVC" (varying coefficients); "wavelet additive coefficient LASSO" |

---

## 3. O que NÃO reabrir

- **Pasta dentro do `wafc-draft` ou repositório próprio para o código.**
  Respondido em D4 (a ratificar) e D7, detalhe em `plano-projeto.md` E0.2:
  protótipo aqui, função no `WaveBased`, compêndio próprio.
- **Se a restrição de identificabilidade precisa de multiplicador ou
  centralização numérica.** Não precisa no caso periódico com `j0 = 0`: a
  constante é a única função de escala e é descartada (D2). Só volta se
  `boundary = "interval"` virar o padrão.

---

## 4. Perguntas em aberto

Ordenadas pelo que bloqueia mais.

1. **D4, D5, D8** (tabela acima). D5 destrava E0.3 e a estrutura de E5a; D4
   destrava E3; D8 é nome e pode esperar E5a.
2. **Índice de nível das wavelets:** `l` (proposta, para deixar `j` à
   covariável linear) ou `j` como no WALL, mudando a covariável para outro
   índice. Decide-se em E1.1 (`notacao.md`, §6).
3. **Aplicação (E6.1):** o autor tem uma base em mente? Os candidatos de
   `plano-projeto.md` E6.1 são genéricos. Decidir cedo evita desenhar a
   simulação longe do caso real.
4. **Sparse group LASSO:** vale a dependência nova (`sparsegl`) no
   `WaveBased`, ou a variante em grupos fica só no protótipo e no artigo
   como comparação? Decide-se em E2.5 com número.
5. **`X` dependente de `U`:** a teoria de E1.4 tenta o caso geral ou o
   artigo assume `X ⊥ U` e discute o geral? Decide-se quando E1.4 mostrar o
   que fecha.

---

## 5. Próximos passos

Em ordem; (a) e (b) são independentes de (c) a (e).

- (a) **Autor:** ratificar D4, D5, D8; criar `michelcias/wafc-draft` no
  GitHub e fazer o push (`CONTINUAR.md`, §1); E0.3 (instruções e template
  da SS); responder a pergunta 3 se já tiver a base.
- (b) **Chat principal:** E1.1 (congelar a notação; `macros.tex`); depois
  catalogar E1.3, E1.4 e E2.1 com os arquivos permitidos definitivos.
- (c) **Chats de tarefa, em paralelo desde já:** L1 (verificação
  bibliográfica) e L2 (busca de novidade), que não dependem da notação.
- (d) **Chats de tarefa, depois de E1.1:** E1.2 e E1.3 (independentes);
  E2.1 (precisa só do enunciado de E1.2).
- (e) **Depois de E1.3 e E1.4:** E1.5, E1.6; E2.2 a E2.5.

## 6. Histórico de sessões

| Data | O que aconteceu |
|---|---|
| 2026-09-18 | Avaliação de viabilidade; criação do repositório e dos documentos de trabalho; template da EJS; plano E0 a E7 |
