# Estado do trabalho, handoff de continuidade

**Última atualização:** 2026-09-18.
**Etapa corrente:** E0 fechada, menos a ratificação de D5; E1.1, L1 e L2
fechadas. Liberadas para chats de tarefa: E1.2, E1.3 e E1.4 (E2.1 espera o
enunciado de E1.2). Nenhum chat de tarefa em curso. Duas decisões do autor
adiadas (D5, D8), seis propostas de L2 a ratificar e sete perguntas na §4.
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
- **O código fica aqui.** O autor decidiu (D4) que o `WaveBased` não recebe
  código por enquanto: o método vive na pasta `wafc/` (`R/`, `tests/`,
  `scripts/`), organizada como pacote sem ser pacote, e o `WaveBased`
  instalado é só dependência para as bases. Empacotar decide-se em E3.3.
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

### 2026-09-18: E0.3 fechada (instruções e template da Statistica Sinica)

O autor capturou as páginas oficiais da revista e baixou os templates; as
instruções estão transcritas em inglês em
[`ss-instrucoes-autores.md`](ss-instrucoes-autores.md) e o
[`alvo-revista.md`](alvo-revista.md) foi corrigido. Os números que fecham a
etapa:

- **Teto de 40 páginas em espaço duplo no template, referências e apêndice
  incluídos** (12pt, margens de 1 polegada), e não as 30 que a busca na web
  tinha dado. A estrutura-alvo soma ~29 páginas, com ~11 de folga não
  alocada.
- **Template obrigatório e conferido na recepção**: manuscrito fora do
  template volta antes de qualquer revisão. Versionado em
  `manuscript/ss-template/`: `SS-template.tex` (bibliografia manual),
  `SS-template-bib.tex` (BibTeX, `chicago`) e `supp-temp_20240820.tex`;
  os três compilam com o TeX Live 2023 local (avisos de `fancyheadings`
  obsoleto e de `\headheight`, sem erro).
- **Triagem:** ~40% das submissões passam dos co-editores; o AE manda a
  referee ~60% do que recebe. Alvo de prazo: 2 a 3 meses na primeira
  rodada.
- **Suplementar:** um único PDF de até 10 MB, mesmo título e mesmos autores,
  revisado junto; seção "Supplementary Material" como última seção do corpo,
  antes dos agradecimentos. A revista *prefere* provas, lemas técnicos e
  detalhes de simulação no suplementar, o que é o que o plano já fazia.
- **Reprodutibilidade é exigência editorial**, não cortesia: "software
  producing the evidence should be available for examination as well as
  pertinent datasets".
- **O que a página não diz**, listado na §7 de `ss-instrucoes-autores.md`: o
  **tipo de revisão** (nem cega simples nem dupla; sem exigência de
  anonimização), taxas, licença, política formal de dados e código, limite
  de resumo e de palavras-chave, e qual dos dois templates (Windows ou Mac)
  é o zip baixado. Um símbolo da seção de fórmulas veio como imagem quebrada
  na captura.

### 2026-09-18: E1.1 fechada (notação congelada)

O autor ratificou os pontos da §6 do [`notacao.md`](notacao.md), que sai do
estado de esboço; [`../derivations/macros.tex`](../derivations/macros.tex)
foi reescrito para implementá-lo e compila (documento de teste com o modelo,
a expansão do bloco e os ambientes). O que ficou:

- **A base fica com `j` e `k`.** `ψ_{jk}` mantém o padrão da literatura de
  wavelets; quem muda de letra são as covariáveis: linear `X_ℓ`,
  `ℓ = 1, …, p`; moduladora `U_m`, `m = 1, …, q`. Isso inverte a proposta do
  esboço, que punha a covariável em `j` e a wavelet em `ψ_{lm}`.
- **Correção de um erro do esboço:** o `notacao.md` dizia que o WALL usa `j`
  para o nível. Não usa: o WALL teórico escreve `ψ_{ℓ,k}` e depois achata o
  par num índice único `b_m` (`ms_theo_1.tex`, §2.2). A notação do WAFC
  mantém o par visível, porque o bloco `(ℓ, m)` é a unidade do desenho e o
  grupo da variante de E2.2.
- **Coeficientes em `θ`**, com `β_ℓ` reservado ao coeficiente funcional; o
  subscrito é `θ_{ℓm,jk}`, bloco antes da wavelet.
- **`U ∈ [0,1]^q` por hipótese populacional**; a transformação monótona da
  prática fica na seção de computação e fora da teoria de E1.3 a E1.6.
- **Colisão registrada:** `ℓ` é índice e é a letra da penalidade `ℓ_1`. A
  convenção que ficou é escrever a penalidade `‖θ‖_1` e nunca colar `ℓ_τ` a
  um índice (`notacao.md`, §1).
- **Ordem das colunas de `Z`** fixada: os `p` termos não penalizados, depois
  os blocos `(ℓ, m)` em ordem lexicográfica, e dentro do bloco `j` crescente
  e `k` crescente. É o que E2.1 tem de implementar.

### 2026-09-18: L1 fechada (verificação bibliográfica)

Chat de tarefa, integrado neste. `docs/referencias-verificadas.bib` tem
**35 entradas**, cada uma conferida no Crossref, e `literatura.md` não tem
mais `[VERIFICAR]` nem `[L1: confirmar]`. Conferido aqui: 35 `@`-entradas,
35 `\bibitem` com `bibtex` e `plain`, sem erro nem aviso; 32 com DOI, e as
três sem DOI são Xue & Yang (2006), Xue & Qu (2012, JMLR) e Haris, Simon &
Shojaie (2018, NeurIPS), que não têm registro.

O que a verificação corrigiu, e que muda o que o artigo cita:

- O arXiv 1903.04631 que `literatura.md` atribuía a Amato, Antoniadis, De
  Feis & Gijbels é de **Haris, Simon & Shojaie** (NeurIPS 2018): aditivos
  com wavelets em desenho irregular. O candidato ao que a linha queria
  citar é Amato et al. (2022, *Stat. Comput.*). As duas entradas ficaram no
  `.bib`, e **as duas entram na varredura de L2**; pelo resumo, nenhuma tem
  o `X_ℓ` multiplicando.
- Antoniadis & Gijbels (2014): o terceiro autor é **Lambert-Lacroix**, não
  Verhasselt.
- Klopp & Pensky é de **2015** (*Ann. Statist.* 43(3) 1273-1299).
- Sardy & Ma (2024) saiu em *Scand. J. Statist.* 51(1) 89-108
  (10.1111/sjos.12680); o WALL cita só o preprint. É o trabalho mais
  próximo, então a citação tem de ser a publicada.
- Donoho & Johnstone (1994 *Biometrika*; 1998 *Ann. Statist.*), Daubechies
  & Lagarias e Bickel, Ritov & Tsybakov estavam marcados `no wall` mas não
  estão no `references_theo_1.bib`; entraram com chave própria.

**Lição de ferramenta** (vale para L2 e para E5b): Crossref e OpenAlex
respondem desta máquina sem chave; o Semantic Scholar limita a ~1
consulta/s; o Project Euclid bloqueia `curl`; a *Statistica Sinica* serve os
PDFs antigos em `statistica/oldpdf/A{vol}n{iss}{art}.pdf`.

### 2026-09-18: L2 fechada (busca de novidade), e o vizinho é outro

Chat de tarefa, integrado neste. `docs/busca-novidade.md` traz as quatro
buscas, a varredura da *Statistica Sinica* e da EJS desde 2015 e o veredito
por contribuição; as linhas novas de `literatura.md` foram coladas aqui,
porque L1 tinha o arquivo na mão quando L2 rodou.

- **A novidade central se sustenta:** nenhum trabalho combina coeficientes
  aditivos em várias moduladoras, wavelets e LASSO. Na varredura por título,
  a interseção "coeficientes variáveis × wavelet" é **zero** na SS e **zero**
  na EJS desde 2015 (SS: 28 `varying coefficient`, 1 `additive coefficient`,
  2 `wavelet`; EJS: 7, 0, 6).
- **Mas o vizinho perigoso não é Sardy & Ma: é Klopp & Pensky (2015,
  *Ann. Statist.* 43(3) 1273-1299).** Para `q = 1` e `X ⊥ U`, eles já têm o
  desenho de produtos com Gram `Ω ⊗ Φ` (eq. 1.8), a concentração da Gram
  empírica restrita (Lema 1), o LASSO em blocos com desigualdade oráculo não
  assintótica, a taxa adaptativa em Besov com `ν < 2` e a cota inferior
  minimax. Ou seja: E1.4 (i), E1.5 e E1.6 estão publicados no caso de uma
  moduladora com desenho independente.
- **O que sobra para o WAFC**, e que passa a ser o eixo do artigo:
  aditividade em `q ≥ 2` moduladoras, com o termo cruzado
  `E[X_ℓ X_{ℓ'} ψ_{jk}(U_m) ψ_{j'k'}(U_{m'})]`, `m ≠ m'`, que não é produto
  de Kronecker; `X` dependente de `U`; LASSO puro com o corolário de
  compressibilidade; identificabilidade no nível da base com constantes não
  penalizadas; e o software, a simulação e a aplicação, que K&P não têm.
- **Sardy & Ma (2024) não é ameaça teórica:** os quatro teoremas deles são
  de otimização (degenerescência do group square-root LASSO em bloco
  ortonormal, limiarização suave em forma fechada, relaxação por blocos,
  SURE), sem teoria estatística e sem sieve.
- **Precursor do próprio grupo, obrigatório citar:** Montoril, Morettin &
  Chiann (2018, *IJWMIP* 16(1) 1850004), wavelets em coeficientes funcionais
  sem penalização.
- **Concorrentes que o referee vai cobrar em E4:** o spline adaptativo de
  Wang, Jiang & Liu (2024, *JCGS*), o block LASSO de K&P no mesmo desenho
  (sai com `grpreg`/`gglasso`) e um não aditivo em várias moduladoras
  (VCBART).
- **Lição de ferramenta:** o Crossref por ISSN com `from-pub-date` varre uma
  revista por título; o OpenAlex dá grafo de citações e resumos; a Wiley
  bloqueia leitura automatizada.

As linhas que L2 acrescentou ao `literatura.md` estão com status `resumo` ou
`[VERIFICAR]` e **não podem entrar no `.bib`** antes de uma verificação como
a de L1.

### Decisões tomadas

| # | Data | Decisão | Razão |
|---|---|---|---|
| D1 | 09-18 | O modelo é `Y = Σ_j β_j(U) X_j + ε`, `β_j(u) = c_j + Σ_k g_{jk}(u_k)`, com `X_1 ≡ 1` permitido (o aditivo puro e o parcialmente linear aditivo são casos particulares) | pedido do autor |
| D2 | 09-18 | Base: wavelets ortonormais de suporte compacto do `WaveBased`; padrão periódico com `j0 = 0` e a função de escala constante descartada (identificabilidade no nível da base, como no `wall()`); `boundary = "interval"` como opção | herda o `wall()`; é o que faz a restrição `∫ g_{jk} = 0` sair de graça |
| D3 | 09-18 | Estimador base: LASSO sobre todos os coeficientes de wavelet, `c_j` não penalizados; sparse group LASSO por par `(j,k)` é a variante a medir em E2 | pedido do autor (LASSO); a variante em grupos é a candidata natural à seleção de estrutura |
| D4 | 09-18 | **O código do método vive na pasta `wafc/` deste repositório** (`R/`, `tests/`, `scripts/`); o `WaveBased` não recebe código por enquanto e é usado só como dependência para as bases (`wbasis()`, `wtable()`), confirmado pelo autor; **as funções criadas ficam privadas por enquanto** (neste repositório privado, sem pacote público, sem `install_github`); empacotar e publicar decide-se em E3.3, com o código testado | decisão do autor, contra a proposta de implementar dentro do `WaveBased` |
| D9 | 09-18 | **Índices:** wavelet `ψ_{jk}` (nível `j`, translação `k`); covariável linear `X_ℓ`, `ℓ = 1, …, p`; moduladora `U_m`, `m = 1, …, q` | decisão do autor: não mexer no padrão da base; a colisão sai das covariáveis |
| D10 | 09-18 | **Coeficientes em `θ`**, com `θ_{ℓm,jk}` (bloco antes da wavelet); `β_ℓ` fica sendo só o coeficiente funcional | `β` não pode ser função e vetor na mesma seção |
| D11 | 09-18 | **`U ∈ [0,1]^q` por hipótese populacional**, com densidade limitada longe de `0` e de `∞`; reescalonamento empírico só na seção de computação | teoria limpa em E1.3 a E1.6; o termo extra não vale o custo agora |
| D12 | 09-18 | Ordem das colunas de `Z`: não penalizados, depois blocos `(ℓ, m)` lexicográficos, dentro do bloco `j` e `k` crescentes | fixa a interface de `wafc_design()` em E2.1 |
| D6 | 09-18 | Documentos de trabalho em português; manuscrito em inglês americano; convenções de git, marcação e continuidade herdadas do `bdm-draft` | pedido do autor ("em linha com o bdm-draft") |
| D7 | 09-18 | Compêndio de simulação e aplicação em repositório próprio, `wafc-studies`, nos moldes do `wall` | o `wall` já resolveu cache, `renv` por commit e proveniência |

**Propostas de L2, a ratificar (2026-09-18).** Nenhuma bloqueia E1.2, E1.3 ou
E1.4; todas mudam documento, não resultado.

| # | Proposta | Onde |
|---|---|---|
| L2a | Frase-tese e contribuição 1 reescritas como extensão de Klopp & Pensky a coeficientes aditivos em várias moduladoras e desenho dependente; "Why not block LASSO?" entra nas perguntas do referee | `alvo-revista.md` §4 |
| L2b | E1.4 (i) vira "recordar K&P (eq. 1.8, Lema 1) e estender ao desenho aditivo"; o entregável central passa a ser o termo cruzado entre moduladoras e a parte (ii) | `plano-projeto.md` E1.4 |
| L2c | O QUT (Giacobino et al. 2017), que Sardy & Ma usam, entra em E2.3 como regra de `λ` sem `σ`, ao lado de BIC/EBIC | `plano-projeto.md` E2.3 |
| L2d | Concorrentes mínimos de E2.4/E4: `mgcv`, spline adaptativo (Wang, Jiang & Liu 2024), block LASSO de K&P no mesmo desenho, VCBART; cenário não homogêneo com as funções de Donoho-Johnstone | `plano-projeto.md` E2.4, E4 |
| L2e | "O mais próximo na teoria é Klopp & Pensky; no método, Sardy & Ma e Amato et al."; citar Montoril, Morettin & Chiann (2018) | `proposta-metodo.md` §5 |
| L2f | E1.7 fica opcional ou vira variante sem teorema de seleção: a ideia já existe em splines (Antoniadis et al. 2014; Ma et al. 2015) e em wavelets sem `X_ℓ` (Amato et al. 2022) | `plano-projeto.md` E1.7 |

**Adiadas pelo autor ("cobre-me depois", 2026-09-18); propostas do assistente:**

| # | Proposta | Razão | Alternativa |
|---|---|---|---|
| D5 | Alvo primário **Statistica Sinica**; reserva EJS | `alvo-revista.md` §2 | EJS como primário se a teoria pesar mais que a aplicação |
| D8 | Nome do método **WAFC** (*wavelet additive functional coefficients*) | curto, ecoa o WALL, e a sigla é a do repositório | "WAVC" (varying coefficients); "wavelet additive coefficient LASSO" |

---

## 3. O que NÃO reabrir

- **Pasta dentro do `wafc-draft`, `WaveBased` ou repositório próprio para o
  código.** Respondido em D4 e D7, detalhe em `plano-projeto.md` E0.2: o
  método em `wafc/` aqui, compêndio próprio. O empacotamento volta como
  decisão só em E3.3, com o código testado.
- **Se a restrição de identificabilidade precisa de multiplicador ou
  centralização numérica.** Não precisa no caso periódico com `j0 = 0`: a
  constante é a única função de escala e é descartada (D2). Só volta se
  `boundary = "interval"` virar o padrão.

---

## 4. Perguntas em aberto

Ordenadas pelo que bloqueia mais.

1. **D5 e D8** (tabela acima), adiadas. E1.1 fechou, que era o marco para
   cobrá-las: D5 fixa o formato de E5a e D8 é o nome que vai no título.
2. **Aplicação (E6.1):** o autor tem uma base em mente? Os candidatos de
   `plano-projeto.md` E6.1 são genéricos. Decidir cedo evita desenhar a
   simulação longe do caso real.
3. **Sparse group LASSO:** a variante em grupos entra em `wafc()` como
   opção (`penalty = "sglasso"`, dependência `sparsegl`) ou fica só no
   piloto e no artigo como comparação? Decide-se em E2.5 com número.
4. **`X` dependente de `U`:** a teoria de E1.4 tenta o caso geral ou o
   artigo assume `X ⊥ U` e discute o geral? Decide-se quando E1.4 mostrar o
   que fecha.
5. **Como o artigo se posiciona diante de Klopp & Pensky** (a mais
   importante desta rodada): o WAFC se apresenta como extensão deles
   (honesto, e o referee da SS reconhece) ou como modelo diferente, com K&P
   citado como caso particular? A primeira muda a frase-tese de
   `alvo-revista.md` §4; a segunda obriga a Seção 2 a mostrar que o desenho
   aditivo não é o deles. Decide-se antes de E5a; E1.4 não depende.
6. **Block LASSO de K&P:** entra em `wafc()` como opção de penalidade, ao
   lado do sparse group LASSO de E2.2, ou fica só como concorrente em E2/E4?
7. **Bibliografia, quatro pontos deixados por L1** (nenhum bloqueia; o
   `.bib` fica como está até a resposta):
   - Amato et al. (2022) ou Haris, Simon & Shojaie (2018) na linha que
     citava o arXiv 1903.04631? Proposta: **as duas**, que são trabalhos
     distintos e ambos interessam a L2.
   - Hastie & Tibshirani (1993): páginas 757-779 (só o artigo) ou 757-796
     (com a discussão)? Proposta: **757-779**, que é o que o Crossref
     registra.
   - Daubechies & Lagarias: as duas partes ou só a parte I, que é a do
     algoritmo do `WaveBased`? Proposta: **só a parte I**, e a parte II
     entra se alguma prova precisar dela.
   - Chaves do `.bib`: padrão único `Autor-Autor-Ano` ou manter as quatro
     herdadas do WALL (`cohen1993wavelets` etc.)? Proposta: **manter as do
     WALL**, porque o `ms_theo_1.tex` é o molde da prova e a citação
     cruzada fica direta.

---

## 5. Próximos passos

(a) e (b) correm em paralelo; (c) espera E1.2.

- (a) **Autor:** ratificar D5 e D8 (a notação já não depende disso);
  responder a pergunta 2 se já tiver a base da aplicação.
- (b) **Chats de tarefa, em paralelo desde já:** E1.2, E1.3 e E1.4, que
  tocam arquivos distintos. E1.4 abre com o achado de L2: a parte (i) é
  citação de K&P, e o entregável é o termo cruzado entre moduladoras e a
  parte (ii).
- (b') **A catalogar quando o autor ratificar as propostas de L2:** a
  verificação bibliográfica das linhas novas de `literatura.md` (uma L1 de
  segunda rodada) e as edições de `alvo-revista.md` e `plano-projeto.md`
  listadas na §6 de `busca-novidade.md`.
- (c) **Chat de tarefa, depois do enunciado de E1.2:** E2.1, que precisa
  saber o que a identificabilidade descarta do desenho.
- (d) **Chat principal:** integrar os handoffs e catalogar E1.5, E1.6 e
  E2.2 a E2.5 quando as dependências fecharem.
- (e) **Depois de E1.3 e E1.4:** E1.5, E1.6; E2.2 a E2.5.

## 6. Histórico de sessões

| Data | O que aconteceu |
|---|---|
| 2026-09-18 | Avaliação de viabilidade; criação do repositório e dos documentos de trabalho; template da EJS; plano E0 a E7 |
| 2026-09-18 | D4 decidida pelo autor (código em `wafc/`, não no `WaveBased`); D5 e D8 adiadas; `prototype/` virou `wafc/`; plano E2 e E3 reescritos; repositório publicado; o autor confirmou o `WaveBased` como dependência e que as funções ficam privadas |
| 2026-09-18 | L2 fechada em chat de tarefa e integrada: novidade central confirmada (interseção zero na SS e na EJS), mas Klopp & Pensky (2015) cobre E1.4 (i), E1.5 e E1.6 para `q = 1` e `X ⊥ U`; seis propostas de mudança de rumo a ratificar |
| 2026-09-18 | L1 fechada em chat de tarefa e integrada: 35 referências verificadas, quatro correções de atribuição, dois trabalhos novos para L2 olhar |
| 2026-09-18 | E1.1 fechada: notação congelada (`ψ_{jk}`, `X_ℓ`, `U_m`, `θ_{ℓm,jk}`, `U ∈ [0,1]^q`), D9 a D12, `macros.tex` reescrito e compilando |
| 2026-09-18 | Máquina nova conferida (R 4.6.1, `WaveBased` 2.6-0, `grpreg`, `gglasso`, `sparsegl` por `apt`); E0.3 fechada: instruções da SS transcritas, templates versionados e compilando, teto corrigido de 30 para 40 páginas |
