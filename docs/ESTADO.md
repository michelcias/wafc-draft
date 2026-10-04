# Estado do trabalho, handoff de continuidade

**Última atualização:** 2026-10-03.
**Etapa corrente:** **E0, E1 (com E1.3b, E1.4c, E1.7a, E1.7c, E1.8,
E1.10 e E1.11), E5a, E5c e E2 inteira fechadas**, mais L1 a L9 e a
sondagem E6.1a. **E2.5 fechou em 2026-10-03 com a decisão de rumo do
autor (D44 a D46): go reposicionado na *Statistica Sinica*, o block LASSO
balanceado seguido do limiar `cv1se` é o WAFC, o LASSO fica como opção, e o
spline de E4 entra com `gam.reml` e `gam.gcv`.** O caminho até ela: o LASSO puro levou
no-go pelo critério literal em E2.5a; E2.5b a E2.5f mediram cinco formas
do block LASSO (a balanceada é a que a teoria cobre e a melhor fora do
suave); E2.5g mostrou que o limiar do Corolário 8 derruba o fator do suave
abaixo de 1,5 contra o `gam.matched`; e **E2.5h, com o `gam` sintonizado
nos termos dele (REML, validação cruzada nas dobras do WAFC e GCV), manteve
o veredito**: o WAFC perde no `smooth` e no `uneven` em toda forma, o
`klopp.balanced+cv` vence no não homogêneo, na `mixed` e no nulo, e o fator
do `wafc.lasso+max` no suave fica no limite de 1,5 (1,49 contra o `gam.cv`
em `n = 1000`). **Em 2026-10-03, relido no chat principal:** o limiar
oráculo dá o teto de qualquer regra de `t`, e só o block LASSO balanceado
passa a perna do não homogêneo com ele; contra o `gam.gcv` (o critério
padrão do `mgcv::gam()`) essa perna some. **E2.5j (2026-10-03):** a regra `cv1se` (o maior
`t` a um erro-padrão do mínimo) chega ao teto no suave, no `uneven`, no
nulo e em `n = 1000`, e paga em `n ≤ 500` no não homogêneo e na `mixed`;
o `c` relativo repete o `+cv`; a porta do QUT resolve o nulo das outras
regras com nível acima do nominal; os cortes de caminho dos dois motores
não tocam o `λ` escolhido. Por D43, nada medido em E2.5 é descartado. **E1.12
fechou (2026-10-03):** a teoria em blocos está provada em `08-blocos.tex`
(Corolário 11 sem logaritmo em `π ≥ 2`) e no adendo do `06` (o risco do
limiarizado); a pergunta 42 foi decidida (D47). **E3.1 fechou
(2026-10-03):** `cv.wafc(x, u, y)` é o WAFC de D44 e D45, a interface espera
a ratificação do autor, dada em D48. **E1.13 e E3.2 fecharam (2026-10-03):**
a taxa lenta em blocos (Corolário 14, que em `π ≥ 2` e `s < 1/2` dá
`n^{−2s/(2s+1)}` sem logaritmo e sem condição de desenho) e os gráficos com
a documentação; pendências decididas em D49. E6.1b está pronta e espera a
ordem do autor para rodar (§5). A decidir: a saída da aplicação (pergunta 2). O
teto de páginas fica para o fim (D21).
**Versão viva do manuscrito:** `k = 3` (`manuscript/ms_3.tex`,
`supp_3.tex`, `references_3.bib`), aberta em 2026-10-01 por decisão do
autor para a notação da seleção (D40); a `k = 2` fica intacta. Os dois
compilam limpos em **34 e 30 páginas** (34 desde a pergunta 40).
**Cor da rodada corrente:** `colR1`, em uso desde `k = 2` (a rodada não foi
aceita, então a marcação de E5c segue em `colR1` na `k = 3`).

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
  `X_ℓ` multiplicando. Detalhe em [`proposta-metodo.md`](proposta-metodo.md)
  §5 e [`literatura.md`](literatura.md); a confirmação é L2.
- **O motor já existe.** A matriz de desenho do WAFC é a do `wall()` com cada
  bloco multiplicado por `X_ℓ`; o ajuste é o mesmo `glmnet`. Inventário em
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

### 2026-09-18: E1.2, E1.3 e E1.4 fechadas (o núcleo da teoria)

Três chats de tarefa, integrados neste. Rodei as três conferências aqui
antes de integrar, e as três imprimem `OK` e reproduzem os números dos
handoffs: `01-identificabilidade.R` (7 s), `03-desenho-produtos.R` (9 s),
`02-aproximacao-besov.R` (53 s).

**E1.2, identificabilidade** (`derivations/01-identificabilidade.md`).
Proposição 1 (das `β_ℓ` a partir de `f` sob `E[XX' | U]` não singular q.c.;
de `(c_ℓ, g_{ℓm})` a partir de `β_ℓ` sob densidade positiva e centralização)
e Lema 1 (a base periódica com `j_0 = 0` impõe a centralização no nível da
base; `(c, θ)` é injetiva e `Σ` definida positiva a cada `J`). Números:
`Z` com 30 colunas tem posto 30 em `n = 200`, `λ_min(Z'Z/n) = 7.8e-03`,
recuperação sem ruído com erro `3e-15`. As deficiências previstas
aparecem exatas: manter `φ_{00}` dá nulidade `4 = p q`; `U_2 = U_1` dá
`14 = p N_J`; `X_2 = 1 + ψ_{10}(U_1)` dá `1`.

**E1.3, aproximação em Besov** (`derivations/02-aproximacao-besov.tex`,
7 páginas). Lema de aproximação em `L_2` e em `L_∞`, corolário para
`f - f_J`, e uma proposição sobre o custo da periodização. Números: em
`sin(2πu)` a queda por nível é `4.00` bits (isto é `2^{−JN}`, `N = 4`); em
`bumps` a queda média em `j = 9..12` é `1.37`, contra `s' = 3/2`; numa
função que não emenda em `0 ≡ 1`, a queda trava em `0.50` bits e o erro
uniforme **não decai** (`0.54` e `0.92` em `J = 9`).

**E1.4, Gram do desenho de produtos**
(`derivations/03-desenho-produtos.tex`, 6 páginas), que entregou mais do que
o catálogo pedia: no caso geral sai o **autovalor mínimo cheio**, não só a
compatibilidade no cone. `λ_min(Σ) ≥ κ_1 c_U` e `λ_max(Σ) ≤ κ_2 C_U`.
Números: sob `X ⊥ U` a fatoração `Π(E[XX'] ⊗ Σ_Ψ)Π'` bate a `7.6e-13`; com
`X` dependente de `U`, `λ_min(Σ) = 0.2747` contra a cota `κ_1 c_U = 0.2503`
(10% de folga); com `κ_1 = 0`, a Gram degenera com `J`, como previsto
(`4e-04`, `<1e-04`, `<1e-04` em `J = 2, 3, 4`).

O que isso muda no rumo:

- **A pergunta 4 da §4 está respondida: o caso geral fecha** (vira D13). A
  independência compra só a fatoração exata, não constante melhor.
- **A linha de risco "a condição de desenho dos produtos não fecha em
  geral" sai** da tabela de `proposta-metodo.md` §6, com número.
- **A taxa `2^{−Js'}` só aparece quando `2^{−J}` fica abaixo da escala da
  função** (em `bumps`, `J ≥ 9`). Os `n` de E4 (`≤ 2000`) com
  `J_n ≍ log_2 n/(2s'+1)` não chegam lá: a simulação vive no regime
  pré-assintótico, que é justamente onde o sieve linear perde e o LASSO
  ganha. Isso é matéria da Seção 5 do artigo, não defeito da conferência.
- **Periodização custa `1/2`:** para qualquer `g` que não emende, a taxa em
  `L_2` trava em `2^{−J/2}` e o erro uniforme não converge. Daí a pergunta 6
  da §4.
- **Em `n` pequeno a compatibilidade é mais saudável que o autovalor cheio:**
  `φ²(S)` fica em 65-75% do populacional em `n = 200`, contra 40-50% de
  `λ_min(Σ̂)`. É `φ²(S)` que E1.5 consome.
- **Achado sobre o `wall()`** (não tocado, D4): com `boundary = "interval"`
  ele mantém as `2^{j_0}` funções de escala de cada covariável e ainda o
  intercepto do `glmnet`, logo a constante entra `d + 1` vezes e a Gram fica
  singular com nulidade `d`. O ajuste não sofre; a leitura dos coeficientes
  de escala, sim. Fica registrado para o autor decidir se avisa lá.

**Numeração global dos resultados** (TAREFA.md §5 tem a tabela): E1.2 fica
com Proposição 1 e Lema 1; E1.3 com Lema 2, Lema 3, Corolário 1 e
Proposição 2; E1.4 com Proposição 3. Os arquivos mantêm os contadores
locais; o mapa é a autoridade quando E5a montar o manuscrito.

### 2026-09-19: E1.5 e E2.1 fechadas (o oráculo e o primeiro código)

Dois chats de tarefa, integrados aqui. Conferências rodadas nesta máquina
antes de integrar: `04-oraculo.R` imprime `OK` em 16 s;
`testthat::test_dir("wafc/tests")` dá **103 passam, 0 falham** em 4,6 s;
`wafc/scripts/01-smoke.R` imprime `OK` em 2,8 s; o `.pdf` de E1.5 tem
9 páginas.

**E1.5, desigualdade oráculo** (`derivations/04-oraculo.tex`). Lemas 4 a 7,
Teorema 1 e Corolários 2 e 3. A estrutura da prova evita a condição de cone
por **perfilagem**: com `B̃ = (I − Π_A)B`, o estimador resolve o problema
residualizado e `λ_min(B̃'B̃/n) ≥ λ_min(Σ̂)`, que é a forma consumível de
E1.4. Números: com `s_0 = 4` e sem viés, a razão `‖f̂ − f‖²_n/(λ² s_0)` fica
entre `1.541` e `1.617` ao variar `n` de `200` a `6400` (máximo sobre
mínimo `1.05`), que é exatamente a escala que o plano mandava medir; nenhuma
violação da cota rápida em 360 ajustes, com folga de `371×` a `187×`.

- **D13 pagou aqui:** é o autovalor cheio de E1.4 que dispensa o cone e dá
  constante que **não depende de `S_0`**.
- **Um termo que Klopp & Pensky não têm: `σ² p / n`**, preço de deixar os
  `c_ℓ` fora da penalidade (D3), já que K&P penalizam a constante. A
  diferença para eles não é só de generalidade, aparece no enunciado.
- **Calibração de `λ`:** a leitura certa de `‖Z‖_max` no plano é
  `σ̂_max = max_a sqrt(Σ̂_aa)`, que é `O_p(1)`; a leitura literal
  `max_{i,a}|Z_{ia}|` é de ordem `2^{J/2}` e **destruiria a taxa de E1.6**.
  Medido: `σ̂_max` vai de `1.02` a `1.14` em `J = 2..6` enquanto
  `max|Z_{ia}|` vai de `3.04` a `11.47`.
- **O Lema 7 dá `‖θ*‖_1` sem dependência de `π`**, o que permite enunciar a
  taxa lenta em `B^s_{π,r}` geral, e não só em `B^s_{2,∞}` como o WALL.
- **Sob Besov, em geral `s_0 = d`:** o enunciado incondicional é a taxa
  lenta, e trocar `s_0` por dimensão efetiva é o corolário de
  compressibilidade de E1.6.

**E2.1, desenho e cenários** (`wafc/R/load.R`, `dgp.R`, `design.R`,
`tests/test-design.R`, `scripts/01-smoke.R`). `wafc_design()` implementa a
ordem de D12 e descarta só o que E1.2 manda (a coluna `φ_{00}` de cada
bloco); a ordem bate com a construção de referência de
`check/03-desenho-produtos.R` a `1e-12`; mínimos quadrados sem ruído com
erro `1.1e-15` e `glmnet` com `λ → 0` a `4.3e-06`. No piloto de fumaça, o
cenário nulo zera as 90 colunas penalizadas (RMSE `0.045` contra
`σ = 0.62`): o desenho não inventa estrutura.

- **Achado que fixa E2.2:** o `glmnet` **descarta colunas de variância
  zero** mesmo com `standardize = FALSE`, logo `X_1 ≡ 1` sai do ajuste. Com
  `intercept = FALSE` perde-se `c_1` inteiro; com `intercept = TRUE` o
  intercepto carrega `c_1` exatamente (`1.00000025` contra `1`). `wafc()`
  ajusta com `intercept = TRUE` e soma o intercepto ao nível da covariável
  constante; `wafc_design()` já devolve `constant` com esses índices. E1.5
  chegou ao mesmo achado por outro caminho e acrescentou uma terceira rota
  conferida (resolver o problema residualizado), além de duas notas: o
  `glmnet` **reescala `penalty.factor` para somar `nvars`**, e `thresh` está
  obsoleto na versão 5.0 (vai em `control`).
- **`lambda.min` deixa lixo nos blocos nulos** (3 e 8 coeficientes nos dois
  blocos de `β_3`, ISE `0.0010` e `0.0080`), e `lambda.1se` reduz sem zerar.
  É o argumento numérico a favor da variante em grupos de E2.2.
- **O reescalonamento herdado do `wall()`** (`eps = 1.9^{−J}`, isto é
  `[0.077, 0.923]` em `J = 4`) faz a componente estimada integrar zero na
  amplitude reescalada, não em `[0,1]`. É deslocamento de nível, não de
  forma, e é matéria da seção de computação.
- **Lição de ferramenta:** `Matrix::cbind2` só aceita dois argumentos, e
  `do.call(Matrix::cbind2, lista)` devolve em silêncio o `cbind` dos dois
  primeiros; usar `Reduce`.

**Numeração global** (mapa na §5 do `TAREFA.md`): E1.5 ficou com Lemas 4 a
7, Teorema 1 e Corolários 2 e 3, e **passou a imprimir o número global** via
`\setcounter`. Adotei a prática; `02` e `03` continuam com contador local e
o mapa é a autoridade até alguém alinhá-los.

### 2026-09-19: E1.6, E2.2 e L3 fechadas; E1 inteira fechada

Três chats de tarefa, integrados aqui. Conferido nesta máquina antes de
integrar: `05-taxas.R` imprime `OK` (24 verificações); a bateria de `wafc/`
dá **231 passam, 0 falham** em 14 s; `01-smoke.R` imprime `OK`; o `.bib` com
57 entradas compila com `bibtex` sem erro nem aviso.

**E1.6, taxas e dimensão efetiva** (`derivations/05-taxas.tex`, 10 páginas).
Lemas 8 e 9, Proposição 4, Teorema 2, Corolários 4 e 5. O eixo é um achado
que muda o enunciado principal do artigo: **a hipótese de Besov de E1.3 já
contém a compressibilidade** (Lema 9(i): weak-`ℓ_τ` com
`τ = (s + 1/2)^{−1}`). O corolário de compressibilidade não é hipótese
extra, é leitura mais fina da mesma hipótese, e é ele que recupera a taxa de
aproximação **não linear** `n^{−2s/(2s+1)}` onde o sieve linear só chega a
`n^{−2s'/(2s'+1)}`; as duas coincidem em `π ≥ 2` e se separam em `π < 2`,
que é a diferença que E1.3 tinha registrado.

- **A perda quadrática dispensa a localização do WALL.** A trilha rápida de
  lá precisa de `(M*)² J 2^J = o(n)` (o que restringe `τ < 2/3`) e de um
  piso em `s`, porque localiza o estimador não modificado. Aqui o Lema 8 é a
  prova de E1.5 relida com outro comparador: `τ ∈ (0,2)` inteiro e piso
  `s' > (2−τ)/4`, que o WALL só obtém para o estimador restrito. É vantagem
  técnica citável na comparação com o companheiro.
- Números: a razão do erro de predição ao alvo `(J_n/n)^{2s'/(2s'+1)}` vai
  de `2.95` a `3.70` com `n` de 200 a 12800 (máx/mín `1.25`; nas
  componentes, `1.72`); a dimensão efetiva `M*` fica de `51×` a `103×`
  abaixo da leitura densa `M = d`.
- **O pré-assintótico, medido de novo:** com `s' = 1/4` a regra leva `J_n` a
  7 só em `n = 12800`, e é só aí que a condição empírica de E1.4 vale em
  95% das réplicas (contra ≤10% em `n ≤ 6400`); a regra fica **dois níveis
  acima** do `J` que minimiza o erro realizado, ao custo de 5% a 11%. A
  Seção 5 do artigo tem de dizer isso.
- **Não há cota inferior**, e construir a de `q ≥ 2` com desenho aditivo é
  trabalho do porte de E1.4 mais E1.5, fora do plano. Ver pergunta 12.

**E2.2, o estimador** (`wafc/R/fit.R`, `reconstruct.R`, `tests/test-fit.R`).
`wafc()` com LASSO e sparse group LASSO (grupo = bloco `(ℓ,m)`), `coef`,
`predict`, `wafc_functions()`, `wafc_blocks()` e `wafc_kkt()`. As KKT fecham
em `74/74` pontos do caminho no LASSO e `100/100` nos grupos.

- **A escala de `λ` ficou fixada por KKT:** o `λ` do objeto é o do objetivo
  de E1.5, e o do `glmnet` é esse dividido por `nvars/npen` (medido:
  `1.016129` contra razão observada `1.016136`). O `sparsegl` não reescala.
- **Resposta com número à pergunta 3 da §4** (a variante em grupos entra em
  `wafc()` ou fica no piloto): **entra**. No mesmo desenho e com ~99 não
  nulos, o LASSO zera `0` dos 6 blocos e o sparse group LASSO zera
  exatamente os `3` blocos nulos verdadeiros; o preço é predição (`RMSE`
  `0.2633` contra `0.2078`, com `σ = 0.7314`) e tempo (`3,7 s` contra
  `0,10 s`). A escolha entre as duas continua sendo de E2.5.
- **Correção do achado de E2.1:** o `glmnet` só descarta a coluna constante
  quando o desenho é **denso**; no esparso ele a mantém e reparte o nível
  com o próprio intercepto. A rota "somar o intercepto ao nível" é a única
  que funciona nas duas formas, e bate a `1e-13` entre elas.
- Dois defeitos de motor contornados: o `sparsegl` devolve o vetor nulo no
  primeiro ponto do caminho (inclusive os `c_ℓ`, que não é solução), e
  `asparse = 1` não é jeito confiável de obter o LASSO por ele (violações de
  2% a 5% de `λ`).

**L3, segunda rodada bibliográfica.** `referencias-verificadas.bib` vai de
35 a **57 entradas**; `literatura.md` não tem mais nenhuma linha em `resumo`
nem `[VERIFICAR]`; os `derivations/` não têm mais marca bibliográfica
pendente. O que ela corrigiu, e que teria ido para o manuscrito errado:

- **Três anos e um autor de L2 estavam errados:** Deshpande et al. é **2026**
  (*Bayesian Anal.* 21(1) 281-308), não 2024; Zhou, Xu & Lin é **2017**, não
  2016; e o terceiro autor de Zhou, Yang & Xiang (2022) é **Xiang**. A
  origem dos dois erros de ano é o prefixo do DOI (`10.1214/24-…` é o ano de
  **aceitação**), e a lição é conferir `published-print` no Crossref, não o
  campo `issued`.
- **Duas citações de teorema do `04-oraculo.tex` estavam trocadas:** em
  Bühlmann & van de Geer (2011) a taxa lenta é o **Corolário 6.1**, não o
  Teorema 6.1, e as coordenadas não penalizadas são a **§6.9** ("The
  weighted Lasso"), não a §6.2.3. A matemática de E1.5 não muda, os rótulos
  sim. O Corolário 6.8 que E1.4 cita está certo. O mapa completo ficou em
  comentário na entrada `buhlmann2011statistics`, para não se refazer a
  busca em E5a.
- Sobra uma única marca `verificar` nos `derivations/`, e é **matemática**,
  não bibliográfica: o esboço da prova da extensão na Proposição 2 de E1.3.

### 2026-09-19: E2.3 e E5a fechadas; o manuscrito existe

Dois chats de tarefa, integrados aqui. Conferido nesta máquina: a bateria de
`wafc/` dá **402 passam, 0 falham** em 25 s; `ms_1.tex` e `supp_1.tex`
compilam com `latexmk` **sem nenhuma referência ou citação indefinida**, em
28 e 26 páginas.

**E5a, o manuscrito nasce** (`manuscript/ms_1.tex`, `supp_1.tex`,
`references_1.bib` com 30 das 57 entradas verificadas). As Seções 1 a 4
estão escritas no template da SS, a Introduction usa o parágrafo de D18 como
está, e a Seção 3 abre pelo Corolário 5 de E1.6, que no manuscrito é o
`Theorem 1` (D16). O mapa completo da numeração global para os rótulos do
LaTeX está no `ms_1.tex` e é o que E5b precisa.

- **O teto ficou para o fim (D21).** Corpo em 24 páginas mais 4 de
  referências, contra as 14,5 que `alvo-revista.md` §5 previa para as
  Seções 1 a 4. Com as 11 planejadas para E5b e as referências crescendo
  para ~6, o total projetado é **39 a 41 páginas**, isto é, no teto ou
  acima. As três saídas estão na §2 do handoff, e a barata é escrever as
  Seções 5 a 7 em menos de 11 páginas, com tabelas ao suplemento, que é o
  que a revista prefere.
- **Sete decisões de redação** foram tomadas para o texto existir, todas
  reversíveis com uma linha, e todas resolvidas na varredura de 2026-09-19
  (D22, D24 e o checklist de submissão): centralização de Lebesgue
  na equação (2.2) (é a única que mexe em hipótese); teoria na base
  periodizada com `interval` como opção de computação; reescalonamento
  descrito sem fixar padrão; `E(·)` e `P(·)` em romano com parênteses e um
  único `B_X` no lugar de `C_X` e `B_X`, exigidos pela §5 das instruções da
  revista; macros copiadas para dentro dos `.tex` em vez de `\input`;
  `\numberwithin{equation}{section}` para corrigir o contador do template;
  e o título "Sparse wavelet estimation of additive functional
  coefficients", com título corrente de 40 caracteres.
- **O `chicago.bst` do template não obedece à própria revista:** abrevia em
  "et al." a partir de três autores, e a §4 das instruções manda listar os
  três. Contornado com `\citet*` onde deu; decidir antes da submissão se
  ajusta o `.bst`.
- **Cinco referências que faltam** e que por isso ficaram sem citação
  (`glmnet` e `WaveBased` aparecem em `\texttt{}` sem `\cite`): Tibshirani
  (1996), Friedman, Hastie & Tibshirani (2010), a citação do R, o
  `sparsegl` e o Johnstone de modelos de sequência.

**E2.3, sintonia** (`wafc/R/tune.R`, `tests/test-tune.R`, mais os dois
drivers salvos aqui em `wafc/scripts/`). `cv.wafc()` sobre `(J, λ)` com
dobras fixas, BIC, EBIC, a regra da teoria e `wafc_tune()` como entrada
única. A comparação rodou nos três cenários, `n ∈ {250, 1000}`, 20 réplicas,
contra o **oráculo da grade** (o par `(J, λ)` que minimiza o erro fora da
amostra, que nenhuma regra enxerga).

- **A validação cruzada em `lambda.min` é a melhor regra de predição**, com
  custo de `1.00` a `1.04` sobre o oráculo, e acerta o `J` do oráculo em
  20/20 réplicas com `n = 1000`. O preço é estrutura: liga os três blocos
  nulos sempre, e no cenário nulo deixa 1 ou 2 coeficientes de lixo, onde
  `lambda.1se`, BIC, EBIC e a regra da teoria zeram tudo.
- **A regra da teoria não é utilizável como regra prática, e o culpado é o
  `λ`, não o `J`:** `λ_n` do Corolário 2 é de **7 a 13 vezes** o `λ` que
  minimiza o erro realizado; no cenário não homogêneo, `J_n` sozinho custa
  2% e 18% enquanto `λ_n` sozinho custa 30% e 54%. É a distinção usual entre
  otimizar a cota e otimizar o erro, e é dela que a seção de computação
  vive.
- **A direção do erro de `J_n` depende de `s'`**, e isso corrige a leitura
  que E1.6 sugeriu: lá, com `s' = 1/4`, a regra ficava dois níveis *acima*;
  aqui, com `s' = 3/2`, fica dois níveis *abaixo*. A frase para o artigo não
  é "a regra fica acima" e sim **"a regra é muito sensível a `s'`, que não
  se conhece"**. A frase de §4.3 do `ms_1.tex` que diz "pode ficar acima"
  tem de virar número ou sair.
- **O EBIC penaliza a resolução, não só o modelo** (`2γ log C(d,k)` cresce
  com `d = pq(2^J − 1)`), e por isso erra o `J` sistematicamente quando `n`
  cresce: fica em `J = 3` em 17/20 e 18/20 réplicas onde o oráculo está em
  4 e 5. Se E2.4 quiser critério de informação, o BIC.
- **`σ` estimado piora tudo mais um passo:** estimado em `J_n`, `σ̂` dá 1.65
  a 1.82 vezes o `σ` verdadeiro, porque ali o viés do sieve está no resíduo.
- **Achado que pede conserto em `design.R`** (nenhuma das duas tarefas podia
  tocar): `wafc_design()` não constrói desenho em `J = 1` com os padrões,
  porque `eps = 1.9^{−J} = 0.526` cai fora do `[0, 0.5)` que
  `wafc_rescale()` exige. A grade de `cv.wafc()` começa em `J = 2` por isso.

### 2026-09-19: L4 fechada (referências de software)

`referencias-verificadas.bib` vai de 57 a **63 entradas**; conferido aqui:
63 `@`-entradas, 63 `\bibitem` e **zero erro e zero aviso** de BibTeX, tanto
com `plain` quanto com o `chicago.bst` do template da revista. As cinco
referências que faltavam existem, e o manuscrito já pode citar o lasso, o R,
o `glmnet`, o `sparsegl`, o `WaveBased` e o Johnstone.

- **Correção herdada do WALL:** a entrada `johnstone2019gaussian` de lá
  atribui o livro à Cambridge University Press, e isso **não se confirma**:
  não há registro no Crossref e a busca em livros da CUP devolveu zero em
  2026-09-19. A página do autor e a folha de rosto dizem "Book Draft,
  version of September 16, 2019", sem editora. Aqui ficou `unpublished`. Se
  o WALL for submetido citando como livro da CUP, a referência sai errada lá
  também.
- **Duas armadilhas de composição**, anotadas no próprio `.bib`: um arroba
  dentro de comentário quebra o BibTeX, e o `chicago.bst` avisa "empty
  organization" num `manual` sem organização, razão de o `WaveBased` ter
  entrado como `misc` com `howpublished`.
- O `glmnet` e o `sparsegl` não têm páginas no Crossref (o JSS não as
  registra); vieram do `CITATION` dos pacotes instalados. O DOI do R não é
  de versão, então o ano é o da versão instalada e a versão vai na `note`.

### 2026-09-19: E1.3b e E2.1b fechadas, e a margem mostra o seu preço

Dois chats de tarefa. Conferido aqui: a bateria de `wafc/` dá **420 passam,
0 falham**; `02-aproximacao-besov.tex` compila em 11 páginas sem referência
indefinida.

**E1.3b, o lema de extensão.** Hipótese 5 (margem fixa e Besov **no
suporte**), Lema 10 em quatro itens e três observações. O que D23 propunha
está provado: com margem, o viés volta a `O(2^{−2Js'})` e o termo `2^{−J}`
da periodização some. Mas a emenda cobra dois preços que não estavam na
conta:

1. **O preço da margem é exato:** `C(eps) ≍ eps^{−(s−1/π)}`, com cota
   superior provada e **cota inferior para qualquer extensão admissível**
   (mergulho de Besov em Hölder). Não é artefato da construção. Medido: a
   norma da extensão de `u − 1/2` vai de `0.92` a `46.7` quando `eps` cai
   de `0.4` a `0.05`, com expoente convergindo para `s − 1/π`.
2. **A margem não pode encolher com `J`.** Com `eps ∝ 2^{−J}` o ganho é
   **zero** (medido: queda de `0.52` a `0.58` bit por nível, contra `0.50`
   sem margem nenhuma); com o `1.9^{−J}` herdado do `wall()` o ganho é
   parcial (`0.96` bit); só margem **fixa** compra a taxa cheia (`0.05` e
   `0.10` dão de 4 a 13 bits por nível até o piso numérico).

**A tensão que isso abre, e que é o achado mais importante da rodada.** A
margem que compra a taxa é a mesma que degenera a Gram restrita:
`λ_min(G_eps) = 0` assim que `V_J` contém função suportada na margem, isto
é `2^J ≥ (L−1)/(2 eps)`. Medido: `λ_min(G_eps)` estabiliza em `0.213` com
`eps = 2^{−(J+1)}` e em `0.0013` com `1.9^{−J}`, e desaba a `10^{−12}` em
`J = 6` com `eps = 0.05`. Não é contradição — as direções que somem são
wavelets invisíveis no suporte, cujas colunas de `Z` se anulam e que não
afetam `f` ali —, mas **o que E1.5 consome é uma constante de
compatibilidade sobre as direções estimáveis, e ela não existe hoje**.
Enquanto não existir, a taxa de E1.6 sob a Hipótese 5 é **conjectura, não
teorema**, e está escrita assim no arquivo. O manuscrito não é afetado: ele
enuncia a teoria na base periodizada com `eps = 0` (D23 já mandava a margem
para a computação).

**E2.1b, o conserto do `eps`.** O padrão deixou de depender de `J`:
`wafc_eps()` passou a receber `q` e `boundary`, e devolve `0` no caso
`"interval"` e a constante fixa **`0.05`, declarada provisória**, no
periódico. `J = 1` deixou de ser inalcançável e virou candidato que a
validação cruzada descarta com número (`mse` `2.446` contra `0.669` em
`J = 3`). Colunas identicamente nulas aparecem, como previsto, a partir de
`J = 8` com `eps = 0.05` (12, 48, 136 colunas em `J = 8, 9, 10`), e nenhuma
das duas penalidades quebra.

**Três linhas de `wafc/R/tune.R` ficaram factualmente falsas** e nenhuma
das duas tarefas podia tocá-las: o Roxygen de `cv.wafc(J = )` e o
comentário de `wafc_J_grid()` justificam a grade começar em `J = 2` por uma
razão que morreu, e a guarda da regra da teoria **para em erro** quando
`J_n = 1`, o que já não é necessário. É conserto de uma linha cada, para
E2.4 ou para uma correção curta. O `inventario-codigo.md:45` também
documenta o padrão antigo.

### 2026-09-19: E1.8 fechada (a rota do intervalo está escrita)

Conferido aqui: `07-rota-intervalo.tex` compila em 11 páginas sem
indefinida, e `check/07-rota-intervalo.R` imprime `OK` (2 min 11 s). O
arquivo **não vai ao manuscrito** (D26); existe para estar pronto.

- **A rota do intervalo fecha, e o preço é todo de amostra finita.** As
  quatro peças transferem: aproximação com a mesma prova (some o operador
  de extensão, o corte, a constante `C(eps)` e o termo `2^{−J}`);
  identificabilidade, que já estava feita na §5 de E1.2; desenho, verbatim,
  com `λ_min(Σ) ≥ κ_1 c_U` sem `λ_min(G_eps)`; oráculo e taxas, que
  consomem só isso. O que é novo e foi provado: **Lema 11**, a cota pontual
  `Σψ² ≤ C_ψ 2^J` na CDV, e **Proposição 5**, o efeito dos coeficientes de
  escala não penalizados, que troca `σ²p/n` por `σ²p_0/n` com
  `p_0 = p{1 + q(2^{j_0} − 1)}`, constante em `n`, de modo que as taxas não
  mudam.
- **O que ela cobra, medido:** `C_ψ^int ≈ 7.55` contra `1.30` na
  periodizada (fator 5.6, que a condição empírica de E1.4(iv) paga em `n`);
  `p_0 = 93` contra `3` com `p = 3`, `q = 2`, `j_0 = 4`, o que dá
  `0.372 σ²` em `n = 250`; e `j_0 ≥ 4`, que tira `J ∈ {1,2,3,4}` da grade —
  **exatamente a faixa que E2.3 achou ótima** nos `n` do piloto.
- **O cruzamento da margem agora é teorema** (Proposição 6): as duas
  exigências medem a mesma largura `(L−1)2^{−J}` e as faixas são disjuntas,
  separadas por **um nível exato**. Manter a Gram viva com `J_n → ∞` obriga
  `eps_n ≍ 2^{−J_n}`, onde o expoente é `1/2`, o mesmo de margem nenhuma. O
  Corolário 7 diz qual enunciado é verdadeiro: amostra finita, `eps` fixo,
  `J` na faixa admissível.
- **Restrição prática para E2.4:** com `eps = 0.05` a faixa admissível
  medida é `J ≤ 5`, e E2.3 achou o ótimo em `J = 3` a `5`. O padrão
  provisório está dentro, mas por pouco: `eps = 0.10` limitaria a `J ≤ 4`.
  É razão contra aumentar a margem.
- **A análise de `s'` (Lema 12 e a tabela das seis componentes)** confirma
  D27 e explica: das seis, só a cúbica paga periodização, e paga na
  **derivada** (`g'(0) = 0.4` contra `g'(1) = 0.6`); as três de Donoho e
  Johnstone são contínuas na emenda. Medido, `smooth` lê `3/2` no regime do
  manuscrito e `4` com margem ou na CDV; `inhomogeneous` lê `1/2` nos três.
  **Na CDV a cúbica é reproduzida exatamente** (grau `3 < N = 4`), com
  coeficientes de wavelet em `1e-10`.
- **Lição de ferramenta:** a base CDV do `wbasis()` perde precisão com o
  nível — o piso de avaliação vai de `7e-08` em `j = 4` a `3e-04` em
  `j = 12`, contra `1e-16` na periodizada. Medir decaimento na CDV acima de
  `j ≈ 7` lê ruído, não coeficiente.

**Registro de um erro meu no versionamento.** O commit `c911cf9`, que era
de plano, foi feito com `git add -A` em vez de caminhos explícitos, contra a
§4 do `instrucoes.md`. Ele levou junto os entregáveis de E1.8, ainda não
integrados, e os de E2.4, **em curso e sem handoff**. Não houve reescrita de
histórico; o commit fica, e esta é a integração que devia tê-lo precedido.

### 2026-09-20: E1.7a, E1.7c e E1.4c fechadas; a seleção de estrutura tem teorema

**E1.7c, a seleção por limiarização (saída (c) de D28), fecha e é barata.**
`06-selecao-limiar.tex`, 7 páginas, com **Lema 13** (determinístico, três
linhas: o sanduíche `{‖g‖ > t + D} ⊆ Ŝ(t) ⊆ {‖g‖ > t − D}`) e **Corolário
8**, que instala nele a taxa do Corolário 4. Nenhuma condição de desenho
nova, nenhuma irrepresentabilidade, e vale para o **LASSO puro** (D3).

- **O achado que não estava previsto: triagem e recuperação se separam.**
  Eliminar os blocos nulos **não custa hipótese nenhuma** (parte (ii) do
  Corolário 8); toda a hipótese de separação é gasta em não perder bloco
  ativo pequeno. É isso que explica o cenário não homogêneo acertar 1.00
  mesmo com a janela do Lema 13 vazia: o viés do sieve encolhe todas as
  normas na mesma direção e não cria falso positivo.
- **A comparação com a variante em grupos mudou de sinal.** Em `n = 1000`,
  `J = 4`, mesmas dobras: o LASSO limiarizado acerta a estrutura em **10 de
  10** réplicas nos dois cenários, e o sparse group LASSO lido em
  `lambda.min` acerta em **0 de 10**, com predição igual (`RMSE` 0.770
  contra 0.776 no suave). A pergunta de E2.5 deixou de ser "grupos ou nada".
- **O platô de limiares é largo:** `P(Ŝ = S) = 1` para `t ∈ [0.062, 0.590]`
  em `n = 2000`, e a regra grosseira `t = 0.15 max‖ĝ‖` acerta 1.00 e 0.90 nos
  dois cenários em `n = 1000`. Ponto de partida para E2.4/E2.5.
- **`wafc_blocks()$norm` já é a estatística do limiar**, exata em `eps = 0`
  (erro relativo `1.4e-09`); com `rescale = TRUE` entra um fator
  `sqrt(scale_m)`, comum aos blocos que dividem `U_m` mas diferente entre
  moduladoras. É calibração para E2.5, não defeito.
- **Quarta aparição do mesmo regime pré-assintótico:** a regra teórica
  `(J_n, λ_n)` zera quase tudo e o melhor limiar vira o piso da grade; em
  `n = 250` no cenário não homogêneo ela acerta **0**.

**E1.7a, a sondagem da saída (a): veredito "vale a pena, com escopo
reduzido, e não antes de E2.5".** A álgebra fechou melhor do que a
conjectura previa: sob `E(XX' | U)` constante e **(BD)**, independência
entre moduladoras, a condição de irrepresentabilidade em grupos de Bach
(2008), já com as colunas não penalizadas perfiladas, vira
`IRR = sqrt(a'V^{(m)}a)` com `a = Ω_{ℓ_0,S}(Ω_{S,S})^{-1}` — **a base, o
nível `J`, `N_J` e as densidades marginais cancelam**.

- Conferido: a redução bate com a `Σ` cheia a `1.0e-17`, e `IRR` é
  **constante em `J`** (`0.846` de `J = 3` a `7`). **A dificuldade de seleção
  no desenho de produtos não vem da alta dimensão em wavelets**; vem de `Ω` e
  do ângulo entre as componentes que dividem uma moduladora. É resposta
  direta à pergunta de referee sobre seleção, e custa meia página.
- **Minha conjectura estava meio certa:** não se reduz a `Ω` sozinha, sobra
  `V^{(m)}`, a matriz de correlações `L_2` das componentes. Mas
  `sqrt(a'Va) ≤ ‖a‖_1`, e `‖a‖_1 < 1` é a condição de Zhao & Yu (2006): a
  conjectura vale como **suficiente**, com folga de até `sqrt(|S_m|)`.
- **O obstáculo real é D13**, `X` dependente de `U`: o erro da redução cresce
  com a dependência (`−0.035`, `−0.073`, `−0.132`). Densidade marginal não
  uniforme não faz nada, porque o complemento de Schur absorve.
- **Margem populacional abaixo de `0.2` não se traduz em seleção nos `n` de
  E2:** com `η = 0.21`, uma réplica em 50 já viola em `n = 250`.
- **Rota alternativa descoberta ao conferir a referência:** Wei & Huang
  (2010) **não usam irrepresentabilidade**; provam seleção para o LASSO
  adaptativo em grupos a partir de condição de Riesz esparsa, que é
  essencialmente o que a Proposição 3 de E1.4 já entrega. Se (a) for
  reaberta depois de E2.5, é por aí, não por Bach.
- **Aviso para E2.4 e E4:** no cenário `smooth` de `dgp.R`,
  `⟨sine, cubic⟩ = 0.969` — as duas componentes de `β_1` são quase
  proporcionais, o que é o pior caso para qualquer medida de acerto de
  estrutura. Não afeta predição.

**E1.4c** traduziu o `03-desenho-produtos.tex` (D29), recompilado em 6
páginas, com `babel` acrescentado e uma citação corrigida (a frase atribuída
ao `proposta-metodo.md` era "no concurvity", e o documento diz "não
colinearidade entre moduladoras").

### 2026-09-20: E2.4 fechada, e o piloto achou um erro de grade

Conferido aqui: **544 testes passam** (eram 420), `01-smoke.R` imprime `OK`,
e uma passagem curta de `04-pilot.R` reproduz o padrão da tabela. O piloto
completo levou ~3 h em 12 núcleos e não pôde ser repetido aqui.

**O achado que decide E2.5: a grade de `J` herdada do `cv.wall()` trunca o
sieve no ótimo.** O oráculo da grade escolhe o `J` do **topo** da grade em
100% das réplicas de toda célula com componente. Alargando de
`⌈log_2 n/2⌉` para `2:8`, no cenário não homogêneo em `n = 1000`, o erro de
predição cai **17,6%** e o ISE **32%**, em **50 de 50 réplicas**; no cenário
suave não muda nada, nem o `J`, nem o `λ`. A razão está em E1.3: `bumps` só
resolve de `J = 9`. O preço é tempo (`5,3×` a `5,6×`).

**Isso inverte o veredito contra os concorrentes.** Com a grade curta, o
WAFC **não ganha** o cenário (b): fica atrás do block LASSO de Klopp &
Pensky (0.96 a 0.98 em toda a faixa) e o **VCBART ganha dele em `n = 1000`
em 50 de 50 réplicas** (0.899). Com a grade profunda, o WAFC passa a
**0.839** contra `0.910` do VCBART, `0.997` do `klopp` e `1.189` do `gam`,
isto é, 7,8% melhor que o melhor concorrente. **E2.5 tem de repetir a parte
`competitors` com a grade larga antes de decidir**: o `jgrid` mediu o WAFC
contra si mesmo.

**No cenário (a) o `gam` ganha do WAFC por 2,2 a 2,6 em ISE** (1,41 a 1,54
em `rmse_f`), acima do fator 1,5 que E2.5 propõe, e a grade profunda não
muda isso. Mas o cenário suave medido é o **antigo**, de curvatura uniforme,
que é o melhor caso do spline penalizado; a componente `C^∞` de **curvatura
desigual** que D30 pede não existe em `dgp.R`, e é ela que testaria a
afirmação de D30. Ordem proposta e que eu endosso: acrescentar a componente,
repetir o cenário (a) com ela e com a grade larga, **e só então** fixar o
fator do critério de saída.

Outros números da rodada:

- **O block LASSO de K&P, no mesmo desenho, prediz melhor que o LASSO puro
  no cenário não homogêneo** (0.96 a 0.98) e pior no nulo (1.30 a 1.50),
  onde penalizar os `c_ℓ` os encolhe. É a resposta numérica à pergunta que
  D18 convida, e o artigo tem de escrevê-la.
- **O QUT domina a regra da teoria uniformemente** e dispensa `σ` e `s'`;
  nenhuma das duas serve como regra de predição (custo 1,3 a 5,0). Mas o
  QUT é a **única regra que nunca liga um bloco falso**, em nenhuma célula:
  é regra de estrutura, não de predição.
- **D20 se sustenta:** `cv.min` custa 1.00 a 1.03 fora do cenário nulo e é a
  melhor regra de predição em todas as células. E **nenhuma variante
  seleciona estrutura com `cv.min`** — o que confirma E2.3 e é exatamente o
  que a limiarização de E1.7c resolve.
- **O padrão `wafc_eps_periodic = 0.05` é dominado:** custa 12% a 15% no
  cenário suave (perdendo em 95% a 100% das réplicas) e não ganha nada no
  não homogêneo (0.992 e 0.994). `eps = 0` ganha de todas as margens no
  suave, perde no máximo 2,4% no não homogêneo, e é o único valor em que
  `λ_min(G_eps) = 1`, isto é, em que a condição de E1.4 vale com a constante
  que a teoria usa.
- **`J = 1` não entra na grade:** escolhido em 0 de 50 réplicas em toda
  célula com componente.
- **D31 não muda conclusão nenhuma e vale 27% do tempo:** mesmo `(J, λ)`,
  `mse` igual na oitava casa, `1.11e-05` de diferença máxima no desenho.
- **Orçamento medido:** a célula `mixed` custa ~23 horas-núcleo para 50
  réplicas contra 3 das outras três juntas, e rodou com **15 réplicas**; as
  demais com as 50 do plano.

### 2026-09-20: E6.1a fechada, com veredito negativo e uma lição que atinge E2.4

**Nenhuma das três candidatas sustenta o argumento do artigo.** O detalhe
está em [`aplicacao-candidatas.md`](aplicacao-candidatas.md); o essencial:

| base | `n` | WAFC | `gam` `k = 10` | `gam` de **dimensão casada** |
|---|---|---|---|---|
| bike | 17 379 | 0.6219 | 0.6480 | **0.6209** |
| beijing | 34 287 | 0.9010 | 0.9009 | **0.9007** |
| housing | 20 640 | 0.2736 | 0.3151 | **0.2713** |

**O ganho aparente era de dimensão, não de base.** Com `k = 10` o WAFC
parecia ganhar de 6,6% a 15,2%; com `k = 2^J` ele empata, e o spline fica
marginalmente à frente nas três. A maior diferença a favor do WAFC que
sobreviveu a um desenho honesto foi **0,4%, numa partição só**.

A estrutura diz o mesmo por outro caminho: a componente do WAFC é **duas
vezes mais rugosa** que a do spline de mesma dimensão, com o **mesmo** erro
de predição, e em bike e housing a energia por nível **não decai**
(`ŝ'` mediano `−0.12` e `0.06`). Rugosidade que não paga em predição é ruído
ajustado, não estrutura.

**A lição que atinge E2.4 e E4, e que é maior que a tarefa.** A comparação
com o spline tem de ser em **dimensão casada**, senão mede a coisa errada.
E `wafc_fit_gam()` tem **`k = 10` como padrão**, que é o que o piloto usou:
um spline de dimensão 10 contra um sieve com 31 funções por bloco em
`J = 5`, ou 127 em `J = 7`. Portanto **a vitória do WAFC sobre o `gam` no
cenário não homogêneo (1.006 a 1.178) é suspeita**, e a derrota no cenário
suave (0.65 a 0.71) é, se algo, conservadora. Isso entra em E2.5 antes de
qualquer veredito: é a segunda coisa, junto com a grade de `J`, que precisa
ser remedida.

- **`mgcv::bam` é duas ordens de grandeza mais rápido que `gam` nestes `n`**
  (6,1 s contra 1453 s para oito suavizadores de `k = 23`). Um estudo que
  sintonize o spline com `gam` vai, por custo, dar a ele base pequena demais
  e medir a diferença errada. É acionável em E2.4b e em E4.2.
- **Partição aleatória em série dependente inverte o veredito, e isso é
  condição de validade, não refinamento.** Em bike, sob partição aleatória a
  CV pedia `J = 8` e o WAFC parecia ganhar 0,8%; retendo **semanas
  inteiras**, a curva ganha mínimo interior em `J = 4` e a diferença cai a
  0,4%. O ganho vinha de o componente fino no índice do dia carregar o nível
  do dia do treino para o teste. Vale para E6 e para qualquer tabela de E5b
  com dado real estruturado no tempo ou no espaço.
- **Moduladora discreta com poucos valores é incompatível com o argumento:**
  com `hr` em 24 valores, o bloco tem posto no máximo 23 e não há onde
  exibir regularidade heterogênea, por maior que seja `J`.
- **O que a próxima busca tem de exigir antes de baixar dado:** razão
  substantiva, documentada na literatura da área, para o efeito de `X_ℓ`
  **saltar** num valor conhecido de `U_m`. As famílias que têm isso: limiar
  administrativo (faixa de imposto, nota de corte, elegibilidade), limiar
  regulatório em dose-resposta, e série com quebra datada.
- **Dois diagnósticos que E4 vai querer:** a energia por nível com reajuste
  sobre o suporte selecionado, e o índice de localização. Eles separam
  "ajustou estrutura" de "ajustou ruído" **sem precisar da verdade**, que é o
  que uma aplicação real nunca tem.

### 2026-09-21: medição no chat principal, dois defeitos de desempenho

Antes de decidir P1, fui ver por que o `bsgl` era 8 vezes mais rápido que o
WAFC na tabela do piloto. **Não é o método, é o nosso ajuste.**

- **A tolerância padrão de `wafc()` é `thresh = 1e-10`**, três ordens abaixo
  do padrão do `glmnet`. Medido em `n = 1000`, `p = 3`, `q = 2`, mesmas
  dobras: a validação cruzada do candidato `J = 5` custa **5,95 s** com
  `1e-10` e **0,21 s** com `1e-7` (**28×**), e a sintonia **não muda**:
  `lambda.min` idêntico, `cvm` mínimo diferente na quinta casa. O `glmnet`
  sozinho, no mesmo desenho, leva `0,01 s`; o resto é a tolerância.
- **Em candidatos de tamanho comparável o WAFC é mais rápido que o `bsgl`:**
  em `J = 4`, que tem as mesmas 93 colunas do maior candidato dele, o ciclo
  completo custa `0,15 s` contra `0,77 s` do `bsgl` inteiro.
- **O `...` de `cv.wafc()` é roteado para `wafc_design()`**, então
  `cv.wafc(..., thresh = 1e-7)` para com "argumento não utilizado": não há
  como afrouxar a tolerância pela interface, o que provavelmente é a razão
  de o defeito ter passado despercebido.

**Consequência para P1:** a objeção de custo à grade profunda encolhe muito.
O piloto estimou 5,3× a 5,6× para alargar a grade; com o ajuste 28× mais
barato no candidato dominante, a grade larga cabe, e a comparação com o
`gam` em dimensão casada (E6.1a) também.

### 2026-09-21: E1.10 fechada (a terminologia)

"Sieve" saiu das seis derivações e do comentário do `macros.tex`, com os seis
`.pdf` recompilados; conferido aqui, os seis compilam sem referência
indefinida e nenhum mudou de tamanho (11, 6, 10, 10, 7 e 11 páginas). O diff
é inteiramente de palavras: nenhum enunciado, hipótese, constante, número ou
`\label` mudou, e nenhum `\ref` foi tocado.

O dicionário, que o manuscrito reaproveita em `k = 2`: o sieve → **o espaço
de aproximação**; sieve linear → **aproximação linear** (para fazer par com
"aproximação não linear", que já era o termo do texto); viés do sieve →
**viés de aproximação**; regime de sieve → **regime de resolução**. Ficou
**uma** menção, entre parênteses, no ponto em que `W_J` é definido.

### 2026-09-21: E2.4b fechada, e a minha leitura da tolerância estava errada

Conferido aqui: **642 testes passam** (eram 544).

**A correção que mais importa, e é sobre o que eu afirmei ontem.** Eu medi
28× entre `thresh = 1e-10` e `1e-7` e concluí que o padrão apertado era
desperdício. A varredura de E2.4b, em 12 células, mostra que **não é**: com
`1e-7` o `wafc_kkt()` **rejeita pontos do caminho em todas as 12 células**
(de 19 a 73 pontos, de caminhos de 58 a 100) e com `1e-8` em dez delas; com
`1e-9` rejeita em uma e com `1e-10` em nenhuma. E a minha frase "a sintonia
não muda" vale **a partir de `1e-8`**, não de `1e-7`: em `smooth`,
`n = 250`, `J = 5`, o `1e-7` escolhe `λ = 0.0914` contra `0.1003` das
demais. O ganho de tempo está quase todo entre `1e-9` e `1e-7`, isto é,
**os 8× restantes se compram entregando a verificação**, não cortando
folga. O padrão novo é `1e-9`, o mais frouxo que mantém o KKT limpo, e o
conserto do `...` é o que deixa um estudo escolher `1e-7` **de propósito**.

**O `...` de `cv.wafc()` está consertado**, com separação por nome e recusa
de nome que não pertença a nenhum destino, de modo que erro de digitação
vira erro em vez de padrão silencioso.

**A tabela de tempo, com a separação que faltava** (não homogêneo, `q = 2`,
medianas; `secs.tune` é a busca, `secs.fit` o ajuste escolhido):

| `n` | `wafc.lasso` | `bsgl` | `gam k=10` | `gam` casado (`bam`) | `klopp` | `vcbart` |
|---|---|---|---|---|---|---|
| 250 | **0,375** (0,371 + 0,005) | 0,424 | 0,855 | 0,094 | 0,901 | 2,482 |
| 1000 | **3,105** (3,096 + 0,010) | 0,694 | 1,955 | 0,589 | 5,641 | 6,266 |

O WAFC deixou de ser o caro da tabela: em `n = 250` ele é mais rápido que
todos menos o `gam` por `bam`, e em `n = 1000` perde só para `bsgl` e para o
`bam`. E a coluna `split` registra o que a tabela do piloto escondia: o
tempo do WAFC **expõe** a busca de sintonia, o dos outros a esconde dentro
de uma chamada.

**O cenário `uneven` de D30 nasceu**, com `gaussians` (razão de escalas
5,6), `chirp` (frequência de 1 a 8 ciclos) e `cosine`, todas `C^∞` por
construção. Medido na maquinaria de E1.8: queda por nível de **3,91** e
**3,63**, contra 3,99 do `sine`, isto é, as duas leem o **teto do filtro**
(`N = 4`), que é o que autoriza `s' = 4`. A janela `C^∞` do Lema 10 não é
enfeite: sem ela a `gaussians` tinha `g(0) = 0.044` contra `g(1) ≈ 3e-11`, e
o salto na emenda derrubava a queda para **0,76**, assinatura de
descontinuidade. E o cenário novo **não é adversário gratuito**:
`⟨gaussians, chirp⟩ = −0,126`, contra `⟨sine, cubic⟩ = 0,969` do `smooth`.

- **`s'` virou atributo** (`wafc_sprime(scenario, regime)`), com o **regime**
  anotado: `smooth` lê `3/2` em `periodic` e `4` com margem ou no intervalo;
  `inhomogeneous` lê `1/2` nos três; `uneven` lê `4` nos três.
- **Quinta aparição do regime pré-assintótico:** no `uneven` a regra da
  teoria **colapsa em `J = 1` e zera tudo**, porque `s' = 4` faz
  `2^J ≥ (n/log n)^{1/9}` pedir `2^J ≥ 1.6`. A coluna `theory` desse cenário
  vai ser degenerada em E2.5, e isso tem de ser dito na tabela.
- **P1 e P2 não foram ratificadas** (o autor respondeu no chat da tarefa):
  a grade continua em `⌈log_2 n/2⌉` e a margem em `0.05`. **E2.5 herda as
  duas**, e a evidência de E2.4 a favor de P1 continua na mesa.
- O `gam` ganhou `k` por moduladora e motor `bam`, e `wafc_k_matched()`.

**Sobre o alarme do `vcbart`: o handoff suspeitou que o piloto tivesse
perdido a coluna, e eu fui conferir.** Ela não se perdeu. O `vcbart` **está**
na tabela de razões da análise de E2.4, com `rr` de `0.999`, `0.992` e
`0.899` no não homogêneo e `0.964`, `0.972`, `0.889` na `mixed` — os números
que o `ESTADO.md` registra. O que aconteceu é outra coisa, e em duas
camadas: ele some da tabela de **medianas absolutas** porque não decompõe em
componentes, então tem `ise = NA` e o `aggregate()` descarta a linha; e a
chamada **hoje** falha, porque o ajuste a D31 passou a mandar
`wavelet.table` a todo concorrente e o `VCBART` não conhece o argumento —
ajuste que, como o próprio handoff de E2.4 registra, **veio depois** de os
números serem produzidos. Ou seja: **os números valem; o script, como está,
perderia a coluna em silêncio na próxima execução.** É conserto obrigatório
antes de E2.5 rodar qualquer coisa.

**Lição de catálogo, segunda vez:** a coluna de arquivos permitidos de E2.4b
não cobria `wafc/R/fit.R` nem o `06-timing.R`, que o próprio texto da tarefa
mandava tocar. Antes disso, E2.4 pôs o QUT em `competitors.R` pela mesma
razão. **Quando eu catalogar, o arquivo vai junto do item.**

### 2026-09-21: a tabela de tempo, medida no chat principal

`06-timing.R`, 5 réplicas, cenário não homogêneo, `p = 3`, `q = 2`, nesta
máquina. Medianas em segundos; `busca` é a sintonia e `ajuste` é o ajuste
que ela seleciona.

| método | 250 | 500 | 1000 |
|---|---|---|---|
| `gam` casado, por `bam` | **0,100** | **0,313** | **0,661** |
| `bsgl` | 0,42 | 0,453 | 0,670 |
| `gam k = 10` | 0,86 | 1,051 | 2,019 |
| **`wafc.lasso`** | **0,361** | **6,856** | **3,204** |
| `klopp` | 0,90 | 6,209 | 5,917 |
| `vcbart` | 2,48 | 3,831 | 6,482 |
| `aspline` | 3,70 | 4,801 | 6,574 |
| `gam` casado, por `gam` | 1,49 | 10,013 | **31,393** |
| `wafc.sglasso` | 6,94 | 39,342 | 37,507 |

- **O ajuste do WAFC custa de 5 a 57 milissegundos**; todo o resto é a busca
  sobre `(J, λ)`, e é a única coluna da tabela que **mostra** a busca: o
  `gam` escolhe suavização por REML dentro da chamada e o `bsgl` faz
  validação cruzada dentro do `cv.grpreg`.
- **`n = 500` é mais caro que `n = 1000` no WAFC** (6,86 contra 3,20), e o
  mesmo aparece no `klopp` (6,2 contra 5,9) e no `sglasso` (39,3 contra
  37,5), que rodam no mesmo desenho: em `n = 500` a grade chega a `J = 5`
  com 31 colunas por bloco e poucas linhas, que é o pior ponto da razão
  entre dimensão e amostra. **É do desenho, não do laço de validação
  cruzada.**
- **A comparação em dimensão casada só é viável pelo `bam`:** `gam()` custa
  **31,4 s** em `n = 1000`, quase cinquenta vezes o `bam` e dez vezes o
  WAFC. É a razão pela qual o piloto rodou em `k = 10`, e é o obstáculo que
  E2.4b removeu.
- **O `sglasso` custa 19× o LASSO em `n = 250` e 12× em `n = 1000`.** Somado
  ao que E1.7c mediu (o LASSO limiarizado acerta a estrutura em 10 de 10
  réplicas onde ele acerta 0 de 10), a variante em grupos chega a E2.5
  perdendo nos dois eixos.

### 2026-09-28: revisão do código de `wafc/` (defeitos, memória, velocidade)

Leitura dos sete arquivos de `wafc/R/` e dos seis scripts, com medição de
tempo e memória no chat principal; aplicado o grupo que **não muda número
nenhum** (nem `J`, nem `λ.min`, nem `λ.1se`). Conferido aqui: **672 testes
passam** (eram 642; 30 novos), e o `04-pilot.R` antigo e o novo, rodados
em 2 réplicas de `smooth` com `n = 250` nas partes `competitors` e
`lambda`, dão **diferença máxima 0** nas 15 colunas numéricas das 30 linhas
comuns.

- **Defeitos corrigidos.** (1) `wafc_competitor()` mandava `wavelet.table`
  ao `VCBART_ind`, que parava com "argumento não utilizado" (pergunta 10):
  agora argumento de `wafc_design()` só chega a quem constrói desenho
  (`klopp`, `oracle`), e o `vcbart` voltou à tabela. (2) O `04-pilot.R`
  tinha um vetor `wafc_sprime` que escondia a função e não tinha `uneven`:
  a célula nova sumiria da parte `lambda`; agora lê `wafc_sprime()` do
  `dgp.R`. (3) Falha de método virava ausência de linha (`try(silent)` +
  `next`), que é como o `vcbart` sumiu sem aviso: agora vira linha sem
  números com a mensagem na coluna `error`, e o `sweep_part` conta as
  falhas na saída. (4) Sem covariável constante, `bsgl` e `klopp` somavam o
  intercepto do `grpreg` a `β_1`, multiplicando-o por `X_1` (RMSE de 7,1
  com nível 5): agora ele fica em `intercept`, fora de `β`. Não afeta o
  piloto, que sempre tem `X_1 ≡ 1`. (5) `wafc_lambda_max()` punha coluna de
  intercepto mesmo quando o ajuste não tem: o topo do caminho saía 0,760 em
  vez de 1,029; agora bate com o do `glmnet` nos dois casos.
- **Memória.** Chamados por `do.call`, `glmnet` e `sparsegl` guardavam no
  `call` o desenho inteiro, e `wafc_design()`, `wafc()` e `cv.wafc()`
  guardavam os dados e o corpo da própria função: era **49% de cada cache**
  do `05-sondagem` (97 a 115 MB), e `print(fit$fit)` imprimia 534 mil
  caracteres. Agora o `call` é compacto. E o `cv.wafc()` guardava o ajuste
  de todo candidato de `J` até o fim; agora só o melhor. Medido: ajuste em
  `J = 6`, `n = 1000`, de **7,35 para 4,35 MB** salvo; `cv.wafc` com grade
  `2:8`, de **11,4 para 3,6 MB**.
- **Velocidade.** A busca de nós do spline adaptativo pontua todos os
  candidatos com uma projeção (acrescentar o nó `κ` acrescenta
  `(v − κ)³₊` ao espaço), em vez de um ajuste por candidato: **mesmos nós
  em 18 de 18 casos**, valores ajustados idênticos, **7 a 11× mais
  rápida** (no piloto, 6,95 s para 0,61 s). A parte `lambda` do piloto
  rodava o mesmo `cv.wafc()` três vezes (`cv.min`, `cv.1se`, `qut`); agora
  uma, com a mesma coluna `time` (cada regra continua carregando o tempo
  da validação cruzada que usa).
- **Lição de ferramenta:** um limiar relativo de `1e-10` para descartar
  candidato degenerado na busca de nós zerava candidatos legítimos (um nó
  entre dois vizinhos acrescenta pouco, mas acrescenta) e mudou os nós em
  5 de 18 casos; `1e-24` separa isso do arredondamento, que é da ordem de
  `1e-32`.

**O que ficou para decidir, com número** (pergunta 29): 97% do tempo do
WAFC está dentro do Fortran do `glmnet`, e a cauda do caminho de `λ`
abaixo de `imin + 20` custa **90% a 98%** das dobras em `J = 5` e `6`
(`J = 6`, `n = 1000`: 11,5 s contra 0,65 s), com `λ.min` entre os pontos
22 e 43 de ~100. Na grade `2:8` com `n = 250`, é ali que o `glmnet` deixa
de convergir (6 avisos). E o `sglasso` herdou o `thresh = 1e-9` medido para
o `glmnet`: com `1e-8`, o padrão do `sparsegl`, ele é **2× mais rápido e
passa o KKT em 200 de 200 pontos** nos dois desenhos medidos (com `1e-7`,
metade falha). Forçar `type.gaussian`, trocar denso por esparso e o
custo de R em volta do motor (< 3%) não compram nada.

### 2026-09-28: P1 e P2 ratificadas (D34, D35)

O autor ratificou as duas propostas de E2.4 e o código já as segue: a
grade padrão de `J` de `cv.wafc()`, `wafc_tune()`, `klopp` e `oracle` passou
a `2:8` (`wafc_J_top` em `tune.R`), e a margem padrão do caso periódico a
`0` (`wafc_eps_periodic` em `design.R`). **672 testes passam**; o
`01-smoke.R` imprime `OK`.

- **O custo medido da grade nova**, cenário não homogêneo, 10 dobras, nesta
  máquina: `cv.wafc` de 0,5 s para **9,1 s** em `n = 250` e de 0,4 s para
  **3,9 s** em `n = 1000`; o `klopp`, que usa a mesma grade por ser o
  concorrente no mesmo desenho, de 0,5 s para 5,7 s e de 1,6 s para
  **20,2 s**. É o que torna a medição da pergunta 29 prioritária antes de
  repetir o piloto: a cauda do caminho é justamente o que cresce com `J`.
- **Números registrados que não se reproduzem mais com os padrões:** os de
  `02-tune.R` e `03-tune-decomp.R` (E2.3), `05-sondagem-aplicacao.R`
  (E6.1a) e `06-timing.R` (2026-09-21) foram produzidos com a grade
  `2:⌈log_2 n/2⌉` e, depois de E2.1b, a margem `0.05`; rodá-los hoje usa
  `2:8` e `0`. Reproduzi-los pede passar os dois argumentos. O `04-pilot.R`
  segue os padrões novos, que é o que a E2.5 quer, e a parte `jgrid` agora
  escreve a regra antiga por extenso, para continuar comparando as duas.
- **Aberto pela ratificação e decidido no mesmo dia** (pergunta 30, D36):
  o `bsgl` escolhia a dimensão entre 4, 8 e 16, que são as de `J ≤ 4`;
  agora entre `2^2` e `2^8`. Em `n = 1000` ele passa a escolher 64 e o erro
  de validação cruzada cai de 2,129 para 1,958, isto é, a grade curta o
  truncava como truncava o WAFC. O preço é ~60 s por ajuste em
  `n = 1000`, contra 0,9 s. **674 testes passam.**

### 2026-09-28: as acelerações da pergunta 29, medidas

`wafc/scripts/07-accel.R`, grade `2:8`, cinco células (as quatro do piloto
mais `uneven` com SNR 4), `n ∈ {250, 500, 1000}`, dobras e dados iguais
entre os braços; 20 réplicas no LASSO (5 na `mixed`) e 10 no sparse group
LASSO (3 na `mixed`); 14 núcleos, 10 min e 48 min. Critério declarado antes:
`J.min`, `λ.min` e `λ.1se` idênticos em **todas** as réplicas.

- **O prefixo é exato**, como a teoria do motor dizia: o `cvm` calculado
  num prefixo do caminho difere do caminho inteiro em no máximo `3.6e-15`
  (LASSO) e `0` (grupos).
- **Validação cruzada em etapas, LASSO: passa, mas ganha pouco.** Seleção
  idêntica em **255 de 255** réplicas; o tempo total cai só **1,5×**
  (4709 s para 3095 s), de 2,6× a 3,4× em `n = 250` e `500` no não
  homogêneo e no nulo, e **fica mais lento** em `n = 1000` em três células
  (suave 10,6 s para 12,6 s), porque o protótipo recalcula o prefixo a cada
  extensão e, em `J = 8`, o mínimo está no ponto ~55 e a busca vai a 80.
  A estimativa de 5 a 10× que eu tinha feito era da cauda em `J = 5` e `6`
  da grade antiga, não da grade `2:8`.
- **Validação cruzada em etapas, grupos: reprovada.** Idêntica em **104 de
  129**; muda a escolha em 18 de 30 réplicas do nulo e em 3 de 3 da `mixed`
  com `n = 1000`, porque nos `J` profundos o mínimo do caminho inteiro
  está no fim (pontos 89 a 100) e a busca para em 40.
- **`thresh = 1e-8` no `sparsegl`: reprovado pelo critério, por uma
  réplica.** Idêntico em **128 de 129** (o `λ.1se` muda numa réplica de
  `uneven` com `n = 500`), caminhos de `λ` iguais nas duas tolerâncias,
  KKT limpo no `λ.min` nas duas; ganho de 1,6×.
- **Achado que não era a pergunta, e que importa para E2.5: a validação
  cruzada do sparse group LASSO se engana no nulo com a grade funda.**
  Escolhe `J = 8` em 6 de 10 réplicas com `n = 500`, com o `λ.min` no fim
  do caminho. Conferido numa réplica: em `J = 8` o `cvm` **desce** até o
  fim (0,4178 para 0,398) num ajuste com 1223 de 1530 coeficientes não
  nulos, todas as dobras devolvem os 100 pontos (não é truncamento), e esse
  ajuste é **pior** contra a verdade (RMSE 0,221 contra 0,050 em `J = 2`).
  O LASSO, no mesmo dado, escolhe em `J = 8` um ajuste com 24 não nulos e
  RMSE 0,092, e o `cvm` dele sobe na cauda (0,614), como deve. A causa não
  está identificada (pergunta 31). **Corrigido pela investigação seguinte
  (bloco abaixo):** não é defeito nem vazamento, e o custo em erro
  verdadeiro é o mesmo do LASSO; a réplica examinada estava no percentil
  90, e "se engana" era forte demais.

### 2026-09-28: a pergunta 31, investigada

`wafc/scripts/08-sgl-null.R`, três partes, cada uma um teste de hipótese;
rodado duas vezes, com resultado idêntico. **Veredito: não há defeito no
código nem vazamento; é a maldição do vencedor numa faixa plana da curva
de validação cruzada, e o preço em predição é o mesmo do LASSO. O que o
sparse group LASSO paga a mais é estrutura, não erro.**

- **O laço de dobras está certo.** O `cv.sparsegl()` do próprio pacote, no
  mesmo desenho, caminho e dobras, dá a mesma curva: diferença máxima de
  `7e-9` em `J = 2` e `3e-6` em `J = 8`, e o mesmo `argmin`.
- **A mecânica.** No nulo, em `J = 8`, o caminho do sparse group LASSO fica
  em zero por ~85 dos 100 pontos e salta para 1223 não nulos nos últimos
  ~15, terminando em `λ = 0,01 λ_max`. A penalidade de grupo vale
  `0,95·√255·λ ≈ 15λ` por bloco; no nulo os gradientes de todos os blocos
  são ruído do mesmo tamanho, então eles entram juntos, densos e muito
  encolhidos (com `asparse = 0.05` a esparsidade dentro do bloco quase não
  age). Isso cria uma faixa de ajustes quase nulos cuja validação cruzada
  difere do nulo menos que o ruído dela, e o mínimo sobre 7 valores de `J`
  e ~20 pontos dessa faixa cai ali por acaso. O LASSO ativa um a um e sua
  curva sobe logo.
- **Não há vazamento.** O excesso de `cvm` sobre o ponto nulo se decompõe
  em `média((η̂ − f)²) − 2·média(ε(η̂ − f))`; sem vazamento o termo
  cruzado tem média zero. Em 30 réplicas do nulo (`n = 500`, `J = 8`):
  termo cruzado `+0,0005` (`t = 0,6`) no fim do caminho do sparse group
  LASSO e `+0,0019` (`t = 1,7`) no LASSO com o mesmo erro de treino; o fim
  do caminho fica abaixo do nulo em só 2 de 30.
- **O custo em erro verdadeiro**, excesso sobre a verdade em fração de
  `σ²` (mediana), `cv.wafc` na grade `2:8`, 20 réplicas por `n`:

  | nulo | LASSO `λ.min` | grupos `λ.min` | os dois `λ.1se` |
  |---|---|---|---|
  | `n = 250` | 0,026 | 0,025 | 0,013 |
  | `n = 500` | 0,013 | 0,012 | 0,005 |
  | `n = 1000` | 0,005 | 0,005 | 0,002 |

  Idêntico entre as variantes. Com `λ.1se` os dois voltam ao modelo nulo
  (mediana de 0 não nulos). No não homogêneo, `n = 500`, as duas também
  empatam (1,24 e 1,26).
- **O que difere é a estrutura:** no nulo, o sparse group LASSO em `λ.min`
  liga **303, 105 e 23** coeficientes (medianas, `n = 250, 500, 1000`),
  contra **3 a 5** do LASSO, e cai nos últimos 20 pontos do caminho em 13,
  7 e 2 de 20 réplicas. E nos dois o mínimo da validação cruzada fica
  abaixo do ponto nulo em 19 ou 20 de 20 réplicas, e o `J` escolhido no nulo
  é essencialmente aleatório: `λ.min` sempre escolhe *alguma* coisa quando
  não há nada.
- **O que isso muda para a E2.5:** a coluna `wafc.sglasso` do piloto pode
  ser lida em predição, inclusive no nulo; em estrutura, `λ.min` não é
  regra de seleção para nenhuma das duas (é o que D20 já dizia), e o
  sparse group LASSO fica pior que o LASSO nisso. Estrutura se lê pela
  limiarização de E1.7c ou por `λ.1se`.

### 2026-09-30: E5c fechada, o manuscrito em `k = 2`

Chat de tarefa, integrado aqui. Conferido nesta máquina: `ms_2.tex` e
`supp_2.tex` compilam com `latexmk` sem erro e **sem referência ou citação
indefinida**, em **33 e 30 páginas** (eram 28 e 26); `references_2.bib` tem
**37 entradas** (eram 30); a tarefa tocou só os arquivos da coluna do
catálogo. As perguntas 15 e 26 fecham.

- **"Sieve" (D33, marcado como D38):** 21 trocas no `ms` e 5 no `supp`, uma
  menção retida na definição de `W_J` ("approximation space (the sieve)");
  os rótulos `sec:sieve` e `thm:sieve` ficaram com o nome antigo, como em
  E1.10. Quatro escolhas fora do dicionário: "sieve oracle" virou "oracle of
  the approximation space" onde é definido e "oracle" no resto; "sieve
  coefficients" do resumo virou "wavelet coefficients"; "approximation error
  of the sieve" virou "of `W_J`"; o título da S3 virou "Approximation by the
  wavelet spaces `W_J`".
- **§4.3, a regra teórica, virou número:** a penalidade teórica é 7 a 13
  vezes a que minimiza o erro realizado e é ela, mais que a resolução, que
  torna a regra cara; `J_n` erra por até dois níveis, acima ou abaixo
  conforme `s'`. Ficou um comentário `% E5b:` para ancorar em tabela.
- **§4.2 alinhada a D23, D26 e D35:** o padrão é `ε = 0`, o cenário da
  Seção 3; a margem fixa aparece como opção de amostra finita, com o preço
  `ε^{−(s−1/π)}` do Lema 10 e as duas razões para ficar fora da teoria. Em
  prosa, sem enunciado nem prova no `supp`.
- **Software citado** (L4): lasso, `glmnet`, R, `WaveBased`, `sparsegl` e
  Johnstone, mais Huang, Horowitz & Wei (2010) junto do novo corolário.
- **Reprodutibilidade:** uma frase no fim do resumo, ainda sem endereço (D4
  mantém `wafc/` privado).
- **O Corolário 8 (D32) está no artigo:** nova §3.6, "Which modulators act
  on which coefficients", com a **Assumption 6** (separação por bloco) e o
  **Corollary 2** (triagem sem hipótese; recuperação sob a Assumption 6),
  dizendo que é estimação seguida de limiar. No `supp`, a **Seção S7** com o
  Lema 13 como **Lemma S7.1** e a prova. O parágrafo da variante estrutural
  da §4.3 remete a ele.
- **O mapa de numeração não existia no cabeçalho do `ms`**, ao contrário do
  que o `TAREFA.md` dizia; E5c o criou no cabeçalho do `ms_2.tex`, com a
  linha nova (Corolário 8 → Corollary 2, Hipótese S → Assumption 6, Lema 13 →
  Lemma S7.1). Nenhum resultado existente mudou de número.
- **Para D21:** o corpo cresceu 5 páginas (§3.6 ~2,5, §4.2 ~1, referências);
  a projeção de 39 a 41 páginas feita no fecho de E5a sobe na mesma medida,
  para ~44 a 46, acima do teto de 40.
- **Três lições de composição**, registradas no §3 do `instrucoes.md`:
  `\textcolor` em título de seção do `supp` quebra a compilação (o cabeçalho
  corrente põe o título em maiúsculas e a cor vira `COLR1`); `\textcolor`
  não atravessa parágrafo nem ambiente de teorema, e inserção longa vai num
  grupo `{\color{colR1} ... }`; cor dentro de `\citep{...}` não funciona, e
  citação acrescentada ao lado de outra vai por `\citetext` com `\citealp`.

As oito perguntas que a tarefa deixou estão na pergunta 32 da §4.

### 2026-09-30: E2.4c fechada, o piloto pronto para a repetição

Chat de tarefa, integrado aqui. Conferido nesta máquina: **674 testes
passam**, 0 falhas; a tarefa tocou só `wafc/scripts/04-pilot.R` e
`wafc/README.md` (o `competitors.R` já repassava `k` e `engine` ao `gam`).

- **O que mudou no `04-pilot.R`:** célula `uneven` (`p = 3`, `q = 2`, SNR 4,
  50 réplicas), acrescentada **por último** para as quatro células de E2.4
  manterem as sementes; método `gam.matched`, com `k = wafc_k_matched(u, J)`
  no `J.min` do `wafc.lasso` da mesma réplica e motor `bam`, ao lado do
  `gam` com `k = 10`; a coluna `J` na tabela de medianas; e `WAFC_TAG`
  como prefixo de todo `.rds` (padrão `e24`).
- **Conserto que não estava pedido:** a semente de cada réplica vinha da
  posição da célula e do `n` nas listas **restritas** pelos argumentos, de
  modo que rodar uma célula ou um `n` sozinho sorteava outros dados, contra
  o que o cabeçalho prometia. Agora vem das listas completas; as rodadas
  completas mantêm as sementes de antes.
- **Fumaça**, 1 réplica, `n = 250`, as cinco células: 50 linhas, **0
  falhas**, 123 s em 8 núcleos. No nulo, `gam` e `gam.matched` coincidem
  na quinta casa, porque com `select = TRUE` o REML leva os suavizadores à
  reta, onde `k` não importa.
- **O risco de E2.5a é o `gam.matched` no topo da grade**, medido em
  ajustes isolados: 0,1 a 1 s até `J = 6`, 5 s em `J = 7`, 72 s em `J = 8`
  com `n = 250`, **522 s** no não homogêneo com `n = 1000` e `J = 8`, e
  **mais de 15 min, com 3 a 4 GB de memória por processo**, na `mixed` com
  `n = 1000` (interrompido). Como o `sweep_part` roda as 15 réplicas da
  `mixed` num mesmo `n` juntas, 12 núcleos pediriam ~48 GB, contra 31 da
  máquina. Daí E2.5a rodar a `mixed` à parte, com poucos núcleos, o que o
  conserto das sementes permite sem mudar os dados.
- **O tempo do `gam.matched` é só o do spline**; a busca que escolheu `J` é
  contada na coluna do `wafc.lasso`, e a tabela de tempo tem de dizer isso.
- A `uneven` entra também nas partes `lambda` e `j1`; `margin` e `jgrid`
  continuam restritas por nome.

### 2026-09-30: E2.5a fechada, e o piloto devolve no-go para o LASSO

Chat de tarefa, integrado aqui. Conferido nesta máquina nos `.rds` de
`wafc/cache/e25a/`: **6 450 linhas, 0 falhas**, e as razões centrais da
tabela abaixo reproduzem. Duas chamadas, como o catálogo mandava: quatro
células em 12 núcleos (4 h 08 min de relógio, 49,1 h de processador) e a
`mixed` em 4 (1 h, 3,75 h). **O risco de memória de E2.4c não aconteceu**
(pico de 1,35 GB por processo), porque o WAFC não escolhe `J = 8` na
`mixed`.

Razões dentro da réplica contra o `wafc.lasso` (mediana; `< 1` é o
concorrente melhor), `rmse_f` / ISE, `n = 1000`:

| célula | `klopp` | `gam.matched` | `gam` `k = 10` | `vcbart` | `bsgl` | `aspline` |
|---|---|---|---|---|---|---|
| smooth | 0.955 / 0.876 | 0.777 / 0.590 | 0.716 / 0.497 | 1.388 / — | 1.114 / 1.255 | 0.835 / 0.688 |
| uneven | 0.946 / 0.879 | 0.791 / 0.626 | 2.766 / 7.477 | 1.262 / — | 1.117 / 1.278 | 0.866 / 0.772 |
| inhomogeneous | 0.956 / 0.919 | 1.000 / 1.034 | 1.434 / 2.024 | 1.086 / — | 1.137 / 1.305 | 1.256 / 1.596 |
| mixed | 1.082 / 1.174 | 0.924 / 0.844 | 1.254 / 1.558 | 0.972 / — | 1.016 / 1.037 | 1.130 / 1.266 |

- **(b), margem no não homogêneo: não há.** O melhor concorrente em todo
  `n` é o **`klopp`**, o block LASSO de K&P no mesmo desenho (0,947 a
  0,956 em `rmse_f`, vencendo em 92% a 96% das réplicas); contra o
  `gam.matched`, empate em `n = 1000` (1,000, 50% para cada lado) e derrota
  pequena em 250 e 500. Na `mixed` (`q = 4`, o eixo do artigo), o
  `gam.matched` vence em **15 de 15** em `n = 1000`. Tomando em cada
  réplica o melhor dos seis concorrentes, o `wafc.lasso` o vence em 0% a
  4% das réplicas de toda célula com componente.
- **(a), o fator no suave: acima de 1,5.** `ISE(wafc.lasso)/ISE(gam.matched)`
  é 1,67 a 1,69 no `smooth` e 1,50 a 1,60 no `uneven`, com o `gam.matched`
  vencendo em 90% a 100% das réplicas. **D30 não se confirma em dimensão
  casada:** no `uneven` o WAFC vence o `gam` de `k = 10` por 2,3 a 7,5 em
  ISE e perde do casado por 1,5 a 1,6; é a lição de E6.1a, agora na
  simulação.
- **O que mudou desde E2.4 foi o `klopp`, não o WAFC.** No não homogêneo
  em `n = 1000` os números absolutos se reproduzem (WAFC 0,841, `vcbart`
  0,910, `gam` 1,189), mas o `klopp` foi de 0,997 a **0,807**, porque D34
  deu a ele a mesma grade. **Lição de leitura:** o "7,8% melhor que o
  melhor concorrente" de E2.4 era o `jgrid`, o WAFC contra si mesmo; um
  padrão que corrige uma truncagem do WAFC vale para o concorrente no mesmo
  desenho, e o veredito só se lê depois de os dois rodarem com ele.
- **Por que o `klopp` ganha** (diagnóstico da tarefa, 10 réplicas,
  `n = 1000`, fora do catálogo, reproduzindo o piloto exatamente): **é o
  agrupamento em pedaços de `⌈log n⌉`**, não a padronização do `grpreg`
  (o `cv.wafc` padronizado e o `klopp` de pedaço unitário coincidem com o
  LASSO) nem a penalização dos níveis. **Com os níveis livres**
  (`wafc_fit_klopp(penalize.levels = FALSE)`, já no código) o ganho se
  mantém ou cresce: 0,940 / 0,882 no `smooth`, 0,913 / 0,840 no `uneven`,
  0,965 / 0,930 no não homogêneo, vencendo em 90% a 100%. Não medido no
  nulo nem na `mixed`. Mesmo assim o `klopp` não passa o fator do suave
  contra o `gam.matched` (1,45 a 1,61 no `smooth`; 1,31 a 1,38 no `uneven`),
  embora vença o `gam.matched` por 4% a 5% no não homogêneo em `n ≥ 500`.
- **O `J` escolhido** não bate mais no topo da grade nas células com
  componente, fora 10 de 50 réplicas do não homogêneo em `n = 1000`;
  no suave fica em 3, no `uneven` em 4 e 5, no não homogêneo sobe de 5 a 7.
- **Tempo:** o `wafc.sglasso` é 42% do processador da rodada (22,2 h); o
  `wafc.lasso` custa de 7 a 26 s por ajuste com a busca; o tempo do
  `gam.matched` **não inclui** a busca de `J`, que é a do `wafc.lasso`.
- **Estrutura, de passagem:** com `cv.min` o `wafc.lasso` liga todos os
  blocos nulos, como D20 já dizia; é assunto da limiarização (pergunta 11).
- **Dois itens de script, sem edição:** o `bsgl` não registra a dimensão
  escolhida (`J` sai `NA`, embora `extra$df` exista), e o cabeçalho do
  `04-pilot.R` ainda descreve o custo da `mixed` com os números da grade
  antiga.

### 2026-09-30: E1.11 fechada, a teoria com blocos transfere e melhora

Chat de tarefa, integrado aqui. `derivations/08a-sondagem-blocos.md`
(sondagem, sem numeração global, como a `06a`) e `check/08a-blocos.R`,
que imprime `OK` nesta máquina em 1 min 30 s. **Veredito: a objeção teórica
da pergunta 33(b) cai.**

- **O oráculo de E1.5 transfere sem cone e com as mesmas constantes:** a
  perfilagem não depende da penalidade, e o passo que dispensa o cone é
  Cauchy-Schwarz sobre o autovalor cheio de E1.4. Muda só a calibração,
  agora pela desigualdade de Hsu, Kakade & Zhang (2012) com união sobre os
  pedaços: com pesos 1, `λ ≍ σ sqrt((b_n + log|𝒢|)/n)`, que é o `δ̂` de K&P
  e o `λ` de Lounici et al. Conferido: cobertura 1,000 em 1 200 réplicas;
  as três cotas do "Teorema 1 em blocos" valem em 360 ajustes, com
  `γ̃ ≥ λ_min(Σ̂)` em todos. Não é o oráculo de Lounici et al. no lugar do de
  E1.5: é o de E1.5 com a calibração em blocos.
- **A compressibilidade vira o risco ideal por pedaços**, e a hipótese de
  Besov de E1.3 ainda o implica sem hipótese nova (o Lema 9 com um Hölder
  dentro do pedaço; o expoente é o do Lema 4 de K&P). Conferido em 864
  sequências, razão máxima 0,786.
- **A taxa ganha `(log n)^{2s'/(2s+1)}`:** fica
  `n^{−2s/(2s+1)}(log n)^{(2/π−1)_+/(2s+1)}` no lugar de
  `(log n/n)^{2s/(2s+1)}`, na mesma janela de `J_n`, e **sem logaritmo em
  `π ≥ 2`**, o que atinge a cota inferior de K&P para `q = 1` e `X ⊥ U`. O
  Teorema 2 perde o logaritmo; o Corolário 8 recebe `ρ_n = n^{−s'/(2s'+1)}`;
  E1.4 não muda. Custo de escrever: ~5 páginas e 2 a 3 dias.
- **O que não fecha é o `klopp` do piloto:** o `grpreg` usa pesos
  `sqrt(|G|)` por padrão, e com eles o `λ` da teoria é ditado pelos pedaços
  unitários do nível 0, o que devolve a cota à ordem do LASSO. A teoria
  cobre pesos 1, os de K&P. No melhor `λ` de uma grade os pesos mudam o erro
  em 10% a 20%, com sinal que depende da verdade (pergunta 34).
- **Nos `n` do piloto o ganho é de constante:** o custo por coordenada de
  um pedaço cheio é 0,35 do LASSO em `n = 500` e cai como `1/log n`; em
  `n ≤ 1000`, `b_n` e `log(|𝒢|/α)` são parecidos (6 a 7).
- **Posicionamento:** o Teorema 2 de K&P não cobre `p` fixo (a (A6) deles
  exige `p ≥ (s + s_0)(1 + log n)`) e controla o viés pela norma dual, o que
  lhes custa `r* ≥ 2` e `L + 1 ≥ n^{1/2}`; a arquitetura de E1.5 dispensa as
  três.
- **Lição de conferência:** neste desenho `X_1 ψ_{jk}(U_m)` e
  `X_2 ψ_{jk}(U_m)` têm correlação 0,65; uma verdade que repete a mesma
  wavelet em dois blocos mede o mau condicionamento, não a penalidade (a
  razão erro/cota variou 2,9 vezes em `n`). Vale para qualquer conferência
  futura de seleção ou de grupos.
- **Duas referências a verificar**, que a tarefa não podia pôr no `.bib`:
  Hsu, Kakade & Zhang (2012, *Electron. Commun. Probab.* 17, DOI
  10.1214/ECP.v17-2079; faltam número do artigo e páginas), a única de que a
  teoria em blocos depende; e Cai (1999, *Ann. Statist.* 27(3), DOI
  10.1214/aos/1018031262), não lida, citada `[VERIFICAR]` como origem dos
  blocos de tamanho `log n`. Vão para a frente curta da pergunta 32(e).
- **Símbolos novos, só se a variante for adotada** (passam pelo
  `notacao.md`): `𝒢`, `G`, `b_n`, `|𝒢|`, `w_G`, `‖θ‖_{𝒢,w}`, `𝒢_0`,
  `W(𝒢_0)`, `Ψ̃_G`, `R_𝒢(θ; η)`, `λ_n^𝒢`; a lista com as razões está na §3
  do handoff, reproduzida no `08a`.

### 2026-09-30: E2.5b fechada, o block LASSO com níveis livres medido

Chat de tarefa, integrado aqui. Conferido nesta máquina: **678 testes
passam** (eram 674); `e25b-joined.rds` tem 8 250 linhas, e as razões do
`klopp.free` contra o `wafc.lasso` em `n = 1000` reproduzem (0,926, 0,939,
0,959, 1,000, 0,916). O `wafc.lasso` refeito nas quatro células bate com
E2.5a nas 600 linhas, e as 15 primeiras réplicas da `mixed` batem nas 450,
**diferença 0 em toda coluna e todo método, inclusive o `vcbart`**. 0
falhas em 2 850 linhas. Para `WAFC_METHODS` não mover nada, o gerador é
restaurado ao estado deixado pelas dobras antes de cada método; é isso que
torna legítimo acrescentar um método a uma tabela antiga sem refazê-la, e
economizou ~40 h de processador aqui. O `bsgl` passou a registrar a
dimensão escolhida.

Mediana da razão dentro da réplica, `rmse_f` / ISE, `klopp.free` contra
os outros (`< 1` é o `klopp.free` melhor), `n = 250`, `500`, `1000`:

| célula | contra `wafc.lasso` | contra `gam.matched` |
|---|---|---|
| smooth | 0.934, 0.957, 0.926 / 0.864, 0.911, 0.856 | 1.215, 1.220, 1.233 / 1.492, 1.545, 1.590 |
| uneven | 0.940, 0.901, 0.939 / 0.882, 0.813, 0.873 | 1.133, 1.139, 1.169 / 1.299, 1.306, 1.358 |
| inhomogeneous | 0.951, 0.947, 0.959 / 0.898, 0.902, 0.922 | 0.994, 0.959, 0.943 / 0.987, 0.904, 0.870 |
| null | 0.996, 1.000, 1.000 / 0.560, 0.822, 0.910 | 0.999, 0.925, 0.975 / 1.183, 0.499, 0.727 |
| mixed (50 réplicas) | 0.944, 0.909, 0.916 / 0.893, 0.826, 0.842 | 1.013, 0.998, 0.995 / 1.023, 1.008, 1.004 |

- **O item (a) da pergunta 33 está respondido.** Onde há componente, o
  `klopp.free` vence o `wafc.lasso` em toda célula e todo `n`: 4% a 10% em
  `rmse_f`, 8% a 19% em ISE, em 78% a 100% das réplicas. No nulo empata
  (cara ou coroa: vence em 54%, 50% e 48%), e **a penalização dos níveis
  explica a derrota inteira do `klopp` ali** (1,31 a 1,43 em E2.5a).
- **Contra o `gam.matched`:** vence no não homogêneo a partir de `n = 500`
  (4% a 6% em `rmse_f`, 10% a 13% em ISE, 70% a 84% das réplicas), empata
  na `mixed`, e perde no suave. O fator de ISE fica em 1,49 a 1,59 no
  `smooth` e 1,30 a 1,36 no `uneven`; com o fator 1,5, passa no `uneven` e
  na `mixed` e reprova no `smooth` em `n ≥ 500`.
- **Contra o `klopp` fiel a K&P**, soltar os níveis compra pouco nas
  células `q = 2` (0,984 a 0,997 em `rmse_f`), fora a `mixed` em
  `n = 1000`, onde o `klopp` escolhe `J = 5` em 29 de 50 réplicas, numa
  distribuição bimodal sem nenhum `J = 6`, e perde por 0,832. A causa não
  está identificada.
- **A `mixed` em 50 réplicas confirma E2.5a:** o `gam.matched` vence o
  `wafc.lasso` em 49 de 50 réplicas em `n = 1000`; tomando o melhor dos
  sete concorrentes em cada réplica, o `wafc.lasso` o vence em 0% das
  réplicas nos três `n`, e, nas cinco células, em uma réplica de 600 com
  componente. Quinze réplicas moveram medianas em até 0,043.
- **Custo:** as quatro células com dois métodos, 42 min em 12 núcleos; a
  `mixed` inteira em 50 réplicas, 2 h 02 min em 8 núcleos (15,7 h de
  processador, pico de 1,49 GB por processo). O `klopp.free` custa ~1,26
  vez o `klopp`.
- **Lição de máquina:** esta tem **8 núcleos físicos** e 16 fios; acima de
  8 processos cada um fica ~25% mais lento, e a estimativa de custo tem de
  contar isso (a `mixed` estimada em 12,5 h levou 15,7 h). Com a grade
  `2:8`, a `mixed` custa por tarefa o mesmo que as outras células, e o
  orçamento de 15 réplicas deixou de ter razão de custo.
- **Pendências de E2.5b** (nenhuma bloqueia): a dimensão do `bsgl` nas
  células `q = 2` não foi medida (~50 min em 12 núcleos, rodando só ele);
  o padrão de 15 réplicas da `mixed` no script (50 daqui em diante?); e o
  `J` bimodal do `klopp` na `mixed`, que só importa se o `klopp` fiel
  continuar como concorrente.

### 2026-09-30: L5 fechada, a bibliografia curta

Chat de tarefa, integrado aqui. Conferido nesta máquina:
`referencias-verificadas.bib` com **69 entradas** (eram 64), 69 `\bibitem`
e zero erro ou aviso de BibTeX. Os seis PDFs de `refs/` foram lidos, e os
nomes seguem D39.

- **Entraram:** Hsu, Kakade & Zhang (2012, *ECP* 17, no. 52, pp. 1–6, do
  cabeçalho do PDF; o Crossref não guarda páginas), Cai (1999, 898–924),
  van de Geer, Bühlmann & Zhou (2011, 688–749), Mallat (2009, 3ª ed.) e
  Restrepo & Leaf (1997); `hardle1998wavelets` passou a "Gerard" e
  "Alexander", como na folha de rosto. As outras três entradas com
  Tsybakov ficam com "Alexandre B.", como nas fontes delas.
- **A pergunta 22 fecha:** a âncora da partição da unidade é **Mallat
  (2009, §7.5.1, Teorema 7.16 e p. 319)**, que cobre também a base
  periodizada com `j_0 = 0`; a prova está no Teorema 7.4 com `k = 0`.
  Härdle et al. só dão a identidade por dedução no cap. 8, não no 5.
  **Restrepo & Leaf não serve de âncora:** a propriedade (5) e o Apêndice
  II afirmam que as funções de escala de todos os níveis formam base
  ortonormal, o que é falso (`⟨φ̃_{0,0}, φ̃_{1,0}⟩ = 2^{−1/2}`); ficou no
  `.bib` com o comentário do que não sustenta.
- **Cai (1999) é a origem dos blocos de tamanho `log n`**, como E1.11 o
  cita, e o Teorema 4 dele dá o mesmo expoente de E1.11,
  `(log n)^{(2/π−1)/(2s+1)}`, com a ressalva `α ≥ 1/p` e no modelo de
  sequência; os blocos em limiarização em geral são de Hall,
  Kerkyacharian & Picard (1999), não verificado.
- **Hsu et al.:** o "Theorem 1" do arXiv que E1.11 leu é o **Teorema 2.1**
  do artigo publicado, que toma `A` quadrada; o `08a` usa `A` retangular, o
  que não muda a matemática (entra só `Σ = A'A`), mas a citação deve ser ao
  Teorema 2.1 com essa observação.
- **van de Geer et al. sustenta a frase de E1.7c** (Teorema 3.2, p. 698),
  com três diferenças que um referee pode apontar: eles reajustam por
  mínimos quadrados depois do limiar, limiarizam coordenadas e evitam
  beta-min.
- **Correções que ficaram fora da coluna de L5** (pergunta 35): a
  citação de Mallat no `01-identificabilidade.md`; o Teorema 2.1 de Hsu
  et al. e as páginas de Cai no `08a-sondagem-blocos.md`; as páginas de van
  de Geer et al. e a ressalva do reajuste no `06-selecao-limiar.tex`; e,
  na próxima rodada do manuscrito, Mallat na prova do Lemma 1 do `supp`.
- **Lições de ferramenta:** o Crossref do *ECP* e da *EJS* não guarda
  páginas, e perde iniciais do meio (Kakade); o cabeçalho do artigo é a
  fonte. O Project Euclid pede verificação humana também no navegador
  embutido. Na extração de texto do Mallat `=` e `≠` viram o mesmo glifo:
  fórmula se lê na página renderizada. Com o `chicago.bst` os prenomes saem
  abreviados, de modo que D39 aparece no `.bib`, não no PDF.

### 2026-09-30: E2.5d fechada, o `gam` autônomo confirma E2.5a

Chat de tarefa, integrado aqui. `wafc/scripts/09-gam-autonomo.R` reproduz
o sorteio da parte `competitors` sem editar o piloto e ajusta o `gam` com
REML (`bam`, `select = TRUE`) e `k` fixo sem olhar o WAFC, 64 e 128 por
moduladora. Conferido nesta máquina: o `gam.matched` refeito bate com
E2.5a nas 600 linhas das células `q = 2` com diferença 0 (a tarefa
conferiu as 645, com a `mixed`), e o fator do suave abaixo reproduz. 645
réplicas, 1 935 ajustes, 0 falhas, 1 h 18 min em 12 núcleos, pico de
1,3 GB por processo.

- **A pergunta 33(d) está respondida: herdar o `J` do WAFC não dava
  vantagem ao `gam`.** O autônomo fica a 1% a 2% do `gam.matched` em quase
  toda célula, e é melhor onde o `J` pequeno do WAFC apertava a base
  (`uneven` em `n = 250`: 0,91 / 0,82).
- **O veredito de E2.5a se mantém com um concorrente autônomo, e o fator
  do suave piora.** `ISE(wafc.lasso)/ISE(gam.k128)`: 1,58, 1,69 e 1,92 no
  `smooth` e 1,81, 1,71 e 1,59 no `uneven`, com o `gam` vencendo em 88% a
  100% das réplicas; para o `klopp.free`, 1,47 a 1,68 e 1,38 a 1,63. Com
  1,5, o `klopp.free` passa só no `smooth` em `n = 250` e no `uneven` em
  `n ≥ 500`.
- **O `klopp.free` é o único método que vence o `gam` autônomo em alguma
  célula:** no não homogêneo em `n ≥ 500` (2% e 5% em `rmse_f` contra o
  `gam.k128`, 5% e 13% em ISE). Empata na `mixed` e perde no suave por 16%
  a 21% em `rmse_f`. Contra o `wafc.lasso`, o `gam` autônomo vence no suave
  (88% a 100% das réplicas), na `mixed` e no não homogêneo em 250 e 500, e
  empata em `n = 1000`; no nulo o WAFC é melhor em 500 e 1000.
- **O motor não é a causa:** `gam` com REML exato, sem discretização, dá o
  mesmo número que o `bam` a 0,5% em 40 ajustes, a 20 a 90 vezes o custo.
- **`k = 128` não aperta nos `n` do piloto** (`edf`/coeficientes do maior
  suavizador ≤ 0,61); `k = 64` aperta no não homogêneo e na `mixed` em
  `n = 1000` (até 0,74). O `gam.matched` apertava no suave (0,78 a 0,92).
- **Achado em E2.5a: o `bam` com `k = 256` e `n = 1000` diverge sem
  erro** em 2 de 38 ajustes (réplicas 12 e 14 do não homogêneo em
  `n = 1000`: `rmse_f` de 451 e 910, `edf` na casa dos milhares). As
  tabelas de E2.5a e E2.5b são medianas e não mudam, mas qualquer leitura
  por média ou gráfico dos `.rds` tem de marcar as duas. **Lição:** "0
  falhas" no piloto quer dizer "nenhum erro", não "nenhum ajuste absurdo";
  concorrente de base grande pede a checagem `edf` ≤ coeficientes.
- **Custo:** o `gam.k64` custa menos que a busca do WAFC nas células
  `q = 2` (2,7 s contra 7 a 26 s); o `gam.k128` custa o mesmo ali e de 8 a
  44 vezes mais na `mixed` (~8 min por ajuste com 12 processos). O tempo
  medido em varredura paralela é cota superior.
- **Pendências** (pergunta 36): qual `gam` entra em E4, a `mixed` em 50
  réplicas para o `gam` autônomo, a guarda de `edf` no `wafc_fit_gam()` e
  a linha do `09-gam-autonomo.R` no `wafc/README.md`, que é arquivo de
  E2.5c e entra na integração dela.

### 2026-09-30: E2.5c fechada, os pesos do block LASSO

Chat de tarefa, integrado aqui. Conferido nesta máquina: **696 testes
passam** (eram 678); `e25c-joined.rds` tem 9 750 linhas, e as razões
contra o `klopp.free` em `n = 1000` reproduzem. O `klopp.free` refeito na
`smooth` bate com E2.5b em toda coluna fora o tempo: é a segunda junção
exata, e o procedimento está validado para acrescentar métodos. Rodada com
a máquina só para ela: 1 500 ajustes, 1 h 14 min em 8 processos, pico de
0,94 GB.

Três formas do block LASSO com níveis livres: `klopp.free` (pesos
`sqrt(|G|)` do `grpreg`), `klopp.unit` (pesos 1, os de K&P e os que E1.11
cobre) e `klopp.merged` (os níveis grossos de cada bloco num pedaço só,
pesos do `grpreg`). Razão contra o `klopp.free`, `rmse_f` / ISE, `n = 250`,
`500`, `1000`:

| célula | `klopp.unit` | `klopp.merged` |
|---|---|---|
| smooth | 1.023, 1.024, 1.046 / 1.052, 1.051, 1.102 | 1.046, 1.043, 1.033 / 1.111, 1.097, 1.068 |
| uneven | 1.036, 1.080, 1.052 / 1.085, 1.165, 1.104 | 1.005, 1.001, 1.003 / 1.006, 1.003, 1.005 |
| inhomogeneous | 1.038, 1.055, 1.045 / 1.066, 1.112, 1.088 | 0.988, 0.992, 0.998 / 0.973, 0.986, 0.994 |
| null | 1.002, 0.995, 0.994 / 0.844, 0.866, 0.826 | 1.004, 1.000, 1.000 / 0.936, 0.986, 0.966 |
| mixed | 1.028, 1.077, 1.063 / 1.065, 1.165, 1.116 | 0.980, 0.988, 0.993 / 0.956, 0.976, 0.985 |

- **A pergunta 34 está respondida em medição: a forma que a teoria cobre
  é a que pior prediz.** Os pesos 1 perdem para o `klopp.free` em toda
  célula com componente, por 2% a 8% em `rmse_f` e 5% a 17% em ISE, e
  escolhem `J` mais raso; no não homogêneo em `n ≥ 500` o ganho de E2.5b
  sobre o `wafc.lasso` some (1,003 e 1,002). O padrão do `grpreg` é o
  melhor no suave; a junção dos níveis grossos é o melhor no não homogêneo
  e na `mixed`, por até 2% (melhor das três em 83% e 86% das réplicas).
- **A escolha entre as formas é de segunda ordem** (1% a 8%) perto da de
  blocos contra LASSO (4% a 10%, E2.5b). Nenhuma forma passa o fator 1,5
  no `smooth` com `n ≥ 500`, nem contra o `gam.matched` nem contra o
  `gam.k128`.
- **Na `mixed`, o `klopp.merged` é o único método abaixo do `gam.matched`
  nos três `n`** (0,980 a 0,990 em `rmse_f`, perto do empate: 56% a 60% das
  réplicas); em `n = 1000`, medianas absolutas de `rmse_f`: `klopp.merged`
  0,854, `klopp.free` 0,863, `gam.matched` 0,863, `klopp.unit` 0,921,
  `wafc.lasso` 0,942.
- **O `J` bimodal do `klopp` na `mixed`** (pendência de E2.5b) não aparece
  em nenhuma forma com níveis livres, que escolhem `J = 6` a 8 e nunca 5:
  indício de que a causa é a penalização dos níveis.
- **Descompasso entre teoria e prática**, que a pergunta 33 tem de pesar:
  se os blocos forem adotados, ou a teoria passa a cobrir os pesos do
  `grpreg` (e, pelo argumento de E1.11, a cota volta à ordem do LASSO por
  causa dos pedaços unitários), ou o artigo roda uma forma que perde 2% a
  8%. A tarefa propõe, como conjectura não verificada, uma quarta forma
  que talvez reconcilie as duas: juntar os níveis grossos e absorver a
  sobra de cada nível no pedaço anterior, de modo que todo pedaço tenha
  tamanho entre `b_n` e `2 b_n` e os pesos `sqrt(|G|)` fiquem a menos de
  `sqrt(2)` dos pesos 1 reescalados. Custo de medir: ~40 min em 8 núcleos
  e uma linha em `wafc_kp_groups()`.

### 2026-10-01: E2.5e fechada, os pedaços balanceados reconciliam teoria e prática em parte

Chat de tarefa, integrado aqui. Conferido nesta máquina: **730 testes
passam** (eram 696); `08a-blocos.R` imprime `OK`; `e25e-joined.rds` tem
10 500 linhas, e as razões do `klopp.balanced` em `n = 1000` reproduzem. A
junção é exata pela terceira vez (`klopp.free` refeito na `smooth`). 750
ajustes, 0 falhas, 36 min em 8 processos, pico de 0,67 GB.

**A forma:** `klopp.balanced` junta num pedaço os níveis com `2^j < b_n`
de cada bloco e absorve a sobra de cada nível fino no pedaço anterior, de
modo que todo pedaço fino tenha entre `b_n` e `2b_n − 1` colunas; níveis
livres e pesos do `grpreg`.

Razão `klopp.balanced` / referência, `rmse_f` / ISE, `n = 250`, `500`,
`1000` (`< 1` é a forma balanceada melhor):

| célula | `klopp.free` | `wafc.lasso` | `gam.matched` |
|---|---|---|---|
| smooth | 1.048, 1.039, 1.020 / 1.111, 1.092, 1.042 | 1.003, 1.009, 0.958 / 0.999, 1.008, 0.901 | 1.267, 1.315, 1.278 / 1.678, 1.733, 1.643 |
| uneven | 0.983, 0.995, 0.995 / 0.967, 0.984, 0.988 | 0.918, 0.900, 0.942 / 0.844, 0.803, 0.895 | 1.097, 1.130, 1.168 / 1.249, 1.296, 1.357 |
| inhomogeneous | 0.977, 0.990, 1.000 / 0.955, 0.984, 1.000 | 0.923, 0.937, 0.960 / 0.865, 0.886, 0.925 | 0.969, 0.949, 0.941 / 0.935, 0.892, 0.866 |
| null | 1.001, 0.982, 0.980 / 0.966, 0.758, 0.743 | 0.972, 0.999, 0.964 / 0.528, 0.825, 0.522 | 0.993, 0.938, 0.929 / 0.879, 0.423, 0.697 |
| mixed (50) | 0.967, 0.982, 0.993 / 0.934, 0.963, 0.986 | 0.911, 0.891, 0.908 / 0.825, 0.791, 0.823 | 0.991, 0.980, 0.982 / 0.976, 0.969, 0.989 |

- **A teoria cobre a forma balanceada com os pesos do `grpreg`** (§11 do
  `08a`, Parte E da conferência): o argumento de E1.11 aceita pesos de
  razão limitada pagando `ρ²` só no termo de estimação, e aqui
  `ρ² ≤ (2b_n − 1)/(b_n − 1) ≤ 3`, isto é 1,67 em `n = 250` e 1,57 em
  `n = 500` e `1000`. A taxa não muda. O `klopp.free` (`ρ² = b_n`) e o
  `klopp.merged` (sobra de uma coluna) não são cobertos. Conferido:
  cobertura de 0,997 a 1,000 em seis esquemas de peso; as três cotas do
  Teorema 1 em blocos com `sqrt(|G|)` valem em 240 ajustes; o risco ideal
  vale com os pedaços balanceados em 1 008 sequências.
- **E prediz bem fora do `smooth`:** é a melhor das quatro formas no
  `uneven`, no não homogêneo em `n ≤ 500` e na `mixed` (0,5% a 3% sobre o
  `klopp.free` em `rmse_f`, 56% a 94% das réplicas); empata no não
  homogêneo em `n = 1000`. Na `mixed`, tem a menor mediana de `rmse_f` nos
  três `n` entre os métodos que não conhecem a estrutura, mas contra o
  `gam.matched` é empate (0,98 a 0,99, 56% a 62%). No não homogêneo vence o
  `gam` autônomo nos três `n` (0,953 a 0,972, 64% a 82%).
- **No `smooth` perde do `klopp.free` por 2% a 5%**, e isso é a junção dos
  níveis grossos, não a absorção da sobra: em `J = 3` o bloco inteiro é um
  pedaço de 7 e ela coincide com o `klopp.merged` (76%, 60% e 32% das
  réplicas). Contra o `wafc.lasso` empata no `smooth` em `n ≤ 500` e vence
  em 1000.
- **O fator do suave não passa:** 1,64 a 1,73 contra o `gam.matched` e
  1,64 a 1,79 contra o `gam.k128` no `smooth`; no `uneven`, 1,25 a 1,36 e
  1,38 a 1,53. Nenhuma das quatro formas passa 1,5 no `smooth` com
  `n ≥ 500`.
- **Um achado que não favorece a leitura óbvia:** nestes `n`, os pesos 1
  ficam mais perto dos de Lounici et al. que os `sqrt(|G|)`; a vantagem de
  predição dos `sqrt(|G|)` com `λ` por validação cruzada não sai da cota.
- **Conjectura da tarefa, não medida** (pergunta 37): deixar livres os
  níveis com `2^j < b_n`, que a forma balanceada penaliza juntos. Na
  teoria iriam para o bloco não penalizado, com `σ²p_0/n`,
  `p_0 ≍ pq b_n`, de ordem `log n/n`, menor que a taxa. Pode recuperar o
  `smooth` e custar no nulo.
- **Lições:** a junção por `WAFC_METHODS` com prova numa célula barata
  funcionou três vezes e está madura para E4; um vigia
  `until ! pgrep -f "<script>"` casa com a própria linha de comando e
  nunca termina, então vigiar por PID ou arquivo de fim.

### 2026-10-01: E2.5f fechada, os níveis grossos livres não pagam

Chat de tarefa, integrado aqui. Conferido nesta máquina: **807 testes
passam** (eram 730); `08a-blocos.R` imprime `OK`; `e25f-joined.rds` tem
11 250 linhas, e as razões contra o `klopp.balanced` reproduzem. Junção
exata pela quarta vez. 750 ajustes, 0 falhas, 46 min em 8 processos.

**A conjectura de E2.5e não se confirma.** O `klopp.freecoarse` (níveis
com `2^j < b_n` sem penalidade, pedaços finos balanceados) é a pior das
cinco formas em quase toda célula e perde do `klopp.balanced` em todas:
em `rmse_f`, 1,06 a 1,22 no `smooth`, 1,05 a 1,10 no `uneven`, 1,01 a 1,09
no não homogêneo, 1,11 a 1,45 na `mixed` e 2,0 a 2,4 no nulo.

- **O mecanismo:** soltar os níveis grossos melhora os blocos ativos no
  `smooth` (ISE ativo 26% menor em `n = 1000`), mas os mesmos níveis ficam
  livres nos blocos inativos, que nunca zeram; o ISE deles sobe 3 a 6
  vezes. No nulo a validação cruzada escolhe `J = 2` em 98% das réplicas
  (mínimos quadrados em 21 colunas) e o erro é exatamente `σ²p_0/n`.
- **A teoria cobre a forma e a taxa fica** (§12 do `08a`, Parte F): o
  Teorema 1 em blocos vale com os níveis grossos no bloco não penalizado,
  `p_0 = p + pq(2^{j*+1} − 1)`, como na Proposição 5 de E1.8; o `ρ²` dos
  pedaços restantes fica abaixo de 2. A cota favorece soltar porque compara
  o custo por bloco ativo; a medição mostra que o `σ²p_0/n` dos blocos
  inativos decide.
- **O fator do suave da forma nova é o pior:** 1,88 a 2,61. Nenhuma das
  cinco formas passa 1,5 no `smooth` com `n ≥ 500`.
- **Lições:** um método que deixa colunas livres em todo bloco se lê com o
  ISE separado em blocos ativos e inativos, não só em `rmse_f`; o `grpreg`
  não aceita problema sem grupo penalizado (a tarefa ajusta por mínimos
  quadrados quando nada sobra).

**Medido no chat principal, sobre os mesmos `.rds`: onde está a perda no
suave.** Separando o ISE do `klopp.balanced` por bloco no `smooth`, os
blocos inativos respondem por ~13% do ISE dele (medianas de 0,0124,
0,0069 e 0,0045), contra ~2% no `gam.matched`, que com `select = TRUE`
leva os inativos a quase zero (0,0030, 0,0005 e 0,0004). **Se os blocos
inativos fossem zerados, o fator do suave do `klopp.balanced` cairia de
1,64 a 1,73 para 1,41, 1,45 e 1,32** (`ISE ativo / ISE do gam.matched`,
mediana dentro da réplica), e no `uneven` de 1,25 a 1,36 para 1,10 a 1,21:
abaixo de 1,5 nos dois. Isso é o limite de um limiar perfeito, não uma
medição de limiar; mas aponta para o Corolário 8, que já está no artigo
(D32), como o passo que falta à comparação: estimação seguida de limiar,
que é também a calibração de `t_n` da pergunta 11.

### 2026-10-01: L6 fechada, as correções de L5 nas derivações

Chat de tarefa, integrado aqui, correndo ao lado da E2.5g sem arquivo em
comum. Conferido nesta máquina: o diff toca só texto de citação (nenhum
enunciado, hipótese, número de resultado nem `\label`), e o
`06-selecao-limiar.pdf` recompila em 7 páginas sem referência indefinida nem
aviso. A tarefa releu cada localização nos PDFs de `refs/`.

- **`01-identificabilidade.md`:** a partição da unidade passa a Mallat
  (2009, §7.5.1, p. 319), com a prova no Teorema 7.4 (p. 284) com `k = 0`;
  a base periodizada com `j_0 = 0` ganha o Teorema 7.16 (p. 318), com a
  nota de que o `j` de Mallat é o `−j` daqui. Daubechies (1992, §9.3) fica
  só para a periodização; Restrepo & Leaf (1997) registrado como âncora que
  não serve.
- **`08a-sondagem-blocos.md`:** Hsu, Kakade & Zhang como Teorema 2.1
  (hipótese (2.1), `μ = 0`, Observação 2.2), com a nota de que a `A`
  retangular não muda nada porque só `Σ = A'A` entra; o "§2.3" que o
  arquivo citava não existe no artigo publicado e saiu. Cai (1999) com o
  Teorema 4 (eq. 5.4, p. 908), a ressalva `s ≥ 1/π` e o modelo de
  sequência; o `[VERIFICAR]` dele sai.
- **`06-selecao-limiar.tex`:** van de Geer, Bühlmann & Zhou (2011) com
  688–749, eq. (1.3) p. 691 e Teorema 3.2 p. 698, e as três diferenças para
  o Corolário 8: eles reajustam por mínimos quadrados, limiarizam
  coordenadas e não pedem beta-min, onde a Hipótese S é condição desse tipo
  por bloco. É a frase que um referee pode pedir na §3.6 do `ms`.
- **Sobram quatro marcas de verificação**, fora de L5 (pergunta 35(e)).

### 2026-10-01: E2.5g fechada, o limiar do Corolário 8 passa o fator do suave

Chat de tarefa, integrado aqui. Conferido nesta máquina: **921 testes
passam** (eram 807), 0 falhas; `e25g-joined.rds` tem **24 750 linhas**
(as 11 250 de E2.5f mais 13 500 limiarizadas), 0 falhas; as 1 500 linhas
de `wafc.lasso` e `klopp.balanced` coincidem com as de E2.5f com diferença
0 em toda coluna fora o tempo (**quinta junção exata**); e as razões
abaixo marcadas reproduzem. Rodada: 750 ajustes, 1 h 59 min em 8
processos, 15,9 h de processador, pico de 0,86 GB.

**O estimador:** `wafc_threshold()` (`wafc/R/threshold.R`) zera os blocos
com `N̂_{ℓm} ≤ t` num ajuste do WAFC ou do `klopp`, sem reajuste (o
Corolário 8) ou com mínimos quadrados no suporte (`+ls`, van de Geer et
al.) ou no bloco inteiro (`+lsb`). Três regras: `max` (`t = c·max N̂`,
`c = 0,15` de E1.7c), `cv` (`t` absoluto, escolhido nas mesmas dobras no
`(J, λ)` já escolhido) e `oracle` (referência). Os candidatos são
exaustivos, porque o ajuste limiarizado é constante entre normas
consecutivas; em `t = 0` o erro de validação cruzada é o `cvm.min` a
`1e-10`.

- **O fator do suave passa com limiar, pela primeira vez.**
  `ISE / ISE(gam.matched)` no `smooth`, `n = 250, 500, 1000`: o
  `wafc.lasso` vai de 1,67, 1,68, 1,69 (reproduzido) a **1,31, 1,37, 1,36**
  com `+max` (reproduzido) e 1,34, 1,40, 1,48 com `+cv`; o `klopp.balanced`
  de 1,68, 1,73, 1,64 a 1,46, 1,45, 1,32 com `+max`. No `uneven`, 1,23 a
  1,38 para as duas formas. Contra o `gam.k128`, só o `+max` fica abaixo de
  1,5 nos três `n` (1,22, 1,37, 1,45 no `wafc.lasso`); o `+cv` chega a 1,68
  em `n = 1000`.
- **O que sobra é ISE ativo.** O limiar sem reajuste não toca os blocos
  mantidos (razão de ISE ativo 1,000, reproduzido): o resto do fator é
  1,38 a 1,41 no `wafc.lasso` e 1,48 a 1,50 no `klopp.balanced`, que é o
  teto de um limiar perfeito; `+max` e `+oracle` já estão nele no suave. Isso
  confirma a estimativa feita aqui sobre os `.rds` de E2.5f.
- **Contra o próprio ajuste**, em `rmse_f`: `+max` ganha 8% a 10% no
  `smooth` (100% das réplicas) e nada no não homogêneo; `+cv` ganha 2% a 6%
  no não homogêneo e na `mixed`.
- **Contra o `gam.matched`**, em `rmse_f`: o `klopp.balanced+cv` vence no
  não homogêneo (0,951, 0,924, **0,906**, reproduzido; 82% a 90% das
  réplicas) e na `mixed` (0,944 a 0,971, 74% a 84%); o `wafc.lasso+cv`
  vence no não homogêneo em `n = 1000` (0,937) e perde na `mixed` (1,044,
  1,051, 1,025, reproduzido). **Com limiar, o LASSO passa à frente dos
  blocos no `smooth` em `n ≤ 500`** (`klopp.balanced+cv / wafc.lasso+cv`
  1,025 e 1,046), e os blocos continuam à frente no não homogêneo, no
  `uneven` e na `mixed` (0,91 a 0,98).
- **Estrutura**, `P(Ŝ = S)`: sem limiar, 0 em toda célula com componente.
  Com `+max`, 0,86, 1, 1 no `smooth` (reproduzido) e 0,52 a 1 no `uneven`,
  mas **0,02 a 0,42 no não homogêneo e 0 a 0,68 na `mixed`**. O `+cv+ls` é a
  melhor regra que não vê a verdade: 0,88 a 1 no suave e no `uneven`, 0,52
  a 0,98 no não homogêneo, 0,44 a 0,98 na `mixed`. O `gam.matched` acerta
  0,18 a 0,34 no suave e 0 na `mixed`.
- **A pergunta 11, medida:** o `c = 0,15` de E1.7c serve ao suave e não ao
  não homogêneo, onde a validação cruzada leva `J` a 5, 6 e 7 e a folga
  entre o menor bloco ativo e o maior nulo cai a 0,26, 0,44 e 0,57 da maior
  norma (0,80 a 0,90 no suave); que a causa seja o `J` mais fundo é
  leitura, não medida. **`c = 0,4` acerta 0,88 a 1 em `n ≥ 500`** nas quatro
  células com componente; o platô do não homogêneo e da `mixed` é estreito e
  fecha dos dois lados. **Nenhuma regra relativa zera o nulo** (o maior
  bloco fica sempre); `t = 0,3` absoluto zera, mas depende da escala do
  sinal e fica como conjectura.
- **Os reajustes:** o `+ls` ajuda o WAFC no suave (0,87 a 0,90 do ajuste) e,
  com `cv`, na `mixed`, e piora o `klopp.balanced` fora do suave; com `max`
  ele é desastroso no nulo (2,2 a 3,1 em `rmse_f`), porque o maior bloco
  nulo é reajustado inteiro. O `+lsb` é instável fora do suave (até 2,7) e
  custa quase todo o acréscimo de processador da rodada (170 a 190 s por
  ajuste na `mixed` em `n = 1000`). Nenhuma forma com reajuste passa o
  fator contra o `gam.k128` em `n = 1000`.
- **Custo do limiar:** `max` e `oracle`, menos de 0,01 s; `cv`, 0,06 a 0,7 s
  com `q = 2` e 2,7 s (`wafc.lasso`) e 10 s (`klopp.balanced`) na `mixed` em
  `n = 1000`, contra 24 e 48 s do ajuste.
- **Lições:** o limiar sem reajuste não mexe no ISE ativo, e a tabela com
  ISE ativo e inativo separados (lição de E2.5f) mostra isso de uma vez; o
  reajuste no bloco inteiro é caro e instável com `q = 4`, e se entrar em
  E4 é só no suporte.

As quatro perguntas da tarefa estão na pergunta 38.

### 2026-10-01: L7 fechada, a numeração conferida nas fontes publicadas

Chat de tarefa, integrado aqui, ao lado da E2.5g. Conferido nesta máquina:
o diff só toca citações, e o `06-selecao-limiar.pdf` recompila em 7
páginas sem referência indefinida. Os cinco PDFs foram postos em `refs/`
pelo autor (`klopp2015`, `lounici2011`, `huang2010`, `donoho1994`,
`buhlmann2011`); o `03` e o `04` não mudaram (o `04` cita a eq. (1.8) e o
Lema 1 de K&P, que conferem no *Annals*, e o `03` não cita K&P).

- **Klopp & Pensky no *Annals* não é o arXiv a menos de revisão no
  Teorema 2:** as condições `r* ≥ 2` e `L + 1 ≥ n^{1/2}` viraram a família
  `L + 1 = n^ς`, `1/2 ≤ ς < 1`, `r* > (2ς)^{−1}` (eq. 3.13), da qual o par
  antigo é o caso `ς = 1/2`; a escolha de `δ` é a (3.14), com `λ = δ̂` e não
  `δ̂/2`; e os Lemas 2 a 4 foram para o suplemento, que não foi lido. A
  ordem `σ sqrt(log n/n)` e as conclusões do `08a` não mudam, e o
  posicionamento de E1.11 continua de pé (K&P precisam de alguma
  regularidade mínima e de `L` polinomial em `n`, e a arquitetura de E1.5
  não), mas **quem citar K&P no manuscrito cita a forma do *Annals***.
- **Lounici et al.:** numeração idêntica; o *Annals* omite a prova do
  Teorema 3.2, e a seção 8 é no modelo multitarefa.
- **Resolvidos:** a VisuShrink de Donoho & Johnstone (1994) na §12 do
  `08a` (Definição 2, p. 445, e §2.4, p. 440); o LASSO limiarizado em
  Bühlmann & van de Geer (2011), Seção 2.9 e Teorema 7.8, com `δ` por
  validação cruzada, que é a regra `cv` de E2.5g; a condição de separação
  de Huang, Horowitz & Wei (2010), a (A1), com a seleção no Teorema 4(i).
- **As ocorrências fora da coluna conferem** (`busca-novidade.md`,
  `plano-projeto.md`, `notacao.md` §7, `ms_2`), fora a frase de
  `busca-novidade.md` l. 45 com as condições do arXiv; e uma marca no
  `05-taxas.tex` não se resolve com número (pergunta 35(f)).

### 2026-10-01: L8 fechada, as referências do `mgcv` e da escolha da dimensão

Chat de tarefa, integrado aqui. Conferido nesta máquina:
`referencias-verificadas.bib` com **82 entradas** (eram 69), 82 `\bibitem`
e zero erro ou aviso de BibTeX com `plain` e com `chicago`; sete linhas
novas no `literatura.md`, numa seção "Splines penalizadas e escolha da
dimensão da base" e no fim da do LASSO. Entraram o `mgcv` 1.9-4 e o que o
`gam` do piloto usa (Wood 2017, 2003, 2011; Wood, Goude & Shaw 2015; Wood,
Li, Shaddick & Augustin 2017; Marra & Wood 2011), a escolha da dimensão
(Ruppert & Carroll 2000; Ruppert 2002; Kauermann & Opsomer 2011; Pya &
Wood 2016, `misc` do arXiv, sem versão publicada) e os graus de liberdade
do LASSO (Zou, Hastie & Tibshirani 2007; Tibshirani & Taylor 2012, com as
páginas tiradas do `journal_ref` das cópias que o IMS depositou no arXiv).

- **Kauermann & Opsomer (2011) usam ML, não REML:** maximizam em `K` a
  log-verossimilhança perfilada do modelo misto (eq. (3), pp. 226–227) e
  não mencionam REML. Quem diz que a busca equivale à REML é Pya & Wood
  (2016, p. 1). O catálogo da E2.5h foi corrigido: o `gam.reml` é o método
  3 de Pya & Wood (§2.3), e o ML de Kauermann & Opsomer é a variante fiel
  ao artigo deles.
- **Ruppert (2002) no modelo aditivo usa um `K` comum a todas as
  componentes, só até 40** (§6, p. 751, `K ∈ {5, 10, 20, 40}`). A grade da
  E2.5h, `k` por moduladora até 120, é extensão da dele, e tem de ser dita
  assim.
- **Pya & Wood (2016)** viram a busca por REML escolher `k` maior, com erro
  maior, no seno com `n = 100` (§3, p. 4): é o sinal a olhar no `gam.reml`.
- **Lições de ferramenta:** as páginas do *Annals* saem do `journal_ref` da
  cópia que o IMS deposita no arXiv, sem o Project Euclid; o Crossref do
  livro da CRC só tem a data do e-book (o ano vem da página de créditos);
  o DOI do CRAN tem `published-print` na primeira versão do pacote (2000 no
  `mgcv`), não na usada; o Crossref abrevia prenomes de Kauermann &
  Opsomer, outro caso de D39.

### 2026-10-01: a notação da seleção (pergunta 32(a), D40), e `k = 3`

No chat principal. O autor aprovou os quatro símbolos e pediu a `k = 3`. A
troca é uniformização de notação por lista fechada de padrões, aplicada por
script **só dentro do modo matemático** e só nos trechos da seleção, e
entra sem marcação (exceção do §3 do `instrucoes.md`). Os padrões, na
ordem: `\widehat{N}_{` → `\widehat{\nu}_{`; `N_{` → `\nu_{`;
`\widehat{S}` → `\widehat{\mathcal{S}}`; `\delta_n` → `\nu_{\min,n}`;
`\delta` → `\nu_{\min}`; `D_n` → `\bar{\Delta}_n`; `S` → `\mathcal{S}`;
`D` → `\bar{\Delta}`, as duas letras isoladas e fora de comando.

| arquivo e trecho | trocas |
|---|---|
| `ms_3.tex`, §3.6 | 31 (2, 5, 4, 7, 0, 0, 8, 5, na ordem acima) |
| `supp_3.tex`, Seção S7 | 101 (8, 13, 10, 12, 13, 7, 12, 26) |
| `06-selecao-limiar.tex`, pelas macros (`\Nhat`, `\Ntrue`, `\Shat`, `\sepmin` e a nova `\sepminn`) mais o modo matemático | 87 (`\sepmin_{n}` 18, `D_{n}` 10, `S_{n}` 3, `S` 21, `D` 35) |

Conferido: o diff do `ms` só toca os símbolos, mais o `\bibliography` e o
comentário que aponta o `supp`; `ms_3` e `supp_3` compilam em 33 e 30
páginas e o `06` em 7, sem referência indefinida. Fora do modo matemático
nada muda, o que protege "D28", "D3" e "Hipótese~S" no `06`. O `N_{ℓm}(ε)`
da S3 do `supp` (outra contagem) ficou fora, por estar fora do trecho, e
agora não colide mais. O `δ` de Bühlmann & van de Geer (2011) na lista de
referências do `06` é outro objeto e ficou.

### 2026-10-01: L9 fechada, Donoho & Johnstone (1998) é sobre corpos de Besov

Chat de tarefa, integrado aqui. Conferido nesta máquina: o diff do
`05-taxas.tex` só tira a entrada de Donoho & Johnstone (1994), que o corpo
não cita e cuja descrição era falsa; o PDF recompila em 10 páginas. A linha
de K&P do `busca-novidade.md` passou às condições do *Annals* (eq. 3.13,
p. 1283; Observação 2, p. 1284). `refs/donoho1998.pdf` entrou.

- **Donoho & Johnstone (1998) não trata de bolas weak-`ℓ_τ`:** o risco
  minimax está nos Teoremas 4 e 5 (§4.1, pp. 890–891), sobre corpos de
  Besov `Θ^α_{p,q}(C)` no modelo de sequência, `≍ C^{2(1−r)} ε^{2r}` com
  `r = 2α/(2α+1)`. Com `α = s` e `τ = (s + 1/2)^{−1}` isso é
  `C^τ ε^{2−τ}`, **a mesma ordem `λ_n^{2−τ}` que a Observação `rem:tres`
  usa**; o que não confere é a classe. A forma sobre bolas `ℓ_p` vem de
  Donoho & Johnstone (1994b, *PTRF* 99, 277–303), fora do `.bib`.
- A marca da entrada de 1998 no `05` **fica** até a redação ser decidida, e
  o `ms_3.tex` l. 577 diz a mesma coisa que a Observação: pergunta 35(h).

### 2026-10-02: E2.5h fechada, o `gam` sintonizado não muda o veredito

Chat de tarefa, integrado aqui. Conferido nesta máquina: **1 006 testes
passam** (eram 921); `e25h-joined.rds` tem **28 350 linhas**, 0 falhas, e
as 24 750 que vêm de E2.5g coincidem com diferença 0 (**sexta junção
exata**); as razões abaixo marcadas reproduzem. Rodada: 10 h 08 min em 8
processos, 79,6 h de processador, pico de 1,20 GB por processo, no escopo
de D41 (grade `{5, 10, 20, 40, 80}` comum, 50 réplicas, `gam.gcv` fora da
`mixed`).

- **O critério de `k` não muda o erro do `gam` com REML:** `gam.reml` e
  `gam.cv` dão razão **1,000 em toda célula e `n`** (reproduzido), embora
  escolham `k` diferentes, e ficam a 2% do `gam.k128` e do `gam.matched`.
  A exceção é o topo da grade: no não homogêneo e na `mixed` em
  `n = 1000`, `k = 80` fica 1,0% a 1,4% acima do `gam.k128` em `rmse_f`.
  **D34 dispara** (mais de 20% no topo) no não homogêneo, na `mixed`, no
  `uneven` e no nulo, mas só aperta de fato naquelas duas células em
  `n = 1000`; a grade não foi estendida.
- **O GCV é outro estimador:** o `gam.gcv` é o melhor `gam` medido no não
  homogêneo (0,93, 0,89 e 0,93 do `gam.matched`, reproduzido; vence em 78%
  a 96%), perde 3% a 4% no `smooth` e 21% a 51% no nulo, e deixa mais
  blocos nulos ligados (o subsuavizamento conhecido do GCV, leitura).
- **Com o mesmo critério e as mesmas dobras dos dois lados (`gam.cv`), o
  veredito de E2.5a a E2.5g se mantém.** Razões contra o `gam.cv`, `rmse_f`,
  `n = 250, 500, 1000`: o `wafc.lasso` perde em 1,25 a 1,38 no `smooth` e
  1,26 a 1,33 no `uneven`; o `klopp.balanced+cv` vence no não homogêneo
  (**0,964, 0,923, 0,895**, reproduzido; 80% a 94% das réplicas) e na
  `mixed` (**0,978, 0,964, 0,941**, reproduzido), e o `wafc.lasso+cv` vence
  no não homogêneo em `n ≥ 500` (0,981 e 0,931). **O fator do suave do
  `wafc.lasso+max`** (`ISE / ISE(gam.cv)`) é **1,27, 1,44, 1,49** no
  `smooth` (reproduzido) e 1,72, 1,50, 1,35 no `uneven`: no limite de 1,5
  no suave, acima dele no `uneven` em `n = 250`. Contra o `gam.reml`, 1,53
  no `smooth` em `n = 1000`.
- **O GCV do WAFC não ganha nada sobre o `cv.min`** (1,00 a 1,02 em
  `rmse_f`, 1,09 a 1,14 no nulo), vai mais fundo em `J`, custa um décimo
  (0,5 a 3,9 s contra 7 a 26 s), e a guarda `df ≥ n/2` decide em 16,5% das
  réplicas, quase todas em `n = 250` nas células fundas. Com o limiar, empata
  com o `cv.min` limiarizado. Contra o `gam.gcv` (GCV dos dois lados) perde
  no suave e no não homogêneo e vence no nulo.
- **Pergunta 39(a), conferida:** `rank(M X_A) = |A|` em 150 de 150
  sintonias refeitas, no `λ` do `wafc.gcv` e no ponto mais denso que a
  guarda admite (`|A|` até 496, menor valor singular relativo 0,051). Nesta
  amostra a contagem "não nulos + `p`" é o `df` de Tibshirani & Taylor
  (2012): condição verificada, não teorema.
- **Estrutura:** o `gam.reml` e o `gam.cv` acertam `P(Ŝ = S)` de 0,16 a
  0,42 no suave, 0,24 a 0,34 no não homogêneo e 0 na `mixed`; o `gam.gcv`,
  de 0 a 0,20. O limiar do WAFC continua sendo a única regra que recupera
  estrutura, com as ressalvas de E2.5g sobre o `c`.
- **Lições:** o escore que o `bam` com `discrete = TRUE` reporta não é
  comparável entre `k` (mudou a escolha em até 10% das réplicas de uma
  célula), e o código usa o REML exato recalculado; escore de motor
  aproximado não serve para comparar modelos diferentes sem conferir contra
  o exato. De passagem: em 10 de 150 sintonias o `glmnet` não convergiu num
  `λ` e devolveu o caminho até o anterior, sem medição de onde (pergunta
  41(d)).
- **Tempo** (mediana por ajuste, busca incluída, 8 processos): `gam.reml`
  3,4 a 8,6 s com `q = 2` e 77 a 86 s na `mixed`; `gam.cv` 29 a 64 s e
  780 a 800 s; `gam.gcv` 103 a 318 s; `wafc.gcv` 0,5 a 3,9 s; `wafc.lasso`
  6,8 a 26 s. O `k = 80` é 75% a 87% do custo da busca do `gam`.

### 2026-10-03: o teto do limiar e o `gam.gcv`, lidos no chat principal

Sobre o `e25h-joined.rds`, sem rodada nova (scripts e saídas em
`wafc/cache/e25h/main-rumo*.R` e `.txt`, não versionados; os números
abaixo são o registro). Mediana da razão dentro da réplica, `n = 250`,
`500`, `1000`, com a fração de réplicas vencidas entre parênteses quando
importa.

- **O limiar oráculo é o teto de qualquer regra de `t`** (ele escolhe `t`
  pelo erro contra a verdade), e acerta a estrutura quase sempre:
  `P(Ŝ = S)` de 0,94 a 1 nos dois ajustes em `n ≥ 500`, nas quatro células
  com componente e no nulo. Contra o `gam.cv`, `rmse_f`:

  | teto | não homogêneo | `mixed` | `smooth`, ISE |
  |---|---|---|---|
  | `wafc.lasso+oracle` | 1,010, 0,973, 0,931 | 1,041, 1,067, 1,014 | 1,26, 1,44, 1,49 |
  | `klopp.balanced+oracle` | 0,945, 0,916, 0,895 (88% a 100%) | 0,956, 0,955, 0,932 | 1,37, 1,50, 1,52 |

  **Nenhuma regra de `t` dá ao LASSO a perna (b) do critério de E2.5**, e
  a diferença de estrutura entre LASSO e blocos medida em E2.5g é da regra,
  não do ajuste. O que uma regra melhor ainda pode comprar nos blocos é a
  distância do `+cv` ao teto: 0,964 → 0,945 no não homogêneo e 0,978 →
  0,956 na `mixed` em `n = 250`; no `smooth`, ISE 1,43, 1,65, 1,72 → 1,37,
  1,50, 1,52; e a estrutura, de 0,40 a 0,88 para 0,94 a 1 em `n ≥ 500`.
- **Contra o `gam.gcv` a perna (b) some:** `klopp.balanced+cv` 1,017,
  1,047, 0,982 em `rmse_f` (38%, 14%, 60%), e o teto não passa de 1,005,
  1,044, 0,982; o `wafc.lasso+cv`, 1,088, 1,123, 1,018. **O `GCV.Cp` é o
  `method` padrão do `mgcv::gam()`** (o do `bam()` é `fREML`; conferido no
  1.9-4). No nulo é o inverso: o `gam.gcv` perde de 0,57 a 0,73 em `rmse_f`
  para os dois ajustes limiarizados.
- **Nenhuma forma passa as duas pernas contra o `gam.cv`:** o
  `klopp.balanced+cv` passa a (b) (0,964, 0,923, 0,895) e reprova a (a) (ISE
  1,43, 1,65, 1,72); o `wafc.lasso+max` passa a (a) (1,27, 1,44, 1,49) e não
  a (b) (1,048, 1,024, 0,973). A mais perto das duas é o
  `klopp.balanced+max`: (b) 0,963, 0,940, 0,917 e (a) 1,43, 1,50, 1,52, na
  linha.
- **A estrutura é a vantagem que não depende do critério do `gam`:**
  `wafc.lasso+cv+ls` acerta 0,92 a 0,98 no não homogêneo e na `mixed` em
  `n ≥ 500` e 0,88 a 1 no suave e no `uneven`; `gam.cv` e `gam.reml`, 0,24 a
  0,34 no não homogêneo e 0 na `mixed`.
- **Contra o melhor dos outros concorrentes na réplica** (`bsgl`,
  `aspline`, `vcbart`, `klopp` fiel, `wafc.sglasso`): o `klopp.balanced+cv`
  vence no não homogêneo (0,974, 0,973, 0,970; 76% a 92%) e na `mixed`
  (0,977, 0,925, 0,882); o `wafc.lasso+cv` perde no não homogêneo em
  `n ≤ 500` (1,047, 1,031).
- **Correção do que o chat principal disse em 2026-10-02:** `boundary =
  "interval"` não é uma opção barata a medir. O `wafc_design()` recusa a
  opção com erro, porque a reparametrização do bloco de escala que mantém a
  identificabilidade não foi feita (`01-identificabilidade.md` §5); o
  `WaveBased::wbasis()` já tem a base. Fica fora da E2.5j e continua sendo o
  teste que o plano manda fazer antes de mudar de rumo (pergunta 33(f)).

### 2026-10-03: E2.5j fechada, a regra de um erro-padrão chega ao teto onde o sinal é forte

Chat de tarefa, integrado aqui. Conferido nesta máquina: **1 133 testes
passam** (eram 1 006); `e25j-joined.rds` tem **37 350 linhas** (as 28 350
de E2.5h mais 9 000), 0 falhas, e as 28 350 antigas coincidem com
diferença 0 (**sétima junção exata**; a rodada refez `wafc.lasso`,
`klopp.balanced` e os seus `+max`, `+cv` e `+oracle`); as razões abaixo
marcadas reproduzem. Rodada: 750 tarefas, `NC = 10`, **1 h 01 min de
relógio**, 9,9 h de processador, pico de 0,77 GB por processo (15,2 GB na
máquina, swap 0): a estimativa de 1,4 a 1,8 h era folgada.

**O que entrou:** em `wafc_threshold()`, as regras `"cvrel"` (o `c`
relativo por validação cruzada nas mesmas dobras) e `"cv1se"` (o maior `t`
cujo erro de validação cruzada fica a um erro-padrão entre dobras do
mínimo, a regra do `lambda.1se`), e a porta `gate = "qut"` antes de
qualquer regra (`wafc_threshold_gate()`; o pivô do QUT saiu para
`wafc_qut_pivot()`, e o `wafc_lambda_qut()` ficou idêntico bit a bit); os
registros de convergência no `wafc()` (`conv`), no `cv.wafc()` (`conv` e
`conv.folds`, lidos por `wafc_cv_convergence()`) e no `wafc_fit_klopp()`
(`extra$conv`); no `04-pilot.R`, seis rótulos por base, o controle
`WAFC_THR_REFITS` e as tabelas laterais `-thr-gate.rds` e `-conv.rds`.
Nenhum padrão de motor mudou.

Mediana da razão dentro da réplica, `n = 250, 500, 1000`.

- **O `cvrel` não recupera nada:** repete o `+cv` em quase toda célula,
  com o `c` que o `+cv` já escolhia (medianas de 0,27 a 0,58, quartis de
  0,05 a 0,82); a escala relativa não corrige os falsos positivos do suave
  (`P(Ŝ = S)` do `klopp.balanced+cvrel` 0,74, 0,56, 0,52, reproduzido).
- **O `cv1se` chega ao teto no suave, no `uneven` e no nulo, nos dois
  ajustes e nos três `n`:** a mediana é a do `+oracle`, e `P(Ŝ = S) = 1`
  no suave e no `uneven` (no nulo, 1 no `klopp.balanced` e 0,96 a 1 no
  `wafc.lasso`), sem porta. Fator do suave (`ISE / gam.cv`, reproduzido):

  | ajuste | `smooth` | `uneven` |
  |---|---|---|
  | `wafc.lasso+cv1se` | 1,26, 1,44, 1,49 | 1,56, 1,49, 1,35 |
  | `klopp.balanced+cv1se` | 1,37, 1,50, 1,52 | 1,33, 1,26, 1,25 |

  Contra o `gam.reml`, no `smooth`: 1,26, 1,48, 1,53 e 1,37, 1,52, 1,52.
- **No não homogêneo e na `mixed` o `cv1se` só chega ao teto em
  `n = 1000`;** em `n ≤ 500` zera bloco ativo fraco (0,78 a 0,90 falsos
  negativos por réplica em `n = 250`) e perde do `+cv`. `rmse_f / gam.cv`
  do `klopp.balanced` (reproduzido):

  | regra | não homogêneo | `mixed` |
  |---|---|---|
  | `+cv1se` | 1,015, 0,930, 0,895 | 1,000, 0,987, 0,932 |
  | `+cv` | 0,964, 0,923, 0,895 | 0,978, 0,964, 0,941 |

  `P(Ŝ = S)` do `klopp.balanced+cv1se`: 0,20, 0,86, 1 no não homogêneo e
  0,30, 0,82, 1 na `mixed`. O `wafc.lasso+cv1se` dá 1,061, 1,017, 0,931 e
  1,073, 1,079, 1,014. Contra o `gam.gcv`, o `klopp.balanced+cv1se` dá
  1,084, 1,047, 0,982 no não homogêneo e 0,42 a 0,48 no nulo.
- **A porta do QUT** zera tudo em 80% a 90% das réplicas nulas, em 6% a
  16% das do não homogêneo e da `mixed` em `n = 250` (0 a 2% em 500, 0 em
  1000) e nunca no suave e no `uneven`; decide igual nos dois ajustes em 723
  de 750 réplicas. Com ela, o `+max`, o `+cv` e o `+cvrel` passam no nulo de
  0,92 a 1,04 para 0,70 a 0,81 em `rmse_f` contra o `gam.cv`, e `P(Ŝ = S)`
  de 0,08 a 0,28 para 0,80 a 0,92 (0,88, 0,92, 0,86 no
  `klopp.balanced+cv+qut`, reproduzido); fora do nulo o custo fica em
  `n = 250` (`klopp.balanced+cv`: 0,964 → 0,966 e 0,978 → 0,984). Ao
  `cv1se` não acrescenta nada. **O nível não é o nominal:** rejeita 12% e
  17% das réplicas nulas contra `α = 0,05`, até 31% onde a validação cruzada
  escolheu `J = 2` ou `3`, e nada em `J = 5`.
- **A pergunta 41(d) está respondida:** `glmnet` com `jerr ≠ 0` em 36 de
  5 250 ajustes `(réplica, J)` (0,7%), quase todos em `J = 7` ou `8` e
  `n ≥ 500`, nove no `J` escolhido, e em nenhum o `lambda.min` encosta no
  corte nem alguma dobra foi cortada; os 4 288 caminhos curtos com
  `jerr = 0` são a parada antecipada do próprio `glmnet`. `grpreg`: nenhum
  caminho esgotou o `max.iter`, com pico de 9 015 de 10 000 na `mixed` em
  `n = 250` e `J = 4`; 3 caminhos validados mais curtos, com o corte em
  `λ ≈ 1e-4` e o `lambda.min` em 0,05 a 0,08.
- **Custo:** as regras novas reaproveitam as dobras do `+cv` (menos de
  0,1 s); a porta, 0,04 s de mediana.
- **Lição:** um teste exato num `J` fixo deixa de ser exato quando o `J` é
  escolhido nos mesmos dados (12% a 17% de rejeição sob o nulo contra 5%);
  é leitura do mecanismo, não medida.

### 2026-10-03: E1.13 fechada, a taxa lenta em blocos ganha um logaritmo

Chat de tarefa, integrado aqui. Conferido nesta máquina: o
`check/08-blocos.R` (com a Parte VII nova), numa cópia congelada, imprime
**`OK` em 209 s**; o `08-blocos.tex` compila em **18 páginas** (eram 14),
sem aviso nem referência indefinida; o diff só acrescenta (três linhas do
cabeçalho e da lista de peças reescritas).

**Corolário 14** (`cor:lenta-blocos`, §5 do `08`): sob as Hipóteses 1 e 2
de E1.5, as de E1.3, a Hipótese B, `A` de posto `p` e
`pq2^{J_n} log(pq2^{J_n})/n → 0`, `‖f̂ − f‖_n² = O_p(λ_n^𝒢 ‖θ*‖_{𝒢,w} +
𝓑²_{J_n} + p/n)`, com `λ_n^𝒢 ‖θ*‖_{𝒢,w} ≤ ϱ̄ λ^𝒢_{n,1} ‖θ*‖_{𝒢,1}`; a
contagem sob Besov, com `ϑ_π = min(1 − 1/π, 1/2)`, dá `‖θ*‖_{𝒢,1} ≤
pqC_g{(1 − 2^{−s'})^{−1} + b_n^{−ϑ_π} Σ_{j*<j<J_n} 2^{j(1/2−s)}}`. (i)
`s > 1/2` e `s' > 1/4`: `(log n/n)^{1/2}`; (ii) `s < 1/2` e `s' > s/2`:
`n^{−2s'/(1−2s+4s')}(log n)^{(2/π−1)_+ 2s'/(1−2s+4s')}`, e **em `π ≥ 2`,
`n^{−2s/(2s+1)}` sem logaritmo**, a taxa do Teorema 4, sem condição de
desenho.

- **A leitura do chat principal estava certa só até `‖θ*‖_1`:** pela
  cadeia `‖θ*‖_{𝒢,1} ≤ ‖θ*‖_1` o corolário reproduz a Proposição 4 com a
  constante `ϱ̄`, e nem o `2b_n − 1` nem o `w_max` mudam o logaritmo. **Mas a
  cadeia joga fora um ganho:** nos níveis finos `‖θ*‖_{𝒢,1}` é menor que
  `‖θ*‖_1` por `b_n^{−ϑ_π}`, o que em (ii) cancela o `sqrt(b_n)` do tamanho
  do pedaço e tira o logaritmo inteiro em `π ≥ 2` (nada em `π = 1`). Em (i)
  nada muda.
- **A conferência (Parte VII):** a cota da contagem contra o supremo exato
  na bola de Besov em 12 pares `(s, π)`, razão máxima 0,999, e a inclinação
  em `b` é `ϑ_π` à terceira casa; os expoentes do logaritmo em `n` de `1e20`
  a `1e300` batem com a teoria nos blocos e no LASSO (por exemplo 0 contra
  0,375 em `(s, π) = (0,3; 2)`); a cota em ajustes (FISTA, `λ` exato, pesos 1
  e `sqrt(|G|)`, inclusive dois desenhos sem a condição de E1.4 e um com
  `d > n`), sem violação, folga mínima de 27 vezes.
- **O desenho só entra por cima:** a calibração em blocos pede
  `λ_max(Σ̂) ≤ Λ` (a metade superior da Proposição 3 de E1.4, sem `c_U` nem
  `κ_1`), além de `σ̂_max`; a Proposição 4 do LASSO só pedia as normas das
  colunas. Nenhuma cota inferior de `Σ̂` é usada.
- **Achado em E1.6:** a Proposição 4 tem duas condições implícitas. Com
  `s' ≤ 1/4` em (i) a escolha `J_n ≥ log₂n/(4s')` dá `2^{J_n} ≥ n`, e com
  `s' ≤ s/2` em (ii) o balanço dá `2^{J_n} ≳ n/log n`; nos dois casos o
  Lema 6 de E1.5 não se aplica. Só afeta `π < 2` (11 pares de uma grade de
  156). O Corolário 14 escreve as duas condições.
- **Lição:** transpor um enunciado por uma cota grosseira de norma
  reproduz a ordem antiga e esconde o ganho; a ordem se confere na classe,
  com o supremo exato, antes de dizer que "é a mesma".

### 2026-10-03: E3.2 fechada, gráficos e documentação

Chat de tarefa, integrado aqui. Conferido nesta máquina: a suíte padrão
passa com **1 202 expectativas em 48 s** (17 puladas) e a inteira com
**1 427 em 104 s**, 0 falhas; as 127 novas são 86 do `test-plot.R` e 41 do
`coef` limiarizado, e o total antigo não mudou. Fora do `plot.R` novo, o
diff em `wafc/R/` é roxygen e a correção do `wafc_thr_object()`;
`wafc/scripts/` intocado.

- **`plot.wafc()`** (`"components"`, um painel por bloco com `ν̂_{ℓm}` e a
  verdade deslocada para a média da estimada, porque o nível é convenção;
  `"path"`, `ν̂_{ℓm}` contra `log λ`) e **`plot.cv.wafc()`** (`"cv"`, uma
  curva por `J` com a faixa de um erro-padrão; `"threshold"`, o erro contra
  `t` com a faixa que o `cv1se` lê; `"components"` limiarizadas; `"path"`),
  em gráficos de base, sem dependência nova; cada um devolve o que desenhou.
- **Roxygen auditado por script** nas 49 funções documentadas, 0 problemas;
  38 de 38 exemplos executáveis rodam. As assinaturas do README batem com o
  código nas 10 entradas.
- **O defeito de D48 corrigido:** o `coef` do objeto limiarizado está na
  parametrização do `coef.wafc()`; os dois testes novos falham contra o
  código anterior (9 falhas) e passam com a correção.
- **O exemplo do README**, de uma sessão limpa, em 16 s: no não homogêneo
  com `n = 500`, `J = 6`, e o `cv1se` mantém exatamente os três blocos
  ativos (normas 0,615, 0,695, 0,848) e zera os três nulos (0,178, 0,167,
  0,051). Figuras em `wafc/man-figures/` (290 KB).
- **Sinal:** a suíte padrão está em 48 s contra o teto de 60 s; teste caro
  novo vai com `skip_slow()`.
- **Lições:** o eixo da curva de validação cruzada tem de parar no erro do
  ajuste nulo, ou os mínimos se achatam; a verdade só se compara à estimada
  a menos de uma constante.

### 2026-10-03: E3.1 fechada, o WAFC numa chamada só

Chat de tarefa, integrado aqui. Conferido nesta máquina: a suíte padrão
passa com **1 075 expectativas em 46 s** (17 puladas) e a inteira
(`WAFC_SLOW_TESTS=1`) com **1 300 em 103 s**, 0 falhas nas duas; a
**oitava junção exata** reproduz (as 3 000 linhas de `wafc.lasso`,
`klopp.balanced` e dos 18 rótulos limiarizados na célula `smooth` iguais
às do `e25j-joined.rds` em toda coluna fora o tempo, 0 erros); e o
`wafc.block` coincide com o `klopp.balanced` nas 150 réplicas (mesmo `J`,
`rmse_f` a 1,1e-16; o `+cv1se` com ISE igual).

- **A interface** (`wafc/README.md`, "Interface (E3.1, proposta para
  ratificação)"): `wafc(penalty = c("block", "lasso", "sglasso"))` e
  `cv.wafc()` com o block LASSO balanceado como padrão, pelo `grpreg` nos
  pedaços de `wafc_kp_groups(penalize.levels = FALSE, balanced = TRUE)`, e
  `block.size` (padrão `⌈log n⌉`); `cv.wafc(..., threshold = c("cv1se",
  "cv", "max", "none"))`, depois do `...` para não casar com `thresh`; o
  limiar fica em `$threshold`, e `coef`, `predict`, `wafc_functions` e
  `wafc_blocks` leem o ajuste limiarizado em `lambda.min` (`thresholded =
  NULL`). `wafc_threshold()` passou a `rule = "cv1se"`. No bloco,
  `intercept = FALSE` e `standardize = TRUE` são erro (o `grpreg` sempre
  ajusta intercepto e ortonormaliza), o `cvm` segue a convenção do
  `cv.grpreg` (média sobre as observações) e `thresh` é `1e-4`. O
  `wafc_tune()`, o `wafc_sigma()` e o `wafc_fit_oracle()` continuam no
  LASSO.
- **Compatibilidade (D43):** 85 chamadas dos testes e os scripts 02 a 08
  fixam o padrão antigo (`penalty = "lasso"`, `threshold = "none"`); o
  `04-pilot.R` ganhou o método `wafc.block` (linhas `wafc.block` e
  `wafc.block+cv1se`). Testes novos em `test-interface.R` (158
  expectativas) e para `wafc_k_matched` e `print.wafc_competitor`, que não
  tinham; `helper-slow.R` com `skip_slow()` em 17 testes.
- **A tolerância do `grpreg`:** em `1e-4` o KKT do bloco falha no fim do
  caminho (49 de 100 pontos, resíduo até 0,79 `λ`; em `1e-8`, 100 de 100),
  mas **não chega à escolha**: em 40 réplicas (`smooth` e não homogêneo,
  `n = 250` e `500`), `1e-8` dá o mesmo `J`, o mesmo `λ` e os mesmos blocos,
  `rmse_f` com razão de 0,9973 a 1,0002, a 4,1 vezes o custo (9,7 contra
  40,0 s por ajuste).
- **Lição:** o `grpreg` soma `1e-5` ao primeiro `λ` de um caminho dado
  quando há grupo livre, e o devolve assim; interpolar o caminho de uma
  dobra em `s = λ_1` mistura dois pontos, e por isso a leitura do bloco é
  pelo ponto mais próximo.
- **Defeito pequeno, não corrigido:** o `wafc_thr_object()` grava
  `coef = c(0, b)` depois de dobrar o intercepto e perde o nível da
  covariável constante; nenhum código lê o campo.

### 2026-10-03: E1.12 fechada, a teoria em blocos provada

Chat de tarefa, integrado aqui. Conferido nesta máquina: o
`check/08-blocos.R`, rodado numa cópia congelada, imprime **`OK` em 207 s**;
o `08-blocos.tex` compila em **14 páginas** e o `06-selecao-limiar.tex` em
**11**, sem referência indefinida nem aviso; o diff do `06` só acrescenta
(o texto anterior não mudou, fora uma linha de comentário do cabeçalho). O
`check` não chama `wafc/R/`, de modo que a E3.1 em paralelo não o afeta.

**Os resultados** (numeração global; mapa no `TAREFA.md` §5):

- `08-blocos.tex`: **Proposição 7** (a forma balanceada: pedaços finos entre
  `b` e `2b − 1`, o grosso com `2^{j*+1} − 1 ∈ [b − 1, 2b − 3]` colunas,
  `|𝒢| ≤ pq(1 + 2^J/b)`, `ρ² ≤ 3` com `sqrt(|G|)`); **Lema 14** (a
  calibração em blocos, Hsu, Kakade & Zhang 2012, Teorema 2.1, e união sobre
  os pedaços); **Teorema 3** (o oráculo sem cone em blocos, escrito por
  inteiro, com comparador arbitrário, as constantes de E1.5 e `W(𝒢_0)` no
  lugar de `s_0`); **Corolário 9** (a forma de risco ideal); **Lema 15** (o
  risco ideal por pedaços sob a Besov de E1.3, dois regimes de `π`, `π` até
  `∞`); **Hipótese B** (o regime em blocos); **Teorema 4** (a taxa do espaço
  de aproximação `n^{−2s'/(2s'+1)}`, sem logaritmo); **Corolário 10** (as
  componentes, `ρ_n^𝒢 = n^{−s'/(2s'+1)}`); **Corolário 11**
  (`n^{−2s/(2s+1)}(log n)^{(2/π−1)_+/(2s+1)}`, sem logaritmo em `π ≥ 2`, na
  janela do Corolário 5). O `Lema 4` (perfilagem) é reimpresso em blocos sem
  número novo.
- `06-selecao-limiar.tex`, adendo: **Lema 16** (o risco do limiarizado,
  determinístico: blocos mantidos não mudam, nulo zerado só reduz, ativo
  zerado custa até `2t² + 2‖ĝ − g‖²`), **Corolário 12** (o Corolário 8 em
  blocos, com `ρ_n^𝒢`) e **Corolário 13** (a taxa do limiarizado nos dois
  ajustes, com a passagem à predição pelo número de condição de `Σ̂`).

**A conferência:**

- Proposição 7 em `b = 2..200`, `J = 1..14`; Lema 15 em 2 280 sequências
  (19 pares `(s, π)`), razão máxima 0,823; os dois expoentes do Corolário 11
  no risco ideal com 60 níveis (`2s/(2s+1)` a `1e−3`; o do logaritmo é
  exatamente `(2/π − 1)/(2s+1)` em `π < 2` e entre −0,12 e 0 em `π ≥ 2`).
- Calibração ao longo de `J_n`, `n = 250` a `4 000`: cobertura 1,000.
- **O realizado já separa as duas taxas:** no caso denso (`π = 2`,
  `s = 1/2`, `n` de 250 a 8 000), a inclinação do erro é **1,027 contra a
  taxa sem logaritmo e 1,194 contra `(log n/n)^{1/2}`**, o que E1.6 não
  conseguia mostrar; cotas do Teorema 3 e do Corolário 10 sem violação.
- **A teoria cobre o estimador exato do código** (`klopp.balanced`, a
  `penalty = "block"` da E3.1): a variante branca da Observação 2 (o `grpreg`
  ortonormaliza cada pedaço), com calibração pivotal e `q_max W` no lugar de
  `W`, conferida em `n = 500` (cobertura 1,000).
- Lema 16 em 14 000 combinações sem violação, com o fator 2 atingido.

**Achados:**

- **Posicionamento:** o Teorema 2 de K&P pede `r* > (2ς)^{−1} > 1/2` (com
  `r* = s'` em `π < 2`), e o Corolário 11 pede só `s' > s/(2s+1)`, que é
  menor que 1/2. **O cenário não homogêneo de D27 (`s' = 1/2`) está no nosso
  enunciado e fora do deles.** E a cota inferior de K&P vale sobre a nossa
  classe em `q = 1` com `X ⊥ U`, o que torna o Corolário 11 ótimo em
  `π ≥ 2` nesse caso.
- **Para D45:** no LASSO, limiarizar no platô de acerto baixa o erro das
  componentes, mas pode subir o de predição na amostra (até +5% em
  `n = 500`), porque o bloco nulo correlacionado com um ativo compensa o
  encolhimento dele; o Corolário 13 cobre isso com o número de condição de
  `Σ̂`.
- **Pré-assintótico:** como em E1.6, o evento de E1.4 só vale bem acima dos
  `n` do piloto (`λ_max(Σ̂)` mediano de 2,0 a 2,7 contra `Λ = 1,458`; em
  `J = 4` fixo, só em `n = 64 000`).
- **Nenhum expoente do `08a` mudou ao escrever a prova;** mudaram
  constantes: `A^𝒢 = 2 + (1 − 2^{1−π(s+1/2)})^{−1}` (`π ≤ 2`) ou
  `2 + (1 − 2^{−2s})^{−1}`, o pedaço máximo `2b_n − 1` em `λ^𝒢_{n,1}`, e
  `(λ_n^𝒢)² W(𝒢_n)` no `Δ̄_n`.
- **Lição de conferência:** um script lido pelo `Rscript` enquanto é
  editado quebra no meio; rodar sempre uma cópia congelada.

**O que muda no manuscrito em `k = 4`** (revisão da tabela da §6 do `08a`):
§2.3 passa à penalidade em blocos e à Proposição 7, com o LASSO como opção;
o Theorem 2 passa ao Teorema 3 (e o Lema 14 ao supp); o Theorem 3 ao
Teorema 4; o Corollary 1 ao Corolário 10; o Lemma 3 ao Lema 15; o
Theorem 1 ao Corolário 11; o Corollary 2 e a Assumption 6 ao Corolário 12
com `ρ_n^𝒢`; a §3.6 ganha o Corolário 13; o §1 e o "Why not block LASSO?"
são reescritos; no supp, S5 (Teorema 3, Lema 14, variante branca), S6
(Lema 15, Teorema 4, Corolários 10 e 11) e S7 (mais Lema 16, Corolários 12
e 13). As Proposições 1 e 2 e os Lemas 1 e 2 do manuscrito não mudam. A
Proposition 3 (taxa lenta) passa ao Corolário 14 (E1.13), no corpo, com as
condições `s' > 1/4` em (i) e `s' > s/2` em (ii) e a frase de que em
`s < 1/2` e `π ≥ 2` a taxa é `n^{−2s/(2s+1)}`, sem logaritmo e sem condição
de desenho, dita também no §3 depois do Theorem `thm:main` (D49(c)); a
contagem de `‖θ*‖_{𝒢,1}` e a prova vão ao supp, em S6. O texto diz que a
única propriedade de `Σ̂` usada é a cota superior (D49(b)).

### Decisões tomadas

| # | Data | Decisão | Razão |
|---|---|---|---|
| D1 | 09-18 | O modelo é `Y = Σ_ℓ β_ℓ(U) X_ℓ + ε`, `β_ℓ(u) = c_ℓ + Σ_m g_{ℓm}(u_m)` (escrito aqui já na notação de D9, congelada depois), com `X_1 ≡ 1` permitido (o aditivo puro e o parcialmente linear aditivo são casos particulares) | pedido do autor |
| D2 | 09-18 | Base: wavelets ortonormais de suporte compacto do `WaveBased`; padrão periódico com `j0 = 0` e a função de escala constante descartada (identificabilidade no nível da base, como no `wall()`); `boundary = "interval"` como opção | herda o `wall()`; é o que faz a restrição `∫ g_{ℓm} = 0` sair de graça (D22) |
| D3 | 09-18 | Estimador base: LASSO sobre todos os coeficientes de wavelet, `c_j` não penalizados; sparse group LASSO por par `(j,k)` é a variante a medir em E2 | pedido do autor (LASSO); a variante em grupos é a candidata natural à seleção de estrutura |
| D4 | 09-18 | **O código do método vive na pasta `wafc/` deste repositório** (`R/`, `tests/`, `scripts/`); o `WaveBased` não recebe código por enquanto e é usado só como dependência para as bases (`wbasis()`, `wtable()`), confirmado pelo autor; **as funções criadas ficam privadas por enquanto** (neste repositório privado, sem pacote público, sem `install_github`); empacotar e publicar decide-se em E3.3, com o código testado | decisão do autor, contra a proposta de implementar dentro do `WaveBased` |
| D9 | 09-18 | **Índices:** wavelet `ψ_{jk}` (nível `j`, translação `k`); covariável linear `X_ℓ`, `ℓ = 1, …, p`; moduladora `U_m`, `m = 1, …, q` | decisão do autor: não mexer no padrão da base; a colisão sai das covariáveis |
| D10 | 09-18 | **Coeficientes em `θ`**, com `θ_{ℓm,jk}` (bloco antes da wavelet); `β_ℓ` fica sendo só o coeficiente funcional | `β` não pode ser função e vetor na mesma seção |
| D11 | 09-18 | **`U ∈ [0,1]^q` por hipótese populacional**, com densidade limitada longe de `0` e de `∞`; reescalonamento empírico só na seção de computação | teoria limpa em E1.3 a E1.6; o termo extra não vale o custo agora |
| D12 | 09-18 | Ordem das colunas de `Z`: não penalizados, depois blocos `(ℓ, m)` lexicográficos, dentro do bloco `j` e `k` crescentes | fixa a interface de `wafc_design()` em E2.1 |
| D13 | 09-18 | **A teoria assume o caso geral** `λ_min(E[XX' \| U]) ≥ κ_1 > 0` q.c.; `X ⊥ U` vira observação (a fatoração de Kronecker) | E1.4 fechou o caso geral com `λ_min(Σ) ≥ κ_1 c_U`, conferido a 10% da verdade; assumir independência custaria generalidade sem comprar constante |
| D14 | 09-18 | **D11 é sobre a densidade conjunta** de `U` em `[0,1]^q`, não sobre as marginais | `U_2 = U_1` tem marginais uniformes e `Σ_Ψ` singular; a prova de E1.4 usa `c_U` da conjunta, e é daí que sai a "não colinearidade entre moduladoras" |
| D15 | 09-19 | Padrões do `wafc_design()`, ratificados do handoff de E2.1: `filter.size = 8` (contra os 20 do `wall()`), nomes de coluna `x2:u1:psi3.5` e de bloco `x2:u1`, erro informativo em `j0 != 0` e `boundary = "interval"`, e cenários com `β_1` aditivo em duas moduladoras, `β_2` em uma e `β_ℓ` constante para `ℓ ≥ 3` | o filtro 8 é o das três conferências de `derivations/check/`; o erro em vez da implementação mantém a pergunta 6 aberta sem código morto; o cenário tem de exibir o termo cruzado, que L2 apontou como o eixo do artigo |
| D16 | 09-19 | **(Emendada por D47 em 2026-10-03: o enunciado principal passa ao Corolário 11 de E1.12.)** **O enunciado principal do artigo é o Corolário 5 de E1.6** (taxa `n^{−2s/(2s+1)}` a menos de logaritmos, sob a mesma hipótese de Besov de E1.3), com o Teorema 2 como a taxa do sieve e a Proposição 4 como enunciado incondicional | ratificada pelo autor em 09-19 junto com D18: sob o posicionamento escolhido, o Corolário 5 é o paralelo direto do Teorema 2 de Klopp & Pensky |
| D18 | 09-19 | **O artigo se apresenta como extensão de Klopp & Pensky (2015)** a coeficientes aditivos em várias moduladoras e a desenho dependente; o parágrafo de posicionamento está escrito e aprovado em `alvo-revista.md` §4, e E5a o usa como está | decisão do autor, escolhendo entre os dois parágrafos redigidos; é a leitura honesta da literatura e a que o referee provável reconhece, ao custo de expor a falta da cota inferior para `q ≥ 2` |
| D17 | 09-19 | **Interface de `wafc()`** (E2.2): o `λ` do objeto é o do objetivo, não o do motor; `intercept` resolvido por presença de covariável constante, com erro informativo nos casos ambíguos; `coef()` dobra o intercepto no nível e `predict()` usa os coeficientes crus; mínimo quadrado escrito no ponto nulo do caminho; `wafc_kkt()` e `wafc_blocks()` públicas | a escala de `λ` é o que liga o código à teoria de E1.5, e as outras quatro saem dos defeitos de motor medidos; ratificada pelo autor em 09-19 |
| D5 | 09-19 | **Alvo primário: *Statistica Sinica***; reserva: *Electronic Journal of Statistics* | ratificada pelo autor; a linhagem do modelo está lá (Xue & Yang 2006; Wei, Huang & Li 2011) e o teto de 40 páginas em espaço duplo é folgado para a estrutura de ~29 planejada. O manuscrito nasce no template da revista (`manuscript/ss-template/`), com provas no suplementar |
| D8 | 09-19 | **O método se chama WAFC**, *wavelet additive functional coefficients* | ratificada pelo autor; é a sigla do repositório, ecoa o WALL e cabe no título. As alternativas "WAVC" e "wavelet additive coefficient LASSO" ficam descartadas |
| D21 | 09-19 | **A decisão sobre o teto de páginas fica para o fim**: E5b escreve sem contar, e quando o corpo estiver completo mede-se o compilado e decide-se entre resumir mais e mandar conteúdo ao suplemento | decisão do autor; a saída mais barata (tabelas ao suplemento) é a que a revista prefere e não exige reescrever prosa, desde que **E5b escreva cada tabela num `\input{}` próprio**, o que torna a mudança de lugar uma linha |
| D22 | 09-19 | **Centralização de Lebesgue:** a restrição de identificabilidade é `∫_0^1 g_{ℓm}(u) du = 0`, e não `E[g_{ℓm}(U_m)] = 0` como o `notacao.md` escrevia; a versão centrada em `P` sai pelo deslocamento `c_ℓ ↦ c_ℓ + Σ_m E{g_{ℓm}(U_m)}` e fica como observação | é o que a base impõe de graça (Lema 1 de E1.2) e o que o estimador estima; K&P (A1) usam base ortonormal em Lebesgue sem centralização, o WALL adota Lebesgue, e Xue & Yang centralizam em `P` mas recentralizam **empiricamente** na (3.4) deles. Adotar `P` obrigaria a mudar o alvo do Corolário 4 de E1.6, que hoje mede `‖ĝ − g‖_{L_2}` |
| D23 | 09-19 | **A margem `eps` do reescalonamento tem função declarada:** periodizar a base num intervalo maior que o suporte dos dados, o que permite trocar `g` por uma extensão que emende em `0 ≡ 1`. Com isso a hipótese não é "densidade limitada por baixo em todo `[0,1]`" e sim sobre o **suporte**, e a constante de E1.4 passa a ser `c_U λ_min(G_eps)`, com `G_eps` a Gram da base restrita ao suporte. `rescale = TRUE` continua padrão e `boundary = "interval"` entra no `wafc()`, mas fica fora do manuscrito por enquanto | argumento do autor em 2026-09-19: tomar `[0,1]` é sem perda de generalidade sobre a escala, não sobre o suporte. O ganho é que a Proposição 2 de E1.3 (periodização trava a taxa em `2^{−J/2}`) deixa de se aplicar quando `eps > 0`; a emenda é E1.3b, e E2.4 mede `λ_min(G_eps)` e o ISE contra `eps` |
| D24 | 09-19 | **As exigências de estilo da revista valem também nas derivações:** `E(·)`, `P(·)` e `Var(·)` em romano e com parênteses (§5 das instruções da SS), e um único `B_X` no lugar de `C_X` (E1.3) e `B_X` (E1.4). O `macros.tex` é ajustado quando E1.3b fechar | decisão do autor; manter dois conjuntos de símbolos para as mesmas quantidades é o que a notação congelada existe para evitar, e a revista não é negociável no ponto |
| D33 | 09-21 | **O termo "sieve" sai do manuscrito e das derivações**, trocado por "espaço de aproximação" / "approximation space" (e "sieve linear" por "aproximação linear", para fazer par com "aproximação não linear"); uma menção fica, entre parênteses, onde o espaço é definido. Cumprida nas derivações por E1.10 em 2026-09-21; no manuscrito, entra em `k = 2` | pedido do autor: é padrão em estatística teórica e em econometria semiparamétrica, e **incomum onde o artigo quer se inserir** — conferido, Klopp & Pensky (2015) e Xue & Yang (2006) não o usam nenhuma vez, e "Approximation space" é **palavra-chave** do artigo de Xue & Yang; o WALL teórico usa 38 vezes, e é de lá que ele entrou no nosso vocabulário |
| D34 | 09-28 | **A grade padrão de `J` é `2:8`**, independente de `n` (proposta P1 de E2.4), para `cv.wafc()`, `wafc_tune()` e os concorrentes no mesmo desenho (`klopp`, `oracle`); `n` maior que o medido pede `J` explícito | a grade `2:⌈log_2 n/2⌉` herdada do `cv.wall()` tinha o oráculo no topo em 100% das réplicas com componente; ir a 8 cortou 17,6% do erro de predição e 32% do ISE no não homogêneo em `n = 1000`, em 50 de 50 réplicas, e inverteu o veredito contra o VCBART (0.839 contra 0.910); ratificada pelo autor |
| D35 | 09-28 | **A margem padrão do reescalonamento é `eps = 0`** também na base periódica (proposta P2 de E2.4); margem positiva continua disponível pelo argumento | em E2.4, `0.05` custou 12% a 15% no suave, perdendo em 95% a 100% das réplicas, e não ganhou no não homogêneo (0.992 e 0.994); `eps = 0` perde no máximo 2,4% ali e é o único valor com `λ_min(G_eps) = 1`, onde a teoria do manuscrito é enunciada (D26); ratificada pelo autor |
| D36 | 09-28 | **O `bsgl` escolhe a dimensão por bloco entre `2^J`, `J = 2, …, 8`**, a grade de D34, em vez de 4, 8 e 16; candidato maior que o número de valores distintos de uma moduladora é descartado | com 4, 8 e 16 (isto é, `J ≤ 4`) a comparação que isola a base voltava a medir dimensão, a lição de E6.1a; medido no não homogêneo em `n = 1000`: a grade estendida escolhe 64 e baixa o erro de validação cruzada de 2,129 para 1,958, ao custo de ~60 s contra 0,9 s; ratificada pelo autor (pergunta 30) |
| D37 | 09-28 | **Nenhuma das duas acelerações da pergunta 29 entra:** o `cv.wafc()` continua calculando o caminho inteiro nas dobras, e o `sparsegl` continua em `thresh = 1e-9` | medido em `07-accel.R`: a validação cruzada em etapas passa no LASSO (255 de 255) mas ganha só 1,5× e fica mais lenta em `n = 1000`, o que não paga truncar o `cvm` devolvido; é reprovada nos grupos (104 de 129); o `1e-8` nos grupos é reprovado pelo critério declarado antes por uma réplica (128 de 129); decisão do autor |
| D38 | 09-30 | **Marcação da troca de "sieve" (D33) no manuscrito:** o termo que substitui vai em `colR1` e o termo antigo é apagado, sem `\sout` nem cinza; vale só para essa troca, e o resto de `k = 2` segue a marcação do §3 do `instrucoes.md` | decisão do autor: são ~31 ocorrências de uma troca de terminologia já ratificada, e o tachado em cada uma polui o PDF sem informar nada que o azul não diga |
| D39 | 09-30 | **Nomes de autor são grafados como estão impressos em cada referência**, sem uniformizar entre entradas: o mesmo autor pode aparecer como "Alexander" numa e "Alexandre" noutra (Tsybakov). A fonte é a folha de rosto ou o cabeçalho do próprio trabalho; o Crossref só na falta dela | decisão do autor: a citação reproduz o trabalho citado, e as grafias variam entre as publicações do mesmo autor |
| D40 | 10-01 | **Notação da seleção de estrutura:** `ν_{ℓm}` e `ν̂_{ℓm}` para as normas dos blocos, `𝒮` e `𝒮̂(t)` para a estrutura, `ν_min` e `ν_{min,n}` para a separação, `Δ̄` e `Δ̄_n` para a cota do erro por bloco (`notacao.md` §9); aplicada no `06` e na `k = 3` sem marcação, por lista fechada de padrões | as quatro antigas colidiam com `N` e `N_J`, com o suporte `S`, com o vetor `δ` e com `D = p + d` (pergunta 32(a)); aprovada pelo autor, que pediu a `k = 3` |
| D41 | 10-01 | **Escopo da E2.5h:** uma grade só, `k ∈ {5, 10, 20, 40, 80}`, comum aos suavizadores, para as três buscas do `gam` (REML, GCV e validação cruzada nas dobras do WAFC); 50 réplicas em todas as células; o `gam.gcv` fora da `mixed`, com a razão na tabela; a E2.5i volta a ser parte da E2.5h | grades diferentes entre os critérios do `gam` deixariam a comparação aberta à pergunta "por que esta grade para este critério?" (objeção do autor); Ruppert (2002, §6) usa `K` comum até 40 no aditivo, e E2.5d viu `k = 64` apertar; o 120 custava de 3 a 4 vezes o resto; o `gam.gcv` na `mixed` custa mais de uma hora por ajuste em `k = 80`, sem discretização, e o GCV mostrou mínimo local ali (medição de uma réplica da E2.5h); a exceção inteira, um método ausente numa célula, é mais limpa que uma grade cortada só ali; decisão do autor |
| D42 | 10-01 | **Quatro convenções bibliográficas de L1** (pergunta 28): citar Amato et al. (2022) e Haris, Simon & Shojaie (2018), que são trabalhos distintos; Hastie & Tibshirani (1993) com as páginas 757–779, sem a discussão; de Daubechies & Lagarias, só a parte I (1991), a do algoritmo de avaliação, e a parte II fica no `.bib` verificado sem ir ao manuscrito; manter as chaves herdadas do WALL (`cohen1993wavelets` e afins) | propostas de L1, ratificadas pelo autor; o `.bib` e o `references_3.bib` já seguiam as quatro |
| D51 | 10-03 | **Duas convenções da `k = 4`:** (a) a teoria do LASSO coordenado sai do corpo e **fica numa seção do supp**, com os enunciados e as provas que já estão no `supp_3`, como a teoria da opção `penalty = "lasso"` (D43), com a Proposição 4 emendada por D49; (b) **remoção de enunciado, prova ou matemática deslocada inteira** se marca com o texto antigo num grupo `{\color{gray} ...}` aberto por "[removed in k = r]", sem tachado, e o novo logo depois em `{\color{colR1} ...}`; trocas pontuais seguem com `\sout` (`instrucoes.md` §3) | decisão do autor, pelas recomendações do chat principal: o supp não conta no teto de 40 páginas, e o `\sout` não atravessa matemática deslocada nem enunciado |
| D50 | 10-03 | **O parágrafo de posicionamento e as contribuições em blocos** (`alvo-revista.md` §4), que substituem o terceiro e o quarto parágrafos da Introduction na `k = 4`: a extensão de K&P com o block LASSO de níveis livres e pesos de razão limitada; a taxa deles, sem logaritmo em `π ≥ 2`, ótima com uma moduladora independente, sob condições mais fracas (`p` fixo, `s' > s/(2s+1)`); **a frase "A lower bound for several modulators remains open." fica na introdução**; **três contribuições** (teoria, estrutura, computação), com a antiga "compressibility costs no assumption" absorvida na primeira e a numérica para a E5b; "Section~5" e o rótulo do Corolário 13 resolvidos na `k = 4` e na E5b | rascunho do chat principal, aprovado pelo autor com as três recomendações; D18 continua, e o texto antigo fica no histórico do git |
| D49 | 10-03 | **Pendências de E1.13 e E3.2 (pergunta 44):** (a) a Proposição 4 do `05-taxas.tex` ganhou as duas condições que a prova usa, `s' > 1/4` em (i) e `s' > s/2` em (ii), no enunciado e na prova (só afetam `π < 2`; o PDF segue em 10 páginas); (b) no manuscrito, o enunciado sem condição de desenho diz que nenhuma cota inferior de `Σ̂` é usada e que a única propriedade de `Σ̂` é `λ_max(Σ̂) ≤ Λ`, sem `c_U` nem `κ_1`, e a variante sem essa cota não é escrita; (c) o Corolário 14(ii) em `π ≥ 2` sobe ao §3 numa frase depois do Theorem `thm:main`; (d) `‖θ‖_{𝒢,1}` e `ϑ_π` no `notacao.md` §10; (e) as assinaturas de `plot.wafc()` e `plot.cv.wafc()` ratificadas; (f) as figuras ficam versionadas em `wafc/man-figures/`; (g) as funções sem roxygen fora de "Internals" vão a E3.3, e a escolha dos blocos desenhados a E6.2 | recomendações do chat principal, aceitas pelo autor |
| D48 | 10-03 | **Interface de E3.1 ratificada (pergunta 43):** `wafc()` e `cv.wafc()` com `penalty = "block"` padrão (o block LASSO balanceado pelo `grpreg`, `block.size = ⌈log n⌉`), `cv.wafc(..., threshold = "cv1se")` padrão, `thresholded = NULL` nos métodos de `cv.wafc`, o `cvm` do bloco na convenção do `cv.grpreg`, `thresh = 1e-4` no bloco, e `wafc_tune()`, `wafc_sigma()` e `wafc_fit_oracle()` no LASSO; o `grpreg` passa a `wafc_depends` (`load.R`, `CONTINUAR.md`); o defeito do `wafc_thr_object()` vai a E3.2 | recomendações do chat principal, aceitas pelo autor: a interface é o estimador medido, reproduzido à última casa na oitava junção exata; `1e-8` não muda a escolha e custa 4,1 vezes |
| D47 | 10-03 | **Pendências de E1.12 (pergunta 42):** (a) **o enunciado principal do artigo passa a ser o Corolário 11** (`08-blocos.tex`), no lugar do Corolário 5, o que emenda D16; (b) a razão dos pesos é **`ϱ`** e a soma dos pesos ao quadrado é **`𝒲(·)`**, trocados no `08` (19 `\rho`, 12 `\bar\rho`, 24 `\bar W`, 9 `W(`) e no `06` (3 `\rho` da razão, 9 `W(`) só no modo matemático e sem tocar a taxa `ρ_n` nem as macros `\rhoG`, `\rhoGt`; a lista de símbolos entrou no `notacao.md` §10, com as macros locais como no `07`; (c) **a taxa lenta (Proposition 3) fica no corpo, em versão de blocos**, como corolário do Teorema 3(i), a escrever no `08` (Corolário 14) e conferir; (d) a variante branca (Observação 2 do `08`) vai ao supp como o enunciado que cobre o estimador exato do código | recomendações do chat principal, aceitas pelo autor; a (c) foi revista antes da aceitação, porque a Proposition 3 é o único enunciado sem condição de desenho (D16) |
| D44 | 10-03 | **Rumo de E2.5 (pergunta 33): go reposicionado, na *Statistica Sinica*.** O WAFC é o **block LASSO na forma balanceada** (níveis livres, pesos `sqrt(\|G\|)` do `grpreg`, pedaços finos entre `b_n` e `2b_n − 1` colunas, níveis com `2^j < b_n` num pedaço só; E2.5e), seguido de limiar nas normas por bloco; o LASSO coordenado de D3 fica como opção (D43). A tese numérica passa a ser: vence o spline sintonizado por REML ou validação cruzada no não homogêneo e na `mixed` em `n ≥ 500`, empata com o sintonizado por GCV, perde no suave perto do fator 1,5, e recupera a estrutura que o spline não recupera. D18 fica (extensão de K&P), com o parágrafo de posicionamento a reescrever sem "penalize coefficient by coefficient rather than in blocks"; a teoria em blocos de E1.11 (`08a`, §§ 11 e 12) vai às derivações numeradas e ao manuscrito em `k = 4` | recomendação do chat principal (33(f)), aceita pelo autor: com o limiar oráculo só os blocos passam a perna do não homogêneo; a estrutura deles chega a 1 com o `cv1se` (E2.5j); a taxa perde o logaritmo em `π ≥ 2`; a SS publica seleção em coeficientes variáveis (Wei, Huang & Li 2011), e o JCGS do plano não é necessário |
| D45 | 10-03 | **O limiar (pergunta 38):** a regra padrão de `t` é o **`cv1se`** (o maior `t` a um erro-padrão do mínimo da validação cruzada, nas dobras do ajuste); o `+cv` fica como opção de predição; a porta do QUT e o reajuste `+ls` ficam fora do padrão e no código; **o limiar é também o estimador de predição**, e a cota de risco do estimador limiarizado tem de ser escrita (bloco ativo zerado tem `ν ≤ t + ‖ĝ − g‖`; os mantidos não mudam) | E2.5j: o `cv1se` chega ao teto do limiar oráculo no suave, no `uneven`, no nulo e em `n = 1000`, e é o que segura o fator do suave e a estrutura; o `+cv` é melhor no sinal fraco em `n ≤ 500`; a porta é redundante com o `cv1se` e tem nível acima do nominal; o `+ls` piora os blocos fora do suave; aceita pelo autor |
| D46 | 10-03 | **O spline de E4 e da tabela do manuscrito (pergunta 41):** duas colunas, **`gam.reml`** e **`gam.gcv`**, com o `gam.cv` em nota (razão 1,000 contra o `gam.reml`); a grade de `k` é a de D41, com o custo do topo registrado; o `wafc.gcv` não vai à tabela de E4 e fica no código; o `k.check()` de Wood (2017, §5.9) entra em E4; a tabela de convergência dos motores é relida em E4 se `n` ou `q` crescerem | o `gam.gcv` é o critério padrão do `mgcv::gam()` e o que tira a vantagem no não homogêneo, e omiti-lo seria a objeção; o `gam.reml` dá o erro do `gam.cv` a um décimo do custo; aceita pelo autor |
| D43 | 10-03 | **Nada medido em E2.5 é descartado.** A decisão de rumo escolhe a variante principal, a regra de `t` e as colunas do `gam`; as outras formas (LASSO puro e limiarizado, `+ls` e `+lsb`, as cinco formas do block LASSO, `wafc.gcv`, `wafc.sglasso`, `gam.matched`, `gam.k128`, `gam.gcv`) ficam no código, nos testes e nas tabelas juntadas, como opção, como coluna ou como resposta ao referee; E3.1 não remove variante sem nova decisão | pedido do autor: os achados podem ser úteis (o LASSO limiarizado é o melhor em estrutura, o `gam.gcv` é a referência que mais aperta, o `wafc.gcv` custa um décimo) |
| D32 | 09-20 | **A seleção de estrutura entra pela limiarização (saída (c)), e a saída (a) não abre agora.** O Corolário 8 de `06-selecao-limiar.tex` é o enunciado do artigo, com a hipótese de separação numerada e dizendo no próprio enunciado que é **estimação seguida de limiar**. A sondagem de E1.7a fica registrada como observação de meia página (a condição em grupos não depende de `J` nem da base), e a saída (a), se voltar depois de E2.5, volta pela rota de Wei & Huang (2010), não pela de Bach | veredito de E1.7a: a redução algébrica tira o risco de a condição falhar por construção, mas não o custo, e exige (BD) mais `E(XX'\|U)` constante, que contraria D13; além disso o valor de (a) continua condicionado a E2.5 escolher a variante em grupos, e E1.7c mostrou o LASSO limiarizado acertando 10 de 10 onde ela acerta 0 de 10 |
| D31 | 09-20 | **Em avaliação numérica repetida, a base é fixada e avaliada por tabela**, não pelo algoritmo de Daubechies-Lagarias a cada ajuste: construir a tabela uma vez com `WaveBased::wtable()` para o par `(family, filter.size)` e passá-la em `wavelet.table` de `wafc_design()`, em vez de deixar a regra `use.table = "auto"` decidir réplica a réplica. Vale para o piloto (E2.4), o compêndio (E4) e a aplicação (E6). **Exceção:** conferência que mede precisão fina (as de `derivations/check/`, que leem decaimento até `1e-11`) continua com avaliação exata ou tabela com `prec.wavelet` alto, porque ali o `3.1e-06` engoliria o que se quer medir | pedido do autor, e a razão é **tempo de execução**: a tabela é o caminho rápido e a aproximação é boa o bastante (o erro de interpolação medido em E2.1 é `3.1e-06`), de modo que a variação entre os dois caminhos não é problema prático. A regra `auto` só dispara em `n q ≥ 2000 L`, isto é `n q ≥ 16000` com `L = 8`, e portanto **não dispara** nos `n` do piloto: deixá-la decidir significa pagar Daubechies-Lagarias em toda a varredura |
| D30 | 09-19 | **O cenário suave é para ganhar, não só para não perder.** A hipótese `C^{p+1}` de Xue & Yang delimita a garantia deles, não o desempenho: com `J` e `λ` por validação cruzada o WAFC pode vencer splines também no suave, e E2.4 e E4 têm de medir isso em condição justa — o concorrente sintonizado nos termos dele (nós por BIC como no artigo deles, e `mgcv` com REML), predição e ISE relatadas em separado, e uma componente **suave de curvatura desigual** (gaussiana estreita ou `doppler` truncado longe da singularidade) acrescentada ao `dgp.R`, que é `C^∞` e portanto dentro da hipótese deles, mas com escala variando ao longo do domínio | argumento do autor; se o ganho aparecer, é ilustração forte para o artigo, e se não aparecer, a paridade no suave já é o que a Seção 5 precisa |
| D29 | 09-19 | **Idioma:** inglês no código de `wafc/` (Roxygen e comentários), português nas derivações e nos documentos de trabalho, inglês americano no manuscrito (D6). O `macros.tex` imprime os ambientes em português; o manuscrito não o usa, porque as macros estão copiadas para dentro dele e o template define os próprios ambientes | decisão do autor; `wafc/` vira pacote em E3.3 e documentação de usuário é em inglês, enquanto a derivação é documento de trabalho e revisar prova em português é mais rápido |
| D28 | 09-19 | **Seleção de estrutura: tentar (a) e (c), com (c) garantida.** A seleção por limiarização (saída (c), corolário do Corolário 4 mais separação) entra no artigo de qualquer forma e cobre o LASSO puro; a propriedade oráculo em grupos (saída (a)) é tentada e entra se fechar. A ordem é sondagem antes de prova: E1.7a decide a viabilidade de (a) com álgebra sob `X ⊥ U` e verificação numérica, e E1.7c faz (c) em paralelo | decisão do autor; as duas não são exclusivas, e o valor de (a) é condicionado a E2.5 escolher a variante em grupos, o que os números de E2.2 e E2.3 não indicam |
| D27 | 09-19 | **Os cenários declaram `s' = 3/2` (suave) e `s' = 1/2` (não homogêneo)**, que é o que E2.3 usou; o número entra em `dgp.R` como atributo do cenário, para E2.4 e E4 não o redescobrirem | a cúbica tem quina na extensão periódica e a teoria é em `eps = 0` (D26); `blocks` e `heavisine` saltam dentro do intervalo, onde base nenhuma ajuda |
| D26 | 09-19 | **A teoria fica em `eps = 0`, na base periodizada; a margem é dispositivo de amostra finita.** O Lema 10 de E1.3b permanece nas derivações como justificativa da margem para `J` na faixa admissível, e não como enunciado de taxa. Ficam **registradas duas rotas de princípio**, para o caso de um referee pedir teoria sem periodicidade: (i) enunciar na base do intervalo (CDV), que dispensa periodicidade, extensão, margem e Gram restrita; (ii) provar a compatibilidade sobre as direções estimáveis, que é, em essência, refazer o trabalho que a CDV já faz | a margem não produz ganho assintótico (o cruzamento das duas exigências é exato); e o manuscrito já enuncia a teoria sem margem, então nada muda nele |
| D25 | 09-19 | **Forma das referências de software:** entram no `.bib` com versão e URL, conferidas no registro oficial do pacote quando não há DOI, numa seção própria do arquivo, e **não** entram em `literatura.md`, que é de trabalho | proposta de L4, ratificável sem custo: software não é literatura a posicionar, é dependência a citar, e a versão é o que a exigência de reprodutibilidade da revista pede |
| D19 | 09-19 | **Interface de `cv.wafc()` e `wafc_tune()`** (E2.3): dobras fixas para toda a grade de `J`, expostas em `foldid`; desenho e caminho de `λ` por candidato construídos na amostra inteira, com as dobras reaproveitando as colunas; empate resolvido pelo menor `J`; `df` do BIC e do EBIC igual a não nulos mais os `p` níveis; `wafc_tune(rule)` como entrada única das cinco regras | segue o `cv.wall()` e é o que torna duas regras comparáveis na mesma réplica; ratificada pelo autor em 09-19 |
| D20 | 09-19 | **O padrão de sintonia do WAFC é `cv.min`**, com `lambda.1se` como variante de estrutura e o BIC como alternativa barata; o EBIC não serve para escolher resolução neste desenho | custo de 1.00 a 1.04 sobre o oráculo da grade contra 1.05 a 1.67 das demais; ratificada pelo autor em 09-19 |
| D6 | 09-18 | Documentos de trabalho em português; manuscrito em inglês americano; convenções de git, marcação e continuidade herdadas do `bdm-draft` | pedido do autor ("em linha com o bdm-draft") |
| D7 | 09-18 | Compêndio de simulação e aplicação em repositório próprio, `wafc-studies`, nos moldes do `wall` | o `wall` já resolveu cache, `renv` por commit e proveniência |

**Propostas de L2.** L2a a L2e foram ratificadas e cumpridas em 2026-09-19;
L2f continua em aberto e é a única que muda o escopo do artigo.

| # | Proposta | Onde |
|---|---|---|
| ~~L2a~~ | **cumprida em 09-19** (D18): frase-tese e contribuições reescritas, "Why not block LASSO?" na lista do referee | `alvo-revista.md` §4 |
| ~~L2b~~ | **cumprida em 09-19:** E1.4 (i) virou "recordar K&P (eq. 1.8, Lema 1) e estender ao desenho aditivo"; o entregável central passa a ser o termo cruzado entre moduladoras e a parte (ii) | `plano-projeto.md` E1.4 |
| ~~L2c~~ | **cumprida em 09-19**, com uma correção: o QUT não é um `type.measure`, é outra regra de `λ`, e por isso entra em **E2.4**, não em E2.3 | `plano-projeto.md` E2.4 |
| ~~L2d~~ | **cumprida em 09-19:** concorrentes mínimos de E2.4/E4: `mgcv`, spline adaptativo (Wang, Jiang & Liu 2024), block LASSO de K&P no mesmo desenho, VCBART; cenário não homogêneo com as funções de Donoho-Johnstone | `plano-projeto.md` E2.4, E4 |
| ~~L2e~~ | **cumprida em 09-19** | `proposta-metodo.md` §5 |
| ~~L2f~~ | **decidida em 09-19 (D28):** tentar (a) e (c). As três saídas estão em [`selecao-estrutura.md`](selecao-estrutura.md)** (2026-09-19): (a) propriedade oráculo para a variante em grupos, que pede irrepresentabilidade no desenho de produtos; (b) variante sem teorema; (c) seleção por limiarização como corolário do Corolário 4 mais uma hipótese de separação, que é a recomendação do assistente | `plano-projeto.md` E1.7 |

**Nenhuma decisão adiada.** D5 e D8, que estavam pendentes desde
2026-09-18, foram ratificadas em 2026-09-19 e estão na tabela acima.

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

1. **~~D5 e D8~~ ratificadas em 2026-09-19.** O alvo é a *Statistica
   Sinica* e o método se chama WAFC; **E5a está destravada**.
2. **Aplicação (E6.1), agora com sondagem feita e negativa.** **Depois de
   D44 (2026-10-03):** a tese passa a estrutura recuperada com predição
   competitiva, o que torna a saída (b) mais forte do que era; **E6.1b
   catalogada** para medir as duas saídas no critério novo, e a escolha
   continua do autor. Texto original: duas saídas,
   e a escolha é do autor: **(a)** abrir uma quarta rodada de busca com o
   critério do §3.3 de E6.1a (salto documentado na literatura da área:
   limiar administrativo, limiar regulatório, quebra datada), o que atrasa
   E6; **(b)** aceitar que a aplicação **empate** com o spline, escolher a
   mais interpretável das três e mover o peso do artigo para a teoria de
   adaptação e para a simulação. A (b) é mais rápida e mais honesta, e
   enfraquece o "método + teoria + simulação + aplicação" que a revista
   espera. Se for (b), housing é a mais forte no papel de aplicação neutra
   (maior `R²`, componentes legíveis, um degrau na longitude da área da
   baía), **mas tem licença não declarada no StatLib**, o que colide com a
   exigência de reprodutibilidade; bike e beijing têm CC BY 4.0 e DOI.
3. **~~Sparse group LASSO em `wafc()`~~ respondida por E2.2:** entrou como
   opção, com o número de seleção de estrutura (3 de 3 blocos nulos contra
   0 do LASSO) e o preço em predição. Qual das duas é a principal continua
   sendo de E2.5.
4. **~~Centralização das componentes~~ decidida (D22):** Lebesgue,
   `∫_0^1 g_{ℓm} = 0`. O `notacao.md` ganhou a §7 com a emenda e com o que
   as fontes fazem; o `ms_1.tex` já estava assim.
5. **~~Periodicidade no enunciado~~ resolvida pelo caminho, e fecho aqui:**
   a proposta (c) de E1.3 é o que está em vigor desde D23 e D26 — teoria na
   base periodizada com `eps = 0`, `boundary` como opção do código —, o
   `ms_1.tex` já foi escrito assim, e E1.8 escreveu a alternativa (b) por
   inteiro, para o caso de ser cobrada.
6. **~~`boundary = "interval"` e o tipo de margem~~ decididos (D23, E1.3b,
   E2.1b).** A margem é **fixa**, porque só a fixa compra a taxa cheia, e já
   está implementada. O que resta é estreito e é de E2.4: **qual valor
   fixo**, medindo ISE e `λ_min(G_eps)` em
   `eps ∈ {0, 2^{−(J+1)}, 1.9^{−J}, 0.02, 0.05, 0.10}`, **mais a família
   `eps = a (L−1) 2^{−J}` com `a` em torno de 1**, que é a escala da faixa
   contaminada: em amostra finita a margem certa acompanha `J`, ainda que
   não compre nada assintoticamente (D26). O padrão atual, `0.05`, está no
   código marcado como provisório, e é largo demais para `J` pequeno e
   estreito demais para `J` grande por esse critério. **Medido em E2.4
   (2026-09-20):** `eps = 0` domina, e `0.05` custa 12% a 15% no cenário
   suave; o padrão passou a `0` com a ratificação da P2 (D35, 2026-09-28).
7. **~~Posicionamento diante de Klopp & Pensky~~ decidido (D18):**
   extensão deles. O parágrafo aprovado e a lista de contribuições reescrita
   estão em `alvo-revista.md` §4, e L2a está cumprida.
8. **~~A comparação com o `gam` foi em dimensão casada?~~ Respondida por
   E2.5a (2026-09-30, §2):** em dimensão casada o `gam` vence no suave
   (1,67 a 1,69 em ISE) e na curvatura desigual (1,50 a 1,60), empata no
   não homogêneo com `q = 2` e vence na `mixed`. O fator do critério
   continua por fixar (pergunta 33). Texto original: não, o piloto usou
   o padrão `k = 10` de `wafc_fit_gam()`, e E6.1a mostrou que é isso que
   separa "o WAFC ganha 6 a 15%" de "empata". **E2.5 tem de remedir o
   cenário não homogêneo com `k` casado ao espaço de aproximação**, o que
   E2.4b tornou barato (`wafc_k_matched()` e o motor `bam`, que custa
   `0,589 s` em `n = 1000` contra `1,955 s` do `k = 10` por `gam`).
9. **~~P1 e P2~~ ratificadas em 2026-09-28 (D34, D35):** grade padrão de
   `J` em `2:8` e margem padrão `0`, já no código.
10. **~~Conserto obrigatório antes de E2.5 rodar~~ feito em 2026-09-28:**
   o `wafc_competitor()` repassava `wavelet.table` a todo concorrente, e o
   `VCBART` parava com "argumento não utilizado" dentro de um `try()`
   silencioso. Agora o argumento só chega a `klopp` e `oracle`, e falha de
   método vira linha com a mensagem em vez de sumir. Os números já
   registrados **não** são afetados (foram produzidos antes do ajuste a
   D31, e estão na tabela de razões da análise de E2.4).
11. **Calibração do limiar `t_n`: ~~catalogada como E2.5g~~ medida por
   E2.5g (2026-10-01, §2):** `c = 0,15` serve ao suave e não ao não
   homogêneo; `c = 0,4` acerta 0,88 a 1 em `n ≥ 500`; nenhuma regra
   relativa zera o nulo. A escolha da regra é a pergunta 38. Texto
   original: (nova, de E1.7c) o Corolário 8 é
   explícito em que `t_n` depende de constantes desconhecidas, e a regra
   grosseira `t = 0.15 max_{ℓm} ‖ĝ_{ℓm}‖` acertou 1.00 e 0.90 nos dois
   cenários em `n = 1000`. Calibrar é de E2.4/E2.5, com a curva de acerto
   contra o limiar como instrumento.
12. **O cenário `smooth` de `dgp.R` fica como está?** `⟨sine, cubic⟩ = 0.969`
   (E1.7a): as duas componentes de `β_1` são quase proporcionais, o pior caso
   para medida de estrutura. Trocar o `cubic` por algo ortogonal ao `sine`
   tornaria o cenário menos adversário gratuito; manter deixa a medida
   conservadora. Decide-se junto com a componente suave de curvatura desigual
   que D30 já pediu.
13. **(BD) vira hipótese nomeada?** Se a observação de E1.7a entrar no artigo,
   a independência entre moduladoras que ela exige é mais forte que D11 e
   precisa de nome e de uma frase dizendo o que exclui.
14. **Block LASSO de K&P:** entra em `wafc()` como opção de penalidade, ao
   lado do sparse group LASSO de E2.2, ou fica só como concorrente em E2/E4?
   **Com número desde E2.5a:** é o agrupamento que compra o ganho sobre o
   LASSO, e com os níveis livres ele domina o LASSO em predição nas três
   células medidas em `n = 1000`. A pergunta virou parte da 33.
15. **~~Edições acumuladas para `k = 2`~~ aplicadas por E5c
   (2026-09-30, §2).** Texto original: (decidido: não tocar em `ms_1`
   avulso). **`k = 2` aberta em 2026-09-30; as edições são a tarefa E5c**,
   que acrescenta à lista a §4.2, cujo `ε` "da ordem de `2^{−J}`" contraria
   D35, e o Corolário 8 de D32 na Seção 3, com a prova no `supp` (decisão
   do autor, 2026-09-30); a troca de "sieve" segue a marcação de D38. Entram de uma vez: **a troca de "sieve" por "approximation
   space" (D33), que são 25 ocorrências no `ms_1.tex` e 5 no `supp_1.tex`**;
   a frase de §4.3
   sobre a regra teórica, que vira número com o que E2.3 mediu (`λ_n` de 7 a
   13 vezes o `λ` ótimo; sensível a `s'`, não sistematicamente acima); as
   citações de software que L4 trouxer; a observação de extensão de E1.3b; e
   a frase de reprodutibilidade no resumo ou na discussão.
16. **~~Idioma do código e das derivações~~ decidido (D29):** inglês no
   código, português nas derivações. O `macros.tex` já imprime os ambientes
   em português e os quatro arquivos em português foram recompilados. Falta
   alinhar o `03-desenho-produtos.tex`, que está em inglês e hoje imprime
   cabeçalhos em português; é tradução, catalogada como E1.4c e **feita
   em 2026-09-20**.
17. **~~`rescale = TRUE` como padrão~~ decidido (D23):** continua padrão, e
   agora com razão declarada, não herdada. O deslocamento de centralização
   que ele causa é um nível, e a leitura correta é que a componente é
   identificada no suporte.
18. **~~Quem escolhe `J`~~ não era pergunta:** o plano sempre disse
   `cv.wafc()` sobre `(J, λ)`, como o `cv.wall()` (`plano-projeto.md` E2.3);
   o handoff de E1.6 leu o plano como se ele só falasse de `λ`. O que fica
   de E1.6 é **matéria de medição para E2.3**, não bloqueio: a regra teórica
   `J_n` do Teorema 2 ficou dois níveis acima do `J` que minimiza o erro
   realizado em toda a varredura, ao custo de 5% a 11% de erro, e as
   escolhas de `J_n` e de `c` dependem de `s'` e `τ`, que ninguém conhece.
   E2.3 compara a regra teórica com a validação cruzada.
19. **~~A compatibilidade sobre as direções estimáveis~~ resolvida por
   decisão de rumo (2026-09-19), e o caminho está em D26.** Registro do que
   se aprendeu, porque é fácil reabrir por engano:

   - A saída que eu havia proposto, restringir a Hipótese 5 a
     `2^J < (L−1)/(2 eps)`, **não é assintótica**: o autor notou que ela
     obriga `eps → 0` para `J` crescer. Com `eps = 2^{−cJ}` a condição vira
     `2^{J(1−c)} < (L−1)/2`, que só admite `J → ∞` se `c ≥ 1`, e em `c = 1`
     o expoente garantido pelo Lema 10 é **exatamente `1/2`**, tanto para
     `π ≥ 2` quanto para `π < 2`: a mesma taxa que a periodização dá sem
     margem nenhuma. No ponto em que a margem passa a ser compatível com
     `J → ∞`, o ganho é zero.
   - **Por que o zero não é acidente.** As duas exigências são a mesma
     faixa, na mesma escala. Excluir a região contaminada pela emenda pede
     `eps ≳ (L−1)2^{−J}`, porque é essa a largura das wavelets que
     atravessam `0 ≡ 1`; preservar a Gram pede `2 eps < (L−1)2^{−J}`, para
     que nenhuma wavelet caiba inteira na margem. A margem teria de ser
     simultaneamente mais larga e mais estreita que `(L−1)2^{−J}`. Daí o
     `1.9^{−J}` do `wall()` ter dado ganho **parcial** (0.96 bit medido em
     E1.3b): ele fica logo acima do cruzamento.
   - **A justificativa do autor para ligar `eps` a `J` continua válida, e é
     numérica:** quanto maior `J`, mais estreita a faixa contaminada, e
     menos margem é necessária. É implicitamente o que motiva a construção
     de Cohen, Daubechies e Vial, que corrige as wavelets da borda em vez
     de negociar a largura da faixa.
20. **Cota inferior, registrada como E1.9 no plano** (mais visível depois de
   D18, porque o artigo se declara extensão de quem tem a dele). Não existe
   aqui, e Klopp & Pensky têm a deles para `q = 1` com `X ⊥ U`. **A decisão
   de abrir fica para depois dos resultados de E2.5 e E4**, por escolha do
   autor: se a evidência numérica sustentar o artigo, a cota é cortesia que
   se responde em revisão; se a contribuição numérica ficar fraca, ela vira
   o peso que falta. As três saídas: (a) citar K&P e dizer que a cota
   superior atinge a referência do modelo de sequência, que é o que o
   `05-taxas.tex` faz hoje; (b) abrir E1.9 e construir a cota para `q ≥ 2`,
   trabalho do porte de E1.4 mais E1.5; (c) restringir a afirmação de
   otimalidade a `q = 1`. O referee da SS pode cobrar a (b).
21. **`p` crescente com `n`.** A teoria fixa `p` e `q`. Se a aplicação de E6
   tiver `p` grande, o termo `σ² p / n` deixa de ser de ordem menor e a
   janela do Corolário 5 estreita; mudar isso mexe em E1.3 e E1.4, não só em
   E1.6.
22. **~~Bibliografia, pontos de L3~~ fechada por L5 (2026-09-30, §2):**
   âncora em Mallat (2009, §7.5.1); Tsybakov por D39. Resta só levar a
   correção de `hardle1998wavelets` ao WALL, se o autor quiser (pergunta
   35). Texto original: (nenhum bloqueia) a identidade
   `Σ_l φ(x − l) ≡ 1`, usada na prova do Lema 1(i) de E1.2, ficou **sem
   âncora** — a expressão não ocorre em Daubechies (1992), a quem estava
   atribuída, e os candidatos a conferir são Härdle et al. (1998, cap. 5) e
   o Mallat; a citação de Meyer (1992, cap. III §11) para a base
   periodizada não foi conferida (a CUP bloqueia acesso automatizado) e foi
   substituída no texto pela de Daubechies (1992, §9.3), que foi; e
   `hardle1998wavelets`, copiada do WALL, grafa "Tsybakov, Alexandre" onde
   o Crossref diz "Alexander" (uniformizar nos dois repositórios ou deixar?).
   **Decidido em 2026-09-30 (D39):** cada entrada grafa o nome como está
   impresso na própria referência; a folha de rosto de Härdle et al.
   (1998) diz "Alexander", e a correção é de L5.
   Proposta: abrir uma frente curta só para a âncora da partição da unidade
   quando E5a precisar dela.
23. **~~O teto de páginas~~ adiado por decisão do autor (D21):** escrever
   sem contar, medir no fim e então decidir entre resumir mais e mandar
   coisa ao suplemento.
24. **~~A rota do intervalo~~ escrita e fechada em E1.8** (2026-09-19),
   em `derivations/07-rota-intervalo.tex`. **É ela a resposta** a um referee
   que peça teoria sem periodicidade, e não a rota (ii) de D26, que continua
   sem existir. Levar a rota ao manuscrito é reescrever as provas de E1.3 a
   E1.6 com a base trocada: redação, não matemática nova. O levantamento
   original: a base CDV garante os quatro pilares sem negociação. Aproximação, porque caracteriza `B^s_{π,r}[0,1]`
   sem periodicidade, com a mesma prova e sem `C(eps)`. Identificabilidade,
   porque a reparametrização `[Φ Q | Ψ]` da §5 de `01-identificabilidade.md`
   dá colunas ortonormais em `L_2[0,1]` e de integral zero, e o bloco fica
   com `2^J − 1` colunas, como no periódico. Desenho, porque a prova de
   E1.4 só usa ortonormalidade em `[0,1]` e as cotas da densidade, de modo
   que `λ_min(Σ) ≥ κ_1 c_U` sai verbatim, sem `λ_min(G_eps)`. Oráculo e
   taxas, porque consomem só essas duas peças. **Duas ressalvas:** os
   `p q (2^{j_0} − 1)` coeficientes de escala não penalizados são constantes
   em `n` mas pesados em amostra finita (`15 p q` com filtro 8, `63 p q`
   com filtro 20) e o `wbasis()` exige `j_0 ≳ log_2(5L)`, o que tira os `J`
   pequenos da grade; e a cota pontual `Σ ψ²_{jk} ≤ C_ψ 2^J`, que E1.4 usa,
   foi conferida na base periódica e **não** na CDV.
25. **~~Qual `s'` cada cenário declara~~ confirmado pelo autor (D27):**
   `s' = 3/2` no `smooth` e `s' = 1/2` no não homogêneo, que é o que E2.3 já
   tinha usado. O `3/2` vale porque a teoria é enunciada em `eps = 0`
   (D26), onde a quina da cúbica na extensão periódica é real; com margem e
   extensão o mesmo cenário leria `s' = 4`, e é uma das coisas que a rota do
   intervalo (E1.8) tem de deixar escritas.
26. **~~Cinco referências~~ fechadas por L4 (2026-09-19)**, e as seis
   citadas no manuscrito por E5c (2026-09-30). Texto original: o que sobra é
   pequeno e vai junto com `k = 2`: copiar as seis chaves novas para
   `references_1.bib` e citá-las onde L4 propôs (o lasso na Introduction e
   na seção do estimador; `glmnet` e R na computação; `WaveBased` onde as
   bases são avaliadas; `sparsegl` na variante em grupos; Johnstone ao lado
   de Donoho & Johnstone 1998). Duas perguntas de forma ficaram: o Johnstone
   é `unpublished` (o que se confirma) ou se mantém a forma do WALL, que não
   se sustenta; e o `citation("WaveBased")` pede também a entrada do método,
   mas os métodos que ele lista não são os que o WAFC usa, então ficou só a
   do pacote. Texto original da pergunta: as mesmas cinco, e que não estão
   verificadas, por isso hoje `glmnet` e `WaveBased` aparecem sem citação:
   Tibshirani (1996), Friedman, Hastie & Tibshirani (2010), a citação do R,
   o `sparsegl`, e o Johnstone de modelos de sequência. Frente curta nos
   moldes de L1 e L3, antes de E5b.
27. **Pendências de acabamento do manuscrito**, todas para a semana da
   submissão e nenhuma bloqueando (estão no checklist de `alvo-revista.md`
   §6): (a) **~~o `chicago.bst`~~ dispensado pelo autor (2026-10-01):** o
   `chicago` do template abrevia em "et al." a partir de três autores e não
   põe o volume em negrito, contra a §4 das instruções, mas a revista não
   tem `.bst` próprio e diagrama o artigo aceito a partir da fonte, então
   fica como está; (b) substituir os `\input` pelas cópias dos arquivos de macro, de
   modo que o pacote enviado seja autocontido, registrando a correspondência
   entre bloco copiado e arquivo de origem; (c) autores, afiliações, e-mails
   na última página e agradecimentos, que hoje são os marcadores do
   template.
28. **~~Bibliografia, quatro pontos deixados por L1~~ ratificados pelo autor
   em 2026-10-01 (D42)**, como propostos; o `.bib` já estava assim. Texto
   original: (nenhum bloqueia; o
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
29. **~~Duas acelerações do WAFC~~ decidida em 2026-09-28 (D37): nenhuma
   entra.** Medidas em 2026-09-28 (§2). A em etapas passa no LASSO (255 de 255) mas ganha só
   1,5× e fica mais lenta em `n = 1000`; é reprovada nos grupos (104 de
   129); o `1e-8` nos grupos é reprovado por uma réplica (128 de 129).
   O autor decidiu que a em etapas não entra (D37). Texto original:
   duas acelerações a medir antes de adotar (revisão de 2026-09-28). (a) **Validação cruzada em etapas:** o
   `glmnet` calcula o caminho em sequência, então as dobras num prefixo do
   caminho dão exatamente as mesmas soluções; calcula-se o prefixo e só se
   estende se o mínimo estiver perto do fim. Ganho estimado de 5 a 10× no
   `cv.wafc()` e fim dos avisos de não convergência da cauda. O `λ.1se`
   nunca muda sozinho (é o primeiro ponto do caminho com
   `cvm ≤ cvm[imin] + cvsd[imin]`, logo antes do `imin`); o `λ.min` muda só
   se a curva tiver um segundo mínimo, mais baixo, além da guarda, e o
   `cvm` devolvido fica truncado na cauda. (b) **`thresh = 1e-8` no
   `sparsegl`:** 2× mais rápido, KKT limpo, mas os `cvm` mudam um pouco e a
   escolha pode pular um ponto da grade em caso apertado; não medido.
   Critério proposto para as duas: `J.min`, `λ.min` e `λ.1se` idênticos em
   todas as réplicas das células do piloto, senão fica o comportamento
   atual. **Depois da P1**, porque o ganho de (a) depende da profundidade
   da grade.
30. **~~O `bsgl` acompanha a grade de D34?~~ Sim, D36 (2026-09-28).** Ele escolhia o número de funções
   B-spline por bloco entre 4, 8 e 16, que são as dimensões de `J ≤ 4`;
   com o WAFC indo a `J = 8` (256), a comparação que isola a base deixa de
   ser em dimensão comparável, que é a lição de E6.1a. Proposta: estender
   os candidatos a `2^(2:8)`, como o `gam` casado de (c). Decide-se antes
   de repetir o piloto.
31. **~~Por que a validação cruzada do sparse group LASSO prefere, no
   cenário nulo e em `J` profundo, o fim do caminho?~~ Respondida em
   2026-09-28 (§2):** maldição do vencedor numa faixa plana da curva, sem
   defeito nem vazamento, com o mesmo custo em predição do LASSO; o preço
   é estrutura (centenas de coeficientes falsos em `λ.min`). Texto
   original: (medição de
   2026-09-28, §2). O `cvm` desce até o último `λ` num ajuste denso que é
   pior contra a verdade, e o LASSO no mesmo dado não faz isso. Enquanto a
   causa não for achada, a coluna `wafc.sglasso` do piloto repetido não
   pode ser lida no nulo, e a escolha entre as duas variantes de E2.5 tem
   de levar isso em conta. Candidatos a examinar: a solução do `sparsegl`
   no fim do caminho (a penalidade de grupo ainda encolhe em bloco, o que
   pode dar um ajuste denso de variância baixa), e se a diferença de `cvm`
   entre `J = 2` e `J = 8` é maior que o ruído das dobras.

32. **Pendências de E5c** (2026-09-30), nenhuma bloqueando:
   - (a) **~~Colisão de símbolos na §3.6~~ resolvida (D40, 2026-10-01,
     §2):** `ν_{ℓm}`, `𝒮`, `ν_{min,n}` e `Δ̄_n`, no `06` e na `k = 3`, mais
     a colisão de `D` com `D = p + d`, que a pergunta não listava. Texto
     original: herdada de `06-selecao-limiar.tex`:
     `N_{ℓm}` (norma do bloco) colide com `N` (momentos nulos) e com `N_J`;
     `S` (conjunto de blocos) com o `S` de suporte genérico; `δ_n` com o
     vetor `δ`. Proposta da tarefa: `ν_{ℓm}` ou a norma por extenso, e
     `\mathcal{S}` para a estrutura. Passa pelo `notacao.md`.
   - (b) **~~Aplicada em 2026-10-01~~:** item (iii) do Corollary 2 no `ms_3`, a janela `Δ̄_n ≤ t < ν_{min,n} − Δ̄_n` num evento de probabilidade `≥ 1 − α − α′ − ω_n`, com `α` o nível de `eq:lambda`, e a prova no `supp_3` (S7) identificando `ω_n`; marcado em `colR1`. Texto original: o Corollary 2 enuncia só triagem e recuperação; a janela não
     assintótica `D_n ≤ t < δ − D_n` ficou na prova. Subir ao enunciado?
   - (c) **~~Aplicada em 2026-10-01~~:** as partes (i) e (ii) valem também na resolução do Theorem 1, com `ρ̃_n = (log n/n)^{s/(2s+1)} ≤ ρ_n` (separação mais fraca com `π < 2`); a (iii) fica no Theorem 3. Parágrafo depois do Corollary 2 no `ms_3`, prova no fim da S7 do `supp_3`, Observação `rem:compressivel` no `06`, `ρ̃_n` no `notacao.md` §9; marcado. Texto original: o Corollary 2 cobre só o regime do Theorem 3 (aproximação linear);
     a mesma prova daria a versão com `ρ_n = (log n/n)^{s/(2s+1)}` sob o
     Theorem 1. Enunciado novo: acrescentar?
   - (d) (**Deixada para a revisão** em 2026-10-01, por recomendação do chat principal: o lema não dá taxa em `ε = 0`, a resposta pronta a um referee é a rota do intervalo de E1.8, e ele convidaria a pergunta pela taxa com margem, cuja resposta é a Proposição 6. Sem decisão do autor.) Levar o Lema 10 ao `supp` como Proposition S3.2, com o expoente
     exato e a observação sobre a margem acoplada a `J`?
   - (e) **~~Catalogado como L5~~ verificado por L5 (2026-09-30):** Teorema
     3.2, pp. 688–749. Texto original: catalogado como L5, com as duas referências de
     E1.11 e a pergunta 22. van de Geer, Bühlmann & Zhou (2011), o LASSO limiarizado, é a
     âncora natural do procedimento em dois passos, mas não está verificado
     (faltam as páginas). Frente curta de bibliografia?
   - (f) A frase de reprodutibilidade ganha endereço quando houver (E3.3,
     E4.1).
   - (g) Os números da §4.3 precisam de tabela em E5b.
   - (h) **~~Aplicada em 2026-10-01~~:** a Introduction do `ms_3` passa de três a quatro contribuições, a quarta sendo a recuperação de estrutura pelo limiar (Corollary 2), sem condição de desenho nova nem irrepresentabilidade, e com a triagem sem separação; marcado. Texto original: a lista de contribuições da Introduction não menciona a
     recuperação de estrutura; acrescentar como quarta ou deixar?

33. **~~O rumo depois do no-go de E2.5a~~ decidido em 2026-10-03 (D44),
   com a recomendação (f).** Texto original: (2026-09-30), a mais importante
   da lista. O plano (`plano-projeto.md` E2.5 e tabela de riscos) manda,
   antes de mudar de rumo, testar `boundary = "interval"`, pesos
   adaptativos e limiarização em blocos, e, se persistir, reposicionar o
   artigo para seleção de estrutura e mudar o alvo para o JCGS. E2.5a
   aponta um candidato anterior a esses: **o block LASSO de K&P com os
   níveis livres**, que já está em `wafc_fit_klopp()` e é a penalidade de
   D3 com o agrupamento de K&P. Quatro pontos a decidir, nesta ordem
   ((a) catalogado como **E2.5b** e a sondagem teórica de (b) como
   **E1.11**, ambos em 2026-09-30, em paralelo; E1.11 fechou no mesmo dia):
   - (a) **~~Medir o `klopp` de níveis livres~~ feito por E2.5b
     (2026-09-30, §2):** vence o `wafc.lasso` em toda célula com componente
     (4% a 10% em `rmse_f`), empata no nulo, empata com o `gam.matched` na
     `mixed`, vence-o no não homogêneo a partir de `n = 500` e perde no
     suave (ISE 1,49 a 1,59 no `smooth`, 1,30 a 1,36 no `uneven`). Falta a
     forma com pesos 1, que é a da teoria (pergunta 34, E2.5c).
   - (b) O `klopp` conta como concorrente ou como variante do WAFC? Mesma
     base, mesmo desenho, outra penalidade. Adotá-lo não contraria D18 (o
     artigo já se diz extensão de K&P), mas muda o que "o WAFC" é.
     **E1.11 (2026-09-30) tirou a objeção teórica:** o oráculo de E1.5
     transfere sem cone com a calibração em blocos, e o Corolário 5 ganha
     `(log n)^{2s'/(2s+1)}`, sem logaritmo em `π ≥ 2`; o custo de escrever é
     de ~5 páginas e 2 a 3 dias, mais `k = 3`. A teoria cobre pesos 1, não os
     `sqrt(|G|)` do `grpreg` (pergunta 34). Se não for adotado, a sondagem
     fica como resposta ao "why not block LASSO?".
   - (c) O fator do suave. Medido contra o `gam.matched`: 1,67 a 1,69
     (`smooth`) e 1,50 a 1,60 (`uneven`) para o LASSO; 1,45 a 1,61 e 1,31
     a 1,38 para o `klopp`. Com 1,5, os dois reprovam no `smooth`. **Contra
     o `gam` autônomo com `k = 128` (E2.5d) o fator piora:** 1,58 a 1,92 e
     1,59 a 1,81 para o LASSO; 1,47 a 1,68 e 1,38 a 1,63 para o
     `klopp.free`.
   - (d) **~~Catalogado como E2.5d~~ respondido por E2.5d (2026-09-30,
     §2):** o `gam` autônomo empata com o `gam.matched` a 1% a 2% e o
     veredito se mantém; a escolha do `gam` de E4 é a pergunta 36. Texto
     original: o `gam.matched` herda o
     `J` da busca do WAFC na mesma réplica, e não é método autônomo; E4 precisa de um `gam` que escolha `k` sozinho
     (REML com `k` generoso é a prática usual, não medida aqui).
   - (e) **O limiar, medido por E2.5g (2026-10-01, §2):** é ele que passa
     o fator do suave (1,31 a 1,37 no `wafc.lasso+max` contra o
     `gam.matched`, 1,22 a 1,45 contra o `gam.k128`); com limiar o LASSO
     volta à frente dos blocos no suave em `n ≤ 500`, e os blocos ficam à
     frente no não homogêneo, no `uneven` e na `mixed`, onde o
     `klopp.balanced+cv` vence o `gam.matched` e o `wafc.lasso+cv` perde
     dele. A escolha entre LASSO e blocos é agora de célula, não de
     método; é o que a decisão de rumo tem de pesar.
   - (f) **Recomendação do chat principal (2026-10-03, números no §2,
     "o teto do limiar e o `gam.gcv`"), sem decisão do autor:** go
     reposicionado, na *Statistica Sinica*, com o **block LASSO balanceado
     como o WAFC** e o LASSO como opção (D43). A razão: com o `t` oráculo
     só os blocos passam a perna (b) (0,90 a 0,95 no não homogêneo, 0,93 a
     0,96 na `mixed` contra o `gam.cv`), a estrutura deles chega a 0,94 a 1
     com uma regra boa, e a taxa perde o logaritmo em `π ≥ 2` (E1.11). A
     tese numérica passa a ser: vence o spline por REML ou validação
     cruzada no não homogêneo e na `mixed`, empata com o spline por GCV,
     perde no suave perto do fator 1,5, e recupera a estrutura que o spline
     não recupera. O JCGS do plano não é necessário: a SS publicou Wei,
     Huang & Li (2011), seleção em coeficientes variáveis. Custo: E1.11 no
     `ms` e no `supp` (~5 páginas, 2 a 3 dias, `k = 4`) e a frase "penalize
     coefficient by coefficient rather than in blocks" do parágrafo de D18.
     Antes da decisão: **E2.5j** (a regra de `t`, pergunta 38). O teste de
     `boundary = "interval"` que o plano manda fazer antes de mudar de rumo
     é trabalho de código (a reparametrização do bloco de escala,
     `01-identificabilidade.md` §5), não opção a ligar; é a única alavanca
     que sobra na perna (a), porque a cúbica tem quina na extensão periódica
     (D27). **Depois de E2.5j (§2):** com o `klopp.balanced+cv1se`, a perna
     (a) fica na linha (ISE 1,37, 1,50, 1,52 contra o `gam.cv`) e a (b)
     passa só em `n ≥ 500` (0,930 e 0,895 no não homogêneo; 0,987 e 0,932
     na `mixed`), com `P(Ŝ = S) = 1` no suave, no `uneven`, no nulo e em
     `n = 1000`; com o `+cv`, a (b) passa nos três `n` e a (a) reprova. A
     recomendação não muda.

34. **~~Os pesos do block LASSO~~ medidos por E2.5c e E2.5e (§2):** a
   forma balanceada (E2.5e) é coberta pela teoria com os pesos do `grpreg`
   (`ρ² ≤ 1,67` no piloto) e é a melhor fora do `smooth`; reconcilia em
   parte. Medição de E2.5c: os
   pesos 1 da teoria são a pior forma na predição; o padrão do `grpreg`
   vence no suave e a junção dos níveis grossos no não homogêneo e na
   `mixed`. A quarta forma que talvez reconcilie teoria e prática é
   conjectura (§2), **catalogada como E2.5e** (2026-09-30), com a parte
   teórica de saber se E1.11 aceita pesos de razão limitada. Texto
   original: (de E1.11, 2026-09-30) o `klopp` e o
   `klopp.free` do código usam os pesos `sqrt(|G|)` do `grpreg`; a teoria
   de E1.11 cobre pesos 1, os de K&P, e com `sqrt(|G|)` a cota volta à
   ordem do LASSO. **E2.5b roda com `sqrt(|G|)`.** Antes de a pergunta 33
   decidir pela variante em blocos, medir pesos 1, e talvez uma terceira
   forma que junta num pedaço só os níveis com `2^j < b_n` de cada bloco
   (tira os pedaços unitários sem mudar os pesos dos cheios), nas células de
   E2.5b. É uma linha em `wafc_kp_groups()` ou `wafc_fit_klopp()`, arquivos
   de E2.5b; abre depois dela. **Catalogada como E2.5c** (2026-09-30).

35. **Pendências de L5** (2026-09-30), nenhuma bloqueando:
   - (a) **~~Fechada por L6~~ (2026-10-01, §2).** Texto original: levar as correções bibliográficas aos `derivations/`: Mallat
     (2009, §7.5.1) na prova do Lema 1(i) do `01-identificabilidade.md`;
     Hsu et al. como Teorema 2.1, com a nota da `A` quadrada, e Cai com
     páginas e Teorema 4 no `08a-sondagem-blocos.md`; páginas e Teorema 3.2
     de van de Geer et al., e a ressalva do reajuste, no
     `06-selecao-limiar.tex` (recompilar o PDF). Tarefa curta, sem conflito
     com as de medição.
   - (b) Na próxima rodada do manuscrito, Mallat (2009) na prova do Lemma 1
     do `supp`; Hsu et al. e Cai só se a variante em blocos for adotada.
   - (c) (**Conferido em 2026-10-02 na folha de rosto da impressão da Springer**, `refs/hardle1998-springer-fm.pdf`, p. iii: "Gerard Kerkyacharian" sem acento e "Alexander Tsybakov", como no `.bib` daqui, que cumpre D39. **O WALL divergia da impressão**; corrigido em 2026-10-02 a pedido do autor na última versão, `references_theo_8.bib` do `wall-manuscript` ("Gérard" → "Gerard"; o "Alexander" já estava certo ali), sem commit lá.) O `hardle1998wavelets` do WALL continua com "Gérard" e
     "Alexandre"; levar D39 ao `references_theo_1.bib` de lá é decisão do
     autor. "Gerard" sem acento veio da versão do seminário, da página da
     Springer e do Crossref; a folha de rosto da impressão não foi vista.
   - (d) Hall, Kerkyacharian & Picard (1999, *Statist. Sinica* 9) e Zhou
     (2010, arXiv:1002.1583) não verificados; entram só se o manuscrito
     precisar.
   - (e) **~~Marcas de verificação que L6 achou~~ resolvidas por L7
     (2026-10-01, §2).** Texto original: Donoho & Johnstone (1994) com `j_0` fixo, no `08a`
     (`[VERIFICAR: não relido]`); no `06-selecao-limiar.tex`, a seção do
     LASSO limiarizado em Bühlmann & van de Geer (2011) e o número do
     resultado de Huang, Horowitz & Wei (2010); e, no `08a`, a numeração do
     *Annals* de Klopp & Pensky e de Lounici et al., lidos no arXiv.
   - (f) **~~Feita por L9~~ (2026-10-01, §2):** a entrada de 1994 saiu do
     `05` e a linha de K&P do `busca-novidade.md` está na forma do *Annals*;
     o que sobrou é a 35(h). Texto original: a marca que L7 achou e não resolve com número: o
     `05-taxas.tex` (l. 835) atribui a Donoho & Johnstone (1994) "a forma
     original" do risco minimax sobre bolas weak-`ℓ_τ`, e o artigo não trata
     disso; o mais próximo é a desigualdade oráculo (Teorema 1, p. 437;
     Corolário 1, p. 440). Ou a frase vira "a desigualdade oráculo que a
     precede", ou a entrada sai. E a frase de `busca-novidade.md` l. 45 com
     as condições do arXiv de K&P, que no *Annals* são `L + 1 = n^ς` e
     `r* > (2ς)^{−1}`. As duas são edição curta; a primeira recompila o PDF.
   - (g) **Suplemento de K&P lido em 2026-10-01** (`refs/klopp2015-supp.pdf`):
     os Lemas 2 a 5 estão lá com os números que o `08a` usa, e o Lema 2(ii)
     é a eq. (5.25), o termo de viés pela norma dual, como o `08a` diz. Dois
     pontos para o `08a` §1.1 (**aplicados em 2026-10-01**, com uma
     correção ao que eu tinha escrito aqui): o suplemento rotula a condição
     do Lema 4 como "(A4)", mas a que ele escreve, eq. (5.30), é a (A3) do
     artigo, e a (A4) de lá é a do ruído, de modo que o "(A3)" do `08a`
     estava certo e o rótulo do suplemento é que destoa; e o lema pede
     `1 ≤ ν < ∞` e `r > min(1/2, 1/ν)`, de modo que a correspondência
     `(r, ν) = (s, π)` só vale com `π ≥ 1`. **O fator 2 do Corolário 1,
     resolvido em 2026-10-02:** é erro do enunciado publicado. A razão entre
     a cota superior (3.19) e a inferior (3.6) do *Annals* é
     `(log n)^{(2−ν)_+/(ν(2r+1))}` vezes `(log p/log n)^{2r/(2r+1)}`, que a
     própria (3.20) registra e que fica limitada sob `n^β ≥ p`; o "2" da
     (3.21) não sai de nenhuma das duas. O arXiv v2 (`refs/klopp2014-arxiv.pdf`)
     não confirma nem desmente, porque tinha outro corolário, com
     `(log p)^{2r/(2r+1)}` no lugar de potências de `log n`: o resultado
     mudou na revisão, e o 2 entrou com ela. Não há errata no Crossref. O
     `08a` §4.3 lia sem o 2, o que está certo, e agora diz por quê; quem
     citar a ótima a menos de log de K&P cita (3.19)–(3.20), não a (3.21). Texto original: dois pontos de K&P que só o arXiv ou o suplemento resolvem: o
     Corolário 1 do *Annals* (3.21) tem um fator 2 no expoente do log que a
     (3.17) e a (3.20) não têm, e o `08a` §4.3 o lê sem o 2 (parece erro de
     impressão, não conferido no arXiv); e o enunciado dos Lemas 2 a 4, que
     estão no suplemento (`10.1214/15-AOS1309SUPP`). Opcional; se o autor
     baixar o suplemento, entra como `refs/klopp2015-supp.pdf`.
   - (h) **~~A classe da referência minimax~~ decidida pelo autor, saída (a), e aplicada (2026-10-01):** no `05`, a Observação `rem:tres` diz "corpos de Besov (que o Lema `lem:weak`(i) põe dentro de bolas weak-`ℓ_τ`)" e cita os Teoremas 4 e 5, e a entrada de 1998 perdeu a marca; no `ms_3` (§3.3, l. 576), "weak `ℓ_τ` balls" virou "Besov bodies" e a citação ganhou "Theorems 4 and 5", marcado em `colR1`; Johnstone (2019) ficou ao lado, sem conferir o resultado. Texto original: (de L9). A ordem `λ_n^{2−τ}`
     que o `05` (Observação `rem:tres`) e o `ms_3` (l. 577) citam está nos
     Teoremas 4 e 5 de Donoho & Johnstone (1998), mas sobre **corpos de
     Besov**, não sobre bolas weak-`ℓ_τ`. Duas saídas: (a) trocar "weak
     `ℓ_τ` balls" por "Besov bodies" nos dois lugares, com "Theorems 4 and
     5", o que casa com a frase seguinte (o Besov já implica a
     compressibilidade, Lema 9) e não pede leitura nova; ou (b) manter
     weak-`ℓ_τ` e apoiar só em Johnstone (2019) ou em Donoho & Johnstone
     (1994b), com o resultado conferido no texto. A redação de (a) para a
     entrada e para a Observação do `05` está no handoff, reproduzida aqui:
     a entrada passa a "Teoremas 4 e 5 (§4.1, pp. 890–891): no modelo de
     sequência gaussiano com ruído `ε`, o risco minimax sobre corpos de
     Besov `Θ^α_{p,q}(C)` é `≍ C^{2(1−r)} ε^{2r}`, com `r = 2α/(2α+1)`; com
     `α = s` e `τ = (s+1/2)^{−1}`, isto é `C^τ ε^{2−τ}`. O artigo não trata
     de bolas weak-`ℓ_τ`."; e a Observação, "referência minimax do modelo
     de sequência gaussiana sobre corpos de Besov (que o Lema `lem:weak`(i)
     põe dentro de bolas weak-`ℓ_τ`) com ruído por coordenada `λ_n`". No
     manuscrito, (a) é edição marcada em `colR1` na `k = 3`.

36. **Pendências de E2.5d** (2026-09-30), nenhuma bloqueando:
   - (a) **Qual `gam` é o concorrente de E4 e da tabela do manuscrito** (medido por E2.5h, 2026-10-02: a decisão está na pergunta 41(a)).
     (**Em curso como E2.5h**, que absorveu a E2.5i em 2026-10-01: o
     `gam` com `k` escolhido por REML, por GCV e por validação cruzada nas
     dobras do WAFC, numa grade única `{5, 10, 20, 40, 80}` e não em
     potências de 2, com o `gam.gcv` fora da `mixed` (D41); ~15,5 h em 8
     núcleos. Wood
     (2017, §5.9) recomenda `k` generoso conferido pelo `k.check()`, não
     busca, o que apoia o `gam.k128`; a busca é a versão simétrica ao `J`
     do WAFC. As referências vão ao `.bib` por L8, catalogada no mesmo
     dia.)
     A tarefa propõe o `gam` autônomo com `k = 128` por moduladora: não
     aperta, não diverge, é o que um usuário rodaria. O `gam.matched` fica
     como a medição de E6.1a (base contra dimensão), não como concorrente.
     O custo na `mixed` é de ~8 min por ajuste.
   - (b) A `mixed` em 50 réplicas para o `gam` autônomo: pede um controle
     como o `WAFC_REPS_MIXED` no `09-gam-autonomo.R` e ~6,5 h de
     processador.
   - (c) Guarda no `wafc_fit_gam()` (livre desde que E2.5c fechou): recusar, ou marcar, ajuste com `edf`
     acima do número de coeficientes, que é a assinatura das duas
     divergências de E2.5a. Toca `competitors.R`, arquivo de E2.5c.
   - (d) **~~Feito na integração de E2.5c.~~** A linha do `09-gam-autonomo.R` no `wafc/README.md` entra na
     integração de E2.5c, que tem o arquivo: "`scripts/09-gam-autonomo.R`
     ✓ | E2.5d | `Rscript wafc/scripts/09-gam-autonomo.R [n_rep] [partes]
     [ncores] [ns] [células] [ks]` (padrão `50`, `fit,report`, `64,128`):
     o `gam` autônomo ao lado do `gam.matched` refeito, no sorteio da parte
     `competitors` do `04-pilot.R` (funções copiadas, o piloto intacto; a
     `mixed` só nas 15 réplicas de E2.5a); `report` confere o `gam.matched`
     contra E2.5a antes de qualquer tabela; lê `WAFC_E25A` e `WAFC_E25B`,
     grava `<WAFC_TAG>-fits.rds`, `-edf.rds` e `-joined.rds` em
     `WAFC_OUT`".
   - (e) E2.5c pode ler o fator do suave contra o `gam.k128` de
     `wafc/cache/e25d/e25d-fits.rds`, que cobre as mesmas réplicas (fora
     as 35 da `mixed` além da 15ª).

37. **~~Uma quinta forma do block LASSO~~ medida por E2.5f (2026-10-01,
   §2): não paga**, e sai da lista de candidatas; a forma balanceada
   continua sendo a que a teoria cobre e a que melhor prediz fora do
   `smooth`. Texto original: (**catalogada como E2.5f** em
   2026-10-01; conjectura de E2.5e, não medida
   nem escrita): deixar livres os níveis com `2^j < b_n` de cada bloco, que
   a forma balanceada penaliza juntos. Na teoria eles iriam para o bloco
   não penalizado `A`, como os coeficientes de escala da Proposição 5 de
   E1.8, com `σ²p_0/n` e `p_0 ≍ pq b_n`, de ordem `log n/n`; os pedaços
   restantes têm `ρ² < 2`. Pode recuperar o `smooth`, onde a perda é a
   penalização conjunta desses níveis, e custar no nulo. Medir: uma opção
   em `wafc_kp_groups()` e ~40 min em 8 núcleos.

38. **~~Pendências de E2.5g~~ decididas em 2026-10-03 (D45).** Resta
   escrever a cota de risco do estimador limiarizado (§5). Texto original:
   (2026-10-01), que entravam na decisão de rumo.
   **(a) e (b) medidas por E2.5j (2026-10-03, §2):** o `cvrel` repete o
   `+cv` e sai da lista; o `cv1se` chega ao teto no suave, no `uneven`, no
   nulo e em `n = 1000`, e perde do `+cv` em `n ≤ 500` no não homogêneo e
   na `mixed` por zerar bloco ativo fraco; a porta do QUT resolve o nulo do
   `+max` e do `+cv`, é redundante com o `cv1se` e tem nível acima do
   nominal (12% a 17%) porque o `J` é escolhido antes. **Recomendação do
   chat principal:** o `cv1se` como padrão (é o que serve à tese de
   estrutura e à perna (a)), o `+cv` como opção de predição, e a porta fora
   do padrão, no código (D43). Não medidos: um meio-termo (meio
   erro-padrão) e a porta seguida do `+cv`. Recomendação do chat principal para (c): sim, o limiar é também o
   estimador de predição, e falta escrever a cota de risco (bloco ativo
   zerado tem `ν ≤ t + ‖ĝ − g‖`, e os mantidos não mudam); para (d): o
   `+ls` piora os blocos fora do suave (1,12 a 1,22 no não homogêneo contra
   o `gam.cv`) e não vai ao padrão deles, mas fica (D43).
   - (a) **Qual regra de `t` vai ao artigo e a E4.** O `max` com
     `c = 0,15` é o melhor no suave e falha no não homogêneo e na `mixed`;
     um `c` de 0,3 a 0,4 teria de ser fixado sem olhar a verdade, ou por
     uma regra mista (o `c` relativo escolhido por validação cruzada), que
     não foi medida. A tarefa escolheu `t` absoluto na regra `cv`; a escala
     relativa talvez corrija os falsos positivos do suave (uma opção em
     `wafc_threshold()` e ~1 h de processador sem o `+lsb`).
   - (b) **O nulo:** nenhuma regra relativa o zera. Acrescentar uma porta
     de "tudo zero" antes do limiar, por exemplo o QUT de
     `wafc_lambda_qut()`, que já mata o nulo?
   - (c) **O limiar entra só como seleção de estrutura (D32) ou também como
     o estimador de predição da variante?** No suave os números pesam para
     a segunda leitura; no manuscrito, o Corollary 2 da §3.6 já é
     "estimação seguida de limiar".
   - (d) Se o `+ls` entrar (é a melhor regra de estrutura que não vê a
     verdade), entra só no suporte, não no bloco inteiro.

39. **Pendências de L8** (2026-10-01), nenhuma bloqueando:
   - (a) **~~PDFs que fecham três `[VERIFICAR]`~~ lidos em 2026-10-01**
     (`refs/zou2007.pdf`, `tibshirani2012.pdf`, `marra2011.pdf`), e o
     `.bib` não tem mais marca. Zou et al. (2007, Teorema 1, p. 2177)
     **exigem `rank(X) = p`**, o que falha no WAFC com `d > n`; quem cobre o
     desenho do WAFC é Tibshirani & Taylor (2012): Teorema 2 (p. 1213),
     `df = E[rank(X_A)]` para qualquer `X`, e a Observação da p. 1215, que
     com colunas não penalizadas dá `df = p + E[rank(M X_A)]`, com `M` a
     projeção que tira os níveis. **A contagem "não nulos + p" do GCV da
     E2.5h é esse `df` quando `M X_A` tem posto coluna cheio**, o que o
     artigo não garante em geral; o manuscrito deve citar Tibshirani &
     Taylor e dizer a condição, e a E2.5h pode conferir `rank(M X_A)`
     contra `|A|` nos `λ` escolhidos. Marra & Wood (2011), §2.1, p. 2374:
     a dupla penalidade é o `select = TRUE`, dito no próprio artigo. Texto
     original: PDFs que fecham três `[VERIFICAR]` no `.bib`: Zou, Hastie &
     Tibshirani (2007) e Tibshirani & Taylor (2012), o teorema que dá os
     graus de liberdade do LASSO e se o caso com colunas não penalizadas
     (os `p` níveis do WAFC) está enunciado. **É a justificativa do `df` do
     GCV da E2.5h**, e tem de estar resolvido antes de o GCV ir ao
     manuscrito; e Marra & Wood (2011), a seção da dupla penalidade. Como
     `refs/zou2007.pdf`, `refs/tibshirani2012.pdf`, `refs/marra2011.pdf`.
     Opcionais, só para D39: `wood2003`, `wood2011`, `wood2015`,
     `wood2017jasa`.
   - (b) (**Wood 2004 entrou no `.bib` em 2026-10-02**, porque o `gam.gcv` da E2.5h é o `bam` com `method = "GCV.Cp"`; Wood, Pya & Säfken 2016 continua fora.) **As duas do `citation("mgcv")` que ficaram fora**, já conferidas
     no Crossref para entrarem sem nova rodada se o handoff da E2.5h disser
     que o código as usa: Wood (2004, *JASA* 99(467), 673–686), se o
     `gam.gcv` for pelo `bam` com `GCV.Cp`; e Wood, Pya & Säfken (2016,
     *JASA* 111(516), 1548–1563 sem a discussão), que hoje não se usa.
   - (c) (Depois da E2.5h: Wood 2008 não se aplica, porque o `gam.gcv` usa o `bam` e não o `gam()`; Li & Wood 2020 continua candidata, junto do `discrete = TRUE` do `gam.reml`.) **Três candidatas que a ajuda do `mgcv` cita**, conferidas só para
     registro: Wood (2008, *JRSS-B* 70(3), 495–518), se o `gam.gcv` cair no
     `gam()`; Li & Wood (2020, *Stat. Comput.* 30(1), 19–25), junto do
     `discrete = TRUE`; e Wood (2025, *Annu. Rev. Stat. Appl.*), não
     conferida.
   - (d) O `chicago.bst` não imprime o `eprint` de um `misc`: na cópia
     para o `references_{k}.bib` do manuscrito, as entradas do arXiv
     (Pya & Wood e as duas de Schnaidt Grez & Vidakovic) precisam de
     `howpublished`. É da próxima rodada do `.tex`, que copia também as 13
     entradas novas.

40. **~~O Lema 9(i) é um mergulho clássico~~ (a) e (b) aplicadas em 2026-10-01:** no `ms_3`, a contribuição 2 da Introduction e uma frase depois do Lemma `lem:weak` citam a Proposição 9.15 de Johnstone e dizem o que o lema acrescenta (as constantes uniformes em `J`); a citação de Johnstone na l. 576 ganhou "Theorem 9.3"; tudo marcado em `colR1`. No `05`, a nota no início da prova do Lema 9 e a entrada de Johnstone, que perdeu o `[verificar]`; no `alvo-revista.md`, a contribuição 2 reescrita. O §3.3 remete ao lema e não repete a citação. Texto original: (2026-10-01, lendo
   Johnstone 2019, `refs/johnstone2019.pdf`). A Proposição 9.15 dele
   (§9.7, p. 273) é `b^α_{p,q} ⊂ wℓ_{p_α}`, `p_α = 2/(2α+1)`, sob
   `p > p_α`, que é exatamente o Lema 9(i) de E1.6 (`τ = (s+1/2)^{−1}`, sob
   `π(s+1/2) > 1`). O `alvo-revista.md` §4 lista como contribuição 2 "a
   leitura de que a compressibilidade não custa hipótese", e o `ms_3`
   (§3.3 e Lemma `lem:weak`) a apresenta sem essa citação. O que o WAFC
   acrescenta é a versão com constantes explícitas e uniformes em `J` para
   o oráculo do desenho de produtos, e o uso dela para dispensar a hipótese
   de compressibilidade do WALL; o mergulho em si não é novo. Decidir: (a)
   citar Johnstone (2019, Proposição 9.15) no Lemma `lem:weak` e no §3.3,
   e reescrever a contribuição 2 como uso do mergulho; (b) também citar o
   Teorema 9.3 dele (§9.4, p. 260), o risco minimax sobre bolas weak-`ℓ_p`,
   ao lado de Donoho & Johnstone (1998) na l. 576, o que daria apoio
   conferido à forma weak-`ℓ_τ` que a 35(h) trocou por Besov. As duas são
   edição marcada na `k = 3`, e (a) toca o `alvo-revista.md`.

41. **~~Pendências de E2.5h~~ decididas em 2026-10-03 (D46).** Texto
   original: (2026-10-02), que entravam na decisão de rumo.
   Recomendação do chat principal (2026-10-03): (a) duas colunas, `gam.reml`
   e `gam.gcv`, com o `gam.cv` em nota (razão 1,000); omitir o `gam.gcv`
   seria a objeção, porque é ele que tira a perna (b) e é o padrão do
   `gam()`; (b) fica como D41; (c) o `wafc.gcv` não vai à tabela de E4 e
   fica no código (D43); (d) **respondida por E2.5j (§2):** os cortes dos
   dois motores existem e não tocam o `λ` escolhido; o `grpreg` chegou a
   9 015 das 10 000 iterações na `mixed` em `n = 250`, e E4, com `n` ou `q`
   maiores, relê a tabela de convergência (`wafc_cv_convergence()`,
   `extra$conv`) antes de confiar no padrão; (e) em E4.
   - (a) **Qual `gam` vai a E4 e à tabela do manuscrito.** A tarefa propõe
     o `gam.cv` (mesmo critério e mesmas dobras do WAFC, e o mesmo erro do
     `gam.reml` e do `gam.k128`); o `gam.reml` dá o mesmo erro a um décimo
     do custo. Nos dois casos, com a ressalva do topo da grade no não
     homogêneo e na `mixed` em `n = 1000` (1,0% a 1,4%). O `gam.gcv` é a
     referência que mais aperta o WAFC no não homogêneo e a que mais o
     favorece no nulo; se entrar, entra como terceira coluna.
   - (b) **O topo da grade:** estender a 120 só onde aperta seria a grade
     diferente por célula que D41 recusou; fica registrado o custo de 1,0% a
     1,4%.
   - (c) **O `wafc.gcv` entra em E4?** Não ganha em erro e perde no nulo; o
     argumento é o custo (um décimo) e o par com o `gam.gcv`. Se entrar, o
     manuscrito cita Tibshirani & Taylor (2012) com a condição de posto,
     conferida em 150 de 150 sintonias.
   - (d) As 10 de 150 sintonias em que o `glmnet` não convergiu num `λ`:
     medir em que `J` e se o corte cai dentro da guarda, antes de E4.
   - (e) O `k.check()` de Wood (2017, §5.9) não foi feito.

42. **~~Pendências de E1.12~~ decididas em 2026-10-03 (D47).** Resta o
   Corolário 14 (a taxa lenta em blocos), §5. Texto original: (2026-10-03),
   que pedem o autor antes da `k = 4`. Recomendação do chat principal em cada item:
   - (a) **D16: o Corolário 11 passa a ser o enunciado principal**, no
     lugar do Corolário 5. Recomendação: sim; é a consequência de D44, e é
     o enunciado que vence o de K&P em `s'` e não tem logaritmo em `π ≥ 2`.
   - (b) **Notação:** a razão dos pesos `ρ` colide com a taxa `ρ_n` (D40), e
     `W(·)` encosta em `W_J`. A tarefa propõe `ϱ` e `𝒲(·)`. Recomendação:
     aceitar os dois, e com eles a lista do handoff entra no `notacao.md`
     §10 e no `macros.tex` (o `08` troca as macros locais).
   - (c) **A Proposition 3 do manuscrito (a taxa lenta, a Proposição 4 de
     E1.6):** o `08` não tem versão em blocos à parte; o Teorema 3(i) a dá
     com `‖θ*‖_{𝒢,w} ≤ w_max‖θ*‖_1`. ~~Recomendação: cortar do corpo e
     deixar uma frase com a consequência no supp, o que poupa página.~~
     **Recomendação revista (2026-10-03):** manter no corpo, em versão de
     blocos, como corolário do Teorema 3(i) com prova de uma linha. É o
     único enunciado do artigo sem condição de desenho (D16 e o próprio
     `ms_3`, "the unconditional statement of this paper"), e cortá-lo
     custaria isso para poupar meia página. A ordem deve ser a mesma da
     Proposition 3 atual, com a constante `ϱ`: com pesos 1,
     `λ^𝒢_{n,1} ≍ sqrt(log n/n)`, a ordem do `λ` do LASSO, e
     `‖θ*‖_{𝒢,1} ≤ ‖θ*‖_1`. Isso é leitura, a conferir ao escrever.
   - (d) **A Observação 2 (variante branca) vai ao supp** como o enunciado
     que cobre o estimador exato do código. Recomendação: sim.

43. **~~Pendências de E3.1~~ decididas em 2026-10-03 (D48).** Texto
   original: (2026-10-03). Recomendação do chat principal:
   - (a) **Ratificar a interface** do `wafc/README.md` (como D17 e D19),
     com os padrões `penalty = "block"` e `threshold = "cv1se"`, o
     `thresholded = NULL` e o `cvm` do bloco na convenção do `cv.grpreg`.
     Recomendação: ratificar; é o estimador medido, reproduzido à última
     casa.
   - (b) **`thresh = 1e-4` no bloco:** manter; `1e-8` não muda a escolha
     em 40 de 40 réplicas e custa 4,1 vezes.
   - (c) **O `grpreg` passa a `wafc_depends`** no `load.R`, com a linha no
     `CONTINUAR.md`: sim, porque o padrão depende dele.
   - (d) O defeito do `wafc_thr_object()`: corrigir em E3.2.

44. **~~Pendências de E1.13 e E3.2~~ decididas em 2026-10-03 (D49).**
   Texto original: (2026-10-03). Recomendação do chat principal:
   - (a) **Emendar a Proposição 4 do `05-taxas.tex`** com as duas condições
     implícitas (`s' > 1/4` em (i), `s' > s/2` em (ii)), uma linha no
     enunciado e uma na prova: sim, porque o LASSO continua como opção
     (D43) e o enunciado hoje afirma mais do que a prova dá.
   - (b) **"No design condition" no manuscrito:** não escrever a variante
     sem cota superior de `Σ̂`; dizer que nenhuma cota inferior é usada e
     que a única propriedade de `Σ̂` é `λ_max(Σ̂) ≤ Λ`, que vale sem `c_U` nem
     `κ_1`.
   - (c) **O Corolário 14(ii) em `π ≥ 2` sobe ao §3** numa frase depois do
     Theorem `thm:main`: a taxa do teorema principal sem condição de desenho
     em `s < 1/2`. Recomendação: sim.
   - (d) **Dois símbolos para o `notacao.md` §10:** `‖θ‖_{𝒢,1}`, a norma de
     grupos com pesos 1, e `ϑ_π = min(1 − 1/π, 1/2)`, o expoente da
     contagem (não há `\vartheta` no `ms_3` nem no `supp_3`). Recomendação:
     aceitar.
   - (e) **Ratificar as duas assinaturas de `plot`** do `wafc/README.md`:
     sim, são as do código, conferidas.
   - (f) **As figuras ficam versionadas** em `wafc/man-figures/` (290 KB):
     sim, o README as mostra.
   - (g) As cinco funções sem roxygen fora de "Internals" e um argumento
     que escolha os blocos desenhados, para `p × q` grande: ficam para E3.3
     e E6.2.

45. **Completar a Proposição 4 de E1.6 e o Corolário 14 de E1.12 nos
   regimes que D49 excluiu** (2026-10-03, pedido do autor: os resultados
   servem no futuro). Os regimes `s' ≤ 1/4` em (i) e `s' ≤ s/2` em (ii)
   (só `π < 2`) ficaram fora porque a escolha de `J_n` passa de `n`, e a
   concentração das normas das colunas (Lema 6 de E1.5, e o análogo em
   blocos) pede `2^{J}log d/n → 0`. **Leitura do chat principal, a
   conferir:** nesses regimes o viés domina, e basta parar `J_n` no maior
   nível admissível, `2^{J_n} ≍ n/(ω_n log n)` com `ω_n → ∞` tão devagar
   quanto se queira; a cota fica `O_p((ω_n log n/n)^{2s'})`, mais lenta que
   a dos regimes de (i) e (ii), mas consistente para todo `s' > 0`. Com isso
   a proposição e o corolário passam a cobrir todo `s' > 0`, com três
   regimes em vez de dois. **Catalogada como E1.14 em 2026-10-03**, com a
   fronteira `s = 1/2`, que nenhum dos dois regimes atuais cobre.

---

## 5. Próximos passos

**Onde parou (2026-10-03).** **E2 fechada.** O autor aceitou as
recomendações das perguntas 33, 38 e 41 (D44 a D46): o WAFC é o block
LASSO balanceado seguido do limiar `cv1se`, o LASSO fica como opção, o
alvo continua a *Statistica Sinica*, e o spline de E4 entra com `gam.reml`
e `gam.gcv`. **E1.12 fechada e integrada.** **E3.1** (a interface) em
**E3.1 e E3.2 fechadas e integradas**, com a interface ratificada (D48) e
as assinaturas de `plot` ratificadas (D49); E3.2 (gráficos e documentação) e E3.3 (empacotamento, do
autor) seguem. **E6.1b** (a aplicação) está pronta e espera a **ordem do autor** para a
rodada de 4 processos: o script 10 cobre seis bases (as três de E6.1a,
`beijing.heat`, `marylebone` e `kelmarsh`, as duas últimas baixadas com a
permissão do autor), a fumaça das seis passou sem falha, e o lançamento
está em `wafc/cache/e61b/run.sh`. A estimativa é de **~24 h de relógio em 4
processos** (~96 h de processador, pico de ~18 GB dos 31 GB da máquina),
com erro possível de um fator 2. **O autor manteve as 20 partições
(2026-10-03)**, pelo erro-padrão; a rodada continua à espera da ordem. O
handoff parcial fica em `docs/handoff-E6.1b.md`. O manuscrito está em
`k = 3`, ainda escrito para o LASSO; a bibliografia sem marca aberta (83
entradas).

**O que D44 a D46 abrem, na ordem sugerida:**

- (a) **Teoria em blocos nas derivações numeradas.** Promover a sondagem
  `08a-sondagem-blocos.md` (§§ 1 a 12) a resultado numerado: o Teorema 1 em
  blocos com a calibração de Hsu, Kakade & Zhang, o risco ideal por pedaços,
  a taxa sem logaritmo em `π ≥ 2` e a cobertura dos pesos `sqrt(|G|)` na
  forma balanceada (`ρ² ≤ 3`). O Corolário 8 passa ao ajuste em blocos, e
  entra a **cota de risco do estimador limiarizado** (D45). D16 diz que o
  enunciado principal é o Corolário 5; com D44 ele passa a ser a versão em
  blocos, o que o autor confirma quando o enunciado existir. Os símbolos de
  E1.11 (`𝒢`, `G`, `b_n`, `|𝒢|`, `w_G`, `‖θ‖_{𝒢,w}`, `𝒢_0`, `W(𝒢_0)`, `Ψ̃_G`,
  `R_𝒢(θ; η)`, `λ_n^𝒢`) passam antes pelo `notacao.md`, com aval do autor.
  **Fechada por E1.12 em 2026-10-03 (§2)**, e a pergunta 42 decidida
  (D47): o Corolário 11 é o enunciado principal, `ϱ` e `𝒲` aplicados, a
  lista no `notacao.md` §10. **Falta o Corolário 14** no `08`, a taxa lenta
  em blocos como corolário do Teorema 3(i) (D47(c)), com a conferência da
  ordem. **Fechado por E1.13 em 2026-10-03 (§2)**, com pendências
  decididas em D49.
- (b) **O parágrafo de posicionamento de D18** (`alvo-revista.md` §4),
  reescrito para a penalidade em blocos, com a frase-tese e as
  contribuições. **Aprovado em 2026-10-03 (D50)** e gravado no
  `alvo-revista.md` §4, com as respostas ao referee revistas; entra no
  `ms` na `k = 4`.
- (c) **E3, a consolidação do código:** o block LASSO balanceado sai de
  `wafc_fit_klopp()` (concorrente) para dentro de `wafc()` e `cv.wafc()`
  como penalidade padrão, com o limiar `cv1se` na interface; o LASSO como
  `penalty = "lasso"`; nada removido (D43). E3.3 (empacotamento) é decisão
  do autor. **E3.1 catalogada em 2026-10-03.**
- (d) **Catalogado como E5d em 2026-10-03**, com D51. **O manuscrito em `k = 4`** (D44): a teoria em blocos no `ms` e no
  `supp` (~5 páginas), o posicionamento de (b), a §3.6 com o limiar em
  blocos e a §4 com a interface de (c). Abre depois de (a) e (b).
- (e) **E4 e E6 em paralelo**, depois de E3.1: o compêndio `wafc-studies`
  com os métodos de D46; e a aplicação, que depende da **pergunta 2**,
  ainda aberta e a maior lacuna para o formato da revista. **E6.1b
  catalogada em 2026-10-03** para dar à pergunta 2 os números do critério
  novo.

**Medições opcionais** (não catalogadas; uma rodada por vez): `n = 2000`
no não homogêneo e na `mixed`, que diz se a vantagem que cresce com `n`
continua e se o `gam.gcv` fica para trás; `boundary = "interval"`, que pede
antes a reparametrização do bloco de escala (`01-identificabilidade.md`
§5) e é a única alavanca que sobra no fator do suave; um meio-termo entre o
`cv1se` e o `+cv` (meio erro-padrão); o `gam.k128` na `mixed` com 50
réplicas (36(b)); a guarda de `edf` no `gam` (36(c)).

**Sem bloquear:** 12 (cenário `smooth`), 13 (nome da (BD)), 20 (cota
inferior, o risco assumido com D18; com os blocos, a taxa atinge a cota de
K&P em `q = 1` e `X ⊥ U`), 21 (`p` crescente); 32(d), o Lema 10 no `supp`.

**Depois:** E5b (Seções 5 a 7, com a tabela da §4.3, 32(g), e o endereço
de reprodutibilidade, 32(f)); E7 (submissão, com o teto de D21 e os itens
27(b) e (c) do checklist).

## 6. Histórico de sessões

| Data | O que aconteceu |
|---|---|
| 2026-10-03 | E5d (o manuscrito em `k = 4`) catalogada com D51 (a teoria do LASSO numa seção do supp; remoções inteiras em bloco cinza sem tachado, no `instrucoes.md` §3); mensagem para rodar a E6.1b com 4 processos e 20 partições passada ao autor |
| 2026-10-03 | E1.14 catalogada: a taxa lenta (Proposição 4 e Corolário 14) em todo `s' > 0` e na fronteira `s = 1/2` |
| 2026-10-03 | Posicionamento aprovado (D50) e gravado no `alvo-revista.md` §4, com frase-tese, contribuições e respostas ao referee revistas; pergunta 45 (completar a Proposição 4 e o Corolário 14 nos regimes excluídos por D49) |
| 2026-10-03 | Pergunta 44 decidida (D49): a Proposição 4 de E1.6 emendada com as duas condições implícitas; o Corolário 14 sobe ao §3 em `k = 4`; dois símbolos no `notacao.md`; as assinaturas de `plot` ratificadas; a rodada da E6.1b fica com 20 partições |
| 2026-10-03 | E1.13 fechada e integrada (`OK` em 209 s, 18 páginas): Corolário 14, a taxa lenta em blocos, que pela norma de grupos ganha um logaritmo sobre a Proposição 4 e em `π ≥ 2` dá a taxa do Teorema 4 sem condição de desenho; achadas duas condições implícitas na Proposição 4 de E1.6. E3.2 fechada e integrada (1 202 testes em 48 s, 1 427 com os lentos): `plot.wafc`, `plot.cv.wafc`, roxygen auditado, o `coef` limiarizado corrigido, exemplo no README. Pergunta 44 |
| 2026-10-03 | E1.13 (a taxa lenta em blocos, Corolário 14) e E3.2 (gráficos, documentação e o defeito de D48) catalogadas, em paralelo com a E6.1b; o parágrafo de posicionamento de D18 começa no chat principal |
| 2026-10-03 | Pergunta 43 decidida (D48): interface de E3.1 ratificada, `thresh = 1e-4` no bloco, `grpreg` em `wafc_depends` e no `CONTINUAR.md` |
| 2026-10-03 | E3.1 fechada e integrada (1 075 testes em 46 s, 1 300 com os lentos, oitava junção exata): `cv.wafc(x, u, y)` é o block LASSO balanceado com limiar `cv1se`, e o `wafc.block` reproduz o `klopp.balanced` à última casa; a tolerância `1e-4` do `grpreg` não muda a escolha; pergunta 43 |
| 2026-10-03 | Pergunta 42 decidida (D47): o Corolário 11 é o enunciado principal (emenda D16); `ϱ` e `𝒲` no `08` e no `06` (14 e 11 páginas, sem aviso); `notacao.md` §10; a taxa lenta fica no corpo em versão de blocos, a escrever |
| 2026-10-03 | E1.12 fechada e integrada (`OK` em 207 s, 14 e 11 páginas): Teorema 3 (oráculo em blocos), Lema 15, Teorema 4 e Corolário 11 sem logaritmo em `π ≥ 2`, a variante branca cobrindo o `grpreg`, e o risco do limiarizado (Lema 16, Corolários 12 e 13); o `s' = 1/2` de D27 fica dentro do nosso enunciado e fora do de K&P; pergunta 42. E6.1b parou antes da rodada, à espera do autor |
| 2026-10-03 | E3.1 (a interface com o block LASSO balanceado e o limiar `cv1se` como padrão) e E6.1b (a aplicação no critério de D44) catalogadas, em paralelo com a E1.12, sem arquivo em comum |
| 2026-10-03 | E1.12 catalogada: a sondagem de E1.11 promovida a `08-blocos.tex` (calibração, oráculo com pesos de razão limitada, risco ideal por pedaços, taxas) e, no `06`, o Corolário 8 em blocos e a cota de risco do estimador limiarizado (D45) |
| 2026-10-03 | O autor aceitou as recomendações das perguntas 33, 38 e 41 (D44 a D46): block LASSO balanceado com limiar `cv1se` como o WAFC, LASSO como opção, *Statistica Sinica* mantida, `gam.reml` e `gam.gcv` em E4; E2.5 e E2 fechadas; §5 reescrita em torno da teoria em blocos, do posicionamento, de E3 e da `k = 4` |
| 2026-10-03 | E2.5j fechada e integrada (1 133 testes, 37 350 linhas, sétima junção exata, 1 h 01 min em 10 processos): o `cv1se` chega ao teto no suave, no `uneven`, no nulo e em `n = 1000`; o `cvrel` repete o `+cv`; a porta do QUT tem nível acima do nominal; os cortes dos motores não tocam o `λ` escolhido |
| 2026-10-03 | Leitura dos `.rds` da E2.5h no chat principal: com o `t` oráculo só os blocos passam a perna do não homogêneo; o `gam.gcv` (padrão do `gam()`) tira essa perna; a estrutura é a vantagem robusta; `boundary = "interval"` não está implementado. Recomendação na 33(f); D43 (nada de E2.5 é descartado); E2.5j catalogada |
| 2026-10-02 | Documentação de continuidade revista para recomeçar em outro chat: §5 do `ESTADO.md` reescrita em torno da decisão de rumo; `CONTINUAR.md` §3 e §4, `plano-projeto.md` E2.5 e `TAREFA.md` alinhados |
| 2026-10-02 | E2.5h fechada e integrada (1 006 testes, 28 350 linhas, sexta junção exata): com o `gam` sintonizado por REML, validação cruzada ou GCV o veredito se mantém; `gam.reml` e `gam.cv` coincidem; o GCV do WAFC não ganha nada; o posto de Tibshirani & Taylor confere; Wood (2004) no `.bib`; pergunta 41 |
| 2026-10-02 | Folha de rosto da impressão de Härdle et al. (1998) conferida: "Gerard" e "Alexander", como no `.bib`; o WALL diverge (35(c)) |
| 2026-10-02 | 35(g) fechada: o fator 2 do Corolário 1 de K&P é erro do enunciado do *Annals*, que (3.19), (3.20) e (3.6) desmentem; o arXiv v2 tinha outro corolário; o `08a` já lia certo e agora diz por quê |
| 2026-10-02 | 32(c) aplicada: o Corollary 2 (i)–(ii) também na resolução do Theorem 1, com `ρ̃_n`; no `ms_3`, no `supp_3` e no `06` |
| 2026-10-01 | 32(b) e 32(h) aplicadas no `ms_3` e no `supp_3` (janela não assintótica no enunciado do Corollary 2; a estrutura como quarta contribuição); 34 e 30 páginas |
| 2026-10-01 | D42 (as quatro convenções de L1) ratificada; o `.bst` próprio dispensado (27(a)); o `08a` corrigido com o suplemento de K&P, e o "(A4)" do suplemento é a (A3) do artigo |
| 2026-10-01 | Pergunta 40 aplicada: a Proposição 9.15 e o Teorema 9.3 de Johnstone (2019) citados no `ms_3` (marcado) e no `05`; contribuição 2 reescrita |
| 2026-10-01 | Johnstone (2019) em `refs/`: a Proposição 9.15 dele é o Lema 9(i), o que tira a novidade do mergulho e pede reescrever a contribuição 2 (pergunta 40); o Teorema 9.3 dá o minimax sobre bolas weak-`ℓ_p` |
| 2026-10-01 | PDFs da 39(a) e da 35(g) lidos: o `df` do GCV do WAFC se apoia em Tibshirani & Taylor (2012, Teorema 2 e a observação das colunas não penalizadas), não em Zou et al., que exigem posto cheio; o Lema 4 de K&P pede `ν ≥ 1` (o "(A4)" do suplemento é a (A3) do artigo) |
| 2026-10-01 | Pergunta 35(h) decidida, saída (a): "Besov bodies" com os Teoremas 4 e 5 de Donoho & Johnstone (1998) no `05` e no `ms_3`, marcado; os dois compilam em 10 e 33 páginas |
| 2026-10-01 | E2.5h parou antes da rodada e mediu uma réplica: o `gam.gcv` na `mixed` custa ~900 h na grade até 120; D41 fixa a grade única `{5, 10, 20, 40, 80}`, tira o `gam.gcv` da `mixed` e junta a E2.5i de volta, ~15,5 h no total |
| 2026-10-01 | L9 fechada e integrada: a entrada falsa de Donoho & Johnstone (1994) saiu do `05`; o risco minimax de 1998 é sobre corpos de Besov, com a mesma ordem, e a classe citada no `05` e no `ms_3` é a pergunta 35(h) |
| 2026-10-01 | Pergunta 32(a) resolvida (D40): quatro símbolos da seleção trocados no `06` e no manuscrito, por padrões no modo matemático (31, 101 e 87 trocas); `k = 3` aberta; L9 catalogada |
| 2026-10-01 | L8 fechada e integrada: `.bib` com 82 entradas; Kauermann & Opsomer usam ML, não REML, e o catálogo da E2.5h foi corrigido; o Ruppert aditivo usa um `K` comum até 40; pergunta 39 |
| 2026-10-01 | E2.5h e E2.5i catalogadas (o `gam` por REML e GCV, e por validação cruzada, na grade de Ruppert; o GCV do WAFC); E2.5i com 50 réplicas em todas as células, por decisão do autor; L8 catalogada |
| 2026-10-01 | E2.5g fechada e integrada (921 testes, 24 750 linhas, quinta junção exata): o limiar do Corolário 8 passa o fator do suave pela primeira vez (1,31 a 1,37 no `wafc.lasso+max`), o que sobra é ISE ativo; `c = 0,15` não generaliza ao não homogêneo, `c = 0,4` sim; pergunta 38 |
| 2026-10-01 | L7 fechada e integrada: a numeração de K&P e Lounici et al. conferida no *Annals*; o Teorema 2 de K&P mudou de condições entre o arXiv e o *Annals*, e os Lemas 2 a 4 foram para o suplemento; perguntas 35(f) e 35(g) |
| 2026-10-01 | L6 catalogada e fechada ao lado da E2.5g: as correções de L5 nas três derivações, sem mudança de matemática; o "§2.3" de Hsu et al. no `08a` não existia; quatro marcas fora de L5 ficam na pergunta 35(e) |
| 2026-10-01 | E2.5f fechada e integrada: os níveis grossos livres melhoram os blocos ativos e pioram os inativos, e perdem em toda célula; medido aqui que, com os blocos inativos zerados, o fator do suave do `klopp.balanced` cairia abaixo de 1,5, o que aponta para o limiar do Corolário 8 |
| 2026-10-01 | E2.5e fechada e integrada: a forma balanceada do block LASSO é coberta pela teoria com os pesos do `grpreg` (`ρ² ≤ 1,67`) e é a melhor das quatro fora do `smooth`; nenhuma forma passa o fator 1,5 no suave; pergunta 37 |
| 2026-09-30 | E2.5c fechada e integrada: os pesos 1 da teoria predizem pior que o padrão do `grpreg` (2% a 8%), a junção dos níveis grossos vence no não homogêneo e na `mixed`; a escolha entre formas é de segunda ordem; catálogo vazio, e a pergunta 33 tem todas as medições que pediu |
| 2026-09-30 | E2.5d fechada e integrada: o `gam` autônomo (REML, `k = 64` e `128`) empata com o `gam.matched`, o veredito de E2.5a se mantém e o fator do suave piora (1,58 a 1,92 para o LASSO); o `klopp.free` é o único que o vence em alguma célula; duas divergências do `bam` em E2.5a; pergunta 36 |
| 2026-09-30 | L5 fechada e integrada: `.bib` com 69 entradas; a partição da unidade ancorada em Mallat (2009); Cai (1999) confere com E1.11; Restrepo & Leaf não serve para a base; D39 aplicada; pergunta 35 |
| 2026-09-30 | E2.5b fechada e integrada: o `klopp.free` vence o `wafc.lasso` em toda célula com componente e empata no nulo, reprodução exata de E2.5a; a `mixed` em 50 réplicas confirma o `gam.matched`; a máquina tem 8 núcleos físicos; E2.5c e E2.5d ganham regras de concorrência |
| 2026-09-30 | E1.11 fechada e integrada: o oráculo de E1.5 transfere para o block LASSO com níveis livres sem cone, e o Corolário 5 ganha `(log n)^{2s'/(2s+1)}`, sem logaritmo em `π ≥ 2`; a teoria cobre pesos 1, e o `grpreg` usa `sqrt(\|G\|)` (pergunta 34) |
| 2026-09-30 | E2.5a fechada e integrada (6 450 linhas, 0 falhas, 5 h de relógio): no-go para a variante LASSO pelo critério literal; o `klopp` vence no não homogêneo porque D34 lhe deu a grade larga, e o ganho é do agrupamento; o `gam` casado vence no suave por 1,5 a 1,7 em ISE; D30 não se confirma; pergunta 33 |
| 2026-09-30 | Documentos alinhados ao estado depois de E5c e E2.4c: frase-tese e respostas ao referee do `alvo-revista.md` passam à limiarização (D32) e à §4.3 do `ms_2`; E2.5 e E5c no `plano-projeto.md`; `notacao.md` sem "sieve" e com o `\E` em romano já feito; E1.9 no lugar de E1.8 na pergunta 20 |
| 2026-09-30 | E2.4c fechada e integrada: `uneven` e `gam.matched` no `04-pilot.R`, `WAFC_TAG`, sementes independentes da restrição de células; o `gam.matched` em `J = 8` custa de 72 s a mais de 15 min e até 4 GB, e a `mixed` roda à parte em E2.5a |
| 2026-09-30 | Este chat vira o orquestrador; documentos alinhados a D5, D8 e D34 a D37; E2.4c e E2.5a catalogadas; `k = 2` aberta; D38; E5c fechada e integrada: Corolário 8 no artigo (§3.6, S7), §4.2 alinhada a D35, "sieve" trocado, 33 e 30 páginas |
| 2026-09-28 | Pergunta 31 investigada (`08-sgl-null.R`, reproduzido): o laço de dobras bate com o `cv.sparsegl`, não há vazamento (termo cruzado `t = 0,6`), e o sparse group LASSO no nulo cai numa faixa plana da curva por maldição do vencedor, com o mesmo custo em predição do LASSO e centenas de coeficientes falsos; D37: nenhuma aceleração entra |
| 2026-09-28 | Pergunta 29 medida (`07-accel.R`, 384 réplicas): o prefixo do caminho é exato; em etapas passa no LASSO (255/255) com ganho de só 1,5×, reprova nos grupos (104/129); `1e-8` nos grupos reprova por uma réplica (128/129); e a validação cruzada do sparse group LASSO se engana no nulo com a grade funda (pergunta 31) |
| 2026-09-28 | P1 e P2 ratificadas (D34, D35): grade padrão `2:8` e margem padrão `0`, já no código; a grade nova custa de 10 a 18 vezes no `cv.wafc` e até 13 no `klopp`, o que torna a pergunta 29 urgente; o `bsgl` passou a acompanhar a grade (D36), escolhendo 64 funções em `n = 1000` onde o teto era 16 |
| 2026-09-28 | Revisão do código de `wafc/`: cinco defeitos corrigidos (o `vcbart` volta à tabela; `uneven` não quebra mais o piloto; falha vira linha), `call` compacto e só o melhor ajuste no `cv.wafc` (objeto salvo de 11,4 para 3,6 MB), busca de nós 7 a 11× mais rápida com os mesmos nós; 672 testes, piloto idêntico ao anterior; duas acelerações que podem mudar a escolha ficam para medir (pergunta 29) |
| 2026-09-28 | Documentos de continuidade alinhados ao estado de 09-21: cabeçalho e §5 do `ESTADO.md`, §3 do `CONTINUAR.md`, catálogo e numeração do `TAREFA.md`, marcas de etapa do `plano-projeto.md`, teto de 40 páginas no `instrucoes.md`; perguntas da §4 renumeradas a partir da 16 (havia duas 15) |
| 2026-09-18 | Avaliação de viabilidade; criação do repositório e dos documentos de trabalho; template da EJS; plano E0 a E7 |
| 2026-09-18 | D4 decidida pelo autor (código em `wafc/`, não no `WaveBased`); D5 e D8 adiadas; `prototype/` virou `wafc/`; plano E2 e E3 reescritos; repositório publicado; o autor confirmou o `WaveBased` como dependência e que as funções ficam privadas |
| 2026-09-19 | E2.3 e E5a fechadas e integradas: `cv.min` é o padrão de sintonia (D20) e a regra da teoria custa de 7 a 13 vezes no `λ`; o manuscrito nasce em `k = 1` com 28 e 26 páginas compilando limpo, e o teto vira a decisão urgente |
| 2026-09-21 | Tabela de tempo medida: o ajuste do WAFC custa milissegundos e o custo é a busca; a comparação em dimensão casada só é viável pelo `bam`; o `sglasso` custa 12 a 19 vezes o LASSO |
| 2026-09-21 | E2.4b fechada: `uneven` nasce e lê `s' = 4`; a tolerância apertada é preço da verificação de KKT e não desperdício, contra o que eu tinha concluído; o WAFC deixa de ser o caro da tabela de tempo |
| 2026-09-21 | E1.10 fechada (terminologia) e o documento das inconsistências do WALL criado; medidos os dois defeitos de desempenho do ajuste |
| 2026-09-20 | E6.1a fechada com veredito negativo: nenhuma das três bases sustenta o argumento, e o ganho aparente sobre o spline era de dimensão, não de base — o que põe em dúvida a comparação com o `gam` no piloto |
| 2026-09-20 | E2.4 fechada: a grade de `J` do `cv.wall()` trunca o sieve no ótimo, e alargá-la corta 17,6% do erro e 32% do ISE no cenário não homogêneo, invertendo o veredito contra o VCBART; a margem `0.05` é dominada por `eps = 0` |
| 2026-09-20 | E1.7c fecha a seleção de estrutura com teorema (Lema 13, Corolário 8) e o LASSO limiarizado vira 10 de 10 contra 0 de 10 da variante em grupos; E1.7a dá veredito de escopo reduzido para a saída (a); E1.4c traduz o `03` |
| 2026-09-19 | E1.8 fechada e integrada: a rota do intervalo está escrita e é a resposta a quem pedir teoria sem periodicidade; o cruzamento da margem virou teorema, com as faixas separadas por um nível exato |
| 2026-09-19 | E1.3b e E2.1b fechadas e integradas: a extensão fecha e o preço da margem é `eps^{−(s−1/π)}`, exato; só margem fixa compra a taxa, e ela degenera a Gram restrita a partir de `2^J ≥ (L−1)/(2eps)`, o que deixa a taxa sob a Hipótese 5 conjectural |
| 2026-09-19 | L4 fechada e integrada: `.bib` com 63 entradas, e a entrada do Johnstone corrigida contra a que o WALL carrega |
| 2026-09-19 | E0.3 fechada de vez: o autor confirmou no SCImago que a *Statistica Sinica* é Q1 em Statistics and Probability |
| 2026-09-19 | D18 decidida (o artigo é extensão de Klopp & Pensky) com o parágrafo de posicionamento escrito; D16, D5 e D8 ratificadas; E5a catalogada e destravada |
| 2026-09-19 | E1.6, E2.2 e L3 fechadas e integradas: **E1 inteira**, com a compressibilidade saindo de graça da hipótese de Besov (D16); `wafc()` com as duas penalidades e KKT fechando no caminho inteiro (D17); `.bib` com 57 entradas e duas citações de teorema corrigidas |
| 2026-09-19 | E1.5 e E2.1 fechadas em dois chats de tarefa e integradas: desigualdade oráculo sem cone, com a razão `‖f̂−f‖²_n/(λ²s_0)` estável a 5% em `n` de 200 a 6400; `wafc/` nasce com 103 testes passando; D15 |
| 2026-09-18 | E1.2, E1.3 e E1.4 fechadas em três chats de tarefa e integradas; as três conferências rodam e imprimem `OK` aqui; D13 e D14; a numeração global dos resultados fixada |
| 2026-09-18 | L2 fechada em chat de tarefa e integrada: novidade central confirmada (interseção zero na SS e na EJS), mas Klopp & Pensky (2015) cobre E1.4 (i), E1.5 e E1.6 para `q = 1` e `X ⊥ U`; seis propostas de mudança de rumo a ratificar |
| 2026-09-18 | L1 fechada em chat de tarefa e integrada: 35 referências verificadas, quatro correções de atribuição, dois trabalhos novos para L2 olhar |
| 2026-09-18 | E1.1 fechada: notação congelada (`ψ_{jk}`, `X_ℓ`, `U_m`, `θ_{ℓm,jk}`, `U ∈ [0,1]^q`), D9 a D12, `macros.tex` reescrito e compilando |
| 2026-09-18 | Máquina nova conferida (R 4.6.1, `WaveBased` 2.6-0, `grpreg`, `gglasso`, `sparsegl` por `apt`); E0.3 fechada: instruções da SS transcritas, templates versionados e compilando, teto corrigido de 30 para 40 páginas |
