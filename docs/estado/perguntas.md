# ESTADO antigo: §4 Perguntas em aberto
Movido de docs/ESTADO.md em 2026-11-16, sem alteração de linhas.
Arquivo morto: não ler inteiro; consultar por grep pelo identificador (D56, E6.1b, pergunta 33) e ler o trecho.

<!-- conteudo movido abaixo -->
## 4. Perguntas em aberto

Ordenadas pelo que bloqueia mais.

1. **~~D5 e D8~~ ratificadas em 2026-09-19.** O alvo é a *Statistica
   Sinica* e o método se chama WAFC; **E5a está destravada**.
2. **Aplicação: decidida em D57 (2026-10-05)**, marylebone como principal e
   beijing.heat como segunda, condicionada à E6.1c. **Aplicação (E6.1), agora com sondagem feita e negativa.** **Depois de
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
20. **Cota inferior: fechada no nível da taxa por E1.9 (2026-10-05, §2);**
   fica aberto só o logaritmo em `π < 2`, a constante e `p` crescente.
   Catalogada como E1.9 em 2026-10-05, com escopo
   reduzido (no nível da taxa, pela subclasse com uma componente; a
   constante em `pq` opcional; o logaritmo em `π < 2` fora). A leitura do
   chat principal: com `p` e `q` fixos, a cota `n^{−2s/(2s+1)}` sai da
   regressão não paramétrica univariada com desenho aleatório, e o
   Corolário 11 seria ótimo em `π ≥ 2` para todo `q`; a frase de D50 sobre
   a cota inferior mudaria. Texto original: **Cota inferior, registrada como E1.9 no plano** (mais visível depois de
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

45. **~~Completar a Proposição 4 de E1.6 e o Corolário 14 de E1.12 nos
   regimes que D49 excluiu~~ feito por E1.14 (2026-10-04, §2).** Texto
   original: Completar a Proposição 4 de E1.6 e o Corolário 14 de E1.12 nos
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

46. **~~Pendências de E1.14~~ decididas em 2026-10-04 (D52).** Texto original: (2026-10-04). Recomendação do chat principal:
   - (a) **O comparador truncado na taxa lenta do LASSO.** O Lema 8(i) com
     `θ*·1{|θ*| > λ}` troca `λ‖θ*‖_1` por `Σ_a min(λ|θ*_a|, θ*_a²)`, e em
     `π < 2` a taxa lenta passaria a `(log n/n)^{2s/(2s+1)}` com
     `s' > s/(2s+1)`, a taxa do Corolário 5(iii) sem condição de desenho; nos
     blocos pede uma contagem nova de `Σ_G min(λ‖θ*_G‖, ‖θ*_G‖²)` sob Besov.
     Só a aritmética foi conferida (cinco pares). Recomendação: catalogar,
     porque o autor quer os resultados completos e é o que faria a taxa lenta
     igualar a rápida também em `π < 2`; não bloqueia nada.
   - (b) Uma cota de `σ̂_max` e `λ_max(Σ̂)` pela contagem de pontos nos
     suportes, muito acima do teto: tiraria parte de um logaritmo em regimes
     que o artigo não usa. Recomendação: não abrir.
   - (c) Os itens (iii) e (iv) no supp, com uma frase no corpo: sim (já está
     na tabela da `k = 4`).
   - (d) **Símbolos locais novos** para o `notacao.md` §10: `μ_J`, `λ_n^+`,
     `c_ψ` (no `05`), `μ^𝒢_J`, `λ_n^{𝒢,+}`, `c_μ`, `Δ = Λ − κ_2C_U` (no `08`).
     Recomendação: aceitar, trocando `Δ` por outra letra, porque `Δ` e `Δ̄_n`
     já são o erro por bloco e a sua cota (D40).

47. **~~Pendências de E5d~~ decididas em 2026-10-04 (D52)**; o resumo (b) aprovado em D53. Texto original: (2026-10-04). Recomendação do chat principal:
   - (a) **Os itens de E1.14 na `k = 4`:** a Proposition 3 e a Proposition
     S8.1 entraram sem o regime (iii) e sem `s = 1/2`; levá-los (a frase no
     corpo, os itens e a prova em S6 e S8, `μ^𝒢_J` e `λ_n^{𝒢,+}`, e a frase de
     D49(b) revista). Recomendação: catalogar uma tarefa curta na `k = 4`.
   - (b) **O resumo** foi trocado pela tarefa, fora do catálogo, para não
     contradizer o corpo (block lasso com pedaços, pesos de razão limitada,
     sem logaritmo em `π ≥ 2`, uma frase sobre o limiar, "group lasso" nas
     palavras-chave). O texto não foi aprovado. Recomendação: o chat
     principal o traz ao autor no formato do §3 do `instrucoes.md`.
   - (c) **"programme"** (grafia britânica, desde a `k = 1`) no texto de
     D50: trocar por "program". Recomendação: sim.
   - (d) **A citação do `grpreg`:** o §4.1 diz "[VERIFICAR: Breheny and
     Huang (2015) and the grpreg manual]", e o §2.3 nomeia o group lasso sem
     citar. Recomendação: catalogar uma frente curta de bibliografia (L10)
     com o `citation("grpreg")`, Breheny & Huang (2015, *Statist. Comput.*) e
     Yuan & Lin (2006, *JRSS-B*), verificados como em L1.
   - (e) **O teto de páginas** (D21): 43 páginas hoje sem o removido, ~54
     projetadas com a E5b, contra 40. A §3 tem 19 páginas, e os candidatos
     são o Theorem 3 e o Corollary 1 ao supp, ou a Proposition 3 e o
     Corollary 3 só com o enunciado no corpo. Recomendação: manter D21
     (decidir no fim), mas registrar que a decisão vai ter de mexer na §3, e
     não só nas tabelas.
   - (f) **Ratificar a leitura de D51 nos deslocamentos do supp:** o material
     do lasso que foi de S5 a S7 para a S8 não foi repetido em cinza no
     lugar antigo (uma linha cinza diz o que saiu e para onde), e a S8
     inteira está em `colR1`; repetir dobraria ~14 páginas. Recomendação:
     ratificar.
   - (g) Símbolos do `ms` fora do `notacao.md` §10 (o comparador `β̄`, `b̄`,
     `v̄`, `𝒢̄`, `𝒲̄`, `M_1`, `Σ̂_{GG}`): levá-los ao §10 junto com os da 46(d).

48. **~~Pendências de E1.15, E5e e L10~~ decididas em 2026-10-05 (D54);
   E5f e L11 catalogadas.** Texto original: (2026-10-05). Recomendação do chat
   principal:
   - (a) **O Corolário 15 e a Proposição 8 na `k = 4`**, na forma (a) da
     E1.15: uma ou duas frases no corpo depois da Proposition 3 (a taxa do
     Theorem 1 sem condição de desenho em `s < 1/2`, e o limite em
     `s ≥ 1/2`), o enunciado e a prova no supp (uma Proposition S6.x a partir
     do Corolário 15 e a versão do lasso na S8); e a frase do §3 que atribui
     a taxa do Theorem 1 à condição de desenho passa a dizer "for the
     components, and for prediction when `s ≥ 1/2`". O rascunho da frase
     está no fim desta pergunta. Recomendação: sim, pelo teto de páginas.
   - (b) **O posicionamento (D50) não muda:** a introdução não diz que a
     condição de desenho é o preço da taxa rápida; a qualificação vai ao §3.
   - (c) A variante branca com o comparador truncado: não abrir; uma linha em
     "what this does not cover" no supp.
   - (d) **As citações de L10 na `k = 4`:** no §4.1, "group descent
     \citep{Breheny-Huang-2015}" e "the grpreg package
     \citep{Breheny-Miller-grpreg-2026}", tirando o `[VERIFICAR]` (nota de
     trabalho, sem cinza); no §2.3, "a group lasso \citep{Yuan-Lin-2006}";
     as três entradas no `references_4.bib`.
   - (e) **As hipóteses da Proposition S8.1:** acrescentar "together with
     the upper bounds in Assumptions 2 and 3", marcado, porque o `λ_n` e a
     prova usam `B_X` e `C_U` (herdado da `k = 3`).
   - (f) **Simon & Tibshirani (2012, *Statist. Sinica* 22, "Standardization
     and the group lasso penalty")** como referência da Remark S5.2 (a
     escala `‖Q_G^{1/2}θ_G‖_2`): é da revista-alvo; verificar antes, numa
     frente curta de bibliografia (L11), com a paginação de Breheny & Huang
     no periódico se o manuscrito citar seção ou página.
   - (g) Breheny & Huang ficam na seção de software do `.bib` (D25), como o
     artigo do `glmnet`.
   - (h) **Símbolos para o `notacao.md` §10:** `ℛ_1(θ; t)`,
     `ℛ_{𝒢,w}(θ; λ)` e `ℛ_{𝒢,1}(θ; t)` (`\mathcal{R}`, sem uso anterior),
     que não se confundem com o risco ideal `R_𝒢(θ; η)` do Lema 15.
   Rascunho da frase de (a), da E1.15: "When $s < 1/2$ the slow rate
   improves further: comparing with the truncated oracle
   $\theta^{*}_{G}\mathbf{1}\{\lVert\theta^{*}_{G}\rVert_{2} > \lambda w_{G}\}$
   instead of $\theta^{*}$ gives, without any design condition, the
   prediction rate of Theorem 1 whenever $s' > s/(2s+1)$ (Proposition S6.x);
   when $s \ge 1/2$ no comparator does better than
   $\lambda^{\mathcal{G}}_{n,1} \asymp (\log n/n)^{1/2}$, and the design
   condition is what the faster rate costs."

49. **~~Pendências de E5f e L11~~ aplicadas em 2026-10-05 (D55).** Texto
   original: três edições pequenas na
   `k = 4`, com o `.tex` proposto na resposta do chat principal de
   2026-10-05: (a) a frase do §3.2 depois do Theorem 1 ("When `s < 1/2` and
   `π ≥ 2`, the rate ... is attained with no design condition at all, at
   another resolution") passa a dizer que, em `s < 1/2`, a taxa de predição
   do Theorem 1 vale sem condição de desenho para todo `π` (Proposition
   S6.2); (b) Simon & Tibshirani (2012) no §4.1, ao lado de Breheny & Huang,
   e não no parágrafo do Theorem 2; (c) a frase de Simon & Tibshirani na
   Remark S5.2 (eq. 1.4 e §2), sem citá-los no argumento pivotal (lá é uma
   esperança sob a nula, aqui uma cota de cauda); a entrada no
   `references_4.bib`. Breheny & Huang se citam por seção, não por página
   inferida. O supp cresceu de 67 para 79 páginas, o que não conta no teto
   (D21).

50. **~~Pendências de E6.1b~~ decididas em 2026-10-05 (D56, D57).** Texto original: (2026-10-05). Recomendação do chat principal:
   - (a) **A pergunta 2, a aplicação:** marylebone, como aplicação de
     estrutura, assumindo os ~2% contra o spline sintonizado em blocos. É a
     única base com um resultado substantivo que a literatura da área
     confirma (o degrau de 2003) e que o spline não dá na forma de decisão.
     Pede a licença limpa (a série MY1 do UK-AIR, OGL v3, e o vento de fonte
     declarada) e, para a data, que não é periódica, pesar
     `boundary = "interval"` (pergunta 6). As alternativas: aplicação neutra
     (bike ou kelmarsh) ou outra busca.
   - (b) ~~**O spline da aplicação:** o `gam.cv` nas mesmas dobras por
     bloco, ao lado dos dois de D46.~~ **Decidido em D56: não.** Cada método
     nos seus termos; contra os dois splines de D46, beijing.heat vence por
     ~9% (`+cv` 0,790 contra `gam.reml` 0,867, 20 de 20 partições) e
     marylebone empata a menos de um erro-padrão (`+cv` 11,79 contra 11,77);
     housing passaria na predição e não tem licença. A escolha da base (a)
     fica entre marylebone (a melhor estrutura, sem licença limpa) e
     beijing.heat (a melhor predição, CC BY 4.0, o limiar de 15/11 em 7 de 20
     partições).
   - (c) **O limiar na aplicação:** o padrão `cv1se` (D45) fica; a aplicação
     relata a predição do `+cv` e a estrutura dos dois, e só afirma ausência
     de efeito onde os dois zeram. Em marylebone, `NOx × vento` é zerado
     pelo `cv1se` e precisa dessa conferência antes de virar resultado.
   - (d) **Truncar `J` por moduladora** no número de valores distintos
     (lição 3): catalogar uma tarefa curta de código.

51. **~~Pendências de E1.9~~ decididas em 2026-10-05 (D57); a E5g leva tudo
   à `k = 4`.** Texto original: (2026-10-05). Recomendação do chat principal:
   - (a) **A introdução:** trocar "so that with a single modulator
     independent of the covariates it is minimax optimal there" por "where
     it is minimax optimal", e a frase de D50 sobre a cota inferior por "The
     lower bound that makes it optimal is proved here for any number of
     modulators and for covariates that depend on them; Klopp and Pensky
     have it with a single modulator independent of the covariates."
   - (b) **A cota no artigo:** no supp, em S6, no lugar da Remark S6.1 (que
     vira uma frase: K&P é o caso `q = 1`, `X ⊥ U`), com o enunciado nas três
     perdas e a prova; no corpo, a terceira leitura do Theorem 1 no texto do
     handoff da E1.9.
   - (c) **A uniformidade do `O_p`:** uma frase na prova do Theorem 1 no
     supp ("the event does not depend on `f` and the bound depends on `f`
     only through `C_g`"), sem a qual "minimax optimal" fica no sentido
     ponto a ponto.
   - (d) O `(pq)²` do Corolário 1 de E1.3 fica como está.
   - (e) Os símbolos do `09` ficam locais; só os que chegarem ao manuscrito
     (`ρ̲_n`, `𝒦`, `K`, `c_A`) entram no `notacao.md` §10.
   - (f) **Catalogar a E5g** para levar (a) a (c) à `k = 4`, com
     `Tsybakov-2009` no `references_4.bib`.

   Texto da terceira leitura do Theorem 1, do handoff da E1.9 (guardado
   aqui porque o handoff foi apagado na integração); "Theorem S6.1" é
   provisório e depende de onde a cota entrar no supp:

   ```latex
   \textcolor{colR1}{Third, \textcolor{gray}{\sout{$n^{-2s/(2s+1)}$ is the minimax rate over Besov bodies
   in the Gaussian sequence model with noise level of order $n^{-1/2}$
   \mbox{\citep[Theorems~4 and~5]{Donoho-Johnstone-1998}}, and with a single modulator,
   independent of the covariates, and Gaussian errors, the minimax lower bound of
   \mbox{\citet[Theorem~1]{Klopp-Pensky-2015}} holds over the class of
   \mbox{Assumption~\ref{ass:besov}} (Remark~S6.1 in the Supplementary Material), so that
   in that case \mbox{Theorem~\ref{thm:main}} is minimax optimal when $\pi \ge 2$.}} the rate
   cannot be improved. With Gaussian errors and any law of $(\bX,\bU)$ that
   satisfies Assumptions~\ref{ass:X} and~\ref{ass:U}, the minimax risk over the
   class of Assumption~\ref{ass:besov}, in the prediction norm and for the
   components, is at least a constant multiple of $pq\,n^{-2s/(2s+1)}$, for every
   $\pi$ (Theorem~S6.1 in the Supplementary Material). This is the regression
   counterpart of the lower bound of \citet[Theorem~5]{Donoho-Johnstone-1998} for
   Besov bodies in the Gaussian sequence model, and with a fixed number of
   covariates it extends that of \citet[Theorem~1]{Klopp-Pensky-2015}, who have a
   single modulator independent of the covariates. Theorem~\ref{thm:main} is
   therefore minimax optimal when $\pi \ge 2$, with any number of modulators and
   with covariates that depend on them. When
   $\pi < 2$ the remaining factor $(\log n)^{(2/\pi-1)/(2s+1)}$ is the one in the
   rate of \citet[eq.~3.19]{Klopp-Pensky-2015} with a fixed number of covariates,
   and in that of block thresholding in the sequence model
   \citep[Theorem~4]{Cai-1999}; \textcolor{gray}{\sout{whether it is necessary for this estimator is not
   known, and no lower bound is proved here for $q \ge 2$.}} the lower bound carries
   no such factor, and whether the factor is necessary for this estimator is not
   known.
   ```

52. **~~Pendências de E3.4 e E5g~~ decididas e aplicadas em 2026-10-05
   (D58).** Texto original: (2026-10-05). Recomendação do chat
   principal:
   - (a) **Ratificar a interface de E3.4** (`cap.J = TRUE` padrão,
     `design$J` efetivo e `J.requested`, `cvfit$J.eff`, o reaproveitamento
     do candidato repetido, a contagem no círculo): sim.
   - (b) Moduladora de um ponto no círculo (dois valores com `eps = 0`):
     aviso, não erro; sim.
   - (c) Moduladora periódica discreta (a hora) como `u/24` com
     `rescale = FALSE`: fica para quando uma aplicação usar uma; nenhuma das
     de D57 usa.
   - (d) "The lower bound that makes it optimal" → "that makes the rate
     optimal" (o "it" fica longe do antecedente): sim.
   - (e) **A cota inferior nas contribuições e no resumo:** na primeira
     contribuição, ", together with a minimax lower bound, for any number of
     modulators, that makes these rates optimal when `π ≥ 2` (Theorem S6.1
     in the Supplementary Material)"; no resumo (D53), depois de "with no
     logarithmic factor when the Besov integrability index is at least
     two", a oração ", which is minimax optimal". Recomendação: as duas, por
     ser o resultado que fecha o posicionamento de D18.
   - (f) A condição `s ≥ 1/π` da taxa de Cai (1999), que saiu com a Remark
     S6.1: uma oração no parágrafo depois do Theorem S6.1 (a condição do
     Theorem 1, `s' > s/(2s+1)`, é mais fraca); sim.
   - (g) O `ω` do hipercubo (S6) e o `ω_n` do Corollary 2(iii) (S7) ficam,
     em seções diferentes.
   - (h) O espaço que falta entre "Lemma 3." e "no condition" na segunda
     leitura do Theorem 1 (E5d): corrigir, sem marca, por ser tipografia.
   - (i) Uma nota no `notacao.md` §10 sobre `|𝒦|` e `|𝒞|` no supp.

53. **~~Pendências de E6.1c~~ decididas em 2026-10-06 (D59); L12
   catalogada.** Texto original: (2026-10-06). Recomendação do chat principal:
   - (a) **O vento na aplicação:** relatar o efeito do vento sobre a fração
     primária como modulação pequena (~17% do degrau), que o `cv1se`
     descarta e o `+cv` e o spline mantêm, numa frase; é a regra de D57 em
     uso. A tese da aplicação fica com o degrau e o ciclo anual.
   - (b) **A borda da data:** a base periodizada cola o fim da série no
     começo, e a componente cai no último mês. Para a figura, cortar a grade
     antes da borda, com uma nota; medir `boundary = "interval"` só se um
     referee pedir (a rota está em `07-rota-intervalo.tex`).
   - (c) **O mergulho do começo de 2002 e o pico do fim de 2003**, que o
     WAFC e o spline mostram nas duas bases: descrever em uma frase, sem
     interpretar, até ler Carslaw (2005) no PDF.
   - (d) **Uma rodada bibliográfica (L12)** com as fontes da aplicação:
     Carslaw (2005), Carslaw & Beevers (2005), Clapp & Jenkin (2001),
     Hersbach et al. (2020), o conjunto ERA5 do C3S, o Open-Meteo, o UK-AIR
     com o relatório de conversão da Defra, e, se o artigo os usar, Nadeau &
     Bengio (2003) e Opsomer, Wang & Yang (2001); Carslaw (2005) lido no PDF
     para conferir o degrau e os dois ressaltos.

54. **~~O desenho do estudo de simulação (E4.2)~~ decidido em 2026-10-06
   (D60), com o topo das grades em `n = 2000` acrescentado e o compêndio
   numa pasta do `wafc-draft`.** Texto original: (2026-10-06), proposta do
   chat principal, à espera do autor; o texto completo está no
   `plano-projeto.md`, E4.2 ("Proposta revista de 2026-10-06"). Os pontos a
   decidir:
   - (a) **Réplicas:** 100 (recomendado; as frações de vitória e as caudas
     pedem mais que as 50 do piloto) ou 50. O núcleo com 100 custa ~290 h de
     processador, ~37 h em 8 núcleos; com 50, a metade.
   - (b) **Os braços**, um fator por vez em `n = 1000` nas células não
     homogênea e `mixed`, ~3 h cada: ruído alto (razão sinal-ruído 1),
     moduladoras correlacionadas (`u_rho = 0,5`), e `X` dependente de `U`
     (pede uma opção nova no `dgp.R`, tarefa curta de código).
   - (c) **O braço de escala:** `(p, q) = (6, 4)` com 6 blocos ativos e a
     grade `2:7`, com o custo medido na E4.3 antes de decidir; o
     `(10, 5)` do plano antigo não cabe (12 750 colunas em `J = 8`).
   - (d) **O VCBART:** é bayesiano (soma de árvores por MCMC, média a
     posteriori), e o autor não sabe se faz sentido compará-lo aos
     frequentistas. **Fica no estudo por curiosidade; a decisão de usá-lo no
     artigo vem depois** (2026-10-06).
   - (e) **O repositório `wafc-studies`:** o autor o cria no GitHub
     (privado) antes da E4.1.

55. **~~Pendências de E4.1b~~ decididas em 2026-10-06 (D62).** Texto original: (2026-10-06). Recomendação do chat principal:
   - (a) **A variância no `x_dist = "uniform"`** salta de 1/3 em `ρ = 0`
     (mantida pela reprodução bit a bit) para 1 em `ρ > 0`. O estudo usa o
     gaussiano, então não pesa: **manter**, com o salto registrado no
     roxygen, como já está.
   - (b) **`B_X`:** a hipótese sobre `X` pede `X` limitada, e o `Z`
     gaussiano não é limitado em célula nenhuma do estudo, não só no braço.
     **Recomendação:** a Seção 5 diz numa frase que as covariáveis
     gaussianas ficam fora da limitação suposta pela teoria e que o braço
     exercita a parte dos autovalores, que é a que D13 generaliza; não
     trocar o gerador.
   - (c) `ρ` negativo fica fora (só inverte `h`).

56. **~~Pendências de E4.1~~ decididas em 2026-10-06 (D63); E4.1c catalogada.** Texto original: (2026-10-06). Recomendação do chat principal:
   - (a) **Os oráculos:** a E4.1 leu "o linear oráculo" como as duas
     colunas do piloto, `linear` (coeficientes constantes, a referência do
     nulo) e `oracle` (o LASSO coordenado restrito aos blocos verdadeiros).
     **Recomendação:** ficam os dois, e o `oracle` passa à forma em blocos
     (o block LASSO balanceado restrito aos blocos verdadeiros), porque
     desde D44 o estimador é o block LASSO e o oráculo coordenado mede o
     custo da seleção de outro estimador; o coordenado fica no código
     (D43). Pede `penalty` no `wafc_fit_oracle()`.
   - (b) **O braço de ruído:** o plano dizia "razão sinal-ruído 1 contra o
     padrão 4", mas o padrão das células não homogênea e `mixed` é 3; a
     E4.1 usou 1 contra 3. **Recomendação:** aceitar e corrigir o plano.
   - (c) **O braço de escala:** o padrão do não homogêneo duas vezes, em
     `(x1, x2) × (u1, u2)` e `(x3, x4) × (u3, u4)`, com `β5` e `β6`
     constantes, `snr = 3`, `gam.gcv` fora; hoje a célula é montada no
     compêndio, porque o `simulate_wafc()` sempre ativa 3 blocos.
     **Recomendação:** aceitar, com um argumento `structure` no
     `simulate_wafc()` e o compêndio só repassando-o.
   - (d) **O lasso no supp:** as duas linhas, `lasso` e `lasso.cv1se`
     (saem do mesmo ajuste, sem custo). **Recomendação:** ficam as duas.
   - (e) **Licença:** `GPL (>= 3)` no `DESCRIPTION` e no README, sem
     arquivo `LICENSE`; o `WaveBased` é `GPL (>= 2)`, compatível. **Decisão
     do autor**; a recomendação é manter e pôr o `LICENSE` na extração.
   - (f) O custo no README entra depois da E4.3.
   - **(a) e (c) juntas são uma tarefa curta de código** (E4.1c: `penalty`
     no `wafc_fit_oracle()`, `structure` no `simulate_wafc()`, com o padrão
     reproduzindo os sorteios, e o compêndio repassando-os), a catalogar se
     aceitas; a E4.3 espera por ela.

57. **~~Pendências de L12~~ decididas em 2026-10-06 (D64).** Texto original: (2026-10-06). Recomendação do chat principal:
   - (a) **A frase de D59 sobre os ressaltos:** o artigo descreve os dois
     numa frase, dizendo que o mergulho do começo de 2002 aparece também
     nas estimativas mensais de Carslaw (2005, Fig. 3(a)) e que a série
     dele termina antes do pico do fim de 2003, sem interpretar nenhum dos
     dois. **Recomendação:** aceitar; fecha a ressalva "até a L12 ler
     Carslaw" de D59.
   - (b) **O número de comparação do degrau:** Marylebone Road (~10 a ~23
     vol%), não a média de Londres; e a inclinação como estimativa da
     fração primária, com a qualificação de Clapp & Jenkin. **Recomendação:**
     aceitar; vale para a E6.2, a E5b e o comentário do script 10 (linhas
     400–405), corrigido pela próxima tarefa que tocar o script; o
     `aplicacao-candidatas.md` ganhou nota nas duas linhas.
   - (c) **A forma do ERA5:** uma entrada só, a de 15 autores (2023), que o
     CDS dá para os dados, com o mesmo título e DOI da entrada do catálogo.
     **Recomendação:** manter.
   - (d) **O ano do Open-Meteo:** 2024, o do registro no DataCite, como no
     `Wood-mgcv-2025`. **Recomendação:** manter.

58. **Pendências de L13** (2026-10-06); **decidida em 2026-10-06 (D65 e D66)**. Recomendação do chat principal:
   - (a) **~~O Opsomer na aplicação~~ recusado (D65):** na `beijing.heat`, a vantagem de ~9%
     do WAFC sobre os splines de D46 some contra o `gam.cv` em blocos
     (E6.1b, §13), e Opsomer et al. dão a razão (o seletor do spline sob
     erro dependente, não a base). **Recomendação:** a segunda frase acima
     entra na aplicação, junto da tabela da `beijing.heat`, para que a
     vantagem seja lida como de sintonia em dado dependente antes que o
     referee a descubra; o `gamm` com erro autorregressivo fica para a
     resposta ao referee (risco de D56). É coerente com D56: o artigo diz
     que nenhum dos métodos modela a dependência e diz o que isso faz com
     o seletor do spline.
   - (b) **As marcas da `beijing.heat`:** ficam nos dias nominais (318 e
     73), com "nominal" na legenda e nos comentários do código (o
     `wafc-studies/R/application_data.R`, linhas 70–74, e o script 10,
     linhas 351–357, dizem "administrative" sem a ressalva).
     **Recomendação:** ficar no nominal; datas efetivas pediriam fonte
     para 2013 a 2016 que nenhum artigo lido dá.

59. **~~Pendências de E6.2~~ decididas em 2026-10-06 (D67).** Texto original: (2026-10-06). Recomendação do chat principal:
   - (a) **As rodadas** (esperam o aviso do autor): (i) a junção do WAFC da
     `beijing.heat`, ~2,6 h num processo, no retrato; (ii) as 20 partições
     e a amostra inteira, 168 unidades, ~70 h de processador, ~24 h em 3
     processos de até ~5 GB, **na árvore limpa** (unidade com mudança não
     commitada não vai ao artigo). **Recomendação:** (i) e (ii) juntas
     agora, 4 processos; a E4.3 depois, sozinha na máquina, porque mede
     tempo e o tempo é métrica do estudo; a E4.4 depois da E4.3.
   - (b) **As partições** do compêndio vêm do mestre do estudo e não são as
     da E6.1b e da E6.1c (o catálogo dizia "as mesmas"); os números do
     artigo vão diferir um pouco dos da §17. **Recomendação:** manter as
     independentes, pela razão da E4.1: marylebone e o `+cv` foram
     escolhidos olhando aquelas partições.
   - (c) **A junção do WAFC** a ~1e-14 em `t` e nos candidatos do limiar,
     exata em toda estimativa e predição. **Recomendação:** aceitar como
     exata.
   - (d) **A figura** com `β_ℓ` remontado e o corte de 60 dias nas duas
     pontas da data (também nas leituras). **Recomendação:** manter.
   - (e) **A borda nas outras moduladoras** (a umidade da `beijing.heat`, o
     vento fraco de marylebone). **Recomendação:** decidir depois da
     rodada, cortando com a mesma nota se a dobra aparecer com o `J`
     escolhido (uma linha em `figure$trim`).
   - (f) **A linha `one` da `beijing.heat`** (o intercepto a 0 °C e 1000
     hPa, que vai de 0 a 5,6 em log PM2.5 ao longo do ano). **Recomendação:**
     decisão de apresentação da E5b; o supp mostra vento, temperatura e
     pressão, e a do nível sai ou vai com uma frase.

60. **~~Pendências de E5h~~ decididas em 2026-10-06 (D68); L14 catalogada.** Texto original: (2026-10-06). Recomendação do chat principal:
   - (a) **"The four rules are compared in the simulation study"** (§4.3):
     o estudo de D60 só roda a validação cruzada do WAFC. **Recomendação:**
     tirar a frase; citar números da regra da teoria pediria que o
     compêndio os reproduzisse (D60(f)).
   - (b) **O braço da cópula está fora da Assumption 3** (densidade
     ilimitada em dois cantos), e o texto diz isso numa oração.
     **Recomendação:** manter a frase e o braço (D60).
   - (c) **O KP não é a forma de Klopp & Pensky** (pedaços dentro do nível,
     pesos do `grpreg`), e o texto diz as duas diferenças.
     **Recomendação:** manter.
   - (d) **`Z_ℓ` na frase do braço** colide com a matriz de desenho `\bZ`.
     **Recomendação:** trocar por `ζ_ℓ` na E5b.
   - (e) **A UCI pede citar Zhang et al. (2017, *Proc. R. Soc. A* 473,
     20170457)**, e nem ele nem o conjunto estão no verificado.
     **Recomendação:** uma L14 curta com as duas entradas.
   - (f) **"independent of the exploratory partitions on which this data
     set was chosen"** (§6.2). **Recomendação:** manter, como a frase das
     sementes do piloto no §5.1.
   - (g) **A linha do nível fora da Figura S1, com uma frase (D67(f)), e a
     legenda da Figure 2 encurtada.** **Recomendação:** manter as duas.
   - (h) **O teto:** 64 páginas sem o removido contra 40, sem a Discussion.
     **Recomendação:** decidir no fim (D21); o que vai ao supp sem
     reescrever é a Table 1 e as Figures 1 e 2 (uma linha cada), e no
     texto encurtam o parágrafo do `s'`, a descrição dos concorrentes e o
     dos braços.
   - (i) **Menores, na E5b:** "cross-validation" com hífen em todo o texto
     (o §4.3 usa sem hífen); o resumo ganha a simulação e a aplicação;
     "Hersbach et al., 2020, 2023" são grupos distintos que o
     `chicago.bst` junta pelo primeiro nome, o que é o padrão do `natbib`
     e fica.

61. **Pendências de E4.3a** (2026-10-06); **(a), (c) e (d) sem objeto depois de D69; (b) e (e) decididas em 2026-10-06 (D70)**. Recomendação do chat principal:
   - (a) **O gatilho de 20%** é lido na **pior célula com componentes**,
     não na fração conjunta: se a grade aperta numa célula, o estimador ali
     está limitado por ela; a leitura mostra as duas. **Recomendação:**
     aceitar.
   - (b) **Os ~10 h da escala (D60(c)) são de relógio em 8 núcleos**, como
     os "~3 h" de cada braço; a estimativa é ~25 a 30 h de processador,
     ~3,5 h de relógio. **Recomendação:** aceitar; a escala entra.
   - (c) **10 réplicas em `n = 2000`**, estendidas a 20 só no método que
     ficar a menos de 0,1 do gatilho (as 10 primeiras vêm do cache).
     **Recomendação:** aceitar.
   - (d) **O custo da grade subida é estimativa** (o fator 2 do `J` é
     conjectura; no `bsgl` pode chegar a 4); se a regra disparar, as
     primeiras unidades da produção em `n = 2000` dão o número.
     **Recomendação:** aceitar, sem medição à parte.
   - (e) **A escala no mestre da produção, réplicas 991 a 993.**
     **Recomendação:** aceitar; a decisão dela é de custo, e os dados ficam
     fora da produção.

62. **~~Pendências de L14~~ decididas em 2026-10-06 (D71).** Texto original: (2026-10-06). Recomendação do chat principal:
   - (a) **A troca de estação do Dongsi vai à S10**, na oração da frase
     acima, só no supp: custa uma linha, não pesa no teto, e antecipa o
     referee que abrir a Tabela S1 de Zhang et al. **Recomendação:**
     aceitar.
   - (b) **"described by" no `wafc-studies/data-raw/sources.yaml`** (a
     frase de licença): a UCI chama o artigo de "Introductory Paper".
     **Recomendação:** passar a "introductory paper:" na próxima vez que
     alguém tocar o arquivo.

63. **~~Pendências de E5i~~ decididas em 2026-10-06 (D72)**, todas para a E5b. Texto original: (2026-10-06). Recomendação do
   chat principal:
   - (a) **A frase que sobra no §4.3** ("the difference between the two is
     a finite sample quantity worth measuring") fica sem quem a meça.
     **Recomendação:** encurtar para "Its interest is as a reference for
     what the cross-validated choice ought to approach."
   - (b) **"$K$ fold"** sem hífen, ao lado de "10-fold". **Recomendação:**
     "$K$-fold", marcado.
   - (c) **"the grids below were settled"** no fim do §5.1 generaliza
     depois de D69. **Recomendação:** "the grids below up to $n = 1000$".

65. **~~Pendências de E5j~~ decididas em 2026-10-07 (D75); E6.2b catalogada.** Texto original: (2026-10-07). Recomendação do chat principal:
   - (a) **Os topos das grades na aplicação:** os dois splines escolhem
     `k = 80`, o topo, nos 21 ajustes de marylebone, e o `k.check()` os
     marca; o WAFC do ajuste inteiro escolhe `J = 8`, o topo. A comparação
     de predição é contra splines presos na grade. Com `n` entre 34 000 e
     61 000, muito acima dos 2 000 em que D69 já estendeu a grade do
     estudo, a grade da aplicação tem o mesmo problema. **Recomendação:**
     uma sondagem curta fora do compêndio (`k = 120` e `160` nos dois
     splines e `J = 9` no WAFC, no ajuste inteiro e em poucas partições das
     duas bases) antes de decidir; se o erro ou a estrutura mudarem, a
     grade da aplicação cresce e as 168 unidades rodam de novo (~1 dia); se
     não mudarem, fica, com uma frase.
   - (b) **As colunas `range_` e `jump_` de (nox, day)** no ajuste inteiro
     do WAFC saem do buraco de 1998 (57,2 contra 20,8 no spline).
     **Recomendação:** o relatório as calcula na grade desenhada, sem
     pontas e sem buracos (muda só o `structure-`).
   - (c) **As leituras incluem os pontos da grade nos buracos;** tirá-los
     muda a subida do ajuste inteiro de 8,88 para 8,96 (o texto passaria a
     9,0). **Recomendação:** excluir também nas leituras, para a figura e
     os números dizerem o mesmo.
   - (d) **As legendas dos floats novos saem pretas** na versão marcada.
     **Recomendação:** `\color{colR1}` dentro de cada `figure` e `table`
     novos da `k = 5`, na E5b.
   - (e) **A Figura S1 sem rug.** **Recomendação:** `rug: {doy: 1, rh: 1}`
     para o mesmo aspecto, na E5b.
   - (f) O comentário no topo do `app-marylebone.tex` ainda tem colchetes
     `[E6.2: ...]` do rascunho. **Recomendação:** limpar na E5b.

66. **~~Pendências de E6.2b~~ decididas em 2026-10-07 (D76).** Texto original: (2026-10-07). Recomendação do chat principal:
   - (a) **Os splines da `beijing.heat`** não foram sondados com `k` maior.
     **Recomendação:** não rodar; em marylebone mais `k` piorou a predição,
     e o seletor no topo sob erro dependente é o mesmo mecanismo.
   - (b) **A frase da data de 90%** na `k = 5` ("by the end of July 2003
     ... in every partition") não se sustenta nos pontos desenhados.
     **Recomendação:** na E5b, "90% by July 2003 in the fit to the whole
     sample and in 18 of the 20 partitions, and by the end of the gap in
     the data of late July and early August 2003 in the other two", com
     9,0 e 13,5 nos lugares de 8,9 e 13,4.
   - (c) **A frase da grade** (acima). **Recomendação:** no supp, junto da
     tabela de escolhas, com uma oração no §6.3, por causa do teto.

67. **~~O piloto e a grade dos splines~~ decidida em 2026-10-07 (D77).** Texto original: (2026-10-07). Recomendação do chat
   principal:
   - (a) **Os splines no topo da grade** no não homogêneo e na `mixed` em
     `n ≥ 1000`, mesmo com `k = 120` em `n = 2000`: a comparação nessas
     células pode medir a grade do spline. Com erro independente, ao
     contrário da aplicação, mais `k` deve ajudar o spline. **Recomendação:**
     uma sondagem de ~1 a 2 h (máquina livre agora), com critério
     declarado antes: o `gam.reml` com `k` até 240 em `n = 1000` e `2000`,
     no não homogêneo e na `mixed`, 5 réplicas do mestre do piloto; se o
     erro (`rmse_f`) cair mais de 1% em média, a grade do `gam.reml` cresce
     nessas `n` (por `n`, não por célula, D41) antes da E4.4; senão fica,
     com uma frase. O `gam.gcv` fica na grade atual pelo custo (expoente
     ~3,7; a razão é a mesma de D41 para tirá-lo da `mixed`), com o topo
     dito no texto.
   - (b) **A produção (E4.4)**: ~59 h de relógio sozinha em 8
     trabalhadores (~2,5 dias), com o aviso do autor, depois de (a).

68. **Pendências de E5k** (2026-10-08). **Respondida em 2026-10-12: o
   autor aceitou as nove recomendações.** A E5k já tinha aplicado (a),
   (b), (c), (e), (f) e a parte de Marylebone do (d) na `k = 5`; (g) e (h)
   ficam; a `k = 6` acrescentou a frase de Dongsi do (d) na S10 e as três
   linhas do tempo por ajuste do (i) na Table S7
   (`tables/app-tuning_6.tex`, da coluna `time` de
   `results/tables/choices-<app>.csv`). Recomendação do chat principal:
   - (a) **"by the end of July 2003"** no lugar do "by July 2003" de D76(b):
     as 18 partições passam de 1 a ~20 de julho. **Recomendação:** aceitar.
   - (b) **O REML "90% of it by the end of that gap"** (12 de agosto, o
     primeiro ponto depois do buraco). **Recomendação:** aceitar.
   - (c) **9,3 → 9,4** no mínimo das partições (9,355). **Recomendação:**
     aceitar.
   - (d) **A frase da grade com o escopo da sondagem** ("At Marylebone
     Road", "in the first three partitions", "in the fit to the whole
     sample") **e uma frase nova na S10:** "At Dongsi, where the splines
     also choose $k = 80$, the grids were not extended." (D76(a)).
     **Recomendação:** aceitar as duas; o referee vai perguntar.
   - (e) **A S6 dividida em S6 e S7.** **Recomendação:** aceitar.
   - (f) **A legenda da Figura S1** no padrão da Figure 2, com os dois
     cortes ditos e o achado no texto. **Recomendação:** aceitar.
   - (g) **A coluna de `R²` na Table S5.** **Recomendação:** fica (é supp e
     situa o erro).
   - (h) **A oração interpretativa da S10** sobre os sinais dos contrastes.
     **Recomendação:** fica; é a leitura literal dos números.
   - (i) **O tempo dos ajustes** da aplicação (o WAFC 1 800 a 8 800 s por
     ajuste contra 7 a 470 s dos splines). **Recomendação:** uma linha na
     Table S7, como a simulação dá na S4; o custo é a primeira pergunta de
     um referee, e as unidades rodam num núcleo (piloto: processador sobre
     tempo 1,00).

69. **Pendências de E4.3b** (2026-10-08). Recomendação do chat principal:
   - (a) **A margem do limite de 80 iterações:** o máximo observado foi 51
     em 523 ajustes; a produção faz ~19 000. **Recomendação:** ficar com
     80 e, no meio da E4.4, ler no `gam_k.csv` a contagem de "not
     converged" por célula; se aparecer fora do `urho` e da `mixed` em `k`
     grande, parar e reabrir antes do fim (subir para 100 custaria ~1 000 s
     por ajuste que deriva).
   - (b) **A E4.4:** pronta para rodar, ~74 h de relógio sozinha em 8
     trabalhadores (~3 dias), com o aviso do autor.

---

