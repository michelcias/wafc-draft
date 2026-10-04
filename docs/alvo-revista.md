# Alvo: a revista

O autor pediu um periódico **Q1 em Statistics and Probability (SCImago)**,
citando *Statistica Sinica* e *Electronic Journal of Statistics* como
exemplos, ou outro de alto nível "que tenha mais a cara da proposta". Este
documento compara os candidatos, propõe um alvo e registra o que ele pede.
A proposta é **decisão D5, a ratificar** (`ESTADO.md`).

**Estado da verificação:** o autor confirmou em 2026-09-19, no próprio
SCImago, que a *Statistica Sinica* é **Q1** em Statistics and Probability,
o que fecha o que restava de E0.3; as páginas do SCImago devolviam 403 desta
máquina. As
páginas da IMS sobre a EJS foram lidas. As instruções da *Statistica Sinica*
foram capturadas pelo autor e estão transcritas em
[`ss-instrucoes-autores.md`](ss-instrucoes-autores.md): **o teto é de 40
páginas em espaço duplo, incluídas as referências**, e não as 30 que a busca
na web tinha dado; o §3 e o §5 abaixo já refletem isso. O tipo de revisão não
consta em lugar nenhum da página oficial.

---

## 1. Candidatos

| Revista | SJR Q (Stat. & Prob.) | Formato típico | Por que sim | Por que não |
|---|---|---|---|---|
| **Statistica Sinica** | **Q1** (confirmado pelo autor no SCImago, 2026-09-19) | método + teoria + simulação + aplicação; teto de **40 páginas em espaço duplo**, referências incluídas; suplementar para provas | é onde a linhagem do modelo mora: Xue & Yang (2006) e Wei, Huang & Li (2011) saíram lá; o leitor conhece additive coefficient models e group LASSO; o padrão de artigo é exatamente o que o plano produz | triagem dura (só ~40% passa dos co-editores); exige aplicação real convincente; a revista prefere provas e detalhes de simulação no suplementar |
| **Electronic Journal of Statistics** | Q1 em 2024 | sem teto de páginas; open access sem taxa (IMS e Bernoulli); "same standard as other IMS journals" | cabe a teoria completa no corpo; sem custo; a revisão tende a ser mais rápida; template `imsart` já versionado aqui | a expectativa de teoria é maior (uma desigualdade oráculo "transposta do WALL" pode parecer pouco); menos peso para o pacote e a aplicação |
| Scandinavian Journal of Statistics | Q2 (2024/2025, a confirmar) | teoria + simulação | é onde Sardy & Ma (2024) e o próprio autor (Montoril, Pinheiro & Vidakovic 2019) publicaram; leitor natural | não atende ao pedido de Q1 |
| Journal of Computational and Graphical Statistics | Q1 (a confirmar) | método computacional + software | o pacote e a comparação numérica são fortes | a teoria ficaria em segundo plano; o WALL aplicado já mira lá, e dois artigos do mesmo grupo no mesmo veículo com o mesmo motor pedem diferenciação |
| Journal of Multivariate Analysis / Journal of Nonparametric Statistics | Q2/Q1 (a confirmar) | teoria | | abaixo do pedido, ou sem a cara aplicada |
| Bernoulli, JASA, Biometrika, JRSS-B, Annals | Q1 topo | | | o nível de novidade teórica exigido é outro; não é o alvo de um primeiro artigo do método |

## 2. Proposta

**Alvo primário: *Statistica Sinica*.** Razões: (i) a proposta é a extensão
natural de dois artigos da revista (Xue & Yang 2006; Wei, Huang & Li 2011) e
o referee provável já está lá; (ii) o formato completo (teoria, simulação,
aplicação, software) é o que o plano produz, e a revista valoriza os quatro;
(iii) o teto de 40 páginas em espaço duplo é folgado para o corpo com as
provas no suplementar.

**Reserva: *Electronic Journal of Statistics*.** Se a teoria ficar mais forte
do que a aplicação (por exemplo, se E1.6 der uma taxa de adaptação limpa e
E6 não achar aplicação convincente), a EJS é o veículo melhor, e a mudança
custa só o template, já versionado em `manuscript/ejs-template/`.

O plano é escrito para o formato da *Statistica Sinica* e mantém a EJS como
saída sem retrabalho: provas em `supp_{k}.tex` desde o começo; corpo com
enunciados completos.

## 3. Requisitos formais

### Statistica Sinica (página oficial lida em 2026-09-18)

Transcrição completa em [`ss-instrucoes-autores.md`](ss-instrucoes-autores.md);
aqui só o que decide o formato do manuscrito.

| Item | Regra | Fonte |
|---|---|---|
| Extensão | "should not exceed 40 double-spaced pages using the Statistica Sinica template", **incluídas referências e apêndice**; espaço duplo, 12pt, margens de 1 polegada | instruções, Editorial Policy e §I.1 |
| Template | obrigatório; a conformidade é conferida na recepção e o manuscrito fora do template volta na hora | Editorial Policy |
| Estrutura | resumo curto, palavras-chave em ordem alfabética, título corrente com menos de 45 caracteres, afiliações e e-mails na última página; seções e subseções numeradas, sem sub-subseções | §I.2 a §I.4 |
| Submissão | ScholarOne (`mc04.manuscriptcentral.com/statisticasinica`); fonte LaTeX só depois do aceite | Submission, Manuscript Preparation |
| Suplementar | um único PDF de até 10 MB, mesmo título e mesmos autores; seção "Supplementary Material" como última seção, antes dos agradecimentos; revisado junto | §IV |
| Código | "software producing the evidence should be available for examination as well as pertinent datasets"; evidência numérica tem de ser reprodutível | Editorial Policy |
| Triagem | co-editores e SAE cortam ~60% antes do AE; o AE manda a referee ~60% do que recebe | Editorial Policy |
| Revisão | **não declarada** na página (nem simples nem dupla cega; sem exigência de anonimização) | §7 de `ss-instrucoes-autores.md` |
| Instruções | https://www3.stat.sinica.edu.tw/statistica/ (aba *Instructions for Authors*) | capturadas em PDF pelo autor |

Template versionado em `manuscript/ss-template/` (`SS-template.tex`,
`SS-template-bib.tex` com `\bibliographystyle{chicago}`, `supp-temp_20240820.tex`
para o suplemento); os três compilam com o `latexmk` local.

### Electronic Journal of Statistics (confirmado em 2026-09-18)

| Item | Regra | Fonte |
|---|---|---|
| Template | `imsart.cls` com opção **`ejsv2`** (2025-03-20; `ejs` ainda funciona), `imsart.sty`, `imsart-nameyear.bst`; "we prefer that no separate macro or .sty files are used" | página da IMS de preparação; `vtex-soft/texsupport.ims_cosponsored-ejs` |
| Extensão | sem teto declarado | página da IMS |
| Taxas | nenhuma obrigatória; contribuição voluntária ao Open Access Fund | página da IMS |
| Licença | CC BY 4.0 | idem |
| Submissão | sistema EJMS da IMS; após aceitação, fonte LaTeX para `ejournals-ejs@vtex.lt` com o número do manuscrito | idem |
| Escopo | "research articles and short notes on theoretical, computational and applied statistics" | idem |

Template versionado em `manuscript/ejs-template/` (clonado em 2026-09-18,
`imsart.cls` 2025/03/18; `ejs-sample.tex` compila com o `latexmk` local).

## 4. Como o artigo se posiciona

**Decidido em 2026-09-19 (D18) e reescrito em 2026-10-03 (D44, D50):** o
artigo se apresenta como **extensão de Klopp & Pensky (2015)**. Com D44 o
WAFC passou a ser o block LASSO na forma balanceada seguido do limiar
`cv1se`, com o LASSO coordenado como opção, e o parágrafo de 2026-09-19
("we penalize coefficient by coefficient rather than in blocks") deixou de
valer; o texto antigo fica no histórico do git. O texto abaixo foi aprovado
pelo autor em 2026-10-03 (D50) e substitui o terceiro e o quarto parágrafos
da Introduction do `ms_3` quando a `k = 4` abrir, com a marcação do §3 do
`instrucoes.md`. Dois rótulos ainda não existem no manuscrito: "Section~5"
(E5b) e o do Corolário 13 de E1.12 (o risco do limiarizado, citado sem
rótulo na segunda contribuição).

```latex
\citet{Klopp-Pensky-2015} study the varying coefficient model in which a
single index modulates every coefficient. They expand each coefficient
function in an orthonormal basis, estimate the resulting array by a block
lasso, and obtain a nonasymptotic oracle inequality, a Besov rate that adapts
to inhomogeneous smoothness, and a matching minimax lower bound, under the
assumption that the covariates are independent of the index. The present paper
carries that programme to the situation an applied problem usually presents,
in which several variables modulate a coefficient at once. We let each
coefficient be additive in the modulators, which keeps the estimator free of
the curse of dimensionality in their number, we allow the covariates to depend
on the modulators, and we keep the block lasso but leave the levels of the
coefficients unpenalized and group the wavelet coefficients in chunks of order
$\log n$ within a resolution level, with weights whose ratio stays bounded, so
that the standardized weights of group lasso software are covered. Additivity
produces cross terms between distinct modulators, so the population Gram
matrix no longer factorizes as a Kronecker product and has to be bounded
directly; dependence between the covariates and the modulators replaces a
marginal second moment by a conditional one; and the unpenalized levels leave
a term in the risk that full block penalization does not incur. Our main
result, Theorem~\ref{thm:main}, states that the estimator attains the rate of
Klopp and Pensky, with no logarithmic factor when the Besov integrability
index is at least two, so that with a single modulator independent of the
covariates it is minimax optimal there. It holds under weaker conditions than
theirs: the number of covariates is fixed, a regime their high dimensional
condition excludes, and the effective smoothness of the components need only
exceed $s/(2s+1)$, which is below one half, where theirs must exceed one half.
Components with effective smoothness one half, such as the inhomogeneous
functions of Section~5, are covered here and not there. A lower bound for
several modulators remains open.

We call the resulting procedure WAFC, for wavelet additive functional
coefficients. Its contributions are of three kinds. The first is the theory of
the additive product design: a design condition for the Gram matrix of the
products between a linear covariate and a wavelet of a modulator, which
carries a cross term between distinct modulators and does not assume
independence (Proposition~\ref{prop:gram}); an oracle inequality that needs no
cone condition, because the unpenalized levels are profiled out, and that
holds for any weights of bounded ratio, including those of the software that
computes the estimator (Theorem~\ref{thm:oracle}); and rates that follow from
the Besov assumption alone, through a bound on the ideal risk over chunks
(Theorem~\ref{thm:main}). The second is the recovery of the structure of the
model, which modulators act on which coefficients, by a threshold on the norms
of the fitted components (Corollary~\ref{cor:threshold}): removing the null
blocks needs no separation condition, and the thresholded estimator keeps the
rate of the fit, so that one estimator serves both for structure and for
prediction. The third is computational: the estimator is one convex program
solved by standard group lasso software, the resolution level and the penalty
are chosen jointly by cross-validation, and the threshold by the one standard
error rule on the same folds.
```

**Frase-tese:** estender o modelo de coeficientes variáveis de Klopp &
Pensky a coeficientes **aditivos em várias moduladoras** e a desenho
**dependente**, mantendo o block LASSO com os níveis livres, produz um
estimador que (i) atinge a taxa deles, sem logaritmo em `π ≥ 2`, sob
condições mais fracas (`p` fixo, regularidade efetiva acima de `s/(2s+1)`;
Corolário 11 de E1.12, D47); (ii) é um único problema convexo resolvido por
software de group LASSO; e (iii) recupera, por limiarização das normas por
bloco, quais moduladoras afetam quais coeficientes, com o mesmo estimador
servindo à predição (Corolários 12 e 13).

**Contribuições a defender (em ordem de força):**

1. A teoria no desenho aditivo de produtos com a penalidade em blocos: a
   condição de desenho com termo cruzado e sem independência (E1.4), o
   oráculo sem cone para pesos de razão limitada, que cobre o estimador
   exato do `grpreg` pela variante branca (E1.12, Teorema 3), e as taxas
   pela hipótese de Besov sozinha (Lema 15, Teorema 4, Corolários 10, 11 e
   14), com o regime `s/(2s+1) < s' ≤ 1/2` fora do Teorema 2 de K&P.
2. A recuperação da estrutura por limiarização (Corolário 12) e a cota de
   risco do limiarizado (Corolário 13): triagem sem separação e um só
   estimador para estrutura e predição.
3. A computação: `cv.wafc(x, u, y)` (D48), `(J, λ)` por validação cruzada e o
   limiar pela regra de um erro-padrão.
4. Evidência numérica (E4, na E5b): a do piloto é vitória sobre o spline
   sintonizado por REML no não homogêneo e na `mixed` em `n ≥ 500` (5% a
   10% em `rmse_f`), empate com o sintonizado por GCV, derrota no suave
   perto do fator 1,5 em ISE, e estrutura recuperada em 0,9 a 1 das réplicas
   em `n ≥ 500` contra 0 a 0,4 do spline (`ESTADO.md` §2, E2.5h a E2.5j).

**O que o referee vai perguntar:**

- **"Isto não é Klopp e Pensky com mais moduladoras?"** Não só: o termo
  cruzado entre moduladoras tira a fatoração de Kronecker da Gram; `X`
  depende de `U`; os níveis são livres; os pesos de razão limitada cobrem o
  software; e o enunciado vale com `p` fixo e com regularidade efetiva entre
  `s/(2s+1)` e `1/2`, que o Teorema 2 deles exclui (`08-blocos.tex`, §7). A
  seleção de estrutura e o software também não estão lá.
- **"Por que blocos e não o LASSO coordenado?"** Mesmo com o limiar oráculo,
  o LASSO não vence o spline no não homogêneo nem na `mixed`, e os blocos
  vencem; a taxa em blocos não tem logaritmo em `π ≥ 2`. O LASSO coordenado
  fica como opção (D43), com a teoria de E1.5 e E1.6.
- **"Onde está a cota inferior para `q ≥ 2`?"** Não existe, e o artigo o diz
  na introdução (D50). Com uma moduladora independente das covariáveis a
  cota de K&P vale sobre a nossa classe, e o Corolário 11 é ótimo em
  `π ≥ 2`.
- **"Por que wavelets e não splines?"** O não homogêneo e a estrutura (item
  4 acima), e a taxa com `π < 2`. O spline por GCV empata no não homogêneo;
  a tabela mostra as duas sintonias do `gam` (D46).
- "Isso não é Sardy & Ma com um `X_ℓ` multiplicando?" Não: a teoria deles é
  de otimização (L2, `busca-novidade.md` §3).
- "Como escolhe `J` e `t`?" Validação cruzada conjunta de `(J, λ)` e a regra
  de um erro-padrão para `t` nas mesmas dobras (D45); a teoria diz a ordem
  de `J_n` e de `t_n`, e a regra teórica de `λ` é cara (§4.3 do `ms`).
- "E a irrepresentabilidade?" Não é pedida: estimação seguida de limiar, com
  triagem sem hipótese e recuperação sob separação por bloco.
- "Aplicação real?" E6 (E6.1b em curso, no critério de estrutura estável
  com predição competitiva).

## 5. Estrutura-alvo do manuscrito (40 páginas, SS)

Páginas em espaço duplo no template da revista. O teto de 40 **inclui as
referências**, ao contrário do que esta tabela supunha.

| Seção | Páginas | Conteúdo |
|---|---|---|
| 1 Introduction | 2,5 | problema, o que existe, o que falta, contribuições |
| 2 Model and wavelet approximation space | 4 | modelo, identificabilidade, base, desenho de produtos, estimador |
| 3 Theory | 6 | hipóteses, aproximação, condição de desenho, oráculo, taxas, adaptação; provas no suplementar |
| 4 Computation and tuning | 2 | `glmnet`, desenho esparso, `(J, λ)` por CV, custo |
| 5 Simulation | 6 | desenho, competidores, métricas, resultados |
| 6 Application | 4 | uma aplicação com interpretação |
| 7 Discussion | 1 | limites, extensões (inferência, grupos, resposta não gaussiana) |
| **Corpo** | **25,5** | |
| Referências e afiliações | 3 a 4 | dentro do teto |
| **Total** | **~29** | contra o teto de 40 |

**Medido em 2026-09-19, com `k = 1`:** as Seções 1 a 4 saíram em 24 páginas,
contra as 14,5 desta tabela, e as referências em 4. Por D21, a conta final é
feita quando o corpo estiver completo, contra o PDF compilado; até lá E5b
escreve sem contar, mas **cada tabela vai num `\input{}` próprio**, para que
mandá-la ao suplemento seja mover uma linha.

**Medido em 2026-09-30, com `k = 2`** (depois de E5c): o `ms` tem 33
páginas e o `supp` 30. A §3.6 (seleção por limiar) e a §4.2 somaram cerca
de 3,5 páginas ao corpo, e a projeção com as Seções 5 a 7 fica em ~44 a 46,
acima do teto; a decisão continua sendo a de D21.

A folga de cerca de 11 páginas não está alocada de propósito: se a teoria
pedir espaço, ela vai para as Seções 3 e 5, nessa ordem. A conferência é
contra o PDF compilado no template, não contra esta tabela.

## 6. Checklist de submissão

- [ ] Instruções oficiais relidas na semana da submissão (`ss-instrucoes-autores.md`).
- [ ] `ms` no template da revista, em espaço duplo, dentro das 40 páginas com
      as referências contadas; provas no suplementar.
- [ ] Título corrente com menos de 45 caracteres; palavras-chave em ordem
      alfabética; afiliações e e-mails na última página; nenhuma sub-subseção.
- [ ] Seção "Supplementary Material" antes dos agradecimentos, e o suplemento
      num único PDF de até 10 MB com o mesmo título e os mesmos autores.
- [ ] Figuras legíveis em preto e branco, sem referência a cor no texto.
- [x] ~~**Estilo de citação:** o `chicago.bst` que o template carrega abrevia
      em "et al." a partir de três autores; a §4 das instruções manda listar
      os três.~~ Dispensado pelo autor em 2026-10-01: a revista não tem
      `.bst` próprio e diagrama o artigo aceito a partir da fonte.
- [ ] **Pacote autocontido:** substituir os `\input` pelas cópias dos
      arquivos de macro, com a correspondência entre bloco copiado e arquivo
      de origem registrada.
- [ ] Autores, afiliações e e-mails na última página; agradecimentos.
- [ ] Todos os números do texto conferem com `results/tables/`.
- [ ] Cada legenda de figura enuncia o achado; unidades nos eixos.
- [ ] Referências verificadas uma a uma no `.bib`; versão oficial, não arXiv.
- [ ] Pacote `WaveBased` com versão etiquetada e citada; compêndio com DOI.
- [ ] Carta de apresentação com a frase-tese e editores associados sugeridos.
