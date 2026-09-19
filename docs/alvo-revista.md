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

**Decidido em 2026-09-19 (D18):** o artigo se apresenta como **extensão de
Klopp & Pensky (2015)**, e não como modelo diferente com eles citados de
passagem. O parágrafo de posicionamento foi escrito e aprovado; E5a o usa
como está, ajustando só o nome do método (D8) e os números dos resultados:

```
Klopp and Pensky (2015) study the varying coefficient model in which a
single index modulates every coefficient. They expand each coefficient
function in an orthonormal basis, estimate the resulting array by a block
LASSO, and obtain a nonasymptotic oracle inequality, a Besov rate that
adapts to inhomogeneous smoothness, and a matching minimax lower bound,
under the assumption that the covariates are independent of the index. The
present paper carries that programme to the situation an applied problem
usually presents, in which several variables modulate a coefficient at once.
We let each coefficient be additive in the modulators, which keeps the
estimator free of the curse of dimensionality in their number, we allow the
covariates to depend on the modulators, and we penalize coefficient by
coefficient rather than in blocks, leaving the levels of the coefficients
unpenalized. Each of these costs something. Additivity produces cross terms
between distinct modulators, so the population Gram matrix no longer
factorizes as a Kronecker product and has to be bounded directly; dependence
between the covariates and the modulators replaces a marginal second moment
by a conditional one; and the unpenalized levels leave a term in the risk
that block penalization does not incur. Our main result, Corollary 5, states
that the estimator attains the rate of nonlinear wavelet approximation up to
a logarithmic factor, under the same Besov assumption that governs the
linear sieve, and it reduces to the rate of Klopp and Pensky when a single
modulator is present.
```

**Frase-tese:** estender o modelo de coeficientes variáveis esparso de Klopp
& Pensky a coeficientes **aditivos em várias moduladoras** e a desenho
**dependente** produz um estimador que (i) atinge a taxa de aproximação não
linear em wavelets a menos de um fator logarítmico, sob a mesma hipótese de
Besov que governa o sieve linear (Corolário 5, D16); (ii) é um único
problema convexo resolvido pelo `glmnet` em segundos; e (iii) seleciona,
pela variante em grupos, quais moduladoras afetam quais coeficientes.

**Contribuições a defender (em ordem de força):**

1. A teoria no desenho aditivo de produtos `X_ℓ ψ_{jk}(U_m)`: a condição de
   desenho com termo cruzado entre moduladoras e sem independência (E1.4), a
   desigualdade oráculo sem condição de cone (E1.5) e a taxa por
   compressibilidade (E1.6). O que Klopp & Pensky já têm para `q = 1` com
   `X ⊥ U` é citado, não reprovado.
2. A leitura de que **a compressibilidade não custa hipótese**: a hipótese de
   Besov já implica weak-`ℓ_τ` (Lema 9), e é isso que separa a taxa da do
   sieve linear quando `π < 2`.
3. Evidência numérica de que a adaptatividade se materializa: contra splines
   (Xue & Yang; `mgcv`), contra o spline adaptativo de Wang, Jiang & Liu
   (2024) e contra o block LASSO de K&P no mesmo desenho (E4).
4. Software: o código de `wafc/`, com a mesma interface do `wall()`, na
   forma de distribuição decidida em E3.3.

**O que o referee vai perguntar:**

- **"Why not block LASSO, as in Klopp and Pensky?"** É a pergunta que o
  posicionamento escolhido convida. A resposta: a penalidade coordenada é um
  único problema do `glmnet`, dá o corolário de compressibilidade sem
  restringir `τ`, e a comparação numérica contra os blocos está em E4. A
  variante em grupos existe em `wafc()` para quem quiser estrutura (E2.2).
- **"Onde está a cota inferior para `q ≥ 2`?"** Não existe (pergunta 12 da
  §4 do `ESTADO.md`), e sob este posicionamento ela fica mais visível: o
  artigo cita a de K&P e afirma otimalidade só onde ela vale. É o risco
  assumido da escolha.
- "Por que wavelets e não splines?" O cenário não homogêneo em E4 mais a
  taxa com `π < 2`. Se E4 não mostrar ganho, o artigo não tem razão de ser.
- "Isso não é Sardy & Ma com um `X_ℓ` multiplicando?" Não: a teoria deles é
  de otimização (L2, `busca-novidade.md` §3). A resposta fica na Seção 2.
- "Como escolhe `J`?" Validação cruzada conjunta com `λ` (E2.3), e a teoria
  diz a ordem; a regra teórica ficou dois níveis acima do melhor `J`
  empírico na conferência de E1.6, e E2.3 mede o custo disso.
- "E a seleção de estrutura?" A variante em grupos zera os blocos nulos
  (E2.2, com número); ou entra com resultado (E1.7) ou o artigo diz que
  seleção não é o objetivo.
- "Aplicação real?" E6, com efeito que varia com covariáveis e interpretação.

## 5. Estrutura-alvo do manuscrito (40 páginas, SS)

Páginas em espaço duplo no template da revista. O teto de 40 **inclui as
referências**, ao contrário do que esta tabela supunha.

| Seção | Páginas | Conteúdo |
|---|---|---|
| 1 Introduction | 2,5 | problema, o que existe, o que falta, contribuições |
| 2 Model and wavelet sieve | 4 | modelo, identificabilidade, base, desenho de produtos, estimador |
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
- [ ] **Estilo de citação:** o `chicago.bst` que o template carrega abrevia
      em "et al." a partir de três autores; a §4 das instruções manda listar
      os três. Corrigir no `.bst`, não no texto.
- [ ] **Pacote autocontido:** substituir os `\input` pelas cópias dos
      arquivos de macro, com a correspondência entre bloco copiado e arquivo
      de origem registrada.
- [ ] Autores, afiliações e e-mails na última página; agradecimentos.
- [ ] Todos os números do texto conferem com `results/tables/`.
- [ ] Cada legenda de figura enuncia o achado; unidades nos eixos.
- [ ] Referências verificadas uma a uma no `.bib`; versão oficial, não arXiv.
- [ ] Pacote `WaveBased` com versão etiquetada e citada; compêndio com DOI.
- [ ] Carta de apresentação com a frase-tese e editores associados sugeridos.
