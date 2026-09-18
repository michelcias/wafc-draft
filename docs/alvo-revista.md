# Alvo: a revista

O autor pediu um periódico **Q1 em Statistics and Probability (SCImago)**,
citando *Statistica Sinica* e *Electronic Journal of Statistics* como
exemplos, ou outro de alto nível "que tenha mais a cara da proposta". Este
documento compara os candidatos, propõe um alvo e registra o que ele pede.
A proposta é **decisão D5, a ratificar** (`ESTADO.md`).

**Estado da verificação (2026-09-18):** as páginas do SCImago devolveram 403 e
a página de instruções da *Statistica Sinica* não resolveu no DNS desta
máquina; o que está abaixo sobre elas vem de buscas na web e precisa ser
confirmado pelo autor nas páginas oficiais (tarefa E0.3). As páginas da IMS
sobre a EJS foram lidas.

---

## 1. Candidatos

| Revista | SJR Q (Stat. & Prob.) | Formato típico | Por que sim | Por que não |
|---|---|---|---|---|
| **Statistica Sinica** | Q1 em 2024 (uma fonte diz Q2; confirmar) | método + teoria + simulação + aplicação; teto de **30 páginas** de manuscrito; suplementar para provas | é onde a linhagem do modelo mora: Xue & Yang (2006) e Wei, Huang & Li (2011) saíram lá; o leitor conhece additive coefficient models e group LASSO; o padrão de artigo é exatamente o que o plano produz | processo lento (12 a 18 meses é comum); exige aplicação real convincente; o teto de 30 páginas obriga provas ao suplementar |
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
(iii) o teto de 30 páginas é folgado para o corpo com as provas no
suplementar.

**Reserva: *Electronic Journal of Statistics*.** Se a teoria ficar mais forte
do que a aplicação (por exemplo, se E1.6 der uma taxa de adaptação limpa e
E6 não achar aplicação convincente), a EJS é o veículo melhor, e a mudança
custa só o template, já versionado em `manuscript/ejs-template/`.

O plano é escrito para o formato da *Statistica Sinica* e mantém a EJS como
saída sem retrabalho: provas em `supp_{k}.tex` desde o começo; corpo com
enunciados completos.

## 3. Requisitos formais

### Statistica Sinica (a confirmar na página oficial, E0.3)

| Item | O que se sabe | Fonte |
|---|---|---|
| Extensão | "should normally not exceed 30 manuscript pages" | busca web, 2026-09-18 |
| Template | LaTeX da revista, disponível no site | idem |
| Submissão | ScholarOne (`mc04.manuscriptcentral.com/statisticasinica`) | idem |
| Suplementar | aceito (provas, tabelas extras) | prática dos artigos citados |
| Código | não é exigido; declarar disponibilidade | prática |
| Revisão | simples cega (a confirmar) | |
| Instruções | https://www3.stat.sinica.edu.tw/statistica/AUTHORS.HTM | inacessível daqui em 2026-09-18 |

Tarefa do autor em E0.3: abrir a página, transcrever as instruções em
`docs/ss-instrucoes-autores.md` e baixar o template para
`manuscript/ss-template/`.

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

**Frase-tese provisória:** em modelos de regressão com coeficientes
funcionais aditivos, expandir cada componente em wavelets e estimar todos os
coeficientes por LASSO produz um estimador que (i) adapta à regularidade
local de cada componente, com taxa quase minimax em espaços de Besov sem
suavidade global; (ii) é um único problema convexo resolvido pelo `glmnet`
em segundos; e (iii) seleciona, via a variante em grupos, quais moduladoras
afetam quais coeficientes.

**Contribuições a defender (em ordem de força, provisória):**

1. O estimador e sua teoria: desigualdade oráculo no desenho de produtos
   `X_j ψ(U_k)`, taxas em Besov e o corolário de adaptação por
   compressibilidade (E1.3 a E1.6).
2. A condição de desenho para produtos (E1.4), que é o que distingue este
   problema do modelo aditivo com wavelets de Sardy & Ma (2024).
3. Evidência numérica de que a adaptatividade se materializa: contra
   splines (Xue & Yang; `mgcv`) em funções não homogêneas, sem perder muito
   nas suaves (E4).
4. Software: `wafc()` no `WaveBased`, com a mesma interface do `wall()`.

**O que o referee vai perguntar:**

- "Por que wavelets e não splines?" A resposta é o cenário não homogêneo em
  E4 mais a taxa em Besov com `π < 2`. Se E4 não mostrar ganho, o artigo não
  tem razão de ser.
- "Isso não é Sardy & Ma com um `X_j` multiplicando?" A resposta tem de estar
  na Seção 2 (identificabilidade e desenho de produtos) e na E1.4.
- "Como escolhe `J`?" Validação cruzada, e a teoria diz qual ordem; E2.3 mede
  se BIC/EBIC serve.
- "E a seleção de estrutura?" Ou a variante em grupos entra com resultado
  (E1.7), ou o artigo diz explicitamente que seleção não é o objetivo.
- "Aplicação real?" E6, com efeito que varia com covariáveis e interpretação.

## 5. Estrutura-alvo do manuscrito (30 páginas, SS)

| Seção | Páginas | Conteúdo |
|---|---|---|
| 1 Introduction | 2,5 | problema, o que existe, o que falta, contribuições |
| 2 Model and wavelet sieve | 4 | modelo, identificabilidade, base, desenho de produtos, estimador |
| 3 Theory | 6 | hipóteses, aproximação, condição de desenho, oráculo, taxas, adaptação; provas no suplementar |
| 4 Computation and tuning | 2 | `glmnet`, desenho esparso, `(J, λ)` por CV, custo |
| 5 Simulation | 6 | desenho, competidores, métricas, resultados |
| 6 Application | 4 | uma aplicação com interpretação |
| 7 Discussion | 1 | limites, extensões (inferência, grupos, resposta não gaussiana) |
| Referências | fora do teto | |

## 6. Checklist de submissão

- [ ] Instruções oficiais relidas na semana da submissão.
- [ ] `ms` no template, dentro do teto, provas no suplementar.
- [ ] Todos os números do texto conferem com `results/tables/`.
- [ ] Cada legenda de figura enuncia o achado; unidades nos eixos.
- [ ] Referências verificadas uma a uma no `.bib`; versão oficial, não arXiv.
- [ ] Pacote `WaveBased` com versão etiquetada e citada; compêndio com DOI.
- [ ] Carta de apresentação com a frase-tese e editores associados sugeridos.
