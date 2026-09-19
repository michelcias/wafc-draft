# Seleção de estrutura: as três saídas para L2f

Documento de decisão, escrito em 2026-09-19 a pedido do autor. A pergunta é
o que fazer com **E1.7**, a etapa condicional que prometia "seleção de
estrutura": dizer quais moduladoras afetam quais coeficientes, isto é,
identificar os blocos `(ℓ, m)` com `g_{ℓm} ≡ 0`.

A pergunta ficou em aberto porque L2 mostrou que a ideia não é nova
(`busca-novidade.md` §5): existe em splines (Wei, Huang & Li 2011; Xue & Qu
2012; Ma, Carroll, Liang & Xu 2015) e em wavelets **sem** o `X_ℓ`
multiplicando (Amato et al. 2022, seleção bi-nível; Yu et al. 2019, sparse
group LASSO). Novo só seria "no desenho de produtos".

---

## O que já está na mão

| Peça | Onde | O que dá |
|---|---|---|
| Taxa por componente | Corolário 4 de `05-taxas.tex` | `Σ_{ℓm} ‖ĝ_{ℓm} − g_{ℓm}‖²_{L_2}` na taxa de E1.6 |
| Norma por bloco, no código | `wafc_blocks()` (E2.2) | não nulos e norma `ℓ_2` de cada bloco `(ℓ, m)` |
| Contraste medido | §2 do handoff de E2.2, integrado no `ESTADO.md` | com ~99 não nulos, o LASSO zera **0** dos 3 blocos nulos e o sparse group LASSO zera **3**; o preço é predição (`RMSE` 0.2633 contra 0.2078) e tempo (3,7 s contra 0,10 s) |
| Regras que compram estrutura | §2 do handoff de E2.3 | `lambda.1se`, BIC e EBIC zeram o cenário nulo inteiro, onde `cv.min` deixa 1 ou 2 coeficientes de lixo |

Nenhuma dessas peças é um teorema de seleção. A teoria de E1.5 e E1.6 é do
**LASSO puro**; a variante em grupos, que é a que seleciona bem, hoje não
tem teoria nenhuma neste projeto.

---

## Saída (a): a propriedade oráculo de seleção, como o plano previa

Provar que o sparse group LASSO acerta o conjunto de blocos ativos com
probabilidade tendendo a 1.

- **O que exige:** uma condição de irrepresentabilidade (ou similar) no
  desenho de produtos **com termo cruzado entre moduladoras**. E1.4 entrega
  cotas de autovalor, não irrepresentabilidade; é hipótese de outra
  natureza, que não se deduz do que está provado.
- **Custo:** comparável a E1.4 mais E1.5 somadas, e mais arriscado, porque
  a condição pode simplesmente falhar no desenho de produtos e só se
  descobre tentando.
- **Novidade:** média. A técnica é padrão uma vez que se tenha a condição;
  o que é novo é o desenho.
- **Quando faz sentido:** se a seleção de estrutura for elevada a
  contribuição principal, o que hoje ela não é (`alvo-revista.md` §4 a lista
  em quarto lugar, e como software).

## Saída (b): variante sem teorema

Manter a variante em grupos como recurso computacional, com os números de
E2.2 e E2.5, e dizer no artigo, em uma frase, que seleção consistente não é
objetivo.

- **O que exige:** nada além do que já existe.
- **Custo:** zero.
- **Novidade:** nenhuma; é a opção de não fazer.
- **Risco:** é a resposta mais fraca à pergunta "e a seleção de estrutura?"
  da lista do referee (`alvo-revista.md` §4), e deixa a variante em grupos
  no artigo sem nenhuma garantia, o que um referee pode achar estranho num
  artigo que é teórico no resto.

## Saída (c): seleção por limiarização, como corolário

Definir `Ŝ = {(ℓ, m) : ‖ĝ_{ℓm}‖_{L_2} > t_n}` e provar `P(Ŝ = S) → 1`.

- **O que exige:** o Corolário 4 (que já está provado) mais uma **condição
  de separação** em nível de bloco: o menor bloco ativo tem norma acima da
  taxa, `min_{(ℓ,m) ∈ S} ‖g_{ℓm}‖ ≫ r_n`, com `t_n` entre `r_n` e essa
  separação. **Não precisa de irrepresentabilidade nem de condição de
  desenho nova.**
- **Custo:** um corolário e uma conferência numérica; ordem de E1.2.
- **Novidade:** o enunciado em si é clássico ("estimar bem e limiarizar
  seleciona"); o que é novo é valer no desenho de produtos e **para o LASSO
  puro**, que é o estimador base (D3).
- **O que ganha de quebra:** a comparação de E2.5 deixa de ser "grupos ou
  nada" e passa a ser "grupos ou LASSO limiarizado", com `wafc_blocks()` já
  fornecendo a estatística do limiar.
- **O que custa admitir:** é mais fraco que (a). Não diz que o estimador
  seleciona sozinho; diz que estimação seguida de limiar seleciona. O
  enunciado tem de dizer isso com todas as letras, e a hipótese de separação
  tem de aparecer como hipótese, não escondida.

---

## Comparação

| | (a) oráculo em grupos | (b) sem teorema | (c) limiarização |
|---|---|---|---|
| Hipótese nova | irrepresentabilidade | nenhuma | separação em nível de bloco |
| Esforço | alto, com risco de não fechar | nenhum | baixo |
| Vale para o estimador base (LASSO) | não | — | sim |
| Responde ao referee | com teorema forte | com ressalva | com teorema condicional |
| Depende de E2.5 escolher a variante em grupos | sim | não | não |

**Recomendação do assistente:** (c). É a única que produz resultado sem
hipótese de desenho nova, cobre o estimador base em vez de só a variante, e
transforma a comparação de E2.5 numa pergunta mais interessante. (a) fica
registrada como extensão natural para um segundo artigo, se a seleção virar
o tema; (b) é o recuo se (c) não fechar.

Se (c) for escolhida, E1.7 é recatalogada como `derivations/06-selecao.tex`
mais `check/06-selecao.R`, medindo nos cenários de E2.1 a curva de acerto de
estrutura contra o limiar, e o LASSO limiarizado contra o sparse group LASSO.
