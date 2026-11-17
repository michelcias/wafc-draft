# Numeração dos resultados
Movido de docs/TAREFA.md §5 em 2026-11-16, sem alteração de linhas.
Referência viva: o mapa da numeração global dos resultados; movido do TAREFA.md §5 em 2026-11-16.

<!-- conteudo movido abaixo -->
## 5. Numeração dos resultados

Os `derivations/` numeram os resultados em sequência global (Proposição 1,
Lema 1, ...), começando em E1.2. O manuscrito terá numeração própria, por
rótulo do LaTeX. Resultado novo numera a partir do último registrado aqui
pelo chat principal, com o mapa abaixo.

| Arquivo | Resultados, na numeração global |
|---|---|
| `01-identificabilidade.md` (E1.2) | Proposição 1, Lema 1 |
| `02-aproximacao-besov.tex` (E1.3) | Lema 2, Lema 3, Corolário 1, Proposição 2 |
| `03-desenho-produtos.tex` (E1.4) | Proposição 3 |
| `04-oraculo.tex` (E1.5) | Lema 4, Lema 5, Lema 6, Lema 7, Teorema 1, Corolário 2, Corolário 3 |
| `05-taxas.tex` (E1.6) | Lema 8, Lema 9, Proposição 4, Teorema 2, Corolário 4, Corolário 5 |
| `02-aproximacao-besov.tex` (emenda E1.3b) | Lema 10 |
| `07-rota-intervalo.tex` (E1.8) | Lema 11, Proposição 5, Corolário 6, Proposição 6, Corolário 7, Lema 12 |
| `06-selecao-limiar.tex` (E1.7c) | Lema 13, Corolário 8 (no manuscrito: Lemma S7.1 e Corollary 2, com a Hipótese S como Assumption 6; mapa no cabeçalho do `ms_4.tex`) |
| `08-blocos.tex` (E1.12) | Proposição 7, Lema 14, Teorema 3, Corolário 9, Lema 15, Teorema 4, Corolário 10, Corolário 11; Hipótese B (letra, como a S); reimprime o Lema 4 em blocos sem número novo |
| `06-selecao-limiar.tex` (adendo de E1.12) | Lema 16, Corolário 12, Corolário 13 |
| `08-blocos.tex` (E1.13) | Corolário 14 (`cor:lenta-blocos`) |
| `05-taxas.tex` (E1.15) | Proposição 8 (`prop:truncada`) |
| `08-blocos.tex` (E1.15) | Corolário 15 (`cor:truncada-blocos`) |
| `09-cota-inferior.tex` (E1.9) | Lema 17, Teorema 5, Corolário 16 (a Proposição 9 ficou livre) |

**Este mapa é a autoridade**, e é ele que E5a usa ao montar o manuscrito. O
`04-oraculo.tex` já imprime o número global (via `\setcounter` no
preâmbulo), prática adotada daqui em diante; `02` e `03` ainda imprimem o
contador local. O próximo resultado novo é **Proposição 9**, **Lema 18**, **Teorema 6**
ou **Corolário 17**, conforme o tipo. O `06-selecao-limiar.tex` imprime `Hipótese S`, com letra
em vez de número, porque é citada lado a lado com a `Hipótese 1` de E1.6;
é desvio local e aceito.

**O enunciado que vai ao resumo do artigo é o Corolário 11** (`08-blocos.tex`;
D47, 2026-10-03, que emenda D16, onde era o Corolário 5).

A notação está congelada (E1.1): `ψ_{jk}` com nível `j` e translação `k`,
covariável linear `X_ℓ`, moduladora `U_m`, coeficiente `θ_{ℓm,jk}`. As
macros de `derivations/macros.tex` implementam isso e são o único lugar onde
uma macro é definida.
