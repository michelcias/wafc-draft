# derivations

Um resultado por arquivo, numerado na ordem de dependência do plano (E1):
`01-identificabilidade.md`, `02-aproximacao-besov.tex` (com a emenda de
E1.3b), `03-desenho-produtos.tex`, `04-oraculo.tex`, `05-taxas.tex`,
`07-rota-intervalo.tex` (E1.8). Cada um com
enunciado, hipóteses, prova, "o que isso não cobre" e o parágrafo da
conferência numérica. Documentos `article` autônomos que fazem
`\input{macros}`; `macros.tex` implementa `docs/notacao.md` e é o único
lugar onde uma macro é definida. Desde D24 ele imprime esperança e
probabilidade em **romano**, como a revista exige.

**A numeração dos resultados é global** e a tabela que manda está na §5 de
[`../docs/TAREFA.md`](../docs/TAREFA.md), não nos contadores de cada
arquivo: `02` e `03` ainda imprimem contador local, `04`, `05` e a emenda de
`02` imprimem o número global. Compilar com `latexmk -pdf <arquivo>.tex`;
o `.pdf` é versionado, os auxiliares não.

`check/` guarda os scripts R de conferência numérica, com o mesmo número do
resultado. A regra `n ≤ 200`, `J ≤ 3` e forma densa vale para conferência de
**desenho**; conferência de **taxa** precisa de `J` grande (até 12 em
`02-aproximacao-besov.R`, por quadratura nível a nível) e de `n` até 12800
(`05-taxas.R`), e é exceção registrada. A conferência roda **antes**
de a prova ser escrita; o script imprime `OK` ou falha. Os scripts chamam o
`WaveBased` instalado para as bases (`wbasis()`), nunca reimplementam a
avaliação de wavelets.

O molde de prova é o WALL teórico
(`../../wall-manuscript/manuscript/theo/ms_theo_1.tex`, Seções 5, 7, 8):
ler o passo correspondente antes de escrever o do WAFC, e registrar no
arquivo o que foi transposto e o que é novo.
