# derivations

Um resultado por arquivo, numerado na ordem de dependência do plano (E1):
`01-identificabilidade.md`, `02-aproximacao-besov.tex`,
`03-desenho-produtos.tex`, `04-oraculo.tex`, `05-taxas.tex`, ... Cada um com
enunciado, hipóteses, prova, "o que isso não cobre" e o parágrafo da
conferência numérica. Documentos `article` autônomos que fazem
`\input{macros}`; `macros.tex` implementa `docs/notacao.md` e é o único
lugar onde uma macro é definida. Compilar com `latexmk -pdf <arquivo>.tex`;
o `.pdf` é versionado, os auxiliares não.

`check/` guarda os scripts R de conferência numérica (`n ≤ 200`, `J ≤ 3`,
forma densa), com o mesmo número do resultado. A conferência roda **antes**
de a prova ser escrita; o script imprime `OK` ou falha. Os scripts chamam o
`WaveBased` instalado para as bases (`wbasis()`), nunca reimplementam a
avaliação de wavelets.

O molde de prova é o WALL teórico
(`../../wall-manuscript/manuscript/theo/ms_theo_1.tex`, Seções 5, 7, 8):
ler o passo correspondente antes de escrever o do WAFC, e registrar no
arquivo o que foi transposto e o que é novo.
