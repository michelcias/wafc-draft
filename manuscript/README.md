# manuscript

Manuscrito do WAFC. **Existe desde 2026-09-19, em `k = 1`** (E5a):
`ms_1.tex` (Seções 1 a 4, 28 páginas), `supp_1.tex` (as provas, 26 páginas)
e `references_1.bib` (30 das 63 entradas verificadas). Os dois compilam com
`latexmk -pdf` sem referência ou citação indefinida. Os três andam juntos com o
mesmo índice `k`; a versão `{k+1}` copia os três e atualiza
`\bibliography{references_{k+1}}` nos dois `.tex` (`docs/instrucoes.md`, §3).
A versão `1` não leva marcação de alteração; o preâmbulo já define `colR1`
e carrega `ulem` para quando `k = 2` existir.

| Pasta | O que é |
|---|---|
| `ejs-template/` | template da *Electronic Journal of Statistics* (`imsart.cls` 2025/03/18, opção `ejsv2`, `imsart.sty`, `imsart-nameyear.bst`, `ejs-sample.tex/.pdf`, LPPL), clonado de `vtex-soft/texsupport.ims_cosponsored-ejs` em 2026-09-18; o `ejs-sample.tex` compila com o `latexmk` local |
| `ss-template/` | template da *Statistica Sinica* (D5), baixado pelo autor em 2026-09-18: `SS-template.tex`, `SS-template-bib.tex` (o que o `ms_1.tex` usa), `supp-temp_20240820.tex` e `carat.pdf`; os três `.tex` compilam localmente |

O alvo é a *Statistica Sinica* (D5, confirmada Q1 no SCImago), e o `ms_1.tex`
já nasceu no template dela. As instruções transcritas estão em
[`../docs/ss-instrucoes-autores.md`](../docs/ss-instrucoes-autores.md): teto
de 40 páginas em espaço duplo **com as referências**, suplementar em PDF
único de até 10 MB, e template conferido na recepção.

As macros estão **copiadas para dentro** dos dois `.tex`, e não em
`\input`, porque o pacote enviado à revista tem de ser autocontido e o
template já define os ambientes de teorema. Substituí-las pelas cópias
geradas dos arquivos de origem é item do checklist de submissão
(`../docs/alvo-revista.md` §6).

## Compilar

```bash
cd manuscript && latexmk -pdf ms_1.tex && latexmk -pdf supp_1.tex && latexmk -c
```
