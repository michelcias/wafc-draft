# manuscript

Manuscrito do WAFC. **Ainda não existe** (nasce em E5a como `k = 1`):
`ms_1.tex`, `supp_1.tex`, `references_1.bib`. Os três andam juntos com o
mesmo índice `k`; a versão `{k+1}` copia os três e atualiza
`\bibliography{references_{k+1}}` nos dois `.tex` (`docs/instrucoes.md`, §3).
A versão `1` não leva marcação de alteração; o preâmbulo já define `colR1`
e carrega `ulem` para quando `k = 2` existir.

| Pasta | O que é |
|---|---|
| `ejs-template/` | template da *Electronic Journal of Statistics* (`imsart.cls` 2025/03/18, opção `ejsv2`, `imsart.sty`, `imsart-nameyear.bst`, `ejs-sample.tex/.pdf`, LPPL), clonado de `vtex-soft/texsupport.ims_cosponsored-ejs` em 2026-09-18; o `ejs-sample.tex` compila com o `latexmk` local |
| `ss-template/` | template da *Statistica Sinica*; **a baixar pelo autor em E0.3** (o site não resolveu daqui) |

O alvo primário proposto é a *Statistica Sinica* (`docs/alvo-revista.md`).
Enquanto o template dela não estiver aqui, o `ms_1.tex` pode nascer em
`article` simples com as provas em `supp_1.tex`, e a mudança de classe é
uma linha.

## Compilar (quando existir)

```bash
cd manuscript && TEXINPUTS=./ejs-template//: BSTINPUTS=./ejs-template//: latexmk -pdf ms_1.tex
```
