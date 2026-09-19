# Continuar de outra máquina

Arquivo de abertura para retomar o projeto num computador novo. A mensagem
de abertura é uma linha:

> Leia `docs/CONTINUAR.md` e me diga o que falta instalar; depois leia
> `docs/ESTADO.md` e continuamos de onde parou.

O chat que ler isto deve: (1) conferir a §2 e dizer o que falta; (2) ler
`CLAUDE.md`, `docs/ESTADO.md`, `docs/instrucoes.md`; (3) esperar o autor
dizer o foco do dia. Não commitar nada por iniciativa própria.

---

## 1. Repositórios e disposição das pastas

Os documentos usam caminhos relativos (`CLAUDE.md`, tabela de vizinhos),
então a disposição tem de ser esta (o nome da pasta-mãe varia entre
máquinas: `~/Documents` nesta, `~/Documentos` em outras; os caminhos
relativos não dependem disso):

```
~/Documents/
├── repo/
│   ├── wafc-draft/       # este repositório (docs, derivations, wafc/ com o código, manuscript)
│   └── bdm-draft/        # projeto irmão, só consulta (origem das convenções)
├── WaveBased/            # o pacote (michelcias/WaveBased); só dependência para as bases (D4)
├── wall-manuscript/      # artigos do WALL (michelcias/wall-manuscript), só leitura; molde da prova
└── wall/                 # compêndio do WALL, só leitura; molde do wafc-studies
```

```bash
mkdir -p ~/Documents/repo && cd ~/Documents/repo
git clone https://github.com/michelcias/wafc-draft.git
cd ~/Documents
git clone https://github.com/michelcias/WaveBased.git
git clone https://github.com/michelcias/wall-manuscript.git
```

O `wall` (compêndio) não tem remoto registrado nesta máquina; se precisar
dele em outra, copiar ou publicar. Só E4.1 depende dele, como molde.

O repositório `michelcias/wafc-draft` é **privado** (publicado em
2026-09-18), e o código em `wafc/` fica privado por enquanto (D4): nada de
`install_github`, pacote público ou cópia para outro repositório sem
decisão em E3.3.

## 2. Ferramentas

| O quê | Para quê | Conferir |
|---|---|---|
| R ≥ 4.1 (aqui 4.6.1) com compilador C | `WaveBased` compila de fonte | `R --version` |
| `latexmk`, `pdflatex`, `bibtex` | `derivations/*.tex`, templates, manuscrito | `latexmk --version` |
| `gh` autenticado como `michelcias` (opcional) | criar repositórios, CI | `gh auth status` |
| pacotes R: `glmnet`, `Matrix`, `mgcv`, `grpreg`, `gglasso`, `bench`, `testthat`, `devtools`, `roxygen2`, `remotes`, `renv` | protótipo, competidores, pacote | comando abaixo |
| `sparsegl` | sparse group LASSO (E2.2), se a variante for adotada | `install.packages("sparsegl")` |
| `quadprog` | constante de compatibilidade exata em `check/03-desenho-produtos.R` (E1.4) | `install.packages("quadprog")` |
| o próprio `WaveBased`, instalado de `../../WaveBased` (ou `remotes::install_github("michelcias/WaveBased")`) | bases de wavelets (`wbasis()`, `wtable()`) chamadas por `wafc/R/design.R`; não recebe código | `cd ~/Documents/WaveBased && R CMD INSTALL .` |

No Ubuntu, `grpreg`, `gglasso` e `sparsegl` também saem do repositório da
distribuição (`sudo apt install r-cran-grpreg r-cran-gglasso r-cran-sparsegl`),
que instala em `/usr/lib/R/site-library`; o `R` acha do mesmo jeito.

```r
pk <- c("glmnet", "Matrix", "mgcv", "grpreg", "gglasso", "bench", "testthat", "devtools", "roxygen2", "remotes", "renv")
install.packages(setdiff(pk, rownames(installed.packages())), repos = "https://cloud.r-project.org")
```

Conferência de que tudo roda (da raiz de `wafc-draft`):

```bash
Rscript -e 'library(WaveBased); w <- wbasis(sort(runif(64)), j0 = 0, J = 3); cat(dim(w), "\n")'   # 64 8
cd manuscript/ejs-template && latexmk -pdf ejs-sample.tex && latexmk -c && cd -   # compila
cd manuscript/ss-template && latexmk -pdf SS-template.tex && latexmk -c && cd -   # compila (9 páginas)
```

As conferências das derivações, que devem imprimir `OK` (tempos desta
máquina):

```bash
Rscript derivations/check/01-identificabilidade.R   # ~7 s
Rscript derivations/check/03-desenho-produtos.R     # ~9 s, precisa de quadprog
Rscript derivations/check/04-oraculo.R              # ~16 s, precisa de glmnet
Rscript derivations/check/02-aproximacao-besov.R    # ~53 s
```

E o código do método, que já existe:

```bash
Rscript -e 'testthat::test_dir("wafc/tests")'   # 103 passam, ~5 s
Rscript wafc/scripts/01-smoke.R                 # imprime OK, ~3 s
```

## 3. Onde o trabalho está (resumo de 2026-09-18; o `ESTADO.md` manda)

- **E0 fechada, menos a ratificação de D5:** D4 decidida (código em `wafc/`);
  D5 e D8 adiadas; repositório publicado; E0.3 feita (instruções da
  *Statistica Sinica* em `docs/ss-instrucoes-autores.md`, templates em
  `manuscript/ss-template/`, teto corrigido para 40 páginas em espaço duplo).
  O tipo de revisão da revista não consta da página oficial.
- **Notação congelada** (E1.1, D9 a D12): `ψ_{jk}`, `X_ℓ`, `U_m`,
  `θ_{ℓm,jk}`, `U ∈ [0,1]^q`; `derivations/macros.tex` implementa e compila.
- **Teoria: E1.2 a E1.5 fechadas**, com conferência numérica rodando; de E1
  só falta E1.6 (taxas e compressibilidade).
- **Código: E2.1 fechada.** `wafc/R/` tem o carregador, os cenários e
  `wafc_design()`; 103 testes passam. Falta `wafc()` (E2.2). **Nada de
  manuscrito ainda.**
- **Tarefas que podem abrir agora em chats de tarefa:** E1.6 e E2.2
  (catálogo em `TAREFA.md`).
- **Nenhum handoff pendente.** Se aparecer um `docs/handoff-*.md`, é de chat
  de tarefa que não foi integrado: o protocolo está na §7 de
  `instrucoes.md`.

## 4. O que não viaja

- A memória local do assistente (`~/.claude/...`) e o scratchpad: nada do
  projeto depende deles; tudo que importa está nos repositórios.
- Resultados de piloto (`.rds`) e `wafc/cache/`: não versionados; só o
  handoff e o `ESTADO.md` registram os números.
- A instalação do `WaveBased`: reinstalar com `R CMD INSTALL .` na pasta
  clonada.

## 5. Regras que valem lá como aqui

`CLAUDE.md`: sem coautoria em commit; respostas curtas. `instrucoes.md`:
commit e push são do autor; chats de tarefa começam por "Leia
`docs/TAREFA.md` e execute a tarefa X" e terminam em handoff; o chat
principal integra. **§8 do `instrucoes.md`**: este arquivo, o `ESTADO.md` e o
`TAREFA.md` são o contrato de continuidade, e mantê-los em dia é parte de
fechar uma etapa.

As mensagens de abertura estão em [`MENSAGENS.md`](MENSAGENS.md).
