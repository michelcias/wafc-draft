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
| `sparsegl` | exigido por `wafc(penalty = "sglasso")`; o resto de `wafc/` roda sem ele, e os testes pulam os blocos de grupo se faltar | `install.packages("sparsegl")` |
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
cd manuscript && latexmk -pdf ms_1.tex && latexmk -pdf supp_1.tex && latexmk -c && cd -   # 28 e 26 páginas, sem indefinida
```

As conferências das derivações, que devem imprimir `OK` (tempos desta
máquina):

```bash
Rscript derivations/check/01-identificabilidade.R   # ~7 s
Rscript derivations/check/03-desenho-produtos.R     # ~9 s, precisa de quadprog
Rscript derivations/check/04-oraculo.R              # ~16 s, precisa de glmnet
Rscript derivations/check/02-aproximacao-besov.R    # ~72 s (com a emenda E1.3b)
Rscript derivations/check/05-taxas.R                # ~7 min
```

E o código do método, que já existe:

```bash
Rscript -e 'testthat::test_dir("wafc/tests")'   # 420 passam, ~20 s
Rscript wafc/scripts/01-smoke.R                 # imprime OK, ~3 s
Rscript wafc/scripts/03-tune-decomp.R           # ~2 min
Rscript wafc/scripts/02-tune.R 20               # a comparação de E2.3, ~31 min
```

## 3. Onde o trabalho está (resumo de 2026-09-19; o `ESTADO.md` manda)

- **E0 fechada.** Código em `wafc/` (D4); alvo *Statistica Sinica* (D5,
  confirmada Q1 no SCImago) e método chamado **WAFC** (D8); instruções da
  revista em `docs/ss-instrucoes-autores.md` e templates em
  `manuscript/ss-template/`; teto de 40 páginas em espaço duplo, referências
  incluídas. O tipo de revisão não consta da página oficial.
- **Notação congelada** (E1.1, D9 a D12): `ψ_{jk}`, `X_ℓ`, `U_m`,
  `θ_{ℓm,jk}`, `U ∈ [0,1]^q`; `derivations/macros.tex` implementa e compila.
- **Teoria: E1 fechada** (E1.2 a E1.6, mais a emenda E1.3b), com as cinco
  conferências numéricas rodando. O enunciado principal é o Corolário 5 de
  `05-taxas.tex` (D16). E1.7 espera a decisão de L2f
  (`docs/selecao-estrutura.md`); E1.8 (rota do intervalo) está aberta e E1.9
  (cota inferior) fica para depois dos resultados.
- **Código: E2.1 a E2.3 fechadas**, mais a emenda E2.1b (margem `eps` fixa,
  desacoplada de `J`). `wafc/R/` tem o carregador, os cenários,
  `wafc_design()`, `wafc()` com as duas penalidades e a sintonia
  (`cv.wafc()`, BIC, EBIC, regra da teoria); **420 testes passam**. E2.4
  (piloto) está aberta; falta o go/no-go (E2.5).
- **Avaliação da base em estudo numérico:** tabela fixa (`wtable()` uma vez,
  passada em `wavelet.table`), nunca a regra `auto` réplica a réplica (D31);
  a exceção são as conferências de `derivations/check/`, que medem precisão
  fina.
- **Manuscrito vivo em `k = 1`:** `manuscript/ms_1.tex` e `supp_1.tex`
  compilam limpos, em 28 e 26 páginas. Alterar pede decidir antes se nasce
  `k = 2`, e há uma lista de edições acumuladas esperando (pergunta 9 do
  `ESTADO.md`).
- **Bibliografia:** 63 entradas verificadas (L1, L3 e L4); nada pendente.
- **Abertas no catálogo:** E1.8 (rota do intervalo, em `derivations/`) e
  E2.4 (piloto, em `wafc/`); não compartilham arquivo.
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
