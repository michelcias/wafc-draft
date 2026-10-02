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
| `VCBART` | exigido por `wafc_competitor("vcbart")` (E2.4); o resto de `wafc/` roda sem ele e os testes pulam o bloco se faltar. Instalado aqui: 1.2.5 | `install.packages("VCBART")` |
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
cd manuscript && latexmk -pdf ms_3.tex && latexmk -pdf supp_3.tex && latexmk -c && cd -   # versão viva k = 3; 34 e 30 páginas, sem indefinida
```

As conferências das derivações, que devem imprimir `OK` (tempos desta
máquina):

```bash
Rscript derivations/check/01-identificabilidade.R   # ~7 s
Rscript derivations/check/03-desenho-produtos.R     # ~9 s, precisa de quadprog
Rscript derivations/check/04-oraculo.R              # ~16 s, precisa de glmnet
Rscript derivations/check/06a-irrepresentabilidade.R # ~4 min
Rscript derivations/check/06-selecao-limiar.R       # ~7 min, precisa de glmnet e sparsegl
Rscript derivations/check/02-aproximacao-besov.R    # ~72 s (com a emenda E1.3b)
Rscript derivations/check/05-taxas.R                # ~7 min
Rscript derivations/check/07-rota-intervalo.R      # ~2 min
Rscript derivations/check/08a-blocos.R             # ~1 min 30 s, a sondagem de E1.11
```

E o código do método, que já existe:

```bash
Rscript -e 'testthat::test_dir("wafc/tests")'   # 1 006 passam, ~70 s
Rscript wafc/scripts/01-smoke.R                 # imprime OK, ~3 s
Rscript wafc/scripts/03-tune-decomp.R           # ~2 min
Rscript wafc/scripts/02-tune.R 20               # a comparação de E2.3, ~31 min com a grade antiga (hoje o padrão é 2:8, D34)
Rscript wafc/scripts/04-pilot.R 5 competitors 8 250 smooth   # fumaça do piloto, ~1 min
Rscript wafc/scripts/04-pilot.R 50 all 12                    # o piloto inteiro de E2.4, ~3 h em 12 núcleos com a grade antiga
WAFC_OUT=wafc/cache/e25a WAFC_TAG=e25a Rscript wafc/scripts/04-pilot.R 50 competitors 12 "" smooth,inhomogeneous,null,uneven   # E2.5a, parte 1: ~4 h em 12 núcleos
# variáveis do 04-pilot.R: WAFC_OUT (pasta), WAFC_TAG (prefixo dos .rds), WAFC_METHODS (só estes métodos, sem mudar dados nem sementes), WAFC_REPS_MIXED (réplicas da mixed); máquina com 8 núcleos físicos
WAFC_OUT=wafc/cache/e25a WAFC_TAG=e25a-mixed Rscript wafc/scripts/04-pilot.R 50 competitors 4 "" mixed                       # E2.5a, parte 2: ~1 h em 4 núcleos, pico de 1,35 GB por processo
Rscript wafc/scripts/05-sondagem-aplicacao.R all           # a sondagem de E6.1a, horas; baixa os dados
Rscript wafc/scripts/06-timing.R                             # só tempo, curto
Rscript wafc/scripts/07-accel.R 20 staged,sgl 14             # pergunta 29: ~10 min + ~48 min em 14 núcleos
Rscript wafc/scripts/08-sgl-null.R curve,cross,cost 14       # pergunta 31: ~25 min em 14 núcleos
WAFC_OUT=wafc/cache/e25d Rscript wafc/scripts/09-gam-autonomo.R 50 fit,report 12   # E2.5d: 1 h 18 min em 12 núcleos, pico de 1,3 GB por processo, ~11 GB no total
WAFC_OUT=wafc/cache/e25g WAFC_TAG=e25g WAFC_METHODS=wafc.lasso,klopp.balanced WAFC_REPS_MIXED=50 Rscript wafc/scripts/04-pilot.R 50 competitors 8   # E2.5g, os limiares: 1 h 59 min em 8 núcleos, pico de 0,86 GB por processo
NC=8 setsid nohup wafc/cache/e25h/run.sh > wafc/cache/e25h/run.out 2>&1 &   # E2.5h, o gam por REML, CV e GCV e o GCV do WAFC (D41): 10 h 08 min em 8 núcleos, pico de 1,2 GB; o run.sh não é versionado, e o comando por baixo é o 04-pilot.R com WAFC_METHODS=gam.reml,gam.cv,gam.gcv,wafc.gcv (gam.gcv fora da mixed)
```

## 3. Onde o trabalho está (resumo de 2026-10-02; o `ESTADO.md` manda)

- **Fechado:** E0, E1 (E1.2 a E1.8, E1.10, E1.11), E2.1 a E2.5h, E5a, E5c, a
  sondagem E6.1a e as rodadas de bibliografia L1 a L9. Decisões D1 a D42 na
  tabela do `ESTADO.md` §2.
- **Teoria:** o enunciado principal é o Corolário 5 de `05-taxas.tex`
  (Theorem 1 do manuscrito, D16); a seleção de estrutura é o Corolário 8 de
  `06-selecao-limiar.tex` (Corollary 2, D32), agora com a janela a `n` fixo
  e a versão no regime do Theorem 1; a rota do intervalo (E1.8) está
  escrita e não vai ao manuscrito (D26); a teoria do block LASSO é a
  sondagem `08a` (E1.11, E2.5e, E2.5f). O mergulho de Besov em weak-`ℓ_τ`
  (Lema 9(i)) é clássico (Johnstone 2019, Proposição 9.15) e está citado
  assim.
- **Código** (`wafc/`, privado, D4): `wafc()`, `cv.wafc()`, os
  concorrentes (`gam` com `k` por REML, CV ou GCV na grade de D41; `bsgl`;
  as cinco formas do `klopp`; `aspline`; `vcbart`; `linear`; `oracle`), o
  limiar `wafc_threshold()` (E2.5g) e o GCV `wafc_gcv()` (E2.5h); **1 006
  testes passam**. A tabela mais recente é
  `wafc/cache/e25h/e25h-joined.rds` (28 350 linhas, não versionada); o
  procedimento para acrescentar um método por `WAFC_METHODS` foi validado
  seis vezes (`TAREFA.md` §3).
- **Avaliação da base em estudo numérico:** tabela fixa (`wtable()` uma vez,
  passada em `wavelet.table`), que é o caminho rápido; a regra `auto` não
  dispara nos `n` do estudo (D31). A exceção são as conferências de
  `derivations/check/`, que medem precisão fina.
- **Manuscrito vivo em `k = 3`** (`ms_3.tex`, `supp_3.tex`,
  `references_3.bib`, aberto em 2026-10-01 para a notação de D40): compilam
  limpos em 34 e 30 páginas; as edições desde a `k = 2` estão marcadas em
  `colR1` (a rodada não foi aceita); as Seções 5 a 7 são de E5b.
- **Bibliografia:** 83 entradas verificadas em
  `docs/referencias-verificadas.bib`, sem marca aberta; as citações das
  derivações conferidas nas fontes publicadas (L6, L7, L9); os PDFs ficam em
  `refs/` (§4).
- **Aplicação:** a sondagem E6.1a deu veredito negativo nas três bases
  (`docs/aplicacao-candidatas.md`); a saída é decisão do autor (pergunta 2).
- **O próximo passo é a decisão de rumo** (perguntas 33, 38 e 41), e a
  §5 do `ESTADO.md` diz com quais números. **Catálogo vazio** e **nenhum
  handoff pendente**; se aparecer um `docs/handoff-*.md`, é de chat de
  tarefa que não foi integrado (§7 de `instrucoes.md`).

## 4. O que não viaja

- A memória local do assistente (`~/.claude/...`) e o scratchpad: nada do
  projeto depende deles; tudo que importa está nos repositórios.
- Resultados de piloto (`.rds`) e `wafc/cache/`: não versionados; só o
  handoff e o `ESTADO.md` registram os números.
- A instalação do `WaveBased`: reinstalar com `R CMD INSTALL .` na pasta
  clonada.
- **Os PDFs de `refs/`** (não versionados, por direito autoral): as fontes
  que as rodadas de bibliografia leram, com os nomes usados nos comentários
  do `.bib` (`buhlmann2011`, `cai1999`, `donoho1994`, `donoho1998`,
  `hardle1998` e `hardle1998-springer-fm`, `hsu2012`, `huang2010`,
  `johnstone2019`, `kauermann2011`, `klopp2014-arxiv`, `klopp2015` e
  `klopp2015-supp`, `lounici2011`, `mallat2009`, `marra2011`, `pya2016`,
  `restrepo1997`, `ruppert2000`, `ruppert2002`, `tibshirani2012`,
  `vandegeer2011`, `wood2017`, `zou2007`). Em outra máquina, copiar a pasta
  à mão; só uma conferência nova de citação precisa deles.

## 5. Regras que valem lá como aqui

`CLAUDE.md`: sem coautoria em commit; respostas curtas. `instrucoes.md`:
commit e push são do autor; chats de tarefa começam por "Leia
`docs/TAREFA.md` e execute a tarefa X" e terminam em handoff; o chat
principal integra. **§8 do `instrucoes.md`**: este arquivo, o `ESTADO.md` e o
`TAREFA.md` são o contrato de continuidade, e mantê-los em dia é parte de
fechar uma etapa.

As mensagens de abertura estão em [`MENSAGENS.md`](MENSAGENS.md).
