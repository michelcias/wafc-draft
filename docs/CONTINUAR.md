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
│   ├── wafc-draft/       # este repositório (docs, derivations, wafc/ com o código, wafc-studies/ com o compêndio, manuscript)
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
| `grpreg` | **exigido pelo padrão** do `wafc()` e do `cv.wafc()`, `penalty = "block"` (D44, D48); está em `wafc_depends` do `load.R` desde 2026-10-03 | já no comando abaixo, ou `sudo apt install r-cran-grpreg` |
| `sparsegl` | exigido por `wafc(penalty = "sglasso")`; o resto de `wafc/` roda sem ele, e os testes pulam os blocos de grupo se faltar | `install.packages("sparsegl")` |
| `quadprog` | constante de compatibilidade exata em `check/03-desenho-produtos.R` (E1.4) | `install.packages("quadprog")` |
| `VCBART` | exigido por `wafc_competitor("vcbart")` (E2.4); o resto de `wafc/` roda sem ele e os testes pulam o bloco se faltar. Instalado aqui: 1.2.5 | `install.packages("VCBART")` |
| `/usr/bin/time` (pacote `time` do Ubuntu) | a rodada do piloto (E4.3) mede o pico de memória e o processador de cada passo com `time -v` | `sudo apt install time` |
| `/usr/bin/time` também na máquina da E4.4 | ver §2b | `sudo apt install time` |
| `yaml` | só o compêndio (`wafc-studies/`), para ler `config/*.yaml`; 2.3.12 no `renv.lock` | `renv::restore()` dentro de `wafc-studies/` instala tudo, com o `WaveBased` do GitHub em `e494b0e` |
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
cd manuscript && latexmk -pdf ms_5.tex && latexmk -pdf supp_5.tex && latexmk -c && cd -   # versão viva k = 5; 73 e 101 páginas, sem indefinida
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
NC=4 setsid nohup wafc/cache/e61b/run.sh > wafc/cache/e61b/run.out 2>&1 &   # E6.1b, fit + report: ~39 h em 3 a 4 processos, pico de 4,9 GB por processo; o run.sh não é versionado, e o comando por baixo é o 10-sondagem-aplicacao-b.R all fit,report 4 20 0; os dados de marylebone e kelmarsh em wafc/cache/data/ (handoff da E6.1b)
E61B_EXTRA=gam.cv Rscript wafc/scripts/10-sondagem-aplicacao-b.R all extra 2 20 0   # E6.1b, o gam.cv em dobras por bloco: ~1 h
Rscript wafc/scripts/10-sondagem-aplicacao-b.R all report 1 20 0                   # E6.1b, as tabelas de aplicacao-candidatas.md §11 a §14
Rscript derivations/check/09-cota-inferior.R                                       # E1.9: ~90 s
Rscript wafc/scripts/10-sondagem-aplicacao-b.R marylebone.ukair fit 3 20 0         # E6.1c, marylebone com dados do UK-AIR e do ERA5: 5 h 27 min em 3 processos, pico de 4,4 GB; os dados em wafc/cache/data/marylebone.ukair/ (fontes e SHA-256 no script e em aplicacao-candidatas.md §17)
E61B_EXTRA=gam.cv Rscript wafc/scripts/10-sondagem-aplicacao-b.R marylebone.ukair extra 3 20 0   # E6.1c, o gam.cv: ~3 min
Rscript wafc/scripts/10-sondagem-aplicacao-b.R marylebone.ukair report 1 20 0      # E6.1c, as tabelas da §17
```

O compêndio (`wafc-studies/`, E4.1), sempre de dentro da pasta, com o
`renv` dela:

```bash
cd wafc-studies
Rscript -e 'renv::restore()'                                                                         # o ambiente do compêndio
Rscript scripts/01_simulate.R --sizes=250 --reps=2 --workers=8 --out=../wafc/cache/e41/smoke-core   # fumaça, ~4,6 min, 0 falhas
Rscript scripts/02_aggregate.R --sizes=250 --reps=2 --out=../wafc/cache/e41/smoke-core             # as tabelas da fumaça
Rscript ../wafc/cache/e41/check-junction.R 6                                                        # décima junção exata, ~5 min; o script não é versionado
Rscript data-raw/fetch.R --from=../wafc/cache/data                                                  # as fontes da aplicação, conferidas pelo SHA-256 (sem --from, baixa; só com permissão)
Rscript data-raw/prepare.R --check                                                                  # reconstrói data/*.rds e compara com os versionados
Rscript scripts/03_application.R --workers=3                                                        # as duas aplicações: 168 unidades, ~70 h de CPU, ~24 h em 3 processos de até 5 GB; só com o aviso do autor
Rscript run_all.R                                                                                   # o estudo inteiro e as aplicações: centenas de horas de CPU, só com o aviso do autor
```

O piloto do estudo (E4.3), da raiz de `wafc-draft`; o `run-pilot.sh` não é
versionado e o texto está no bloco da E4.3a do `ESTADO.md`:

```bash
bash wafc/cache/e43/run-pilot.sh smoke                                                              # fumaça: ~4 min, 1 trabalhador
setsid nohup bash wafc/cache/e43/run-pilot.sh > wafc/cache/e43/run-pilot.out 2>&1 < /dev/null &    # a rodada: ~3 h em 8 trabalhadores, sozinha na máquina, só com o aviso do autor
Rscript wafc/scripts/11-pilot-read.R                                                                # a leitura; na fumaça, --dir=wafc/cache/e43/smoke --top-n=250 --arms-n=250
```

## 2b. A produção do estudo (E4.4) noutra máquina

A E4.4 roda numa máquina só, porque o tempo de cada unidade é métrica do
estudo (D81). A máquina prevista é a desktop do autor (Intel i9-10900KF,
10 núcleos físicos com hyperthreading, ~50 GB): **10 trabalhadores**, um
por núcleo físico, até ~38 GB no pico (3,8 GB por processo na `mixed` em
`n = 2000`). O passo a passo, nela:

1. **Ferramentas:** R **4.6.1** (a versão do `renv.lock`), compilador C e
   Fortran (`sudo apt install build-essential gfortran`), o `/usr/bin/time`
   (`sudo apt install time`), `git` e acesso ao `michelcias/wafc-draft`
   (privado; `gh auth login` ou um token). **A BLAS tem de ser a de
   referência**, como aqui (não instalar `libopenblas`):
   `Rscript -e 'extSoftVersion()["BLAS"]'` deve mostrar
   `.../blas/libblas.so.3...`.
2. **O repositório e o ambiente:**

   ```bash
   mkdir -p ~/Documentos/repo && cd ~/Documentos/repo
   git clone https://github.com/michelcias/wafc-draft.git
   cd wafc-draft/wafc-studies
   Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
   Rscript -e 'renv::restore()'      # instala os 19 pacotes, o WaveBased do GitHub em e494b0e
   ```

3. **A junção entre máquinas** (minutos; tem de dar `JUNCTION OK` antes da
   produção; as unidades de referência foram ajustadas aqui em 2026-10-08,
   com o código `6aa4a3a`):

   ```bash
   cd ~/Documentos/repo/wafc-draft
   mkdir -p wafc/cache/e44 && tar xzf wafc/scripts/14-junction-ref.tar.gz -C wafc/cache/e44
   cd wafc-studies
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript scripts/01_simulate.R --config=../wafc/scripts/11-pilot.yaml --out=../wafc/cache/e44/junction-other --workers=10 --cells=smooth,mixed,inhomogeneous.xu --sizes=250,1000 --reps=1
   cd .. && Rscript wafc/scripts/14-junction.R wafc/cache/e44/junction-ref wafc/cache/e44/junction-other
   ```

   Se der `JUNCTION FAILED`, parar e trazer a saída: a causa provável é a
   BLAS ou a versão de algum pacote.
4. **A produção** (~592 h de processador; em 10 núcleos mais rápidos que
   os daqui, ~45 a 60 h), da raiz de `wafc-draft`, com a árvore limpa:

   ```bash
   setsid nohup bash wafc/scripts/14-production.sh 10 > wafc-studies/outputs/production.out 2>&1 < /dev/null &
   tail wafc-studies/outputs/production.log      # as marcas; "production end" no fim
   ```

   O script recusa começar com mudança não commitada no código ou na
   configuração, grava a máquina em `outputs/machine.txt` e retoma do
   cache se for parado e relançado.
5. **No meio da rodada** (pergunta 69(a)): contar os candidatos de `k` que
   não convergiram, por célula, nas unidades já gravadas:

   ```bash
   cd ~/Documentos/repo/wafc-draft/wafc-studies
   Rscript -e 'f <- list.files("outputs", "rds$", recursive = TRUE, full.names = TRUE); f <- f[grepl("/gam\\.reml/", f)]; e <- do.call(rbind, lapply(f, function(x) { g <- readRDS(x)$side$gam_k; if (is.null(g)) NULL else data.frame(cell = g$cell, n = g$n, k = g$k, nc = grepl("not converged", g$error)) })); print(aggregate(nc ~ cell + n + k, e, sum))'
   ```

   Fora do `urho` e da `mixed` em `k` grande, um "not converged" é motivo
   para parar e avisar.
6. **Depois da rodada:** as tabelas vão a `results/` na própria desktop
   (com origem e SHA-256 no `results/README.md`). Se for preciso levar as
   saídas para outra máquina, empacotar:

   ```bash
   cd ~/Documentos/repo/wafc-draft/wafc-studies
   tar czf ~/e44-outputs.tar.gz outputs/core outputs/arms outputs/scale outputs/machine.txt outputs/production.log outputs/steps
   ```

   e, no destino, `tar xzf e44-outputs.tar.gz -C wafc-studies/`.

   **Para levar do notebook à desktop** (antes de começar lá; não
   versionado): `tar czf ~/wafc-extra.tar.gz refs wafc-studies/outputs/application`
   na raiz de `wafc-draft` no notebook, e `tar xzf wafc-extra.tar.gz` na
   raiz de `wafc-draft` na desktop.

## 3. Onde o trabalho está (resumo de 2026-10-08; o `ESTADO.md` manda)

- **Fechado:** E0 a E3 inteiras, E4.1 a E4.3b, E5a, E5c a E5k, E6.1a a
  E6.2b, L1 a L14. Decisões D1 a D81 na tabela do `ESTADO.md` §2.
  Catálogo vazio.
- **O método (D44, D45):** o block LASSO na forma balanceada seguido do
  limiar `cv1se`; `cv.wafc(x, u, y)` é o estimador (D48).
- **Código** (`wafc/`, D4): a suíte padrão passa em ~1 min (1 770 testes;
  2 003 com `WAFC_SLOW_TESTS=1`); o último commit de `wafc/R` é
  `6aa4a3a` (árvore `34f2e04`).
- **O compêndio** (`wafc-studies/`): o estudo (`core`, `arms`, `scale`) e
  as duas aplicações; as aplicações rodaram (168 unidades) e as tabelas e
  figuras estão em `results/`; o estudo espera a E4.4.
- **Manuscrito vivo em `k = 5`**: 73 e 101 páginas com a marcação (66 e
  100 sem o removido); a Seção 6 com os números finais, a 5 à espera da
  E4.4, sem a Seção 7.
- **Bibliografia:** 101 entradas verificadas em
  `docs/referencias-verificadas.bib`.
- **O próximo passo é a E4.4 na desktop do autor** (§2b; D81).
- **Se aparecer um `docs/handoff-*.md`,** é de chat de tarefa que não foi
  integrado (§7 de `instrucoes.md`).

## 4. O que não viaja

- A memória local do assistente (`~/.claude/...`) e o scratchpad: nada do
  projeto depende deles; tudo que importa está nos repositórios.
- Resultados de piloto (`.rds`) e `wafc/cache/` (1,1 GB): não
  versionados; só o `ESTADO.md` registra os números. Ficam na máquina onde
  rodaram.
- **As unidades e as saídas do compêndio** (`wafc-studies/outputs/`): não
  versionadas. As da aplicação (3,5 MB) estão no notebook; para refazer o
  relatório das figuras noutra máquina sem reajustar, copiar a pasta
  `wafc-studies/outputs/application/`. As da E4.4 nascem na desktop e ficam
  lá; o que vai ao repositório são as tabelas copiadas para `results/`.
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
