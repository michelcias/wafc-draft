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
Rscript run_all.R                                                                                   # o estudo inteiro: centenas de horas de CPU, só com o aviso do autor
```

## 3. Onde o trabalho está (resumo de 2026-10-06; o `ESTADO.md` manda)

- **Fechado:** E0, E1 (E1.2 a E1.15, E1.9 no nível da taxa), E2 inteira
  (E2.1 a E2.5j; o rumo em D44 a D46), E3.1 a E3.4, E5a, E5c, E5d a E5g, as
  sondagens E6.1a a E6.1c, e as rodadas de bibliografia L1 a L12; no
  catálogo, a E4.1c (o oráculo em blocos e o `structure`). Decisões
  D1 a D62 na tabela do `ESTADO.md` §2; E4.1 e E4.1b fechadas, o compêndio na pasta `wafc-studies/` (D60).
- **O método (D44, D45):** o block LASSO na forma balanceada, com os níveis
  livres e os pesos do `grpreg`, seguido do limiar `cv1se`; o LASSO
  coordenado fica como opção (D43). `cv.wafc(x, u, y)` é o estimador (D48).
- **Teoria:** o enunciado principal é o Corolário 11 de `08-blocos.tex`
  (Theorem 1 do manuscrito, D47); a taxa lenta sem condição de desenho
  cobre todo `s' > 0` (Corolário 14, E1.13 e E1.14) e, pelo comparador
  truncado, alcança a do Theorem 1 em `s < 1/2` (Corolário 15, E1.15); a
  cota inferior no nível da taxa (`09-cota-inferior.tex`, E1.9) faz o
  Corolário 11 ótimo em `π ≥ 2` para todo `q`; a teoria do LASSO é a de
  E1.5 e E1.6, com a Proposição 4 completada.
- **Código** (`wafc/`, privado, D4): a interface de D48 com os gráficos de
  E3.2; a suíte padrão passa em menos de 60 s (1 202 testes) e a inteira
  com `WAFC_SLOW_TESTS=1` (1 427). A tabela de simulação mais recente é
  `wafc/cache/e25j/e25j-joined.rds` (37 350 linhas, não versionada).
- **Manuscrito vivo em `k = 4`** (`ms_4.tex`, `supp_4.tex`,
  `references_4.bib`): a teoria em blocos, 52 e 80 páginas com a marcação
  (44 e 79 sem o removido); as Seções 5 a 7 são da E5b; o teto de 40 páginas
  fica para o fim (D21).
- **Bibliografia:** 88 entradas verificadas em
  `docs/referencias-verificadas.bib`; os PDFs em `refs/` (§4), com o padrão
  sobrenome e ano.
- **Aplicação (D57, D59):** marylebone com dados do UK-AIR e do ERA5 (a
  base `marylebone.ukair`, E6.1c) como principal e beijing.heat como
  segunda, cada método nos seus termos (D56)
  (`docs/aplicacao-candidatas.md` §9 a §17).
- **O próximo passo é o estudo de simulação (E4)**, com a proposta de
  desenho na pergunta 54 do `ESTADO.md`.
- **Se aparecer um `docs/handoff-*.md`,** é de chat de tarefa que não foi
  integrado (§7 de `instrucoes.md`).

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
