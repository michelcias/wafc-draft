# wafc-draft

Repositório de trabalho do manuscrito sobre **regressão com coeficientes
funcionais aditivos estimados por LASSO em bases de wavelets** (WAFC, nome
provisório). O modelo é

```
Y = Σ_j β_j(U) X_j + ε,      β_j(u) = c_j + Σ_k g_{jk}(u_k),
```

cada `g_{jk}` expandida numa base ortonormal de wavelets de suporte compacto
e o vetor de coeficientes estimado por LASSO. É a generalização, para
regressão com coeficientes funcionais, do que o WALL (`wall()` no pacote
`WaveBased`) faz para o modelo aditivo logístico.

Alvo de publicação: periódico Q1 em Statistics and Probability (SCImago);
proposta **Statistica Sinica**, reserva *Electronic Journal of Statistics*
([`docs/alvo-revista.md`](docs/alvo-revista.md)).

> **Comece por [`CLAUDE.md`](CLAUDE.md)** (as regras que valem desde a primeira
> mensagem) **e depois [`docs/ESTADO.md`](docs/ESTADO.md)** (onde o trabalho
> parou, o que foi decidido e o que vem a seguir).

---

## Mapa do repositório

```
wafc-draft/
├── CLAUDE.md                    # regras da sessão (leia primeiro)
├── README.md
├── docs/                        # documentos de trabalho
│   ├── ESTADO.md                #   handoff: estado atual e próximos passos
│   ├── instrucoes.md            #   normativo: fluxo, escrita, marcação, git, chats paralelos
│   ├── TAREFA.md                #   abertura de um chat de tarefa: o que ler, catálogo, regras
│   ├── CONTINUAR.md             #   retomar de outra máquina: pastas, ferramentas, o que não viaja
│   ├── MENSAGENS.md             #   cola do autor: como abrir conversa principal, de tarefa e em outra máquina
│   ├── plano-projeto.md         #   etapas E0 a E7 (+ L1, L2), dependências, entregáveis, critérios de saída
│   ├── proposta-metodo.md       #   a ideia inicial: modelo, estimador, teoria esperada, riscos
│   ├── alvo-revista.md          #   escolha da revista e o que ela pede
│   ├── literatura.md            #   leituras, o que cada uma resolve, status de verificação
│   ├── inventario-codigo.md     #   o que o WaveBased, o wall e o wall-manuscript têm de reaproveitável
│   ├── notacao.md               #   símbolos (esboço; congela em E1.1)
│   └── referencias-verificadas.bib
├── manuscript/                  # ms_{k}.tex, supp_{k}.tex, references_{k}.bib (nascem em E5a)
│   └── ejs-template/            #   template da EJS (imsart.cls, opção ejsv2), compila
├── derivations/                 # rascunhos matemáticos, um resultado por arquivo
│   └── check/                   #   conferência numérica em R (n pequeno, forma densa)
├── prototype/                   # protótipo do estimador em R sobre o WaveBased (E2)
└── results/                     # cópia de referência de figuras e tabelas (origem: wafc-studies)
```

## Convenção de versionamento do manuscrito

Os três arquivos da versão viva andam juntos com o mesmo índice `k`:
`manuscript/ms_{k}.tex`, `manuscript/supp_{k}.tex`,
`manuscript/references_{k}.bib`. Antes de qualquer alteração pergunta-se se a
versão `{k+1}` deve ser criada. Regra completa em
[`docs/instrucoes.md`](docs/instrucoes.md).

## Convenção de marcação

Nenhuma edição em `.tex` é silenciosa a partir de `k = 2`.

| Finalidade | Comando | Cor |
|---|---|---|
| Remoção sugerida | `\textcolor{gray}{\sout{...}}` | cinza tachado |
| Inserção da rodada corrente | `\textcolor{colR1}{...}` | azul |

O índice da cor (`colR1`, `colR2`, ...) avança a cada rodada aceita; a rodada
corrente é registrada em `docs/ESTADO.md`.

## Git

Commit e push só quando pedidos, **e nunca com coautoria**: sem
`Co-Authored-By:`, sem "Generated with", sem trailer que nomeie o assistente.
Regra completa em [`docs/instrucoes.md`](docs/instrucoes.md).

## Repositórios relacionados

| Repositório | Papel |
|---|---|
| [`WaveBased`](https://github.com/michelcias/WaveBased) | pacote R/C do autor: bases de wavelets (Daubechies–Lagarias, tabelas, CDV), `wall()`; recebe `wafc()` em E3 |
| [`wall-manuscript`](https://github.com/michelcias/wall-manuscript) | artigos do WALL; o teórico é o molde da prova |
| `wall` (local) | compêndio de benchmark do WALL; molde do `wafc-studies` |
| [`bdm-draft`](https://github.com/michelcias/bdm-draft) | projeto irmão, origem destas convenções |
| `wafc-studies` | compêndio de simulação e aplicação; a criar em E4.1 |
