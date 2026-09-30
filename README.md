# wafc-draft

Repositório de trabalho do manuscrito sobre **regressão com coeficientes
funcionais aditivos estimados por LASSO em bases de wavelets** (WAFC,
*wavelet additive functional coefficients*, D8). O modelo é

```
Y = Σ_ℓ β_ℓ(U) X_ℓ + ε,      β_ℓ(u) = c_ℓ + Σ_m g_{ℓm}(u_m),      ∫_0^1 g_{ℓm} = 0,
```

cada `g_{ℓm}` expandida numa base ortonormal de wavelets `ψ_{jk}` de suporte
compacto e o vetor de coeficientes `θ_{ℓm,jk}` estimado por LASSO (notação
de [`docs/notacao.md`](docs/notacao.md)). É a generalização, para
regressão com coeficientes funcionais, do que o WALL (`wall()` no pacote
`WaveBased`) faz para o modelo aditivo logístico. O código do método vive
na pasta [`wafc/`](wafc/) deste repositório.

Alvo de publicação: periódico Q1 em Statistics and Probability (SCImago);
alvo **Statistica Sinica** (D5), reserva *Electronic Journal of Statistics*
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
│   ├── ss-instrucoes-autores.md #   instruções da Statistica Sinica, transcritas
│   ├── literatura.md            #   leituras, o que cada uma resolve, status de verificação
│   ├── busca-novidade.md        #   L2: busca de novidade e veredito por contribuição
│   ├── selecao-estrutura.md     #   as três saídas para a seleção de estrutura (D28, D32)
│   ├── aplicacao-candidatas.md  #   E6.1a: sondagem das bases de aplicação, com veredito
│   ├── inventario-codigo.md     #   o que o WaveBased, o wall e o wall-manuscript têm de reaproveitável
│   ├── notacao.md               #   símbolos (congelada em E1.1)
│   └── referencias-verificadas.bib
├── manuscript/                  # ms_{k}.tex, supp_{k}.tex, references_{k}.bib (vivo em k = 1)
│   ├── ss-template/             #   template da Statistica Sinica, compila
│   └── ejs-template/            #   template da EJS (imsart.cls, opção ejsv2), compila
├── derivations/                 # rascunhos matemáticos, um resultado por arquivo
│   └── check/                   #   conferência numérica em R (n pequeno, forma densa)
├── wafc/                        # o código do método em R (E2): R/, tests/, scripts/
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
| [`WaveBased`](https://github.com/michelcias/WaveBased) | pacote R/C do autor: bases de wavelets (Daubechies–Lagarias, tabelas, CDV); dependência do código em `wafc/`; `wall()` é referência de leitura |
| [`wall-manuscript`](https://github.com/michelcias/wall-manuscript) | artigos do WALL; o teórico é o molde da prova |
| `wall` (local) | compêndio de benchmark do WALL; molde do `wafc-studies` |
| [`bdm-draft`](https://github.com/michelcias/bdm-draft) | projeto irmão, origem destas convenções |
| `wafc-studies` | compêndio de simulação e aplicação; a criar em E4.1 |
