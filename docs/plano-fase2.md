# Plano da fase 2: do parecer à submissão

**Data:** 2026-11-16. **Estado:** ratificado pelo autor em 2026-11-16
(perguntas 70 a 77, §8; decisões D83 a D90 em `docs/estado/decisoes.md`). O
andamento de cada tarefa está no quadro da cabeça do `ESTADO.md`, não aqui.
**Origem:** `docs/parecer-wafc.html` (commit `cb58595` do branch
`origin/ccr-53b4afa8-g7xu8z`, ainda fora do `main`), lido contra
`main @ 67e45f5`.
**Alcance:** as etapas E0 a E6 do `plano-projeto.md` ficam como registro da
fase 1. Daqui em diante a ordem de execução é a deste arquivo (E8 a E12), e
a E7 (submissão) passa a depender da E12.
**Restrição de execução:** orquestrador e executores são no máximo Opus.
Cada tarefa traz modelo e esforço (§4), e o plano é escrito para que o
julgamento difícil esteja nos critérios declarados aqui, não na hora da
execução.

---

## 1. Leitura do parecer

### 1.1 O que procede

Conferido nesta rodada contra o repositório, não só lido:

| Objeção | Fato conferido | Resposta no plano |
|---|---|---|
| M1, novidade teórica estreita (regime mais fácil que o de K&P) | julgamento do parecer; as frases "weaker conditions" e "no cone because profiled" estão no `ms_6` | E10.3 (reposicionar), E10.5 (constantes em `p`, `q`); Rota T fora (§3) |
| M2, nenhum teorema cobre o procedimento usado | `J` médio 6,9 na `mixed` em `n = 1000`: 2 032 colunas penalizadas | E10.3 (remark da taxa lenta, §1.2(d)); L15 |
| M3, estrutura confundida com o limiar | `lasso.cv1se` marca 1,00 nas cinco células em `n = 1000`, igual ao `wafc.cv1se` (`results/tables/core/summary.csv`) | **G1** |
| M4, ganho pequeno, numa faixa de `n` | razões WAFC/`gam.reml` reproduzidas: 0,993, 0,926, 0,919, 1,015 (`inhomogeneous`); 1,23 a 1,26 e 1,18 (`smooth`) | tese recentrada (§1.3); **G4** |
| M5, toda função de teste vale o mesmo nas duas pontas | `wafc_design()` para em `boundary = "interval"` (`wafc/R/design.R:261`); as funções, a reconferir no cartão da E9.2a | **G2** |
| M6, as aplicações não exercem `q ≥ 2`; o "8%" não resiste ao `gam.cv` | registrado na E6.1b e em D56 | **G3**; perguntas 71 e 72 |
| M7, faltam os concorrentes localmente adaptativos | fato | **G4**; E10.7 |
| M8, 78 páginas, tempo em `n` grande, sem incerteza, deslizes | `ms_6` em 78 páginas; "9% and 4%" contra 8% e 4% | E12, E11.3, E11.4, E10.3, E10.6 |

### 1.2 O que o parecer não viu e muda o plano

(a) **A base do intervalo já existe na dependência.** `WaveBased::wbasis()`
aceita `boundary = "interval"` (Cohen, Daubechies e Vial, `2^J` colunas como
no caso periódico), e `derivations/07-rota-intervalo.tex` (E1.8) já tem a
teoria. O remédio de fronteira com menos prova nova é ligar isso no
`wafc_design()`; o termo polinomial do parecer fica como segunda opção.
Volta a pergunta de identificabilidade do `ESTADO.md` §3 (a constante deixa
de ser uma função de escala única), que o cartão da E9.2a resolve.

(b) **As unidades em cache guardam linhas, não ajustes**
(`wafc-studies/R/run.R:118`). "Reaproveita os ajustes já feitos" não vale:
linha nova é ajuste novo. O que vale: as sementes são determinísticas por
célula, `n` e réplica, então método novo entra como rótulo novo, nas mesmas
réplicas, pareado com as linhas existentes e sem tocar as ~500 h do `gam`.
Regra do plano: **nenhuma tarefa muda a configuração de um método já
rodado.**

(c) **Onde a estrutura discrimina é a `scale` e o `snr1`, não o núcleo.** Em
`n = 1000`: `scale`, `wafc.cv1se` 0,65 contra `lasso.cv1se` 0,47 e `gam.reml`
0,00; `snr1`, `wafc.cv` 0,62 e 0,47 contra 0,06 a 0,08 de todo `cv1se`. O
G1 mede ali.

(d) **A M2 tem meia resposta pronta.** A taxa lenta sem condição de desenho
(Corolários 14 e 15) cobre o ajuste com mais colunas que observações. Ela
pode ir ao suplemento na E12, mas o corpo guarda um remark que diga isso.

(e) **O custo fixo de contexto é o maior desperdício do projeto.** O
`ESTADO.md` tem 404 KB (da ordem de 110 a 130 mil tokens) e o protocolo
manda todo chat lê-lo. Com executores abaixo do Fable isso também derruba a
qualidade. A E8 vem antes de tudo por isso.

### 1.3 Estratégia

Correção de rumo, não recomeço, em três movimentos:

1. **Quatro medições com critério declarado antes de rodar** (G1 a G4, §2),
   em paralelo com as correções que valem em qualquer rota (E10).
2. **A rota é decidida pelo autor com os quatro números na mão** (§3).
   Recomendação condicional: Rota A, com a tese recentrada em **seleção de
   estrutura com muitos blocos, a uma fração do custo do `gam`, com garantia
   minimax**, e a predição como segundo argumento, qualificado por `n` e por
   critério. É o que os dados já sustentam (`scale`: 78 s contra 440 s,
   estrutura 0,65 contra 0,00) e o que o parecer chama de "único cenário em
   que ele já vence".
3. **Reescrita para 40 páginas a partir de texto limpo** (E12), com os
   números entrando por macro gerada, não digitados.

A Rota T (teoria com `p`, `q` crescentes e esparsidade de blocos) fica fora
deste plano: prova longa e nova é o pior uso de um executor abaixo do
Fable, e o próprio parecer pede leitura humana linha a linha do que já
existe. Entra só a versão barata, a E10.5.

---

## 2. As quatro medições (gates)

Regras comuns: 30 réplicas; nas células já rodadas, as réplicas 1 a 30 da
produção (pareamento com as linhas existentes); nas células novas, chaves
novas ao fim de `seeds.keys`, que é só de acréscimo; cada gate tem o seu
`config/gate<i>.yaml`, sem tocar `core`, `arms` e `scale`; a produção de 100
réplicas só roda depois da rota, e só das linhas que forem ao artigo. Os
custos abaixo são estimativas, a medir no passo de preparação.

### G1. A estrutura resiste ao limiar igualado? (M3)

- **Linhas novas:** `wafc.raw` (sem limiar), `gam.reml` com limiar,
  `bsgl` com limiar, `klopp` com limiar, `bsgl.mcp` (grupo-MCP com BIC, pelo
  `grpreg`), `gam.select` (penalidade dupla do `mgcv`, no `k` escolhido pelo
  `gam.reml`).
- **O mesmo limiar:** a semântica de `wafc_thr_cv` e `wafc_thr_pick`
  (`wafc/R/threshold.R:581` e `:650`) sobre a norma `L_2` de cada componente
  ajustada; o cartão fixa o funcional e confere que, no WAFC, a regra nova
  reproduz o `wafc.cv1se` em 30 de 30.
- **Células:** `scale`, `inhomogeneous.snr1`, `mixed.snr1` (as difíceis),
  `inhomogeneous` e `mixed` (controle), `n = 1000`.
- **Regra por célula, declarada:** `cv1se` na `scale`; `cv` nos dois `snr1`
  (onde o `cv1se` cai a 6% a 8% em todo método).
- **Passa:** em duas das três células difíceis, a estrutura exata do WAFC
  supera a do melhor concorrente com limiar em 10 pontos ou mais, com
  diferença pareada acima de 2 erros-padrão.
- **Passa por custo:** empate (menos de 10 pontos) e tempo mediano do
  concorrente com limiar de 3 vezes ou mais o do WAFC na `scale`.
- **Não passa:** algum concorrente iguala a estrutura a menos de 3 vezes o
  tempo. A estrutura passa a ser atribuída à regra de limiar, não ao
  estimador.
- **Custo:** dominado pelas dobras do `gam`; ordem de 50 a 150 h de
  processador. **Máquina:** desktop, sozinha (o tempo é leitura do gate).

### G2. A vantagem sobrevive a componentes não periódicas? (M5)

- **Células novas:** `trend` (componentes suaves e monótonas, valores
  diferentes nas pontas) e `trend.inhom` (as não homogêneas do núcleo
  somadas a uma tendência).
- **Linhas:** `wafc.cv1se` (periódica), `wafc.int` (base do intervalo),
  `wafc.poly` (termo polinomial não penalizado por bloco), `gam.reml`,
  `gam.gcv`, `oracle`. `n` em 500 e 1000.
- **Passa**, para o melhor dos dois remédios: (i) na `trend`, a razão contra
  o `gam.reml` não passa da razão da `smooth` no mesmo `n` mais 0,05 (ou
  seja, 1,29 e 1,31); (ii) na `trend.inhom`, a razão é no máximo 1,00 num
  `n` e 1,03 no outro; (iii) nas células do núcleo, o remédio não piora o
  WAFC periódico em mais de 2% (se piorar, entra como opção, não como
  padrão).
- **Empate entre remédios:** fica a base do intervalo, que já tem teoria
  escrita.
- **Não passa:** o artigo restringe a tese a moduladoras periódicas (hora,
  dia do ano, direção), e o G3 só conta com elas.
- **Custo:** menos de 20 h. **Máquina:** notebook.

### G3. Existe aplicação que exerça o método? (M6)

- **(a) O spline certo nas duas bases atuais:** `bam` com erro AR(1); o
  `gam.cv` em dobras por bloco já está medido (E6.1b).
- **(b) Marylebone com moduladoras periódicas:** direção do vento e hora do
  dia. A especificação (resposta, lineares, moduladoras, grades, partições)
  é escrita no cartão **antes** de qualquer ajuste; 20 partições por bloco;
  o último bloco temporal fica intocado até a versão final.
- **Passa**, contra o melhor entre REML, GCV, `gam.cv` em blocos e AR(1): o
  critério da E6.1b (`docs/aplicacao-candidatas.md` §9: estrutura estável e
  interpretável; predição a menos de um erro-padrão do melhor spline), com
  "estável" fixado em 16 de 20 partições, e mais `q ≥ 2` exercido (duas
  moduladoras mantidas em 16 de 20). **Passa forte** se a razão média de
  predição for 0,99 ou menos.
- **Não passa em (b):** a E9.3c procura até três bases novas com moduladoras
  periódicas e mudança abrupta documentada, licença primeiro. Se nenhuma
  passar, não há Rota A na *Statistica Sinica*.
- **Custo:** 10 a 20 h de relógio. **Máquina:** notebook.

### G4. O WAFC resiste a um concorrente localmente adaptativo? (M7)

- **Linhas:** `gam` com `bs = "ad"` na célula `inhomogeneous` (`n` em 500 e
  1000); com `p = 1`, trend filtering aditivo. O primeiro passo do cartão é
  de viabilidade (se o `bs = "ad"` aceita `by`; que pacote dá o trend
  filtering aditivo).
- **Passa:** razão WAFC/adaptativo de 1,05 ou menos no erro de predição.
- **Não passa:** a predição sai da lista de argumentos; ficam estrutura e
  custo, e o tempo do adaptativo entra na tabela.
- **Custo:** desconhecido, a medir com 3 réplicas. **Máquina:** notebook.

---

## 3. A decisão de rota

Do autor, numa rodada só, com os quatro resultados:

| O que o gate decide | Se passa | Se não passa |
|---|---|---|
| G1, o primeiro argumento | estrutura com muitos blocos (ou "a uma fração do custo") abre o artigo | o limiar vira dispositivo geral; o artigo abre por predição e custo |
| G2, o escopo | qualquer moduladora, fronteira resolvida | só moduladoras periódicas, dito no título da seção e no resumo |
| G3, o alvo | *Statistica Sinica* fica (D5) | sem aplicação convincente, alvo abaixo (quartil a conferir na L16) |
| G4, a predição | fica como segundo argumento, qualificada por `n` | sai dos argumentos |

- **Rota A** (reforço aplicado, E11 inteira): G3 passa e pelo menos um entre
  G1 e G4 passa.
- **Rota 0** (só as correções, alvo abaixo): G3 não passa, ou G1 e G4 não
  passam. Da E11 ficam E11.1, E11.2 e E11.6.
- **Parar e reavaliar:** nenhum dos quatro passa.

---

## 4. Tarefas

Formato da coluna de execução: `modelo/esforço`. "Autor" é passo que só o
autor faz. Os arquivos da última coluna são os únicos que a tarefa escreve.

### E8. Infraestrutura (antes de tudo; depende das perguntas 70, 73 e 74)

| Tarefa | Entrega | Execução | Depende | Escreve |
|---|---|---|---|---|
| E8.1 | `ESTADO.md` e `TAREFA.md` divididos em cabeça e arquivo (§6), por script que move linhas sem reescrever; conferência: o multiconjunto de linhas antes e depois é igual | sonnet/medium | P74 | `docs/ESTADO.md`, `docs/estado/*`, `docs/TAREFA.md` |
| E8.2 | protocolo novo: `TAREFA.md` curto, modelo de cartão e de handoff (§5), `.claude/agents/` com modelo e esforço, ajustes em `instrucoes.md` §1, §7 e §8, `MENSAGENS.md` e `CLAUDE.md` | opus/medium (o orquestrador) | E8.1 | os citados |
| E8.3 | `wafc/scripts/15-digest.R`: de qualquer `summary.csv` ou `paired.csv`, uma tabela de até 40 linhas por célula, `n` e método | sonnet/medium | nada | o script |
| E8.4 | `k = 7` aberta de texto limpo (rodada aceita), uma seção por arquivo; compila; páginas contadas | sonnet/medium | P73 | `manuscript/ms_7.tex`, `supp_7.tex`, `references_7.bib`, `manuscript/k7/*` |

### E9. As medições

Cada medição tem três passos: `a` prepara e para (código, teste, fumaça,
custo medido); `b` roda por script, sem modelo, depois do aviso do autor;
`c` lê o resumo contra o critério. As três preparações que tocam
`wafc-studies/R/methods.R` e `wafc/R/load.R` correm em série, nesta ordem:
E9.1a, E9.2a, E9.4a. O resto é paralelo.

| Tarefa | Entrega | Execução | Depende | Escreve |
|---|---|---|---|---|
| E9.1a | G1 preparado: limiar para `gam` e `bsgl`, as seis linhas, testes, `gate1.yaml`, custo em 3 réplicas | opus/high | E8 | `wafc/R/threshold-competitors.R` (novo), `wafc/tests/test-threshold-competitors.R` (novo), `wafc/R/load.R`, `wafc-studies/R/methods.R`, `wafc-studies/config/gate1.yaml` |
| E9.1b | G1 rodado | script; aviso do autor | E9.1a | `wafc-studies/outputs/` (desktop) |
| E9.1c | G1 lido: tabela do critério, diferenças pareadas, tempos | sonnet/medium | E9.1b | handoff |
| E9.2a | G2 preparado: `boundary = "interval"` pela `wbasis()`, opção `poly`, células `trend` e `trend.inhom`, recuperação exata, o padrão reproduz o `HEAD`, a ordem do viés conferida numericamente | opus/high | E9.1a | `wafc/R/design.R`, `fit.R`, `reconstruct.R`, `dgp.R`, os testes deles, `methods.R`, `config/gate2.yaml`, `derivations/check/10-fronteira.R` |
| E9.2b | G2 rodado | script; aviso do autor | E9.2a | `outputs/` (notebook) |
| E9.2c | G2 lido | sonnet/medium | E9.2b | handoff |
| E9.3a | `bam` com AR(1) nas duas bases, nas partições da E6.2 | sonnet/high | E8 | `wafc/scripts/16-app-ar1.R` (novo) |
| E9.3b | marylebone com moduladoras periódicas: especificação no cartão (opus/high, o orquestrador), execução e leitura | sonnet/high | E9.3a | `wafc/scripts/17-app-periodic.R` (novo), `docs/aplicacao-candidatas.md` (seção nova) |
| E9.3c | só se E9.3b não passa: até três bases novas, licença primeiro; download só com permissão | sonnet/medium | E9.3b | `docs/aplicacao-candidatas.md` (seção nova) |
| E9.4a | G4 preparado: viabilidade, linhas, `gate4.yaml`, custo em 3 réplicas | sonnet/high | E9.2a | `wafc/R/competitors-adaptive.R` (novo), teste novo, `load.R`, `methods.R`, `config/gate4.yaml` |
| E9.4b, E9.4c | G4 rodado e lido | script; sonnet/medium | E9.4a | `outputs/`; handoff |

### E10. Correções que valem em qualquer rota (paralelas à E9)

| Tarefa | Entrega | Execução | Depende | Escreve |
|---|---|---|---|---|
| L15 | citações do parecer verificadas no PDF ou na página do periódico: Stone (1985), Meinshausen e Yu (2009), van de Geer, Bühlmann e Zhou (2011), Huang, Horowitz e Wei (2010), Ma, Carroll, Liang e Xu (2015), Sadhanala e Tibshirani (2019), Petersen, Witten e Simon (2016), Haris, Simon e Shojaie (2018), Oh, Naveau e Lee (2001); e duas minhas, de memória: Marra e Wood (2011, penalidade dupla) e Chetverikov, Liao e Chernozhukov (2021, lasso com validação cruzada). Todas `[VERIFICAR]` até lá | sonnet/medium | nada | `docs/referencias-verificadas.bib`, `docs/literatura.md` |
| L16 | por revista candidata (SS, EJS, JCGS, Technometrics, CSDA, Statistics and Computing): quartil no SCImago, teto de páginas, política de IA generativa e a declaração pedida | sonnet/low; autor nas páginas inacessíveis | nada | `docs/alvo-revista.md` |
| E10.3 | proposta, em bloco `.tex` marcado no handoff, das frases que não dependem de número novo: "different regime" no lugar de "weaker conditions"; o cone atribuído à condição de desenho; "levels" num sentido só; 8% e 4%; o "8% better in Beijing" fora do resumo; o rótulo "Adaptive spline" sem o nome de Wang, Jiang e Liu; as hipóteses que os experimentos não satisfazem, ditas; o remark de §1.2(d). Aplicação só depois do aval | opus/high (proposta); sonnet/medium (aplicação) | E8.4 | handoff; depois `manuscript/k7/*` |
| E10.4 | auditoria adversarial, em contexto novo, das provas que o parecer não conferiu (Teorema 2, Lema 3, S6 a S8), uma seção por tarefa; saída: dúvidas por linha, em ordem de gravidade, para o coautor. Não substitui a leitura humana, ordena-a | opus/xhigh | E8.4 | `docs/auditoria-provas.md` |
| E10.5 | o que, na cadeia do Teorema 1 ao 2, depende de `p` e `q`; se o enunciado admite `p_n q_n 2^J log(·)/n → 0` com as constantes rastreadas. Resultado afirmativo é conjectura até ter conferência numérica e prova | opus/high | nada | `derivations/08b-constantes-pq.md` (novo) |
| E10.6 | sensibilidade: outras famílias e filtros, outro tamanho de pedaço; só linhas do WAFC, núcleo em `n = 1000` | sonnet/medium | E9.1a | `config/sens.yaml`, handoff |
| E10.7 | VCBART no suplemento (os números existem) | sonnet/low | E8.4 | `manuscript/k7/supp-*`, tabela |

### E11. Reforço (Rota A; depois da decisão do §3)

| Tarefa | Entrega | Execução | Depende | Escreve |
|---|---|---|---|---|
| E11.1 | fronteira no produto: o remédio do G2 como padrão ou opção; teoria no suplemento (intervalo: adaptar o `07`; polinomial: lema novo, com conferência numérica antes e auditoria em contexto novo depois) | opus/high (código); opus/xhigh (lema) | G2 | `wafc/R/*`, `derivations/07` ou `10`, `manuscript/k7/supp-*` |
| E11.2 | a aplicação que passou o G3 no compêndio, com os quatro splines; tabelas e figura | sonnet/high | G3 | `wafc-studies/R/application*.R`, `config/application.yaml`, `results/` |
| E11.3 | velocidade em `n` grande. Perfil primeiro. Hipótese a medir: laço de grupo sobre a Gram, com a Gram de cada dobra obtida por subtração, o que tira `n` do custo de cada ajuste. Critério: a mesma escolha de `(J, λ, t)` e predição a 0,5% em 20 de 20 partições; meta de menos de 10 min por ajuste em `n ≈ 60 000`. Reabre D37 (pergunta 76) | opus/high (desenho); sonnet/high (código) | rota | `wafc/R/fit-gram.R` (novo), testes |
| E11.4 | bandas: reajuste pós-seleção com bootstrap em blocos, sem teorema, com a cobertura medida. Critério de entrada: cobertura pontual média de 90% ou mais a 95% nominais, núcleo, `n ≥ 1000`; abaixo disso não entram e a lacuna fica na Discussão | opus/high (desenho); sonnet/high (código) | rota | `wafc/R/bands.R` (novo), testes, `config/bands.yaml` |
| E11.5 | escala: o braço `scale` em mais `n` e um braço maior (dezenas de blocos), com as linhas que passaram o G1 | sonnet/medium; script | G1, rota | `config/scale2.yaml` |
| E11.6 | produção de 100 réplicas das linhas e células novas que vão ao artigo | script; aviso do autor | E11.1, E11.5 | `outputs/` (desktop), `results/` |

### E12. O manuscrito em 40 páginas

| Tarefa | Entrega | Execução | Depende | Escreve |
|---|---|---|---|---|
| E12.1 | esqueleto: o que fica no corpo (hipóteses, Teorema 1, Proposição 2, Teorema 2, Corolário 2), o que vai ao suplemento, orçamento de páginas por seção. O autor aprova antes da E12.2 | opus/high | E11 | `docs/esqueleto-k8.md` |
| E12.2 | uma tarefa por seção, cada uma lendo só a sua seção, o esqueleto e as tabelas: resumo, §1, §3 e §6 | opus/high | E12.1 | o arquivo da seção |
| | §2, §4, §5 e o suplemento | sonnet/high | E12.1 | o arquivo da seção |
| | §7 | opus/medium | as outras | o arquivo da seção |
| E12.3 | todo número do texto vira macro gerada dos CSV (`numbers_<k>.tex`), e um script confere texto contra tabela | sonnet/medium | E11.6 | `manuscript/tables/make-numbers.R`, `numbers_<k>.tex` |
| E12.4 | compilação e contagem de páginas | haiku/low | cada E12.2 | nada |
| | passada de coerência: referências cruzadas, notação, afirmações do resumo contra as tabelas | opus/medium | E12.3 | handoff |
| E12.5 | leitura das provas, linha a linha, com a lista da E10.4 | autor e coautor | E10.4, E12.2 | |

Depois, a E7 como está no `plano-projeto.md`, mais a declaração da L16.

---

## 5. Orquestração

**Orquestrador:** `opus/medium`, um por vez; `opus/high` só ao escrever os
cartões da E9, da E11 e da E12.1 e ao julgar um gate. É o chat principal do
`instrucoes.md` §7: o único que edita a cabeça do `ESTADO.md` e os cartões.

**Uma sessão do orquestrador é um lote, e termina:**

1. lê `CLAUDE.md`, a cabeça do `ESTADO.md` e os handoffs pendentes; mais nada;
2. integra os handoffs (uma edição por arquivo, conferida);
3. compara cada gate fechado com o critério do §2 e, se for o caso, leva a
   decisão ao autor;
4. escreve os cartões das tarefas destravadas;
5. despacha: tarefa que cabe num turno e não pede aviso do autor vai a um
   subagente, com o modelo e o esforço do cartão; as outras voltam ao autor
   como uma linha de abertura de chat;
6. atualiza o quadro, lembra do commit e para. A sessão seguinte começa de
   contexto novo.

**Cartão (`docs/tarefas/<ID>.md`, até 60 linhas), a única coisa que o
executor lê além do `CLAUDE.md`:**

```
# <ID>: <título>
modelo: sonnet   esforço: high   máquina: notebook
objetivo: <uma frase>
lê: <caminho:linhas>, ... (lista fechada; nada fora dela)
escreve: <lista fechada>
passos: <numerados, com os valores e nomes já decididos>
critério de saída: <verificável por comando; o comando>
não faz: commit, push, rodada longa sem aviso, arquivo fora da lista
handoff: docs/handoff-<ID>.md
```

**Handoff (até 60 linhas):** o que foi feito (arquivos); os números numa
tabela, com o caminho do CSV completo; o que entra no quadro; perguntas. A
mensagem final do executor tem até 10 linhas: caminho do handoff e veredito.

**Escalada:** tarefa que falha o critério de saída duas vezes no seu nível
recebe cartão novo um nível acima (esforço, depois modelo), com a causa da
falha escrita; não se repete no mesmo nível. Tarefa que fechou sem atrito
rebaixa a próxima do mesmo tipo.

**Verificação:** por script sempre que existir (testes, compilação, junção
exata, conferência numérica, texto contra tabela). Segundo modelo só para
matemática e para as afirmações do resumo, em contexto novo.

**Rodadas longas:** script com PID gravado e marcador de fim; nenhum modelo
fica à espera nem consultando; rodada cujo tempo é métrica roda na desktop,
sozinha; começam só com o aviso do autor.

---

## 6. Organização dos arquivos

```
docs/
  ESTADO.md              cabeça, até 150 linhas: onde está, quadro da fase,
                         perguntas abertas, próximo passo, o que não reabrir
  estado/
    decisoes.md          a tabela D1 em diante; consulta por grep '^| D56 '
    perguntas.md         a §4 de hoje (fechadas); as abertas vivem na cabeça
    diario-2026-09.md    os blocos datados do §2
    diario-2026-10.md
    historico.md         a §6
    etapas-fechadas.md   a tabela §2 do TAREFA.md
  plano-fase2.md         este arquivo; só o orquestrador lê
  plano-projeto.md       fase 1, congelado
  TAREFA.md              protocolo do executor, até 60 linhas
  tarefas/<ID>.md        um cartão por tarefa aberta; apagado ao integrar
  handoff-<ID>.md        como hoje
  parecer-wafc.html      referência do autor; executor não lê
.claude/agents/          um arquivo por papel, com modelo e esforço
manuscript/
  ms_7.tex, supp_7.tex   preâmbulo e \input
  k7/                    uma seção por arquivo
  numbers_7.tex          macros geradas dos CSV
wafc/R/                  código novo em arquivo novo, para tarefas paralelas
                         não disputarem arquivo
wafc-studies/config/     um yaml por gate e por braço novo
```

O que não muda: `ms_{k}`, `supp_{k}` e `references_{k}` andam juntos (a
pasta `k7/` leva o mesmo índice); as versões 1 a 6 ficam intactas; a
marcação recomeça em `colR2` sobre o texto limpo.

---

## 7. Economia de tokens

1. **Custo fixo por chat.** Hoje: `ESTADO.md` (404 KB), `instrucoes.md`
   (16 KB) e `TAREFA.md` (28 KB), da ordem de 125 mil tokens antes da
   primeira ação. Depois da E8: `CLAUDE.md`, o protocolo e o cartão, abaixo
   de 10 mil.
2. **Arquivo morto se consulta, não se lê.** `grep` pelo identificador
   (`D56`, `E6.1b`) e leitura do trecho por intervalo de linhas.
3. **Lista fechada de leitura no cartão**, com intervalos de linhas. O
   manuscrito por seção custa de 4 a 10 mil tokens por tarefa em vez dos 35
   mil do `ms_6` ou 55 mil do `supp_6` inteiros.
4. **Nunca ler:** versões `k` anteriores, PDFs compilados, `wafc/cache/`,
   CSV de réplicas, HTML. Regras de negação de leitura no
   `.claude/settings.json` reforçam (a sintaxe é conferida na E8.2).
5. **Saída de comando no contexto, até 50 linhas.** Log vai para arquivo;
   lê-se o fim ou o resumo da E8.3.
6. **Número não se digita.** Macro gerada e conferência por script (E12.3):
   elimina a classe de erro do "9% contra 8%" e a releitura de tabelas.
7. **O esforço caro vai ao cartão, não à execução.** Um cartão preciso do
   Opus deixa o Sonnet executar; cartão ambíguo custa mais que ele.
8. **Modelo pelo tipo de trabalho:** `opus/xhigh` só para prova e auditoria
   de prova; `opus/high` para desenho de código do estimador, tese e texto
   de posicionamento; `sonnet` para código especificado, rodadas, leitura de
   resultado, bibliografia e edição já aprovada; `haiku` para contagem,
   compilação e busca no arquivo morto.
9. **Sessões curtas do orquestrador**, uma por lote, sem histórico
   acumulado; ele não lê saída bruta de tarefa, só o handoff.
10. **Sem modelo à espera.** Rodada longa é script; quem lê o resultado é
    uma tarefa nova, de contexto novo.
11. **Paralelo só com arquivos disjuntos**, para não gastar rodada em
    conflito; o que é comum corre em série (§4, E9).

---

## 8. Perguntas ao autor

Todas respondidas em 2026-11-16, como recomendado: 70 a 77 são D83 a D90,
na ordem.

| # | Pergunta | Recomendação |
|---|---|---|
| 70 | Ratifica o plano e os critérios dos quatro gates (§2), antes de qualquer rodada? | sim; critério mudado depois de ver o número não vale |
| 71 | D56 reabre: o `gam.cv` em blocos e o `bam` com AR(1) entram como concorrentes na aplicação? | sim; é o primeiro pedido de um referee |
| 72 | O "8% better in Beijing" sai do resumo já, sem esperar o G3? | sim |
| 73 | A `k = 7` nasce de texto limpo (rodada `colR1` aceita) e com uma seção por arquivo? | sim; a reescrita para 40 páginas não se faz sobre 78 marcadas |
| 74 | O `ESTADO.md` e o `TAREFA.md` podem ser divididos como no §6? | sim; conteúdo movido, não reescrito |
| 75 | O parecer entra no `main`? O commit `cb58595` tem "Claude" como autor, o que fere a regra de autoria. Saída sem reescrever histórico: `git checkout origin/ccr-53b4afa8-g7xu8z -- docs/parecer-wafc.html` e commit seu | trazer o arquivo assim; o branch remoto fica ou é apagado por você |
| 76 | D37 (nenhuma aceleração) reabre na Rota A, sob o critério da E11.3? | sim, só se a rota for a A |
| 77 | Quem faz a leitura linha a linha das provas (E12.5), e quando? | um coautor, depois da E12.2, com a lista da E10.4 |

D5 (alvo) não reabre agora: decide-se no §3, com a L16 feita.

---

## 9. Riscos

| Risco | O que se faz |
|---|---|
| Nenhum gate passa | para-se antes de gastar a E11; o §3 já prevê a saída |
| Busca de aplicação vira caminho de bifurcações | especificação e critério no cartão antes do ajuste; último bloco intocado; toda base tentada é relatada |
| Prova nova com erro | nenhuma prova longa nova no plano; lema da E11.1 só com conferência numérica antes e auditoria em contexto novo depois; E12.5 é humana |
| A versão do `WaveBased` com a base do intervalo difere da fixada no `renv.lock` | antes da E9.2a, a junção exata das 58 unidades de referência (`wafc/scripts/14-junction.R`) com a versão nova |
| Tempos de máquinas diferentes na mesma tabela | tempo que vai a tabela vem só da desktop, sozinha |
| O parecer é a leitura de um modelo | o que os gates usam foi reconferido aqui (§1.1); o resto, no primeiro passo de cada cartão |
| Critério de gate mal calibrado | a pergunta 70 é o momento de ajustar; depois da rodada, não |
