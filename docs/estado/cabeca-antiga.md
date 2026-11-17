# ESTADO antigo: cabeçalho, §1, §2 (abertura), §3 e §5
Movido de docs/ESTADO.md em 2026-11-16, sem alteração de linhas.
Arquivo morto: não ler inteiro; consultar por grep pelo identificador (D56, E6.1b, pergunta 33) e ler o trecho.

<!-- conteudo movido abaixo -->
# Estado do trabalho, handoff de continuidade

**Última atualização:** 2026-11-16.
**Fase 2 (2026-11-16):** parecer interno no papel de referee recebido
(`docs/parecer-wafc.html`, commit `cb58595` do branch
`origin/ccr-53b4afa8-g7xu8z`, fora do `main`); o plano que responde a ele
está **proposto** em [`plano-fase2.md`](plano-fase2.md) e espera a
ratificação do autor (perguntas 70 a 77, §8 de lá). Nada da fase 2 foi
executado; a E7 passa a depender dela.
**Etapa corrente:** **E0 a E3 inteiras; E4.1 a E4.4; E5a, E5c a E5k; E6.1a
a E6.2b; L1 a L14 fechadas. A produção do estudo de simulação (E4.4)
rodou na desktop do autor (D81, D82), e as tabelas estão em
`results/tables/{core,arms,scale}/`. O próximo passo é a E5b, os números
do estudo no manuscrito.** O WAFC é o block LASSO na forma balanceada seguido do
limiar `cv1se` (D44, D45); o alvo é a *Statistica Sinica* (D5). A teoria
em blocos está provada (o Corolário 11 de `08-blocos.tex` é o Theorem 1,
D47), com a taxa lenta sem condição de desenho e a cota inferior que a faz
ótima em `π ≥ 2` (E1.9). **O compêndio** é a pasta `wafc-studies/` (D60),
que vira repositório só no fim, por cópia; o código do método fica em
`wafc/` e vai ao `WaveBased` no fim (D61). **As duas aplicações rodaram**
(168 unidades, 0 falhas; marylebone no corpo, beijing.heat no supp; os
números na `k = 5`). **O desenho do estudo** está fixado (D60, D62, D63,
D69, D70, D77 a D80) e o piloto rodou (E4.3: a escala entra; a E4.4 custa
~592 h de processador, ~74 h em 8 trabalhadores, ~45 a 60 h nos 10 da
desktop). O teto de páginas fica para o fim (D21).
**Versão viva do manuscrito:** `k = 6` (`manuscript/ms_6.tex`,
`supp_6.tex`, `references_6.bib`), aberta pela E5b em 2026-10-12 como cópia
da `k = 5`; os `\input` de `tables/` e `figures/` que a `k = 6` muda ganham
o sufixo `_6` (o primeiro: `tables/app-tuning_6.tex`), para a `k = 5`
compilar como estava. As Seções 5 e 6 estão escritas, a 6 com os números
finais e a 5 com marcadores `[E4.4]`, que a E5b preenche com as tabelas da
E4.4; a Discussion (Seção 7) não existe. Os dois compilam limpos em **73 e
101 páginas** com a marcação, contra o teto de 40 (D21). As versões `k = 1`
a `k = 5` ficam intactas.
**Cor da rodada corrente:** `colR1`, em uso desde `k = 2` (nenhuma rodada
foi aceita, então toda a marcação segue em `colR1`).

Este documento é o ponto de partida de cada sessão. Ele diz onde o trabalho
parou, o que já foi decidido (para não reabrir) e o que vem a seguir.

---

## 1. Para começar uma sessão

1. Ler `CLAUDE.md` (as duas regras de abertura: sem coautoria, respostas
   curtas).
2. Ler este arquivo.
3. Ler [`instrucoes.md`](instrucoes.md). É a regra, não a sugestão.
4. Ler a etapa em foco em [`plano-projeto.md`](plano-projeto.md). Em E1,
   também [`proposta-metodo.md`](proposta-metodo.md).
5. Se for escrever fórmula: [`notacao.md`](notacao.md) antes, e as Seções 2,
   3 e 5 de `../wall-manuscript/manuscript/theo/ms_theo_1.tex`, que é o
   molde da prova.

### Mensagem de abertura sugerida

Chat principal (o que edita este arquivo e commita):

> Leia `CLAUDE.md`, `docs/ESTADO.md` e `docs/instrucoes.md` para se situar.
> Hoje vamos <foco do dia>.

Outra máquina (instalação e disposição das pastas em `CONTINUAR.md`):

> Leia `docs/CONTINUAR.md` e me diga o que falta instalar; depois leia
> `docs/ESTADO.md` e continuamos de onde parou.

Chat de tarefa (uma etapa, arquivos restritos, termina em handoff):

> Leia `docs/TAREFA.md` e execute a tarefa <etapa>.

### Git

O autor commita e faz push. O assistente lê o estado e lembra de commitar nos
momentos apropriados. **Sem coautoria nas mensagens**, sem exceção.

---

## 2. Onde o trabalho está

## 3. O que NÃO reabrir

- **Pasta dentro do `wafc-draft`, `WaveBased` ou repositório próprio para o
  código.** Respondido em D4 e D7, detalhe em `plano-projeto.md` E0.2: o
  método em `wafc/` aqui, compêndio próprio. O empacotamento volta como
  decisão só em E3.3, com o código testado.
- **Se a restrição de identificabilidade precisa de multiplicador ou
  centralização numérica.** Não precisa no caso periódico com `j0 = 0`: a
  constante é a única função de escala e é descartada (D2). Só volta se
  `boundary = "interval"` virar o padrão.

---

## 5. Próximos passos

**Fase 2 (2026-11-16), antes do que segue.** O parecer interno
(`docs/parecer-wafc.html`, no branch `origin/ccr-53b4afa8-g7xu8z`) conclui
rejeição provável na *Statistica Sinica* no estado da `k = 6` e pede
correção de rumo. Conferido aqui contra `results/tables/`: as razões
WAFC/`gam.reml` (0,919 a 1,256) e o `lasso.cv1se` com estrutura exata 1,00
nas cinco células em `n = 1000`. O plano está em
[`plano-fase2.md`](plano-fase2.md): quatro medições com critério declarado
(G1, o limiar igualado; G2, componentes não periódicas; G3, a aplicação;
G4, o concorrente adaptativo), a rota decidida depois delas, modelo e
esforço por tarefa (no máximo Opus) e a reorganização dos arquivos (E8).
**Pede resposta do autor:** as perguntas 70 a 77, no §8 do plano. O teto
de páginas (D21) passa à E12 de lá. Achados que mudam o plano: o
`WaveBased::wbasis()` já tem `boundary = "interval"`; as unidades em cache
guardam linhas, não ajustes; este arquivo tem 404 KB e é o maior custo fixo
de contexto.

**Onde parou (2026-10-12).** Nenhum handoff pendente; **catálogo
vazio**; nenhum processo rodando. **A E4.4 fechou** (bloco de 2026-10-12
no §2; D82): 26 200 unidades, 0 falhas, 647,9 h de processador; as
tabelas em `results/tables/{core,arms,scale}/`. A pergunta 69 fica
respondida por D82; a 68 foi aceita como recomendada em 2026-10-12 e
aplicada na `k = 6`. **Nenhuma pergunta em aberto pede resposta do
autor.**

**A E5b está em curso na `k = 6`** (2026-10-12). Feito: a Table 1
(predição e estrutura exata em `n = 1000`; o ISE foi à S1 porque, em
espaço duplo, três painéis passavam da página), as Tabelas S1 a S11 (uma
por célula com RMSE, ISE e estrutura em todo `n`; nulos mantidos; pares;
braços e modelo maior; sintonia e tempo; inclinações), geradas por
`manuscript/tables/make-sim-main.R` e `make-sim-supp.R`; a leitura da
Seção 5.4 e da S9, a frase da grade do REML no 5.2 e a subseção "The grids
of the splines" na S9 (D82), o processador. O enquadramento combinado com o
autor: a estrutura exata e a predição lado a lado, as concessões ditas e
explicadas (células suaves, `n = 2000`, `snr1`, o REML no topo em
`n = 500`). A Figure 1: as medianas ponto a ponto das 100 réplicas do
WAFC e do GAM (REML) em bumps e blocks (o WAFC alcança de 14% a 44% mais
da altura de cada um dos 11 picos; os dois abaixo da verdade em [0,3; 0,4]),
e a Figure S1 com as faixas de 5% a 95% (largura média quase igual), as
duas de `manuscript/figures/make-sim-components.R`, que lê os `curves.csv`
da desktop; a figura de Dongsi passa a S2. A Seção 7 (Discussion: o que o
método entrega, onde não ganha e por que, o que fica aberto), a frase da
simulação e da aplicação no resumo (D68(i)), o roteiro da Introdução, e a
oração do limiar de menor erro de validação cruzada no `snr1` (estrutura
em 62% e 47%, à frente do REML por 4% e 3%) na 5.4. 78 e 114 páginas.
Falta: D21, o teto de páginas.

**O que era o plano da E5b**, os números do estudo no manuscrito, a
partir de `results/tables/`:

1. Os marcadores `[E4.4]` das Seções 5 e S9: a Table 1, a Figure 1 (dos
   `curves.csv`, que estão só na desktop), as Tabelas S1 a S4, o
   processador e o tempo (647,9 h, a máquina de D81).
2. A frase sobre o topo das grades dos splines (D78, D80, D82): o
   `gam.reml` no topo admissível na `mixed` e em `n = 2000`, o 240 fora
   de `n = 1000` pelo custo, o `gam.gcv` preso pelo custo.
3. A Seção 7 (Discussion), o resumo com a simulação e a aplicação
   (D68(i)), os itens de texto pendentes (D76; a pergunta 68 se aceita).
4. **A decisão do teto de páginas** (D21): 66 páginas sem o removido
   contra 40, sem a Discussion; o que vai ao supp sem reescrever está na
   pergunta 60(h).
5. **E7**, a submissão: o `wafc/` ao `WaveBased` (D61), a extração do
   compêndio (D60(f), com o `LICENSE` e a lista das árvores de código de
   cada rodada no `PROVENANCE.md`), os DOIs, o checklist do
   `alvo-revista.md` §6.

**Ainda falta trazer do notebook para a desktop** (não versionado; o
comando no `CONTINUAR.md` §2b, passo 6): os PDFs de `refs/` (114 MB, só
para conferir citação) e as unidades da aplicação,
`wafc-studies/outputs/application/` (3,5 MB, para refazer o relatório das
figuras na E5b sem reajustar). Os `curves.csv` do estudo estão só na
desktop, então a Figure 1 da E5b se faz lá. O `wafc/cache/` do notebook
(1,1 GB, pilotos e sondagens) fica lá; só os números registrados neste
arquivo importam.

**Sem bloquear:** 12 (o cenário `smooth`), 13 (o nome da (BD)), 14 (o
KP como opção do `wafc()`), 21 (`p` crescente), 27 (o acabamento da
semana de submissão, E7), 32 (as pendências de E5c, entre elas o Lema 10
no supp), 35, 36 e 39 (pendências de L5, E2.5d e L8); as opções da figura da
`beijing.heat` no compêndio (sem a linha do nível, proporção mais baixa;
bloco da E5k), para a E5b; `boundary = "interval"` só se um referee pedir.

