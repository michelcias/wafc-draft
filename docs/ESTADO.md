# Estado do trabalho

**Última atualização:** 2026-11-16.
**Fase:** 2, a resposta ao parecer interno ([`parecer-wafc.html`](parecer-wafc.html)).
O plano é o [`plano-fase2.md`](plano-fase2.md), ratificado (D83 a D90).
**Etapa corrente:** a E8 (infraestrutura) fechou em 2026-11-16. O próximo
lote abre a E9.1a (preparar o G1) e, em paralelo, E9.3a, L15, L16 e E10.5.
**Versão viva do manuscrito:** `k = 7` (`manuscript/ms_7.tex`, `supp_7.tex`,
`references_7.bib`; o texto em `manuscript/k7/`, uma seção por arquivo). É o
texto limpo da `k = 6`, com a rodada `colR1` aceita (D86): **71 e 111
páginas**, contra o teto de 40. As versões 1 a 6 ficam intactas.
**Cor da rodada corrente:** `colR2`.

Este arquivo é a cabeça: diz onde o trabalho está, não como chegou lá. Tem
teto de 150 linhas; o que deixa de ser corrente vai ao diário (§7).

---

## 1. Para começar

- **Orquestrador** (Opus, esforço médio; um por vez): este arquivo, os
  `docs/handoff-*.md` pendentes e, do plano, só a seção da tarefa em foco.
- **Executor:** [`TAREFA.md`](TAREFA.md) e o cartão (`docs/tarefas/<ID>.md`).
- As regras estão no [`instrucoes.md`](instrucoes.md), §7 e §8; as mensagens
  de abertura, no [`MENSAGENS.md`](MENSAGENS.md).

## 2. Onde o trabalho está

- **O método** (D44, D45): o block LASSO na forma balanceada seguido do
  limiar `cv1se`; `cv.wafc(x, u, y)` em `wafc/` (1 770 testes).
- **A teoria** (E1): o Theorem 1 é o Corolário 11 do `08-blocos.tex` (D47);
  taxa lenta sem condição de desenho; cota inferior que a faz ótima em
  `π ≥ 2`. O mapa dos resultados está em `numeracao-resultados.md`.
- **O estudo** (E4.4, D82): 26 200 unidades, 0 falhas, 647,9 h de
  processador; tabelas em `results/tables/{core,arms,scale}/`. O cache das
  unidades e os `curves.csv` estão só na desktop.
- **As aplicações** (E6.2): marylebone no corpo, beijing.heat no supp; 168
  unidades.
- **O parecer** (2026-11-16): rejeição provável na *Statistica Sinica* no
  estado da `k = 6`; as objeções M1 a M8 estão resumidas no plano, §1.
- **A fase 2:** quatro medições com critério fixado antes de rodar (G1, o
  limiar igualado; G2, componentes não periódicas; G3, a aplicação; G4, o
  concorrente adaptativo), depois a rota (plano §3), depois E11, E12 e E7.
- **O que orienta a fase:** `WaveBased::wbasis()` já tem
  `boundary = "interval"`; o cache guarda linhas, não ajustes, então método
  novo é rótulo novo nas mesmas réplicas; a estrutura só discrimina na
  `scale` e no `snr1`; nenhuma tarefa muda a configuração de um método já
  rodado.

## 3. Quadro da fase 2

Modelo e esforço como no plano, §4. "Pronta" quer dizer sem dependência
aberta e ainda sem cartão.

| Tarefa | Estado | Execução | Espera |
|---|---|---|---|
| E8.1 a E8.4, infraestrutura | **fechadas** em 2026-11-16 | | |
| E9.1a, preparar o G1 | pronta | opus/high | |
| E9.1b, E9.1c, rodar e ler o G1 | | script (desktop); sonnet/medium | E9.1a e o aviso |
| E9.2a, preparar o G2 | | opus/high | E9.1a (arquivos comuns) |
| E9.3a, `bam` com AR(1) | pronta | sonnet/high | |
| E9.3b, marylebone periódica | | sonnet/high | E9.3a |
| E9.4a, preparar o G4 | | sonnet/high | E9.2a |
| L15, citações do parecer | pronta | sonnet/medium | |
| L16, revistas e política de IA | pronta | sonnet/low; autor | |
| E10.3, frases sem número novo | pronta | opus/high | |
| E10.4, auditoria das provas | pronta | opus/xhigh | |
| E10.5, constantes em `p` e `q` | pronta | opus/high | |
| E10.6, sensibilidade | | sonnet/medium | E9.1a |
| E10.7, VCBART no supp | pronta | sonnet/low | |
| a rota | | autor | G1 a G4 |
| E11, E12, E7 | | | a rota |

## 4. Perguntas que pedem o autor

Nenhuma. As últimas (70 a 77) viraram D83 a D90. A próxima pergunta nova é
a 78; a próxima decisão, D91.

## 5. Próximo passo

1. **Autor:** commitar a E8. Está tudo na árvore de trabalho, nada no
   índice.
2. **Orquestrador, lote seguinte:** escrever os cartões de E9.1a, E9.3a,
   L15, L16 e E10.5 e despachar. A E9.1a define o custo do G1, que roda na
   desktop.
3. **Antes da E9.2a:** a junção exata das 58 unidades de referência
   (`wafc/scripts/14-junction.R`) com a versão do `WaveBased` que tem a base
   do intervalo.
4. **Ao regerar tabela ou figura para a `k = 7`:** os geradores
   (`manuscript/tables/make-*.R`, `figures/make-*.R`) ainda escrevem
   `\color{colR1}` no float; trocar por `colR2` ou tirar (E12.3).
5. **Sem bloquear, herdadas da fase 1:** as perguntas 12, 13, 14, 21, 27,
   32, 35, 36 e 39, em `estado/perguntas.md`.

## 6. O que não reabrir

- **Pasta dentro do `wafc-draft`, `WaveBased` ou repositório próprio para o
  código.** Respondido em D4, D7 e D61: o método em `wafc/` aqui até o fim,
  o compêndio na pasta `wafc-studies/`.
- **Se a restrição de identificabilidade precisa de multiplicador ou
  centralização numérica.** Não precisa no caso periódico com `j0 = 0` (D2).
  Volta na E9.2a, e só lá, porque a base do intervalo muda a função de
  escala.
- **Os critérios dos quatro gates** (D83): fixados antes de rodar.
- O resto está na tabela de decisões, D1 a D90.

## 7. Arquivo

Consulta por `grep` do identificador e leitura do trecho; nunca inteiro.

| O que | Onde |
|---|---|
| decisões, D1 em diante (tabela viva: decisão nova entra ao fim) | `estado/decisoes.md` |
| o que cada etapa fechada mediu, um bloco por etapa | `estado/diario-2026-09.md`, `-10.md`, `-11.md` |
| perguntas fechadas da fase 1 (1 a 69) | `estado/perguntas.md` |
| etapas fechadas da fase 1, em tabela | `estado/etapas-fechadas.md` |
| a cabeça e o protocolo antigos, o histórico de sessões (congelados) | `estado/cabeca-antiga.md`, `tarefa-antiga.md`, `historico.md` |
