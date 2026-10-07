# results

Cópia de referência do que o compêndio `wafc-studies` gera. Nunca gerado
aqui. Quando divergir do compêndio, o compêndio é a verdade.

## Origem

### As duas aplicações (E6.2, E5j, E6.2b; 2026-10-07)

| Arquivo | Origem no compêndio | SHA-256 |
|---|---|---|
| `figures/marylebone.ukair.pdf` | `outputs/application/figures/marylebone.ukair.pdf` (Figure 2 do `ms_5`) | `fee33368f9bf65380645b1aff550ba6487e1f310946bc3a9040ab8c92be269b7` |
| `tables/prediction-marylebone.ukair.csv` | `outputs/application/tables/` | `0cbbdf99c8d742ee585700cbed5d07fa43216f2db19a52211701411d1124afcb` |
| `tables/structure-marylebone.ukair.csv` | `outputs/application/tables/` | `f99ed30eebb4bb5ffac0eae3e927c52f7e33158533be1536710385fd4ee17e16` |
| `tables/choices-marylebone.ukair.csv` | `outputs/application/tables/` | `58636030a49f7b827ebb2fe82a2326d28ef5ea9a12dc172d70252930613afcbe` |
| `tables/readings-marylebone.ukair.csv` | `outputs/application/tables/` | `c320d2fe7995785c115d0af9b1ed660e3bbb57e006a03166d95e305319aeb0e1` |
| `figures/beijing.heat.pdf` | `outputs/application/figures/beijing.heat.pdf` (Figure S1 do `supp_5`) | `341bd9301a36a986fa8f69293b2fe9b07239e59ecd92058252415293aca65501` |
| `tables/prediction-beijing.heat.csv` | `outputs/application/tables/` | `d5f11c2ea3817d0205ea6758adc33abd8a60a9dd254522721b1a057705a7bd80` |
| `tables/structure-beijing.heat.csv` | `outputs/application/tables/` | `567996c4ccc47524ffc493b48df5069f9b385c38048bc92b351f1cdb70ba8b3a` |
| `tables/choices-beijing.heat.csv` | `outputs/application/tables/` | `f8e79f3a42313ee509b766d3874e563bf50154b53034f52918d51c6fe763ce11` |
| `tables/readings-beijing.heat.csv` | `outputs/application/tables/` | `eaf6c755babe4da96086669ec2ba3aa362b3b5d1172a09f4de475b7c6d579e48` |

- **As unidades:** as 168 da aplicação (`split00` a `split20` dos métodos
  `wafc`, `gam.reml`, `gam.gcv` e `linear` nas duas bases; 0 falhas):
  marylebone terminada em 2026-10-07 às 02h35, a `beijing.heat` às 16h31,
  com o código do método na árvore `e262faa` de `wafc/R`
  (`sessionInfo.txt` da rodada), o mestre `20261006` e a configuração
  `config/application.yaml`.
- **O relatório:** o da cadeia, em 2026-10-07 às 16h31, com
  `Rscript scripts/03_application.R --parts=report` e o
  `R/application_report.R` da E6.2b (D73: curvas e faixa interrompidas nos
  buracos da data mais longos que `2^{-J}` do período, ~10,7 dias em
  `J = 8`, e o rug; D75(b), (c): as colunas `range_`, `jump_` e as leituras
  só nos pontos desenhados). Em marylebone mudaram, contra a cópia da E5j,
  só `structure-` e `readings-` (a subida do ajuste inteiro de 8,876 a
  8,963).
- **A Figura S1** refeita em 2026-10-07 depois do piloto, só o relatório da `beijing.heat`: o rug ao longo do dia do ano e da umidade (D75(e)) e a umidade cortada a ~15% e ~95% (`trim: rh: [6, 4]`), onde a base periodizada dobra (D67(e)); muda só o `structure-beijing.heat.csv` (as colunas `range_` e `jump_` ao longo da umidade).
- O `manuscript/figures/marylebone.ukair.pdf` é cópia deste
  `figures/marylebone.ukair.pdf`, para o manuscrito compilar autocontido.
- O `.pdf` traz a data de criação, então a soma muda a cada relatório; as
  tabelas, não.
