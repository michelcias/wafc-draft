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
- O `manuscript/figures/marylebone.ukair.pdf` e o
  `manuscript/figures/beijing.heat.pdf` são cópias dos daqui, para o
  manuscrito compilar autocontido.
- O `.pdf` traz a data de criação, então a soma muda a cada relatório; as
  tabelas, não.

### O estudo de simulação (E4.4; 2026-10-08 a 2026-10-11)

| Arquivo | Origem no compêndio | SHA-256 |
|---|---|---|
| `tables/core/components.csv` | `outputs/core/tables/components.csv` | `062f07048841be402ee6442f92d22df6ee35cde40b7bc021753bdd2ec64d351a` |
| `tables/core/convergence.csv` | `outputs/core/tables/convergence.csv` | `3cb786b416f4857029a5a3773f11c4dace2e82376f6524e2439421ca8e03ba0e` |
| `tables/core/gam_k.csv` | `outputs/core/tables/gam_k.csv` | `f0c6622cc59c36fc7a8d91e6ef458050d0986cc5c882ab77e93044d9b4162794` |
| `tables/core/kcheck.csv` | `outputs/core/tables/kcheck.csv` | `b9c6319a1ba0e9fb590ad8f390142eae6de8a12f8671032b557499fc9de35e17` |
| `tables/core/paired.csv` | `outputs/core/tables/paired.csv` | `af81fa638af3e282569b4718d72b3c68c3de53ca308626931c61265f51afc697` |
| `tables/core/replicates.csv` | `outputs/core/tables/replicates.csv` | `46885b4c11c38ef363273e5fc57068abb77feae40e7c3ac158dbbe6dee14ec1c` |
| `tables/core/slopes.csv` | `outputs/core/tables/slopes.csv` | `777146af441c1a12a4f8f35ab5da5fcc9a08184db901878305d339b0d8b315d8` |
| `tables/core/summary.csv` | `outputs/core/tables/summary.csv` | `ff1c9694288af912af4ff742cc92baf51cdf5b4fb86972394971d2c408a10056` |
| `tables/core/threshold.csv` | `outputs/core/tables/threshold.csv` | `845ff7c96889e83a233b168ce9944513866fc568b987cd07cb5e5c5930a38820` |
| `tables/arms/components.csv` | `outputs/arms/tables/components.csv` | `c892db16f8812016bd4d27c900d874c105fc1f39d08caee0e43d2e271a80aea3` |
| `tables/arms/convergence.csv` | `outputs/arms/tables/convergence.csv` | `927421e89221c548f9c9474490036175f3ff50a98f368607e9dda2986030aa99` |
| `tables/arms/gam_k.csv` | `outputs/arms/tables/gam_k.csv` | `765732dbe11939031af7ce6d50dad6e03d960b1f061574f23a594114988591a6` |
| `tables/arms/kcheck.csv` | `outputs/arms/tables/kcheck.csv` | `c2a7cb18a83afaa4d8c40ba5674ecc292adb3a012e6f93268906ecb4bd42a5b7` |
| `tables/arms/paired.csv` | `outputs/arms/tables/paired.csv` | `1df83bcec0775c6916e899ab614ea86f8245e1effcbf09569ad8600b55ca908b` |
| `tables/arms/replicates.csv` | `outputs/arms/tables/replicates.csv` | `21a2e1429b9b177680417d09fc3a68305ef765d212ff0c3f5393617cd19780b8` |
| `tables/arms/summary.csv` | `outputs/arms/tables/summary.csv` | `4b4d1859055882d5378da8797298d6eebb03256aba4985b03494f44d721e50d8` |
| `tables/arms/threshold.csv` | `outputs/arms/tables/threshold.csv` | `25a2efac1e3a47b37076a3faf8b1150183b0b4dcf38aa0a476aa0a51db85a37d` |
| `tables/scale/components.csv` | `outputs/scale/tables/components.csv` | `9912b436ef584da5f6e123e9389425d339bbecc0b887084b6c5939273bfcbf99` |
| `tables/scale/convergence.csv` | `outputs/scale/tables/convergence.csv` | `2aaedb65d497aa9ac37475c76d5aebde3428c8a0837f19c7793b4feef529cea4` |
| `tables/scale/gam_k.csv` | `outputs/scale/tables/gam_k.csv` | `6ba892ffc9c50797c4f95aa4c5b8550eba5232999bafbc4187ece6bb732aa478` |
| `tables/scale/kcheck.csv` | `outputs/scale/tables/kcheck.csv` | `d650055f48ccf624eacf01309d775fedd15bf42bf20d57401daf55934136b836` |
| `tables/scale/paired.csv` | `outputs/scale/tables/paired.csv` | `cd9a0778dd328828787b497d71f940ce1361d7bd51f6a40768da99220b257b26` |
| `tables/scale/replicates.csv` | `outputs/scale/tables/replicates.csv` | `63efb8a6eaac710e7adce223c33fc46b359932e7ad0923551e553e2594061cf3` |
| `tables/scale/summary.csv` | `outputs/scale/tables/summary.csv` | `1b440a7c362e141c18b41baf6ea0707daa1b508db782da84d7586c3ca322b50f` |
| `tables/scale/threshold.csv` | `outputs/scale/tables/threshold.csv` | `22a9b58f265029e52f2efcfcd1f702a27cb4869434b4e13043a107af6a4651f4` |

- **Os `curves.csv` não vêm** (85 a 256 MB cada, 21 a 58 MB comprimidos):
  ficam em `wafc-studies/outputs/<desenho>/tables/` na desktop do autor;
  o que deles for ao manuscrito (a Figure 1) vem como figura. As somas:
  `core` `25201674e3cfa5919d79f3ca90565040853ab509a1b96cfbc05d92f39034ba00` (147 MB);
  `arms` `ab178bf32cd63ccc3d4d8c1972419eeb814f1af7b2fd3f2a0ccc52e361681c89` (256 MB);
  `scale` `4bc0b1676231b5b65235547b0dbeafa475c5411a6b95daaece57754a9a083efe` (85 MB).
- **As unidades:** 19 600 do `core`, 5 700 dos braços e 900 da `scale`,
  0 falhas, na desktop do autor (i9-10900KF, 10 trabalhadores, BLAS de
  referência em uma thread; `outputs/machine.txt`), com o código do
  método na árvore `34f2e04` de `wafc/R` (`sessionInfo.txt` de cada
  desenho) e a configuração do commit `aa11e76` (D82: o `gam.reml` em
  `n = 1000` até `k = 160`). A junção com o notebook deu 58 de 58
  unidades idênticas fora o tempo (D81).
- **A rodada:** de 2026-10-08 17h34 a 2026-10-11 ~04h00 (o `core` até
  19 598 unidades, parado por D82), e de 2026-10-11 04h11 a 20h51 (as
  500 unidades do `gam.reml` em `n = 1000` reajustadas, os braços e a
  `scale`). Processador: 498,6 h no `core`, 126,6 h nos braços e 22,7 h
  na `scale`, 647,9 h no total, sem as 34,5 h das 498 unidades
  descartadas por D82 nem as duas que não terminaram.
- **A agregação:** `scripts/02_aggregate.R` de cada desenho, chamado pelo
  `wafc/scripts/14-production.sh` logo depois do ajuste.
