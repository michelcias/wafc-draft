# results

Cópia de referência do que o compêndio `wafc-studies` gera. Nunca gerado
aqui. Quando divergir do compêndio, o compêndio é a verdade.

## Origem

### Aplicação `marylebone.ukair` (E6.2, E5j; 2026-10-07)

| Arquivo | Origem no compêndio | SHA-256 |
|---|---|---|
| `figures/marylebone.ukair.pdf` | `outputs/application/figures/marylebone.ukair.pdf` (Figure 2 do `ms_5`) | `4a41b38c831552377dcdb9b90e4ca659b5240a504faed188549a10ea7fa4c3ad` |
| `tables/prediction-marylebone.ukair.csv` | `outputs/application/tables/` | `0cbbdf99c8d742ee585700cbed5d07fa43216f2db19a52211701411d1124afcb` |
| `tables/structure-marylebone.ukair.csv` | idem | `6c4aa43bce52592a7d2765f14639c00d24f3c66f097e945e3398331122eafcd0` |
| `tables/choices-marylebone.ukair.csv` | idem | `58636030a49f7b827ebb2fe82a2326d28ef5ea9a12dc172d70252930613afcbe` |
| `tables/readings-marylebone.ukair.csv` | idem | `fae0134f44da37062be74af86b258971e6396834a0c29cc90fd5c678eb5d805f` |

- **As unidades:** a rodada de marylebone (84 unidades, `split00` a
  `split20` dos métodos `wafc`, `gam.reml`, `gam.gcv` e `linear`; 0 falhas),
  terminada em 2026-10-07 às 02h35, com o código do método na árvore
  `e262faa` de `wafc/R` (`sessionInfo.txt` da rodada), o mestre `20261006` e
  a configuração `config/application.yaml`.
- **O relatório:** refeito em 2026-10-07, sem ajuste, com
  `Rscript scripts/03_application.R --parts=report --bases=marylebone.ukair`,
  depois da mudança da figura de D73 (curvas e faixa interrompidas nos
  buracos da data mais longos que `2^{-J}` do período, ~10,7 dias em
  `J = 8`, e o rug). As quatro tabelas são idênticas, bit a bit, às da
  rodada de 02h35; só a figura mudou. O `.png` da figura (não copiado) é
  idêntico ao do esboço aprovado pelo autor.
- O `manuscript/figures/marylebone.ukair.pdf` é cópia deste
  `figures/marylebone.ukair.pdf`, para o manuscrito compilar autocontido.
- O `.pdf` traz a data de criação, então a soma muda a cada relatório; as
  tabelas, não.
