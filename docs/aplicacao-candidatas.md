# Candidatas para a aplicação (E6.1a)

Sondagem com veredito, não escolha final. A tarefa é decidir **qual base
sustenta a aplicação de E6**, e o critério não é o usual.

**Última atualização:** 2026-09-20. Números produzidos por
`wafc/scripts/05-sondagem-aplicacao.R`, semente `20260920`.

---

## 1. O critério, e por que ele não é o usual

A escolha de uma base de aplicação costuma pesar interpretabilidade, tamanho
e disponibilidade. Aqui há uma exigência anterior a essas três: **a
modulação tem de ser não homogênea**.

A razão está em `proposta-metodo.md` §6 e em `busca-novidade.md`. Xue & Yang
(2006) provam a taxa ótima univariada para o modelo de coeficientes aditivos
por splines polinomiais **sob `α_ls ∈ C^{p+1}[0,1]`**, isto é, na classe em
que splines são ótimos. Se as `β_ℓ(u)` da aplicação forem suaves nesse
sentido, o spline não é um competidor a ser batido: ele é o estimador ótimo,
e o WAFC no máximo empata. O argumento do artigo, a adaptação a regularidade
espacialmente heterogênea que a base de wavelets compra e a base de splines
de nós fixos não compra, só tem onde morder se a modulação
verdadeira tiver quina, salto ou mudança brusca de escala.

Daí o critério de aceitação desta sondagem, em três partes, com o peso nesta
ordem:

1. **Estrutura**: a componente estimada tem energia em níveis finos e
   variação localizada, e não a forma que um spline de poucos nós daria.
2. **Predição**: o WAFC não perde para o `mgcv::gam` com
   `s(u_m, by = x_ℓ)` fora da amostra. Empatar em predição e ganhar em
   estrutura é resultado bom; perder em predição é fatal, porque o referee
   lê a tabela antes da figura.
3. **Interpretação**: `Y`, `X_ℓ` e `U_m` que um leitor de estatística geral
   reconheça sem parágrafo de contexto.

Uma base que falhe (1) é **favorável ao spline** e não deve ser a aplicação,
por melhor que seja nos outros quesitos.

---

## 2. As três candidatas

Todas públicas, citáveis, com `n` na casa das dezenas de milhares. Os dados
**não são versionados** (`.gitignore` tem `wafc/cache/`): o script os baixa
para `wafc/cache/data/` na primeira execução e confere a soma SHA-256 de
cada arquivo contra a registrada abaixo e no próprio script.

### A. Bike sharing, Capital bikeshare, horário

| Item | |
|---|---|
| Fonte | UCI Machine Learning Repository, conjunto 275, `bike+sharing+dataset.zip`, arquivo `hour.csv` |
| URL | `https://archive.ics.uci.edu/static/public/275/bike+sharing+dataset.zip` |
| DOI do conjunto | `10.24432/C5W894` |
| SHA-256 do `.zip` | `b70182d0d0508e9abbb79306ce5c0cec34869000f8220175ac83d11dbe845401` |
| Licença | Creative Commons Attribution 4.0 International (CC BY 4.0), declarada na página do conjunto |
| `n` | 17 379 horas (2011-01-01 a 2012-12-31), sem faltantes |
| Citação | Fanaee-T & Gama (2014), conferida (§7) |

Mapeamento:

| Papel | Variável |
|---|---|
| `Y` | `log(cnt)`, log da contagem horária de aluguéis |
| `X_1` | constante (o nível, `β_1(u)` é a superfície de base) |
| `X_2` | `temp`, temperatura normalizada |
| `X_3` | `hum`, umidade relativa |
| `X_4` | `windspeed`, velocidade do vento |
| `U_1` | `hr`, hora do dia (0 a 23) |
| `U_2` | dias decorridos desde 2011-01-01 (0 a 730) |

A leitura substantiva é `β_temp(hora, dia) = c + g_{temp,hora}(hora) +
g_{temp,dia}(dia)`: o efeito da temperatura sobre a demanda, modulado pela
hora do dia e pela posição na série. A quina esperada está nos horários de
início e fim da jornada.

**Ressalva de desenho, registrada antes dos números:** `hr` tem 24 valores
distintos, de modo que o bloco `(ℓ, 1)` tem posto no máximo 23 qualquer que
seja `J`, e para `2^J − 1 > 23` o bloco é sobreparametrizado por
construção. Não invalida o ajuste (o LASSO lida com isso), mas significa
que a "não homogeneidade em hora" que se pode detectar é a de uma função
sobre 24 pontos, não a de uma função sobre um intervalo. O mesmo limite
restringe o competidor: o `mgcv` recusa `k` maior que o número de valores
distintos.

### B. Beijing multi-site air quality, sítio Dongsi, horário

| Item | |
|---|---|
| Fonte | UCI Machine Learning Repository, conjunto 501, `beijing+multi+site+air+quality+data.zip`, arquivo `PRSA_Data_Dongsi_20130301-20170228.csv` |
| URL | `https://archive.ics.uci.edu/static/public/501/beijing+multi+site+air+quality+data.zip` |
| DOI do conjunto | `10.24432/C5RK5G` |
| SHA-256 do `.zip` | `b04da438b2f331ac0ffd45aebdfec0d20d2367feb5f6948c4b1f7ce1191e33c4` |
| Licença | Creative Commons Attribution 4.0 International (CC BY 4.0), declarada na página do conjunto |
| `n` | 35 064 horas no sítio Dongsi; 34 294 casos completos nas variáveis usadas (420 768 no conjunto inteiro, 12 sítios) |
| Citação | Zhang, Guo, Dong, He, Xu & Chen (2017), conferida (§7) |

Mapeamento:

| Papel | Variável |
|---|---|
| `Y` | `log(PM2.5)` |
| `X_1` | constante |
| `X_2` | `WSPM`, velocidade do vento (m/s) |
| `X_3` | `TEMP`, temperatura (°C) |
| `X_4` | `PRES − 1000`, pressão centrada (hPa) |
| `U_1` | `hour`, hora do dia (0 a 23) |
| `U_2` | umidade relativa, calculada de `TEMP` e `DEWP` pela forma de Magnus |

A leitura substantiva é o efeito de dispersão do vento sobre a concentração
de material particulado, modulado pela hora (camada limite noturna contra
diurna) e pela umidade (higroscopia do aerossol).

**Ressalvas:** a série é horária, de modo que os resíduos são
autocorrelacionados e a partição treino/teste aleatória usada aqui é
otimista para todos os métodos igualmente; a umidade relativa é derivada,
não medida; e apenas um dos doze sítios entra, para que `n` fique na casa
da dezena de milhar sem empilhar sítios com níveis diferentes.

### C. California block groups, censo de 1990

| Item | |
|---|---|
| Fonte | StatLib (Carnegie Mellon), `datasets/houses.zip`, arquivo `cadata.txt` |
| URL | `http://lib.stat.cmu.edu/datasets/houses.zip` |
| SHA-256 do `.zip` | `8b18f0a01cf9c99a65174d18fa582aa31971dfe55a26ad794f3299937c3708d7` |
| Licença | **sem licença explícita.** Arquivo público do StatLib, redistribuído por pacotes de uso geral. É o ponto fraco da candidata diante da exigência editorial de reprodutibilidade (`ss-instrucoes-autores.md`) |
| `n` | 20 640 grupos de quarteirão |
| Citação | Pace & Barry (1997), conferida (§7) |

Mapeamento:

| Papel | Variável |
|---|---|
| `Y` | `log(median house value)` |
| `X_1` | constante |
| `X_2` | `median income` |
| `X_3` | `log(total rooms / households)`, cômodos por domicílio |
| `X_4` | `housing median age / 10` |
| `U_1` | latitude |
| `U_2` | longitude |

A leitura substantiva é o efeito da renda e do tamanho do imóvel sobre o
preço, modulado pela posição geográfica; a quina esperada está em fronteiras
de mercado (litoral, limite metropolitano).

**Ressalvas, ambas substantivas:** a resposta é censurada à direita em
500 001 dólares, e a **aditividade em latitude e longitude é hipótese forte**
sobre um campo espacial. Um efeito que varie ao longo da diagonal da costa
da Califórnia não é aditivo nas duas coordenadas, e um modelo aditivo o
representa mal, não porque falte regularidade mas porque falta o termo de
interação. Se esta candidata for a escolhida, a discussão da hipótese
aditiva tem de estar na seção da aplicação, e não numa nota.

---

## 3. Como a sondagem mediu

`wafc/scripts/05-sondagem-aplicacao.R`, uma partição 70/30 fixada pela
semente, `cv.wafc()` com os padrões atuais (grade `J = 2:⌈log₂ n / 2⌉`,
10 dobras, LASSO), e o `mgcv::gam` de `wafc_competitor("gam", ...)` com
`s(u_m, by = x_ℓ)`, `select = TRUE` e REML. Uma tabela de avaliação da base
(`WaveBased::wtable()`, Daublets 8) é construída uma vez por candidata e
passada explicitamente, que é D31.

O `gam` roda em **três tamanhos de base**, e a razão é de honestidade da
comparação: uma vitória do WAFC sobre um spline de dez nós seria vitória
sobre o tamanho da base do competidor, não sobre a classe de suavidade dele,
e é exatamente essa a confusão que a §4 de `plano-projeto.md` (E4.2) manda
evitar.

| tamanho | qual | como é ajustado |
|---|---|---|
| `k = 10` | o padrão de `wafc_fit_gam()` | `mgcv::gam`, REML, `select = TRUE` |
| `k` grande | o maior que os dados admitem, até 30 | idem |
| **dimensão casada** | `k_m = min(2^J, valores distintos de U_m − 1)`, com `J` o nível que a validação cruzada escolheu | `mgcv::bam`, `fREML`, `discrete = TRUE`, `select = TRUE` |

A terceira linha é a que decide, e duas coisas sobre ela precisam estar
ditas. Primeiro, `wafc_fit_gam()` recebe **um** `k` para todos os
suavizadores, o que trava cada moduladora na mais grosseira delas (com `hr`
em 24 valores, nenhum suavizador daquele ajuste passa de `k = 23`); a
dimensão casada é portanto ajustada por código próprio do script, com `k`
por moduladora. Segundo, ela usa `bam` e não `gam`, porque nestas dimensões
o `gam` não termina em tempo útil (§4.3). É troca de **algoritmo de ajuste**
mais uma discretização das covariáveis em caixas, não troca de estimador, e
está rotulada onde quer que o número apareça. Que a aproximação é benigna se
vê em beijing, onde o `bam` com `k = (16, 16)` dá 0.90070 e o `gam` com
`k = 23` dá 0.90090: a quarta casa decimal.

Três leituras, na ordem de peso do §1:

1. **Energia por nível.** Para cada bloco `(ℓ,m)`, `‖θ̂_{ℓm,j·}‖_2` nível a
   nível. Na convenção `B^{s'}_{2,2}` essa norma é de ordem `2^{-j s'}`, de
   modo que a inclinação de `log₂‖θ̂_{j·}‖_2` contra `j` estima `−s'`.
   Reportam-se `fine.share`, a fração da energia nos dois níveis mais finos,
   e `ŝ'`. O LASSO encolhe, e encolhe mais os níveis finos, o que enviesa
   `ŝ'` **para cima**, ou seja, a favor da leitura "suave" que derrubaria o
   argumento; por isso a leitura principal é o reajuste de mínimos quadrados
   sobre o suporte selecionado, e a leitura encolhida entra ao lado.
2. **Índice de localização.** Fração da variação total de `ĝ_{ℓm}` carregada
   pelos 5% maiores incrementos na grade de 512 pontos. A escala vem das
   próprias formas de `wafc_component()`, que é o vocabulário de D27:

   | forma | `sine` | `cubic` | `heavisine` | `bumps` | `blocks` |
   |---|---|---|---|---|---|
   | índice | 0.080 | 0.135 | 0.177 | 0.583 | 1.000 |

   O índice sobe com ruído de estimação, e não só com quina verdadeira. Por
   isso o valor da componente do `mgcv` na **mesma** base é o piso interno:
   ela é suave por construção, e a comparação que informa é WAFC contra
   `mgcv`, não WAFC contra a tabela acima.
3. **Predição fora da amostra**, RMSE no terço retido, com o WAFC em
   `lambda.min` e `lambda.1se`.

---

## 4. Resultados

### 4.1 Predição fora da amostra

RMSE no terço retido; `R²` é `1 − (RMSE/sd(Y))²`, para dar escala. A coluna
que decide é a última.

| base | `n` | `J` escolhido | `sd(Y)` | WAFC `λ.min` | WAFC `λ.1se` | `gam` `k = 10` | `gam` `k` grande | `gam` dimensão casada |
|---|---|---|---|---|---|---|---|---|
| bike | 17 379 | **7, topo da grade** | 1.486 | 0.6219 | 0.6216 | 0.6480 | 0.6307 (`k = 23`) | **0.6209** (`k = 23, 128`) |
| beijing | 34 287 | 4, interior | 1.160 | 0.9010 | 0.9051 | 0.9009 | 0.9009 (`k = 23`) | **0.9007** (`k = 16, 16`) |
| housing | 20 640 | **7, topo da grade** | 0.569 | 0.2736 | 0.2760 | 0.3151 | 0.2916 (`k = 30`) | **0.2713** (`k = 128, 128`) |

Em `R²`: bike 0.825 (WAFC) contra 0.826 (spline casado); beijing 0.397
contra 0.397; housing 0.769 contra 0.773.

O resultado que importa está na diferença entre as três últimas colunas. O
WAFC **bate** o `mgcv` nos dois tamanhos de base convencionais, e a margem é
grande em housing (RMSE 15.2% menor que `k = 10`, 6.6% menor que `k = 30`).
**Quando o spline recebe a mesma dimensão que a peneira do WAFC, a vantagem
desaparece e inverte**, nas três bases: 0.1% em bike, 0.03% em beijing, 0.9%
em housing, sempre a favor do spline. São diferenças pequenas demais para
serem sinal numa única partição, e é exatamente esse o ponto: **empate**.

O ganho aparente sobre `k = 10` e `k = 30` era ganho de **dimensão**, não de
base.

### 4.2 Estrutura das componentes

Mediana sobre os oito blocos, com a faixa entre parênteses.

| base | localização WAFC | localização `gam` casado | `fine.share` | `ŝ'` (reajustado) | `ŝ'` (encolhido) |
|---|---|---|---|---|---|
| bike | 0.27 (0.21 a 0.35) | 0.14 (0.07 a 0.23) | 0.49 (0.11 a 0.86) | −0.12 (−0.63 a 0.37) | −0.02 |
| beijing | 0.15 (0.12 a 0.24) | 0.11 (0.08 a 0.20) | 0.26 (0.03 a 0.41) | 0.37 (−0.06 a 1.04) | 0.63 |
| housing | 0.24 (0.17 a 0.35) | 0.14 (0.10 a 0.16) | 0.25 (0.08 a 0.37) | 0.06 (−0.25 a 0.33) | 0.21 |

O índice de localização do WAFC é sistematicamente o dobro do índice do
spline de mesma dimensão, nas três bases. **Isso não é evidência a favor da
base de wavelets, e sim contra**: se a rugosidade a mais fosse estrutura
verdadeira, ela apareceria também na predição, e não aparece. Duas
componentes estimadas com o mesmo número de colunas, uma rugosa e outra
suave, com o mesmo erro fora da amostra, dizem que a rugosa está ajustando
ruído.

A leitura de `ŝ'` aponta no mesmo sentido, por um caminho diferente. Em bike
e housing a energia por nível **não decai** (`ŝ'` mediana de −0.12 e 0.06):
um expoente de Besov em torno de zero, ou negativo, não descreve função
alguma com regularidade positiva. Descreve ruído branco espalhado pelos
níveis finos. É o que se espera quando a verdade é suave e o LASSO gasta
coeficientes finos em flutuação amostral. Em beijing a decaimento existe
(`ŝ'` mediana 0.37, e 0.63 na leitura encolhida), mas a energia total dos
blocos de `temp` e `pres` é de ordem `10⁻³`, isto é, esses efeitos
praticamente não se modulam.

As figuras confirmam à vista (`wafc/cache/05-bike.png`,
`05-housing.png`, `05-beijing.png`; a curva do spline nelas é a de `k`
grande, não a de dimensão casada, porque o estágio casado não redesenha):
a componente do WAFC é a curva do spline mais oscilação de alta frequência,
sem quina que o spline tenha deixado passar. Em bike,
`g[hum, day]` e `g[wind, day]` são quase inteiramente oscilação, com norma
de bloco alta (0.73 e 0.85) e nenhuma contrapartida no spline. A única
estrutura que se parece com salto verdadeiro está em `g[one, lon]` de
housing, com um degrau perto de `lon = −122.5` (área da baía) e um pico em
`−118.5` (Los Angeles); é uma componente entre dezesseis.

### 4.3 A grade de `J` foi binding, e o que acontece acima dela

Em bike e em housing a validação cruzada escolheu o maior `J` que a grade
`2:⌈log₂ n / 2⌉` de `cv.wafc()` oferece, com a curva ainda caindo (bike,
`mse` 0.3893 em `J = 6` contra 0.3807 em `J = 7`; housing, 0.0793 contra
0.0738). O teto é a regra herdada de `cv.wall()`, e em `n` da casa de 10⁴
ela para antes do que os dados pedem. A comparação do §4.1 é portanto de um
WAFC que não foi deixado descer tão fundo quanto queria, e isso tem de ser
medido antes de qualquer veredito.

Medido: com a grade estendida em um nível, `J = 8` e cinco dobras, e o
spline acompanhando a subida para `k = 2^8`.

| base | WAFC `J = 7` | WAFC `J = 8` | `gam` casado a `J = 7` | `gam` casado a `J = 8` |
|---|---|---|---|---|
| bike | 0.6219 | **0.6137** | 0.6209 (`k = 23, 128`) | 0.6185 (`k = 23, 256`) |
| housing | 0.2736 | 0.2704 | 0.2713 (`k = 128, 128`) | **0.2677** (`k = 256, 256`) |

Em housing nada muda: os dois melhoram e o spline continua à frente por 1.0%.
Em bike o WAFC **passa à frente** por 0.8%. É o único sinal positivo da
sondagem inteira, e o §4.4 o examina, porque ele não sobrevive ao exame.

### 4.4 O sinal de bike não é adaptação, é vazamento

O ganho de bike em `J = 8` está inteiro nos blocos modulados pelo índice do
dia: `fine.share` de 0.47, 0.74, 0.71 e 0.75 em `g[one, day]`,
`g[temp, day]`, `g[hum, day]` e `g[wind, day]`, com `ŝ'` de −0.16 a −0.64,
isto é, **energia crescendo com o nível**. Um componente com 255 funções de
base sobre 731 dias resolve cerca de três dias; sob partição aleatória de
uma série horária, as horas de teste estão nos mesmos dias das horas de
treino, e um componente fino no índice do dia carrega o nível daquele dia de
uma para a outra. Isso é memorização, não adaptação.

O teste que separa as duas coisas é reter **semanas inteiras**: as horas de
teste nunca dividem um dia com as de treino, a faixa do índice do dia
continua coberta, e nada é extrapolado. Com 105 semanas, 32 retidas, e as
dobras também bloqueadas por semana:

| `J` | 2 | 3 | **4** | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| `mse` da validação cruzada | 0.7770 | 0.5052 | **0.4120** | 0.4187 | 0.4239 | 0.4284 | 0.4281 |

**A curva passa a ter mínimo interior em `J = 4`**, e os níveis finos que a
partição aleatória comprava deixam de pagar. No conjunto retido:

| método | RMSE | relativo |
|---|---|---|
| WAFC `λ.min` (`J = 4`) | **0.6469** | 1.000 |
| `gam` casado a `J = 4` (`k = 16, 16`) | 0.6493 | 1.004 |
| `gam` casado a `J = 8` (`k = 23, 256`) | 0.6559 | 1.014 |
| WAFC `λ.1se` | 0.6635 | 1.026 |

O WAFC fica à frente por 0.4%, que numa partição só não é sinal. E o spline
de dimensão grande (`k = 23, 256`) fica **atrás** do de `k = 16`: sob
partição honesta, resolução alta no índice do dia prejudica os dois métodos.
Ou seja, o ganho de 0.8% do §4.3 era do desenho do experimento, não do
estimador.

O leitor tira daqui duas coisas. Sobre bike: a modulação por hora é suave e
a por dia não tem estrutura fina real; empate. Sobre método: **qualquer
aplicação do WAFC a série temporal tem de usar partição por bloco**, ou a
tabela mede vazamento.

### 4.5 O `mgcv::gam` não é sintonizável nestes `n`, e o `mgcv::bam` é

O `gam`
com oito suavizadores de `k = 23` levou 1453 s em bike e 1875 s em beijing;
o `bam` com `fREML` e covariáveis discretizadas ajustou `k = (23, 128)` em
**6.1 s** e `k = (128, 128)` em **16.3 s**. Duas ordens de grandeza. Um
estudo que sintonize o spline com `gam` vai, por custo, dar a ele uma base
pequena demais, e vai medir a diferença errada.

---

## 5. Veredito por base

O critério do §1 tem três partes e a primeira é eliminatória.

### bike: **neutra**, depois de um falso positivo

É a candidata que deu mais trabalho, porque passou por três leituras
diferentes. Na grade padrão, empata com o spline de dimensão casada
(0.6216 contra 0.6209). Um nível acima do teto da grade, passa à frente por
0.8% (0.6137 contra 0.6185), o único sinal positivo da sondagem. Com semanas
inteiras retidas, o sinal evapora: a validação cruzada volta a escolher
`J = 4`, os níveis finos param de pagar, e o WAFC fica 0.4% à frente do
spline de mesma dimensão, o que numa partição só não é sinal.

A substância é que a modulação por hora é suave: o perfil diário de
`g[one, hour]` é exatamente a forma que um spline penalizado reproduz, e a
quina dos horários de pico, que era a aposta, existe como subida e descida
de curvatura moderada, não como quina. E a modulação pelo índice do dia não
tem estrutura fina real: o que parecia estrutura era o nível do dia passando
do treino para o teste. Some-se o limite estrutural de `hr` ter 24 valores,
onde não há onde exibir regularidade espacialmente heterogênea.

### beijing: **neutra, e fraca por outras razões**

Aqui a base é honesta sobre si mesma: a validação cruzada escolheu `J = 4`
no interior da grade, o decaimento por nível existe, e os três métodos
empatam na terceira casa decimal. Não há o que o WAFC mostre que o spline
não mostre. Além disso, `R² ≈ 0.40` com `sd(Y) = 1.16`: quatro covariáveis
meteorológicas não explicam PM2.5 horário, e a modulação dos efeitos de
temperatura e pressão tem norma de bloco de ordem `10⁻³`, isto é, é nula
para fins práticos. Uma aplicação sobre esta base seria uma aplicação em que
o modelo aditivo de coeficientes não tem muito o que dizer.

### housing: **favorece o spline, mas é a menos ruim**

É a candidata com o maior ganho aparente sobre o `mgcv` convencional (6.6%
sobre `k = 30`, 15.2% sobre `k = 10`) e a única com estrutura que se parece
com salto verdadeiro, o degrau de `g[one, lon]` na longitude da área da
baía. Mas o ganho é de dimensão: com `k = 128` o spline passa à frente
(0.2713 contra 0.2736), e um nível acima o quadro não muda (0.2677 contra
0.2704). E a hipótese aditiva em latitude e longitude, já
registrada no §2 como forte, é provavelmente a razão pela qual nenhum dos
dois modelos vai longe: um campo espacial de preços não é aditivo nas
coordenadas, e os dois métodos pagam o mesmo preço por supor que é.

### Veredito da tarefa

**Nenhuma das três sustenta o argumento do artigo.** Não há aqui uma base em
que o WAFC ganhe do spline por adaptação a regularidade heterogênea; há três
bases em que ele empata quando o spline é sintonizado com honestidade, e
ganha quando não é. A maior diferença a favor do WAFC que sobreviveu a um
desenho de comparação honesto foi de 0.4%, numa partição só.

O resultado não é um fracasso da sondagem: é a resposta que ela existia para
dar, e chega antes de E6 gastar semanas numa aplicação que o referee derruba
com um `k` maior. O que ele diz sobre a próxima rodada de busca está no §3
de `handoff-E6.1a.md`: o critério que falta às três é **quina verdadeira na
moduladora**, e as três foram escolhidas por terem moduladora com
interpretação forte (hora, coordenada), não por haver evidência prévia de
descontinuidade. As famílias que têm essa evidência são outras (regressão
em descontinuidade, com limiar administrativo, faixa de imposto ou nota de
corte; limiar regulatório em dose-resposta; séries com quebra datada), e são
elas que a próxima sondagem deve levantar.

---

## 6. O que esta sondagem não cobre

- **Uma partição, uma semente.** As diferenças de RMSE não vêm com erro
  padrão de reamostragem. Diferenças de poucos milésimos não são sinal, e
  isso inclui a de 0.4% que sobrou a favor do WAFC em bike sob partição
  por blocos.
- **Dependência serial, medida em bike e não em beijing.** Duas das três
  candidatas são séries horárias, e a partição aleatória trata vizinhos
  temporais como independentes. Em bike isso foi medido e importou muito
  (§4.4). Em **beijing não foi medido**: os números dela são de partição
  aleatória, e como os três métodos empatam ali de qualquer forma, o que a
  partição por bloco mudaria é o valor absoluto do erro, não a ordem. Em
  housing o análogo é a dependência espacial, e o bloqueio correspondente
  seria por região; **também não foi feito**, e a comparação do §4.1 para
  housing herda essa ressalva.
- **Sintonia do competidor.** O `mgcv` roda em REML com `select = TRUE` em
  três dimensões, a maior delas por `bam`. Não entraram aqui o spline adaptativo de Wang,
  Jiang & Liu (2024) nem o B-spline com group LASSO, que são de E2.4 e de
  E4; se a base escolhida sobreviver a eles é pergunta de E6.2.
- **`p` e `q` fixos.** Quatro covariáveis lineares e duas moduladoras em
  todas as candidatas, para que as três sejam comparáveis entre si. A
  aplicação final pode querer outro conjunto.
- **Nada aqui decide a aplicação.** O veredito do §5 diz quais bases
  sustentam o argumento; escolher entre as que sustentam é decisão do autor,
  e envolve interpretação e licença tanto quanto número.

---

## 7. Entradas de bibliografia, conferidas

Conferidas no Crossref em 2026-09-20 (título, autores, ano, veículo, volume,
número, páginas, DOI), no padrão de L1, L3 e L4. **Não foram escritas em
`docs/referencias-verificadas.bib`**, que não é arquivo desta tarefa; entram
quando E6.1 fechar e a base for escolhida.

```bibtex
% Crossref registra o fascículo impresso em 2014 (vol. 2, n. 2-3) e a
% publicacao eletronica em 2013-11-26. O ano de citacao e 2014.
@article{Fanaee-T-Gama-2014,
  author  = {Fanaee-T, Hadi and Gama, Jo\~{a}o},
  title   = {Event labeling combining ensemble detectors and background knowledge},
  journal = {Progress in Artificial Intelligence},
  year    = {2014},
  volume  = {2},
  number  = {2-3},
  pages   = {113--127},
  doi     = {10.1007/s13748-013-0040-3}
}

@article{Zhang-Guo-Dong-He-Xu-Chen-2017,
  author  = {Zhang, Shuyi and Guo, Bin and Dong, Anlan and He, Jing and
             Xu, Ziping and Chen, Song Xi},
  title   = {Cautionary tales on air-quality improvement in {Beijing}},
  journal = {Proceedings of the Royal Society A: Mathematical, Physical and
             Engineering Sciences},
  year    = {2017},
  volume  = {473},
  number  = {2205},
  pages   = {20170457},
  doi     = {10.1098/rspa.2017.0457}
}

% O Crossref grafa o primeiro autor como "Kelley Pace, R."; a forma usada
% pelo proprio arquivo do StatLib e pela literatura e "Pace, R. Kelley".
@article{Pace-Barry-1997,
  author  = {Pace, R. Kelley and Barry, Ronald},
  title   = {Sparse spatial autoregressions},
  journal = {Statistics \& Probability Letters},
  year    = {1997},
  volume  = {33},
  number  = {3},
  pages   = {291--297},
  doi     = {10.1016/S0167-7152(96)00140-X}
}
```

Os dois conjuntos do UCI também têm DOI próprio de conjunto de dados
(`10.24432/C5W894` e `10.24432/C5RK5G`), que é o que se cita ao lado do
artigo quando a revista exige a fonte dos dados e não só o trabalho que os
descreve.

---

## 8. Reprodução

Da raiz do repositório:

```
Rscript wafc/scripts/05-sondagem-aplicacao.R all
```

Sem argumentos, roda as três candidatas sem subamostragem, com a semente
`20260920`. Um candidato isolado: `... 05-sondagem-aplicacao.R bike`. Um
teto de amostra, para inspeção rápida: `... 05-sondagem-aplicacao.R bike 3000`.

Os argumentos são `[candidatas] [n_max] [semente] [estágio]`, e os estágios
são quatro:

| estágio | o que faz | custo |
|---|---|---|
| `full` (padrão) | a sondagem inteira do §4.1 e §4.2 | horas (bike 2265 s, beijing 6627 s, housing 12 329 s só na validação cruzada) |
| `matched` | acrescenta o spline de dimensão casada a uma corrida pronta | segundos |
| `deep` | o nível acima do teto da grade, mais o spline que o acompanha (§4.3) | bike 528 s, housing 3713 s |
| `blocked` | a partição por semanas inteiras, só para bike (§4.4) | 950 s |

Os três últimos leem `05-sondagem.rds` e precisam de uma corrida `full`
antes. Tudo é guardado em `wafc/cache/05-cache-*.rds` por unidade, de modo
que repetir um estágio não repete o que já custou; apagar esses arquivos
força o recálculo.

O script baixa os arquivos na primeira execução, confere as somas SHA-256
acima e avisa (sem parar) se algum arquivo mudou. Saídas em `wafc/cache/`:
`05-sondagem.rds` com todos os números, um `.log` por estágio e
`05-<candidata>.png` com as componentes estimadas. **Atenção ao ler as
figuras:** a curva do spline nelas é a de `k` grande (`k = 23` ou `k = 30`),
porque só o estágio `full` desenha, e é justamente essa a dimensão que o
§4.1 mostra ser pequena demais. A figura serve para ver a forma da
componente do WAFC, não para julgar a comparação.

---

## 9. E6.1b: o critério de D44 e como a sondagem mediu

**Atualização de 2026-10-05.** As seções 1 a 8 acima são de E6.1a e ficam
como estão. Daqui em diante, os números são de
`wafc/scripts/10-sondagem-aplicacao-b.R`, com o código do retrato
`b0ea096` (`wafc/cache/e61b/snap/`, conferido arquivo a arquivo contra o
commit a cada início).

**O critério mudou com D44.** A tese deixou de ser "o WAFC ganha do spline
por adaptação" e passou a ser **estrutura recuperada com predição
competitiva**. Uma base passa se: (1) a estrutura que o WAFC devolve é
estável entre partições e interpretável; e (2) a predição fica a menos de um
erro-padrão do melhor `gam`.

**O estimador é o de D44 e D45**, com todos os argumentos escritos na
chamada:
- o ajuste é o block LASSO balanceado (`wafc_fit_klopp` com níveis livres e
  pesos do `grpreg`), com `J` em `2:8` e 10 dobras;
- em seguida vem o limiar nas normas por bloco, com a regra `cv1se` (o padrão
  de D45) e também com a `cv`;
- as dobras do limiar são ajustadas uma vez para as duas regras.

**Os concorrentes:**
- os dois splines de D46, `gam.reml` e `gam.gcv` (`bam`, grade
  `5, 10, 20, 40, 80` de D41);
- o linear;
- depois da primeira leitura, um terceiro spline, `gam.cv` (§13): `k`
  escolhido na mesma grade pelo erro de validação cruzada nas mesmas dobras
  por bloco do WAFC, com REML dentro de cada dobra (`k.select = "cv"` de
  E2.5h).

**A partição é por bloco** onde há dependência, a lição de §4.4:
- os blocos são a semana nas séries e o ladrilho de 0,25° em housing;
- 30% dos blocos vão ao teste;
- as dobras internas também são formadas por blocos inteiros;
- são **20 partições**, com semente por base e partição.

O erro-padrão de cada método é o desvio-padrão entre partições dividido por
`sqrt(20)`. A comparação com o melhor spline é pareada (a diferença de RMSE
na mesma partição), e vem com dois erros-padrão:
- o ingênuo, `sd/sqrt(20)`;
- o corrigido para treinos sobrepostos, pelo fator de Nadeau & Bengio
  (2003), `sqrt((1/K + n_teste/n_treino)/(1/K))`, que dá ~3,1 aqui.

O critério (2) usa o ingênuo. O corrigido mostra o quanto a conclusão depende
dessa escolha.

**A estabilidade da estrutura** é medida assim:
- o treino de cada partição é uma subamostra de 70% dos blocos;
- para cada bloco `(ℓ, m)`, conta-se a fração das 20 em que ele fica no ajuste
  limiarizado;
- ao lado, a mesma fração para os suavizadores do `gam` com `edf > 0,1` (o
  `select = TRUE` encolhe, mas nunca zera exatamente), e o `edf` mediano.

**O custo:**
- 120 unidades (6 bases × 20 partições);
- lançada em 2026-10-03 às 22h16 com 4 processos;
- parada às 07h31 de 2026-10-04, quando a swap passou de 1 GB (duas
  conferências da E1.14 começaram às 07h24), com 4 unidades perdidas;
- relançada às 07h33 com 3 processos, retomando do que estava gravado;
- fim em 2026-10-05 às 13h48: 30 h 15 min de relógio desde o relançamento, 39 h 32 min desde o primeiro lançamento.

Mediana por unidade:

| base | WAFC (ajuste, as duas regras de limiar) | `gam.reml` | `gam.gcv` | `gam.cv` |
|---|---|---|---|---|
| bike | 1 691 s | 6 s | 85 s | 40 s |
| beijing | 3 322 s | 8 s | 134 s | 42 s |
| beijing.heat | 7 750 s | 18 s | 392 s | 108 s |
| housing | 5 009 s | 12 s | 237 s | 77 s |
| marylebone | 2 530 s | 7 s | 108 s | 25 s |
| kelmarsh | 598 s | 6 s | 128 s | 35 s |

O pico foi de 4,9 GB num processo e 16 GB somados; a swap chegou a 2,2 GB na fase de 3 processos, partindo do resíduo de 1,57 GB deixado pela parada.

---

## 10. As candidatas novas (frente ii)

O critério de busca é o do §5: salto ou limiar documentado na literatura da
área, dado público e citável, licença declarada, `n` na casa dos milhares,
duas ou mais moduladoras. Três passaram e foram sondadas.

| base | dado e licença | o salto documentado | mapeamento |
|---|---|---|---|
| `beijing.heat` | o Dongsi da §2.B, CC BY 4.0, DOI 10.24432/C5RK5G; mesmo arquivo, sem download | a temporada de aquecimento de Pequim, de 15 de novembro a 15 de março (Liang et al. 2015, Proc. R. Soc. A 471: 20150257) | `Y = log PM2.5`; `X = (1, vento, temp, pressão − 1000)`; `U = (dia do ano, umidade relativa)`. O dia do ano é periódico, como a base |
| `marylebone` | `mydata` do `openair` 3.1.0 (Carslaw & Ropkins 2012): 65 533 horas, de 1998-01-01 a 2005-06-23, 60 780 completas | a fração primária NO2/NOx do tráfego de Londres subiu de ~5–6% (1997) a ~17% (2003), mudança ligada aos filtros dos ônibus (Carslaw 2005) | `Y = NO2 + O3` (o oxidante); `X = (1, NOx/100)`; `U = (data em dias, velocidade do vento)`. É a relação de Clapp & Jenkin (2001): `β_NOx(u)` é a fração primária |
| `kelmarsh` | SCADA de 10 min de 2017 da turbina 1 (Senvion MM92, 2 050 kW); Plumley (2022), Zenodo 10.5281/zenodo.5841834, CC BY 4.0; 51 185 intervalos sem parada nem corte, de 52 560 | a velocidade nominal da curva de potência: a densidade do ar aumenta a potência só abaixo dela (normalização da IEC 61400-12-1; Lee, Ding, Genton & Xie 2015, JASA) | `Y = potência (MW)`; `X = (1, temperatura/10)`; `U = (velocidade, direção do vento)`. `β_temp(u)` deve ser negativo abaixo da nominal e voltar ao nível de fora da faixa acima dela |

Detalhes de cada base:
- **Licença do `mydata`.** O pacote é MIT. As medidas de poluição vêm do
  London Air Quality Archive, e a London Air Quality Network declara a Open
  Government Licence v2 na página da sua API. Ficam dois `[VERIFICAR]`: a
  página não diz se a licença cobre o arquivo histórico, e o `mydata.Rd` não
  diz de onde vêm o vento e a direção. A saída limpa, se a base for escolhida,
  é a série de Marylebone Road (MY1) do UK-AIR (Defra, OGL v3), com vento de
  fonte declarada.
- **Direção do vento em marylebone.** Não é moduladora: tem 38 valores, em
  passos de 10°, o caso discreto que o §5 exclui.
- **Ano de Kelmarsh: 2017.** É um ano civil inteiro, com 97% dos intervalos
  passando no filtro, e o menor zip depois do primeiro ano (174,6 MB).
- **Turbulência em kelmarsh.** O desvio-padrão do vento, que daria a
  intensidade de turbulência de Lee et al., falta em 73% dos intervalos e
  ficou fora.

---

## 11. Predição

RMSE no teste, média sobre as 20 partições, com o erro-padrão entre
partições. "Razão" é a média, por partição, da razão ao melhor spline da
mesma partição (o melhor de `gam.reml`, `gam.gcv` e `gam.cv`). "Diferença" é
a média da diferença pareada ao spline de menor RMSE médio, com o erro-padrão
ingênuo e o corrigido. "Vence" é a fração de partições em que o método bate
esse spline.

| base | método | RMSE | razão | diferença ao melhor `gam` | vence |
|---|---|---|---|---|---|
| bike | WAFC `+cv1se` | 0,6680 ± 0,0057 | 1,036 | 0,0225 ± 0,0037 (0,0114) | 1/20 |
| | WAFC `+cv` | 0,6642 ± 0,0092 | 1,030 | 0,0187 ± 0,0067 (0,0208) | 1/20 |
| | `gam.cv` (melhor) | 0,6454 ± 0,0041 | | | |
| | `gam.reml` / `gam.gcv` / linear | 0,6563 / 0,7437 / 1,2939 | | | |
| beijing | WAFC `+cv1se` | 0,9211 ± 0,0050 | 1,010 | 0,0087 ± 0,0011 (0,0034) | 0/20 |
| | WAFC `+cv` | 0,9174 ± 0,0051 | 1,006 | 0,0050 ± 0,0005 (0,0016) | 0/20 |
| | `gam.cv` (melhor) | 0,9124 ± 0,0050 | | | |
| | `gam.reml` / `gam.gcv` / linear | 0,9125 / 0,9131 / 1,0543 | | | |
| beijing.heat | WAFC `+cv1se` | 0,7945 ± 0,0041 | 1,026 | 0,0199 ± 0,0027 (0,0085) | 0/20 |
| | WAFC `+cv` | 0,7902 ± 0,0037 | 1,020 | 0,0156 ± 0,0019 (0,0059) | 0/20 |
| | `gam.cv` (melhor) | 0,7746 ± 0,0042 | | | |
| | `gam.reml` / `gam.gcv` / linear | 0,8666 / 0,8844 / 1,0505 | | | |
| housing | WAFC `+cv1se` | 0,3600 ± 0,0072 | 1,066 | 0,0176 ± 0,0051 (0,0155) | 4/20 |
| | WAFC `+cv` | 0,3517 ± 0,0075 | 1,041 | 0,0093 ± 0,0049 (0,0150) | 6/20 |
| | `gam.cv` (melhor) | 0,3424 ± 0,0060 | | | |
| | `gam.reml` / `gam.gcv` / linear | 0,3483 / 0,3705 / 0,4236 | | | |
| marylebone | WAFC `+cv1se` | 12,02 ± 0,10 | 1,040 | 0,444 ± 0,062 (0,192) | 1/20 |
| | WAFC `+cv` | 11,79 ± 0,10 | 1,019 | 0,212 ± 0,040 (0,125) | 1/20 |
| | `gam.cv` (melhor) | 11,58 ± 0,10 | | | |
| | `gam.reml` / `gam.gcv` / linear | 11,77 / 11,80 / 15,69 | | | |
| kelmarsh | WAFC `+cv1se` | 0,06129 ± 0,00076 | 1,072 | 0,0041 ± 0,0008 (0,0024) | 0/20 |
| | WAFC `+cv` | 0,05901 ± 0,00075 | 1,033 | 0,0018 ± 0,0008 (0,0026) | 0/20 |
| | `gam.cv` / `gam.reml` (melhores, iguais) | 0,05723 ± 0,00053 | | | |
| | `gam.gcv` / linear | 0,05724 / 0,5988 | | | |

O `R²` no teste fica em ~0,80 em bike, ~0,39 em beijing, ~0,55 em beijing.heat, ~0,63 em
housing, ~0,69 em marylebone e ~0,99 em kelmarsh.

Três leituras:
1. **Nenhuma base passa a perna da predição contra o melhor spline.** O WAFC
   fica de 0,6% a 7% atrás, em mais de um erro-padrão ingênuo em toda base.
   Com o erro-padrão corrigido, só o `+cv` de bike, housing e kelmarsh fica
   dentro.
2. **Contra os dois splines de D46, o quadro seria outro.**
   - Em beijing.heat o WAFC venceria por 8% em 20 de 20 partições (razão
     0,917 contra o `gam.reml`).
   - Em marylebone e housing o `+cv` empataria (diferença de 0,018 ± 0,043 e
     0,0034 ± 0,0074 contra o `gam.reml`).
   - A diferença vem toda do `gam.cv`; a §13 diz por quê.
3. **O `cv1se` custa predição em toda base**, de 0,4% (beijing) a 3,9%
   (kelmarsh) sobre o `+cv`. Ele zera blocos de efeito pequeno, mas real, que
   o spline e o `+cv` usam: a direção do vento em kelmarsh, a velocidade do
   vento em marylebone.

---

## 12. Estrutura, base por base

Fração das 20 partições em que o bloco fica. Para o WAFC, a decisão do
limiar; para os splines, `edf > 0,1`, com o `edf` mediano do `gam.cv` entre
parênteses. "Forma" resume a componente mediana do WAFC `+cv1se`.

**bike** (`U = hora, dia`).

| bloco | `+cv1se` | `+cv` | `gam.reml` | `gam.cv` (`edf`) |
|---|---|---|---|---|
| one × hora | 1 | 1 | 1 | 1 (18) |
| one × dia | 0,95 | 1 | 1 | 1 (12) |
| temp × hora | 0,85 | 1 | 1 | 1 (8,6) |
| temp × dia | 0,70 | 1 | 1 | 1 (15) |
| umidade × hora | 0 | 0,20 | 1 | 1 (4,9) |
| umidade × dia | 0,40 | 0,95 | 1 | 1 (17) |
| vento × hora | 0,55 | 1 | 0,55 | 0,50 (0,1) |
| vento × dia | 0,60 | 0,95 | 0,75 | 0,85 (7,9) |

**A estrutura não é estável:**
- metade dos blocos tem decisão entre 0,40 e 0,85 no `+cv1se`;
- o único zero estável é umidade × hora, que o spline mantém com `edf` 4,9;
- os quatro blocos da hora têm a mesma ressalva de beijing (24 valores,
  forma e norma não identificadas com `J ≥ 5`, que a validação cruzada
  escolheu em 15 de 20 partições);
- o efeito do vento é o único que o spline também trata como fraco
  (`edf` 0,1 na hora).

**beijing** (`U = hora, umidade`).

| bloco | `+cv1se` | `+cv` | `gam.reml` | `gam.cv` (`edf`) |
|---|---|---|---|---|
| one × hora | 1 | 1 | 1 | 1 (8,8) |
| one × umidade | 1 | 1 | 1 | 1 (4,9) |
| vento × hora | 0,95 | 1 | 1 | 1 (4,5) |
| vento × umidade | 1 | 1 | 1 | 1 (7,6) |
| temp × hora | 0,15 | 0,85 | 0,75 | 0,80 (2,6) |
| temp × umidade | 0,05 | 0,65 | 1 | 1 (3,9) |
| pressão × hora | 0,70 | 1 | 1 | 1 (7,2) |
| pressão × umidade | 1 | 1 | 1 | 1 (6,7) |

O WAFC zera, de modo estável, a modulação do efeito da temperatura nas duas
moduladoras: o efeito é constante. Isso confirma §4.2 (normas de ordem
`10⁻³`). Os três blocos da hora têm uma ressalva que os invalida como
leitura de forma:
- a hora tem 24 valores, e com `J ≥ 5` o bloco tem 31 ou mais colunas sobre
  24 pontos;
- a curva entre as horas inteiras e a norma do bloco não são identificadas
  (a norma mediana de `one × hora` é 639, contra 0,29 de `one × umidade`);
- a decisão do limiar nesses blocos depende de direções sem dado
  (pressão × hora em 0,70), e o script não reporta a forma deles.

**beijing.heat** (`U = dia do ano, umidade`). Todo bloco fica em todo
método (pressão × umidade em 0,85 no `+cv1se`). Não há seleção, então a
estabilidade é trivial. As formas são legíveis:
- o nível tem o "U" sazonal, alto no inverno, mas a amplitude dele (~6 na
  escala do log) não é o efeito sazonal: os termos da temperatura e da
  pressão também variam com a estação, e só a soma é;
- o efeito do vento é mais negativo no inverno (−0,21 por m/s em relação ao
  verão): a dispersão pesa mais na temporada de aquecimento;
- o efeito da temperatura é mais negativo no inverno (−0,056 por °C).

**O limiar administrativo não aparece como salto estável.**
- Na componente mediana do WAFC, a maior inclinação de `g[temp, dia do ano]`
  está exatamente no dia 318 (15 de novembro), a primeira de 255 posições.
- Partição a partição, porém, ela cai a até 10 dias de 15 de novembro em 7
  de 20. No `gam.cv` são 5 de 20 (entre os dias 310 e 314), e no `gam.reml`,
  1 de 20.
- Em 15 de março, nada.
- Com `J = 4` (17 de 20 partições), o nível mais fino do WAFC tem suporte de
  quase um ano. Ele não resolve um degrau de dias; a validação cruzada por
  semana preferiu suavidade.

**housing** (`U = latitude, longitude`). Todo bloco fica, menos renda ×
latitude (0,55 no `+cv1se`, 0,95 no `gam.cv`) e renda × longitude (0,90 e
1). Esses dois blocos são os únicos com decisão instável.

A forma das componentes é que não é estável. As faixas entre partições são
largas, porque o ladrilho retido leva consigo trechos inteiros da costa. O
pico de `g[one, longitude]` perto de −118,5 (Los Angeles) aparece; o degrau
da área da baía de E6.1a não se distingue. A ressalva de §2.C sobre
aditividade num campo espacial continua, e a licença também.

**marylebone** (`U = data, velocidade do vento`).

| bloco | `+cv1se` | `+cv` | `gam.reml` | `gam.cv` (`edf`) |
|---|---|---|---|---|
| one × data | 1 | 1 | 1 | 1 (37) |
| one × vento | 0,60 | 1 | 1 | 1 (4,9) |
| NOx × data | 1 | 1 | 1 | 1 (36) |
| NOx × vento | 0 | 1 | 1 | 1 (10) |

**É a estrutura mais forte da sondagem, e coincide com a literatura.**
- `g[NOx, data]`, a fração primária ao longo do tempo, fica plana até 2002 e
  sobe em degrau:
  - no WAFC, ~10,3 por 100 ppb de NOx (entre 9,7 e 11,2 nas 20 partições),
    de −3,2 a +7,2;
  - no `gam.cv`, 11,2 (entre 10,9 e 12,3);
  - os 10% da subida caem entre março e outubro de 2002 e os 90%, entre
    julho e setembro de 2003, em toda partição.
- São ~10 pontos percentuais de fração primária no período que Carslaw
  (2005) dá para a passagem de ~5–6% a ~17%.
- `g[one, data]` é o ciclo anual do oxidante de fundo, sete ciclos.
- O `+cv1se` diz, em 20 de 20 partições, que a fração primária não depende
  da velocidade do vento, e zera o vento no fundo em 8 de 20. O spline mantém
  os dois blocos (`edf` 10 e 4,9).
- Essa estrutura mais parcimoniosa custa 2,0% de predição contra o `+cv`, e
  4% contra o `gam.cv`.
- A queda de `g[NOx, data]` no último mês da série é artefato de borda da
  base periodizada (data não é periódica).

**kelmarsh** (`U = velocidade, direção do vento`).

| bloco | `+cv1se` | `+cv` | `gam.reml` | `gam.cv` (`edf`) |
|---|---|---|---|---|
| one × velocidade | 1 | 1 | 1 | 1 (28) |
| one × direção | 0 | 1 | 1 | 1 (36) |
| temp × velocidade | 1 | 1 | 1 | 1 (16) |
| temp × direção | 0 | 1 | 1 | 1 (31) |

A forma de `g[temp, velocidade]` é a da física, e o WAFC e o `gam.cv` a dão
iguais:
- +0,056 MW por 10 °C abaixo da entrada (3 m/s);
- mínimo de −0,090 em 9,7 m/s;
- +0,045 acima de 12,5 m/s.

O efeito da temperatura é negativo na carga parcial e volta ao nível de fora
da faixa acima da nominal. A amplitude, ~0,14 MW por 10 °C na carga parcial,
é cerca de 2,5 vezes a que a densidade sozinha daria (3,5% de 1,6 MW). A
temperatura carrega também a estabilidade atmosférica sazonal; é leitura,
não medida.

Sobre a direção:
- o `+cv1se` zera os dois blocos da direção em 20 de 20 partições;
- o spline os mantém, com `edf` de 31 a 36 e amplitude pequena (0,17 e
  0,08 MW);
- o `+cv` os mantém e prediz 3,9% melhor que o `+cv1se`.

Então o efeito das esteiras existe e é pequeno, e o `cv1se` o descarta.
A borda da base periodizada aparece em `g[one, velocidade]` abaixo de 2 m/s.

---

## 13. A conferência dos splines: o `gam` sintonizado nas dobras por bloco

Na primeira leitura, beijing.heat dava ao WAFC 8% de vantagem em 20 de 20
partições. Mas o `k` do `gam.reml` estava no topo da grade (80) em toda
partição, com `edf` de 61 a 71 nos suavizadores do dia do ano. Duas leituras
cabiam, e diziam coisas opostas sobre o WAFC:
- (a) a grade de D41 limita o spline, a lição de dimensão de E6.1a;
- (b) o REML e o GCV supõem erros independentes e escolhem pouca suavização
  numa série horária de resíduos correlacionados (Opsomer, Wang & Yang
  2001), enquanto o WAFC escolhe `(J, λ)` por validação cruzada em semanas
  inteiras.

O `gam.cv` separa as duas: ele escolhe `k` na grade de D41 pelo erro nas
dobras por semana do próprio WAFC, ~1 a 2 min por partição.

| base | `k` do `gam.cv` (partições) | `gam.cv` contra `gam.reml` |
|---|---|---|
| beijing.heat | 5 (13), 10 (7); nunca o topo | 0,7746 contra 0,8666: −10,6% |
| marylebone | 20 (9), 40 (11) | 11,58 contra 11,77: −1,6% |
| housing | 10 (7), 20 (12), 80 (1) | 0,3424 contra 0,3483: −1,7% |
| beijing | 10 (10), 20 (6), 23/40 (2), 23/80 (2) | igual (0,9124 contra 0,9125) |
| kelmarsh | 40 (12), 80 (8) | igual |
| bike | 20 (13), 23/40 (7) | 0,6454 contra 0,6563: −1,7% |

**A leitura (b) é a que vale.**
- Com `k` escolhido nas dobras por bloco, o spline de beijing.heat desce a
  `k = 5` ou 10 (`edf` ~4 por suavizador) e passa à frente do WAFC por 2,0%
  a 2,6%.
- O topo da grade não estava limitando: o REML queria mais `k` porque
  suavizava de menos.
- Na fumaça com 3 000 pontos, o REML numa grade até 320 subiu a 320 e
  predisse pior; não foi medido no tamanho cheio.

**A lição vale para E4 e para o manuscrito, além da aplicação.** Em dado
dependente, o spline tem de escolher a dimensão pela mesma validação cruzada
por bloco que o WAFC usa. Senão a tabela mede o critério de sintonia, não o
estimador, e mede a favor do WAFC.

---

## 14. Veredito por base no critério de D44

| base | estrutura estável e interpretável? | predição a menos de 1 EP do melhor `gam`? | licença | veredito |
|---|---|---|---|---|
| bike | não: decisões entre 0,40 e 0,85 em metade dos blocos; os blocos da hora não identificados | não (+3,0% a +3,6%, 1/20) | CC BY 4.0 | **não passa** |
| beijing | em parte: o efeito da temperatura sem modulação, estável; os blocos da hora não identificados | não (+0,6% a +1,0%, 0/20) | CC BY 4.0 | **não passa** |
| beijing.heat | estável só porque nada é zerado; legível (sazonalidade dos efeitos); o limiar de 15 de novembro em 7/20 | não (+2,0% a +2,6% contra o `gam.cv`, 0/20); ganharia 8% contra os splines de D46, por causa do critério de sintonia deles | CC BY 4.0 | **não passa** |
| housing | a decisão é estável fora dos dois blocos da renda, a forma não; aditividade espacial duvidosa | não (+4,1% a +6,6%; o `+cv` fica dentro só pelo EP corrigido) | não declarada | **não passa** |
| marylebone | **sim**: o degrau documentado da fração primária em 2002–2003, estável em 20/20, mais a ausência de modulação pelo vento | não (+1,9% a +4,0%, 1/20); o `+cv` empataria com os splines de D46 | MIT; OGL v2 nas medidas; vento `[VERIFICAR]` | **não passa na predição; é a melhor estrutura** |
| kelmarsh | sim na velocidade (curva de potência; efeito da temperatura que some acima da nominal); a direção, zerada, faz falta na predição | não (+3,3% a +7,2%, 0/20) | CC BY 4.0 | **não passa** |

**Nenhuma das seis passa o critério de D44 contra o spline sintonizado nas
mesmas dobras.**

**Se o critério for contra os dois splines de D46**, a escolha que o
`TAREFA.md` e o D46 fixam:
- marylebone passa (a estrutura documentada mais o `+cv` empatado);
- housing passa na predição e falha na licença;
- beijing.heat venceria.

A sondagem recomenda não usar essa leitura. A vantagem vem do critério de
sintonia do spline, e um referee que rode `k` por validação cruzada em
blocos a derruba.

**A base mais forte para uma aplicação de estrutura é marylebone.** Ela
mostra o que a tese de D44 promete: uma componente com um degrau datado,
documentado na literatura da área, recuperado em toda partição. Também mostra
o custo: o WAFC prediz ~2% pior que um spline que vê o mesmo degrau. A
escolha continua do autor (pergunta 2 do `ESTADO.md`).

---

## 15. O que esta sondagem não cobre

- **O erro-padrão é entre partições, e as partições compartilham ~70% do
  treino.** A correção de Nadeau & Bengio supõe partições aleatórias de
  observações, não de blocos; as duas colunas cercam o valor certo, sem
  dá-lo.
- **O `gam.cv` escolhe `k`, não os parâmetros de suavização**, que continuam
  por REML dentro de cada dobra. Um spline com suavização escolhida por
  validação cruzada por bloco, ou com erro correlacionado (`gamm`, `bam` com
  `rho`), pode ir mais longe.
- **A grade de REML acima de 80 não foi medida no tamanho cheio.**
- **Moduladora discreta.** Com a hora (24 valores) e `J ≥ 5`, o bloco do WAFC
  não é identificado fora das horas inteiras, e a norma que o limiar lê
  também não. Truncar `J` por moduladora no número de valores distintos,
  como `wafc_k_matched()` faz com o `k` do spline, resolveria; é proposta
  para o código, não feita aqui.
- **Borda da base periodizada.** Data, velocidade do vento e coordenadas não
  são periódicas, e as componentes têm artefatos na primeira e na última
  fração da grade. Os índices de localização e os "maiores saltos" nas
  bordas não são leitura.
- **Uma turbina, um sítio, um ano.** Kelmarsh usa a turbina 1 de 2017;
  beijing.heat, o Dongsi; marylebone, um sítio.
- **As datas efetivas da temporada de aquecimento** de cada ano de
  2013–2017 não foram conferidas; as marcas são 15 de novembro e 15 de março
  em todo ano.

---

## 16. Reprodução

Da raiz do repositório:

```
NC=4 setsid nohup wafc/cache/e61b/run.sh > wafc/cache/e61b/run.out 2>&1 &
E61B_EXTRA=gam.cv Rscript wafc/scripts/10-sondagem-aplicacao-b.R all extra 4 20 0
Rscript wafc/scripts/10-sondagem-aplicacao-b.R all report 1 20 0
```

- **Ordem:** o `run.sh` (não versionado) roda as partes `fit` e `report` das
  seis bases e escreve `run end` no `call1.log` ao terminar; a parte `extra`
  acrescenta o `gam.cv`; o `report` final junta as duas.
- **Custo:** ~39 h de relógio em 3 a 4 processos na parte `fit` (30 h sem a parada); ~1 h em 2
  processos na `extra`.
- **Retomada:** uma unidade já gravada em `units/` ou `units-extra/` não é
  refeita.
- **Saídas,** em `wafc/cache/e61b/`, não versionadas: `units/`,
  `units-extra/`, `e61b-summary.rds`, `e61b-<base>.png`, `call1.log`,
  `extra-call.log` e `mem.log`.
- **Dados:** as somas SHA-256 dos arquivos estão no script; os de E6.1a,
  também na §2.
