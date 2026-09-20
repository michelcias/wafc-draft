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
