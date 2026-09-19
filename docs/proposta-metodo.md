# Proposta de método: a ideia inicial

Registro do que foi discutido em 2026-09-18. É um ponto de partida, não uma
especificação: tudo aqui pode ser melhorado, e o que for provado passa para
`derivations/`. Notação em [`notacao.md`](notacao.md), congelada em E1.1
(2026-09-18); os símbolos deste arquivo foram trocados para ela na mesma
rodada, e o conteúdo é o de 2026-09-18.

---

## 1. O modelo

Para uma resposta escalar `Y`, `p` covariáveis "lineares" `X = (X_1, …, X_p)`
e `q` covariáveis "moduladoras" `U = (U_1, …, U_q)`, com `U_m` a valores num
intervalo (reescalado para `[0, 1]`):

```
Y = Σ_{ℓ=1}^{p} β_ℓ(U) X_ℓ + ε,           E[ε | X, U] = 0
β_ℓ(u) = c_ℓ + Σ_{m=1}^{q} g_{ℓm}(u_m),   ∫_0^1 g_{ℓm}(u) du = 0
```

- **Modelo de coeficientes funcionais** (varying-coefficient) no sentido de
  Hastie & Tibshirani (1993) e Cai, Fan & Yao (2000): o efeito de `X_ℓ` é uma
  função de `U`.
- **Estrutura aditiva em cada coeficiente**: cada `β_ℓ` é soma de funções
  univariadas, uma por moduladora. É o *additive coefficient model* de Xue &
  Yang (2006), que escapa da maldição da dimensionalidade em `q`.
- **Casos particulares que o método tem de recuperar**: `X_1 ≡ 1` dá um
  intercepto aditivo `c_1 + Σ_m g_{1m}(U_m)`, logo o modelo aditivo puro
  (`p = 1`) e o modelo parcialmente linear aditivo (`p > 1`, só `β_1`
  funcional) são casos particulares; `q = 1` é o varying-coefficient
  clássico com um índice; `g_{ℓm} ≡ 0` para todo `m` devolve a regressão
  linear.

A restrição `∫_0^1 g_{ℓm} = 0` (centralização de Lebesgue, D22; era
`E[g_{ℓm}(U_m)] = 0` no registro de 2026-09-18) é a de identificabilidade
das componentes
aditivas; `c_ℓ` absorve o nível. Identificar `(c_ℓ, g_{ℓm})` a partir da
distribuição de `(Y, X, U)` exige que `E[X X' | U]` seja não singular
(quase certamente) e que as `U_m` não sejam funções uma da outra; isso vira
hipótese em E1.2.

## 2. A aproximação por wavelets

Cada `g_{ℓm}` é projetada no espaço de multirresolução `V_J` de uma base
ortonormal de wavelets de suporte compacto em `[0, 1]` (Daublets, Symmlets ou
Coiflets, avaliadas por Daubechies–Lagarias no `WaveBased`):

```
g_{ℓm}(u) ≈ Σ_k α_{ℓm,k} φ_{j0,k}(u) + Σ_{j=j0}^{J−1} Σ_k θ_{ℓm,jk} ψ_{jk}(u)
```

Duas escolhas herdadas do desenho do `wall()` do `WaveBased` (reimplementadas
em `wafc/`, D4):

- **Periodização com `j0 = 0`**: a única função de escala é a constante, e
  ela é **descartada**, absorvida em `c_ℓ`. Isso impõe a restrição de
  identificabilidade `∫ g_{ℓm} = 0` no nível da base, sem restrição
  numérica no ajuste: cada `g_{ℓm}` fica só com as `2^J − 1` wavelets
  `ψ_{jk}`, todas de integral zero. A alternativa `boundary = "interval"`
  (Cohen, Daubechies & Vial 1993) evita o artefato de periodização à custa
  de `j0` maior e de uma restrição explícita.
- **Reescalonamento das moduladoras** para `[ε, 1 − ε]` (`eps` do `wall`),
  que afasta os dados da faixa de borda onde a periodização deixa artefato.

O que importa para a teoria: a coluna da base é o **produto** `X_ℓ ψ_{jk}(U_m)`.
Com `p` lineares, `q` moduladoras e `2^J − 1` wavelets por par, a matriz de
desenho tem `p + p q (2^J − 1)` colunas (as `p` primeiras são os `X_ℓ`,
não penalizados). É um sieve cuja dimensão cresce com `n` via `J = J_n`.

## 3. O estimador

```
(ĉ, θ̂) = argmin  (1/2n) Σ_i (Y_i − Σ_ℓ X_iℓ [c_ℓ + Σ_m Σ_{jk} θ_{ℓm,jk} ψ_{jk}(U_im)])²
                 + λ Σ_{ℓ,m,j,k} w_{ℓm,jk} |θ_{ℓm,jk}|
```

- **LASSO** sobre todos os coeficientes de wavelet, com `c_ℓ` não penalizado
  (`penalty.factor = 0` no `glmnet`), exatamente como o `wall()`. A esparsidade
  no domínio das wavelets é o que dá **adaptatividade espacial**: funções
  `g_{ℓm}` com picos, saltos ou regularidade não homogênea (Besov `B^s_{π,r}`
  com `π < 2`) são bem aproximadas com poucos coeficientes grandes, o que
  splines e kernels com suavidade global não capturam (Donoho & Johnstone
  1994, 1998).
- **Variante em grupos** (a decidir em E2): sparse group LASSO com um grupo
  por par `(ℓ, m)` (Simon, Friedman, Hastie & Tibshirani 2013), que elimina
  uma `g_{ℓm}` inteira e recupera a seleção de "quais moduladoras afetam
  quais coeficientes", no espírito de Huang, Horowitz & Wei (2010) e Wei,
  Huang & Li (2011). O LASSO puro é mais simples, tem teoria mais direta e
  já basta para a adaptatividade; a variante em grupos ganha se a seleção
  de estrutura for uma contribuição do artigo.
- **Sintonia**: `(J, λ)` por validação cruzada (como `cv.wall()`), com BIC ou
  EBIC como alternativa barata a medir em E2.

Pesos `w` iguais a 1 no LASSO base; pesos adaptativos (Zou 2006) são uma
extensão para a propriedade oráculo de seleção, se E1.7 abrir.

## 4. A teoria que se espera provar (e o molde)

O artigo teórico do WALL (`wall-manuscript/manuscript/theo/`) já tem a
arquitetura para o caso logístico aditivo, em cinco passos: viés de
aproximação em Besov; processo empírico dos incrementos da perda;
desigualdade oráculo de taxa lenta; passagem para o risco; e a trilha rápida
com compatibilidade populacional e corolário de compressibilidade
(weak-`ℓ_τ`; o WALL escreve `q`, que aqui é o número de moduladoras).
O WAFC tem **perda quadrática**, o que simplifica os passos 2 a 4, e um
**desenho de produtos**, o que complica o passo da compatibilidade. Resultados
esperados, em ordem de dependência:

1. **Aproximação (E1.3).** Se `g_{ℓm} ∈ B^s_{π,r}[0,1]` com `s > 1/π`, o erro
   de projeção em `V_J` em `L_2(P_{U_m})` é `O(2^{−Js'})` com `s' = s − (1/π −
   1/2)_+`, uniforme em `(ℓ, m)`, quando a densidade de `U_m` é limitada.
2. **Desenho (E1.4).** A matriz de Gram populacional do desenho de produtos,
   `Σ = E[Z Z']` com `Z = (X_ℓ ψ_{jk}(U_m))`, tem autovalor mínimo restrito
   afastado de zero sob: densidade conjunta de `U` limitada por baixo e por
   cima; `E[X X' | U]` com autovalores em `[κ_1, κ_2]`; e uma condição de
   "não colinearidade entre moduladoras" (as `ψ_{jk}(U_m)` para `m`
   distintos não são combinações lineares). A versão empírica segue por
   concentração, com `p q 2^J log(pq 2^J)/n → 0` (ou a condição mais fraca
   de compatibilidade só no cone).
3. **Desigualdade oráculo (E1.5).** Para `λ ≍ σ sqrt(log(p q 2^J)/n)` e erro
   sub-gaussiano: `‖Ẑθ̂ − Zθ*‖²_n + λ‖θ̂ − θ*‖_1 ≲ λ² s_0/φ² + viés²`, com
   `s_0` o número de coeficientes não nulos do oráculo (Bühlmann & van de
   Geer 2011, cap. 6).
4. **Taxas (E1.6).** Com `J_n ≍ log_2 n / (2s' + 1)` a taxa de predição é
   `n^{−2s'/(2s'+1)}` a menos de logaritmos, e o mesmo para cada `ĝ_{ℓm}` em
   `L_2` pela condição de desenho; sob compressibilidade dos coeficientes
   (weak-`ℓ_τ`, `τ < 2`) o LASSO **adapta** à dimensão efetiva, com taxa
   melhor que a do sieve cheio: é o corolário do WALL transposto.
5. **(Opcional, E1.7) Seleção de estrutura.** Com a variante em grupos e pesos
   adaptativos, `P(ĝ_{ℓm} ≡ 0 para todo (ℓ,m) nulo) → 1`.

O que **não** se pretende: inferência (intervalos) para `g_{ℓm}`, taxa minimax
inferior no modelo completo (citar Donoho & Johnstone para o caso univariado
e deixar como observação), e o caso `p` crescendo mais rápido que `n`.

## 5. Por que isso é novo (a defender)

- Xue & Yang (2006) e sucessores usam splines polinomiais, que exigem
  suavidade global (Hölder) e não adaptam a regularidade local.
- Zhou & You (2004) usam wavelets em coeficientes variáveis, mas sem
  estrutura aditiva, sem LASSO e com uma só moduladora.
- **O mais próximo na teoria é Klopp & Pensky (2015)** (L2, `busca-novidade.md`
  §1): para uma moduladora e `X ⊥ U`, eles já têm o desenho de produtos com
  Gram de Kronecker, a concentração da Gram empírica restrita, a desigualdade
  oráculo não assintótica, a taxa adaptativa em Besov e a cota inferior
  minimax. O artigo se apresenta como extensão deles (D18).
- **No método, os mais próximos são Sardy & Ma (2024) e Amato et al. (2022)**:
  aditivos esparsos com wavelets e penalidade `‖·‖_1`, mas **sem o `X_ℓ`
  multiplicando**; o desenho é `ψ(U_m)`, não `X_ℓ ψ(U_m)`. A teoria de Sardy &
  Ma é de otimização, não estatística.
- **Montoril, Morettin & Chiann (2018)** estimam coeficientes funcionais com
  wavelets, sem penalização e com uma moduladora: precursor do próprio grupo,
  obrigatório citar.
- Wei, Huang & Li (2011) e Xue & Qu (2012) fazem seleção em coeficientes
  variáveis por group LASSO com B-splines, com uma moduladora.

A lacuna: coeficientes funcionais **aditivos em várias moduladoras**, base de
**wavelets** (adaptatividade espacial), estimação por **LASSO** com desenho de
produtos, teoria de taxa e adaptação, e software. A busca de novidade (L2)
tem de confirmar que não existe.

## 6. Riscos

| Risco | Sinal | Resposta |
|---|---|---|
| ~~A condição de desenho dos produtos não fecha em geral~~ **fechado em E1.4 (2026-09-18)** | | sob `λ_min(E[XX' \| U]) ≥ κ_1` e densidade conjunta de `U` limitada por baixo, `λ_min(Σ) ≥ κ_1 c_U`, sem cone (D13); a fatoração sob `X ⊥ U` vira observação. O que `κ_1 > 0` exclui é uma combinação linear de `X` ser função de `U` |
| O LASSO puro seleciona coeficientes isolados e as `ĝ_{ℓm}` ficam "espinhosas" | ISE ruim nas funções suaves em E2.4 | sparse group LASSO, ou pós-processamento por limiarização em blocos; medir |
| Periodização deixa artefato de borda nas `g_{ℓm}` | erro concentrado em `u` perto de 0 e 1 | `eps` do `wall`; `boundary = "interval"`; reportar o erro no interior |
| Custo: `p q (2^J − 1)` colunas com `J` grande | `glmnet` lento com `n = 2000`, `p q = 50`, `J = 6` (`~3000` colunas) | desenho esparso do `wall` (`sparse = "auto"`), tabela `wtable()`; o `glmnet` aceita `dgCMatrix` |
| A revista pede aplicação real convincente | não há base natural | E6.1 escolhe entre candidatos com efeito que varia com covariáveis (renda × idade/escolaridade; ambiente × tempo/temperatura); decidir cedo |
| Novidade insuficiente | L2 encontra "additive coefficient + wavelet + lasso" | o peso vai para a teoria de adaptação e para o pacote |
