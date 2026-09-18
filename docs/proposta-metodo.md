# Proposta de método: a ideia inicial

Registro do que foi discutido em 2026-09-18. É um ponto de partida, não uma
especificação: tudo aqui pode ser melhorado, e o que for provado passa para
`derivations/`. Notação (esboço) em [`notacao.md`](notacao.md).

---

## 1. O modelo

Para uma resposta escalar `Y`, `p` covariáveis "lineares" `X = (X_1, …, X_p)`
e `q` covariáveis "moduladoras" `U = (U_1, …, U_q)`, com `U_k` a valores num
intervalo (reescalado para `[0, 1]`):

```
Y = Σ_{j=1}^{p} β_j(U) X_j + ε,          E[ε | X, U] = 0
β_j(u) = c_j + Σ_{k=1}^{q} g_{jk}(u_k),   E[g_{jk}(U_k)] = 0
```

- **Modelo de coeficientes funcionais** (varying-coefficient) no sentido de
  Hastie & Tibshirani (1993) e Cai, Fan & Yao (2000): o efeito de `X_j` é uma
  função de `U`.
- **Estrutura aditiva em cada coeficiente**: cada `β_j` é soma de funções
  univariadas, uma por moduladora. É o *additive coefficient model* de Xue &
  Yang (2006), que escapa da maldição da dimensionalidade em `q`.
- **Casos particulares que o método tem de recuperar**: `X_1 ≡ 1` dá um
  intercepto aditivo `c_1 + Σ_k g_{1k}(U_k)`, logo o modelo aditivo puro
  (`p = 1`) e o modelo parcialmente linear aditivo (`p > 1`, só `β_1`
  funcional) são casos particulares; `q = 1` é o varying-coefficient
  clássico com um índice; `g_{jk} ≡ 0` para todo `k` devolve a regressão
  linear.

A restrição `E[g_{jk}(U_k)] = 0` é a de identificabilidade das componentes
aditivas; `c_j` absorve o nível. Identificar `(c_j, g_{jk})` a partir da
distribuição de `(Y, X, U)` exige que `E[X X' | U]` seja não singular
(quase certamente) e que as `U_k` não sejam funções uma da outra; isso vira
hipótese em E1.2.

## 2. A aproximação por wavelets

Cada `g_{jk}` é projetada no espaço de multirresolução `V_J` de uma base
ortonormal de wavelets de suporte compacto em `[0, 1]` (Daublets, Symmlets ou
Coiflets, avaliadas por Daubechies–Lagarias no `WaveBased`):

```
g_{jk}(u) ≈ Σ_{m} α_{jk,m} φ_{j0,m}(u) + Σ_{l=j0}^{J-1} Σ_{m} θ_{jk,lm} ψ_{lm}(u)
```

Duas escolhas herdadas do `wall()` do `WaveBased`:

- **Periodização com `j0 = 0`**: a única função de escala é a constante, e
  ela é **descartada**, absorvida em `c_j`. Isso impõe a restrição de
  identificabilidade `∫ g_{jk} = 0` no nível da base, sem restrição
  numérica no ajuste: cada `g_{jk}` fica só com as `2^J − 1` wavelets
  `ψ_{lm}`, todas de integral zero. A alternativa `boundary = "interval"`
  (Cohen, Daubechies & Vial 1993) evita o artefato de periodização à custa
  de `j0` maior e de uma restrição explícita.
- **Reescalonamento das moduladoras** para `[ε, 1 − ε]` (`eps` do `wall`),
  que afasta os dados da faixa de borda onde a periodização deixa artefato.

O que importa para a teoria: a coluna da base é o **produto** `X_j ψ_{lm}(U_k)`.
Com `p` lineares, `q` moduladoras e `2^J − 1` wavelets por par, a matriz de
desenho tem `p + p q (2^J − 1)` colunas (as `p` primeiras são os `X_j`,
não penalizados). É um sieve cuja dimensão cresce com `n` via `J = J_n`.

## 3. O estimador

```
(ĉ, θ̂) = argmin  (1/2n) Σ_i (Y_i − Σ_j X_ij [c_j + Σ_k Σ_{lm} θ_{jk,lm} ψ_{lm}(U_ik)])²
                 + λ Σ_{j,k,l,m} w_{jk,lm} |θ_{jk,lm}|
```

- **LASSO** sobre todos os coeficientes de wavelet, com `c_j` não penalizado
  (`penalty.factor = 0` no `glmnet`), exatamente como o `wall()`. A esparsidade
  no domínio das wavelets é o que dá **adaptatividade espacial**: funções
  `g_{jk}` com picos, saltos ou regularidade não homogênea (Besov `B^s_{π,r}`
  com `π < 2`) são bem aproximadas com poucos coeficientes grandes, o que
  splines e kernels com suavidade global não capturam (Donoho & Johnstone
  1994, 1998).
- **Variante em grupos** (a decidir em E2): sparse group LASSO com um grupo
  por par `(j, k)` (Simon, Friedman, Hastie & Tibshirani 2013), que elimina
  uma `g_{jk}` inteira e recupera a seleção de "quais moduladoras afetam
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
com compatibilidade populacional e corolário de compressibilidade (weak-`ℓ_q`).
O WAFC tem **perda quadrática**, o que simplifica os passos 2 a 4, e um
**desenho de produtos**, o que complica o passo da compatibilidade. Resultados
esperados, em ordem de dependência:

1. **Aproximação (E1.3).** Se `g_{jk} ∈ B^s_{π,r}[0,1]` com `s > 1/π`, o erro
   de projeção em `V_J` em `L_2(P_{U_k})` é `O(2^{−Js'})` com `s' = s − (1/π −
   1/2)_+`, uniforme em `(j, k)`, quando a densidade de `U_k` é limitada.
2. **Desenho (E1.4).** A matriz de Gram populacional do desenho de produtos,
   `Σ = E[Z Z']` com `Z = (X_j ψ_{lm}(U_k))`, tem autovalor mínimo restrito
   afastado de zero sob: densidade conjunta de `U` limitada por baixo e por
   cima; `E[X X' | U]` com autovalores em `[κ_1, κ_2]`; e uma condição de
   "não colinearidade entre moduladoras" (as `ψ_{lm}(U_k)` para `k`
   distintos não são combinações lineares). A versão empírica segue por
   concentração, com `p q 2^J log(pq 2^J)/n → 0` (ou a condição mais fraca
   de compatibilidade só no cone).
3. **Desigualdade oráculo (E1.5).** Para `λ ≍ σ sqrt(log(p q 2^J)/n)` e erro
   sub-gaussiano: `‖Ẑθ̂ − Zθ*‖²_n + λ‖θ̂ − θ*‖_1 ≲ λ² s_0/φ² + viés²`, com
   `s_0` o número de coeficientes não nulos do oráculo (Bühlmann & van de
   Geer 2011, cap. 6).
4. **Taxas (E1.6).** Com `J_n ≍ log_2 n / (2s' + 1)` a taxa de predição é
   `n^{−2s'/(2s'+1)}` a menos de logaritmos, e o mesmo para cada `ĝ_{jk}` em
   `L_2` pela condição de desenho; sob compressibilidade dos coeficientes
   (weak-`ℓ_q`, `q < 2`) o LASSO **adapta** à dimensão efetiva, com taxa
   melhor que a do sieve cheio: é o corolário do WALL transposto.
5. **(Opcional, E1.7) Seleção de estrutura.** Com a variante em grupos e pesos
   adaptativos, `P(ĝ_{jk} ≡ 0 para todo (j,k) nulo) → 1`.

O que **não** se pretende: inferência (intervalos) para `g_{jk}`, taxa minimax
inferior no modelo completo (citar Donoho & Johnstone para o caso univariado
e deixar como observação), e o caso `p` crescendo mais rápido que `n`.

## 5. Por que isso é novo (a defender)

- Xue & Yang (2006) e sucessores usam splines polinomiais, que exigem
  suavidade global (Hölder) e não adaptam a regularidade local.
- Zhou & You (2004) usam wavelets em coeficientes variáveis, mas sem
  estrutura aditiva, sem LASSO e com uma só moduladora.
- Sardy & Tseng (2004) e Sardy & Ma (2024) fazem modelos aditivos esparsos com
  wavelets e penalidade `ℓ_1`, mas **sem o `X_j` multiplicando**: o desenho é
  `ψ(U_k)`, não `X_j ψ(U_k)`, e as condições de desenho e a identificabilidade
  mudam. É o trabalho mais próximo e o que o referee vai citar.
- Wei, Huang & Li (2011) e Xue & Qu (2012) fazem seleção em coeficientes
  variáveis por group LASSO com B-splines, com uma moduladora.

A lacuna: coeficientes funcionais **aditivos em várias moduladoras**, base de
**wavelets** (adaptatividade espacial), estimação por **LASSO** com desenho de
produtos, teoria de taxa e adaptação, e software. A busca de novidade (L2)
tem de confirmar que não existe.

## 6. Riscos

| Risco | Sinal | Resposta |
|---|---|---|
| A condição de desenho dos produtos não fecha em geral | contra-exemplo com `X` e `U` dependentes | enunciar sob `X ⊥ U` primeiro (a Gram fatora em `E[XX'] ⊗ E[ψψ']`) e generalizar com hipótese explícita |
| O LASSO puro seleciona coeficientes isolados e as `ĝ_{jk}` ficam "espinhosas" | ISE ruim nas funções suaves em E2.4 | sparse group LASSO, ou pós-processamento por limiarização em blocos; medir |
| Periodização deixa artefato de borda nas `g_{jk}` | erro concentrado em `u` perto de 0 e 1 | `eps` do `wall`; `boundary = "interval"`; reportar o erro no interior |
| Custo: `p q (2^J − 1)` colunas com `J` grande | `glmnet` lento com `n = 2000`, `p q = 50`, `J = 6` (`~3000` colunas) | desenho esparso do `wall` (`sparse = "auto"`), tabela `wtable()`; o `glmnet` aceita `dgCMatrix` |
| A revista pede aplicação real convincente | não há base natural | E6.1 escolhe entre candidatos com efeito que varia com covariáveis (renda × idade/escolaridade; ambiente × tempo/temperatura); decidir cedo |
| Novidade insuficiente | L2 encontra "additive coefficient + wavelet + lasso" | o peso vai para a teoria de adaptação e para o pacote |
