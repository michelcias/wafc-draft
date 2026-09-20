# Handoff E1.8: a rota do intervalo, o preço da margem e o `s'` dos cenários

Chat de tarefa, 2026-09-19. Tocou só os arquivos do catálogo.

## 1. O que foi feito

| Arquivo | O que é |
|---|---|
| `derivations/07-rota-intervalo.tex` | novo; 11 páginas, compila com `latexmk` sem aviso e sem referência indefinida |
| `derivations/07-rota-intervalo.pdf` | o compilado |
| `derivations/check/07-rota-intervalo.R` | novo; imprime `OK` em 2 min 11 s |
| `docs/handoff-E1.8.md` | este arquivo |

Nenhum arquivo de etapa fechada foi editado. O `.tex` cita E1.2 a E1.6 e
E1.3b sem mexer neles, e **não vai ao manuscrito** (D26).

**Numeração global** (`TAREFA.md` §5): o arquivo continuou de Proposição 4,
Lema 10, Teorema 2 e Corolário 5, e imprime o número global via
`\setcounter`, como `04` e `05`. Ficou com

> **Lema 11, Proposição 5, Corolário 6, Proposição 6, Corolário 7, Lema 12**,
> mais uma observação numerada localmente.

O próximo resultado numera a partir de Proposição 6, Lema 12, Teorema 2 e
Corolário 7.

O arquivo tem três partes independentes, que é a divisão que o catálogo pedia.

**Parte A, a rota do intervalo.** O que transfere da base periodizada para a
base de Cohen, Daubechies e Vial, item a item e com a razão:

- *Aproximação* (Lemas 2 e 3, Corolário 1 de E1.3): vale com
  `B^s_{π,r}[0,1]` no lugar da classe periódica e **a mesma prova**, porque a
  CDV caracteriza o espaço do intervalo pelos coeficientes. Somem o operador
  de extensão, o corte `χ_eps`, a constante `C(eps) ≍ eps^{−(s−1/π)}` e o
  termo `2^{−J}` da Proposição 2, que deixa de ter análogo.
- *Identificabilidade*: já estava feita na §5 de `01-identificabilidade.md`.
  O Lema 1(ii) vale literalmente com `[Φ Q | Ψ]` porque a prova usa só duas
  propriedades das colunas do bloco, ortonormalidade em `L_2[0,1]` e integral
  zero, e a reparametrização entrega as duas.
- *Desenho* (Proposição 3 de E1.4): verbatim, com `C_ψ` trocado por
  `C_ψ^int`. `λ_min(Σ) ≥ κ_1 c_U` sai sem `λ_min(G_eps)` e sem hipótese de
  suporte próprio.
- *Oráculo e taxas*: consomem só essas peças (Corolário 6).

O que **não** transfere de graça, e que é o conteúdo novo da parte:

- **Lema 11**, a cota pontual `Σ ψ²_{jk}(u) ≤ C_ψ 2^J` na CDV, provada pela
  contagem de sobreposição com as `L/2` funções de borda por extremo, que são
  um conjunto fixo de perfis dilatados. A reparametrização `[Φ Q | Ψ]` não
  aumenta a soma, porque `QQ'` é projeção ortogonal.
- **Proposição 5**, o efeito dos `p q (2^{j_0} − 1)` coeficientes de escala
  não penalizados: eles entram no bloco `A` do Lema 4 de E1.5, e o termo
  `σ² p/n` do Teorema 1(iii) vira `σ² p_0/n` com
  `p_0 = p{1 + q(2^{j_0} − 1)}`. Como `p_0` é constante em `n`, as taxas do
  Teorema 2 e dos Corolários 4 e 5 não mudam.

**Parte B, o caso `eps > 0`.** **Proposição 6** dá o cruzamento exato, com
`T(eps) = (L−1)/(2 eps)`: `λ_min(G_eps) > 0` exige `2^J < T(eps)`, e excluir
a faixa contaminada é garantido por `2^J ≥ 2 T(eps)`. As duas faixas são
disjuntas e separadas por **exatamente um nível de resolução**, e manter a
Gram viva ao longo de `J_n → ∞` obriga `eps_n ≍ 2^{−J_n}`, onde o expoente
garantido é exatamente `1/2`, o mesmo de nenhuma margem. **Corolário 7** diz
qual enunciado é verdadeiro: amostra finita, `eps` fixo, `J` na faixa
admissível, com `γ = κ_1 c_U λ_min(G_eps)` e `C̃_g ≍ eps^{−(s−1/π)}`; não há
enunciado de taxa, porque para `eps` fixo todo `J` admissível é `O(1)`.

**Parte C, a análise de `s'`.** **Lema 12**: se a função é `C^N` de cada lado
de um ponto e a primeira derivada que não casa é a de ordem `r < N`, os
coeficientes das `O(1)` translações que cobrem o ponto são de ordem exata
`2^{−j(r+1/2)}`, logo `‖θ_{j·}‖_2 ≍ 2^{−j(r+1/2)}` e a regularidade efetiva
imposta pela singularidade é `r + 1/2`; se não há singularidade, a queda é
`N`, o número de momentos nulos, e não a suavidade da função. A leitura que
organiza a tabela é que, na base periodizada, a emenda `0 ≡ 1` é um ponto
interior para a extensão periódica, com `Δ_r = g^{(r)}(0) − g^{(r)}(1)`.

## 2. O que a conferência numérica mostrou

`Rscript derivations/check/07-rota-intervalo.R` imprime `OK` em 2 min 11 s.

**Parte A.**

- A base CDV é ortonormal em `L_2[0,1]`: `max|G − I|` de `1.8e-06` (`J = 4`)
  a `3.9e-04` (`J = 8`) na grade de `2^14`, que é erro de quadratura.
- A reparametrização fecha: `μ'μ = 1`, `μ'Q = 0` a `1e-10`,
  `(ΦQ)'(ΦQ) = I`, `ΦQ ⊥ Ψ`, e `15 + 48 = 63 = 2^6 − 1` colunas no bloco,
  como no periódico. Toda coluna tem integral zero.
- **`C_ψ^int ≈ 7.55`** (estável em `J = 4..9`, com leve queda até `6.85`),
  contra `1.30` a `1.36` na periódica: **fator 5.6**. A constante entra
  linearmente em `R_J`, logo a condição empírica de E1.4(iv) pede cerca de
  seis vezes o `n` na rota do intervalo. É preço de constante, não de taxa.
- **`p_0 = 93` contra `3`** com `p = 3`, `q = 2`, `j_0 = 4`. Medido em
  `n = 1000`, `σ = 1`, 200 réplicas: `E‖P_A ε‖²_n = 0.0922` contra os
  `0.0930` previstos (e `0.00295` contra `0.00300` no periódico). O termo é
  `0.093 σ²` em `n = 1000` e `0.372 σ²` em `n = 250`.

**Parte B.**

- Família `eps = a (L−1) 2^{−J}`, `J = 5..9`, `L − 1 = 7`. `λ_min(G_eps)`:
  `1` em `a = 0`; `3.5e-04` em `a = 0.25`, estável em `J`; `1.3e-14` em
  `a = 0.5`; nulo em `a = 1`. Queda por nível do erro restrito de `exp(u)`:
  `0.50`, `0.50`, `0.08`, `−0.39`. **Nenhum `a` dá ao mesmo tempo Gram viva e
  queda acima de `1/2` bit.**
- O erro restrito em `J = 9` vai de `2.3e-02` (`eps = 0`) a `2.3e-04`
  (`a = 0.25`, ganho parcial, mas ainda `1/2` bit por nível) e a `1.3e-07` e
  `7.5e-08` (`a = 0.5` e `1`, onde a Gram já morreu). O ganho e a
  degenerescência são o mesmo fenômeno medido de dois lados.
- Com `eps` fixo, o maior `J` com `λ_min(G_eps) > 1e-8` é **`7`, `5` e `4`**
  para `eps = 0.02, 0.05, 0.10`, contra o limiar suficiente `7`, `6`, `5` e
  contra o `9`, `8`, `7` que excluiria a faixa contaminada.

**Parte C.** Nenhuma das seis componentes de `dgp.R` salta na emenda (salto
relativo `≤ 2.5e-05`); a única com quina é a cúbica (`g'(0) = 0.4`,
`g'(1) = 0.6`, salto relativo `0.334`). Queda média por nível medida:

| | `sine` | `cosine` | `cubic` | `bumps` | `blocks` | `heavisine` |
|---|---|---|---|---|---|---|
| (a) periódico `eps = 0` | 3.99 | 3.95 | 1.45 | 1.37 | 0.48 | 0.51 |
| (b) margem e extensão | 3.85 | 3.85 | 3.85 | 1.30 | 0.49 | 0.60 |
| (c) intervalo (CDV) | 5.11 | 4.25 | no piso | 1.37 | 0.49 | 0.58 |

contra a tabela declarada `4, 4, 3/2, 3/2, 1/2, 1/2` em (a) e
`4, 4, 4, 3/2, 1/2, 1/2` em (b) e (c). Daí `s'` do cenário, que é o mínimo
sobre as componentes:

| | `smooth` | `inhomogeneous` |
|---|---|---|
| (a) periódico `eps = 0` | **3/2** | **1/2** |
| (b) margem e extensão | 4 | 1/2 |
| (c) intervalo (CDV) | 4 | 1/2 |

**D27 está confirmada no regime que a teoria usa**, e agora com a razão
escrita.

## 3. O que deve entrar no `ESTADO.md`

### Achados

1. **A rota do intervalo fecha, e o preço é todo de amostra finita.** As
   quatro peças transferem; o que ela cobra é `C_ψ` seis vezes maior,
   `p_0 = 93` coeficientes não penalizados em vez de `3` (com `p = 3`,
   `q = 2`, `L = 8`) e `j_0 ≥ 4`, o que tira `J ∈ {1,2,3,4}` da grade de
   `cv.wafc()` — **exatamente a faixa que E2.3 encontrou como ótima** nos `n`
   do piloto. Assintoticamente não custa nada; nos `n` de E2.4, custa.
2. **O cruzamento da margem é de um nível exato, e é provado.** As duas
   exigências medem a mesma largura, `(L−1)2^{−J}`, e pedem que a margem seja
   simultaneamente mais larga e mais estreita que ela. O `1.9^{−J}` do
   `wall()` fica logo acima do cruzamento, o que explica o ganho parcial de
   `0.96` bit que E1.3b mediu.
3. **A faixa admissível de `J` com `eps = 0.05` é `J ≤ 5`**, e E2.3 achou o
   ótimo da validação cruzada em `J = 3` a `5`. A margem provisória de E2.1b
   está dentro da faixa, mas **por pouco**: `eps = 0.10` já limitaria a
   `J ≤ 4`. É restrição que a varredura de E2.4 tem de respeitar, e é uma
   razão contra aumentar `eps`.
4. **O `1/2` do cenário não homogêneo é independente de regime; o `3/2` do
   cenário suave é inteiramente artefato de periodização.** As três funções
   de Donoho e Johnstone, como `wafc_component()` as escreve, são todas
   contínuas na emenda (em `blocks` as alturas somam zero; em `heavisine` os
   dois termos de sinal se cancelam nos dois extremos; em `bumps` os núcleos
   de borda já decaíram a `1e-05`). A única que paga periodização é a cúbica,
   e paga na **derivada**, não no valor: `g(0) = g(1) = 0` porque
   `1 − 1.4 + 0.4 = 0`.
5. **A cúbica é reproduzida exatamente na base do intervalo:** grau `3 < N = 4`
   e a CDV reproduz polinômios de grau `< N`, de modo que todos os
   coeficientes de wavelet se anulam e a componente inteira mora em
   `V_{j_0}`, na parte não penalizada. Medido: coeficientes em `1e-10` a
   `5e-09`, abaixo do próprio piso de avaliação da base.
6. **O teto de `s'` é `N`, o número de momentos nulos, e não a suavidade da
   função.** Trocar `filter.size` muda o `4` da tabela; D15 fixou `L = 8`.
7. **Lição de ferramenta, útil a qualquer conferência futura:** a base CDV do
   `wbasis()` perde precisão com o nível. O piso de avaliação, medido pela
   norma dos coeficientes que têm de ser exatamente nulos
   (`⟨1, ψ_{jk}⟩`), vai de `7e-08` em `j = 4` a `3e-04` em `j = 12`, contra
   `1e-16` na base periodizada em qualquer nível. Medir taxa de decaimento na
   CDV acima de `j ≈ 7` lê ruído de avaliação, não coeficiente. O piso é
   barato de medir e a conferência o usa como critério de janela.

### Propostas de decisão

- **E1.8a.** Registrar que **a rota do intervalo está escrita e fechada** em
  `derivations/07-rota-intervalo.tex`, e que a resposta a um referee que peça
  teoria sem periodicidade é ela, e não a rota (ii) de D26 (provar a
  compatibilidade sobre as direções estimáveis), que continua sem existir.
  Custo de levar a rota ao manuscrito: reescrever as provas de E1.3 a E1.6
  com a base trocada, o que é redação, não matemática nova.
- **E1.8b.** Registrar a **faixa admissível de `J` como restrição de E2.4**:
  a varredura de `eps` mede `λ_min(G_eps)`, e um `eps` que force `J ≤ 4`
  elimina do desenho os níveis que E2.3 mediu como ótimos. A escolha de
  `wafc_eps_periodic` é um compromisso entre viés de borda e faixa de `J`,
  sem ótimo assintótico.
- **E1.8d.** Uma linha de `derivations/README.md` ficou desatualizada e não
  podia ser tocada por esta tarefa: ela já lista `07-rota-intervalo.tex`
  (E1.8), mas continua dizendo que "`04`, `05` e a emenda de `02` imprimem o
  número global". O `07` também imprime. É conserto de três palavras.
- **E1.8c.** Registrar que o atributo `s'` de `dgp.R` (D27) é o número do
  **regime (a)**, que é o do manuscrito, e que o piloto roda no **regime
  (b)**; ver a pergunta 1 abaixo.

## 4. Perguntas em aberto

1. **O atributo `s'` de `dgp.R` deve declarar o regime.** D27 manda gravar
   `3/2` no cenário suave, que é o número do regime `eps = 0`, o do
   manuscrito. Mas o piloto roda com `rescale = TRUE` e `eps = 0.05`, isto é,
   no regime (b), onde o mesmo cenário lê `4`. A escolha continua defensável
   (o atributo documenta a teoria enunciada, e a Parte B mostra que a margem
   não compra taxa), mas o número deveria vir com o regime escrito ao lado,
   sob pena de E4 comparar o `J_n` teórico com o cenário errado. **Proposta:**
   manter `3/2`, e acrescentar ao atributo um campo com o regime e o valor no
   outro regime. É edição em `wafc/R/dgp.R`, que é de E2.4, não desta tarefa.
2. **A regra teórica de `J_n` muda de um nível com o regime declarado.** Com
   `2^{J_n} ≥ (n/log n)^{1/(2s'+1)}` e `n = 1000`: `s' = 3/2` dá `J_n = 2` e
   `s' = 4` dá `J_n = 1`. E2.3 já media a regra dois níveis abaixo do oráculo
   com `s' = 3/2`; com `4` seriam três. Isso **reforça** a conclusão de E2.3
   ("a regra é muito sensível a `s'`, que não se conhece"), e agora com o
   agravante de que `s'` depende até da convenção de borda. Vale uma frase na
   Seção 5 do artigo? É decisão de E5b.
3. **A cota do Lema 11 não é afirmada justa.** O `C_ψ^int ≈ 7.55` medido é
   para Daublets de filtro 8 e `j_0 = 4`; nada diz que a dependência em `j_0`
   seja a melhor possível. Não bloqueia nada.
4. **O limiar de degenerescência da Proposição 6(i) é suficiente, não
   necessário.** A conferência acha a primeira direção nula até um nível antes
   dele, porque uma combinação de funções de escala pode se anular fora da
   margem sem que nenhuma delas caiba nela. A faixa admissível prática é a
   medida, não a do limiar; quem quiser o limiar justo tem trabalho pela
   frente, e nada depende disso hoje.
5. **`[VERIFICAR]` herdado, não resolvido aqui.** O `02-aproximacao-besov.tex`
   ainda carrega duas marcas de localização interna não conferida (Triebel
   1983 para o ingrediente (E2), Härdle et al. 1998 para o mergulho de Besov
   em Hölder). Este arquivo não depende de nenhuma das duas: a Parte A usa
   Cohen (2003, §3.9) e a Parte B usa só as observações já provadas de E1.3b.
