# Teoria da Credibilidade — Laboratório 2, Parte 1 (2026/1)

**Previsão de Mortalidade em Fundo de Pensão via Bühlmann-Straub**

> Disciplina: Teoria da Credibilidade (2026/1) — DME/IM-UFRJ  
> Autor: Arthur Pontes Motta  
> Professora: Viviana G. R. Lobo

---

## Relatório interativo

O trabalho está publicado como página web com código, gráficos interativos e tabelas completas:

### **[https://arthurpmotta02.github.io/credibilidade-mortalidade-efpc/](https://arthurpmotta02.github.io/credibilidade-mortalidade-efpc/)**

---

## Sobre o trabalho

Análise de mortalidade de uma **Entidade Fechada de Previdência Complementar (EFPC)** brasileira com dados de expostos e óbitos por idade para o triênio 2012–2014. O objetivo central é estimar o número esperado de sinistros $\hat{N}_i^*$ por classe de risco $i$ para o próximo ano, com base em três abordagens de crescente complexidade probabilística.

Os dados cobrem 116 idades ($x = 0, \ldots, 115$), $v = 417{.}964$ pessoas-ano de exposição e $N = 3{.}658$ óbitos ao longo de $T = 3$ anos. As classes de risco são cinco faixas etárias: $[0,18)$, $[18,40)$, $[40,60)$, $[60,80)$ e $80{+}$.

---

## Estrutura do trabalho

| Questão | Conteúdo |
|---|---|
| (i) | Análise exploratória: pirâmide etária, distribuição de expostos e óbitos, $\log\hat{q}_x$ por idade e ano, heatmap de Lexis (superfície discreta), métricas de maturidade do plano |
| (ii) | Construção das $m = 5$ faixas de risco com critérios de homogeneidade interna, exposição $v_{it} = E_{it}$ suficiente e relevância atuarial |
| (iii) | Seleção de tábuas de referência (AT-2000 e BR-EMSsb-m v.2021), razão A/E por faixa e 11 testes formais de aderência via `mortalityAdherence` |
| (iv) | Definição dos perfis de risco e exposições $v_{it}$ para o modelo de credibilidade |
| (v) | **Bühlmann-Straub:** estimação iterativa de $\lambda_0$ e $\tau^2$, fatores $\hat{\omega}_i$ e previsão $\hat{N}_i^* = v_i^* \cdot \hat{\lambda}_i^H$ |
| (vi) | **Poisson-Gama** (equivalência analítica com B-S) e **Poisson hierárquico bayesiano** via `brms`/Stan: diagnóstico MCMC, PPC, priori vs. posteriori de $\sigma_u$ e IC 95% preditivo do total $\sum_i \hat{N}_i^*$ |

---

## Modelos implementados

### 1. Bühlmann-Straub (questão v)

Estimador linear de credibilidade para a frequência de sinistros. Para cada grupo $i = 1,\ldots,5$ e ano $t = 1,2,3$, define-se $F_{it} = N_{it}/v_{it}$, com $E[N_{it}|\theta_i] = \text{Var}[N_{it}|\theta_i] = v_{it}\lambda_0\theta_i$ (hipótese Poisson condicional).

O prêmio de credibilidade homogêneo é:

$$\hat{\lambda}_i^H = \hat{\omega}_i \bar{F}_i^v + (1 - \hat{\omega}_i)\hat{\lambda}_0, \qquad \hat{\omega}_i = \frac{v_i}{v_i + \hat{\kappa}}, \qquad \hat{\kappa} = \frac{\hat{\lambda}_0}{\hat{\tau}^2}$$

onde $v_i = \sum_t v_{it}$ e $\bar{F}_i^v = \sum_t (v_{it}/v_i) F_{it}$. Os parâmetros estruturais $\hat{\lambda}_0$ e $\hat{\tau}^2$ são estimados por momentos não-viciados via algoritmo iterativo (convergência em 5 iterações).

**Resultados:**

| Parâmetro | Valor |
|---|---|
| $\hat{\lambda}_0$ | $0{,}017921$ |
| $\hat{\tau}^2$ | $0{,}000415$ |
| $\hat{\kappa} = \hat{\lambda}_0/\hat{\tau}^2$ | $43{,}18$ |
| $\hat{\omega}_i$ (grupos com óbitos) | $> 0{,}999$ |
| $\hat{N}^* = \sum_i \hat{N}_i^*$ | $\approx 1{.}267$ sinistros |

---

### 2. Poisson-Gama paramétrico (questão vi — parte 1)

Exploração da equivalência analítica entre o prêmio de Bühlmann-Straub e o prêmio de Bayes do modelo conjugado com priori $\theta_i \sim \text{Gama}(\alpha, \alpha)$:

$$\theta_i \mid \mathbf{N}_i \sim \text{Gama}(\alpha + N_i,\; \alpha + v_i\lambda_0)$$

O prêmio de Bayes é:

$$\hat{\lambda}_i^{\text{Bayes}} = \frac{\alpha}{\alpha + v_i\lambda_0}\,\lambda_0 + \frac{v_i\lambda_0}{\alpha + v_i\lambda_0}\,\bar{F}_i^v$$

que coincide com $\hat{\lambda}_i^H$ quando $\kappa = \alpha$. Inclui visualização das distribuições $p(\theta_i)$ e $p(\theta_i|\mathbf{N}_i)$ para cada faixa etária.

---

### 3. Poisson hierárquico bayesiano (questão vi — parte 2)

Modelo Poisson-LogNormal com efeitos aleatórios por grupo, ajustado por MCMC (NUTS via `brms`):

$$N_{it} \sim \text{Poisson}(\mu_{it}), \qquad \log\mu_{it} = \beta_0 + u_i + \log v_{it}$$

$$\beta_0 \sim \mathcal{N}(-4{,}7,\; 1), \qquad u_i \sim \mathcal{N}(0,\; \sigma_u^2), \qquad \sigma_u \sim \text{Exponencial}(1)$$

A taxa estimada para o grupo $i$ é $\hat{\lambda}_i = e^{\hat{\beta}_0 + \hat{u}_i}$.

**Resultados (4 cadeias, 2.000 amostras pós-aquecimento, `adapt_delta = 0.99`):**

| Parâmetro | Estimativa | IC 95% | $\hat{R}$ | ESS |
|---|---|---|---|---|
| $\beta_0$ | $-5{,}21$ | $(-6{,}60;\; -3{,}81)$ | $1{,}00$ | $1{.}739$ |
| $\sigma_u$ | $2{,}25$ | $(1{,}21;\; 4{,}10)$ | $1{,}00$ | $1{.}689$ |

**Previsão:** mediana $\approx 1{.}266$ sinistros, IC 95% preditivo $[1{.}196;\; 1{.}341]$.  
O P95 $\approx 1{.}341$ é o carregamento de segurança para provisionamento de reservas.

A concordância entre os três modelos confirma empiricamente a equivalência teórica: o estimador linear de Bühlmann-Straub é o prêmio de Bayes ótimo sob perda quadrática no modelo Poisson-Gama.

---

## Tecnologias

| Pacote | Uso |
|---|---|
| `brms` + Stan | Modelo hierárquico bayesiano com MCMC/NUTS |
| `bayesplot` | Diagnóstico MCMC (`mcmc_trace`, `pp_check`) |
| `mortalityAdherence` | 11 testes formais de aderência às tábuas |
| `scico` | Paletas de cor perceptualmente uniformes (Fabio Crameri) |
| `ggiraph` | Gráficos ggplot2 interativos com tooltips |
| `gt` | Tabelas publicáveis com formatação LaTeX |
| `patchwork` | Composição de múltiplos gráficos |

---

## Arquivos

```
.
├── index.qmd   # Documento-fonte (Quarto)
├── referencias.bib             # Referências bibliográficas (ABNT)
├── dadosfundopensao.csv        # Dados do fundo: E_{it} e N_{it} por idade, 2012–2014
├── at2000_iba.csv              # Tábua AT-2000 masculina
├── brems2021_iba.csv           # Tábua BR-EMSsb-m v.2021
└── README.md
```

---

## Reprodução local

### Dependências R

```r
install.packages(c(
  "tidyverse", "gt", "scales", "patchwork",
  "scico", "ggiraph", "kableExtra",
  "brms", "bayesplot"
))

# mortalityAdherence (GitHub)
remotes::install_url(
  "https://github.com/eduardoflm/mortalityAdherence/archive/refs/heads/main.zip"
)
```

### Renderização

```bash
quarto render AtividadeLab2_parte1.qmd
```

O chunk `bayes-fit` usa `#| cache: true` — a primeira renderização compila o modelo Stan (~1–2 min) e as seguintes usam o cache. Para forçar recompilação: `quarto render --cache-refresh`.

---

## Referências principais

- Bühlmann, H.; Gisler, A. *A Course in Credibility Theory and Its Applications*. Springer, 2005.
- Klugman, S. A.; Panjer, H. H.; Willmot, G. E. *Loss Models: From Data to Decisions*. 4. ed. Wiley, 2012.
- Bürkner, P.-C. brms: An R Package for Bayesian Multilevel Models Using Stan. *Journal of Statistical Software*, 80(1), 2017.
- Gelman, A. et al. *Bayesian Data Analysis*. 3. ed. CRC Press, 2014.
- de Melo, E. F. L.; Graziadei, H.; Targino, R. `mortalityAdherence`: Mortality Table Adherence Tests. GitHub, 2026.
- Rau, R. et al. *Visualizing Mortality Dynamics in the Lexis Diagram*. Springer, 2018.
- Lord, D.; Miranda-Moreno, L. F. Effects of low sample mean values and small sample size on the estimation of the fixed dispersion parameter of Poisson-gamma models. *Safety Science*, 46, 2008.
- Park, B.-J.; Lord, D.; Hart, J. D. Bias properties of Bayesian statistics in finite mixture of negative binomial regression models in crash data analysis. *Accident Analysis & Prevention*, 42(2), 2010.