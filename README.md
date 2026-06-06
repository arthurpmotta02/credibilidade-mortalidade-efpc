# Teoria da Credibilidade — Laboratório 2, Parte 1 (2026/1)

**Previsão de Mortalidade em Fundo de Pensão via Bühlmann-Straub**

> Disciplina: Teoria da Credibilidade (2026/1) — DME/IM-UFRJ
> Autor: Arthur Pontes Motta

---

## Relatório interativo

O trabalho está publicado como página web com código, gráficos e tabelas completas:

**[https://arthurpmotta02.github.io/credibilidade-mortalidade-efpc/](https://arthurpmotta02.github.io/credibilidade-mortalidade-efpc/)**

---

## Sobre o trabalho

Análise de mortalidade de uma **Entidade Fechada de Previdência Complementar (EFPC)** brasileira com dados de expostos e óbitos por idade para o triênio 2012–2014. O objetivo central é estimar o número esperado de sinistros por classe de risco para o próximo ano.

Os dados cobrem 116 idades (0–115), 417.964 pessoas-ano de exposição e 3.658 óbitos ao longo de três anos. As faixas etárias definidas como grupos de risco são: `[0,18)`, `[18,40)`, `[40,60)`, `[60,80)` e `80+`.

---

## Estrutura do trabalho

| Questão | Conteúdo |
|---|---|
| (i) | Análise exploratória: pirâmide etária, distribuição de expostos e óbitos, log-taxa de mortalidade, heatmap de Lexis e métricas de maturidade do plano |
| (ii) | Construção das faixas etárias de risco com critérios de homogeneidade interna, exposição suficiente e relevância atuarial |
| (iii) | Seleção e avaliação de tábuas de referência (AT-2000 e BR-EMSsb-m v.2021) com razão A/E e 11 testes formais de aderência via `mortalityAdherence` |
| (iv) | Definição dos perfis de risco para o modelo de credibilidade |
| (v) | **Bühlmann-Straub:** estimação iterativa de λ₀ e τ², fatores de credibilidade ω̂ᵢ e previsão de sinistros por faixa |
| (vi) | **Poisson-Gama** (equivalência analítica com B-S) e **Poisson hierárquico bayesiano** via `brms`/Stan com diagnóstico MCMC, PPC, priori vs. posteriori de σᵤ e IC 95% preditivo |

---

## Modelos implementados

### 1. Bühlmann-Straub (questão v)

Estimador linear de credibilidade com pesos de volume v_{it} = E_{it} (exposição anual). Os parâmetros estruturais λ₀ e τ² são estimados pelos estimadores de momentos não-viciados via algoritmo iterativo.

Resultados principais:
- κ̂ ≈ 43,2 (razão λ₀/τ²)
- ω̂ᵢ > 0,999 para todos os grupos com óbitos positivos
- Previsão: **≈ 1.267 sinistros** (com exposição de 2014)

### 2. Poisson-Gama paramétrico (questão vi — parte 1)

Exploração da equivalência analítica entre o prêmio de Bühlmann-Straub e o prêmio de Bayes do modelo Poisson-Gama conjugado (θᵢ ~ Gama(α, α), com κ = α). Inclui visualização das distribuições a priori e a posteriori de θᵢ para cada faixa etária.

### 3. Poisson hierárquico bayesiano (questão vi — parte 2)

Modelo: N_{it} ~ Poisson(v_{it} · exp(β₀ + uᵢ)), uᵢ ~ N(0, σ²ᵤ), ajustado com NUTS via `brms` (4 cadeias, 2.000 amostras pós-aquecimento, adapt_delta = 0.99).

Resultados:
- β̂₀ = −5,21 (IC 95%: −6,60 ; −3,81) — R̂ = 1,00
- σ̂ᵤ = 2,25 (IC 95%: 1,21 ; 4,10) — ESS > 1.500
- Previsão: **≈ 1.266 sinistros**, IC 95% preditivo: **[1.196 ; 1.341]**

A concordância entre os modelos é esperada e confirma empiricamente a equivalência teórica. O diferencial do modelo bayesiano é o IC 95% preditivo do total da carteira — o P95 ≈ 1.341 sinistros representa o carregamento de segurança para provisionamento de reservas.

---

## Arquivos

```
.
├── AtividadeLab2_parte1.qmd   # Documento-fonte (Quarto)
├── referencias.bib             # Referências bibliográficas
├── dadosfundopensao.csv        # Dados do fundo (expostos e óbitos por idade, 2012-2014)
├── at2000_iba.csv              # Tábua AT-2000 masculina
├── brems2021_iba.csv           # Tábua BR-EMSsb-m v.2021
└── README.md
```

---

## Reprodução local

### Dependências R

```r
install.packages(c(
  "tidyverse", "gt", "scales", "patchwork", "viridisLite",
  "kableExtra", "brms", "bayesplot"
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

O chunk `bayes-fit` usa `#| cache: true` — a primeira renderização compila o modelo Stan (~1–2 min); as seguintes usam o cache automaticamente.

---

## Referências principais

- Bühlmann, H.; Gisler, A. *A Course in Credibility Theory and Its Applications*. Springer, 2005.
- Klugman, S. A.; Panjer, H. H.; Willmot, G. E. *Loss Models: From Data to Decisions*. 4. ed. Wiley, 2012.
- Bürkner, P.-C. brms: An R Package for Bayesian Multilevel Models Using Stan. *Journal of Statistical Software*, 80(1), 2017.
- Gelman, A. et al. *Bayesian Data Analysis*. 3. ed. CRC Press, 2014.
- de Melo, E. F. L.; Graziadei, H.; Targino, R. `mortalityAdherence`: Mortality Table Adherence Tests. GitHub, 2026.
- Rau, R. et al. *Visualizing Mortality Dynamics in the Lexis Diagram*. Springer, 2018.
- Lord, D.; Miranda-Moreno, L. F. Effects of low sample mean values and small sample size on the estimation of the fixed dispersion parameter of Poisson-gamma models. *Safety Science*, 46, 2008.