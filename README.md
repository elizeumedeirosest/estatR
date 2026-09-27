# estatR <img src="man/figures/logo.png" align="right" height="120" alt="" />

> **Uma camada estatística intuitiva e didática para o R.**  
> Saídas em português, formatação limpa e gráficos automáticos integrados com o pacote [`metaR`](https://github.com/elizeumedeirosest/metaR).

---

## Visão Geral

O `estatR` é um pacote R que oferece uma interface amigável e didática sobre as principais funções estatísticas do R base. Pensado para **ensino e análise de dados**, ele produz saídas no console com formatação profissional em português e gráficos automáticos com um tema visual consistente.

---

## Instalação

```r
# Instale o devtools se ainda não tiver
install.packages("devtools")

# Instale o estatR diretamente do GitHub
devtools::install_github("elizeumedeirosest/estatR")
```

> **Dependência recomendada:** Para os gráficos funcionarem com o tema visual padrão, instale também o [`metaR`](https://github.com/elizeumedeirosest/metaR):
> ```r
> devtools::install_github("elizeumedeirosest/metaR")
> ```

---

## Módulos Disponíveis

### 📊 Estatística Descritiva

```r
# Visão geral (diagnóstico) de um banco ou variável
diagnostico(mtcars)

# Tabela de frequência e contingência
tabela_contingencia(mtcars, cyl)
tabela_contingencia(mtcars, cyl, am)

# Estatísticas resumidas numéricas
descrever(mtcars, mpg)
descrever(mtcars, mpg, grupo = cyl)
```

---

### 🎲 Amostragem

```r
# Calcular tamanho de amostra
tamanho_amostra(populacao = 50000, erro = 0.03)

# Amostragem estratificada (proporcional, uniforme ou ótima de Neyman)
amostra_estratificada(iris, estrato = Species, n = 30)
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "otima", variavel = Sepal.Length)
```

---

### 📈 Regressão Linear

```r
# Regressão simples
modelo <- lm(mpg ~ wt, data = mtcars)
regressao_linear(modelo)

# Regressão múltipla
modelo2 <- lm(mpg ~ wt + hp, data = mtcars)
regressao_linear(modelo2)

# Métricas de ajuste
metricas(modelo)

# Diagnóstico de resíduos (tabela + painel 4 gráficos)
analise_residual(modelo)
```

---

### 🔗 Correlação

```r
# Teste de correlação bivariada (Pearson, Spearman ou Kendall)
teste_correlacao(mtcars, wt, mpg)
teste_correlacao(mtcars, wt, mpg, metodo = "spearman")

# Correlograma completo do banco de dados
matriz_correlacao(mtcars)
matriz_correlacao(mtcars, c("mpg", "wt", "hp", "disp"), metodo = "kendall")
```

---

### 🎰 Probabilidade

```r
# Probabilidade da distribuição Normal (com gráfico de área sombreada)
prob_normal(media = 100, dp = 15, q1 = 120, tipo = "maior")
prob_normal(media = 100, dp = 15, q1 = 85, q2 = 115, tipo = "entre")

# Probabilidade da distribuição Binomial
prob_binomial(ensaios = 10, prob = 0.5, q = 5, tipo = "exato")

# Probabilidade da distribuição Poisson
prob_poisson(lambda = 2, q = 2, tipo = "menor")

# Gerar amostras aleatórias com comparação teórica
gerar_amostra(n = 1000, distribuicao = "normal", media = 50, dp = 5)
gerar_amostra(n = 500,  distribuicao = "binomial", ensaios = 10, prob = 0.5)
gerar_amostra(n = 300,  distribuicao = "poisson",  lambda = 3)
gerar_amostra(n = 1000, distribuicao = "exponencial", taxa = 0.5)
```

---

## Roadmap

- [ ] Módulo de Simulações (Bootstrap, Teorema do Limite Central)
- [ ] Módulo de Testes de Hipótese (Teste-t, Qui-Quadrado, Wilcoxon)
- [ ] Módulo de ANOVA e Delineamento Experimental (Tukey, Duncan)

---

## Autor

**Elizeu Medeiros**  
[GitHub](https://github.com/elizeumedeirosest)

---

## Licença

MIT © Elizeu Medeiros
