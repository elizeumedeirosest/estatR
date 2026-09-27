# estatR <img src="man/figures/logo.png" align="right" height="120" alt="" />

> **Uma camada estatística intuitiva e didática para o R.**  
> Saídas em português, formatação limpa e gráficos automáticos.

---

## Visão Geral

O estatR é um pacote R que oferece uma interface amigável e didática sobre as principais funções estatísticas do R base. Pensado para **ensino e análise de dados**, ele produz saídas no console com formatação em português e gráficos automáticos, de modo à representar os resultados.

---

## Instalação

```
# Instale o devtools se ainda não tiver
install.packages("devtools")

# Instale o estatR diretamente do GitHub
devtools::install_github("elizeumedeirosest/estatR")
```
---

## Módulos Disponíveis

### 📊 Estatística Descritiva

```
# Visão geral (diagnóstico) de um banco ou variável
diagnostico(mtcars)
diagnostico(mtcars$mpg)

# Estatísticas resumidas e tabelas de frequência (Polimórfica)
descrever(mtcars)                        # Resumo completo do banco
descrever(mtcars$mpg)                    # Estatísticas da variável numérica
descrever(mtcars$cyl)                    # Tabela de frequência da variável categórica
descrever(mtcars$mpg, por = mtcars$cyl)  # Numérica agrupada

# Tabela de contingência (Cruzamento de duas variáveis)
tabela_contingencia(mtcars, cyl, am)
tabela_contingencia(mtcars, cyl, am, proporcao = "linha")

# Medidas Específicas (Tendência Central, Dispersão, Forma)
med_tend_central(mtcars$mpg)
med_dispersao(mtcars$mpg)
med_forma(mtcars$mpg)

# Medidas de Posição (Separatrizes)
quartis(mtcars$mpg)
quintis(mtcars$mpg)
decis(mtcars$mpg)
```

---

### 📐 Amostragem

```
# Calcular tamanho de amostra (proporção ou média)
tamanho_amostra(populacao = 50000, erro = 0.03)
tamanho_amostra(desvio_padrao = 15, erro = 2, confianca = 0.95)

# Amostragem Aleatória Simples
amostra_aleatoria(iris, n = 30)

# Amostragem Sistemática
amostra_sistematica(iris, n = 30)

# Amostragem Estratificada (proporcional, uniforme ou ótima)
amostra_estratificada(iris, estrato = Species, n = 30)
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "otima", variavel = Sepal.Length)

# Amostragem por Conglomerados
amostra_conglomerados(iris, conglomerado = Species, n_conglomerados = 2)
```

---

### 📈 Regressão Linear

```
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

```
# Teste de correlação bivariada (Pearson, Spearman ou Kendall)
teste_correlacao(mtcars, wt, mpg)
teste_correlacao(mtcars, wt, mpg, metodo = "spearman")

# Correlograma completo do banco de dados
matriz_correlacao(mtcars)
matriz_correlacao(mtcars, c("mpg", "wt", "hp", "disp"), metodo = "kendall")
```

---

### 🎲 Probabilidade

```
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

## Autor

**Elizeu S. de Medeiros, Estatístico.**  
[GitHub](https://github.com/elizeumedeirosest)

---

## Licença

MIT © Elizeu Medeiros
