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
---

## Módulos Disponíveis

### 🔍 Diagnóstico de Dados

```r
# Painel completo (Estrutura dos dados + Valores Ausentes)
diagnostico(mtcars)
diagnostico(mtcars$mpg)

# Funções fragmentadas (caso queira visualizar separadamente)
estrutura(mtcars)
valores_ausentes(mtcars)
```

---

### 📊 Estatística Descritiva

```r
# Estatísticas resumidas completas (Polimórfica)
descrever(mtcars)                        # Resumo completo do banco
descrever(mtcars$mpg)                    # Estatísticas da variável numérica
descrever(mtcars$cyl)                    # Tabela de frequência da variável categórica
descrever(mtcars$mpg, por = mtcars$cyl)  # Numérica agrupada

# Tabela de contingência (Cruzamento de duas variáveis)
tabela_contingencia(mtcars, cyl, am)
tabela_contingencia(mtcars, cyl, am, proporcao = "linha")

# Medidas Específicas
med_tend_central(mtcars$mpg)  # Média, Mediana, Moda
med_dispersao(mtcars$mpg)     # Variância, DP, Erro Padrão, CV, Amplitude
med_forma(mtcars$mpg)         # Assimetria e Curtose

# Separatrizes (Medidas de Posição)
quartis(mtcars$mpg)
quintis(mtcars$mpg)
decis(mtcars$mpg)
percentis(mtcars$mpg, p = c(5, 10, 50, 90, 95))

# Detecção de Anomalias
outliers(mtcars$mpg)          # Limites, quantidade de outliers e lista (método IQR)
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
modelo <- regressao_linear(mpg ~ wt, dados = mtcars)
predicao(modelo, data.frame(wt = c(2.5, 3.0)))

# Regressão múltipla
modelo2 <- regressao_linear(mpg ~ wt + hp, dados = mtcars)
predicao(modelo2, data.frame(wt = 2.5, hp = 110), intervalo = "predicao")

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

### 📏 Intervalos de Confiança

```r
# IC para a média de uma amostra (distribuição t)
ic_media(mtcars$mpg)
ic_media(mtcars$mpg, confianca = 0.99)

# IC para uma proporção (aproximação Normal / Wald)
ic_proporcao(x = 87, n = 300)
ic_proporcao(x = 87, n = 300, confianca = 0.99)

# IC para a variância e desvio padrão (distribuição Qui-Quadrado)
ic_variancia(mtcars$mpg)

# IC para a diferença de médias entre dois grupos
ic_diferenca_medias(mtcars$mpg, mtcars$am)
```

---

### 🎯 Estimação de Parâmetros

```r
# Estimador de Máxima Verossimilhança (EMV)
# Plota a curva de log-verossimilhança com o EMV destacado
estimar_verossimilhanca(mtcars$mpg,  distribuicao = "normal")
estimar_verossimilhanca(mtcars$carb, distribuicao = "poisson")
estimar_verossimilhanca(mtcars$wt,   distribuicao = "gamma")

# Método dos Momentos
# Plota histograma com a curva teórica ajustada
estimar_momentos(mtcars$mpg, distribuicao = "normal")
estimar_momentos(mtcars$wt,  distribuicao = "gamma")
```

---

### 🎲 Simulação Estatística

```r
# Bootstrap: IC e variabilidade de qualquer estatística por reamostragem
bootstrap(mtcars$mpg, estatistica = "media",  repeticoes = 2000, semente = 42)
bootstrap(mtcars$mpg, estatistica = "mediana", repeticoes = 2000, semente = 42)
bootstrap(mtcars$mpg, estatistica = function(x) quantile(x, 0.9), semente = 42)

# Teorema Central do Limite: painel visual distribuicao original vs. medias
simular_tcl(mtcars$mpg, n_amostra = 5,  repeticoes = 1000)
simular_tcl(mtcars$mpg, n_amostra = 30, repeticoes = 1000)

# Monte Carlo: estimar probabilidade de qualquer evento por simulacao
monte_carlo(experimento = function() mean(rnorm(30)) > 0.3, repeticoes = 10000, semente = 42)
```

---

### 🧪 Modelagem

```r
# Divisao treino/teste — modo aleatorio (padrao)
partes <- dividir_dados(mtcars, proporcao = 0.8, semente = 42)
modelo <- regressao_linear(mpg ~ wt + hp, dados = partes$treino)
predicao(modelo, partes$teste)

# Divisao estratificada (mantem proporcoes dos grupos)
dividir_dados(mtcars, proporcao = 0.8, estrato = cyl, semente = 42)

# Divisao temporal (corte cronologico sem embaralhamento — series temporais)
dividir_dados(dados_serie, proporcao = 0.8, temporal = TRUE)
```

## Autor

**Elizeu Medeiros**  
[GitHub](https://github.com/elizeumedeirosest)

---

## Licença

MIT © Elizeu Medeiros
