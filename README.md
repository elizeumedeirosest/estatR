# estatR <img src="man/figures/logo.png" align="right" height="120" alt="" />

> **Uma camada estatística intuitiva e didática para o R.**  
> Saídas em português, formatação limpa e gráficos prontos para uso.

---

## Instalação

```r
install.packages("devtools")
devtools::install_github("elizeumedeirosest/estatR")
```

---

## Módulos Disponíveis

<details>
<summary><strong>► Análise Descritiva e Exploratória de Dados</strong></summary>
<br>

**Diagnóstico de Dados**

```r
diagnostico(mtcars)          # Painel completo: estrutura + valores ausentes
diagnostico(mtcars$mpg)      # Diagnóstico de uma variável

estrutura(mtcars)            # Estrutura do banco de dados
valores_ausentes(mtcars)     # Mapa de valores ausentes
```

**Estatística Descritiva**

```r
# Resumo polimórfico
descrever(mtcars)                        # Resumo completo do banco
descrever(mtcars$mpg)                    # Estatísticas de variável numérica
descrever(mtcars$cyl)                    # Tabela de frequência de variável categórica
descrever(mtcars$mpg, por = mtcars$cyl)  # Numérica agrupada por categoria

# Tabela de contingência (cruzamento de duas variáveis)
tabela_contingencia(mtcars, cyl, am)
tabela_contingencia(mtcars, cyl, am, proporcao = "linha")

# Medidas específicas
med_tend_central(mtcars$mpg)  # Média, mediana, moda
med_dispersao(mtcars$mpg)     # Variância, DP, erro padrão, CV, amplitude
med_forma(mtcars$mpg)         # Assimetria e curtose

# Separatrizes
quartis(mtcars$mpg)
percentis(mtcars$mpg, p = c(5, 25, 50, 75, 95))

# Detecção de anomalias
outliers(mtcars$mpg)          # Limites, quantidade e lista (método IQR)
```

</details>

---

<details>
<summary><strong>► Inferência Estatística</strong></summary>
<br>

**Testes de Hipótese**

```r
# Teste t
teste_t_uma_amostra(mtcars$mpg, mu = 20)
teste_t_duas_amostras(mtcars$mpg, mtcars$am)
teste_t_pareado(dados$antes, dados$depois)

# Teste de proporção
teste_proporcao(x = 87, n = 300, p0 = 0.25)

# Testes não-paramétricos
teste_wilcoxon(mtcars$mpg, mu = 20)
teste_mann_whitney(mtcars$mpg, mtcars$am)

# Qui-quadrado e homogeneidade de variâncias
teste_qui_quadrado(mtcars, cyl, am)
teste_levene(mtcars$mpg, mtcars$am)
```

**Intervalos de Confiança**

```r
ic_media(mtcars$mpg)
ic_media(mtcars$mpg, confianca = 0.99)
ic_proporcao(x = 87, n = 300)
ic_variancia(mtcars$mpg)
ic_diferenca_medias(mtcars$mpg, mtcars$am)
```

**Estimação de Parâmetros**

```r
# Máxima verossimilhança (plota curva de log-verossimilhança)
estimar_verossimilhanca(mtcars$mpg, distribuicao = "normal")
estimar_verossimilhanca(mtcars$wt,  distribuicao = "gamma")

# Método dos momentos (plota histograma com curva teórica)
estimar_momentos(mtcars$mpg, distribuicao = "normal")
estimar_momentos(mtcars$wt,  distribuicao = "gamma")
```

**Simulação Estatística**

```r
# Bootstrap
bootstrap(mtcars$mpg, estatistica = "media",  repeticoes = 2000)
bootstrap(mtcars$mpg, estatistica = "mediana", repeticoes = 2000)

# Teorema Central do Limite
simular_tcl(mtcars$mpg, n_amostra = 30, repeticoes = 1000)

# Monte Carlo
monte_carlo(experimento = function() mean(rnorm(30)) > 0.3, repeticoes = 10000)
```

</details>

---

<details>
<summary><strong>► Probabilidade e Distribuições</strong></summary>
<br>

```r
# Distribuição Normal (com gráfico de área sombreada)
prob_normal(media = 100, dp = 15, q1 = 120, tipo = "maior")
prob_normal(media = 100, dp = 15, q1 = 85, q2 = 115, tipo = "entre")

# Distribuição Binomial
prob_binomial(ensaios = 10, prob = 0.5, q = 5, tipo = "exato")

# Distribuição Poisson
prob_poisson(lambda = 2, q = 2, tipo = "menor")

# Gerar amostras aleatórias com comparação teórica
gerar_amostra(n = 1000, distribuicao = "normal",      media = 50, dp = 5)
gerar_amostra(n = 500,  distribuicao = "binomial",    ensaios = 10, prob = 0.5)
gerar_amostra(n = 300,  distribuicao = "poisson",     lambda = 3)
gerar_amostra(n = 1000, distribuicao = "exponencial", taxa = 0.5)
```

</details>

---

<details>
<summary><strong>► Modelagem e Regressão</strong></summary>
<br>

**Regressão Linear**

```r
# Regressão simples e múltipla
modelo  <- regressao_linear(mpg ~ wt, dados = mtcars)
modelo2 <- regressao_linear(mpg ~ wt + hp, dados = mtcars)

# Predição
predicao(modelo,  data.frame(wt = c(2.5, 3.0)))
predicao(modelo2, data.frame(wt = 2.5, hp = 110), intervalo = "predicao")

# Métricas e diagnóstico de resíduos
metricas(modelo)
analise_residual(modelo)
```

**Correlação**

```r
# Teste de correlação bivariada
teste_correlacao(mtcars, wt, mpg)
teste_correlacao(mtcars, wt, mpg, metodo = "spearman")

# Correlograma completo
matriz_correlacao(mtcars)
matriz_correlacao(mtcars, c("mpg", "wt", "hp", "disp"), metodo = "kendall")
```

**Divisão de Dados**

```r
# Divisão aleatória treino/teste
partes <- dividir_dados(mtcars, proporcao = 0.8, semente = 42)
modelo <- regressao_linear(mpg ~ wt + hp, dados = partes$treino)
predicao(modelo, partes$teste)

# Divisão estratificada
dividir_dados(mtcars, proporcao = 0.8, estrato = cyl, semente = 42)

# Divisão temporal (para séries temporais)
dividir_dados(dados_serie, proporcao = 0.8, temporal = TRUE)
```

</details>

---

<details>
<summary><strong>► Amostragem</strong></summary>
<br>

```r
# Tamanho de amostra
tamanho_amostra(populacao = 50000, erro = 0.03)
tamanho_amostra(populacao = Inf, erro = 0.05, confianca = 0.99)

# Amostragem aleatória simples
amostra_aleatoria(mtcars, n = 10)
amostra_aleatoria(mtcars, proporcao = 0.3, semente = 42)

# Amostragem sistemática
amostra_sistematica(mtcars, n = 10)

# Amostragem estratificada
amostra_estratificada(iris, estrato = Species, n = 30)
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "uniforme")
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "otima", variavel = Sepal.Length)

# Amostragem por conglomerados
amostra_conglomerados(mtcars, conglomerado = cyl, n_conglomerados = 2)
```

</details>

---

<details>
<summary><strong>► Visualização e Gráficos</strong></summary>
<br>

**Gráficos Estatísticos Prontos**

```r
# Boxplot básico
grafico_boxplot(mtcars, cyl, mpg)
grafico_boxplot(mtcars, cyl, mpg, violino = TRUE)

# Boxplot agrupado com dispersão e paleta estatR
grafico_boxplot(mtcars, cyl, mpg, grupo = am,
                dispersao_pts = TRUE, paleta = "vibrant")

# Destaque de categorias e ordenamento por mediana
grafico_boxplot(mtcars, cyl, mpg, ordenar = TRUE, destaque = "8")
```

**Tema e Paletas para ggplot2**

```r
library(ggplot2)

# Aplicar tema padronizado
ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
  geom_boxplot() +
  tema_estatR(estilo = 2) +     # 1 = classic, 2 = minimal com grade, 3 = limpo
  paleta_estatR("academic")     # Paleta de cores padronizada

# Outros ajustes de tema
tema_estatR(modo = "dark")      # Tema escuro
tema_estatR(inclinar = 45)      # Inclina rótulos do eixo X

# Visualizar paletas disponíveis
paleta_estatR()                 # Painel com todas as paletas
paleta_estatR("vibrant")        # Ver uma paleta específica
```

</details>

---

## Autor

**Elizeu Medeiros**  
[GitHub](https://github.com/elizeumedeirosest)

---

## Licença

MIT © Elizeu Medeiros
