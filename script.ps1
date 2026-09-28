
$text = Get-Content 'C:\Users\Elize\OneDrive\Documentos\estatR\README.md' -Raw -Encoding UTF8
$text = $text -replace '(?s)### 📐 Amostragem.*?`', '### 📐 Amostragem

`
# Calcular tamanho de amostra
tamanho_amostra(populacao = 50000, erro = 0.03)

# Amostragem Aleatória Simples
amostra_aleatoria(iris, n = 30)

# Amostragem Sistemática
amostra_sistematica(iris, n = 30)

# Amostragem Estratificada
amostra_estratificada(iris, estrato = Species, n = 30)
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "otima", variavel = Sepal.Length)

# Amostragem por Conglomerados
amostra_conglomerados(iris, conglomerado = Species, n_conglomerados = 2)
`'
$text = $text -replace '(?s)### 📈 Regressão Linear.*?`', '### 📈 Regressão Linear

`
# Construir o modelo, exibir análise e gráficos automaticamente
modelo_simples <- regressao_linear(dados = mtcars, mpg ~ wt)

# Regressão múltipla com transformações
modelo_multiplo <- regressao_linear(dados = mtcars, mpg ~ log(wt) + hp)

# Métricas de ajuste
metricas(modelo_simples)

# Diagnóstico de resíduos
analise_residual(modelo_simples)
`'
Set-Content -Path 'C:\Users\Elize\OneDrive\Documentos\estatR\README.md' -Value $text -Encoding UTF8

