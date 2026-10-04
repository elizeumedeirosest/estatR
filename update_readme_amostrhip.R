readme <- readLines("C:/Users/Elize/OneDrive/Documentos/estatR/README.md", encoding = "UTF-8")

# Encontrar a secao de Amostragem e substituir
idx_ini <- grep("### .* Amostragem", readme)
idx_fim <- grep("### .* Regress", readme)

secao_amostragem <- c(
  "### \U0001F3B2 Amostragem",
  "",
  "```r",
  "# Calcular tamanho \u00f3timo de amostra",
  "tamanho_amostra(populacao = 50000, erro = 0.03)",
  "tamanho_amostra(populacao = Inf, erro = 0.05, confianca = 0.99)",
  "",
  "# Amostragem Aleat\u00f3ria Simples",
  "amostra_aleatoria(mtcars, n = 10)",
  "amostra_aleatoria(mtcars, proporcao = 0.3, semente = 42)",
  "",
  "# Amostragem Sistem\u00e1tica",
  "amostra_sistematica(mtcars, n = 10)",
  "amostra_sistematica(mtcars, salto = 3, semente = 1)",
  "",
  "# Amostragem Estratificada (proporcional, uniforme ou \u00f3tima de Neyman)",
  "amostra_estratificada(iris, estrato = Species, n = 30)",
  "amostra_estratificada(iris, estrato = Species, n = 30, alocacao = \"uniforme\")",
  "amostra_estratificada(iris, estrato = Species, n = 30, alocacao = \"otima\", variavel = Sepal.Length)",
  "",
  "# Amostragem por Conglomerados",
  "amostra_conglomerados(mtcars, conglomerado = cyl, n_conglomerados = 2)",
  "```",
  "",
  "---",
  ""
)

secao_hipotese <- c(
  "### \U0001F52C Testes de Hip\u00f3tese",
  "",
  "```r",
  "# Teste T (uma amostra)",
  "teste_t_uma_amostra(mtcars$mpg, mu = 20)",
  "teste_t_uma_amostra(mtcars$mpg, mu = 20, hipotese = \"maior\")",
  "",
  "# Teste T (duas amostras independentes)",
  "teste_t_duas_amostras(mtcars$mpg, mtcars$am)",
  "teste_t_duas_amostras(mtcars$mpg, mtcars$am, hipotese = \"maior\")",
  "",
  "# Teste T Pareado",
  "teste_t_pareado(dados$antes, dados$depois)",
  "",
  "# Teste de Propor\u00e7\u00e3o",
  "teste_proporcao(x = 87, n = 300, p0 = 0.25)",
  "teste_proporcao(x = 87, n = 300, p0 = 0.25, hipotese = \"maior\")",
  "",
  "# Testes N\u00e3o-Param\u00e9tricos",
  "teste_wilcoxon(mtcars$mpg, mu = 20)            # Wilcoxon (uma amostra)",
  "teste_mann_whitney(mtcars$mpg, mtcars$am)      # Mann-Whitney (duas amostras)",
  "",
  "# Teste Qui-Quadrado de Independ\u00eancia",
  "teste_qui_quadrado(mtcars, cyl, am)",
  "",
  "# Teste de Levene (homogeneidade de vari\u00e2ncias)",
  "teste_levene(mtcars$mpg, mtcars$am)",
  "```",
  "",
  "---",
  ""
)

# Montar o novo README
readme_novo <- c(
  readme[1:(idx_ini[1] - 1)],
  secao_amostragem,
  secao_hipotese,
  readme[idx_fim[1]:length(readme)]
)

writeLines(readme_novo, "C:/Users/Elize/OneDrive/Documentos/estatR/README.md", useBytes = TRUE)
cat("README atualizado!\n")
