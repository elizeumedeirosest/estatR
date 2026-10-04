readme <- readLines("C:/Users/Elize/OneDrive/Documentos/estatR/README.md", encoding = "UTF-8")

# Encontrar a secao de gghelpers para incluir os graficos nela
idx_gghelpers <- grep("### \U0001F3A8 Visualiza\u00e7\u00e3o e Gr\u00e1ficos", readme)

# Encontrar a próxima seção (que começa com ### ou ##)
idx_proxima <- grep("^##?#? ", readme[(idx_gghelpers + 1):length(readme)])[1]
idx_proxima_absoluto <- idx_gghelpers + idx_proxima

secao_boxplot <- c(
  "",
  "### \U0001F4CA Gr\u00e1ficos Estat\u00edsticos Prontos",
  "",
  "O `estatR` tamb\u00e9m disponibiliza fun\u00e7\u00f5es gr\u00e1ficas completas com sintaxe limpa (sem aspas):",
  "",
  "```r",
  "# Boxplot B\u00e1sico",
  "grafico_boxplot(mtcars, cyl, mpg)",
  "",
  "# Boxplot Avan\u00e7ado (Agrupado, com dispers\u00e3o e paleta inteligente)",
  "grafico_boxplot(mtcars, cyl, mpg, grupo = am, ",
  "                dispersao_pts = TRUE, paleta = \"vibrant\")",
  "",
  "# Destaque de categorias e ordenamento por mediana",
  "grafico_boxplot(mtcars, cyl, mpg, ",
  "                ordenar = TRUE, destaque = c(\"8\"))",
  "```",
  "",
  "---",
  ""
)

# Inserir antes da próxima seção (que seria Autor)
readme_novo <- c(
  readme[1:(idx_proxima_absoluto - 1)],
  secao_boxplot,
  readme[idx_proxima_absoluto:length(readme)]
)

writeLines(readme_novo, "C:/Users/Elize/OneDrive/Documentos/estatR/README.md", useBytes = TRUE)
cat("README atualizado com grafico_boxplot!\n")
