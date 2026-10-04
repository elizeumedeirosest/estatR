readme <- readLines("C:/Users/Elize/OneDrive/Documentos/estatR/README.md", encoding = "UTF-8")

# Encontrar a secao Autor para inserir antes dela
idx_autor <- grep("## Autor", readme)

secao_gghelpers <- c(
  "### \U0001F3A8 Visualiza\u00e7\u00e3o e Gr\u00e1ficos",
  "",
  "O `estatR` oferece helpers para estilizar gr\u00e1ficos do `ggplot2` de forma r\u00e1pida e limpa.",
  "",
  "```r",
  "library(ggplot2)",
  "",
  "# Visualizar todas as paletas dispon\u00edveis",
  "paleta_estatR()",
  "",
  "# Aplicar tema e paleta padronizados em um gr\u00e1fico",
  "ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +",
  "  geom_boxplot() +",
  "  tema_estatR(estilo = 2) +       # Tema minimalista com grade inteligente",
  "  paleta_estatR(\"academic\")       # Aplica a paleta 'academic'",
  "",
  "# Outros ajustes do tema",
  "tema_estatR(modo = \"dark\")        # Tema escuro",
  "tema_estatR(inclinar = 45)        # Inclina os r\u00f3tulos do eixo X",
  "```",
  "",
  "---",
  ""
)

# Inserir
readme_novo <- c(
  readme[1:(idx_autor[1] - 1)],
  secao_gghelpers,
  readme[idx_autor[1]:length(readme)]
)

writeLines(readme_novo, "C:/Users/Elize/OneDrive/Documentos/estatR/README.md", useBytes = TRUE)
cat("README atualizado com gghelpers!\n")
