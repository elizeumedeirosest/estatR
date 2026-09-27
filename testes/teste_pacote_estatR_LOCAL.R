# ════════════════════════════════════════════
# SCRIPT: Testes das Funções do Pacote estatR
# AUTOR: Elizeu S. De Medeiros
# DATA: 27/09/2026
# ════════════════════════════════════════════


# ── Pacote ──────────────────
devtools::install_github("elizeumedeirosest/estatR")
library(estatR)



# ════════════════════════════════════════════
# 1. ESTATÍSTICA DESCRITIVA
# ════════════════════════════════════════════

# ────────────────────────────────────────────
# 1.1 Diagnóstico de Dados
# ────────────────────────────────────────────

# Visão geral de um data frame
diagnostico(mtcars)

# Visão geral de um vetor
diagnostico(mtcars$mpg)


# ────────────────────────────────────────────
# 1.2 Estatísticas e Frequências (Polimórfica)
# ────────────────────────────────────────────

# Resumo completo do banco
descrever(mtcars)

# Variável numérica em um data frame
descrever(mtcars$mpg)

# Descrição numérica agrupada por uma variável categórica
descrever(mtcars$mpg, por = mtcars$cyl)

# Tabela de frequência simples (vetor categórico)
descrever(mtcars$cyl)


# ────────────────────────────────────────────
# 1.3 Tabela de Contingência
# ────────────────────────────────────────────

# Cruzamento de duas variáveis categóricas
tabela_contingencia(mtcars, cyl, am)

# Contingência com porcentagem por linha
tabela_contingencia(mtcars, cyl, am, proporcao = "linha")



# ════════════════════════════════════════════
# 2. AMOSTRAGEM
# ════════════════════════════════════════════

# ────────────────────────────────────────────
# 2.1 Tamanho de Amostra
# ────────────────────────────────────────────

# Para proporção — cenário conservador
tamanho_amostra(populacao = 50000, erro = 0.03, confianca = 0.95)

# Para média — com desvio-padrão conhecido
tamanho_amostra(desvio_padrao = 15, erro = 2, confianca = 0.95)


# ────────────────────────────────────────────
# 2.2 Amostragem Estratificada
# ────────────────────────────────────────────

# Alocação proporcional (padrão)
amostra_estratificada(iris, estrato = Species, n = 30, semente = 42)

# Alocação uniforme (igual por estrato)
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "uniforme", semente = 42)

# Alocação ótima de Neyman (usa variabilidade interna)
amostra_estratificada(iris, estrato = Species, n = 30,
                      alocacao = "otima", variavel = Sepal.Length, semente = 42)




# ════════════════════════════════════════════
# 3. REGRESSÃO LINEAR
# ════════════════════════════════════════════

# ────────────────────────────────────────────
# 3.1 Regressão Simples
# ────────────────────────────────────────────

modelo_simples <- lm(mpg ~ wt, data = mtcars)

regressao_linear(modelo_simples)
metricas(modelo_simples)
analise_residual(modelo_simples)


# ────────────────────────────────────────────
# 3.2 Regressão Múltipla
# ────────────────────────────────────────────

modelo_multiplo <- lm(mpg ~ wt + hp + cyl, data = mtcars)

regressao_linear(modelo_multiplo)
metricas(modelo_multiplo)
analise_residual(modelo_multiplo)




# ════════════════════════════════════════════
# 4. CORRELAÇÃO
# ════════════════════════════════════════════

# ────────────────────────────────────────────
# 4.1 Teste de Correlação Bivariada
# ────────────────────────────────────────────

# Pearson (padrão) — variáveis contínuas
teste_correlacao(mtcars, wt, mpg)

# Spearman — mais robusto a outliers
teste_correlacao(iris, Sepal.Length, Petal.Length, metodo = "spearman")


# ────────────────────────────────────────────
# 4.2 Matriz de Correlação
# ────────────────────────────────────────────

# Todas as variáveis numéricas do banco
matriz_correlacao(mtcars)

# Selecionando variáveis específicas
matriz_correlacao(iris, c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"))

# Com Spearman
matriz_correlacao(mtcars, metodo = "spearman")




# ════════════════════════════════════════════
# 5. PROBABILIDADE
# ════════════════════════════════════════════

# ────────────────────────────────────────────
# 5.1 Distribuição Normal
# ────────────────────────────────────────────

# P(X < 1.96) — Normal padrão
prob_normal(media = 0, dp = 1, q1 = 1.96, tipo = "menor")

# P(X > 120) — QI, média 100, DP 15
prob_normal(media = 100, dp = 15, q1 = 120, tipo = "maior")

# P(85 < X < 115) — intervalo central
prob_normal(media = 100, dp = 15, q1 = 85, q2 = 115, tipo = "entre")


# ────────────────────────────────────────────
# 5.2 Distribuição Binomial
# ────────────────────────────────────────────

# P(X = 5) — 10 lançamentos, p = 0.5
prob_binomial(ensaios = 10, prob = 0.5, q = 5, tipo = "exato")

# P(X <= 3) — probabilidade acumulada
prob_binomial(ensaios = 10, prob = 0.3, q = 3, tipo = "menor")


# ────────────────────────────────────────────
# 5.3 Distribuição Poisson
# ────────────────────────────────────────────

# P(X = 2) — média de 3 eventos por hora
prob_poisson(lambda = 3, q = 2, tipo = "exato")

# P(X >= 5)
prob_poisson(lambda = 3, q = 5, tipo = "maior")


# ────────────────────────────────────────────
# 5.4 Geração de Amostras
# ────────────────────────────────────────────

# Normal
gerar_amostra(n = 1000, distribuicao = "normal", media = 50, dp = 5, semente = 42)

# Binomial
gerar_amostra(n = 500, distribuicao = "binomial", ensaios = 10, prob = 0.5, semente = 42)

# Poisson
gerar_amostra(n = 500, distribuicao = "poisson", lambda = 3, semente = 42)

# Exponencial
gerar_amostra(n = 1000, distribuicao = "exponencial", taxa = 0.5, semente = 42)
