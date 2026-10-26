<div align="center">

<img src="man/figures/logo.svg" alt="estatR Logo" width="150" height="150" />

# estatR

### Estatística com Elegância, Precisão e Produtividade em R

[![Documentação Oficial](https://img.shields.io/badge/docs-oficial-0071e3?style=for-the-badge&logo=google-chrome&logoColor=white)](https://elizeumedeirosest.github.io/estatR) 
[![Licença MIT](https://img.shields.io/badge/licença-MIT-10b981?style=for-the-badge)](LICENSE) 
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D_4.0-38bdf8?style=for-the-badge&logo=R&logoColor=white)](https://cran.r-project.org/) 
[![Funções](https://img.shields.io/badge/funções-120+-818cf8?style=for-the-badge)](https://elizeumedeirosest.github.io/estatR/reference/index.html)

<p align="center">
  <b>Uma camada estatística moderna, didática e unificada para o ecossistema R.</b><br>
  Saídas estruturadas em português, formatação ANSI de alto padrão, tabelas prontas para publicação e gráficos editoriais em um único pacote.
</p>

[**Explorar Documentação Interativa →**](https://elizeumedeirosest.github.io/estatR) &nbsp;•&nbsp; [**Referência de Funções**](https://elizeumedeirosest.github.io/estatR/reference/index.html) &nbsp;•&nbsp; [**Instalação Rápida**](#instalação)

</div>

---

### ✨ Pilares de Design do estatR

| 🎨 Padrão Editorial Automático | ⚡ Sintaxe Fluida & Polimórfica | 📋 Saídas Estruturadas & ANSI | 🔬 Rigor Matemático Completo |
| :--- | :--- | :--- | :--- |
| Gráficos com paleta corporativa e tema dual grid (`grade = 'dupla'`). Prontos para teses, artigos e relatórios corporativos sem necessidade de formatação manual. | Funções inteligentes que detectam o tipo de entrada (`descrever()`, `metricas()`, `analise_residual()`) e entregam diagnósticos sob medida. | Cabeçalhos visualmente estilizados no console, alinhamento rigoroso de colunas, casas decimais padronizadas e notas interpretativas explicativas. | Do diagnóstico preliminar à modelagem linear generalizada (GLM/Youden), séries temporais com machine learning, sobrevivência e multivariada. |

---

## Instalação

```r
install.packages("devtools")
devtools::install_github("elizeumedeirosest/estatR")
```

---

## Módulos Disponíveis

<details>
<summary><strong> Análise Descritiva e Exploratória de Dados</strong></summary>
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
<summary><strong> Inferência Estatística</strong></summary>
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
# Bootstrap Simples (vetor independente)
bootstrap(mtcars$mpg, estatistica = "media",  repeticoes = 2000)

# Bootstrap Estratificado (reamostragem balanceada dentro de subgrupos)
bootstrap(mtcars$mpg, estatistica = "media", tipo = "estratificado", estrato = mtcars$am)

# Bootstrap Residual para Regressão (distribuições e ICs robustos para coeficientes e R²)
bootstrap(modelo = modelo2, tipo = "residual", repeticoes = 1000)

# Bootstrap de Blocos Móveis (Moving Block Bootstrap para séries temporais)
bootstrap(serie_temp, tipo = "bloco", repeticoes = 1000)

# Teorema Central do Limite
simular_tcl(mtcars$mpg, n_amostra = 30, repeticoes = 1000)

# Monte Carlo
monte_carlo(experimento = function() mean(rnorm(30)) > 0.3, repeticoes = 10000)
```

</details>

---

<details>
<summary><strong> Probabilidade e Distribuições</strong></summary>
<br>

**Distribuições (Parâmetros, Momentos Teóricos e Quantis)**

```r
# Visão Geral Didática de Distribuições
distribuicao("normal", media = 100, dp = 15)
distribuicao("poisson", lambda = 3)
distribuicao("t", gl = 5)
distribuicao("qui_quadrado", gl = 3)
distribuicao("f", gl1 = 2, gl2 = 15)
distribuicao("uniforme", min = 0, max = 10)
distribuicao("gama", forma = 2, escala = 2)

# Painel Visual de Simulação (Lei dos Grandes Números / Convergência)
distribuicao("exponencial", taxa = 0.5, simulacao = TRUE)
```

**Cálculo de Probabilidade Específica**

```r
# Cálculo 
prob_normal(media = 100, dp = 15, q1 = 85, q2 = 115, tipo = "entre")
prob_binomial(ensaios = 10, prob = 0.5, q = 5, tipo = "exato")
prob_poisson(lambda = 2, q = 2, tipo = "menor")
prob_exponencial(taxa = 0.5, q1 = 2, tipo = "maior")
```

**Amostras Aleatórias de Distribuiçoes Específicas**
```r
# Gerar amostras aleatórias com comparação teórica
gerar_amostra(n = 1000, distribuicao = "normal", media = 50, dp = 5)
gerar_amostra(n = 500,  distribuicao = "binomial", ensaios = 10, prob = 0.5)
gerar_amostra(n = 1000, distribuicao = "exponencial", taxa = 0.5)
```
</details>

---

<details>
<summary><strong> Modelagem e Regressão</strong></summary>
<br>

**Modelos Lineares Generalizados (MLG / GLM)**

```r
# Função
# Suporta: "binomial", "poisson", "quasipoisson", "gama", "gaussian", "binomial_negativa", "beta"
mod_mlg <- regressao_mlg(mtcars, mpg ~ wt + hp, familia = "gaussian", teste = "ambos")

# Gráfico de Coeficientes estilo Forest Plot
grafico_coeficientes(mod_mlg)
```

**Regressão Linear**

```r
# Regressão simples e múltipla
modelo1  <- regressao_linear(mpg ~ wt, dados = mtcars)
modelo2 <- regressao_linear(mpg ~ wt + hp, dados = mtcars)

# Predição
predicao(modelo1,  data.frame(wt = c(2.5, 3.0)))
predicao(modelo2, data.frame(wt = 2.5, hp = 110), intervalo = "predicao")

# Métricas e diagnóstico de resíduos polimórfico
metricas(modelo1)
analise_residual(modelo1)
grafico_coeficientes(modelo1)

# Diagnóstico de Alavancagem e Pontos Influentes (Cook, Leverage, DFFITS e reajuste comparativo)
pontos_influentes(modelo1)
pontos_influentes(modelo2)

# Diagnóstico de Multicolinearidade (VIF para numéricas, GVIF para fatores e Número de Condição)
multicolinearidade(modelo2)
```

**Regressão Logística (Binomial, Multinomial e Ordinal)**

```r
# Logística Binária com definição de categoria de referência
mod_bin <- regressao_logistica(vs ~ wt + hp, dados = mtcars, referencia = 0)

# Logística Multinomial (3+ níveis não ordenados)
mod_multi <- regressao_logistica(Species ~ Sepal.Length + Petal.Length, 
                                 dados = iris, tipo = "nominal")

# Logística Ordinal (3+ níveis ordenados, Odds Proporcionais)
# (Requer que a variável resposta seja factor ordenado - ordered)
mod_ord <- regressao_logistica(nota ~ horas + esforco, 
                               dados = dados_alunos, tipo = "ordinal")

# Métricas, diagnóstico de resíduos e Best Subsets
metricas(mod_bin)
analise_residual(mod_bin)
grafico_coeficientes(mod_bin) # Plota Odds Ratios (OR) automaticamente
selecao_modelos(vs ~ wt + hp + qsec, dados = mtcars)
```

**Regressão de Poisson**

```r
mod_pois <- regressao_poisson(breaks ~ wool + tension, dados = warpbreaks)

# Funções polimórficas (adaptam-se automaticamente ao modelo de Poisson)
metricas(mod_pois)
analise_residual(mod_pois)
grafico_coeficientes(mod_pois) # Plota Razões de Taxas de Incidência (IRR) automaticamente
selecao_modelos(breaks ~ wool + tension, dados = warpbreaks)
```

**Regressão Binomial Negativa (Contagens com Sobredispersão)**

```r
mod_nb <- regressao_binomial_negativa(breaks ~ wool + tension, dados = warpbreaks)

# Métricas globais, análise de resíduos e seleção
metricas(mod_nb)
analise_residual(mod_nb)
grafico_coeficientes(mod_nb)
selecao_modelos(breaks ~ wool + tension, dados = warpbreaks)
```

**Regressão Beta (Proporções e Taxas)**

```r
# A variável resposta deve estar no intervalo (0, 1) exclusivo.
mod_beta <- regressao_beta(taxa_sucesso ~ metodo + experiencia, dados = dados_taxas)

# Métricas, diagnóstico de resíduos e seleção polimórficos
metricas(mod_beta)
analise_residual(mod_beta)
grafico_coeficientes(mod_beta)
selecao_modelos(taxa_sucesso ~ metodo + experiencia + idade, dados = dados_taxas)
```

**Regressão Gama**

```r
# A variável resposta deve ser estritamente positiva (> 0). Ideal para custos e tempos.
mod_gama <- regressao_gama(custos ~ tratamento + idade, dados = dados_saude)

# Métricas, resíduos e seleção polimórficos
metricas(mod_gama)
analise_residual(mod_gama)
grafico_coeficientes(mod_gama)
selecao_modelos(custos ~ tratamento + idade + bmi, dados = dados_saude)
```

**Correlação**

```r
# Teste de correlação bivariada
teste_correlacao(mtcars, wt, mpg)
teste_correlacao(mtcars, wt, mpg, metodo = "spearman")

# Matriz completa de correlação com rankings TOP (positivas e negativas) e correlograma:
matriz_correlacao(mtcars, top = 5)

# Matriz e gráfico filtrados apenas para as variáveis com maiores correlações:
matriz_correlacao(mtcars, top = 3, filtrar_matriz = TRUE)

# Seleção manual de variáveis específicas
matriz_correlacao(mtcars, variaveis = c("mpg", "wt", "hp", "disp"), metodo = "kendall")
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
<summary><strong> Estatística Multivariada</strong></summary>
<br>

**Análise de Componentes Principais (ACP / PCA)**

```r
# Ajuste com diagnóstico completo (KMO, Bartlett, Kaiser, Horn)
modelo_acp <- acp(mtcars)

# Acesso direto aos elementos do modelo
modelo_acp$autovalores     # Variância explicada e acumulada por componente
modelo_acp$cargas          # Cargas dos componentes retidos (loadings)
modelo_acp$escores         # Coordenadas/escores das observações
modelo_acp$kmo             # KMO global e individual por variável
modelo_acp$bartlett        # Teste de Esfericidade de Bartlett
modelo_acp$criterios       # Comparação de retenção (Kaiser vs Análise Paralela)

# Visualizações gráficas dedicadas
biplot_acp(modelo_acp)
screeplot_acp(modelo_acp)

# Biplot com grupos e elipses de confiança
biplot_acp(modelo_acp, grupo = mtcars$cyl, elipses = TRUE)
```

**Análise Fatorial Exploratória (AFE)**

```r
# Ajuste fatorial com rotação Varimax e método dos Eixos Principais
modelo_afe <- analise_fatorial(mtcars, rotacao = "varimax", metodo = "pa")

# Acesso direto aos elementos do modelo
modelo_afe$cargas          # Cargas fatoriais, Comunalidades (h²) e Unicidades (u²)
modelo_afe$comunalidades   # Vetor de comunalidades por variável
modelo_afe$unicidades      # Vetor de unicidades por variável
modelo_afe$variancia       # Variância explicada por fator rotacionado
modelo_afe$escores         # Escores fatoriais individuais
modelo_afe$kmo             # KMO global e individual por variável
modelo_afe$bartlett        # Teste de Esfericidade de Bartlett

# Gráficos dedicados
cargasplot_fatorial(modelo_afe, corte = 0.30)   # Gráfico de barras das cargas
screeplot_fatorial(modelo_afe)                  # Scree plot com Horn
biplot_fatorial(modelo_afe, grupo = mtcars$cyl) # Biplot 2D com escores e cargas
```

**Análise de Agrupamentos (Clustering: K-Means e Hierárquico)**

```r
# K-Means com otimização automática de nstart e busca de k ideal
modelo_km <- cluster_kmeans(mtcars, k = 3)

# Acesso aos elementos do modelo K-Means
modelo_km$clusters           # Atribuição de cada observação ao cluster
modelo_km$dados_cluster      # Dados originais com a coluna $Cluster anexada
modelo_km$centroides         # Médias das variáveis por cluster (escala original)
modelo_km$resumo             # Tamanho (n), proporção (%) e silhueta de cada grupo
modelo_km$criterios_k        # WSS, Silhueta e Calinski-Harabasz de k = 2 a k_max
modelo_km$diagnostico_nstart # Teste de estabilização de variância do nstart

# Gráficos do K-Means
plot_clusters(modelo_km, elipses = TRUE) # Projeção 2D com centroidas
cotovelo_plot(modelo_km)                 # Método do Cotovelo e Silhueta por k

# Cluster Hierárquico com corte em k e Dendrograma
modelo_hc <- cluster_hierarquico(mtcars, k = 3, metodo_ligacao = "ward.D2")
dendrograma(modelo_hc)                   # Dendrograma colorido
```

**Análise de Variância Multivariada (MANOVA)**

```r
# Ajuste da MANOVA
modelo_manova <- analise_manova(c(Sepal.Length, Sepal.Width, Petal.Length), grupo = Species, dados = iris)

# Acesso direto aos elementos do modelo
modelo_manova$box_m               # Teste M de Box (Homogeneidade de Covariâncias)
modelo_manova$tabela_multivariada # 4 Critérios (Pillai, Wilks, Hotelling-Lawley, Roy)
modelo_manova$anovas_univariadas  # ANOVAs univariadas com correção de Bonferroni
modelo_manova$medias_grupos       # Médias observadas por grupo
modelo_manova$modelo_manova       # Objeto nativo manova do R

# Gráfico de Perfis Médios Multivariados (Escores Z padronizados)
plot_manova(modelo_manova)
```

**Análise Discriminante Linear (LDA)**

```r
# Sintaxe
modelo_lda <- analise_discriminante(variaveis = c("Sepal.Length", "Sepal.Width"), grupo = "Species", dados = iris)

# Acesso direto aos elementos do modelo
modelo_lda$funcoes_autovalores   # Autovalores, variância explicada (%) e correlação canônica
modelo_lda$scaling               # Coeficientes das funções discriminantes lineares
modelo_lda$matriz_confusao       # Matriz de confusão (resubstituição)
modelo_lda$metricas              # Acurácia, Erro, Kappa, Sensibilidade, Precisão e F1-Score (Treino vs CV)
modelo_lda$box_m                 # Teste M de Box (Homogeneidade de Covariâncias)
modelo_lda$escores               # Coordenadas/escores individuais projetados nas funções discriminantes
modelo_lda$modelo_lda            # Objeto nativo MASS::lda

# Biplot Discriminante com elipses e vetores direcionais (tema 2 com grade dupla)
plot_discriminante(modelo_lda)
```

</details>

---

<details>
<summary><strong> Planejamento e Análise de Experimentos</strong></summary>
<br>

**Delineamento Inteiramente Casualizado (DIC)**

```r
# Sintaxe
modelo_dic <- analise_dic("Produtividade", "Adubo", dados = soja, teste = "scott_knott")

# Teste de Scheffé com contrastes em linguagem natural:
analise_dic(Produtividade, Adubo, dados = soja, teste = "scheffe",
    contrastes = c("A vs c(B, C)", "c(A, B) vs c(C, D)"))

# Teste de Dunnett contra testemunha/controle:
analise_dic(Produtividade, Adubo, dados = soja, teste = "dunnett", controle = "Testemunha")

# Fator Quantitativo (Doses/Níveis numéricos — Regressão Polinomial com Ponto Ótimo):
modelo_quanti <- analise_dic(Produtividade, Dose_N, dados = milho)

# Acesso direto aos elementos do modelo
modelo_dic$anova                 # Quadro completo da ANOVA (Tratamento, Resíduo, Total)
modelo_dic$parametros            # Média geral, CV experimental (%) e R² do modelo
modelo_dic$tabela_grupos         # Médias, desvios padrões e letras do agrupamento (a, b, ab)
modelo_dic$comparacoes           # Diferenças pareadas com intervalos de confiança e p-valores
modelo_dic$contrastes            # Estimativas e testes dos contrastes de Scheffé
modelo_dic$regressao             # ANOVA polinomial (Linear, Quadrático, Cúbico) e ponto ótimo
modelo_dic$modelo_aov            # Objeto nativo aov do R

# Visualizações gráficas dedicadas (tema 2 com grade dupla)
grafico_dic(modelo_dic)          # Painel lado a lado: Médias com letras + ICs das diferenças
grafico_medias(modelo_dic)       # Gráfico de barras com erro e letras do teste sobrepostas
grafico_ic(modelo_dic)           # Gráfico horizontal de ICs com linha vermelha no zero e destaque
grafico_regressao(modelo_quanti) # Curva suave ajustada com pontos reais e dose ótima destacada
```

**Delineamento em Blocos Casualizados (DBC)**

```r
# Sintaxe
modelo_dbc <- analise_dbc("Produtividade", "Tratamento", "Bloco", dados = dados, teste = "scheffe",
                  contrastes = c("T1 vs c(T2, T3)", "c(T1, T2) vs c(T3, T4)"))

# Fator Quantitativo com Blocos (Regressão Polinomial controlada por blocos):
modelo_dbc_quanti <- analise_dbc(Produtividade, Dose, Bloco, dados = dados_quanti)

# Acesso direto aos elementos do modelo DBC
modelo_dbc$anova                 # Quadro da ANOVA (Bloco, Tratamento, Resíduo, Total)
modelo_dbc$parametros            # Média geral, CV experimental (%), R² e Eficiência da Blocagem (%)
modelo_dbc$tabela_grupos         # Médias com letras do teste de agrupamento
modelo_dbc$comparacoes           # Intervalos de confiança das diferenças pareadas
modelo_dbc$regressao             # Regressão polinomial controlada por blocos e dose ótima
modelo_dbc$modelo_aov            # Objeto aov nativo

# Visualização gráfica dedicada em painel
grafico_dbc(modelo_dbc)          # Painel automático de Médias + Intervalos de Confiança
```

**Delineamento em Quadrado Latino (DQL)**

```r
# Sintaxe
modelo_dql <- analise_dql("Produtividade", "Tratamento", "Linha", "Coluna", dados = dados, teste = "scott_knott")

# Teste de Scheffé com contrastes personalizados:
analise_dql(Produtividade, Tratamento, Linha, Coluna, dados = dados, teste = "scheffe",
    contrastes = c("A vs c(B, C)", "c(A, B) vs c(C, D)"))

# Acesso direto aos elementos do modelo DQL
modelo_dql$anova                 # Quadro da ANOVA (Linha, Coluna, Tratamento, Resíduo, Total)
modelo_dql$parametros            # Média geral, CV (%), R² e Eficiências Relativas (vs Linhas, vs Colunas, vs DIC)
modelo_dql$tabela_grupos         # Médias observadas com letras do teste de médias
modelo_dql$comparacoes           # Intervalos de confiança das diferenças pareadas

# Visualização gráfica dedicada em painel
grafico_dql(modelo_dql)          # Painel com Gráfico de Médias + ICs das diferenças
```

**Esquema Fatorial Duplo (A × B — DIC e DBC)**

```r
# Sintaxe - DIC
modelo_fat <- analise_fatorial_duplo(Produtividade, Cultivar, Adubacao, dados = dados, teste = "tukey")

# Em DBC
modelo_fat_dbc <- analise_fatorial_duplo(Produtividade, Cultivar, Adubacao, bloco = Bloco, dados = dados)

# Desdobramento automático:
# - Se a interação for significativa (p < 0.05): realiza o desdobramento dos fatores
#   e exibe a tabela matricial combinada de médias com letras minúsculas (na linha) e maiúsculas (na coluna).
# - Se a interação não for significativa: avalia e apresenta os efeitos principais independentes.

# Acesso direto aos elementos do modelo
modelo_fat$anova                   # Quadro da ANOVA com interação formatada (ex: Cultivar × Adubacao)
modelo_fat$parametros              # Média geral, CV experimental (%) e R² do modelo
modelo_fat$medias_comb             # Matriz de médias por combinação dos níveis
modelo_fat$letras_min              # Letras minúsculas do desdobramento (linha: Fator 1 dentro de Fator 2)
modelo_fat$letras_mai              # Letras maiúsculas do desdobramento (coluna: Fator 2 dentro de Fator 1)
modelo_fat$desdobramento_f1_em_f2  # Quadro de ANOVA do Fator 1 desdobrado em cada nível do Fator 2
modelo_fat$desdobramento_f2_em_f1  # Quadro de ANOVA do Fator 2 desdobrado em cada nível do Fator 1
modelo_fat$modelo_aov              # Objeto aov nativo

# Visualização gráfica do perfil de interação
grafico_fatorial(modelo_fat)       # Perfil de interação com barras de erro e letras sobrepostas
```

**Parcelas Subdivididas (Split-Plot — DIC e DBC)**

```r
# Em DBC (com Bloco):
modelo_psub <- analise_parcela_subdividida(Produtividade, Parcela = Manejo, Subparcela = Adubacao, bloco = Bloco, dados = dados)

# Em DIC (sem Bloco):
modelo_psub_dic <- analise_parcela_subdividida(Produtividade, Parcela = Manejo, Subparcela = Adubacao, dados = dados)

# Acesso aos resultados
modelo_psub$anova          # ANOVA com Erro (a) para Parcelas e Erro (b) para Subparcelas
modelo_psub$parametros     # CV(a) das parcelas, CV(b) das subparcelas e QM combinado
modelo_psub$letras_min     # Letras minúsculas: Subparcela dentro de Parcela (Erro b)
modelo_psub$letras_mai     # Letras maiúsculas: Parcela dentro de Subparcela (Erro combinado Satterthwaite)

# Diagnóstico de resíduos integrado
analise_residual(modelo_psub)
```

**Esquema Fatorial com Testemunha Adicional (A × B + 1)**

```r
# Fatorial 2x2 + 1 Testemunha isolada
modelo_fat_ad <- analise_fatorial_adicional(
  resposta   = Produtividade,
  fator1     = Cultivar,
  fator2     = Adubo,
  tratamento = Tratamento,
  testemunha = "Testemunha",
  bloco      = Bloco,
  dados      = dados
)

# Acesso aos resultados
modelo_fat_ad$anova        # Decomposição ortogonal: Fatorial, F1, F2, F1xF2 e Fatorial vs Testemunha
modelo_fat_ad$dunnett      # Teste de Dunnett comparando cada combinação diretamente contra a Testemunha

# Diagnóstico de resíduos integrado
analise_residual(modelo_fat_ad)
```

</details>

---

<details>
<summary><strong> Análise de Sobrevivência</strong></summary>
<br>

**Estimador de Kaplan-Meier**

```r
# Curva de Kaplan-Meier com grupos e teste de Log-Rank:
modelo_km <- analise_sobrevivencia(tempo, evento, grupo = tratamento, dados = dados)

# Com estratificação (Log-Rank Estratificado para controle de variáveis):
modelo_km_est <- analise_sobrevivencia(tempo, evento, grupo = tratamento, estrato = sexo, dados = dados)

# Estimador de Nelson-Aalen e Fleming-Harrington:
modelo_na <- analise_sobrevivencia(tempo, evento, grupo = tratamento, dados = dados, metodo = "nelson-aalen")

# Curva global (sem grupos):
km_global <- analise_sobrevivencia(tempo, evento, dados = dados)

# Acesso direto aos resultados:
modelo_km$tabela           # Resumo por grupo: Total, Eventos, Censuras (%), Mediana e IC 95%
modelo_km$curva            # Tabela pontual completa: S(t), Risco Acumulado H(t), Taxa de Risco h(t) e Densidade f(t)
modelo_km$teste_logrank    # Estatística Chi2, GL, p-valor (com indicador de estratificação)
modelo_km$fit_km           # Objeto nativo survfit do pacote survival

# Gráficos Analíticos com Grade Dupla (4 funções da Sobrevivência):
grafico_sobrevivencia(modelo_km, tipo = "sobrevivencia")     # Curva clássica S(t) com Risk Table
grafico_sobrevivencia(modelo_km, tipo = "risco_acumulado")   # Risco acumulado H(t) de Nelson-Aalen
grafico_sobrevivencia(modelo_km, tipo = "risco_instantaneo") # Taxa de risco h(t) com tendência LOESS
grafico_sobrevivencia(modelo_km, tipo = "densidade")         # Função de densidade de falhas f(t)
```

**Modelo de Riscos Proporcionais de Cox**

```r
# Ajuste do modelo semi-paramétrico de Cox
modelo_cox <- regressao_cox(Surv(tempo, evento) ~ idade + tratamento + estagio, dados = dados)

# Acesso direto aos resultados
modelo_cox$tabela          # Coeficientes (Beta), EP, Hazard Ratio (HR), IC 95%, z e p-valor
modelo_cox$c_index         # Índice de Concordância (C-index / AUC temporal) e Erro Padrão
modelo_cox$lrt             # Teste da Razão de Verossimilhança global
modelo_cox$wald            # Teste de Wald global
modelo_cox$diagnostico     # Teste de Schoenfeld para o pressuposto de Riscos Proporcionais

# Forest Plot dos Hazard Ratios com intervalos de confiança de 95%:
grafico_cox(modelo_cox)

# Métricas polimórficas (C-index, Pseudo-R2, AIC, BIC, LRT, Wald):
metricas(modelo_cox)

# Seleção polimórfica de modelos (Best Subsets por AIC e C-index):
selecao_modelos(Surv(tempo, evento) ~ idade + tratamento + estagio, dados = dados)

# Diagnóstico e Análise Residual integrada (Schoenfeld, Deviance e Martingale):
analise_residual(modelo_cox)
```

</details>

---

<details>
<summary><strong> Amostragem</strong></summary>
<br>

**Tamanho ideal da amostra**

```r
# Tamanho de amostra
tamanho_amostra(populacao = 50000, erro = 0.03)
tamanho_amostra(populacao = Inf, erro = 0.05, confianca = 0.99)
```
**Amostragem Aleatória Simples**

```r
# AAS
amostra_aleatoria(mtcars, n = 10)
amostra_aleatoria(mtcars, proporcao = 0.3, semente = 42)
```

**Amostragem Sistemática**

```r
# AS
amostra_sistematica(mtcars, n = 10)
```

**Amostragem Estratificada**

```r
# AE
amostra_estratificada(iris, estrato = Species, n = 30)
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "uniforme")
amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "otima", variavel = Sepal.Length)
```

**Amostragem por Conglomerados**

```r
# AC
amostra_conglomerados(mtcars, conglomerado = cyl, n_conglomerados = 2)
```

</details>

---

<details>
<summary><strong> Visualização e Gráficos</strong></summary>
<br>

**Gráficos**

```r
# Boxplot básico
grafico_boxplot(mtcars, cyl, mpg)
grafico_boxplot(mtcars, cyl, mpg, violino = TRUE)

# Boxplot agrupado com dispersão e paleta estatR
grafico_boxplot(mtcars, cyl, mpg, grupo = am,
                dispersao_pts = TRUE, paleta = "vibrant")

# Destaque de categorias e ordenamento por mediana
grafico_boxplot(mtcars, cyl, mpg, ordenar = TRUE, destaque = "8")

# Forest Plot de coeficientes de modelos
grafico_coeficientes(modelo)
```

**Tema e Paletas para ggplot2**

```r
# Aplicar tema padronizado
ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
  geom_boxplot() +
  tema_estatR(estilo = 2, grade = "dupla") + # 1 = classic, 2 = minimal com grade, 3 = limpo
  paleta_estatR("academic")                 # Paleta de cores padronizada

# Outros ajustes de tema
tema_estatR(modo = "dark")                  # Tema escuro
tema_estatR(inclinar = 45)                  # Inclina rótulos do eixo X
tema_estatR(fundo_transparente = TRUE)      # Fundo 100% transparente (para slides/apresentações)

# Visualizar paletas disponíveis
paleta_estatR()                             # Painel com todas as paletas
paleta_estatR("vibrant")                    # Ver uma paleta específica
```

</details>

---

<details>
<summary><strong> Séries Temporais</strong></summary>
<br>

**Diagnóstico e Decomposição**

```r
# Análise completa da série: estatísticas, força dos componentes e painel gráfico (4 painéis)
analise_serie(AirPassengers, tipo = "multiplicativa")
analise_serie(dados, y = vendas, data = mes, frequencia = 12)
```

**Estacionariedade**

```r
# Testes ADF e KPSS com veredito
teste_estacionariedade(AirPassengers)
teste_estacionariedade(AirPassengers, diferencas = 1)
```

**Autocorrelação (FAC e FACP)**

```r
# Gráfico de barras (lollipop) com faixas de Bartlett + Ljung-Box
grafico_autocorrelacao(AirPassengers, defasagens = 24)
grafico_acf(AirPassengers)  # alias amigável
```

**Previsão com Suavização Exponencial (Fan Chart)**

```r
# Seleção automática (menor AICc) ou método explícito
previsao_serie(AirPassengers, horizonte = 12, metodo = "auto")
previsao_serie(AirPassengers, horizonte = 12, metodo = "holt_winters", tipo = "multiplicativa")
```

**Métricas de Acurácia**

```r
# RMSE, MAE, MAPE, sMAPE, MASE, UMBRAE, AIC, AICc, BIC
metricas_serie(modelo_previsao)          # Objeto de previsão
metricas(modelo_previsao)           
```

**Diagnóstico de Resíduos (Polimórfico)**

```r
# Testes formais (Média Zero, Normalidade, Ljung-Box) + Painel com 4 gráficos
analise_residual(modelo_previsao)             # Objeto de previsão estatR
analise_residual(modelo_arima)                # Modelos ARIMA / SARIMA
```

**Seleção de Modelos (Out-of-Sample Ranking)**

```r
# Compara famílias de modelos por desempenho fora da amostra (holdout)
selecao_modelos_serie(AirPassengers, h = 12)
selecao_modelos_serie(AirPassengers, familia = "suavizacao", top = 5, ordenar_por = "sMAPE")
selecao_modelos_serie(AirPassengers, familia = c("estocastica"))  # ARIMA/SARIMA

# Famílias disponíveis: "todas", "suavizacao", "estocastica", "decomposicao"
```

</details>

---

<details>
<summary><strong> Estatística Não-Paramétrica e Concordância</strong></summary>
<br>

**Comparações de Múltiplos Grupos Não-Paramétricos**

```r
# Kruskal-Wallis com pós-teste automático de Dunn e ajuste de p-valor (Holm/Bonferroni)
teste_kruskal(airquality, Ozone ~ Month)
teste_kruskal(dados, y = rendimento, grupo = tratamento)

# Teste de Friedman para medidas repetidas / blocos casualizados com Wilcoxon pareado
teste_friedman(dados, valor ~ tempo | sujeito)
teste_friedman(dados, y = resposta, grupo = periodo, bloco = paciente)
```

**Concordância e Confiabilidade entre Avaliadores**

```r
# Kappa de Cohen (2 avaliadores - simples, linear ou quadrático)
concordancia_kappa(dados, avaliador1 = Medico_A, avaliador2 = Medico_B)
concordancia_kappa(dados, Medico_A, Medico_B, ponderado = "quadratico")

# Kappa de Fleiss (3 ou mais avaliadores)
concordancia_kappa(matriz_avaliadores)  # Matriz de sujeitos x categorias
```

**Concordância Contínua de Métodos (Bland-Altman)**

```r
# Viés médio, limites de concordância de 95%, ICs e gráfico com grade dupla
grafico_bland_altman(dados, metodo1 = Aparelho_A, metodo2 = Aparelho_B)
```

</details>

---

<details>
<summary><strong> Aprendizado de Máquina (Machine Learning)</strong></summary>
<br>

**Previsão Numérica (Regressão)**

```r
# Random Forest com auto-tuning de árvores por estabilização OOB
mod_prev <- previsao_aprend_maq(mpg ~ wt + hp + qsec, dados = mtcars, metodo = "rf")

# Parametrização manual em português
mod_prev <- previsao_aprend_maq(mpg ~ wt + hp, dados = mtcars, metodo = "rf",
                                n_arvores = 500, vars_por_no = 2)

# Outros algoritmos suportados: "xgboost", "svm", "arvore"
mod_xgb <- previsao_aprend_maq(mpg ~ ., dados = mtcars, metodo = "xgboost",
                               taxa_aprendizado = 0.05, profundidade = 4)

# Métricas de desempenho (Treino vs Teste)
metricas(mod_prev)

# Diagnóstico de resíduos formal (viés sistemático, normalidade e outliers)
analise_residual(mod_prev)
```

**Classificação Categórica e Binária**

```r
# Classificação com Random Forest e matriz de confusão
mod_cls <- classificacao_aprend_maq(Species ~ ., dados = iris, metodo = "rf")

# Métodos disponíveis: "rf", "xgboost", "svm", "knn", "arvore"
mod_knn <- classificacao_aprend_maq(Species ~ ., dados = iris, metodo = "knn", k = 5)

# Métricas polimórficas (Acurácia, Kappa de Cohen)
metricas(mod_cls)
```

**Seleção e Comparação de Modelos (Estatística vs Aprendizado de Máquina)**

```r
# Compara modelos estatísticos clássicos (OLS/GLM) com modelos de ML sob Validação Cruzada (5-fold CV)
selecao_modelos(mpg ~ wt + hp + qsec, dados = mtcars,
                comparar_modelos = c("estatistico", "aprend_maq"), top = 5)
```

</details>

---

## Autor

**Elizeu Medeiros**  
[GitHub](https://github.com/elizeumedeirosest)

---

## Licença

MIT © Elizeu Medeiros
