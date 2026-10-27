#' estatR: Uma Camada Estatística Intuitiva e Didática para o R
#'
#' O \code{estatR} oferece uma interface amigável e didática sobre as principais
#' análises estatísticas do R. Pensado para ensino, pesquisa acadêmica e ciência de dados, ele
#' produz saídas no console com formatação profissional em português, gráficos
#' editoriais automáticos com grade dupla e retorno estruturado para manipulação direta.
#'
#' @section Módulos Disponíveis (14 Módulos Especializados):
#'
#' \strong{1. 📊 Diagnóstico e Estatística Descritiva}
#' \itemize{
#'   \item \code{\link{diagnostico}} — Visão geral de integridade, tipos e anomalias de um banco de dados.
#'   \item \code{\link{descrever}} — Estatísticas descritivas completas e tabelas de frequência.
#'   \item \code{\link{tabela_contingencia}} — Tabela cruzada entre duas variáveis categóricas.
#'   \item \code{\link{med_tend_central}}, \code{\link{med_dispersao}}, \code{\link{med_forma}} — Medidas de síntese numérica.
#'   \item \code{\link{quartis}}, \code{\link{quintis}}, \code{\link{decis}}, \code{\link{percentis}} — Separatrizes.
#'   \item \code{\link{outliers}} — Detecção e identificação de anomalias pelo método IQR.
#' }
#'
#' \strong{2. 🎯 Testes de Hipótese e Inferência}
#' \itemize{
#'   \item \code{\link{teste_t_uma_amostra}}, \code{\link{teste_t_duas_amostras}}, \code{\link{teste_t_pareado}} — Testes t de Student e Welch com tamanho de efeito (d de Cohen).
#'   \item \code{\link{teste_proporcao}} — Teste Z de proporção amostral.
#'   \item \code{\link{teste_equivalencia}} — Procedimento de Dois Testes Unicaudais (TOST) para equivalência clínica/estatística.
#'   \item \code{\link{teste_mcnemar}} — Teste de McNemar para dados categóricos binários pareados com Odds Ratio pareado.
#'   \item \code{\link{teste_wilcoxon}}, \code{\link{teste_mann_whitney}} — Alternativas não-paramétricas para 1 e 2 amostras.
#'   \item \code{\link{teste_kruskal}}, \code{\link{teste_friedman}} — Comparações não-paramétricas para múltiplos grupos ou medidas repetidas.
#'   \item \code{\link{teste_qui_quadrado}} — Qui-quadrado de independência com tamanho de efeito (V de Cramér).
#'   \item \code{\link{teste_levene}} — Teste de homogeneidade de variâncias.
#'   \item \code{\link{concordancia_kappa}}, \code{\link{grafico_bland_altman}} — Confiabilidade entre avaliadores e concordância contínua.
#' }
#'
#' \strong{3. 📏 Intervalos de Confiança}
#' \itemize{
#'   \item \code{\link{ic_media}}, \code{\link{ic_proporcao}}, \code{\link{ic_variancia}}, \code{\link{ic_diferenca_medias}} — Estimativas intervalares formais.
#' }
#'
#' \strong{4. 🎲 Probabilidade e Simulação}
#' \itemize{
#'   \item \code{\link{distribuicao}} — Painel completo de distribuições (Normal, t, F, Chi², Poisson, Binomial, Weibull, Gama, etc.).
#'   \item \code{\link{prob_normal}}, \code{\link{prob_binomial}}, \code{\link{prob_poisson}}, \code{\link{prob_exponencial}}, \code{\link{prob_weibull}} — Cálculo exato de probabilidades com curvas sombreadas.
#'   \item \code{\link{gerar_amostra}} — Geração de amostras aleatórias teóricas.
#'   \item \code{\link{estimar_verossimilhanca}}, \code{\link{estimar_momentos}} — Estimação pontual com gráficos de suporte.
#'   \item \code{\link{bootstrap}}, \code{\link{simular_tcl}}, \code{\link{monte_carlo}} — Reamostragem computacional e simulação.
#' }
#'
#' \strong{5. 🎛️ Transformações de Variáveis}
#' \itemize{
#'   \item \code{\link{transformar}} — Família de potências de Box-Cox com intervalo de confiança de lambda, Yeo-Johnson, Log, Raiz, Logit, Arco-Seno, Padronização e seleção automática ("auto") por Shapiro-Wilk.
#' }
#'
#' \strong{6. 📈 Modelos Lineares e Regressão}
#' \itemize{
#'   \item \code{\link{regressao_linear}} — Regressão linear simples e múltipla com ANOVA, coeficientes, VIF e interpretação.
#'   \item \code{\link{coeficientes}}, \code{\link{residuos}}, \code{\link{valores_ajustados}}, \code{\link{resumo}} — Extratores universais de modelos em português.
#'   \item \code{\link{analise_residual}} — Diagnóstico residual formal (Normalidade, Homocedasticidade, Independência).
#'   \item \code{\link{metricas}}, \code{\link{selecao_modelos}}, \code{\link{predicao}} — Avaliação, Best Subsets e intervalos de predição.
#'   \item \code{\link{multicolinearidade}}, \code{\link{pontos_influentes}}, \code{\link{grafico_coeficientes}} — Diagnóstico de alavancagem e Forest Plot.
#' }
#'
#' \strong{7. 📊 Modelos Lineares Generalizados (MLG)}
#' \itemize{
#'   \item \code{\link{regressao_mlg}} — Ajuste flexível para Binomial, Poisson, Quasipoisson, Gama, Gaussiana, Binomial Negativa e Beta.
#'   \item \code{\link{regressao_logistica}}, \code{\link{odds_ratio}} — Modelos logísticos (binário, multinomial, ordinal) e Razões de Chances.
#'   \item \code{\link{regressao_poisson}}, \code{\link{regressao_binomial_negativa}}, \code{\link{irr}} — Modelos de contagem e Razões de Taxas de Incidência.
#'   \item \code{\link{regressao_gama}}, \code{\link{regressao_beta}} — Modelos para dados assimétricos estritos e taxas no intervalo (0, 1).
#' }
#'
#' \strong{8. 📉 Séries Temporais}
#' \itemize{
#'   \item \code{\link{analise_serie}}, \code{\link{teste_estacionariedade}} — Decomposição temporal e testes ADF/KPSS.
#'   \item \code{\link{previsao_serie}}, \code{\link{previsao_serie_aprend_maq}} — Suavização exponencial e previsão com Machine Learning.
#'   \item \code{\link{grafico_autocorrelacao}}, \code{\link{selecao_modelos_serie}} — FAC/FACP e ranking fora da amostra.
#' }
#'
#' \strong{9. 🧬 Análise Multivariada e Agrupamento}
#' \itemize{
#'   \item \code{\link{acp}}, \code{\link{biplot_acp}}, \code{\link{screeplot_acp}} — Análise de Componentes Principais com KMO e Bartlett.
#'   \item \code{\link{analise_fatorial}}, \code{\link{biplot_fatorial}} — Análise Fatorial Exploratória (AFE).
#'   \item \code{\link{cluster_kmeans}}, \code{\link{cluster_hierarquico}}, \code{\link{dendrograma}} — Agrupamentos com silhueta e cotovelo.
#'   \item \code{\link{analise_discriminante}}, \code{\link{analise_manova}} — LDA/QDA e MANOVA com teste M de Box.
#' }
#'
#' \strong{10. 🧪 Delineamentos Experimentais}
#' \itemize{
#'   \item \code{\link{analise_dic}}, \code{\link{analise_dbc}}, \code{\link{analise_dql}} — Delineamentos clássicos com Tukey, Scott-Knott, Dunnett e regressão polinomial.
#'   \item \code{\link{analise_fatorial_duplo}}, \code{\link{analise_fatorial_adicional}}, \code{\link{analise_parcela_subdividida}} — Esquemas fatoriais e Split-Plot.
#'   \item \code{\link{analise_medidas_repetidas}} — ANOVA de medidas repetidas com teste de Mauchly e correções GG/HF.
#'   \item \code{\link{teste_dunnett}}, \code{\link{grafico_dunnett}} — Comparações múltiplas contra controle/testemunha.
#' }
#'
#' \strong{11. 🛡️ Controle Estatístico da Qualidade (CEQ & Six Sigma)}
#' \itemize{
#'   \item \code{\link{grafico_controle}} — Cartas de Shewhart (I-MR, Xbar-R, Xbar-S, p, np, c, u) com regras de Western Electric / Nelson.
#'   \item \code{\link{capacidade_processo}} — Análise de capacidade e desempenho (Cp, Cpk, Cpm, Pp, Ppk, PPM e nível Sigma).
#' }
#'
#' \strong{12. 🤖 Aprendizado de Máquina (Machine Learning)}
#' \itemize{
#'   \item \code{\link{classificacao_aprend_maq}}, \code{\link{previsao_aprend_maq}} — Treinamento supervisionado com Random Forest e métricas completas.
#'   \item \code{\link{validacao_cruzada}} — Validação cruzada K-Fold para regressão e classificação com avaliação fora da amostra.
#'   \item \code{\link{analise_residual_aprend_maq}} — Diagnóstico residual para modelos de Machine Learning.
#' }
#'
#' \strong{13. ⏱️ Análise de Sobrevivência}
#' \itemize{
#'   \item \code{\link{analise_sobrevivencia}}, \code{\link{grafico_sobrevivencia}} — Curvas de Kaplan-Meier, risco de Nelson-Aalen e teste de Log-Rank.
#'   \item \code{\link{regressao_cox}}, \code{\link{metricas_cox}}, \code{\link{selecao_modelos_cox}}, \code{\link{grafico_cox}} — Modelo de Cox, resíduos de Schoenfeld e Forest Plot de Hazard Ratios.
#' }
#'
#' \strong{14. 👥 Amostragem Probabilística}
#' \itemize{
#'   \item \code{\link{amostra_aleatoria}}, \code{\link{amostra_estratificada}}, \code{\link{amostra_sistematica}}, \code{\link{amostra_conglomerados}} — Planos amostrais probabilísticos.
#'   \item \code{\link{tamanho_amostra}} — Dimensionamento ótimo de amostra.
#'   \item \code{\link{poder_estatistico}} — Cálculo de poder (1 - beta) e curvas de poder para testes estatísticos.
#' }
#'
#' @author Elizeu Medeiros \email{elizeumedeirosest@@github.com}
#' @docType package
#' @name estatR
"_PACKAGE"
