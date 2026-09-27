#' @title Probabilidade da Distribuição Normal
#' @description Calcula a probabilidade para uma distribuição Normal e gera o gráfico com a área sombreada.
#' @param media Média teórica da distribuição.
#' @param dp Desvio-padrão teórico.
#' @param q1 Quantil de referência (limite de integração).
#' @param q2 Segundo quantil de referência (apenas para tipo = "entre").
#' @param tipo Tipo de probabilidade desejada: "menor" (P(X < q1)), "maior" (P(X > q1)) ou "entre" (P(q1 < X < q2)).
#' @param grafico Lógico. Se TRUE, gera o gráfico da distribuição com a área hachurada.
#' @export
prob_normal <- function(media = 0, dp = 1, q1 = NULL, q2 = NULL, tipo = c("menor", "maior", "entre"), grafico = TRUE) {
  tipo <- match.arg(tipo)
  
  if (is.null(q1)) stop("Você precisa informar o 'q1'.")
  if (dp <= 0) stop("O desvio-padrão (dp) deve ser maior que zero.")
  
  if (tipo == "entre") {
    if (is.null(q2)) stop("Para o tipo 'entre', você precisa informar o 'q2'.")
    if (q1 > q2) {
      tmp <- q1
      q1 <- q2
      q2 <- tmp
    }
  }
  
  # Cálculos estatísticos
  if (tipo == "menor") {
    prob <- pnorm(q1, mean = media, sd = dp)
    z1 <- (q1 - media) / dp
    inf_lim <- media - 4 * dp
    sup_lim <- q1
    titulo <- sprintf("P(X < %g)", q1)
  } else if (tipo == "maior") {
    prob <- 1 - pnorm(q1, mean = media, sd = dp)
    z1 <- (q1 - media) / dp
    inf_lim <- q1
    sup_lim <- media + 4 * dp
    titulo <- sprintf("P(X > %g)", q1)
  } else {
    prob <- pnorm(q2, mean = media, sd = dp) - pnorm(q1, mean = media, sd = dp)
    z1 <- (q1 - media) / dp
    z2 <- (q2 - media) / dp
    inf_lim <- q1
    sup_lim <- q2
    titulo <- sprintf("P(%g < X < %g)", q1, q2)
  }
  
  linha <- paste0("\n", strrep("\u2500", 54), "\n\n")
  
  # Saída Console
  cat("\nCÁLCULO DE PROBABILIDADE (Distribuição Normal)\n\n")
  cat(sprintf("Média (\u03bc):          %g\n", media))
  cat(sprintf("Desvio Padrão (\u03c3):  %g\n", dp))
  cat(linha)
  cat("RESULTADO\n\n")
  
  if (tipo == "entre") {
    cat(sprintf("Quantis (X):       %g e %g\n", q1, q2))
    cat(sprintf("Escores (Z):       %.2f e %.2f\n", z1, z2))
  } else {
    cat(sprintf("Quantil (X):       %g\n", q1))
    cat(sprintf("Escore (Z):        %.2f\n", z1))
  }
  
  cat(sprintf("Probabilidade:     %.2f%%\n", prob * 100))
  cat(linha)
  
  # Gráfico
  if (grafico) {
    df_curve <- data.frame(x = seq(media - 4*dp, media + 4*dp, length.out = 500))
    df_curve$y <- dnorm(df_curve$x, mean = media, sd = dp)
    df_poly <- df_curve[df_curve$x >= inf_lim & df_curve$x <= sup_lim, ]
    
    # Garantir que o polígono desça até o zero no eixo Y
    if (nrow(df_poly) > 0) {
      df_poly <- rbind(
        data.frame(x = df_poly$x[1], y = 0),
        df_poly,
        data.frame(x = df_poly$x[nrow(df_poly)], y = 0)
      )
    }
    
    p <- ggplot2::ggplot(df_curve, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(color = "#333333", linewidth = 1) +
      ggplot2::geom_polygon(data = df_poly, ggplot2::aes(x = x, y = y), fill = "dodgerblue1", alpha = 0.6) +
      ggplot2::labs(title = titulo, 
                    subtitle = sprintf("Probabilidade: %.2f%%", prob * 100),
                    x = "Valor (X)", y = "Densidade",
                    caption = "estatR")
    
    if (exists("meu_tema")) {
      p <- p + meu_tema(estilo = 1)
    } else {
      p <- p + ggplot2::theme_minimal()
    }
    
    suppressMessages(print(p))
  }
  
  invisible(prob)
}


#' @title Gerador de Amostras de Distribuições
#' @description Gera um conjunto de dados aleatórios baseado em distribuições teóricas e plota o histograma.
#' @param n Tamanho da amostra (quantidade de observações geradas).
#' @param distribuicao Tipo da distribuição ("normal", "binomial", "poisson", "exponencial").
#' @param media Parâmetro da Normal (média).
#' @param dp Parâmetro da Normal (desvio-padrão).
#' @param ensaios Parâmetro da Binomial (número de tentativas/ensaios).
#' @param prob Parâmetro da Binomial (probabilidade de sucesso).
#' @param lambda Parâmetro da Poisson (taxa média de ocorrência).
#' @param taxa Parâmetro da Exponencial (taxa lambda).
#' @param grafico Lógico. Se TRUE, exibe histograma dos dados gerados sobreposto com a curva teórica.
#' @param semente Semente opcional para reprodutibilidade (set.seed).
#' @export
gerar_amostra <- function(n, distribuicao = c("normal", "binomial", "poisson", "exponencial"),
                          media = 0, dp = 1, ensaios = 10, prob = 0.5, lambda = 1, taxa = 1,
                          grafico = TRUE, semente = NULL) {
  
  distribuicao <- match.arg(distribuicao)
  if (!is.null(semente)) set.seed(semente)
  
  if (n <= 0) stop("O tamanho da amostra (n) deve ser positivo.")
  
  # Geração e cálculos teóricos
  if (distribuicao == "normal") {
    amostra <- rnorm(n, mean = media, sd = dp)
    teo_media <- media
    teo_dp <- dp
    tit_graf <- sprintf("Amostra Normal (\u03bc=%g, \u03c3=%g)", media, dp)
  } else if (distribuicao == "binomial") {
    amostra <- rbinom(n, size = ensaios, prob = prob)
    teo_media <- ensaios * prob
    teo_dp <- sqrt(ensaios * prob * (1 - prob))
    tit_graf <- sprintf("Amostra Binomial (ensaios=%d, p=%g)", ensaios, prob)
  } else if (distribuicao == "poisson") {
    amostra <- rpois(n, lambda = lambda)
    teo_media <- lambda
    teo_dp <- sqrt(lambda)
    tit_graf <- sprintf("Amostra Poisson (\u03bb=%g)", lambda)
  } else if (distribuicao == "exponencial") {
    amostra <- rexp(n, rate = taxa)
    teo_media <- 1 / taxa
    teo_dp <- 1 / taxa
    tit_graf <- sprintf("Amostra Exponencial (\u03bb=%g)", taxa)
  }
  
  emp_media <- mean(amostra)
  emp_dp <- sd(amostra)
  
  linha <- paste0("\n", strrep("\u2500", 54), "\n\n")
  
  # Saída Console
  cat(sprintf("\nGERAÇÃO DE DADOS ALEATÓRIOS (%s)\n\n", tools::toTitleCase(distribuicao)))
  cat(sprintf("Tamanho da amostra (n): %d\n", n))
  cat(linha)
  cat(sprintf("%-20s %15s %15s\n", "Estatística", "Amostra Gerada", "Valor Teórico"))
  cat(paste(rep("-", 54), collapse = ""), "\n")
  cat(sprintf("%-20s %15.2f %15.2f\n", "Média", emp_media, teo_media))
  cat(sprintf("%-20s %15.2f %15.2f\n", "Desvio Padrão", emp_dp, teo_dp))
  cat(linha)
  
  if (grafico) {
    df <- data.frame(X = amostra)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = X))
    
    # Diferencia gráficos de distribuições contínuas e discretas
    if (distribuicao %in% c("normal", "exponencial")) {
      p <- p + 
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), 
                                bins = max(10, min(30, round(n/10))),
                                fill = "#DDDDDD", color = "black", alpha = 0.7)
      
      if (distribuicao == "normal") {
        p <- p + ggplot2::stat_function(fun = dnorm, args = list(mean = media, sd = dp), color = "#1F77B4", linewidth = 1.2)
      } else {
        p <- p + ggplot2::stat_function(fun = dexp, args = list(rate = taxa), color = "#1F77B4", linewidth = 1.2)
      }
    } else {
      # Binomial e Poisson (Discretas)
      p <- p + 
        ggplot2::geom_bar(ggplot2::aes(y = ggplot2::after_stat(prop)), 
                          fill = "#DDDDDD", color = "black", alpha = 0.7)
      
      # Cálculo manual das proporções teóricas para plotar pontos/linhas teóricas
      x_vals <- min(amostra):max(amostra)
      if (distribuicao == "binomial") {
        y_teo <- dbinom(x_vals, size = ensaios, prob = prob)
      } else {
        y_teo <- dpois(x_vals, lambda = lambda)
      }
      df_teo <- data.frame(x = x_vals, y = y_teo)
      
      p <- p + 
        ggplot2::geom_point(data = df_teo, ggplot2::aes(x = x, y = y), color = "#1F77B4", size = 3) +
        ggplot2::geom_segment(data = df_teo, ggplot2::aes(x = x, xend = x, y = 0, yend = y), color = "#1F77B4", linewidth = 1)
    }
    
    # Define o subtítulo com base no tipo de distribuição
    if (distribuicao %in% c("normal", "exponencial")) {
      texto_sub <- "Barras (Amostra Gerada) vs Curva (Distribuição Teórica)"
    } else {
      texto_sub <- "Barras (Amostra Gerada) vs Pontos (Distribuição Teórica)"
    }
    
    p <- p + ggplot2::labs(title = tit_graf,
                           subtitle = texto_sub,
                           y = ifelse(distribuicao %in% c("normal", "exponencial"), "Densidade", "Proporção"),
                           caption = "estatR")
                           
    if (exists("meu_tema")) {
      p <- p + meu_tema(estilo = 1)
    } else {
      p <- p + ggplot2::theme_minimal()
    }
    
    suppressMessages(print(p))
  }
  
  invisible(amostra)
}

#' @title Probabilidade da Distribuição Binomial
#' @description Calcula a probabilidade para uma distribuição Binomial e plota as barras.
#' @param ensaios Número de tentativas (n).
#' @param prob Probabilidade de sucesso (p).
#' @param q Valor de interesse (número de sucessos).
#' @param tipo "exato" (P(X=q)), "menor" (P(X<=q)), "maior" (P(X>=q)).
#' @param grafico Lógico. Se TRUE, gera o gráfico.
#' @export
prob_binomial <- function(ensaios, prob, q, tipo = c("exato", "menor", "maior"), grafico = TRUE) {
  tipo <- match.arg(tipo)
  if (q < 0 || q > ensaios) stop("'q' deve estar entre 0 e o número de ensaios.")
  
  if (tipo == "exato") {
    res <- dbinom(q, ensaios, prob)
    tit <- sprintf("P(X = %d)", q)
    destaque <- q
  } else if (tipo == "menor") {
    res <- pbinom(q, ensaios, prob)
    tit <- sprintf("P(X \u2264 %d)", q)
    destaque <- 0:q
  } else {
    res <- 1 - pbinom(q - 1, ensaios, prob)
    tit <- sprintf("P(X \u2265 %d)", q)
    destaque <- q:ensaios
  }
  
  linha <- paste0("\n", strrep("\u2500", 54), "\n\n")
  cat("\nCÁLCULO DE PROBABILIDADE (Distribuição Binomial)\n\n")
  cat(sprintf("Ensaios (n):           %d\n", ensaios))
  cat(sprintf("Probabilidade (p):     %g\n", prob))
  cat(linha)
  cat(sprintf("%-22s %d\n", ifelse(tipo=="exato", "Valor (X = x):", ifelse(tipo=="menor", "Valores (X \u2264 x):", "Valores (X \u2265 x):")), q))
  cat(sprintf("%-22s %.2f%%\n", "Probabilidade:", res * 100))
  cat(linha)
  
  if (grafico) {
    df <- data.frame(x = 0:ensaios)
    df$y <- dbinom(df$x, ensaios, prob)
    df$cor <- ifelse(df$x %in% destaque, "dodgerblue1", "#DDDDDD")
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = as.factor(x), y = y, fill = cor)) +
      ggplot2::geom_col(color = "black", alpha = 0.8) +
      ggplot2::scale_fill_identity() +
      ggplot2::labs(title = tit, subtitle = sprintf("Probabilidade: %.2f%%", res * 100),
                    x = "Número de Sucessos (X)", y = "Probabilidade", caption = "estatR")
    
    if (exists("meu_tema")) p <- p + meu_tema(estilo = 1) else p <- p + ggplot2::theme_minimal()
    suppressMessages(print(p))
  }
  invisible(res)
}

#' @title Probabilidade da Distribuição Poisson
#' @description Calcula a probabilidade para uma distribuição Poisson e plota as barras.
#' @param lambda Taxa de ocorrência.
#' @param q Valor de interesse.
#' @param tipo "exato" (P(X=q)), "menor" (P(X<=q)), "maior" (P(X>=q)).
#' @param grafico Lógico. Se TRUE, gera o gráfico.
#' @export
prob_poisson <- function(lambda, q, tipo = c("exato", "menor", "maior"), grafico = TRUE) {
  tipo <- match.arg(tipo)
  if (q < 0) stop("'q' deve ser \u2265 0.")
  
  if (tipo == "exato") {
    res <- dpois(q, lambda)
    tit <- sprintf("P(X = %d)", q)
    destaque <- q
  } else if (tipo == "menor") {
    res <- ppois(q, lambda)
    tit <- sprintf("P(X \u2264 %d)", q)
    destaque <- 0:q
  } else {
    res <- 1 - ppois(q - 1, lambda)
    tit <- sprintf("P(X \u2265 %d)", q)
    destaque <- q:max(q + 10, ceiling(lambda + 4*sqrt(lambda)))
  }
  
  linha <- paste0("\n", strrep("\u2500", 54), "\n\n")
  cat("\nCÁLCULO DE PROBABILIDADE (Distribuição Poisson)\n\n")
  cat(sprintf("Taxa (\u03bb):               %g\n", lambda))
  cat(linha)
  cat(sprintf("%-22s %d\n", ifelse(tipo=="exato", "Valor (X = x):", ifelse(tipo=="menor", "Valores (X \u2264 x):", "Valores (X \u2265 x):")), q))
  cat(sprintf("%-22s %.2f%%\n", "Probabilidade:", res * 100))
  cat(linha)
  
  if (grafico) {
    xmax <- max(q + 5, ceiling(lambda + 4*sqrt(lambda)))
    df <- data.frame(x = 0:xmax)
    df$y <- dpois(df$x, lambda)
    df$cor <- ifelse(df$x %in% destaque, "dodgerblue1", "#DDDDDD")
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = as.factor(x), y = y, fill = cor)) +
      ggplot2::geom_col(color = "black", alpha = 0.8) +
      ggplot2::scale_fill_identity() +
      ggplot2::labs(title = tit, subtitle = sprintf("Probabilidade: %.2f%%", res * 100),
                    x = "Número de Ocorrências (X)", y = "Probabilidade", caption = "estatR")
    
    if (exists("meu_tema")) p <- p + meu_tema(estilo = 1) else p <- p + ggplot2::theme_minimal()
    suppressMessages(print(p))
  }
  invisible(res)
}
