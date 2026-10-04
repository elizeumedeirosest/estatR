# ==============================================================================
# MÓDULO: GRÁFICO BOXPLOT
# ==============================================================================

#' @title Gráfico de Boxplot Avançado
#' @description Cria um boxplot completo com opções de violino, dispersão,
#' detecção de outliers, médias e personalização visual.
#' 
#' @param dados Data frame.
#' @param x Variável categórica (eixo X). (Passe sem aspas)
#' @param y Variável contínua/numérica (eixo Y). (Passe sem aspas)
#' @param grupo Variável opcional para agrupamento lado a lado. (Passe sem aspas)
#' @param bigode Lógico. Desenha barras horizontais nas pontas dos bigodes.
#' @param ponto_media Lógico. Desenha um diamante na média.
#' @param violino Lógico. Plota um gráfico de violino ao fundo.
#' @param dispersao_pts Lógico. Plota os pontos reais (jitter) ao fundo.
#' @param arredondar_borda_caixa Lógico. Suaviza as bordas do boxplot (disponível visualmente via tema).
#' @param ligacao_media Lógico. Conecta as médias com linha tracejada.
#' @param ligacao_mediana Lógico. Conecta as medianas com linha tracejada.
#' @param outlier Lógico. Exibe os outliers. Se FALSE, oculta.
#' @param nomes_outliers Lógico ou string. Se TRUE usa rownames. Se string, usa o nome da coluna para rotular outliers.
#' @param fator_iqr Numérico. Multiplicador do IQR para limite de outliers (padrão 1.5).
#' @param tam_dispersao_pts Numérico. Tamanho dos pontos (outliers e dispersão).
#' @param tam_texto_outliers Numérico. Tamanho da fonte dos rótulos dos outliers.
#' @param paleta Numérico (1, 2, 3...) para paletas básicas, ou String ("academic", "vibrant") para paleta_estatR.
#' @param cor Cor fixa para pintar todos os boxplots ignorando a paleta.
#' @param destaque Vetor com categorias a serem destacadas (o resto fica em cinza).
#' @param ordenar Lógico. Ordena o eixo X pelas medianas (decrescente).
#' @param ordem_grupo Vetor manual definindo a ordem dos grupos.
#' @param nomes_grupo Vetor manual para renomear os grupos.
#' @param ordem_eixo Vetor manual definindo a ordem do eixo X.
#' @param nomes_eixo Vetor manual para renomear as categorias do eixo X.
#' 
#' @import ggplot2
#' @importFrom rlang enquo as_name quo_is_null sym !! :=
#' @export
grafico_boxplot <- function(
    dados, x, y, grupo = NULL,
    bigode = TRUE, ponto_media = TRUE, violino = FALSE, 
    dispersao_pts = FALSE, arredondar_borda_caixa = FALSE,
    ligacao_media = FALSE, ligacao_mediana = FALSE,
    outlier = TRUE, nomes_outliers = FALSE, fator_iqr = 1.5, 
    tam_dispersao_pts = 1, tam_texto_outliers = 1,
    paleta = 1, cor = NULL, destaque = NULL,
    ordenar = FALSE, ordem_grupo = NULL, nomes_grupo = NULL, 
    ordem_eixo = NULL, nomes_eixo = NULL
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Pacote ggplot2 necess\u00e1rio.")
  
  # Captura as variáveis sem aspas (NSE)
  q_x <- rlang::enquo(x)
  q_y <- rlang::enquo(y)
  q_grupo <- rlang::enquo(grupo)
  
  nome_x <- rlang::as_name(q_x)
  nome_y <- rlang::as_name(q_y)
  tem_grupo <- !rlang::quo_is_null(q_grupo)
  nome_grupo <- if (tem_grupo) rlang::as_name(q_grupo) else NULL
  
  if (!nome_x %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", nome_x))
  if (!nome_y %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", nome_y))
  if (tem_grupo && !nome_grupo %in% names(dados)) stop(sprintf("Vari\u00e1vel '%s' n\u00e3o encontrada.", nome_grupo))
  
  colunas_uso <- c(nome_x, nome_y)
  if (tem_grupo) colunas_uso <- c(colunas_uso, nome_grupo)
  
  nome_rotulo <- NULL
  if (is.character(nomes_outliers)) {
    if (!nomes_outliers %in% names(dados)) stop("Coluna para 'nomes_outliers' n\u00e3o encontrada.")
    colunas_uso <- c(colunas_uso, nomes_outliers)
    nome_rotulo <- nomes_outliers
  } else if (isTRUE(nomes_outliers)) {
    dados$.rotulo_outlier <- rownames(dados)
    colunas_uso <- c(colunas_uso, ".rotulo_outlier")
    nome_rotulo <- ".rotulo_outlier"
  }
  
  # Filtra NAs
  dados <- dados[stats::complete.cases(dados[, colunas_uso, drop = FALSE]), ]
  
  # Ordenação do Eixo X
  if (ordenar) {
    meds <- tapply(dados[[nome_y]], dados[[nome_x]], stats::median)
    ordem_eixo <- names(sort(meds, decreasing = TRUE))
  }
  if (!is.null(ordem_eixo)) {
    dados[[nome_x]] <- factor(dados[[nome_x]], levels = ordem_eixo)
  } else {
    dados[[nome_x]] <- factor(dados[[nome_x]])
  }
  if (!is.null(nomes_eixo)) levels(dados[[nome_x]]) <- nomes_eixo
  
  # Ordenação do Grupo
  if (tem_grupo) {
    if (!is.null(ordem_grupo)) dados[[nome_grupo]] <- factor(dados[[nome_grupo]], levels = ordem_grupo)
    else dados[[nome_grupo]] <- factor(dados[[nome_grupo]])
    if (!is.null(nomes_grupo)) levels(dados[[nome_grupo]]) <- nomes_grupo
  }
  
  fill_var <- if (tem_grupo) q_grupo else q_x
  
  # Cores Base
  cor_unica <- NULL
  cores_vetor <- NULL
  
  if (tem_grupo) {
    if (is.character(paleta)) {
      # Será resolvido pelo gghelpers
    } else {
      # Com grupo e paleta numerica: default é academic
      paleta <- "academic"
    }
  } else {
    if (!is.null(cor)) {
      cor_unica <- cor
    } else if (is.character(paleta)) {
      if (exists(".paletas_estatR")) cor_unica <- .paletas_estatR[[paleta]][1] else cor_unica <- "steelblue"
    } else {
      # Sem grupo, paleta numerica -> usa paleta monocromatica do metaR
      paletas_mono <- list(
        "1" = "#CCCCCC", "2" = "#D9E2EC", "3" = "#D6F0FF", 
        "4" = "#E6ECF5", "5" = "#EAF3FB"
      )
      idx <- as.character(paleta)
      cor_unica <- if (idx %in% names(paletas_mono)) paletas_mono[[idx]] else paletas_mono[["1"]]
    }
  }
  
  # ---------------------------------------------------------
  # CONSTRUÇÃO DO PLOT
  # ---------------------------------------------------------
  p <- ggplot2::ggplot(dados, ggplot2::aes(x = !!q_x, y = !!q_y, fill = !!fill_var))
  
  # Ocultar legenda se nao tiver grupo
  hide_leg <- !tem_grupo
  
  if (violino) {
    p <- p + ggplot2::geom_violin(alpha = 0.5, color = NA, trim = FALSE, show.legend = !hide_leg)
  }
  
  if (dispersao_pts) {
    p <- p + ggplot2::geom_jitter(
      ggplot2::aes(color = !!fill_var), 
      width = 0.2, size = tam_dispersao_pts, alpha = 0.4, show.legend = FALSE
    )
  }
  
  if (bigode) {
    p <- p + ggplot2::stat_boxplot(
      geom = "errorbar", width = 0.2, color = "gray20",
      position = if (tem_grupo) ggplot2::position_dodge(0.75) else "identity"
    )
  }
  
  mostrar_outlier_padrao <- if (outlier && !is.character(nome_rotulo) && !dispersao_pts) TRUE else FALSE
  
  geom_bx_args <- list(
    width = if (violino) 0.3 else 0.6, 
    alpha = 1,  # Alpha 1 para nao deixar o bigode vazar por tras da caixa
    outlier.size = tam_dispersao_pts,
    show.legend = !hide_leg
  )
  
  if (tem_grupo) geom_bx_args$position <- ggplot2::position_dodge(0.75)
  if (!mostrar_outlier_padrao) geom_bx_args$outlier.shape <- NA
  if (!is.null(cor_unica)) {
    geom_bx_args$fill <- cor_unica
    geom_bx_args$color <- "black"
  }
  if (arredondar_borda_caixa) {
    geom_bx_args$linejoin <- "round"
  }
  
  p <- p + do.call(ggplot2::geom_boxplot, geom_bx_args)
  
  # ---------------------------------------------------------
  # Outliers Customizados (com rótulos)
  # ---------------------------------------------------------
  if (outlier && is.character(nome_rotulo) && requireNamespace("ggrepel", quietly = TRUE)) {
    grupos_calc <- if (tem_grupo) c(nome_x, nome_grupo) else nome_x
    calc_out <- function(d) {
      q1 <- stats::quantile(d[[nome_y]], 0.25, na.rm = TRUE)
      q3 <- stats::quantile(d[[nome_y]], 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      d$is_outlier <- d[[nome_y]] < (q1 - fator_iqr * iqr) | d[[nome_y]] > (q3 + fator_iqr * iqr)
      d
    }
    dados_out <- do.call(rbind, lapply(split(dados, dados[, grupos_calc, drop=FALSE]), calc_out))
    dados_out <- dados_out[dados_out$is_outlier, ]
    
    if (nrow(dados_out) > 0) {
      pos_out <- if (tem_grupo) ggplot2::position_dodge(0.75) else "identity"
      if (!dispersao_pts) {
        p <- p + ggplot2::geom_point(
          data = dados_out, ggplot2::aes(x = !!q_x, y = !!q_y, group = !!fill_var),
          position = pos_out, size = tam_dispersao_pts, color = "black", alpha = 0.8,
          show.legend = FALSE
        )
      }
      p <- p + ggrepel::geom_text_repel(
        data = dados_out,
        ggplot2::aes(x = !!q_x, y = !!q_y, label = !!rlang::sym(nome_rotulo), group = !!fill_var),
        position = pos_out, size = 3.5 * tam_texto_outliers, color = "black",
        box.padding = 0.5, point.padding = 0.2, min.segment.length = 0, show.legend = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------
  # Médias e Linhas
  # ---------------------------------------------------------
  pos_sum <- if (tem_grupo) ggplot2::position_dodge(0.75) else "identity"
  if (ponto_media) {
    p <- p + ggplot2::stat_summary(fun = "mean", geom = "point", shape = 18, size = 3.5, color = "darkred", position = pos_sum, show.legend = FALSE)
  }
  if (ligacao_media) {
    if (tem_grupo) {
      p <- p + ggplot2::stat_summary(fun = "mean", geom = "line", ggplot2::aes(group = !!q_grupo, color = !!q_grupo), linetype = "dashed", linewidth = 0.8, position = pos_sum, show.legend = FALSE)
    } else {
      p <- p + ggplot2::stat_summary(fun = "mean", geom = "line", ggplot2::aes(group = 1), linetype = "dashed", linewidth = 0.8, color = "gray40", show.legend = FALSE)
    }
  }
  if (ligacao_mediana) {
    if (tem_grupo) {
      p <- p + ggplot2::stat_summary(fun = "median", geom = "line", ggplot2::aes(group = !!q_grupo, color = !!q_grupo), linetype = "dotted", linewidth = 0.8, position = pos_sum, show.legend = FALSE)
    } else {
      p <- p + ggplot2::stat_summary(fun = "median", geom = "line", ggplot2::aes(group = 1), linetype = "dotted", linewidth = 0.8, color = "black", show.legend = FALSE)
    }
  }
  
  # ---------------------------------------------------------
  # Cores e Paleta Final
  # ---------------------------------------------------------
  if (!is.null(destaque)) {
    nivs <- levels(dados[[rlang::as_name(fill_var)]])
    cores_destaque <- rep("gray85", length(nivs))
    cor_dest <- if (is.character(paleta)) {
      if (exists(".paletas_estatR")) .paletas_estatR[[paleta]][1] else "#E31A1C"
    } else { "#E31A1C" }
    
    cores_destaque[nivs %in% destaque] <- cor_dest
    p <- p + ggplot2::scale_fill_manual(values = cores_destaque) + ggplot2::scale_color_manual(values = cores_destaque)
    
  } else if (tem_grupo) {
    # Com grupo -> aplica paleta_estatR() (academic se não especificada)
    if (exists("paleta_estatR", mode = "function")) {
       p <- p + paleta_estatR(paleta)
    }
  }
  
  p
}
