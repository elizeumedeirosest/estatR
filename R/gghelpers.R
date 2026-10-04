# ==============================================================================
# MÓDULO: GGHELPERS
# Funções auxiliares para gráficos ggplot2: tema e paletas do estatR
# ==============================================================================

# ============================================
# BANCO DE CORES E PALETAS
# ============================================

.paletas_estatR <- list(
  academic = c("#377EB8", "#E31A1C", "seagreen1", "darkorchid1", "#666666"),
  vibrant  = c("orangered1", "deepskyblue1", "gold", "magenta", "#33A02C"),
  contrast = c("firebrick1", "royalblue1", "limegreen", "darkorchid1", "cyan1"),
  nature   = c("seagreen1", "sienna1", "royalblue1", "gold", "darkorchid1"),
  neon     = c("dodgerblue1", "orangered1", "limegreen", "magenta", "gold"),
  electric = c("dodgerblue1", "firebrick1", "cyan1", "limegreen", "gold"),
  aurora   = c("lavenderblush2", "plum3", "mediumpurple3", "royalblue3", "midnightblue"),
  ember    = c("cornsilk2",      "burlywood3", "peru",      "sienna3",   "maroon4"),
  dusk     = c("steelblue3", "darksalmon", "mediumseagreen", "mediumpurple3", "cadetblue4"),
  terra    = c("peru", "cadetblue3", "darkseagreen3", "rosybrown3", "slateblue3")
)

.formatar_label_paleta <- function(x) {
  x <- gsub("([a-z])([0-9])", "\\1 \\2", x)
  x <- gsub("_|\\.", " ", x)
  # str_to_title sem depender do stringr
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}

.remover_escalas <- function(plot, aesthetics) {
  plot$scales$scales <- Filter(function(sc) {
    !any(aesthetics %in% sc$aesthetics)
  }, plot$scales$scales)
  plot
}

# Configura as ações ao somar a paleta a um ggplot
#' @export
ggplot_add.paleta_estatR_layer <- function(object, plot, object_name) {
  pal <- .paletas_estatR[[object$nome]]
  
  palette_fn <- function(n) {
    if (n <= length(pal)) unname(pal[1:n])
    else grDevices::colorRampPalette(pal)(n)
  }
  
  plot <- .remover_escalas(plot, c("fill", "colour", "color"))
  
  has_fill_map <- !is.null(plot$mapping[["fill"]]) || any(vapply(plot$layers, function(l) !is.null(l$mapping[["fill"]]), logical(1)))
  has_colour_map <- !is.null(plot$mapping[["colour"]]) || !is.null(plot$mapping[["color"]]) || any(vapply(plot$layers, function(l) !is.null(l$mapping[["colour"]]) || !is.null(l$mapping[["color"]]), logical(1)))
  
  if (has_fill_map) {
    plot <- plot + ggplot2::discrete_scale("fill", palette = palette_fn, labels = .formatar_label_paleta)
  }
  if (has_colour_map) {
    plot <- plot + ggplot2::discrete_scale("colour", palette = palette_fn, labels = .formatar_label_paleta, guide = "none")
  }
  
  # Recolore geometrias estáticas
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    geom_class <- class(layer$geom)
    
    has_fill_mapping <- !is.null(plot$mapping$fill) || !is.null(layer$mapping$fill)
    has_color_mapping <- !is.null(plot$mapping$color) || !is.null(plot$mapping$colour) || !is.null(layer$mapping$color) || !is.null(layer$mapping$colour)
    
    if ("GeomPoint" %in% geom_class) {
      if (!has_color_mapping && !is.null(layer$aes_params$colour)) layer$aes_params$colour <- pal[1]
      if (!has_fill_mapping && !is.null(layer$aes_params$fill)) layer$aes_params$fill <- pal[1]
    } else if ("GeomBar" %in% geom_class || "GeomCol" %in% geom_class || "GeomRect" %in% geom_class) {
      if (!has_fill_mapping && !is.null(layer$aes_params$fill)) layer$aes_params$fill <- pal[1]
    } else if ("GeomLine" %in% geom_class) {
      if (!has_color_mapping && !is.null(layer$aes_params$colour)) layer$aes_params$colour <- pal[1]
    } else if ("GeomPolygon" %in% geom_class || "GeomArea" %in% geom_class) {
      if (!has_fill_mapping && !is.null(layer$aes_params$fill)) layer$aes_params$fill <- pal[1]
    }
  }
  
  plot
}

# ============================================
# FUNÇÕES EXPORTADAS
# ============================================

#' @title Paleta de Cores do estatR
#' @description Aplica paletas de cores padronizadas aos gráficos ggplot2.
#' @param nome Nome da paleta (ex: "academic", "vibrant", "nature"). Se NULL, plota todas as paletas disponíveis.
#' @return Um objeto ggplot_add para adicionar ao gráfico, ou (se chamada sem argumentos) 
#' plota um demonstrativo visual de todas as paletas na aba Plots.
#' @export
paleta_estatR <- function(nome = NULL) {
  nomes_validos <- names(.paletas_estatR)
  
  # Se o usuário não passou nome, mostra o plot demonstrativo de TODAS as paletas
  if (is.null(nome)) {
    cores_df <- data.frame()
    for (i in seq_along(nomes_validos)) {
      pal_nome <- nomes_validos[i]
      pal_cores <- .paletas_estatR[[pal_nome]]
      # Cria um grid para plotar
      tmp <- data.frame(
        x = seq_along(pal_cores),
        y = length(nomes_validos) - i + 1,
        cor = pal_cores,
        paleta = sprintf("%d. %s", i, pal_nome)
      )
      cores_df <- rbind(cores_df, tmp)
    }
    
    p <- ggplot2::ggplot(cores_df, ggplot2::aes(x = x, y = as.factor(y), fill = I(cor))) +
      ggplot2::geom_tile(color = "white", linewidth = 1) +
      ggplot2::scale_y_discrete(labels = rev(unique(cores_df$paleta))) +
      ggplot2::theme_void() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(hjust = 1, margin = ggplot2::margin(r = 10), face = "bold", size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 20)),
        plot.margin = ggplot2::margin(20, 20, 20, 20)
      ) +
      ggplot2::labs(title = "Paletas de Cores do estatR")
    
    print(p)
    return(invisible(NULL))
  }
  
  if (!nome %in% nomes_validos) {
    stop("Paleta n\u00e3o encontrada. Op\u00e7\u00f5es: ", paste(nomes_validos, collapse = ", "))
  }
  
  # Se o usuário passou nome MAS usou sozinho (não somou a um ggplot)
  # A classe set_cores_layer será avaliada pelo R console e podemos interceptar com um método print
  obj <- structure(list(nome = nome), class = c("paleta_estatR_layer", "ggproto"))
  return(obj)
}

#' @export
print.paleta_estatR_layer <- function(x, ...) {
  # Se o objeto foi impresso sozinho (não somado), mostra as cores daquela paleta específica
  pal_cores <- .paletas_estatR[[x$nome]]
  df <- data.frame(x = seq_along(pal_cores), y = 1, cor = pal_cores)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = as.factor(x), y = y, fill = I(cor))) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = cor), color = ifelse(colSums(grDevices::col2rgb(pal_cores)) > 450, "black", "white")) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 20))
    ) +
    ggplot2::labs(title = sprintf("Paleta: %s", x$nome))
  
  print(p)
  invisible(x)
}


# ============================================
# TEMA
# ============================================

.detectar_eixo_grade <- function(plot) {
  obter_coluna <- function(aes_role) {
    for (l in plot$layers) {
      m <- l$mapping[[aes_role]]
      if (is.null(m)) m <- plot$mapping[[aes_role]]
      if (is.null(m)) next
      d <- if (inherits(l$data, "waiver") || is.null(l$data)) plot$data else l$data
      if (is.function(d)) d <- tryCatch(d(plot$data), error = function(e) NULL)
      nm <- tryCatch(rlang::as_label(m), error = function(e) NULL)
      if (!is.null(d) && !is.null(nm) && nm %in% names(d)) return(d[[nm]])
    }
    NULL
  }
  
  x_continua <- is.numeric(obter_coluna("x"))
  y_continua <- is.numeric(obter_coluna("y"))
  
  eixo_valor <- if (y_continua && !x_continua) "y" else if (x_continua && !y_continua) "x" else "y"
  tem_flip <- inherits(plot$coordinates, "CoordFlip")
  if (tem_flip) (if (eixo_valor == "x") "y" else "x") else eixo_valor
}

#' @title Tema Padronizado do estatR
#' @description Aplica formata\u00e7\u00e3o limpa e profissional aos gr\u00e1ficos do ggplot2.
#' @param fonte Nome da fonte do sistema (padr\u00e3o: "sans", carrega Arial/Helvetica nativamente).
#' @param estilo 1 (Classic limpo), 2 (Minimal com grades), 3 (Limpo absoluto).
#' @param escala Escala de tamanho dos textos.
#' @param inclinar Grau de inclina\u00e7\u00e3o dos r\u00f3tulos do eixo X (ex: 45). Se TRUE, assume 45.
#' @param modo "light" (padr\u00e3o) ou "dark" (fundo escuro).
#' @param grade Posi\u00e7\u00e3o da grade no estilo 2 ("auto", "x", "y", "dupla", "nenhuma").
#' @return Objeto de tema para adicionar ao ggplot.
#' @export
tema_estatR <- function(fonte = "sans", estilo = 1, escala = 1, inclinar = FALSE, modo = "light", grade = "auto") {
  if (is.numeric(fonte)) {
    estilo <- fonte
    fonte  <- "sans"
  }
  structure(list(fonte = fonte, estilo = estilo, escala = escala, inclinar = inclinar, modo = modo, grade = grade), class = "tema_estatR_layer")
}

#' @export
ggplot_add.tema_estatR_layer <- function(object, plot, object_name) {
  fonte <- object$fonte
  estilo <- object$estilo
  escala <- object$escala
  inclinar <- object$inclinar
  modo <- object$modo
  grade_opt <- object$grade
  
  is_dark <- modo == "dark"
  cor_fundo <- if (is_dark) "#1A1A1A" else "white"
  cor_texto <- if (is_dark) "#E5E5E5" else "black"
  cor_sub <- if (is_dark) "#A0A0A0" else "gray40"
  cor_linha <- if (is_dark) "#E5E5E5" else "black"
  
  angulo_inclinacao <- if (is.numeric(inclinar)) inclinar else if (isTRUE(inclinar)) 45 else NULL
  inclinacao_text <- if (!is.null(angulo_inclinacao)) {
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angulo_inclinacao, hjust = 1, vjust = 1))
  } else {
    ggplot2::theme()
  }
  
  textos_comuns <- ggplot2::theme(
    text              = ggplot2::element_text(family = fonte, size = 16 * escala, color = cor_texto),
    axis.title.x      = ggplot2::element_text(size = 15 * escala, color = cor_texto, margin = ggplot2::margin(t = 15)),
    axis.title.y      = ggplot2::element_text(size = 15 * escala, color = cor_texto, margin = ggplot2::margin(r = 15)),
    axis.text         = ggplot2::element_text(size = 14 * escala, color = cor_texto),
    plot.title        = ggplot2::element_text(face = "bold", size = 16 * escala, color = cor_texto),
    plot.subtitle     = ggplot2::element_text(color = cor_sub, size = 13 * escala),
    legend.title      = ggplot2::element_text(size = 14 * escala, color = cor_texto, margin = ggplot2::margin(b = 8)),
    legend.text       = ggplot2::element_text(size = 13 * escala, color = cor_texto),
    legend.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo),
    legend.key        = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo),
    strip.text        = ggplot2::element_text(size = 14 * escala, face = "bold", color = cor_texto),
    strip.background  = ggplot2::element_rect(fill = if(is_dark) "#333333" else "gray90", color = if(is_dark) "#333333" else "gray90"),
    plot.background   = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo, linewidth = 1)
  )
  
  tema_base <- if (estilo == 2) {
    cor_grade_estilo2 <- if (is_dark) "gray30" else "gray92"
    eixo_grade <- if (identical(grade_opt, "auto")) .detectar_eixo_grade(plot) else grade_opt
    
    grid_theme <- if (is.na(eixo_grade) || eixo_grade == "nenhuma") {
      ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    } else if (eixo_grade == "dupla") {
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5), panel.grid.major.x = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5))
    } else if (eixo_grade == "y") {
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5), panel.grid.major.x = ggplot2::element_blank())
    } else {
      ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5), panel.grid.major.y = ggplot2::element_blank())
    }
    
    ggplot2::theme_minimal(base_family = fonte) + textos_comuns + grid_theme +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + inclinacao_text
  } else if (estilo == 3) {
    ggplot2::theme_minimal(base_family = fonte) + textos_comuns +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + inclinacao_text
  } else {
    ggplot2::theme_classic(base_family = fonte) + textos_comuns +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo), axis.line = ggplot2::element_line(color = cor_linha, linewidth = 0.35), axis.ticks = ggplot2::element_line(color = cor_linha, linewidth = 0.35), axis.ticks.length = ggplot2::unit(4, "pt"), axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)), axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 7))) + inclinacao_text
  }
  
  if (is_dark) {
    for (i in seq_along(plot$layers)) {
      layer <- plot$layers[[i]]
      if (!is.null(layer$aes_params$color) && layer$aes_params$color == "white") layer$aes_params$color <- "#1A1A1A"
      if (!is.null(layer$aes_params$fill) && layer$aes_params$fill == "white") layer$aes_params$fill <- "#1A1A1A"
    }
  }
  
  plot + tema_base
}
