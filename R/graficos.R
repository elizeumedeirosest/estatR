# ============================================
# MÓDULO INTERNO: GRÁFICOS (Adaptado do metaR)
# ============================================
# Este arquivo contém as funções de tema e visualização essenciais
# para tornar o estatR auto-suficiente, sem depender do metaR.

#' Tema padronizado para gráficos do estatR
#'
#' @param fonte Nome da família da fonte (ex: "Segoe UI", "Arial").
#' @param estilo Inteiro: 1 (classic limpo), 2 (minimal com grade no valor), 3 (branco puro).
#' @param escala Multiplicador de tamanho para fontes e espaçamentos.
#' @param inclinar Lógico ou numérico. Se TRUE, inclina o eixo X a 45 graus.
#' @param modo "light" (padrão) ou "dark".
#' @param grade "auto" (padrão), "x", "y", "dupla", "nenhuma".
#' @export
meu_tema <- function(fonte = "Segoe UI", estilo = 2, escala = 1, inclinar = FALSE,
                     modo = "light", grade = "auto") {

  if (is.numeric(fonte)) {
    estilo <- fonte
    fonte  <- "sans"
  }

  estrutura <- list(fonte = fonte, estilo = estilo, escala = escala,
                     inclinar = inclinar, modo = modo, grade = grade)
  class(estrutura) <- "meu_tema_layer"
  estrutura
}

#' @keywords internal
.detectar_eixo_grade <- function(plot) {
  obter_coluna <- function(aes_role) {
    for (l in plot$layers) {
      m <- l$mapping[[aes_role]] %||% plot$mapping[[aes_role]]
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

  eixo_valor <- if (y_continua && !x_continua) "y"
                else if (x_continua && !y_continua) "x"
                else return("y")

  tem_flip <- inherits(plot$coordinates, "CoordFlip")
  if (tem_flip) {
    if (eixo_valor == "x") "y" else "x"
  } else {
    eixo_valor
  }
}

#' @keywords internal
#' @export
ggplot_add.meu_tema_layer <- function(object, plot, object_name) {
  fonte     <- object$fonte
  estilo    <- object$estilo
  escala    <- object$escala
  inclinar  <- object$inclinar
  modo      <- object$modo
  grade_opt <- object$grade %||% "auto"

  is_dark <- modo == "dark"
  cor_fundo  <- if (is_dark) "#1A1A1A" else "white"
  cor_texto  <- if (is_dark) "#E5E5E5" else "black"
  cor_sub    <- if (is_dark) "#A0A0A0" else "gray40"
  cor_linha  <- if (is_dark) "#E5E5E5" else "black"

  angulo_inclinacao <- if (is.numeric(inclinar)) inclinar else if (isTRUE(inclinar)) 45 else NULL

  inclinacao_text <- if (!is.null(angulo_inclinacao)) {
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angulo_inclinacao, hjust = 1, vjust = 1))
  } else {
    ggplot2::theme()
  }

  leg_pos_salva <- attr(plot, ".legenda_posicao")
  pos_legenda_atual <- if (isTRUE(attr(plot, ".legenda_aplicada")) && !is.null(leg_pos_salva)) leg_pos_salva else "right"

  textos_comuns <- ggplot2::theme(
    text              = ggplot2::element_text(family = fonte, size = 16 * escala, color = cor_texto),
    axis.title.x      = ggplot2::element_text(size = 15 * escala, color = cor_texto, margin = ggplot2::margin(t = 15)),
    axis.title.y      = ggplot2::element_text(size = 15 * escala, color = cor_texto, margin = ggplot2::margin(r = 15)),
    axis.text         = ggplot2::element_text(size = 14 * escala, color = cor_texto),
    plot.title        = ggplot2::element_text(face = "bold", size = 16 * escala, color = cor_texto),
    plot.subtitle     = ggplot2::element_text(color = cor_sub, size = 13 * escala),
    legend.position   = pos_legenda_atual,
    legend.title      = ggplot2::element_text(size = 14 * escala, color = cor_texto, margin = ggplot2::margin(b = 8)),
    legend.text       = ggplot2::element_text(size = 13 * escala, color = cor_texto),
    legend.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo),
    legend.key        = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo),
    legend.spacing.y  = ggplot2::unit(2 * escala, "pt"),
    strip.text        = ggplot2::element_text(size = 14 * escala, face = "bold", color = cor_texto),
    strip.background  = ggplot2::element_rect(fill = if(is_dark) "#333333" else "gray90", color = if(is_dark) "#333333" else "gray90"),
    plot.background   = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo, linewidth = 1)
  )

  tema_base <- if (estilo == 2) {
    cor_grade_estilo2 <- if (is_dark) "gray30" else "gray92"
    eixo_grade <- if (identical(grade_opt, "auto")) .detectar_eixo_grade(plot) else if (grade_opt %in% c("x", "y", "dupla")) grade_opt else NA_character_
    
    grid_theme <- if (is.na(eixo_grade)) {
      ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    } else if (eixo_grade == "dupla") {
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5), panel.grid.major.x = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5))
    } else if (eixo_grade == "y") {
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5), panel.grid.major.x = ggplot2::element_blank())
    } else {
      ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = cor_grade_estilo2, linewidth = 0.5), panel.grid.major.y = ggplot2::element_blank())
    }

    ggplot2::theme_minimal(base_family = fonte) + textos_comuns + grid_theme +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo), panel.ontop = FALSE, axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + inclinacao_text
  } else if (estilo == 3) {
    ggplot2::theme_minimal(base_family = fonte) + textos_comuns +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + inclinacao_text
  } else {
    ggplot2::theme_classic(base_family = fonte) + textos_comuns +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = cor_fundo, color = cor_fundo), axis.line = ggplot2::element_line(color = cor_linha, linewidth = 0.35), axis.ticks = ggplot2::element_line(color = cor_linha, linewidth = 0.35), axis.ticks.length = ggplot2::unit(4, "pt"), axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)), axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 7))) + inclinacao_text
  }

  plot + tema_base
}

#' @keywords internal
.registrar_meu_tema_layer <- function() {
  registerS3method("ggplot_add", "meu_tema_layer", ggplot_add.meu_tema_layer, envir = asNamespace("ggplot2"))
}

# Auto-registro quando o pacote é carregado
if (isNamespaceLoaded("ggplot2")) {
  .registrar_meu_tema_layer()
} else {
  setHook(packageEvent("ggplot2", "onLoad"), function(...) .registrar_meu_tema_layer())
}

# Inicializador nativo do pacote para registrar o S3
.onLoad <- function(libname, pkgname) {
  .registrar_meu_tema_layer()
}
