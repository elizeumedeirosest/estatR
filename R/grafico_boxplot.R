# ==============================================================================
# MÓDULO: GRÁFICO BOXPLOT
# ==============================================================================

#' @title Gráfico de Boxplot Avançado
#' @description Cria um boxplot completo com opções de violino, dispersão,
#' detecção de outliers, médias e personalização visual.
#'
#' @param dados Data frame.
#' @param x Variável categórica (eixo X). Passe sem aspas.
#' @param y Variável contínua/numérica (eixo Y). Passe sem aspas.
#' @param grupo Variável opcional para agrupamento lado a lado. Passe sem aspas.
#' @param bigode Lógico. Desenha barras horizontais nas pontas dos bigodes.
#'   Quando \code{violino = TRUE} e este parâmetro não for informado, o padrão
#'   muda automaticamente para \code{FALSE}.
#' @param ponto_media Lógico. Desenha um diamante preto na média de cada caixa.
#' @param violino Lógico. Plota um violino ao fundo das caixas.
#' @param dispersao_pts Lógico. Plota os pontos reais (jitter) por cima da caixa,
#'   excluindo outliers.
#' @param arredondar_borda_caixa Lógico. Ativa arredondamento suave dos cantos da caixa.
#' @param ligacao_media Lógico. Conecta as médias dos grupos com linha tracejada.
#' @param ligacao_mediana Lógico. Conecta as medianas dos grupos com linha pontilhada.
#' @param outlier Lógico. Exibe outliers. Se \code{FALSE}, oculta completamente.
#' @param nomes_outliers Lógico ou string. Se \code{TRUE} usa rownames. Se string,
#'   usa o nome de uma coluna do data frame para rotular os outliers.
#' @param fator_iqr Numérico. Multiplicador do IQR para detectar outliers (padrão 1.5).
#' @param tam_dispersao_pts Numérico. Multiplicador do tamanho dos pontos.
#' @param tam_texto_outliers Numérico. Multiplicador do tamanho dos rótulos dos outliers.
#' @param paleta Numérico (1–5) para paletas internas do boxplot, ou string
#'   (\code{"academic"}, \code{"vibrant"}, etc.) para usar \code{paleta_estatR()}.
#' @param cor Cor fixa (hex ou nome R) para todos os boxplots, ignora a paleta.
#'   Pode ser um vetor com uma cor por grupo (na ordem dos níveis do eixo X).
#' @param destaque Vetor de nomes de categorias a destacar (as demais ficam em cinza).
#' @param ordenar Lógico. Se \code{TRUE}, ordena o eixo X pelas medianas (decrescente).
#' @param ordem_grupo Vetor manual com a ordem dos níveis do grupo.
#' @param nomes_grupo Vetor manual para renomear os níveis do grupo.
#' @param ordem_eixo Vetor manual com a ordem dos níveis do eixo X.
#' @param nomes_eixo Vetor manual para renomear as categorias do eixo X.
#'
#' @return Um objeto \code{ggplot}.
#'
#' @import ggplot2
#' @importFrom rlang enquo as_name quo_is_null sym
#' @export
grafico_boxplot <- function(
    dados, x, y, grupo = NULL,
    bigode         = TRUE,
    ponto_media    = TRUE,
    violino        = FALSE,
    dispersao_pts  = FALSE,
    arredondar_borda_caixa = FALSE,
    ligacao_media  = FALSE,
    ligacao_mediana = FALSE,
    outlier        = TRUE,
    nomes_outliers = FALSE,
    fator_iqr      = 1.5,
    tam_dispersao_pts   = 1,
    tam_texto_outliers  = 1,
    paleta  = 1,
    cor     = NULL,
    destaque = NULL,
    ordenar      = FALSE,
    ordem_grupo  = NULL,
    nomes_grupo  = NULL,
    ordem_eixo   = NULL,
    nomes_eixo   = NULL
) {

  # Quando violino = TRUE e bigode não foi especificado, desativa o bigode
  if (missing(bigode) && isTRUE(violino)) bigode <- FALSE

  # ---- captura NSE ---------------------------------------------------------
  q_x     <- rlang::enquo(x)
  q_y     <- rlang::enquo(y)
  q_grupo <- rlang::enquo(grupo)

  nome_x    <- rlang::as_name(q_x)
  nome_y    <- rlang::as_name(q_y)
  tem_grupo <- !rlang::quo_is_null(q_grupo)
  nome_grupo <- if (tem_grupo) rlang::as_name(q_grupo) else NULL

  # ---- validações ----------------------------------------------------------
  if (!nome_x %in% names(dados))
    stop(sprintf("grafico_boxplot: variavel '%s' nao encontrada.", nome_x))
  if (!nome_y %in% names(dados))
    stop(sprintf("grafico_boxplot: variavel '%s' nao encontrada.", nome_y))
  if (tem_grupo && !nome_grupo %in% names(dados))
    stop(sprintf("grafico_boxplot: variavel de grupo '%s' nao encontrada.", nome_grupo))

  # ---- rótulo de outliers --------------------------------------------------
  colunas_uso <- c(nome_x, nome_y)
  if (tem_grupo) colunas_uso <- c(colunas_uso, nome_grupo)

  nome_rotulo <- NULL
  if (is.character(nomes_outliers)) {
    if (!nomes_outliers %in% names(dados))
      stop("grafico_boxplot: coluna para 'nomes_outliers' nao encontrada.")
    colunas_uso <- c(colunas_uso, nomes_outliers)
    nome_rotulo <- nomes_outliers
  } else if (isTRUE(nomes_outliers)) {
    dados$.rotulo_outlier <- rownames(dados)
    colunas_uso  <- c(colunas_uso, ".rotulo_outlier")
    nome_rotulo  <- ".rotulo_outlier"
  }

  # ---- limpeza e ordenação -------------------------------------------------
  dados <- dados[stats::complete.cases(dados[, colunas_uso, drop = FALSE]), ]

  if (isTRUE(ordenar)) {
    meds       <- tapply(dados[[nome_y]], dados[[nome_x]], stats::median)
    ordem_eixo <- names(sort(meds, decreasing = TRUE))
  }
  dados[[nome_x]] <- if (!is.null(ordem_eixo)) {
    factor(dados[[nome_x]], levels = ordem_eixo)
  } else {
    factor(dados[[nome_x]])
  }
  if (!is.null(nomes_eixo)) levels(dados[[nome_x]]) <- nomes_eixo

  if (tem_grupo) {
    dados[[nome_grupo]] <- if (!is.null(ordem_grupo)) {
      factor(dados[[nome_grupo]], levels = ordem_grupo)
    } else {
      factor(dados[[nome_grupo]])
    }
    if (!is.null(nomes_grupo)) levels(dados[[nome_grupo]]) <- nomes_grupo
  }

  # ---- paletas internas do boxplot (mono / básicas) -------------------------
  .paletas_box <- list(
    "1" = list(fill = "#CCCCCC", borda = "#555555", out = "#D90429"),
    "2" = list(fill = "#D9E2EC", borda = "#243B53", out = "#D4A017"),
    "3" = list(fill = "#D6F0FF", borda = "#0096C7", out = "#FF6B6B"),
    "4" = list(fill = "#E6ECF5", borda = "midnightblue", out = "gold"),
    "5" = list(fill = "#EAF3FB", borda = "#377EB8",      out = "#E31A1C")
  )
  .paletas_box_multi <- list(
    "1" = c("#555555", "#888888", "#AAAAAA", "#CCCCCC", "#EEEEEE"),
    "2" = c("#243B53", "#486581", "#6E9DC9", "#BCCCDC", "#D9E2EC"),
    "3" = c("#0096C7", "#48CAE4", "#90E0EF", "#ADE8F4", "#CEEDFD"),
    "4" = c("midnightblue", "#1B4F8A", "#2E75B6", "#5B9BD5", "#9DC3E6"),
    "5" = c("#377EB8", "#5E9ED6", "#82B6E0", "#AACEE8", "#C7DFEF")
  )

  usa_paleta_box_num <- is.numeric(paleta)
  idx_pal            <- as.character(if (usa_paleta_box_num) paleta else 1)
  pal_box            <- if (idx_pal %in% names(.paletas_box)) .paletas_box[[idx_pal]] else .paletas_box[["1"]]

  # ---- flag de visibilidade da legenda -------------------------------------
  hide_leg <- !tem_grupo

  # ---- cálculo antecipado dos outliers ------------------------------------
  grupos_iqr <- if (tem_grupo) c(nome_x, nome_grupo) else nome_x
  .marca_outliers <- function(d) {
    if (nrow(d) == 0) { d$.is_out <- logical(0); return(d) }
    yv <- d[[nome_y]]
    q1 <- stats::quantile(yv, 0.25, na.rm = TRUE)
    q3 <- stats::quantile(yv, 0.75, na.rm = TRUE)
    iq <- q3 - q1
    d$.is_out <- yv < (q1 - fator_iqr * iq) | yv > (q3 + fator_iqr * iq)
    d
  }
  split_key <- if (length(grupos_iqr) > 1) {
    interaction(dados[, grupos_iqr], drop = TRUE)
  } else {
    dados[[grupos_iqr]]
  }
  dados <- do.call(rbind, lapply(split(dados, split_key), .marca_outliers))

  # ---- fill_var (o que mapeia cor das caixas) ------------------------------
  fill_var <- if (tem_grupo) q_grupo else q_x

  # ===========================================================================
  # CONSTRUÇÃO DO MAPA DE CORES (named vector: nível → cor)
  # Isso centraliza toda a lógica de cor em um único lugar e elimina bugs.
  # ===========================================================================
  # A variável de preenchimento mapeada no aes
  nome_fill_var <- rlang::as_name(fill_var)
  nivs_fill     <- levels(dados[[nome_fill_var]])
  n_fill        <- length(nivs_fill)

  # Inicializa com a paleta padrão
  if (usa_paleta_box_num) {
    if (tem_grupo) {
      multi   <- if (idx_pal %in% names(.paletas_box_multi)) .paletas_box_multi[[idx_pal]] else .paletas_box_multi[["1"]]
      base_cores <- if (n_fill <= length(multi)) multi[seq_len(n_fill)] else grDevices::colorRampPalette(multi)(n_fill)
    } else {
      base_cores <- rep(pal_box$fill, n_fill)
    }
  } else {
    # paleta string: começa com cinza neutro, será sobrescrita depois via paleta_estatR()
    base_cores <- rep("#CCCCCC", n_fill)
  }
  mapa_cores <- stats::setNames(base_cores, nivs_fill)

  # Aplica `cor` (vetor ou escalar) sobrescrevendo as cores base
  usa_cor_manual <- FALSE
  if (!is.null(cor)) {
    usa_cor_manual <- TRUE
    cors_expandidas <- if (length(cor) == 1) {
      rep(cor, n_fill)
    } else {
      # cor é vetor: mapeia posicionalmente aos níveis do fill_var
      if (length(cor) < n_fill) rep(cor, length.out = n_fill) else cor[seq_len(n_fill)]
    }
    mapa_cores <- stats::setNames(cors_expandidas, nivs_fill)
  }

  # Aplica `destaque`: destacados ficam com sua cor, demais ficam cinza
  usa_destaque <- FALSE
  if (!is.null(destaque)) {
    usa_destaque <- TRUE
    destaque_chr <- as.character(destaque)

    # Cor para grupos NÃO destacados
    cinza_dest <- "gray85"

    # Cor para grupos DESTACADOS:
    # Prioridade: (1) cor[i] correspondente ao nível se cor for vetor, (2) paleta string[1], (3) pal_box$fill
    cor_para_destaque <- function(nivel) {
      if (usa_cor_manual && length(cor) > 1) {
        # Usa a cor correspondente ao índice desse nível no fill_var
        idx <- which(nivs_fill == nivel)
        if (length(idx) > 0 && idx <= length(cor)) return(cor[idx])
      }
      if (usa_cor_manual && length(cor) == 1) return(cor)
      if (is.character(paleta) && exists(".paletas_estatR", envir = parent.env(environment()), inherits = TRUE)) {
        pals <- tryCatch(get(".paletas_estatR", envir = parent.env(environment()), inherits = TRUE), error = function(e) NULL)
        if (!is.null(pals) && !is.null(pals[[paleta]])) return(pals[[paleta]][1])
      }
      pal_box$fill
    }

    for (nv in nivs_fill) {
      if (nv %in% destaque_chr) {
        mapa_cores[nv] <- cor_para_destaque(nv)
      } else {
        mapa_cores[nv] <- cinza_dest
      }
    }
  }

  # ===========================================================================
  # MONTAGEM DAS CAMADAS
  # ===========================================================================

  p <- ggplot2::ggplot(dados, ggplot2::aes(x = !!q_x, y = !!q_y, fill = !!fill_var))

  # 0. Violino (fundo, antes da caixa)
  if (isTRUE(violino)) {
    p <- p + ggplot2::geom_violin(
      alpha      = 0.22,
      color      = "black",
      linewidth  = 0.4,
      trim       = FALSE,
      show.legend = !hide_leg
    )
  }

  # 1. Cap de bigode (antes da caixa para ficar por baixo)
  if (isTRUE(bigode)) {
    p <- p + ggplot2::stat_boxplot(
      geom     = "errorbar",
      width    = 0.22,
      color    = "gray30",
      position = if (tem_grupo) ggplot2::position_dodge(0.75) else "identity"
    )
  }

  # 2. Caixa principal
  geom_bx_args <- list(
    width         = if (isTRUE(violino)) 0.28 else 0.60,
    alpha         = 1,
    outlier.shape = NA,
    show.legend   = !hide_leg
  )
  if (tem_grupo)                      geom_bx_args$position  <- ggplot2::position_dodge(0.75)
  if (isTRUE(arredondar_borda_caixa)) geom_bx_args$linejoin  <- "round"

  p <- p + do.call(ggplot2::geom_boxplot, geom_bx_args)

  # 3. Outliers — cor sempre via aes(colour = fill_var) + scale_color_manual
  #    para ficar 100% sincronizada com a cor da caixa correspondente.
  if (isTRUE(outlier)) {
    dados_out <- dados[dados$.is_out, ]
    if (nrow(dados_out) > 0) {
      pos_out <- if (tem_grupo) ggplot2::position_dodge(0.75) else "identity"

      p <- p + ggplot2::geom_point(
        data        = dados_out,
        mapping     = ggplot2::aes(
          x      = !!q_x,
          y      = !!q_y,
          group  = !!fill_var,
          colour = !!fill_var
        ),
        position    = pos_out,
        size        = 2.2 * tam_dispersao_pts,
        alpha       = 0.9,
        inherit.aes = FALSE,
        show.legend = FALSE
      )

      # Rótulos de outliers (ggrepel)
      if (!is.null(nome_rotulo) && requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data    = dados_out,
          mapping = ggplot2::aes(
            x     = !!q_x,
            y     = !!q_y,
            label = !!rlang::sym(nome_rotulo),
            group = !!fill_var
          ),
          position           = pos_out,
          size               = 3.5 * tam_texto_outliers,
          color              = "black",
          box.padding        = 0.5,
          point.padding      = 0.2,
          min.segment.length = 0,
          show.legend        = FALSE
        )
      }
    }
  }

  # 4. Dispersão por cima (excluindo outliers)
  if (isTRUE(dispersao_pts)) {
    dados_norm <- dados[!dados$.is_out, ]
    pos_jitter <- if (tem_grupo) {
      ggplot2::position_jitterdodge(jitter.width = 0.18, dodge.width = 0.75)
    } else {
      ggplot2::position_jitter(width = 0.18, seed = 42)
    }

    p <- p + ggplot2::geom_point(
      data        = dados_norm,
      mapping     = ggplot2::aes(x = !!q_x, y = !!q_y, group = !!fill_var),
      position    = pos_jitter,
      color       = "gray30",
      size        = tam_dispersao_pts,
      alpha       = 0.45,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  # 5. Ponto da média — sempre PRETO (param fixo, imune a qualquer scale)
  pos_sum <- if (tem_grupo) ggplot2::position_dodge(0.75) else "identity"

  if (isTRUE(ponto_media)) {
    p <- p + ggplot2::stat_summary(
      mapping      = ggplot2::aes(group = !!fill_var),
      fun          = mean,
      geom         = "point",
      shape        = 18,
      size         = 3.8,
      colour       = "black",
      fill         = "black",
      position     = pos_sum,
      show.legend  = FALSE
    )
  }

  # 6. Linhas de ligação (média / mediana)
  if (isTRUE(ligacao_media)) {
    if (tem_grupo) {
      p <- p + ggplot2::stat_summary(
        ggplot2::aes(group = !!q_grupo),
        fun = mean, geom = "line",
        linetype = "dashed", linewidth = 0.8, color = "gray30",
        position = pos_sum, show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::stat_summary(
        ggplot2::aes(group = 1),
        fun = mean, geom = "line",
        linetype = "dashed", linewidth = 0.8, color = "gray30",
        show.legend = FALSE
      )
    }
  }

  if (isTRUE(ligacao_mediana)) {
    if (tem_grupo) {
      p <- p + ggplot2::stat_summary(
        ggplot2::aes(group = !!q_grupo),
        fun = stats::median, geom = "line",
        linetype = "dotted", linewidth = 0.8, color = "black",
        position = pos_sum, show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::stat_summary(
        ggplot2::aes(group = 1),
        fun = stats::median, geom = "line",
        linetype = "dotted", linewidth = 0.8, color = "black",
        show.legend = FALSE
      )
    }
  }

  # ===========================================================================
  # ESCALAS DE COR — um único ponto de aplicação
  # ===========================================================================

  # Quando a paleta é numérica OU há cor/destaque manual: usa scale_*_manual
  # com o mapa_cores construído acima (fill e colour apontam para o mesmo vetor).
  # Quando paleta é string e não há cor/destaque: delega ao sistema paleta_estatR().

  usa_scale_manual <- usa_paleta_box_num || usa_cor_manual || usa_destaque

  if (usa_scale_manual) {
    # Cor do outlier para paleta numérica sem cor/destaque manual
    cor_outlier_scale <- if (!usa_cor_manual && !usa_destaque) {
      stats::setNames(rep(pal_box$out, n_fill), nivs_fill)
    } else {
      # Outlier herda a cor do fill (mapa_cores) — perfeita sincronização
      mapa_cores
    }

    p <- p +
      ggplot2::scale_fill_manual(
        values = mapa_cores,
        guide  = if (hide_leg) "none" else ggplot2::guide_legend()
      ) +
      ggplot2::scale_color_manual(
        values = cor_outlier_scale,
        guide  = "none"
      )

  } else if (is.character(paleta)) {
    # Paleta estatR (string): delega ao sistema de paletas
    if (exists("paleta_estatR", mode = "function")) p <- p + paleta_estatR(paleta)
    if (hide_leg) p <- p + ggplot2::guides(fill = "none", color = "none")
  }

  p
}
