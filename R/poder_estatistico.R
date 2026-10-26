# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
# M\u00d3DULO: PODER ESTAT\u00cdSTICO
# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500

# ── Helpers internos ──────────────────────────────────────────────────────────

# Converte r\u00f3tulo de efeito em valor num\u00e9rico conforme o tipo de teste
.pe_resolver_efeito <- function(efeito, teste) {
  if (is.numeric(efeito)) return(efeito)
  efeito <- tolower(trimws(efeito))
  tab <- list(
    t_uma_amostra  = c(pequeno = 0.20, medio = 0.50, grande = 0.80),
    t_duas_amostras = c(pequeno = 0.20, medio = 0.50, grande = 0.80),
    t_pareado      = c(pequeno = 0.20, medio = 0.50, grande = 0.80),
    anova          = c(pequeno = 0.10, medio = 0.25, grande = 0.40),
    correlacao     = c(pequeno = 0.10, medio = 0.30, grande = 0.50),
    proporcao      = c(pequeno = 0.20, medio = 0.50, grande = 0.80)
  )
  vals <- tab[[teste]]
  efeito_norm <- switch(efeito,
    "pequeno" = "pequeno", "small"  = "pequeno",
    "medio"   = "medio",   "medium" = "medio",   "m\u00e9dio" = "medio",
    "grande"  = "grande",  "large"  = "grande",
    stop(sprintf("Efeito '%s' n\u00e3o reconhecido. Use num\u00e9rico ou 'pequeno', 'medio', 'grande'.", efeito), call. = FALSE)
  )
  unname(vals[efeito_norm])
}

.pe_label_efeito <- function(efeito_raw) {
  if (is.numeric(efeito_raw)) return(NULL)
  tolower(trimws(efeito_raw))
}

# Calcular poder dado n (retorna escalar entre 0 e 1)
.pe_calcular_poder <- function(teste, efeito_val, n, alpha, grupos, bilateral) {
  switch(teste,
    t_uma_amostra = {
      df  <- n - 1
      tc  <- qt(1 - alpha / ifelse(bilateral, 2, 1), df = df)
      lam <- efeito_val * sqrt(n)
      if (bilateral)
        1 - pt(tc, df = df, ncp = lam) + pt(-tc, df = df, ncp = lam)
      else
        1 - pt(tc, df = df, ncp = lam)
    },
    t_duas_amostras = {
      df  <- 2 * n - 2
      tc  <- qt(1 - alpha / ifelse(bilateral, 2, 1), df = df)
      lam <- efeito_val * sqrt(n / 2)
      if (bilateral)
        1 - pt(tc, df = df, ncp = lam) + pt(-tc, df = df, ncp = lam)
      else
        1 - pt(tc, df = df, ncp = lam)
    },
    t_pareado = {
      df  <- n - 1
      tc  <- qt(1 - alpha / ifelse(bilateral, 2, 1), df = df)
      lam <- efeito_val * sqrt(n)
      if (bilateral)
        1 - pt(tc, df = df, ncp = lam) + pt(-tc, df = df, ncp = lam)
      else
        1 - pt(tc, df = df, ncp = lam)
    },
    anova = {
      lam <- efeito_val^2 * n * grupos
      fc  <- qf(1 - alpha, df1 = grupos - 1, df2 = grupos * (n - 1))
      1 - pf(fc, df1 = grupos - 1, df2 = grupos * (n - 1), ncp = lam)
    },
    correlacao = {
      # Teste t para correlacao: t = r*sqrt(n-2)/sqrt(1-r^2)
      df  <- n - 2
      if (df <= 0) return(0)
      tc  <- qt(1 - alpha / ifelse(bilateral, 2, 1), df = df)
      r   <- efeito_val
      lam <- r * sqrt(n - 2) / sqrt(1 - r^2)
      if (bilateral)
        1 - pt(tc, df = df, ncp = lam) + pt(-tc, df = df, ncp = lam)
      else
        1 - pt(tc, df = df, ncp = lam)
    },
    proporcao = {
      # Usando aprox normal (z-teste)
      # d = efeito em d de Cohen para propor\u00e7\u00e3o
      zc  <- qnorm(1 - alpha / ifelse(bilateral, 2, 1))
      lam <- efeito_val * sqrt(n)
      if (bilateral)
        pnorm(lam - zc) + pnorm(-lam - zc)
      else
        pnorm(lam - zc)
    },
    stop("Teste n\u00e3o reconhecido.", call. = FALSE)
  )
}

# Busca bin\u00e1ria do n m\u00ednimo para atingir poder desejado
.pe_calcular_n <- function(teste, efeito_val, poder_alvo, alpha, grupos, bilateral) {
  for (n_try in 2:10000) {
    p_try <- .pe_calcular_poder(teste, efeito_val, n_try, alpha, grupos, bilateral)
    if (!is.na(p_try) && p_try >= poder_alvo) return(n_try)
  }
  return(NA_integer_)
}

# ─────────────────────────────────────────────────────────────────────────────

#' @title An\u00e1lise de Poder Estat\u00edstico
#' @description Calcula o poder de um teste estat\u00edstico dado o tamanho amostral,
#'   ou determina o n necess\u00e1rio para atingir um poder alvo. Suporta testes T,
#'   ANOVA, correla\u00e7\u00e3o e propor\u00e7\u00e3o sem depend\u00eancias externas.
#' @param teste Tipo de teste: \code{"t_uma_amostra"}, \code{"t_duas_amostras"},
#'   \code{"t_pareado"}, \code{"anova"}, \code{"correlacao"} ou
#'   \code{"proporcao"}.
#' @param efeito Tamanho de efeito num\u00e9rico ou texto (\code{"pequeno"},
#'   \code{"medio"}, \code{"grande"}).
#' @param n Tamanho amostral por grupo. Forne\u00e7a \code{n} para calcular o poder,
#'   ou \code{poder} para calcular o \code{n} necess\u00e1rio.
#' @param poder Poder desejado (0 a 1). Forne\u00e7a \code{poder} para calcular o
#'   \code{n} necess\u00e1rio.
#' @param alpha N\u00edvel de signific\u00e2ncia (padr\u00e3o \code{0.05}).
#' @param grupos N\u00famero de grupos para ANOVA (padr\u00e3o \code{2}).
#' @param hipotese \code{"bilateral"} (padr\u00e3o) ou \code{"unilateral"}.
#' @param grafico L\u00f3gico. Se \code{TRUE}, exibe curva de poder vs n.
#' @return Retorna invisivelmente uma lista com poder calculado, n necess\u00e1rio e
#'   par\u00e2metros utilizados.
#' @export
poder_estatistico <- function(
    teste    = c("t_uma_amostra", "t_duas_amostras", "t_pareado",
                 "anova", "correlacao", "proporcao"),
    efeito   = NULL,
    n        = NULL,
    poder    = NULL,
    alpha    = 0.05,
    grupos   = 2,
    hipotese = c("bilateral", "unilateral"),
    grafico  = TRUE
) {
  teste    <- match.arg(teste)
  hipotese <- match.arg(hipotese)

  if (is.null(efeito))
    stop("Forne\u00e7a o tamanho de efeito em 'efeito' (num\u00e9rico ou 'pequeno'/'medio'/'grande').", call. = FALSE)
  if (is.null(n) && is.null(poder))
    stop("Forne\u00e7a 'n' (para calcular poder) ou 'poder' (para calcular n necess\u00e1rio).", call. = FALSE)
  if (!is.null(n) && !is.null(poder))
    stop("Forne\u00e7a apenas um de 'n' ou 'poder', n\u00e3o ambos.", call. = FALSE)

  efeito_val  <- .pe_resolver_efeito(efeito, teste)
  efeito_lbl  <- .pe_label_efeito(efeito)
  bilateral   <- hipotese == "bilateral"

  # Nomes amig\u00e1veis
  teste_nome <- switch(teste,
    t_uma_amostra   = "Teste T para Uma Amostra",
    t_duas_amostras = "Teste T para Duas Amostras Independentes",
    t_pareado       = "Teste T Pareado",
    anova           = sprintf("ANOVA com %d grupos", grupos),
    correlacao      = "Teste de Correla\u00e7\u00e3o (Pearson)",
    proporcao       = "Teste de Propor\u00e7\u00e3o (Z)"
  )

  metrica_efeito <- switch(teste,
    anova        = "f",
    correlacao   = "r",
    "d"
  )

  efeito_str <- if (!is.null(efeito_lbl)) {
    cap <- paste0(toupper(substr(efeito_lbl, 1, 1)), substr(efeito_lbl, 2, nchar(efeito_lbl)))
    sprintf("%s = %.2f (%s)", metrica_efeito, efeito_val, cap)
  } else {
    sprintf("%s = %.2f", metrica_efeito, efeito_val)
  }

  # Modo de c\u00e1lculo
  if (!is.null(n)) {
    # Calcular PODER dado n
    poder_calc <- .pe_calcular_poder(teste, efeito_val, n, alpha, grupos, bilateral)
    n_calc     <- n
    n_total    <- n * ifelse(teste %in% c("t_duas_amostras", "anova"), grupos, 1)
  } else {
    # Calcular N dado poder desejado
    if (poder <= 0 || poder >= 1)
      stop("'poder' deve estar entre 0 e 1 (exclusive).", call. = FALSE)
    n_calc     <- .pe_calcular_n(teste, efeito_val, poder, alpha, grupos, bilateral)
    if (is.na(n_calc))
      stop("N\u00e3o foi poss\u00edvel determinar o n para o poder desejado (n > 10000).", call. = FALSE)
    poder_calc <- .pe_calcular_poder(teste, efeito_val, n_calc, alpha, grupos, bilateral)
    n_total    <- n_calc * ifelse(teste %in% c("t_duas_amostras", "anova"), grupos, 1)
  }

  # ── Imprimir resultados ──────────────────────────────────────────────────────
  .print_titulo("AN\u00c1LISE DE PODER ESTAT\u00cdSTICO")

  cat(sprintf("  Teste: %s\n", teste_nome))
  cat(sprintf(
    "  Hip\u00f3tese: %s   |   N\u00edvel \u03b1: %.0f%%   |   Tamanho de efeito: %s\n\n",
    if (bilateral) "Bilateral" else "Unilateral",
    alpha * 100,
    efeito_str
  ))

  .print_topico("RESULTADO PRINCIPAL")

  if (!is.null(n)) {
    cat(sprintf("  n fornecido (por grupo):          %d observa\u00e7\u00f5es\n", n_calc))
    cat(sprintf("  Poder calculado:                  %.1f%%\n",               poder_calc * 100))
    cat(sprintf("  n total do estudo:                %d observa\u00e7\u00f5es\n", n_total))
  } else {
    cat(sprintf("  Poder desejado:                   %.1f%%\n",               poder * 100))
    cat(sprintf("  n por grupo necess\u00e1rio:           %d observa\u00e7\u00f5es\n", n_calc))
    cat(sprintf("  Poder atingido:                   %.1f%%\n",               poder_calc * 100))
    cat(sprintf("  n total do estudo:                %d observa\u00e7\u00f5es\n", n_total))
  }
  cat("\n")

  .print_topico("INTERPRETA\u00c7\u00c3O")

  if (!is.null(n)) {
    adv <- if (poder_calc >= 0.80) "adequado" else "INSUFICIENTE (abaixo de 80%)"
    cat(sprintf(
      "  Com %d observa\u00e7\u00f5es (n = %d por grupo), o estudo tem %.1f%% de poder\n",
      n_total, n_calc, poder_calc * 100
    ))
    cat(sprintf(
      "  para detectar um efeito %s ao n\u00edvel de %.0f%% de signific\u00e2ncia %s.\n",
      efeito_str, alpha * 100, hipotese
    ))
    cat(sprintf("  O poder \u00e9 %s.\n\n", adv))
  } else {
    cat(sprintf(
      "  Com %d observa\u00e7\u00f5es por grupo (%d no total), o estudo tem %.1f%% de\n",
      n_calc, n_total, poder_calc * 100
    ))
    cat(sprintf(
      "  poder para detectar um efeito %s, ao n\u00edvel de\n",
      efeito_str
    ))
    cat(sprintf(
      "  %.0f%% de signific\u00e2ncia %s. Se o efeito realmente existir,\n",
      alpha * 100, hipotese
    ))
    cat(sprintf(
      "  o teste ter\u00e1 %.1f%% de chance de detect\u00e1-lo.\n\n",
      poder_calc * 100
    ))
  }

  .print_rodape()

  # ── Curva de poder ──────────────────────────────────────────────────────────
  if (grafico) {
    tryCatch({
      n_max_graf <- max(n_calc * 3, 50)
      n_seq      <- seq(2, n_max_graf, by = max(1, floor(n_max_graf / 200)))
      pod_seq    <- sapply(n_seq, function(ni)
        .pe_calcular_poder(teste, efeito_val, ni, alpha, grupos, bilateral))

      df_curva <- data.frame(n = n_seq, poder = pod_seq)

      tema_fn <- if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      p_graf <- ggplot2::ggplot(df_curva, ggplot2::aes(x = n, y = poder)) +
        ggplot2::geom_line(colour = "#377EB8", linewidth = 1) +
        ggplot2::geom_hline(
          yintercept = 0.80,
          linetype   = "dashed",
          colour     = "firebrick",
          linewidth  = 0.7
        ) +
        ggplot2::geom_vline(
          xintercept = n_calc,
          linetype   = "dotted",
          colour     = "darkgreen",
          linewidth  = 0.7
        ) +
        ggplot2::annotate(
          "text",
          x     = n_calc,
          y     = 0.05,
          label = sprintf("n = %d", n_calc),
          hjust = -0.15,
          colour = "darkgreen",
          size   = 4
        ) +
        ggplot2::scale_y_continuous(
          labels = scales::percent_format(accuracy = 1),
          limits = c(0, 1)
        ) +
        ggplot2::labs(
          title    = "Curva de Poder Estat\u00edstico",
          subtitle = sprintf("%s | %s | \u03b1 = %.0f%%", teste_nome, efeito_str, alpha * 100),
          x        = "n por grupo",
          y        = "Poder (1 - \u03b2)"
        ) +
        tema_fn

      print(p_graf)
    }, error = function(e) {
      # Fallback sem scales
      tryCatch({
        n_max_graf <- max(n_calc * 3, 50)
        n_seq      <- seq(2, n_max_graf, by = max(1, floor(n_max_graf / 200)))
        pod_seq    <- sapply(n_seq, function(ni)
          .pe_calcular_poder(teste, efeito_val, ni, alpha, grupos, bilateral))
        df_curva   <- data.frame(n = n_seq, poder = pod_seq)

        tema_fn <- if (exists("tema_estatR", mode = "function")) {
          tema_estatR(estilo = 2, grade = "dupla")
        } else {
          ggplot2::theme_minimal()
        }

        p_graf <- ggplot2::ggplot(df_curva, ggplot2::aes(x = n, y = poder)) +
          ggplot2::geom_line(colour = "#377EB8", linewidth = 1) +
          ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed", colour = "firebrick") +
          ggplot2::geom_vline(xintercept = n_calc, linetype = "dotted", colour = "darkgreen") +
          ggplot2::labs(
            title = "Curva de Poder Estat\u00edstico",
            x     = "n por grupo",
            y     = "Poder (1 - beta)"
          ) +
          tema_fn
        print(p_graf)
      }, error = function(e2) {
        message("  [aviso] N\u00e3o foi poss\u00edvel gerar o gr\u00e1fico de poder: ", e2$message)
      })
    })
  }

  invisible(list(
    teste       = teste,
    efeito      = efeito_val,
    n           = n_calc,
    n_total     = n_total,
    poder       = poder_calc,
    alpha       = alpha,
    grupos      = grupos,
    hipotese    = hipotese
  ))
}
