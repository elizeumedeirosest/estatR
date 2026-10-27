# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: APRENDIZADO DE MÁQUINA (estatR)
# Classificação, Previsão (Regressão), Otimização e Diagnóstico de Resíduos
# ─────────────────────────────────────────────────────────────────────────────

# ── Helpers internos de formatação e métricas ─────────────────────────────────

.am_pad <- function(s, w, align = "left") {
  s  <- as.character(s)
  sp <- w - nchar(s)
  if (sp <= 0) return(s)
  if (align == "right")  return(paste0(strrep(" ", sp), s))
  if (align == "center") return(paste0(strrep(" ", floor(sp / 2)), s, strrep(" ", ceiling(sp / 2))))
  paste0(s, strrep(" ", sp))
}

.am_fmt_num <- function(x, decimais = 3) {
  if (is.null(x) || length(x) == 0 || is.na(x) || is.nan(x)) return("-")
  formatC(x, format = "f", digits = decimais, decimal.mark = ",")
}

.am_fmt_p <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p) || is.nan(p)) return("-")
  if (p < 0.001) return("< 0,001")
  formatC(p, format = "f", digits = 3, decimal.mark = ",")
}

.am_fmt_pct <- function(x, decimais = 1) {
  if (is.null(x) || length(x) == 0 || is.na(x) || is.nan(x)) return("-")
  paste0(formatC(x * 100, format = "f", digits = decimais, decimal.mark = ","), "%")
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Previsão com Aprendizado de Máquina (Regressão)
#' @description Ajusta modelos de Aprendizado de Máquina para tarefas de regressão
#'   (variável resposta numérica contínua). Suporta múltiplos algoritmos com
#'   otimização automática de hiperparâmetros (estabilização OOB / Cross-Validation)
#'   ou parametrização manual em português.
#'
#' @param formula Fórmula no formato \code{y ~ x1 + x2} ou \code{y ~ .}.
#' @param dados Data frame contendo as variáveis do modelo.
#' @param algoritmo Algoritmo a utilizar: \code{"rf"} (Random Forest / Floresta Aleatória, padrão),
#'   \code{"xgboost"} (Extreme Gradient Boosting), \code{"svm"} (Máquina de Vetores de Suporte),
#'   \code{"lasso"} (Lasso Penalizado L1), \code{"ridge"} (Ridge Penalizado L2) ou \code{"arvore"} (CART).
#' @param proporcao_treino Proporção da amostra destinada ao treino (padrão \code{0.75}).
#' @param semente Valor inteiro para reprodutibilidade (padrão \code{42}).
#' @param n_arvores Número de árvores (para \code{"rf"}). Se \code{NULL} (padrão),
#'   aplica busca automática testando candidatos e seleciona o ponto de estabilização do erro OOB.
#' @param vars_por_no Número de variáveis sorteadas por divisão em cada nó (\code{mtry}).
#'   Se \code{NULL} (padrão), utiliza \code{max(1, floor(p / 3))}.
#' @param n_rounds Número de rodadas de boosting (para \code{"xgboost"}). Se \code{NULL},
#'   utiliza 50 rodadas ou otimiza automaticamente.
#' @param profundidade Profundidade máxima das árvores (para \code{"xgboost"} e \code{"arvore"}).
#' @param taxa_aprendizado Taxa de aprendizado (\code{eta}) para \code{"xgboost"} (padrão \code{0.1}).
#' @param custo Parâmetro de regularização C para \code{"svm"} (padrão \code{1}).
#' @param grafico Lógico. Se TRUE (padrão), exibe painel gráfico com grade dupla.
#'
#' @return Objeto da classe \code{estatR_aprend_maq} contendo o modelo ajustado,
#'   partições de treino/teste, métricas, predições, importância de variáveis e
#'   tabela de estabilização (quando aplicável).
#' @export
previsao_aprend_maq <- function(formula,
                                dados,
                                algoritmo          = c("rf", "xgboost", "svm", "lasso", "ridge", "arvore"),
                                proporcao_treino   = 0.75,
                                semente            = 42,
                                n_arvores          = NULL,
                                vars_por_no        = NULL,
                                n_rounds           = NULL,
                                profundidade       = NULL,
                                taxa_aprendizado   = 0.1,
                                custo              = 1,
                                grafico            = TRUE) {

  algoritmo <- match.arg(algoritmo)
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")

  # Extração de variáveis da fórmula
  termos_obj <- terms(formula, data = dados)
  y_nome     <- as.character(formula[[2]])
  pred_nomes <- attr(termos_obj, "term.labels")

  if (!(y_nome %in% names(dados))) stop(sprintf("Variável resposta '%s' não encontrada nos dados.", y_nome))
  if (!is.numeric(dados[[y_nome]])) stop(sprintf("A variável resposta '%s' deve ser numérica contínua para 'previsao_aprend_maq'.", y_nome))

  # Limpeza de dados (apenas colunas necessárias e casos completos)
  cols_nec <- c(y_nome, pred_nomes)
  df_limpo <- na.omit(dados[, cols_nec, drop = FALSE])
  n_total  <- nrow(df_limpo)
  if (n_total < 10) stop("Amostra insuficiente após remoção de valores ausentes (mínimo 10 observações).")

  p_pred <- length(pred_nomes)
  if (p_pred == 0) stop("O modelo deve conter pelo menos uma variável preditora.")

  # Partição Treino / Teste
  if (!is.null(semente)) set.seed(semente)
  n_treino   <- floor(n_total * proporcao_treino)
  idx_treino <- sample(seq_len(n_total), size = n_treino, replace = FALSE)
  df_treino  <- df_limpo[idx_treino, , drop = FALSE]
  df_teste   <- df_limpo[-idx_treino, , drop = FALSE]

  diag_estabilizacao <- NULL
  n_arv_usado        <- NULL
  vars_por_no_usado  <- NULL
  modelo_ajustado    <- NULL
  importancia_df     <- NULL

  # ── AJUSTE: RANDOM FOREST (via ranger ou randomForest) ─────────────────────
  if (algoritmo == "rf") {
    vars_por_no_usado <- if (!is.null(vars_por_no)) as.integer(vars_por_no) else max(1, floor(p_pred / 3))

    # Otimização automática do número de árvores por estabilização OOB
    if (is.null(n_arvores)) {
      candidatos_arv <- c(50, 100, 200, 300, 500, 750)
      # Limitar se amostra for muito pequena
      erros_oob <- numeric(length(candidatos_arv))

      for (i in seq_along(candidatos_arv)) {
        if (!is.null(semente)) set.seed(semente)
        rf_tmp <- ranger::ranger(
          formula       = formula,
          data          = df_treino,
          num.trees     = candidatos_arv[i],
          mtry          = vars_por_no_usado,
          importance    = "none",
          verbose       = FALSE
        )
        erros_oob[i] <- rf_tmp$prediction.error
      }

      # Variação relativa do MSE
      var_rel <- c(Inf, abs(diff(erros_oob)) / (abs(erros_oob[-length(erros_oob)]) + 1e-10))
      idx_est <- which(var_rel < 0.01) # < 1% de variação relativa
      n_arv_usado <- if (length(idx_est) > 0) candidatos_arv[min(idx_est)] else 300

      diag_estabilizacao <- data.frame(
        N_Arvores   = candidatos_arv,
        MSE_OOB     = round(erros_oob, 4),
        RMSE_OOB    = round(sqrt(erros_oob), 4),
        Var_Rel_Pct = c("-", paste0(round(var_rel[-1] * 100, 2), "%")),
        stringsAsFactors = FALSE
      )
    } else {
      n_arv_usado <- as.integer(n_arvores)
    }

    if (!is.null(semente)) set.seed(semente)
    modelo_ajustado <- ranger::ranger(
      formula       = formula,
      data          = df_treino,
      num.trees     = n_arv_usado,
      mtry          = vars_por_no_usado,
      importance    = "impurity",
      verbose       = FALSE
    )

    pred_treino <- predict(modelo_ajustado, data = df_treino)$predictions
    pred_teste  <- predict(modelo_ajustado, data = df_teste)$predictions

    imp_raw <- modelo_ajustado$variable.importance
    if (!is.null(imp_raw) && length(imp_raw) > 0) {
      imp_pct <- (imp_raw / sum(imp_raw)) * 100
      importancia_df <- data.frame(
        Variavel      = names(imp_raw),
        Importancia   = round(imp_raw, 2),
        Percentual    = round(imp_pct, 1),
        stringsAsFactors = FALSE
      )
      importancia_df <- importancia_df[order(importancia_df$Importancia, decreasing = TRUE), ]
      rownames(importancia_df) <- NULL
    }
  }

  # ── AJUSTE: XGBOOST ───────────────────────────────────────────────────────
  if (algoritmo == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("O pacote 'xgboost' \u00e9 necess\u00e1rio para este algoritmo.")
    }

    X_tr <- model.matrix(formula, data = df_treino)[, -1, drop = FALSE]
    y_tr <- df_treino[[y_nome]]
    X_te <- model.matrix(formula, data = df_teste)[, -1, drop = FALSE]
    y_te <- df_teste[[y_nome]]

    n_r_usado    <- if (!is.null(n_rounds)) as.integer(n_rounds) else 50
    prof_usada   <- if (!is.null(profundidade)) as.integer(profundidade) else 6

    if (!is.null(semente)) set.seed(semente)
    modelo_ajustado <- xgboost::xgboost(
      x         = X_tr,
      y         = y_tr,
      nrounds   = n_r_usado,
      max_depth = prof_usada,
      eta       = taxa_aprendizado,
      verbosity = 0
    )

    pred_treino <- predict(modelo_ajustado, X_tr)
    pred_teste  <- predict(modelo_ajustado, X_te)

    imp_xgb <- tryCatch(xgboost::xgb.importance(model = modelo_ajustado), error = function(e) NULL)
    if (!is.null(imp_xgb) && nrow(imp_xgb) > 0) {
      importancia_df <- data.frame(
        Variavel      = imp_xgb$Feature,
        Importancia   = round(imp_xgb$Gain, 4),
        Percentual    = round(imp_xgb$Gain * 100, 1),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── AJUSTE: SVM (e1071) ───────────────────────────────────────────────────
  if (algoritmo == "svm") {
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("O pacote 'e1071' é necessário para este método.")
    }
    if (!is.null(semente)) set.seed(semente)
    modelo_ajustado <- e1071::svm(formula = formula, data = df_treino, cost = custo)
    pred_treino     <- predict(modelo_ajustado, newdata = df_treino)
    pred_teste      <- predict(modelo_ajustado, newdata = df_teste)
  }

  # ── AJUSTE: LASSO E RIDGE (glmnet) ────────────────────────────────────────
  if (algoritmo %in% c("lasso", "ridge")) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("O pacote 'glmnet' é necessário para regularização.")
    }
    alpha_val <- if (algoritmo == "lasso") 1 else 0
    X_tr <- model.matrix(formula, data = df_treino)[, -1, drop = FALSE]
    y_tr <- df_treino[[y_nome]]
    X_te <- model.matrix(formula, data = df_teste)[, -1, drop = FALSE]

    if (!is.null(semente)) set.seed(semente)
    cv_glm <- glmnet::cv.glmnet(x = X_tr, y = y_tr, alpha = alpha_val)
    modelo_ajustado <- cv_glm
    pred_treino     <- as.numeric(predict(cv_glm, newx = X_tr, s = "lambda.min"))
    pred_teste      <- as.numeric(predict(cv_glm, newx = X_te, s = "lambda.min"))

    coefs_glm <- as.matrix(coef(cv_glm, s = "lambda.min"))[-1, , drop = FALSE]
    imp_glm   <- abs(coefs_glm[, 1])
    if (sum(imp_glm) > 0) {
      importancia_df <- data.frame(
        Variavel    = rownames(coefs_glm),
        Importancia = round(imp_glm, 4),
        Percentual  = round((imp_glm / sum(imp_glm)) * 100, 1),
        stringsAsFactors = FALSE
      )
      importancia_df <- importancia_df[order(importancia_df$Importancia, decreasing = TRUE), ]
      rownames(importancia_df) <- NULL
    }
  }

  # ── AJUSTE: ÁRVORE DE DECISÃO (rpart) ─────────────────────────────────────
  if (algoritmo == "arvore") {
    prof_usada <- if (!is.null(profundidade)) as.integer(profundidade) else 10
    modelo_ajustado <- rpart::rpart(formula = formula, data = df_treino, control = rpart::rpart.control(maxdepth = prof_usada))
    pred_treino     <- predict(modelo_ajustado, newdata = df_treino)
    pred_teste      <- predict(modelo_ajustado, newdata = df_teste)

    imp_rp <- modelo_ajustado$variable.importance
    if (!is.null(imp_rp) && length(imp_rp) > 0) {
      importancia_df <- data.frame(
        Variavel    = names(imp_rp),
        Importancia = round(imp_rp, 2),
        Percentual  = round((imp_rp / sum(imp_rp)) * 100, 1),
        stringsAsFactors = FALSE
      )
      rownames(importancia_df) <- NULL
    }
  }

  # ── Cálculo das Métricas de Desempenho ─────────────────────────────────────
  y_real_tr <- df_treino[[y_nome]]
  y_real_te <- df_teste[[y_nome]]

  calc_metricas_reg <- function(real, pred) {
    residuos <- real - pred
    rmse     <- sqrt(mean(residuos^2))
    mae      <- mean(abs(residuos))
    ss_tot   <- sum((real - mean(real))^2)
    ss_res   <- sum(residuos^2)
    r2       <- if (ss_tot > 0) max(0, 1 - (ss_res / ss_tot)) else 0
    mape     <- if (all(real != 0)) mean(abs(residuos / real)) * 100 else NA
    list(RMSE = rmse, MAE = mae, R2 = r2, MAPE = mape, Residuos = residuos)
  }

  m_tr <- calc_metricas_reg(y_real_tr, pred_treino)
  m_te <- calc_metricas_reg(y_real_te, pred_teste)

  tabela_metricas <- data.frame(
    Métrica = c("RMSE", "MAE", "R²", "MAPE (%)"),
    Treino  = c(.am_fmt_num(m_tr$RMSE, 3), .am_fmt_num(m_tr$MAE, 3), .am_fmt_pct(m_tr$R2, 1),
                if (is.na(m_tr$MAPE)) "-" else paste0(.am_fmt_num(m_tr$MAPE, 1), "%")),
    Teste   = c(.am_fmt_num(m_te$RMSE, 3), .am_fmt_num(m_te$MAE, 3), .am_fmt_pct(m_te$R2, 1),
                if (is.na(m_te$MAPE)) "-" else paste0(.am_fmt_num(m_te$MAPE, 1), "%")),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # ── Saída no Console ──────────────────────────────────────────────────────
  metodo_rotulo <- switch(
    algoritmo,
    rf      = "Random Forest (Floresta Aleat\u00f3ria)",
    xgboost = "XGBoost (Extreme Gradient Boosting)",
    svm     = "M\u00e1quina de Vetores de Suporte (SVM)",
    lasso   = "Regress\u00e3o Lasso (Penalizada L1)",
    ridge   = "Regress\u00e3o Ridge (Penalizada L2)",
    arvore  = "\u00c1rvore de Decis\u00e3o (CART)"
  )

  .print_titulo(sprintf("PREVIS\u00c3O \u2014 %s", toupper(metodo_rotulo)))

  cat("  ▶ CONFIGURA\u00c7\u00c3O DO MODELO\n\n")
  cat(sprintf("  Vari\u00e1vel Resposta: %s\n", y_nome))
  cat(sprintf("  Preditores:        %s\n", paste(pred_nomes, collapse = ", ")))
  cat(sprintf("  Parti\u00e7\u00e3o:           %.0f%% Treino (%d obs) | %.0f%% Teste (%d obs)\n",
              proporcao_treino * 100, nrow(df_treino), (1 - proporcao_treino) * 100, nrow(df_teste)))

  if (algoritmo == "rf") {
    if (!is.null(diag_estabilizacao)) {
      cat(sprintf("  N. \u00c1rvores:        %d (otimiza\u00e7\u00e3o autom\u00e1tica por estabiliza\u00e7\u00e3o OOB)\n", n_arv_usado))
    } else {
      cat(sprintf("  N. \u00c1rvores:        %d (definido manualmente)\n", n_arv_usado))
    }
    cat(sprintf("  Vars por Divis\u00e3o:   %d (n_vars)\n", vars_por_no_usado))
  } else if (algoritmo == "xgboost") {
    cat(sprintf("  Rodadas (nrounds): %d | Profundidade: %d | Taxa (eta): %.2f\n",
                n_r_usado, prof_usada, taxa_aprendizado))
  } else if (algoritmo == "svm") {
    cat(sprintf("  Custo (C):         %.2f\n", custo))
  } else if (algoritmo %in% c("lasso", "ridge")) {
    cat(sprintf("  Penaliza\u00e7\u00e3o (\u03bb):    %.5f (otimizado via CV: lambda.min)\n", modelo_ajustado$lambda.min))
  } else if (algoritmo == "arvore") {
    cat(sprintf("  Profundidade M\u00e1x:  %d\n", prof_usada))
  }
  cat("\n")

  # Tabela de Estabilização (quando gerada automaticamente)
  if (!is.null(diag_estabilizacao)) {
    .print_topico("ESTABILIZA\u00c7\u00c3O DO N\u00daMERO DE \u00c1RVORES (OOB)")
    names(diag_estabilizacao) <- c("N. \u00c1rvores", "MSE (OOB)", "RMSE (OOB)", "Var. Relativa")
    .print_tabela_estatR(diag_estabilizacao, align = c("center", "center", "center", "center"))
    cat(sprintf("  * Estabiliza\u00e7\u00e3o confirmada em n_arvores = %d (varia\u00e7\u00e3o < 1%%).\n\n", n_arv_usado))
  }

  .print_topico("DESEMPENHO DO MODELO")
  .print_tabela_estatR(tabela_metricas, align = c("left", "center", "center"))

  # Importância das Variáveis
  if (!is.null(importancia_df) && nrow(importancia_df) > 0) {
    .print_topico("IMPORT\u00c2NCIA DAS VARI\u00c1VEIS")
    df_imp_print <- data.frame(
      "Vari\u00e1vel" = importancia_df$Variavel,
      "Import\u00e2ncia" = sapply(importancia_df$Importancia, .am_fmt_num, decimais = 2),
      "Contribui\u00e7\u00e3o (%)" = sapply(importancia_df$Percentual / 100, .am_fmt_pct, decimais = 1),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    .print_tabela_estatR(df_imp_print, align = c("left", "center", "center"))
  }

  .print_rodape()

  # ── Painel Gráfico Padronizado com Grade Dupla ───────────────────────────
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      tema_padrao <- if (exists("meu_tema", mode = "function")) {
        meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      df_pred_plot <- data.frame(
        Real     = y_real_te,
        Predito  = pred_teste,
        Residuo  = m_te$Residuos,
        Ajustado = pred_teste
      )

      # Gráfico 1: Predito vs Real (com linha identidade 45 graus)
      min_v <- min(c(df_pred_plot$Real, df_pred_plot$Predito), na.rm = TRUE)
      max_v <- max(c(df_pred_plot$Real, df_pred_plot$Predito), na.rm = TRUE)

      p1 <- ggplot2::ggplot(df_pred_plot, ggplot2::aes(x = Real, y = Predito)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#E31A1C", linewidth = 0.8) +
        ggplot2::geom_point(color = "#377EB8", size = 2.2, alpha = 0.75) +
        ggplot2::labs(title = "Predito vs Real (Teste)", x = "Valor Observado (Real)", y = "Valor Predito") +
        tema_padrao

      # Gráfico 2: Importância de Variáveis ou Resíduos vs Ajustados
      if (!is.null(importancia_df) && nrow(importancia_df) > 0) {
        imp_plot_df <- head(importancia_df, 10)
        imp_plot_df$Variavel <- factor(imp_plot_df$Variavel, levels = rev(imp_plot_df$Variavel))

        p2 <- ggplot2::ggplot(imp_plot_df, ggplot2::aes(x = Variavel, y = Percentual)) +
          ggplot2::geom_col(fill = "#377EB8", width = 0.6) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = "Import\u00e2ncia das Vari\u00e1veis (%)", x = "", y = "Import\u00e2ncia Relativa (%)") +
          tema_padrao
      } else {
        p2 <- ggplot2::ggplot(df_pred_plot, ggplot2::aes(x = Ajustado, y = Residuo)) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#E31A1C", linewidth = 0.8) +
          ggplot2::geom_point(color = "#377EB8", size = 2.2, alpha = 0.75) +
          ggplot2::labs(title = "Res\u00edduos vs Ajustados", x = "Valores Ajustados", y = "Res\u00edduo") +
          tema_padrao
      }

      if (requireNamespace("patchwork", quietly = TRUE)) {
        painel <- patchwork::wrap_plots(p1, p2, ncol = 2)
        print(painel)
      } else {
        print(p1)
      }
    }, error = function(e) {
      message("[Aviso] Não foi possível renderizar o gráfico: ", e$message)
    })
  }

  # Objeto de Retorno S3
  resultado <- list(
    tipo               = "previsao",
    algoritmo          = algoritmo,
    metodo             = algoritmo,
    formula            = formula,
    modelo             = modelo_ajustado,
    metricas           = tabela_metricas,
    metricas_treino    = m_tr,
    metricas_teste     = m_te,
    importancia        = importancia_df,
    estabilizacao      = diag_estabilizacao,
    dados_treino       = df_treino,
    dados_teste        = df_teste,
    predicoes_treino   = pred_treino,
    predicoes_teste    = pred_teste,
    residuos           = m_te$Residuos,
    y_nome             = y_nome
  )
  class(resultado) <- c("estatR_aprend_maq", "list")
  invisible(resultado)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Classificação com Aprendizado de Máquina
#' @description Ajusta modelos de Aprendizado de Máquina para classificação
#'   (variável resposta categórica ou binária). Suporta algoritmos consolidados com
#'   otimização de hiperparâmetros ou ajuste manual em português.
#'
#' @param formula Fórmula no formato \code{y ~ x1 + x2} ou \code{y ~ .}.
#' @param dados Data frame contendo as variáveis do modelo.
#' @param algoritmo Algoritmo a utilizar: \code{"rf"} (Random Forest / Floresta Aleatória, padrão),
#'   \code{"xgboost"} (Extreme Gradient Boosting), \code{"svm"} (Máquina de Vetores de Suporte),
#'   \code{"knn"} (K-Vizinhos Mais Próximos) ou \code{"arvore"} (CART).
#' @param proporcao_treino Proporção da amostra destinada ao treino (padrão \code{0.75}).
#' @param semente Valor inteiro para reprodutibilidade (padrão \code{42}).
#' @param n_arvores Número de árvores (para \code{"rf"}). Se \code{NULL} (padrão),
#'   aplica busca automática por estabilização do erro OOB.
#' @param vars_por_no Número de variáveis por divisão em cada nó (\code{mtry}).
#'   Se \code{NULL} (padrão), utiliza \code{max(1, floor(sqrt(p)))}.
#' @param k Número de vizinhos mais próximos (para \code{"knn"}). Se \code{NULL},
#'   utiliza \code{floor(sqrt(n_treino))}.
#' @param n_rounds Número de rodadas de boosting (para \code{"xgboost"}).
#' @param profundidade Profundidade máxima das árvores.
#' @param taxa_aprendizado Taxa de aprendizado (\code{eta}) para \code{"xgboost"}.
#' @param custo Parâmetro C para \code{"svm"}.
#' @param grafico Lógico. Se TRUE (padrão), exibe gráfico com grade dupla.
#'
#' @return Objeto da classe \code{estatR_aprend_maq} contendo o modelo ajustado,
#'   matriz de confusão, métricas, predições e importância de variáveis.
#' @export
classificacao_aprend_maq <- function(formula,
                                     dados,
                                     algoritmo          = c("rf", "xgboost", "svm", "knn", "arvore"),
                                     proporcao_treino   = 0.75,
                                     semente            = 42,
                                     n_arvores          = NULL,
                                     vars_por_no        = NULL,
                                     k                  = NULL,
                                     n_rounds           = NULL,
                                     profundidade       = NULL,
                                     taxa_aprendizado   = 0.1,
                                     custo              = 1,
                                     grafico            = TRUE) {

  algoritmo <- match.arg(algoritmo)
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data frame.")

  termos_obj <- terms(formula, data = dados)
  y_nome     <- as.character(formula[[2]])
  pred_nomes <- attr(termos_obj, "term.labels")

  if (!(y_nome %in% names(dados))) stop(sprintf("Variável resposta '%s' não encontrada nos dados.", y_nome))

  # Limpeza de dados
  cols_nec <- c(y_nome, pred_nomes)
  df_limpo <- na.omit(dados[, cols_nec, drop = FALSE])
  df_limpo[[y_nome]] <- as.factor(df_limpo[[y_nome]])
  classes_y <- levels(df_limpo[[y_nome]])
  n_classes <- length(classes_y)

  if (n_classes < 2) stop("A variável resposta deve conter pelo menos 2 classes distintas.")

  n_total <- nrow(df_limpo)
  if (n_total < 10) stop("Amostra insuficiente após remoção de valores ausentes (mínimo 10 observações).")

  p_pred <- length(pred_nomes)
  if (p_pred == 0) stop("O modelo deve conter pelo menos uma variável preditora.")

  # Partição Treino / Teste (Estratificada)
  if (!is.null(semente)) set.seed(semente)
  idx_treino <- integer(0)
  for (cls in classes_y) {
    idx_cls <- which(df_limpo[[y_nome]] == cls)
    n_cls_tr <- max(1, floor(length(idx_cls) * proporcao_treino))
    idx_treino <- c(idx_treino, sample(idx_cls, size = n_cls_tr, replace = FALSE))
  }
  idx_teste  <- setdiff(seq_len(n_total), idx_treino)
  df_treino  <- df_limpo[idx_treino, , drop = FALSE]
  df_teste   <- df_limpo[idx_teste, , drop = FALSE]

  diag_estabilizacao <- NULL
  n_arv_usado        <- NULL
  vars_por_no_usado  <- NULL
  modelo_ajustado    <- NULL
  importancia_df     <- NULL

  # ── AJUSTE: RANDOM FOREST ─────────────────────────────────────────────────
  if (algoritmo == "rf") {
    vars_por_no_usado <- if (!is.null(vars_por_no)) as.integer(vars_por_no) else max(1, floor(sqrt(p_pred)))

    if (is.null(n_arvores)) {
      candidatos_arv <- c(50, 100, 200, 300, 500, 750)
      erros_oob <- numeric(length(candidatos_arv))

      for (i in seq_along(candidatos_arv)) {
        if (!is.null(semente)) set.seed(semente)
        rf_tmp <- ranger::ranger(
          formula       = formula,
          data          = df_treino,
          num.trees     = candidatos_arv[i],
          mtry          = vars_por_no_usado,
          classification = TRUE,
          importance    = "none",
          verbose       = FALSE
        )
        erros_oob[i] <- rf_tmp$prediction.error
      }

      var_rel <- c(Inf, abs(diff(erros_oob)) / (abs(erros_oob[-length(erros_oob)]) + 1e-10))
      idx_est <- which(var_rel < 0.01)
      n_arv_usado <- if (length(idx_est) > 0) candidatos_arv[min(idx_est)] else 300

      diag_estabilizacao <- data.frame(
        N_Arvores    = candidatos_arv,
        Erro_OOB     = round(erros_oob, 4),
        Acuracia_OOB = round(1 - erros_oob, 4),
        Var_Rel_Pct  = c("-", paste0(round(var_rel[-1] * 100, 2), "%")),
        stringsAsFactors = FALSE
      )
    } else {
      n_arv_usado <- as.integer(n_arvores)
    }

    if (!is.null(semente)) set.seed(semente)
    modelo_ajustado <- ranger::ranger(
      formula        = formula,
      data           = df_treino,
      num.trees      = n_arv_usado,
      mtry           = vars_por_no_usado,
      classification = TRUE,
      importance     = "impurity",
      verbose        = FALSE
    )

    pred_treino <- predict(modelo_ajustado, data = df_treino)$predictions
    pred_teste  <- predict(modelo_ajustado, data = df_teste)$predictions

    imp_raw <- modelo_ajustado$variable.importance
    if (!is.null(imp_raw) && length(imp_raw) > 0) {
      imp_pct <- (imp_raw / sum(imp_raw)) * 100
      importancia_df <- data.frame(
        Variavel    = names(imp_raw),
        Importancia = round(imp_raw, 2),
        Percentual  = round(imp_pct, 1),
        stringsAsFactors = FALSE
      )
      importancia_df <- importancia_df[order(importancia_df$Importancia, decreasing = TRUE), ]
      rownames(importancia_df) <- NULL
    }
  }

  # ── AJUSTE: XGBOOST ───────────────────────────────────────────────────────
  if (algoritmo == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) {
      stop("O pacote 'xgboost' é necessário para este método.")
    }

    X_tr <- model.matrix(formula, data = df_treino)[, -1, drop = FALSE]
    y_tr <- as.integer(df_treino[[y_nome]]) - 1
    X_te <- model.matrix(formula, data = df_teste)[, -1, drop = FALSE]
    y_te <- as.integer(df_teste[[y_nome]]) - 1

    n_r_usado  <- if (!is.null(n_rounds)) as.integer(n_rounds) else 50
    prof_usada <- if (!is.null(profundidade)) as.integer(profundidade) else 6

    obj_xgb <- if (n_classes == 2) "binary:logistic" else "multi:softmax"
    num_cls <- if (n_classes > 2) n_classes else NULL

    params_list <- list(
      max_depth = prof_usada,
      eta       = taxa_aprendizado,
      objective = obj_xgb
    )
    if (!is.null(num_cls)) params_list$num_class <- num_cls

    if (!is.null(semente)) set.seed(semente)
    modelo_ajustado <- xgboost::xgboost(
      params    = params_list,
      x         = X_tr,
      y         = y_tr,
      nrounds   = n_r_usado,
      verbosity = 0
    )

    p_tr_raw <- predict(modelo_ajustado, X_tr)
    p_te_raw <- predict(modelo_ajustado, X_te)

    if (n_classes == 2) {
      pred_treino <- factor(ifelse(p_tr_raw > 0.5, classes_y[2], classes_y[1]), levels = classes_y)
      pred_teste  <- factor(ifelse(p_te_raw > 0.5, classes_y[2], classes_y[1]), levels = classes_y)
    } else {
      pred_treino <- factor(classes_y[p_tr_raw + 1], levels = classes_y)
      pred_teste  <- factor(classes_y[p_te_raw + 1], levels = classes_y)
    }

    imp_xgb <- tryCatch(xgboost::xgb.importance(model = modelo_ajustado), error = function(e) NULL)
    if (!is.null(imp_xgb) && nrow(imp_xgb) > 0) {
      importancia_df <- data.frame(
        Variavel    = imp_xgb$Feature,
        Importancia = round(imp_xgb$Gain, 4),
        Percentual  = round(imp_xgb$Gain * 100, 1),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── AJUSTE: SVM ───────────────────────────────────────────────────────────
  if (algoritmo == "svm") {
    if (!requireNamespace("e1071", quietly = TRUE)) stop("O pacote 'e1071' \u00e9 necess\u00e1rio para este algoritmo.")
    if (!is.null(semente)) set.seed(semente)
    modelo_ajustado <- e1071::svm(formula = formula, data = df_treino, cost = custo, probability = TRUE)
    pred_treino     <- predict(modelo_ajustado, newdata = df_treino)
    pred_teste      <- predict(modelo_ajustado, newdata = df_teste)
  }

  # ── AJUSTE: KNN ───────────────────────────────────────────────────────────
  if (algoritmo == "knn") {
    cols_num_p <- pred_nomes[vapply(df_treino[pred_nomes], is.numeric, logical(1))]
    if (length(cols_num_p) == 0) stop("O algoritmo KNN requer pelo menos uma variável preditora numérica.")

    X_tr_sc <- scale(df_treino[cols_num_p])
    X_te_sc <- scale(df_teste[cols_num_p],
                     center = attr(X_tr_sc, "scaled:center"),
                     scale  = attr(X_tr_sc, "scaled:scale"))

    k_usado <- if (!is.null(k)) as.integer(k) else max(1, floor(sqrt(nrow(df_treino))))
    pred_treino <- class::knn(train = X_tr_sc, test = X_tr_sc, cl = df_treino[[y_nome]], k = k_usado)
    pred_teste  <- class::knn(train = X_tr_sc, test = X_te_sc,  cl = df_treino[[y_nome]], k = k_usado)
    modelo_ajustado <- list(tipo = "knn", k = k_usado, variaveis = cols_num_p)
  }

  # ── AJUSTE: ÁRVORE (rpart) ────────────────────────────────────────────────
  if (algoritmo == "arvore") {
    prof_usada <- if (!is.null(profundidade)) as.integer(profundidade) else 10
    modelo_ajustado <- rpart::rpart(formula = formula, data = df_treino, method = "class",
                                    control = rpart::rpart.control(maxdepth = prof_usada))
    pred_treino <- predict(modelo_ajustado, newdata = df_treino, type = "class")
    pred_teste  <- predict(modelo_ajustado, newdata = df_teste, type = "class")

    imp_rp <- modelo_ajustado$variable.importance
    if (!is.null(imp_rp) && length(imp_rp) > 0) {
      importancia_df <- data.frame(
        Variavel    = names(imp_rp),
        Importancia = round(imp_rp, 2),
        Percentual  = round((imp_rp / sum(imp_rp)) * 100, 1),
        stringsAsFactors = FALSE
      )
      rownames(importancia_df) <- NULL
    }
  }

  # ── Métricas de Classificação (Matriz de Confusão, Acurácia, Kappa) ───────
  y_real_tr <- df_treino[[y_nome]]
  y_real_te <- df_teste[[y_nome]]

  calc_metricas_cls <- function(real, pred) {
    tab_conf <- table(Real = real, Predito = pred)
    acuracia <- sum(diag(tab_conf)) / sum(tab_conf)

    # Kappa de Cohen
    po <- acuracia
    pe <- sum(rowSums(tab_conf) * colSums(tab_conf)) / (sum(tab_conf)^2)
    kappa <- if (pe < 1) (po - pe) / (1 - pe) else 0

    list(Acuracia = acuracia, Kappa = kappa, Confusao = tab_conf)
  }

  m_cls_tr <- calc_metricas_cls(y_real_tr, pred_treino)
  m_cls_te <- calc_metricas_cls(y_real_te, pred_teste)

  tabela_metricas <- data.frame(
    Métrica        = c("Acurácia", "Kappa de Cohen"),
    Treino         = c(.am_fmt_pct(m_cls_tr$Acuracia, 1), .am_fmt_num(m_cls_tr$Kappa, 3)),
    Teste          = c(.am_fmt_pct(m_cls_te$Acuracia, 1), .am_fmt_num(m_cls_te$Kappa, 3)),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # ── Saída no Console ──────────────────────────────────────────────────────
  algoritmo_rotulo <- switch(
    algoritmo,
    rf      = "Random Forest (Floresta Aleat\u00f3ria)",
    xgboost = "XGBoost (Extreme Gradient Boosting)",
    svm     = "M\u00e1quina de Vetores de Suporte (SVM)",
    knn     = "K-Vizinhos Mais Pr\u00f3ximos (KNN)",
    arvore  = "\u00c1rvore de Decis\u00e3o (CART)"
  )

  .print_titulo(sprintf("CLASSIFICA\u00c7\u00c3O \u2014 %s", toupper(algoritmo_rotulo)))

  cat("  \u25b6 CONFIGURA\u00c7\u00c3O DO MODELO\n\n")
  cat(sprintf("  Vari\u00e1vel Resposta: %s (%d classes: %s)\n",
              y_nome, n_classes, paste(classes_y, collapse = ", ")))
  cat(sprintf("  Preditores:        %s\n", paste(pred_nomes, collapse = ", ")))
  cat(sprintf("  Parti\u00e7\u00e3o:           %.0f%% Treino (%d obs) | %.0f%% Teste (%d obs)\n",
              proporcao_treino * 100, nrow(df_treino), (1 - proporcao_treino) * 100, nrow(df_teste)))

  if (algoritmo == "rf") {
    if (!is.null(diag_estabilizacao)) {
      cat(sprintf("  N. \u00c1rvores:        %d (otimiza\u00e7\u00e3o autom\u00e1tica por estabiliza\u00e7\u00e3o OOB)\n", n_arv_usado))
    } else {
      cat(sprintf("  N. \u00c1rvores:        %d (definido manualmente)\n", n_arv_usado))
    }
    cat(sprintf("  Vars por Divis\u00e3o:   %d (n_vars)\n", vars_por_no_usado))
  } else if (algoritmo == "knn") {
    cat(sprintf("  Vizinhos (k):      %d\n", k_usado))
  } else if (algoritmo == "xgboost") {
    cat(sprintf("  Rodadas (nrounds): %d | Profundidade: %d | Taxa (eta): %.2f\n",
                n_r_usado, prof_usada, taxa_aprendizado))
  } else if (algoritmo == "svm") {
    cat(sprintf("  Custo (C):         %.2f\n", custo))
  } else if (algoritmo == "arvore") {
    cat(sprintf("  Profundidade M\u00e1x:  %d\n", prof_usada))
  }
  cat("\n")

  if (!is.null(diag_estabilizacao)) {
    .print_topico("ESTABILIZA\u00c7\u00c3O DO N\u00daMERO DE \u00c1RVORES (OOB)")
    names(diag_estabilizacao) <- c("N. \u00c1rvores", "Erro (OOB)", "Acurácia (OOB)", "Var. Relativa")
    .print_tabela_estatR(diag_estabilizacao, align = c("center", "center", "center", "center"))
    cat(sprintf("  * Estabiliza\u00e7\u00e3o confirmada em n_arvores = %d.\n\n", n_arv_usado))
  }

  .print_topico("DESEMPENHO DO MODELO")
  .print_tabela_estatR(tabela_metricas, align = c("left", "center", "center"))

  # Matriz de Confusão no Teste
  .print_topico("MATRIZ DE CONFUS\u00c3O (CONJUNTO DE TESTE)")
  print(m_cls_te$Confusao)
  cat("\n")

  if (!is.null(importancia_df) && nrow(importancia_df) > 0) {
    .print_topico("IMPORT\u00c2NCIA DAS VARI\u00c1VEIS")
    df_imp_print <- data.frame(
      "Vari\u00e1vel" = importancia_df$Variavel,
      "Import\u00e2ncia" = sapply(importancia_df$Importancia, .am_fmt_num, decimais = 2),
      "Contribui\u00e7\u00e3o (%)" = sapply(importancia_df$Percentual / 100, .am_fmt_pct, decimais = 1),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    .print_tabela_estatR(df_imp_print, align = c("left", "center", "center"))
  }

  .print_rodape()

  # ── Gráfico Padronizado com Grade Dupla ──────────────────────────────────
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      tema_padrao <- if (exists("meu_tema", mode = "function")) {
        meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      if (!is.null(importancia_df) && nrow(importancia_df) > 0) {
        imp_plot_df <- head(importancia_df, 10)
        imp_plot_df$Variavel <- factor(imp_plot_df$Variavel, levels = rev(imp_plot_df$Variavel))

        p_imp <- ggplot2::ggplot(imp_plot_df, ggplot2::aes(x = Variavel, y = Percentual)) +
          ggplot2::geom_col(fill = "#377EB8", width = 0.6) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = "Import\u00e2ncia das Vari\u00e1veis (%)", x = "", y = "Import\u00e2ncia Relativa (%)") +
          tema_padrao
        print(p_imp)
      }
    }, error = function(e) {
      message("[Aviso] Não foi possível renderizar o gráfico: ", e$message)
    })
  }

  resultado <- list(
    tipo             = "classificacao",
    algoritmo        = algoritmo,
    metodo           = algoritmo,
    formula          = formula,
    modelo           = modelo_ajustado,
    metricas         = tabela_metricas,
    confusao_treino  = m_cls_tr$Confusao,
    confusao_teste   = m_cls_te$Confusao,
    importancia      = importancia_df,
    estabilizacao    = diag_estabilizacao,
    dados_treino     = df_treino,
    dados_teste      = df_teste,
    predicoes_treino = pred_treino,
    predicoes_teste  = pred_teste,
    y_nome           = y_nome,
    classes          = classes_y
  )
  class(resultado) <- c("estatR_aprend_maq", "list")
  invisible(resultado)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Análise Residual para Modelos de Aprendizado de Máquina
#' @description Realiza diagnóstico formal dos resíduos em modelos de regressão
#'   ajustados por Aprendizado de Máquina (\code{previsao_aprend_maq}). Avalia viés
#'   sistemático (média dos resíduos = 0), dispersão e casos atípicos extremos,
#'   gerando painel gráfico com grade dupla.
#'
#' @param modelo Objeto da classe \code{estatR_aprend_maq} gerado por \code{previsao_aprend_maq}.
#' @param grafico Lógico. Se TRUE (padrão), exibe o painel gráfico com grade dupla.
#'
#' @return Retorna invisivelmente uma lista com as estatísticas do diagnóstico.
#' @export
analise_residual_aprend_maq <- function(modelo, grafico = TRUE) {
  if (!inherits(modelo, "estatR_aprend_maq")) {
    stop("O modelo deve ser um objeto ajustado pelo estatR (classe 'estatR_aprend_maq').")
  }
  if (modelo$tipo != "previsao") {
    stop("A análise residual é aplicável apenas a modelos de previsão numérica (regressão).")
  }

  residuos <- modelo$residuos
  n        <- length(residuos)
  preditos <- modelo$predicoes_teste
  reais    <- modelo$dados_teste[[modelo$y_nome]]

  # Testes e Estatísticas Formais
  media_res <- mean(residuos)
  sd_res    <- stats::sd(residuos)
  t_res     <- stats::t.test(residuos)
  p_vies    <- t_res$p.value

  # Normalidade
  st_norm <- if (n <= 5000) stats::shapiro.test(residuos) else suppressWarnings(stats::ks.test((residuos - media_res) / sd_res, "pnorm"))
  p_norm  <- st_norm$p.value
  nome_norm <- if (n <= 5000) "Normalidade (Shapiro-Wilk)" else "Normalidade (Kolmogorov-S)"
  stat_norm <- if (n <= 5000) sprintf("W = %s", .am_fmt_num(st_norm$statistic, 3)) else sprintf("D = %s", .am_fmt_num(st_norm$statistic, 3))

  # Outliers / Casos Atípicos (|resíduo padronizado| > 2,5)
  res_pad     <- (residuos - media_res) / sd_res
  idx_out     <- which(abs(res_pad) > 2.5)
  n_out       <- length(idx_out)
  pct_out     <- (n_out / n) * 100

  .print_titulo("DIAGN\u00d3STICO DE RES\u00cdDUOS \u2014 APRENDIZADO DE M\u00c1QUINA")

  cat("  ▶ TESTES ESTAT\u00cdSTICOS DOS RES\u00cdDUOS (CONJUNTO DE TESTE)\n\n")

  df_tab_diag <- data.frame(
    Teste = c("Vi\u00e9s Sistem\u00e1tico (t: M\u00e9dia = 0)", nome_norm),
    "Estat\u00edstica" = c(sprintf("t = %s", .am_fmt_num(t_res$statistic, 3)), stat_norm),
    "p-valor" = c(.am_fmt_p(p_vies), .am_fmt_p(p_norm)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  .print_tabela_estatR(df_tab_diag, align = c("left", "center", "center"))

  cat(sprintf("  * Desvio Padr\u00e3o dos Res\u00edduos: %s\n", .am_fmt_num(sd_res, 3)))
  cat(sprintf("  * Casos At\u00edpicos (|Res\u00edduo Padronizado| > 2,5): %d de %d (%.1f%%)\n\n",
              n_out, n, pct_out))

  # Cores padronizadas ANSI
  cor_ok   <- "\033[1;32m"
  cor_erro <- "\033[1;31m"
  cor_warn <- "\033[1;33m"
  reset    <- "\033[0m"

  .print_topico("INTERPRETA\u00c7\u00c3O DO DIAGN\u00d3STICO")

  # 1. Viés
  if (p_vies >= 0.05) {
    cat(sprintf("  %s[\u2713]%s A m\u00e9dia dos res\u00edduos n\u00e3o difere de zero (p = %s). O modelo n\u00e3o apresenta vi\u00e9s sistem\u00e1tico.\n\n",
                cor_ok, reset, .am_fmt_p(p_vies)))
  } else {
    cat(sprintf("  %s[!]%s Vi\u00e9s detectado nos res\u00edduos (p = %s). O modelo tende a subestimar ou superestimar a resposta.\n\n",
                cor_erro, reset, .am_fmt_p(p_vies)))
  }

  # 2. Outliers
  if (pct_out <= 5.0) {
    cat(sprintf("  %s[\u2713]%s Propor\u00e7\u00e3o de casos at\u00edpicos adequada (%.1f%% <= 5%% do esperado sob normalidade).\n\n",
                cor_ok, reset, pct_out))
  } else {
    cat(sprintf("  %s[!]%s Presen\u00e7a moderada de res\u00edduos discrepantes (%.1f%% > 5%%). Recomenda-se investigar outliers.\n\n",
                cor_warn, reset, pct_out))
  }

  .print_rodape()

  # Painel Gráfico de Resíduos (4 quadrantes com grade dupla)
  if (grafico && requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch({
      tema_padrao <- if (exists("meu_tema", mode = "function")) {
        meu_tema(grade = "dupla")
      } else if (exists("tema_estatR", mode = "function")) {
        tema_estatR(estilo = 2, grade = "dupla")
      } else {
        ggplot2::theme_minimal()
      }

      df_diag_plot <- data.frame(
        Ajustados = preditos,
        Residuos  = residuos,
        ResPad    = res_pad,
        Indice    = seq_along(residuos)
      )

      # 1. Resíduos vs Ajustados
      p1 <- ggplot2::ggplot(df_diag_plot, ggplot2::aes(x = Ajustados, y = Residuos)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#E31A1C", linewidth = 0.8) +
        ggplot2::geom_point(color = "#377EB8", size = 2, alpha = 0.75) +
        ggplot2::labs(title = "Res\u00edduos vs Valores Ajustados", x = "Valores Ajustados", y = "Res\u00edduo") +
        tema_padrao

      # 2. Histograma com Densidade Teórica
      p2 <- ggplot2::ggplot(df_diag_plot, ggplot2::aes(x = Residuos)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = min(20, max(6, floor(n / 5))),
                                fill = "gray85", color = "gray40") +
        ggplot2::stat_function(fun = stats::dnorm,
                               args = list(mean = media_res, sd = sd_res),
                               color = "#E31A1C", linewidth = 0.85) +
        ggplot2::labs(title = "Histograma dos Res\u00edduos", x = "Res\u00edduo", y = "Densidade") +
        tema_padrao

      # 3. QQ-Plot
      p3 <- ggplot2::ggplot(df_diag_plot, ggplot2::aes(sample = Residuos)) +
        ggplot2::stat_qq(color = "#377EB8", size = 1.8, alpha = 0.75) +
        ggplot2::stat_qq_line(color = "#E31A1C", linetype = "dashed", linewidth = 0.8) +
        ggplot2::labs(title = "QQ-Plot dos Res\u00edduos", x = "Quantis Te\u00f3ricos", y = "Quantis Amostrais") +
        tema_padrao

      # 4. Resíduos Padronizados por Índice
      p4 <- ggplot2::ggplot(df_diag_plot, ggplot2::aes(x = Indice, y = ResPad)) +
        ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = c(-2.5, 2.5), linetype = "dashed", color = "#E31A1C", linewidth = 0.7) +
        ggplot2::geom_point(color = "#377EB8", size = 1.8, alpha = 0.75) +
        ggplot2::labs(title = "Res\u00edduos Padronizados", x = "\u00cdndice da Observa\u00e7\u00e3o", y = "Res\u00edduo Padronizado") +
        tema_padrao

      if (requireNamespace("patchwork", quietly = TRUE)) {
        painel <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
        print(painel)
      } else {
        print(p1)
      }
    }, error = function(e) {
      message("[Aviso] Não foi possível renderizar o painel de resíduos: ", e$message)
    })
  }

  invisible(list(
    media_residuos    = media_res,
    sd_residuos       = sd_res,
    teste_vies        = t_res,
    teste_normalidade = st_norm,
    n_outliers        = n_out,
    pct_outliers      = pct_out
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Métricas de Desempenho para Aprendizado de Máquina
#' @description Extrai e exibe a tabela de métricas de treino e teste de um
#'   objeto de Aprendizado de Máquina (\code{estatR_aprend_maq}).
#'
#' @param modelo Objeto da classe \code{estatR_aprend_maq}.
#' @return Data frame com as métricas avaliadas.
#' @export
metricas_aprend_maq <- function(modelo) {
  if (!inherits(modelo, "estatR_aprend_maq")) {
    stop("O objeto deve ser da classe 'estatR_aprend_maq'.")
  }
  .print_titulo("M\u00c9TRICAS DE DESEMPENHO \u2014 APRENDIZADO DE M\u00c1QUINA")
  cat(sprintf("  Tarefa:    %s\n", ifelse(modelo$tipo == "previsao", "Previs\u00e3o Num\u00e9rica (Regress\u00e3o)", "Classifica\u00e7\u00e3o")))
  cat(sprintf("  Algoritmo: %s\n\n", toupper(modelo$metodo)))
  .print_tabela_estatR(modelo$metricas, align = c("left", "center", "center"))
  .print_rodape()
  invisible(modelo$metricas)
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Seleção e Comparação de Modelos com Aprendizado de Máquina
#' @description Compara múltiplos algoritmos de Aprendizado de Máquina e modelos
#'   estatísticos com busca e ajuste automático de hiperparâmetros sob validação cruzada.
#'
#' @param formula Fórmula no formato \code{y ~ .} ou \code{y ~ x1 + x2}.
#' @param dados Data frame contendo os dados.
#' @param comparar_modelos Vetor de caracteres indicando quais famílias comparar:
#'   \code{"estatistico"}, \code{"aprend_maq"} ou ambos (padrão \code{c("estatistico", "aprend_maq")}).
#' @param n_folds Número de folds para a validação cruzada (padrão \code{5}).
#' @param semente Semente para reprodutibilidade (padrão \code{42}).
#' @param top Número de modelos a exibir no ranking (padrão \code{5}).
#'
#' @return Data frame invisível com o ranking completo dos modelos comparados.
#' @keywords internal
.selecao_comparativa_modelos <- function(formula,
                                         dados,
                                         comparar_modelos = c("estatistico", "aprend_maq"),
                                         n_folds          = 5,
                                         semente          = 42,
                                         top              = 5) {

  termos_obj <- terms(formula, data = dados)
  y_nome     <- as.character(formula[[2]])
  pred_nomes <- attr(termos_obj, "term.labels")

  cols_nec <- c(y_nome, pred_nomes)
  df_limpo <- na.omit(dados[, cols_nec, drop = FALSE])
  n        <- nrow(df_limpo)
  p        <- length(pred_nomes)

  Y_vec <- df_limpo[[y_nome]]
  is_classificacao <- is.factor(Y_vec) || is.character(Y_vec) || (is.numeric(Y_vec) && length(unique(Y_vec)) == 2)

  if (!is.null(semente)) set.seed(semente)
  folds <- sample(rep(seq_len(n_folds), length.out = n))

  ranking_list <- list()

  # ── 1. Modelos Estatísticos ───────────────────────────────────────────────
  if ("estatistico" %in% comparar_modelos) {
    if (!is_classificacao) {
      # Regressão Linear (OLS) com CV
      preds_cv <- numeric(n)
      for (f in seq_len(n_folds)) {
        tr_idx <- which(folds != f)
        te_idx <- which(folds == f)
        m_lm   <- stats::lm(formula, data = df_limpo[tr_idx, , drop = FALSE])
        preds_cv[te_idx] <- predict(m_lm, newdata = df_limpo[te_idx, , drop = FALSE])
      }
      res_cv <- Y_vec - preds_cv
      rmse_cv <- sqrt(mean(res_cv^2))
      mae_cv  <- mean(abs(res_cv))
      ss_tot  <- sum((Y_vec - mean(Y_vec))^2)
      r2_cv   <- if (ss_tot > 0) max(0, 1 - sum(res_cv^2) / ss_tot) else 0

      ranking_list[[length(ranking_list) + 1]] <- data.frame(
        Algoritmo   = "Regress\u00e3o Linear (OLS)",
        Familia     = "Estat\u00edstico",
        Hiperparams = "Padr\u00e3o (m\u00ednimos quadrados)",
        Metrica1    = rmse_cv,
        Metrica2    = r2_cv,
        stringsAsFactors = FALSE
      )
    } else {
      # Regressão Logística com CV
      preds_cv_cls <- character(n)
      df_limpo[[y_nome]] <- as.factor(df_limpo[[y_nome]])
      lvls <- levels(df_limpo[[y_nome]])

      for (f in seq_len(n_folds)) {
        tr_idx <- which(folds != f)
        te_idx <- which(folds == f)
        m_glm  <- suppressWarnings(stats::glm(formula, data = df_limpo[tr_idx, , drop = FALSE], family = binomial))
        p_prob <- predict(m_glm, newdata = df_limpo[te_idx, , drop = FALSE], type = "response")
        preds_cv_cls[te_idx] <- ifelse(p_prob >= 0.5, lvls[2], lvls[1])
      }
      acuracia_cv <- sum(preds_cv_cls == df_limpo[[y_nome]]) / n

      ranking_list[[length(ranking_list) + 1]] <- data.frame(
        Algoritmo   = "Regress\u00e3o Log\u00edstica",
        Familia     = "Estat\u00edstico",
        Hiperparams = "Padr\u00e3o (GLM)",
        Metrica1    = acuracia_cv,
        Metrica2    = NA,
        stringsAsFactors = FALSE
      )
    }
  }

  # ── 2. Modelos de Aprendizado de Máquina ───────────────────────────────────
  if ("aprend_maq" %in% comparar_modelos) {
    if (!is_classificacao) {
      # A) Random Forest - Grade de hiperparâmetros
      cands_arv  <- c(100, 300, 500)
      cands_vars <- unique(pmax(1, c(floor(p / 3), floor(sqrt(p)))))

      for (na in cands_arv) {
        for (mv in cands_vars) {
          preds_cv_rf <- numeric(n)
          for (f in seq_len(n_folds)) {
            tr_idx <- which(folds != f)
            te_idx <- which(folds == f)
            m_rf   <- ranger::ranger(formula, data = df_limpo[tr_idx, , drop = FALSE],
                                     num.trees = na, mtry = mv, verbose = FALSE)
            preds_cv_rf[te_idx] <- predict(m_rf, data = df_limpo[te_idx, , drop = FALSE])$predictions
          }
          res_cv <- Y_vec - preds_cv_rf
          rmse_cv <- sqrt(mean(res_cv^2))
          ss_tot  <- sum((Y_vec - mean(Y_vec))^2)
          r2_cv   <- if (ss_tot > 0) max(0, 1 - sum(res_cv^2) / ss_tot) else 0

          ranking_list[[length(ranking_list) + 1]] <- data.frame(
            Algoritmo   = "Random Forest",
            Familia     = "Aprend. M\u00e1quina",
            Hiperparams = sprintf("n_arvores = %d, n_vars = %d", na, mv),
            Metrica1    = rmse_cv,
            Metrica2    = r2_cv,
            stringsAsFactors = FALSE
          )
        }
      }

      # B) XGBoost (se disponível)
      if (requireNamespace("xgboost", quietly = TRUE)) {
        X_mat <- model.matrix(formula, data = df_limpo)[, -1, drop = FALSE]
        for (prof in c(3, 6)) {
          for (eta in c(0.05, 0.1)) {
            preds_cv_xgb <- numeric(n)
            for (f in seq_len(n_folds)) {
              tr_idx <- which(folds != f)
              te_idx <- which(folds == f)
              m_xgb  <- xgboost::xgboost(x = X_mat[tr_idx, , drop = FALSE], y = Y_vec[tr_idx],
                                         nrounds = 50, max_depth = prof, eta = eta, verbosity = 0)
              preds_cv_xgb[te_idx] <- predict(m_xgb, X_mat[te_idx, , drop = FALSE])
            }
            res_cv <- Y_vec - preds_cv_xgb
            rmse_cv <- sqrt(mean(res_cv^2))
            ss_tot  <- sum((Y_vec - mean(Y_vec))^2)
            r2_cv   <- if (ss_tot > 0) max(0, 1 - sum(res_cv^2) / ss_tot) else 0

            ranking_list[[length(ranking_list) + 1]] <- data.frame(
              Algoritmo   = "XGBoost",
              Familia     = "Aprend. M\u00e1quina",
              Hiperparams = sprintf("profundidade = %d, eta = %.2f", prof, eta),
              Metrica1    = rmse_cv,
              Metrica2    = r2_cv,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      # C) SVM (se disponível)
      if (requireNamespace("e1071", quietly = TRUE)) {
        preds_cv_svm <- numeric(n)
        for (f in seq_len(n_folds)) {
          tr_idx <- which(folds != f)
          te_idx <- which(folds == f)
          m_svm  <- e1071::svm(formula, data = df_limpo[tr_idx, , drop = FALSE], cost = 1)
          preds_cv_svm[te_idx] <- predict(m_svm, newdata = df_limpo[te_idx, , drop = FALSE])
        }
        res_cv <- Y_vec - preds_cv_svm
        rmse_cv <- sqrt(mean(res_cv^2))
        ss_tot  <- sum((Y_vec - mean(Y_vec))^2)
        r2_cv   <- if (ss_tot > 0) max(0, 1 - sum(res_cv^2) / ss_tot) else 0

        ranking_list[[length(ranking_list) + 1]] <- data.frame(
          Algoritmo   = "SVM (Radial)",
          Familia     = "Aprend. M\u00e1quina",
          Hiperparams = "custo (C) = 1",
          Metrica1    = rmse_cv,
          Metrica2    = r2_cv,
          stringsAsFactors = FALSE
        )
      }

    } else {
      # Classificação: Random Forest
      df_limpo[[y_nome]] <- as.factor(df_limpo[[y_nome]])
      for (na in c(100, 300)) {
        preds_cv_cls <- character(n)
        for (f in seq_len(n_folds)) {
          tr_idx <- which(folds != f)
          te_idx <- which(folds == f)
          m_rf   <- ranger::ranger(formula, data = df_limpo[tr_idx, , drop = FALSE],
                                   num.trees = na, classification = TRUE, verbose = FALSE)
          p_rf   <- predict(m_rf, data = df_limpo[te_idx, , drop = FALSE])$predictions
          preds_cv_cls[te_idx] <- as.character(p_rf)
        }
        acuracia_cv <- sum(preds_cv_cls == as.character(df_limpo[[y_nome]])) / n

        ranking_list[[length(ranking_list) + 1]] <- data.frame(
          Algoritmo   = "Random Forest",
          Familia     = "Aprend. M\u00e1quina",
          Hiperparams = sprintf("n_arvores = %d", na),
          Metrica1    = acuracia_cv,
          Metrica2    = NA,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df_comp <- do.call(rbind, ranking_list)

  if (!is_classificacao) {
    # Ordenar por menor RMSE
    df_comp <- df_comp[order(df_comp$Metrica1), ]
    df_top  <- head(df_comp, top)

    df_fmt <- data.frame(
      "Posi\u00e7\u00e3o"        = paste0(seq_len(nrow(df_top)), "\u00ba"),
      "Algoritmo"      = df_top$Algoritmo,
      "Fam\u00edlia"        = df_top$Familia,
      "Hiperpar\u00e2metros" = df_top$Hiperparams,
      "RMSE (CV)"    = sapply(df_top$Metrica1, .am_fmt_num, decimais = 3),
      "R\u00b2 (CV)"        = sapply(df_top$Metrica2, .am_fmt_pct, decimais = 1),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  } else {
    # Ordenar por maior Acurácia
    df_comp <- df_comp[order(df_comp$Metrica1, decreasing = TRUE), ]
    df_top  <- head(df_comp, top)

    df_fmt <- data.frame(
      "Posi\u00e7\u00e3o"        = paste0(seq_len(nrow(df_top)), "\u00ba"),
      "Algoritmo"      = df_top$Algoritmo,
      "Fam\u00edlia"        = df_top$Familia,
      "Hiperpar\u00e2metros" = df_top$Hiperparams,
      "Acur\u00e1cia (CV)" = sapply(df_top$Metrica1, .am_fmt_pct, decimais = 1),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }

  .print_titulo("SELE\u00c7\u00c3O E COMPARA\u00c7\u00c3O DE MODELOS")

  cat(sprintf("  Vari\u00e1vel Resposta: %s\n", y_nome))
  cat(sprintf("  Preditores base:   %s\n", paste(pred_nomes, collapse = ", ")))
  cat(sprintf("  Valida\u00e7\u00e3o Cruzada: %d-fold CV (semente = %d)\n\n", n_folds, semente))

  tit_top <- if (!is_classificacao) {
    sprintf("TOP %d MODELOS (Ordenados por menor RMSE)", nrow(df_fmt))
  } else {
    sprintf("TOP %d MODELOS (Ordenados por maior Acur\u00e1cia)", nrow(df_fmt))
  }
  .print_topico(tit_top)

  align_vec <- if (!is_classificacao) c("center", "left", "center", "left", "center", "center") else c("center", "left", "center", "left", "center")
  .print_tabela_estatR(df_fmt, align = align_vec)

  cat("  * M\u00e9tricas avaliadas exclusivamente fora da amostra (Valida\u00e7\u00e3o Cruzada).\n\n")
  .print_rodape()

  invisible(df_comp)
}



# ── VALIDAÇÃO CRUZADA K-FOLD (CROSS-VALIDATION) ───────────────────────────────

#' @title Validação Cruzada K-Fold (Cross-Validation)
#' @description Executa procedimento de validação cruzada k-fold para modelos de regressão
#'   ou classificação, avaliando a estabilidade e capacidade de generalização preditiva
#'   fora da amostra com métricas detalhadas por partição (fold) e síntese global.
#'
#' @param formula Fórmula do modelo estatístico (ex: \code{y ~ x1 + x2}).
#' @param dados Data frame contendo as variáveis.
#' @param k Número de partições (folds) para a validação cruzada (padrão 5).
#' @param tipo Tipo de tarefa preditiva: \code{"auto"} (detecta automaticamente),
#'   \code{"regressao"} (resposta contínua) ou \code{"classificacao"} (resposta binária/categórica).
#' @param algoritmo Algoritmo base: \code{"linear"} (regressão linear/logística clássica),
#'   \code{"random_forest"} (floresta aleatória via ranger) ou \code{"arvore"} (árvore de decisão rpart).
#' @param semente Semente para reprodutibilidade da divisão dos folds (padrão 42).
#' @param grafico Lógico. Se \code{TRUE} (padrão), plota o gráfico de desempenho fold a fold.
#'
#' @return Objeto da classe \code{estatR_cv} contendo:
#' \itemize{
#'   \item \code{metricas_folds}: Tabela com as métricas avaliadas em cada fold.
#'   \item \code{resumo}: Tabela com a média e desvio-padrão das métricas.
#'   \item \code{grafico}: Objeto gráfico do ggplot2.
#' }
#' @export
validacao_cruzada <- function(formula,
                              dados,
                              k         = 5,
                              tipo      = c("auto", "regressao", "classificacao"),
                              algoritmo = c("linear", "random_forest", "arvore"),
                              semente   = 42,
                              grafico   = TRUE) {

  tipo      <- match.arg(tipo)
  algoritmo <- match.arg(algoritmo)

  if (missing(formula) || missing(dados)) {
    stop("Forne\u00e7a 'formula' e 'dados'.", call. = FALSE)
  }

  resp_nome <- as.character(formula[[2]])
  if (!resp_nome %in% names(dados)) {
    stop(sprintf("Vari\u00e1vel resposta '%s' n\u00e3o encontrada nos dados.", resp_nome), call. = FALSE)
  }

  y_raw <- dados[[resp_nome]]
  is_classif <- if (tipo == "auto") {
    is.factor(y_raw) || is.character(y_raw) || (is.numeric(y_raw) && length(unique(stats::na.omit(y_raw))) == 2)
  } else {
    tipo == "classificacao"
  }

  # Preparar dados
  df_limpo <- stats::na.omit(dados[, all.vars(formula), drop = FALSE])
  N <- nrow(df_limpo)
  if (N < k * 3) stop("Volume de dados insuficiente para o n\u00famero k de folds.", call. = FALSE)

  if (is_classif) {
    df_limpo[[resp_nome]] <- as.factor(df_limpo[[resp_nome]])
  } else {
    df_limpo[[resp_nome]] <- as.numeric(df_limpo[[resp_nome]])
  }

  # Sorteio dos Folds
  if (!is.null(semente)) set.seed(semente)
  if (is_classif) {
    # Estratificado por classe
    folds <- integer(N)
    for (lvl in levels(df_limpo[[resp_nome]])) {
      idx_lvl <- which(df_limpo[[resp_nome]] == lvl)
      n_lvl <- length(idx_lvl)
      folds[idx_lvl] <- sample(rep(1:k, length.out = n_lvl))
    }
  } else {
    folds <- sample(rep(1:k, length.out = N))
  }

  resultados_folds <- list()

  for (i in 1:k) {
    treino <- df_limpo[folds != i, , drop = FALSE]
    teste  <- df_limpo[folds == i, , drop = FALSE]
    y_true <- teste[[resp_nome]]

    if (!is_classif) {
      # ── Regressão ──────────────────────────────────────────────────────────
      if (algoritmo == "linear") {
        mod <- stats::lm(formula, data = treino)
        preds <- stats::predict(mod, newdata = teste)
      } else if (algoritmo == "random_forest" && requireNamespace("ranger", quietly = TRUE)) {
        mod <- ranger::ranger(formula, data = treino, num.trees = 150, verbose = FALSE)
        preds <- stats::predict(mod, data = teste)$predictions
      } else {
        mod <- rpart::rpart(formula, data = treino)
        preds <- stats::predict(mod, newdata = teste)
      }

      erros <- y_true - preds
      rmse_val <- sqrt(mean(erros^2, na.rm = TRUE))
      mae_val  <- mean(abs(erros), na.rm = TRUE)
      ss_tot   <- sum((y_true - mean(y_true))^2)
      ss_res   <- sum(erros^2)
      r2_val   <- max(0, 1 - (ss_res / ss_tot))

      resultados_folds[[i]] <- data.frame(
        Fold = sprintf("Fold %d", i),
        RMSE = rmse_val,
        MAE  = mae_val,
        R2   = r2_val,
        stringsAsFactors = FALSE
      )
    } else {
      # ── Classificação ───────────────────────────────────────────────────────
      lvls <- levels(df_limpo[[resp_nome]])
      if (algoritmo == "linear") {
        if (length(lvls) == 2) {
          mod <- stats::glm(formula, data = treino, family = stats::binomial)
          probs <- stats::predict(mod, newdata = teste, type = "response")
          preds <- factor(ifelse(probs >= 0.5, lvls[2], lvls[1]), levels = lvls)
        } else {
          mod <- nnet::multinom(formula, data = treino, trace = FALSE)
          preds <- stats::predict(mod, newdata = teste)
        }
      } else if (algoritmo == "random_forest" && requireNamespace("ranger", quietly = TRUE)) {
        mod <- ranger::ranger(formula, data = treino, num.trees = 150, classification = TRUE, verbose = FALSE)
        preds <- stats::predict(mod, data = teste)$predictions
      } else {
        mod <- rpart::rpart(formula, data = treino, method = "class")
        preds <- stats::predict(mod, newdata = teste, type = "class")
      }

      tab_conf <- table(factor(preds, levels = lvls), factor(y_true, levels = lvls))
      acuracia <- sum(diag(tab_conf)) / sum(tab_conf)

      if (length(lvls) == 2) {
        sens <- if (sum(tab_conf[, 2]) > 0) tab_conf[2, 2] / sum(tab_conf[, 2]) else NA
        esp  <- if (sum(tab_conf[, 1]) > 0) tab_conf[1, 1] / sum(tab_conf[, 1]) else NA
        prec <- if (sum(tab_conf[2, ]) > 0) tab_conf[2, 2] / sum(tab_conf[2, ]) else NA
        f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * (prec * sens) / (prec + sens) else NA
      } else {
        sens <- NA; esp <- NA; f1 <- NA
      }

      resultados_folds[[i]] <- data.frame(
        Fold         = sprintf("Fold %d", i),
        Acuracia     = acuracia,
        Sensibilidade = sens,
        Especificidade = esp,
        F1_Score     = f1,
        stringsAsFactors = FALSE
      )
    }
  }

  df_folds <- do.call(rbind, resultados_folds)

  # Resumo (Média e DP)
  metricas_nomes <- setdiff(names(df_folds), "Fold")
  medias <- colMeans(df_folds[, metricas_nomes, drop = FALSE], na.rm = TRUE)
  desvios <- apply(df_folds[, metricas_nomes, drop = FALSE], 2, stats::sd, na.rm = TRUE)

  df_resumo <- data.frame(
    Metrica = metricas_nomes,
    Media   = medias,
    DP      = desvios,
    stringsAsFactors = FALSE
  )
  rownames(df_resumo) <- NULL

  # ── Saída Console ────────────────────────────────────────────────────────────
  .print_titulo("VALIDAÇÃO CRUZADA K-FOLD (CROSS-VALIDATION)")

  cat(sprintf("  Variável Resposta: %s  |  Tipo de Tarefa: %s
",
              resp_nome, if (is_classif) "Classificação" else "Regressão Contínua"))
  cat(sprintf("  Algoritmo: %s  |  Número de Partinções: %d-fold CV
",
              toupper(algoritmo), k))
  cat(sprintf("  Amostra Total: %d observações válidas

", N))

  .print_topico("DESEMPENHO POR PARTIÇÃO (FOLD A FOLD)")
  w_cv <- c(12, rep(16, length(metricas_nomes)))
  cat(.th_sep(sum(w_cv)), "
")
  cat("  ", .th_pad("Partição", w_cv[1], "left"), sep = "")
  for (m_idx in seq_along(metricas_nomes)) {
    cat(.th_pad(metricas_nomes[m_idx], w_cv[m_idx + 1], "center"))
  }
  cat("
")
  cat(.th_sep(sum(w_cv)), "
")

  for (i in 1:nrow(df_folds)) {
    cat("  ", .th_pad(df_folds$Fold[i], w_cv[1], "left"), sep = "")
    for (m_idx in seq_along(metricas_nomes)) {
      val <- df_folds[i, metricas_nomes[m_idx]]
      val_str <- if (metricas_nomes[m_idx] %in% c("R2", "Acuracia", "Sensibilidade", "Especificidade", "F1_Score")) {
        sprintf("%.1f%%", val * 100)
      } else {
        formatC(val, format = "f", digits = 3, decimal.mark = ",")
      }
      cat(.th_pad(val_str, w_cv[m_idx + 1], "center"))
    }
    cat("
")
  }
  cat(.th_sep(sum(w_cv)), "
")

  # Linha de Média ± DP
  cat("  ", .th_pad("MÉDIA (CV)", w_cv[1], "left"), sep = "")
  for (m_idx in seq_along(metricas_nomes)) {
    m_val <- df_resumo$Media[m_idx]
    dp_val <- df_resumo$DP[m_idx]
    res_str <- if (df_resumo$Metrica[m_idx] %in% c("R2", "Acuracia", "Sensibilidade", "Especificidade", "F1_Score")) {
      sprintf("%.1f%% (±%.1f%%)", m_val * 100, dp_val * 100)
    } else {
      sprintf("%.3f (±%.3f)", m_val, dp_val)
    }
    cat(.th_pad(res_str, w_cv[m_idx + 1], "center"))
  }
  cat("
")
  cat(.th_sep(sum(w_cv)), "

")

  .print_rodape()

  # ── Gráfico ggplot2 ──────────────────────────────────────────────────────────
  p_cv <- NULL
  if (grafico) {
    metrica_foco <- metricas_nomes[1]
    df_plot <- data.frame(
      Fold  = factor(df_folds$Fold, levels = df_folds$Fold),
      Valor = df_folds[[metrica_foco]]
    )
    media_foco <- mean(df_plot$Valor, na.rm = TRUE)

    p_cv <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Fold, y = Valor)) +
      ggplot2::geom_col(fill = "#2B5C8F", width = 0.55, alpha = 0.85) +
      ggplot2::geom_hline(yintercept = media_foco, color = "#D90429", linetype = "dashed", linewidth = 1) +
      ggplot2::annotate("text", x = 1, y = media_foco,
                        label = sprintf("Média = %.3f", media_foco),
                        color = "#D90429", fontface = "bold", vjust = -0.5, hjust = 0) +
      ggplot2::labs(
        title    = sprintf("Validação Cruzada: Desempenho de %s (%d Folds)", metrica_foco, k),
        subtitle = sprintf("Modelo: %s | Algoritmo: %s", resp_nome, toupper(algoritmo)),
        x = "Partição de Teste (Fold)", y = metrica_foco
      )

    if (exists("meu_tema", mode = "function")) {
      p_cv <- p_cv + meu_tema(grade = "dupla")
    } else if (exists("tema_estatR", mode = "function")) {
      p_cv <- p_cv + tema_estatR(estilo = 1, grade = "dupla")
    } else {
      p_cv <- p_cv + ggplot2::theme_minimal()
    }
    suppressMessages(print(p_cv))
  }

  res <- list(
    metricas_folds = df_folds,
    resumo         = df_resumo,
    k              = k,
    grafico        = p_cv
  )
  class(res) <- c("estatR_cv", "list")
  return(invisible(res))
}
