#' @title Resumo Estatistico de Variaveis ou Banco de Dados
#' @description Gera um resumo estatistico completo.
#' @param x Um data frame ou vetor.
#' @param por Variavel categorica para agrupar o resumo.
#' @param numericas Logico. Se TRUE, exibe numericas.
#' @param categoricas Logico. Se TRUE, exibe categoricas.
#' @param decimais Casas decimais (padrao 2).
#' @export
descrever <- function(x, por = NULL, numericas = TRUE, categoricas = TRUE, decimais = 2) {
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)

  # Captura a expressão de 'por' sem avaliar prematuramente (permite usar com ou sem aspas)
  por_sub <- substitute(por)

  if (!missing(por) && !is.null(por_sub)) {
    por_expr <- deparse(por_sub)
    por_nome <- sub(".*\\$", "", por_expr)

    if (is.data.frame(x)) {
      # 1. Verifica se foi passado o nome da coluna diretamente sem aspas (ex: por = Species)
      if (por_nome %in% names(x)) {
        grupo_vec <- x[[por_nome]]
        nome_g <- por_nome
        return(.descrever_dataframe_por_grupo(x, grupo = grupo_vec, grupo_nome = nome_g, numericas = numericas, decimais = decimais))
      }

      # 2. Se não encontrou o nome direto, avalia no ambiente pai (caso seja uma string ou variável externa)
      por_val <- tryCatch(eval(por_sub, parent.frame()), error = function(e) NULL)
      if (!is.null(por_val)) {
        if (is.character(por_val) && length(por_val) == 1 && por_val %in% names(x)) {
          grupo_vec <- x[[por_val]]
          nome_g <- por_val
        } else if (length(por_val) == nrow(x)) {
          grupo_vec <- por_val
          nome_g <- por_nome
        } else {
          stop(sprintf("A coluna '%s' não foi encontrada no banco de dados.", por_nome))
        }
        return(.descrever_dataframe_por_grupo(x, grupo = grupo_vec, grupo_nome = nome_g, numericas = numericas, decimais = decimais))
      } else {
        stop(sprintf("A coluna '%s' não foi encontrada no banco de dados.", por_nome))
      }

    } else {
      # Se x for um vetor (ex: iris$Sepal.Length, por = iris$Species)
      por_val <- eval(por_sub, parent.frame())
      if (is.numeric(x)) {
        if (length(x) != length(por_val)) stop("Os vetores 'x' e 'por' devem ter o mesmo comprimento.")
        return(.descrever_vetor_numerico_por_grupo(x, grupo = por_val, nome_var = var_nome, nome_grupo = por_nome, decimais = decimais))
      } else if (is.character(x) || is.factor(x)) {
        return(tabela_contingencia(x, por_val, decimais = decimais))
      }
    }
  }

  # Caso padrão (sem agrupamento)
  if (is.data.frame(x)) {
    return(.descrever_dataframe(x, numericas = numericas, categoricas = categoricas, decimais = decimais))
  } else if (is.numeric(x)) {
    return(.descrever_vetor_numerico(x, nome = var_nome, decimais = decimais))
  } else if (is.character(x) || is.factor(x)) {
    return(.descrever_vetor_categorico(x, nome = var_nome, decimais = decimais))
  } else {
    stop("O tipo do objeto não é suportado pelo descrever(). Use um vetor numérico, categórico ou data.frame.")
  }
}

#' Função interna para descrever vetor numérico por grupo
#' @noRd
.descrever_vetor_numerico_por_grupo <- function(x, grupo, nome_var = "Variável", nome_grupo = "Grupo", decimais = 2, dentro_df = FALSE) {
  valido <- !is.na(x) & !is.na(grupo)
  x_c <- x[valido]
  g_c <- as.factor(grupo[valido])
  
  niveis <- levels(g_c)
  
  if (!dentro_df) {
    tit <- sprintf("ESTATÍSTICA DESCRITIVA POR GRUPO (%s por %s)", nome_var, nome_grupo)
    .print_header(tit)
  } else {
    cat(sprintf("\033[1m▶ %s (por %s)\033[0m\n\n", nome_var, nome_grupo))
  }
  
  w_grp <- max(14, max(nchar(niveis)), nchar(nome_grupo), nchar("Total")) + 2
  w_num <- c(media = 7, dp = 6, cv = 7, min = 6, q1 = 6, q2 = 6, q3 = 6, max = 6, assim = 6, curt = 6)
  
  col_titulos <- c(
    .pad_string(nome_grupo, w_grp, "left"),
    .pad_string("Média", w_num["media"], "center"),
    .pad_string("DP", w_num["dp"], "center"),
    .pad_string("CV(%)", w_num["cv"], "center"),
    .pad_string("Mín", w_num["min"], "center"),
    .pad_string("Q1", w_num["q1"], "center"),
    .pad_string("Q2", w_num["q2"], "center"),
    .pad_string("Q3", w_num["q3"], "center"),
    .pad_string("Máx", w_num["max"], "center"),
    .pad_string("Assim", w_num["assim"], "center"),
    .pad_string("Curt", w_num["curt"], "center")
  )
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  calc_linha <- function(vals, rotulo) {
    if (length(vals) == 0) return(NULL)
    media <- mean(vals)
    sd_val <- sd(vals)
    cv <- if (media != 0) (sd_val / media) * 100 else NA
    min_val <- min(vals)
    q1 <- quantile(vals, 0.25)
    med <- median(vals)
    q3 <- quantile(vals, 0.75)
    max_val <- max(vals)
    assim <- .calcular_assimetria(vals)
    curt <- .calcular_curtose(vals)
    
    rot_trunc <- if (nchar(rotulo) > w_grp) paste0(substr(rotulo, 1, w_grp - 2), "..") else rotulo
    cv_str <- if (!is.na(cv)) .fmt_num(cv, 1) else "-"
    
    c(
      .pad_string(rot_trunc, w_grp, "left"),
      .pad_string(.fmt_num(media, decimais), w_num["media"], "center"),
      .pad_string(.fmt_num(sd_val, decimais), w_num["dp"], "center"),
      .pad_string(cv_str, w_num["cv"], "center"),
      .pad_string(.fmt_num(min_val, decimais), w_num["min"], "center"),
      .pad_string(.fmt_num(q1, decimais), w_num["q1"], "center"),
      .pad_string(.fmt_num(med, decimais), w_num["q2"], "center"),
      .pad_string(.fmt_num(q3, decimais), w_num["q3"], "center"),
      .pad_string(.fmt_num(max_val, decimais), w_num["max"], "center"),
      .pad_string(.fmt_num(assim, decimais), w_num["assim"], "center"),
      .pad_string(.fmt_num(curt, decimais), w_num["curt"], "center")
    )
  }
  
  for (niv in niveis) {
    vals_niv <- x_c[g_c == niv]
    l_vals <- calc_linha(vals_niv, niv)
    if (!is.null(l_vals)) {
      cat("  ", paste(l_vals, collapse = " "), "\n", sep = "")
    }
  }
  
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  l_tot <- calc_linha(x_c, "Total")
  cat("  ", paste(l_tot, collapse = " "), "\n", sep = "")
  
  if (!dentro_df) {
    .print_footer()
  } else {
    cat("\n")
  }
  invisible(x)
}

#' Função interna para descrever dataframe por grupo
#' @noRd
.descrever_dataframe_por_grupo <- function(dados, grupo, grupo_nome, numericas = TRUE, decimais = 2) {
  cols_num <- names(dados)[sapply(dados, is.numeric)]
  
  if (!numericas || length(cols_num) == 0) {
    cat("Nenhuma variável numérica encontrada no banco de dados para agrupamento.\n")
    return(invisible(dados))
  }
  
  tit <- sprintf("ESTATÍSTICA DESCRITIVA POR GRUPO (%s)", grupo_nome)
  .print_header(tit)
  
  for (col in cols_num) {
    v <- dados[[col]]
    .descrever_vetor_numerico_por_grupo(v, grupo = grupo, nome_var = col, nome_grupo = grupo_nome, decimais = decimais, dentro_df = TRUE)
  }
  
  .print_footer()
  invisible(dados)
}


#' Função interna para descrever data.frame
#' @noRd
.descrever_dataframe <- function(dados, numericas = TRUE, categoricas = TRUE, decimais = 2) {
  .print_header("RESUMO ESTATÍSTICO DO BANCO DE DADOS")
  
  # Identifica tipos
  cols_num <- names(dados)[sapply(dados, is.numeric)]
  cols_cat <- names(dados)[sapply(dados, function(v) is.character(v) || is.factor(v))]
  
  # 1. VARIÁVEIS NUMÉRICAS
  if (numericas && length(cols_num) > 0) {
    cat("\033[1m▶ VARIÁVEIS NUMÉRICAS\033[0m\n\n")
    
    # Larguras compactas (sem N e sem NA)
    w_num <- c(var = 16, media = 7, dp = 6, cv = 7, min = 6, q1 = 6, q2 = 6, q3 = 6, max = 6, assim = 6, curt = 6)
    
    # Cabeçalho da Tabela (valores centralizados sob o cabeçalho)
    col_titulos <- c(
      .pad_string("Variável", w_num["var"], "left"),
      .pad_string("Média", w_num["media"], "center"),
      .pad_string("DP", w_num["dp"], "center"),
      .pad_string("CV(%)", w_num["cv"], "center"),
      .pad_string("Mín", w_num["min"], "center"),
      .pad_string("Q1", w_num["q1"], "center"),
      .pad_string("Q2", w_num["q2"], "center"),
      .pad_string("Q3", w_num["q3"], "center"),
      .pad_string("Máx", w_num["max"], "center"),
      .pad_string("Assim", w_num["assim"], "center"),
      .pad_string("Curt", w_num["curt"], "center")
    )
    hdr <- paste(col_titulos, collapse = " ")
    cat("  ", hdr, "\n", sep = "")
    cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
    
    for (col in cols_num) {
      v <- dados[[col]]
      v_clean <- v[!is.na(v)]
      
      media <- mean(v_clean)
      sd_val <- sd(v_clean)
      cv <- if (media != 0) (sd_val / media) * 100 else NA
      min_val <- min(v_clean)
      q1 <- quantile(v_clean, 0.25)
      med <- median(v_clean)
      q3 <- quantile(v_clean, 0.75)
      max_val <- max(v_clean)
      assim <- .calcular_assimetria(v_clean)
      curt <- .calcular_curtose(v_clean)
      
      nome_trunc <- if (nchar(col) > w_num["var"]) paste0(substr(col, 1, w_num["var"] - 2), "..") else col
      cv_str <- if (!is.na(cv)) .fmt_num(cv, 1) else "-"
      
      col_valores <- c(
        .pad_string(nome_trunc, w_num["var"], "left"),
        .pad_string(.fmt_num(media, decimais), w_num["media"], "center"),
        .pad_string(.fmt_num(sd_val, decimais), w_num["dp"], "center"),
        .pad_string(cv_str, w_num["cv"], "center"),
        .pad_string(.fmt_num(min_val, decimais), w_num["min"], "center"),
        .pad_string(.fmt_num(q1, decimais), w_num["q1"], "center"),
        .pad_string(.fmt_num(med, decimais), w_num["q2"], "center"),
        .pad_string(.fmt_num(q3, decimais), w_num["q3"], "center"),
        .pad_string(.fmt_num(max_val, decimais), w_num["max"], "center"),
        .pad_string(.fmt_num(assim, decimais), w_num["assim"], "center"),
        .pad_string(.fmt_num(curt, decimais), w_num["curt"], "center")
      )
      cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
    }
    cat("\n")
  }
  
  # 2. VARIÁVEIS CATEGÓRICAS
  if (categoricas && length(cols_cat) > 0) {
    cat("\033[1m▶ VARIÁVEIS CATEGÓRICAS\033[0m\n\n")
    
    w_cat <- c(var = 16, niveis = 8, moda = 20, freq = 10)
    
    col_titulos_cat <- c(
      .pad_string("Variável", w_cat["var"], "left"),
      .pad_string("Níveis", w_cat["niveis"], "center"),
      .pad_string("Moda", w_cat["moda"], "center"),
      .pad_string("Freq (%)", w_cat["freq"], "center")
    )
    hdr_cat <- paste(col_titulos_cat, collapse = " ")
    cat("  ", hdr_cat, "\n", sep = "")
    cat("  ", paste(rep("─", nchar(hdr_cat)), collapse = ""), "\n", sep = "")
    
    for (col in cols_cat) {
      v <- dados[[col]]
      v_clean <- v[!is.na(v)]
      n_obs <- length(v_clean)
      n_levels <- length(unique(v_clean))
      
      moda_val <- .calcular_moda(v_clean)
      tab <- table(v_clean)
      max_freq <- if(length(tab) > 0) max(tab) else 0
      pct_moda <- if(n_obs > 0) max_freq / n_obs else 0
      
      nome_trunc <- if (nchar(col) > w_cat["var"]) paste0(substr(col, 1, w_cat["var"] - 2), "..") else col
      moda_trunc <- if (nchar(moda_val) > w_cat["moda"]) paste0(substr(moda_val, 1, w_cat["moda"] - 2), "..") else moda_val
      
      col_valores_cat <- c(
        .pad_string(nome_trunc, w_cat["var"], "left"),
        .pad_string(as.character(n_levels), w_cat["niveis"], "center"),
        .pad_string(moda_trunc, w_cat["moda"], "center"),
        .pad_string(.fmt_pct(pct_moda, 1), w_cat["freq"], "center")
      )
      cat("  ", paste(col_valores_cat, collapse = " "), "\n", sep = "")
    }
    cat("\n")
  }
  
  .print_footer()
  invisible(dados)
}


#' Função interna para descrever um único vetor numérico
#' @noRd
.descrever_vetor_numerico <- function(x, nome = NULL, decimais = 2) {
  x_clean <- x[!is.na(x)]
  n_obs <- length(x_clean)
  
  media <- mean(x_clean)
  mediana <- median(x_clean)
  
  # Identifica se a variável é discreta ou contínua para exibir Moda ou Classe Modal
  e_discreta <- is.integer(x_clean) || all(x_clean == round(x_clean))
  
  if (e_discreta) {
    tab <- table(x_clean)
    max_freq <- max(tab)
    pct_moda <- if (n_obs > 0) max_freq / n_obs else 0
    modas <- names(tab)[tab == max_freq]
    
    if (max_freq == 1 && length(tab) > 1) {
      str_moda <- "Amodal"
    } else if (length(modas) > 3) {
      str_moda <- sprintf("Multimodal (%s)", .fmt_pct(pct_moda, 1))
    } else {
      str_moda <- sprintf("%s (%s)", paste(modas, collapse = ", "), .fmt_pct(pct_moda, 1))
    }
    linha_modal <- sprintf("  Moda:         %s\n\n", str_moda)
  } else {
    h <- hist(x_clean, plot = FALSE)
    idx_max <- which.max(h$counts)
    max_count <- h$counts[idx_max]
    pct_classe <- if (n_obs > 0) max_count / n_obs else 0
    lim_inf <- .fmt_num(h$breaks[idx_max], decimais)
    lim_sup <- .fmt_num(h$breaks[idx_max + 1], decimais)
    str_classe_modal <- sprintf("[%s; %s) (%s)", lim_inf, lim_sup, .fmt_pct(pct_classe, 1))
    linha_modal <- sprintf("  Classe Modal: %s\n\n", str_classe_modal)
  }
  
  dp <- sd(x_clean)
  var_val <- var(x_clean)
  amp <- max(x_clean) - min(x_clean)
  cv <- (dp / media) * 100
  
  q_vals <- quantile(x_clean, probs = c(0, 0.25, 0.5, 0.75, 1))
  
  assim <- .calcular_assimetria(x_clean)
  curt <- .calcular_curtose(x_clean)
  
  tit <- if (!is.null(nome) && nzchar(nome) && !identical(nome, "x")) {
    sprintf("DESCRITIVA: VARIÁVEL NUMÉRICA (%s)", nome)
  } else {
    "DESCRITIVA: VARIÁVEL NUMÉRICA"
  }
  .print_header(tit)
  
  cat("▶ Tendência Central\n")
  cat(sprintf("  Média:        %s\n", .fmt_num(media, decimais)))
  cat(sprintf("  Mediana:      %s\n", .fmt_num(mediana, decimais)))
  cat(linha_modal)
  
  cat("▶ Dispersão e Escala\n")
  cat(sprintf("  Desvio-padrão: %s\n", .fmt_num(dp, decimais)))
  cat(sprintf("  Variância:     %s\n", .fmt_num(var_val, decimais)))
  cat(sprintf("  Amplitude:     %s\n", .fmt_num(amp, decimais)))
  cat(sprintf("  CV%%:           %s%%\n\n", .fmt_num(cv, decimais)))
  
  cat("▶ Posição (Quartis)\n")
  cat(sprintf("  Mín: %s  |  Q1: %s  |  Q2: %s  |  Q3: %s  |  Máx: %s\n\n",
              .fmt_num(q_vals[1], decimais), .fmt_num(q_vals[2], decimais),
              .fmt_num(q_vals[3], decimais), .fmt_num(q_vals[4], decimais),
              .fmt_num(q_vals[5], decimais)))
  
  cat("▶ Forma da Distribuição\n")
  cat(sprintf("  Assimetria: %s\n", .fmt_num(assim, decimais)))
  cat(sprintf("  Curtose:    %s\n", .fmt_num(curt, decimais)))
  
  .print_footer()
  invisible(x)
}


#' Função interna para descrever um vetor categórico (Frequência)
#' @noRd
.descrever_vetor_categorico <- function(x, nome = NULL, decimais = 2) {
  x_clean <- x[!is.na(x)]
  n_obs <- length(x_clean)
  
  tab <- table(x_clean)
  df_freq <- as.data.frame(tab)
  colnames(df_freq) <- c("Categoria", "Frequencia")
  df_freq$Relativa <- (df_freq$Frequencia / n_obs)
  df_freq$Acumulada <- cumsum(df_freq$Relativa)
  
  tit <- if (!is.null(nome) && nzchar(nome) && !identical(nome, "x")) {
    sprintf("DESCRITIVA: VARIÁVEL CATEGÓRICA (%s)", nome)
  } else {
    "DESCRITIVA: VARIÁVEL CATEGÓRICA"
  }
  .print_header(tit)
  
  w_cat_ind <- c(cat = 22, n = 8, rel = 14, acum = 14)
  col_titulos <- c(
    .pad_string("Categoria", w_cat_ind["cat"], "left"),
    .pad_string("N", w_cat_ind["n"], "center"),
    .pad_string("Relativa (%)", w_cat_ind["rel"], "center"),
    .pad_string("Acumulada (%)", w_cat_ind["acum"], "center")
  )
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  for (i in 1:nrow(df_freq)) {
    cat_name <- as.character(df_freq$Categoria[i])
    if (nchar(cat_name) > w_cat_ind["cat"]) cat_name <- paste0(substr(cat_name, 1, w_cat_ind["cat"] - 2), "..")
    
    col_valores <- c(
      .pad_string(cat_name, w_cat_ind["cat"], "left"),
      .pad_string(as.character(df_freq$Frequencia[i]), w_cat_ind["n"], "center"),
      .pad_string(.fmt_pct(df_freq$Relativa[i], 1), w_cat_ind["rel"], "center"),
      .pad_string(.fmt_pct(df_freq$Acumulada[i], 1), w_cat_ind["acum"], "center")
    )
    cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
  }
  
  # Linha de Total
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  col_total <- c(
    .pad_string("Total", w_cat_ind["cat"], "left"),
    .pad_string(as.character(n_obs), w_cat_ind["n"], "center"),
    .pad_string(.fmt_pct(1, 1), w_cat_ind["rel"], "center"),
    .pad_string("-", w_cat_ind["acum"], "center")
  )
  cat("  ", paste(col_total, collapse = " "), "\n", sep = "")
  
  .print_footer()
  invisible(x)
}


#' Tabela de Contingência (Cruzamento Bivariado)
#'
#' Cria uma tabela de dupla entrada para cruzar duas variáveis categóricas,
#' incluindo totais marginais e opções de proporção por linha, coluna ou total geral.
#'
#' @param x Vetor categórico (linhas).
#' @param y Vetor categórico (colunas).
#' @param proporcao Tipo de proporção a exibir: "nenhuma" (apenas N), "linha" (%), "coluna" (%) ou "total" (%).
#' @param decimais Casas decimais para as porcentagens (padrão 1).
#' @export
tabela_contingencia <- function(dados, var_x = NULL, var_y = NULL, proporcao = c("nenhuma", "linha", "coluna", "total"), decimais = 1) {
  proporcao <- match.arg(proporcao)
  
  sub_dados <- substitute(dados)
  sub_x <- substitute(var_x)
  sub_y <- substitute(var_y)
  
  if (is.data.frame(dados)) {
    if (missing(var_x) || missing(var_y)) stop("Informe as variaveis x e y. Ex: tabela_contingencia(banco, var1, var2)")
    
    x_nome <- deparse(sub_x)
    y_nome <- deparse(sub_y)
    
    x <- eval(sub_x, dados, parent.frame())
    y <- eval(sub_y, dados, parent.frame())
    
    if (is.null(x)) stop(sprintf("Variavel '%s' nao encontrada.", x_nome))
    if (is.null(y)) stop(sprintf("Variavel '%s' nao encontrada.", y_nome))
  } else {
    # Comportamento antigo: vetores diretos
    x <- dados
    y <- var_x
    x_nome <- sub(".*\\$", "", deparse(sub_dados))
    y_nome <- sub(".*\\$", "", deparse(sub_x))
    if (missing(var_x) || is.null(y)) stop("Informe as duas variaveis validas.")
  }
  
  valido <- !is.na(x) & !is.na(y)
  x_c <- x[valido]
  y_c <- y[valido]
  
  tab <- table(x_c, y_c)
  tot_linha <- rowSums(tab)
  tot_col <- colSums(tab)
  tot_geral <- sum(tab)
  
  categorias_x <- rownames(tab)
  categorias_y <- colnames(tab)
  
  tit <- sprintf("TABELA DE CONTINGÊNCIA (%s x %s)", x_nome, y_nome)
  .print_header(tit)
  
  w_lin <- max(14, max(nchar(categorias_x)), nchar(x_nome), nchar("Total")) + 2
  w_col <- max(12, max(nchar(categorias_y)), nchar("Total")) + 2
  
  # Rótulo da variável de coluna (eixo Y) centralizado acima das categorias
  largura_cols_y <- length(categorias_y) * (w_col + 1) - 1
  cat("  ", .pad_string("", w_lin, "left"), " ", .pad_string(y_nome, largura_cols_y, "center"), "\n", sep = "")
  
  # Cabeçalho da tabela com o nome da variável de linha (eixo X)
  col_titulos <- c(
    .pad_string(x_nome, w_lin, "left"),
    sapply(categorias_y, function(cat_y) .pad_string(cat_y, w_col, "center"), USE.NAMES = FALSE),
    .pad_string("Total", w_col, "center")
  )
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  fmt_celula <- function(n_val, tot_ref) {
    if (proporcao == "nenhuma") {
      as.character(n_val)
    } else {
      pct <- if (tot_ref > 0) n_val / tot_ref else 0
      sprintf("%d (%s)", n_val, .fmt_pct(pct, decimais))
    }
  }
  
  for (i in seq_along(categorias_x)) {
    valores_linha <- sapply(seq_along(categorias_y), function(j) {
      n_ij <- tab[i, j]
      tot_ref <- switch(proporcao,
                        "linha" = tot_linha[i],
                        "coluna" = tot_col[j],
                        "total" = tot_geral,
                        1)
      .pad_string(fmt_celula(n_ij, tot_ref), w_col, "center")
    }, USE.NAMES = FALSE)
    
    tot_l_str <- if (proporcao == "linha") {
      sprintf("%d (100,0%%)", tot_linha[i])
    } else if (proporcao == "coluna") {
      as.character(tot_linha[i])
    } else if (proporcao == "total") {
      sprintf("%d (%s)", tot_linha[i], .fmt_pct(tot_linha[i] / tot_geral, decimais))
    } else {
      as.character(tot_linha[i])
    }
    
    linha_txt <- c(
      .pad_string(categorias_x[i], w_lin, "left"),
      valores_linha,
      .pad_string(tot_l_str, w_col, "center")
    )
    cat("  ", paste(linha_txt, collapse = " "), "\n", sep = "")
  }
  
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  valores_tot_col <- sapply(seq_along(categorias_y), function(j) {
    n_j <- tot_col[j]
    str_val <- if (proporcao == "coluna") {
      sprintf("%d (100,0%%)", n_j)
    } else if (proporcao == "total") {
      sprintf("%d (%s)", n_j, .fmt_pct(n_j / tot_geral, decimais))
    } else {
      as.character(n_j)
    }
    .pad_string(str_val, w_col, "center")
  }, USE.NAMES = FALSE)
  
  str_tot_geral <- if (proporcao %in% c("total", "linha", "coluna")) {
    sprintf("%d (100,0%%)", tot_geral)
  } else {
    as.character(tot_geral)
  }
  
  linha_tot <- c(
    .pad_string("Total", w_lin, "left"),
    valores_tot_col,
    .pad_string(str_tot_geral, w_col, "center")
  )
  cat("  ", paste(linha_tot, collapse = " "), "\n", sep = "")
  
  .print_footer()
  invisible(tab)
}


#' Medidas de Tendência Central
#'
#' @param x Vetor numérico.
#' @param decimais Casas decimais (padrão 2).
#' @export
med_tend_central <- function(x, decimais = 2) {
  if (!is.numeric(x)) stop("O argumento 'x' deve ser numérico.")
  
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)
  
  x_c <- x[!is.na(x)]
  n_obs <- length(x_c)
  
  media <- mean(x_c)
  mediana <- median(x_c)
  
  e_discreta <- is.integer(x_c) || all(x_c == round(x_c))
  
  if (e_discreta) {
    tab <- table(x_c)
    max_freq <- max(tab)
    pct_moda <- if (n_obs > 0) max_freq / n_obs else 0
    modas <- names(tab)[tab == max_freq]
    
    if (max_freq == 1 && length(tab) > 1) {
      str_moda <- "Amodal"
    } else if (length(modas) > 3) {
      str_moda <- sprintf("Multimodal (%s)", .fmt_pct(pct_moda, 1))
    } else {
      str_moda <- sprintf("%s (%s)", paste(modas, collapse = ", "), .fmt_pct(pct_moda, 1))
    }
    rotulo_modal <- "Moda"
  } else {
    h <- hist(x_c, plot = FALSE)
    idx_max <- which.max(h$counts)
    max_count <- h$counts[idx_max]
    pct_classe <- if (n_obs > 0) max_count / n_obs else 0
    lim_inf <- .fmt_num(h$breaks[idx_max], decimais)
    lim_sup <- .fmt_num(h$breaks[idx_max + 1], decimais)
    str_moda <- sprintf("[%s; %s) (%s)", lim_inf, lim_sup, .fmt_pct(pct_classe, 1))
    rotulo_modal <- "Classe Modal"
  }
  
  tit <- if (!is.null(var_nome) && nzchar(var_nome) && !identical(var_nome, "x")) {
    sprintf("MEDIDAS DE TENDÊNCIA CENTRAL (%s)", var_nome)
  } else {
    "MEDIDAS DE TENDÊNCIA CENTRAL"
  }
  .print_header(tit)
  
  w <- c(media = 10, mediana = 10, modal = max(18, nchar(str_moda) + 2))
  col_titulos <- c(
    .pad_string("Média", w["media"], "center"),
    .pad_string("Mediana", w["mediana"], "center"),
    .pad_string(rotulo_modal, w["modal"], "center")
  )
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  col_valores <- c(
    .pad_string(.fmt_num(media, decimais), w["media"], "center"),
    .pad_string(.fmt_num(mediana, decimais), w["mediana"], "center"),
    .pad_string(str_moda, w["modal"], "center")
  )
  cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
  
  .print_footer()
  invisible(list(media = media, mediana = mediana, moda_ou_classe = str_moda))
}


#' Medidas de Dispersão e Escala
#'
#' @param x Vetor numérico.
#' @param decimais Casas decimais (padrão 2).
#' @export
med_dispersao <- function(x, decimais = 2) {
  if (!is.numeric(x)) stop("O argumento 'x' deve ser numérico.")
  
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)
  
  x_c <- x[!is.na(x)]
  
  dp <- sd(x_c)
  var_val <- var(x_c)
  amp <- max(x_c) - min(x_c)
  cv <- (dp / mean(x_c)) * 100
  
  tit <- if (!is.null(var_nome) && nzchar(var_nome) && !identical(var_nome, "x")) {
    sprintf("MEDIDAS DE DISPERSÃO E ESCALA (%s)", var_nome)
  } else {
    "MEDIDAS DE DISPERSÃO E ESCALA"
  }
  .print_header(tit)
  
  w <- c(dp = 8, var = 12, amp = 11, cv = 10)
  col_titulos <- c(
    .pad_string("DP", w["dp"], "center"),
    .pad_string("Variância", w["var"], "center"),
    .pad_string("Amplitude", w["amp"], "center"),
    .pad_string("CV (%)", w["cv"], "center")
  )
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  col_valores <- c(
    .pad_string(.fmt_num(dp, decimais), w["dp"], "center"),
    .pad_string(.fmt_num(var_val, decimais), w["var"], "center"),
    .pad_string(.fmt_num(amp, decimais), w["amp"], "center"),
    .pad_string(paste0(.fmt_num(cv, decimais), "%"), w["cv"], "center")
  )
  cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
  
  .print_footer()
  invisible(list(desvio_padrao = dp, variancia = var_val, amplitude = amp, cv = cv))
}


#' Medidas de Forma da Distribuição
#'
#' @param x Vetor numérico.
#' @param decimais Casas decimais (padrão 2).
#' @export
med_forma <- function(x, decimais = 2) {
  if (!is.numeric(x)) stop("O argumento 'x' deve ser numérico.")
  
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)
  
  x_c <- x[!is.na(x)]
  
  assim <- .calcular_assimetria(x_c)
  curt <- .calcular_curtose(x_c)
  
  interp_assim <- if (is.na(assim)) {
    "-"
  } else if (assim > 0.5) {
    "Assimetria positiva / Cauda à direita"
  } else if (assim < -0.5) {
    "Assimetria negativa / Cauda à esquerda"
  } else {
    "Aproximadamente simétrica"
  }
  
  interp_curt <- if (is.na(curt)) {
    "-"
  } else if (curt > 0.5) {
    "Leptocúrtica / Mais pontiaguda que a normal"
  } else if (curt < -0.5) {
    "Platicúrtica / Mais achatada que a normal"
  } else {
    "Mesocúrtica / Similar à distribuição normal"
  }
  
  tit <- if (!is.null(var_nome) && nzchar(var_nome) && !identical(var_nome, "x")) {
    sprintf("MEDIDAS DE FORMA (%s)", var_nome)
  } else {
    "MEDIDAS DE FORMA"
  }
  .print_header(tit)
  
  cat(sprintf("  Assimetria: %s  (%s)\n", .fmt_num(assim, decimais), interp_assim))
  cat(sprintf("  Curtose:    %s  (%s)\n", .fmt_num(curt, decimais), interp_curt))
  
  .print_footer()
  invisible(list(assimetria = assim, curtose = curt))
}


#' Quartis da Distribuição
#'
#' @param x Vetor numérico.
#' @param decimais Casas decimais (padrão 2).
#' @export
quartis <- function(x, decimais = 2) {
  if (!is.numeric(x)) stop("O argumento 'x' deve ser numérico.")
  
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)
  
  x_c <- x[!is.na(x)]
  q_vals <- quantile(x_c, probs = c(0, 0.25, 0.5, 0.75, 1))
  
  tit <- if (!is.null(var_nome) && nzchar(var_nome) && !identical(var_nome, "x")) {
    sprintf("QUARTIS (%s)", var_nome)
  } else {
    "QUARTIS"
  }
  .print_header(tit)
  
  w <- c(min = 8, q1 = 8, q2 = 8, q3 = 8, max = 8)
  col_titulos <- c(
    .pad_string("Mín", w["min"], "center"),
    .pad_string("Q1", w["q1"], "center"),
    .pad_string("Q2", w["q2"], "center"),
    .pad_string("Q3", w["q3"], "center"),
    .pad_string("Máx", w["max"], "center")
  )
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  col_valores <- c(
    .pad_string(.fmt_num(q_vals[1], decimais), w["min"], "center"),
    .pad_string(.fmt_num(q_vals[2], decimais), w["q1"], "center"),
    .pad_string(.fmt_num(q_vals[3], decimais), w["q2"], "center"),
    .pad_string(.fmt_num(q_vals[4], decimais), w["q3"], "center"),
    .pad_string(.fmt_num(q_vals[5], decimais), w["max"], "center")
  )
  cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
  
  .print_footer()
  invisible(q_vals)
}


#' Quintis da Distribuição
#'
#' @param x Vetor numérico.
#' @param decimais Casas decimais (padrão 2).
#' @export
quintis <- function(x, decimais = 2) {
  if (!is.numeric(x)) stop("O argumento 'x' deve ser numérico.")
  
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)
  
  x_c <- x[!is.na(x)]
  probs <- seq(0.2, 0.8, by = 0.2)
  q_vals <- quantile(x_c, probs = probs)
  
  tit <- if (!is.null(var_nome) && nzchar(var_nome) && !identical(var_nome, "x")) {
    sprintf("QUINTIS (%s)", var_nome)
  } else {
    "QUINTIS"
  }
  .print_header(tit)
  
  w_col <- 8
  nomes_col <- c("20%", "40%", "60%", "80%")
  col_titulos <- sapply(nomes_col, function(nm) .pad_string(nm, w_col, "center"), USE.NAMES = FALSE)
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  col_valores <- sapply(q_vals, function(v) .pad_string(.fmt_num(v, decimais), w_col, "center"), USE.NAMES = FALSE)
  cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
  
  .print_footer()
  invisible(q_vals)
}


#' Decis da Distribuição
#'
#' @param x Vetor numérico.
#' @param decimais Casas decimais (padrão 2).
#' @export
decis <- function(x, decimais = 2) {
  if (!is.numeric(x)) stop("O argumento 'x' deve ser numérico.")
  
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)
  
  x_c <- x[!is.na(x)]
  probs <- seq(0.1, 0.9, by = 0.1)
  d_vals <- quantile(x_c, probs = probs)
  
  tit <- if (!is.null(var_nome) && nzchar(var_nome) && !identical(var_nome, "x")) {
    sprintf("DECIS (%s)", var_nome)
  } else {
    "DECIS"
  }
  .print_header(tit)
  
  w_col <- 7
  nomes_col <- paste0("D", 1:9)
  col_titulos <- sapply(nomes_col, function(nm) .pad_string(nm, w_col, "center"), USE.NAMES = FALSE)
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("─", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  col_valores <- sapply(d_vals, function(v) .pad_string(.fmt_num(v, decimais), w_col, "center"), USE.NAMES = FALSE)
  cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
  
  .print_footer()
  invisible(d_vals)
}


#' Análise de Outliers pelo Método IQR
#'
#' Detecta valores discrepantes (outliers) em um vetor numérico ou em todas as
#' variáveis numéricas de um data.frame, usando o critério do Intervalo
#' Interquartil (IQR): valores abaixo de Q1 - 1,5×IQR ou acima de Q3 + 1,5×IQR.
#'
#' @param x Um vetor numérico ou data.frame.
#' @param decimais Casas decimais para exibição dos limites e valores (padrão 2).
#' @return Retorna invisivelmente um data.frame com as colunas: variavel,
#'   n_outliers, pct_outliers, valores_abaixo, valores_acima.
#' @export
outliers <- function(x, decimais = 2) {
  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)

  .outliers_vetor <- function(vec, nome) {
    vec_c <- vec[!is.na(vec)]
    n     <- length(vec)
    q1    <- quantile(vec_c, 0.25)
    q3    <- quantile(vec_c, 0.75)
    iqr   <- q3 - q1
    li    <- q1 - 1.5 * iqr
    ls    <- q3 + 1.5 * iqr

    abaixo <- sort(vec_c[vec_c < li])
    acima  <- sort(vec_c[vec_c > ls])
    n_out  <- length(abaixo) + length(acima)
    n_sem  <- n - n_out
    pct_out <- n_out / n
    pct_sem <- n_sem / n

    tit <- sprintf("ANÁLISE DE OUTLIERS (%s)", nome)
    .print_header(tit)
    cat(sprintf("  Limite inferior (LI): %s  |  Limite superior (LS): %s\n\n",
                .fmt_num(li, decimais), .fmt_num(ls, decimais)))

    if (n_out == 0) {
      cat("  Nenhum outlier detectado pelo método IQR.\n\n")
    } else {
      w  <- c(res = 24, n = 12, pct = 14)
      hdr <- paste(
        .pad_string("Resultado",   w["res"], "left"),
        .pad_string("Contagem",    w["n"],   "center"),
        .pad_string("Proporção",   w["pct"], "center")
      )
      sep <- paste(rep("─", nchar(hdr)), collapse = "")
      cat("  ", hdr, "\n", sep = "")
      cat("  ", sep, "\n", sep = "")
      cat("  ", paste(.pad_string("Sem outliers",         w["res"], "left"),
                      .pad_string(as.character(n_sem),    w["n"],   "center"),
                      .pad_string(.fmt_pct(pct_sem, 1),   w["pct"], "center")), "\n", sep = "")
      cat("  ", paste(.pad_string("Outliers detectados",  w["res"], "left"),
                      .pad_string(as.character(n_out),    w["n"],   "center"),
                      .pad_string(.fmt_pct(pct_out, 1),   w["pct"], "center")), "\n", sep = "")
      cat("  ", sep, "\n", sep = "")
      cat("  ", paste(.pad_string("Total",                w["res"], "left"),
                      .pad_string(as.character(n),        w["n"],   "center"),
                      .pad_string("100,0%",               w["pct"], "center")), "\n", sep = "")

      cat("\n  Valores identificados como outliers:\n")
      w_col <- 40
      hdr_v <- paste(.pad_string("Abaixo do LI", w_col, "left"),
                     .pad_string("Acima do LS",  w_col, "right"))
      sep_v <- paste(rep("─", nchar(hdr_v)), collapse = "")
      cat("  ", sep_v, "\n", sep = "")
      cat("  ", hdr_v, "\n", sep = "")
      cat("  ", sep_v, "\n", sep = "")
      n_rows <- max(length(abaixo), length(acima))
      for (i in seq_len(n_rows)) {
        val_abaixo <- if (i <= length(abaixo)) .fmt_num(abaixo[i], decimais) else ""
        val_acima  <- if (i <= length(acima))  .fmt_num(acima[i],  decimais) else ""
        cat("  ", paste(.pad_string(val_abaixo, w_col, "left"),
                        .pad_string(val_acima,  w_col, "right")), "\n", sep = "")
      }
      cat("\n")
    }

    .print_footer()
    invisible(data.frame(
      variavel       = nome,
      n_outliers     = n_out,
      pct_outliers   = round(pct_out * 100, 1),
      valores_abaixo = I(list(abaixo)),
      valores_acima  = I(list(acima))
    ))
  }

  if (is.data.frame(x)) {
    nums <- names(x)[sapply(x, is.numeric)]
    if (length(nums) == 0) stop("Nenhuma variável numérica encontrada no data.frame.")
    resultado <- lapply(nums, function(nm) .outliers_vetor(x[[nm]], nm))
    invisible(do.call(rbind, resultado))
  } else if (is.numeric(x)) {
    .outliers_vetor(x, var_nome)
  } else {
    stop("O argumento 'x' deve ser um vetor numérico ou um data.frame.")
  }
}


#' Percentis da Distribuição
#'
#' Calcula e exibe os percentis de um vetor numérico. Por padrão exibe
#' P5, P10, P25, P50, P75, P90 e P95, mas o usuário pode definir quais
#' percentis quer calcular pelo parâmetro \code{p}.
#'
#' @param x Vetor numérico.
#' @param p Vetor numérico com os percentis desejados (entre 0 e 100). Padrão: c(5, 10, 25, 50, 75, 90, 95).
#' @param decimais Casas decimais (padrão 2).
#' @return Retorna invisivelmente um vetor nomeado com os valores dos percentis.
#' @export
percentis <- function(x, p = c(5, 10, 25, 50, 75, 90, 95), decimais = 2) {
  if (!is.numeric(x)) stop("O argumento 'x' deve ser numérico.")
  if (any(p < 0 | p > 100)) stop("Os valores de 'p' devem estar entre 0 e 100.")

  var_expr <- deparse(substitute(x))
  var_nome <- sub(".*\\$", "", var_expr)

  x_c    <- x[!is.na(x)]
  probs  <- p / 100
  p_vals <- quantile(x_c, probs = probs)
  nomes  <- paste0("P", p)

  tit <- if (!is.null(var_nome) && nzchar(var_nome) && !identical(var_nome, "x")) {
    sprintf("PERCENTIS (%s)", var_nome)
  } else {
    "PERCENTIS"
  }
  .print_header(tit)

  w_label <- max(nchar("Percentil"), max(nchar(nomes))) + 2
  w_val   <- max(nchar("Valor"), 8) + 2

  hdr <- paste(.pad_string("Percentil", w_label, "left"),
               .pad_string("Valor",     w_val,   "right"))
  sep <- paste(rep("─", nchar(hdr)), collapse = "")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", sep, "\n", sep = "")

  for (i in seq_along(p_vals)) {
    linha <- paste(.pad_string(nomes[i],                   w_label, "left"),
                   .pad_string(.fmt_num(p_vals[i], decimais), w_val, "right"))
    cat("  ", linha, "\n", sep = "")
  }
  cat("  ", sep, "\n", sep = "")

  .print_footer()
  names(p_vals) <- nomes
  invisible(p_vals)
}
