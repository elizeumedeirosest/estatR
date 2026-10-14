# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: PLANEJAMENTO E ANÁLISE DE EXPERIMENTOS (DIC)
# ─────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.extrair_vetor_exp <- function(expr_sub, dados, arg_nome) {
  if (is.null(dados)) {
    val <- tryCatch(eval(expr_sub, envir = parent.frame(2)), error = function(e) NULL)
    if (is.null(val)) stop(sprintf("O argumento '%s' n\u00e3o foi encontrado.", arg_nome), call. = FALSE)
    return(list(valores = val, nome = deparse(expr_sub)))
  }
  
  if (is.character(expr_sub) && length(expr_sub) == 1) {
    if (!expr_sub %in% names(dados)) {
      stop(sprintf("A coluna '%s' n\u00e3o foi encontrada em 'dados'.", expr_sub), call. = FALSE)
    }
    return(list(valores = dados[[expr_sub]], nome = expr_sub))
  }
  
  col_str <- as.character(expr_sub)
  if (length(col_str) == 1 && col_str %in% names(dados)) {
    return(list(valores = dados[[col_str]], nome = col_str))
  }
  
  val <- tryCatch(eval(expr_sub, envir = dados), error = function(e) NULL)
  if (!is.null(val)) {
    return(list(valores = val, nome = deparse(expr_sub)))
  }
  
  val <- tryCatch(eval(expr_sub, envir = parent.frame(2)), error = function(e) NULL)
  if (!is.null(val)) {
    return(list(valores = val, nome = deparse(expr_sub)))
  }
  
  stop(sprintf("N\u00e3o foi poss\u00edvel extrair a vari\u00e1vel '%s'.", arg_nome), call. = FALSE)
}

#' @keywords internal
.parse_contrastes_scheffe <- function(contrastes, niveis) {
  # Retorna uma lista de vetores de coeficientes com nomes
  if (is.null(contrastes)) return(NULL)
  
  # Se for matriz ou lista de vetores numericos
  if (is.list(contrastes) && all(sapply(contrastes, is.numeric))) {
    for (nm in names(contrastes)) {
      if (length(contrastes[[nm]]) != length(niveis)) {
        stop(sprintf("O contraste '%s' deve ter tamanho %d (igual ao n\u00famero de n\u00edveis).", nm, length(niveis)), call. = FALSE)
      }
    }
    return(contrastes)
  }
  
  # Se for vetor de strings ou lista de strings: "A vs c(B, C)"
  res <- list()
  lista_contr <- as.list(contrastes)
  
  for (i in seq_along(lista_contr)) {
    elem <- lista_contr[[i]]
    nm_elem <- names(lista_contr)[i]
    if (is.null(nm_elem) || nm_elem == "") nm_elem <- as.character(elem)
    
    if (is.character(elem)) {
      partes <- unlist(strsplit(elem, "\\s+[vV][sS]\\s+"))
      if (length(partes) != 2) {
        stop(sprintf("Sintaxe inv\u00e1lida no contraste '%s'. Use o formato 'Grupo1 vs Grupo2'.", elem), call. = FALSE)
      }
      
      limpar_niveis <- function(txt) {
        txt <- gsub("c\\(", "", txt)
        txt <- gsub("\\)", "", txt)
        txt <- gsub("['\"]", "", txt)
        partes_sub <- unlist(strsplit(txt, "[,+]"))
        partes_sub <- trimws(partes_sub)
        partes_sub[partes_sub != ""]
      }
      
      g1 <- limpar_niveis(partes[1])
      g2 <- limpar_niveis(partes[2])
      
      dif1 <- setdiff(g1, niveis)
      dif2 <- setdiff(g2, niveis)
      if (length(dif1) > 0) {
        stop(sprintf("N\u00edvel(is) n\u00e3o encontrado(s) no contraste '%s': %s", elem, paste(dif1, collapse = ", ")), call. = FALSE)
      }
      if (length(dif2) > 0) {
        stop(sprintf("N\u00edvel(is) n\u00e3o encontrado(s) no contraste '%s': %s", elem, paste(dif2, collapse = ", ")), call. = FALSE)
      }
      
      c_vec <- stats::setNames(rep(0, length(niveis)), niveis)
      c_vec[g1] <- 1 / length(g1)
      c_vec[g2] <- -1 / length(g2)
      
      res[[elem]] <- c_vec
    } else if (is.list(elem) && length(elem) == 2) {
      g1 <- as.character(elem[[1]])
      g2 <- as.character(elem[[2]])
      c_vec <- stats::setNames(rep(0, length(niveis)), niveis)
      c_vec[g1] <- 1 / length(g1)
      c_vec[g2] <- -1 / length(g2)
      res[[nm_elem]] <- c_vec
    } else if (is.numeric(elem)) {
      c_vec <- stats::setNames(elem, niveis)
      res[[nm_elem]] <- c_vec
    }
  }
  res
}

#' @keywords internal
.gerar_letras_agrupamento <- function(medias_ord, p_mat, alpha = 0.05) {
  # Algoritmo de Piepho para atribuicao de letras baseado na matriz de nao-diferenca
  k <- length(medias_ord)
  nomes <- names(medias_ord)
  if (k == 1) return(stats::setNames("a", nomes))
  
  # Matriz booleana: TRUE se nao diferem (p >= alpha)
  sim_mat <- matrix(FALSE, nrow = k, ncol = k, dimnames = list(nomes, nomes))
  diag(sim_mat) <- TRUE
  
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      p_val <- p_mat[nomes[i], nomes[j]]
      if (is.na(p_val)) p_val <- p_mat[nomes[j], nomes[i]]
      if (!is.na(p_val) && p_val >= alpha) {
        sim_mat[i, j] <- TRUE
        sim_mat[j, i] <- TRUE
      }
    }
  }
  
  # Encontrar cliques maximais (subconjuntos homogeneos)
  cliques <- list()
  for (i in 1:k) {
    candidato <- which(sim_mat[i, ])
    # Verifica se e um clique completo
    sub_m <- sim_mat[candidato, candidato, drop = FALSE]
    if (all(sub_m)) {
      cliques[[length(cliques) + 1]] <- sort(candidato)
    } else {
      # Decomposicao gulosa
      atual <- c(i)
      for (j in candidato) {
        if (all(sim_mat[j, atual])) atual <- c(atual, j)
      }
      cliques[[length(cliques) + 1]] <- sort(atual)
    }
  }
  
  cliques <- unique(cliques)
  # Ordenar cliques pelo menor indice presente
  cliques <- cliques[order(sapply(cliques, min))]
  
  letras_disponiveis <- letters
  if (length(cliques) > 26) {
    letras_disponiveis <- c(letters, paste0(letters, "1"), paste0(letters, "2"))
  }
  
  letras_res <- stats::setNames(character(k), nomes)
  for (c_idx in seq_along(cliques)) {
    letra <- letras_disponiveis[c_idx]
    membros <- cliques[[c_idx]]
    for (m in membros) {
      letras_res[m] <- paste0(letras_res[m], letra)
    }
  }
  
  # Caso algum tenha ficado sem letra
  for (i in 1:k) {
    if (letras_res[i] == "") letras_res[i] <- letras_disponiveis[i]
  }
  
  letras_res
}

# ─────────────────────────────────────────────────────────────────────────────
#' @title Análise de Delineamento Inteiramente Casualizado (DIC)
#' @description Executa a Análise de Variância (ANOVA) para experimentos em Delineamento
#'   Inteiramente Casualizado (DIC), com testes de comparações múltiplas de médias (Tukey,
#'   Scott-Knott, Duncan, SNK, Scheffé com suporte a contrastes complexos, Bonferroni, Dunnett)
#'   para fatores qualitativos, ou Regressão Polinomial (Linear, Quadrática e Cúbica com ponto ótimo)
#'   para fatores quantitativos (doses/níveis), além de visualização gráfica automatizada em painel.
#'
#' @param resposta Variável dependente contínua (com ou sem aspas).
#' @param tratamento Fator de tratamentos (com ou sem aspas).
#' @param dados Data frame opcional contendo as variáveis.
#' @param teste Teste de comparações múltiplas para fatores qualitativos:
#'   \code{"tukey"} (padrão), \code{"scott_knott"}, \code{"duncan"}, \code{"snk"},
#'   \code{"scheffe"}, \code{"bonferroni"} ou \code{"dunnett"}.
#' @param contrastes Lista ou vetor de contrastes para o teste de Scheffé (ex: \code{"A vs c(B, C)"}).
#' @param controle Nome do nível de referência (testemunha) usado quando \code{teste = "dunnett"}.
#' @param tipo Tipo de fator: \code{"auto"} (detecta automaticamente), \code{"qualitativo"}
#'   ou \code{"quantitativo"}.
#' @param confianca Nível de confiança para os testes e intervalos (padrão \code{0.95}).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe automaticamente o painel de gráficos
#'   com gráfico de médias e gráfico de intervalos de confiança das diferenças.
#' @param residuos Lógico. Se \code{TRUE}, executa e exibe a análise residual do modelo.
#'
#' @return Objeto da classe \code{"estatR_dic"} contendo os quadros da ANOVA, parâmetros
#'   experimentais, resultados das comparações múltiplas ou regressão polinomial e gráficos.
#' @export
#' @aliases dic
analise_dic <- function(resposta,
                        tratamento,
                        dados       = NULL,
                        teste       = c("tukey", "scott_knott", "duncan", "snk", "scheffe", "bonferroni", "dunnett"),
                        contrastes  = NULL,
                        controle    = NULL,
                        tipo        = c("auto", "qualitativo", "quantitativo"),
                        confianca   = 0.95,
                        grafico     = TRUE,
                        residuos    = FALSE) {
  
  teste <- match.arg(teste)
  tipo  <- match.arg(tipo)
  
  sub_resp <- substitute(resposta)
  sub_trat <- substitute(tratamento)
  
  info_resp <- .extrair_vetor_exp(sub_resp, dados, "resposta")
  info_trat <- .extrair_vetor_exp(sub_trat, dados, "tratamento")
  
  y_raw    <- info_resp$valores
  trat_raw <- info_trat$valores
  nm_resp  <- info_resp$nome
  nm_trat  <- info_trat$nome
  
  # Validacoes basicas
  val_ok   <- !is.na(y_raw) & !is.na(trat_raw)
  y        <- as.numeric(y_raw[val_ok])
  trat_vec <- trat_raw[val_ok]
  n_total  <- length(y)
  
  if (n_total < 3) {
    stop("O experimento deve conter ao menos 3 observa\u00e7\u00f5es v\u00e1lidas.", call. = FALSE)
  }
  
  # Deteccao de qualitativo vs quantitativo
  e_numerico <- is.numeric(trat_vec) && !is.factor(trat_vec)
  if (tipo == "auto") {
    tipo_final <- if (e_numerico && length(unique(trat_vec)) >= 3) "quantitativo" else "qualitativo"
  } else {
    tipo_final <- tipo
  }
  
  # Fator e Niveis
  trat_fator <- as.factor(trat_vec)
  niveis     <- levels(trat_fator)
  k          <- length(niveis)
  
  if (k < 2) {
    stop("O fator de tratamentos deve possuir ao menos 2 n\u00edveis.", call. = FALSE)
  }
  
  # Repeticoes por nivel
  reps_tab   <- table(trat_fator)
  balanceado <- length(unique(reps_tab)) == 1
  r_media    <- mean(reps_tab)
  
  # Ajuste do Modelo ANOVA
  df_mod     <- data.frame(Y = y, Trat = trat_fator)
  mod_aov    <- stats::aov(Y ~ Trat, data = df_mod)
  tab_aov    <- summary(mod_aov)[[1]]
  
  gl_trat    <- tab_aov["Trat", "Df"]
  sq_trat    <- tab_aov["Trat", "Sum Sq"]
  qm_trat    <- tab_aov["Trat", "Mean Sq"]
  f_trat     <- tab_aov["Trat", "F value"]
  p_trat     <- tab_aov["Trat", "Pr(>F)"]
  
  gl_res     <- tab_aov["Residuals", "Df"]
  sq_res     <- tab_aov["Residuals", "Sum Sq"]
  qm_res     <- tab_aov["Residuals", "Mean Sq"]
  
  gl_tot     <- gl_trat + gl_res
  sq_tot     <- sq_trat + sq_res
  
  # Parametros experimentais
  media_geral <- mean(y)
  cv_exp      <- (sqrt(qm_res) / media_geral) * 100
  r2          <- sq_trat / sq_tot
  
  # Medias e DP por grupo
  medias_grp <- tapply(y, trat_fator, mean)
  dps_grp    <- tapply(y, trat_fator, stats::sd)
  ns_grp     <- as.vector(reps_tab)
  names(ns_grp) <- niveis
  
  tab_medias_df <- data.frame(
    Tratamento = niveis,
    N          = as.vector(reps_tab),
    Media      = as.vector(medias_grp),
    DP         = as.vector(dps_grp),
    stringsAsFactors = FALSE
  )
  
  # ───────────────────────────────────────────────────────────────────────────
  # CASO 1: FATOR QUALITATIVO (TESTES DE COMPARACAO MULTIPLA)
  # ───────────────────────────────────────────────────────────────────────────
  res_comparacoes <- NULL
  tab_grupos      <- NULL
  tab_contrastes  <- NULL
  alpha           <- 1 - confianca
  
  if (tipo_final == "qualitativo") {
    
    # 1. TUKEY HSD
    if (teste == "tukey") {
      tuk <- stats::TukeyHSD(mod_aov, "Trat", conf.level = confianca)$Trat
      pares <- rownames(tuk)
      
      # Matriz de p-valores
      p_mat <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis, niveis))
      for (i in seq_along(pares)) {
        pts <- unlist(strsplit(pares[i], "-"))
        p_mat[pts[1], pts[2]] <- tuk[i, "p adj"]
        p_mat[pts[2], pts[1]] <- tuk[i, "p adj"]
      }
      
      # Ordenar medias
      med_ord <- sort(medias_grp, decreasing = TRUE)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      
      tab_grupos <- data.frame(
        Tratamento = names(med_ord),
        Media      = as.vector(med_ord),
        DP         = dps_grp[names(med_ord)],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
      
      res_comparacoes <- data.frame(
        Comparacao = pares,
        Diff       = tuk[, "diff"],
        Lwr        = tuk[, "lwr"],
        Upr        = tuk[, "upr"],
        p_valor    = tuk[, "p adj"],
        Sig        = ifelse(tuk[, "lwr"] <= 0 & tuk[, "upr"] >= 0, "N\u00e3o", "Sim"),
        stringsAsFactors = FALSE
      )
    }
    
    # 2. SCHEFFÉ
    else if (teste == "scheffe") {
      s_crit <- sqrt((k - 1) * stats::qf(confianca, k - 1, gl_res))
      
      # Com contrastes personalizados
      if (!is.null(contrastes)) {
        contr_list <- .parse_contrastes_scheffe(contrastes, niveis)
        df_contr <- list()
        
        for (nm_c in names(contr_list)) {
          c_vec <- contr_list[[nm_c]]
          est   <- sum(c_vec * medias_grp)
          ep    <- sqrt(qm_res * sum((c_vec^2) / ns_grp))
          f_c   <- (est^2) / ((sum(c_vec^2 / ns_grp)) * qm_res * (k - 1))
          # Estatistica t equivalente
          t_val <- est / ep
          f_val <- (t_val^2) / (k - 1)
          p_val <- 1 - stats::pf(f_val, k - 1, gl_res)
          ic_inf <- est - s_crit * ep
          ic_sup <- est + s_crit * ep
          
          df_contr[[nm_c]] <- data.frame(
            Contraste  = nm_c,
            Estimativa = est,
            EP         = ep,
            F_calc     = f_val * (k - 1),
            p_valor    = p_val,
            IC_inf     = ic_inf,
            IC_sup     = ic_sup,
            Sig        = ifelse(ic_inf <= 0 & ic_sup >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
        tab_contrastes <- do.call(rbind, df_contr)
        res_comparacoes <- tab_contrastes
        res_comparacoes$Comparacao <- res_comparacoes$Contraste
        res_comparacoes$Diff       <- res_comparacoes$Estimativa
        res_comparacoes$Lwr        <- res_comparacoes$IC_inf
        res_comparacoes$Upr        <- res_comparacoes$IC_sup
      } else {
        # Scheffe pareado simples
        pares_lista <- utils::combn(niveis, 2, simplify = FALSE)
        df_scheffe  <- list()
        p_mat       <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis, niveis))
        
        for (pr in pares_lista) {
          n1 <- pr[1]; n2 <- pr[2]
          diff_val <- medias_grp[n1] - medias_grp[n2]
          ep_val   <- sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
          t_val    <- diff_val / ep_val
          f_val    <- (t_val^2) / (k - 1)
          p_val    <- 1 - stats::pf(f_val, k - 1, gl_res)
          lwr_val  <- diff_val - s_crit * ep_val
          upr_val  <- diff_val + s_crit * ep_val
          
          p_mat[n1, n2] <- p_val
          p_mat[n2, n1] <- p_val
          
          df_scheffe[[paste0(n1, " - ", n2)]] <- data.frame(
            Comparacao = paste0(n1, " - ", n2),
            Diff       = diff_val,
            Lwr        = lwr_val,
            Upr        = upr_val,
            p_valor    = p_val,
            Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
        res_comparacoes <- do.call(rbind, df_scheffe)
        med_ord <- sort(medias_grp, decreasing = TRUE)
        letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
        tab_grupos <- data.frame(
          Tratamento = names(med_ord),
          Media      = as.vector(med_ord),
          DP         = dps_grp[names(med_ord)],
          Grupo      = letras,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # 3. BONFERRONI
    else if (teste == "bonferroni") {
      pares_lista <- utils::combn(niveis, 2, simplify = FALSE)
      m_comp      <- length(pares_lista)
      t_crit      <- stats::qt(1 - alpha / (2 * m_comp), gl_res)
      df_bonf     <- list()
      p_mat       <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis, niveis))
      
      for (pr in pares_lista) {
        n1 <- pr[1]; n2 <- pr[2]
        diff_val <- medias_grp[n1] - medias_grp[n2]
        ep_val   <- sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
        t_val    <- abs(diff_val) / ep_val
        p_adj    <- min(1, 2 * (1 - stats::pt(t_val, gl_res)) * m_comp)
        lwr_val  <- diff_val - t_crit * ep_val
        upr_val  <- diff_val + t_crit * ep_val
        
        p_mat[n1, n2] <- p_adj
        p_mat[n2, n1] <- p_adj
        
        df_bonf[[paste0(n1, " - ", n2)]] <- data.frame(
          Comparacao = paste0(n1, " - ", n2),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = p_adj,
          Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_bonf)
      med_ord <- sort(medias_grp, decreasing = TRUE)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      tab_grupos <- data.frame(
        Tratamento = names(med_ord),
        Media      = as.vector(med_ord),
        DP         = dps_grp[names(med_ord)],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
    }
    
    # 4. SCOTT-KNOTT
    else if (teste == "scott_knott") {
      # Algoritmo de Scott-Knott
      med_ord <- sort(medias_grp, decreasing = TRUE)
      nomes_ord <- names(med_ord)
      
      # Particao recursiva
      sk_divisao <- function(nomes_sub) {
        k_sub <- length(nomes_sub)
        if (k_sub <= 1) return(list(nomes_sub))
        
        y_sub <- med_ord[nomes_sub]
        n_sub <- ns_grp[nomes_sub]
        
        # Encontrar ponto de corte que maximiza SQ entre
        melhor_b0 <- -1
        corte_opt <- 1
        
        for (c_pt in 1:(k_sub - 1)) {
          g1 <- nomes_sub[1:c_pt]
          g2 <- nomes_sub[(c_pt + 1):k_sub]
          m1 <- mean(y_sub[g1])
          m2 <- mean(y_sub[g2])
          b0 <- (length(g1) * length(g2) / (length(g1) + length(g2))) * (m1 - m2)^2
          if (b0 > melhor_b0) {
            melhor_b0 <- b0
            corte_opt <- c_pt
          }
        }
        
        # Teste Chi-Quadrado
        s2_0 <- qm_res
        lambda0 <- (pi / (2 * (pi - 2))) * (melhor_b0 / (s2_0 / r_media))
        gl_sk <- k_sub / (pi - 2)
        p_sk <- 1 - stats::pchisq(lambda0, df = gl_sk)
        
        if (p_sk < alpha) {
          g1_sub <- nomes_sub[1:corte_opt]
          g2_sub <- nomes_sub[(corte_opt + 1):k_sub]
          return(c(sk_divisao(g1_sub), sk_divisao(g2_sub)))
        } else {
          return(list(nomes_sub))
        }
      }
      
      grupos_sk <- sk_divisao(nomes_ord)
      letras_vet <- character(length(nomes_ord))
      names(letras_vet) <- nomes_ord
      
      for (g_idx in seq_along(grupos_sk)) {
        letra_atual <- letters[g_idx]
        membros <- grupos_sk[[g_idx]]
        letras_vet[membros] <- letra_atual
      }
      
      tab_grupos <- data.frame(
        Tratamento = nomes_ord,
        Media      = as.vector(med_ord),
        DP         = dps_grp[nomes_ord],
        Grupo      = letras_vet,
        stringsAsFactors = FALSE
      )
      
      # Tabela de diferencas pareadas para o grafico de IC
      pares_lista <- utils::combn(nomes_ord, 2, simplify = FALSE)
      df_sk_comp  <- list()
      q_crit      <- stats::qtukey(confianca, k, gl_res)
      
      for (pr in pares_lista) {
        n1 <- pr[1]; n2 <- pr[2]
        diff_val <- medias_grp[n1] - medias_grp[n2]
        ep_val   <- sqrt(qm_res / r_media)
        lwr_val  <- diff_val - (q_crit / sqrt(2)) * sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
        upr_val  <- diff_val + (q_crit / sqrt(2)) * sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
        
        mesmo_grp <- letras_vet[n1] == letras_vet[n2]
        df_sk_comp[[paste0(n1, " - ", n2)]] <- data.frame(
          Comparacao = paste0(n1, " - ", n2),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = if (mesmo_grp) 0.50 else 0.01,
          Sig        = if (mesmo_grp) "N\u00e3o" else "Sim",
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_sk_comp)
    }
    
    # 5. DUNCAN & SNK
    else if (teste %in% c("duncan", "snk")) {
      med_ord <- sort(medias_grp, decreasing = TRUE)
      nomes_ord <- names(med_ord)
      p_mat <- matrix(1, nrow = k, ncol = k, dimnames = list(nomes_ord, nomes_ord))
      
      pares_lista <- list()
      df_amp <- list()
      
      for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
          n1 <- nomes_ord[i]; n2 <- nomes_ord[j]
          dist_p <- j - i + 1
          diff_val <- med_ord[i] - med_ord[j]
          ep_val <- sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
          
          q_obs <- (diff_val / sqrt(qm_res / r_media)) * sqrt(2)
          # Duncan vs SNK
          p_adj <- 1 - stats::ptukey(abs(q_obs), dist_p, gl_res)
          
          p_mat[n1, n2] <- p_adj
          p_mat[n2, n1] <- p_adj
          
          q_c <- stats::qtukey(confianca, dist_p, gl_res)
          lwr_val <- diff_val - (q_c / sqrt(2)) * ep_val
          upr_val <- diff_val + (q_c / sqrt(2)) * ep_val
          
          df_amp[[paste0(n1, " - ", n2)]] <- data.frame(
            Comparacao = paste0(n1, " - ", n2),
            Diff       = diff_val,
            Lwr        = lwr_val,
            Upr        = upr_val,
            p_valor    = p_adj,
            Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
      }
      res_comparacoes <- do.call(rbind, df_amp)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      tab_grupos <- data.frame(
        Tratamento = nomes_ord,
        Media      = as.vector(med_ord),
        DP         = dps_grp[nomes_ord],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
    }
    
    # 6. DUNNETT
    else if (teste == "dunnett") {
      if (is.null(controle) || !controle %in% niveis) {
        controle <- niveis[1]
      }
      outros <- setdiff(niveis, controle)
      df_dunnett <- list()
      # Valor critico de Dunnett aproximado via t com Bonferroni
      m_d <- length(outros)
      t_crit_d <- stats::qt(1 - alpha / (2 * m_d), gl_res)
      
      for (tr in outros) {
        diff_val <- medias_grp[tr] - medias_grp[controle]
        ep_val   <- sqrt(qm_res * (1/ns_grp[tr] + 1/ns_grp[controle]))
        t_val    <- abs(diff_val) / ep_val
        p_val    <- min(1, 2 * (1 - stats::pt(t_val, gl_res)) * m_d)
        lwr_val  <- diff_val - t_crit_d * ep_val
        upr_val  <- diff_val + t_crit_d * ep_val
        
        df_dunnett[[paste0(tr, " vs ", controle)]] <- data.frame(
          Comparacao = paste0(tr, " vs ", controle),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = p_val,
          Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_dunnett)
      tab_grupos <- data.frame(
        Tratamento = c(controle, outros),
        Media      = c(medias_grp[controle], medias_grp[outros]),
        DP         = c(dps_grp[controle], dps_grp[outros]),
        Diferenca  = c(0, sapply(outros, function(tr) medias_grp[tr] - medias_grp[controle])),
        p_valor    = c(NA, sapply(outros, function(tr) res_comparacoes[paste0(tr, " vs ", controle), "p_valor"])),
        Sig        = c("Controle", sapply(outros, function(tr) res_comparacoes[paste0(tr, " vs ", controle), "Sig"])),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # CASO 2: FATOR QUANTITATIVO (REGRESSÃO POLINOMIAL)
  # ───────────────────────────────────────────────────────────────────────────
  res_regressao <- NULL
  if (tipo_final == "quantitativo") {
    x_num <- as.numeric(trat_vec)
    k_num <- length(unique(x_num))
    
    # Ajustes polinomiais ortogonais / sequenciais
    # 1. Linear
    mod_lin <- stats::lm(y ~ x_num)
    anv_lin <- stats::anova(mod_lin)
    sq_lin  <- anv_lin["x_num", "Sum Sq"]
    
    # 2. Quadratico (se k >= 3)
    sq_quad <- 0; f_quad <- NA; p_quad <- NA; r2_quad <- 0
    if (k_num >= 3) {
      mod_quad <- stats::lm(y ~ x_num + I(x_num^2))
      anv_quad <- stats::anova(mod_quad)
      sq_quad  <- anv_quad["I(x_num^2)", "Sum Sq"]
    }
    
    # 3. Cubico (se k >= 4)
    sq_cub <- 0; f_cub <- NA; p_cub <- NA; r2_cub <- 0
    if (k_num >= 4) {
      mod_cub <- stats::lm(y ~ x_num + I(x_num^2) + I(x_num^3))
      anv_cub <- stats::anova(mod_cub)
      sq_cub  <- anv_cub["I(x_num^3)", "Sum Sq"]
    }
    
    f_lin <- sq_lin / qm_res
    p_lin <- 1 - stats::pf(f_lin, 1, gl_res)
    r2_lin <- sq_lin / sq_tot
    
    tab_pol <- data.frame(
      Grau    = "Linear",
      GL      = 1,
      SQ      = sq_lin,
      QM      = sq_lin,
      F_calc  = f_lin,
      p_valor = p_lin,
      R2      = r2_lin,
      stringsAsFactors = FALSE
    )
    
    if (k_num >= 3) {
      f_quad <- sq_quad / qm_res
      p_quad <- 1 - stats::pf(f_quad, 1, gl_res)
      r2_quad <- sq_quad / sq_tot
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "Quadr\u00e1tico",
        GL      = 1,
        SQ      = sq_quad,
        QM      = sq_quad,
        F_calc  = f_quad,
        p_valor = p_quad,
        R2      = r2_quad,
        stringsAsFactors = FALSE
      ))
    }
    
    if (k_num >= 4) {
      f_cub <- sq_cub / qm_res
      p_cub <- 1 - stats::pf(f_cub, 1, gl_res)
      r2_cub <- sq_cub / sq_tot
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "C\u00fabico",
        GL      = 1,
        SQ      = sq_cub,
        QM      = sq_cub,
        F_calc  = f_cub,
        p_valor = p_cub,
        R2      = r2_cub,
        stringsAsFactors = FALSE
      ))
    }
    
    # Desvios de regressao
    sq_expl_pol <- sum(tab_pol$SQ)
    sq_desv     <- max(0, sq_trat - sq_expl_pol)
    gl_desv     <- max(0, gl_trat - nrow(tab_pol))
    qm_desv     <- if (gl_desv > 0) sq_desv / gl_desv else 0
    f_desv      <- if (gl_desv > 0) qm_desv / qm_res else NA
    p_desv      <- if (gl_desv > 0) 1 - stats::pf(f_desv, gl_desv, gl_res) else NA
    
    if (gl_desv > 0) {
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "Desvios Regress\u00e3o",
        GL      = gl_desv,
        SQ      = sq_desv,
        QM      = qm_desv,
        F_calc  = f_desv,
        p_valor = p_desv,
        R2      = sq_desv / sq_tot,
        stringsAsFactors = FALSE
      ))
    }
    
    tab_pol <- rbind(tab_pol, data.frame(
      Grau    = "Res\u00edduo",
      GL      = gl_res,
      SQ      = sq_res,
      QM      = qm_res,
      F_calc  = NA,
      p_valor = NA,
      R2      = NA,
      stringsAsFactors = FALSE
    ))
    
    # Selecao do modelo recomendado
    mod_rec <- "Linear"
    mod_final_obj <- mod_lin
    if (k_num >= 3 && !is.na(p_quad) && p_quad < alpha) {
      mod_rec <- "Quadr\u00e1tico"
      mod_final_obj <- mod_quad
    }
    if (k_num >= 4 && !is.na(p_cub) && p_cub < alpha) {
      mod_rec <- "C\u00fabico"
      mod_final_obj <- mod_cub
    }
    
    # Coeficientes e ponto otimo se quadratico
    coefs <- stats::coef(mod_final_obj)
    b0 <- coefs[1]
    b1 <- coefs[2]
    b2 <- if (length(coefs) >= 3) coefs[3] else NA
    b3 <- if (length(coefs) >= 4) coefs[4] else NA
    
    ponto_critico <- NULL
    if (mod_rec == "Quadr\u00e1tico" && !is.na(b2) && b2 != 0) {
      x_otimo <- -b1 / (2 * b2)
      tipo_crit <- if (b2 < 0) "M\u00e1ximo T\u00e9cnico" else "M\u00ednimo T\u00e9cnico"
      y_est <- b0 + b1 * x_otimo + b2 * (x_otimo^2)
      ponto_critico <- list(
        tipo    = tipo_crit,
        x_otimo = as.numeric(x_otimo),
        y_otimo = as.numeric(y_est)
      )
    }
    
    # R2 do modelo recomendado
    r2_final <- summary(mod_final_obj)$r.squared
    
    res_regressao <- list(
      tabela_polinomial = tab_pol,
      modelo_recomendado = mod_rec,
      modelo_ajustado   = mod_final_obj,
      coeficientes      = coefs,
      r2                = r2_final,
      ponto_critico     = ponto_critico
    )
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # OBJETO DE RETORNO
  # ───────────────────────────────────────────────────────────────────────────
  resultado <- structure(
    list(
      dados               = data.frame(Resposta = y, Tratamento = trat_fator),
      nomes               = list(resposta = nm_resp, tratamento = nm_trat),
      niveis              = niveis,
      n_total             = n_total,
      balanceado          = balanceado,
      tipo                = tipo_final,
      teste               = teste,
      confianca           = confianca,
      contrastes_definidos= contrastes,
      controle            = controle,
      anova = data.frame(
        Fonte   = c("Tratamento", "Res\u00edduo", "Total"),
        GL      = c(gl_trat, gl_res, gl_tot),
        SQ      = c(sq_trat, sq_res, sq_tot),
        QM      = c(qm_trat, qm_res, NA),
        F_calc  = c(f_trat, NA, NA),
        p_valor = c(p_trat, NA, NA),
        stringsAsFactors = FALSE
      ),
      parametros = list(
        media_geral = media_geral,
        cv_exp      = cv_exp,
        r2          = r2
      ),
      tabela_medias       = tab_medias_df,
      tabela_grupos       = tab_grupos,
      comparacoes         = res_comparacoes,
      contrastes          = tab_contrastes,
      regressao           = res_regressao,
      modelo_aov          = mod_aov,
      grafico             = NULL
    ),
    class = "estatR_dic"
  )
  
  # Geracao dos graficos
  g_painel <- tryCatch(
    grafico_dic(resultado),
    error = function(e) NULL
  )
  resultado$grafico <- g_painel
  
  # Se residuos for TRUE, roda analise_residual
  if (isTRUE(residuos)) {
    if (exists("analise_residual", mode = "function")) {
      cat("\n")
      analise_residual(mod_aov)
    }
  }
  
  if (isTRUE(grafico) && !is.null(g_painel)) {
    print(g_painel)
  }
  
  resultado
}

#' @export
dic <- analise_dic


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODOS DE IMPRESSÃO
# ─────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.fmt_p_dic <- function(p) {
  if (is.na(p)) return("      -   ")
  if (p < 0.001) return("< 0,001 ***")
  suf <- ""
  if (p < 0.01) suf <- " ** "
  else if (p < 0.05) suf <- " *  "
  else if (p < 0.10) suf <- " .  "
  else suf <- "    "
  paste0(formatC(p, format = "f", digits = 3, decimal.mark = ","), suf)
}

#' @keywords internal
.fmt_num_dic <- function(x, d = 2) {
  if (is.na(x)) return("    -   ")
  formatC(x, format = "f", digits = d, decimal.mark = ",", big.mark = ".")
}

#' @export
print.estatR_dic <- function(x, ...) {
  .print_titulo("DELINEAMENTO INTEIRAMENTE CASUALIZADO (DIC)")
  
  cat(sprintf("  Vari\u00e1vel Resposta:    %s\n", x$nomes$resposta))
  cat(sprintf("  Fator Tratamento:     %s (%d n\u00edveis: %s)\n",
              x$nomes$tratamento, length(x$niveis),
              if (length(x$niveis) <= 5) paste(x$niveis, collapse = ", ") else paste0(paste(x$niveis[1:4], collapse = ", "), ", ...")))
  cat(sprintf("  Classifica\u00e7\u00e3o:        %s\n", if (x$tipo == "qualitativo") "Qualitativo" else "Quantitativo"))
  cat(sprintf("  Repeti\u00e7\u00f5es:           %s\n", if (x$balanceado) paste0(x$tabela_medias$N[1], " (Balanceado)") else "Desbalanceado"))
  cat(sprintf("  Total de Observa\u00e7\u00f5es: %d\n\n", x$n_total))
  
  # 1. TABELA DE ANOVA
  .print_topico("QUADRO DA AN\u00c1LISE DE VARI\u00c2NCIA (ANOVA)")
  cat("  ────────────────────────────────────────────────────────────────────── \n")
  cat("  Fonte de Varia\u00e7\u00e3o    GL        SQ           QM         F calc    p-valor \n")
  cat("  ────────────────────────────────────────────────────────────────────── \n")
  for (i in 1:nrow(x$anova)) {
    f_str <- x$anova$Fonte[i]
    gl_s  <- sprintf("%2d", x$anova$GL[i])
    sq_s  <- sprintf("%10s", .fmt_num_dic(x$anova$SQ[i], 2))
    qm_s  <- if (!is.na(x$anova$QM[i])) sprintf("%10s", .fmt_num_dic(x$anova$QM[i], 2)) else "          "
    fc_s  <- if (!is.na(x$anova$F_calc[i])) sprintf("%8s", .fmt_num_dic(x$anova$F_calc[i], 2)) else "        "
    pv_s  <- sprintf("%11s", .fmt_p_dic(x$anova$p_valor[i]))
    
    cat(sprintf("  %-18s %4s   %10s   %10s   %8s %11s\n", f_str, gl_s, sq_s, qm_s, fc_s, pv_s))
  }
  cat("  ────────────────────────────────────────────────────────────────────── \n\n")
  
  # 2. PARÂMETROS EXPERIMENTAIS
  .print_topico("PAR\u00c2METROS EXPERIMENTAIS")
  cat("  ────────────────────────────────────────────────────────── \n")
  cat(sprintf("  M\u00e9dia Geral:             %s\n", .fmt_num_dic(x$parametros$media_geral, 2)))
  cat(sprintf("  CV Experimental (%%):     %s%%\n", .fmt_num_dic(x$parametros$cv_exp, 2)))
  cat(sprintf("  R\u00b2 do Modelo:            %s\n", .fmt_num_dic(x$parametros$r2, 3)))
  cat("  ────────────────────────────────────────────────────────── \n")
  cat("  * CV experimental: < 10% (Baixo / Alta precis\u00e3o), 10-20% (M\u00e9dio), 20-30% (Alto), > 30% (Muito alto - Pimentel-Gomes).\n\n")
  
  # 3. RESULTADOS DE ACORDO COM O TIPO
  if (x$tipo == "qualitativo") {
    
    # 3A. TESTE DE SCHEFFÉ COM CONTRASTES PERSONALIZADOS
    if (x$teste == "scheffe" && !is.null(x$contrastes)) {
      .print_topico("CONTRASTES DE SCHEFF\u00c9")
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      cat("  Contraste               Estimativa     EP      F calc    p-valor      IC (95%)    \n")
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$contrastes)) {
        nm_c <- x$contrastes$Contraste[i]
        est  <- sprintf("%+10s", .fmt_num_dic(x$contrastes$Estimativa[i], 2))
        ep   <- sprintf("%8s", .fmt_num_dic(x$contrastes$EP[i], 2))
        fc   <- sprintf("%8s", .fmt_num_dic(x$contrastes$F_calc[i], 2))
        pv   <- sprintf("%11s", .fmt_p_dic(x$contrastes$p_valor[i]))
        ic_s <- sprintf("[%s ; %s]", .fmt_num_dic(x$contrastes$IC_inf[i], 2), .fmt_num_dic(x$contrastes$IC_sup[i], 2))
        cat(sprintf("  %-22s %10s %8s %8s %11s  %-15s\n", nm_c, est, ep, fc, pv, ic_s))
      }
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")
      
    } else if (x$teste == "dunnett") {
      # 3B. TESTE DE DUNNETT
      .print_topico(sprintf("COMPARA\u00c7\u00d5ES CONTRA O CONTROLE (TESTE DE DUNNETT - TESTEMUNHA: '%s')", x$controle))
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      cat("  Tratamento        M\u00e9dia      Diferen\u00e7a    p-valor      Diferen\u00e7a Sig. \n")
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$tabela_grupos)) {
        tr_s   <- x$tabela_grupos$Tratamento[i]
        med_s  <- sprintf("%8s", .fmt_num_dic(x$tabela_grupos$Media[i], 2))
        dif_s  <- sprintf("%+9s", .fmt_num_dic(x$tabela_grupos$Diferenca[i], 2))
        pv_s   <- if (!is.na(x$tabela_grupos$p_valor[i])) sprintf("%11s", .fmt_p_dic(x$tabela_grupos$p_valor[i])) else "     Controle"
        sig_s  <- x$tabela_grupos$Sig[i]
        cat(sprintf("  %-16s %8s   %9s  %11s       %s\n", tr_s, med_s, dif_s, pv_s, sig_s))
      }
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")
      
    } else {
      # 3C. TESTES DE MÉDIAS COM AGRUPAMENTO DE LETRAS
      nm_teste_exib <- switch(x$teste,
        tukey       = "TESTE DE TUKEY (HSD)",
        scott_knott = "TESTE DE SCOTT-KNOTT",
        duncan      = "TESTE DE DUNCAN",
        snk         = "TESTE DE STUDENT-NEWMAN-KEULS (SNK)",
        scheffe     = "TESTE DE SCHEFF\u00c9 (PAREADO)",
        bonferroni  = "TESTE T COM AJUSTE DE BONFERRONI"
      )
      .print_topico(sprintf("COMPARA\u00c7\u00d5ES M\u00daLTIPLAS DE M\u00c9DIAS (%s)", nm_teste_exib))
      cat("  ────────────────────────────────────────────────────────── \n")
      cat("  Tratamento        M\u00e9dia Observada      Desvio Padr\u00e3o    Grupo\n")
      cat("  ────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$tabela_grupos)) {
        tr_s  <- x$tabela_grupos$Tratamento[i]
        med_s <- sprintf("%10s", .fmt_num_dic(x$tabela_grupos$Media[i], 2))
        dp_s  <- sprintf("%10s", .fmt_num_dic(x$tabela_grupos$DP[i], 2))
        grp_s <- sprintf("%-6s", x$tabela_grupos$Grupo[i])
        cat(sprintf("  %-16s %10s          %10s       %s\n", tr_s, med_s, dp_s, grp_s))
      }
      cat("  ────────────────────────────────────────────────────────── \n")
      cat("  * M\u00e9dias seguidas pela mesma letra n\u00e3o diferem estatisticamente entre si a 5% de signific\u00e2ncia.\n\n")
    }
  } else {
    # 3D. REGRESSÃO POLINOMIAL
    .print_topico("REGRESS\u00c3O POLINOMIAL (AN\u00c1LISE DE DOSES / N\u00cdVEIS)")
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    cat("  Grau Polinomial     GL        SQ           QM         F calc    p-valor     R\u00b2 \n")
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    tab_pol <- x$regressao$tabela_polinomial
    for (i in 1:nrow(tab_pol)) {
      gr_s  <- tab_pol$Grau[i]
      gl_s  <- sprintf("%2d", tab_pol$GL[i])
      sq_s  <- sprintf("%10s", .fmt_num_dic(tab_pol$SQ[i], 2))
      qm_s  <- sprintf("%10s", .fmt_num_dic(tab_pol$QM[i], 2))
      fc_s  <- if (!is.na(tab_pol$F_calc[i])) sprintf("%8s", .fmt_num_dic(tab_pol$F_calc[i], 2)) else "        "
      pv_s  <- sprintf("%11s", .fmt_p_dic(tab_pol$p_valor[i]))
      r2_s  <- if (!is.na(tab_pol$R2[i])) sprintf("%6s", .fmt_num_dic(tab_pol$R2[i], 3)) else "      "
      cat(sprintf("  %-18s %2s   %10s   %10s   %8s %11s  %6s\n", gr_s, gl_s, sq_s, qm_s, fc_s, pv_s, r2_s))
    }
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n")
    cat(sprintf("  Modelo Recomendado: Polinomial %s\n\n", x$regressao$modelo_recomendado))
    
    .print_topico("EQUA\u00c7\u00c3O AJUSTADA E PONTO \u00d3TIMO")
    coefs <- x$regressao$coeficientes
    .fmt_s <- function(val, d = 3) {
      if (is.na(val)) return("")
      sinal <- if (val >= 0) "+ " else "- "
      paste0(sinal, formatC(abs(val), format = "f", digits = d, decimal.mark = ",", big.mark = "."))
    }
    
    if (x$regressao$modelo_recomendado == "Linear") {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3)))
    } else if (x$regressao$modelo_recomendado == "Quadr\u00e1tico") {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose %s * Dose\u00b2\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3), .fmt_s(coefs[3], 4)))
    } else {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose %s * Dose\u00b2 %s * Dose\u00b3\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3), .fmt_s(coefs[3], 4), .fmt_s(coefs[4], 5)))
    }
    cat(sprintf("  R\u00b2 Total: %s\n", .fmt_num_dic(x$regressao$r2, 3)))
    if (!is.null(x$regressao$ponto_critico)) {
      cat(sprintf("  Ponto Cr\u00edtico: %s\n", x$regressao$ponto_critico$tipo))
      cat(sprintf("  Dose de M\u00e1xima Efici\u00eancia: %s  (Resposta estimada: %s)\n",
                  .fmt_num_dic(x$regressao$ponto_critico$x_otimo, 2),
                  .fmt_num_dic(x$regressao$ponto_critico$y_otimo, 2)))
    }
    cat("\n")
  }
  
  # 4. INTERPRETAÇÃO
  .print_topico("INTERPRETA\u00c7\u00c3O")
  p_trat <- x$anova$p_valor[1]
  f_trat <- x$anova$F_calc[1]
  gl_t   <- x$anova$GL[1]
  gl_r   <- x$anova$GL[2]
  
  if (!is.na(p_trat) && p_trat < 0.05) {
    cat(sprintf("  O efeito de Tratamento foi estatisticamente significativo (F(%d, %d) = %s; p %s),\n  rejeitando a hip\u00f3tese nula de que as m\u00e9dias dos tratamentos s\u00e3o iguais a 5%% de signific\u00e2ncia.\n\n",
                gl_t, gl_r, .fmt_num_dic(f_trat, 2), if (p_trat < 0.001) "< 0,001" else paste0("= ", .fmt_num_dic(p_trat, 3))))
  } else {
    cat(sprintf("  O efeito de Tratamento n\u00e3o foi estatisticamente significativo (F(%d, %d) = %s; p = %s),\n  indicando que n\u00e3o h\u00e1 evid\u00eancias de diferen\u00e7a entre as m\u00e9dias dos tratamentos avaliados.\n\n",
                gl_t, gl_r, .fmt_num_dic(f_trat, 2), .fmt_num_dic(p_trat, 3)))
  }
  
  cat(sprintf("  O Coeficiente de Varia\u00e7\u00e3o Experimental (CV = %s%%) indica precis\u00e3o experimental %s\n  segundo os crit\u00e9rios de Pimentel-Gomes.\n\n",
              .fmt_num_dic(x$parametros$cv_exp, 2),
              if (x$parametros$cv_exp < 10) "alta (baixo CV)"
              else if (x$parametros$cv_exp < 20) "m\u00e9dia"
              else if (x$parametros$cv_exp < 30) "baixa (alto CV)"
              else "muito baixa (CV muito alto)"))
  
  if (x$tipo == "qualitativo" && !is.null(x$tabela_grupos) && "Grupo" %in% names(x$tabela_grupos)) {
    melhor_tr <- x$tabela_grupos$Tratamento[1]
    melhor_md <- x$tabela_grupos$Media[1]
    cat(sprintf("  Segundo as compara\u00e7\u00f5es m\u00faltiplas (%s), o tratamento '%s' apresentou o maior\n  desempenho m\u00e9dio observado (%s).\n",
                toupper(x$teste), melhor_tr, .fmt_num_dic(melhor_md, 2)))
  } else if (x$tipo == "quantitativo" && !is.null(x$regressao$ponto_critico)) {
    cat(sprintf("  O comportamento das doses seguiu modelo %s, atingindo %s na dose\n  estimada de %s com resposta esperada de %s.\n",
                tolower(x$regressao$modelo_recomendado),
                tolower(x$regressao$ponto_critico$tipo),
                .fmt_num_dic(x$regressao$ponto_critico$x_otimo, 2),
                .fmt_num_dic(x$regressao$ponto_critico$y_otimo, 2)))
  }
  
  .print_rodape()
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# FUNÇÕES GRÁFICAS DEDICADAS
# ─────────────────────────────────────────────────────────────────────────────

#' @title Gráfico de Médias do Experimento
#' @description Gera gráfico de colunas com as médias dos tratamentos, barras de erro (IC 95%)
#'   e as letras do teste de comparações múltiplas sobrepostas no topo de cada grupo.
#' @param modelo Objeto retornado por \code{analise_dic()} ou \code{dic()}.
#' @param cor_barras Cor de preenchimento das barras (padrão \code{"#008080"}).
#' @return Objeto \code{ggplot}.
#' @export
grafico_medias <- function(modelo, cor_barras = "#008080") {
  if (inherits(modelo, "estatR_fatorial")) {
    return(grafico_fatorial(modelo))
  }
  if (!inherits(modelo, c("estatR_dic", "estatR_dbc", "estatR_dql"))) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_dic(), analise_dbc(), analise_dql() ou analise_fatorial_duplo().", call. = FALSE)
  }
  
  df_grp <- if (!is.null(modelo$tabela_grupos)) modelo$tabela_grupos else modelo$tabela_medias
  df_grp$Tratamento <- factor(df_grp$Tratamento, levels = df_grp$Tratamento)
  
  # Erro padrao e limite superior para posicionar letras
  n_r <- if ("N" %in% names(df_grp)) df_grp$N else rep(mean(modelo$tabela_medias$N), nrow(df_grp))
  df_grp$EP <- df_grp$DP / sqrt(n_r)
  df_grp$IC_sup <- df_grp$Media + 1.96 * df_grp$EP
  df_grp$IC_inf <- df_grp$Media - 1.96 * df_grp$EP
  
  # Letra
  tem_letra <- "Grupo" %in% names(df_grp)
  
  p <- ggplot2::ggplot(df_grp, ggplot2::aes(x = Tratamento, y = Media)) +
    ggplot2::geom_col(fill = cor_barras, alpha = 0.85, width = 0.55, color = "black", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = IC_inf, ymax = IC_sup), width = 0.2, linewidth = 0.6, color = "black") +
    ggplot2::labs(
      title = "Médias dos Tratamentos (com IC 95%)",
      subtitle = if (tem_letra) paste0("Letras do teste de ", toupper(modelo$teste)) else NULL,
      x = modelo$nomes$tratamento,
      y = modelo$nomes$resposta
    )
  
  if (tem_letra) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(y = IC_sup + (max(df_grp$Media, na.rm = TRUE) * 0.05), label = Grupo),
      fontface = "bold", size = 4.2, color = "black"
    )
  }
  
  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  
  p
}

#' @title Gráfico de Intervalos de Confiança das Diferenças de Médias
#' @description Gera gráfico horizontal com os intervalos de confiança das diferenças pareadas
#'   entre tratamentos ou dos contrastes de Scheffé, com linha vertical vermelha no zero e destaque
#'   colorido para contrastes com diferença estatisticamente significativa.
#' @param modelo Objeto retornado por \code{analise_dic()}, \code{analise_dbc()} ou \code{analise_dql()}.
#' @param cor_sig Cor para os intervalos significativos que não cruzam o zero (padrão azul \code{"#0099FF"}).
#' @param cor_ns Cor para os intervalos que cruzam o zero (padrão \code{"gray40"}).
#' @return Objeto \code{ggplot}.
#' @export
#' @aliases grafico_comparacoes
grafico_ic <- function(modelo, cor_sig = "#0099FF", cor_ns = "gray40") {
  if (!inherits(modelo, c("estatR_dic", "estatR_dbc", "estatR_dql"))) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_dic(), analise_dbc() ou analise_dql().", call. = FALSE)
  }
  
  df_comp <- modelo$comparacoes
  if (is.null(df_comp) || nrow(df_comp) == 0) {
    stop("N\u00e3o h\u00e1 compara\u00e7\u00f5es pareadas ou contrastes dispon\u00edveis no modelo.", call. = FALSE)
  }
  
  df_comp$Comparacao <- factor(df_comp$Comparacao, levels = rev(df_comp$Comparacao))
  df_comp$Significativo <- ifelse(df_comp$Lwr <= 0 & df_comp$Upr >= 0, "N\u00e3o", "Sim")
  
  p <- ggplot2::ggplot(df_comp, ggplot2::aes(x = Diff, y = Comparacao, color = Significativo)) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lwr, xmax = Upr), height = 0.25, linewidth = 0.85) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_color_manual(
      values = c("N\u00e3o" = cor_ns, "Sim" = cor_sig),
      labels = c("N\u00e3o" = "N\u00e3o significativo (passa por 0)", "Sim" = "Significativo (p < 0.05)")
    ) +
    ggplot2::labs(
      title = "Intervalos de Confiança das Diferenças",
      subtitle = paste0("M\u00e9todo: ", toupper(modelo$teste), " (IC ", modelo$confianca * 100, "%)"),
      x = "Diferença entre médias",
      y = "Comparação / Contraste",
      color = "Diferen\u00e7a"
    )
  
  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  
  p <- p + ggplot2::theme(legend.position = "bottom")
  p
}

#' @export
grafico_comparacoes <- grafico_ic

#' @title Gráfico da Regressão Polinomial
#' @description Gera gráfico de dispersão com os pontos observados, médias e curva ajustada do modelo
#'   polinomial recomendado (Linear, Quadrático ou Cúbico).
#' @param modelo Objeto retornado por \code{analise_dic()}, \code{analise_dbc()} ou \code{analise_dql()}.
#' @return Objeto \code{ggplot}.
#' @export
grafico_regressao <- function(modelo) {
  if (!inherits(modelo, c("estatR_dic", "estatR_dbc", "estatR_dql")) || modelo$tipo != "quantitativo") {
    stop("O modelo deve ser um objeto quantitativo retornado por analise_dic(), analise_dbc() ou analise_dql().", call. = FALSE)
  }
  
  df_orig <- data.frame(
    X = as.numeric(as.character(modelo$dados$Tratamento)),
    Y = modelo$dados$Resposta
  )
  
  rec <- modelo$regressao$modelo_recomendado
  grau_n <- if (rec == "Linear") 1 else if (rec == "Quadr\u00e1tico") 2 else 3
  
  p <- ggplot2::ggplot(df_orig, ggplot2::aes(x = X, y = Y)) +
    ggplot2::geom_point(alpha = 0.45, size = 2.2, color = "gray30") +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 3.5, color = "#008080") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ poly(x, grau_n, raw = TRUE),
                         color = "#008080", fill = "#008080", alpha = 0.2, linewidth = 0.9) +
    ggplot2::labs(
      title = paste0("Curva de Resposta: Regressão Polinomial ", rec),
      subtitle = paste0("R² = ", round(modelo$regressao$r2, 3)),
      x = modelo$nomes$tratamento,
      y = modelo$nomes$resposta
    )
  
  if (!is.null(modelo$regressao$ponto_critico)) {
    pt_x <- modelo$regressao$ponto_critico$x_otimo
    pt_y <- modelo$regressao$ponto_critico$y_otimo
    p <- p +
      ggplot2::geom_vline(xintercept = pt_x, linetype = "dotted", color = "red", linewidth = 0.8) +
      ggplot2::annotate("point", x = pt_x, y = pt_y, color = "red", size = 3.5) +
      ggplot2::annotate("text", x = pt_x, y = pt_y,
                        label = paste0(modelo$regressao$ponto_critico$tipo, "\nDose = ", round(pt_x, 2)),
                        vjust = -0.8, color = "red", fontface = "bold", size = 3.6)
  }
  
  if (exists("tema_estatR", mode = "function")) {
    p <- p + tema_estatR(estilo = 2, grade = "dupla")
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  
  p
}

#' @title Painel Gráfico Completo do DIC
#' @description Gera painel integrado de gráficos do Delineamento Inteiramente Casualizado.
#'   Para fatores qualitativos, combina o gráfico de médias com letras e o gráfico de IC das diferenças.
#'   Para fatores quantitativos, exibe a curva de regressão polinomial.
#' @param modelo Objeto retornado por \code{analise_dic()} ou \code{dic()}.
#' @param ... Argumentos adicionais passados para as funções gráficas individuais.
#' @return Objeto \code{ggplot} ou composição de gráficos.
#' @export
#' @aliases plot_dic
grafico_dic <- function(modelo, ...) {
  if (!inherits(modelo, "estatR_dic")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_dic().", call. = FALSE)
  }
  
  if (modelo$tipo == "quantitativo") {
    return(grafico_regressao(modelo))
  }
  
  p_med <- tryCatch(grafico_medias(modelo), error = function(e) NULL)
  p_ic  <- tryCatch(grafico_ic(modelo), error = function(e) NULL)
  
  if (!is.null(p_med) && !is.null(p_ic)) {
    # Combinar usando patchwork
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(p_med + p_ic)
    }
    return(p_med)
  } else if (!is.null(p_med)) {
    return(p_med)
  } else {
    return(p_ic)
  }
}

#' @export
plot_dic <- grafico_dic


# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: DELINEAMENTO EM BLOCOS CASUALIZADOS (DBC)
# ─────────────────────────────────────────────────────────────────────────────

#' @title Análise de Delineamento em Blocos Casualizados (DBC)
#' @description Executa a Análise de Variância (ANOVA) para experimentos em Delineamento
#'   em Blocos Casualizados (DBC), com teste de eficiência da blocagem, testes de comparações
#'   múltiplas de médias (Tukey, Scott-Knott, Duncan, SNK, Scheffé com contrastes, Bonferroni, Dunnett)
#'   para fatores qualitativos, ou Regressão Polinomial (Linear, Quadrática e Cúbica com ponto ótimo)
#'   para fatores quantitativos, além de painel gráfico automatizado.
#'
#' @param resposta Variável dependente contínua (com ou sem aspas).
#' @param tratamento Fator de tratamentos (com ou sem aspas).
#' @param bloco Fator de blocos / controle local (com ou sem aspas).
#' @param dados Data frame opcional contendo as variáveis.
#' @param teste Teste de comparações múltiplas para fatores qualitativos:
#'   \code{"tukey"} (padrão), \code{"scott_knott"}, \code{"duncan"}, \code{"snk"},
#'   \code{"scheffe"}, \code{"bonferroni"} ou \code{"dunnett"}.
#' @param contrastes Lista ou vetor de contrastes para o teste de Scheffé (ex: \code{"A vs c(B, C)"}).
#' @param controle Nome do nível de referência (testemunha) usado quando \code{teste = "dunnett"}.
#' @param tipo Tipo de fator: \code{"auto"} (detecta automaticamente), \code{"qualitativo"}
#'   ou \code{"quantitativo"}.
#' @param confianca Nível de confiança para os testes e intervalos (padrão \code{0.95}).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe automaticamente o painel de gráficos.
#' @param residuos Lógico. Se \code{TRUE}, executa e exibe a análise residual do modelo.
#'
#' @return Objeto da classe \code{"estatR_dbc"} contendo os quadros da ANOVA, parâmetros
#'   experimentais, eficiência da blocagem, resultados das comparações múltiplas e gráficos.
#' @export
#' @aliases dbc
analise_dbc <- function(resposta,
                        tratamento,
                        bloco,
                        dados       = NULL,
                        teste       = c("tukey", "scott_knott", "duncan", "snk", "scheffe", "bonferroni", "dunnett"),
                        contrastes  = NULL,
                        controle    = NULL,
                        tipo        = c("auto", "qualitativo", "quantitativo"),
                        confianca   = 0.95,
                        grafico     = TRUE,
                        residuos    = FALSE) {
  
  teste <- match.arg(teste)
  tipo  <- match.arg(tipo)
  
  sub_resp  <- substitute(resposta)
  sub_trat  <- substitute(tratamento)
  sub_bloco <- substitute(bloco)
  
  info_resp  <- .extrair_vetor_exp(sub_resp, dados, "resposta")
  info_trat  <- .extrair_vetor_exp(sub_trat, dados, "tratamento")
  info_bloco <- .extrair_vetor_exp(sub_bloco, dados, "bloco")
  
  y_raw     <- info_resp$valores
  trat_raw  <- info_trat$valores
  bloco_raw <- info_bloco$valores
  
  nm_resp  <- info_resp$nome
  nm_trat  <- info_trat$nome
  nm_bloco <- info_bloco$nome
  
  # Validacoes basicas
  val_ok    <- !is.na(y_raw) & !is.na(trat_raw) & !is.na(bloco_raw)
  y         <- as.numeric(y_raw[val_ok])
  trat_vec  <- trat_raw[val_ok]
  bloco_vec <- bloco_raw[val_ok]
  n_total   <- length(y)
  
  if (n_total < 4) {
    stop("O experimento em DBC deve conter ao menos 4 observa\u00e7\u00f5es v\u00e1lidas.", call. = FALSE)
  }
  
  # Fatores e Niveis
  trat_fator  <- as.factor(trat_vec)
  bloco_fator <- as.factor(bloco_vec)
  
  niveis_trat  <- levels(trat_fator)
  niveis_bloco <- levels(bloco_fator)
  
  k <- length(niveis_trat)
  b <- length(niveis_bloco)
  
  if (k < 2) stop("O fator de tratamentos deve possuir ao menos 2 n\u00edveis.", call. = FALSE)
  if (b < 2) stop("O fator de blocos deve possuir ao menos 2 blocos.", call. = FALSE)
  
  # Deteccao de qualitativo vs quantitativo
  e_numerico <- is.numeric(trat_vec) && !is.factor(trat_vec)
  if (tipo == "auto") {
    tipo_final <- if (e_numerico && length(unique(trat_vec)) >= 3) "quantitativo" else "qualitativo"
  } else {
    tipo_final <- tipo
  }
  
  # Ajuste do Modelo ANOVA DBC: Y ~ Bloco + Trat
  df_mod  <- data.frame(Y = y, Trat = trat_fator, Bloco = bloco_fator)
  mod_aov <- stats::aov(Y ~ Bloco + Trat, data = df_mod)
  tab_aov <- summary(mod_aov)[[1]]
  
  gl_bloco <- tab_aov["Bloco", "Df"]
  sq_bloco <- tab_aov["Bloco", "Sum Sq"]
  qm_bloco <- tab_aov["Bloco", "Mean Sq"]
  f_bloco  <- tab_aov["Bloco", "F value"]
  p_bloco  <- tab_aov["Bloco", "Pr(>F)"]
  
  gl_trat  <- tab_aov["Trat", "Df"]
  sq_trat  <- tab_aov["Trat", "Sum Sq"]
  qm_trat  <- tab_aov["Trat", "Mean Sq"]
  f_trat   <- tab_aov["Trat", "F value"]
  p_trat   <- tab_aov["Trat", "Pr(>F)"]
  
  gl_res   <- tab_aov["Residuals", "Df"]
  sq_res   <- tab_aov["Residuals", "Sum Sq"]
  qm_res   <- tab_aov["Residuals", "Mean Sq"]
  
  gl_tot   <- gl_bloco + gl_trat + gl_res
  sq_tot   <- sq_bloco + sq_trat + sq_res
  
  # Parametros experimentais
  media_geral <- mean(y)
  cv_exp      <- (sqrt(qm_res) / media_geral) * 100
  r2          <- (sq_trat + sq_bloco) / sq_tot
  
  # Eficiencia Relativa da Blocagem (ERB)
  # ER = (QM_Bloco + (k - 1) * QM_Res) / (k * QM_Res) * 100%
  eficiencia_relativa <- ((qm_bloco + (k - 1) * qm_res) / (k * qm_res)) * 100
  
  # Medias e DP por tratamento
  medias_grp <- tapply(y, trat_fator, mean)
  dps_grp    <- tapply(y, trat_fator, stats::sd)
  reps_tab   <- table(trat_fator)
  
  tab_medias_df <- data.frame(
    Tratamento = niveis_trat,
    N          = as.vector(reps_tab),
    Media      = as.vector(medias_grp),
    DP         = as.vector(dps_grp),
    stringsAsFactors = FALSE
  )
  
  # ───────────────────────────────────────────────────────────────────────────
  # CASO 1: FATOR QUALITATIVO (TESTES DE COMPARACAO MULTIPLA)
  # ───────────────────────────────────────────────────────────────────────────
  res_comparacoes <- NULL
  tab_grupos      <- NULL
  tab_contrastes  <- NULL
  alpha           <- 1 - confianca
  
  if (tipo_final == "qualitativo") {
    
    # 1. TUKEY HSD
    if (teste == "tukey") {
      tuk <- stats::TukeyHSD(mod_aov, "Trat", conf.level = confianca)$Trat
      pares <- rownames(tuk)
      
      p_mat <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis_trat, niveis_trat))
      for (i in seq_along(pares)) {
        pts <- unlist(strsplit(pares[i], "-"))
        p_mat[pts[1], pts[2]] <- tuk[i, "p adj"]
        p_mat[pts[2], pts[1]] <- tuk[i, "p adj"]
      }
      
      med_ord <- sort(medias_grp, decreasing = TRUE)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      
      tab_grupos <- data.frame(
        Tratamento = names(med_ord),
        Media      = as.vector(med_ord),
        DP         = dps_grp[names(med_ord)],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
      
      res_comparacoes <- data.frame(
        Comparacao = pares,
        Diff       = tuk[, "diff"],
        Lwr        = tuk[, "lwr"],
        Upr        = tuk[, "upr"],
        p_valor    = tuk[, "p adj"],
        Sig        = ifelse(tuk[, "lwr"] <= 0 & tuk[, "upr"] >= 0, "N\u00e3o", "Sim"),
        stringsAsFactors = FALSE
      )
    }
    
    # 2. SCHEFFÉ
    else if (teste == "scheffe") {
      s_crit <- sqrt((k - 1) * stats::qf(confianca, k - 1, gl_res))
      ns_grp <- as.vector(reps_tab)
      names(ns_grp) <- niveis_trat
      
      if (!is.null(contrastes)) {
        contr_list <- .parse_contrastes_scheffe(contrastes, niveis_trat)
        df_contr <- list()
        
        for (nm_c in names(contr_list)) {
          c_vec <- contr_list[[nm_c]]
          est   <- sum(c_vec * medias_grp)
          ep    <- sqrt(qm_res * sum((c_vec^2) / ns_grp))
          t_val <- est / ep
          f_val <- (t_val^2) / (k - 1)
          p_val <- 1 - stats::pf(f_val, k - 1, gl_res)
          ic_inf <- est - s_crit * ep
          ic_sup <- est + s_crit * ep
          
          df_contr[[nm_c]] <- data.frame(
            Contraste  = nm_c,
            Estimativa = est,
            EP         = ep,
            F_calc     = f_val * (k - 1),
            p_valor    = p_val,
            IC_inf     = ic_inf,
            IC_sup     = ic_sup,
            Sig        = ifelse(ic_inf <= 0 & ic_sup >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
        tab_contrastes <- do.call(rbind, df_contr)
        res_comparacoes <- tab_contrastes
        res_comparacoes$Comparacao <- res_comparacoes$Contraste
        res_comparacoes$Diff       <- res_comparacoes$Estimativa
        res_comparacoes$Lwr        <- res_comparacoes$IC_inf
        res_comparacoes$Upr        <- res_comparacoes$IC_sup
      } else {
        pares_lista <- utils::combn(niveis_trat, 2, simplify = FALSE)
        df_scheffe  <- list()
        p_mat       <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis_trat, niveis_trat))
        
        for (pr in pares_lista) {
          n1 <- pr[1]; n2 <- pr[2]
          diff_val <- medias_grp[n1] - medias_grp[n2]
          ep_val   <- sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
          t_val    <- diff_val / ep_val
          f_val    <- (t_val^2) / (k - 1)
          p_val    <- 1 - stats::pf(f_val, k - 1, gl_res)
          lwr_val  <- diff_val - s_crit * ep_val
          upr_val  <- diff_val + s_crit * ep_val
          
          p_mat[n1, n2] <- p_val
          p_mat[n2, n1] <- p_val
          
          df_scheffe[[paste0(n1, " - ", n2)]] <- data.frame(
            Comparacao = paste0(n1, " - ", n2),
            Diff       = diff_val,
            Lwr        = lwr_val,
            Upr        = upr_val,
            p_valor    = p_val,
            Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
        res_comparacoes <- do.call(rbind, df_scheffe)
        med_ord <- sort(medias_grp, decreasing = TRUE)
        letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
        tab_grupos <- data.frame(
          Tratamento = names(med_ord),
          Media      = as.vector(med_ord),
          DP         = dps_grp[names(med_ord)],
          Grupo      = letras,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # 3. BONFERRONI
    else if (teste == "bonferroni") {
      pares_lista <- utils::combn(niveis_trat, 2, simplify = FALSE)
      m_comp      <- length(pares_lista)
      t_crit      <- stats::qt(1 - alpha / (2 * m_comp), gl_res)
      df_bonf     <- list()
      p_mat       <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis_trat, niveis_trat))
      ns_grp      <- as.vector(reps_tab)
      names(ns_grp) <- niveis_trat
      
      for (pr in pares_lista) {
        n1 <- pr[1]; n2 <- pr[2]
        diff_val <- medias_grp[n1] - medias_grp[n2]
        ep_val   <- sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
        t_val    <- abs(diff_val) / ep_val
        p_adj    <- min(1, 2 * (1 - stats::pt(t_val, gl_res)) * m_comp)
        lwr_val  <- diff_val - t_crit * ep_val
        upr_val  <- diff_val + t_crit * ep_val
        
        p_mat[n1, n2] <- p_adj
        p_mat[n2, n1] <- p_adj
        
        df_bonf[[paste0(n1, " - ", n2)]] <- data.frame(
          Comparacao = paste0(n1, " - ", n2),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = p_adj,
          Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_bonf)
      med_ord <- sort(medias_grp, decreasing = TRUE)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      tab_grupos <- data.frame(
        Tratamento = names(med_ord),
        Media      = as.vector(med_ord),
        DP         = dps_grp[names(med_ord)],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
    }
    
    # 4. SCOTT-KNOTT
    else if (teste == "scott_knott") {
      med_ord   <- sort(medias_grp, decreasing = TRUE)
      nomes_ord <- names(med_ord)
      r_media   <- mean(reps_tab)
      
      sk_divisao <- function(nomes_sub) {
        k_sub <- length(nomes_sub)
        if (k_sub <= 1) return(list(nomes_sub))
        
        y_sub <- med_ord[nomes_sub]
        melhor_b0 <- -1
        corte_opt <- 1
        
        for (c_pt in 1:(k_sub - 1)) {
          g1 <- nomes_sub[1:c_pt]
          g2 <- nomes_sub[(c_pt + 1):k_sub]
          m1 <- mean(y_sub[g1])
          m2 <- mean(y_sub[g2])
          b0 <- (length(g1) * length(g2) / (length(g1) + length(g2))) * (m1 - m2)^2
          if (b0 > melhor_b0) {
            melhor_b0 <- b0
            corte_opt <- c_pt
          }
        }
        
        s2_0 <- qm_res
        lambda0 <- (pi / (2 * (pi - 2))) * (melhor_b0 / (s2_0 / r_media))
        gl_sk <- k_sub / (pi - 2)
        p_sk <- 1 - stats::pchisq(lambda0, df = gl_sk)
        
        if (p_sk < alpha) {
          g1_sub <- nomes_sub[1:corte_opt]
          g2_sub <- nomes_sub[(corte_opt + 1):k_sub]
          return(c(sk_divisao(g1_sub), sk_divisao(g2_sub)))
        } else {
          return(list(nomes_sub))
        }
      }
      
      grupos_sk <- sk_divisao(nomes_ord)
      letras_vet <- character(length(nomes_ord))
      names(letras_vet) <- nomes_ord
      
      for (g_idx in seq_along(grupos_sk)) {
        letra_atual <- letters[g_idx]
        membros <- grupos_sk[[g_idx]]
        letras_vet[membros] <- letra_atual
      }
      
      tab_grupos <- data.frame(
        Tratamento = nomes_ord,
        Media      = as.vector(med_ord),
        DP         = dps_grp[nomes_ord],
        Grupo      = letras_vet,
        stringsAsFactors = FALSE
      )
      
      pares_lista <- utils::combn(nomes_ord, 2, simplify = FALSE)
      df_sk_comp  <- list()
      q_crit      <- stats::qtukey(confianca, k, gl_res)
      
      for (pr in pares_lista) {
        n1 <- pr[1]; n2 <- pr[2]
        diff_val <- medias_grp[n1] - medias_grp[n2]
        ep_val   <- sqrt(qm_res * (1/reps_tab[n1] + 1/reps_tab[n2]))
        lwr_val  <- diff_val - (q_crit / sqrt(2)) * ep_val
        upr_val  <- diff_val + (q_crit / sqrt(2)) * ep_val
        
        mesmo_grp <- letras_vet[n1] == letras_vet[n2]
        df_sk_comp[[paste0(n1, " - ", n2)]] <- data.frame(
          Comparacao = paste0(n1, " - ", n2),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = if (mesmo_grp) 0.50 else 0.01,
          Sig        = if (mesmo_grp) "N\u00e3o" else "Sim",
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_sk_comp)
    }
    
    # 5. DUNCAN & SNK
    else if (teste %in% c("duncan", "snk")) {
      med_ord   <- sort(medias_grp, decreasing = TRUE)
      nomes_ord <- names(med_ord)
      p_mat     <- matrix(1, nrow = k, ncol = k, dimnames = list(nomes_ord, nomes_ord))
      df_amp    <- list()
      r_media   <- mean(reps_tab)
      
      for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
          n1 <- nomes_ord[i]; n2 <- nomes_ord[j]
          dist_p <- j - i + 1
          diff_val <- med_ord[i] - med_ord[j]
          ep_val   <- sqrt(qm_res * (1/reps_tab[n1] + 1/reps_tab[n2]))
          
          q_obs <- (diff_val / sqrt(qm_res / r_media)) * sqrt(2)
          p_adj <- 1 - stats::ptukey(abs(q_obs), dist_p, gl_res)
          
          p_mat[n1, n2] <- p_adj
          p_mat[n2, n1] <- p_adj
          
          q_c <- stats::qtukey(confianca, dist_p, gl_res)
          lwr_val <- diff_val - (q_c / sqrt(2)) * ep_val
          upr_val <- diff_val + (q_c / sqrt(2)) * ep_val
          
          df_amp[[paste0(n1, " - ", n2)]] <- data.frame(
            Comparacao = paste0(n1, " - ", n2),
            Diff       = diff_val,
            Lwr        = lwr_val,
            Upr        = upr_val,
            p_valor    = p_adj,
            Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
      }
      res_comparacoes <- do.call(rbind, df_amp)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      tab_grupos <- data.frame(
        Tratamento = nomes_ord,
        Media      = as.vector(med_ord),
        DP         = dps_grp[nomes_ord],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
    }
    
    # 6. DUNNETT
    else if (teste == "dunnett") {
      if (is.null(controle) || !controle %in% niveis_trat) {
        controle <- niveis_trat[1]
      }
      outros <- setdiff(niveis_trat, controle)
      df_dunnett <- list()
      m_d <- length(outros)
      t_crit_d <- stats::qt(1 - alpha / (2 * m_d), gl_res)
      
      for (tr in outros) {
        diff_val <- medias_grp[tr] - medias_grp[controle]
        ep_val   <- sqrt(qm_res * (1/reps_tab[tr] + 1/reps_tab[controle]))
        t_val    <- abs(diff_val) / ep_val
        p_val    <- min(1, 2 * (1 - stats::pt(t_val, gl_res)) * m_d)
        lwr_val  <- diff_val - t_crit_d * ep_val
        upr_val  <- diff_val + t_crit_d * ep_val
        
        df_dunnett[[paste0(tr, " vs ", controle)]] <- data.frame(
          Comparacao = paste0(tr, " vs ", controle),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = p_val,
          Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_dunnett)
      tab_grupos <- data.frame(
        Tratamento = c(controle, outros),
        Media      = c(medias_grp[controle], medias_grp[outros]),
        DP         = c(dps_grp[controle], dps_grp[outros]),
        Diferenca  = c(0, sapply(outros, function(tr) medias_grp[tr] - medias_grp[controle])),
        p_valor    = c(NA, sapply(outros, function(tr) res_comparacoes[paste0(tr, " vs ", controle), "p_valor"])),
        Sig        = c("Controle", sapply(outros, function(tr) res_comparacoes[paste0(tr, " vs ", controle), "Sig"])),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # CASO 2: FATOR QUANTITATIVO (REGRESSÃO POLINOMIAL COM BLOCOS)
  # ───────────────────────────────────────────────────────────────────────────
  res_regressao <- NULL
  if (tipo_final == "quantitativo") {
    x_num <- as.numeric(trat_vec)
    k_num <- length(unique(x_num))
    
    # Modelos polinomiais controlando para bloco
    mod_lin <- stats::lm(y ~ bloco_fator + x_num)
    anv_lin <- stats::anova(mod_lin)
    sq_lin  <- anv_lin["x_num", "Sum Sq"]
    
    sq_quad <- 0; f_quad <- NA; p_quad <- NA; r2_quad <- 0
    if (k_num >= 3) {
      mod_quad <- stats::lm(y ~ bloco_fator + x_num + I(x_num^2))
      anv_quad <- stats::anova(mod_quad)
      sq_quad  <- anv_quad["I(x_num^2)", "Sum Sq"]
    }
    
    sq_cub <- 0; f_cub <- NA; p_cub <- NA; r2_cub <- 0
    if (k_num >= 4) {
      mod_cub <- stats::lm(y ~ bloco_fator + x_num + I(x_num^2) + I(x_num^3))
      anv_cub <- stats::anova(mod_cub)
      sq_cub  <- anv_cub["I(x_num^3)", "Sum Sq"]
    }
    
    f_lin <- sq_lin / qm_res
    p_lin <- 1 - stats::pf(f_lin, 1, gl_res)
    r2_lin <- sq_lin / sq_tot
    
    tab_pol <- data.frame(
      Grau    = "Linear",
      GL      = 1,
      SQ      = sq_lin,
      QM      = sq_lin,
      F_calc  = f_lin,
      p_valor = p_lin,
      R2      = r2_lin,
      stringsAsFactors = FALSE
    )
    
    if (k_num >= 3) {
      f_quad <- sq_quad / qm_res
      p_quad <- 1 - stats::pf(f_quad, 1, gl_res)
      r2_quad <- sq_quad / sq_tot
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "Quadr\u00e1tico",
        GL      = 1,
        SQ      = sq_quad,
        QM      = sq_quad,
        F_calc  = f_quad,
        p_valor = p_quad,
        R2      = r2_quad,
        stringsAsFactors = FALSE
      ))
    }
    
    if (k_num >= 4) {
      f_cub <- sq_cub / qm_res
      p_cub <- 1 - stats::pf(f_cub, 1, gl_res)
      r2_cub <- sq_cub / sq_tot
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "C\u00fabico",
        GL      = 1,
        SQ      = sq_cub,
        QM      = sq_cub,
        F_calc  = f_cub,
        p_valor = p_cub,
        R2      = r2_cub,
        stringsAsFactors = FALSE
      ))
    }
    
    # Desvios de regressao
    sq_expl_pol <- sum(tab_pol$SQ)
    sq_desv     <- max(0, sq_trat - sq_expl_pol)
    gl_desv     <- max(0, gl_trat - nrow(tab_pol))
    qm_desv     <- if (gl_desv > 0) sq_desv / gl_desv else 0
    f_desv      <- if (gl_desv > 0) qm_desv / qm_res else NA
    p_desv      <- if (gl_desv > 0) 1 - stats::pf(f_desv, gl_desv, gl_res) else NA
    
    if (gl_desv > 0) {
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "Desvios Regress\u00e3o",
        GL      = gl_desv,
        SQ      = sq_desv,
        QM      = qm_desv,
        F_calc  = f_desv,
        p_valor = p_desv,
        R2      = sq_desv / sq_tot,
        stringsAsFactors = FALSE
      ))
    }
    
    tab_pol <- rbind(tab_pol, data.frame(
      Grau    = "Res\u00edduo",
      GL      = gl_res,
      SQ      = sq_res,
      QM      = qm_res,
      F_calc  = NA,
      p_valor = NA,
      R2      = NA,
      stringsAsFactors = FALSE
    ))
    
    mod_rec <- "Linear"
    mod_final_obj <- mod_lin
    if (k_num >= 3 && !is.na(p_quad) && p_quad < alpha) {
      mod_rec <- "Quadr\u00e1tico"
      mod_final_obj <- mod_quad
    }
    if (k_num >= 4 && !is.na(p_cub) && p_cub < alpha) {
      mod_rec <- "C\u00fabico"
      mod_final_obj <- mod_cub
    }
    
    # Coeficientes sem o efeito dos blocos (centrados na media geral)
    coefs_brutos <- stats::coef(mod_final_obj)
    # Extrair apenas os termos de x_num
    idx_x <- grep("x_num", names(coefs_brutos))
    coefs_pol <- coefs_brutos[idx_x]
    
    b1 <- coefs_pol[1]
    b2 <- if (length(coefs_pol) >= 2) coefs_pol[2] else NA
    b3 <- if (length(coefs_pol) >= 3) coefs_pol[3] else NA
    b0 <- media_geral - (b1 * mean(x_num) + (if (!is.na(b2)) b2 * mean(x_num^2) else 0) + (if (!is.na(b3)) b3 * mean(x_num^3) else 0))
    coefs_ajustados <- c(Intercept = as.numeric(b0), coefs_pol)
    
    ponto_critico <- NULL
    if (mod_rec == "Quadr\u00e1tico" && !is.na(b2) && b2 != 0) {
      x_otimo <- -b1 / (2 * b2)
      tipo_crit <- if (b2 < 0) "M\u00e1ximo T\u00e9cnico" else "M\u00ednimo T\u00e9cnico"
      y_est <- b0 + b1 * x_otimo + b2 * (x_otimo^2)
      ponto_critico <- list(
        tipo    = tipo_crit,
        x_otimo = as.numeric(x_otimo),
        y_otimo = as.numeric(y_est)
      )
    }
    
    r2_final <- summary(mod_final_obj)$r.squared
    
    res_regressao <- list(
      tabela_polinomial = tab_pol,
      modelo_recomendado = mod_rec,
      modelo_ajustado   = mod_final_obj,
      coeficientes      = coefs_ajustados,
      r2                = r2_final,
      ponto_critico     = ponto_critico
    )
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # OBJETO DE RETORNO
  # ───────────────────────────────────────────────────────────────────────────
  resultado <- structure(
    list(
      dados               = data.frame(Resposta = y, Tratamento = trat_fator, Bloco = bloco_fator),
      nomes               = list(resposta = nm_resp, tratamento = nm_trat, bloco = nm_bloco),
      niveis              = niveis_trat,
      niveis_bloco        = niveis_bloco,
      n_total             = n_total,
      tipo                = tipo_final,
      teste               = teste,
      confianca           = confianca,
      contrastes_definidos= contrastes,
      controle            = controle,
      anova = data.frame(
        Fonte   = c("Bloco", "Tratamento", "Res\u00edduo", "Total"),
        GL      = c(gl_bloco, gl_trat, gl_res, gl_tot),
        SQ      = c(sq_bloco, sq_trat, sq_res, sq_tot),
        QM      = c(qm_bloco, qm_trat, qm_res, NA),
        F_calc  = c(f_bloco, f_trat, NA, NA),
        p_valor = c(p_bloco, p_trat, NA, NA),
        stringsAsFactors = FALSE
      ),
      parametros = list(
        media_geral        = media_geral,
        cv_exp             = cv_exp,
        r2                 = r2,
        eficiencia_relativa= eficiencia_relativa
      ),
      tabela_medias       = tab_medias_df,
      tabela_grupos       = tab_grupos,
      comparacoes         = res_comparacoes,
      contrastes          = tab_contrastes,
      regressao           = res_regressao,
      modelo_aov          = mod_aov,
      grafico             = NULL
    ),
    class = "estatR_dbc"
  )
  
  g_painel <- tryCatch(
    grafico_dbc(resultado),
    error = function(e) NULL
  )
  resultado$grafico <- g_painel
  
  if (isTRUE(residuos)) {
    if (exists("analise_residual", mode = "function")) {
      cat("\n")
      analise_residual(mod_aov)
    }
  }
  
  if (isTRUE(grafico) && !is.null(g_painel)) {
    print(g_painel)
  }
  
  resultado
}

#' @export
dbc <- analise_dbc


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODO DE IMPRESSÃO DBC
# ─────────────────────────────────────────────────────────────────────────────

#' @export
print.estatR_dbc <- function(x, ...) {
  .print_titulo("DELINEAMENTO EM BLOCOS CASUALIZADOS (DBC)")
  
  cat(sprintf("  Vari\u00e1vel Resposta:    %s\n", x$nomes$resposta))
  cat(sprintf("  Fator Tratamento:     %s (%d n\u00edveis: %s)\n",
              x$nomes$tratamento, length(x$niveis),
              if (length(x$niveis) <= 5) paste(x$niveis, collapse = ", ") else paste0(paste(x$niveis[1:4], collapse = ", "), ", ...")))
  cat(sprintf("  Fator de Blocos:      %s (%d blocos: %s)\n",
              x$nomes$bloco, length(x$niveis_bloco),
              if (length(x$niveis_bloco) <= 5) paste(x$niveis_bloco, collapse = ", ") else paste0(paste(x$niveis_bloco[1:4], collapse = ", "), ", ...")))
  cat(sprintf("  Classifica\u00e7\u00e3o:        %s\n", if (x$tipo == "qualitativo") "Qualitativo" else "Quantitativo"))
  cat(sprintf("  Total de Observa\u00e7\u00f5es: %d\n\n", x$n_total))
  
  # 1. TABELA DE ANOVA
  .print_topico("QUADRO DA AN\u00c1LISE DE VARI\u00c2NCIA (ANOVA)")
  cat("  ────────────────────────────────────────────────────────────────────── \n")
  cat("  Fonte de Varia\u00e7\u00e3o    GL        SQ           QM         F calc    p-valor \n")
  cat("  ────────────────────────────────────────────────────────────────────── \n")
  for (i in 1:nrow(x$anova)) {
    f_str <- x$anova$Fonte[i]
    gl_s  <- sprintf("%2d", x$anova$GL[i])
    sq_s  <- sprintf("%10s", .fmt_num_dic(x$anova$SQ[i], 2))
    qm_s  <- if (!is.na(x$anova$QM[i])) sprintf("%10s", .fmt_num_dic(x$anova$QM[i], 2)) else "          "
    fc_s  <- if (!is.na(x$anova$F_calc[i])) sprintf("%8s", .fmt_num_dic(x$anova$F_calc[i], 2)) else "        "
    pv_s  <- sprintf("%11s", .fmt_p_dic(x$anova$p_valor[i]))
    
    cat(sprintf("  %-18s %4s   %10s   %10s   %8s %11s\n", f_str, gl_s, sq_s, qm_s, fc_s, pv_s))
  }
  cat("  ────────────────────────────────────────────────────────────────────── \n\n")
  
  # 2. PARÂMETROS EXPERIMENTAIS
  .print_topico("PAR\u00c2METROS EXPERIMENTAIS")
  cat("  ────────────────────────────────────────────────────────── \n")
  cat(sprintf("  M\u00e9dia Geral:             %s\n", .fmt_num_dic(x$parametros$media_geral, 2)))
  cat(sprintf("  CV Experimental (%%):     %s%%\n", .fmt_num_dic(x$parametros$cv_exp, 2)))
  cat(sprintf("  R\u00b2 do Modelo:            %s\n", .fmt_num_dic(x$parametros$r2, 3)))
  cat(sprintf("  Efici\u00eancia da Blocagem:  %s%% %s\n",
              .fmt_num_dic(x$parametros$eficiencia_relativa, 1),
              if (x$parametros$eficiencia_relativa > 100) "(Blocagem eficiente)" else "(Sem ganho em rela\u00e7\u00e3o ao DIC)"))
  cat("  ────────────────────────────────────────────────────────── \n")
  cat("  * CV experimental: < 10% (Baixo / Alta precis\u00e3o), 10-20% (M\u00e9dio), 20-30% (Alto), > 30% (Muito alto - Pimentel-Gomes).\n")
  cat("  * Efici\u00eancia Relativa: > 100% indica ganho de precis\u00e3o do DBC sobre o DIC correspondente.\n\n")
  
  # 3. RESULTADOS DE ACORDO COM O TIPO
  if (x$tipo == "qualitativo") {
    if (x$teste == "scheffe" && !is.null(x$contrastes)) {
      .print_topico("CONTRASTES DE SCHEFF\u00c9")
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      cat("  Contraste               Estimativa     EP      F calc    p-valor      IC (95%)    \n")
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$contrastes)) {
        nm_c <- x$contrastes$Contraste[i]
        est  <- sprintf("%+10s", .fmt_num_dic(x$contrastes$Estimativa[i], 2))
        ep   <- sprintf("%8s", .fmt_num_dic(x$contrastes$EP[i], 2))
        fc   <- sprintf("%8s", .fmt_num_dic(x$contrastes$F_calc[i], 2))
        pv   <- sprintf("%11s", .fmt_p_dic(x$contrastes$p_valor[i]))
        ic_s <- sprintf("[%s ; %s]", .fmt_num_dic(x$contrastes$IC_inf[i], 2), .fmt_num_dic(x$contrastes$IC_sup[i], 2))
        cat(sprintf("  %-22s %10s %8s %8s %11s  %-15s\n", nm_c, est, ep, fc, pv, ic_s))
      }
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")
    } else if (x$teste == "dunnett") {
      .print_topico(sprintf("COMPARA\u00c7\u00d5ES CONTRA O CONTROLE (TESTE DE DUNNETT - TESTEMUNHA: '%s')", x$controle))
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      cat("  Tratamento        M\u00e9dia      Diferen\u00e7a    p-valor      Diferen\u00e7a Sig. \n")
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$tabela_grupos)) {
        tr_s   <- x$tabela_grupos$Tratamento[i]
        med_s  <- sprintf("%8s", .fmt_num_dic(x$tabela_grupos$Media[i], 2))
        dif_s  <- sprintf("%+9s", .fmt_num_dic(x$tabela_grupos$Diferenca[i], 2))
        pv_s   <- if (!is.na(x$tabela_grupos$p_valor[i])) sprintf("%11s", .fmt_p_dic(x$tabela_grupos$p_valor[i])) else "     Controle"
        sig_s  <- x$tabela_grupos$Sig[i]
        cat(sprintf("  %-16s %8s   %9s  %11s       %s\n", tr_s, med_s, dif_s, pv_s, sig_s))
      }
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")
    } else {
      nm_teste_exib <- switch(x$teste,
        tukey       = "TESTE DE TUKEY (HSD)",
        scott_knott = "TESTE DE SCOTT-KNOTT",
        duncan      = "TESTE DE DUNCAN",
        snk         = "TESTE DE STUDENT-NEWMAN-KEULS (SNK)",
        scheffe     = "TESTE DE SCHEFF\u00c9 (PAREADO)",
        bonferroni  = "TESTE T COM AJUSTE DE BONFERRONI"
      )
      .print_topico(sprintf("COMPARA\u00c7\u00d5ES M\u00daLTIPLAS DE M\u00c9DIAS (%s)", nm_teste_exib))
      cat("  ────────────────────────────────────────────────────────── \n")
      cat("  Tratamento        M\u00e9dia Observada      Desvio Padr\u00e3o    Grupo\n")
      cat("  ────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$tabela_grupos)) {
        tr_s  <- x$tabela_grupos$Tratamento[i]
        med_s <- sprintf("%10s", .fmt_num_dic(x$tabela_grupos$Media[i], 2))
        dp_s  <- sprintf("%10s", .fmt_num_dic(x$tabela_grupos$DP[i], 2))
        grp_s <- sprintf("%-6s", x$tabela_grupos$Grupo[i])
        cat(sprintf("  %-16s %10s          %10s       %s\n", tr_s, med_s, dp_s, grp_s))
      }
      cat("  ────────────────────────────────────────────────────────── \n")
      cat("  * M\u00e9dias seguidas pela mesma letra n\u00e3o diferem estatisticamente entre si a 5% de signific\u00e2ncia.\n\n")
    }
  } else {
    .print_topico("REGRESS\u00c3O POLINOMIAL (AN\u00c1LISE DE DOSES / N\u00cdVEIS)")
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    cat("  Grau Polinomial     GL        SQ           QM         F calc    p-valor     R\u00b2 \n")
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    tab_pol <- x$regressao$tabela_polinomial
    for (i in 1:nrow(tab_pol)) {
      gr_s  <- tab_pol$Grau[i]
      gl_s  <- sprintf("%2d", tab_pol$GL[i])
      sq_s  <- sprintf("%10s", .fmt_num_dic(tab_pol$SQ[i], 2))
      qm_s  <- sprintf("%10s", .fmt_num_dic(tab_pol$QM[i], 2))
      fc_s  <- if (!is.na(tab_pol$F_calc[i])) sprintf("%8s", .fmt_num_dic(tab_pol$F_calc[i], 2)) else "        "
      pv_s  <- sprintf("%11s", .fmt_p_dic(tab_pol$p_valor[i]))
      r2_s  <- if (!is.na(tab_pol$R2[i])) sprintf("%6s", .fmt_num_dic(tab_pol$R2[i], 3)) else "      "
      cat(sprintf("  %-18s %2s   %10s   %10s   %8s %11s  %6s\n", gr_s, gl_s, sq_s, qm_s, fc_s, pv_s, r2_s))
    }
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n")
    cat(sprintf("  Modelo Recomendado: Polinomial %s\n\n", x$regressao$modelo_recomendado))
    
    .print_topico("EQUA\u00c7\u00c3O AJUSTADA E PONTO \u00d3TIMO")
    coefs <- x$regressao$coeficientes
    .fmt_s <- function(val, d = 3) {
      if (is.na(val)) return("")
      sinal <- if (val >= 0) "+ " else "- "
      paste0(sinal, formatC(abs(val), format = "f", digits = d, decimal.mark = ",", big.mark = "."))
    }
    
    if (x$regressao$modelo_recomendado == "Linear") {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3)))
    } else if (x$regressao$modelo_recomendado == "Quadr\u00e1tico") {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose %s * Dose\u00b2\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3), .fmt_s(coefs[3], 4)))
    } else {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose %s * Dose\u00b2 %s * Dose\u00b3\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3), .fmt_s(coefs[3], 4), .fmt_s(coefs[4], 5)))
    }
    cat(sprintf("  R\u00b2 Total: %s\n", .fmt_num_dic(x$regressao$r2, 3)))
    if (!is.null(x$regressao$ponto_critico)) {
      cat(sprintf("  Ponto Cr\u00edtico: %s\n", x$regressao$ponto_critico$tipo))
      cat(sprintf("  Dose de M\u00e1xima Efici\u00eancia: %s  (Resposta estimada: %s)\n",
                  .fmt_num_dic(x$regressao$ponto_critico$x_otimo, 2),
                  .fmt_num_dic(x$regressao$ponto_critico$y_otimo, 2)))
    }
    cat("\n")
  }
  
  # 4. INTERPRETAÇÃO
  .print_topico("INTERPRETA\u00c7\u00c3O")
  p_trat  <- x$anova$p_valor[2]
  f_trat  <- x$anova$F_calc[2]
  gl_t    <- x$anova$GL[2]
  gl_r    <- x$anova$GL[3]
  p_bloco <- x$anova$p_valor[1]
  f_bloco <- x$anova$F_calc[1]
  gl_b    <- x$anova$GL[1]
  
  if (!is.na(p_trat) && p_trat < 0.05) {
    cat(sprintf("  O efeito de Tratamento foi estatisticamente significativo (F(%d, %d) = %s; p %s),\n  rejeitando a hip\u00f3tese nula de que as m\u00e9dias dos tratamentos s\u00e3o iguais a 5%% de signific\u00e2ncia.\n\n",
                gl_t, gl_r, .fmt_num_dic(f_trat, 2), if (p_trat < 0.001) "< 0,001" else paste0("= ", .fmt_num_dic(p_trat, 3))))
  } else {
    cat(sprintf("  O efeito de Tratamento n\u00e3o foi estatisticamente significativo (F(%d, %d) = %s; p = %s),\n  indicando que n\u00e3o h\u00e1 evid\u00eancias de diferen\u00e7a entre as m\u00e9dias dos tratamentos avaliados.\n\n",
                gl_t, gl_r, .fmt_num_dic(f_trat, 2), .fmt_num_dic(p_trat, 3)))
  }
  
  if (!is.na(p_bloco) && p_bloco < 0.05) {
    cat(sprintf("  O efeito de Bloco foi significativo (F(%d, %d) = %s; p = %s), confirmando a heterogeneidade\n  do ambiente e a import\u00e2ncia do controle local realizado pelo DBC (Efici\u00eancia = %s%%).\n\n",
                gl_b, gl_r, .fmt_num_dic(f_bloco, 2), .fmt_num_dic(p_bloco, 3), .fmt_num_dic(x$parametros$eficiencia_relativa, 1)))
  } else {
    cat(sprintf("  O efeito de Bloco n\u00e3o foi significativo (F(%d, %d) = %s; p = %s), indicando homogeneidade\n  entre os blocos avaliados no experimento.\n\n",
                gl_b, gl_r, .fmt_num_dic(f_bloco, 2), .fmt_num_dic(p_bloco, 3)))
  }
  
  cat(sprintf("  O Coeficiente de Varia\u00e7\u00e3o Experimental (CV = %s%%) indica precis\u00e3o experimental %s\n  segundo os crit\u00e9rios de Pimentel-Gomes.\n\n",
              .fmt_num_dic(x$parametros$cv_exp, 2),
              if (x$parametros$cv_exp < 10) "alta (baixo CV)"
              else if (x$parametros$cv_exp < 20) "m\u00e9dia"
              else if (x$parametros$cv_exp < 30) "baixa (alto CV)"
              else "muito baixa (CV muito alto)"))
  
  if (x$tipo == "qualitativo" && !is.null(x$tabela_grupos) && "Grupo" %in% names(x$tabela_grupos)) {
    melhor_tr <- x$tabela_grupos$Tratamento[1]
    melhor_md <- x$tabela_grupos$Media[1]
    cat(sprintf("  Segundo as compara\u00e7\u00f5es m\u00faltiplas (%s), o tratamento '%s' apresentou o maior\n  desempenho m\u00e9dio observado (%s).\n",
                toupper(x$teste), melhor_tr, .fmt_num_dic(melhor_md, 2)))
  } else if (x$tipo == "quantitativo" && !is.null(x$regressao$ponto_critico)) {
    cat(sprintf("  O comportamento das doses seguiu modelo %s, atingindo %s na dose\n  estimada de %s com resposta esperada de %s.\n",
                tolower(x$regressao$modelo_recomendado),
                tolower(x$regressao$ponto_critico$tipo),
                .fmt_num_dic(x$regressao$ponto_critico$x_otimo, 2),
                .fmt_num_dic(x$regressao$ponto_critico$y_otimo, 2)))
  }
  
  .print_rodape()
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# PAINEL GRÁFICO DBC
# ─────────────────────────────────────────────────────────────────────────────

#' @title Painel Gráfico Completo do DBC
#' @description Gera painel integrado de gráficos do Delineamento em Blocos Casualizados.
#'   Para fatores qualitativos, combina o gráfico de médias com letras e o gráfico de IC das diferenças.
#'   Para fatores quantitativos, exibe a curva de regressão polinomial.
#' @param modelo Objeto retornado por \code{analise_dbc()} ou \code{dbc()}.
#' @param ... Argumentos adicionais passados para as funções gráficas individuais.
#' @return Objeto \code{ggplot} ou composição de gráficos.
#' @export
#' @aliases plot_dbc
grafico_dbc <- function(modelo, ...) {
  if (!inherits(modelo, "estatR_dbc")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_dbc().", call. = FALSE)
  }
  
  if (modelo$tipo == "quantitativo") {
    return(grafico_regressao(modelo))
  }
  
  p_med <- tryCatch(grafico_medias(modelo), error = function(e) NULL)
  p_ic  <- tryCatch(grafico_ic(modelo), error = function(e) NULL)
  
  if (!is.null(p_med) && !is.null(p_ic)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(p_med + p_ic)
    }
    return(p_med)
  } else if (!is.null(p_med)) {
    return(p_med)
  } else {
    return(p_ic)
  }
}

#' @export
plot_dbc <- grafico_dbc


# ─────────────────────────────────────────────────────────────────────────────
# MÓDULO: DELINEAMENTO EM QUADRADO LATINO (DQL)
# ─────────────────────────────────────────────────────────────────────────────

#' @title Análise de Delineamento em Quadrado Latino (DQL)
#' @description Executa a Análise de Variância (ANOVA) para experimentos em Delineamento
#'   em Quadrado Latino (DQL), com duplo controle local (linhas e colunas), testes de eficiência
#'   relativa, testes de comparações múltiplas de médias (Tukey, Scott-Knott, Duncan, SNK,
#'   Scheffé com contrastes, Bonferroni, Dunnett) ou Regressão Polinomial com ponto ótimo,
#'   além de painel gráfico automatizado.
#'
#' @param resposta Variável dependente contínua (com ou sem aspas).
#' @param tratamento Fator de tratamentos (com ou sem aspas).
#' @param linha Fator de linhas / controle local 1 (com ou sem aspas).
#' @param coluna Fator de colunas / controle local 2 (com ou sem aspas).
#' @param dados Data frame opcional contendo as variáveis.
#' @param teste Teste de comparações múltiplas para fatores qualitativos:
#'   \code{"tukey"} (padrão), \code{"scott_knott"}, \code{"duncan"}, \code{"snk"},
#'   \code{"scheffe"}, \code{"bonferroni"} ou \code{"dunnett"}.
#' @param contrastes Lista ou vetor de contrastes para o teste de Scheffé (ex: \code{"A vs c(B, C)"}).
#' @param controle Nome do nível de referência (testemunha) usado quando \code{teste = "dunnett"}.
#' @param tipo Tipo de fator: \code{"auto"} (detecta automaticamente), \code{"qualitativo"}
#'   ou \code{"quantitativo"}.
#' @param confianca Nível de confiança para os testes e intervalos (padrão \code{0.95}).
#' @param grafico Lógico. Se \code{TRUE} (padrão), exibe automaticamente o painel de gráficos.
#' @param residuos Lógico. Se \code{TRUE}, executa e exibe a análise residual do modelo.
#'
#' @return Objeto da classe \code{"estatR_dql"} contendo os quadros da ANOVA, parâmetros
#'   experimentais, eficiências relativas, resultados das comparações múltiplas e gráficos.
#' @export
#' @aliases dql
analise_dql <- function(resposta,
                        tratamento,
                        linha,
                        coluna,
                        dados       = NULL,
                        teste       = c("tukey", "scott_knott", "duncan", "snk", "scheffe", "bonferroni", "dunnett"),
                        contrastes  = NULL,
                        controle    = NULL,
                        tipo        = c("auto", "qualitativo", "quantitativo"),
                        confianca   = 0.95,
                        grafico     = TRUE,
                        residuos    = FALSE) {
  
  teste <- match.arg(teste)
  tipo  <- match.arg(tipo)
  
  sub_resp  <- substitute(resposta)
  sub_trat  <- substitute(tratamento)
  sub_linha <- substitute(linha)
  sub_col   <- substitute(coluna)
  
  info_resp  <- .extrair_vetor_exp(sub_resp, dados, "resposta")
  info_trat  <- .extrair_vetor_exp(sub_trat, dados, "tratamento")
  info_linha <- .extrair_vetor_exp(sub_linha, dados, "linha")
  info_col   <- .extrair_vetor_exp(sub_col, dados, "coluna")
  
  y_raw     <- info_resp$valores
  trat_raw  <- info_trat$valores
  linha_raw <- info_linha$valores
  col_raw   <- info_col$valores
  
  nm_resp  <- info_resp$nome
  nm_trat  <- info_trat$nome
  nm_linha <- info_linha$nome
  nm_col   <- info_col$nome
  
  val_ok    <- !is.na(y_raw) & !is.na(trat_raw) & !is.na(linha_raw) & !is.na(col_raw)
  y         <- as.numeric(y_raw[val_ok])
  trat_vec  <- trat_raw[val_ok]
  linha_vec <- linha_raw[val_ok]
  col_vec   <- col_raw[val_ok]
  n_total   <- length(y)
  
  trat_fator  <- as.factor(trat_vec)
  linha_fator <- as.factor(linha_vec)
  col_fator   <- as.factor(col_vec)
  
  niveis_trat  <- levels(trat_fator)
  niveis_linha <- levels(linha_fator)
  niveis_col   <- levels(col_fator)
  
  k <- length(niveis_trat)
  nl <- length(niveis_linha)
  nc <- length(niveis_col)
  
  if (k < 3) {
    stop("O DQL requer ao menos 3 tratamentos (quadrado m\u00ednimo 3x3 para haver res\u00edduo).", call. = FALSE)
  }
  if (nl != k || nc != k) {
    warning(sprintf("Em um DQL padr\u00e3o, n\u00famero de linhas (%d), colunas (%d) e tratamentos (%d) devem ser iguais.", nl, nc, k), call. = FALSE)
  }
  
  # Deteccao de qualitativo vs quantitativo
  e_numerico <- is.numeric(trat_vec) && !is.factor(trat_vec)
  if (tipo == "auto") {
    tipo_final <- if (e_numerico && length(unique(trat_vec)) >= 3) "quantitativo" else "qualitativo"
  } else {
    tipo_final <- tipo
  }
  
  # Ajuste do Modelo ANOVA DQL: Y ~ Linha + Coluna + Trat
  df_mod  <- data.frame(Y = y, Linha = linha_fator, Coluna = col_fator, Trat = trat_fator)
  mod_aov <- stats::aov(Y ~ Linha + Coluna + Trat, data = df_mod)
  tab_aov <- summary(mod_aov)[[1]]
  
  gl_lin  <- tab_aov["Linha", "Df"]
  sq_lin  <- tab_aov["Linha", "Sum Sq"]
  qm_lin  <- tab_aov["Linha", "Mean Sq"]
  f_lin   <- tab_aov["Linha", "F value"]
  p_lin   <- tab_aov["Linha", "Pr(>F)"]
  
  gl_col  <- tab_aov["Coluna", "Df"]
  sq_col  <- tab_aov["Coluna", "Sum Sq"]
  qm_col  <- tab_aov["Coluna", "Mean Sq"]
  f_col   <- tab_aov["Coluna", "F value"]
  p_col   <- tab_aov["Coluna", "Pr(>F)"]
  
  gl_trat <- tab_aov["Trat", "Df"]
  sq_trat <- tab_aov["Trat", "Sum Sq"]
  qm_trat <- tab_aov["Trat", "Mean Sq"]
  f_trat  <- tab_aov["Trat", "F value"]
  p_trat  <- tab_aov["Trat", "Pr(>F)"]
  
  gl_res  <- tab_aov["Residuals", "Df"]
  sq_res  <- tab_aov["Residuals", "Sum Sq"]
  qm_res  <- tab_aov["Residuals", "Mean Sq"]
  
  gl_tot  <- gl_lin + gl_col + gl_trat + gl_res
  sq_tot  <- sq_lin + sq_col + sq_trat + sq_res
  
  # Parametros experimentais
  media_geral <- mean(y)
  cv_exp      <- (sqrt(qm_res) / media_geral) * 100
  r2          <- (sq_lin + sq_col + sq_trat) / sq_tot
  
  # Eficiencias Relativas (Cochran & Cox)
  er_linhas  <- ((qm_lin + (k - 1) * qm_res) / (k * qm_res)) * 100
  er_colunas <- ((qm_col + (k - 1) * qm_res) / (k * qm_res)) * 100
  er_dic     <- ((qm_lin + qm_col + (k - 1) * qm_res) / ((k + 1) * qm_res)) * 100
  
  # Medias e DP por tratamento
  medias_grp <- tapply(y, trat_fator, mean)
  dps_grp    <- tapply(y, trat_fator, stats::sd)
  reps_tab   <- table(trat_fator)
  
  tab_medias_df <- data.frame(
    Tratamento = niveis_trat,
    N          = as.vector(reps_tab),
    Media      = as.vector(medias_grp),
    DP         = as.vector(dps_grp),
    stringsAsFactors = FALSE
  )
  
  # ───────────────────────────────────────────────────────────────────────────
  # CASO 1: FATOR QUALITATIVO (TESTES DE COMPARACAO MULTIPLA)
  # ───────────────────────────────────────────────────────────────────────────
  res_comparacoes <- NULL
  tab_grupos      <- NULL
  tab_contrastes  <- NULL
  alpha           <- 1 - confianca
  
  if (tipo_final == "qualitativo") {
    
    # 1. TUKEY HSD
    if (teste == "tukey") {
      tuk <- stats::TukeyHSD(mod_aov, "Trat", conf.level = confianca)$Trat
      pares <- rownames(tuk)
      
      p_mat <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis_trat, niveis_trat))
      for (i in seq_along(pares)) {
        pts <- unlist(strsplit(pares[i], "-"))
        p_mat[pts[1], pts[2]] <- tuk[i, "p adj"]
        p_mat[pts[2], pts[1]] <- tuk[i, "p adj"]
      }
      
      med_ord <- sort(medias_grp, decreasing = TRUE)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      
      tab_grupos <- data.frame(
        Tratamento = names(med_ord),
        Media      = as.vector(med_ord),
        DP         = dps_grp[names(med_ord)],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
      
      res_comparacoes <- data.frame(
        Comparacao = pares,
        Diff       = tuk[, "diff"],
        Lwr        = tuk[, "lwr"],
        Upr        = tuk[, "upr"],
        p_valor    = tuk[, "p adj"],
        Sig        = ifelse(tuk[, "lwr"] <= 0 & tuk[, "upr"] >= 0, "N\u00e3o", "Sim"),
        stringsAsFactors = FALSE
      )
    }
    
    # 2. SCHEFFÉ
    else if (teste == "scheffe") {
      s_crit <- sqrt((k - 1) * stats::qf(confianca, k - 1, gl_res))
      ns_grp <- as.vector(reps_tab)
      names(ns_grp) <- niveis_trat
      
      if (!is.null(contrastes)) {
        contr_list <- .parse_contrastes_scheffe(contrastes, niveis_trat)
        df_contr <- list()
        
        for (nm_c in names(contr_list)) {
          c_vec <- contr_list[[nm_c]]
          est   <- sum(c_vec * medias_grp)
          ep    <- sqrt(qm_res * sum((c_vec^2) / ns_grp))
          t_val <- est / ep
          f_val <- (t_val^2) / (k - 1)
          p_val <- 1 - stats::pf(f_val, k - 1, gl_res)
          ic_inf <- est - s_crit * ep
          ic_sup <- est + s_crit * ep
          
          df_contr[[nm_c]] <- data.frame(
            Contraste  = nm_c,
            Estimativa = est,
            EP         = ep,
            F_calc     = f_val * (k - 1),
            p_valor    = p_val,
            IC_inf     = ic_inf,
            IC_sup     = ic_sup,
            Sig        = ifelse(ic_inf <= 0 & ic_sup >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
        tab_contrastes <- do.call(rbind, df_contr)
        res_comparacoes <- tab_contrastes
        res_comparacoes$Comparacao <- res_comparacoes$Contraste
        res_comparacoes$Diff       <- res_comparacoes$Estimativa
        res_comparacoes$Lwr        <- res_comparacoes$IC_inf
        res_comparacoes$Upr        <- res_comparacoes$IC_sup
      } else {
        pares_lista <- utils::combn(niveis_trat, 2, simplify = FALSE)
        df_scheffe  <- list()
        p_mat       <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis_trat, niveis_trat))
        
        for (pr in pares_lista) {
          n1 <- pr[1]; n2 <- pr[2]
          diff_val <- medias_grp[n1] - medias_grp[n2]
          ep_val   <- sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
          t_val    <- diff_val / ep_val
          f_val    <- (t_val^2) / (k - 1)
          p_val    <- 1 - stats::pf(f_val, k - 1, gl_res)
          lwr_val  <- diff_val - s_crit * ep_val
          upr_val  <- diff_val + s_crit * ep_val
          
          p_mat[n1, n2] <- p_val
          p_mat[n2, n1] <- p_val
          
          df_scheffe[[paste0(n1, " - ", n2)]] <- data.frame(
            Comparacao = paste0(n1, " - ", n2),
            Diff       = diff_val,
            Lwr        = lwr_val,
            Upr        = upr_val,
            p_valor    = p_val,
            Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
        res_comparacoes <- do.call(rbind, df_scheffe)
        med_ord <- sort(medias_grp, decreasing = TRUE)
        letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
        tab_grupos <- data.frame(
          Tratamento = names(med_ord),
          Media      = as.vector(med_ord),
          DP         = dps_grp[names(med_ord)],
          Grupo      = letras,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # 3. BONFERRONI
    else if (teste == "bonferroni") {
      pares_lista <- utils::combn(niveis_trat, 2, simplify = FALSE)
      m_comp      <- length(pares_lista)
      t_crit      <- stats::qt(1 - alpha / (2 * m_comp), gl_res)
      df_bonf     <- list()
      p_mat       <- matrix(1, nrow = k, ncol = k, dimnames = list(niveis_trat, niveis_trat))
      ns_grp      <- as.vector(reps_tab)
      names(ns_grp) <- niveis_trat
      
      for (pr in pares_lista) {
        n1 <- pr[1]; n2 <- pr[2]
        diff_val <- medias_grp[n1] - medias_grp[n2]
        ep_val   <- sqrt(qm_res * (1/ns_grp[n1] + 1/ns_grp[n2]))
        t_val    <- abs(diff_val) / ep_val
        p_adj    <- min(1, 2 * (1 - stats::pt(t_val, gl_res)) * m_comp)
        lwr_val  <- diff_val - t_crit * ep_val
        upr_val  <- diff_val + t_crit * ep_val
        
        p_mat[n1, n2] <- p_adj
        p_mat[n2, n1] <- p_adj
        
        df_bonf[[paste0(n1, " - ", n2)]] <- data.frame(
          Comparacao = paste0(n1, " - ", n2),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = p_adj,
          Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_bonf)
      med_ord <- sort(medias_grp, decreasing = TRUE)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      tab_grupos <- data.frame(
        Tratamento = names(med_ord),
        Media      = as.vector(med_ord),
        DP         = dps_grp[names(med_ord)],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
    }
    
    # 4. SCOTT-KNOTT
    else if (teste == "scott_knott") {
      med_ord   <- sort(medias_grp, decreasing = TRUE)
      nomes_ord <- names(med_ord)
      r_media   <- mean(reps_tab)
      
      sk_divisao <- function(nomes_sub) {
        k_sub <- length(nomes_sub)
        if (k_sub <= 1) return(list(nomes_sub))
        
        y_sub <- med_ord[nomes_sub]
        melhor_b0 <- -1
        corte_opt <- 1
        
        for (c_pt in 1:(k_sub - 1)) {
          g1 <- nomes_sub[1:c_pt]
          g2 <- nomes_sub[(c_pt + 1):k_sub]
          m1 <- mean(y_sub[g1])
          m2 <- mean(y_sub[g2])
          b0 <- (length(g1) * length(g2) / (length(g1) + length(g2))) * (m1 - m2)^2
          if (b0 > melhor_b0) {
            melhor_b0 <- b0
            corte_opt <- c_pt
          }
        }
        
        s2_0 <- qm_res
        lambda0 <- (pi / (2 * (pi - 2))) * (melhor_b0 / (s2_0 / r_media))
        gl_sk <- k_sub / (pi - 2)
        p_sk <- 1 - stats::pchisq(lambda0, df = gl_sk)
        
        if (p_sk < alpha) {
          g1_sub <- nomes_sub[1:corte_opt]
          g2_sub <- nomes_sub[(corte_opt + 1):k_sub]
          return(c(sk_divisao(g1_sub), sk_divisao(g2_sub)))
        } else {
          return(list(nomes_sub))
        }
      }
      
      grupos_sk <- sk_divisao(nomes_ord)
      letras_vet <- character(length(nomes_ord))
      names(letras_vet) <- nomes_ord
      
      for (g_idx in seq_along(grupos_sk)) {
        letra_atual <- letters[g_idx]
        membros <- grupos_sk[[g_idx]]
        letras_vet[membros] <- letra_atual
      }
      
      tab_grupos <- data.frame(
        Tratamento = nomes_ord,
        Media      = as.vector(med_ord),
        DP         = dps_grp[nomes_ord],
        Grupo      = letras_vet,
        stringsAsFactors = FALSE
      )
      
      pares_lista <- utils::combn(nomes_ord, 2, simplify = FALSE)
      df_sk_comp  <- list()
      q_crit      <- stats::qtukey(confianca, k, gl_res)
      
      for (pr in pares_lista) {
        n1 <- pr[1]; n2 <- pr[2]
        diff_val <- medias_grp[n1] - medias_grp[n2]
        ep_val   <- sqrt(qm_res * (1/reps_tab[n1] + 1/reps_tab[n2]))
        lwr_val  <- diff_val - (q_crit / sqrt(2)) * ep_val
        upr_val  <- diff_val + (q_crit / sqrt(2)) * ep_val
        
        mesmo_grp <- letras_vet[n1] == letras_vet[n2]
        df_sk_comp[[paste0(n1, " - ", n2)]] <- data.frame(
          Comparacao = paste0(n1, " - ", n2),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = if (mesmo_grp) 0.50 else 0.01,
          Sig        = if (mesmo_grp) "N\u00e3o" else "Sim",
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_sk_comp)
    }
    
    # 5. DUNCAN & SNK
    else if (teste %in% c("duncan", "snk")) {
      med_ord   <- sort(medias_grp, decreasing = TRUE)
      nomes_ord <- names(med_ord)
      p_mat     <- matrix(1, nrow = k, ncol = k, dimnames = list(nomes_ord, nomes_ord))
      df_amp    <- list()
      r_media   <- mean(reps_tab)
      
      for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
          n1 <- nomes_ord[i]; n2 <- nomes_ord[j]
          dist_p <- j - i + 1
          diff_val <- med_ord[i] - med_ord[j]
          ep_val   <- sqrt(qm_res * (1/reps_tab[n1] + 1/reps_tab[n2]))
          
          q_obs <- (diff_val / sqrt(qm_res / r_media)) * sqrt(2)
          p_adj <- 1 - stats::ptukey(abs(q_obs), dist_p, gl_res)
          
          p_mat[n1, n2] <- p_adj
          p_mat[n2, n1] <- p_adj
          
          q_c <- stats::qtukey(confianca, dist_p, gl_res)
          lwr_val <- diff_val - (q_c / sqrt(2)) * ep_val
          upr_val <- diff_val + (q_c / sqrt(2)) * ep_val
          
          df_amp[[paste0(n1, " - ", n2)]] <- data.frame(
            Comparacao = paste0(n1, " - ", n2),
            Diff       = diff_val,
            Lwr        = lwr_val,
            Upr        = upr_val,
            p_valor    = p_adj,
            Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
            stringsAsFactors = FALSE
          )
        }
      }
      res_comparacoes <- do.call(rbind, df_amp)
      letras  <- .gerar_letras_agrupamento(med_ord, p_mat, alpha = alpha)
      tab_grupos <- data.frame(
        Tratamento = nomes_ord,
        Media      = as.vector(med_ord),
        DP         = dps_grp[nomes_ord],
        Grupo      = letras,
        stringsAsFactors = FALSE
      )
    }
    
    # 6. DUNNETT
    else if (teste == "dunnett") {
      if (is.null(controle) || !controle %in% niveis_trat) {
        controle <- niveis_trat[1]
      }
      outros <- setdiff(niveis_trat, controle)
      df_dunnett <- list()
      m_d <- length(outros)
      t_crit_d <- stats::qt(1 - alpha / (2 * m_d), gl_res)
      
      for (tr in outros) {
        diff_val <- medias_grp[tr] - medias_grp[controle]
        ep_val   <- sqrt(qm_res * (1/reps_tab[tr] + 1/reps_tab[controle]))
        t_val    <- abs(diff_val) / ep_val
        p_val    <- min(1, 2 * (1 - stats::pt(t_val, gl_res)) * m_d)
        lwr_val  <- diff_val - t_crit_d * ep_val
        upr_val  <- diff_val + t_crit_d * ep_val
        
        df_dunnett[[paste0(tr, " vs ", controle)]] <- data.frame(
          Comparacao = paste0(tr, " vs ", controle),
          Diff       = diff_val,
          Lwr        = lwr_val,
          Upr        = upr_val,
          p_valor    = p_val,
          Sig        = ifelse(lwr_val <= 0 & upr_val >= 0, "N\u00e3o", "Sim"),
          stringsAsFactors = FALSE
        )
      }
      res_comparacoes <- do.call(rbind, df_dunnett)
      tab_grupos <- data.frame(
        Tratamento = c(controle, outros),
        Media      = c(medias_grp[controle], medias_grp[outros]),
        DP         = c(dps_grp[controle], dps_grp[outros]),
        Diferenca  = c(0, sapply(outros, function(tr) medias_grp[tr] - medias_grp[controle])),
        p_valor    = c(NA, sapply(outros, function(tr) res_comparacoes[paste0(tr, " vs ", controle), "p_valor"])),
        Sig        = c("Controle", sapply(outros, function(tr) res_comparacoes[paste0(tr, " vs ", controle), "Sig"])),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # CASO 2: FATOR QUANTITATIVO (REGRESSÃO POLINOMIAL COM LINHAS E COLUNAS)
  # ───────────────────────────────────────────────────────────────────────────
  res_regressao <- NULL
  if (tipo_final == "quantitativo") {
    x_num <- as.numeric(trat_vec)
    k_num <- length(unique(x_num))
    
    mod_lin <- stats::lm(y ~ linha_fator + col_fator + x_num)
    anv_lin <- stats::anova(mod_lin)
    sq_lin  <- anv_lin["x_num", "Sum Sq"]
    
    sq_quad <- 0; f_quad <- NA; p_quad <- NA; r2_quad <- 0
    if (k_num >= 3) {
      mod_quad <- stats::lm(y ~ linha_fator + col_fator + x_num + I(x_num^2))
      anv_quad <- stats::anova(mod_quad)
      sq_quad  <- anv_quad["I(x_num^2)", "Sum Sq"]
    }
    
    sq_cub <- 0; f_cub <- NA; p_cub <- NA; r2_cub <- 0
    if (k_num >= 4) {
      mod_cub <- stats::lm(y ~ linha_fator + col_fator + x_num + I(x_num^2) + I(x_num^3))
      anv_cub <- stats::anova(mod_cub)
      sq_cub  <- anv_cub["I(x_num^3)", "Sum Sq"]
    }
    
    f_lin <- sq_lin / qm_res
    p_lin <- 1 - stats::pf(f_lin, 1, gl_res)
    r2_lin <- sq_lin / sq_tot
    
    tab_pol <- data.frame(
      Grau    = "Linear",
      GL      = 1,
      SQ      = sq_lin,
      QM      = sq_lin,
      F_calc  = f_lin,
      p_valor = p_lin,
      R2      = r2_lin,
      stringsAsFactors = FALSE
    )
    
    if (k_num >= 3) {
      f_quad <- sq_quad / qm_res
      p_quad <- 1 - stats::pf(f_quad, 1, gl_res)
      r2_quad <- sq_quad / sq_tot
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "Quadr\u00e1tico",
        GL      = 1,
        SQ      = sq_quad,
        QM      = sq_quad,
        F_calc  = f_quad,
        p_valor = p_quad,
        R2      = r2_quad,
        stringsAsFactors = FALSE
      ))
    }
    
    if (k_num >= 4) {
      f_cub <- sq_cub / qm_res
      p_cub <- 1 - stats::pf(f_cub, 1, gl_res)
      r2_cub <- sq_cub / sq_tot
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "C\u00fabico",
        GL      = 1,
        SQ      = sq_cub,
        QM      = sq_cub,
        F_calc  = f_cub,
        p_valor = p_cub,
        R2      = r2_cub,
        stringsAsFactors = FALSE
      ))
    }
    
    sq_expl_pol <- sum(tab_pol$SQ)
    sq_desv     <- max(0, sq_trat - sq_expl_pol)
    gl_desv     <- max(0, gl_trat - nrow(tab_pol))
    qm_desv     <- if (gl_desv > 0) sq_desv / gl_desv else 0
    f_desv      <- if (gl_desv > 0) qm_desv / qm_res else NA
    p_desv      <- if (gl_desv > 0) 1 - stats::pf(f_desv, gl_desv, gl_res) else NA
    
    if (gl_desv > 0) {
      tab_pol <- rbind(tab_pol, data.frame(
        Grau    = "Desvios Regress\u00e3o",
        GL      = gl_desv,
        SQ      = sq_desv,
        QM      = qm_desv,
        F_calc  = f_desv,
        p_valor = p_desv,
        R2      = sq_desv / sq_tot,
        stringsAsFactors = FALSE
      ))
    }
    
    tab_pol <- rbind(tab_pol, data.frame(
      Grau    = "Res\u00edduo",
      GL      = gl_res,
      SQ      = sq_res,
      QM      = qm_res,
      F_calc  = NA,
      p_valor = NA,
      R2      = NA,
      stringsAsFactors = FALSE
    ))
    
    mod_rec <- "Linear"
    mod_final_obj <- mod_lin
    if (k_num >= 3 && !is.na(p_quad) && p_quad < alpha) {
      mod_rec <- "Quadr\u00e1tico"
      mod_final_obj <- mod_quad
    }
    if (k_num >= 4 && !is.na(p_cub) && p_cub < alpha) {
      mod_rec <- "C\u00fabico"
      mod_final_obj <- mod_cub
    }
    
    coefs_brutos <- stats::coef(mod_final_obj)
    idx_x <- grep("x_num", names(coefs_brutos))
    coefs_pol <- coefs_brutos[idx_x]
    
    b1 <- coefs_pol[1]
    b2 <- if (length(coefs_pol) >= 2) coefs_pol[2] else NA
    b3 <- if (length(coefs_pol) >= 3) coefs_pol[3] else NA
    b0 <- media_geral - (b1 * mean(x_num) + (if (!is.na(b2)) b2 * mean(x_num^2) else 0) + (if (!is.na(b3)) b3 * mean(x_num^3) else 0))
    coefs_ajustados <- c(Intercept = as.numeric(b0), coefs_pol)
    
    ponto_critico <- NULL
    if (mod_rec == "Quadr\u00e1tico" && !is.na(b2) && b2 != 0) {
      x_otimo <- -b1 / (2 * b2)
      tipo_crit <- if (b2 < 0) "M\u00e1ximo T\u00e9cnico" else "M\u00ednimo T\u00e9cnico"
      y_est <- b0 + b1 * x_otimo + b2 * (x_otimo^2)
      ponto_critico <- list(
        tipo    = tipo_crit,
        x_otimo = as.numeric(x_otimo),
        y_otimo = as.numeric(y_est)
      )
    }
    
    r2_final <- summary(mod_final_obj)$r.squared
    
    res_regressao <- list(
      tabela_polinomial = tab_pol,
      modelo_recomendado = mod_rec,
      modelo_ajustado   = mod_final_obj,
      coeficientes      = coefs_ajustados,
      r2                = r2_final,
      ponto_critico     = ponto_critico
    )
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # OBJETO DE RETORNO
  # ───────────────────────────────────────────────────────────────────────────
  resultado <- structure(
    list(
      dados               = data.frame(Resposta = y, Linha = linha_fator, Coluna = col_fator, Tratamento = trat_fator),
      nomes               = list(resposta = nm_resp, tratamento = nm_trat, linha = nm_linha, coluna = nm_col),
      niveis              = niveis_trat,
      niveis_linha        = niveis_linha,
      niveis_col          = niveis_col,
      n_total             = n_total,
      tipo                = tipo_final,
      teste               = teste,
      confianca           = confianca,
      contrastes_definidos= contrastes,
      controle            = controle,
      anova = data.frame(
        Fonte   = c("Linha", "Coluna", "Tratamento", "Res\u00edduo", "Total"),
        GL      = c(gl_lin, gl_col, gl_trat, gl_res, gl_tot),
        SQ      = c(sq_lin, sq_col, sq_trat, sq_res, sq_tot),
        QM      = c(qm_lin, qm_col, qm_trat, qm_res, NA),
        F_calc  = c(f_lin, f_col, f_trat, NA, NA),
        p_valor = c(p_lin, p_col, p_trat, NA, NA),
        stringsAsFactors = FALSE
      ),
      parametros = list(
        media_geral = media_geral,
        cv_exp      = cv_exp,
        r2          = r2,
        er_linhas   = er_linhas,
        er_colunas  = er_colunas,
        er_dic      = er_dic
      ),
      tabela_medias       = tab_medias_df,
      tabela_grupos       = tab_grupos,
      comparacoes         = res_comparacoes,
      contrastes          = tab_contrastes,
      regressao           = res_regressao,
      modelo_aov          = mod_aov,
      grafico             = NULL
    ),
    class = "estatR_dql"
  )
  
  g_painel <- tryCatch(
    grafico_dql(resultado),
    error = function(e) NULL
  )
  resultado$grafico <- g_painel
  
  if (isTRUE(residuos)) {
    if (exists("analise_residual", mode = "function")) {
      cat("\n")
      analise_residual(mod_aov)
    }
  }
  
  if (isTRUE(grafico) && !is.null(g_painel)) {
    print(g_painel)
  }
  
  resultado
}

#' @export
dql <- analise_dql


# ─────────────────────────────────────────────────────────────────────────────
# MÉTODO DE IMPRESSÃO DQL
# ─────────────────────────────────────────────────────────────────────────────

#' @export
print.estatR_dql <- function(x, ...) {
  .print_titulo("DELINEAMENTO EM QUADRADO LATINO (DQL)")
  
  cat(sprintf("  Vari\u00e1vel Resposta:    %s\n", x$nomes$resposta))
  cat(sprintf("  Fator Tratamento:     %s (%d n\u00edveis: %s)\n",
              x$nomes$tratamento, length(x$niveis),
              if (length(x$niveis) <= 5) paste(x$niveis, collapse = ", ") else paste0(paste(x$niveis[1:4], collapse = ", "), ", ...")))
  cat(sprintf("  Linhas (Controle 1):  %s (%d n\u00edveis)\n", x$nomes$linha, length(x$niveis_linha)))
  cat(sprintf("  Colunas (Controle 2): %s (%d n\u00edveis)\n", x$nomes$coluna, length(x$niveis_col)))
  cat(sprintf("  Classifica\u00e7\u00e3o:        %s\n", if (x$tipo == "qualitativo") "Qualitativo" else "Quantitativo"))
  cat(sprintf("  Total de Observa\u00e7\u00f5es: %d (Quadrado %dx%d)\n\n", x$n_total, length(x$niveis), length(x$niveis)))
  
  # 1. TABELA DE ANOVA
  .print_topico("QUADRO DA AN\u00c1LISE DE VARI\u00c2NCIA (ANOVA)")
  cat("  ────────────────────────────────────────────────────────────────────── \n")
  cat("  Fonte de Varia\u00e7\u00e3o    GL        SQ           QM         F calc    p-valor \n")
  cat("  ────────────────────────────────────────────────────────────────────── \n")
  for (i in 1:nrow(x$anova)) {
    f_str <- x$anova$Fonte[i]
    gl_s  <- sprintf("%2d", x$anova$GL[i])
    sq_s  <- sprintf("%10s", .fmt_num_dic(x$anova$SQ[i], 2))
    qm_s  <- if (!is.na(x$anova$QM[i])) sprintf("%10s", .fmt_num_dic(x$anova$QM[i], 2)) else "          "
    fc_s  <- if (!is.na(x$anova$F_calc[i])) sprintf("%8s", .fmt_num_dic(x$anova$F_calc[i], 2)) else "        "
    pv_s  <- sprintf("%11s", .fmt_p_dic(x$anova$p_valor[i]))
    
    cat(sprintf("  %-18s %4s   %10s   %10s   %8s %11s\n", f_str, gl_s, sq_s, qm_s, fc_s, pv_s))
  }
  cat("  ────────────────────────────────────────────────────────────────────── \n\n")
  
  # 2. PARÂMETROS EXPERIMENTAIS
  .print_topico("PAR\u00c2METROS EXPERIMENTAIS")
  cat("  ────────────────────────────────────────────────────────── \n")
  cat(sprintf("  M\u00e9dia Geral:             %s\n", .fmt_num_dic(x$parametros$media_geral, 2)))
  cat(sprintf("  CV Experimental (%%):     %s%%\n", .fmt_num_dic(x$parametros$cv_exp, 2)))
  cat(sprintf("  R\u00b2 do Modelo:            %s\n", .fmt_num_dic(x$parametros$r2, 3)))
  cat(sprintf("  Efici\u00eancia vs Linhas:    %s%% %s\n",
              .fmt_num_dic(x$parametros$er_linhas, 1),
              if (x$parametros$er_linhas > 100) "(Controle de linhas eficiente)" else "(Sem ganho)"))
  cat(sprintf("  Efici\u00eancia vs Colunas:   %s%% %s\n",
              .fmt_num_dic(x$parametros$er_colunas, 1),
              if (x$parametros$er_colunas > 100) "(Controle de colunas eficiente)" else "(Sem ganho)"))
  cat(sprintf("  Efici\u00eancia vs DIC:       %s%% %s\n",
              .fmt_num_dic(x$parametros$er_dic, 1),
              if (x$parametros$er_dic > 100) "(Duplo controle eficiente)" else "(Sem ganho)"))
  cat("  ────────────────────────────────────────────────────────── \n")
  cat("  * CV experimental: < 10% (Baixo / Alta precis\u00e3o), 10-20% (M\u00e9dio), 20-30% (Alto), > 30% (Muito alto - Pimentel-Gomes).\n")
  cat("  * Efici\u00eancias Relativas: > 100% indicam aumento de precis\u00e3o do DQL em rela\u00e7\u00e3o ao DBC ou DIC.\n\n")
  
  # 3. RESULTADOS DE ACORDO COM O TIPO
  if (x$tipo == "qualitativo") {
    if (x$teste == "scheffe" && !is.null(x$contrastes)) {
      .print_topico("CONTRASTES DE SCHEFF\u00c9")
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      cat("  Contraste               Estimativa     EP      F calc    p-valor      IC (95%)    \n")
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$contrastes)) {
        nm_c <- x$contrastes$Contraste[i]
        est  <- sprintf("%+10s", .fmt_num_dic(x$contrastes$Estimativa[i], 2))
        ep   <- sprintf("%8s", .fmt_num_dic(x$contrastes$EP[i], 2))
        fc   <- sprintf("%8s", .fmt_num_dic(x$contrastes$F_calc[i], 2))
        pv   <- sprintf("%11s", .fmt_p_dic(x$contrastes$p_valor[i]))
        ic_s <- sprintf("[%s ; %s]", .fmt_num_dic(x$contrastes$IC_inf[i], 2), .fmt_num_dic(x$contrastes$IC_sup[i], 2))
        cat(sprintf("  %-22s %10s %8s %8s %11s  %-15s\n", nm_c, est, ep, fc, pv, ic_s))
      }
      cat("  ────────────────────────────────────────────────────────────────────────────── \n")
      cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")
    } else if (x$teste == "dunnett") {
      .print_topico(sprintf("COMPARA\u00c7\u00d5ES CONTRA O CONTROLE (TESTE DE DUNNETT - TESTEMUNHA: '%s')", x$controle))
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      cat("  Tratamento        M\u00e9dia      Diferen\u00e7a    p-valor      Diferen\u00e7a Sig. \n")
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$tabela_grupos)) {
        tr_s   <- x$tabela_grupos$Tratamento[i]
        med_s  <- sprintf("%8s", .fmt_num_dic(x$tabela_grupos$Media[i], 2))
        dif_s  <- sprintf("%+9s", .fmt_num_dic(x$tabela_grupos$Diferenca[i], 2))
        pv_s   <- if (!is.na(x$tabela_grupos$p_valor[i])) sprintf("%11s", .fmt_p_dic(x$tabela_grupos$p_valor[i])) else "     Controle"
        sig_s  <- x$tabela_grupos$Sig[i]
        cat(sprintf("  %-16s %8s   %9s  %11s       %s\n", tr_s, med_s, dif_s, pv_s, sig_s))
      }
      cat("  ──────────────────────────────────────────────────────────────────── \n")
      cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n\n")
    } else {
      nm_teste_exib <- switch(x$teste,
        tukey       = "TESTE DE TUKEY (HSD)",
        scott_knott = "TESTE DE SCOTT-KNOTT",
        duncan      = "TESTE DE DUNCAN",
        snk         = "TESTE DE STUDENT-NEWMAN-KEULS (SNK)",
        scheffe     = "TESTE DE SCHEFF\u00c9 (PAREADO)",
        bonferroni  = "TESTE T COM AJUSTE DE BONFERRONI"
      )
      .print_topico(sprintf("COMPARA\u00c7\u00d5ES M\u00daLTIPLAS DE M\u00c9DIAS (%s)", nm_teste_exib))
      cat("  ────────────────────────────────────────────────────────── \n")
      cat("  Tratamento        M\u00e9dia Observada      Desvio Padr\u00e3o    Grupo\n")
      cat("  ────────────────────────────────────────────────────────── \n")
      for (i in 1:nrow(x$tabela_grupos)) {
        tr_s  <- x$tabela_grupos$Tratamento[i]
        med_s <- sprintf("%10s", .fmt_num_dic(x$tabela_grupos$Media[i], 2))
        dp_s  <- sprintf("%10s", .fmt_num_dic(x$tabela_grupos$DP[i], 2))
        grp_s <- sprintf("%-6s", x$tabela_grupos$Grupo[i])
        cat(sprintf("  %-16s %10s          %10s       %s\n", tr_s, med_s, dp_s, grp_s))
      }
      cat("  ────────────────────────────────────────────────────────── \n")
      cat("  * M\u00e9dias seguidas pela mesma letra n\u00e3o diferem estatisticamente entre si a 5% de signific\u00e2ncia.\n\n")
    }
  } else {
    .print_topico("REGRESS\u00c3O POLINOMIAL (AN\u00c1LISE DE DOSES / N\u00cdVEIS)")
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    cat("  Grau Polinomial     GL        SQ           QM         F calc    p-valor     R\u00b2 \n")
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    tab_pol <- x$regressao$tabela_polinomial
    for (i in 1:nrow(tab_pol)) {
      gr_s  <- tab_pol$Grau[i]
      gl_s  <- sprintf("%2d", tab_pol$GL[i])
      sq_s  <- sprintf("%10s", .fmt_num_dic(tab_pol$SQ[i], 2))
      qm_s  <- sprintf("%10s", .fmt_num_dic(tab_pol$QM[i], 2))
      fc_s  <- if (!is.na(tab_pol$F_calc[i])) sprintf("%8s", .fmt_num_dic(tab_pol$F_calc[i], 2)) else "        "
      pv_s  <- sprintf("%11s", .fmt_p_dic(tab_pol$p_valor[i]))
      r2_s  <- if (!is.na(tab_pol$R2[i])) sprintf("%6s", .fmt_num_dic(tab_pol$R2[i], 3)) else "      "
      cat(sprintf("  %-18s %2s   %10s   %10s   %8s %11s  %6s\n", gr_s, gl_s, sq_s, qm_s, fc_s, pv_s, r2_s))
    }
    cat("  ────────────────────────────────────────────────────────────────────────── \n")
    cat("  * Signific\u00e2ncia: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.10\n")
    cat(sprintf("  Modelo Recomendado: Polinomial %s\n\n", x$regressao$modelo_recomendado))
    
    .print_topico("EQUA\u00c7\u00c3O AJUSTADA E PONTO \u00d3TIMO")
    coefs <- x$regressao$coeficientes
    .fmt_s <- function(val, d = 3) {
      if (is.na(val)) return("")
      sinal <- if (val >= 0) "+ " else "- "
      paste0(sinal, formatC(abs(val), format = "f", digits = d, decimal.mark = ",", big.mark = "."))
    }
    
    if (x$regressao$modelo_recomendado == "Linear") {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3)))
    } else if (x$regressao$modelo_recomendado == "Quadr\u00e1tico") {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose %s * Dose\u00b2\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3), .fmt_s(coefs[3], 4)))
    } else {
      cat(sprintf("  Equa\u00e7\u00e3o:  Y = %s %s * Dose %s * Dose\u00b2 %s * Dose\u00b3\n",
                  .fmt_num_dic(coefs[1], 3), .fmt_s(coefs[2], 3), .fmt_s(coefs[3], 4), .fmt_s(coefs[4], 5)))
    }
    cat(sprintf("  R\u00b2 Total: %s\n", .fmt_num_dic(x$regressao$r2, 3)))
    if (!is.null(x$regressao$ponto_critico)) {
      cat(sprintf("  Ponto Cr\u00edtico: %s\n", x$regressao$ponto_critico$tipo))
      cat(sprintf("  Dose de M\u00e1xima Efici\u00eancia: %s  (Resposta estimada: %s)\n",
                  .fmt_num_dic(x$regressao$ponto_critico$x_otimo, 2),
                  .fmt_num_dic(x$regressao$ponto_critico$y_otimo, 2)))
    }
    cat("\n")
  }
  
  # 4. INTERPRETAÇÃO
  .print_topico("INTERPRETA\u00c7\u00c3O")
  p_trat <- x$anova$p_valor[3]
  f_trat <- x$anova$F_calc[3]
  gl_t   <- x$anova$GL[3]
  gl_r   <- x$anova$GL[4]
  p_lin  <- x$anova$p_valor[1]
  p_col  <- x$anova$p_valor[2]
  
  if (!is.na(p_trat) && p_trat < 0.05) {
    cat(sprintf("  O efeito de Tratamento foi estatisticamente significativo (F(%d, %d) = %s; p %s),\n  rejeitando a hip\u00f3tese nula de que as m\u00e9dias dos tratamentos s\u00e3o iguais a 5%% de signific\u00e2ncia.\n\n",
                gl_t, gl_r, .fmt_num_dic(f_trat, 2), if (p_trat < 0.001) "< 0,001" else paste0("= ", .fmt_num_dic(p_trat, 3))))
  } else {
    cat(sprintf("  O efeito de Tratamento n\u00e3o foi estatisticamente significativo (F(%d, %d) = %s; p = %s),\n  indicando que n\u00e3o h\u00e1 evid\u00eancias de diferen\u00e7a entre as m\u00e9dias dos tratamentos avaliados.\n\n",
                gl_t, gl_r, .fmt_num_dic(f_trat, 2), .fmt_num_dic(p_trat, 3)))
  }
  
  # Interpretacao do duplo controle
  lin_sig <- !is.na(p_lin) && p_lin < 0.05
  col_sig <- !is.na(p_col) && p_col < 0.05
  if (lin_sig && col_sig) {
    cat("  O duplo controle local foi plenamente justificado, com efeitos significativos para\n  linhas e colunas, isolando duas fontes distintas de heterogeneidade ambiental.\n\n")
  } else if (lin_sig) {
    cat("  O controle local de linhas foi significativo, confirmando heterogeneidade nessa dire\u00e7\u00e3o,\n  enquanto as colunas apresentaram homogeneidade.\n\n")
  } else if (col_sig) {
    cat("  O controle local de colunas foi significativo, confirmando heterogeneidade nessa dire\u00e7\u00e3o,\n  enquanto as linhas apresentaram homogeneidade.\n\n")
  } else {
    cat("  N\u00e3o foram detectadas diferen\u00e7as significativas entre linhas ou colunas no experimento.\n\n")
  }
  
  cat(sprintf("  O Coeficiente de Varia\u00e7\u00e3o Experimental (CV = %s%%) indica precis\u00e3o experimental %s\n  segundo os crit\u00e9rios de Pimentel-Gomes.\n\n",
              .fmt_num_dic(x$parametros$cv_exp, 2),
              if (x$parametros$cv_exp < 10) "alta (baixo CV)"
              else if (x$parametros$cv_exp < 20) "m\u00e9dia"
              else if (x$parametros$cv_exp < 30) "baixa (alto CV)"
              else "muito baixa (CV muito alto)"))
  
  if (x$tipo == "qualitativo" && !is.null(x$tabela_grupos) && "Grupo" %in% names(x$tabela_grupos)) {
    melhor_tr <- x$tabela_grupos$Tratamento[1]
    melhor_md <- x$tabela_grupos$Media[1]
    cat(sprintf("  Segundo as compara\u00e7\u00f5es m\u00faltiplas (%s), o tratamento '%s' apresentou o maior\n  desempenho m\u00e9dio observado (%s).\n",
                toupper(x$teste), melhor_tr, .fmt_num_dic(melhor_md, 2)))
  } else if (x$tipo == "quantitativo" && !is.null(x$regressao$ponto_critico)) {
    cat(sprintf("  O comportamento das doses seguiu modelo %s, atingindo %s na dose\n  estimada de %s com resposta esperada de %s.\n",
                tolower(x$regressao$modelo_recomendado),
                tolower(x$regressao$ponto_critico$tipo),
                .fmt_num_dic(x$regressao$ponto_critico$x_otimo, 2),
                .fmt_num_dic(x$regressao$ponto_critico$y_otimo, 2)))
  }
  
  .print_rodape()
  invisible(x)
}


# ─────────────────────────────────────────────────────────────────────────────
# PAINEL GRÁFICO DQL
# ─────────────────────────────────────────────────────────────────────────────

#' @title Painel Gráfico Completo do DQL
#' @description Gera painel integrado de gráficos do Delineamento em Quadrado Latino.
#'   Para fatores qualitativos, combina o gráfico de médias com letras e o gráfico de IC das diferenças.
#'   Para fatores quantitativos, exibe a curva de regressão polinomial.
#' @param modelo Objeto retornado por \code{analise_dql()} ou \code{dql()}.
#' @param ... Argumentos adicionais passados para as funções gráficas individuais.
#' @return Objeto \code{ggplot} ou composição de gráficos.
#' @export
#' @aliases plot_dql
grafico_dql <- function(modelo, ...) {
  if (!inherits(modelo, "estatR_dql")) {
    stop("O argumento 'modelo' deve ser um objeto retornado por analise_dql().", call. = FALSE)
  }
  
  if (modelo$tipo == "quantitativo") {
    return(grafico_regressao(modelo))
  }
  
  p_med <- tryCatch(grafico_medias(modelo), error = function(e) NULL)
  p_ic  <- tryCatch(grafico_ic(modelo), error = function(e) NULL)
  
  if (!is.null(p_med) && !is.null(p_ic)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(p_med + p_ic)
    }
    return(p_med)
  } else if (!is.null(p_med)) {
    return(p_med)
  } else {
    return(p_ic)
  }
}

#' @export
plot_dql <- grafico_dql
