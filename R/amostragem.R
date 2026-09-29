#' CÃ¡lculo de Tamanho Amostral Ã“timo
#'
#' Calcula o tamanho mÃ­nimo de amostra necessÃ¡rio para estimar uma proporÃ§Ã£o
#' ou uma mÃ©dia com nÃ­vel de confianÃ§a e margem de erro especificados.
#'
#' @param populacao Tamanho da populaÃ§Ã£o (N). Se NULL, considera populaÃ§Ã£o infinita.
#' @param erro Margem de erro mÃ¡xima permitida (ex: 0.05 para 5% ou valor na mesma unidade da mÃ©dia).
#' @param confianca NÃ­vel de confianÃ§a da estimativa (padrÃ£o 0.95 para 95%).
#' @param proporcao ProporÃ§Ã£o esperada do evento (padrÃ£o 0.5 para cenÃ¡rio conservador de variÃ¢ncia mÃ¡xima).
#' @param desvio_padrao Desvio-padrÃ£o populacional ou preliminar. Se informado, calcula para mÃ©dia em vez de proporÃ§Ã£o.
#' @return Retorna o nÃºmero inteiro do tamanho da amostra (n) invisivelmente.
#' @export
#' @examples
#' # Para proporÃ§Ã£o com populaÃ§Ã£o de 50.000, 3% de erro e 95% de confianÃ§a:
#' tamanho_amostra(populacao = 50000, erro = 0.03, confianca = 0.95)
#'
#' # Para mÃ©dia com desvio-padrÃ£o de 15, erro de 2 unidades e confianÃ§a de 95%:
#' tamanho_amostra(desvio_padrao = 15, erro = 2, confianca = 0.95)
tamanho_amostra <- function(populacao = NULL, erro = 0.05, confianca = 0.95,
                            proporcao = 0.5, desvio_padrao = NULL) {
  if (erro <= 0) stop("A margem de erro deve ser maior que zero.")
  if (confianca <= 0 || confianca >= 1) stop("O nÃ­vel de confianÃ§a deve estar entre 0 e 1 (ex: 0.95).")
  
  alfa <- 1 - confianca
  z <- qnorm(1 - alfa / 2)
  
  se_media <- !is.null(desvio_padrao)
  
  if (se_media) {
    if (desvio_padrao <= 0) stop("O desvio-padrÃ£o deve ser maior que zero.")
    n0 <- (z^2 * desvio_padrao^2) / (erro^2)
    tipo_param <- sprintf("MÃ©dia (DP = %s)", .fmt_num(desvio_padrao, 2))
    var_info <- sprintf("Margem de Erro:          %s unidades", .fmt_num(erro, 2))
  } else {
    if (proporcao <= 0 || proporcao >= 1) stop("A proporÃ§Ã£o esperada deve estar entre 0 e 1 (ex: 0.5).")
    n0 <- (z^2 * proporcao * (1 - proporcao)) / (erro^2)
    tipo_conservador <- if (proporcao == 0.5) " (CenÃ¡rio Conservador)" else ""
    tipo_param <- sprintf("ProporÃ§Ã£o (p = %s%s)", .fmt_pct(proporcao, 1), tipo_conservador)
    var_info <- sprintf("Margem de Erro:          %s", .fmt_pct(erro, 1))
  }
  
  if (!is.null(populacao)) {
    if (populacao <= 1) stop("O tamanho da populaÃ§Ã£o deve ser maior que 1.")
    n_final <- (populacao * n0) / (populacao + n0 - 1)
    pop_str <- sprintf("%s (Finita)", format(populacao, big.mark = ".", decimal.mark = ","))
    frac_amostral <- (n_final / populacao)
    frac_str <- sprintf("  FraÃ§Ã£o Amostral:             %s da populaÃ§Ã£o\n", .fmt_pct(frac_amostral, 2))
  } else {
    n_final <- n0
    pop_str <- "Infinita (ou nÃ£o informada)"
    frac_str <- ""
  }
  
  n_otimo <- ceiling(n_final)
  
  .print_titulo("DIMENSIONAMENTO DE AMOSTRA (TAMANHO Ã“TIMO)")
  cat("ParÃ¢metros do Estudo:\n")
  cat(sprintf("  â€¢ Tipo de ParÃ¢metro:       %s\n", tipo_param))
  cat(sprintf("  â€¢ NÃ­vel de ConfianÃ§a:      %s (Z = %s)\n", .fmt_pct(confianca, 1), .fmt_num(z, 2)))
  cat(sprintf("  â€¢ %s\n", var_info))
  cat(sprintf("  â€¢ PopulaÃ§Ã£o (N):           %s\n\n", pop_str))
  
  cat("\033[1mâ–¶ Resultado do Dimensionamento\033[0m\n")
  cat(sprintf("  Tamanho Amostral MÃ­nimo (n): %s observaÃ§Ãµes\n", format(n_otimo, big.mark = ".", decimal.mark = ",")))
  cat(frac_str)
  
  .print_rodape()
  
  invisible(n_otimo)
}


#' Amostragem AleatÃ³ria Simples (AAS)
#'
#' Seleciona uma amostra aleatÃ³ria simples de um conjunto de dados ou vetor,
#' com ou sem reposiÃ§Ã£o.
#'
#' @param dados Um data.frame ou vetor representando a populaÃ§Ã£o.
#' @param n Quantidade exata de observaÃ§Ãµes a sortear.
#' @param proporcao ProporÃ§Ã£o da populaÃ§Ã£o a sortear (ex: 0.20 para 20%).
#' @param reposicao LÃ³gico. Se TRUE, permite sortear o mesmo elemento mais de uma vez.
#' @param semente NÃºmero inteiro para semente aleatÃ³ria (garante reprodutibilidade).
#' @return Retorna o subconjunto de dados sorteado.
#' @export
#' @examples
#' amostra_aleatoria(iris, n = 30)
#' amostra_aleatoria(iris, proporcao = 0.2, semente = 123)
amostra_aleatoria <- function(dados, n = NULL, proporcao = NULL, reposicao = FALSE, semente = NULL) {
  if (!is.null(semente)) set.seed(semente)
  
  e_df <- is.data.frame(dados)
  N <- if (e_df) nrow(dados) else length(dados)
  
  if (is.null(n) && is.null(proporcao)) {
    stop("Informe o tamanho desejado da amostra ('n') ou a 'proporcao' a sortear.")
  }
  
  if (!is.null(proporcao)) {
    if (proporcao <= 0 || proporcao > 1) stop("A 'proporcao' deve estar entre 0 e 1 (ex: 0.2 para 20%).")
    n_sorteio <- round(N * proporcao)
  } else {
    n_sorteio <- n
  }
  
  if (!reposicao && n_sorteio > N) {
    stop("O tamanho da amostra (n) nÃ£o pode ser maior que a populaÃ§Ã£o em sorteios sem reposiÃ§Ã£o.")
  }
  
  indices <- sample(1:N, size = n_sorteio, replace = reposicao)
  amostra <- if (e_df) dados[indices, , drop = FALSE] else dados[indices]
  
  pct_amostrada <- (n_sorteio / N) * 100
  tipo_rep <- if (reposicao) "Com reposiÃ§Ã£o" else "Sem reposiÃ§Ã£o"
  
  .print_titulo("AMOSTRAGEM ALEATÃ“RIA SIMPLES")
  cat(sprintf("PopulaÃ§Ã£o (N):       %s observaÃ§Ãµes\n", format(N, big.mark = ".", decimal.mark = ",")))
  cat(sprintf("Amostra Sorteada (n):%s observaÃ§Ãµes (%s%% da populaÃ§Ã£o)\n",
              format(n_sorteio, big.mark = ".", decimal.mark = ","),
              .fmt_num(pct_amostrada, 1)))
  cat(sprintf("Tipo de Sorteio:     %s\n", tipo_rep))
  .print_rodape()
  
  invisible(amostra)
}


#' Amostragem SistemÃ¡tica
#'
#' Seleciona uma amostra sistemÃ¡tica a cada salto de k elementos a partir de um
#' ponto inicial sorteado aleatoriamente.
#'
#' @param dados Um data.frame ou vetor representando a populaÃ§Ã£o.
#' @param n Quantidade desejada de elementos na amostra.
#' @param salto Intervalo de salto (k). Se informado, seleciona 1 elemento a cada k.
#' @param semente NÃºmero inteiro para semente aleatÃ³ria.
#' @return Retorna o subconjunto de dados sorteado.
#' @export
#' @examples
#' amostra_sistematica(iris, n = 30)
#' amostra_sistematica(iris, salto = 5, semente = 42)
amostra_sistematica <- function(dados, n = NULL, salto = NULL, semente = NULL) {
  if (!is.null(semente)) set.seed(semente)
  
  e_df <- is.data.frame(dados)
  N <- if (e_df) nrow(dados) else length(dados)
  
  if (is.null(n) && is.null(salto)) {
    stop("Informe o tamanho desejado da amostra ('n') ou o intervalo de 'salto' (k).")
  }
  
  if (!is.null(salto)) {
    if (salto <= 0 || salto > N) stop("O 'salto' deve ser maior que zero e menor que o tamanho da populaÃ§Ã£o.")
    k <- round(salto)
    n_calculado <- floor(N / k)
  } else {
    if (n <= 0 || n > N) stop("O tamanho da amostra (n) deve ser maior que zero e menor ou igual Ã  populaÃ§Ã£o.")
    k <- floor(N / n)
    if (k < 1) k <- 1
    n_calculado <- n
  }
  
  ponto_partida <- sample(1:k, 1)
  indices <- seq(ponto_partida, N, by = k)
  
  amostra <- if (e_df) dados[indices, , drop = FALSE] else dados[indices]
  n_real <- length(indices)
  pct_amostrada <- (n_real / N) * 100
  
  .print_titulo("AMOSTRAGEM SISTEMÃTICA")
  cat(sprintf("PopulaÃ§Ã£o (N):       %s observaÃ§Ãµes\n", format(N, big.mark = ".", decimal.mark = ",")))
  cat(sprintf("Amostra Sorteada (n):%s observaÃ§Ãµes (%s%% da populaÃ§Ã£o)\n",
              format(n_real, big.mark = ".", decimal.mark = ","),
              .fmt_num(pct_amostrada, 1)))
  cat(sprintf("Salto (Intervalo k): A cada %d elementos\n", k))
  cat(sprintf("Ponto de Partida:    Elemento %d (sorteado aleatoriamente entre 1 e %d)\n", ponto_partida, k))
  .print_rodape()
  
  invisible(amostra)
}


#' Amostragem Estratificada
#'
#' Seleciona uma amostra estratificada de um data.frame,
#' com suporte a trÃªs tipos de alocaÃ§Ã£o.
#'
#' @param dados Um data.frame representando a populaÃ§Ã£o.
#' @param estrato Nome da coluna do estrato (sem aspas ou como string).
#' @param n Quantidade total de observaÃ§Ãµes desejadas na amostra.
#' @param proporcao ProporÃ§Ã£o da populaÃ§Ã£o total a sortear (ex: 0.20 para 20%).
#' @param alocacao Tipo de alocaÃ§Ã£o: "proporcional" (padrÃ£o), "uniforme" ou "otima" (Neyman).
#' @param variavel Nome da coluna numÃ©rica para cÃ¡lculo do desvio-padrÃ£o por estrato (obrigatÃ³rio para alocaÃ§Ã£o Ã³tima).
#' @param semente NÃºmero inteiro para semente aleatÃ³ria.
#' @return Retorna o data.frame com a amostra estratificada selecionada.
#' @export
#' @examples
#' amostra_estratificada(iris, estrato = Species, n = 30)
#' amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "uniforme")
#' amostra_estratificada(iris, estrato = Species, n = 30, alocacao = "otima", variavel = Sepal.Length)
amostra_estratificada <- function(dados, estrato, n = NULL, proporcao = NULL, 
                                   alocacao = c("proporcional", "uniforme", "otima"),
                                   variavel = NULL, semente = NULL) {
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data.frame.")
  if (!is.null(semente)) set.seed(semente)
  
  alocacao <- match.arg(alocacao)
  
  N <- nrow(dados)
  
  # Captura o estrato via NSE ou string
  estrato_sub <- substitute(estrato)
  estrato_nome <- deparse(estrato_sub)
  estrato_nome <- sub(".*\\$", "", estrato_nome)
  
  if (estrato_nome %in% names(dados)) {
    estrato_vec <- dados[[estrato_nome]]
    nome_exibicao <- estrato_nome
  } else {
    estrato_val <- tryCatch(eval(estrato_sub, parent.frame()), error = function(e) NULL)
    if (is.character(estrato_val) && length(estrato_val) == 1 && estrato_val %in% names(dados)) {
      estrato_vec <- dados[[estrato_val]]
      nome_exibicao <- estrato_val
    } else if (length(estrato_val) == N) {
      estrato_vec <- estrato_val
      nome_exibicao <- estrato_nome
    } else {
      stop(sprintf("A coluna de estrato '%s' nÃ£o foi encontrada no banco de dados.", estrato_nome))
    }
  }
  
  # Captura variÃ¡vel para alocaÃ§Ã£o Ã³tima (NSE)
  if (alocacao == "otima") {
    var_sub <- substitute(variavel)
    if (is.null(var_sub) || identical(var_sub, quote(expr=))) {
      stop("Para alocaÃ§Ã£o Ã³tima (Neyman), informe o parÃ¢metro 'variavel' (coluna numÃ©rica para o cÃ¡lculo do desvio-padrÃ£o).")
    }
    var_nome <- deparse(var_sub)
    var_nome <- sub(".*\\$", "", var_nome)
    
    if (var_nome %in% names(dados)) {
      var_vec <- dados[[var_nome]]
    } else {
      var_val <- tryCatch(eval(var_sub, parent.frame()), error = function(e) NULL)
      if (is.character(var_val) && length(var_val) == 1 && var_val %in% names(dados)) {
        var_vec <- dados[[var_val]]
        var_nome <- var_val
      } else {
        stop(sprintf("A coluna '%s' nÃ£o foi encontrada no banco de dados.", var_nome))
      }
    }
    if (!is.numeric(var_vec)) stop(sprintf("A coluna '%s' deve ser numÃ©rica para alocaÃ§Ã£o Ã³tima.", var_nome))
  }
  
  if (is.null(n) && is.null(proporcao)) {
    stop("Informe o tamanho total desejado da amostra ('n') ou a 'proporcao' a sortear.")
  }
  
  if (!is.null(proporcao)) {
    if (proporcao <= 0 || proporcao > 1) stop("A 'proporcao' deve estar entre 0 e 1.")
    n_total <- round(N * proporcao)
  } else {
    n_total <- n
  }
  
  fator_estrato <- as.factor(estrato_vec)
  niveis <- levels(fator_estrato)
  L <- length(niveis)
  
  tab_pop <- table(fator_estrato)
  indices_selecionados <- c()
  tab_amostra_contagem <- integer(L)
  names(tab_amostra_contagem) <- niveis
  
  # Calcula Nh e Sh por estrato (para Ã³tima)
  Nh <- as.numeric(tab_pop)
  
  if (alocacao == "otima") {
    Sh <- numeric(L)
    for (i in seq_along(niveis)) {
      idx_h <- which(fator_estrato == niveis[i])
      Sh[i] <- sd(var_vec[idx_h])
      if (is.na(Sh[i])) Sh[i] <- 0
    }
    NhSh <- Nh * Sh
    soma_NhSh <- sum(NhSh)
    if (soma_NhSh == 0) {
      warning("Todos os estratos possuem desvio-padrÃ£o zero. Usando alocaÃ§Ã£o proporcional.")
      alocacao <- "proporcional"
    }
  }
  
  for (i in seq_along(niveis)) {
    niv <- niveis[i]
    idx_estrato <- which(fator_estrato == niv)
    n_h_pop <- length(idx_estrato)
    
    if (alocacao == "proporcional") {
      n_h_amostra <- round(n_total * (n_h_pop / N))
    } else if (alocacao == "uniforme") {
      n_h_amostra <- round(n_total / L)
    } else if (alocacao == "otima") {
      n_h_amostra <- round(n_total * (NhSh[i] / soma_NhSh))
    }
    
    if (n_h_amostra < 1 && n_h_pop > 0) n_h_amostra <- 1
    if (n_h_amostra > n_h_pop) n_h_amostra <- n_h_pop
    
    sorteados_h <- sample(idx_estrato, size = n_h_amostra, replace = FALSE)
    indices_selecionados <- c(indices_selecionados, sorteados_h)
    tab_amostra_contagem[niv] <- n_h_amostra
  }
  
  amostra_final <- dados[indices_selecionados, , drop = FALSE]
  n_real_total <- length(indices_selecionados)
  
  # RÃ³tulo da alocaÃ§Ã£o
  rotulo_alocacao <- switch(alocacao,
                             proporcional = "Proporcional",
                             uniforme = "Uniforme (igual)",
                             otima = sprintf("Ã“tima de Neyman (var: %s)", var_nome))
  
  # ImpressÃ£o do painel e tabela comparativa
  tit <- sprintf("AMOSTRAGEM ESTRATIFICADA (Por: %s)", nome_exibicao)
  .print_titulo(tit)
  
  cat(sprintf("PopulaÃ§Ã£o (N): %s | Amostra Sorteada (n): %s (%s%%)\n",
              format(N, big.mark = ".", decimal.mark = ","),
              format(n_real_total, big.mark = ".", decimal.mark = ","),
              .fmt_num((n_real_total / N) * 100, 1)))
  cat(sprintf("AlocaÃ§Ã£o:      %s\n\n", rotulo_alocacao))
  
  w_est <- max(14, max(nchar(niveis)), nchar("Estrato"), nchar("Total")) + 2
  w_col <- 15
  
  col_titulos <- c(
    .pad_string("Estrato", w_est, "left"),
    .pad_string("PopulaÃ§Ã£o (%)", w_col, "center"),
    .pad_string("Amostra (N)", w_col, "center"),
    .pad_string("Amostra (%)", w_col, "center")
  )
  hdr <- paste(col_titulos, collapse = " ")
  cat("  ", hdr, "\n", sep = "")
  cat("  ", paste(rep("â”€", nchar(hdr)), collapse = ""), "\n", sep = "")
  
  for (i in seq_along(niveis)) {
    niv <- niveis[i]
    pop_n <- as.numeric(tab_pop[niv])
    pop_pct <- pop_n / N
    am_n <- as.numeric(tab_amostra_contagem[niv])
    am_pct <- if (n_real_total > 0) am_n / n_real_total else 0
    
    str_pop <- sprintf("%d (%s)", pop_n, .fmt_pct(pop_pct, 1))
    
    col_valores <- c(
      .pad_string(niv, w_est, "left"),
      .pad_string(str_pop, w_col, "center"),
      .pad_string(as.character(am_n), w_col, "center"),
      .pad_string(.fmt_pct(am_pct, 1), w_col, "center")
    )
    cat("  ", paste(col_valores, collapse = " "), "\n", sep = "")
  }
  
  cat("  ", paste(rep("â”€", nchar(hdr)), collapse = ""), "\n", sep = "")
  col_totais <- c(
    .pad_string("Total", w_est, "left"),
    .pad_string(sprintf("%d (100,0%%)", N), w_col, "center"),
    .pad_string(as.character(n_real_total), w_col, "center"),
    .pad_string("100,0%", w_col, "center")
  )
  cat("  ", paste(col_totais, collapse = " "), "\n", sep = "")
  
  nota <- switch(alocacao,
                  proporcional = "As proporÃ§Ãµes populacionais de cada estrato foram preservadas.",
                  uniforme = "Cada estrato recebeu o mesmo nÃºmero de observaÃ§Ãµes, independente do seu tamanho.",
                  otima = sprintf("A alocaÃ§Ã£o priorizou estratos com maior variabilidade em '%s' (Neyman, 1934).", var_nome))
  cat(sprintf("\n  * %s\n", nota))
  .print_rodape()
  
  invisible(amostra_final)
}


#' Amostragem por Conglomerados (Clusters)
#'
#' Seleciona aleatoriamente conglomerados (grupos inteiros) de um data.frame,
#' retornando todos os indivÃ­duos pertencentes aos conglomerados sorteados.
#'
#' @param dados Um data.frame representando a populaÃ§Ã£o.
#' @param conglomerado Nome da coluna do conglomerado (sem aspas ou como string).
#' @param n_conglomerados Quantidade de conglomerados inteiros a sortear.
#' @param proporcao ProporÃ§Ã£o dos conglomerados a sortear (ex: 0.20 para 20%).
#' @param semente NÃºmero inteiro para semente aleatÃ³ria.
#' @return Retorna o data.frame contendo todos os dados dos conglomerados sorteados.
#' @export
#' @examples
#' # Sorteia 2 espÃ©cies inteiras no banco iris:
#' amostra_conglomerados(iris, conglomerado = Species, n_conglomerados = 2)
amostra_conglomerados <- function(dados, conglomerado, n_conglomerados = NULL, proporcao = NULL, semente = NULL) {
  if (!is.data.frame(dados)) stop("O argumento 'dados' deve ser um data.frame.")
  if (!is.null(semente)) set.seed(semente)
  
  N_total <- nrow(dados)
  
  # Captura conglomerado via NSE ou string
  cong_sub <- substitute(conglomerado)
  cong_nome <- deparse(cong_sub)
  cong_nome <- sub(".*\\$", "", cong_nome)
  
  if (cong_nome %in% names(dados)) {
    cong_vec <- dados[[cong_nome]]
    nome_exibicao <- cong_nome
  } else {
    cong_val <- tryCatch(eval(cong_sub, parent.frame()), error = function(e) NULL)
    if (is.character(cong_val) && length(cong_val) == 1 && cong_val %in% names(dados)) {
      cong_vec <- dados[[cong_val]]
      nome_exibicao <- cong_val
    } else if (length(cong_val) == N_total) {
      cong_vec <- cong_val
      nome_exibicao <- cong_nome
    } else {
      stop(sprintf("A coluna de conglomerado '%s' nÃ£o foi encontrada no banco de dados.", cong_nome))
    }
  }
  
  todos_congs <- unique(cong_vec)
  M_congs <- length(todos_congs)
  
  if (is.null(n_conglomerados) && is.null(proporcao)) {
    stop("Informe a quantidade de conglomerados ('n_conglomerados') ou a 'proporcao' a sortear.")
  }
  
  if (!is.null(proporcao)) {
    if (proporcao <= 0 || proporcao > 1) stop("A 'proporcao' deve estar entre 0 e 1.")
    m_sorteio <- round(M_congs * proporcao)
  } else {
    m_sorteio <- n_conglomerados
  }
  
  if (m_sorteio <= 0 || m_sorteio > M_congs) {
    stop("A quantidade de conglomerados a sortear deve estar entre 1 e o total existente.")
  }
  
  congs_sorteados <- sample(todos_congs, size = m_sorteio, replace = FALSE)
  indices <- which(cong_vec %in% congs_sorteados)
  amostra_final <- dados[indices, , drop = FALSE]
  
  n_individuos <- nrow(amostra_final)
  pct_individuos <- (n_individuos / N_total) * 100
  pct_congs <- (m_sorteio / M_congs) * 100
  
  str_congs_sorteados <- paste(as.character(congs_sorteados), collapse = ", ")
  if (nchar(str_congs_sorteados) > 50) {
    str_congs_sorteados <- paste0(substr(str_congs_sorteados, 1, 47), "...")
  }
  
  tit <- sprintf("AMOSTRAGEM POR CONGLOMERADOS (Por: %s)", nome_exibicao)
  .print_titulo(tit)
  cat(sprintf("Total de Conglomerados na PopulaÃ§Ã£o: %d\n", M_congs))
  cat(sprintf("Conglomerados Sorteados:             %d (%s%%)\n", m_sorteio, .fmt_num(pct_congs, 1)))
  cat(sprintf("Grupos Selecionados:                 %s\n", str_congs_sorteados))
  cat(sprintf("Total de IndivÃ­duos Amostrados:      %s observaÃ§Ãµes (%s%% da populaÃ§Ã£o)\n",
              format(n_individuos, big.mark = ".", decimal.mark = ","),
              .fmt_num(pct_individuos, 1)))
  .print_rodape()
  
  invisible(amostra_final)
}
