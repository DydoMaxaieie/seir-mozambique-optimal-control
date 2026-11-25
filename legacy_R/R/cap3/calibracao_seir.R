# ================================================
# calibracao_seir.R — versão final completa
# Calibração SEIR — Moçambique Jan-Mar 2021
# ================================================

# ------------------------------------------------
# 1. VARIÁVEIS GLOBAIS
# ------------------------------------------------
N_moz       <- 31e6
R0_moz      <- 18642    # recuperados até Dez 2020
sigma_fixo  <- 1/5.2    # fixo da literatura

# ------------------------------------------------
# 2. DOWNLOAD DADOS OWID
# ------------------------------------------------
url <- paste0(
  "https://raw.githubusercontent.com/",
  "owid/covid-19-data/master/public/data/",
  "owid-covid-data.csv"
)

cat("A descarregar dados OWID...\n")
tryCatch({
  owid <- read.csv(url)
  cat("Sucesso! Linhas:", nrow(owid), "\n")
}, error = function(e) {
  stop("Erro no download: ", conditionMessage(e))
})

# ------------------------------------------------
# 3. FILTRAR MOÇAMBIQUE — Jan-Mar 2021
# ------------------------------------------------
owid$date <- as.Date(owid$date)
moz       <- owid[owid$iso_code == "MOZ", ]

q1_2021 <- moz[
  moz$date >= as.Date("2021-01-01") &
  moz$date <= as.Date("2021-03-31"),
  c("date", "new_cases", "new_cases_smoothed",
    "total_cases", "new_deaths", "total_deaths")
]

cat("Número de dias:", nrow(q1_2021), "\n")
cat("Valores em falta:\n")
print(colSums(is.na(q1_2021)))

cat("\nTotais mensais — verificação:\n")
cat("Janeiro:",
    sum(q1_2021$new_cases[
      q1_2021$date <= as.Date("2021-01-31")],
      na.rm=TRUE), "\n")
cat("Fevereiro:",
    sum(q1_2021$new_cases[
      q1_2021$date >= as.Date("2021-02-01") &
      q1_2021$date <= as.Date("2021-02-28")],
      na.rm=TRUE), "\n")
cat("Março:",
    sum(q1_2021$new_cases[
      q1_2021$date >= as.Date("2021-03-01")],
      na.rm=TRUE), "\n")

# ------------------------------------------------
# 4. PREPARAR DADOS PARA CALIBRAÇÃO
# ------------------------------------------------
novos_obs <- q1_2021$new_cases_smoothed
datas     <- q1_2021$date
dias      <- as.numeric(datas - datas[1])

# Condições iniciais
I0_init <- round(novos_obs[1] * 10)
E0_init <- round(I0_init * 1.5)
S0_init <- N_moz - E0_init - I0_init - R0_moz

cat("\nCondições iniciais base:\n")
cat("S(0) =", S0_init, "\n")
cat("E(0) =", E0_init, "\n")
cat("I(0) =", I0_init, "\n")
cat("R(0) =", R0_moz,  "\n")
cat("novos_obs max =", round(max(novos_obs)),
    "no dia", dias[which.max(novos_obs)], "\n")

# ------------------------------------------------
# 5. INTEGRADOR SEIR RK4
# ------------------------------------------------
seir_rk4 <- function(beta, gamma, S0, E0,
                     I0, R0, N, T_dias) {
  sigma <- sigma_fixo
  n     <- T_dias + 1

  S  <- numeric(n); S[1]  <- S0
  E  <- numeric(n); E[1]  <- E0
  Iv <- numeric(n); Iv[1] <- I0
  R  <- numeric(n); R[1]  <- R0
  NC <- numeric(n); NC[1] <- beta * S0 * I0 / N

  for (i in 1:T_dias) {
    s <- S[i]; e <- E[i]
    ii <- Iv[i]; r <- R[i]

    force <- function(ss, ii) beta * ss * ii / N

    k1s <- -force(s, ii)
    k1e <-  force(s, ii) - sigma * e
    k1i <-  sigma * e    - gamma * ii
    k1r <-  gamma * ii

    s2 <- s+0.5*k1s; e2 <- e+0.5*k1e
    i2 <- ii+0.5*k1i

    k2s <- -force(s2, i2)
    k2e <-  force(s2, i2) - sigma * e2
    k2i <-  sigma * e2    - gamma * i2
    k2r <-  gamma * i2

    s3 <- s+0.5*k2s; e3 <- e+0.5*k2e
    i3 <- ii+0.5*k2i

    k3s <- -force(s3, i3)
    k3e <-  force(s3, i3) - sigma * e3
    k3i <-  sigma * e3    - gamma * i3
    k3r <-  gamma * i3

    s4 <- s+k3s; e4 <- e+k3e
    i4 <- ii+k3i

    k4s <- -force(s4, i4)
    k4e <-  force(s4, i4) - sigma * e4
    k4i <-  sigma * e4    - gamma * i4
    k4r <-  gamma * i4

    S[i+1]  <- s  + (k1s+2*k2s+2*k3s+k4s)/6
    E[i+1]  <- e  + (k1e+2*k2e+2*k3e+k4e)/6
    Iv[i+1] <- ii + (k1i+2*k2i+2*k3i+k4i)/6
    R[i+1]  <- r  + (k1r+2*k2r+2*k3r+k4r)/6
    NC[i+1] <- force(S[i+1], Iv[i+1])
  }

  data.frame(S=S, E=E, I=Iv, R=R, NC=NC)
}

# ------------------------------------------------
# 6. FUNÇÃO DE CUSTO
# ------------------------------------------------
custo_mq <- function(params) {
  beta  <- params[1]
  gamma <- params[2]
  I0    <- params[3]

  # Restrições biológicas
  if (beta  < 0.10 || beta  > 1.50) return(1e12)
  if (gamma < 0.05 || gamma > 0.33) return(1e12)
  if (I0    < 100  || I0    > 5e4)  return(1e12)

  I0_cal <- round(I0)
  E0_cal <- round(I0_cal * 1.5)
  S0_cal <- N_moz - E0_cal - I0_cal - R0_moz
  if (S0_cal <= 0) return(1e12)

  sol <- tryCatch(
    seir_rk4(beta, gamma,
             S0_cal, E0_cal, I0_cal, R0_moz,
             N_moz, 89),
    error = function(e) NULL
  )
  if (is.null(sol)) return(1e12)

  # RMSE relativo sobre novos casos diários
  sqrt(mean((sol$NC[-1] - novos_obs[-1])^2)) /
    mean(novos_obs)
}

# ------------------------------------------------
# 7. OPTIMIZAÇÃO — GRELHA DE PONTOS INICIAIS
# ------------------------------------------------
cat("\nA calibrar parâmetros...\n")

beta_grid  <- c(0.20, 0.30, 0.40, 0.50, 0.60)
gamma_grid <- c(0.07, 0.10, 0.14)
I0_grid    <- c(500, 1000, 2000, 5000)

melhor_custo  <- Inf
melhor_params <- NULL
total_iter    <- 0

for (b in beta_grid) {
  for (g in gamma_grid) {
    for (i0 in I0_grid) {
      res <- optim(
        par     = c(b, g, i0),
        fn      = custo_mq,
        method  = "Nelder-Mead",
        control = list(maxit  = 10000,
                       reltol = 1e-10)
      )
      total_iter <- total_iter + 1
      if (res$value < melhor_custo) {
        melhor_custo  <- res$value
        melhor_params <- res$par
        cat("Novo mínimo [", total_iter, "]:",
            "beta=", round(res$par[1], 4),
            "gamma=", round(res$par[2], 4),
            "I0=", round(res$par[3]),
            "RMSE%=", round(res$value*100, 2),
            "\n")
      }
    }
  }
}

# ------------------------------------------------
# 8. RESULTADOS FINAIS
# ------------------------------------------------
beta_cal  <- melhor_params[1]
gamma_cal <- melhor_params[2]
I0_cal    <- round(melhor_params[3])
E0_cal    <- round(I0_cal * 1.5)
S0_cal    <- N_moz - E0_cal - I0_cal - R0_moz
R0_num    <- beta_cal / gamma_cal

cat("\n=== PARÂMETROS CALIBRADOS ===\n")
cat("beta  =", round(beta_cal,  4), "\n")
cat("sigma =", round(sigma_fixo, 4),
    "(fixo — WHO/literatura)\n")
cat("gamma =", round(gamma_cal, 4), "\n")
cat("I(0)  =", I0_cal, "\n")
cat("R0    =", round(R0_num, 3), "\n")
cat("RMSE% =", round(melhor_custo*100, 2), "%\n")

cat("\n=== CONDIÇÕES INICIAIS CALIBRADAS ===\n")
cat("S(0) =", S0_cal, "\n")
cat("E(0) =", E0_cal, "\n")
cat("I(0) =", I0_cal, "\n")
cat("R(0) =", R0_moz,  "\n")

cat("\n=== VALIDAÇÃO BIOLÓGICA ===\n")
cat("R0 esperado 501.V2: 1.5 — 2.5\n")
cat("R0 obtido:         ", round(R0_num, 3), "\n")
if (R0_num >= 1.2 && R0_num <= 3.5) {
  cat("STATUS: PLAUSÍVEL\n")
} else {
  cat("STATUS: REVER\n")
}

# Guardar parâmetros para uso nos gráficos
saveRDS(list(
  beta      = beta_cal,
  gamma     = gamma_cal,
  sigma     = sigma_fixo,
  I0        = I0_cal,
  E0        = E0_cal,
  S0        = S0_cal,
  R0_moz    = R0_moz,
  N_moz     = N_moz,
  novos_obs = novos_obs,
  datas     = datas,
  R0_num    = R0_num
), "params_calibrados.rds")

cat("\nParâmetros guardados em params_calibrados.rds\n")