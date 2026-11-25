# ================================================
# seir_meta.R — Modelo Metapopulacional Completo
# SEIR-V Corredor N1, Moçambique Jan-Mar 2021
# Versão com FBS diferenciado por faixa etária
# ================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

dir.create("Figuras",    showWarnings = FALSE)
dir.create("resultados", showWarnings = FALSE)

# ------------------------------------------------
# 1. PARÂMETROS GLOBAIS
# ------------------------------------------------
cat("A definir parâmetros...\n")

K <- 6
G <- 9

nos_nomes <- c("Maputo","Xai-Xai","Maxixe",
               "Caia","Inchope","Namacurra")

grupos <- c("0-4","5-14","15-24","25-34",
            "35-44","45-54","55-64",
            "65-74","75+")

N_nos <- c(1600,980,750,610,700,840) * 1000

prop_etaria <- c(0.148,0.229,0.197,0.152,
                 0.103,0.063,0.038,0.021,0.011)
prop_etaria <- prop_etaria / sum(prop_etaria)

N_kg    <- outer(N_nos, prop_etaria)
N_total <- sum(N_kg)
cat("Pop. total corredor N1:",
    format(round(N_total), big.mark=","), "\n")

sigma    <- 1/5.2
gamma    <- 0.08
T_V      <- 14
T_Q      <- 7
T_Hw     <- 8
T_Hc     <- 12
e_vac    <- 0.70
delta    <- 0.50
p_x      <- 0.23
mu_decay <- 0.017

p_h_vec  <- c(0.010,0.005,0.010,0.020,
              0.040,0.080,0.130,0.200,0.300)
p_c_vec  <- c(0.050,0.030,0.040,0.060,
              0.080,0.120,0.180,0.250,0.350)
mu_g_vec <- c(0.0002,0.0001,0.0002,0.0005,
              0.0010,0.0030,0.0080,0.0200,
              0.0500)
kappa <- c(1.0,1.2,1.4,1.8,1.6,1.5)

p_h_kg <- matrix(0, K, G)
p_c_kg <- matrix(0, K, G)
mu_kg  <- matrix(0, K, G)
for (k in 1:K) {
  p_h_kg[k,] <- p_h_vec
  p_c_kg[k,] <- p_c_vec
  mu_kg[k,]  <- mu_g_vec * kappa[k]
}

# Pesos etários A_{kg}
# Reflectem vulnerabilidade clínica
# documentada nos boletins MISAU
omega_g <- c(
  0.5,  # 0-4
  0.3,  # 5-14
  0.4,  # 15-24
  0.6,  # 25-34
  0.8,  # 35-44
  1.0,  # 45-54 — referência
  1.5,  # 55-64
  3.0,  # 65-74
  5.0   # 75+
)

# Matriz A_kg base: mortalidade * peso etário
A_kg_base <- matrix(0, K, G)
for (k in 1:K) {
  A_kg_base[k,] <- mu_g_vec *
                   kappa[k] * omega_g
}
# Normalizar para que média = 1
A_kg_base <- A_kg_base /
             mean(A_kg_base)

cat("Pesos A_kg (média por faixa):\n")
for (g in 1:G) {
  cat(sprintf("  %-8s %.3f\n",
      grupos[g], mean(A_kg_base[,g])))
}

# Matriz de mobilidade Θ
Theta <- matrix(1e-4, K, K)
diag(Theta) <- c(0.880,0.872,0.881,
                 0.862,0.884,0.926)
Theta[1,2] <- 0.072; Theta[2,1] <- 0.065
Theta[2,3] <- 0.048; Theta[3,2] <- 0.042
Theta[3,4] <- 0.058; Theta[4,3] <- 0.051
Theta[4,5] <- 0.042; Theta[5,4] <- 0.038
Theta[5,6] <- 0.061; Theta[6,5] <- 0.054
Theta[1,4] <- 0.030; Theta[4,1] <- 0.025
for (i in 1:K) {
  Theta[i,] <- Theta[i,] / sum(Theta[i,])
}
cat("Somas Theta:\n")
print(round(rowSums(Theta), 6))

# Matriz de contacto C (Prem 2017)
C_mat <- matrix(c(
  1.92,0.83,0.41,1.21,0.97,0.28,0.12,0.08,0.04,
  0.61,4.87,1.23,0.58,0.72,0.31,0.15,0.07,0.03,
  0.38,1.45,5.12,1.87,0.64,0.42,0.18,0.09,0.04,
  1.18,0.72,1.93,3.84,1.42,0.61,0.24,0.11,0.05,
  1.02,0.98,0.71,1.58,3.21,0.87,0.31,0.14,0.06,
  0.29,0.41,0.58,0.74,1.12,2.87,0.72,0.28,0.09,
  0.14,0.22,0.29,0.38,0.51,0.84,2.14,0.48,0.15,
  0.09,0.12,0.16,0.19,0.27,0.41,0.57,1.62,0.31,
  0.05,0.07,0.09,0.12,0.15,0.18,0.24,0.38,1.12
), nrow=G, byrow=TRUE)

params <- list(
  K=K, G=G,
  sigma=sigma, gamma=gamma,
  T_V=T_V, T_Q=T_Q,
  T_Hw=T_Hw, T_Hc=T_Hc,
  e=e_vac, delta=delta,
  p_h=p_h_kg, p_c=p_c_kg,
  mu_kg=mu_kg,
  Theta=Theta, C_mat=C_mat,
  N_kg=N_kg
)

# ------------------------------------------------
# 2. CARREGAR DADOS MISAU
# ------------------------------------------------
cat("\nA carregar dados MISAU...\n")

dados <- read.csv("dados/covid_moz_Q1_2021.csv")
dados$Date <- as.Date(dados$Date)

dias_boletins <- c(0, 33, 64)
casos_misau   <- matrix(0, K, 3)

for (j in 1:3) {
  d_j <- as.Date("2021-01-01") +
         dias_boletins[j]
  sub <- dados[dados$Date == d_j, ]
  for (k in 1:K) {
    no_k <- paste0("k", k)
    idx  <- which(sub$Node == no_k)
    if (length(idx) > 0) {
      casos_misau[k,j] <-
        sub$Total_Cases[idx[1]]
    }
  }
}

cat("Casos por nó (Jan→Mar):\n")
for (k in 1:K) {
  cat(sprintf("  %-12s: %5d → %5d\n",
      nos_nomes[k],
      casos_misau[k,1],
      casos_misau[k,3]))
}

# ------------------------------------------------
# 3. CONDIÇÕES INICIAIS
# ------------------------------------------------
cat("\nA calcular condições iniciais...\n")

subtotal_n1   <- sum(casos_misau[,1])
rec_total_jan <- 16680
obi_total_jan <- 167
total_jan_nac <- 18794

rec_n1 <- round(rec_total_jan *
                subtotal_n1/total_jan_nac)
obi_n1 <- round(obi_total_jan *
                subtotal_n1/total_jan_nac)

I0_nos <- numeric(K)
R0_nos <- numeric(K)

for (k in 1:K) {
  prop_k    <- casos_misau[k,1]/subtotal_n1
  R0_nos[k] <- round(rec_n1 * prop_k)
  obi_k     <- round(obi_n1 * prop_k)
  activos_k <- casos_misau[k,1] -
               R0_nos[k] - obi_k
  I0_nos[k] <- max(10, activos_k)
}

cat("Condições iniciais por nó:\n")
for (k in 1:K) {
  cat(sprintf("  %-12s I0=%4d R0=%5d\n",
      nos_nomes[k], I0_nos[k], R0_nos[k]))
}

# ------------------------------------------------
# 4. CALIBRAÇÃO β_k
# ------------------------------------------------
cat("\nA calibrar beta_k...\n")

seir_simples <- function(beta0, S0, I0,
                         R0, N, T_dias) {
  n    <- T_dias + 1
  S    <- numeric(n); S[1]  <- S0
  E    <- numeric(n); E[1]  <- round(I0*1.5)
  Iv   <- numeric(n); Iv[1] <- I0
  R    <- numeric(n); R[1]  <- R0
  C    <- numeric(n); C[1]  <- I0 + R0
  for (i in 1:T_dias) {
    bt     <- beta0 * exp(-mu_decay*(i-1))
    f      <- bt * S[i] * Iv[i] / N
    S[i+1] <- max(0, S[i]  - f)
    E[i+1] <- max(0, E[i]  + f -
                  sigma*E[i])
    Iv[i+1]<- max(0, Iv[i] +
                  sigma*E[i] - gamma*Iv[i])
    R[i+1] <- max(0, R[i]  + gamma*Iv[i])
    C[i+1] <- Iv[i+1] + R[i+1]
  }
  C
}

beta_k   <- numeric(K)
rmse_pct <- numeric(K)

for (k in 1:K) {
  C_alvo <- casos_misau[k,]
  S0_k   <- N_nos[k] - I0_nos[k] - R0_nos[k]

  custo_k <- function(b) {
    if (b <= 0 || b > 1.5) return(1e12)
    C_sim <- seir_simples(
      b, S0_k, I0_nos[k],
      R0_nos[k], N_nos[k], 64)
    sum((C_sim[dias_boletins+1] -
         C_alvo)^2) / mean(C_alvo)^2
  }

  res        <- optimize(custo_k,
                         c(0.05,1.0),
                         tol=1e-10)
  beta_k[k]  <- res$minimum
  C_sim      <- seir_simples(
    beta_k[k], S0_k, I0_nos[k],
    R0_nos[k], N_nos[k], 64)
  rmse       <- sqrt(mean(
    (C_sim[dias_boletins+1]-C_alvo)^2))
  rmse_pct[k] <- rmse/mean(C_alvo)*100

  cat(sprintf(
    "  %-12s beta=%.4f R0=%.3f RMSE%%=%.1f\n",
    nos_nomes[k], beta_k[k],
    beta_k[k]/gamma, rmse_pct[k]))
}
cat("RMSE% global:",
    round(mean(rmse_pct),2),"%\n")

# ------------------------------------------------
# 5. FORÇA DE INFECÇÃO
# ------------------------------------------------
calc_lambda <- function(I_kg, Iv_kg,
                        u_kg, beta_k, p) {
  ensure_mat <- function(x) {
    matrix(as.numeric(x), p$K, p$G)
  }
  I_kg  <- ensure_mat(I_kg)
  Iv_kg <- ensure_mat(Iv_kg)
  u_kg  <- ensure_mat(u_kg)

  N_hat <- numeric(p$K)
  for (m in 1:p$K) {
    for (l in 1:p$K) {
      N_hat[m] <- N_hat[m] +
        sum(p$Theta[l,m] * p$N_kg[l,])
    }
  }
  N_hat <- pmax(N_hat, 1)

  lam <- matrix(0, p$K, p$G)
  for (k in 1:p$K) {
    for (g in 1:p$G) {
      soma <- 0
      for (m in 1:p$K) {
        inf_m <- numeric(p$G)
        for (l in 1:p$K) {
          inf_m <- inf_m +
            p$Theta[l,m] *
            (I_kg[l,] +
             p$delta * Iv_kg[l,])
        }
        soma <- soma +
          p$Theta[k,m] *
          sum(p$C_mat[g,] *
              inf_m / N_hat[m])
      }
      lam[k,g] <- beta_k[k] *
                  (1-u_kg[k,g]) * soma
    }
  }
  lam
}

# ------------------------------------------------
# 6. INTEGRADOR RK4
# ------------------------------------------------
rk4_meta <- function(x, v_kg, u_kg,
                     beta_k, p) {
  K <- p$K; G <- p$G

  f_dx <- function(x, v_kg, u_kg) {
    lam <- calc_lambda(
      x$I, x$Iv, u_kg, beta_k, p)
    list(
      Sx = -lam * x$Sx,
      Su = -(lam + v_kg) * x$Su,
      Sv =  v_kg * x$Su -
            (lam + 1/p$T_V) * x$Sv,
      Sp = (1-p$e)*x$Sv/p$T_V -
           (1-p$e)*(1-u_kg)*lam*x$Sp,
      E  =  lam*(x$Sx+x$Su+x$Sv) -
            p$sigma*x$E,
      Ev = (1-p$e)*(1-u_kg)*lam*x$Sp -
            p$sigma*x$Ev,
      I  =  p$sigma*x$E  - p$gamma*x$I,
      Iv =  p$sigma*x$Ev - p$gamma*x$Iv,
      Q  = (1-p$p_h)*p$gamma*(x$I+x$Iv) -
            x$Q/p$T_Q,
      Hw =  p$p_h*p$gamma*(x$I+x$Iv) -
            (1/p$T_Hw+p$mu_kg)*x$Hw,
      Hc =  p$p_c*x$Hw/p$T_Hw -
            (1/p$T_Hc+p$mu_kg)*x$Hc,
      R  =  p$e*x$Sv/p$T_V +
            x$Q/p$T_Q,
      RH = (1-p$p_c)*x$Hw/p$T_Hw +
            x$Hc/p$T_Hc,
      D  =  p$mu_kg*x$Q,
      DH =  p$mu_kg*(x$Hw+x$Hc),
      Vw =  v_kg*x$Su
    )
  }

  add_x <- function(a, b, s=1) {
    nms <- names(a)
    setNames(
      lapply(nms, function(n)
        matrix(as.numeric(a[[n]]) +
               s*as.numeric(b[[n]]),
               K, G)), nms)
  }

  k1 <- f_dx(x, v_kg, u_kg)
  k2 <- f_dx(add_x(x,k1,.5), v_kg, u_kg)
  k3 <- f_dx(add_x(x,k2,.5), v_kg, u_kg)
  k4 <- f_dx(add_x(x,k3,1.), v_kg, u_kg)

  nms <- names(x)
  setNames(
    lapply(nms, function(n) {
      val <- as.numeric(x[[n]]) +
             (as.numeric(k1[[n]]) +
              2*as.numeric(k2[[n]]) +
              2*as.numeric(k3[[n]]) +
              as.numeric(k4[[n]])) / 6
      matrix(pmax(0, val), K, G)
    }), nms)
}

# ------------------------------------------------
# 7. CONDIÇÕES INICIAIS DO SISTEMA
# ------------------------------------------------
cat("\nA inicializar estado...\n")

nomes_comp <- c("Sx","Su","Sv","Sp",
                "E","Ev","I","Iv",
                "Q","Hw","Hc","R",
                "RH","D","DH","Vw")

init_estado <- function() {
  x0 <- setNames(
    lapply(nomes_comp, function(n)
      matrix(0, K, G)),
    nomes_comp)
  for (k in 1:K) {
    I0g <- round(I0_nos[k] * prop_etaria)
    R0g <- round(R0_nos[k] * prop_etaria)
    E0g <- round(I0g * 1.5)
    S0g <- pmax(0,
                N_kg[k,]-I0g-R0g-E0g)
    x0$Sx[k,] <- round(p_x     * S0g)
    x0$Su[k,] <- round((1-p_x) * S0g)
    x0$E[k,]  <- E0g
    x0$I[k,]  <- I0g
    x0$R[k,]  <- R0g
  }
  x0
}

x0     <- init_estado()
T_dias <- 89
tvec   <- 0:T_dias

cat("I(0):", sum(x0$I),
    "| S(0):", sum(x0$Su)+sum(x0$Sx), "\n")

# ------------------------------------------------
# 8. FUNÇÃO DE SIMULAÇÃO
# ------------------------------------------------
simular <- function(x0, v_kg, u_kg,
                    beta_k, p,
                    T_dias, nome) {
  cat("  Cenário:", nome, "\n")
  N_t  <- T_dias + 1
  hist <- vector("list", N_t)
  hist[[1]] <- x0

  for (t in 1:T_dias) {
    hist[[t+1]] <- rk4_meta(
      hist[[t]], v_kg, u_kg, beta_k, p)
  }

  # Séries temporais globais
  I_t  <- numeric(N_t)
  D_t  <- numeric(N_t)
  Hw_t <- numeric(N_t)
  R_t  <- numeric(N_t)

  for (ti in 1:N_t) {
    h        <- hist[[ti]]
    I_t[ti]  <- sum(h$I)  + sum(h$Iv)
    D_t[ti]  <- sum(h$D)  + sum(h$DH)
    Hw_t[ti] <- sum(h$Hw)
    R_t[ti]  <- sum(h$R)  + sum(h$RH)
  }

  # Por nó
  I_nos  <- matrix(0, K, N_t)
  D_nos  <- matrix(0, K, N_t)
  Rt_nos <- matrix(0, K, N_t)

  for (ti in 1:N_t) {
    h <- hist[[ti]]
    for (k in 1:K) {
      I_nos[k,ti] <- sum(h$I[k,]) +
                     sum(h$Iv[k,])
      D_nos[k,ti] <- sum(h$D[k,]) +
                     sum(h$DH[k,])
      S_k <- sum(h$Sx[k,]) +
             sum(h$Su[k,]) +
             sum(h$Sv[k,])
      Rt_nos[k,ti] <- beta_k[k]*S_k /
        (p$gamma*sum(p$N_kg[k,]))
    }
  }

  # Por faixa etária
  I_grupos <- matrix(0, G, N_t)
  D_grupos <- matrix(0, G, N_t)

  for (ti in 1:N_t) {
    h <- hist[[ti]]
    for (g in 1:G) {
      I_grupos[g,ti] <- sum(h$I[,g]) +
                        sum(h$Iv[,g])
      D_grupos[g,ti] <- sum(h$D[,g]) +
                        sum(h$DH[,g])
    }
  }

  list(
    nome=nome, I_t=I_t, D_t=D_t,
    Hw_t=Hw_t, R_t=R_t,
    I_nos=I_nos, D_nos=D_nos,
    Rt_nos=Rt_nos,
    I_grupos=I_grupos,
    D_grupos=D_grupos,
    hist=hist,
    v_kg=v_kg, u_kg=u_kg)
}

# ------------------------------------------------
# 9. FBS METAPOPULACIONAL
#    COM ADJUNTA DIFERENCIADA POR FAIXA
# ------------------------------------------------
fbs_meta <- function(x0, beta_k, p,
                     T_dias, A_kg,
                     B_peso, C_peso,
                     v_max, u_max,
                     nome) {
  cat("  FBS:", nome, "\n")
  K <- p$K; G <- p$G
  v_kg  <- matrix(0, K, G)
  u_kg  <- matrix(0, K, G)
  omega <- 0.5
  tol   <- 1e-6
  err   <- Inf
  it    <- 0
  N_t   <- T_dias + 1

  while (err > tol && it < 500) {
    v_old <- v_kg
    u_old <- u_kg

    # Forward
    hist <- vector("list", N_t)
    hist[[1]] <- x0
    for (t in 1:T_dias) {
      hist[[t+1]] <- rk4_meta(
        hist[[t]], v_kg, u_kg,
        beta_k, p)
    }

    # Backward — adjunta K×G
    # q_I_{kg}(t) = A_{kg}/gamma *
    #   (1 - exp(-gamma*(T-t)))
    v_num <- matrix(0, K, G)
    u_num <- matrix(0, K, G)
    count <- 0

    for (ti in 1:N_t) {
      tau <- T_dias - (ti - 1)

      # Adjunta diferenciada por (k,g)
      q_I_kg <- A_kg / p$gamma *
                (1 - exp(-p$gamma * tau))

      x_t <- hist[[ti]]
      lam <- calc_lambda(
        x_t$I, x_t$Iv,
        u_kg, beta_k, p)

      # Condição optimalidade vacinação
      # v* = (q_Su - q_Sv) * Su / B
      q_Su_mat <- lam * q_I_kg
      q_Sv_mat <- q_Su_mat * p$e
      vh <- (q_Su_mat - q_Sv_mat) *
             x_t$Su / B_peso
      vh <- matrix(pmax(pmin(
        as.numeric(vh),
        as.numeric(v_max)), 0), K, G)

      # Condição optimalidade distanciamento
      # u* = Phi_{kg} / C
      S_tot <- x_t$Sx + x_t$Su + x_t$Sv
      uh <- lam * q_I_kg * S_tot / C_peso
      uh <- matrix(pmax(pmin(
        as.numeric(uh), u_max), 0), K, G)

      v_num <- v_num + vh
      u_num <- u_num + uh
      count <- count + 1
    }

    v_hat <- v_num / count
    u_hat <- u_num / count
    v_kg  <- omega*v_hat + (1-omega)*v_old
    u_kg  <- omega*u_hat + (1-omega)*u_old

    err <- max(
      max(abs(v_kg-v_old)) /
        (max(abs(v_kg))+1e-10),
      max(abs(u_kg-u_old)) /
        (max(abs(u_kg))+1e-10))
    it <- it + 1
  }

  cat(sprintf(
    "    Iter:%d Erro:%.2e u_med:%.4f\n",
    it, err, mean(u_kg)))

  res      <- simular(x0, v_kg, u_kg,
                       beta_k, p,
                       T_dias, nome)
  res$v_kg <- v_kg
  res$u_kg <- u_kg
  res
}

# ------------------------------------------------
# 10. SIMULAR TODOS OS CENÁRIOS
# ------------------------------------------------
cat("\nA simular cenários...\n")

v_zero <- matrix(0,    K, G)
u_zero <- matrix(0,    K, G)
v_max  <- matrix(0.05, K, G)
u_max  <- 0.60
B_peso <- 100
C_peso <- 200

# Cenário 0 — sem intervenção
res0 <- simular(x0, v_zero, u_zero,
                 beta_k, params,
                 T_dias, "0-Base")

cat("Infectados máx:",
    round(max(res0$I_t)),
    "no dia", which.max(res0$I_t), "\n")
cat("Óbitos totais:",
    round(res0$D_t[T_dias+1]), "\n")

# Cenário A — prioridade etária 65+
# Amplificar pesos grupos 8 e 9
A_kg_A <- A_kg_base
A_kg_A[,8:9] <- A_kg_A[,8:9] * 5

resA <- fbs_meta(x0, beta_k, params,
                  T_dias, A_kg_A,
                  B_peso, C_peso,
                  v_max, u_max,
                  "A-Prior.Etaria")

# Cenário B — prioridade regional Maputo
# Reduzir custo de vacinação em k=1
A_kg_B <- A_kg_base
B_kg_B <- matrix(B_peso, K, G)
B_kg_B[1,] <- B_peso / 10

resB <- fbs_meta(x0, beta_k, params,
                  T_dias, A_kg_base,
                  mean(B_kg_B), C_peso,
                  v_max, u_max,
                  "B-Prior.Regional")

# Cenário C — controlo óptimo simultâneo
resC <- fbs_meta(x0, beta_k, params,
                  T_dias, A_kg_base,
                  B_peso, C_peso,
                  v_max, u_max,
                  "C-Ctrl.Optimo")

# Cenário D — só distanciamento
resD <- fbs_meta(x0, beta_k, params,
                  T_dias, A_kg_base,
                  B_peso, C_peso,
                  v_zero, u_max,
                  "D-So.Distanciamento")

# ------------------------------------------------
# 11. TABELA DE RESULTADOS
# ------------------------------------------------
D0 <- res0$D_t[T_dias+1]

cat("\n=== TABELA DE RESULTADOS ===\n")
cat(sprintf("%-28s %8s %10s %8s %8s\n",
    "Cenário","Óbitos",
    "Evitados","Red%","Rt<1"))
cat(strrep("-",66),"\n")
cat(sprintf("%-28s %8.0f %10s %8s %8s\n",
    "0-Base", D0, "---","---","---"))

for (res in list(resA,resB,resC,resD)) {
  D          <- res$D_t[T_dias+1]
  ev         <- D0 - D
  pct        <- ev/D0*100
  dia_limiar <- which(
    colMeans(res$Rt_nos) < 1)[1]
  dia_str    <- ifelse(
    is.na(dia_limiar),">89",
    as.character(dia_limiar))
  cat(sprintf(
    "%-28s %8.0f %10.0f %7.1f%% %8s\n",
    res$nome, D, ev, pct, dia_str))
}

# ------------------------------------------------
# 12. TABELA ÓBITOS POR FAIXA ETÁRIA
# ------------------------------------------------
cat("\n=== ÓBITOS POR FAIXA ETÁRIA ===\n")
cat(sprintf("%-8s %8s %8s %8s %8s %8s %8s\n",
    "Faixa","Base","Cen.A",
    "Cen.B","Cen.C","Cen.D","Red%C"))
cat(strrep("-",60),"\n")

tab_etaria <- data.frame()
for (g in 1:G) {
  D0_g  <- res0$D_grupos[g,T_dias+1]
  DA_g  <- resA$D_grupos[g,T_dias+1]
  DB_g  <- resB$D_grupos[g,T_dias+1]
  DC_g  <- resC$D_grupos[g,T_dias+1]
  DD_g  <- resD$D_grupos[g,T_dias+1]
  red_C <- ifelse(D0_g > 0,
                  (D0_g-DC_g)/D0_g*100, 0)

  cat(sprintf(
    "%-8s %8.1f %8.1f %8.1f %8.1f %8.1f %7.1f%%\n",
    grupos[g], D0_g, DA_g,
    DB_g, DC_g, DD_g, red_C))

  tab_etaria <- rbind(tab_etaria,
    data.frame(
      Faixa   = grupos[g],
      Base    = D0_g,
      CenA    = DA_g,
      CenB    = DB_g,
      CenC    = DC_g,
      CenD    = DD_g,
      Red_pct = red_C))
}

# ------------------------------------------------
# 13. TABELA ACE
# ------------------------------------------------
cat("\n=== TABELA ACE ===\n")
tab_ace <- data.frame(
  Cenario = c("A-Prior.Etaria",
              "B-Prior.Regional",
              "C-Ctrl.Optimo",
              "D-So.Distanciamento"),
  Obitos  = c(resA$D_t[T_dias+1],
              resB$D_t[T_dias+1],
              resC$D_t[T_dias+1],
              resD$D_t[T_dias+1])
)
tab_ace$Evitados <- D0 - tab_ace$Obitos
tab_ace$Red_pct  <- tab_ace$Evitados/D0*100
tab_ace$ACER     <- tab_ace$Obitos/
                    tab_ace$Evitados

cat(sprintf("%-28s %8s %10s %8s %8s\n",
    "Cenário","Óbitos",
    "Evitados","Red%","ACER"))
cat(strrep("-",66),"\n")
for (i in 1:nrow(tab_ace)) {
  cat(sprintf(
    "%-28s %8.1f %10.1f %7.1f%% %8.4f\n",
    tab_ace$Cenario[i],
    tab_ace$Obitos[i],
    tab_ace$Evitados[i],
    tab_ace$Red_pct[i],
    tab_ace$ACER[i]))
}

# ------------------------------------------------
# 14. GUARDAR RESULTADOS
# ------------------------------------------------
resultados <- list(
  res0=res0, resA=resA,
  resB=resB, resC=resC, resD=resD,
  beta_k=beta_k, rmse_pct=rmse_pct,
  params=params, x0=x0,
  T_dias=T_dias, tvec=tvec,
  nos_nomes=nos_nomes, grupos=grupos,
  N_kg=N_kg, N_nos=N_nos,
  casos_misau=casos_misau,
  dias_boletins=dias_boletins,
  D0=D0, tab_ace=tab_ace,
  tab_etaria=tab_etaria,
  prop_etaria=prop_etaria,
  I0_nos=I0_nos, R0_nos=R0_nos,
  subtotal_n1=subtotal_n1,
  A_kg_base=A_kg_base,
  omega_g=omega_g
)

saveRDS(resultados,
        "resultados/resultados_meta.rds")
cat("\nResultados guardados em",
    "resultados/resultados_meta.rds\n")