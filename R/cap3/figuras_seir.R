# ================================================
# figuras_seir.R — versão final
# SEIR Moçambique Jan-Mar 2021
# ================================================
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(patchwork)

dir.create("figuras", showWarnings = FALSE)

# ------------------------------------------------
# 1. CARREGAR PARÂMETROS
# ------------------------------------------------
p         <- readRDS("params_finais.rds")
beta0     <- p$beta0
mu        <- p$mu
gamma     <- p$gamma
sigma     <- p$sigma
N_moz     <- p$N_moz
S0        <- p$S0
E0        <- p$E0
I0        <- p$I0
R0_moz    <- p$R0_moz
C_obs_adj <- p$C_obs_adj
datas     <- p$datas
R0_num    <- p$R0_num
T_dias    <- 89
N_fbs     <- T_dias + 1
tvec      <- 0:T_dias
A_peso    <- 1
omega     <- 0.5
tol_fbs   <- 1e-6

# ------------------------------------------------
# 2. INTEGRADORES
# ------------------------------------------------

# SEIR com beta variável no tempo
seir_rk4_beta_var <- function(beta0, mu,
                               gamma, S0, E0,
                               I0, R0, N,
                               T_dias) {
  sigma <- 1/5.2
  n     <- T_dias + 1
  S  <- numeric(n); S[1]  <- S0
  E  <- numeric(n); E[1]  <- E0
  Iv <- numeric(n); Iv[1] <- I0
  R  <- numeric(n); R[1]  <- R0
  C  <- numeric(n); C[1]  <- I0 + R0

  for (i in 1:T_dias) {
    s  <- S[i]; e <- E[i]
    ii <- Iv[i]; r <- R[i]
    bt  <- beta0 * exp(-mu*(i-1))
    bt2 <- beta0 * exp(-mu*(i-0.5))
    bt4 <- beta0 * exp(-mu*i)

    f  <- function(ss,ii) bt *ss*ii/N
    f2 <- function(ss,ii) bt2*ss*ii/N
    f4 <- function(ss,ii) bt4*ss*ii/N

    k1s <- -f(s,ii)
    k1e <-  f(s,ii)  - sigma*e
    k1i <-  sigma*e  - gamma*ii
    k1r <-  gamma*ii

    s2 <- s+0.5*k1s; e2 <- e+0.5*k1e
    i2 <- ii+0.5*k1i
    k2s <- -f2(s2,i2)
    k2e <-  f2(s2,i2) - sigma*e2
    k2i <-  sigma*e2  - gamma*i2
    k2r <-  gamma*i2

    s3 <- s+0.5*k2s; e3 <- e+0.5*k2e
    i3 <- ii+0.5*k2i
    k3s <- -f2(s3,i3)
    k3e <-  f2(s3,i3) - sigma*e3
    k3i <-  sigma*e3  - gamma*i3
    k3r <-  gamma*i3

    s4 <- s+k3s; e4 <- e+k3e; i4 <- ii+k3i
    k4s <- -f4(s4,i4)
    k4e <-  f4(s4,i4) - sigma*e4
    k4i <-  sigma*e4  - gamma*i4
    k4r <-  gamma*i4

    S[i+1]  <- s +(k1s+2*k2s+2*k3s+k4s)/6
    E[i+1]  <- e +(k1e+2*k2e+2*k3e+k4e)/6
    Iv[i+1] <- ii+(k1i+2*k2i+2*k3i+k4i)/6
    R[i+1]  <- r +(k1r+2*k2r+2*k3r+k4r)/6
    C[i+1]  <- Iv[i+1] + R[i+1]
  }
  data.frame(S=S, E=E, I=Iv, R=R, C=C)
}

# SEIR com controlo — beta variável + u(t)
seir_ctrl <- function(beta0, mu, gamma,
                      S0, E0, I0, R0,
                      N, u_vec, T_dias) {
  sigma <- 1/5.2
  n     <- T_dias + 1
  S  <- numeric(n); S[1]  <- S0
  E  <- numeric(n); E[1]  <- E0
  Iv <- numeric(n); Iv[1] <- I0
  R  <- numeric(n); R[1]  <- R0

  for (i in 1:T_dias) {
    s <- S[i]; e <- E[i]
    ii <- Iv[i]; r <- R[i]
    bt <- beta0*exp(-mu*(i-1))*(1-u_vec[i])

    k1s <- -bt*s*ii/N
    k1e <-  bt*s*ii/N  - sigma*e
    k1i <-  sigma*e    - gamma*ii
    k1r <-  gamma*ii

    bt2 <- beta0*exp(-mu*(i-0.5))*
           (1-u_vec[i])
    s2 <- s+0.5*k1s; e2 <- e+0.5*k1e
    i2 <- ii+0.5*k1i
    k2s <- -bt2*s2*i2/N
    k2e <-  bt2*s2*i2/N - sigma*e2
    k2i <-  sigma*e2    - gamma*i2
    k2r <-  gamma*i2

    s3 <- s+0.5*k2s; e3 <- e+0.5*k2e
    i3 <- ii+0.5*k2i
    k3s <- -bt2*s3*i3/N
    k3e <-  bt2*s3*i3/N - sigma*e3
    k3i <-  sigma*e3    - gamma*i3
    k3r <-  gamma*i3

    bt4 <- beta0*exp(-mu*i)*(1-u_vec[i])
    s4 <- s+k3s; e4 <- e+k3e; i4 <- ii+k3i
    k4s <- -bt4*s4*i4/N
    k4e <-  bt4*s4*i4/N - sigma*e4
    k4i <-  sigma*e4    - gamma*i4
    k4r <-  gamma*i4

    S[i+1]  <- s +(k1s+2*k2s+2*k3s+k4s)/6
    E[i+1]  <- e +(k1e+2*k2e+2*k3e+k4e)/6
    Iv[i+1] <- ii+(k1i+2*k2i+2*k3i+k4i)/6
    R[i+1]  <- r +(k1r+2*k2r+2*k3r+k4r)/6
  }
  data.frame(S=S, E=E, I=Iv, R=R)
}

# Passo RK4 adjunta
rk4_adj <- function(q, x, u, beta0,
                    mu, gamma, A, N,
                    dia, h) {
  sigma  <- 1/5.2
  bt     <- beta0*exp(-mu*(dia-1))*(1-u)

  dq_dt <- function(q, x) {
    S <- x[1]; I <- x[3]
    c(
      (q[1]-q[2])*bt*I/N,
      sigma*(q[2]-q[3]),
      -A + (q[1]-q[2])*bt*S/N +
        gamma*(q[3]-q[4]),
      0
    )
  }

  k1 <- dq_dt(q, x)
  k2 <- dq_dt(q+h/2*k1, x)
  k3 <- dq_dt(q+h/2*k2, x)
  k4 <- dq_dt(q+h*k3,   x)
  q + (h/6)*(k1+2*k2+2*k3+k4)
}

# ------------------------------------------------
# 3. FBS — função reutilizável
# ------------------------------------------------
fbs_seir <- function(B_peso) {
  x_m <- matrix(0, N_fbs, 4)
  x_m[1,] <- c(S0, E0, I0, R0_moz)
  q_m <- matrix(0, N_fbs, 4)
  u   <- rep(0, N_fbs)
  err <- Inf; it <- 0

  while (err > tol_fbs && it < 3000) {
    u_old <- u

    # Forward
    for (i in 1:T_dias) {
      x_m[i+1,] <- as.numeric(
        seir_ctrl(
          beta0, mu, gamma,
          x_m[i,1], 0, x_m[i,3],
          x_m[i,4], N_moz,
          c(u[i], 0), 1
        )[2,]
      )
      # corrigir E
      bt <- beta0*exp(-mu*(i-1))*(1-u[i])
      dE <- bt*x_m[i,1]*x_m[i,3]/N_moz -
            sigma*x_m[i,2]
      x_m[i+1, 2] <- x_m[i,2] + dE
    }

    # Backward via tau
    q_tau    <- matrix(0, N_fbs, 4)
    q_tau[1,] <- c(0,0,0,0)

    for (i in 1:(N_fbs-1)) {
      idx   <- N_fbs - i + 1
      q_tau[i+1,] <- rk4_adj(
        q_tau[i,], x_m[idx,],
        u[idx], beta0, mu, gamma,
        A_peso, N_moz, idx, 1
      )
    }

    for (j in 1:4) {
      q_m[,j] <- rev(q_tau[,j])
    }

    # Actualizar u
    u_hat <- (q_m[,2]-q_m[,1]) *
             beta0 *
             exp(-mu*(tvec)) *
             x_m[,1] * x_m[,3] /
             (B_peso * N_moz)
    u_hat <- pmax(pmin(u_hat, 1), 0)
    u     <- omega*u_hat + (1-omega)*u_old
    err   <- max(abs(u-u_old)) /
             (max(abs(u))+1)
    it    <- it + 1
  }

  cat("B =", B_peso, "| iter =", it,
      "| err =", format(err,
      scientific=TRUE, digits=2), "\n")
  list(x=x_m, q=q_m, u=u,
       iter=it, err=err)
}

# ------------------------------------------------
# 4. SIMULAR — sem controlo
# ------------------------------------------------
cat("A simular cenário livre...\n")
sol_livre <- seir_rk4_beta_var(
  beta0, mu, gamma,
  S0, E0, I0, R0_moz,
  N_moz, T_dias
)

# ------------------------------------------------
# 5. SIMULAR — com controlo (B=100)
# ------------------------------------------------
cat("A executar FBS (B=100)...\n")
res_ctrl <- fbs_seir(B_peso = 100)

# ------------------------------------------------
# 6. FIGURA 1 — Ajuste modelo vs dados reais
# ------------------------------------------------
cat("A gerar Figura 1...\n")

df_ajuste <- data.frame(
  dia     = rep(tvec, 2),
  Valor   = c(sol_livre$C,
              as.numeric(C_obs_adj)),
  Tipo    = rep(c("Modelo SEIR",
                  "Dados OWID/MISAU"),
                each = N_fbs)
)

p1 <- ggplot(df_ajuste,
             aes(x = dia, y = Valor/1000,
                 color = Tipo,
                 linetype = Tipo)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = c("Modelo SEIR"      = "#2C3E50",
               "Dados OWID/MISAU" = "#C0392B")
  ) +
  scale_linetype_manual(
    values = c("Modelo SEIR"      = "solid",
               "Dados OWID/MISAU" = "dashed")
  ) +
  annotate("text", x = 70, y = 55,
           label = paste0(
             "RMSE% = 3.92%\n",
             "R\u2080 = ", round(R0_num, 2)),
           color = "#2C3E50",
           size  = 3.5,
           hjust = 0) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Ajuste do Modelo SEIR aos Dados Reais",
    subtitle = "Moçambique — Janeiro a Março 2021",
    x        = "Dias (a partir de 1 Jan 2021)",
    y        = "Casos Acumulados (milhares)",
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("figuras/seir_ajuste.pdf",
       plot = p1, width = 8, height = 5,
       device = cairo_pdf)
cat("Figura 1 guardada.\n")

# ------------------------------------------------
# 7. FIGURA 2 — Compartimentos S, E, I, R
# ------------------------------------------------
cat("A gerar Figura 2...\n")

df_comp <- data.frame(
  dia     = rep(tvec, 8),
  Valor   = c(
    sol_livre$S,     sol_livre$E,
    sol_livre$I,     sol_livre$R,
    res_ctrl$x[,1], res_ctrl$x[,2],
    res_ctrl$x[,3], res_ctrl$x[,4]
  ),
  Comp    = rep(rep(c("S","E","I","R"),
                    each = N_fbs), 2),
  Cenario = rep(
    c("Sem Controlo","Com Controlo"),
    each = 4 * N_fbs)
)

df_comp$Comp <- factor(
  df_comp$Comp,
  levels = c("S","E","I","R"),
  labels = c("Susceptíveis (S)",
             "Expostos (E)",
             "Infectados (I)",
             "Recuperados (R)")
)

p2 <- ggplot(df_comp,
             aes(x = dia,
                 y = Valor / 1000,
                 color    = Cenario,
                 linetype = Cenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Comp,
             scales   = "free_y",
             ncol     = 2) +
  scale_color_manual(
    values = c("Sem Controlo" = "#C0392B",
               "Com Controlo" = "#27AE60")
  ) +
  scale_linetype_manual(
    values = c("Sem Controlo" = "dashed",
               "Com Controlo" = "solid")
  ) +
  theme_minimal(base_size = 11) +
  labs(
    title    = paste0(
      "Dinâmica SEIR — Moçambique ",
      "Jan-Mar 2021"),
    subtitle = paste0(
      "Comparação: evolução natural ",
      "vs controlo óptimo (B=100)"),
    x        = "Dias",
    y        = "Indivíduos (milhares)",
    color    = NULL,
    linetype = NULL
  ) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(
                         face = "bold")
  )

ggsave("figuras/seir_compartimentos.pdf",
       plot = p2, width = 10, height = 7,
       device = cairo_pdf)
cat("Figura 2 guardada.\n")

# ------------------------------------------------
# 8. FIGURA 3 — Controlo óptimo u(t)
# ------------------------------------------------
cat("A gerar Figura 3...\n")

df_u <- data.frame(
  dia = tvec,
  u   = res_ctrl$u
)

p3 <- ggplot(df_u,
             aes(x = dia, y = u)) +
  geom_line(color     = "#2980B9",
            linewidth = 1.2) +
  geom_area(alpha = 0.12,
            fill  = "#2980B9") +
  geom_hline(yintercept = 0,
             linetype   = "dotted",
             color      = "grey50") +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Estratégia de Controlo Óptimo",
    subtitle = paste0(
      "Taxa de mitigação u*(t) ",
      "derivada pelo PMP (B=100)"),
    x        = "Dias",
    y        = "Intensidade u*(t)"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave("figuras/seir_controlo.pdf",
       plot = p3, width = 8, height = 5,
       device = cairo_pdf)
cat("Figura 3 guardada.\n")

# ------------------------------------------------
# 9. FIGURA 4 — R_t com e sem controlo
# ------------------------------------------------
cat("A gerar Figura 4...\n")

Rt_livre <- sapply(1:N_fbs, function(i) {
  bt <- beta0 * exp(-mu*(i-1))
  bt * sol_livre$S[i] / (gamma * N_moz)
})

Rt_ctrl <- sapply(1:N_fbs, function(i) {
  bt <- beta0 * exp(-mu*(i-1)) *
        (1 - res_ctrl$u[i])
  bt * res_ctrl$x[i,1] / (gamma * N_moz)
})

df_rt <- data.frame(
  dia     = rep(tvec, 2),
  Rt      = c(Rt_livre, Rt_ctrl),
  Cenario = rep(
    c("Sem Controlo", "Com Controlo"),
    each = N_fbs)
)

# Dia em que Rt cruza 1 sem controlo
dia_limiar <- which(Rt_livre < 1)[1]

p4 <- ggplot(df_rt,
             aes(x = dia, y = Rt,
                 color    = Cenario,
                 linetype = Cenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1,
             linetype   = "dotted",
             color      = "grey30",
             linewidth  = 0.8) +
  geom_vline(xintercept = dia_limiar,
             linetype   = "dashed",
             color      = "#7F8C8D",
             linewidth  = 0.6) +
  annotate("text",
           x     = dia_limiar + 2,
           y     = 1.8,
           label = paste0("R_t = 1\n",
                          "dia ", dia_limiar),
           color = "#7F8C8D",
           size  = 3.2,
           hjust = 0) +
  annotate("text",
           x     = 75,
           y     = 1.08,
           label = "Limiar R_t = 1",
           color = "grey30",
           size  = 3.2) +
  scale_color_manual(
    values = c("Sem Controlo" = "#C0392B",
               "Com Controlo" = "#27AE60")
  ) +
  scale_linetype_manual(
    values = c("Sem Controlo" = "dashed",
               "Com Controlo" = "solid")
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Número de Reprodução Efectivo",
    subtitle = paste0(
      "R_t cruza o limiar unitário ",
      "no dia ", dia_limiar,
      " sem intervenção"),
    x        = "Dias",
    y        = expression(R[t]),
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("figuras/seir_rt.pdf",
       plot = p4, width = 8, height = 5,
       device = cairo_pdf)
cat("Figura 4 guardada.\n")

# ------------------------------------------------
# 10. FIGURA 5 — Variáveis adjuntas
# ------------------------------------------------
cat("A gerar Figura 5...\n")

df_adj <- data.frame(
  dia      = rep(tvec, 4),
  Valor    = c(res_ctrl$q[,1],
               res_ctrl$q[,2],
               res_ctrl$q[,3],
               res_ctrl$q[,4]),
  Variavel = rep(
    c("q[S](t)", "q[E](t)",
      "q[I](t)", "q[R](t)"),
    each = N_fbs)
)

p5 <- ggplot(df_adj,
             aes(x        = dia,
                 y        = Valor,
                 color    = Variavel,
                 linetype = Variavel)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0,
             linetype   = "dotted",
             color      = "grey50") +
  scale_color_manual(
    values = c(
      "q[S](t)" = "#2C3E50",
      "q[E](t)" = "#E67E22",
      "q[I](t)" = "#C0392B",
      "q[R](t)" = "#27AE60"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "q[S](t)" = "solid",
      "q[E](t)" = "dashed",
      "q[I](t)" = "dotdash",
      "q[R](t)" = "dotted"
    )
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Variáveis Adjuntas do Sistema SEIR",
    subtitle = paste0(
      "Condição de transversalidade: ",
      "q(T) = 0 verificada"),
    x        = "Dias",
    y        = "Valor Sombra (Custo Marginal)",
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("figuras/seir_adjuntas.pdf",
       plot = p5, width = 8, height = 5,
       device = cairo_pdf)
cat("Figura 5 guardada.\n")

# ------------------------------------------------
# 11. FIGURA 6 — Sensibilidade ao parâmetro B
# ------------------------------------------------
cat("A executar análise de sensibilidade...\n")

B_cenarios <- c(
  "B=10 (Prioridade Sanitária)"   = 10,
  "B=100 (Equilíbrio)"            = 100,
  "B=500 (Moderado)"              = 500,
  "B=1000 (Prioridade Económica)" = 1000
)

cores_B <- c(
  "B=10 (Prioridade Sanitária)"   = "#C0392B",
  "B=100 (Equilíbrio)"            = "#2980B9",
  "B=500 (Moderado)"              = "#27AE60",
  "B=1000 (Prioridade Económica)" = "#E67E22"
)

df_B_u  <- data.frame()
df_B_I  <- data.frame()
df_B_Rt <- data.frame()
tab_ace <- data.frame()

for (nm in names(B_cenarios)) {
  B_t  <- B_cenarios[nm]
  res_B <- fbs_seir(B_peso = B_t)

  # Controlo
  df_B_u <- rbind(df_B_u, data.frame(
    dia     = tvec,
    u       = res_B$u,
    Cenario = nm
  ))

  # Infectados
  df_B_I <- rbind(df_B_I, data.frame(
    dia     = tvec,
    I       = res_B$x[,3],
    Cenario = nm
  ))

  # Rt
  Rt_B <- sapply(1:N_fbs, function(i) {
    bt <- beta0*exp(-mu*(i-1))*
          (1-res_B$u[i])
    bt*res_B$x[i,1]/(gamma*N_moz)
  })
  df_B_Rt <- rbind(df_B_Rt, data.frame(
    dia     = tvec,
    Rt      = Rt_B,
    Cenario = nm
  ))

  # Métricas para tabela ACE
  C_ctrl  <- res_B$x[,3] + res_B$x[,4]
  C_livre <- sol_livre$I  + sol_livre$R
  eficacia <- sum(C_livre - C_ctrl)
  custo    <- sum(0.5 * B_t * res_B$u^2)
  acer     <- custo / max(eficacia, 1)
  u_medio  <- mean(res_B$u)
  I_pico   <- round(max(res_B$x[,3]))
  red_pct  <- round(
    (1 - I_pico/round(max(sol_livre$I))
    ) * 100, 1)

  tab_ace <- rbind(tab_ace, data.frame(
    Cenario    = nm,
    B          = B_t,
    Custo      = round(custo, 1),
    Eficacia   = round(eficacia),
    ACER       = round(acer, 6),
    u_medio    = round(u_medio, 3),
    I_pico     = I_pico,
    Reducao_pct = red_pct
  ))
}

# --- 6a. Trajectória do controlo u(t) ---
p6a <- ggplot(df_B_u,
              aes(x        = dia,
                  y        = u,
                  color    = Cenario,
                  linetype = Cenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cores_B) +
  scale_linetype_manual(
    values = c("solid","dashed",
               "dotdash","dotted")) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  theme_minimal(base_size = 11) +
  labs(
    title    = "Trajectória do Controlo u*(t)",
    subtitle = "Por cenário de investimento",
    x        = "Dias",
    y        = "Intensidade u*(t)",
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        legend.text      = element_text(
                             size = 8))

# --- 6b. Infectados por cenário ---
p6b <- ggplot(df_B_I,
              aes(x        = dia,
                  y        = I / 1000,
                  color    = Cenario,
                  linetype = Cenario)) +
  geom_line(data = data.frame(
    dia     = tvec,
    I       = sol_livre$I,
    Cenario = "Sem Controlo"
  ), aes(x=dia, y=I/1000),
  color    = "grey50",
  linetype = "dotted",
  linewidth = 0.8,
  inherit.aes = FALSE) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cores_B) +
  scale_linetype_manual(
    values = c("solid","dashed",
               "dotdash","dotted")) +
  theme_minimal(base_size = 11) +
  labs(
    title    = "Infectados Activos I(t)",
    subtitle = "Impacto por cenário de B",
    x        = "Dias",
    y        = "Infectados (milhares)",
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        legend.text      = element_text(
                             size = 8))

# --- 6c. R_t por cenário ---
p6c <- ggplot(df_B_Rt,
              aes(x        = dia,
                  y        = Rt,
                  color    = Cenario,
                  linetype = Cenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 1,
             linetype   = "dotted",
             color      = "grey30",
             linewidth  = 0.7) +
  annotate("text", x = 75, y = 1.08,
           label = "R_t = 1",
           color = "grey30",
           size  = 3.2) +
  scale_color_manual(values = cores_B) +
  scale_linetype_manual(
    values = c("solid","dashed",
               "dotdash","dotted")) +
  theme_minimal(base_size = 11) +
  labs(
    title    = expression(R[t]*" por Cenário"),
    subtitle = "Efeito do parâmetro B sobre R_t",
    x        = "Dias",
    y        = expression(R[t]),
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        legend.text      = element_text(
                             size = 8))

# --- 6d. Fronteira de eficiência ---
p6d <- ggplot(tab_ace,
              aes(x     = Eficacia / 1000,
                  y     = Custo,
                  color = Cenario,
                  label = paste0("B=", B))) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8,
            hjust =  0.5,
            size  =  3.5) +
  geom_path(color     = "grey50",
            linetype  = "dashed",
            linewidth = 0.6) +
  scale_color_manual(values = cores_B) +
  theme_minimal(base_size = 11) +
  labs(
    title    = "Fronteira de Eficiência",
    subtitle = "Plano Custo-Eficácia (ACER)",
    x        = "Eficácia (casos evitados, milhares)",
    y        = "Custo Total da Intervenção",
    color    = NULL
  ) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank())

# --- Painel combinado 2x2 ---
p6_painel <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(
    title = paste0(
      "Análise de Sensibilidade — ",
      "Parâmetro B"),
    subtitle = paste0(
      "Impacto do peso do custo de ",
      "intervenção na estratégia óptima"),
    theme = theme(
      plot.title    = element_text(
                        face = "bold",
                        size = 13),
      plot.subtitle = element_text(
                        size = 10,
                        color = "grey40")
    )
  )

ggsave("figuras/seir_sensibilidade.pdf",
       plot   = p6_painel,
       width  = 14,
       height = 10,
       device = cairo_pdf)
cat("Figura 6 (painel) guardada.\n")

# --- Tabela ACE no console ---
cat("\n=== TABELA CUSTO-EFICÁCIA ===\n")
cat(sprintf(
  "%-35s %6s %10s %10s %8s %8s\n",
  "Cenário", "B", "Custo",
  "Eficácia", "ACER", "Red%"
))
cat(strrep("-", 80), "\n")
for (i in 1:nrow(tab_ace)) {
  cat(sprintf(
    "%-35s %6d %10.1f %10.0f %8.6f %7.1f%%\n",
    tab_ace$Cenario[i],
    tab_ace$B[i],
    tab_ace$Custo[i],
    tab_ace$Eficacia[i],
    tab_ace$ACER[i],
    tab_ace$Reducao_pct[i]
  ))
}

# ------------------------------------------------
# 12. CONFIRMAÇÃO FINAL
# ------------------------------------------------
cat("\n=== FICHEIROS GERADOS ===\n")
print(list.files("figuras/")
  )