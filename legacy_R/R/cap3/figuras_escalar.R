# ================================================
# figuras_escalar.R — versão 4 (backward via tau)
# ================================================
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

dir.create("figuras", showWarnings = FALSE)

# ------------------------------------------------
# 1. PARÂMETROS
# ------------------------------------------------
T_f  <- 1.0
dt   <- 0.01
tvec <- seq(0, T_f, by = dt)
N    <- length(tvec)
x0   <- 0.5

# ------------------------------------------------
# 2. INTEGRADOR RK4 — passo simples
# ------------------------------------------------
rk4_step <- function(y, fn, t_n, h) {
  k1 <- fn(y,             t_n)
  k2 <- fn(y + h/2 * k1, t_n + h/2)
  k3 <- fn(y + h/2 * k2, t_n + h/2)
  k4 <- fn(y + h   * k3, t_n + h)
  y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
}

# ------------------------------------------------
# 3. SOLUÇÃO ANALÍTICA
# ------------------------------------------------
A_mat <- matrix(
  c(1, 1,
    (1 - sqrt(3)) * exp(sqrt(3)),
    (1 + sqrt(3)) * exp(-sqrt(3))),
  nrow = 2, byrow = TRUE
)
coefs       <- solve(A_mat, c(0.5, 0))
c1          <- coefs[1]
c2          <- coefs[2]
x_analitico <- c1 * exp(sqrt(3) * tvec) +
               c2 * exp(-sqrt(3) * tvec)

cat("=== SOLUÇÃO ANALÍTICA ===\n")
cat("c1 =", round(c1, 7), "\n")
cat("c2 =", round(c2, 7), "\n")
cat("x(0) =", round(x_analitico[1], 6),
    "— deve ser 0.5\n")

# ------------------------------------------------
# 4. LADO DIREITO DAS EDOs
# ------------------------------------------------

# Forward: dx/dt = x + u
dx_dt <- function(x, u) x + u

# Backward via mudança de variável tau = T_f - t
# dq/dtau = +(2*x(T_f - tau) + q)
# x_tau é o valor de x no instante T_f - tau
dq_dtau <- function(q, x_tau) 2 * x_tau + q

# ------------------------------------------------
# 5. FBS COM RK4 — abordagem tau
# ------------------------------------------------
x_num <- numeric(N); x_num[1] <- x0
q_num <- numeric(N)
u     <- rep(0, N)

tolerancia <- 1e-7
omega      <- 0.5
erro_fbs   <- Inf
iter       <- 0

while (erro_fbs > tolerancia && iter < 2000) {
  u_old <- u

  # --- Passo Forward: t de 0 a T_f ---
  for (i in 1:(N - 1)) {
    ui         <- u[i]
    x_num[i+1] <- rk4_step(
      x_num[i],
      function(y, s) dx_dt(y, ui),
      tvec[i], dt
    )
  }

  # --- Passo Backward via tau = T_f - t ---
  # tau=0 corresponde a t=T_f: q(tau=0) = q(T_f) = 0
  # tau=T_f corresponde a t=0
  # x no instante t = T_f - tau é x_num[N - i]
  q_tau    <- numeric(N)
  q_tau[1] <- 0   # condição terminal: q(T_f) = 0

  for (i in 1:(N - 1)) {
    # índice do instante t = T_f - tau_i
    # tau_i = (i-1)*dt → t = T_f - (i-1)*dt → índice N-i+1
    idx_t      <- N - i + 1
    x_tau_i    <- x_num[idx_t]

    q_tau[i+1] <- rk4_step(
      q_tau[i],
      function(y, s) dq_dtau(y, x_tau_i),
      (i-1) * dt, dt
    )
  }

  # Converter de volta: q(t) = q_tau(T_f - t)
  # tau = T_f - t → índice em tau = N - índice em t + 1
  q_num <- rev(q_tau)

  # --- Controlo óptimo com projecção ---
  u_hat <- pmax(pmin(-q_num, 1), -1)
  u     <- omega * u_hat + (1 - omega) * u_old

  # --- Critério L-infinito ---
  erro_fbs <- max(abs(u - u_old)) / (max(abs(u)) + 1)
  iter     <- iter + 1
}

cat("\n=== VERIFICAÇÃO FBS ===\n")
cat("Convergiu em", iter, "iterações.\n")
cat("Erro FBS:   ", format(erro_fbs, scientific=TRUE, digits=3), "\n")
cat("x_num[1]   =", round(x_num[1], 6), "— deve ser 0.5\n")
cat("x_num[N]   =", round(x_num[N], 6), "\n")
cat("q_num[N]   =", round(q_num[N], 8), "— deve ser ~0\n")
cat("u min      =", round(min(u), 6), "— deve ser >= -1\n")
cat("u max      =", round(max(u), 6), "— deve ser <= 1\n")

cat("\n=== ERRO RESIDUAL ===\n")
erro_abs    <- abs(x_num - x_analitico)
erro_maximo <- max(erro_abs)
erro_medio  <- mean(erro_abs)
cat("Erro máximo:", format(erro_maximo, scientific=TRUE, digits=3), "\n")
cat("Erro médio: ", format(erro_medio,  scientific=TRUE, digits=3), "\n")

# ------------------------------------------------
# 6. FIGURA 1 — Estado: FBS vs Analítico
# ------------------------------------------------
df_estado <- data.frame(
  tvec   = rep(tvec, 2),
  Valor  = c(x_num, x_analitico),
  Metodo = rep(c("Numérico (FBS + RK4)",
                 "Analítico"), each = N)
)

p1 <- ggplot(df_estado,
             aes(x = tvec, y = Valor,
                 color = Metodo, linetype = Metodo)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Analítico"            = "#C0392B",
               "Numérico (FBS + RK4)" = "#2C3E50")
  ) +
  scale_linetype_manual(
    values = c("Analítico"            = "dashed",
               "Numérico (FBS + RK4)" = "solid")
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Trajectória de Estado: Analítica vs. Numérica",
    subtitle = "Minimização de J com restrições no controlo",
    x        = "Tempo (t)",
    y        = "Estado x(t)",
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("figuras/escalar_estado.pdf",
       plot = p1, width = 8, height = 5,
       device = cairo_pdf)
cat("Figura 1 guardada.\n")

# ------------------------------------------------
# 7. FIGURA 2 — Erro residual
# ------------------------------------------------
df_erro <- data.frame(
  tvec = tvec,
  Erro = erro_abs + 1e-16
)

p2 <- ggplot(df_erro, aes(x = tvec, y = Erro)) +
  geom_line(color = "#C0392B", linewidth = 0.8) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Erro Residual: |x_analitico - x_FBS|",
    subtitle = "Escala logarítmica — precisão do integrador RK4",
    x        = "Tempo (t)",
    y        = "Erro Absoluto"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave("figuras/escalar_erro.pdf",
       plot = p2, width = 8, height = 4,
       device = cairo_pdf)
cat("Figura 2 guardada.\n")

# ------------------------------------------------
# 8. FIGURA 3 — Controlo óptimo e variável adjunta
# ------------------------------------------------
df_adjunta <- data.frame(
  tvec     = rep(tvec, 2),
  Valor    = c(u, q_num),
  Variavel = rep(c("u(t)", "q(t)"), each = N)
)

p3 <- ggplot(df_adjunta,
             aes(x = tvec, y = Valor,
                 color    = Variavel,
                 linetype = Variavel)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0,
             linetype = "dotted",
             color    = "grey50") +
  scale_color_manual(
    values = c("u(t)" = "#E67E22",
               "q(t)" = "#2980B9"),
    labels = c("u*(t)", "q(t)")
  ) +
  scale_linetype_manual(
    values = c("u(t)" = "solid",
               "q(t)" = "dashed"),
    labels = c("u*(t)", "q(t)")
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Controlo Optimo e Variavel Adjunta",
    subtitle = "Relacao u*(t) = -q(t) com saturacao em [-1, 1]",
    x        = "Tempo (t)",
    y        = "Amplitude",
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("figuras/escalar_adjunta.pdf",
       plot = p3, width = 8, height = 5,
       device = cairo_pdf)
cat("Figura 3 guardada.\n")

# ------------------------------------------------
# 9. FIGURA 4 — Sistema controlado vs livre
# ------------------------------------------------
x_livre <- 0.5 * exp(tvec)

df_comp <- data.frame(
  tvec    = rep(tvec, 2),
  Valor   = c(x_num, x_livre),
  Cenario = rep(c("Com Controlo Optimo",
                  "Sem Controlo (u=0)"), each = N)
)

p4 <- ggplot(df_comp,
             aes(x = tvec, y = Valor,
                 color    = Cenario,
                 linetype = Cenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Com Controlo Optimo" = "#27AE60",
               "Sem Controlo (u=0)"  = "#C0392B")
  ) +
  scale_linetype_manual(
    values = c("Com Controlo Optimo" = "solid",
               "Sem Controlo (u=0)"  = "dotted")
  ) +
  annotate("text",
           x     = 0.78,
           y     = 1.18,
           label = "Crescimento\nExponencial",
           color = "#C0392B",
           size  = 3.5) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Sistema Controlado vs. Sistema Livre",
    subtitle = "Esforco de minimizacao de J(u)",
    x        = "Tempo (t)",
    y        = "Estado x(t)",
    color    = NULL,
    linetype = NULL
  ) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave("figuras/escalar_impacto.pdf",
       plot = p4, width = 8, height = 5,
       device = cairo_pdf)
cat("Figura 4 guardada.\n")

cat("\nFicheiros em figuras/:\n")
print(list.files("figuras/"))