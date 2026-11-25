# ================================================
# sensibilidade_R0.R — Versão final
# Dois mapas de calor:
# A — impacto absoluto por nó
# B — sensibilidade por faixa etária
# ================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

dir.create("Figuras", showWarnings = FALSE)

# ------------------------------------------------
# 1. PARÂMETROS BASE
# ------------------------------------------------
K     <- 6
G     <- 9
gamma <- 0.08
sigma <- 1/5.2
delta <- 0.50
e_vac <- 0.70
p_x   <- 0.23
T_V   <- 14

nos_nomes <- c("Maputo","Xai-Xai",
               "Maxixe","Caia",
               "Inchope","Namacurra")

grupos <- c("0-4","5-14","15-24",
            "25-34","35-44","45-54",
            "55-64","65-74","75+")

beta_k <- c(0.1484, 0.1637, 0.2487,
            0.1400, 0.1943, 0.2249)

# Parâmetros clínicos por faixa etária
p_h_g <- c(0.010,0.005,0.010,0.020,
           0.040,0.080,0.130,0.200,0.300)
p_c_g <- c(0.050,0.030,0.040,0.060,
           0.080,0.120,0.180,0.250,0.350)
mu_g  <- c(0.0002,0.0001,0.0002,0.0005,
           0.0010,0.0030,0.0080,0.0200,
           0.0500)
omega_g <- c(0.5,0.3,0.4,0.6,0.8,
             1.0,1.5,3.0,5.0)
kappa_k <- c(1.0,1.2,1.4,1.8,1.6,1.5)

# ------------------------------------------------
# 2. FUNÇÃO R0
# ------------------------------------------------
calc_R0 <- function(beta, gamma,
                    delta, e, p_x) {
  S_eff <- p_x +
           (1-p_x) +
           (1-p_x)*(1-e)*delta
  beta / gamma * S_eff
}

# R0 base por nó
R0_base <- sapply(beta_k, function(b)
  calc_R0(b, gamma, delta, e_vac, p_x))

cat("=== R0 BASE POR NÓ ===\n")
for (k in 1:K) {
  cat(sprintf("  %-12s R0=%.3f\n",
      nos_nomes[k], R0_base[k]))
}

# ------------------------------------------------
# 3A. MAPA A — IMPACTO ABSOLUTO
#     Variação de R0 face a
#     perturbação de 10% em
#     cada parâmetro, por nó
# ------------------------------------------------
cat("\nA calcular Mapa A...\n")

perturb <- 0.10  # 10%

params_A <- list(
  "beta[k]"   = list(
    nodal=TRUE,
    fn=function(b,g,d,e,px,h)
      calc_R0(b*(1+h),g,d,e,px)),
  "gamma"     = list(
    nodal=FALSE,
    fn=function(b,g,d,e,px,h)
      calc_R0(b,g*(1+h),d,e,px)),
  "italic(e)" = list(
    nodal=FALSE,
    fn=function(b,g,d,e,px,h)
      calc_R0(b,g,d,e*(1+h),px)),
  "delta"     = list(
    nodal=FALSE,
    fn=function(b,g,d,e,px,h)
      calc_R0(b,g,d*(1+h),e,px)),
  "T[V]"      = list(
    nodal=FALSE,
    fn=function(b,g,d,e,px,h)
      calc_R0(b,g,d,e*(1-h*0.1),px)),
  "p[x]"      = list(
    nodal=FALSE,
    fn=function(b,g,d,e,px,h)
      calc_R0(b,g,d,e,px*(1+h)))
)

res_A <- data.frame()
for (nome_p in names(params_A)) {
  par <- params_A[[nome_p]]
  for (k in 1:K) {
    R0_p <- par$fn(
      beta_k[k], gamma, delta,
      e_vac, p_x, perturb)
    delta_R0 <- R0_p - R0_base[k]
    res_A <- rbind(res_A,
      data.frame(
        No       = nos_nomes[k],
        Param    = nome_p,
        Delta_R0 = delta_R0
      ))
  }
}

# Ordenar por impacto médio absoluto
ordem_A <- res_A |>
  group_by(Param) |>
  summarise(m=mean(abs(Delta_R0))) |>
  arrange(desc(m)) |>
  pull(Param)

res_A$Param <- factor(res_A$Param,
                      levels=ordem_A)
res_A$No    <- factor(res_A$No,
                      levels=nos_nomes)

lim_A <- max(abs(res_A$Delta_R0))*1.05

pA <- ggplot(res_A,
             aes(x=Param,
                 y=No,
                 fill=Delta_R0)) +
  geom_tile(color="white",
            linewidth=0.7) +
  geom_text(
    aes(label=sprintf("%+.3f",
                      Delta_R0),
        color=ifelse(
          abs(Delta_R0)>lim_A*0.5,
          "white","grey20")),
    size=3.2, fontface="bold") +
  scale_color_identity() +
  scale_fill_gradient2(
    low="#2980B9", mid="white",
    high="#E74C3C", midpoint=0,
    limits=c(-lim_A, lim_A),
    name=expression(
      Delta*R[0]*
      " (perturb. 10%)")
  ) +
  scale_x_discrete(labels=c(
    "beta[k]"   =expression(beta[k]),
    "gamma"     =expression(gamma),
    "italic(e)" =expression(italic(e)),
    "delta"     =expression(delta),
    "T[V]"      =expression(T[V]),
    "p[x]"     =expression(p[x]))
  ) +
  theme_minimal(base_size=10) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(
      size=11, color="#2C3E50"),
    axis.text.y     = element_text(
      size=9, color="#2C3E50"),
    legend.position = "bottom",
    legend.key.width= unit(2,"cm"),
    legend.key.height=unit(0.4,"cm"),
    plot.title      = element_text(
      face="bold", size=10,
      color="#2C3E50"),
    plot.subtitle   = element_text(
      size=8, color="grey40")
  ) +
  labs(
    title    = paste0(
      "(a) Variação absoluta de ",
      "\u211c\u2080 por nó"),
    subtitle = "Perturbação de 10% em cada parâmetro",
    x = "Parâmetro",
    y = "Nó"
  )

# ------------------------------------------------
# 3B. MAPA B — SENSIBILIDADE
#     POR FAIXA ETÁRIA
#     Elasticidade de R0_efectivo
#     aos parâmetros clínicos g
# ------------------------------------------------
cat("A calcular Mapa B...\n")

# R0 efectivo por faixa etária:
# R0_g = R0_base × (1 - p_h_g×IFR_g)
# onde IFR_g ∝ mu_g × omega_g
# Sensibilidade ao parâmetro clínico

# Parâmetros clínicos por faixa
params_B <- list(
  "p[h]^g"   = p_h_g,
  "p[c]^g"   = p_c_g,
  "mu[g]"    = mu_g,
  "omega[g]" = omega_g
)

# R0 efectivo por faixa:
# IFR_g = mu_g * (p_h_g + p_c_g)
# Mortalidade ponderada pelo peso
calc_IFR <- function(mu, p_h, p_c,
                     omega, kappa) {
  mu * (p_h + p_c) * omega * kappa
}

res_B <- data.frame()
for (nome_p in names(params_B)) {
  vals <- params_B[[nome_p]]
  for (g in 1:G) {
    # IFR base grupo g (média nós)
    IFR_b <- mean(sapply(1:K,
      function(k)
        calc_IFR(mu_g[g], p_h_g[g],
                 p_c_g[g], omega_g[g],
                 kappa_k[k])))

    # IFR perturbado
    vals_p <- vals
    vals_p[g] <- vals[g] * 1.10

    IFR_p <- mean(sapply(1:K,
      function(k) {
        mu_  <- if(nome_p=="mu[g]")
                  vals_p[g] else mu_g[g]
        ph_  <- if(nome_p=="p[h]^g")
                  vals_p[g] else p_h_g[g]
        pc_  <- if(nome_p=="p[c]^g")
                  vals_p[g] else p_c_g[g]
        om_  <- if(nome_p=="omega[g]")
                  vals_p[g] else omega_g[g]
        calc_IFR(mu_, ph_, pc_,
                 om_, kappa_k[k])
      }))

    # Elasticidade
    elast <- (IFR_p - IFR_b) /
             IFR_b / 0.10

    res_B <- rbind(res_B,
      data.frame(
        Grupo    = grupos[g],
        Param    = nome_p,
        Elast    = elast
      ))
  }
}

ordem_B <- res_B |>
  group_by(Param) |>
  summarise(m=mean(abs(Elast))) |>
  arrange(desc(m)) |>
  pull(Param)

res_B$Param <- factor(res_B$Param,
                      levels=ordem_B)
res_B$Grupo <- factor(res_B$Grupo,
                      levels=grupos)

lim_B <- max(abs(res_B$Elast))*1.05

pB <- ggplot(res_B,
             aes(x=Param,
                 y=Grupo,
                 fill=Elast)) +
  geom_tile(color="white",
            linewidth=0.7) +
  geom_text(
    aes(label=sprintf("%.2f", Elast),
        color=ifelse(
          abs(Elast)>lim_B*0.5,
          "white","grey20")),
    size=3.2, fontface="bold") +
  scale_color_identity() +
  scale_fill_gradient2(
    low="#2980B9", mid="white",
    high="#E74C3C", midpoint=0,
    limits=c(-lim_B, lim_B),
    name="Elasticidade IFR"
  ) +
  scale_x_discrete(labels=c(
    "p[h]^g"   =expression(p[h]^g),
    "p[c]^g"   =expression(p[c]^g),
    "mu[g]"    =expression(mu[g]),
    "omega[g]" =expression(omega[g]))
  ) +
  theme_minimal(base_size=10) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(
      size=11, color="#2C3E50"),
    axis.text.y     = element_text(
      size=9, color="#2C3E50"),
    legend.position = "bottom",
    legend.key.width= unit(2,"cm"),
    legend.key.height=unit(0.4,"cm"),
    plot.title      = element_text(
      face="bold", size=10,
      color="#2C3E50"),
    plot.subtitle   = element_text(
      size=8, color="grey40")
  ) +
  labs(
    title    = paste0(
      "(b) Elasticidade da mortalidade ",
      "IFR por faixa etária"),
    subtitle = paste0(
      "Sensibilidade aos parâmetros ",
      "clínicos por grupo g"),
    x = "Parâmetro clínico",
    y = "Faixa etária"
  )

# ------------------------------------------------
# 4. COMBINAR E GUARDAR
# ------------------------------------------------
cat("A combinar painéis...\n")

mapa_final <- (pA | pB) +
  plot_annotation(
    title    = paste0(
      "Análise de Sensibilidade de ",
      "\u211c\u2080 e IFR aos ",
      "Parâmetros do Modelo"),
    subtitle = paste0(
      "Painel (a): variação absoluta ",
      "de \u211c\u2080 por nó | ",
      "Painel (b): elasticidade da ",
      "mortalidade por faixa etária"),
    theme = theme(
      plot.title = element_text(
        face="bold", size=11,
        color="#2C3E50"),
      plot.subtitle = element_text(
        size=8.5, color="grey40"))
  )

ggsave(
  "Figuras/sensibilidade_R0.pdf",
  plot   = mapa_final,
  width  = 16,
  height = 9,
  units  = "cm",
  device = cairo_pdf,
  dpi    = 300
)

cat("Mapa guardado em",
    "Figuras/sensibilidade_R0.pdf\n")