# ================================================
# figuras_meta.R — Figuras do Capítulo 5
# Versão com diferenciação etária
# ================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(ggrepel)

dir.create("Figuras", showWarnings = FALSE)

# ------------------------------------------------
# 1. CARREGAR RESULTADOS
# ------------------------------------------------
cat("A carregar resultados...\n")

r   <- readRDS("resultados/resultados_meta.rds")
res0 <- r$res0; resA <- r$resA
resB <- r$resB; resC <- r$resC
resD <- r$resD
beta_k      <- r$beta_k
nos_nomes   <- r$nos_nomes
grupos      <- r$grupos
T_dias      <- r$T_dias
tvec        <- r$tvec
N_kg        <- r$N_kg
N_nos       <- r$N_nos
params      <- r$params
casos_misau <- r$casos_misau
dias_boletins <- r$dias_boletins
D0          <- r$D0
tab_ace     <- r$tab_ace
tab_etaria  <- r$tab_etaria
prop_etaria <- r$prop_etaria
I0_nos      <- r$I0_nos
R0_nos      <- r$R0_nos
subtotal_n1 <- r$subtotal_n1
omega_g     <- r$omega_g
K <- params$K; G <- params$G

sigma    <- params$sigma
gamma    <- params$gamma
mu_decay <- 0.017

# ------------------------------------------------
# 2. CONFIGURAÇÃO VISUAL
# ------------------------------------------------
cores_cen <- c(
  "0-Base"               = "grey50",
  "A-Prior.Etaria"       = "#E74C3C",
  "B-Prior.Regional"     = "#E67E22",
  "C-Ctrl.Optimo"        = "#27AE60",
  "D-So.Distanciamento"  = "#2980B9"
)

linhas_cen <- c(
  "0-Base"               = "dotted",
  "A-Prior.Etaria"       = "dashed",
  "B-Prior.Regional"     = "dotdash",
  "C-Ctrl.Optimo"        = "solid",
  "D-So.Distanciamento"  = "longdash"
)

cores_grupos <- c(
  "0-4"   = "#2ECC71",
  "5-14"  = "#27AE60",
  "15-24" = "#F1C40F",
  "25-34" = "#E67E22",
  "35-44" = "#E74C3C",
  "45-54" = "#C0392B",
  "55-64" = "#9B59B6",
  "65-74" = "#2980B9",
  "75+"   = "#1A5276"
)

tema_base <- theme_minimal(base_size=11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    plot.title       = element_text(
      face="bold", size=11),
    plot.subtitle    = element_text(
      size=9, color="grey40"),
    strip.text       = element_text(
      face="bold", size=9)
  )

# ------------------------------------------------
# 3. FIG 1 — Ajuste calibração por nó
# ------------------------------------------------
cat("A gerar Figura 1 — Ajuste...\n")

seir_traj <- function(beta0, S0, I0,
                      R0, N, T_dias) {
  n    <- T_dias + 1
  S    <- numeric(n); S[1]  <- S0
  E    <- numeric(n); E[1]  <- round(I0*1.5)
  Iv   <- numeric(n); Iv[1] <- I0
  R    <- numeric(n); R[1]  <- R0
  C    <- numeric(n); C[1]  <- I0 + R0
  for (i in 1:T_dias) {
    bt     <- beta0*exp(-mu_decay*(i-1))
    f      <- bt*S[i]*Iv[i]/N
    S[i+1] <- max(0, S[i]-f)
    E[i+1] <- max(0, E[i]+f-sigma*E[i])
    Iv[i+1]<- max(0, Iv[i]+
                  sigma*E[i]-gamma*Iv[i])
    R[i+1] <- max(0, R[i]+gamma*Iv[i])
    C[i+1] <- Iv[i+1]+R[i+1]
  }
  C
}

df_ajuste <- data.frame()
for (k in 1:K) {
  S0_k  <- N_nos[k]-I0_nos[k]-R0_nos[k]
  C_sim <- seir_traj(
    beta_k[k], S0_k, I0_nos[k],
    R0_nos[k], N_nos[k], T_dias)
  df_ajuste <- rbind(df_ajuste,
    data.frame(dia=tvec, C_sim=C_sim,
               No=nos_nomes[k]))
}

df_obs <- data.frame()
for (k in 1:K) {
  df_obs <- rbind(df_obs,
    data.frame(dia=dias_boletins,
               C_obs=casos_misau[k,],
               No=nos_nomes[k]))
}

df_ajuste$No <- factor(df_ajuste$No,
                       levels=nos_nomes)
df_obs$No    <- factor(df_obs$No,
                       levels=nos_nomes)

p1 <- ggplot() +
  geom_line(data=df_ajuste,
            aes(x=dia, y=C_sim/1000),
            color="#2C3E50",
            linewidth=0.9) +
  geom_point(data=df_obs,
             aes(x=dia, y=C_obs/1000),
             color="#E74C3C",
             size=3, shape=16) +
  facet_wrap(~No, scales="free_y",
             ncol=3) +
  tema_base +
  labs(
    title    = paste0(
      "Ajuste do Modelo SEIR ",
      "aos Dados MISAU por Nó"),
    subtitle = paste0(
      "Linha: modelo calibrado | ",
      "Pontos: Boletins 290, 323 e 354 | ",
      "RMSE% global = ",
      round(mean(r$rmse_pct),1),"%"),
    x = "Dias (a partir de 1 Jan 2021)",
    y = "Casos acumulados (milhares)"
  )

ggsave("Figuras/calibracao_nos.pdf",
       p1, width=10, height=7,
       device=cairo_pdf)
cat("  Figura 1 guardada.\n")

# ------------------------------------------------
# 4. FIG 2 — Cenário Base por nó
# ------------------------------------------------
cat("A gerar Figura 2 — Cenário Base...\n")

df_base <- data.frame()
for (k in 1:K) {
  df_base <- rbind(df_base,
    data.frame(
      dia = tvec,
      I   = res0$I_nos[k,]/1000,
      D   = res0$D_nos[k,],
      No  = nos_nomes[k]))
}
df_base$No <- factor(df_base$No,
                     levels=nos_nomes)

p2a <- ggplot(df_base,
              aes(x=dia, y=I,
                  color=No)) +
  geom_line(linewidth=0.9) +
  scale_color_brewer(palette="Dark2") +
  tema_base +
  labs(title="Infectados activos I(t)",
       x="Dias", y="Infectados (mil)")

p2b <- ggplot(df_base,
              aes(x=dia, y=D,
                  color=No)) +
  geom_line(linewidth=0.9) +
  scale_color_brewer(palette="Dark2") +
  tema_base +
  labs(title="Óbitos acumulados D(t)",
       x="Dias", y="Óbitos")

p2 <- (p2a/p2b) +
  plot_annotation(
    title    = "Cenário 0 — Sem Intervenção",
    subtitle = paste0(
      "Corredor N1, Jan-Mar 2021. ",
      "Maputo lidera o pico com ",
      "desfasamento de ~10 dias ",
      "para Namacurra."),
    theme=theme(
      plot.title=element_text(
        face="bold",size=12))
  )

ggsave("Figuras/cenario_base.pdf",
       p2, width=10, height=8,
       device=cairo_pdf)
cat("  Figura 2 guardada.\n")

# ------------------------------------------------
# 5. FIG 3 — Infectados: todos os cenários
# ------------------------------------------------
cat("A gerar Figura 3 — Comparação I(t)...\n")

df_It <- data.frame(
  dia     = rep(tvec, 5),
  I       = c(res0$I_t, resA$I_t,
              resB$I_t, resC$I_t,
              resD$I_t)/1000,
  Cenario = rep(
    c("0-Base","A-Prior.Etaria",
      "B-Prior.Regional",
      "C-Ctrl.Optimo",
      "D-So.Distanciamento"),
    each=T_dias+1)
)
df_It$Cenario <- factor(df_It$Cenario,
                        levels=names(cores_cen))

p3 <- ggplot(df_It,
             aes(x=dia, y=I,
                 color=Cenario,
                 linetype=Cenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=cores_cen) +
  scale_linetype_manual(values=linhas_cen) +
  tema_base +
  labs(
    title    = paste0(
      "Infectados Activos — ",
      "Comparação de Cenários"),
    subtitle = paste0(
      "Corredor N1, Jan-Mar 2021. ",
      "Cenário A reduz pico em ",
      round(tab_ace$Red_pct[1],1),"%."),
    x = "Dias",
    y = "Infectados activos (mil)"
  )

ggsave("Figuras/cenario_comparacao_I.pdf",
       p3, width=9, height=5,
       device=cairo_pdf)
cat("  Figura 3 guardada.\n")

# ------------------------------------------------
# 6. FIG 4 — Mortalidade acumulada
# ------------------------------------------------
cat("A gerar Figura 4 — Mortalidade...\n")

df_Dt <- data.frame(
  dia     = rep(tvec, 5),
  D       = c(res0$D_t, resA$D_t,
              resB$D_t, resC$D_t,
              resD$D_t),
  Cenario = rep(
    c("0-Base","A-Prior.Etaria",
      "B-Prior.Regional",
      "C-Ctrl.Optimo",
      "D-So.Distanciamento"),
    each=T_dias+1)
)
df_Dt$Cenario <- factor(df_Dt$Cenario,
                        levels=names(cores_cen))

p4 <- ggplot(df_Dt,
             aes(x=dia, y=D,
                 color=Cenario,
                 linetype=Cenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=cores_cen) +
  scale_linetype_manual(values=linhas_cen) +
  tema_base +
  labs(
    title    = paste0(
      "Mortalidade Acumulada — ",
      "Comparação de Cenários"),
    subtitle = paste0(
      "Cenário A evita ",
      round(tab_ace$Red_pct[1],1),
      "% dos óbitos relativamente ",
      "ao cenário base."),
    x = "Dias",
    y = "Óbitos acumulados"
  )

ggsave("Figuras/mortalidade_comparacao.pdf",
       p4, width=9, height=5,
       device=cairo_pdf)
cat("  Figura 4 guardada.\n")

# ------------------------------------------------
# 7. FIG 5 — Controlos óptimos Cenário C
# ------------------------------------------------
cat("A gerar Figura 5 — Controlos...\n")

# u*(t) e v*(t) por nó e faixa etária
# Mostrar variação por grupo etário
df_u_grupo <- data.frame()
df_v_grupo <- data.frame()
for (g in 1:G) {
  df_u_grupo <- rbind(df_u_grupo,
    data.frame(
      dia    = tvec,
      u      = rep(mean(resC$u_kg[,g]),
                   T_dias+1),
      Faixa  = grupos[g]))
  df_v_grupo <- rbind(df_v_grupo,
    data.frame(
      dia    = tvec,
      v      = rep(mean(resC$v_kg[,g]),
                   T_dias+1),
      Faixa  = grupos[g]))
}
df_u_grupo$Faixa <- factor(
  df_u_grupo$Faixa, levels=grupos)
df_v_grupo$Faixa <- factor(
  df_v_grupo$Faixa, levels=grupos)

p5a <- ggplot(df_u_grupo,
              aes(x=dia, y=u,
                  color=Faixa)) +
  geom_hline(yintercept=0.60,
             linetype="dotted",
             color="grey60",
             linewidth=0.6) +
  geom_line(linewidth=0.9) +
  scale_color_manual(
    values=cores_grupos) +
  scale_y_continuous(
    labels=scales::percent,
    limits=c(0,0.65)) +
  annotate("text", x=80, y=0.62,
           label="u_max=60%",
           size=3, color="grey50") +
  tema_base +
  labs(
    title = "Distanciamento u*(t) por faixa",
    x = "Dias", y = "Intensidade u*(t)")

p5b <- ggplot(df_v_grupo,
              aes(x=dia, y=v,
                  color=Faixa)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(
    values=cores_grupos) +
  tema_base +
  labs(
    title = "Vacinação v*(t) por faixa",
    x = "Dias", y = "Taxa v*(t)")

p5 <- (p5a|p5b) +
  plot_annotation(
    title    = paste0(
      "Cenário C — Controlos Óptimos ",
      "por Faixa Etária"),
    subtitle = paste0(
      "FBS diferenciado por (k,g). ",
      "Grupos 65+ recebem maior ",
      "intensidade de vacinação."),
    theme=theme(
      plot.title=element_text(
        face="bold",size=12))
  )

ggsave("Figuras/cenario_C_controlos.pdf",
       p5, width=12, height=5,
       device=cairo_pdf)
cat("  Figura 5 guardada.\n")

# ------------------------------------------------
# 8. FIG 6 — Compartimentos C vs Base
# ------------------------------------------------
cat("A gerar Figura 6 — Compartimentos...\n")

df_comp <- data.frame(
  dia     = rep(tvec, 4),
  Valor   = c(
    res0$I_t/1000, resC$I_t/1000,
    res0$D_t,      resC$D_t),
  Tipo    = rep(
    c("Infectados (mil)",
      "Infectados (mil)",
      "Óbitos acumulados",
      "Óbitos acumulados"),
    each=T_dias+1),
  Cenario = rep(
    c("0-Base","C-Ctrl.Optimo",
      "0-Base","C-Ctrl.Optimo"),
    each=T_dias+1)
)
df_comp$Cenario <- factor(
  df_comp$Cenario,
  levels=c("0-Base","C-Ctrl.Optimo"))
df_comp$Tipo <- factor(
  df_comp$Tipo,
  levels=c("Infectados (mil)",
           "Óbitos acumulados"))

p6 <- ggplot(df_comp,
             aes(x=dia, y=Valor,
                 color=Cenario,
                 linetype=Cenario)) +
  geom_line(linewidth=1.0) +
  facet_wrap(~Tipo, scales="free_y") +
  scale_color_manual(
    values=c(
      "0-Base"        = "grey50",
      "C-Ctrl.Optimo" = "#27AE60")) +
  scale_linetype_manual(
    values=c(
      "0-Base"        = "dashed",
      "C-Ctrl.Optimo" = "solid")) +
  tema_base +
  labs(
    title    = paste0(
      "Cenário C vs Base — ",
      "Infectados e Óbitos"),
    subtitle = paste0(
      "Controlo óptimo simultâneo ",
      "reduz mortalidade em ",
      round(tab_ace$Red_pct[3],1),"%."),
    x = "Dias", y = "Valor"
  )

ggsave(
  "Figuras/cenario_C_compartimentos.pdf",
  p6, width=10, height=5,
  device=cairo_pdf)
cat("  Figura 6 guardada.\n")

# ------------------------------------------------
# 9. FIG 7 — R_t por cenário
# ------------------------------------------------
cat("A gerar Figura 7 — Rt...\n")

df_Rt <- data.frame()
for (res in list(res0,resA,resB,
                 resC,resD)) {
  df_Rt <- rbind(df_Rt, data.frame(
    dia     = tvec,
    Rt      = colMeans(res$Rt_nos),
    Cenario = res$nome))
}
df_Rt$Cenario <- factor(
  df_Rt$Cenario, levels=names(cores_cen))

p7 <- ggplot(df_Rt,
             aes(x=dia, y=Rt,
                 color=Cenario,
                 linetype=Cenario)) +
  geom_line(linewidth=1.0) +
  geom_hline(yintercept=1,
             linetype="dotted",
             color="grey30",
             linewidth=0.8) +
  annotate("text", x=75, y=1.12,
           label=expression(R[t]*"=1"),
           size=3.2, color="grey30") +
  scale_color_manual(values=cores_cen) +
  scale_linetype_manual(
    values=linhas_cen) +
  tema_base +
  labs(
    title    = paste0(
      "Número de Reprodução ",
      "Efectivo R_t"),
    subtitle = paste0(
      "Cenários A e B cruzam ",
      "R_t=1 no dia 71; ",
      "Cenário D não cruza até dia 89."),
    x = "Dias",
    y = expression(R[t])
  )

ggsave("Figuras/seir_rt.pdf",
       p7, width=9, height=5,
       device=cairo_pdf)

# Cenário C vs Base apenas
df_Rt_C <- df_Rt[df_Rt$Cenario %in%
  c("0-Base","C-Ctrl.Optimo"),]

p7b <- ggplot(df_Rt_C,
              aes(x=dia, y=Rt,
                  color=Cenario,
                  linetype=Cenario)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=1,
             linetype="dotted",
             color="grey30",
             linewidth=0.8) +
  annotate("text", x=75, y=1.08,
           label=expression(R[t]*"=1"),
           size=3.2, color="grey30") +
  scale_color_manual(
    values=cores_cen[
      c("0-Base","C-Ctrl.Optimo")]) +
  scale_linetype_manual(
    values=linhas_cen[
      c("0-Base","C-Ctrl.Optimo")]) +
  tema_base +
  labs(
    title    = paste0(
      "R_t — Cenário C vs Base"),
    subtitle = paste0(
      "Controlo óptimo abrevia ",
      "a fase de crescimento epidémico."),
    x = "Dias",
    y = expression(R[t])
  )

ggsave("Figuras/cenario_C_rt.pdf",
       p7b, width=9, height=5,
       device=cairo_pdf)
cat("  Figuras 7 e cenario_C_rt guardadas.\n")

# ------------------------------------------------
# 10. FIG 8 — Fronteira de eficiência
# ------------------------------------------------
cat("A gerar Figura 8 — Fronteira...\n")

tab_ace$Custo <- c(
  sum(0.5*100*resA$u_kg^2) +
    sum(0.5*100*resA$v_kg^2),
  sum(0.5*10 *resB$u_kg^2) +
    sum(0.5*10 *resB$v_kg^2),
  sum(0.5*100*resC$u_kg^2) +
    sum(0.5*100*resC$v_kg^2),
  sum(0.5*100*resD$u_kg^2))

tab_ace$Label <- c(
  "Prior.\nEtária",
  "Prior.\nRegional",
  "Ctrl.\nÓptimo",
  "Só\nDistanc.")

p8 <- ggplot(tab_ace,
             aes(x=Evitados,
                 y=Custo,
                 color=Cenario)) +
  geom_point(size=5) +
  geom_path(
    data=tab_ace[
      order(tab_ace$Evitados),],
    aes(group=1),
    color="grey60",
    linetype="dashed",
    linewidth=0.6) +
  geom_label_repel(
    aes(label=paste0(
      Label,"\n(",
      round(Red_pct,1),"%)")),
    size=3, fill="white",
    label.size=0.2,
    box.padding=0.6,
    point.padding=0.3,
    segment.color="grey60") +
  scale_color_manual(
    values=cores_cen[
      c("A-Prior.Etaria",
        "B-Prior.Regional",
        "C-Ctrl.Optimo",
        "D-So.Distanciamento")]) +
  scale_x_continuous(
    labels=scales::comma) +
  tema_base +
  theme(legend.position="none") +
  labs(
    title    = paste0(
      "Fronteira de Eficiência — ",
      "Plano Custo-Eficácia"),
    subtitle = paste0(
      "Percentagem = redução de ",
      "mortalidade vs cenário base."),
    x = "Óbitos evitados",
    y = "Custo total da intervenção"
  )

ggsave("Figuras/fronteira_eficiencia.pdf",
       p8, width=8, height=6,
       device=cairo_pdf)
cat("  Figura 8 guardada.\n")

# ------------------------------------------------
# 11. FIG 9 — Óbitos por faixa etária
# ------------------------------------------------
cat("A gerar Figura 9 — Faixas etárias...\n")

df_etaria_long <- tab_etaria %>%
  select(Faixa, Base, CenA,
         CenB, CenC, CenD) %>%
  pivot_longer(
    cols = c(Base,CenA,CenB,
             CenC,CenD),
    names_to  = "Cenario",
    values_to = "Obitos") %>%
  mutate(
    Cenario = recode(Cenario,
      "Base" = "0-Base",
      "CenA" = "A-Prior.Etaria",
      "CenB" = "B-Prior.Regional",
      "CenC" = "C-Ctrl.Optimo",
      "CenD" = "D-So.Distanciamento"),
    Faixa = factor(Faixa,
                   levels=grupos))

df_etaria_long$Cenario <- factor(
  df_etaria_long$Cenario,
  levels=names(cores_cen))

p9a <- ggplot(df_etaria_long,
              aes(x=Faixa,
                  y=Obitos,
                  fill=Cenario)) +
  geom_col(position="dodge",
           width=0.7) +
  scale_fill_manual(
    values=cores_cen) +
  tema_base +
  theme(
    axis.text.x=element_text(
      angle=45, hjust=1)) +
  labs(
    title = "Óbitos por Faixa Etária",
    x = "Faixa etária",
    y = "Óbitos acumulados"
  )

# Redução % do Cenário C por faixa
tab_etaria$Faixa <- factor(
  tab_etaria$Faixa, levels=grupos)

p9b <- ggplot(tab_etaria,
              aes(x=Faixa,
                  y=Red_pct,
                  fill=Faixa)) +
  geom_col(width=0.7) +
  geom_text(
    aes(label=paste0(
      round(Red_pct,1),"%")),
    vjust=-0.5, size=3,
    fontface="bold") +
  scale_fill_manual(
    values=cores_grupos) +
  scale_y_continuous(
    limits=c(0,100),
    labels=function(x)
      paste0(x,"%")) +
  tema_base +
  theme(
    legend.position="none",
    axis.text.x=element_text(
      angle=45, hjust=1)) +
  labs(
    title = paste0(
      "Redução de Óbitos — ",
      "Cenário C por Faixa Etária"),
    x = "Faixa etária",
    y = "Redução %"
  )

p9 <- (p9a/p9b) +
  plot_annotation(
    title    = paste0(
      "Impacto Etário do ",
      "Controlo Óptimo"),
    subtitle = paste0(
      "Painel superior: óbitos ",
      "absolutos por cenário. ",
      "Painel inferior: redução ",
      "percentual do Cenário C ",
      "relativamente ao base."),
    theme=theme(
      plot.title=element_text(
        face="bold",size=12))
  )

ggsave("Figuras/obitos_etarios.pdf",
       p9, width=10, height=10,
       device=cairo_pdf)
cat("  Figura 9 guardada.\n")

# ------------------------------------------------
# 12. FIG 10 — Infectados por faixa etária
# ------------------------------------------------
cat("A gerar Figura 10 — I(t) etário...\n")

df_Ig <- data.frame()
for (g in 1:G) {
  df_Ig <- rbind(df_Ig, data.frame(
    dia   = rep(tvec, 2),
    I     = c(res0$I_grupos[g,]/1000,
              resC$I_grupos[g,]/1000),
    Cen   = rep(
      c("0-Base","C-Ctrl.Optimo"),
      each=T_dias+1),
    Faixa = grupos[g]))
}
df_Ig$Faixa <- factor(df_Ig$Faixa,
                      levels=grupos)
df_Ig$Cen   <- factor(df_Ig$Cen,
  levels=c("0-Base","C-Ctrl.Optimo"))

p10 <- ggplot(df_Ig,
              aes(x=dia, y=I,
                  color=Cen,
                  linetype=Cen)) +
  geom_line(linewidth=0.8) +
  facet_wrap(~Faixa,
             scales="free_y",
             ncol=3) +
  scale_color_manual(
    values=c(
      "0-Base"        = "grey50",
      "C-Ctrl.Optimo" = "#27AE60")) +
  scale_linetype_manual(
    values=c(
      "0-Base"        = "dashed",
      "C-Ctrl.Optimo" = "solid")) +
  tema_base +
  labs(
    title    = paste0(
      "Infectados I(t) por ",
      "Faixa Etária — ",
      "Base vs Cenário C"),
    subtitle = paste0(
      "Cada painel é uma faixa etária. ",
      "Escala Y livre entre painéis."),
    x = "Dias",
    y = "Infectados (mil)"
  )

ggsave("Figuras/infectados_etarios.pdf",
       p10, width=12, height=10,
       device=cairo_pdf)
cat("  Figura 10 guardada.\n")

# ------------------------------------------------
# 13. FIG 11 — Sensibilidade pesos B e C
# ------------------------------------------------
cat("A gerar Figura 11 — Sensibilidade...\n")

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
        sum(p$Theta[l,m]*p$N_kg[l,])
    }
  }
  N_hat <- pmax(N_hat,1)
  lam   <- matrix(0, p$K, p$G)
  for (k in 1:p$K) {
    for (g in 1:p$G) {
      soma <- 0
      for (m in 1:p$K) {
        inf_m <- numeric(p$G)
        for (l in 1:p$K) {
          inf_m <- inf_m +
            p$Theta[l,m] *
            (I_kg[l,]+p$delta*Iv_kg[l,])
        }
        soma <- soma +
          p$Theta[k,m] *
          sum(p$C_mat[g,]*
              inf_m/N_hat[m])
      }
      lam[k,g] <- beta_k[k]*
                  (1-u_kg[k,g])*soma
    }
  }
  lam
}

rk4_meta <- function(x, v_kg, u_kg,
                     beta_k, p) {
  K <- p$K; G <- p$G
  f_dx <- function(x, v_kg, u_kg) {
    lam <- calc_lambda(
      x$I, x$Iv, u_kg, beta_k, p)
    list(
      Sx=-lam*x$Sx,
      Su=-(lam+v_kg)*x$Su,
      Sv=v_kg*x$Su-(lam+1/p$T_V)*x$Sv,
      Sp=(1-p$e)*x$Sv/p$T_V-
         (1-p$e)*(1-u_kg)*lam*x$Sp,
      E=lam*(x$Sx+x$Su+x$Sv)-
        p$sigma*x$E,
      Ev=(1-p$e)*(1-u_kg)*lam*x$Sp-
          p$sigma*x$Ev,
      I=p$sigma*x$E-p$gamma*x$I,
      Iv=p$sigma*x$Ev-p$gamma*x$Iv,
      Q=(1-p$p_h)*p$gamma*(x$I+x$Iv)-
        x$Q/p$T_Q,
      Hw=p$p_h*p$gamma*(x$I+x$Iv)-
         (1/p$T_Hw+p$mu_kg)*x$Hw,
      Hc=p$p_c*x$Hw/p$T_Hw-
         (1/p$T_Hc+p$mu_kg)*x$Hc,
      R=p$e*x$Sv/p$T_V+x$Q/p$T_Q,
      RH=(1-p$p_c)*x$Hw/p$T_Hw+
          x$Hc/p$T_Hc,
      D=p$mu_kg*x$Q,
      DH=p$mu_kg*(x$Hw+x$Hc),
      Vw=v_kg*x$Su)
  }
  add_x <- function(a,b,s=1) {
    nms <- names(a)
    setNames(lapply(nms, function(n)
      matrix(as.numeric(a[[n]])+
             s*as.numeric(b[[n]]),
             K,G)),nms)
  }
  k1 <- f_dx(x,v_kg,u_kg)
  k2 <- f_dx(add_x(x,k1,.5),v_kg,u_kg)
  k3 <- f_dx(add_x(x,k2,.5),v_kg,u_kg)
  k4 <- f_dx(add_x(x,k3,1.),v_kg,u_kg)
  nms <- names(x)
  setNames(lapply(nms, function(n) {
    val <- as.numeric(x[[n]])+
           (as.numeric(k1[[n]])+
            2*as.numeric(k2[[n]])+
            2*as.numeric(k3[[n]])+
            as.numeric(k4[[n]]))/6
    matrix(pmax(0,val),K,G)
  }),nms)
}

sim_rapida <- function(x0, v_kg, u_kg,
                       beta_k, p,
                       T_dias, nome) {
  N_t  <- T_dias+1
  hist <- vector("list",N_t)
  hist[[1]] <- x0
  for (t in 1:T_dias) {
    hist[[t+1]] <- rk4_meta(
      hist[[t]],v_kg,u_kg,beta_k,p)
  }
  I_t <- numeric(N_t)
  D_t <- numeric(N_t)
  Rt_nos <- matrix(0,K,N_t)
  for (ti in 1:N_t) {
    h       <- hist[[ti]]
    I_t[ti] <- sum(h$I)+sum(h$Iv)
    D_t[ti] <- sum(h$D)+sum(h$DH)
    for (k in 1:K) {
      S_k <- sum(h$Sx[k,])+
             sum(h$Su[k,])+
             sum(h$Sv[k,])
      Rt_nos[k,ti] <- beta_k[k]*S_k/
        (p$gamma*sum(p$N_kg[k,]))
    }
  }
  list(nome=nome,I_t=I_t,D_t=D_t,
       Rt_nos=Rt_nos,
       v_kg=v_kg,u_kg=u_kg)
}

fbs_sens <- function(x0,beta_k,p,
                     T_dias,A_kg,
                     B_peso,C_peso,
                     nome) {
  cat("  FBS:",nome,"\n")
  K <- p$K; G <- p$G
  v_kg  <- matrix(0,K,G)
  u_kg  <- matrix(0,K,G)
  v_max <- matrix(0.05,K,G)
  u_max <- 0.60
  omega <- 0.5; tol <- 1e-6
  err   <- Inf; it <- 0
  N_t   <- T_dias+1
  while (err>tol && it<500) {
    v_old <- v_kg; u_old <- u_kg
    hist <- vector("list",N_t)
    hist[[1]] <- x0
    for (t in 1:T_dias) {
      hist[[t+1]] <- rk4_meta(
        hist[[t]],v_kg,u_kg,beta_k,p)
    }
    v_num <- matrix(0,K,G)
    u_num <- matrix(0,K,G)
    count <- 0
    for (ti in 1:N_t) {
      tau     <- T_dias-(ti-1)
      q_I_kg  <- A_kg/p$gamma*
                 (1-exp(-p$gamma*tau))
      x_t <- hist[[ti]]
      lam <- calc_lambda(
        x_t$I,x_t$Iv,u_kg,beta_k,p)
      q_Su <- lam*q_I_kg
      q_Sv <- q_Su*p$e
      vh   <- (q_Su-q_Sv)*
               x_t$Su/B_peso
      vh   <- matrix(pmax(pmin(
        as.numeric(vh),
        as.numeric(v_max)),0),K,G)
      S_tot <- x_t$Sx+x_t$Su+x_t$Sv
      uh    <- lam*q_I_kg*S_tot/C_peso
      uh    <- matrix(pmax(pmin(
        as.numeric(uh),u_max),0),K,G)
      v_num <- v_num+vh
      u_num <- u_num+uh
      count <- count+1
    }
    v_hat <- v_num/count
    u_hat <- u_num/count
    v_kg  <- omega*v_hat+(1-omega)*v_old
    u_kg  <- omega*u_hat+(1-omega)*u_old
    err   <- max(
      max(abs(v_kg-v_old))/
        (max(abs(v_kg))+1e-10),
      max(abs(u_kg-u_old))/
        (max(abs(u_kg))+1e-10))
    it <- it+1
  }
  cat(sprintf(
    "    Iter:%d Erro:%.2e\n",it,err))
  res <- sim_rapida(x0,v_kg,u_kg,
                    beta_k,p,T_dias,nome)
  res$v_kg <- v_kg
  res$u_kg <- u_kg
  res
}

A_kg_base <- r$A_kg_base
B_vals  <- c(100,  500,  1000)
C_vals  <- c(200, 1000,  2000)
nomes_s <- c("B=100 (Sanitário)",
             "B=500 (Moderado)",
             "B=1000 (Económico)")

df_s_I  <- data.frame()
df_s_Rt <- data.frame()
tab_s   <- data.frame()
x0_s    <- r$x0

for (i in seq_along(B_vals)) {
  res_s <- fbs_sens(
    x0_s, beta_k, params,
    T_dias, A_kg_base,
    B_vals[i], C_vals[i],
    nomes_s[i])

  df_s_I <- rbind(df_s_I, data.frame(
    dia     = tvec,
    I       = res_s$I_t/1000,
    Cenario = nomes_s[i]))

  df_s_Rt <- rbind(df_s_Rt, data.frame(
    dia     = tvec,
    Rt      = colMeans(res_s$Rt_nos),
    Cenario = nomes_s[i]))

  D_s <- res_s$D_t[T_dias+1]
  tab_s <- rbind(tab_s, data.frame(
    Cenario = nomes_s[i],
    B       = B_vals[i],
    Obitos  = D_s,
    Red_pct = (D0-D_s)/D0*100,
    u_medio = mean(res_s$u_kg)))
}

cat("\n=== SENSIBILIDADE AOS PESOS ===\n")
cat(sprintf("%-25s %6s %8s %7s %7s\n",
    "Cenário","B","Óbitos",
    "Red%","u_med"))
cat(strrep("-",57),"\n")
for (i in 1:nrow(tab_s)) {
  cat(sprintf(
    "%-25s %6d %8.0f %6.1f%% %7.4f\n",
    tab_s$Cenario[i],tab_s$B[i],
    tab_s$Obitos[i],tab_s$Red_pct[i],
    tab_s$u_medio[i]))
}

cores_s <- c(
  "B=100 (Sanitário)" = "#C0392B",
  "B=500 (Moderado)"  = "#2980B9",
  "B=1000 (Económico)"= "#27AE60"
)
df_s_I$Cenario  <- factor(
  df_s_I$Cenario,  levels=nomes_s)
df_s_Rt$Cenario <- factor(
  df_s_Rt$Cenario, levels=nomes_s)

p11a <- ggplot(df_s_I,
               aes(x=dia,y=I,
                   color=Cenario,
                   linetype=Cenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=cores_s) +
  tema_base +
  labs(title="Infectados I(t)",
       x="Dias",y="Infectados (mil)")

p11b <- ggplot(df_s_Rt,
               aes(x=dia,y=Rt,
                   color=Cenario,
                   linetype=Cenario)) +
  geom_line(linewidth=1.0) +
  geom_hline(yintercept=1,
             linetype="dotted",
             color="grey40",
             linewidth=0.6) +
  scale_color_manual(values=cores_s) +
  tema_base +
  labs(title=expression(
         R[t]*" por cenário de B"),
       x="Dias",y=expression(R[t]))

amplitude <- max(tab_s$Red_pct)-
             min(tab_s$Red_pct)

p11 <- (p11a|p11b) +
  plot_annotation(
    title    = paste0(
      "Análise de Sensibilidade — ",
      "Parâmetro B"),
    subtitle = paste0(
      "A eficácia varia ",
      round(amplitude,1),
      " p.p. entre extremos; ",
      "o custo varia por factor ",
      max(B_vals)/min(B_vals),"x."),
    theme=theme(
      plot.title=element_text(
        face="bold",size=12))
  )

ggsave("Figuras/sensibilidade_pesos.pdf",
       p11, width=12, height=5,
       device=cairo_pdf)
cat("  Figura 11 guardada.\n")

# ------------------------------------------------
# 14. FIGURAS INDIVIDUAIS POR CENÁRIO
# ------------------------------------------------
cat("A gerar figuras por cenário...\n")

for (res in list(resA,resB,resD)) {
  nome_fig <- switch(res$nome,
    "A-Prior.Etaria"      = "cenario_A",
    "B-Prior.Regional"    = "cenario_B",
    "D-So.Distanciamento" = "cenario_D")

  titulo <- switch(res$nome,
    "A-Prior.Etaria"      =
      "Cenário A — Prioridade Etária",
    "B-Prior.Regional"    =
      "Cenário B — Prioridade Regional",
    "D-So.Distanciamento" =
      "Cenário D — Só Distanciamento")

  red_pct <- tab_ace$Red_pct[
    tab_ace$Cenario == res$nome]

  df_tmp <- data.frame()
  for (k in 1:K) {
    df_tmp <- rbind(df_tmp, data.frame(
      dia = tvec,
      I   = res$I_nos[k,]/1000,
      No  = nos_nomes[k]))
  }
  df_tmp$No <- factor(df_tmp$No,
                      levels=nos_nomes)

  p_tmp <- ggplot(df_tmp,
                  aes(x=dia,y=I,
                      color=No)) +
    geom_line(
      data=data.frame(
        dia=tvec,
        I=res0$I_t/1000,
        No="0-Base"),
      aes(x=dia,y=I),
      color="grey50",
      linetype="dotted",
      linewidth=0.8,
      inherit.aes=FALSE) +
    geom_line(linewidth=0.9) +
    scale_color_brewer(palette="Dark2") +
    tema_base +
    labs(
      title    = titulo,
      subtitle = paste0(
        "Infectados por nó. ",
        "Linha cinzenta: Cenário Base. ",
        "Redução global: ",
        round(red_pct,1),"%"),
      x = "Dias",
      y = "Infectados activos (mil)"
    )

  ggsave(paste0("Figuras/",
                nome_fig,".pdf"),
         p_tmp, width=9, height=5,
         device=cairo_pdf)
  cat(" ", nome_fig, ".pdf guardado.\n")
}

# Cenário C adjuntas
q_I_traj <- 1/params$gamma *
            (1-exp(-params$gamma*
                   rev(tvec)))

df_adj <- data.frame(
  dia      = rep(tvec, 2),
  Valor    = c(q_I_traj,
               colMeans(resC$Rt_nos)),
  Variavel = rep(
    c("q_I(t) — Custo marginal",
      "R_t médio — Cenário C"),
    each=T_dias+1)
)

pC_adj <- ggplot(df_adj,
                 aes(x=dia,y=Valor,
                     color=Variavel,
                     linetype=Variavel)) +
  geom_line(linewidth=1.0) +
  geom_hline(
    data=data.frame(
      Variavel="R_t médio — Cenário C",
      y=1),
    aes(yintercept=y),
    linetype="dotted",
    color="grey40",
    linewidth=0.6,
    inherit.aes=FALSE) +
  scale_color_manual(
    values=c(
      "q_I(t) — Custo marginal"="#C0392B",
      "R_t médio — Cenário C"="#27AE60")) +
  scale_linetype_manual(
    values=c(
      "q_I(t) — Custo marginal"="solid",
      "R_t médio — Cenário C"="dashed")) +
  facet_wrap(~Variavel,
             scales="free_y") +
  tema_base +
  theme(legend.position="none") +
  labs(
    title    = paste0(
      "Cenário C — Variável ",
      "Adjunta e R_t"),
    subtitle = paste0(
      "q_I(T)=0: transversalidade ",
      "verificada. R_t cruza limiar ",
      "no dia ",
      which(colMeans(resC$Rt_nos)<1)[1],"."),
    x = "Dias", y = "Valor"
  )

ggsave("Figuras/cenario_C_adjuntas.pdf",
       pC_adj, width=10, height=5,
       device=cairo_pdf)
cat("  cenario_C_adjuntas.pdf guardado.\n")

# ------------------------------------------------
# 15. CONFIRMAÇÃO FINAL
# ------------------------------------------------
cat("\n=== FIGURAS GERADAS ===\n")
figs <- list.files("Figuras/",
                   pattern="\\.pdf$")
for (f in sort(figs)) cat(" ",f,"\n")
cat("Total:",length(figs),
    "ficheiros PDF\n")