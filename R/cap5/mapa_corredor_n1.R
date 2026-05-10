# ================================================
# mapa_corredor_n1.R — Versão final
# + Rio Zambeze com dados reais
# + Linhas de fluxo Θ_{km}
# + Gradiente β_k nos nós
# + Inset África no Painel A
# ================================================
library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggspatial)
library(cowplot)
library(ggrepel)

dir.create("Figuras", showWarnings = FALSE)

# ------------------------------------------------
# 1. DADOS GEOGRÁFICOS
# ------------------------------------------------
cat("A carregar dados geográficos...\n")

africa <- ne_countries(
  scale       = "medium",
  continent   = "Africa",
  returnclass = "sf"
)

mundo <- ne_countries(
  scale       = "small",
  returnclass = "sf"
)

moz <- africa[
  africa$name == "Mozambique", ]

vizinhos <- africa[
  africa$name %in% c(
    "Tanzania", "Malawi", "Zambia",
    "Zimbabwe", "South Africa",
    "Eswatini", "Mozambique"), ]

cat("Dados base carregados.\n")

# ------------------------------------------------
# 2. RIO ZAMBEZE — dados reais
# ------------------------------------------------
cat("A carregar Rio Zambeze...\n")

rios <- ne_download(
  scale       = 10,
  type        = "rivers_lake_centerlines",
  category    = "physical",
  returnclass = "sf"
)

zambeze_sf <- rios[
  grepl("Zambezi|Zambeze|Zambèze",
        rios$name,
        ignore.case = TRUE), ]

# Recortar ao domínio do mapa
bbox_moz <- st_bbox(c(
  xmin=31.5, xmax=37.8,
  ymin=-26.8, ymax=-16.0),
  crs=st_crs(zambeze_sf))

zambeze_sf <- st_crop(
  zambeze_sf, bbox_moz)

cat("Zambeze carregado:",
    nrow(zambeze_sf), "linhas\n")

# ------------------------------------------------
# 3. NÓS DO CORREDOR N1
# ------------------------------------------------
nos <- data.frame(
  k    = 1:6,
  nome = c("Maputo","Xai-Xai",
           "Maxixe","Caia",
           "Inchope","Namacurra"),
  lon  = c(32.573, 33.644, 35.347,
           35.338, 33.961, 36.883),
  lat  = c(-25.966,-25.052,-23.860,
           -17.857,-19.133,-17.083),
  tipo = c("Hub principal",
           "Nó regular","Nó regular",
           "Nó crítico",
           "Nó regular","Nó regular"),
  beta = c(0.1484, 0.1637, 0.2487,
           0.1400, 0.1943, 0.2249),
  R0   = c(1.855, 2.046, 3.108,
           1.751, 2.429, 2.812)
)

# Capitais provinciais
capitais <- data.frame(
  nome = c("Pemba","Lichinga",
           "Nampula","Quelimane",
           "Tete","Chimoio",
           "Beira","Inhambane"),
  lon  = c(40.517, 35.233,
           39.261, 36.888,
           33.589, 33.463,
           34.838, 35.383),
  lat  = c(-12.972,-13.313,
           -15.116,-17.878,
           -16.156,-19.117,
           -19.844,-23.865)
)

# ------------------------------------------------
# 4. LINHAS DE FLUXO Θ_{km}
# ------------------------------------------------
fluxos <- data.frame(
  k_orig = c(1,2, 2,3, 3,4,
             4,5, 5,6, 1,4),
  k_dest = c(2,1, 3,2, 4,3,
             5,4, 6,5, 4,1),
  theta  = c(0.072, 0.065,
             0.048, 0.042,
             0.058, 0.051,
             0.042, 0.038,
             0.061, 0.054,
             0.030, 0.025)
)

linhas_list <- vector("list",
                      nrow(fluxos))
for (i in seq_len(nrow(fluxos))) {
  k1 <- fluxos$k_orig[i]
  k2 <- fluxos$k_dest[i]
  coords <- matrix(c(
    nos$lon[k1], nos$lat[k1],
    nos$lon[k2], nos$lat[k2]
  ), ncol=2, byrow=TRUE)
  linhas_list[[i]] <-
    st_linestring(coords)
}

fluxos_sf <- st_sf(
  theta    = fluxos$theta,
  geometry = st_sfc(
    linhas_list, crs=4326)
)

# ------------------------------------------------
# 5. TRAÇADO DA N1
# ------------------------------------------------
n1_coords <- matrix(c(
  32.573,-25.966,
  32.900,-25.500,
  33.100,-25.200,
  33.644,-25.052,
  34.200,-24.500,
  34.800,-24.000,
  35.100,-23.900,
  35.347,-23.860,
  35.400,-23.000,
  35.300,-22.000,
  35.200,-21.000,
  35.100,-20.000,
  34.800,-19.500,
  34.200,-19.200,
  33.961,-19.133,
  34.200,-18.500,
  34.800,-18.200,
  35.100,-18.000,
  35.338,-17.857,
  35.500,-17.500,
  36.000,-17.200,
  36.500,-17.100,
  36.883,-17.083
), ncol=2, byrow=TRUE)

n1_sf <- st_sf(
  geometry = st_sfc(
    st_linestring(n1_coords),
    crs=4326))

# ------------------------------------------------
# 6. TEMA BASE
# ------------------------------------------------
tema_mapa <- theme_void(base_size=10) +
  theme(
    plot.background = element_rect(
      fill      = "white",
      color     = "grey80",
      linewidth = 0.5),
    plot.margin     = margin(6,6,6,6),
    legend.position = "bottom",
    legend.title    = element_text(
      size=8, face="bold"),
    legend.text     = element_text(
      size=7.5),
    legend.key.size = unit(0.45,"cm"),
    plot.title      = element_text(
      size=9, face="bold",
      hjust=0.5, color="#2C3E50"),
    plot.subtitle   = element_text(
      size=7.5, hjust=0.5,
      color="grey50")
  )

# ------------------------------------------------
# 7. INSET — ÁFRICA COMPLETA
# ------------------------------------------------
cat("A construir inset África...\n")

inset_africa <- ggplot() +

  geom_sf(data      = mundo,
          fill      = "grey88",
          color     = "white",
          linewidth = 0.2) +

  geom_sf(data      = africa,
          fill      = "grey75",
          color     = "white",
          linewidth = 0.2) +

  geom_sf(data      = moz,
          fill      = "#E74C3C",
          color     = "white",
          linewidth = 0.3) +

  coord_sf(xlim   = c(-20, 52),
           ylim   = c(-36, 38),
           expand = FALSE) +

  theme_void() +
  theme(
    plot.background = element_rect(
      fill      = "#D6EAF8",
      color     = "grey60",
      linewidth = 0.4),
    plot.margin = margin(2,2,2,2)
  )

# ------------------------------------------------
# 8. PAINEL A — CONTEXTO REGIONAL
# ------------------------------------------------
cat("A construir Painel A...\n")

pA_base <- ggplot() +

  geom_sf(data      = vizinhos,
          fill      = "grey92",
          color     = "white",
          linewidth = 0.3) +

  geom_sf(data      = moz,
          fill      = "#D6EAF8",
          color     = "#2C3E50",
          linewidth = 0.6) +

  annotate("rect",
           xmin=31.5, xmax=37.8,
           ymin=-26.8, ymax=-16.0,
           fill=NA, color="#E74C3C",
           linewidth=0.9,
           linetype="dashed") +

  annotate("text", x=38.5, y=-12.5,
           label="Tanzânia",
           size=2.8, color="grey40",
           fontface="italic") +
  annotate("text", x=33.0, y=-13.8,
           label="Malawi",
           size=2.8, color="grey40",
           fontface="italic") +
  annotate("text", x=30.2, y=-15.5,
           label="Zâmbia",
           size=2.8, color="grey40",
           fontface="italic") +
  annotate("text", x=30.0, y=-20.5,
           label="Zimbabwe",
           size=2.8, color="grey40",
           fontface="italic") +
  annotate("text", x=28.2, y=-25.5,
           label="África\ndo Sul",
           size=2.8, color="grey40",
           fontface="italic") +
  annotate("text", x=35.0, y=-21.0,
           label="MOÇAMBIQUE",
           size=3.2, fontface="bold",
           color="#2C3E50") +
  annotate("text", x=36.8, y=-15.3,
           label="Corredor N1",
           size=2.6, color="#E74C3C",
           fontface="italic") +
  annotate("text", x=40.8, y=-19.5,
           label="Oceano\nÍndico",
           size=2.8, color="#5DADE2",
           fontface="italic") +

  coord_sf(xlim   = c(27, 42),
           ylim   = c(-27, -10),
           expand = FALSE) +

  annotation_north_arrow(
    location    = "tl",
    which_north = "true",
    height = unit(0.8,"cm"),
    width  = unit(0.6,"cm"),
    style  = north_arrow_fancy_orienteering(
      text_size=6)) +

  annotation_scale(
    location   = "bl",
    width_hint = 0.3,
    text_cex   = 0.6) +

  labs(
    title    = "(a) Localização Regional",
    subtitle = "Moçambique na África Austral"
  ) +

  tema_mapa

# Combinar Painel A com inset
pA <- ggdraw(pA_base) +
  draw_plot(
    inset_africa,
    x=0.62, y=0.62,
    width=0.36, height=0.36
  )

# ------------------------------------------------
# 9. PAINEL B — CORREDOR N1
# ------------------------------------------------
cat("A construir Painel B...\n")

pB <- ggplot() +

  # Fundo Moçambique
  geom_sf(data      = moz,
          fill      = "#EBF5FB",
          color     = "#2C3E50",
          linewidth = 0.6) +

  # Linhas de fluxo Θ_{km}
  geom_sf(data        = fluxos_sf,
          aes(linewidth = theta),
          color       = "#8E44AD",
          alpha       = 0.55,
          show.legend = TRUE) +

  scale_linewidth_continuous(
    name   = expression(Theta[km]),
    range  = c(0.4, 2.8),
    breaks = c(0.025, 0.050, 0.072),
    labels = c("0.025","0.050","0.072")
  ) +

  # Rio Zambeze — dados reais
  geom_sf(data      = zambeze_sf,
          color     = "#3498DB",
          linewidth = 1.4,
          inherit.aes = FALSE) +

  annotate("text",
           x=33.5, y=-17.25,
           label="Rio Zambeze",
           size=2.8, color="#2980B9",
           fontface="italic",
           angle=-8) +

  # Estrada N1
  geom_sf(data      = n1_sf,
          color     = "#E74C3C",
          linewidth = 0.8,
          linetype  = "dashed") +

  annotate("text",
           x=32.05, y=-23.2,
           label="N1",
           size=3.0, color="#E74C3C",
           fontface="bold",
           angle=72) +

  # Capitais provinciais
  geom_point(
    data   = capitais,
    aes(x=lon, y=lat),
    shape  = 4,
    size   = 2.0,
    color  = "grey40",
    stroke = 0.8) +

  geom_text_repel(
    data          = capitais,
    aes(x=lon, y=lat, label=nome),
    size          = 2.2,
    color         = "grey35",
    fontface      = "italic",
    box.padding   = 0.25,
    point.padding = 0.2,
    segment.color = "grey60",
    segment.size  = 0.3,
    max.overlaps  = 20) +

  # Nós com gradiente β_k
  geom_point(
    data   = nos,
    aes(x=lon, y=lat,
        fill=beta,
        shape=tipo),
    size   = 5.5,
    color  = "#2C3E50",
    stroke = 0.9) +

  scale_fill_gradient(
    name   = expression(
      beta[k]*" (dia"^{-1}*")"),
    low    = "#AED6F1",
    high   = "#1A5276",
    breaks = c(0.14, 0.18,
               0.22, 0.25),
    labels = c("0.14","0.18",
               "0.22","0.25")
  ) +

  scale_shape_manual(
    name   = "Tipo de nó",
    values = c(
      "Hub principal" = 23,
      "Nó crítico"    = 24,
      "Nó regular"    = 21)
  ) +

  # Labels dos nós com R0
  geom_label_repel(
    data          = nos,
    aes(x=lon, y=lat,
        label=paste0(
          "k=",k," ",nome,
          "\nR\u2080=",R0)),
    size          = 2.4,
    fontface      = "bold",
    fill          = "white",
    color         = "#2C3E50",
    box.padding   = 0.55,
    point.padding = 0.4,
    label.size    = 0.25,
    segment.color = "grey50",
    segment.size  = 0.4,
    max.overlaps  = 20) +

  coord_sf(xlim   = c(31.5, 37.8),
           ylim   = c(-26.8, -16.0),
           expand = FALSE) +

  annotation_scale(
    location   = "bl",
    width_hint = 0.25,
    text_cex   = 0.6) +

  annotation_north_arrow(
    location    = "tl",
    which_north = "true",
    height = unit(0.8,"cm"),
    width  = unit(0.6,"cm"),
    style  = north_arrow_fancy_orienteering(
      text_size=6)) +

  labs(
    title    = paste0(
      "(b) Corredor N1 — ",
      "Grafo de Mobilidade"),
    subtitle = paste0(
      "Gradiente azul: \u03b2\u2096 ",
      "calibrado | ",
      "Espessura: fluxo \u0398\u2096\u2098 | ",
      "\u2715 Capitais provinciais")
  ) +

  tema_mapa +
  theme(
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.spacing.x = unit(0.3,"cm")
  ) +

  guides(
    fill      = guide_colorbar(
      barwidth       = 4,
      barheight      = 0.5,
      title.position = "top"),
    linewidth = guide_legend(
      title.position = "top"),
    shape     = guide_legend(
      title.position = "top",
      override.aes   = list(
        fill="grey70", size=3))
  )

# ------------------------------------------------
# 10. PAINEL COMBINADO
# ------------------------------------------------
cat("A combinar painéis...\n")

mapa_final <- plot_grid(
  pA, pB,
  ncol       = 2,
  rel_widths = c(0.42, 0.58)
)

# ------------------------------------------------
# 11. GUARDAR
# ------------------------------------------------
ggsave(
  "Figuras/mapa_corredor_n1.pdf",
  plot   = mapa_final,
  width  = 18,
  height = 12,
  units  = "cm",
  device = cairo_pdf,
  dpi    = 300
)

cat("\nMapa guardado em",
    "Figuras/mapa_corredor_n1.pdf\n")