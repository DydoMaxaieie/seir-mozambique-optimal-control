# SEIR Mozambique Optimal Control

R implementation of a metapopulation
SEIR-V epidemiological model for
COVID-19 in Mozambique with optimal
control of vaccination and social
distancing.

## Description

This repository contains the R code
for the Master's dissertation:

> Maxaieie, D. R. J. (2026).
> *Modelação Matemático-Epidemiológica
> da COVID-19 em Moçambique: Um Modelo
> SEIR Multirregional e Etário com
> Controlo Óptimo*.
> Universidade Eduardo Mondlane,
> Maputo, Moçambique.

## Repository Structure

```
R/
├── cap3/
│   ├── figuras_escalar.R
│   ├── calibracao_seir.R
│   └── figuras_seir.R
└── cap5/
    ├── seir_meta.R
    ├── figuras_meta.R
    ├── mapa_corredor_n1.R
    └── sensibilidade_R0.R
dados/
├── covid_moz_Q1_2021.csv
└── n1_real.gpkg
```

## Requirements

R >= 4.0 and the following packages:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr",
  "scales", "patchwork", "sf",
  "rnaturalearth", "ggrepel",
  "ggspatial", "cowplot",
  "httr", "jsonlite"
))
```

## Reproduction

Run scripts in this order:

1. `R/cap5/seir_meta.R` — calibration
   and simulation (~15 min)
2. `R/cap5/figuras_meta.R` — figures
   (~8 min)
3. `R/cap5/mapa_corredor_n1.R` — map
   (~3 min)
4. `R/cap5/sensibilidade_R0.R` —
   sensitivity analysis (~2 min)

## License

MIT License — see `LICENSE` file.

## Author

Dário Rogério Júlio Maxaieie  
Universidade Eduardo Mondlane  
Maputo, Moçambique
