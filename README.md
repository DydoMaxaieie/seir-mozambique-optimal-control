# seir-mozambique-optimal-control

Código da dissertação de mestrado *Modelação Matemático-Epidemiológica da COVID-19 em Moçambique: modelo SEIR-V multirregional estruturado por idade e controlo óptimo* (Dário R. J. Maxaieie, UEM, Mestrado em Matemática e Aplicações, Biomatemática).

Modelo: 4 macrorregiões x 5 grupos etários x 13 compartimentos (260 EDO), mobilidade lagrangiana, matriz de contactos de Prem et al. (2021). Controlo óptimo por Princípio do Máximo de Pontryagin, resolvido por varrimento progressivo-regressivo com relaxação 0,5 e RK4 (h = 0,5 d).

## Instalação

    python -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt

## Ordem de execução

| Script | Conteúdo | Capítulo | Tempo |
|---|---|---|---|
| `scripts/01_verificacao_cap3.py` | FBS no problema escalar e no SEIR homogéneo | 3 (Tab. 3.1, 3.2) | segundos |
| `scripts/02_calibracao.py` | Calibração por mínimos quadrados e R0 por bloco | 5 (Tab. 5.5 a 5.7) | ~2 min |
| `scripts/03_controlo.py` | Quatro estratégias e sensibilidades | 6 (Tab. 6.2 a 6.5) | ~10 min (`--rapido`: ~1 min) |
| `scripts/04_verificacoes_cap6.py` | Gradiente vs diferenças finitas, 40 perturbações, convergência em h | 6 e Apêndice B | ~2 min |
| `scripts/05a_figuras_cap3.py` | Figuras do Capítulo 3 | 3 | segundos |
| `scripts/06_louvain.py` | Agrupamento de Louvain (200 sementes) | 4 | segundos |
| `scripts/07_calibracao_idade.py`, `08_controlo_idade.py`, `09_verificacao_idade.py`, `11_sensibilidade_idade.py` | Calibração, controlo, verificação e sensibilidade com o denominador N̂_mh | Resultados (normalização) | minutos |
| `scripts/10_horizonte.py` | Horizontes de 120 a 300 dias e óbitos aos 600 dias | Resultados (horizonte) | ~10 min |
| `scripts/12_validacao_previsao.py` | Validação fora da amostra a 14 dias | Resultados (validação) | ~2 min |
| `scripts/13_vacinacao_estratos.py` | Doses, cobertura e distanciamento por região e idade | Resultados (estratos) | segundos |
| `scripts/14_robustez_controlo.py ic` ou `clin` | Condições iniciais do FBS; eficácia vacinal e parâmetros clínicos | Resultados (robustez) | ~1 min |
| `scripts/15_matriz_reciproca.py cal` ou `ctl` | Matriz de contactos recíproca: recalibração e controlo | Resultados (robustez) | ~3 min |

Teste de regressão: `pytest -q`.

## Resultados de referência (caso base, rho = 0,1)

Sem controlo 17 673 óbitos; só vacinação 14 647; só distanciamento 4 614; vacinação e distanciamento 3 921. Convergência em h: 3920,6 / 3921,4 / 3921,6 (h = 1; 0,5; 0,25).

## Dados

`data/` contém as séries de casos por macrorregião, a matriz de fluxos, a população por idade e as matrizes de contacto agregadas. Os parâmetros calibrados estão em `results/calibracao/`; `scripts/02_calibracao.py` recalcula-os de raiz, e o optimizador local pode terminar em mínimos ligeiramente diferentes (R² dentro de 0,03 dos publicados), pelo que a dissertação usa os valores guardados.

## Licença

MIT.
