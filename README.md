# benchmark_cams_forecast

Previsões globais do CAMS para MP2,5 na estação CETESB 73 (Congonhas), em médias diárias
por data de emissão e lead. São a entrada da trilha opcional de previsão do
[pm25bench](https://github.com/GabrielDSAlves/Benchmark) (`data.cams_forecast_dir`).

## Fonte

| | |
|---|---|
| Dataset | CAMS global atmospheric composition forecasts, Atmosphere Data Store (ADS), `cams-global-atmospheric-composition-forecasts`, DOI [10.24381/04a0b097](https://doi.org/10.24381/04a0b097) |
| Variável | `particulate_matter_2.5um` (no NetCDF: `pm2p5`, kg/m³) |
| Rodada | 00 UTC (`time: 00:00`, `type: forecast`) |
| Horas de previsão | 4 a 99, de hora em hora, as únicas usadas nos leads 0 a 3 (o ADS oferece passo horário para o MP2,5 desde 2015-01-01). Os blocos de 2015 foram baixados com as horas 0 a 102; as horas extras são ignoradas no preparo, e as médias diárias são idênticas |
| Área | caixa `[N, W, S, E] = [-22.8, -47.5, -24.5, -45.8]`, com pelo menos 2 células de 0,4° de margem em volta da estação |
| Grade | 0,4° × 0,4° (grade regular entregue pelo ADS) |
| Formato do pedido | `netcdf_zip`, um pedido por bloco de dias (1–10, 11–20, 21–fim de cada mês) |

As chaves e os valores do pedido foram conferidos no formulário e nas restrições do dataset
no ADS em 2026-11-05.

## Pastas de `data/`

| pasta | conteúdo |
|---|---|
| `data/73/` | previsão interpolada na estação 73 (Congonhas) |
| `data/media_sp5/` | média, hora a hora, da previsão interpolada em 5 estações do município de São Paulo (98, 270, 73, 95, 83), antes da média diária |

```bash
uv run python scripts/preparar_cams.py --estacao 98 270 73 95 83 --nome media_sp5
```

## Colunas de `data/<pasta>/{ano}.csv`

- **Rodada usada:** a previsão emitida às **00 UTC do dia t** (21h de t−1 no horário
  local). Ela fica disponível na manhã do dia t, coerente com uma emissão no dia t.
- **Dia CETESB:** o dia D vai de 01:00 de D até 24:00 de D, em horário local UTC−3 (a hora
  24 é registrada como 00:00 de D+1). Em UTC, isso é de **04 UTC de D até 03 UTC de D+1**,
  inclusive.
- **Prazos necessários a partir da rodada 00 UTC de t:** o dia t+k corresponde às horas de
  previsão **4 + 24k até 27 + 24k**, inclusive. Para k = 0…3, isso dá as horas 4 a 99.
- **Coluna de previsão:** `fc_pm2_5_l{k}` na linha da data de emissão t = média diária (dia
  CETESB) da previsão de MP2,5 emitida em t para o dia t+k, em µg/m³, na estação. Variáveis
  adicionais seguem o mesmo padrão (`fc_pm10_l{k}`, `fc_no2_l{k}`…). Todo nome de coluna de
  previsão começa com `fc_`.
- **Unidade:** o CAMS entrega kg/m³; os valores são multiplicados por 10⁹ para µg/m³.
- `date`: data de emissão t (`AAAA-MM-DD`). Uma linha por data, em ordem.

O valor na estação vem de **interpolação bilinear** do campo para as coordenadas da estação
73, lidas de `pm25bench/data/stations.csv` do benchmark. A média diária exige pelo menos 18
das 24 horas válidas (6 dos 8 passos, se os dados forem 3-horários), a mesma regra das
observações; abaixo disso, `NaN`.

## Cobertura

Datas de emissão de **2015-01-01 a 2025-12-31** (4.018 dias), todas com as quatro colunas
`fc_pm2_5_l0` a `fc_pm2_5_l3` válidas (pelo menos 18 das 24 horas em cada dia alvo).

| ano | dias | completos | faltantes |
|---|---|---|---|
| 2015 | 365 | 365 | 0 |
| 2016 | 366 | 366 | 0 |
| 2017 | 365 | 365 | 0 |
| 2018 | 365 | 365 | 0 |
| 2019 | 365 | 365 | 0 |
| 2020 | 366 | 366 | 0 |
| 2021 | 365 | 365 | 0 |
| 2022 | 365 | 365 | 0 |
| 2023 | 365 | 365 | 0 |
| 2024 | 366 | 366 | 0 |
| 2025 | 365 | 365 | 0 |

Download feito entre 2026-11-10 e 2026-11-19 (2025 baixado por último). Janeiro de 2016 falhou nas primeiras tentativas
(pedidos presos por horas no ADS) e veio numa nova execução; nenhum bloco ficou faltando.

## Aviso: série não estacionária

O sistema de previsão do CAMS é atualizado aproximadamente uma vez por ano (novos ciclos do
IFS, mudanças em emissões, assimilação e química). Em julho de 2019, por exemplo, o modelo
passou de 60 para 137 níveis verticais. O viés da previsão pode mudar a cada atualização, e
a série 2015–2024 não é estacionária. Um modelo de correção ajustado num período pode não
valer em outro. Isso afeta em particular a divisão treino 2015–2019 / teste 2022–2024 do
benchmark.

Nesta estação o efeito é grande. Média anual de `fc_pm2_5_l0` comparada com a média
observada em Congonhas (CETESB, dia com pelo menos 18 h válidas):

| ano | previsão CAMS (µg/m³) | observado (µg/m³) |
|---|---|---|
| 2015 | 9,5 | 19,8 |
| 2016 | 9,1 | 19,2 |
| 2017 | 40,8 | 18,0 |
| 2018 | 48,8 | 18,6 |
| 2019 | 41,4 | 17,8 |
| 2020 | 20,7 | 16,6 |
| 2021 | 31,3 | 16,6 |
| 2022 | 31,0 | 15,8 |
| 2023 | 27,1 | 16,7 |
| 2024 | 22,9 | 19,2 |
| 2025 | 21,1 | 15,2 |

O salto de 2017 acontece em **2017-01-25**: até 24/01 a previsão de lead 0 fica entre 3 e
10 µg/m³ e a partir de 25/01 entre 15 e 27 µg/m³. A data coincide com uma quebra de versão do
sistema listada nas restrições do dataset no ADS (2017-01-24). A queda na virada de 2019
para 2020 é gradual e não coincide com nenhuma data de troca listada; a causa não foi
identificada. O período de treino do benchmark (2015–2019) mistura, portanto, dois regimes
de viés muito diferentes, e o teste está num terceiro.

## Como reproduzir

```bash
uv sync
# credenciais do ADS em ~/.cdsapirc (nunca versionar); para outra conta:
#   export CDSAPI_RC=~/.cdsapirc_outra_conta
uv run python scripts/baixar_cams.py --inicio 2015-01-01 --fim 2024-12-31
uv run python scripts/preparar_cams.py --raw raw --saida data --estacao 73
uv run pytest
```

Os pedidos são feitos um bloco de ~10 dias por vez, em sequência. Um pedido de um mês
inteiro (30 datas × 103 horas) foi derrubado pelo ADS depois de 6 h rodando; blocos menores
pedem as mesmas horas e geram exatamente as mesmas médias diárias. Se o ADS rejeitar ou
apagar um pedido, o bloco é pedido de novo depois de 15 min, dobrando a espera a cada falha
(até 4 h, 6 tentativas). Um pedido que passa de 3 h **rodando** sem terminar (`--limite-h`; a
espera na fila não conta, porque passa de 2 h em dias movimentados) é cancelado no ADS e
refeito; a fila tem um teto separado de 24 h (`--limite-fila-h`). Ao reiniciar, o script
retoma um pedido idêntico que ainda esteja na fila, sem perder o lugar; se isso acontecer duas vezes, ou se o MARS responder que a fita está
indisponível, o bloco é anotado em `raw/indisponiveis.txt` e o download segue. Rodar o
script de novo tenta esses blocos outra vez. O ADS limita o número de pedidos na
fila por usuário e por dataset e rejeita os excedentes ("Maximum number queued requests for
this dataset has been reached. Please configure your scripts to request data
sequentially"). Dados com mais de 30 dias vêm do arquivo em fita e podem esperar bastante na
fila. O download é retomável: os blocos já baixados (pasta `raw/<bloco>/.completo`) são
pulados. Os NetCDF ficam em `raw/`, fora do git.

## Licença e atribuição

Os dados do CAMS são distribuídos pelo ADS sob a licença
[CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/), que exige crédito à fonte,
link para a licença e indicação das modificações. Atribuição:

> Contém informações modificadas do Copernicus Atmosphere Monitoring Service (CAMS),
> 2015–2025. Dados: CAMS global atmospheric composition forecasts, ECMWF, Atmosphere Data
> Store, doi:10.24381/04a0b097. Modificações: recorte espacial, interpolação bilinear para a
> estação CETESB 73, conversão para µg/m³ e médias diárias por lead. Nem a Comissão Europeia
> nem o ECMWF se responsabilizam pelo uso destas informações.
