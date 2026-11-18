"""Converte os NetCDF do CAMS em médias diárias (dia CETESB) por data de emissão e lead.

Para cada data de emissão t (rodada ``HH`` UTC) e cada lead k = 0 a 3, calcula a média das
horas de previsão que caem no dia CETESB t + k, interpolando bilinearmente o campo para a
estação e convertendo kg/m³ em µg/m³. Grava ``data/<estacao>/<ano>.csv``.

Uso::

    uv run python scripts/preparar_cams.py --raw raw --saida data --estacao 73
"""

from __future__ import annotations

# argparse: leitura dos parâmetros de linha de comando
import argparse

# math: arredondamento do mínimo de horas válidas
import math

# Path: caminhos de arquivo
from pathlib import Path

# numpy: vetores
import numpy as np

# pandas: tabelas de saída
import pandas as pd

# xarray: leitura dos NetCDF
import xarray as xr

# coordenadas da estação, lidas do benchmark
from baixar_cams import coordenadas_estacao

#: Nome da variável no NetCDF do ADS -> nome usado nas colunas ``fc_<nome>_l{k}``.
NOMES: dict[str, str] = {
    # MP2,5
    "pm2p5": "pm2_5",
    # MP10
    "pm10": "pm10",
    # dióxido de nitrogênio
    "no2": "no2",
    # ozônio
    "go3": "o3",
    # monóxido de carbono
    "co": "co",
    # dióxido de enxofre
    "so2": "so2",
}

#: kg/m³ -> µg/m³.
KG_M3_PARA_UG_M3 = 1e9

#: Leads diários gerados.
LEADS: list[int] = [0, 1, 2, 3]

#: Fração mínima de valores válidos no dia (18 de 24 horas, a regra das observações).
FRACAO_MINIMA = 18 / 24


def horas_previsao(lead: int, hora_rodada: int = 0) -> np.ndarray:
    """Horas de previsão (a partir da rodada) que caem no dia CETESB t + ``lead``.

    O dia CETESB D vai de 01:00 a 24:00 local (UTC-3), ou seja, de 04 UTC de D a 03 UTC de
    D+1, inclusive. Com a rodada às ``hora_rodada`` UTC do dia t, o dia t + k corresponde às
    horas ``4 + 24k - hora_rodada`` a ``27 + 24k - hora_rodada``. Para a rodada 00 UTC: 4 a 27
    no lead 0 e 76 a 99 no lead 3. Horas negativas (antes da rodada) são descartadas.

    Parameters
    ----------
    lead : int
        Dias depois da data de emissão.
    hora_rodada : int
        Hora UTC da rodada (0 ou 12).

    Returns
    -------
    numpy.ndarray
        Horas de previsão inteiras, em ordem crescente.
    """
    # primeira hora do dia alvo, contada da rodada
    inicio = 4 + 24 * lead - hora_rodada
    # 24 horas consecutivas, inclusive a última (27 + 24k - hora_rodada)
    horas = np.arange(inicio, inicio + 24)
    # só horas a partir da rodada existem na previsão
    return horas[horas >= 0]


def minimo_valores(passo_h: int) -> int:
    """Mínimo de valores válidos no dia: 18 de 24 se horário, 6 de 8 se 3-horário."""
    # valores esperados num dia com esse passo
    esperados = 24 // passo_h
    # 75% deles, arredondado para cima
    return math.ceil(FRACAO_MINIMA * esperados)


def _horas(fp: xr.DataArray) -> np.ndarray:
    """Converte a coordenada ``forecast_period`` em horas inteiras (float ou timedelta)."""
    # valores da coordenada
    v = fp.values
    # timedelta64: converte para horas
    if np.issubdtype(v.dtype, np.timedelta64):
        return (v / np.timedelta64(1, "h")).round().astype(int)
    # número: já está em horas
    return np.asarray(v, dtype="float64").round().astype(int)


def abrir(raw: Path, hora_rodada: int) -> xr.Dataset:
    """Abre e junta todos os NetCDF de ``raw`` (recursivo), só da rodada pedida."""
    # todos os NetCDF extraídos, em ordem
    arquivos = sorted(raw.rglob("*.nc"))
    # sem arquivos não há o que preparar
    if not arquivos:
        raise FileNotFoundError(f"nenhum NetCDF em {raw}")
    # abre cada arquivo e mantém só as variáveis conhecidas
    partes = []
    # percorre os arquivos
    for a in arquivos:
        # carrega o arquivo na memória (as caixas são pequenas)
        ds = xr.open_dataset(a).load()
        # variáveis de interesse presentes neste arquivo
        vs = [v for v in ds.data_vars if v in NOMES]
        # ignora arquivos sem variáveis conhecidas
        if vs:
            partes.append(ds[vs])
    # junta ao longo das datas de emissão; join="outer" porque blocos baixados com horas
    # diferentes (0-102 em 2015, 4-99 depois) ficam com NaN nas horas que não pediram,
    # e essas horas nunca entram nas médias diárias
    ds = xr.combine_by_coords(partes, combine_attrs="drop_conflicts", join="outer")
    # garante Dataset (combine_by_coords pode devolver DataArray)
    assert isinstance(ds, xr.Dataset)
    # horas UTC de cada rodada
    horas_ref = pd.DatetimeIndex(ds["forecast_reference_time"].values).hour
    # mantém só a rodada pedida
    return ds.isel(forecast_reference_time=np.flatnonzero(horas_ref == hora_rodada))


def interpolar(campo: xr.DataArray, lat: float, lon: float) -> xr.DataArray:
    """Interpolação bilinear do campo para o ponto (lat, lon).

    Em uma grade regular, a interpolação linear em latitude seguida da linear em longitude
    é a interpolação bilinear: o valor é a média dos quatro vizinhos ponderada pelas áreas
    opostas, ``f = Σ w_i·f_i`` com ``w_i = (1 - |Δx_i|/h)·(1 - |Δy_i|/h)``.
    """
    # dimensões que precisam ficar, mesmo com tamanho 1 (uma única rodada, por exemplo)
    manter = {"forecast_reference_time", "forecast_period", "latitude", "longitude"}
    # dimensões extras de tamanho 1 (por exemplo, nível de pressão) saem
    campo = campo.squeeze([d for d in campo.dims if d not in manter], drop=True)
    # latitudes e longitudes da grade
    lats, lons = campo["latitude"].values, campo["longitude"].values
    # o ponto precisa estar estritamente dentro da caixa para ter os quatro vizinhos
    if not (lats.min() < lat < lats.max() and lons.min() < lon < lons.max()):
        raise ValueError(f"ponto ({lat}, {lon}) fora da grade baixada; amplie --area")
    # interpolação linear nas duas direções = bilinear
    return campo.interp(latitude=lat, longitude=lon, method="linear")


def medias_diarias(serie: xr.DataArray, hora_rodada: int = 0) -> pd.DataFrame:
    """Média diária por data de emissão e lead, com a regra de completude.

    Parameters
    ----------
    serie : xarray.DataArray
        Valores na estação, em µg/m³, com dimensões ``forecast_reference_time`` e
        ``forecast_period``.
    hora_rodada : int
        Hora UTC da rodada.

    Returns
    -------
    pandas.DataFrame
        Índice = data de emissão, colunas ``l0`` ... ``l3``; ``NaN`` onde há menos valores
        válidos que o mínimo.
    """
    # horas de previsão disponíveis no arquivo
    horas = _horas(serie["forecast_period"])
    # passo da grade de horas (1 h ou 3 h)
    passo = int(np.min(np.diff(np.unique(horas)))) if len(horas) > 1 else 1
    # mínimo de valores válidos no dia para esse passo
    minimo = minimo_valores(passo)
    # matriz (emissão × horas de previsão)
    valores = serie.transpose("forecast_reference_time", "forecast_period").values
    # datas de emissão (dia da rodada)
    datas = pd.DatetimeIndex(serie["forecast_reference_time"].values).normalize()
    # resultado, uma coluna por lead
    saida = {}
    # percorre os leads
    for k in LEADS:
        # colunas cujas horas caem no dia t + k
        sel = np.isin(horas, horas_previsao(k, hora_rodada))
        # valores dessas horas
        janela = valores[:, sel]
        # quantos são válidos em cada data de emissão
        n_validos = np.sum(~np.isnan(janela), axis=1)
        # média ignorando NaN (sem aviso quando a linha está toda vazia)
        with np.errstate(invalid="ignore"):
            # soma dos válidos dividida pela contagem
            media = np.nansum(janela, axis=1) / np.where(n_validos > 0, n_validos, np.nan)
        # dias com poucos valores ficam NaN
        saida[f"l{k}"] = np.where(n_validos >= minimo, media, np.nan)
    # tabela final
    return pd.DataFrame(saida, index=pd.Index(datas, name="date"))


def preparar(raw: Path, lat: float, lon: float, hora_rodada: int = 0) -> pd.DataFrame:
    """Tabela ``date`` × ``fc_<var>_l{k}`` a partir dos NetCDF de ``raw``."""
    # NetCDF da rodada pedida
    ds = abrir(raw, hora_rodada)
    # uma tabela por variável
    tabelas = []
    # percorre as variáveis presentes
    for var in sorted(str(v) for v in ds.data_vars):
        # valor na estação, convertido para µg/m³
        serie = interpolar(ds[var], lat, lon) * KG_M3_PARA_UG_M3
        # médias diárias por lead
        t = medias_diarias(serie, hora_rodada)
        # renomeia l{k} -> fc_<nome>_l{k}
        t.columns = [f"fc_{NOMES[var]}_{c}" for c in t.columns]
        # guarda
        tabelas.append(t)
    # junta as variáveis lado a lado, uma linha por data de emissão
    out = pd.concat(tabelas, axis=1).sort_index()
    # remove datas duplicadas (meses sobrepostos), mantendo a primeira
    return out[~out.index.duplicated()]


def gravar(tabela: pd.DataFrame, pasta: Path) -> list[Path]:
    """Grava um CSV por ano (``<pasta>/<ano>.csv``) e devolve os caminhos."""
    # garante a pasta
    pasta.mkdir(parents=True, exist_ok=True)
    # caminhos gravados
    caminhos = []
    # agrupa pelo ano da data de emissão
    for ano, g in tabela.groupby(tabela.index.year):
        # arquivo do ano
        caminho = pasta / f"{ano}.csv"
        # data como AAAA-MM-DD e valores com 4 casas
        g.to_csv(caminho, date_format="%Y-%m-%d", float_format="%.4f")
        # guarda o caminho
        caminhos.append(caminho)
    # devolve a lista
    return caminhos


def resumo(tabela: pd.DataFrame) -> pd.DataFrame:
    """Por ano: dias no calendário, datas de emissão com o lead 0 a 3 completo e faltantes."""
    # calendário completo do período coberto
    dias = pd.date_range(tabela.index.min(), tabela.index.max(), freq="D")
    # datas com todos os leads de MP2,5 válidos
    cols = [c for c in tabela.columns if c.startswith("fc_pm2_5_")]
    # máscara de completude reindexada ao calendário
    ok = tabela[cols].notna().all(axis=1).reindex(dias, fill_value=False)
    # contagens por ano
    r = ok.groupby(ok.index.year).agg(dias="size", completos="sum")
    # faltantes = dias do calendário sem todos os leads
    r["faltantes"] = r["dias"] - r["completos"]
    # nome do índice
    r.index.name = "ano"
    # devolve o resumo
    return r


def main() -> None:
    """Lê os argumentos, prepara e grava as médias diárias."""
    # analisador de argumentos
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    # pasta com os NetCDF baixados
    p.add_argument("--raw", type=Path, default=Path("raw"))
    # pasta de saída (uma subpasta por estação)
    p.add_argument("--saida", type=Path, default=Path("data"))
    # estação CETESB
    p.add_argument("--estacao", type=int, default=73)
    # hora da rodada
    p.add_argument("--rodada", default="00:00", choices=["00:00", "12:00"])
    # lê a linha de comando
    a = p.parse_args()
    # hora da rodada como inteiro
    hora = int(a.rodada.split(":")[0])
    # coordenadas da estação, do benchmark
    lat, lon = coordenadas_estacao(a.estacao)
    # tabela diária
    tabela = preparar(a.raw, lat, lon, hora)
    # grava por ano
    caminhos = gravar(tabela, a.saida / str(a.estacao))
    # relatório
    print(f"estação {a.estacao} ({lat}, {lon}); {len(caminhos)} arquivos em {a.saida}")
    # resumo de cobertura para o README
    print(resumo(tabela).to_string())


# executa só quando chamado como script
if __name__ == "__main__":
    main()
