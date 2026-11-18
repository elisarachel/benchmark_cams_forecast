"""Testes de preparar_cams com NetCDF sintéticos no formato do ADS."""

from __future__ import annotations

# zipfile: ZIP falso do download
import zipfile

# Path: pasta temporária dos testes
from pathlib import Path

# numpy: valores sintéticos
import numpy as np

# pandas: datas
import pandas as pd

# pytest: verificação de exceções
import pytest

# xarray: criação dos NetCDF sintéticos
import xarray as xr

# funções testadas
from baixar_cams import (
    FITA_INDISPONIVEL,
    area_padrao,
    baixar_bloco,
    blocos,
    coordenadas_estacao,
    meses,
    pedido,
)
from preparar_cams import horas_previsao, medias_diarias, minimo_valores, preparar

# coordenadas reais da estação 73, lidas do benchmark
LAT, LON = coordenadas_estacao(73)


def _netcdf(
    pasta: Path, valor: float, horas: np.ndarray, emissoes: list[str], nome: str = "a.nc"
) -> np.ndarray:
    """Grava um NetCDF com a estrutura do ADS e devolve o array de dados (editável antes)."""
    # caixa padrão em volta da estação
    n, w, s, e = area_padrao(LAT, LON)
    # latitudes decrescentes, como no ADS
    lats = np.round(np.arange(n, s - 1e-9, -0.4), 2)
    # longitudes crescentes
    lons = np.round(np.arange(w, e + 1e-9, 0.4), 2)
    # campo constante (kg/m³) com as dimensões do ADS
    dados = np.full((len(horas), len(emissoes), len(lats), len(lons)), valor, dtype="float32")
    # Dataset no formato entregue pelo ADS
    ds = xr.Dataset(
        {"pm2p5": (("forecast_period", "forecast_reference_time", "latitude", "longitude"), dados)},
        coords={
            # horas de previsão como número
            "forecast_period": horas.astype("float64"),
            # rodadas
            "forecast_reference_time": pd.to_datetime(emissoes),
            # grade
            "latitude": lats,
            "longitude": lons,
        },
    )
    # grava
    ds.to_netcdf(pasta / nome)
    # devolve os dados
    return dados


def test_horas_previsao() -> None:
    """Rodada 00 UTC: lead 0 = horas 4 a 27; lead 3 = horas 76 a 99."""
    # lead 0
    assert list(horas_previsao(0)) == list(range(4, 28))
    # lead 3
    assert list(horas_previsao(3)) == list(range(76, 100))
    # rodada 12 UTC: o dia t começa antes da rodada, então sobram só as horas >= 0
    assert list(horas_previsao(0, 12)) == list(range(0, 16))


def test_minimo_valores() -> None:
    """18 de 24 horas se horário, 6 de 8 se 3-horário."""
    # passo horário
    assert minimo_valores(1) == 18
    # passo de 3 h
    assert minimo_valores(3) == 6


def test_valor_constante(tmp_path: Path) -> None:
    """Campo constante de 2e-8 kg/m³ vira 20 µg/m³ em todos os leads e datas."""
    # duas rodadas 00 UTC, horas 0 a 102
    _netcdf(tmp_path, 2e-8, np.arange(0, 103), ["2020-01-01", "2020-01-02"])
    # prepara na estação 73
    t = preparar(tmp_path, LAT, LON)
    # colunas e datas esperadas
    assert list(t.columns) == [f"fc_pm2_5_l{k}" for k in range(4)]
    # uma linha por data de emissão
    assert list(t.index) == list(pd.to_datetime(["2020-01-01", "2020-01-02"]))
    # todos os valores iguais à constante convertida (float32 -> tolerância relativa)
    np.testing.assert_allclose(t.to_numpy(), 20.0, rtol=1e-5)


def test_poucas_horas_vira_nan(tmp_path: Path) -> None:
    """Com 7 das 24 horas do lead 1 faltando (17 válidas < 18), o lead 1 vira NaN."""
    # horas 0 a 102
    horas = np.arange(0, 103)
    # campo constante
    dados = _netcdf(tmp_path, 1e-8, horas, ["2020-01-01"])
    # apaga 7 horas do dia t+1 (horas 28 a 51)
    dados[np.isin(horas, range(28, 35))] = np.nan
    # NetCDF a partir dos dados editados
    ds = xr.open_dataset(tmp_path / "a.nc").load()
    # substitui os valores
    ds["pm2p5"].values = dados
    # regrava por cima
    ds.to_netcdf(tmp_path / "b.nc")
    # remove o original
    (tmp_path / "a.nc").unlink()
    # prepara
    t = preparar(tmp_path, LAT, LON)
    # lead 1 sem horas suficientes
    assert np.isnan(t["fc_pm2_5_l1"].iloc[0])
    # os outros leads continuam válidos
    assert t[["fc_pm2_5_l0", "fc_pm2_5_l2", "fc_pm2_5_l3"]].notna().all(axis=None)


def test_passo_tres_horas() -> None:
    """Com passos de 3 h, a média usa os 8 passos do dia e exige 6."""
    # horas 0, 3, ..., 102
    horas = np.arange(0, 103, 3)
    # série na estação: valor = hora (para conferir quais horas entram)
    serie = xr.DataArray(
        horas[None, :].astype("float64"),
        dims=("forecast_reference_time", "forecast_period"),
        coords={
            "forecast_reference_time": pd.to_datetime(["2020-01-01"]),
            "forecast_period": horas.astype("float64"),
        },
    )
    # médias diárias
    t = medias_diarias(serie)
    # lead 0: passos 6, 9, ..., 27 (8 valores), média 16,5
    assert t["l0"].iloc[0] == pytest.approx(np.mean(np.arange(6, 28, 3)))


def test_meses_e_area() -> None:
    """Meses recortados pelo intervalo e caixa com margem em volta da estação."""
    # intervalo que começa e termina no meio do mês
    m = meses(pd.Timestamp("2020-01-15").date(), pd.Timestamp("2020-03-10").date())
    # três meses, recortados
    assert [(a.isoformat(), b.isoformat()) for a, b in m] == [
        ("2020-01-15", "2020-01-31"),
        ("2020-02-01", "2020-02-29"),
        ("2020-03-01", "2020-03-10"),
    ]
    # caixa [N, W, S, E]
    n, w, s, e = area_padrao(LAT, LON)
    # pelo menos uma célula de 0,4° de margem em cada lado
    assert n - LAT >= 0.4 and LAT - s >= 0.4 and LON - w >= 0.4 and e - LON >= 0.4


def test_blocos_fixos_no_calendario() -> None:
    """Cada mês vira 1-10, 11-20 e 21-fim; o recorte não muda os limites dos outros blocos."""
    # datas como date
    d = lambda x: pd.Timestamp(x).date()  # noqa: E731
    # fevereiro de 2016 (bissexto) inteiro
    b = blocos(d("2016-02-01"), d("2016-02-29"))
    # três blocos, o último até o dia 29
    assert [(x.isoformat(), y.isoformat()) for x, y in b] == [
        ("2016-02-01", "2016-02-10"),
        ("2016-02-11", "2016-02-20"),
        ("2016-02-21", "2016-02-29"),
    ]
    # começando no dia 15: o primeiro bloco é recortado, os outros ficam iguais
    b = blocos(d("2016-02-15"), d("2016-03-05"))
    # blocos esperados
    assert [(x.isoformat(), y.isoformat()) for x, y in b] == [
        ("2016-02-15", "2016-02-20"),
        ("2016-02-21", "2016-02-29"),
        ("2016-03-01", "2016-03-05"),
    ]


class RemotoFalso:
    """Simula o pedido remoto do ADS: devolve os status da lista, repetindo o último."""

    def __init__(self, status: list[str], erro: str = "", request: dict | None = None) -> None:
        # sequência de status das consultas
        self._status = list(status)
        # mensagem de erro de get_results
        self.erro = erro
        # conteúdo do pedido
        self.request = request or {}
        # se delete foi chamado
        self.apagado = False
        # identificador fictício
        self.request_id = "falso"

    @property
    def status(self) -> str:
        # consome a sequência, repetindo o último status
        return self._status.pop(0) if len(self._status) > 1 else self._status[0]

    def get_results(self) -> None:
        # como o ADS, levanta o erro do pedido
        raise RuntimeError(self.erro)

    def delete(self) -> None:
        # registra o cancelamento
        self.apagado = True

    def download(self, destino: str) -> None:
        # grava um ZIP com um NetCDF vazio de mentira
        with zipfile.ZipFile(destino, "w") as z:
            z.writestr("data_sfc.nc", b"")


class _Jobs:
    """Lista de pedidos abertos, como a devolvida por get_jobs."""

    def __init__(self, ids: list[str]) -> None:
        # identificadores
        self.request_ids = ids


class ClienteFalso:
    """Simula o cliente do ADS: todo pedido devolve um RemotoFalso com os status dados."""

    def __init__(self, status: list[str], erro: str = "") -> None:
        # status e erro dos pedidos
        self.status, self.erro = status, erro
        # pedidos feitos
        self.remotos: list[RemotoFalso] = []
        # pedidos já abertos no ADS (para testar a retomada)
        self.abertos: dict[str, RemotoFalso] = {}

    def submit(self, _: str, req: dict) -> RemotoFalso:
        # novo pedido
        r = RemotoFalso(self.status, self.erro, req)
        # registra
        self.remotos.append(r)
        # devolve
        return r

    def get_jobs(self, status: object = None) -> _Jobs:
        # identificadores dos pedidos abertos
        return _Jobs(list(self.abertos))

    def get_remote(self, i: str) -> RemotoFalso:
        # pedido aberto pelo identificador
        return self.abertos[i]


#: Datas de um bloco de teste.
INI, FIM = pd.Timestamp("2015-01-01").date(), pd.Timestamp("2015-01-10").date()


def _baixar(cliente: ClienteFalso, pasta: Path, **kw: float) -> Path | None:
    """Chama baixar_bloco sem esperas reais."""
    # padrões: limite rodando minúsculo, fila folgada e consultas imediatas
    opcoes = {"limite_h": 1e-6, "limite_fila_h": 1.0, "intervalo_s": 0.0, **kw}
    # baixa o bloco
    return baixar_bloco(
        cliente,
        INI,
        FIM,
        "00:00",
        ["particulate_matter_2.5um"],
        [0, 0, 0, 0],
        range(3),
        pasta,
        tentativas=5,
        espera_min=0,
        **opcoes,
    )


def test_fita_indisponivel_pula_sem_repetir(tmp_path: Path) -> None:
    """Erro de fita indisponível: uma única tentativa, bloco anotado e pulado."""
    # pedido que falha com o texto real do MARS
    cliente = ClienteFalso(["failed"], f"AccessError: The {FITA_INDISPONIVEL}. Please visit ...")
    # bloco pulado
    assert _baixar(cliente, tmp_path) is None
    # sem insistir na mesma fita
    assert len(cliente.remotos) == 1
    # anotado para a próxima execução
    assert (tmp_path / "indisponiveis.txt").read_text().startswith("2015-01-01 2015-01-10")


def test_pedido_rodando_demais_cancela_e_pula(tmp_path: Path) -> None:
    """Pedido que roda sem terminar: cancelado no ADS, tentado 2 vezes, depois pulado."""
    # pedido eternamente rodando
    cliente = ClienteFalso(["running"])
    # bloco pulado
    assert _baixar(cliente, tmp_path) is None
    # duas tentativas (tentativas_demora padrão)
    assert len(cliente.remotos) == 2
    # cada pedido demorado foi cancelado no ADS
    assert all(r.apagado for r in cliente.remotos)
    # motivo anotado
    assert "demora" in (tmp_path / "indisponiveis.txt").read_text()


def test_fila_longa_nao_conta_no_limite_rodando(tmp_path: Path) -> None:
    """Muito tempo na fila não estoura o limite de processamento; o bloco é baixado."""
    # várias consultas na fila e depois sucesso, com o limite rodando minúsculo
    cliente = ClienteFalso(["accepted"] * 50 + ["successful"])
    # baixa normalmente
    destino = _baixar(cliente, tmp_path)
    # bloco completo
    assert destino is not None and (destino / ".completo").exists()
    # um único pedido, nunca cancelado
    assert len(cliente.remotos) == 1 and not cliente.remotos[0].apagado


def test_retoma_pedido_identico_na_fila(tmp_path: Path) -> None:
    """Se já existe um pedido idêntico na fila, ele é retomado em vez de criar outro."""
    # cliente cujos pedidos novos dariam sucesso
    cliente = ClienteFalso(["successful"])
    # pedido igual ao que o bloco faria, já aberto no ADS (área com inteiros, como no JSON)
    req = pedido(INI, FIM, "00:00", ["particulate_matter_2.5um"], [0, 0, 0, 0], range(3))
    # registra o pedido aberto
    cliente.abertos["x"] = RemotoFalso(["running", "successful"], request=dict(req))
    # baixa
    assert _baixar(cliente, tmp_path, limite_h=1.0) is not None
    # nenhum pedido novo foi criado
    assert cliente.remotos == []


def test_outro_erro_esgota_tentativas_e_para(tmp_path: Path) -> None:
    """Erro desconhecido (por exemplo, credencial): tenta e depois para o download."""
    # pedido rejeitado por outro motivo
    cliente = ClienteFalso(["rejected"], "licença não aceita")
    # o erro sobe depois das tentativas
    with pytest.raises(RuntimeError, match="licença"):
        _baixar(cliente, tmp_path)
    # usou todas as tentativas
    assert len(cliente.remotos) == 5


def test_media_de_pontos_em_campo_linear(tmp_path: Path) -> None:
    """Num campo que varia linearmente com a longitude, a média dos pontos é o valor no centro."""
    # campo constante em todas as horas, para gravar a grade
    _netcdf(tmp_path, 1e-8, np.arange(0, 103), ["2020-01-01"])
    # reabre e substitui por um campo linear em longitude: valor = (lon + 50) * 1e-9 kg/m³
    ds = xr.open_dataset(tmp_path / "a.nc").load()
    # campo linear
    ds["pm2p5"].values[:] = ((ds["longitude"] + 50) * 1e-9).values
    # regrava
    ds.to_netcdf(tmp_path / "b.nc")
    # remove o original
    (tmp_path / "a.nc").unlink()
    # dois pontos simétricos em volta da estação
    t = preparar(tmp_path, [LAT, LAT], [LON - 0.1, LON + 0.1])
    # a interpolação bilinear é exata num campo linear; a média é o valor no centro
    np.testing.assert_allclose(t["fc_pm2_5_l0"].iloc[0], LON + 50, rtol=1e-4)
