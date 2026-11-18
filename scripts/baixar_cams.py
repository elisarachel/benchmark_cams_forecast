"""Baixa as previsões globais do CAMS do Atmosphere Data Store (ADS), um bloco de dias por vez.

Dataset: ``cams-global-atmospheric-composition-forecasts``. Os nomes de chaves e valores
abaixo foram conferidos no formulário e nas restrições (constraints) do dataset no ADS em
2026-11-05: ``variable`` (por exemplo ``particulate_matter_2.5um``), ``date`` no formato
``AAAA-MM-DD/AAAA-MM-DD``, ``time`` (``00:00`` ou ``12:00``), ``leadtime_hour`` (0 a 120,
horário para o MP2,5 desde 2015-01-01), ``type`` (``forecast``), ``data_format``
(``grib`` ou ``netcdf_zip``) e ``area`` como ``[norte, oeste, sul, leste]``.

Uso::

    uv run python scripts/baixar_cams.py --inicio 2015-01-01 --fim 2024-12-31

Por padrão são pedidas só as horas de previsão 4 a 99, as únicas usadas nas médias diárias
dos leads 0 a 3 (ver ``preparar_cams.horas_previsao``).

As credenciais ficam em ``~/.cdsapirc`` (nunca versionar). Para usar outra conta, aponte a
variável de ambiente ``CDSAPI_RC`` para outro arquivo de credenciais. O script é retomável:
blocos já baixados são pulados.
"""

from __future__ import annotations

# argparse: leitura dos parâmetros de linha de comando
import argparse

# calendar: número de dias de cada mês
import calendar

# csv: leitura da tabela de estações do benchmark
import csv

# json: comparação de pedidos (forma canônica)
import json

# math: arredondamento da caixa para fora
import math

# os: variável de ambiente CDSAPI_RC (arquivo de credenciais)
import os

# time: espera entre tentativas
import time

# zipfile: o ADS entrega o NetCDF dentro de um ZIP
import zipfile

# date: datas de início e fim
from datetime import date

# Path: caminhos de arquivo
from pathlib import Path

#: Identificador do dataset no ADS.
DATASET = "cams-global-atmospheric-composition-forecasts"

#: Tabela de estações do benchmark (fonte única das coordenadas).
ESTACOES = Path(__file__).resolve().parents[2] / "Benchmark/src/pm25bench/data/stations.csv"

#: Nome curto de cada variável do ADS, usado no nome do arquivo.
CURTOS: dict[str, str] = {
    # MP2,5
    "particulate_matter_2.5um": "pm2p5",
    # MP10
    "particulate_matter_10um": "pm10",
    # dióxido de nitrogênio
    "nitrogen_dioxide": "no2",
    # ozônio
    "ozone": "go3",
    # monóxido de carbono
    "carbon_monoxide": "co",
    # dióxido de enxofre
    "sulphur_dioxide": "so2",
}

#: Trecho do erro do MARS quando a fita com os dados está fora do ar (não adianta repetir).
FITA_INDISPONIVEL = "tape on which the data reside is unavailable"

#: Resolução da grade do CAMS entregue pelo ADS, em graus.
RESOLUCAO = 0.4


def coordenadas_estacao(codigo: int, tabela: Path = ESTACOES) -> tuple[float, float]:
    """Latitude e longitude de uma estação, lidas de ``stations.csv`` do benchmark."""
    # abre a tabela de estações
    with open(tabela, encoding="utf-8") as f:
        # percorre as linhas como dicionários
        for linha in csv.DictReader(f):
            # encontrou a estação pedida
            if int(linha["codigo"]) == codigo:
                # devolve (lat, lon) como float
                return float(linha["lat"]), float(linha["lon"])
    # estação ausente é erro do usuário
    raise KeyError(f"estação {codigo} não está em {tabela}")


def area_padrao(lat: float, lon: float, celulas: int = 2) -> list[float]:
    """Caixa ``[N, W, S, E]`` em volta do ponto com ``celulas`` células de margem por lado.

    Com 2 células de 0,4° (0,8°) de margem, os quatro vizinhos da interpolação bilinear
    ficam dentro da caixa com folga, qualquer que seja a posição do ponto na grade. Os
    limites são arredondados para fora, em décimos de grau.
    """
    # margem em graus
    m = celulas * RESOLUCAO
    # arredonda para cima (norte e leste) e para baixo (sul e oeste)
    return [
        # norte
        math.ceil((lat + m) * 10) / 10,
        # oeste
        math.floor((lon - m) * 10) / 10,
        # sul
        math.floor((lat - m) * 10) / 10,
        # leste
        math.ceil((lon + m) * 10) / 10,
    ]


def meses(inicio: date, fim: date) -> list[tuple[date, date]]:
    """Lista de (primeiro dia, último dia) de cada mês entre ``inicio`` e ``fim``, recortada."""
    # resultado
    saida = []
    # começa no primeiro dia do mês de início
    atual = date(inicio.year, inicio.month, 1)
    # percorre mês a mês
    while atual <= fim:
        # primeiro dia do mês seguinte
        proximo = date(atual.year + atual.month // 12, atual.month % 12 + 1, 1)
        # último dia do mês atual
        ultimo = date.fromordinal(proximo.toordinal() - 1)
        # recorta pelo intervalo pedido
        saida.append((max(atual, inicio), min(ultimo, fim)))
        # avança
        atual = proximo
    # devolve a lista
    return saida


def blocos(inicio: date, fim: date) -> list[tuple[date, date]]:
    """Divide cada mês em blocos fixos de dias (1-10, 11-20, 21-fim), recortados.

    Pedidos menores têm mais chance de terminar no ADS: dados com mais de 30 dias vêm do
    arquivo em fita, e um pedido de um mês inteiro com 103 horas foi derrubado depois de 6 h
    rodando. Os limites são fixos no calendário, para que o nome de cada bloco não dependa
    de ``--inicio`` e o download continue retomável.
    """
    # resultado
    saida = []
    # percorre os meses do intervalo (só o primeiro dia de cada mês interessa)
    for ini_mes, _ in meses(inicio, fim):
        # número de dias do mês (considera anos bissextos)
        dias_mes = calendar.monthrange(ini_mes.year, ini_mes.month)[1]
        # limites fixos: 1-10, 11-20, 21-fim do mês
        for d0, d1 in [(1, 10), (11, 20), (21, dias_mes)]:
            # bloco recortado pelo intervalo pedido
            a = max(date(ini_mes.year, ini_mes.month, d0), inicio)
            # fim do bloco recortado
            b = min(date(ini_mes.year, ini_mes.month, d1), fim)
            # guarda só blocos não vazios
            if a <= b:
                saida.append((a, b))
    # devolve a lista
    return saida


def nome_arquivo(inicio: date, fim: date, rodada: str, variaveis: list[str]) -> str:
    """Nome base de um bloco: ``cams_fc_AAAA-MM-DD_AAAA-MM-DD_HHUTC_<variaveis>``."""
    # hora da rodada sem os minutos
    hora = rodada.split(":")[0]
    # nomes curtos das variáveis, em ordem fixa
    curtos = "+".join(sorted(CURTOS.get(v, v) for v in variaveis))
    # junta as partes
    return f"cams_fc_{inicio:%Y-%m-%d}_{fim:%Y-%m-%d}_{hora}UTC_{curtos}"


def pedido(
    inicio: date, fim: date, rodada: str, variaveis: list[str], area: list[float], horas: range
) -> dict[str, object]:
    """Monta o dicionário do pedido ao ADS (chaves conferidas no 'Show API request')."""
    # dicionário no formato aceito pelo cdsapi
    return {
        # variáveis pedidas
        "variable": variaveis,
        # intervalo de datas de emissão (data base do modelo)
        "date": [f"{inicio:%Y-%m-%d}/{fim:%Y-%m-%d}"],
        # hora da rodada
        "time": [rodada],
        # horas de previsão como texto
        "leadtime_hour": [str(h) for h in horas],
        # campos de previsão (não análise)
        "type": ["forecast"],
        # NetCDF compactado
        "data_format": "netcdf_zip",
        # caixa [N, W, S, E]
        "area": area,
    }


def ler_credenciais() -> tuple[str, str]:
    """URL e chave do ADS, lidas de ``$CDSAPI_RC`` (ou ``~/.cdsapirc``), no formato do cdsapi."""
    # arquivo de credenciais: o da variável de ambiente ou o padrão
    arquivo = Path(os.environ.get("CDSAPI_RC", Path.home() / ".cdsapirc"))
    # pares "chave: valor" de cada linha
    campos = {}
    # lê o arquivo linha a linha
    for linha in arquivo.read_text(encoding="utf-8").splitlines():
        # só linhas com ":" (separa no primeiro, a URL também tem ":")
        if ":" in linha:
            # chave e valor sem espaços
            k, v = linha.split(":", 1)
            # guarda
            campos[k.strip()] = v.strip()
    # devolve url e chave
    return campos["url"], campos["key"]


class PedidoDemorado(Exception):
    """O pedido passou do limite de tempo no ADS e foi cancelado."""


def baixar_bloco(
    cliente: object,
    inicio: date,
    fim: date,
    rodada: str,
    variaveis: list[str],
    area: list[float],
    horas: range,
    saida: Path,
    tentativas: int = 6,
    espera_min: float = 15.0,
    limite_h: float = 3.0,
    tentativas_demora: int = 2,
    intervalo_s: float = 60.0,
    limite_fila_h: float = 24.0,
) -> Path | None:
    """Baixa um bloco, tentando de novo com espera crescente se o ADS falhar.

    O ADS pode rejeitar, derrubar ou apagar um pedido (por exemplo, "job not found" depois
    de horas rodando). Em vez de parar o download inteiro, o bloco é pedido de novo depois
    de ``espera_min`` minutos, dobrando a espera a cada falha (até 4 h). Os pedidos
    continuam sequenciais: nunca há mais de um pedido deste script na fila.

    Dois casos não insistem indefinidamente. Se o MARS responder que a fita com os dados
    está indisponível, tentar de novo não adianta (é a mesma fita). Se o pedido passar de
    ``limite_h`` horas sem terminar ``tentativas_demora`` vezes, o bloco provavelmente está
    numa fita com problema (``limite_h`` conta só o tempo rodando, não a fila). Nos dois
    casos o bloco é anotado em ``saida/indisponiveis.txt``,
    o download segue para o próximo e a função devolve ``None``. Como o bloco não ganha o
    marcador ``.completo``, uma nova execução do script tenta esse bloco outra vez.
    """
    # espera atual, em minutos
    espera = espera_min
    # quantas vezes o pedido passou do limite de tempo
    demoras = 0
    # percorre as tentativas
    for n in range(1, tentativas + 1):
        # tenta baixar
        try:
            # sucesso: devolve a pasta
            return baixar_mes(
                cliente,
                inicio,
                fim,
                rodada,
                variaveis,
                area,
                horas,
                saida,
                limite_h,
                intervalo_s,
                limite_fila_h,
            )
        # qualquer erro do ADS ou da rede
        except Exception as erro:
            # conta as demoras excessivas
            demoras += isinstance(erro, PedidoDemorado)
            # fita indisponível, ou demoras repetidas: não insiste, anota e segue
            if FITA_INDISPONIVEL in str(erro) or demoras >= tentativas_demora:
                # motivo curto
                motivo = "fita indisponível" if demoras < tentativas_demora else "demora"
                # mensagem
                print(f"{inicio}..{fim}: {motivo} no ECMWF; bloco pulado")
                # anota o bloco para a próxima execução e para o README
                with open(saida / "indisponiveis.txt", "a", encoding="utf-8") as f:
                    # uma linha por ocorrência: início, fim, motivo e quando
                    f.write(f"{inicio} {fim} {motivo} {time.strftime('%Y-%m-%dT%H:%M')}\n")
                # segue sem esse bloco
                return None
            # sem mais tentativas: repassa o erro e o laço principal para
            if n == tentativas:
                raise
            # registra a falha
            print(f"{inicio}..{fim}: tentativa {n} falhou ({erro}); nova em {espera:.0f} min")
            # espera antes de pedir de novo
            time.sleep(espera * 60)
            # dobra a espera, com teto de 4 h
            espera = min(espera * 2, 240.0)
    # inalcançável: o laço devolve ou levanta
    raise RuntimeError("tentativas esgotadas")


def esperar_pedido(
    remoto: object, limite_h: float, intervalo_s: float, limite_fila_h: float = 24.0
) -> None:
    """Espera o pedido terminar; cancela e levanta ``PedidoDemorado`` se passar do limite.

    São dois limites separados. ``limite_h`` conta só o tempo **rodando** (status
    ``running``): é o que indica um problema na leitura da fita, já que um bloco normal roda
    em 10 a 30 min. ``limite_fila_h`` conta o tempo esperando na fila (``accepted``), que
    depende do movimento do ADS e pode passar de 2 h sem nada de errado com o pedido.

    Levanta ``RuntimeError`` com a mensagem do ADS (que inclui o erro do MARS) se o pedido
    falhar, for rejeitado ou apagado.
    """
    # início da espera na fila
    t_fila = time.monotonic()
    # início do processamento (None enquanto não começar a rodar)
    t_rodando: float | None = None
    # consulta o status até um estado final
    while True:
        # status atual do pedido ("accepted", "running", "successful", "failed", ...)
        status = remoto.status  # type: ignore[attr-defined]
        # pronto para baixar
        if status == "successful":
            return
        # estados de erro: busca a mensagem completa do ADS
        if status in ("failed", "rejected", "dismissed", "deleted"):
            # get_results levanta um erro com o motivo (por exemplo, a fita indisponível)
            try:
                remoto.get_results()  # type: ignore[attr-defined]
            # repassa a mensagem como texto
            except Exception as erro:
                raise RuntimeError(f"pedido {status}: {erro}") from erro
            # sem mensagem detalhada
            raise RuntimeError(f"pedido {status}")
        # agora
        agora = time.monotonic()
        # marca o início do processamento na primeira vez que o pedido aparece rodando
        if status == "running" and t_rodando is None:
            t_rodando = agora
        # o pedido pode voltar para a fila (o ADS às vezes o recoloca): zera o relógio
        if status == "accepted" and t_rodando is not None:
            t_rodando, t_fila = None, agora
        # tempo rodando acima do limite, ou fila acima do teto
        rodando_demais = t_rodando is not None and agora - t_rodando > limite_h * 3600
        # fila acima do teto
        fila_demais = t_rodando is None and agora - t_fila > limite_fila_h * 3600
        # passou de algum limite: cancela no ADS para não ocupar a fila
        if rodando_demais or fila_demais:
            # cancela o pedido
            remoto.delete()  # type: ignore[attr-defined]
            # descrição do limite que estourou
            qual = (
                f"mais de {limite_h:g} h rodando"
                if rodando_demais
                else (f"mais de {limite_fila_h:g} h na fila")
            )
            # sinaliza a demora
            raise PedidoDemorado(f"{qual}; pedido cancelado")
        # espera antes de consultar de novo
        time.sleep(intervalo_s)


def _canonico(req: dict[str, object]) -> str:
    """Forma canônica de um pedido, para comparar pedidos sem depender da ordem das chaves."""
    # JSON com chaves ordenadas e números como float
    return json.dumps(
        {k: [float(x) for x in v] if k == "area" else v for k, v in req.items()}, sort_keys=True
    )


def pedido_na_fila(cliente: object, req: dict[str, object]) -> object | None:
    """Pedido idêntico a ``req`` que já esteja na fila ou rodando no ADS, ou ``None``.

    Permite reiniciar o script sem perder o lugar na fila: em vez de criar um pedido novo
    (e voltar para o fim da fila), o script retoma o que já existe.
    """
    # pedidos ainda não terminados da conta
    ids = cliente.get_jobs(status=["accepted", "running"]).request_ids  # type: ignore[attr-defined]
    # forma canônica do pedido procurado
    alvo = _canonico(req)
    # percorre os pedidos abertos
    for i in ids:
        # pedido remoto
        r = cliente.get_remote(i)  # type: ignore[attr-defined]
        # mesmo conteúdo: retoma este
        if _canonico(r.request) == alvo:
            return r
    # nenhum pedido igual na fila
    return None


def baixar_mes(
    cliente: object,
    inicio: date,
    fim: date,
    rodada: str,
    variaveis: list[str],
    area: list[float],
    horas: range,
    saida: Path,
    limite_h: float = 3.0,
    intervalo_s: float = 60.0,
    limite_fila_h: float = 24.0,
) -> Path:
    """Baixa um intervalo de datas e extrai os NetCDF em ``saida/<nome>/``. Pula se já existir."""
    # nome base do bloco
    nome = nome_arquivo(inicio, fim, rodada, variaveis)
    # pasta final com os NetCDF
    destino = saida / nome
    # marcador de conclusão: só existe se o download e a extração terminaram
    pronto = destino / ".completo"
    # retomável: bloco já baixado
    if pronto.exists():
        # avisa e devolve
        print(f"{nome}: já baixado")
        return destino
    # arquivo temporário do download
    zip_parcial = saida / f"{nome}.zip.part"
    # mostra o pedido
    print(f"{nome}: pedindo {inicio} a {fim}")
    # conteúdo do pedido
    req = pedido(inicio, fim, rodada, variaveis, area, horas)
    # retoma um pedido idêntico já na fila (por exemplo, depois de reiniciar o script)
    remoto = pedido_na_fila(cliente, req)
    # sem pedido na fila: envia um novo
    if remoto is None:
        # novo pedido na fila do ADS
        remoto = cliente.submit(DATASET, req)  # type: ignore[attr-defined]
        # identifica o pedido no log
        print(f"{nome}: pedido {remoto.request_id}")  # type: ignore[attr-defined]
    else:
        # identifica o pedido retomado no log
        print(f"{nome}: retomando pedido {remoto.request_id} já na fila")  # type: ignore[attr-defined]
    # espera terminar, com limites de tempo
    esperar_pedido(remoto, limite_h, intervalo_s, limite_fila_h)
    # baixa o arquivo pronto
    remoto.download(str(zip_parcial))
    # cria a pasta do bloco
    destino.mkdir(parents=True, exist_ok=True)
    # extrai os NetCDF
    with zipfile.ZipFile(zip_parcial) as z:
        # todos os membros do ZIP
        z.extractall(destino)
    # remove o ZIP depois da extração
    zip_parcial.unlink()
    # marca o bloco como completo
    pronto.touch()
    # print de conclusão
    print(f"{nome}: concluído")
    # devolve a pasta
    return destino


def main() -> None:
    """Lê os argumentos e baixa bloco a bloco."""
    # analisador de argumentos
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    # data inicial (data de emissão)
    p.add_argument("--inicio", type=date.fromisoformat, required=True)
    # data final (data de emissão)
    p.add_argument("--fim", type=date.fromisoformat, required=True)
    # hora da rodada
    p.add_argument("--rodada", default="00:00", choices=["00:00", "12:00"])
    # variáveis do ADS
    p.add_argument("--variaveis", nargs="+", default=["particulate_matter_2.5um"])
    # estação usada para a caixa padrão
    p.add_argument("--estacao", type=int, default=73)
    # caixa explícita, sobrescreve a padrão
    p.add_argument("--area", type=float, nargs=4, metavar=("N", "W", "S", "E"))
    # primeira hora de previsão pedida (a primeira usada no lead 0 é a 4)
    p.add_argument("--hora-min", type=int, default=4)
    # última hora de previsão pedida (a última usada no lead 3 é a 99)
    p.add_argument("--hora-max", type=int, default=99)
    # pasta de saída
    p.add_argument("--saida", type=Path, default=Path("raw"))
    # tentativas por bloco antes de desistir
    p.add_argument("--tentativas", type=int, default=6)
    # espera inicial entre tentativas, em minutos (dobra a cada falha)
    p.add_argument("--espera-min", type=float, default=15.0)
    # tempo máximo de um pedido rodando no ADS, em horas (a fila não conta)
    p.add_argument("--limite-h", type=float, default=3.0)
    # tempo máximo de um pedido esperando na fila do ADS, em horas
    p.add_argument("--limite-fila-h", type=float, default=24.0)
    # lê a linha de comando
    a = p.parse_args()
    # caixa: a informada ou a padrão em volta da estação
    area = a.area or area_padrao(*coordenadas_estacao(a.estacao))
    # horas de previsão pedidas
    horas = range(a.hora_min, a.hora_max + 1)
    # garante a pasta de saída
    a.saida.mkdir(parents=True, exist_ok=True)
    # importação tardia: os testes não precisam do cliente do ADS
    from ecmwf.datastores import Client

    # cliente com as credenciais de $CDSAPI_RC (ou ~/.cdsapirc)
    url, chave = ler_credenciais()
    # cliente do ADS (submit/status/delete permitem o limite de tempo)
    cliente = Client(url=url, key=chave)
    # mostra a caixa e as horas usadas
    print(f"área [N, W, S, E] = {area}; horas {a.hora_min} a {a.hora_max}")
    # um pedido por bloco de ~10 dias, em sequência: o ADS rejeita pedidos acima do limite de
    # fila por usuário e pede explicitamente scripts sequenciais
    for ini, fim in blocos(a.inicio, a.fim):
        # baixa (ou pula) o bloco, com novas tentativas se o ADS falhar
        baixar_bloco(
            cliente,
            ini,
            fim,
            a.rodada,
            a.variaveis,
            area,
            horas,
            a.saida,
            a.tentativas,
            a.espera_min,
            a.limite_h,
            limite_fila_h=a.limite_fila_h,
        )


# executa só quando chamado como script
if __name__ == "__main__":
    main()
