"""UMI 分子共识去噪 API（FastAPI 入口）。"""

from __future__ import annotations

import json
from typing import Any

from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from . import db, service
from .errors import APIError
from .parsing import parse_submission


@asynccontextmanager
async def lifespan(_: FastAPI):
    db.init_db()
    yield


app = FastAPI(
    title="UMI 分子共识去噪 API",
    description=(
        "按 UMI 频次/字典序分簇、汉明距离 <=1 归并，簇内以 Phred 权重 "
        "10^(Q/10) 逐位投票生成不可变版本化共识。"
    ),
    version="1.0.0",
    lifespan=lifespan,
)


@app.exception_handler(APIError)
async def _api_error_handler(_: Request, exc: APIError) -> JSONResponse:
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": {
                "code": exc.code,
                "message": exc.message,
                "details": exc.details,
            }
        },
    )


@app.get("/health", tags=["meta"])
def health() -> dict:
    return {"status": "ok"}


@app.post("/consensus", tags=["consensus"])
async def create_consensus(request: Request) -> JSONResponse:
    """提交读数并生成一个不可变共识版本。

    同一 request_id + 同一载荷为幂等重放（200）；不同载荷 409；
    同一 read_id 内容冲突 409；首次提交 201。
    """
    raw = await request.body()
    if not raw:
        body: Any = None
    else:
        try:
            body = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            return JSONResponse(
                status_code=400,
                content={
                    "error": {
                        "code": "invalid_json",
                        "message": "请求体不是合法的 JSON",
                        "details": {},
                    }
                },
            )

    parsed = parse_submission(body)
    result, replayed = service.submit(parsed)
    return JSONResponse(
        status_code=200 if replayed else 201,
        content=result,
        headers={"X-Idempotent-Replay": "true" if replayed else "false"},
    )


@app.get("/consensus/{version}", tags=["consensus"])
def read_consensus(version: int) -> dict:
    """按版本号读取不可变共识快照。"""
    return service.get_version(version)


@app.get("/versions", tags=["consensus"])
def list_versions(request_id: str | None = None) -> dict:
    """列出全部（或指定 request_id 的）版本元数据。"""
    return service.list_versions(request_id)
