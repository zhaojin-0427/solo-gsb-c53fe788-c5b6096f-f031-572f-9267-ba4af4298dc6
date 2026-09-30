"""提交与版本读取的业务编排。"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Optional

from . import db
from .consensus import build_consensus_result
from .errors import NotFound, PayloadConflictError
from .parsing import ParsedSubmission, canonical_hash


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _read_content_tuple(row) -> tuple:
    return (
        row["locus"],
        row["umi"],
        row["sequence"],
        tuple(json.loads(row["qualities_json"])),
    )


def submit(parsed: ParsedSubmission) -> tuple[dict, bool]:
    """处理一次提交，返回 (版本结果, 是否为幂等重放)。

    冲突约定：
    - 同一 request_id 但不同载荷 -> 409 payload_conflict；
    - 同一 read_id 但内容（locus/umi/sequence/qualities）不同 -> 409 read_conflict；
    - 同 request_id 同载荷 -> 幂等，返回已存在的同一版本（200）。
    """
    payload_hash = canonical_hash(parsed.canonical_payload)

    with db.write_transaction() as conn:
        existing = db.get_submission(conn, parsed.request_id)
        if existing is not None:
            if existing["payload_hash"] != payload_hash:
                raise PayloadConflictError(
                    "payload_conflict",
                    f"request_id {parsed.request_id} 已存在但请求载荷不同",
                    {
                        "request_id": parsed.request_id,
                        "existing_version": existing["version"],
                    },
                )
            row = db.get_version_row(conn, existing["version"])
            return json.loads(row["result_json"]), True

        # 新 request_id：先检查 read_id 是否与历史内容冲突。
        for read in parsed.reads:
            prior = db.find_read(conn, read.read_id)
            if prior is None:
                continue
            if _read_content_tuple(prior) != (
                read.locus,
                read.umi,
                read.sequence,
                tuple(read.qualities),
            ):
                raise PayloadConflictError(
                    "read_conflict",
                    f"read_id {read.read_id} 已存在但内容与本次提交冲突",
                    {
                        "read_id": read.read_id,
                        "existing_locus": prior["locus"],
                        "incoming_locus": read.locus,
                    },
                )

        result = build_consensus_result(
            parsed.request_id, parsed.sample_id, parsed.reads
        )
        created_at = utc_now_iso()
        version = db.next_version(conn)
        content_hash = canonical_hash(
            {
                "request_id": result["request_id"],
                "sample_id": result["sample_id"],
                "consensus": result["consensus"],
                "clusters": result["clusters"],
                "eliminated": result["eliminated"],
                "stats": result["stats"],
            }
        )
        result["version"] = version
        result["content_hash"] = content_hash
        result["created_at"] = created_at

        result_json = json.dumps(result, ensure_ascii=False, separators=(",", ":"))

        db.insert_submission(
            conn,
            parsed.request_id,
            parsed.sample_id,
            payload_hash,
            version,
            created_at,
        )
        for read in parsed.reads:
            db.upsert_read(conn, read, created_at)
            db.link_submission_read(conn, parsed.request_id, read.read_id)
        db.insert_version(
            conn,
            version,
            parsed.request_id,
            parsed.sample_id,
            payload_hash,
            content_hash,
            result_json,
            created_at,
        )
        return result, False


def get_version(version: int) -> dict:
    with db.read_connection() as conn:
        row = db.get_version_row(conn, version)
        if row is None:
            raise NotFound(
                "version_not_found",
                f"版本 {version} 不存在",
                {"version": version},
            )
        return json.loads(row["result_json"])


def list_versions(request_id: Optional[str] = None) -> dict:
    with db.read_connection() as conn:
        rows = db.list_version_rows(conn, request_id)
        latest = db.latest_version(conn)
    return {
        "latest_version": latest,
        "count": len(rows),
        "versions": [
            {
                "version": r["version"],
                "request_id": r["request_id"],
                "sample_id": r["sample_id"],
                "payload_hash": r["payload_hash"],
                "content_hash": r["content_hash"],
                "created_at": r["created_at"],
            }
            for r in rows
        ],
    }
