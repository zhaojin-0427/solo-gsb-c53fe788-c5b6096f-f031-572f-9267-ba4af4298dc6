"""SQLite 存储：提交、读数、不可变版本快照。

写路径在同一把进程锁内用 BEGIN IMMEDIATE 串行化；读路径为单条只读
SELECT，配合 WAL，可与写入并发并始终读到已提交的一致快照。
"""

from __future__ import annotations

import json
import os
import sqlite3
import threading
from contextlib import contextmanager
from typing import Iterator, Optional

_SCHEMA = """
CREATE TABLE IF NOT EXISTS submissions (
    request_id    TEXT PRIMARY KEY,
    sample_id     TEXT NOT NULL,
    payload_hash  TEXT NOT NULL,
    version       INTEGER NOT NULL,
    created_at    TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS read_records (
    read_id        TEXT PRIMARY KEY,
    locus          TEXT NOT NULL,
    umi            TEXT NOT NULL,
    sequence       TEXT NOT NULL,
    qualities_json TEXT NOT NULL,
    created_at     TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS submission_reads (
    request_id TEXT NOT NULL,
    read_id    TEXT NOT NULL,
    PRIMARY KEY (request_id, read_id),
    FOREIGN KEY (request_id) REFERENCES submissions(request_id),
    FOREIGN KEY (read_id) REFERENCES read_records(read_id)
);

CREATE TABLE IF NOT EXISTS versions (
    version       INTEGER PRIMARY KEY AUTOINCREMENT,
    request_id    TEXT NOT NULL,
    sample_id     TEXT NOT NULL,
    payload_hash  TEXT NOT NULL,
    content_hash  TEXT NOT NULL,
    result_json   TEXT NOT NULL,
    created_at    TEXT NOT NULL
);
"""

_write_lock = threading.RLock()


def db_path() -> str:
    return os.environ.get("UMI_DB_PATH", os.path.join("data", "umi.db"))


@contextmanager
def _connect(write: bool = False) -> Iterator[sqlite3.Connection]:
    path = db_path()
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    conn = sqlite3.connect(path, timeout=30, isolation_level=None)
    try:
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("PRAGMA busy_timeout=30000")
        if write:
            conn.execute("BEGIN IMMEDIATE")
        yield conn
        if write:
            conn.execute("COMMIT")
    except Exception:
        if write:
            conn.execute("ROLLBACK")
        raise
    finally:
        conn.close()


def init_db() -> None:
    # DDL 在自动提交模式下执行，避免 executescript 隐式提交事务。
    with _connect(write=False) as conn:
        conn.executescript(_SCHEMA)


@contextmanager
def write_transaction() -> Iterator[sqlite3.Connection]:
    """进程级串行化的写事务。"""
    with _write_lock:
        with _connect(write=True) as conn:
            yield conn


@contextmanager
def read_connection() -> Iterator[sqlite3.Connection]:
    with _connect(write=False) as conn:
        yield conn


def get_submission(conn: sqlite3.Connection, request_id: str) -> Optional[sqlite3.Row]:
    cur = conn.execute(
        "SELECT * FROM submissions WHERE request_id = ?", (request_id,)
    )
    return cur.fetchone()


def find_read(conn: sqlite3.Connection, read_id: str) -> Optional[sqlite3.Row]:
    cur = conn.execute(
        "SELECT * FROM read_records WHERE read_id = ?", (read_id,)
    )
    return cur.fetchone()


def insert_submission(
    conn: sqlite3.Connection,
    request_id: str,
    sample_id: str,
    payload_hash: str,
    version: int,
    created_at: str,
) -> None:
    conn.execute(
        "INSERT INTO submissions(request_id, sample_id, payload_hash, version, created_at)"
        " VALUES (?, ?, ?, ?, ?)",
        (request_id, sample_id, payload_hash, version, created_at),
    )


def upsert_read(conn: sqlite3.Connection, read, created_at: str) -> None:
    conn.execute(
        "INSERT OR IGNORE INTO read_records"
        "(read_id, locus, umi, sequence, qualities_json, created_at)"
        " VALUES (?, ?, ?, ?, ?, ?)",
        (
            read.read_id,
            read.locus,
            read.umi,
            read.sequence,
            json.dumps(list(read.qualities), separators=(",", ":")),
            created_at,
        ),
    )


def link_submission_read(conn: sqlite3.Connection, request_id: str, read_id: str) -> None:
    conn.execute(
        "INSERT OR IGNORE INTO submission_reads(request_id, read_id) VALUES (?, ?)",
        (request_id, read_id),
    )


def next_version(conn: sqlite3.Connection) -> int:
    row = conn.execute("SELECT COALESCE(MAX(version), 0) + 1 AS v FROM versions").fetchone()
    return int(row["v"])


def insert_version(
    conn: sqlite3.Connection,
    version: int,
    request_id: str,
    sample_id: str,
    payload_hash: str,
    content_hash: str,
    result_json: str,
    created_at: str,
) -> None:
    conn.execute(
        "INSERT INTO versions(version, request_id, sample_id, payload_hash,"
        " content_hash, result_json, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            version,
            request_id,
            sample_id,
            payload_hash,
            content_hash,
            result_json,
            created_at,
        ),
    )


def get_version_row(conn: sqlite3.Connection, version: int) -> Optional[sqlite3.Row]:
    cur = conn.execute("SELECT * FROM versions WHERE version = ?", (version,))
    return cur.fetchone()


def list_version_rows(conn: sqlite3.Connection, request_id: Optional[str] = None):
    if request_id is None:
        cur = conn.execute(
            "SELECT version, request_id, sample_id, payload_hash, content_hash,"
            " created_at FROM versions ORDER BY version"
        )
    else:
        cur = conn.execute(
            "SELECT version, request_id, sample_id, payload_hash, content_hash,"
            " created_at FROM versions WHERE request_id = ? ORDER BY version",
            (request_id,),
        )
    return cur.fetchall()


def latest_version(conn: sqlite3.Connection) -> Optional[int]:
    row = conn.execute("SELECT MAX(version) AS v FROM versions").fetchone()
    v = row["v"]
    return int(v) if v is not None else None
