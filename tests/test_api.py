"""API 层测试：校验、幂等、409 冲突、版本不可变、读取一致。"""

import os
import tempfile

import pytest
from fastapi.testclient import TestClient


@pytest.fixture()
def client():
    tmp = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
    tmp.close()
    os.environ["UMI_DB_PATH"] = tmp.name

    from app import db
    from app.main import app

    db.init_db()
    with TestClient(app) as c:
        yield c
    for suffix in ("", "-wal", "-shm"):
        try:
            os.unlink(tmp.name + suffix)
        except OSError:
            pass


def _rec(rid, locus="L1", umi="AAAA", seq="AAAA", quals=None):
    return {
        "read_id": rid,
        "locus": locus,
        "umi": umi,
        "sequence": seq,
        "qualities": [40] * len(seq) if quals is None else quals,
    }


def _envelope(records, request_id="req1", sample_id="s1"):
    return {
        "request_id": request_id,
        "sample_id": sample_id,
        "reads": records,
    }


def test_create_and_get_version(client):
    payload = _envelope(
        [_rec("r1", seq="AAAA"), _rec("r2", seq="AAAA"), _rec("r3", seq="AAAA")]
    )
    r = client.post("/consensus", json=payload)
    assert r.status_code == 201
    body = r.json()
    assert body["version"] == 1
    assert body["consensus"][0]["consensus"] == "AAAA"
    assert r.headers["X-Idempotent-Replay"] == "false"

    g = client.get("/consensus/1")
    assert g.status_code == 200
    assert g.json() == body


def test_idempotent_same_payload_same_version(client):
    payload = _envelope([_rec("r1"), _rec("r2"), _rec("r3")])
    first = client.post("/consensus", json=payload)
    assert first.status_code == 201

    # 载荷顺序不同但内容相同 -> 规范化后仍为同一幂等键
    reordered = _envelope(list(reversed(payload["reads"])))
    second = client.post("/consensus", json=reordered)
    assert second.status_code == 200
    assert second.headers["X-Idempotent-Replay"] == "true"
    assert second.json()["version"] == first.json()["version"]

    assert client.get("/versions").json()["count"] == 1


def test_different_payload_same_request_id_conflicts(client):
    p1 = _envelope([_rec("r1"), _rec("r2"), _rec("r3")])
    client.post("/consensus", json=p1)
    p2 = _envelope([_rec("r1", seq="AAAT"), _rec("r2"), _rec("r3")])
    r = client.post("/consensus", json=p2)
    assert r.status_code == 409
    assert r.json()["error"]["code"] == "payload_conflict"


def test_same_read_id_different_content_conflicts(client):
    p1 = _envelope([_rec("r1"), _rec("r2"), _rec("r3")], request_id="reqA")
    client.post("/consensus", json=p1)
    # 不同 request_id，但 read_id r1 的 sequence 与历史冲突
    p2 = _envelope(
        [_rec("r1", seq="CCCC"), _rec("r9"), _rec("r8")], request_id="reqB"
    )
    r = client.post("/consensus", json=p2)
    assert r.status_code == 409
    assert r.json()["error"]["code"] == "read_conflict"


def test_same_read_id_same_content_distinct_request_ok(client):
    p1 = _envelope([_rec("r1"), _rec("r2"), _rec("r3")], request_id="reqA")
    client.post("/consensus", json=p1)
    p2 = _envelope([_rec("r1"), _rec("r4"), _rec("r5")], request_id="reqB")
    r = client.post("/consensus", json=p2)
    assert r.status_code == 201
    assert r.json()["version"] == 2


def test_versions_are_immutable(client):
    p1 = _envelope(
        [_rec("r1", seq="AAAA"), _rec("r2", seq="AAAA"), _rec("r3", seq="AAAA")],
        request_id="reqA",
    )
    client.post("/consensus", json=p1)
    p2 = _envelope(
        [_rec("r4", seq="CCCC"), _rec("r5", seq="CCCC"), _rec("r6", seq="CCCC")],
        request_id="reqB",
    )
    client.post("/consensus", json=p2)

    v1_first = client.get("/consensus/1").json()
    v1_again = client.get("/consensus/1").json()
    assert v1_first == v1_again
    assert v1_first["consensus"][0]["consensus"] == "AAAA"
    assert client.get("/consensus/2").json()["consensus"][0]["consensus"] == "CCCC"


def test_get_missing_version_404(client):
    r = client.get("/consensus/99")
    assert r.status_code == 404
    assert r.json()["error"]["code"] == "version_not_found"


def test_validation_bad_bases(client):
    rec = _rec("r1", seq="AXAA")
    r = client.post("/consensus", json=_envelope([rec, _rec("r2"), _rec("r3")]))
    assert r.status_code == 422
    codes = [e["code"] for e in r.json()["error"]["details"]["errors"]]
    assert "invalid_value" in codes


def test_validation_quality_range_and_length(client):
    bad_q = _rec("r1", quals=[41, 40, 40, 40])
    r = client.post("/consensus", json=_envelope([bad_q, _rec("r2"), _rec("r3")]))
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "validation_failed"

    bad_len = _rec("r1", quals=[40, 40])
    r2 = client.post("/consensus", json=_envelope([bad_len, _rec("r2"), _rec("r3")]))
    assert r2.status_code == 422
    paths = [e["path"] for e in r2.json()["error"]["details"]["errors"]]
    assert any("qualities" in p for p in paths)


def test_flat_single_record_accepted(client):
    rec = {**_rec("solo"), "request_id": "req1", "sample_id": "s1"}
    r = client.post("/consensus", json=rec)
    # 单条 -> 覆盖不足被淘汰，但仍创建版本
    assert r.status_code == 201
    assert r.json()["eliminated"][0]["reason"] == "coverage_below_3"


def test_bare_array_accepted(client):
    recs = [
        {**_rec("r1"), "request_id": "req1", "sample_id": "s1"},
        {**_rec("r2"), "request_id": "req1", "sample_id": "s1"},
        {**_rec("r3"), "request_id": "req1", "sample_id": "s1"},
    ]
    r = client.post("/consensus", json=recs)
    assert r.status_code == 201
    assert r.json()["version"] == 1


def test_in_batch_duplicate_read_id_conflicting_content_409(client):
    recs = [
        _rec("dup", seq="AAAA"),
        _rec("dup", seq="CCCC"),
        _rec("r3"),
    ]
    r = client.post("/consensus", json=_envelope(recs))
    assert r.status_code == 422
    codes = [e["code"] for e in r.json()["error"]["details"]["errors"]]
    assert "read_conflict" in codes


def test_in_batch_identical_duplicate_is_deduplicated(client):
    recs = [_rec("r1"), _rec("r1"), _rec("r2"), _rec("r3")]
    r = client.post("/consensus", json=_envelope(recs))
    assert r.status_code == 201
    assert r.json()["stats"]["total_reads"] == 3


def test_invalid_json_returns_400(client):
    r = client.post("/consensus", content=b"{not json", headers={"Content-Type": "application/json"})
    assert r.status_code == 400
    assert r.json()["error"]["code"] == "invalid_json"


def test_empty_reads_422(client):
    r = client.post("/consensus", json=_envelope([]))
    assert r.status_code == 422


def test_empty_qualities_length_mismatch_422(client):
    rec = _rec("r1", quals=[])
    r = client.post("/consensus", json=_envelope([rec, _rec("r2"), _rec("r3")]))
    assert r.status_code == 422
    paths = [e["path"] for e in r.json()["error"]["details"]["errors"]]
    assert any(p.endswith("qualities") for p in paths)


def test_ambiguous_cluster_eliminated_end_to_end(client):
    recs = [
        _rec("x1", umi="AAAA", seq="AAAA"),
        _rec("x2", umi="AAAA", seq="AAAA"),
        _rec("x3", umi="AAAA", seq="AAAA"),
        _rec("y1", umi="ACCA", seq="CCCC"),
        _rec("y2", umi="ACCA", seq="CCCC"),
        _rec("y3", umi="ACCA", seq="CCCC"),
        _rec("z1", umi="ACAA", seq="AAAA"),
    ]
    r = client.post("/consensus", json=_envelope(recs))
    assert r.status_code == 201
    body = r.json()
    assert body["consensus"] == []
    assert len(body["eliminated"]) == 2
    assert all(e["reason"] == "ambiguous" for e in body["eliminated"])
    assert body["stats"]["ambiguous_clusters"] == 2


def test_concurrent_distinct_requests_get_distinct_versions(client):
    import concurrent.futures

    def submit(i):
        payload = _envelope(
            [_rec(f"r{i}a"), _rec(f"r{i}b"), _rec(f"r{i}c")],
            request_id=f"req-{i}",
            sample_id="s1",
        )
        r = client.post("/consensus", json=payload)
        return r.status_code, r.json()["version"]

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(submit, range(12)))

    assert all(code == 201 for code, _ in results)
    versions = sorted(v for _, v in results)
    assert versions == list(range(1, 13))
