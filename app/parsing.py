"""POST /consensus 请求解析与校验，并生成幂等用的规范化载荷。"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

from .consensus import Read
from .errors import ValidationError

ALLOWED_BASES = frozenset("ACGT")


@dataclass(frozen=True)
class ParsedSubmission:
    request_id: str
    sample_id: str
    reads: Tuple[Read, ...]
    canonical_payload: dict


def _err(errors: List[dict], path: str, code: str, message: str) -> None:
    errors.append({"path": path, "code": code, "message": message})


def _str_field(
    obj: dict,
    key: str,
    path: str,
    errors: List[dict],
    *,
    nonempty: bool = True,
) -> str:
    if key not in obj:
        _err(errors, f"{path}.{key}", "missing_field", f"缺少必填字段 {key}")
        return ""
    value = obj[key]
    # bool 是 int 的子类，但这里要求 str
    if not isinstance(value, str):
        _err(errors, f"{path}.{key}", "invalid_type", f"字段 {key} 必须是字符串")
        return ""
    if nonempty and value.strip() == "":
        _err(errors, f"{path}.{key}", "empty_field", f"字段 {key} 不能为空")
        return ""
    return value.strip()


def _parse_record(
    raw: Any,
    index: int,
    errors: List[dict],
    top_request_id: str,
    top_sample_id: str,
) -> Tuple[str, str, str, str, str, str, Tuple[int, ...]] | None:
    path = f"reads[{index}]"
    if not isinstance(raw, dict):
        _err(errors, path, "invalid_type", "每条读数必须是对象")
        return None
    error_start = len(errors)

    read_id = _str_field(raw, "read_id", path, errors)
    locus = _str_field(raw, "locus", path, errors)

    umi_raw = _str_field(raw, "umi", path, errors)
    umi = umi_raw.upper()
    if umi and not set(umi) <= ALLOWED_BASES:
        _err(errors, f"{path}.umi", "invalid_value", "umi 只能包含 A/C/G/T")

    seq_raw = _str_field(raw, "sequence", path, errors)
    sequence = seq_raw.upper()
    if sequence and not set(sequence) <= ALLOWED_BASES:
        _err(
            errors,
            f"{path}.sequence",
            "invalid_value",
            "sequence 只能包含 A/C/G/T",
        )

    # qualities：必须是 [0,40] 的整数数组，且与 sequence 等长
    qualities: Tuple[int, ...] = ()
    if "qualities" not in raw:
        _err(errors, f"{path}.qualities", "missing_field", "缺少必填字段 qualities")
    else:
        quals = raw["qualities"]
        if not isinstance(quals, list):
            _err(
                errors,
                f"{path}.qualities",
                "invalid_type",
                "qualities 必须是整数数组",
            )
        else:
            parsed_quals: List[int] = []
            bad = False
            for qi, q in enumerate(quals):
                if isinstance(q, bool) or not isinstance(q, int):
                    _err(
                        errors,
                        f"{path}.qualities[{qi}]",
                        "invalid_type",
                        "质量值必须是整数",
                    )
                    bad = True
                    continue
                if not 0 <= q <= 40:
                    _err(
                        errors,
                        f"{path}.qualities[{qi}]",
                        "invalid_value",
                        "质量值必须在 0~40 之间",
                    )
                    bad = True
                parsed_quals.append(q)
            if not bad:
                qualities = tuple(parsed_quals)
                if sequence and len(qualities) != len(sequence):
                    _err(
                        errors,
                        f"{path}.qualities",
                        "length_mismatch",
                        "qualities 与 sequence 必须等长",
                    )

    # 记录级 request_id / sample_id 若提供必须与顶层一致
    rec_request_id = top_request_id
    rec_sample_id = top_sample_id
    if "request_id" in raw:
        v = raw["request_id"]
        if not isinstance(v, str) or not v.strip():
            _err(
                errors,
                f"{path}.request_id",
                "invalid_type",
                "request_id 必须是非空字符串",
            )
        elif top_request_id and v.strip() != top_request_id:
            _err(
                errors,
                f"{path}.request_id",
                "field_mismatch",
                "记录级 request_id 与顶层不一致",
            )
        else:
            rec_request_id = v.strip()
    if "sample_id" in raw:
        v = raw["sample_id"]
        if not isinstance(v, str) or not v.strip():
            _err(
                errors,
                f"{path}.sample_id",
                "invalid_type",
                "sample_id 必须是非空字符串",
            )
        elif top_sample_id and v.strip() != top_sample_id:
            _err(
                errors,
                f"{path}.sample_id",
                "field_mismatch",
                "记录级 sample_id 与顶层不一致",
            )
        else:
            rec_sample_id = v.strip()

    if len(errors) > error_start:
        return None
    return rec_request_id, rec_sample_id, read_id, locus, umi, sequence, qualities


def parse_submission(body: Any) -> ParsedSubmission:
    """解析三种受支持的请求体：

    - 信封: {"request_id", "sample_id", "reads": [...]}
    - 批量数组: [ {读数, request_id, sample_id}, ... ]
    - 单条扁平读数: {request_id, sample_id, read_id, locus, umi, sequence, qualities}
    """
    errors: List[dict] = []

    if isinstance(body, dict) and "reads" in body:
        request_id = _str_field(body, "request_id", "$", errors)
        sample_id = _str_field(body, "sample_id", "$", errors)
        raw_reads = body["reads"]
        if not isinstance(raw_reads, list):
            _err(errors, "$.reads", "invalid_type", "reads 必须是数组")
            raw_reads = []
        records_raw = list(enumerate(raw_reads))
    elif isinstance(body, list):
        request_id = ""
        sample_id = ""
        records_raw = list(enumerate(body))
    elif isinstance(body, dict):
        request_id = ""
        sample_id = ""
        records_raw = [(0, body)]
    else:
        raise ValidationError(
            [{"path": "$", "code": "invalid_type", "message": "请求体必须是对象或数组"}]
        )

    if not records_raw:
        _err(errors, "$.reads", "empty", "至少需要一条读数")

    parsed: List[Tuple[str, str, str, str, str, str, Tuple[int, ...]]] = []
    seen_ids: Dict[str, Tuple[str, str, str, str, Tuple[int, ...]]] = {}

    for index, raw in records_raw:
        before = len(errors)
        rec = _parse_record(raw, index, errors, request_id, sample_id)
        if rec is None:
            continue
        rec_req, rec_sample, read_id, locus, umi, sequence, qualities = rec

        # 数组形式下，request_id/sample_id 必须全部一致。
        if not request_id:
            if not rec_req:
                _err(
                    errors,
                    f"reads[{index}].request_id",
                    "missing_field",
                    "缺少 request_id",
                )
            elif not request_id:
                request_id = rec_req
            elif rec_req != request_id:
                _err(
                    errors,
                    f"reads[{index}].request_id",
                    "inconsistent",
                    "同一请求中的 request_id 必须一致",
                )
        if not sample_id:
            if not rec_sample:
                _err(
                    errors,
                    f"reads[{index}].sample_id",
                    "missing_field",
                    "缺少 sample_id",
                )
            elif not sample_id:
                sample_id = rec_sample
            elif rec_sample != sample_id:
                _err(
                    errors,
                    f"reads[{index}].sample_id",
                    "inconsistent",
                    "同一请求中的 sample_id 必须一致",
                )

        content = (locus, umi, sequence, qualities)
        if read_id and read_id in seen_ids:
            if seen_ids[read_id] != content:
                _err(
                    errors,
                    f"reads[{index}].read_id",
                    "read_conflict",
                    f"read_id {read_id} 在同一请求中对应了不同的内容",
                )
            # 完全相同的重复读数 -> 静默去重
            continue
        if read_id:
            seen_ids[read_id] = content
        parsed.append(rec)

    if errors:
        raise ValidationError(errors[:50])

    reads = tuple(
        Read(
            read_id=p[2],
            locus=p[3],
            umi=p[4],
            sequence=p[5],
            qualities=p[6],
        )
        for p in parsed
    )

    # 规范化载荷：与提交顺序无关，用于幂等哈希。
    canonical = {
        "request_id": request_id,
        "sample_id": sample_id,
        "reads": sorted(
            (
                {
                    "read_id": r.read_id,
                    "locus": r.locus,
                    "umi": r.umi,
                    "sequence": r.sequence,
                    "qualities": list(r.qualities),
                }
                for r in reads
            ),
            key=lambda x: (x["locus"], x["read_id"]),
        ),
    }
    return ParsedSubmission(
        request_id=request_id,
        sample_id=sample_id,
        reads=reads,
        canonical_payload=canonical,
    )


def canonical_hash(canonical: dict) -> str:
    import hashlib

    blob = json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()
