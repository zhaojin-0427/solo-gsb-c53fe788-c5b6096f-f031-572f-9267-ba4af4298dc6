"""UMI 分簇与 Phred 权重共识的纯函数核心。

流程：
1. 在每个 locus 内，按 UMI 频次降序、频次相同按 UMI 字典序处理；
2. 每个 UMI 与已有簇代表比较汉明距离，距离 <= 1 即归入该簇；
   - 命中 0 个代表：自成一簇；
   - 命中 1 个代表：归入该簇；
   - 命中 >= 2 个代表：所有命中簇整体标记为 ambiguous（桥接）；
3. 簇内逐位以 Phred 权重 10^(Q/10) 对 A/C/G/T 投票，
   最高权重必须 >= 其余总和的 2 倍，且该位覆盖读数 >= 3，
   否则该位记为 N；
4. 含 N 的共识（以及 ambiguous、覆盖不足的簇）被淘汰。
"""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

BASES = ("A", "C", "G", "T")


@dataclass(frozen=True)
class Read:
    """一条去重后的读数。"""

    read_id: str
    locus: str
    umi: str
    sequence: str
    qualities: Tuple[int, ...]


@dataclass
class _Cluster:
    """分簇过程中的中间结构。"""

    index: int
    representative: str
    # umi -> 归入该 UMI 的读数（构建完成后再排序输出）
    members: Dict[str, List[Read]] = field(default_factory=dict)
    ambiguous: bool = False
    bridges: List[dict] = field(default_factory=list)

    @property
    def all_reads(self) -> List[Read]:
        out: List[Read] = []
        for reads in self.members.values():
            out.extend(reads)
        return out


def hamming_distance(a: str, b: str) -> int:
    """等长字符串的汉明距离；不等长视为不可能归入同簇。"""
    if len(a) != len(b):
        return 10**9
    return sum(1 for x, y in zip(a, b) if x != y)


def cluster_umis(reads: Sequence[Read]) -> List[_Cluster]:
    """按 UMI 频次降序、字典序分簇。"""
    by_umi: Dict[str, List[Read]] = defaultdict(list)
    for r in reads:
        by_umi[r.umi].append(r)

    ordered_umis = sorted(by_umi.keys(), key=lambda u: (-len(by_umi[u]), u))

    clusters: List[_Cluster] = []
    for umi in ordered_umis:
        hits = [
            c
            for c in clusters
            if hamming_distance(umi, c.representative) <= 1
        ]
        if not hits:
            cluster = _Cluster(index=len(clusters) + 1, representative=umi)
            cluster.members[umi] = list(by_umi[umi])
            clusters.append(cluster)
            continue

        chosen = hits[0]  # hits 按代表处理顺序（即频次/字典序优先级）排列
        chosen.members.setdefault(umi, []).extend(by_umi[umi])

        if len(hits) >= 2:
            bridge = {
                "umi": umi,
                "representatives": [h.representative for h in hits],
                "cluster_indices": [h.index for h in hits],
                "assigned_to": chosen.index,
            }
            for h in hits:
                h.ambiguous = True
                h.bridges.append(bridge)
    return clusters


def weighted_consensus(
    reads: Sequence[Read],
) -> Tuple[str, List[dict]]:
    """簇内逐位 Phred 权重投票。

    返回 (共识序列, 逐位明细)。长度取簇内最长序列；长度不足的读数
    不参与超过其长度的位置（覆盖按实际参与该位的读数计）。
    """
    length = max((len(r.sequence) for r in reads), default=0)
    consensus: List[str] = []
    positions: List[dict] = []

    for pos in range(length):
        weights: Dict[str, float] = {b: 0.0 for b in BASES}
        participants = 0
        for r in reads:
            if pos >= len(r.sequence):
                continue
            weights[r.sequence[pos]] += 10.0 ** (r.qualities[pos] / 10.0)
            participants += 1

        # 确定性打破并列：取权重最高者；权重相同按 A,C,G,T 字典序。
        winner = max(BASES, key=lambda b: (weights[b], -BASES.index(b)))
        best = weights[winner]
        rest = sum(weights[b] for b in BASES if b != winner)
        enough_coverage = participants >= 3
        dominant = best >= 2.0 * rest
        call = winner if enough_coverage and dominant else "N"
        consensus.append(call)

        reasons: List[str] = []
        if not enough_coverage:
            reasons.append("coverage_below_3")
        if not dominant:
            reasons.append("weight_not_dominant")

        positions.append(
            {
                "position": pos + 1,
                "call": call,
                "winner": winner if dominant else None,
                "coverage": participants,
                "weights": {b: round(weights[b], 4) for b in BASES},
                "reasons": reasons,
            }
        )

    return "".join(consensus), positions


def _umi_rows(cluster: _Cluster) -> List[dict]:
    return [
        {
            "umi": umi,
            "read_count": len(cluster.members[umi]),
            "read_ids": sorted(r.read_id for r in cluster.members[umi]),
        }
        for umi in sorted(cluster.members.keys())
    ]


def _reads_by_locus(reads: Sequence[Read]) -> OrderedDict:
    out: OrderedDict = OrderedDict()
    for r in reads:
        out.setdefault(r.locus, []).append(r)
    return out


def build_consensus_result(
    request_id: str,
    sample_id: str,
    reads: Sequence[Read],
) -> dict:
    """对一次提交的全部读数构建完整的共识结果。"""
    by_locus = _reads_by_locus(reads)

    consensus_rows: List[dict] = []
    mapping_rows: List[dict] = []
    eliminated_rows: List[dict] = []
    stats: Dict[str, object] = {
        "total_reads": len(reads),
        "total_umis": 0,
        "total_clusters": 0,
        "ambiguous_clusters": 0,
        "eliminated_clusters": 0,
        "kept_clusters": 0,
        "eliminated_by_reason": {
            "ambiguous": 0,
            "coverage_below_3": 0,
            "n_in_consensus": 0,
        },
        "n_positions": 0,
    }

    for locus in sorted(by_locus.keys()):
        locus_reads = by_locus[locus]
        clusters = cluster_umis(locus_reads)
        stats["total_umis"] += len({r.umi for r in locus_reads})
        stats["total_clusters"] += len(clusters)

        for cluster in clusters:
            cluster_id = f"{locus}:c{cluster.index}"
            member_reads = cluster.all_reads
            read_ids = sorted(r.read_id for r in member_reads)
            row_base = {
                "cluster_id": cluster_id,
                "locus": locus,
                "cluster_index": cluster.index,
                "representative": cluster.representative,
                "read_count": len(member_reads),
                "read_ids": read_ids,
                "umis": _umi_rows(cluster),
                "ambiguous": cluster.ambiguous,
            }

            # 1) 多代表命中 -> 整簇 ambiguous，淘汰且不再投票。
            if cluster.ambiguous:
                reason = "ambiguous"
                stats["ambiguous_clusters"] += 1
                stats["eliminated_clusters"] += 1
                stats["eliminated_by_reason"][reason] += 1
                eliminated_rows.append(
                    {
                        **row_base,
                        "reason": reason,
                        "bridges": cluster.bridges,
                        "consensus": None,
                        "positions": [],
                        "n_positions": [],
                    }
                )
                mapping_rows.append(
                    {
                        **row_base,
                        "status": "eliminated",
                        "reason": reason,
                        "consensus": None,
                        "bridges": [b["umi"] for b in cluster.bridges],
                    }
                )
                continue

            sequence, positions = weighted_consensus(member_reads)
            n_positions = [
                {
                    "position": p["position"],
                    "coverage": p["coverage"],
                    "reasons": p["reasons"],
                    "weights": p["weights"],
                }
                for p in positions
                if p["call"] == "N"
            ]
            stats["n_positions"] += len(n_positions)

            # 2) 覆盖不足（读数 < 3）优先归因；3) 其次为含 N。
            if len(member_reads) < 3:
                reason = "coverage_below_3"
            elif "N" in sequence:
                reason = "n_in_consensus"
            else:
                reason = None

            if reason is not None:
                stats["eliminated_clusters"] += 1
                stats["eliminated_by_reason"][reason] += 1
                eliminated_rows.append(
                    {
                        **row_base,
                        "reason": reason,
                        "bridges": [],
                        "consensus": sequence,
                        "positions": positions,
                        "n_positions": n_positions,
                    }
                )
                mapping_rows.append(
                    {
                        **row_base,
                        "status": "eliminated",
                        "reason": reason,
                        "consensus": sequence,
                        "bridges": [],
                    }
                )
                continue

            stats["kept_clusters"] += 1
            mean_cov = round(
                sum(p["coverage"] for p in positions) / max(len(positions), 1),
                2,
            )
            consensus_rows.append(
                {
                    **row_base,
                    "consensus": sequence,
                    "length": len(sequence),
                    "mean_coverage": mean_cov,
                }
            )
            mapping_rows.append(
                {
                    **row_base,
                    "status": "kept",
                    "reason": None,
                    "consensus": sequence,
                    "bridges": [],
                }
            )

    return {
        "request_id": request_id,
        "sample_id": sample_id,
        "consensus": consensus_rows,
        "clusters": mapping_rows,
        "eliminated": eliminated_rows,
        "stats": stats,
    }
