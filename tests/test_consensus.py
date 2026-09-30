"""算法核心测试：分簇、权重投票、淘汰原因、排序。"""

from app.consensus import (
    Read,
    build_consensus_result,
    hamming_distance,
    weighted_consensus,
)


def R(read_id, locus, umi, sequence, qualities):
    return Read(read_id, locus, umi, sequence, tuple(qualities))


def test_hamming():
    assert hamming_distance("ACGT", "ACGT") == 0
    assert hamming_distance("ACGT", "ACGA") == 1
    assert hamming_distance("ACGT", "TGCA") == 4
    assert hamming_distance("ACG", "ACGT") > 1


def test_clean_cluster_consensus_and_stats():
    reads = [
        R("r1", "L1", "ACGT", "AAAA", [40, 40, 40, 40]),
        R("r2", "L1", "ACGT", "AAAA", [30, 30, 30, 30]),
        R("r3", "L1", "ACGT", "AAAA", [20, 20, 20, 20]),
    ]
    res = build_consensus_result("req", "s1", reads)
    assert res["consensus"][0]["consensus"] == "AAAA"
    assert res["stats"]["kept_clusters"] == 1
    assert res["stats"]["eliminated_clusters"] == 0


def test_weight_threshold_boundary():
    # 3 票：2 票 A（10000 each），1 票 C（10000）-> 20000 >= 2*10000 -> A
    reads = [
        R("r1", "L", "U", "A", [40]),
        R("r2", "L", "U", "A", [40]),
        R("r3", "L", "U", "C", [40]),
    ]
    seq, positions = weighted_consensus(reads)
    assert seq == "A"

    # 2 票 A（10000），2 票 C（10000）-> 10000*2 < 2*20000 -> N（并列）
    reads2 = [
        R("r1", "L", "U", "A", [40]),
        R("r2", "L", "U", "A", [40]),
        R("r3", "L", "U", "C", [40]),
        R("r4", "L", "U", "C", [40]),
    ]
    seq2, positions2 = weighted_consensus(reads2)
    assert seq2 == "N"
    assert "weight_not_dominant" in positions2[0]["reasons"]


def test_coverage_below_three_yields_n():
    reads = [
        R("r1", "L", "U", "A", [40]),
        R("r2", "L", "U", "A", [40]),
    ]
    seq, positions = weighted_consensus(reads)
    assert seq == "N"
    assert positions[0]["coverage"] == 2
    assert "coverage_below_3" in positions[0]["reasons"]


def test_cluster_ordering_frequency_then_lexicographic():
    reads = []
    # BBB: 3 次，AAA: 2 次 -> BBB 频次更高
    for i in range(3):
        reads.append(R(f"b{i}", "L", "BBB", "AAAA", [40] * 4))
    for i in range(2):
        reads.append(R(f"a{i}", "L", "AAA", "CCCC", [40] * 4))
    res = build_consensus_result("req", "s", reads)
    # AAA vs BBB 汉明距离 3，各自成簇
    assert res["clusters"][0]["representative"] == "BBB"
    assert res["clusters"][1]["representative"] == "AAA"


def test_hamming_one_merges_into_representative():
    reads = [
        R("r1", "L", "AAAA", "AAAA", [40] * 4),
        R("r2", "L", "AAAA", "AAAA", [40] * 4),
        R("r3", "L", "AAAA", "AAAA", [40] * 4),
        # AAAT 与 AAAA 距离 1，归入同簇（单代表命中，非 ambiguous）
        R("r4", "L", "AAAT", "AAAA", [40] * 4),
    ]
    res = build_consensus_result("req", "s", reads)
    assert len(res["clusters"]) == 1
    assert res["clusters"][0]["ambiguous"] is False
    assert res["clusters"][0]["read_count"] == 4
    assert {u["umi"] for u in res["clusters"][0]["umis"]} == {"AAAA", "AAAT"}


def test_multiple_representative_hits_marks_all_ambiguous():
    # AAAA freq 3，ACCA freq 3（距离 2），ACAA freq 1 与两者距离都为 1
    reads = []
    for i in range(3):
        reads.append(R(f"x{i}", "L", "AAAA", "AAAA", [40] * 4))
    for i in range(3):
        reads.append(R(f"y{i}", "L", "ACCA", "CCCC", [40] * 4))
    reads.append(R("z0", "L", "ACAA", "AAAA", [40] * 4))
    res = build_consensus_result("req", "s", reads)
    assert len(res["clusters"]) == 2
    assert all(c["ambiguous"] for c in res["clusters"])
    assert res["stats"]["ambiguous_clusters"] == 2
    assert all(e["reason"] == "ambiguous" for e in res["eliminated"])
    # 桥接 UMI 分配给优先级更高（字典序 AAAA）的簇
    assert res["clusters"][0]["representative"] == "AAAA"
    assert res["clusters"][0]["read_count"] == 4
    assert res["clusters"][1]["read_count"] == 3
    bridges = {b for c in res["clusters"] for b in c["bridges"]}
    assert bridges == {"ACAA"}


def test_n_consensus_eliminated():
    # 3 票各为 A/C/G -> 无人占优 -> N -> 整簇淘汰
    reads = [
        R("r1", "L", "U", "A", [40]),
        R("r2", "L", "U", "C", [40]),
        R("r3", "L", "U", "G", [40]),
    ]
    res = build_consensus_result("req", "s", reads)
    assert res["consensus"] == []
    assert res["eliminated"][0]["reason"] == "n_in_consensus"
    assert res["stats"]["eliminated_by_reason"]["n_in_consensus"] == 1
    assert res["stats"]["n_positions"] == 1


def test_locus_sorting_in_output():
    reads = [
        R("r1", "L2", "U", "AAAA", [40] * 4),
        R("r2", "L2", "U", "AAAA", [40] * 4),
        R("r3", "L2", "U", "AAAA", [40] * 4),
        R("r4", "L1", "U", "CCCC", [40] * 4),
        R("r5", "L1", "U", "CCCC", [40] * 4),
        R("r6", "L1", "U", "CCCC", [40] * 4),
    ]
    res = build_consensus_result("req", "s", reads)
    loci = [c["locus"] for c in res["consensus"]]
    assert loci == ["L1", "L2"]
