# UMI 分子共识去噪 API

基于 **Python 3.11 + FastAPI + SQLite** 实现的 UMI（Unique Molecular Identifier）分子共识去噪服务。
它把同一批测序读数按 UMI 分簇、归并误差 UMI，再用 Phred 质量加权投票生成**不可变、版本化**的一致性序列。

- 接口：`POST /consensus`、`GET /consensus/{version}`（另附 `GET /versions`、`GET /health`）
- 存储：SQLite（WAL 模式），写入串行化，读取始终命中已提交快照
- 一键启动：Docker / Docker Compose

---

## 1. 算法规则

### 1.1 分簇（在每个 `locus` 内独立进行）

1. 按 UMI **出现频次降序**排序，频次相同按 UMI **字典序升序**处理；
2. 对每个 UMI，与所有已建簇的代表比较汉明距离：
   - 命中 **0** 个距离 ≤ 1 的代表 → 自成一簇；
   - 命中 **1** 个 → 归入该簇（代表保持不变）；
   - 命中 **≥2** 个 → 这些被命中的簇**整体标记为 `ambiguous`**（该 UMI 是桥接分子），全部淘汰；
   - 长度不同的 UMI 汉明距离视为无穷大，不会归并。

### 1.2 簇内共识（逐位加权投票）

- 对每个位点，碱基 `b ∈ {A,C,G,T}` 的权重为该位所有支持读数的 **Phred 权重之和**：
  `weight(b) = Σ 10^(Q/10)`；
- 该位点被判定为最高权重碱基，当且仅当：
  - **覆盖读数 ≥ 3**，且
  - **最高权重 ≥ 其余碱基权重总和的 2 倍**（`best ≥ 2 × rest`）；
- 任一条件不满足，该位记为 **`N`**（并列权重天然无法满足 2 倍条件，确定性地记为 N）。

### 1.3 淘汰规则

簇被淘汰的原因（`eliminated[].reason`，按优先级）：

| reason | 含义 |
| --- | --- |
| `ambiguous` | 桥接 UMI 同时命中多个代表，整簇 ambiguous（不进行投票） |
| `coverage_below_3` | 簇内读数不足 3 条，无法保证覆盖 |
| `n_in_consensus` | 共识序列含 `N`（权重不占优或个别位置覆盖不足） |

最终 `consensus` 只保留通过全部检查的簇，并按 **`locus` 字典序**排序。

---

## 2. 请求与幂等 / 冲突语义

### 请求体（三种等价写法）

推荐使用**信封**格式：

```json
{
  "request_id": "demo-001",
  "sample_id": "sample-A",
  "reads": [
    {
      "read_id": "l1-01",
      "locus": "L1",
      "umi": "AAAA",
      "sequence": "AAAA",
      "qualities": [40, 40, 40, 40]
    }
  ]
}
```

也支持：

- **裸数组**：`[ { ...read, request_id, sample_id}, ... ]`（每条带 `request_id/sample_id` 且须一致）；
- **单条扁平读数**：`{ request_id, sample_id, read_id, locus, umi, sequence, qualities }`。

字段校验：

- `request_id`、`sample_id`、`read_id`、`locus`：非空字符串；
- `umi`：非空，仅含 `A/C/G/T`（小写会自动转大写）；
- `sequence`：仅含 `A/C/G/T`（小写自动转大写）；
- `qualities`：**整数数组**，取值 `0~40`，且与 `sequence` **等长**；
- 校验失败返回 `422`，错误体带逐项 `path/code/message`。

### 幂等与冲突

幂等键为 `(request_id, 规范化载荷 SHA-256)`。规范化会对读数按
`(locus, read_id)` 排序，因此**改变读数提交顺序不影响幂等性**。

| 情形 | 结果 |
| --- | --- |
| 同 `request_id` + 同载荷重复提交 | `200` 幂等重放，返回**同一个 version**，响应头 `X-Idempotent-Replay: true` |
| 同 `request_id` + **不同载荷** | `409 payload_conflict` |
| **同 `read_id` 但内容不同**（locus/umi/sequence/qualities 任一不同，跨 request 也算） | `409 read_conflict` |
| 同 `read_id` 同内容（即使在不同 request 中复用） | 正常，不冲突 |
| 全新提交 | `201`，`X-Idempotent-Replay: false` |
| 同一请求内完全相同的重复读数 | 静默去重；同 `read_id` 不同内容则 `422` |

### 版本不可变、读取一致

- 每次**新**提交生成一个全局单调递增的 `version`，结果整体存入 SQLite；
- `GET /consensus/{version}` 返回该版本的**不可变快照**（含 `version/content_hash/created_at`），之后的提交不会改变它；
- SQLite 运行在 WAL 模式，写事务 `BEGIN IMMEDIATE` 并由进程锁串行化，读请求永远读到已提交的一致版本。

---

## 3. Docker Compose 一键启动

前置：已安装 Docker 与 Docker Compose。

```bash
docker compose up --build
```

启动后：

- API 地址：<http://localhost:8000>
- 健康检查：<http://localhost:8000/health>
- Swagger 交互文档：<http://localhost:8000/docs>
- SQLite 数据持久化在命名卷 `umi-data`（容器内 `/data/umi.db`）

后台运行 / 停止：

```bash
docker compose up --build -d
docker compose down            # 停止；数据保留在卷中
docker compose down -v         # 同时删除数据卷
```

---

## 4. 实际请求 / 响应示例

仓库内置了一个覆盖全部淘汰原因的样例数据集（19 条读数，3 个 locus）：

- `L1`：1 个干净簇，距离 1 的 UMI 变体并入，共识 `AAAA`（保留）；
- `L2`：`AAAA` 与 `ACCA` 距离 2，`ACAA` 同时距离二者为 1 → 两簇均 `ambiguous`；
- `L3`：一个 2v2 等权并列产生 `N` 的簇（`n_in_consensus`），一个仅 2 条读数的簇（`coverage_below_3`）。

完整请求见 [`examples/request.json`](examples/request.json)，完整响应见
[`examples/response.json`](examples/response.json)（下方为同一真实运行结果的节选）。

### 4.1 提交

```bash
curl -i -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' \
  --data @examples/request.json
```

首次返回 **`201 Created`**，响应头：

```
x-idempotent-replay: false
```

响应体（节选）：

```json
{
  "request_id": "demo-001",
  "sample_id": "sample-A",
  "consensus": [
    {
      "cluster_id": "L1:c1",
      "locus": "L1",
      "cluster_index": 1,
      "representative": "AAAA",
      "read_count": 6,
      "read_ids": ["l1-01", "l1-02", "l1-03", "l1-04", "l1-05", "l1-06"],
      "umis": [
        { "umi": "AAAA", "read_count": 4, "read_ids": ["l1-01", "l1-02", "l1-03", "l1-04"] },
        { "umi": "AAAT", "read_count": 2, "read_ids": ["l1-05", "l1-06"] }
      ],
      "ambiguous": false,
      "consensus": "AAAA",
      "length": 4,
      "mean_coverage": 6.0
    }
  ],
  "clusters": [
    { "cluster_id": "L1:c1", "locus": "L1", "representative": "AAAA", "read_count": 6, "ambiguous": false, "status": "kept",     "reason": null,              "consensus": "AAAA", "bridges": [] },
    { "cluster_id": "L2:c1", "locus": "L2", "representative": "AAAA", "read_count": 4, "ambiguous": true,  "status": "eliminated", "reason": "ambiguous",       "consensus": null,   "bridges": ["ACAA"] },
    { "cluster_id": "L2:c2", "locus": "L2", "representative": "ACCA", "read_count": 3, "ambiguous": true,  "status": "eliminated", "reason": "ambiguous",       "consensus": null,   "bridges": ["ACAA"] },
    { "cluster_id": "L3:c1", "locus": "L3", "representative": "TTTT", "read_count": 4, "ambiguous": false, "status": "eliminated", "reason": "n_in_consensus",  "consensus": "NTTT", "bridges": [] },
    { "cluster_id": "L3:c2", "locus": "L3", "representative": "GGGG", "read_count": 2, "ambiguous": false, "status": "eliminated", "reason": "coverage_below_3","consensus": "NNNN", "bridges": [] }
  ],
  "eliminated": [ "…（含每个淘汰簇的逐位明细，见 examples/response.json）…" ],
  "stats": {
    "total_reads": 19,
    "total_umis": 7,
    "total_clusters": 5,
    "ambiguous_clusters": 2,
    "eliminated_clusters": 4,
    "kept_clusters": 1,
    "eliminated_by_reason": {
      "ambiguous": 2,
      "coverage_below_3": 1,
      "n_in_consensus": 1
    },
    "n_positions": 5
  },
  "version": 1,
  "content_hash": "b1ec239ad171ac5b76eb308aa7147fc00af8fad088b7f48459774142d325248b",
  "created_at": "2026-09-30T23:26:19.542578Z"
}
```

`eliminated` 中还包含每个 N 位点的 `coverage / weights / reasons`。例如 `L3:c1` 第 1 位
A、T 各 2 票且等权（各 `20000.0`），不满足 2 倍条件：

```json
{
  "position": 1,
  "call": "N",
  "winner": null,
  "coverage": 4,
  "weights": { "A": 20000.0, "C": 0.0, "G": 0.0, "T": 20000.0 },
  "reasons": ["weight_not_dominant"]
}
```

### 4.2 幂等重放

原样（或仅打乱 `reads` 顺序）再次提交：

```bash
curl -i -X POST http://localhost:8000/consensus \
  -H 'Content-Type: application/json' \
  --data @examples/request.json
```

返回 **`200 OK`**，仍是 **version 1**，响应头 `x-idempotent-replay: true`。

### 4.3 读取不可变版本

```bash
curl http://localhost:8000/consensus/1
```

返回与首次提交完全一致的快照（含相同的 `content_hash`）。不存在的版本返回 `404`：

```json
{ "error": { "code": "version_not_found", "message": "版本 99 不存在", "details": { "version": 99 } } }
```

### 4.4 冲突与校验示例

同一 `request_id` 但改一个碱基再提交 → **`409`**：

```json
{ "error": { "code": "payload_conflict",
  "message": "request_id demo-001 已存在但请求载荷不同",
  "details": { "request_id": "demo-001", "existing_version": 1 } } }
```

非法 UMI / 越界质量值 → **`422`**：

```bash
curl -X POST http://localhost:8000/consensus -H 'Content-Type: application/json' \
  -d '{"request_id":"bad","sample_id":"s","reads":[{"read_id":"r1","locus":"L","umi":"AX","sequence":"AC","qualities":[41,9]}]}'
```

```json
{ "error": { "code": "validation_failed", "message": "请求载荷校验失败", "details": { "errors": [
  { "path": "reads[0].umi", "code": "invalid_value", "message": "umi 只能包含 A/C/G/T" },
  { "path": "reads[0].qualities[0]", "code": "invalid_value", "message": "质量值必须在 0~40 之间" }
] } } }
```

### 4.5 列出版本

```bash
curl http://localhost:8000/versions
curl "http://localhost:8000/versions?request_id=demo-001"
```

---

## 5. 本地开发（不使用 Docker）

```bash
python3 -m venv --without-pip .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/uvicorn app.main:app --reload --port 8000
```

数据库路径可用环境变量 `UMI_DB_PATH` 覆盖（默认 `./data/umi.db`）。

### 运行测试

```bash
.venv/bin/python -m pip install pytest httpx
.venv/bin/python -m pytest -q
```

测试覆盖：汉明距离与归并、频次/字典序、多代表 ambiguous、Phred 权重 2 倍边界、
覆盖 <3、N 淘汰、locus 排序，以及 API 层的幂等、409（载荷/read_id）、422 校验、
404、版本不可变与读取一致。

## 6. 项目结构

```
app/
  consensus.py   # 分簇 + Phred 权重投票纯函数核心
  parsing.py     # 请求解析/校验 + 规范化幂等哈希
  db.py          # SQLite schema、串行写事务、读连接
  service.py     # 幂等/冲突/版本编排
  errors.py      # 统一错误类型
  main.py        # FastAPI 路由
examples/        # 真实请求与响应样例
tests/           # pytest 测试
Dockerfile
docker-compose.yml
```
