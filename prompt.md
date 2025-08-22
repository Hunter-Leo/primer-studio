## 目的说明

本项目旨在开发一个**批量引物设计工具**，实现以下目标：

1. 从用户提供的 FASTA 文件中读取多个 DNA 序列。
2. 自动为每条序列设计正向和反向引物。
3. 输出每条引物的序列、Tm、GC 含量、长度等信息。
4. 支持命令行调用，也可作为 Python 包导入使用。
5. 支持配置经验参数，方便不熟悉 PCR 引物设计的用户使用。
6. 高效运行，尽量减少内存消耗和计算资源占用。
7. 采用面向对象设计和结构化代码，便于扩展和维护。

---

## 需求说明

* **输入**

  * 支持 FASTA 文件，含一条或多条 DNA 序列。
  * 可通过命令行参数或 Python 接口传入。

* **输出**

  * 每条序列对应的引物信息表，包含：

    * 正向引物序列
    * 反向引物序列
    * Tm
    * GC 含量
    * 长度
  * 输出格式可为 CSV 或 JSON，可配置。

* **功能**

  1. 批量处理多条序列
  2. 提供经验参数，用户可快速使用
  3. 使用 log 模块记录运行信息，而非 print
  4. 支持面向对象设计，易于扩展
  5. 高性能和低资源消耗
  6. 使用 `pydantic v2` 进行数据和配置校验
  7. 命令行入口（CLI）设计，支持 uv 项目管理启动

* **核心依赖**

  * `primer3-py`
  * `pydantic v2`
  * `Biopython`（读取 FASTA）
  * `logging` 日志
  * `uv`（项目管理工具）

* **编码要求**
  * 代码应进行注释，包括函数、类、模块等
  * 参数应尽量做类型注解和注释
  * 注释使用 Google 风格，全英文注释
---

## 代码层级和结构设计（src 布局）

### 项目目录示例

```
primer_designer_project/
│
├─ src/
│   └─ primer_designer/
│       ├─ __init__.py
│       ├─ config.py          # PrimerDesignerConfig 定义及常用实例
│       ├─ models.py          # Pydantic 数据模型（引物信息、序列信息）
│       ├─ primer.py          # PrimerDesigner 类封装 primer3-py 调用
│       ├─ fasta_reader.py    # FASTA 文件读取工具
│       ├─ batch_designer.py  # 批量处理逻辑
│       ├─ utils.py           # 工具函数（GC计算、Tm计算、日志初始化）
│       └─ cli.py             # 命令行入口（使用 Typer 或 argparse）
│
├─ tests/
│   ├─ test_primer.py
│   ├─ test_batch.py
│   └─ test_fasta.py
│
├─ pyproject.toml             # 使用 uv 进行项目管理
└─ README.md
```

---

### 类和模块设计

#### 1. `PrimerDesignerConfig`（config.py）

* 基于 `pydantic.BaseModel`
* 定义引物设计参数，如：

  * `PRIMER_OPT_SIZE`, `PRIMER_MIN_SIZE`, `PRIMER_MAX_SIZE`
  * `PRIMER_OPT_TM`, `PRIMER_MIN_TM`, `PRIMER_MAX_TM`
  * `PRIMER_MAX_SELF_ANY`, `PRIMER_MAX_POLY_X`
* 提供几个常用经验参数实例，如：

  * `STANDARD_CONFIG`
  * `GC_RICH_CONFIG`
  * `GC_POOR_CONFIG`

#### 2. `Primer`（models.py）

* Pydantic 模型
* 属性：

  * `sequence_id`
  * `forward_primer`
  * `reverse_primer`
  * `tm_forward`
  * `tm_reverse`
  * `gc_forward`
  * `gc_reverse`
  * `length_forward`
  * `length_reverse`

#### 3. `PrimerDesigner`（primer.py）

* 面向对象类，封装 primer3-py 调用
* **通过 `PrimerDesignerConfig` 初始化**
* 方法：

  * `design_primer(sequence: str) -> Primer`
  * 内部使用初始化时传入的 `PrimerDesignerConfig` 对象

#### 4. `FastaReader`（fasta\_reader.py）

* 读取 FASTA 文件
* 返回序列列表，或者封装为 Pydantic 模型

#### 5. `BatchDesigner`（batch\_designer.py）

* 接收多条序列，批量生成引物
* 支持并行处理（可选，使用 multiprocessing 或线程池）
* 输出 CSV/JSON

#### 6. `CLI`（cli.py）

* 命令行入口
* 功能：

  * 输入 FASTA 文件路径
  * 输出文件路径和格式
  * 可覆盖默认经验参数（选择不同 `PrimerDesignerConfig`）

#### 7. `utils.py`

* GC 含量计算
* Tm 计算（可用 primer3 自带函数）
* 日志初始化

---

### 数据流示意

```
FASTA 文件 → FastaReader → BatchDesigner → PrimerDesigner(PrimerDesignerConfig) → Primer 模型 → 输出 CSV/JSON
```

---

### 性能优化建议

1. 批量处理时，可使用生成器减少内存消耗。
2. 并行化设计引物（序列数量大时可用线程池或进程池）。
3. 避免重复计算 Tm 和 GC，可直接用 primer3 返回的结果。
4. 日志按级别输出，避免打印大量中间信息。

---

这个 `prompt.md` 文件可以直接用于 AI 编程工具生成：

* 支持 CLI 和 Python 包使用
* 面向对象
* 使用 `primer3-py`
* `PrimerDesigner` 通过 `PrimerDesignerConfig` 初始化
* 包含常用经验参数实例
* 基于 `src/primer_designer` 布局
