
<center>
 <img src="resources\spacerscope_pixel_icon.png" style="zoom: 5%;" />
</center>
# **<center>SpacerScope 用户手册（中文版)</center>**





## 1. 简介

### 什么是 SpacerScope？

**SpacerScope** 是一款面向 **CRISPR/Cas 基因组编辑**应用的高效脱靶分析工具，旨在解决大规模参考基因组，尤其是复杂植物基因组场景下传统方法计算开销高、响应慢、难以兼顾速度与精度的问题。

SpacerScope 的核心思想是将脱靶搜索拆分为**位点提取、快速预筛和精确验证**三个阶段。程序首先从查询序列和参考基因组中提取满足 PAM 约束的候选 spacer 位点；随后针对不同搜索任务采用不同的高效编码策略：在**错配搜索**中，使用紧凑的 **2-bit 序列编码**与批量化扫描策略进行快速汉明距离判断；在**插入/缺失搜索**中，先将序列编码为多通道二进制特征，利用**通道过滤与批量掩码融合**大幅压缩候选集合，再对保留下来的少量候选执行精确比对。对于 indel 验证，SpacerScope 采用的是**右端锚定的编辑距离动态规划**，而不是简单的局部比对，从而更适合 CRISPR spacer 场景下插入、缺失与末端约束共存的序列比较需求。

因此，SpacerScope 并不是单纯依赖传统全序列比对，而是通过**二进制编码驱动的高效预过滤**与**精确编辑距离验证**相结合，在保证结果完整性和准确性的前提下，显著降低整体搜索成本。实际测试表明，SpacerScope 能够在复杂基因组背景下实现快速、稳定且高精度的 gRNA 筛选与脱靶分析。

###  SpacerScope的核心特点

1. **完美的脱靶检出能力**
   在给定 PAM、错配数和 indel 数等约束条件下，SpacerScope 能够对满足条件的候选位点进行系统性搜索，并在精确验证阶段尽可能保证结果完整性与判定准确性。
   
2. **同时支持错配与 indel 脱靶分析**
   SpacerScope 不仅支持常规的**碱基替换（mismatch）**检测，还支持**插入（insertion）**与**缺失（deletion）**的脱靶识别，能够覆盖更常见也更复杂的 CRISPR/Cas 潜在脱靶类型。
3. **高效的二进制编码与预过滤机制**
   程序针对不同任务采用专门设计的二进制表示与预筛策略：
   - 在错配搜索中，使用紧凑的 **2-bit 编码**和批量扫描加速大规模汉明距离计算；
   - 在 indel 搜索中，使用多通道编码、二进制过滤和批量掩码融合显著降低需要进入精确比对的候选规模。
     这种“先快速过滤、后精确验证”的设计显著提升了整体运行效率。
4. **适用于大规模基因组分析**
   通过针对内存访问、并行调度和候选筛选流程的优化，SpacerScope 能够在较大的参考基因组上保持较好的计算效率，尤其适合复杂植物基因组等高负载场景。
5. **多模式输入与分析流程**
   SpacerScope 同时支持：
   - 直接输入 FASTA 查询序列；
   - 基于注释文件提取目标序列进行分析；
   - 常规模式与低内存 cut 模式。
     其中 cut 模式可在内存资源有限的设备上按染色体流式处理参考基因组。
6. **跨平台兼容**
   提供 **Windows** 与 **Linux** 平台支持，既可作为命令行工具运行，也可通过 spacerscope-web 以浏览器方式使用。
7. **并行与底层优化支持**
   程序支持多线程并行执行，并可在支持的硬件与编译环境下利用 **AVX2 / SSE4.1** 等指令优化，以进一步提升实际运行性能。
8. **适配多种 CRISPR 系统**
   支持用户自定义 **PAM 序列**及其相对方向，可适配不同核酸酶系统和不同实验设计需求。

---

## 2. 程序组成

当前版本发布时通常包含两个可执行文件：

- **`spacerscope`**：主程序
- **`spacerscope-web`**：Web 图形壳，内部通过子进程调用 `spacerscope`

两者关系如下：

- `spacerscope` 负责所有实际计算，可单独运行。
- `spacerscope-web` 通过浏览器实现GUI，也方便网络访问。spacerscope-web 不是独立搜索引擎，必须与 `spacerscope` 搭配使用

---

## 3. 环境要求

### 支持的操作系统

- Windows x64
- Linux x86_64

### CPU 指令集

推荐 CPU 支持：

- AVX2
- SSE4.1

Linux 可参考：

```bash
lscpu | grep -E 'avx2|sse4_1' | sed 's/^[[:space:]]*//'
```

程序即使在不支持 AVX2 的机器上也可运行，但性能可能下降。

### 运行依赖

核心依赖为：

- OpenMP
- zlib

如使用 `spacerscope-web`，无需额外安装浏览器以外的软件；只需保证同目录或指定路径下能找到 `spacerscope` 主程序。

---

## 4. 运行模式

SpacerScope 当前支持 4 种主模式：

| 模式 | 输入方式 | 内存策略 | 说明 |
|------|----------|----------|------|
| `fasta` | 查询 FASTA + 参考 FASTA | 常规模式 | 直接使用 FASTA 查询序列 |
| `fasta-cut` | 查询 FASTA + 参考 FASTA | 低内存模式 | 按染色体流式处理参考序列 |
| `anno` | 注释文件 + gene ID + 参考 FASTA | 常规模式 | 从注释提取查询序列后运行 |
| `anno-cut` | 注释文件 + gene ID + 参考 FASTA | 低内存模式 | 注释提取 + 按染色体低内存处理 |

### 普通模式与 cut 模式的区别

- **普通模式**：整份参考基因组加载进内存，速度通常更高。

- **cut 模式**：参考基因组按染色体逐条处理，峰值内存更低，更适合内存受限机器。

### 内存峰值估算

程序启动时会自动根据参考 FASTA 估算峰值内存，并要求用户确认；如需跳过，可使用：

```bash
-y
```
或
```bash
--yes
```
或
```bash
--skip-memory-check
```

- 普通模式峰值主要与**总基因组大小**相关
- cut 模式峰值主要与**最长染色体大小**相关

---

## 5. 命令行主程序：`spacerscope`

### 通用用法

```bash
spacerscope fasta <options>
spacerscope fasta-cut <options>
spacerscope anno <options>
spacerscope anno-cut <options>
spacerscope --report-tsv <file> --html-output <file>
```

查看帮助：

```bash
spacerscope --help
spacerscope help fasta
spacerscope help anno
```

查看版本：

```bash
spacerscope --version
```

### 通用参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-R, --ref <file>` | 参考基因组 FASTA 文件 | 必填 |
| `-o, --output <file>` | 输出 TSV 文件 | 必填 |
| `--pam <string>` | 主 PAM 序列 | `NGG` |
| `--alt-pam <string>` | 备用 PAM 序列 | `NAG` |
| `--spacer-len <int>` | spacer 长度 | `20` |
| `--direction <up|down>` | spacer 相对于 PAM 的方向 | `up` |
| `-m, --mismatch <int>` | 最大错配数 | `4` |
| `-i, --indel <int>` | 最大 indel 数 | `2` |
| `--dust-threshold <int>` | DUST 低复杂度过滤阈值 | `12` |
| `-t, --threads <int>` | 线程数 | 逻辑核心数 |
| `-y, --yes` | 跳过启动前内存确认 | 关闭 |
| `--skip-memory-check` | 跳过启动前内存确认 | 关闭 |
| `--raw-output` | 输出 7 列未评分 TSV | 关闭 |
| `--progress-format <text|jsonl>` | 进度输出格式 | `text` |
| `--html <0|1>` | 是否生成 `<output>.html` | `1` |
| `-ht <0|1>` | `--html` 的别名 | `1` |
| `--html-output <file>` | 指定 HTML 输出文件 | 空 |
| `--report-tsv <file>` | 对已有 scored TSV 单独生成 HTML | 空 |

### `fasta / fasta-cut` 模式参数

| 参数 | 说明 |
|------|------|
| `-I, --query <file>` | 查询 FASTA 文件 |

#### 示例

```bash
spacerscope fasta --query q.fa --ref ref.fa -o out.tsv
```

```bash
spacerscope fasta-cut --query q.fa --ref ref.fa -o out.tsv --yes
```

### `anno / anno-cut` 模式参数

| 参数 | 说明 |
|------|------|
| `-A, --annotation <file>` | 注释文件（GFF/GFF3） |
| `--geneID <id>` | 注释文件中的**精确 gene ID** |
| `--seq-type <gene|mrna|cds>` | 要提取的序列层级，默认 `gene` |

#### 示例

```bash
spacerscope anno --annotation genes.gff3 --geneID gene:Gene1 --seq-type cds --ref ref.fa -o out.tsv
```

```bash
spacerscope anno-cut --annotation genes.gff3 --geneID gene:Gene1 --seq-type gene --ref ref.fa -o out.tsv --yes
```

###  `anno` 模式的重要限制

当前版本的 `anno / anno-cut` 已经并入 C++ 主程序，但注释兼容策略是**按旧约定发布**，不是任意 GFF3 的通用解析器。

请注意：

1. `--geneID` 必须与注释文件中的 **gene ID 完全一致**，不是自由 gene symbol。
2. `mrna/cds` 提取遵循旧约定：转录本关联依赖 **`GeneX-mRNA` 命名关系**。
3. 对于通用 `Parent=transcript:...` 风格的 GFF3，当前版本 **不承诺** 在 `mrna/cds` 模式下完全兼容。

如果你的注释文件格式复杂或不统一，建议：

- 先手动提取目标基因序列
- 再使用 `fasta` 或 `fasta-cut` 模式运行

### `--raw-output` 与 HTML 报告

当前版本中：

- 正常输出：**9 列 scored TSV**
- `--raw-output`：输出 **7 列未评分合并 TSV**

`--raw-output` 生成的文件列为：

| 列名 |
|------|
| QueryName |
| QuerySeq |
| TargetName |
| TargetSeq |
| Distance |
| SearchType |
| GenomicFrequency |

注意：

- `--raw-output` 与 `--html-output` **不能同时使用**
- HTML 报告只支持基于 scored TSV 生成

### 仅根据已有 TSV 生成 HTML

如果已经有 scored TSV，可单独生成 HTML：

```bash
spacerscope --report-tsv out.tsv --html-output out.html
```

---

## 6. Web 版本（GUI）：`spacerscope-web`

### 适用场景

`spacerscope-web` 是当前版本提供的本地网页壳，适合：

- 无桌面环境的 Linux 服务器
- 希望通过浏览器配置参数和查看结果的用户
- 通过 SSH 端口转发访问本机服务的场景

### 启动方式

```bash
spacerscope-web
```

默认监听：

```text
127.0.0.1:8787
```

常用参数：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--bind <host>` | 绑定地址 | `127.0.0.1` |
| `--port <int>` | 监听端口 | `8787` |
| `--workspace-root <dir>` | Web 工作目录 | 系统临时目录下 `spacerscope-web` |
| `--spacerscope-bin <path>` | 指定 `spacerscope` 主程序路径 | 同目录自动查找 |
| `-h, --help` | 显示帮助 | - |
| `-V, --version` | 显示版本 | - |

#### 示例

```bash
spacerscope-web --bind 127.0.0.1 --port 8787
spacerscope-web 
```

### 访问方式

浏览器打开：

```text
http://127.0.0.1:8787
```

如果在远程 Linux 服务器上使用，推荐通过 SSH 端口转发访问：

```bash
ssh -L 8787:127.0.0.1:8787 user@server
```

然后在本地浏览器打开：

```text
http://127.0.0.1:8787
```

### Web 版本说明

Web 页面支持：

- `fasta / fasta-cut / anno / anno-cut`
- 服务器路径输入模式
- 文件上传模式
- 中英文切换
- 实时日志
- 实时进度
- TSV / HTML / stdout / stderr 下载

### 6.5 安全提示

当前 `spacerscope-web` v1 **无认证机制**。

因此：

- 默认只应绑定 `127.0.0.1`
- 若绑定到 `0.0.0.0` 或其他局域网地址，只应在**可信网络**中使用
- 不建议直接暴露到公网

---

## 7. 输出文件说明

### 默认主输出（scored TSV）

默认输出为 **9 列 scored TSV**：

| 列名 | 含义 |
|------|------|
| `QueryName` | 查询位点名称 |
| `QuerySeq` | 查询位点序列 |
| `TargetName` | 脱靶位点名称 |
| `TargetSeq` | 脱靶位点序列 |
| `Distance` | 编辑距离 |
| `SearchType` | 匹配类型：`Exact` / `Mismatch` / `Indel` |
| `GenomicFrequency` | 该序列在基因组中的出现次数 |
| `Score` | 脱靶评分 |
| `Details` | 比对细节 |

### HTML 报告

如果启用 HTML 输出，程序会额外生成：

```text
<output>.html
```

或你通过 `--html-output` 指定的文件。

HTML 报告包含：

- 每个 query 的统计汇总
- 完全匹配 / 错配 / indel 的数量分布
- 详细结果表格

### raw TSV

当使用 `--raw-output` 时，程序输出的是**7 列未评分 TSV**，不包含：

- `Score`
- `Details`

同时不会生成 HTML 报告。

---

## 8. 参数与运行建议

### 一般建议

- 大多数情况下，默认参数已经可用。
- 如无明确需要，不建议随意放大 `--indel`。
- 推荐先从默认参数开始，再根据需求微调。

### 关于 `--indel`

`--indel` 的增长会明显增加计算量。

通常建议：

- `--indel 0`：做完全匹配 + 错配
- `--indel 1~2`：常规推荐范围
- 更高的 `--indel`：仅在确有必要时尝试

### 关于 `--threads`

程序支持多线程，但线程数并非越大越好。实际速度会受到内存带宽、NUMA、机器架构等等的影响。

建议：

- 先用默认线程数或机器的常用高性能线程数测试
- 在服务器上实测确定最优值

### 8.4 关于 cut 模式

cut 模式的目标是**降低峰值内存**，不一定总是更快。内存足够时，普通模式通常更直接；内存受限时，cut 模式更安全。

---

## 9. 注意事项与已知限制

1. 建议使用时全部相关文件在**英文路径**下。
2. 不建议将整个基因组作为查询 FASTA 直接运行。
3. 建议使用**完整的连续基因序列**作为查询文件（例如采用 `fasta` 模式）。直接使用 mRNA 或 CDS 层级作为查询序列可能因剪接（Splicing）问题导致设计的 gRNA 序列在基因组中不存在，表现为完全匹配位点数量为 0 或未位于目标基因上。
4. 对注释文件：
   - 当前 `anno` 模式不是任意 GFF3 的通用解析器
   - `mrna/cds` 模式尤其依赖旧注释约定
5. 由于不同的基因注释文件中格式中不统一，建议提取出**基因序列**使用**fasta模式**。
6. `spacerscope-web` v1 无认证，不建议直接公网暴露。
7. `--progress-format jsonl` 主要用于 Web 壳和程序化调用，普通用户通常使用默认的 `text` 即可。
8. 选择最终的 gRNA 时，请确认所有完全匹配位点是否均位于**目标基因**上。
9. 在使用软件时，请给予软件足够的权限。

---

## 10. 常用命令速查

### FASTA 模式

```bash
spacerscope fasta --query query.fa --ref genome.fa -o result.tsv
```

### FASTA-CUT模式

```bash
spacerscope fasta-cut --query query.fa --ref genome.fa -o result.tsv --yes
```

### ANNO 模式

```bash
spacerscope anno --annotation genes.gff3 --geneID gene:Gene1 --seq-type gene --ref genome.fa -o result.tsv
```

### ANNO-CUT模式

```bash
spacerscope anno-cut --annotation genes.gff3 --geneID gene:Gene1 --seq-type cds --ref genome.fa -o result.tsv --yes
```

### 仅生成 HTML 报告

```bash
spacerscope --report-tsv result.tsv --html-output result.html
```

### 启动 Web 壳

```bash
spacerscope-web --bind 127.0.0.1 --port 8787
```

