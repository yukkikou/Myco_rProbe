# Myco_rProbe

**真菌 rRNA 探针设计与评估流程**

Myco_rProbe 是一个综合性生物信息学流程，用于设计、评估和优化物种特异性及泛真菌 rRNA 靶向杂交探针。它支持单物种和多物种探针设计、环状 RNA (circRNA) 鉴定以及基于 RNA-seq 数据的计算机模拟探针评估。

---

## 功能特点

- **单物种探针设计**：为单个真菌物种设计 rRNA 探针
- **多物种（泛真菌）探针设计**：设计靶向多个真菌物种保守 rRNA 区域的探针，并进行系统发育分组
- **环状 RNA (circRNA) 鉴定**：使用 CIRIquant 鉴定和定量 circRNA
- **计算机模拟探针评估**：利用 RNA-seq 比对数据评估探针特异性和覆盖度
- **分类决策树**：生成跨分类学术平的探针覆盖层级报告
- **优化探针选择**：基于覆盖度和数量评分自动选择最佳探针配置
- **核心 vs. 扩展探针分组**：区分靶向保守区域（核心）和可变区域（扩展）的探针

---

## 流程概览

本流程包含三个主要工作流，各步骤详解如下。

### A. 单物种探针设计 (`myco_rprobe_single.sh`)

```
基因组 FASTA ──► Barrnap (rRNA 预测) ──► Silva 数据库 (seqkit grep 物种注释)
    │
    ▼
rRNA 序列合并去重 ──► 滑动窗口 (探针片段生成)
    │
    ▼
CD-HIT-EST (冗余聚类, -c 阈值) ──► 二聚体过滤 (primer3, ΔG 阈值) ──► 最终探针
```

1. **rRNA 预测** — Barrnap 从输入基因组中鉴定 5S、5.8S、18S 和 28S rRNA 序列。
2. **Silva 数据库整合** — 通过 seqkit grep 将预测的 rRNA 与 Silva 参考数据库比对，向 FASTA 头部添加物种谱系信息。
3. **去重** — 基于序列一致性去除冗余（seqkit rmdup）。
4. **滑动窗口** — 以长度 `-l`、步长 `-w` 生成重叠探针片段。
5. **CD-HIT-EST 聚类** — 在阈值 `-c` 下去除冗余探针序列。
6. **二聚体过滤** — Primer3 评估二聚化自由能（ΔG），低于阈值的探针由 `remove_dimer_para.py` 去除。
7. **反向互补** — 最终探针序列反向互补为实验可用形式。

### B. 多物种（泛真菌）探针设计 (`myco_rprobe_multi.sh`)

```
基因组列表 (TSV) ──► 每个菌株: Barrnap + Silva DB ──► 合并所有 rRNA
    │
    ▼
按类型分离 rRNA (5S / 5.8S / 18S / 28S) 通过 seqkit grep
    │
    ▼
对每个 rRNA 配置行 (type, step, identity, width):
    │
    ├── 滑动窗口 (步长 = 配置步长, 宽度 = 配置宽度)
    ├── CD-HIT-EST 聚类 (阈值 = 配置一致性)
    ├── 内部探针去除 (remove_inner.py — 去除重叠探针)
    ├── BLAST 比对 Silva 参考 (覆盖度评估)
    ├── 覆盖度解析 (parse_coverage.py → rRNA_coverage_summary.txt)
    ├── 探针系统发育分组 (group_probes.py / group_probes_plus.py)
    │   └── 核心 (保守) vs 扩展 (可变) 分组
    └── 按分类级别的覆盖度统计 (grouped_coverage.py)
    │
    ▼
得分 = (1 / probe_num) × q3 ──► 自动选择最佳配置
    │
    ▼
决策树 (decision_tree.py) — 跨所有分类级别的层级报告:
    亚界 → 门 → 亚门 → 纲 → 亚纲 → 目 → 亚目 → 科 → 亚科
```

1. **多菌株 rRNA 提取** — 对基因组列表中每个菌株执行 Barrnap 和 Silva 注释，合并所有 rRNA 序列。
2. **按类型分离** — `all_rrna.fasta` 分割为 `sep_5s.fasta`、`sep_5_8s.fasta`、`sep_18s.fasta`、`sep_28s.fasta`。
3. **参数扫描** — rRNA 配置 (`-r`) 定义多组探针设计参数组合（如 `5S_40_0.8_40`、`28S_59_0.8_59`），每组在独立子目录中处理。
4. **探针生成与聚类** — 与单物种相同的滑动窗口、CD-HIT-EST 和内部探针去除步骤。
5. **覆盖度评估** — BLASTN 比对 Silva 参考，计算每个靶序列的覆盖度。
6. **系统发育分组** — 基于 BLAST 命中分布将探针分配至分类级别。核心探针靶向保守（高命中率）区域，扩展探针靶向可变（低命中率）区域。
7. **配置评分** — 每个配置评分为 `(1 / probe_num) × q3`，自动选取每类 rRNA 的最佳配置。
8. **决策树** — 生成跨所有分类级别的层级探针覆盖报告。

### C. 探针评估 (`myco_rprobe_evaluate.sh`)

```
RNA-seq (R1/R2) ──► fastp (质控) ──► 清洁读段
    │
    ├── Bowtie2 → rRNA 参考 (每类: 5S/5.8S/18S/28S) → 深度统计
    │       └── rRNA 比对率 → rrna_state.tsv
    │
    ├── Bowtie2 → 全基因组 → 排序 BAM
    │       └── 基因组比对率 → genome_state.tsv
    │
    └── 覆盖度图 (R: stat_coverage.R)
```

1. **质控** — fastp 去除接头序列，过滤低质量读段。
2. **rRNA 比对** — 清洁读段比对至各 rRNA 参考（Bowtie2），记录每碱基深度和比对率。
3. **基因组比对** — 读段比对至全基因组作为背景评估。
4. **覆盖度统计** — 汇总 rRNA/基因组比对率至 `rrna_state.tsv` 和 `genome_state.tsv`。
5. **可视化** — R 脚本生成每个菌株的覆盖度图。

### D. circRNA 鉴定 (`myco_rprobe_circiden.sh`)

```
非 rRNA 残留读段 ──► BWA 索引 + HISAT2 索引 ──► CIRIquant ──► circRNA GTF / BAM
```

### E. circRNA 序列提取 (`myco_rprobe_circseq.sh`)

```
circRNA GTF + 基因组 FASTA ──► bedtools getfasta ──► circRNA FASTA
```

---

## 环境要求

### 软件
- [Singularity](https://docs.sylabs.io/) (v3.x) - 容器化工具执行
- 所需容器：
  - `singularity/pan_fungi.sif` - 整合镜像，包含所有生物信息学工具（barrnap, seqkit, cd-hit-est, BLAST, Bowtie2, Samtools, fastp, HISAT2, BWA，以及 Python 环境含 Biopython/primer3-py/pandas/matplotlib/numpy，R 含 tidyverse/purrr）

### Python 依赖（容器内）
- Biopython
- primer3-py
- pandas
- matplotlib
- numpy

### R 依赖（容器内）
- tidyverse
- purrr

---

## 安装

```bash
git clone https://github.com/your-org/Myco_rProbe.git
cd Myco_rProbe
```

确保 Singularity 镜像文件位于指定路径。

---

## 使用方法

### 1. 单物种探针设计

```bash
bash myco_rprobe_single.sh \
    -s Aspergillus_fumigatus \
    -g example/Aspergillus_fumigatus.ASM265v1.dna.toplevel.fa.gz \
    -t 20 \
    -o single_af \
    -m singularity/pan_fungi.sif

| 参数 | 说明 |
|------|------|
| `-s` | 物种名（下划线分隔，如 `Aspergillus_fumigatus`） |
| `-g` | 基因组 FASTA 文件路径（支持 `.gz`） |
| `-t` | CPU 线程数 |
| `-o` | 输出目录 |
| `-m` | Singularity 镜像路径 |
| `-l` | 探针长度（默认：40 nt） |
| `-w` | 滑动窗口步长（默认：40 nt） |
| `-c` | CD-HIT-EST 序列一致性阈值（默认：0.8） |

### 2. 多物种泛真菌探针设计

```bash
bash myco_rprobe_multi.sh \
    -g example/fungi_genome.list \
    -t 10 \
    -o pathgen_fungi \
    -m singularity/pan_fungi.sif \
    -r data/rRNA_list.tsv
```

**基因组列表格式**（TSV，无表头）：
```
<基因组路径>  <属名>  <物种标识符>
```

**rRNA 配置格式**（TSV）：
```
<rRNA类型>  <滑动步长>  <一致性阈值>  <宽度>
```

示例：
```
5S    40    0.8    40
5S    50    0.85   50
18S   59    0.8    59
28S   59    0.8    59
```

### 3. 探针评估（使用 RNA-seq 数据）

```bash
bash myco_rprobe_evaluate.sh \
    -g example/batch_lib_config.tsv \
    -t 10 \
    -o batch_mapping \
    -m singularity/pan_fungi.sif \
    -r pathgen_fungi
```

**文库配置格式**（TSV）：
```
<菌株名>  <前缀>  <R1路径>  <R2路径>  <基因组路径>  <GTF路径>
```

### 4. circRNA 鉴定

```bash
bash myco_rprobe_circiden.sh \
    -g example/batch_lib_config.tsv \
    -t 10 \
    -o batch_mapping \
    -m singularity/pan_fungi.sif
```

### 5. circRNA 序列提取

```bash
bash myco_rprobe_circseq.sh \
    -g example/circ_seq.tsv \
    -t 5 \
    -o output_dir \
    -m singularity/pan_fungi.sif
```

---

## 输出结构

```
output_dir/
├── sep_5s.fasta              # 分离的 5S rRNA 序列
├── sep_5_8s.fasta            # 分离的 5.8S rRNA 序列
├── sep_18s.fasta             # 分离的 18S rRNA 序列
├── sep_28s.fasta             # 分离的 28S rRNA 序列
├── all_rrna.fasta            # 合并的所有 rRNA 序列
├── coverage_summary_results.txt  # 覆盖度汇总
│
├── {rna_type}_{step}_{identity}_{width}/
│   ├── sliding.fasta         # 滑动窗口探针
│   ├── cdhit_merge.fasta     # CD-HIT-EST 合并探针
│   ├── cdhit_final.fasta     # 过滤后的最终探针
│   ├── cdhit_final.alig.out  # BLAST 比对结果
│   ├── rrna_probes.fasta     # 探针序列（反向互补，RNA转DNA）
│   ├── rRNA_coverage_summary.txt  # 每个靶标的覆盖度
│   ├── probe_groups.tsv      # 探针系统发育分组
│   ├── probe_filtered.tsv    # 过滤后的最高水平分组
│   ├── fungal_tree_report.txt    # 分类决策树
│   └── {rna_type}_{level}_grouped.tsv  # 各分类水平的覆盖度
│
├── coverage/                 # 每个菌株的覆盖度数据
├── index/                    # Bowtie2/BLAST 索引
│   ├── rrna/                 # rRNA 参考索引
│   └── genome/               # 基因组索引
├── rrna_state.tsv            # rRNA 比对统计
└── genome_state.tsv          # 基因组比对统计
```

---

## 主要脚本

| 脚本 | 说明 |
|------|------|
| `reassemble.sh` | 从分卷重建 `singularity/pan_fungi.sif` |
| `myco_rprobe_single.sh` | 单物种探针设计 |
| `myco_rprobe_multi.sh` | 多物种泛真菌探针设计 |
| `myco_rprobe_evaluate.sh` | 基于 RNA-seq 的计算机模拟探针评估 |
| `myco_rprobe_circiden.sh` | 环状 RNA 鉴定 |
| `myco_rprobe_circseq.sh` | circRNA 序列提取 |
| `src/circ_blast.sh` | BLAST 数组作业脚本 |
| `src/plus_design/s01_run_mpile.sh` | 用于 plus 设计的 pileup 生成 |
| `src/plus_design/s02_mismatch_pileup.py` | 错配 pileup 分析 |

### Python 分析脚本 (src/)

| 脚本 | 说明 |
|------|------|
| `remove_dimer_para.py` | 使用 primer3 去除二聚体倾向探针（多线程） |
| `remove_inner.py` | 去除重叠的内部探针 |
| `parse_coverage.py` | 解析 BLAST 结果并计算覆盖度统计 |
| `group_probes.py` | 按系统发育水平分组探针 |
| `group_probes_plus.py` | 高级探针分组（含 NCBI 分类学） |
| `grouped_coverage.py` | 计算每个系统发育组的覆盖度 |
| `grouped_coverage_plus.py` | 高级覆盖度计算（含 NCBI 分类学） |
| `decision_tree.py` | 生成分层分类决策树报告 |
| `mismatch_pileup.py` | 分析 pileup 以生成共有序列和探针 |

---

---

## 引用

如果在研究中使用了 Myco_rProbe，请引用：

> (引用信息待添加)

---

## 许可证

本项目基于 Apache License 2.0 许可 - 详见 LICENSE 文件。

---

## 联系方式

如有问题或支持需求，请在 GitHub 上提交 issue。
